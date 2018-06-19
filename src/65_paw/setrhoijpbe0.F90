!{\src2tex{textfont=tt}}
!!****f* ABINIT/setrhoijpbe0
!! NAME
!! setrhoijpbe0
!!
!! FUNCTION
!! PAW local exact exchange only:
!! Impose value of rhoij for f electrons using an auxiliairy file
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FJ,MT,MD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  istep=index of the number of steps in the routine scfcv
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_comm_read=MPI communicator containing all the processes reading the PBE0 file
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell
!!  ntypat=number of types of atoms in unit cell
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type integer for each atom in cell
!!
!! SIDE EFFECTS
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Only valid for f electrons !!!
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,pawio_print_ij,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setrhoijpbe0(dtset,initialized,istep,istep_mix,&
&                       mpi_comm_read,my_natom,natom,ntypat,pawrhoij,pawtab,typat,& 
&                       mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_io_tools,   only : open_file
 use m_pawtab,     only : pawtab_type
 use m_pawrhoij,   only : pawrhoij_type
 use m_paw_io,     only : pawio_print_ij
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setrhoijpbe0'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: initialized,istep,mpi_comm_read,my_natom,natom,ntypat
 integer,intent(inout) :: istep_mix
 integer,optional,intent(in) :: comm_atom
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: ll=3
 integer :: iatom,iatom_tot,ierr,ii,ios,iread,irhoij,ispden,itypat,jj,klmn,my_comm_atom,my_rank,nselect
 integer :: nstep1,nstep1_abs,rhoijshft,rhoijsz
 logical :: my_atmtab_allocated,paral_atom,test0
 character(len=9),parameter :: filnam='rhoijpbe0'
 character(len=9),parameter :: dspin(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: strg, message
!arrays
 integer, allocatable :: nspden_tmp(:)
 integer,pointer :: my_atmtab(:)
 real(dp),allocatable :: rhoijtmp(:,:),rhoijtmp1(:,:),rhoijtmp2(:,:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Test existence of file and open it
 inquire(file=filnam,iostat=ios,exist=test0)
 if(.not.test0) return

!Look for parallelisation over atomic sites
 paral_atom=(present(comm_atom).and.(my_natom/=natom))

!Test if exact-exch. is on f electrons
 test0=.false.
 do itypat=1,ntypat
   if (pawtab(itypat)%useexexch>0.and.pawtab(itypat)%lexexch/=ll) test0=.true.
 end do
 if (test0) then
   write(message, '(3a,i1,a)' ) &
&   ' Local exact exchange: occ. matrix can only be imposed for l=',ll,' !'
   MSG_ERROR(message)
 end if

!============================================================
!===== First case: no parallelisation over atomic sites =====
!============================================================

 if (.not.paral_atom) then

!  Open file
   if (open_file(filnam,message,unit=77,form='formatted') /= 0) then
     MSG_ERROR(message)
   end if

!  Read step number and eventually exit
   nstep1=0;test0=.false.
   do while (.not.test0)
     read(77,'(A)') strg
     test0=(strg(1:1)/="#")
     if (test0) read(unit=strg,fmt=*) nstep1
   end do
   nstep1_abs=abs(nstep1)
   if (nstep1_abs==0.or.istep>nstep1_abs.or.(nstep1>0.and.initialized/=0)) then
     close(77)
!    Reinitalize mixing when rhoij is allowed to change; for experimental purpose...
     if (dtset%userib==1234.and.istep==1+nstep1_abs.and.(nstep1<0.or.initialized==0)) istep_mix=1
     return
   end if

!  Loop on atoms
   do iatom=1,natom
     itypat=typat(iatom)
     if (pawtab(itypat)%useexexch>0) then

!      Set sizes depending on ll
       rhoijsz=4*ll+2
       rhoijshft=2*ll*ll
       
!      Uncompress rhoij
       ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
         end do
       end do
!      Read rhoij from file
       ABI_ALLOCATE(rhoijtmp1,(rhoijsz,rhoijsz))
       do ispden=1,pawrhoij(iatom)%nspden
         do ii=1,rhoijsz
           test0=.false.
           do while (.not.test0)
             read(77,'(A)') strg
             test0=(strg(1:1)/="#")
             if (test0)  read(unit=strg,fmt=*) (rhoijtmp1(ii,jj), jj=1,rhoijsz)
           end do
         end do

!        Impose rhoij
         do jj=1,rhoijsz
           do ii=1,jj
             rhoijtmp((jj+rhoijshft)*((jj+rhoijshft)-1)/2+ii+rhoijshft,ispden)=rhoijtmp1(ii,jj)
           end do
         end do

       end do
       ABI_DEALLOCATE(rhoijtmp1)

!      Compress rhoij
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
           end do
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)

!      Print new rhoij
       do ispden=1,pawrhoij(iatom)%nspden
         write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,&
&         ' == Imposed occupation matrix'
         if (pawrhoij(iatom)%nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
         if (pawrhoij(iatom)%nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
         if (pawrhoij(iatom)%nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&         trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4)))," =="
         call wrtout(std_out,message,'COLL')
         call pawio_print_ij(std_out,pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&         pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,ll,&
&         pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size),&
&         1,-1,pawrhoij(iatom)%rhoijselect(:),-1.d0,1,mode_paral='COLL')
       end do

!      End loop on atoms
     end if
   end do

!  Close file
   close (77)

 else

!  ============================================================
!  ====== 2nd case: no parallelisation over atomic sites =====
!  ============================================================
   
   my_rank=xmpi_comm_rank(mpi_comm_read)

!  Read step number and eventually exit
   iread=0
   if (my_rank==0) then
     if (open_file(filnam,message,unit=77,form='formatted') /=0 ) then
       MSG_ERROR(message)
     end if
     nstep1=0;test0=.false.
     do while (.not.test0)
       read(77,'(A)') strg
       test0=(strg(1:1)/="#")
       if (test0) read(unit=strg,fmt=*) nstep1
     end do
     nstep1_abs=abs(nstep1)
     if (nstep1_abs==0.or.istep>nstep1_abs.or.(nstep1>0.and.initialized/=0)) then
       close(77)
!      Reinitalize mixing when rhoij is allowed to change; for experimental purpose...
       if (dtset%userib==1234.and.istep==1+nstep1_abs.and.(nstep1<0.or.initialized==0)) istep_mix=1
       iread=1
     end if
   end if
   call xmpi_sum(iread,mpi_comm_read,ierr)
   if (iread/=0) return

!  Set up parallelism over atoms
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!  Store number of component for rhoij
   ABI_ALLOCATE(nspden_tmp,(natom))
   nspden_tmp(:)=zero
   do iatom=1,my_natom
     iatom_tot=my_atmtab(iatom)
     nspden_tmp(iatom_tot)=pawrhoij(iatom)%nspden
   end do
   call xmpi_sum(nspden_tmp,mpi_comm_read,ierr)
   
!  To be improve if too much memory
   ABI_ALLOCATE(rhoijtmp2,(natom,4,rhoijsz,rhoijsz))
   rhoijtmp2=zero 

!  Read rhoij from file
   if (my_rank==0) then
     do iatom=1,natom
       itypat=typat(iatom)
       if (pawtab(itypat)%useexexch>0) then
         rhoijsz=4*ll+2
         do ispden=1,nspden_tmp(iatom)
           do ii=1,rhoijsz
             test0=.false.
             do while (.not.test0)
               read(77,'(A)') strg
               test0=(strg(1:1)/="#")
               if (test0)  read(unit=strg,fmt=*) (rhoijtmp2(iatom,ispden,ii,jj),jj=1,rhoijsz)
             end do
           end do
         end do
       end if
     end do
   end if
   call xmpi_sum(rhoijtmp2,mpi_comm_read,ierr)

!  Now, distribute rhoij
   do iatom=1,my_natom
     iatom_tot=my_atmtab(iatom)
     itypat=pawrhoij(iatom)%itypat

     if (pawtab(itypat)%useexexch>0) then

!      Set sizes depending on ll
       rhoijsz=4*ll+2
       rhoijshft=2*ll*ll

!      Uncompress rhoij
       ABI_ALLOCATE(rhoijtmp,(pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=pawrhoij(iatom)%rhoijselect(irhoij)
           rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
         end do

!        Impose rhoij
         do jj=1,rhoijsz
           do ii=1,jj
             rhoijtmp((jj+rhoijshft)*((jj+rhoijshft)-1)/2+ii+rhoijshft,ispden)=rhoijtmp2(iatom_tot,ispden,ii,jj)
           end do
         end do

       end do

!      Compress rhoij
       nselect=0
       do klmn=1,pawrhoij(iatom)%lmn2_size
         if (any(abs(rhoijtmp(klmn,:))>tol10)) then
           nselect=nselect+1
           do ispden=1,pawrhoij(iatom)%nspden
             pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
           end do
           pawrhoij(iatom)%rhoijselect(nselect)=klmn
         end if
       end do
       pawrhoij(iatom)%nrhoijsel=nselect
       ABI_DEALLOCATE(rhoijtmp)

     end if ! useexexch>0

!    Print new rhoij
     do ispden=1,pawrhoij(iatom)%nspden
       write(message,'(2a,i3,a)') ch10,'== Atom ',iatom,' == Imposed occupation matrix'
       if (pawrhoij(iatom)%nspden==1) write(message,fmt='(2a)')     trim(message)," for spin up =="
       if (pawrhoij(iatom)%nspden==2) write(message,fmt='(2a,i3,a)')trim(message)," for spin ",ispden," =="
       if (pawrhoij(iatom)%nspden==4) write(message,fmt='(4a)')     trim(message)," for component ", &
&       trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4)))," =="
       call wrtout(std_out,message,'PERS')
       call pawio_print_ij(std_out,pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&       pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,ll,&
&       pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size),&
&       1,-1,pawrhoij(iatom)%rhoijselect(:),-1.d0,1,mode_paral='PERS')
     end do

!    end loop on atoms
   end do

   ABI_DEALLOCATE(nspden_tmp)
   ABI_DEALLOCATE(rhoijtmp2)

!  Destroy atom table used for parallelism
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)

!  ============================================================
 end if ! paral_atom

 DBG_EXIT("COLL")

end subroutine setrhoijpbe0
!!***
