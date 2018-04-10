!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_all_gkq
!!
!! NAME
!! get_all_gkq
!!
!! FUNCTION
!! This routine determines what to do with the initial qspace
!!   matrix elements of the electron phonon coupling (to disk or in memory),
!!   then reads those given in the gkk file and completes them
!!   (for kpts, then perturbations)
!!   01/2010: removed completion on qpoints here (MJV)
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   Cryst<crystal_t>=Info on the unit cell and on its symmetries.
!!   Ifc<ifc_type>=Object containing the interatomic force constants.
!!   Bst<ebands_t>=GS energies, occupancies and Fermi level.
!!   FSfullpqtofull = mapping of k+q to another k
!!   kphon_full2full = mapping of FS kpoints under symops
!!   kpt_phon = fermi surface kpoints
!!   %k_phon%wtk = integration weights for bands and kpoints near the FS
!!   gkk_flag = flag to
!!   nband = number of bands
!!   n1wf = number of file headers from perturbation calculations
!!      which are present in the initial gkk input file.
!!   onegkksize = size of one record of the new gkk output file, in bytes
!!   qpttoqpt = mapping of qpoints onto each other under symmetries
!!   unitgkk = fortran unit for initial gkk input file
!!   xred = reduced coordinates of atoms
!!
!! OUTPUT
!!   elph_ds%gkq = recip space elphon matrix elements.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      complete_gkk,int2char4,read_gkk,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_all_gkq (elph_ds,Cryst,ifc,Bst,FSfullpqtofull,nband,n1wf,onegkksize,&
&    qpttoqpt,ep_prt_yambo,unitgkk,ifltransport)

 use defs_basis
 use defs_datatypes
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_io_tools
 use m_xmpi

 use m_fstrings, only : int2char4
 use m_crystal,  only : crystal_t
 use m_ifc,      only : ifc_type
 use m_iogkk,    only : read_gkk

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_all_gkq'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => get_all_gkq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1wf,nband,onegkksize,unitgkk,ep_prt_yambo,ifltransport
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: Bst
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 integer,intent(in) :: qpttoqpt(2,Cryst%nsym,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: iost,ierr,me,sz2,sz3,sz4,sz5,sz6
 character(len=10) :: procnum
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer,allocatable :: gkk_flag(:,:,:,:,:)

! *************************************************************************

!attribute file unit number
 elph_ds%unitgkq = get_unit()

!============================================
!save gkk for all qpts in memory or to disk
!============================================

!DEBUG
!write(std_out,*) ' 4 bytes / ??'
!write(std_out,*) ' kind(real) = ', kind(one)
!write(std_out,*) ' elph_ds%ngkkband = ', elph_ds%ngkkband, '^2'
!write(std_out,*) ' elph_ds%nbranch = ', elph_ds%nbranch, '^2'
!write(std_out,*) ' elph_ds%k_phon%nkpt = ', elph_ds%k_phon%nkpt
!write(std_out,*) ' elph_ds%nsppol = ', elph_ds%nsppol
!write(std_out,*) ' elph_ds%nqptirred ', elph_ds%nqptirred
!ENDDEBUG

 write(message,'(a,f14.4,a)')&
& ' get_all_gkq : gkq file/array size = ',&
 4.0*dble(onegkksize)*dble(elph_ds%k_phon%my_nkpt)*dble(elph_ds%nqptirred)/1024.0_dp/1024.0_dp/1024.0_dp,' Gb'
 call wrtout(std_out,message,'COLL')

 if (elph_ds%gkqwrite == 0) then !calculate gkk(q) keeping all in memory

   call wrtout(std_out,' get_all_gkq : keep gkk(q) in memory ','COLL')

   sz2=elph_ds%ngkkband*elph_ds%ngkkband
   sz3=elph_ds%nbranch*elph_ds%nbranch
   sz4=elph_ds%k_phon%my_nkpt
   sz5=elph_ds%nsppol
   if (ifltransport == 3) then
     sz6=elph_ds%nqpt_full
   else
     sz6=elph_ds%nqptirred
   end if
   ABI_STAT_ALLOCATE(elph_ds%gkk_qpt,(2,sz2,sz3,sz4,sz5,sz6), ierr)
   ABI_CHECK(ierr==0, 'Trying to allocate array elph_ds%gkk_qpt')

   elph_ds%gkk_qpt = zero

 else if (elph_ds%gkqwrite == 1) then !calculate gkk(q) and write to file
   me = xmpi_comm_rank(xmpi_world)
   call int2char4(me,procnum)
   ABI_CHECK((procnum(1:1)/='#'),'Bug: string length too short!')
   fname=trim(elph_ds%elph_base_name) // "_P" // trim(procnum) // '_GKKQ'

   iost=open_file(file=fname,iomsg=message,newunit=elph_ds%unitgkq,access='direct',&
&   recl=onegkksize,form='unformatted')
   if (iost /= 0) then
     write (message,'(2a)')' get_all_gkq : ERROR- opening file ',trim(fname)
     MSG_ERROR(message)
   end if

   write (message,'(5a)')&
&   ' get_all_gkq : gkq matrix elements  will be written to file : ',trim(fname),ch10,&
&   ' Nothing is in files yet',ch10
   call wrtout(std_out,message,'COLL')

 else
   write(message,'(a,i0)')' gkqwrite must be 0 or 1 while it is : ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if !if gkqwrite

!=====================================================
!read in g_kk matrix elements for all bands, kpoints,
!and calculated qpoints
!=====================================================
 call wrtout(std_out,' get_all_gkq : calling read_gkk to read in the g_kk matrix elements',"COLL")

 sz2=elph_ds%nbranch;sz3=elph_ds%k_phon%my_nkpt
 sz4=elph_ds%nsppol;sz5=elph_ds%nqpt_full
 ABI_STAT_ALLOCATE(gkk_flag,(sz2,sz2,sz3,sz4,sz5), ierr)
 ABI_CHECK(ierr==0, "allocating gkk_flag")

 call read_gkk(elph_ds,Cryst,ifc,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,ep_prt_yambo,unitgkk)

!if (elph_ds%symgkq ==1) then
!MJV 01/2010 removed the completion on qpt here: it should be done after FS integration
!so that everything is lighter in memory etc... (only irred qpt)
! if (0==1) then
 if (ifltransport == 3) then !  bxu, complete gkk is necessary

!  ==============================================================
!  complete gkk matrices for other qpoints on the full grid qpt_full
!  inspired and cannibalized from symdm9.f
!  FIXME: should add the possibility to copy over to other qpoints,
!  without full symmetrization, for testing purposes.
!  ==============================================================

   write(message,'(4a)')ch10,&
&   ' get_all_gkq : calling complete_gkk to complete ',ch10,&
&   ' gkk matrices for other qpoints on the full grid'
   call wrtout(std_out,message,'COLL')

   call complete_gkk(elph_ds,gkk_flag,Cryst%gprimd,Cryst%indsym,&
&   Cryst%natom,Cryst%nsym,qpttoqpt,Cryst%rprimd,Cryst%symrec,Cryst%symrel)

   call wrtout(std_out,' get_all_gkq : out of complete_gkk','COLL')

 end if !symgkq

!TODO Do we need gkk_flag in elphon?
 ABI_DEALLOCATE(gkk_flag)

end subroutine get_all_gkq
!!***
