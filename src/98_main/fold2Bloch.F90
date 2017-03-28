!{\src2tex{textfont=tt}}
!!****p* ABINIT/fold2Bloch
!! NAME
!! fold2Bloch
!!
!! FUNCTION
!! Main routine for the unfolding of the wavefuntion.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2017 ABINIT group (AB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (main program)
!!
!! NOTES
!! folds= Array of folds in X,Y, and Z directions
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,getargs,newk,progress,prompt
!!      sortc,wfk_close,wfk_open_read,wfk_read_band_block,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program fold2Bloch

 use defs_basis
 use defs_abitypes
 use m_distribfft
 use m_errors
 use m_profiling_abi
 use m_splines
 use m_wfk
 use m_xmpi
 use m_nctk
 use m_hdr
 use m_fold2block

 use m_io_tools, only : get_unit, iomode_from_fname, open_file, prompt, file_exists
 use m_pptools,  only : print_fofr_ri, print_fofr_xyzri , print_fofr_cube

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fold2Bloch'
!End of the abilint section

implicit none

!Arguments --------------------------------------------------------------

!Local variables-------------------------------
type(wfk_t) :: Wfk
integer :: ikpt, iband,nspinor,nsppol,mband,nkpt,mcg,csppol, cspinor 
integer :: comm, my_rank, nargs, iomode
integer :: folds(3),cg_b, count, outfile, outfile1, outfile2, lwcg, hicg, pos !,unitfi
real(dp), allocatable :: cg(:,:), eig(:),kpts(:,:), weights(:),coefc(:,:), nkval(:,:)
integer, allocatable :: kg(:,:),nband(:), npwarr(:)
character(fnlen) :: fname, outname,seedname
character(len=500) :: msg

!*************************************************************************

!0) Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(xmpi_world)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 if (xmpi_comm_size(comm) /= 1) then
   MSG_ERROR("fold2bloch not programmed for parallel execution.")
 end if

 ! Test if the netcdf library supports MPI-IO
 !call nctk_test_mpiio()

 nargs = command_argument_count()

 if (nargs == 0) then
   call prompt("Enter WFK file name:", fname)
   call prompt("Enter x y z integers giving the multiplicity:", folds)
 else
   call getargs(folds, fname) !Process command line arguments
 end if

 if (nctk_try_fort_or_ncfile(fname, msg) /= 0) then
   MSG_ERROR(msg)
 end if

 pos=INDEX(fname,"_") 
 !read(fname(1:pos-1),*) seedname ! File name root
 write(seedname,'(a)') fname(1:pos-1) 
!folds=(/1,2,3/)

 write(std_out,*) '         '//achar(27)//'[97m ***********************' !print program header in pearl white
 write(std_out,*) '          ** Fold2Bloch V 1.1  **'
 write(std_out,*) '          **Build  Mar 16, 2015**'
 write(std_out,*) '          ***********************'//achar(27)//'[0m'

 iomode = iomode_from_fname(fname)
 call wfk_open_read(wfk,fname,0,iomode,get_unit(),comm)

 nkpt=wfk%hdr%nkpt
 ABI_ALLOCATE(npwarr,(nkpt))
 ABI_ALLOCATE(nband,(nkpt))
 ABI_ALLOCATE(kpts,(3,nkpt)) 

 nsppol=wfk%hdr%nsppol
 nspinor=wfk%hdr%nspinor
 npwarr=wfk%hdr%npwarr
 kpts=wfk%hdr%kptns
 nband=wfk%hdr%nband
 mband=maxval(nband)
 mcg=maxval(npwarr)*nspinor*mband

 do csppol=1, nsppol
   if (nsppol==1) then !Determine spin polarization for output file
     outname=trim(seedname)//".f2b"
   elseif ((nsppol==2).and.(csppol==1)) then
     outname=trim(seedname)//"_UP.f2b"
     write(std_out,*) "     ===================="
     write(std_out,*) "     SPIN POLARIZATION UP"
     write(std_out,*) "     ===================="
   elseif ((nsppol==2).and.(csppol==2)) then
     outname=trim(seedname)//"_DOWN.f2b"
     write(std_out,*) "     ======================"
     write(std_out,*) "     SPIN POLARIZATION DOWN"
     write(std_out,*) "     ======================"
   end if
   if (nspinor==2) then
     !open output file
     if (open_file(trim(seedname)//"_SPOR_1.f2b", msg, newunit=outfile1, form="formatted", status="unknown") /= 0) then
       MSG_ERROR(msg)
     end if 
     if (open_file(trim(seedname)//"_SPOR_2.f2b", msg, newunit=outfile2, form="formatted", status="unknown") /= 0) then
       MSG_ERROR(msg)
     end if
   else
     if (open_file(outname, msg, newunit=outfile1,form="formatted", status="unknown") /= 0) then
       MSG_ERROR(msg)
     end if 
   end if

   do ikpt=1, nkpt !For each K point
     ABI_ALLOCATE(cg,(2,mcg))
     ABI_ALLOCATE(eig,((2*mband)**0*mband))
     ABI_ALLOCATE(kg,(3,npwarr(ikpt)))
     ABI_ALLOCATE(coefc,(2,nspinor*npwarr(ikpt)))
     ABI_ALLOCATE(weights,(product(folds)))
     ABI_ALLOCATE(nkval,(3,product(folds)))
     call progress(ikpt,nkpt,kpts(:,ikpt)) !Write progress information

     !Read a block of data
     call wfk_read_band_block(wfk, (/1,nband(ikpt)/), ikpt, csppol, xmpio_single, kg_k=kg, cg_k=cg, eig_k=eig)

     !Determine unfolded K point states
     call newk(kpts(1,ikpt),kpts(2,ikpt),kpts(3,ikpt),folds(1),folds(2),folds(3),nkval)

     cg_b=1
     do iband=1, nband(ikpt) !Foe each Eigenvalue
       coefc=cg(:,cg_b:(cg_b+nspinor*npwarr(ikpt)-1)) !Split coefficients per eigen value according to the number of "kg" 
       do cspinor=1,nspinor
         if (cspinor==1) then
           outfile=outfile1
           lwcg=1
           hicg=npwarr(ikpt)
         else
          ! Move coefficient span to spinor 2
           outfile=outfile2
           lwcg=npwarr(ikpt)+1
           hicg=npwarr(ikpt)*nspinor
         end if
         call sortc(folds(1),folds(2),folds(3),kg,coefc(:,lwcg:hicg),npwarr(ikpt),weights)  
        !Write out results, format: new k states(x, y, and z), eigenvalue, weight
         do count=1, product(folds) 
           write(outfile,50) nkval(1,count),nkval(2,count),nkval(3,count),eig(iband),weights(count)
           50 format(f11.6, f11.6, f11.6, f11.6, f11.6)
         end do
       end do
       cg_b=cg_b+nspinor*npwarr(ikpt) !shift coefficient pointer for next eigenvalue
     end do
     ABI_DEALLOCATE(cg)
     ABI_DEALLOCATE(eig)
     ABI_DEALLOCATE(kg)
     ABI_DEALLOCATE(coefc)
     ABI_DEALLOCATE(weights)
     ABI_DEALLOCATE(nkval)
   end do
   if (nspinor==2) then
     close(outfile1) !close output file
     close(outfile2)
   else 
     close(outfile1)
   end if
 end do
 call wfk_close(wfk)

 ABI_FREE(kpts)
 ABI_FREE(nband)
 ABI_FREE(npwarr)

! Print summary
 write(std_out,*) '    '//achar(27)//'[97m Number of K points processed:', nkpt
 if (nsppol==2) then
   write(std_out,*) '     Data was written to: ', trim(seedname)//"_UP.f2b", " & ", trim(seedname)//"_DOWN.f2b"
 else
   if (nspinor==2) then
     write(std_out,*) '     Data was written to: ', trim(seedname)//"_SPOR_1.f2b", " & ", trim(seedname)//"_SPOR_2.f2b"
   else
     write(std_out,*) '     Data was written to: ', trim(seedname)//".f2b"
   end if
 end if    
 write(std_out,*) '     Data format: KX, KY, KZ, Eigenvalue(Ha), Weight'//achar(27)//'[0m'

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__fold2bloch")

 call xmpi_end()

 end program fold2Bloch
!!***
