!!****p* ABINIT/fold2Bloch
!! NAME
!! fold2Bloch
!!
!! FUNCTION
!! Main routine for the unfolding of the wavefuntion.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2020 ABINIT group (AB)
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
!!      abi_io_redirect,abimem_init,abinit_doctor,crystal_free,crystal_from_hdr
!!      ebands_free,getargs,newk,progress,prompt,sortc,wfk_close,wfk_open_read
!!      wfk_read_band_block,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program fold2Bloch

 use defs_basis
 use m_errors
 use m_abicore
 use m_wfk
 use m_xmpi
 use m_nctk
 use m_hdr
 use m_crystal
 use m_ebands
 use m_fold2block
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,       only : strcat
 use m_io_tools,       only : get_unit, iomode_from_fname, open_file, prompt
 use defs_datatypes,   only : ebands_t
implicit none

!Arguments --------------------------------------------------------------

!Local variables-------------------------------
!scalars
integer :: ikpt, iband,nspinor,nsppol,mband,nkpt,mcg,csppol, cspinor, nfold, iss, ii
integer :: comm, my_rank, nargs, iomode, ncid, ncerr, fform, timrev, kunf_varid, weights_varid, eigunf_varid
integer :: cg_b, count, outfile, outfile1, outfile2, lwcg, hicg, pos
character(fnlen) :: fname, outname,seedname
character(len=500) :: msg
type(wfk_t) :: wfk
type(crystal_t) :: cryst
type(ebands_t) :: ebands
!arrays
integer :: folds(3),fold_matrix(3,3)
integer, allocatable :: kg(:,:),nband(:), npwarr(:)
real(dp), allocatable :: cg(:,:), eig(:),kpts(:,:), weights(:),coefc(:,:), nkval(:,:)

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

 nargs = command_argument_count()

 if (nargs == 0) then
   call prompt("Enter WFK file name:", fname)
   call prompt("Enter x y z integers giving the multiplicity:", folds)
 else
   call getargs(folds, fname) !Process command line arguments
 end if
 ! Use fold_matrix instead of folds(1:3) to prepare possible generalization.
 fold_matrix = 0
 do ii=1,3
   fold_matrix(ii,ii) = folds(ii)
 end do

 ! Test if the netcdf library supports MPI-IO
 !call nctk_test_mpiio()

 if (nctk_try_fort_or_ncfile(fname, msg) /= 0) then
   MSG_ERROR(msg)
 end if

 pos=INDEX(fname, "_")
 write(seedname,'(a)') fname(1:pos-1)

 write(std_out,*) '         '//achar(27)//'[97m ***********************' !print program header in pearl white
 write(std_out,*) '          ** Fold2Bloch V 1.1  **'
 write(std_out,*) '          **Build  Mar 16, 2015**'
 write(std_out,*) '          ***********************'//achar(27)//'[0m'

 ebands = wfk_read_ebands(fname, xmpi_comm_self)
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
 nfold = product(folds)

#ifdef HAVE_NETCDF
 timrev = 2; if (any(wfk%hdr%kptopt == [3, 4])) timrev = 1
 cryst = wfk%hdr%get_crystal(timrev)

 NCF_CHECK(nctk_open_create(ncid, strcat(seedname, "_FOLD2BLOCH.nc"), xmpi_comm_self))
 fform = fform_from_ext("FOLD2BLOCH.nc")
 NCF_CHECK(wfk%hdr%ncwrite(ncid, fform, nc_define=.True.))
 NCF_CHECK(cryst%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(ebands, ncid))

 ncerr = nctk_def_dims(ncid, [ &
 nctkdim_t("nk_unfolded", nkpt * nfold), &
 nctkdim_t("nsppol_times_nspinor", wfk%hdr%nsppol * wfk%hdr%nspinor)], defmode=.True.)
 NCF_CHECK(ncerr)
 ncerr = nctk_def_arrays(ncid, [ &
 nctkarr_t("fold_matrix", "int", "number_of_reduced_dimensions, number_of_reduced_dimensions"), &
 nctkarr_t("reduced_coordinates_of_unfolded_kpoints", "dp", "number_of_reduced_dimensions, nk_unfolded"), &
 nctkarr_t("unfolded_eigenvalues", "dp", "max_number_of_states, nk_unfolded, number_of_spins"), &
 nctkarr_t("spectral_weights", "dp", "max_number_of_states, nk_unfolded, nsppol_times_nspinor") &
 ])
 NCF_CHECK(ncerr)
 NCF_CHECK(nf90_inq_varid(ncid, "reduced_coordinates_of_unfolded_kpoints", kunf_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "unfolded_eigenvalues", eigunf_varid))
 NCF_CHECK(nf90_inq_varid(ncid, "spectral_weights", weights_varid))
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "fold_matrix"), fold_matrix))
 call cryst%free()
#endif

 call ebands_free(ebands)

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
     ABI_ALLOCATE(weights, (nfold))
     ABI_ALLOCATE(nkval,(3, nfold))
     call progress(ikpt,nkpt,kpts(:,ikpt)) !Write progress information

     !Read a block of data
     call wfk%read_band_block([1, nband(ikpt)], ikpt, csppol, xmpio_single, kg_k=kg, cg_k=cg, eig_k=eig)

     !Determine unfolded K point states
     call newk(kpts(1,ikpt),kpts(2,ikpt),kpts(3,ikpt),folds(1),folds(2),folds(3),nkval)
#ifdef HAVE_NETCDF
     if (csppol == 1) then
       NCF_CHECK(nf90_put_var(ncid, kunf_varid, nkval, start=[1, 1 + (ikpt-1) * nfold], count=[3, nfold]))
     end if
#endif

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
         ! Write out results, format: new k states(x, y, and z), eigenvalue, weight
         do count=1, nfold
           write(outfile,50) nkval(1,count),nkval(2,count),nkval(3,count),eig(iband),weights(count)
           50 format(f11.6, f11.6, f11.6, f11.6, f11.6)
         end do
#ifdef HAVE_NETCDF
         iss = csppol; if (nspinor == 2) iss = cspinor
         ncerr = nf90_put_var(ncid, weights_varid, weights, start=[iband, 1 + (ikpt-1) * nfold, iss], &
         stride=[mband, 1, 1], count=[1, nfold, 1])
         NCF_CHECK(ncerr)
         if (cspinor == 1) then
           weights = eig(iband) ! Use weights as workspace array.
           ncerr = nf90_put_var(ncid, eigunf_varid, weights, start=[iband, 1 + (ikpt-1) * nfold, csppol], &
           stride=[mband, 1, 1], count=[1, nfold, 1])
             !count=[1, nfold, 1])
           NCF_CHECK(ncerr)
         end if
#endif
       end do ! cspinor
       cg_b=cg_b+nspinor*npwarr(ikpt) !shift coefficient pointer for next eigenvalue
     end do ! iband

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
 call wfk%close()

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

#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_close(ncid))
#endif

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__fold2bloch")

 call xmpi_end()

 end program fold2Bloch
!!***
