!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrgdv
!! NAME
!! mrgdv
!!
!! FUNCTION
!! This program merges DFPT potentials for different q-vectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! DVDB file format:
!!   version (integer)
!!   number of potentials (integer)
!!   for each potential:
!!     Abinit header with info on the perturbation and the FFT mesh
!!     potential on the FFT mesh 
!!
!! CHILDREN
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,dvdb_free,dvdb_ftinterp_qpt
!!      dvdb_ftinterp_setup,dvdb_init,dvdb_list_perts,dvdb_merge_files
!!      dvdb_open_read,dvdb_print,dvdb_readsym_allv1,dvdb_test_symmetries
!!      get_command_argument,herald,ngfft_seq,prompt,vdiff_print,wrtout
!!      xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgdv

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_build_info
 use m_profiling_abi
 use m_dvdb

 use m_fstrings,        only : sjoin, itoa, ktoa
 use m_numeric_tools,   only : vdiff_eval, vdiff_print
 use m_io_tools,        only : file_exists, prompt
 use m_fftcore,         only : ngfft_seq

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgdv'
 use interfaces_14_hidewrite
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 integer :: ii,nargs,nfiles,comm,prtvol,nfft,ifft
 integer :: iq,cplex,mu,ispden,my_rank
 character(len=24) :: codename
 character(len=500) :: command,arg !msg,
 character(len=fnlen) :: db_path
 type(dvdb_t) :: db
!arrays
 integer :: ngfft(18),ngqpt(3)
 real(dp),allocatable :: file_v1r(:,:,:,:),intp_v1r(:,:,:,:),tmp_v1r(:,:,:,:)
 character(len=fnlen),allocatable :: v1files(:)

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world
 my_rank = xmpi_comm_rank(comm)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 ! write greating,read the file names, etc.
 codename='MRGDV'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

 ABI_CHECK(xmpi_comm_size(comm) == 1, "Not programmed for parallel execution")
 prtvol = 0

 nargs = command_argument_count()

 if (nargs == 0) then
   ! We are reading from stdin
   call prompt("Enter name of output file:", db_path)
   call prompt("Enter total number of DFPT POT files:", nfiles)
   ABI_MALLOC(v1files, (nfiles))
   do ii=1,nfiles
     call prompt(sjoin("Enter name of POT file",itoa(ii),":"), v1files(ii))
   end do
   call dvdb_merge_files(nfiles, v1files, db_path, prtvol)
   ABI_FREE(v1files)

 else
   ! Command line options.
   do ii=1,command_argument_count()
     call get_command_argument(ii, arg)
     if (arg == "-v" .or. arg == "--version") then
       write(std_out,"(a)") trim(abinit_version); goto 100

     else if (arg == "-h" .or. arg == "--help") then
       ! Document the options.
       write(std_out,*)"-v, --version              Show version number and exit."
       write(std_out,*)"merge out_DVDB POT1 POT2   Merge list of POT files, produce out_DVDB file."
       write(std_out,*)"info out_DVDB              Print information on DVDB file"
       write(std_out,*)"-h, --help                 Show this help and exit."
       goto 100
     end if
   end do

   call get_command_argument(1, command)

   select case (command)
   case ("merge")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, db_path)
     if (file_exists(db_path)) then
       MSG_ERROR(sjoin("Cannot overwrite existing file:", db_path))
     end if

     nfiles = nargs - 2
     ABI_MALLOC(v1files, (nfiles))
     do ii=1,nfiles
       call get_command_argument(ii+2, v1files(ii))
     end do

     call dvdb_merge_files(nfiles, v1files, db_path, prtvol)
     ABI_FREE(v1files)

   case ("info")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, db_path)

     call dvdb_init(db, db_path, comm)

     call dvdb_print(db)
     call dvdb_list_perts(db, [-1,-1,-1])
     !call dvdb_list_perts(db, [2, 2, 2])

     !call ngfft_seq(ngfft, [12, 12, 12])
     !nfft = product(ngfft(1:3))
     !call dvdb_open_read(db, ngfft, xmpi_comm_self)

     call dvdb_free(db)

   case ("test_symmetries")
     call wrtout(std_out," Testing symmetries",'COLL')
     call get_command_argument(2, db_path)

     call dvdb_init(db, db_path, comm)

     call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
     nfft = product(ngfft(1:3))
     call dvdb_open_read(db, ngfft, xmpi_comm_self)

     call dvdb_test_symmetries(db)
     !call dvdb_check_v1sym(db)

     call dvdb_free(db)

   case ("test_ftinterp")
     call get_command_argument(2, db_path)
     if (nargs >= 3) then
       call get_command_argument(3, arg)
       !call replace_char(arg, ",", " ")
       !read(arg,"(3(i0,a))")ngqpt
     else
       ngqpt = [2,2,2]
       ngqpt = [4,4,4]
       !ngqpt = [8,8,8]
     end if

     call wrtout(std_out," Calling dvdb_ftinterp_setup",'COLL')
     call dvdb_init(db, db_path, comm)

     call ngfft_seq(ngfft, db%ngfft3_v1(:,1))
     nfft = product(ngfft(1:3))
     call dvdb_open_read(db, ngfft, xmpi_comm_self)

     ABI_MALLOC(intp_v1r, (2,nfft,db%nspden,db%natom3))
     ABI_MALLOC(file_v1r, (2,nfft,db%nspden,db%natom3))
     
     call dvdb_ftinterp_setup(db,ngqpt,1,[zero,zero,zero],nfft,ngfft,comm)

     do iq=1,db%nqpt
       ! Read data from file
       call dvdb_readsym_allv1(db, dvdb_findq(db, db%qpts(:,iq)), cplex, nfft, ngfft, tmp_v1r, comm)
       if (cplex == 1) then
         file_v1r(1,:,:,:) = tmp_v1r(1,:,:,:)
         file_v1r(2,:,:,:) = zero
       else
         file_v1r = tmp_v1r
       end if
       ABI_FREE(tmp_v1r)

       ! Interpolate data on the same q-point
       call dvdb_ftinterp_qpt(db, db%qpts(:,iq), nfft, ngfft, intp_v1r, comm)

       write(std_out,*)sjoin("For q-point:", ktoa(db%qpts(:,iq)))
       do mu=1,db%natom3
         do ispden=1,db%nspden
           call vdiff_print(vdiff_eval(2,nfft,file_v1r(:,:,ispden,mu),intp_v1r(:,:,ispden,mu),db%cryst%ucvol))

           do ifft=1,2
           !do ifft=1,nfft
             write(std_out,*)file_v1r(1,ifft,ispden,mu),intp_v1r(1,ifft,ispden,mu),&
             file_v1r(2,ifft,ispden,mu),intp_v1r(2,ifft,ispden,mu) 
           end do
         end do
       end do
       write(std_out,*)""
     end do ! iq

     call dvdb_free(db)
     ABI_FREE(intp_v1r)
     ABI_FREE(file_v1r)

   case default
     MSG_ERROR(sjoin("Unknown command:", command))
   end select

 end if

 call wrtout(std_out," Done",'COLL')

 call abinit_doctor("__mrgdv")

 100 call xmpi_end()

 end program mrgdv
!!***
