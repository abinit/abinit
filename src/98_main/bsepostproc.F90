!{\src2tex{textfont=tt}}
!!****p* ABINIT/bsepostproc
!! NAME
!! bsepostproc
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2013-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
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

program bsepostproc

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_errors
 use m_hdr
 use m_ebands
 use m_crystal
 use m_kpts
 use m_tetrahedron


 use defs_datatypes,   only : ebands_t
 use m_fstrings,       only : sjoin, strcat, basename
 use m_specialmsg,     only : herald
 use m_argparse,       only : get_arg, get_arg_list

 implicit none

!Arguments ----------------------------
!Local variables-----------------------
!scalars
 integer,parameter :: master = 0
 integer :: ii, nargs, comm, my_rank, nprocs, prtvol, fform, rdwr, prtebands
 integer :: kptopt, nshiftk, new_nshiftk, chksymbreak, nkibz, nkbz, ierr, occopt, intmeth
 real(dp) :: spinmagntarget, tsmear, extrael, step, broad
 character(len=500) :: command, arg, msg
 character(len=fnlen) :: path !, prefix
 type(hdr_type) :: hdr
 type(ebands_t) :: ebands !, ebands_kpath
 type(edos_t) :: edos
 type(crystal_t) :: cryst
 type(t_tetrahedron) :: tetra
!arrays
 integer :: kptrlatt(3,3), new_kptrlatt(3,3)
 !real(dp):: params(4)
 real(dp),allocatable :: shiftk(:,:), new_shiftk(:,:), wtk(:), kibz(:,:), kbz(:,:)

!*******************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 nargs = command_argument_count()
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)

 !if (nargs == 0) then
 !else
 ! Command line options.
 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     write(std_out,"(a)")" --version              Show version number and exit."
     write(std_out,"(a)")" -h, --help                 Show this help and exit."
     write(std_out,"(a)")" -v                         Increase verbosity level"
     write(std_out,"(2a)")ch10,"=== HEADER ==="
     write(std_out,"(a)")"hdr FILE                   Print ABINIT header."
     write(std_out,"(2a)")ch10,"=== KPOINTS ==="
     write(std_out,"(a)")"ibz FILE --kptopt 1 --kptrlatt or --ngkpt --shiftk 0.5 0.5, 0.5 --chksymbreak 1"
     write(std_out,"(2a)")ch10,"=== CRYSTAL ==="
     write(std_out,"(a)")"crystal_print FILE                   Print info on crystalline structure."
     write(std_out,"(2a)")ch10,"=== ELECTRONS ==="
     write(std_out,"(a)")"ebands_print FILE                    Print info on electron band structure."
     write(std_out,"(a)")"ebands_xmgrace FILE                  Produce XMGRACE file with bands."
     write(std_out,"(a)")"ebands_gnuplot FILE                  Produce GNUPLOT file with bands."
     write(std_out,"(a)")"ebands_dos FILE --intmeth, --step, --broad  Compute electron DOS."
     !write(std_out,"(a)")"ebands_jdos FILE --intmeth, --step, --broad  Compute electron DOS."
     write(std_out,"(a)")"ebands_bxsf FILE                     Produce BXSF file for Xcrysden."
     !write(std_out,"(a)")"ebands_skw_path FILE                     Produce BXSF file for Xcrysden."
     write(std_out,"(a)")"ebands_extrael FILE --occopt --tsmear --extrael  Change number of electron, compute new Fermi level."
     write(std_out,"(2a)")ch10,"=== DEVELOPERS ==="
     write(std_out,"(a)")"test_mjv                             Old tetrahedron routine"
     write(std_out,"(a)")"test_unit_tests                      Run unit tests for tetrahedron routines."
     goto 100
   end if
 end do

 call get_command_argument(1, command)

 select case (command)
 case ("hdr_print")
   ABI_CHECK(nargs > 1, "FILE argument is required.")
   call get_command_argument(2, path)
   call hdr_read_from_fname(hdr, path, fform, comm)
   ABI_CHECK(fform /= 0, "fform == 0")
   !rdwr = 3; if (prtvol > 0) rdwr = 4
   rdwr = 4
   call hdr_echo(hdr, fform, rdwr, unit=std_out)
   
 case ("ibz")
   ! Print list of kpoints in the IBZ with the corresponding weights
   call get_path_cryst(path, cryst, comm)
   call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   ABI_CHECK(any(kptrlatt /= 0), "kptrlatt or ngkpt must be specified")

   call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
      new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk) !, bz2ibz)  ! Optional
   new_nshiftk = size(new_shiftk, dim=2)

   write(std_out, "(/, a)")" Input_kptrlatt | New_kptrlatt"
   do ii=1,3
     write(std_out, "(2(a, 3(i0,1x)), a)")" [", kptrlatt(ii, :), "],       [ ", new_kptrlatt(ii, :), "]"
   end do
   write(std_out, "(/, a)")" Input_shiftk   |   New_shiftk"
   do ii=1,max(nshiftk, new_nshiftk)
     if (ii <= nshiftk .and. ii <= new_nshiftk) then
       write(std_out, "(2(a, 3(f3.1, 1x)), a)")" [", shiftk(:, ii), "],   [", new_shiftk(:, ii), "]"
     else if (ii <= nshiftk .and. ii > new_nshiftk) then
       write(std_out, "((a, 3(f3.1, 1x)), a)")" [", shiftk(:, ii), "]      ......"
     end if
   end do
   write(std_out, "(/, a, i0)")" nkibz: ", nkibz

 !case ("testkgrid")
   !call get_path_cryst(path, cryst, comm)
   !call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   !call testkgrid(bravais, iout, kptrlatt, kptrlen, msym, nshiftk, nsym, prtkpt, rprimd, shiftk, symafm, symrel, vacuum)

 case ("crystal_print")
    call get_path_cryst(path, cryst, comm)
    call crystal_print(cryst, unit=std_out, prtvol=prtvol)
 
 case ("ebands_print", "ebands_xmgrace", "ebands_gnuplot")
   call get_path_ebands(path, ebands, comm)
   if (command == "ebands_print") then
     call ebands_print(ebands, unit=std_out, prtvol=prtvol)
   else
     prtebands = 1; if (command == "ebands_gnuplot") prtebands = 2
     call ebands_write(ebands, prtebands, basename(path))
   end if

 case ("ebands_dos", "ebands_jdos")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   ABI_CHECK(get_arg("intmeth", intmeth, msg, default=1) == 0, msg)
   ABI_CHECK(get_arg("step", step, msg, default=0.02 * eV_Ha) == 0, msg)
   ABI_CHECK(get_arg("broad", broad, msg, default=0.04 * ev_Ha) == 0, msg)

   if (command == "ebands_dos") then
     edos = ebands_get_edos(ebands, cryst, intmeth, step, broad, comm)
     call edos%write(strcat(basename(path), "_EDOS"))
     call edos%free()
   else if (command == "ebands_jdos") then
     !jdos = ebands_get_jdos(ebands, cryst, intmeth, step, broad, comm, ierr)
     !call jdos%write(strcat(basename(path), "_EJDOS"))
     !call jdos%free()
   end if

 case ("ebands_bxsf")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   if (ebands_write_bxsf(ebands, cryst, strcat(basename(path), "_BXSF")) /= 0)  then
     MSG_ERROR("Cannot produce file for Fermi surface, check log file for more info")
   end if

 !case ("nesting")
   !call get_path_ebands_cryst(path, ebands, cryst, comm)
   !ierr = ebands_write_nesting(ebands,cryst,filepath,prtnest,tsmear,fermie_nest, qpath_vertices,errmsg)

 case ("ebands_skw_kpath")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   ! Generate k-path
   !ABI_CHECK(get_arg("ndivsm", ndivsm, msg, default=20) == 0, msg)
   !nbounds = dtset%nkpath
   !if (nbounds <= 0) then
   !  MSG_COMMENT("Using hard-coded k-path because nkpath not present in input file.")
   !  nbounds = 5
   !  ABI_MALLOC(bounds, (3, 5))
   !  bounds = reshape([zero, zero, zero, half, zero, zero, zero, half, zero, zero, zero, zero, zero, zero, half], [3,5])
   !else
   !  call alloc_copy(dtset%kptbounds, bounds)
   !end if
   !kpath = kpath_new(bounds, cryst%gprimd, ndivsm)
   !call kpath_print(kpath, header="Interpolating energies on k-path", unit=std_out)
   ! Interpolate band energies with star-functions
   !params = 0; params(1) = 1; params(2) = 5
   !if (nint(dtset%einterp(1)) == 1) params = dtset%einterp
   !ebands_kpath = ebands_interp_kpath(ebands, cryst, kpath, params, [1, ebands%mband], comm)
   !call kpath%free()
   !call ebands_free(ebands_kpath)
   !ABI_CHECK(get_arg("prtebands", prtebands, msg, default=1) == 0, msg)
   !call ebands_write(ebands_kpath, prtebands, basename(path))
   !call ebands_kpath(ebands_kpath)

   !ebands_kmesh = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
   !     band_block, comm, out_prefix)
   !call ebands_free(ebands_kmesh)

 case ("ebands_extrael")
   call get_path_ebands(path, ebands, comm)
   ABI_CHECK(get_arg("occopt", occopt, msg, default=3) == 0, msg)
   ABI_CHECK(get_arg("tsmear", tsmear, msg, default=tol2) == 0, msg)
   ABI_CHECK(get_arg("extrael", extrael, msg) == 0, msg)
   ABI_CHECK(get_arg("spinmagntarget", spinmagntarget, msg, default=-99.99_dp) == 0, msg)

   call ebands_set_scheme(ebands, occopt, tsmear, spinmagntarget, prtvol)
   call ebands_set_nelect(ebands, ebands%nelect + extrael, spinmagntarget, msg)
   write(std_out, "(a)") msg
   call ebands_update_occ(ebands, spinmagntarget, prtvol=prtvol)
   call ebands_print(ebands, prtvol=prtvol)

 case ("tetra_mjv")
   call get_path_cryst(path, cryst, comm)
   call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   ABI_CHECK(any(kptrlatt /= 0), "ngkpt or kptrlatt must be specified")

   call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
      new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk) !, bz2ibz)  ! Optional

   new_nshiftk = size(new_shiftk, dim=2)
   tetra = tetra_from_kptrlatt(cryst, kptopt, new_kptrlatt, new_nshiftk, new_shiftk, nkibz, kibz, comm, msg, ierr)
   ABI_CHECK(ierr == 0, msg)
   call tetra_write(tetra, nkibz, kibz, strcat(basename(path), "_TETRA"))
   call destroy_tetra(tetra)

 case ("tetra_unit_tests")
   !call phdos_unittests(comm)

 case default
   MSG_ERROR(sjoin("Unknown command:", command))
 end select

 ! Deallocate memory to make memcheck happy.
 call hdr_free(hdr)
 call cryst%free()
 call ebands_free(ebands)

 ABI_SFREE(kibz)
 ABI_SFREE(wtk)
 ABI_SFREE(kbz)
 ABI_SFREE(new_shiftk)

 call abinit_doctor("__bsepostproc")

 100 call xmpi_end()

contains 
!!***

!!****f* ABINIT/get_ebands_path_cryst
!! NAME
!! get_ebands_path_cryst
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


subroutine get_path_ebands_cryst(path, ebands, cryst, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(out) :: path
 type(ebands_t),intent(out) :: ebands
 type(crystal_t),intent(out) :: cryst
 integer,intent(in) :: comm

!Arguments ----------------------------
!Local variables-----------------------
 integer :: nargs

 nargs = command_argument_count()
 ABI_CHECK(nargs > 1, "FILE argument is required.")
 call get_command_argument(2, path)
 cryst = crystal_from_file(path, comm)
 ebands = ebands_from_file(path, comm)

end subroutine get_path_ebands_cryst
!!***

!!****f* ABINIT/get_path_ebands
!! NAME
!! get_path_ebands
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_path_ebands(path, ebands, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(out) :: path
 type(ebands_t),intent(out) :: ebands
 integer,intent(in) :: comm

!Arguments ----------------------------
!Local variables-----------------------
 integer :: nargs

 nargs = command_argument_count()
 ABI_CHECK(nargs > 1, "FILE argument is required.")
 call get_command_argument(2, path)
 ebands = ebands_from_file(path, comm)

end subroutine get_path_ebands
!!***

!!****f* ABINIT/get_path_cryst
!! NAME
!! get_path_cryst
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_path_cryst(path, cryst, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(out) :: path
 type(crystal_t),intent(out) :: cryst
 integer,intent(in) :: comm

!Arguments ----------------------------
!Local variables-----------------------
 integer :: nargs

 nargs = command_argument_count()
 ABI_CHECK(nargs > 1, "FILE argument is required.")
 call get_command_argument(2, path)
 cryst = crystal_from_file(path, comm)

end subroutine get_path_cryst
!!***

!!****f* ABINIT/parse_kargs
!! NAME
!! parse_kargs
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)

 use m_fstrings,       only : atoi

!Arguments ------------------------------------
 integer,intent(out) :: kptopt, nshiftk, chksymbreak
 integer,intent(out) :: kptrlatt(3,3)
 real(dp),allocatable,intent(out) :: shiftk(:,:)

!Local variables-------------------------------
 integer :: ii, lenr
 character(len=500) :: msg
 integer :: ivec9(9), ngkpt(3)
 real(dp) :: my_shiftk(3 * MAX_NSHIFTK)

! *************************************************************************
 
 ABI_CHECK(get_arg("kptopt", kptopt, msg, default=1) == 0, msg)
 ABI_CHECK(get_arg("chksymbreak", chksymbreak, msg, default=1) == 0, msg)
 ABI_CHECK(get_arg_list("ngkpt", ngkpt, lenr, msg, exclude="kptrlatt", want_len=3) == 0, msg)
 if (lenr == 3) then
   kptrlatt = 0
   do ii=1,3
     kptrlatt(ii, ii) = ngkpt(ii)
   end do
 end if
 ABI_CHECK(get_arg_list("kptrlatt", ivec9, lenr, msg, exclude="ngkpt", want_len=9) == 0, msg)
 if (lenr == 9) kptrlatt = transpose(reshape(ivec9, [3, 3]))

 ! Init default
 ABI_CHECK(get_arg_list("shiftk", my_shiftk, lenr, msg) == 0, msg)
 if (lenr /= 0) then
   ABI_CHECK(mod(lenr, 3) == 0, "Expecting 3 * nshift array")
   nshiftk = lenr / 3
   ABI_MALLOC(shiftk, (3, nshiftk))
   shiftk = reshape(my_shiftk(1:lenr), [3, nshiftk])
 else
   nshiftk = 1
   ABI_MALLOC(shiftk, (3, nshiftk))
   shiftk(:, 1) = [half, half, half]
 end if

 !write(std_out, *)"kptopt = ", kptopt, ", chksymbreak = ", chksymbreak, ", nshiftk = ", nshiftk, ", kptrlatt = ", kptrlatt

end subroutine parse_kargs
!!***

!!****f* ABINIT/crystal_from_file
!! NAME
!! crystal_from_file
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


type(crystal_t) function crystal_from_file(path, comm) result(new)

 use m_fstrings,       only : endswith

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer :: fform, timrev
 type(hdr_type) :: hdr

! *************************************************************************

 if (endswith(path, ".abi") .or. endswith(path, ".in")) then
   NOT_IMPLEMENTED_ERROR()
   !new = crystal_from_abinit_input(path, comm)
   !call iofn1(filnam, filstat, comm)
   !call get_dtsets_pspheads(path, ndtset, lenstr, string, timopt, dtsets, pspheads, mx, dmatpuflag, comm)
 else
 !if (endswith(path, ".nc")) then
    !new = crystal_from_ncfile(path, comm)
    call hdr_read_from_fname(hdr, path, fform, comm)
    ABI_CHECK(fform /= 0, "fform == 0")
    timrev = 2 !; (if kpts_timrev_from_kptopt(hdr%kptopt) == 0) timrev = 1
    new = hdr_get_crystal(hdr, timrev)
    call hdr_free(hdr)
 !else
 !  MSG_ERROR(sjoin("Don't know how to construct crystal structure from:", path, ch10, "Supported extensions: .nc, .abi., .in"))
 end if

end function crystal_from_file
!!***

!!****f* ABINIT/ebands_from_file
!! NAME
!! ebands_from_file
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


type(ebands_t) function ebands_from_file(path, comm) result(new)

 use m_fstrings,       only : endswith

 use m_wfk

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer :: fform
 type(hdr_type) :: hdr
!arrays
 real(dp),pointer :: gs_eigen(:,:,:)

! *************************************************************************

 if (endswith(path, "_WFK") .or. endswith(path, "_WFK.nc")) then
   call hdr_read_from_fname(hdr, path, fform, comm)
   ABI_CHECK(fform /= 0, "fform == 0")
   call wfk_read_eigenvalues(path, gs_eigen, hdr, comm)
   new = ebands_from_hdr(hdr, maxval(hdr%nband), gs_eigen)
   call hdr_free(hdr)
   ABI_FREE(gs_eigen)
 else if (endswith(path, ".nc")) then
   NOT_IMPLEMENTED_ERROR()
   !new = ebands_ncread(path, comm)
 else
   MSG_ERROR(sjoin("Don't know how to construct crystal structure from: ", path, ch10, "Supported extensions: _WFK or .nc"))
 end if

end function ebands_from_file
!!***

end program bsepostproc
!!***
