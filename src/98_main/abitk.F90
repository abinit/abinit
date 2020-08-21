!!****p* ABINIT/abitk
!! NAME
!! abitk
!!
!! FUNCTION
!!  Command line interface to perform very basic post-processing of output files (mainly netcdf files).
!!  Use `abitk --help` to get list of possible commands.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (main routine)
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

program abitk

 use defs_basis
 use m_abicore
 use m_build_info
 use m_xmpi
 use m_errors
 use m_hdr
 use m_ebands
 use m_crystal
 use m_kpts
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk

 use defs_datatypes,   only : ebands_t
 use m_fstrings,       only : sjoin, strcat, basename
 use m_specialmsg,     only : herald
 use m_symtk,          only : matr3inv
 use m_bz_mesh,        only : kpath_new, kpath_t
 use m_unittests,      only : tetra_unittests, kptrank_unittests
 use m_argparse,       only : get_arg, get_arg_list, parse_kargs
 use m_common,         only : ebands_from_file, crystal_from_file
 use m_parser,         only : geo_t, geo_from_poscar_path

 implicit none

!Arguments ----------------------------
!Local variables-----------------------
!scalars
 integer,parameter :: master = 0
 integer :: ii, nargs, comm, my_rank, nprocs, prtvol, fform, rdwr, prtebands
 integer :: kptopt, nshiftk, new_nshiftk, chksymbreak, nkibz, nkbz, occopt, intmeth, lenr
 integer :: ndivsm, abimem_level, ierr !, spin
 real(dp) :: spinmagntarget, tsmear, extrael, step, broad, abimem_limit_mb !, tolsym
 character(len=500) :: command, arg, msg, ptgroup
 character(len=fnlen) :: path, other_path !, prefix
 type(hdr_type) :: hdr
 type(ebands_t) :: ebands, ebands_kpath, other_ebands
 type(edos_t) :: edos
 type(jdos_t) :: jdos
 type(crystal_t) :: cryst, other_cryst
 type(geo_t) :: geo
 type(kpath_t) :: kpath
!arrays
 integer :: kptrlatt(3,3), new_kptrlatt(3,3), ngqpt(3)
 real(dp) :: skw_params(4)
 real(dp),allocatable :: bounds(:,:)
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
 ABI_CHECK(get_arg("abimem-level", abimem_level, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("abimem-limit-mb", abimem_limit_mb, msg, default=20.0_dp) == 0, msg)
#ifdef HAVE_MEM_PROFILING
 call abimem_init(abimem_level, limit_mb=abimem_limit_mb)
#endif

 nargs = command_argument_count()
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)

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
     write(std_out,"(a)")"ibz FILE --ngkpt 2 2 2 or --kptrlatt [--kptopt 1] [--shiftk 0.5 0.5 0.5] [--chksymbreak 1]"

     write(std_out,"(2a)")ch10,"=== CRYSTAL ==="
     write(std_out,"(a)")"crystal_print FILE                   Print info on crystalline structure."
     write(std_out,"(a)")"from_poscar POSCAR_FILE              Read POSCAR file, print abinit variables."

     write(std_out,"(2a)")ch10,"=== ELECTRONS ==="
     write(std_out,"(a)")"ebands_print FILE                    Print info on electron band structure."
     write(std_out,"(a)")"ebands_xmgrace FILE                  Produce XMGRACE file with electron bands."
     write(std_out,"(a)")"ebands_gnuplot FILE                  Produce GNUPLOT file with electron bands."
     write(std_out,"(a)")"ebands_dos FILE --intmeth, --step, --broad  Compute electron DOS."
     write(std_out,"(a)")"ebands_bxsf FILE                     Produce BXSF file for Xcrysden."
     write(std_out,"(a)")"ebands_extrael FILE --occopt --tsmear --extrael  Change number of electron, compute new Fermi level."
     write(std_out,"(a)")"ebands_gaps FILE                     Print info on gaps"
     !write(std_out,"(a)")"ebands_jdos FILE --intmeth, --step, --broad  Compute electron DOS."
     !write(std_out,"(a)")"skw_path FILE                       Interpolate band structure along a k-path."
     write(std_out,"(a)")"skw_compare IBZ_WFK KPATH_WFK        Use e_nk from IBZ_WFK to interpolate on the k-path in KPATH_WFK."

     write(std_out,"(2a)")ch10,"=== DEVELOPERS ==="
     write(std_out,"(a)")"tetra_unit_tests                      Run unit tests for tetrahedron routines."
     write(std_out,"(a)")"kptrank_unit_tests                    Run unit tests for kptrank routines."
     goto 100
   end if
 end do

 call get_command_argument(1, command)

 select case (command)
 case ("from_poscar")
   call get_command_argument(2, path)
   geo = geo_from_poscar_path(path, comm)
   call geo%print_abivars(std_out)
   call geo%free()

 !case ("to_poscar")
 !  call get_command_argument(2, path)
 !  call get_path_cryst(path, cryst, comm)
 !  call prtposcar(fcart, fnameradix, natom, ntypat, rprimd, typat, ucvol, xred, znucl)

 !case ("to_cif")
 !  call get_command_argument(2, path)
 !  call get_path_cryst(path, cryst, comm)
 !  call prt_cif(brvltt, ciffname, natom, nsym, ntypat, rprimd, spgaxor, spgroup, spgorig, symrel, tnon, typat, xred, znucl)

 case ("hdr_print")
   ABI_CHECK(nargs > 1, "FILE argument is required.")
   call get_command_argument(2, path)
   call hdr_read_from_fname(hdr, path, fform, comm)
   ABI_CHECK(fform /= 0, "fform == 0")
   rdwr = 3; if (prtvol > 0) rdwr = 4
   !rdwr = 4
   call hdr%echo(fform, rdwr, unit=std_out)

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
   write(std_out,"(/,a)")"# List of kpoints in the IBZ with the corresponding weights:"
   write(std_out,"(a,i0)")" kptopt ", kptopt
   write(std_out,"(a,i0,/,a)")" nkpt ", nkibz, " kpt"
   do ii=1,nkibz
     write(std_out, "(3(es11.4),a,es11.4)") kibz(:, ii), " # wtk ", wtk(ii)
   end do

 !case ("testkgrid")
   !call get_path_cryst(path, cryst, comm)
   !call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   !call testkgrid(bravais, iout, kptrlatt, kptrlen, msym, nshiftk, nsym, prtkpt, rprimd, shiftk, symafm, symrel, vacuum)

 case ("crystal_print")
    call get_path_cryst(path, cryst, comm)
    call cryst%print(unit=std_out, prtvol=prtvol)

 case ("ebands_print", "ebands_xmgrace", "ebands_gnuplot")
   call get_path_ebands(path, ebands, comm)
   if (command == "ebands_print") then
     call ebands_print(ebands, unit=std_out, prtvol=prtvol)
   else
     prtebands = 1; if (command == "ebands_gnuplot") prtebands = 2
     call ebands_write(ebands, prtebands, basename(path))
   end if

 case ("ebands_gaps")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   call ebands_print_gaps(ebands, std_out)

 case ("ebands_dos", "ebands_jdos")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   ABI_CHECK(get_arg("intmeth", intmeth, msg, default=2) == 0, msg)
   ABI_CHECK(get_arg("step", step, msg, default=0.02 * eV_Ha) == 0, msg)
   ABI_CHECK(get_arg("broad", broad, msg, default=0.04 * ev_Ha) == 0, msg)

   if (command == "ebands_dos") then
     edos = ebands_get_edos(ebands, cryst, intmeth, step, broad, comm)
     call edos%write(strcat(basename(path), "_EDOS"))
     call edos%free()

   else if (command == "ebands_jdos") then
     NOT_IMPLEMENTED_ERROR()
     jdos = ebands_get_jdos(ebands, cryst, intmeth, step, broad, comm, ierr)
     !call jdos%write(strcat(basename(path), "_EJDOS"))
     call jdos%free()
   end if

 case ("ebands_bxsf")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   if (ebands_write_bxsf(ebands, cryst, strcat(basename(path), "_BXSF")) /= 0)  then
     MSG_ERROR("Cannot produce file for Fermi surface in BXSF format. Check log file for info.")
   end if

 !case ("ebands_nesting")
   !call get_path_ebands_cryst(path, ebands, cryst, comm)
   !if (ebands_write_nesting(ebands, cryst, filepath, prtnest, tsmear, fermie_nest, qpath_vertices, errmsg) /= 0) then
   !  MSG_ERROR("Cannot produce file for Fermi surface in BXSF format. Check log file for info.")
   !end if

 case ("skw_kpath")
   ! Get energies on the IBZ from path file
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ! Generate k-path
   ABI_CHECK(get_arg("ndivsm", ndivsm, msg, default=20) == 0, msg)
   !MSG_COMMENT("Using hard-coded k-path because nkpath not present in input file.")
   !nbounds = 5
   ABI_MALLOC(bounds, (3, 5))
   bounds = reshape([zero, zero, zero, half, zero, zero, zero, half, zero, zero, zero, zero, zero, zero, half], [3,5])
   !ABI_CHECK(get_args_from_string("bounds", bounds, msg, default="[0,0,0,1/2,0,0,0,1/2]") == 0, msg)
   kpath = kpath_new(bounds, cryst%gprimd, ndivsm)
   ABI_FREE(bounds)

   call kpath%print(header="Interpolating energies on k-path", unit=std_out)

   ! Interpolate band energies with star-functions
   call parse_skw_params(skw_params)
   !ABI_CHECK(get_arg_list("einterp", ivec9, lenr, msg, exclude="ngkpt", want_len=9) == 0, msg)
   !if (nint(dtset%einterp(1)) == 1) skw_params = dtset%einterp
   ebands_kpath = ebands_interp_kpath(ebands, cryst, kpath, skw_params, [1, ebands%mband], comm)

   !call wrtout(std_out, sjoin(" Writing interpolated bands to:",  path)
   ABI_CHECK(get_arg("prtebands", prtebands, msg, default=2) == 0, msg)
   call ebands_write(ebands_kpath, prtebands, path)

   !ebands_kmesh = ebands_interp_kmesh(ebands, cryst, skw_params, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
   !     band_block, comm, out_prefix)
   !call ebands_free(ebands_kmesh)

 case ("skw_compare")
   ! Get energies on the IBZ from filepath
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ! Get ab-initio energies for the second file (assume k-path!)
   call get_path_ebands_cryst(other_path, other_ebands, other_cryst, comm, argpos=3)

   ! Interpolate band energies on the path with star-functions
   kpath = kpath_new(other_ebands%kptns, other_cryst%gprimd, -1)

   call parse_skw_params(skw_params)
   ebands_kpath = ebands_interp_kpath(ebands, cryst, kpath, skw_params, [1, ebands%mband], comm)

   ! Compare gaps
   call ebands_print_gaps(other_ebands, std_out, header="Ab-initio gaps")
   call ebands_print_gaps(ebands_kpath, std_out, header="SKW interpolated gaps")

   !ABI_CHECK(get_arg("prtebands", prtebands, msg, default=2) == 0, msg)
   !call ebands_write(ebands_kpath, prtebands, path)

   ! Write EBANDS file
#ifdef HAVE_NETCDF
   NCF_CHECK(ebands_ncwrite_path(other_ebands, cryst, "abinitio_EBANDS.nc"))
   NCF_CHECK(ebands_ncwrite_path(ebands_kpath, other_cryst, "skw_EBANDS.nc"))
#endif

   call wrtout(std_out, &
     ch10//" Use `abicomp.py ebands abinitio_EBANDS.nc skw_EBANDS.nc -p combiplot` to compare the bands with AbiPy.", &
     newlines=2)

 case ("ebands_extrael")

   ! Get energies on the IBZ from filepath
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ABI_CHECK(get_arg("extrael", extrael, msg) == 0, msg)
   ABI_CHECK(get_arg("occopt", occopt, msg, default=3) == 0, msg)
   ABI_CHECK(get_arg("tsmear", tsmear, msg, default=tol2) == 0, msg)
   ABI_CHECK(get_arg("spinmagntarget", spinmagntarget, msg, default=-99.99_dp) == 0, msg)

   !extrael = - dprarr(1) * cryst%ucvol * (Bohr_meter * 100) ** 3
   !ebands%nelect = ebands%nelect + extrael
   !call ebands_print_gaps(ebands, std_out, header="KS gaps")
   !integer,intent(in) :: ntemp
   !real(dp),intent(in) :: kTmesh(ntemp)
   !real(dp),intent(out) :: mu_e(ntemp)
   !call ebands_get_muT_with_fd(ebands, ntemp, kTmesh, spinmagntarget, prtvol, mu_e, comm)

   call ebands_set_scheme(ebands, occopt, tsmear, spinmagntarget, prtvol)
   call ebands_set_nelect(ebands, ebands%nelect + extrael, spinmagntarget, msg)
   write(std_out, "(a)") msg
   call ebands_update_occ(ebands, spinmagntarget, prtvol=prtvol)
   call ebands_print(ebands, prtvol=prtvol)

 ! ====================
 ! Tools for developers
 ! ====================
 !case ("denpot_f2nc")
 !   call get_command_argument(2, path)
 !   call get_command_argument(3, other_path)
 !   call denpot_f2nc(path, other_path)

   ! Get important dimensions from the first header and rewind the file.
   !call hdr_fort_read(new%hdr_ref, unt, fform)
   !if (dvdb_check_fform(fform, "read_dvdb", msg) /= 0) then
   !  MSG_ERROR(sjoin("While reading:", path, ch10, msg))
   !end if
   !! Fortran IO
   !do ispden=1,hdr1%nspden
   !  read(units(jj), err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
   !  write(ount, err=10, iomsg=msg) (v1(ifft), ifft=1,cplex*nfft)
   !end do
   !if (new%debug) call new%hdr_ref%echo(fform, 4, unit=std_out)
   !call hdr1%free()

   !! first order potentials are always written because the eph code requires them
   !! the files are small (much much smaller that 1WFK, actually we should avoid writing 1WFK)
   !rdwrpaw=0
   !call appdig(pertcase,dtfil%fnameabo_pot,fi1o)
   !! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
   !call fftdatar_write_from_hdr("first_order_potential",fi1o,dtset%iomode,hdr,&
   !ngfftf,cplex,nfftf,dtset%nspden,vtrial1,mpi_enreg)

 ! ===========
 ! Unit tests
 ! ===========

 case ("tetra_unit_tests")
   ABI_CHECK(get_arg("ptgroup", ptgroup, msg, default="m-3m") == 0, msg)
   ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default_list=[20, 20, 20]) == 0, msg)
   call tetra_unittests(ptgroup, ngqpt, comm)

 case ("kptrank_unit_tests")
   ABI_CHECK(get_arg("ptgroup", ptgroup, msg, default="m-3m") == 0, msg)
   ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default_list=[100, 100, 100]) == 0, msg)
   call kptrank_unittests(ptgroup, ngqpt, comm)

 case default
   MSG_ERROR(sjoin("Unknown command:", command))
 end select

 ! Deallocate memory to make memcheck happy.
 call hdr%free()
 call cryst%free()
 call other_cryst%free()
 call kpath%free()
 call ebands_free(ebands)
 call ebands_free(ebands_kpath)
 call ebands_free(other_ebands)

 ABI_SFREE(kibz)
 ABI_SFREE(wtk)
 ABI_SFREE(kbz)
 ABI_SFREE(new_shiftk)

 call abinit_doctor("__abitk")

 100 call xmpi_end()

contains
!!***

!!****f* abitk/get_path_ebands_cryst
!! NAME
!! get_path_ebands_cryst
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abitk
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_path_ebands_cryst(path, ebands, cryst, comm, argpos)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(out) :: path
 type(ebands_t),intent(out) :: ebands
 type(crystal_t),intent(out) :: cryst
 integer,intent(in) :: comm
 integer,intent(in),optional :: argpos

!Arguments ----------------------------
!Local variables-----------------------
 integer :: nargs, apos

 apos = 2; if (present(argpos)) apos = argpos
 nargs = command_argument_count()
 ABI_CHECK(nargs >= apos, "FILE argument is required.")
 call get_command_argument(apos, path)
 cryst = crystal_from_file(path, comm)
 ebands = ebands_from_file(path, comm)

end subroutine get_path_ebands_cryst
!!***

!!****f* abitk/get_path_ebands
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
!!      abitk
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

!!****f* abitk/get_path_cryst
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
!!      abitk
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

!!****f* abitk/parse_skw_params
!! NAME
!! parse_skw_params
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abitk
!!
!! CHILDREN
!!
!! SOURCE

subroutine parse_skw_params(params)

!Arguments ----------------------------
!Local variables-----------------------
 real(dp),intent(out) :: params(4)

 integer :: lpratio
 real(dp) :: rcut, rsigma

!*******************************************************

 ABI_CHECK(get_arg("lpratio", lpratio, msg, default=5) == 0, msg)
 ABI_CHECK(get_arg("rcut", rcut, msg, default=zero) == 0, msg)
 ABI_CHECK(get_arg("rsigma", rsigma, msg, default=zero) == 0, msg)
 params = [one, dble(lpratio), rcut, rsigma]

end subroutine parse_skw_params
!!***

end program abitk
!!***
