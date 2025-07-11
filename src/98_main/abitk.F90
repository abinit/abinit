!!****p* ABINIT/abitk
!! NAME
!! abitk
!!
!! FUNCTION
!!  Command line interface to perform very post-processing of output files (mainly netcdf files).
!!  Use `abitk --help` to get list of possible commands.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2025 ABINIT group (MG)
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program abitk

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hdr
 use m_ebands
 use m_crystal
 use m_kpts
 use netcdf
 use m_nctk

 use m_build_info,     only : abinit_version
 use m_fstrings,       only : sjoin, strcat, basename, itoa
 use m_io_tools,       only : open_file, enforce_fortran_io
 use m_specialmsg,     only : herald
 use m_matrix,         only : matr3inv
 use m_numeric_tools,  only : arth
 use m_bz_mesh,        only : kpath_new, kpath_t
 use m_unittests,      only : tetra_unittests, kptrank_unittests, tetra_zinv_convergence
 use m_argparse,       only : get_arg, get_arg_list, parse_kargs
 use m_common,         only : ebands_from_file, crystal_from_file
 use m_parser,         only : geo_t, geo_from_poscar_path
 use m_phgamma,        only : find_ewin

 implicit none

!Arguments ----------------------------
!Local variables-----------------------
!scalars
 integer,parameter :: master = 0
 integer :: ii, nargs, comm, my_rank, nprocs, prtvol, fform, rdwr, prtebands, spin, ncid, ncerr
 integer :: kptopt, nshiftk, new_nshiftk, chksymbreak, nkibz, nkbz, intmeth, lenr !occopt,
 integer :: ndivsm, abimem_level, ierr, ntemp, ios, itemp, use_symmetries, ltetra
 real(dp) :: spinmagntarget, extrael, doping, step, broad, abimem_limit_mb, fs_ewin !, tolsym, tsmear
 logical :: is_metal
 character(len=500) :: command, arg, msg, ptgroup
 character(len=fnlen) :: path, other_path, out_path !, prefix
 type(hdr_type) :: hdr
 type(ebands_t) :: ebands, ebands_kpath, other_ebands
 type(edos_t) :: edos
 type(gaps_t) :: gaps
 type(jdos_t) :: jdos
 type(crystal_t) :: cryst, other_cryst
 type(geo_t) :: geo
 type(kpath_t) :: kpath
!arrays
 integer :: kptrlatt(3,3), new_kptrlatt(3,3), ngqpt(3)
 integer, allocatable :: bz2ibz(:,:)
 real(dp) :: skw_params(4), tmesh(3)
 real(dp),allocatable :: bounds(:,:), kTmesh(:), mu_e(:), n_ehst(:,:,:)
 real(dp),allocatable :: shiftk(:,:), new_shiftk(:,:), wtk(:), kibz(:,:), kbz(:,:)
!*******************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 if (my_rank == master) then
#ifdef FC_NAG
   open(unit=ab_out,file="abitk.out",form='formatted',status='unknown', action="write", recl=ABI_RECL, iomsg=msg, iostat=ios)
#else
   open(unit=ab_out,file="abitk.out",form='formatted',status='unknown', action="write", iomsg=msg, iostat=ios)
#endif
   ABI_CHECK(ios == 0, msg)
   rewind (unit=ab_out)
 else
   close(std_out)
   if (open_file(NULL_FILE, msg, unit=std_out, action="write") /= 0) then
     ABI_ERROR(msg)
   end if
 end if

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

   else if (arg == "--enforce-fortran-io") then
      call enforce_fortran_io(.True.)

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     call abitk_show_help()
     goto 100
   end if
 end do

 call get_command_argument(1, command)
 call get_command_argument(2, path)

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
   call hdr%from_fname(path, fform, comm)
   ABI_CHECK(fform /= 0, "fform == 0")
   rdwr = 3; if (prtvol > 0) rdwr = 4
   call hdr%echo(fform, rdwr, unit=std_out)

 case ("ibz")
   ! Print list of k-points in the IBZ with the corresponding weights
   call get_path_cryst(path, cryst, comm)

   call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   ABI_CHECK(any(kptrlatt /= 0), "kptrlatt or ngkpt must be specified")

   call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, & ! in
                               nkibz, kibz, wtk, nkbz, kbz,              & ! out
                               new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk, bz2ibz=bz2ibz) ! out
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

   NCF_CHECK(nctk_open_create(ncid, "KMESH.nc", comm))
   NCF_CHECK(cryst%ncwrite(ncid))

   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nkibz", nkibz), nctkdim_t("nkbz", nkbz), &
     nctkdim_t("nshiftk", nshiftk), nctkdim_t("new_nshiftk", new_nshiftk) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "kptopt" &
   ])
   NCF_CHECK(ncerr)

   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "wr_step", &
   !])
   !NCF_CHECK(ncerr)

   ! Define arrays with results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("kptrlatt", "int", "three, three"), &
     nctkarr_t("new_kptrlatt", "int", "three, three"), &
     nctkarr_t("shiftk", "dp", "three, nshiftk"), &
     nctkarr_t("new_shiftk", "dp", "three, new_nshiftk"), &
     nctkarr_t("kibz", "dp", "three, nkibz"), &
     nctkarr_t("bz2ibz", "int", "six, nkbz"), &
     nctkarr_t("wtk", "dp", "nkibz"), &
     nctkarr_t("kbz", "dp", "three, nkbz") &
   ])
   NCF_CHECK(ncerr)

   ! ==========
   ! Write data
   ! ==========
   NCF_CHECK(nctk_set_datamode(ncid))
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "kptopt"], &
     [kptopt])
   NCF_CHECK(ncerr)

   !ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "wr_step", "ecuteps", "ecut", "ecutsigx", "gwr_boxcutmin", &
   !  [gwr%wr_step, dtset%ecuteps, dtset%ecut, dtset%ecutsigx, dtset%gwr_boxcutmin, &
   !  ])
   !NCF_CHECK(ncerr)

   NCF_CHECK(nf90_put_var(ncid, vid("kptrlatt"), kptrlatt))
   NCF_CHECK(nf90_put_var(ncid, vid("new_kptrlatt"), new_kptrlatt))
   NCF_CHECK(nf90_put_var(ncid, vid("shiftk"), shiftk))
   NCF_CHECK(nf90_put_var(ncid, vid("new_shiftk"), new_shiftk))
   NCF_CHECK(nf90_put_var(ncid, vid("kibz"), kibz))
   NCF_CHECK(nf90_put_var(ncid, vid("wtk"), wtk))
   NCF_CHECK(nf90_put_var(ncid, vid("kbz"), kbz))
   NCF_CHECK(nf90_put_var(ncid, vid("bz2ibz"), bz2ibz))
   NCF_CHECK(nf90_close(ncid))

   ABI_FREE(bz2ibz)

 !case ("testkgrid")
   !call get_path_cryst(path, cryst, comm)
   !call parse_kargs(kptopt, kptrlatt, nshiftk, shiftk, chksymbreak)
   !call testkgrid(bravais, iout, kptrlatt, kptrlen, msym, nshiftk, nsym, prtkpt, rprimd, shiftk, symafm, symrel, vacuum)

 case ("crystal_print")
    call get_path_cryst(path, cryst, comm)
    call cryst%print(unit=std_out, prtvol=prtvol)

 case ("crystal_abivars")
    call get_path_cryst(path, cryst, comm)
    call cryst%print_abivars(std_out)

 case ("ebands_print", "ebands_xmgrace", "ebands_gnuplot")
   call get_path_ebands(path, ebands, comm)
   if (command == "ebands_print") then
     call ebands%print([std_out], prtvol=prtvol)
   else
     prtebands = 1; if (command == "ebands_gnuplot") prtebands = 2
     call ebands%write(prtebands, basename(path))
   end if

 case ("ebands_gaps")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   call ebands%print_gaps([std_out])

 case ("ebands_edos", "ebands_jdos")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   ABI_CHECK(get_arg("intmeth", intmeth, msg, default=2) == 0, msg)
   ABI_CHECK(get_arg("step", step, msg, default=0.02 * eV_Ha) == 0, msg)
   ABI_CHECK(get_arg("broad", broad, msg, default=0.06 * eV_Ha) == 0, msg)

   if (command == "ebands_edos") then
     edos = ebands%get_edos(cryst, intmeth, step, broad, comm)
     call edos%print(std_out, header="Electron DOS")
     out_path = strcat(basename(path), "_EDOS")
     call wrtout(std_out, sjoin("Writing electron DOS to file:", out_path))
     call edos%write(out_path)

   else if (command == "ebands_jdos") then
     NOT_IMPLEMENTED_ERROR()
     jdos = ebands%get_jdos(cryst, intmeth, step, broad, comm, ierr)
     !call jdos%write(strcat(basename(path), "_EJDOS"))
   end if

 case ("ebands_bxsf")
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   if (ebands%write_bxsf(cryst, strcat(basename(path), "_BXSF")) /= 0)  then
     ABI_ERROR("Cannot produce file for Fermi surface in BXSF format. Check log file for info.")
   end if

 !case ("ebands_nesting")
   !call get_path_ebands_cryst(path, ebands, cryst, comm)
   !if (ebands_write_nesting(ebands, cryst, filepath, prtnest, tsmear, fermie_nest, qpath_vertices, errmsg) /= 0) then
   !  ABI_ERROR("Cannot produce file for nesting factor. Check log file for info.")
   !end if

 case ("skw_kpath")
   ! Get energies on the IBZ from path file
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ! Generate k-path
   ABI_CHECK(get_arg("ndivsm", ndivsm, msg, default=20) == 0, msg)
   !ABI_COMMENT("Using hard-coded k-path because nkpath not present in input file.")
   !nbounds = 5
   ABI_MALLOC(bounds, (3, 5))
   bounds = reshape([zero, zero, zero, half, zero, zero, zero, half, zero, zero, zero, zero, zero, zero, half], [3,5])
   !ABI_CHECK(get_args_from_string("bounds", bounds, msg, default="[0,0,0,1/2,0,0,0,1/2]") == 0, msg)
   kpath = kpath_new(bounds, cryst%gprimd, ndivsm)
   ABI_FREE(bounds)

   call kpath%print([std_out], header="Interpolating energies on k-path")

   ! Interpolate band energies with star-functions
   call parse_skw_params(skw_params)
   !ABI_CHECK(get_arg_list("einterp", ivec9, lenr, msg, exclude="ngkpt", want_len=9) == 0, msg)
   !if (nint(dtset%einterp(1)) == 1) skw_params = dtset%einterp
   ebands_kpath = ebands%interp_kpath(cryst, kpath, skw_params, [1, ebands%mband], comm)

   !call wrtout(std_out, sjoin(" Writing interpolated bands to:",  path)
   ABI_CHECK(get_arg("prtebands", prtebands, msg, default=2) == 0, msg)
   call ebands_kpath%write(prtebands, path)

   !ebands_kmesh = ebands%interp_kmesh(cryst, skw_params, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
   !     band_block, comm, out_prefix)
   !call ebands_kmesh%free()

 case ("skw_compare")
   ! Get energies on the IBZ from path
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ! Get ab-initio energies for the second file (assume k-path!)
   call get_path_ebands_cryst(other_path, other_ebands, other_cryst, comm, argpos=3)

   ! Interpolate band energies on the path with star-functions
   kpath = kpath_new(other_ebands%kptns, other_cryst%gprimd, -1)

   call parse_skw_params(skw_params)
   ebands_kpath = ebands%interp_kpath(cryst, kpath, skw_params, [1, ebands%mband], comm)

   ! Compare gaps
   ABI_CHECK(get_arg("is-metal", is_metal, msg, default=.False.) == 0, msg)
   if (.not. is_metal) then
     write(std_out, "(2a)")" Will try to compare gaps. Use --is-metal option to skip this check.",ch10
     call other_ebands%print_gaps([std_out], header="Ab-initio gaps")
     call ebands_kpath%print_gaps([std_out], header="SKW interpolated gaps")
   end if

   !ABI_CHECK(get_arg("prtebands", prtebands, msg, default=2) == 0, msg)
   !call ebands_write(ebands_kpath, prtebands, path)

   ! Write EBANDS file
   NCF_CHECK(other_ebands%ncwrite_path(cryst, "abinitio_EBANDS.nc"))
   NCF_CHECK(ebands_kpath%ncwrite_path(other_cryst, "skw_EBANDS.nc"))

   call wrtout(std_out, &
     ch10//" Use `abicomp.py ebands abinitio_EBANDS.nc skw_EBANDS.nc -p combiplot` to compare the bands with AbiPy.", &
     newlines=2)

 case ("ebands_mu_T")
   ! Get energies on the IBZ from filepath
   call get_path_ebands_cryst(path, ebands, cryst, comm)

   ABI_CHECK(get_arg("extrael", extrael, msg, default=zero) == 0, msg)
   if (extrael == zero) then
     ABI_CHECK(get_arg("doping", doping, msg, default=zero) == 0, msg)
     ! Units of eph_doping is e_charge / cm^3
     extrael = - doping * cryst%ucvol * (Bohr_meter * 100) ** 3
     write(std_out, *)"Adding doping", doping
   end if

   !ABI_CHECK(get_arg("occopt", occopt, msg, default=3) == 0, msg)
   !ABI_CHECK(get_arg("tsmear", tsmear, msg, default=tol2) == 0, msg)
   ABI_CHECK(get_arg("spinmagntarget", spinmagntarget, msg, default=-99.99_dp) == 0, msg)

   ebands%nelect = ebands%nelect + extrael
   gaps = ebands%get_gaps(ierr)

   ABI_CHECK(get_arg_list("tmesh", tmesh, lenr, msg, default_list=[5._dp, 59._dp, 6._dp] ) == 0, msg)
   ntemp = nint(tmesh(3))
   ABI_CHECK_IGEQ(ntemp, 1, "ntemp <= 0")
   ABI_MALLOC(kTmesh, (ntemp))
   kTmesh = arth(tmesh(1), tmesh(2), ntemp) * kb_HaK

   ABI_MALLOC(mu_e, (ntemp))
   call ebands%get_muT_with_fd(ntemp, kTmesh, spinmagntarget, prtvol, mu_e, comm)
   !mu_e = 6.715 * eV_Ha

   call gaps%print([std_out], header="KS gaps", kTmesh=kTmesh, mu_e=mu_e)
   !stop

   ABI_MALLOC(n_ehst, (2, ebands%nsppol, ntemp))
   call ebands%get_carriers(ntemp, kTmesh, mu_e, n_ehst)

   !write(msg, "(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
   do spin=1,ebands%nsppol
     do itemp=1,ntemp
       write(std_out, "(a, 2f16.2, 2e16.2)")&
         " T (K), mu_e (eV), nh, ne", kTmesh(itemp) / kb_HaK, mu_e(itemp) * Ha_eV, &
        n_ehst(2,spin,itemp) / cryst%ucvol / Bohr_cm**3, &
        n_ehst(1,spin,itemp) / cryst%ucvol / Bohr_cm**3
     end do
   end do

   ABI_CHECK(get_arg("intmeth", intmeth, msg, default=2) == 0, msg)
   ABI_CHECK(get_arg("step", step, msg, default=0.02 * eV_Ha) == 0, msg)
   ABI_CHECK(get_arg("broad", broad, msg, default=0.06 * eV_Ha) == 0, msg)

   edos = ebands%get_edos(cryst, intmeth, step, broad, comm)
   call edos%print(std_out, header="Electron DOS")
   call edos%get_carriers(ntemp, kTmesh, mu_e, n_ehst)

   !write(msg, "(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
   do spin=1,ebands%nsppol
     do itemp=1,ntemp
       write(std_out, "(a, 2f16.2, 2e16.2)")&
        " T (K), mu_e (eV), nh, ne", kTmesh(itemp) / kb_HaK, mu_e(itemp) * Ha_eV, &
        n_ehst(2, itemp, spin) / cryst%ucvol / Bohr_cm**3, &
        n_ehst(1, itemp, spin) / cryst%ucvol / Bohr_cm**3
     end do
   end do

   ABI_FREE(kTmesh)
   ABI_FREE(mu_e)
   ABI_FREE(n_ehst)

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
   !  ABI_ERROR(sjoin("While reading:", path, ch10, msg))
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

 case ("find_ewin")
   ltetra = 2
   ABI_CHECK(get_arg("ltetra", ltetra, msg, default=2) == 0, msg)
   call get_path_ebands_cryst(path, ebands, cryst, comm)
   ! Use Q-mesh == K-mesh for double delta.
   call find_ewin(ebands%nkpt, ebands%kptns, cryst, ebands, ltetra, fs_ewin, comm)

 ! ===========
 ! Unit tests
 ! ===========

 case ("tetra_unit_tests")
   ABI_CHECK(get_arg("ptgroup", ptgroup, msg, default="m-3m") == 0, msg)
   ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default_list=[20, 20, 20]) == 0, msg)
   ABI_CHECK(get_arg("use-symmetries", use_symmetries, msg, default=1) == 0, msg)
   call tetra_unittests(ptgroup, ngqpt, use_symmetries, prtvol, comm)

 case ("tetra_zinv_convergence")
   ABI_CHECK(get_arg("ptgroup", ptgroup, msg, default="m-3m") == 0, msg)
   !ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default_list=[20, 20, 20]) == 0, msg)
   ABI_CHECK(get_arg("use-symmetries", use_symmetries, msg, default=1) == 0, msg)
   call tetra_zinv_convergence(ptgroup, use_symmetries, comm)

 case ("kptrank_unit_tests")
   ABI_CHECK(get_arg("ptgroup", ptgroup, msg, default="m-3m") == 0, msg)
   ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default_list=[20, 20, 20]) == 0, msg)
   ABI_CHECK(get_arg("use-symmetries", use_symmetries, msg, default=1) == 0, msg)
   call kptrank_unittests(ptgroup, ngqpt, use_symmetries, comm)

 case default
   call abitk_show_help()
   ABI_ERROR(sjoin("Invalid command:", command))
 end select

 ! Deallocate memory to make memcheck happy.
 call hdr%free(); call cryst%free(); call other_cryst%free(); call kpath%free()
 call ebands%free(); call ebands_kpath%free(); call other_ebands%free()
 call edos%free(); call jdos%free(); call gaps%free()

 ABI_SFREE(kibz)
 ABI_SFREE(wtk)
 ABI_SFREE(kbz)
 ABI_SFREE(new_shiftk)
 ABI_SFREE(shiftk)

 call abinit_doctor("__abitk")

 100 call xmpi_end()

contains
!!***

!!****f* abitk/abitk_show_help
!! NAME
!! abitk_show_help
!!
!! FUNCTION
!!  Show command line help
!!
!! SOURCE

subroutine abitk_show_help()

  write(std_out,"(a)")" --version              Show version number and exit."
  write(std_out,"(a)")" -h, --help             Show this help and exit."
  write(std_out,"(a)")" -v                     Increase verbosity level"

  write(std_out,"(2a)")ch10,"=== HEADER ==="
  write(std_out,"(a)")"hdr_print FILE [--prtvol 0]  Print ABINIT header."

  write(std_out,"(2a)")ch10,"=== KPOINTS ==="
  write(std_out,"(a)")"ibz FILE --ngkpt 2 2 2 or --kptrlatt [--kptopt 1] [--shiftk 0.5 0.5 0.5] [--chksymbreak 1]"

  write(std_out,"(2a)")ch10,"=== CRYSTAL ==="
  write(std_out,"(a)")"crystal_print FILE                   Print info on crystalline structure."
  write(std_out,"(a)")"from_poscar POSCAR_FILE              Read POSCAR file, print abinit variables."

  write(std_out,"(2a)")ch10,"=== ELECTRONS ==="
  write(std_out,"(a)")"ebands_print FILE                    Print info on electron band structure."
  write(std_out,"(a)")"ebands_xmgrace FILE                  Produce XMGRACE file with electron bands."
  write(std_out,"(a)")"ebands_gnuplot FILE                  Produce GNUPLOT file with electron bands."
  write(std_out,"(a)")"ebands_edos FILE --intmeth, --step, --broad  Compute electron DOS."
  write(std_out,"(a)")"ebands_bxsf FILE                     Produce BXSF file for Xcrysden."
  !write(std_out,"(a)")"ebands_mu_t FILE --occopt --tsmear --extrael  Change number of electron, compute new Fermi level."
  write(std_out,"(a)")"ebands_gaps FILE                     Print info on gaps"
  !write(std_out,"(a)")"ebands_jdos FILE --intmeth, --step, --broad  Compute electron DOS."
  write(std_out,"(a)")"skw_kpath FILE                       Interpolate band structure along a (hardcoded) k-path."
  write(std_out,"(a)")"skw_compare IBZ_WFK KPATH_WFK        Use e_nk from IBZ_WFK to interpolate on the k-path in KPATH_WFK."

  write(std_out,"(2a)")ch10,"=== DEVELOPERS ==="
  write(std_out,"(a)")"tetra_unit_tests                      Run unit tests for tetrahedron routines."
  write(std_out,"(a)")"kptrank_unit_tests                    Run unit tests for kptrank routines."

end subroutine abitk_show_help
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
!*******************************************************

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
!*******************************************************

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
!!  Build a crystal object from a filepath.
!!
!! INPUTS
!!
!! OUTPUT
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
!*******************************************************

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
!! SOURCE

subroutine parse_skw_params(params)

!Arguments ----------------------------
 real(dp),intent(out) :: params(4)

!Local variables-----------------------
 integer :: lpratio
 real(dp) :: rcut, rsigma
!*******************************************************

 ABI_CHECK(get_arg("lpratio", lpratio, msg, default=5) == 0, msg)
 ABI_CHECK(get_arg("rcut", rcut, msg, default=zero) == 0, msg)
 ABI_CHECK(get_arg("rsigma", rsigma, msg, default=zero) == 0, msg)
 params = [one, dble(lpratio), rcut, rsigma]

end subroutine parse_skw_params
!!***

integer function vid(vname)
  character(len=*),intent(in) :: vname
  vid = nctk_idname(ncid, vname)
end function vid

end program abitk
!!***
