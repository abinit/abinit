!!****p*ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2026 ABINIT group (XG,DCA,JCC,CL,XW,GA,MR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program anaddb

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_xmpi
 use m_xomp
 use m_abicore
 use m_errors
 use m_argparse
 use m_nctk
 use netcdf

 use m_build_info,     only : abinit_version
 use m_io_tools,       only : open_file, flush_unit
 use m_fstrings,       only : int2char4, itoa, sjoin, strcat, inupper
 use m_specialmsg,     only : specialmsg_getcount, herald
 use m_time,           only : asctime, timein, timab, cwtime, cwtime_report
 use m_dtfil,          only : isfile
 use m_crystal,        only : crystal_t
 use m_ddb,            only : ddb_type, asrq0_t, ddb_lw_copy
 use m_ddb_hdr,        only : ddb_hdr_type
 use m_ifc,            only : ifc_type
 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_anaddb_driver, only : anaddb_driver_type
 use m_ddb_interpolate, only : ddb_interpolate
 use m_elphon,         only : elphon
 use m_thmeig,         only : thmeig
 use m_phonons,        only : mkphbs
 use m_gruneisen,      only : gruns_anaddb

 implicit none

!Local variables-------------------------------
 integer, parameter:: master = 0
 integer:: comm, ii, ierr
 integer:: nproc, my_rank, ana_ncid
 logical:: iam_master
 integer:: units(2)
 real(dp):: tcpu, tcpui, twall, twalli !,cpu, wall, gflops
 real(dp):: tsec(2)
 character(len=10):: procstr
 character(len=24):: codename, start_datetime
 character(len=fnlen):: worker_logfile
 character(len=500):: msg
 type(args_t):: args
 type(anaddb_dataset_type):: dtset
 type(anaddb_driver_type):: driver
 type(crystal_t):: crystal
 type(ifc_type):: Ifc
 type(ddb_type):: ddb
 type(ddb_type):: ddb_lw
 type(ddb_hdr_type):: ddb_hdr
 type(asrq0_t):: asrq0

! ========================================================================== !

! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm = xmpi_world)

! These units are defined in defs_basis
 units = [std_out, ab_out]

! Initialize MPI
 call xmpi_init()

! MPI variables
 comm = xmpi_world; nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

! Initialize memory profiling if activated at configure time.
! if a full report is desired, set the argument of abimem_init to "2" instead of "0" via the command line.
! note that the file can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level, limit_mb = args%abimem_limit_mb)
#endif

! Initialisation of the timing
 call timein(tcpui, twalli)

 if (iam_master) then
   codename='ANADDB'//repeat(' ',18)
   call herald(codename, abinit_version, std_out)
 end if

 start_datetime = asctime()

! Zero out all accumulators of time and init timers
 call timab(1, 0, tsec)

! Initialise the code: write heading, and read names of files.
 if (iam_master) then
   call dtset%init(args%input_path)
 end if

! Broadcast file names
 call dtset%bcast_files(comm)

! make log file for non-master procs
 if (.not. iam_master) then
   call int2char4(my_rank, procstr)
   ABI_CHECK((procstr(1:1)/='#'), 'Bug: string length too short!')
   worker_logfile = trim(dtset%filename_output) // "_LOG_P" // trim(procstr)
   if (open_file(worker_logfile, msg, unit = std_out, form="formatted", action="write") /= 0) then
     ABI_ERROR(msg)
   end if
 end if

! ========================================================================== !
! Read input variables
 call dtset%read_input(comm)

 if (args%dry_run /= 0) then
   call wrtout(std_out, "Dry run mode. Exiting after have read the input")
   call dtset%free()
   goto 100
 end if

! ========================================================================== !
! Open output file
 if (iam_master) then
   call isfile(dtset%filename_output, 'new')
   if (open_file(dtset%filename_output, msg, unit=ab_out, form='formatted', status='new') /= 0) then
     ABI_ERROR(msg)
   end if
   rewind (unit = ab_out)
   call herald(codename, abinit_version, ab_out)

   ! Echo the inputs to console and main output file
   call dtset%outvars(std_out)
   call dtset%outvars(ab_out)
 else
   ab_out = dev_null
 end if

! =========================================================================== !

! Initialize driver
 call driver%init(dtset)

! Read the DDB information and symmetrize partially the DDB
 write(msg, '(a, a)' )' read the DDB information and perform some checks',ch10
 call wrtout(units, msg)

 call ddb%from_file(dtset%filename_ddb, ddb_hdr, crystal, comm, prtvol=dtset%prtvol)

! Change the bravais lattice if needed
 call ddb%set_brav(dtset%brav)

! Copy the long-wave ddb
 if (ddb_hdr%has_d3E_lw) then
   call ddb_lw_copy(ddb, ddb_lw, ddb_hdr)
 end if

! Acoustic Sum Rule
 call asrq0%init(ddb, dtset%asr, dtset%rfmeth, crystal%xcart)

! Open netcdf output and write basic quantities
 call driver%open_write_nc(ana_ncid, dtset, crystal, comm)

! =========================================================================== !

! Compute dielectric tensor, Born effective charges, and quadrupoles.
 if (driver%do_electric_tensors) then
   call driver%electric_tensors(dtset, crystal, ddb, ddb_lw, ddb_hdr, ana_ncid, comm)
 end if

! Structural response at fixed polarization
 if (dtset%polflag == 1) then
   call driver%structural_response(dtset, crystal, ddb)
 end if

! Compute non-linear optical susceptibilities
! and first-order change in the linear dielectric susceptibility
 if (dtset%nlflag > 0) then
   call driver%susceptibilities(dtset, ddb, ana_ncid, comm)
 end if

! Interatomic force constants
 if (driver%do_ifc) then
   call driver%interatomic_force_constants(Ifc, dtset, crystal, ddb, ana_ncid, comm)
 end if

! Phonon density of states
 if (driver%do_phonon_dos) then
   call driver%phdos(dtset, crystal, Ifc, comm)
 end if

! Phonon density of states and thermodynamical properties calculation
 if (dtset%ifcflag == 1 .and. any(dtset%thmflag==[1, 2])) then
   call driver%harmonic_thermo(dtset, crystal, Ifc, ddb, comm)
 end if

! Phonon band structure
 if (driver%do_phonon_bs) then
   call mkphbs(Ifc, crystal, dtset, ddb, asrq0, dtset%prefix_outdata, comm)
 end if

! DDB interpolation
 if (dtset%prtddb == 1 .and. dtset%ifcflag == 1) then
   call ddb_interpolate(Ifc, crystal, dtset, ddb, ddb_hdr, asrq0, comm)
 end if

! =========================================================================== !
! Electron-phonon section

 if (dtset%elphflag == 1) then
   call elphon(dtset, crystal, Ifc, comm)
 end if

 ! Thermal supercell calculation
 if (sum(abs(dtset%thermal_supercell))>0 .and. dtset%ifcflag == 1) then
   call driver%thermal_supercell(dtset, crystal, ifc)
 end if

! Thermal corrections to eigenvalues (old)
 if (dtset%thmflag >= 3 .and. dtset%thmflag <= 8) then
   call thmeig(dtset, ddb, crystal, ab_out, crystal%natom, dtset%mpert, dtset%msize, asrq0%d2asr, comm)
 end if

! =========================================================================== !

! Compute the dielectric function and oscillator strength.
 if (driver%do_dielectric_q0) then
   call driver%dielectric_q0(dtset, crystal, ifc, ddb, asrq0, ana_ncid, comm)
 end if

! Non-linear response: electrooptic and Raman (q = Gamma, TO modes only)
 if (dtset%nlflag == 1) then
   call driver%nonlinear_response(dtset, crystal, ana_ncid, comm)
 end if

! Non-analyticity in the dynamical matrix
 if (driver%do_dielectric_nonana) then
   call driver%dielectric_nonana(dtset, crystal, ddb, ana_ncid, comm)
 end if

! =========================================================================== !
! Linear response with strain

! Internal strain (needed for the other linear response functions)
 if (dtset%instrflag /= 0) then
   call driver%internal_strain(dtset, ddb, asrq0)
 end if

! Elastic tensor
 if (dtset%elaflag /= 0) then
   call driver%elastic_tensor(dtset, crystal, ddb, asrq0, ana_ncid)
 end if

! Piezoelectric tensor
 if (dtset%piezoflag /= 0 .or. dtset%dieflag == 4 .or. dtset%elaflag == 4) then
   call driver%piezoelectric_tensor(dtset, crystal, ddb, ana_ncid)
 end if

! Flexoelectric tensor
 if (dtset%flexoflag /= 0) then
   call driver%flexoelectric_tensor(dtset, crystal, ddb, ddb_lw, ddb_hdr, asrq0)
 end if

! =========================================================================== !

 ! Gruneisen parameters
 if (dtset%gruns_nddbs /= 0) then
   call gruns_anaddb(dtset, comm)
 end if

 ! Lattice Wannier functions
 if (dtset%ifcflag == 1 .and. dtset%lwfflag > 0 ) then
   call driver%lattice_wannier(dtset, crystal, Ifc, comm)
 endif

 ! Output phonon frequencies for BoltzTrap
 if (iam_master .and. dtset%ifcflag == 1 .and. dtset%outboltztrap == 1) then
   call ifc%outphbtrap(crystal, dtset%ng2qpt, 1, dtset%q2shft, dtset%prefix_outdata)
 end if

! =========================================================================== !
! Close netcdf file
 if (iam_master) then
   NCF_CHECK(nf90_close(ana_ncid))
 end if

! =========================================================================== !
! Free memory
 call asrq0%free()
 call ifc%free()
 call crystal%free()
 call ddb%free()
 call ddb_hdr%free()
 call ddb_lw%free()
 call driver%free()
 call dtset%free()

! =========================================================================== !
! Output timing and memory reports, then close output files

 call timein(tcpu, twall)
 tsec(1)=tcpu-tcpui; tsec(2)=twall-twalli
 write(msg, '(a, i4, a, f13.1, a, f13.1)' )' Proc.',my_rank, ' individual time (sec): cpu=',tsec(1), '  wall=',tsec(2)
 call wrtout(std_out, msg)

 if (iam_master) then
   write(ab_out, '(a, a, a, i4, a, f13.1, a, f13.1)' )'-',ch10, &
    '- Proc.',my_rank, ' individual time (sec): cpu=',tsec(1), '  wall=',tsec(2)
 end if

 call xmpi_sum(tsec, comm, ierr)

 write(msg, '(a, (80a), a, a, a, f11.3, a, f11.3, a, a, a, a)' ) ch10, &
  ('=',ii = 1, 80), ch10, ch10, &
   '+Total cpu time',tsec(1), '  and wall time',tsec(2), ' sec',ch10, ch10, &
   ' anaddb : the run completed successfully.'
 call wrtout(units, msg)

 if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
   write(std_out, "(a)")"--- !FinalSummary"
   write(std_out, "(a)")"program: anaddb"
   write(std_out, "(2a)")"version: ",trim(abinit_version)
   write(std_out, "(2a)")"start_datetime: ",start_datetime
   write(std_out, "(2a)")"end_datetime: ",asctime()
   write(std_out, "(a, f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out, "(a, f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out, "(a, i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out, "(a, i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out, "(a, i0)")"num_warnings: ",nwarning
   !write(std_out, "(a, i0)")"num_comments: ",ncomment
   write(std_out, "(a)")"..."
   call flush_unit(std_out)
 end if

 ! Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor(dtset%filename_output)

 call flush_unit(ab_out)
 call flush_unit(std_out)

 if (iam_master) close(ab_out)

 100 call xmpi_end()

 end program anaddb
!!***
