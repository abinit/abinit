!!****p*ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group (XG,DCA,JCC,CL,XW,GA,MR)
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
 use m_ifc
 use m_ddb
 use m_ddb_hdr
 use m_phonons
 use m_gruneisen
 use m_supercell
 use m_nctk
 use netcdf

 use m_build_info,     only : abinit_version
 use m_io_tools,       only : open_file, flush_unit
 use m_fstrings,       only : int2char4, itoa, sjoin, strcat, inupper
 use m_specialmsg,     only : specialmsg_getcount, herald
 use m_time,           only : asctime, timein, timab, cwtime, cwtime_report
 use m_parser,         only : instrng
 use m_dtfil,          only : isfile
 use m_anaddb_dataset, only : anaddb_init, anaddb_dataset_type, anaddb_dtset_free, outvars_anaddb, invars9
 use m_ddb_interpolate, only : ddb_interpolate
 use m_crystal,        only : crystal_t
 use m_dynmat,         only : gtdyn9, dfpt_phfrq, dfpt_prtph
 use m_elphon,         only : elphon
 use m_harmonic_thermo, only : harmonic_thermo
 use m_thmeig,         only : thmeig
 use m_symfind,        only : symanal
 use m_raman,          only : ramansus, electrooptic
 use m_ddb_diel,       only : ddb_diel
 use m_relaxpol,       only : relaxpol
 use m_ddb_elast,      only : ddb_elast
 use m_ddb_piezo,      only : ddb_piezo
 use m_ddb_internalstr, only : ddb_internalstr
 use m_ddb_flexo
 use m_lwf,            only : run_lattice_wannier

 implicit none

!Local variables-------------------------------
 integer:: msym !  msym = maximum number of symmetry elements in space group
!Define input and output unit numbers (some are defined in defs_basis-all should be there ...):
 integer, parameter:: master = 0
 integer:: comm, iblok, iblok_stress, iblok_epsinf, iblock_quadrupoles, ii
 integer:: ierr, iphl2, lenstr, lwsym, mtyp, mpert, msize, natom
 integer:: nsym, ntypat, usepaw, nproc, my_rank, ana_ncid, prt_internalstr
 integer:: phdos_ncid, ncerr
 logical:: iam_master
 integer:: rfelfd(4), rfphon(4), rfstrs(4), ngqpt_coarse(3), units(2)
 integer:: count_wminmax(2)
 integer, allocatable:: d2flg(:)
 real(dp):: etotal, tcpu, tcpui, twall, twalli !,cpu, wall, gflops
 real(dp):: epsinf(3, 3), dielt_rlx(3, 3)
 real(dp):: compl(6, 6), compl_clamped(6, 6), compl_stress(6, 6)
 real(dp):: elast(6, 6), elast_clamped(6, 6), elast_stress(6, 6)
 real(dp):: red_ptot(3), pel(3)
 real(dp):: piezo(6, 3), qphnrm(3), qphon(3, 3), strten(6), tsec(2)
 real(dp):: wminmax(2)
 real(dp), allocatable:: d2cart(:,:), dchide(:,:,:), lst(:)
 real(dp), allocatable:: dchidt(:,:,:,:), displ(:), eigval(:,:)
 real(dp), allocatable:: eigvec(:,:,:,:,:), fact_oscstr(:,:,:), instrain(:,:)
 real(dp), allocatable:: gred(:,:), phfrq(:)
 real(dp), allocatable:: rsus(:,:,:)
 real(dp), allocatable:: zeff(:,:,:)
 real(dp), allocatable:: qdrp_cart(:,:,:,:)
 character(len = 10):: procstr
 character(len = 24):: codename, start_datetime
 character(len = strlen):: string, raw_string
 character(len = fnlen):: filnam(8), elph_base_name, tmpfilename, phibz_prefix
 character(len = 500):: msg
 type(args_t):: args
 type(anaddb_dataset_type):: inp
 type(phdos_t):: Phdos
 type(ifc_type):: Ifc, Ifc_coarse
 type(ddb_type):: ddb
 type(ddb_type):: ddb_lw
 type(ddb_hdr_type):: ddb_hdr
 type(asrq0_t):: asrq0
 type(crystal_t):: crystal
 type(supercell_type), allocatable:: thm_scells(:)

!******************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm = xmpi_world)
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
 if (iam_master) call anaddb_init(args%input_path, filnam)
 call xmpi_bcast(filnam, master, comm, ierr)

 ! make log file for non-master procs
 if (.not. iam_master) then
   call int2char4(my_rank, procstr)
   ABI_CHECK((procstr(1:1)/='#'), 'Bug: string length too short!')
   tmpfilename = trim(filnam(2)) // "_LOG_P" // trim(procstr)
   if (open_file(tmpfilename, msg, unit = std_out, form="formatted", action="write") /= 0) then
     ABI_ERROR(msg)
   end if
 end if

!******************************************************************

 ! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 call ddb_hdr%open_read(filnam(3), comm = comm, dimonly = 1)

 usepaw = ddb_hdr%usepaw
 natom = ddb_hdr%natom
 ntypat = ddb_hdr%ntypat
 mtyp = ddb_hdr%mblktyp
 mpert = ddb_hdr%mpert
 msize = ddb_hdr%msize
 call ddb_hdr%free()

 ! Read the input file, and store the information in a long string of characters
 ! strlen from defs_basis module
 if (iam_master) then
   call instrng(filnam(1), lenstr, 1, strlen, string, raw_string)
   ! To make case-insensitive, map characters to upper case.
   call inupper(string(1:lenstr))
 end if

 call xmpi_bcast(string, master, comm, ierr)
 call xmpi_bcast(raw_string, master, comm, ierr)
 call xmpi_bcast(lenstr, master, comm, ierr)

 ! Save input string in global variable so that we can access it in ntck_open_create
 ABI_MALLOC_TYPE_SCALAR(character(len=len_trim(raw_string)), INPUT_STRING)
 INPUT_STRING = trim(raw_string)

 ! Read the inputs
 call invars9(inp, lenstr, natom, string)

 if (args%dry_run /= 0) then
   call wrtout(std_out, "Dry run mode. Exiting after have read the input")
   goto 100
 end if

 ! Open output files and ab_out (might change its name if needed)
 ! MJV 1/2010 : now output file is open, but filnam(2) continues unmodified

 ! so the other output files are overwritten instead of accumulating.
 if (iam_master) then
   tmpfilename = filnam(2)
   call isfile(tmpfilename, 'new')
   if (open_file(tmpfilename, msg, unit = ab_out, form='formatted',status='new') /= 0) then
     ABI_ERROR(msg)
   end if
   rewind (unit = ab_out)
   call herald(codename, abinit_version, ab_out)

   ! Echo the inputs to console and main output file
   call outvars_anaddb(inp, std_out)
   call outvars_anaddb(inp, ab_out)
 else
   ab_out = dev_null
 end if

 ! Check the value and transform the meaning of atifc (1 and 0 only)
 call chkin9(inp%atifc,inp%natifc,natom)

!******************************************************************
!******************************************************************
! Read the DDB information, also perform some checks, and symmetrize partially the DDB
 write(msg, '(a, a)' )' read the DDB information and perform some checks',ch10
 call wrtout(units, msg)

 call ddb%from_file(filnam(3), ddb_hdr, crystal, comm, prtvol = inp%prtvol)
 call ddb_hdr%free()

 ! Change the bravais lattice if needed
 call ddb%set_brav(inp%brav)

 nsym = crystal%nsym

 ! MR: a new ddb is necessary for the longwave quantities due to incompability of it with automatic reshapes
 ! that ddb%val and ddb%flg experience when passed as arguments of some routines
 ! GA: Should replace with ddb_hdr%with_d3E_lw
 if (mtyp == BLKTYP_d3E_lw) then
   call ddb_lw_copy(ddb, ddb_lw, mpert, natom, ntypat)
 end if

 ! Acoustic Sum Rule
 ! In case the interatomic forces are not calculated, the
 ! ASR-correction (asrq0%d2asr) has to be determined here from the Dynamical matrix at Gamma.
 asrq0 = ddb%get_asrq0(inp%asr, inp%rfmeth, crystal%xcart)

 ! TODO: This is to maintain the previous behaviour in which all the arrays were initialized to zero.
 ! In the new version asrq0%d2asr is always computed if the Gamma block is present
 ! and this causes changes in [v5][t28]
 if (.not. (inp%ifcflag == 0 .or. inp%instrflag /= 0 .or. inp%elaflag /= 0)) then
   asrq0%d2asr = zero
   if (asrq0%asr == 3 .or. asrq0%asr == 4) then
     asrq0%singular = zero; asrq0%uinvers = zero; asrq0%vtinvers = zero
   end if
 end if

 ! Open the netcdf file that will contain the anaddb results
 ana_ncid = nctk_noid
 if (iam_master) then
   NCF_CHECK_MSG(nctk_open_create(ana_ncid, trim(filnam(8))//"_anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
   ncerr = nctk_def_dims(ana_ncid, [ &
       nctkdim_t('number_of_atoms', natom), &
       nctkdim_t('natom3', 3*natom), &
       nctkdim_t('number_of_phonon_modes', 3*natom), &
       nctkdim_t('anaddb_input_len', lenstr) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_arrays(ana_ncid, [ &
     nctkarr_t("anaddb_input_string", "char", "anaddb_input_len") &
   ])
   NCF_CHECK(ncerr)
   !NCF_CHECK(nctk_defnwrite_ivars(ana_ncid, ["anaddb_version"], [1]))
   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "anaddb_input_string"), string(:lenstr)))
   NCF_CHECK(crystal%ncwrite(ana_ncid))
 end if

 ! Calculation of Gruneisen parameters.
 if (inp%gruns_nddbs /= 0) then
   call gruns_anaddb(inp, filnam(2), comm)
   goto 50
 end if

 ABI_MALLOC(instrain, (3*natom, 6))
 ABI_MALLOC(d2cart, (2, msize))
 ABI_MALLOC(displ, (2*3*natom*3*natom))
 ABI_MALLOC(eigval, (3, natom))
 ABI_MALLOC(eigvec, (2, 3, natom, 3, natom))
 ABI_MALLOC(phfrq, (3*natom))
 ABI_MALLOC(zeff, (3, 3, natom))
 ABI_MALLOC(qdrp_cart, (3, 3, 3, natom))

!**********************************************************************
!**********************************************************************

 ! Get Quadrupole tensor
 iblock_quadrupoles = 0
 qdrp_cart = zero
 if (mtyp == BLKTYP_d3E_lw) then
   write(msg, '(2a, (80a), 2a)') ch10, ('=',ii = 1, 80)
   call wrtout(units, msg)
   lwsym = 1
   iblock_quadrupoles = ddb_lw%get_quadrupoles(ddb_hdr%ddb_version, lwsym, BLKTYP_d3E_lw, qdrp_cart)
 end if

 ! The default value is 1. Here we set the flags to zero if Q*is not available.
 if (iblock_quadrupoles == 0) then
   inp%dipquad = 0
   inp%quadquad = 0
 end if

 ! Get the electronic dielectric tensor (epsinf) and Born effective charges (zeff)
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblok = ddb%get_dielt_zeff(crystal, inp%rfmeth, inp%chneut, inp%selectz, epsinf, zeff)

 ! Try to get epsinf, in case just the DDE are present
 if (iblok == 0) then
   iblok_epsinf = ddb%get_dielt(inp%rfmeth, epsinf)
 end if

 !if (iblok_epsinf == 0) then
 if (epsinf(1, 1)==one .and. epsinf(2, 2)==one .and. epsinf(3, 3)==one) then
   write(msg, '(5a)') ch10, &
     ' The DDB file does not contain the derivatives w.r.t. electric field perturbation. ',ch10, &
     ' The program will continue by setting the electronic dielectric tensor to 1. ',ch10
  ! call wrtout([ab_out], msg)
   ABI_WARNING(msg)
 end if

 if (my_rank == master) then
   ncerr = nctk_def_arrays(ana_ncid, [&
   nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'), &
   nctkarr_t('quadrupoles_cart', "dp", 'three, three, three, number_of_atoms'), &
   nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")],&
   defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ana_ncid, [character(len = nctk_slen) :: &
       "asr", "chneut", "dipdip", "symdynmat", "dipquad", "quadquad"])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'emacro_cart'), epsinf))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'quadrupoles_cart'), qdrp_cart))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'becs_cart'), zeff))
   ncerr = nctk_write_iscalars(ana_ncid, [character(len = nctk_slen) :: &
     "asr", "chneut", "dipdip", "symdynmat", "dipquad", "quadquad"], &
     [inp%asr, inp%chneut, inp%dipdip, inp%symdynmat, inp%dipquad, inp%quadquad])
   NCF_CHECK(ncerr)
 end if

!***************************************************************************
! Structural response at fixed polarization
!***************************************************************************
 if (inp%polflag == 1) then
   ABI_MALLOC(d2flg, (msize))

   if(iblok /= 0)then
     ! Save the second-order derivatives
     d2cart(1:2, 1:msize) = ddb%val(1:2, 1:msize, iblok)
     d2flg(1:msize) = ddb%flg(1:msize, iblok)

   else
     ! the gamma blok has not been found
     if (inp%relaxat == 0 .and. inp%relaxstr == 0) then
       ! The gamma blok is not needed
       d2cart(1:2, 1:msize)=zero
       d2flg(1:msize)=1
     else
       ! There is a problem !
       write(msg, '(7a)' )&
        'The dynamical matrix at Gamma is needed, in order to perform ',ch10, &
        "relaxation at constant polarisation (Na Sai's method)",ch10, &
        'However, this was not found in the DDB.',ch10, &
        'Action: complete your DDB with the dynamical matrix at Gamma.'
       ABI_ERROR(msg)
     end if
   end if  ! iblok not found

   ! Extract the block with the total energy
   if (ddb%get_etotal(etotal) == 0) then
     ABI_ERROR("DDB file does not contain GS etotal")
   end if

   ! Extract the polarizability
   iblok = ddb%get_pel(pel, inp%relaxat, inp%relaxstr)

   ! Extract the forces
   iblok = ddb%get_gred(gred, inp%relaxat, inp%relaxstr)

   ! Extract the stress tensor
   iblok = ddb%get_strten(strten, inp%relaxat, inp%relaxstr)

   ! when called from here red_ptot is not set  ! So set it to zero
   red_ptot(:)=zero

   call relaxpol(crystal, d2flg, d2cart, etotal, gred, inp%iatfix, &
&   ab_out, inp%istrfix, mpert, msize, inp%natfix, natom, &
&   inp%nstrfix, pel, red_ptot, inp%relaxat, inp%relaxstr, &
&   strten, inp%targetpol, usepaw)

   ABI_SFREE(gred)
   ABI_FREE(d2flg)
 end if

!***************************************************************************
! Compute non-linear optical susceptibilities and, if inp%nlflag < 3,
! First-order change in the linear dielectric susceptibility induced by an atomic displacement
!***************************************************************************
 if (inp%nlflag > 0) then
   ABI_MALLOC(dchide, (3, 3, 3))
   ABI_MALLOC(dchidt, (natom, 3, 3, 3))
   ABI_MALLOC(rsus, (3*natom, 3, 3))

   if (ddb%get_dchidet(inp%ramansr, inp%nlflag, dchide, dchidt) == 0) then
     ABI_ERROR("Cannot find block corresponding to non-linear optical susceptibilities in DDB file")
   end if

   ! Save to the netcdf
   if (my_rank == master) then
     ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchide", "dp", "three, three, three")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ana_ncid))
     NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchide"), dchide))

     ! dchidt only present if nlflag == 1 or 2
     if (inp%nlflag < 3) then
       ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchidt", "dp", &
       "number_of_atoms, three, three, three")], defmode=.True.)
       NCF_CHECK(ncerr)
       NCF_CHECK(nctk_set_datamode(ana_ncid))
       NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchidt"), dchidt))
     end if
   end if

 end if  ! nlflag

!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
 if (inp%ifcflag == 1) then
   ! ifc to be calculated for interpolation
   write(msg, '(a, a, (80a), a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, &
    ' Calculation of the interatomic forces ',ch10
   call wrtout(units, msg)

! TODO : check if this wrtout should be removed in latest merge 17 feb 2017
   call timein(tcpu, twall)
   write(msg, '(a, f11.3, a, f11.3, a)' )'-begin at tcpu',tcpu-tcpui, '  and twall',twall-twalli, ' sec'
   call wrtout(units, msg)

   if (any(inp%qrefine(:) > 1)) then
     ! Gaal-Nagy's algorithm in PRB 73 014117 [[cite:GaalNagy2006]]
     ! Build the IFCs using the coarse q-mesh.
     do ii = 1, 3
       ngqpt_coarse(ii) = inp%ngqpt(ii) / inp%qrefine(ii)
     end do
     call Ifc_coarse%init(crystal, ddb, &
       inp%brav, inp%asr, inp%symdynmat, inp%dipdip, inp%rfmeth, ngqpt_coarse, inp%nqshft, inp%q1shft, epsinf, zeff, qdrp_cart, &
       inp%nsphere, inp%rifcsph, inp%prtsrlr, inp%enunit, comm, dipquad = inp%dipquad, quadquad = inp%quadquad)

     ! Now use the coarse q-mesh to fill the entries in dynmat(q)
     ! on the dense q-mesh that cannot be obtained from the DDB file.
     call ifc%init(crystal, ddb, &
      inp%brav, inp%asr, inp%symdynmat, inp%dipdip, inp%rfmeth, &
      inp%ngqpt(1:3), inp%nqshft, inp%q1shft, epsinf, zeff, qdrp_cart, &
      inp%nsphere, inp%rifcsph, inp%prtsrlr, inp%enunit, comm, &
      Ifc_coarse = Ifc_coarse, dipquad = inp%dipquad, quadquad = inp%quadquad)
     call Ifc_coarse%free()

   else
     call ifc%init(crystal, ddb, &
       inp%brav, inp%asr, inp%symdynmat, inp%dipdip, inp%rfmeth, &
       inp%ngqpt(1:3), inp%nqshft, inp%q1shft, epsinf, zeff, qdrp_cart, &
       inp%nsphere, inp%rifcsph, inp%prtsrlr, inp%enunit, comm, dipquad = inp%dipquad, quadquad = inp%quadquad)
   end if

   call ifc%print(unit = std_out)

   ! Compute speed of sound.
   if (inp%vs_qrad_tolkms(1) > zero) call ifc%speedofsound(crystal, inp%vs_qrad_tolkms, ana_ncid, comm)

   ! Print analysis of the real-space interatomic force constants
   ! TODO: ifc_out should not have side effects
   if (my_rank == master .and. inp%ifcout /= 0) then
     call ifc%write(inp%ifcana, inp%atifc, inp%ifcout, inp%prt_ifc, ana_ncid, filnam(8))
   end if
 end if

!***************************************************************************
! Electron-phonon section
!***************************************************************************
 if (inp%elphflag == 1) then
   call elphon(inp, crystal, Ifc, filnam, comm)
 end if

!***************************************************************************

 if (sum(abs(inp%thermal_supercell))>0 .and. inp%ifcflag == 1) then
   ABI_MALLOC(thm_scells, (inp%ntemper))
   call zacharias_supercell_make(crystal, Ifc, inp%ntemper, inp%thermal_supercell, inp%tempermin, inp%temperinc, thm_scells)
   call zacharias_supercell_print(filnam(8), inp%ntemper, inp%tempermin, inp%temperinc, thm_scells)
 end if

!***************************************************************************
! Phonon density of states calculation, Start if interatomic forces have been calculated
!***************************************************************************
 if (inp%ifcflag == 1 .and. any(inp%prtdos==[1, 2])) then
   write(msg, '(a, (80a), 4a)')ch10, ('=',ii = 1, 80), ch10, ch10, ' Calculation of phonon density of states ',ch10
   call wrtout(units, msg)

   ! Only 1 shift in q-mesh
   wminmax = zero
   phibz_prefix = ""
   !phibz_prefix = "freq_displ" ! Uncomment this line to activate output of PHIBZ
   do
     call Phdos%init(crystal, Ifc, inp%prtdos, inp%dosdeltae, inp%dossmear, inp%ng2qpt, 1, inp%q2shft, &
                     phibz_prefix, wminmax, count_wminmax, comm, dos_maxmode = inp%dos_maxmode)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05; wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(msg, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
   end do

   if (iam_master) then
     call phdos%print_msqd(filnam(8), inp%ntemper, inp%tempermin, inp%temperinc)
     call phdos%print(strcat(filnam(8), "_PHDOS"))
     call phdos%print_debye(crystal%ucvol)
     call phdos%print_thermo(strcat(filnam(8), "_THERMO"), inp%ntemper, inp%tempermin, inp%temperinc)

     ncerr = nctk_open_create(phdos_ncid, strcat(filnam(8), "_PHDOS.nc"), xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, "Creating PHDOS.nc file")
     NCF_CHECK(crystal%ncwrite(phdos_ncid))
     call phdos%ncwrite(phdos_ncid)
     NCF_CHECK(nf90_close(phdos_ncid))
   end if

   call phdos%free()
 end if

 if (iam_master .and. inp%ifcflag == 1 .and. inp%outboltztrap == 1) then
   call ifc%outphbtrap(crystal, inp%ng2qpt, 1, inp%q2shft, filnam(8))
 end if

 ! Phonon density of states and thermodynamical properties calculation
 ! Start if interatomic forces and thermal flags are on
 if (inp%ifcflag == 1 .and. any(inp%thmflag==[1, 2])) then

   write(msg, '(a, (80a), a, a, a, a, a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, &
    ' Calculation of phonon density of states, ',ch10, &
    '    thermodynamical properties, ',ch10, &
    '    and Debye-Waller factors.',ch10
   call wrtout(units, msg)

   if (inp%thmflag == 1) then
     call harmonic_thermo(Ifc, crystal, ddb%amu, inp, ab_out, filnam(8), comm)

   else if (inp%thmflag == 2) then
     write(msg, '(a, (80a), a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, ' Entering thm9 routine with thmflag = 2 ',ch10
     call wrtout(units, msg)
     call harmonic_thermo(Ifc, crystal, ddb%amu, inp, ab_out, filnam(8), comm, thmflag = inp%thmflag)
   end if
 end if

!***************************************************************************
! Lattice Wannier function section.
! Compute the Dynamical matrix for a dense Q-mesh
! Input the eigenvectors and eigenvalues to the Lattcie Wannier module
! Construct the Lattice Wannier functions
!***************************************************************************
 if (inp%ifcflag == 1 .and. inp%lwfflag > 0 ) then
   write(msg, '(a, (80a), 4a)')ch10, ('=',ii = 1, 80), ch10, ch10, ' Calculation of lattice Wannier functions ',ch10
   call wrtout(units, msg)
   call run_lattice_wannier(ifc = Ifc, crystal = crystal, dtset = inp, prefix = filnam(8), comm = comm)
   write(msg, '(a, (80a))')ch10, ('=',ii = 1, 80)
   call wrtout(units, msg)
 endif



!**********************************************************************
! Treat the first list of vectors (without non-analyticities)
! (print the phonon freq. at each qpt (and eigenvectors if asked by the user)
! and the mode characters (Gamma only as of 12.04.2020)
!**********************************************************************
 call mkphbs(Ifc, crystal, inp, ddb, asrq0, filnam(8), comm)


 ! Interpolate the DDB onto the first list of vectors and write the file.
 if (inp%prtddb == 1 .and. inp%ifcflag == 1) then
   call ddb_hdr%open_read(filnam(3), comm)
   call ddb_hdr%close()
   ddb_hdr%crystal%space_group = crystal%space_group  ! GA: the space group is not written in the DDB text file.

   call ddb_interpolate(Ifc, crystal, inp, ddb, ddb_hdr, asrq0, filnam(8), comm)
   call ddb_hdr%free()
 end if


!***********************************************************************
! Thermal flags
!***********************************************************************

 if (inp%thmflag >= 3 .and. inp%thmflag <= 8) then
   elph_base_name = trim(filnam(8))//"_ep"
   call thmeig(inp, ddb, crystal, elph_base_name, filnam(5), ab_out, natom, mpert, msize, asrq0%d2asr, comm)
 end if


!**********************************************************************
! q = Gamma quantities (without non-analycities):
! - Dielectric constant calculations (diefalg options)
! and related properties: mode effective charges, oscillator strength
! - Raman tensor (at q = 0 with only TO modes) and EO coef. (nlflag)
!**********************************************************************

 ABI_MALLOC(fact_oscstr, (2, 3, 3*natom))
 ABI_MALLOC(lst, (inp%nph2l+1))
 lst(:)=zero

 ! Print the electronic contribution to the dielectric tensor
 ! It can be extracted directly from the DDB if perturbation with E-field is present
 if ((inp%dieflag /= 0 .and. inp%dieflag /= 2) .or. inp%nph2l /= 0 .or. inp%nlflag == 1) then

  !***************************************************************
  ! Generates the dynamical matrix at Gamma
  ! TODO: Check if we can avoid recomputing the phonon freq and eigendispla at Gamma becasue
  ! it is already done before in this routine. (EB)
  ! The problem is that it is done through mkphbs, which has only printing and does not return anything as out... (EB)

  qphon(:,1)=zero; qphnrm(1)=zero
  ! Generation of the dynamical matrix in cartesian coordinates
  if (inp%ifcflag == 1) then
    ! Get d2cart using the interatomic forces and the
    ! long-range coulomb interaction through Ewald summation
    call gtdyn9(ddb%acell, Ifc%atmfrc, epsinf, Ifc%dipdip, &
      Ifc%dyewq0, d2cart, crystal%gmet, ddb%gprim, mpert, natom, &
      Ifc%nrpt, qphnrm(1), qphon, crystal%rmet, ddb%rprim, Ifc%rpt, &
      Ifc%trans, crystal%ucvol, Ifc%wghatm, crystal%xred, zeff, qdrp_cart, Ifc%ewald_option, xmpi_comm_self, &
      dipquad = Ifc%dipquad, quadquad = Ifc%quadquad)

  else if (inp%ifcflag == 0) then
    ! Look for the information in the DDB
    rfphon(1:2)=1; rfelfd(1:2)=2; rfstrs(1:2)=0
    call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, inp%rfmeth)
    ! Copy the dynamical matrix in d2cart
    d2cart(:,1:msize)=ddb%val(:,:,iblok)
    ! Eventually impose the acoustic sum rule
    call asrq0%apply(natom, mpert, msize, crystal%xcart, d2cart)

  end if  ! end of the generation of the dynamical matrix at gamma.
  !***************************************************************

  ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
  call dfpt_phfrq(ddb%amu, displ, d2cart, eigval, eigvec, crystal%indsym, &
    mpert, msym, natom, nsym, ntypat, phfrq, qphnrm(1), qphon, &
    crystal%rprimd, inp%symdynmat, crystal%symrel, crystal%symafm, crystal%typat, crystal%ucvol)

  ! calculation of the oscillator strengths, mode effective charge and
  ! dielectric tensor, frequency dependent dielectric tensor (dieflag)
  ! and mode by mode decomposition of epsilon if dieflag == 3
  if (inp%dieflag /= 0) then
    !if (iblok_epsinf == 0) then
    if (epsinf(1, 1)==one .and. epsinf(2, 2)==one .and. epsinf(3, 3)==one) then
      write(msg, '(7a)') ch10, &
       ' The DDB file does not contain the derivatives w.r.t. electric field perturbation. ',ch10, &
       ' This is mandatory to calculate the dielectric constant, ',ch10, &
       ' Please check your DDB file or use dieflag = 0.',ch10
      ABI_ERROR(msg)
    end if

    write(msg, '(a, (80a), a)' ) ch10, ('=',ii = 1, 80), ch10
    call wrtout(units, msg)

    call ddb_diel(crystal, ddb%amu, inp, dielt_rlx, displ, d2cart, epsinf, fact_oscstr, &
      ab_out, lst, mpert, natom, 0, phfrq, comm, ana_ncid)
  end if

end if  ! dieflag!=0 or inp%nph2l /= 0


!**********************************************************************
! Non-linear response: electrooptic and Raman (q = Gamma, TO modes only)
!**********************************************************************
! if (inp%nlflag > 0) then
!   ABI_MALLOC(rsus, (3*natom, 3, 3))
! end if

 if (inp%nlflag == 1) then
   ! Raman susceptibilities for the 1st list (only TO  modes at q = Gamma)
   call ramansus(d2cart, dchide, dchidt, displ, mpert, natom, phfrq, qphon, qphnrm(1), rsus, crystal%ucvol)

   if (my_rank == master) then
     call nctk_defwrite_raman_terms(ana_ncid, natom, rsus, phfrq)
   end if

   ! EO coef:
   call electrooptic(dchide, inp%dieflag, epsinf, fact_oscstr, natom, phfrq, inp%prtmbm, rsus, crystal%ucvol)
end if  ! condition on nlflag

!**********************************************************************
! Calculation of properties associated to the second list of wv: nph2l /= 0
! (can include non-analyticities in the DM)
!**********************************************************************

 if (inp%nph2l /= 0) then

   write(msg, '(a, (80a), a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, ' Treat the second list of vectors ',ch10
   call wrtout(units, msg)

   if (my_rank == master) then
     iphl2 = 0
     call nctk_defwrite_nonana_terms(ana_ncid, iphl2, inp%nph2l, inp%qph2l, natom, phfrq, displ, "define")
     if (inp%nlflag == 1) then
       call nctk_defwrite_nonana_raman_terms(ana_ncid, iphl2, inp%nph2l, natom, rsus, "define")
     end if
   end if
!  Get the log of product of the square of the phonon frequencies without non-analyticities (q = 0)
!  For the Lyddane-Sachs-Teller relation, it is stored in lst(nph2+1)
   if (inp%dieflag /= 2 .and. inp%dieflag /= 0) then
     do ii = 4, 3*natom
       lst(inp%nph2l+1)=lst(inp%nph2l+1)+2*log(phfrq(ii))
     end do
   end if

   ! Examine every wavevector of this list
   do iphl2 = 1, inp%nph2l

     ! Initialisation of the phonon wavevector
     qphon(:,1)=inp%qph2l(:,iphl2)
     qphnrm(1)=inp%qnrml2(iphl2)

     !TODO: Quadrupole interactions need to be incorporated here (MR)

     ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
     ! for the second list of wv (can include non-analyticities if q /= 0)
     call dfpt_phfrq(ddb%amu, displ, d2cart, eigval, eigvec, crystal%indsym, &
       mpert, msym, natom, nsym, ntypat, phfrq, qphnrm(1), qphon, crystal%rprimd, inp%symdynmat, &
       crystal%symrel, crystal%symafm, crystal%typat, crystal%ucvol)

     ! Write the phonon frequencies for the second list of wv (can include non-analyticities if q /= 0)
     call dfpt_prtph(displ, inp%eivec, inp%enunit, ab_out, natom, phfrq, qphnrm(1), qphon)
     ! TODO: Mode effective charge could be printed here for LO modes (EB)

     if (my_rank == master) then
       ! Loop is not MPI-parallelized--> no need for MPI-IO API.
       call nctk_defwrite_nonana_terms(ana_ncid, iphl2, inp%nph2l, inp%qph2l, natom, phfrq, displ, "write")
     end if

     !  Get the log of product of the square of the phonon frequencies with non-analyticities (q-->0)
     !  For the Lyddane-Sachs-Teller relation
     ! The fourth mode should have positive frequency otherwise there is an instability: LST relationship should not be evaluated
     ! Isn't it tested somewhere else (i.e. stop of the code if there are imaginary freq.)? (EB)
     if (inp%dieflag /= 2 .and. inp%dieflag /= 0) then
       do ii = 4, 3*natom
         lst(iphl2)=lst(iphl2)+2*log(phfrq(ii))
       end do
     end if

     ! Write Raman susceptibilities for the 2nd list (can includes LO modes if q /= 0 0 0)
     if (inp%nlflag == 1) then
       call ramansus(d2cart, dchide, dchidt, displ, mpert, natom, phfrq, qphon, qphnrm(1), rsus, crystal%ucvol)
       if (my_rank == master) then
         call nctk_defwrite_nonana_raman_terms(ana_ncid, iphl2, inp%nph2l, natom, rsus, "write")
       end if
     end if  ! nlflag = 1 (Raman suscep for the 2nd list of wv.)
   end do  ! iphl2

   ! Lyddane-Sachs-Teller relation:
   if (inp%dieflag /= 2 .and. inp%dieflag /= 0) then
     call ddb_diel(crystal, ddb%amu, inp, dielt_rlx, displ, d2cart, epsinf, fact_oscstr, &
       ab_out, lst, mpert, natom, inp%nph2l, phfrq, comm, ana_ncid)
   end if
 end if  ! nph2l /= 0
 ! End of second list of wv stuff (nph2l /= 0)


 ABI_FREE(fact_oscstr)
 ABI_FREE(lst)
 if (inp%nlflag > 0) then
   ABI_FREE(rsus)
   ABI_FREE(dchide)
   ABI_FREE(dchidt)
 end if

!**********************************************************************
! Linear response with strain: elastic, piezo, etc
!**********************************************************************

 if (inp%instrflag /= 0) then
   ! Here treating the internal strain tensors at Gamma point
   write(msg, '(a, a, (80a), a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
    ' Calculation of the internal-strain  tensor',ch10
   call wrtout(units, msg)

   if (inp%instrflag == 1) then
     call wrtout(std_out, 'instrflag = 1, so extract the internal strain constant from the 2DTE')

     ! looking after the no. of blok that contains the internal strain tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, inp%rfmeth)
     if (iblok == 0) then
       ABI_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
     end if

     ! then print the internal stain tensor
     prt_internalstr = 2
     call ddb_internalstr(inp%asr, ddb%val, asrq0%d2asr, iblok, instrain, ab_out, mpert, natom, ddb%nblok, prt_internalstr)
   end if
 end if  ! internal strain

!**********************************************************************

 if (inp%elaflag /= 0) then
   ! here treating the elastic tensors at Gamma Point
   write(msg, '(a, a, (80a), a, a, a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
    ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
   call wrtout(units, msg)

   if (any(inp%elaflag == [1, 2, 3, 4, 5])) then
     call wrtout(std_out, 'so extract the elastic constant from the 2DTE')

     ! look after the blok no. that contains the stress tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=0

     call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, BLKTYP_d1E_xx)
     iblok_stress = iblok

     ! look after the blok no.iblok that contains the elastic tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     ! for both diagonal and shear parts
     call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, inp%rfmeth)
     if (iblok == 0) then
       ABI_ERROR("DDB file must contain both uniaxial and shear strain when elaflag != 0, Check your calculations")
     end if

     ! print the elastic tensor
     call ddb_elast(inp, crystal, ddb%val, compl, compl_clamped, compl_stress, asrq0%d2asr, &
       elast, elast_clamped, elast_stress, iblok, iblok_stress, &
       instrain, ab_out, mpert, natom, ddb%nblok, ana_ncid)
   end if
 end if  ! elastic tensors

!**********************************************************************

 if (inp%piezoflag /= 0 .or. inp%dieflag == 4 .or. inp%elaflag == 4) then
   ! here treating the piezoelectric tensor at Gamma Point
   write(msg, '(a, a, (80a), a, a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
   ' Calculation of the tensor related to piezoelectric effetc',ch10, &
   '  (Elastic indices in Voigt notation)',ch10
   call wrtout(units, msg)

   if (any(inp%piezoflag == [1, 2, 3, 4, 5, 6, 7]) .or. inp%dieflag == 4 .or. inp%elaflag == 4) then
     call wrtout(std_out, 'extract the piezoelectric constant from the 2DTE')

     ! looking for the gamma point block
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     ! for both diagonal and shear parts
     call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, inp%rfmeth)
     if (iblok == 0) then
       ABI_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
     end if

     ! then print out the piezoelectric constants
     call ddb_piezo(inp, ddb%val, dielt_rlx, elast, iblok, &
         & instrain, ab_out, mpert, natom, ddb%nblok, piezo, &
         & crystal%ucvol, ana_ncid)
   end if
 end if

!**********************************************************************
! Flexoelectric response
!**********************************************************************

 if (inp%flexoflag /= 0 ) then
   ! Here treating the flexoelectric tensor
   write(msg, '(a, a, (80a), a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
   ' Calculation of the tensors related to flexoelectric effect',ch10
   call wrtout(units, msg)

   ! Compute and print the contributions to the flexoelectric tensor
   call ddb_flexo(inp%asr, asrq0%d2asr, ddb, ddb_lw, ddb_hdr%ddb_version, crystal, &
       & filnam(3), inp%flexoflag, inp%prtvol, zeff)
 end if

!**********************************************************************

 ! Free memory
 ABI_FREE(displ)
 ABI_FREE(d2cart)
 ABI_FREE(eigval)
 ABI_FREE(eigvec)
 ABI_FREE(phfrq)
 ABI_FREE(zeff)
 ABI_FREE(qdrp_cart)
 ABI_FREE(instrain)

 50 continue

 call asrq0%free()
 call ifc%free()
 call crystal%free()
 call ddb%free()
 call anaddb_dtset_free(inp)
 call thermal_supercell_free(inp%ntemper, thm_scells)
 call ddb_lw%free()

 if (sum(abs(inp%thermal_supercell))>0 .and. inp%ifcflag == 1) then
   ABI_FREE(thm_scells)
 end if

 ! Close files
 if (iam_master) then
   NCF_CHECK(nf90_close(ana_ncid))
 end if

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
   ' anaddb : the run completed succesfully.'
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
 call abinit_doctor(filnam(2))

 call flush_unit(ab_out)
 call flush_unit(std_out)

 if (iam_master) close(ab_out)

 100 call xmpi_end()

 end program anaddb
!!***
