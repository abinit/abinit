!{\src2tex{textfont=tt}}
!!****p* ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (XG,DCA,JCC,CL,XW,GA,MR)
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
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,anaddb_dtset_free,anaddb_init
!!      asrq0_apply,asrq0_free,crystal_free,ddb_diel,ddb_elast,ddb_free
!!      ddb_from_file,ddb_hdr_free,ddb_hdr_open_read,ddb_internalstr
!!      ddb_interpolate,ddb_piezo,dfpt_phfrq,dfpt_prtph,dfpt_symph
!!      elast_ncwrite,electrooptic,elphon,flush_unit,gruns_anaddb,gtblk9,gtdyn9
!!      harmonic_thermo,herald,ifc_free,ifc_init,ifc_outphbtrap,ifc_print
!!      ifc_speedofsound,ifc_write,instrng,int2char4,inupper,invars9,isfile
!!      mkphbs,mkphdos,nctk_defwrite_nonana_terms,outvars_anaddb,phdos_free
!!      phdos_ncwrite,phdos_print,phdos_print_debye,phdos_print_msqd
!!      phdos_print_thermo,ramansus,relaxpol,thermal_supercell_free
!!      thermal_supercell_make,thermal_supercell_print,thmeig,timein,wrtout
!!      xmpi_bcast,xmpi_end,xmpi_init,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program anaddb

 use defs_basis
 use m_build_info
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
 use iso_c_binding
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

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
 use m_harmonic_thermo,only : harmonic_thermo
 use m_thmeig,         only : thmeig
 use m_raman,          only : ramansus, electrooptic
 use m_ddb_diel,       only : ddb_diel
 use m_relaxpol,       only : relaxpol
 use m_ddb_elast,      only : ddb_elast
 use m_ddb_piezo,      only : ddb_piezo
 use m_ddb_internalstr, only : ddb_internalstr
#ifdef MR_DEV
 use m_ddb_flexo
#endif

 implicit none

!Local variables-------------------------------
 integer :: msym !  msym =maximum number of symmetry elements in space group
!Define input and output unit numbers (some are defined in defs_basis -all should be there ...):
 integer,parameter :: ddbun=2,master=0 ! FIXME: these should not be reserved unit numbers!
 integer,parameter :: rftyp4=4
 integer :: comm,iatom,iblok,iblok_stress,iblok_tmp,idir,ii,index
 integer :: ierr,iphl2,lenstr,lwsym,mtyp,mpert,msize,natom
 integer :: nsym,ntypat,option,usepaw,nproc,my_rank,ana_ncid,prt_internalstr
 logical :: iam_master
 integer :: rfelfd(4),rfphon(4),rfstrs(4),ngqpt_coarse(3)
 integer :: count_wminmax(2)
 integer,allocatable :: d2flg(:)
 real(dp) :: etotal,tcpu,tcpui,twall,twalli !,cpu, wall, gflops
 real(dp) :: dielt(3,3)
 real(dp) :: compl(6,6),compl_clamped(6,6),compl_stress(6,6)
 real(dp) :: dielt_rlx(3,3),elast(6,6),elast_clamped(6,6),elast_stress(6,6)
 real(dp) :: epsinf(3,3),red_ptot(3),pel(3)
 real(dp) :: piezo(6,3),qphnrm(3),qphon(3,3),strten(6),tsec(2)
 real(dp) :: wminmax(2)
 real(dp),allocatable :: d2cart(:,:),dchide(:,:,:)
 real(dp),allocatable :: dchidt(:,:,:,:),displ(:),eigval(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:,:),fact_oscstr(:,:,:),instrain(:,:)
 real(dp),allocatable :: fred(:,:),lst(:),phfrq(:)
 real(dp),allocatable :: rsus(:,:,:)
 real(dp),allocatable :: zeff(:,:,:)
 real(dp),allocatable :: qdrp_cart(:,:,:,:)
 character(len=10) :: procstr
 character(len=24) :: codename
 character(len=24) :: start_datetime
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(7),elph_base_name,tmpfilename, phibz_prefix
 character(len=500) :: msg
 type(args_t) :: args
 type(anaddb_dataset_type) :: inp
 type(phonon_dos_type) :: Phdos
 type(ifc_type) :: Ifc,Ifc_coarse
 type(ddb_type) :: ddb
#ifdef MR_DEV
 type(ddb_type) :: ddb_lw
#endif
 type(ddb_hdr_type) :: ddb_hdr
 type(asrq0_t) :: asrq0
 type(crystal_t) :: Crystal
 type(supercell_type), allocatable :: thm_scells(:)
#ifdef HAVE_NETCDF
 integer :: phdos_ncid, ncerr
#endif

!******************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

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
 call abimem_init(args%abimem_level, limit_mb=args%abimem_limit_mb)
#endif

 ! Initialisation of the timing
 call timein(tcpui,twalli)

 if (iam_master) then
   codename='ANADDB'//repeat(' ',18)
   call herald(codename,abinit_version,std_out)
 end if

 start_datetime = asctime()

 ! Zero out all accumulators of time and init timers
 call timab(1, 0, tsec)

 ! Initialise the code: write heading, and read names of files.
 if (iam_master) call anaddb_init(filnam)
 call xmpi_bcast (filnam, master, comm, ierr)

 ! make log file for non-master procs
 if (.not. iam_master) then
   call int2char4(my_rank, procstr)
   ABI_CHECK((procstr(1:1)/='#'),'Bug: string length too short!')
   tmpfilename = trim(filnam(2)) // "_LOG_P" // trim(procstr)
   if (open_file(tmpfilename, msg, unit=std_out, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if
 end if

!******************************************************************

 ! Must read natom from the DDB before being able to allocate some arrays needed for invars9

 call ddb_hdr_open_read(ddb_hdr,filnam(3),ddbun,DDB_VERSION,comm=comm, dimonly=1)

 natom = ddb_hdr%natom
 ntypat = ddb_hdr%ntypat
 mtyp = ddb_hdr%mblktyp
 usepaw = ddb_hdr%usepaw

 call ddb_hdr_free(ddb_hdr)

 mpert=natom+MPERT_MAX
 msize=3*mpert*3*mpert; if (mtyp==3) msize=msize*3*mpert

 ! Read the input file, and store the information in a long string of characters
 ! strlen from defs_basis module
 option=1
 if (iam_master) then
   call instrng (filnam(1),lenstr,option,strlen,string)

   ! To make case-insensitive, map characters to upper case.
   call inupper(string(1:lenstr))
 end if

 call xmpi_bcast(string, master, comm, ierr)
 call xmpi_bcast(lenstr, master, comm, ierr)

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
   call isfile(tmpfilename,'new')
   if (open_file(tmpfilename,msg,unit=ab_out,form='formatted',status='new') /= 0) then
     MSG_ERROR(msg)
   end if
   rewind (unit=ab_out)
   call herald(codename,abinit_version,ab_out)

   ! Echo the inputs to console and main output file
   call outvars_anaddb(inp,std_out)
   call outvars_anaddb(inp,ab_out)
 else
   ab_out = dev_null
 end if

!******************************************************************

 ! Read the DDB information, also perform some checks, and symmetrize partially the DDB
 write(msg, '(a,a)' )' read the DDB information and perform some checks',ch10
 call wrtout([std_out, ab_out], msg)

 call ddb_from_file(ddb,filnam(3),inp%brav,natom,inp%natifc,inp%atifc,Crystal,comm, prtvol=inp%prtvol)
 nsym = Crystal%nsym

#ifdef MR_DEV
 ! A new ddb is necessary for the quadrupoles due to incompability of it with authomatic reshapes
 ! that ddb%val and ddb%flg experience when passed as arguments of some routines
 if (mtyp==33) then
   call ddb_lw_copy(ddb,ddb_lw,mpert,natom,ntypat)
 end if
#endif

 ! Acoustic Sum Rule
 ! In case the interatomic forces are not calculated, the
 ! ASR-correction (asrq0%d2asr) has to be determined here from the Dynamical matrix at Gamma.
 asrq0 = ddb%get_asrq0(inp%asr, inp%rfmeth, crystal%xcart)

 ! TODO: This is to maintain the previous behaviour in which all the arrays were initialized to zero.
 ! In the new version asrq0%d2asr is always computed if the Gamma block is present
 ! and this causes changes in [v5][t28]
 if (.not. (inp%ifcflag==0 .or. inp%instrflag/=0 .or. inp%elaflag/=0)) then
   asrq0%d2asr = zero
   if (asrq0%asr==3.or.asrq0%asr==4) then
     asrq0%singular = zero; asrq0%uinvers = zero; asrq0%vtinvers = zero
   end if
 end if

 ! Open the netcdf file that will contain the anaddb results
 ana_ncid = nctk_noid
 if (iam_master) then
#ifdef HAVE_NETCDF
   NCF_CHECK_MSG(nctk_open_create(ana_ncid, "anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
   ncerr = nctk_def_dims(ana_ncid, [ &
       nctkdim_t('number_of_atoms', natom), &
       nctkdim_t('natom3', 3 * natom), &
       nctkdim_t('number_of_phonon_modes', 3 * natom), &
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
#endif
 end if

 ! Calculation of Grunesein parameters.
 if (inp%gruns_nddbs /= 0) then
   call gruns_anaddb(inp, filnam(2), comm)
   goto 50
 end if

 ABI_MALLOC(instrain,(3*natom,6))
 ABI_MALLOC(d2cart,(2,msize))
 ABI_MALLOC(displ,(2*3*natom*3*natom))
 ABI_MALLOC(eigval,(3,natom))
 ABI_MALLOC(eigvec,(2,3,natom,3,natom))
 ABI_MALLOC(phfrq,(3*natom))
 ABI_MALLOC(zeff,(3,3,natom))
 ABI_MALLOC(qdrp_cart,(3,3,3,natom))
 ABI_MALLOC(lst,(inp%nph2l))

!**********************************************************************
!**********************************************************************

 ! Get Quadrupole tensor
 write(msg,'(2a,(80a),2a)') ch10,('=',ii=1,80)
 call wrtout([ab_out,std_out],msg,'COLL')
 lwsym=1
 iblok = ddb_lw%get_quadrupoles(crystal,lwsym,33,qdrp_cart)
 if ((inp%dipquad==1.or.inp%quadquad==1).and.iblok == 0) then
   call wrtout(std_out, "--- !WARNING")
   call wrtout(std_out, sjoin("- Cannot find Dynamical Quadrupoles tensor in DDB file:", filnam(3)))
   call wrtout(std_out, "  dipquad=1 or quadquad=1 requires the DDB file to include the corresponding long wave 3rd derivatives")
 end if

 ! Get Dielectric Tensor and Effective Charges
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblok = ddb%get_dielt_zeff(crystal,inp%rfmeth,inp%chneut,inp%selectz,dielt,zeff)

 ! Try to get dielt, in case just the DDE are present
 if (iblok == 0) then
   iblok_tmp = ddb%get_dielt(inp%rfmeth,dielt)
 end if

 !if (iblok == 0) then
 !  call wrtout(std_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", filnam(3)))
 !  call wrtout(std_out, "Values initialized with zeros")
 !else
 !  call wrtout(std_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", filnam(3)))
 !end if

 if (my_rank == master) then
#ifdef HAVE_NETCDF
   ncerr = nctk_def_arrays(ana_ncid, [&
   nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'),&
   nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")],&
   defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ana_ncid, [character(len=nctk_slen) :: &
       "asr", "chneut", "dipdip", "symdynmat"])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'emacro_cart'), dielt))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'becs_cart'), zeff))
   ncerr = nctk_write_iscalars(ana_ncid, [character(len=nctk_slen) :: &
     "asr", "chneut", "dipdip", "symdynmat"], &
     [inp%asr, inp%chneut, inp%dipdip, inp%symdynmat])
   NCF_CHECK(ncerr)
#endif
 end if

 ! Structural response at fixed polarization
 if (inp%polflag == 1) then
   ABI_MALLOC(d2flg, (msize))

   if(iblok/=0)then
     ! Save the second-order derivatives
     d2cart(1:2,1:msize) = ddb%val(1:2,1:msize,iblok)
     d2flg(1:msize) = ddb%flg(1:msize,iblok)

   else
     ! the gamma blok has not been found
     if (inp%relaxat==0 .and. inp%relaxstr==0) then
       ! The gamma blok is not needed
       d2cart(1:2,1:msize)=zero
       d2flg(1:msize)=1
     else
       ! There is a problem !
       write(msg, '(7a)' )&
        'The dynamical matrix at Gamma is needed, in order to perform ',ch10,&
        "relaxation at constant polarisation (Na Sai's method)",ch10,&
        'However, this was not found in the DDB.',ch10,&
        'Action: complete your DDB with the dynamical matrix at Gamma.'
       MSG_ERROR(msg)
     end if
   end if ! iblok not found

   ! Extract the block with the total energy
   if (ddb%get_etotal(etotal) == 0) then
     MSG_ERROR("DDB file does not contain GS etotal")
   end if

   ! Extract the block with the gradients
   ABI_MALLOC(fred, (3,natom))
   qphon(:,:) = zero; qphnrm(:) = zero
   rfphon(:) = 0; rfstrs(:) = 0; rfelfd(:) = 2
   if (inp%relaxat == 1) rfphon(:) = 1
   if (inp%relaxstr == 1) rfstrs(:) = 3

   call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp4)

   if (inp%relaxat == 1) then
     index = 0
     do iatom = 1, natom
       do idir = 1, 3
         index = index + 1
         fred(idir,iatom) = ddb%val(1,index,iblok)
       end do
     end do
   end if

   pel(1:3) = ddb%val(1,3*natom+4:3*natom+6,iblok)

   if (inp%relaxstr == 1) then
     index = 3*natom + 6
     do ii = 1, 6
       index = index + 1
       strten(ii) = ddb%val(1,index,iblok)
     end do
   end if

   ! when called from here red_ptot is not set! So set it to zero
   red_ptot(:)=zero

   call relaxpol(Crystal,d2flg,d2cart,etotal,fred,inp%iatfix,&
&   ab_out,inp%istrfix,mpert,msize,inp%natfix,natom,&
&   inp%nstrfix,pel,red_ptot,inp%relaxat,inp%relaxstr,&
&   strten,inp%targetpol,usepaw)

   ABI_FREE(fred)
   ABI_FREE(d2flg)
 end if

!***************************************************************************

 ! Compute non-linear optical susceptibilities and, if inp%nlflag < 3,
 ! First-order change in the linear dielectric susceptibility induced by an atomic displacement
 if (inp%nlflag > 0) then
   ABI_MALLOC(dchide, (3,3,3))
   ABI_MALLOC(dchidt, (natom,3,3,3))

   if (ddb%get_dchidet(inp%ramansr,inp%nlflag,dchide,dchidt) == 0) then
     MSG_ERROR("Cannot find block corresponding to non-linear optical susceptibilities in DDB file")
   end if

   ! Save to the netcdf
   if (my_rank == master) then
#ifdef HAVE_NETCDF
     ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchide", "dp", "three, three, three")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ana_ncid))
     NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchide"), dchide))

     ! dchidt only present if nlflag==1 or 2
     if (inp%nlflag < 3) then
       ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchidt", "dp", &
       "number_of_atoms, three, three, three")], defmode=.True.)
       NCF_CHECK(ncerr)
       NCF_CHECK(nctk_set_datamode(ana_ncid))
       NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchidt"), dchidt))
     end if
#endif
   end if

 end if ! nlflag

!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
 if (inp%ifcflag ==1) then
   ! ifc to be calculated for interpolation
   write(msg, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
    ' Calculation of the interatomic forces ',ch10
   call wrtout([std_out, ab_out], msg)

! TODO : check if this wrtout should be removed in latest merge 17 feb 2017
   call timein(tcpu,twall)
   write(msg, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout([std_out, ab_out], msg)

#ifdef MR_DEV
   if (any(inp%qrefine(:) > 1)) then
     ! Gaal-Nagy's algorithm in PRB 73 014117 [[cite:GaalNagy2006]]
     ! Build the IFCs using the coarse q-mesh.
     do ii = 1, 3
       ngqpt_coarse(ii) = inp%ngqpt(ii) / inp%qrefine(ii)
     end do
     call ifc_init(Ifc_coarse,Crystal,ddb,&
       inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,ngqpt_coarse,inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
       inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm,dipquad=inp%dipquad,quadquad=inp%quadquad)

     ! Now use the coarse q-mesh to fill the entries in dynmat(q)
     ! on the dense q-mesh that cannot be obtained from the DDB file.
     call ifc_init(Ifc,Crystal,ddb,&
      inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
      inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm,Ifc_coarse=Ifc_coarse,dipquad=inp%dipquad,quadquad=inp%quadquad)
     call Ifc_coarse%free()

   else
     call ifc_init(Ifc,Crystal,ddb,&
       inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
       inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm,dipquad=inp%dipquad,quadquad=inp%quadquad)
   end if
#else
   if (any(inp%qrefine(:) > 1)) then
     ! Gaal-Nagy's algorithm in PRB 73 014117 [[cite:GaalNagy2006]]
     ! Build the IFCs using the coarse q-mesh.
     do ii = 1, 3
       ngqpt_coarse(ii) = inp%ngqpt(ii) / inp%qrefine(ii)
     end do
     call ifc_init(Ifc_coarse,Crystal,ddb,&
       inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,ngqpt_coarse,inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
       inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm)

     ! Now use the coarse q-mesh to fill the entries in dynmat(q)
     ! on the dense q-mesh that cannot be obtained from the DDB file.
     call ifc_init(Ifc,Crystal,ddb,&
      inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
      inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm,Ifc_coarse=Ifc_coarse)
     call Ifc_coarse%free()

   else
     call ifc_init(Ifc,Crystal,ddb,&
       inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,qdrp_cart,&
       inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit,comm)
   end if
#endif

   call ifc%print(unit=std_out)

   ! Compute speed of sound.
   if (inp%vs_qrad_tolkms(1) > zero) then
     call ifc%speedofsound(crystal, inp%vs_qrad_tolkms, ana_ncid, comm)
   end if

   ! Print analysis of the real-space interatomic force constants
   ! TODO: ifc_out should not have side effects
   if (my_rank == master .and. inp%ifcout/=0) then
     call ifc%write(inp%ifcana,inp%atifc,inp%ifcout,inp%prt_ifc,ana_ncid)
   end if
 end if

!**********************************************************************

!Electron-phonon section
 if (inp%elphflag == 1) then
   call elphon(inp,Crystal,Ifc,filnam,comm)
 end if

!**********************************************************************

 if (sum(abs(inp%thermal_supercell))>0 .and. inp%ifcflag==1) then
   ABI_MALLOC(thm_scells, (inp%ntemper))
   call zacharias_supercell_make(Crystal, Ifc, inp%ntemper, inp%thermal_supercell, inp%tempermin, inp%temperinc, thm_scells)
   call zacharias_supercell_print(filnam(2), inp%ntemper, inp%tempermin, inp%temperinc, thm_scells)
 end if

 ! Phonon density of states calculation, Start if interatomic forces have been calculated
 if (inp%ifcflag==1 .and. any(inp%prtdos==[1, 2])) then
   write(msg,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,' Calculation of phonon density of states ',ch10
   call wrtout([std_out, ab_out], msg)

   ! Only 1 shift in q-mesh
   wminmax = zero
   phibz_prefix = ""
   !phibz_prefix = "freq_displ" ! Uncomment this line to activate output of PHIBZ
   do
     call mkphdos(Phdos, Crystal, Ifc, inp%prtdos, inp%dosdeltae, inp%dossmear, inp%ng2qpt, 1, inp%q2shft, &
      phibz_prefix, wminmax, count_wminmax, comm)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(msg, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
   end do

   if (iam_master) then
     call phdos%print_msqd(filnam(2), inp%ntemper, inp%tempermin, inp%temperinc)
     call phdos%print(strcat(filnam(2), "_PHDOS"))
     call phdos%print_debye(Crystal%ucvol)
     call phdos%print_thermo(strcat(filnam(2), "_THERMO"), inp%ntemper, inp%tempermin, inp%temperinc)

#ifdef HAVE_NETCDF
     ncerr = nctk_open_create(phdos_ncid, strcat(filnam(2), "_PHDOS.nc"), xmpi_comm_self)
     NCF_CHECK_MSG(ncerr, "Creating PHDOS.nc file")
     NCF_CHECK(Crystal%ncwrite(phdos_ncid))
     call phdos%ncwrite(phdos_ncid)
     NCF_CHECK(nf90_close(phdos_ncid))
#endif
   end if

   call phdos%free()
 end if

 if (iam_master .and. inp%ifcflag==1 .and. inp%outboltztrap==1) then
   call ifc%outphbtrap(Crystal, inp%ng2qpt, 1, inp%q2shft, filnam(2))
 end if

 ! Phonon density of states and thermodynamical properties calculation
 ! Start if interatomic forces and thermal flags are on
 if (inp%ifcflag==1 .and. any(inp%thmflag==[1,2])) then

   write(msg, '(a,(80a),a,a,a,a,a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
    ' Calculation of phonon density of states, ',ch10,&
    '    thermodynamical properties, ',ch10,&
    '    and Debye-Waller factors.',ch10
   call wrtout([std_out, ab_out], msg)

   if (inp%thmflag==1) then
     call harmonic_thermo(Ifc,Crystal,ddb%amu,inp,ab_out,filnam(2),comm)

   else if (inp%thmflag==2) then
     write(msg, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,' Entering thm9 routine with thmflag=2 ',ch10
     call wrtout([std_out, ab_out], msg)
     call harmonic_thermo(Ifc,Crystal,ddb%amu,inp,ab_out,filnam(2),comm,thmflag=inp%thmflag)
   end if
 end if

!**********************************************************************

 ! Now treat the first list of vectors (without non-analyticities)
 call mkphbs(Ifc,Crystal,inp,ddb,asrq0,filnam(2),comm)

!***********************************************************************

 if (inp%prtddb == 1 .and. inp%ifcflag == 1) then
   ! Interpolate the DDB onto the first list of vectors and write the file.
   call ddb_hdr_open_read(ddb_hdr,filnam(3),ddbun,DDB_VERSION)
   close(ddbun)
   call ddb_interpolate(Ifc,Crystal,inp,ddb,ddb_hdr,asrq0,filnam(2),comm)
   call ddb_hdr_free(ddb_hdr)
 end if

!***********************************************************************

 if (inp%thmflag>=3 .and. inp%thmflag<=8) then
   elph_base_name=trim(filnam(2))//"_ep"
   call thmeig(inp,ddb,Crystal,elph_base_name,filnam(5),ddbun,ab_out,natom,mpert,msize,asrq0%d2asr,comm)
 end if

!**********************************************************************

 ! Now treat the second list of vectors (only at the Gamma point,
 ! but can include non-analyticities), as well as the frequency-dependent dielectric tensor

 if (inp%nlflag > 0) then
   ABI_MALLOC(rsus, (3*natom,3,3))
 end if
 ABI_MALLOC(fact_oscstr, (2,3,3*natom))

 if (inp%nph2l/=0 .or. inp%dieflag==1) then

   write(msg, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,' Treat the second list of vectors ',ch10
   call wrtout([std_out, ab_out], msg)

   ! Before examining every direction or the dielectric tensor, generates the dynamical matrix at gamma
   qphon(:,1)=zero; qphnrm(1)=zero

   ! Generation of the dynamical matrix in cartesian coordinates
   if (inp%ifcflag==1) then

     ! Get d2cart using the interatomic forces and the
     ! long-range coulomb interaction through Ewald summation
#ifdef MR_DEV
     call gtdyn9(ddb%acell,Ifc%atmfrc,dielt,Ifc%dipdip,&
       Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,mpert,natom,&
       Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,&
       Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,zeff,qdrp_cart,Ifc%ewald_option,xmpi_comm_self,&
       dipquad=Ifc%dipquad,quadquad=Ifc%quadquad)
#else
     call gtdyn9(ddb%acell,Ifc%atmfrc,dielt,inp%dipdip,&
       Ifc%dyewq0,d2cart,Crystal%gmet,ddb%gprim,mpert,natom,&
       Ifc%nrpt,qphnrm(1),qphon,Crystal%rmet,ddb%rprim,Ifc%rpt,&
       Ifc%trans,Crystal%ucvol,Ifc%wghatm,Crystal%xred,zeff,qdrp_cart,Ifc%ewald_option,xmpi_comm_self)
#endif

   else if (inp%ifcflag==0) then

     ! Look after the information in the DDB
     rfphon(1:2)=1; rfelfd(1:2)=2; rfstrs(1:2)=0
     call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,inp%rfmeth)

     ! Copy the dynamical matrix in d2cart
     d2cart(:,1:msize)=ddb%val(:,:,iblok)

     ! Eventually impose the acoustic sum rule
     call asrq0%apply(natom, mpert, msize, crystal%xcart, d2cart)
   end if ! end of the generation of the dynamical matrix at gamma.

   if (inp%nph2l/=0) then

     if (my_rank == master) then
#ifdef HAVE_NETCDF
       iphl2 = 0
       call nctk_defwrite_nonana_terms(ana_ncid, iphl2, inp%nph2l, inp%qph2l, natom, phfrq, displ, "define")
#endif
     end if

     ! Examine every wavevector of this list
     do iphl2=1,inp%nph2l

       ! Initialisation of the phonon wavevector
       qphon(:,1)=inp%qph2l(:,iphl2)
       qphnrm(1)=inp%qnrml2(iphl2)

       !TODO: Quadrupole interactions need to be incorporated here (MR)
       ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
       call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
         mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,Crystal%rprimd,inp%symdynmat,&
         Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

       ! Write the phonon frequencies
       call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

       if (my_rank == master) then
#ifdef HAVE_NETCDF
         ! Loop is not MPI-parallelized --> no need for MPI-IO API.
         call nctk_defwrite_nonana_terms(ana_ncid, iphl2, inp%nph2l, inp%qph2l, natom, phfrq, displ, "write")
#endif
       end if

       ! Determine the symmetries of the phonon modes at Gamma
       if (sum(abs(qphon(:,1)))<DDB_QTOL) then
         call dfpt_symph(ab_out,ddb%acell,eigvec,Crystal%indsym,natom,nsym,phfrq,ddb%rprim,Crystal%symrel)
       end if

       ! Write Raman susceptibilities
       if (inp%nlflag == 1) then
         call ramansus(d2cart,dchide,dchidt,displ,mpert,natom,phfrq,qphon,qphnrm(1),rsus,Crystal%ucvol)
       end if

       ! Prepare the evaluation of the Lyddane-Sachs-Teller relation
       if(inp%dieflag==1 .and. natom>1)then
         lst(iphl2)=zero
         ! The fourth mode should have positive frequency, otherwise,
         ! there is an instability, and the LST relationship should not be evaluated
         if(phfrq(4)>tol6)then
           do ii=4,3*natom
             lst(iphl2)=lst(iphl2)+2*log(phfrq(ii))
           end do
         end if
       end if

     end do ! iphl2
   end if ! nph2l/=0

   ! The frequency-dependent dielectric tensor (and oscillator strength).
   if (inp%dieflag==1)then
     write(msg, '(6a)' )&
      ' the frequency-dependent dielectric tensor (and also once more',ch10,&
      ' the phonons at gamma - without non-analytic part )',ch10,ch10,&
      ' The frequency-dependent dielectric tensor'
     call wrtout(std_out, msg)

     ! Initialisation of the phonon wavevector
     qphon(:,1)=zero; qphnrm(1)=zero

     ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
     call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
       mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
       Crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)

     ! Write the phonon frequencies (not to ab_out, however)
     call dfpt_prtph(displ,0,inp%enunit,-1,natom,phfrq,qphnrm(1),qphon)

     ! Evaluation of the oscillator strengths and frequency-dependent dielectric tensor.
     call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
       ab_out,lst,mpert,natom,inp%nph2l,phfrq,comm,ana_ncid)
     ! write(std_out,*)'after ddb_diel, dielt_rlx(:,:)=',dielt_rlx(:,:)
   end if

   ! If the electronic dielectric tensor only is needed...
   if (inp%dieflag==2.or.inp%dieflag==3.or. inp%dieflag==4) then
     ! Everything is already in place...
     call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
       ab_out,lst,mpert,natom,inp%nph2l,phfrq,comm,ana_ncid)
   end if

 end if ! either nph2l/=0  or  dieflag==1

!**********************************************************************

!In case nph2l was equal to 0, the electronic dielectric tensor has to be computed independently.

 if (inp%dieflag==2 .and. inp%nph2l==0) then
   call wrtout(std_out, 'nph2l=0, so compute the electronic dielectric tensor independently')
   ! Look after the second derivative matrix at gamma in the DDB
   ! Note that the information on the dielectric tensor is completely
   ! independent of the interatomic force constant calculation
   qphon(:,1)=zero; qphnrm(1)=zero
   rfphon(1:2)=0; rfelfd(1:2)=2; rfstrs(1:2)=0
   call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,inp%rfmeth)
   if (iblok == 0) then
     MSG_ERROR("DDB file must contain both derivatives wrt electric field, Check your calculations")
   end if

   d2cart(:,1:msize)=ddb%val(:,:,iblok)

   ! Print the electronic dielectric tensor
   call ddb_diel(Crystal,ddb%amu,inp,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
   ab_out,lst,mpert,natom,inp%nph2l,phfrq,comm,ana_ncid)
 end if

!**********************************************************************

 if (inp%nlflag == 1) then
   ! Compute the electrooptic tensor
   ! In case dieflag = 2, recompute phonon frequencies and eigenvectors without non-analyticity
   if (inp%dieflag == 2) then
     qphon(:,1)=zero; qphnrm(1)=zero
     call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,Crystal%indsym,&
       mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,&
       Crystal%rprimd,inp%symdynmat,Crystal%symrel,Crystal%symafm,Crystal%typat,Crystal%ucvol)
   end if

   rsus = zero
   call ramansus(d2cart,dchide,dchidt,displ,mpert,natom,phfrq(1),qphon,qphnrm(1),rsus,Crystal%ucvol)

   call electrooptic(dchide,inp%dieflag,epsinf,fact_oscstr,natom,phfrq,inp%prtmbm,rsus,Crystal%ucvol)
 end if ! condition on nlflag

 ABI_FREE(fact_oscstr)
 if (inp%nlflag > 0) then
   ABI_FREE(dchide)
   ABI_FREE(rsus)
   ABI_FREE(dchidt)
 end if

!**********************************************************************

 if (inp%instrflag/=0) then
   ! Here treating the internal strain tensors at Gamma point
   write(msg, '(a,a,(80a),a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
    ' Calculation of the internal-strain  tensor',ch10
   call wrtout([std_out, ab_out], msg)

   if (inp%instrflag==1) then
     call wrtout(std_out, 'instrflag=1, so extract the internal strain constant from the 2DTE')

     ! looking after the no. of blok that contains the internal strain tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,inp%rfmeth)
     if (iblok == 0) then
       MSG_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
     end if

     ! then print the internal stain tensor
     prt_internalstr=2
     call ddb_internalstr(inp%asr,ddb%val,asrq0%d2asr,iblok,instrain,ab_out,mpert,natom,ddb%nblok,prt_internalstr)
   end if
 end if ! internal strain

!**********************************************************************

 if (inp%elaflag/=0) then
   ! here treating the elastic tensors at Gamma Point
   write(msg, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
    ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
   call wrtout([std_out, ab_out], msg)

   if (any(inp%elaflag == [1,2,3,4,5])) then
     call wrtout(std_out, 'so extract the elastic constant from the 2DTE')

     ! look after the blok no. that contains the stress tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=0

     call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp4)
     iblok_stress=iblok

     ! look after the blok no.iblok that contains the elastic tensor
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     ! for both diagonal and shear parts
     call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,inp%rfmeth)
     if (iblok == 0) then
       MSG_ERROR("DDB file must contain both uniaxial and shear strain when elaflag != 0, Check your calculations")
     end if

     ! print the elastic tensor
     call ddb_elast(inp,crystal,ddb%val,compl,compl_clamped,compl_stress,asrq0%d2asr,&
       elast,elast_clamped,elast_stress,iblok,iblok_stress,&
       instrain,ab_out,mpert,natom,ddb%nblok,ana_ncid)
   end if
 end if ! elastic tensors

!**********************************************************************

 if (inp%piezoflag/=0 .or. inp%dieflag==4 .or. inp%elaflag==4) then
   ! here treating the piezoelectric tensor at Gamma Point
   write(msg, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
   ' Calculation of the tensor related to piezoelectric effect',ch10,&
   '  (Elastic indices in Voigt notation)',ch10
   call wrtout([std_out, ab_out], msg)

   if (any(inp%piezoflag == [1,2,3,4,5,6,7]) .or. inp%dieflag==4 .or.inp%elaflag==4) then
     call wrtout(std_out, 'extract the piezoelectric constant from the 2DTE')

     ! looking for the gamma point block
     qphon(:,1)=zero; qphnrm(1)=zero
     rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

     ! for both diagonal and shear parts
     call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,inp%rfmeth)
     if (iblok == 0) then
       MSG_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
     end if

     ! then print out the piezoelectric constants
     call ddb_piezo(inp,ddb%val,dielt_rlx,elast,iblok,instrain,ab_out,mpert,natom,ddb%nblok,piezo,Crystal%ucvol,ana_ncid)
   end if
 end if

!**********************************************************************

#ifdef MR_DEV
 if (inp%flexoflag/=0 ) then
   ! Here treating the flexoelectric tensor
   write(msg, '(a,a,(80a),a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
   ' Calculation of the tensors related to flexoelectric effect',ch10
   call wrtout([std_out, ab_out], msg)

   ! Compute and print the contributions to the flexoelectric tensor
   call ddb_flexo(inp%asr,asrq0%d2asr,ddb,ddb_lw,crystal,filnam(3),inp%flexoflag,zeff)
 end if
#endif

!**********************************************************************

 ! Free memory
 ABI_FREE(displ)
 ABI_FREE(d2cart)
 ABI_FREE(eigval)
 ABI_FREE(eigvec)
 ABI_FREE(lst)
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
#ifdef MR_DEV
 call ddb_lw%free()
#endif

 if (sum(abs(inp%thermal_supercell))>0 .and. inp%ifcflag==1) then
   ABI_FREE(thm_scells)
 end if

 ! Close files
 if (iam_master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ana_ncid))
#endif
 end if

 call timein(tcpu,twall)
 tsec(1)=tcpu-tcpui; tsec(2)=twall-twalli
 write(msg, '(a,i4,a,f13.1,a,f13.1)' )' Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 call wrtout(std_out, msg)

 if (iam_master) then
   write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' )'-',ch10,&
    '- Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if

 call xmpi_sum(tsec,comm,ierr)

 write(msg, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
  ('=',ii=1,80),ch10,ch10,&
   '+Total cpu time',tsec(1),'  and wall time',tsec(2),' sec',ch10,ch10,&
   ' anaddb : the run completed succesfully.'
 call wrtout([std_out, ab_out], msg)

 if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
   write(std_out,"(a)")"--- !FinalSummary"
   write(std_out,"(a)")"program: anaddb"
   write(std_out,"(2a)")"version: ",trim(abinit_version)
   write(std_out,"(2a)")"start_datetime: ",start_datetime
   write(std_out,"(2a)")"end_datetime: ",asctime()
   write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
   write(std_out,"(a)")"..."
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
