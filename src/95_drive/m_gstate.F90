!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gstate
!! NAME
!!  m_gstate
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ, MB, DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_gstate

 use defs_basis
 use defs_rectypes
 use m_errors
 use m_xmpi
 use m_abicore
 use libxc_functionals
 use m_exit
 use m_crystal
 use m_scf_history
 use m_abimover
 use m_wffile
 use m_rec
 use m_efield
 use m_ddb
 use m_bandfft_kpt
 use m_invovl
 use m_gemm_nonlop
 use m_wfk
 use m_nctk
 use m_hdr
 use m_ebands
 use m_dtfil

 use defs_datatypes,     only : pseudopotential_type, ebands_t
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_symtk,            only : matr3inv
 use m_io_tools,         only : open_file
 use m_occ,              only : newocc, getnel
 use m_ddb_hdr,          only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_fstrings,         only : strcat, sjoin
 use m_geometry,         only : fixsym, mkradim, metric
 use m_kpts,             only : tetra_from_kptrlatt
 use m_kg,               only : kpgio, getph
 use m_fft,              only : fourdp
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,          only : pawcprj_type,pawcprj_free,pawcprj_alloc, pawcprj_getdim
 use m_pawfgr,           only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_abi2big,          only : wvl_occ_abi2big, wvl_setngfft, wvl_setBoxGeometry
 use m_energies,         only : energies_type, energies_init
 use m_args_gs,          only : args_gs_type
 use m_results_gs,       only : results_gs_type
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_copy, pawrhoij_free
 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,print_sc_dmft,paw_dmft_type,readocc_dmft
 use m_paw_sphharm,      only : setsym_ylm
 use m_paw_occupancies,  only : initrhoij
 use m_paw_init,         only : pawinit,paw_gencond
 use m_paw_correlations, only : pawpuxinit
 use m_orbmag,           only : initorbmag,destroy_orbmag,orbmag_type
 use m_paw_uj,           only : pawuj_ini,pawuj_free,pawuj_det, macro_uj_type
 use m_data4entropyDMFT, only : data4entropyDMFT_t, data4entropyDMFT_init, data4entropyDMFT_destroy
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron, &
                                electronpositron_calctype
 use m_scfcv,            only : scfcv_t, scfcv_init, scfcv_destroy, scfcv_run
 use m_dtfil,            only : dtfil_init_time
 use m_jellium,          only : jellium
 use m_iowf,             only : outwf, outresid
 use m_outqmc,           only : outqmc
 use m_ioarr,            only : ioarr,read_rhor
 use m_inwffil,          only : inwffil
 use m_spacepar,         only : setsym
 use m_mkrho,            only : mkrho, initro, prtrhomxmn
 use m_initylmg,         only : initylmg
 use m_pspini,           only : pspini
 use m_mover,            only : mover
 use m_mpinfo,           only : proc_distrb_cycle
 use m_common,           only : setup1, prteigrs, prtene
 use m_fourier_interpol, only : transgrid
 use m_psolver,          only : psolver_kernel
 use m_paw2wvl,          only : paw2wvl, wvl_paw_free
 use m_berryphase_new,   only : init_e_field_vars,prtefield
 use m_wvl_wfs,          only : wvl_wfs_set, wvl_wfs_free, wvl_wfs_lr_copy
 use m_wvl_rho,          only : wvl_initro, wvl_mkrho
 use m_wvl_descr_psp,    only : wvl_descr_psp_set, wvl_descr_free, wvl_descr_atoms_set, wvl_descr_atoms_set_sym
 use m_wvl_denspot,      only : wvl_denspot_set, wvl_denspot_free
 use m_wvl_projectors,   only : wvl_projectors_set, wvl_projectors_free

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

 use defs_wvltypes,      only : wvl_data,coulomb_operator,wvl_wf_type
#if defined HAVE_BIGDFT
 use BigDFT_API,         only : wvl_timing => timing,xc_init,xc_end,XC_MIXED,XC_ABINIT,&
                                local_potential_dimensions,nullify_gaussian_basis, &
                                copy_coulomb_operator,deallocate_coulomb_operator
#else
 use defs_wvltypes,      only : coulomb_operator
#endif

#if defined HAVE_LOTF
 use defs_param_lotf,    only : lotfparam_init
#endif

 implicit none

 private
!!***

 public :: gstate
!!***

contains
!!***

!!****f* m_gstate/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! INPUTS
!!  args_gs<type(args_gs_type)>=various input arguments for the GS calculation
!!                              Possibly different from dtset
!!  codvsn=code version
!!  cpui=initial CPU time
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =number of k points treated by this processor (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  initialized= 0 for the first GS calculation (not initialized), else 1
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  vel(3,natom)=value of velocity
!!  vel_cell(3,3)=value of cell parameters velocity
!!  wvl <type(wvl_data)>=all wavelets data
!!  xred(3,natom) = reduced atomic coordinates
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      gstateimg
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine gstate(args_gs,acell,codvsn,cpui,dtfil,dtset,iexit,initialized,&
&                 mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,&
&                 psps,results_gs,rprim,scf_history,vel,vel_cell,wvl,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit,initialized
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(args_gs_type),intent(in) :: args_gs
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),target,intent(inout) :: scf_history
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: acell(3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom),vel_cell(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=101,response=0 ,cplex1=1, master=0
 integer :: ndtpawuj=0  ! Cannot use parameter because scfargs points to this! Have to get rid of pointers to scalars!
#if defined HAVE_BIGDFT
 integer :: icoulomb
#endif
 integer :: accessfil,ask_accurate,bantot,choice,comm_psp,fform
 integer :: gnt_option,gscase,iatom,idir,ierr,ii,indx,jj,kk,ios,itypat
 integer :: ixfh,izero,mband_cprj,mcg,mcprj,me,mgfftf,mpert,msize,mu,my_natom,my_nspinor
 integer :: nblok,ncpgr,nfftf,nfftot,npwmin
 integer :: openexit,option,optorth,psp_gencond,conv_retcode
 integer :: pwind_alloc,rdwrpaw,comm,tim_mkrho,use_sc_dmft
 integer :: cnt,spin,band,ikpt,usecg,usecprj,ylm_option
 real(dp) :: cpus,ecore,ecut_eff,ecutdg_eff,etot,fermie
 real(dp) :: gsqcut_eff,gsqcut_shp,gsqcutc_eff,hyb_range_fock,residm,ucvol
 logical :: read_wf_or_den,has_to_init,call_pawinit,write_wfk
 logical :: is_dfpt=.false.,wvlbigdft=.false.,wvl_debug=.false.
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt,filnam,wfkfull_path
 real(dp) :: fatvshift
 type(crystal_t) :: cryst
 type(ebands_t) :: bstruct,ebands
 type(efield_type) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type) :: hdr,hdr_den
 type(macro_uj_type) :: dtpawuj(0)
 type(orbmag_type) :: dtorbmag
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(recursion_type) ::rec_set
 type(wffile_type) :: wff1,wffnew,wffnow
 type(ab_xfh_type) :: ab_xfh
 type(ddb_type) :: ddb
 type(ddb_hdr_type) :: ddb_hdr
 type(scfcv_t) :: scfcv_args
!arrays
 integer :: ngfft(18),ngfftf(18)
 integer,allocatable :: atindx(:),atindx1(:),indsym(:,:,:),dimcprj_srt(:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),symrec(:,:,:)
 integer,allocatable,target :: npwarr(:)
 integer,pointer :: npwarr_(:),pwind(:,:,:)
 real(dp) :: efield_band(3),gmet(3,3),gmet_for_kg(3,3),gprimd(3,3),gprimd_for_kg(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_for_kg(3,3),tsec(2)
 real(dp),allocatable :: doccde(:)
 real(dp),allocatable :: ph1df(:,:),phnons(:,:,:),resid(:),rhowfg(:,:)
 real(dp),allocatable :: rhowfr(:,:),spinat_dum(:,:),start(:,:),work(:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: cg(:,:),eigen(:),pwnsfac(:,:),rhog(:,:),rhor(:,:)
 real(dp),pointer :: taug(:,:),taur(:,:),xred_old(:,:)
 type(pawrhoij_type),pointer :: pawrhoij(:)
 type(coulomb_operator) :: kernel_dummy
 type(pawcprj_type),allocatable :: cprj(:,:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(32,1,tsec)
 call timab(33,3,tsec)

!###########################################################
!### 01. Initializations XML, MPI, WVL, etc

!Init MPI data
 comm=mpi_enreg%comm_cell; me=xmpi_comm_rank(comm)

!Set up MPI information from the dataset
 my_natom=mpi_enreg%my_natom

!Set up information when wavelets are in use
 if (dtset%usewvl == 1) then

!  If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
   wvlbigdft=(dtset%wvl_bigdft_comp==1)

!  Default value, to be set-up elsewhere.
   wvl%descr%h(:) = dtset%wvl_hgrid

#if defined HAVE_BIGDFT
   wvl%descr%paw%usepaw=psps%usepaw
   wvl%descr%paw%natom=dtset%natom
#endif

!  We set the atom-related internal wvl variables.
   call wvl_descr_atoms_set(acell, dtset%icoulomb, dtset%natom, &
&   dtset%ntypat, dtset%typat, wvl%descr)
   if(dtset%usepaw==0) then
!    nullify PAW proj_G in NC case:
#if defined HAVE_BIGDFT
     ABI_DATATYPE_ALLOCATE(wvl%projectors%G,(dtset%ntypat))
     do itypat=1,dtset%ntypat
       call nullify_gaussian_basis(wvl%projectors%G(itypat))
     end do
#endif
   end if

   wvl%descr%exctxpar = "OP2P"
 end if

 if (me == master .and. dtset%prtxml == 1) then
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  <dataSet>'
 end if

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,' gstate : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!###########################################################
!### 02. Calls setup1, kpgio, initylmg

 ecore=zero
 results_gs%pel(1:3)   =zero
 results_gs%grchempottn(:,:)=zero
 results_gs%grewtn(:,:)=zero
!MT Feb 2012: I dont know why but grvdw has to be allocated
!when using BigDFT to ensure success on inca_gcc44_sdebug
 if ((dtset%vdw_xc>=5.and.dtset%vdw_xc<=7).or.dtset%usewvl==1) then
   results_gs%ngrvdw=dtset%natom
   if (allocated(results_gs%grvdw)) then
     ABI_DEALLOCATE(results_gs%grvdw)
   end if
   ABI_ALLOCATE(results_gs%grvdw,(3,dtset%natom))
   results_gs%grvdw(:,:)=zero
 end if
 call energies_init(results_gs%energies)

!Set up for iterations
 call setup1(acell,bantot,dtset,&
  ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
  ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
  response,rmet,rprim,rprimd,ucvol,psps%usepaw)

!In some cases (e.g. getcell/=0), the plane wave vectors have
! to be generated from the original simulation cell
 rprimd_for_kg=rprimd
 if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=args_gs%rprimd_orig
 call matr3inv(rprimd_for_kg,gprimd_for_kg)
 gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet_for_kg,dtset%istwfk,kg, &
&   dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&   dtset%mpw,npwarr,npwtot,dtset%nsppol)
   call bandfft_kpt_init1(bandfft_kpt,dtset%istwfk,kg,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol)
 else
   ABI_ALLOCATE(kg,(0,0))
   npwarr(:) = 0
   npwtot(:) = 0
 end if

 if(dtset%wfoptalg == 1 .and. psps%usepaw == 1) then
   call init_invovl(dtset%nkpt)
 end if

 if(dtset%use_gemm_nonlop == 1 .and. dtset%use_gpu_cuda/=1) then
   ! set global variable
   gemm_nonlop_use_gemm = .true.
   call init_gemm_nonlop(dtset%nkpt)
 else
   gemm_nonlop_use_gemm = .false.
 end if

!Set up the Ylm for each k point
 if ( dtset%tfkinfunc /= 2) then
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     ylm_option=0
     if (dtset%prtstm==0.and.dtset%iscf>0.and.dtset%positron/=1) ylm_option=1 ! compute gradients of YLM
     if (dtset%berryopt==4 .and. dtset%optstress /= 0 .and. psps%usepaw==1) ylm_option = 1 ! compute gradients of YLM
     if ((dtset%orbmag.NE.0) .AND. (psps%usepaw==1)) ylm_option = 1 ! compute gradients of YLM
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,ylm_option,rprimd,ylm,ylmgr)
  end if
 else
   ABI_ALLOCATE(ylm,(0,0))
   ABI_ALLOCATE(ylmgr,(0,0,0))
 end if

!SCF history management (allocate it at first call)
 if (initialized==0) then
!  This call has to be done before any use of SCF history
   usecg=0
   if(dtset%extrapwf>0 .or. dtset%imgwfstor==1)usecg=1
   call scf_history_init(dtset,mpi_enreg,usecg,scf_history)
 end if
 has_to_init=(initialized==0.or.scf_history%history_size<0)

 call timab(33,2,tsec)
 call timab(701,3,tsec)

!###########################################################
!### 03. Calls pspini

!Open and read pseudopotential files
 comm_psp=mpi_enreg%comm_cell;if (dtset%usewvl==1) comm_psp=mpi_enreg%comm_wvl
 if (dtset%nimage>1) psps%mixalch(:,:)=args_gs%mixalch(:,:) ! mixalch can evolve for some image algos
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,pawrad,pawtab,psps,rprimd,comm_mpi=comm_psp)
 call timab(701,2,tsec)
 call timab(33,3,tsec)

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy as the ion-ion interaction.
 if (dtset%icoulomb == 1) ecore = zero

!WVL - Now that psp data are available, we compute rprimd, acell... from the atomic positions.
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_set(trim(dtfil%filnam_ds(3))//"_OCCUP",dtset%nsppol,psps,dtset%spinat,wvl%descr)
   call wvl_setBoxGeometry(dtset%prtvol, psps%gth_params%radii_cf, rprimd, xred, &
&   wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult)
   call mkradim(acell,rprim,rprimd)
   rprimd_for_kg=rprimd
   call wvl_denspot_set(wvl%den, psps%gth_params, dtset%ixc, dtset%natom, dtset%nsppol, rprimd, &
&   wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, mpi_enreg%comm_wvl, xred)
!  TODO: to be moved in a routine.
#if defined HAVE_BIGDFT
   if (wvl%descr%atoms%astruct%geocode == "F") then
     icoulomb = 1
   else if (wvl%descr%atoms%astruct%geocode == "S") then
     icoulomb = 2
   else
     icoulomb = 0
   end if
!  calculation of the Poisson kernel anticipated to reduce memory peak for small systems
   call psolver_kernel( wvl%den%denspot%dpbox%hgrids, 1, icoulomb, mpi_enreg%me_wvl, wvl%den%denspot%pkernel , &
&   mpi_enreg%comm_wvl, wvl%den%denspot%dpbox%ndims, mpi_enreg%nproc_wvl, dtset%nscforder)
   nullify(wvl%den%denspot%pkernelseq%kernel)
   !call copy_coulomb_operator(wvl%den%denspot%pkernel,wvl%den%denspot%pkernelseq, "gstate")
!  Associate the denspot distribution into mpi_enreg.
   mpi_enreg%nscatterarr  => wvl%den%denspot%dpbox%nscatterarr
   mpi_enreg%ngatherarr   => wvl%den%denspot%dpbox%ngatherarr
   mpi_enreg%ngfft3_ionic =  wvl%den%denspot%dpbox%n3pi
   call wvl_setngfft(mpi_enreg%me_wvl, dtset%mgfft, dtset%nfft, &
&   dtset%ngfft, mpi_enreg%nproc_wvl, wvl%den%denspot%dpbox%ndims(1), &
&   wvl%den%denspot%dpbox%ndims(2),wvl%den%denspot%dpbox%ndims(3),&
&   wvl%den%denspot%dpbox%nscatterarr(mpi_enreg%me_wvl, 1))
#endif
   nfftf     = dtset%nfft
   mgfftf    = dtset%mgfft
   ngfftf(:) = dtset%ngfft(:)
!  Recalculate gprimd
   call matr3inv(rprimd,gprimd)
!  PAW section
   if(psps%usepaw==1) then
!    Reinitialize Pawfgr with new values of ngfft
     call pawfgr_destroy(pawfgr)
     call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)
!    fill wvl objects from paw objects
!    wvl%descr%paw%usepaw=dtset%usepaw
     call paw2wvl(pawtab,wvl%projectors,wvl%descr)
   end if
 else if (dtset%icoulomb /= 0) then
#if defined HAVE_BIGDFT
   if (dtset%ixc < 0) then
     call xc_init(wvl%den%denspot%xc, dtset%ixc, XC_MIXED, dtset%nsppol)
   else
     call xc_init(wvl%den%denspot%xc, dtset%ixc, XC_ABINIT, dtset%nsppol)
   end if
#endif
 end if

!Initialize band structure datatype
 if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
   ABI_ALLOCATE(npwarr_,(dtset%nkpt))
   npwarr_(:)=npwarr(:)
   call xmpi_sum(npwarr_,mpi_enreg%comm_bandfft,ierr)
 else
   npwarr_ => npwarr
 end if

 bstruct = ebands_from_dtset(dtset, npwarr_)

 if (dtset%paral_kgb/=0)  then
   ABI_DEALLOCATE(npwarr_)
 end if
 nullify(npwarr_)

!Initialize PAW atomic occupancies
 if (scf_history%history_size>=0) then
   pawrhoij => scf_history%pawrhoij_last
 else
   ABI_DATATYPE_ALLOCATE(pawrhoij,(my_natom*psps%usepaw))
 end if
 if (psps%usepaw==1.and.has_to_init) then
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,&
&   dtset%lpawu,my_natom,dtset%natom,dtset%nspden,dtset%nspinor,dtset%nsppol,&
&   dtset%ntypat,pawrhoij,dtset%pawspnorb,pawtab,cplex1,dtset%spinat,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr%update(bantot,etot,fermie,&
   residm,rprimd,occ,pawrhoij,xred,args_gs%amu,&
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 ! PW basis set: test if the problem is ill-defined.
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
   npwmin=minval(hdr%npwarr(:))
   if (dtset%mband > npwmin) then
     ! No way we can solve the problem. Abort now!
     write(message,"(2(a,i0),4a)")&
&     "Number of bands nband= ",dtset%mband," > number of planewaves npw= ",npwmin,ch10,&
&     "The number of eigenvectors cannot be greater that the size of the Hamiltonian!",ch10,&
&     "Action: decrease nband or, alternatively, increase ecut"
     if (dtset%ionmov/=23) then
       MSG_ERROR(message)
     else
       MSG_WARNING(message)
     end if

   else if (dtset%mband >= 0.9 * npwmin) then
     ! Warn the user
     write(message,"(a,i0,a,f6.1,4a)")&
&     "Number of bands nband= ",dtset%mband," >= 0.9 * maximum number of planewaves= ",0.9*npwmin,ch10,&
&     "The problem is ill-defined and the GS algorithm will show numerical instabilities!",ch10,&
&     "Assume experienced user. Execution will continue."
     MSG_WARNING(message)
   end if
 end if

!###########################################################
!### 04. Symmetry operations when nsym>1

!Do symmetry stuff only for nsym>1
 if (dtset%usewvl == 0) then
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 else
#if defined HAVE_BIGDFT
   nfftot=product(wvl%den%denspot%dpbox%ndims)
#else
   BIGDFT_NOTENABLED_ERROR()
#endif
 end if
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,dtset%natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon(:,:,:)=0
 phnons(:,:,:)=zero
 indsym(:,:,:)=0
 symrec(:,:,:)=0

 if (dtset%nsym>1) then
   call setsym(indsym,irrzon,dtset%iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!  Make sure dtset%iatfix does not break symmetry
   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)
 else
!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1
 end if
 if (dtset%usewvl == 1) then
   call wvl_descr_atoms_set_sym(wvl%descr, dtset%efield, irrzon, dtset%nsppol, &
&   dtset%nsym, phnons, dtset%symafm, dtset%symrel, dtset%tnons, dtset%tolsym)
#if defined HAVE_BIGDFT
   wvl%den%symObj = wvl%descr%atoms%astruct%sym%symObj
#endif
 end if

!###########################################################
!### 05. Calls inwffil

 ! if paral_kgb == 0, it may happen that some processors are idle (no entry in proc_distrb)
 ! but mkmem == nkpt and this can cause integer overflow in mcg or allocation error.
 ! Here we count the number of states treated by the proc. if cnt == 0, mcg is then set to 0.
 cnt = 0
 do spin=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     do band=1,dtset%nband(ikpt + (spin-1) * dtset%nkpt)
       if (.not. proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, band, band, spin, mpi_enreg%me_kpt)) cnt = cnt + 1
     end do
   end do
 end do

 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcg=dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 if (cnt == 0) then
   mcg = 0
   write(message,"(2(a,i0))")"rank: ",mpi_enreg%me, "does not have wavefunctions to treat. Setting mcg to: ",mcg
   MSG_WARNING(message)
 end if

 if (dtset%usewvl == 0 .and. dtset%mpw > 0 .and. cnt /= 0)then
   if (my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol > floor(real(HUGE(0))/real(dtset%mpw) )) then
     write (message,'(9a)')&
&     "Default integer is not wide enough to store the size of the wavefunction array (mcg).",ch10,&
&     "This usually happens when paral_kgb == 0 and there are not enough procs to distribute kpts and spins",ch10,&
&     "Action: if paral_kgb == 0, use nprocs = nkpt * nsppol to reduce the memory per node.",ch10,&
&     "If this does not solve the problem, use paral_kgb 1 with nprocs > nkpt * nsppol and use npfft/npband/npspinor",ch10,&
&     "to decrease the memory requirements. Consider also OpenMP threads."
     MSG_ERROR_NOSTOP(message,ii)
     write (message,'(5(a,i0), 2a)')&
&     "my_nspinor: ",my_nspinor, ", mpw: ",dtset%mpw, ", mband: ",dtset%mband,&
&     ", mkmem: ",dtset%mkmem, ", nsppol: ",dtset%nsppol,ch10,&
&     'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS...) compiled in int64 mode'
     MSG_ERROR(message)
   end if
 end if

 if (dtset%imgwfstor==1) then
   cg => scf_history%cg(:,:,1)
   eigen => scf_history%eigen(:,1)
 else
   ABI_MALLOC_OR_DIE(cg,(2,mcg), ierr)
   ABI_ALLOCATE(eigen,(dtset%mband*dtset%nkpt*dtset%nsppol))
 end if

 ABI_ALLOCATE(resid,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen(:)=zero ; resid(:)=zero
!mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0

!WVL - Branching, allocating wavefunctions as wavelets.
 if (dtset%usewvl == 1) then
   call wvl_wfs_lr_copy(wvl%wfs, wvl%descr)
!  Create access arrays for wavefunctions and allocate wvl%wfs%psi (other arrays are left unallocated).
   call wvl_wfs_set(dtset%strprecon,dtset%spinmagntarget, dtset%kpt, mpi_enreg%me_wvl,&
&   dtset%natom, sum(dtset%nband), &
&   dtset%nkpt, mpi_enreg%nproc_wvl, dtset%nspinor, dtset%nsppol, dtset%nwfshist, occ, &
&   psps, rprimd, wvl%wfs, dtset%wtk, wvl%descr, dtset%wvl_crmult, dtset%wvl_frmult, &
&   xred)
!  We transfer wavelets information to the hdr structure.
#if defined HAVE_BIGDFT
   call local_potential_dimensions(mpi_enreg%me_wvl,wvl%wfs%ks%lzd,wvl%wfs%ks%orbs,wvl%den%denspot%xc,&
&   wvl%den%denspot%dpbox%ngatherarr(0,1))
   hdr%nwvlarr(1) = wvl%wfs%ks%lzd%Glr%wfd%nvctr_c
   hdr%nwvlarr(2) = 7 * wvl%wfs%ks%lzd%Glr%wfd%nvctr_f
#endif
!  Create access arrays for projectors and allocate them.
!  Compute projectors from each atom.
   call wvl_projectors_set(mpi_enreg%me_wvl, dtset%natom, wvl%projectors, psps, rprimd, &
&   wvl%wfs, wvl%descr, dtset%wvl_frmult, xred)
 end if

 read_wf_or_den=(dtset%iscf<=0.or.dtfil%ireadwf/=0.or.(dtfil%ireadden/=0.and.dtset%positron<=0))
 read_wf_or_den=(read_wf_or_den.and.has_to_init)

!RECURSION -  initialization
 if(has_to_init .and. dtset%userec==1) then
   call InitRec(dtset,mpi_enreg,rec_set,rmet,maxval(psps%indlmn(3,:,:)))
 end if

!LOTF - initialization
#if defined HAVE_LOTF
 if(has_to_init .and. dtset%ionmov==23) then
   call lotfparam_init(dtset%natom,dtset%lotf_version,1,&
&   dtset%lotf_nitex,dtset%lotf_nneigx,&
&   dtset%lotf_classic,1,1)
 end if
#endif

!Initialize wavefunctions.
 if(dtset%imgwfstor==1 .and. initialized==1)then
   cg(:,:)=scf_history%cg(:,:,1)
   eigen(:)=scf_history%eigen(:,1)
 else if(dtset%tfkinfunc /=2) then
!if(dtset%tfkinfunc /=2) then
   wff1%unwff=dtfil%unwff1
   optorth=1   !if (psps%usepaw==1) optorth=0
   if(psps%usepaw==1 .and. dtfil%ireadwf==1)optorth=0
   hdr%rprimd=rprimd_for_kg ! We need the rprimd that was used to generate de G vectors
   call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,&
&   dtset%exchn2n3d,formeig,hdr,dtfil%ireadwf,dtset%istwfk,kg,&
&   dtset%kptns,dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,ngfft,dtset%nkpt,npwarr,&
&   dtset%nsppol,dtset%nsym,occ,optorth,dtset%symafm,&
&   dtset%symrel,dtset%tnons,dtfil%unkg,wff1,wffnow,dtfil%unwff1,&
&   dtfil%fnamewffk,wvl)
   hdr%rprimd=rprimd
 end if

 if (psps%usepaw==1.and.dtfil%ireadwf==1)then
   call pawrhoij_copy(hdr%pawrhoij,pawrhoij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
!  Has to update header again (because pawrhoij has changed)  -  MT 2007-10-22: Why ?
!  call hdr%update(bantot,etot,fermie,residm,rprimd,occ,pawrhoij,xred,args_gs%amu, &
!  &               comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!###########################################################
!### 06. Operations related to restartxf (Old version)

!TODO: Remove ab_xfh
!Initialize xf history (should be put in inwffil)
 ab_xfh%nxfh=0
 if(dtset%restartxf>=1 .and. dtfil%ireadwf==1)then

!  Should exchange the data about history in parallel localrdwf==0
   if(xmpi_paral==1 .and. dtset%localrdwf==0)then
     write(message, '(a,a,a)' )&
&     'It is not yet possible to use non-zero restartxf,',ch10,&
&     'in parallel, when localrdwf=0. Sorry for this ...'
     MSG_BUG(message)
   end if

   ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,0))
   call outxfhist(ab_xfh,dtset%natom,2,wff1,ios)
   ABI_DEALLOCATE(ab_xfh%xfhist)

   if(ios>0)then
     write(message,'(a,a,a)')&
&     'An error occurred reading the input wavefunction file,',ch10,&
&     'with restartxf=1.'
     MSG_ERROR(message)
   else if(ios==0)then
     write(message, '(a,a,i4,a)' )ch10,&
&     ' gstate : reading',ab_xfh%nxfh,' (x,f) history pairs from input wf file.'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if
!  WARNING : should check that restartxf is not negative
!  WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by dtset%ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 ab_xfh%mxfh=(ab_xfh%nxfh-dtset%restartxf+1)+dtset%ntime+5
 ABI_ALLOCATE(ab_xfh%xfhist,(3,dtset%natom+4,2,ab_xfh%mxfh))
 ab_xfh%xfhist(:,:,:,:) = zero
!WARNING : should check that the number of atoms in the wf file and natom are the same

!Initialize the xf history array
 if(ab_xfh%nxfh>=dtset%restartxf .and. ab_xfh%nxfh>0)then
!  Eventually skip some of the previous history
   if(dtset%restartxf>=2)then
     do ixfh=1,dtset%restartxf-1
       call WffReadSkipRec(ios,1,wff1)
     end do
   end if

!  Read and store the relevant history
   ab_xfh%nxfhr=ab_xfh%nxfh-dtset%restartxf+1
   call outxfhist(ab_xfh,dtset%natom,3,wff1,ios)
 end if

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
   call WffClose(wff1,ierr)
 end if

!###########################################################
!### 07. Calls setup2

!Further setup
 ABI_ALLOCATE(start,(3,dtset%natom))
 call setup2(dtset,npwtot,start,wvl%wfs,xred)

!Allocation of previous atomic positions
 if (scf_history%history_size>=0) then
   xred_old => scf_history%xred_last
 else
   ABI_ALLOCATE(xred_old,(3,dtset%natom))
 end if
 if (has_to_init) xred_old=xred

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,3,tsec)

!###########################################################
!### 08. Compute new occupation numbers

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
& (dtset%occopt>=3.and.dtset%occopt<=8) .and. &
& (dtset%iscf>0 .or. dtset%iscf==-3) .and. dtset%positron/=1 ) then

   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
!  Do not take into account the possible STM bias
   call newocc(doccde,eigen,results_gs%energies%entropy,&
&   results_gs%energies%e_fermie,&
&   dtset%spinmagntarget,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,&
&   dtset%occopt,dtset%prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
   if (dtset%dmftcheck>=0.and.dtset%usedmft>=1.and.(sum(args_gs%upawu(:))>=tol8.or.  &
&   sum(args_gs%jpawu(:))>tol8).and.dtset%dmft_entropy==0) results_gs%energies%entropy=zero
   ABI_DEALLOCATE(doccde)

!  Transfer occupations to bigdft object:
   if(dtset%usewvl==1 .and. .not. wvlbigdft) then
     call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,1,wvl%wfs)
!    call wvl_energies_abi2big(results_gs%energies,wvl%wfs,2)
   end if

 else
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
   results_gs%energies%entropy=zero
 end if

!###########################################################
!### 09. Generate an index table of atoms

!Definition of atindx array
!Generate an index table of atoms, in order for them to be used type after type.
 ABI_ALLOCATE(atindx,(dtset%natom))
 ABI_ALLOCATE(atindx1,(dtset%natom))
 ABI_ALLOCATE(nattyp,(psps%ntypat))
 indx=1
 do itypat=1,psps%ntypat
   nattyp(itypat)=0
   do iatom=1,dtset%natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Compute structure factor phases for current atomic pos:
 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init)) then
   ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
   call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 end if

!Here allocation of GPU for vtorho calculations
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,npwarr,2,psps,dtset%use_gpu_cuda)
 end if
#endif

!###########################################################
!### 10. PAW related operations

!Initialize paw_dmft, even if neither dmft not paw are used
!write(std_out,*) "dtset%usedmft",dtset%usedmft
 use_sc_dmft=dtset%usedmft
! if(dtset%paral_kgb>0) use_sc_dmft=0
 call init_sc_dmft(dtset%nbandkss,dtset%dmftbandi,dtset%dmftbandf,dtset%dmft_read_occnd,dtset%mband,&
& dtset%nband,dtset%nkpt,dtset%nspden, &
& dtset%nspinor,dtset%nsppol,occ,dtset%usedmft,paw_dmft,use_sc_dmft,dtset%dmft_solv,mpi_enreg)
 if (paw_dmft%use_dmft==1.and.me==0) then
   call readocc_dmft(paw_dmft,dtfil%filnam_ds(3),dtfil%filnam_ds(4))
 end if
 !Should be done inside init_sc_dmft
 if ( dtset%usedmft /= 0 ) then
   call data4entropyDMFT_init(paw_dmft%forentropyDMFT,&
   dtset%natom,&
   dtset%typat,&
   dtset%lpawu,&
   dtset%dmft_t2g==1, &
   args_gs%upawu,&
   args_gs%jpawu)
 end if
!write(std_out,*) "paw_dmft%use_dmft",paw_dmft%use_dmft

!PAW: 1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
!  1-
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

!  Test if we have to call pawinit
!  Some gen-cond have to be added...
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or.call_pawinit) then
     call timab(553,1,tsec)
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     hyb_range_fock=zero;if (dtset%ixc<0) call libxc_functionals_get_hybridparams(hyb_range=hyb_range_fock)
     call pawinit(dtset%effmass_free,gnt_option,gsqcut_shp,hyb_range_fock,dtset%pawlcutd,dtset%pawlmix,&
&     psps%mpsang,dtset%pawnphi,dtset%nsym,dtset%pawntheta,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%xclevel,dtset%usekden,dtset%usepotzero)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit)
     call timab(553,2,tsec)
#if defined HAVE_BIGDFT
!    In the PAW+WVL case, copy sij:
     if(dtset%usewvl==1) then
       do itypat=1,dtset%ntypat
         wvl%descr%paw%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
       end do
     end if
#endif
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsym_ylm(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,rprimd,symrec,pawang%zarot)

!  2-Initialize and compute data for LDA+U, EXX, or LDA+DMFT
   if(paw_dmft%use_dmft==1) call print_sc_dmft(paw_dmft,dtset%pawprtvol)
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
&     is_dfpt,args_gs%jpawu,dtset%lexexch,dtset%lpawu,dtset%ntypat,pawang,dtset%pawprtvol,&
&     pawrad,pawtab,args_gs%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu,ucrpa=dtset%ucrpa)
 end if

!###########################################################
!### 11. Initialize (eventually) electron-positron data and
!###     electric and magnetic field data

!Initialize (eventually) electron-positron data
 nullify (electronpositron)
 if (dtset%positron/=0) then
   call init_electronpositron(dtfil%ireadwf,dtset,electronpositron,mpi_enreg,nfftf,pawrhoij,pawtab)
 end if

!###########################################################
! Initialisation of cprj
 usecprj=0; mcprj=0;mband_cprj=0
 if (dtset%usepaw==1) then
   if (associated(electronpositron)) then
     if (dtset%positron/=0.and.electronpositron%dimcprj>0) usecprj=1
   end if
   if (dtset%prtnabla>0) usecprj=1
   if (dtset%extrapwf>0) usecprj=1
   if (dtset%pawfatbnd>0)usecprj=1
   if (dtset%prtdos==3)  usecprj=1
   if (dtset%usewvl==1)  usecprj=1
   if (dtset%nstep==0) usecprj=0
   if (dtset%usefock==1)  usecprj=1
 end if
 if (usecprj==0) then
   ABI_DATATYPE_ALLOCATE(cprj,(0,0))
 end if
 if (usecprj==1) then
   mband_cprj=dtset%mband;if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
   mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
!Was allocated above for valgrind sake so should always be true (safety)
   if (allocated(cprj)) then
     call pawcprj_free(cprj)
     ABI_DATATYPE_DEALLOCATE(cprj)
   end if
   ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
   ncpgr=0
   if (dtset%usefock==1) then
     if (dtset%optforces == 1) then
       ncpgr = 3
     end if
!       if (dtset%optstress /= 0) then
!         ncpgr = 6 ; ctocprj_choice = 3
!       end if
   end if
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(dimcprj_srt,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
   call pawcprj_alloc(cprj,ncpgr,dimcprj_srt)
 end if


!###########################################################
!### 12. Operations dependent of iscf value

!Get starting charge density : rhor as well as rhog
!Also initialize the kinetic energy density
 if (scf_history%history_size>=0) then
   rhor => scf_history%rhor_last
   taur => scf_history%taur_last
 else
   ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))
   ABI_ALLOCATE(taur,(nfftf,dtset%nspden*dtset%usekden))
 end if
 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(taug,(2,nfftf*dtset%usekden))

 if (has_to_init) then
   if (dtset%iscf>0 .or. (dtset%iscf==0 .and. dtset%usewvl==1 )) then ! .and. dtset%usepaw==1)) then

     if(dtfil%ireadden/=0.and.dtset%positron<=0)then
       ! Read density
       rdwrpaw=psps%usepaw; if(dtfil%ireadwf/=0) rdwrpaw=0
       if (dtset%usewvl==0) then
         call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, &
         mpi_enreg, rhor, hdr_den, pawrhoij, comm, check_hdr=hdr, allow_interp=.True.)
         results_gs%etotal = hdr_den%etot; call hdr_den%free()
       else
         fform=52 ; accessfil=0
         if (dtset%iomode == IO_MODE_MPI ) accessfil=4
         if (dtset%iomode == IO_MODE_ETSF) accessfil=3
         call ioarr(accessfil,rhor,dtset,results_gs%etotal,fform,dtfil%fildensin,hdr,&
&         mpi_enreg,ngfftf,cplex1,nfftf,pawrhoij,1,rdwrpaw,wvl%den)
       end if

       if (rdwrpaw/=0) then
         call hdr%update(bantot,etot,fermie,residm,&
           rprimd,occ,pawrhoij,xred,args_gs%amu,&
           comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       end if

       ! Read kinetic energy density
       if(dtfil%ireadkden/=0 .and. dtset%usekden==1 )then
         call read_rhor(dtfil%filkdensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, &
         mpi_enreg, taur, hdr_den, pawrhoij, comm, check_hdr=hdr, allow_interp=.True.)
         call hdr_den%free()
       end if

!      Compute up+down rho(G) by fft
       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
       if(dtset%usekden==1)then
         call fourdp(1,taug,taur(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
       end if

     else if(dtfil%ireadwf/=0)then
       izero=0
!      Obtain the charge density from wfs that were read previously
!      Be careful: in PAW, rho does not include the compensation
!      density (to be added in scfcv.F90) !
!      tim_mkrho=1 ; mpi_enreg%paralbd=0
       tim_mkrho=1
       if (psps%usepaw==1) then
         ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
         ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
!        write(std_out,*) "mkrhogstate"
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
         call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
         if(dtset%usekden==1)then
           call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&           mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
           call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,taug,rhowfr,taur)
         end if
         ABI_DEALLOCATE(rhowfg)
         ABI_DEALLOCATE(rhowfr)
       else
         call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&         mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
         if(dtset%usekden==1)then
           call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&           mpi_enreg,npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
         end if

       end if

     else if(dtfil%ireadwf==0.and.dtset%positron/=1)then

!      Crude, but realistic initialisation of the density
!      There is not point to compute it from random wavefunctions.
       if (dtset%usewvl == 0) then
         call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,&
&         mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,nfftf,&
&         ngfftf,dtset%nspden,psps%ntypat,psps,pawtab,ph1df,&
&         psps%qgrid_vl,rhog,rhor,dtset%spinat,ucvol,psps%usepaw,&
&         dtset%ziontypat,dtset%znucl)
!        Update initialized density taking into account jellium slab
         if(dtset%jellslab/=0) then
           option=2
           ABI_ALLOCATE(work,(nfftf))
           call jellium(gmet,gsqcut_eff,mpi_enreg,nfftf,ngfftf,dtset%nspden,&
&           option,dtset%slabwsrad,rhog,rhor,rprimd,work,dtset%slabzbeg,dtset%slabzend)
           ABI_DEALLOCATE(work)
         end if ! of usejell
!        Kinetic energy density initialized to zero (used only in metaGGAs ... )
         if(dtset%usekden==1)then
           taur=zero ; taug=zero
         end if
       else if (dtset%usewvl/=0) then
!        skip for the moment for wvl+paw, not yet coded
         if(dtset%usepaw==0) then
           !Wavelet density corresponds exactly to the wavefunctions,
           !since wavefunctions are taken from diagonalisation of LCAO.
           call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%den)
         else !usepaw
#if defined HAVE_BIGDFT
           call wvl_initro(atindx1,wvl%descr%atoms%astruct%geocode,wvl%descr%h,mpi_enreg%me_wvl,&
&           dtset%natom,nattyp,nfftf,dtset%nspden,psps%ntypat,&
&           wvl%descr%Glr%d%n1,wvl%descr%Glr%d%n1i,&
&           wvl%descr%Glr%d%n2,wvl%descr%Glr%d%n2i,&
&           wvl%descr%Glr%d%n3,&
&           pawrad,pawtab,psps%gth_params%psppar,rhor,rprimd,&
&           dtset%spinat,wvl%den,dtset%xc_denpos,xred,dtset%ziontypat)
           call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%den)
#endif
         end if
       end if

     end if

   else if ((dtset%iscf==-1.or.dtset%iscf==-2.or.dtset%iscf==-3).and.dtset%positron<=0) then

!    Read rho(r) from a disk file
     rdwrpaw=psps%usepaw
!    Note : results_gs%etotal is read here,
!    and might serve in the tddft routine, but it is contrary to the intended use of results_gs ...
!    Warning : should check the use of results_gs%e_fermie
!    Warning : should check the use of results_gs%residm
!    One might make them separate variables.

    ! Read density and get Fermi level from hdr_den
     call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, &
     mpi_enreg, rhor, hdr_den, pawrhoij, comm, check_hdr=hdr)
     results_gs%etotal = hdr_den%etot; results_gs%energies%e_fermie = hdr_den%fermie
     call hdr_den%free()

     if(dtfil%ireadkden/=0 .and. dtset%usekden==1)then
       call read_rhor(dtfil%filkdensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, &
       mpi_enreg, taur, hdr_den, pawrhoij, comm, check_hdr=hdr)
       call hdr_den%free()
     end if

!    Compute up+down rho(G) by fft
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
     if(dtset%usekden==1)then
       call fourdp(1,taug,taur(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
     end if

   end if
 end if ! has_to_init

!###########################################################
!### 13. If needed, initialize SCF history variables

!If needed, initialize atomic density in SCF history
 if (scf_history%history_size>0.and.has_to_init) then
!  If rhor is an atomic density, just store it in history
   if (.not.read_wf_or_den) then
     scf_history%atmrho_last(:)=rhor(:,1)
   else
!    If rhor is not an atomic density, has to compute rho_at(r)
     ABI_ALLOCATE(rhowfg,(2,nfftf))
     ABI_ALLOCATE(rhowfr,(nfftf,1))
     ABI_ALLOCATE(spinat_dum,(3,dtset%natom))
     spinat_dum=zero
     call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,&
&     psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,1,psps%ntypat,psps,pawtab,&
&     ph1df,psps%qgrid_vl,rhowfg,rhowfr,spinat_dum,ucvol,&
&     psps%usepaw,dtset%ziontypat,dtset%znucl)
     scf_history%atmrho_last(:)=rhowfr(:,1)
     ABI_DEALLOCATE(rhowfg)
     ABI_DEALLOCATE(rhowfr)
     ABI_DEALLOCATE(spinat_dum)
   end if
 end if

 if ((.not.read_wf_or_den).or.(scf_history%history_size>0.and.has_to_init))  then
   ABI_DEALLOCATE(ph1df)
 end if

!!Electric field: initialization stage
!!further initialization and updates happen in scfcv.F90
 call init_e_field_vars(dtefield,dtset,gmet,gprimd,kg,&
& mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
& pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)

 !! orbital magnetization initialization
 dtorbmag%orbmag = dtset%orbmag
 if (dtorbmag%orbmag .NE. 0) then
    call initorbmag(dtorbmag,dtset,gmet,gprimd,kg,mpi_enreg,npwarr,occ,&
&                   pawtab,psps,pwind,pwind_alloc,pwnsfac,&
&                   rprimd,symrec,xred)
 end if


 fatvshift=one

 if (dtset%usewvl == 1 .and. wvl_debug) then
#if defined HAVE_BIGDFT
   call wvl_timing(me,'INIT','PR')
#endif
 end if

!Check whether exiting was required by the user. If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call exit_check(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg%comm_cell,openexit)

!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) eigen(:)=zero

 call timab(34,2,tsec)

 conv_retcode = 0

 if (iexit==0) then

!  ###########################################################
!  ### 14. Move atoms and acell according to ionmov value


!  Eventually symmetrize atomic coordinates over space group elements:
!  call symmetrize_xred(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!  If atoms are not being moved and U should not be determined,
!  use scfcv directly; else
!  call move, pawuj_drive or brdmin which in turn calls scfcv.

   call timab(35,3,tsec)

   call scfcv_init(scfcv_args,atindx,atindx1,cg,cprj,cpus,&
&   args_gs%dmatpawu,dtefield,dtfil,dtorbmag,dtpawuj,dtset,ecore,eigen,hdr,&
&   indsym,initialized,irrzon,kg,mcg,mcprj,mpi_enreg,my_natom,nattyp,ndtpawuj,&
&   nfftf,npwarr,occ,pawang,pawfgr,pawrad,pawrhoij,&
&   pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,&
&   resid,results_gs,scf_history,fatvshift,&
&   symrec,taug,taur,wvl,ylm,ylmgr,paw_dmft,wffnew,wffnow)

   call dtfil_init_time(dtfil,0)

   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if (dtset%ionmov==0 .or. dtset%imgmov==6) then

!    Should merge this call with the call for dtset%ionmov==4 and 5

     if (dtset%macro_uj==0) then
       call scfcv_run(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
     else
!      Conduct determination of U
       call pawuj_drive(scfcv_args,dtset,electronpositron,rhog,rhor,rprimd,xred,xred_old)
     end if

!    ========================================
!    New structure for geometry optimization
!    ========================================
   else if (dtset%ionmov>50.or.dtset%ionmov<=27) then

     ! TODO: return conv_retcode
     call mover(scfcv_args,ab_xfh,acell,args_gs%amu,dtfil,&
&     electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old)

!    Compute rprim from rprimd and acell
     do kk=1,3
       do jj=1,3
         rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
       end do
     end do

!    =========================================
!    New structure for geometry optimization
!    =========================================

   else ! Not an allowed option
     write(message, '(a,i0,2a)' )&
&     'Disallowed value for ionmov=',dtset%ionmov,ch10,&
&     'Allowed values are: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,21,22,23,24 and 30'
     MSG_BUG(message)
   end if

   call scfcv_destroy(scfcv_args)

   call timab(35,2,tsec)

!  ###########################################################
!  ### 15. Final operations and output for gstate

 end if !  End of the check of hasty exit

 call timab(36,3,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ----iterations are completed or convergence reached----',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Mark this GS computation as done
 initialized=1
 if (dtset%usewvl == 1 .and. wvl_debug) then
#if defined HAVE_BIGDFT
   call wvl_timing(me,'WFN_OPT','PR')
#endif
 end if

!Will be put here later.
!! ! WVL - maybe compute the tail corrections to energy
!! if (dtset%tl_radius > zero) then
!!    ! Store xcart for each atom
!!    allocate(xcart(3, dtset%natom))
!!    call xred2xcart(dtset%natom, rprimd, xcart, xred)
!!    ! Use the tails to improve energy precision.
!!    call wvl_tail_corrections(dtset, results_gs%energies, results_gs%etotal, &
!!         & mpi_enreg, occ, psps, vtrial, wvl, xcart)
!!    deallocate(xcart)
!! end if

!Update the header, before using it
 call hdr%update(bantot,results_gs%etotal,results_gs%energies%e_fermie,&
   results_gs%residm,rprimd,occ,pawrhoij,xred,args_gs%amu,&
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde=zero

 call ebands_init(bantot,ebands,dtset%nelect,doccde,eigen,hdr%istwfk,hdr%kptns,hdr%nband,&
& hdr%nkpt,hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk,&
& hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, &
& hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

 ebands%fermie = results_gs%energies%e_fermie
 ABI_DEALLOCATE(doccde)
 !write(std_out,*)"efermi after ebands_init",ebands%fermie

 ! Compute and print the gaps.
 call ebands_report_gap(ebands,header="Gap info",unit=std_out,mode_paral="COLL",gaps=results_gs%gaps)

 if(dtset%nqpt==0)filnam=dtfil%fnameabo_wfk
 if(dtset%nqpt==1)filnam=dtfil%fnameabo_wfq

 ! Write wavefunctions file only if convergence was not achieved.
 !write(std_out,*)"conv_retcode", conv_retcode
 write_wfk = .True.
 if (dtset%prtwf==-1 .and. conv_retcode == 0) then
   write_wfk = .False.
   message = "GS calculation converged with prtwf=-1 --> Skipping WFK file output"
   call wrtout(ab_out, message, "COLL")
   MSG_COMMENT(message)
 end if

!To print out the WFs, need the rprimd that was used to generate the G vectors
 hdr%rprimd=rprimd_for_kg

 if (write_wfk) then
   call outresid(dtset,dtset%kptns,dtset%mband,&
&                dtset%nband,dtset%nkpt,&
&                dtset%nsppol,resid)

   call outwf(cg,dtset,psps,eigen,filnam,hdr,kg,dtset%kptns,&
    dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,&
    dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,&
    occ,response,dtfil%unwff2,wvl%wfs,wvl%descr)

   ! Generate WFK with k-mesh from WFK containing list of k-points inside pockets.
   if (dtset%getkerange_path /= ABI_NOFILE) then
     call wfk_klist2mesh(dtfil%fnameabo_wfk, dtset%getkerange_path, dtset, comm)
   end if

   !SPr: add input variable managing the .vtk file OUTPUT (Please don't remove the next commented line)
   !call printmagvtk(mpi_enreg,cplex1,dtset%nspden,nfftf,ngfftf,rhor,rprimd,'DEN')
 end if

 if (dtset%prtwf==2) call outqmc(cg,dtset,eigen,gprimd,hdr,kg,mcg,mpi_enreg,npwarr,occ,psps,results_gs)

!Restore the original rprimd in hdr
 hdr%rprimd=rprimd

 ! Generate WFK in full BZ (needed by LOBSTER)
 if (me == master .and. dtset%prtwf == 1 .and. dtset%prtwf_full == 1 .and. dtset%nqpt == 0) then
   wfkfull_path = strcat(dtfil%filnam_ds(4), "_FULL_WFK")
   if (dtset%iomode == IO_MODE_ETSF) wfkfull_path = nctk_ncify(wfkfull_path)
   call wfk_tofullbz(filnam, dtset, psps, pawtab, wfkfull_path)
   call cryst%free()
 end if

 call clnup1(acell,dtset,eigen,results_gs%energies%e_fermie,&
& dtfil%fnameabo_dos,dtfil%fnameabo_eig,results_gs%fred,&
& mpi_enreg,nfftf,ngfftf,occ,dtset%optforces,&
& resid,rhor,rprimd,results_gs%vxcavg,xred)

 if ( (dtset%iscf>=0 .or. dtset%iscf==-3) .and. dtset%prtstm==0) then
   call prtene(dtset,results_gs%energies,ab_out,psps%usepaw)
 end if

!write final electric field components HONG

 if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt ==7 .or.  &
& dtset%berryopt == 14 .or. dtset%berryopt == 16 .or. dtset%berryopt ==17 ) then ! output final elctric field data    !!HONG
   if (dtset%berryopt == 4) then
     write(message,'(a,a)')   ch10, 'Constant unreduced E calculation  - final values:'
   else if (dtset%berryopt == 6 ) then
     write(message,'(a,a)')   ch10, 'Constant unreduced D calculation  - final values:'
   else if (dtset%berryopt == 14) then
     write(message,'(a,a)')   ch10, 'Constant reduced ebar calculation  - final values:'
   else if (dtset%berryopt == 16 ) then
     write(message,'(a,a)')   ch10, 'Constant reduced d calculation  - final values:'
   else if (dtset%berryopt == 17) then
     write(message,'(a,a)')   ch10, 'Constant reduced ebar and d calculation  - final values:'
   end if

   call wrtout(ab_out,message,'COLL')
   call prtefield(dtset,dtefield,ab_out,rprimd)

   call wrtout(std_out,message,'COLL')
   call prtefield(dtset,dtefield,std_out,rprimd)

!  To check if the final electric field is below the critical field
   do kk = 1, 3
     efield_band(kk) = abs(dtset%red_efieldbar(kk))*dtefield%nkstr(kk)
   end do
!  eg = maxval(eg_dir)
!  eg_ev = eg*Ha_eV
   write(message,'(a,a,a,a,a,a,a,a,f7.2,a,a)')ch10,&
&   ' Please check: COMMENT - ',ch10,&
&   '  As a rough estimate,',ch10,&
&   '  to be below the critical field, the bandgap of your system',ch10,&
&   '  should be larger than ',maxval(efield_band)*Ha_eV,' eV.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message,'(a)')  '--------------------------------------------------------------------------------'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 end if

!Open the formatted derivative database file, and write the preliminary information
!In the // case, only one processor writes the energy and the gradients to the DDB

 if (me==0.and.dtset%nimage==1.and.((dtset%iscf > 0).or.&
& (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

   if (dtset%iscf > 0) then
     nblok = 2          ! 1st blok = energy, 2nd blok = gradients
   else
     nblok = 1
   end if

   dscrpt=' Note : temporary (transfer) database '
   ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'

   call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&   nblok,xred=xred,occ=occ,ngfft=ngfft)

   call ddb_hdr_open_write(ddb_hdr,ddbnm,dtfil%unddb,fullinit=0)

   call ddb_hdr_free(ddb_hdr)

   choice=2
   mpert = dtset%natom + 6 ; msize = 3*mpert

!  create a ddb structure with just one blok
   call ddb%malloc(msize,1,dtset%natom,dtset%ntypat)

   ddb%flg = 0
   ddb%qpt = zero
   ddb%nrm = one
   ddb%val = zero

!  Write total energy to the DDB
   if (dtset%iscf > 0) then
     ddb%typ(1) = 0
     ddb%val(1,1,1) = results_gs%etotal
     ddb%flg(1,1) = 1
     call ddb%write_block(1,choice,dtset%mband,mpert,msize,dtset%nkpt,dtfil%unddb)
   end if

!  Write gradients to the DDB
   ddb%typ = 4
   ddb%flg = 0
   ddb%val = zero
   indx = 0
   if (dtset%iscf > 0) then
     do iatom = 1, dtset%natom
       do idir = 1, 3
         indx = indx + 1
         ddb%flg(indx,1) = 1
         ddb%val(1,indx,1) = results_gs%fred(idir,iatom)
       end do
     end do
   end if

   indx = 3*dtset%natom + 3
   if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
     do idir = 1, 3
       indx = indx + 1
       if (dtset%rfdir(idir) == 1) then
         ddb%flg(indx,1) = 1
         ddb%val(1,indx,1) = results_gs%pel(idir)
       end if
     end do
   end if

   indx = 3*dtset%natom + 6
   if (dtset%iscf > 0) then
     ddb%flg(indx+1:indx+6,1) = 1
     ddb%val(1,indx+1:indx+6,1) = results_gs%strten(1:6)
   end if

   call ddb%write_block(1,choice,dtset%mband,mpert,msize,dtset%nkpt,dtfil%unddb)
   call ddb%free()

!  Close DDB
   close(dtfil%unddb)
 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0 .and. dtset%positron/=1) then
   call clnup2(psps%n1xccc,results_gs%fred,results_gs%grchempottn,results_gs%gresid,&
&   results_gs%grewtn,results_gs%grvdw,results_gs%grxc,dtset%iscf,dtset%natom,&
&   results_gs%ngrvdw,dtset%optforces,dtset%optstress,dtset%prtvol,start,&
&   results_gs%strten,results_gs%synlgr,xred)
 end if

 if(dtset%imgwfstor==1)then
   scf_history%cg(:,:,1)=cg(:,:)
   scf_history%eigen(:,1)=eigen(:)
 endif

!Deallocate arrays
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(resid)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(start)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(taug)
 ABI_DEALLOCATE(ab_xfh%xfhist)
 call pawfgr_destroy(pawfgr)

 if(dtset%imgwfstor==0)then
   ABI_DEALLOCATE(cg)
   ABI_DEALLOCATE(eigen)
 else
   nullify(cg,eigen)
 endif

 if (dtset%usewvl == 0 .or. dtset%nsym <= 1) then
!  In wavelet case, irrzon and phnons are deallocated by wavelet object.
   ABI_DEALLOCATE(irrzon)
   ABI_DEALLOCATE(phnons)
 end if

 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)

 if (scf_history%history_size<0) then
   if (psps%usepaw==1) then
     call pawrhoij_free(pawrhoij)
   end if
   ABI_DEALLOCATE(rhor)
   ABI_DEALLOCATE(taur)
   ABI_DATATYPE_DEALLOCATE(pawrhoij)
   ABI_DEALLOCATE(xred_old)
 else
   nullify(rhor,taur,pawrhoij,xred_old)
 end if

!PAW+DMFT
 call destroy_sc_dmft(paw_dmft)
 ! This call should be done inside destroy_sc_dmft
 if ( dtset%usedmft /= 0 ) then
   call data4entropyDMFT_destroy(paw_dmft%forentropyDMFT)
 end if

!Destroy electronpositron datastructure
 if (dtset%positron/=0) then
   call destroy_electronpositron(electronpositron)
 end if

!Deallocating the basis set.
 if (dtset%usewvl == 1) then
   call wvl_projectors_free(wvl%projectors)
   call wvl_wfs_free(wvl%wfs)
   call wvl_descr_free(wvl%descr)
   call wvl_denspot_free(wvl%den)
   if(dtset%usepaw == 1) then
     call wvl_paw_free(wvl%descr)
   end if
 end if

 ABI_DEALLOCATE(kg)

 if (dtset%icoulomb /= 0) then
   call psolver_kernel((/ 0._dp, 0._dp, 0._dp /), 0, dtset%icoulomb, 0, kernel_dummy, &
&   0, dtset%ngfft, 1, dtset%nscforder)
 end if

 if (associated(pwind)) then
   ABI_DEALLOCATE(pwind)
 end if
 if (associated(pwnsfac)) then
   ABI_DEALLOCATE(pwnsfac)
 end if
 if ((dtset%berryopt<0).or.&
& (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17)) then
   if (xmpi_paral == 1) then
     ABI_DEALLOCATE(mpi_enreg%kptdstrb)
     if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
&     dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
       ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
     end if
   end if
   if (allocated(mpi_enreg%kpt_loc2ibz_sp))  then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   end if
   if (allocated(mpi_enreg%kpt_loc2fbz_sp)) then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2fbz_sp)
   end if
   if (allocated(mpi_enreg%mkmem)) then
     ABI_DEALLOCATE(mpi_enreg%mkmem)
   end if
 end if
 ! deallocate cprj
 if(usecprj==1) then
   ABI_DEALLOCATE(dimcprj_srt)
   call pawcprj_free(cprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cprj)

 ! deallocate efield
 call destroy_efield(dtefield)

 ! deallocate orbmag
 call destroy_orbmag(dtorbmag)

!deallocate Recursion
 if (dtset%userec == 1) then
   call CleanRec(rec_set)
 end if

 call hdr%free()
 call ebands_free(ebands)

 if (me == master .and. dtset%prtxml == 1) then
!  The dataset given in argument has been treated, then we output its variables.
!  call outvarsXML()
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  </dataSet>'
 end if

 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2 .and. dtset%optdriver /= RUNL_GWLS) then
!  Plane-wave case
   call bandfft_kpt_destroy_array(bandfft_kpt,mpi_enreg)
 end if

 if(dtset%wfoptalg == 1 .and. psps%usepaw == 1) then
   call destroy_invovl(dtset%nkpt)
 end if

 if(gemm_nonlop_use_gemm) then
   call destroy_gemm_nonlop(dtset%nkpt)
   gemm_nonlop_use_gemm = .false.
 end if

!Eventually clean cuda runtime
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call dealloc_hamilt_gpu(2,dtset%use_gpu_cuda)
 end if
#endif

 call timab(36,2,tsec)
 call timab(32,2,tsec)

 DBG_EXIT("COLL")
end subroutine gstate
!!***

!!****f* m_gstate/setup2
!!
!! NAME
!! setup2
!!
!! FUNCTION
!! Call within main routine for setup of various arrays.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ecut=kinetic energy cutoff for planewave basis (hartree)
!!   | natom=number of atoms in unit cell
!!   | nkpt=number of k points
!!   | wtk(nkpt)=integration weight associated with each k point
!!   | iscf=parameter controlling scf or non-scf choice
!!  npwtot(nkpt)=number of planewaves in basis and boundary at each k point
!!  xred(3,natom)=starting reduced atomic coordinates
!!
!! OUTPUT
!!  start(3,natom)=copy of starting xred
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine setup2(dtset,npwtot,start,wfs,xred)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 integer,intent(in) :: npwtot(dtset%nkpt)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: start(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ikpt,npw
 real(dp) :: arith,geom,wtknrm
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' setup2 : enter '
!ENDDEBUG

   if (dtset%iscf>=0) then

!  Copy coordinates into array start
     start(:,:)=xred(:,:)

     if (dtset%usewvl == 0) then
!    Get average number of planewaves per k point:
!    both arithmetic and GEOMETRIC averages are desired--
!    need geometric average to use method of Francis and Payne,
!    J. Phys.: Condens. Matter 2, 4395-4404 (1990) [[cite:Francis1990]].
!    Also note: force k point wts to sum to 1 for this averaging.
!    (wtk is not forced to add to 1 in a case with occopt=2)
       arith=zero
       geom=one
       wtknrm=zero
       do ikpt=1,dtset%nkpt
         npw=npwtot(ikpt)
         wtknrm=wtknrm+dtset%wtk(ikpt)
         arith=arith+npw*dtset%wtk(ikpt)
         geom=geom*npw**dtset%wtk(ikpt)
       end do

!    Enforce normalization of weights to 1
       arith=arith/wtknrm
       geom=geom**(1.0_dp/wtknrm)

     end if

!  Ensure portability of output thanks to tol8
     if (dtset%usewvl == 0) then
       write(message, '(a,2f12.3)' ) &
&       '_setup2: Arith. and geom. avg. npw (full set) are',arith+tol8,geom
     else
#if defined HAVE_BIGDFT
       write(message, '(a,2I8)' ) ' setup2: nwvl coarse and fine are', &
&       wfs%ks%lzd%Glr%wfd%nvctr_c, wfs%ks%lzd%Glr%wfd%nvctr_f
#endif
     end if
     call wrtout(ab_out,  message,'COLL')
     call wrtout(std_out, message,'COLL')

   end if

#if !defined HAVE_BIGDFT
   if (.false.) write(std_out,*) wfs%ks
#endif

 end subroutine setup2
!!***

!!****f* m_gstate/clnup1
!! NAME
!! clnup1
!!
!! FUNCTION
!! Perform "cleanup" at end of execution of gstate routine.
!!
!! INPUTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  dosdeltae=DOS delta of Energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree) for all bands
!!                           at each k point
!!  enunit=choice for units of output eigenvalues: 0=>hartree,
!!   1=> eV, 2=> hartree and eV
!!  fermie=fermi energy (Hartree)
!!  fnameabo_dos=filename of output DOS file
!!  fnameabo_eig=filename of output EIG file
!!  fred(3,natom)=d(E)/d(xred) (hartree)
!!  iatfix(3,natom)=0 if not fixed along specified direction,
!!                  1 if fixed
!!  iscf=parameter controlling scf or non-scf choice
!!  kptopt=option for the generation of k points
!!  kptns(3,nkpt)=k points in terms of recip primitive translations
!!  mband=maximum number of bands
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in unit cell
!!  nband(nkpt*nsppol)=number of bands
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nstep=desired number of electron iteration steps
!!  occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!!  occopt=option for occupancies
!!  prtdos= if == 1, will print the density of states
!!  prtfor= if >0, will print the forces
!!  prtstm= input variable prtstm
!!  prtvol=control print volume and debugging
!!  resid(mband*nkpt*nsppol)=squared residuals for each band and k point where
!!                     resid(n,k)=|<C(n,k)|(H-e(n,k))|C(n,k)>|^2
!!  rhor(nfft,nspden)=electron density (electrons/bohr^3)
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  tphysel="physical" electronic temperature with FD occupations
!!  tsmear=smearing energy or temperature (if metal)
!!  vxcavg=average of vxc potential
!!  wtk(nkpt)=real(dp) array of k-point weights
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  (only print and write to disk)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      getnel,metric,prteigrs,prtrhomxmn,prtxf,write_eig,wrtout
!!
!! SOURCE

subroutine clnup1(acell,dtset,eigen,fermie,&
  & fnameabo_dos,fnameabo_eig,fred,&
  & mpi_enreg,nfft,ngfft,occ,prtfor,&
  & resid,rhor,rprimd,vxcavg,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 integer,intent(in) :: prtfor
 real(dp),intent(in) :: fermie
 real(dp),intent(in) :: vxcavg
 character(len=*),intent(in) :: fnameabo_dos,fnameabo_eig
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: acell(3)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: fred(3,dtset%natom)
 real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: comm,iatom,ii,iscf_dum,iwfrc,me,nnonsc,option,unitdos
 real(dp) :: entropy,grmax,grsum,maxocc,nelect,tolwf,ucvol
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 character(len=500) :: message
 character(len=fnlen) filename
!arrays
 real(dp),allocatable :: doccde(:)

! ****************************************************************

 comm=mpi_enreg%comm_cell; me=xmpi_comm_rank(comm)

 if(dtset%prtstm==0)then ! Write reduced coordinates xred
   write(message, '(a,i5,a)' )' reduced coordinates (array xred) for',dtset%natom,' atoms'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,dtset%natom
     write(message, '(1x,3f20.12)' ) xred(:,iatom)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

!Write reduced gradients if iscf > 0 and dtset%nstep>0 and
 if (dtset%iscf>=0.and.dtset%nstep>0.and.dtset%prtstm==0) then

!  Compute absolute maximum and root mean square value of gradients
   grmax=0.0_dp
   grsum=0.0_dp
   do iatom=1,dtset%natom
     do ii=1,3
!      To be activated in v5.5
!      grmax=max(grmax,abs(fred(ii,iatom)))
       grmax=max(grmax,fred(ii,iatom))
       grsum=grsum+fred(ii,iatom)**2
     end do
   end do
   grsum=sqrt(grsum/dble(3*dtset%natom))

   write(message, '(1x,a,1p,e12.4,a,e12.4,a)' )'rms dE/dt=',grsum,'; max dE/dt=',grmax,'; dE/dt below (all hartree)'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,dtset%natom
     write(message, '(i5,1x,3f20.12)' ) iatom,fred(1:3,iatom)
     call wrtout(ab_out,message,'COLL')
   end do

 end if

 if(dtset%prtstm==0)then

!  Compute and write out dimensional cartesian coords and forces:
   call wrtout(ab_out,' ','COLL')

!  (only write forces if iscf > 0 and dtset%nstep>0)
   if (dtset%iscf<0.or.dtset%nstep<=0.or.prtfor==0) then
     iwfrc=0
   else
     iwfrc=1
   end if

   call prtxf(fred,dtset%iatfix,ab_out,iwfrc,dtset%natom,rprimd,xred)

!  Write length scales
   write(message, '(1x,a,3f16.12,a)' )'length scales=',acell,' bohr'
   call wrtout(ab_out,message,'COLL')
   write(message, '(14x,a,3f16.12,a)' )'=',Bohr_Ang*acell(1:3),' angstroms'
   call wrtout(ab_out,message,'COLL')

 end if

 option=1; nnonsc=0; tolwf=0.0_dp

 if(dtset%iscf<0 .and. dtset%iscf/=-3)option=3
 iscf_dum=dtset%iscf
 if(dtset%nstep==0)iscf_dum=-1

 if(dtset%tfkinfunc==0)then
   if (me == master) then
     call prteigrs(eigen,dtset%enunit,fermie,fnameabo_eig,ab_out,&
&     iscf_dum,dtset%kptns,dtset%kptopt,dtset%mband,&
&     dtset%nband,dtset%nkpt,nnonsc,dtset%nsppol,occ,&
&     dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwf,&
&     vxcavg,dtset%wtk)
     call prteigrs(eigen,dtset%enunit,fermie,fnameabo_eig,std_out,&
&     iscf_dum,dtset%kptns,dtset%kptopt,dtset%mband,&
&     dtset%nband,dtset%nkpt,nnonsc,dtset%nsppol,occ,&
&     dtset%occopt,option,dtset%prteig,dtset%prtvol,resid,tolwf,&
&     vxcavg,dtset%wtk)
   end if

#if defined HAVE_NETCDF
   if (dtset%prteig==1 .and. me == master) then
     filename=trim(fnameabo_eig)//'.nc'
     call write_eig(eigen,filename,dtset%kptns,dtset%mband,dtset%nband,dtset%nkpt,dtset%nsppol)
   end if
#endif
 end if

!Compute and print location of maximal and minimal density
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 call prtrhomxmn(std_out,mpi_enreg,nfft,ngfft,dtset%nspden,2,rhor,ucvol=ucvol)
 if( dtset%prtvol>1)then
   call prtrhomxmn(ab_out,mpi_enreg,nfft,ngfft,dtset%nspden,2,rhor,ucvol=ucvol)
 end if

!If needed, print DOS (unitdos is closed in getnel, occ is not changed if option == 2
 if (dtset%prtdos==1 .and. me == master) then
   if (open_file(fnameabo_dos,message, newunit=unitdos, status='unknown', action="write", form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unitdos)
   maxocc=two/(dtset%nspinor*dtset%nsppol)  ! Will not work in the fixed moment case
   option=2
   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
   call getnel(doccde,dtset%dosdeltae,eigen,entropy,fermie,&
&   maxocc,dtset%mband,dtset%nband,nelect,dtset%nkpt,&
&   dtset%nsppol,occ,dtset%occopt,option,dtset%tphysel,&
&   dtset%tsmear,unitdos,dtset%wtk)
   ABI_DEALLOCATE(doccde)
 end if

end subroutine clnup1
!!***

!!****f* m_gstate/prtxf
!! NAME
!! prtxf
!!
!! FUNCTION
!! Compute and print out dimensional cartesian coordinates and forces.
!! Note: for x=cartesian coordinates, t=reduced coordinates (xred),
!! =>
!!  $ x= R t $
!! =>
!!  $ x(1)=rprimd(1,1) t(1)+rprimd(2,1) t(2)+rprimd(3,1) t(3)$
!!  etc. Also $ t = (R^{-1}) x$ .
!!  To convert gradients, $d(E)/dx(n) = [d(E)/dt(m)] [dt(m)/dx(n)]$
!!  and $ dt(m)/dx(n) = (R^{-1})_{mn} = G_{nm}$ because G is the
!!  inverse transpose of R.  Finally then
!!  $d(E)/dx(n) = G_{nm} [d(E)/dt(m)]$.
!!  The vector $d(E)/dt(m)$ for each atom is input in fred
!!  (grad. wrt xred).
!!
!! INPUTS
!!  fred(3,natom)=gradients of Etot (hartree) wrt xred(3,natom)
!!  iatfix(3,natom)=1 for each fixed atom along specified
!!  direction, else 0
!!  iout=unit number for output file
!!  iwfrc=controls force output: 0=> no forces output,
!!                               1=>forces out in eV/A and Ha/bohr,
!!                               2=>forces out in Ha/bohr
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=relative coordinates of atoms (in terms of prim. transl.)
!!
!! OUTPUT
!!  (data written to unit iout)
!!
!! PARENTS
!!      clnup1
!!
!! CHILDREN
!!      matr3inv,wrtout
!!
!! SOURCE

subroutine prtxf(fred,iatfix,iout,iwfrc,natom,rprimd,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iwfrc,natom
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fred(3,natom),rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu,unfixd
 real(dp) :: convt,fmax,frms
 character(len=15) :: format_line21
 character(len=15) :: format_line25
 character(len=15) :: format_line
 character(len=500) :: message
!arrays
 real(dp) :: favg(3),favg_out(3),ff(3),gprimd(3,3),xx(3)

! ****************************************************************

 format_line21='(i5,1x,3f21.14)'
 format_line25='(i5,1x,3f25.14)'

!Write cartesian coordinates in angstroms
 call wrtout(iout,' cartesian coordinates (angstrom) at end:','COLL')
 do iatom=1,natom
   format_line=format_line21
   do mu=1,3
     xx(mu)=(rprimd(mu,1)*xred(1,iatom)+&
&     rprimd(mu,2)*xred(2,iatom)+&
&     rprimd(mu,3)*xred(3,iatom))*Bohr_Ang
     if(xx(mu)>99999 .or. xx(mu)<-9999)format_line=format_line25
   end do
   write(message,format_line) iatom,xx
   call wrtout(iout,message,'COLL')
 end do

!Optionally write cartesian forces in eV/Angstrom (also provide same in hartree/bohr)
 if (iwfrc/=0) then
!  First, provide results in hartree/bohr
   write(message, '(a,a)' ) ch10,' cartesian forces (hartree/bohr) at end:'
   call wrtout(iout,message,'COLL')
   frms=zero
   fmax=zero
   favg(1)=zero
   favg(2)=zero
   favg(3)=zero
!  To get cartesian forces from input gradients with respect to
!  dimensionless coordinates xred, multiply by G and negate
!  (see notes at top of this subroutine)
   call matr3inv(rprimd,gprimd)
!  First compute (spurious) average force favg
   do iatom=1,natom
     do mu=1,3
       ff(mu)=-(gprimd(mu,1)*fred(1,iatom)+&
&       gprimd(mu,2)*fred(2,iatom)+&
&       gprimd(mu,3)*fred(3,iatom))
       favg(mu)=favg(mu)+ff(mu)
     end do
   end do
   favg(1) = favg(1)/dble(natom)
   favg(2) = favg(2)/dble(natom)
   favg(3) = favg(3)/dble(natom)

!  Subtract off average force in what follows
!  (avg is also subtracted off in carfor, called by loopcv,
!  called by grad)
   unfixd=0
   do iatom=1,natom
     format_line=format_line21
     do mu=1,3
       ff(mu)=-(gprimd(mu,1)*fred(1,iatom)+&
&       gprimd(mu,2)*fred(2,iatom)+&
&       gprimd(mu,3)*fred(3,iatom))-favg(mu)
       if(ff(mu)>99999 .or. ff(mu)<-9999)format_line=format_line25
!      For rms and max force, include only unfixed components
       if (iatfix(mu,iatom) /= 1) then
         unfixd=unfixd+1
         frms=frms+ff(mu)**2
         fmax=max(fmax,abs(ff(mu)))
       end if
     end do
     write(message, format_line) iatom,ff
     call wrtout(iout,message,'COLL')
   end do
   if ( unfixd /= 0 ) frms = sqrt(frms/dble(unfixd))

!  The average force is obtained from the cancellation of numbers
!  of typical size unity, so an absolute value lower
!  than tol14 is meaningless for the output file.
   favg_out(:)=favg(:)
   if(abs(favg_out(1))<tol14)favg_out(1)=zero
   if(abs(favg_out(2))<tol14)favg_out(2)=zero
   if(abs(favg_out(3))<tol14)favg_out(3)=zero

   write(message, '(a,1p,2e14.7,1x,3e11.3,a)' )' frms,max,avg=',frms,fmax,favg_out(1:3),' h/b'
   call wrtout(iout,message,'COLL')

   if (iwfrc==1) then

     write(message, '(a,a)' )ch10,' cartesian forces (eV/Angstrom) at end:'
     call wrtout(iout,message,'COLL')
     convt=Ha_eV/Bohr_Ang

!    Note: subtract off average force
     do iatom=1,natom
       format_line=format_line21
       do mu=1,3
         ff(mu)=(-(gprimd(mu,1)*fred(1,iatom)+&
&         gprimd(mu,2)*fred(2,iatom)+&
&         gprimd(mu,3)*fred(3,iatom))-favg(mu))*convt
         if(ff(mu)>99999 .or. ff(mu)<-9999)format_line=format_line25
       end do
       write(message, format_line) iatom,ff
       call wrtout(iout,message,'COLL')
     end do
     write(message, '(a,1p,2e14.7,1x,3e11.3,a)' )' frms,max,avg=',convt*frms,convt*fmax,convt*favg_out(1:3),' e/A'
     call wrtout(iout,message,'COLL')

   end if
 end if

end subroutine prtxf
!!***

!!****f* m_gstate/clnup2
!! NAME
!! clnup2
!!
!! FUNCTION
!! Perform more "cleanup" after completion of iterations.
!! This subroutine prints out more breakdown of force
!! information, shifts of atomic positions, and stresses.
!!
!! INPUTS
!!  fred(3,natom)=d(E_total)/d(xred) derivatives (hartree)
!!  grchempottn(3,natom)=d(E_chempot)/d(xred) derivatives (hartree)
!!  grewtn(3,natom)=d(E_Ewald)/d(xred) derivatives (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D2 dispersion (hartree)
!!  grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!  iscf=parameter controlling scf or non-scf iterations
!!  natom=number of atoms in unit cell
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  prtfor= >0 if forces have to be printed (0 otherwise)
!!  prtstr= >0 if stresses have to be printed (0 otherwise)
!!  prtvol=control print volume and debugging output
!!  start(3,natom)=starting coordinates in terms of real space
!!   primitive translations
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  synlgr(3,natom)=d(E_nlpsp)/d(xred) derivatives (hartree)
!!  xred(3,natom)=final coordinates in terms of primitive translations
!!
!! OUTPUT
!!  (only print)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine clnup2(n1xccc,fred,grchempottn,gresid,grewtn,grvdw,grxc,iscf,natom,ngrvdw,&
&                 prtfor,prtstr,prtvol,start,strten,synlgr,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,n1xccc,natom,ngrvdw,prtfor,prtstr,prtvol
!arrays
 real(dp),intent(in) :: fred(3,natom),grchempottn(3,natom),gresid(3,natom)
 real(dp),intent(in) :: grewtn(3,natom),grvdw(3,ngrvdw)
 real(dp),intent(in) :: grxc(3,natom),start(3,natom),strten(6),synlgr(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
 character(len=*), parameter :: format01020 ="(i5,1x,3f20.12)"
!scalars
 integer :: iatom,mu
 real(dp) :: devsqr,grchempot2
 character(len=500) :: message

! *************************************************************************
!
!DEBUG
!write(std_out,*)' clnup2 : enter '
!ENDDEBUG

!Only print additional info for scf calculations
 if (iscf>=0) then

   if((prtvol>=10).and.(prtfor>0))then

     write(message, '(a,10x,a)' ) ch10,&
&     '===> extra information on forces <==='
     call wrtout(ab_out,message,'COLL')

     write(message, '(a)' ) ' ewald contribution to reduced grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(grewtn(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

     grchempot2=sum(grchempottn(:,:)**2)
     if(grchempot2>tol16)then
       write(message, '(a)' ) ' chemical potential contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grchempottn(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     write(message, '(a)' ) ' nonlocal contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(synlgr(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

     write(message, '(a)' ) ' local psp contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     if (n1xccc/=0) then
       do iatom=1,natom
         write(message,format01020) iatom,fred(:,iatom)-&
&         (grewtn(:,iatom)+grchempottn(:,iatom)+synlgr(:,iatom)+grxc(:,iatom)+gresid(:,iatom))
         call wrtout(ab_out,message,'COLL')
       end do
     else
       do iatom=1,natom
         write(message,format01020) iatom,fred(:,iatom)-&
&         (grewtn(:,iatom)+grchempottn(:,iatom)+synlgr(:,iatom)+gresid(:,iatom))
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     if (n1xccc/=0) then
       write(message, '(a)' ) ' core charge xc contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grxc(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     if (ngrvdw==natom) then
       write(message, '(a)' ) ' Van der Waals DFT-D contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grvdw(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     write(message, '(a)' ) ' residual contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(gresid(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

   end if

!  Compute mean squared deviation from starting coords
   devsqr=0.0_dp
   do iatom=1,natom
     do mu=1,3
       devsqr=devsqr+(xred(mu,iatom)-start(mu,iatom))**2
     end do
   end do

!  When shift is nonnegligible then print values
   if (devsqr>1.d-14) then
     write(message, '(a,1p,e12.4,3x,a)' ) &
&     ' rms coord change=',sqrt(devsqr/dble(3*natom)),&
&     'atom, delta coord (reduced):'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message, '(1x,i5,2x,3f20.12)' ) iatom,&
&       (xred(mu,iatom)-start(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do
   end if

!  Write out stress results
   if (prtstr>0) then
     write(message, '(a,a)' ) ch10,&
&     ' Cartesian components of stress tensor (hartree/bohr^3)'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!    Also output the pressure (minus one third the trace of the stress
!    tensor.
     write(message, '(a,a,es12.4,a)' ) ch10,&
&     '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&     -(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa]'

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(1 1)=',strten(1)*HaBohr3_GPa,&
&     '  sigma(3 2)=',strten(4)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(2 2)=',strten(2)*HaBohr3_GPa,&
&     '  sigma(3 1)=',strten(5)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(3 3)=',strten(3)*HaBohr3_GPa,&
&     '  sigma(2 1)=',strten(6)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

!  Last end if above refers to iscf > 0
 end if

!DEBUG
!write(std_out,*)' clnup2 : exit '
!ENDDEBUG

end subroutine clnup2
!!***

!!****f* m_gstate/pawuj_drive
!! NAME
!! pawuj_drive
!!
!! FUNCTION
!!  Drive for automatic determination of U
!!  Relevant only in PAW+U context
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points treated by this node.
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions.
!!  dtefield <type(efield_type)> = variables related to Berry phase calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      pawuj_det,pawuj_free,pawuj_ini,scfcv_run
!!
!! SOURCE

subroutine pawuj_drive(scfcv_args, dtset,electronpositron,rhog,rhor,rprimd, xred,xred_old)

!Arguments ------------------------------------
!scalars
 type(scfcv_t), intent(inout) :: scfcv_args
 type(dataset_type),intent(inout) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 !type(wffile_type),intent(inout) :: wffnew,wffnow
!arrays
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)

!Local variables -------------------------
!scalars
 integer,target :: ndtpawuj=4
 integer :: iuj,conv_retcode
 real(dp) :: ures
 !character(len=500) :: message
!arrays
 real(dp),allocatable :: cgstart(:,:)
 type(macro_uj_type),allocatable,target :: dtpawuj(:)
! *********************************************************************

 DBG_ENTER("COLL")

 if (dtset%macro_uj==0) then
   MSG_BUG('Macro_uj must be set !')
 end if

 ABI_DATATYPE_ALLOCATE(dtpawuj,(0:ndtpawuj))
 ABI_ALLOCATE(cgstart,(2,scfcv_args%mcg))

!DEBUG
!write(std_out,*)'pawuj_drive: before ini dtpawuj(:)%iuj ', dtpawuj(:)%iuj
!END DEBUG
 call pawuj_ini(dtpawuj,ndtpawuj)

 cgstart=scfcv_args%cg
 do iuj=1,ndtpawuj
!  allocate(dtpawuj(iuj)%rprimd(3,3)) ! this has already been done in pawuj_ini
   dtpawuj(iuj)%macro_uj=dtset%macro_uj
   dtpawuj(iuj)%pawprtvol=dtset%pawprtvol
   dtpawuj(iuj)%diemix=dtset%diemix
   dtpawuj(iuj)%pawujat=dtset%pawujat
   dtpawuj(iuj)%nspden=dtset%nspden
   dtpawuj(iuj)%rprimd=dtset%rprimd_orig(1:3,1:3,1)
 end do

!allocate(dtpawuj(0)%vsh(0,0),dtpawuj(0)%occ(0,0))

 do iuj=1,2
   if (iuj>1) scfcv_args%cg(:,:)=cgstart(:,:)

!  DEBUG
!  write(std_out,*)'drive_pawuj before count dtpawuj(:)%iuj ', dtpawuj(:)%iuj
!  END DEBUG

   dtpawuj(iuj*2-1)%iuj=iuj*2-1

   scfcv_args%ndtpawuj=>ndtpawuj
   scfcv_args%dtpawuj=>dtpawuj

   !call scfcv_new(ab_scfcv_in,ab_scfcv_inout,dtset,electronpositron,&
!&   paw_dmft,rhog,rhor,rprimd,wffnew,wffnow,xred,xred_old,conv_retcode)
   call scfcv_run(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)

   scfcv_args%fatvshift=scfcv_args%fatvshift*(-one)
 end do

!Calculate Hubbard U (or J)
 call pawuj_det(dtpawuj,ndtpawuj,trim(scfcv_args%dtfil%filnam_ds(4))//"_UJDET.nc",ures)
 dtset%upawu(dtset%typat(dtset%pawujat),1)=ures/Ha_eV

!Deallocations
 do iuj=0,ndtpawuj
   call pawuj_free(dtpawuj(iuj))
 end do

 ABI_DATATYPE_DEALLOCATE(dtpawuj)
 ABI_DEALLOCATE(cgstart)

 DBG_EXIT("COLL")

end subroutine pawuj_drive
!!***

!!****f* ABINIT/outxfhist
!! NAME
!! outxfhist
!!
!! FUNCTION
!!  read/write xfhist
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  option =
!!   1: write
!!   2: read only nxfh
!!   3: read xfhist
!!  response =
!!   0: GS wavefunctions
!!   1: RF wavefunctions
!!  natom = number of atoms in unit cell
!!  mxfh = last dimension of the xfhist array
!!
!! OUTPUT
!!  ios = error code returned by read operations
!!
!! SIDE EFFECTS
!!  nxfh = actual number of (x,f) history pairs, see xfhist array
!!  wff2 = structured info for wavefunctions
!!  xfhist(3,natom+4,2,ab_xfh%mxfh) = (x,f) history array, also including
!!   rprim and stress
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      xderiveread,xderiverrecend,xderiverrecinit,xderivewrecend
!!      xderivewrecinit,xderivewrite
!!
!! SOURCE

subroutine outxfhist(ab_xfh,natom,option,wff2,ios)

 use defs_basis
 use m_abicore
 use m_abimover
 use m_xmpi
 use m_wffile
 use m_errors
#if defined HAVE_NETCDF
 use netcdf
#endif

!Arguments ------------------------------------
 integer          ,intent(in)    :: natom,option
 integer          ,intent(out)   :: ios
 type(wffile_type),intent(inout)    :: wff2
 type(ab_xfh_type),intent(inout) :: ab_xfh

!Local variables-------------------------------
 integer :: ierr,ixfh,ncid_hdr,spaceComm,xfdim2
 real(dp),allocatable :: xfhist_tmp(:)
 character(len=500) :: message
!no_abirules
#if defined HAVE_NETCDF
 integer :: ncerr
 integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
 integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp
#endif


! *************************************************************************

!DEBUG
!write(std_out,*)'outxfhist  : enter, option = ', option
!ENDDEBUG
 ncid_hdr = wff2%unwff
 xfdim2 = natom+4

 ios = 0

!### (Option=1) Write out content of all iterations
!#####################################################################
 if ( option == 1 ) then

!  Write the (x,f) history
   if (wff2%iomode == IO_MODE_FORTRAN) then
     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     write(unit=wff2%unwff)ab_xfh%nxfh
     do ixfh=1,ab_xfh%nxfh
       write(unit=wff2%unwff)ab_xfh%xfhist(:,:,:,ixfh)
     end do

!    insert mpi broadcast here

   else if(wff2%iomode==IO_MODE_MPI)then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     spaceComm=xmpi_comm_self
     call xderiveWRecInit(wff2,ierr)
     call xderiveWrite(wff2,ab_xfh%nxfh,ierr)
     call xderiveWRecEnd(wff2,ierr)
     do ixfh=1,ab_xfh%nxfh
       xfhist_tmp(:)=reshape(ab_xfh%xfhist(:,:,:,ixfh),(/3*(natom+4)*2/))
       call xderiveWRecInit(wff2,ierr)
       call xderiveWrite(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveWRecEnd(wff2,ierr)
     end do
     ABI_DEALLOCATE(xfhist_tmp)

#if defined HAVE_NETCDF
   else if (wff2%iomode == IO_MODE_NETCDF) then
!    check if nxfh and xfhist are defined
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

     if (ncerr /= NF90_NOERR) then
!      need to define everything
       ncerr = nf90_redef (ncid=ncid_hdr)
       NCF_CHECK_MSG(ncerr," outxfhist : going to define mode ")

       ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define dim2inout")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=ab_xfh%mxfh,dimid=mxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define mxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=ab_xfh%nxfh,dimid=nxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define nxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define xfdim2")

       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire dimr3")

!      ab_xfh%xfhist(3,natom+4,2,ab_xfh%mxfh)
       ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
&       dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
       NCF_CHECK_MSG(ncerr," outxfhist : define xfhist")

!      End define mode and go to data mode
       ncerr = nf90_enddef(ncid=ncid_hdr)
       NCF_CHECK_MSG(ncerr," outxfhist : enddef call ")
     else
!      check that the dimensions are correct
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&       len=nxfh_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire xfdim2")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
&       len=xfdim2_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get xfdim2")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire mxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
&       len=mxfh_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get mxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire dim2inout")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
&       len=dim2inout_tmp)
       NCF_CHECK_MSG(ncerr,"  outxfhist : get dim2inout")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
       NCF_CHECK_MSG(ncerr," outxfhist : inquire xfhist")

       if (mxfh_tmp /= ab_xfh%mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
         write (message,"(A)") 'outxfhist : ERROR xfhist has bad dimensions in NetCDF file. Can not re-write it.'
         MSG_ERROR(message)
       end if

     end if

!    Now fill the data
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist)
     NCF_CHECK_MSG(ncerr," outxfhist : fill xfhist")

!    end NETCDF definition ifdef
#endif
   end if  ! end iomode if

!  ### (Option=2) Read in number of iterations
!  #####################################################################
 else if ( option == 2 ) then

   if (wff2%iomode == IO_MODE_FORTRAN) then
     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     read(unit=wff2%unwff,iostat=ios)ab_xfh%nxfh

   else if (wff2%iomode == IO_MODE_MPI) then
     call xderiveRRecInit(wff2,ierr)
     call xderiveRead(wff2,ab_xfh%nxfh,ierr)
     call xderiveRRecEnd(wff2,ierr)

#if defined HAVE_NETCDF
   else if (wff2%iomode == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfh)
     NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")
#endif
   end if

!  ### (Option=3) Read in iteration content
!  #####################################################################
 else if ( option == 3 ) then
   if (wff2%iomode == IO_MODE_FORTRAN) then
     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do
   else if (wff2%iomode == IO_MODE_FORTRAN_MASTER) then
!    FIXME: should copy the xfhist to other processors, and check that we are on the master to read in this case
!    if node is master
     write(message, "(A,A,A,A)") ch10, " outxfhist: ERROR -", ch10, &
&     'iomode == -1 (localrdwf ) has not been coded yet for xfhist rereading.'
     MSG_ERROR(message)

     do ixfh=1,ab_xfh%nxfhr
       read(unit=wff2%unwff,iostat=ios)ab_xfh%xfhist(:,:,:,ixfh)
     end do

   else if (wff2%iomode == IO_MODE_MPI) then
     ABI_ALLOCATE(xfhist_tmp,(3*(natom+4)*2))
     spaceComm=xmpi_comm_self
     do ixfh=1,ab_xfh%nxfhr
       call xderiveRRecInit(wff2,ierr)
       call xderiveRead(wff2,xfhist_tmp,3*(natom+4)*2,spaceComm,ierr)
       call xderiveRRecEnd(wff2,ierr)
       xfhist_tmp(:)=xfhist_tmp(:)
     end do
     ABI_DEALLOCATE(xfhist_tmp)
   end if

!  FIXME: should this be inside the if not mpi as above for options 1 and 2?
!  it is placed here because the netcdf read is a single operation
#if defined HAVE_NETCDF
   if (wff2%iomode == IO_MODE_NETCDF) then
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
     NCF_CHECK_MSG(ncerr," outxfhist : inquire nxfh")
     ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
&     len=ab_xfh%nxfhr)
     NCF_CHECK_MSG(ncerr,"  outxfhist : get nxfh")

     ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=xfhist_id,name="xfhist")
     NCF_CHECK_MSG(ncerr," outxfhist : inquire xfhist")
     ncerr = nf90_get_var(ncid=ncid_hdr,varid=xfhist_id,values=ab_xfh%xfhist,&
&     start=(/1,1,1,1/),count=(/3,natom+4,2,ab_xfh%nxfhr/))
     NCF_CHECK_MSG(ncerr," outxfhist : read xfhist")
   end if
#endif

 else
!  write(std_out,*)' outxfhist : option ', option , ' not available '
   write(message, "(A,A,A,A,I3,A)") ch10, "outxfhist: ERROR -", ch10, &
&   "option ", option, " not available."
   MSG_ERROR(message)
 end if

end subroutine outxfhist
!!***

end module m_gstate
!!***
