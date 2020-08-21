!!****m* ABINIT/m_dtset
!! NAME
!!  m_dtset
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1992-2020 ABINIT group (XG, MG, FJ, DCA, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dtset

 use defs_basis
 use m_abicore
 use m_copy
 use m_errors
 use m_xmpi

 use m_fstrings,     only : inupper
 use m_symtk,        only : mati3inv, littlegroup_q, symatm
 use m_symkpt,       only : symkpt
 use m_geometry,     only : mkrdim, metric, littlegroup_pert, irreducible_set_pert
 use m_parser,       only : intagm, chkvars_in_string

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_dtset/dataset_type
!! NAME
!! dataset_type
!!
!! FUNCTION
!! The dataset_type structured datatype gather all the input variables,
!! except those that are labelled NOT INTERNAL.
!! For one dataset, it is initialized in driver.f, and will not change
!! at all during the treatment of the dataset.
!! The "evolving" input variables are also stored, with their
!! name appended with _orig, to make clear that this is the original
!! value, decided by the user, and not a possibly modified, intermediate value.
!! The following input variables are NOT INTERNAL, that is, they
!! are input variables used to determine other input variables,
!! after suitable processing, and do not appear anymore afterwards
!! (so, they do not appear as components of a dataset_type variable) :
!! cpuh,cpum(but cpus is present),fband,kptbounds,ndivk,ndism,nobj,
!! objaat,objbat,objaax,objbax,objan,objbn,objarf,objbrf,objaro,objbro
!! objatr,objbtr,vaclst,vacuum
!!
!! WARNING: if you modify this datatype, please check whether there might be
!! creation/destruction/copy routines, declared in another part of ABINIT,
!! that might need to take into account your modification.
!!
!! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.
!!
!! Since all these input variables are described in the abinit_help.html and
!! associated html files they are not described in length here ...
!!
!! SOURCE

type, public :: dataset_type

! Integer
 integer :: iomode
 integer :: accuracy
 integer :: adpimd
 integer :: autoparal
 integer :: auxc_ixc
 integer :: awtr = 1
 integer :: bandpp
 integer :: bdeigrf
 integer :: berryopt
 integer :: berrysav
 integer :: berrystep
 integer :: brav = 1
 integer :: brvltt
 integer :: bs_nstates
 integer :: bs_hayd_term = 0
 integer :: builtintest
 integer :: cd_full_grid
 integer :: cd_frqim_method
 integer :: cd_customnimfrqs
 integer :: chkdilatmx
 integer :: chkexit
 integer :: chkprim
 integer :: chksymbreak
 integer :: cineb_start
 integer :: delayperm
 integer :: densfor_pred
 integer :: diismemory
 integer :: dmatpuopt
 integer :: dmatudiag
 integer :: dmft_dc
 integer :: dmft_entropy
 integer :: dmft_iter
 integer :: dmft_kspectralfunc
 integer :: dmft_nlambda
 integer :: dmft_nwli
 integer :: dmft_nwlo
 integer :: dmft_occnd_imag
 integer :: dmft_rslf
 integer :: dmft_read_occnd
 integer :: dmft_solv
 integer :: dmft_t2g
!integer :: dmft_x2my2d
 integer :: dmftbandi
 integer :: dmftbandf
 integer :: dmftcheck
 integer :: dmftctqmc_basis
 integer :: dmftctqmc_check
 integer :: dmftctqmc_correl
 integer :: dmftctqmc_gmove
 integer :: dmftctqmc_grnns
 integer :: dmftctqmc_meas
 integer :: dmftctqmc_mov
 integer :: dmftctqmc_mrka
 integer :: dmftctqmc_order
 integer :: dmftctqmc_triqs_nleg
 integer :: dmftqmc_l
 integer :: dmftqmc_seed
 integer :: dmftqmc_therm
 integer :: d3e_pert1_elfd
 integer :: d3e_pert1_phon
 integer :: d3e_pert2_elfd
 integer :: d3e_pert2_phon
 integer :: d3e_pert2_strs
 integer :: d3e_pert3_elfd
 integer :: d3e_pert3_phon
 integer :: efmas
 integer :: efmas_calc_dirs
 integer :: efmas_deg
 integer :: efmas_dim
 integer :: efmas_n_dirs
 integer :: efmas_ntheta
 integer :: enunit
 integer :: eph_restart = 0
 integer :: eph_task = 1
 integer :: exchn2n3d
 integer :: extrapwf
 integer :: fftgw = 21
 integer :: fockoptmix
 integer :: fock_icutcoul
 integer :: frzfermi
 integer :: ga_algor
 integer :: ga_fitness
 integer :: ga_n_rules
 integer :: getcell = 0
 integer :: getddb = 0
 integer :: getdvdb = 0
 integer :: getddk = 0
 integer :: getdelfd = 0
 integer :: getdkdk = 0
 integer :: getdkde = 0
 integer :: getden = 0
 integer :: getefmas = 0
 integer :: getgam_eig2nkq = 0
 integer :: getocc = 0
 integer :: getpawden = 0
 integer :: getqps = 0
 integer :: getscr = 0
 integer :: getsuscep = 0
 integer :: getvel = 0
 integer :: getwfk = 0
 integer :: getwfkfine = 0
 integer :: getwfq = 0
 integer :: getxcart = 0
 integer :: getxred = 0
 integer :: get1den = 0
 integer :: get1wf = 0
 integer :: getbseig = 0
 integer :: getbsreso = 0
 integer :: getbscoup = 0
 integer :: gethaydock = 0
 integer :: goprecon
 integer :: gwaclowrank = 0
 integer :: gwcalctyp = 0
 integer :: gwcomp = 0
 integer :: gwgamma = 0
 integer :: gwrpacorr = 0
 integer :: gw_customnfreqsp
 integer :: gw_icutcoul
 integer :: gw_invalid_freq
 integer :: gw_qprange
 integer :: gw_nqlwl
 ! TODO: REMOVE?
 integer :: gw_nstep = 30
 integer :: gw_sigxcore = 0

 ! GWLS
 integer :: gwls_stern_kmax             ! number of Lanczos steps taken by the gw_sternheimer routine
 integer :: gwls_npt_gauss_quad         ! number of points used in Gaussian quadrature in gw_sternheimer routine
 integer :: gwls_diel_model             ! switch to determine which dielectic model should be used in integration
 integer :: gwls_print_debug            ! switch to determine what to print out for debugging
 integer :: gwls_nseeds                 ! number of seeds in the Lanczos description of the dielectric matrix
 integer :: gwls_n_proj_freq            ! Number of projection frequencies to be used for the construction of the sternheimer basis
 integer :: gwls_kmax_complement        ! number of Lanczos steps taken in the complement space
 integer :: gwls_kmax_poles             ! number of Lanczos steps taken to compute Poles contribution
 integer :: gwls_kmax_analytic          ! number of Lanczos steps taken to compute the analytic contribution
 integer :: gwls_kmax_numeric           ! number of Lanczos steps taken to compute the numeric contribution
 integer :: gwls_band_index             ! band index of the state to be corrected
 integer :: gwls_exchange               ! Flag to determine if Exchange energy will be computed
 integer :: gwls_correlation            ! Flag to determine if Correlation energy will be computed
 integer :: gwls_first_seed             ! index of the first seed used in the Lanczos algorithm;
                                        ! seeds will go from first_seed to first_seed+nseeds
 !integer :: gwls_n_ext_freq            ! The number of frequencies to be read in gwls_ext_freq
 integer :: gwls_recycle                ! Recycle the sternheimer solutions computed to obtain the static dielectric matrix
                                        ! and add them to the other solutions requested.
                                        ! 0 : don't recycle. 1 : store in RAM. 2 : Store on disk.
 integer :: gw_frqim_inzgrid = 0
 integer :: gw_frqre_inzgrid = 0
 integer :: gw_frqre_tangrid = 0
 integer :: gw_sctype
 integer :: gwmem = 11
 integer :: gwpara = 2
 integer :: hmcsst
 integer :: hmctt
 integer :: iboxcut
 integer :: icoulomb
 integer :: icutcoul
 integer :: ieig2rf
 integer :: imgmov
 integer :: imgwfstor
 integer :: inclvkb = 2
 integer :: intxc
 integer :: ionmov
 integer :: iprcel
 integer :: iprcfc
 integer :: irandom
 integer :: irdddb = 0
 integer :: irddvdb = 0
 integer :: irdddk = 0
 integer :: irdden = 0
 integer :: irdefmas = 0
 integer :: irdhaydock = 0
 integer :: irdpawden = 0
 integer :: irdqps = 0
 integer :: irdscr = 0
 integer :: irdsuscep = 0
 integer :: irdvdw = 0
 integer :: irdwfk = 0
 integer :: irdwfkfine = 0
 integer :: irdwfq = 0
 integer :: ird1den = 0
 integer :: ird1wf = 0
 integer :: irdbseig = 0
 integer :: irdbsreso = 0
 integer :: irdbscoup = 0
 integer :: iscf
 integer :: isecur
 integer :: istatimg
 integer :: istatr
 integer :: istatshft
 integer :: ixc
 integer :: ixc_sigma
 integer :: ixcpositron
 integer :: ixcrot
 integer :: jdtset !  jdtset contains the current dataset number
 integer :: jellslab
 integer :: kptopt
 integer :: kssform = 1
 integer :: localrdwf = 1
 integer :: lotf_classic
 integer :: lotf_nitex
 integer :: lotf_nneigx
 integer :: lotf_version
 integer :: lw_flexo
 integer :: lw_qdrpl
 integer :: magconon
 integer :: maxnsym
 integer :: max_ncpus = 0
 integer :: mband
 integer :: mep_solver
 integer :: mem_test = 1
 integer :: mffmem
 integer :: mgfft
 integer :: mgfftdg
 integer :: mkmem
 integer :: mkqmem
 integer :: mk1mem
 integer :: nnos
 integer :: mpw
 integer :: mqgrid
 integer :: mqgriddg
 integer :: natom
 integer :: natpawu
 integer :: natrd
 integer :: natsph
 integer :: natsph_extra
 integer :: natvshift
 integer :: nbandhf
 integer :: nbandkss = 0
 integer :: nbdblock
 integer :: nbdbuf
 integer :: nberry
 integer :: nc_xccc_gspace = 0
 integer :: nconeq
 integer :: nctime
 integer :: ndtset
 integer :: ndynimage
 integer :: neb_algo
 integer :: nfft
 integer :: nfftdg
 integer :: nfreqim = -1
 integer :: nfreqre = -1
 integer :: nfreqsp = 0
 integer :: nimage
 integer :: nkpt
 integer :: nkptgw
 integer :: nkpthf
 integer :: nline
 integer :: nnsclo
 integer :: nnsclohf
 integer :: nomegasf = 100
 integer :: nomegasi = 12
 integer :: nomegasrd = 9
 integer :: nonlinear_info
 integer :: npband
 integer :: npfft
 integer :: nphf
 integer :: npimage
 integer :: npkpt
 integer :: nppert
 integer :: npspinor
 integer :: npsp
 integer :: npspalch
 integer :: npulayit
 integer :: npvel
 integer :: npweps = 0
 integer :: npwkss = 0
 integer :: npwsigx = 0
 integer :: npwwfn = 0
 integer :: np_slk
 integer :: nqpt
 integer :: nqptdm
 integer :: nscforder
 integer :: nshiftk
 integer :: nshiftk_orig  ! original number of shifts given in input (changed in inkpts, the actual value is nshiftk)
 integer :: nspden
 integer :: nspinor
 integer :: nsppol
 integer :: nstep
 integer :: nsym
 integer :: ntime
 integer :: ntimimage
 integer :: ntypalch
 integer :: ntypat
 integer :: ntyppure
 integer :: nwfshist
 integer :: nzchempot
 integer :: occopt
 integer :: optcell
 integer :: optdriver
 integer :: optforces
 integer :: optnlxccc
 integer :: optstress
 integer :: orbmag
 integer :: ortalg
 integer :: paral_atom
 integer :: paral_kgb
 integer :: paral_rf
 integer :: pawcpxocc
 integer :: pawcross
 integer :: pawfatbnd
 integer :: pawlcutd
 integer :: pawlmix
 integer :: pawmixdg
 integer :: pawnhatxc
 integer :: pawnphi
 integer :: pawntheta
 integer :: pawnzlm
 integer :: pawoptmix
 integer :: pawoptosc
 integer :: pawprtdos
 integer :: pawprtvol
 integer :: pawprtwf
 integer :: pawprt_k
 integer :: pawprt_b
 integer :: pawspnorb
 integer :: pawstgylm
 integer :: pawsushat
 integer :: pawusecp
 integer :: macro_uj
 integer :: pawujat
 integer :: pawxcdev
 integer :: pimd_constraint
 integer :: pitransform
 integer :: plowan_bandi
 integer :: plowan_bandf
 integer :: plowan_compute
 integer :: plowan_natom
 integer :: plowan_nt
 integer :: plowan_realspace
 integer :: posdoppler
 integer :: positron
 integer :: posnstep
 integer :: ppmodel = 1
 integer :: prepalw
 integer :: prepanl
 integer :: prepgkk = 0
 integer :: prtbbb = 0
 integer :: prtbltztrp = 0
 integer :: prtcif = 0
 integer :: prtden
 integer :: prtdensph = 1
 integer :: prtdipole = 0
 integer :: prtdos = 0
 integer :: prtdosm = 0
 integer :: prtebands = 1
 integer :: prtefg = 0
 integer :: prtefmas = 1
 integer :: prteig
 integer :: prtelf = 0
 integer :: prtfc = 0
 integer :: prtfull1wf = 0
 integer :: prtfsurf = 0
 integer :: prtgsr = 1
 integer :: prtgden = 0
 integer :: prtgeo = 0
 integer :: prtgkk = 0
 integer :: prtkden = 0
 integer :: prtkpt
 integer :: prtlden = 0
 integer :: prtnabla = 0
 integer :: prtnest = 0
 integer :: prtpmp
 integer :: prtposcar = 0
 integer :: prtprocar = 0
 integer :: prtphdos = 1
 integer :: prtphbands = 1
 integer :: prtphsurf = 0
 integer :: prtpot = 0
 integer :: prtpsps = 0
 integer :: prtspcur = 0
 integer :: prtstm = 0
 integer :: prtsuscep = 0
 integer :: prtvclmb = 0
 integer :: prtvdw = 0
 integer :: prtvha = 0
 integer :: prtvhxc = 0
 integer :: prtkbff = 0
 integer :: prtvol = 0
 integer :: prtvolimg = 0
 integer :: prtvpsp = 0
 integer :: prtvxc = 0
 integer :: prtwant = 0
 integer :: prtwf
 integer :: prtwf_full = 0
 integer :: prtxml = 0
 integer :: prt1dm = 0
 integer :: ptgroupma
 integer :: qptopt
 integer :: random_atpos
 integer :: recgratio
 integer :: recnpath
 integer :: recnrec
 integer :: recptrott
 integer :: rectesteg
 integer :: restartxf
 integer :: rfasr
 integer :: rfddk
 integer :: rfelfd
 integer :: rfmagn
 integer :: rfmeth
 integer :: rfphon
 integer :: rfstrs
 integer :: rfuser
 integer :: rf2_dkdk
 integer :: rf2_dkde
 integer :: signperm
 integer :: slk_rankpp
 integer :: smdelta
 integer :: spgaxor
 integer :: spgorig
 integer :: spgroup
 integer :: spmeth = 0
 integer :: string_algo
 integer :: symmorphi = 1
 integer :: symchi = 1
 integer :: symsigma = 1
 integer :: td_mexcit
 integer :: tfkinfunc
 integer :: tim1rev
 integer :: timopt
 integer :: tl_nprccg
 integer :: ucrpa
 integer :: use_gpu_cuda
 integer :: usedmatpu
 integer :: usedmft
 integer :: useexexch
 integer :: usefock
 integer :: usekden
 integer :: use_gemm_nonlop
 integer :: use_nonscf_gkk
 integer :: usepaw
 integer :: usepawu
 integer :: usepead
 integer :: usepotzero
 integer :: userec
 integer :: useria = 0
 integer :: userib = 0
 integer :: useric = 0
 integer :: userid = 0
 integer :: userie = 0
 integer :: usewvl
 integer :: usexcnhat_orig
 integer :: useylm
 integer :: use_yaml = 0
 integer :: use_slk
 integer :: vacnum
 integer :: vdw_nfrag
 integer :: vdw_df_ndpts
 integer :: vdw_df_ngpts
 integer :: vdw_df_nqpts
 integer :: vdw_df_nrpts
 integer :: vdw_df_nsmooth
 integer :: vdw_df_tweaks
 integer :: vdw_xc
 integer :: wfoptalg
 integer :: wfk_task
 integer :: wvl_bigdft_comp
 integer :: wvl_nprccg
 integer :: w90iniprj
 integer :: w90prtunk
 integer :: xclevel

!Integer arrays
 integer :: bdberry(4)
 integer :: bravais(11)
 integer :: cd_subset_freq(2)
 integer :: d3e_pert1_atpol(2)
 integer :: d3e_pert1_dir(3)
 integer :: d3e_pert2_atpol(2)
 integer :: d3e_pert2_dir(3)
 integer :: d3e_pert3_atpol(2)
 integer :: d3e_pert3_dir(3)
 integer :: fockdownsampling(3)
 integer :: jfielddir(3)
 integer :: kptrlatt(3,3)
 integer :: kptrlatt_orig(3,3)=0
 integer :: qptrlatt(3,3)
 integer :: ga_rules(30)
 integer :: gpu_devices(5)
 integer :: ngfft(18)
 integer :: ngfftdg(18)
 integer :: nloalg(3)
 integer :: ngkpt(3)   ! Number of division for MP sampling.
 integer :: qprtrb(3)
 integer :: rfatpol(2)
 integer :: rfdir(3)
 integer :: rf2_pert1_dir(3)
 integer :: rf2_pert2_dir(3)
 integer :: supercell_latt(3)
 integer :: ucrpa_bands(2)
 integer :: vdw_supercell(3)
 integer :: vdw_typfrag(100)
 integer :: wvl_ngauss(2)

!Integer allocatables
 integer, allocatable ::  algalch(:)         ! algalch(ntypalch)
 integer, allocatable ::  bdgw(:,:,:)        ! bdgw(2,nkptgw,nsppol)
 integer, allocatable ::  constraint_kind(:) ! constraint_kind(ntypat)
 integer, allocatable ::  dynimage(:)        ! dynimage(nimage or mxnimage)
 integer, allocatable ::  efmas_bands(:,:)   ! efmas_bands(2,nkptgw)
 integer, allocatable ::  iatfix(:,:)        ! iatfix(3,natom)
 integer, allocatable ::  iatsph(:)          ! iatsph(natsph)
 integer, allocatable ::  istwfk(:)          ! istwfk(nkpt)
 integer, allocatable ::  kberry(:,:)        ! kberry(3,nberry)
 integer, allocatable ::  lexexch(:)         ! lexexch(ntypat)
 integer, allocatable ::  ldaminushalf(:)    ! ldaminushalf(ntypat)
 integer, allocatable ::  lpawu(:)           ! lpawu(ntypat)
 integer, allocatable ::  nband(:)           ! nband(nkpt*nsppol)
 integer, allocatable ::  plowan_iatom(:)    ! plowan_iatom(plowan_natom)
 integer, allocatable ::  plowan_it(:)       ! plowan_it(plowan_nt*3)
 integer, allocatable ::  plowan_lcalc(:)    ! plowan_lcalc(\sum_iatom plowan_nbl)
 integer, allocatable ::  plowan_nbl(:)      ! plowan_nbl(plowan_natom)
 integer, allocatable ::  plowan_projcalc(:) ! plowan_projcalc(\sum_iatom plowan_nbl)
 integer, allocatable ::  prtatlist(:)       ! prtatlist(natom)
 integer, allocatable ::  so_psp(:)          ! so_psp(npsp)
 integer, allocatable ::  symafm(:)          ! symafm(nsym)
 integer, allocatable ::  symrel(:,:,:)      ! symrel(3,3,nsym)
 integer, allocatable ::  typat(:)           ! typat(natom)

!Real
 real(dp) :: adpimd_gamma
 real(dp) :: auxc_scal
 real(dp) :: bmass
 real(dp) :: boxcutmin
 real(dp) :: bxctmindg
 real(dp) :: cd_halfway_freq
 real(dp) :: cd_max_freq
 real(dp) :: charge
 real(dp) :: cpus
 real(dp) :: ddamp
 real(dp) :: dfpt_sciss
 real(dp) :: diecut
 real(dp) :: diegap
 real(dp) :: dielam
 real(dp) :: dielng
 real(dp) :: diemac
 real(dp) :: diemix
 real(dp) :: diemixmag
 real(dp) :: dilatmx
 real(dp) :: dmft_charge_prec
 real(dp) :: dmft_mxsf
 real(dp) :: dmft_tolfreq
 real(dp) :: dmft_tollc
 real(dp) :: dmftqmc_n
 real(dp) :: dosdeltae
 real(dp) :: dtion
 real(dp) :: dvdb_qcache_mb = 1024.0_dp
 real(dp) :: dvdb_qdamp = 0.1_dp
 real(dp) :: ecut
 real(dp) :: ecuteps
 real(dp) :: ecutsigx
 real(dp) :: ecutsm
 real(dp) :: ecutwfn
 real(dp) :: effmass_free
 real(dp) :: efmas_deg_tol
 real(dp) :: elph2_imagden
 real(dp) :: eshift
 real(dp) :: esmear
 real(dp) :: exchmix
 real(dp) :: fband
 real(dp) :: fermie_nest = zero
 real(dp) :: focktoldfe
 real(dp) :: freqim_alpha
 real(dp) :: freqremin = zero
 real(dp) :: freqremax = zero
 real(dp) :: freqspmin = zero
 real(dp) :: freqspmax = zero
 real(dp) :: friction
 real(dp) :: fxcartfactor
 real(dp) :: ga_opt_percent
 real(dp) :: gwencomp = 2.0_dp
 real(dp) :: gwls_model_parameter         ! Parameter used in dielectric function model
 real(dp) :: gw_toldfeig
 real(dp) :: hyb_mixing
 real(dp) :: hyb_mixing_sr
 real(dp) :: hyb_range_dft
 real(dp) :: hyb_range_fock
 real(dp) :: kptnrm
 real(dp) :: kptrlen
 real(dp) :: magcon_lambda
 real(dp) :: maxestep
 real(dp) :: mbpt_sciss = zero
 real(dp) :: mdf_epsinf = zero
 real(dp) :: mdwall
 real(dp) :: mep_mxstep
 real(dp) :: nelect
 real(dp) :: noseinert
 real(dp) :: omegasimax = 50/Ha_eV
 real(dp) :: omegasrdmax = 1.0_dp/Ha_eV  ! = 1eV
 real(dp) :: pawecutdg
 real(dp) :: pawovlp
 real(dp) :: pawujrad
 real(dp) :: pawujv
 real(dp) :: posocc
 real(dp) :: postoldfe
 real(dp) :: postoldff
 real(dp) :: ppmfrq = zero
 real(dp) :: pw_unbal_thresh
 real(dp) :: ratsm
 real(dp) :: ratsph_extra
 real(dp) :: recrcut
 real(dp) :: recefermi
 real(dp) :: rectolden
 real(dp) :: rhoqpmix
 real(dp) :: rcut
 real(dp) :: slabwsrad
 real(dp) :: slabzbeg
 real(dp) :: slabzend
 real(dp) :: spbroad
 real(dp) :: spinmagntarget
 real(dp) :: spnorbscl
 real(dp) :: stmbias
 real(dp) :: strfact
 real(dp) :: strprecon
 real(dp) :: td_maxene
 real(dp) :: tfw_toldfe
 real(dp) :: tl_radius
 real(dp) :: toldfe
 real(dp) :: tolmxde
 real(dp) :: toldff
 real(dp) :: tolimg
 real(dp) :: tolmxf
 real(dp) :: tolrde
 real(dp) :: tolrff
 real(dp) :: tolsym
 real(dp) :: tolvrs
 real(dp) :: tolwfr
 real(dp) :: tphysel
 real(dp) :: tsmear
 real(dp) :: userra = zero
 real(dp) :: userrb = zero
 real(dp) :: userrc = zero
 real(dp) :: userrd = zero
 real(dp) :: userre = zero
 real(dp) :: vacwidth
 real(dp) :: vdw_tol
 real(dp) :: vdw_tol_3bt
 real(dp) :: vdw_df_acutmin
 real(dp) :: vdw_df_aratio
 real(dp) :: vdw_df_damax
 real(dp) :: vdw_df_damin
 real(dp) :: vdw_df_dcut
 real(dp) :: vdw_df_dratio
 real(dp) :: vdw_df_dsoft
 real(dp) :: vdw_df_gcut
 real(dp) :: vdw_df_phisoft
 real(dp) :: vdw_df_qcut
 real(dp) :: vdw_df_qratio
 real(dp) :: vdw_df_rcut
 real(dp) :: vdw_df_rsoft
 real(dp) :: vdw_df_threshold = zero
 real(dp) :: vdw_df_tolerance
 real(dp) :: vdw_df_zab
 real(dp) :: vis
 real(dp) :: wfmix
 real(dp) :: wtq
 real(dp) :: wvl_hgrid
 real(dp) :: wvl_crmult
 real(dp) :: wvl_frmult
 real(dp) :: xc_denpos
 real(dp) :: xc_tb09_c
 real(dp) :: zcut

!Real arrays
 real(dp) :: boxcenter(3)
 real(dp) :: bfield(3)
 real(dp) :: dfield(3)
 real(dp) :: efield(3)
 real(dp) :: genafm(3)
 real(dp) :: goprecprm(3)
 real(dp) :: neb_spring(2)
 real(dp) :: pol(3)
 real(dp) :: polcen(3)
 real(dp) :: pvelmax(3)
 real(dp) :: qptn(3)
 real(dp) :: red_efield(3)
 real(dp) :: red_dfield(3)
 real(dp) :: red_efieldbar(3)
 real(dp) :: strtarget(6)
 real(dp) :: ucrpa_window(2)
 real(dp) :: vcutgeo(3)
 real(dp) :: vprtrb(2)
 real(dp) :: zeemanfield(3)
 real(dp) :: mdtemp(2)

!Real allocatables
 real(dp), allocatable :: acell_orig(:,:)   ! acell_orig(3,nimage)
 real(dp), allocatable :: amu_orig(:,:)     ! amu(ntypat,nimage)
 real(dp), allocatable :: atvshift(:,:,:)   ! atvshift(16,nsppol,natom)
 real(dp), allocatable :: cd_imfrqs(:)      ! cd_imfrqs(cd_customnimfrqs)
 real(dp), allocatable :: chrgat(:)         ! chrgat(natom)
 real(dp), allocatable :: chempot(:,:,:)    ! chempot(3,nzchempot,ntypat)
 real(dp), allocatable :: corecs(:)         ! corecs(ntypat)
 real(dp), allocatable :: densty(:,:)       ! densty(ntypat,4)
 real(dp), allocatable :: dmatpawu(:,:,:,:,:) ! dmatpawu(2*lpawu+1,2*lpawu+1,nsppol*nspinor,natpu,nimage)
                                              ! where natpu=number of atoms with lpawu/=1
 real(dp), allocatable :: efmas_dirs(:,:)   ! efmas_dirs(3,efmas_n_dirs)
 real(dp), allocatable :: f4of2_sla(:)      ! f4of2_sla(ntypat)
 real(dp), allocatable :: f6of2_sla(:)      ! f6of2_sla(ntypat)
 real(dp), allocatable :: gw_qlwl(:,:)      ! gw_qlwl(3,gw_nqlwl)
 real(dp), allocatable :: gw_freqsp(:)      ! gw_freqsp(gw_customnfreqsp)
 real(dp), allocatable :: gwls_list_proj_freq(:)      ! gwls_list_proj_freq(gwls_n_proj_freq)
 real(dp), allocatable :: jpawu(:,:)        ! jpawu(ntypat,nimage)
 real(dp), allocatable :: kpt(:,:)          ! kpt(3,nkpt)
 real(dp), allocatable :: kptgw(:,:)        ! kptgw(3,nkptgw)
 real(dp), allocatable :: kptns(:,:)        ! kptns(3,nkpt) k-points renormalized and shifted.
                                            !  The ones that should be used inside the code.
 real(dp), allocatable :: kptns_hf(:,:)     ! kpthf(3,nkptns_hf)

 real(dp), allocatable :: mixalch_orig(:,:,:) ! mixalch_orig(npspalch,ntypalch,nimage)
 real(dp), allocatable :: mixesimgf(:)        ! mixesimgf(nimage)
 real(dp), allocatable :: nucdipmom(:,:)      ! nucdipmom(3,natom)
 real(dp), allocatable :: occ_orig(:,:)       ! occ_orig(mband*nkpt*nsppol,nimage)
 real(dp), allocatable :: pimass(:)           ! pimass(ntypat)
 real(dp), allocatable :: ptcharge(:)         ! ptcharge(ntypat)
 real(dp), allocatable :: qmass(:)            ! qmass(nnos)
 real(dp), allocatable :: qptdm(:,:)          ! qptdm(3,nqptdm)
 real(dp), allocatable :: quadmom(:)          ! quadmom(ntypat)
 real(dp), allocatable :: ratsph(:)           ! ratsph(ntypat)
 real(dp), allocatable :: rprim_orig(:,:,:)   ! rprim_orig(3,3,nimage)
 real(dp), allocatable :: rprimd_orig(:,:,:)  ! rprimd_orig(3,3,nimage)
 real(dp), allocatable :: shiftk(:,:)         ! shiftk(3,nshiftk)
 real(dp) :: shiftk_orig(3,MAX_NSHIFTK)       ! original shifts given in input (changed in inkpts).

 real(dp), allocatable :: spinat(:,:)         ! spinat(3,natom)
 real(dp), allocatable :: tnons(:,:)          ! tnons(3,nsym)
 real(dp), allocatable :: upawu(:,:)          ! upawu(ntypat,nimage)
 real(dp), allocatable :: vel_cell_orig(:,:,:)! vel_cell_orig(3,3,nimage)
 real(dp), allocatable :: vel_orig(:,:,:)     ! vel_orig(3,natom,nimage)
 real(dp), allocatable :: wtatcon(:,:,:)      ! wtatcon(3,natom,nconeq)
 real(dp), allocatable :: wtk(:)              ! wtk(nkpt)
 real(dp), allocatable :: xred_orig(:,:,:)    ! xred_orig(3,natom,nimage)
 real(dp), allocatable :: xredsph_extra(:,:)  ! xredsph_extra(3,natsph_extra)
 real(dp), allocatable :: ziontypat(:)        ! ziontypat(ntypat)
 real(dp), allocatable :: znucl(:)            ! znucl(npsp)

!BEGIN VARIABLES FOR @Bethe-Salpeter
 integer :: bs_algorithm = 2
 integer :: bs_haydock_niter = 100
 integer :: bs_exchange_term = 1
 integer :: bs_coulomb_term = 11
 integer :: bs_calctype = 1
 integer :: bs_coupling = 0
 integer :: bs_interp_mode = 0 ! No interpolation
 integer :: bs_interp_prep = 0 ! Do not prepare interp
 integer :: bs_interp_method = 1 ! YG interpolation
 integer :: bs_interp_rl_nb = 1

 real(dp) :: bs_interp_m3_width = one

 integer  :: bs_interp_kmult(3) = 0
 real(dp) :: bs_haydock_tol(2) = [0.02_dp, zero]

 integer,allocatable :: bs_loband(:)

!  Take big absolute value numbers, but the the biggest ones, otherwise overflow can happen
 real(dp) :: bs_eh_cutoff(2) = [smallest_real*tol6, greatest_real*tol6]
 real(dp) :: bs_freq_mesh(3) = [zero,zero, 0.01_dp/Ha_eV]
!END VARIABLES FOR @Bethe-Salpeter.

 integer :: gpu_linalg_limit

!EPH variables
! ifc variables
 integer :: asr = 1
 integer :: dipdip = 1
 integer :: chneut = 1
 integer :: symdynmat = 1

! Phonon variables.
 integer :: ph_ndivsm = 20
 integer :: ph_nqpath = 0
 integer :: ph_ngqpt(3) = 20
 integer :: ph_nqshift = 1

 real(dp),allocatable :: ph_qshift(:,:)
  ! ph_qshift(3, ph_nqshift)
 real(dp),allocatable :: ph_qpath(:,:)
  ! ph_qpath(3, nqpath)

! e-ph variables
 real(dp) :: eph_mustar = 0.1_dp
 integer :: eph_intmeth = 2
 real(dp) :: eph_extrael = zero
 real(dp) :: eph_fermie = zero
 integer :: eph_frohlichm = 0
 real(dp) :: eph_fsmear = 0.01_dp
 real(dp) :: eph_fsewin = 0.04_dp
 real(dp) :: eph_ecutosc = zero
 !real(dp) :: eph_alpha_gmin = zero !sqrt(5)
 real(dp) :: eph_tols_idelta(2) = [tol12, tol12]
 integer :: eph_phrange(2) = 0

 integer :: eph_ngqpt_fine(3) = 0
 integer :: eph_np_pqbks(5) = 0

 integer :: eph_stern = 0
 integer :: eph_transport = 0
 integer :: eph_use_ftinterp = 0

 integer :: ph_intmeth = 2
 integer :: prteliash = 0
 real(dp) :: ph_wstep = 0.1_dp / Ha_meV
 real(dp) :: ph_smear = 0.00002_dp
 integer :: ddb_ngqpt(3) = 0
 real(dp) :: ddb_shiftq(3) = zero

 integer :: mixprec = 0
 integer :: symv1scf = 0
 integer :: dvdb_add_lr = 1
 integer :: dvdb_rspace_cell = 0

 integer :: sigma_bsum_range(2) = 0

 real(dp) :: sigma_erange(2) = zero

 integer :: transport_ngkpt(3) = 0
 ! K-mesh for Transport calculation.

 integer :: sigma_ngkpt(3) = 0
 ! K-mesh for Sigma_{nk} (only IBZ points). Alternative to kptgw.

 integer :: sigma_nshiftk = 1
 ! Number of shifts in k-mesh for Sigma_{nk}.

 real(dp),allocatable :: sigma_shiftk(:,:)
 ! sigma_shiftk(3, sigma_nshiftk)
 ! shifts in k-mesh for Sigma_{nk}.
!END EPH

 integer :: ndivsm = 0
 integer :: nkpath = 0
 real(dp) :: einterp(4) = zero
 real(dp),allocatable :: kptbounds(:,:)
 real(dp) :: tmesh(3) = [5._dp, 59._dp, 6._dp]

 character(len=fnlen) :: getddb_filepath = ABI_NOFILE
 character(len=fnlen) :: getden_filepath = ABI_NOFILE
 character(len=fnlen) :: getdvdb_filepath = ABI_NOFILE
 character(len=fnlen) :: getwfk_filepath = ABI_NOFILE
 character(len=fnlen) :: getwfkfine_filepath = ABI_NOFILE
 character(len=fnlen) :: getwfq_filepath = ABI_NOFILE
 character(len=fnlen) :: getkerange_filepath = ABI_NOFILE
 character(len=fnlen) :: getpot_filepath = ABI_NOFILE
 character(len=fnlen) :: getscr_filepath = ABI_NOFILE
 character(len=fnlen) :: getsigeph_filepath = ABI_NOFILE

 contains

 procedure :: chkneu => dtset_chkneu
   ! Check neutrality of system based on band occupancies and valence charges of pseudo-atoms.

 procedure :: copy => dtset_copy
   ! Copy object.

 procedure :: free => dtset_free
   ! Free dynamic memory.

 procedure :: free_nkpt_arrays => dtset_free_nkpt_arrays
   ! Free arrays that depend on input nkpt (used in EPH code)

 procedure :: get_npert_rbz => dtset_get_npert_rbz
   ! Get the number of effective pertubation done in looper3, nkpt_rbz, nband_rbz

 procedure :: testsusmat => dtset_testsusmat
   ! Test wether a new susceptibility matrix and/or a new dielectric matrix must be computed

 end type dataset_type
!!***

 public :: find_getdtset           ! Find the number of the dataset (iget) for a given value of a "get" variable (getvalue)
 public :: macroin                 ! Treat "macro" input variables
 public :: macroin2
 public :: chkvars                 ! Examines the input string, to check whether all names are allowed.

CONTAINS  !==============================================================================
!!***

!!****f* m_dtset/dtset_chkneu
!! NAME
!! dtset_chkneu
!!
!! FUNCTION
!! Check neutrality of system based on band occupancies and valence charges of pseudo-atoms.
!! Eventually initialize occ if occopt==1 or 3...8
!! Also return nelect, the number of valence electron per unit cell
!!
!! INPUTS
!!  charge=number of electrons missing (+) or added (-) to system (usually 0)
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | iscf= if>0, SCF calculation ; if<=0, non SCF calculation (wtk might
!!   |  not be defined)
!!   | natom=number of atoms in unit cell
!!   | nband(nkpt*nsppol)=number of bands at each k point
!!   | nkpt=number of k points
!!   | nspinor=number of spinorial components of the wavefunctions
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | ntypat=number of pseudopotentials
!!   | positron=0 if electron GS calculation
!!   |          1 if positron GS calculation
!!   |          2 if electron GS calcultaion in presence of the positron
!!   | typat(natom)=atom type (integer) for each atom
!!   | wtk(nkpt)=k point weights (defined if iscf>0 or iscf==-3)
!!   | ziontypat(ntypat)=ionic charge of each pseudoatom
!!  occopt=option for occupancies
!!
!! OUTPUT
!!  Writes warning and/or aborts if error condition exists
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | nelect=number of valence electrons per unit cell
!!   |  (from counting valence electrons in psps, and taking into
!!   |   account the input variable "charge")
!!
!! SIDE EFFECTS
!! Input/Output :
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | occ_orig(dtset%nband(1)*nkpt*nsppol,nimage)=occupation numbers for each band and k point
!!   |   must be input for occopt==0 or 2,
!!   |   will be an output for occopt==1 or 3 ... 8
!!
!! PARENTS
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine dtset_chkneu(dtset, charge, occopt)

!Arguments ------------------------------------
!scalars
 class(dataset_type),intent(inout) :: dtset
 integer,intent(in) :: occopt
 real(dp),intent(in) :: charge

!Local variables-------------------------------
!scalars
 integer :: bantot,iatom,iband,ii,iimage,ikpt,isppol,nocc
 real(dp) :: maxocc,nelect_occ,nelect_spin,occlast,sign_spin,zval
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: tmpocc(:)

! *************************************************************************

!(1) count nominal valence electrons according to ziontypat
 zval=zero
 do iatom=1,dtset%natom
   zval=zval+dtset%ziontypat(dtset%typat(iatom))
 end do
 if (dtset%positron/=1) then
   dtset%nelect=zval-charge
 else
   dtset%nelect=one
 end if

! write(std_out,*)ch10,' chkneu : enter, dtset%nelect=',dtset%nelect
! write(std_out,*)' occopt,dtset%nsppol,dtset%nspden=',occopt,dtset%nsppol,dtset%nspden

!(2) Optionally initialize occ with semiconductor occupancies
!(even for a metal : at this stage, the eigenenergies are unknown)
!Note that nband(1)=nband(2) in this section, as occopt=2 is avoided.
 if(occopt==1 .or. (occopt>=3 .and. occopt<=8) )then
!  Here, initialize a real(dp) variable giving the
!  maximum occupation number per band
   maxocc=2.0_dp/real(dtset%nsppol*dtset%nspinor,dp)

!  Determine the number of bands fully or partially occupied
   nocc=int((dtset%nelect-1.0d-8)/maxocc) + 1
!  Occupation number of the highest level
   occlast=dtset%nelect-maxocc*(nocc-1)
   !write(std_out,*)' maxocc,nocc,occlast=',maxocc,nocc,occlast

!  The number of allowed bands must be sufficiently large
   if( nocc<=dtset%nband(1)*dtset%nsppol .or. dtset%iscf==-2) then

     if(dtset%iscf==-2 .and. nocc>dtset%nband(1)*dtset%nsppol)nocc=dtset%nband(1)*dtset%nsppol

!    First treat the case where the spin magnetization is not imposed, is zero with nspden==1, or has sufficient flexibility
!    for a target not to be matched at the initialisation, but later
     if(abs(dtset%spinmagntarget+99.99_dp)<tol8 .or. (dtset%nspden==4) .or. &
       (abs(dtset%spinmagntarget)<tol8.and.dtset%nspden==1))then

!      Use a temporary array for defining occupation numbers
       ABI_MALLOC(tmpocc,(dtset%nband(1)*dtset%nsppol))
!      First do it for fully occupied bands
       if (1<nocc) tmpocc(1:nocc-1)=maxocc
!      Then, do it for highest occupied band
       if (1<=nocc) tmpocc(nocc)=occlast
!      Finally do it for eventual unoccupied bands
       if ( nocc<dtset%nband(1)*dtset%nsppol ) tmpocc(nocc+1:dtset%nband(1)*dtset%nsppol)=0.0_dp

!      Now copy the tmpocc array in the occ array, taking into account the spin
       if(dtset%nsppol==1)then
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             dtset%occ_orig(iband+(ikpt-1)*dtset%nband(1),:)=tmpocc(iband)
           enddo
         end do
       else
         do ikpt=1,dtset%nkpt
           do iband=1,dtset%nband(1)
             do isppol=1,dtset%nsppol
               dtset%occ_orig(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)),:) =  &
&               tmpocc(isppol+dtset%nsppol*(iband-1))
             end do
           end do
         end do
       end if
       ABI_FREE(tmpocc)

     ! Second, treat the case in which one imposes the spin magnetization (only possible for nspden==2)
     ! Also treat antiferromagnetic case (nsppol==1, nspden==2), although spinmagntarget must be zero
     else if(abs(dtset%spinmagntarget+99.99_dp)>1.0d-10)then
       do isppol=1,dtset%nsppol
         sign_spin=real(3-2*isppol,dp)
         nelect_spin=half*(dtset%nelect*maxocc+sign_spin*dtset%spinmagntarget)

         !write(std_out,*)' isppol,sign_spin,nelect_spin=',isppol,sign_spin,nelect_spin
         ! Determines the last state, and its occupation
         if(abs(nint(nelect_spin)-nelect_spin)<tol10)then
           nocc=nint(nelect_spin/maxocc)
           occlast=maxocc
         else
           nocc=ceiling(nelect_spin/maxocc)
           occlast=nelect_spin-(real(nocc,dp)-one)*maxocc
         end if
         !write(std_out,*)' dtset%nband(1),maxocc,occlast=',dtset%nband(1),maxocc,occlast
         if(dtset%nband(1)*nint(maxocc)<nocc)then
           write(msg, '(a,i0,a, a,2i0,a, a,es16.6,a, a,es16.6,6a)' )&
           'Initialization of occ, with nspden = ',dtset%nspden,ch10,&
           'number of bands = ',dtset%nband(1:2),ch10,&
           'number of electrons = ',dtset%nelect,ch10,&
           'and spinmagntarget = ',dtset%spinmagntarget,ch10,&
           'This combination is not possible, because of a lack of bands.',ch10,&
           'Action: modify input file ... ',ch10,&
           '(you should likely increase nband, but also check nspden, nspinor, nsppol, and spinmagntarget)'
           MSG_ERROR(msg)
         end if
         do ikpt=1,dtset%nkpt
           ! Fill all bands, except the upper one
           if(dtset%nband(1)>1)then
             do iband=1,nocc-1
               dtset%occ_orig(iband+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)),:)=maxocc
             end do
           end if
           ! Fill the upper occupied band
           dtset%occ_orig(nocc+dtset%nband(1)*(ikpt-1+dtset%nkpt*(isppol-1)),:)=occlast
         end do
       end do

     else
       write(msg, '(a,i0,a,a,es16.6,6a)' )&
       'Initialization of occ, with nspden = ',dtset%nspden,ch10,&
       'and spinmagntarget = ',dtset%spinmagntarget,ch10,&
       'This combination is not possible.',ch10,&
       'Action: modify input file ... ',ch10,&
       '(check nspden, nspinor, nsppol and spinmagntarget)'
       MSG_ERROR(msg)
     end if

     ! Now print the values (only the first image, since they are all the same)
     if(dtset%nsppol==1)then
       if (dtset%prtvol > 0) then
         write(msg, '(a,i0,a,a)' ) &
          ' chkneu: initialized the occupation numbers for occopt= ',occopt,', spin-unpolarized or antiferromagnetic case:'
         call wrtout(std_out,msg)
         do ii=0,(dtset%nband(1)-1)/12
           write(msg,'(12f6.2)') dtset%occ_orig( 1+ii*12 : min(12+ii*12,dtset%nband(1)),1 )
           call wrtout(std_out,msg)
         end do
       end if
     else
       write(msg, '(a,i0,2a)' ) &
        ' dtset_chkneu: initialized the occupation numbers for occopt= ',occopt,ch10,'    spin up   values:'
       call wrtout(std_out,msg)
       if (dtset%prtvol > 0) then
         do ii=0,(dtset%nband(1)-1)/12
           write(msg,'(12f6.2)') dtset%occ_orig( 1+ii*12 : min(12+ii*12,dtset%nband(1)),1 )
           call wrtout(std_out,msg)
         end do
         call wrtout(std_out,'    spin down values:')
         do ii=0,(dtset%nband(1)-1)/12
           write(msg,'(12f6.2)') &
             dtset%occ_orig( 1+ii*12+dtset%nkpt*dtset%nband(1) : min(12+ii*12,dtset%nband(1))+dtset%nkpt*dtset%nband(1) ,1)
           call wrtout(std_out,msg)
         end do
       end if

     end if

!    Here, treat the case when the number of allowed bands is not large enough
   else
     write(msg, '(a,i0,8a)' )&
     'Initialization of occ, with occopt: ',occopt,ch10,&
     'There are not enough bands to get charge balance right',ch10,&
     'Action: modify input file ... ',ch10,&
     '(check the pseudopotential charges, the variable charge,',ch10,&
     'and the declared number of bands, nband)'
     MSG_ERROR(msg)
   end if
 end if

!The remaining of the routine is for SCF runs and special options
 if(dtset%iscf>0 .or. dtset%iscf==-1 .or. dtset%iscf==-3)then

   do iimage=1,dtset%nimage

!    (3) count electrons in bands (note : in case occ has just been
!    initialized, point (3) and (4) is a trivial test
     nelect_occ=0.0_dp
     bantot=0
     do isppol=1,dtset%nsppol
       do ikpt=1,dtset%nkpt
         do iband=1,dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
           bantot=bantot+1
           nelect_occ=nelect_occ+dtset%wtk(ikpt)*dtset%occ_orig(bantot,iimage)
         end do
       end do
     end do

!    (4) if dtset%iscf/=-3, dtset%nelect must equal nelect_occ
!    if discrepancy exceeds tol11, give warning;  tol8, stop with error

     if (abs(nelect_occ-dtset%nelect)>tol11 .and. dtset%iscf/=-3) then

!      There is a discrepancy
       write(msg, &
       '(a,a,e16.8,a,e16.8,a,a,a,e22.14,a,a,a,i5,a,a,a,a)' ) ch10,&
       ' chkneu: nelect_occ=',nelect_occ,', zval=',zval,',',ch10,&
       '         and input value of charge=',charge,',',ch10,&
       '   nelec_occ is computed from occ and wtk, iimage=',iimage,ch10,&
       '   zval is nominal charge of all nuclei, computed from zion (read in psp),',ch10,&
       '   charge is an input variable (usually 0).'
       call wrtout(std_out,msg)

       if (abs(nelect_occ-dtset%nelect)>tol8) then
!        The discrepancy is severe
         write(msg,'(a,a,e9.2,a,a)')ch10,&
         'These must obey zval-nelect_occ=charge to better than ',tol8,ch10,&
         ' This is not the case. '
       else
!        The discrepancy is not so severe
         write(msg, '(2a,e9.2)' )ch10,'These should obey zval-nelect_occ=charge to better than: ',tol11
       end if
       MSG_WARNING(msg)

       write(msg, '(6a)' ) &
       'Action: check input file for occ,wtk, and charge.',ch10,&
       'Note that wtk is NOT automatically normalized when occopt=2,',ch10,&
       'but IS automatically normalized otherwise.',ch10
       call wrtout(std_out,msg)

       ! If the discrepancy is severe, stop
       if (abs(nelect_occ-dtset%nelect)>tol8)then
         MSG_ERROR(msg)
       end if

     end if
   end do

 end if ! condition dtset%iscf>0 or -1 or -3 .

end subroutine dtset_chkneu
!!***

!----------------------------------------------------------------------

!!****f* m_dtset/dtset_copy
!! NAME
!! dtset_copy
!!
!! FUNCTION
!! Copy all values of dataset dtin to dataset dtout. allocatables of dtout are
!! allocated if required. Use dtset_free() to free a dataset after use.
!!
!! INPUTS
!!  dtin <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!  dtout <type(dataset_type)>
!!
!! PARENTS
!!      chkinp,dfpt_looppert,driver,gwls_hamiltonian,m_io_kss
!!
!! CHILDREN
!!
!! SOURCE

type(dataset_type) function dtset_copy(dtin) result(dtout)

!Arguments ------------------------------------
!scalars
 class(dataset_type),intent(in) :: dtin

! *************************************************************************

 DBG_ENTER("COLL")

 !@dataset_type

!BEGIN VARIABLES FOR @Bethe-Salpeter
 dtout%bs_algorithm     = dtin%bs_algorithm
 dtout%bs_haydock_niter = dtin%bs_haydock_niter
 dtout%bs_exchange_term = dtin%bs_exchange_term
 dtout%bs_coulomb_term  = dtin%bs_coulomb_term
 dtout%bs_calctype      = dtin%bs_calctype
 dtout%bs_coupling      = dtin%bs_coupling
 dtout%bs_haydock_tol   = dtin%bs_haydock_tol
 dtout%bs_hayd_term     = dtin%bs_hayd_term
 dtout%bs_interp_m3_width = dtin%bs_interp_m3_width
 dtout%bs_interp_method = dtin%bs_interp_method
 dtout%bs_interp_mode   = dtin%bs_interp_mode
 dtout%bs_interp_prep   = dtin%bs_interp_prep
 dtout%bs_interp_rl_nb  = dtin%bs_interp_rl_nb
 dtout%bs_interp_kmult(:) = dtin%bs_interp_kmult(:)
 dtout%bs_eh_cutoff(:) = dtin%bs_eh_cutoff(:)
 dtout%bs_freq_mesh(:) = dtin%bs_freq_mesh(:)
!END VARIABLES FOR @Bethe-Salpeter.

!Copy integers from dtin to dtout
 dtout%iomode             = dtin%iomode
 dtout%accuracy           = dtin%accuracy
 dtout%adpimd             = dtin%adpimd
 dtout%autoparal          = dtin%autoparal
 dtout%auxc_ixc           = dtin%auxc_ixc
 dtout%auxc_scal          = dtin%auxc_scal
 dtout%awtr               = dtin%awtr
 dtout%bandpp             = dtin%bandpp
 dtout%bdeigrf            = dtin%bdeigrf
 dtout%berryopt           = dtin%berryopt
 dtout%berrysav           = dtin%berrysav
 dtout%berrystep          = dtin%berrystep
 dtout%brav               = dtin%brav
 dtout%brvltt             = dtin%brvltt
 dtout%bs_nstates         = dtin%bs_nstates
 dtout%builtintest        = dtin%builtintest
 dtout%cd_customnimfrqs   = dtin%cd_customnimfrqs
 dtout%cd_frqim_method    = dtin%cd_frqim_method
 dtout%cd_full_grid       = dtin%cd_full_grid
 dtout%chkdilatmx         = dtin%chkdilatmx
 dtout%chkexit            = dtin%chkexit
 dtout%chkprim            = dtin%chkprim
 dtout%chksymbreak        = dtin%chksymbreak
 dtout%cineb_start        = dtin%cineb_start
 dtout%delayperm          = dtin%delayperm
 dtout%diismemory         = dtin%diismemory
 dtout%dmatpuopt          = dtin%dmatpuopt
 dtout%dmatudiag          = dtin%dmatudiag
 dtout%dmft_dc            = dtin%dmft_dc
 dtout%dmft_entropy       = dtin%dmft_entropy
 dtout%dmft_charge_prec   = dtin%dmft_charge_prec
 dtout%dmft_iter          = dtin%dmft_iter
 dtout%dmft_kspectralfunc = dtin%dmft_kspectralfunc
 dtout%dmft_nlambda       = dtin%dmft_nlambda
 dtout%dmft_mxsf          = dtin%dmft_mxsf
 dtout%dmft_nwlo          = dtin%dmft_nwlo
 dtout%dmft_nwli          = dtin%dmft_nwli
 dtout%dmft_occnd_imag    = dtin%dmft_occnd_imag
 dtout%dmft_read_occnd    = dtin%dmft_read_occnd
 dtout%dmft_rslf          = dtin%dmft_rslf
 dtout%dmft_solv          = dtin%dmft_solv
 dtout%dmft_t2g           = dtin%dmft_t2g
!dtout%dmft_x2my2d        = dtin%dmft_x2my2d
 dtout%dmft_tolfreq       = dtin%dmft_tolfreq
 dtout%dmft_tollc         = dtin%dmft_tollc
 dtout%dmftbandi          = dtin%dmftbandi
 dtout%dmftbandf          = dtin%dmftbandf
 dtout%dmftcheck          = dtin%dmftcheck
 dtout%dmftctqmc_basis    = dtin%dmftctqmc_basis
 dtout%dmftctqmc_check    = dtin%dmftctqmc_check
 dtout%dmftctqmc_correl   = dtin%dmftctqmc_correl
 dtout%dmftctqmc_gmove    = dtin%dmftctqmc_gmove
 dtout%dmftctqmc_grnns    = dtin%dmftctqmc_grnns
 dtout%dmftctqmc_meas     = dtin%dmftctqmc_meas
 dtout%dmftctqmc_mrka     = dtin%dmftctqmc_mrka
 dtout%dmftctqmc_mov      = dtin%dmftctqmc_mov
 dtout%dmftctqmc_order    = dtin%dmftctqmc_order
 dtout%dmftctqmc_triqs_nleg = dtin%dmftctqmc_triqs_nleg
 dtout%dmftqmc_n          = dtin%dmftqmc_n
 dtout%dmftqmc_l          = dtin%dmftqmc_l
 dtout%dmftqmc_seed       = dtin%dmftqmc_seed
 dtout%dmftqmc_therm      = dtin%dmftqmc_therm
 dtout%d3e_pert1_elfd     = dtin%d3e_pert1_elfd
 dtout%d3e_pert1_phon     = dtin%d3e_pert1_phon
 dtout%d3e_pert2_elfd     = dtin%d3e_pert2_elfd
 dtout%d3e_pert2_phon     = dtin%d3e_pert2_phon
 dtout%d3e_pert2_strs     = dtin%d3e_pert2_strs
 dtout%d3e_pert3_elfd     = dtin%d3e_pert3_elfd
 dtout%d3e_pert3_phon     = dtin%d3e_pert3_phon
 dtout%efmas              = dtin%efmas
 dtout%efmas_calc_dirs    = dtin%efmas_calc_dirs
 dtout%efmas_deg          = dtin%efmas_deg
 dtout%efmas_dim          = dtin%efmas_dim
 dtout%efmas_n_dirs       = dtin%efmas_n_dirs
 dtout%efmas_ntheta       = dtin%efmas_ntheta
 dtout%enunit             = dtin%enunit

! begin eph variables
 dtout%asr                = dtin%asr
 dtout%dipdip             = dtin%dipdip
 dtout%chneut             = dtin%chneut

 dtout%eph_mustar         = dtin%eph_mustar
 dtout%eph_intmeth        = dtin%eph_intmeth
 dtout%eph_tols_idelta    = dtin%eph_tols_idelta
 dtout%eph_phrange        = dtin%eph_phrange
 dtout%eph_extrael        = dtin%eph_extrael
 dtout%eph_fermie         = dtin%eph_fermie
 dtout%eph_frohlichm      = dtin%eph_frohlichm
 dtout%eph_fsmear         = dtin%eph_fsmear
 dtout%eph_fsewin         = dtin%eph_fsewin
 dtout%eph_ecutosc        = dtin%eph_ecutosc
 !dtout%eph_alpha_gmin     = dtin%eph_alpha_gmin
 dtout%eph_ngqpt_fine     = dtin%eph_ngqpt_fine
 dtout%eph_np_pqbks       = dtin%eph_np_pqbks
 dtout%eph_restart        = dtin%eph_restart
 dtout%eph_task           = dtin%eph_task
 dtout%eph_stern          = dtin%eph_stern
 dtout%eph_use_ftinterp   = dtin%eph_use_ftinterp
 dtout%eph_transport      = dtin%eph_transport

 dtout%ph_wstep          = dtin%ph_wstep
 dtout%ph_intmeth        = dtin%ph_intmeth
 dtout%symdynmat         = dtin%symdynmat
 dtout%symv1scf          = dtin%symv1scf
 dtout%ph_nqshift        = dtin%ph_nqshift
 if (allocated(dtin%ph_qshift)) call alloc_copy(dtin%ph_qshift, dtout%ph_qshift)
 dtout%ph_smear          = dtin%ph_smear
 dtout%ddb_ngqpt         = dtin%ddb_ngqpt
 dtout%ddb_shiftq        = dtin%ddb_shiftq
 dtout%dvdb_qcache_mb    = dtin%dvdb_qcache_mb
 dtout%dvdb_qdamp        = dtin%dvdb_qdamp
 dtout%dvdb_add_lr       = dtin%dvdb_add_lr
 dtout%dvdb_rspace_cell  = dtin%dvdb_rspace_cell

 dtout%sigma_bsum_range = dtin%sigma_bsum_range
 dtout%sigma_erange = dtin%sigma_erange
 dtout%sigma_ngkpt = dtin%sigma_ngkpt
 dtout%sigma_nshiftk = dtin%sigma_nshiftk
 if (allocated(dtin%sigma_shiftk)) call alloc_copy(dtin%sigma_shiftk, dtout%sigma_shiftk)

 dtout%transport_ngkpt = dtin%transport_ngkpt

 dtout%ph_ndivsm          = dtin%ph_ndivsm
 dtout%ph_nqpath          = dtin%ph_nqpath
 dtout%ph_ngqpt           = dtin%ph_ngqpt
 if (allocated(dtin%ph_qpath)) call alloc_copy(dtin%ph_qpath, dtout%ph_qpath)
! end eph variables

 dtout%exchn2n3d          = dtin%exchn2n3d
 dtout%extrapwf           = dtin%extrapwf
 dtout%pawfatbnd          = dtin%pawfatbnd
 dtout%fermie_nest        = dtin%fermie_nest
 dtout%fftgw              = dtin%fftgw
 dtout%fockdownsampling   = dtin%fockdownsampling
 dtout%fockoptmix         = dtin%fockoptmix
 dtout%fock_icutcoul      = dtin%fock_icutcoul
 dtout%freqim_alpha       = dtin%freqim_alpha
 dtout%freqremin          = dtin%freqremin
 dtout%freqremax          = dtin%freqremax
 dtout%freqspmin          = dtin%freqspmin
 dtout%freqspmax          = dtin%freqspmax
 dtout%frzfermi           = dtin%frzfermi
 dtout%ga_algor           = dtin%ga_algor
 dtout%ga_fitness         = dtin%ga_fitness
 dtout%ga_n_rules         = dtin%ga_n_rules
 dtout%getbseig           = dtin%getbseig
 dtout%getbsreso          = dtin%getbsreso
 dtout%getbscoup          = dtin%getbscoup
 dtout%getcell            = dtin%getcell
 dtout%getddb             = dtin%getddb
 dtout%getdvdb            = dtin%getdvdb
 dtout%getddk             = dtin%getddk
 dtout%getdelfd           = dtin%getdelfd
 dtout%getdkdk            = dtin%getdkdk
 dtout%getdkde            = dtin%getdkde
 dtout%getden             = dtin%getden
 dtout%getefmas           = dtin%getefmas
 dtout%getgam_eig2nkq     = dtin%getgam_eig2nkq
 dtout%gethaydock         = dtin%gethaydock
 dtout%getocc             = dtin%getocc
 dtout%getpawden          = dtin%getpawden
 dtout%getddb_filepath        = dtin%getddb_filepath
 dtout%getden_filepath        = dtin%getden_filepath
 dtout%getdvdb_filepath       = dtin%getdvdb_filepath
 dtout%getpot_filepath        = dtin%getpot_filepath
 dtout%getsigeph_filepath     = dtin%getsigeph_filepath
 dtout%getscr_filepath        = dtin%getscr_filepath
 dtout%getwfk_filepath        = dtin%getwfk_filepath
 dtout%getwfkfine_filepath    = dtin%getwfkfine_filepath
 dtout%getwfq_filepath        = dtin%getwfq_filepath
 dtout%getqps             = dtin%getqps
 dtout%getscr             = dtin%getscr
 dtout%getsuscep          = dtin%getsuscep
 dtout%getvel             = dtin%getvel
 dtout%getwfk             = dtin%getwfk
 dtout%getwfkfine         = dtin%getwfkfine
 dtout%getwfq             = dtin%getwfq
 dtout%getxcart           = dtin%getxcart
 dtout%getxred            = dtin%getxred
 dtout%get1den            = dtin%get1den
 dtout%get1wf             = dtin%get1wf
 dtout%goprecon           = dtin%goprecon
 dtout%gpu_linalg_limit   = dtin%gpu_linalg_limit
 dtout%gwaclowrank        = dtin%gwaclowrank
 dtout%gwcalctyp          = dtin%gwcalctyp
 dtout%gwcomp             = dtin%gwcomp
 dtout%gwencomp           = dtin%gwencomp
 dtout%gwmem              = dtin%gwmem
 dtout%gwpara             = dtin%gwpara
 dtout%gwgamma            = dtin%gwgamma
 dtout%gwrpacorr          = dtin%gwrpacorr
 dtout%gw_customnfreqsp   = dtin%gw_customnfreqsp
 dtout%gw_icutcoul        = dtin%gw_icutcoul
 dtout%gw_nqlwl           = dtin%gw_nqlwl
 dtout%gw_nstep           = dtin%gw_nstep
 dtout%gw_frqim_inzgrid   = dtin%gw_frqim_inzgrid
 dtout%gw_frqre_inzgrid   = dtin%gw_frqre_inzgrid
 dtout%gw_frqre_tangrid   = dtin%gw_frqre_tangrid
 dtout%gw_invalid_freq    = dtin%gw_invalid_freq
 dtout%gw_qprange         = dtin%gw_qprange
 dtout%gw_sctype          = dtin%gw_sctype
 dtout%gw_sigxcore        = dtin%gw_sigxcore
 dtout%gw_toldfeig        = dtin%gw_toldfeig
 dtout%gwls_stern_kmax= dtin%gwls_stern_kmax
 dtout%gwls_npt_gauss_quad  = dtin%gwls_npt_gauss_quad
 dtout%gwls_diel_model= dtin%gwls_diel_model
 dtout%gwls_print_debug     = dtin%gwls_print_debug
 dtout%gwls_nseeds          = dtin%gwls_nseeds
 dtout%gwls_n_proj_freq     = dtin%gwls_n_proj_freq
 dtout%gwls_kmax_complement = dtin%gwls_kmax_complement
 dtout%gwls_kmax_poles      = dtin%gwls_kmax_poles
 dtout%gwls_kmax_analytic   = dtin%gwls_kmax_analytic
 dtout%gwls_kmax_numeric    = dtin%gwls_kmax_numeric
 dtout%gwls_band_index      = dtin%gwls_band_index
 dtout%gwls_exchange        = dtin%gwls_exchange
 dtout%gwls_correlation     = dtin%gwls_correlation
 dtout%gwls_first_seed      = dtin%gwls_first_seed
 dtout%gwls_recycle         = dtin%gwls_recycle
 dtout%hyb_mixing      = dtin%hyb_mixing
 dtout%hyb_mixing_sr   = dtin%hyb_mixing_sr
 dtout%hyb_range_dft   = dtin%hyb_range_dft
 dtout%hyb_range_fock  = dtin%hyb_range_fock
 dtout%hmcsst             = dtin%hmcsst
 dtout%hmctt              = dtin%hmctt
 dtout%iboxcut            = dtin%iboxcut
 dtout%icoulomb           = dtin%icoulomb
 dtout%icutcoul           = dtin%icutcoul
 dtout%ieig2rf            = dtin%ieig2rf
 dtout%imgmov             = dtin%imgmov
 dtout%imgwfstor          = dtin%imgwfstor
 dtout%inclvkb            = dtin%inclvkb
 dtout%intxc              = dtin%intxc
 dtout%ionmov             = dtin%ionmov
 dtout%densfor_pred       = dtin%densfor_pred
 dtout%iprcel             = dtin%iprcel
 dtout%iprcfc             = dtin%iprcfc
 dtout%irandom            = dtin%irandom
 dtout%irdbseig           = dtin%irdbseig
 dtout%irdbsreso          = dtin%irdbsreso
 dtout%irdbscoup          = dtin%irdbscoup
 dtout%irdddb             = dtin%irdddb
 dtout%irddvdb            = dtin%irddvdb
 dtout%irdddk             = dtin%irdddk
 dtout%irdden             = dtin%irdden
 dtout%irdefmas           = dtin%irdefmas
 dtout%irdhaydock         = dtin%irdhaydock
 dtout%irdpawden          = dtin%irdpawden
 dtout%irdqps             = dtin%irdqps
 dtout%irdscr             = dtin%irdscr
 dtout%irdsuscep          = dtin%irdsuscep
 dtout%irdvdw             = dtin%irdvdw
 dtout%irdwfk             = dtin%irdwfk
 dtout%irdwfkfine         = dtin%irdwfkfine
 dtout%irdwfq             = dtin%irdwfq
 dtout%ird1den            = dtin%ird1den
 dtout%ird1wf             = dtin%ird1wf
 dtout%iscf               = dtin%iscf
 dtout%isecur             = dtin%isecur
 dtout%istatimg           = dtin%istatimg
 dtout%istatr             = dtin%istatr
 dtout%istatshft          = dtin%istatshft
 dtout%ixc                = dtin%ixc
 dtout%ixc_sigma          = dtin%ixc_sigma
 dtout%ixcpositron        = dtin%ixcpositron
 dtout%ixcrot             = dtin%ixcrot
 dtout%jdtset             = dtin%jdtset
 dtout%jellslab           = dtin%jellslab
 dtout%kptopt             = dtin%kptopt
 dtout%kssform            = dtin%kssform
 dtout%localrdwf          = dtin%localrdwf
#if defined HAVE_LOTF
 dtout%lotf_classic       = dtin%lotf_classic
 dtout%lotf_nitex         = dtin%lotf_nitex
 dtout%lotf_nneigx        = dtin%lotf_nneigx
 dtout%lotf_version       = dtin%lotf_version
#endif
 dtout%lw_flexo           = dtin%lw_flexo
 dtout%lw_qdrpl           = dtin%lw_qdrpl
 dtout%magconon           = dtin%magconon
 dtout%maxnsym            = dtin%maxnsym
 dtout%max_ncpus          = dtin%max_ncpus
 dtout%mband              = dtin%mband
 dtout%mdf_epsinf         = dtin%mdf_epsinf
 dtout%mep_solver         = dtin%mep_solver
 dtout%mem_test           = dtin%mem_test
 dtout%mixprec            = dtin%mixprec
 dtout%mffmem             = dtin%mffmem
 dtout%mgfft              = dtin%mgfft
 dtout%mgfftdg            = dtin%mgfftdg
 dtout%mkmem              = dtin%mkmem
 dtout%mkqmem             = dtin%mkqmem
 dtout%mk1mem             = dtin%mk1mem
 dtout%mpw                = dtin%mpw
 dtout%mqgrid             = dtin%mqgrid
 dtout%mqgriddg           = dtin%mqgriddg
 dtout%natom              = dtin%natom
 dtout%natrd              = dtin%natrd
 dtout%natsph             = dtin%natsph
 dtout%natsph_extra       = dtin%natsph_extra
 dtout%natpawu            = dtin%natpawu
 dtout%natvshift          = dtin%natvshift
 dtout%nbdblock           = dtin%nbdblock
 dtout%nbdbuf             = dtin%nbdbuf
 dtout%nbandhf            = dtin%nbandhf
 dtout%nberry             = dtin%nberry
 dtout%nc_xccc_gspace     = dtin%nc_xccc_gspace
 dtout%nbandkss           = dtin%nbandkss
 dtout%nconeq             = dtin%nconeq
 dtout%nctime             = dtin%nctime
 dtout%ndtset             = dtin%ndtset
 dtout%ndynimage          = dtin%ndynimage
 dtout%neb_algo           = dtin%neb_algo
 dtout%nfft               = dtin%nfft
 dtout%nfftdg             = dtin%nfftdg
 dtout%nfreqim            = dtin%nfreqim
 dtout%nfreqre            = dtin%nfreqre
 dtout%nfreqsp            = dtin%nfreqsp
 dtout%nimage             = dtin%nimage
 dtout%nkpt               = dtin%nkpt
 dtout%nkpthf             = dtin%nkpthf
 dtout%nkptgw             = dtin%nkptgw
 dtout%nonlinear_info     = dtin%nonlinear_info
 dtout%nline              = dtin%nline
 dtout%nnsclo             = dtin%nnsclo
 dtout%nnsclohf           = dtin%nnsclohf
 dtout%nomegasf           = dtin%nomegasf
 dtout%nomegasi           = dtin%nomegasi
 dtout%nomegasrd          = dtin%nomegasrd
 dtout%npband             = dtin%npband
 dtout%npfft              = dtin%npfft
 dtout%nphf               = dtin%nphf
 dtout%npimage            = dtin%npimage
 dtout%npkpt              = dtin%npkpt
 dtout%nppert             = dtin%nppert
 dtout%npspinor           = dtin%npspinor
 dtout%npsp               = dtin%npsp
 dtout%npspalch           = dtin%npspalch
 dtout%npulayit           = dtin%npulayit
 dtout%npvel              = dtin%npvel
 dtout%npweps             = dtin%npweps
 dtout%npwkss             = dtin%npwkss
 dtout%npwsigx            = dtin%npwsigx
 dtout%npwwfn             = dtin%npwwfn
 dtout%np_slk             = dtin%np_slk
 dtout%nqpt               = dtin%nqpt
 dtout%nqptdm             = dtin%nqptdm
 dtout%nscforder          = dtin%nscforder
 dtout%nshiftk            = dtin%nshiftk
 dtout%nshiftk_orig       = dtin%nshiftk_orig
 dtout%nspden             = dtin%nspden
 dtout%nspinor            = dtin%nspinor
 dtout%nsppol             = dtin%nsppol
 dtout%nstep              = dtin%nstep
 dtout%nsym               = dtin%nsym
 dtout%ntime              = dtin%ntime
 dtout%ntimimage          = dtin%ntimimage
 dtout%ntypalch           = dtin%ntypalch
 dtout%ntypat             = dtin%ntypat
 dtout%ntyppure           = dtin%ntyppure
 dtout%nwfshist           = dtin%nwfshist
 dtout%nzchempot          = dtin%nzchempot
 dtout%occopt             = dtin%occopt
 dtout%optcell            = dtin%optcell
 dtout%optdriver          = dtin%optdriver
 dtout%optforces          = dtin%optforces
 dtout%optnlxccc          = dtin%optnlxccc
 dtout%optstress          = dtin%optstress
 dtout%orbmag             = dtin%orbmag
 dtout%ortalg             = dtin%ortalg
 dtout%paral_atom         = dtin%paral_atom
 dtout%paral_kgb          = dtin%paral_kgb
 dtout%paral_rf           = dtin%paral_rf
 dtout%pawcpxocc          = dtin%pawcpxocc
 dtout%pawcross           = dtin%pawcross
 dtout%pawlcutd           = dtin%pawlcutd
 dtout%pawlmix            = dtin%pawlmix
 dtout%pawmixdg           = dtin%pawmixdg
 dtout%pawnhatxc          = dtin%pawnhatxc
 dtout%pawnphi            = dtin%pawnphi
 dtout%pawntheta          = dtin%pawntheta
 dtout%pawnzlm            = dtin%pawnzlm
 dtout%pawoptmix          = dtin%pawoptmix
 dtout%pawoptosc          = dtin%pawoptosc
 dtout%pawprtdos          = dtin%pawprtdos
 dtout%pawprtvol          = dtin%pawprtvol
 dtout%pawprtwf           = dtin%pawprtwf
 dtout%pawprt_k           = dtin%pawprt_k
 dtout%pawprt_b           = dtin%pawprt_b
 dtout%pawspnorb          = dtin%pawspnorb
 dtout%pawstgylm          = dtin%pawstgylm
 dtout%pawsushat          = dtin%pawsushat
 dtout%pawusecp           = dtin%pawusecp
 dtout%pawujat            = dtin%pawujat
 dtout%macro_uj           = dtin%macro_uj
 dtout%pawujrad           = dtin%pawujrad
 dtout%pawujv             = dtin%pawujv
 dtout%pawxcdev           = dtin%pawxcdev
 dtout%pimd_constraint    = dtin%pimd_constraint
 dtout%pitransform        = dtin%pitransform
 dtout%plowan_compute     = dtin%plowan_compute
 dtout%plowan_bandi       = dtin%plowan_bandi
 dtout%plowan_bandf       = dtin%plowan_bandf
 dtout%plowan_natom       = dtin%plowan_natom
 dtout%plowan_nt          = dtin%plowan_nt
 dtout%plowan_realspace   = dtin%plowan_realspace
 dtout%posdoppler         = dtin%posdoppler
 dtout%positron           = dtin%positron
 dtout%posnstep           = dtin%posnstep
 dtout%ppmodel            = dtin%ppmodel
 dtout%prepalw            = dtin%prepalw
 dtout%prepanl            = dtin%prepanl
 dtout%prepgkk            = dtin%prepgkk
 dtout%prtbbb             = dtin%prtbbb
 dtout%prtbltztrp         = dtin%prtbltztrp
 dtout%prtcif             = dtin%prtcif
 dtout%prtden             = dtin%prtden
 dtout%prtdensph          = dtin%prtdensph
 dtout%prtdipole          = dtin%prtdipole
 dtout%prtdos             = dtin%prtdos
 dtout%prtdosm            = dtin%prtdosm
 dtout%prtebands          = dtin%prtebands    ! TODO prteig could be replaced by prtebands...
 dtout%prtefg             = dtin%prtefg
 dtout%prtefmas           = dtin%prtefmas
 dtout%prteig             = dtin%prteig
 dtout%prtelf             = dtin%prtelf
 dtout%prteliash          = dtin%prteliash
 dtout%prtfc              = dtin%prtfc
 dtout%prtfull1wf         = dtin%prtfull1wf
 dtout%prtfsurf           = dtin%prtfsurf
 dtout%prtgsr             = dtin%prtgsr
 dtout%prtgden            = dtin%prtgden
 dtout%prtgeo             = dtin%prtgeo
 dtout%prtgkk             = dtin%prtgkk
 dtout%prtkden            = dtin%prtkden
 dtout%prtkpt             = dtin%prtkpt
 dtout%prtlden            = dtin%prtlden
 dtout%prtnabla           = dtin%prtnabla
 dtout%prtnest            = dtin%prtnest
 dtout%prtphbands         = dtin%prtphbands
 dtout%prtphdos           = dtin%prtphdos
 dtout%prtphsurf          = dtin%prtphsurf
 dtout%prtposcar          = dtin%prtposcar
 dtout%prtprocar          = dtin%prtprocar
 dtout%prtpot             = dtin%prtpot
 dtout%prtpsps            = dtin%prtpsps
 dtout%prtspcur           = dtin%prtspcur
 dtout%prtsuscep          = dtin%prtsuscep
 dtout%prtstm             = dtin%prtstm
 dtout%prtvclmb           = dtin%prtvclmb
 dtout%prtvdw             = dtin%prtvdw
 dtout%prtvha             = dtin%prtvha
 dtout%prtvhxc            = dtin%prtvhxc
 dtout%prtkbff            = dtin%prtkbff
 dtout%prtvol             = dtin%prtvol
 dtout%prtvolimg          = dtin%prtvolimg
 dtout%prtvpsp            = dtin%prtvpsp
 dtout%prtvxc             = dtin%prtvxc
 dtout%prtwant            = dtin%prtwant
 dtout%prtwf              = dtin%prtwf
 dtout%prtwf_full         = dtin%prtwf_full
 dtout%prtxml             = dtin%prtxml
 dtout%prt1dm             = dtin%prt1dm
 dtout%ptgroupma          = dtin%ptgroupma
 dtout%qptopt             = dtin%qptopt
 dtout%random_atpos       = dtin%random_atpos
 dtout%recgratio          = dtin%recgratio
 dtout%recnpath           = dtin%recnpath
 dtout%recnrec            = dtin%recnrec
 dtout%recptrott          = dtin%recptrott
 dtout%rectesteg          = dtin%rectesteg
 dtout%rcut               = dtin%rcut
 dtout%restartxf          = dtin%restartxf
 dtout%rfasr              = dtin%rfasr
 dtout%rfddk              = dtin%rfddk
 dtout%rfelfd             = dtin%rfelfd
 dtout%rfmagn             = dtin%rfmagn
 dtout%rfmeth             = dtin%rfmeth
 dtout%rfphon             = dtin%rfphon
 dtout%rfstrs             = dtin%rfstrs
 dtout%rfuser             = dtin%rfuser
 dtout%rf2_dkdk           = dtin%rf2_dkdk
 dtout%rf2_dkde           = dtin%rf2_dkde
 dtout%rhoqpmix           = dtin%rhoqpmix
 dtout%signperm           = dtin%signperm
 dtout%slabwsrad          = dtin%slabwsrad
 dtout%slabzbeg           = dtin%slabzbeg
 dtout%slabzend           = dtin%slabzend
 dtout%slk_rankpp         = dtin%slk_rankpp
 dtout%smdelta            = dtin%smdelta
 dtout%spgaxor            = dtin%spgaxor
 dtout%spgorig            = dtin%spgorig
 dtout%spgroup            = dtin%spgroup
 dtout%spmeth             = dtin%spmeth
 dtout%string_algo        = dtin%string_algo
 dtout%symchi             = dtin%symchi
 dtout%symmorphi          = dtin%symmorphi
 dtout%symsigma           = dtin%symsigma
 dtout%td_mexcit          = dtin%td_mexcit
 dtout%tfkinfunc          = dtin%tfkinfunc
 dtout%tim1rev            = dtin%tim1rev
 dtout%timopt             = dtin%timopt
 dtout%use_gemm_nonlop    = dtin%use_gemm_nonlop
 dtout%use_gpu_cuda       = dtin%use_gpu_cuda
 dtout%use_yaml           = dtin%use_yaml   ! This variable activates the Yaml output for testing purposes
                                            ! It will be removed when Yaml output enters production.
 dtout%use_slk            = dtin%use_slk
 dtout%usedmatpu          = dtin%usedmatpu
 dtout%usedmft            = dtin%usedmft
 dtout%useexexch          = dtin%useexexch
 dtout%usefock            = dtin%usefock
 dtout%usekden            = dtin%usekden
 dtout%use_nonscf_gkk     = dtin%use_nonscf_gkk
 dtout%usepaw             = dtin%usepaw
 dtout%usepawu            = dtin%usepawu
 dtout%usepead            = dtin%usepead
 dtout%usepotzero         = dtin%usepotzero
 dtout%userec             = dtin%userec
 dtout%useria             = dtin%useria
 dtout%userib             = dtin%userib
 dtout%useric             = dtin%useric
 dtout%userid             = dtin%userid
 dtout%userie             = dtin%userie
 dtout%usewvl             = dtin%usewvl
 dtout%usexcnhat_orig     = dtin%usexcnhat_orig
 dtout%useylm             = dtin%useylm
 dtout%vacnum             = dtin%vacnum
 dtout%vdw_df_acutmin     = dtin%vdw_df_acutmin
 dtout%vdw_df_aratio      = dtin%vdw_df_aratio
 dtout%vdw_df_damax       = dtin%vdw_df_damax
 dtout%vdw_df_damin       = dtin%vdw_df_damin
 dtout%vdw_df_dcut        = dtin%vdw_df_dcut
 dtout%vdw_df_dratio      = dtin%vdw_df_dratio
 dtout%vdw_df_dsoft       = dtin%vdw_df_dsoft
 dtout%vdw_df_gcut        = dtin%vdw_df_gcut
 dtout%vdw_df_ndpts       = dtin%vdw_df_ndpts
 dtout%vdw_df_ngpts       = dtin%vdw_df_ngpts
 dtout%vdw_df_nqpts       = dtin%vdw_df_nqpts
 dtout%vdw_df_nrpts       = dtin%vdw_df_nrpts
 dtout%vdw_df_nsmooth     = dtin%vdw_df_nsmooth
 dtout%vdw_df_phisoft     = dtin%vdw_df_phisoft
 dtout%vdw_df_qcut        = dtin%vdw_df_qcut
 dtout%vdw_df_qratio      = dtin%vdw_df_qratio
 dtout%vdw_df_rcut        = dtin%vdw_df_rcut
 dtout%vdw_df_rsoft       = dtin%vdw_df_rsoft
 dtout%vdw_df_tolerance   = dtin%vdw_df_tolerance
 dtout%vdw_df_threshold   = dtin%vdw_df_threshold
 dtout%vdw_df_tweaks      = dtin%vdw_df_tweaks
 dtout%vdw_df_zab         = dtin%vdw_df_zab
 dtout%vdw_nfrag          = dtin%vdw_nfrag
 dtout%vdw_xc             = dtin%vdw_xc
 dtout%wfoptalg           = dtin%wfoptalg
 dtout%wvl_bigdft_comp    = dtin%wvl_bigdft_comp
 dtout%w90iniprj          = dtin%w90iniprj
 dtout%w90prtunk          = dtin%w90prtunk
 dtout%xclevel            = dtin%xclevel
 dtout%xc_denpos          = dtin%xc_denpos

!Copy allocated integer arrays from dtin to dtout
 dtout%bdberry(:)         = dtin%bdberry(:)
 dtout%cd_subset_freq(:)  = dtin%cd_subset_freq(:)
 dtout%d3e_pert1_atpol(:) = dtin%d3e_pert1_atpol(:)
 dtout%d3e_pert1_dir(:)   = dtin%d3e_pert1_dir(:)
 dtout%d3e_pert2_atpol(:) = dtin%d3e_pert2_atpol(:)
 dtout%d3e_pert2_dir(:)   = dtin%d3e_pert2_dir(:)
 dtout%d3e_pert3_atpol(:) = dtin%d3e_pert3_atpol(:)
 dtout%d3e_pert3_dir(:)   = dtin%d3e_pert3_dir(:)
 dtout%ga_rules(:)        = dtin%ga_rules(:)
 dtout%gpu_devices(:)     = dtin%gpu_devices(:)
 dtout%jfielddir(:)       = dtin%jfielddir(:)
 dtout%kptrlatt(:,:)      = dtin%kptrlatt(:,:)
 dtout%kptrlatt_orig      = dtin%kptrlatt_orig
 dtout%qptrlatt(:,:)      = dtin%qptrlatt(:,:)
 dtout%ngfft(:)           = dtin%ngfft(:)
 dtout%ngfftdg(:)         = dtin%ngfftdg(:)
 dtout%nloalg(:)          = dtin%nloalg(:)
 dtout%ngkpt(:)           = dtin%ngkpt(:)
 dtout%qprtrb(:)          = dtin%qprtrb(:)
 dtout%rfatpol(:)         = dtin%rfatpol(:)
 dtout%rfdir(:)           = dtin%rfdir(:)
 dtout%rf2_pert1_dir(:)   = dtin%rf2_pert1_dir(:)
 dtout%rf2_pert2_dir(:)   = dtin%rf2_pert2_dir(:)
 dtout%supercell_latt(:)= dtin%supercell_latt(:)
 dtout%ucrpa_bands(:)     = dtin%ucrpa_bands(:)
 dtout%vdw_supercell(:)   = dtin%vdw_supercell(:)
 dtout%vdw_typfrag(:)     = dtin%vdw_typfrag(:)
 dtout%wvl_ngauss(:)      = dtin%wvl_ngauss(:)

!Copy reals from dtin to dtout
 dtout%adpimd_gamma       = dtin%adpimd_gamma
 dtout%boxcutmin          = dtin%boxcutmin
 dtout%bxctmindg          = dtin%bxctmindg
 dtout%cd_halfway_freq    = dtin%cd_halfway_freq
 dtout%cd_max_freq        = dtin%cd_max_freq
 dtout%charge             = dtin%charge
 dtout%cpus               = dtin%cpus
 dtout%ddamp              = dtin%ddamp
 dtout%diecut             = dtin%diecut
 dtout%diegap             = dtin%diegap
 dtout%dielam             = dtin%dielam
 dtout%dielng             = dtin%dielng
 dtout%diemac             = dtin%diemac
 dtout%diemix             = dtin%diemix
 dtout%diemixmag          = dtin%diemixmag
 dtout%dilatmx            = dtin%dilatmx
 dtout%dosdeltae          = dtin%dosdeltae
 dtout%dtion              = dtin%dtion
 dtout%ecut               = dtin%ecut
 dtout%ecuteps            = dtin%ecuteps
 dtout%ecutsigx           = dtin%ecutsigx
 dtout%ecutsm             = dtin%ecutsm
 dtout%ecutwfn            = dtin%ecutwfn
 dtout%effmass_free       = dtin%effmass_free
 dtout%efmas_deg_tol      = dtin%efmas_deg_tol
 dtout%elph2_imagden      = dtin%elph2_imagden
 dtout%eshift             = dtin%eshift
 dtout%esmear             = dtin%esmear
 dtout%exchmix            = dtin%exchmix
 dtout%fband              = dtin%fband
 dtout%focktoldfe         = dtin%focktoldfe
 dtout%friction           = dtin%friction
 dtout%fxcartfactor       = dtin%fxcartfactor
 dtout%ga_opt_percent     = dtin%ga_opt_percent
 dtout%gwls_model_parameter = dtin%gwls_model_parameter
 dtout%kptnrm             = dtin%kptnrm
 dtout%kptrlen            = dtin%kptrlen
 dtout%maxestep           = dtin%maxestep
 dtout%bmass              = dtin%bmass
 dtout%magcon_lambda      = dtin%magcon_lambda
 dtout%mdwall             = dtin%mdwall
 dtout%mep_mxstep         = dtin%mep_mxstep
 dtout%nelect             = dtin%nelect
 dtout%nnos               = dtin%nnos
 dtout%noseinert          = dtin%noseinert
 dtout%omegasimax         = dtin%omegasimax
 dtout%omegasrdmax        = dtin%omegasrdmax
 dtout%pawecutdg          = dtin%pawecutdg
 dtout%pawovlp            = dtin%pawovlp
 dtout%posocc             = dtin%posocc
 dtout%postoldfe          = dtin%postoldfe
 dtout%postoldff          = dtin%postoldff
 dtout%ppmfrq             = dtin%ppmfrq
 dtout%pw_unbal_thresh    = dtin%pw_unbal_thresh
 dtout%ratsm              = dtin%ratsm
 dtout%ratsph_extra       = dtin%ratsph_extra
 dtout%recrcut            = dtin%recrcut
 dtout%recefermi          = dtin%recefermi
 dtout%rectolden          = dtin%rectolden
 dtout%dfpt_sciss         = dtin%dfpt_sciss
 dtout%mbpt_sciss         = dtin%mbpt_sciss
 dtout%spinmagntarget     = dtin%spinmagntarget
 dtout%spbroad            = dtin%spbroad
 dtout%spnorbscl          = dtin%spnorbscl
 dtout%stmbias            = dtin%stmbias
 dtout%strfact            = dtin%strfact
 dtout%strprecon          = dtin%strprecon
 dtout%tfw_toldfe         = dtin%tfw_toldfe
 dtout%tl_radius          = dtin%tl_radius
 dtout%tl_nprccg          = dtin%tl_nprccg
 dtout%td_maxene          = dtin%td_maxene
 dtout%toldfe             = dtin%toldfe
 dtout%tolmxde            = dtin%tolmxde
 dtout%toldff             = dtin%toldff
 dtout%tolimg             = dtin%tolimg
 dtout%tolmxf             = dtin%tolmxf
 dtout%tolrde             = dtin%tolrde
 dtout%tolrff             = dtin%tolrff
 dtout%tolsym             = dtin%tolsym
 dtout%tolvrs             = dtin%tolvrs
 dtout%tolwfr             = dtin%tolwfr
 dtout%tphysel            = dtin%tphysel
 dtout%tsmear             = dtin%tsmear
 dtout%ucrpa              = dtin%ucrpa
 dtout%userra             = dtin%userra
 dtout%userrb             = dtin%userrb
 dtout%userrc             = dtin%userrc
 dtout%userrd             = dtin%userrd
 dtout%userre             = dtin%userre
 dtout%vacwidth           = dtin%vacwidth
 dtout%vdw_tol            = dtin%vdw_tol
 dtout%vdw_tol_3bt        = dtin%vdw_tol_3bt
 dtout%vis                = dtin%vis
 dtout%wfmix              = dtin%wfmix
 dtout%wfk_task           = dtin%wfk_task
 dtout%wtq                = dtin%wtq
 dtout%wvl_hgrid          = dtin%wvl_hgrid
 dtout%wvl_crmult         = dtin%wvl_crmult
 dtout%wvl_frmult         = dtin%wvl_frmult
 dtout%wvl_nprccg         = dtin%wvl_nprccg
 dtout%xc_tb09_c          = dtin%xc_tb09_c
 dtout%zcut               = dtin%zcut

!Copy allocated real arrays from dtin to dtout
 dtout%boxcenter(:)       = dtin%boxcenter(:)
 dtout%bfield(:)          = dtin%bfield(:)
 dtout%dfield(:)          = dtin%dfield(:)
 dtout%efield(:)          = dtin%efield(:)
 dtout%genafm(:)          = dtin%genafm(:)
 dtout%goprecprm(:)       = dtin%goprecprm(:)
 dtout%mdtemp(:)          = dtin%mdtemp(:)
 dtout%neb_spring(:)      = dtin%neb_spring(:)
 dtout%polcen(:)          = dtin%polcen(:)
 dtout%qptn(:)            = dtin%qptn(:)
 dtout%pvelmax(:)         = dtin%pvelmax(:)
 dtout%red_efield(:)      = dtin%red_efield(:)
 dtout%red_dfield(:)      = dtin%red_dfield(:)
 dtout%red_efieldbar(:)   = dtin%red_efieldbar(:)
 dtout%shiftk_orig        = dtin%shiftk_orig
 dtout%strtarget(:)       = dtin%strtarget(:)
 dtout%ucrpa_window(:)    = dtin%ucrpa_window(:)
 dtout%vcutgeo(:)         = dtin%vcutgeo(:)
 dtout%vprtrb(:)          = dtin%vprtrb(:)
 dtout%zeemanfield(:)     = dtin%zeemanfield(:)

!Use alloc_copy to allocate and copy the allocatable arrays

!integer allocatables
 call alloc_copy(dtin%algalch, dtout%algalch)
 call alloc_copy(dtin%bdgw, dtout%bdgw)
 call alloc_copy(dtin%bs_loband, dtout%bs_loband)
 call alloc_copy(dtin%constraint_kind, dtout%constraint_kind)
 call alloc_copy(dtin%dynimage, dtout%dynimage)
 call alloc_copy(dtin%efmas_bands, dtout%efmas_bands)
 call alloc_copy(dtin%iatfix, dtout%iatfix)
 call alloc_copy(dtin%iatsph, dtout%iatsph)
 call alloc_copy(dtin%istwfk, dtout%istwfk)
 call alloc_copy(dtin%kberry, dtout%kberry)
 call alloc_copy(dtin%lexexch, dtout%lexexch)
 call alloc_copy(dtin%ldaminushalf, dtout%ldaminushalf)
 call alloc_copy(dtin%lpawu, dtout%lpawu)
 call alloc_copy(dtin%nband, dtout%nband)
 call alloc_copy(dtin%plowan_iatom, dtout%plowan_iatom)
 call alloc_copy(dtin%plowan_it, dtout%plowan_it)
 call alloc_copy(dtin%plowan_nbl, dtout%plowan_nbl)
 call alloc_copy(dtin%plowan_lcalc, dtout%plowan_lcalc)
 call alloc_copy(dtin%plowan_projcalc, dtout%plowan_projcalc)
 call alloc_copy(dtin%prtatlist, dtout%prtatlist)
 call alloc_copy(dtin%so_psp, dtout%so_psp)
 call alloc_copy(dtin%symafm, dtout%symafm)
 call alloc_copy(dtin%symrel, dtout%symrel)
 call alloc_copy(dtin%typat, dtout%typat)

!Allocate and copy real allocatable
 call alloc_copy(dtin%acell_orig, dtout%acell_orig)
 call alloc_copy(dtin%amu_orig, dtout%amu_orig)
 call alloc_copy(dtin%atvshift, dtout%atvshift)
 call alloc_copy(dtin%chrgat, dtout%chrgat)
 call alloc_copy(dtin%cd_imfrqs, dtout%cd_imfrqs)
 call alloc_copy(dtin%chempot, dtout%chempot)
 call alloc_copy(dtin%corecs, dtout%corecs)
 call alloc_copy(dtin%densty, dtout%densty)
 call alloc_copy(dtin%dmatpawu, dtout%dmatpawu)
 call alloc_copy(dtin%efmas_dirs, dtout%efmas_dirs)
 call alloc_copy(dtin%f4of2_sla, dtout%f4of2_sla)
 call alloc_copy(dtin%f6of2_sla, dtout%f6of2_sla)
 call alloc_copy(dtin%gw_qlwl, dtout%gw_qlwl)
 call alloc_copy(dtin%gw_freqsp, dtout%gw_freqsp)
 call alloc_copy(dtin%gwls_list_proj_freq, dtout%gwls_list_proj_freq)
 call alloc_copy(dtin%jpawu, dtout%jpawu)
 call alloc_copy(dtin%kpt, dtout%kpt)
 call alloc_copy(dtin%kptgw, dtout%kptgw)
 call alloc_copy(dtin%kptns, dtout%kptns)
 call alloc_copy(dtin%kptns_hf, dtout%kptns_hf)
 call alloc_copy(dtin%mixalch_orig, dtout%mixalch_orig)
 call alloc_copy(dtin%mixesimgf, dtout%mixesimgf)
 call alloc_copy(dtin%nucdipmom, dtout%nucdipmom)
 call alloc_copy(dtin%occ_orig, dtout%occ_orig)
 call alloc_copy(dtin%pimass, dtout%pimass)
 call alloc_copy(dtin%ptcharge, dtout%ptcharge)
 call alloc_copy(dtin%qmass, dtout%qmass)
 call alloc_copy(dtin%qptdm, dtout%qptdm)
 call alloc_copy(dtin%quadmom, dtout%quadmom)
 call alloc_copy(dtin%ratsph, dtout%ratsph)
 call alloc_copy(dtin%rprim_orig, dtout%rprim_orig)
 call alloc_copy(dtin%rprimd_orig, dtout%rprimd_orig)
 call alloc_copy(dtin%shiftk, dtout%shiftk)
 call alloc_copy(dtin%spinat, dtout%spinat)
 call alloc_copy(dtin%tnons, dtout%tnons)
 call alloc_copy(dtin%upawu, dtout%upawu)
 call alloc_copy(dtin%vel_orig, dtout%vel_orig)
 call alloc_copy(dtin%vel_cell_orig, dtout%vel_cell_orig)
 call alloc_copy(dtin%wtatcon, dtout%wtatcon)
 call alloc_copy(dtin%wtk, dtout%wtk)
 call alloc_copy(dtin%xred_orig, dtout%xred_orig)
 call alloc_copy(dtin%xredsph_extra, dtout%xredsph_extra)
 call alloc_copy(dtin%ziontypat, dtout%ziontypat)
 call alloc_copy(dtin%znucl, dtout%znucl)

 dtout%ndivsm = dtin%ndivsm
 dtout%nkpath = dtin%nkpath
 dtout%einterp = dtin%einterp
 call alloc_copy(dtin%kptbounds, dtout%kptbounds)
 dtout%tmesh = dtin%tmesh
 dtout%getkerange_filepath = dtin%getkerange_filepath

 DBG_EXIT("COLL")

end function dtset_copy
!!***

!----------------------------------------------------------------------

!!****f* m_dtset/dtset_free
!! NAME
!! dtset_free
!!
!! FUNCTION
!! Free a dataset after use.
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=free all allocated allocatable.
!!
!! PARENTS
!!      m_xchybrid
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine dtset_free(dtset)

!Arguments ------------------------------------
!scalars
 class(dataset_type),intent(inout) :: dtset

! *************************************************************************

!please, use the same order as the one used in the declaration of the type (see defs_abitypes).

 !@dataset_type
!integer allocatable
 ABI_SFREE(dtset%algalch)
 ABI_SFREE(dtset%bdgw)
 ABI_SFREE(dtset%bs_loband)
 ABI_SFREE(dtset%constraint_kind)
 ABI_SFREE(dtset%dynimage)
 ABI_SFREE(dtset%efmas_bands)
 ABI_SFREE(dtset%iatfix)
 ABI_SFREE(dtset%iatsph)
 ABI_SFREE(dtset%istwfk)
 ABI_SFREE(dtset%kberry)
 ABI_SFREE(dtset%lexexch)
 ABI_SFREE(dtset%ldaminushalf)
 ABI_SFREE(dtset%lpawu)
 ABI_SFREE(dtset%nband)
 ABI_SFREE(dtset%ph_qpath)
 ABI_SFREE(dtset%ph_qshift)
 ABI_SFREE(dtset%plowan_iatom)
 ABI_SFREE(dtset%plowan_it)
 ABI_SFREE(dtset%plowan_lcalc)
 ABI_SFREE(dtset%plowan_nbl)
 ABI_SFREE(dtset%plowan_projcalc)
 ABI_SFREE(dtset%prtatlist)
 ABI_SFREE(dtset%so_psp)
 ABI_SFREE(dtset%symafm)
 ABI_SFREE(dtset%symrel)
 ABI_SFREE(dtset%typat)

!real allocatable
 ABI_SFREE(dtset%acell_orig)
 ABI_SFREE(dtset%amu_orig)
 ABI_SFREE(dtset%atvshift)
 ABI_SFREE(dtset%cd_imfrqs)
 ABI_SFREE(dtset%chrgat)
 ABI_SFREE(dtset%chempot)
 ABI_SFREE(dtset%corecs)
 ABI_SFREE(dtset%densty)
 ABI_SFREE(dtset%dmatpawu)
 ABI_SFREE(dtset%efmas_dirs)
 ABI_SFREE(dtset%gw_qlwl)
 ABI_SFREE(dtset%gw_freqsp)
 ABI_SFREE(dtset%gwls_list_proj_freq)
 ABI_SFREE(dtset%f4of2_sla)
 ABI_SFREE(dtset%f6of2_sla)
 ABI_SFREE(dtset%jpawu)
 ABI_SFREE(dtset%kpt)
 ABI_SFREE(dtset%kptbounds)
 ABI_SFREE(dtset%kptgw)
 ABI_SFREE(dtset%kptns)
 ABI_SFREE(dtset%kptns_hf)
 ABI_SFREE(dtset%mixalch_orig)
 ABI_SFREE(dtset%mixesimgf)
 ABI_SFREE(dtset%nucdipmom)
 ABI_SFREE(dtset%occ_orig)
 ABI_SFREE(dtset%pimass)
 ABI_SFREE(dtset%ptcharge)
 ABI_SFREE(dtset%qmass)
 ABI_SFREE(dtset%qptdm)
 ABI_SFREE(dtset%quadmom)
 ABI_SFREE(dtset%ratsph)
 ABI_SFREE(dtset%rprim_orig)
 ABI_SFREE(dtset%rprimd_orig)
 ABI_SFREE(dtset%shiftk)
 ABI_SFREE(dtset%spinat)
 ABI_SFREE(dtset%tnons)
 ABI_SFREE(dtset%sigma_shiftk)
 ABI_SFREE(dtset%upawu)
 ABI_SFREE(dtset%vel_orig)
 ABI_SFREE(dtset%vel_cell_orig)
 ABI_SFREE(dtset%wtatcon)
 ABI_SFREE(dtset%wtk)
 ABI_SFREE(dtset%xred_orig)
 ABI_SFREE(dtset%xredsph_extra)
 ABI_SFREE(dtset%ziontypat)
 ABI_SFREE(dtset%znucl)

end subroutine dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_dtset/dtset_free_nkpt_arrays
!! NAME
!! dtset_free_nkpt_arrays
!!
!! FUNCTION
!!  Free arrays that depend on input nkpt (used in EPH code, because EPH has its own
!!  treatment of BZ sampling and we don't want to waste memory with large and useless arrays
!!  especially if very dense k-meshes are used.
!!
!! PARENTS
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine dtset_free_nkpt_arrays(dtset)

!Arguments ------------------------------------
!scalars
 class(dataset_type),intent(inout) :: dtset

! *************************************************************************

 ABI_SFREE(dtset%istwfk)
 !ABI_SFREE(dtset%nband)
 ABI_SFREE(dtset%kpt)
 ABI_SFREE(dtset%kptns)
 ABI_SFREE(dtset%occ_orig)
 ABI_SFREE(dtset%wtk)
 ! Free HF k-points as well.
 ABI_SFREE(dtset%kptns_hf)

end subroutine dtset_free_nkpt_arrays
!!***

!!****f* m_dtset/find_getdtset
!! NAME
!! find_getdtset
!!
!! FUNCTION
!! Find the number of the dataset (iget) for a given value of a "get" variable (getvalue)
!! of name getname, given the number of the current dataset (idtset).
!! Also find the coefficients of mixing of the images of the old dataset, to initialize the new dataset images
!! (use a linear interpolation)
!!
!! INPUTS
!! dtsets(0:ndtset_alloc)=<type datasets_type>contains all input variables
!! getvalue=value of the get variable
!! getname=name of the get variable
!! idtset=number of the current dataset
!! mxnimage=dimension of miximage
!! ndtset_alloc=dimension of dtsets
!!
!! OUTPUT
!! iget=number of the dataset from which the value must be get, 0 if the data should not be got from another dataset
!! miximage(mxnimage,mxnimage)=coefficients of mixing of the images of the old dataset, to initialize the new dataset images
!!
!! PARENTS
!!      m_driver
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine find_getdtset(dtsets,getvalue,getname,idtset,iget,miximage,mxnimage,ndtset_alloc)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: getvalue,idtset,mxnimage,ndtset_alloc
 integer, intent(out) :: iget
 real(dp), intent(out) :: miximage(mxnimage,mxnimage)
 character(len=*),intent(in) :: getname
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables-------------------------------
 integer :: iimage
 real(dp) :: newimage_get,ratio
 character(len=500) :: msg

! *************************************************************************

 iget=0
 if(getvalue>0 .or. (getvalue<0 .and. idtset+getvalue>0) )then
!  In case getvalue is a negative number (so must add to idtset)
   if(getvalue<0 .and. idtset+getvalue>0) iget=idtset+getvalue
   if(getvalue>0)then
     do iget=1,idtset
       if( dtsets(iget)%jdtset==getvalue )exit
     end do
     if(iget==idtset)then
!      The index of the dataset, from which the data ought to be taken,
!      does not correspond to a previous dataset.
       write(msg, '(a,i0,4a,i3,7a)' )&
        'The component number ',idtset,' of the input variable ',trim(getname),',',' equal to ',getvalue,',',ch10,&
        'does not correspond to an existing index.',ch10,&
        'Action: correct ',trim(getname),' or jdtset in your input file.'
       MSG_ERROR(msg)
     end if
   end if
   write(msg, '(3a,i3,2a)' )&
&   ' find_getdtset : ',trim(getname),'/=0, take data from output of dataset with index',dtsets(iget)%jdtset,'.',ch10
   call wrtout([std_out, ab_out], msg)
 end if

!For the time being, uses a simple interpolation when the images do not match. If only one image, take the first get image.
 miximage(:,:)=zero
 if(dtsets(idtset)%nimage==1)then
   miximage(1,1)=one
 else
   do iimage=1,dtsets(idtset)%nimage
     ratio=(iimage-one)/real(dtsets(idtset)%nimage-one)
     newimage_get=one+ratio*(dtsets(iget)%nimage-one)
     if(abs(newimage_get-nint(newimage_get))<tol8)then
       miximage(iimage,nint(newimage_get))=one
     else
       miximage(iimage,floor(newimage_get))=one-(newimage_get-floor(newimage_get))
       miximage(iimage,ceiling(newimage_get))=one-miximage(iimage,floor(newimage_get))
     end if
   end do
 end if

end subroutine find_getdtset
!!***

!!****f* m_dtset/dtset_get_npert_rbz
!! NAME
!! dtset_get_npert_rbz
!!
!! FUNCTION
!! Get the number of effective pertubation done in looper3, nkpt_rbz, nband_rbz
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!!  npert=number of effective pertubation done in looper3
!!  nkpt_rbz= nkpt in the reduced brillouin zone
!!  nband_rbz= nband in the reduced brillouin zone
!!
!! PARENTS
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine dtset_get_npert_rbz(dtset, nband_rbz, nkpt_rbz, npert)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: npert
!arrays
 integer,pointer :: nkpt_rbz(:)
 real(dp),pointer :: nband_rbz(:,:)
 class(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
 integer :: icase,idir,ikpt,ikpt1,ipert,isppol,isym,maxidir,mpert,nband_k,nsym1,timrev,timrev_pert
 integer :: to_compute_this_pert
 real(dp) :: tolsym8,ucvol
 character(len=500) :: msg
!arrays
 integer :: rfdir(9),rf2dir(9),rf2_dir1(3),rf2_dir2(3)
 integer,allocatable :: indkpt1(:,:),indsym(:,:,:),pertsy(:,:),rfpert(:),symq(:,:,:),symrec(:,:,:)
 integer, allocatable :: pert_tmp(:,:), pert_calc(:,:)
 integer,allocatable :: symaf1(:),symrc1(:,:,:),symrl1(:,:,:),symrl1_tmp(:,:,:), bz2ibz_smap(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: tnons1_tmp(:,:),wtk_folded(:)

! *************************************************************************

!Define the set of admitted perturbations
 mpert=dtset%natom+6
 if(dtset%natom+10/=0.or.dtset%natom+11/=0) mpert=dtset%natom+11

 ABI_MALLOC(symrec,(3,3,dtset%nsym))
!Get the symmetry matrices in terms of reciprocal basis
 do isym=1,dtset%nsym
   call mati3inv(dtset%symrel(:,:,isym),symrec(:,:,isym))
 end do

 ABI_MALLOC(indsym,(4,dtset%nsym,dtset%natom))
!Obtain a list of rotated atom labels:
 tolsym8=tol8
 call symatm(indsym,dtset%natom,dtset%nsym,symrec,dtset%tnons,tolsym8,dtset%typat,dtset%xred_orig, print_indsym=50)

 ABI_MALLOC(symq,(4,2,dtset%nsym))
 timrev=1
 call littlegroup_q(dtset%nsym,dtset%qptn,symq,symrec,dtset%symafm,timrev)

!Initialize the list of perturbations rfpert
 ABI_MALLOC(rfpert,(mpert))
 rfpert(:)=0
 if(dtset%rfphon==1)rfpert(dtset%rfatpol(1):dtset%rfatpol(2))=1

 if(dtset%rfddk==1)rfpert(dtset%natom+1)=1
 if(dtset%rfddk==2)rfpert(dtset%natom+6)=1

 if(dtset%rf2_dkdk/=0) rfpert(dtset%natom+10)=1
 if(dtset%rf2_dkde/=0) rfpert(dtset%natom+11)=1

 if(dtset%rfelfd==1.or.dtset%rfelfd==2)rfpert(dtset%natom+1)=1
 if(dtset%rfelfd==1.or.dtset%rfelfd==3)rfpert(dtset%natom+2)=1

 if(dtset%rfstrs==1.or.dtset%rfstrs==3)rfpert(dtset%natom+3)=1
 if(dtset%rfstrs==2.or.dtset%rfstrs==3)rfpert(dtset%natom+4)=1

 if(dtset%rfuser==1.or.dtset%rfuser==3)rfpert(dtset%natom+6)=1
 if(dtset%rfuser==2.or.dtset%rfuser==3)rfpert(dtset%natom+7)=1

 if(dtset%rfmagn==1) rfpert(dtset%natom+5)=1

 ABI_MALLOC(pertsy,(3,mpert))
 call irreducible_set_pert(indsym,mpert,dtset%natom,dtset%nsym,pertsy,dtset%rfdir,rfpert,symq,symrec,dtset%symrel)

!MR: Desactivate perturbation symmetries for a longwave calculation (TODO)
 if (dtset%prepalw==1) then
   do ipert=1,dtset%natom+6
     do idir=1,3
       if( pertsy(idir,ipert)==-1 ) pertsy(idir,ipert)=1
     end do
   end do
 end if

 npert=0
! ABI_MALLOC(pert_tmp,(3*mpert))

! do ipert=1,mpert
!   do idir=1,3
!     if( rfpert(ipert)==1 .and. dtset%rfdir(idir) == 1 )then
!       if (pertsy(idir,ipert)==1.or.&
!&       (dtset%prepanl == 1.and.ipert == dtset%natom+2).or.&
!&       (dtset%prepgkk == 1.and.ipert <= dtset%natom)  ) then
!         npert = npert+1;
!         pert_tmp(npert) = idir+(ipert-1)*3;
!       else
!         write(msg, '(a,a,i0,a,i0,a,a,a,a,a,a)' )ch10,&
!&         'The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
!&         'symmetric of a previously calculated perturbation.',ch10,&
!&         'So, its SCF calculation is not needed.',ch10
!         call wrtout(std_out,msg,'COLL')
!       end if ! Test of existence of symmetry of perturbation
!     end if ! Test of existence of perturbation
!   end do
! end do

!Initialize rf2dir :
 rf2_dir1(1:3)=dtset%rf2_pert1_dir(1:3)
 rf2_dir2(1:3)=dtset%rf2_pert2_dir(1:3)
!Diagonal terms :
 rf2dir(1) = rf2_dir1(1)*rf2_dir2(1)
 rf2dir(2) = rf2_dir1(2)*rf2_dir2(2)
 rf2dir(3) = rf2_dir1(3)*rf2_dir2(3)
!Upper triangular terms :
 rf2dir(4) = rf2_dir1(2)*rf2_dir2(3)
 rf2dir(5) = rf2_dir1(1)*rf2_dir2(3)
 rf2dir(6) = rf2_dir1(1)*rf2_dir2(2)
!Lower triangular terms :
 rf2dir(7) = rf2_dir1(3)*rf2_dir2(2)
 rf2dir(8) = rf2_dir1(3)*rf2_dir2(1)
 rf2dir(9) = rf2_dir1(2)*rf2_dir2(1)

!Determine existence of pertubations and of pertubation symmetries
!Create array with pertubations which have to be calculated
 ABI_MALLOC(pert_tmp,(2,3*(dtset%natom+6)+18))

 do ipert=1,mpert
   if (ipert<dtset%natom+10) then
     maxidir = 3
     rfdir(1:3) = dtset%rfdir(:)
     rfdir(4:9) = 0
   else
     maxidir = 9
     rfdir(1:9) = rf2dir(:)
   end if
   do idir=1,maxidir
     if( rfpert(ipert)==1 .and. rfdir(idir) == 1 )then
       to_compute_this_pert = 0
       if (ipert>=dtset%natom+10) then
         to_compute_this_pert = 1
       else if ((pertsy(idir,ipert)==1).or.&
          ((dtset%prepanl == 1).and.(ipert == dtset%natom+2)).or.&
          ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
         to_compute_this_pert = 1
       end if
       if (to_compute_this_pert /= 0) then
         npert = npert+1;
         pert_tmp(1,npert) = ipert
         pert_tmp(2,npert) = idir
       else
         write(msg, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
           ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
           ' symmetric of a previously calculated perturbation.',ch10,&
           ' So, its SCF calculation is not needed.',ch10
         call wrtout(std_out,msg)
       end if ! Test of existence of symmetry of perturbation
     end if ! Test of existence of perturbation
   end do
 end do
 ABI_MALLOC(pert_calc,(2,npert))
 do icase=1,npert
   pert_calc(:,icase)=pert_tmp(:,icase)
 end do
 ABI_FREE(pert_tmp)
 ABI_FREE(pertsy)
 ABI_FREE(rfpert)

! Write YAML doc with the list of irreducible perturbations. Example.
!
!--- !IrredPerts
!# List of irreducible perturbations
!irred_perts:
!  - qpt: [ 0.0000000000000000,  0.0000000000000000,  0.0000000000000000]
!    ipert : 1
!    idir  : 1
!  - qpt: [ 0.0000000000000000,  0.0000000000000000,  0.0000000000000000]
!    ipert : 2
!    idir  : 1
!..
 write(std_out,'(a)')"--- !IrredPerts"
 write(std_out,'(a)')'# List of irreducible perturbations'
 write(std_out,'(a)')'irred_perts:'

 do icase=1,npert
!   pert = pert_tmp(icase)

!   if (pert <= dtset%natom*3) then
!     idir = mod(pert, 3)
!     if (idir==0) idir=3
!     ipert=((pert-idir) / 3 + 1)
!   else
!     idir = mod(pert, 3)
!     if (idir==0) idir=3
!     ipert = dtset%natom + ((pert - 3*dtset%natom - 1) / 3) + 1
!   end if
   ipert = pert_calc(1,icase)
   idir = pert_calc(2,icase)

   write(std_out,'(a,3(f20.16,a))')"   - qpt: [ ",dtset%qptn(1),", ", dtset%qptn(2),", ", dtset%qptn(3),"]"
   write(std_out,'(a,i0)')"     ipert: ",ipert
   write(std_out,'(a,i0)')"     idir: ",idir
 end do

 write(std_out,'(a)')"..."

! ABI_MALLOC(pert_calc,(npert))
! do icase=1,npert
!   pert_calc(icase) = pert_tmp(icase)
! end do

 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)

 ABI_MALLOC(nkpt_rbz,(npert))
 ABI_MALLOC(indkpt1,(dtset%nkpt,npert))
 indkpt1=0

 do icase=1,npert
!   if (pert_calc(icase) <= dtset%natom*3) then
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!     ipert=( (pert_calc(icase)-idir) / 3 + 1)
!   else
!     ipert = dtset%natom + ((pert_calc(icase) - 3*dtset%natom - 1) / 3) + 1
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!   end if
   ipert = pert_calc(1,icase)
   idir = pert_calc(2,icase)

   ABI_MALLOC(symrl1_tmp,(3,3,dtset%nsym))
   ABI_MALLOC(symaf1,(dtset%nsym))
   ABI_MALLOC(tnons1_tmp,(3,dtset%nsym))
!  MJV TODO: check whether prepgkk should be used here
   if (dtset%prepanl /= 1 .and. dtset%berryopt /=4 .and. dtset%berryopt /=6 .and. dtset%berryopt /=7 .and. &
&   dtset%berryopt /=14 .and. dtset%berryopt /=16 .and. dtset%berryopt /=17) then   !!HONG
     call littlegroup_pert(gprimd,idir,indsym,std_out,ipert,dtset%natom,dtset%nsym,nsym1,2,&
&     dtset%symafm,symaf1,symq,symrec,&
&     dtset%symrel,symrl1_tmp,0,dtset%tnons,tnons1_tmp)
   else
     nsym1 = 1
   end if
   ABI_FREE(tnons1_tmp)
   ABI_FREE(symaf1)

   ABI_MALLOC(symrc1,(3,3,nsym1))
   ABI_MALLOC(symrl1,(3,3,nsym1))
   symrl1(:,:,1:nsym1)=symrl1_tmp(:,:,1:nsym1)
   ABI_FREE(symrl1_tmp)
   do isym=1,nsym1
     call mati3inv(symrl1(:,:,isym),symrc1(:,:,isym))
   end do
   ABI_FREE(symrl1)

   ABI_MALLOC(wtk_folded,(dtset%nkpt))
   ABI_MALLOC(bz2ibz_smap, (6, dtset%nkpt))
   timrev_pert=timrev
   if(dtset%ieig2rf>0) then
     call symkpt(0,gmet,indkpt1(:,icase),std_out,dtset%kptns,dtset%nkpt,nkpt_rbz(icase),&
&     1,symrc1,0,dtset%wtk,wtk_folded, bz2ibz_smap, xmpi_comm_self)
   else
!    For the time being, the time reversal symmetry is not used
!    for ddk, elfd, mgfd perturbations.
     if(ipert==dtset%natom+1 .or. ipert==dtset%natom+10 .or. ipert==dtset%natom+11 .or. &
        ipert==dtset%natom+2 .or. dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7  &
        .or. dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17 )timrev_pert=0  !!HONG
     call symkpt(0,gmet,indkpt1(:,icase),std_out,dtset%kptns,dtset%nkpt,nkpt_rbz(icase),&
     nsym1,symrc1,timrev_pert,dtset%wtk,wtk_folded, bz2ibz_smap, xmpi_comm_self)
   end if
   ABI_FREE(bz2ibz_smap)
   ABI_FREE(wtk_folded)
   ABI_FREE(symrc1)
 end do

 ABI_MALLOC(nband_rbz,(maxval(nkpt_rbz)*dtset%nsppol,npert))
 nband_rbz=zero
 do icase=1,npert
   do isppol=1,dtset%nsppol
     ikpt1=1
     do ikpt=1,dtset%nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
!      Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
       if(ikpt1/=nkpt_rbz(icase)+1)then
         if(ikpt==indkpt1(ikpt1,icase))then
           nband_rbz(ikpt1+(isppol-1)*nkpt_rbz(icase),icase)=nband_k
         end if
       end if
     end do
   end do

 end do

 ABI_FREE(indkpt1)
 ABI_FREE(symq)
 ABI_FREE(symrec)
 ABI_FREE(indsym)
 ABI_FREE(pert_calc)

end subroutine dtset_get_npert_rbz
!!***

!!****f* m_dtset/dtset_testsusmat
!! NAME
!! dtset_testsusmat
!!
!! FUNCTION
!! Test wether a new susceptibility matrix and/or a new dielectric matrix must be computed
!! and return the logical result
!!
!! INPUTS
!! dielop: option for the computation of the dielectric matrix
!! dtset:
!! istep: number of the current SCF cycle
!!
!! OUTPUT
!! compute:
!!  * if dielop >= 1 and istep == 1 => TRUE
!!  * if dielop >= 2 and istep == dtset%dielstrt => TRUE
!!  * if (dtset%iprcel >= 140 and <=170) depends on the periodicity modulo 10 of istep and iprcel
!!  * otherwise FALSE
!!
!! PARENTS
!!      prcref,prcref_PMA,vtorho
!!
!! CHILDREN
!!
!! SOURCE

logical function dtset_testsusmat(dtset, dielop, dielstrt, istep) result(compute)

!Arguments-------------------------------
!scalars
 integer,intent(in) :: dielop,dielstrt,istep
 !logical,intent(out) :: compute
 class(dataset_type),intent(in) :: dtset

! *********************************************************************

 compute=.FALSE.
 if((dtset%iprcel >= 140).and.(dtset%iprcel<=170)) then
   if(modulo(dtset%iprcel,10).ne.0) then
     compute=(modulo(istep,modulo(dtset%iprcel,10))==0)
   else
     compute=(modulo(istep,10)==0)
   end if
 end if
 if (istep==1 .and. dielop>=2) compute=.TRUE.
 if (istep==dielstrt .and. dielop>=1) compute=.TRUE.

end function dtset_testsusmat
!!***

!!****f* m_dtset/macroin
!! NAME
!! macroin
!!
!! FUNCTION
!! Treat "macro" input variables, that can:
!!
!!      - initialize several other input variables for one given dataset
!!      - initialize several other input variables for a set of datasets.
!1
!! Note that the treatment of these different types of macro input variables is different.
!! Documentation of such input variables is very important, including the
!! proper echo, in the output file, of what such input variables have done.
!!
!! Important information : all the "macro" input variables should be properly
!! identifiable to be so, and it is proposed to make them start with the string "macro".
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!  ecut_tmp(3,2,10)= possible ecut values as read in psp files
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=contains all input variables, some of which are given a value here.
!!   The dataset with number 0 should NOT be modified in the present routine.
!!
!! PARENTS
!!      m_common
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine macroin(dtsets,ecut_tmp,lenstr,ndtset_alloc,string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc,lenstr
 character(len=*),intent(inout) :: string
!arrays
 real(dp),intent(in) :: ecut_tmp(3,2,10)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i ziontypat

!Local variables -------------------------------
!scalars
 integer :: idtset,iatom,jdtset,marr,tread
!!arrays
 integer,allocatable :: intarr(:)
 real(dp) :: ecutmax(3),ecutdgmax(3)
 real(dp),allocatable :: dprarr(:)
 character(len=500) :: msg
!******************************************************************

 do idtset=1,ndtset_alloc
   jdtset=dtsets(idtset)%jdtset
   if (dtsets(idtset)%macro_uj>0) then
     dtsets(idtset)%irdwfk   = 1        ! preconverged wave function compulsory
!    dtsets(idtset)%nline    = maxval((/ int(dtsets(idtset)%natom/2) , 6 /))   ! using default value: \DeltaU< 1%
!    dtsets(idtset)%nnsclo   = 4        ! using default value: \DeltaU< 1%
     dtsets(idtset)%tolvrs   = 10d-8    ! convergence on the potential; 10d-8^= 10d-5 on occupation
     dtsets(idtset)%diemix   = 0.45_dp  ! fastest convergence: dn= E^(-istep * 0.229 )
     dtsets(idtset)%dmatpuopt= 3        ! normalization of the occupation operator
!    dtsets(idtset)%nstep    = 255      ! expected convergence after 10 \pm 3, 30 as in default normally suficient
!    dtsets(idtset)%iscf     = 17       ! mixing on potential, 17: default for PAW
   end if ! macro_uj

  !Read parameters
   marr=dtsets(idtset)%npsp;if (dtsets(idtset)%npsp<3) marr=3
   marr=max(marr,dtsets(idtset)%nimage)
   ABI_MALLOC(intarr,(marr))
   ABI_MALLOC(dprarr,(marr))

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),"accuracy",tread,'INT')

   ecutmax=-one
   ecutdgmax=-one
   do iatom=1,dtsets(idtset)%natom
     ecutmax(:)=max(ecutmax(:),ecut_tmp(:,1,dtsets(idtset)%typat(iatom)))
     ecutdgmax(:)=max(ecutdgmax(:),ecut_tmp(:,2,dtsets(idtset)%typat(iatom)))
   end do

   if(tread==1) then
     dtsets(idtset)%accuracy=intarr(1)
     if (dtsets(idtset)%accuracy==1) then
       if (ecutmax(1)>zero) dtsets(idtset)%ecut=ecutmax(1)
       if (ecutdgmax(1)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(1)
       dtsets(idtset)%boxcutmin=1.5_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.5_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=10
         dtsets(idtset)%pawnhatxc=0
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol3
       dtsets(idtset)%tolmxf=1.0d-3
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=1
       dtsets(idtset)%timopt=0
       dtsets(idtset)%npulayit=4
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=0
       dtsets(idtset)%prtden=0
     else if (dtsets(idtset)%accuracy==2) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=1.8_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.8_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=7
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol5
       dtsets(idtset)%tolmxf=5.0d-4
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=1
       dtsets(idtset)%timopt=0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=0
       dtsets(idtset)%prtden=0
     else if (dtsets(idtset)%accuracy==3) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=1.8_dp
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=1.8_dp
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=7
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol7
       dtsets(idtset)%tolmxf=1.0d-4
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==4) then
       if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
       if (ecutdgmax(3)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(3)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=1
         dtsets(idtset)%pawmixdg=0
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol9
       dtsets(idtset)%tolmxf=5.0d-5
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=7
       dtsets(idtset)%nstep=30
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==5) then
       if (ecutmax(2)>zero) dtsets(idtset)%ecut=ecutmax(2)
       if (ecutdgmax(2)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(2)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=2
         dtsets(idtset)%pawmixdg=1
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol10
       dtsets(idtset)%tolmxf=1.0d-6
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=15
       dtsets(idtset)%nstep=50
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     else if (dtsets(idtset)%accuracy==6) then
       if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
       if (ecutdgmax(3)>zero.and.dtsets(idtset)%usepaw==1) dtsets(idtset)%pawecutdg=ecutdgmax(3)
       dtsets(idtset)%boxcutmin=two
       if (dtsets(idtset)%usepaw==1) then
         dtsets(idtset)%bxctmindg=two
         dtsets(idtset)%pawxcdev=2
         dtsets(idtset)%pawmixdg=1
         dtsets(idtset)%pawovlp=5
         dtsets(idtset)%pawnhatxc=1
         dtsets(idtset)%mqgriddg=0
       end if
       dtsets(idtset)%mqgrid=0
       dtsets(idtset)%tolimg=5.0d-5
       dtsets(idtset)%tolvrs=tol12
       dtsets(idtset)%tolmxf=1.0d-6
       dtsets(idtset)%toldff=zero
       dtsets(idtset)%optforces=2
       dtsets(idtset)%timopt=1
       if(xmpi_paral==1) dtsets(idtset)%timopt = 0
       dtsets(idtset)%npulayit=15
       dtsets(idtset)%nstep=50
       dtsets(idtset)%prteig=1
       dtsets(idtset)%prtden=1
     elseif(dtsets(idtset)%accuracy>6)then
       write(msg, '(a,a,a)' )&
         'accuracy >6 is forbidden !',ch10,&
         'Action: check your input data file.'
       MSG_ERROR(msg)
     end if
   else
     if (ecutmax(3)>zero) dtsets(idtset)%ecut=ecutmax(3)
   end if
   ABI_FREE(intarr)
   ABI_FREE(dprarr)
 end do

end subroutine macroin
!!***

!!****f* m_dtset/macroin2
!! NAME
!! macroin2
!!
!! FUNCTION
!! Treat "macro" input variables, that can :
!! - initialize several other input variables for one given dataset
!! - initialize several other input variables for a set of datasets.
!! Note that the treatment of these different types of macro input variables is different.
!! Documentation of such input variables is very important, including the
!! proper echo, in the output file, of what such input variables have done.
!!
!! Important information : all the "macro" input variables should be properly
!! identifiable to be so, and it is proposed to make them start with the string "macro".
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=contains all input variables, some of which are given a value here.
!!   The dataset with number 0 should NOT be modified in the present routine.
!!
!! PARENTS
!!      m_common
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine macroin2(dtsets,ndtset_alloc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------------
!scalars
 integer :: idtset,pawujat

!******************************************************************

 do idtset=1,ndtset_alloc
!  Set first PAW+U atom to perform atomic level shift
   if (dtsets(idtset)%typat(1)==0) cycle
   pawujat=dtsets(idtset)%pawujat
   pawujat=pawujat-count(dtsets(idtset)%lpawu( dtsets(idtset)%typat( 1:pawujat ))<0)

   if (dtsets(idtset)%macro_uj>0) then
     ! Level shift atom with amplitude pawujv
     dtsets(idtset)%atvshift(:,:,pawujat)=dtsets(idtset)%pawujv
     ! Case level shift only on one spin channel
     if ((dtsets(idtset)%macro_uj==2.or.dtsets(idtset)%macro_uj==3).and.dtsets(idtset)%nsppol==2) then
       dtsets(idtset)%atvshift(:,2,pawujat)=0_dp
     end if
   end if ! macro_uj
 end do

end subroutine macroin2
!!***

!!****f* ABINIT/chkvars
!! NAME
!! chkvars
!!
!! FUNCTION
!!  Examines the input string, to check whether all names are allowed.
!!
!! INPUTS
!!  string*(*)=string of character
!!   the string (with upper case) from the input file, to which the XYZ data is (possibly) appended
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit,m_multibinit_driver,m_multibinit_manager
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine chkvars(string)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string

!Local variables-------------------------------
!scalars
 integer,parameter :: protocol1=1
 character(len=100) :: list_logicals,list_strings
 character(len=10000) :: list_vars

!************************************************************************


!Here, list all admitted variable names (max 10 per line, to fix the ideas)
!<ABINIT_VARS>
!A
 list_vars=                 ' accuracy acell adpimd adpimd_gamma'
 list_vars=trim(list_vars)//' algalch amu analyze_anh_pot angdeg asr atvshift autoparal'
 list_vars=trim(list_vars)//' auxc_ixc auxc_scal awtr'
!B
 list_vars=trim(list_vars)//' bandpp bdberry bdeigrf bdgw berryopt berrysav berrystep bfield bmass'
 list_vars=trim(list_vars)//' boxcenter boxcutmin brav brvltt builtintest'
 list_vars=trim(list_vars)//' bound_SPCoupling bound_anhaStrain bound_cell bound_cutoff'
 list_vars=trim(list_vars)//' bound_maxCoeff bound_model bound_rangePower bound_step bound_temp'
 list_vars=trim(list_vars)//' bs_algorithm bs_calctype bs_coulomb_term bs_coupling'
 list_vars=trim(list_vars)//' bs_interp_kmult bs_interp_m3_width bs_interp_method bs_interp_mode bs_interp_prep'
 list_vars=trim(list_vars)//' bs_interp_rl_nb bs_eh_cutoff bs_exchange_term bs_freq_mesh'
 list_vars=trim(list_vars)//' bs_haydock_niter bs_haydock_tol bs_hayd_term'
 list_vars=trim(list_vars)//' bs_loband bs_nstates'
 list_vars=trim(list_vars)//' bxctmindg'
!C
 list_vars=trim(list_vars)//' cd_customnimfrqs cd_frqim_method cd_full_grid cd_imfrqs'
 list_vars=trim(list_vars)//' cd_halfway_freq cd_max_freq cd_subset_freq'
 list_vars=trim(list_vars)//' chrgat charge chempot chkdilatmx chkexit chkprim'
 list_vars=trim(list_vars)//' chksymbreak chneut cineb_start coefficients constraint_kind cpus cpum cpuh'
!D
 list_vars=trim(list_vars)//' ddamp ddb_ngqpt ddb_shiftq'
 list_vars=trim(list_vars)//' delayperm densfor_pred densty dfield'
 list_vars=trim(list_vars)//' dfpt_sciss diecut diegap dielam dielng diemac'
 list_vars=trim(list_vars)//' diemix diemixmag diismemory dilatmx dipdip dipdip_prt dipdip_range'
 list_vars=trim(list_vars)//' dmatpawu dmatpuopt dmatudiag'
 list_vars=trim(list_vars)//' dmftbandi dmftbandf dmftctqmc_basis'
 list_vars=trim(list_vars)//' dmftctqmc_check dmftctqmc_correl dmftctqmc_gmove'
 list_vars=trim(list_vars)//' dmftctqmc_grnns dmftctqmc_meas dmftctqmc_mrka'
 list_vars=trim(list_vars)//' dmftctqmc_mov dmftctqmc_order dmftctqmc_triqs_nleg'
 list_vars=trim(list_vars)//' dmftcheck dmftqmc_l dmftqmc_n dmftqmc_seed dmftqmc_therm'
 list_vars=trim(list_vars)//' dmft_charge_prec dmft_dc dmft_entropy dmft_iter dmft_kspectralfunc'
 list_vars=trim(list_vars)//' dmft_mxsf dmft_nlambda dmft_nwli dmft_nwlo'
 list_vars=trim(list_vars)//' dmft_occnd_imag dmft_read_occnd dmft_rslf dmft_solv'
!list_vars=trim(list_vars)//' dmft_tolfreq dmft_tollc dmft_t2g dmft_x2my2d'
 list_vars=trim(list_vars)//' dmft_tolfreq dmft_tollc dmft_t2g'
 list_vars=trim(list_vars)//' dosdeltae dtion dynamics dynimage'
 list_vars=trim(list_vars)//' dvdb_add_lr dvdb_ngqpt dvdb_qcache_mb dvdb_qdamp dvdb_rspace_cell'
 list_vars=trim(list_vars)//' dyn_chksym dyn_tolsym'
 list_vars=trim(list_vars)//' d3e_pert1_atpol d3e_pert1_dir d3e_pert1_elfd d3e_pert1_phon'
 list_vars=trim(list_vars)//' d3e_pert2_atpol d3e_pert2_dir d3e_pert2_elfd d3e_pert2_phon'
 list_vars=trim(list_vars)//' d3e_pert2_strs'
 list_vars=trim(list_vars)//' d3e_pert3_atpol d3e_pert3_dir d3e_pert3_elfd d3e_pert3_phon'
!E
 list_vars=trim(list_vars)//' ecut ecuteps ecutsigx ecutsm ecutwfn effmass_free efmas'
 list_vars=trim(list_vars)//' efmas_bands efmas_calc_dirs efmas_deg efmas_deg_tol'
 list_vars=trim(list_vars)//' efmas_dim efmas_dirs efmas_n_dirs efmas_ntheta'
 list_vars=trim(list_vars)//' efield einterp elph2_imagden energy_reference enunit'
 list_vars=trim(list_vars)//' eph_doping eph_ecutosc eph_extrael eph_fermie eph_frohlich eph_frohlichm eph_fsewin eph_fsmear '
 list_vars=trim(list_vars)//' eph_np_pqbks'
  ! XG20200321, please provide testing for eph_np_pqbks
  ! Well, eph_np_pqbks cannot be tested with the present infrastructure because it's a MPI-related variable
  ! and all the tests in the paral and mpiio directory are done with a single input file
  ! whereas EPH requires GS + DFPT + MRGDV + MRGDDB + TESTS_MULTIPLES_PROCS
 list_vars=trim(list_vars)//' eph_intmeth eph_mustar eph_ngqpt_fine'
 list_vars=trim(list_vars)//' eph_phrange eph_tols_idelta '
 list_vars=trim(list_vars)//' eph_restart eph_stern eph_task eph_transport eph_use_ftinterp'
 list_vars=trim(list_vars)//' eshift esmear exchmix exchn2n3d extrapwf'
!F
 list_vars=trim(list_vars)//' fband fermie_nest'
 list_vars=trim(list_vars)//' fftalg fftcache fftgw'
 list_vars=trim(list_vars)//' fit_anhaStrain fit_bancoeff fit_coeff fit_cutoff fit_fixcoeff'
 list_vars=trim(list_vars)//' fit_EFS'
 list_vars=trim(list_vars)//' fit_generateCoeff fit_iatom fit_initializeData fit_nbancoeff fit_ncoeff fit_nfixcoeff'
 list_vars=trim(list_vars)//' fit_rangePower fit_SPCoupling fit_SPC_maxS fit_tolMSDE fit_tolMSDF fit_tolMSDFS fit_tolMSDS'
 list_vars=trim(list_vars)//' fockoptmix focktoldfe fockdownsampling fock_icutcoul'
 list_vars=trim(list_vars)//' freqim_alpha freqremax freqremin freqspmax'
 list_vars=trim(list_vars)//' freqspmin friction frzfermi fxcartfactor'
 list_vars=trim(list_vars)//' freqspmin friction frzfermi fxcartfactor'
 list_vars=trim(list_vars)//' f4of2_sla f6of2_sla'
!G
 list_vars=trim(list_vars)//' ga_algor ga_fitness ga_n_rules ga_opt_percent ga_rules'
 list_vars=trim(list_vars)//' genafm getbscoup getbseig getbsreso getcell'
 list_vars=trim(list_vars)//' getddb getddb_filepath getden_filepath getddk'
 list_vars=trim(list_vars)//' getdelfd getdkdk getdkde getden getdvdb getdvdb_filepath'
 list_vars=trim(list_vars)//' getefmas getkerange_filepath getgam_eig2nkq'
 list_vars=trim(list_vars)//' gethaydock getocc getpawden getpot_filepath getsigeph_filepath getqps getscr getscr_filepath'
 list_vars=trim(list_vars)//' getwfkfine getwfkfine_filepath getsuscep'
 list_vars=trim(list_vars)//' getvel getwfk getwfk_filepath getwfq getwfq_filepath getxcart getxred'
 list_vars=trim(list_vars)//' get1den get1wf goprecon goprecprm'
 list_vars=trim(list_vars)//' gpu_devices gpu_linalg_limit gwaclowrank gwcalctyp gwcomp gwencomp gwgamma gwmem'
 list_vars=trim(list_vars)//' gwpara gwrpacorr gw_customnfreqsp'
 list_vars=trim(list_vars)//' gw_frqim_inzgrid gw_frqre_inzgrid gw_frqre_tangrid gw_freqsp'
 list_vars=trim(list_vars)//' gw_invalid_freq'
 list_vars=trim(list_vars)//' gw_icutcoul'
 list_vars=trim(list_vars)//' gw_qprange gw_nqlwl gw_nstep gw_qlwl'
 list_vars=trim(list_vars)//' gw_sctype gw_sigxcore gw_toldfeig'
 list_vars=trim(list_vars)//' gwls_stern_kmax gwls_kmax_complement gwls_kmax_poles'
 list_vars=trim(list_vars)//' gwls_kmax_analytic gwls_kmax_numeric'
 list_vars=trim(list_vars)//' gwls_list_proj_freq gwls_nseeds gwls_n_proj_freq gwls_recycle'
 list_vars=trim(list_vars)//' gwls_first_seed gwls_model_parameter gwls_npt_gauss_quad'
 list_vars=trim(list_vars)//' gwls_diel_model gwls_print_debug gwls_band_index gwls_exchange gwls_correlation'
!H
 list_vars=trim(list_vars)//' hmcsst hmctt hyb_mixing hyb_mixing_sr hyb_range_dft hyb_range_fock'
!I
 list_vars=trim(list_vars)//' iatcon iatfix iatfixx iatfixy iatfixz iatsph'
 list_vars=trim(list_vars)//' iboxcut icoulomb icutcoul ieig2rf'
 list_vars=trim(list_vars)//' imgmov imgwfstor inclvkb indata_prefix intxc iomode ionmov iqpt'
 list_vars=trim(list_vars)//' iprcel iprcfc irandom irdbscoup'
 list_vars=trim(list_vars)//' irdbseig irdbsreso irdddb irdddk irdden irddvdb irdefmas'
 list_vars=trim(list_vars)//' irdhaydock irdpawden irdqps'
 list_vars=trim(list_vars)//' irdscr irdsuscep irdwfk irdwfq ird1den'
 list_vars=trim(list_vars)//' irdwfkfine'
 list_vars=trim(list_vars)//' ird1wf iscf isecur istatimg istatr'
 list_vars=trim(list_vars)//' istatshft istwfk ixc ixc_sigma ixcpositron ixcrot'
 list_vars=trim(list_vars)//' irdvdw'
!J
 list_vars=trim(list_vars)//' jdtset jellslab jfielddir jpawu'
!K
 list_vars=trim(list_vars)//' kberry kpt kptbounds kptgw'
 list_vars=trim(list_vars)//' kptnrm kptopt kptrlatt kptrlen kssform'
!L
 list_vars=trim(list_vars)//' latt_friction latt_taut latt_taup latt_compressibility latt_mask'
 list_vars=trim(list_vars)//' ldaminushalf lexexch localrdwf lpawu'
 list_vars=trim(list_vars)//' lotf_classic lotf_nitex lotf_nneigx lotf_version'
 list_vars=trim(list_vars)//' lw_flexo lw_qdrpl'
!M
 list_vars=trim(list_vars)//' max_ncpus macro_uj maxestep maxnsym mdf_epsinf mdtemp mdwall'
 list_vars=trim(list_vars)//' magconon magcon_lambda mbpt_sciss'
 list_vars=trim(list_vars)//' mep_mxstep mep_solver mem_test mixalch mixprec mixesimgf'
 list_vars=trim(list_vars)//' mqgrid mqgriddg'
!N
 list_vars=trim(list_vars)//' natcon natfix natfixx natfixy natfixz'
 list_vars=trim(list_vars)//' natom natrd natsph natsph_extra natvshift nband nbandkss nbandhf'
 list_vars=trim(list_vars)//' ncell ncoeff nbdblock nbdbuf nberry nconeq nc_xccc_gspace'
 list_vars=trim(list_vars)//' nctime ndivk ndivsm ndtset neb_algo neb_spring'
 list_vars=trim(list_vars)//' nfreqim nfreqre nfreqsp ngfft ngfftdg'
 list_vars=trim(list_vars)//' ngkpt ngqpt nimage nkpath nkpt nkptgw nkpthf'
 list_vars=trim(list_vars)//' nline nloc_alg nloc_mem nnos nnsclo nnsclohf'
 list_vars=trim(list_vars)//' nobj nomegasf nomegasi nomegasrd nonlinear_info noseinert npband'
 list_vars=trim(list_vars)//' npfft nphf nph1l npimage npkpt nppert npsp npspinor'
 list_vars=trim(list_vars)//' npulayit npvel npwkss'
 list_vars=trim(list_vars)//' np_slk nqpt nqptdm nscforder nshiftk nshiftq nqshft'
 list_vars=trim(list_vars)//' nspden nspinor nsppol nstep nsym'
 list_vars=trim(list_vars)//' ntime ntimimage ntypalch ntypat nucdipmom nwfshist nzchempot'
!O
 list_vars=trim(list_vars)//' objaat objbat objaax objbax objan objbn objarf'
 list_vars=trim(list_vars)//' objbrf objaro objbro objatr objbtr occ'
 list_vars=trim(list_vars)//' occopt omegasimax omegasrdmax optcell optdriver optforces'
 list_vars=trim(list_vars)//' optnlxccc optstress orbmag ortalg'
 list_vars=trim(list_vars)//' opt_effpot opt_ncoeff opt_coeff output_file outdata_prefix'
!P
 list_vars=trim(list_vars)//' papiopt paral_atom paral_kgb paral_rf pawcpxocc pawcross'
 list_vars=trim(list_vars)//' pawecutdg pawfatbnd pawlcutd pawlmix'
 list_vars=trim(list_vars)//' pawmixdg pawnhatxc pawnphi pawntheta pawnzlm pawoptmix pawoptosc pawovlp'
 list_vars=trim(list_vars)//' pawprtdos pawprtvol pawprtwf pawprt_b pawprt_k pawspnorb pawstgylm'
 list_vars=trim(list_vars)//' pawsushat pawujat pawujrad pawujv'
 list_vars=trim(list_vars)//' pawusecp pawxcdev pimass pimd_constraint'
 list_vars=trim(list_vars)//' ph_intmeth ph_ndivsm ph_ngqpt ph_nqpath ph_nqshift ph_qpath'
 list_vars=trim(list_vars)//' ph_qshift ph_smear ph_wstep pitransform'
 list_vars=trim(list_vars)//' plowan_bandi plowan_bandf plowan_compute plowan_iatom plowan_it plowan_lcalc'
 list_vars=trim(list_vars)//' plowan_natom plowan_nbl plowan_nt plowan_projcalc plowan_realspace'
 list_vars=trim(list_vars)//' polcen posdoppler positron posnstep posocc postoldfe postoldff'
 list_vars=trim(list_vars)//' ppmfrq ppmodel pp_dirpath'
 list_vars=trim(list_vars)//' prepalw prepanl prepgkk'
 list_vars=trim(list_vars)//' prtatlist prtbbb prtbltztrp prtcif prtden'
 list_vars=trim(list_vars)//' prtdensph prtdipole prtdos prtdosm prtebands prtefg prtefmas prteig prteliash prtelf'
 list_vars=trim(list_vars)//' prtfc prtfull1wf prtfsurf prtgden prtgeo prtgsr prtgkk prtkden prtkpt prtlden'
 list_vars=trim(list_vars)//' prt_model prtnabla prtnest prtphbands prtphdos prtphsurf prtposcar'
 list_vars=trim(list_vars)//' prtprocar prtpot prtpsps'
 list_vars=trim(list_vars)//' prtspcur prtstm prtsuscep prtvclmb prtvha prtvdw prtvhxc prtkbff'
 list_vars=trim(list_vars)//' prtvol prtvpsp prtvxc prtwant prtwf prtwf_full prtxml prt1dm'
 list_vars=trim(list_vars)//' pseudos ptcharge'
 list_vars=trim(list_vars)//' pvelmax pw_unbal_thresh'
!Q
 list_vars=trim(list_vars)//' q1shft qmass qprtrb qpt qptdm qptnrm qph1l'
 list_vars=trim(list_vars)//' qptopt qptrlatt quadmom'
!R
 list_vars=trim(list_vars)//' random_atpos ratsm ratsph ratsph_extra rcut'
 list_vars=trim(list_vars)//' recefermi recgratio recnpath recnrec recptrott recrcut rectesteg rectolden'
 list_vars=trim(list_vars)//' red_dfield red_efield red_efieldbar restartxf rfasr'
 list_vars=trim(list_vars)//' rfatpol rfddk rfdir rfelfd rfmagn rfmeth rfphon'
 list_vars=trim(list_vars)//' rfstrs rfuser rf2_dkdk rf2_dkde rf2_pert1_dir rf2_pert2_dir rhoqpmix rprim'
 !These input parameters are obsolete (keep them for compatibility)
 list_vars=trim(list_vars)//' rf1atpol rf1dir rf1elfd rf1phon'
 list_vars=trim(list_vars)//' rf2atpol rf2dir rf2elfd rf2phon rf2strs'
 list_vars=trim(list_vars)//' rf3atpol rf3dir rf3elfd rf3phon'
!S
 list_vars=trim(list_vars)//' scalecart shiftk shiftq signperm'
 list_vars=trim(list_vars)//' sel_EFS'
 list_vars=trim(list_vars)//' sigma_bsum_range sigma_erange sigma_ngkpt sigma_nshiftk sigma_shiftk'
!MS Variables for SCALE-UP
!This is only for the developer version, not for the production version. So, was commented.
! @Marcus: simply uncomment these lines in v9.1 (not v9.0 !), and continue to develop without worrying.
!list_vars=trim(list_vars)//' scup_elec_model scup_ksamp scup_tcharge scup_initorbocc scup_ismagnetic'
!list_vars=trim(list_vars)//' scup_istddft scup_printbands scup_printgeom scup_printeigv scup_printeltic '
!list_vars=trim(list_vars)//' scup_printorbocc scup_printniter scup_nspeck scup_speck scup_ndivsm'
!list_vars=trim(list_vars)//' scup_scfmixing scup_scfthresh scup_startpulay scup_maxscfstep'
!list_vars=trim(list_vars)//' scup_smearing scup_freezden'
!End SCALE-UP variables
 list_vars=trim(list_vars)//' slabwsrad slabzbeg slabzend slk_rankpp smdelta so_psp'
 list_vars=trim(list_vars)//' slc_coupling'
 list_vars=trim(list_vars)//' spbroad spgaxor spgorig spgroup spgroupma'
 list_vars=trim(list_vars)//' spin_calc_correlation_obs spin_calc_thermo_obs spin_calc_traj_obs'
 list_vars=trim(list_vars)//' spin_damping'
 list_vars=trim(list_vars)//' spin_dipdip spin_dt spin_dynamics '
 list_vars=trim(list_vars)//' spin_init_orientation spin_init_qpoint spin_init_rotate_axis spin_init_state'
 list_vars=trim(list_vars)//' spin_mag_field spin_nctime spin_ntime spin_ntime_pre'
 list_vars=trim(list_vars)//' spin_n1l spin_n2l spin_projection_qpoint'
 list_vars=trim(list_vars)//' spin_sia_add spin_sia_k1amp spin_sia_k1dir'
 list_vars=trim(list_vars)//' spin_temperature spin_temperature_end'
 list_vars=trim(list_vars)//' spin_temperature_nstep spin_temperature_start spin_tolavg'
 list_vars=trim(list_vars)//' spin_tolvar spin_var_temperature spin_write_traj'
 list_vars=trim(list_vars)//' spinat spinmagntarget spmeth'
 list_vars=trim(list_vars)//' spnorbscl stmbias strfact string_algo strprecon strtarget'
 list_vars=trim(list_vars)//' supercell_latt symafm symchi symdynmat symmorphi symrel symsigma symv1scf'
 list_vars=trim(list_vars)//' structure '
!T
 list_vars=trim(list_vars)//' td_maxene td_mexcit tfkinfunc temperature test_effpot test_prt_ph tfw_toldfe tim1rev timopt'
 list_vars=trim(list_vars)//' tmesh tmpdata_prefix transport_ngkpt'
 list_vars=trim(list_vars)//' tl_nprccg tl_radius tnons toldfe tolmxde toldff tolimg tolmxf tolrde tolrff tolsym'
 list_vars=trim(list_vars)//' tolvrs tolwfr tphysel ts_option tsmear typat'
!U
 list_vars=trim(list_vars)//' ucrpa ucrpa_bands ucrpa_window udtset upawu usepead usedmatpu '
 list_vars=trim(list_vars)//' usedmft useexexch usekden use_nonscf_gkk usepawu usepotzero'
 list_vars=trim(list_vars)//' useria userib useric userid userie'
 list_vars=trim(list_vars)//' userra userrb userrc userrd userre'
 list_vars=trim(list_vars)//' usewvl usexcnhat useylm use_gemm_nonlop use_gpu_cuda use_slk use_yaml'
!V
 list_vars=trim(list_vars)//' vaclst vacnum vacuum vacwidth vcutgeo'
 list_vars=trim(list_vars)//' vdw_nfrag vdw_supercell'
 list_vars=trim(list_vars)//' vdw_tol vdw_tol_3bt vdw_typfrag vdw_xc'
 list_vars=trim(list_vars)//' vdw_df_acutmin vdw_df_aratio vdw_df_damax'
 list_vars=trim(list_vars)//' vdw_df_damin vdw_df_dcut vdw_df_dratio'
 list_vars=trim(list_vars)//' vdw_df_dsoft vdw_df_gcut'
 list_vars=trim(list_vars)//' vdw_df_ndpts vdw_df_ngpts vdw_df_nqpts'
 list_vars=trim(list_vars)//' vdw_df_nrpts vdw_df_nsmooth vdw_df_phisoft vdw_df_qcut'
 list_vars=trim(list_vars)//' vdw_df_qratio vdw_df_rcut vdw_df_rsoft'
 list_vars=trim(list_vars)//' vdw_df_threshold vdw_df_tolerance'
 list_vars=trim(list_vars)//' vdw_df_tweaks vdw_df_zab'
 list_vars=trim(list_vars)//' vel vel_cell vis vprtrb'
!W
 list_vars=trim(list_vars)//' wfmix wfoptalg wtatcon wtk wtq'
 list_vars=trim(list_vars)//' wvl_bigdft_comp wvl_crmult wvl_frmult wvl_hgrid wvl_ngauss wvl_nprccg'
 list_vars=trim(list_vars)//' w90iniprj w90prtunk'
!X
 list_vars=trim(list_vars)//' xcart xc_denpos xc_tb09_c xred xredsph_extra xyzfile'
!Y
!Z
 list_vars=trim(list_vars)//' zcut zeemanfield znucl'

!Logical input variables
 list_logicals=' SpinPolarized'

!String input variables
 list_strings=' XCname wfk_task'
!</ABINIT_VARS>

!Extra token, also admitted:
!<ABINIT_UNITS>
 list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha'
 list_vars=trim(list_vars)//' Hartree Hartrees K nm Ry Rydberg Rydbergs S Sec Second T Tesla'
!</ABINIT_UNITS>

!<ABINIT_OPERATORS>
 list_vars=trim(list_vars)//' sqrt end'
!</ABINIT_OPERATORS>

 ! Transform to upper case
 call inupper(list_vars)
 call inupper(list_logicals)
 call inupper(list_strings)

 call chkvars_in_string(protocol1, list_vars, list_logicals, list_strings, string)

end subroutine chkvars
!!***

end module m_dtset
!!***
