!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_abitypes
!! NAME
!! defs_abitypes
!!
!! FUNCTION
!! This module contains definitions of high-level structured datatypes for the ABINIT package.
!!
!! If you are sure a new high-level structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are well documented, it is a shame ...).
!! Do not forget: you will likely be the major winner if you document properly.
!!
!! Proper documentation of a structured datatype means:
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!  (5) Declare variables on separated lines in order to reduce the occurence of git conflicts.
!!
!! List of datatypes:
!! * aim_dataset_type: the "dataset" for aim
!! * bandfft_kpt_type: the "dataset" for triple band-fft-kpt parallelization
!! * datafiles_type: gather all the variables related to files
!! * dataset_type: the "dataset" for the main abinit code
!! * MPI_type: the data related to MPI parallelization
!! * hdr_type: the header of wf, den and pot files
!! * macro_uj_type: TO BE COMPLETED
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_abitypes

 use defs_basis
 use m_abicore
 use m_distribfft

 use m_pawrhoij, only : pawrhoij_type

 implicit none
!!***

!!****t* defs_abitypes/aim_dataset_type
!! NAME
!! aim_dataset_type
!!
!! FUNCTION
!! The aim_dataset_type structured datatype
!! gathers all the input variables for the aim code
!!
!! SOURCE

 type aim_dataset_type

! Since all these input variables are described in the aim_help.html
! file, they are not described in length here ...

! Integer
  integer :: crit
  integer :: denout
  integer :: dltyp
  integer :: gpsurf
  integer :: irho
  integer :: ivol
  integer :: lapout
  integer :: nsa
  integer :: nsb
  integer :: nsc

  integer :: batom  ! Warning : corresponds to the input variable atom
  integer :: foll   ! Warning : corresponds to the input variable follow
  integer :: isurf  ! Warning : corresponds to the input variable surf
  integer :: irsur  ! Warning : corresponds to the input variable rsurf
  integer :: nph    ! Warning : corresponds to the input variable nphi
  integer :: npt    ! Warning : corresponds to the input variable inpt
  integer :: nth    ! Warning : corresponds to the input variable ntheta
  integer :: plden  ! Warning : not documented in help file ?!

  integer :: ngrid(3)

! Real
  real(dp) :: atrad
  real(dp) :: coff1
  real(dp) :: coff2
  real(dp) :: dpclim
  real(dp) :: folstp
  real(dp) :: lgrad
  real(dp) :: lgrad2
  real(dp) :: lstep
  real(dp) :: lstep2
  real(dp) :: maxatd
  real(dp) :: maxcpd
  real(dp) :: phimax
  real(dp) :: phimin

  real(dp) :: dr0    ! Warning : correspond to the input variable radstp
  real(dp) :: phi0   ! Warning : correspond to the input variable rsurdir(2)
  real(dp) :: rmin   ! Warning : correspond to the input variable ratmin
  real(dp) :: th0    ! Warning : correspond to the input variable rsurdir(1)
  real(dp) :: themax ! Warning : correspond to the input variable thetamax
  real(dp) :: themin ! Warning : correspond to the input variable thetamin

  real(dp) :: foldep(3)
  real(dp) :: scal(3)
  real(dp) :: vpts(3,4)

 end type aim_dataset_type
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/dataset_type
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
!
!! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.
!
!! Since all these input variables are described in the abinit_help.html and
!! associated html files they are not described in length here ...
!!
!! SOURCE

type dataset_type

! Integer
 integer :: iomode
 integer :: accuracy
 integer :: adpimd
 integer :: autoparal
 integer :: auxc_ixc
 integer :: awtr
 integer :: bandpp
 integer :: bdeigrf
 integer :: berryopt
 integer :: berrysav
 integer :: berrystep
 integer :: brvltt
 integer :: bs_nstates
 integer :: bs_hayd_term
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
 integer :: dmft_x2my2d
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
 integer :: eph_task
 integer :: exchn2n3d
 integer :: extrapwf
 integer :: fftgw
 integer :: fockoptmix
 integer :: frzfermi
 integer :: ga_algor
 integer :: ga_fitness
 integer :: ga_n_rules
 integer :: getcell
 integer :: getddb
 integer :: getdvdb = 0
 integer :: getddk
 integer :: getdelfd
 integer :: getdkdk
 integer :: getdkde
 integer :: getden
 integer :: getefmas
 integer :: getgam_eig2nkq
 integer :: getocc
 integer :: getpawden
 integer :: getqps
 integer :: getscr
 integer :: getsuscep
 integer :: getvel
 integer :: getwfk
 integer :: getwfkfine
 integer :: getwfq
 integer :: getxcart
 integer :: getxred
 integer :: get1den
 integer :: get1wf
 integer :: getbseig
 integer :: getbsreso
 integer :: getbscoup
 integer :: gethaydock
 integer :: goprecon
 integer :: gwcalctyp
 integer :: gwcomp
 integer :: gwgamma
 integer :: gwrpacorr
 integer :: gw_customnfreqsp
 integer :: gw_invalid_freq
 integer :: gw_qprange
 integer :: gw_nqlwl
 integer :: gw_nstep
 integer :: gw_sigxcore

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
 integer :: gw_frqim_inzgrid
 integer :: gw_frqre_inzgrid
 integer :: gw_frqre_tangrid
 integer :: gw_sctype
 integer :: gwmem
 integer :: gwpara
 integer :: hmcsst
 integer :: hmctt
 integer :: iboxcut
 integer :: icoulomb
 integer :: icutcoul
 integer :: ieig2rf
 integer :: imgmov
 integer :: imgwfstor
 integer :: inclvkb
 integer :: intxc
 integer :: ionmov
 integer :: iprcel
 integer :: iprcfc
 integer :: irandom
 integer :: irdddb
 integer :: irddvdb = 0
 integer :: irdddk
 integer :: irdden
 integer :: irdefmas
 integer :: irdhaydock
 integer :: irdpawden
 integer :: irdqps
 integer :: irdscr
 integer :: irdsuscep
 integer :: irdvdw
 integer :: irdwfk
 integer :: irdwfkfine
 integer :: irdwfq
 integer :: ird1den
 integer :: ird1wf
 integer :: irdbseig
 integer :: irdbsreso
 integer :: irdbscoup
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
 integer :: kssform
 integer :: localrdwf
 integer :: lotf_classic
 integer :: lotf_nitex
 integer :: lotf_nneigx
 integer :: lotf_version
 integer :: magconon
 integer :: maxnsym
 integer :: max_ncpus
 integer :: mband
 integer :: mep_solver
 integer :: mem_test=1
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
 integer :: nbandkss
 integer :: nbdblock
 integer :: nbdbuf
 integer :: nberry
 integer :: nc_xccc_gspace=0
 integer :: nconeq
 integer :: nctime
 integer :: ndtset
 integer :: ndynimage
 integer :: neb_algo
 integer :: nfft
 integer :: nfftdg
 integer :: nfreqim
 integer :: nfreqre
 integer :: nfreqsp
 integer :: nimage
 integer :: nkpt
 integer :: nkptgw
 integer :: nkpthf
 integer :: nline
 integer :: nnsclo
 integer :: nnsclohf
 integer :: nomegasf
 integer :: nomegasi
 integer :: nomegasrd
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
 integer :: npweps
 integer :: npwkss
 integer :: npwsigx
 integer :: npwwfn
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
 integer :: ppmodel
 integer :: prepanl
 integer :: prepgkk
 integer :: prtbbb
 integer :: prtbltztrp
 integer :: prtcif
 integer :: prtden
 integer :: prtdensph
 integer :: prtdipole
 integer :: prtdos
 integer :: prtdosm
 integer :: prtebands=1
 integer :: prtefg
 integer :: prtefmas
 integer :: prteig
 integer :: prtelf
 integer :: prtfc
 integer :: prtfull1wf
 integer :: prtfsurf
 integer :: prtgsr=1
 integer :: prtgden
 integer :: prtgeo
 integer :: prtgkk
 integer :: prtkden
 integer :: prtkpt
 integer :: prtlden
 integer :: prtnabla
 integer :: prtnest
 integer :: prtpmp
 integer :: prtposcar
 integer :: prtphdos
 integer :: prtphbands=1
 integer :: prtphsurf=0
 integer :: prtpot
 integer :: prtpsps=0
 integer :: prtspcur
 integer :: prtstm
 integer :: prtsuscep
 integer :: prtvclmb
 integer :: prtvdw
 integer :: prtvha
 integer :: prtvhxc
 integer :: prtkbff=0
 integer :: prtvol
 integer :: prtvolimg
 integer :: prtvpsp
 integer :: prtvxc
 integer :: prtwant
 integer :: prtwf
 integer :: prtwf_full
 integer :: prtxml
 integer :: prt1dm
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
 integer :: spmeth
 integer :: string_algo
 integer :: symmorphi
 integer :: symchi
 integer :: symsigma
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
 integer :: useria=0
 integer :: userib=0
 integer :: useric=0
 integer :: userid=0
 integer :: userie=0
 integer :: usewvl
 integer :: usexcnhat_orig
 integer :: useylm
 integer :: use_slk
 integer :: use_yaml
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
 integer :: supercell_latt(3,3)
 integer :: ucrpa_bands(2)
 integer :: vdw_supercell(3)
 integer :: vdw_typfrag(100)
 integer :: wvl_ngauss(2)

!Integer allocatables
 integer, allocatable ::  algalch(:)    ! algalch(ntypalch)
 integer, allocatable ::  bdgw(:,:,:)   ! bdgw(2,nkptgw,nsppol)
 integer, allocatable ::  dynimage(:)   ! dynimage(nimage or mxnimage)
 integer, allocatable ::  efmas_bands(:,:) ! efmas_bands(2,nkptgw)
 integer, allocatable ::  iatfix(:,:)   ! iatfix(3,natom)
 integer, allocatable ::  iatsph(:)     ! iatsph(natsph)
 integer, allocatable ::  istwfk(:)     ! istwfk(nkpt)
 integer, allocatable ::  kberry(:,:)   ! kberry(3,nberry)
 integer, allocatable ::  lexexch(:)    ! lexexch(ntypat)
 integer, allocatable ::  ldaminushalf(:) !lminushalf(ntypat)
 integer, allocatable ::  lpawu(:)      ! lpawu(ntypat)
 integer, allocatable ::  nband(:)      ! nband(nkpt*nsppol)
 integer, allocatable ::  plowan_iatom(:)    ! plowan_iatom(plowan_natom)
 integer, allocatable ::  plowan_it(:)     ! plowan_it(plowan_nt*3)
 integer, allocatable ::  plowan_lcalc(:)    ! plowan_lcalc(\sum_iatom plowan_nbl)
 integer, allocatable ::  plowan_nbl(:)     ! plowan_nbl(plowan_natom)
 integer, allocatable ::  plowan_projcalc(:) ! plowan_projcalc(\sum_iatom plowan_nbl)
 integer, allocatable ::  prtatlist(:)  ! prtatlist(natom)
 integer, allocatable ::  so_psp(:)     ! so_psp(npsp)
 integer, allocatable ::  symafm(:)     ! symafm(nsym)
 integer, allocatable ::  symrel(:,:,:) ! symrel(3,3,nsym)
 integer, allocatable ::  typat(:)      ! typat(natom)

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
 real(dp) :: fermie_nest
 real(dp) :: focktoldfe
 real(dp) :: freqim_alpha
 real(dp) :: freqremin
 real(dp) :: freqremax
 real(dp) :: freqspmin
 real(dp) :: freqspmax
 real(dp) :: friction
 real(dp) :: fxcartfactor
 real(dp) :: ga_opt_percent
 real(dp) :: gwencomp
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
 real(dp) :: mbpt_sciss
 real(dp) :: mdf_epsinf
 real(dp) :: mdwall
 real(dp) :: mep_mxstep
 real(dp) :: nelect
 real(dp) :: noseinert
 real(dp) :: omegasimax
 real(dp) :: omegasrdmax
 real(dp) :: pawecutdg
 real(dp) :: pawovlp
 real(dp) :: pawujrad
 real(dp) :: pawujv
 real(dp) :: posocc
 real(dp) :: postoldfe
 real(dp) :: postoldff
 real(dp) :: ppmfrq
 real(dp) :: pw_unbal_thresh
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
 real(dp) :: userra=zero
 real(dp) :: userrb=zero
 real(dp) :: userrc=zero
 real(dp) :: userrd=zero
 real(dp) :: userre=zero
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
 real(dp) :: vdw_df_threshold
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
 integer :: bs_algorithm
 integer :: bs_haydock_niter
 integer :: bs_exchange_term
 integer :: bs_coulomb_term
 integer :: bs_calctype
 integer :: bs_coupling
 integer :: bs_interp_mode
 integer :: bs_interp_prep
 integer :: bs_interp_method
 integer :: bs_interp_rl_nb

 real(dp) :: bs_interp_m3_width

 integer  :: bs_interp_kmult(3)
 real(dp) :: bs_haydock_tol(2)

 integer,allocatable :: bs_loband(:)

 real(dp) :: bs_eh_cutoff(2)
 real(dp) :: bs_freq_mesh(3)

!END VARIABLES FOR @Bethe-Salpeter.

 integer :: gpu_linalg_limit ! MT sept2012: had to add this keyword at the end
                             ! to get BigDFT automatic tests work on shiva and littlebuda (why????)

!EPH variables
! ifc variables
 integer :: asr
 integer :: dipdip
 integer :: chneut
 integer :: symdynmat

! Phonon variables.
 integer :: ph_freez_disp_addStrain
 integer :: ph_freez_disp_option
 integer :: ph_freez_disp_nampl
 integer :: ph_ndivsm    ! =20
 integer :: ph_nqpath    !=0
 integer :: ph_ngqpt(3)  !0
 integer :: ph_nqshift

 real(dp),allocatable :: ph_freez_disp_ampl(:,:)
  ! ph_freez_disp_ampl(5,ph_freez_disp_nampl)
 real(dp),allocatable :: ph_qshift(:,:)
  ! ph_qshift(3, ph_nqshift)
 real(dp),allocatable :: ph_qpath(:,:)
  ! ph_qpath(3, nqpath)

! e-ph variables
 real(dp) :: eph_mustar
 integer :: eph_intmeth ! = 1
 real(dp) :: eph_extrael != zero
 real(dp) :: eph_fermie != huge(one)
 integer :: eph_frohlichm != 0
 real(dp) :: eph_fsmear != 0.01
 real(dp) :: eph_fsewin != 0.04
 real(dp) :: eph_tols_idelta(2) = [tol12, tol12]
 integer :: eph_phrange(2) = 0

 integer :: eph_ngqpt_fine(3)
 integer :: eph_np_pqbks(5) = 0

 integer :: eph_stern = 0
 integer :: eph_transport

 integer :: ph_intmeth
 integer :: prteliash = 0
 real(dp) :: ph_wstep
 real(dp) :: ph_smear
 integer :: ddb_ngqpt(3)
 real(dp) :: ddb_shiftq(3)

 integer :: mixprec = 0
 integer :: symv1scf = 0
 integer :: dvdb_add_lr = 1

 integer :: sigma_bsum_range(2) = 0

 real(dp) :: sigma_erange(2) = -one

 integer :: sigma_ngkpt(3) = 0
 ! K-mesh for Sigma_{nk} (only IBZ points). Alternative to kptgw.

 integer :: sigma_nshiftk = 1
 ! Number of shifts in k-mesh for Sigma_{nk}.

 real(dp) :: frohl_params(4) = zero

 real(dp),allocatable :: sigma_shiftk(:,:)
 ! sigma_shiftk(3, sigma_nshiftk)
 ! shifts in k-mesh for Sigma_{nk}.
!END EPH

 integer :: ndivsm=0
 integer :: nkpath=0
 real(dp) :: einterp(4)=zero
 real(dp),allocatable :: kptbounds(:,:)
 real(dp) :: tmesh(3) ! = [5._dp, 59._dp, 6._dp] This triggers a bug in the bindings

 character(len=fnlen) :: getkerange_path = ABI_NOFILE
 !character(len=fnlen) :: getpot_path = ABI_NOFILE
 !character(len=fnlen) :: getsigeph_path = ABI_NOFILE

 end type dataset_type
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/MPI_type
!! NAME
!! MPI_type
!!
!! FUNCTION
!! The MPI_type structured datatype gather different information
!! about the MPI parallelisation: number of processors,
!! the index of my processor, the different groups of processors, etc ...
!!
!! SOURCE

 type MPI_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.
! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.

! *****************************************************************************************
! Please make sure that initmpi_seq is changed so that any variable or any flag in MPI_type
! is initialized with the value used for sequential executions.
! In particular any MPI communicator should be set to MPI_COMM_SELF
! ************************************************************************************

  ! Set of variables for parallelism, that do NOT depend on input variables.
  ! These are defined for each dataset

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main variables for parallelisation
  integer :: comm_world
  ! number of the world communicator MPI COMM WORLD

  integer :: me
  ! number of my processor in the group of all processors

  integer :: nproc
  ! number of processors

  integer :: me_g0
  ! if set to 1, means that the current processor is taking care of the G(0 0 0) planewave.

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over atoms (PAW)
   integer :: comm_atom
   ! Communicator over atoms

   integer :: nproc_atom
   ! Size of the communicator over atoms

   integer :: my_natom
   ! Number of atoms treated by current proc

   integer,pointer :: my_atmtab(:) => null()
   ! Indexes of the atoms treated by current processor

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over perturbations
   integer :: paral_pert
   ! to activate parallelisation over perturbations for linear response

   integer :: comm_pert
   ! communicator for calculating perturbations

   integer :: comm_cell_pert
   ! general communicator over all processors treating the same cell

   integer :: me_pert
   ! number of my processor in my group of perturbations

   integer :: nproc_pert
   ! number of processors in my group of perturbations

   integer, allocatable :: distrb_pert(:)
   ! distrb_pert(1:npert)
   ! index of processor treating each perturbation

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over images
   integer :: paral_img
   ! Flag activated if parallelization over image is on

   integer :: my_nimage
   ! Number of images of the cell treated by current proc (i.e. local nimage)

   integer :: comm_img
   ! Communicator over all images

   integer :: me_img
   ! Index of my processor in the comm. over all images

   integer :: nproc_img
   ! Size of the communicator over all images

   integer,allocatable :: distrb_img(:)
   ! distrb_img(1:dtset%nimage)
   ! index of processor treating each image (in comm_img communicator)

   integer,allocatable :: my_imgtab(:)
   ! index_img(1:my_nimage) indexes of images treated by current proc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over the cell
   integer :: comm_cell
   ! local Communicator over all processors treating the same cell

   integer :: me_cell
   ! Index of my processor in the comm. over one cell

   integer :: nproc_cell
   ! Size of the communicator over one cell

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over fft
   integer :: comm_fft
   ! Communicator over fft

   integer :: me_fft
   ! Rank of my processor in my group of FFT

   integer :: nproc_fft
   ! number of processors in my group of FFT

   type(distribfft_type),pointer :: distribfft  => null()
   ! Contains all the information related to the FFT distribution

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over bands
   integer :: paralbd
    ! paralbd=0 : (no //ization on bands)
    ! paralbd=1 : (//ization on bands)

   integer :: comm_band
   ! Communicator over bands

   integer :: me_band
   ! Rank of my proc in my group of bands

   integer :: nproc_band
   ! Number of procs on which we distribute bands

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the spinor parallelisation
   integer :: paral_spinor
   ! Flag: activation of parallelization over spinors

   integer :: comm_spinor
   ! Communicator over spinors

   integer :: me_spinor
   ! Rank of my proc in the communicator over spinors
   ! Note: me_spinor is related to the index treated by current proc
   ! (nspinor_index= mpi_enreg%me_spinor + 1)

   integer :: nproc_spinor
   ! Number of procs on which we distribute spinors

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the kpt/nsppol parallelisation
   integer :: comm_kpt
   ! Communicator over kpt

   integer :: me_kpt
   ! Rank of my proc in the communicator over kpt

   integer :: nproc_kpt
   ! Number of procs on which we distribute kpt

   integer, allocatable :: proc_distrb(:,:,:)
    ! proc_distrb(nkpt,mband,nsppol)
    ! number of the processor that will treat
    ! each band in each k point.

   integer :: my_isppoltab(2)
    ! my_isppoltab(2) contains the flags telling which value of isppol is treated by current proc
    ! in sequential, its value is (1,0) when nsppol=1 and (1,1) when nsppol=2
    ! in parallel,   its value is (1,0) when nsppol=1 and (1,0) when nsppol=2 and up-spin is treated
    !                                                  or (0,1) when nsppol=2 and dn-spin is treated
    !                                                  or (1,1) when nsppol=2 and both spins are treated

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the band-FFT-kpt-spinor parallelisation
   integer :: paral_kgb
   ! Flag: activation of parallelization over kpt/band/fft

   integer :: bandpp
   ! # of Bands Per Processor

   integer :: comm_bandspinorfft
   ! Cartesian communicator over band-fft-spinor

   integer :: comm_bandfft
   ! Cartesian communicator over the band-fft

   integer :: comm_kptband
   ! Communicator over kpt-band subspace

   integer :: comm_spinorfft
   ! Communicator over fft-spinors subspace

   integer :: comm_bandspinor
   ! Communicator over band-spinors subspace

   integer, allocatable :: my_kgtab(:,:)
    ! (mpw, mkmem)
    ! Indexes of kg treated by current proc
    ! i.e. mapping betwee the G-vector stored by this proc and the list of G-vectors
    ! one would have in the sequential version. See kpgsph in m_fftcore.

   integer, allocatable :: my_kpttab(:)
    ! Indicates the correspondence between the ikpt and ikpt_this_proc

   real(dp) :: pw_unbal_thresh
    !Threshold (in %) activating the plane-wave load balancing process (see kpgsph routine)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over kpt/nsppol in the Berry Phase case
   integer, allocatable :: kptdstrb(:,:,:)
    ! kptdstrb(me,ineigh,ikptloc)
    ! tab of processors required for dfptnl_mv.f and berryphase_new.f

   integer, allocatable :: kpt_loc2fbz_sp(:,:,:)
    ! kpt_loc2fbz_sp(nproc, dtefield%fmkmem_max, 2)
    ! K-PoinT LOCal TO Full Brilloin Zone and Spin Polarization
    ! given a processor and the local number of the k-point for this proc,
    ! give the number of the k-point in the FBZ and the isppol;
    ! necessary for synchronisation in berryphase_new
    ! kpt_loc2fbz(iproc, ikpt_loc,1) = ikpt
    ! kpt_loc2fbz(iproc, ikpt_loc,2) = isppol

   integer, allocatable :: kpt_loc2ibz_sp(:,:,:)

   ! TODO: Is is still used?
   integer, allocatable :: mkmem(:)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for Hartree-Fock's parallelisation
   integer :: paral_hf
   ! Flag: activation of parallelization for Hartree-Fock

   integer :: comm_hf
   ! Communicator over the k-points and bands of occupied states for Hartree-Fock

   integer :: me_hf
   ! Rank of my proc in the communicator for Hartree-Fock

   integer :: nproc_hf
   ! Number of procs on which we distribute the occupied states for Hartree-Fock

   integer, allocatable :: distrb_hf(:,:,:)
    ! distrb_hf(nkpthf,nbandhf,1)
    ! index of processor treating each occupied states for Hartree Fock.
    ! No spin-dependence because only the correct spin is treated (in parallel) or both spins are considered (sequential)
    ! but we keep the third dimension (always equal to one) to be able to use the same routines as the one for proc_distrb

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the wavelet parallelisation
   integer :: comm_wvl
   ! Communicator over real space grids for WVLs

   integer :: me_wvl
   ! Rank of my proc for WVLs

   integer :: nproc_wvl
   ! Number of procs for WVLs
   ! Array to store the description of the scaterring in real space of
   ! the potentials and density. It is allocated to (0:nproc-1,4).
   ! The four values are:
   ! - the density size in z direction ( = ngfft(3)) ;
   ! - the potential size in z direction ( <= ngfft(3)) ;
   ! - the position of the first value in the complete array ;
   ! - the shift for the potential in the array.
   ! This array is a pointer on a BigDFT handled one.

   integer, pointer :: nscatterarr(:,:) => null()
   ! Array to store the total size (of this proc) of the potentails arrays when
   ! the memory is distributed following nscatterarr.
   ! This array is a pointer on a BigDFT handled one.

   integer, pointer :: ngatherarr(:,:) => null()
   ! Store the ionic potential size in z direction.
   ! This array is a pointer on a BigDFT handled one.

   integer :: ngfft3_ionic
   ! End wavelet additions

 end type MPI_type
!!***

!----------------------------------------------------------------------

!!****t* defs_datatypes/datafiles_type
!! NAME
!! datafiles_type
!!
!! FUNCTION
!! The datafiles_type structures datatype gather all the variables
!! related to files, such as filename, and file units.
!! For one dataset, it is initialized in 95_drive/dtfil_init1.F90,
!! and will not change at all during the treatment of the dataset.
!!
!! SOURCE

 type datafiles_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! These keywords are only used in algorithms using images of the cell
  integer :: getwfk_from_image
   ! index of image from which read WFK file (0 if standard WFK)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getden_from_image
   ! index of image from which read DEN file (0 if standard DEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getpawden_from_image
   ! index of image from which read PAWDEN file (0 if standard PAWDEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: ireadddb
   ! ireadddb non-zero  if the ddb file must be read

  integer :: ireadden
   ! ireadden non-zero  if the den file must be read

  integer :: ireadkden
   ! ireadkden non-zero  if the kden file must be read

  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffk file must be read
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1 file must be read
   !   (if ird1wf non-zero or get1wf non-zero)

  integer :: unchi0  ! unit number for chi0 files
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: unddk   ! unit number for ddk 1WF file
  integer :: undkdk  ! unit number for 2WF file (dkdk)
  integer :: undkde  ! unit number for 2WF file (dkde)
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unkss   ! unit number for KSS file
  integer :: unqps   ! unit number for QPS file
  integer :: unscr   ! unit number for SCR file
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwff3  ! unit number for wavefunctions, number three
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwft3  ! unit number for wavefunctions, temporary three
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: unpaw   ! unit number for temporary PAW data (for ex. rhoij_nk) (Paw only)
  integer :: unpaw1  ! unit number for temporary PAW first-order cprj1=<c1_k,q|p>(1) data
  integer :: unpawq  ! unit number for temporary PAW cprjq=<c+_k+q|p> at k+qdata
  integer :: unpos   ! unit number for restart molecular dynamics

  character(len=fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input :
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character(len=fnlen) :: filddbsin
   ! if no dataset mode             : abi//'DDB'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DDB'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DDB'

  character(len=fnlen) :: fildensin
   ! if no dataset mode             : abi//'DEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character(len=fnlen) :: fildvdbin
   ! if no dataset mode             : abi//'DVDB'
   ! if dataset mode, and getdvdb==0 : abi//'_DS'//trim(jdtset)//'DVDB'
   ! if dataset mode, and getdvdb/=0 : abo//'_DS'//trim(jgetden)//'DVDB'

  character(len=fnlen) :: filkdensin
   ! if no dataset mode             : abi//'KDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'KDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'KDEN'

  character(len=fnlen) :: filpawdensin
   ! if no dataset mode             : abi//'PAWDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'PAWDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'PAWDEN'

! character(len=fnlen) :: filpsp(ntypat)
   ! the filenames of the pseudopotential files, from the standard input.

  character(len=fnlen) :: filstat
   ! tmp//'_STATUS'

  character(len=fnlen) :: fnamewffk
   ! the name of the ground-state wavefunction file to be read (see driver.F90)

  character(len=fnlen) :: fnamewffq
   ! the name of the k+q ground-state wavefunction file to be read (see driver.F90)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdelfd
   ! the generic name of the electric field response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkdk
   ! the generic name of the 2nd order dkdk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkde
   ! the generic name of the 2nd order dkde response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewff1
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fildens1in   ! to be described by MVeithen

  character(len=fnlen) :: fname_tdwf

  character(len=fnlen) :: fname_w90

  character(len=fnlen) :: fnametmp_wf1
  character(len=fnlen) :: fnametmp_wf2
  character(len=fnlen) :: fnametmp_1wf1
  character(len=fnlen) :: fnametmp_1wf2
  character(len=fnlen) :: fnametmp_wfgs
  character(len=fnlen) :: fnametmp_wfkq
   ! Set of filemanes formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _WF1, _WF2 ...
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_kg
  character(len=fnlen) :: fnametmp_kgq
  character(len=fnlen) :: fnametmp_kg1
  character(len=fnlen) :: fnametmp_dum
  character(len=fnlen) :: fnametmp_ylm
  character(len=fnlen) :: fnametmp_ylm1
  character(len=fnlen) :: fnametmp_paw
  character(len=fnlen) :: fnametmp_paw1
  character(len=fnlen) :: fnametmp_pawq
   ! Set of filemanes formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _KG, _DUM, followed
   ! by the index of the processor.
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_cg
  character(len=fnlen) :: fnametmp_cprj
  character(len=fnlen) :: fnametmp_eig
  character(len=fnlen) :: fnametmp_1wf1_eig
  character(len=fnlen) :: fnametmp_fft
  character(len=fnlen) :: fnametmp_kgs
  character(len=fnlen) :: fnametmp_sustr
  character(len=fnlen) :: fnametmp_tdexcit
  character(len=fnlen) :: fnametmp_tdwf

!@Bethe-Salpeter
! New files introduced for the Bethe-Salpeter part.

   character(len=fnlen) :: fnameabi_bsham_reso
    ! if no dataset mode             : abi//'BSR'
    ! if dataset mode, and getbsreso==0 : abi//'_DS'//trim(jdtset)//'BSR'
    ! if dataset mode, and getbsreso/=0 : abo//'_DS'//trim(jget_reso_bsham)//'BSR'

   character(len=fnlen) :: fnameabi_bsham_coup
    ! if no dataset mode             : abi//'BSC'
    ! if dataset mode, and getbscoup==0 : abi//'_DS'//trim(jdtset)//'BSC'
    ! if dataset mode, and getbscoup/=0 : abo//'_DS'//trim(jget_coup_bsham)//'BSC'

  character(len=fnlen) :: fnameabi_bseig
   ! The name of the file containing the eigenstates and eigenvalues of the Bethe-Salpeter Hamiltonian
   ! if no dataset mode             : abi//'BS_EIG'
   ! if dataset mode, and getbseig==0 : abi//'_DS'//trim(jdtset)//'BS_EIG'
   ! if dataset mode, and getbseig/=0 : abo//'_DS'//trim(jget_bseig)//'BS_EIG'

   character(len=fnlen) :: fnameabi_haydock
   ! The prefix used to construct the names of the files containing the coefficients of the
   ! continued fractions produced by the Haydock iterative algorithm.
   ! if no dataset mode             : abi//'HAYDOCK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'HAYDOCK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'HAYDOCK'

   character(len=fnlen) :: fnameabi_wfkfine
   ! The name of the file containing the wavefunctions on a fine grid
   ! if no dataset mode             : abi//'WFK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'WFK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'WFK'

!END @BEthe-Salpeter

!The following filenames do not depend on itimimage, iimage and itime loops.
!Note the following convention:
!  fnameabo_* are filenames used for ouput results (writing)
!  fnameabi_* are filenames used for data that should be read by the code.
!  fnametmp_* are filenames used for temporary files that should be erased at the end of each dataset.
!
!If a file does not have the corresponding "abi" or the corresponding "abo" name, that means that
!that particular file is only used for writing or for reading results, respectively.

  character(len=fnlen) :: fnameabi_efmas
  character(len=fnlen) :: fnameabi_hes
  character(len=fnlen) :: fnameabi_phfrq
  character(len=fnlen) :: fnameabi_phvec
  character(len=fnlen) :: fnameabi_qps
  character(len=fnlen) :: fnameabi_scr            ! SCReening file (symmetrized inverse dielectric matrix)
  character(len=fnlen) :: fnameabi_sus            ! KS independent-particle polarizability file
  character(len=fnlen) :: fnameabo_ddb
  character(len=fnlen) :: fnameabo_den
  character(len=fnlen) :: fnameabo_dos
  character(len=fnlen) :: fnameabo_dvdb
  character(len=fnlen) :: fnameabo_eelf
  character(len=fnlen) :: fnameabo_eig
  character(len=fnlen) :: fnameabo_eigi2d
  character(len=fnlen) :: fnameabo_eigr2d
  character(len=fnlen) :: fnameabo_em1
  character(len=fnlen) :: fnameabo_em1_lf
  character(len=fnlen) :: fnameabo_em1_nlf
  character(len=fnlen) :: fnameabo_fan
  character(len=fnlen) :: fnameabo_gkk
  character(len=fnlen) :: fnameabo_gw
  character(len=fnlen) :: fnameabo_gw_nlf_mdf
  character(len=fnlen) :: fnameabo_kss
  character(len=fnlen) :: fnameabo_moldyn
  character(len=fnlen) :: fnameabo_pot
  character(len=fnlen) :: fnameabo_qps            ! Quasi-Particle band structure file.
  character(len=fnlen) :: fnameabo_qp_den
  character(len=fnlen) :: fnameabo_qp_pawden      ! Full QP density
  character(len=fnlen) :: fnameabo_qp_dos
  character(len=fnlen) :: fnameabo_qp_eig
  character(len=fnlen) :: fnameabo_rpa
  character(len=fnlen) :: fnameabo_rpa_nlf_mdf
  character(len=fnlen) :: fnameabo_scr
  character(len=fnlen) :: fnameabo_sgm
  character(len=fnlen) :: fnameabo_sgr
  character(len=fnlen) :: fnameabo_sig
  character(len=fnlen) :: fnameabo_spcur
  character(len=fnlen) :: fnameabo_sus
  character(len=fnlen) :: fnameabo_vha
  character(len=fnlen) :: fnameabo_vpsp
  character(len=fnlen) :: fnameabo_vso
  character(len=fnlen) :: fnameabo_vxc
  character(len=fnlen) :: fnameabo_wan
  character(len=fnlen) :: fnameabo_wfk
  character(len=fnlen) :: fnameabo_wfq
  character(len=fnlen) :: fnameabo_w90
  character(len=fnlen) :: fnameabo_1wf
  character(len=fnlen) :: fnameabo_gwdiag
  character(len=fnlen) :: fnameabo_nlcc_derivs
  character(len=fnlen) :: fnameabo_pspdata

!The following filenames are initialized only iniside itimimage, iimage and itime loops,
!and are appended with the adequate specifier 'app'.

  character(len=fnlen) :: fnameabo_app
  character(len=fnlen) :: fnameabo_app_atmden_core
  character(len=fnlen) :: fnameabo_app_atmden_full
  character(len=fnlen) :: fnameabo_app_atmden_val
  character(len=fnlen) :: fnameabo_app_n_tilde
  character(len=fnlen) :: fnameabo_app_n_one
  character(len=fnlen) :: fnameabo_app_nt_one
  character(len=fnlen) :: fnameabo_app_bxsf
  character(len=fnlen) :: fnameabo_app_cif
  character(len=fnlen) :: fnameabo_app_den
  character(len=fnlen) :: fnameabo_app_dos
  character(len=fnlen) :: fnameabo_app_elf
  character(len=fnlen) :: fnameabo_app_elf_down
  character(len=fnlen) :: fnameabo_app_elf_up
  character(len=fnlen) :: fnameabo_app_eig
  character(len=fnlen) :: fnameabo_app_fatbands
  character(len=fnlen) :: fnameabo_app_gden1
  character(len=fnlen) :: fnameabo_app_gden2
  character(len=fnlen) :: fnameabo_app_gden3
  character(len=fnlen) :: fnameabo_app_geo
  character(len=fnlen) :: fnameabo_app_kden
  character(len=fnlen) :: fnameabo_app_lden
  character(len=fnlen) :: fnameabo_app_nesting
  character(len=fnlen) :: fnameabo_app_pawden
  character(len=fnlen) :: fnameabo_app_pot
  character(len=fnlen) :: fnameabo_app_opt
  character(len=fnlen) :: fnameabo_app_opt2
  character(len=fnlen) :: fnameabo_app_stm
  character(len=fnlen) :: fnameabo_app_vclmb
  character(len=fnlen) :: fnameabo_app_vha
  character(len=fnlen) :: fnameabo_app_vhxc
  character(len=fnlen) :: fnameabo_app_vhpsp
  character(len=fnlen) :: fnameabo_app_vpsp
  character(len=fnlen) :: fnameabo_app_vxc
  character(len=fnlen) :: fnameabo_app_wfk
  character(len=fnlen) :: fnameabo_app_1dm
  character(len=fnlen) :: fnameabo_app_vha_1dm
  character(len=fnlen) :: fnameabo_app_vclmb_1dm
  character(len=fnlen) :: fnametmp_app_den
  character(len=fnlen) :: fnametmp_app_kden

 end type datafiles_type
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/hdr_type
!! NAME
!! hdr_type
!!
!! FUNCTION
!! It contains all the information needed to write a header for a wf, den or pot file.
!! The structure of the header is explained in the abinit_help.html and other associated html files.
!! The datatype is considered as an object, to which are attached a whole
!! set of "methods", actually, different subroutines.
!! A few of these subroutines are: hdr_init, hdr_update, hdr_free, hdr_check, hdr_io, hdr_skip.
!!
!! SOURCE

 type hdr_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
  integer :: date          ! starting date
  integer :: headform      ! format of the header
  integer :: intxc         ! input variable
  integer :: ixc           ! input variable
  integer :: mband         ! maxval(hdr%nband)
  integer :: natom         ! input variable
  integer :: nkpt          ! input variable
  integer :: npsp          ! input variable
  integer :: nspden        ! input variable
  integer :: nspinor       ! input variable
  integer :: nsppol        ! input variable
  integer :: nsym          ! input variable
  integer :: ntypat        ! input variable
  integer :: occopt        ! input variable
  integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
  integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
  integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)

  integer :: kptopt          ! input variable (defines symmetries used for k-point sampling)
  integer :: pawcpxocc       ! input variable
  integer :: nshiftk_orig=1  ! original number of shifts given in input (changed in inkpts, the actual value is nshiftk)
  integer :: nshiftk=1       ! number of shifts after inkpts.
  integer :: icoulomb        ! input variable.

  real(dp) :: ecut         ! input variable
  real(dp) :: ecutdg       ! input variable (ecut for NC psps, pawecutdg for paw)
  real(dp) :: ecutsm       ! input variable
  real(dp) :: ecut_eff     ! ecut*dilatmx**2 (dilatmx is an input variable)
  real(dp) :: etot         ! EVOLVING variable
  real(dp) :: fermie       ! EVOLVING variable
  real(dp) :: residm       ! EVOLVING variable
  real(dp) :: stmbias      ! input variable
  real(dp) :: tphysel      ! input variable
  real(dp) :: tsmear       ! input variable
  real(dp) :: nelect       ! number of electrons (computed from pseudos and charge)
  real(dp) :: charge       ! input variable

! This record is not a part of the hdr_type, although it is present in the
! header of the files. This is because it depends on the kind of file
! that is written, while all other information does not depend on it.
! It was preferred to let it be initialized or defined outside of hdr_type.
! integer :: fform         ! file format

  real(dp) :: qptn(3)
 ! the wavevector, in case of a perturbation

  real(dp) :: rprimd(3,3)
  ! EVOLVING variables

  integer :: ngfft(3)
  ! input variable

  integer :: nwvlarr(2)
  ! nwvlarr(2) array holding the number of wavelets for each resolution.

  integer :: kptrlatt_orig(3,3)
  ! Original kptrlatt

  integer :: kptrlatt(3,3)
  ! kptrlatt after inkpts.

  integer, allocatable :: istwfk(:)
  ! input variable istwfk(nkpt)

  integer, allocatable :: lmn_size(:)
  ! lmn_size(npsp) from psps

  integer, allocatable :: nband(:)
  ! input variable nband(nkpt*nsppol)

  integer, allocatable :: npwarr(:)
  ! npwarr(nkpt) array holding npw for each k point

  integer, allocatable :: pspcod(:)
  ! pscod(npsp) from psps

  integer, allocatable :: pspdat(:)
  ! psdat(npsp) from psps

  integer, allocatable :: pspso(:)
  ! pspso(npsp) from psps

  integer, allocatable :: pspxc(:)
  ! pspxc(npsp) from psps

  integer, allocatable :: so_psp(:)
  ! input variable so_psp(npsp)

  integer, allocatable :: symafm(:)
  ! input variable symafm(nsym)

  integer, allocatable :: symrel(:,:,:)
  ! input variable symrel(3,3,nsym)

  integer, allocatable :: typat(:)
  ! input variable typat(natom)

  real(dp), allocatable :: kptns(:,:)
  ! input variable kptns(3,nkpt)

  real(dp), allocatable :: occ(:)
  ! EVOLVING variable occ(bantot)

  real(dp), allocatable :: tnons(:,:)
  ! input variable tnons(3,nsym)

  real(dp), allocatable :: wtk(:)
  ! weight of kpoints wtk(nkpt)

  real(dp),allocatable :: shiftk_orig(:,:)
  ! original shifts given in input (changed in inkpts).

  real(dp),allocatable :: shiftk(:,:)
  ! shiftk(3,nshiftk), shiftks after inkpts

  real(dp),allocatable :: amu(:)
  ! amu(ntypat) ! EVOLVING variable

  real(dp), allocatable :: xred(:,:)
  ! EVOLVING variable xred(3,natom)

  real(dp), allocatable :: zionpsp(:)
  ! zionpsp(npsp) from psps

  real(dp), allocatable :: znuclpsp(:)
  ! znuclpsp(npsp) from psps
  ! Note the difference between (znucl|znucltypat) and znuclpsp !

  real(dp), allocatable :: znucltypat(:)
  ! znucltypat(ntypat) from alchemy

  character(len=6) :: codvsn
  ! version of the code

  character(len=132), allocatable :: title(:)
  ! title(npsp) from psps

  character(len=md5_slen),allocatable :: md5_pseudos(:)
  ! md5pseudos(npsp)
  ! md5 checksums associated to pseudos (read from file)

  ! EVOLVING variable, only for paw
  type(pawrhoij_type), allocatable :: pawrhoij(:)

 end type hdr_type
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/ab_dimensions
!! NAME
!! ab_dimensions
!!
!! FUNCTION
!! One record for each dimension of arrays used in ABINIT.
!! Will be used to e.g.:
!! - contain the maximum size attained over all datasets (mxvals)
!! - indicate whether this dimension is the same for all datasets or not (multivals).
!! Used for example inside outvars
!!
!! SOURCE

 type ab_dimensions

    integer :: ga_n_rules   ! maximal value of input ga_n_rules for all the datasets
    integer :: gw_nqlwl     ! maximal value of input gw_nqlwl for all the datasets
    integer :: lpawu        ! maximal value of input lpawu for all the datasets
    integer :: mband
    integer :: mband_upper ! maximal value of input nband for all the datasets
                           ! Maybe this one could be removed
    integer :: natom
    integer :: natpawu     ! maximal value of number of atoms on which +U is applied for all the datasets
    integer :: natsph      ! maximal value of input natsph for all the datasets
    integer :: natsph_extra  ! maximal value of input natsph_extra for all the datasets
    integer :: natvshift   ! maximal value of input natvshift for all the datasets
    integer :: nberry = 20 ! This is presently a fixed value. Should be changed.
    integer :: nbandhf
    integer :: nconeq      ! maximal value of input nconeq for all the datasets
    integer :: n_efmas_dirs
    integer :: nfreqsp
    integer :: n_projection_frequencies
    integer :: nimage
    integer :: nimfrqs
    integer :: nkpt       ! maximal value of input nkpt for all the datasets
    integer :: nkptgw     ! maximal value of input nkptgw for all the datasets
    integer :: nkpthf     ! maximal value of input nkpthf for all the datasets
    integer :: nnos       ! maximal value of input nnos for all the datasets
    integer :: nqptdm     ! maximal value of input nqptdm for all the datasets
    integer :: nshiftk
    integer :: nsp
    integer :: nspinor    ! maximal value of input nspinor for all the datasets
    integer :: nsppol     ! maximal value of input nsppol for all the datasets
    integer :: nsym       ! maximum number of symmetries
    integer :: ntypalch
    integer :: ntypat     ! maximum number of types of atoms
    integer :: nzchempot  ! maximal value of input nzchempot for all the datasets

 end type ab_dimensions
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/macro_uj_type
!! NAME
!! dtmacro_uj
!!
!! FUNCTION
!! This data type contains the potential shifts and the occupations
!! for the determination of U and J for the DFT+U calculations.
!! iuj=1,2: non-selfconsistent calculations. iuj=3,4 selfconsistent calculations.
!! iuj=2,4  => pawujsh<0 ; iuj=1,3 => pawujsh >0,
!!
!! SOURCE

 type macro_uj_type

! Integer
  integer :: iuj        ! dataset treated
  integer :: nat        ! number of atoms U (J) is determined on
  integer :: ndtset     ! total number of datasets
  integer :: nspden     ! number of densities treated
  integer :: macro_uj   ! which mode the determination runs in
  integer :: pawujat    ! which atom U (J) is determined on
  integer :: pawprtvol  ! controlling amount of output
  integer :: option     ! controls the determination of U (1 with compensating charge bath, 2 without)
  integer :: dmatpuopt  ! controls the renormalisation of the PAW projectors

! Real
  real(dp) :: diemix    ! mixing parameter
  real(dp) :: mdist     ! maximal distance of ions
  real(dp) :: pawujga   ! gamma for inversion of singular matrices
  real(dp) :: ph0phiint ! integral of phi(:,1)*phi(:,1)
  real(dp) :: pawujrad  ! radius to which U should be extrapolated.
  real(dp) :: pawrad    ! radius of the paw atomic data

! Integer arrays
  integer , allocatable  :: scdim(:)
  ! size of supercell

! Real arrays
  real(dp) , allocatable :: occ(:,:)
  ! occupancies after a potential shift: occ(ispden,nat)

  real(dp) , allocatable :: rprimd(:,:)
  ! unit cell for symmetrization

  real(dp) , allocatable :: vsh(:,:)
  ! potential shifts on atoms, dimensions: nspden,nat

  real(dp) , allocatable :: xred(:,:)
  ! atomic position for symmetrization

  real(dp) , allocatable :: wfchr(:)
  ! wfchr(1:3): zion, n and l of atom on which projection is done
  ! wfchr(4:6): coefficients ai of a0+a1*r+a2*r^2, fit to the wfc for r< r_paw

  real(dp), allocatable :: zioneff(:)
  ! zioneff(ij_proj), "Effective charge"*n "seen" at r_paw, deduced from Phi at r_paw, n:
  ! pricipal quantum number; good approximation to model wave function outside PAW-sphere through

 end type macro_uj_type
!!***

end module defs_abitypes
!!***
