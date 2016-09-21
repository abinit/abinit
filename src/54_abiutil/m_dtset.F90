!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dtset
!! NAME
!!  m_dtset
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1992-2016 ABINIT group (XG, MG)
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
 use m_profiling_abi
 use m_copy
 use m_errors

 use defs_abitypes, only : dataset_type

 implicit none

 private

 public :: dtset_copy
 public :: dtset_free

CONTAINS  !==============================================================================
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
!!      calc_vhxc_me,chkinp,dfpt_looppert,driver,gwls_hamiltonian,hybrid_corr
!!      m_io_kss,m_kxc,xchybrid_ncpp_cc
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtset_copy(dtout, dtin)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtset_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtin
 type(dataset_type),intent(out) :: dtout

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
 dtout%iomode          = dtin%iomode
 dtout%accuracy           = dtin%accuracy
 dtout%adpimd             = dtin%adpimd
 dtout%autoparal          = dtin%autoparal
 dtout%awtr               = dtin%awtr
 dtout%bandpp             = dtin%bandpp
 dtout%bdeigrf            = dtin%bdeigrf
 dtout%berryopt           = dtin%berryopt
 dtout%berrysav           = dtin%berrysav
 dtout%berrystep          = dtin%berrystep
 dtout%brvltt             = dtin%brvltt
 dtout%bs_nstates         = dtin%bs_nstates
 dtout%builtintest        = dtin%builtintest
 dtout%cd_customnimfrqs   = dtin%cd_customnimfrqs
 dtout%cd_frqim_method    = dtin%cd_frqim_method
 dtout%cd_full_grid       = dtin%cd_full_grid
 dtout%chkexit            = dtin%chkexit
 dtout%chkprim            = dtin%chkprim
 dtout%chksymbreak        = dtin%chksymbreak
 dtout%cineb_start        = dtin%cineb_start
 dtout%cgtyphf            = dtin%cgtyphf
 dtout%delayperm          = dtin%delayperm
 dtout%diismemory         = dtin%diismemory
 dtout%dmatpuopt          = dtin%dmatpuopt
 dtout%dmatudiag          = dtin%dmatudiag
 dtout%dmft_dc            = dtin%dmft_dc
 dtout%dmft_entropy       = dtin%dmft_entropy
 dtout%dmft_iter          = dtin%dmft_iter
 dtout%dmft_nlambda       = dtin%dmft_nlambda
 dtout%dmft_mxsf          = dtin%dmft_mxsf
 dtout%dmft_nwlo          = dtin%dmft_nwlo
 dtout%dmft_nwli          = dtin%dmft_nwli
 dtout%dmft_read_occnd    = dtin%dmft_read_occnd
 dtout%dmft_rslf          = dtin%dmft_rslf
 dtout%dmft_solv          = dtin%dmft_solv
 dtout%dmft_t2g           = dtin%dmft_t2g
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
 dtout%dmftqmc_n          = dtin%dmftqmc_n
 dtout%dmftqmc_l          = dtin%dmftqmc_l
 dtout%dmftqmc_seed       = dtin%dmftqmc_seed
 dtout%dmftqmc_therm      = dtin%dmftqmc_therm
 dtout%d3e_pert1_elfd     = dtin%d3e_pert1_elfd
 dtout%d3e_pert1_phon     = dtin%d3e_pert1_phon
 dtout%d3e_pert2_elfd     = dtin%d3e_pert2_elfd
 dtout%d3e_pert2_phon     = dtin%d3e_pert2_phon
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
 dtout%eph_extrael        = dtin%eph_extrael
 dtout%eph_fermie         = dtin%eph_fermie
 dtout%eph_fsmear         = dtin%eph_fsmear
 dtout%eph_fsewin         = dtin%eph_fsewin
 dtout%eph_ngqpt_fine     = dtin%eph_ngqpt_fine
 dtout%eph_task           = dtin%eph_task

 dtout%ph_wstep          = dtin%ph_wstep
 dtout%ph_intmeth        = dtin%ph_intmeth
 dtout%symdynmat         = dtin%symdynmat
 dtout%ph_nqshift        = dtin%ph_nqshift
 if (allocated(dtin%ph_qshift)) call alloc_copy(dtin%ph_qshift, dtout%ph_qshift)
 dtout%ph_smear          = dtin%ph_smear
 dtout%ddb_ngqpt         = dtin%ddb_ngqpt
 dtout%ddb_shiftq        = dtin%ddb_shiftq

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
 dtout%getddk             = dtin%getddk
 dtout%getdelfd           = dtin%getdelfd
 dtout%getdkdk            = dtin%getdkdk
 dtout%getdkde            = dtin%getdkde
 dtout%getden             = dtin%getden
 dtout%getgam_eig2nkq     = dtin%getgam_eig2nkq
 dtout%gethaydock         = dtin%gethaydock
 dtout%getocc             = dtin%getocc
 dtout%getpawden          = dtin%getpawden
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
 dtout%gwcalctyp          = dtin%gwcalctyp
 dtout%gwcomp             = dtin%gwcomp
 dtout%gwencomp           = dtin%gwencomp
 dtout%gwmem              = dtin%gwmem
 dtout%gwpara             = dtin%gwpara
 dtout%gwgamma            = dtin%gwgamma
 dtout%gwrpacorr          = dtin%gwrpacorr
 dtout%gw_customnfreqsp   = dtin%gw_customnfreqsp
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
 dtout%gwls_sternheimer_kmax= dtin%gwls_sternheimer_kmax
 dtout%gwls_npt_gauss_quad  = dtin%gwls_npt_gauss_quad
 dtout%gwls_dielectric_model= dtin%gwls_dielectric_model
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
 dtout%iboxcut            = dtin%iboxcut
 dtout%icoulomb           = dtin%icoulomb
 dtout%icutcoul           = dtin%icutcoul
 dtout%ieig2rf            = dtin%ieig2rf
 dtout%imgmov             = dtin%imgmov
 dtout%inclvkb            = dtin%inclvkb
 dtout%intxc              = dtin%intxc
 dtout%ionmov             = dtin%ionmov
 dtout%densfor_pred             = dtin%densfor_pred
 dtout%iprcel             = dtin%iprcel
 dtout%iprcfc             = dtin%iprcfc
 dtout%irandom            = dtin%irandom
 dtout%irdbseig           = dtin%irdbseig
 dtout%irdbsreso          = dtin%irdbsreso
 dtout%irdbscoup          = dtin%irdbscoup
 dtout%irdddb             = dtin%irdddb
 dtout%irdddk             = dtin%irdddk
 dtout%irdden             = dtin%irdden
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
 dtout%ixcpositron        = dtin%ixcpositron
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
 dtout%magconon           = dtin%magconon
 dtout%maxnsym            = dtin%maxnsym
 dtout%max_ncpus          = dtin%max_ncpus
 dtout%mband              = dtin%mband
 dtout%mdf_epsinf         = dtin%mdf_epsinf
 dtout%mep_solver         = dtin%mep_solver
 dtout%mem_test           = dtin%mem_test
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
 dtout%occopt             = dtin%occopt
 dtout%optcell            = dtin%optcell
 dtout%optdriver          = dtin%optdriver
 dtout%optforces          = dtin%optforces
 dtout%optnlxccc          = dtin%optnlxccc
 dtout%optstress          = dtin%optstress
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
 dtout%prtefg             = dtin%prtefg
 dtout%prteig             = dtin%prteig
 dtout%prtelf             = dtin%prtelf
 dtout%prtfc              = dtin%prtfc
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
 dtout%prtphdos           = dtin%prtphdos
 dtout%prtposcar          = dtin%prtposcar
 dtout%prtpot             = dtin%prtpot
 dtout%prtpsps            = dtin%prtpsps
 dtout%prtspcur           = dtin%prtspcur
 dtout%prtsuscep          = dtin%prtsuscep
 dtout%prtstm             = dtin%prtstm
 dtout%prtvclmb           = dtin%prtvclmb  
 dtout%prtvdw             = dtin%prtvdw
 dtout%prtvha             = dtin%prtvha
 dtout%prtvhxc            = dtin%prtvhxc
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
 dtout%timopt             = dtin%timopt
 dtout%use_gemm_nonlop    = dtin%use_gemm_nonlop
 dtout%use_gpu_cuda       = dtin%use_gpu_cuda
 dtout%use_slk            = dtin%use_slk
 dtout%usedmatpu          = dtin%usedmatpu
 dtout%usedmft            = dtin%usedmft
 dtout%useexexch          = dtin%useexexch
 dtout%usefock            = dtin%usefock
 dtout%usekden            = dtin%usekden
 dtout%use_nonscf_gkk     = dtin%use_nonscf_gkk
 dtout%usepaw             = dtin%usepaw
 dtout%usepawu            = dtin%usepawu
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
 dtout%supercell(:)       = dtin%supercell(:)
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
 dtout%effmass            = dtin%effmass
 dtout%efmas_deg_tol      = dtin%efmas_deg_tol
 dtout%elph2_imagden      = dtin%elph2_imagden
 dtout%eshift             = dtin%eshift
 dtout%esmear             = dtin%esmear
 dtout%exchmix            = dtin%exchmix
 dtout%fband              = dtin%fband
 dtout%gwls_model_parameter = dtin%gwls_model_parameter
 dtout%gwls_second_model_parameter = dtin%gwls_second_model_parameter
 dtout%spinmagntarget     = dtin%spinmagntarget
 dtout%friction           = dtin%friction
 dtout%fxcartfactor       = dtin%fxcartfactor
 dtout%ga_opt_percent     = dtin%ga_opt_percent
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
 dtout%ratsph_extra       = dtin%ratsph_extra
 dtout%recrcut            = dtin%recrcut
 dtout%recefermi          = dtin%recefermi
 dtout%rectolden          = dtin%rectolden
 dtout%dfpt_sciss         = dtin%dfpt_sciss
 dtout%mbpt_sciss         = dtin%mbpt_sciss
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
 dtout%tolmxde          = dtin%tolmxde
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
 call alloc_copy( dtin%algalch, dtout%algalch)

 call alloc_copy( dtin%bdgw, dtout%bdgw)

 call alloc_copy( dtin%bs_loband, dtout%bs_loband)

 call alloc_copy( dtin%dynimage, dtout%dynimage)

 call alloc_copy( dtin%efmas_bands, dtout%efmas_bands)

 call alloc_copy( dtin%iatfix, dtout%iatfix)

 call alloc_copy( dtin%iatsph, dtout%iatsph)

 call alloc_copy( dtin%istwfk, dtout%istwfk)

 call alloc_copy( dtin%kberry, dtout%kberry)

 call alloc_copy( dtin%lexexch, dtout%lexexch)

 call alloc_copy( dtin%lpawu, dtout%lpawu)

 call alloc_copy( dtin%nband, dtout%nband)

 call alloc_copy( dtin%plowan_iatom, dtout%plowan_iatom)

 call alloc_copy( dtin%plowan_it, dtout%plowan_it)

 call alloc_copy( dtin%plowan_nbl, dtout%plowan_nbl)

 call alloc_copy( dtin%plowan_lcalc, dtout%plowan_lcalc)

 call alloc_copy( dtin%plowan_projcalc, dtout%plowan_projcalc)

 call alloc_copy( dtin%prtatlist, dtout%prtatlist)

 call alloc_copy( dtin%so_psp, dtout%so_psp)

 call alloc_copy( dtin%symafm, dtout%symafm)

 call alloc_copy( dtin%symrel, dtout%symrel)

 call alloc_copy( dtin%typat, dtout%typat)

!Allocate and copy real allocatable
 call alloc_copy( dtin%acell_orig, dtout%acell_orig)

 call alloc_copy( dtin%amu_orig, dtout%amu_orig)

 call alloc_copy( dtin%atvshift, dtout%atvshift)

 call alloc_copy( dtin%cd_imfrqs, dtout%cd_imfrqs)

 call alloc_copy( dtin%corecs, dtout%corecs)

 call alloc_copy( dtin%densty, dtout%densty)

 call alloc_copy( dtin%dmatpawu, dtout%dmatpawu)

 call alloc_copy( dtin%efmas_dirs, dtout%efmas_dirs)

 call alloc_copy( dtin%f4of2_sla, dtout%f4of2_sla)

 call alloc_copy( dtin%f6of2_sla, dtout%f6of2_sla)

 call alloc_copy( dtin%gw_qlwl, dtout%gw_qlwl)

 call alloc_copy( dtin%gw_freqsp, dtout%gw_freqsp)

 call alloc_copy( dtin%gwls_list_proj_freq, dtout%gwls_list_proj_freq) 

 call alloc_copy( dtin%jpawu, dtout%jpawu)

 call alloc_copy( dtin%kpt, dtout%kpt)

 call alloc_copy( dtin%kptgw, dtout%kptgw)

 call alloc_copy( dtin%kptns, dtout%kptns)

 call alloc_copy( dtin%mixalch_orig, dtout%mixalch_orig)

 call alloc_copy( dtin%nucdipmom, dtout%nucdipmom)

 call alloc_copy( dtin%occ_orig, dtout%occ_orig)

 call alloc_copy( dtin%pimass, dtout%pimass)

 call alloc_copy( dtin%ptcharge, dtout%ptcharge)

 call alloc_copy( dtin%qmass, dtout%qmass)

 call alloc_copy( dtin%qptdm, dtout%qptdm)

 call alloc_copy( dtin%quadmom, dtout%quadmom)

 call alloc_copy( dtin%ratsph, dtout%ratsph)

 call alloc_copy( dtin%rprim_orig, dtout%rprim_orig)

 call alloc_copy( dtin%rprimd_orig, dtout%rprimd_orig)

 call alloc_copy( dtin%shiftk, dtout%shiftk)

 call alloc_copy( dtin%spinat, dtout%spinat)

 call alloc_copy( dtin%tnons, dtout%tnons)

 call alloc_copy( dtin%upawu, dtout%upawu)

 call alloc_copy( dtin%vel_orig, dtout%vel_orig)

 call alloc_copy( dtin%vel_cell_orig, dtout%vel_cell_orig)

 call alloc_copy( dtin%wtatcon, dtout%wtatcon)

 call alloc_copy( dtin%wtk, dtout%wtk)

 call alloc_copy( dtin%xred_orig, dtout%xred_orig)

 call alloc_copy( dtin%xredsph_extra, dtout%xredsph_extra)

 call alloc_copy( dtin%ziontypat, dtout%ziontypat)

 call alloc_copy( dtin%znucl, dtout%znucl)

 DBG_EXIT("COLL")

end subroutine dtset_copy
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
!!      calc_vhxc_me,chkinp,dfpt_looppert,driver,gwls_hamiltonian,hybrid_corr
!!      m_ab7_invars_f90,m_io_kss,m_kxc,xchybrid_ncpp_cc
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtset_free(dtset)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtset_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset

! *************************************************************************

!please, use the same order as the one used in the declaration of the type (see defs_abitypes).

 !@dataset_type
!integer allocatable
 if (allocated(dtset%algalch))     then
   ABI_DEALLOCATE(dtset%algalch)
 end if
 if (allocated(dtset%bdgw))        then
   ABI_DEALLOCATE(dtset%bdgw)
 end if
  if (allocated(dtset%bs_loband))  then
    ABI_DEALLOCATE(dtset%bs_loband)
  end if

 if (allocated(dtset%dynimage))    then
   ABI_DEALLOCATE(dtset%dynimage)
 end if
 if (allocated(dtset%efmas_bands))        then
   ABI_DEALLOCATE(dtset%efmas_bands)
 end if
 if (allocated(dtset%iatfix))      then
   ABI_DEALLOCATE(dtset%iatfix)
 end if
 if (allocated(dtset%iatsph))      then
   ABI_DEALLOCATE(dtset%iatsph)
 end if
 if (allocated(dtset%istwfk))      then
   ABI_DEALLOCATE(dtset%istwfk)
 end if
 if (allocated(dtset%kberry))      then
   ABI_DEALLOCATE(dtset%kberry)
 end if
 if (allocated(dtset%lexexch))     then
   ABI_DEALLOCATE(dtset%lexexch)
 end if
 if (allocated(dtset%lpawu))       then
   ABI_DEALLOCATE(dtset%lpawu)
 end if
 if (allocated(dtset%nband))       then
   ABI_DEALLOCATE(dtset%nband)
 end if
 if (allocated(dtset%ph_qpath)) then
   ABI_DEALLOCATE(dtset%ph_qpath)
 end if
 if (allocated(dtset%ph_qshift)) then
   ABI_DEALLOCATE(dtset%ph_qshift)
 end if
 if (allocated(dtset%plowan_iatom))       then
   ABI_DEALLOCATE(dtset%plowan_iatom)
 end if
 if (allocated(dtset%plowan_it))       then
   ABI_DEALLOCATE(dtset%plowan_it)
 end if
 if (allocated(dtset%plowan_lcalc))       then
   ABI_DEALLOCATE(dtset%plowan_lcalc)
 end if
 if (allocated(dtset%plowan_nbl))       then
   ABI_DEALLOCATE(dtset%plowan_nbl)
 end if
 if (allocated(dtset%plowan_projcalc))       then
   ABI_DEALLOCATE(dtset%plowan_projcalc)
 end if
 if (allocated(dtset%prtatlist))   then
   ABI_DEALLOCATE(dtset%prtatlist)
 end if
 if (allocated(dtset%so_psp))      then
   ABI_DEALLOCATE(dtset%so_psp)
 end if
 if (allocated(dtset%symafm))      then
   ABI_DEALLOCATE(dtset%symafm)
 end if
 if (allocated(dtset%symrel))      then
   ABI_DEALLOCATE(dtset%symrel)
 end if
 if (allocated(dtset%typat))       then
   ABI_DEALLOCATE(dtset%typat)
 end if

!real allocatable
 if (allocated(dtset%acell_orig))  then
   ABI_DEALLOCATE(dtset%acell_orig)
 end if
 if (allocated(dtset%amu_orig))    then
   ABI_DEALLOCATE(dtset%amu_orig)
 end if
 if (allocated(dtset%atvshift))    then
   ABI_DEALLOCATE(dtset%atvshift)
 end if
 if (allocated(dtset%cd_imfrqs))   then
   ABI_DEALLOCATE(dtset%cd_imfrqs)
 end if
 if (allocated(dtset%corecs))      then
   ABI_DEALLOCATE(dtset%corecs)
 end if
 if (allocated(dtset%densty))      then
   ABI_DEALLOCATE(dtset%densty)
 end if
 if (allocated(dtset%dmatpawu))    then
   ABI_DEALLOCATE(dtset%dmatpawu)
 end if
 if (allocated(dtset%efmas_dirs))        then
   ABI_DEALLOCATE(dtset%efmas_dirs)
 end if
 if (allocated(dtset%gw_qlwl))     then
   ABI_DEALLOCATE(dtset%gw_qlwl)
 end if
 if (allocated(dtset%gw_freqsp))   then
   ABI_DEALLOCATE(dtset%gw_freqsp)
 end if
 if (allocated(dtset%gwls_list_proj_freq))   then
   ABI_DEALLOCATE(dtset%gwls_list_proj_freq)
 end if
 if (allocated(dtset%f4of2_sla))   then
   ABI_DEALLOCATE(dtset%f4of2_sla)
 end if
 if (allocated(dtset%f6of2_sla))   then
   ABI_DEALLOCATE(dtset%f6of2_sla)
 end if
 if (allocated(dtset%jpawu))       then
   ABI_DEALLOCATE(dtset%jpawu)
 end if
 if (allocated(dtset%kpt))         then
   ABI_DEALLOCATE(dtset%kpt)
 end if
 if (allocated(dtset%kptgw))       then
   ABI_DEALLOCATE(dtset%kptgw)
 end if
 if (allocated(dtset%kptns))       then
   ABI_DEALLOCATE(dtset%kptns)
 end if
 if (allocated(dtset%mixalch_orig))     then
   ABI_DEALLOCATE(dtset%mixalch_orig)
 end if
 if (allocated(dtset%nucdipmom))      then
   ABI_DEALLOCATE(dtset%nucdipmom)
 end if
 if (allocated(dtset%occ_orig))    then
   ABI_DEALLOCATE(dtset%occ_orig)
 end if
 if (allocated(dtset%pimass))      then
   ABI_DEALLOCATE(dtset%pimass)
 end if
 if (allocated(dtset%ptcharge))    then
   ABI_DEALLOCATE(dtset%ptcharge)
 end if
 if (allocated(dtset%qmass))       then
   ABI_DEALLOCATE(dtset%qmass)
 end if
 if (allocated(dtset%qptdm))       then
   ABI_DEALLOCATE(dtset%qptdm)
 end if
 if (allocated(dtset%quadmom))     then
   ABI_DEALLOCATE(dtset%quadmom)
 end if
 if (allocated(dtset%ratsph))      then
   ABI_DEALLOCATE(dtset%ratsph)
 end if
 if (allocated(dtset%rprim_orig))  then
   ABI_DEALLOCATE(dtset%rprim_orig)
 end if
 if (allocated(dtset%rprimd_orig)) then
   ABI_DEALLOCATE(dtset%rprimd_orig)
 end if
 if (allocated(dtset%shiftk))      then
   ABI_DEALLOCATE(dtset%shiftk)
 end if
 if (allocated(dtset%spinat))      then
   ABI_DEALLOCATE(dtset%spinat)
 end if
 if (allocated(dtset%tnons))       then
   ABI_DEALLOCATE(dtset%tnons)
 end if
 if (allocated(dtset%upawu))       then
   ABI_DEALLOCATE(dtset%upawu)
 end if
 if (allocated(dtset%vel_orig))    then
   ABI_DEALLOCATE(dtset%vel_orig)
 end if
 if (allocated(dtset%vel_cell_orig))    then
   ABI_DEALLOCATE(dtset%vel_cell_orig)
 end if
 if (allocated(dtset%wtatcon))     then
   ABI_DEALLOCATE(dtset%wtatcon)
 end if
 if (allocated(dtset%wtk))         then
   ABI_DEALLOCATE(dtset%wtk)
 end if
 if (allocated(dtset%xred_orig))   then
   ABI_DEALLOCATE(dtset%xred_orig)
 end if
 if (allocated(dtset%xredsph_extra))   then
   ABI_DEALLOCATE(dtset%xredsph_extra)
 end if
 if (allocated(dtset%ziontypat))   then
   ABI_DEALLOCATE(dtset%ziontypat)
 end if
 if (allocated(dtset%znucl))       then
   ABI_DEALLOCATE(dtset%znucl)
 end if

end subroutine dtset_free
!!***

END MODULE m_dtset
!!***
