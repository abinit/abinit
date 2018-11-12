! This file has been automatically generated, do not modify.

  subroutine ab7_invars_get_integer(dtsetsId, value, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    integer, intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB7_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB7_ERROR_INVARS_ID
      return
    end if
    
    errno = AB7_NO_ERROR
    select case (id)
    case (ab7_invars_iomode)
      value = token%dtsets(idtset)%iomode
    case (ab7_invars_accuracy)
      value = token%dtsets(idtset)%accuracy
    case (ab7_invars_adpimd)
      value = token%dtsets(idtset)%adpimd
    case (ab7_invars_autoparal)
      value = token%dtsets(idtset)%autoparal
    case (ab7_invars_auxc_ixc)
      value = token%dtsets(idtset)%auxc_ixc
    case (ab7_invars_awtr)
      value = token%dtsets(idtset)%awtr
    case (ab7_invars_bandpp)
      value = token%dtsets(idtset)%bandpp
    case (ab7_invars_bdeigrf)
      value = token%dtsets(idtset)%bdeigrf
    case (ab7_invars_berryopt)
      value = token%dtsets(idtset)%berryopt
    case (ab7_invars_berrysav)
      value = token%dtsets(idtset)%berrysav
    case (ab7_invars_berrystep)
      value = token%dtsets(idtset)%berrystep
    case (ab7_invars_brvltt)
      value = token%dtsets(idtset)%brvltt
    case (ab7_invars_bs_nstates)
      value = token%dtsets(idtset)%bs_nstates
    case (ab7_invars_bs_hayd_term)
      value = token%dtsets(idtset)%bs_hayd_term
    case (ab7_invars_builtintest)
      value = token%dtsets(idtset)%builtintest
    case (ab7_invars_cd_full_grid)
      value = token%dtsets(idtset)%cd_full_grid
    case (ab7_invars_cd_frqim_method)
      value = token%dtsets(idtset)%cd_frqim_method
    case (ab7_invars_cd_customnimfrqs)
      value = token%dtsets(idtset)%cd_customnimfrqs
    case (ab7_invars_chkdilatmx)
      value = token%dtsets(idtset)%chkdilatmx
    case (ab7_invars_chkexit)
      value = token%dtsets(idtset)%chkexit
    case (ab7_invars_chkprim)
      value = token%dtsets(idtset)%chkprim
    case (ab7_invars_chksymbreak)
      value = token%dtsets(idtset)%chksymbreak
    case (ab7_invars_cineb_start)
      value = token%dtsets(idtset)%cineb_start
    case (ab7_invars_delayperm)
      value = token%dtsets(idtset)%delayperm
    case (ab7_invars_densfor_pred)
      value = token%dtsets(idtset)%densfor_pred
    case (ab7_invars_diismemory)
      value = token%dtsets(idtset)%diismemory
    case (ab7_invars_dmatpuopt)
      value = token%dtsets(idtset)%dmatpuopt
    case (ab7_invars_dmatudiag)
      value = token%dtsets(idtset)%dmatudiag
    case (ab7_invars_dmft_dc)
      value = token%dtsets(idtset)%dmft_dc
    case (ab7_invars_dmft_entropy)
      value = token%dtsets(idtset)%dmft_entropy
    case (ab7_invars_dmft_iter)
      value = token%dtsets(idtset)%dmft_iter
    case (ab7_invars_dmft_nlambda)
      value = token%dtsets(idtset)%dmft_nlambda
    case (ab7_invars_dmft_nwli)
      value = token%dtsets(idtset)%dmft_nwli
    case (ab7_invars_dmft_nwlo)
      value = token%dtsets(idtset)%dmft_nwlo
    case (ab7_invars_dmft_rslf)
      value = token%dtsets(idtset)%dmft_rslf
    case (ab7_invars_dmft_read_occnd)
      value = token%dtsets(idtset)%dmft_read_occnd
    case (ab7_invars_dmft_solv)
      value = token%dtsets(idtset)%dmft_solv
    case (ab7_invars_dmft_t2g)
      value = token%dtsets(idtset)%dmft_t2g
    case (ab7_invars_dmftbandi)
      value = token%dtsets(idtset)%dmftbandi
    case (ab7_invars_dmftbandf)
      value = token%dtsets(idtset)%dmftbandf
    case (ab7_invars_dmftcheck)
      value = token%dtsets(idtset)%dmftcheck
    case (ab7_invars_dmftctqmc_basis)
      value = token%dtsets(idtset)%dmftctqmc_basis
    case (ab7_invars_dmftctqmc_check)
      value = token%dtsets(idtset)%dmftctqmc_check
    case (ab7_invars_dmftctqmc_correl)
      value = token%dtsets(idtset)%dmftctqmc_correl
    case (ab7_invars_dmftctqmc_gmove)
      value = token%dtsets(idtset)%dmftctqmc_gmove
    case (ab7_invars_dmftctqmc_grnns)
      value = token%dtsets(idtset)%dmftctqmc_grnns
    case (ab7_invars_dmftctqmc_meas)
      value = token%dtsets(idtset)%dmftctqmc_meas
    case (ab7_invars_dmftctqmc_mov)
      value = token%dtsets(idtset)%dmftctqmc_mov
    case (ab7_invars_dmftctqmc_mrka)
      value = token%dtsets(idtset)%dmftctqmc_mrka
    case (ab7_invars_dmftctqmc_order)
      value = token%dtsets(idtset)%dmftctqmc_order
    case (ab7_invars_dmftctqmc_triqs_nleg)
      value = token%dtsets(idtset)%dmftctqmc_triqs_nleg
    case (ab7_invars_dmftqmc_l)
      value = token%dtsets(idtset)%dmftqmc_l
    case (ab7_invars_dmftqmc_seed)
      value = token%dtsets(idtset)%dmftqmc_seed
    case (ab7_invars_dmftqmc_therm)
      value = token%dtsets(idtset)%dmftqmc_therm
    case (ab7_invars_d3e_pert1_elfd)
      value = token%dtsets(idtset)%d3e_pert1_elfd
    case (ab7_invars_d3e_pert1_phon)
      value = token%dtsets(idtset)%d3e_pert1_phon
    case (ab7_invars_d3e_pert2_elfd)
      value = token%dtsets(idtset)%d3e_pert2_elfd
    case (ab7_invars_d3e_pert2_phon)
      value = token%dtsets(idtset)%d3e_pert2_phon
    case (ab7_invars_d3e_pert3_elfd)
      value = token%dtsets(idtset)%d3e_pert3_elfd
    case (ab7_invars_d3e_pert3_phon)
      value = token%dtsets(idtset)%d3e_pert3_phon
    case (ab7_invars_efmas)
      value = token%dtsets(idtset)%efmas
    case (ab7_invars_efmas_calc_dirs)
      value = token%dtsets(idtset)%efmas_calc_dirs
    case (ab7_invars_efmas_deg)
      value = token%dtsets(idtset)%efmas_deg
    case (ab7_invars_efmas_dim)
      value = token%dtsets(idtset)%efmas_dim
    case (ab7_invars_efmas_n_dirs)
      value = token%dtsets(idtset)%efmas_n_dirs
    case (ab7_invars_efmas_ntheta)
      value = token%dtsets(idtset)%efmas_ntheta
    case (ab7_invars_enunit)
      value = token%dtsets(idtset)%enunit
    case (ab7_invars_eph_task)
      value = token%dtsets(idtset)%eph_task
    case (ab7_invars_exchn2n3d)
      value = token%dtsets(idtset)%exchn2n3d
    case (ab7_invars_extrapwf)
      value = token%dtsets(idtset)%extrapwf
    case (ab7_invars_fftgw)
      value = token%dtsets(idtset)%fftgw
    case (ab7_invars_fockoptmix)
      value = token%dtsets(idtset)%fockoptmix
    case (ab7_invars_frzfermi)
      value = token%dtsets(idtset)%frzfermi
    case (ab7_invars_ga_algor)
      value = token%dtsets(idtset)%ga_algor
    case (ab7_invars_ga_fitness)
      value = token%dtsets(idtset)%ga_fitness
    case (ab7_invars_ga_n_rules)
      value = token%dtsets(idtset)%ga_n_rules
    case (ab7_invars_getcell)
      value = token%dtsets(idtset)%getcell
    case (ab7_invars_getddb)
      value = token%dtsets(idtset)%getddb
    case (ab7_invars_getddk)
      value = token%dtsets(idtset)%getddk
    case (ab7_invars_getdelfd)
      value = token%dtsets(idtset)%getdelfd
    case (ab7_invars_getdkdk)
      value = token%dtsets(idtset)%getdkdk
    case (ab7_invars_getdkde)
      value = token%dtsets(idtset)%getdkde
    case (ab7_invars_getden)
      value = token%dtsets(idtset)%getden
    case (ab7_invars_getefmas)
      value = token%dtsets(idtset)%getefmas
    case (ab7_invars_getgam_eig2nkq)
      value = token%dtsets(idtset)%getgam_eig2nkq
    case (ab7_invars_getocc)
      value = token%dtsets(idtset)%getocc
    case (ab7_invars_getpawden)
      value = token%dtsets(idtset)%getpawden
    case (ab7_invars_getqps)
      value = token%dtsets(idtset)%getqps
    case (ab7_invars_getscr)
      value = token%dtsets(idtset)%getscr
    case (ab7_invars_getsuscep)
      value = token%dtsets(idtset)%getsuscep
    case (ab7_invars_getvel)
      value = token%dtsets(idtset)%getvel
    case (ab7_invars_getwfk)
      value = token%dtsets(idtset)%getwfk
    case (ab7_invars_getwfkfine)
      value = token%dtsets(idtset)%getwfkfine
    case (ab7_invars_getwfq)
      value = token%dtsets(idtset)%getwfq
    case (ab7_invars_getxcart)
      value = token%dtsets(idtset)%getxcart
    case (ab7_invars_getxred)
      value = token%dtsets(idtset)%getxred
    case (ab7_invars_get1den)
      value = token%dtsets(idtset)%get1den
    case (ab7_invars_get1wf)
      value = token%dtsets(idtset)%get1wf
    case (ab7_invars_getbseig)
      value = token%dtsets(idtset)%getbseig
    case (ab7_invars_getbsreso)
      value = token%dtsets(idtset)%getbsreso
    case (ab7_invars_getbscoup)
      value = token%dtsets(idtset)%getbscoup
    case (ab7_invars_gethaydock)
      value = token%dtsets(idtset)%gethaydock
    case (ab7_invars_goprecon)
      value = token%dtsets(idtset)%goprecon
    case (ab7_invars_gwcalctyp)
      value = token%dtsets(idtset)%gwcalctyp
    case (ab7_invars_gwcomp)
      value = token%dtsets(idtset)%gwcomp
    case (ab7_invars_gwgamma)
      value = token%dtsets(idtset)%gwgamma
    case (ab7_invars_gwrpacorr)
      value = token%dtsets(idtset)%gwrpacorr
    case (ab7_invars_gw_customnfreqsp)
      value = token%dtsets(idtset)%gw_customnfreqsp
    case (ab7_invars_gw_invalid_freq)
      value = token%dtsets(idtset)%gw_invalid_freq
    case (ab7_invars_gw_qprange)
      value = token%dtsets(idtset)%gw_qprange
    case (ab7_invars_gw_nqlwl)
      value = token%dtsets(idtset)%gw_nqlwl
    case (ab7_invars_gw_nstep)
      value = token%dtsets(idtset)%gw_nstep
    case (ab7_invars_gw_sigxcore)
      value = token%dtsets(idtset)%gw_sigxcore
    case (ab7_invars_gwls_stern_kmax)
      value = token%dtsets(idtset)%gwls_stern_kmax
    case (ab7_invars_gwls_npt_gauss_quad)
      value = token%dtsets(idtset)%gwls_npt_gauss_quad
    case (ab7_invars_gwls_diel_model)
      value = token%dtsets(idtset)%gwls_diel_model
    case (ab7_invars_gwls_print_debug)
      value = token%dtsets(idtset)%gwls_print_debug
    case (ab7_invars_gwls_nseeds)
      value = token%dtsets(idtset)%gwls_nseeds
    case (ab7_invars_gwls_n_proj_freq)
      value = token%dtsets(idtset)%gwls_n_proj_freq
    case (ab7_invars_gwls_kmax_complement)
      value = token%dtsets(idtset)%gwls_kmax_complement
    case (ab7_invars_gwls_kmax_poles)
      value = token%dtsets(idtset)%gwls_kmax_poles
    case (ab7_invars_gwls_kmax_analytic)
      value = token%dtsets(idtset)%gwls_kmax_analytic
    case (ab7_invars_gwls_kmax_numeric)
      value = token%dtsets(idtset)%gwls_kmax_numeric
    case (ab7_invars_gwls_band_index)
      value = token%dtsets(idtset)%gwls_band_index
    case (ab7_invars_gwls_exchange)
      value = token%dtsets(idtset)%gwls_exchange
    case (ab7_invars_gwls_correlation)
      value = token%dtsets(idtset)%gwls_correlation
    case (ab7_invars_gwls_first_seed)
      value = token%dtsets(idtset)%gwls_first_seed
    case (ab7_invars_gwls_recycle)
      value = token%dtsets(idtset)%gwls_recycle
    case (ab7_invars_gw_frqim_inzgrid)
      value = token%dtsets(idtset)%gw_frqim_inzgrid
    case (ab7_invars_gw_frqre_inzgrid)
      value = token%dtsets(idtset)%gw_frqre_inzgrid
    case (ab7_invars_gw_frqre_tangrid)
      value = token%dtsets(idtset)%gw_frqre_tangrid
    case (ab7_invars_gw_sctype)
      value = token%dtsets(idtset)%gw_sctype
    case (ab7_invars_gwmem)
      value = token%dtsets(idtset)%gwmem
    case (ab7_invars_gwpara)
      value = token%dtsets(idtset)%gwpara
    case (ab7_invars_hmcsst)
      value = token%dtsets(idtset)%hmcsst
    case (ab7_invars_hmctt)
      value = token%dtsets(idtset)%hmctt
    case (ab7_invars_iboxcut)
      value = token%dtsets(idtset)%iboxcut
    case (ab7_invars_icoulomb)
      value = token%dtsets(idtset)%icoulomb
    case (ab7_invars_icutcoul)
      value = token%dtsets(idtset)%icutcoul
    case (ab7_invars_ieig2rf)
      value = token%dtsets(idtset)%ieig2rf
    case (ab7_invars_imgmov)
      value = token%dtsets(idtset)%imgmov
    case (ab7_invars_imgwfstor)
      value = token%dtsets(idtset)%imgwfstor
    case (ab7_invars_inclvkb)
      value = token%dtsets(idtset)%inclvkb
    case (ab7_invars_intxc)
      value = token%dtsets(idtset)%intxc
    case (ab7_invars_ionmov)
      value = token%dtsets(idtset)%ionmov
    case (ab7_invars_iprcel)
      value = token%dtsets(idtset)%iprcel
    case (ab7_invars_iprcfc)
      value = token%dtsets(idtset)%iprcfc
    case (ab7_invars_irandom)
      value = token%dtsets(idtset)%irandom
    case (ab7_invars_irdddb)
      value = token%dtsets(idtset)%irdddb
    case (ab7_invars_irdddk)
      value = token%dtsets(idtset)%irdddk
    case (ab7_invars_irdden)
      value = token%dtsets(idtset)%irdden
    case (ab7_invars_irdefmas)
      value = token%dtsets(idtset)%irdefmas
    case (ab7_invars_irdhaydock)
      value = token%dtsets(idtset)%irdhaydock
    case (ab7_invars_irdpawden)
      value = token%dtsets(idtset)%irdpawden
    case (ab7_invars_irdqps)
      value = token%dtsets(idtset)%irdqps
    case (ab7_invars_irdscr)
      value = token%dtsets(idtset)%irdscr
    case (ab7_invars_irdsuscep)
      value = token%dtsets(idtset)%irdsuscep
    case (ab7_invars_irdvdw)
      value = token%dtsets(idtset)%irdvdw
    case (ab7_invars_irdwfk)
      value = token%dtsets(idtset)%irdwfk
    case (ab7_invars_irdwfkfine)
      value = token%dtsets(idtset)%irdwfkfine
    case (ab7_invars_irdwfq)
      value = token%dtsets(idtset)%irdwfq
    case (ab7_invars_ird1den)
      value = token%dtsets(idtset)%ird1den
    case (ab7_invars_ird1wf)
      value = token%dtsets(idtset)%ird1wf
    case (ab7_invars_irdbseig)
      value = token%dtsets(idtset)%irdbseig
    case (ab7_invars_irdbsreso)
      value = token%dtsets(idtset)%irdbsreso
    case (ab7_invars_irdbscoup)
      value = token%dtsets(idtset)%irdbscoup
    case (ab7_invars_iscf)
      value = token%dtsets(idtset)%iscf
    case (ab7_invars_isecur)
      value = token%dtsets(idtset)%isecur
    case (ab7_invars_istatimg)
      value = token%dtsets(idtset)%istatimg
    case (ab7_invars_istatr)
      value = token%dtsets(idtset)%istatr
    case (ab7_invars_istatshft)
      value = token%dtsets(idtset)%istatshft
    case (ab7_invars_ixc)
      value = token%dtsets(idtset)%ixc
    case (ab7_invars_ixc_sigma)
      value = token%dtsets(idtset)%ixc_sigma
    case (ab7_invars_ixcpositron)
      value = token%dtsets(idtset)%ixcpositron
    case (ab7_invars_ixcrot)
      value = token%dtsets(idtset)%ixcrot
    case (ab7_invars_jdtset)
      value = token%dtsets(idtset)%jdtset
    case (ab7_invars_jellslab)
      value = token%dtsets(idtset)%jellslab
    case (ab7_invars_kptopt)
      value = token%dtsets(idtset)%kptopt
    case (ab7_invars_kssform)
      value = token%dtsets(idtset)%kssform
    case (ab7_invars_localrdwf)
      value = token%dtsets(idtset)%localrdwf
    case (ab7_invars_lotf_classic)
      value = token%dtsets(idtset)%lotf_classic
    case (ab7_invars_lotf_nitex)
      value = token%dtsets(idtset)%lotf_nitex
    case (ab7_invars_lotf_nneigx)
      value = token%dtsets(idtset)%lotf_nneigx
    case (ab7_invars_lotf_version)
      value = token%dtsets(idtset)%lotf_version
    case (ab7_invars_magconon)
      value = token%dtsets(idtset)%magconon
    case (ab7_invars_maxnsym)
      value = token%dtsets(idtset)%maxnsym
    case (ab7_invars_max_ncpus)
      value = token%dtsets(idtset)%max_ncpus
    case (ab7_invars_mband)
      value = token%dtsets(idtset)%mband
    case (ab7_invars_mep_solver)
      value = token%dtsets(idtset)%mep_solver
    case (ab7_invars_mem_test)
      value = token%dtsets(idtset)%mem_test
    case (ab7_invars_mffmem)
      value = token%dtsets(idtset)%mffmem
    case (ab7_invars_mgfft)
      value = token%dtsets(idtset)%mgfft
    case (ab7_invars_mgfftdg)
      value = token%dtsets(idtset)%mgfftdg
    case (ab7_invars_mkmem)
      value = token%dtsets(idtset)%mkmem
    case (ab7_invars_mkqmem)
      value = token%dtsets(idtset)%mkqmem
    case (ab7_invars_mk1mem)
      value = token%dtsets(idtset)%mk1mem
    case (ab7_invars_nnos)
      value = token%dtsets(idtset)%nnos
    case (ab7_invars_mpw)
      value = token%dtsets(idtset)%mpw
    case (ab7_invars_mqgrid)
      value = token%dtsets(idtset)%mqgrid
    case (ab7_invars_mqgriddg)
      value = token%dtsets(idtset)%mqgriddg
    case (ab7_invars_natom)
      value = token%dtsets(idtset)%natom
    case (ab7_invars_natpawu)
      value = token%dtsets(idtset)%natpawu
    case (ab7_invars_natrd)
      value = token%dtsets(idtset)%natrd
    case (ab7_invars_natsph)
      value = token%dtsets(idtset)%natsph
    case (ab7_invars_natsph_extra)
      value = token%dtsets(idtset)%natsph_extra
    case (ab7_invars_natvshift)
      value = token%dtsets(idtset)%natvshift
    case (ab7_invars_nbandhf)
      value = token%dtsets(idtset)%nbandhf
    case (ab7_invars_nbandkss)
      value = token%dtsets(idtset)%nbandkss
    case (ab7_invars_nbdblock)
      value = token%dtsets(idtset)%nbdblock
    case (ab7_invars_nbdbuf)
      value = token%dtsets(idtset)%nbdbuf
    case (ab7_invars_nberry)
      value = token%dtsets(idtset)%nberry
    case (ab7_invars_nc_xccc_gspace)
      value = token%dtsets(idtset)%nc_xccc_gspace
    case (ab7_invars_nconeq)
      value = token%dtsets(idtset)%nconeq
    case (ab7_invars_nctime)
      value = token%dtsets(idtset)%nctime
    case (ab7_invars_ndtset)
      value = token%dtsets(idtset)%ndtset
    case (ab7_invars_ndynimage)
      value = token%dtsets(idtset)%ndynimage
    case (ab7_invars_neb_algo)
      value = token%dtsets(idtset)%neb_algo
    case (ab7_invars_nfft)
      value = token%dtsets(idtset)%nfft
    case (ab7_invars_nfftdg)
      value = token%dtsets(idtset)%nfftdg
    case (ab7_invars_nfreqim)
      value = token%dtsets(idtset)%nfreqim
    case (ab7_invars_nfreqre)
      value = token%dtsets(idtset)%nfreqre
    case (ab7_invars_nfreqsp)
      value = token%dtsets(idtset)%nfreqsp
    case (ab7_invars_nimage)
      value = token%dtsets(idtset)%nimage
    case (ab7_invars_nkpt)
      value = token%dtsets(idtset)%nkpt
    case (ab7_invars_nkptgw)
      value = token%dtsets(idtset)%nkptgw
    case (ab7_invars_nkpthf)
      value = token%dtsets(idtset)%nkpthf
    case (ab7_invars_nline)
      value = token%dtsets(idtset)%nline
    case (ab7_invars_nnsclo)
      value = token%dtsets(idtset)%nnsclo
    case (ab7_invars_nnsclohf)
      value = token%dtsets(idtset)%nnsclohf
    case (ab7_invars_nomegasf)
      value = token%dtsets(idtset)%nomegasf
    case (ab7_invars_nomegasi)
      value = token%dtsets(idtset)%nomegasi
    case (ab7_invars_nomegasrd)
      value = token%dtsets(idtset)%nomegasrd
    case (ab7_invars_nonlinear_info)
      value = token%dtsets(idtset)%nonlinear_info
    case (ab7_invars_npband)
      value = token%dtsets(idtset)%npband
    case (ab7_invars_npfft)
      value = token%dtsets(idtset)%npfft
    case (ab7_invars_nphf)
      value = token%dtsets(idtset)%nphf
    case (ab7_invars_npimage)
      value = token%dtsets(idtset)%npimage
    case (ab7_invars_npkpt)
      value = token%dtsets(idtset)%npkpt
    case (ab7_invars_nppert)
      value = token%dtsets(idtset)%nppert
    case (ab7_invars_npspinor)
      value = token%dtsets(idtset)%npspinor
    case (ab7_invars_npsp)
      value = token%dtsets(idtset)%npsp
    case (ab7_invars_npspalch)
      value = token%dtsets(idtset)%npspalch
    case (ab7_invars_npulayit)
      value = token%dtsets(idtset)%npulayit
    case (ab7_invars_npvel)
      value = token%dtsets(idtset)%npvel
    case (ab7_invars_npweps)
      value = token%dtsets(idtset)%npweps
    case (ab7_invars_npwkss)
      value = token%dtsets(idtset)%npwkss
    case (ab7_invars_npwsigx)
      value = token%dtsets(idtset)%npwsigx
    case (ab7_invars_npwwfn)
      value = token%dtsets(idtset)%npwwfn
    case (ab7_invars_np_slk)
      value = token%dtsets(idtset)%np_slk
    case (ab7_invars_nqpt)
      value = token%dtsets(idtset)%nqpt
    case (ab7_invars_nqptdm)
      value = token%dtsets(idtset)%nqptdm
    case (ab7_invars_nscforder)
      value = token%dtsets(idtset)%nscforder
    case (ab7_invars_nshiftk)
      value = token%dtsets(idtset)%nshiftk
    case (ab7_invars_nshiftk_orig)
      value = token%dtsets(idtset)%nshiftk_orig
    case (ab7_invars_nspden)
      value = token%dtsets(idtset)%nspden
    case (ab7_invars_nspinor)
      value = token%dtsets(idtset)%nspinor
    case (ab7_invars_nsppol)
      value = token%dtsets(idtset)%nsppol
    case (ab7_invars_nstep)
      value = token%dtsets(idtset)%nstep
    case (ab7_invars_nsym)
      value = token%dtsets(idtset)%nsym
    case (ab7_invars_ntime)
      value = token%dtsets(idtset)%ntime
    case (ab7_invars_ntimimage)
      value = token%dtsets(idtset)%ntimimage
    case (ab7_invars_ntypalch)
      value = token%dtsets(idtset)%ntypalch
    case (ab7_invars_ntypat)
      value = token%dtsets(idtset)%ntypat
    case (ab7_invars_ntyppure)
      value = token%dtsets(idtset)%ntyppure
    case (ab7_invars_nwfshist)
      value = token%dtsets(idtset)%nwfshist
    case (ab7_invars_nzchempot)
      value = token%dtsets(idtset)%nzchempot
    case (ab7_invars_occopt)
      value = token%dtsets(idtset)%occopt
    case (ab7_invars_optcell)
      value = token%dtsets(idtset)%optcell
    case (ab7_invars_optdriver)
      value = token%dtsets(idtset)%optdriver
    case (ab7_invars_optforces)
      value = token%dtsets(idtset)%optforces
    case (ab7_invars_optnlxccc)
      value = token%dtsets(idtset)%optnlxccc
    case (ab7_invars_optstress)
      value = token%dtsets(idtset)%optstress
    case (ab7_invars_orbmag)
      value = token%dtsets(idtset)%orbmag
    case (ab7_invars_ortalg)
      value = token%dtsets(idtset)%ortalg
    case (ab7_invars_paral_atom)
      value = token%dtsets(idtset)%paral_atom
    case (ab7_invars_paral_kgb)
      value = token%dtsets(idtset)%paral_kgb
    case (ab7_invars_paral_rf)
      value = token%dtsets(idtset)%paral_rf
    case (ab7_invars_pawcpxocc)
      value = token%dtsets(idtset)%pawcpxocc
    case (ab7_invars_pawcross)
      value = token%dtsets(idtset)%pawcross
    case (ab7_invars_pawfatbnd)
      value = token%dtsets(idtset)%pawfatbnd
    case (ab7_invars_pawlcutd)
      value = token%dtsets(idtset)%pawlcutd
    case (ab7_invars_pawlmix)
      value = token%dtsets(idtset)%pawlmix
    case (ab7_invars_pawmixdg)
      value = token%dtsets(idtset)%pawmixdg
    case (ab7_invars_pawnhatxc)
      value = token%dtsets(idtset)%pawnhatxc
    case (ab7_invars_pawnphi)
      value = token%dtsets(idtset)%pawnphi
    case (ab7_invars_pawntheta)
      value = token%dtsets(idtset)%pawntheta
    case (ab7_invars_pawnzlm)
      value = token%dtsets(idtset)%pawnzlm
    case (ab7_invars_pawoptmix)
      value = token%dtsets(idtset)%pawoptmix
    case (ab7_invars_pawoptosc)
      value = token%dtsets(idtset)%pawoptosc
    case (ab7_invars_pawprtdos)
      value = token%dtsets(idtset)%pawprtdos
    case (ab7_invars_pawprtvol)
      value = token%dtsets(idtset)%pawprtvol
    case (ab7_invars_pawprtwf)
      value = token%dtsets(idtset)%pawprtwf
    case (ab7_invars_pawprt_k)
      value = token%dtsets(idtset)%pawprt_k
    case (ab7_invars_pawprt_b)
      value = token%dtsets(idtset)%pawprt_b
    case (ab7_invars_pawspnorb)
      value = token%dtsets(idtset)%pawspnorb
    case (ab7_invars_pawstgylm)
      value = token%dtsets(idtset)%pawstgylm
    case (ab7_invars_pawsushat)
      value = token%dtsets(idtset)%pawsushat
    case (ab7_invars_pawusecp)
      value = token%dtsets(idtset)%pawusecp
    case (ab7_invars_macro_uj)
      value = token%dtsets(idtset)%macro_uj
    case (ab7_invars_pawujat)
      value = token%dtsets(idtset)%pawujat
    case (ab7_invars_pawxcdev)
      value = token%dtsets(idtset)%pawxcdev
    case (ab7_invars_pimd_constraint)
      value = token%dtsets(idtset)%pimd_constraint
    case (ab7_invars_pitransform)
      value = token%dtsets(idtset)%pitransform
    case (ab7_invars_plowan_bandi)
      value = token%dtsets(idtset)%plowan_bandi
    case (ab7_invars_plowan_bandf)
      value = token%dtsets(idtset)%plowan_bandf
    case (ab7_invars_plowan_compute)
      value = token%dtsets(idtset)%plowan_compute
    case (ab7_invars_plowan_natom)
      value = token%dtsets(idtset)%plowan_natom
    case (ab7_invars_plowan_nt)
      value = token%dtsets(idtset)%plowan_nt
    case (ab7_invars_plowan_realspace)
      value = token%dtsets(idtset)%plowan_realspace
    case (ab7_invars_posdoppler)
      value = token%dtsets(idtset)%posdoppler
    case (ab7_invars_positron)
      value = token%dtsets(idtset)%positron
    case (ab7_invars_posnstep)
      value = token%dtsets(idtset)%posnstep
    case (ab7_invars_ppmodel)
      value = token%dtsets(idtset)%ppmodel
    case (ab7_invars_prepanl)
      value = token%dtsets(idtset)%prepanl
    case (ab7_invars_prepgkk)
      value = token%dtsets(idtset)%prepgkk
    case (ab7_invars_prtbbb)
      value = token%dtsets(idtset)%prtbbb
    case (ab7_invars_prtbltztrp)
      value = token%dtsets(idtset)%prtbltztrp
    case (ab7_invars_prtcif)
      value = token%dtsets(idtset)%prtcif
    case (ab7_invars_prtden)
      value = token%dtsets(idtset)%prtden
    case (ab7_invars_prtdensph)
      value = token%dtsets(idtset)%prtdensph
    case (ab7_invars_prtdipole)
      value = token%dtsets(idtset)%prtdipole
    case (ab7_invars_prtdos)
      value = token%dtsets(idtset)%prtdos
    case (ab7_invars_prtdosm)
      value = token%dtsets(idtset)%prtdosm
    case (ab7_invars_prtebands)
      value = token%dtsets(idtset)%prtebands
    case (ab7_invars_prtefg)
      value = token%dtsets(idtset)%prtefg
    case (ab7_invars_prtefmas)
      value = token%dtsets(idtset)%prtefmas
    case (ab7_invars_prteig)
      value = token%dtsets(idtset)%prteig
    case (ab7_invars_prtelf)
      value = token%dtsets(idtset)%prtelf
    case (ab7_invars_prtfc)
      value = token%dtsets(idtset)%prtfc
    case (ab7_invars_prtfull1wf)
      value = token%dtsets(idtset)%prtfull1wf
    case (ab7_invars_prtfsurf)
      value = token%dtsets(idtset)%prtfsurf
    case (ab7_invars_prtgsr)
      value = token%dtsets(idtset)%prtgsr
    case (ab7_invars_prtgden)
      value = token%dtsets(idtset)%prtgden
    case (ab7_invars_prtgeo)
      value = token%dtsets(idtset)%prtgeo
    case (ab7_invars_prtgkk)
      value = token%dtsets(idtset)%prtgkk
    case (ab7_invars_prtkden)
      value = token%dtsets(idtset)%prtkden
    case (ab7_invars_prtkpt)
      value = token%dtsets(idtset)%prtkpt
    case (ab7_invars_prtlden)
      value = token%dtsets(idtset)%prtlden
    case (ab7_invars_prtnabla)
      value = token%dtsets(idtset)%prtnabla
    case (ab7_invars_prtnest)
      value = token%dtsets(idtset)%prtnest
    case (ab7_invars_prtpmp)
      value = token%dtsets(idtset)%prtpmp
    case (ab7_invars_prtposcar)
      value = token%dtsets(idtset)%prtposcar
    case (ab7_invars_prtphdos)
      value = token%dtsets(idtset)%prtphdos
    case (ab7_invars_prtphbands)
      value = token%dtsets(idtset)%prtphbands
    case (ab7_invars_prtphsurf)
      value = token%dtsets(idtset)%prtphsurf
    case (ab7_invars_prtpot)
      value = token%dtsets(idtset)%prtpot
    case (ab7_invars_prtpsps)
      value = token%dtsets(idtset)%prtpsps
    case (ab7_invars_prtspcur)
      value = token%dtsets(idtset)%prtspcur
    case (ab7_invars_prtstm)
      value = token%dtsets(idtset)%prtstm
    case (ab7_invars_prtsuscep)
      value = token%dtsets(idtset)%prtsuscep
    case (ab7_invars_prtvclmb)
      value = token%dtsets(idtset)%prtvclmb
    case (ab7_invars_prtvdw)
      value = token%dtsets(idtset)%prtvdw
    case (ab7_invars_prtvha)
      value = token%dtsets(idtset)%prtvha
    case (ab7_invars_prtvhxc)
      value = token%dtsets(idtset)%prtvhxc
    case (ab7_invars_prtkbff)
      value = token%dtsets(idtset)%prtkbff
    case (ab7_invars_prtvol)
      value = token%dtsets(idtset)%prtvol
    case (ab7_invars_prtvolimg)
      value = token%dtsets(idtset)%prtvolimg
    case (ab7_invars_prtvpsp)
      value = token%dtsets(idtset)%prtvpsp
    case (ab7_invars_prtvxc)
      value = token%dtsets(idtset)%prtvxc
    case (ab7_invars_prtwant)
      value = token%dtsets(idtset)%prtwant
    case (ab7_invars_prtwf)
      value = token%dtsets(idtset)%prtwf
    case (ab7_invars_prtwf_full)
      value = token%dtsets(idtset)%prtwf_full
    case (ab7_invars_prtxml)
      value = token%dtsets(idtset)%prtxml
    case (ab7_invars_prt1dm)
      value = token%dtsets(idtset)%prt1dm
    case (ab7_invars_ptgroupma)
      value = token%dtsets(idtset)%ptgroupma
    case (ab7_invars_qptopt)
      value = token%dtsets(idtset)%qptopt
    case (ab7_invars_random_atpos)
      value = token%dtsets(idtset)%random_atpos
    case (ab7_invars_recgratio)
      value = token%dtsets(idtset)%recgratio
    case (ab7_invars_recnpath)
      value = token%dtsets(idtset)%recnpath
    case (ab7_invars_recnrec)
      value = token%dtsets(idtset)%recnrec
    case (ab7_invars_recptrott)
      value = token%dtsets(idtset)%recptrott
    case (ab7_invars_rectesteg)
      value = token%dtsets(idtset)%rectesteg
    case (ab7_invars_restartxf)
      value = token%dtsets(idtset)%restartxf
    case (ab7_invars_rfasr)
      value = token%dtsets(idtset)%rfasr
    case (ab7_invars_rfddk)
      value = token%dtsets(idtset)%rfddk
    case (ab7_invars_rfelfd)
      value = token%dtsets(idtset)%rfelfd
    case (ab7_invars_rfmagn)
      value = token%dtsets(idtset)%rfmagn
    case (ab7_invars_rfmeth)
      value = token%dtsets(idtset)%rfmeth
    case (ab7_invars_rfphon)
      value = token%dtsets(idtset)%rfphon
    case (ab7_invars_rfstrs)
      value = token%dtsets(idtset)%rfstrs
    case (ab7_invars_rfuser)
      value = token%dtsets(idtset)%rfuser
    case (ab7_invars_rf2_dkdk)
      value = token%dtsets(idtset)%rf2_dkdk
    case (ab7_invars_rf2_dkde)
      value = token%dtsets(idtset)%rf2_dkde
    case (ab7_invars_signperm)
      value = token%dtsets(idtset)%signperm
    case (ab7_invars_slk_rankpp)
      value = token%dtsets(idtset)%slk_rankpp
    case (ab7_invars_smdelta)
      value = token%dtsets(idtset)%smdelta
    case (ab7_invars_spgaxor)
      value = token%dtsets(idtset)%spgaxor
    case (ab7_invars_spgorig)
      value = token%dtsets(idtset)%spgorig
    case (ab7_invars_spgroup)
      value = token%dtsets(idtset)%spgroup
    case (ab7_invars_spmeth)
      value = token%dtsets(idtset)%spmeth
    case (ab7_invars_string_algo)
      value = token%dtsets(idtset)%string_algo
    case (ab7_invars_symmorphi)
      value = token%dtsets(idtset)%symmorphi
    case (ab7_invars_symchi)
      value = token%dtsets(idtset)%symchi
    case (ab7_invars_symsigma)
      value = token%dtsets(idtset)%symsigma
    case (ab7_invars_td_mexcit)
      value = token%dtsets(idtset)%td_mexcit
    case (ab7_invars_tfkinfunc)
      value = token%dtsets(idtset)%tfkinfunc
    case (ab7_invars_tim1rev)
      value = token%dtsets(idtset)%tim1rev
    case (ab7_invars_timopt)
      value = token%dtsets(idtset)%timopt
    case (ab7_invars_tl_nprccg)
      value = token%dtsets(idtset)%tl_nprccg
    case (ab7_invars_ucrpa)
      value = token%dtsets(idtset)%ucrpa
    case (ab7_invars_use_gpu_cuda)
      value = token%dtsets(idtset)%use_gpu_cuda
    case (ab7_invars_usedmatpu)
      value = token%dtsets(idtset)%usedmatpu
    case (ab7_invars_usedmft)
      value = token%dtsets(idtset)%usedmft
    case (ab7_invars_useexexch)
      value = token%dtsets(idtset)%useexexch
    case (ab7_invars_usefock)
      value = token%dtsets(idtset)%usefock
    case (ab7_invars_usekden)
      value = token%dtsets(idtset)%usekden
    case (ab7_invars_use_gemm_nonlop)
      value = token%dtsets(idtset)%use_gemm_nonlop
    case (ab7_invars_use_nonscf_gkk)
      value = token%dtsets(idtset)%use_nonscf_gkk
    case (ab7_invars_usepaw)
      value = token%dtsets(idtset)%usepaw
    case (ab7_invars_usepawu)
      value = token%dtsets(idtset)%usepawu
    case (ab7_invars_usepead)
      value = token%dtsets(idtset)%usepead
    case (ab7_invars_usepotzero)
      value = token%dtsets(idtset)%usepotzero
    case (ab7_invars_userec)
      value = token%dtsets(idtset)%userec
    case (ab7_invars_useria)
      value = token%dtsets(idtset)%useria
    case (ab7_invars_userib)
      value = token%dtsets(idtset)%userib
    case (ab7_invars_useric)
      value = token%dtsets(idtset)%useric
    case (ab7_invars_userid)
      value = token%dtsets(idtset)%userid
    case (ab7_invars_userie)
      value = token%dtsets(idtset)%userie
    case (ab7_invars_usewvl)
      value = token%dtsets(idtset)%usewvl
    case (ab7_invars_usexcnhat_orig)
      value = token%dtsets(idtset)%usexcnhat_orig
    case (ab7_invars_useylm)
      value = token%dtsets(idtset)%useylm
    case (ab7_invars_use_slk)
      value = token%dtsets(idtset)%use_slk
    case (ab7_invars_vacnum)
      value = token%dtsets(idtset)%vacnum
    case (ab7_invars_vdw_nfrag)
      value = token%dtsets(idtset)%vdw_nfrag
    case (ab7_invars_vdw_df_ndpts)
      value = token%dtsets(idtset)%vdw_df_ndpts
    case (ab7_invars_vdw_df_ngpts)
      value = token%dtsets(idtset)%vdw_df_ngpts
    case (ab7_invars_vdw_df_nqpts)
      value = token%dtsets(idtset)%vdw_df_nqpts
    case (ab7_invars_vdw_df_nrpts)
      value = token%dtsets(idtset)%vdw_df_nrpts
    case (ab7_invars_vdw_df_nsmooth)
      value = token%dtsets(idtset)%vdw_df_nsmooth
    case (ab7_invars_vdw_df_tweaks)
      value = token%dtsets(idtset)%vdw_df_tweaks
    case (ab7_invars_vdw_xc)
      value = token%dtsets(idtset)%vdw_xc
    case (ab7_invars_wfoptalg)
      value = token%dtsets(idtset)%wfoptalg
    case (ab7_invars_wfk_task)
      value = token%dtsets(idtset)%wfk_task
    case (ab7_invars_wvl_bigdft_comp)
      value = token%dtsets(idtset)%wvl_bigdft_comp
    case (ab7_invars_wvl_nprccg)
      value = token%dtsets(idtset)%wvl_nprccg
    case (ab7_invars_w90iniprj)
      value = token%dtsets(idtset)%w90iniprj
    case (ab7_invars_w90prtunk)
      value = token%dtsets(idtset)%w90prtunk
    case (ab7_invars_xclevel)
      value = token%dtsets(idtset)%xclevel
    case (ab7_invars_bs_algorithm)
      value = token%dtsets(idtset)%bs_algorithm
    case (ab7_invars_bs_haydock_niter)
      value = token%dtsets(idtset)%bs_haydock_niter
    case (ab7_invars_bs_exchange_term)
      value = token%dtsets(idtset)%bs_exchange_term
    case (ab7_invars_bs_coulomb_term)
      value = token%dtsets(idtset)%bs_coulomb_term
    case (ab7_invars_bs_calctype)
      value = token%dtsets(idtset)%bs_calctype
    case (ab7_invars_bs_coupling)
      value = token%dtsets(idtset)%bs_coupling
    case (ab7_invars_bs_interp_mode)
      value = token%dtsets(idtset)%bs_interp_mode
    case (ab7_invars_bs_interp_prep)
      value = token%dtsets(idtset)%bs_interp_prep
    case (ab7_invars_bs_interp_method)
      value = token%dtsets(idtset)%bs_interp_method
    case (ab7_invars_bs_interp_rl_nb)
      value = token%dtsets(idtset)%bs_interp_rl_nb
    case (ab7_invars_gpu_linalg_limit)
      value = token%dtsets(idtset)%gpu_linalg_limit
    case (ab7_invars_asr)
      value = token%dtsets(idtset)%asr
    case (ab7_invars_dipdip)
      value = token%dtsets(idtset)%dipdip
    case (ab7_invars_chneut)
      value = token%dtsets(idtset)%chneut
    case (ab7_invars_symdynmat)
      value = token%dtsets(idtset)%symdynmat
    case (ab7_invars_ph_freez_disp_addstrain)
      value = token%dtsets(idtset)%ph_freez_disp_addstrain
    case (ab7_invars_ph_freez_disp_option)
      value = token%dtsets(idtset)%ph_freez_disp_option
    case (ab7_invars_ph_freez_disp_nampl)
      value = token%dtsets(idtset)%ph_freez_disp_nampl
    case (ab7_invars_ph_ndivsm)
      value = token%dtsets(idtset)%ph_ndivsm
    case (ab7_invars_ph_nqpath)
      value = token%dtsets(idtset)%ph_nqpath
    case (ab7_invars_ph_nqshift)
      value = token%dtsets(idtset)%ph_nqshift
    case (ab7_invars_ph_)
      value = token%dtsets(idtset)%ph_
    case (ab7_invars_eph_intmeth)
      value = token%dtsets(idtset)%eph_intmeth
    case (ab7_invars_eph_frohlichm)
      value = token%dtsets(idtset)%eph_frohlichm
    case (ab7_invars_eph_transport)
      value = token%dtsets(idtset)%eph_transport
    case (ab7_invars_ph_intmeth)
      value = token%dtsets(idtset)%ph_intmeth
    case (ab7_invars_ndivsm)
      value = token%dtsets(idtset)%ndivsm
    case (ab7_invars_nkpath)
      value = token%dtsets(idtset)%nkpath

    case default
      errno = AB7_ERROR_INVARS_ATT
    end select
  end subroutine ab7_invars_get_integer

  subroutine ab7_invars_get_real(dtsetsId, value, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    real(dp), intent(out) :: value
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB7_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB7_ERROR_INVARS_ID
      return
    end if
    
    errno = AB7_NO_ERROR
    select case (id)
    case (ab7_invars_adpimd_gamma)
      value = token%dtsets(idtset)%adpimd_gamma
    case (ab7_invars_auxc_scal)
      value = token%dtsets(idtset)%auxc_scal
    case (ab7_invars_bmass)
      value = token%dtsets(idtset)%bmass
    case (ab7_invars_boxcutmin)
      value = token%dtsets(idtset)%boxcutmin
    case (ab7_invars_bxctmindg)
      value = token%dtsets(idtset)%bxctmindg
    case (ab7_invars_cd_halfway_freq)
      value = token%dtsets(idtset)%cd_halfway_freq
    case (ab7_invars_cd_max_freq)
      value = token%dtsets(idtset)%cd_max_freq
    case (ab7_invars_charge)
      value = token%dtsets(idtset)%charge
    case (ab7_invars_cpus)
      value = token%dtsets(idtset)%cpus
    case (ab7_invars_ddamp)
      value = token%dtsets(idtset)%ddamp
    case (ab7_invars_dfpt_sciss)
      value = token%dtsets(idtset)%dfpt_sciss
    case (ab7_invars_diecut)
      value = token%dtsets(idtset)%diecut
    case (ab7_invars_diegap)
      value = token%dtsets(idtset)%diegap
    case (ab7_invars_dielam)
      value = token%dtsets(idtset)%dielam
    case (ab7_invars_dielng)
      value = token%dtsets(idtset)%dielng
    case (ab7_invars_diemac)
      value = token%dtsets(idtset)%diemac
    case (ab7_invars_diemix)
      value = token%dtsets(idtset)%diemix
    case (ab7_invars_diemixmag)
      value = token%dtsets(idtset)%diemixmag
    case (ab7_invars_dilatmx)
      value = token%dtsets(idtset)%dilatmx
    case (ab7_invars_dmft_mxsf)
      value = token%dtsets(idtset)%dmft_mxsf
    case (ab7_invars_dmft_tolfreq)
      value = token%dtsets(idtset)%dmft_tolfreq
    case (ab7_invars_dmft_tollc)
      value = token%dtsets(idtset)%dmft_tollc
    case (ab7_invars_dmftqmc_n)
      value = token%dtsets(idtset)%dmftqmc_n
    case (ab7_invars_dosdeltae)
      value = token%dtsets(idtset)%dosdeltae
    case (ab7_invars_dtion)
      value = token%dtsets(idtset)%dtion
    case (ab7_invars_ecut)
      value = token%dtsets(idtset)%ecut
    case (ab7_invars_ecuteps)
      value = token%dtsets(idtset)%ecuteps
    case (ab7_invars_ecutsigx)
      value = token%dtsets(idtset)%ecutsigx
    case (ab7_invars_ecutsm)
      value = token%dtsets(idtset)%ecutsm
    case (ab7_invars_ecutwfn)
      value = token%dtsets(idtset)%ecutwfn
    case (ab7_invars_effmass_free)
      value = token%dtsets(idtset)%effmass_free
    case (ab7_invars_efmas_deg_tol)
      value = token%dtsets(idtset)%efmas_deg_tol
    case (ab7_invars_elph2_imagden)
      value = token%dtsets(idtset)%elph2_imagden
    case (ab7_invars_eshift)
      value = token%dtsets(idtset)%eshift
    case (ab7_invars_esmear)
      value = token%dtsets(idtset)%esmear
    case (ab7_invars_exchmix)
      value = token%dtsets(idtset)%exchmix
    case (ab7_invars_fband)
      value = token%dtsets(idtset)%fband
    case (ab7_invars_fermie_nest)
      value = token%dtsets(idtset)%fermie_nest
    case (ab7_invars_focktoldfe)
      value = token%dtsets(idtset)%focktoldfe
    case (ab7_invars_freqim_alpha)
      value = token%dtsets(idtset)%freqim_alpha
    case (ab7_invars_freqremin)
      value = token%dtsets(idtset)%freqremin
    case (ab7_invars_freqremax)
      value = token%dtsets(idtset)%freqremax
    case (ab7_invars_freqspmin)
      value = token%dtsets(idtset)%freqspmin
    case (ab7_invars_freqspmax)
      value = token%dtsets(idtset)%freqspmax
    case (ab7_invars_friction)
      value = token%dtsets(idtset)%friction
    case (ab7_invars_fxcartfactor)
      value = token%dtsets(idtset)%fxcartfactor
    case (ab7_invars_ga_opt_percent)
      value = token%dtsets(idtset)%ga_opt_percent
    case (ab7_invars_gwencomp)
      value = token%dtsets(idtset)%gwencomp
    case (ab7_invars_gwls_model_parameter)
      value = token%dtsets(idtset)%gwls_model_parameter
    case (ab7_invars_gw_toldfeig)
      value = token%dtsets(idtset)%gw_toldfeig
    case (ab7_invars_hyb_mixing)
      value = token%dtsets(idtset)%hyb_mixing
    case (ab7_invars_hyb_mixing_sr)
      value = token%dtsets(idtset)%hyb_mixing_sr
    case (ab7_invars_hyb_range_dft)
      value = token%dtsets(idtset)%hyb_range_dft
    case (ab7_invars_hyb_range_fock)
      value = token%dtsets(idtset)%hyb_range_fock
    case (ab7_invars_kptnrm)
      value = token%dtsets(idtset)%kptnrm
    case (ab7_invars_kptrlen)
      value = token%dtsets(idtset)%kptrlen
    case (ab7_invars_magcon_lambda)
      value = token%dtsets(idtset)%magcon_lambda
    case (ab7_invars_maxestep)
      value = token%dtsets(idtset)%maxestep
    case (ab7_invars_mbpt_sciss)
      value = token%dtsets(idtset)%mbpt_sciss
    case (ab7_invars_mdf_epsinf)
      value = token%dtsets(idtset)%mdf_epsinf
    case (ab7_invars_mdwall)
      value = token%dtsets(idtset)%mdwall
    case (ab7_invars_mep_mxstep)
      value = token%dtsets(idtset)%mep_mxstep
    case (ab7_invars_nelect)
      value = token%dtsets(idtset)%nelect
    case (ab7_invars_noseinert)
      value = token%dtsets(idtset)%noseinert
    case (ab7_invars_omegasimax)
      value = token%dtsets(idtset)%omegasimax
    case (ab7_invars_omegasrdmax)
      value = token%dtsets(idtset)%omegasrdmax
    case (ab7_invars_pawecutdg)
      value = token%dtsets(idtset)%pawecutdg
    case (ab7_invars_pawovlp)
      value = token%dtsets(idtset)%pawovlp
    case (ab7_invars_pawujrad)
      value = token%dtsets(idtset)%pawujrad
    case (ab7_invars_pawujv)
      value = token%dtsets(idtset)%pawujv
    case (ab7_invars_posocc)
      value = token%dtsets(idtset)%posocc
    case (ab7_invars_postoldfe)
      value = token%dtsets(idtset)%postoldfe
    case (ab7_invars_postoldff)
      value = token%dtsets(idtset)%postoldff
    case (ab7_invars_ppmfrq)
      value = token%dtsets(idtset)%ppmfrq
    case (ab7_invars_pw_unbal_thresh)
      value = token%dtsets(idtset)%pw_unbal_thresh
    case (ab7_invars_ratsph_extra)
      value = token%dtsets(idtset)%ratsph_extra
    case (ab7_invars_recrcut)
      value = token%dtsets(idtset)%recrcut
    case (ab7_invars_recefermi)
      value = token%dtsets(idtset)%recefermi
    case (ab7_invars_rectolden)
      value = token%dtsets(idtset)%rectolden
    case (ab7_invars_rhoqpmix)
      value = token%dtsets(idtset)%rhoqpmix
    case (ab7_invars_rcut)
      value = token%dtsets(idtset)%rcut
    case (ab7_invars_slabwsrad)
      value = token%dtsets(idtset)%slabwsrad
    case (ab7_invars_slabzbeg)
      value = token%dtsets(idtset)%slabzbeg
    case (ab7_invars_slabzend)
      value = token%dtsets(idtset)%slabzend
    case (ab7_invars_spbroad)
      value = token%dtsets(idtset)%spbroad
    case (ab7_invars_spinmagntarget)
      value = token%dtsets(idtset)%spinmagntarget
    case (ab7_invars_spnorbscl)
      value = token%dtsets(idtset)%spnorbscl
    case (ab7_invars_stmbias)
      value = token%dtsets(idtset)%stmbias
    case (ab7_invars_strfact)
      value = token%dtsets(idtset)%strfact
    case (ab7_invars_strprecon)
      value = token%dtsets(idtset)%strprecon
    case (ab7_invars_td_maxene)
      value = token%dtsets(idtset)%td_maxene
    case (ab7_invars_tfw_toldfe)
      value = token%dtsets(idtset)%tfw_toldfe
    case (ab7_invars_tl_radius)
      value = token%dtsets(idtset)%tl_radius
    case (ab7_invars_toldfe)
      value = token%dtsets(idtset)%toldfe
    case (ab7_invars_tolmxde)
      value = token%dtsets(idtset)%tolmxde
    case (ab7_invars_toldff)
      value = token%dtsets(idtset)%toldff
    case (ab7_invars_tolimg)
      value = token%dtsets(idtset)%tolimg
    case (ab7_invars_tolmxf)
      value = token%dtsets(idtset)%tolmxf
    case (ab7_invars_tolrde)
      value = token%dtsets(idtset)%tolrde
    case (ab7_invars_tolrff)
      value = token%dtsets(idtset)%tolrff
    case (ab7_invars_tolsym)
      value = token%dtsets(idtset)%tolsym
    case (ab7_invars_tolvrs)
      value = token%dtsets(idtset)%tolvrs
    case (ab7_invars_tolwfr)
      value = token%dtsets(idtset)%tolwfr
    case (ab7_invars_tphysel)
      value = token%dtsets(idtset)%tphysel
    case (ab7_invars_tsmear)
      value = token%dtsets(idtset)%tsmear
    case (ab7_invars_userra)
      value = token%dtsets(idtset)%userra
    case (ab7_invars_userrb)
      value = token%dtsets(idtset)%userrb
    case (ab7_invars_userrc)
      value = token%dtsets(idtset)%userrc
    case (ab7_invars_userrd)
      value = token%dtsets(idtset)%userrd
    case (ab7_invars_userre)
      value = token%dtsets(idtset)%userre
    case (ab7_invars_vacwidth)
      value = token%dtsets(idtset)%vacwidth
    case (ab7_invars_vdw_tol)
      value = token%dtsets(idtset)%vdw_tol
    case (ab7_invars_vdw_tol_3bt)
      value = token%dtsets(idtset)%vdw_tol_3bt
    case (ab7_invars_vdw_df_acutmin)
      value = token%dtsets(idtset)%vdw_df_acutmin
    case (ab7_invars_vdw_df_aratio)
      value = token%dtsets(idtset)%vdw_df_aratio
    case (ab7_invars_vdw_df_damax)
      value = token%dtsets(idtset)%vdw_df_damax
    case (ab7_invars_vdw_df_damin)
      value = token%dtsets(idtset)%vdw_df_damin
    case (ab7_invars_vdw_df_dcut)
      value = token%dtsets(idtset)%vdw_df_dcut
    case (ab7_invars_vdw_df_dratio)
      value = token%dtsets(idtset)%vdw_df_dratio
    case (ab7_invars_vdw_df_dsoft)
      value = token%dtsets(idtset)%vdw_df_dsoft
    case (ab7_invars_vdw_df_gcut)
      value = token%dtsets(idtset)%vdw_df_gcut
    case (ab7_invars_vdw_df_phisoft)
      value = token%dtsets(idtset)%vdw_df_phisoft
    case (ab7_invars_vdw_df_qcut)
      value = token%dtsets(idtset)%vdw_df_qcut
    case (ab7_invars_vdw_df_qratio)
      value = token%dtsets(idtset)%vdw_df_qratio
    case (ab7_invars_vdw_df_rcut)
      value = token%dtsets(idtset)%vdw_df_rcut
    case (ab7_invars_vdw_df_rsoft)
      value = token%dtsets(idtset)%vdw_df_rsoft
    case (ab7_invars_vdw_df_threshold)
      value = token%dtsets(idtset)%vdw_df_threshold
    case (ab7_invars_vdw_df_tolerance)
      value = token%dtsets(idtset)%vdw_df_tolerance
    case (ab7_invars_vdw_df_zab)
      value = token%dtsets(idtset)%vdw_df_zab
    case (ab7_invars_vis)
      value = token%dtsets(idtset)%vis
    case (ab7_invars_wfmix)
      value = token%dtsets(idtset)%wfmix
    case (ab7_invars_wtq)
      value = token%dtsets(idtset)%wtq
    case (ab7_invars_wvl_hgrid)
      value = token%dtsets(idtset)%wvl_hgrid
    case (ab7_invars_wvl_crmult)
      value = token%dtsets(idtset)%wvl_crmult
    case (ab7_invars_wvl_frmult)
      value = token%dtsets(idtset)%wvl_frmult
    case (ab7_invars_xc_denpos)
      value = token%dtsets(idtset)%xc_denpos
    case (ab7_invars_xc_tb09_c)
      value = token%dtsets(idtset)%xc_tb09_c
    case (ab7_invars_zcut)
      value = token%dtsets(idtset)%zcut
    case (ab7_invars_bs_interp_m3_width)
      value = token%dtsets(idtset)%bs_interp_m3_width
    case (ab7_invars_eph_mustar)
      value = token%dtsets(idtset)%eph_mustar
    case (ab7_invars_eph_extrael)
      value = token%dtsets(idtset)%eph_extrael
    case (ab7_invars_eph_fermie)
      value = token%dtsets(idtset)%eph_fermie
    case (ab7_invars_eph_fsmear)
      value = token%dtsets(idtset)%eph_fsmear
    case (ab7_invars_eph_fsewin)
      value = token%dtsets(idtset)%eph_fsewin
    case (ab7_invars_ph_wstep)
      value = token%dtsets(idtset)%ph_wstep
    case (ab7_invars_ph_smear)
      value = token%dtsets(idtset)%ph_smear

    case default
      errno = AB7_ERROR_INVARS_ATT
    end select
  end subroutine ab7_invars_get_real

  subroutine ab7_invars_get_integer_array(dtsetsId, values, n, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset, n
    integer, intent(out) :: values(n)
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token
    integer :: n_dt

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB7_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB7_ERROR_INVARS_ID
      return
    end if
    
    errno = AB7_NO_ERROR
    select case (id)
    case (ab7_invars_bdberry)
      n_dt = product(shape(token%dtsets(idtset)%bdberry))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bdberry, (/ n_dt /))
      end if
    case (ab7_invars_bravais)
      n_dt = product(shape(token%dtsets(idtset)%bravais))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bravais, (/ n_dt /))
      end if
    case (ab7_invars_cd_subset_freq)
      n_dt = product(shape(token%dtsets(idtset)%cd_subset_freq))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%cd_subset_freq, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert1_atpol)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert1_atpol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert1_atpol, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert1_dir)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert1_dir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert1_dir, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert2_atpol)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert2_atpol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert2_atpol, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert2_dir)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert2_dir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert2_dir, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert3_atpol)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert3_atpol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert3_atpol, (/ n_dt /))
      end if
    case (ab7_invars_d3e_pert3_dir)
      n_dt = product(shape(token%dtsets(idtset)%d3e_pert3_dir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%d3e_pert3_dir, (/ n_dt /))
      end if
    case (ab7_invars_fockdownsampling)
      n_dt = product(shape(token%dtsets(idtset)%fockdownsampling))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%fockdownsampling, (/ n_dt /))
      end if
    case (ab7_invars_jfielddir)
      n_dt = product(shape(token%dtsets(idtset)%jfielddir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%jfielddir, (/ n_dt /))
      end if
    case (ab7_invars_kptrlatt)
      n_dt = product(shape(token%dtsets(idtset)%kptrlatt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptrlatt, (/ n_dt /))
      end if
    case (ab7_invars_kptrlatt_orig)
      n_dt = product(shape(token%dtsets(idtset)%kptrlatt_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptrlatt_orig, (/ n_dt /))
      end if
    case (ab7_invars_qptrlatt)
      n_dt = product(shape(token%dtsets(idtset)%qptrlatt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qptrlatt, (/ n_dt /))
      end if
    case (ab7_invars_ga_rules)
      n_dt = product(shape(token%dtsets(idtset)%ga_rules))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ga_rules, (/ n_dt /))
      end if
    case (ab7_invars_gpu_devices)
      n_dt = product(shape(token%dtsets(idtset)%gpu_devices))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%gpu_devices, (/ n_dt /))
      end if
    case (ab7_invars_ngfft)
      n_dt = product(shape(token%dtsets(idtset)%ngfft))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ngfft, (/ n_dt /))
      end if
    case (ab7_invars_ngfftdg)
      n_dt = product(shape(token%dtsets(idtset)%ngfftdg))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ngfftdg, (/ n_dt /))
      end if
    case (ab7_invars_nloalg)
      n_dt = product(shape(token%dtsets(idtset)%nloalg))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%nloalg, (/ n_dt /))
      end if
    case (ab7_invars_ngkpt)
      n_dt = product(shape(token%dtsets(idtset)%ngkpt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ngkpt, (/ n_dt /))
      end if
    case (ab7_invars_qprtrb)
      n_dt = product(shape(token%dtsets(idtset)%qprtrb))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qprtrb, (/ n_dt /))
      end if
    case (ab7_invars_rfatpol)
      n_dt = product(shape(token%dtsets(idtset)%rfatpol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rfatpol, (/ n_dt /))
      end if
    case (ab7_invars_rfdir)
      n_dt = product(shape(token%dtsets(idtset)%rfdir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rfdir, (/ n_dt /))
      end if
    case (ab7_invars_rf2_pert1_dir)
      n_dt = product(shape(token%dtsets(idtset)%rf2_pert1_dir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf2_pert1_dir, (/ n_dt /))
      end if
    case (ab7_invars_rf2_pert2_dir)
      n_dt = product(shape(token%dtsets(idtset)%rf2_pert2_dir))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rf2_pert2_dir, (/ n_dt /))
      end if
    case (ab7_invars_supercell_latt)
      n_dt = product(shape(token%dtsets(idtset)%supercell_latt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%supercell_latt, (/ n_dt /))
      end if
    case (ab7_invars_ucrpa_bands)
      n_dt = product(shape(token%dtsets(idtset)%ucrpa_bands))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ucrpa_bands, (/ n_dt /))
      end if
    case (ab7_invars_vdw_supercell)
      n_dt = product(shape(token%dtsets(idtset)%vdw_supercell))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vdw_supercell, (/ n_dt /))
      end if
    case (ab7_invars_vdw_typfrag)
      n_dt = product(shape(token%dtsets(idtset)%vdw_typfrag))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vdw_typfrag, (/ n_dt /))
      end if
    case (ab7_invars_wvl_ngauss)
      n_dt = product(shape(token%dtsets(idtset)%wvl_ngauss))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%wvl_ngauss, (/ n_dt /))
      end if
    case (ab7_invars_algalch)
      n_dt = product(shape(token%dtsets(idtset)%algalch))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%algalch, (/ n_dt /))
      end if
    case (ab7_invars_bdgw)
      n_dt = product(shape(token%dtsets(idtset)%bdgw))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bdgw, (/ n_dt /))
      end if
    case (ab7_invars_dynimage)
      n_dt = product(shape(token%dtsets(idtset)%dynimage))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%dynimage, (/ n_dt /))
      end if
    case (ab7_invars_efmas_bands)
      n_dt = product(shape(token%dtsets(idtset)%efmas_bands))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%efmas_bands, (/ n_dt /))
      end if
    case (ab7_invars_iatfix)
      n_dt = product(shape(token%dtsets(idtset)%iatfix))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%iatfix, (/ n_dt /))
      end if
    case (ab7_invars_iatsph)
      n_dt = product(shape(token%dtsets(idtset)%iatsph))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%iatsph, (/ n_dt /))
      end if
    case (ab7_invars_istwfk)
      n_dt = product(shape(token%dtsets(idtset)%istwfk))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%istwfk, (/ n_dt /))
      end if
    case (ab7_invars_kberry)
      n_dt = product(shape(token%dtsets(idtset)%kberry))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kberry, (/ n_dt /))
      end if
    case (ab7_invars_lexexch)
      n_dt = product(shape(token%dtsets(idtset)%lexexch))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%lexexch, (/ n_dt /))
      end if
    case (ab7_invars_ldaminushalf)
      n_dt = product(shape(token%dtsets(idtset)%ldaminushalf))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ldaminushalf, (/ n_dt /))
      end if
    case (ab7_invars_lpawu)
      n_dt = product(shape(token%dtsets(idtset)%lpawu))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%lpawu, (/ n_dt /))
      end if
    case (ab7_invars_nband)
      n_dt = product(shape(token%dtsets(idtset)%nband))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%nband, (/ n_dt /))
      end if
    case (ab7_invars_plowan_iatom)
      n_dt = product(shape(token%dtsets(idtset)%plowan_iatom))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%plowan_iatom, (/ n_dt /))
      end if
    case (ab7_invars_plowan_it)
      n_dt = product(shape(token%dtsets(idtset)%plowan_it))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%plowan_it, (/ n_dt /))
      end if
    case (ab7_invars_plowan_lcalc)
      n_dt = product(shape(token%dtsets(idtset)%plowan_lcalc))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%plowan_lcalc, (/ n_dt /))
      end if
    case (ab7_invars_plowan_nbl)
      n_dt = product(shape(token%dtsets(idtset)%plowan_nbl))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%plowan_nbl, (/ n_dt /))
      end if
    case (ab7_invars_plowan_projcalc)
      n_dt = product(shape(token%dtsets(idtset)%plowan_projcalc))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%plowan_projcalc, (/ n_dt /))
      end if
    case (ab7_invars_prtatlist)
      n_dt = product(shape(token%dtsets(idtset)%prtatlist))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%prtatlist, (/ n_dt /))
      end if
    case (ab7_invars_so_psp)
      n_dt = product(shape(token%dtsets(idtset)%so_psp))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%so_psp, (/ n_dt /))
      end if
    case (ab7_invars_symafm)
      n_dt = product(shape(token%dtsets(idtset)%symafm))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%symafm, (/ n_dt /))
      end if
    case (ab7_invars_symrel)
      n_dt = product(shape(token%dtsets(idtset)%symrel))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%symrel, (/ n_dt /))
      end if
    case (ab7_invars_typat)
      n_dt = product(shape(token%dtsets(idtset)%typat))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%typat, (/ n_dt /))
      end if
    case (ab7_invars_bs_interp_kmult)
      n_dt = product(shape(token%dtsets(idtset)%bs_interp_kmult))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_interp_kmult, (/ n_dt /))
      end if
    case (ab7_invars_bs_loband)
      n_dt = product(shape(token%dtsets(idtset)%bs_loband))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_loband, (/ n_dt /))
      end if
    case (ab7_invars_ph_ngqpt)
      n_dt = product(shape(token%dtsets(idtset)%ph_ngqpt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ph_ngqpt, (/ n_dt /))
      end if
    case (ab7_invars_eph_ngqpt_fine)
      n_dt = product(shape(token%dtsets(idtset)%eph_ngqpt_fine))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%eph_ngqpt_fine, (/ n_dt /))
      end if
    case (ab7_invars_ddb_ngqpt)
      n_dt = product(shape(token%dtsets(idtset)%ddb_ngqpt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ddb_ngqpt, (/ n_dt /))
      end if

    case default
      errno = AB7_ERROR_INVARS_ATT
    end select
  end subroutine ab7_invars_get_integer_array

  subroutine ab7_invars_get_real_array(dtsetsId, values, n, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset, n
    real(dp), intent(out) :: values(n)
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token
    integer :: n_dt

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB7_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB7_ERROR_INVARS_ID
      return
    end if
    
    errno = AB7_NO_ERROR
    select case (id)
    case (ab7_invars_boxcenter)
      n_dt = product(shape(token%dtsets(idtset)%boxcenter))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%boxcenter, (/ n_dt /))
      end if
    case (ab7_invars_bfield)
      n_dt = product(shape(token%dtsets(idtset)%bfield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bfield, (/ n_dt /))
      end if
    case (ab7_invars_dfield)
      n_dt = product(shape(token%dtsets(idtset)%dfield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%dfield, (/ n_dt /))
      end if
    case (ab7_invars_efield)
      n_dt = product(shape(token%dtsets(idtset)%efield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%efield, (/ n_dt /))
      end if
    case (ab7_invars_genafm)
      n_dt = product(shape(token%dtsets(idtset)%genafm))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%genafm, (/ n_dt /))
      end if
    case (ab7_invars_goprecprm)
      n_dt = product(shape(token%dtsets(idtset)%goprecprm))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%goprecprm, (/ n_dt /))
      end if
    case (ab7_invars_neb_spring)
      n_dt = product(shape(token%dtsets(idtset)%neb_spring))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%neb_spring, (/ n_dt /))
      end if
    case (ab7_invars_pol)
      n_dt = product(shape(token%dtsets(idtset)%pol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%pol, (/ n_dt /))
      end if
    case (ab7_invars_polcen)
      n_dt = product(shape(token%dtsets(idtset)%polcen))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%polcen, (/ n_dt /))
      end if
    case (ab7_invars_pvelmax)
      n_dt = product(shape(token%dtsets(idtset)%pvelmax))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%pvelmax, (/ n_dt /))
      end if
    case (ab7_invars_qptn)
      n_dt = product(shape(token%dtsets(idtset)%qptn))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qptn, (/ n_dt /))
      end if
    case (ab7_invars_red_efield)
      n_dt = product(shape(token%dtsets(idtset)%red_efield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%red_efield, (/ n_dt /))
      end if
    case (ab7_invars_red_dfield)
      n_dt = product(shape(token%dtsets(idtset)%red_dfield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%red_dfield, (/ n_dt /))
      end if
    case (ab7_invars_red_efieldbar)
      n_dt = product(shape(token%dtsets(idtset)%red_efieldbar))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%red_efieldbar, (/ n_dt /))
      end if
    case (ab7_invars_strtarget)
      n_dt = product(shape(token%dtsets(idtset)%strtarget))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%strtarget, (/ n_dt /))
      end if
    case (ab7_invars_ucrpa_window)
      n_dt = product(shape(token%dtsets(idtset)%ucrpa_window))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ucrpa_window, (/ n_dt /))
      end if
    case (ab7_invars_vcutgeo)
      n_dt = product(shape(token%dtsets(idtset)%vcutgeo))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vcutgeo, (/ n_dt /))
      end if
    case (ab7_invars_vprtrb)
      n_dt = product(shape(token%dtsets(idtset)%vprtrb))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vprtrb, (/ n_dt /))
      end if
    case (ab7_invars_zeemanfield)
      n_dt = product(shape(token%dtsets(idtset)%zeemanfield))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%zeemanfield, (/ n_dt /))
      end if
    case (ab7_invars_mdtemp)
      n_dt = product(shape(token%dtsets(idtset)%mdtemp))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%mdtemp, (/ n_dt /))
      end if
    case (ab7_invars_acell_orig)
      n_dt = product(shape(token%dtsets(idtset)%acell_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%acell_orig, (/ n_dt /))
      end if
    case (ab7_invars_amu_orig)
      n_dt = product(shape(token%dtsets(idtset)%amu_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%amu_orig, (/ n_dt /))
      end if
    case (ab7_invars_atvshift)
      n_dt = product(shape(token%dtsets(idtset)%atvshift))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%atvshift, (/ n_dt /))
      end if
    case (ab7_invars_cd_imfrqs)
      n_dt = product(shape(token%dtsets(idtset)%cd_imfrqs))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%cd_imfrqs, (/ n_dt /))
      end if
    case (ab7_invars_chempot)
      n_dt = product(shape(token%dtsets(idtset)%chempot))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%chempot, (/ n_dt /))
      end if
    case (ab7_invars_corecs)
      n_dt = product(shape(token%dtsets(idtset)%corecs))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%corecs, (/ n_dt /))
      end if
    case (ab7_invars_densty)
      n_dt = product(shape(token%dtsets(idtset)%densty))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%densty, (/ n_dt /))
      end if
    case (ab7_invars_dmatpawu)
      n_dt = product(shape(token%dtsets(idtset)%dmatpawu))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%dmatpawu, (/ n_dt /))
      end if
    case (ab7_invars_efmas_dirs)
      n_dt = product(shape(token%dtsets(idtset)%efmas_dirs))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%efmas_dirs, (/ n_dt /))
      end if
    case (ab7_invars_f4of2_sla)
      n_dt = product(shape(token%dtsets(idtset)%f4of2_sla))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%f4of2_sla, (/ n_dt /))
      end if
    case (ab7_invars_f6of2_sla)
      n_dt = product(shape(token%dtsets(idtset)%f6of2_sla))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%f6of2_sla, (/ n_dt /))
      end if
    case (ab7_invars_gw_qlwl)
      n_dt = product(shape(token%dtsets(idtset)%gw_qlwl))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%gw_qlwl, (/ n_dt /))
      end if
    case (ab7_invars_gw_freqsp)
      n_dt = product(shape(token%dtsets(idtset)%gw_freqsp))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%gw_freqsp, (/ n_dt /))
      end if
    case (ab7_invars_gwls_list_proj_freq)
      n_dt = product(shape(token%dtsets(idtset)%gwls_list_proj_freq))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%gwls_list_proj_freq, (/ n_dt /))
      end if
    case (ab7_invars_jpawu)
      n_dt = product(shape(token%dtsets(idtset)%jpawu))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%jpawu, (/ n_dt /))
      end if
    case (ab7_invars_kpt)
      n_dt = product(shape(token%dtsets(idtset)%kpt))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kpt, (/ n_dt /))
      end if
    case (ab7_invars_kptgw)
      n_dt = product(shape(token%dtsets(idtset)%kptgw))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptgw, (/ n_dt /))
      end if
    case (ab7_invars_kptns)
      n_dt = product(shape(token%dtsets(idtset)%kptns))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptns, (/ n_dt /))
      end if
    case (ab7_invars_kptns_hf)
      n_dt = product(shape(token%dtsets(idtset)%kptns_hf))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptns_hf, (/ n_dt /))
      end if
    case (ab7_invars_mixalch_orig)
      n_dt = product(shape(token%dtsets(idtset)%mixalch_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%mixalch_orig, (/ n_dt /))
      end if
    case (ab7_invars_mixesimgf)
      n_dt = product(shape(token%dtsets(idtset)%mixesimgf))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%mixesimgf, (/ n_dt /))
      end if
    case (ab7_invars_nucdipmom)
      n_dt = product(shape(token%dtsets(idtset)%nucdipmom))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%nucdipmom, (/ n_dt /))
      end if
    case (ab7_invars_occ_orig)
      n_dt = product(shape(token%dtsets(idtset)%occ_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%occ_orig, (/ n_dt /))
      end if
    case (ab7_invars_pimass)
      n_dt = product(shape(token%dtsets(idtset)%pimass))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%pimass, (/ n_dt /))
      end if
    case (ab7_invars_ptcharge)
      n_dt = product(shape(token%dtsets(idtset)%ptcharge))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ptcharge, (/ n_dt /))
      end if
    case (ab7_invars_qmass)
      n_dt = product(shape(token%dtsets(idtset)%qmass))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qmass, (/ n_dt /))
      end if
    case (ab7_invars_qptdm)
      n_dt = product(shape(token%dtsets(idtset)%qptdm))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%qptdm, (/ n_dt /))
      end if
    case (ab7_invars_quadmom)
      n_dt = product(shape(token%dtsets(idtset)%quadmom))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%quadmom, (/ n_dt /))
      end if
    case (ab7_invars_ratsph)
      n_dt = product(shape(token%dtsets(idtset)%ratsph))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ratsph, (/ n_dt /))
      end if
    case (ab7_invars_rprim_orig)
      n_dt = product(shape(token%dtsets(idtset)%rprim_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rprim_orig, (/ n_dt /))
      end if
    case (ab7_invars_rprimd_orig)
      n_dt = product(shape(token%dtsets(idtset)%rprimd_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%rprimd_orig, (/ n_dt /))
      end if
    case (ab7_invars_shiftk)
      n_dt = product(shape(token%dtsets(idtset)%shiftk))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%shiftk, (/ n_dt /))
      end if
    case (ab7_invars_shiftk_orig)
      n_dt = product(shape(token%dtsets(idtset)%shiftk_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%shiftk_orig, (/ n_dt /))
      end if
    case (ab7_invars_spinat)
      n_dt = product(shape(token%dtsets(idtset)%spinat))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%spinat, (/ n_dt /))
      end if
    case (ab7_invars_tnons)
      n_dt = product(shape(token%dtsets(idtset)%tnons))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%tnons, (/ n_dt /))
      end if
    case (ab7_invars_upawu)
      n_dt = product(shape(token%dtsets(idtset)%upawu))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%upawu, (/ n_dt /))
      end if
    case (ab7_invars_vel_cell_orig)
      n_dt = product(shape(token%dtsets(idtset)%vel_cell_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vel_cell_orig, (/ n_dt /))
      end if
    case (ab7_invars_vel_orig)
      n_dt = product(shape(token%dtsets(idtset)%vel_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%vel_orig, (/ n_dt /))
      end if
    case (ab7_invars_wtatcon)
      n_dt = product(shape(token%dtsets(idtset)%wtatcon))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%wtatcon, (/ n_dt /))
      end if
    case (ab7_invars_wtk)
      n_dt = product(shape(token%dtsets(idtset)%wtk))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%wtk, (/ n_dt /))
      end if
    case (ab7_invars_xred_orig)
      n_dt = product(shape(token%dtsets(idtset)%xred_orig))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%xred_orig, (/ n_dt /))
      end if
    case (ab7_invars_xredsph_extra)
      n_dt = product(shape(token%dtsets(idtset)%xredsph_extra))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%xredsph_extra, (/ n_dt /))
      end if
    case (ab7_invars_ziontypat)
      n_dt = product(shape(token%dtsets(idtset)%ziontypat))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ziontypat, (/ n_dt /))
      end if
    case (ab7_invars_znucl)
      n_dt = product(shape(token%dtsets(idtset)%znucl))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%znucl, (/ n_dt /))
      end if
    case (ab7_invars_bs_haydock_tol)
      n_dt = product(shape(token%dtsets(idtset)%bs_haydock_tol))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_haydock_tol, (/ n_dt /))
      end if
    case (ab7_invars_bs_eh_cutoff)
      n_dt = product(shape(token%dtsets(idtset)%bs_eh_cutoff))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_eh_cutoff, (/ n_dt /))
      end if
    case (ab7_invars_bs_freq_mesh)
      n_dt = product(shape(token%dtsets(idtset)%bs_freq_mesh))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%bs_freq_mesh, (/ n_dt /))
      end if
    case (ab7_invars_ph_freez_disp_ampl)
      n_dt = product(shape(token%dtsets(idtset)%ph_freez_disp_ampl))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ph_freez_disp_ampl, (/ n_dt /))
      end if
    case (ab7_invars_ph_qshift)
      n_dt = product(shape(token%dtsets(idtset)%ph_qshift))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ph_qshift, (/ n_dt /))
      end if
    case (ab7_invars_ph_qpath)
      n_dt = product(shape(token%dtsets(idtset)%ph_qpath))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ph_qpath, (/ n_dt /))
      end if
    case (ab7_invars_ddb_shiftq)
      n_dt = product(shape(token%dtsets(idtset)%ddb_shiftq))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%ddb_shiftq, (/ n_dt /))
      end if
    case (ab7_invars_einterp)
      n_dt = product(shape(token%dtsets(idtset)%einterp))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%einterp, (/ n_dt /))
      end if
    case (ab7_invars_kptbounds)
      n_dt = product(shape(token%dtsets(idtset)%kptbounds))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%kptbounds, (/ n_dt /))
      end if
    case (ab7_invars_tmesh)
      n_dt = product(shape(token%dtsets(idtset)%tmesh))
      if (n_dt /= n) then
        errno = AB7_ERROR_INVARS_SIZE
      else
        values = reshape(token%dtsets(idtset)%tmesh, (/ n_dt /))
      end if

    case default
      errno = AB7_ERROR_INVARS_ATT
    end select
  end subroutine ab7_invars_get_real_array

  subroutine ab7_invars_get_shape(dtsetsId, dims, ndim, id, idtset, errno)
    integer, intent(in) :: dtsetsId
    integer, intent(in) :: id, idtset
    integer, intent(out) :: dims(7), ndim
    integer, intent(out) :: errno

    type(dtsets_list), pointer :: token

    call get_token(token, dtsetsId)
    if (.not. associated(token)) then
      errno = AB7_ERROR_OBJ
      return
    end if
    if (idtset < 0 .or. idtset > size(token%dtsets)) then
      errno = AB7_ERROR_INVARS_ID
      return
    end if
    
    errno = AB7_NO_ERROR
    select case (id)
    case (ab7_invars_bdberry)
      ndim = size(shape(token%dtsets(idtset)%bdberry))
      dims(1:ndim) = shape(token%dtsets(idtset)%bdberry)
    case (ab7_invars_bravais)
      ndim = size(shape(token%dtsets(idtset)%bravais))
      dims(1:ndim) = shape(token%dtsets(idtset)%bravais)
    case (ab7_invars_cd_subset_freq)
      ndim = size(shape(token%dtsets(idtset)%cd_subset_freq))
      dims(1:ndim) = shape(token%dtsets(idtset)%cd_subset_freq)
    case (ab7_invars_d3e_pert1_atpol)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert1_atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert1_atpol)
    case (ab7_invars_d3e_pert1_dir)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert1_dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert1_dir)
    case (ab7_invars_d3e_pert2_atpol)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert2_atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert2_atpol)
    case (ab7_invars_d3e_pert2_dir)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert2_dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert2_dir)
    case (ab7_invars_d3e_pert3_atpol)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert3_atpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert3_atpol)
    case (ab7_invars_d3e_pert3_dir)
      ndim = size(shape(token%dtsets(idtset)%d3e_pert3_dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%d3e_pert3_dir)
    case (ab7_invars_fockdownsampling)
      ndim = size(shape(token%dtsets(idtset)%fockdownsampling))
      dims(1:ndim) = shape(token%dtsets(idtset)%fockdownsampling)
    case (ab7_invars_jfielddir)
      ndim = size(shape(token%dtsets(idtset)%jfielddir))
      dims(1:ndim) = shape(token%dtsets(idtset)%jfielddir)
    case (ab7_invars_kptrlatt)
      ndim = size(shape(token%dtsets(idtset)%kptrlatt))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptrlatt)
    case (ab7_invars_kptrlatt_orig)
      ndim = size(shape(token%dtsets(idtset)%kptrlatt_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptrlatt_orig)
    case (ab7_invars_qptrlatt)
      ndim = size(shape(token%dtsets(idtset)%qptrlatt))
      dims(1:ndim) = shape(token%dtsets(idtset)%qptrlatt)
    case (ab7_invars_ga_rules)
      ndim = size(shape(token%dtsets(idtset)%ga_rules))
      dims(1:ndim) = shape(token%dtsets(idtset)%ga_rules)
    case (ab7_invars_gpu_devices)
      ndim = size(shape(token%dtsets(idtset)%gpu_devices))
      dims(1:ndim) = shape(token%dtsets(idtset)%gpu_devices)
    case (ab7_invars_ngfft)
      ndim = size(shape(token%dtsets(idtset)%ngfft))
      dims(1:ndim) = shape(token%dtsets(idtset)%ngfft)
    case (ab7_invars_ngfftdg)
      ndim = size(shape(token%dtsets(idtset)%ngfftdg))
      dims(1:ndim) = shape(token%dtsets(idtset)%ngfftdg)
    case (ab7_invars_nloalg)
      ndim = size(shape(token%dtsets(idtset)%nloalg))
      dims(1:ndim) = shape(token%dtsets(idtset)%nloalg)
    case (ab7_invars_ngkpt)
      ndim = size(shape(token%dtsets(idtset)%ngkpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%ngkpt)
    case (ab7_invars_qprtrb)
      ndim = size(shape(token%dtsets(idtset)%qprtrb))
      dims(1:ndim) = shape(token%dtsets(idtset)%qprtrb)
    case (ab7_invars_rfatpol)
      ndim = size(shape(token%dtsets(idtset)%rfatpol))
      dims(1:ndim) = shape(token%dtsets(idtset)%rfatpol)
    case (ab7_invars_rfdir)
      ndim = size(shape(token%dtsets(idtset)%rfdir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rfdir)
    case (ab7_invars_rf2_pert1_dir)
      ndim = size(shape(token%dtsets(idtset)%rf2_pert1_dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf2_pert1_dir)
    case (ab7_invars_rf2_pert2_dir)
      ndim = size(shape(token%dtsets(idtset)%rf2_pert2_dir))
      dims(1:ndim) = shape(token%dtsets(idtset)%rf2_pert2_dir)
    case (ab7_invars_supercell_latt)
      ndim = size(shape(token%dtsets(idtset)%supercell_latt))
      dims(1:ndim) = shape(token%dtsets(idtset)%supercell_latt)
    case (ab7_invars_ucrpa_bands)
      ndim = size(shape(token%dtsets(idtset)%ucrpa_bands))
      dims(1:ndim) = shape(token%dtsets(idtset)%ucrpa_bands)
    case (ab7_invars_vdw_supercell)
      ndim = size(shape(token%dtsets(idtset)%vdw_supercell))
      dims(1:ndim) = shape(token%dtsets(idtset)%vdw_supercell)
    case (ab7_invars_vdw_typfrag)
      ndim = size(shape(token%dtsets(idtset)%vdw_typfrag))
      dims(1:ndim) = shape(token%dtsets(idtset)%vdw_typfrag)
    case (ab7_invars_wvl_ngauss)
      ndim = size(shape(token%dtsets(idtset)%wvl_ngauss))
      dims(1:ndim) = shape(token%dtsets(idtset)%wvl_ngauss)
    case (ab7_invars_algalch)
      ndim = size(shape(token%dtsets(idtset)%algalch))
      dims(1:ndim) = shape(token%dtsets(idtset)%algalch)
    case (ab7_invars_bdgw)
      ndim = size(shape(token%dtsets(idtset)%bdgw))
      dims(1:ndim) = shape(token%dtsets(idtset)%bdgw)
    case (ab7_invars_dynimage)
      ndim = size(shape(token%dtsets(idtset)%dynimage))
      dims(1:ndim) = shape(token%dtsets(idtset)%dynimage)
    case (ab7_invars_efmas_bands)
      ndim = size(shape(token%dtsets(idtset)%efmas_bands))
      dims(1:ndim) = shape(token%dtsets(idtset)%efmas_bands)
    case (ab7_invars_iatfix)
      ndim = size(shape(token%dtsets(idtset)%iatfix))
      dims(1:ndim) = shape(token%dtsets(idtset)%iatfix)
    case (ab7_invars_iatsph)
      ndim = size(shape(token%dtsets(idtset)%iatsph))
      dims(1:ndim) = shape(token%dtsets(idtset)%iatsph)
    case (ab7_invars_istwfk)
      ndim = size(shape(token%dtsets(idtset)%istwfk))
      dims(1:ndim) = shape(token%dtsets(idtset)%istwfk)
    case (ab7_invars_kberry)
      ndim = size(shape(token%dtsets(idtset)%kberry))
      dims(1:ndim) = shape(token%dtsets(idtset)%kberry)
    case (ab7_invars_lexexch)
      ndim = size(shape(token%dtsets(idtset)%lexexch))
      dims(1:ndim) = shape(token%dtsets(idtset)%lexexch)
    case (ab7_invars_ldaminushalf)
      ndim = size(shape(token%dtsets(idtset)%ldaminushalf))
      dims(1:ndim) = shape(token%dtsets(idtset)%ldaminushalf)
    case (ab7_invars_lpawu)
      ndim = size(shape(token%dtsets(idtset)%lpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%lpawu)
    case (ab7_invars_nband)
      ndim = size(shape(token%dtsets(idtset)%nband))
      dims(1:ndim) = shape(token%dtsets(idtset)%nband)
    case (ab7_invars_plowan_iatom)
      ndim = size(shape(token%dtsets(idtset)%plowan_iatom))
      dims(1:ndim) = shape(token%dtsets(idtset)%plowan_iatom)
    case (ab7_invars_plowan_it)
      ndim = size(shape(token%dtsets(idtset)%plowan_it))
      dims(1:ndim) = shape(token%dtsets(idtset)%plowan_it)
    case (ab7_invars_plowan_lcalc)
      ndim = size(shape(token%dtsets(idtset)%plowan_lcalc))
      dims(1:ndim) = shape(token%dtsets(idtset)%plowan_lcalc)
    case (ab7_invars_plowan_nbl)
      ndim = size(shape(token%dtsets(idtset)%plowan_nbl))
      dims(1:ndim) = shape(token%dtsets(idtset)%plowan_nbl)
    case (ab7_invars_plowan_projcalc)
      ndim = size(shape(token%dtsets(idtset)%plowan_projcalc))
      dims(1:ndim) = shape(token%dtsets(idtset)%plowan_projcalc)
    case (ab7_invars_prtatlist)
      ndim = size(shape(token%dtsets(idtset)%prtatlist))
      dims(1:ndim) = shape(token%dtsets(idtset)%prtatlist)
    case (ab7_invars_so_psp)
      ndim = size(shape(token%dtsets(idtset)%so_psp))
      dims(1:ndim) = shape(token%dtsets(idtset)%so_psp)
    case (ab7_invars_symafm)
      ndim = size(shape(token%dtsets(idtset)%symafm))
      dims(1:ndim) = shape(token%dtsets(idtset)%symafm)
    case (ab7_invars_symrel)
      ndim = size(shape(token%dtsets(idtset)%symrel))
      dims(1:ndim) = shape(token%dtsets(idtset)%symrel)
    case (ab7_invars_typat)
      ndim = size(shape(token%dtsets(idtset)%typat))
      dims(1:ndim) = shape(token%dtsets(idtset)%typat)
    case (ab7_invars_boxcenter)
      ndim = size(shape(token%dtsets(idtset)%boxcenter))
      dims(1:ndim) = shape(token%dtsets(idtset)%boxcenter)
    case (ab7_invars_bfield)
      ndim = size(shape(token%dtsets(idtset)%bfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%bfield)
    case (ab7_invars_dfield)
      ndim = size(shape(token%dtsets(idtset)%dfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%dfield)
    case (ab7_invars_efield)
      ndim = size(shape(token%dtsets(idtset)%efield))
      dims(1:ndim) = shape(token%dtsets(idtset)%efield)
    case (ab7_invars_genafm)
      ndim = size(shape(token%dtsets(idtset)%genafm))
      dims(1:ndim) = shape(token%dtsets(idtset)%genafm)
    case (ab7_invars_goprecprm)
      ndim = size(shape(token%dtsets(idtset)%goprecprm))
      dims(1:ndim) = shape(token%dtsets(idtset)%goprecprm)
    case (ab7_invars_neb_spring)
      ndim = size(shape(token%dtsets(idtset)%neb_spring))
      dims(1:ndim) = shape(token%dtsets(idtset)%neb_spring)
    case (ab7_invars_pol)
      ndim = size(shape(token%dtsets(idtset)%pol))
      dims(1:ndim) = shape(token%dtsets(idtset)%pol)
    case (ab7_invars_polcen)
      ndim = size(shape(token%dtsets(idtset)%polcen))
      dims(1:ndim) = shape(token%dtsets(idtset)%polcen)
    case (ab7_invars_pvelmax)
      ndim = size(shape(token%dtsets(idtset)%pvelmax))
      dims(1:ndim) = shape(token%dtsets(idtset)%pvelmax)
    case (ab7_invars_qptn)
      ndim = size(shape(token%dtsets(idtset)%qptn))
      dims(1:ndim) = shape(token%dtsets(idtset)%qptn)
    case (ab7_invars_red_efield)
      ndim = size(shape(token%dtsets(idtset)%red_efield))
      dims(1:ndim) = shape(token%dtsets(idtset)%red_efield)
    case (ab7_invars_red_dfield)
      ndim = size(shape(token%dtsets(idtset)%red_dfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%red_dfield)
    case (ab7_invars_red_efieldbar)
      ndim = size(shape(token%dtsets(idtset)%red_efieldbar))
      dims(1:ndim) = shape(token%dtsets(idtset)%red_efieldbar)
    case (ab7_invars_strtarget)
      ndim = size(shape(token%dtsets(idtset)%strtarget))
      dims(1:ndim) = shape(token%dtsets(idtset)%strtarget)
    case (ab7_invars_ucrpa_window)
      ndim = size(shape(token%dtsets(idtset)%ucrpa_window))
      dims(1:ndim) = shape(token%dtsets(idtset)%ucrpa_window)
    case (ab7_invars_vcutgeo)
      ndim = size(shape(token%dtsets(idtset)%vcutgeo))
      dims(1:ndim) = shape(token%dtsets(idtset)%vcutgeo)
    case (ab7_invars_vprtrb)
      ndim = size(shape(token%dtsets(idtset)%vprtrb))
      dims(1:ndim) = shape(token%dtsets(idtset)%vprtrb)
    case (ab7_invars_zeemanfield)
      ndim = size(shape(token%dtsets(idtset)%zeemanfield))
      dims(1:ndim) = shape(token%dtsets(idtset)%zeemanfield)
    case (ab7_invars_mdtemp)
      ndim = size(shape(token%dtsets(idtset)%mdtemp))
      dims(1:ndim) = shape(token%dtsets(idtset)%mdtemp)
    case (ab7_invars_acell_orig)
      ndim = size(shape(token%dtsets(idtset)%acell_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%acell_orig)
    case (ab7_invars_amu_orig)
      ndim = size(shape(token%dtsets(idtset)%amu_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%amu_orig)
    case (ab7_invars_atvshift)
      ndim = size(shape(token%dtsets(idtset)%atvshift))
      dims(1:ndim) = shape(token%dtsets(idtset)%atvshift)
    case (ab7_invars_cd_imfrqs)
      ndim = size(shape(token%dtsets(idtset)%cd_imfrqs))
      dims(1:ndim) = shape(token%dtsets(idtset)%cd_imfrqs)
    case (ab7_invars_chempot)
      ndim = size(shape(token%dtsets(idtset)%chempot))
      dims(1:ndim) = shape(token%dtsets(idtset)%chempot)
    case (ab7_invars_corecs)
      ndim = size(shape(token%dtsets(idtset)%corecs))
      dims(1:ndim) = shape(token%dtsets(idtset)%corecs)
    case (ab7_invars_densty)
      ndim = size(shape(token%dtsets(idtset)%densty))
      dims(1:ndim) = shape(token%dtsets(idtset)%densty)
    case (ab7_invars_dmatpawu)
      ndim = size(shape(token%dtsets(idtset)%dmatpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%dmatpawu)
    case (ab7_invars_efmas_dirs)
      ndim = size(shape(token%dtsets(idtset)%efmas_dirs))
      dims(1:ndim) = shape(token%dtsets(idtset)%efmas_dirs)
    case (ab7_invars_f4of2_sla)
      ndim = size(shape(token%dtsets(idtset)%f4of2_sla))
      dims(1:ndim) = shape(token%dtsets(idtset)%f4of2_sla)
    case (ab7_invars_f6of2_sla)
      ndim = size(shape(token%dtsets(idtset)%f6of2_sla))
      dims(1:ndim) = shape(token%dtsets(idtset)%f6of2_sla)
    case (ab7_invars_gw_qlwl)
      ndim = size(shape(token%dtsets(idtset)%gw_qlwl))
      dims(1:ndim) = shape(token%dtsets(idtset)%gw_qlwl)
    case (ab7_invars_gw_freqsp)
      ndim = size(shape(token%dtsets(idtset)%gw_freqsp))
      dims(1:ndim) = shape(token%dtsets(idtset)%gw_freqsp)
    case (ab7_invars_gwls_list_proj_freq)
      ndim = size(shape(token%dtsets(idtset)%gwls_list_proj_freq))
      dims(1:ndim) = shape(token%dtsets(idtset)%gwls_list_proj_freq)
    case (ab7_invars_jpawu)
      ndim = size(shape(token%dtsets(idtset)%jpawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%jpawu)
    case (ab7_invars_kpt)
      ndim = size(shape(token%dtsets(idtset)%kpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%kpt)
    case (ab7_invars_kptgw)
      ndim = size(shape(token%dtsets(idtset)%kptgw))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptgw)
    case (ab7_invars_kptns)
      ndim = size(shape(token%dtsets(idtset)%kptns))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptns)
    case (ab7_invars_kptns_hf)
      ndim = size(shape(token%dtsets(idtset)%kptns_hf))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptns_hf)
    case (ab7_invars_mixalch_orig)
      ndim = size(shape(token%dtsets(idtset)%mixalch_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%mixalch_orig)
    case (ab7_invars_mixesimgf)
      ndim = size(shape(token%dtsets(idtset)%mixesimgf))
      dims(1:ndim) = shape(token%dtsets(idtset)%mixesimgf)
    case (ab7_invars_nucdipmom)
      ndim = size(shape(token%dtsets(idtset)%nucdipmom))
      dims(1:ndim) = shape(token%dtsets(idtset)%nucdipmom)
    case (ab7_invars_occ_orig)
      ndim = size(shape(token%dtsets(idtset)%occ_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%occ_orig)
    case (ab7_invars_pimass)
      ndim = size(shape(token%dtsets(idtset)%pimass))
      dims(1:ndim) = shape(token%dtsets(idtset)%pimass)
    case (ab7_invars_ptcharge)
      ndim = size(shape(token%dtsets(idtset)%ptcharge))
      dims(1:ndim) = shape(token%dtsets(idtset)%ptcharge)
    case (ab7_invars_qmass)
      ndim = size(shape(token%dtsets(idtset)%qmass))
      dims(1:ndim) = shape(token%dtsets(idtset)%qmass)
    case (ab7_invars_qptdm)
      ndim = size(shape(token%dtsets(idtset)%qptdm))
      dims(1:ndim) = shape(token%dtsets(idtset)%qptdm)
    case (ab7_invars_quadmom)
      ndim = size(shape(token%dtsets(idtset)%quadmom))
      dims(1:ndim) = shape(token%dtsets(idtset)%quadmom)
    case (ab7_invars_ratsph)
      ndim = size(shape(token%dtsets(idtset)%ratsph))
      dims(1:ndim) = shape(token%dtsets(idtset)%ratsph)
    case (ab7_invars_rprim_orig)
      ndim = size(shape(token%dtsets(idtset)%rprim_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%rprim_orig)
    case (ab7_invars_rprimd_orig)
      ndim = size(shape(token%dtsets(idtset)%rprimd_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%rprimd_orig)
    case (ab7_invars_shiftk)
      ndim = size(shape(token%dtsets(idtset)%shiftk))
      dims(1:ndim) = shape(token%dtsets(idtset)%shiftk)
    case (ab7_invars_shiftk_orig)
      ndim = size(shape(token%dtsets(idtset)%shiftk_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%shiftk_orig)
    case (ab7_invars_spinat)
      ndim = size(shape(token%dtsets(idtset)%spinat))
      dims(1:ndim) = shape(token%dtsets(idtset)%spinat)
    case (ab7_invars_tnons)
      ndim = size(shape(token%dtsets(idtset)%tnons))
      dims(1:ndim) = shape(token%dtsets(idtset)%tnons)
    case (ab7_invars_upawu)
      ndim = size(shape(token%dtsets(idtset)%upawu))
      dims(1:ndim) = shape(token%dtsets(idtset)%upawu)
    case (ab7_invars_vel_cell_orig)
      ndim = size(shape(token%dtsets(idtset)%vel_cell_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%vel_cell_orig)
    case (ab7_invars_vel_orig)
      ndim = size(shape(token%dtsets(idtset)%vel_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%vel_orig)
    case (ab7_invars_wtatcon)
      ndim = size(shape(token%dtsets(idtset)%wtatcon))
      dims(1:ndim) = shape(token%dtsets(idtset)%wtatcon)
    case (ab7_invars_wtk)
      ndim = size(shape(token%dtsets(idtset)%wtk))
      dims(1:ndim) = shape(token%dtsets(idtset)%wtk)
    case (ab7_invars_xred_orig)
      ndim = size(shape(token%dtsets(idtset)%xred_orig))
      dims(1:ndim) = shape(token%dtsets(idtset)%xred_orig)
    case (ab7_invars_xredsph_extra)
      ndim = size(shape(token%dtsets(idtset)%xredsph_extra))
      dims(1:ndim) = shape(token%dtsets(idtset)%xredsph_extra)
    case (ab7_invars_ziontypat)
      ndim = size(shape(token%dtsets(idtset)%ziontypat))
      dims(1:ndim) = shape(token%dtsets(idtset)%ziontypat)
    case (ab7_invars_znucl)
      ndim = size(shape(token%dtsets(idtset)%znucl))
      dims(1:ndim) = shape(token%dtsets(idtset)%znucl)
    case (ab7_invars_bs_interp_kmult)
      ndim = size(shape(token%dtsets(idtset)%bs_interp_kmult))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_interp_kmult)
    case (ab7_invars_bs_haydock_tol)
      ndim = size(shape(token%dtsets(idtset)%bs_haydock_tol))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_haydock_tol)
    case (ab7_invars_bs_loband)
      ndim = size(shape(token%dtsets(idtset)%bs_loband))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_loband)
    case (ab7_invars_bs_eh_cutoff)
      ndim = size(shape(token%dtsets(idtset)%bs_eh_cutoff))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_eh_cutoff)
    case (ab7_invars_bs_freq_mesh)
      ndim = size(shape(token%dtsets(idtset)%bs_freq_mesh))
      dims(1:ndim) = shape(token%dtsets(idtset)%bs_freq_mesh)
    case (ab7_invars_ph_ngqpt)
      ndim = size(shape(token%dtsets(idtset)%ph_ngqpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%ph_ngqpt)
    case (ab7_invars_ph_freez_disp_ampl)
      ndim = size(shape(token%dtsets(idtset)%ph_freez_disp_ampl))
      dims(1:ndim) = shape(token%dtsets(idtset)%ph_freez_disp_ampl)
    case (ab7_invars_ph_qshift)
      ndim = size(shape(token%dtsets(idtset)%ph_qshift))
      dims(1:ndim) = shape(token%dtsets(idtset)%ph_qshift)
    case (ab7_invars_ph_qpath)
      ndim = size(shape(token%dtsets(idtset)%ph_qpath))
      dims(1:ndim) = shape(token%dtsets(idtset)%ph_qpath)
    case (ab7_invars_eph_ngqpt_fine)
      ndim = size(shape(token%dtsets(idtset)%eph_ngqpt_fine))
      dims(1:ndim) = shape(token%dtsets(idtset)%eph_ngqpt_fine)
    case (ab7_invars_ddb_ngqpt)
      ndim = size(shape(token%dtsets(idtset)%ddb_ngqpt))
      dims(1:ndim) = shape(token%dtsets(idtset)%ddb_ngqpt)
    case (ab7_invars_ddb_shiftq)
      ndim = size(shape(token%dtsets(idtset)%ddb_shiftq))
      dims(1:ndim) = shape(token%dtsets(idtset)%ddb_shiftq)
    case (ab7_invars_einterp)
      ndim = size(shape(token%dtsets(idtset)%einterp))
      dims(1:ndim) = shape(token%dtsets(idtset)%einterp)
    case (ab7_invars_kptbounds)
      ndim = size(shape(token%dtsets(idtset)%kptbounds))
      dims(1:ndim) = shape(token%dtsets(idtset)%kptbounds)
    case (ab7_invars_tmesh)
      ndim = size(shape(token%dtsets(idtset)%tmesh))
      dims(1:ndim) = shape(token%dtsets(idtset)%tmesh)

    case default
      errno = AB7_ERROR_INVARS_ATT
    end select
  end subroutine ab7_invars_get_shape
