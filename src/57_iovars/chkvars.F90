!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkvars
!! NAME
!! chkvars
!!
!! FUNCTION
!!  Examines the input string, to check whether all names are allowed.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2017 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string*(*)=string of character
!!   the string (with upper case) from the input file, to which the XYZ data is (possibly) appended
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkvars (string)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkvars'
 use interfaces_32_util
 use interfaces_57_iovars, except_this_one => chkvars
!End of the abilint section

 implicit none

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
!Note: Do not use "double quotation mark" for the string since it triggers a bug in docchk.py (abirules script)
!<ABINIT_VARS>
!A
 list_vars=                 ' accuracy acell adpimd adpimd_gamma'
 list_vars=trim(list_vars)//' algalch amu angdeg asr atvshift autoparal awtr'
!B
 list_vars=trim(list_vars)//' bandpp bdberry bdeigrf bdgw berryopt berrysav berrystep bfield bmass'
 list_vars=trim(list_vars)//' boxcenter boxcutmin brvltt builtintest'
 list_vars=trim(list_vars)//' bs_algorithm bs_calctype bs_coulomb_term bs_coupling'
 list_vars=trim(list_vars)//' bs_interp_kmult bs_interp_m3_width bs_interp_method bs_interp_mode bs_interp_prep'
 list_vars=trim(list_vars)//' bs_interp_rl_nb bs_eh_cutoff bs_exchange_term bs_freq_mesh'
 list_vars=trim(list_vars)//' bs_haydock_niter bs_haydock_tol bs_hayd_term'
 list_vars=trim(list_vars)//' bs_loband bs_nstates'
 list_vars=trim(list_vars)//' bxctmindg'
!C
 list_vars=trim(list_vars)//' cd_customnimfrqs cd_frqim_method cd_full_grid cd_imfrqs'
 list_vars=trim(list_vars)//' cd_halfway_freq cd_max_freq cd_subset_freq'
 list_vars=trim(list_vars)//' charge chempot chkexit chkprim chksymbreak chneut cineb_start cpus cpum cpuh cgtyphf'
!D
 list_vars=trim(list_vars)//' ddamp ddb_ngqpt ddb_shiftq delayperm densfor_pred densty dfield'
 list_vars=trim(list_vars)//' dfpt_sciss diecut diegap dielam dielng diemac'
 list_vars=trim(list_vars)//' diemix diemixmag diismemory dilatmx dipdip dmatpawu dmatpuopt dmatudiag'
 list_vars=trim(list_vars)//' dmft_entropy dmft_nlambda'
 list_vars=trim(list_vars)//' dmft_dc dmft_iter dmft_mxsf dmft_nwli dmft_nwlo'
 list_vars=trim(list_vars)//' dmft_read_occnd dmft_rslf dmft_solv dmft_t2g'
 list_vars=trim(list_vars)//' dmft_tolfreq dmft_tollc dmftbandi dmftbandf dmftctqmc_basis'
 list_vars=trim(list_vars)//' dmftctqmc_check dmftctqmc_correl dmftctqmc_gmove'
 list_vars=trim(list_vars)//' dmftctqmc_grnns dmftctqmc_meas dmftctqmc_mrka' 
 list_vars=trim(list_vars)//' dmftctqmc_mov dmftctqmc_order dmftctqmc_triqs_nleg'
 list_vars=trim(list_vars)//' dmftcheck dmftqmc_l dmftqmc_n dmftqmc_seed dmftqmc_therm dosdeltae dtion dynimage'
 list_vars=trim(list_vars)//' d3e_pert1_atpol d3e_pert1_dir d3e_pert1_elfd d3e_pert1_phon'
 list_vars=trim(list_vars)//' d3e_pert2_atpol d3e_pert2_dir d3e_pert2_elfd d3e_pert2_phon'
 list_vars=trim(list_vars)//' d3e_pert3_atpol d3e_pert3_dir d3e_pert3_elfd d3e_pert3_phon'
!E
 list_vars=trim(list_vars)//' ecut ecuteps ecutsigx ecutsm ecutwfn effmass efmas'
 list_vars=trim(list_vars)//' efmas_bands efmas_calc_dirs efmas_deg efmas_deg_tol'
 list_vars=trim(list_vars)//' efmas_dim efmas_dirs efmas_n_dirs efmas_ntheta'
 list_vars=trim(list_vars)//' efield einterp elph2_imagden enunit eshift esmear exchmix exchn2n3d extrapwf'
 list_vars=trim(list_vars)//' eph_intmeth eph_extrael eph_fermie eph_fsmear eph_fsewin eph_mustar '
 list_vars=trim(list_vars)//' eph_ngqpt_fine eph_task eph_transport'
!F
 list_vars=trim(list_vars)//' fband fermie_nest'
 list_vars=trim(list_vars)//' fftalg fftcache fftgw'
 list_vars=trim(list_vars)//' freqim_alpha freqremax freqremin freqspmax'
 list_vars=trim(list_vars)//' freqspmin'
 list_vars=trim(list_vars)//' friction frzfermi fxcartfactor '
 list_vars=trim(list_vars)//' f4of2_sla f6of2_sla'
!G
 list_vars=trim(list_vars)//' ga_algor ga_fitness ga_n_rules ga_opt_percent ga_rules'
 list_vars=trim(list_vars)//' genafm getbscoup getbseig getbsreso getcell'
 list_vars=trim(list_vars)//' getddb getddk getdelfd getdkdk getdkde getden getgam_eig2nkq'
 list_vars=trim(list_vars)//' gethaydock getocc getpawden getqps getscr'
 list_vars=trim(list_vars)//' getwfkfine'
 list_vars=trim(list_vars)//' getsuscep '
 list_vars=trim(list_vars)//' getvel getwfk getwfq getxcart getxred get1den get1wf goprecon goprecprm'
 list_vars=trim(list_vars)//' gpu_devices gpu_linalg_limit gwcalctyp gwcomp gwencomp gwgamma gwmem'
 list_vars=trim(list_vars)//' gwpara gwrpacorr gw_customnfreqsp'
 list_vars=trim(list_vars)//' gw_frqim_inzgrid gw_frqre_inzgrid gw_frqre_tangrid gw_freqsp'
 list_vars=trim(list_vars)//' gw_invalid_freq '
 list_vars=trim(list_vars)//' gw_qprange gw_nqlwl gw_nstep gw_qlwl'
 list_vars=trim(list_vars)//' gw_sctype gw_sigxcore gw_toldfeig'
 list_vars=trim(list_vars)//' gwfockmix'
 list_vars=trim(list_vars)//' gwls_sternheimer_kmax gwls_kmax_complement gwls_kmax_poles'
 list_vars=trim(list_vars)//' gwls_kmax_analytic gwls_kmax_numeric'
 list_vars=trim(list_vars)//' gwls_list_proj_freq gwls_nseeds gwls_n_proj_freq gwls_recycle'
 list_vars=trim(list_vars)//' gwls_first_seed gwls_model_parameter gwls_npt_gauss_quad'
 list_vars=trim(list_vars)//' gwls_diel_model gwls_print_debug gwls_band_index gwls_exchange gwls_correlation'
!I
 list_vars=trim(list_vars)//' iatcon iatfix iatfixx iatfixy iatfixz iatsph'
 list_vars=trim(list_vars)//' iboxcut icoulomb icutcoul ieig2rf'
 list_vars=trim(list_vars)//' imgmov inclvkb intxc iomode ionmov iqpt'
 list_vars=trim(list_vars)//' iprcel iprcfc irandom irdbscoup'
 list_vars=trim(list_vars)//' irdbseig irdbsreso irdddb irdddk irdden'
 list_vars=trim(list_vars)//' irdhaydock irdpawden irdqps'
 list_vars=trim(list_vars)//' irdscr irdsuscep irdwfk irdwfq ird1den'
 list_vars=trim(list_vars)//' irdwfkfine'
 list_vars=trim(list_vars)//' ird1wf iscf isecur istatimg istatr istatshft istwfk ixc ixcpositron'
 list_vars=trim(list_vars)//' irdvdw'
!J
 list_vars=trim(list_vars)//' jdtset jellslab jfielddir jpawu'
!K
 list_vars=trim(list_vars)//' kberry kpt kptbounds kptgw'
 list_vars=trim(list_vars)//' kptnrm kptopt kptrlatt kptrlen kssform'
!L
 list_vars=trim(list_vars)//' lexexch localrdwf lpawu'
 list_vars=trim(list_vars)//' lotf_classic lotf_nitex lotf_nneigx lotf_version'
!M
 list_vars=trim(list_vars)//' max_ncpus macro_uj maxestep maxnsym mdf_epsinf mdtemp mdwall'
 list_vars=trim(list_vars)//' magconon magcon_lambda mbpt_sciss'
 list_vars=trim(list_vars)//' mep_mxstep mep_solver mem_test mixalch'
 list_vars=trim(list_vars)//' mqgrid mqgriddg'
!N
 list_vars=trim(list_vars)//' natcon natfix natfixx natfixy natfixz'
 list_vars=trim(list_vars)//' natom natrd natsph natsph_extra natvshift nband nbandkss nbandhf'
 list_vars=trim(list_vars)//' nbdblock nbdbuf nberry nconeq nc_xccc_gspace'
 list_vars=trim(list_vars)//' nctime ndivk ndivsm ndtset neb_algo neb_spring'
 list_vars=trim(list_vars)//' nfreqim nfreqre nfreqsp ngfft ngfftdg'
 list_vars=trim(list_vars)//' ngkpt ngqpt nimage nkpath nkpt nkptgw nkpthf'
 list_vars=trim(list_vars)//' nline nloc_alg nloc_mem nnos nnsclo nnsclohf'
 list_vars=trim(list_vars)//' nobj nomegasf nomegasi nomegasrd noseinert npband'
 list_vars=trim(list_vars)//' npfft nphf npimage npkpt nppert npsp npspinor'
 list_vars=trim(list_vars)//' npulayit npvel npwkss'
 list_vars=trim(list_vars)//' np_slk nqpt nqptdm nscforder nshiftk nshiftq'
 list_vars=trim(list_vars)//' nspden nspinor nsppol nstep nsym'
 list_vars=trim(list_vars)//' ntime ntimimage ntypalch ntypat nucdipmom nwfshist nzchempot'
!O
 list_vars=trim(list_vars)//' objaat objbat objaax objbax objan objbn objarf'
 list_vars=trim(list_vars)//' objbrf objaro objbro objatr objbtr occ'
 list_vars=trim(list_vars)//' occopt omegasimax omegasrdmax optcell optdriver optforces'
 list_vars=trim(list_vars)//' optnlxccc optstress ortalg'
!P
 list_vars=trim(list_vars)//' paral_atom paral_kgb paral_rf pawcpxocc pawcross'
 list_vars=trim(list_vars)//' pawecutdg pawfatbnd pawlcutd pawlmix'
 list_vars=trim(list_vars)//' pawmixdg pawnhatxc pawnphi pawntheta pawnzlm pawoptmix pawoptosc pawovlp'
 list_vars=trim(list_vars)//' pawprtdos pawprtvol pawprtwf pawprt_b pawprt_k pawspnorb pawstgylm'
 list_vars=trim(list_vars)//' pawsushat pawujat pawujrad pawujv pawusecp pawxcdev pimass pimd_constraint'
 list_vars=trim(list_vars)//' pitransform ph_ndivsm ph_nqpath ph_qpath ph_ngqpt'
 list_vars=trim(list_vars)//' ph_wstep ph_intmeth ph_smear ph_nqshift ph_qshift'
 list_vars=trim(list_vars)//' plowan_bandi plowan_bandf plowan_compute plowan_iatom plowan_it plowan_lcalc'
 list_vars=trim(list_vars)//' plowan_natom plowan_nbl plowan_nt plowan_projcalc plowan_realspace'
 list_vars=trim(list_vars)//' polcen posdoppler positron posnstep posocc postoldfe postoldff ppmfrq ppmodel'
 list_vars=trim(list_vars)//' prepanl prepgkk papiopt'
 list_vars=trim(list_vars)//' prtatlist prtbbb prtbltztrp prtcif prtden'
 list_vars=trim(list_vars)//' prtdensph prtdipole prtdos prtdosm prtebands prtefg prteig prtelf'
 list_vars=trim(list_vars)//' prtfc prtfsurf prtgden prtgeo prtgsr prtgkk prtkden prtkpt prtlden'
 list_vars=trim(list_vars)//' prtnabla prtnest prtphbands prtphdos prtphsurf prtposcar prtpot prtpsps'
 list_vars=trim(list_vars)//' prtspcur prtstm prtsuscep prtvclmb prtvha prtvdw prtvhxc'
 list_vars=trim(list_vars)//' prtvol prtvpsp prtvxc prtwant prtwf prtwf_full prtxml prt1dm ptcharge'
 list_vars=trim(list_vars)//' pvelmax pw_unbal_thresh'
!Q
 list_vars=trim(list_vars)//' qmass qprtrb qpt qptdm qptnrm '
 list_vars=trim(list_vars)//' qptopt qptrlatt quadmom'
!R
 list_vars=trim(list_vars)//' random_atpos ratsph ratsph_extra rcut'
 list_vars=trim(list_vars)//' recefermi recgratio recnpath recnrec recptrott recrcut rectesteg rectolden'
 list_vars=trim(list_vars)//' red_dfield red_efield red_efieldbar restartxf rfasr'
 list_vars=trim(list_vars)//' rfatpol rfddk rfdir rfelfd rfmagn rfmeth rfphon'
 list_vars=trim(list_vars)//' rfstrs rfuser rf2_dkdk rf2_dkde rf2_pert1_dir rf2_pert2_dir rhoqpmix rprim'
 !These input parameters are obsolete (keep them for compatibility)
 list_vars=trim(list_vars)//' rf1atpol rf1dir rf1elfd rf1phon'
 list_vars=trim(list_vars)//' rf2atpol rf2dir rf2elfd rf2phon'
 list_vars=trim(list_vars)//' rf3atpol rf3dir rf3elfd rf3phon'
!S
 list_vars=trim(list_vars)//' scalecart shiftk shiftq signperm'
 list_vars=trim(list_vars)//' slabwsrad slabzbeg slabzend smdelta so_psp'
 list_vars=trim(list_vars)//' spbroad spgaxor spgorig spgroup spgroupma spinat spinmagntarget spmeth '
 list_vars=trim(list_vars)//' spnorbscl stmbias strfact string_algo strprecon strtarget supercell'
 list_vars=trim(list_vars)//' symafm symchi symdynmat symmorphi symrel symsigma'
!T
 list_vars=trim(list_vars)//' td_maxene td_mexcit tfkinfunc tfw_toldfe timopt tl_nprccg tl_radius'
 list_vars=trim(list_vars)//' tnons toldfe tolmxde toldff tolimg tolmxf tolrde tolrff tolsym'
 list_vars=trim(list_vars)//' tolvrs tolwfr tphysel tsmear typat'
!U
 list_vars=trim(list_vars)//' ucrpa ucrpa_bands ucrpa_window udtset upawu usedmatpu'
 list_vars=trim(list_vars)//' usedmft useexexch usekden use_nonscf_gkk usepawu usepotzero'
 list_vars=trim(list_vars)//' useria userib useric userid userie'
 list_vars=trim(list_vars)//' userra userrb userrc userrd userre'
 list_vars=trim(list_vars)//' usewvl usexcnhat useylm use_gemm_nonlop use_gpu_cuda use_slk'
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
 list_vars=trim(list_vars)//' wfoptalg wtatcon wtk wtq'
 list_vars=trim(list_vars)//' wvl_bigdft_comp wvl_crmult wvl_frmult wvl_hgrid wvl_ngauss wvl_nprccg'
 list_vars=trim(list_vars)//' w90iniprj w90prtunk'
!X
 list_vars=trim(list_vars)//' xangst xcart xc_denpos xc_tb09_c xred xredsph_extra xyzfile'
!Y
!Z
 list_vars=trim(list_vars)//' zcut zeemanfield znucl'

!Logical input variables
 list_logicals=' SpinPolarized'

!String input variables
 list_strings=' XCname wfk_task'
!</ABINIT_VARS>

!Extra token, also admitted :
!<ABINIT_UNITS>
 list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha'
 list_vars=trim(list_vars)//' Hartree Hartrees K Ry Rydberg Rydbergs T Tesla'
!</ABINIT_UNITS>

!<ABINIT_OPERATORS>
 list_vars=trim(list_vars)//' sqrt end'
!</ABINIT_OPERATORS>

!Transform to upper case
 call inupper(list_vars)
 call inupper(list_logicals)
 call inupper(list_strings)

 call chkvars_in_string(protocol1, list_vars, list_logicals, list_strings, string)

end subroutine chkvars
!!***

!!****f* ABINIT/chkvars_in_string
!! NAME
!!  chkvars_in_string
!!
!! FUNCTION
!!  Analyze variable names in string. Abort if name is not recognized.
!!
!! INPUTS
!!  protocol=
!!    0 if parser does not accept multiple datasets and +* syntax (e.g. anaddb)
!!    1 if parser accepts multiple datasets and +* syntax (e.g. abinit)
!!
!!  list_vars(len=*)=string with the (upper case) names of the variables (excluding logicals and chars).
!!  list_logicals(len=*)=string with the (upper case) names of the logical variables.
!!  list_strings(len=*)=string with the (upper case) names of the character variables.
!!  string(len=*)=string (with upper case) from the input file.
!!
!! OUTPUT
!!  Abort if variable name is not recognized.
!!
!! PARENTS
!!      chkvars,m_anaddb_dataset
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkvars_in_string(protocol, list_vars, list_logicals, list_strings, string)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkvars_in_string'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: protocol
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: list_logicals,list_strings,list_vars

!Local variables-------------------------------
 character,parameter :: blank=' '
!scalars
 integer :: index_blank,index_current,index_endword,index_endwordnow,index_list_vars
 character(len=500) :: message

!************************************************************************

 index_current=1
 do ! Infinite do-loop, to identify the presence of each potential variable names

   if(len_trim(string)<=index_current)exit
   index_blank=index(string(index_current:),blank)+index_current-1

   if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_current:index_current))/=0)then

     index_endword = index_blank -1
     if (protocol == 1) then
       ! Skip characters like : + or the digits at the end of the word
       ! Start from the blank that follows the end of the word
       do index_endword=index_blank-1,index_current,-1
         if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
       end do
     end if
     !write(std_out,*)"Will analyze:", string(index_current:index_endword)

     ! Find the index of the potential variable name in the list of variables
     index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)

     ! Treat the complications due to the possibility of images
     if (index_list_vars==0 .and. protocol==1) then

       ! Treat possible LASTIMG appendix
       if(index_endword-6>=1)then
         if(string(index_endword-6:index_endword)=='LASTIMG')index_endword=index_endword-7
       end if

       ! Treat possible IMG appendix
       if(index_endword-2>=1)then
         if(string(index_endword-2:index_endword)=='IMG')index_endword=index_endword-3
       end if

       index_endwordnow=index_endword

       ! Again skip characters like : + or the digits before IMG
       ! Start from the blank that follows the end of the word
       do index_endword=index_endwordnow,index_current,-1
         if(index('ABCDEFGHIJKLMNOPQRSTUVWXYZ',string(index_endword:index_endword))/=0)exit
       end do

       ! Find the index of the potential variable name in the list of variables
       index_list_vars=index(list_vars,blank//string(index_current:index_endword)//blank)
     end if

     if(index_list_vars==0)then

       ! Treat possible logical input variables
       if(index(list_logicals,blank//string(index_current:index_endword)//blank)/=0)then
         !write(std_out,*)"Found logical variable: ",string(index_current:index_endword)
         index_blank=index(string(index_current:),blank)+index_current-1
         if(index(' F T ',string(index_blank:index_blank+2))==0)then
           write(message, '(8a)' )&
&           'Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&           'This variable should be given a logical value (T or F), but the following string was found :',&
&           string(index_blank:index_blank+2),ch10,&
&           'Action: check your input file. You likely misused the input variable.'
           MSG_ERROR(message)
         else
           index_blank=index_blank+2
         end if
!        Treat possible string input variables
       else if(index(list_strings,blank//string(index_current:index_endword)//blank)/=0)then
!        Every following string is accepted
         !write(std_out,*)"Found string variable: ",string(index_current:index_endword)
         !write(std_out,*)"in string: ",trim(string(index_current:))
         index_current=index(string(index_current:),blank)+index_current
         index_blank=index(string(index_current:),blank)+index_current-1
         !write(std_out,*)"next:: ",string(index_current:index_endword)

!        If still not admitted, then there is a problem
       else
         write(message, '(7a)' )&
&         'Found the token ',string(index_current:index_endword),' in the input file.',ch10,&
&         'This name is not one of the registered input variable names (see http://www.abinit.org/doc).',ch10,&
&         'Action: check your input file. You likely mistyped the input variable.'
         MSG_ERROR(message)
       end if
     end if
   end if

   index_current=index_blank+1
 end do

end subroutine chkvars_in_string
!!***
