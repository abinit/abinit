!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout
!! NAME
!!  wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!   The routine uses optional arguments, therefore the interface must be explicit.
!!   Be careful when writing CPP macros that use wrtout since abilint won't see the call
!!   and no interface will be added to the source file.
!!
!! PARENTS
!!      abinit,afterscfloop,anaddb,append_xyz,atm2fft,atomden,berryphase
!!      berryphase_new,bethe_salpeter,bonds_lgth_angles,bsepostproc,calc_efg
!!      calc_fc,calc_optical_mels,calc_rpa_functional,calc_sigc_me,calc_sigx_me
!!      calc_ubare,calc_ucrpa,calc_vhxc_me,calcdensph,cchi0,cchi0q0
!!      cchi0q0_intraband,cgwf,chebfi,chern_number,chiscwrt,chkdpr,chkinp
!!      chkint_prt,chkpawovlp,clnup1,clnup2,cohsex_me,compute_anharmonics
!!      compute_kgb_indicator,compute_levels,constrf,cut3d,d2frnl,datafordmft
!!      ddb_diel,ddb_elast,ddb_internalstr,ddb_interpolate,ddb_io_out,ddb_piezo
!!      debug_tools,deloc2xcart,denfgr,dfpt_cgwf,dfpt_looppert,dfpt_mkrho
!!      dfpt_nselt,dfpt_nstdy,dfpt_nstpaw,dfpt_prtene,dfpt_prtph,dfpt_scfcv
!!      dfpt_symph,dfpt_vtowfk,dfpt_wfkfermi,dfptff_initberry,dfptnl_doutput
!!      dfptnl_loop,dfptnl_mv,dielmt,dieltcel,dmft_solve,dos_hdr_write,driver
!!      dyson,echo_xc_name,eig2stern,eliashberg_1d,elph2_fanddw,elphon
!!      elpolariz,elt_ewald,entropyrec,ep_fs_weights,ep_setupqpt,eph,erlxconv
!!      evdw_wannier,exc_build_block,exc_build_ham,exc_den,exc_plot,extraprho
!!      fconv,fermi_green,fermisolverec,fftprof,find_getdtset,finddistrproc
!!      findmin,findminscf,first_rec,fred2fdeloc,fresidrsp,fsumrule,get_all_gkq
!!      get_fs_bands,get_npert_rbz,get_nv_fs_en,get_nv_fs_temp,get_tau_k
!!      getcgqphase,getcut,getdim_nloc,getghc,getmpw,getnel,getng,getshell
!!      gran_potrec,green_kernel,gstate,gstateimg,gwls_polarisability
!!      harmonic_thermo,hdr_vs_dtset,herald,hermit,hubbard_one,importxyz
!!      impurity_solve,inarray,ingeo,ingeobld,initberry,initmpi_grid,initorbmag
!!      initro,initwf,inkpts,inpgkk,inpspheads,instrng,intagm,integrate_gamma
!!      integrate_gamma_alt,integrate_gamma_tr,integrate_gamma_tr_lova,invars1
!!      invars2,inwffil,ioniondist,ioprof,irrzg,kpgio,kramerskronig,ks_ddiago
!!      lapackprof,ldau_self,leave_new,linemin,littlegroup_pert,littlegroup_q
!!      lobpcgwf,local_ks_green,m_abi_etsf,m_abihist,m_abilasi,m_anaddb_dataset
!!      m_argparse,m_atom,m_bfgs,m_bs_defs,m_bse_io,m_bz_mesh,m_cgtools,m_chi0
!!      m_crystal,m_ddb,m_ddb_hdr,m_ddk,m_double_grid,m_dtset,m_dvdb,m_dynmat
!!      m_dyson_solver,m_ebands,m_effective_potential
!!      m_effective_potential_file,m_energy,m_entropyDMFT,m_epjdos,m_errors
!!      m_esymm,m_eval_lotf,m_exc_diago,m_exc_itdiago,m_exc_spectra,m_exit
!!      m_fft,m_fft_mesh,m_fft_prof,m_fftcore,m_fftw3,m_fit_polynomial_coeff
!!      m_fock,m_fstab,m_geometry,m_gkk,m_gpu_detect,m_green,m_gruneisen
!!      m_gsphere,m_hamiltonian,m_harmonics_terms,m_haydock,m_hdr,m_hexc
!!      m_hidecudarec,m_hu,m_ifc,m_initcuda,m_invovl,m_io_kss,m_io_screening
!!      m_ioarr,m_iowf,m_iterators,m_kxc,m_libpaw_libxc,m_libxc_functionals
!!      m_lotf,m_matlu,m_matrix,m_melemts,m_mep,m_multibinit_dataset,m_nctk
!!      m_numeric_tools,m_oper,m_optic_tools,m_paw_an,m_paw_dmft,m_paw_gaussfit
!!      m_paw_ij,m_paw_io,m_paw_slater,m_paw_sphharm,m_pawang,m_pawdij,m_pawfgr
!!      m_pawfgrtab,m_pawpsp,m_pawpwij,m_pawrad,m_pawrhoij,m_pawtab,m_phgamma
!!      m_phonons,m_phpi,m_pimd,m_plowannier,m_polynomial_coeff,m_ppmodel
!!      m_pptools,m_pred_lotf,m_pretty_rec,m_psps,m_psxml2ab,m_ptgroups
!!      m_qparticles,m_rec,m_rf2,m_screen,m_screening,m_self,m_shirley,m_sigma
!!      m_sigmaph,m_skw,m_slk,m_strain,m_vcoul,m_wfd,m_wfk,m_work_var_lotf
!!      m_xc_vdw,m_xpapi,mag_constr_e,mblktyp1,mblktyp5,memana,memorf,memory
!!      metric,mka2f,mka2fQgrid,mka2f_tr,mka2f_tr_lova,mkcore_paw,mkcore_wvl
!!      mkfilename,mkfskgrid,mklocl_recipspace,mklocl_wavelets,mknormpath
!!      mkph_linwid,mkqptequiv,mkrho,mlwfovlp,mlwfovlp_proj,mlwfovlp_projpaw
!!      mlwfovlp_pw,mlwfovlp_qp,mlwfovlp_seedname,mlwfovlp_setup,mover
!!      mover_effpot,mpi_setup,mrgddb,mrgdv,mrggkk,mrgscr,multibinit
!!      multipoles_out,my_calc_wfwfg,newfermie1,newkpt,newocc,newton
!!      nlenergyrec,nonlinear,normsq_gkq,optic,out1dm,out_acknowl,outelph
!!      outgkk,outkss,outqmc,outscfcv,outvars,outwant,partial_dos_fractions
!!      paw_mknewh0,paw_qpscgw,pawdenpot,pawdensities,pawmkaewf,pawmkrhoij
!!      pawprt,pawpuxinit,pawuenergy,pawuj_det,pawuj_red,pawuj_utils,pawxenergy
!!      pimd_nosehoover_nvt,polcart,posdoppler,poslifetime,pred_delocint
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_nose,pred_verlet
!!      predictimg,prep_calc_ucrpa,prtefield,prteigrs,prtene,prtimg,prtrhomxmn
!!      prtspgroup,prtxf,prtxfase,prtxvf,psichi_renormalization,psolver_hartree
!!      psolver_kernel,psolver_rhohxc,psp10in,psp1in,psp2in,psp2lo,psp3in
!!      psp5in,psp6in,psp8in,psp9in,pspatm,pspini,pspnl_hgh_rec
!!      pspnl_operat_rec,qmc_prep_ctqmc,randac,random_stopping_power
!!      rayleigh_ritz,read_gkk,read_plowannier,recursion_nl,remove_inversion
!!      respfn,rf2_init,rotmat,scfcge,scfcv,scfeig,scfopt,scprqt,screening
!!      setnoccmmp,setrhoijpbe0,setsymrhoij,setup1,setup_bse,setup_bse_interp
!!      setup_positron,setup_screening,setup_sigma,shellstruct,sigma,smpbz
!!      spectral_function,stress,sumrule,symanal,symatm,symaxes,symcharac
!!      symfind,symkchk,symkpt,symlatt,symplanes,symspgr,tddft,testkgrid,thmeig
!!      timana,uderiv,ujdet,update_e_field_vars,vdw_dftd2,vdw_dftd3
!!      vdw_kernelgen,vtorho,vtorhorec,vtorhotf,vtowfk,wfconv,wfd_mkrho
!!      wfd_pawrhoij,wfk_analyze,wfsinp,wrt_moldyn_netcdf,wrtloctens
!!      wvl_denspot_set,wvl_descr_atoms_set_sym,wvl_descr_psp_set,wvl_hpsitopsi
!!      wvl_initro,wvl_memory,wvl_mkrho,wvl_nl_gradient,wvl_projectors_set
!!      wvl_psitohpsi,wvl_rwwf,wvl_setboxgeometry,wvl_setngfft
!!      wvl_tail_corrections,wvl_wfs_set,wvl_wfsinp_disk,wvl_wfsinp_reformat
!!      wvl_wfsinp_scratch
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wrtout(unit,msg,mode_paral,do_flush)

 use defs_basis

 use m_xmpi,      only : xmpi_world, xmpi_comm_rank, xmpi_comm_size
 use m_io_tools,  only : flush_unit, write_lines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout'
 use interfaces_14_hidewrite, except_this_one => wrtout
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: do_flush

!Local variables-------------------------------
 integer :: comm,me,nproc
 integer,save :: master=0
 logical :: my_flush
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if ((unit == std_out).and.(.not.do_write_log)) RETURN
 if (unit == dev_null) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral
 my_flush = .false.; if (PRESENT(do_flush)) my_flush = do_flush

!Communicator is xmpi_world by default, except for the parallelization over images
 if (abinit_comm_output/=-1) then
   comm=abinit_comm_output
 else
   comm=xmpi_world
 end if

!Determine who I am in COMM_WORLD
 nproc = xmpi_comm_size(comm)
 me    = xmpi_comm_rank(comm)

 if( (my_mode_paral=='COLL') .or. (nproc==1) ) then
   if (me==master) then
     call wrtout_myproc(unit, msg, do_flush=my_flush)
   end if

 else if (my_mode_paral=='PERS') then
   call write_lines(unit,msg)

   ! Flush unit
   if (my_flush) then
     call flush_unit(unit)
   end if

 else if (my_mode_paral=='INIT') then
   master=unit

 else
   write(string,'(7a)')ch10,&
&   'wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)
 end if

end subroutine wrtout
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/wrtout_myproc
!! NAME
!!  wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc. For parallel or sequential output use wrtout()
!!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!!  Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! INPUTS
!!  unit=unit number for writing
!!  message=(character(len=*)) message to be written
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      wrtout
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wrtout_myproc(unit,message,do_flush) ! optional argument

 use defs_basis
 use m_profiling_abi

 use m_xmpi,       only : xmpi_sum
 use m_specialmsg, only : specialmsg_setcount
 use m_io_tools,   only : flush_unit, write_lines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtout_myproc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: message
 logical,optional,intent(in) :: do_flush

!Local variables-------------------------------
!scalars
 integer :: i_one=1
 logical :: print_std_err

!******************************************************************

 print_std_err = (unit == std_out .and. std_out /= std_err .and. &
& (index(trim(message),'BUG')/=0.or.index(trim(message),'ERROR')/=0))

!Print message
 call write_lines(unit,message)
 if (print_std_err) call write_lines(std_err,message)

!Append "Contact Abinit group" to BUG messages
 if( index(trim(message),'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action: contact ABINIT group (please attach the output of `abinit -b`)'
   write(unit,*)
   if (print_std_err) then
     write(std_err, '(a)' ) '  Action: contact ABINIT group (please attach the output of `abinit -b`)'
     write(std_err,*)
   end if
 end if

!Count the number of warnings and comments. Only take into
!account unit std_out, in order not to duplicate these numbers.
 if( index(trim(message),'WARNING') /= 0 .and. unit==std_out )then
   call specialmsg_setcount(n_add_warning=i_one)
 end if
 if( index(trim(message),'COMMENT') /= 0 .and. unit==std_out )then
   call specialmsg_setcount(n_add_comment=i_one)
 end if
 if( index(trim(message),'Exit') /= 0 )then
   call specialmsg_setcount(n_add_exit=i_one)
 end if

!Flush unit
 if (present(do_flush)) then
   if (do_flush) call flush_unit(unit)
 end if
#ifdef DEBUG_MODE
 call flush_unit(unit)
 if (print_std_err) call flush_unit(std_err)
#endif

end subroutine wrtout_myproc
!!***
