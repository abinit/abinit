!!****m* ABINIT/m_mover
!! NAME
!!  m_mover
!!
!! FUNCTION
!! Move ion or change acell according to forces and stresses
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, SE)
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

module m_mover

 use defs_basis
 use m_abicore
 use m_errors
 use m_profiling_abi
 use m_abimover
 use m_abihist
 use m_dtset
 use m_xmpi
 use m_nctk
 use m_dtfil
 use m_yaml
#ifdef HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_LOTF
 use lotfpath
 use m_pred_lotf
#endif

 use defs_abitypes,        only : MPI_type
 use m_fstrings,           only : strcat, sjoin, indent
 use m_symtk,              only : matr3inv, symmetrize_xred
 use m_geometry,           only : fcart2fred, chkdilatmx, xred2xcart
 use m_time,               only : abi_wtime, sec2str
 use m_exit,               only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_electronpositron,   only : electronpositron_type
 use m_scfcv,              only : scfcv_t, scfcv_run
 use m_effective_potential,only : effective_potential_type, effective_potential_evaluate
 use m_dtfil,              only : dtfil_init_time
 use m_initylmg,           only : initylmg
 use m_xfpack,             only : xfh_update
 use m_precpred_1geo,      only : precpred_1geo
 use m_pred_delocint,      only : pred_delocint
 use m_pred_bfgs,          only : pred_bfgs, pred_lbfgs
 use m_pred_fire,          only : pred_fire
 use m_pred_isokinetic,    only : pred_isokinetic
 use m_pred_diisrelax,     only : pred_diisrelax
 use m_pred_nose,          only : pred_nose
 use m_pred_srkhna14,      only : pred_srkna14
 use m_pred_isothermal,    only : pred_isothermal
 use m_pred_verlet,        only : pred_verlet
 use m_pred_velverlet,     only : pred_velverlet
 use m_pred_moldyn,        only : pred_moldyn
 use m_pred_langevin,      only : pred_langevin
 use m_pred_steepdesc,     only : pred_steepdesc
 use m_pred_simple,        only : pred_simple, prec_simple
 use m_pred_hmc,           only : pred_hmc
 use m_generate_training_set, only : generate_training_set
 use m_wvl_wfsinp, only : wvl_wfsinp_reformat
 use m_wvl_rho,      only : wvl_mkrho
 use m_effective_potential_file, only : effective_potential_file_mapHistToRef
#if defined DEV_MS_SCALEUP
 use scup_global, only : global_set_parent_iter,global_set_print_parameters
#endif
 use m_scup_dataset
 implicit none

 private
!!***

 public :: mover
!!***

contains
!!***

!!****f* ABINIT/mover
!! NAME
!! mover
!!
!! FUNCTION
!! Move ion or change acell acording to forces and stresses
!!
!! INPUTS
!!  amu_curr(ntypat)=mass of each atom for the current image
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points treated by this node
!!   |  angular momentum for nonlocal pseudopotential
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in unit cell
!!   |  except on first call (hartree/bohr); updated on output
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=informations about MPI parallelization
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  npwarr(nkpt)=number of planewaves in basis and boundary at this k point.
!!  nattyp(ntypat)= # atoms of each type.
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  scup_dtset <type(scup_dtset_type) = derived datatype holding all options
 !!            for the evaluation of an effective electronic model using SCALE UP
!!
!! OUTPUT
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!
!! SIDE EFFECTS
!! Rest of i/o is related to lda
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mcg)=array for planewave coefficients of wavefunctions.
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  initialized= if 0 the initialisation of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  occ(mband*nkpt*nsppol=occupation number for each band (usually 2) at each k point.
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  vel(3,natom)=old value of velocity; updated on output
!!  vel_cell(3,3)=old value of cell parameters velocity; updated on output
!!  xred(3,natom)=reduced dimensionless atomic coordinates; updated on output
!!  xred_old(3,natom)=work space for old xred
!!  eff_pot<type(effective_potential_type)> = optional,effective_potential datatype
!!  verbose = optional, default is true, flag to disable the verbose mode
!!  write_HIST = optional, default is true, flag to disble the write of the HIST file
!!
!! NOTES
!! This subroutine uses the arguments natom, xred, vel, amu_curr,
!! vis, and dtion (the last two contained in dtset) to make
!! molecular dynamics updates.  The rest of the lengthy
!! argument list supports the underlying lda computation
!! of forces, returned from subroutine scfcv
!!
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
!!
!! PARENTS
!!      gstate,mover_effpot
!!
!! CHILDREN
!!      abiforstr_fin,abiforstr_ini,abihist_bcast,abihist_compare_and_copy
!!      abihist_free,abihist_init,abimover_fin,abimover_ini
!!      chkdilatmx,crystal_init,dtfil_init_time
!!      effective_potential_evaluate,erlxconv,fcart2fred,fconv,hist2var
!!      initylmg,matr3inv,mttk_fin,mttk_ini,prec_simple,pred_bfgs,pred_delocint
!!      pred_diisrelax,pred_hmc,pred_isokinetic,pred_isothermal,pred_langevin
!!      pred_lbfgs,pred_lotf,pred_moldyn,pred_nose,pred_simple,pred_srkna14
!!      pred_steepdesc,pred_velverlet,pred_verlet,prtxfase,read_md_hist
!!      scfcv_run,status,symmetrize_xred,var2hist,vel2hist,write_md_hist
!!      wrt_moldyn_netcdf,wrtout,wvl_mkrho,wvl_wfsinp_reformat,xfh_update
!!      xmpi_barrier,xmpi_isum,xmpi_wait
!!
!! SOURCE

subroutine mover(scfcv_args,ab_xfh,acell,amu_curr,dtfil,&
& electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old,&
& effective_potential,filename_ddb,verbose,writeHIST,scup_dtset)

!Arguments ------------------------------------
!scalars
type(scfcv_t),intent(inout) :: scfcv_args
type(datafiles_type),intent(inout),target :: dtfil
type(electronpositron_type),pointer :: electronpositron
type(ab_xfh_type),intent(inout) :: ab_xfh
type(effective_potential_type),optional,intent(inout) :: effective_potential
logical,optional,intent(in) :: verbose
logical,optional,intent(in) :: writeHIST
character(len=fnlen),optional,intent(in) :: filename_ddb
!arrays
real(dp),intent(inout) :: acell(3)
real(dp), intent(in),target :: amu_curr(:) !(scfcv%dtset%ntypat)
real(dp), pointer :: rhog(:,:),rhor(:,:)
real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom),xred_old(3,scfcv_args%dtset%natom)
real(dp), intent(inout) :: vel(3,scfcv_args%dtset%natom),vel_cell(3,3),rprimd(3,3)
type(scup_dtset_type),optional, intent(inout) :: scup_dtset

!Local variables-------------------------------
!scalars
integer,parameter :: level=102,master=0
type(abihist) :: hist,hist_prev
type(abimover) :: ab_mover
type(abimover_specs) :: specs
type(abiforstr) :: preconforstr ! Preconditioned forces and stress
type(delocint) :: deloc
type(mttk_type) :: mttk_vars
integer :: itime,icycle,itime_hist,iexit=0,ifirst,ihist_prev,ihist_prev2,timelimit_exit,ncycle,nhisttot,kk,jj,me
integer :: ntime,option,comm
integer :: nerr_dilatmx,my_quit,ierr,quitsum_request
integer ABI_ASYNC :: quitsum_async
character(len=500) :: message
!character(len=500) :: dilatmx_errmsg
character(len=8) :: stat4xml
character(len=35) :: fmt
character(len=fnlen) :: filename,fname_ddb, name_file
character(len=500) :: MY_NAME = "mover"
real(dp) :: favg
logical :: DEBUG=.FALSE., need_verbose=.TRUE.,need_writeHIST=.TRUE.
logical :: need_scfcv_cycle = .TRUE., need_elec_eval = .FALSE.
logical :: change,useprtxfase
logical :: skipcycle
integer :: minIndex,ii,similar,conv_retcode
integer :: iapp
real(dp) :: minE,wtime_step,now,prev
logical :: file_exists
!arrays
real(dp) :: gprimd(3,3),rprim(3,3),rprimd_prev(3,3)
real(dp),allocatable :: fred_corrected(:,:),xred_prev(:,:)
! ***************************************************************
 need_verbose=.TRUE.
 if(present(verbose)) need_verbose = verbose

 need_writeHIST=.TRUE.
 if(present(writeHIST)) need_writeHIST = writeHIST

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(MY_NAME) == MY_NAME) then
   if (need_verbose) then
     write(std_out,*)"Enabling timelimit check in function: ",trim(MY_NAME)," with timelimit: ",trim(sec2str(get_timelimit()))
   end if
 end if

!Table of contents
!(=>) Refers to an important call (scfcv,pred_*)
!
!01. Initialization of indexes and allocations of arrays
!02. Particularities of each predictor
!03. Set the number of iterations ntime
!04. Try to read history of previous calculations
!05. Allocate the hist structure
!06. First output before any itime or icycle
!07. Fill the history of the first SCFCV
!08. Loop for itime (From 1 to ntime)
!09. Loop for icycle (From 1 to ncycle)
!10. Output for each icycle (and itime)
!11. Symmetrize atomic coordinates over space group elements
!12. => Call to SCFCV routine and fill history with forces
!13. Write the history into the _HIST file
!14. Output after SCFCV
!15. => Test Convergence of forces and stresses
!16. => Precondition forces, stress and energy
!17. => Call to each predictor
!18. Use the history  to extract the new values
!19. End loop icycle
!20. End loop itime
!21. Set the final values of xred
!22. XML Output at the end
!23. Deallocate hist and ab_mover datatypes
!
 call abimover_ini(ab_mover,amu_curr,dtfil,scfcv_args%dtset,specs)

 if(ab_mover%ionmov==10 .or. ab_mover%ionmov==11)then
   call delocint_ini(deloc)
 end if

 if (ab_mover%ionmov==13)then
   call mttk_ini(mttk_vars,ab_mover%nnos)
 end if

!###########################################################
!### 03. Set the number of iterations ntime
!###     By default ntime==1 but if the user enter a lower
!###     value mover will execute at least one iteration

 if (scfcv_args%dtset%ntime<1)then
   ntime=1
 else
   ntime=scfcv_args%dtset%ntime
 end if

!###########################################################
!### 04. Try to read history of previous calculations
!###     It requires access to the NetCDF library

!Init MPI data
 comm=scfcv_args%mpi_enreg%comm_cell
 me=xmpi_comm_rank(comm)


#if defined HAVE_NETCDF
 filename=trim(ab_mover%filnam_ds(4))//'_HIST.nc'

 if (ab_mover%restartxf<0)then
!  Read history from file (and broadcast if MPI)
   if (me==master) then
     call read_md_hist(filename,hist_prev,specs%isVused,specs%isARused,ab_mover%restartxf==-3)
   end if
   call abihist_bcast(hist_prev,master,comm)

!  If restartxf specifies to reconstruct the history
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-1)then
     ntime=ntime+hist_prev%mxhist
   end if

!  If restartxf specifies to start from the lowest energy
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-2)then
     minE=hist_prev%etot(1)
     minIndex=1
     do ii=1,hist_prev%mxhist
       if(need_verbose) write(std_out,*) 'Iteration:',ii,' Total Energy:',hist_prev%etot(ii)
       if (minE>hist_prev%etot(ii))then
         minE=hist_prev%etot(ii)
         minIndex=ii
       end if
     end do
     if(need_verbose)write(std_out,*) 'The lowest energy occurs at iteration:',minIndex,'etotal=',minE
     acell(:)   =hist_prev%acell(:,minIndex)
     rprimd(:,:)=hist_prev%rprimd(:,:,minIndex)
     xred(:,:)  =hist_prev%xred(:,:,minIndex)
     call abihist_free(hist_prev)
   end if
!  If restarxf specifies to start to the last iteration
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-3)then
     acell(:)   =hist_prev%acell(:,hist_prev%mxhist)
     rprimd(:,:)=hist_prev%rprimd(:,:,hist_prev%mxhist)
     xred(:,:)  =hist_prev%xred(:,:,hist_prev%mxhist)
     call abihist_free(hist_prev)
   end if

 end if !if (ab_mover%restartxf<=0)
#endif

!###########################################################
!### 05. Allocate the hist structure

 iexit=0; timelimit_exit=0
 ncycle=specs%ncycle

 if(ab_mover%ionmov==25.and.scfcv_args%dtset%hmctt>=0)then
   ncycle=scfcv_args%dtset%hmctt
   if(scfcv_args%dtset%hmcsst>0.and.ab_mover%optcell/=0)then
     ncycle=ncycle+scfcv_args%dtset%hmcsst
   endif
 endif

 nhisttot=ncycle*ntime;if (scfcv_args%dtset%nctime>0) nhisttot=nhisttot+1
!AM_2017 New version of the hist, we just store the needed history step not all of them...
 if(specs%nhist/=-1)then
  nhisttot = specs%nhist! We don't need to store all the history
 endif

 call abihist_init(hist,ab_mover%natom,nhisttot,specs%isVused,specs%isARused)
 call abiforstr_ini(preconforstr,ab_mover%natom)

!###########################################################
!### 06. First output before any itime or icycle

!If effective potential is present,
!  forces will be compute with it
 if (present(effective_potential))then
   need_scfcv_cycle = .FALSE.
   if(need_verbose)then
     write(message,'(2a,i2,5a,80a)')&
&     ch10,'=== [ionmov=',ab_mover%ionmov,'] ',trim(specs%method),' with effective potential',&
&     ch10,('=',kk=1,80)
     call wrtout([std_out, ab_out], message)
   end if
   need_elec_eval = .FALSE.
   if(present(scup_dtset))then
     need_elec_eval = scup_dtset%scup_elec_model
   endif
 else
   if(need_verbose)then
     write(message,'(a,a,i2,a,a,a,80a)')&
&     ch10,'=== [ionmov=',ab_mover%ionmov,'] ',specs%method,&
&     ch10,('=',kk=1,80)
     call wrtout([std_out, ab_out], message)
   end if
 end if

!Format for printing on each cycle
 write(fmt,'(a6,i2,a4,i2,a4,i2,a4,i2,a9)')&
& '(a,a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,a,80a)'

!###########################################################
!### 07. Fill the history of the first SCFCV

 if (ab_mover%ionmov==26)then
!Tdep call need to merge with adewandre branch
 else if (ab_mover%ionmov==27)then
   if(present(filename_ddb))then
     fname_ddb = trim(filename_ddb)
   else
     fname_ddb = trim(ab_mover%filnam_ds(3))//'_DDB'
   end if
   INQUIRE(FILE=filename, EXIST=file_exists)

   call generate_training_set(acell,ab_mover%ph_freez_disp_addStrain==1,ab_mover%ph_freez_disp_ampl,&
&                             fname_ddb,hist,ab_mover%natom,ab_mover%ph_freez_disp_nampl,ntime,&
&                             ab_mover%ph_ngqpt,ab_mover%ph_nqshift,ab_mover%ph_freez_disp_option,&
&                             ab_mover%ph_qshift,scfcv_args%dtset%supercell_latt,&
&                             rprimd,ab_mover%mdtemp(2),xred,comm,DEBUG)


   !Fill history with the values of xred, acell and rprimd of the first configuration
   acell(:)   =hist%acell(:,1)
   rprimd(:,:)=hist%rprimd(:,:,1)
   xred(:,:)  =hist%xred(:,:,1)

 else

   !Fill history with the values of xred, acell and rprimd
   call var2hist(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)

   !Fill velocities and ionic kinetic energy
   call vel2hist(ab_mover%amass,hist,vel,vel_cell)
   hist%time(hist%ihist)=zero

 end if
!Decide if prtxfase will be called
 useprtxfase=.FALSE.
 do ii=1,ab_mover%natom
   if (ab_mover%prtatlist(ii)/=0)then
     useprtxfase=.TRUE.
     exit
   end if
 end do

!At beginning no error
 nerr_dilatmx = 0

 ABI_ALLOCATE(xred_prev,(3,scfcv_args%dtset%natom))

!###########################################################
!### 08. Loop for itime (From 1 to ntime)
 quitsum_request = xmpi_request_null

 do itime=1,ntime

   call yaml_iterstart("itime", itime, dev_null, scfcv_args%dtset%use_yaml)

   ! Handle time limit condition.
   if (itime == 1) prev = abi_wtime()
   if (itime  > 1) then
     now = abi_wtime()
     wtime_step = now - prev
     prev = now
     write(message,*)sjoin("mover: previous time step took ",sec2str(wtime_step))
     if(need_verbose)call wrtout(std_out, message)
     if (have_timelimit_in(MY_NAME)) then
       if (itime > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then
           write(message,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),&
&           ". Will exit itime loop in mover."
           if(need_verbose)MSG_COMMENT(message)
           if(need_verbose)call wrtout(ab_out, message)
           timelimit_exit = 1
           exit
         end if
       end if

       my_quit = 0; if (now - get_start_time() + 2.15 * wtime_step > get_timelimit()) my_quit = 1
       call xmpi_isum(my_quit,quitsum_async,comm,quitsum_request,ierr)
     end if
   end if

   skipcycle=.FALSE.
#if defined HAVE_LOTF
   if(ab_mover%ionmov==23 .and. .not. lotf_extrapolation(itime)) skipcycle=.True.
#endif

!  ###########################################################
!  ### 09. Loop for icycle (From 1 to ncycle)
   do icycle=1,ncycle

     call yaml_iterstart("icycle", icycle, dev_null, scfcv_args%dtset%use_yaml)

     itime_hist = (itime-1)*ncycle + icycle ! Store the time step in the history

!    ###########################################################
!    ### 10. Output for each icycle (and itime)
     if(need_verbose)then
       write(message,fmt)&
&       ch10,'--- Iteration: (',itime,'/',ntime,') Internal Cycle: (',icycle,'/',ncycle,')',ch10,('-',kk=1,80)
        call wrtout([std_out, ab_out], message)
     end if
     if (useprtxfase) then
       call prtxfase(ab_mover,hist,itime_hist,std_out,mover_BEFORE)
     end if

     xred_prev(:,:)=xred(:,:)
     rprimd_prev(:,:)=rprimd(:,:)

!    ###########################################################
!    ### 11. Symmetrize atomic coordinates over space group elements

     call symmetrize_xred(scfcv_args%indsym,ab_mover%natom,&
&     scfcv_args%dtset%nsym,scfcv_args%dtset%symrel,scfcv_args%dtset%tnons,xred)

     change=any(xred(:,:)/=xred_prev(:,:))
     if (change)then
       hist%xred(:,:,hist%ihist)=xred(:,:)
       if(need_verbose) then
         write(std_out,*) 'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED'
         write(std_out,*) 'DIFFERENCES:'
         do kk=1,ab_mover%natom
           write(std_out,*) xred(:,kk)-xred_prev(:,kk)
         end do
       end if
       xred_prev(:,:)=xred(:,:)
     end if

!    ###########################################################
!    ### 12. => Call to SCFCV routine and fill history with forces
     if (need_verbose) then
       if (need_scfcv_cycle) then
         write(message,'(a,3a,33a,44a)')&
          ch10,('-',kk=1,3),'SELF-CONSISTENT-FIELD CONVERGENCE',('-',kk=1,44)
       else
         write(message,'(a,3a,33a,44a)')&
          ch10,('-',kk=1,3),'EFFECTIVE POTENTIAL CALCULATION',('-',kk=1,44)
       end if
       call wrtout([std_out, ab_out], message)
     end if

     if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then

       call abihist_compare_and_copy(hist_prev,hist,ab_mover%natom,similar,tol8,specs%nhist==nhisttot)
       hist_prev%ihist=hist_prev%ihist+1

     else
       scfcv_args%ndtpawuj=0
       iapp=itime
       if(icycle>1.and.icycle/=ncycle) iapp=-1
       if(itime==1 .and. icycle/=ncycle ) iapp=-icycle-1
       if (ab_mover%ionmov==14.and.(icycle<ncycle)) iapp=-1

#if defined HAVE_LOTF
       if (ab_mover%ionmov/=23 .or.(lotf_extrapolation(itime).and.(icycle/=1.or.itime==1)))then
#endif
         !call scfcv_new2(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)

         !WVL - reformat the wavefunctions in the case of xred != xred_old
         if (scfcv_args%dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
          !Before running scfcv, on non-first geometry step iterations,
          ! we need to reformat the wavefunctions, taking into acount the new
          ! coordinates. We prepare to change rhog (to be removed) and rhor.
           ABI_DEALLOCATE(rhog)
           ABI_DEALLOCATE(rhor)
           call wvl_wfsinp_reformat(scfcv_args%dtset, scfcv_args%mpi_enreg,&
&           scfcv_args%psps, rprimd, scfcv_args%wvl, xred, xred_old)
           scfcv_args%nfftf = scfcv_args%dtset%nfft
           ABI_ALLOCATE(rhog,(2, scfcv_args%dtset%nfft))
           ABI_ALLOCATE(rhor,(2, scfcv_args%dtset%nfft))
           call wvl_mkrho(scfcv_args%dtset, scfcv_args%irrzon, scfcv_args%mpi_enreg,&
&           scfcv_args%phnons, rhor,scfcv_args%wvl%wfs,scfcv_args%wvl%den)
         end if

!        MAIN CALL TO SELF-CONSISTENT FIELD ROUTINE
         if (need_scfcv_cycle) then

           call dtfil_init_time(dtfil,iapp)
           call scfcv_run(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
           if (conv_retcode == -1) then
               message = "Scf cycle returned conv_retcode == -1 (timelimit is approaching), this should not happen inside mover"
               MSG_WARNING(message)
           end if

         else
!          For monte carlo don't need to recompute energy here
!          (done in pred_montecarlo)
           name_file='MD_anharmonic_terms_energy.dat'
             if(itime == 1 .and. ab_mover%restartxf==-3)then
               call effective_potential_file_mapHistToRef(effective_potential,hist,comm,scfcv_args%dtset%iatfix,need_verbose) ! Map Hist to Ref to order atoms
               xred(:,:) = hist%xred(:,:,1) ! Fill xred with new ordering
               hist%ihist = 1
             end if

#if defined DEV_MS_SCALEUP
           !If we a SCALE UP effective electron model give the iteration and set print-options
           if(need_elec_eval)then
              call global_set_parent_iter(itime)
              ! Set all print options to false.
              call global_set_print_parameters(geom=.FALSE.,eigvals=.FALSE.,eltic=.FALSE.,&
&                      orbocc=.FALSE.,bands=.FALSE.)
              if(itime == 1 .or. modulo(itime,scup_dtset%scup_printniter) == 0)then
                 call global_set_print_parameters(scup_dtset%scup_printgeom,scup_dtset%scup_printeigv,scup_dtset%scup_printeltic,&
&                         scup_dtset%scup_printorbocc,scup_dtset%scup_printbands)
              end if
           end if
#endif

           call effective_potential_evaluate( &
&           effective_potential,scfcv_args%results_gs%etotal,scfcv_args%results_gs%fcart,scfcv_args%results_gs%fred,&
&           scfcv_args%results_gs%strten,ab_mover%natom,rprimd,xred=xred,verbose=need_verbose,&
&           filename=name_file,elec_eval=need_elec_eval)

!          Check if the simulation did not diverge...
           if(itime > 3 .and.ABS(scfcv_args%results_gs%etotal - hist%etot(1)) > 1E5)then
!            We set to false the flag corresponding to the bound
             effective_potential%anharmonics_terms%bounded = .FALSE.
             if(need_verbose.and.me==master)then
               MSG_WARNING("The simulation is diverging, please check your effective potential")
             end if
!            Set the flag to finish the simulation
             iexit=1
             stat4xml="Failed"
           else
!            We set to true the flag corresponding to the bound
             effective_potential%anharmonics_terms%bounded = .TRUE.
           end if
         end if
#if defined HAVE_LOTF
       end if
#endif
!      ANOMALOUS SITUATION
!      This is the only case where rprimd could change inside scfcv
!      It generates an weird condition, we start with a certain
!      value for rprimd before scfcv and after we finish with
!      a different value.
!      Notice that normally scfcv should not change rprimd
!      And even worse if optcell==0
!      The solution here is to recompute acell and store these value
!      in the present record even if initially it was not exactly
!      the value entering in scfcv
!      One test case with these condition is bigdft/t10
       if (any(rprimd(:,:)/=rprimd_prev(:,:))) then
         hist%acell(:,hist%ihist)=acell(:)
         hist%rprimd(:,:,hist%ihist)=rprimd(:,:)
       end if

!      ANOMALOUS SITUATIONS
!      * In ionmov 4 & 5 xred could change inside SCFCV
!      So we need to take the values from the output
!
!      * Inside scfcv_core.F90 there is a call to symmetrize_xred.F90
!      for the first SCF cycle symmetrize_xred could change xred
       if (ab_mover%ionmov<10)then
         change=any(xred(:,:)/=xred_prev(:,:))
         if (change)then
           hist%xred(:,:,hist%ihist)=xred(:,:)
           if(need_verbose)then
             write(std_out,*) 'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED AFTER SCFCV'
             write(std_out,*) 'DIFFERENCES:'
             do kk=1,ab_mover%natom
               write(std_out,*) xred(:,kk)-xred_prev(:,kk)
             end do
           end if
         end if !if (change)
       end if

!      Fill velocities and ionic kinetic energy
       call vel2hist(ab_mover%amass,hist,vel,vel_cell)

       hist%fcart(:,:,hist%ihist)=scfcv_args%results_gs%fcart(:,:)
       hist%strten(:,hist%ihist) =scfcv_args%results_gs%strten(:)
       hist%etot(hist%ihist)     =scfcv_args%results_gs%etotal
       hist%entropy(hist%ihist)  =scfcv_args%results_gs%energies%entropy
       hist%time(hist%ihist)     =real(itime,kind=dp)

!      !######################################################################
     end if ! if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then

!    Store trajectory in xfh file
     if((ab_xfh%nxfh==0.or.itime/=1)) then
       ABI_ALLOCATE(fred_corrected,(3,scfcv_args%dtset%natom))
       call fcart2fred(hist%fcart(:,:,hist%ihist),fred_corrected,rprimd,ab_mover%natom)
!      Get rid of mean force on whole unit cell,
!       but only if no generalized constraints are in effect
       if (ab_mover%nconeq==0)then
         do ii=1,3
           if (ii/=3.or.ab_mover%jellslab==0) then
             favg=sum(fred_corrected(ii,:))/dble(ab_mover%natom)
             fred_corrected(ii,:)=fred_corrected(ii,:)-favg
           end if
         end do
       end if
       if (ncycle<10.and.ab_mover%restartxf>=0) then
         do ii=1,3
           rprim(ii,1:3)=rprimd(ii,1:3)/acell(1:3)
         end do

!        The size of ab_xfh%xfhist is to big for very large supercell.
!        Call it only for specific ionmov
         if(any((/2,3,10,11,22/)==ab_mover%ionmov)) then
           call xfh_update(ab_xfh,acell,fred_corrected,ab_mover%natom,rprim,hist%strten(:,hist%ihist),xred)
         end if
       end if
       ABI_DEALLOCATE(fred_corrected)
     end if

!    ###########################################################
!    ### 13. Write the history into the _HIST file
!    ###

#if defined HAVE_NETCDF
     if (need_writeHIST.and.me==master) then
       ifirst=merge(0,1,(itime>1.or.icycle>1))
       call write_md_hist(hist,filename,ifirst,itime_hist,ab_mover%natom,scfcv_args%dtset%nctime,&
&       ab_mover%ntypat,ab_mover%typat,amu_curr,ab_mover%znucl,ab_mover%dtion,scfcv_args%dtset%mdtemp)
     end if
#endif

!    ###########################################################
!    ### 14. Output after SCFCV
     if(need_verbose.and.need_scfcv_cycle)then
       write(message,'(a,3a,a,72a)')ch10,('-',kk=1,3),'OUTPUT',('-',kk=1,71)
       call wrtout([std_out, ab_out], message)
     end if
     if (useprtxfase) then
       call prtxfase(ab_mover,hist,itime_hist,ab_out,mover_AFTER)
       call prtxfase(ab_mover,hist,itime_hist,std_out,mover_AFTER)
     end if

!    ###########################################################
!    ### 15. => Test Convergence of forces and stresses

     if (itime==ntime.and.icycle==ncycle)then
       iexit=1
       stat4xml="Failed"
     else
       stat4xml="Succeded"
     end if

!    Only if convergence is needed
     if(specs%isFconv)then
       if ((ab_mover%ionmov/=4.and.ab_mover%ionmov/=5).or.mod(itime,2)==1)then
         if (scfcv_args%dtset%tolmxf/=0)then
           call fconv(hist%fcart(:,:,hist%ihist),&
&           scfcv_args%dtset%iatfix, &
&           iexit,itime,&
&           ab_mover%natom,&
&           ntime,&
&           ab_mover%optcell,&
&           scfcv_args%dtset%strfact,&
&           scfcv_args%dtset%strtarget,&
&           hist%strten(:,hist%ihist),&
&           scfcv_args%dtset%tolmxf)
         else
           call erlxconv(hist,iexit,itime,itime_hist,ntime,scfcv_args%dtset%tolmxde)
         end if
       end if
     end if

     if (itime==ntime.and.icycle==ncycle) iexit=1

!    ###########################################################
!    ### 16. => Precondition forces, stress and energy
!    ### 17. => Call to each predictor

     call precpred_1geo(ab_mover,ab_xfh,amu_curr,deloc,&
&     scfcv_args%dtset%chkdilatmx,&
&     scfcv_args%mpi_enreg%comm_cell,&
&     scfcv_args%dtset%dilatmx,dtfil%filnam_ds(4),&
&     hist,scfcv_args%dtset%hmctt,&
&     icycle,iexit,itime,mttk_vars,&
&     scfcv_args%dtset%nctime,ncycle,nerr_dilatmx,scfcv_args%dtset%npsp,ntime,&
&     scfcv_args%dtset%rprimd_orig,skipcycle,&
&     scfcv_args%dtset%usewvl)

!    Write MOLDYN netcdf and POSABIN files (done every dtset%nctime time step)
!    This file is not created for multibinit run
     if(need_scfcv_cycle .and. (ab_mover%ionmov/=23 .or. icycle==1))then
       if (scfcv_args%dtset%nctime>0) then
         jj=itime; if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1) jj=jj-hist_prev%mxhist
         if (jj>0) then
           option=3
           ihist_prev = abihist_findIndex(hist,-1)
           call wrt_moldyn_netcdf(ab_mover%amass,scfcv_args%dtset,jj,option,dtfil%fnameabo_moldyn,&
&           scfcv_args%mpi_enreg,scfcv_args%results_gs,&
&           hist%rprimd(:,:,ihist_prev),dtfil%unpos,hist%vel(:,:,hist%ihist),&
&           hist%xred(:,:,ihist_prev))
         end if
         if (iexit==1) hist%ihist=ihist_prev
       end if
     end if
     if(iexit/=0) exit

!    ###########################################################
!    ### 18. Use the history  to extract the new values
!    ###     acell, rprimd and xred

     call hist2var(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)

     if(ab_mover%optcell/=0)then

       call matr3inv(rprimd,gprimd)

!      If metric has changed since the initialization, update the Ylm's
       if (scfcv_args%psps%useylm==1)then
         option=0;
         if (scfcv_args%dtset%iscf>0) option=1
         call initylmg(gprimd,&
&         scfcv_args%kg,&
&         scfcv_args%dtset%kptns,&
&         scfcv_args%dtset%mkmem,&
&         scfcv_args%mpi_enreg,&
&         scfcv_args%psps%mpsang,&
&         scfcv_args%dtset%mpw,&
&         scfcv_args%dtset%nband,&
&         scfcv_args%dtset%nkpt,&
&         scfcv_args%npwarr,&
&         scfcv_args%dtset%nsppol,&
&         option,rprimd,&
&         scfcv_args%ylm,&
&         scfcv_args%ylmgr)
       end if

     end if

     vel(:,:)=hist%vel(:,:,hist%ihist)

!    vel_cell(3,3)= velocities of cell parameters
!    Not yet used here but compute it for consistency
     vel_cell(:,:)=zero
     if (ab_mover%ionmov==13 .and. hist%mxhist >= 2) then
       if (itime_hist>2) then
         ihist_prev2 = abihist_findIndex(hist,-2)
         vel_cell(:,:)=(hist%rprimd(:,:,hist%ihist)- hist%rprimd(:,:,ihist_prev2))/(two*ab_mover%dtion)
       else if (itime_hist>1) then
         ihist_prev = abihist_findIndex(hist,-1)
         vel_cell(:,:)=(hist%rprimd(:,:,hist%ihist)-hist%rprimd(:,:,ihist_prev))/(ab_mover%dtion)
       end if
     end if

!    This is needed for some compilers such as
!    pathscale, g95, xlf that do not exit
!    from a loop if you change the upper limit
!    inside
     if (icycle>=ncycle .and. scfcv_args%mpi_enreg%me == 0) then
       if(need_verbose)write(std_out,*) 'EXIT:',icycle,ncycle
       exit
     end if


     if (need_verbose) then
       write(message,*) 'ICYCLE',icycle,skipcycle
       call wrtout(std_out,message)
       write(message,*) 'NCYCLE',ncycle
       call wrtout(std_out,message)
     end if
     if (skipcycle) exit

     !write(std_out,*)' m_mover : will call precpred_1geo'


!    ###########################################################
!    ### 19. End loop icycle

   end do ! do icycle=1,ncycle

   if(iexit/=0)exit

!  ###########################################################
!  ### 20. End loop itime

 end do ! do itime=1,ntime

 ! Call fconv here if we exited due to wall time limit.
 if (timelimit_exit==1 .and. specs%isFconv) then
   iexit = timelimit_exit
   ntime = itime-1
   ihist_prev = abihist_findIndex(hist,-1)
   if ((ab_mover%ionmov/=4.and.ab_mover%ionmov/=5)) then
     if (scfcv_args%dtset%tolmxf/=0)then
       call fconv(hist%fcart(:,:,ihist_prev),&
&       scfcv_args%dtset%iatfix, &
&       iexit, itime,&
&       ab_mover%natom,&
&       ntime,&
&       ab_mover%optcell,&
&       scfcv_args%dtset%strfact,&
&       scfcv_args%dtset%strtarget,&
&       hist%strten(:,ihist_prev),&
&       scfcv_args%dtset%tolmxf)
     else
       call erlxconv(hist,iexit,itime,itime_hist,ntime,scfcv_args%dtset%tolmxde)
     end if
   end if
 end if

 ! Avoid pending requests if itime == ntime.
 call xmpi_wait(quitsum_request,ierr)

!###########################################################
!### 21. Set the final values of xred with the last
!###     computed values (not the last predicted)

 call hist2var(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)
 vel(:,:)=hist%vel(:,:,hist%ihist)

 if (DEBUG .and. ab_mover%ionmov==1)then
   write (std_out,*) 'vel'
   do kk=1,ab_mover%natom
     write (std_out,*) hist%vel(:,kk,hist%ihist)
   end do
 end if

!###########################################################
!### 22. XML Output at the end

!XML output of the status
 if (scfcv_args%mpi_enreg%me == 0 .and. scfcv_args%dtset%prtxml == 1) then
   write(ab_xml_out, "(3a)") '    <geometryMinimisation type="',trim(specs%type4xml),'">'
   write(ab_xml_out, "(5a)") '      <status cvState="',trim(stat4xml) ,'" stop-criterion="',trim(specs%crit4xml),'" />'
   write(ab_xml_out, "(3a)") '    </geometryMinimisation>'
 end if

!###########################################################
!### 23. Deallocate hist and ab_mover datatypes

!This call is needed to free an internal matrix. However, this is not optimal ...
!One should instead have a datastructure associated with the preconditioner...
 if (ab_mover%goprecon>0)then
   call prec_simple(ab_mover,preconforstr,hist,1,1,1)
 end if

 if (ab_mover%ionmov==13)then
   call mttk_fin(mttk_vars)
 end if

 if (ab_mover%ionmov==10 .or. ab_mover%ionmov==11)then
   call delocint_fin(deloc)
 end if

 ABI_DEALLOCATE(xred_prev)

 call abihist_free(hist)
 call abihist_free(hist_prev)

 call abimover_destroy(ab_mover)
 call abiforstr_fin(preconforstr)

contains
!!***

!!****f* ABINIT/fconv
!!
!! NAME
!! fconv
!!
!! FUNCTION
!! Check maximal absolute value of force (hartree/bohr) against
!! input tolerance; if below tolerance, return iexit=1.
!! Takes into account the fact that the Broyden (or moldyn) step
!! might be the last one (last itime), to print eventually modified message.
!! Stresses are also included in the check, provided that optcell/=0.
!! If optcell=1, takes only the trace into account
!!    optcell=2, takes all components into account
!!    optcell=3, takes traceless stress into account
!!    optcell=4, takes sigma(1 1) into account
!!    optcell=5, takes sigma(2 2) into account
!!    optcell=6, takes sigma(3 3) into account
!!    optcell=7, takes sigma(2,2),(2,3) and (3 3) into account
!!    optcell=8, takes sigma(1,1),(1,3) and (3 3) into account
!!    optcell=9, takes sigma(1,1),(1,2) and (2 2) into account
!! In the case of stresses, target the tensor strtarget, and
!! take into account the factor strfact
!!
!! INPUTS
!!  fcart(3,natom)= forces on atoms in hartree/bohr in cartesian coordinates
!!  iatfix(3,natom)=1 for frozen atom, 0 for unfrozen
!!  itime=current number of Broyden/Moldyn iterations
!!  natom=number of atoms in unit cell
!!  ntime=maximum number of Broyden/Moldyn iterations allowed
!!  optcell=option for taking stresses into account (see above)
!!  strfact=factor that multiplies the stresses when they are compared to forces.
!!  strtarget(6)=components of the target stress tensor (hartree/bohr^3)
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  tolmxf=tolerance on maximal absolute value of components of forces
!!
!! OUTPUT
!!  writes to unit std_out and to ab_out, and returns
!!
!! SIDE EFFECTS
!! Input/Output
!!  at input  : iexit=  0 if not the last itime,  1 if the last itime
!!  at output : iexit=  0 if not below tolerance, 1 if below tolerance
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine fconv(fcart,iatfix,iexit,itime,natom,ntime,optcell,strfact,strtarget,strten,tolmxf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ntime,optcell
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: strfact,tolmxf
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fcart(3,natom),strtarget(6),strten(6)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,istr
 real(dp) :: fmax,strdiag
 character(len=500) :: message
!arrays
 real(dp) :: dstr(6)

! *************************************************************************

!Compute maximal component of forces, EXCLUDING any fixed components
 fmax=zero
 do iatom=1,natom
   do idir=1,3
     if (iatfix(idir,iatom) /= 1) then
       if( abs(fcart(idir,iatom)) >= fmax ) fmax=abs(fcart(idir,iatom))
     end if
   end do
 end do

 dstr(:)=strten(:)-strtarget(:)

!Eventually take into account the stress
 if(optcell==1)then
   strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   if(abs(strdiag)*strfact >= fmax ) fmax=abs(strdiag)*strfact
 else if(optcell==2)then
   do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
   end do
 else if(optcell==3)then
!  Must take away the trace from diagonal elements
   strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   do istr=1,3
     if(abs(dstr(istr)-strdiag)*strfact >= fmax ) fmax=abs(dstr(istr)-strdiag)*strfact
   end do
   do istr=4,6
     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
   end do
 else if(optcell==4 .or. optcell==5 .or. optcell==6)then
   if(abs(dstr(optcell-3))*strfact >= fmax ) fmax=abs(dstr(optcell-3))*strfact
 else if(optcell==7)then
   if(abs(dstr(2))*strfact >= fmax ) fmax=abs(dstr(2))*strfact
   if(abs(dstr(3))*strfact >= fmax ) fmax=abs(dstr(3))*strfact
   if(abs(dstr(4))*strfact >= fmax ) fmax=abs(dstr(4))*strfact
 else if(optcell==8)then
   if(abs(dstr(1))*strfact >= fmax ) fmax=abs(dstr(1))*strfact
   if(abs(dstr(3))*strfact >= fmax ) fmax=abs(dstr(3))*strfact
   if(abs(dstr(5))*strfact >= fmax ) fmax=abs(dstr(5))*strfact
 else if(optcell==9)then
   if(abs(dstr(1))*strfact >= fmax ) fmax=abs(dstr(1))*strfact
   if(abs(dstr(2))*strfact >= fmax ) fmax=abs(dstr(2))*strfact
   if(abs(dstr(6))*strfact >= fmax ) fmax=abs(dstr(6))*strfact
 end if

 if (fmax<tolmxf) then
   write(message, '(a,a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&   ' At Broyd/MD step',itime,', gradients are converged : ',ch10,&
&   '  max grad (force/stress) =',fmax,' < tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
   call wrtout([std_out, ab_out], message)
   iexit=1
 else
   if(iexit==1)then
     write(message, '(a,a,a,a,i5,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&     ' fconv : WARNING -',ch10,&
&     '  ntime=',ntime,' was not enough Broyd/MD steps to converge gradients: ',ch10,&
&     '  max grad (force/stress) =',fmax,' > tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
     call wrtout([std_out, ab_out], message)

     write(std_out,"(8a)")ch10,&
&     "--- !RelaxConvergenceWarning",ch10,&
&     "message: | ",ch10,TRIM(indent(message)),ch10,&
&     "..."

   else
     write(message, '(a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) &
&     ' fconv : at Broyd/MD step',itime,', gradients have not converged yet. ',ch10,&
&     '  max grad (force/stress) =',fmax,' > tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
     call wrtout(std_out,message,'COLL')
   end if
   iexit=0
 end if

end subroutine fconv
!!***

!!****f* ABINIT/erlxconv
!! NAME
!!  erlxconv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine erlxconv(hist,iexit,itime,itime_hist,ntime,tolmxde)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,itime_hist,ntime
 integer,intent(inout) :: iexit
 real(dp), intent(in) :: tolmxde
!arrays
 type(abihist),intent(inout) :: hist

!Local variables-------------------------------
 integer :: ihist,ihist_prev,ihist_prev2
 real(dp) :: ediff1,ediff2,maxediff
 character(len=500) :: message
! *************************************************************************

 if (itime_hist<3) then
   write(message, '(a,a,a)' ) ch10,&
&   ' erlxconv : minimum 3 Broyd/MD steps to check convergence of energy in relaxations',ch10
   call wrtout(std_out,message,'COLL')
 else
   ihist = hist%ihist
   ihist_prev  = abihist_findIndex(hist,-1)
   ihist_prev2 = abihist_findIndex(hist,-2)
   ediff1 = hist%etot(ihist) - hist%etot(ihist_prev)
   ediff2 = hist%etot(ihist) - hist%etot(ihist_prev2)
   if ((abs(ediff1)<tolmxde).and.(abs(ediff2)<tolmxde)) then
     write(message, '(a,a,i4,a,a,a,a,a,es11.4,a,a)' ) ch10,&
&     ' At Broyd/MD step',itime,', energy is converged : ',ch10,&
&     '  the difference in energy with respect to the two ',ch10,&
&     '  previous steps is < tolmxde=',tolmxde,' ha',ch10
     call wrtout([std_out, ab_out], message)
     iexit=1
   else
     maxediff = max(abs(ediff1),abs(ediff2))
     if(iexit==1)then
       write(message, '(a,a,a,a,i5,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&       ' erlxconv : WARNING -',ch10,&
&       '  ntime=',ntime,' was not enough Broyd/MD steps to converge energy: ',ch10,&
&       '  max difference in energy =',maxediff,' > tolmxde=',tolmxde,' ha',ch10
       call wrtout([std_out, ab_out], message)

       write(std_out,"(8a)")ch10,&
&       "--- !RelaxConvergenceWarning",ch10,&
&       "message: | ",ch10,TRIM(indent(message)),ch10,&
&       "..."
     else
       write(message, '(a,a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&       ' erlxconv : at Broyd/MD step',itime,', energy has not converged yet. ',ch10,&
&       '  max difference in energy=',maxediff,' > tolmxde=',tolmxde,' ha',ch10
       call wrtout(std_out,message,'COLL')
     end if
   end if
 end if

end subroutine erlxconv
!!***

end subroutine mover
!!***

!!****f* ABINIT/prtxfase
!!
!! NAME
!! prtxfase
!!
!! FUNCTION
!! Print the values of xcart (X), forces (F) acell (A), Stresses (S), and energy (E)
!! All values come from the history hist
!! Also compute and print max and rms forces.
!! Also compute absolute and relative differences with previous calculation
!!
!! INPUTS
!! ab_mover<type abimover>=Subset of dtset only related with
!!          |                 movement of ions and acell, contains:
!!          | dtion:  Time step
!!          ! natom:  Number of atoms
!!          | vis:    viscosity
!!          | iatfix: Index of atoms and directions fixed
!!          | amass:  Mass of ions
!! hist<type abihist>=Historical record of positions, forces,
!!                               stresses, cell and energies,
!! itime= time step
!! iout=unit number for printing
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE

subroutine prtxfase(ab_mover,hist,itime,iout,pos)

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(in),target :: hist
 integer,intent(in) :: itime,iout
 integer,intent(in) :: pos
!arrays

!Local variables-------------------------------
!scalars
 integer :: jj,kk,unfixd,iprt
 real(dp) :: val_max,val_rms,ucvol ! Values maximal and RMS, Volume of Unitary cell
 real(dp) :: dEabs,dErel ! Diff of energy absolute and relative
 real(dp) :: ekin
 real(dp) :: angle(3),rmet(3,3)
!character(len=80*(max(ab_mover%natom,3)+1)) :: msg
!MGNAG: This is not very safe. One should use line-based output istead of appending chars
! and then outputting everything! For the time being I use this temporary hack to solve the problem with NAG
 character(len=max(80*(max(ab_mover%natom,3)+1),50000)) :: msg
 character(len=18)   :: fmt1
 logical :: prtallatoms
!arrays
 logical :: atlist(ab_mover%natom)
 real(dp),allocatable :: fred(:,:),xcart(:,:)
 real(dp),pointer :: acell(:),fcart(:,:),rprimd(:,:),strten(:),vel(:,:),xred(:,:)

! ***********************************************************

 fmt1='(a,a,1p,3e22.14)'

!##########################################################
!### 1. Organize list of atoms to print

 prtallatoms=.TRUE.
 do kk=1,ab_mover%natom
   if (ab_mover%prtatlist(kk)/=kk) prtallatoms=.FALSE.
 end do

 atlist(:)=.FALSE.
 do iprt=1,ab_mover%natom
   if (ab_mover%prtatlist(iprt)>0.and.ab_mover%prtatlist(iprt)<=ab_mover%natom) atlist(ab_mover%prtatlist(iprt))=.TRUE.
 end do

 acell  => hist%acell(:,hist%ihist)
 rprimd => hist%rprimd(:,:,hist%ihist)
 xred   => hist%xred(:,:,hist%ihist)
 fcart  => hist%fcart(:,:,hist%ihist)
 strten => hist%strten(:,hist%ihist)
 vel    => hist%vel(:,:,hist%ihist)

!###########################################################
!### 1. Positions

 ABI_ALLOCATE(xcart,(3,ab_mover%natom))
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

 write(msg, '(a,a)' )ch10,' Cartesian coordinates (xcart) [bohr]'
 call prtnatom(atlist,iout,msg,ab_mover%natom,prtallatoms,xcart)

 write(msg, '(a)' )' Reduced coordinates (xred)'
 call prtnatom(atlist,iout,msg,ab_mover%natom,prtallatoms,xred)

 ABI_DEALLOCATE(xcart)

!###########################################################
!### 2. Forces

 if(pos==mover_AFTER)then

   ABI_ALLOCATE(fred,(3,ab_mover%natom))
   call fcart2fred(fcart,fred,rprimd,ab_mover%natom)

!  Compute max |f| and rms f,
!  EXCLUDING the components determined by iatfix
   val_max=0.0_dp
   val_rms=0.0_dp
   unfixd=0
   do kk=1,ab_mover%natom
     do jj=1,3
       if (ab_mover%iatfix(jj,kk) /= 1) then
         unfixd=unfixd+1
         val_rms=val_rms+fcart(jj,kk)**2
         val_max=max(val_max,abs(fcart(jj,kk)**2))
       end if
     end do
   end do
   if ( unfixd /= 0 ) val_rms=sqrt(val_rms/dble(unfixd))

   write(msg, '(a,1p,2e12.5,a)' ) ' Cartesian forces (fcart) [Ha/bohr]; max,rms=',sqrt(val_max),val_rms,' (free atoms)'
   call prtnatom(atlist,iout,msg,ab_mover%natom,prtallatoms,fcart)

   write(msg, '(a)' )' Reduced forces (fred)'
   call prtnatom(atlist,iout,msg,ab_mover%natom,prtallatoms,fred)
   ABI_DEALLOCATE(fred)
 end if

!###########################################################
!### 3. Velocities

!Only if the velocities are being used
 if (hist%isVused)then
!  Only if velocities are recorded in a history
   if (allocated(hist%vel))then
!    Compute max |v| and rms v,
!    EXCLUDING the components determined by iatfix
     val_max=0.0_dp
     val_rms=0.0_dp
     unfixd=0
     do kk=1,ab_mover%natom
       do jj=1,3
         if (ab_mover%iatfix(jj,kk) /= 1) then
           unfixd=unfixd+1
           val_rms=val_rms+vel(jj,kk)**2
           val_max=max(val_max,abs(vel(jj,kk)**2))
         end if
       end do
     end do
     if ( unfixd /= 0 ) val_rms=sqrt(val_rms/dble(unfixd))

     write(msg, '(a,1p,2e12.5,a)' ) &
&     ' Cartesian velocities (vel) [bohr*Ha/hbar]; max,rms=',sqrt(val_max),val_rms,' (free atoms)'
     call prtnatom(atlist,iout,msg,ab_mover%natom,prtallatoms,vel)

!    Compute the ionic kinetic energy (no cell shape kinetic energy yet)
     ekin=0.0_dp
     do kk=1,ab_mover%natom
       do jj=1,3
!        Warning : the fixing of atoms is implemented in reduced
!        coordinates, so that this expression is wrong
         if (ab_mover%iatfix(jj,kk) == 0) then
           ekin=ekin+0.5_dp*ab_mover%amass(kk)*vel(jj,kk)**2
         end if
       end do
     end do
     write(msg, '(a,1p,e22.14,a)' )' Kinetic energy of ions (ekin) [Ha]=',ekin
     call wrtout(iout,msg,'COLL')
   end if
 end if

!###########################################################
!### 3. ACELL

!Only if the acell is being used
 if (hist%isARused)then
!  Only if acell is recorded in a history
   if (allocated(hist%acell))then
     write(msg, '(a)' )' Scale of Primitive Cell (acell) [bohr]'
     write(msg,fmt1)TRIM(msg),ch10,acell(:)
     call wrtout(iout,msg,'COLL')
   end if
 end if

!###########################################################
!### 4. RPRIMD

!Only if the acell is being used
 if (hist%isARused)then
!  Only if rprimd is recorded in a history
   if (allocated(hist%rprimd))then
     write(msg, '(a)' )' Real space primitive translations (rprimd) [bohr]'
     do kk=1,3
       write(msg,fmt1)TRIM(msg),ch10,rprimd(:,kk)
     end do
     call wrtout(iout,msg,'COLL')
   end if
 end if

!###########################################################
!### 5. Unitary cell volume

 if (ab_mover%optcell/=0)then

   ucvol=&
&   rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
&   rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
&   rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

   write(msg, '(a,1p,e22.14)' )' Unitary Cell Volume (ucvol) [Bohr^3]=',ucvol
   call wrtout(iout,msg,'COLL')

!  ###########################################################
!  ### 5. Angles and lengths

!  Compute real space metric.
   rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0

   write(msg, '(a)' )' Angles (23,13,12)= [degrees]'
   write(msg,fmt1)TRIM(msg),ch10,angle(:)
   call wrtout(iout,msg,'COLL')

   write(msg, '(a)' ) ' Lengths [Bohr]'
   write(msg,fmt1)TRIM(msg),ch10,sqrt(rmet(1,1)),sqrt(rmet(2,2)),sqrt(rmet(3,3))
   call wrtout(iout,msg,'COLL')

!  ###########################################################
!  ### 5. Stress Tensor

   if(pos==mover_AFTER)then
!    Only if strten is recorded in a history
     if (allocated(hist%strten))then

       write(msg, '(a)' ) ' Stress tensor in cartesian coordinates (strten) [Ha/bohr^3]'

       write(msg,fmt1)TRIM(msg),ch10,strten(1),strten(6),strten(5)
       write(msg,fmt1)TRIM(msg),ch10,strten(6),strten(2),strten(4)
       write(msg,fmt1)TRIM(msg),ch10,strten(5),strten(4),strten(3)
       call wrtout(iout,msg,'COLL')
     end if
   end if
 end if

!###########################################################
!### 6. Energy

 if(pos==mover_AFTER)then
   write(msg, '(a,1p,e22.14)' )' Total energy (etotal) [Ha]=',hist%etot(hist%ihist)

   if (itime>1)then
     jj = abihist_findIndex(hist,-1)
     dEabs=hist%etot(hist%ihist)-hist%etot(jj)
     dErel=2*dEabs/(abs(hist%etot(hist%ihist))+abs(hist%etot(jj)))
     write(msg, '(a,a,a,a)' )TRIM(msg),ch10,ch10,' Difference of energy with previous step (new-old):'
     write(msg, '(a,a,10a,a,1p,e12.5,a,10a,a,1p,e12.5)')&
&     TRIM(msg),ch10,&
&     (' ',jj=1,10),' Absolute (Ha)=',dEabs,ch10,&
&     (' ',jj=1,10),' Relative     =',dErel
   end if
   call wrtout(iout,msg,'COLL')
 end if

 contains
!!***

!!****f* ABINIT/gettag
!!
!! NAME
!! gettag
!!
!! FUNCTION
!! Set the tag associated to each atom,
!!
!! INPUTS
!! prtallatoms = Logical for PRTint ALL ATOMS
!! atlist      = ATom LIST
!! index       = index for each atom
!! natom       = Number of ATOMs
!!
!! OUTPUT
!!  tag = The string to put for each atom
!!
!! PARENTS
!!      prtxfase
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE

subroutine gettag(atlist,index,natom,prtallatoms,tag)

!Arguments ------------------------------------
!scalars
  logical,intent(in) :: prtallatoms
  integer,intent(in) :: natom
  logical,intent(in) :: atlist(natom)
  integer,intent(in) :: index
  character(len=7),intent(out)   :: tag

!Local variables -------------------------

! *********************************************************************
!The numbering will be from (1) to (9999)

 if (prtallatoms)then
   tag=''
 elseif (atlist(index)) then
   if (natom<10) then
     write(tag, '(a,I1.1,a)') ' (',index,')'
   elseif (natom<100) then
     write(tag, '(a,I2.2,a)') ' (',index,')'
   elseif (natom<1000) then
     write(tag, '(a,I3.3,a)') ' (',index,')'
   elseif (natom<10000) then
     write(tag, '(a,I4.4,a)') ' (',index,')'
   end if
 end if

 end subroutine gettag
!!***

!!****f* ABINIT/prtnatom
!!
!! NAME
!! prtnatom
!!
!! FUNCTION
!! Print information for N atoms
!!
!!
!! INPUTS
!! prtallatoms = Logical for PRTint ALL ATOMS
!! atlist      = ATom LIST
!! index       = index for each atom
!! natom       = Number of ATOMs
!!
!! OUTPUT
!!  tag = The string to put for aech atom
!!
!! PARENTS
!!      prtxfase
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE


subroutine prtnatom(atlist,iout,message,natom,prtallatoms,thearray)

!Arguments ------------------------------------
!scalars
  logical,intent(in) :: prtallatoms
  integer,intent(in) :: natom
  logical,intent(in) :: atlist(natom)
  integer,intent(in) :: iout
  character(len=*),intent(inout) :: message
!arrays
  real(dp) :: thearray(3,natom)

!Local variables-------------------------------
!scalars
  integer :: kk
  character(len=7)   :: tag ! Maximal ' (9999)'
  character(len=18)   :: fmt

! *********************************************************************

 fmt='(a,a,1p,3e22.14,a)'

 do kk=1,natom
   if (atlist(kk)) then
     call gettag(atlist,kk,natom,prtallatoms,tag)
     write(message,fmt)TRIM(message),ch10,thearray(:,kk),tag
   end if
 end do
 !MGNAG Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output
 call wrtout(iout,message,'COLL')

 end subroutine prtnatom
!!***

end subroutine prtxfase
!!***

!!****f* ABINIT/wrt_moldyn_netcdf
!! NAME
!! wrt_moldyn_netcdf
!!
!! FUNCTION
!! Write two files for later molecular dynamics analysis:
!!  - MOLDYN.nc (netcdf format) : evolution of key quantities with time (pressure, energy, ...)
!!  - POSABIN : values of coordinates and velocities for the next time step
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (FLambert,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  itime=time step index
!!  option=1: write MOLDYN.nc file (netcdf format)
!!         2: write POSABIN file
!!         3: write both
!!  moldyn_file=name of the MD netcdf file
!!  mpi_enreg=informations about MPI parallelization
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!  rprimd(3,3)=real space primitive translations
!!  unpos=unit number for POSABIN file
!!  vel(3,natom)=velocities of atoms
!!  xred(3,natom)=reduced coordinates of atoms
!!
!! OUTPUT
!!  -- only printing --
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      metric,wrtout,xcart2xred,xred2xcart
!!
!! SOURCE

subroutine wrt_moldyn_netcdf(amass,dtset,itime,option,moldyn_file,mpi_enreg,&
&                            results_gs,rprimd,unpos,vel,xred)

 use defs_basis
 use defs_abitypes
 use m_results_gs
 use m_abicore
 use m_errors
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,   only : open_file, get_unit
 use m_geometry,   only : xcart2xred, xred2xcart, metric

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,option,unpos
 character(fnlen),intent(in) :: moldyn_file
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(results_gs_type),intent(in) :: results_gs
!arrays
 real(dp),intent(in) :: amass(dtset%natom),rprimd(3,3)
 real(dp),intent(in),target :: vel(3,dtset%natom)
 real(dp),intent(in) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer,save :: ipos=0
 integer :: iatom,ii
 character(len=500) :: message
#if defined HAVE_NETCDF
 integer :: AtomNumDimid,AtomNumId,CelId,CellVolumeId,DimCoordid,DimScalarid,DimVectorid
 integer :: EkinDimid,EkinId,EpotDimid,EpotId,EntropyDimid,EntropyId,MassDimid,MassId,NbAtomsid
 integer :: ncerr,ncid,PosId,StressDimid,StressId,TensorSymDimid
 integer :: TimeDimid,TimestepDimid,TimestepId
 logical :: atom_fix
 real(dp) :: ekin,ucvol
 character(len=fnlen) :: ficname
 character(len=16) :: chain
#endif
!arrays
#if defined HAVE_NETCDF
 integer :: PrimVectId(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable ::  xcart(:,:)
 real(dp),pointer :: vcart(:,:),vred(:,:),vtmp(:,:)
#endif

! *************************************************************************

!Only done by master processor, every nctime step
 if (mpi_enreg%me==0.and.dtset%nctime>0) then

!  Netcdf file name
#if defined HAVE_NETCDF
   ficname = trim(moldyn_file)//'.nc'
#endif

!  Xcart from Xred
#if defined HAVE_NETCDF
   ABI_ALLOCATE(xcart,(3,dtset%natom))
   call xred2xcart(dtset%natom,rprimd,xcart,xred)
#endif

!  ==========================================================================
!  First time step: write header of netcdf file
!  ==========================================================================
   if (itime==1.and.(option==1.or.option==3)) then

     ipos=0

#if defined HAVE_NETCDF
!    Write message
     write(message,'(4a)')ch10,' Open file ',trim(ficname),' to store molecular dynamics information.'
     call wrtout(std_out,message,'COLL')

!    Create netcdf file
     ncerr = nf90_create(ficname, NF90_CLOBBER , ncid)
     NCF_CHECK_MSG(ncerr,'nf90_create')

!    Dimension time for netcdf (time dim is unlimited)
     ncerr = nf90_def_dim(ncid, "time", nf90_unlimited, TimeDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Symetric Tensor Dimension
     ncerr = nf90_def_dim(ncid, "DimTensor", size(results_gs%strten), TensorSymDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Coordinates Dimension
     ncerr = nf90_def_dim(ncid, "DimCoord", size(xcart,1), DimCoordid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Atoms Dimensions
     ncerr = nf90_def_dim(ncid, "NbAtoms", dtset%natom, NbAtomsid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Vector Dimension
     ncerr = nf90_def_dim(ncid, "DimVector", 3 , DimVectorid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Scalar Dimension
     ncerr = nf90_def_dim(ncid, "DimScalar", 1 , DimScalarid)
     NCF_CHECK_MSG(ncerr,'nf90_def_dim')

!    Time step and time unit
     ncerr = nf90_def_var(ncid, "Time_step", nf90_double , DimScalarid, TimestepDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, TimestepDimid, "units", "atomic time unit")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Ionic masses
     ncerr = nf90_def_var(ncid, "Ionic_Mass", nf90_double , NbAtomsid, MassDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, MassDimid, "units", "atomic mass unit")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Ionic atomic numbers
     ncerr = nf90_def_var(ncid, "Ionic_Atomic_Number", nf90_double , NbAtomsid, AtomNumDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')

!    E_pot
     ncerr = nf90_def_var(ncid, "E_pot", nf90_double , TimeDimid, EpotDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, EpotDimid, "units", "hartree")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    E_kin
     ncerr = nf90_def_var(ncid, "E_kin", nf90_double , TimeDimid, EkinDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, EkinDimid, "units", "hartree")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Entropy
     ncerr = nf90_def_var(ncid, "Entropy", nf90_double , TimeDimid, EntropyDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, EntropyDimid, "units", "")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Stress tensor
     ncerr = nf90_def_var(ncid, "Stress", nf90_double , (/TensorSymDimid,TimeDimid/), StressDimid)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, StressDimid, "units", "hartree/bohr^3")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Positions
     ncerr = nf90_def_var(ncid, "Position", nf90_double ,(/DimCoordid,NbAtomsid,TimeDimid/), PosId)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, PosId, "units", "bohr")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    Celerities
     ncerr = nf90_def_var(ncid, "Celerity", nf90_double ,(/DimCoordid,NbAtomsid,TimeDimid/), CelId)
     NCF_CHECK_MSG(ncerr,'nf90_def_var')
     ncerr = nf90_put_att(ncid, CelId, "units", "bohr/(atomic time unit)")
     NCF_CHECK_MSG(ncerr,'nf90_put_att')

!    In case of volume cell constant
     if (dtset%optcell==0) then
!      Primitive vectors
       do ii = 1,3
         write(unit=chain,fmt='(a15,i1)') "PrimitiveVector",ii
         ncerr = nf90_def_var(ncid, trim(chain), nf90_double , DimVectorid, PrimVectId(ii))
         NCF_CHECK_MSG(ncerr,'nf90_def_var')
       end do
!      Cell Volume
       ncerr = nf90_def_var(ncid, "Cell_Volume", nf90_double , DimScalarid, CellVolumeId)
       NCF_CHECK_MSG(ncerr,'nf90_def_var')
       ncerr = nf90_put_att(ncid, CellVolumeId, "units", "bohr^3")
       NCF_CHECK_MSG(ncerr,'nf90_put_att')
     end if

!    Leave define mode and close file
     ncerr = nf90_enddef(ncid)
     NCF_CHECK_MSG(ncerr,'nf90_enddef')
     ncerr = nf90_close(ncid)
     NCF_CHECK_MSG(ncerr,'nf90_close')
#endif
   end if

!  ==========================================================================
!  Write data to netcdf file (every nctime time step)
!  ==========================================================================
   if (mod(itime, dtset%nctime)==0.and.(option==1.or.option==3)) then

     ipos=ipos+1

#if defined HAVE_NETCDF
!    Write message
     write(message,'(3a)')ch10,' Store molecular dynamics information in file ',trim(ficname)
     call wrtout(std_out,message,'COLL')

!    Open netcdf file
     ncerr = nf90_open(ficname, nf90_write, ncid)
     NCF_CHECK_MSG(ncerr,'nf90_open')

!    Time step
     ncerr = nf90_inq_varid(ncid, "Time_step", TimestepId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, TimestepId, dtset%dtion)
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Ionic masses
     ncerr = nf90_inq_varid(ncid, "Ionic_Mass", MassId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, MassId, amass, start = (/ 1 /), count=(/dtset%natom/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Ionic atomic numbers
     ncerr = nf90_inq_varid(ncid, "Ionic_Atomic_Number", AtomNumId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, AtomNumId, dtset%znucl(dtset%typat(:)),start=(/1/),count=(/dtset%natom/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Epot
     ncerr = nf90_inq_varid(ncid, "E_pot", EpotId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, EpotId, (/results_gs%etotal/), start=(/ipos/),count=(/1/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Ekin
     ekin=zero;atom_fix=(maxval(dtset%iatfix)>0)
     if (dtset%ionmov==1.or.(.not.atom_fix)) then
       vcart => vel
     else
       ABI_ALLOCATE(vcart,(3,dtset%natom))
       ABI_ALLOCATE(vred,(3,dtset%natom))
       vtmp => vel
       call xcart2xred(dtset%natom,rprimd,vtmp,vred)
       do iatom=1,dtset%natom
         do ii=1,3
           if (dtset%iatfix(ii,iatom)==1) vred(ii,iatom)=zero
         end do
       end do
       call xred2xcart(dtset%natom,rprimd,vcart,vred)
       ABI_DEALLOCATE(vred)
     end if
     do iatom=1,dtset%natom
       do ii=1,3
         ekin=ekin+half*amass(iatom)*vcart(ii,iatom)**2
       end do
     end do
     if (dtset%ionmov/=1.and.atom_fix)  then
       ABI_DEALLOCATE(vcart)
     end if
     ncerr = nf90_inq_varid(ncid, "E_kin", EkinId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, EkinId, (/ekin/), start = (/ipos/),count=(/1/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    EntropyDimid
     ncerr = nf90_inq_varid(ncid, "Entropy", EntropyId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, EntropyId, (/results_gs%energies%entropy/),start = (/ipos/),count=(/1/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Stress tensor
     ncerr = nf90_inq_varid(ncid, "Stress", StressId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, StressId, results_gs%strten, &
&     start=(/1,ipos/),count=(/size(results_gs%strten)/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Positions
     ncerr = nf90_inq_varid(ncid, "Position", PosId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, PosId, xcart, start=(/1,1,ipos/), &
&     count=(/size(xcart,1),dtset%natom,1/))
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    Celerities
     ncerr = nf90_inq_varid(ncid, "Celerity", CelId)
     NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
     ncerr = nf90_put_var(ncid, CelId, vel, start=(/1,1,ipos/), &
&     count=(/size(vel,1),dtset%natom,1/)  )
     NCF_CHECK_MSG(ncerr,'nf90_put_var')

!    In case of volume cell constant
     if (dtset%optcell==0.and.ipos==1) then
!      Primitive vectors
       do ii = 1,3
         write(unit=chain,fmt='(a15,i1)') "PrimitiveVector",ii
         ncerr = nf90_inq_varid(ncid, trim(chain), PrimVectId(ii) )
         NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
         ncerr = nf90_put_var(ncid, PrimVectId(ii), rprimd(:,ii))
         NCF_CHECK_MSG(ncerr,'nf90_put_var')
       end do
!      Cell Volume
       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
       ncerr = nf90_inq_varid(ncid, "Cell_Volume" , CellVolumeId)
       NCF_CHECK_MSG(ncerr,'nf90_inq_varid')
       ncerr = nf90_put_var(ncid, CellVolumeId, ucvol)
       NCF_CHECK_MSG(ncerr,'nf90_put_var')
     end if

!    Close file
     ncerr = nf90_close(ncid)
     NCF_CHECK_MSG(ncerr,'nf90_close')
#endif
   end if

!  ==========================================================================
!  Write data to POSABIN file (every nctime time step if option=3)
!  ==========================================================================
   if ((mod(itime, dtset%nctime)==0.and.option==3).or.(option==2)) then

!    Open file for writing
     if (open_file('POSABIN',message,unit=unpos,status='replace',form='formatted') /= 0 ) then
       MSG_ERROR(message)
     end if

!    Write Positions
     if (dtset%natom>=1) write(unpos,'(a7,3d18.5)') 'xred  ',(xred(ii,1),ii=1,3)
     if (dtset%natom>1) then
       do iatom=2,dtset%natom
         write(unpos,'(7x,3d18.5)') (xred(ii,iatom),ii=1,3)
       end do
     end if

!    Write Velocities
     if (dtset%natom>=1) write(unpos,'(a7,3d18.5)') 'vel  ',(vel(ii,1),ii=1,3)
     if (dtset%natom>1) then
       do iatom=2,dtset%natom
         write(unpos,'(7x,3d18.5)') (vel(ii,iatom),ii=1,3)
       end do
     end if

!    Close file
     close(unpos)
   end if

#if defined HAVE_NETCDF
   ABI_DEALLOCATE(xcart)
#endif

!  ==========================================================================
!  End if master proc
 end if

!Fake lines
#if !defined HAVE_NETCDF
 if (.false.) write(std_out,*) moldyn_file,results_gs%etotal,rprimd(1,1)
#endif

end subroutine wrt_moldyn_netcdf
!!***

end module m_mover
!!***
