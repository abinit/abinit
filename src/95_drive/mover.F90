!{\src2tex{textfont=tt}}
!!****f* ABINIT/mover
!! NAME
!! mover
!!
!! FUNCTION
!! Move ion or change acell acording to forces and stresses
 !!
 !! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
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
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
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
!! This subroutine uses the arguments natom, xred, vel, amass,
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
!!      abihist_free,abihist_init,abimover_fin,abimover_ini,abimover_nullify
!!      chkdilatmx,crystal_free,crystal_init,dtfil_init_time
!!      effective_potential_evaluate,erlxconv,fcart2fred,fconv,hist2var
!!      initylmg,matr3inv,monte_carlo_step,mttk_fin,mttk_ini,prec_simple
!!      pred_bfgs,pred_delocint,pred_diisrelax,pred_hmc,pred_isokinetic
!!      pred_isothermal,pred_langevin,pred_lbfgs,pred_lotf,pred_moldyn
!!      pred_nose,pred_simple,pred_srkna14,pred_steepdesc,pred_velverlet
!!      pred_verlet,prtxfase,read_md_hist,scfcv_run,status,symmetrize_xred
!!      var2hist,vel2hist,write_md_hist,wrt_moldyn_netcdf,wrtout,wvl_mkrho
!!      wvl_wfsinp_reformat,xfh_update,xmpi_barrier,xmpi_isum,xmpi_wait
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mover(scfcv_args,ab_xfh,acell,amass,dtfil,&
& electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old,&
& effective_potential,verbose,writeHIST)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abimover
 use m_abihist
 use m_xmpi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_LOTF
 use lotfpath
 use m_pred_lotf
#endif

 use m_fstrings,           only : strcat, sjoin
 use m_crystal,            only : crystal_init, crystal_free, crystal_t
 use m_crystal_io,         only : crystal_ncwrite_path
 use m_time,               only : abi_wtime, sec2str
 use m_exit,               only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_electronpositron,   only : electronpositron_type
 use m_scfcv,              only : scfcv_t, scfcv_run
 use m_effective_potential,only : effective_potential_type,effective_potential_evaluate

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mover'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_45_geomoptim
 use interfaces_56_recipspace
 use interfaces_59_ionetcdf
 use interfaces_67_common
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => mover
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
type(scfcv_t),intent(inout) :: scfcv_args
type(datafiles_type),intent(inout),target :: dtfil
type(electronpositron_type),pointer :: electronpositron
type(ab_xfh_type),intent(inout) :: ab_xfh
type(effective_potential_type),optional,intent(inout) :: effective_potential
logical,optional,intent(in) :: verbose
logical,optional,intent(in) :: writeHIST
!arrays
real(dp),intent(inout) :: acell(3)
real(dp), intent(in),target :: amass(:) !(scfcv%dtset%natom) cause segfault of g95 on yquem_g95 A check of the dim has been added
real(dp), pointer :: rhog(:,:),rhor(:,:)
real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom),xred_old(3,scfcv_args%dtset%natom)
real(dp), intent(inout) :: vel(3,scfcv_args%dtset%natom),vel_cell(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
integer,parameter :: level=102,master=0
type(abihist) :: hist,hist_prev
type(abimover) :: ab_mover
type(abimover_specs) :: specs
type(abiforstr) :: preconforstr ! Preconditioned forces and stress
type(mttk_type) :: mttk_vars
integer :: itime,icycle,itime_hist,iexit=0,ifirst,timelimit_exit,ncycle,nhisttot,kk,jj,me
integer :: nloop,ntime,option,comm
integer :: nerr_dilatmx,my_quit,ierr,quitsum_request
integer ABI_ASYNC :: quitsum_async
character(len=500) :: message
character(len=500) :: dilatmx_errmsg
character(len=8) :: stat4xml
character(len=35) :: fmt
character(len=fnlen) :: filename
real(dp) :: favg
logical :: DEBUG=.FALSE., need_verbose=.TRUE.,need_writeHIST=.TRUE.
logical :: need_scfcv_cycle = .TRUE.
logical :: change,useprtxfase
logical :: skipcycle
integer :: minIndex,ii,similar,conv_retcode
integer :: iapp
real(dp) :: minE,wtime_step,now,prev
type(crystal_t) :: crystal
!arrays
real(dp) :: gprimd(3,3),rprim(3,3),rprimd_prev(3,3)
real(dp),allocatable :: amu(:),fred_corrected(:,:),xred_prev(:,:)

! ***************************************************************
 need_verbose=.TRUE.
 if(present(verbose)) need_verbose = verbose

 need_writeHIST=.TRUE.
 if(present(writeHIST)) need_writeHIST = writeHIST
 
 call status(0,dtfil%filstat,iexit,level,'init          ')

 if ( size(amass) /= scfcv_args%dtset%natom ) then
   MSG_BUG("amass does not have the proper size.")
 end if

 ! enable time limit handler if not done in callers.
 if (need_verbose .and. enable_timelimit_in(ABI_FUNC) == ABI_FUNC) then
   write(std_out,*)"Enabling timelimit check in function: ",trim(ABI_FUNC)," with timelimit: ",trim(sec2str(get_timelimit()))
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
!09. Loop for icycle (From 1 to ncycles)
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
!IONMOV values:
!
!1.  Molecular dynamics without viscosity (vis=0)
!1.  Molecular dynamics with viscosity (vis/=0)
!2.  Broyden-Fletcher-Goldfard-Shanno method (forces)
!3.  Broyden-Fletcher-Goldfard-Shanno method (forces,Tot energy)
!4.  Conjugate gradient of potential and ionic degrees of freedom
!5.  Simple relaxation of ionic positions
!6.  Verlet algorithm for molecular dynamics
!7.  Verlet algorithm blocking every atom where dot(vel,force)<0
!8.  Verlet algorithm with a nose-hoover thermostat
!9.  Langevin molecular dynamics
!10. BFGS with delocalized internal coordinates
!11. Conjugate gradient algorithm
!12. Isokinetic ensemble molecular dynamics
!13. Isothermal/isenthalpic ensemble molecular dynamics
!14. Symplectic algorithm Runge-Kutta-Nystrom SRKNa14
!20. Ionic positions relaxation using DIIS
!21. Steepest descent algorithm
!23. LOTF method
!24. Velocity verlet molecular dynamics

 call abimover_nullify(ab_mover)

 call abimover_ini(ab_mover,specs,&
& scfcv_args%dtset%delayperm,&
& scfcv_args%dtset%diismemory,&
& scfcv_args%dtset%goprecon,&
& scfcv_args%dtset%jellslab,&
& scfcv_args%dtset%natom,&
& scfcv_args%dtset%nconeq,&
& scfcv_args%dtset%nnos,&
& scfcv_args%dtset%nsym,&
& scfcv_args%dtset%ntypat,&
& scfcv_args%dtset%optcell,&
& scfcv_args%dtset%restartxf,&
& scfcv_args%dtset%signperm,&
& scfcv_args%dtset%ionmov,&
& scfcv_args%dtset%bmass,&
& scfcv_args%dtset%dtion,&
& scfcv_args%dtset%friction,&
& scfcv_args%dtset%mdwall,&
& scfcv_args%dtset%noseinert,&
& scfcv_args%dtset%strprecon,&
& scfcv_args%dtset%vis,&
& scfcv_args%dtset%iatfix,&
& scfcv_args%dtset%symrel,&
& scfcv_args%dtset%typat,&
& scfcv_args%dtset%prtatlist,&
& amass,&
& scfcv_args%dtset%goprecprm,&
& scfcv_args%dtset%mdtemp,&
& scfcv_args%dtset%strtarget,&
& scfcv_args%dtset%qmass,&
& scfcv_args%dtset%znucl,&
& dtfil%fnameabi_hes,&
& dtfil%filnam_ds)

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

 nhisttot=ncycle*ntime;if (scfcv_args%dtset%nctime>0) nhisttot=nhisttot+1

!AM_2017 New version of the hist, we just store the needed history step not all of them...
 if(specs%nhist/=-1) nhisttot = specs%nhist! We don't need to store all the history
 call abihist_init(hist,ab_mover%natom,nhisttot,specs%isVused,specs%isARused)
 call abiforstr_ini(preconforstr,ab_mover%natom)

!###########################################################
!### 06. First output before any itime or icycle

!If effective potential is present,
!  forces will be compute with it
 if (present(effective_potential)) then
   need_scfcv_cycle = .FALSE.
   if(need_verbose)then
     write(message,'(2a,i2,5a,80a)')&
&   ch10,'=== [ionmov=',ab_mover%ionmov,'] ',trim(specs%method),' with effective potential',&
&   ch10,('=',kk=1,80)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
 else
   if(need_verbose)then
     write(message,'(a,a,i2,a,a,a,80a)')&
&     ch10,'=== [ionmov=',ab_mover%ionmov,'] ',specs%method,&
&     ch10,('=',kk=1,80)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
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

!Compute atomic masses (in a.u.)
 ABI_ALLOCATE(amu,(scfcv_args%dtset%ntypat))
 do kk=1,ab_mover%ntypat
   do jj=1,ab_mover%natom
     if (kk==ab_mover%typat(jj)) then
       amu(kk)=amass(jj)/amu_emass
       exit
     end if
   end do
 end do

!Fill history with the values of xred, acell and rprimd
 call var2hist(acell,hist,ab_mover%natom,rprimd,xred,DEBUG)

!Fill velocities and ionic kinetic energy
 call vel2hist(ab_mover%amass,hist,vel,vel_cell)
 hist%time(hist%ihist)=zero

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

   ! Handle time limit condition.
   if (itime == 1) prev = abi_wtime()
   if (itime  > 1) then
     now = abi_wtime()
     wtime_step = now - prev
     prev = now
     write(message,*)sjoin("mover: previous time step took ",sec2str(wtime_step))
     if(need_verbose)call wrtout(std_out, message, "COLL")
     if (have_timelimit_in(ABI_FUNC)) then
       if (itime > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then 
           write(message,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),&
&           ". Will exit itime loop in mover."
           if(need_verbose)MSG_COMMENT(message)
           if(need_verbose)call wrtout(ab_out, message, "COLL")
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
!  ### 09. Loop for icycle (From 1 to ncycles)
   do icycle=1,ncycle
     itime_hist = (itime-1)*ncycle + icycle ! Store the time step in of the history
     
!    ###########################################################
!    ### 10. Output for each icycle (and itime)
     if(need_verbose)then
       write(message,fmt)&
&     ch10,'--- Iteration: (',itime,'/',ntime,') Internal Cycle: (',icycle,'/',ncycle,')',ch10,('-',kk=1,80)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
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
&       ch10,('-',kk=1,3),&
&       'SELF-CONSISTENT-FIELD CONVERGENCE',('-',kk=1,44)
       else
         write(message,'(a,3a,33a,44a)')&
&       ch10,('-',kk=1,3),&
&       'EFFECTIVE POTENTIAL CALCULATION',('-',kk=1,44)
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
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
           call effective_potential_evaluate( &
&           effective_potential,scfcv_args%results_gs%etotal,&
&           scfcv_args%results_gs%fcart,scfcv_args%results_gs%fred,&
&           scfcv_args%results_gs%strten,ab_mover%natom,rprimd,xred=xred,verbose=need_verbose)

!          Check if the simulation does not diverged...
           if(itime > 10 .and.ABS(scfcv_args%results_gs%etotal - hist%etot(1)) > 1E4)then
!            We set to false the flag corresponding to the bound 
             effective_potential%anharmonics_terms%bounded = .FALSE.
             if(need_verbose.and.me==master)then
               message = "The simulation is diverging, please check your effective potential"
               MSG_WARNING(message)
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
!      * Inside scfcv.F90 there is a call to symmetrize_xred.F90
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
           call xfh_update(ab_xfh,acell,fred_corrected,ab_mover%natom,rprim,&
&           hist%strten(:,hist%ihist),xred)
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
       call write_md_hist(hist,filename,ifirst,itime_hist,ab_mover%natom,ab_mover%ntypat,&
&       ab_mover%typat,amu,ab_mover%znucl,ab_mover%dtion,scfcv_args%dtset%mdtemp)
     end if
#endif

!    ###########################################################
!    ### 14. Output after SCFCV
     if(need_verbose.and.need_scfcv_cycle)then
       write(message,'(a,3a,a,72a)')&
&       ch10,('-',kk=1,3),'OUTPUT',('-',kk=1,71)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
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
           call erlxconv(hist%etot(:),iexit,hist%ihist,&
&           itime,hist%mxhist,ntime,scfcv_args%dtset%tolmxde)
         end if
       end if
     end if

     if (itime==ntime.and.icycle==ncycle) iexit=1

!    ###########################################################
!    ### 16. => Precondition forces, stress and energy

     write(message,*) 'Geometry Optimization Precondition:',ab_mover%goprecon
     if(need_verbose)call wrtout(std_out,message,'COLL')
     if (ab_mover%goprecon>0)then
       call prec_simple(ab_mover,preconforstr,hist,icycle,itime,0)
     end if

!    ###########################################################
!    ### 17. => Call to each predictor

!    MT->GAF: dirty trick to predict vel(t)
!    do a double loop: 1- compute vel, 2- exit
     nloop=1


     if (scfcv_args%dtset%nctime>0.and.iexit==1) then
       iexit=0;nloop=2
     end if

     do ii=1,nloop
       if (ii==2) iexit=1

       select case (ab_mover%ionmov)
       case (1)
         call pred_moldyn(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit)
       case (2,3)
         call pred_bfgs(ab_mover,ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
       case (4,5)
         call pred_simple(ab_mover,hist,iexit)
       case (6,7)
         call pred_verlet(ab_mover,hist,ab_mover%ionmov,itime,ntime,DEBUG,iexit)
       case (8)
         call pred_nose(ab_mover,hist,itime,ntime,DEBUG,iexit)
       case (9)
         call pred_langevin(ab_mover,hist,icycle,itime,ncycle,ntime,DEBUG,iexit,skipcycle)
       case (10,11)
         call pred_delocint(ab_mover,ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
       case (12)
         call pred_isokinetic(ab_mover,hist,itime,ntime,DEBUG,iexit)
       case (13)
         call pred_isothermal(ab_mover,hist,itime,mttk_vars,ntime,DEBUG,iexit)
       case (14)
         call pred_srkna14(ab_mover,hist,icycle,DEBUG,iexit,skipcycle)
       case (20)
         call pred_diisrelax(ab_mover,hist,itime,ntime,DEBUG,iexit)
       case (21)
         call pred_steepdesc(ab_mover,preconforstr,hist,itime,DEBUG,iexit)
       case (22)
         call pred_lbfgs(ab_mover,ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
#if defined HAVE_LOTF
       case (23)
         call pred_lotf(ab_mover,hist,itime,icycle,DEBUG,iexit)
#endif
       case (24)
         call pred_velverlet(ab_mover,hist,itime,ntime,DEBUG,iexit)
       case (25)
         call pred_hmc(ab_mover,hist,itime,icycle,ntime,ncycle,DEBUG,iexit)

       case default
         write(message,"(a,i0)") "Wrong value of ionmov: ",ab_mover%ionmov
         MSG_ERROR(message)
       end select

     end do

     ! check dilatmx here and correct if necessary
     if (scfcv_args%dtset%usewvl == 0) then
       call chkdilatmx(scfcv_args%dtset%dilatmx,rprimd,scfcv_args%dtset%rprimd_orig(1:3,1:3,1),&
&       dilatmx_errmsg)
       _IBM6("dilatxm_errmsg: "//TRIM(dilatmx_errmsg))
       if (LEN_TRIM(dilatmx_errmsg) /= 0) then
         MSG_WARNING(dilatmx_errmsg)
         nerr_dilatmx = nerr_dilatmx+1
         if (nerr_dilatmx > 3) then
           ! Write last structure before aborting, so that we can restart from it.
           ! zion is not available, but it's not useful here.
           if (me == master) then
             ! Init crystal
             call crystal_init(scfcv_args%dtset%amu_orig(:,1),crystal,0,ab_mover%natom,&
&             scfcv_args%dtset%npsp,ab_mover%ntypat,scfcv_args%dtset%nsym,rprimd,ab_mover%typat,xred,&
&             [(-one, ii=1,ab_mover%ntypat)],ab_mover%znucl,2,.False.,.False.,"dilatmx_structure",&
&             symrel=scfcv_args%dtset%symrel,tnons=scfcv_args%dtset%tnons,symafm=scfcv_args%dtset%symafm)

#ifdef HAVE_NETCDF
             ! Write netcdf file
             filename = strcat(dtfil%filnam_ds(4), "_DILATMX_STRUCT.nc")
             NCF_CHECK(crystal_ncwrite_path(crystal, filename))
#endif
             call crystal_free(crystal)
           end if
           call xmpi_barrier(comm)
           write (dilatmx_errmsg, '(a,i0,3a)') &
&           'Dilatmx has been exceeded too many times (', nerr_dilatmx, ')',ch10, &
&           'Restart your calculation from larger lattice vectors and/or a larger dilatmx'
           MSG_ERROR_CLASS(dilatmx_errmsg, "DilatmxError")
         end if
       end if
     end if

!    Write MOLDYN netcdf and POSABIN files (done every dtset%nctime time step)
     if(ab_mover%ionmov/=23 .or. icycle==1)then
       if (scfcv_args%dtset%nctime>0) then
         jj=itime; if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1) jj=jj-hist_prev%mxhist
         if (jj>0) then
           option=3
           call wrt_moldyn_netcdf(amass,scfcv_args%dtset,jj,option,dtfil%fnameabo_moldyn,&
&           scfcv_args%mpi_enreg,scfcv_args%results_gs,&
&           hist%rprimd(:,:,hist%ihist-1),dtfil%unpos,hist%vel(:,:,hist%ihist),&
&           hist%xred(:,:,hist%ihist-1))
         end if
         if (iexit==1) hist%ihist=hist%ihist-1
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
     if (ab_mover%ionmov==13) then
       if (hist%ihist>2) then
         ii = abihist_findIndex(hist,-2)
         vel_cell(:,:)=(hist%rprimd(:,:,hist%ihist)- hist%rprimd(:,:,ii)/(two*ab_mover%dtion))
       else if (hist%ihist>1) then
         ii = abihist_findIndex(hist,-1)
         vel_cell(:,:)=(hist%rprimd(:,:,hist%ihist)-hist%rprimd(:,:,ii))/(ab_mover%dtion)
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

     write(message,*) 'ICYCLE',icycle,skipcycle
     if(need_verbose)call wrtout(std_out,message,'COLL')
     write(message,*) 'NCYCLE',ncycle
     if(need_verbose)call wrtout(std_out,message,'COLL')
     if (skipcycle) exit

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
   if ((ab_mover%ionmov/=4.and.ab_mover%ionmov/=5)) then
     if (scfcv_args%dtset%tolmxf/=0)then
       call fconv(hist%fcart(:,:,hist%ihist-1),&
&       scfcv_args%dtset%iatfix, &
&       iexit, itime,&
&       ab_mover%natom,&
&       ntime,&
&       ab_mover%optcell,&
&       scfcv_args%dtset%strfact,&
&       scfcv_args%dtset%strtarget,&
&       hist%strten(:,hist%ihist-1),&
&       scfcv_args%dtset%tolmxf)
     else
       call erlxconv(hist%etot(:),iexit,hist%ihist,&
&       itime,hist%mxhist,ntime,scfcv_args%dtset%tolmxde)
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

 if (ab_mover%goprecon>0)then
   call prec_simple(ab_mover,preconforstr,hist,1,1,1)
 end if

 if (ab_mover%ionmov==13)then
   call mttk_fin(mttk_vars)
 end if

 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(xred_prev)

 call abihist_free(hist)
 call abihist_free(hist_prev)

 call abimover_fin(ab_mover)
 call abiforstr_fin(preconforstr)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine mover
!!***
