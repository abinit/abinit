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
!!
!! NOTES
!! This subroutine uses the arguments natom, xred, vel, fcart, amass,
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
!!      abiforstr_fin,abiforstr_ini,abihist_bcast,abihist_compare,abihist_fin
!!      abihist_ini,abimover_fin,abimover_ini,abimover_nullify,chkdilatmx
!!      crystal_free,crystal_init,dtfil_init_time,effective_potential_evaluate
!!      erlxconv,fcart2fred,fconv,fred2fcart,hist2var,initylmg,metric,mkradim
!!      monte_carlo_step,mttk_fin,mttk_ini,prec_simple,pred_bfgs,pred_delocint
!!      pred_diisrelax,pred_isokinetic,pred_isothermal,pred_langevin,pred_lbfgs
!!      pred_lotf,pred_moldyn,pred_nose,pred_simple,pred_srkna14,pred_steepdesc
!!      pred_verlet,prtxfase,read_md_hist,scfcv_run,status,symmetrize_xred
!!      var2hist,vel2hist,write_md_hist,wrt_moldyn_netcdf,wrtout,wvl_mkrho
!!      wvl_wfsinp_reformat,xfh_update,xmpi_barrier,xmpi_isum,xmpi_wait
!!      xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mover(scfcv_args,ab_xfh,acell,amass,dtfil,&
& electronpositron,rhog,rhor,rprimd,vel,vel_cell,xred,xred_old,effective_potential)

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

 use m_fstrings,         only : strcat, sjoin
 use m_crystal,          only : crystal_init, crystal_free, crystal_t
 use m_crystal_io,       only : crystal_ncwrite_path
 use m_time,             only : abi_wtime, sec2str
 use m_exit,             only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_electronpositron, only : electronpositron_type
 use m_scfcv,            only : scfcv_t, scfcv_run
 use m_effective_potential
 use m_monte_carlo

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
!type(dataset_type),intent(inout),target :: dtset
type(electronpositron_type),pointer :: electronpositron
!type(paw_dmft_type) :: paw_dmft
!type(ab_scfcv_args_in),intent(inout) :: ab_scfcv_in
!type(ab_scfcv_args_inout),intent(inout) :: ab_scfcv_inout
type(ab_xfh_type),intent(inout) :: ab_xfh
!arrays
!no_abirules
real(dp),intent(inout) :: acell(3)
real(dp), intent(in),target :: amass(:) !(scfcv%dtset%natom) cause segfault of g95 on yquem_g95 A check of the dim has been added
real(dp), pointer :: rhog(:,:),rhor(:,:)
real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom),xred_old(3,scfcv_args%dtset%natom)
real(dp), intent(inout) :: vel(3,scfcv_args%dtset%natom),vel_cell(3,3),rprimd(3,3)
type(effective_potential_type),optional,intent(in):: effective_potential
!Local variables-------------------------------
!scalars
integer,parameter :: level=102,master=0
type(abihist) :: hist,hist_prev
type(abimover) :: ab_mover
type(abimover_specs) :: specs
type(abiforstr) :: preconforstr ! Preconditioned forces and stress
type(mttk_type) :: mttk_vars
integer :: itime,icycle,iexit=0,timelimit_exit,ncycle,kk,jj,me,nloop,ntime,option,comm
integer :: nerr_dilatmx,my_quit,ierr,quitsum_request
integer ABI_ASYNC :: quitsum_async
character(len=500) :: message
character(len=500) :: dilatmx_errmsg
character(len=8) :: stat4xml
character(len=35) :: fmt
character(len=fnlen) :: filename
real(dp) :: ucvol,favg
logical :: DEBUG=.FALSE.
logical :: need_scfcv_cycle = .TRUE.
logical :: change,useprtxfase
logical :: skipcycle
integer :: minIndex,ii,similar,conv_retcode
integer :: iapp
real(dp) :: minE,tolerance,wtime_step,now,prev
type(crystal_t) :: crystal
!arrays
real(dp) :: xred_tmp(3,scfcv_args%dtset%natom)
real(dp) :: xcart(3,scfcv_args%dtset%natom)
real(dp) :: fcart2(3,scfcv_args%dtset%natom)
real(dp) :: fred2(3,scfcv_args%dtset%natom)
real(dp) :: favg2(3)
real(dp) :: fred_corrected(3,scfcv_args%dtset%natom)
real(dp) :: rprim(3,3)
real(dp) :: gprimd(3,3)
real(dp) :: gmet(3,3)
real(dp) :: rmet(3,3)

! ***************************************************************

 call status(0,dtfil%filstat,iexit,level,'init          ')

 if ( size(amass) /= scfcv_args%dtset%natom ) then
   MSG_BUG("amass does not have the proper size.")
 end if

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(ABI_FUNC) == ABI_FUNC) then
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
!07. Compute xcart and fill the history of the first SCFCV
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
!21. Set the final values of xcart and xred
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


 if (ab_mover%ionmov==31.and..not.present(effective_potential)) then
   MSG_BUG("effective_potential is not present.")
 end if

 if (ab_mover%ionmov==13)then
   call mttk_ini(mttk_vars,ab_mover%nnos)
 end if

!write(std_out,*) 'mover 03'
!###########################################################
!### 03. Set the number of iterations ntime
!###     By default ntime==1 but if the user enter a lower
!###     value mover will execute at least one iteration

 if (scfcv_args%dtset%ntime<1)then
   ntime=1
 else
   ntime=scfcv_args%dtset%ntime
 end if

!write(std_out,*) 'mover 04'
!###########################################################
!### 04. Try to read history of previous calculations
!###     It requires access to the NetCDF library

!Init MPI data
 comm=scfcv_args%mpi_enreg%comm_cell
 me=xmpi_comm_rank(comm)

#if defined HAVE_NETCDF
 filename=trim(ab_mover%filnam_ds(4))//'_HIST.nc'

 if (ab_mover%restartxf<=0)then
!  Read history from file (and broadcast if MPI)
   if (me==master) then
     call read_md_hist(filename,hist_prev)
   end if
   if (xmpi_comm_size(comm) > 1) then
     call abihist_bcast(hist_prev,master,comm)
   end if

!  If restartxf specifies to reconstruct the history
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-1)then
     ntime=ntime+hist_prev%mxhist
   end if

!  If restartxf specifies to start from the lowest energy
   if (hist_prev%mxhist>0.and.ab_mover%restartxf==-2)then
     minE=hist_prev%histE(1)
     minIndex=1
     do ii=1,hist_prev%mxhist
       write(std_out,*) 'Iteration:',ii,' Total Energy:',hist_prev%histE(ii)
       if (minE>hist_prev%histE(ii))then
         minE=hist_prev%histE(ii)
         minIndex=ii
       end if
     end do
     write(std_out,*) 'The lowest energy occurs at iteration:',minIndex,'etotal=',minE
     acell(:)   =hist_prev%histA(:,minIndex)
     xcart(:,:) =hist_prev%histXF(:,:,1,minIndex)
     xred(:,:)  =hist_prev%histXF(:,:,2,minIndex)
     rprimd(:,:)=hist_prev%histR(:,:,minIndex)
   end if

 end if !if (ab_mover%restartxf<=0)

#endif

!write(std_out,*) 'mover 05'
!###########################################################
!### 05. Allocate the hist structure
 iexit=0; timelimit_exit=0
 hist%isVused=specs%isVused
 hist%isARused=specs%isARused
 ncycle=specs%ncycle

 if (scfcv_args%dtset%nctime>0) then
   call abihist_ini(hist,ab_mover%natom,ncycle*ntime+1)
 else
   call abihist_ini(hist,ab_mover%natom,ncycle*ntime)
 end if

 call abiforstr_ini(preconforstr,ab_mover%natom)

!write(std_out,*) 'mover 06'
!###########################################################
!### 06. First output before any itime or icycle

 ! if effective potential is present,
 ! forces will be compute with it
 if (present(effective_potential)) then
   need_scfcv_cycle = .FALSE.
   write(message,'(a,a,i2,a,a,a,a,80a)')&
&   ch10,'=== [ionmov=',ab_mover%ionmov,'] ',trim(specs%method),&
&   ' with effective potential',ch10,('=',kk=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 else
   write(message,'(a,a,i2,a,a,a,80a)')&
&   ch10,'=== [ionmov=',ab_mover%ionmov,'] ',specs%method,&
&   ch10,('=',kk=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Format for printing on each cycle
 write(fmt,'(a6,i2,a4,i2,a4,i2,a4,i2,a9)')&
& '(a,a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ntime))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,i',int(log10(real(ncycle))+1),&
& ',a,a,80a)'
!write(std_out,*) fmt
!write(std_out,*) HUGE(1)
!write(std_out,*) int(log10(real(HUGE(1)))+1)

!write(std_out,*) 'mover 07'
!###########################################################


!Compute xcart from xred, and rprimd
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

!Compute rprim from rprimd and acell
 do kk=1,3
   do jj=1,3
     rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
   end do
 end do

!Fill history with the values of xred,xcart,acell and rprimd
 call var2hist(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)

!Fill velocities and ionic kinetic energy
 call vel2hist(ab_mover%amass,hist,ab_mover%natom,vel)
 hist%histT(hist%ihist)=0.0

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

!write(std_out,*) 'mover 08'
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
     call wrtout(std_out, message, "COLL")
     if (have_timelimit_in(ABI_FUNC)) then
       if (itime > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then 
           write(message,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),&
&           ". Will exit itime loop in mover."
           MSG_COMMENT(message)
           call wrtout(ab_out, message, "COLL")
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
!  write(std_out,*) 'mover 09'
!  ###########################################################
!  ### 09. Loop for icycle (From 1 to ncycles)
   do icycle=1,ncycle

!    write(std_out,*) 'mover 10'
!    ###########################################################
!    ### 10. Output for each icycle (and itime)

     write(message,fmt)&
&     ch10,'--- Iteration: (',itime,'/',ntime,') Internal Cycle: (',icycle,'/',ncycle,')',ch10,('-',kk=1,80)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     if (useprtxfase) then
       call prtxfase(ab_mover,hist,std_out,mover_BEFORE)
     end if

!    write(std_out,*) 'mover 11'
!    ###########################################################
!    ### 11. Symmetrize atomic coordinates over space group elements
     change=.FALSE.
     xred_tmp(:,:)=xred(:,:)

     call symmetrize_xred(scfcv_args%indsym,ab_mover%natom,&
&     scfcv_args%dtset%nsym,scfcv_args%dtset%symrel,scfcv_args%dtset%tnons,xred)
     do kk=1,ab_mover%natom
       do jj=1,3
         if (xred(jj,kk)/=xred_tmp(jj,kk)) change=.TRUE.
       end do
     end do
     
     if (change)then
       call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
       hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
       hist%histXF(:,:,2,hist%ihist)=xred(:,:)
       write(std_out,*) 'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED'
       write(std_out,*) 'DIFFERENCES:'

       do kk=1,ab_mover%natom
         write(std_out,*) xred(:,kk)-xred_tmp(:,kk)
       end do
       
       xred_tmp(:,:)=xred(:,:)
     end if

!    write(std_out,*) 'mover 12'
!    ###########################################################
!    ### 12. => Call to SCFCV routine and fill history with forces
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

     if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then
       
       call abihist_compare(hist_prev,hist,ab_mover%natom,similar,tolerance)
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
         if (need_scfcv_cycle.and.scfcv_args%dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
         !  WVL - Before running scfcv, on non-first geometry step iterations,
         !  we need to reformat the wavefunctions, taking into acount the new
         !  coordinates.
         !  We prepare to change rhog (to be removed) and rhor.
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
         if (need_scfcv_cycle) then

           call dtfil_init_time(dtfil,iapp)
           call scfcv_run(scfcv_args,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
           
           if (conv_retcode == -1) then
             message = "Scf cycle returned conv_retcode == -1 (timelimit is approaching), this should not happen inside mover"
             MSG_WARNING(message)
           end if
!        For monte carlo don't need ton recompute energy here
!        every is done in pred_montecarlo
         else if(ab_mover%ionmov /= 31) then
           call effective_potential_evaluate(effective_potential,&
&           scfcv_args%results_gs%etotal,&
&           scfcv_args%results_gs%fcart,scfcv_args%results_gs%fred,&
&           scfcv_args%results_gs%strten,ab_mover%natom,rprimd,&
&           xcart)
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
!      The solution here is to recompute acell and rprim
!      and store those values in the present record
!      even if initially those values were not exactly
!      the values entering in scfcv
!      
!      One test case with these condition is bigdft/t10
       if (need_scfcv_cycle.and.scfcv_args%dtset%usewvl == 1) then
         call mkradim(acell,rprim,rprimd)
         hist%histA(:,hist%ihist)=acell(:)
         hist%histR(:,:,hist%ihist)=rprimd(:,:)
       end if
       
!      ANOMALOUS SITUATIONS
!      * In ionmov 4 & 5 xred could change inside SCFCV
!      So we need to take the values from the output
!      
!      * Inside scfcv.F90 there is a call to symmetrize_xred.F90
!      for the first SCF cycle symmetrize_xred could change xred
!      so we need always take xred convert to xcart and
!      store in the history

       if (ab_mover%ionmov<10)then
         change=.FALSE.
         
         do kk=1,ab_mover%natom
           do jj=1,3
             if (xred(jj,kk)/=xred_tmp(jj,kk)) change=.TRUE.
           end do
         end do

         
         if (change)then
           call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
           hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
           hist%histXF(:,:,2,hist%ihist)=xred(:,:)
           call wrtout(std_out,'WARNING: ATOMIC COORDINATES WERE SYMMETRIZED AFTER SCFCV','COLL') 
           call wrtout(std_out,'DIFFERENCES:','COLL') 
           
           do kk=1,ab_mover%natom
             write(std_out,*) xred(:,kk)-xred_tmp(:,kk)
           end do
           
         end if !if (change)

!        call xred2xcart(ab_mover%natom,rprimd,xcart,xred)
!        hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
!        hist%histXF(:,:,2,hist%ihist)=xred(:,:)
       end if

!        Fill velocities and ionic kinetic energy
       call vel2hist(ab_mover%amass,hist,ab_mover%natom,vel)
       hist%histE(hist%ihist)       =scfcv_args%results_gs%etotal
       hist%histEnt(hist%ihist)     =scfcv_args%results_gs%energies%entropy
       hist%histXF(:,:,3,hist%ihist)=scfcv_args%results_gs%fcart(:,:)
       hist%histXF(:,:,4,hist%ihist)=scfcv_args%results_gs%fred(:,:)
       hist%histS(:,hist%ihist)     =scfcv_args%results_gs%strten(:)
       hist%histT(hist%ihist)       =itime

!      !######################################################################
!      ! Test of convertion fcart and fred

       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
       
       call fred2fcart(favg2,fcart2,scfcv_args%results_gs%fred,gprimd,ab_mover%jellslab,ab_mover%natom)
       call fcart2fred(scfcv_args%results_gs%fcart,fred2,rprimd,ab_mover%natom)

!      write (std_out,*) 'FCART'
!      do kk=1,ab_mover%natom
!      write (std_out,*) scfcv_args%results_gs%fcart(:,kk)
!      end do
!      write (std_out,*) 'FCART converted from fred'
!      do kk=1,ab_mover%natom
!      write (std_out,*) fcart2(:,kk)
!      end do
!      write (std_out,*) 'FRED'
!      do kk=1,ab_mover%natom
!      write (std_out,*) scfcv_args%results_gs%fred(:,kk)
!      end do
!      write (std_out,*) 'FRED converted from fcart'
!      do kk=1,ab_mover%natom
!      write (std_out,*) fred2(:,kk)
!      end do

!      !######################################################################

     end if ! if(hist_prev%mxhist>0.and.ab_mover%restartxf==-1.and.hist_prev%ihist<=hist_prev%mxhist)then

     if((ab_xfh%nxfh==0.or.itime/=1)) then
       !Not yet needing for effective potential
       call mkradim(acell,rprim,rprimd)
!      Get rid of mean force on whole unit cell, but only if no
!      generalized constraints are in effect
!      hist%histXF(:,:,4,hist%ihist) are reduced forces
       if(ab_mover%nconeq==0)then
         do kk=1,3
           favg=sum(hist%histXF(kk,:,4,hist%ihist))/dble(ab_mover%natom)
           fred_corrected(kk,:)=hist%histXF(kk,:,4,hist%ihist)-favg
           if(ab_mover%jellslab/=0.and.kk==3)&
&           fred_corrected(kk,:)=hist%histXF(kk,:,4,hist%ihist)
         end do
       else
         fred_corrected(:,:)=hist%histXF(:,:,4,hist%ihist)
       end if

       if (ncycle<10.and.ab_mover%restartxf>=0) then
         call xfh_update(ab_xfh,acell,fred_corrected,ab_mover%natom,rprim,&
&         hist%histS(:,hist%ihist),xred)
       end if

     end if

!    write(std_out,*) 'mover 13'
!    ###########################################################
!    ### 13. Write the history into the _HIST file
!    ###

#if defined HAVE_NETCDF
     if (me==master) then
       call write_md_hist(hist,ab_mover,filename,icycle,itime)
     end if
#endif

!    write(std_out,*) 'mover 14'
!    ###########################################################
!    ### 14. Output after SCFCV
     write(message,'(a,3a,a,72a)')&
&     ch10,('-',kk=1,3),'OUTPUT',('-',kk=1,71)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')

     if (useprtxfase) then
       call prtxfase(ab_mover,hist,ab_out,mover_AFTER)
       call prtxfase(ab_mover,hist,std_out,mover_AFTER)
     end if

!    !    DEBUG (XRA AFTER SCFCV)
!    if(DEBUG)then
!    write (std_out,*) '---XRA AFTER SCFCV---'
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

!    write(std_out,*) 'mover 15'
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
           call fconv(hist%histXF(:,:,3,hist%ihist),&
&           scfcv_args%dtset%iatfix, &
&           iexit, itime,&
&           ab_mover%natom,&
&           ntime,&
&           ab_mover%optcell,&
&           scfcv_args%dtset%strfact,&
&           scfcv_args%dtset%strtarget,&
&           hist%histS(:,hist%ihist),&
&           scfcv_args%dtset%tolmxf)
         else
           call erlxconv(hist%histE(:),iexit,hist%ihist,&
&           itime,hist%mxhist,ntime,scfcv_args%dtset%tolmxde)
         end if
       end if
     end if

     if (itime==ntime.and.icycle==ncycle) iexit=1

!    write(std_out,*) 'mover 16'
!    ###########################################################
!    ### 16. => Precondition forces, stress and energy

     write(message,*) 'Geometry Optimization Precondition:',ab_mover%goprecon
     call wrtout(std_out,message,'COLL')
     if (ab_mover%goprecon>0)then
       call prec_simple(ab_mover,preconforstr,hist,icycle,itime,0)
     end if

!    write(std_out,*) 'mover 16'
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

       case (31)         
         call monte_carlo_step(ab_mover,effective_potential,hist,itime,ntime,DEBUG,iexit)
         write(std_out,*) "Developpement Monte carlo"
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
&           hist%histR(:,:,hist%ihist-1),dtfil%unpos,hist%histV(:,:,hist%ihist),&
&           hist%histXF(:,:,1,hist%ihist-1),hist%histXF(:,:,2,hist%ihist-1))
         end if
         if (iexit==1) hist%ihist=hist%ihist-1
       end if
     end if
     if(iexit/=0) exit

!    write(std_out,*) 'mover 17'
!    ###########################################################
!    ### 18. Use the history  to extract the new values
!    ###     acell, rprimd, xcart and xred

!    !    DEBUG (XRA BEFORE USE OF HISTORY)
!    if(DEBUG)then
!    write (std_out,*) '---XRA BEFORE USE OF HISTORY---',hist%ihist
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

     call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)

     if(ab_mover%optcell/=0)then

       call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

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

     vel(:,:)=hist%histV(:,:,hist%ihist)

!    vell_cell(3,3)= velocities of cell parameters
!    Not yet used here but compute it for consistency
     vel_cell(:,:)=zero
     if (ab_mover%ionmov==13) then
       if (hist%ihist>2) then
         vel_cell(:,:)=(hist%histR(:,:,hist%ihist)-hist%histR(:,:,hist%ihist-2))/(two*ab_mover%dtion)
       else if (hist%ihist>1) then
         vel_cell(:,:)=(hist%histR(:,:,hist%ihist)-hist%histR(:,:,hist%ihist-1))/(ab_mover%dtion)
       end if
     end if

!    !    DEBUG (XRA AFTER USE OF HISTORY)
!    if(DEBUG)then
!    write (std_out,*) '---XRA AFTER USE OF HISTORY---',hist%ihist
!    write (std_out,*) 'XCART'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xcart(:,kk)
!    end do
!    write (std_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (std_out,*) xred(:,kk)
!    end do
!    if (ab_mover%ionmov==1)then
!    write (std_out,*) 'VEL'
!    do kk=1,ab_mover%natom
!    write (std_out,*) hist%histV(:,kk,hist%ihist)
!    end do
!    end if
!    write(std_out,*) 'RPRIMD'
!    do kk=1,3
!    write(std_out,*) rprimd(:,kk)
!    end do
!    write(std_out,*) 'ACELL'
!    write(std_out,*) acell(:)
!    end if

!    This is needed for some compilers such as
!    pathscale, g95, xlf that do not exit
!    from a loop if you change the upper limit
!    inside
     if (icycle>=ncycle .and. scfcv_args%mpi_enreg%me == 0) then
       write(std_out,*) 'EXIT:',icycle,ncycle
       exit
     end if

     write(message,*) 'ICYCLE',icycle,skipcycle
     call wrtout(std_out,message,'COLL')
     write(message,*) 'NCYCLE',ncycle
     call wrtout(std_out,message,'COLL')
     if (skipcycle) exit

!    write(std_out,*) 'mover 18'
!    ###########################################################
!    ### 19. End loop icycle
   end do ! do icycle=1,ncycle

   if(iexit/=0)exit
!  write(std_out,*) 'mover 19'
!  ###########################################################
!  ### 20. End loop itime
 end do ! do itime=1,ntime


 ! Call fconv here if we exited due to wall time limit.
 if (timelimit_exit==1 .and. specs%isFconv) then
   iexit = timelimit_exit
   ntime = itime-1
   if ((ab_mover%ionmov/=4.and.ab_mover%ionmov/=5)) then
     !write(std_out,*)"hist%ihist",hist%ihist
     if (scfcv_args%dtset%tolmxf/=0)then
       call fconv(hist%histXF(:,:,3,hist%ihist-1),&
&       scfcv_args%dtset%iatfix, &
&       iexit, itime,&
&       ab_mover%natom,&
&       ntime,&
&       ab_mover%optcell,&
&       scfcv_args%dtset%strfact,&
&       scfcv_args%dtset%strtarget,&
&       hist%histS(:,hist%ihist-1),&
&       scfcv_args%dtset%tolmxf)
     else
       call erlxconv(hist%histE(:),iexit,hist%ihist,&
&       itime,hist%mxhist,ntime,scfcv_args%dtset%tolmxde)
     end if
   end if
 end if

 ! Avoid pending requests if itime == ntime.
 call xmpi_wait(quitsum_request,ierr)

!###########################################################
!### 21. Set the final values of xcart and xred with the last
!###     computed values (not the last predicted)

 call hist2var(acell,hist,ab_mover%natom,rprim,rprimd,xcart,xred,DEBUG)
 vel(:,:)=hist%histV(:,:,hist%ihist)

 if (DEBUG .and. ab_mover%ionmov==1)then
   write (std_out,*) 'vel'
   do kk=1,ab_mover%natom
     write (std_out,*) hist%histV(:,kk,hist%ihist)
   end do
 end if

!write(std_out,*) 'mover 21'
!###########################################################
!### 22. XML Output at the end

!XML output of the status
 if (scfcv_args%mpi_enreg%me == 0 .and. scfcv_args%dtset%prtxml == 1) then
   write(ab_xml_out, "(3a)") '    <geometryMinimisation type="',trim(specs%type4xml),'">'
   write(ab_xml_out, "(5a)") '      <status cvState="',trim(stat4xml) ,'" stop-criterion="',trim(specs%crit4xml),'" />'
   write(ab_xml_out, "(3a)") '    </geometryMinimisation>'
 end if

!write(std_out,*) 'mover 22'
!###########################################################
!### 23. Deallocate hist and ab_mover datatypes

 if (ab_mover%goprecon>0)then
   call prec_simple(ab_mover,preconforstr,hist,1,1,1)
 end if

 if (ab_mover%ionmov==13)then
   call mttk_fin(mttk_vars)
 end if

 call abihist_fin(hist)
 call abihist_fin(hist_prev)
 call abimover_fin(ab_mover)
 call abiforstr_fin(preconforstr)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

end subroutine mover
!!***
