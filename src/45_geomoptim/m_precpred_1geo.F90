!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_precpred_1geo
!! NAME
!!  m_precpred_1geo
!!
!! FUNCTION
!! Single geometry : apply force and stress preconditioner followed by geometry predictor. 
!! Choose among the whole set of geometry predictors defined by iomov.
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (DCA, XG, GMR, SE)
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

module m_precpred_1geo

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

 use m_fstrings,           only : strcat, sjoin, indent
 use m_symtk,              only : matr3inv, symmetrize_xred
 use m_geometry,           only : fcart2fred, chkdilatmx, xred2xcart
 use m_crystal,            only : crystal_init, crystal_free, crystal_t
 use m_crystal_io,         only : crystal_ncwrite_path
 use m_time,               only : abi_wtime, sec2str
 use m_exit,               only : get_start_time
 use m_scfcv,              only : scfcv_t, scfcv_run
 use m_dtfil,              only : dtfil_init_time, status
 use m_initylmg,           only : initylmg
 use m_xfpack,             only : xfh_update
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

 implicit none

 private
!!***

 public :: precpred_1geo
!!***

contains
!!***

!!****f* ABINIT/precpred_1geo
!! NAME
!! mover
!!
!! FUNCTION
!! Single geometry : apply force and stress preconditioner followed by geometry predictor.
!! Choose among the whole set of geometry predictors defined by iomov.
!!
!! INPUTS
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
!!  cg(2,mcg)=array for planewave coefficients of wavefunctions.
!!  initialized= if 0 the initialisation of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  occ(mband*nkpt*nsppol=occupation number for each band (usually 2) at each k point.
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  xred(3,natom)=reduced dimensionless atomic coordinates; updated on output
!!  verbose = optional, default is true, flag to disable the verbose mode
!!  write_HIST = optional, default is true, flag to disble the write of the HIST file
!!
!! NOTES
!! This subroutine uses the arguments natom, xred, 
!! vis, and dtion (the last two contained in dtset) to make
!! molecular dynamics updates.  
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine precpred_1geo(scfcv_args,ab_xfh,dtfil,rprimd,xred,verbose)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mover'
 use interfaces_14_hidewrite
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
type(scfcv_t),intent(inout) :: scfcv_args
type(datafiles_type),intent(inout),target :: dtfil
type(ab_xfh_type),intent(inout) :: ab_xfh
logical,optional,intent(in) :: verbose
!arrays
real(dp), intent(inout) :: xred(3,scfcv_args%dtset%natom)
real(dp), intent(inout) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
integer,parameter :: master=0
type(abihist) :: hist
type(abimover) :: ab_mover
type(abiforstr) :: preconforstr ! Preconditioned forces and stress
type(mttk_type) :: mttk_vars
integer :: itime,icycle,iexit,ncycle,me
integer :: nloop,ntime,comm
integer :: nerr_dilatmx
character(len=500) :: message
character(len=500) :: dilatmx_errmsg
character(len=fnlen) :: filename
logical :: DEBUG=.FALSE., need_verbose=.TRUE.
logical :: skipcycle
integer :: ii
type(crystal_t) :: crystal
!arrays
real(dp) :: rprim(3,3)
real(dp),allocatable :: amu(:)
! ***************************************************************
 need_verbose=.TRUE.
 if(present(verbose)) need_verbose = verbose

 comm=scfcv_args%mpi_enreg%comm_cell
 me=xmpi_comm_rank(comm)

!Precondition forces, stress and energy
 call abiforstr_ini(preconforstr,ab_mover%natom)
 write(message,'(2a)') 'Geometry Optimization Precondition:',ab_mover%goprecon
 if(need_verbose)call wrtout(std_out,message,'COLL')
 if (ab_mover%goprecon>0)then
   call prec_simple(ab_mover,preconforstr,hist,icycle,itime,0)
 end if

!Call to each predictor
!MT->GAF: dirty trick to predict vel(t)
!do a double loop: 1- compute vel, 2- exit
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
   case (15)
     call pred_fire(ab_mover, ab_xfh,preconforstr,hist,ab_mover%ionmov,itime,DEBUG,iexit)
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
     call pred_hmc(ab_mover,hist,itime,icycle,ntime,scfcv_args%dtset%hmctt,DEBUG,iexit)
   case (27)
     !In case of ionmov 27, all the atomic configurations have been computed at the
     !begining of the routine in generate_training_set, thus we just need to increase the indexes
     !in the hist
     hist%ihist = abihist_findIndex(hist,+1)

   case default
     write(message,"(a,i0)") "Wrong value of ionmov: ",ab_mover%ionmov
     MSG_ERROR(message)
   end select

 end do

 ! check dilatmx here and correct if necessary
 if (scfcv_args%dtset%usewvl == 0) then
   call chkdilatmx(scfcv_args%dtset%chkdilatmx,scfcv_args%dtset%dilatmx,&
    rprimd,scfcv_args%dtset%rprimd_orig(1:3,1:3,1),dilatmx_errmsg)
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
&         scfcv_args%dtset%npsp,ab_mover%ntypat,scfcv_args%dtset%nsym,rprimd,ab_mover%typat,xred,&
&         [(-one, ii=1,ab_mover%ntypat)],ab_mover%znucl,2,.False.,.False.,"dilatmx_structure",&
&         symrel=scfcv_args%dtset%symrel,tnons=scfcv_args%dtset%tnons,symafm=scfcv_args%dtset%symafm)

#ifdef HAVE_NETCDF
         ! Write netcdf file
         filename = strcat(dtfil%filnam_ds(4), "_DILATMX_STRUCT.nc")
         NCF_CHECK(crystal_ncwrite_path(crystal, filename))
#endif
         call crystal_free(crystal)
       end if
       call xmpi_barrier(comm)
       write (dilatmx_errmsg, '(a,i0,3a)') &
        'Dilatmx has been exceeded too many times (', nerr_dilatmx, ')',ch10, &
        'Restart your calculation from larger lattice vectors and/or a larger dilatmx'
       MSG_ERROR_CLASS(dilatmx_errmsg, "DilatmxError")
     end if
   end if
 end if

 call abiforstr_fin(preconforstr)

end subroutine precpred_1geo
!!***

end module m_precpred_1geo
!!***
