!!****m* ABINIT/m_dfpt_cgwf
!! NAME
!!  m_dfpt_cgwf
!!
!! FUNCTION
!! Update one single wavefunction (cwavef), non self-consistently.
!! Uses a conjugate-gradient algorithm.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG,DRH,XW,FJ,MT,LB)
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

module m_dfpt_cgwf

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_cgtools
 use m_rf2

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_set_zero, pawcprj_axpby
 use m_hamiltonian, only : gs_hamiltonian_type, rf_hamiltonian_type, KPRIME_H_KPRIME
 use m_getghc,      only : getghc
 use m_getgh1c,     only : getgh1c, getdc1

 implicit none

 private
!!***

 public :: dfpt_cgwf
!!***

contains
!!***

!!****f* ABINIT/dfpt_cgwf
!! NAME
!! dfpt_cgwf
!!
!! FUNCTION
!! Update one single wavefunction (cwavef), non self-consistently.
!! Uses a conjugate-gradient algorithm.
!! Try to keep close to the formulas in PRB55, 10337 (1997) [[cite:Gonze1997]], for the
!! non-self-consistent case, except that we are computing here
!! the second-derivative of the total energy, and not E(2). There
!! is a factor of 2 between the two quantities.
!! The wavefunction that is generated is always orthogonal to cgq.
!! It is orthogonal to the active Hilbert space, and will be complemented
!! by contributions from the active space in the calling routine, if needed.
!!
!! INPUTS
!!  band=which particular band we are converging.
!!  band_me=cpu-local index of band which we are converging.
!!  berryopt=option for Berry phase
!!  cgq(2,mcgq)=wavefunction coefficients for MY bands at k+Q
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwaveprj0(natom,nspinor*usecprj)=GS wave function at k projected with nl projectors
!!  band_procs(nband)=tags for processors which have the other bands for cgq below
!!  eig0_k=0-order eigenvalues for the present wavefunction at k
!!  eig0_kq(nband)=GS eigenvalues at k+Q (hartree)
!!  grad_berry(2,mpw1,dtefield%mband_occ) = the gradient of the Berry phase term
!!  gscq(2,mgscq)=<g|S|Cnk+q> coefficients for MY bands (PAW) at k+Q
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+Q
!!  icgq=shift to be applied on the location of data in the array cgq
!!  igscq=shift to be applied on the location of data in the array gscq
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  mcgq=second dimension of the cgq array
!!  mgscq=second dimension of gscq, with only mband_mem bands
!!  mpi_enreg=information about MPI parallelization
!!  mpw1=maximum number of planewave for first-order wavefunctions
!!  natom=number of atoms in cell.
!!  nband=number of bands.
!!  nband_me=number of bands on this cpu
!!  nbdbuf=number of buffer bands for the minimisation
!!  nline=number of line minimizations per band.
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+Q
!!  nspinor=number of spinorial components of the wavefunctions
!!  opt_gvnlx1=option controlling the use of gvnlx1 array:
!!            0: used as an output
!!            1: used as an input:    - used only for ipert=natom+2
!!                 NCPP: contains the ddk 1-st order WF
!!                 PAW: contains frozen part of 1st-order hamiltonian
!!            2: used as input/ouput:    - used only for PAW and ipert=natom+2
!!                 At input: contains the ddk 1-st order WF (times i)
!!                 At output: contains frozen part of 1st-order hamiltonian
!!  prtvol=control print volume and debugging output
!!  quit= if 1, proceeds to smooth ending of the job.
!!  dfpt_sciss=scissor shift (Ha)
!!  tolrde=tolerance on the ratio of differences of energies (for the line minimisation)
!!  tolwfr=tolerance on largest wf residual
!!  usedcwavef=flag controlling the use of dcwavef array (PAW only):
!!             0: not used (not allocated)
!!             1: used as input
!!             2: used as output
!!  wfoptalg=govern the choice of algorithm for wf optimisation (0 or 10, at present)
!!
!! OUTPUT
!!  eig1_k(2*nband**2)=matrix of first-order eigenvalues (hartree)
!!                     eig1(:,ii,jj)=<C0 ii|H1|C0 jj> for norm-conserving psps
!!                     eig1(:,ii,jj)=<C0 ii|H1-(eig0_k+eig0_k+q)/2.S(1)|C0 jj> for PAW
!!  ghc(2,npw1*nspinor)=<G|H0-eig0_k.I|C1 band,k> (NCPP) or <G|H0-eig0_k.S0|C1 band,k> (PAW)
!!  gvnlxc(2,npw1*nspinor)=<G|Vnl+VFockACE|C1 band,k>
!!  gvnlx1(2,npw1*nspinor)=  part of <G|K1+Vnl1|C0 band,k> not depending on VHxc1           (NCPP)
!!                       or part of <G|K1+Vnl1-eig0k.S1|C0 band,k> not depending on VHxc1 (PAW)
!!  resid=wf residual for current band
!!  gh1c_n= <G|H1|C0 band,k> (NCPP) or <G|H1-eig0k.S1|C0 band,k> (PAW).
!!          This vector is not projected on the subspace orthogonal to the WF.
!!  === if gs_hamkq%usepaw==1 ===
!!  gsc(2,npw1*nspinor*usepaw)=<G|S0|C1 band,k>
!!
!! SIDE EFFECTS
!!  Input/Output:
!!  cwavef(2,npw1*nspinor)=first-order  wavefunction at k,q, in reciprocal space (updated)
!!
!!  === if gs_hamkq%usepaw==1 ===
!!  cwaveprj(natom,nspinor)= wave functions at k projected with nl projectors
!!
!!  === if also usedcwavef>0 ===
!!  dcwavef(2,npw1*nspinor)=change of wavefunction due to change of overlap:
!!         dcwavef is delta_Psi(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!!         see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42)
!!         input if usedcwavef=1, output if usedcwavef=2
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_precon,cg_zaxpy,cg_zcopy,dotprod_g,getdc1,getgh1c,getghc
!!      pawcprj_alloc,pawcprj_axpby,pawcprj_free,pawcprj_set_zero,projbd
!!      sqnorm_g,timab,wrtout
!!
!! SOURCE

subroutine dfpt_cgwf(band,band_me,band_procs,berryopt,cgq,cwavef,cwave0,cwaveprj,cwaveprj0,&
& rf2,dcwavef,&
& eig0_k,eig0_kq,eig1_k,ghc,gh1c_n,grad_berry,gsc,gscq,&
& gs_hamkq,gvnlxc,gvnlx1,icgq,idir,ipert,igscq,&
& mcgq,mgscq,mpi_enreg,mpw1,natom,nband,nband_me,nbdbuf,nline_in,npw,npw1,nspinor,&
& opt_gvnlx1,prtvol,quit,resid,rf_hamkq,dfpt_sciss,tolrde,tolwfr,&
& usedcwavef,wfoptalg,nlines_done)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,berryopt
 integer,intent(in) :: band_me, nband_me
 integer,intent(in) :: icgq,idir,igscq,ipert,mcgq,mgscq,mpw1,natom,nband
 integer,intent(in) :: nbdbuf,nline_in,npw,npw1,nspinor,opt_gvnlx1
 integer,intent(in) :: prtvol,quit,usedcwavef,wfoptalg
 integer,intent(inout) :: nlines_done
 real(dp),intent(in) :: dfpt_sciss,tolrde,tolwfr
 real(dp),intent(out) :: resid
 type(MPI_type),intent(in) :: mpi_enreg
 type(rf2_t), intent(in) :: rf2
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
!arrays
 integer,intent(in) :: band_procs(nband)
 real(dp),intent(in) :: cgq(2,mcgq),eig0_kq(nband)
 real(dp),intent(in) :: eig0_k(nband)
 real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband),gscq(2,mgscq)
 real(dp),intent(inout) :: cwave0(2,npw*nspinor),cwavef(2,npw1*nspinor)
 real(dp),intent(inout) :: dcwavef(2,npw1*nspinor*((usedcwavef+1)/2))
 real(dp),intent(inout) :: eig1_k(2*nband**2)
 real(dp),intent(out) :: gh1c_n(2,npw1*nspinor)
 real(dp),intent(out) :: ghc(2,npw1*nspinor)
 real(dp),intent(out) :: gsc(2,npw1*nspinor*gs_hamkq%usepaw)
 real(dp),intent(inout) :: gvnlx1(2,npw1*nspinor),gvnlxc(2,npw1*nspinor)
 type(pawcprj_type),intent(inout) :: cwaveprj(natom,nspinor)
 type(pawcprj_type),intent(inout) :: cwaveprj0(natom,nspinor*gs_hamkq%usecprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=15,tim_getgh1c=1,tim_getghc=2,tim_projbd=2
 integer,save :: nskip=0
 integer :: cpopt,iband,igs,iline,indx_cgq,ipw,me_g0,comm_fft
 integer :: iband_me, jband_me, ierr, me_band, band_off, unit_me
 integer :: ipws,ispinor,istwf_k,jband,nline,optlocal,optnl,dc_shift_band,sij_opt
 integer :: test_is_ok,useoverlap,usepaw,usevnl,usetolrde
 real(dp) :: d2edt2,d2te,d2teold,dedt,deltae,deold,dotgg
 real(dp) :: dotgp,doti,dotr,eshift,eshiftkq,gamma,optekin,prod1,prod2
 real(dp) :: theta,tol_restart,u1h0me0u1
 logical :: gen_eigenpb
 integer :: skipme, bands_skipped_now(nband)
 character(len=500) :: msg
!arrays
 integer :: bands_treated_now (nband)
 real(dp) :: dummy(0,0),tsec(2)
 real(dp) :: eig1_k_loc(2,nband,nband)
 real(dp),allocatable :: conjgr(:,:),cwaveq(:,:),cwwork(:,:),direc(:,:)
 real(dp),allocatable :: gberry(:,:),gh1c(:,:),gh_direc(:,:),gresid(:,:),gvnlx1_saved(:,:)
 real(dp),allocatable :: gs1c(:,:),gvnlx_direc(:,:),pcon(:),sconjgr(:,:)
 real(dp),allocatable :: scprod(:,:),work(:,:),work1(:,:),work2(:,:)
 real(dp),pointer :: kinpw1(:)
 type(pawcprj_type),allocatable :: conjgrprj(:,:)
 type(pawcprj_type) :: cprj_dummy(0,0)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(122,1,tsec)

 !======================================================================
 !========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
 !====================================================================

 nline = nline_in
 usetolrde = 1

 ! LB-23/04/17:
 ! For ipert=natom+10 or ipert=natom+11, the Sternheimer equation is non-self-consistent, so we have
 ! to solve a true linear problem (A.X = B) for each kpoint and band. In this case, the conjugate
 ! gradient algorithm can find the exact solution (within the numerical precision) with only ONE call
 ! of dfpt_cgwf (per kpoint and band...). The solution is found with at most N iterations, N being the dimension of X.
 ! In order to avoid useless scfcv loops (and many calls of rf2_init, which can be time consuming),
 ! we want to leave this routine only if 'tolwfr' is reached. Consequently, 'tolrde' is not used and 'nline' is set to 100.
 ! One could use nline=npw1*nspinor (>> 100 !) instead, but when the method cannot converge (i.e when tolwfr is lower than the
 ! numerical noise) the program could be stuck here for a very long time.
 ! NOTE : This is also true for ipert==natom+1, but a lot of references in the test suite have to be changed...
 if(ipert==natom+10.or.ipert==natom+11) then
   nline = 100 ! The default value is only 4... This should be sufficient to converge with nstep=1 or 2
   if (nline_in > 100) nline = nline_in ! Keep the possibility to increase nline
   usetolrde = 0  ! see below
 end if

 if (prtvol>=10) then
   !Tell us what is going on:
   write(msg,'(a,i6,2x,a,i3,a)')' --- dfpt_cgwf is called for band',band,'for',nline,' lines'
   call wrtout(std_out,msg)
 end if

 me_g0 = mpi_enreg%me_g0
 comm_fft = mpi_enreg%comm_fft

 me_band = mpi_enreg%me_band
 !unit_me = 300+band
 unit_me = 6
 bands_skipped_now = 0

 bands_treated_now = 0
 bands_treated_now(band) = 1
 call xmpi_sum(bands_treated_now,mpi_enreg%comm_band,ierr)
#ifdef DEV_MJV
print *, 'bands_treated_now ', bands_treated_now
print *, 'bands_skipped_now ', bands_skipped_now
#endif


 skipme = 0

 ! if PAW, one has to solve a generalized eigenproblem
 usepaw=gs_hamkq%usepaw
 gen_eigenpb=(usepaw==1)
 useoverlap=0;if (gen_eigenpb) useoverlap=1

 ! Use scissor shift on 0-order eigenvalue
 eshift=eig0_k(band)-dfpt_sciss

 ! Additional initializations
 istwf_k=gs_hamkq%istwf_k
 optekin=0;if (wfoptalg>=10) optekin=1
 tol_restart=tol12;if (gen_eigenpb) tol_restart=tol8
 if (ipert == natom+10 .or. ipert == natom+11) tol_restart = tol7

 kinpw1 => gs_hamkq%kinpw_kp

 ! Memory allocations
 ABI_ALLOCATE(gh1c,(2,npw1*nspinor))
 ABI_ALLOCATE(pcon,(npw1))
 ABI_ALLOCATE(scprod,(2,nband_me))

 if (berryopt== 4.or.berryopt== 6.or.berryopt== 7.or. berryopt==14.or.berryopt==16.or.berryopt==17) then
   ABI_ALLOCATE(gberry,(2,npw1*nspinor))
   gberry(:,1:npw1*nspinor)=grad_berry(:,1:npw1*nspinor,band)
#ifdef DEV_MJV
print *, 'gberry 294 ', gberry(:,1:10)
#endif
 else
   ABI_ALLOCATE(gberry,(0,0))
 end if

!TODO MJV: this should probably be adjusted as well for natom+10/11 perts: set to band_me instead of band
 dc_shift_band=(band-1)*npw1*nspinor

! this is used many times - no use de and re allocating
 ABI_ALLOCATE(work,(2,npw1*nspinor))

#ifdef DEV_MJV
print *, 'bands_treated_now 308 ', bands_treated_now
print *, 'bands_skipped_now 308 ', bands_skipped_now
#endif
!DEBUG!! Several checking statements
 if (prtvol==-level.or.prtvol==-19.or.prtvol==-20) then
   write(msg,'(a)') " ** cgwf3 : debugging mode, tests will be done"
   ! Search CGWF3_WARNING in the log file to find errors (if any)
   call wrtout(std_out,msg)
   ABI_ALLOCATE(work1,(2,npw1*nspinor))
   !  ===== Check <Psi_k+q^(0)|S(0)|Psi_k+q^(0)>=delta_{ij}
   if (.not.gen_eigenpb) work1(:,:)=cgq(:,1+npw1*nspinor*(band_me-1)+icgq:npw1*nspinor*band_me+icgq)
   if (     gen_eigenpb) work1(:,:)=gscq(:,1+npw1*nspinor*(band_me-1)+igscq:npw1*nspinor*band_me+igscq)
! NB: this loop is not band-block diagonal
!  the present logic does a lot of communication: 1 for each band.
!  Could be grouped outside the jband loop into 1 big one, but if the wf are big this is a waste. Tradeoffs...
   jband_me = 0
   do jband=1,nband
     if (bands_treated_now(jband)-bands_skipped_now(jband) == 0) cycle
     if (band_procs(jband) == me_band) then
       jband_me = jband_me + 1
       work(:,:)=cgq(:,1+npw1*nspinor*(jband_me-1)+icgq:npw1*nspinor*jband_me+icgq)
     end if
! send to everyone else, who is also working on jband right now
     call xmpi_bcast(work,band_procs(jband),mpi_enreg%comm_band,ierr)  

     call dotprod_g(dotr,doti,istwf_k,npw1*nspinor,2,work1,work,me_g0,mpi_enreg%comm_spinorfft)
     test_is_ok=1
     if(jband==band) then
       if(abs(dotr-one)>tol12) test_is_ok=0
     else
       if(abs(dotr)>tol12) test_is_ok=0
     end if
     if(abs(doti)>tol12) test_is_ok=0
     if(test_is_ok/=1) then
       write(msg,'(a,i3,a,2es22.15)') "CGWF3_WARNING : <Psi_k+q,i^(0)|S(0)|Psi_k+q,j^(0)> for band j=",jband," is ",dotr,doti
       call wrtout(std_out,msg)
     end if
   end do

   !  ===== Check Pc.Psi_k+q^(0)=0
   ! each jband is checked by everybody, against the bands attributed to present cpu
   ! NB - this does not depend on the "band" input
   jband_me = 0
   do jband=1,nband
     if (bands_treated_now(jband)-bands_skipped_now(jband) == 0) cycle
     if (band_procs(jband) == me_band) then
       jband_me = jband_me + 1
       work(:,:)=cgq(:,1+npw1*nspinor*(jband_me-1)+icgq:npw1*nspinor*jband_me+icgq)
     end if
! send to everyone else, who is also working on jband right now
     call xmpi_bcast(work,band_procs(jband),mpi_enreg%comm_band,ierr)
     work1 = work

     call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
       gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)

!
! if bands are parallelized, I have only projected against bands on my cpu
!   Pc|work>  = |work> - Sum_l <psi_{k+q, l}|work> |psi_{k+q, l}>
!             = Sum_nproc_band (|work> - Sum_{my l} <psi_{k+q, l}|work> |psi_{k+q, l}>) - (nproc_band-1) |work>
!
     if (mpi_enreg%nproc_band > 1) then
       call xmpi_sum(work,mpi_enreg%comm_band,ierr)
!TODO: make this a blas call? zaxpy
       work = work - (mpi_enreg%nproc_band-1)*work1
     end if

     call sqnorm_g(dotr,istwf_k,npw1*nspinor,work,me_g0,comm_fft)
     if(sqrt(dotr)>tol12) then
       write(msg,'(a,i3,a,es22.15)') "CGWF3_WARNING : Norm of Pc.Psi_k+q_j^(0) for band j=",jband," is ",sqrt(dotr)
       call wrtout(std_out,msg)
     end if
   end do

   ! ===== Check Pc.Psi_k^(0)=0
   ! NB: this _does_ depend on the input band "band" stored in cwave0
   do iband = 1, nband
     if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
     !if (band_procs(iband) == me_band) then ! these 2 conditions should be the same
     if (iband == band) then
       work(:,:)=cwave0(:,:)
     end if
! send to everyone else, who is also working on jband right now
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
     work1 = work
     call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
       gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)

     if (mpi_enreg%nproc_band > 1) then
       call xmpi_sum(work,mpi_enreg%comm_band,ierr)
!TODO: make this a blas call? zaxpy
       work = work - (mpi_enreg%nproc_band-1)*work1
     end if

     call sqnorm_g(dotr,istwf_k,npw1*nspinor,work,me_g0,comm_fft)
     if(sqrt(dotr)>tol12) then
       write(msg,'(a,i3,a,es22.15)') "CGWF3_WARNING : Norm of Pc.Psi_k^(0) for band ",band," is ",sqrt(dotr)
       call wrtout(std_out,msg)
     end if
   end do

 
   ! ===== Check Pc.dcwavef=0 (for 2nd order only)
   if(ipert==natom+10.or.ipert==natom+11) then
     do iband = 1, nband
       if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
       !if (band_procs(iband) == me_band) then ! these 2 conditions should be the same
       if (iband == band) then
         work(:,:)=rf2%dcwavef(:,1+dc_shift_band:npw1*nspinor+dc_shift_band)
       end if
! send to everyone else, who is also working on jband right now
       call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
       work1 = work
       call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
         gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)

       if (mpi_enreg%nproc_band > 1) then
         call xmpi_sum(work,mpi_enreg%comm_band,ierr)
!TODO: make this a blas call? zaxpy
         work = work - (mpi_enreg%nproc_band-1)*work1
       end if

       call sqnorm_g(dotr,istwf_k,npw1*nspinor,work,me_g0,comm_fft)
       if(sqrt(dotr)>tol10) then
         write(msg,'(a,i3,a,es22.15)') "CGWF3_WARNING : Norm of Pc.dcwavef for band ",band," is ",sqrt(dotr)
         call wrtout(std_out,msg)
       end if
     end do
   end if

   ! ===== Check Pc^*.S(0).Psi_k+q^(0)=0
   ! NB: here again does not depend on input "band"
   if (gen_eigenpb) then
     jband_me = 0
     do jband=1,nband
       if (bands_treated_now(jband)-bands_skipped_now(jband) == 0) cycle
       if (band_procs(jband) == me_band) then
         jband_me = jband_me + 1
         work(:,:)=gscq(:,1+npw1*nspinor*(jband_me-1)+igscq:npw1*nspinor*jband_me+igscq)
       end if
! send to everyone else, who is also working on jband right now
       call xmpi_bcast(work,band_procs(jband),mpi_enreg%comm_band,ierr)
       work1 = work

       call projbd(gscq,work,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
         cgq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)

       if (mpi_enreg%nproc_band > 1) then
         call xmpi_sum(work,mpi_enreg%comm_band,ierr)
!TODO: make this a blas call? zaxpy
         work = work - (mpi_enreg%nproc_band-1)*work1
       end if

       call sqnorm_g(dotr,istwf_k,npw1*nspinor,work,me_g0,comm_fft)
       if(sqrt(dotr)>tol12) then
         write(msg,'(a,i3,a,es22.15)') "CGWF3_WARNING : Norm of Pc^*.S(0).Psi_k+q_j^(0) for band j=",jband," is ",sqrt(dotr)
         call wrtout(std_out,msg)
       end if
     end do
   end if
   ABI_DEALLOCATE(work1)
 end if
!ENDDEBUG!! Several checking statements


 !======================================================================
 !========== INITIALISATION OF MINIMIZATION ITERATIONS =================
 !======================================================================

 if(ipert/=natom+10.and.ipert/=natom+11) then
   !  The following is needed for first order perturbations only
   !  Otherwise, the work is already done in rf2_init (called in dfpt_vtowfk.F90)

   ! Compute H(1) applied to GS wavefunction Psi(0)
   if (gen_eigenpb) then
     sij_opt=1
     ABI_ALLOCATE(gs1c,(2,npw1*nspinor))
   else
     ABI_ALLOCATE(gs1c,(0,0))
     sij_opt=0
   end if
   usevnl=1; optlocal=1; optnl=2
   if (prtvol==-level.or.prtvol==-19) then
     ABI_ALLOCATE(gvnlx1_saved,(2,npw1*nspinor))
     gvnlx1_saved(:,:) = gvnlx1(:,:)
   end if
   call getgh1c(berryopt,cwave0,cwaveprj0,gh1c,gberry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,&
     mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

#ifdef DEV_MJV
print *, 'gh1c 490 ', gh1c(:,1:10)
#endif
   if (gen_eigenpb) then
     if (ipert/=natom+2) then  ! S^(1) is zero for ipert=natom+2
!$OMP PARALLEL
!$OMP DO
       do ipw=1,npw1*nspinor
         gh1c (1:2,ipw)=gh1c (1:2,ipw)-eshift*gs1c(1:2,ipw)
       end do
!$OMP END DO NOWAIT
       if (opt_gvnlx1/=1) then
!$OMP DO
         do ipw=1,npw1*nspinor
           gvnlx1(1:2,ipw)=gvnlx1(1:2,ipw)-eshift*gs1c(1:2,ipw)
         end do
!$OMP END DO NOWAIT
       end if
!$OMP END PARALLEL
     end if

     ! If generalized eigenPb and dcwavef requested, compute it:
     ! dcwavef is delta_Psi(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
     ! see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42)
     if (usedcwavef==2) then
       call getdc1(band,band_procs,bands_treated_now,cgq,cprj_dummy,dcwavef,cprj_dummy,&
&           0,icgq,istwf_k,mcgq,0,&
&           mpi_enreg,natom,nband,nband_me,npw1,nspinor,0,gs1c)
#ifdef DEV_MJV
print *, 'dcwavef 490 ', dcwavef(:,1:10)
#endif
     end if
   end if ! gen_eigenpb

 else
   ! 2nd order case (wrt k perturbation)
   ! Copy RHS_Stern(:,:) of the given band in gh1c
   gh1c(:,:)=rf2%RHS_Stern(:,1+dc_shift_band:npw1*nspinor+dc_shift_band)
 end if

 if (prtvol==-level.and.usedcwavef==2) then
   !Check that Pc^*.(H^(0)-E.S^(0)).delta_Psi^(1) is zero ! This is a consequence of P_c delta_Psi^(1) = 0
   ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
   do iband = 1, nband
     if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
     if (iband == band) then
       cwwork=dcwavef
       !  - Apply H^(0)-E.S^(0) to delta_Psi^(1)
       sij_opt=0;if (gen_eigenpb) sij_opt=-1
       cpopt=-1
       ABI_ALLOCATE(work1,(2,npw1*nspinor*((sij_opt+1)/2)))
       ABI_ALLOCATE(work2,(2,npw1*nspinor))
       call getghc(cpopt,cwwork,conjgrprj,work,work1,gs_hamkq,work2,eshift,mpi_enreg,&
         1,prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
       ABI_DEALLOCATE(work1)
       ABI_DEALLOCATE(work2)
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
     cwwork=work

   ! -Apply Pc^*
     call projbd(gscq,cwwork,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
       cgq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
  
     call xmpi_sum(cwwork,mpi_enreg%comm_band,ierr)

     if (mpi_enreg%nproc_band > 1) then
!TODO: make this a blas call? zaxpy
       cwwork = cwwork - (mpi_enreg%nproc_band-1)*work
     end if

     call sqnorm_g(dotr,istwf_k,npw1*nspinor,cwwork,me_g0,comm_fft)
     if(sqrt(dotr)>tol12) then
       write(msg,'(a,i3,a,es22.15)') 'CGWF3_WARNING : |Pc^*.(H^(0)-E.S^(0)).delta_Psi^(1)| (band ',band,')=',sqrt(dotr)
       call wrtout(std_out,msg)
     end if
   end do
   ABI_DEALLOCATE(cwwork)
 end if ! prtvol==-level.and.usedcwavef==2

 call cg_zcopy(npw1*nspinor,gh1c,gh1c_n)

#ifdef DEV_MJV
print *, 'gh1c 578 ', gh1c(:,1:10)
#endif
 ! Projecting out all bands
 ! While we could avoid calculating all the eig1_k to obtain the perturbed density,
 ! we do need all of the matrix elements when outputing the full 1st-order wfn.
 ! Note the subtlety:
 ! -For the generalized eigenPb, S|cgq> is used in place of |cgq>,
 ! in order to apply P_c+ projector (see PRB 73, 235101 (2006) [[cite:Audouze2006]], Eq. (71), (72))
 eig1_k_loc = zero
 do iband = 1, nband
#ifdef DEV_MJV
print *, 'iband, bands_treated_now(iband), bands_skipped_now(iband) ', iband, bands_treated_now(iband), bands_skipped_now(iband)
#endif
   if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
   work = zero
   if (iband == band) then
     work = gh1c
   end if
#ifdef DEV_MJV
print *, ' band_procs(iband), mpi_enreg%comm_band,mpi_enreg%nproc_band ', &
           band_procs(iband), mpi_enreg%comm_band, mpi_enreg%nproc_band
#endif
   call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)

   if(gen_eigenpb)then
     call projbd(gscq,work,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
       cgq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
   else
     call projbd(cgq,work,-1,icgq,0,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
       dummy,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
   end if

#ifdef DEV_MJV
print *, 'work1 611 ', work(:,1:10)
#endif
! sum projections against all bands k+q
   call xmpi_sum(work,mpi_enreg%comm_band,ierr)
#ifdef DEV_MJV
print *, 'work2 616 ', work(:,1:10)
#endif

   ! scprod now contains scalar products of band iband (runs over all bands in current queue) with local bands j
   jband_me = 0
   do jband=1,nband
     if (band_procs(jband) /= me_band) cycle
     jband_me = jband_me + 1
     eig1_k_loc(:,jband,iband)=scprod(:,jband_me)
   end do

! save this for me only
   if (iband == band) then
!TODO: make this a blas call? zaxpy
     gh1c = work - (mpi_enreg%nproc_band-1)*gh1c
   end if

#ifdef DEV_MJV
print *, 'eig1_k_loc 619 line for iband ', iband, eig1_k_loc(:,:,iband)
#endif
 end do !iband

#ifdef DEV_MJV
print *, 'gh1c 619 ', gh1c(:,1:10)
print *, 'gscq 619 ', gscq(:,1:10)
#endif

 if(ipert/=natom+10.and.ipert/=natom+11) then
   ! For ipert=natom+10 or natom+11, this is done in rf2_init

   ! The array eig1_k contains:
   ! <u_(jband,k+q)^(0)|H_(k+q,k)^(1)|u_(iband,k)^(0)>                              (NC psps)
   ! or <u_(jband,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(iband,k)^(0)> (PAW)
   ! so in case of PAW need to add the overlap term below
   !
   ! NB: 2019 11 15: MJV: I swapped the names of jband and iband to be more consistent with other loops above
   if (gen_eigenpb) then
     do iband=1,nband
       if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
       work = zero
       if (iband==band) then
         work = gs1c
       end if
     ! for iband on this proc, bcast to all others to get full line of iband,jband pairs
       call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr) 

     ! add PAW overlap correction term to present iband (all procs) and local jband elements
       indx_cgq=icgq
       do jband=1,nband
         if (band_procs(jband) /= me_band) cycle
     
         eshiftkq=half*(eig0_kq(jband)-eig0_k(iband))
         call dotprod_g(dotr,doti,istwf_k,npw1*nspinor,2,cgq(:,indx_cgq+1:indx_cgq+npw1*nspinor),work,&
           me_g0,mpi_enreg%comm_spinorfft)
         eig1_k_loc(1,jband,iband)=eig1_k_loc(1,jband,iband)-eshiftkq*dotr
         eig1_k_loc(2,jband,iband)=eig1_k_loc(2,jband,iband)-eshiftkq*doti
         indx_cgq=indx_cgq+npw1*nspinor
       end do
#ifdef DEV_MJV
print *, 'me ', me_band, ' eig1_k_loc 672 line for iband ', iband, eig1_k_loc(:,:,iband)
#endif
     end do ! iband
   end if ! PAW and generalized eigenproblem

   ! No more need of gs1c
   ABI_DEALLOCATE(gs1c)

   ! reduce over band procs to fill in the matrix for all jband (distributed over procs)
   ! must only do this once for eig1_k_loc: now have all jband for current ibands on all procs
   call xmpi_sum(eig1_k_loc,mpi_enreg%comm_band,ierr)
   
! TODO: I think this is just a reshape
   do iband=1,nband
     if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
     band_off=(iband-1)*2*nband
     do jband=1,nband
       eig1_k(2*jband-1+band_off)=eig1_k_loc(1,jband,iband)
       eig1_k(2*jband  +band_off)=eig1_k_loc(2,jband,iband)
     end do
#ifdef DEV_MJV
print *, 'me ', me_band, ' eig1_k 693 line for iband ', iband, eig1_k(band_off+1:band_off+2*nband)
#endif
   end do
 end if ! ipert/=natom+10.and.ipert/=natom+11

 ! Filter the wavefunctions for large modified kinetic energy (see routine mkkin.f)
! TODO: should this also be applied to cwaveq for the preconditioning with kinpw1 below?
 do ispinor=1,nspinor
   ipws=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(cwavef,kinpw1,ipws,npw1)
   do ipw=1+ipws,npw1+ipws
     if(kinpw1(ipw-ipws)>huge(zero)*1.d-11)then
       cwavef(1:2,ipw)=zero
     end if
   end do
 end do
#ifdef DEV_MJV
print *, 'cwavef 670 ', cwavef(:,1:5)
print *, 'mpi_enreg%nproc_band ', mpi_enreg%nproc_band
#endif

 ! Apply the orthogonality condition: <C1 k,q|C0 k+q>=0 (NCPP) or <C1 k,q|S0|C0 k+q>=0 (PAW)
 ! Project out all bands from cwavef, i.e. apply P_c projector on cwavef
 ! (this is needed when there are some partially or unoccupied states)
 do iband = 1, nband
   if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
   if (iband == band) then
     work = cwavef
   end if
   call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
 
   call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
     gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)

   call xmpi_sum(work,mpi_enreg%comm_band,ierr)

! save this for me_band only
   if (iband == band) then
!TODO: make this a blas call? zaxpy
     cwavef = work - (mpi_enreg%nproc_band-1)*cwavef
   end if
 end do
#ifdef DEV_MJV
print *, 'cwavef 692 ', cwavef(:,1:5)
#endif


 if(ipert/=natom+10.and.ipert/=natom+11) then
   ! If PAW, the orthogonality condition is <C1 k,q|S0|C0 k+q>+1/2<C0 k|S1|C0 k+q>=0
   if (usepaw==1.and.usedcwavef>0) then
!$OMP PARALLEL DO
     do ipw=1,npw1*nspinor
       cwavef(1:2,ipw)=cwavef(1:2,ipw)+dcwavef(1:2,ipw)
     end do
#ifdef DEV_MJV
print *, 'dcwavef 739 ', dcwavef(:,1:5)
#endif
   end if
 else
   ! In 2nd order case, dcwavef/=0 even in NC, and it is already computed in rf2_init (called in dfpt_vtowfk.F90)
   do ipw=1,npw1*nspinor
     cwavef(:,ipw)=cwavef(:,ipw)+rf2%dcwavef(:,ipw+dc_shift_band)
   end do
 end if

#ifdef DEV_MJV
print *, 'cwavef 739 ', cwavef(:,1:5)
#endif
 if(band>max(1,nband-nbdbuf))then
   ! Treat the case of buffer bands
   cwavef=zero
   ghc   =zero
   gvnlxc =zero
   if (gen_eigenpb) gsc=zero
   if (usedcwavef==2) dcwavef=zero
   if (usepaw==1) then
     call pawcprj_set_zero(cwaveprj)
   end if
   if (usedcwavef==2) dcwavef=zero
   ! Number of one-way 3D ffts skipped
   nskip=nskip+nline

   ! At the end of the treatment of a set of bands, write the number of one-way 3D ffts skipped
   if (xmpi_paral==1 .and. band==nband .and. prtvol>=10) then
     write(msg,'(a,i0)')' dfpt_cgwf: number of one-way 3D ffts skipped in cgwf3 until now =',nskip
     call wrtout(std_out,msg)
   end if
  
   skipme = 1
#ifdef DEV_MJV
print *, 'skipping for buffer band'
#endif
 end if

 ! If not a buffer band, perform the optimisation

 ABI_ALLOCATE(conjgr,(2,npw1*nspinor))
 ABI_ALLOCATE(direc,(2,npw1*nspinor))
 ABI_ALLOCATE(gresid,(2,npw1*nspinor))
 ABI_ALLOCATE(cwaveq,(2,npw1*nspinor))
 if (usepaw==1) then
   ABI_DATATYPE_ALLOCATE(conjgrprj,(natom,nspinor))
   call pawcprj_alloc(conjgrprj,0,gs_hamkq%dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(conjgrprj,(0,0))
 end if

 cwaveq(:,:)=cgq(:,1+npw1*nspinor*(band_me-1)+icgq:npw1*nspinor*band_me+icgq)
 dotgp=one

#ifdef DEV_MJV
if (gen_eigenpb) then
print *, 'cwave 807 ', cwaveprj(1,1)%cp(:,:)
end if
#endif
 ! Here apply H(0) at k+q to input orthogonalized 1st-order wfs
 sij_opt=0;if (gen_eigenpb) sij_opt=1
 cpopt=-1+usepaw
 call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamkq,gvnlxc,eshift,mpi_enreg,1,&
   prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

 ! ghc also includes the eigenvalue shift
 if (gen_eigenpb) then
   call cg_zaxpy(npw1*nspinor, [-eshift, zero], gsc,ghc)
 else
   call cg_zaxpy(npw1*nspinor, [-eshift, zero], cwavef,ghc)
 end if
#ifdef DEV_MJV
print *, 'cwavef 763 ', cwavef(:,1:10)
print *, 'gsc 778 ', band, gsc(:,1:10)
print *, 'ghc 778 ',band,  ghc(:,1:10)
print *, 'gh1c 778 ',band,  gh1c(:,1:10)
#endif

 ! Initialize resid, in case of nline==0
 resid=zero



 ! ======================================================================
 ! ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
 ! ======================================================================
 do iline=1,nline

   bands_skipped_now = 0

   ! ======================================================================
   ! ================= COMPUTE THE RESIDUAL ===============================
   ! ======================================================================
   ! Note that gresid (=steepest-descent vector, Eq.(26) of PRB 55, 10337 (1996) [[cite:Gonze1997]])
   ! is precomputed to garantee cancellation of errors
   ! and allow residuals to reach values as small as 1.0d-24 or better.
   if (berryopt== 4.or.berryopt== 6.or.berryopt== 7.or. berryopt==14.or.berryopt==16.or.berryopt==17) then
     if (ipert==natom+2) then
       if (opt_gvnlx1/=1) gvnlx1=zero
!$OMP PARALLEL DO
       do ipw=1,npw1*nspinor
         gresid(1:2,ipw)=-ghc(1:2,ipw)-gh1c(1:2,ipw)
       end do
     else
!$OMP PARALLEL DO
       do ipw=1,npw1*nspinor
         gresid(1,ipw)=-ghc(1,ipw)-gh1c(1,ipw)+gberry(2,ipw)
         gresid(2,ipw)=-ghc(2,ipw)-gh1c(2,ipw)-gberry(1,ipw)
       end do
     end if
   else
!$OMP PARALLEL DO
     do ipw=1,npw1*nspinor
       gresid(1:2,ipw)=-ghc(1:2,ipw)-gh1c(1:2,ipw)
     end do
   end if

#ifdef DEV_MJV
print *, 'band gresid 821 ', band, gresid(:,1:5)
#endif
   ! ======================================================================
   ! =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
   ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
   ! ======================================================================
   ! Project all bands from gresid into direc:
   ! The following projection over the subspace orthogonal to occupied bands
   ! is not optional in the RF case, unlike the GS case.
   ! However, the order of operations could be changed, so that
   ! as to make it only applied at the beginning, to H(1) psi(0),
   ! so, THIS IS TO BE REEXAMINED
   ! Note the subtlety:
   ! -For the generalized eigenPb, S|cgq> is used in place of |cgq>,
   ! in order to apply P_c+ projector (see PRB 73, 235101 (2006) [[cite:Audouze2006]], Eq. (71), (72)
   do iband = 1, nband
     work = zero
     if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
     if (iband == band) then
       work = gresid
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
  
     if(gen_eigenpb)then
       call projbd(gscq,work,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
         cgq,  scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     else
       call projbd( cgq,work,-1, icgq,   0,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
         dummy,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     end if
  
       call xmpi_sum(work,mpi_enreg%comm_band,ierr)
    
! save this for me_band only
     if (iband == band) then 
!TODO: make this a blas call? zaxpy
       gresid = work - (mpi_enreg%nproc_band-1)*gresid
     end if
   end do

   call cg_zcopy(npw1*nspinor,gresid,direc)
#ifdef DEV_MJV
print *, ' cgwf band gresid 906 ', band, gresid (:,1:5)
#endif

   ! ======================================================================
   ! ============== CHECK FOR CONVERGENCE CRITERIA ========================
   ! ======================================================================

   ! Compute second-order derivative of the energy using a variational expression
   call dotprod_g(prod1,doti,istwf_k,npw1*nspinor,1,cwavef,gresid,me_g0,mpi_enreg%comm_spinorfft)
   call dotprod_g(prod2,doti,istwf_k,npw1*nspinor,1,cwavef,gh1c,me_g0,mpi_enreg%comm_spinorfft)
   d2te=two*(-prod1+prod2)
   ! write(std_out,'(a,f14.6,a,f14.6)') 'prod1 = ',prod1,' prod2 = ',prod2 ! Keep this debugging feature!

   ! Compute residual (squared) norm
   call sqnorm_g(resid,istwf_k,npw1*nspinor,gresid,me_g0,comm_fft)
#ifdef DEV_MJV
print *, ' cgwf band, resid ', band, resid
#endif
   if (prtvol==-level.or.prtvol==-19)then
     write(msg,'(a,a,i3,f14.6,a,a,4es12.4)') ch10,&
      ' dfpt_cgwf : iline,eshift     =',iline,eshift,ch10,&
      '         resid,prod1,prod2,d2te=',resid,prod1,prod2,d2te
     call wrtout(std_out,msg)
   end if

   ! Compute <u_m(1)|H(0)-e_m(0)|u_m(1)>
   ! (<u_m(1)|H(0)-e_m(0).S|u_m(1)> if gen. eigenPb),
   ! that should be positive,
   ! except when the eigenvalue eig_mk(0) is higher than
   ! the lowest non-treated eig_mk+q(0). For insulators, this
   ! has no influence, but for metallic occupations,
   ! the conjugate gradient algorithm breaks down. The solution adopted here
   ! is very crude, and rely upon the fact that occupancies of such
   ! levels should be smaller and smaller with increasing nband, so that
   ! a convergence study will give the right result.
   ! The same trick is also used later.
   u1h0me0u1=-prod1-prod2

   ! Some tolerance is allowed, to account for very small numerical inaccuracies and cancellations.
   if(u1h0me0u1<-tol_restart .and. skipme == 0)then
     if (prtvol==-level.or.prtvol==-19) then
       write(msg,'(a,es22.13e3)') '  cgwf3: u1h0me0u1 = ',u1h0me0u1
       call wrtout(std_out,msg)
     end if
     cwavef =zero
     ghc    =zero
     gvnlxc  =zero
     if (gen_eigenpb) gsc(:,:)=zero
     if (usepaw==1) then
       call pawcprj_set_zero(cwaveprj)
     end if
     ! A negative residual will be the signal of this problem ...
     resid=-one
     if (prtvol > 0) call wrtout(std_out,' dfpt_cgwf: problem of minimisation (likely metallic), set resid to -1')
     ! Number of one-way 3D ffts skipped
     nskip=nskip+(nline-iline+1)

!DEBUG     exit ! Exit from the loop on iline
     skipme = 1
#ifdef DEV_MJV
print *, 'skipping for u1h0me0u1'
#endif
   end if

   ! If residual sufficiently small stop line minimizations
   if (resid<tolwfr .and. skipme == 0) then
     if(prtvol>=10)then
       write(msg,'(a,i4,a,i2,a,es12.4)')' dfpt_cgwf: band',band,' converged after ',iline,' line minimizations: resid = ',resid
       call wrtout(std_out,msg)
     end if
     nskip=nskip+(nline-iline+1)  ! Number of two-way 3D ffts skipped
!DEBUG     exit                         ! Exit from the loop on iline
     skipme = 1
#ifdef DEV_MJV
print *, 'skipping for resid'
#endif
   end if

   ! If user require exiting the job, stop line minimisations
   if (quit==1) then
     write(msg,'(a,i0)')' dfpt_cgwf: user require exiting => skip update of band ',band
     call wrtout(std_out,msg)
     nskip=nskip+(nline-iline+1)  ! Number of two-way 3D ffts skipped
     exit                         ! Exit from the loop on iline
   end if

   ! Check that d2te is decreasing on succeeding lines:
   if (iline/=1) then
     if (d2te>d2teold+tol6) then
       write(msg,'(a,i0,a,e14.6,a,e14.6)')'New trial energy at line ',iline,'=',d2te,'is higher than former:',d2teold
       MSG_WARNING(msg)
     end if
   end if
   d2teold=d2te

   !DEBUG Keep this debugging feature !
   !call sqnorm_g(dotr,istwf_k,npw1*nspinor,direc,me_g0,comm_fft)
   !write(std_out,*)' dfpt_cgwf : before precon, direc**2=',dotr
   !if (gen_eigenpb) then
   !call dotprod_g(dotr,doti,istwf_k,npw1*nspinor,1,cwaveq,&
   !&                 gscq(:,1+npw1*nspinor*(band-1)+igscq:npw1*nspinor*band+igscq),me_g0,mpi_enreg%comm_spinorfft)
   !else
   !call sqnorm_g(dotr,istwf_k,npw1*nspinor,cwaveq,me_g0,comm_fft)
   !end if
   !write(std_out,*)' dfpt_cgwf : before precon, cwaveq**2=',dotr
   !ENDDEBUG

   ! ======================================================================
   ! ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
   ! ======================================================================

   ! If wfoptalg>=10, the precondition matrix is kept constant
   ! during iteration; otherwise it is recomputed
   if (wfoptalg<10.or.iline==1) then
!print *, 'cwaveq cg_precon ', cwaveq
!print *, 'direc cg_precon ', direc
     call cg_precon(cwaveq,zero,istwf_k,kinpw1,npw1,nspinor,me_g0,0,pcon,direc,mpi_enreg%comm_fft)
   else
     do ispinor=1,nspinor
       igs=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(igs,npw,direc,pcon)
       do ipw=1+igs,npw1+igs
         direc(1:2,ipw)=direc(1:2,ipw)*pcon(ipw-igs)
       end do
     end do
   end if

   !DEBUG Keep this debugging feature !
   !call sqnorm_g(dotr,istwf_k,npw1*nspinor,direc,me_g0,comm_fft)
   !write(std_out,*)' dfpt_cgwf : after precon, direc**2=',dotr
   !ENDDEBUG

   ! ======================================================================
   ! ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
   ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
   ! ======================================================================

   ! Projecting again out all bands:
   ! -For the simple eigenPb, gscq is used as dummy argument
   do iband = 1, nband
     if (bands_treated_now(iband)-bands_skipped_now(iband) == 0) cycle
     if (iband == band) then
       work = direc
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
  
     call projbd(cgq,work,-1,icgq,igscq,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
       gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
  
     call xmpi_sum(work,mpi_enreg%comm_band,ierr)
  
! save this for me_band only
     if (iband == band) then 
!TODO: make this a blas call? zaxpy
       direc = work - (mpi_enreg%nproc_band-1)*direc
     end if
   end do


   !DEBUG Keep this debugging feature !
   !call sqnorm_g(dotr,istwf_k,npw1*nspinor,direc,me_g0,comm_fft)
   !write(std_out,*)' dfpt_cgwf : after projbd, direc**2=',dotr
   !ENDDEBUG

   ! ======================================================================
   ! ================= COMPUTE THE CONJUGATE-GRADIENT =====================
   ! ======================================================================

   ! get dot of direction vector with residual vector
   call dotprod_g(dotgg,doti,istwf_k,npw1*nspinor,1,direc,gresid,me_g0,mpi_enreg%comm_spinorfft)

   if (iline==1) then
     ! At first iteration, gamma is set to zero
     gamma=zero
     dotgp=dotgg
     call cg_zcopy(npw1*nspinor,direc,conjgr)
   else
     ! At next iterations, h = g + gamma * h
     gamma=dotgg/dotgp
     dotgp=dotgg
     if (prtvol==-level.or.prtvol==-19)then
       write(msg,'(a,2es16.6)') 'dfpt_cgwf: dotgg,gamma = ',dotgg,gamma
       call wrtout(std_out,msg)
     end if
!$OMP PARALLEL DO
     do ipw=1,npw1*nspinor
       conjgr(1:2,ipw)=direc(1:2,ipw)+gamma*conjgr(1:2,ipw)
     end do
     if (prtvol==-level.or.prtvol==-19) call wrtout(std_out,'dfpt_cgwf: conjugate direction has been found')
   end if

   ! ======================================================================
   ! ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
   ! ======================================================================
   ! ...along the search direction

   ! Compute dedt, Eq.(29) of of PRB55, 10337 (1997) [[cite:Gonze1997]],
   ! with an additional factor of 2 for the difference between E(2) and the 2DTE
   dedt = zero
   call dotprod_g(dedt,doti,istwf_k,npw1*nspinor,1,conjgr,gresid,me_g0,mpi_enreg%comm_spinorfft)
   dedt=-two*two*dedt

#ifdef DEV_MJV
print *, 'gresid 781 ', gresid(1:2,1:5)
#endif
   if((prtvol==-level.or.prtvol==-19.or.prtvol==-20).and.dedt-tol14>0) call wrtout(std_out,' CGWF3_WARNING : dedt>0')
   ABI_ALLOCATE(gvnlx_direc,(2,npw1*nspinor))
   ABI_ALLOCATE(gh_direc,(2,npw1*nspinor))
   if (gen_eigenpb)  then
     ABI_ALLOCATE(sconjgr,(2,npw1*nspinor))
   else
     ABI_ALLOCATE(sconjgr,(0,0))
   end if
   sij_opt=0;if (gen_eigenpb) sij_opt=1
   cpopt=-1+usepaw
   call getghc(cpopt,conjgr,conjgrprj,gh_direc,sconjgr,gs_hamkq,gvnlx_direc,&
     eshift,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

   ! ghc also includes the eigenvalue shift
   if (gen_eigenpb) then
!$OMP PARALLEL DO
     do ipw=1,npw1*nspinor
       gh_direc(1:2,ipw)=gh_direc(1:2,ipw)-eshift*sconjgr(1:2,ipw)
     end do
   else
!$OMP PARALLEL DO
     do ipw=1,npw1*nspinor
       gh_direc(1:2,ipw)=gh_direc(1:2,ipw)-eshift*conjgr(1:2,ipw)
     end do
   end if
#ifdef DEV_MJV
print *, 'gh_direc ', gh_direc(1:2,1:5)
#endif

   ! compute d2edt2, Eq.(30) of of PRB55, 10337 (1997) [[cite:Gonze1997]],
   ! with an additional factor of 2 for the difference
   ! between E(2) and the 2DTE, and neglect of local fields (SC terms)
   d2edt2 = zero
   call dotprod_g(d2edt2,doti,istwf_k,npw1*nspinor,1,conjgr,gh_direc,me_g0,mpi_enreg%comm_spinorfft)
   d2edt2=two*two*d2edt2
   if(prtvol==-level.or.prtvol==-19)then
     write(msg,'(a,2es14.6)') 'dfpt_cgwf: dedt,d2edt2=',dedt,d2edt2
     call wrtout(std_out,msg)
   end if

   ! ======================================================================
   ! ======= COMPUTE MIXING FACTOR - CHECK FOR CONVERGENCE ===============
   ! ======================================================================

#ifdef DEV_MJV
print *, 'dedt, d2edt2, theta 1120 ', dedt, d2edt2, theta
#endif

   ! see Eq.(31) of PRB55, 10337 (1997) [[cite:Gonze1997]]
   !
   if(d2edt2<-tol_restart .and. skipme == 0)then
     ! This may happen when the eigenvalue eig_mk(0) is higher than
     ! the lowest non-treated eig_mk+q(0). The solution adopted here
     ! is very crude, and rely upon the fact that occupancies of such
     ! levels should be smaller and smaller with increasing nband, so that
     ! a convergence study will give the right result.
     theta=zero

     cwavef=zero
     ghc   =zero
     gvnlxc =zero
     if (gen_eigenpb) gsc=zero
     if (usepaw==1) then
       call pawcprj_set_zero(cwaveprj)
     end if
     ! A negative residual will be the signal of this problem ...
     resid=-two
     if (prtvol > 0) call wrtout(std_out,' dfpt_cgwf: problem of minimisation (likely metallic), set resid to -2')
   else if (d2edt2 > 1.d-40) then
     ! Here, the value of theta that gives the minimum
     theta=-dedt/d2edt2
     !write(std_out,*)' dfpt_cgwf: dedt,d2edt2=',dedt,d2edt2
   else
     write(msg,'(a)') 'DFPT_CGWF WARNING : d2edt2 is zero, skipping update'
     call wrtout(std_out,msg,'COLL')
     theta=zero
#ifdef DEV_MJV
print *, 'dedt, d2edt2, theta 1152 ', dedt, d2edt2, theta
#endif
   end if

   ! Check that result is above machine precision
   if (one+theta==one) then
     if (prtvol > 0) then
       write(msg, '(a,es16.4)' ) ' dfpt_cgwf: converged with theta=',theta
       call wrtout(std_out,msg)
     end if
     nskip=nskip+2*(nline-iline) ! Number of one-way 3D ffts skipped
     skipme = 1
#ifdef DEV_MJV
print *, 'skipping for theta below machine prec'
#endif
!DEBUG     exit                        ! Exit from the loop on iline
   end if

   ! ======================================================================
   ! ================ GENERATE NEW |wf>, H|wf>, Vnl|Wf ... ================
   ! ======================================================================
#ifdef DEV_MJV
print *, 'skipme ', skipme
print *, ' cwavef 1174 ', cwavef(:,1:5)
#endif

   if (skipme == 0) then
     call cg_zaxpy(npw1*nspinor, [theta, zero], conjgr,cwavef)
 ! Filter the wavefunctions for large modified kinetic energy (see routine mkkin.f)
     do ispinor=1,nspinor
       ipws=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(cwavef,kinpw1,ipws,npw1)
       do ipw=1+ipws,npw1+ipws
         if(kinpw1(ipw-ipws)>huge(zero)*1.d-11)then
           cwavef(1:2,ipw)=zero
         end if
       end do
     end do

     call cg_zaxpy(npw1*nspinor, [theta, zero], gh_direc,ghc)
     call cg_zaxpy(npw1*nspinor, [theta, zero], gvnlx_direc,gvnlxc)
  
     if (gen_eigenpb) then
       call cg_zaxpy(npw1*nspinor, [theta, zero], sconjgr, gsc)
     end if
     if (usepaw==1) then
       call pawcprj_axpby(theta,one,conjgrprj,cwaveprj)
     end if
   end if

   ABI_DEALLOCATE(gh_direc)
   ABI_DEALLOCATE(gvnlx_direc)
   ABI_DEALLOCATE(sconjgr)

   ! ======================================================================
   ! =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
   ! ======================================================================

   if(usetolrde/=0) then
     ! Check reduction in trial energy deltae, Eq.(28) of PRB55, 10337 (1997) [[cite:Gonze1997]]
     deltae=half*d2edt2*theta**2+theta*dedt

     if (iline==1) then
       deold=deltae
       ! The extra factor of two should be removed !
     else if (abs(deltae)<tolrde*two*abs(deold) .and. iline/=nline) then
       if(prtvol>=10.or.prtvol==-level.or.prtvol==-19)then
         write(msg, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
          ' dfpt_cgwf: line',iline,' deltae=',deltae,' < tolrde*',deold,' =>skip lines'
         call wrtout(std_out,msg)
       end if
       nskip=nskip+2*(nline-iline) ! Number of one-way 3D ffts skipped
       skipme = 1
#ifdef DEV_MJV
print *, 'skipping for deltae diff'
#endif
!DEBUG       exit                        ! Exit from the loop on iline
     end if
   end if

! if all bands are skippable, we can exit the iline loop for good.
!   otherwise, all procs are needed for the projbd and other operations,
!   even if the present band will not be updated
   bands_skipped_now(band) = skipme
   call xmpi_sum(bands_skipped_now,mpi_enreg%comm_band,ierr)
#ifdef DEV_MJV
print *, 'bands_skipped_now 1 ', bands_skipped_now
#endif
   bands_skipped_now = bands_skipped_now - bands_treated_now
#ifdef DEV_MJV
print *, 'bands_skipped_now 2 ', bands_skipped_now
#endif
   if (sum(abs(bands_skipped_now)) == 0) exit

   ! ======================================================================
   ! ================== END LOOP FOR GIVEN BAND ===========================
   ! ======================================================================

   ! Note that there are five "exit" instruction inside the loop.
   nlines_done = nlines_done + 1
#ifdef DEV_MJV
print *, 'cgwf  cwavef  1180 ',cwavef(:,1:5), cwavef(:,23)
#endif
 end do ! iline
#ifdef DEV_MJV
print *, 'shape ', shape(cwavef)
print *, 'cgwf band,  cwavef', band, cwavef(:,1:5)
print *, 'cgwf  cwavef  1183 ',cwavef(:,23)
print *, 'cgwf band,  ghc', band, ghc(:,1:5)
#endif

!--------------------------------------------------------------------------
!             DEBUG
!--------------------------------------------------------------------------
 ! Check that final cwavef (Psi^(1)) satisfies the orthogonality condition
 if (prtvol==-level.or.prtvol==-19) then
   sij_opt=0 ; usevnl=1 ; optlocal=1 ; optnl=2 ; if (gen_eigenpb)  sij_opt=1
   ABI_ALLOCATE(work2,(2,npw1*nspinor*sij_opt))
   iband_me = 0
   do iband=1,nband
     if (bands_treated_now(iband) == 0) cycle
     if (band_procs(iband)==me_band) then 
       iband_me = iband_me+1
       if (gen_eigenpb) then
         work(:,:)=gscq(:,1+npw1*nspinor*(iband-1)+igscq:npw1*nspinor*iband+igscq)
       else
         work(:,:)=cgq(:,1+npw1*nspinor*(iband-1)+icgq:npw1*nspinor*iband+icgq)
       end if
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)

     ! Compute: <Psi^(0)_i,k+q|Psi^(1)_j,k,q>
     call dotprod_g(prod1,prod2,istwf_k,npw1*nspinor,2,work,cwavef,me_g0,mpi_enreg%comm_spinorfft)

     if (ipert/=natom+10.and.ipert/=natom+11) then
       if (gen_eigenpb) then
         call getgh1c(berryopt,cwave0,cwaveprj0,work1,gberry,work2,gs_hamkq,gvnlx1_saved,idir,ipert,eshift,&
           mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

         if (band_procs(iband)==me_band) then 
           work(:,:)=cgq(:,1+npw1*nspinor*(iband_me-1)+icgq:npw1*nspinor*iband_me+icgq)
         end if
         call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
         call dotprod_g(dotr,doti,istwf_k,npw1*nspinor,2,work,work2,me_g0,mpi_enreg%comm_spinorfft)
       else
         dotr=zero; doti=zero
       end if
       dotr=prod1+half*dotr
       doti=prod2+half*doti
     else if(prtvol==-19) then ! 2nd order case
       dotr=prod1+half*rf2%amn(1,iband+(band-1)*nband)
       doti=prod2+half*rf2%amn(2,iband+(band-1)*nband)
     else
       write(msg,'(a)') 'CGWF3_WARNING : Use prtvol=-19 to test orthogonality for ipert=natom+10 or +11'
       call wrtout(std_out,msg,'COLL')
     end if
     dotr=sqrt(dotr**2+doti**2)
     if(dotr>tol10) then
!         if (gen_eigenpb) then
!           write(msg,'(2a,i3,a,2es22.15)') 'CGWF3_WARNING : <Psi^(1)_i,k,q|S^(0)|Psi^(0)_j,k+q>',&
!             '+ 1/2<Psi^(0)_i,k|S^(1)|Psi^(0)_j,k+q>, for j= ',iband,' is ',dotr,doti
!         else
       write(msg,'(a,i3,a,es22.15)') 'CGWF3_WARNING : |<Psi^(0)_i,k+q|Psi^(1)_j,k,q>+amn(i,j)/2|, for j= ',iband,' is ',dotr
!         end if
       call wrtout(std_out,msg)
     end if
   end do
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   if (ipert/=natom+10.and.ipert/=natom+11) then
     ABI_DEALLOCATE(gvnlx1_saved)
   end if
 end if ! prtvol==-level.or.prtvol==-19

 if (prtvol==-level.or.prtvol==-19)then
   !  Check that final cwavef Psi^(1) is Pc.Psi^(1)+delta_Psi^(1)
   ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
   ! -Apply Pc to Psi^(1)
   do iband = 1, nband
     if (bands_treated_now(iband) == 0) cycle
     if (iband == band) then
       cwwork=cwavef
     end if
     call xmpi_bcast(cwwork,band_procs(iband),mpi_enreg%comm_band,ierr)

     if(gen_eigenpb)then
       call projbd(cgq,cwwork,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
         cgq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     else
       call projbd(cgq,cwwork,-1,icgq,0,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
         dummy,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     end if

     call xmpi_sum(cwwork,mpi_enreg%comm_band,ierr)

! save this for me_band only
     if (iband == band) then 
!TODO: make this a blas call? zaxpy
       cwwork = cwwork - (mpi_enreg%nproc_band-1)*cwavef
     end if
   end do

   ! -Add delta_Psi^(1)
   if (usedcwavef>0) cwwork=cwwork+dcwavef
   if(ipert==natom+10.or.ipert==natom+11) cwwork=cwwork+rf2%dcwavef(:,1+dc_shift_band:npw1*nspinor+dc_shift_band)
   ! -Compare to Psi^(1)
   cwwork=cwwork-cwavef
   call sqnorm_g(dotr,istwf_k,npw1*nspinor,cwwork,me_g0,comm_fft)
   ABI_DEALLOCATE(cwwork)
   if(sqrt(dotr)>tol10) then
!       if (gen_eigenpb) then
!         write(msg,'(a,i3,a,es22.15)') &
!         'CGWF3_WARNING : |(Pc.Psi^(1)_i,k,q + delta_Psi^(1)_i,k) - Psi^(1)_i,k,q|^2 (band ',band,')=',dotr
!       else
     write(msg,'(a,es22.15)') 'CGWF3_WARNING : |(Pc.Psi^(1)_i,k,q + delta_Psi^(1)_i,k) - Psi^(1)_i,k,q| = ',sqrt(dotr)
!       end if
     call wrtout(std_out,msg)
   end if
 end if  ! prtvol==-level.or.prtvol==-19

 if(prtvol==-level.or.prtvol==-19.or.prtvol==-20)then
   ! Check that final cwavef (Psi^(1)) solves the Sternheimer equation
   ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
   ! -Apply Pc to Psi^(1)
   do iband = 1, nband
     if (bands_treated_now(iband) == 0) cycle
     if (iband == band) then
       cwwork=cwavef
     end if
     call xmpi_bcast(cwwork,band_procs(iband),mpi_enreg%comm_band,ierr)

     if(gen_eigenpb)then
       call projbd(cgq,cwwork,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
         gscq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     else
       call projbd(cgq,cwwork,-1,icgq,0,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
         dummy,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     end if

     call xmpi_sum(cwwork,mpi_enreg%comm_band,ierr)

! save this for me_band only
     if (iband == band) then 
!TODO: make this a blas call? zaxpy
      cwwork = cwwork - (mpi_enreg%nproc_band-1)*cwavef
     end if
   end do

   ! - Apply H^(0)-E.S^(0)
   sij_opt=0;if (gen_eigenpb) sij_opt=1
   cpopt=-1
   ABI_ALLOCATE(work1,(2,npw1*nspinor*((sij_opt+1)/2)))
   ABI_ALLOCATE(work2,(2,npw1*nspinor))
   call getghc(cpopt,cwwork,conjgrprj,work,work1,gs_hamkq,work2,eshift,&
     mpi_enreg,1,prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   if (gen_eigenpb) then
     cwwork=work-eshift*work1
   else
     cwwork=work-eshift*cwwork
   end if
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
  
   ! The following is not mandatory, as Pc has been already applied to Psi^(1)
   ! and Pc^* H^(0) Pc = Pc^* H^(0) = H^(0) Pc (same for S^(0)).
   ! However, in PAW, to apply Pc^* here seems to reduce the numerical error
   ! -Apply Pc^*
   do iband = 1, nband
     if (bands_treated_now(iband) == 0) cycle
     if (iband == band) then
       work=cwwork
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)
  
     if(gen_eigenpb)then
       call projbd(gscq,  work,-1,igscq,icgq,istwf_k,mgscq,mcgq,nband_me,npw1,nspinor,&
         cgq,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     else
       call projbd(cgq,work,-1,icgq,0,istwf_k,mcgq,mgscq,nband_me,npw1,nspinor,&
         dummy,scprod,0,tim_projbd,useoverlap,me_g0,comm_fft)
     end if
  
     call xmpi_sum(work,mpi_enreg%comm_band,ierr)

! save this for me_band only
     if (iband == band) then 
!TODO: make this a blas call? zaxpy
       cwwork = cwwork - (mpi_enreg%nproc_band-1)*work
     end if
   end do

   ! - Add Pc^*(H^(1)-E.S^(1)).Psi^(0)
   cwwork=cwwork+gh1c
   call sqnorm_g(dotr,istwf_k,npw1*nspinor,cwwork,me_g0,comm_fft)
   ABI_DEALLOCATE(cwwork)
   write(msg,'(a,i3,a,es22.15,2a,i4)') &
     '*** CGWF3 Sternheimer equation test for band ',band,'=',sqrt(dotr),ch10,&
     'It should go to zero for large nline : nlines_done = ',nlines_done
   call wrtout(std_out,msg)
 end if ! prtvol==-level.or.prtvol==-19.or.prtvol==-20

 if(prtvol==-level.or.prtvol==-19.or.prtvol==-20)then
   ! Check that < Psi^(0) | ( H^(0)-eps^(0) S^(0) ) | Psi^(1) > is in agreement with eig^(1)
   ABI_ALLOCATE(cwwork,(2,npw1*nspinor))
   cwwork=cwavef
   ! - Apply H^(0)-E.S^(0)
   sij_opt=0;if (gen_eigenpb) sij_opt=1
   cpopt=-1
   ABI_ALLOCATE(work1,(2,npw1*nspinor*((sij_opt+1)/2)))
   ABI_ALLOCATE(work2,(2,npw1*nspinor))
   call getghc(cpopt,cwwork,conjgrprj,work,work1,gs_hamkq,work2,eshift,&
     mpi_enreg,1,prtvol,sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   if (gen_eigenpb) then
     cwwork=work-eshift*work1
   else
     cwwork=work-eshift*cwwork
   end if
   ABI_DEALLOCATE(work1)
   ABI_DEALLOCATE(work2)
   cwwork=cwwork+gh1c_n
   jband=(band-1)*2*nband
   iband_me = 0
   do iband=1,nband
     if (bands_treated_now(iband) == 0) cycle
     if (band_procs(iband)==me_band) then
       iband_me = iband_me+1
       work(:,:)=cgq(:,1+npw1*nspinor*(iband_me-1)+icgq:npw1*nspinor*iband_me+icgq)
     end if
     call xmpi_bcast(work,band_procs(iband),mpi_enreg%comm_band,ierr)

     call dotprod_g(dotr,doti,istwf_k,npw1*nspinor,2,work,cwwork,me_g0,mpi_enreg%comm_spinorfft)
     dotr = dotr - eig1_k(2*iband-1+jband)
     doti = doti - eig1_k(2*iband  +jband)
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol8) then
       write(msg,'(2(a,i3),a,es22.15)') &
         'CGWF3_WARNING < Psi^(0) | ( H^(0)-eps^(0) S^(0) ) | Psi^(1) > for i=',iband,' j=',band,&
       ' : ',sqrt(dotr**2+doti**2)
       call wrtout(std_out,msg)
     end if
   end do
   !write(std_out,'(a)') '< Psi^(0) | ( H^(0)-eps^(0) S^(0) ) | Psi^(1) > is done.'
   ABI_DEALLOCATE(cwwork)
 end if ! prtvol==-level.or.prtvol==-19.or.prtvol==-20
!--------------------------------------------------------------------------
!            END DEBUG
!--------------------------------------------------------------------------

 if (allocated(gh_direc))  then
   ABI_DEALLOCATE(gh_direc)
 end if
 if (allocated(gvnlx_direc))  then
   ABI_DEALLOCATE(gvnlx_direc)
 end if
 ABI_DEALLOCATE(conjgr)
 ABI_DEALLOCATE(cwaveq)
 ABI_DEALLOCATE(direc)
 ABI_DEALLOCATE(gresid)
 if (usepaw==1) then
   call pawcprj_free(conjgrprj)
 end if
 ABI_DATATYPE_DEALLOCATE(conjgrprj)

#ifdef DEV_MJV
print *, 'cgwf band, resid 1553   ', band, resid
#endif
 if(band>max(1,nband-nbdbuf))then
   ! A small negative residual will be associated with these
   ! in the present algorithm all bands need to be in the loops over cgq etc... for the parallelization
   resid=-0.1_dp
 end if

 ! At the end of the treatment of a set of bands, write the number of one-way 3D ffts skipped
 if (xmpi_paral==1 .and. band==nband .and. prtvol>=10) then
   write(msg,'(a,i0)')' dfpt_cgwf: number of one-way 3D ffts skipped in cgwf3 until now =',nskip
   call wrtout(std_out,msg)
 end if

 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(gh1c)
 ABI_DEALLOCATE(pcon)
 ABI_DEALLOCATE(scprod)
 ABI_DEALLOCATE(gberry)

 call timab(122,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_cgwf
!!***

end module m_dfpt_cgwf
!!***
