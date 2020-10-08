!!****m* ABINIT/m_cgwfnew
!! NAME
!!  m_cgwfnew
!!
!! FUNCTION
!!  Conjugate-gradient eigensolver.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (DCA, XG, GMR, MT, MVeithen, ISouza, JIniguez)
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

module m_cgwfnew

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_cgtools
 use m_efield

 use defs_abitypes,   only : MPI_type
 use m_time,          only : timab
 use m_numeric_tools, only : rhophi
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_put, pawcprj_copy, &
                             pawcprj_get, pawcprj_mpi_allgather, pawcprj_free, pawcprj_symkn
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_fock,          only : fock_set_ieigen,fock_set_getghc_call
 use m_getchc,        only : getchc
 use m_getghc,        only : getghc
 use m_berrytk,       only : smatrix
 use m_nonlop,        only : nonlop
 use m_paw_overlap,   only : smatrix_k_paw
 use m_cgprj,         only : getcprj
 use m_cgwf,          only : mksubham
 !LTEST
 use testing
 !LTEST

 implicit none

 private
!!***

 public :: cgwfnew

!!***

contains
!!***

!!****f* m_cgwf/cgwfnew
!! NAME
!! cgwfnew
!!
!! FUNCTION
!! Update all wavefunction |C>, non self-consistently.
!! also compute the corresponding H|C> and Vnl|C> (and S|C> if paw).
!! Uses a conjugate-gradient algorithm.
!! In case of paw, resolves a generalized eigenproblem using an
!!  overlap matrix (not used for norm conserving psps).
!!
!! INPUTS
!!  berryopt == 4/14: electric field is on;
!           6/7/16/7: electric displacement field is on;
!!              all other values, no field is present
!!  chkexit= if non-zero, check whether the user wishes to exit
!!  cpus = CPU time limit
!!  filnam_ds1=name of input file (used for exit checking)
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  ikpt=number of the k-point
!!  inonsc=index of non self-consistent loop
!!  isppol=spin polarization currently treated
!!  mband =maximum number of bands
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mgsc=second dimension of the gsc array
!!  mkgq = second dimension of pwnsfacq
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw
!!  nband=number of bands.
!!  nbdblock=number of bands in a block
!!  nkpt=number of k points
!!  nline=number of line minimizations per band.
!!  npw=number of planewaves in basis sphere at given k.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  prtvol=control print volume and debugging output
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!  tolrde=tolerance on the ratio of differences of energies (for the line minimisation)
!!  tolwfr=tolerance on largest wf residual
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  use_subvnlx=1 if subvnlx has to be computed
!!  wfoptalg=govern the choice of algorithm for wf optimisation
!!   (0, 1, 10 and 11 : in the present routine, usual CG algorithm ;
!!   (2 and 3 : use shifted square Hamiltonian)
!!  zshift(nband)=in case wfoptalg is 2 or 3, shift of the Hamiltonian
!!
!! OUTPUT
!!  dphase_k(3) = change in Zak phase for the current k-point in case berryopt = 4/14,6/16,7/17 (electric (displacement) field)
!!  resid(nband)=wf residual for new states=|(H-e)|C>|^2 (hartree^2)
!!  subham(nband*(nband+1))=Hamiltonian expressed in the WFs subspace
!!  subovl(nband*(nband+1)*use_subovl)=overlap matrix expressed in sthe WFs subspace
!!  subvnlx(nband*(nband+1)*use_subvnlx))=non-local Hamiltonian (if NCPP)  plus Fock ACE operator (if usefock_ACE)
!!   expressed in the WFs subspace
!!
!! SIDE EFFECTS
!!  cg(2,mcg)
!!    at input =wavefunction <G|C band,k> coefficients for ALL bands
!!    at output same as input except that
!!      the current band, with number 'band' has been updated
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  quit= if 1, proceeds to smooth ending of the job.
!!  if(gs_hamk%usepaw==1)
!!   gsc(2,mgsc)=<G|S|C band,k> coefficients for ALL bands
!!               where S is the overlap matrix (used only for paw)
!!
!! NOTES
!!  1) cg should not be filtered and normalized : it should already be OK at input !
!!  2) Not sure that that the generalized eigenproblem (when gs_hamk%usepaw=1)
!!     is compatible with wfoptalg=2 or 3 (use of shifted square  Hamiltonian) - to be verified
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      cg_precon,cg_zaxpy,cg_zcopy,cg_zscal,dotprod_g,etheta,fock_set_ieigen
!!      getcprj,getghc,linemin,make_grad_berry,mksubham,pawcprj_alloc
!!      pawcprj_copy,pawcprj_free,pawcprj_get,pawcprj_mpi_allgather,pawcprj_put
!!      pawcprj_symkn,projbd,smatrix,smatrix_k_paw,sqnorm_g,timab,wrtout
!!      xmpi_allgather
!!
!! SOURCE

subroutine cgwfnew(berryopt,cg,cgq,chkexit,cpus,dphase_k,dtefield,&
&                filnam_ds1,gsc,gs_hamk,icg,igsc,ikpt,inonsc,&
&                isppol,mband,mcg,mcgq,mgsc,mkgq,mpi_enreg,&
&                mpw,nband,nbdblock,nkpt,nline,npw,npwarr,&
&                nspinor,nsppol,ortalg,prtvol,pwind,&
&                pwind_alloc,pwnsfac,pwnsfacq,quit,resid,subham,subovl,&
&                subvnlx,tolrde,tolwfr,use_subovl,use_subvnlx,wfoptalg,zshift)

!Arguments ------------------------------------
 integer,intent(in) :: berryopt,chkexit,icg,igsc,ikpt,inonsc,isppol
 integer,intent(in) :: mband,mcg,mcgq,mgsc,mkgq,mpw,nband,nbdblock,nkpt,nline
 integer,intent(in) :: npw,nspinor,nsppol,ortalg,prtvol,pwind_alloc
 integer,intent(in) :: use_subovl,use_subvnlx,wfoptalg
 integer,intent(in) :: quit
 real(dp),intent(in) :: cpus,tolrde,tolwfr
 character(len=*),intent(in) :: filnam_ds1
 type(MPI_type),intent(in) :: mpi_enreg
 type(efield_type),intent(inout) :: dtefield
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp),intent(in) :: cgq(2,mcgq)
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq),zshift(nband)
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)
 real(dp),intent(inout) :: dphase_k(3)
 real(dp),intent(out) :: subham(nband*(nband+1)),subovl(nband*(nband+1)*use_subovl)
 real(dp),intent(out) :: subvnlx(nband*(nband+1)*use_subvnlx)
 real(dp),intent(out) :: resid(nband)

!Local variables-------------------------------
 integer,parameter :: level=113,tim_getghc=1,tim_projbd=1
 integer,save :: nskip=0
 integer :: choice,counter,cpopt,ddkflag,dimenlc1,dimenlr1,dimenl2,iat,iatom,itypat
 integer :: iband,ibandmin,ibandmax,me_g0
 integer :: ibdblock,iblock,icg1,icg_shift,icp1,icp2,idir,idum1,ierr,ifor,igs,igsc_shift,ii,ikgf
 integer :: ikpt2,ikpt2f,ikptf,iline,iproc,ipw,ispinor,istwf_k,isubh,isubo,itrs
 integer :: job,mcg_q,me_distrb,natom,ncpgr,nblock,nproc_distrb,npw_k2
 integer :: optekin,paw_opt,signs,shiftbd,sij_opt,spaceComm_distrb
 integer :: useoverlap,wfopta10
 real(dp) :: chc,costh,deltae,deold,dhc,dhd,diff,dotgg,dotgp,doti,dotr
 real(dp) :: dphase_aux2,e0,e0_old,e1,e1_old,eval,gamma
 real(dp) :: lam0,lamold,root,sinth,sintn,swap,tan2th,theta,thetam
 real(dp) :: xnorm
 logical :: gen_eigenpb
 character(len=500) :: message
 integer :: hel(2,3)
 integer,allocatable :: dimlmn(:),dimlmn_srt(:),ikptf_recv(:),pwind_k(:),sflag_k(:)
 real(dp) :: bcut(2,3),dphase_aux1(3),dtm_k(2),phase_end(3)
 real(dp) :: phase_init(3),tsec(2)
 real(dp),allocatable :: cg1_k(:,:),cgq_k(:,:),conjgr(:,:),cwavef(:,:)
 real(dp),allocatable :: detovc(:,:,:),detovd(:,:,:),direc(:,:),direc_tmp(:,:)
 real(dp),allocatable :: gh_direc(:,:),gh_direcws(:,:),ghc(:,:),ghc_all(:,:),ghcws(:,:)
 real(dp),allocatable :: grad_berry(:,:),grad_total(:,:),gs_direc(:,:)
 real(dp) :: gsc_dummy(0,0)
 real(dp),allocatable :: gvnlxc(:,:),gvnlx_direc(:,:),gvnlx_dummy(:,:)
 real(dp),allocatable :: pcon(:),pwnsfac_k(:,:),scprod(:,:),scwavef(:,:)
 real(dp),allocatable :: smat_inv(:,:,:),smat_k(:,:,:),smat_k_paw(:,:,:),swork(:,:),vresid(:,:),work(:,:)
 real(dp),pointer :: kinpw(:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:)
 type(pawcprj_type),allocatable :: cprj_direc(:,:),cprj_band_srt(:,:),cprj_gat(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Starting the routine
 call timab(22,1,tsec)

!Touching chkexit, cpus,filnam_ds to avoid warning for abirules. This is dirty...
 if(chkexit<0)then
   MSG_BUG('chkexit should be positive!')
 end if

 if(cpus<0 .and. filnam_ds1=='a')then
   MSG_BUG('cpus should be positive!')
 end if

!======================================================================
!========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
!======================================================================

!MPI data
 spaceComm_distrb=mpi_enreg%comm_cell
 nproc_distrb=xmpi_comm_size(spaceComm_distrb)
 me_distrb=mpi_enreg%me_kpt
 me_g0 = mpi_enreg%me_g0

!if PAW, one has to solve a generalized eigenproblem (H|Psi>=Lambda.S|Psi>)
!else,   one has to solve a classical eigenproblem   (H|Psi>=Lambda.|Psi>)
 gen_eigenpb=(gs_hamk%usepaw==1)
 useoverlap=0;if (gen_eigenpb) useoverlap=1
 !LTEST
 useoverlap=0
 !LTEST

!Initializations and allocations
 isubh=1;isubo=1
 nblock=(nband-1)/nbdblock+1
 istwf_k=gs_hamk%istwf_k
 wfopta10=mod(wfoptalg,10)
 if (wfopta10/=0) then
   MSG_ERROR('cgwf new is implemented only for wfopta10==0')
 end if
 !LTEST
 call writeout(999,'wfoptalg',wfoptalg)
 call writeout(999,'wfopta10',wfopta10)
 call writeout(999,'ortalg',ortalg)
 xnorm=0.0_DP
 !LTEST
 optekin=0;if (wfoptalg>=10) optekin=1
 natom=gs_hamk%natom
 cpopt=-1
 kinpw => gs_hamk%kinpw_k

 ABI_ALLOCATE(pcon,(npw))
 ABI_ALLOCATE(ghc,(2,npw*nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw*nspinor))
 ABI_ALLOCATE(conjgr,(2,npw*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw*nspinor))
 ABI_ALLOCATE(direc,(2,npw*nspinor))
 ABI_ALLOCATE(scprod,(2,nband))

 ABI_ALLOCATE(gh_direc,(2,npw*nspinor))
 ABI_ALLOCATE(gs_direc,(2,0))
 ABI_ALLOCATE(gvnlx_direc,(2,npw*nspinor))
 ABI_ALLOCATE(vresid,(2,npw*nspinor))

 if (gen_eigenpb) then
   ABI_ALLOCATE(direc_tmp,(2,npw*nspinor))
 end if

 ABI_ALLOCATE(swork,(0,0))
 ABI_ALLOCATE(gvnlx_dummy,(0,0))

 ! Loop over blocks of bands. In the standard band-sequential algorithm, nblock=nband.
 do iblock=1,nblock
   counter=100*iblock*nbdblock+inonsc

   ! Loop over bands in a block
   ! This loop can be MPI-parallelized, over processors attached to the same k point
   ibandmin=1+(iblock-1)*nbdblock
   ibandmax=min(iblock*nbdblock,nband)

   ! Big iband loop
   do iband=ibandmin,ibandmax
     ibdblock=iband-(iblock-1)*nbdblock
     counter=100*iband+inonsc
     icg_shift=npw*nspinor*(iband-1)+icg
     igsc_shift=npw*nspinor*(iband-1)+igsc

     ! ======================================================================
     ! ========== INITIALISATION OF MINIMIZATION ITERATIONS =================
     ! ======================================================================

     if (prtvol>=10) then ! Tell us what is going on:
       write(message, '(a,i6,2x,a,i3,a)' )' --- cgwf is called for band',iband,'for',nline,' lines'
       call wrtout(std_out,message,'PERS')
     end if

     dotgp=one

     ! Extraction of the vector that is iteratively updated
     call cg_zcopy(npw*nspinor,cg(1,1+icg_shift),cwavef)

     ! Normalize incoming wf (and S.wf, if generalized eigenproblem):
     ! WARNING : It might be interesting to skip the following operation.
     ! The associated routines should be reexamined to see whether cwavef is not already normalized.
!LTEST
!     if (gen_eigenpb) then
!       call dotprod_g(dotr,doti,istwf_k,npw*nspinor,2,cwavef,scwavef,me_g0,mpi_enreg%comm_spinorfft)
!       dotr=sqrt(dotr**2+doti**2); xnorm=one/sqrt(dotr)
!       call cg_zscal(npw*nspinor,(/xnorm,zero/),cwavef)
!       call cg_zscal(npw*nspinor,(/xnorm,zero/),scwavef)
!     else
!       call sqnorm_g(dotr,istwf_k,npw*nspinor,cwavef,me_g0,mpi_enreg%comm_fft)
!       xnorm=one/sqrt(abs(dotr))
!       call cg_zscal(npw*nspinor,(/xnorm,zero/),cwavef)
!     end if
!
!     if (prtvol==-level) then
!       write(message,'(a,f14.6)')' cgwf: xnorm = ',xnorm
!       call wrtout(std_out,message,'PERS')
!     end if
!
     xnorm = one
!LTEST

     ! ======================================================================
     ! ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
     ! ======================================================================
     if(nline/=0)then
       do iline=1,nline

         ! === COMPUTE THE RESIDUAL ===

         ! Compute lambda = <C|H|C> or <C|(H-zshift)**2|C>
         sij_opt = 0
         call getchc(chc,doti,cpopt,cwavef,cwavef,cprj_dum,gs_hamk,&
&         eval,mpi_enreg,1,npw,nspinor,prtvol,sij_opt,tim_getghc,0)
         lam0=chc

         ! Check that lam0 is decreasing on succeeding lines:
         if (iline==1) then
           lamold=lam0
         else
           if (lam0 > lamold+tol12) then
             write(message, '(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&             'New trial energy at line ',iline,' = ',lam0,ch10,&
&             'is higher than former =',lamold,ch10
             MSG_WARNING(message)
           end if
           lamold=lam0
         end if

         ! Compute residual vector:
         ! Note that vresid is precomputed to garantee cancellation of errors
         ! and allow residuals to reach values as small as 1.0d-24 or better.

         eval=chc
         !TO DO : to compute ( H - epsilon S ) with getghc
         vresid(:,:) = zero
         sij_opt = 0
         if (gen_eigenpb) sij_opt = -1
         call getghc(cpopt,cwavef,cprj_dum,ghc,scwavef,gs_hamk,gvnlxc,&
&         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)

         ! Compute residual (squared) norm
         call sqnorm_g(resid(iband),istwf_k,npw*nspinor,vresid,me_g0,mpi_enreg%comm_fft)

         if (prtvol==-level) then
           write(message,'(a,i0,2f14.6)')' cgwf: iline,eval,resid = ',iline,eval,resid(iband)
           call wrtout(std_out,message,'PERS')
         end if

         ! ======================================================================
         ! ============== CHECK FOR CONVERGENCE CRITERIA ========================
         ! ======================================================================

         ! If residual sufficiently small stop line minimizations
         if (resid(iband)<tolwfr) then
           if (prtvol>=10) then
             write(message, '(a,i4,a,i2,a,es12.4)' ) &
&             ' cgwf: band ',iband,' converged after ',iline,' line minimizations: resid =',resid(iband)
             call wrtout(std_out,message,'PERS')
           end if
           nskip=nskip+(nline-iline+1)  ! Number of two-way 3D ffts skipped
           exit                         ! Exit from the loop on iline
         end if

         ! If user require exiting the job, stop line minimisations
         if (quit==1) then
           write(message, '(a,i0)' )' cgwf: user require exiting => skip update of band ',iband
           call wrtout(std_out,message,'PERS')

           nskip=nskip+(nline-iline+1)  ! Number of two-way 3D ffts skipped
           exit                         ! Exit from the loop on iline
         end if

         ! ======================================================================
         ! =========== COMPUTE THE STEEPEST DESCENT DIRECTION ===================
         ! ======================================================================

         ! Compute the steepest descent direction
         if (gen_eigenpb) then
           call cg_zcopy(npw*nspinor,vresid,direc)  ! Store <G|H-lambda.S|C> in direc
         else
           call cg_zcopy(npw*nspinor,ghc,direc)     ! Store <G|H|C> in direc
         end if

         ! =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
         ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================

         ! The following projection over the subspace orthogonal to occupied bands
         ! is optional. It is a bit more accurate, but doubles the number of N^3 ops.
         ! It is done only if ortalg>=0.

         ! Project the steepest descent direction:
         ! direc(2,npw)=<G|H|Cnk> - \sum_{(i<=n)} <G|H|Cik> , normalized.

         if(ortalg>=0)then
           if (gen_eigenpb) then
             call projbd(cg,direc,iband,icg,igsc,istwf_k,mcg,mgsc,nband,npw,nspinor,&
&             gsc,scprod,0,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
           else
             call projbd(cg,direc,-1   ,icg,igsc,istwf_k,mcg,mgsc,nband,npw,nspinor,&
&             gsc,scprod,0,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
           end if
         else
!          For negative ortalg must still project current band out of conjugate vector (unneeded if gen_eigenpb)
           if (.not.gen_eigenpb) then
             call dotprod_g(dotr,doti,istwf_k,npw*nspinor,3,cwavef,direc,me_g0,mpi_enreg%comm_spinorfft)
             if(istwf_k==1)then
               call cg_zaxpy(npw*nspinor,-(/dotr,doti/),cwavef,direc)
             else
               call cg_zaxpy(npw*nspinor,(/-dotr,zero/),cwavef,direc)
             end if
           end if
         end if

         ! For a generalized eigenpb, store the steepest descent direction
         if (gen_eigenpb) direc_tmp=direc

         ! ======================================================================
         ! ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
         ! ======================================================================

         ! If wfoptalg>=10, the precondition matrix is kept constant during iteration ; otherwise it is recomputed
         if (wfoptalg<10.or.iline==1) then
           call cg_precon(cwavef,zero,istwf_k,kinpw,npw,nspinor,me_g0,optekin,pcon,direc,mpi_enreg%comm_fft)
         else
           do ispinor=1,nspinor
             igs=(ispinor-1)*npw
!$OMP PARALLEL DO
             do ipw=1+igs,npw+igs
               direc(1,ipw)=direc(1,ipw)*pcon(ipw-igs)
               direc(2,ipw)=direc(2,ipw)*pcon(ipw-igs)
             end do
           end do
         end if

         ! ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
         ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
         ! Projecting again out all bands (not normalized).
         call projbd(cg,direc,-1,icg,igsc,istwf_k,mcg,mgsc,nband,npw,nspinor,&
&         gsc,scprod,0,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)

         ! ======================================================================
         ! ================= COMPUTE THE CONJUGATE-GRADIENT =====================
         ! ======================================================================

         if (gen_eigenpb) then
           call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,direc_tmp,me_g0,mpi_enreg%comm_spinorfft)
         else
           call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,ghc,me_g0,mpi_enreg%comm_spinorfft)
         end if

         ! MJV: added 5 Feb 2012 - causes divide by 0 on next iteration of iline
         if (abs(dotgg) < TINY(0.0_dp)*1.e50_dp) dotgg = TINY(0.0_dp)*1.e50_dp

         ! At first iteration, gamma is set to zero
         if (iline==1) then
           gamma=zero
           dotgp=dotgg
           call cg_zcopy(npw*nspinor,direc,conjgr)

         else
           gamma=dotgg/dotgp
           dotgp=dotgg

           if (prtvol==-level)then
             write(message,'(a,2es16.6)')' cgwf: dotgg,gamma = ',dotgg,gamma
             call wrtout(std_out,message,'PERS')
           end if

           ! Note: another way to compute gamma: Polak, Ribiere no real improvement ; to be more carrefully tested
           ! call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direc,direc_tmp)
           ! !direcp must be set to zero at the beginning
           ! direcp=direc-direcp
           ! call dotprod_g(dotgmg,doti,istwf_k,mpi_enreg,npw*nspinor,1,direcp,direc_tmp)
           ! direcp=direc;gamma=dotgmg/dotgp;dotgp=dotgmg

!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             conjgr(1,ipw)=direc(1,ipw)+gamma*conjgr(1,ipw)
             conjgr(2,ipw)=direc(2,ipw)+gamma*conjgr(2,ipw)
           end do
           !call cg_zaxpby(npw*nspinor,cg_one,direc,(/gamma,zero/),conjgr)
         end if

         ! ======================================================================
         ! ============ PROJECTION OF THE CONJUGATED GRADIENT ===================
         ! ======================================================================

         call dotprod_g(dotr,doti,istwf_k,npw*nspinor,3,cwavef,conjgr,me_g0,mpi_enreg%comm_spinorfft)

         ! Project the conjugated gradient onto the current band
         ! MG: TODO: this is an hot spot that could be rewritten with BLAS! provided
         ! that direc --> conjgr
         if(istwf_k==1)then

!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             direc(1,ipw)=conjgr(1,ipw)-(dotr*cwavef(1,ipw)-doti*cwavef(2,ipw))
             direc(2,ipw)=conjgr(2,ipw)-(dotr*cwavef(2,ipw)+doti*cwavef(1,ipw))
           end do
         else
!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             direc(1,ipw)=conjgr(1,ipw)-dotr*cwavef(1,ipw)
             direc(2,ipw)=conjgr(2,ipw)-dotr*cwavef(2,ipw)
           end do
         end if

         ! In case of generalized eigenproblem, normalization of direction vector
         ! cannot be done here (because S|D> is not known here).
         if (.not.gen_eigenpb) then
           call sqnorm_g(dotr,istwf_k,npw*nspinor,direc,me_g0,mpi_enreg%comm_fft)
           xnorm=one/sqrt(abs(dotr))
           call cg_zscal(npw*nspinor,(/xnorm,zero/),direc)
           xnorm=one
         end if

         ! ======================================================================
         ! ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
         ! ======================================================================

         ! Compute gh_direc = <G|H|D> and eventually gs_direc = <G|S|D>
         sij_opt=0;if (gen_eigenpb) sij_opt=1
         !LTEST
         sij_opt=0
         !LTEST

         call getghc(cpopt,direc,cprj_dum,gh_direc,gs_direc,gs_hamk,gvnlx_direc,&
&         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)

         ! In case of generalized eigenproblem, compute now the norm of the conjugated gradient
         !LTEST
!         if (gen_eigenpb) then
!           call dotprod_g(dotr,doti,istwf_k,npw*nspinor,1,direc,gs_direc,me_g0,mpi_enreg%comm_spinorfft)
!           xnorm=one/sqrt(abs(dotr))
!         end i
         !LTEST

         ! Compute dhc = Re{<D|H|C>}
         call dotprod_g(dhc,doti,istwf_k,npw*nspinor,1,direc,ghc,me_g0,mpi_enreg%comm_spinorfft)
         dhc=dhc*xnorm

         ! Compute <D|H|D> or <D|(H-zshift)^2|D>
         call dotprod_g(dhd,doti,istwf_k,npw*nspinor,1,direc,gh_direc,me_g0,mpi_enreg%comm_spinorfft)
         dhd=dhd*xnorm**2

         if(prtvol==-level)then
           write(message,'(a,3f14.6)') 'cgwf: chc,dhc,dhd=',chc,dhc,dhd
           call wrtout(std_out,message,'PERS')
         end if

         ! ======================================================================
         ! ======= COMPUTE MIXING FACTORS - CHECK FOR CONVERGENCE ===============
         ! ======================================================================

         ! Compute tan(2 theta),sin(theta) and cos(theta)
         tan2th=2.0_dp*dhc/(chc-dhd)

         if (abs(tan2th)<1.d-05) then
           costh=1.0_dp-0.125_dp*tan2th**2
           sinth=0.5_dp*tan2th*(1.0_dp-0.375_dp*tan2th**2)

           ! Check that result is above machine precision
           if (abs(sinth)<epsilon(0._dp)) then
             write(message, '(a,es16.4)' ) ' cgwf: converged with tan2th=',tan2th
             call wrtout(std_out,message,'PERS')
             ! Number of one-way 3D ffts skipped
             nskip=nskip+2*(nline-iline)
             exit ! Exit from the loop on iline
           end if

         else
           root=sqrt(1.0_dp+tan2th**2)
           costh=sqrt(0.5_dp+0.5_dp/root)
           sinth=sign(sqrt(0.5_dp-0.5_dp/root),tan2th)
         end if

         ! Check for lower of two possible roots (same sign as curvature at theta where slope is zero)
         diff=(chc-dhd)
         ! Swap c and d if value of diff is positive
         if (diff>zero) then
           swap=costh
           costh=-sinth
           sinth=swap
           if(prtvol<0 .or. prtvol>=10)then
             write(message,*)'   Note: swap roots, iline,diff=',iline,diff
             call wrtout(std_out,message,'PERS')
           end if
         end if

         ! ======================================================================
         ! =========== GENERATE NEW |wf>, H|wf>, Vnl|Wf>, S|Wf> ... =============
         ! ======================================================================

         sintn=sinth*xnorm

!$OMP PARALLEL DO
         do ipw=1,npw*nspinor
           cwavef(1,ipw)=cwavef(1,ipw)*costh+direc(1,ipw)*sintn
           cwavef(2,ipw)=cwavef(2,ipw)*costh+direc(2,ipw)*sintn
         end do

!        call cg_zaxpby(npw*nspinor,(/sintn,zero/),direc,(/costh,zero/),cwavef)
         call cg_zcopy(npw*nspinor,cwavef,cg(1,1+icg_shift))

         do ipw=1,npw*nspinor
           ghc(1,ipw)  =ghc(1,ipw)*costh + gh_direc(1,ipw)*sintn
           ghc(2,ipw)  =ghc(2,ipw)*costh + gh_direc(2,ipw)*sintn
         end do


         if (use_subvnlx==1) then
!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             gvnlxc(1,ipw)=gvnlxc(1,ipw)*costh + gvnlx_direc(1,ipw)*sintn
             gvnlxc(2,ipw)=gvnlxc(2,ipw)*costh + gvnlx_direc(2,ipw)*sintn
           end do
!          call cg_zaxpby(npw*nspinor,(/sintn,zero/),gvnlx_direc,(/costh,zero/),gvnlxc)
         end if

         ! ======================================================================
         ! =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
         ! ======================================================================

         ! Compute delta(E)
         deltae=chc*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc

!        Check convergence and eventually exit
         if (iline==1) then
           deold=deltae
         else if (abs(deltae)<tolrde*abs(deold) .and. iline/=nline .and. wfopta10<2)then
           if(prtvol>=10)then
             write(message, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
&             ' cgwf: line',iline,&
&             ' deltae=',deltae,' < tolrde*',deold,' =>skip lines'
             call wrtout(std_out,message,'PERS')
           end if
           nskip=nskip+2*(nline-iline)  ! Number of one-way 3D ffts skipped
           exit                         ! Exit from the loop on iline
         end if

       end do ! END LOOP FOR A GIVEN BAND Note that there are three "exit" instructions inside

     else ! nline==0 , needs to provide a residual
       resid(iband)=-one
     end if ! End nline==0 case

     ! ======================================================================
     ! =============== END OF CURRENT BAND: CLEANING ========================
     ! ======================================================================

     ! At the end of the treatment of a set of bands, write the number of one-way 3D ffts skipped
     if (xmpi_paral==0 .and. mpi_enreg%paral_kgb==0 .and. iband==nband .and. prtvol/=0) then
       write(message,'(a,i0)')' cgwf: number of one-way 3D ffts skipped in cgwf until now =',nskip
       call wrtout(std_out,message,'PERS')
     end if

   end do !  End big iband loop. iband in a block

   !  ======================================================================
   !  ============= COMPUTE HAMILTONIAN IN WFs SUBSPACE ====================
   !  ======================================================================
   !LTEST
!   call writeout(999,'iblock',iblock)
!   call writeout(999,'ghc-cwavef',ghc-cwavef)
!   call writeout(999,'gsc-cwavef',scwavef-cwavef)
   !LTEST
   call mksubham(cg,ghc,gsc,gvnlxc,iblock,icg,igsc,istwf_k,&
&   isubh,isubo,mcg,mgsc,nband,nbdblock,npw,&
&   nspinor,subham,subovl,subvnlx,use_subovl,use_subvnlx,me_g0)
   !LTEST
   call sqnorm_g(dotr,istwf_k,npw*nspinor,cwavef,me_g0,mpi_enreg%comm_fft)
   xnorm=xnorm+dotr
!   call writeout(999,'subham',subham(isubh-2)-dotr)
   !LTEST

 end do ! iblock End loop over block of bands
 !LTEST
! call writeout(999,'xnorm',xnorm)
 !LTEST

 if (allocated(dimlmn_srt)) then
   ABI_DEALLOCATE(dimlmn_srt)
 end if

 ! Debugging ouputs
 if(prtvol==-level)then
   isubh=1
   if (use_subvnlx==1) write(message,'(a)') ' cgwf : isubh  subham(isubh:isubh+1)  subvnlx(isubh:isubh+1)'
   if (use_subvnlx==0) write(message,'(a)') ' cgwf : isubh  subham(isubh:isubh+1)'
   do iband=1,nband
     do ii=1,iband
       if (use_subvnlx==1) then
         write(message,'(i5,4es16.6)')isubh,subham(isubh:isubh+1),subvnlx(isubh:isubh+1)
       else
         write(message,'(i5,2es16.6)')isubh,subham(isubh:isubh+1)
       end if
       call wrtout(std_out,message,'PERS')
       isubh=isubh+2
     end do
   end do
 end if

 ! ===================
 ! FINAL DEALLOCATIONS
 ! ===================
 ABI_DEALLOCATE(conjgr)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(direc)
 ABI_DEALLOCATE(pcon)
 ABI_DEALLOCATE(scprod)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(gh_direc)
 ABI_DEALLOCATE(gvnlx_direc)
 ABI_DEALLOCATE(vresid)
 ABI_DEALLOCATE(gs_direc)

 if (gen_eigenpb)  then
   ABI_DEALLOCATE(direc_tmp)
 end if

 ABI_DEALLOCATE(gvnlx_dummy)

 ABI_DEALLOCATE(swork)

! Do not delete this line, needed to run with open MP
 write(unit=message,fmt=*) resid(1)

 call timab(22,2,tsec)

 DBG_EXIT("COLL")

end subroutine cgwfnew
!!***

end module m_cgwfnew
!!***
