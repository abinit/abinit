!!****m* ABINIT/m_cgwf_paw
!! NAME
!!  m_cgwf_paw
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

module m_cgwf_paw

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_cgtools

 use defs_abitypes,   only : MPI_type
 use m_time,          only : timab
 use m_numeric_tools, only : rhophi
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_put, pawcprj_copy, &
                             pawcprj_get, pawcprj_mpi_allgather, pawcprj_free, pawcprj_symkn
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getchc,        only : getchc,getcsc
 use m_getghc,        only : getghc
 use m_nonlop,        only : nonlop
 use m_paw_overlap,   only : smatrix_k_paw
 use m_cgprj,         only : getcprj,cprj_axpby
 use m_fft,           only : fourwf

 implicit none

 private
!!***

 public :: cgwf_paw

!!***

contains
!!***

!!****f* m_cgwf/cgwf_paw
!! NAME
!! cgwf_paw
!!
!! FUNCTION
!! Update all wavefunction |C>, non self-consistently.
!! also compute the corresponding H|C> and Vnl|C> (and S|C> if paw).
!! Uses a conjugate-gradient algorithm.
!! In case of paw, resolves a generalized eigenproblem using an
!!  overlap matrix (not used for norm conserving psps).
!!
!! INPUTS
!!  chkexit= if non-zero, check whether the user wishes to exit
!!  cpus = CPU time limit
!!  filnam_ds1=name of input file (used for exit checking)
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  inonsc=index of non self-consistent loop
!!  isppol=spin polarization currently treated
!!  mcg=second dimension of the cg array
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands.
!!  nbdblock=number of bands in a block
!!  nkpt=number of k points
!!  nline=number of line minimizations per band.
!!  npw=number of planewaves in basis sphere at given k.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  prtvol=control print volume and debugging output
!!  tolrde=tolerance on the ratio of differences of energies (for the line minimisation)
!!  tolwfr=tolerance on largest wf residual
!!  wfoptalg=govern the choice of algorithm for wf optimisation
!!   (0, 1, 10 and 11 : in the present routine, usual CG algorithm ;
!!   (2 and 3 : use shifted square Hamiltonian)
!!
!! OUTPUT
!!  resid(nband)=wf residual for new states=|(H-e)|C>|^2 (hartree^2)
!!
!! SIDE EFFECTS
!!  cg(2,mcg)
!!    at input =wavefunction <G|C band,k> coefficients for ALL bands
!!    at output same as input except that
!!      the current band, with number 'band' has been updated
!!  quit= if 1, proceeds to smooth ending of the job.
!!
!! NOTES
!!  1) cg should not be filtered and normalized : it should already be OK at input !
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

subroutine cgwf_paw(cg,chkexit,cpus,eig,&
&                filnam_ds1,gs_hamk,icg,ikpt,inonsc,&
&                isppol,mcg,mpi_enreg,&
&                nband,nbdblock,nkpt,nline,npw,&
&                nspinor,nsppol,ortalg,prtvol,quit,resid,subham,&
&                tolrde,tolwfr,wfoptalg)
!Arguments ------------------------------------
 integer,intent(in) :: chkexit,icg,ikpt,inonsc,isppol
 integer,intent(in) :: mcg,nband,nbdblock,nkpt,nline
 integer,intent(in) :: npw,nspinor,nsppol,ortalg,prtvol
 integer,intent(in) :: wfoptalg
 integer,intent(in) :: quit
 real(dp),intent(in) :: cpus,tolrde,tolwfr
 character(len=*),intent(in) :: filnam_ds1
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 real(dp),intent(inout),target :: cg(2,mcg)
 real(dp), intent(inout) :: eig(nband)
 real(dp),intent(out) :: resid(nband),subham(nband*(nband+1))

!Local variables-------------------------------
integer,parameter :: level=113,tim_getghc=1,tim_projbd=1,type_calc=0
integer,parameter :: tim_getchc=0,tim_getcsc=3,tim_getcsc_band=4,tim_fourwf=40
 integer,save :: nskip=0
 integer :: counter,cpopt,itypat
 integer :: i1,i2,i3,ia,iband,ibandmin,ibandmax,isubh,jband,me_g0
 integer :: ibdblock,iblock,igs
 integer :: iline,ipw,ispinor,istwf_k
 integer :: n4,n5,n6,natom,ncpgr,nblock
 integer :: optekin,sij_opt
 integer :: useoverlap,wfopta10
 real(dp) :: chc,costh,deltae,deold,dhc,dhd,diff,dotgg,dotgp,doti,dotr,eval,gamma
 real(dp) :: lam0,lamold,root,sinth,sintn,swap,tan2th,theta,weight_fft,xnorm
 real(dp) :: dot(2)
 character(len=500) :: message
 integer,allocatable :: dimlmn(:)
 real(dp) :: tsec(2),dotrr(1),dotii(1)
 real(dp),allocatable :: conjgr(:,:),gvnlxc(:,:),denpot_dum(:,:,:),fofgout_dum(:,:)
 real(dp), pointer :: cwavef(:,:),cwavef_left(:,:),cwavef_bands(:,:)
 real(dp), pointer :: cwavef_r(:,:,:,:),cwavef_r_left(:,:,:,:)
 real(dp),allocatable,target :: cwavef_r_bands(:,:,:,:,:)
 real(dp),allocatable :: direc(:,:),direc_tmp(:,:),pcon(:),scprod(:,:),scwavef_dum(:,:),direc_r(:,:,:,:)
 real(dp),pointer :: scprod_csc(:),kinpw(:)
 real(dp) :: z_tmp(2),z_tmp2(2)
 type(pawcprj_type),allocatable,target :: cprj_cwavef_bands(:,:)
 type(pawcprj_type),pointer :: cprj_cwavef(:,:),cprj_cwavef_left(:,:)
 type(pawcprj_type),allocatable :: cprj_direc(:,:),cprj_conjgr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Starting the routine
 call timab(1300,1,tsec)

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
 me_g0 = mpi_enreg%me_g0

 useoverlap=0

 cpopt = 2

 weight_fft = one

!Initializations and allocations
 nblock=(nband-1)/nbdblock+1
 if(nbdblock>1)then
   MSG_BUG('cgwf_paw not implemented for nbdblock>1!')
 end if

 istwf_k=gs_hamk%istwf_k
 wfopta10=mod(wfoptalg,10)
 if (wfopta10/=0) then
   MSG_BUG('cgwf_paw is implemented only for wfopta10==0')
 end if
 if (ortalg>=0) then
   MSG_BUG('cgwf_paw tested only for ortalg<0')
 end if
 optekin=0;if (wfoptalg>=10) optekin=1
 natom=gs_hamk%natom
 kinpw => gs_hamk%kinpw_k

 ABI_ALLOCATE(pcon,(npw))
 ABI_ALLOCATE(conjgr,(2,npw*nspinor))
 ABI_ALLOCATE(scwavef_dum,(0,0))
 ABI_ALLOCATE(direc,(2,npw*nspinor))
 ABI_ALLOCATE(scprod,(2,nband))
 ABI_ALLOCATE(scprod_csc,(2*nband))
 ABI_ALLOCATE(direc_tmp,(2,npw*nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw*nspinor))

 ABI_DATATYPE_ALLOCATE(cprj_cwavef_bands,(natom,nband))
 ABI_DATATYPE_ALLOCATE(cprj_direc ,(natom,nbdblock))
 ABI_DATATYPE_ALLOCATE(cprj_conjgr ,(natom,nbdblock))

 n4=gs_hamk%ngfft(4);n5=gs_hamk%ngfft(5);n6=gs_hamk%ngfft(6)
 ABI_ALLOCATE(cwavef_r_bands,(2,n4,n5,n6,nband))
 ABI_ALLOCATE(direc_r, (2,n4,n5,n6))
 ABI_ALLOCATE(denpot_dum, (0,0,0))
 ABI_ALLOCATE(fofgout_dum, (0,0))

 ncpgr = 0 ! no need of gradients here
!Dimensioning and allocation of <p_i|Cnk>
 ABI_ALLOCATE(dimlmn,(natom))
 dimlmn=0  ! Type-sorted cprj
 ia=0
 do itypat=1,gs_hamk%ntypat
   dimlmn(ia+1:ia+gs_hamk%nattyp(itypat))=count(gs_hamk%indlmn(3,:,itypat)>0)
   ia=ia+gs_hamk%nattyp(itypat)
 end do
 call pawcprj_alloc(cprj_cwavef_bands,ncpgr,dimlmn)
 call pawcprj_alloc(cprj_direc,ncpgr,dimlmn)
 call pawcprj_alloc(cprj_conjgr,ncpgr,dimlmn)

 cwavef_bands => cg(:,1+icg:nblock*npw*nspinor+icg)

 do iblock=1,nblock
   ibandmin=1+(iblock-1)*nbdblock
   ibandmax=min(iblock*nbdblock,nband)
   do iband=ibandmin,ibandmax
     ibdblock=iband-(iblock-1)*nbdblock

     cwavef => cwavef_bands(:,1+(iband-1)*npw*nspinor:iband*npw*nspinor)
     cprj_cwavef => cprj_cwavef_bands(:,iband:iband)
     cwavef_r => cwavef_r_bands(:,:,:,:,iband)

     call getcprj(1,0,cwavef,cprj_cwavef,&
&      gs_hamk%ffnl_k,0,gs_hamk%indlmn,istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&      gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,natom,gs_hamk%nattyp,&
&      gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&      gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)

     ! Compute wavefunction in real space
     call fourwf(0,denpot_dum,cwavef,fofgout_dum,cwavef_r,gs_hamk%gbound_k,gs_hamk%gbound_k,istwf_k,&
&      gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,gs_hamk%npw_fft_k,gs_hamk%npw_fft_k,&
&      n4,n5,n6,0,tim_fourwf,weight_fft,weight_fft)

   end do
 end do

 isubh = 1

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

     ! ======================================================================
     ! ========== INITIALISATION OF MINIMIZATION ITERATIONS =================
     ! ======================================================================

     if (prtvol>=10) then ! Tell us what is going on:
       write(message, '(a,i6,2x,a,i3,a)' )' --- cgwf is called for band',iband,'for',nline,' lines'
       call wrtout(std_out,message,'PERS')
     end if

     dotgp=one

     ! Extraction of the vector that is iteratively updated
     cwavef => cwavef_bands(:,1+(iband-1)*npw*nspinor:iband*npw*nspinor)
     cprj_cwavef => cprj_cwavef_bands(:,iband:iband)
     cwavef_r => cwavef_r_bands(:,:,:,:,iband)

     ! Normalize incoming wf (and S.wf, if generalized eigenproblem):
     ! WARNING : It might be interesting to skip the following operation.
     ! The associated routines should be reexamined to see whether cwavef is not already normalized.
     call getcsc(dot,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,&
&     gs_hamk,mpi_enreg,1,prtvol,tim_getcsc)
     xnorm=one/sqrt(dot(1))
     z_tmp = (/xnorm,zero/)
     call cg_zscal(npw*nspinor,z_tmp,cwavef)
     z_tmp2 = (/zero,zero/)
     call cprj_axpby(cprj_cwavef,cprj_cwavef,cprj_cwavef,z_tmp,z_tmp2,&
&             gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&             natom,gs_hamk%nattyp,1,nspinor,gs_hamk%ntypat)
     cwavef_r=cwavef_r*xnorm

     if (prtvol==-level) then
       write(message,'(a,f14.6)')' cgwf: xnorm = ',xnorm
       call wrtout(std_out,message,'PERS')
     end if

     ! ======================================================================
     ! ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
     ! ======================================================================
     if(nline/=0)then
       do iline=1,nline

         ! === COMPUTE THE RESIDUAL ===
         ! Compute lambda = <C|H|C>
         sij_opt = 0
         call getchc(chc,doti,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
           &          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
         lam0=chc
         eval=chc
         eig(iband)=chc

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

         ! ======================================================================
         ! =========== COMPUTE THE STEEPEST DESCENT DIRECTION ===================
         ! ======================================================================

         sij_opt = -1
         call getghc(cpopt,cwavef,cprj_cwavef,direc,scwavef_dum,gs_hamk,gvnlxc,&
           &         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,type_calc)

         ! Compute residual (squared) norm
         call sqnorm_g(resid(iband),istwf_k,npw*nspinor,direc,me_g0,mpi_enreg%comm_fft)

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

         ! =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
         ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================

         ! The following projection over the subspace orthogonal to occupied bands
         ! is optional. It is a bit more accurate, but doubles the number of N^3 ops.
         ! It is done only if ortalg>=0.

         ! Project the steepest descent direction:
         ! direc(2,npw)=<G|H|Cnk> - \sum_{(i<=n)} <G|H|Cik> , normalized.
         if(ortalg>=0)then
           call getcprj(1,0,direc,cprj_direc,&
&           gs_hamk%ffnl_k,0,gs_hamk%indlmn,istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&           gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,natom,gs_hamk%nattyp,&
&           gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&           gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
           call getcsc(scprod_csc,cpopt,direc,cwavef_bands,cprj_direc,cprj_cwavef_bands,&
&           gs_hamk,mpi_enreg,nband,prtvol,tim_getcsc_band)
           scprod = reshape(scprod_csc,(/2,nband/))
           call projbd(cg,direc,iband,icg,icg,istwf_k,mcg,mcg,nband,npw,nspinor,&
&           direc,scprod,1,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
         end if

         ! For a generalized eigenpb, store the steepest descent direction
         direc_tmp=direc

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
         call getcprj(1,0,direc,cprj_direc,&
&         gs_hamk%ffnl_k,0,gs_hamk%indlmn,istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&         gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,natom,gs_hamk%nattyp,&
&         gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&         gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
         call getcsc(scprod_csc,cpopt,direc,cwavef_bands,cprj_direc,cprj_cwavef_bands,&
&         gs_hamk,mpi_enreg,nband,prtvol,tim_getcsc_band)
         ! Projecting again out all bands (not normalized).
         scprod = reshape(scprod_csc,(/2,nband/))
         call projbd(cg,direc,-1,icg,icg,istwf_k,mcg,mcg,nband,npw,nspinor,&
&         direc,scprod,1,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
         ! Apply projbd to cprj_direc
         z_tmp  = (/one,zero/)
         scprod_csc = -scprod_csc
         call cprj_axpby(cprj_direc,cprj_direc,cprj_cwavef_bands,z_tmp,scprod_csc,&
&                 gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&                 natom,gs_hamk%nattyp,nband,nspinor,gs_hamk%ntypat)

         ! ======================================================================
         ! ================= COMPUTE THE CONJUGATE-GRADIENT =====================
         ! ======================================================================

         call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,direc_tmp,me_g0,mpi_enreg%comm_spinorfft)

         ! MJV: added 5 Feb 2012 - causes divide by 0 on next iteration of iline
         if (abs(dotgg) < TINY(0.0_dp)*1.e50_dp) dotgg = TINY(0.0_dp)*1.e50_dp

         ! At first iteration, gamma is set to zero
         if (iline==1) then
           gamma=zero
           dotgp=dotgg
           call cg_zcopy(npw*nspinor,direc,conjgr)
           z_tmp  = (/one,zero/)
           z_tmp2 = (/zero,zero/)
           call cprj_axpby(cprj_conjgr,cprj_direc,cprj_direc,z_tmp,z_tmp2,&
&                   gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&                   natom,gs_hamk%nattyp,1,nspinor,gs_hamk%ntypat)

         else
           gamma=dotgg/dotgp
           dotgp=dotgg

           if (prtvol==-level)then
             write(message,'(a,2es16.6)')' cgwf: dotgg,gamma = ',dotgg,gamma
             call wrtout(std_out,message,'PERS')
           end if

!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             conjgr(1,ipw)=direc(1,ipw)+gamma*conjgr(1,ipw)
             conjgr(2,ipw)=direc(2,ipw)+gamma*conjgr(2,ipw)
           end do
           z_tmp   = (/one,zero/)
           z_tmp2  = (/gamma,zero/)
           call cprj_axpby(cprj_conjgr,cprj_direc,cprj_conjgr,z_tmp,z_tmp2,&
&                   gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&                   natom,gs_hamk%nattyp,1,nspinor,gs_hamk%ntypat)
         end if

         ! ======================================================================
         ! ============ PROJECTION OF THE CONJUGATED GRADIENT ===================
         ! ======================================================================

         call getcsc(dot,cpopt,conjgr,cwavef,cprj_conjgr,cprj_cwavef,&
&         gs_hamk,mpi_enreg,1,prtvol,tim_getcsc)
         dotr=dot(1)

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
         z_tmp   = (/one,zero/)
         z_tmp2  = (/-dotr,zero/)
         call cprj_axpby(cprj_direc,cprj_conjgr,cprj_cwavef,z_tmp,z_tmp2,&
&         gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&         natom,gs_hamk%nattyp,1,nspinor,gs_hamk%ntypat)

         ! ======================================================================
         ! ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
         ! ======================================================================

         ! Compute norm of direc
         call getcsc(dot,cpopt,direc,direc,cprj_direc,cprj_direc,&
&         gs_hamk,mpi_enreg,1,prtvol,tim_getcsc)
         xnorm=one/sqrt(abs(dot(1)))

         sij_opt=0
         ! Compute dhc = Re{<D|H|C>}
         ! Compute direc in real space
         call fourwf(0,denpot_dum,direc,fofgout_dum,direc_r,gs_hamk%gbound_k,gs_hamk%gbound_k,istwf_k,&
&          gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,gs_hamk%npw_fft_k,gs_hamk%npw_fft_k,&
&          n4,n5,n6,0,tim_fourwf,weight_fft,weight_fft)
         call getchc(dhc,doti,cpopt,cwavef,direc,cprj_cwavef,cprj_direc,cwavef_r,direc_r,&
           &          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
         dhc=dhc*xnorm

         ! Compute <D|H|D> or <D|(H-zshift)^2|D>
         call getchc(dhd,doti,cpopt,direc,direc,cprj_direc,cprj_direc,direc_r,direc_r,&
&          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
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
         z_tmp  = (/costh,zero/)
         z_tmp2 = (/sintn,zero/)
         call cprj_axpby(cprj_cwavef,cprj_cwavef,cprj_direc,z_tmp,z_tmp2,&
&         gs_hamk%indlmn,istwf_k,gs_hamk%lmnmax,mpi_enreg,&
&         natom,gs_hamk%nattyp,1,nspinor,gs_hamk%ntypat)
         do i3=1,gs_hamk%n6
           do i2=1,gs_hamk%n5
             do i1=1,gs_hamk%n4
               cwavef_r(1,i1,i2,i3)=cwavef_r(1,i1,i2,i3)*costh+direc_r(1,i1,i2,i3)*sintn
               cwavef_r(2,i1,i2,i3)=cwavef_r(2,i1,i2,i3)*costh+direc_r(2,i1,i2,i3)*sintn
             end do
           end do
         end do

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
           ! Update chc before exit
           call getchc(chc,doti,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
             &          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
           eig(iband)=chc

           nskip=nskip+2*(nline-iline)  ! Number of one-way 3D ffts skipped
           exit                         ! Exit from the loop on iline
         end if

         ! Update chc only if last iteration, otherwise it will be done at the beginning of the next one
         if (iline==nline) then
           call getchc(chc,doti,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
             &          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
           eig(iband)=chc
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

   do iband=ibandmin,ibandmax
     cwavef => cwavef_bands(:,1+(iband-1)*npw*nspinor:iband*npw*nspinor)
     cprj_cwavef => cprj_cwavef_bands(:,iband:iband)
     cwavef_r => cwavef_r_bands(:,:,:,:,iband)
     do jband=1,iband
       cwavef_left => cwavef_bands(:,1+(jband-1)*npw*nspinor:jband*npw*nspinor)
       cprj_cwavef_left => cprj_cwavef_bands(:,jband:jband)
       cwavef_r_left => cwavef_r_bands(:,:,:,:,jband)
       call getchc(subham(isubh),subham(isubh+1),cpopt,cwavef,cwavef_left,&
         &          cprj_cwavef,cprj_cwavef_left,cwavef_r,cwavef_r_left,&
         &          gs_hamk,zero,mpi_enreg,1,prtvol,sij_opt,tim_getchc,type_calc)
       isubh=isubh+2
     end do
   end do

 end do ! iblock End loop over block of bands

 ! Debugging ouputs
 if(prtvol==-level)then
   isubh=1
   do iband=1,nband
     do jband=1,iband
       write(message,'(i5,2es16.6)')isubh,subham(isubh:isubh+1)
       call wrtout(std_out,message,'PERS')
       isubh=isubh+2
     end do
   end do
 end if

 ! ===================
 ! FINAL DEALLOCATIONS
 ! ===================

 nullify(cprj_cwavef)
 call pawcprj_free(cprj_cwavef_bands)
 call pawcprj_free(cprj_direc)
 call pawcprj_free(cprj_conjgr)
 ABI_DATATYPE_DEALLOCATE(cprj_cwavef_bands)
 ABI_DATATYPE_DEALLOCATE(cprj_conjgr)
 ABI_DEALLOCATE(dimlmn)
 ABI_DEALLOCATE(conjgr)
 ABI_DEALLOCATE(scwavef_dum)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(direc)
 ABI_DEALLOCATE(pcon)
 ABI_DEALLOCATE(scprod)

 ABI_DEALLOCATE(direc_tmp)
 ABI_DEALLOCATE(cwavef_r_bands)
 ABI_DEALLOCATE(direc_r)
 ABI_DEALLOCATE(denpot_dum)
 ABI_DEALLOCATE(fofgout_dum)

! Do not delete this line, needed to run with open MP
 write(unit=message,fmt=*) resid(1)

 call timab(1300,2,tsec)

 DBG_EXIT("COLL")

end subroutine cgwf_paw
!!***

end module m_cgwf_paw
!!***
