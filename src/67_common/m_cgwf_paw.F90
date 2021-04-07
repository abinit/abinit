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
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_copy, &
                             pawcprj_free,&
                             pawcprj_axpby,pawcprj_zaxpby,pawcprj_projbd
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getchc,        only : getchc,getcsc
 use m_getghc,        only : getghc
 use m_nonlop,        only : nonlop
 use m_paw_overlap,   only : smatrix_k_paw
 use m_cgprj,         only : getcprj
 use m_fft,           only : fourwf
 use m_dtset,         only : dataset_type

 implicit none

 private
!!***

 public :: cgwf_paw
 public :: mksubovl
 public :: cprj_update
 public :: cprj_update_oneband
 public :: cprj_check
 public :: get_cprj_id
 public :: enable_cgwf_paw

!!***

contains
!!***

!!****f* m_cgwf/cgwf_paw
!! NAME
!! cgwf_paw
!!
!! FUNCTION
!! Update all wavefunction |C>, non self-consistently.
!! Uses a conjugate-gradient algorithm.
!! This version is available for PAW only, as it needs the cprj in memory.
!!
!! INPUTS
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  mcg=second dimension of the cg array
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands.
!!  nline=number of line minimizations per band.
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  prtvol=control print volume and debugging output
!!  tolrde=tolerance on the ratio of differences of energies (for the line minimisation)
!!  tolwfr=tolerance on largest wf residual
!!  wfoptalg=govern the choice of algorithm for wf optimisation
!!   (0, 1, 10 and 11 : in the present routine, usual CG algorithm ;
!!
!! OUTPUT
!!  resid(nband)=wf residual for new states=|(H-e)|C>|^2 (hartree^2)
!!
!! SIDE EFFECTS
!!  cg(2,mcg)
!!    at input =wavefunction <G|C band,k> coefficients for ALL bands
!!    at output same as input except that
!!      the current band, with number 'band' has been updated
!!  cprj_cwavef_bands=<p_i|c_n> coefficients for all bands n, they are updated when the WFs change.
!!  quit= if 1, proceeds to smooth ending of the job.
!!
!! NOTES
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

subroutine cgwf_paw(cg,cprj_cwavef_bands,cprj_update_lvl,eig,&
&                gs_hamk,icg,mcg,mpi_enreg,nband,nline,ortalg,prtvol,quit,resid,subham,&
&                tolrde,tolwfr,wfoptalg)
!Arguments ------------------------------------
 integer,intent(in) :: cprj_update_lvl,icg
 integer,intent(in) :: mcg,nband,nline,ortalg,prtvol
 integer,intent(in) :: wfoptalg,quit
 real(dp),intent(in) :: tolrde,tolwfr
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 real(dp),intent(inout),target :: cg(2,mcg)
 real(dp), intent(inout) :: eig(:)
 real(dp),intent(out) :: resid(:),subham(:)
 type(pawcprj_type),intent(inout),target :: cprj_cwavef_bands(:,:)

!Local variables-------------------------------
integer,parameter :: level=113,tim_getghc=1,tim_projbd=1,type_calc=0
integer,parameter :: useoverlap=0,tim_getcsc=3
 integer,save :: nskip=0
 integer :: cpopt,i1,i2,i3,iband,isubh,isubh0,jband,me_g0,igs
 integer :: iline,ipw,ispinor,istwf_k
 integer :: n4,n5,n6,natom,ncpgr,npw,nspinor
 integer :: optekin,sij_opt,wfopta10
 real(dp) :: chc,costh,deltae,deold,dhc,dhd,diff,dotgg,dotgp,doti,dotr,eval,gamma
 real(dp) :: lam0,lamold,root,sinth,sintn,swap,tan2th,xnorm,xnormd,xnormd_previous
 character(len=500) :: message
!arrays
 real(dp) :: dot(2)
 real(dp) :: tsec(2)
 real(dp),allocatable :: conjgr(:,:),gvnlxc(:,:)
 real(dp), pointer :: cwavef(:,:),cwavef_left(:,:),cwavef_bands(:,:)
 real(dp), allocatable :: cwavef_r(:,:,:,:,:)
 real(dp),allocatable,target :: direc(:,:)
 real(dp),allocatable :: direc_tmp(:,:),pcon(:),scprod(:,:),scwavef_dum(:,:)
 real(dp),allocatable :: direc_r(:,:,:,:,:),scprod_csc(:)
 real(dp) :: z_tmp(2),z_tmp2(2)
 type(pawcprj_type),pointer :: cprj_cwavef(:,:),cprj_cwavef_left(:,:)
 type(pawcprj_type),allocatable :: cprj_direc(:,:),cprj_conjgr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Starting the routine
 call timab(1300,1,tsec)

!Some checks
!Only wfoptalg==10 for now
 if (wfoptalg/=10) then
   ABI_ERROR("cgwf_paw is implemented for wfoptalg==10 only")
 end if

!======================================================================
!========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
!======================================================================

!MPI data
 me_g0 = mpi_enreg%me_g0

! cprj are already in memory
 cpopt = 2

!Initializations and allocations
 istwf_k=gs_hamk%istwf_k
 wfopta10=mod(wfoptalg,10)

 optekin=0;if (wfoptalg>=10) optekin=1
 natom=gs_hamk%natom
 npw=gs_hamk%npw_k
 nspinor=gs_hamk%nspinor

 ABI_MALLOC(pcon,(npw))
 ABI_MALLOC(conjgr,(2,npw*nspinor))
 ABI_MALLOC(scwavef_dum,(0,0))
 ABI_MALLOC(direc,(2,npw*nspinor))
 ABI_MALLOC(scprod,(2,nband))
 ABI_MALLOC(scprod_csc,(2*nband))
 ABI_MALLOC(direc_tmp,(2,npw*nspinor))
 ABI_MALLOC(gvnlxc,(2,npw*nspinor))

 ABI_MALLOC(cprj_direc ,(natom,nspinor))
 ABI_MALLOC(cprj_conjgr ,(natom,nspinor))

 n4=gs_hamk%ngfft(4);n5=gs_hamk%ngfft(5);n6=gs_hamk%ngfft(6)
 ABI_MALLOC(cwavef_r,(2,n4,n5,n6,nspinor))
 ABI_MALLOC(direc_r, (2,n4,n5,n6,nspinor))

 ncpgr  = 0 ! no need of gradients here...
 call pawcprj_alloc(cprj_direc,ncpgr,gs_hamk%dimcprj)
 call pawcprj_alloc(cprj_conjgr,ncpgr,gs_hamk%dimcprj)

 cwavef_bands => cg(:,1+icg:nband*npw*nspinor+icg)

 if (cprj_update_lvl==-1) then
   call timab(1203,1,tsec)
   call cprj_update(cg,cprj_cwavef_bands,gs_hamk,icg,nband,mpi_enreg)
   call timab(1203,2,tsec)
 end if

 ! Big iband loop
 do iband=1,nband

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
   cprj_cwavef => cprj_cwavef_bands(:,nspinor*(iband-1)+1:nspinor*iband)

   ! Normalize incoming wf:
   call getcsc(dot,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,&
&   gs_hamk,mpi_enreg,1,tim_getcsc)
   xnorm=one/sqrt(dot(1))
   z_tmp = (/xnorm,zero/)
   ! cwavef = xnorm * cwavef
   call timab(1305,1,tsec)
   call cg_zscal(npw*nspinor,z_tmp,cwavef)
   call timab(1305,2,tsec)
   ! cprj = xnorm * cprj
   call timab(1302,1,tsec)
   call pawcprj_axpby(zero,xnorm,cprj_cwavef,cprj_cwavef)
   call timab(1302,2,tsec)

   ! Compute wavefunction in real space
   call get_cwavefr(cwavef,cwavef_r,gs_hamk,mpi_enreg)

   if (prtvol==-level) then
     write(message,'(a,f26.14)')' cgwf: xnorm = ',xnorm
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
       call getchc(z_tmp,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
         &          gs_hamk,zero,mpi_enreg,1,sij_opt,type_calc)
       chc=z_tmp(1)
       lam0=chc
       eval=chc
       eig(iband)=chc

       ! Check that lam0 is decreasing on succeeding lines:
       if (iline==1) then
         lamold=lam0
       else
         if (lam0 > lamold+tol12) then
           write(message, '(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&           'New trial energy at line ',iline,' = ',lam0,ch10,&
&           'is higher than former =',lamold,ch10
           ABI_WARNING(message)
         end if
         lamold=lam0
       end if

       ! ======================================================================
       ! =========== COMPUTE THE STEEPEST DESCENT DIRECTION ===================
       ! ======================================================================

       sij_opt = -1
       call getghc(cpopt,cwavef,cprj_cwavef,direc,scwavef_dum,gs_hamk,gvnlxc,&
         &         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,type_calc)

       ! Compute residual (squared) norm
       call timab(1305,1,tsec)
       call sqnorm_g(resid(iband),istwf_k,npw*nspinor,direc,me_g0,mpi_enreg%comm_fft)
       call timab(1305,2,tsec)

       if (prtvol==-level) then
         write(message,'(a,i0,f26.14,es27.14e3)')' cgwf: iline,eval,resid = ',iline,eval,resid(iband)
         call wrtout(std_out,message,'PERS')
       end if
       xnormd=1/sqrt(resid(iband))
       z_tmp = (/xnormd,zero/)
       call cg_zscal(npw*nspinor,z_tmp,direc)

       ! ======================================================================
       ! ============== CHECK FOR CONVERGENCE CRITERIA ========================
       ! ======================================================================

       ! If residual sufficiently small stop line minimizations
       if (resid(iband)<tolwfr) then
         if (prtvol>=10) then
           write(message, '(a,i4,a,i2,a,es12.4)' ) &
&           ' cgwf: band ',iband,' converged after ',iline,' line minimizations: resid =',resid(iband)
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
       ! direc(2,npw)=<G|H|Cnk> - \sum_{(i/=n)} <G|H|Cik>
       if(ortalg>=0)then
         call timab(1203,1,tsec)
         call cprj_update_oneband(direc,cprj_direc,gs_hamk,mpi_enreg)
         call timab(1203,2,tsec)
         call getcsc(scprod_csc,cpopt,direc,cwavef_bands,cprj_direc,cprj_cwavef_bands,&
&         gs_hamk,mpi_enreg,nband,tim_getcsc)
         scprod = reshape(scprod_csc,(/2,nband/))
         ! Note that the current band (|C_iband>) is not used here
         call projbd(cg,direc,iband,icg,icg,istwf_k,mcg,mcg,nband,npw,nspinor,&
&         direc,scprod,1,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
       end if

       ! For a generalized eigenpb, store the steepest descent direction
       direc_tmp=direc

       ! ======================================================================
       ! ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
       ! ======================================================================

       ! If wfoptalg>=10, the precondition matrix is kept constant during iteration ; otherwise it is recomputed
       call timab(1305,1,tsec)
       if (wfoptalg<10.or.iline==1) then
         call cg_precon(cwavef,zero,istwf_k,gs_hamk%kinpw_k,npw,nspinor,me_g0,optekin,pcon,direc,mpi_enreg%comm_fft)
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
       call timab(1305,2,tsec)

       ! ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
       ! ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
       call timab(1203,1,tsec)
       call cprj_update_oneband(direc,cprj_direc,gs_hamk,mpi_enreg)
       call timab(1203,2,tsec)
       call getcsc(scprod_csc,cpopt,direc,cwavef_bands,cprj_direc,cprj_cwavef_bands,&
&       gs_hamk,mpi_enreg,nband,tim_getcsc)
!       if (abs(xnorm-one)>tol10) then ! True if iline==1 and if input WFs are random
       if (iline==1) then
         ! We compensate the normalization of the current band
         scprod_csc(2*iband-1:2*iband) = scprod_csc(2*iband-1:2*iband)/xnorm
       end if
       ! Projecting again out all bands (not normalized).
       scprod = reshape(scprod_csc,(/2,nband/))
       call projbd(cg,direc,-1,icg,icg,istwf_k,mcg,mcg,nband,npw,nspinor,&
&       direc,scprod,1,tim_projbd,useoverlap,me_g0,mpi_enreg%comm_fft)
       if (iline==1) then
         ! Again we have to compensate the normalization of the current band.
         ! Indeed, by calling projbd we compute:
         ! |direc'_i> = |direc_i> - \sum_j <c'_j|S|direc_i>|c'_j>
         !            = |direc_i> - \sum_{j/=i} <c_j|S|direc_i>|c_j> - <c'_i|S|direc_i>|c'_i>
         ! where |c'_j> = |c_j> for j/=i and |c'_i> = xnorm.|c_i>
         ! As we compensated "scprod" before the call of projbd we actually computed:
         ! |direc'_i> = |direc_i> - \sum_{j/=i} <c_j|S|direc_i>|c_j> - <c_i|S|direc_i>|c'_i>
         !            = |direc_i> - \sum_{j/=i} <c_j|S|direc_i>|c_j> - xnorm.<c_i|S|direc_i>|c_i>
         ! The correct projected direction should be:
         ! |projdirec_i> = |direc_i> - \sum_j <c_j|S|direc_i>|c_j>
         !               = |direc_i> - \sum_{j/=i} <c_j|S|direc_i>|c_j> - <c_i|S|direc_i>|c_i>
         ! So:
         ! |projdirec_i> = |direc'_i> - <c_i|S|direc_i>|c_i> + xnorm.<c_i|S|direc_i>|c_i>
         !               = |direc'_i> - (1-xnorm) <c_i|S|direc_i>|c_i>
         !               = |direc'_i> - (1-xnorm)/xnorm <c_i|S|direc_i>|c'_i>
         !
         z_tmp = -scprod_csc(2*iband-1:2*iband)*(one-xnorm)/xnorm
         call timab(1305,1,tsec)
         do ispinor=1,nspinor
           igs=(ispinor-1)*npw
           do ipw=1+igs,npw+igs
             direc(1,ipw)=direc(1,ipw) + z_tmp(1)*cwavef(1,ipw) - z_tmp(2)*cwavef(2,ipw)
             direc(2,ipw)=direc(2,ipw) + z_tmp(1)*cwavef(2,ipw) + z_tmp(2)*cwavef(1,ipw)
           end do
         end do
         call timab(1305,2,tsec)
       end if
       ! Apply projbd to cprj_direc
       scprod=-scprod
       call timab(1303,1,tsec)
       call pawcprj_projbd(scprod,cprj_cwavef_bands,cprj_direc)
       call timab(1303,2,tsec)
       if (iline==1) then
         ! Same correction than for WFs
         z_tmp  = -scprod_csc(2*iband-1:2*iband)*(one-xnorm)/xnorm
         z_tmp2 = (/one,zero/)
         ! cprj = z_tmp*cprjx + z_tmp2*cprj
         call timab(1302,1,tsec)
         call pawcprj_zaxpby(z_tmp,z_tmp2,cprj_cwavef,cprj_direc)
         call timab(1302,2,tsec)
       end if

       ! ======================================================================
       ! ================= COMPUTE THE CONJUGATE-GRADIENT =====================
       ! ======================================================================

       call timab(1305,1,tsec)
       call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,direc_tmp,me_g0,mpi_enreg%comm_spinorfft)
       call timab(1305,2,tsec)

       dotgg=dotgg/xnormd**2

       ! MJV: added 5 Feb 2012 - causes divide by 0 on next iteration of iline
       if (abs(dotgg) < TINY(0.0_dp)*1.e50_dp) dotgg = TINY(0.0_dp)*1.e50_dp

       ! At first iteration, gamma is set to zero
       if (iline==1) then
         gamma=zero
         dotgp=dotgg
         call timab(1305,1,tsec)
         call cg_zcopy(npw*nspinor,direc,conjgr)
         call timab(1305,2,tsec)
         call timab(1302,1,tsec)
         call pawcprj_copy(cprj_direc,cprj_conjgr)
         call timab(1302,2,tsec)
         if (prtvol==-level)then
           write(message,'(a,es27.14e3)')' cgwf: dotgg = ',dotgg
           call wrtout(std_out,message,'PERS')
         end if
       else
         gamma=dotgg/dotgp
         dotgp=dotgg

         if (prtvol==-level)then
           write(message,'(a,2es27.14e3)')' cgwf: dotgg,gamma = ',dotgg,gamma
           call wrtout(std_out,message,'PERS')
         end if

         gamma=gamma*xnormd/xnormd_previous

         call timab(1305,1,tsec)
!$OMP PARALLEL DO
         do ipw=1,npw*nspinor
           conjgr(1,ipw)=direc(1,ipw)+gamma*conjgr(1,ipw)
           conjgr(2,ipw)=direc(2,ipw)+gamma*conjgr(2,ipw)
         end do
         call timab(1305,2,tsec)
         call timab(1302,1,tsec)
         call pawcprj_axpby(one,gamma,cprj_direc,cprj_conjgr)
         call timab(1302,2,tsec)
       end if
       call timab(1305,1,tsec)

       xnormd_previous=xnormd

       call getcsc(dot,cpopt,conjgr,conjgr,cprj_conjgr,cprj_conjgr,&
&       gs_hamk,mpi_enreg,1,tim_getcsc)

       ! ======================================================================
       ! ============ PROJECTION OF THE CONJUGATED GRADIENT ===================
       ! ======================================================================

       call getcsc(dot,cpopt,conjgr,cwavef,cprj_conjgr,cprj_cwavef,&
&       gs_hamk,mpi_enreg,1,tim_getcsc)
       dotr=dot(1)
       doti=dot(2)

       ! Project the conjugated gradient onto the current band
       call timab(1305,1,tsec)
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
       call timab(1305,2,tsec)
       call timab(1302,1,tsec)
       call pawcprj_copy(cprj_conjgr,cprj_direc)
       z_tmp  = (/-dotr,-doti/)
       z_tmp2 = (/one,zero/)
       call pawcprj_zaxpby(z_tmp,z_tmp2,cprj_cwavef,cprj_direc)
       call timab(1302,2,tsec)

       ! ======================================================================
       ! ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
       ! ======================================================================

       ! Compute norm of direc
       call getcsc(dot,cpopt,direc,direc,cprj_direc,cprj_direc,&
&       gs_hamk,mpi_enreg,1,tim_getcsc)
       xnormd=one/sqrt(abs(dot(1)))

       sij_opt=0
       ! Compute direc in real space
       call get_cwavefr(direc,direc_r,gs_hamk,mpi_enreg)
       ! Compute dhc = Re{<D|H|C>}
       call getchc(z_tmp,cpopt,cwavef,direc,cprj_cwavef,cprj_direc,cwavef_r,direc_r,&
         &          gs_hamk,zero,mpi_enreg,1,sij_opt,type_calc)
       dhc=z_tmp(1)
       dhc=dhc*xnormd

       ! Compute <D|H|D> or <D|(H-zshift)^2|D>
       call getchc(z_tmp,cpopt,direc,direc,cprj_direc,cprj_direc,direc_r,direc_r,&
&        gs_hamk,zero,mpi_enreg,1,sij_opt,type_calc)
       dhd=z_tmp(1)
       dhd=dhd*xnormd**2

       if (prtvol==-level) then
         write(message,'(a,3f26.14)') 'cgwf: chc,dhc,dhd=',chc,dhc,dhd
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

       ! ============================================================
       ! =========== GENERATE NEW wf(G), wf(r) and cprj =============
       ! ============================================================

       sintn=sinth*xnormd

       call timab(1305,1,tsec)
!$OMP PARALLEL DO
       do ipw=1,npw*nspinor
         cwavef(1,ipw)=cwavef(1,ipw)*costh+direc(1,ipw)*sintn
         cwavef(2,ipw)=cwavef(2,ipw)*costh+direc(2,ipw)*sintn
       end do
       do ispinor=1,nspinor
         do i3=1,gs_hamk%n6
           do i2=1,gs_hamk%n5
             do i1=1,gs_hamk%n4
               cwavef_r(1,i1,i2,i3,ispinor)=cwavef_r(1,i1,i2,i3,ispinor)*costh+direc_r(1,i1,i2,i3,ispinor)*sintn
               cwavef_r(2,i1,i2,i3,ispinor)=cwavef_r(2,i1,i2,i3,ispinor)*costh+direc_r(2,i1,i2,i3,ispinor)*sintn
             end do
           end do
         end do
       end do
       call timab(1305,2,tsec)
       if (cprj_update_lvl<=0.and.cprj_update_lvl>=-1) then
         call timab(1203,1,tsec)
         call cprj_update_oneband(cwavef,cprj_cwavef,gs_hamk,mpi_enreg)
         call timab(1203,1,tsec)
       else
         call timab(1302,1,tsec)
         call pawcprj_axpby(sintn,costh,cprj_direc,cprj_cwavef)
         call timab(1302,2,tsec)
       end if

       ! ======================================================================
       ! =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
       ! ======================================================================

       ! Compute delta(E)
       deltae=chc*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc

!      Check convergence and eventually exit
       if (iline==1) then
         deold=deltae
       else if (abs(deltae)<tolrde*abs(deold) .and. iline/=nline .and. wfopta10<2)then
         if(prtvol>=10)then
           write(message, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
&           ' cgwf: line',iline,&
&           ' deltae=',deltae,' < tolrde*',deold,' =>skip lines'
           call wrtout(std_out,message,'PERS')
         end if
         ! Update chc before exit
         call getchc(z_tmp,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
           &          gs_hamk,zero,mpi_enreg,1,sij_opt,type_calc)
         eig(iband)=z_tmp(1)
         nskip=nskip+2*(nline-iline)  ! Number of one-way 3D ffts skipped
         exit                         ! Exit from the loop on iline
       end if

       ! Update chc only if last iteration, otherwise it will be done at the beginning of the next one
       if (iline==nline) then
         call getchc(z_tmp,cpopt,cwavef,cwavef,cprj_cwavef,cprj_cwavef,cwavef_r,cwavef_r,&
           &          gs_hamk,zero,mpi_enreg,1,sij_opt,type_calc)
         eig(iband)=z_tmp(1)
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

 end do !  End big iband loop.

 !  ======================================================================
 !  ============= COMPUTE HAMILTONIAN IN WFs SUBSPACE ====================
 !  ======================================================================

 if (cprj_update_lvl<=2.and.cprj_update_lvl>=-1) then
   call timab(1203,1,tsec)
   call cprj_update(cg,cprj_cwavef_bands,gs_hamk,icg,nband,mpi_enreg)
   call timab(1203,2,tsec)
 end if

 sij_opt=0

 isubh=1
 isubh0=1
 do iband=1,nband
   cwavef => cwavef_bands(:,1+(iband-1)*npw*nspinor:iband*npw*nspinor)
   cprj_cwavef => cprj_cwavef_bands(:,nspinor*(iband-1)+1:nspinor*iband)
   ! Compute local+kinetic part
   call getghc(cpopt,cwavef,cprj_cwavef,direc,scwavef_dum,gs_hamk,gvnlxc,&
     &         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,3)
   isubh=isubh0
   call timab(1304,1,tsec)
   do jband=1,iband
     cwavef_left => cwavef_bands(:,1+(jband-1)*npw*nspinor:jband*npw*nspinor)
     call dotprod_g(subham(isubh),subham(isubh+1),istwf_k,npw*nspinor,2,cwavef_left,direc,me_g0,mpi_enreg%comm_spinorfft)
     isubh=isubh+2
   end do
   call timab(1304,2,tsec)
   cprj_cwavef_left => cprj_cwavef_bands(:,1:nspinor*iband)
   ! Add the nonlocal part
   call getchc(subham(isubh0:isubh0+2*iband-1),cpopt,cwavef,cwavef,&
     &          cprj_cwavef,cprj_cwavef_left,cwavef_r,cwavef_r,&
     &          gs_hamk,zero,mpi_enreg,iband,sij_opt,4)
   isubh0=isubh
 end do

 ! Debugging ouputs
 if(prtvol==-level)then
   isubh=1
   do iband=1,nband
     do jband=1,iband
       if (jband<=10) then
         write(message,'(i7,2f26.14)')isubh,subham(isubh:isubh+1)
         call wrtout(std_out,message,'PERS')
       end if
       isubh=isubh+2
     end do
   end do
 end if

 ! ===================
 ! FINAL DEALLOCATIONS
 ! ===================

 nullify(cwavef_left)
 nullify(cwavef_bands)
 nullify(cwavef)
 nullify(cprj_cwavef)
 nullify(cprj_cwavef_left)
 call pawcprj_free(cprj_direc)
 call pawcprj_free(cprj_conjgr)
 ABI_FREE(cprj_direc)
 ABI_FREE(cprj_conjgr)
 ABI_FREE(conjgr)
 ABI_FREE(scwavef_dum)
 ABI_FREE(gvnlxc)
 ABI_FREE(direc)
 ABI_FREE(pcon)
 ABI_FREE(scprod)
 ABI_FREE(scprod_csc)

 ABI_FREE(direc_tmp)
 ABI_FREE(cwavef_r)
 ABI_FREE(direc_r)

! Do not delete this line, needed to run with open MP
 write(unit=message,fmt=*) resid(1)

 call timab(1300,2,tsec)

 DBG_EXIT("COLL")

end subroutine cgwf_paw
!!***

!!****f* m_cgwf/mksubovl
!! NAME
!! mksubovl
!!
!! FUNCTION
!!   Compute the triangular matrix <c_m|S|c_n> for m<=n
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction <G|C band,k> coefficients for ALL bands
!!  cprj_cwavef_bands=<p_i|c_n> coefficients for all bands n
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  nband=number of bands.
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  subovl(nband*(nband+1)) = <c_m|S|c_n> for m<=n
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine mksubovl(cg,cprj_cwavef_bands,gs_hamk,icg,nband,subovl,mpi_enreg)
!Arguments ------------------------------------
 integer,intent(in) :: icg,nband
!arrays
 real(dp),intent(out) :: subovl(nband*(nband+1))
 type(pawcprj_type),intent(in),target :: cprj_cwavef_bands(:,:)
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(in),target :: cg(:,:)

!Local variables-------------------------------
 integer,parameter :: tim_getcsc=4
 integer :: cpopt,iband,isubh,nspinor,wfsize
 real(dp),pointer :: cwavef(:,:),cwavef_left(:,:)
 real(dp),pointer :: cwavef_bands(:,:)
 type(pawcprj_type),pointer :: cprj_cwavef(:,:),cprj_cwavef_left(:,:)

! *********************************************************************

 nspinor = gs_hamk%nspinor
 wfsize=gs_hamk%npw_k*nspinor
 cpopt=2

 cwavef_bands => cg(:,1+icg:nband*wfsize+icg)

 isubh=1
 do iband=1,nband
   cwavef => cwavef_bands(:,1+(iband-1)*wfsize:iband*wfsize)
   cprj_cwavef => cprj_cwavef_bands(:,nspinor*(iband-1)+1:nspinor*iband)
   cwavef_left => cwavef_bands(:,1:iband*wfsize)
   cprj_cwavef_left => cprj_cwavef_bands(:,1:nspinor*iband)
   ! Compute csc matrix
   call getcsc(subovl(isubh:isubh+2*iband-1),cpopt,cwavef,cwavef_left,&
     &          cprj_cwavef,cprj_cwavef_left,&
     &          gs_hamk,mpi_enreg,iband,tim_getcsc)
   isubh=isubh+2*iband
 end do

end subroutine mksubovl
!!***

!!****f* m_cgwf/cprj_update
!! NAME
!! cprj_update
!!
!! FUNCTION
!!   Compute the cprj = <p_i|c_n> from input cg (at fixed k point).
!!
!! INPUTS
!!  cg(2,mcg)=wavefunction <G|C band,k> coefficients for ALL bands
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  nband=number of bands.
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  cprj_cwavef_bands=<p_i|c_n> coefficients for all bands n
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine cprj_update(cg,cprj_cwavef_bands,gs_hamk,icg,nband,mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: icg,nband
!arrays
 type(pawcprj_type),intent(in),target :: cprj_cwavef_bands(:,:)
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(inout),target :: cg(:,:)

!Local variables-------------------------------
 integer :: choice,iband,wfsize
 real(dp),pointer :: cwavef(:,:),cwavef_bands(:,:)
 type(pawcprj_type),pointer :: cprj_cwavef(:,:)

 wfsize=gs_hamk%npw_k*gs_hamk%nspinor
 cwavef_bands => cg(:,1+icg:nband*wfsize+icg)

 choice = 1
 if (cprj_cwavef_bands(1,1)%ncpgr==3) then
   choice = 2
 end if

 do iband=1,nband
   cwavef => cwavef_bands(:,1+(iband-1)*wfsize:iband*wfsize)
   cprj_cwavef => cprj_cwavef_bands(:,gs_hamk%nspinor*(iband-1)+1:gs_hamk%nspinor*iband)

   call getcprj(choice,0,cwavef,cprj_cwavef,&
&    gs_hamk%ffnl_k,0,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&    gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,gs_hamk%natom,gs_hamk%nattyp,&
&    gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&    gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
 end do

end subroutine cprj_update
!!***

!!****f* m_cgwf/cprj_update_oneband
!! NAME
!! cprj_update_oneband
!!
!! FUNCTION
!!   Compute the cprj = <p_i|c_n> from input wavefunction (one band only).
!!
!! INPUTS
!!  cwavef(2,npw)=a wavefunction <G|C band,k>
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  cprj_cwavef=<p_i|c_n> coefficients of the input WF
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine cprj_update_oneband(cwavef,cprj_cwavef,gs_hamk,mpi_enreg)

!Arguments ------------------------------------
!arrays
 type(pawcprj_type),intent(inout) :: cprj_cwavef(:,:)
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(inout) :: cwavef(:,:)

!Local variables-------------------------------
 integer :: choice,wfsize

 wfsize=gs_hamk%npw_k*gs_hamk%nspinor

 choice = 1
 if (cprj_cwavef(1,1)%ncpgr==3) then
   choice = 2
 end if

 call getcprj(choice,0,cwavef,cprj_cwavef,&
&  gs_hamk%ffnl_k,0,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&  gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,gs_hamk%natom,gs_hamk%nattyp,&
&  gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&  gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)

end subroutine cprj_update_oneband
!!***

!!****f* m_cgwf/cprj_check
!! NAME
!! cprj_check
!!
!! FUNCTION
!!   Check if the cprj are consistent with the bands contained in cg.
!!   Useful for debugging only.
!!
!! INPUTS
!!  cwavef(2,npw)=a wavefunction <G|C band,k>
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  message=string to specify the context of the test
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine cprj_check(cg,cprj_cwavef_bands,gs_hamk,icg,nband,message,mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: icg,nband
!arrays
 type(pawcprj_type),intent(in),target :: cprj_cwavef_bands(:,:)
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(inout),target :: cg(:,:)
 character(len=*),intent(in) :: message

!Local variables-------------------------------
 integer :: choice,iband,ispinor,ncpgr,wfsize
 real(dp),pointer :: cwavef(:,:),cwavef_bands(:,:)
 integer :: iatom
 real(dp) :: re,ratio
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)

 write(std_out,'(a)') ''
 write(std_out,'(2a)') 'cprj_check : ',message

 choice = 1
 ncpgr  = 0
 if (cprj_cwavef_bands(1,1)%ncpgr==3) then
   choice = 2
   ncpgr  = 3
 end if
 write(std_out,'(a,i3)') 'ncpgr : ',ncpgr

 ABI_MALLOC(cprj_tmp,(gs_hamk%natom,gs_hamk%nspinor))
 call pawcprj_alloc(cprj_tmp,ncpgr,gs_hamk%dimcprj)

 wfsize=gs_hamk%npw_k*gs_hamk%nspinor
 cwavef_bands => cg(:,1+icg:nband*wfsize+icg)

 do iband=1,nband
   cwavef => cwavef_bands(:,1+(iband-1)*wfsize:iband*wfsize)
   call getcprj(choice,0,cwavef,cprj_tmp,&
&    gs_hamk%ffnl_k,0,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&    gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,gs_hamk%natom,gs_hamk%nattyp,&
&    gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&    gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
   do ispinor=1,gs_hamk%nspinor
     do iatom=1,gs_hamk%natom
       re = sum(abs(cprj_tmp(iatom,ispinor)%cp-cprj_cwavef_bands(iatom,gs_hamk%nspinor*(iband-1)+ispinor)%cp))
       if (re>tol6) then
         ratio = sum(abs(cprj_tmp(iatom,ispinor)%cp))/sum(abs(cprj_cwavef_bands(iatom,gs_hamk%nspinor*(iband-1)+ispinor)%cp))
         write(std_out,'(a)') 'cprj_check:'
         write(std_out,'(a,2i5,2es11.3e3)') 'iband,iatom:',iband,iatom,re,ratio
         write(std_out,'(a)') ''
         flush(std_out)
         ABI_ERROR('dif too large')
       end if
       if (ncpgr>0) then
         re = sum(abs(cprj_tmp(iatom,ispinor)%dcp-cprj_cwavef_bands(iatom,gs_hamk%nspinor*(iband-1)+ispinor)%dcp))
         if (re>tol6) then
           write(std_out,'(a)') 'cprj_check:'
           write(std_out,'(a,2i5,es11.3e3)') 'iband,iatom:',iband,iatom,re
           write(std_out,'(a)') ''
           flush(std_out)
           ABI_ERROR('dif too large (dcp)')
         end if
       end if
     end do
   end do
 end do

 write(std_out,'(a)') 'cprj_check : ok'
 flush(std_out)

 call pawcprj_free(cprj_tmp)
 ABI_FREE(cprj_tmp)

end subroutine cprj_check

!!****f* m_cgwf/cprj_check_oneband
!! NAME
!! cprj_check_oneband
!!
!! FUNCTION
!!   Check if the input cprj is consistent with the input WF.
!!   Useful for debugging only.
!!
!! INPUTS
!!  cwavef(2,npw)=a wavefunction <G|C band,k>
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  message=string to specify the context of the test
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine cprj_check_oneband(cwavef,cprj_cwavef,gs_hamk,message,mpi_enreg)

!Arguments ------------------------------------
!arrays
 type(pawcprj_type),intent(in),target :: cprj_cwavef(:,:)
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(inout) :: cwavef(:,:)
 character(len=*),intent(in) :: message

!Local variables-------------------------------
 integer :: choice,ncpgr,wfsize
 integer :: iatom
 real(dp) :: re
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)

 write(std_out,'(a)') ''
 write(std_out,'(2a)') 'cprj_check (oneband): ',message

 choice = 1
 ncpgr  = 0
 if (cprj_cwavef(1,1)%ncpgr==3) then
   choice = 2
   ncpgr  = 3
 end if
 write(std_out,'(a,i3)') 'ncpgr : ',ncpgr

 ABI_MALLOC(cprj_tmp,(gs_hamk%natom,1))
 call pawcprj_alloc(cprj_tmp,ncpgr,gs_hamk%dimcprj)

 wfsize=gs_hamk%npw_k*gs_hamk%nspinor

 call getcprj(choice,0,cwavef,cprj_tmp,&
&  gs_hamk%ffnl_k,0,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&  gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,gs_hamk%natom,gs_hamk%nattyp,&
&  gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&  gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
 do iatom=1,gs_hamk%natom
   re = sum(abs(cprj_tmp(iatom,1)%cp-cprj_cwavef(iatom,1)%cp))
   if (re>tol6) then
     write(std_out,'(a)') 'cprj_check (oneband):'
     write(std_out,'(a,i5,es11.3e3)') 'iatom:',iatom,re
     write(std_out,'(a)') ''
     flush(std_out)
     ABI_ERROR('dif too large')
   end if
   if (ncpgr>0) then
     re = sum(abs(cprj_tmp(iatom,1)%dcp-cprj_cwavef(iatom,1)%dcp))
     if (re>tol6) then
       write(std_out,'(a)') 'cprj_check (oneband):'
       write(std_out,'(a,i5,es11.3e3)') 'iatom:',iatom,re
       write(std_out,'(a)') ''
       flush(std_out)
       ABI_ERROR('dif too large (dcp)')
     end if
   end if
 end do

 write(std_out,'(a)') 'cprj_check : ok'
 flush(std_out)

 call pawcprj_free(cprj_tmp)
 ABI_FREE(cprj_tmp)

end subroutine cprj_check_oneband
!!***

!!****f* m_cgwf/get_cprj_id
!! NAME
!! get_cprj_id
!!
!! FUNCTION
!!   Get an id from cprj array.
!!   Useful for debugging only.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
real(dp) function get_cprj_id(cprj)

!Arguments ------------------------------------
!arrays
 type(pawcprj_type),intent(in) :: cprj(:,:)

!Local variables-------------------------------
 integer :: i1,i2
 real(dp) :: id_tmp

 id_tmp=zero
 do i2=1,size(cprj,2)
   do i1=1,size(cprj,1)
     id_tmp = id_tmp + sum(abs(cprj(i1,i2)%cp))
   end do
 end do
 get_cprj_id=id_tmp

end function get_cprj_id
!!***

!!****f* m_cgwf/get_cwavefr
!! NAME
!! get_cwavefr
!!
!! FUNCTION
!!   Compute the real space wavefunction by FFT of the input WF in reciprocal space.
!!
!! INPUTS
!!  cwavef(2,npw)=a wavefunction <G|C band,k>
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!  cwavef_r(2,n4,n5,n6,nspinor) = real space coefficients of the WF
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine get_cwavefr(cwavef,cwavef_r,gs_hamk,mpi_enreg)

!Arguments ------------------------------------
!arrays
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(in),target :: cwavef(:,:)
 real(dp),intent(out) :: cwavef_r(:,:,:,:,:)

!Local variables-------------------------------
 integer,parameter :: tim_fourwf=40
 integer :: n4,n5,n6,npw
 real(dp),parameter :: weight_fft = one
 real(dp), pointer :: cwavef_fft1(:,:),cwavef_fft2(:,:)
 real(dp),allocatable :: denpot_dum(:,:,:),fofgout_dum(:,:)

 n4=gs_hamk%ngfft(4);n5=gs_hamk%ngfft(5);n6=gs_hamk%ngfft(6)
 npw=gs_hamk%npw_k

 ABI_MALLOC(denpot_dum, (0,0,0))
 ABI_MALLOC(fofgout_dum, (0,0))

 if (gs_hamk%nspinor==1) then
   cwavef_fft1 => cwavef
 else
   cwavef_fft1 => cwavef(:,1:npw)
   cwavef_fft2 => cwavef(:,1+npw:2*npw)
 end if
 call fourwf(0,denpot_dum,cwavef_fft1,fofgout_dum,cwavef_r(:,:,:,:,1),gs_hamk%gbound_k,gs_hamk%gbound_k,gs_hamk%istwf_k,&
&  gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,gs_hamk%npw_fft_k,gs_hamk%npw_fft_k,&
&  n4,n5,n6,0,tim_fourwf,weight_fft,weight_fft)
 if (gs_hamk%nspinor==2) then
   call fourwf(0,denpot_dum,cwavef_fft2,fofgout_dum,cwavef_r(:,:,:,:,2),gs_hamk%gbound_k,gs_hamk%gbound_k,gs_hamk%istwf_k,&
&    gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,gs_hamk%npw_fft_k,gs_hamk%npw_fft_k,&
&    n4,n5,n6,0,tim_fourwf,weight_fft,weight_fft)
 end if

 ABI_FREE(denpot_dum)
 ABI_FREE(fofgout_dum)

end subroutine get_cwavefr
!!***

!!****f* ABINIT/enable_cgwf_paw
!!
!! NAME
!! enable_cgwf_paw
!!
!! FUNCTION
!! Return a logical value determining if the "cgwf_paw" implementation can be used or not
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!!  enable_cgwf_paw : true if "cgwf_paw" can be used
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

logical function enable_cgwf_paw(dtset)

   type(dataset_type), intent(in) :: dtset

   enable_cgwf_paw = .false.
   ! cgwf_paw is supposed to work for the following cases :
   if (dtset%usepaw==1.and.dtset%wfoptalg==10.and.dtset%paral_kgb==0) then
     enable_cgwf_paw = .true.
   end if
   ! but the following cases are not implemented yet:
   if (dtset%berryopt/=0.or.dtset%usefock/=0.or.sum(abs(dtset%nucdipmom))>tol16) then
     enable_cgwf_paw = .false.
   end if

end function enable_cgwf_paw
!!***

end module m_cgwf_paw
!!***
