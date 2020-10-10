!!****m* ABINIT/m_cgwf
!! NAME
!!  m_cgwf
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

module m_cgwf

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
 use m_getghc,        only : getghc
 use m_berrytk,       only : smatrix
 use m_nonlop,      only : nonlop
 use m_paw_overlap,   only : smatrix_k_paw
 use m_cgprj,         only : getcprj

 implicit none

 private
!!***

 public :: cgwf

!!***

contains
!!***

!!****f* m_cgwf/cgwf
!! NAME
!! cgwf
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
!!  subovl(nband*(nband+1)*use_subovl)=overlap matrix expressed in the WFs subspace
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
!!  1) cg should not be filtered and normalized: it should already be OK at input !
!!  2) Not sure that that the generalized eigenproblem (when gs_hamk%usepaw=1)
!!     is compatible with wfoptalg=2 or 3 (use of shifted square Hamiltonian) - to be verified
!!
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

subroutine cgwf(berryopt,cg,cgq,chkexit,cpus,dphase_k,dtefield,&
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
 logical :: gen_eigenpb, finite_field
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

!Initializations and allocations
 isubh=1;isubo=1
 nblock=(nband-1)/nbdblock+1
 istwf_k=gs_hamk%istwf_k
 wfopta10=mod(wfoptalg,10)
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
 ABI_ALLOCATE(gvnlx_direc,(2,npw*nspinor))
 ABI_ALLOCATE(vresid,(2,npw*nspinor))

 if (gen_eigenpb)  then
   ABI_ALLOCATE(gs_direc,(2,npw*nspinor))
 else
   ABI_ALLOCATE(gs_direc,(0,0))
 end if

 if (gen_eigenpb) then
   ABI_ALLOCATE(scwavef,(2,npw*nspinor))
   ABI_ALLOCATE(direc_tmp,(2,npw*nspinor))
 end if

 if (gen_eigenpb.and.(inonsc==1))  then
   ABI_MALLOC_OR_DIE(ghc_all,(2,nband*npw*nspinor), ierr)
 end if

 if (wfopta10==2.or.wfopta10==3)  then
   ABI_ALLOCATE(work,(2,npw*nspinor))
 end if

 if (gen_eigenpb.and.(wfopta10==2.or.wfopta10==3))  then
   ABI_ALLOCATE(swork,(2,npw*nspinor))
 else
   ABI_ALLOCATE(swork,(0,0))
 end if

 if (wfopta10==2 .or. wfopta10==3) then
   ABI_ALLOCATE(ghcws,(2,npw*nspinor))
   ABI_ALLOCATE(gh_direcws,(2,npw*nspinor))
   ABI_ALLOCATE(gvnlx_dummy,(2,npw*nspinor))
 else
   ABI_ALLOCATE(gvnlx_dummy,(0,0))
 end if

!if "generalized eigenproblem", not sure of wfoptalg=2,3 algorithms
 if ((gen_eigenpb).and.(wfopta10==2.or.wfopta10==3)) then
   write(message, '(a,a,a,a,a)' ) &
&   'Conjugate gradient algorithm not tested with',ch10,&
&   'wfoptalg=2 or 3 and usepaw==1 !',ch10,&
&   'Program will continue at your own risk...'
   MSG_WARNING(message)
 end if

!Electric field: definition of local variables:
!detovc(1:2,ifor,idir) determinant of the overlap matrix
!S_{nm}(k,k+dk)=<u_{n,k}|u_{m,k+dk}>, with the states at
!k as bras (neighbor is specified by ifor and idir)
!detovd                same as detovc but with <u_{n,k}| replaced by
!<D| (search direction) in the band-th line
!grad_berry(1:2,ipw)   Berry phase term contribution to the gradient vector
!hel(ifor,idir)        helicity of the ellipse associated w/ (ifor,idir)
!bcut(ifor,idir)       branch cut of the ellipse associated w/ (ifor,idir)
!theta_min             optimal angle theta in line_minimization when electric
!field is on
!grad_total(1:2,ipw)   total gradient (zero field term plus Berry phase term)

 finite_field = ( (berryopt == 4) .or. (berryopt == 6) .or. (berryopt == 7) .or.   &
& (berryopt == 14) .or. (berryopt == 16) .or. (berryopt == 17) )
 ncpgr = 0 ! do not think the cprj's here need gradients (no force computation in cgwf)

 if (finite_field) then

!  ji : These could be a couple of new input variables (but it is OK to define them here)
   ikptf = dtefield%i2fbz(ikpt)
   ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind
   mcg_q = mpw*mband*nspinor
   ABI_ALLOCATE(detovc,(2,2,3))
   ABI_ALLOCATE(detovd,(2,2,3))
   ABI_ALLOCATE(grad_berry,(2,npw*nspinor))
   ABI_ALLOCATE(cg1_k,(2,mpw))
   ABI_ALLOCATE(cgq_k,(2,mcg_q))
   ABI_ALLOCATE(grad_total,(2,npw*nspinor))
   ABI_ALLOCATE(sflag_k,(mband))
   ABI_ALLOCATE(pwind_k,(mpw))
   ABI_ALLOCATE(pwnsfac_k,(4,mpw))
   ABI_ALLOCATE(smat_k,(2,mband,mband))
   ABI_ALLOCATE(smat_inv,(2,mband,mband))
!  now the special features if using PAW
   if (gs_hamk%usepaw /= 0) then
     ABI_ALLOCATE(smat_k_paw,(2,gs_hamk%usepaw*mband,gs_hamk%usepaw*mband))
     smat_k_paw(:,:,:) = zero
!    the following are arguments to nonlop used to apply the on-site dipole to direc vector
     choice = 1 ! only apply projectors
     paw_opt = 1 ! only apply Dij
     signs = 2 ! apply nonlop to vector in k space
!    following two items are the nonlocal potential strength dij due to the on-site dipoles
     dimenlc1 = 2*gs_hamk%lmnmax*(gs_hamk%lmnmax+1)/2
     dimenlr1 = gs_hamk%lmnmax*(gs_hamk%lmnmax+1)/2
     dimenl2 = natom
!    cprj structures for finite_field case
     ABI_ALLOCATE(dimlmn,(natom))
     do iatom = 1, natom
       itypat = gs_hamk%typat(iatom)
       dimlmn(iatom)=dtefield%lmn_size(itypat)
     end do
     ABI_ALLOCATE(dimlmn_srt,(natom))
     iatom = 0
     do itypat = 1, gs_hamk%ntypat
       do iat = 1, gs_hamk%nattyp(itypat)
         iatom = iatom + 1
         dimlmn_srt(iatom) = dtefield%lmn_size(itypat)
       end do
     end do
     ABI_ALLOCATE(ikptf_recv,(nproc_distrb))
     ABI_DATATYPE_ALLOCATE(cprj_k,(natom,mband*nspinor))
     ABI_DATATYPE_ALLOCATE(cprj_kb,(natom,mband*nspinor))
     ABI_DATATYPE_ALLOCATE(cprj_direc,(natom,mband*nspinor))
     ABI_DATATYPE_ALLOCATE(cprj_band_srt,(natom,nspinor))
     ABI_DATATYPE_ALLOCATE(cprj_gat,(natom,mband*nspinor*nproc_distrb))
     call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_direc,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_band_srt,ncpgr,dimlmn_srt)
     call pawcprj_alloc(cprj_gat,ncpgr,dimlmn)
     if (nkpt /= dtefield%fnkpt) then
       ABI_DATATYPE_ALLOCATE(cprj_fkn,(natom,mband*nspinor))
       ABI_DATATYPE_ALLOCATE(cprj_ikn,(natom,mband*nspinor))
       call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
       call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
     end if
   else
     ABI_ALLOCATE(smat_k_paw,(0,0,0))
     ABI_ALLOCATE(dimlmn,(0))
     ABI_ALLOCATE(dimlmn_srt,(0))
     ABI_ALLOCATE(ikptf_recv,(0))
     ABI_DATATYPE_ALLOCATE(cprj_k,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_kb,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_direc,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_band_srt,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_gat,(0,0))
   end if
 end if ! finite_field

 ! ======================================================================
 ! If generalized eigenproblem: has to know <g|S|c> for all
 ! bands (for orthogonalization purpose); take benefit of this
 ! calculation to compute <g|H|c> at the same time, and cprj_k if finite_field
 ! ======================================================================
 if (gen_eigenpb.and.inonsc==1) then
   do iblock=1,nblock
     ibandmin=1+(iblock-1)*nbdblock
     ibandmax=min(iblock*nbdblock,nband)
     do iband=ibandmin,ibandmax
       ibdblock=iband-(iblock-1)*nbdblock
       icg_shift=npw*nspinor*(iband-1)+icg
       igsc_shift=npw*nspinor*(iband-1)+igsc

       call cg_zcopy(npw*nspinor,cg(1,1+icg_shift),cwavef)

!      Compute <g|H|c>
!      By setting ieigen to iband, Fock contrib. of this iband to the energy will be calculated
       call fock_set_ieigen(gs_hamk%fockcommon,iband)
       sij_opt=1
       if (finite_field .and. gs_hamk%usepaw == 1) then
         call getghc(0,cwavef,cprj_band_srt,ghc,scwavef,gs_hamk,gvnlxc,&
&         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
         call pawcprj_put(gs_hamk%atindx,cprj_band_srt,cprj_k,natom,iband,0,ikpt,&
&         1,isppol,mband,1,natom,1,mband,dimlmn,nspinor,nsppol,0,&
&         mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)
       else
         call getghc(cpopt,cwavef,cprj_dum,ghc,scwavef,gs_hamk,gvnlxc,&
&         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
       end if

       call cg_zcopy(npw*nspinor,ghc,ghc_all(1,1+icg_shift-icg))
       call cg_zcopy(npw*nspinor,scwavef,gsc(1,1+igsc_shift))

     end do
   end do
 end if

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
     if (finite_field) then
       detovc(:,:,:) = zero ; detovd(:,:,:) = zero
       phase_init(:) = zero
       dphase_aux1(:) = zero
       phase_end(:) = zero
       bcut(:,:) = zero
       hel(:,:) = 0
     end if

     ! Extraction of the vector that is iteratively updated
     call cg_zcopy(npw*nspinor,cg(1,1+icg_shift),cwavef)

     ! If generalized eigenproblem: extraction of the overlap information
     if (gen_eigenpb) then
       call cg_zcopy(npw*nspinor,gsc(1,1+igsc_shift),scwavef)
     end if

     ! Normalize incoming wf (and S.wf, if generalized eigenproblem):
     ! WARNING : It might be interesting to skip the following operation.
     ! The associated routines should be reexamined to see whether cwavef is not already normalized.
     if (gen_eigenpb) then
       call dotprod_g(dotr,doti,istwf_k,npw*nspinor,2,cwavef,scwavef,me_g0,mpi_enreg%comm_spinorfft)
       dotr=sqrt(dotr**2+doti**2); xnorm=one/sqrt(dotr)
       call cg_zscal(npw*nspinor,(/xnorm,zero/),cwavef)
       call cg_zscal(npw*nspinor,(/xnorm,zero/),scwavef)
     else
       call sqnorm_g(dotr,istwf_k,npw*nspinor,cwavef,me_g0,mpi_enreg%comm_fft)
       xnorm=one/sqrt(abs(dotr))
       call cg_zscal(npw*nspinor,(/xnorm,zero/),cwavef)
     end if

     if (prtvol==-level) then
       write(message,'(a,f14.6)')' cgwf: xnorm = ',xnorm
       call wrtout(std_out,message,'PERS')
     end if

     ! Compute (or extract) <g|H|c>
     if (gen_eigenpb.and.(inonsc==1)) then

!$OMP PARALLEL DO PRIVATE(ipw)
       do ipw=1,npw*nspinor
         ghc(1,ipw)=xnorm*ghc_all(1,ipw+icg_shift-icg)
         ghc(2,ipw)=xnorm*ghc_all(2,ipw+icg_shift-icg)
       end do

     else
!      By setting ieigen to iband, Fock contrib. of this iband to the energy will be calculated
       call fock_set_ieigen(gs_hamk%fockcommon,iband)
       sij_opt=0
       call getghc(cpopt,cwavef,cprj_dum,ghc,gsc_dummy,gs_hamk,gvnlxc,&
&       eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
     end if


     ! Minimisation of the residual: compute <G|(H-zshift)^2|C iband,k>
     if(wfopta10==2 .or. wfopta10==3) then
       ghcws(:,:)=ghc(:,:)
       if (gen_eigenpb) then
         sij_opt=1
         work(:,:)=ghc(:,:)-zshift(iband)*scwavef(:,:)
       else
         sij_opt=0
         work(:,:)=ghc(:,:)-zshift(iband)*cwavef(:,:)
       end if
       call getghc(cpopt,work,cprj_dum,ghc,swork,gs_hamk,gvnlx_dummy,&
&       eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
       if (gen_eigenpb) then
         ghc(:,:)=ghc(:,:)-zshift(iband)*swork(:,:)
       else
         ghc(:,:)=ghc(:,:)-zshift(iband)*work(:,:)
       end if
     end if

     ! ======================================================================
     ! ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
     ! ======================================================================
     if(nline/=0)then
       do iline=1,nline

         ! === COMPUTE THE RESIDUAL ===

         ! Compute lambda = <C|H|C> or <C|(H-zshift)**2|C>
         call dotprod_g(chc,doti,istwf_k,npw*nspinor,1,cwavef,ghc,me_g0,mpi_enreg%comm_spinorfft)
         lam0=chc

         ! Check that lam0 is decreasing on succeeding lines:
         if (.not.finite_field) then
           if (iline==1) then
             lamold=lam0
           else
             if (lam0 > lamold+tol12) then
               write(message, '(a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&               'New trial energy at line ',iline,' = ',lam0,ch10,&
&               'is higher than former =',lamold,ch10
               MSG_WARNING(message)
             end if
             lamold=lam0
           end if
         end if

         ! Compute residual vector:
         ! Note that vresid is precomputed to garantee cancellation of errors
         ! and allow residuals to reach values as small as 1.0d-24 or better.

         if (wfopta10<=1) then
           eval=chc
           if (gen_eigenpb) then

!$OMP PARALLEL DO
             do ipw=1,npw*nspinor
               vresid(1,ipw)=ghc(1,ipw)-chc*scwavef(1,ipw)
               vresid(2,ipw)=ghc(2,ipw)-chc*scwavef(2,ipw)
             end do
           else
!$OMP PARALLEL DO
             do ipw=1,npw*nspinor
               vresid(1,ipw)=ghc(1,ipw)-chc*cwavef(1,ipw)
               vresid(2,ipw)=ghc(2,ipw)-chc*cwavef(2,ipw)
             end do
           end if
         else
           call dotprod_g(eval,doti,istwf_k,npw*nspinor,1,cwavef,ghcws,me_g0,mpi_enreg%comm_spinorfft)
           if (gen_eigenpb) then
!$OMP PARALLEL DO
             do ipw=1,npw*nspinor
               vresid(1,ipw)=ghcws(1,ipw)-eval*scwavef(1,ipw)
               vresid(2,ipw)=ghcws(2,ipw)-eval*scwavef(2,ipw)
             end do
           else
!$OMP PARALLEL DO
             do ipw=1,npw*nspinor
               vresid(1,ipw)=ghcws(1,ipw)-eval*cwavef(1,ipw)
               vresid(2,ipw)=ghcws(2,ipw)-eval*cwavef(2,ipw)
             end do
           end if
         end if

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

         ! Electric field: compute the gradient of the Berry phase part of the energy functional.
         ! See PRL 89, 117602 (2002) [[cite:Souza2002]], grad_berry(:,:) is the second term of Eq. (4)
         if (finite_field) then

           call make_grad_berry(cg,cgq,cprj_k,detovc,dimlmn,dimlmn_srt,direc,dtefield,grad_berry,&
&           gs_hamk,iband,icg,ikpt,isppol,mband,mcg,mcgq,mkgq,mpi_enreg,mpw,natom,nkpt,npw,npwarr,&
&           nspinor,nsppol,pwind,pwind_alloc,pwnsfac,pwnsfacq)

!          Add grad_berry to direc and store original gradient
           direc(:,:) = direc(:,:) + grad_berry(:,:)
           grad_total(:,:) = direc(:,:)
!          DEBUG: check that grad_berry is orthogonal to the occupied manifold at k
!          do jband = 1, dtefield%mband_occ
!          dotr = zero  ;  doti = zero
!          do ipw = 1, npw*nspinor
!          if(.not.gen_eigenpb) then
!          dotr = dotr + cg(1,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(1,ipw) + &
!          &                 cg(2,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(2,ipw)
!          doti = doti + cg(1,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(2,ipw) - &
!          &                 cg(2,icg + (jband-1)*npw*nspinor + ipw)*grad_berry(1,ipw)
!          end if
!          end do
!          if ((abs(dotr) > tol12).or.(abs(doti) > tol12)) then
!          write(std_out,'(a)')'cgwf-berry : ERROR (orthogonality)'
!          write(std_out,'(3(2x,i3),2(5x,e16.9))')ikpt,iband,jband,dotr,doti
!          stop
!          end if
!          end do
!          ENDDEBUG
         end if   ! finite_field

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

           ! Minimisation of the residual: must precondition twice
           ! (might make only one call, with modified precon routine - might also make a shift !!!)
           if(wfopta10==2 .or. wfopta10==3)then
             call cg_precon(cwavef,zero,istwf_k,kinpw,npw,nspinor,me_g0,optekin,pcon,direc,mpi_enreg%comm_fft)
             if(iline==1)then
!$OMP PARALLEL DO
               do ipw=1,npw
                 pcon(ipw)=pcon(ipw)**2
                 pcon(ipw)=pcon(ipw)**2
               end do
             end if
           end if
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

         if (finite_field) then
           call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,grad_total,me_g0,mpi_enreg%comm_spinorfft)
           !DEBUG (electric field)
           !check that the dotproduct is real
           !if (abs(doti) > tol8) then
           !write(std_out,*) ' cgwf-berry: ERROR'
           !write(std_out,*) ' doti = ',doti
           !stop
           !end if
           !ENDDEBUG
         else
           if (gen_eigenpb) then
             call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,direc_tmp,me_g0,mpi_enreg%comm_spinorfft)
           else
             call dotprod_g(dotgg,doti,istwf_k,npw*nspinor,1,direc,ghc,me_g0,mpi_enreg%comm_spinorfft)
           end if
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

         if (gen_eigenpb) then
           call dotprod_g(dotr,doti,istwf_k,npw*nspinor,3,scwavef,conjgr,me_g0,mpi_enreg%comm_spinorfft)
         else
           call dotprod_g(dotr,doti,istwf_k,npw*nspinor,3,cwavef,conjgr,me_g0,mpi_enreg%comm_spinorfft)
         end if

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

         call getghc(cpopt,direc,cprj_dum,gh_direc,gs_direc,gs_hamk,gvnlx_direc,&
&         eval,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)

         if(wfopta10==2 .or. wfopta10==3)then
           ! Minimisation of the residual, so compute <G|(H-zshift)^2|D>
           gh_direcws(:,:)=gh_direc(:,:)
           if (gen_eigenpb) then
             sij_opt=1
             work(:,:)=gh_direc(:,:)-zshift(iband)*gs_direc(:,:)
           else
             sij_opt=0
             work(:,:)=gh_direc(:,:)-zshift(iband)*direc(:,:)
           end if

           call getghc(cpopt,work,cprj_dum,gh_direc,swork,gs_hamk,gvnlx_dummy,&
&           eval,mpi_enreg,1,prtvol,0,tim_getghc,0)

           if (gen_eigenpb) then
             gh_direc(:,:)=gh_direc(:,:)-zshift(iband)*swork(:,:)
           else
             gh_direc(:,:)=gh_direc(:,:)-zshift(iband)*work(:,:)
           end if
         end if

         ! In case of generalized eigenproblem, compute now the norm of the conjugated gradient
         if (gen_eigenpb) then
           call dotprod_g(dotr,doti,istwf_k,npw*nspinor,1,direc,gs_direc,me_g0,mpi_enreg%comm_spinorfft)
           xnorm=one/sqrt(abs(dotr))
         end if

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

         if (.not.finite_field) then
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

         else
           ! In case the eletric field is on, the line minimization has to be done numerically

           ! Compute determinant of the overlap matrix where in the band-th line
           ! the wavefunction is replaced by the search direction
           job = 10 ; shiftbd = 0
           do idir = 1, 3
             ! do not do this for efield_dot(idir)=0
             if (abs(dtefield%efield_dot(idir)) < tol12) cycle
             do ifor = 1, 2
               ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
               if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
                 itrs = 10
               else
                 itrs = 0
               end if
               ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
               npw_k2 = npwarr(ikpt2)
               pwind_k(1:npw) = pwind(ikgf+1:ikgf+npw,ifor,idir)
               pwnsfac_k(1:2,1:npw) = pwnsfac(1:2,ikgf+1:ikgf+npw)
               sflag_k(:) = dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir)
               smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir)

               if (mpi_enreg%nproc_cell > 1) then
                 icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
                 cgq_k(:,1:dtefield%mband_occ*nspinor*npw_k2) = &
&                 cgq(:,icg1+1:icg1+dtefield%mband_occ*nspinor*npw_k2)
                 idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
                 pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
               else
                 icg1 = dtefield%cgindex(ikpt2,isppol)
                 cgq_k(:,1:dtefield%mband_occ*nspinor*npw_k2) = &
&                 cg(:,icg1+1:icg1+dtefield%mband_occ*nspinor*npw_k2)
                 idum1=dtefield%fkgindex(ikpt2f)
                 pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
               end if

               icg1 = 0 ; ddkflag = 0
               if (gen_eigenpb) then
!$OMP PARALLEL DO
                 do ipw=1,npw*nspinor
                   direc_tmp(1,ipw)=direc(1,ipw)*xnorm
                   direc_tmp(2,ipw)=direc(2,ipw)*xnorm
                 end do
                 ! need cprj corresponding to direc_tmp in order to make smat_k_paw properly
                 call getcprj(1,0,direc_tmp,cprj_band_srt,&
&                 gs_hamk%ffnl_k,0,gs_hamk%indlmn,gs_hamk%istwf_k,gs_hamk%kg_k,&
&                 gs_hamk%kpg_k,gs_hamk%kpt_k,gs_hamk%lmnmax,gs_hamk%mgfft,&
&                 mpi_enreg,gs_hamk%natom,gs_hamk%nattyp,gs_hamk%ngfft,gs_hamk%nloalg,&
&                 gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,gs_hamk%phkxred,gs_hamk%ph1d,&
&                 gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)

                 call pawcprj_copy(cprj_k,cprj_direc)
                 call pawcprj_put(gs_hamk%atindx,cprj_band_srt,cprj_direc,gs_hamk%natom,iband,0,&
&                 ikpt,1,isppol,mband,1,gs_hamk%natom,1,mband,dimlmn,gs_hamk%nspinor,nsppol,0,&
&                 mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)

!                icp1=dtefield%mband_occ*(ikptf-1)
                 icp2=mband*nspinor*(ikpt2-1)
                 call pawcprj_get(gs_hamk%atindx,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,isppol,&
&                 mband,dtefield%fnkpt,natom,mband,mband,nspinor,nsppol,0,&
&                 mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)

                 if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
                   call pawcprj_copy(cprj_kb,cprj_ikn)
                   call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
&                   dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
&                   dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
&                   dtefield%lmax,dtefield%lmnmax,mband,natom,dtefield%mband_occ,nspinor,&
&                   dtefield%nsym,gs_hamk%ntypat,gs_hamk%typat,dtefield%zarot)
                   call pawcprj_copy(cprj_fkn,cprj_kb)
                 end if

                 call smatrix_k_paw(cprj_direc,cprj_kb,dtefield,idir,ifor,mband,natom,smat_k_paw,gs_hamk%typat)

                 call smatrix(direc_tmp,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&                 itrs,job,iband,npw*nspinor,mcg_q,mpw,iband,&
&                 mpw,dtefield%mband_occ,dtefield%nband_occ(isppol),&
&                 npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&                 shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
               else
                 call smatrix(direc,cgq_k,cg1_k,ddkflag,dtm_k,icg1,icg1,&
&                 itrs,job,iband,npw*nspinor,mcg_q,mpw,iband,&
&                 mpw,dtefield%mband_occ,dtefield%nband_occ(isppol),&
&                 npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
&                 shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
               end if
               detovd(:,ifor,idir) = dtm_k(:) ! Store the determinant of the overlap
!              matrix (required to compute theta_min)
!              DEBUG
!              write(std_out,*)'cgwf-berry: detovc and detovd'
!              write(std_out,*)detovc(:,ifor,idir)
!              write(std_out,*)detovd(:,ifor,idir)
!              write(std_out,*)'smat_k'
!              do jband = 1, 4
!              write(std_out,'(4(2x,e14.6))')smat_k(1,jband,:)
!              write(std_out,'(4(2x,e14.6))')smat_k(2,jband,:)
!              write(std_out,*)
!              end do
!              ENDDEBUG
             end do  ! ifor
           end do    ! idir

           call linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,&
&           dphase_aux1,dtefield%efield_dot,iline,&
&           dtefield%fnkpt,dtefield%nstr,hel,phase_end,&
&           phase_init,dtefield%sdeg,sinth,thetam)
!          DEBUG
!          if (mpi_enreg%me == 1) then
!          write(std_out,*)'after linemin '
!          write(std_out,'(a,3(2x,f16.9))')'phase_init  = ',phase_init(:)
!          write(std_out,'(a,3(2x,f16.9))')'phase_end   = ',phase_end(:)
!          write(std_out,'(a,3(2x,f16.9))')'dphase_aux1 = ',dphase_aux1(:)
!          write(std_out,*) 'thetam',thetam
!          end if
!          ENDDEBUG
         end if  ! finite_field

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

         if (gen_eigenpb) then
!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             scwavef(1,ipw)=scwavef(1,ipw)*costh+gs_direc(1,ipw)*sintn
             scwavef(2,ipw)=scwavef(2,ipw)*costh+gs_direc(2,ipw)*sintn
!            gsc(1,ipw+igsc_shift)=scwavef(1,ipw)
!            gsc(2,ipw+igsc_shift)=scwavef(2,ipw)
           end do
!          call cg_zaxpby(npw*nspinor,(/sintn,zero/),gs_direc,(/costh,zero/),scwavef)
           call cg_zcopy(npw*nspinor,scwavef,gsc(1,1+igsc_shift))

           if (finite_field) then  ! must update cprj for the new wavefunction
             call getcprj(1,0,cwavef,cprj_band_srt,&
&             gs_hamk%ffnl_k,0,gs_hamk%indlmn,istwf_k,gs_hamk%kg_k,gs_hamk%kpg_k,gs_hamk%kpt_k,&
&             gs_hamk%lmnmax,gs_hamk%mgfft,mpi_enreg,natom,gs_hamk%nattyp,&
&             gs_hamk%ngfft,gs_hamk%nloalg,gs_hamk%npw_k,gs_hamk%nspinor,gs_hamk%ntypat,&
&             gs_hamk%phkxred,gs_hamk%ph1d,gs_hamk%ph3d_k,gs_hamk%ucvol,gs_hamk%useylm)
             call pawcprj_put(gs_hamk%atindx,cprj_band_srt,cprj_k,gs_hamk%natom,iband,0,ikpt,&
&             1,isppol,mband,1,gs_hamk%natom,1,mband,dimlmn,gs_hamk%nspinor,nsppol,0,&
&             mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)
           end if
         end if

         if(wfopta10==2 .or. wfopta10==3)then
           ! Need to keep track of ghcws, in order to avoid recomputing it
!$OMP PARALLEL DO
           do ipw=1,npw*nspinor
             ghcws(1,ipw)=ghcws(1,ipw)*costh + gh_direcws(1,ipw)*sintn
             ghcws(2,ipw)=ghcws(2,ipw)*costh + gh_direcws(2,ipw)*sintn
           end do
!          call cg_zaxpby(npw*nspinor,(/sintn,zero/),gh_direcws,(/costh,zero/),ghcws)
         end if

         ! ======================================================================
         ! =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
         ! ======================================================================

         ! Compute delta(E)
         if (.not.finite_field) then
           deltae=chc*(costh**2-1._dp)+dhd*sinth**2+2._dp*costh*sinth*dhc
         else
           ! Compute deltae
           call etheta(bcut,chc,detovc,detovd,dhc,dhd,dtefield%efield_dot,e0,e1,&
&           hel,dtefield%fnkpt,dtefield%nstr,dtefield%sdeg,thetam)
           theta = zero

           call etheta(bcut,chc,detovc,detovd,dhc,dhd,&
&           dtefield%efield_dot,e0_old,e1_old,&
&           hel,dtefield%fnkpt,dtefield%nstr,dtefield%sdeg,theta)
           deltae = e0 - e0_old
!          DEBUG
!          write(std_out,*) 'e0, e0_old, deltae', e0, e0_old, deltae
!          ENDDEBUG
!          Check that e0 is decreasing on succeeding lines:
!          if (deltae > zero) then
           if (deltae > tol12) then ! exploring different checks for finit_field
             write(message, '(3a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&             '  (electric field)',ch10,&
&             '  New trial energy at line',iline,' = ',e0,ch10,&
&             '  is higher than former:',e0_old,ch10
             MSG_WARNING(message)
           end if
         end if         ! finite_field

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

       ! Additional computations in case of electric field
       if (finite_field) then
         ! Bring present contribution to dphasek(idir) into [-pi,pi]
         do idir = 1, 3
           dphase_aux2 = mod(phase_end(idir) - phase_init(idir) + 100*two_pi,two_pi)
           if (dphase_aux2 > pi) dphase_aux2 = dphase_aux2 - two_pi
           ! DEBUG
           ! dphase_aux1(idir)=mod(dphase_aux1(idir)+100*two_pi,two_pi)
           ! if(dphase_aux1(idir)>pi) dphase_aux1(idir)=dphase_aux1(idir)-two_pi
           ! diff = dphase_aux2 - dphase_aux1(idir)
           ! if (abs(diff) > tol10) then
           ! write(std_out,*)'cgwf-berry: ERROR'
           ! write(std_out,'(a,3(2x,i3),f16.9)')'ikpt,iband,idir,diff',ikpt,iband,idir,diff
           ! stop
           ! end if
           ! write(100,*) idir, dphase_aux2
           ! ENDDEBUG
           dphase_k(idir) = dphase_k(idir) + dphase_aux2
           ! DEBUG
           ! write(std_out,*) 'idir,phase_init,phase_end,dphase_k'
           ! write(std_out,*) idir,phase_init(idir),phase_end(idir),dphase_k(idir)
           ! ENDDEBUG
         end do
       end if   ! finite_field

     else ! nline==0 , needs to provide a residual
       resid(iband)=-one
     end if ! End nline==0 case

     ! ======================================================================
     ! =============== END OF CURRENT BAND: CLEANING ========================
     ! ======================================================================

     ! It was checked that getghc is NOT needed here : equivalent results with the copy below.
     if(wfopta10==2 .or. wfopta10==3) ghc(:,:)=ghcws(:,:)

     if (finite_field) dtefield%sflag(:,ikpt + (isppol-1)*nkpt,:,:) = 0

     ! At the end of the treatment of a set of bands, write the number of one-way 3D ffts skipped
     if (xmpi_paral==0 .and. mpi_enreg%paral_kgb==0 .and. iband==nband .and. prtvol/=0) then
       write(message,'(a,i0)')' cgwf: number of one-way 3D ffts skipped in cgwf until now =',nskip
       call wrtout(std_out,message,'PERS')
     end if

   end do !  End big iband loop. iband in a block

   !  ======================================================================
   !  ============= COMPUTE HAMILTONIAN IN WFs SUBSPACE ====================
   !  ======================================================================
   call mksubham(cg,ghc,gsc,gvnlxc,iblock,icg,igsc,istwf_k,&
&   isubh,isubo,mcg,mgsc,nband,nbdblock,npw,&
&   nspinor,subham,subovl,subvnlx,use_subovl,use_subvnlx,me_g0)

 end do ! iblock End loop over block of bands

 if (allocated(dimlmn_srt)) then
   ABI_DEALLOCATE(dimlmn_srt)
 end if

 if (finite_field .and. gs_hamk%usepaw == 1) then ! store updated cprjs for this kpt
   ! switch from ikptf to ikpt
   ikptf = ikpt
   call xmpi_allgather(ikptf,ikptf_recv,spaceComm_distrb,ierr)
   call pawcprj_mpi_allgather(cprj_k,cprj_gat,natom,nspinor*mband,1,dimlmn,ncpgr,nproc_distrb,&
&   spaceComm_distrb,ierr,rank_ordered=.true.)
   do iproc = 1, nproc_distrb
     icp2=nspinor*mband*(iproc-1)
     call pawcprj_get(gs_hamk%atindx1,cprj_k,cprj_gat,natom,1,icp2,ikpt,0,isppol,&
&     mband,nproc_distrb,natom,mband,mband,nspinor,nsppol,0,&
&     mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)
!    ikptf = ikptf_recv(iproc)
     icp1 = nspinor*mband*(ikptf_recv(iproc)-1)
     call pawcprj_put(gs_hamk%atindx1,cprj_k,dtefield%cprj,natom,1,icp1,ikpt,0,isppol,&
&     mband,nkpt,natom,mband,mband,dimlmn,nspinor,nsppol,0,&
&     mpicomm=spaceComm_distrb,proc_distrb=mpi_enreg%proc_distrb)
   end do
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
   ABI_DEALLOCATE(scwavef)
   ABI_DEALLOCATE(direc_tmp)
 end if

 if (gen_eigenpb.and.(inonsc==1))  then
   ABI_DEALLOCATE(ghc_all)
 end if

 if(wfopta10==2.or.wfopta10==3)  then
   ABI_DEALLOCATE(ghcws)
   ABI_DEALLOCATE(gh_direcws)
 end if
 ABI_DEALLOCATE(gvnlx_dummy)

 if(wfopta10==2.or.wfopta10==3)  then
   ABI_DEALLOCATE(work)
 end if

 ABI_DEALLOCATE(swork)

 if (finite_field) then
   ABI_DEALLOCATE(cg1_k)
   ABI_DEALLOCATE(cgq_k)
   ABI_DEALLOCATE(detovc)
   ABI_DEALLOCATE(detovd)
   ABI_DEALLOCATE(grad_berry)
   ABI_DEALLOCATE(sflag_k)
   ABI_DEALLOCATE(smat_inv)
   ABI_DEALLOCATE(smat_k)
   ABI_DEALLOCATE(pwind_k)
   ABI_DEALLOCATE(pwnsfac_k)
   ABI_DEALLOCATE(grad_total)
   if (gs_hamk%usepaw /= 0) then
     call pawcprj_free(cprj_k)
     call pawcprj_free(cprj_kb)
     call pawcprj_free(cprj_direc)
     call pawcprj_free(cprj_band_srt)
     call pawcprj_free(cprj_gat)
     if (nkpt /= dtefield%fnkpt) then
       call pawcprj_free(cprj_fkn)
       call pawcprj_free(cprj_ikn)
       ABI_DATATYPE_DEALLOCATE(cprj_fkn)
       ABI_DATATYPE_DEALLOCATE(cprj_ikn)
     end if
   end if
   ABI_DEALLOCATE(smat_k_paw)
   ABI_DEALLOCATE(dimlmn)
   ABI_DATATYPE_DEALLOCATE(cprj_k)
   ABI_DATATYPE_DEALLOCATE(cprj_kb)
   ABI_DATATYPE_DEALLOCATE(cprj_direc)
   ABI_DATATYPE_DEALLOCATE(cprj_gat)
   ABI_DEALLOCATE(ikptf_recv)
   ABI_DATATYPE_DEALLOCATE(cprj_band_srt)
 end if

! Do not delete this line, needed to run with open MP
 write(unit=message,fmt=*) resid(1)

 call timab(22,2,tsec)

 DBG_EXIT("COLL")

end subroutine cgwf
!!***

!!****f* m_cgwf/linemin
!! NAME
!! linemin
!!
!! FUNCTION
!! Performs the "line minimization" w.r.t. the angle theta on a unit circle
!! to update the wavefunction associated with the current k-point and
!! band label.
!! This routine is used only when the electric field is on (otherwise it could
!! in principle also be used, but there is a simpler procedure, as originally
!! coded in abinit).
!!
!! INPUTS
!! chc = <C|H_0|C> where |C> is the wavefunction of the current band
!! detovc = determinant of the overlap matrix S
!! detovd = determinant of the overlap matrix where for the band
!!          that is being updated <C| is replaced by <D| (search direction)
!! dhc = Re[<D|H_0|C>]
!! dhd = <D|H_0|D>
!! efield_dot = reciprocal lattice coordinates of the electric field
!! iline = index of the current line minimization
!! nkpt = number of k-points
!! nstr(idir) = number of strings along the idir-th direction
!! sdeg = spin degeneracy
!!
!! OUTPUT
!! bcut(ifor,idir) = branch cut of the ellipse associated with (ifor,idir)
!! costh = cos(thetam)
!! hel(ifor,idir) = helicity of the ellipse associated with (ifor,idir)
!! phase_end = total change in Zak phase, must be equal to
!!             dphase_aux1 + n*two_pi
!! sinth = sin(thetam)
!! thetam = optimal angle theta in line_minimization
!!
!!
!! SIDE EFFECTS
!! Input/Output
!! dphase_aux1 = change in Zak phase accumulated during the loop over iline
!!               (can be used for debugging in cgwf.f)
!! phase_init = initial Zak phase (before doing the first line minimization)
!!
!! NOTES
!! We are making the "frozen Hamiltonian approximation", i.e., the
!! Hamiltonian does not change with theta (we are neglecting the dependence
!! of the Hartree and exchange-correlation terms on theta; the original
!! abinit routine does the same)
!!
!! PARENTS
!!      m_cgwf
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

subroutine linemin(bcut,chc,costh,detovc,detovd,dhc,dhd,dphase_aux1,&
&  efield_dot,iline,nkpt,nstr,hel,phase_end,phase_init,sdeg,sinth,thetam)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iline,nkpt
 real(dp),intent(in) :: chc,dhc,dhd,sdeg
 real(dp),intent(out) :: costh,sinth,thetam
!arrays
 integer,intent(in) :: nstr(3)
 integer,intent(out) :: hel(2,3)
 real(dp),intent(in) :: detovc(2,2,3),detovd(2,2,3),efield_dot(3)
 real(dp),intent(inout) :: dphase_aux1(3),phase_init(3)
 real(dp),intent(out) :: bcut(2,3),phase_end(3)

!Local variables -------------------------
!scalars
 integer :: idir,ifor,igrid,iter,maxiter,ngrid
 real(dp) :: aa,angle,bb,big_axis,cc,cphi_0,delta_theta,e0,e1
 real(dp) :: excentr,iab,phase0,phase_min,phi_0,rdum,sgn,small_axis,sphi_0
 real(dp) :: theta,theta_0,val
 logical :: flag_neg
 character(len=500) :: message
!arrays
 real(dp) :: g_theta(2),theta_min(2),theta_try(2)
 real(dp) :: esave(251),e1save(251)   !!REC

! ***********************************************************************

!Compute the helicity and the branch cut of the ellipse in the complex
!plane associated with the overlap between a k-point and one of its neighbours

 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   do ifor = 1, 2

     aa = half*(detovc(1,ifor,idir)*detovc(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovc(2,ifor,idir) + &
&     detovd(1,ifor,idir)*detovd(1,ifor,idir) + &
&     detovd(2,ifor,idir)*detovd(2,ifor,idir))

     bb = half*(detovc(1,ifor,idir)*detovc(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovc(2,ifor,idir) - &
&     detovd(1,ifor,idir)*detovd(1,ifor,idir) - &
&     detovd(2,ifor,idir)*detovd(2,ifor,idir))

     cc = detovc(1,ifor,idir)*detovd(1,ifor,idir) + &
&     detovc(2,ifor,idir)*detovd(2,ifor,idir)

     iab = detovc(1,ifor,idir)*detovd(2,ifor,idir) - &
&     detovc(2,ifor,idir)*detovd(1,ifor,idir)

     if (iab >= zero) then
       hel(ifor,idir) = 1
     else
       hel(ifor,idir) = -1
     end if

     if (abs(bb) > tol8) then
       theta_0 = half*atan(cc/bb)
     else
       theta_0 = quarter*pi
     end if

     if (bb < zero) theta_0 = theta_0 + pi*half

     g_theta(:) = cos(theta_0)*detovc(:,ifor,idir) + &
&     sin(theta_0)*detovd(:,ifor,idir)
!    DEBUG
!    write(std_out,*)'before rhophi, g_theta =',g_theta
!    ENDDEBUG
     call rhophi(g_theta,phi_0,rdum)
!    DEBUG
!    write(std_out,*)'after rhophi, phi_0 = ',phi_0
!    ENDDEBUG

     cphi_0 = cos(phi_0)
     sphi_0 = sin(phi_0)

     rdum = aa - sqrt(bb*bb + cc*cc)
     if (rdum < zero) rdum = zero
     small_axis = sqrt(rdum)
     big_axis = sqrt(aa + sqrt(bb*bb + cc*cc))
     excentr = hel(ifor,idir)*small_axis/big_axis

!    Find angle for which phi = pi
     if (abs(excentr) > tol8) then
       angle = atan(tan(pi-phi_0)/excentr)
     else
       if (tan(pi-phi_0)*hel(ifor,idir) > zero) then
         angle = half*pi
       else
         angle = -0.5_dp*pi
       end if
     end if
     bcut(ifor,idir) = angle + theta_0


!    Compute the branch-cut angle
     if (hel(ifor,idir) == 1) then
       if ((sphi_0 > 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) + pi
       if ((sphi_0 < 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) - pi
     else
       if ((sphi_0 > 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) - pi
       if ((sphi_0 < 0).and.(cphi_0 > 0)) bcut(ifor,idir) = bcut(ifor,idir) + pi
     end if

     if (bcut(ifor,idir) > pi) bcut(ifor,idir) = bcut(ifor,idir) - two_pi
     if (bcut(ifor,idir) < -1_dp*pi) bcut(ifor,idir) = bcut(ifor,idir) + two_pi

!    DEBUG
!    write(std_out,'(a,2x,i3,2x,i3,5x,f16.9,5x,i2)')'linemin: ifor,idir,bcut,hel',&
!    &   ifor,idir,bcut(ifor,idir),hel(ifor,idir)
!    write(std_out,'(a,5x,f16.9,5x,f16.9)')'linemin: big_axis,small_axis ',&
!    &     big_axis,small_axis
!    ENDDEBUG

   end do   ! ifor
 end do   ! idir

!---------------------------------------------------------------------------

!Perform the "line minimization" w.r.t. the angle theta on a unit circle
!to update the wavefunction associated with the current k-point and band label.

 ngrid = 250   ! initial number of subdivisions in [-pi/2,pi/2]
!for finding extrema
 maxiter = 100
 delta_theta = pi/ngrid

!DEBUG
!write(std_out,*)'linemin: theta, e0, e1, e1fdiff'
!ENDDEBUG


!Get the interval where the absolute minimum of E(theta) is located

 val = huge(one)             ! large number
 flag_neg=.false.
 theta_min(:) = ten
 do igrid = 1, ngrid+1

   theta = (igrid - 1)*delta_theta - pi*half
   call etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&   hel,nkpt,nstr,sdeg,theta)

   esave(igrid)=e0      !!REC
   e1save(igrid)=e1     !!REC

!  It is important to detect when the slope changes from negative to positive
!  Moreover, a slope being extremely close to zero must be ignored

!  DEBUG
!  write(std_out,*)' igrid,e0,e1,val,theta_min(:)=',igrid,theta,e0,e1,val,theta_min(:)
!  ENDDEBUG

!  Store e1 and theta if negative ...
   if(e1 < -tol10)then
     theta_try(1)=theta
     flag_neg=.true.
   end if
!  A change of sign is just happening
   if(e1 > tol10 .and. flag_neg)then
     theta_try(2)=theta
     flag_neg=.false.
!    Still, must be better than the previous minimum in order to succeed
     if (e0 < val-tol10) then
       val=e0
       theta_min(:)=theta_try(:)
     end if
   end if
 end do

!In case the minimum was not found

 if (abs(theta_min(1) - ten) < tol10) then
!  REC start
   write(message,'(a,a)')ch10,&
&   ' linemin: ERROR- cannot find theta_min.'
   call wrtout(std_out,message,'COLL')
   write(message,'(a,a)')ch10,&
&   ' igrid      theta          esave(igrid)    e1save(igrid) '
   call wrtout(std_out,message,'COLL')
   do igrid = 1, ngrid+1
     theta = (igrid - 1)*delta_theta - pi*half
     write(std_out,'(i6,3f16.9)')igrid,theta,esave(igrid),e1save(igrid)
     !write(101,'(i6,3f16.9)')igrid,theta,esave(igrid),e1save(igrid)
   end do
   write(message,'(6a)')ch10,&
&   ' linemin: ERROR - ',ch10,&
&   '  Cannot find theta_min. No minimum exists: the field is too strong ! ',ch10,&
&   '  Try decreasing difference between D and 4 Pi P by changing structure or D (only for fixed D calculation)'
   call wrtout(std_out,message,'COLL')
   MSG_ERROR('linemin cannot find theta_min')
 end if

!Compute the mimum of E(theta)


 iter = 0
 do while ((delta_theta > tol8).and.(iter < maxiter))
   delta_theta = half*(theta_min(2) - theta_min(1))
   theta = theta_min(1) + delta_theta
   call etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&   hel,nkpt,nstr,sdeg,theta)
   if (e1 > zero) then
     theta_min(2) = theta
   else
     theta_min(1) = theta
   end if
   iter = iter + 1

!  DEBUG
!  write(std_out,'(a,2x,i3,2(2x,f16.9))')'iter,e0,e1 = ',iter,e0,e1
!  ENDDEBUG

 end do

 costh = cos(theta)
 sinth = sin(theta)

 thetam = theta

!DEBUG
!write(std_out,*)'linemin : thetam = ',thetam
!ENDDEBUG

!---------------------------------------------------------------------------

!Compute and store the change in electronic polarization

 sgn = one
 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   phase_end(idir) = zero
   do ifor = 1, 2

     g_theta(:) = detovc(:,ifor,idir)
!    DEBUG
!    write(std_out,*)'before rhophi (2nd call), g_theta =',g_theta
!    ENDDEBUG
     call rhophi(g_theta,phase0,rdum)
!    DEBUG
!    write(std_out,*)'after rhophi, phase0 = ',phase0
!    ENDDEBUG

     if(iline == 1) phase_init(idir) = phase_init(idir) + sgn*phase0

     g_theta(:) = costh*detovc(:,ifor,idir) + sinth*detovd(:,ifor,idir)
     call rhophi(g_theta,phase_min,rdum)

     phase_end(idir) = phase_end(idir) + sgn*phase_min

!    Correct for branch cuts (remove jumps)
     if (bcut(ifor,idir) <= zero) phase0 = phase0 + hel(ifor,idir)*two_pi
     if(thetam >= bcut(ifor,idir)) phase_min = phase_min + hel(ifor,idir)*two_pi

     dphase_aux1(idir) = dphase_aux1(idir) + sgn*(phase_min - phase0)

     sgn = -1_dp*sgn

   end do   ! idir
 end do    ! ifor

!DEBUG
!write(std_out,'(a,3(2x,f16.9))')'dphase_aux1 = ',(dphase_aux1(idir),idir = 1, 3)
!write(std_out,*)' linemin: debug, exit.'
!ENDDEBUG

end subroutine linemin
!!***

!!****f* m_cgwf/etheta
!! NAME
!! etheta
!!
!! FUNCTION
!! Computes the energy per unit cell and its first derivative
!! for a given angle theta. More precisely, computes only the part of
!! the energy that changes with theta.
!!
!! INPUTS
!! bcut(ifor,idir) = branch cut of the ellipse associated with (ifor,idir)
!! chc = <C|H_0|C> where |C> is the wavefunction of the current band
!! detovc = determinant of the overlap matrix S
!! detovd = determinant of the overlap matrix where for the band
!!          that is being updated <C| is replaced by <D| (search direction)
!! dhc = Re[<D|H_0|C>]
!! dhd = <D|H_0|D>
!! efield_dot = reciprocal lattice coordinates of the electric field
!! hel(ifor,idir) = helicity of the ellipse associated with (ifor,idir)
!! nkpt = number of k-points
!! nsppol = 1 for unpolarized, 2 for spin-polarized
!! nstr(idir) = number of strings along the idir-th direction
!! sdeg = spin degeneracy
!! theta = value of the angle for which the energy (e0) and its
!!         derivative (e1) are computed
!!
!! OUTPUT
!! e0 = energy for the given value of theta
!! e1 = derivative of the energy with respect to theta
!!
!! PARENTS
!!      m_cgwf
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

subroutine etheta(bcut,chc,detovc,detovd,dhc,dhd,efield_dot,e0,e1,&
&    hel,nkpt,nstr,sdeg,theta)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 real(dp),intent(in) :: chc,dhc,dhd,sdeg,theta
 real(dp),intent(out) :: e0,e1
!arrays
 integer,intent(in) :: hel(2,3),nstr(3)
 real(dp),intent(in) :: bcut(2,3),detovc(2,2,3),detovd(2,2,3),efield_dot(3)

!Local variables -------------------------
!scalars
 integer :: idir,ifor
 real(dp) :: c2theta,ctheta,dphase,gnorm,phase,rho,s2theta,sgn,stheta
!arrays
 real(dp) :: dg_theta(2),g_theta(2)

! ***********************************************************************

 e0 = zero ; e1 = zero

 ctheta = cos(theta)
 stheta = sin(theta)
 c2theta = ctheta*ctheta - stheta*stheta   ! cos(2*theta)
 s2theta = two*ctheta*stheta               ! sin(2*theta)

 e0 = chc*ctheta*ctheta + dhd*stheta*stheta + dhc*s2theta
 e0 = e0*sdeg/nkpt

!DEBUG
!e0 = zero
!ENDDEBUG

 e1 = (dhd - chc)*s2theta + two*dhc*c2theta
 e1 = e1*sdeg/nkpt

 sgn = -1_dp
 do idir = 1, 3

   if (abs(efield_dot(idir)) < tol12) cycle

   do ifor = 1, 2

     g_theta(:)  = ctheta*detovc(:,ifor,idir) + &
&     stheta*detovd(:,ifor,idir)
     dg_theta(:) = -1_dp*stheta*detovc(:,ifor,idir) + &
&     ctheta*detovd(:,ifor,idir)

!    Compute E(theta)

     call rhophi(g_theta,phase,rho)
     if (theta >= bcut(ifor,idir)) phase = phase + hel(ifor,idir)*two_pi

!    DEBUG
!    unit = 100 + 10*idir + ifor
!    write(unit,'(4(f16.9))')theta,g_theta(:),phase
!    ENDDEBUG

     e0 = e0 + sgn*sdeg*efield_dot(idir)*phase/(two_pi*nstr(idir))


!    Compute dE/dtheta

!    imaginary part of the derivative of ln(g_theta)
     gnorm = g_theta(1)*g_theta(1) + g_theta(2)*g_theta(2)
     dphase = (dg_theta(2)*g_theta(1) - dg_theta(1)*g_theta(2))/gnorm

     e1 = e1 + sgn*sdeg*efield_dot(idir)*dphase/(two_pi*nstr(idir))

     sgn = -1_dp*sgn

   end do
 end do

end subroutine etheta
!!***

!!****f* m_cgwf/mksubham
!! NAME
!! mksubham
!!
!! FUNCTION
!! Build the Hamiltonian matrix in the eigenfunctions subspace,
!! for one given band (or for one given block of bands)
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!  iblock=index of block of bands
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array cg
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock=number of bands in a block
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  use_subvnlx= 1 if <C band,k|H|C band_prime,k> has to be computed
!!  me_g0=1 if this processors has G=0, 0 otherwise
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ghc(2,npw_k*nspinor)=<G|H|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  gvnlxc(2,npw_k*nspinor)=<G|Vnl|C band,k> for the current state
!!                       This is an input in non-blocked algorithm
!!                               an output in blocked algorithm
!!  isubh=index of current state in array subham
!!  isubo=index of current state in array subovl
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in the WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in the WFs subspace
!!  subvnlx(nband_k*(nband_k+1)*use_subvnlx)=non-local Hamiltonian (if NCPP)  plus Fock ACE operator (if usefock_ACE)
!!   expressed in the WFs subspace
!!
!! PARENTS
!!      m_cgwf
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

subroutine mksubham(cg,ghc,gsc,gvnlxc,iblock,icg,igsc,istwf_k,&
&                    isubh,isubo,mcg,mgsc,nband_k,nbdblock,npw_k,&
&                    nspinor,subham,subovl,subvnlx,use_subovl,use_subvnlx,me_g0)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iblock,icg,igsc,istwf_k,mcg,mgsc,nband_k
 integer,intent(in) :: nbdblock,npw_k,nspinor,use_subovl,use_subvnlx,me_g0
 integer,intent(inout) :: isubh,isubo
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: gsc(2,mgsc)
 real(dp),intent(inout) :: ghc(2,npw_k*nspinor),gvnlxc(2,npw_k*nspinor)
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1))
 real(dp),intent(inout) :: subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(inout) :: subvnlx(nband_k*(nband_k+1)*use_subvnlx)

!Local variables-------------------------------
!scalars
 integer :: iband,ibdblock,ii,ipw,ipw1,isp,iwavef,jwavef
 real(dp) :: cgimipw,cgreipw,chcim,chcre,cscim,cscre,cvcim,cvcre
!real(dp) :: chc(2),cvc(2),csc(2)

! *********************************************************************

!Loop over bands in a block This loop can be parallelized
 do iband=1+(iblock-1)*nbdblock,min(iblock*nbdblock,nband_k)
   ibdblock=iband-(iblock-1)*nbdblock

!  Compute elements of subspace Hamiltonian <C(i)|H|C(n)> and <C(i)|Vnl|C(n)>
   if(istwf_k==1)then

     do ii=1,iband
       iwavef=(ii-1)*npw_k*nspinor+icg
       chcre=zero ; chcim=zero
       if (use_subvnlx==0) then
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
         end do
!        chc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),ghc)
       else
#if 1
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           chcim=chcim+cgreipw*ghc(2,ipw)-cgimipw*ghc(1,ipw)
         end do
         cvcre=zero ; cvcim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           cvcre=cvcre+cgreipw*gvnlxc(1,ipw)+cgimipw*gvnlxc(2,ipw)
           cvcim=cvcim+cgreipw*gvnlxc(2,ipw)-cgimipw*gvnlxc(1,ipw)
         end do
         subvnlx(isubh  )=cvcre
         subvnlx(isubh+1)=cvcim
#else
!        New version with BLAS1, will require some update of the refs.
         cvc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),gvnlxc)
         subvnlx(isubh  )=cvc(1)
         subvnlx(isubh+1)=cvc(2)
         chc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),ghc)
         chcre = chc(1)
         chcim = chc(2)
#endif
!        Store real and imag parts in Hermitian storage mode:
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=chcim
!      subham(isubh  )=chc(1)
!      subham(isubh+1)=chc(2)
       isubh=isubh+2
     end do

   else if(istwf_k>=2)then
     do ii=1,iband
       iwavef=(ii-1)*npw_k+icg
!      Use the time-reversal symmetry, but should not double-count G=0
       if(istwf_k==2 .and. me_g0==1) then
         chcre = half*cg(1,1+iwavef)*ghc(1,1)
         if (use_subvnlx==1) cvcre=half*cg(1,1+iwavef)*gvnlxc(1,1)
         ipw1=2
       else
         chcre=zero; ipw1=1
         if (use_subvnlx==1) cvcre=zero
       end if
       if (use_subvnlx==0) then
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
           end do
         end do
         chcre=two*chcre
       else
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             chcre=chcre+cgreipw*ghc(1,ipw)+cgimipw*ghc(2,ipw)
             cvcre=cvcre+cgreipw*gvnlxc(1,ipw)+cgimipw*gvnlxc(2,ipw)
           end do
         end do
         chcre=two*chcre
         cvcre=two*cvcre
!        Store real and imag parts in Hermitian storage mode:
         subvnlx(isubh  )=cvcre
         subvnlx(isubh+1)=zero
       end if
       subham(isubh  )=chcre
       subham(isubh+1)=zero
       isubh=isubh+2
     end do
   end if

!  Compute elements of subspace <C(i)|S|C(n)> (S=overlap matrix)
!  <C(i)|S|C(n)> should be closed to Identity.
   if (use_subovl==1) then
     jwavef=(iband-1)*npw_k*nspinor+igsc
     if(istwf_k==1)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         cscre=zero ; cscim=zero
         do ipw=1,npw_k*nspinor
           cgreipw=cg(1,ipw+iwavef)
           cgimipw=cg(2,ipw+iwavef)
           cscre=cscre+cgreipw*gsc(1,ipw+jwavef)+cgimipw*gsc(2,ipw+jwavef)
           cscim=cscim+cgreipw*gsc(2,ipw+jwavef)-cgimipw*gsc(1,ipw+jwavef)
         end do
!        csc = cg_zdotc(npw_k*nspinor,cg(1,1+iwavef),gsc)
!        subovl(isubo  )=csc(1)
!        subovl(isubo+1)=csc(2)
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=cscim
         isubo=isubo+2
       end do
     else if(istwf_k>=2)then
       do ii=1,iband
         iwavef=(ii-1)*npw_k*nspinor+icg
         if(istwf_k==2 .and. me_g0==1)then
           cscre=half*cg(1,1+iwavef)*gsc(1,1+jwavef)
           ipw1=2
         else
           cscre=zero; ipw1=1
         end if
         do isp=1,nspinor
           do ipw=ipw1+(isp-1)*npw_k,npw_k*isp
             cgreipw=cg(1,ipw+iwavef)
             cgimipw=cg(2,ipw+iwavef)
             cscre=cscre+cg(1,ipw+iwavef)*gsc(1,ipw+jwavef)+cg(2,ipw+iwavef)*gsc(2,ipw+jwavef)
           end do
         end do
         cscre=two*cscre
!        Store real and imag parts in Hermitian storage mode:
         subovl(isubo  )=cscre
         subovl(isubo+1)=zero
         isubo=isubo+2
       end do
     end if
   end if

 end do ! iband in a block

end subroutine mksubham
!!***

!!****f* ABINIT/make_grad_berry
!! NAME
!! make_grad_berry
!!
!! FUNCTION
!! compute gradient contribution from berry phase in finite
!! electric field case
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)=input wavefunctions
!!  cgq(2,mcgq) = wavefunctions at neighboring k points
!!  cprj_k(natom,nband_k*usepaw)=cprj at this k point
!!  dimlmn(natom)=lmn_size for each atom in input order
!!  dimlmn_srt(natom)=lmn_size for each atom sorted by type
!!  direc(2,npw*nspinor)=gradient vector
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  iband=index of band currently being treated
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point currently being treated
!!  isppol=spin polarization currently treated
!!  natom=number of atoms in cell.
!!  mband =maximum number of bands
!!  mpw=maximum dimensioned size of npw
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array
!!  mkgq = second dimension of pwnsfacq
!!  nkpt=number of k points
!!  mpi_enreg=information about MPI parallelization
!!  npw=number of planewaves in basis sphere at given k.
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=number of spin polarizations
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
!!                     current k-point (electric field, MPI //)
!!
!! OUTPUT
!! grad_berry(2,npw*nspinor) :: contribution to gradient in finite electric field case
!!
!! SIDE EFFECTS
!!  dtefield <type(efield_type)> = variables related to Berry phase calculations (see initberry.f)
!!
!! NOTES
!!
!! PARENTS
!!      m_cgwf
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawcprj_symkn,smatrix,smatrix_k_paw
!!
!! SOURCE

subroutine make_grad_berry(cg,cgq,cprj_k,detovc,dimlmn,dimlmn_srt,direc,dtefield,grad_berry,&
&                          gs_hamk,iband,icg,ikpt,isppol,mband,mcg,mcgq,mkgq,mpi_enreg,mpw,natom,nkpt,&
&                          npw,npwarr,nspinor,nsppol,pwind,pwind_alloc,pwnsfac,pwnsfacq)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iband,icg,ikpt,isppol,mband,mcg,mcgq
  integer,intent(in) :: mkgq,mpw,natom,nkpt,npw,nspinor,nsppol,pwind_alloc
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(efield_type),intent(inout) :: dtefield
  type(MPI_type),intent(in) :: mpi_enreg

  !arrays
  integer,intent(in) :: dimlmn(natom),dimlmn_srt(natom)
  integer,intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: cg(2,mcg),cgq(2,mcgq)
  real(dp),intent(inout) :: direc(2,npw*nspinor)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
  real(dp),intent(out) :: detovc(2,2,3),grad_berry(2,npw*nspinor)
  type(pawcprj_type),intent(in) :: cprj_k(natom,dtefield%mband_occ*gs_hamk%usepaw*dtefield%nspinor)

  !Local variables-------------------------------
  !scalars
  integer :: choice,cpopt,ddkflag,dimenlr1,iatom,icg1,icp2,idum1
  integer :: idir,ifor,ikgf,ikptf,ikpt2,ikpt2f,ipw,i_paw_band,ispinor,itrs,itypat,job
  integer :: klmn,mcg1_k,mcg_q,nbo,npw_k2,nspinortot,paw_opt,shiftbd,signs
  real(dp) :: fac
  character(len=500) :: message
  !arrays
  integer :: pwind_k(npw),sflag_k(dtefield%mband_occ)
  real(dp) :: cg1_k(2,npw*nspinor),dtm_k(2),pwnsfac_k(4,mpw)
  real(dp) :: smat_k(2,dtefield%mband_occ,dtefield%mband_occ)
  real(dp) :: smat_inv(2,dtefield%mband_occ,dtefield%mband_occ),svectout_dum(2,0)
  real(dp) :: dummy_enlout(0)
  real(dp),allocatable :: cgq_k(:,:),enl_rij(:,:,:,:),grad_berry_ev(:,:)
  real(dp),allocatable :: qijbkk(:,:,:,:),smat_k_paw(:,:,:)
  ! type(pawcprj_type) :: cprj_dum(1,1) ! was used in on-site dipole, now suppressed
  ! 15 June 2012 J Zwanziger
  type(pawcprj_type),allocatable :: cprj_kb(:,:),cprj_band_srt(:,:)
  type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)


  ! *********************************************************************

  !DBG_ENTER("COLL")

  nbo = dtefield%mband_occ

  !allocations

  !Electric field: compute the gradient of the Berry phase part of the energy functional.
  !See PRL 89, 117602 (2002) [[cite:Souza2002]], grad_berry(:,:) is the second term of Eq. (4)
  grad_berry(:,:) = zero
  job = 11 ; shiftbd = 1
  mcg_q = mpw*mband*nspinor
  mcg1_k = npw*nspinor

  if (gs_hamk%usepaw /= 0) then
     dimenlr1 = gs_hamk%lmnmax*(gs_hamk%lmnmax+1)/2
     ABI_ALLOCATE(qijbkk,(dimenlr1,natom,nspinor**2,2))
     ABI_ALLOCATE(enl_rij,(nspinor*dimenlr1,natom,nspinor**2,1))
     ABI_ALLOCATE(smat_k_paw,(2,nbo,nbo))
     ABI_ALLOCATE(grad_berry_ev,(2,npw*nspinor))
     enl_rij = zero
     qijbkk = zero
     smat_k_paw = zero
     ABI_DATATYPE_ALLOCATE(cprj_kb,(natom,nbo*nspinor))
     call pawcprj_alloc(cprj_kb,0,dimlmn)
     ABI_DATATYPE_ALLOCATE(cprj_band_srt,(natom,nspinor))
     call pawcprj_alloc(cprj_band_srt,0,dimlmn_srt)
     if (nkpt /= dtefield%fnkpt) then
        ABI_DATATYPE_ALLOCATE(cprj_fkn,(natom,nbo*nspinor))
        ABI_DATATYPE_ALLOCATE(cprj_ikn,(natom,nbo*nspinor))
        call pawcprj_alloc(cprj_fkn,0,dimlmn)
        call pawcprj_alloc(cprj_ikn,0,dimlmn)
     else
        ABI_DATATYPE_ALLOCATE(cprj_fkn,(0,0))
        ABI_DATATYPE_ALLOCATE(cprj_ikn,(0,0))
     end if
  else
     ABI_ALLOCATE(qijbkk,(0,0,0,0))
     ABI_ALLOCATE(enl_rij,(0,0,0,0))
     ABI_ALLOCATE(smat_k_paw,(0,0,0))
     ABI_ALLOCATE(grad_berry_ev,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_kb,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_band_srt,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(0,0))
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(0,0))
  end if

  ikptf = dtefield%i2fbz(ikpt)
  ikgf = dtefield%fkgindex(ikptf)  ! this is the shift for pwind

  do idir = 1, 3
     !  skip idir values for which efield_dot(idir)=0
     if (abs(dtefield%efield_dot(idir)) < tol12) cycle
     !  Implicitly, we use the gradient multiplied by the number of k points in the FBZ
     fac = dtefield%efield_dot(idir)*dble(dtefield%fnkpt)/&
          &   (dble(dtefield%nstr(idir))*four_pi)
     do ifor = 1, 2
        !    Handle dtefield%i2fbz properly and ask whether t.r.s. is used
        ikpt2f = dtefield%ikpt_dk(ikptf,ifor,idir)
        if (dtefield%indkk_f2ibz(ikpt2f,6) == 1) then
           itrs = 10
        else
           itrs = 0
        end if
        ikpt2 = dtefield%indkk_f2ibz(ikpt2f,1)
        npw_k2 = npwarr(ikpt2)
        ABI_ALLOCATE(cgq_k,(2,nbo*nspinor*npw_k2))
        pwind_k(1:npw) = pwind(ikgf+1:ikgf+npw,ifor,idir)
        pwnsfac_k(1:2,1:npw) = pwnsfac(1:2,ikgf+1:ikgf+npw)
        sflag_k(:) = dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir)
        smat_k(:,:,:) = dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir)
        if (mpi_enreg%nproc_cell > 1) then
           icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
           cgq_k(:,1:nbo*nspinor*npw_k2) = &
                &       cgq(:,icg1+1:icg1+nbo*nspinor*npw_k2)
           idum1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*nkpt)
           pwnsfac_k(3:4,1:npw_k2) = pwnsfacq(1:2,idum1+1:idum1+npw_k2)
        else
           icg1 = dtefield%cgindex(ikpt2,isppol)
           cgq_k(:,1:nbo*nspinor*npw_k2) = &
                &       cg(:,icg1+1:icg1+nbo*nspinor*npw_k2)
           idum1 = dtefield%fkgindex(ikpt2f)
           pwnsfac_k(3:4,1:npw_k2) = pwnsfac(1:2,idum1+1:idum1+npw_k2)
        end if
        if (gs_hamk%usepaw == 1) then
           icp2=nbo*(ikpt2-1)*nspinor
           call pawcprj_get(gs_hamk%atindx1,cprj_kb,dtefield%cprj,natom,1,icp2,ikpt,0,isppol,&
                &       nbo,dtefield%fnkpt,natom,nbo,nbo,nspinor,nsppol,0,&
                &       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           if (ikpt2 /= ikpt2f) then ! construct cprj_kb by symmetry
              call pawcprj_copy(cprj_kb,cprj_ikn)
              call pawcprj_symkn(cprj_fkn,cprj_ikn,dtefield%atom_indsym,dimlmn,-1,gs_hamk%indlmn,&
                   &         dtefield%indkk_f2ibz(ikpt2f,2),dtefield%indkk_f2ibz(ikpt2f,6),&
                   &         dtefield%fkptns(:,dtefield%i2fbz(ikpt2)),&
                   &         dtefield%lmax,dtefield%lmnmax,mband,natom,nbo,nspinor,&
                   &         dtefield%nsym,gs_hamk%ntypat,gs_hamk%typat,dtefield%zarot)
              call pawcprj_copy(cprj_fkn,cprj_kb)
           end if
           call smatrix_k_paw(cprj_k,cprj_kb,dtefield,idir,ifor,mband,natom,smat_k_paw,gs_hamk%typat)
        end if

        icg1 = 0 ; ddkflag = 1
        call smatrix(cg,cgq_k,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,&
             &     job,iband,mcg,mcg_q,mcg1_k,iband,mpw,nbo,dtefield%nband_occ(isppol),&
             &     npw,npw_k2,nspinor,pwind_k,pwnsfac_k,sflag_k,&
             &     shiftbd,smat_inv,smat_k,smat_k_paw,gs_hamk%usepaw)
        ABI_DEALLOCATE(cgq_k)
        detovc(:,ifor,idir) = dtm_k(:) !store the determinant of the overlap
        if (sqrt(dtm_k(1)*dtm_k(1) + dtm_k(2)*dtm_k(2)) < tol12) then
           write(message,'(3a,i5,a,i3,a,a,a)') &
                &       '  (electric field)',ch10,&
                &       '  For k-point #',ikpt,' and band # ',iband,',',ch10,&
                &       '  the determinant of the overlap matrix is found to be 0. Fixing...'
           !      REC try this:
           write(std_out,*)message,dtm_k(1:2)
           if(abs(dtm_k(1))<=1d-12)dtm_k(1)=1d-12
           if(abs(dtm_k(2))<=1d-12)dtm_k(2)=1d-12
           write(std_out,*)' Changing to:',dtm_k(1:2)
           !      REC       MSG_BUG(message)
        end if

        if (gs_hamk%usepaw == 1) then
           !      this loop applies discretized derivative of projectors
           !      note that qijb_kk is sorted by input atom order, but nonlop wants it sorted by type
           do iatom = 1, natom
              itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
              do klmn = 1, dtefield%lmn2_size(itypat)
                 !          note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the qijb is diagonal
                 !          in spin space so only the first two are nonzero and they are equal
                 do ispinor = 1, nspinor
                    qijbkk(klmn,iatom,ispinor,1) = dtefield%qijb_kk(1,klmn,gs_hamk%atindx1(iatom),idir)
                    qijbkk(klmn,  iatom,ispinor,2) = dtefield%qijb_kk(2,klmn,gs_hamk%atindx1(iatom),idir)
                    if (ifor > 1) qijbkk(klmn,iatom,ispinor,2) = -qijbkk(klmn,iatom,ispinor,2)
                 end do
              end do ! end loop over lmn2_size
           end do ! end loop over natom

           choice = 1
           signs = 2
           paw_opt = 1
           cpopt = 2 ! use cprj_kb in memory
           nspinortot=min(2,nspinor*(1+mpi_enreg%paral_spinor))
           do i_paw_band = 1, nbo

              call pawcprj_get(gs_hamk%atindx,cprj_band_srt,cprj_kb,natom,i_paw_band,0,ikpt,1,&
                   &         isppol,nbo,1,natom,1,nbo,nspinor,nsppol,0,&
                   &         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

              ! Pass dummy_enlout to avoid aliasing (enl, enlout)
              call nonlop(choice,cpopt,cprj_band_srt,dummy_enlout,gs_hamk,idir,(/zero/),mpi_enreg,1,0,&
                   &         paw_opt,signs,svectout_dum,0,direc,grad_berry_ev,enl=qijbkk)

              !        Add i*fac*smat_inv(i_paw_band,iband)*grad_berry_ev to the gradient
              do ipw = 1, npw*nspinor

                 grad_berry(1,ipw) = grad_berry(1,ipw) - &
                      &           fac*(smat_inv(2,i_paw_band,iband)*grad_berry_ev(1,ipw) + &
                      &           smat_inv(1,i_paw_band,iband)*grad_berry_ev(2,ipw))

                 grad_berry(2,ipw) = grad_berry(2,ipw) + &
                      &           fac*(smat_inv(1,i_paw_band,iband)*grad_berry_ev(1,ipw) - &
                      &           smat_inv(2,i_paw_band,iband)*grad_berry_ev(2,ipw))

              end do
           end do
        end if ! end if PAW

        !    Add i*fac*cg1_k to the gradient
        do ipw = 1, npw*nspinor
           grad_berry(1,ipw) = grad_berry(1,ipw) - fac*cg1_k(2,ipw)
           grad_berry(2,ipw) = grad_berry(2,ipw) + fac*cg1_k(1,ipw)
        end do
        fac = -1._dp*fac
        dtefield%sflag(:,ikpt+(isppol-1)*nkpt,ifor,idir) = sflag_k(:)
        dtefield%sflag(iband,ikpt+(isppol-1)*nkpt,ifor,idir) = 0
        dtefield%smat(:,:,:,ikpt+(isppol-1)*nkpt,ifor,idir) = smat_k(:,:,:)
     end do  ! ifor

     !  if (gs_hamk%usepaw == 1) then
     !  !    call nonlop to apply on-site dipole <EV> part to direc
     !  !    note that rij is sorted by input atom order, but nonlop wants it sorted by type
     !  do iatom = 1, natom
     !  itypat = gs_hamk%typat(gs_hamk%atindx1(iatom))
     !  do klmn = 1, dtefield%lmn2_size(itypat)
     !  !        note: D_ij-like terms have 4 spinor components: 11, 22, 12, and 21. Here the enl_rij is diagonal
     !  !        in spin space so only the first two are nonzero and they are equal
     !  do ispinor = 1, nspinor
     !  if (nspinor == 1) then
     !  enl_rij(klmn,iatom,ispinor) = dtefield%rij(klmn,itypat,idir)
     !  else
     !  enl_rij(2*klmn-1,iatom,ispinor) = dtefield%rij(klmn,itypat,idir)
     !  end if
     !  end do
     !  end do ! end loop over lmn2_size
     !  end do ! end loop over natom
     !  cpopt = -1 ! compute cprj inside nonlop because we do not have them for direc
     !  call nonlop(choice,cpopt,cprj_dum,dummy_enlout,gs_hamk,idir,zero,mpi_enreg,1,0,&
     !  &           paw_opt,signs,svectout_dum,0,direc,grad_berry_ev,enl=enl_rij)
     !  grad_berry(:,:) = grad_berry(:,:) - dtefield%efield_dot(idir)*grad_berry_ev(:,:)/two_pi
     !  end if

  end do    ! idir

  !deallocations
  if(gs_hamk%usepaw /= 0) then
     call pawcprj_free(cprj_kb)
     call pawcprj_free(cprj_band_srt)
     if (nkpt /= dtefield%fnkpt) then
        call pawcprj_free(cprj_fkn)
        call pawcprj_free(cprj_ikn)
     end if
  end if
  ABI_DEALLOCATE(grad_berry_ev)
  ABI_DEALLOCATE(qijbkk)
  ABI_DEALLOCATE(enl_rij)
  ABI_DEALLOCATE(smat_k_paw)
  ABI_DATATYPE_DEALLOCATE(cprj_kb)
  ABI_DATATYPE_DEALLOCATE(cprj_band_srt)
  ABI_DATATYPE_DEALLOCATE(cprj_fkn)
  ABI_DATATYPE_DEALLOCATE(cprj_ikn)

  !DBG_EXIT("COLL")

end subroutine make_grad_berry
!!***


end module m_cgwf
!!***
