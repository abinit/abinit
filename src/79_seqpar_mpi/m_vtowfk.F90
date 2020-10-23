!!****m* ABINIT/m_vtowfk
!! NAME
!!  m_vtowfk
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MT)
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

module m_vtowfk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_efield
 use m_linalg_interfaces
 use m_cgtools
 use m_dtset
 use m_dtfil

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab, cwtime, sec2str
 use m_fstrings,    only : sjoin, itoa, ftoa
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_paw_dmft,    only : paw_dmft_type
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_put,pawcprj_copy
 use m_paw_dmft,    only : paw_dmft_type
 use m_gwls_hamiltonian, only : build_H
 use m_fftcore,     only : fftcore_set_mixprec, fftcore_mixprec
 use m_cgwf,        only : cgwf
 use m_cgwf_paw,    only : cgwf_paw
 use m_lobpcgwf_old,only : lobpcgwf
 use m_lobpcgwf,    only : lobpcgwf2
 use m_spacepar,    only : meanvalue_g
 use m_chebfi,      only : chebfi
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_nonlop, prep_fourwf
 use m_fft,         only : fourwf

 implicit none

 private
!!***

 public :: vtowfk
!!***

contains
!!***

!!****f* ABINIT/vtowfk
!! NAME
!! vtowfk
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed Hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point
!!
!! INPUTS
!!  cgq = array that holds the WF of the nearest neighbours of
!!        the current k-point (electric field, MPI //)
!!  cpus= cpu time limit in seconds
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fixed_occ=true if electronic occupations are fixed (occopt<3)
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cprj
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=number of the k-point
!!  iscf=(<= 0  =>non-SCF), >0 => SCF
!!  isppol isppol=1 for unpolarized, 2 for spin-polarized
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw_k)=(modified) kinetic energy for each plane wave (Hartree)
!!  mcg=second dimension of the cg array
!!  mcgq=second dimension of the cgq array (electric field, MPI //)
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mkgq = second dimension of pwnsfacq
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nkpt=number of k points.
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!             (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  optforces=option for the computation of forces
!!  prtvol=control print volume and debugging output
!!  pwind(pwind_alloc,2,3)= array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc= first dimension of pwind
!!  pwnsfac(2,pwind_alloc)= phase factors for non-symmorphic translations
!!                          (see initberry.f)
!!  pwnsfacq(2,mkgq)= phase factors for the nearest neighbours of the
!!                    current k-point (electric field, MPI //)
!!  usebanfft=flag for band-fft parallelism
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  wtk=weight assigned to the k point.
!!  zshift(nband_k)=energy shifts for the squared shifted hamiltonian algorithm
!!
!! OUTPUT
!!  dphase_k(3)=change in Zak phase for the current k-point
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  ek_k(nband_k)=contribution from each band to kinetic energy, at this k-point
!!  ek_k_nd(2,nband_k,nband_k*use_dmft)=contribution to kinetic energy, including non-diagonal terms, at this k-point (usefull if use_dmft)
!!  resid_k(nband_k)=residuals for each band over all k points,
!!                   BEFORE the band rotation.
!!  ==== if optforces>0 ====
!!    grnl_k(3*natom,nband_k)=nonlocal gradients, at this k-point
!!  ==== if (gs_hamk%usepaw==0) ====
!!    enlx_k(nband_k)=contribution from each band to nonlocal pseudopotential + potential Fock-type part of total energy, at this k-point
!!  ==== if (gs_hamk%usepaw==1) ====
!!    cprj(natom,mcprj*usecprj)= wave functions projected with non-local projectors:
!!                               cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!  rhoaug(n4,n5,n6,nvloc)= density in electrons/bohr**3, on the augmented fft grid.
!!                    (cumulative, so input as well as output). Update only
!!                    for occopt<3 (fixed occupation numbers)
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      build_h,cgwf,chebfi,dsymm,fourwf,fxphas,lobpcgwf,lobpcgwf2,meanvalue_g
!!      nonlop,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_put,prep_fourwf
!!      prep_nonlop,pw_orthon,subdiago,timab,wrtout,xmpi_sum,zhemm
!!
!! NOTES
!!  The cprj are distributed over band and spinors processors.
!!  One processor doesn't know all the cprj.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
!!  are stored on each proc.
!!
!! SOURCE

subroutine vtowfk(cg,cgq,cprj,cpus,dphase_k,dtefield,dtfil,dtset,&
& eig_k,ek_k,ek_k_nd,enlx_k,fixed_occ,grnl_k,gs_hamk,&
& ibg,icg,ikpt,iscf,isppol,kg_k,kinpw,mband_cprj,mcg,mcgq,mcprj,mkgq,mpi_enreg,&
& mpw,natom,nband_k,nkpt,nnsclo_now,npw_k,npwarr,occ_k,optforces,prtvol,&
& pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,rhoaug,paw_dmft,wtk,zshift)

!Arguments ------------------------------------
 integer, intent(in) :: ibg,icg,ikpt,iscf,isppol,mband_cprj,mcg,mcgq,mcprj,mkgq,mpw
 integer, intent(in) :: natom,nband_k,nkpt,nnsclo_now,npw_k,optforces
 integer, intent(in) :: prtvol,pwind_alloc
 logical,intent(in) :: fixed_occ
 real(dp), intent(in) :: cpus,wtk
 type(datafiles_type), intent(in) :: dtfil
 type(efield_type), intent(inout) :: dtefield
 type(dataset_type), intent(in) :: dtset
 type(gs_hamiltonian_type), intent(inout) :: gs_hamk
 type(MPI_type), intent(inout) :: mpi_enreg
 type(paw_dmft_type), intent(in)  :: paw_dmft
 integer, intent(in) :: kg_k(3,npw_k)
 integer, intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp), intent(in) :: cgq(2,mcgq),kinpw(npw_k),occ_k(nband_k)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp), intent(in) :: zshift(nband_k)
 real(dp), intent(out) :: eig_k(nband_k),ek_k(nband_k),dphase_k(3),ek_k_nd(2,nband_k,nband_k*paw_dmft%use_dmft)
 real(dp), intent(out) :: enlx_k(nband_k)
 real(dp), intent(out) :: grnl_k(3*natom,nband_k*optforces)
 real(dp), intent(out) :: resid_k(nband_k)
 real(dp), intent(inout) :: cg(2,mcg),rhoaug(gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,gs_hamk%nvloc)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*gs_hamk%usecprj)

!Local variables-------------------------------
 logical :: has_fock,newlobpcg,enable_cgwf_paw
 integer,parameter :: level=112,tim_fourwf=2,tim_nonlop_prep=11,enough=3
 integer,save :: nskip=0
!     Flag use_subovl: 1 if "subovl" array is computed (see below)
!     subovl should be Identity (in that case we should use use_subovl=0)
!     But this is true only if conjugate gradient algo. converges
 integer :: use_subovl=0
 integer :: use_subvnlx=0
 integer :: use_totvnlx=0
 integer :: bandpp_cprj,blocksize,choice,cpopt,iband,iband1
 integer :: iblock,iblocksize,ibs,idir,ierr,igs,igsc,ii,pidx,inonsc
 integer :: iorder_cprj,ipw,ispinor,ispinor_index,istwf_k,iwavef,jj,mgsc,my_nspinor,n1,n2,n3 !kk
 integer :: nband_k_cprj,nblockbd,ncpgr,ndat,nkpt_max,nnlout,ortalgo
 integer :: paw_opt,quit,signs,spaceComm,tim_nonlop,wfoptalg,wfopta10
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 real(dp) :: ar,ar_im,eshift,occblock
 real(dp) :: res,residk,weight,cpu,wall,gflops
 character(len=500) :: message
 real(dp) :: dummy(2,1),nonlop_dum(1,1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef1(:,:),cwavef_x(:,:),cwavef_y(:,:),cwavefb(:,:,:)
 real(dp),allocatable :: eig_save(:),enlout(:),evec(:,:),evec_loc(:,:),gsc(:,:)
 real(dp),allocatable :: mat_loc(:,:),mat1(:,:,:),matvnl(:,:,:)
 real(dp),allocatable :: subham(:),subovl(:),subvnlx(:),totvnlx(:,:),wfraug(:,:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! **********************************************************************

 DBG_ENTER("COLL")

 call timab(28,1,tsec) ! Keep track of total time spent in "vtowfk"

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'vtowfk: enter'
   call wrtout(std_out,message,'PERS')
 end if

!=========================================================================
!============= INITIALIZATIONS AND ALLOCATIONS ===========================
!=========================================================================

 nkpt_max=50; if(xmpi_paral==1)nkpt_max=-1

 wfoptalg=mod(dtset%wfoptalg,100); wfopta10=mod(wfoptalg,10)
 newlobpcg = (dtset%wfoptalg == 114 .and. dtset%use_gpu_cuda == 0)
 istwf_k=gs_hamk%istwf_k
 has_fock=(associated(gs_hamk%fockcommon))
 quit=0

!Parallelization over spinors management
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 if (mpi_enreg%paral_spinor==0) then
   ispinor_index=1
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(dtset%nspinor==2)
 else
   ispinor_index=mpi_enreg%me_spinor+1
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

!Parallelism over FFT and/or bands: define sizes and tabs
 !if (mpi_enreg%paral_kgb==1) then
 nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
 !else
 !  nblockbd=nband_k/mpi_enreg%nproc_fft
 !  if (nband_k/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
 !end if
 blocksize=nband_k/nblockbd

!Save eshift
 if(wfoptalg==3)then
   eshift=zshift(1)
   ABI_ALLOCATE(eig_save,(nband_k))
   eig_save(:)=eshift
 end if

 n1=gs_hamk%ngfft(1); n2=gs_hamk%ngfft(2); n3=gs_hamk%ngfft(3)

 enable_cgwf_paw = (gs_hamk%usepaw==1).and.(wfopta10==0)
 mgsc=0
 if ( .not. newlobpcg .and. .not. enable_cgwf_paw) then
   igsc=0
   mgsc=nband_k*npw_k*my_nspinor*gs_hamk%usepaw

   ABI_MALLOC_OR_DIE(gsc,(2,mgsc), ierr)
   gsc=zero
 else
   ABI_ALLOCATE(gsc,(0,0))
 end if

 if(wfopta10 /= 1 .and. .not. newlobpcg ) then
   !chebfi already does this stuff inside
   ABI_ALLOCATE(evec,(2*nband_k,nband_k))
   ABI_ALLOCATE(subham,(nband_k*(nband_k+1)))

   ABI_ALLOCATE(subvnlx,(0))
   ABI_ALLOCATE(totvnlx,(0,0))
   if (wfopta10==4) then
!    Later, will have to generalize to Fock case, like when wfopta10/=4
     if (gs_hamk%usepaw==0) then
       ABI_DEALLOCATE(totvnlx)
       if (istwf_k==1) then
         ABI_ALLOCATE(totvnlx,(2*nband_k,nband_k))
       else if (istwf_k==2) then
         ABI_ALLOCATE(totvnlx,(nband_k,nband_k))
       end if
       use_totvnlx=1
     endif
   else
     if (gs_hamk%usepaw==0 .or. has_fock) then
       ABI_DEALLOCATE(subvnlx)
       ABI_ALLOCATE(subvnlx,(nband_k*(nband_k+1)))
       use_subvnlx=1
     end if
   end if

   if (use_subovl==1) then
     ABI_ALLOCATE(subovl,(nband_k*(nband_k+1)))
   else
     ABI_ALLOCATE(subovl,(0))
   end if
 end if

 ! Carry out UP TO dtset%nline steps, or until resid for every band is < dtset%tolwfr
 if (prtvol>2 .or. ikpt<=nkpt_max) then
   write(message,'(a,i5,2x,a,3f9.5,2x,a)')' non-scf iterations; kpt # ',ikpt,', k= (',gs_hamk%kpt_k,'), band residuals:'
   call wrtout(std_out,message,'PERS')
 end if

!Electric field: initialize dphase_k
 dphase_k(:) = zero

!=========================================================================
!==================== NON-SELF-CONSISTENT LOOP ===========================
!=========================================================================

!nnsclo_now=number of non-self-consistent loops for the current vtrial
!(often 1 for SCF calculation, =nstep for non-SCF calculations)
 call timab(39,1,tsec) ! "vtowfk (loop)"

 do inonsc=1,nnsclo_now
   if (iscf < 0 .and. (inonsc <= enough .or. mod(inonsc, 10) == 0)) call cwtime(cpu, wall, gflops, "start")

   ! This initialisation is needed for the MPI-parallelisation (gathering using sum)
   if(wfopta10 /= 1 .and. .not. newlobpcg) then
     subham(:)=zero
     if (gs_hamk%usepaw==0) then
       if (wfopta10==4) then
         totvnlx(:,:)=zero
       else
         subvnlx(:)=zero
       end if
     end if
     if (use_subovl==1)subovl(:)=zero
   end if
   resid_k(:)=zero

!  Filter the WFs when modified kinetic energy is too large (see routine mkkin.f)
!  !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(igs,iwavef)
   do ispinor=1,my_nspinor
     do iband=1,nband_k
       igs=(ispinor-1)*npw_k
       iwavef=(iband-1)*npw_k*my_nspinor+icg
       do ipw=1+igs,npw_k+igs
         if(kinpw(ipw-igs)>huge(zero)*1.d-11)then
           cg(1,ipw+iwavef)=zero
           cg(2,ipw+iwavef)=zero
         end if
       end do
     end do
   end do

   ! JLJ 17/10/2014: If it is a GWLS calculation, construct the hamiltonian
   ! as in a usual GS calc., but skip any minimisation procedure.
   ! This would be equivalent to nstep=0, if the latter did work.
   if(dtset%optdriver/=RUNL_GWLS) then

     if(wfopta10==4.or.wfopta10==1) then
       if (istwf_k/=1.and.istwf_k/=2) then !no way to use lobpcg
         write(message,'(3a)')&
&         'Only istwfk=1 or 2 are allowed with wfoptalg=4/14 !',ch10,&
&         'Action: put istwfk to 1 or remove k points with half integer coordinates.'
         MSG_ERROR(message)
       end if

!    =========================================================================
!    ============ MINIMIZATION OF BANDS: LOBPCG ==============================
!    =========================================================================
       if (wfopta10==4) then
         if ( .not. newlobpcg ) then
           call lobpcgwf(cg,dtset,gs_hamk,gsc,icg,igsc,kinpw,mcg,mgsc,mpi_enreg,&
&           nband_k,nblockbd,npw_k,prtvol,resid_k,subham,totvnlx,use_totvnlx)
!          In case of FFT parallelism, exchange subspace arrays
           spaceComm=mpi_enreg%comm_bandspinorfft
           call xmpi_sum(subham,spaceComm,ierr)
           if (gs_hamk%usepaw==0) then
             if (wfopta10==4) then
               call xmpi_sum(totvnlx,spaceComm,ierr)
             else
               call xmpi_sum(subvnlx,spaceComm,ierr)
             end if
           end if
           if (use_subovl==1) call xmpi_sum(subovl,spaceComm,ierr)
         else
           call lobpcgwf2(cg(:,icg+1:),dtset,eig_k,enlx_k,gs_hamk,kinpw,mpi_enreg,&
&           nband_k,npw_k,my_nspinor,prtvol,resid_k)
         end if
!        In case of FFT parallelism, exchange subspace arrays

!    =========================================================================
!    ============ MINIMIZATION OF BANDS: CHEBYSHEV FILTERING =================
!    =========================================================================
       else if (wfopta10 == 1) then
         call chebfi(cg(:, icg+1:),dtset,eig_k,enlx_k,gs_hamk,gsc,kinpw,&
&         mpi_enreg,nband_k,npw_k,my_nspinor,prtvol,resid_k)
       end if

!      =========================================================================
!      ======== MINIMIZATION OF BANDS: CONJUGATE GRADIENT (Teter et al.) =======
!      =========================================================================
     else
!      use_subvnlx=0; if (gs_hamk%usepaw==0 .or. associated(gs_hamk%fockcommon)) use_subvnlx=1
!      use_subvnlx=0; if (gs_hamk%usepaw==0) use_subvnlx=1
       if (enable_cgwf_paw) then
         call cgwf_paw(cg,dtset%chkexit,cpus,eig_k,dtfil%filnam_ds(1),&
&         gs_hamk,icg,ikpt,inonsc,isppol,mcg,&
&         mpi_enreg,nband_k,dtset%nbdblock,nkpt,dtset%nline,npw_k,my_nspinor,&
&         dtset%nsppol,dtset%ortalg,prtvol,quit,resid_k,subham,dtset%tolrde,dtset%tolwfr,wfoptalg)
       else
         call cgwf(dtset%berryopt,cg,cgq,dtset%chkexit,cpus,dphase_k,dtefield,dtfil%filnam_ds(1),&
&         gsc,gs_hamk,icg,igsc,ikpt,inonsc,isppol,dtset%mband,mcg,mcgq,mgsc,mkgq,&
&         mpi_enreg,mpw,nband_k,dtset%nbdblock,nkpt,dtset%nline,npw_k,npwarr,my_nspinor,&
&         dtset%nsppol,dtset%ortalg,prtvol,pwind,pwind_alloc,pwnsfac,pwnsfacq,quit,resid_k,&
&         subham,subovl,subvnlx,dtset%tolrde,dtset%tolwfr,use_subovl,use_subvnlx,wfoptalg,zshift)
       end if
     end if
   end if

!  =========================================================================
!  ===================== FIND LARGEST RESIDUAL =============================
!  =========================================================================

!  Find largest resid over bands at this k point
!  Note that this operation is done BEFORE rotation of bands:
!  it would be time-consuming to recompute the residuals after.
   residk=maxval(resid_k(1:max(1,nband_k-dtset%nbdbuf)))

!  Print residuals
   if(prtvol>2 .or. ikpt<=nkpt_max)then
     do ii=0,(nband_k-1)/8
       write(message,'(a,8es10.2)')' res:',(resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
       call wrtout(std_out,message,'PERS')
     end do
   end if

!  =========================================================================
!  ========== DIAGONALIZATION OF HAMILTONIAN IN WFs SUBSPACE ===============
!  =========================================================================

   if( .not. wfopta10 == 1 .and. .not. newlobpcg ) then
     call timab(585,1,tsec) !"vtowfk(subdiago)"
     call subdiago(cg,eig_k,evec,gsc,icg,igsc,istwf_k,&
&     mcg,mgsc,nband_k,npw_k,my_nspinor,dtset%paral_kgb,&
&     subham,subovl,use_subovl,gs_hamk%usepaw,mpi_enreg%me_g0)
     call timab(585,2,tsec)
   end if

!  Print energies
   if(prtvol>2 .or. ikpt<=nkpt_max)then
     do ii=0,(nband_k-1)/8
       write(message, '(a,8es10.2)' )' ene:',(eig_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
       call wrtout(std_out,message,'PERS')
     end do
   end if

!  THIS CHANGE OF SHIFT DOES NOT WORK WELL
!  Update zshift in the case of wfoptalg==3
!  if(wfoptalg==3 .and. inonsc/=1)then
!  do iband=1,nband_k
!  if(eig_k(iband)<eshift .and. eig_save(iband)<eshift)then
!  zshift(iband)=max(eig_k(iband),eig_save(iband))
!  end if
!  if(eig_k(iband)>eshift .and. eig_save(iband)>eshift)then
!  zshift(iband)=min(eig_k(iband),eig_save(iband))
!  end if
!  end do
!  eig_save(:)=eig_k(:)
!  end if

!  =========================================================================
!  =============== ORTHOGONALIZATION OF WFs (if needed) ====================
!  =========================================================================

!  Re-orthonormalize the wavefunctions at this k point--
!  this is redundant but is performed to combat rounding error in wavefunction orthogonality

   call timab(583,1,tsec) ! "vtowfk(pw_orthon)"
   ortalgo=mpi_enreg%paral_kgb
   if ((wfoptalg/=14 .and. wfoptalg /= 1 .and. .not. enable_cgwf_paw).or.dtset%ortalg>0) then
     call pw_orthon(icg,igsc,istwf_k,mcg,mgsc,npw_k*my_nspinor,nband_k,ortalgo,gsc,gs_hamk%usepaw,cg,&
&     mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
   end if
   call timab(583,2,tsec)

!  DEBUG seq==par comment next block
!  Fix phases of all bands
   if (xmpi_paral/=1 .or. mpi_enreg%paral_kgb/=1) then
     if ( .not. newlobpcg .and. .not. enable_cgwf_paw ) then
       call fxphas(cg,gsc,icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nband_k,npw_k*my_nspinor,gs_hamk%usepaw)
     else
       ! GSC is local to vtowfk and is completely useless since everything
       ! is calcultated in my lobpcg, we don't care about the phase of gsc !
       call fxphas(cg,gsc,icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nband_k,npw_k*my_nspinor,0)
     end if
   end if

   if (iscf < 0) then
     if (residk > dtset%tolwfr .and. residk < tol7) then
       if (fftcore_mixprec == 1) call wrtout(std_out, " Approaching NSCF convergence. Activating FFT in double-precision")
       ii = fftcore_set_mixprec(0)
     end if

     ! Print residual and wall-time required by NSCF iteration.
     if (inonsc <= enough .or. mod(inonsc, 20) == 0) then
       call cwtime(cpu, wall, gflops, "stop")
       if (inonsc == 1) call wrtout(std_out, sjoin(" k-point: [", itoa(ikpt), "/", itoa(nkpt), "], spin:", itoa(isppol)))
       call wrtout(std_out, sjoin("   Max resid =", ftoa(residk, fmt="es13.5"), &
         " (exclude nbdbuf bands). One NSCF iteration cpu-time:", &
         sec2str(cpu), ", wall-time:", sec2str(wall)), do_flush=.True.)
       if (inonsc == enough) call wrtout(std_out, "   Printing residuals every mod(20) iteration ...")
     end if
   end if

   ! Exit loop over inonsc if converged
   if (residk < dtset%tolwfr) then
     if (iscf < 0) call wrtout(std_out, sjoin("   NSCF loop completed after", itoa(inonsc), "iterations"))
     exit
   end if
 end do ! End loop over inonsc (NON SELF-CONSISTENT LOOP)

 call timab(39,2,tsec)
 call timab(30,1,tsec) ! "vtowfk  (afterloop)"

!###################################################################

!Compute kinetic energy and non-local energy for each band, and in the SCF
!case, contribution to forces, and eventually accumulate rhoaug

 ndat=1;if (mpi_enreg%paral_kgb==1) ndat=mpi_enreg%bandpp
 if(iscf>0 .and. fixed_occ)  then
   ABI_ALLOCATE(wfraug,(2,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6*ndat))
 end if

!"nonlop" routine input parameters
 nnlout=3*natom*optforces
 signs=1;idir=0
 if (gs_hamk%usepaw==0) then
   choice=1+optforces
   paw_opt=0;cpopt=-1;tim_nonlop=2
 else
   choice=2*optforces
   paw_opt=2;cpopt=0;tim_nonlop=10-8*optforces
   if (dtset%usefock==1) then
!     if (dtset%optforces/= 0) then
     if (optforces/= 0) then
       choice=2;cpopt=1; nnlout=3*natom
     end if
   end if
 end if

 ABI_ALLOCATE(enlout,(nnlout*blocksize))

!Allocation of memory space for one WF
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
 if (gs_hamk%usepaw==1.and.(iscf>0.or.gs_hamk%usecprj==1)) then
   iorder_cprj=0
   nband_k_cprj=nband_k*(mband_cprj/dtset%mband)
   bandpp_cprj=mpi_enreg%bandpp
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor*bandpp_cprj))
   ncpgr=0;if (cpopt==1) ncpgr=cprj(1,1)%ncpgr
   call pawcprj_alloc(cwaveprj,ncpgr,gs_hamk%dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 end if

#undef DEV_NEW_CODE
!#define DEV_NEW_CODE

!The code below is more efficient if paral_kgb==1 (less MPI communications)
!however OMP is not compatible with paral_kgb since we should define
!which threads performs the call to MPI_ALL_REDUCE.
!This problem can be easily solved by removing MPI_enreg from meanvalue_g so that
!the MPI call is done only once outside the OMP parallel region.

#ifdef DEV_NEW_CODE
!Loop over bands or blocks of bands. Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
!!$OMP PARALLEL DO PRIVATE(iband,ar)
 do iblock=1,nblockbd

!  Compute kinetic energy of each band
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize

     call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),0)

     ek_k(iband)=ar

     if(paw_dmft%use_dmft==1) then
       do iband1=1,nband_k
         call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
&         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),paw_dmft%use_dmft,ar_im=ar_im)
         ek_k_nd(1,iband,iband1)=ar
         ek_k_nd(2,iband,iband1)=ar_im
       end do
     end if
!    if(use_dmft==1) then
!    do iband1=1,nband_k
!    call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
!    &         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
!    &         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),use_dmft)
!    ek_k_nd(iband,iband1)=ar
!    end do
!    end if

   end do
 end do
!TODO: xmpi_sum is missing but I have to understand the logic used to deal with the different
!MPI options and communicators.
#endif


!Loop over bands or blocks of bands. Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
 do iblock=1,nblockbd
   occblock=maxval(occ_k(1+(iblock-1)*blocksize:iblock*blocksize))
   cwavef(:,:)=cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)

#ifndef DEV_NEW_CODE
!  Compute kinetic energy of each band
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize

     call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),&
&     cg(:,1+(iband-1)*npw_k*my_nspinor+icg:iband*npw_k*my_nspinor+icg),0)

     ek_k(iband)=ar

     if(paw_dmft%use_dmft==1) then
       do iband1=1,nband_k
         call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&         cg(:,1+(iband -1)*npw_k*my_nspinor+icg:iband *npw_k*my_nspinor+icg),&
&         cg(:,1+(iband1-1)*npw_k*my_nspinor+icg:iband1*npw_k*my_nspinor+icg),paw_dmft%use_dmft,ar_im=ar_im)
         ek_k_nd(1,iband,iband1)=ar
         ek_k_nd(2,iband,iband1)=ar_im
       end do
     end if
   end do
#endif

   if(iscf>0)then
     ! In case of fixed occupation numbers, accumulates the partial density
     if (fixed_occ .and. mpi_enreg%paral_kgb/=1) then
       if (abs(occ_k(iblock))>=tol8) then
         weight=occ_k(iblock)*wtk/gs_hamk%ucvol
!        Accumulate charge density in real space in array rhoaug

!        The same section of code is also found in mkrho.F90 : should be rationalized !
         call fourwf(1,rhoaug(:,:,:,1),cwavef,dummy,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&         istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&         gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&         tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

         if(dtset%nspinor==2)then
           ABI_ALLOCATE(cwavef1,(2,npw_k))
           cwavef1(:,:)=cwavef(:,1+npw_k:2*npw_k) ! EB FR spin dn part and used for m_z component (cwavef_z)

           if(dtset%nspden==1) then

             call fourwf(1,rhoaug(:,:,:,1),cwavef1,dummy,wfraug,&
&             gs_hamk%gbound_k,gs_hamk%gbound_k,&
&             istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

           else if(dtset%nspden==4) then
             ! Build the four components of rho. We use only norm quantities and, so fourwf.
             ! $\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
             ! $\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
             ! $\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
             ABI_ALLOCATE(cwavef_x,(2,npw_k))
             ABI_ALLOCATE(cwavef_y,(2,npw_k))
             !$(\Psi^{1}+\Psi^{2})$
             cwavef_x(:,:)=cwavef(:,1:npw_k)+cwavef1(:,1:npw_k)
             !$(\Psi^{1}-i \Psi^{2})$
             cwavef_y(1,:)=cwavef(1,1:npw_k)+cwavef1(2,1:npw_k)
             cwavef_y(2,:)=cwavef(2,1:npw_k)-cwavef1(1,1:npw_k)
             ! z component
             call fourwf(1,rhoaug(:,:,:,4),cwavef1,dummy,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&             istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)
             ! x component
             call fourwf(1,rhoaug(:,:,:,2),cwavef_x,dummy,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&             istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)
             ! y component
             call fourwf(1,rhoaug(:,:,:,3),cwavef_y,dummy,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&             istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,gs_hamk%mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,&
&             gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,1,&
&             tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

             ABI_DEALLOCATE(cwavef_x)
             ABI_DEALLOCATE(cwavef_y)

           end if ! dtset%nspden/=4
           ABI_DEALLOCATE(cwavef1)
         end if
       else
         nskip=nskip+1
       end if

!      In case of fixed occupation numbers,in bandFFT mode accumulates the partial density
     else if (fixed_occ .and. mpi_enreg%paral_kgb==1) then

       if (dtset%nspinor==1) then
         call timab(537,1,tsec) ! "prep_fourwf%vtow"
         call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavef,wfraug,iblock,istwf_k,&
&         gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,npw_k,&
&         gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,&
&         1,gs_hamk%ucvol,wtk,use_gpu_cuda=dtset%use_gpu_cuda)
         call timab(537,2,tsec)
       else if (dtset%nspinor==2) then
         ABI_ALLOCATE(cwavefb,(2,npw_k*blocksize,2))
         ibs=(iblock-1)*npw_k*my_nspinor*blocksize+icg
!        --- No parallelization over spinors ---
         if (mpi_enreg%paral_spinor==0) then
           do iband=1,blocksize
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,1)=cg(:,1+(2*iband-2)*npw_k+ibs:(iband*2-1)*npw_k+ibs)
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,2)=cg(:,1+(2*iband-1)*npw_k+ibs:iband*2*npw_k+ibs)
           end do
         else
!          --- Parallelization over spinors ---
!          (split the work between 2 procs)
           cwavefb(:,:,3-ispinor_index)=zero
           do iband=1,blocksize
             cwavefb(:,(iband-1)*npw_k+1:iband*npw_k,ispinor_index) = cg(:,1+(iband-1)*npw_k+ibs:iband*npw_k+ibs)
           end do
           call xmpi_sum(cwavefb,mpi_enreg%comm_spinor,ierr)
         end if

         call timab(537,1,tsec) !"prep_fourwf%vtow"
         if (nspinor1TreatedByThisProc) then
           call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavefb(:,:,1),wfraug,iblock,&
&           istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,npw_k,&
&           gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&           use_gpu_cuda=dtset%use_gpu_cuda)
         end if
         if(dtset%nspden==1) then
           if (nspinor2TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,1),blocksize,cwavefb(:,:,2),wfraug,&
&             iblock,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,&
&             gs_hamk%ngfft,npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,&
&             gs_hamk%ucvol,wtk,use_gpu_cuda=dtset%use_gpu_cuda)
           end if
         else if(dtset%nspden==4) then
           ABI_ALLOCATE(cwavef_x,(2,npw_k*blocksize))
           ABI_ALLOCATE(cwavef_y,(2,npw_k*blocksize))
           cwavef_x(:,:)=cwavefb(:,1:npw_k*blocksize,1)+cwavefb(:,:,2)
           cwavef_y(1,:)=cwavefb(1,1:npw_k*blocksize,1)+cwavefb(2,:,2)
           cwavef_y(2,:)=cwavefb(2,:,1)-cwavefb(1,:,2)
           if (nspinor1TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,4),blocksize,cwavefb(:,:,2),wfraug,&
&             iblock,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&             use_gpu_cuda=dtset%use_gpu_cuda)
           end if
           if (nspinor2TreatedByThisProc) then
             call prep_fourwf(rhoaug(:,:,:,2),blocksize,cwavef_x,wfraug,&
&             iblock,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&             use_gpu_cuda=dtset%use_gpu_cuda)
             call prep_fourwf(rhoaug(:,:,:,3),blocksize,cwavef_y,wfraug,&
&             iblock,istwf_k,gs_hamk%mgfft,mpi_enreg,nband_k,ndat,gs_hamk%ngfft,&
&             npw_k,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,occ_k,1,gs_hamk%ucvol,wtk,&
&             use_gpu_cuda=dtset%use_gpu_cuda)
           end if
           ABI_DEALLOCATE(cwavef_x) 
           ABI_DEALLOCATE(cwavef_y) 
         end if                     
         call timab(537,2,tsec)
         ABI_DEALLOCATE(cwavefb)
       end if
     end if
   end if ! End of SCF calculation

!    Call to nonlocal operator:
!    - Compute nonlocal forces from most recent wfs
!    - PAW: compute projections of WF onto NL projectors (cprj)
   if(iscf>0.or.gs_hamk%usecprj==1)then
     if (gs_hamk%usepaw==1.or.optforces/=0) then
!      Treat all wavefunctions in case of varying occupation numbers or PAW
!      Only treat occupied bands in case of fixed occupation numbers and NCPP
       if(fixed_occ.and.abs(occblock)<=tol8.and.gs_hamk%usepaw==0) then
         if (optforces>0) grnl_k(:,(iblock-1)*blocksize+1:iblock*blocksize)=zero
       else
         if(gs_hamk%usepaw==1) then
           call timab(554,1,tsec)  ! "vtowfk:rhoij"
         end if
         if(cpopt==1) then
           iband=1+(iblock-1)*bandpp_cprj
           call pawcprj_copy(cprj(:,1+(iblock-1)*my_nspinor*blocksize+ibg:iblock*my_nspinor*blocksize+ibg),cwaveprj)
         end if
         if (mpi_enreg%paral_kgb==1) then
           call timab(572,1,tsec) ! 'prep_nonlop%vtowfk'
           call prep_nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir, &
&           eig_k(1+(iblock-1)*blocksize:iblock*blocksize),blocksize,&
&           mpi_enreg,nnlout,paw_opt,signs,nonlop_dum,tim_nonlop_prep,cwavef,cwavef)
           call timab(572,2,tsec)
         else
           call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,eig_k(1+(iblock-1)*blocksize:iblock*blocksize),&
&           mpi_enreg,blocksize,nnlout,&
&           paw_opt,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
         end if
         if(gs_hamk%usepaw==1) then
           call timab(554,2,tsec)
         end if
!        Acccumulate forces
         if (optforces>0) then
           iband=(iblock-1)*blocksize
           do iblocksize=1,blocksize
             ii=0
             if (nnlout>3*natom) ii=6
             iband=iband+1;ibs=ii+nnlout*(iblocksize-1)
             grnl_k(1:nnlout,iband)=enlout(ibs+1:ibs+nnlout)
           end do
         end if
!        Store cprj (<Pnl|Psi>)
         if (gs_hamk%usepaw==1.and.gs_hamk%usecprj==1) then
           iband=1+(iblock-1)*bandpp_cprj
           call pawcprj_put(gs_hamk%atindx,cwaveprj,cprj,natom,iband,ibg,ikpt,iorder_cprj,isppol,&
&           mband_cprj,dtset%mkmem,natom,bandpp_cprj,nband_k_cprj,gs_hamk%dimcprj,my_nspinor,&
&           dtset%nsppol,dtfil%unpaw,mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if
       end if
     end if ! PAW or forces
   end if ! iscf>0 or iscf=-3

 end do !  End of loop on blocks

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(enlout)

 if (gs_hamk%usepaw==1.and.(iscf>0.or.gs_hamk%usecprj==1)) then
   call pawcprj_free(cwaveprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 if (fixed_occ.and.iscf>0) then
   ABI_DEALLOCATE(wfraug)
 end if

!Write the number of one-way 3D ffts skipped until now (in case of fixed occupation numbers
 if(iscf>0 .and. fixed_occ .and. (prtvol>2 .or. ikpt<=nkpt_max) )then
   write(message,'(a,i0)')' vtowfk: number of one-way 3D ffts skipped in vtowfk until now =',nskip
   call wrtout(std_out,message,'PERS')
 end if

!Norm-conserving or FockACE: Compute nonlocal+FockACE part of total energy : rotate subvnlx
 if (gs_hamk%usepaw==0 .and. wfopta10 /= 1 .and. .not. newlobpcg ) then
   call timab(586,1,tsec)   ! 'vtowfk(nonlocalpart)'
   ABI_ALLOCATE(matvnl,(2,nband_k,nband_k))
   ABI_ALLOCATE(mat1,(2,nband_k,nband_k))
   mat1=zero

   if (wfopta10==4) then
     enlx_k(1:nband_k)=zero

     if (istwf_k==1) then
       call zhemm('l','l',nband_k,nband_k,cone,totvnlx,nband_k,evec,nband_k,czero,mat1,nband_k)
       do iband=1,nband_k
         res = cg_real_zdotc(nband_k,evec(:,iband),mat1(:,:,iband))
         enlx_k(iband)= res
       end do
     else if (istwf_k==2) then
       ABI_ALLOCATE(evec_loc,(nband_k,nband_k))
       ABI_ALLOCATE(mat_loc,(nband_k,nband_k))
       do iband=1,nband_k
         do jj=1,nband_k
           evec_loc(iband,jj)=evec(2*iband-1,jj)
         end do
       end do
       call dsymm('l','l',nband_k,nband_k,one,totvnlx,nband_k,evec_loc,nband_k,zero,mat_loc,nband_k)
       do iband=1,nband_k
         enlx_k(iband)=ddot(nband_k,evec_loc(:,iband),1,mat_loc(:,iband),1)
       end do
       ABI_DEALLOCATE(evec_loc)
       ABI_DEALLOCATE(mat_loc)
     end if

   else
!    MG: This version is much faster with good OMP scalability.
!    Construct upper triangle of matvnl from subvnlx using full storage mode.
     pidx=0
     do jj=1,nband_k
       do ii=1,jj
         pidx=pidx+1
         matvnl(1,ii,jj)=subvnlx(2*pidx-1)
         matvnl(2,ii,jj)=subvnlx(2*pidx  )
       end do
     end do

     call zhemm('L','U',nband_k,nband_k,cone,matvnl,nband_k,evec,nband_k,czero,mat1,nband_k)

!$OMP PARALLEL DO PRIVATE(res)
     do iband=1,nband_k
       res = cg_real_zdotc(nband_k,evec(:,iband),mat1(:,:,iband))
       enlx_k(iband) = res
     end do
   end if

   ABI_DEALLOCATE(matvnl)
   ABI_DEALLOCATE(mat1)
   call timab(586,2,tsec)
 end if

!###################################################################

 if (iscf<=0 .and. residk > dtset%tolwfr) then
   write(message,'(a,2(i0,1x),a,es13.5)')&
&   'Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,' max resid= ',residk
   MSG_WARNING(message)
 end if

!Print out eigenvalues (hartree)
 if (prtvol>2 .or. ikpt<=nkpt_max) then
   write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&   'eigenvalues (hartree) for',nband_k,'bands',ch10,&
&   '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
   call wrtout(std_out,message,'PERS')
   do ii=0,(nband_k-1)/6
     write(message, '(1p,6e12.4)' ) (eig_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
     call wrtout(std_out,message,'PERS')
   end do
 else if(ikpt==nkpt_max+1)then
   call wrtout(std_out,' vtowfk : prtvol=0 or 1, do not print more k-points.','PERS')
 end if

!Print out decomposition of eigenvalues in the non-selfconsistent case or if prtvol>=10
 if( (iscf<0 .and. (prtvol>2 .or. ikpt<=nkpt_max)) .or. prtvol>=10)then
   write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&   ' mean kinetic energy (hartree) for ',nband_k,' bands',ch10,&
&   '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
   call wrtout(std_out,message,'PERS')

   do ii=0,(nband_k-1)/6
     write(message, '(1p,6e12.4)' ) (ek_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
     call wrtout(std_out,message,'PERS')
   end do

   if (gs_hamk%usepaw==0) then
     write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&     ' mean NL+Fock-type energy (hartree) for ',nband_k,' bands',ch10,&
&     '              after ',inonsc,' non-SCF iterations with ',dtset%nline,' CG line minimizations'
     call wrtout(std_out,message,'PERS')

     do ii=0,(nband_k-1)/6
       write(message,'(1p,6e12.4)') (enlx_k(iband),iband=1+6*ii,min(6+6*ii,nband_k))
       call wrtout(std_out,message,'PERS')
     end do
   end if
 end if

 !Hamiltonian constructor for gwls_sternheimer
 if(dtset%optdriver==RUNL_GWLS) then
   call build_H(dtset,mpi_enreg,cpopt,cg,gs_hamk,kg_k,kinpw)
 end if

 if(wfopta10 /= 1 .and. .not. newlobpcg) then
   ABI_DEALLOCATE(evec)
   ABI_DEALLOCATE(subham)
   ABI_DEALLOCATE(totvnlx)
   ABI_DEALLOCATE(subvnlx)
   ABI_DEALLOCATE(subovl)
 end if
! if ( .not. newlobpcg .and. .not. enable_cgwf_paw) then
 ABI_DEALLOCATE(gsc)
! end if

 if(wfoptalg==3) then
   ABI_DEALLOCATE(eig_save)
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(message,'(a,a,a,i0,a)')' vtowfk : exit ',ch10,'  prtvol=-',level,', debugging mode => stop '
   MSG_ERROR(message)
 end if

 call timab(30,2,tsec)
 call timab(28,2,tsec)

 DBG_EXIT("COLL")

end subroutine vtowfk
!!***

!!****f* ABINIT/fxphas
!!
!! NAME
!! fxphas
!!
!! FUNCTION
!! Fix phase of all bands. Keep normalization but maximize real part
!! (minimize imag part). Also fix the sign of real part
!! by setting the first non-zero element to be positive.
!! See also fxphas_seq in m_cgtools
!!
!! INPUTS
!!  cg(2,mcg)= contains the wavefunction |c> coefficients.
!!  gsc(2,mgsc)= if useoverlap==1, contains the S|c> coefficients,
!!               where S is an overlap matrix.
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwfk=input option parameter that describes the storage of wfs
!!    (set to 1 if usual complex vectors)
!!  mcg=size of second dimension of cg
!!  mgsc=size of second dimension of gsc
!!  mpi_enreg=information about MPI parallelization
!!  nband_k=number of bands
!!  npw_k=number of planewaves
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identi0,ty_matrix)
!!               1: wavefunctions are overlapping
!!
!! OUTPUT
!!  cg(2,mcg)=same array with altered phase.
!!  gsc(2,mgsc)= same array with altered phase.
!!
!! NOTES
!! When the sign of the real part was fixed (modif v3.1.3g.6), the
!! test Tv3#5 , dataset 5, behaved differently than previously.
!! This should be cleared up.
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine fxphas(cg,gsc,icg,igsc,istwfk,mcg,mgsc,mpi_enreg,nband_k,npw_k,useoverlap)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc*useoverlap)

!Local variables-------------------------------
!scalars
 integer :: iband,ierr,ii,indx
 real(dp) :: cim,cre,gscim,gscre,quotient,root1,root2,saa,sab,sbb,theta
 real(dp) :: thppi,xx,yy
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(nband_k,2),buffer3(nband_k,3),tsec(2)
 real(dp),allocatable :: cimb(:),creb(:),saab(:),sabb(:),sbbb(:) !,sarr(:,:)

! *************************************************************************

!The general case, where a complex phase indeterminacy is present
 if(istwfk==1)then

   ABI_ALLOCATE(cimb,(nband_k))
   ABI_ALLOCATE(creb,(nband_k))
   ABI_ALLOCATE(saab,(nband_k))
   ABI_ALLOCATE(sabb,(nband_k))
   ABI_ALLOCATE(sbbb,(nband_k))
   cimb(:)=zero ; creb(:)=zero

!  Loop over bands
!  TODO: MG store saa arrays in sarr(3,nband_k) to reduce false sharing.
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nband_k,icg,npw_k,cg,saab,sbbb,sabb)
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

!    Compute several sums over Re, Im parts of c
     saa=zero; sbb=zero; sab=zero
     do ii=1+indx,npw_k+indx
       saa=saa+cg(1,ii)*cg(1,ii)
       sbb=sbb+cg(2,ii)*cg(2,ii)
       sab=sab+cg(1,ii)*cg(2,ii)
     end do
     saab(iband)=saa
     sbbb(iband)=sbb
     sabb(iband)=sab
   end do

!  XG030513 : MPIWF : should transmit saab,sbbb,sabb from the procs
!  of the WF group to the master processor of the WF group
   if (mpi_enreg%paral_kgb == 1) then
     buffer3(:,1)=saab(:)
     buffer3(:,2)=sbbb(:)
     buffer3(:,3)=sabb(:)
     call timab(48,1,tsec)
     call xmpi_sum(buffer3,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(buffer3,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
     saab(:)=buffer3(:,1)
     sbbb(:)=buffer3(:,2)
     sabb(:)=buffer3(:,3)
   end if

!  XG030513 : MPIWF this loop should only be executed by the master of the WF group

   if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then
     do iband=1,nband_k
       indx=icg+(iband-1)*npw_k

       saa=saab(iband)
       sbb=sbbb(iband)
       sab=sabb(iband)

!      Get phase angle theta
       if (sbb+saa>tol8)then
         if(abs(sbb-saa)>tol8*(sbb+saa) .or. 2*abs(sab)>tol8*(sbb+saa))then
           if (abs(sbb-saa)>tol8*abs(sab)) then
             quotient=sab/(sbb-saa)
             theta=0.5_dp*atan(2.0_dp*quotient)
           else
!            Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
             theta=0.25_dp*(pi-(sbb-saa)/sab)
           end if
!          Check roots to get theta for max Re part
           root1=cos(theta)**2*saa+sin(theta)**2*sbb-2.0_dp*cos(theta)*sin(theta)*sab
           thppi=theta+0.5_dp*pi
           root2=cos(thppi)**2*saa+sin(thppi)**2*sbb-2.0_dp*cos(thppi)*sin(thppi)*sab
           if (root2>root1) theta=thppi
         else
!          The real part vector and the imaginary part vector are orthogonal, and of same norm. Strong indeterminacy.
!          Will determine the first non-zero coefficient, and fix its phase
!          Hypothesis : there is at least one non-zero element on the master node ...
           do ii=1+indx,npw_k+indx
             cre=cg(1,ii)
             cim=cg(2,ii)
             if(cre**2+cim**2>tol8**2*(saa+sbb))then
               if(cre**2>tol8**2**cim**2)then
                 theta=atan(cim/cre)
               else
!                Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
                 theta=pi/2-cre/cim
               end if
               exit
             end if
           end do
         end if
       else
         write(message,'(a,i0,5a)')&
&         'The eigenvector with band ',iband,' has zero norm.',ch10,&
&         'This usually happens when the number of bands (nband) is comparable to the number of planewaves (mpw)',ch10,&
&         'Action: Check the parameters of the calculation. If nband ~ mpw, then decrease nband or, alternatively, increase ecut'
         MSG_ERROR(message)
       end if

       xx=cos(theta)
       yy=sin(theta)

!      Here, set the first non-zero element to be positive
!      Comment the next nine lines to recover the behaviour of pre v3.1.3g
!      Hypothesis : there is at least one non-zero element on the master node ...
       do ii=1+indx,npw_k+indx
         cre=cg(1,ii)
         cim=cg(2,ii)
         cre=xx*cre-yy*cim
         if(abs(cre)>tol8)exit
       end do
       if(cre<zero)then
         xx=-xx ; yy=-yy
       end if

       creb(iband)=xx
       cimb(iband)=yy

     end do
   end if

!  XG030513 : MPIWF : should transmit creb(:),cimb(:) of the master
!  processor of the WF group to the others procs of the WF group
   if (mpi_enreg%paral_kgb == 1) then
     call timab(48,1,tsec)
     buffer2(:,1)=creb(:)
     buffer2(:,2)=cimb(:)
     call xmpi_sum(buffer2,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(buffer2,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
     creb(:)=buffer2(:,1)
     cimb(:)=buffer2(:,2)
   end if

!  MG TODO: Scaling can be done with zscal
!$OMP PARALLEL DO PRIVATE(indx,xx,yy,cre,cim,gscre,gscim)
   do iband=1,nband_k
     indx=icg+(iband-1)*npw_k

     xx=creb(iband)
     yy=cimb(iband)
!    Alter phase of array |cg>
     do ii=1+indx,npw_k+indx
       cre=cg(1,ii)
       cim=cg(2,ii)
       cg(1,ii)=xx*cre-yy*cim
       cg(2,ii)=xx*cim+yy*cre
     end do

!    Alter phase of array S|cg>
     if (useoverlap==1) then
       indx=igsc+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         gscre=gsc(1,ii)
         gscim=gsc(2,ii)
         gsc(1,ii)=xx*gscre-yy*gscim
         gsc(2,ii)=xx*gscim+yy*gscre
       end do
     end if
   end do ! iband

   ABI_DEALLOCATE(cimb)
   ABI_DEALLOCATE(creb)
   ABI_DEALLOCATE(saab)
   ABI_DEALLOCATE(sabb)
   ABI_DEALLOCATE(sbbb)

 else  ! if istwfk/=1.  Storages that take into account the time-reversal symmetry : the freedom is only a sign freedom

   ABI_ALLOCATE(creb,(nband_k))
   creb(:)=zero
!  XG030513 : MPIWF : this loop should be done only by the master processor of the WF group

   if (mpi_enreg%paral_kgb==0.or.mpi_enreg%me_fft==0) then

!    Loop over bands
     do iband=1,nband_k

       indx=icg+(iband-1)*npw_k

!      Here, set the first non-zero real element to be positive
       do ii=1+indx,npw_k+indx
         cre=cg(1,ii)
         if(abs(cre)>tol8)exit
       end do
       creb(iband)=cre

     end do ! iband

   end if
!  XG030513 : MPIWF : should transmit cre(:) of the master processor of the WF group to the others
   if (mpi_enreg%paral_kgb == 1) then
     call timab(48,1,tsec)
     call xmpi_sum(creb,mpi_enreg%comm_fft,ierr)
     if (mpi_enreg%paral_spinor==1) then
       call xmpi_sum(creb,mpi_enreg%comm_spinor,ierr)
     end if
     call timab(48,2,tsec)
   end if

   do iband=1,nband_k
     cre=creb(iband)
     if(cre<zero)then
       indx=icg+(iband-1)*npw_k
       do ii=1+indx,npw_k+indx
         cg(1,ii)=-cg(1,ii)
         cg(2,ii)=-cg(2,ii)
       end do
       if(useoverlap==1)then
         indx=igsc+(iband-1)*npw_k
         do ii=1+indx,npw_k+indx
           gsc(1,ii)=-gsc(1,ii)
           gsc(2,ii)=-gsc(2,ii)
         end do
       end if
     end if
   end do ! iband

   ABI_DEALLOCATE(creb)
 end if ! istwfk

end subroutine fxphas
!!***

end module m_vtowfk
!!***
