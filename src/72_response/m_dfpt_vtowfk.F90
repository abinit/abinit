!!****m* ABINIT/m_dfpt_vtowfk
!! NAME
!!  m_dfpt_vtowfk
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG, AR, DRH, MB, MVer,XW, MT, GKA)
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

module m_dfpt_vtowfk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_mpinfo
 use m_cgtools
 use m_wfk
 use m_rf2
 use m_dtset
 use m_dtfil


 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,  only : MPI_type
 use m_rf2_init,     only : rf2_init
 use m_time,         only : timab
 use m_pawrhoij,     only : pawrhoij_type
 use m_pawcprj
! use m_pawcprj,      only : pawcprj_type, pawcprj_alloc, pawcprj_put, pawcprj_free, pawcprj_get, pawcprj_copy, pawcprj_zaxpby, pawcprj_set_zero
 use m_hamiltonian,  only : gs_hamiltonian_type, rf_hamiltonian_type, KPRIME_H_KPRIME
 use m_spacepar,     only : meanvalue_g
 use m_dfpt_mkrho,   only : dfpt_accrho
 use m_dfpt_cgwf,    only : dfpt_cgwf
 use m_getghc,       only : getgsc

 implicit none

 private
!!***

 public :: dfpt_vtowfk
!!***

contains
!!***

!!****f* ABINIT/dfpt_vtowfk
!! NAME
!! dfpt_vtowfk
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed potential (vlocal1).
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex=1 if rhoaug1 is real, 2 if rhoaug1 is complex
!TODO MJV: PAW mband_mem
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  dim_eig2rf = dimension for the second order eigenvalues
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig0_k(nband_k)=GS eigenvalues at k (hartree)
!!  eig0_kq(nband_k)=GS eigenvalues at k+Q (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  grad_berry(2,mpw1,dtefield%mband_occ) = the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  ibg1=shift to be applied on the location of data in the array cprj1
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=k-point index number
!!  ipert=type of the perturbation
!!  isppol=1 index of current spin component
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mkmem =number of k points trated by this node (GS data).
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  n4,n5,n6 used for dimensioning real space arrays
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q
!!  rf_hamk_dir2 <type(rf_hamiltonian_type)>= (used only when ipert=natom+11, so q=0)
!!    same as rf_hamkq, but the direction of the perturbation is different
!!  rhoaug1(cplex*n4,n5,n6,nspden)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this ratio has been attributed to the band n (second argument), zero otherwise
!!  ddk<wfk_t>=struct info for DDK file.
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the occupied states.
!!  cg1_active(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the active. Only needed for ieigrf/=0
!!  edocc_k(nband_k)=correction to 2nd-order total energy coming
!!      from changes of occupation
!!  eeig0_k(nband_k)=zero-order eigenvalues contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  ek0_k(nband_k)=0-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  ek1_k(nband_k)=1st-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eloc0_k(nband_k)=zero-order local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl0_k(nband_k)=zero-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl1_k(nband_k)=first-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  gh1c_set(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(1)}|nK>
!!  gh0c1_set(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(0)}k+q-eig^{(0)}nk|\Psi^{(1)}kq>
!!      The wavefunction is orthogonal to the active space (for metals). It is not coherent with cg1.
!!  resid_k(nband_k)=residuals for each band over all k points,
!!  rhoaug1(cplex*n4,n5,n6,nspden)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!  ==== if (gs_hamkq%usepaw==1) ====
!TODO MJV: PAW mband_mem
!!    cprj1(natom,nspinor*mband*mk1mem*nsppol*usecprj)=
!!              1st-order wave functions at k,q projected with non-local projectors:
!!                       cprj1=<p_i|C1nk,q> where p_i is a non-local projector
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                            (cumulative, so input as well as output)
!!
!! PARENTS
!!      dfpt_vtorho
!!
!! CHILDREN
!!      cg_zcopy,corrmetalwf1,dfpt_accrho,dfpt_cgwf,dotprod_g,getgsc
!!      matrixelmt_g,meanvalue_g,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawcprj_get,pawcprj_put,rf2_destroy,rf2_init,sqnorm_g,status,timab
!!      wfk_read_bks,wrtout
!!
!! SOURCE

subroutine dfpt_vtowfk(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,&
& dim_eig2rf,dtfil,dtset,&
& edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,&
& ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&
& fermie1,ffnl1,ffnl1_test,gh0c1_set,gh1c_set,grad_berry,gs_hamkq,&
& ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,ipert,&
& isppol,mband,mband_mem,mcgq,mcprjq,mkmem,mk1mem,&
& mpi_enreg,mpw,mpw1,natom,nband_k,ncpgr,&
& nnsclo_now,npw_k,npw1_k,nspinor,nsppol,&
& n4,n5,n6,occ_k,pawrhoij1,prtvol,psps,resid_k,rf_hamkq,rf_hamk_dir2,rhoaug1,rocceig,&
& ddk_f,wtk_k,nlines_done,cg1_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dim_eig2rf,ibg
 integer,intent(in) :: ibg1,ibgq,icg,icg1,icgq,idir,ikpt,ipert,isppol
 integer,intent(in) :: mband,mcgq,mcprjq,mk1mem,mkmem
 integer,intent(in) :: mband_mem
 integer,intent(in) :: mpw,mpw1,n4,n5,n6,natom,ncpgr
 integer,intent(in) :: nnsclo_now,nspinor,nsppol,prtvol
 integer,optional,intent(in) :: cg1_out
 integer,intent(in) :: nband_k,npw1_k,npw_k
 integer,intent(inout) :: nlines_done
 real(dp),intent(in) :: fermie1,wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq,rf_hamk_dir2
 type(pseudopotential_type),intent(in) :: psps
!arrays
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: eig0_k(nband_k),eig0_kq(nband_k)
 real(dp),intent(in) :: ffnl1(:,:,:,:),ffnl1_test(:,:,:,:)
 real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband_k)
 real(dp),intent(in) :: occ_k(nband_k),rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)
 real(dp),intent(inout) :: rhoaug1(cplex*n4,n5,n6,gs_hamkq%nvloc)
 real(dp),intent(inout) :: cg1_active(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(inout) :: gh1c_set(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(inout) :: gh0c1_set(2,mpw1*nspinor*mband_mem*mk1mem*nsppol*dim_eig2rf)
 real(dp),intent(inout) :: edocc_k(nband_k),eeig0_k(nband_k),eig1_k(2*nband_k**2)
 real(dp),intent(out) :: ek0_k(nband_k),eloc0_k(nband_k)
 real(dp),intent(inout) :: ek1_k(nband_k)
 real(dp),intent(out) :: enl0_k(nband_k),enl1_k(nband_k)
 real(dp),intent(out) :: resid_k(nband_k)
!TODO: PAW distrib bands mband_mem
 type(pawcprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*gs_hamkq%usecprj)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq)
 type(pawcprj_type),intent(inout) :: cprj1(natom,nspinor*mband*mk1mem*nsppol*gs_hamkq%usecprj)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(natom*gs_hamkq%usepaw)
 type(wfk_t),intent(inout) :: ddk_f(4)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=14,tim_fourwf=5
 integer,save :: nskip=0
 integer :: iband,idir0,ierr,igs,igscq,ii,dim_dcwf,inonsc
 integer :: iband_me,nband_me, unit_me
 integer :: iorder_cprj,iorder_cprj1,ipw,iscf_mod,ispinor,me,mgscq,nkpt_max
 integer :: option,opt_gvnlx1,quit,test_ddk
 integer :: tocceig,usedcwavef,ptr,shift_band
 real(dp) :: aa,ai,ar,eig0nk,resid,residk,scprod,energy_factor
 character(len=500) :: message
 type(rf2_t) :: rf2
!arrays
 logical,allocatable :: cycle_bands(:)
 integer :: band_procs(nband_k)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwave0(:,:),cwave1(:,:),cwavef(:,:)
 real(dp),allocatable :: dcwavef(:,:),gh1c_n(:,:),gh0c1(:,:)
 real(dp),allocatable :: gsc(:,:),gscq(:,:),gvnlx1(:,:),gvnlxc(:,:)
 real(dp),pointer :: kinpw1(:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:),cwaveprj0(:,:),cwaveprj1(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Keep track of total time spent in dfpt_vtowfk
 call timab(128,1,tsec)

 nkpt_max=50; if (xmpi_paral==1) nkpt_max=-1

 if(prtvol>2 .or. ikpt<=nkpt_max)then
   write(message,'(2a,i5,2x,a,3f9.5,2x,a)')ch10,' Non-SCF iterations; k pt #',ikpt,'k=',&
&   gs_hamkq%kpt_k(:),'band residuals:'
   call wrtout(std_out,message,'PERS')
 end if

!Initializations and allocations
 me=mpi_enreg%me_kpt
 quit=0

!The value of iscf must be modified if ddk perturbation
 iscf_mod=dtset%iscf;if(ipert==natom+1.or.ipert==natom+10.or.ipert==natom+11) iscf_mod=-3

 kinpw1 => gs_hamkq%kinpw_kp
 ABI_ALLOCATE(gh0c1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
 ABI_ALLOCATE(cwave1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gh1c_n,(2,npw1_k*nspinor))
 if (gs_hamkq%usepaw==1) then
   ABI_ALLOCATE(gsc,(2,npw1_k*nspinor))
 else
   ABI_ALLOCATE(gsc,(0,0))
 end if

!Read the npw and kg records of wf files
 test_ddk=0
 if ((ipert==natom+2.and.sum((dtset%qptn(1:3))**2)<1.0d-7.and.&
& (dtset%berryopt/= 4.and.dtset%berryopt/= 6.and.dtset%berryopt/= 7.and.&
& dtset%berryopt/=14.and.dtset%berryopt/=16.and.dtset%berryopt/=17)).or.&
& ipert==natom+10.or.ipert==natom+11) then
   test_ddk=1
   if(ipert==natom+10.or.ipert==natom+11) test_ddk=0
 end if

!Additional stuff for PAW
 ABI_DATATYPE_ALLOCATE(cwaveprj0,(0,0))
 if (gs_hamkq%usepaw==1) then
!  1-Compute all <g|S|Cnk+q>
   igscq=0
!TODO MJV: PAW mband_mem
   mgscq=mpw1*nspinor*mband
   ABI_MALLOC_OR_DIE(gscq,(2,mgscq), ierr)

   call getgsc(cgq,cprjq,gs_hamkq,gscq,ibgq,icgq,igscq,ikpt,isppol,mcgq,mcprjq,&
&   mgscq,mpi_enreg,natom,nband_k,npw1_k,dtset%nspinor,select_k=KPRIME_H_KPRIME)
!  2-Initialize additional scalars/arrays
   iorder_cprj=0;iorder_cprj1=0
   dim_dcwf=npw1_k*nspinor;if (ipert==natom+2.or.ipert==natom+10.or.ipert==natom+11) dim_dcwf=0
   ABI_ALLOCATE(dcwavef,(2,dim_dcwf))
   if (gs_hamkq%usecprj==1) then
     ABI_DATATYPE_DEALLOCATE(cwaveprj0)
     ABI_DATATYPE_ALLOCATE(cwaveprj0,(natom,nspinor))
     call pawcprj_alloc(cwaveprj0,1,gs_hamkq%dimcprj)
   end if
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,nspinor))
   ABI_DATATYPE_ALLOCATE(cwaveprj1,(natom,nspinor))
   call pawcprj_alloc(cwaveprj ,0,gs_hamkq%dimcprj)
   call pawcprj_alloc(cwaveprj1,0,gs_hamkq%dimcprj)
 else
   igscq=0;mgscq=0;dim_dcwf=0
   ABI_ALLOCATE(gscq,(0,0))
   ABI_ALLOCATE(dcwavef,(0,0))
   ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
   ABI_DATATYPE_ALLOCATE(cwaveprj1,(0,0))
 end if

 energy_factor=two
 if(ipert==natom+10.or.ipert==natom+11) energy_factor=six

!For rf2 perturbation :
 if(ipert==natom+10.or.ipert==natom+11) then
   call rf2_init(cg,cprj,rf2,dtset,dtfil,eig0_k,eig1_k,ffnl1,ffnl1_test,gs_hamkq,ibg,icg,idir,ikpt,ipert,isppol,mkmem,&
   mpi_enreg,mpw,nband_k,nsppol,rf_hamkq,rf_hamk_dir2,occ_k,rocceig,ddk_f)
 end if

 call timab(139,1,tsec)

!======================================================================
!==================  LOOP OVER BANDS ==================================
!======================================================================

 call proc_distrb_band(band_procs,mpi_enreg%proc_distrb,ikpt,isppol,mband,&
&  mpi_enreg%me_band,mpi_enreg%me_kpt,mpi_enreg%comm_band)
#ifdef DEV_MJV
print *, 'band_procs  ', band_procs
#endif

 iband_me = 0
 do iband=1,nband_k

!  Skip bands not treated by current proc
   if( (mpi_enreg%proc_distrb(ikpt, iband,isppol)/=me)) cycle
   iband_me = iband_me + 1
 
!unit_me = 300+iband
unit_me = 6
!  Get ground-state wavefunctions
   ptr = 1+(iband_me-1)*npw_k*nspinor+icg
   call cg_zcopy(npw_k*nspinor,cg(1,ptr),cwave0)

!  Get PAW ground state projected WF (cprj)
   if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1.and.ipert/=natom+10.and.ipert/=natom+11) then
     idir0 = idir
     if(ipert==natom+3.or.ipert==natom+4) idir0 =1
!TODO MJV: PAW distribute cprj mband_mem
     call pawcprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,natom,iband,ibg,ikpt,iorder_cprj,&
&     isppol,mband,mkmem,natom,1,nband_k,nspinor,nsppol,dtfil%unpaw,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,&
&     icpgr=idir0,ncpgr=ncpgr)
   end if

!  Get first-order wavefunctions
   ptr = 1+(iband_me-1)*npw1_k*nspinor+icg1
   call cg_zcopy(npw1_k*nspinor,cg1(1,ptr),cwavef)

!  Read PAW projected 1st-order WF (cprj)
!  Unuseful for the time being (will be recomputed in dfpt_cgwf)
!  if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!TODO MJV: PAW
!  call pawcprj_get(gs_hamkq%atindx1,cwaveprj,cprj1,natom,iband,ibg1,ikpt,iorder_cprj1,&
!  &    isppol,mband,mk1mem,natom,1,nband_k,nspinor,nsppol,dtfil%unpaw1,
!  &    mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
!  end if

!  Filter the wavefunctions for large modified kinetic energy
!  The GS wavefunctions should already be non-zero
   do ispinor=1,nspinor
     igs=(ispinor-1)*npw1_k
#ifdef DEV_MJV
print *, ' igs ', igs, '  huge 0, cutoff = ', huge(zero), huge(zero)*1.d-11
#endif
     do ipw=1+igs,npw1_k+igs
#ifdef DEV_MJV
print *, 'ipw, igs kinpw1 ', ipw, igs, kinpw1(ipw-igs)
#endif
       if(kinpw1(ipw-igs)>huge(zero)*1.d-20)then
#ifdef DEV_MJV
print *, 'zeroing ipw ', ipw
#endif
         cwavef(1,ipw)=zero
         cwavef(2,ipw)=zero
       end if
     end do
   end do
#ifdef DEV_MJV
print *, ' cwavef elements filtered ', cwavef(:,23), cwavef(:,49), cwavef(:,50)
#endif


!  If electric field, the derivative of the wf should be read, and multiplied by i.
   if(test_ddk==1) then
     ii = ddk_f(1)%findk(gs_hamkq%kpt_k)
     ABI_CHECK(ii == ikpt, "ii != ikpt, something is wrong with k-point, check kptopt/ngkpt, etc")
!TODO MJV: check if this iband should be _me
     call ddk_f(1)%read_bks(iband, ikpt, isppol, xmpio_single, cg_bks=gvnlx1)

!    Multiplication by -i
!    MVeithen 021212 : use + i instead,
!    See X. Gonze, Phys. Rev. B 55, 10337 (1997) [[cite:Gonze1997]] Eq. (79)
!    the operator used to compute the first-order derivative
!    of the wavefunctions with respect to an electric field
!    is $+i \frac{d}{dk}$
!    This change will affect the computation of the 2dtes from non
!    stationary expressions, see dfpt_nstdy.f and dfpt_nstwf.f
     do ipw=1,npw1_k*nspinor
!      aa=gvnlx1(1,ipw)
!      gvnlx1(1,ipw)=gvnlx1(2,ipw)
!      gvnlx1(2,ipw)=-aa
       aa=gvnlx1(1,ipw)
       gvnlx1(1,ipw)=-gvnlx1(2,ipw)
       gvnlx1(2,ipw)=aa
     end do
   end if

!  Unlike in GS calculations, the inonsc loop is inside the band loop
!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!  (often 1 for SCF calculation, =nstep for non-SCF calculations)
   do inonsc=1,nnsclo_now

!    Note that the following translation occurs in the called routine :
!    iband->band, nband_k->nband, npw_k->npw, npw1_k->npw1
     eig0nk=eig0_k(iband)
     usedcwavef=gs_hamkq%usepaw;if (dim_dcwf==0) usedcwavef=0
     if (inonsc==1) usedcwavef=2*usedcwavef
     opt_gvnlx1=0;if (ipert==natom+2) opt_gvnlx1=1
     if (ipert==natom+2.and.gs_hamkq%usepaw==1.and.inonsc==1) opt_gvnlx1=2

     if ( (ipert/=natom+10 .and. ipert/=natom+11) .or. abs(occ_k(iband))>tol8 ) then
       nband_me = proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)

#ifdef DEV_MJV
print *, ' vtowfk isppol.ikpt, nband_me ', isppol, ikpt, nband_me, iband, iband_me 
print *, ' occ_k, eig0nk,eig0_kq ', occ_k(iband), eig0nk,eig0_kq
print *, 'cgwf cwavef 444 ', cwavef(:,23)
print *, 'cgwf cwave0 444 ', cwave0(:,1:5)
print *, 'cgwf cgq 444 ', cgq(:,1:5)
#endif
       call dfpt_cgwf(iband,iband_me,band_procs,dtset%berryopt,cgq,cwavef,cwave0,cwaveprj,cwaveprj0,&
&       rf2,dcwavef,&
&       eig0nk,eig0_kq,eig1_k,gh0c1,gh1c_n,grad_berry,gsc,gscq,gs_hamkq,gvnlxc,gvnlx1,icgq,&
&       idir,ipert,igscq,mcgq,mgscq,mpi_enreg,mpw1,natom,nband_k,nband_me,dtset%nbdbuf,dtset%nline,&
&       npw_k,npw1_k,nspinor,opt_gvnlx1,prtvol,quit,resid,rf_hamkq,dtset%dfpt_sciss,dtset%tolrde,&
&       dtset%tolwfr,usedcwavef,dtset%wfoptalg,nlines_done)
       resid_k(iband)=resid
     else
       resid_k(iband)=zero
     end if
#ifdef DEV_MJV
print *, 'cgwf cwavef 454 ', cwavef(:,23)
#endif

     if (ipert/=natom+10 .and. ipert/= natom+11) then
!    At this stage, the 1st order function cwavef is orthogonal to cgq (unlike
!    when it is input to dfpt_cgwf). Here, restore the "active space" content
!    of the first-order wavefunction, to give cwave1.
!    PAW: note that dcwavef (1st-order change of WF due to overlap change)
!         remains in the subspace orthogonal to cgq
       call proc_distrb_cycle_bands(cycle_bands, mpi_enreg%proc_distrb,ikpt,isppol,me)
#ifdef DEV_MJV
print *, 'prtfull1wf ', dtset%prtfull1wf
print *, 'vtowfk after corrmetal iband, cwavef ', iband, cwavef(:,1:5)
print *, 'cgwf cwavef 465 ', cwavef(:,23)
print *, 'vtowfk after corrmetal iband, cwave1 ', iband, cwave1(:,1:5)
#endif
       if (dtset%prtfull1wf>0) then
         call full_active_wf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,cycle_bands,eig1_k,fermie1,&
&         eig0nk,eig0_kq,dtset%elph2_imagden,iband,ibgq,icgq,mcgq,mcprjq,mpi_enreg,natom,nband_k,npw1_k,nspinor,&
&         0,gs_hamkq%usepaw)
         edocc_k=zero
         tocceig=1
       else
         call corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,cycle_bands,edocc_k,eig1_k,fermie1,gh0c1,&
&         iband,ibgq,icgq,gs_hamkq%istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband_k,npw1_k,nspinor,&
&         occ_k,rocceig,0,gs_hamkq%usepaw,tocceig)
       end if
#ifdef DEV_MJV
print *,'vtowfk after corrmetal iband, cwave1 ', iband, cwave1(:,1:5)
ii = (iband-1)*nband_k*2
print *, 'vtowfk iband, eig1_k ', iband, eig1_k(ii+1:ii+min(10,nband_k))
#endif
       ABI_DEALLOCATE (cycle_bands)
     else
       tocceig=0
       call cg_zcopy(npw1_k*nspinor,cwavef,cwave1)
       if (gs_hamkq%usepaw==1) then
         call pawcprj_copy(cwaveprj,cwaveprj1)
       end if
     end if

     if (abs(occ_k(iband))<= tol8) then
       ek0_k(iband)=zero
       ek1_k(iband)=zero
       eeig0_k(iband)=zero
       enl0_k(iband)=zero
       enl1_k(iband)=zero
       eloc0_k(iband)=zero
       nskip=nskip+1
     else
!      Compute the 0-order kinetic operator contribution (with cwavef)
       call meanvalue_g(ar,kinpw1,0,gs_hamkq%istwf_k,mpi_enreg,npw1_k,nspinor,cwavef,cwavef,0)
#ifdef DEV_MJV
print *, ' ik isppol iband ar kinpw1 ', ikpt, isppol, iband, ar, kinpw1(1:30) 
print *, ' cwavef elements ', cwavef(:,23), cwavef(:,49), cwavef(:,50)
print *, ' rf_hamkq%dkinpw_k w1 w0 ', rf_hamkq%dkinpw_k(23), cwave1(:,23),cwave0(:,23)
#endif
!      There is an additional factor of 2 with respect to the bare matrix element
       ek0_k(iband)=energy_factor*ar
!      Compute the 1-order kinetic operator contribution (with cwave1 and cwave0), if needed.
!      Note that this is called only for ddk or strain, so that npw1_k=npw_k
       if(ipert==natom+1 .or. ipert==natom+3 .or. ipert==natom+4)then
         call matrixelmt_g(ai,ar,rf_hamkq%dkinpw_k,gs_hamkq%istwf_k,0,npw_k,nspinor,cwave1,cwave0,&
&         mpi_enreg%me_g0, mpi_enreg%comm_fft)
!        There is an additional factor of 4 with respect to the bare matrix element
         ek1_k(iband)=two*energy_factor*ar
       end if

!      Compute eigenvalue part of total energy (with cwavef)
       if (gs_hamkq%usepaw==1) then
         call dotprod_g(scprod,ai,gs_hamkq%istwf_k,npw1_k*nspinor,1,cwavef,gsc,mpi_enreg%me_g0,&
&         mpi_enreg%comm_spinorfft)
       else
         call sqnorm_g(scprod,gs_hamkq%istwf_k,npw1_k*nspinor,cwavef,mpi_enreg%me_g0,&
&         mpi_enreg%comm_fft)
       end if
       eeig0_k(iband)=-energy_factor*(eig0_k(iband)- (dtset%dfpt_sciss) )*scprod

!      Compute nonlocal psp contributions to nonlocal energy:
!      <G|Vnl+VFockACE|C1nk(perp)> is contained in gvnlxc (with cwavef)
       call dotprod_g(scprod,ai,gs_hamkq%istwf_k,npw1_k*nspinor,1,cwavef,gvnlxc,mpi_enreg%me_g0,&
&       mpi_enreg%comm_spinorfft)
       enl0_k(iband)=energy_factor*scprod

       if(ipert/=natom+10.and.ipert/=natom+11) then
!        <G|Vnl1|Cnk> is contained in gvnlx1 (with cwave1)
         call dotprod_g(scprod,ai,gs_hamkq%istwf_k,npw1_k*nspinor,1,cwave1,gvnlx1,mpi_enreg%me_g0,&
&         mpi_enreg%comm_spinorfft)
         enl1_k(iband)=two*energy_factor*scprod
       end if

!      Removal of the 1st-order kinetic energy from the 1st-order non-local part.
       if(ipert==natom+1 .or. ipert==natom+3 .or. ipert==natom+4) then
         enl1_k(iband)=enl1_k(iband)-ek1_k(iband)
       end if

!      Accumulate 1st-order density (only at the last inonsc)
!      Accumulate zero-order potential part of the 2nd-order total energy
!   BUGFIX from Max Stengel: need to initialize eloc at each inonsc iteration, in case nnonsc > 1
       eloc0_k(iband) = zero
       option=2;if (iscf_mod>0.and.inonsc==nnsclo_now) option=3
       call dfpt_accrho(cplex,cwave0,cwave1,cwavef,cwaveprj0,cwaveprj1,eloc0_k(iband),&
&       gs_hamkq,iband,idir,ipert,isppol,dtset%kptopt,mpi_enreg,natom,nband_k,ncpgr,&
&       npw_k,npw1_k,nspinor,occ_k,option,pawrhoij1,rhoaug1,tim_fourwf,tocceig,wtk_k)
       if(ipert==natom+10.or.ipert==natom+11) eloc0_k(iband)=energy_factor*eloc0_k(iband)/two

       if(ipert==natom+10.or.ipert==natom+11) then
         shift_band=(iband-1)*npw1_k*nspinor
         call dotprod_g(scprod,ai,gs_hamkq%istwf_k,npw1_k*nspinor,1,cwave1,&
&         rf2%RHS_Stern(:,1+shift_band:npw1_k*nspinor+shift_band),mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
         ek1_k(iband)=two*energy_factor*scprod
       end if

     end if ! End of non-zero occupation

!    Exit loop over inonsc if converged and if non-self-consistent
     if (iscf_mod<0 .and. resid<dtset%tolwfr) exit

   end do ! End loop over inonsc

!  Get first-order eigenvalues and wavefunctions
   ptr = 1+(iband_me-1)*npw1_k*nspinor+icg1
   if (.not. present(cg1_out)) then
     call cg_zcopy(npw1_k*nspinor,cwave1,cg1(1,ptr))
   end if
   if(dim_eig2rf > 0) then
     if (.not. present(cg1_out)) then
       cg1_active(:,1+(iband_me-1)*npw1_k*nspinor+icg1:iband_me*npw1_k*nspinor+icg1)=cwavef(:,:)
     end if
     gh1c_set(:,1+(iband_me-1)*npw1_k*nspinor+icg1:iband_me*npw1_k*nspinor+icg1)=gh1c_n(:,:)
     gh0c1_set(:,1+(iband_me-1)*npw1_k*nspinor+icg1:iband_me*npw1_k*nspinor+icg1)=gh0c1(:,:)
   end if

!  PAW: write first-order projected wavefunctions
   if (psps%usepaw==1.and.gs_hamkq%usecprj==1) then
!TODO MJV: PAW distribute cprj and cprj1 over bands? mband_mem
     call pawcprj_put(gs_hamkq%atindx,cwaveprj,cprj1,natom,iband,ibg1,ikpt,iorder_cprj1,isppol,&
&     mband,mk1mem,natom,1,nband_k,gs_hamkq%dimcprj,nspinor,nsppol,dtfil%unpaw1,&
&     mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,to_be_gathered=.true.)
   end if

write (unit_me, *) 'vtowfk iband cg1 : ', iband, cwave1(:,1:5)
 end do

!======================================================================
!==================  END LOOP OVER BANDS ==============================
!======================================================================

!For rf2 perturbation
 if(ipert==natom+10.or.ipert==natom+11) call rf2_destroy(rf2)

!Find largest resid over bands at this k point
 residk=maxval(resid_k(:))
 if (prtvol>2 .or. ikpt<=nkpt_max) then
   do ii=0,(nband_k-1)/8
     write(message,'(1p,8e10.2)')(resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
     call wrtout(std_out,message,'PERS')
   end do
 end if

 call timab(139,2,tsec)
 call timab(130,1,tsec)

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(cwave1)
 ABI_DEALLOCATE(gh0c1)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(gvnlx1)
 ABI_DEALLOCATE(gh1c_n)

 if (gs_hamkq%usepaw==1) then
   call pawcprj_free(cwaveprj)
   call pawcprj_free(cwaveprj1)
   if (gs_hamkq%usecprj==1) then
     call pawcprj_free(cwaveprj0)
   end if
 end if
 ABI_DEALLOCATE(dcwavef)
 ABI_DEALLOCATE(gscq)
 ABI_DEALLOCATE(gsc)
 ABI_DATATYPE_DEALLOCATE(cwaveprj0)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj1)


!###################################################################

 ! Write the number of one-way 3D ffts skipped until now (in case of fixed occupation numbers)
 if (iscf_mod>0 .and. (prtvol>2 .or. ikpt<=nkpt_max)) then
   write(message,'(a,i0)')' dfpt_vtowfk : number of one-way 3D ffts skipped in vtowfk3 until now =',nskip
   call wrtout(std_out,message,'PERS')
 end if

 if (prtvol<=2 .and. ikpt==nkpt_max+1) then
   write(message,'(3a)') ch10,' dfpt_vtowfk : prtvol=0, 1 or 2, do not print more k-points.',ch10
   call wrtout(std_out,message,'PERS')
 end if

 if (residk>dtset%tolwfr .and. iscf_mod<=0 .and. iscf_mod/=-3) then
   write(message,'(a,2i0,a,es13.5)')'Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,' max resid=',residk
   MSG_WARNING(message)
 end if

 call timab(130,2,tsec)
 call timab(128,2,tsec)

 DBG_EXIT('COLL')

end subroutine dfpt_vtowfk
!!***

!!****f* ABINIT/full_active_wf1
!!
!! NAME
!! full_active_wf1
!!
!! FUNCTION
!! Response function calculation only:
!! Restore the full "active space" contribution to the 1st-order wavefunctions.
!! The 1st-order WF corrected in this way will no longer be othogonal to the other occupied states.
!! This routine will be only used in a non self-consistent calculation of the
!! 1st-order WF for post-processing purposes. Therefore, it does not compute
!! the contribution of the 2DTE coming from the change of occupations.
!!
!! INPUTS
!!  cgq(2,mcgq)=planewave coefficients of wavefunctions at k+q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors
!!  cwavef(2,npw1*nspinor)= 1st-order wave-function before correction
!!  cwaveprj(natom,nspinor)= 1st-order wave-function before correction
!!                           projected on NL projectors (PAW)
!!  cycle_bands(nband)=array of logicals for bands we have on this cpu
!!  eig1(2*nband**2)=first-order eigenvalues (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  eig0nk=energy of the band at k being corrected
!!  eig0_kq(nband)=energies of the bands at k+q
!!  elph2_imagden=imaginary parameter to broaden the energy denominators
!!  iband=index of current band
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  npw1=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  timcount=index used to accumulate timing (0 from dfpt_vtowfk, 1 from dfpt_nstwf)
!!  usepaw=flag for PAW
!!
!! OUTPUT
!!  cwave1(2,npw1*nspinor)= 1st-order wave-function after correction
!!  cwaveprj1(natom,nspinor)= 1st-order wave-function after correction
!!                            projected on NL projectors (PAW)
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_zcopy,dotprod_g,pawcprj_copy,pawcprj_zaxpby,timab
!!
!! SOURCE

subroutine full_active_wf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,cycle_bands,eig1,&
&               fermie1,eig0nk,eig0_kq,elph2_imagden,&
&               iband,ibgq,icgq,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,&
&               nspinor,timcount,usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ibgq,icgq,mcgq,mcprjq,natom,nband,npw1,nspinor,timcount,usepaw
 real(dp),intent(in) :: fermie1, eig0nk
 real(dp),intent(in) :: elph2_imagden
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 logical,intent(in)  :: cycle_bands(nband)
 real(dp),intent(in) :: cgq(2,mcgq),cwavef(2,npw1*nspinor)
 real(dp),intent(in) :: eig0_kq(nband)
 real(dp),intent(in) :: eig1(2*nband**2)
 real(dp),intent(out) :: cwave1(2,npw1*nspinor)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq),cwaveprj(natom,nspinor*usepaw)
 type(pawcprj_type),intent(inout) :: cwaveprj1(natom,nspinor*usepaw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ibandkq,index_cgq,index_cprjq,index_eig1,ii
 integer :: ibandkq_me, ierr
 real(dp) :: facti,factr,eta,delta_E,inv_delta_E,gkkr
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

 DBG_ENTER("COLL")
 ABI_UNUSED(mpi_enreg%comm_cell)

 call timab(214+timcount,1,tsec)

!At this stage, the 1st order function cwavef is orthogonal to cgq (unlike when it is input to dfpt_cgwf).
!Here, restore the "active space" content of the 1st-order wavefunction, to give cwave1 .

! New logic 11/11/2019: accumulate correction in cwave1 and cwaveprj1, 
!   then add it to cwavef at the end with a modified blas call
 cwave1 = zero
 
 if (usepaw==1) then
   call pawcprj_set_zero(cwaveprj1)
 end if

 eta = elph2_imagden

!Loop over WF at k+q subspace
 ibandkq_me = 0
 do ibandkq=1,nband

!TODO MJV: here we have an issue - the cgq are no longer present for all bands!
!   we only have diagonal terms for iband iband1 and ibandq in same set of bands
! 1) filter with distrb
   if(cycle_bands(ibandkq)) cycle
   ibandkq_me = ibandkq_me + 1

! 2) get contributions for correction factors of cgq from bands present on this cpu

   delta_E = eig0nk - eig0_kq(ibandkq)
   inv_delta_E = delta_E / ( delta_E ** 2 + eta ** 2)

   index_eig1=2*ibandkq-1+(iband-1)*2*nband
   index_cgq=npw1*nspinor*(ibandkq_me-1)+icgq

   if(ibandkq==iband) then
     gkkr = eig1(index_eig1) - fermie1
   else
     gkkr = eig1(index_eig1)
   end if
   factr = inv_delta_E * gkkr
   facti = inv_delta_E * eig1(index_eig1+1)

!  Apply correction to 1st-order WF
!$OMP PARALLEL DO PRIVATE(ii) SHARED(cgq,cwave1,facti,factr,index_cgq,npw1,nspinor)
   do ii=1,npw1*nspinor
     cwave1(1,ii)=cwave1(1,ii)+(factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq))
     cwave1(2,ii)=cwave1(2,ii)+(facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq))
   end do

!  In the PAW case, also apply correction to projected WF
   if (usepaw==1) then
     index_cprjq=nspinor*(ibandkq-1)+ibgq
!TODO : PAW distrib cprj over bands 
     !index_cprjq=nspinor*(ibandkq_me-1)+ibgq
     call pawcprj_zaxpby((/factr,facti/),(/one,zero/),cprjq(:,index_cprjq+1:index_cprjq+nspinor),cwaveprj1)
   end if

 end do ! Loop over k+q subspace

! 3) reduce over bands to get all contributions to correction
! need MPI reduce over band communicator only
 call xmpi_sum(cwave1,mpi_enreg%comm_band,ierr)
 if (usepaw==1) then
   call pawcprj_mpi_sum(cwaveprj1, mpi_enreg%comm_band, ierr)
 end if

! 4) add correction to the cwave1
!Now add on input WF into output WF
 call cg_zaxpy(npw1*nspinor,(/one,zero/),cwavef,cwave1)

!Idem for cprj
 if (usepaw==1) then
   call pawcprj_zaxpby((/one,zero/),(/one,zero/),cwaveprj,cwaveprj1)
 end if

 call timab(214+timcount,2,tsec)

 DBG_EXIT("COLL")

end subroutine full_active_wf1
!!***

!!****f* ABINIT/corrmetalwf1
!!
!! NAME
!! corrmetalwf1
!!
!! FUNCTION
!! Response function calculation only:
!! Correct 1st-order wave-function, taking into account "metallic" occupations.
!! 1st-order WF orthogonal to C_n,k+q, restore the "active space" content of the first-order WF.
!! receives a single band at k as input, and works on all bands at k+q
!!
!! INPUTS
!!  cgq(2,mcgq)=planewave coefficients of wavefunctions at k+q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors
!!  cwavef(2,npw1*nspinor)= 1st-order wave-function before correction
!!  cwaveprj(natom,nspinor)= 1st-order wave-function before correction
!!                           projected on NL projectors (PAW)
!!  cycle_bands(nband)=array of logicals for bands we have on this cpu
!!  eig1(2*nband**2)=first-order eigenvalues (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  ghc(2,npw1*nspinor)=<G|H0-eig0_k.I|C1 band,k> (NCPP) or <G|H0-eig0_k.S0|C1 band,k> (PAW)
!!                      (C1 before correction)
!!  iband=index of current band
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  npw1=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ(nband)=occupation number for each band for each k.
!!  rocceig(nband,nband)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this ratio has been attributed to the band n (second argument), zero otherwise
!!  timcount=index used to accumulate timing (0 from dfpt_vtowfk, 1 from dfpt_nstwf)
!!  usepaw=flag for PAW
!!
!! OUTPUT
!!  cwave1(2,npw1*nspinor)= 1st-order wave-function after correction
!!  cwaveprj1(natom,nspinor)= 1st-order wave-function after correction
!!                            projected on NL projectors (PAW)
!!  edocc(nband)=correction to 2nd-order total energy coming from changes of occupations
!!  wf_corrected=flag put to 1 if input cwave1 is effectively different from output cwavef
!!
!! NOTES
!!  Was part of dfpt_vtowfk before.
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_zcopy,dotprod_g,pawcprj_copy,pawcprj_zaxpby,timab
!!
!! SOURCE

subroutine corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,cycle_bands,edocc,eig1,fermie1,ghc,iband, &
&          ibgq,icgq,istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,nspinor,occ,rocceig,timcount,&
&          usepaw,wf_corrected)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ibgq,icgq,istwf_k,mcgq,mcprjq,natom,nband,npw1,nspinor,timcount,usepaw
 integer,intent(out) :: wf_corrected
 real(dp),intent(in) :: fermie1
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 logical,intent(in) :: cycle_bands(nband)
 real(dp),intent(in) :: cgq(2,mcgq),cwavef(2,npw1*nspinor)
 real(dp),intent(in) :: eig1(2*nband**2),ghc(2,npw1*nspinor),occ(nband),rocceig(nband,nband)
 real(dp),intent(out) :: cwave1(2,npw1*nspinor),edocc(nband)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq),cwaveprj(natom,nspinor*usepaw)
 type(pawcprj_type),intent(inout) :: cwaveprj1(natom,nspinor*usepaw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ibandkq,index_cgq,index_cprjq,index_eig1,ii
 integer :: ibandkq_me, ierr, iband_
 real(dp) :: facti,factr,invocc
 real(dp) :: edocc_tmp
!arrays
 integer :: bands_treated_now(nband)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwcorr(:,:)
 type(pawcprj_type) :: cwaveprj1_corr(natom,nspinor*usepaw)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(214+timcount,1,tsec)

 bands_treated_now = 0
 bands_treated_now(iband) = 1
 call xmpi_sum(bands_treated_now,mpi_enreg%comm_band,ierr)

!At this stage, the 1st order function cwavef is orthogonal to cgq (unlike when it is input to dfpt_cgwf).
!Here, restore the "active space" content of the 1st-order wavefunction, to give cwave1 .

 ABI_ALLOCATE(cwcorr,(2,npw1*nspinor))

 wf_corrected=0

 cwave1 = zero

 if (usepaw==1) then
   call pawcprj_copy(cwaveprj,cwaveprj1)
   call pawcprj_alloc(cwaveprj1_corr, cwaveprj1(1,1)%ncpgr, cwaveprj1(:,1)%nlmn)
 end if

 edocc(iband)=zero

! loop iband_ over all bands being treated for the moment
! all procs in pool of bands should be working on the same iband_ at a given time
! I will save in _my_ array cwave1, if iband==iband_
 do iband_ = 1, nband
   if (bands_treated_now(iband_) == 0) cycle
#ifdef DEV_MJV
print *, 'iband_ ', iband_, ' nband ', nband
#endif

   cwcorr = zero
   if (usepaw==1) then
     call pawcprj_set_zero(cwaveprj1_corr)
   end if
   edocc_tmp = zero

!Correct WF only for occupied states
   if (abs(occ(iband_)) > tol8) then
     invocc=one/occ(iband_)


!    Loop over WF at k+q subspace
     ibandkq_me = 0
     do ibandkq=1,nband
       if(cycle_bands(ibandkq)) cycle
       ibandkq_me = ibandkq_me + 1
#ifdef DEV_MJV
print *, 'ibandkq ', ibandkq, ibandkq_me, rocceig(ibandkq,iband_)
#endif

!      Select bands with variable occupation
       if (abs(rocceig(ibandkq,iband_))>tol8) then

         if (iband_ == iband) wf_corrected=1

         index_eig1=2*ibandkq-1+(iband_-1)*2*nband
         index_cgq=npw1*nspinor*(ibandkq_me-1)+icgq

         if(ibandkq==iband_) then
           factr=rocceig(ibandkq,iband_)*invocc*(eig1(index_eig1)-fermie1)
         else
           factr=rocceig(ibandkq,iband_)*invocc*eig1(index_eig1)
         end if
         facti = rocceig(ibandkq,iband_)*invocc*eig1(index_eig1+1)
#ifdef DEV_MJV
print *, 'factr, facti ', factr, facti
#endif

!        Apply correction to 1st-order WF
!$OMP PARALLEL DO PRIVATE(ii) SHARED(cgq,cwcorr,facti,factr,index_cgq,npw1,nspinor)
         do ii=1,npw1*nspinor
           cwcorr(1,ii)=cwcorr(1,ii)+(factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq))
           cwcorr(2,ii)=cwcorr(2,ii)+(facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq))
         end do

!        In the PAW case, also apply correction to projected WF
         if (usepaw==1) then
           index_cprjq=nspinor*(ibandkq-1)+ibgq
!TODO : distribution of PAW cprj over bands 
           !index_cprjq=nspinor*(ibandkq_me-1)+ibgq
           call pawcprj_zaxpby((/factr,facti/),(/one,zero/),cprjq(:,index_cprjq+1:index_cprjq+nspinor),cwaveprj1_corr)
         end if

!        The factor of two is needed because we compute the 2DTE, and not E(2)
         edocc_tmp = edocc_tmp-two*(factr*eig1(index_eig1)+facti*eig1(index_eig1+1))

       end if ! Variable occupations
     end do ! Loop over k+q subspace
   end if ! if occupied states

! 3) reduce over bands to get all contributions to correction
! need MPI reduce over band communicator only
   call xmpi_sum(cwcorr, mpi_enreg%comm_band, ierr)
   if (usepaw==1) then
     call pawcprj_mpi_sum(cwaveprj1_corr, mpi_enreg%comm_band, ierr)
   end if

! this sums over the k+q contributions to the present iband_
   call xmpi_sum(edocc_tmp, mpi_enreg%comm_band, ierr)


! 4) add correction to the cwave1
! if I have iband_, correct my cwave1
   if (iband_==iband) then
     edocc(iband) = edocc_tmp
     cwave1 = cwcorr
     call cg_zaxpy(npw1*nspinor,(/one,zero/),cwavef,cwave1)
!Idem for cprj
     if (usepaw==1) then
       call pawcprj_zaxpby((/one,zero/),(/one,zero/),cwaveprj1_corr,cwaveprj1)
     end if
   end if
 end do ! loop over all bands presently running in parallel

#ifdef DEV_MJV
print *, 'edocc ', edocc
#endif

!In the PAW case, compute <Psi^(1)_ortho|H-Eig0_k.S|Psi^(1)_parallel> contribution to 2DTE
 if (usepaw==1.and.wf_corrected==1) then
!$OMP WORKSHARE
   cwcorr(:,:)=cwave1(:,:)-cwavef(:,:)
!$OMP END WORKSHARE
   call dotprod_g(factr,facti,istwf_k,npw1*nspinor,1,cwcorr,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   edocc(iband)=edocc(iband)+four*factr
 end if

 ABI_DEALLOCATE(cwcorr)
 if (usepaw==1) then
   call pawcprj_free(cwaveprj1_corr)
 end if

 call timab(214+timcount,2,tsec)

 DBG_EXIT("COLL")

end subroutine corrmetalwf1
!!***

end module m_dfpt_vtowfk
!!***
