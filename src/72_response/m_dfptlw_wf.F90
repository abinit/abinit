!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_wf
!! NAME
!!  m_dfptlw_wf
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_dfptlw_wf
    
 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_dtset
 use m_errors
 use m_profiling_abi
 use m_hamiltonian
 use m_cgtools
 use m_pawcprj
 use m_pawfgr
 use m_wfk
 use m_xmpi
 use m_getgh1c
 use m_mklocl

 use m_fstrings, only : itoa, sjoin
 use m_io_tools, only : file_exists
 use m_time, only : cwtime
 use m_kg, only : mkkpg

 implicit none

 public :: dfpt_1wf

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_wf/dfpt_1wf
!! NAME
!!  dfpt_1wf
!!
!! FUNCTION
!!  Compute the spin, band and kpt resolved contributions
!!  to the spatial-dispersion third-order energy derivatives
!!  that depend on first-order response functions.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cplex: if 1, several magnitudes are REAL, if 2, COMPLEX
!!  ddk_f = wf files
!!  d2_dkdk_f = wf files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  gsqcut=large sphere cut-off
!!  icg=shift to be applied on the location of data in the array cg
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  ikpt=number of the k-point
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kpt(3)=reduced coordinates of k point
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  natom= number of atoms in the unit cell
!!  natpert=number of atomic displacement perturbations
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nylmgr=second dimension of ylmgr_k
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)=1-dimensional phases
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  ucvol=unit cell volume in bohr**3.
!!  useylmgr= if 1 use the derivative of spherical harmonics
!!  vpsp1_i1pertdq(cplex*nfft,nspden,n1dq)= local potential of first-order
!!          gradient Hamiltonian for i1pert
!!  vpsp1_i2pertdq(cplex*nfft,nspden,n2dq)= local potential of first-order
!!          gradient Hamiltonian for i2pert
!!  vtrial1_i1pert(cplex*nfft,nspden)=firs-order local potential
!!  vtrial1_i2pert(cplex*nfft,nspden)=firs-order local potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)=real spherical harmonics for the k point
!!  ylmgr_k(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)= k-gradients of real spherical
!!                                                                      harmonics for the k point
!!
!! OUTPUT
!! d3etot_t(1-5)_k= stationary 1wf contributions to the third-order energy
!!                  derivatives for kpt
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

subroutine dfpt_1wf(atindx,cg,cplex,ddk_f,d2_dkdk_f,&
     & d3etot_t1_k,d3etot_t2_k,d3etot_t3_k,&
     & d3etot_t4_k,d3etot_t5_k,dtset,gs_hamkq,gsqcut,icg,&
     & i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,ikpt,isppol,istwf_k,&
     & kg_k,kpt,kxc,mkmem,mpi_enreg,mpw,natom,nattyp,nband_k,&
     & n1dq,n2dq,nfft,ngfft,nkxc,npw_k,nspden,nsppol,nylmgr,occ_k,&
     & ph1d,psps,rhog,rhor,rmet,ucvol,useylmgr,&
     & vpsp1_i1pertdq,vpsp1_i2pertdq,vtrial1_i1pert,vtrial1_i2pert,&
     & xred,ylm_k,ylmgr_k)
    
 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer,intent(in) :: icg,ikpt,isppol,istwf_k
 integer,intent(in) :: mkmem,mpw,natom,nband_k,n1dq,n2dq,nfft
 integer,intent(in) :: nkxc,npw_k,nspden,nsppol,nylmgr
 integer,intent(in) :: useylmgr
 real(dp),intent(in) :: gsqcut,ucvol
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wfk_t),intent(inout) :: ddk_f,d2_dkdk_f

!arrays
 integer,intent(in) :: atindx(natom)
 integer,intent(in) :: kg_k(3,npw_k),nattyp(dtset%ntypat),ngfft(18)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(out) :: d3etot_t1_k(2)
 real(dp),intent(out) :: d3etot_t2_k(2)
 real(dp),intent(out) :: d3etot_t3_k(2)
 real(dp),intent(out) :: d3etot_t4_k(2,n2dq)
 real(dp),intent(out) :: d3etot_t5_k(2,n1dq)
 real(dp),intent(in) :: kpt(3),occ_k(nband_k),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,nspden),rmet(3,3)
 real(dp),intent(in) :: vpsp1_i1pertdq(2*nfft,nspden,n1dq)
 real(dp),intent(in) :: vpsp1_i2pertdq(2*nfft,nspden,n2dq)
 real(dp),intent(in) :: vtrial1_i1pert(cplex*nfft,nspden)
 real(dp),intent(in) :: vtrial1_i2pert(cplex*nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr_k(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)

!Local variables-------------------------------
!scalars
 integer :: berryopt,iband,jband,nkpg,nkpg1
 integer :: opt_gvnl1,optlocal,optnl,sij_opt,tim_getgh1c,usevnl,useylmgr1
 real(dp) :: dum_lambda

!arrays
 real(dp),allocatable :: cwave0i(:,:),cwave0j(:,:)
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: dkinpw(:),gv1c(:,:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg_k(:,:),kpg1_k(:,:)
 real(dp),allocatable :: part_ylmgr_k(:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: dum_vlocal(:,:,:,:),vlocal1(:,:,:,:),dum_vpsp(:)
 real(dp),allocatable :: vpsp1(:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)
 type(rf_hamiltonian_type),allocatable :: rf_hamkq(:)

 
! *************************************************************************

 DBG_ENTER("COLL")
 
!Additional definitions
 tim_getgh1c=0

!Additional allocations
 ABI_MALLOC(cwave0i,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(cwave0j,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(cwavef1,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(cwavef2,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(gv1c,(2,npw_k*dtset%nspinor))
 ABI_MALLOC(vlocal1,(cplex*ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
 ABI_MALLOC(dum_vpsp,(nfft))
 ABI_MALLOC(dum_vlocal,(ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
 ABI_MALLOC(vpsp1,(cplex*nfft))
 ABI_MALLOC(dum_cwaveprj,(0,0))
 ABI_MALLOC(part_ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
 part_ylmgr_k(:,:,:)=ylmgr_k(:,1:3,:)

!------------------------------------T1------------------------------------------------
!q1-gradient of gs Hamiltonian:
! < u_{i,k}^{\lambda1}} | \partial_{gamma} H^{(0)} | u_{i,k}^{\lambda2} >
!--------------------------------------------------------------------------------------

!Specific definitions
 vlocal1=zero
 useylmgr1=1
 dum_lambda=zero
 berryopt=0;optlocal=0;optnl=1;usevnl=0;opt_gvnl1=0;sij_opt=0

!Initialize rf Hamiltonian (the k-dependent part is prepared in getgh1c_setup)
 ABI_MALLOC(rf_hamkq,(1))
 call init_rf_hamiltonian(cplex,gs_hamkq,natom+1,rf_hamkq(1),& 
 & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
 & mpi_spintab=mpi_enreg%my_isppoltab)

 call rf_hamkq(1)%load_spin(isppol,vlocal1=vlocal1,with_nonlocal=.true.)

 !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
 call getgh1c_setup(gs_hamkq,rf_hamkq(1),dtset,psps,&                     ! In
 kpt,kpt,i3dir,natom+1,natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,& ! In
 npw_k,npw_k,useylmgr1,kg_k,ylm_k,kg_k,ylm_k,part_ylmgr_k,&               ! In
 dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)            ! Out


!Clean rf_hamiltonian
 call rf_hamkq(1)%free()
 ABI_FREE(rf_hamkq)

 !Deallocations
 ABI_FREE(kpg_k)
 ABI_FREE(kpg1_k)
 ABI_FREE(dkinpw)
 ABI_FREE(kinpw1)
 ABI_FREE(ffnlk)
 ABI_FREE(ffnl1)
 ABI_FREE(ph3d)


!Deallocations
 ABI_FREE(cwave0i)
 ABI_FREE(cwave0j)
 ABI_FREE(cwavef1)
 ABI_FREE(cwavef2)
 ABI_FREE(gv1c)
 ABI_FREE(vlocal1)
 ABI_FREE(dum_vpsp)
 ABI_FREE(dum_vlocal)
 ABI_FREE(vpsp1)
 ABI_FREE(dum_cwaveprj)
 ABI_FREE(part_ylmgr_k)


 DBG_EXIT("COLL")

end subroutine dfpt_1wf
!!***

end module m_dfptlw_wf
!!***
