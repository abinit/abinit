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
!!  dimffnl= third dimension of ffnl_k
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig1_k(2*nband_k**2)=1st-order eigenvalues at k for i1pert,i1dir
!!  eig2_k(2*nband_k**2)=1st-order eigenvalues at k for i2pert,i2dir
!!  ffnl_k(dtset%mpw,dimffnl,psps%lmnmax,psps%ntypat)= Nonlocal projectors and their derivatives for this k point
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  cg1 = first derivative of cg with respect the perturbation i1pert
!!  cg2 = first derivative of cg with respect the perturbation i2pert
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
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)=1-dimensional phases
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  samepert= .true. if i1pert=i2pert and i1dir=i2dir
!!  ucvol=unit cell volume in bohr**3.
!!  useylmgr= if 1 use the derivative of spherical harmonics
!!  vpsp1_i1pertdq(cplex*nfft,nspden,n1dq)= local potential of first-order
!!          gradient Hamiltonian for i1pert
!!  vpsp1_i2pertdq(cplex*nfft,nspden,n2dq)= local potential of first-order
!!          gradient Hamiltonian for i2pert
!!  wtk_k=weight assigned to the k point.
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

subroutine dfpt_1wf(atindx,cg,cg1,cg2,cplex,ddk_f,d2_dkdk_f,&
     & d3etot_t1_k,d3etot_t2_k,d3etot_t3_k,&
     & d3etot_t4_k,d3etot_t5_k,dimffnl,dtset,eig1_k,eig2_k,ffnl_k,gs_hamkq,gsqcut,icg,&
     & i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,ikpt,isppol,istwf_k,&
     & kg_k,kpt,kxc,mkmem,mpi_enreg,mpw,natom,nattyp,nband_k,&
     & n1dq,n2dq,nfft,ngfft,nkxc,npw_k,nspden,nsppol,nylmgr,occ_k,&
     & pawfgr,ph1d,psps,rhog,rhor,rmet,rprimd,samepert,ucvol,useylmgr,&
     & vpsp1_i1pertdq,vpsp1_i2pertdq,&
     & wtk_k,xred,ylm_k,ylmgr_k)
    
 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dimffnl,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer,intent(in) :: icg,ikpt,isppol,istwf_k
 integer,intent(in) :: mkmem,mpw,natom,nband_k,n1dq,n2dq,nfft
 integer,intent(in) :: nkxc,npw_k,nspden,nsppol,nylmgr
 integer,intent(in) :: useylmgr
 real(dp),intent(in) :: gsqcut,ucvol,wtk_k
 logical,intent(in) :: samepert
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wfk_t),intent(inout) :: ddk_f,d2_dkdk_f
 type(pawfgr_type),intent(in) :: pawfgr

!arrays
 integer,intent(in) :: atindx(natom)
 integer,intent(in) :: kg_k(3,npw_k),nattyp(dtset%ntypat),ngfft(18)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(out) :: d3etot_t1_k(2)
 real(dp),intent(out) :: d3etot_t2_k(2)
 real(dp),intent(out) :: d3etot_t3_k(2)
 real(dp),intent(out) :: d3etot_t4_k(2,n2dq)
 real(dp),intent(out) :: d3etot_t5_k(2,n1dq)
 real(dp),intent(in) :: eig1_k(2*nband_k**2),eig2_k(2*nband_k**2)
 real(dp),intent(in) :: ffnl_k(npw_k,dimffnl,psps%lmnmax,psps%ntypat)
 real(dp),intent(in) :: kpt(3),occ_k(nband_k),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: vpsp1_i1pertdq(2*nfft,nspden,n1dq)
 real(dp),intent(in) :: vpsp1_i2pertdq(2*nfft,nspden,n2dq)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr_k(npw_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)

!Local variables-------------------------------
!scalars
 integer :: berryopt,dimffnlk,dimffnl1,iband,idir,idq,ii,jband,nkpg,nkpg1,nylmgrtmp
 integer :: offset_cgi,offset_cgj,opt_gvnl1,optlocal,optnl,reuse_ffnlk,reuse_ffnl1,sij_opt
 integer :: size_wf,tim_getgh1c,usepaw,usevnl,useylmgr1
 real(dp) :: cprodi,cprodr,doti,dotr,dum_lambda,fac,tmpim,tmpre
 real(dp) :: cpu,wall,gflops
 logical :: with_nonlocal_i1pert,with_nonlocal_i2pert
 type(rf_hamiltonian_type) :: rf_hamkq

!arrays
 real(dp),allocatable :: cg1_aux(:,:),cg1_ddk(:,:,:),cwave0i(:,:),cwave0j(:,:)
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: dkinpw(:),gv1c(:,:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable :: gvloc1dqc(:,:),gvnl1dqc(:,:)
 real(dp) :: cj_h1_ci(2),dum_grad_berry(1,1),dum_gs1(1,1),dum_gvnl1(1,1)
 real(dp),allocatable :: kinpw1(:),kpg_k(:,:),kpg1_k(:,:)
 real(dp),allocatable :: part_ylmgr_k(:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: dum_vlocal(:,:,:,:),vlocal1(:,:,:,:),dum_vpsp(:)
 real(dp),allocatable :: vpsp1(:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)

 !AZ_test_ini**************************************************************************
 character(40) :: i1dir_text, i2dir_text, i3dir_text, iband_text, ikpt_text
 character(20) :: kpt_1_text, kpt_2_text, kpt_3_text, jband_text
 character(100) :: file_name
 real(dp) :: AZ_sum_re, AZ_sum_im
 !AZ_test_fin**************************************************************************
 
! *************************************************************************

 DBG_ENTER("COLL")
 
!Additional definitions
 tim_getgh1c=0
 useylmgr1=useylmgr
 usepaw=dtset%usepaw
 size_wf= dtset%nspinor*npw_k
 with_nonlocal_i1pert=.true. ; if (i1pert==natom+2) with_nonlocal_i1pert=.false.
 with_nonlocal_i2pert=.true. ; if (i2pert==natom+2) with_nonlocal_i2pert=.false.
 reuse_ffnlk=1 ; if (dtset%ffnl_lw==1) reuse_ffnlk=0
 reuse_ffnl1=1 ; if (dtset%ffnl_lw==1) reuse_ffnl1=0

!Additional allocations
 ABI_MALLOC(cwave0i,(2,size_wf))
 ABI_MALLOC(cwave0j,(2,size_wf))
 ABI_MALLOC(cwavef1,(2,size_wf))
 ABI_MALLOC(cwavef2,(2,size_wf))
 ABI_MALLOC(cg1_aux,(2,size_wf))
 ABI_MALLOC(gv1c,(2,size_wf))
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

 call cwtime(cpu, wall, gflops, "start")

!Specific definitions
 d3etot_t1_k=zero
 vlocal1=zero
 dum_lambda=zero
 berryopt=0;optlocal=0;optnl=1;usevnl=0;opt_gvnl1=0;sij_opt=0

!Initialize rf Hamiltonian (the k-dependent part is prepared in getgh1c_setup)
 call init_rf_hamiltonian(cplex,gs_hamkq,natom+1,rf_hamkq,& 
 & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
 & mpi_spintab=mpi_enreg%my_isppoltab)

 call rf_hamkq%load_spin(isppol,vlocal1=vlocal1,with_nonlocal=.true.)

 !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
 if (dtset%ffnl_lw==0) then
   ABI_MALLOC(ffnlk,(npw_k,0,psps%lmnmax,psps%ntypat))
   ABI_MALLOC(ffnl1,(npw_k,2,psps%lmnmax,psps%ntypat))
   ffnl1(:,1,:,:)=ffnl_k(:,1,:,:)
   ffnl1(:,2,:,:)=ffnl_k(:,1+i3dir,:,:)
 end if
 call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,&                        ! In
 kpt,kpt,i3dir,natom+1,natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,& ! In
 npw_k,npw_k,useylmgr1,kg_k,ylm_k,kg_k,ylm_k,part_ylmgr_k,&               ! In
 dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&           ! Out
 reuse_ffnlk=reuse_ffnlk, reuse_ffnl1=reuse_ffnl1)                        ! Optional

 !LOOP OVER BANDS
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle
   
   !Select bks wf1
   offset_cgi = (iband-1)*size_wf+icg
   cwavef1(:,:)= cg1(:,1+offset_cgi:size_wf+offset_cgi)
   cwavef2(:,:)= cg2(:,1+offset_cgi:size_wf+offset_cgi)

   !AZ_try_ini*****************************************************
   ! Print 1WF fro electric field perturbation
   !write(kpt_1_text,'(f6.3)') kpt(1)
   !write(kpt_2_text,'(f6.3)') kpt(2)
   !write(kpt_3_text,'(f6.3)') kpt(3)
   !write(iband_text,'(i5)') iband 
   !write(i1dir_text,'(i5)') i1dir
   !write(i2dir_text,'(i5)') i2dir

   !! Write cg1 wavefunction
   !file_name = 'AZ_cg1_iband_'//trim(adjustl(iband_text))//&
  !& '_iq2grad_'//trim(adjustl(i1dir_text))//&
  !& '_x_'//trim(adjustl(kpt_1_text))//&
  !& '_y_'//trim(adjustl(kpt_2_text))//&
  !& '_z_'//trim(adjustl(kpt_3_text))//'.dat'   

   !open(unit=999,file=file_name,action='write',status='replace')
   !do ii=1,size(cwavef1(1,:))
   !  write(999,'(i10,2f12.6)') ii, cwavef1(1,ii), cwavef1(2,ii)
   !enddo
   !close(999)

   ! Write cg2 wavefunction
   !file_name = 'AZ_cg2_iband_'//trim(adjustl(iband_text))//&
  !& '_iq2grad_'//trim(adjustl(i2dir_text))//&
  !& '_x_'//trim(adjustl(kpt_1_text))//&
  !& '_y_'//trim(adjustl(kpt_2_text))//&
  !& '_z_'//trim(adjustl(kpt_3_text))//'.dat'
  
   !open(unit=999,file=file_name,action='write',status='replace')
   !do ii=1,size(cwavef2(1,:))
   !  write(999,'(i10,2f12.6)') ii, cwavef2(1,ii), cwavef2(2,ii)
   !enddo
   !close(999)

   !AZ_try_fin*****************************************************
   
   !Compute < g |\partial_{gamma} H^{(0)} | u_{i,k}^{\lambda2} >
   call getgh1c(berryopt,cwavef2,dum_cwaveprj,gv1c,dum_grad_berry,&
 & dum_gs1,gs_hamkq,dum_gvnl1,i3dir,natom+1,dum_lambda,mpi_enreg,optlocal,&
 & optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
    
   !Apply the dot product with the ket wf (take into account occupation here)
   ! < u_{i,k}^{\lambda1}} | \partial_{gamma} H^{(0)} | u_{i,k}^{lambda2}} >
   call dotprod_g(dotr,doti,istwf_k,size_wf,2,cwavef1,gv1c, &
 & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   d3etot_t1_k(1)=d3etot_t1_k(1)+occ_k(iband)*dotr
   d3etot_t1_k(2)=d3etot_t1_k(2)+occ_k(iband)*doti

   !AZ_try_ini********************************************************
   write(i1dir_text,'(i8)') i1dir
   write(i2dir_text,'(i8)') i2dir
   write(i3dir_text,'(i8)') i3dir
   write(iband_text,'(i8)') iband
   write(kpt_1_text,'(f6.3)') kpt(1)
   write(kpt_2_text,'(f6.3)') kpt(2)
   write(kpt_3_text,'(f6.3)') kpt(3)
   file_name = 'AZ_T1_iband_'//trim(adjustl(iband_text))//&
 & '_iq1grad_'//trim(adjustl(i1dir_text))//&
 & '_iq2grad_'//trim(adjustl(i2dir_text))//&
 & '_iq3grad_'//trim(adjustl(i3dir_text))//&
 & '_x_'//trim(adjustl(kpt_1_text))//&
 & '_y_'//trim(adjustl(kpt_2_text))//&
 & '_z_'//trim(adjustl(kpt_3_text))//'.dat'
   open(unit=888,file=file_name,action='write',status='replace')
   write(888,'(2f12.6)') dotr*occ_k(iband)*wtk_k*eight*pi, doti*occ_k(iband)*wtk_k*eight*pi
   close(888)
   !AZ_try_fin********************************************************

 end do !iband

!Clean rf_hamiltonian
 call rf_hamkq%free()

 !Deallocations
 ABI_FREE(kpg_k)
 ABI_FREE(kpg1_k)
 ABI_FREE(dkinpw)
 ABI_FREE(kinpw1)
 ABI_FREE(ffnlk)
 ABI_FREE(ffnl1)
 ABI_FREE(ph3d)

 call cwtime(cpu, wall, gflops, "stop")

!------------------------------------T2------------------------------------------------
!q-gradient of CB projector x rf Hamiltonian lambda 2:
! < u_{i,k}^{\lambda1}} | \partial_{gamma} Q_k H^{\lambda2} | u_{i,k}^{(0)} >
!--------------------------------------------------------------------------------------

!Create array for ddk 1wf from file
 ABI_MALLOC(cg1_ddk,(2,size_wf,nband_k))
 do iband=1,nband_k
   call ddk_f%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=cg1_aux)
   cg1_ddk(:,:,iband)=cg1_aux(:,:)
 end do

 call cwtime(cpu, wall, gflops, "start")

!For \lambda1=\lambda2 T2 is inferred from the cc of T3
if (.not.samepert) then

  !Specific definitions
  d3etot_t2_k=zero

  !LOOP OVER BANDS
  do iband=1,nband_k
  
    if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle
  
    !Select bks wfs
    offset_cgi = (iband-1)*size_wf+icg
    cwavef1(:,:)= cg1(:,1+offset_cgi:size_wf+offset_cgi)
  
    !LOOP OVER BANDS
    do jband=1,nband_k
  
      !Select ddk wf1
      cg1_aux(:,:)=cg1_ddk(:,:,jband)

      !Load < u_{j,k}^{(0) | H^{\lambda2}+V^{\lambda2}} | u_{i,k}^{(0)} >
      ii=2*jband-1+(iband-1)*2*nband_k
      cj_h1_ci(1)=eig2_k(ii)
      cj_h1_ci(2)=eig2_k(ii+1)

      !Calculate: < u_{i,k}^{lambda1}} | u_{j,k}^{k_{\gamma}} >
      call dotprod_g(dotr,doti,istwf_k,size_wf,2,cwavef1,cg1_aux, &
    & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
  
      !Calculate the contribution to T2
      cprodr=dotr*cj_h1_ci(1)-doti*cj_h1_ci(2)
      cprodi=dotr*cj_h1_ci(2)+doti*cj_h1_ci(1)
      d3etot_t2_k(1)=d3etot_t2_k(1)-cprodr*occ_k(iband)
      d3etot_t2_k(2)=d3etot_t2_k(2)-cprodi*occ_k(iband)
  
    end do !jband 
  
  end do !iband
  
end if !samepert

 call cwtime(cpu, wall, gflops, "stop")

!------------------------------------T3------------------------------------------------
!rf Hamiltonian lambda 1 x q-gradient of CB projector
! < u_{i,k}^{(0) | (H^{\lambda1})^{\dagger} \partial_{gamma} Q_k | u_{i,k}^{\lambda2}}  >
!--------------------------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")
!Specific definitions
 d3etot_t3_k=zero

 !LOOP OVER BANDS
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle

   !Select bks wfs
   offset_cgi = (iband-1)*size_wf+icg
   cwavef2(:,:)= cg2(:,1+offset_cgi:size_wf+offset_cgi)

   !LOOP OVER BANDS
   do jband=1,nband_k

     !Select ddk wf1
     cg1_aux(:,:)=cg1_ddk(:,:,jband)

     !Load (< u_{j,k}^{(0) | H^{\lambda1}+V^{\lambda1}} | u_{i,k}^{(0)} >)^*
     ii=2*jband-1+(iband-1)*2*nband_k
     cj_h1_ci(1)=eig1_k(ii)
     cj_h1_ci(2)=-eig1_k(ii+1)

     !Calculate: < u_{j,k}^{k_{\gamma}} | u_{i,k}^{lambda2}} >
     call dotprod_g(dotr,doti,istwf_k,size_wf,2,cg1_aux,cwavef2, &
   & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

     !Calculate the contribution to T3
     cprodr=dotr*cj_h1_ci(1)-doti*cj_h1_ci(2)
     cprodi=dotr*cj_h1_ci(2)+doti*cj_h1_ci(1)
     d3etot_t3_k(1)=d3etot_t3_k(1)-cprodr*occ_k(iband)
     d3etot_t3_k(2)=d3etot_t3_k(2)-cprodi*occ_k(iband)

     !AZ_try_ini********************************************************
      write(i1dir_text,'(i8)') i1dir
      write(i2dir_text,'(i8)') i2dir
      write(i3dir_text,'(i8)') i3dir
      write(iband_text,'(i8)') iband
      write(jband_text,'(i8)') jband
      write(kpt_1_text,'(f6.3)') kpt(1)
      write(kpt_2_text,'(f6.3)') kpt(2)
      write(kpt_3_text,'(f6.3)') kpt(3)
      file_name = 'AZ_T3_iband_'//trim(adjustl(iband_text))//&
    & '_jband_'//trim(adjustl(jband_text))//&
    & '_iq1grad_'//trim(adjustl(i1dir_text))//&
    & '_iq2grad_'//trim(adjustl(i2dir_text))//&
    & '_iq3grad_'//trim(adjustl(i3dir_text))//&
    & '_x_'//trim(adjustl(kpt_1_text))//&
    & '_y_'//trim(adjustl(kpt_2_text))//&
    & '_z_'//trim(adjustl(kpt_3_text))//'.dat'
      open(unit=888,file=file_name,action='write',status='replace')
      write(888,'(2f12.6)') -cprodr*occ_k(iband)*wtk_k*eight*pi, -cprodi*occ_k(iband)*wtk_k*eight*pi
      close(888)
      !AZ_try_fin********************************************************

   end do !jband

 end do !iband

 ABI_FREE(cg1_ddk)
 ABI_FREE(cg1_aux)
 ABI_FREE(vpsp1)
 ABI_FREE(vlocal1)

if (samepert) then
  d3etot_t2_k(1)=d3etot_t3_k(1)
  d3etot_t2_k(2)=-d3etot_t3_k(2)
end if

 call cwtime(cpu, wall, gflops, "stop")
!------------------------------------T4------------------------------------------------
!q-gradient of rf Hamiltonian lambda 2 
! < u_{i,k}^{\lambda1} | H^{\lambda2}_{gamma} | u_{i,k}^{(0)} >
!--------------------------------------------------------------------------------------

!For \lambda1=\lambda2 T4 is inferred from the cc of T5
if (.not.samepert) then

  !Specific definitions and allocations
   d3etot_t4_k=zero
   optlocal=1;optnl=1
   dimffnlk=0
   if (i2pert/=natom+2) then
     ABI_MALLOC(vlocal1,(2*ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
     ABI_MALLOC(vpsp1,(2*nfft))
     ABI_MALLOC(gvloc1dqc,(2,size_wf))
     ABI_MALLOC(gvnl1dqc,(2,size_wf))
   end if
   if (i2pert<=natom) fac=-one
   if (i2pert==natom+2) fac=one
   if (i2pert==natom+3.or.i2pert==natom+4) fac=-half
   if (i2pert<=natom) then
     nylmgrtmp=3
     dimffnlk=1
     dimffnl1=2
   else if (i2pert==natom+3.or.i2pert==natom+4) then
     nylmgrtmp=nylmgr
     dimffnl1=10
     ABI_FREE(part_ylmgr_k)
     ABI_MALLOC(part_ylmgr_k,(npw_k,nylmgrtmp,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
     part_ylmgr_k(:,:,:)=ylmgr_k(:,:,:)
   end if
  
  !Do loop to compute both extradiagonal shear-strain components
   do idq=1,n2dq
  
 call cwtime(cpu, wall, gflops, "start")

     if (i2pert/=natom+2) then
       idir=i2dir; if (i2pert==natom+4) idir=idq*3+i2dir
       !Initialize rf Hamiltonian (the k-dependent part is prepared in getgh1c_setup)
       call init_rf_hamiltonian(2,gs_hamkq,i2pert,rf_hamkq,& 
       & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
       & mpi_spintab=mpi_enreg%my_isppoltab)
  
       !Set up local potentials with proper dimensioning
       !and load the spin-dependent part of the Hamiltonians
       vpsp1=vpsp1_i2pertdq(:,isppol,idq)
       call rf_transgrid_and_pack(isppol,nspden,usepaw,2,nfft,nfft,ngfft,&
       & gs_hamkq%nvloc,pawfgr,mpi_enreg,dum_vpsp,vpsp1,dum_vlocal,vlocal1)
       call rf_hamkq%load_spin(isppol,vlocal1=vlocal1,& 
       & with_nonlocal=with_nonlocal_i2pert)
  
       !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
       if (dtset%ffnl_lw==0) then
         ABI_MALLOC(ffnlk,(npw_k,dimffnlk,psps%lmnmax,psps%ntypat))
         if (dimffnlk==1) ffnlk(:,1,:,:)=ffnl_k(:,1,:,:)
         ABI_MALLOC(ffnl1,(npw_k,dimffnl1,psps%lmnmax,psps%ntypat))
         if (dimffnl1==2) then
           ffnl1(:,1,:,:)=ffnl_k(:,1,:,:)
           ffnl1(:,2,:,:)=ffnl_k(:,1+i3dir,:,:)
         else
           ffnl1(:,1:dimffnl1,:,:)=ffnl_k(:,1:dimffnl1,:,:)
         end if
       end if
       call getgh1dqc_setup(gs_hamkq,rf_hamkq,dtset,psps,kpt,kpt,idir,i2pert,i3dir, &
     & dtset%natom,rmet,rprimd,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,npw_k,npw_k,nylmgrtmp,useylmgr1,kg_k, &
     & ylm_k,kg_k,ylm_k,part_ylmgr_k,nkpg,nkpg1,kpg_k,kpg1_k,dkinpw,kinpw1,ffnlk,ffnl1,ph3d,ph3d1, &
     & reuse_ffnlk=reuse_ffnlk,reuse_ffnl1=reuse_ffnl1)
  
     end if
  
     !LOOP OVER BANDS
     do iband=1,nband_k
  
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle
  
       !Select bks wfs
       offset_cgi = (iband-1)*size_wf+icg
       cwavef1(:,:)= cg1(:,1+offset_cgi:size_wf+offset_cgi)
  
       !Perturbation-specific part
       if (i2pert==natom+2) then
         call d2_dkdk_f%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=gv1c)
       else
         cwave0i(:,:)= cg(:,1+offset_cgi:size_wf+offset_cgi)
  
         !Compute < g |H^{\lambda2}}_{\gamma} | u_{i,k}^{(0)} >
         call getgh1dqc(cwave0i,dum_cwaveprj,gv1c,gvloc1dqc,gvnl1dqc,gs_hamkq, &
         & idir,i2pert,mpi_enreg,optlocal,optnl,i3dir,rf_hamkq)
       end if
       
       !Calculate: < u_{j,k}^{\lambda1} | |H^{\lambda2}}_{\gamma} | u_{i,k}^{(0)} >
       call dotprod_g(dotr,doti,istwf_k,size_wf,2,cwavef1,gv1c, &
     & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
  
       !Calculate the contribution to T4
       d3etot_t4_k(1,idq)=d3etot_t4_k(1,idq)+dotr*occ_k(iband)
       d3etot_t4_k(2,idq)=d3etot_t4_k(2,idq)+doti*occ_k(iband)
  
     end do !iband
  
     if (i2pert/=natom+2) then
  
       !Clean the rf_hamiltonian
       call rf_hamkq%free()
  
       !Deallocations
       ABI_FREE(kpg_k)
       ABI_FREE(kpg1_k)
       ABI_FREE(dkinpw)
       ABI_FREE(kinpw1)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
  
     end if

 call cwtime(cpu, wall, gflops, "stop")
  
     !Apply the perturbation-dependent prefactors on T4
     tmpre=d3etot_t4_k(1,idq); tmpim=d3etot_t4_k(2,idq)
     if (i2pert<=natom.or.i2pert==natom+2) then
       d3etot_t4_k(1,idq)=-tmpim
       d3etot_t4_k(2,idq)=tmpre
     end if 
     d3etot_t4_k(:,idq)=d3etot_t4_k(:,idq)*fac
  
   end do !idq
  
   if (i2pert/=natom+2) then
     ABI_FREE(gvloc1dqc)
     ABI_FREE(gvnl1dqc)
     ABI_FREE(vlocal1)
     ABI_FREE(vpsp1)
   end if

 end if !samepert

!------------------------------------T5------------------------------------------------
!q-gradient of rf Hamiltonian lambda 1 
! < u_{i,k}^{(0)} | (H^{\lambda1}_{gamma})^{\dagger} | u_{i,k}^{\lambda2} >
!--------------------------------------------------------------------------------------

!Specific definitions and allocations
 d3etot_t5_k=zero
 optlocal=1;optnl=1
 dimffnlk=0
 if (i1pert/=natom+2) then
   ABI_MALLOC(vlocal1,(2*ngfft(4),ngfft(5),ngfft(6),gs_hamkq%nvloc))
   ABI_MALLOC(vpsp1,(2*nfft))
   ABI_MALLOC(gvloc1dqc,(2,size_wf))
   ABI_MALLOC(gvnl1dqc,(2,size_wf))
 end if
 if (i1pert<=natom) fac=-one
 if (i1pert==natom+2) fac=one
 if (i1pert==natom+3.or.i1pert==natom+4) fac=-half
 if (i1pert<=natom) then
   nylmgrtmp=3
   dimffnlk=1
   dimffnl1=2
   ABI_FREE(part_ylmgr_k)
   ABI_MALLOC(part_ylmgr_k,(npw_k,nylmgrtmp,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
   part_ylmgr_k(:,:,:)=ylmgr_k(:,1:3,:)
 else if (i1pert==natom+3.or.i1pert==natom+4) then
   nylmgrtmp=nylmgr
   dimffnl1=10
   ABI_FREE(part_ylmgr_k)
   ABI_MALLOC(part_ylmgr_k,(npw_k,nylmgrtmp,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
   part_ylmgr_k(:,:,:)=ylmgr_k(:,:,:)
 end if

!Do loop to compute both extradiagonal shear-strain components
 do idq=1,n1dq

 call cwtime(cpu, wall, gflops, "start")
   if (i1pert/=natom+2) then
     idir=i1dir; if (i1pert==natom+4) idir=idq*3+i1dir
     !Initialize rf Hamiltonian (the k-dependent part is prepared in getgh1c_setup)
     call init_rf_hamiltonian(2,gs_hamkq,i1pert,rf_hamkq,& 
     & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
     & mpi_spintab=mpi_enreg%my_isppoltab)

     !Set up local potentials with proper dimensioning
     !and load the spin-dependent part of the Hamiltonians
     vpsp1=vpsp1_i1pertdq(:,isppol,idq)
     call rf_transgrid_and_pack(isppol,nspden,usepaw,2,nfft,nfft,ngfft,&
     & gs_hamkq%nvloc,pawfgr,mpi_enreg,dum_vpsp,vpsp1,dum_vlocal,vlocal1)
     call rf_hamkq%load_spin(isppol,vlocal1=vlocal1,& 
     & with_nonlocal=with_nonlocal_i1pert)

     !Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
     if (dtset%ffnl_lw==0) then
       ABI_MALLOC(ffnlk,(npw_k,dimffnlk,psps%lmnmax,psps%ntypat))
       if (dimffnlk==1) ffnlk(:,1,:,:)=ffnl_k(:,1,:,:)
       ABI_MALLOC(ffnl1,(npw_k,dimffnl1,psps%lmnmax,psps%ntypat))
       if (dimffnl1==2) then
         ffnl1(:,1,:,:)=ffnl_k(:,1,:,:)
         ffnl1(:,2,:,:)=ffnl_k(:,1+i3dir,:,:)
       else
         ffnl1(:,1:dimffnl1,:,:)=ffnl_k(:,1:dimffnl1,:,:)
       end if
     end if
     call getgh1dqc_setup(gs_hamkq,rf_hamkq,dtset,psps,kpt,kpt,idir,i1pert,i3dir, &
   & dtset%natom,rmet,rprimd,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,npw_k,npw_k,nylmgrtmp,useylmgr1,kg_k, &
   & ylm_k,kg_k,ylm_k,part_ylmgr_k,nkpg,nkpg1,kpg_k,kpg1_k,dkinpw,kinpw1,ffnlk,ffnl1,ph3d,ph3d1, &
   & reuse_ffnlk=reuse_ffnlk,reuse_ffnl1=reuse_ffnl1)

   end if

   !LOOP OVER BANDS
   do iband=1,nband_k

     if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle

     !Select bks wfs
     offset_cgi = (iband-1)*size_wf+icg
     cwavef2(:,:)= cg2(:,1+offset_cgi:size_wf+offset_cgi)

     !Perturbation-specific part
     if (i1pert==natom+2) then
       call d2_dkdk_f%read_bks(iband,ikpt,isppol,xmpio_single,cg_bks=gv1c)
     else
       cwave0i(:,:)= cg(:,1+offset_cgi:size_wf+offset_cgi)

       !Compute < g |H^{\lambda1}}_{\gamma} | u_{i,k}^{(0)} >
       call getgh1dqc(cwave0i,dum_cwaveprj,gv1c,gvloc1dqc,gvnl1dqc,gs_hamkq, &
       & idir,i1pert,mpi_enreg,optlocal,optnl,i3dir,rf_hamkq)
     end if
     
     !Calculate: < u_{j,k}^{\lambda2} | |H^{\lambda1}}_{\gamma} | u_{i,k}^{(0)} >
     call dotprod_g(dotr,doti,istwf_k,size_wf,2,cwavef2,gv1c, &
   & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

     !Calculate the contribution to T5:
     d3etot_t5_k(1,idq)=d3etot_t5_k(1,idq)+dotr*occ_k(iband)
     d3etot_t5_k(2,idq)=d3etot_t5_k(2,idq)+doti*occ_k(iband)

     !AZ_try_ini********************************************************
     write(i1dir_text,'(i8)') i1dir
     write(i2dir_text,'(i8)') i2dir
     write(i3dir_text,'(i8)') i3dir
     write(iband_text,'(i8)') iband
     write(kpt_1_text,'(f6.3)') kpt(1)
     write(kpt_2_text,'(f6.3)') kpt(2)
     write(kpt_3_text,'(f6.3)') kpt(3)
     file_name = 'AZ_T5_iband_'//trim(adjustl(iband_text))//&
   & '_iq1grad_'//trim(adjustl(i1dir_text))//&
   & '_iq2grad_'//trim(adjustl(i2dir_text))//&
   & '_iq3grad_'//trim(adjustl(i3dir_text))//&
   & '_x_'//trim(adjustl(kpt_1_text))//&
   & '_y_'//trim(adjustl(kpt_2_text))//&
   & '_z_'//trim(adjustl(kpt_3_text))//'.dat'
     open(unit=888,file=file_name,action='write',status='replace')
     ! I don't care about the real part...
     write(888,'(2f12.6)') 0.0_dp, -dotr*occ_k(iband)*wtk_k*eight*pi
     close(888)
     !AZ_try_fin********************************************************


   end do !iband

   if (i1pert/=natom+2) then

     !Clean the rf_hamiltonian
     call rf_hamkq%free()

     !Deallocations
     ABI_FREE(kpg_k)
     ABI_FREE(kpg1_k)
     ABI_FREE(dkinpw)
     ABI_FREE(kinpw1)
     ABI_FREE(ffnlk)
     ABI_FREE(ffnl1)
     ABI_FREE(ph3d)

   end if

 call cwtime(cpu, wall, gflops, "stop")

   !Apply the perturbation-dependent prefactors on T5
   tmpre=d3etot_t5_k(1,idq); tmpim=d3etot_t5_k(2,idq)
   if (i1pert<=natom.or.i1pert==natom+2) then
     d3etot_t5_k(1,idq)=-tmpim
     d3etot_t5_k(2,idq)=tmpre
   end if 
   d3etot_t5_k(:,idq)=d3etot_t5_k(:,idq)*fac

   !Apply now the conjugate complex:
   !(< u_{j,k}^{\lambda2} | |H^{\lambda1}}_{\gamma} | u_{i,k}^{(0)} >)* 
   tmpim=d3etot_t5_k(2,idq)
   d3etot_t5_k(2,idq)=-tmpim

 end do !idq


 if (i1pert/=natom+2) then
   ABI_FREE(gvloc1dqc)
   ABI_FREE(gvnl1dqc)
   ABI_FREE(vlocal1)
   ABI_FREE(vpsp1)
 end if

if (samepert) then
  d3etot_t4_k(1,:)=d3etot_t5_k(1,:)
  d3etot_t4_k(2,:)=-d3etot_t5_k(2,:)
end if

!Scale d3etot_k contributions by the kpt weight
d3etot_t1_k(:)=d3etot_t1_k(:)*wtk_k
d3etot_t2_k(:)=d3etot_t2_k(:)*wtk_k
d3etot_t3_k(:)=d3etot_t3_k(:)*wtk_k
d3etot_t4_k(:,:)=d3etot_t4_k(:,:)*wtk_k
d3etot_t5_k(:,:)=d3etot_t5_k(:,:)*wtk_k

!Deallocations
 ABI_FREE(cwave0i)
 ABI_FREE(cwave0j)
 ABI_FREE(cwavef1)
 ABI_FREE(cwavef2)
 ABI_FREE(gv1c)
 ABI_FREE(dum_vpsp)
 ABI_FREE(dum_vlocal)
 ABI_FREE(dum_cwaveprj)
 ABI_FREE(part_ylmgr_k)


 DBG_EXIT("COLL")

end subroutine dfpt_1wf
!!***

end module m_dfptlw_wf
!!***
