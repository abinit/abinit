!!****m* ABINIT/m_fock_getghc
!! NAME
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2025 ABINIT group (CMartins, FJ, MT, XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_fock_getghc

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_fock
 use m_pawcprj
 !use m_cgtools
 use, intrinsic :: iso_c_binding, only: c_size_t

 use defs_abitypes, only : mpi_type
 use defs_datatypes, only : pseudopotential_type
 use m_time,         only : timab, time_accu
 use m_matrix,       only : matr3inv
 use m_cgtools,      only : dotprod_g,dotprod_g_batch_half
 use m_kg,           only : mkkpg
 use m_fftcore,      only : sphereboundary
 use m_fft,          only : fftpac, fourwf, fourdp
 use m_fstrings,     only : sjoin, itoa
 use m_hamiltonian,  only : gs_hamiltonian_type, K_H_KPRIME
 use m_paw_nhat,     only : pawmknhat_psipsi, pawdijhat_ndat
 use m_spacepar,     only : hartre
 use m_nonlop,       only : nonlop
 use m_bandfft_kpt,      only : bandfft_kpt, bandfft_kpt_type, bandfft_kpt_savetabs,bandfft_kpt_restoretabs, &
                                prep_bandfft_tabs
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type
 use m_mkffnl,           only : mkffnl
 use m_mpinfo,           only : proc_distrb_cycle
 use m_abi_linalg

 use, intrinsic :: iso_c_binding, only: c_loc

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

 implicit none

 private
!!***

 public :: fock_getghc
 public :: fock2ACE
 public :: fock_ACE_getghc
!!***

contains
!!***

subroutine select_ndat_occ_for_gpu(ndat_occ,nband_k,ndat,npw,cplex_fock,nfftf,&
    n4,n5,n6,natom,nspinor,lmn2_size,usepaw,cprj,ieigen,need_ghc,optfor,optstr)

!Arguments ------------------------------------
! Scalars
 integer,intent(in)     :: nband_k,ndat,npw,cplex_fock,nfftf,n4,n5,n6
 integer,intent(in)     :: natom,nspinor,lmn2_size,usepaw,ieigen
 logical,intent(in)     :: optfor,optstr,need_ghc
 integer,intent(out)    :: ndat_occ
 type(pawcprj_type),intent(in) :: cprj(natom,nspinor*ndat)

!Local variables-------------------------------
 integer :: i,ider,nprojs
 integer(kind=c_size_t) :: sum_mem,free_mem

! *************************************************************************

 ider=0
 if (usepaw==1) then
   nprojs=0
   do i = 1,natom
     nprojs = nprojs + cprj(i, 1)%nlmn
   end do
 end if

#ifdef HAVE_GPU
 call gpu_get_max_mem(free_mem)
 free_mem = 0.95 * free_mem ! Cutting 5% out to be safe
#endif

 do i=1,nband_k
   if(modulo(nband_k,i)/=0) cycle
   ndat_occ=nband_k/i
   sum_mem = 0

   ! cwavef
   sum_mem = sum_mem + INT(2,c_size_t)*n4*n5*n6*ndat
   if(need_ghc) then
     ! ghc
     sum_mem = sum_mem + INT(2,c_size_t)*npw*ndat
     ! ghc1
     sum_mem = sum_mem + INT(2,c_size_t)*npw*ndat
     ! ghc2
     sum_mem = sum_mem + INT(2,c_size_t)*npw*ndat
   end if

   ! cwavef_r
   sum_mem = sum_mem + INT(2,c_size_t)*n4*n5*n6*ndat
   ! vlocpsi_r
   sum_mem = sum_mem + INT(cplex_fock,c_size_t)*nfftf*ndat
   ! work (ompgpu_fourwf internal array)
   sum_mem = sum_mem + INT(2,c_size_t)*n4*n5*n6*ndat

   ! rhor_munu
   sum_mem = sum_mem + INT(cplex_fock,c_size_t)*nfftf*ndat_occ*ndat
   ! rhog_munu
   sum_mem = sum_mem + INT(2,c_size_t)*nfftf*ndat_occ*ndat
   ! vfock
   sum_mem = sum_mem + INT(cplex_fock,c_size_t)*nfftf*ndat_occ*ndat
   ! occ
   sum_mem = sum_mem + INT(ndat_occ,c_size_t)

  !*Additional arrays in case of paw
   if (usepaw==1) then
     if ((optfor).and.(ieigen/=0)) then
       ider=3
       ! forout
       sum_mem = sum_mem + INT(2,c_size_t)*npw*nspinor*ndat_occ

       ! dprojs (gemm_nonlop_ompgpu internal work array)
       sum_mem = sum_mem + INT(2,c_size_t)*nprojs*npw*3
     end if

     if (optstr.and.(ieigen/=0)) then
       ider=3
       ! strout
       sum_mem = sum_mem + INT(2,c_size_t)*npw*nspinor*ndat_occ

       ! dprojs (gemm_nonlop_ompgpu internal work array)
       sum_mem = sum_mem + INT(2,c_size_t)*nprojs*npw*6
     end if
     ! grnhat_12
     sum_mem = sum_mem + INT(2,c_size_t)*nfftf*nspinor**2*3*natom*(ider/3)*ndat_occ*ndat
     ! gvnlxc
     sum_mem = sum_mem + INT(2,c_size_t)*npw*nspinor*ndat_occ
     ! rho12
     sum_mem = sum_mem + INT(2,c_size_t)*nfftf*nspinor**2*ndat_occ*ndat

     ! rho12 (paw_psipsi internal work array)
     sum_mem = sum_mem + INT(2,c_size_t)*nfftf*nspinor**2*ndat_occ*ndat*natom
     ! cprj1 (paw_psipsi internal work array)
     sum_mem = sum_mem + INT(2,c_size_t)*nprojs*nspinor*ndat
     ! cprj2 (paw_psipsi internal work array)
     sum_mem = sum_mem + INT(2,c_size_t)*nprojs*nspinor*ndat
     ! cpf  (paw_psipsi internal work array)
     sum_mem = sum_mem + INT(2,c_size_t)*lmn2_size*nspinor*ndat*ndat_occ*natom

     ! dijhat
     sum_mem = sum_mem + INT(2,c_size_t)*lmn2_size*nspinor*ndat*ndat_occ*natom

     ! projs (gemm_nonlop_ompgpu internal work array)
     sum_mem = sum_mem + INT(2,c_size_t)*nprojs*npw
   end if

   ! cwaveocc_r
   sum_mem = sum_mem + INT(2,c_size_t)*n4*n5*n6*ndat_occ

   sum_mem = sum_mem*dp

   if(sum_mem < free_mem) exit
 end do

 if(sum_mem > free_mem) then
   ABI_WARNING("Test case doesn't fit in GPU memory. Try to lower bandpp.")
 end if
!#ifdef DEBUG_VERBOSE
 write(std_out,*) "-----------DEBUG fock_getghc%select_ndat_occ_for_gpu : "
 write(std_out,'(A,F10.3,1x,A)') "Free GPU memory                : ", real(free_mem,dp)/(1024*1024), "MiB"
 write(std_out,'(A,I4)')         "selected ndat_occ              : ", ndat_occ
 write(std_out,'(A,F10.3,1x,A)') "Forecasted consumed GPU memory : ", real(sum_mem,dp)/(1024*1024), "MiB"
 write(std_out,*) "-----------END DEBUG fock_getghc%select_ndat_occ_for_gpu : "
!#endif

 end subroutine select_ndat_occ_for_gpu

!!****f* ABINIT/fock_getghc
!! NAME
!!  fock_getghc
!!
!! FUNCTION
!!  Compute the matrix elements <G|Vx|psi> of the Fock operator.
!!
!! INPUTS
!!  cwavef(2,npw*nspinor*ndat)= planewave coefficients of wavefunctions on which Fock operator is applied.
!!  cwaveprj <type(pawcprj_type> = <cwavevf|proj>
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!!  mpi_enreg= information about MPI parallelization
!!
!! SIDE EFFECTS
!!  ghc(2,npw*ndat)= matrix elements <G|H|C> or <G|H-lambda.S|C> (if sij_opt>=0 or =-1 in getghc)
!!                   contains the fock exchange term for cwavef at the end.
!!
!! NOTES
!!  The current version assumes that:
!!   * nspinor = 1
!!   * no "my_nspinor"
!!   * no restriction to the value of istwfk_bz (but must be tested in all case)
!!   * all the data for the occupied states (cgocc_bz) are the same as those for the current states (cg)
!!
!! SOURCE

subroutine fock_getghc(cwavef,cwaveprj,ghc,gs_ham,mpi_enreg,ndat)

!Arguments ------------------------------------
! Scalars
 integer,intent(in)     :: ndat
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)
 real(dp),intent(inout) :: cwavef(:,:)!,ghc(2,gs_ham%npw_k)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 integer,parameter :: tim_fourwf_fock_getghc=10,tim_fourdp_fock_getghc=10
 integer :: bdtot_jindex,choice,cplex_fock,cplex_dij,cpopt,i1,i2,i3,ia,iatom,idat,idat_occ,iatm
 integer :: iband_cprj,ider,idir,idir1,ier,ii,ind,ipw,ifft,itypat,izero,jband,jbg,jcg,jkg
 integer :: jkpt,my_jsppol,jstwfk,lmn2_size,mgfftf,mpw,n1,n2,n3,n4,n5,n6,ndat_occ
 integer :: n1f,n2f,n3f,n4f,n5f,n6f,natom,nband_k,ndij,nfft,nfftf,nfftotf,nhat12_grdim,nnlout
 integer :: npw,npwj,nspden_fock,nspinor,nkpg,paw_opt,signs,tim_nonlop,gpu_option
 integer, save :: ncount=0
 logical :: need_ghc,qeq0
 real(dp),parameter :: weight1=one
 real(dp) :: doti,eigen,imcwf,imcwocc,imvloc,invucvol,recwf,recwocc,revloc,wtk
 complex(dp) :: cinvucvol,cucvol
 type(fock_common_type),pointer :: fockcommon
 type(fock_BZ_type),pointer :: fockbz
! Arrays
 integer :: ngfft(18),ngfftf(18)
 integer,pointer :: gboundf(:,:),kg_occ(:,:),gbound_kp(:,:)
 real(dp) :: fockstr(6),qphon(3),qvec_j(3),tsec(2),gsc_dum(2,0),rhodum(2,1)
 real(dp) :: rhodum0(0,1,1)
 real(dp), allocatable :: dummytab(:,:),dijhat(:,:,:,:,:,:),dijhat_tmp(:,:,:),ffnl_kp_dum(:,:,:,:),kpg_kp(:,:),occ(:)
 real(dp), allocatable, target :: gvnlxc(:,:),ghc1(:,:),ghc2(:,:),grnhat12(:,:,:,:,:,:),grnhat_12(:,:,:,:,:,:,:),forikpt(:,:,:)
 real(dp), allocatable :: rho12(:,:,:,:,:),rhog_munu(:,:,:,:),rhor_munu(:,:,:,:),vlocpsi_r(:,:),strdat(:,:,:,:)
 real(dp), allocatable :: vfock(:,:,:),psilocal(:,:,:),enlout_dum(:),vectin_dum(:,:),vqg(:),forout(:,:),strout(:,:),for1(:,:,:,:)
 real(dp), allocatable,target ::cwavef_r(:,:,:,:),vdotr(:,:,:,:),vdoti(:),vfockstr(:,:,:)
 real(dp), ABI_CONTIGUOUS  pointer :: cwaveocc_r(:,:,:,:,:)
 type(pawcprj_type),pointer :: cwaveocc_prj(:,:)

 real(dp) :: rprimd(3,3),for12(3),esum
 integer,  ABI_CONTIGUOUS pointer :: atom_ifftsph(:,:),atom_nfgd(:)
 real(dp), ABI_CONTIGUOUS pointer :: stress_ikpt(:,:),atom_rfgd(:,:,:)
 integer   :: ieigen


! *************************************************************************
!return

 ncount=ncount+1

 call timab(1504,1,tsec) ; call timab(1505,-1,tsec) ; call timab(1515,-1,tsec) ; call timab(1541,-1,tsec)

 ABI_CHECK(associated(gs_ham%fockcommon),"fock_common must be associated!")
 fockcommon => gs_ham%fockcommon
 ABI_CHECK(associated(gs_ham%fockbz),"fock_bz must be associated!")
 fockbz => gs_ham%fockbz

 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")
 if (fockcommon%usepaw==1) then
   ABI_CHECK((size(cwaveprj,1)==gs_ham%natom.and.size(cwaveprj,2)==gs_ham%nspinor*ndat),"error on cwaveprj dims")
 end if
 need_ghc=(size(ghc,2)>0)

!Some constants
 invucvol=1.d0/sqrt(gs_ham%ucvol)
 cinvucvol=dcmplx(invucvol,0.0_dp)
 cucvol=dcmplx(sqrt(gs_ham%ucvol),0.0_dp)
 call matr3inv(gs_ham%gprimd,rprimd)
 cplex_fock=2;nspden_fock=1
 natom=fockcommon%natom
 nspinor=gs_ham%nspinor
 mpw=maxval(fockbz%npwarr)
 npw=gs_ham%npw_k
 gpu_option=gs_ham%gpu_option
 ider=0;izero=0
 if (fockcommon%usepaw==1) then
   nfft =fockcommon%pawfgr%nfftc ; ngfft =fockcommon%pawfgr%ngfftc
   nfftf=fockcommon%pawfgr%nfft  ; ngfftf=fockcommon%pawfgr%ngfft
   mgfftf=fockcommon%pawfgr%mgfft
 else
   nfft =gs_ham%nfft  ; nfftf =nfft
   ngfft=gs_ham%ngfft ; ngfftf=ngfft
   mgfftf=gs_ham%mgfft
 end if
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 n4=ngfft(4);n5=ngfft(5);n6=ngfft(6)
 n1f=ngfftf(1);n2f=ngfftf(2);n3f=ngfftf(3)
 n4f=ngfftf(4);n5f=ngfftf(5);n6f=ngfftf(6)

! ===========================
! === Initialize arrays   ===
! ===========================
! transient optfor and optstress
! fockcommon%optfor=.false.
! fockcommon%optstr=.false.
!*Initialization of local pointers
!*Initialization of the array cwavef_r
!*cwavef_r = current wavefunction in r-space
 ABI_MALLOC(cwavef_r,(2,n4f,n5f,n6f*ndat))
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(alloc:cwavef_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
!*dummytab = variables for fourwf
 ABI_MALLOC(dummytab,(2,nfft*ndat))
!*vqg = 4pi/(G+q)**2
 ABI_MALLOC(vqg,(nfftf))

 if(need_ghc) then
  !*Initialization of the array ghc1
  !*ghc1 will contain the exact exchange contribution to the Hamiltonian
   ABI_MALLOC(ghc1,(2,npw*ndat))
   ABI_MALLOC(ghc2,(2,npw*ndat))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:ghc1,ghc2,vqg) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   if(gpu_option==ABI_GPU_DISABLED) then
     ghc1=zero
     ghc2=zero
   else if(gpu_option==ABI_GPU_OPENMP) then
     call gpu_set_to_zero(ghc1, int(2,c_size_t)*npw*ndat)
     call gpu_set_to_zero(ghc2, int(2,c_size_t)*npw*ndat)
   end if
 end if
!*Initialization of the array vlocpsi_r
!*vlocpsi_r = partial local Fock operator applied to cwavef in r-space and summed over all occupied (jkpt,mu)
 ABI_MALLOC(vlocpsi_r,(cplex_fock*nfftf,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(alloc:vlocpsi_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
 if(gpu_option==ABI_GPU_DISABLED) then
   vlocpsi_r=zero
 else if(gpu_option==ABI_GPU_OPENMP) then
   call gpu_set_to_zero(vlocpsi_r, int(cplex_fock,c_size_t)*nfftf*ndat)
 end if

!*Additional arrays in case of paw
 if (fockcommon%usepaw==1) then
   nhat12_grdim=0
 end if

 if (fockcommon%usepaw==1) then
   if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
     ABI_MALLOC(forikpt,(3,natom,ndat))
     forikpt=zero
   end if
 end if
 if (fockcommon%usepaw==1.or.fockcommon%optstr) then
   ABI_MALLOC(gboundf,(2*mgfftf+8,2))
   call sphereboundary(gboundf,gs_ham%istwf_k,gs_ham%kg_k,mgfftf,npw)
 else
   gboundf=>gs_ham%gbound_k
 end if
! ==========================================
! === Get cwavef in real space using FFT ===
! ==========================================
 if(gpu_option==ABI_GPU_DISABLED) then
   cwavef_r=zero
 else if(gpu_option==ABI_GPU_OPENMP) then
   call gpu_set_to_zero(cwavef_r, int(2,c_size_t)*n4f*n5f*n6f*ndat)
 end if
 call timab(1515,2,tsec) ; call timab(1541,-2,tsec) ; call timab(1512,-1,tsec)
 call fourwf(0,rhodum0,cwavef,rhodum,cwavef_r,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
& mgfftf,mpi_enreg,ndat,ngfftf,npw,1,n4f,n5f,n6f,0,tim_fourwf_fock_getghc,weight1,weight1,&
& gpu_option=gs_ham%gpu_option)
 call timab(1512,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1541,-1,tsec)
 if(gpu_option==ABI_GPU_DISABLED) then
   cwavef_r=cwavef_r*invucvol
 else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_ADDR(cwavef_r)
   call abi_gpu_xscal(2,n4f*n5f*n6f*ndat,cinvucvol,c_loc(cwavef_r),1)
   !$OMP END TARGET DATA
#endif
 end if

! =====================================================
! === Select the states in cgocc_bz with the same spin ===
! =====================================================
!* Initialization of the indices/shifts, according to the value of isppol
!* bdtot_jindex = shift to be applied on the location of data in the array occ_bz ?
 bdtot_jindex=0
!* jbg = shift to be applied on the location of data in the array cprj/occ
 jbg=0;jcg=0
 my_jsppol=fockcommon%isppol
 if((fockcommon%isppol==2).and.(mpi_enreg%nproc_spkpt/=1)) my_jsppol=1

!===================================
!=== Loop on the k-points in IBZ ===
!===================================
 jkg=0

 if (associated(gs_ham%ph3d_kp)) then
   nullify (gs_ham%ph3d_kp)
 end if

 call timab(1505,2,tsec) ; call timab(1506,-1,tsec) ; call timab(1541,-2,tsec)

 do jkpt=1,fockbz%mkpt

   if(fockbz%nbandocc_bz(jkpt,my_jsppol)==0) cycle

   call timab(1521,1,tsec)

!* nband_k = number of bands at point k_j
   nband_k=fockbz%nbandocc_bz(jkpt,my_jsppol)

   ! Select ndat_occ value :
   ! GPU : check GPU memory available and compute maximum value for ndat_occ
   ! CPU : 4, because it seems to be optimal from my observations
   if(gpu_option/=ABI_GPU_DISABLED) then
     lmn2_size=0
     if(fockcommon%usepaw==1) lmn2_size=fockcommon%pawtab(1)%lmn2_size
     call select_ndat_occ_for_gpu(ndat_occ,nband_k,ndat,npw,cplex_fock,&
&        nfftf,n4f,n5f,n6f,natom,nspinor,lmn2_size,&
&        fockcommon%usepaw,cwaveprj,fockcommon%ieigen,need_ghc,fockcommon%optfor,fockcommon%optstr)
   else
     ndat_occ=min(nband_k,4)
     do ii=1,nband_k
       if(modulo(nband_k,ndat_occ)==0) exit
       ndat_occ=ndat_occ-1
     end do
   end if

!* wtk = weight in BZ of this k point
   wtk=fockbz%wtk_bz(jkpt) !*sqrt(gs_ham%ucvol)
!* jstwfk= how is stored the wavefunction
   jstwfk=fockbz%istwfk_bz(jkpt)
!* npwj= number of plane wave in basis for the wavefunction
   npwj=fockbz%npwarr(jkpt)
!* Basis sphere of G vectors
   if (allocated(fockbz%cgocc)) then
     gbound_kp => fockbz%gbound_bz(:,:,jkpt)
     kg_occ => fockbz%kg_bz(:,1+jkg:npwj+jkg)
   end if

  !*rhormunu = overlap matrix between cwavef and (jkpt,mu) in R-space
   ABI_MALLOC(rhor_munu,(cplex_fock,nfftf,ndat_occ,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:rhor_munu) IF(gpu_option==ABI_GPU_OPENMP)
#endif
  !*rhogmunu = overlap matrix between cwavef and (jkpt,mu) in G-space
   ABI_MALLOC(rhog_munu,(2,nfftf,ndat_occ,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:rhog_munu) IF(gpu_option==ABI_GPU_OPENMP)
#endif
  !*vfock = Fock potential
   ABI_MALLOC(vfock,(cplex_fock*nfftf,ndat_occ,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:vfock) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_MALLOC(occ,(ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:occ) IF(gpu_option==ABI_GPU_OPENMP)
#endif

  !*Additional arrays in case of paw
   if (fockcommon%usepaw==1) then
     if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
       ider=3
       ABI_MALLOC(forout,(2,npw*nspinor*ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET ENTER DATA MAP(alloc:forout) IF(gpu_option==ABI_GPU_OPENMP)
#endif
     end if

     if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
       ider=3
       ABI_MALLOC(strout,(2,npw*nspinor*ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET ENTER DATA MAP(alloc:strout) IF(gpu_option==ABI_GPU_OPENMP)
#endif
     end if
     ABI_MALLOC(grnhat_12,(2,nfftf,nspinor**2,3,natom*(ider/3),ndat_occ,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET ENTER DATA MAP(alloc:grnhat_12) IF(gpu_option==ABI_GPU_OPENMP .and. ider==3)
#endif
     ABI_MALLOC(gvnlxc,(2,npw*nspinor*ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET ENTER DATA MAP(alloc:gvnlxc) IF(gpu_option==ABI_GPU_OPENMP)
#endif
     ABI_MALLOC(grnhat12,(2,nfftf,nspinor**2,3*nhat12_grdim,ndat_occ,ndat))
     ABI_MALLOC(rho12,(2,nfftf,nspinor**2,ndat_occ,ndat))
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET ENTER DATA MAP(alloc:rho12) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   end if

!* Load k^prime hamiltonian in the gs_ham datastructure
!  Note: ffnl_kp / ph3d_kp / gbound_kp are not used

   if (.not. associated(gs_ham%ph3d_kp)) then
     ABI_MALLOC(gs_ham%ph3d_kp,(2,npwj,gs_ham%matblk))
   end if

   call gs_ham%load_kprime(kpt_kp=fockbz%kptns_bz(:,jkpt),&
&   istwf_kp=jstwfk,npw_kp=npwj,kg_kp=fockbz%kg_bz(:,1+jkg:npwj+jkg))
!* Some temporary allocations needed for PAW
   if (fockcommon%usepaw==1) then
     ABI_MALLOC(enlout_dum,(ndat_occ))
     ABI_MALLOC(vectin_dum,(2,npwj*nspinor*ndat_occ))
     vectin_dum=zero
     ABI_MALLOC(ffnl_kp_dum,(npwj,1,gs_ham%lmnmax,gs_ham%ntypat))
     nkpg=size(gs_ham%kpg_k,2)
     ABI_MALLOC(kpg_kp,(npwj,nkpg))
     if (nkpg>0) then
       call mkkpg(gs_ham%kg_kp,kpg_kp,gs_ham%kpt_kp,nkpg,npwj)
     end if
     call gs_ham%load_kprime(ffnl_kp=ffnl_kp_dum,kpg_kp=kpg_kp)
   end if

! ======================================
! === Calculate the vector q=k_i-k_j ===
! ======================================
!* Evaluation of kpoint_j, the considered k-point in reduced coordinates
!     kpoint_j(:)=fockbz%kptns_bz(:,jkpt)
!* the vector qvec is expressed in reduced coordinates.
!     qvec(:)=kpoint_i(:)-kpoint_j(:)
   qvec_j(:)=gs_ham%kpt_k(:)-fockbz%kptns_bz(:,jkpt)
   qeq0=(qvec_j(1)**2+qvec_j(2)**2+qvec_j(3)**2<1.d-15)

   ! Get the Coulomb interaction in reciprocal space
   call bare_vqg(qvec_j,fockcommon,gs_ham%gmet,nfftf,fockbz%nkpt_bz,ngfftf,gs_ham%ucvol,vqg)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE TO(vqg) IF(gpu_option==ABI_GPU_OPENMP)
#endif

   call timab(1521,2,tsec)

! =================================================
! === Loop on the band indices jband of cgocc_k ===
! =================================================
   do jband=1,nband_k,ndat_occ

!*   occ = occupancy of jband at this k point
     occ(1:ndat_occ)=fockbz%occ_bz(jband+bdtot_jindex:jband+ndat_occ-1+bdtot_jindex,my_jsppol)
     if(occ(1)<tol8) cycle
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE TO(occ) IF(gpu_option==ABI_GPU_OPENMP)
#endif

!    This timing is placed after the cycle ...
     call timab(1522,1,tsec) ; call timab(1542,-1,tsec)

! ==============================================
! === Get cwaveocc_r in real space using FFT ===
! ==============================================
     if (allocated(fockbz%cwaveocc_bz)) then
       cwaveocc_r => fockbz%cwaveocc_bz(:,:,:,:,jband+jbg:jband+jbg+ndat_occ-1,my_jsppol)
     else
       ABI_MALLOC(cwaveocc_r,(2,n4f,n5f,n6f,ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET ENTER DATA MAP(alloc:cwaveocc_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       if(gpu_option==ABI_GPU_DISABLED) then
         cwaveocc_r=zero
       else if(gpu_option==ABI_GPU_OPENMP) then
         call gpu_set_to_zero(cwaveocc_r, int(2,c_size_t)*n4f*n5f*n6f*ndat_occ)
       end if
       call timab(1515,2,tsec) ; call timab(1512,-1,tsec) ; call timab(1542,-2,tsec)
       call fourwf(1,rhodum0,fockbz%cgocc(:,1+jcg+npwj*(jband-1):jcg+(jband+ndat_occ-1)*npwj,my_jsppol),rhodum,cwaveocc_r, &
&       gbound_kp,gbound_kp,jstwfk,kg_occ,kg_occ,mgfftf,mpi_enreg,ndat_occ,ngfftf,&
&       npwj,1,n4f,n5f,n6f,0,tim_fourwf_fock_getghc,weight1,weight1,gpu_option=gs_ham%gpu_option)
       call timab(1512,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1542,-1,tsec)
       if(gpu_option==ABI_GPU_DISABLED) then
         cwaveocc_r=cwaveocc_r*invucvol
       else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET DATA USE_DEVICE_ADDR(cwaveocc_r)
         call abi_gpu_xscal(2,n4f*n5f*n6f*ndat_occ,cinvucvol,c_loc(cwaveocc_r),1)
         !$OMP END TARGET DATA
#endif
       end if
     end if

! ================================================
! === Get the overlap density matrix rhor_munu ===
! ================================================
!* Calculate the overlap density matrix in real space = conj(cwaveocc_r)*cwavef_r
!* rhor_munu will contain the overlap density matrix.
! vfock=-int{conj(cwaveocc_r)*cwavef_r*dr'/|r-r'|}

     call timab(1522,2,tsec) ; call timab(1542,-2,tsec) ; call timab(1523,-1,tsec)

     if(gpu_option==ABI_GPU_DISABLED) then
       !$OMP PARALLEL DO COLLAPSE(2) &
       !$OMP& PRIVATE(ind,imcwf,recwf,recwocc,imcwocc)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
         do i3=1,n3f
           do i2=1,n2f
             do i1=1,n1f
               ind=i1+(i2-1)*n1f+(i3-1)*n2f*n1f
               recwf  =cwavef_r(1,i1,i2,(idat-1)*n3f+i3)
               imcwf  =cwavef_r(2,i1,i2,(idat-1)*n3f+i3)
               recwocc=cwaveocc_r(1,i1,i2,i3,idat_occ)
               imcwocc=cwaveocc_r(2,i1,i2,i3,idat_occ)
               rhor_munu(1,ind,idat_occ,idat)= recwocc*recwf+imcwocc*imcwf
               rhor_munu(2,ind,idat_occ,idat)= recwocc*imcwf-imcwocc*recwf
             end do ! i1
           end do ! i2
         end do ! i3
       end do ! idat_occ
       end do ! idat
     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& MAP(to:rhor_munu,cwavef_r,cwaveocc_r) PRIVATE(idat,idat_occ)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
         !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ind,i3,i2,i1)
         do i3=1,n3f
           do i2=1,n2f
             do i1=1,n1f
               ind=i1+(i2-1)*n1f+(i3-1)*n2f*n1f
               rhor_munu(1,ind,idat_occ,idat)= cwaveocc_r(1,i1,i2,i3,idat_occ)*cwavef_r(1,i1,i2,(idat-1)*n3f+i3)&
                                              +cwaveocc_r(2,i1,i2,i3,idat_occ)*cwavef_r(2,i1,i2,(idat-1)*n3f+i3)
               rhor_munu(2,ind,idat_occ,idat)= cwaveocc_r(1,i1,i2,i3,idat_occ)*cwavef_r(2,i1,i2,(idat-1)*n3f+i3)&
                                              -cwaveocc_r(2,i1,i2,i3,idat_occ)*cwavef_r(1,i1,i2,(idat-1)*n3f+i3)
             end do ! i1
           end do ! i2
         end do ! i3
       end do ! idat_occ
       end do ! idat
#endif
     end if ! gpu_option

     call timab(1523,2,tsec)

! =======================================================
! === Add compensation charge density in the PAW case ===
! =======================================================

     call timab(1524,-1,tsec) ; call timab(1544,-1,tsec)

     if (fockcommon%usepaw==1) then

       iband_cprj=(my_jsppol-1)*fockbz%mkptband+jbg+jband
       cwaveocc_prj=>fockbz%cwaveocc_prj(:,iband_cprj:iband_cprj+ndat_occ*nspinor-1)

       call pawmknhat_psipsi(cwaveprj(:,:),cwaveocc_prj(:,:),&
&       ider,izero,natom,natom,nfftf,ngfftf,&
&       nhat12_grdim,nspinor,fockcommon%ntypat,ndat,ndat_occ,fockbz%pawang,fockcommon%pawfgrtab,grnhat12,&
&       rho12,&
&       fockcommon%pawtab,gprimd=gs_ham%gprimd,grnhat_12=grnhat_12,qphon=qvec_j,&
&       xred=gs_ham%xred,atindx=gs_ham%atindx,gpu_option=gpu_option,nattyp=gs_ham%nattyp)

       if(gpu_option==ABI_GPU_DISABLED) then
         !$OMP PARALLEL DO COLLAPSE(2) &
         !$OMP& PRIVATE(idat,idat_occ)
         do idat=1,ndat
         do idat_occ=1,ndat_occ
           rhor_munu(1,:,idat_occ,idat)=rhor_munu(1,:,idat_occ,idat)+rho12(1,:,nspinor,idat_occ,idat)
           rhor_munu(2,:,idat_occ,idat)=rhor_munu(2,:,idat_occ,idat)-rho12(2,:,nspinor,idat_occ,idat)
         end do ! idat_occ
         end do ! idat
       else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
         !$OMP& MAP(to:rhor_munu,rho12) PRIVATE(idat,idat_occ)
         do idat=1,ndat
         do idat_occ=1,ndat_occ
         !$OMP PARALLEL DO PRIVATE(ifft)
         do ifft=1,nfftf
           rhor_munu(1,ifft,idat_occ,idat)=rhor_munu(1,ifft,idat_occ,idat)+rho12(1,ifft,nspinor,idat_occ,idat)
           rhor_munu(2,ifft,idat_occ,idat)=rhor_munu(2,ifft,idat_occ,idat)-rho12(2,ifft,nspinor,idat_occ,idat)
         end do ! ifft
         end do ! idat_occ
         end do ! idat
#endif
       end if
     end if


     call timab(1515,2,tsec) ; call timab(1513,-1,tsec) ; call timab(1544,-2,tsec)
     ! Perform an FFT using fourwf to get rhog_munu = FFT^-1(rhor_munu)
     call fourdp(cplex_fock,rhog_munu,rhor_munu,-1,mpi_enreg,nfftf,ndat*ndat_occ,&
&         ngfftf,tim_fourdp_fock_getghc,gpu_option=gpu_option)
     call timab(1513,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1544,-1,tsec)

     if(fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
       ABI_MALLOC(vfockstr, (6,ndat_occ,ndat))
       call strfock(fockcommon,gs_ham%gprimd,vfockstr,&
&                   mpi_enreg,nfftf,ngfftf,fockbz%nkpt_bz,ndat*ndat_occ,rhog_munu,gs_ham%ucvol,&
&                   qvec_j,gpu_option=gpu_option)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
         fockcommon%stress_ikpt(:,fockcommon%ieigen+idat-1)=fockcommon%stress_ikpt(:,fockcommon%ieigen+idat-1)+vfockstr(:,idat_occ,idat)*occ(idat_occ)*wtk
       end do ! idat_occ
       end do ! idat
       ABI_FREE(vfockstr)
       if (fockcommon%usepaw==0.and.(.not.need_ghc)) then
         if (allocated(fockbz%cgocc)) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET EXIT DATA MAP(delete:cwaveocc_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
           ABI_FREE(cwaveocc_r)
         end if
         call timab(1524,2,tsec) ; call timab(1544,-2,tsec)
         cycle
       end if
     end if
     call timab(1524,2,tsec) ; call timab(1544,-2,tsec)

! ===================================================
! === Calculate the local potential vfockloc_munu ===
! ===================================================
!* Apply the Poisson solver to "rhog_munu" while taking into account the effect of the vector "qvec"
!* This is precisely what is done in the subroutine hartre, with option cplex=2.
!* vfock will contain the local Fock potential, the result of hartre routine.
!* vfock = FFT( rhog_munu/|g+qvec|^2 )
     call timab(1525,-1,tsec) ; call timab(1545,-1,tsec)
#if 0

     do idat=1,ndat
     do idat_occ=1,ndat_occ
     call timab(1515,-2,tsec) ; call timab(1513,-1,tsec)
     call hartre(cplex_fock,fockcommon%gsqcut,fockcommon%usepaw,mpi_enreg,nfftf,ngfftf,&
&     mpi_enreg%paral_kgb,rhog_munu(:,:,idat_occ,idat),rprimd,vfock(:,idat_occ,idat),divgq0=fock%divgq0,qpt=qvec_j)
     call timab(1513,2,tsec) ; call timab(1515,-1,tsec)
     end do ! idat_occ
     end do ! idat

#else
     if(gpu_option==ABI_GPU_DISABLED) then
       !$OMP PARALLEL DO COLLAPSE(2)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
       do ifft=1,nfftf
         rhog_munu(1,ifft,idat_occ,idat) = rhog_munu(1,ifft,idat_occ,idat) * vqg(ifft)
         rhog_munu(2,ifft,idat_occ,idat) = rhog_munu(2,ifft,idat_occ,idat) * vqg(ifft)
       end do
       end do ! idat_occ
       end do ! idat
     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
       !$OMP& MAP(to:rhog_munu,vqg) PRIVATE(idat,idat_occ)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
       !$OMP PARALLEL DO PRIVATE(ifft)
       do ifft=1,nfftf
         rhog_munu(1,ifft,idat_occ,idat) = rhog_munu(1,ifft,idat_occ,idat) * vqg(ifft)
         rhog_munu(2,ifft,idat_occ,idat) = rhog_munu(2,ifft,idat_occ,idat) * vqg(ifft)
       end do
       end do ! idat_occ
       end do ! idat
#endif
     end if ! gpu_option

     call timab(1515,2,tsec) ; call timab(1513,-1,tsec) ; call timab(1545,-2,tsec)
     call fourdp(cplex_fock,rhog_munu,vfock,+1,mpi_enreg,nfftf,ndat*ndat_occ,ngfftf,tim_fourdp_fock_getghc,gpu_option=gpu_option)
     call timab(1513,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1545,-1,tsec)
#endif
     call timab(1525,-2,tsec) ; call timab(1545,-2,tsec)

!===============================================================
!======== Calculate Dij_Fock_hat contribution in case of PAW ===
!===============================================================

     call timab(1526,-1,tsec) ; call timab(1546,-1,tsec)

     if (fockcommon%usepaw==1) then
       qphon=qvec_j;nfftotf=product(ngfftf(1:3))
       ndij=nspden_fock
       ! dimekb1 is dimensioned as cplex_dij*lmnmax*(lmnmax+1)/2
       cplex_dij=2*gs_ham%dimekb1/(gs_ham%lmnmax*(gs_ham%lmnmax+1))
       ABI_MALLOC(dijhat,(gs_ham%dimekb1,natom,ndij,ndat_occ,cplex_fock,ndat))
       dijhat=zero

#ifdef HAVE_OPENMP_OFFLOAD
       !!$OMP TARGET UPDATE FROM(vfock) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       iatm=0
       do itypat=1,gs_ham%ntypat
         lmn2_size=fockcommon%pawtab(itypat)%lmn2_size
         ABI_MALLOC(dijhat_tmp,(cplex_fock*cplex_dij*lmn2_size,ndij*ndat_occ*ndat,gs_ham%nattyp(itypat)))
         call pawdijhat_ndat(dijhat_tmp,cplex_dij,cplex_fock,gs_ham%gprimd,iatm,&
&           natom,ndij,nfftf,nfftotf,nspden_fock,nspden_fock,ndat_occ*ndat,&
&           gs_ham%nattyp(itypat),fockbz%pawang,fockcommon%pawfgrtab,&
&           fockcommon%pawtab(itypat),vfock,qphon,gs_ham%ucvol,gs_ham%xred,&
&           gpu_option=gpu_option)
         do ia=1,gs_ham%nattyp(itypat)
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               do ii=1,cplex_fock
                 iatom=iatm+ia
                 ind=(ii-1)*lmn2_size*cplex_dij
                 dijhat(1:cplex_dij*lmn2_size,iatom,:,idat_occ,ii,idat)=&
  &                dijhat_tmp(ind+1:ind+cplex_dij*lmn2_size,&
                              1+(idat_occ-1)*ndij+(idat-1)*ndat_occ*ndij:idat_occ*ndij+(idat-1)*ndat_occ*ndij,ia)
               end do
             end do ! idat_occ
           end do ! idat
         end do ! ia
         ABI_FREE(dijhat_tmp)
         iatm=iatm+gs_ham%nattyp(itypat)
       end do
       signs=2; cpopt=2;idir=0; paw_opt=1;nnlout=1;tim_nonlop=17

       if(need_ghc) then
         choice=1
         call timab(1515,2,tsec) ; call timab(1514,-1,tsec) ; call timab(1546,-2,tsec)
         do idat=1,ndat
           call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),&
&               mpi_enreg,ndat_occ,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
&               gvnlxc,enl_ndat=dijhat(:,:,:,:,:,idat),&
&               select_k=K_H_KPRIME)

           if(gpu_option==ABI_GPU_DISABLED) then
             do idat_occ=1,ndat_occ
               ghc2(:,1+(idat-1)*npw*nspinor:idat*npw*nspinor)=ghc2(:,1+(idat-1)*npw*nspinor:idat*npw*nspinor)&
  &               -gvnlxc(:,1+(idat_occ-1)*npw*nspinor:idat_occ*npw*nspinor)*occ(idat_occ)*wtk
             end do ! idat_occ
           else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
             !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO &
             !$OMP& MAP(to:gvnlxc,ghc2,occ) PRIVATE(idat_occ,ipw)
             do ipw=1,npw
               do idat_occ=1,ndat_occ
                 ghc2(1:2,ipw+(idat-1)*npw*nspinor)=ghc2(1:2,ipw+(idat-1)*npw*nspinor)&
    &               -gvnlxc(1:2,ipw+(idat_occ-1)*npw*nspinor)*occ(idat_occ)*wtk
               end do
             end do ! idat_occ
#endif
           end if
         end do ! idat
         call timab(1514,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1546,-1,tsec)
       end if

! Forces calculation

       if (fockcommon%optfor.and.(fockcommon%ieigen/=0)) then
         ABI_MALLOC(vdotr,(ndat_occ,3,natom,ndat))
         ABI_MALLOC(vdoti,(ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET ENTER DATA MAP(alloc:vdotr,vdoti) IF(gpu_option==ABI_GPU_OPENMP)
#endif
         ABI_MALLOC(for1,(ndat_occ,3,natom,ndat))
         ABI_MALLOC(atom_nfgd,    (natom))
         do iatom=1,natom
           atom_nfgd(iatom) =      fockcommon%pawfgrtab(iatom)%nfgd
         end do
         ABI_MALLOC(atom_ifftsph, (maxval(atom_nfgd), natom))
         ABI_MALLOC(atom_rfgd,    (3, maxval(atom_nfgd), natom))
         do iatom=1,natom
           atom_ifftsph(1:atom_nfgd(iatom),iatom) = fockcommon%pawfgrtab(iatom)%ifftsph(1:atom_nfgd(iatom))
           atom_rfgd(:,1:atom_nfgd(iatom),iatom) =  fockcommon%pawfgrtab(iatom)%rfgd(:,1:atom_nfgd(iatom))
         end do
         choice=2; vdotr=zero;doti=zero;cpopt=4;tim_nonlop=17
         do idir=1,3
           do idat=1,ndat
             do iatom=1,natom
               call timab(1515,2,tsec) ; call timab(1514,-1,tsec) ; call timab(1546,-2,tsec)
               call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
  &             ndat_occ,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
  &             forout,enl_ndat=dijhat(:,:,:,:,:,idat),iatom_only=iatom,&
  &             select_k=K_H_KPRIME)
               call timab(1514,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1546,-1,tsec)
               call dotprod_g_batch_half(vdotr(:,idir,iatom,idat),vdoti,gs_ham%istwf_k,npw,ndat_occ,2,&
                 cwavef(:,npw*nspinor*(idat-1)+1:npw*nspinor*idat),&
                 forout,mpi_enreg%me_g0,mpi_enreg%comm_fft,gpu_option=gpu_option)
             end do ! iatom
           end do ! idat
         end do ! idir

         if(gpu_option==ABI_GPU_DISABLED) then
           do idat=1,ndat
             do iatom=1,natom
               do idir=1,3
                 do idat_occ=1,ndat_occ
                   esum=0
                   do ifft=1,atom_nfgd(iatom)
                     ind=atom_ifftsph(ifft,iatom)
                     esum=esum &
    &                + vfock(2*ind-1,idat_occ,idat)*grnhat_12(1,ind,1,idir,iatom,idat_occ,idat) &
    &                - vfock(2*ind,idat_occ,idat)*grnhat_12(2,ind,1,idir,iatom,idat_occ,idat)
                   end do
                   for1(idat_occ,idir,iatom,idat)=esum
                 end do ! idat_occ
               end do ! idir
             end do ! iatom
           end do ! idat
         else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) MAP(tofrom:for1) &
           !$OMP& MAP(to:vfock,grnhat_12,atom_nfgd,atom_rfgd,atom_ifftsph) &
           !$OMP& PRIVATE(ifft,ind,iatom) PRIVATE(idat_occ,idir,esum)
           do idat=1,ndat
             do iatom=1,natom
               !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(esum,ind,idir,idat_occ,ifft)
               do idir=1,3
                 do idat_occ=1,ndat_occ
                   esum=0
                   do ifft=1,atom_nfgd(iatom)
                     ind=atom_ifftsph(ifft,iatom)
                     esum=esum &
    &                + vfock(2*ind-1,idat_occ,idat)*grnhat_12(1,ind,1,idir,iatom,idat_occ,idat) &
    &                - vfock(2*ind,idat_occ,idat)*grnhat_12(2,ind,1,idir,iatom,idat_occ,idat)
                   end do
                   for1(idat_occ,idir,iatom,idat)=esum
                 end do ! idat_occ
               end do ! idir
             end do ! iatom
           end do ! idat
#endif
         end if

         if(.true.) then
         !if(gpu_option==ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET UPDATE FROM(vdotr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
           do idat=1,ndat
             do iatom=1,natom
               do idir=1,3
                 do idat_occ=1,ndat_occ
                   for12(idir)=rprimd(1,idir)*for1(idat_occ,1,iatom,idat)&
                   &          +rprimd(2,idir)*for1(idat_occ,2,iatom,idat)&
                   &          +rprimd(3,idir)*for1(idat_occ,3,iatom,idat)
                   forikpt(idir,iatom,idat)=forikpt(idir,iatom,idat)&
                       -(for12(idir)*gs_ham%ucvol/nfftf+vdotr(idat_occ,idir,iatom,idat))*occ(idat_occ)*wtk
                 end do ! idat_occ
               end do ! idir
             end do ! iatom
           end do ! idat
         end if
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET EXIT DATA MAP(delete:vdotr,vdoti) IF(gpu_option==ABI_GPU_OPENMP)
#endif
         ABI_FREE(vdotr)
         ABI_FREE(vdoti)
         ABI_FREE(for1)
         ABI_FREE(atom_ifftsph)
         ABI_FREE(atom_nfgd)
         ABI_FREE(atom_rfgd)
       end if

! Stresses calculation
       if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
         signs=2;choice=3;cpopt=4;tim_nonlop=17

       ! first contribution
         ABI_MALLOC(vdotr,(ndat_occ,6,1,1))
         ABI_MALLOC(vdoti,(ndat_occ))
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET ENTER DATA MAP(alloc:vdotr,vdoti) IF(gpu_option==ABI_GPU_OPENMP)
#endif
         do idir=1,6
           do idat=1,ndat
             call timab(1515,2,tsec) ; call timab(1514,-1,tsec) ; call timab(1546,-2,tsec)
             call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
  &           ndat_occ,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
  &           strout,enl_ndat=dijhat(:,:,:,:,:,idat),select_k=K_H_KPRIME)
             call timab(1514,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1546,-1,tsec)
             call dotprod_g_batch_half(vdotr(:,idir,1,1),vdoti,gs_ham%istwf_k,npw,ndat_occ,2,&
                   cwavef(:,npw*nspinor*(idat-1)+1:npw*nspinor*idat),strout,&
                   mpi_enreg%me_g0,mpi_enreg%comm_fft,gpu_option=gpu_option)
#ifdef HAVE_OPENMP_OFFLOAD
             !$OMP TARGET UPDATE FROM(vdotr) IF(gpu_option==ABI_GPU_OPENMP)
#endif
             do idat_occ=1,ndat_occ
               fockcommon%stress_ikpt(idir,fockcommon%ieigen+idat-1)=fockcommon%stress_ikpt(idir,fockcommon%ieigen+idat-1)-&
    &             vdotr(idat_occ,idir,1,1)*occ(idat_occ)*wtk/gs_ham%ucvol
             end do ! idat_occ
           end do ! idat
         end do ! idir
#ifdef HAVE_OPENMP_OFFLOAD
         !$OMP TARGET EXIT DATA MAP(delete:vdotr,vdoti) IF(gpu_option==ABI_GPU_OPENMP)
#endif
         ABI_FREE(vdotr)
         ABI_FREE(vdoti)

         ABI_MALLOC(atom_nfgd,    (natom))
         do iatom=1,natom
           atom_nfgd(iatom) =      fockcommon%pawfgrtab(iatom)%nfgd
         end do
         ABI_MALLOC(atom_ifftsph, (maxval(atom_nfgd), natom))
         ABI_MALLOC(atom_rfgd,    (3, maxval(atom_nfgd), natom))
         do iatom=1,natom
           atom_ifftsph(1:atom_nfgd(iatom),iatom) = fockcommon%pawfgrtab(iatom)%ifftsph(1:atom_nfgd(iatom))
           atom_rfgd(:,1:atom_nfgd(iatom),iatom) =  fockcommon%pawfgrtab(iatom)%rfgd(:,1:atom_nfgd(iatom))
         end do
         stress_ikpt =>  fockcommon%stress_ikpt
         ieigen = fockcommon%ieigen

       ! second contribution
         !if(.true.) then
         if(gpu_option==ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
           !$OMP TARGET UPDATE FROM(vfock,grnhat_12) IF(gpu_option==ABI_GPU_OPENMP)
#endif
           ABI_MALLOC(strdat, (3,3,ndat_occ,ndat))
           strdat=zero
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               do idir=1,3
                 do idir1=1,3
                   esum=0
                   do iatom=1,natom
                     do ifft=1,atom_nfgd(iatom)
                       !ind=fockcommon%pawfgrtab(iatom)%ifftsph(ifft)
                       ind=atom_ifftsph(ifft,iatom)
                       !strdat(idir,idir1,idat_occ,idat)=strdat(idir,idir1,idat_occ,idat)+(vfock(2*ind-1,idat_occ,idat)*grnhat_12(1,ind,1,idir,iatom,idat_occ,idat)-&
                       esum=esum+(vfock(2*ind-1,idat_occ,idat)*grnhat_12(1,ind,1,idir,iatom,idat_occ,idat)-&
                          vfock(2*ind,idat_occ,idat)*grnhat_12(2,ind,1,idir,iatom,idat_occ,idat))*&
                          atom_rfgd(idir1,ifft,iatom)
                          !fockcommon%pawfgrtab(iatom)%rfgd(idir1,ifft)
                      !   vfock(2*ind,idat_occ,idat)*grnhat_12(2,ind,1,idir,iatom,idat_occ,idat))*atom_rfgd(idir1,ifft,iatom)
                     end do
                   end do
                   strdat(idir,idir1,idat_occ,idat)=esum
                 end do
               end do
             end do ! idat_occ
           end do ! idat
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               do idir=1,3
                 fockstr(idir)=strdat(idir,idir,idat_occ,idat)
               end do
               fockstr(4)=(strdat(3,2,idat_occ,idat)+strdat(2,3,idat_occ,idat))*half
               fockstr(5)=(strdat(3,1,idat_occ,idat)+strdat(1,3,idat_occ,idat))*half
               fockstr(6)=(strdat(1,2,idat_occ,idat)+strdat(2,1,idat_occ,idat))*half
               do idir=1,6
                 fockcommon%stress_ikpt(idir,fockcommon%ieigen+idat-1)=fockcommon%stress_ikpt(idir,fockcommon%ieigen+idat-1)+&
        &           fockstr(idir)/nfftf*occ(idat_occ)*wtk
               end do
             end do ! idat_occ
           end do ! idat
           ABI_FREE(strdat)
         else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           ABI_MALLOC(strdat, (3,3,ndat_occ,ndat))
           strdat=zero
           !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(4) MAP(tofrom:strdat) &
           !$OMP& MAP(to:vfock,grnhat_12,atom_nfgd,atom_rfgd,atom_ifftsph) &
           !$OMP& PRIVATE(idat,idat_occ,idir,idir1,esum,iatom)
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               do idir=1,3
                 do idir1=1,3
                   do iatom=1,natom
                     esum=0
                     !$OMP PARALLEL DO PRIVATE(ifft,ind) REDUCTION(+:esum)
                     do ifft=1,atom_nfgd(iatom)
                       ind=atom_ifftsph(ifft,iatom)
                       esum=esum+(vfock(2*ind-1,idat_occ,idat)*grnhat_12(1,ind,1,idir,iatom,idat_occ,idat)-&
      &                   vfock(2*ind,idat_occ,idat)*grnhat_12(2,ind,1,idir,iatom,idat_occ,idat))*atom_rfgd(idir1,ifft,iatom)
                     end do
                     strdat(idir,idir1,idat_occ,idat)=strdat(idir,idir1,idat_occ,idat)+esum
                   end do
                 end do
               end do
             end do ! idat_occ
           end do ! idat
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               do idir=1,3
                 fockstr(idir)=strdat(idir,idir,idat_occ,idat)
               end do
               fockstr(4)=(strdat(3,2,idat_occ,idat)+strdat(2,3,idat_occ,idat))*half
               fockstr(5)=(strdat(3,1,idat_occ,idat)+strdat(1,3,idat_occ,idat))*half
               fockstr(6)=(strdat(1,2,idat_occ,idat)+strdat(2,1,idat_occ,idat))*half
               do idir=1,6
                 stress_ikpt(idir,ieigen+idat-1)=stress_ikpt(idir,ieigen+idat-1)+&
        &           fockstr(idir)/nfftf*occ(idat_occ)*wtk
               end do
             end do ! idat_occ
           end do ! idat
           ABI_FREE(strdat)
#endif
         end if
         ABI_FREE(atom_ifftsph)
         ABI_FREE(atom_nfgd)
         ABI_FREE(atom_rfgd)

       ! third contribution
         if(gpu_option==ABI_GPU_DISABLED) then
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               doti=zero
               do ifft=1,nfftf
                 doti=doti+vfock(2*ifft-1,idat_occ,idat)*rho12(1,ifft,nspinor,idat_occ,idat)-vfock(2*ifft,idat_occ,idat)*rho12(2,ifft,nspinor,idat_occ,idat)
               end do
               fockcommon%stress_ikpt(1:3,fockcommon%ieigen+idat-1)=fockcommon%stress_ikpt(1:3,fockcommon%ieigen+idat-1)-doti/nfftf*occ(idat_occ)*wtk
             end do ! idat_occ
           end do ! idat
         else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
           ABI_MALLOC(strdat, (1,1,ndat_occ,ndat))
           strdat=zero
           !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(2) &
           !$OMP& MAP(tofrom:strdat) MAP(to:vfock,rho12) &
           !$OMP& PRIVATE(idat,idat_occ,idir,idir1,doti)
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               doti=zero
               !$OMP PARALLEL DO PRIVATE(ifft,ind,iatom) REDUCTION(+:doti)
               do ifft=1,nfftf
                 doti=doti+vfock(2*ifft-1,idat_occ,idat)*rho12(1,ifft,nspinor,idat_occ,idat)&
                          -vfock(2*ifft,idat_occ,idat)  *rho12(2,ifft,nspinor,idat_occ,idat)
               end do
               strdat(1,1,idat_occ,idat)=doti
             end do ! idat_occ
           end do ! idat
           do idat=1,ndat
             do idat_occ=1,ndat_occ
               fockcommon%stress_ikpt(1:3,fockcommon%ieigen+idat-1)=fockcommon%stress_ikpt(1:3,fockcommon%ieigen+idat-1) &
                 -strdat(1,1,idat_occ,idat)/nfftf*occ(idat_occ)*wtk
             end do ! idat_occ
           end do ! idat
           ABI_FREE(strdat)
#endif
         end if
       end if ! end stresses

       ABI_FREE(dijhat)
     end if !end PAW
     call timab(1526,2,tsec) ; call timab(1546,-2,tsec)

! =============================================================
! === Apply the local potential vfockloc_munu to cwaveocc_r ===
! =============================================================
     call timab(1527,-1,tsec)
     if(gpu_option==ABI_GPU_DISABLED) then
       !$OMP PARALLEL DO &
       !$OMP& PRIVATE(ind,recwocc,imcwocc,revloc,imvloc)
       do idat=1,ndat
       do idat_occ=1,ndat_occ
       do i3=1,ngfftf(3)
         do i2=1,ngfftf(2)
           do i1=1,ngfftf(1)
             ind=i1+(i2-1)*ngfftf(1)+(i3-1)*ngfftf(2)*ngfftf(1)
             revloc=vfock(2*ind-1,idat_occ,idat) ; imvloc=vfock(2*ind,idat_occ,idat)
             recwocc=cwaveocc_r(1,i1,i2,i3,idat_occ)
             imcwocc=cwaveocc_r(2,i1,i2,i3,idat_occ)
             vlocpsi_r(2*ind-1,idat)=vlocpsi_r(2*ind-1,idat)-(revloc*recwocc-imvloc*imcwocc)*occ(idat_occ)*wtk
             vlocpsi_r(2*ind  ,idat)=vlocpsi_r(2*ind  ,idat)-(revloc*imcwocc+imvloc*recwocc)*occ(idat_occ)*wtk
           end do
         end do
       end do
       end do ! idat_occ
       end do ! idat
     else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
       do idat_occ=1,ndat_occ
         !$OMP TARGET TEAMS DISTRIBUTE &
         !$OMP& MAP(to:vlocpsi_r) MAP(to:cwaveocc_r,occ,vfock) PRIVATE(idat)
         do idat=1,ndat
           !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ind,recwocc,imcwocc,revloc,imvloc,i3,i2,i1)
           do i3=1,ngfftf(3)
             do i2=1,ngfftf(2)
               do i1=1,ngfftf(1)
                 ind=i1+(i2-1)*ngfftf(1)+(i3-1)*ngfftf(2)*ngfftf(1)
                 revloc=vfock(2*ind-1,idat_occ,idat) ; imvloc=vfock(2*ind,idat_occ,idat)
                 recwocc=cwaveocc_r(1,i1,i2,i3,idat_occ)
                 imcwocc=cwaveocc_r(2,i1,i2,i3,idat_occ)
                 vlocpsi_r(2*ind-1,idat)=vlocpsi_r(2*ind-1,idat)-(revloc*recwocc-imvloc*imcwocc)*occ(idat_occ)*wtk
                 vlocpsi_r(2*ind  ,idat)=vlocpsi_r(2*ind  ,idat)-(revloc*imcwocc+imvloc*recwocc)*occ(idat_occ)*wtk
               end do
             end do
           end do
         end do ! idat
       end do ! idat_occ
#endif
     end if
     if (allocated(fockbz%cgocc)) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET EXIT DATA MAP(delete:cwaveocc_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       ABI_FREE(cwaveocc_r)
     end if
     call timab(1527,2,tsec)
   end do ! jband

! ========================================================
! === End of loop : update of shifts and deallocations ===
! ==============================:=========================
!* Update of the shifts to be applied (reminder : mkmem is not 0, nspinor=1)
   call timab(1528,1,tsec)
   jcg=jcg+npwj*nband_k
   jbg=jbg+nband_k
   bdtot_jindex=bdtot_jindex+nband_k
   jkg=jkg+npwj
   if (fockcommon%usepaw==1) then
     ABI_FREE(enlout_dum)
     ABI_FREE(vectin_dum)
     ABI_FREE(ffnl_kp_dum)
     ABI_FREE(kpg_kp)
   end if
   if (associated(gs_ham%ph3d_kp)) then
     ABI_FREE(gs_ham%ph3d_kp)
   end if
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:rhor_munu,rhog_munu,vfock,occ,rho12,gvnlxc) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(rhor_munu)
   ABI_FREE(rhog_munu)
   ABI_FREE(vfock)
   ABI_FREE(occ)
  !*Additional arrays in case of paw
   if (fockcommon%usepaw==1) then
     if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET EXIT DATA MAP(delete:forout) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       ABI_FREE(forout)
     end if
     if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
#ifdef HAVE_OPENMP_OFFLOAD
       !$OMP TARGET EXIT DATA MAP(delete:strout) IF(gpu_option==ABI_GPU_OPENMP)
#endif
       ABI_FREE(strout)
     end if
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET EXIT DATA MAP(delete:grnhat_12) IF(gpu_option==ABI_GPU_OPENMP .and. ider==3)
#endif
     ABI_FREE(grnhat_12)
     ABI_FREE(gvnlxc)
     ABI_FREE(grnhat12)
     ABI_FREE(rho12)
   end if

   call timab(1528,2,tsec)

 end do ! jkpt

! ========================================================
! === After loop                                       ===
! ========================================================

 call timab(1506,2,tsec) ; call timab(1507,1,tsec) ; call timab(1547,-1,tsec)

 if (fockcommon%usepaw==1) then
   if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
     call timab(1547,2,tsec) ; call timab(1548,-1,tsec)
     call xmpi_sum(forikpt,mpi_enreg%comm_hf,ier)
     call timab(1548,2,tsec) ; call timab(1547,-1,tsec)
     do idat=1,ndat
       do iatom=1,natom !Loop over atom
         ia=gs_ham%atindx(iatom)
         fockcommon%forces_ikpt(:,ia,fockcommon%ieigen+idat-1)=forikpt(:,iatom,idat)
       end do
     end do
   end if
 end if
 if(fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
   call timab(1547,2,tsec) ; call timab(1548,-1,tsec)
   call xmpi_sum(fockcommon%stress_ikpt,mpi_enreg%comm_hf,ier)
   call timab(1548,2,tsec) ; call timab(1547,-1,tsec)
 end if

 if (.not.need_ghc) then

! ===============================
! === Deallocate local arrays ===
! ===============================
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:cwavef_r,vqg,vlocpsi_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(cwavef_r)
   ABI_FREE(vlocpsi_r)
   ABI_FREE(dummytab)
   ABI_FREE(vqg)
   if(fockcommon%usepaw==1.or.fockcommon%optstr) then
     ABI_FREE(gboundf)
   end if
   if (fockcommon%usepaw==1) then
     if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
       ABI_FREE(forikpt)
     end if
   end if
!*Restore gs_ham datastructure

   if (associated(gs_ham%ph3d_kp)) then
     ABI_MALLOC(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
   end if
   call gs_ham%load_kprime(kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
&   npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!   if (fockcommon%ieigen/=0) fockcommon%ieigen=0

 else

!  *Restore gs_ham datastructure

   if (associated(gs_ham%ph3d_kp)) then
     ABI_MALLOC(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
   end if
   call gs_ham%load_kprime(kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
&   npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!  * Perform an FFT using fourwf to get ghc1 = FFT^-1(vlocpsi_r)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET UPDATE FROM(vlocpsi_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_MALLOC(psilocal,(cplex_fock*n4f,n5f,n6f*ndat))
   do idat=1,ndat
     call fftpac(1,mpi_enreg,nspden_fock,cplex_fock*n1f,n2f,n3f,&
       cplex_fock*n4f,n5f,n6f,ngfft,vlocpsi_r(:,idat),psilocal(:,:,1+(idat-1)*n6f:idat*n6f),2)
   end do ! idat

   call timab(1515,2,tsec) ; call timab(1512,-1,tsec) ; call timab(1547,-2,tsec)
   call fourwf(0,rhodum0,rhodum,ghc1,psilocal,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
&   mgfftf,mpi_enreg,ndat,ngfftf,1,npw,n4f,n5f,n6f,3,tim_fourwf_fock_getghc,weight1,weight1,&
&   gpu_option=gs_ham%gpu_option)
   call timab(1512,2,tsec) ; call timab(1515,-1,tsec) ; call timab(1547,-1,tsec)
   ABI_FREE(psilocal)

   if(gpu_option==ABI_GPU_DISABLED) then
     ghc1=ghc1*sqrt(gs_ham%ucvol)+ghc2
   else if(gpu_option==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(ghc1,ghc2)
     call abi_gpu_xaxpy(2,npw*ndat,cucvol,c_loc(ghc1),1,c_loc(ghc2),1)
     !$OMP END TARGET DATA
     call gpu_copy(ghc1,ghc2,int(2,c_size_t)*npw*ndat)
     !$OMP TARGET UPDATE FROM(ghc1)
#endif
   end if

!  * If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
   ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_hf + ghc1(:,:)

   call timab(1547,2,tsec) ; call timab(1548,-1,tsec)
   call xmpi_sum(ghc,mpi_enreg%comm_hf,ier)
   call timab(1548,2,tsec) ; call timab(1547,-1,tsec)

!   ===============================
!   === Deallocate local PAW arrays ===
!   ===============================
   if (fockcommon%usepaw==1) then
     if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
       ABI_FREE(forikpt)
     end if
   end if
   if(fockcommon%usepaw==1.or.fockcommon%optstr) then
     ABI_FREE(gboundf)
   end if
!   ============================================
!   === Calculate the contribution to energy ===
!   ============================================
!  * Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!  * at each convergence step. (consistent definition with the definition of hartree energy)
   if (fockcommon%ieigen/=0) then
     do idat=1,ndat
       eigen=zero
  !  * Dot product of cwavef and ghc
  !  * inspired from the routine 54_spacepar/meanvalue_g but without the reference to parallelism and filtering
       if(gs_ham%istwf_k==2) then
         eigen=half*cwavef(1,1+(idat-1)*npw)*ghc1(1,1+(idat-1)*npw)
       else
         eigen=cwavef(1,1+(idat-1)*npw)*ghc1(1,1+(idat-1)*npw)+cwavef(2,1+(idat-1)*npw)*ghc1(2,1+(idat-1)*npw)
       end if
       do ipw=2,npw
         eigen=eigen+cwavef(1,ipw+(idat-1)*npw)*ghc1(1,ipw+(idat-1)*npw)+cwavef(2,ipw+(idat-1)*npw)*ghc1(2,ipw+(idat-1)*npw)
       end do
       if(gs_ham%istwf_k>=2) eigen=two*eigen
       call timab(1547,2,tsec) ; call timab(1548,-1,tsec)
       call xmpi_sum(eigen,mpi_enreg%comm_hf,ier)
       call timab(1548,2,tsec) ; call timab(1547,-1,tsec)
       fockcommon%eigen_ikpt(fockcommon%ieigen+idat-1)= eigen
       if(fockcommon%use_ACE==0) fockcommon%ieigen = 0
     end do ! idat
   end if

!   ===============================
!   === Deallocate local arrays ===
!   ===============================
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:cwavef_r,ghc1,ghc2,vqg,vlocpsi_r) IF(gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(cwavef_r)
   ABI_FREE(ghc1)
   ABI_FREE(ghc2)
   ABI_FREE(vlocpsi_r)
   ABI_FREE(dummytab)
   ABI_FREE(vqg)

 endif

 call timab(1504,2,tsec) ; call timab(1507,-2,tsec) ; call timab(1515,-2,tsec) ; call timab(1547,-2,tsec)

 end subroutine fock_getghc
!!***

!!****f* ABINIT/fock2ACE
!! NAME
!! fock2ACE
!!
!! FUNCTION
!! Compute nonlocal contribution to the Fock part of the hamiltonian in the ACE formalism.
!! optionally contribution to Fock forces
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions (may be read from disk file)
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced coordinates (integers) of G vecs in basis
!!  kpt(3,nkpt)=k points in reduced coordinates
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpsang=
!!  mpw= maximum number of plane waves
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nband(nkpt)=number of bands at each k point
!!  nfft=number of FFT grid points
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points in Brillouin zone
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves in basis and boundary at each k
!!  nspden=Number of spin Density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms
!!  occ(mband*nkpt*nsppol)=occupation numbers for each band over all k points
!!  optfor=1 if computation of forces is required
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  usecprj=1 if cprj datastructure has been allocated
!!  gpu_option= GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!  wtk(nkpt)=weight associated with each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!
!! fock%fockACE(ikpt,isppol)%xi
!! if optfor=1, fock%fock_common%forces
!!
!! SOURCE

subroutine fock2ACE(cg,cprj,fock,istwfk,kg,kpt,mband,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,my_natom,natom,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,&
&  ntypat,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,typat,usecprj,gpu_option,wtk,xred,ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfft,mkmem,mpsang,mpw,my_natom,natom,nfft,nkpt
 integer,intent(in) :: nspden,nsppol,nspinor,ntypat,optfor
 integer,intent(in) :: usecprj,gpu_option
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(3),npwarr(nkpt)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
 type(fock_type),pointer, intent(inout) :: fock
!Local variables-------------------------------
!scalars
 integer :: bandpp,bdtot_index,dimffnl,iband,iband_cprj,iband_last,ibg,icg,ider
 integer :: ierr,info,idir,ikg,ikpt,ilm,ipw,isppol,istwf_k,kk,ll
 integer :: mband_cprj,me_distrb,my_ikpt,my_nspinor,nband_k,nband_cprj_k,nkpg
 integer :: npw_k,spaceComm
 integer :: use_ACE_old
 integer :: blocksize,iblock,jblock,iblocksize,jblocksize,nblockbd
!integer, save :: counter=0
 type(gs_hamiltonian_type) :: gs_hamk
 logical :: compute_gbound
 character(len=500) :: msg
 type(fock_common_type),pointer :: fockcommon
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: bb(:,:,:),cwavef(:,:),cwavefk(:,:),ffnl_sav(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:),kpg_k_sav(:,:)
 real(dp),allocatable :: mkl(:,:,:),occblock(:),ph3d(:,:,:),ph3d_sav(:,:,:)
 real(dp),allocatable :: wi(:,:,:),weight(:),ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),allocatable,target :: ffnl(:,:,:,:)
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 type(pawcprj_type),target,allocatable :: cwaveprj(:,:)

!*************************************************************************

 call timab(1560,1,tsec)
 call timab(1561,1,tsec)

!DEBUG
!if(counter>0)return
!counter=counter+1
!ENDDEBUG

!Init mpicomm and me
 if(mpi_enreg%paral_kgb==1)then
   spaceComm=mpi_enreg%comm_kpt
   me_distrb=mpi_enreg%me_kpt
 else
!* In case of HF calculation
   if (mpi_enreg%paral_hf==1) then
     spaceComm=mpi_enreg%comm_kpt
     me_distrb=mpi_enreg%me_kpt
   else
     spaceComm=mpi_enreg%comm_cell
     me_distrb=mpi_enreg%me_cell
   end if
 end if

!Some initializations
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
 compute_gbound=.true.
 fockcommon => fock%fock_common
 use_ACE_old=fockcommon%use_ACE
 fockcommon%use_ACE=0
 fockcommon%e_fock0=zero

!Initialize Hamiltonian (k- and spin-independent terms)

 call gs_hamk%init(psps,pawtab,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,usecprj=usecprj,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& paw_ij=paw_ij,ph1d=ph1d,fock=fock,&
& gpu_option=gpu_option)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)
 fockcommon%use_ACE=use_ACE_old

!need to reorder cprj=<p_lmn|Cnk> (from unsorted to atom-sorted)
 if (psps%usepaw==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx)
 end if

!LOOP OVER SPINS
 bdtot_index=0;ibg=0;icg=0
 call timab(1561,2,tsec) ; call timab(1562,-1,tsec)

 do isppol=1,nsppol
   fockcommon%isppol=isppol
!  Continue to initialize the Hamiltonian (PAW DIJ coefficients)
   call gs_hamk%load_spin(isppol,with_nonlocal=.true.)

!  Loop over k points
   ikg=0
   do ikpt=1,nkpt
     fockcommon%ikpt=ikpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     npw_k=npwarr(ikpt)
     kpoint(:)=kpt(:,ikpt)
     istwf_k=istwfk(ikpt)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

!    Parallelism over FFT and/or bands: define sizes and tabs
     if (mpi_enreg%paral_kgb==1) then
       my_ikpt=mpi_enreg%my_kpttab(ikpt)
       nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
       bandpp=mpi_enreg%bandpp
       my_bandfft_kpt => bandfft_kpt(my_ikpt)
     else
       my_ikpt=ikpt
       bandpp=mpi_enreg%bandpp
       nblockbd=nband_k/bandpp
     end if
     blocksize=nband_k/nblockbd
     mband_cprj=mband/mpi_enreg%nproc_band
     nband_cprj_k=nband_k/mpi_enreg%nproc_band

     ABI_MALLOC(cwavef,(2,npw_k*my_nspinor*blocksize))
     if (psps%usepaw==1) then
       ABI_MALLOC(cwaveprj,(natom,my_nspinor*bandpp))
       call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
     else
       ABI_MALLOC(cwaveprj,(0,0))
     end if

     ABI_MALLOC(kg_k,(3,mpw))
!$OMP PARALLEL DO
     do ipw=1,npw_k
       kg_k(:,ipw)=kg(:,ipw+ikg)
     end do

     ABI_MALLOC(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_MALLOC(ylmgr_k,(0,0,0))
     if (psps%useylm==1) then
!$OMP PARALLEL DO COLLAPSE(2)
       do ilm=1,mpsang*mpsang
         do ipw=1,npw_k
           ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
         end do
       end do
     end if

!    Compute (k+G) vectors
     nkpg=3*nloalg(3)
     ABI_MALLOC(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if


!    Compute nonlocal form factors ffnl at all (k+G)
     ider=0;idir=0;dimffnl=1
     ABI_MALLOC(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&     ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&     nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!    Load k-dependent part in the Hamiltonian datastructure
!     - Compute 3D phase factors
!     - Prepare various tabs in case of band-FFT parallelism
!     - Load k-dependent quantities in the Hamiltonian

     ABI_MALLOC(ph3d,(2,npw_k,gs_hamk%matblk))
     call gs_hamk%load_k(kpt_k=kpoint,istwf_k=istwf_k,npw_k=npw_k,&
&     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_gbound=compute_gbound,compute_ph3d=.true.)

!    Load band-FFT tabs (transposed k-dependent arrays)
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_savetabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kpg=kpg_k_sav)
       call prep_bandfft_tabs(gs_hamk,ikpt,mkmem,mpi_enreg)
       call gs_hamk%load_k(npw_fft_k=my_bandfft_kpt%ndatarecv, &
&       kg_k     =my_bandfft_kpt%kg_k_gather, &
&       kpg_k    =my_bandfft_kpt%kpg_k_gather, &
       ffnl_k   =my_bandfft_kpt%ffnl_gather, &
       ph3d_k   =my_bandfft_kpt%ph3d_gather,compute_gbound=compute_gbound)
     end if

!    The following is now wrong. In sequential, nblockbd=nband_k/bandpp
!    blocksize= bandpp (JB 2016/04/16)
!    Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
!
     ABI_MALLOC(occblock,(blocksize))
     ABI_MALLOC(weight,(blocksize))
     occblock=zero;weight=zero

     if (fockcommon%optfor) then
       fockcommon%forces_ikpt=zero
     end if

     ABI_MALLOC(wi,(2,npw_k*my_nspinor*blocksize,nblockbd))
     wi=zero
     ABI_MALLOC(mkl,(2,nband_k,nband_k))
     mkl=zero
! Calculate all the Wi for the current k-point

     do iblock=1,nblockbd

       iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband_k)
       iband_cprj=(iblock-1)*bandpp+1
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,me_distrb)) cycle

!      Select occupied bandsddk
       occblock(:)=occ(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
       weight(:)=wtk(ikpt)*occblock(:)

!        Load contribution from n,k
       cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&       cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
       if (psps%usepaw==1) then
         call pawcprj_get(gs_hamk%atindx1,cwaveprj,cprj,natom,iband_cprj,ibg,ikpt,0,isppol,&
&         mband_cprj,mkmem,natom,bandpp,nband_cprj_k,my_nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if

       if (mpi_enreg%paral_kgb==1) then
         msg='fock2ACE: Paral_kgb is not yet implemented for fock calculations'
         ABI_BUG(msg)
       end if

       fockcommon%ieigen=(iblock-1)*blocksize+1
       fockcommon%iband=(iblock-1)*blocksize+1
       call timab(1562,2,tsec) ; call timab(1563,-1,tsec)
       call fock_getghc(cwavef,&
             cwaveprj,&
             wi(:,:,iblock),&
             gs_hamk,mpi_enreg,blocksize)


       call timab(1563,2,tsec) ; call timab(1562,-1,tsec)
       do iblocksize=1,blocksize
         mkl(1,fockcommon%ieigen+iblocksize-1,fockcommon%ieigen+iblocksize-1)=fockcommon%eigen_ikpt(fockcommon%ieigen+iblocksize-1)
         fockcommon%e_fock0=fockcommon%e_fock0+half*weight(iblocksize)*fockcommon%eigen_ikpt(fockcommon%ieigen+iblocksize-1)
         if (fockcommon%optfor) then
           fockcommon%forces(:,:)=fockcommon%forces(:,:)+weight(iblocksize)*fockcommon%forces_ikpt(:,:,fockcommon%ieigen+iblocksize-1)
         end if
       end do


     end do ! End of loop on block of bands

! Calculate Mkl for the current k-point
     ABI_MALLOC(cwavefk,(2,npw_k*my_nspinor))
     do iblock=1,nblockbd
       cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&       cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
       do iblocksize=1,blocksize
         kk=(iblock-1)*blocksize+iblocksize
         cwavefk(:,:)=cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor)
         do jblock=1,iblock
           do jblocksize=1,blocksize
             ll=(jblock-1)*blocksize+jblocksize
             if (ll<kk) then
               call dotprod_g(mkl(1,kk,ll),mkl(2,kk,ll),gs_hamk%istwf_k,npw_k,2,wi(:,1+(jblocksize-1)*npw_k*my_nspinor:&
&               jblocksize*npw_k*my_nspinor,jblock),cwavefk,mpi_enreg%me_g0,mpi_enreg%comm_fft)
             end if
           end do
         end do
       end do
     end do ! End of loop on block of bands

     ABI_FREE(cwavefk)
     mkl=-mkl

! Cholesky factorisation of -mkl=Lx(trans(L)*. On output mkl=L
     call zpotrf("L",nband_k,mkl,nband_k,info)

! calculate trans(L-1)
     ABI_MALLOC(bb,(2,nband_k,nband_k))
     bb=zero
     do kk=1,nband_k
       bb(1,kk,kk)=one
     end do
     call ztrtrs("L","T","N",nband_k,nband_k,mkl,nband_k,bb,nband_k,info)
     fock%fockACE(ikpt,isppol)%xi=zero

! Calculate ksi
     do kk=1,nband_k
       do jblock=1,nblockbd
         do jblocksize=1,blocksize
           ll=(jblock-1)*blocksize+jblocksize
           fock%fockACE(ikpt,isppol)%xi(1,:,kk)=fock%fockACE(ikpt,isppol)%xi(1,:,kk)+bb(1,ll,kk)*wi(1,1+(jblocksize-1)*&
&           npw_k*my_nspinor:jblocksize*npw_k*my_nspinor,jblock)-&
&           bb(2,ll,kk)*wi(2,1+(jblocksize-1)*npw_k*my_nspinor:jblocksize*npw_k*my_nspinor,jblock)
           fock%fockACE(ikpt,isppol)%xi(2,:,kk)=fock%fockACE(ikpt,isppol)%xi(2,:,kk)+bb(1,ll,kk)*wi(2,1+(jblocksize-1)*&
           npw_k*my_nspinor:jblocksize*npw_k*my_nspinor,jblock)+&
&           bb(2,ll,kk)*wi(1,1+(jblocksize-1)*npw_k*my_nspinor:jblocksize*npw_k*my_nspinor,jblock)
         end do
       end do
     end do

!    DEBUG
!    fock%fockACE(ikpt,isppol)%xi=zero
!    ENDDEBUG

     ABI_FREE(wi)
     ABI_FREE(mkl)

!    Restore the bandfft tabs
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_restoretabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kpg=kpg_k_sav)
     end if

!    Increment indices
     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       ibg=ibg+my_nspinor*nband_cprj_k
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if

     if (psps%usepaw==1) then
       call pawcprj_free(cwaveprj)
     end if
     ABI_FREE(cwaveprj)
     ABI_FREE(cwavef)
     ABI_FREE(bb)
     ABI_FREE(occblock)
     ABI_FREE(weight)
     ABI_FREE(ffnl)
     ABI_FREE(kg_k)
     ABI_FREE(kpg_k)
     ABI_FREE(ylm_k)
     ABI_FREE(ylmgr_k)
     ABI_FREE(ph3d)
   end do ! End k point loop
 end do ! End loop over spins

 call timab(1562,2,tsec)
 call timab(1565,1,tsec)

!Parallel case: accumulate (n,k) contributions
 if (xmpi_paral==1) then
   call xmpi_sum(fockcommon%e_fock0,spaceComm,ierr)
!  Forces
   if (optfor==1) then
     if (psps%usepaw==1) then
       call xmpi_sum(fockcommon%forces,spaceComm,ierr)
     end if
   end if
 end if

!need to reorder cprj=<p_lmn|Cnk> (from atom-sorted to unsorted)
 if (psps%usepaw==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx1)
 end if
!Deallocate temporary space
 call gs_hamk%free()

 call timab(1565,2,tsec)
 call timab(1560,2,tsec)

end subroutine fock2ACE
!!***

!!****f* ABINIT/fock_ACE_getghc
!! NAME
!!  fock_ACE_getghc
!!
!! FUNCTION
!!  Compute the matrix elements <G|Vx|psi> of the Fock operator in the ACE context.
!!
!! INPUTS
!!  cwavef(2,npw*nspinor*ndat)= planewave coefficients of wavefunctions on which Fock operator is applied.
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!!  mpi_enreg= information about MPI parallelization
!!
!! SIDE EFFECTS
!!  ghc(2,npw*ndat)= matrix elements <G|H|C> or <G|H-lambda.S|C> (if sij_opt>=0 or =-1 in getghc)
!!                   contains the fock exchange term for cwavef at the end.
!!
!! NOTES
!!  The current version assumes that :
!!   * nspinor = 1
!!   * no "my_nspinor"
!!   * no restriction to the value of istwfk_bz (but must be tested in all case)
!!   * all the data for the occupied states (cgocc_bz) are the same as those for the current states (cg)
!!
!! SOURCE

subroutine fock_ACE_getghc(cwavef,ghc,gs_ham,mpi_enreg,ndat,gpu_option)

!Arguments ------------------------------------
! Scalars
 integer :: ndat
 integer,optional :: gpu_option
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 real(dp),intent(inout) :: cwavef(:,:)!,ghc(2,gs_ham%npw_k*ndat)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 complex(dp), parameter :: cminusone  = (-1._dp,0._dp)
 integer :: iband,ikpt,ipw,my_nspinor,nband_k,npw,idat,gpu_option_
 real(dp) :: eigen
 type(fock_common_type),pointer :: fockcommon
! Arrays
 real(dp) :: tsec(2)
 real(dp), allocatable :: mat(:,:,:),ghc1(:,:),vdotr(:),vdoti(:)
 real(dp), ABI_CONTIGUOUS pointer :: xi(:,:,:)

! *************************************************************************

 call timab(1580,1,tsec)

 ABI_CHECK(associated(gs_ham%fockcommon),"fock must be associated!")
 fockcommon => gs_ham%fockcommon

 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")

 ikpt=fockcommon%ikpt
 npw=gs_ham%npw_k
 nband_k=fockcommon%nband(ikpt)
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
 gpu_option_=ABI_GPU_DISABLED; if(present(gpu_option)) gpu_option_ = gpu_option
!*Initialization of the array ghc1
!*ghc1 will contain the exact exchange contribution to the Hamiltonian
 ABI_MALLOC(ghc1,(2,npw*my_nspinor*ndat))
 ghc1=zero

 xi => gs_ham%fockACE_k%xi(:,:,:)

 if(gpu_option_==ABI_GPU_DISABLED) then
   if(gs_ham%istwf_k==1) then
     ABI_MALLOC(mat,(2,nband_k,ndat))
     call abi_zgemm_2r('C', 'N', nband_k, ndat, npw, cone, &
                       xi, npw, &
                       cwavef, npw, &
                       czero, &
                       mat, nband_k)
     call abi_zgemm_2r('N', 'N', npw, ndat, nband_k, cminusone, &
                       xi, npw, &
                       mat, nband_k, &
                       czero, &
                       ghc1, npw)
     ABI_FREE(mat)
   else
     ABI_MALLOC(vdotr,(nband_k))
     ABI_MALLOC(vdoti,(nband_k))
     do idat=1,ndat
       call dotprod_g_batch_half(vdotr,vdoti,gs_ham%istwf_k,npw*my_nspinor,nband_k,2,&
       &    cwavef(:,1+(idat-1)*npw:idat*npw),xi(:,:,:),mpi_enreg%me_g0,mpi_enreg%comm_fft)

       do iband=1, nband_k
         ghc1(1,1+(idat-1)*npw:idat*npw)=ghc1(1,1+(idat-1)*npw:idat*npw)-vdotr(iband)*xi(1,:,iband)
         ghc1(2,1+(idat-1)*npw:idat*npw)=ghc1(2,1+(idat-1)*npw:idat*npw)-vdotr(iband)*xi(2,:,iband)
       end do
     end do
     ABI_FREE(vdotr)
     ABI_FREE(vdoti)
   end if

   !* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
   ! ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_spkpt + ghc1(:,:)
   ghc(:,:)=ghc(:,:) + ghc1(:,:)

   ! call xmpi_sum(ghc,mpi_enreg%comm_kpt,ier)

 else if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(alloc:xi,ghc1)
   if(gs_ham%istwf_k==1) then
     ABI_MALLOC(mat,(2,nband_k,ndat))
     !$OMP TARGET ENTER DATA MAP(alloc:mat)
     !$OMP TARGET UPDATE TO(xi)
     !$OMP TARGET DATA USE_DEVICE_ADDR(xi,cwavef,mat,ghc1)
     call abi_gpu_xgemm(2, 'C', 'N', nband_k, ndat, npw, cone, &
                       c_loc(xi), npw, &
                       c_loc(cwavef), npw, &
                       czero, &
                       c_loc(mat), nband_k)
     call abi_gpu_xgemm(2, 'N', 'N', npw, ndat, nband_k, cminusone, &
                       c_loc(xi), npw, &
                       c_loc(mat), nband_k, &
                       czero, &
                       c_loc(ghc1), npw)
     !$OMP END TARGET DATA

     !$OMP TARGET EXIT DATA MAP(delete:mat)
     ABI_FREE(mat)
   else
     ABI_MALLOC(vdotr,(nband_k))
     ABI_MALLOC(vdoti,(nband_k))
     !$OMP TARGET ENTER DATA MAP(alloc:vdotr,vdoti)
     do idat=1,ndat
       call dotprod_g_batch_half(vdotr,vdoti,gs_ham%istwf_k,npw*my_nspinor,nband_k,2,&
       &    cwavef(:,1+(idat-1)*npw:idat*npw),xi(:,:,:),mpi_enreg%me_g0,mpi_enreg%comm_fft,&
       &    gpu_option=gpu_option_)

       !$OMP TARGET TEAMS DISTRIBUTE PRIVATE(iband) MAP(ghc1,xi,vdotr)
       do iband=1, nband_k
         !$OMP PARALLEL DO PRIVATE(ipw)
         do ipw=1,npw
           ghc1(1,ipw+(idat-1)*npw)=ghc1(1,ipw+(idat-1)*npw)-vdotr(iband)*xi(1,ipw,iband)
           ghc1(2,ipw+(idat-1)*npw)=ghc1(2,ipw+(idat-1)*npw)-vdotr(iband)*xi(2,ipw,iband)
         end do
       end do
     end do
     !$OMP TARGET EXIT DATA MAP(delete:vdotr,vdoti)
     ABI_FREE(vdotr)
     ABI_FREE(vdoti)
   end if

   !* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
   ! ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_spkpt + ghc1(:,:)

   !$OMP TARGET DATA USE_DEVICE_ADDR(ghc1,ghc)
   call abi_gpu_xaxpy(2,npw*ndat,cone,c_loc(ghc1),1,c_loc(ghc),1)
   !$OMP END TARGET DATA

   ! call xmpi_sum(ghc,mpi_enreg%comm_kpt,ier)
   !$OMP TARGET UPDATE FROM(ghc1)

   !$OMP TARGET EXIT DATA MAP(delete:xi,ghc1)
#endif
 end if


! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!* at each convergence step. (consistent definition with the definition of hartree energy)
 if (fockcommon%ieigen/=0) then
   do idat=1,ndat
     eigen=zero
!   * Dot product of cwavef and ghc
!   * inspired from the routine 54_spacepar/meanvalue_g but without the reference to parallelism and filtering
     if(gs_ham%istwf_k==2) then
       eigen=half*cwavef(1,1+(idat-1)*npw)*ghc1(1,1+(idat-1)*npw)
     else
       eigen=cwavef(1,1+(idat-1)*npw)*ghc1(1,1+(idat-1)*npw)+cwavef(2,1+(idat-1)*npw)*ghc1(2,1+(idat-1)*npw)
     end if
     do ipw=2,npw
       eigen=eigen+cwavef(1,ipw+(idat-1)*npw)*ghc1(1,ipw+(idat-1)*npw)+cwavef(2,ipw+(idat-1)*npw)*ghc1(2,ipw+(idat-1)*npw)
     end do
     if(gs_ham%istwf_k>=2) eigen=two*eigen
!    call xmpi_sum(eigen,mpi_enreg%comm_kpt,ier)
     fockcommon%eigen_ikpt(fockcommon%ieigen+idat-1)= eigen
     fockcommon%ieigen = 0
   end do ! idat
 end if

! ===============================
! === Deallocate local arrays ===
! ===============================

 ABI_FREE(ghc1)
 call timab(1580,2,tsec)

end subroutine fock_ACE_getghc
!!***

end module m_fock_getghc
!!***
