!!****m* ABINIT/m_fock_getghc
!! NAME
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (CMartins, FJ, MT, XG)
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

module m_fock_getghc

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_fock
 use m_pawcprj
 !use m_cgtools

 use defs_abitypes, only : mpi_type
 use defs_datatypes, only : pseudopotential_type
 use m_time,         only : timab
 use m_symtk,        only : matr3inv
 use m_cgtools,      only : dotprod_g
 use m_kg,           only : mkkpg
 use m_fftcore,      only : sphereboundary
 use m_fft,          only : fftpac, fourwf, fourdp
 use m_fstrings,     only : sjoin, itoa
 use m_hamiltonian,  only : gs_hamiltonian_type, K_H_KPRIME, init_hamiltonian
 use m_pawdij,       only : pawdijhat
 use m_paw_nhat,     only : pawmknhat_psipsi
 use m_spacepar,     only : hartre
 use m_nonlop,       only : nonlop
 use m_bandfft_kpt,      only : bandfft_kpt, bandfft_kpt_type, bandfft_kpt_savetabs,bandfft_kpt_restoretabs, &
                                prep_bandfft_tabs
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type
 use m_mkffnl,           only : mkffnl
 use m_mpinfo,           only : proc_distrb_cycle

 implicit none

 private
!!***

 public :: fock_getghc
 public :: fock2ACE
 public :: fock_ACE_getghc
!!***

contains
!!***

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
!! PARENTS
!!      m_fock_getghc,m_forstr,m_getghc
!!
!! CHILDREN
!!      dotprod_g
!!
!! SOURCE

subroutine fock_getghc(cwavef,cwaveprj,ghc,gs_ham,mpi_enreg)

!Arguments ------------------------------------
! Scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)
 real(dp),intent(inout) :: cwavef(:,:)!,ghc(2,gs_ham%npw_k)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 integer,parameter :: tim_fourwf0=0,tim_fourdp0=0,ndat1=1
 integer :: bdtot_jindex,choice,cplex_fock,cplex_dij,cpopt,i1,i2,i3,ia,iatom
 integer :: iband_cprj,ider,idir,idir1,ier,ii,ind,ipw,ifft,itypat,izero,jband,jbg,jcg,jkg
 integer :: jkpt,my_jsppol,jstwfk,lmn2_size,mgfftf,mpw,n1,n2,n3,n4,n5,n6
 integer :: n1f,n2f,n3f,n4f,n5f,n6f,natom,nband_k,ndij,nfft,nfftf,nfftotf,nhat12_grdim,nnlout
 integer :: npw,npwj,nspden_fock,nspinor,paw_opt,signs,tim_nonlop
 logical :: need_ghc,qeq0
 real(dp),parameter :: weight1=one
 real(dp) :: doti,eigen,imcwf,imcwocc,imvloc,invucvol,recwf,recwocc,revloc,occ,wtk
 type(fock_common_type),pointer :: fockcommon
 type(fock_BZ_type),pointer :: fockbz
! Arrays
 integer :: ngfft(18),ngfftf(18)
 integer,pointer :: gboundf(:,:),kg_occ(:,:),gbound_kp(:,:)
 real(dp) :: enlout_dum(1),dotr(6),fockstr(6),for1(3),qphon(3),qvec_j(3),tsec(2),gsc_dum(2,0),rhodum(2,1)
 real(dp) :: rhodum0(0,1,1),str(3,3)
 real(dp), allocatable :: dummytab(:,:),dijhat(:,:,:,:),dijhat_tmp(:,:),ffnl_kp_dum(:,:,:,:)
 real(dp), allocatable :: gvnlxc(:,:),ghc1(:,:),ghc2(:,:),grnhat12(:,:,:,:),grnhat_12(:,:,:,:,:),forikpt(:,:)
 real(dp), allocatable :: rho12(:,:,:),rhog_munu(:,:),rhor_munu(:,:),vlocpsi_r(:)
 real(dp), allocatable :: vfock(:),psilocal(:,:,:),vectin_dum(:,:),vqg(:),forout(:,:),strout(:,:)
 real(dp), allocatable,target ::cwavef_r(:,:,:,:)
 real(dp), ABI_CONTIGUOUS  pointer :: cwaveocc_r(:,:,:,:)
 type(pawcprj_type),pointer :: cwaveocc_prj(:,:)

 real(dp) :: rprimd(3,3),for12(3)

! *************************************************************************
!return
 call timab(1504,1,tsec)
 call timab(1505,1,tsec)
 call timab(1515,1,tsec)

 ABI_CHECK(associated(gs_ham%fockcommon),"fock_common must be associated!")
 fockcommon => gs_ham%fockcommon
 ABI_CHECK(associated(gs_ham%fockbz),"fock_bz must be associated!")
 fockbz => gs_ham%fockbz

 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")
 if (fockcommon%usepaw==1) then
   ABI_CHECK((size(cwaveprj,1)==gs_ham%natom.and.size(cwaveprj,2)==gs_ham%nspinor),"error on cwaveprj dims")
 end if
 need_ghc=(size(ghc,2)>0)

!Some constants
 invucvol=1.d0/sqrt(gs_ham%ucvol)
 call matr3inv(gs_ham%gprimd,rprimd)
 cplex_fock=2;nspden_fock=1
 natom=fockcommon%natom
 nspinor=gs_ham%nspinor
 mpw=maxval(fockbz%npwarr)
 npw=gs_ham%npw_k
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
 ABI_ALLOCATE(cwavef_r,(2,n4f,n5f,n6f))
!*rhormunu = overlap matrix between cwavef and (jkpt,mu) in R-space
 ABI_ALLOCATE(rhor_munu,(cplex_fock,nfftf))
!*rhogmunu = overlap matrix between cwavef and (jkpt,mu) in G-space
 ABI_ALLOCATE(rhog_munu,(2,nfftf))
!*dummytab = variables for fourwf
 ABI_ALLOCATE(dummytab,(2,nfft))
!*vfock = Fock potential
 ABI_ALLOCATE(vfock,(cplex_fock*nfftf))
!*vqg = 4pi/(G+q)**2
 ABI_ALLOCATE(vqg,(nfftf))

!*Initialization of the array ghc1
!*ghc1 will contain the exact exchange contribution to the Hamiltonian
 ABI_ALLOCATE(ghc1,(2,npw))
 ABI_ALLOCATE(ghc2,(2,npw))
 ghc1=zero
 ghc2=zero
!*Initialization of the array vlocpsi_r
!*vlocpsi_r = partial local Fock operator applied to cwavef in r-space and summed over all occupied (jkpt,mu)
 ABI_ALLOCATE(vlocpsi_r,(cplex_fock*nfftf))
 vlocpsi_r=zero

!*Additional arrays in case of paw
 if (fockcommon%usepaw==1) then
   nhat12_grdim=0
   if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
     ider=3
     ABI_ALLOCATE(strout,(2,npw*nspinor))
   end if
   if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
     ider=3
     ABI_ALLOCATE(forout,(2,npw*nspinor))
     ABI_ALLOCATE(forikpt,(3,natom))
     forikpt=zero
   end if
   ABI_ALLOCATE(grnhat_12,(2,nfftf,nspinor**2,3,natom*(ider/3)))
   ABI_ALLOCATE(gvnlxc,(2,npw*nspinor))
   ABI_ALLOCATE(grnhat12,(2,nfftf,nspinor**2,3*nhat12_grdim))
 end if

 if (fockcommon%usepaw==1.or.fockcommon%optstr) then
   ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
   call sphereboundary(gboundf,gs_ham%istwf_k,gs_ham%kg_k,mgfftf,npw)
 else
   gboundf=>gs_ham%gbound_k
 end if
 call timab(1515,2,tsec)
! ==========================================
! === Get cwavef in real space using FFT ===
! ==========================================
 call timab(840+tim_fourwf0,1,tsec)
 cwavef_r=zero
 call fourwf(0,rhodum0,cwavef,rhodum,cwavef_r,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
& mgfftf,mpi_enreg,ndat1,ngfftf,npw,1,n4f,n5f,n6f,0,tim_fourwf0,weight1,weight1,&
& use_gpu_cuda=gs_ham%use_gpu_cuda)
 cwavef_r=cwavef_r*invucvol
 call timab(840+tim_fourwf0,2,tsec)

! =====================================================
! === Select the states in cgocc_bz with the same spin ===
! =====================================================
!* Initialization of the indices/shifts, according to the value of isppol
!* bdtot_jindex = shift to be applied on the location of data in the array occ_bz ?
 call timab(1515,1,tsec)
 bdtot_jindex=0
!* jbg = shift to be applied on the location of data in the array cprj/occ
 jbg=0;jcg=0
 my_jsppol=fockcommon%isppol
 if((fockcommon%isppol==2).and.(mpi_enreg%nproc_kpt/=1)) my_jsppol=1

 call timab(1505,2,tsec)
 call timab(1506,1,tsec)
 call timab(1515,2,tsec)
!===================================
!=== Loop on the k-points in IBZ ===
!===================================
 call timab(1515,1,tsec)
 jkg=0

 if (associated(gs_ham%ph3d_kp)) then
   nullify (gs_ham%ph3d_kp)
 end if

 do jkpt=1,fockbz%mkpt
!* nband_k = number of bands at point k_j
   nband_k=fockbz%nbandocc_bz(jkpt,my_jsppol)
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
!* Load k^prime hamiltonian in the gs_ham datastructure
!  Note: ffnl_kp / ph3d_kp / gbound_kp are not used

   if (associated(gs_ham%ph3d_kp)) then
     ABI_ALLOCATE(gs_ham%ph3d_kp,(2,npwj,gs_ham%matblk))
   end if

   call gs_ham%load_kprime(kpt_kp=fockbz%kptns_bz(:,jkpt),&
&   istwf_kp=jstwfk,npw_kp=npwj,kg_kp=fockbz%kg_bz(:,1+jkg:npwj+jkg))
!* Some temporary allocations needed for PAW
   if (fockcommon%usepaw==1) then
     ABI_ALLOCATE(vectin_dum,(2,npwj*nspinor))
     vectin_dum=zero
     ABI_ALLOCATE(ffnl_kp_dum,(npwj,0,gs_ham%lmnmax,gs_ham%ntypat))
     call gs_ham%load_kprime(ffnl_kp=ffnl_kp_dum)
   end if
   call timab(1515,2,tsec)

! ======================================
! === Calculate the vector q=k_i-k_j ===
! ======================================
!* Evaluation of kpoint_j, the considered k-point in reduced coordinates
!     kpoint_j(:)=fockbz%kptns_bz(:,jkpt)
!* the vector qvec is expressed in reduced coordinates.
!     qvec(:)=kpoint_i(:)-kpoint_j(:)
   call timab(1515,1,tsec)
   qvec_j(:)=gs_ham%kpt_k(:)-fockbz%kptns_bz(:,jkpt)
   qeq0=(qvec_j(1)**2+qvec_j(2)**2+qvec_j(3)**2<1.d-15)
   call bare_vqg(qvec_j,fockcommon%gsqcut,gs_ham%gmet,fockcommon%usepaw,fockcommon%hyb_mixing,&
&   fockcommon%hyb_mixing_sr,fockcommon%hyb_range_fock,nfftf,fockbz%nkpt_bz,ngfftf,gs_ham%ucvol,vqg)
   call timab(1515,2,tsec)

! =================================================
! === Loop on the band indices jband of cgocc_k ===
! =================================================
   call timab(1515,1,tsec)
   do jband=1,nband_k

!*   occ = occupancy of jband at this k point
     occ=fockbz%occ_bz(jband+bdtot_jindex,my_jsppol)
     if(occ<tol8) cycle
   call timab(1515,2,tsec)

! ==============================================
! === Get cwaveocc_r in real space using FFT ===
! ==============================================
     if (allocated(fockbz%cwaveocc_bz)) then
       cwaveocc_r => fockbz%cwaveocc_bz(:,:,:,:,jband+jbg,my_jsppol)
     else
       ABI_ALLOCATE(cwaveocc_r,(2,n4f,n5f,n6f))
       cwaveocc_r=zero
       call timab(840+tim_fourwf0,1,tsec)
       call fourwf(1,rhodum0,fockbz%cgocc(:,1+jcg+npwj*(jband-1):jcg+jband*npwj,my_jsppol),rhodum,cwaveocc_r, &
&       gbound_kp,gbound_kp,jstwfk,kg_occ,kg_occ,mgfftf,mpi_enreg,ndat1,ngfftf,&
&       npwj,1,n4f,n5f,n6f,tim_fourwf0,0,weight1,weight1,use_gpu_cuda=gs_ham%use_gpu_cuda)
       cwaveocc_r=cwaveocc_r*invucvol
       call timab(840+tim_fourwf0,2,tsec)
     end if

! ================================================
! === Get the overlap density matrix rhor_munu ===
! ================================================
!* Calculate the overlap density matrix in real space = conj(cwaveocc_r)*cwavef_r
!* rhor_munu will contain the overlap density matrix.
! vfock=-int{conj(cwaveocc_r)*cwavef_r*dr'/|r-r'|}

     call timab(1508,1,tsec)
     call timab(1515,1,tsec)
     ind=0
     do i3=1,n3f
       do i2=1,n2f
         do i1=1,n1f
           ind=ind+1
           recwf  =cwavef_r(1,i1,i2,i3)   ; imcwf  =cwavef_r(2,i1,i2,i3)
           recwocc=cwaveocc_r(1,i1,i2,i3) ; imcwocc=cwaveocc_r(2,i1,i2,i3)
           rhor_munu(1,ind)= recwocc*recwf+imcwocc*imcwf
           rhor_munu(2,ind)= recwocc*imcwf-imcwocc*recwf
         end do ! i1
       end do ! i2
     end do ! i3
     call timab(1508,2,tsec)

! =======================================================
! === Add compensation charge density in the PAW case ===
! === Get the overlap density matrix rhor_munu        ===
! =======================================================
     call timab(1509,1,tsec)
     if (fockcommon%usepaw==1) then
       iband_cprj=(my_jsppol-1)*fockbz%mkptband+jbg+jband

       ABI_ALLOCATE(rho12,(2,nfftf,nspinor**2))

       cwaveocc_prj=>fockbz%cwaveocc_prj(:,iband_cprj:iband_cprj+nspinor-1)

       call pawmknhat_psipsi(cwaveprj,cwaveocc_prj,ider,izero,natom,natom,nfftf,ngfftf,&
&       nhat12_grdim,nspinor,fockcommon%ntypat,fockbz%pawang,fockcommon%pawfgrtab,grnhat12,rho12,&
&       fockcommon%pawtab,gprimd=gs_ham%gprimd,grnhat_12=grnhat_12,qphon=qvec_j,xred=gs_ham%xred,atindx=gs_ham%atindx)

       rhor_munu(1,:)=rhor_munu(1,:)+rho12(1,:,nspinor)
       rhor_munu(2,:)=rhor_munu(2,:)-rho12(2,:,nspinor)
     end if
     call timab(1515,2,tsec)
     ! Perform an FFT using fourwf to get rhog_munu = FFT^-1(rhor_munu)
     call timab(260+tim_fourdp0,1,tsec)
     call fourdp(cplex_fock,rhog_munu,rhor_munu,-1,mpi_enreg,nfftf,1,ngfftf,tim_fourdp0)
     call timab(260+tim_fourdp0,1,tsec)
     call timab(1509,2,tsec)

     if(fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
       call strfock(gs_ham%gprimd,fockcommon%gsqcut,fockstr,fockcommon%hyb_mixing,fockcommon%hyb_mixing_sr,&
&       fockcommon%hyb_range_fock,mpi_enreg,nfftf,ngfftf,fockbz%nkpt_bz,rhog_munu,gs_ham%ucvol,qvec_j)
       fockcommon%stress_ikpt(:,fockcommon%ieigen)=fockcommon%stress_ikpt(:,fockcommon%ieigen)+fockstr(:)*occ*wtk
       if (fockcommon%usepaw==0.and.(.not.need_ghc)) then
         if (allocated(fockbz%cgocc)) then
           ABI_DEALLOCATE(cwaveocc_r)
         end if
         cycle
       end if
     end if

! ===================================================
! === Calculate the local potential vfockloc_munu ===
! ===================================================
!* Apply the Poisson solver to "rhog_munu" while taking into account the effect of the vector "qvec"
!* This is precisely what is done in the subroutine hartre, with option cplex=2.
!* vfock will contain the local Fock potential, the result of hartre routine.
!* vfock = FFT( rhog_munu/|g+qvec|^2 )
     call timab(1510,1,tsec)
     call timab(1515,1,tsec)
#if 0

     call hartre(cplex_fock,fockcommon%gsqcut,fockcommon%usepaw,mpi_enreg,nfftf,ngfftf,&
&     mpi_enreg%paral_kgb,rhog_munu,rprimd,vfock,divgq0=fock%divgq0,qpt=qvec_j)

#else
     do ifft=1,nfftf
       rhog_munu(1,ifft) = rhog_munu(1,ifft) * vqg(ifft)
       rhog_munu(2,ifft) = rhog_munu(2,ifft) * vqg(ifft)
     end do

     call timab(1515,2,tsec)
     call timab(260+tim_fourdp0,1,tsec)
     call fourdp(cplex_fock,rhog_munu,vfock,+1,mpi_enreg,nfftf,1,ngfftf,tim_fourdp0)
     call timab(260+tim_fourdp0,2,tsec)

#endif
     call timab(1510,2,tsec)

!===============================================================
!======== Calculate Dij_Fock_hat contribution in case of PAW ===
!===============================================================

     if (fockcommon%usepaw==1) then
       call timab(1515,1,tsec)
       qphon=qvec_j;nfftotf=product(ngfftf(1:3))
       cplex_dij=1;ndij=nspden_fock
       ABI_ALLOCATE(dijhat,(cplex_dij*gs_ham%dimekb1,natom,ndij,cplex_fock))
       dijhat=zero
       do iatom=1,natom
         itypat=gs_ham%typat(iatom)
         lmn2_size=fockcommon%pawtab(itypat)%lmn2_size
         ABI_ALLOCATE(dijhat_tmp,(cplex_fock*cplex_dij*lmn2_size,ndij))
         dijhat_tmp=zero
         call pawdijhat(dijhat_tmp,cplex_dij,cplex_fock,gs_ham%gprimd,iatom,&
&         natom,ndij,nfftf,nfftotf,nspden_fock,nspden_fock,fockbz%pawang,fockcommon%pawfgrtab(iatom),&
&         fockcommon%pawtab(itypat),vfock,qphon,gs_ham%ucvol,gs_ham%xred)
         do ii=1,cplex_fock
           ind=(ii-1)*lmn2_size
           dijhat(1:cplex_dij*lmn2_size,iatom,:,ii)=dijhat_tmp(ind+1:ind+cplex_dij*lmn2_size,:)
         end do
         ABI_DEALLOCATE(dijhat_tmp)
       end do
       signs=2; cpopt=2;idir=0; paw_opt=1;nnlout=1;tim_nonlop=1
       
       if(need_ghc) then
         choice=1
         call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&         ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,gvnlxc,enl=dijhat,&
&         select_k=K_H_KPRIME)
         ghc2=ghc2-gvnlxc*occ*wtk       
       end if
       call timab(1515,2,tsec)

! Forces calculation

       if (fockcommon%optfor.and.(fockcommon%ieigen/=0)) then
         call timab(1515,1,tsec)
         choice=2; dotr=zero;doti=zero;cpopt=4
         do iatom=1,natom
           do idir=1,3
             call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&             ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
&             forout,enl=dijhat,iatom_only=iatom,&
&             select_k=K_H_KPRIME)
             call dotprod_g(dotr(idir),doti,gs_ham%istwf_k,npw,2,cwavef,forout,mpi_enreg%me_g0,mpi_enreg%comm_fft)
             for1(idir)=zero
             do ifft=1,fockcommon%pawfgrtab(iatom)%nfgd
               ind=fockcommon%pawfgrtab(iatom)%ifftsph(ifft)
               for1(idir)=for1(idir)+vfock(2*ind-1)*grnhat_12(1,ind,1,idir,iatom)-&
&               vfock(2*ind)*grnhat_12(2,ind,1,idir,iatom)
             end do
           end do
           do idir=1,3
             for12(idir)=rprimd(1,idir)*for1(1)+rprimd(2,idir)*for1(2)+rprimd(3,idir)*for1(3)
             forikpt(idir,iatom)=forikpt(idir,iatom)-(for12(idir)*gs_ham%ucvol/nfftf+dotr(idir))*occ*wtk
           end do
         end do
         call timab(1515,2,tsec)
       end if

! Stresses calculation
       if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
         signs=2;choice=3;cpopt=4

       ! first contribution
         dotr=zero
         do idir=1,6
           call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&           ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
&           strout,enl=dijhat,select_k=K_H_KPRIME)
           call dotprod_g(dotr(idir),doti,gs_ham%istwf_k,npw,2,cwavef,strout,mpi_enreg%me_g0,mpi_enreg%comm_fft)
           fockcommon%stress_ikpt(idir,fockcommon%ieigen)=fockcommon%stress_ikpt(idir,fockcommon%ieigen)-&
&           dotr(idir)*occ*wtk/gs_ham%ucvol
           call timab(1515,2,tsec)
         end do
       ! second contribution
         call timab(1515,1,tsec)
         str=zero
         do iatom=1,natom
           do idir=1,3
             do idir1=1,3
               do ifft=1,fockcommon%pawfgrtab(iatom)%nfgd
                 ind=fockcommon%pawfgrtab(iatom)%ifftsph(ifft)
                 str(idir,idir1)=str(idir,idir1)+(vfock(2*ind-1)*grnhat_12(1,ind,1,idir,iatom)-&
&                 vfock(2*ind)*grnhat_12(2,ind,1,idir,iatom))*fockcommon%pawfgrtab(iatom)%rfgd(idir1,ifft)

               end do
             end do
           end do
         end do
         do idir=1,3
           fockstr(idir)=str(idir,idir)
         end do
         fockstr(4)=(str(3,2)+str(2,3))*half
         fockstr(5)=(str(3,1)+str(1,3))*half
         fockstr(6)=(str(1,2)+str(2,1))*half
         do idir=1,6
           fockcommon%stress_ikpt(idir,fockcommon%ieigen)=fockcommon%stress_ikpt(idir,fockcommon%ieigen)+&
&           fockstr(idir)/nfftf*occ*wtk
         end do

       ! third contribution
         doti=zero
         do ifft=1,nfftf
           doti=doti+vfock(2*ifft-1)*rho12(1,ifft,nspinor)-vfock(2*ifft)*rho12(2,ifft,nspinor)
         end do
         fockcommon%stress_ikpt(1:3,fockcommon%ieigen)=fockcommon%stress_ikpt(1:3,fockcommon%ieigen)-doti/nfftf*occ*wtk
         call timab(1515,2,tsec)
!         doti=zero
!         do ifft=1,nfftf
!           doti=doti+vfock(2*ifft-1)*rhor_munu(1,ifft)-vfock(2*ifft)*rhor_munu(2,ifft)
!         end do
!         fockcommon%stress_ikpt(1:3,fockcommon%ieigen)=fockcommon%stress_ikpt(1:3,fockcommon%ieigen)+doti/nfftf*occ*wtk*half
       end if ! end stresses

       ABI_DEALLOCATE(dijhat)
       ABI_DEALLOCATE(rho12)
     end if !end PAW

! =============================================================
! === Apply the local potential vfockloc_munu to cwaveocc_r ===
! =============================================================
     call timab(1507,1,tsec)
     call timab(1515,1,tsec)
     ind=0
     do i3=1,ngfftf(3)
       do i2=1,ngfftf(2)
         do i1=1,ngfftf(1)
           ind=ind+1
!          ind=i1+ngfftf(1)*(i2-1+ngfftf(2)*(i3-1))
           revloc=vfock(2*ind-1) ; imvloc=vfock(2*ind)
           recwocc=cwaveocc_r(1,i1,i2,i3) ; imcwocc=cwaveocc_r(2,i1,i2,i3)
           vlocpsi_r(2*ind-1)=vlocpsi_r(2*ind-1)-(revloc*recwocc-imvloc*imcwocc)*occ*wtk
           vlocpsi_r(2*ind  )=vlocpsi_r(2*ind  )-(revloc*imcwocc+imvloc*recwocc)*occ*wtk
         end do
       end do
     end do
     call timab(1507,2,tsec)
     if (allocated(fockbz%cgocc)) then
       ABI_DEALLOCATE(cwaveocc_r)
     end if

   end do ! jband

! ================================================
! === End : update of shifts and deallocations ===
! ================================================
!* Update of the shifts to be applied (reminder : mkmem is not 0, nspinor=1)
   jcg=jcg+npwj*nband_k
   jbg=jbg+nband_k
   bdtot_jindex=bdtot_jindex+nband_k
   jkg=jkg+npwj
   if (fockcommon%usepaw==1) then
     ABI_DEALLOCATE(vectin_dum)
     ABI_DEALLOCATE(ffnl_kp_dum)
   end if
   if (associated(gs_ham%ph3d_kp)) then
     ABI_DEALLOCATE(gs_ham%ph3d_kp)
   end if
 end do ! jkpt

 if (fockcommon%usepaw==1) then
   if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
     call xmpi_sum(forikpt,mpi_enreg%comm_hf,ier)
     do iatom=1,natom !Loop over atom
       ia=gs_ham%atindx(iatom)
       fockcommon%forces_ikpt(:,ia,fockcommon%ieigen)=forikpt(:,iatom)
     end do
   end if
 end if
 if(fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
   call xmpi_sum(fockcommon%stress_ikpt,mpi_enreg%comm_hf,ier)
 end if

 if (.not.need_ghc) then

! ===============================
! === Deallocate local arrays ===
! ===============================
   ABI_DEALLOCATE(cwavef_r)
   ABI_DEALLOCATE(ghc1)
   ABI_DEALLOCATE(ghc2)
   ABI_DEALLOCATE(rhor_munu)
   ABI_DEALLOCATE(rhog_munu)
   ABI_DEALLOCATE(vlocpsi_r)
   ABI_DEALLOCATE(dummytab)
   ABI_DEALLOCATE(vfock)
   ABI_DEALLOCATE(vqg)
   if (fockcommon%usepaw==1) then
     ABI_DEALLOCATE(gvnlxc)
     ABI_DEALLOCATE(grnhat12)
     if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
       ABI_DEALLOCATE(forikpt)
       ABI_DEALLOCATE(forout)
     end if
     if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
       ABI_DEALLOCATE(strout)
     end if
     ABI_DEALLOCATE(grnhat_12)
   end if
   if(fockcommon%usepaw==1.or.fockcommon%optstr) then
     ABI_DEALLOCATE(gboundf)
   end if
!*Restore gs_ham datastructure

   if (associated(gs_ham%ph3d_kp)) then
     ABI_ALLOCATE(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
   end if
   call gs_ham%load_kprime(kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
&   npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!   if (fockcommon%ieigen/=0) fockcommon%ieigen=0
   return
 end if

 call timab(1506,2,tsec)
 call timab(1511,1,tsec)

!*Restore gs_ham datastructure

 if (associated(gs_ham%ph3d_kp)) then
   ABI_ALLOCATE(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
 end if
 call gs_ham%load_kprime(kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
& npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!* Perform an FFT using fourwf to get ghc1 = FFT^-1(vlocpsi_r)
 ABI_ALLOCATE(psilocal,(cplex_fock*n4f,n5f,n6f))
 call fftpac(1,mpi_enreg,nspden_fock,cplex_fock*n1f,n2f,n3f,cplex_fock*n4f,n5f,n6f,ngfft,vlocpsi_r,psilocal,2)

 call timab(1515,2,tsec)

 call timab(840+tim_fourwf0,1,tsec)
 call fourwf(0,rhodum0,rhodum,ghc1,psilocal,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
& mgfftf,mpi_enreg,ndat1,ngfftf,1,npw,n4f,n5f,n6f,3,tim_fourwf0,weight1,weight1,&
& use_gpu_cuda=gs_ham%use_gpu_cuda)
 call timab(840+tim_fourwf0,2,tsec)
 ABI_DEALLOCATE(psilocal)

 ghc1=ghc1*sqrt(gs_ham%ucvol)+ghc2

!* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
 ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_hf + ghc1(:,:)

 call xmpi_sum(ghc,mpi_enreg%comm_hf,ier)

 call timab(1511,2,tsec)


! ===============================
! === Deallocate local PAW arrays ===
! ===============================
 call timab(1515,1,tsec)
 if (fockcommon%usepaw==1) then
   ABI_DEALLOCATE(gvnlxc)
   ABI_DEALLOCATE(grnhat12)
   if ((fockcommon%optfor).and.(fockcommon%ieigen/=0)) then
     ABI_DEALLOCATE(forikpt)
     ABI_DEALLOCATE(forout)
   end if
   if (fockcommon%optstr.and.(fockcommon%ieigen/=0)) then
     ABI_DEALLOCATE(strout)
   end if
   ABI_DEALLOCATE(grnhat_12)
 end if
 if(fockcommon%usepaw==1.or.fockcommon%optstr) then
   ABI_DEALLOCATE(gboundf)
 end if
 call timab(1515,2,tsec)
! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!* at each convergence step. (consistent definition with the definition of hartree energy)
 call timab(1515,1,tsec)
 if (fockcommon%ieigen/=0) then
   eigen=zero
!* Dot product of cwavef and ghc
!* inspired from the routine 54_spacepar/meanvalue_g but without the reference to parallelism and filtering
   if(gs_ham%istwf_k==2) then
     eigen=half*cwavef(1,1)*ghc1(1,1)
   else
     eigen=cwavef(1,1)*ghc1(1,1)+cwavef(2,1)*ghc1(2,1)
   end if
   do ipw=2,npw
     eigen=eigen+cwavef(1,ipw)*ghc1(1,ipw)+cwavef(2,ipw)*ghc1(2,ipw)
   end do
   if(gs_ham%istwf_k>=2) eigen=two*eigen
   call xmpi_sum(eigen,mpi_enreg%comm_hf,ier)
   fockcommon%eigen_ikpt(fockcommon%ieigen)= eigen
   if(fockcommon%use_ACE==0) fockcommon%ieigen = 0
 end if

! ===============================
! === Deallocate local arrays ===
! ===============================
 ABI_DEALLOCATE(cwavef_r)
 ABI_DEALLOCATE(ghc1)
 ABI_DEALLOCATE(ghc2)
 ABI_DEALLOCATE(rhor_munu)
 ABI_DEALLOCATE(rhog_munu)
 ABI_DEALLOCATE(vlocpsi_r)
 ABI_DEALLOCATE(dummytab)
 ABI_DEALLOCATE(vfock)
 ABI_DEALLOCATE(vqg)

 call timab(1504,2,tsec)
 call timab(1515,2,tsec)

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
!!  use_gpu_cuda= 0 or 1 to know if we use cuda for nonlop call
!!  wtk(nkpt)=weight associated with each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!
!! fock%fockACE(ikpt,isppol)%xi
!! if optfor=1, fock%fock_common%forces
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      dotprod_g
!!
!! SOURCE

subroutine fock2ACE(cg,cprj,fock,istwfk,kg,kpt,mband,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,my_natom,natom,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,&
&  ntypat,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,typat,usecprj,use_gpu_cuda,wtk,xred,ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfft,mkmem,mpsang,mpw,my_natom,natom,nfft,nkpt
 integer,intent(in) :: nspden,nsppol,nspinor,ntypat,optfor
 integer,intent(in) :: usecprj,use_gpu_cuda
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
 integer :: mband_cprj,me_distrb,my_ikpt,my_nspinor,nband_k,nband_cprj_k,ndat,nkpg
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
 type(pawcprj_type),pointer :: cwaveprj_idat(:,:)

!*************************************************************************

 call timab(920,1,tsec)
 call timab(921,1,tsec)

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

 call init_hamiltonian(gs_hamk,psps,pawtab,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,usecprj=usecprj,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& paw_ij=paw_ij,ph1d=ph1d,fock=fock,&
& use_gpu_cuda=use_gpu_cuda)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)
 fockcommon%use_ACE=use_ACE_old
 call timab(921,2,tsec)

!need to reorder cprj=<p_lmn|Cnk> (from unsorted to atom-sorted)
 if (psps%usepaw==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx)
 end if

!LOOP OVER SPINS
 bdtot_index=0;ibg=0;icg=0

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

     call timab(922,1,tsec)

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

     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
     if (psps%usepaw==1) then
       ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor*bandpp))
       call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
     else
       ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
     end if

     ABI_ALLOCATE(kg_k,(3,mpw))
!$OMP PARALLEL DO
     do ipw=1,npw_k
       kg_k(:,ipw)=kg(:,ipw+ikg)
     end do

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr_k,(0,0,0))
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
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if


!    Compute nonlocal form factors ffnl at all (k+G)
     ider=0;idir=0;dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&     ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&     nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!    Load k-dependent part in the Hamiltonian datastructure
!     - Compute 3D phase factors
!     - Prepare various tabs in case of band-FFT parallelism
!     - Load k-dependent quantities in the Hamiltonian

     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
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

     call timab(922,2,tsec)

!    The following is now wrong. In sequential, nblockbd=nband_k/bandpp
!    blocksize= bandpp (JB 2016/04/16)
!    Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
!
     ABI_ALLOCATE(occblock,(blocksize))
     ABI_ALLOCATE(weight,(blocksize))
     occblock=zero;weight=zero

     if (fockcommon%optfor) then
       fockcommon%forces_ikpt=zero
     end if

     ABI_ALLOCATE(wi,(2,npw_k*my_nspinor*blocksize,nblockbd))
     wi=zero
     ABI_ALLOCATE(mkl,(2,nband_k,nband_k))
     mkl=zero
! Calculate all the Wi for the current k-point

     do iblock=1,nblockbd

       iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband_k)
       iband_cprj=(iblock-1)*bandpp+1
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,me_distrb)) cycle

!      Select occupied bandsddk
       occblock(:)=occ(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
       call timab(926,1,tsec)
       weight(:)=wtk(ikpt)*occblock(:)

!        Load contribution from n,k
       cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&       cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
       if (psps%usepaw==1) then
         call pawcprj_get(gs_hamk%atindx1,cwaveprj,cprj,natom,iband_cprj,ibg,ikpt,0,isppol,&
&         mband_cprj,mkmem,natom,bandpp,nband_cprj_k,my_nspinor,nsppol,0,&
&         mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
       end if

       call timab(926,2,tsec)

       if (mpi_enreg%paral_kgb==1) then
         msg='fock2ACE: Paral_kgb is not yet implemented for fock calculations'
         MSG_BUG(msg)
       end if
       ndat=mpi_enreg%bandpp
       if (gs_hamk%usepaw==0) cwaveprj_idat => cwaveprj

       do iblocksize=1,blocksize
         fockcommon%ieigen=(iblock-1)*blocksize+iblocksize
         fockcommon%iband=(iblock-1)*blocksize+iblocksize
         if (gs_hamk%usepaw==1) then
           cwaveprj_idat => cwaveprj(:,(iblocksize-1)*my_nspinor+1:iblocksize*my_nspinor)
         end if
         call fock_getghc(cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),cwaveprj_idat,&
&         wi(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor,iblock),gs_hamk,mpi_enreg)
         mkl(1,fockcommon%ieigen,fockcommon%ieigen)=fockcommon%eigen_ikpt(fockcommon%ieigen)
         fockcommon%e_fock0=fockcommon%e_fock0+half*weight(iblocksize)*fockcommon%eigen_ikpt(fockcommon%ieigen)
         if (fockcommon%optfor) then
           fockcommon%forces(:,:)=fockcommon%forces(:,:)+weight(iblocksize)*fockcommon%forces_ikpt(:,:,fockcommon%ieigen)
         end if
       end do


     end do ! End of loop on block of bands

! Calculate Mkl for the current k-point
     ABI_ALLOCATE(cwavefk,(2,npw_k*my_nspinor))
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

     ABI_DEALLOCATE(cwavefk)
     mkl=-mkl

! Cholesky factorisation of -mkl=Lx(trans(L)*. On output mkl=L
     call zpotrf("L",nband_k,mkl,nband_k,info)

! calculate trans(L-1)
     ABI_ALLOCATE(bb,(2,nband_k,nband_k))
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

     ABI_DEALLOCATE(wi)
     ABI_DEALLOCATE(mkl)

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
     ABI_DATATYPE_DEALLOCATE(cwaveprj)
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(bb)
     ABI_DEALLOCATE(occblock)
     ABI_DEALLOCATE(weight)
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(ph3d)
   end do ! End k point loop
 end do ! End loop over spins

!Parallel case: accumulate (n,k) contributions
 if (xmpi_paral==1) then
   call xmpi_sum(fockcommon%e_fock0,spaceComm,ierr)
!  Forces
   if (optfor==1) then
     call timab(65,2,tsec)
     if (psps%usepaw==1) then
       call xmpi_sum(fockcommon%forces,spaceComm,ierr)
     end if
   end if
 end if

 call timab(925,1,tsec)

!need to reorder cprj=<p_lmn|Cnk> (from atom-sorted to unsorted)
 if (psps%usepaw==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx1)
 end if
!Deallocate temporary space
 call gs_hamk%free()

 call timab(925,2,tsec)
 call timab(920,2,tsec)

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
!! PARENTS
!!      m_getghc
!!
!! CHILDREN
!!      dotprod_g
!!
!! SOURCE

subroutine fock_ACE_getghc(cwavef,ghc,gs_ham,mpi_enreg)

!Arguments ------------------------------------
! Scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 real(dp),intent(inout) :: cwavef(:,:)!,ghc(2,gs_ham%npw_k)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 integer :: iband,ikpt,ipw,my_nspinor,nband_k,npw
 real(dp) :: doti,dotr,eigen
 type(fock_common_type),pointer :: fockcommon
! Arrays
 real(dp), allocatable :: ghc1(:,:),xi(:,:)

! *************************************************************************

 ABI_CHECK(associated(gs_ham%fockcommon),"fock must be associated!")
 fockcommon => gs_ham%fockcommon

 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")

 ikpt=fockcommon%ikpt
 npw=gs_ham%npw_k
 nband_k=fockcommon%nband(ikpt)
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
!*Initialization of the array ghc1
!*ghc1 will contain the exact exchange contribution to the Hamiltonian
 ABI_ALLOCATE(ghc1,(2,npw*my_nspinor))
 ghc1=zero
 ABI_ALLOCATE(xi,(2,npw*my_nspinor))

 do iband=1, nband_k
   xi(1,:)=gs_ham%fockACE_k%xi(1,:,iband)
   xi(2,:)=gs_ham%fockACE_k%xi(2,:,iband)

   call dotprod_g(dotr,doti,gs_ham%istwf_k,npw*my_nspinor,2,xi,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_fft)

   ghc1(1,:)=ghc1(1,:)-(dotr*gs_ham%fockACE_k%xi(1,:,iband)-doti*gs_ham%fockACE_k%xi(2,:,iband))
   ghc1(2,:)=ghc1(2,:)-(dotr*gs_ham%fockACE_k%xi(2,:,iband)+doti*gs_ham%fockACE_k%xi(1,:,iband))
 end do
 ABI_DEALLOCATE(xi)

!* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
! ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_kpt + ghc1(:,:)
 ghc(:,:)=ghc(:,:) + ghc1(:,:)

! call xmpi_sum(ghc,mpi_enreg%comm_kpt,ier)

! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!* at each convergence step. (consistent definition with the definition of hartree energy)
 if (fockcommon%ieigen/=0) then
   eigen=zero
!* Dot product of cwavef and ghc
!* inspired from the routine 54_spacepar/meanvalue_g but without the reference to parallelism and filtering
   if(gs_ham%istwf_k==2) then
     eigen=half*cwavef(1,1)*ghc1(1,1)
   else
     eigen=cwavef(1,1)*ghc1(1,1)+cwavef(2,1)*ghc1(2,1)
   end if
   do ipw=2,npw
     eigen=eigen+cwavef(1,ipw)*ghc1(1,ipw)+cwavef(2,ipw)*ghc1(2,ipw)
   end do
   if(gs_ham%istwf_k>=2) eigen=two*eigen
!   call xmpi_sum(eigen,mpi_enreg%comm_kpt,ier)
   fockcommon%eigen_ikpt(fockcommon%ieigen)= eigen
   fockcommon%ieigen = 0
 end if

! ===============================
! === Deallocate local arrays ===
! ===============================

 ABI_DEALLOCATE(ghc1)

end subroutine fock_ACE_getghc
!!***

end module m_fock_getghc
!!***
