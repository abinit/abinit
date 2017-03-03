!{\src2tex{textfont=tt}}
!!****f* ABINIT/fock_getghc
!! NAME
!!  fock_getghc
!!
!! FUNCTION
!!  Compute the matrix elements <G|Vx|psi> of the Fock operator.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2017 ABINIT group (CMartins,FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
!!  The current version assumes that :
!!   * nspinor = 1
!!   * no "my_nspinor"
!!   * no restriction to the value of istwfk_bz (but must be tested in all case)
!!   * all the data for the occupied states (cgocc_bz) are the same as those for the current states (cg)
!!
!! PARENTS
!!      forstrnps,getghc
!!
!! CHILDREN
!!      bare_vqg,dotprod_g,fftpac,fourdp,fourwf,hartre,load_kprime_hamiltonian
!!      nonlop,pawdijhat,pawmknhat_psipsi,sphereboundary,strfock,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fock_getghc(cwavef,cwaveprj,ghc,gs_ham,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi
 use m_cgtools, only :dotprod_g
 use m_fock
 use m_hamiltonian, only : gs_hamiltonian_type,load_kprime_hamiltonian,K_H_KPRIME

 use m_pawcprj
 use m_pawdij, only : pawdijhat

 use defs_datatypes, only: pseudopotential_type
 use m_pawrhoij,     only : pawrhoij_type, pawrhoij_free, pawrhoij_alloc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fock_getghc'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_65_paw
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! Scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),target,intent(inout) :: gs_ham
! Arrays
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)
 real(dp),intent(inout) :: cwavef(2,gs_ham%npw_k)!,ghc(2,gs_ham%npw_k)
 real(dp),intent(inout) :: ghc(:,:)

!Local variables-------------------------------
! Scalars
 integer,parameter :: tim_fourwf0=0,tim_fourdp0=0,ndat1=1
 integer :: bdtot_jindex,choice,cplex_fock,cplex_dij,cpopt,i1,i2,i3,ia,iatom
 integer :: iband_cprj,ider,idir,ier,ind,ipert,ipw,ifft,itypat,izero,jband,jbg,jcg,jkg
 integer :: jkpt,my_jsppol,jstwfk,lmn2_size,mgfftf,mpw,n1,n2,n3,n4,n5,n6
 integer :: n1f,n2f,n3f,n4f,n5f,n6f,natom,nband_k,ndij,nfft,nfftf,nfftotf,nhat12_grdim,nnlout
 integer :: npw,npwj,nspden_fock,nspinor,paw_opt,signs,tim_nonlop
 logical :: need_ghc
 real(dp),parameter :: weight1=one
 real(dp) :: doti,efock,eigen,imcwf,imcwocc,imvloc,invucvol,recwf,recwocc,revloc,occ,wtk
 type(fock_type),pointer :: fock
! Arrays
 integer :: ngfft(18),ngfftf(18)
 integer,pointer :: gboundf(:,:),kg_occ(:,:),gbound_kp(:,:)
 real(dp) ::enlout_dum(1),dotr(3),fockstr(6),for1(3),qphon(3),qvec_j(3),tsec(2),gsc_dum(2,0),rhodum(2,1),rhodum0(0,1,1)
 real(dp), allocatable :: dummytab(:,:),dijhat(:,:,:),dijhat_tmp(:,:),dijhat_tmp1(:,:,:),ffnl_kp_dum(:,:,:,:)
 real(dp), allocatable :: gvnlc(:,:),ghc1(:,:),ghc2(:,:),grnhat12(:,:,:,:),grnhat_12(:,:,:,:,:),forikpt(:,:)
 real(dp), allocatable :: rho12(:,:,:),rhog_munu(:,:),rhor_munu(:,:),vlocpsi_r(:)
 real(dp), allocatable :: vfock(:),psilocal(:,:,:),vectin_dum(:,:),vqg(:),forout(:,:),strout(:,:)
 real(dp), allocatable,target ::cwavef_r(:,:,:,:)
 real(dp), ABI_CONTIGUOUS  pointer :: cwaveocc_r(:,:,:,:)
 type(pawcprj_type),pointer :: cwaveocc_prj(:,:)

integer :: optgr, optgr2,optstr,optstr2,dumint(0)
real(dp) :: dummy(0),nhat_dum(0,0),rprimd(3,3),for11(3),for12(3),eigen1
 real(dp), allocatable ::grnl(:),vxc(:,:), vtrial(:,:)
type(pseudopotential_type) :: psps
 type(pawrhoij_type),allocatable  ::  pawrhoij(:)
 logical :: testtrue
! *************************************************************************
!return
 call timab(1504,1,tsec)
 call timab(1505,1,tsec)

 ABI_CHECK(associated(gs_ham%fock),"fock must be associated!")
 fock => gs_ham%fock

 testtrue=.false.
 if(testtrue)then
   ABI_DATATYPE_ALLOCATE(pawrhoij,(fock%natom))
   call pawrhoij_alloc(pawrhoij,2,1,gs_ham%nspinor,1,gs_ham%typat,pawtab=fock%pawtab,use_rhoijp=1)
   do iatom=1,fock%natom
     pawrhoij(iatom)%nrhoijsel=1
     pawrhoij(iatom)%rhoijselect(1)=1
   end do
 end if
 ABI_CHECK(gs_ham%nspinor==1,"only allowed for nspinor=1!")
 ABI_CHECK(gs_ham%npw_k==gs_ham%npw_kp,"only allowed for npw_k=npw_kp (ground state)!")
 if (fock%usepaw==1) then
   ABI_CHECK((size(cwaveprj,1)==gs_ham%natom.and.size(cwaveprj,2)==gs_ham%nspinor),"error on cwaveprj dims")
 end if
 need_ghc=(size(ghc,2)>0)

!Some constants
 invucvol=1.d0/sqrt(gs_ham%ucvol)
 cplex_fock=2;nspden_fock=1
 natom=fock%natom
 nspinor=gs_ham%nspinor
 mpw=maxval(fock%npwarr)
 npw=gs_ham%npw_k
 ider=0;izero=0
 if (fock%usepaw==1) then
   nfft =fock%pawfgr%nfftc ; ngfft =fock%pawfgr%ngfftc
   nfftf=fock%pawfgr%nfft  ; ngfftf=fock%pawfgr%ngfft
   mgfftf=fock%pawfgr%mgfft
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
! fock%optfor=.false.
 fock%optstr=.false.
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
 if (fock%usepaw==1) then
   nhat12_grdim=0
   if (((fock%optfor).and.(fock%ieigen/=0)).or.fock%optstr) ider=3
   ABI_ALLOCATE(grnhat_12,(2,nfftf,nspinor**2,3,natom*(ider/3)))
   ABI_ALLOCATE(gvnlc,(2,npw*nspinor))
   ABI_ALLOCATE(grnhat12,(2,nfftf,nspinor**2,3*nhat12_grdim))
   if ((fock%optfor).and.(fock%ieigen/=0)) then
     ABI_ALLOCATE(forout,(2,npw*nspinor))
     ABI_ALLOCATE(forikpt,(3,natom))
     forikpt=zero
   end if
   if (fock%optstr) then
     ABI_ALLOCATE(strout,(2,npw*nspinor))
   end if
 end if 
! fock%optstr=.true.
 if (fock%usepaw==1.or.fock%optstr) then
   ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
   call sphereboundary(gboundf,gs_ham%istwf_k,gs_ham%kg_k,mgfftf,npw)
 else
   gboundf=>gs_ham%gbound_k
 end if

! ==========================================
! === Get cwavef in real space using FFT ===
! ==========================================
 cwavef_r=zero

 call fourwf(0,rhodum0,cwavef,rhodum,cwavef_r,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
& mgfftf,mpi_enreg,ndat1,ngfftf,npw,1,n4f,n5f,n6f,0,mpi_enreg%paral_kgb,tim_fourwf0,weight1,weight1,&
& use_gpu_cuda=gs_ham%use_gpu_cuda)
 cwavef_r=cwavef_r*invucvol

! =====================================================
! === Select the states in cgocc_bz with the same spin ===
! =====================================================
!* Initialization of the indices/shifts, according to the value of isppol
!* bdtot_jindex = shift to be applied on the location of data in the array occ_bz ?
 bdtot_jindex=0
!* jbg = shift to be applied on the location of data in the array cprj/occ
 jbg=0;jcg=0
 my_jsppol=fock%isppol
 if((fock%isppol==2).and.(mpi_enreg%nproc_kpt/=1)) my_jsppol=1

 call timab(1505,2,tsec)
 call timab(1506,1,tsec)

 if (fock%optstr) then
   efock=zero
 end if




!===================================
!=== Loop on the k-points in IBZ ===
!===================================
 jkg=0

 if (associated(gs_ham%ph3d_kp)) then
   nullify (gs_ham%ph3d_kp)
 end if

 do jkpt=1,fock%mkpt
!* nband_k = number of bands at point k_j
   nband_k=fock%nbandocc_bz(jkpt,my_jsppol)
!* wtk = weight in BZ of this k point
   wtk=fock%wtk_bz(jkpt) !*sqrt(gs_ham%ucvol)
!* jstwfk= how is stored the wavefunction
   jstwfk=fock%istwfk_bz(jkpt)
!* npwj= number of plane wave in basis for the wavefunction
   npwj=fock%npwarr(jkpt)
!* Basis sphere of G vectors
   if (allocated(fock%cgocc)) then
     gbound_kp => fock%gbound_bz(:,:,jkpt)
     kg_occ => fock%kg_bz(:,1+jkg:npwj+jkg)
   end if
!* Load k^prime hamiltonian in the gs_ham datastructure
!  Note: ffnl_kp / ph3d_kp / gbound_kp are not used

   if (associated(gs_ham%ph3d_kp)) then
     ABI_ALLOCATE(gs_ham%ph3d_kp,(2,npwj,gs_ham%matblk))
   end if
   call load_kprime_hamiltonian(gs_ham,kpt_kp=fock%kptns_bz(:,jkpt),&
&   istwf_kp=jstwfk,npw_kp=npwj,kg_kp=fock%kg_bz(:,1+jkg:npwj+jkg))
!* Some temporary allocations needed for PAW
   if (fock%usepaw==1) then
     ABI_ALLOCATE(vectin_dum,(2,npwj*nspinor))
     ABI_ALLOCATE(ffnl_kp_dum,(npwj,0,gs_ham%lmnmax,gs_ham%ntypat))
     call load_kprime_hamiltonian(gs_ham,ffnl_kp=ffnl_kp_dum)
   end if

! ======================================
! === Calculate the vector q=k_i-k_j ===
! ======================================
!* Evaluation of kpoint_j, the considered k-point in reduced coordinates
!     kpoint_j(:)=fock%kptns_bz(:,jkpt)
!* the vector qvec is expressed in reduced coordinates.
!     qvec(:)=kpoint_i(:)-kpoint_j(:)
   qvec_j(:)=gs_ham%kpt_k(:)-fock%kptns_bz(:,jkpt)
   call bare_vqg(qvec_j,fock%gsqcut,gs_ham%gmet,fock%usepaw,fock%hybrid_mixing,&
&   fock%hybrid_mixing_sr,fock%hybrid_range,nfftf,fock%nkpt_bz,ngfftf,gs_ham%ucvol,vqg)


! =================================================
! === Loop on the band indices jband of cgocc_k ===
! =================================================
   do jband=1,nband_k

! ==============================================
! === Get cwaveocc_r in real space using FFT ===
! ==============================================
     if (allocated(fock%cwaveocc_bz)) then
       cwaveocc_r => fock%cwaveocc_bz(:,:,:,:,jband+jbg,my_jsppol)
     else
       ABI_ALLOCATE(cwaveocc_r,(2,n4f,n5f,n6f))
       cwaveocc_r=zero
       call fourwf(0,rhodum0,fock%cgocc(:,1+jcg+npwj*(jband-1):jcg+jband*npwj,my_jsppol),rhodum,cwaveocc_r, &
&       gbound_kp,gbound_kp,gs_ham%istwf_k,kg_occ,kg_occ,mgfftf,mpi_enreg,ndat1,ngfftf,&
&       npwj,1,n4f,n5f,n6f,tim_fourwf0,mpi_enreg%paral_kgb,0,weight1,weight1,use_gpu_cuda=gs_ham%use_gpu_cuda)
       cwaveocc_r=cwaveocc_r*invucvol
     end if

!*   occ = occupancy of jband at this k point
     occ=fock%occ_bz(jband+bdtot_jindex,my_jsppol)

! If((jkpt/=fock%ikpt).or.(jband/=fock%ieigen)) cycle
! ================================================
! === Get the overlap density matrix rhor_munu ===
! ================================================
!* Calculate the overlap density matrix in real space = conj(cwaveocc_r)*cwavef_r
!* rhor_munu will contain the overlap density matrix.
! vfock=-int{conj(cwaveocc_r)*cwavef_r*dr'/|r-r'|}
     call timab(1508,1,tsec)
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
     if (fock%usepaw==1) then
       ABI_ALLOCATE(rho12,(2,nfftf,nspinor**2))
       iband_cprj=(my_jsppol-1)*fock%mkptband+jbg+jband
       cwaveocc_prj=>fock%cwaveocc_prj(:,iband_cprj:iband_cprj+nspinor-1)

       if(testtrue)then
         do iatom=1,fock%natom
           pawrhoij(iatom)%rhoijp (1,1)=(cwaveprj(iatom,1)%cp(1,1)*cwaveocc_prj(iatom,1)%cp(1,1)+&
&                                  cwaveprj(iatom,1)%cp(2,1)*cwaveocc_prj(iatom,1)%cp(2,1))
!          pawrhoij(iatom)%rhoijp (2,1)=(-cwaveprj(iatom,1)%cp(1,1)*cwaveocc_prj(iatom,1)%cp(2,1)+&
!&                                  cwaveprj(iatom,1)%cp(2,1)*cwaveocc_prj(iatom,1)%cp(1,1))
           pawrhoij(iatom)%rhoijp (2,1)=zero
         end do
       end if

       call pawmknhat_psipsi(cwaveprj,cwaveocc_prj,ider,izero,natom,natom,nfftf,ngfftf,&
&       nhat12_grdim,nspinor,fock%ntypat,fock%pawang,fock%pawfgrtab,grnhat12,rho12,&
&       fock%pawtab,gprimd=gs_ham%gprimd,grnhat_12=grnhat_12,qphon=qvec_j,xred=gs_ham%xred,atindx=gs_ham%atindx)
       rhor_munu(1,:)=rhor_munu(1,:)+rho12(1,:,nspinor)
       rhor_munu(2,:)=rhor_munu(2,:)-rho12(2,:,nspinor)
       ABI_DEALLOCATE(rho12)
     end if

    !Perform an FFT using fourwf to get rhog_munu = FFT^-1(rhor_munu)
     call fourdp(cplex_fock,rhog_munu,rhor_munu,-1,mpi_enreg,nfftf,ngfftf,mpi_enreg%paral_kgb,tim_fourdp0)
     call timab(1509,2,tsec)
     if(fock%optstr.and.(fock%ieigen/=0)) then
!     if(fock%optstr) then
       call strfock(efock,gs_ham%gprimd,fock%gsqcut,fockstr,fock%hybrid_mixing,fock%hybrid_mixing_sr,&
&       fock%hybrid_range,mpi_enreg,nfftf,ngfftf,fock%nkpt_bz,rhog_munu,gs_ham%ucvol,qvec_j)
       fock%stress_ikpt(:,fock%ieigen)=fock%stress_ikpt(:,fock%ieigen)+fockstr(:)*occ*wtk
!write(80,*) "fock-getghc",occ*wtk,fock%ieigen
       if (fock%usepaw==0.and.(.not.need_ghc)) cycle
     end if

! ===================================================
! === Calculate the local potential vfockloc_munu ===
! ===================================================
!* Apply the Poisson solver to "rhog_munu" while taking into account the effect of the vector "qvec"
!* This is precisely what is done in the subroutine hartre, with option cplex=2.
!* vfock will contain the local Fock potential, the result of hartre routine.
!* vfock = FFT( rhog_munu/|g+qvec|^2 )
     call timab(1510,1,tsec)
#if 0
!    rhog_munu=rhog_munu*fock%wtk_bz(jkpt)

     call hartre(cplex_fock,gs_ham%gmet,fock%gsqcut,fock%usepaw,mpi_enreg,nfftf,ngfftf,&
&     mpi_enreg%paral_kgb,qvec_j,rhog_munu,vfock,divgq0=fock%divgq0)

#else
     do ifft=1,nfftf
       rhog_munu(1,ifft) = rhog_munu(1,ifft) * vqg(ifft)
       rhog_munu(2,ifft) = rhog_munu(2,ifft) * vqg(ifft)
     end do


     call fourdp(cplex_fock,rhog_munu,vfock,+1,mpi_enreg,nfftf,ngfftf,mpi_enreg%paral_kgb,tim_fourdp0)

#endif
     call timab(1510,2,tsec)


!===============================================================
!======== Calculate Dij_Fock_hat contribution in case of PAW ===
!===============================================================

     if (fock%usepaw==1) then
       qphon=qvec_j;nfftotf=product(ngfftf(1:3))
       cplex_dij=cplex_fock;ndij=nspden_fock
       ABI_ALLOCATE(dijhat,(cplex_dij*gs_ham%dimekb1,natom,ndij))
       dijhat=zero
       do iatom=1,natom
         ipert=iatom
         itypat=gs_ham%typat(iatom)
         lmn2_size=fock%pawtab(itypat)%lmn2_size
         ABI_ALLOCATE(dijhat_tmp,(cplex_dij*lmn2_size,ndij))
         dijhat_tmp=zero
         call pawdijhat(cplex_fock,cplex_dij,dijhat_tmp,gs_ham%gprimd,iatom,ipert,&
&         natom,ndij,nfftf,nfftotf,nspden_fock,my_jsppol,fock%pawang,fock%pawfgrtab(iatom),&
&         fock%pawtab(itypat),vfock,qphon,gs_ham%ucvol,gs_ham%xred)
         dijhat(1:cplex_dij*lmn2_size,iatom,:)=dijhat_tmp(1:cplex_dij*lmn2_size,:)
         ABI_DEALLOCATE(dijhat_tmp)
       end do
       signs=2; cpopt=2;idir=0; paw_opt=1;nnlout=1;tim_nonlop=1
       if(need_ghc) then
         choice=1
         call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&         ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,gvnlc,enl=dijhat,&
&         select_k=K_H_KPRIME)
          ghc2=ghc2-gvnlc*occ*wtk
       end if

! Forces calculation
       if (fock%optfor.and.(fock%ieigen/=0)) then
         call matr3inv(gs_ham%gprimd,rprimd)
         choice=2; dotr=zero;doti=zero;cpopt=4

         if (testtrue) then
           optgr=1;optgr2=0;optstr=0;optstr2=0
           ABI_ALLOCATE(grnl,(3*natom))
           ABI_ALLOCATE(vxc,(nfftf,cplex_fock))
           ABI_ALLOCATE(vtrial,(nfftf,1))
           vxc=zero
           do ifft=1,nfftf
             vtrial(ifft,1)=vfock(2*ifft-1)
           end do
           grnl=zero
           call pawgrnl(gs_ham%atindx1,0,dummy,1,dummy,grnl,fock%gsqcut,mgfftf,natom,natom,&
&            gs_ham%nattyp,nfftf,ngfftf,nhat_dum,dummy,1,0,fock%ntypat,optgr,optgr2,optstr,optstr2,&
&            fock%pawang,fock%pawfgrtab,pawrhoij,fock%pawtab,gs_ham%ph1d,psps,qphon,rprimd,dumint,&
&            gs_ham%typat,vtrial,vxc,gs_ham%xred)
           ABI_DEALLOCATE(vxc)
           ABI_DEALLOCATE(vtrial)
           ABI_DEALLOCATE(grnl)
         end if

         do iatom=1,natom
           do idir=1,3
             call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&             ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,&
&             forout,enl=dijhat,iatom_only=iatom,&
&             select_k=K_H_KPRIME)
             call dotprod_g(dotr(idir),doti,gs_ham%istwf_k,npw,2,cwavef,forout,mpi_enreg%me_g0,mpi_enreg%comm_fft)
             for1(idir)=zero
             do ifft=1,fock%pawfgrtab(iatom)%nfgd
               ind=fock%pawfgrtab(iatom)%ifftsph(ifft)
               for1(idir)=for1(idir)+vfock(2*ind-1)*grnhat_12(1,ind,1,idir,iatom)-&
&               vfock(2*ind)*grnhat_12(2,ind,1,idir,iatom)
             end do
           end do

           do idir=1,3
             for12(idir)=rprimd(1,idir)*for1(1)+rprimd(2,idir)*for1(2)+rprimd(3,idir)*for1(3)
             forikpt(idir,iatom)=forikpt(idir,iatom)-(for12(idir)*gs_ham%ucvol/nfftf+dotr(idir))*occ*wtk
           end do
         end do

       end if
! Stresses calculation
       if (fock%optstr) then
         choice=3;cpopt=4
         do idir=1,6
           call nonlop(choice,cpopt,cwaveocc_prj,enlout_dum,gs_ham,idir,(/zero/),mpi_enreg,&
&           ndat1,nnlout,paw_opt,signs,gsc_dum,tim_nonlop,vectin_dum,strout,enl=dijhat,&
&           select_k=K_H_KPRIME)
           call dotprod_g(dotr(1),doti,gs_ham%istwf_k,npw,2,cwavef,strout,mpi_enreg%me_g0,mpi_enreg%comm_fft)
           fock%stress_ikpt(idir,1)=fock%stress_ikpt(idir,1)-dotr(1)*occ*wtk
         end do
       end if

       ABI_DEALLOCATE(dijhat)
!       if (.not.need_ghc) cycle
     end if

! =============================================================
! === Apply the local potential vfockloc_munu to cwaveocc_r ===
! =============================================================
     call timab(1507,1,tsec)
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
     if (allocated(fock%cgocc)) then
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
   if (fock%usepaw==1) then
     ABI_DEALLOCATE(vectin_dum)
     ABI_DEALLOCATE(ffnl_kp_dum)
   end if
   if (associated(gs_ham%ph3d_kp)) then
     ABI_DEALLOCATE(gs_ham%ph3d_kp)
   end if
 end do ! jkpt

 if (fock%usepaw==1) then
   if ((fock%optfor).and.(fock%ieigen/=0)) then
     call xmpi_sum(forikpt,mpi_enreg%comm_hf,ier)
!     call sygradfock(fock%forces_ikpt(:,:,fock%ieigen),natom,forikpt,fock%pawang%nsym,fock%symrec,fock%indsym)
     do iatom=1,natom !Loop over atom
       ia=gs_ham%atindx(iatom)
       fock%forces_ikpt(:,ia,fock%ieigen)=forikpt(:,iatom)
     end do
   end if
 end if
 if(fock%optstr) then
   call xmpi_sum(fock%stress_ikpt,mpi_enreg%comm_hf,ier)
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
   if (fock%usepaw==1) then
     ABI_DEALLOCATE(gvnlc)
     ABI_DEALLOCATE(grnhat12)
     if ((fock%optfor).and.(fock%ieigen/=0)) then
       ABI_DEALLOCATE(forikpt)
       ABI_DEALLOCATE(forout)
     end if
     if (fock%optstr) then
       ABI_DEALLOCATE(strout)
     end if
     ABI_DEALLOCATE(grnhat_12)
   end if
   if(fock%usepaw==1.or.fock%optstr) then
     ABI_DEALLOCATE(gboundf)
   end if
!*Restore gs_ham datastructure

   if (associated(gs_ham%ph3d_kp)) then
     ABI_ALLOCATE(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
   end if
   call load_kprime_hamiltonian(gs_ham,kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
&   npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!   if (fock%ieigen/=0) fock%ieigen=0
   return
 end if
 call timab(1506,2,tsec)
 call timab(1511,1,tsec)

!*Restore gs_ham datastructure

 if (associated(gs_ham%ph3d_kp)) then
   ABI_ALLOCATE(gs_ham%ph3d_kp,(2,gs_ham%npw_k,gs_ham%matblk))
 end if
 call load_kprime_hamiltonian(gs_ham,kpt_kp=gs_ham%kpt_k,istwf_kp=gs_ham%istwf_k,&
& npw_kp=gs_ham%npw_k,kg_kp=gs_ham%kg_k,ffnl_kp=gs_ham%ffnl_k,ph3d_kp=gs_ham%ph3d_k)

!* Perform an FFT using fourwf to get ghc1 = FFT^-1(vlocpsi_r)
 ABI_ALLOCATE(psilocal,(cplex_fock*n4f,n5f,n6f))
 call fftpac(1,mpi_enreg,nspden_fock,cplex_fock*n1f,n2f,n3f,cplex_fock*n4f,n5f,n6f,ngfft,vlocpsi_r,psilocal,2)

 call fourwf(0,rhodum0,rhodum,ghc1,psilocal,gboundf,gboundf,gs_ham%istwf_k,gs_ham%kg_k,gs_ham%kg_k,&
& mgfftf,mpi_enreg,ndat1,ngfftf,1,npw,n4f,n5f,n6f,3,mpi_enreg%paral_kgb,tim_fourwf0,weight1,weight1,&
& use_gpu_cuda=gs_ham%use_gpu_cuda)
 ABI_DEALLOCATE(psilocal)

 ghc1=ghc1*sqrt(gs_ham%ucvol)+ghc2

!* If the calculation is parallelized, perform an MPI_allreduce to sum all the contributions in the array ghc
 ghc(:,:)=ghc(:,:)/mpi_enreg%nproc_hf + ghc1(:,:)



 call xmpi_sum(ghc,mpi_enreg%comm_hf,ier)

 call timab(1511,2,tsec)


! ===============================
! === Deallocate local PAW arrays ===
! ===============================

 if (fock%usepaw==1) then
   ABI_DEALLOCATE(gvnlc)
   ABI_DEALLOCATE(grnhat12)
   if ((fock%optfor).and.(fock%ieigen/=0)) then
     ABI_DEALLOCATE(forikpt)
     ABI_DEALLOCATE(forout)
   end if
   if (fock%optstr) then
     ABI_DEALLOCATE(strout)
   end if
   ABI_DEALLOCATE(grnhat_12)
 end if
 if(fock%usepaw==1.or.fock%optstr) then
   ABI_DEALLOCATE(gboundf)
 end if


! ============================================
! === Calculate the contribution to energy ===
! ============================================
!* Only the contribution when cwavef=cgocc_bz are calculated, in order to cancel exactly the self-interaction
!* at each convergence step. (consistent definition with the defintion of hartree energy)
 if (fock%ieigen/=0) then
   eigen=zero
!* Dot product of cwavef and ghc
!* inspired from the routine 53_spacepar/meanvalue_g but without the reference to parallelism and filtering
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
   fock%eigen_ikpt(fock%ieigen)= eigen
   fock%ieigen = 0
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
 if(testtrue)then
   call pawrhoij_free(pawrhoij)
   ABI_DATATYPE_DEALLOCATE(pawrhoij)
 end if
 call timab(1504,2,tsec)

end subroutine fock_getghc
!!***
