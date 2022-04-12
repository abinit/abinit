!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2022 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

#define MATPACK(row,col) (MAX(row,col)*(MAX(row,col)-1)/2 + MIN(row,col))

module m_orbmag

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_cgprj,            only : getcprj
  use m_fft,              only : fftpac
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type
  use m_kg,               only : getph,mkkin,mkkpg,ph1d3d
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle
  use m_nonlop,           only : nonlop
  use m_pawfgr,           only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen
  use m_paw_sphharm,      only : setsym_ylm,slxyzs,realgaunt
  use m_pawtab,           only : pawtab_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get, pawcprj_put
  use m_spacepar,         only : make_vectornd
  use m_time,             only : timab

  implicit none

  ! antisymmetric unit tensor, for doing the crossproduct summations
  real(dp),parameter :: eijk(3,3,3) = reshape((/zero,zero,zero,& !{1..3}11
                                               &zero,zero,-one,& !{1..3}21
                                               &zero,one,zero,& !{1..3}31
                                               &zero,zero,one,& !{1..3}12
                                               &zero,zero,zero,& !{1..3}22
                                               &-one,zero,zero,& !{1..3}32
                                               &zero,-one,zero,& !{1..3}13
                                               &one,zero,zero,& !{1..3}23
                                               &zero,zero,zero/),& !{1..3}33
                                               &(/3,3,3/))

  ! Bound methods:

  public :: orbmag_tt

  private :: gs_eigenvalues
  private :: d2lr_p2
  private :: d2lr_Anp
  private :: dl_Anp
  private :: dl_p2
  private :: dl_q
  private :: dl_vhnzc
  private :: dl_vha
  private :: apply_pw_k
  private :: apply_onsite_d0_k
  private :: apply_d2lr_term_k
  private :: apply_dl_term_k
  private :: lamb_core
  private :: orbmag_tt_output
  private :: local_rhoij
  
CONTAINS  !========================================================================================
!!***

!!****f* ABINIT/orbmag_tt
!! NAME
!! orbmag_tt
!!
!! FUNCTION
!! This routine computes the orbital magnetization and Berry curvature based on input 
!! wavefunctions and DDK wavefuntions. It is based on using the modern theory of
!! orbital magnetism directly.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2022 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006) [[cite:Ceresoli2006]],
!! DDK wavefunctions are used for the derivatives.
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine orbmag_tt(cg,cg1,cprj,dtset,gsqcut,kg,mcg,mcg1,mcprj,mpi_enreg,&
    & nfftf,ngfftf,npwarr,occ,paw_ij,pawfgr,pawrad,pawtab,psps,rprimd,vtrial,&
    & xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcprj,mcg,mcg1,nfftf
 real(dp),intent(in) :: gsqcut
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps

 !arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),ngfftf(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3),rprimd(3,3),xred(3,dtset%natom)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,buff_size,choice,cpopt,dimffnl,exchn2n3d,iat,iatom,icg,icmplx,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat,lmn2max
 integer :: klm
 integer :: me,mcgk,my_lmax,my_nspinor,nband_k,ncpgr,ndat,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 integer,parameter :: ibcpw=1,ibcdpdp=2,ibcdpu=3,ibcdudu=4,ibcd=5
 integer,parameter :: iompw=6,iomdpdp=7,iomdpu=8,iomdudu=9,iomlr=10,iomlmb=11,iomanp=12
 integer,parameter :: iomdlanp=13,iomdlp2=14,iomdlvhnzc=15,iomdlvh=16,iomdlq=17
 integer,parameter :: nterms=17
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 complex(dpc) :: cberry,com
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bc_k(:,:,:),bco_k(:,:,:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: Enk(:),ffnl_k(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: om_k(:,:,:),omo_k(:,:,:,:),orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:),rhoij(:,:,:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 complex(dpc),allocatable :: dlanp(:,:,:),dlp2(:,:,:),dlq(:,:,:),dlvh(:,:,:),dlvhnzc(:,:,:)
 complex(dpc),allocatable :: mp2(:,:,:),mpan(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj1_k(:,:,:),cwaveprj(:,:)

 !----------------------------------------------

 ! set up basic FFT parameters
 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 nband_k = dtset%mband
 istwf_k = 1
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me = mpi_enreg%me_kpt
 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0; ikg1 = 0

 !Definition of atindx array
 !Generate an index table of atoms, in order for them to be used type after type.
 ABI_MALLOC(atindx,(dtset%natom))
 ABI_MALLOC(atindx1,(dtset%natom))
 ABI_MALLOC(nattyp,(psps%ntypat))
 indx=1
 do itypat=1,psps%ntypat
   nattyp(itypat)=0
   do iatom=1,dtset%natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do
 ABI_MALLOC(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
 call getph(atindx,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

 ABI_MALLOC(kg_k,(3,dtset%mpw))
 ABI_MALLOC(kinpw,(dtset%mpw))

 ncpgr = 3
 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
 ABI_MALLOC(cprj_k,(dtset%natom,dtset%mband))
 call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
 ABI_MALLOC(cprj1_k,(dtset%natom,dtset%mband,3))
 do adir = 1, 3
   call pawcprj_alloc(cprj1_k(:,:,adir),ncpgr,dimlmn)
 end do
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
      & paw_ij=paw_ij)

 !========= construct local potential ==================
 ! nspden=1 is essentially hard-coded in the following line
 ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
 call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
 ABI_MALLOC(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
 call fftpac(isppol,mpi_enreg,dtset%nspden,&
      & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)
 ABI_FREE(cgrvtrial)
 call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)
 
 !========  compute nuclear dipole vector potential (may be zero) ==========
 with_vectornd=0
 has_nucdip = ANY( ABS(dtset%nucdipmom) .GT. tol8 )
 if (has_nucdip) with_vectornd=1
 ABI_MALLOC(vectornd,(with_vectornd*nfftf,3))
 vectornd = zero
 if(has_nucdip) then
   call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,&
     & dtset%nucdipmom,rprimd,vectornd,xred)
   ABI_MALLOC(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
   ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
   do idir = 1, 3
     call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
     call fftpac(isppol,mpi_enreg,dtset%nspden,&
       & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
   end do
   ABI_FREE(cgrvtrial)
   call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
 end if

 !========  compute onsite terms ==========

 lmn2max = psps%lmnmax*(psps%lmnmax+1)/2

 ABI_MALLOC(rhoij,(2,lmn2max,dtset%natom))
 write(std_out,'(a)')' Calling local_rhoij...'
 call local_rhoij(atindx,dtset,cprj,lmn2max,mcprj,mpi_enreg,rhoij,pawtab,ucvol)

 ABI_MALLOC(dlq,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling dl_q...'
 call dl_q(dlq,dtset,gprimd,lmn2max,pawtab)

 ABI_MALLOC(mp2,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling d2lr_p2...'
 call d2lr_p2(dtset,gprimd,lmn2max,mp2,pawrad,pawtab)

 ABI_MALLOC(mpan,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling d2lr_Anp...'
 call d2lr_Anp(dtset,gntselect,gprimd,lmn2max,mpan,my_lmax,pawrad,pawtab,realgnt)

 ABI_MALLOC(dlanp,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling dl_Anp...'
 call dl_Anp(dlanp,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)
 
 ABI_MALLOC(dlvhnzc,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling dl_vhnzc...'
 call dl_vhnzc(dlvhnzc,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)
 
 ABI_MALLOC(dlvh,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling dl_vha...'
 call dl_vha(atindx,dlvh,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,psps,realgnt,rhoij)
 
 ABI_MALLOC(dlp2,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling dl_p2...'
 call dl_p2(dlp2,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)

 icg = 0
 ikg = 0
 icprj = 0

 ABI_MALLOC(orbmag_terms,(2,nband_k,3,nterms))
 orbmag_terms = zero 

 write(std_out,'(a)')' entering kpt loop...'

 !============= BIG FAT KPT LOOP :) ===========================
 do ikpt = 1, dtset%nkpt

   ! if the current kpt is not on the current processor, cycle
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

   ! trace norm: assume occupation of two for each band and weight by kpts.
   ! division by ucvol arises from integration over BZ (so multiplication by 
   ! \Omega_{BZ} or division by \Omega
   trnrm = two*dtset%wtk(ikpt)/ucvol

   kpoint(:)=dtset%kptns(:,ikpt)
   npw_k = npwarr(ikpt)

   ! retrieve kg_k at this k point
   kg_k(1:3,1:npw_k) = kg(1:3,ikg+1:ikg+npw_k)

   ! retrieve ylm at this k point
   ABI_MALLOC(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
   ABI_MALLOC(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
   do ilm=1,psps%mpsang*psps%mpsang
     ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
   end do

   ! Compute kinetic energy at kpt
   kinpw(:) = zero
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

   ! Compute k+G at this k point
   nkpg = 3
   ABI_MALLOC(kpg_k,(npw_k,nkpg))
   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

   ! Make 3d phase factors
   ABI_MALLOC(phkxred,(2,dtset%natom))
   do iat = 1, dtset%natom
     iatom = atindx(iat)
     arg=two_pi*DOT_PRODUCT(kpoint,xred(:,iat))
     phkxred(1,iatom)=DCOS(arg);phkxred(2,iatom)=DSIN(arg)
   end do
   ABI_MALLOC(ph3d,(2,npw_k,dtset%natom))
   call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,&
     & npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)

   ! Compute nonlocal form factors ffnl at all (k+G):
   ider=1 ! ffnl and 1st derivatives
   idir=4 ! ignored when ider = 0; idir=0 means d ffnl/ dk in reduced units referenced 
          ! to reciprocal translations
          ! idir=4 meand d ffnl / dk in reduced units referenced to real space
          ! translations. rfddk = 1 wavefunctions are computed using this convention.
   dimffnl=4 ! 1 + number of derivatives
   ABI_MALLOC(ffnl_k,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,&
     & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
     & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
     & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
     & psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
   !  - Load k-dependent quantities in the Hamiltonian
   call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
     & kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,&
     & compute_gbound=.TRUE.)

   ! retrieve ground state wavefunctions at this k point
   mcgk = npw_k*nband_k
   ABI_MALLOC(cg_k,(2,mcgk))
   cg_k = cg(1:2,icg+1:icg+mcgk)

   ! retrieve first order wavefunctions at this k point
   ABI_MALLOC(cg1_k,(2,mcgk,3))
   cg1_k = cg1(1:2,icg+1:icg+mcgk,1:3)

   ! retrieve cprj_k
   call pawcprj_get(atindx,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
     & dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

   ! compute <p|cg1> cprjs
   ABI_MALLOC(cwavef,(2,npw_k))
   choice = 5
   cpopt = 0
   idir = 0
   do nn = 1, nband_k
     cwavef(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)
     call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl_k,idir,psps%indlmn,istwf_k,&
       & kg_k,kpg_k,kpoint,psps%lmnmax,dtset%mgfft,mpi_enreg,dtset%natom,nattyp,dtset%ngfft,&
       & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)
     call pawcprj_put(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
       & dtset%mkmem,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)
     do adir = 1, 3
       cwavef(1:2,1:npw_k) = cg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir)
       call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl_k,idir,psps%indlmn,istwf_k,&
         & kg_k,kpg_k,kpoint,psps%lmnmax,dtset%mgfft,mpi_enreg,dtset%natom,nattyp,dtset%ngfft,&
         & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)
       call pawcprj_put(atindx,cwaveprj,cprj1_k(:,:,adir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
         & dtset%mkmem,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)
     end do
   end do
   ABI_FREE(cwavef)

   ! compute GS eigenvalues
   ABI_MALLOC(Enk,(nband_k))
   call gs_eigenvalues(atindx,cg_k,cprj_k,dimlmn,dtset,Enk,gs_hamk,ikpt,mcgk,mpi_enreg,nband_k,npw_k)

   !--------------------------------------------------------------------------------
   ! Finally ready to compute contributions to orbital magnetism and Berry curvature
   !--------------------------------------------------------------------------------

   ! overall pre-factors
   cberry = j_dpc/two_pi
   com = -half*j_dpc
   
   !--------------------------------------------------------------------------------
   ! planewave contributions to Berry curvature and orbital magnetism
   !--------------------------------------------------------------------------------
   ABI_MALLOC(bc_k,(2,nband_k,3))
   ABI_MALLOC(bco_k,(2,nband_k,3,3))
   ABI_MALLOC(om_k,(2,nband_k,3))
   ABI_MALLOC(omo_k,(2,nband_k,3,3))

   !call bc_pw_k(cg1_k,cberry,dimlmn,gs_hamk,mcgk,mpi_enreg,dtset%natom,nband_k,npw_k,dtset%nspinor,bcpw_k)
   call apply_pw_k(bc_k,cberry,cg1_k,dimlmn,Enk,gs_hamk,mcgk,mpi_enreg,dtset%natom,&
     & nband_k,npw_k,dtset%nspinor,om_k,com)

   orbmag_terms(:,:,:,ibcpw) = orbmag_terms(:,:,:,ibcpw) + trnrm*bc_k
   orbmag_terms(:,:,:,iompw) = orbmag_terms(:,:,:,iompw) + trnrm*om_k
   
   !--------------------------------------------------------------------------------
   ! additional onsite terms to both Berry and orbmag (no |phi> derivatives) 
   !--------------------------------------------------------------------------------

   call apply_onsite_d0_k(atindx,bco_k,cberry,cprj_k,cprj1_k,Enk,&
     & dtset%natom,nband_k,dtset%ntypat,omo_k,com,paw_ij,pawtab,dtset%typat)

   orbmag_terms(:,:,:,ibcdpdp) = orbmag_terms(:,:,:,ibcdpdp) + trnrm*bco_k(:,:,:,1)
   orbmag_terms(:,:,:,ibcdpu) = orbmag_terms(:,:,:,ibcdpu) + trnrm*bco_k(:,:,:,2)
   orbmag_terms(:,:,:,ibcdudu) = orbmag_terms(:,:,:,ibcdudu) + trnrm*bco_k(:,:,:,3)

   orbmag_terms(:,:,:,iomdpdp) = orbmag_terms(:,:,:,iomdpdp) + trnrm*omo_k(:,:,:,1)
   orbmag_terms(:,:,:,iomdpu) = orbmag_terms(:,:,:,iomdpu) + trnrm*omo_k(:,:,:,2)
   orbmag_terms(:,:,:,iomdudu) = orbmag_terms(:,:,:,iomdudu) + trnrm*omo_k(:,:,:,3)

   !--------------------------------------------------------------------------------
   ! apply dlq term to Berry curvature
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlq,bc_k,cberry,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat)
   orbmag_terms(:,:,:,ibcd) = orbmag_terms(:,:,:,ibcd) + trnrm*bc_k

   !--------------------------------------------------------------------------------
   ! apply dl_vh term to orb mag
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlvh,om_k,com,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat)
   orbmag_terms(:,:,:,iomdlvh) = orbmag_terms(:,:,:,iomdlvh) + trnrm*om_k

   !--------------------------------------------------------------------------------
   ! apply dl_vhnzc term to orb mag
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlvhnzc,om_k,com,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat)
   orbmag_terms(:,:,:,iomdlvhnzc) = orbmag_terms(:,:,:,iomdlvhnzc) + trnrm*om_k

   !--------------------------------------------------------------------------------
   ! apply dl_Anp term to orb mag
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlanp,om_k,com,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat)
   orbmag_terms(:,:,:,iomdlanp) = orbmag_terms(:,:,:,iomdlanp) + trnrm*om_k

   !--------------------------------------------------------------------------------
   ! apply dl_Eqij term to orb mag
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlq,om_k,com,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat,Enk=Enk)
   orbmag_terms(:,:,:,iomdlq) = orbmag_terms(:,:,:,iomdlq) + trnrm*om_k


   !--------------------------------------------------------------------------------
   ! apply dl_p2 term to orb mag
   !--------------------------------------------------------------------------------

   call apply_dl_term_k(atindx,cprj_k,cprj1_k,dlp2,om_k,com,lmn2max,dtset%natom,&
     & nband_k,dtset%ntypat,pawtab,dtset%typat)
   orbmag_terms(:,:,:,iomdlp2) = orbmag_terms(:,:,:,iomdlp2) + trnrm*om_k

   !--------------------------------------------------------------------------------
   ! onsite <phi|r_b p^2/2 r_g>
   !--------------------------------------------------------------------------------
   call apply_d2lr_term_k(atindx,com,cprj_k,dtset,lmn2max,mp2,nband_k,om_k,pawtab)
   orbmag_terms(:,:,:,iomlr) = orbmag_terms(:,:,:,iomlr) + trnrm*om_k
   
   !--------------------------------------------------------------------------------
   ! onsite <phi|r_b p.A0 r_g> 
   !--------------------------------------------------------------------------------
   call apply_d2lr_term_k(atindx,com,cprj_k,dtset,lmn2max,mpan,nband_k,om_k,pawtab)
   orbmag_terms(:,:,:,iomanp) = orbmag_terms(:,:,:,iomanp) + trnrm*om_k

   ABI_FREE(bc_k)
   ABI_FREE(bco_k)
   ABI_FREE(om_k)
   ABI_FREE(omo_k)

   icg = icg + npw_k*nband_k
   ikg = ikg + npw_k
   icprj = icprj + nband_k

   ABI_FREE(cg_k)
   ABI_FREE(cg1_k)
   ABI_FREE(Enk)
   ABI_FREE(ylm_k)
   ABI_FREE(ylmgr_k)
   ABI_FREE(kpg_k)
   ABI_FREE(ffnl_k)
   ABI_FREE(ph3d)
   ABI_FREE(phkxred)

 end do ! end loop over kpts

 !! collect orbmag_terms if distributed over different processes
 if (nproc > 1) then
   buff_size=size(orbmag_terms)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = reshape(orbmag_terms,(/2*nband_k*3*nterms/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   orbmag_terms(1:2,1:nband_k,1:3,1:nterms)=reshape(buffer2,(/2,nband_k,3,nterms/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

 !! convert to cartesian frame from reduced triclinic
 ! general results: [rprimd]*r_red = r_cart
 !                  [gprimd]*k_red = k_cart
 !                  [gprimd]*grad_(r_red) = grad_(r_cart)
 !                  [rprimd]*grad_(k_red) = grad_(k_cart)
 ! most terms in orbital magnetism look like
 ! \grad_(k_cart) x \grad_(k_cart) but are calculated in reduced k
 ! so [rprimd]*grad_(k_red) x [rprimd]*grad_(k_red) = det[rprimd]*[rprimd]^{-1,T} grad_(k_red) x grad_(k_red)
 ! 
 do iterm = 1, nterms
   if(iterm.EQ.iomlmb) cycle
   do nn = 1, nband_k
     do icmplx = 1, 2
       if((iterm.EQ.iomlr).OR.(iterm.EQ.iomanp)) then
         orbmag_terms(icmplx,nn,1:3,iterm) = MATMUL(rprimd,orbmag_terms(icmplx,nn,1:3,iterm))
       else
         orbmag_terms(icmplx,nn,1:3,iterm) = ucvol*MATMUL(gprimd,orbmag_terms(icmplx,nn,1:3,iterm))
       end if
     end do
   end do
 end do

 ! convert orbmag magnetization to orbital moment
 ! Berry curvature terms are ignored
 ! Lamb term ignored
 do iterm = 1, nterms
   if ((iterm.EQ.ibcpw).OR.&
     & (iterm.EQ.ibcdpdp).OR.&
     & (iterm.EQ.ibcdpu).OR.&
     & (iterm.EQ.ibcdudu).OR.&
     & (iterm.EQ.ibcd).OR.&
     & (iterm.EQ.iomlmb)) cycle
   orbmag_terms(:,:,:,iterm) = ucvol*orbmag_terms(:,:,:,iterm)
 end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace = orbmag_trace + orbmag_terms(:,nn,:,:)
 end do

 ! get the Lamb term
 call lamb_core(atindx,dtset,omlamb,pawtab)
 orbmag_trace(:,:,iomlmb) = omlamb

 call orbmag_tt_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)

!---------------------------------------------------
! deallocate memory
!---------------------------------------------------

 call gs_hamk%free()

 ABI_FREE(vlocal)
 ABI_FREE(vectornd)
 if(has_nucdip) then
   ABI_FREE(vectornd_pac)
 end if
 ABI_FREE(kg_k)
 ABI_FREE(kinpw)
 ABI_FREE(ph1d)

 ABI_FREE(realgnt)
 ABI_FREE(gntselect)

 ABI_FREE(atindx)
 ABI_FREE(atindx1)
 ABI_FREE(nattyp)

 ABI_FREE(dimlmn)
 call pawcprj_free(cprj_k)
 ABI_FREE(cprj_k)
 do adir = 1, 3
   call pawcprj_free(cprj1_k(:,:,adir))
 end do
 ABI_FREE(cprj1_k)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 ABI_FREE(rhoij)
 ABI_FREE(dlq)
 ABI_FREE(mp2)
 ABI_FREE(mpan)
 ABI_FREE(dlanp)
 ABI_FREE(dlvh)
 ABI_FREE(dlvhnzc)
 ABI_FREE(dlp2)

 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

end subroutine orbmag_tt
!!***

!!****f* ABINIT/apply_pw_k
!! NAME
!! om_pw_k
!!
!! FUNCTION
!! compute pure planewave part of orbital magnetization and Berry curvature
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!
!! SOURCE

subroutine apply_pw_k(bc_k,bcpre,cg1_k,dimlmn,Enk,gs_hamk,mcgk,mpi_enreg,&
    & natom,nband_k,npw_k,nspinor,om_k,ompre)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcgk,natom,nband_k,npw_k,nspinor
  complex(dpc),intent(in) :: bcpre,ompre
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: dimlmn(natom)
  real(dp),intent(in) :: cg1_k(2,mcgk,3),Enk(nband_k)
  real(dp),intent(out) :: bc_k(2,nband_k,3),om_k(2,nband_k,3)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,gdir,isppol,my_ndat,my_nspinor,my_timer
  integer :: nbeg,nend,nn,prtvol,sij_opt,type_calc
  real(dp) :: c2,doti,dotr,eabg,lambda
  complex(dpc) :: cme
  !arrays
  real(dp),allocatable :: bcc(:,:),bwavef(:,:),ghc(:,:),gsc(:,:)
  real(dp),allocatable :: gvnlc(:,:),gwavef(:,:),omc(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

  ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
  ! term arises to properly normalize the derivatives (there are two in 
  ! the following terms, either <du|S|du> or <du|dS|u>
  c2=1.0d0/(two_pi*two_pi)

  ABI_MALLOC(cwaveprj,(natom,1))
  call pawcprj_alloc(cwaveprj,0,dimlmn)
  ABI_MALLOC(bwavef,(2,npw_k))
  ABI_MALLOC(gwavef,(2,npw_k))
  ABI_MALLOC(ghc,(2,npw_k))
  ABI_MALLOC(gsc,(2,npw_k))
  ABI_MALLOC(gvnlc,(2,npw_k))
  ABI_MALLOC(omc,(2,npw_k))
  ABI_MALLOC(bcc,(2,npw_k))

  my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
  isppol = 1
  cpopt=-1 ! no cprj
  lambda=zero 
  sij_opt = 1 ! retain gsc
  type_calc = 0 ! use full H
  my_timer = 0
  my_ndat = 1
  prtvol=0

  bc_k = zero
  om_k = zero
  do nn = 1, nband_k
    nbeg = (nn-1)*npw_k+1
    nend = nn*npw_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < tol8) cycle

          bwavef(1:2,1:npw_k) = cg1_k(1:2,nbeg:nend,bdir)
          gwavef(1:2,1:npw_k) = cg1_k(1:2,nbeg:nend,gdir)

          ! compute H|g>, S|g>, Vnl|g>
          call getghc(cpopt,gwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,my_ndat,&
            &           prtvol,sij_opt,my_timer,type_calc)
     
          ! orbital mag needs <du/dk|H + E|du/dk> , no nonlocal parts here, do that elsewhere
          omc = (ghc - gvnlc) + Enk(nn)*gwavef
          
          !omc = ghc + Enk(nn)*gsc
       
          dotr = DOT_PRODUCT(bwavef(1,:),omc(1,:))+DOT_PRODUCT(bwavef(2,:),omc(2,:))
          doti = DOT_PRODUCT(bwavef(1,:),omc(2,:))-DOT_PRODUCT(bwavef(2,:),omc(1,:))
          cme = cmplx(dotr,doti,KIND=dpc)
          om_k(1,nn,adir) = om_k(1,nn,adir) + real(ompre*eabg*c2*cme)
          om_k(2,nn,adir) = om_k(2,nn,adir) + aimag(ompre*eabg*c2*cme)

          ! berry curvature needs <du/dk|S|du/dk>, no nonlocal parts here, do that elsewhere
          bcc = gwavef
          
          !bcc = gsc

          dotr = DOT_PRODUCT(bwavef(1,:),bcc(1,:))+DOT_PRODUCT(bwavef(2,:),bcc(2,:))
          doti = DOT_PRODUCT(bwavef(1,:),bcc(2,:))-DOT_PRODUCT(bwavef(2,:),bcc(1,:))
          cme = cmplx(dotr,doti,KIND=dpc)
          bc_k(1,nn,adir) = bc_k(1,nn,adir) + real(bcpre*eabg*c2*cme)
          bc_k(2,nn,adir) = bc_k(2,nn,adir) + aimag(bcpre*eabg*c2*cme)

        end do ! gdir
      end do ! bdir
    end do ! adir
  end do !loop over bands

  ABI_FREE(bwavef)
  ABI_FREE(gwavef)
  ABI_FREE(ghc)
  ABI_FREE(gsc)
  ABI_FREE(gvnlc)
  ABI_FREE(omc)
  ABI_FREE(bcc)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)
 
end subroutine apply_pw_k
!!***

!!****f* ABINIT/apply_dl_term_k
!! NAME
!! apply_dl_term_k
!!
!! FUNCTION
!! apply left derivative term arising from <d/db phi|H|phi> =
!! i<phi|r_b H|phi>. Assume given in crystal coords
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine apply_dl_term_k(atindx,cprj_k,cprj1_k,dlij,dl_k,dlpre,&
    & lmn2max,natom,nband_k,ntypat,pawtab,typat,&
    & Enk) ! Enk is an optional argument

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: lmn2max,natom,nband_k,ntypat
 complex(dpc),intent(in) :: dlpre

 !arrays
 integer,intent(in) :: atindx(natom),typat(natom)
 real(dp),optional,intent(in) :: Enk(nband_k)
 real(dp),intent(out) :: dl_k(2,nband_k,3)
 complex(dpc),intent(in) :: dlij(lmn2max,natom,3)
 type(pawcprj_type),intent(inout) :: cprj_k(natom,nband_k)
 type(pawcprj_type),intent(inout) :: cprj1_k(natom,nband_k,3)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
 real(dp) :: c2,eabg
 complex(dpc) :: cpi,cpj,dlme

 !arrays
 real(dp),allocatable :: ewt(:)
 complex(dpc) :: dup(2,2),udp(2,2)

 !-----------------------------------------------------------------------

 ABI_MALLOC(ewt,(nband_k))
 ewt = one
 if(present(Enk)) then
   ewt = Enk
 end if

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in the Chern number,
 ! one for each wavefunction derivative) 
 c2=1.0d0/(two_pi*two_pi)

 dl_k = czero
 do nn = 1, nband_k
   do iat=1,natom
     iatom=atindx(iat)
     itypat=typat(iat)
     do klmn=1,pawtab(itypat)%lmn2_size
       ilmn = pawtab(itypat)%indklmn(7,klmn)
       jlmn = pawtab(itypat)%indklmn(8,klmn)

       do adir = 1, 3
         dlme = czero
         do bdir = 1, 3
           do gdir = 1, 3
             eabg = eijk(adir,bdir,gdir)
             if ( abs(eabg) < tol8 ) cycle

             cpi = cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn))
             cpj = cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn))
             ! <u|d_m p_n> : udp(m,n) means direction m, state m
             udp(1,1) = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
             udp(1,2) = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,jlmn),cprj_k(iatom,nn)%dcp(2,bdir,jlmn),KIND=dpc)
             udp(2,1) = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,ilmn),cprj_k(iatom,nn)%dcp(2,gdir,ilmn),KIND=dpc)     
             udp(2,2) = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)     

             ! <d_m u|p_n> : dup(m,n) means direction m, state m
             dup(1,1) = cmplx(cprj1_k(iatom,nn,bdir)%cp(1,ilmn),cprj1_k(iatom,nn,bdir)%cp(2,ilmn),KIND=dpc)
             dup(1,2) = cmplx(cprj1_k(iatom,nn,bdir)%cp(1,jlmn),cprj1_k(iatom,nn,bdir)%cp(2,jlmn),KIND=dpc)
             dup(2,1) = cmplx(cprj1_k(iatom,nn,gdir)%cp(1,ilmn),cprj1_k(iatom,nn,gdir)%cp(2,ilmn),KIND=dpc)     
             dup(2,2) = cmplx(cprj1_k(iatom,nn,gdir)%cp(1,jlmn),cprj1_k(iatom,nn,gdir)%cp(2,jlmn),KIND=dpc)     

             dlme = dlme + dlpre*c2*eabg*ewt(nn)*dlij(klmn,iat,bdir)*conjg(cpi)*(dup(2,2)+udp(2,2))
             if (ilmn /= jlmn) then
               dlme = dlme - dlpre*c2*eabg*ewt(nn)*conjg(dlij(klmn,iat,bdir))*conjg(cpj)*(dup(2,1)+udp(2,1))
             end if

             ! add on the right derivative term
             dlme = dlme + dlpre*c2*eabg*ewt(nn)*(-dlij(klmn,iat,gdir))*conjg(dup(1,1)+udp(1,1))*cpj
             if (ilmn /= jlmn) then
               dlme = dlme - dlpre*c2*eabg*ewt(nn)*conjg(-dlij(klmn,iat,gdir))*conjg(dup(1,2)+udp(1,2))*cpi
             end if

           end do ! end loop over gdir
         end do ! end loop over bdir
         dl_k(1,nn,adir) = dl_k(1,nn,adir) + REAL(dlme)
         dl_k(2,nn,adir) = dl_k(2,nn,adir) + AIMAG(dlme)

       end do ! end loop over adir

     end do ! end loop over klmn
   end do ! end loop over atoms
   
 end do ! end loop over nn

 ABI_FREE(ewt)

end subroutine apply_dl_term_k
!!***

!!****f* ABINIT/apply_onsite_d0_k
!! NAME
!! apply_onsite_d0_k
!!
!! FUNCTION
!! Compute the Berry curvature due to onsite terms 
!! <du|p> and <u|dp> but no derivatives of |phi>
!! Uses the DDK wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine apply_onsite_d0_k(atindx,bc_k,bcpre,cprj_k,cprj1_k,Enk,&
    & natom,nband_k,ntypat,om_k,ompre,paw_ij,pawtab,typat)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: natom,nband_k,ntypat
 complex(dpc),intent(in) :: bcpre,ompre

 !arrays
 integer,intent(in) :: atindx(natom),typat(natom)
 real(dp),intent(in) :: Enk(nband_k)
 real(dp),intent(out) :: bc_k(2,nband_k,3,3),om_k(2,nband_k,3,3)
 type(pawcprj_type),intent(inout) :: cprj_k(natom,nband_k)
 type(pawcprj_type),intent(inout) :: cprj1_k(natom,nband_k,3)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,gdir,iat,iatom,ilmn,iterm,itypat,jlmn,klmn,nn
 integer,parameter :: idpdp=1,idpu=2,idudu=3
 real(dp) :: c2,eabg,qij
 complex(dpc) :: dij,uij,uji
 logical :: cplex_dij

 !arrays
 complex(dpc) :: bcme(3),dup(2,2),udp(2,2),omme(6)

 !-----------------------------------------------------------------------

 cplex_dij = (paw_ij(1)%cplex_dij .EQ. 2)

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in the Chern number,
 ! one for each wavefunction derivative) 
 c2=1.0d0/(two_pi*two_pi)

 bc_k = czero
 om_k = czero
 do nn = 1, nband_k
   do iat=1,natom
     iatom=atindx(iat)
     itypat=typat(iat)
     do klmn=1,pawtab(itypat)%lmn2_size
       ilmn = pawtab(itypat)%indklmn(7,klmn)
       jlmn = pawtab(itypat)%indklmn(8,klmn)

       do adir = 1, 3

         bcme = czero
         omme = czero
         do bdir = 1, 3
           do gdir = 1, 3
             eabg = eijk(adir,bdir,gdir)
             if ( abs(eabg) < tol8 ) cycle

             ! <u|d_m p_n> : udp(m,n) means direction m, state n
             udp(1,1) = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
             udp(1,2) = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,jlmn),cprj_k(iatom,nn)%dcp(2,bdir,jlmn),KIND=dpc)
             udp(2,1) = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,ilmn),cprj_k(iatom,nn)%dcp(2,gdir,ilmn),KIND=dpc)     
             udp(2,2) = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)     

             ! <d_m u|p_n> : dup(m,n) means direction m, state n
             dup(1,1) = cmplx(cprj1_k(iatom,nn,bdir)%cp(1,ilmn),cprj1_k(iatom,nn,bdir)%cp(2,ilmn),KIND=dpc)
             dup(1,2) = cmplx(cprj1_k(iatom,nn,bdir)%cp(1,jlmn),cprj1_k(iatom,nn,bdir)%cp(2,jlmn),KIND=dpc)
             dup(2,1) = cmplx(cprj1_k(iatom,nn,gdir)%cp(1,ilmn),cprj1_k(iatom,nn,gdir)%cp(2,ilmn),KIND=dpc)     
             dup(2,2) = cmplx(cprj1_k(iatom,nn,gdir)%cp(1,jlmn),cprj1_k(iatom,nn,gdir)%cp(2,jlmn),KIND=dpc)     

             qij = pawtab(itypat)%sij(klmn)

             if (cplex_dij) then
               dij = cmplx(paw_ij(iat)%dij(2*klmn-1,1),paw_ij(iat)%dij(2*klmn,1),KIND=dpc)
             else
               dij = cmplx(paw_ij(iat)%dij(klmn,1),zero,KIND=dpc)
             end if

             ! uij is d/d_beta (<u|p_i> * d/d_gamma (<p_j|u>) 
             !uij = conjg(udp(1,1)+dup(1,1))*(udp(2,2)+dup(2,2))
             !uji = conjg(udp(1,2)+dup(1,2))*(udp(2,1)+dup(2,1))

             !uij = conjg(udp(1,1))*udp(2,2)+conjg(udp(1,1))*dup(2,2)+conjg(dup(1,1))*udp(2,2)
             !uji = conjg(udp(1,2))*udp(2,1)+conjg(udp(1,2))*dup(2,1)+conjg(dup(1,2))*udp(2,1)

             !uij = conjg(dup(1,1))*dup(2,2)
             !uji = conjg(dup(1,2))*dup(2,1)

             !uij = conjg(udp(1,1))*udp(2,2)
             !uji = conjg(udp(1,2))*udp(2,1)

             !bcme = bcme + bcpre*c2*epsabg*qij*uij
             bcme(idpdp) = bcme(idpdp) + bcpre*c2*eabg*qij*conjg(udp(1,1))*udp(2,2)
             bcme(idpu) =  bcme(idpu)  + bcpre*c2*eabg*qij*(conjg(udp(1,1))*dup(2,2)+conjg(dup(1,1))*udp(2,2))
             bcme(idudu) = bcme(idudu) + bcpre*c2*eabg*qij*conjg(dup(1,1))*dup(2,2)

             omme(idpdp) = omme(idpdp) + ompre*c2*eabg*(dij+Enk(nn)*qij)*conjg(udp(1,1))*udp(2,2)
             omme(idpu) =  omme(idpu)  + ompre*c2*eabg*(dij+Enk(nn)*qij)*(conjg(udp(1,1))*dup(2,2)+conjg(dup(1,1))*udp(2,2))
             omme(idudu) = omme(idudu) + ompre*c2*eabg*(dij+Enk(nn)*qij)*conjg(dup(1,1))*dup(2,2)

             !omme = omme + ompre*c2*epsabg*(dij+Enk(nn)*qij)*uij

             if (ilmn /= jlmn) then
               bcme(idpdp) = bcme(idpdp) + bcpre*c2*eabg*qij*conjg(udp(1,2))*udp(2,1)
               bcme(idpu) =  bcme(idpu)  + bcpre*c2*eabg*qij*(conjg(udp(1,2))*dup(2,1)+conjg(dup(1,2))*udp(2,1))
               bcme(idudu) = bcme(idudu) + bcpre*c2*eabg*qij*conjg(dup(1,2))*dup(2,1)

               omme(idpdp) = omme(idpdp) + ompre*c2*eabg*conjg(dij+Enk(nn)*qij)*conjg(udp(1,2))*udp(2,1)
               omme(idpu) =  omme(idpu)  + ompre*c2*eabg*conjg(dij+Enk(nn)*qij)*(conjg(udp(1,2))*dup(2,1)+conjg(dup(1,2))*udp(2,1))
               omme(idudu) = omme(idudu) + ompre*c2*eabg*conjg(dij+Enk(nn)*qij)*conjg(dup(1,2))*dup(2,1)
               !bcme = bcme + bcpre*c2*epsabg*qij*uji
               !omme = omme + ompre*c2*epsabg*conjg(dij+Enk(nn)*qij)*uji
             end if
           end do ! gdir
         end do ! bdir

         do iterm = idpdp,idudu
           bc_k(1,nn,adir,iterm) = bc_k(1,nn,adir,iterm) + REAL(bcme(iterm))
           bc_k(2,nn,adir,iterm) = bc_k(2,nn,adir,iterm) + AIMAG(bcme(iterm))
           om_k(1,nn,adir,iterm) = om_k(1,nn,adir,iterm) + REAL(omme(iterm))
           om_k(2,nn,adir,iterm) = om_k(2,nn,adir,iterm) + AIMAG(omme(iterm))
         end do

       end do ! end loop over adir

     end do ! end loop over klmn
   end do ! end loop over atoms
   
 end do ! end loop over nn

end subroutine apply_onsite_d0_k
!!***

!!****f* ABINIT/gs_eigenvalues
!! NAME
!! gs_eigenvalues
!!
!! FUNCTION
!! compute GS eigenvalues from GS wavefunctions
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!
!! SOURCE

subroutine gs_eigenvalues(atindx,cg_k,cprj_k,dimlmn,dtset,Enk,&
    & gs_hamk,ikpt,mcgk,mpi_enreg,nband_k,npw_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,mcgk,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk)
  real(dp),intent(out) :: Enk(nband_k)
  type(pawcprj_type),intent(inout) :: cprj_k(dtset%natom,dtset%mband)

  !Local variables -------------------------
  !scalars
  integer :: cpopt,isppol,my_ndat,my_nspinor,my_timer,nn,prtvol,sij_opt,type_calc
  real(dp) :: lambda
  !arrays
  real(dp),allocatable :: cwavef(:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

  ABI_MALLOC(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)
  ABI_MALLOC(cwavef,(2,npw_k))
  ABI_MALLOC(ghc,(2,npw_k))
  ABI_MALLOC(gsc,(2,npw_k))
  ABI_MALLOC(gvnlc,(2,npw_k))

  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  isppol = 1
  cpopt=4 ! use cprj in memory, not using lambda
  lambda=zero 
  sij_opt = 0 ! no gsc saved
  type_calc = 0 ! use all of Hamiltonian: kinetic, local, nonlocal
  my_timer = 0
  my_ndat = 1
  prtvol=0

  do nn = 1, nband_k
    cwavef(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)

    ! compute Enk and S|u_nk>
    call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
         &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
    call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,my_ndat,&
       &           prtvol,sij_opt,my_timer,type_calc)
    Enk(nn) = DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))

  end do

  ABI_FREE(cwavef)
  ABI_FREE(ghc)
  ABI_FREE(gsc)
  ABI_FREE(gvnlc)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)
 
end subroutine gs_eigenvalues
!!***


!!****f* ABINIT/lamb_core
!! NAME
!! lamb_core
!!
!! FUNCTION
!! add core electron contribution to the orbital magnetic moment
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine lamb_core(atindx,dtset,omlamb,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: omlamb(2,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,itypat

!--------------------------------------------------------------------

  omlamb = zero 
  do adir = 1, 3
    do iat=1,dtset%natom
      iatom = atindx(iat)
      itypat = dtset%typat(iat)
      omlamb(1,adir) = omlamb(1,adir) - dtset%lambsig(itypat)*dtset%nucdipmom(adir,iat)
    end do ! end loop over atoms
  end do ! end loop over adir
 
end subroutine lamb_core
!!***

!!****f* ABINIT/dl_vh
!! NAME
!! dl_vh
!!
!! FUNCTION
!! Compute onsite r_b vhartree
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computed in Cart directions then transformed to xtal coords
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine dl_vha(atindx,dlvh,dtset,gntselect,gprimd,lmn2max,my_lmax,&
    & pawrad,pawtab,psps,realgnt,rhoij)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  real(dp),intent(in) :: rhoij(2,lmn2max,dtset%natom)
  complex(dpc),intent(out) :: dlvh(lmn2max,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bigl,biglm,bigm,gnt1,gnt2,gnt3,iat,iatom,il,ilm,ilmn,im,ir,itypat
  integer :: jlm,jlmn,klmn,kln
  integer :: ll,llmm,llmmjlm,lm1b,max_ij_size,mm,msz,olmn,olm,oln
  real(dp) :: cdij,cl,intg1,rg1,rg2,rg3,rr
  complex(dpc) :: cpre,rhomn

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:),fmnl(:),intg(:,:,:,:),tfmnl(:)
  complex(dpc) :: dij_cart(3),dij_red(3)

!--------------------------------------------------------------------

 ! dij prefactor
 cdij = sqrt(four_pi/three)
 cpre = j_dpc

 max_ij_size = psps%lnmax*(psps%lnmax+1)/2 
 ! radial integrals
 ABI_MALLOC(intg,(max_ij_size,max_ij_size,2*psps%mpsang,dtset%ntypat))
 do itypat = 1, dtset%ntypat
   msz = pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(msz))
   ABI_MALLOC(fmnl,(msz))
   ABI_MALLOC(tfmnl,(msz))
   do oln = 1, pawtab(itypat)%ij_size
     do bigl = 0, pawtab(itypat)%l_size - 1
       do ir=2,msz

         rr=pawrad(itypat)%rad(ir)
         ff(2:ir) = pawtab(itypat)%phiphj(2:ir,oln)*pawrad(itypat)%rad(2:ir)**bigl/rr**(bigl+1)
         ff(ir:msz) = pawtab(itypat)%phiphj(ir:msz,oln)*rr**bigl/pawrad(itypat)%rad(ir:msz)**(bigl+1)
         call pawrad_deducer0(ff,msz,pawrad(itypat))
         call simp_gen(intg1,ff,pawrad(itypat))
         fmnl(ir) = intg1

         ff(2:ir) = pawtab(itypat)%tphitphj(2:ir,oln)*pawrad(itypat)%rad(2:ir)**bigl/rr**(bigl+1)
         ff(ir:msz) = pawtab(itypat)%tphitphj(ir:msz,oln)*rr**bigl/pawrad(itypat)%rad(ir:msz)**(bigl+1)
         call pawrad_deducer0(ff,msz,pawrad(itypat))
         call simp_gen(intg1,ff,pawrad(itypat))
         tfmnl(ir) = intg1

       end do ! end loop over ir

       do kln = 1, pawtab(itypat)%ij_size
         ff(2:msz) = pawrad(itypat)%rad(2:msz)*pawtab(itypat)%phiphj(2:msz,kln)*fmnl(2:msz)
         ff(2:msz) = ff(2:msz) - pawrad(itypat)%rad(2:msz)*pawtab(itypat)%tphitphj(2:msz,kln)*tfmnl(2:msz)
         call pawrad_deducer0(ff,msz,pawrad(itypat))
         call simp_gen(intg1,ff,pawrad(itypat))
         intg(kln,oln,bigl+1,itypat) = intg1
       end do ! end loop over kln
     end do ! end loop over bigl
   end do ! end loop over oln
   ABI_FREE(ff)
   ABI_FREE(fmnl)
   ABI_FREE(tfmnl)
 end do ! end loop over itypat

 dlvh = czero
 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   do klmn = 1, pawtab(itypat)%lmn2_size
     kln = pawtab(itypat)%indklmn(2,klmn)
     ilmn = pawtab(itypat)%indklmn(7,klmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn)
     il = pawtab(itypat)%indlmn(1,ilmn)
     im = pawtab(itypat)%indlmn(2,ilmn)
     ilm = pawtab(itypat)%indlmn(4,ilmn)
     jlm = pawtab(itypat)%indlmn(4,jlmn)

     dij_cart = czero
     do adir = 1, 3
       lm1b = MATPACK(ilm,adir_to_sij(adir))
       do ll = abs(il-1),il+1
         do mm = -ll,ll
           llmm = ll*ll+ll+1+mm
           llmmjlm = MATPACK(llmm,jlm)
           gnt2 = gntselect(llmm,lm1b)
           rg2 = zero
           if (gnt2 /= 0) rg2 = realgnt(gnt2)
           do olmn = 1, pawtab(itypat)%lmn2_size
             olm = pawtab(itypat)%indklmn(1,olmn)
             oln = pawtab(itypat)%indklmn(2,olmn)
             rhomn = cmplx(rhoij(1,olmn,iatom),rhoij(2,olmn,iatom),KIND=dpc)
             do bigl = pawtab(itypat)%indklmn(3,olmn),pawtab(itypat)%indklmn(4,olmn)
               cl = four_pi/(two*bigl + one)
               do bigm = -bigl, bigl
                 biglm = bigl*bigl+bigl+1+bigm
                 gnt1 = gntselect(biglm,olm)
                 rg1 = zero
                 if(gnt1 /= 0) rg1 = realgnt(gnt1)
                 gnt3 = gntselect(biglm,llmmjlm)
                 rg3=zero
                 if(gnt3 /= 0) rg3 = realgnt(gnt3)
                 dij_cart(adir) = dij_cart(adir) + cpre*cdij*cl*rhomn*rg1*rg2*rg3*intg(kln,oln,bigl+1,itypat)
               end do ! end loop over bigm
             end do ! end loop over bigl
           end do ! end loop over olmn
         end do ! end loop over mm
       end do ! end loop over ll
     end do ! end loop over adir
       
     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

     dlvh(klmn,iat,1:3) = dij_red(1:3)
   end do ! end loop over klmn
 end do ! end loop over iat

 ABI_FREE(intg) 
end subroutine dl_vha


!!****f* ABINIT/dl_vhnzc
!! NAME
!! dl_vhnzc
!!
!! FUNCTION
!! Compute onsite r_b vhnzc
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computed in Cart directions then transformed to xtal coords
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine dl_vhnzc(dlvhnzc,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: dlvhnzc(lmn2max,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,iat,itypat,klm,klmn,kln,lm1b,mesh_size
  real(dp) :: cdij,intg
  complex(dpc) :: cpre

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

 ! dij prefactor
 cdij = sqrt(four_pi/three)
 cpre = j_dpc
 dlvhnzc=czero
 do itypat = 1, dtset%ntypat
   mesh_size=pawrad(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))

   do klmn=1, pawtab(itypat)%lmn2_size

     klm = pawtab(itypat)%indklmn(1,klmn)
     kln = pawtab(itypat)%indklmn(2,klmn)

     ff(2:mesh_size) = pawtab(itypat)%phiphj(2:mesh_size,kln)*pawtab(itypat)%vhnzc(2:mesh_size)
     ff(2:mesh_size) = ff(2:mesh_size) - pawtab(itypat)%tphitphj(2:mesh_size,kln)*pawtab(itypat)%vhtnzc(2:mesh_size)
     ff(2:mesh_size) = ff(2:mesh_size)*pawrad(itypat)%rad(2:mesh_size)
     call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
     call simp_gen(intg,ff,pawrad(itypat))

     dij_cart = czero
     do adir = 1, 3
       lm1b = adir_to_sij(adir)
       gint = gntselect(lm1b,klm)
       if (gint == 0) cycle
       dij_cart(adir) = cdij*realgnt(gint)*intg
     end do
     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
     do iat = 1, dtset%natom
       if(dtset%typat(iat) == itypat) then
         dlvhnzc(klmn,iat,1:3) = cpre*dij_red(1:3)
       end if
     end do ! end loop over atoms
   end do ! end loop over klmn
   ABI_FREE(ff)
 end do ! end loop over types
 
end subroutine dl_vhnzc

!!****f* ABINIT/dl_q
!! NAME
!! dl_q
!!
!! FUNCTION
!! Compute onsite r_b Qij
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computed in Cart directions from pawtab%qijl(dir,klmn)
!! then transformed to xtal coords
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine dl_q(dlq,dtset,gprimd,lmn2max,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max
  type(dataset_type),intent(in) :: dtset

  !arrays
  real(dp),intent(in) :: gprimd(3,3)
  complex(dpc),intent(out) :: dlq(lmn2max,dtset%natom,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,itypat,klmn
  real(dp) :: cdij
  complex(dpc) :: cpre

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)

!--------------------------------------------------------------------

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)
 cpre = j_dpc
 dlq=czero
 do itypat = 1, dtset%ntypat
   do klmn=1, pawtab(itypat)%lmn2_size
     do adir = 1, 3
       dij_cart(adir) = cdij*pawtab(itypat)%qijl(adir_to_sij(adir),klmn)
     end do
     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
     do iat = 1, dtset%natom
       if(dtset%typat(iat) == itypat) then
         dlq(klmn,iat,1:3) = cpre*dij_red(1:3)
       end if
     end do ! end loop over atoms
   end do ! end loop over klmn
 end do ! end loop over types
 
end subroutine dl_q

!!****f* ABINIT/dl_p2
!! NAME
!! dl_p2
!!
!! FUNCTION
!! Compute onsite r_b p^2/2
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine dl_p2(dlp2,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: dlp2(lmn2max,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,iat,itypat,ilmn,iln,jl,jlmn,jln,klm,klmn,lm1b,mesh_size,pmesh_size
  real(dp) :: c1m,intg

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:),uj(:),ujder(:),uj2der(:)
  complex(dpc) :: dij_cart(3)

!--------------------------------------------------------------------

  c1m = sqrt(four_pi/three)
  dlp2=czero
  do iat=1,dtset%natom
    itypat=dtset%typat(iat)
    mesh_size=pawrad(itypat)%mesh_size
    pmesh_size=pawtab(itypat)%partialwave_mesh_size
    ABI_MALLOC(ff,(mesh_size))
    ABI_MALLOC(uj,(mesh_size))
    ABI_MALLOC(ujder,(mesh_size))
    ABI_MALLOC(uj2der,(mesh_size))

    do klmn=1, pawtab(itypat)%lmn2_size

      klm = pawtab(itypat)%indklmn(1,klmn)
      ilmn = pawtab(itypat)%indklmn(7,klmn)
      jlmn = pawtab(itypat)%indklmn(8,klmn)
      iln = pawtab(itypat)%indlmn(5,ilmn)
      jln = pawtab(itypat)%indlmn(5,jlmn)
      jl = pawtab(itypat)%indlmn(1,jlmn)

      uj = zero
      uj(1:pmesh_size) = pawtab(itypat)%phi(1:pmesh_size,jln)
      call nderiv_gen(ujder,uj,pawrad(itypat),uj2der)
      ff(2:mesh_size) = pawrad(itypat)%rad(2:mesh_size)*pawtab(itypat)%phi(2:mesh_size,iln)*(uj2der(2:mesh_size)&
        & - jl*(jl+1)*pawtab(itypat)%phi(2:mesh_size,jln)/pawrad(itypat)%rad(2:mesh_size)**2)
      uj = zero
      uj(1:pmesh_size) = pawtab(itypat)%tphi(1:pmesh_size,jln)
      call nderiv_gen(ujder,uj,pawrad(itypat),uj2der)
      ff(2:mesh_size) = ff(2:mesh_size) - &
        & (pawrad(itypat)%rad(2:mesh_size)*pawtab(itypat)%tphi(2:mesh_size,iln)*(uj2der(2:mesh_size)&
        & - jl*(jl+1)*pawtab(itypat)%tphi(2:mesh_size,jln)/pawrad(itypat)%rad(2:mesh_size)**2))
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))

      dij_cart = czero 
      do adir = 1, 3
        lm1b = adir_to_sij(adir)
        gint = gntselect(lm1b,klm)
        if (gint == 0) cycle
        dij_cart(adir) = cmplx(-half*c1m*realgnt(gint)*intg,zero)
      end do ! end loop over adir

      dlp2(klmn,iat,1:3) = MATMUL(TRANSPOSE(gprimd),dij_cart)
    end do ! end loop over klmn
    ABI_FREE(ff)
    ABI_FREE(uj)
    ABI_FREE(ujder)
    ABI_FREE(uj2der)
  end do ! end loop over atoms
 
end subroutine dl_p2

!!****f* ABINIT/dl_Anp
!! NAME
!! dl_Anp
!!
!! FUNCTION
!! Compute onsite r_b Anp 
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! effectively a vector in Cart directions, similar to pawtab%qijl(dir,klmn)
!! but for each atom rather than each type of atom, due to L.m term
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine dl_Anp(dlanp,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: dlanp(lmn2max,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bigl,biglm,bigm,gint,iat,ilmn,il,ilm,ilml1b,itypat
  integer :: jl,jlmn,jm,klmn,kln,l1b,ldir,mesh_size
  real(dp) :: a2,intg
  complex(dpc) :: cme,orbl_me

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:)
  complex(dpc) :: dij_cart(3),dij_red(3)

!--------------------------------------------------------------------

  a2 = FineStructureConstant2

  dlanp=czero
  do iat=1,dtset%natom
    itypat=dtset%typat(iat)
    mesh_size=pawtab(itypat)%mesh_size
    ABI_MALLOC(ff,(mesh_size))
    do klmn=1, pawtab(itypat)%lmn2_size
      kln = pawtab(itypat)%indklmn(2,klmn) 

      ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r^2
      ff(2:mesh_size) = (pawtab(itypat)%phiphj(2:mesh_size,kln)-&
        &                pawtab(itypat)%tphitphj(2:mesh_size,kln))/&
        &               pawrad(itypat)%rad(2:mesh_size)**2
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))
 
      ilmn = pawtab(itypat)%indklmn(7,klmn)
      ilm=pawtab(itypat)%indklmn(5,klmn)
      il=pawtab(itypat)%indlmn(1,ilmn)

      jlmn = pawtab(itypat)%indklmn(8,klmn)
      jl=pawtab(itypat)%indlmn(1,jlmn)
      jm=pawtab(itypat)%indlmn(2,jlmn)

      dij_cart = czero
      do adir = 1, 3
        cme = czero
        l1b = adir_to_sij(adir)
        ilml1b = MATPACK(ilm,l1b)

        do bigl = abs(il-1),il+1
          do bigm=-bigl,bigl

            biglm = bigl*bigl + bigl + 1 + bigm
            gint = gntselect(biglm,ilml1b)
            if (gint > 0) then
              do ldir = 1, 3
                if (abs(dtset%nucdipmom(ldir,iat)) .LT. tol8) cycle
                ! compute <L_dir>
                call slxyzs(bigl,bigm,ldir,jl,jm,orbl_me)
                cme = cme + j_dpc*a2*realgnt(gint)*orbl_me*dtset%nucdipmom(ldir,iat)*intg
              end do ! end loop over ldir in L.m
            end if ! end check for nonzero Gaunt coeff
          end do ! end loop over big m
        end do ! end loop over big L
        dij_cart(adir) = cme
      end do ! end loop over adir
      dlanp(klmn,iat,1:3) = MATMUL(TRANSPOSE(gprimd),dij_cart)
    end do ! end loop over klmn
    ABI_FREE(ff)
  end do ! end loop over atoms
 
end subroutine dl_Anp


!!****f* ABINIT/d2lr_p2
!! NAME
!! d2lr_p2
!!
!! FUNCTION
!! Compute onsite r_b p^2/2 r_g - r_g p^2/2 r_b == -i <LR> 
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine d2lr_p2(dtset,gprimd,lmn2max,lr,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max
  type(dataset_type),intent(in) :: dtset

  !arrays
  complex(dpc),intent(out) :: lr(lmn2max,dtset%natom,3)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,ilmn,il,im,itypat,jlmn,jl,jm,klmn,kln,mesh_size
  real(dp) :: intg
  complex(dpc) :: cme,orbl_me

  !arrays
  complex(dpc) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

  lr=czero
  do itypat=1,dtset%ntypat
    mesh_size=pawtab(itypat)%mesh_size
    ABI_MALLOC(ff,(mesh_size))
    do klmn=1, pawtab(itypat)%lmn2_size

      ilmn = pawtab(itypat)%indklmn(7,klmn)
      il=pawtab(itypat)%indlmn(1,ilmn)
      im=pawtab(itypat)%indlmn(2,ilmn)

      jlmn = pawtab(itypat)%indklmn(8,klmn)
      jl=pawtab(itypat)%indlmn(1,jlmn)
      jm=pawtab(itypat)%indlmn(2,jlmn)

      if ( il /= jl ) cycle ! <l'm'|L|lm> = 0 if l' /= l
      if ( il == 0 ) cycle ! <00|L|00> = 0

      kln = pawtab(itypat)%indklmn(2,klmn) 
      ff(1:mesh_size) = pawtab(itypat)%phiphj(1:mesh_size,kln)-pawtab(itypat)%tphitphj(1:mesh_size,kln)
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))

      do adir = 1, 3
      ! compute <L_dir>
        call slxyzs(il,im,adir,jl,jm,orbl_me)
        dij_cart(adir) = -j_dpc*orbl_me*intg
      end do ! end loop over adir

      ! convert to crystal frame
      dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

      do iat=1,dtset%natom
          if(dtset%typat(iat) .EQ. itypat) then
            lr(klmn,iat,1:3) = dij_red(1:3)
          end if
      end do
    end do ! end loop over klmn
    ABI_FREE(ff)
  end do ! end loop over atoms
 
end subroutine d2lr_p2

!!****f* ABINIT/apply_d2lr_term_k
!! NAME
!! apply_d2lr_term_k
!!
!! FUNCTION
!! Compute contributions due to onsite r_b H r_g - r_g H r_b to orbital magnetization at given k point
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine apply_d2lr_term_k(atindx,cpre,cprj_k,dtset,lmn2max,mterm,nband_k,omm,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  complex(dpc),intent(in) :: cpre
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: omm(2,nband_k,3)
  complex(dpc),intent(in) :: mterm(lmn2max,dtset%natom,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  complex(dpc) :: cpb,cpk,cme

  !arrays

!--------------------------------------------------------------------

  omm = zero 
  do nn = 1, nband_k
    do adir = 1, 3
      do iat=1,dtset%natom
        iatom = atindx(iat)
        itypat = dtset%typat(iat)
        do klmn=1,pawtab(itypat)%lmn2_size
          ilmn = pawtab(itypat)%indklmn(7,klmn)
          jlmn = pawtab(itypat)%indklmn(8,klmn)
          cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
          cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)

          cme = conjg(cpb)*mterm(klmn,iat,adir)*cpk
          ! note the mterms are not Hermitian, they have the form j_dpc*Dij where Dij itself 
          ! is Hermitian, hence the "-" sign in the following line
          if(ilmn /= jlmn) then
            cme = cme - cpb*conjg(mterm(klmn,iat,adir))*conjg(cpk)
          end if
          omm(1,nn,adir)=omm(1,nn,adir) + real(cpre*cme)
          omm(2,nn,adir)=omm(2,nn,adir) + aimag(cpre*cme)

        end do ! end loop over klmn
      end do ! end loop over atoms
    end do ! end loop over adir
  end do ! end loop over nn
 
end subroutine apply_d2lr_term_k

!!****f* ABINIT/d2lr_Anp
!! NAME
!! d2lr_Anp
!!
!! FUNCTION
!! Compute onsite term arising from <phi|r_b A0.p r_g|phi>
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,xmpi_sum
!!
!! SOURCE

subroutine d2lr_Anp(dtset,gntselect,gprimd,lmn2_max,mpan,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2_max,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: mpan(lmn2_max,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)


  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,gint_b,gint_g,iat,il,ilm,ilmn,itypat
  integer :: jl,jlmn,jlm,klmn,kln,lb,ldir,lg,ll,llp,lm1b,lm1g,llmm,llmmp
  integer :: mesh_size,mm,mmp
  real(dp) :: a2,eabg,intg,sij
  complex(dpc) :: cme,orbl_me

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:)
  complex(dpc) :: dij_cart(3),dij_red(3)

  ! ***********************************************************************

  a2 = FineStructureConstant2
  sij = four_pi/three

  mpan= czero

  do iat=1,dtset%natom
    itypat=dtset%typat(iat)
    mesh_size=pawtab(itypat)%mesh_size
    ABI_MALLOC(ff,(mesh_size))
    do klmn=1,pawtab(itypat)%lmn2_size
      kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
      ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r
      ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln) - &
           &           pawtab(itypat)%tphitphj(2:mesh_size,kln)) / &
           &           pawrad(itypat)%rad(2:mesh_size)
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))

      ilmn = pawtab(itypat)%indklmn(7,klmn)
      il = pawtab(itypat)%indlmn(1,ilmn)
      ilm = pawtab(itypat)%indlmn(4,ilmn)
      jlmn = pawtab(itypat)%indklmn(8,klmn)
      jl = pawtab(itypat)%indlmn(1,jlmn)
      jlm = pawtab(itypat)%indlmn(4,jlmn)

      do adir = 1, 3
        cme = czero
        do bdir = 1, 3
          do gdir = 1, 3
            eabg = eijk(adir,bdir,gdir)
            if ( abs(eabg) < tol8 ) cycle

            lb = adir_to_sij(bdir)
            lm1b = MATPACK(ilm,lb)

            lg = adir_to_sij(gdir)
            lm1g = MATPACK(jlm,lg)

            do ll = abs(il-1),il+1
              do mm = -ll,ll
                llmm = ll*ll + ll + 1 + mm
                gint_b = gntselect(llmm,lm1b)
                if(gint_b == 0) cycle

                do llp = abs(jl-1),jl+1
                  do mmp = -llp,llp
                    llmmp = llp*llp + llp + 1 + mmp
                    gint_g = gntselect(llmmp,lm1g)
                    if(gint_g == 0) cycle

                    do ldir = 1, 3
                      call slxyzs(ll,mm,ldir,llp,mmp,orbl_me)
                      cme = cme + a2*sij*eabg*intg*realgnt(gint_b)*realgnt(gint_g)*&
                        & orbl_me*dtset%nucdipmom(ldir,iat)
                    end do ! end loop over ldir
                  end do
                end do ! end loop over llp
              end do 
            end do ! end loop over ll
          end do ! end  gdir
        end do ! bdir
        dij_cart(adir) = cme
      end do ! end loop over adir
      dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
      mpan(klmn,iat,1:3) = dij_red(1:3)

    end do ! end loop over klmn
    ABI_FREE(ff)
  end do ! end loop over atom

end subroutine d2lr_Anp
!!***

!!****f* ABINIT/orbmag_tt_output
!! NAME
!! orbmag_tt_output
!!
!! FUNCTION
!! This routine outputs orbmag terms tailored for the tt ddk routine
!!
!! COPYRIGHT
!! Copyright (C) 2003-2022 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine orbmag_tt_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)


 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nband_k,nterms
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: orbmag_terms(2,nband_k,3,nterms),orbmag_trace(2,3,nterms)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,ikpt,iterms
 integer,parameter :: ibcpw=1,ibcdpdp=2,ibcdpu=3,ibcdudu=4,ibcd=5
 integer,parameter :: iompw=6,iomdpdp=7,iomdpu=8,iomdudu=9,iomlr=10,iomlmb=11,iomanp=12
 integer,parameter :: iomdlanp=13,iomdlp2=14,iomdlvhnzc=15,iomdlvh=16,iomdlq=17
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(2,nband_k,3),berry_total(2,3),orbmag_bb(2,nband_k,3),orbmag_total(2,3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = ibcd+1, nterms
   orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do iband=1, nband_k
     if (iterms == iomlmb) cycle
     orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
   end do
 end do
 berry_bb=zero;berry_total=zero
 do iterms = ibcpw,ibcd
   berry_total(1:2,1:3)=berry_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do iband=1, nband_k
     berry_bb(1:2,iband,1:3) = berry_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
   end do
 end do

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a,a)')' Orbital magnetic moment computed with DFPT derivative wavefunctions ',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a)')' Orbital magnetic moment, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (orbmag_total(1,adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')
 write(message,'(a)')ch10
 call wrtout(ab_out,message,'COLL')
 write(message,'(a)')' Integral of Berry curvature, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (berry_total(1,adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .GE. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '       <du/dk|H+E|du/dk> : ',(orbmag_trace(1,adir,iompw),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     onsite <dp|D+Eq|dp> : ',(orbmag_trace(1,adir,iomdpdp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '      onsite <dp|D+Eq|u> : ',(orbmag_trace(1,adir,iomdpu),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     onsite <du|D+Eq|du> : ',(orbmag_trace(1,adir,iomdudu),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                   dl Eq : ',(orbmag_trace(1,adir,iomdlq),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                   dl vh : ',(orbmag_trace(1,adir,iomdlvh),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                dl vhnzc : ',(orbmag_trace(1,adir,iomdlvhnzc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                 dl An.p : ',(orbmag_trace(1,adir,iomdlanp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                dl p^2/2 : ',(orbmag_trace(1,adir,iomdlp2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <u|p>r_b p^2/2 r_g<p|u> : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <u|p>r_b A0.p r_g <p|u> : ',(orbmag_trace(1,adir,iomanp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '               Lamb core : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    <du/dk|du/dk> : ',(orbmag_trace(1,adir,ibcpw),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  <u|dp><dp|u>Qij : ',(orbmag_trace(1,adir,ibcdpdp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  <du|p><dp|u>Qij : ',(orbmag_trace(1,adir,ibcdpu),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  <du|p><p|du>Qij : ',(orbmag_trace(1,adir,ibcdudu),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     <du/dk|u>dij : ',(orbmag_trace(1,adir,ibcd),adir=1,3)
   call wrtout(ab_out,message,'COLL')
 end if

 if(abs(dtset%orbmag) .EQ. 3) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Term-by-term breakdowns for each band : '
   call wrtout(ab_out,message,'COLL')
   do iband = 1, nband_k
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,i2,a,i2)') ' band ',iband,' of ',nband_k
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Orbital magnetic moment : ',(orbmag_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '       <du/dk|H+E|du/dk> : ',(orbmag_terms(1,iband,adir,iompw),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    onsite <dp|D+Eq|dp>  : ',(orbmag_terms(1,iband,adir,iomdpdp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      onsite <dp|D+Eq|u> : ',(orbmag_terms(1,iband,adir,iomdpu),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '     onsite <du|D+Eq|du> : ',(orbmag_terms(1,iband,adir,iomdudu),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                   dl Eq : ',(orbmag_terms(1,iband,adir,iomdlq),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                   dl vh : ',(orbmag_terms(1,iband,adir,iomdlvh),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                dl vhnzc : ',(orbmag_terms(1,iband,adir,iomdlvhnzc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                 dl An.p : ',(orbmag_terms(1,iband,adir,iomdlanp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                dl p^2/2 : ',(orbmag_terms(1,iband,adir,iomdlp2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' <u|p>r_b p^2/2 r_g<p|u> : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' <u|p>r_b A0.p r_g <p|u> : ',(orbmag_terms(1,iband,adir,iomanp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    <du/dk|du/dk> : ',(orbmag_terms(1,iband,adir,ibcpw),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  <u|dp><dp|u>Qij : ',(orbmag_terms(1,iband,adir,ibcdpdp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  <du|p><p|du>Qij : ',(orbmag_terms(1,iband,adir,ibcdpu),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  <du|p><p|du>Qij : ',(orbmag_terms(1,iband,adir,ibcdudu),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '     <du/dk|u>dij : ',(orbmag_terms(1,iband,adir,ibcd),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_tt_output
!!***

!****f* ABINIT/local_rhoij
!!! NAME
!! local_rhoij
!!
!! FUNCTION
!! make local rhoij 
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine local_rhoij(atindx,dtset,cprj,lmn2max,mcprj,mpi_enreg,rhoij,pawtab,ucvol)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: lmn2max,mcprj
 real(dp),intent(in) :: ucvol
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg

 !arrays
 integer,intent(in) :: atindx(dtset%natom)
 real(dp),intent(out) :: rhoij(2,lmn2max,dtset%natom)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

 !Local
 !scalars
 integer :: buff_size,iat,iatom,icprj,ierr,ikpt,ilmn,itypat,jlmn,klmn
 integer :: me,nband_k,nn,nproc,spaceComm
 real(dp) :: trnrm
 complex(dpc) :: cpi,cpj,cterm

 !arrays
 real(dp),allocatable :: buffer1(:),buffer2(:)

 !----------------------------------------------

 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me = mpi_enreg%me_kpt
 nband_k = dtset%mband
 icprj = 0 
 
 rhoij = zero

 do ikpt = 1, dtset%nkpt

   ! if the current kpt is not on the current processor, cycle
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

   ! trace norm: assume occupation of two for each band and weight by kpts.
   ! division by ucvol arises from integration over BZ (so multiplication by 
   ! \Omega_{BZ} or division by \Omega
   trnrm = two*dtset%wtk(ikpt)/ucvol

   do iat=1,dtset%natom
     iatom=atindx(iat)
     itypat=dtset%typat(iat)
     do nn = 1, nband_k
       do klmn = 1, pawtab(itypat)%lmn2_size
         ilmn = pawtab(itypat)%indklmn(7,klmn)
         jlmn = pawtab(itypat)%indklmn(8,klmn)
         cpi = cmplx(cprj(iatom,icprj+nn)%cp(1,ilmn),cprj(iatom,icprj+nn)%cp(2,ilmn),KIND=dpc)
         cpj = cmplx(cprj(iatom,icprj+nn)%cp(1,jlmn),cprj(iatom,icprj+nn)%cp(2,jlmn),KIND=dpc)     
         cterm = trnrm*conjg(cpi)*cpj
         rhoij(1,klmn,iatom) = rhoij(1,klmn,iatom) + real(cterm)
         rhoij(2,klmn,iatom) = rhoij(2,klmn,iatom) + aimag(cterm)
       end do ! end loop over klmn
     end do ! end loop over nn
   end do ! end loop over atoms
   
   icprj = icprj + nband_k

 end do ! end loop over k points

 if (nproc > 1) then
   buff_size=size(rhoij)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = reshape(rhoij,(/2*lmn2max*dtset%natom/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   rhoij(1:2,1:lmn2max,1:dtset%natom)=reshape(buffer2,(/2,lmn2max,dtset%natom/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

end subroutine local_rhoij


end module m_orbmag
