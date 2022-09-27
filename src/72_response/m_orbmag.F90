!*** ABINIT/m_orbmag
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

#define MATPACK(row,col) (MAX(row,col)*(MAX(row,col)-1)/2 + MIN(row,col))
#define LMPACK(l,m) (l*l+l+1+m)

module m_orbmag

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : projbd
  use m_dfpt_nstwf,       only : gaugetransfo
  use m_fft,              only : fftpac
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type
  use m_kg,               only : getph,mkkin,mkkpg,ph1d3d
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle
  use m_nonlop,           only : nonlop
  use m_paw_an,           only : paw_an_type
  use m_pawang,           only : pawang_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get, pawcprj_put
  use m_paw_dmft,         only : paw_dmft_type
  use m_pawfgr,           only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_paw_denpot,       only : pawdensities
  use m_paw_occupancies,  only : pawmkrhoij
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen
  use m_pawrhoij,         only : pawrhoij_alloc, pawrhoij_free, pawrhoij_type
  use m_paw_sphharm,      only : setsym_ylm,slxyzs,realgaunt
  use m_pawtab,           only : pawtab_type
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
 
  ! these paramters name the various output terms                                             
  integer,parameter :: if1=1,if2=2,if3=3,if4=4
  integer,parameter :: iomlr=5,iomanp=6,iomlmb=7
  integer,parameter :: ig1=8,ig4=9,ih1=10,ih2=11,ih3=12,ih4=13
  integer,parameter :: ig2a_dij=14,ig2b_p2=15,ig2b_ap=16,ig2b_vha=17,ig2b_vhnzc=18
  integer,parameter :: ig3a_dij=19,ig3b_p2=20,ig3b_ap=21,ig3b_vha=22,ig3b_vhnzc=23,ig3b_vhnhat=24,ig2b_vhnhat=25
  integer,parameter :: nterms=25

  ! these parameters name the various d terms
  integer,parameter :: idp2=1,idpa=2,idvhnzc=3,idvha=4,idvhnhat=5,idsij=6
  integer,parameter :: ndterms=6

  ! these parameters are constants used repeatedly
  real(dp),parameter :: c2 = one/(two_pi*two_pi) ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  complex(dpc),parameter :: cbc = j_dpc/two_pi ! Berry curvature pre-factor
  complex(dpc),parameter :: com = -half*j_dpc  ! Orbital magnetism pre-factor

  ! Bound methods:

  public :: orbmag_ptpaw

  private :: make_d
  private :: make_ddir_sij
  private :: make_ddir_vhnzc
  private :: make_ddir_p2
  private :: make_ddir_ap
  private :: make_ddir_vha
  private :: make_ddir_vha2
  private :: make_ddir_vhnhat
  private :: make_ddir_vha3
  private :: make_fgh1
  private :: make_fh2
  private :: make_fh3
  private :: make_fgh4
  private :: make_g2a
  private :: make_g2b
  private :: make_g3a
  private :: make_g3b

  private :: d2lr_p2
  private :: d2lr_Anp

  private :: tdt_me
  private :: apply_d2lr_term_k
  private :: lamb_core
  private :: make_pcg1
  private :: pack_pawrhoij
  private :: orbmag_ptpaw_output
  
CONTAINS  !========================================================================================
!!***

!!****f* ABINIT/orbmag_ptpaw
!! NAME
!! orbmag_ptpaw
!!
!! FUNCTION
!! This routine computes the orbital magnetization and Berry curvature based on input 
!! wavefunctions and DDK wavefuntions. It is based on using the modern theory of
!! orbital magnetism directly.
!! It is assumed that only completely filled bands are present.
!! It uses the perturbation theory followed by PAW approach, as in Umari, Gonze, and Pasquarello,
!! PRB 69, 235102 (2004)
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
!! SOURCE

subroutine orbmag_ptpaw(cg,cg1,cprj,dtset,eigen0,gsqcut,kg,mcg,mcg1,mcprj,mpi_enreg,&
    & nfftf,ngfftf,npwarr,occ,paw_ij,paw_an,pawang,pawfgr,pawrad,pawtab,psps,rprimd,vtrial,xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcprj,mcg,mcg1,nfftf
 real(dp),intent(in) :: gsqcut
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps

 !arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),ngfftf(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,buff_size,choice,cpopt,dimffnl,exchn2n3d,iat,iatom,icg,icmplx,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat,klmn,lmn2max
 integer :: me,mcgk,my_lmax,my_nspinor,nband_k,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: eig_k(:),ff_k(:,:,:),ffnl_k(:,:,:,:),gg_k(:,:,:),hh_k(:,:,:) 
 real(dp),allocatable :: kinpw(:),kpg_k(:,:),orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 complex(dpc),allocatable :: dterms(:,:,:,:,:)
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

 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
 ABI_MALLOC(cprj_k,(dtset%natom,dtset%mband))
 call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimlmn)
 ABI_MALLOC(cprj1_k,(dtset%natom,dtset%mband,3))
 do adir = 1, 3
   call pawcprj_alloc(cprj1_k(:,:,adir),0,dimlmn)
 end do
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,0,dimlmn)

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

 ABI_MALLOC(dterms,(ndterms,psps%lmnmax,psps%lmnmax,dtset%natom,3))
 call make_d(atindx,atindx1,cprj,dimlmn,dterms,dtset,gprimd,mcprj,&
   & mpi_enreg,occ,paw_an,pawang,pawrad,pawtab,psps)

 ABI_MALLOC(mp2,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling d2lr_p2...'
 call d2lr_p2(dtset,gprimd,lmn2max,mp2,pawrad,pawtab)

 ABI_MALLOC(mpan,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling d2lr_Anp...'
 call d2lr_Anp(dtset,gntselect,gprimd,lmn2max,mpan,my_lmax,pawrad,pawtab,realgnt)

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

   ! retrieve zeroth order eigenvalues
   ABI_MALLOC(eig_k,(nband_k))
   eig_k(1:nband_k) = eigen0(nband_k*(ikpt-1)+1:nband_k*ikpt)

   ! retrieve cprj_k
   call pawcprj_get(atindx,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
     & dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

   ! compute P_c|cg1>
   ABI_MALLOC(pcg1_k,(2,mcgk,3))
   call make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,ikpt,isppol,&
     & mcgk,mpi_enreg,nband_k,npw_k,my_nspinor,pcg1_k)

   ! compute <p|Pc cg1> cprjs
   ABI_MALLOC(cwavef,(2,npw_k))
   choice = 1
   cpopt = 0
   idir = 0
   do nn = 1, nband_k
     do adir = 1, 3
       cwavef(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir)
       call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl_k,idir,psps%indlmn,istwf_k,&
         & kg_k,kpg_k,kpoint,psps%lmnmax,dtset%mgfft,mpi_enreg,dtset%natom,nattyp,dtset%ngfft,&
         & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)
       call pawcprj_put(atindx,cwaveprj,cprj1_k(:,:,adir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
         & dtset%mkmem,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)
     end do
   end do
   ABI_FREE(cwavef)

   !--------------------------------------------------------------------------------
   ! Finally ready to compute contributions to orbital magnetism and Berry curvature
   !--------------------------------------------------------------------------------
   
   ABI_MALLOC(ff_k,(2,nband_k,3))
   ABI_MALLOC(gg_k,(2,nband_k,3))
   ABI_MALLOC(hh_k,(2,nband_k,3))

   call make_fgh1(atindx,cprj1_k,dimlmn,dtset,eig_k,ff_k,gg_k,gs_hamk,hh_k,ikpt,isppol,mcgk,&
     & mpi_enreg,my_nspinor,nband_k,npw_k,pcg1_k)
   orbmag_terms(:,:,:,if1) = orbmag_terms(:,:,:,if1) + trnrm*ff_k
   orbmag_terms(:,:,:,ih1) = orbmag_terms(:,:,:,ih1) + trnrm*hh_k
   orbmag_terms(:,:,:,ig1) = orbmag_terms(:,:,:,ig1) + trnrm*gg_k
  
   call make_fh2(atindx,cprj_k,cprj1_k,dterms(idsij,:,:,:,:),dtset,eig_k,ff_k,hh_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,if2) = orbmag_terms(:,:,:,if2) + trnrm*ff_k
   orbmag_terms(:,:,:,ih2) = orbmag_terms(:,:,:,ih2) + trnrm*hh_k
 
   call make_fh3(atindx,cprj_k,dterms(idsij,:,:,:,:),dtset,eig_k,ff_k,hh_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,if3) = orbmag_terms(:,:,:,if3) + trnrm*ff_k
   orbmag_terms(:,:,:,ih3) = orbmag_terms(:,:,:,ih3) + trnrm*hh_k
   
   call make_fgh4(atindx,cprj_k,dterms(idsij,:,:,:,:),dtset,eig_k,ff_k,gg_k,hh_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,if4) = orbmag_terms(:,:,:,if4) + trnrm*ff_k
   orbmag_terms(:,:,:,ih4) = orbmag_terms(:,:,:,ih4) + trnrm*hh_k
   orbmag_terms(:,:,:,ig4) = orbmag_terms(:,:,:,ig4) + trnrm*gg_k
   
   call make_g2a(atindx,cprj_k,cprj1_k,dtset,gg_k,lmn2max,nband_k,paw_ij,pawtab)
   orbmag_terms(:,:,:,ig2a_dij) = orbmag_terms(:,:,:,ig2a_dij) + trnrm*gg_k
   
   call make_g2b(atindx,cprj_k,cprj1_k,dterms(idp2,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig2b_p2) = orbmag_terms(:,:,:,ig2b_p2) + trnrm*gg_k

   call make_g2b(atindx,cprj_k,cprj1_k,dterms(idpa,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig2b_ap) = orbmag_terms(:,:,:,ig2b_ap) + trnrm*gg_k

   call make_g2b(atindx,cprj_k,cprj1_k,dterms(idvha,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig2b_vha) = orbmag_terms(:,:,:,ig2b_vha) + trnrm*gg_k
   
   call make_g2b(atindx,cprj_k,cprj1_k,dterms(idvhnhat,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig2b_vhnhat) = orbmag_terms(:,:,:,ig2b_vhnhat) + trnrm*gg_k

   call make_g2b(atindx,cprj_k,cprj1_k,dterms(idvhnzc,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig2b_vhnzc) = orbmag_terms(:,:,:,ig2b_vhnzc) + trnrm*gg_k

   call make_g3a(atindx,cprj_k,dtset,gg_k,lmn2max,nband_k,paw_ij,pawtab)
   orbmag_terms(:,:,:,ig3a_dij) = orbmag_terms(:,:,:,ig3a_dij) + trnrm*gg_k
   
   call make_g3b(atindx,cprj_k,dterms(idp2,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig3b_p2) = orbmag_terms(:,:,:,ig3b_p2) + trnrm*gg_k
   
   call make_g3b(atindx,cprj_k,dterms(idpa,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig3b_ap) = orbmag_terms(:,:,:,ig3b_ap) + trnrm*gg_k
   
   call make_g3b(atindx,cprj_k,dterms(idvha,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig3b_vha) = orbmag_terms(:,:,:,ig3b_vha) + trnrm*gg_k
   
   call make_g3b(atindx,cprj_k,dterms(idvhnhat,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig3b_vhnhat) = orbmag_terms(:,:,:,ig3b_vhnhat) + trnrm*gg_k
   
   call make_g3b(atindx,cprj_k,dterms(idvhnzc,:,:,:,:),dtset,gg_k,psps%lmnmax,nband_k,pawtab)
   orbmag_terms(:,:,:,ig3b_vhnzc) = orbmag_terms(:,:,:,ig3b_vhnzc) + trnrm*gg_k
   
   !--------------------------------------------------------------------------------
   ! onsite <phi|r_b p^2/2 r_g>
   !--------------------------------------------------------------------------------
   call apply_d2lr_term_k(atindx,cprj_k,dtset,iomlr,lmn2max,mp2,nband_k,gg_k,pawtab)
   orbmag_terms(:,:,:,iomlr) = orbmag_terms(:,:,:,iomlr) + trnrm*gg_k
   
   !--------------------------------------------------------------------------------
   ! onsite <phi|r_b p.A0 r_g> 
   !--------------------------------------------------------------------------------
   call apply_d2lr_term_k(atindx,cprj_k,dtset,iomanp,lmn2max,mpan,nband_k,gg_k,pawtab)
   orbmag_terms(:,:,:,iomanp) = orbmag_terms(:,:,:,iomanp) + trnrm*gg_k

   ABI_FREE(ff_k)
   ABI_FREE(gg_k)
   ABI_FREE(hh_k)

   icg = icg + npw_k*nband_k
   ikg = ikg + npw_k
   icprj = icprj + nband_k

   ABI_FREE(cg_k)
   ABI_FREE(cg1_k)
   ABI_FREE(pcg1_k)
   ABI_FREE(eig_k)
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
   if (iterm .EQ. iomlmb) cycle
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
 do iterm = iomlr, nterms
   if (iterm .EQ. iomlmb) cycle
   orbmag_terms(:,:,:,iterm) = ucvol*orbmag_terms(:,:,:,iterm)
 end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace = orbmag_trace + orbmag_terms(:,nn,:,:)
 end do

! get the Lamb term
 call lamb_core(atindx,dtset,omlamb)
 orbmag_trace(:,:,iomlmb) = omlamb
!
 call orbmag_ptpaw_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)

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

 ABI_FREE(dterms)
 ABI_FREE(mp2)
 ABI_FREE(mpan)
 
 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

end subroutine orbmag_ptpaw
!!***

!!****f* ABINIT/make_pcg1
!! NAME
!! make_pcg1
!!
!! FUNCTION
!! compute Pc|cg1> from |cg1> and |cg>
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
!! see Audouze et al PRB 78, 035105 (2008) Eq. 40
!!
!! SOURCE

subroutine make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,&
    & ikpt,isppol,mcgk,mpi_enreg,nband_k,npw_k,my_nspinor,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),cg1_k(2,mcgk,3)
  real(dp),intent(out) :: pcg1_k(2,mcgk,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)

  !Local variables -------------------------
  !scalars
  integer :: adir,choice,cpopt,iband,jband
  integer :: ndat,nnlout,paw_opt,signs,tim_nonlop
  real(dp) :: doti,dotr,s0,s1
  real(dp) :: cg1,pcg1,vcg
  !arrays
  real(dp) :: lambda(1)
  real(dp),allocatable :: cwavef(:,:),enlout(:),svectout(:,:)
  real(dp),allocatable :: vcg1(:,:),vectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

  choice = 5 ! dS/dk in nonlop
  cpopt = 4 ! cprj and derivatives already in memory
  paw_opt = 3
  signs = 2
  tim_nonlop = 0
  lambda = zero
  nnlout = 0
  ndat = 1

  ABI_MALLOC(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,3,dimlmn)
  ABI_MALLOC(cwavef,(2,npw_k))
  ABI_MALLOC(vectout,(2,npw_k))
  ABI_MALLOC(svectout,(2,npw_k))
  ABI_MALLOC(vcg1,(2,npw_k))

  pcg1_k = zero

  do adir = 1, 3

    do iband = 1, nband_k

      cwavef(1:2,1:npw_k)=cg_k(1:2,(iband-1)*npw_k+1:iband*npw_k)
      call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,iband,0,ikpt,0,isppol,dtset%mband,&
        & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
      call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,adir,lambda,mpi_enreg,ndat,&
        & nnlout,paw_opt,signs,svectout,tim_nonlop,cwavef,vectout)

      ! form vcg1 = -1/2 \sum |u_j^0><u_j^0|S^1|u_i^0>, the valence band part of cg1
      vcg1 = zero
      do jband = 1, nband_k
        cwavef(1:2,1:npw_k)=cg_k(1:2,(jband-1)*npw_k+1:jband*npw_k)
        dotr = DOT_PRODUCT(cwavef(1,:),svectout(1,:))+DOT_PRODUCT(cwavef(2,:),svectout(2,:))
        doti = DOT_PRODUCT(cwavef(1,:),svectout(2,:))-DOT_PRODUCT(cwavef(2,:),svectout(1,:))
        vcg1(1,:) = vcg1(1,:) - half*( dotr*cwavef(1,:) - doti*cwavef(2,:))
        vcg1(2,:) = vcg1(2,:) - half*( dotr*cwavef(2,:) + doti*cwavef(1,:))
      end do
    
      ! subtract vcg1 from cg1_k to obtain pcg1, the conduction band part of cg1 
      pcg1_k(1:2,(iband-1)*npw_k+1:iband*npw_k,adir) =cg1_k(1:2,(iband-1)*npw_k+1:iband*npw_k,adir)-&
       &  vcg1(1:2,1:npw_k)
      ! here's the unprojected version, for testing purposes
      ! pcg1_k(1:2,(iband-1)*npw_k+1:iband*npw_k,adir) =cg1_k(1:2,(iband-1)*npw_k+1:iband*npw_k,adir)

    end do
  end do

  ABI_FREE(cwavef)
  ABI_FREE(vectout)
  ABI_FREE(vcg1)
  ABI_FREE(svectout)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)
 
end subroutine make_pcg1
!!***

!!****f* ABINIT/tdt_me
!! NAME
!! tdt_me
!!
!! FUNCTION
!! compute the matrix element <bra|T^\dag (\partial T)|ket>
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx
!! bra_cprj_k : cprj for bra wavefunctions
!! bra_n : index of band to use on bra side
!! ddir : set of <phi|r|phi> - <t phi|r|t phi>
!! dtset
!! idir : direction (1..3) for derivatives
!! ket_cprj_k : cprj for ket wavefunctions
!! ket_n index of band to use on ket side
!! lmn2max
!! nband_k
!! pawtab
!! tdtme : output matrix element
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! SOURCE

subroutine tdt_me(atindx,bra_cprj_k,bra_n,ddir,dtset,idir,ket_cprj_k,&
    & ket_n,lmnmax,nband_k,pawtab,tdtme)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bra_n,idir,ket_n,lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  complex(dpc),intent(out) :: tdtme
  type(pawcprj_type),intent(in) :: bra_cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: ket_cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,ilmn,itypat,jlmn,klmn
  real(dp) :: qij
  complex(dpc) :: cpdj,cpi,cpj
  
  !arrays

!--------------------------------------------------------------------

  tdtme = czero
  do iat=1,dtset%natom
    iatom=atindx(iat)
    itypat=dtset%typat(iat)
    do ilmn = 1, pawtab(itypat)%lmn_size
      cpi = CMPLX(bra_cprj_k(iatom,bra_n)%cp(1,ilmn),bra_cprj_k(iatom,bra_n)%cp(2,ilmn),KIND=dpc)
      do jlmn = 1, pawtab(itypat)%lmn_size
        cpj = CMPLX(ket_cprj_k(iatom,ket_n)%cp(1,jlmn),ket_cprj_k(iatom,ket_n)%cp(2,jlmn),KIND=dpc)
        klmn = MATPACK(ilmn,jlmn)
        qij = pawtab(itypat)%sij(klmn)
        cpdj = CMPLX(ket_cprj_k(iatom,ket_n)%dcp(1,idir,jlmn),ket_cprj_k(iatom,ket_n)%dcp(2,idir,jlmn),KIND=dpc)

        tdtme = tdtme + CONJG(cpi)*qij*cpdj + &
          & CONJG(cpi)*(-j_dpc)*ddir(ilmn,jlmn,iat,idir)*cpj

      end do ! jlmn
    end do ! ilmn
  end do ! iat

end subroutine tdt_me
!!***


!!****f* ABINIT/make_fh2
!! NAME
!! make_fh2
!!
!! FUNCTION
!! compute term 2 contribution to f and h functions
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
!! SOURCE

subroutine make_fh2(atindx,cprj_k,cprj1_k,ddir,dtset,eig_k,ff_k,hh_k,lmnmax,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  real(dp),intent(out) :: ff_k(2,nband_k,3),hh_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: eabg,qij
  complex(dpc) :: cpdb,cpdg,cpi,cpj,ct,ct2,ct3,pcbi,pcgj
  
  !arrays

!--------------------------------------------------------------------

  ff_k = zero; hh_k = zero
  do nn = 1, nband_k
    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          call tdt_me(atindx,cprj1_k(:,:,bdir),nn,ddir,dtset,gdir,cprj_k,&
              & nn,lmnmax,nband_k,pawtab,ct2)
          call tdt_me(atindx,cprj1_k(:,:,gdir),nn,ddir,dtset,bdir,cprj_k,&
              & nn,lmnmax,nband_k,pawtab,ct3)

          ! Chern number (f function) weighted by 1
          ct = c2*cbc*eabg*(ct2+CONJG(ct3))
          ff_k(1,nn,adir) = ff_k(1,nn,adir) + REAL(ct)
          ff_k(2,nn,adir) = ff_k(2,nn,adir) + AIMAG(ct)
          
          ! Mag E (h function) weighted by E_nk
          ct = c2*com*eabg*eig_k(nn)*(ct2+CONJG(ct3))
          hh_k(1,nn,adir) = hh_k(1,nn,adir) + REAL(ct)
          hh_k(2,nn,adir) = hh_k(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_fh2
!!***

!!****f* ABINIT/make_g3b
!! NAME
!! make_g3b
!!
!! FUNCTION
!! compute contribution to g3 due to r*H term
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
!! SOURCE

subroutine make_g3b(atindx,cprj_k,ddir,dtset,gg,lmnmax,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: gg(2,nband_k,3)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: eabg
  complex(dpc) :: cpi,cpj,ct,ctermdij,dbi,dgj
  
  !arrays

!--------------------------------------------------------------------

  gg = zero
  do nn = 1, nband_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          ctermdij = czero
          do iat=1,dtset%natom
            iatom=atindx(iat)
            itypat=dtset%typat(iat)
            do ilmn = 1, pawtab(itypat)%lmn_size
              do jlmn = 1, pawtab(itypat)%lmn_size
                klmn = MATPACK(ilmn,jlmn)

                cpi = CMPLX(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
                cpj = CMPLX(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)

                dbi = CMPLX(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
                dgj = CMPLX(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)
             
                ctermdij = ctermdij + (-j_dpc)*ddir(ilmn,jlmn,iat,gdir)*CONJG(dbi)*cpj + &
                  & (+j_dpc)*ddir(ilmn,jlmn,iat,bdir)*CONJG(cpi)*dgj
            
              end do ! jlmn
            end do ! ilmn
          end do ! iat
        
          ct = c2*com*eabg*ctermdij 
          gg(1,nn,adir) = gg(1,nn,adir) + REAL(ct)
          gg(2,nn,adir) = gg(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_g3b
!!***


!!****f* ABINIT/make_g2b
!! NAME
!! make_g2b
!!
!! FUNCTION
!! compute contribution to g2 due to r*H term
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
!! SOURCE

subroutine make_g2b(atindx,cprj_k,cprj1_k,ddir,dtset,gg,lmnmax,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: gg(2,nband_k,3)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: eabg
  complex(dpc) :: cpi,cpj,ct,ctermdij,dbi,dgj
  
  !arrays

!--------------------------------------------------------------------

  gg = zero
  do nn = 1, nband_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          ctermdij = czero
          do iat=1,dtset%natom
            iatom=atindx(iat)
            itypat=dtset%typat(iat)
            do ilmn = 1, pawtab(itypat)%lmn_size
              do jlmn = 1, pawtab(itypat)%lmn_size
                klmn = MATPACK(ilmn,jlmn)

                cpi = CMPLX(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
                cpj = CMPLX(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)

                dbi = CMPLX(cprj1_k(iatom,nn,bdir)%cp(1,ilmn),cprj1_k(iatom,nn,bdir)%cp(2,ilmn),KIND=dpc)
                dgj = CMPLX(cprj1_k(iatom,nn,gdir)%cp(1,jlmn),cprj1_k(iatom,nn,gdir)%cp(2,jlmn),KIND=dpc)
             
                ctermdij = ctermdij + (-j_dpc)*ddir(ilmn,jlmn,iat,gdir)*CONJG(dbi)*cpj + &
                  & (+j_dpc)*ddir(ilmn,jlmn,iat,bdir)*CONJG(cpi)*dgj
            
              end do ! jlmn
            end do ! ilmn
          end do ! iat
        
          ct = c2*com*eabg*ctermdij 
          gg(1,nn,adir) = gg(1,nn,adir) + REAL(ct)
          gg(2,nn,adir) = gg(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_g2b
!!***

!!****f* ABINIT/make_g2a
!! NAME
!! make_g2a
!!
!! FUNCTION
!! compute contribution to g2 due to simple Dij terms only
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
!! SOURCE

subroutine make_g2a(atindx,cprj_k,cprj1_k,dtset,gg,lmn2max,nband_k,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: gg(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: eabg
  logical :: cplex_dij
  complex(dpc) :: ct,ctermdij,dij
  
  !arrays
  complex(dpc) :: dup(2,2),udp(2,2)

!--------------------------------------------------------------------

  cplex_dij = (paw_ij(1)%cplex_dij .EQ. 2)

  gg = zero
  do nn = 1, nband_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          ctermdij = czero
          do iat=1,dtset%natom
            iatom=atindx(iat)
            itypat=dtset%typat(iat)
            do klmn = 1, pawtab(itypat)%lmn2_size
              ilmn = pawtab(itypat)%indklmn(7,klmn)
              jlmn = pawtab(itypat)%indklmn(8,klmn)
             
              if (cplex_dij) then
                dij = cmplx(paw_ij(iat)%dij(2*klmn-1,1),paw_ij(iat)%dij(2*klmn,1),KIND=dpc)
              else
                dij = cmplx(paw_ij(iat)%dij(klmn,1),zero,KIND=dpc)
              end if

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

              ctermdij = ctermdij + dij*(CONJG(udp(1,1))*dup(2,2)+CONJG(dup(1,1))*udp(2,2))
              if (ilmn /= jlmn) then
                ctermdij = ctermdij + CONJG(dij)*(CONJG(udp(1,2))*dup(2,1)+CONJG(dup(1,2))*udp(2,1))
              end if
          
            end do ! klmn
          end do ! iat
        
          ct = c2*com*eabg*ctermdij 
          gg(1,nn,adir) = gg(1,nn,adir) + REAL(ct)
          gg(2,nn,adir) = gg(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_g2a
!!***


!!****f* ABINIT/make_g3a
!! NAME
!! make_g3a
!!
!! FUNCTION
!! compute term 3 contribution to g due to simple Dij terms only
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
!! SOURCE

subroutine make_g3a(atindx,cprj_k,dtset,gg,lmn2max,nband_k,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: gg(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: eabg
  logical :: cplex_dij
  complex(dpc) :: cpdb,cpdg,ct,ctermdij,dij
  
  !arrays

!--------------------------------------------------------------------

  cplex_dij = (paw_ij(1)%cplex_dij .EQ. 2)

  gg = zero
  do nn = 1, nband_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          ctermdij = czero
          do iat=1,dtset%natom
            iatom=atindx(iat)
            itypat=dtset%typat(iat)
            do klmn = 1, pawtab(itypat)%lmn2_size
              ilmn = pawtab(itypat)%indklmn(7,klmn)
              jlmn = pawtab(itypat)%indklmn(8,klmn)
            
              if (cplex_dij) then
                dij = cmplx(paw_ij(iat)%dij(2*klmn-1,1),paw_ij(iat)%dij(2*klmn,1),KIND=dpc)
              else
                dij = cmplx(paw_ij(iat)%dij(klmn,1),zero,KIND=dpc)
              end if
          
              cpdb = CMPLX(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
              cpdg = CMPLX(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)
              ctermdij = ctermdij + CONJG(cpdb)*dij*cpdg
              
              if (ilmn /= jlmn) then
                cpdb = CMPLX(cprj_k(iatom,nn)%dcp(1,bdir,jlmn),cprj_k(iatom,nn)%dcp(2,bdir,jlmn),KIND=dpc)
                cpdg = CMPLX(cprj_k(iatom,nn)%dcp(1,gdir,ilmn),cprj_k(iatom,nn)%dcp(2,gdir,ilmn),KIND=dpc)
                ctermdij = ctermdij + CONJG(cpdb)*CONJG(dij)*cpdg
              end if
          
            end do ! klmn
          end do ! iat
        
          ct = c2*com*eabg*ctermdij 
          gg(1,nn,adir) = gg(1,nn,adir) + REAL(ct)
          gg(2,nn,adir) = gg(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_g3a
!!***

!!****f* ABINIT/make_fh3
!! NAME
!! make_fh3
!!
!! FUNCTION
!! compute term 3 contribution to f and h functions
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
!! SOURCE

subroutine make_fh3(atindx,cprj_k,ddir,dtset,eig_k,ff,hh,lmnmax,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  real(dp),intent(out) :: ff(2,nband_k,3),hh(2,nband_k,3)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn,np
  real(dp) :: eabg,qij
  complex(dpc) :: cpdb,cpdg,cpi,cpj,ctermdb,ctermdg,ctermq
  
  !arrays

!--------------------------------------------------------------------

  ff = zero; hh = zero
  do nn = 1, nband_k

    do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle

          ctermq = czero; ctermdb = czero; ctermdg = czero
          do iat=1,dtset%natom
            iatom=atindx(iat)
            itypat=dtset%typat(iat)
            do ilmn = 1, pawtab(itypat)%lmn_size
              cpi = CMPLX(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
              do jlmn = 1, pawtab(itypat)%lmn_size
                cpj = CMPLX(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)
                klmn = MATPACK(ilmn,jlmn)
                qij = pawtab(itypat)%sij(klmn)
                cpdb = CMPLX(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
                cpdg = CMPLX(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)
                ctermq = ctermq + CONJG(cpdb)*qij*cpdg
                ctermdb = ctermdb + CONJG(cpi)*(j_dpc)*ddir(ilmn,jlmn,iat,bdir)*cpdg
                ctermdg = ctermdg + CONJG(cpdb)*(-j_dpc)*ddir(ilmn,jlmn,iat,gdir)*cpj
              end do ! jlmn
            end do ! ilmn
          end do ! iat
         
          ff(1,nn,adir) = ff(1,nn,adir) + REAL(c2*cbc*eabg*(ctermq+ctermdb+ctermdg))
          ff(2,nn,adir) = ff(2,nn,adir) + AIMAG(c2*cbc*eabg*(ctermq+ctermdb+ctermdg))
          hh(1,nn,adir) = hh(1,nn,adir) + REAL(c2*com*eabg*eig_k(nn)*(ctermq+ctermdb+ctermdg))
          hh(2,nn,adir) = hh(2,nn,adir) + AIMAG(c2*com*eabg*eig_k(nn)*(ctermq+ctermdb+ctermdg))
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine make_fh3
!!***

!!****f* ABINIT/make_fgh4
!! NAME
!! make_fgh4
!!
!! FUNCTION
!! compute term 4 contribution to f, g, and h functions
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
!! SOURCE

subroutine make_fgh4(atindx,cprj_k,ddir,dtset,eig_k,ff_k,gg_k,hh_k,lmnmax,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  complex(dpc),intent(in) :: ddir(lmnmax,lmnmax,dtset%natom,3)
  real(dp),intent(out) :: ff_k(2,nband_k,3),gg_k(2,nband_k,3),hh_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,nn,np
  real(dp) :: eabg
  complex(dpc) :: cnp,cnpb,cnpg
  
  !arrays

!--------------------------------------------------------------------
  
  ff_k = zero; gg_k = zero; hh_k = zero
  do nn = 1, nband_k
   do adir = 1, 3
      do bdir = 1, 3
        do gdir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle
  
          do np = 1, nband_k
    

            call tdt_me(atindx,cprj_k,np,ddir,dtset,bdir,cprj_k,&
              & nn,lmnmax,nband_k,pawtab,cnpb)
            call tdt_me(atindx,cprj_k,np,ddir,dtset,gdir,cprj_k,&
              & nn,lmnmax,nband_k,pawtab,cnpg)

            ! Chern (term ff) weighted by 1
            cnp = c2*cbc*eabg*CONJG(cnpb)*cnpg
            ff_k(1,nn,adir) = ff_k(1,nn,adir) - REAL(cnp)
            ff_k(2,nn,adir) = ff_k(2,nn,adir) - AIMAG(cnp)
            
            ! Mag E (term hh) weighted by E_nk
            cnp = c2*com*eabg*eig_k(nn)*CONJG(cnpb)*cnpg
            hh_k(1,nn,adir) = hh_k(1,nn,adir) - REAL(cnp)
            hh_k(2,nn,adir) = hh_k(2,nn,adir) - AIMAG(cnp)
            
            ! Mag H (term gg) weighted by E_n'k
            cnp = c2*com*eabg*eig_k(np)*CONJG(cnpb)*cnpg
            gg_k(1,nn,adir) = gg_k(1,nn,adir) - REAL(cnp)
            gg_k(2,nn,adir) = gg_k(2,nn,adir) - AIMAG(cnp)

          end do ! loop over np

        end do ! gdir
      end do ! bdir
    end do ! adir
  end do !loop over bands

end subroutine make_fgh4
!!***

!!****f* ABINIT/make_fgh1
!! NAME
!! make_fgh1
!!
!! FUNCTION
!! compute term 1 contribution to f, g, and h functions
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
!! SOURCE

subroutine make_fgh1(atindx,cprj1_k,dimlmn,dtset,eig_k,ff_k,gg_k,gs_hamk,&
    & hh_k,ikpt,isppol,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: ff_k(2,nband_k,3),gg_k(2,nband_k,3),hh_k(2,nband_k,3)
  type(pawcprj_type),intent(in) ::  cprj1_k(dtset%natom,nband_k,3)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,gdir,ndat,nn,sij_opt,tim_getghc,type_calc
  real(dp) :: doti,dotr,eabg,lambda
  complex(dpc) :: cme,ct
  
  !arrays
  real(dp),allocatable :: bwavef(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),gwavef(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

  ABI_MALLOC(bwavef,(2,npw_k))
  ABI_MALLOC(gwavef,(2,npw_k))
  ABI_MALLOC(ghc,(2,npw_k))
  ABI_MALLOC(gsc,(2,npw_k))
  ABI_MALLOC(gvnlxc,(2,npw_k))
  ABI_MALLOC(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,0,dimlmn)

  cpopt = 2 ! cprj in memory
  lambda = zero
  ndat = 1
  sij_opt = 1 ! save S overlap in gsc
  tim_getghc = 0
  type_calc = 0 ! use full PAW Hamiltonian
  
  ff_k = zero; gg_k = zero; hh_k = zero
  do nn = 1, nband_k
    do gdir = 1, 3
      gwavef(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
      call pawcprj_get(atindx,cwaveprj,cprj1_k(:,:,gdir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
        & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)

      call getghc(cpopt,gwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,lambda,mpi_enreg,&
        & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)

      do bdir = 1, 3
        bwavef(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
        do adir = 1, 3
          eabg = eijk(adir,bdir,gdir)
          if (abs(eabg) < half) cycle
  
          dotr = DOT_PRODUCT(bwavef(1,:),gsc(1,:))+DOT_PRODUCT(bwavef(2,:),gsc(2,:))
          doti = DOT_PRODUCT(bwavef(1,:),gsc(2,:))-DOT_PRODUCT(bwavef(2,:),gsc(1,:))
          cme = CMPLX(dotr,doti,KIND=dpc)

          ! Chern number (function ff) weighted by 1
          ct = cbc*c2*eabg*cme
          ff_k(1,nn,adir) = ff_k(1,nn,adir) + real(ct)
          ff_k(2,nn,adir) = ff_k(2,nn,adir) + aimag(ct)
          
          ! Mag E  (function hh) weighted by E_nk
          ct = com*c2*eabg*cme*eig_k(nn)
          hh_k(1,nn,adir) = hh_k(1,nn,adir) + real(ct)
          hh_k(2,nn,adir) = hh_k(2,nn,adir) + aimag(ct)

          ! Mag H  (function gg) weighted by H
          dotr = DOT_PRODUCT(bwavef(1,:),ghc(1,:))+DOT_PRODUCT(bwavef(2,:),ghc(2,:))
          doti = DOT_PRODUCT(bwavef(1,:),ghc(2,:))-DOT_PRODUCT(bwavef(2,:),ghc(1,:))
          cme = CMPLX(dotr,doti,KIND=dpc)

          ct = com*c2*eabg*cme
          gg_k(1,nn,adir) = gg_k(1,nn,adir) + real(ct)
          gg_k(2,nn,adir) = gg_k(2,nn,adir) + aimag(ct)
 
        end do ! adir
      end do ! bdir
    end do ! gdir
  end do !loop over bands

  ABI_FREE(bwavef)
  ABI_FREE(gwavef)
  ABI_FREE(ghc)
  ABI_FREE(gsc)
  ABI_FREE(gvnlxc)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)

end subroutine make_fgh1
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
!! SOURCE

subroutine lamb_core(atindx,dtset,omlamb)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: omlamb(2,3)

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

!!****f* ABINIT/make_ddir_sij
!! NAME
!! make_ddir_sij
!!
!! FUNCTION
!! Compute onsite <phi_i|r_dir|phi_j> - <\tilde{phi}_i|r_dir|\tilde{phi}_j>
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
!! SOURCE

subroutine make_ddir_sij(ddir_sij,dtset,gprimd,lmnmax,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  real(dp),intent(in) :: gprimd(3,3)
  complex(dpc),intent(out) :: ddir_sij(lmnmax,lmnmax,dtset%natom,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,itypat,ilmn,jlmn,klmn
  real(dp) :: cdij

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: ddir_cart(3),ddir_red(3)

!--------------------------------------------------------------------

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)
 ddir_sij=czero
 do iat = 1, dtset%natom
   itypat = dtset%typat(iat)
   do ilmn=1, pawtab(itypat)%lmn_size
     do jlmn=1, pawtab(itypat)%lmn_size
       klmn=MATPACK(ilmn,jlmn)
       do adir = 1, 3
         ! note adir 1,2,3 corresponds to cartesian x,y,z. In the abinit real spherical 
         ! harmonics, x~m=1 so LM=4; y~m=-1 so LM=2; z~m=0 so LM=3
         ! adir_to_sij encodes this map
         ddir_cart(adir) = cdij*pawtab(itypat)%qijl(adir_to_sij(adir),klmn)
       end do
       ! now have ddir in cart coords, convert to crystal coords where ddk wavefunctions are
       ddir_red = MATMUL(TRANSPOSE(gprimd),ddir_cart)
       do adir = 1, 3
         ddir_sij(ilmn,jlmn,iat,adir) = CMPLX(ddir_red(adir),zero,KIND=dpc)
       end do
     end do !jlmn
   end do ! ilmn
 end do ! end loop over types
 
end subroutine make_ddir_sij
!!***

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
  integer :: adir,iat,ilmn,il,im,itypat,jlmn,jl,jm,klmn,kln,mesh_size,pwave_size
  real(dp) :: intg
  complex(dpc) :: orbl_me

  !arrays
  complex(dpc) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

  lr=czero
  do itypat=1,dtset%ntypat
    mesh_size=pawtab(itypat)%mesh_size
    pwave_size=size(pawtab(itypat)%phiphj(:,1))
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
      ff=0
      ff(2:pwave_size) = pawtab(itypat)%phiphj(2:pwave_size,kln)-pawtab(itypat)%tphitphj(2:pwave_size,kln)
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
!!***

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
!! SOURCE

subroutine apply_d2lr_term_k(atindx,cprj_k,dtset,iterm,lmn2max,mterm,nband_k,omm,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iterm,lmn2max,nband_k
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
  complex(dpc) :: cpb,cpk,cme,cpre

  !arrays

!--------------------------------------------------------------------
! iterm <= if4 is berry type, others are orb mag type
  if ( iterm <= if4 ) then
    cpre = cbc
  else
    cpre = com
  end if

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
!!***

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
  integer :: mesh_size,mm,mmp,pwave_size
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
    pwave_size=size(pawtab(itypat)%phiphj(:,1))
    ABI_MALLOC(ff,(mesh_size))
    do klmn=1,pawtab(itypat)%lmn2_size
      kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
      ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r
      ff = zero
      ff(2:pwave_size)=(pawtab(itypat)%phiphj(2:pwave_size,kln) - &
           &           pawtab(itypat)%tphitphj(2:pwave_size,kln)) / &
           &           pawrad(itypat)%rad(2:pwave_size)
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
                llmm = LMPACK(ll,mm)
                gint_b = gntselect(llmm,lm1b)
                if(gint_b == 0) cycle

                do llp = abs(jl-1),jl+1
                  do mmp = -llp,llp
                    llmmp = LMPACK(llp,mmp)
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

!!****f* ABINIT/orbmag_ptpaw_output
!! NAME
!! orbmag_ptpaw_output
!!
!! FUNCTION
!! This routine outputs orbmag terms tailored for the ptpaw ddk routine
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
!! SOURCE

subroutine orbmag_ptpaw_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)


 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nband_k,nterms
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: orbmag_terms(2,nband_k,3,nterms),orbmag_trace(2,3,nterms)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,iterms
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(2,nband_k,3),berry_total(2,3),orbmag_bb(2,nband_k,3),orbmag_total(2,3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = iomlr, nterms
   orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do iband=1, nband_k
     orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
   end do
 end do
 berry_bb=zero;berry_total=zero
 do iterms = if1,if4
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
   write(message,'(a,3es16.8)')'        <Pc du|THT|Pc_du> : ',(orbmag_trace(1,adir,ig1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'      <Pc du|T(dT)|u>*Dij : ',(orbmag_trace(1,adir,ig2a_dij),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     <Pc du|T(dT)|u>d(p2) : ',(orbmag_trace(1,adir,ig2b_p2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'    <Pc du|T(dT)|u>d(A.p) : ',(orbmag_trace(1,adir,ig2b_ap),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     <Pc du|T(dT)|u>d(vH) : ',(orbmag_trace(1,adir,ig2b_vha),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'  <Pc du|T(dT)|u>d(vHnzc) : ',(orbmag_trace(1,adir,ig2b_vhnzc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')' <Pc du|T(dT)|u>d(vHnhat) : ',(orbmag_trace(1,adir,ig2b_vhnhat),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'       <u|(dT)(dT)|u>*Dij : ',(orbmag_trace(1,adir,ig3a_dij),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'      <u|(dT)(dT)|u>d(p2) : ',(orbmag_trace(1,adir,ig3b_p2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     <u|(dT)(dT)|u>d(A.p) : ',(orbmag_trace(1,adir,ig3b_ap),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'      <u|(dT)(dT)|u>d(vH) : ',(orbmag_trace(1,adir,ig3b_vha),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   <u|(dT)(dT)|u>d(vHnzc) : ',(orbmag_trace(1,adir,ig3b_vhnzc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'  <u|(dT)(dT)|u>d(vHnhat) : ',(orbmag_trace(1,adir,ig3b_vhnhat),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'  <u|(dT)(dT)|u>dd(p2)(L) : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'    <u|(dT)(dT)|u>dd(A.p) : ',(orbmag_trace(1,adir,iomanp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   -<u|(dT)Pv H Pv(dT)|u> : ',(orbmag_trace(1,adir,ig4),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     Enk*<Pc du|TT|Pc_du> : ',(orbmag_trace(1,adir,ih1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')' Enk*<Pc du|T(dT)|u>+cnjg : ',(orbmag_trace(1,adir,ih2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'       Enk*<u|(dT)(dT)|u> : ',(orbmag_trace(1,adir,ih3),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'    -Enk*<u|(dT)Pv(dT)|u> : ',(orbmag_trace(1,adir,ih4),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'                     Lamb : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '      <Pc du|TT|Pc_du> : ',(orbmag_trace(1,adir,if1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <Pc du|T(dT)|u>+CONJG : ',(orbmag_trace(1,adir,if2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '        <u|(dT)(dT)|u> : ',(orbmag_trace(1,adir,if3),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     -<u|(dT)Pv(dT)|u> : ',(orbmag_trace(1,adir,if4),adir=1,3)
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
     write(message,'(a,3es16.8)') '        g1 : ',(orbmag_terms(1,iband,adir,ig1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   g2a_dij : ',(orbmag_terms(1,iband,adir,ig2a_dij),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    g2b_p2 : ',(orbmag_terms(1,iband,adir,ig2b_p2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    g2b_ap : ',(orbmag_terms(1,iband,adir,ig2b_ap),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   g2b_vha : ',(orbmag_terms(1,iband,adir,ig2b_vha),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' g2b_vhnzc : ',(orbmag_terms(1,iband,adir,ig2b_vhnzc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   g3a_dij : ',(orbmag_terms(1,iband,adir,ig3a_dij),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    g3b_p2 : ',(orbmag_terms(1,iband,adir,ig3b_p2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    g3b_ap : ',(orbmag_terms(1,iband,adir,ig3b_ap),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   g3b_vha : ',(orbmag_terms(1,iband,adir,ig3b_vha),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' g3b_vhnzc : ',(orbmag_terms(1,iband,adir,ig3b_vhnzc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   g3d L : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' g3d A.p : ',(orbmag_terms(1,iband,adir,iomanp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      g4 : ',(orbmag_terms(1,iband,adir,ig4),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      h1 : ',(orbmag_terms(1,iband,adir,ih1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      h2 : ',(orbmag_terms(1,iband,adir,ih2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      h3 : ',(orbmag_terms(1,iband,adir,ih3),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '      h4 : ',(orbmag_terms(1,iband,adir,ih4),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  f1 : ',(orbmag_terms(1,iband,adir,if1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  f2 : ',(orbmag_terms(1,iband,adir,if2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  f3 : ',(orbmag_terms(1,iband,adir,if3),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  f4 : ',(orbmag_terms(1,iband,adir,if4),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_ptpaw_output
!!***

!!****f* ABINIT/make_ddir_vhnzc
!! NAME
!! make_ddir_vhnzc
!!
!! FUNCTION
!! Compute onsite r*vhnzc
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

subroutine make_ddir_vhnzc(ddir_vhnzc,dtset,gntselect,gprimd,lmnmax,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vhnzc(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,iat,ilmn,itypat,jlmn,klm,klmn,kln,lm1b,mesh_size,pwave_size
  real(dp) :: cdij,intg

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

 ! dij prefactor
 cdij = sqrt(four_pi/three)
 ddir_vhnzc=czero
 do itypat = 1, dtset%ntypat
   mesh_size=pawtab(itypat)%mesh_size
   pwave_size=size(pawtab(itypat)%phiphj(:,1))
   ABI_MALLOC(ff,(mesh_size))

   do ilmn = 1, pawtab(itypat)%lmn_size
     do jlmn = 1, pawtab(itypat)%lmn_size
       klmn=MATPACK(ilmn,jlmn)

       klm = pawtab(itypat)%indklmn(1,klmn)
       kln = pawtab(itypat)%indklmn(2,klmn)

       ff = zero
       ff(2:pwave_size) = pawtab(itypat)%phiphj(2:pwave_size,kln)*pawtab(itypat)%vhnzc(2:pwave_size)
       ff(2:pwave_size) = ff(2:pwave_size) - pawtab(itypat)%tphitphj(2:pwave_size,kln)*pawtab(itypat)%vhtnzc(2:pwave_size)
       ff(2:pwave_size) = ff(2:pwave_size)*pawrad(itypat)%rad(2:pwave_size)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(intg,ff,pawrad(itypat))

       dij_cart = zero
       do adir = 1, 3
         lm1b = adir_to_sij(adir)
         gint = gntselect(lm1b,klm)
         if (gint == 0) cycle
         dij_cart(adir) = cdij*realgnt(gint)*intg
       end do
       dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
       do adir = 1, 3
         do iat = 1, dtset%natom
           if(dtset%typat(iat) == itypat) then
             ddir_vhnzc(ilmn,jlmn,iat,adir) = CMPLX(dij_red(adir),zero,KIND=dpc)
           end if
         end do ! iat
       end do ! adir
     end do ! jlmn
   end do ! ilmn
   ABI_FREE(ff)
 end do ! end loop over types
 
end subroutine make_ddir_vhnzc

!!****f* ABINIT/make_ddir_vha3
!! NAME
!! make_ddir_vha3
!!
!! FUNCTION
!! Compute onsite r*vhartree
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
!!
!! SOURCE

subroutine make_ddir_vha3(atindx,ddir_vha,dtset,gntselect,gprimd,lmnmax,my_lmax,pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vha(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,cplex,g1,g2,g3,iat,iatom,ilmn,ij_klm,ij_klmn,ij_kln,imesh,jlmn,jmesh
  integer :: kl_klm,kl_klmn,kl_kln,klmn,ll,llmm,llmn,lm1a,lp,lpmp,itypat,mesh_size,mm,mp
  real(dp) :: c1,c2,nl1,nlt1,rfac,rp,rr,vhaint
  complex(dpc) :: rhokl

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  complex(dpc) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),ff1(:),fft1(:)

!--------------------------------------------------------------------

 c1 = sqrt(four_pi/three)

 ddir_vha = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   cplex = pawrhoij(iatom)%qphase
   mesh_size=pawtab(itypat)%mesh_size
   
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(ff1,(mesh_size))
   ABI_MALLOC(fft1,(mesh_size))

   do ilmn = 1, pawtab(itypat)%lmn_size
     do jlmn = 1, pawtab(itypat)%lmn_size
       ij_klmn = MATPACK(ilmn,jlmn)
       ij_klm = pawtab(itypat)%indklmn(1,ij_klmn)
       ij_kln = pawtab(itypat)%indklmn(2,ij_klmn)

       dij_cart = czero

       do klmn = 1, pawtab(itypat)%lmn_size
         do llmn = 1, pawtab(itypat)%lmn_size
           kl_klmn = MATPACK(klmn,llmn)
           rhokl = CMPLX(pawrhoij(iatom)%rhoij_(2*kl_klmn-1,1),pawrhoij(iatom)%rhoij_(2*kl_klmn,1),KIND=dpc)
           if (klmn .GT. llmn) rhokl = CONJG(rhokl)
           kl_klm = pawtab(itypat)%indklmn(1,kl_klmn)
           kl_kln = pawtab(itypat)%indklmn(2,kl_klmn)

           do ll = pawtab(itypat)%indklmn(3,kl_klmn),pawtab(itypat)%indklmn(4,kl_klmn),2
             c2 = four_pi/(two*ll + one)

             do imesh = 2, mesh_size
               rr = pawrad(itypat)%rad(imesh)
               do jmesh = 2, imesh
                 rp = pawrad(itypat)%rad(jmesh)
                 rfac = (rp**ll)/(rr**(ll+1))
                 ff1(jmesh) = rfac*pawtab(itypat)%phiphj(jmesh,kl_kln)
                 fft1(jmesh) = rfac*pawtab(itypat)%tphitphj(jmesh,kl_kln)
               end do ! jmesh
               do jmesh = imesh+1,mesh_size
                 rp = pawrad(itypat)%rad(jmesh)
                 rfac = (rr**ll)/(rp**(ll+1))
                 ff1(jmesh) = rfac*pawtab(itypat)%phiphj(jmesh,kl_kln)
                 fft1(jmesh) = rfac*pawtab(itypat)%tphitphj(jmesh,kl_kln)
               end do ! jmesh
               call pawrad_deducer0(ff1,mesh_size,pawrad(itypat))
               call simp_gen(nl1,ff1,pawrad(itypat))
               call pawrad_deducer0(fft1,mesh_size,pawrad(itypat))
               call simp_gen(nlt1,fft1,pawrad(itypat))

               ff(imesh) = rr*pawtab(itypat)%phiphj(imesh,ij_kln)*nl1 - &
                 rr*pawtab(itypat)%tphitphj(imesh,ij_kln)

             end do ! imesh
             call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
             call simp_gen(vhaint,ff,pawrad(itypat))

             do mm = -ll, ll
               llmm = LMPACK(ll,mm)
               g3 = gntselect(llmm,kl_klm)
               if (g3 .EQ. 0) cycle
               
               do lp = pawtab(itypat)%indklmn(3,ij_klmn),pawtab(itypat)%indklmn(4,ij_klmn),2
                 do mp = -lp, lp
                   lpmp = LMPACK(lp,mp)
                   g1 = gntselect(lpmp,ij_klm)
                   if (g1 .EQ. 0) cycle

                   do adir = 1, 3
                     lm1a = MATPACK(lpmp,adir_to_sij(adir))
                     g2 = gntselect(llmm,lm1a)
                     if (g2 .EQ. 0) cycle

                     dij_cart(adir) = dij_cart(adir)+rhokl*c1*c2*vhaint*&
                       & realgnt(g1)*realgnt(g2)*realgnt(g3)
                   end do !adir
                 end do !mp
               end do !lp
             end do !mm
           end do !ll
         end do !llmn
       end do !klmn

       dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
       ddir_vha(ilmn,jlmn,iat,1:3) = dij_red(1:3)

     end do !jlmn
   end do !ilmn
  
   ABI_FREE(ff)
   ABI_FREE(ff1)
   ABI_FREE(fft1)

 end do ! end loop over iat
 
end subroutine make_ddir_vha3

!!****f* ABINIT/make_ddir_vnhat
!! NAME
!! make_ddir_vnhat
!!
!! FUNCTION
!! Compute onsite r*vnhat
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
!!
!! SOURCE

subroutine make_ddir_vhnhat(atindx,ddir_vhnhat,dtset,gntselect,gprimd,lmnmax,my_lmax,paw_an,&
    & pawang,pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vhnhat(lmnmax,lmnmax,dtset%natom,3)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,cplex,gint,g2int,iat,iatom,ilm,ilmn,ilmp,imesh,itypat,isel
  integer :: jlmn,jmesh,klmadir,klmn,klm,kln,ll,lm_size,lmax,lmin,lp
  integer :: mesh_size,nzlmopt
  integer :: opt_compch,opt_dens,opt_l,opt_print
  real(dp) :: cdij,compch_sph,nl1,nlt1,rfac,rr,rp,vhaint,rrhokl

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),ff1(:),fft1(:),nhat1(:,:,:)
  real(dp),allocatable :: rho1(:,:,:),trho1(:,:,:)
  logical,allocatable :: lmselectin(:),lmselectout(:)

!--------------------------------------------------------------------

 nzlmopt = -1
 opt_compch = 0
 opt_dens = 0
 opt_l = -1
 opt_print = 0

 cdij = four_pi*sqrt(four_pi/3.0D0)
 ddir_vhnhat = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   cplex = pawrhoij(iatom)%qphase
   mesh_size=pawtab(itypat)%mesh_size
   lm_size = paw_an(iatom)%lm_size
   ABI_MALLOC(lmselectin,(lm_size))
   ABI_MALLOC(lmselectout,(lm_size))
   ABI_MALLOC(nhat1,(cplex*mesh_size,lm_size,dtset%nspden*(1-((opt_dens+1)/2))))
   ABI_MALLOC(rho1,(cplex*mesh_size,lm_size,dtset%nspden))
   ABI_MALLOC(trho1,(cplex*mesh_size,lm_size,dtset%nspden*(1-(opt_dens/2))))
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(ff1,(mesh_size))
   ABI_MALLOC(fft1,(mesh_size))
  
   lmselectin = .TRUE. 
   call pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,&
     & lm_size,nhat1,dtset%nspden,nzlmopt,opt_compch,opt_dens,opt_l,&
     & opt_print,pawang,dtset%pawprtvol,pawrad(itypat),pawrhoij(iatom),&
     & pawtab(itypat),rho1,trho1)

   do klmn = 1, pawtab(itypat)%lmn2_size
     klm = pawtab(itypat)%indklmn(1,klmn)
     kln = pawtab(itypat)%indklmn(2,klmn)
     lmin = pawtab(itypat)%indklmn(3,klmn)
     lmax = pawtab(itypat)%indklmn(4,klmn)
     ilmn = pawtab(itypat)%indklmn(7,klmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn)

     dij_cart = zero

     do ll = lmin,lmax,2
       do ilm = ll**2+1,(ll+1)**2
         gint = gntselect(ilm,klm)
         if (gint .EQ. 0) cycle

         do adir = 1, 3

           do lp = abs(ll-1), ll+1, 2
             do ilmp = lp**2+1,(lp+1)**2
               klmadir = MATPACK(ilm,adir_to_sij(adir))
               g2int = gntselect(ilmp,klmadir)
               if (g2int .EQ. 0) cycle

               ! construct integrand for rho1 vHa, trho1
               do imesh = 2, mesh_size
                 rr = pawrad(itypat)%rad(imesh)
                   
                 ! for this mesh point, do the interior nonlocal integral over Hartree potential
                 do jmesh = 2, imesh
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rp**lp)/(rr**(lp+1))
                   ff1(jmesh) = nhat1(jmesh,ilmp,1)*rfac
                 end do
                 do jmesh=imesh+1, mesh_size
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rr**lp)/(rp**(lp+1))
                   ff1(jmesh) = nhat1(jmesh,ilmp,1)*rfac
                 end do

                 call pawrad_deducer0(ff1,mesh_size,pawrad(itypat))
                 call simp_gen(nl1,ff1,pawrad(itypat))
                 
                 ff(imesh)=-rr*pawtab(itypat)%tphitphj(imesh,kln)*nl1
               
               end do ! end loop over imesh

               call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
               call simp_gen(vhaint,ff,pawrad(itypat))

               dij_cart(adir) = dij_cart(adir) + &
                 & cdij*realgnt(gint)*realgnt(g2int)*vhaint/(two*ll+one)
               
             end do ! end loop over ilmp
           end do ! end loop over lp
         end do ! end loop over adir
       end do ! end loop over ilm
     end do ! end loop over ll

     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

     do adir = 1, 3
       ddir_vhnhat(ilmn,jlmn,iat,adir) = CMPLX(dij_red(adir),zero,KIND=dpc)
       ddir_vhnhat(jlmn,ilmn,iat,adir) = ddir_vhnhat(ilmn,jlmn,iat,adir)
     end do

   end do ! end loop over klmn

   ABI_FREE(lmselectin)
   ABI_FREE(lmselectout)
   ABI_FREE(nhat1)
   ABI_FREE(rho1)
   ABI_FREE(trho1)
   ABI_FREE(ff)
   ABI_FREE(ff1)
   ABI_FREE(fft1)

 end do ! end loop over iat
 
end subroutine make_ddir_vhnhat

!!****f* ABINIT/make_ddir_vxc
!! NAME
!! make_ddir_vxc
!!
!! FUNCTION
!! Compute onsite r*vxc
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
!!
!! SOURCE

subroutine make_ddir_vxc(atindx,ddir_vxc,dtset,gntselect,gprimd,lmnmax,my_lmax,paw_an,&
    & pawang,pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vxc(lmnmax,lmnmax,dtset%natom,3)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,cplex,gint,g2int,iat,iatom,ilm,ilmn,ilmp,imesh,itypat,isel
  integer :: jlmn,jmesh,klmadir,klmn,klm,kln,ll,lm_size,lmin,lmax,lp
  integer :: mesh_size,nzlmopt
  integer :: opt_compch,opt_dens,opt_l,opt_print
  real(dp) :: cdij,compch_sph,nl1,nlt1,rfac,rr,rp,vhaint,rrhokl
  real(dp) :: nterm

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: nhat1(:,:,:),rho1(:,:,:),trho1(:,:,:)
  logical,allocatable :: lmselectin(:),lmselectout(:)

!--------------------------------------------------------------------

 nzlmopt = -1
 opt_compch = 0
 opt_dens = 1
 opt_l = -1
 opt_print = 0

 cdij = sqrt(four_pi/3.0D0)
 ddir_vxc = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   cplex = pawrhoij(iatom)%qphase
   mesh_size=pawtab(itypat)%mesh_size
   lm_size = paw_an(iatom)%lm_size
   ABI_MALLOC(lmselectin,(lm_size))
   ABI_MALLOC(lmselectout,(lm_size))
   ABI_MALLOC(nhat1,(cplex*mesh_size,lm_size,dtset%nspden*(1-((opt_dens+1)/2))))
   ABI_MALLOC(rho1,(cplex*mesh_size,lm_size,dtset%nspden))
   ABI_MALLOC(trho1,(cplex*mesh_size,lm_size,dtset%nspden*(1-(opt_dens/2))))
  
   lmselectin = .TRUE. 
   call pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,&
     & lm_size,nhat1,dtset%nspden,nzlmopt,opt_compch,opt_dens,opt_l,&
     & opt_print,pawang,dtset%pawprtvol,pawrad(itypat),pawrhoij(iatom),&
     & pawtab(itypat),rho1,trho1)

   ABI_FREE(lmselectin)
   ABI_FREE(lmselectout)
   ABI_FREE(nhat1)
   ABI_FREE(rho1)
   ABI_FREE(trho1)

 end do ! end loop over iat
 
end subroutine make_ddir_vxc


!!****f* ABINIT/make_ddir_vha2
!! NAME
!! make_ddir_vha2
!!
!! FUNCTION
!! Compute onsite r*vhartree
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
!!
!! SOURCE

subroutine make_ddir_vha2(atindx,ddir_vha,dtset,gntselect,gprimd,lmnmax,my_lmax,paw_an,&
    & pawang,pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vha(lmnmax,lmnmax,dtset%natom,3)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,cplex,gint,g2int,iat,iatom,ilm,ilmn,ilmp,imesh,itypat,isel
  integer :: jlmn,jmesh,klmadir,klmn,klm,kln,ll,lm_size,lmin,lmax,lp
  integer :: mesh_size,nzlmopt
  integer :: opt_compch,opt_dens,opt_l,opt_print
  real(dp) :: cdij,compch_sph,nl1,nlt1,rfac,rr,rp,vhaint,rrhokl
  real(dp) :: nterm

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),ff1(:),fft1(:),nhat1(:,:,:)
  real(dp),allocatable :: my_rho1(:,:,:),rho1(:,:,:),trho1(:,:,:)
  logical,allocatable :: lmselectin(:),lmselectout(:)

!--------------------------------------------------------------------

 nzlmopt = -1
 opt_compch = 0
 opt_dens = 1
 opt_l = -1
 opt_print = 0

 cdij = four_pi*sqrt(four_pi/3.0D0)
 ddir_vha = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   cplex = pawrhoij(iatom)%qphase
   mesh_size=pawtab(itypat)%mesh_size
   lm_size = paw_an(iatom)%lm_size
   ABI_MALLOC(lmselectin,(lm_size))
   ABI_MALLOC(lmselectout,(lm_size))
   ABI_MALLOC(nhat1,(cplex*mesh_size,lm_size,dtset%nspden*(1-((opt_dens+1)/2))))
   ABI_MALLOC(rho1,(cplex*mesh_size,lm_size,dtset%nspden))
   !ABI_MALLOC(my_rho1,(cplex*mesh_size,lm_size,dtset%nspden))
   ABI_MALLOC(trho1,(cplex*mesh_size,lm_size,dtset%nspden*(1-(opt_dens/2))))
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(ff1,(mesh_size))
   ABI_MALLOC(fft1,(mesh_size))
  
   lmselectin = .TRUE. 
   call pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,&
     & lm_size,nhat1,dtset%nspden,nzlmopt,opt_compch,opt_dens,opt_l,&
     & opt_print,pawang,dtset%pawprtvol,pawrad(itypat),pawrhoij(iatom),&
     & pawtab(itypat),rho1,trho1)

   !!! check rho1 density
   !write(std_out,'(a)') 'JWZ debug call my pawdensities '
   !my_rho1 = zero
   !do isel = 1, pawrhoij(iatom)%nrhoijsel
   !  klmn = pawrhoij(iatom)%rhoijselect(isel)
   !  klm = pawtab(itypat)%indklmn(1,klmn)
   !  kln = pawtab(itypat)%indklmn(2,klmn)
   !  lmin = pawtab(itypat)%indklmn(3,klmn)
   !  lmax = pawtab(itypat)%indklmn(4,klmn)
   !  rrhokl = pawrhoij(iatom)%rhoijp(2*isel-1,1)
   !  rrhokl = rrhokl*pawtab(itypat)%dltij(klmn)
   !  do ll=lmin,lmax,2
   !    do ilm = ll**2+1,(ll+1)**2
   !      gint=gntselect(ilm,klm)
   !      if (gint .NE. 0) then
   !        do imesh = 2, mesh_size
   !          nterm = rrhokl*pawtab(itypat)%phiphj(imesh,kln)*realgnt(gint)/&
   !            & (pawrad(itypat)%rad(imesh)**2)
   !          my_rho1(imesh,ilm,1) = my_rho1(imesh,ilm,1) + nterm
   !        end do !imesh
   !      end if ! gint .NE. 0
   !    end do !ilm
   !  end do ! ll
   !end do !isel
   !
   !do ll=0,pawtab(itypat)%lcut_size-1
   !  do ilm = ll**2+1,(ll+1)**2
   !    if (lmselectout(ilm)) then
   !      call pawrad_deducer0(my_rho1(:,ilm,1),mesh_size,pawrad(itypat))
   !      do imesh = 1, 20
   !        write(std_out,'(a,2i4,2es16.8)')'JWZ debug ilm imesh myrho1 rho1 ',ilm,&
   !          & imesh,my_rho1(imesh,ilm,1),rho1(imesh,ilm,1)
   !      end do !imesh
   !    end if
   !  end do ! ilm
   !end do ! ll
   
   do klmn = 1, pawtab(itypat)%lmn2_size
     klm = pawtab(itypat)%indklmn(1,klmn)
     kln = pawtab(itypat)%indklmn(2,klmn)
     lmin = pawtab(itypat)%indklmn(3,klmn)
     lmax = pawtab(itypat)%indklmn(4,klmn)
     ilmn = pawtab(itypat)%indklmn(7,klmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn)

     dij_cart = zero

     do ll = lmin,lmax,2
       do ilm = ll**2+1,(ll+1)**2
         gint = gntselect(ilm,klm)
         if (gint .EQ. 0) cycle

         do adir = 1, 3

           do lp = abs(ll-1), ll+1, 2
             do ilmp=lp**2+1,(lp+1)**2
               klmadir = MATPACK(ilm,adir_to_sij(adir))
               g2int = gntselect(ilmp,klmadir)
               if (g2int .EQ. 0) cycle

               ! construct integrand for rho1 vHa, trho1
               ff = zero
               do imesh = 2, mesh_size
                 rr = pawrad(itypat)%rad(imesh)
                   
                 ! for this mesh point, do the interior nonlocal integral over Hartree potential
                 ff1 = zero; fft1 = zero
                 do jmesh = 2, imesh
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rp**lp)/(rr**(lp+1))
                   ff1(jmesh) = rho1(jmesh,ilmp,1)*rfac
                   fft1(jmesh) = trho1(jmesh,ilmp,1)*rfac
                 end do
                 do jmesh=imesh+1, mesh_size
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rr**lp)/(rp**(lp+1))
                   ff1(jmesh) = rho1(jmesh,ilmp,1)*rfac
                   fft1(jmesh) = trho1(jmesh,ilmp,1)*rfac
                 end do

                 call pawrad_deducer0(ff1,mesh_size,pawrad(itypat))
                 call simp_gen(nl1,ff1,pawrad(itypat))
                 call pawrad_deducer0(fft1,mesh_size,pawrad(itypat))
                 call simp_gen(nlt1,fft1,pawrad(itypat))
                 
                 ff(imesh)=rr*pawtab(itypat)%phiphj(imesh,kln)*nl1 - &
                   & rr*pawtab(itypat)%tphitphj(imesh,kln)*nlt1
               
               end do ! end loop over imesh

               call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
               call simp_gen(vhaint,ff,pawrad(itypat))

               dij_cart(adir) = dij_cart(adir) + &
                 & cdij*realgnt(gint)*realgnt(g2int)*vhaint/(two*ll+one)
               
             end do ! end loop over ilmp
           end do ! end loop over lp
         end do ! end loop over adir
       end do ! end loop over ilm
     end do ! end loop over ll

     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

     do adir = 1, 3
       ddir_vha(ilmn,jlmn,iat,adir) = CMPLX(dij_red(adir),zero,KIND=dpc)
       ddir_vha(jlmn,ilmn,iat,adir) = ddir_vha(ilmn,jlmn,iat,adir)
     end do

   end do ! end loop over klmn

   ABI_FREE(lmselectin)
   ABI_FREE(lmselectout)
   ABI_FREE(nhat1)
   !ABI_FREE(my_rho1)
   ABI_FREE(rho1)
   ABI_FREE(trho1)
   ABI_FREE(ff)
   ABI_FREE(ff1)
   ABI_FREE(fft1)

 end do ! end loop over iat
 
end subroutine make_ddir_vha2

!!****f* ABINIT/make_ddir_vha
!! NAME
!! make_ddir_vha
!!
!! FUNCTION
!! Compute onsite r*vhartree
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

subroutine make_ddir_vha(atindx,ddir_vha,dtset,gntselect,gprimd,lmnmax,lmn2max,my_lmax,&
    & pawrad,pawrhoij,pawtab,psps,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,lmn2max,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vha(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bigl,biglm,bigm,gnt1,gnt2,gnt3,iat,iatom,il,ilm,ilmn,im,ir,itypat
  integer :: jlm,jlmn,klmn,kln
  integer :: ll,llmm,llmmjlm,lm1b,max_ij_size,mm,msz,olmn,olmni,olmnj,olm,oln,psz
  real(dp) :: cdij,cl,intg1,rg1,rg2,rg3,rr
  complex(dpc) :: rhomn

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:),fmnl(:),intg(:,:,:,:),tfmnl(:)
  complex(dpc) :: dij_cart(3),dij_red(3)

!--------------------------------------------------------------------

 ! dij prefactor
 cdij = sqrt(four_pi/three)

 max_ij_size = psps%lnmax*(psps%lnmax+1)/2 
 ! radial integrals
 ABI_MALLOC(intg,(max_ij_size,max_ij_size,2*psps%mpsang,dtset%ntypat))
 do itypat = 1, dtset%ntypat
   msz = pawtab(itypat)%mesh_size
   psz = size(pawtab(itypat)%phiphj(:,1))
   ABI_MALLOC(ff,(msz))
   ABI_MALLOC(fmnl,(msz))
   ABI_MALLOC(tfmnl,(msz))
   do oln = 1, pawtab(itypat)%ij_size
     do bigl = 0, pawtab(itypat)%l_size - 1
       do ir=2,psz

         ff=zero
         rr=pawrad(itypat)%rad(ir)
         ff(2:ir) = pawtab(itypat)%phiphj(2:ir,oln)*pawrad(itypat)%rad(2:ir)**bigl/rr**(bigl+1)
         ff(ir:psz) = pawtab(itypat)%phiphj(ir:psz,oln)*rr**bigl/pawrad(itypat)%rad(ir:psz)**(bigl+1)
         call pawrad_deducer0(ff,msz,pawrad(itypat))
         call simp_gen(intg1,ff,pawrad(itypat))
         fmnl(ir) = intg1

         ff=zero
         ff(2:ir) = pawtab(itypat)%tphitphj(2:ir,oln)*pawrad(itypat)%rad(2:ir)**bigl/rr**(bigl+1)
         ff(ir:psz) = pawtab(itypat)%tphitphj(ir:psz,oln)*rr**bigl/pawrad(itypat)%rad(ir:psz)**(bigl+1)
         call pawrad_deducer0(ff,msz,pawrad(itypat))
         call simp_gen(intg1,ff,pawrad(itypat))
         tfmnl(ir) = intg1

       end do ! end loop over ir

       do kln = 1, pawtab(itypat)%ij_size
         ff=zero
         ff(2:psz) = pawrad(itypat)%rad(2:psz)*pawtab(itypat)%phiphj(2:psz,kln)*fmnl(2:psz)
         ff(2:psz) = ff(2:psz) - pawrad(itypat)%rad(2:psz)*pawtab(itypat)%tphitphj(2:psz,kln)*tfmnl(2:psz)
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

 ddir_vha = czero
 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   do ilmn = 1, pawtab(itypat)%lmn_size
     do jlmn = 1, pawtab(itypat)%lmn_size
       klmn = MATPACK(ilmn,jlmn)
       kln = pawtab(itypat)%indklmn(2,klmn)
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
               olmni = pawtab(itypat)%indklmn(7,olmn)
               olmnj = pawtab(itypat)%indklmn(8,olmn)
               olm = pawtab(itypat)%indklmn(1,olmn)
               oln = pawtab(itypat)%indklmn(2,olmn)
               rhomn=CMPLX(pawrhoij(iatom)%rhoij_(2*olmn-1,1),pawrhoij(iatom)%rhoij_(2*olmn,1),KIND=dpc)
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
                   dij_cart(adir) = dij_cart(adir) + cdij*cl*rhomn*rg1*rg2*rg3*intg(kln,oln,bigl+1,itypat)
                   if(olmni /= olmnj) then
                     dij_cart(adir) = dij_cart(adir) + cdij*cl*CONJG(rhomn)*rg1*rg2*rg3*intg(kln,oln,bigl+1,itypat)
                   end if ! end loop over bigm
                 end do ! end bigm
               end do ! end loop over bigl
             end do ! end loop over olmn
           end do ! end loop over mm
         end do ! end loop over ll
       end do ! end loop over adir
         
       dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

       ddir_vha(ilmn,jlmn,iat,1:3) = dij_red(1:3)
     end do ! jlmn 
   end do ! ilmn
 end do ! end loop over iat

 ABI_FREE(intg) 
end subroutine make_ddir_vha

!!****f* ABINIT/make_ddir_p2
!! NAME
!! make_ddir_p2
!!
!! FUNCTION
!! Compute onsite r*p^2/2
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

subroutine make_ddir_p2(ddir_p2,dtset,gntselect,gprimd,lmnmax,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_p2(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,iat,itypat,ilmn,iln,jl,jlmn,jln,klm,klmn,lm1b,mesh_size,pwave_size
  real(dp) :: c1m,intg

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:),uj(:),ujder(:),uj2der(:)
  complex(dpc) :: dij_cart(3)

!--------------------------------------------------------------------

  c1m = sqrt(four_pi/three)
  ddir_p2=czero
  do iat=1,dtset%natom
    itypat=dtset%typat(iat)
    mesh_size=pawrad(itypat)%mesh_size
    pwave_size=size(pawtab(itypat)%phiphj(:,1))
    ABI_MALLOC(ff,(mesh_size))
    ABI_MALLOC(uj,(mesh_size))
    ABI_MALLOC(ujder,(mesh_size))
    ABI_MALLOC(uj2der,(mesh_size))

    do ilmn=1, pawtab(itypat)%lmn_size
      do jlmn=1, pawtab(itypat)%lmn_size
        klmn = MATPACK(ilmn,jlmn)

        klm = pawtab(itypat)%indklmn(1,klmn)
        iln = pawtab(itypat)%indlmn(5,ilmn)
        jln = pawtab(itypat)%indlmn(5,jlmn)
        jl = pawtab(itypat)%indlmn(1,jlmn)

        ff = zero
        uj = zero
        uj(1:pwave_size) = pawtab(itypat)%phi(1:pwave_size,jln)
        call nderiv_gen(ujder,uj,pawrad(itypat),uj2der)
        ff(2:pwave_size) = pawrad(itypat)%rad(2:pwave_size)*pawtab(itypat)%phi(2:pwave_size,iln)*(uj2der(2:pwave_size)&
          & - jl*(jl+1)*pawtab(itypat)%phi(2:pwave_size,jln)/pawrad(itypat)%rad(2:pwave_size)**2)
        uj = zero
        uj(1:pwave_size) = pawtab(itypat)%tphi(1:pwave_size,jln)
        call nderiv_gen(ujder,uj,pawrad(itypat),uj2der)
        ff(2:pwave_size) = ff(2:pwave_size) - &
          & (pawrad(itypat)%rad(2:pwave_size)*pawtab(itypat)%tphi(2:pwave_size,iln)*(uj2der(2:pwave_size)&
          & - jl*(jl+1)*pawtab(itypat)%tphi(2:pwave_size,jln)/pawrad(itypat)%rad(2:pwave_size)**2))
        call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
        call simp_gen(intg,ff,pawrad(itypat))

        dij_cart = czero 
        do adir = 1, 3
          lm1b = adir_to_sij(adir)
          gint = gntselect(lm1b,klm)
          if (gint == 0) cycle
          dij_cart(adir) = cmplx(-half*c1m*realgnt(gint)*intg,zero)
        end do ! end loop over adir

        ddir_p2(ilmn,jlmn,iat,1:3) = MATMUL(TRANSPOSE(gprimd),dij_cart)
      end do ! jlmn
    end do ! ilmn
    ABI_FREE(ff)
    ABI_FREE(uj)
    ABI_FREE(ujder)
    ABI_FREE(uj2der)
  end do ! end loop over atoms
 
end subroutine make_ddir_p2
!!***

!!****f* ABINIT/make_ddir_ap
!! NAME
!! make_ddir_ap
!!
!! FUNCTION
!! Compute onsite r*A.p 
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

subroutine make_ddir_ap(ddir_ap,dtset,gntselect,gprimd,lmnmax,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_ap(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bigl,biglm,bigm,gint,iat,ilmn,il,ilm,ilml1b,itypat
  integer :: jl,jlmn,jm,klmn,kln,l1b,ldir,mesh_size,pwave_size
  real(dp) :: a2,intg
  complex(dpc) :: cme,orbl_me

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp),allocatable :: ff(:)
  complex(dpc) :: dij_cart(3)

!--------------------------------------------------------------------

  a2 = FineStructureConstant2

  ddir_ap=czero
  do iat=1,dtset%natom
    itypat=dtset%typat(iat)
    mesh_size=pawtab(itypat)%mesh_size
    pwave_size=size(pawtab(itypat)%phiphj(:,1))
    ABI_MALLOC(ff,(mesh_size))
    do ilmn = 1, pawtab(itypat)%lmn_size
      do jlmn = 1, pawtab(itypat)%lmn_size
        klmn = MATPACK(ilmn,jlmn)
        kln = pawtab(itypat)%indklmn(2,klmn) 

        ff=zero

        ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r^2
        ff(2:pwave_size) = (pawtab(itypat)%phiphj(2:pwave_size,kln)-&
          &                pawtab(itypat)%tphitphj(2:pwave_size,kln))/&
          &               pawrad(itypat)%rad(2:pwave_size)**2
        call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
        call simp_gen(intg,ff,pawrad(itypat))
 
        ilm=pawtab(itypat)%indklmn(5,klmn)
        il=pawtab(itypat)%indlmn(1,ilmn)

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
        ddir_ap(ilmn,jlmn,iat,1:3) = MATMUL(TRANSPOSE(gprimd),dij_cart)
      end do ! jlmn
    end do ! ilmn
    ABI_FREE(ff)
  end do ! end loop over atoms
 
end subroutine make_ddir_ap
!!***

!!****f* ABINIT/pack_pawrhoij
!! NAME
!! pack_pawrhoij
!!
!! FUNCTION
!! for a pawhoij data structure, add packed data from unpacked
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

subroutine pack_pawrhoij(dtset,pawrhoij)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),intent(in) :: dtset
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom)

  !arrays

  !Local variables -------------------------
  !scalars
  integer :: iatom,isel,klmn
  real(dp) :: ri,rr
  logical :: cprho
  character(len=500) :: msg
 
  !arrays
!--------------------------------------------------------------------

  do iatom = 1, dtset%natom
    if (pawrhoij(iatom)%use_rhoij_ .EQ. 0) then
      msg='Unpacked rhoij elements not available.'
      ABI_BUG(msg)
    end if

    if(allocated(pawrhoij(iatom)%rhoijselect)) then
      deallocate(pawrhoij(iatom)%rhoijselect)
    end if
    if(allocated(pawrhoij(iatom)%rhoijp)) then
      deallocate(pawrhoij(iatom)%rhoijp)
    end if

    ! count number of nonzero pawrhoij elements
    isel = 0
    cprho = (pawrhoij(iatom)%cplex_rhoij .EQ. 2)
    do klmn = 1, pawrhoij(iatom)%lmn2_size
      if (cprho) then
        rr = pawrhoij(iatom)%rhoij_(2*klmn-1,1)
        ri = pawrhoij(iatom)%rhoij_(2*klmn,1)
      else
        rr = pawrhoij(iatom)%rhoij_(klmn,1)
        ri = zero
      end if
      !if ( (ABS(rr) .GT. tol16) .OR. (ABS(ri) .GT. tol16) ) then
      !  isel = isel + 1
      !end if
      isel = isel+1
    end do ! end loop over klmn
    pawrhoij(iatom)%nrhoijsel = isel

    if (cprho) then
      allocate(pawrhoij(iatom)%rhoijp(1:2*isel,dtset%nspden))
    else
      allocate(pawrhoij(iatom)%rhoijp(1:isel,dtset%nspden))
    end if
    allocate(pawrhoij(iatom)%rhoijselect(1:isel))

    pawrhoij(iatom)%rhoijp = zero
    pawrhoij(iatom)%rhoijselect = 0

    isel = 0 
    do klmn = 1, pawrhoij(iatom)%lmn2_size
      if (cprho) then
        rr = pawrhoij(iatom)%rhoij_(2*klmn-1,1)
        ri = pawrhoij(iatom)%rhoij_(2*klmn,1)
      else
        rr = pawrhoij(iatom)%rhoij_(klmn,1)
        ri = zero
      end if
      !if ( (ABS(rr) .GT. tol16) .OR. (ABS(ri) .GT. tol16) ) then
        isel = isel + 1
        pawrhoij(iatom)%rhoijselect(isel) = klmn
        if (cprho) then
          pawrhoij(iatom)%rhoijp(2*isel-1,1) = rr
          pawrhoij(iatom)%rhoijp(2*isel,1) = ri
        else
          pawrhoij(iatom)%rhoijp(isel,1) = rr
        end if ! if cprho
      !end if ! if rhoij nonzero
    end do ! end loop over klmn

    pawrhoij(iatom)%use_rhoijp=1

  end do ! end loop over iatom
 
end subroutine pack_pawrhoij
!!***

!!****f* ABINIT/make_d
!! NAME
!! make_d
!!
!! FUNCTION
!! this is a driver to compute all the different d terms
!! <phi_i|d*V|phi_j> - <tphi_i|d*V|tphi_j>
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

subroutine make_d(atindx,atindx1,cprj,dimlmn,dterms,dtset,&
    & gprimd,mcprj,mpi_enreg,occ,paw_an,pawang,pawrad,pawtab,psps)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcprj
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
  integer,intent(in) :: dimlmn(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  complex(dpc),intent(out) :: dterms(ndterms,psps%lmnmax,psps%lmnmax,dtset%natom,3)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iband,isel,itypat,ilmn,jlmn,klmn
  integer :: my_lmax,ngnt
  complex(dpc) :: cpi,cpj,rhoij
  type(paw_dmft_type) :: paw_dmft
 
  !arrays
  integer,allocatable :: gntselect(:,:)
  real(dp),allocatable :: realgnt(:)
  type(pawrhoij_type),allocatable :: pawrhoij(:)
!--------------------------------------------------------------------

 !make pawrhoij
 ABI_MALLOC(pawrhoij,(dtset%natom))
 call pawrhoij_alloc(pawrhoij,dtset%pawcpxocc,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%typat,&
   & pawtab=pawtab)
 paw_dmft%use_sc_dmft=0
 paw_dmft%use_dmft=0
 call pawmkrhoij(atindx,atindx1,cprj,dimlmn,dtset%istwfk,dtset%kptopt,dtset%mband,dtset%mband,&
   & mcprj,dtset%mkmem,mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
   & occ,dtset%paral_kgb,paw_dmft,pawrhoij,0,dtset%usewvl,dtset%wtk)
 call pack_pawrhoij(dtset,pawrhoij)

 !do iat = 1, dtset%natom
 !  itypat = dtset%typat(iat)
 !  write(std_out,'(a,2i4)')'JWZ debug pawrhoij cplex_rhoij qphase ',pawrhoij(iat)%cplex_rhoij,pawrhoij(iat)%qphase
 !  do klmn = 1, pawrhoij(iat)%lmn2_size
 !    ilmn = pawtab(itypat)%indklmn(7,klmn)
 !    jlmn = pawtab(itypat)%indklmn(8,klmn)
 !    rhoij = czero
 !    do iband = 1, dtset%mband
 !      cpi = CMPLX(cprj(1,iband)%cp(1,ilmn),cprj(1,iband)%cp(2,ilmn))
 !      cpj = CMPLX(cprj(1,iband)%cp(1,jlmn),cprj(1,iband)%cp(2,jlmn))
 !      rhoij = rhoij + two*CONJG(cpi)*cpj
 !    end do
 !    do isel = 1, pawrhoij(iat)%nrhoijsel
 !      if (pawrhoij(iat)%rhoijselect(isel) .EQ. klmn) then
 !        write(std_out,'(a,4es16.8)')' JWZ debug pawrhoij rhoij ',pawrhoij(iat)%rhoijp(2*isel-1,1),pawrhoij(iat)%rhoijp(2*isel,1),&
 !          & real(rhoij),aimag(rhoij)
 !      end if
 !    end do
 !  end do
 !end do
 
 ! make Gaunt integrals
 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 ! generate d terms

 ! term idsij arises from <phi|r|phi> - <tphi|r|rphi>, is the d analog of pawtab(itypat)%sij
 call make_ddir_sij(dterms(idsij,:,:,:,:),dtset,gprimd,psps%lmnmax,pawtab)

 ! term idp2 due to onsite p^2/2, corresponds to term 1 of Torrent PAW roadmap paper appendix E
 ! Comp. Mat. Sci. 42, 337-351 (2008)
 call make_ddir_p2(dterms(idp2,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)

 ! term idpa due to onsite A.p, not part of the roadmap paper but similar to term 1
 call make_ddir_ap(dterms(idpa,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)
 
 ! term idvha due to v_H[n1], corresponds to term 2a of roadmap paper
 call make_ddir_vha2(atindx,dterms(idvha,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,paw_an,&
   & pawang,pawrad,pawrhoij,pawtab,realgnt)

 ! term idvhnzc due to v_H[nZ], corresponds to term 2b of roadmap paper
 call make_ddir_vhnzc(dterms(idvhnzc,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)

 ! term idvhnhat due to v_H[nhat], corresponds to term 2c of roadmap paper
 call make_ddir_vhnhat(atindx,dterms(idvhnhat,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,paw_an,&
   & pawang,pawrad,pawrhoij,pawtab,realgnt)

 call pawrhoij_free(pawrhoij)
 ABI_FREE(pawrhoij)
 ABI_FREE(realgnt)
 ABI_FREE(gntselect)
 
end subroutine make_d
!!***


end module m_orbmag
