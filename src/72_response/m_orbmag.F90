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

  integer,parameter :: ibcc=1,ibcv=2,ibvva=3,ibvvb=4
  integer,parameter :: iomlr=5,iomanp=6,iomlmb=7
  integer,parameter :: imcch=8,imvvbh=9,imcce=10,imcve=11,imvvae=12,imvvbe=13
  integer,parameter :: nterms=13
  real(dp),parameter :: c2 = one/(two_pi*two_pi) ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  complex(dpc),parameter :: cbc = j_dpc/two_pi ! Berry curvature pre-factor
  complex(dpc),parameter :: com = -half*j_dpc  ! Orbital magnetism pre-factor

  ! Bound methods:

  public :: orbmag_ptpaw

  private :: make_ddir
  private :: cc_fgh
  private :: cv_fh
  private :: vva_fh
  private :: vvb_fgh

  private :: d2lr_p2
  private :: d2lr_Anp

  private :: tdt_me
  private :: apply_d2lr_term_k
  private :: lamb_core
  private :: make_pcg1
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
    & nfftf,ngfftf,npwarr,paw_ij,pawfgr,pawrad,pawtab,psps,rprimd,vtrial,xred,ylm,ylmgr)

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
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
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
 integer :: me,mcgk,my_lmax,my_nspinor,nband_k,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:),ddir(:,:,:)
 real(dp),allocatable :: eig_k(:),ff_k(:,:,:),ffnl_k(:,:,:,:),gg_k(:,:,:),hh_k(:,:,:) 
 real(dp),allocatable :: kinpw(:),kpg_k(:,:),orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 complex(dp),allocatable :: mp2(:,:,:),mpan(:,:,:)
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

 !ABI_MALLOC(rhoij,(2,lmn2max,dtset%natom))
 !write(std_out,'(a)')' Calling local_rhoij...'
 !call local_rhoij(atindx,dtset,cprj,lmn2max,mcprj,mpi_enreg,rhoij,pawtab,ucvol)

 !ABI_MALLOC(dlq,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling dl_q...'
 !call dl_q(dlq,dtset,gprimd,lmn2max,pawtab)

 ABI_MALLOC(ddir,(lmn2max,dtset%ntypat,3))
 call make_ddir(ddir,dtset,gprimd,lmn2max,pawtab)

 ABI_MALLOC(mp2,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling d2lr_p2...'
 call d2lr_p2(dtset,gprimd,lmn2max,mp2,pawrad,pawtab)

 ABI_MALLOC(mpan,(lmn2max,dtset%natom,3))
 write(std_out,'(a)')' Calling d2lr_Anp...'
 call d2lr_Anp(dtset,gntselect,gprimd,lmn2max,mpan,my_lmax,pawrad,pawtab,realgnt)

 !ABI_MALLOC(dlanp,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling dl_Anp...'
 !call dl_Anp(dlanp,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)
 !
 !ABI_MALLOC(dlvhnzc,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling dl_vhnzc...'
 !call dl_vhnzc(dlvhnzc,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)
 !
 !ABI_MALLOC(dlvh,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling dl_vha...'
 !call dl_vha(atindx,dlvh,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,psps,realgnt,rhoij)
 !
 !ABI_MALLOC(dlp2,(lmn2max,dtset%natom,3))
 !write(std_out,'(a)')' Calling dl_p2...'
 !call dl_p2(dlp2,dtset,gntselect,gprimd,lmn2max,my_lmax,pawrad,pawtab,realgnt)

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

   call cc_fgh(atindx,cprj1_k,dimlmn,dtset,eig_k,ff_k,gg_k,gs_hamk,hh_k,ikpt,isppol,mcgk,&
     & mpi_enreg,my_nspinor,nband_k,npw_k,pcg1_k)
   orbmag_terms(:,:,:,ibcc) = orbmag_terms(:,:,:,ibcc) + trnrm*ff_k
   orbmag_terms(:,:,:,imcce) = orbmag_terms(:,:,:,imcce) + trnrm*hh_k
   orbmag_terms(:,:,:,imcch) = orbmag_terms(:,:,:,imcch) + trnrm*gg_k
  
   call vva_fh(atindx,cprj_k,ddir,dtset,eig_k,ff_k,hh_k,lmn2max,nband_k,pawtab)
   orbmag_terms(:,:,:,ibvva) = orbmag_terms(:,:,:,ibvva) + trnrm*ff_k
   orbmag_terms(:,:,:,imvvae) = orbmag_terms(:,:,:,imvvae) + trnrm*hh_k
   
   call vvb_fgh(atindx,cprj_k,ddir,dtset,eig_k,ff_k,gg_k,hh_k,lmn2max,nband_k,pawtab)
   orbmag_terms(:,:,:,ibvvb) = orbmag_terms(:,:,:,ibvvb) + trnrm*ff_k
   orbmag_terms(:,:,:,imvvbe) = orbmag_terms(:,:,:,imvvbe) + trnrm*hh_k
   orbmag_terms(:,:,:,imvvbh) = orbmag_terms(:,:,:,imvvbh) + trnrm*gg_k
   
   call cv_fh(atindx,cprj_k,cprj1_k,ddir,dtset,eig_k,ff_k,hh_k,lmn2max,nband_k,pawtab)
   orbmag_terms(:,:,:,ibcv) = orbmag_terms(:,:,:,ibcv) + trnrm*ff_k
   orbmag_terms(:,:,:,imcve) = orbmag_terms(:,:,:,imcve) + trnrm*hh_k
   
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

 ABI_FREE(ddir)
 ABI_FREE(mp2)
 ABI_FREE(mpan)
 
 !ABI_FREE(rhoij)
 !ABI_FREE(dlq)
 !ABI_FREE(mp2)
 !ABI_FREE(mpan)
 !ABI_FREE(dlanp)
 !ABI_FREE(dlvh)
 !ABI_FREE(dlvhnzc)
 !ABI_FREE(dlp2)

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
    & ket_n,lmn2max,nband_k,pawtab,tdtme)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bra_n,idir,ket_n,lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: ddir(lmn2max,dtset%ntypat,3)
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
          & CONJG(cpi)*(-j_dpc)*ddir(klmn,itypat,idir)*cpj

      end do ! jlmn
    end do ! ilmn
  end do ! iatom

end subroutine tdt_me
!!***


!!****f* ABINIT/cv_fh
!! NAME
!! cv_fh
!!
!! FUNCTION
!! compute CV contribution to f and h functions
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

subroutine cv_fh(atindx,cprj_k,cprj1_k,ddir,dtset,eig_k,ff_k,hh_k,lmn2max,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: ddir(lmn2max,dtset%ntypat,3),eig_k(nband_k)
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
              & nn,lmn2max,nband_k,pawtab,ct2)
          call tdt_me(atindx,cprj1_k(:,:,gdir),nn,ddir,dtset,bdir,cprj_k,&
              & nn,lmn2max,nband_k,pawtab,ct3)

          ct = c2*cbc*eabg*(ct2+CONJG(ct3))
          ff_k(1,nn,adir) = ff_k(1,nn,adir) + REAL(ct)
          ff_k(2,nn,adir) = ff_k(2,nn,adir) + AIMAG(ct)
          
          ct = c2*com*eabg*eig_k(nn)*(ct2+CONJG(ct3))
          hh_k(1,nn,adir) = hh_k(1,nn,adir) + REAL(ct)
          hh_k(2,nn,adir) = hh_k(2,nn,adir) + AIMAG(ct)
        
        end do ! gdir
      end do ! bdir
    end do ! adir

  end do !loop over bands

end subroutine cv_fh
!!***


!!****f* ABINIT/vva_fh
!! NAME
!! vva_fh
!!
!! FUNCTION
!! compute VVa contribution to f and h functions
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

subroutine vva_fh(atindx,cprj_k,ddir,dtset,eig_k,ff,hh,lmn2max,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: ddir(lmn2max,dtset%ntypat,3),eig_k(nband_k)
  real(dp),intent(out) :: ff(2,nband_k,3),hh(2,nband_k,3)
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
                ctermdb = ctermdb + CONJG(cpi)*(j_dpc)*ddir(klmn,itypat,bdir)*cpdg
                ctermdg = ctermdg + CONJG(cpdb)*(-j_dpc)*ddir(klmn,itypat,gdir)*cpj
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

end subroutine vva_fh
!!***

!!****f* ABINIT/vvb_fgh
!! NAME
!! vvb_fgh
!!
!! FUNCTION
!! compute VVb contribution to f, g, and h functions
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

subroutine vvb_fgh(atindx,cprj_k,ddir,dtset,eig_k,ff_k,gg_k,hh_k,lmn2max,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: ddir(lmn2max,dtset%ntypat,3),eig_k(nband_k)
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
              & nn,lmn2max,nband_k,pawtab,cnpb)
            call tdt_me(atindx,cprj_k,np,ddir,dtset,gdir,cprj_k,&
              & nn,lmn2max,nband_k,pawtab,cnpg)

            cnp = c2*cbc*eabg*CONJG(cnpb)*cnpg
            ff_k(1,nn,adir) = ff_k(1,nn,adir) - REAL(cnp)
            ff_k(2,nn,adir) = ff_k(2,nn,adir) - AIMAG(cnp)
            
            cnp = c2*com*eabg*eig_k(nn)*CONJG(cnpb)*cnpg
            hh_k(1,nn,adir) = hh_k(1,nn,adir) - REAL(cnp)
            hh_k(2,nn,adir) = hh_k(2,nn,adir) - AIMAG(cnp)
            
            cnp = c2*com*eabg*eig_k(np)*CONJG(cnpb)*cnpg
            gg_k(1,nn,adir) = gg_k(1,nn,adir) - REAL(cnp)
            gg_k(2,nn,adir) = gg_k(2,nn,adir) - AIMAG(cnp)

          end do ! loop over np

        end do ! gdir
      end do ! bdir
    end do ! adir
  end do !loop over bands

end subroutine vvb_fgh
!!***

!!****f* ABINIT/cc_fgh
!! NAME
!! cc_be
!!
!! FUNCTION
!! compute CC contribution to f, g, and h functions
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

subroutine cc_fgh(atindx,cprj1_k,dimlmn,dtset,eig_k,ff_k,gg_k,gs_hamk,&
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

          ct = cbc*c2*eabg*cme
          ff_k(1,nn,adir) = ff_k(1,nn,adir) + real(ct)
          ff_k(2,nn,adir) = ff_k(2,nn,adir) + aimag(ct)
          
          ct = com*c2*eabg*cme*eig_k(nn)
          hh_k(1,nn,adir) = hh_k(1,nn,adir) + real(ct)
          hh_k(2,nn,adir) = hh_k(2,nn,adir) + aimag(ct)

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

end subroutine cc_fgh
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

!!****f* ABINIT/make_ddir
!! NAME
!! make_ddir
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

subroutine make_ddir(ddir,dtset,gprimd,lmn2max,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max
  type(dataset_type),intent(in) :: dtset

  !arrays
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out) :: ddir(lmn2max,dtset%ntypat,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,itypat,klmn
  real(dp) :: cdij

  !arrays
  integer,dimension(3) :: adir_to_sij = (/4,2,3/)
  real(dp) :: ddir_cart(3)

!--------------------------------------------------------------------

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)
 ddir=zero
 do itypat = 1, dtset%ntypat
   do klmn=1, pawtab(itypat)%lmn2_size
     do adir = 1, 3
       ! note adir 1,2,3 corresponds to cartesian x,y,z. In the abinit real spherical 
       ! harmonics, x~m=1 so LM=4; y~m=-1 so LM=2; z~m=0 so LM=3
       ! adir_to_sij encodes this map
       ddir_cart(adir) = cdij*pawtab(itypat)%qijl(adir_to_sij(adir),klmn)
     end do
     ! now have ddir in cart coords, convert to crystal coords where ddk wavefunctions are
     ddir(klmn,itypat,1:3) = MATMUL(TRANSPOSE(gprimd),ddir_cart)
   end do ! end loop over klmn
 end do ! end loop over types
 
end subroutine make_ddir
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
  ! iterm <= ibcd is berry type, others are orb mag type
  if ( iterm <= ibvvb ) then
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
 do iterms = ibcc,ibvvb
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
   write(message,'(a,3es16.8)')'    CC: H : ',(orbmag_trace(1,adir,imcch),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   VVb: H : ',(orbmag_trace(1,adir,imvvbh),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     VV L : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   VV A.p : ',(orbmag_trace(1,adir,iomanp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'    CC: E : ',(orbmag_trace(1,adir,imcce),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'    CV: E : ',(orbmag_trace(1,adir,imcve),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   VVa: E : ',(orbmag_trace(1,adir,imvvae),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'   VVb: E : ',(orbmag_trace(1,adir,imvvbe),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)')'     Lamb : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  C-C : ',(orbmag_trace(1,adir,ibcc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  C-V : ',(orbmag_trace(1,adir,ibcv),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' V-Va : ',(orbmag_trace(1,adir,ibvva),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' V-Vb : ',(orbmag_trace(1,adir,ibvvb),adir=1,3)
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
     write(message,'(a,3es16.8)') '    CC: H : ',(orbmag_terms(1,iband,adir,imcch),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   VVb: H : ',(orbmag_terms(1,iband,adir,imvvbh),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   VV L : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' VV A.p : ',(orbmag_terms(1,iband,adir,iomanp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    CC: E : ',(orbmag_terms(1,iband,adir,imcce),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    CV: E : ',(orbmag_terms(1,iband,adir,imcve),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   VVa: E : ',(orbmag_terms(1,iband,adir,imvvae),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   VVb: E : ',(orbmag_terms(1,iband,adir,imvvbe),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  C-C : ',(orbmag_terms(1,iband,adir,ibcc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  C-V : ',(orbmag_terms(1,iband,adir,ibcv),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' V-Va : ',(orbmag_terms(1,iband,adir,ibvva),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' V-Vb : ',(orbmag_terms(1,iband,adir,ibvvb),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_ptpaw_output
!!***

end module m_orbmag
