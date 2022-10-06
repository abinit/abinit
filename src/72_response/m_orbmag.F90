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
  use m_pawxc,            only : pawxc
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


  ! map from adir = 1, 2, 3 to S_1,adir packed for spherical harmonic calls
  ! x = r*S_1(1), y = r*S_1(-1), z = r*S_1(0)  
  integer,parameter :: pack1a(3) = (/4,2,3/)

  ! these parameters name the various output terms                                             
  integer,parameter :: ichern1=1,ichern2=2
  integer,parameter :: nterms=2

  ! these parameters name the various d terms
  integer,parameter :: idsij=1
  !integer,parameter :: idp2=1,idpa=2,idvhnzc=3,idvha=4,idvhnhat=5,idsij=6,idvxc=7
  integer,parameter :: ndterms=1

  ! these parameters are constants used repeatedly
  
  ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  real(dp),parameter :: c2=one/(two_pi*two_pi) 
  complex(dpc),parameter :: cbc=j_dpc/two_pi ! Berry curvature pre-factor
  complex(dpc),parameter :: com=-half*j_dpc  ! Orbital magnetism pre-factor

  ! local datatype for d_\alpha terms
  type,private :: dterm_type
    ! scalars
    integer :: lmn2max
    integer :: natom

    ! real arrays

    ! not really a d term but convenient to have it with the others in the same format
    ! <phi|phi> - <tphi|tphi>
    ! qij(natom,lmn2max)
    complex(dpc),allocatable :: qij(:,:)
    
    ! <phi|r_alpha|phi> - <tphi|r_alpha|tphi>
    ! dqij(natom,lmn2max,3)
    complex(dpc),allocatable :: dqij(:,:,:)

  end type dterm_type

  ! local datatype for pawdensities
  type,private :: paw_sph_den_type
    ! integer scalars
    integer :: lm_size
    integer :: mesh_size
    integer :: cplex_density
    integer :: nspden

    ! logical arrays
    logical,allocatable :: lmselectin(:)
    logical,allocatable :: lmselectout(:)

    ! real arrays
    real(dp),allocatable :: rho1(:,:,:)
    real(dp),allocatable :: trho1(:,:,:)
    real(dp),allocatable :: nhat1(:,:,:)

  end type paw_sph_den_type

  ! Bound methods:

  public :: orbmag_ptpaw

  private :: make_chern_v
  private :: make_d
  private :: make_ddir_sij
  private :: make_ddir_vhnzc
  private :: make_ddir_p2
  private :: make_ddir_ap
  private :: make_ddir_vha
  private :: make_ddir_vhnhat
  private :: make_ddir_vxc
  private :: make_fgh_c
  private :: tdt_me
  private :: dtdt_me

  private :: d2lr_p2
  private :: d2lr_Anp

  private :: apply_d2lr_term_k
  private :: lamb_core
  private :: make_pcg1
  private :: pack_pawrhoij
  private :: orbmag_ptpaw_output
  private :: paw_sph_den_alloc
  private :: paw_sph_den_free
  private :: dterm_alloc
  private :: dterm_free
  
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
 type(dterm_type) :: dterm
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),ct_k(:,:,:),cwavef(:,:)
 real(dp),allocatable :: eig_k(:),ffnl_k(:,:,:,:),gg_k(:,:,:),hh_k(:,:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:),orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
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

 ! note: in make_d, terms will be filled as iatom using atindx
 call dterm_alloc(dterm,lmn2max,dtset%natom)

 call make_d(atindx,atindx1,cprj,dimlmn,dterm,dtset,gprimd,mcprj,&
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

 ABI_MALLOC(ct_k,(2,nband_k,3))
 ABI_MALLOC(gg_k,(2,nband_k,3))
 ABI_MALLOC(hh_k,(2,nband_k,3))

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

   call make_fgh_c(atindx,cprj1_k,dimlmn,dtset,eig_k,ct_k,gg_k,gs_hamk,hh_k,ikpt,isppol,&
     & mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pcg1_k)
   orbmag_terms(:,:,:,ichern1) = orbmag_terms(:,:,:,ichern1) + trnrm*ct_k(:,:,:)

   call make_chern_v(atindx,cprj_k,cprj1_k,ct_k,dterm,dtset,nband_k,pawtab)
   orbmag_terms(:,:,:,ichern2) = orbmag_terms(:,:,:,ichern2) + trnrm*ct_k(:,:,:)
   
   !!--------------------------------------------------------------------------------
   !! onsite <phi|r_b p^2/2 r_g>
   !!--------------------------------------------------------------------------------
   !call apply_d2lr_term_k(atindx,cprj_k,dtset,iomlr,lmn2max,mp2,nband_k,gg_k,pawtab)
   !orbmag_terms(:,:,:,iomlr) = orbmag_terms(:,:,:,iomlr) + trnrm*gg_k
   !
   !!--------------------------------------------------------------------------------
   !! onsite <phi|r_b p.A0 r_g> 
   !!--------------------------------------------------------------------------------
   !call apply_d2lr_term_k(atindx,cprj_k,dtset,iomanp,lmn2max,mpan,nband_k,gg_k,pawtab)
   !orbmag_terms(:,:,:,iomanp) = orbmag_terms(:,:,:,iomanp) + trnrm*gg_k

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
   do nn = 1, nband_k
     do icmplx = 1, 2
       orbmag_terms(icmplx,nn,1:3,iterm) = ucvol*MATMUL(gprimd,orbmag_terms(icmplx,nn,1:3,iterm))
     end do
   end do
 end do

 !! convert orbmag magnetization to orbital moment
 !! Berry curvature terms are ignored
 !! Lamb term ignored
 !do iterm = iomlr, nterms
 !  if (iterm .EQ. iomlmb) cycle
 !  orbmag_terms(:,:,:,iterm) = ucvol*orbmag_terms(:,:,:,iterm)
 !end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace = orbmag_trace + orbmag_terms(:,nn,:,:)
 end do

!! get the Lamb term
! call lamb_core(atindx,dtset,omlamb)
! orbmag_trace(:,:,iomlmb) = omlamb
 
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

 call dterm_free(dterm)

 ABI_FREE(mp2)
 ABI_FREE(mpan)
 
 ABI_FREE(orbmag_terms)
 ABI_FREE(ct_k)
 ABI_FREE(gg_k)
 ABI_FREE(hh_k)
 ABI_FREE(orbmag_trace)

end subroutine orbmag_ptpaw
!!***

!!****f* ABINIT/make_chern_v
!! NAME
!! make_chern_v
!!
!! FUNCTION
!! compute berry curvature originating from all terms but Pc/Pc
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

subroutine make_chern_v(atindx,cprj_k,cprj1_k,ct,dterm,dtset,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: ct(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn,np
  real(dp) :: dotr,doti,epsabg
  complex(dpc) :: ci,cj,cpi,cpj,cterm,dtdt,tdt,tdtp
  !arrays

!--------------------------------------------------------------------

 ct = zero

 do adir = 1, 3
   do bdir = 1, 3
     do gdir = 1, 3
       epsabg = eijk(adir,bdir,gdir)
       if (ABS(epsabg) .LT. tol8) cycle
       do nn = 1, nband_k
         
         cterm = czero

         call tdt_me(dterm%qij,atindx,cprj1_k(:,nn,bdir),dterm%dqij,dtset,gdir,cprj_k(:,nn),&
           & dterm%lmn2max,pawtab,tdt)
         cterm = cterm + tdt
         
         call tdt_me(dterm%qij,atindx,cprj1_k(:,nn,gdir),dterm%dqij,dtset,bdir,cprj_k(:,nn),&
           & dterm%lmn2max,pawtab,tdt)
         cterm = cterm + CONJG(tdt)
         
         call dtdt_me(dterm%qij,atindx,cprj_k(:,nn),bdir,dterm%dqij,dtdt,dtset,&
           & gdir,cprj_k(:,nn),dterm%lmn2max,pawtab)
         cterm = cterm + dtdt
         
         cterm=cbc*c2*epsabg*cterm
         ct(1,nn,adir) = ct(1,nn,adir) +  REAL(cterm)
         ct(2,nn,adir) = ct(2,nn,adir) + AIMAG(cterm)

         cterm = czero
         do np = 1, nband_k
           call tdt_me(dterm%qij,atindx,cprj_k(:,np),dterm%dqij,dtset,gdir,cprj_k(:,nn),&
             & dterm%lmn2max,pawtab,tdt)
           
           call tdt_me(dterm%qij,atindx,cprj_k(:,np),dterm%dqij,dtset,bdir,cprj_k(:,nn),&
             & dterm%lmn2max,pawtab,tdtp)
           cterm = cterm + CONJG(tdtp)*tdt
         end do !np
         cterm=cbc*c2*epsabg*cterm
         ct(1,nn,adir) = ct(1,nn,adir) -  REAL(cterm)
         ct(2,nn,adir) = ct(2,nn,adir) - AIMAG(cterm)

       end do !nn
     end do !gdir
   end do !bdir
 end do !adir
 
end subroutine make_chern_v
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

!!****f* ABINIT/dtdt_me
!! NAME
!! dtdt_me
!!
!! FUNCTION
!! Matrix element <u_n|(\partial_bdir T)^\dag A (\partial_gdir T)|u_m>, for one atom
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! aij :: scalar ij couplings
!! bcp :: bra side cprj
!! bdir :: bra side deriv direc
!! daij :: vector ij couplings
!! kcp :: ket side cprj
!! kdir :: ket side deriv direc
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! not for use with momentum terms, we assume here the double deriv coupling is zero
!! as is the case for Hartree like terms
!!
!! SOURCE

subroutine dtdt_me(aij,atindx,bcp,bdir,daij,dtdt,dtset,gdir,kcp,lmn2max,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,gdir,lmn2max
  complex(dpc),intent(out) :: dtdt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max),daij(dtset%natom,lmn2max,3)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom),kcp(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,jlmn,klmn
  complex(dpc) :: cpi,cpj,dcpi,dcpj

  !arrays

!--------------------------------------------------------------------

  dtdt = czero
  do iat = 1, dtset%natom
    itypat=dtset%typat(iat)
    iatom = atindx(iat)
    do ilmn = 1, pawtab(itypat)%lmn_size
      do jlmn = 1, pawtab(itypat)%lmn_size
        klmn = MATPACK(ilmn,jlmn)
        cpi =  CMPLX(bcp(iatom)%cp(1,ilmn),bcp(iatom)%cp(2,ilmn))
        cpj =  CMPLX(kcp(iatom)%cp(1,jlmn),kcp(iatom)%cp(2,jlmn))
        dcpi = CMPLX(bcp(iatom)%dcp(1,bdir,ilmn),bcp(iatom)%dcp(2,bdir,ilmn))
        dcpj = CMPLX(kcp(iatom)%dcp(1,gdir,jlmn),kcp(iatom)%dcp(2,gdir,jlmn))
        dtdt = dtdt + CONJG(dcpi)*dcpj*aij(iatom,klmn)
        dtdt = dtdt + CONJG(dcpi)*cpj*daij(iatom,klmn,gdir)
        dtdt = dtdt + CONJG(cpi)*dcpj*CONJG(daij(iatom,klmn,bdir))
      end do !jlmn
    end do !ilmn
  end do !iat

end subroutine dtdt_me
!!***


!!****f* ABINIT/tdt_me
!! NAME
!! tdt_me
!!
!! FUNCTION
!! Matrix element <u_n|T^\dag A (\partial_idir T)|u_m>
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! aij :: scalar ij couplings
!! bcp :: bra side cprj
!! daij :: vector ij couplings
!! idir :: direction of interest for derivatives
!! kcp :: ket side cprj
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

subroutine tdt_me(aij,atindx,bcp,daij,dtset,idir,kcp,lmn2max,pawtab,tdt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,lmn2max
  complex(dpc),intent(out) :: tdt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max),daij(dtset%natom,lmn2max,3)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom),kcp(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,jlmn,klmn
  complex(dpc) :: cpi,cpj,dcpj

  !arrays

!--------------------------------------------------------------------

  tdt = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat=dtset%typat(iat)
    do ilmn = 1, pawtab(itypat)%lmn_size
      do jlmn = 1, pawtab(itypat)%lmn_size
        klmn = MATPACK(ilmn,jlmn)
        cpi =  CMPLX(bcp(iatom)%cp(1,ilmn),bcp(iatom)%cp(2,ilmn))
        cpj =  CMPLX(kcp(iatom)%cp(1,jlmn),kcp(iatom)%cp(2,jlmn))
        dcpj = CMPLX(kcp(iatom)%dcp(1,idir,jlmn),kcp(iatom)%dcp(2,idir,jlmn))
        tdt = tdt + CONJG(cpi)*dcpj*aij(iatom,klmn) + &
          & CONJG(cpi)*cpj*daij(iatom,klmn,idir)
      end do !jlmn
    end do !ilmn
  end do !iat

end subroutine tdt_me
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

subroutine make_ddir_sij(atindx,dterm,dtset,gprimd,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,itypat,ilmn,jlmn,klmn
  real(dp) :: cdij

  !arrays
  real(dp) :: ddir_cart(3),ddir_red(3)

!--------------------------------------------------------------------

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)
 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   do klmn = 1, pawtab(itypat)%lmn2_size
     
     ! qij is not really a d term but it's convenient to have it 
     ! in the same format as the d terms 
     dterm%qij(iatom,klmn)=CMPLX(pawtab(itypat)%sij(klmn),zero)

     do adir = 1, 3
       ddir_cart(adir) = cdij*pawtab(itypat)%qijl(pack1a(adir),klmn)
     end do !adir
    
     ! now have ddir in cart coords, convert to crystal coords where ddk wavefunctions are
     ddir_red = MATMUL(TRANSPOSE(gprimd),ddir_cart)
     dterm%dqij(iatom,klmn,1:3) = -j_dpc*ddir_red(1:3)
   end do ! klmn
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
  if ( iterm <= ichern2 ) then
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

            lb = pack1a(bdir)
            lm1b = MATPACK(ilm,lb)

            lg = pack1a(gdir)
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
 !do iterms = iomlr, nterms
 !  orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
 !  do iband=1, nband_k
 !    orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
 !  end do
 !end do
 berry_bb=zero;berry_total=zero
 do iterms = ichern1,ichern2
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
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Pc terms : ',(orbmag_trace(1,adir,ichern1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Pv terms : ',(orbmag_trace(1,adir,ichern2),adir=1,3)
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
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Pc terms : ',(orbmag_terms(1,iband,adir,ichern1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Pv terms : ',(orbmag_terms(1,iband,adir,ichern2),adir=1,3)
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
         lm1b = pack1a(adir)
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

subroutine make_ddir_vhnhat(atindx,ddir_vhnhat,dtset,gntselect,gprimd,lmnmax,my_lmax,&
    & pawrad,pawsphden,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  complex(dpc),intent(out) :: ddir_vhnhat(lmnmax,lmnmax,dtset%natom,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(paw_sph_den_type),intent(in) :: pawsphden(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,g2int,iat,iatom,ilm,ilmn,ilmp,imesh,itypat
  integer :: jlmn,jmesh,klmadir,klmn,klm,kln,ll,lmax,lmin,lp
  integer :: mesh_size
  real(dp) :: cdij,nlt1,rfac,rr,rp,vhaint

  !arrays
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),fft1(:)

!--------------------------------------------------------------------

 cdij = four_pi*sqrt(four_pi/3.0D0)
 ddir_vhnhat = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(fft1,(mesh_size))
  
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
               klmadir = MATPACK(ilm,pack1a(adir))
               g2int = gntselect(ilmp,klmadir)
               if (g2int .EQ. 0) cycle

               ! construct integrand for rho1 vHa, trho1
               do imesh = 2, mesh_size
                 rr = pawrad(itypat)%rad(imesh)
                   
                 ! for this mesh point, do the interior nonlocal integral over Hartree potential
                 do jmesh = 2, imesh
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rp**lp)/(rr**(lp+1))
                   fft1(jmesh) = pawsphden(iatom)%nhat1(jmesh,ilmp,1)*rfac
                 end do
                 do jmesh=imesh+1, mesh_size
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rr**lp)/(rp**(lp+1))
                   fft1(jmesh) = pawsphden(iatom)%nhat1(jmesh,ilmp,1)*rfac
                 end do

                 call pawrad_deducer0(fft1,mesh_size,pawrad(itypat))
                 call simp_gen(nlt1,fft1,pawrad(itypat))
                 !write(std_out,'(a,i4,3es16.8)')' JWZ debug imesh rr nlt1 tphitphj ',imesh,rr,nlt1,&
                 !  & pawtab(itypat)%tphitphj(imesh,kln)
                 
                 ff(imesh)=-rr*pawtab(itypat)%tphitphj(imesh,kln)*nlt1
               
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

   ABI_FREE(ff)
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

subroutine make_ddir_vxc(atindx,ddir_vxc,dtset,gntselect,gprimd,lmnmax,my_lmax,&
    & pawang,pawrad,pawsphden,pawtab,realgnt)

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
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(paw_sph_den_type),intent(in) :: pawsphden(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,ilm,ilmn,imesh,ipt,itypat,jlm,jlmn
  integer :: klmn,kln,mesh_size,nkxc,nk3xc,usecore,usexcnhat,xc_option
  real(dp) :: eexc,eexcdc,hyb_mixing,rr,xcint
  logical :: non_magnetic_xc

  !arrays
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),kxc(:,:,:),k3xc(:,:,:)
  real(dp),allocatable :: vxc(:,:,:),tvxc(:,:,:)

!--------------------------------------------------------------------

 hyb_mixing = zero
 nkxc = 0
 nk3xc = 0
 non_magnetic_xc = .FALSE.
 xc_option = 0
 usecore = 1
 usexcnhat = 0

 ddir_vxc = zero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   mesh_size = pawsphden(iatom)%mesh_size

   ABI_MALLOC(vxc,(mesh_size,pawang%angl_size,dtset%nspden))
   ABI_MALLOC(tvxc,(mesh_size,pawang%angl_size,dtset%nspden))
   ABI_MALLOC(ff,(mesh_size))
         
   call pawxc(pawtab(itypat)%coredens,eexc,eexcdc,hyb_mixing,dtset%ixc,kxc,k3xc,&
     & pawsphden(iatom)%lm_size,pawsphden(iatom)%lmselectout,pawsphden(iatom)%nhat1,&
     & nkxc,nk3xc,non_magnetic_xc,mesh_size,pawsphden(iatom)%nspden,&
     & xc_option,pawang,pawrad(itypat),pawsphden(iatom)%rho1,usecore,usexcnhat,&
     & vxc,dtset%xclevel,dtset%xc_denpos)
         
   call pawxc(pawtab(itypat)%tcoredens(:,1),eexc,eexcdc,hyb_mixing,dtset%ixc,kxc,k3xc,&
     & pawsphden(iatom)%lm_size,pawsphden(iatom)%lmselectout,pawsphden(iatom)%nhat1,&
     & nkxc,nk3xc,non_magnetic_xc,mesh_size,pawsphden(iatom)%nspden,&
     & xc_option,pawang,pawrad(itypat),pawsphden(iatom)%trho1,usecore,usexcnhat,&
     & tvxc,dtset%xclevel,dtset%xc_denpos)

   do klmn = 1, pawtab(itypat)%lmn2_size
     ilmn = pawtab(itypat)%indklmn(7,klmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn)
     
     ilm  = pawtab(itypat)%indklmn(5,klmn)
     jlm  = pawtab(itypat)%indklmn(6,klmn)
     
     kln  = pawtab(itypat)%indklmn(2,klmn)

     dij_cart = zero
     do ipt = 1, pawang%angl_size
       do imesh = 2, mesh_size
         rr = pawrad(itypat)%rad(imesh)
         ff(imesh) = rr*vxc(imesh,ipt,1)*pawtab(itypat)%phiphj(imesh,kln) - &
           & rr*tvxc(imesh,ipt,1)*pawtab(itypat)%tphitphj(imesh,kln)
       end do !imesh
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(xcint,ff,pawrad(itypat))

       do adir = 1, 3
         dij_cart(adir) = dij_cart(adir) + &
           & pawang%angwgth(ipt)*pawang%ylmr(ilm,ipt)*pawang%ylmr(jlm,ipt)*&
           & pawang%ylmr(pack1a(adir),ipt)*xcint
       end do ! adir
     end do ! ipt

     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

     ddir_vxc(ilmn,jlmn,iat,1:3) = dij_red(1:3)
     ddir_vxc(jlmn,ilmn,iat,1:3) = dij_red(1:3)

   end do ! klmn

   ABI_FREE(vxc)
   ABI_FREE(tvxc)
   ABI_FREE(ff)
 
 end do !iat

end subroutine make_ddir_vxc


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
!!
!! SOURCE

subroutine make_ddir_vha(atindx,ddir_vha,dtset,gntselect,gprimd,lmnmax,my_lmax,&
    & pawrad,pawsphden,pawtab,realgnt)

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
  type(paw_sph_den_type),intent(in) :: pawsphden(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,g2int,iat,iatom,ilm,ilmn,ilmp,imesh,itypat
  integer :: jlmn,jmesh,klmadir,klmn,klm,kln,ll,lmin,lmax,lp
  integer :: mesh_size
  real(dp) :: cdij,nl1,nlt1,rfac,rr,rp,vhaint,rrhokl
  real(dp) :: nterm

  !arrays
  real(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:),ff1(:),fft1(:)

!--------------------------------------------------------------------

 cdij = four_pi*sqrt(four_pi/3.0D0)
 ddir_vha = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(ff1,(mesh_size))
   ABI_MALLOC(fft1,(mesh_size))
  
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
               klmadir = MATPACK(ilm,pack1a(adir))
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
                   ff1(jmesh) = pawsphden(iatom)%rho1(jmesh,ilmp,1)*rfac
                   fft1(jmesh) = pawsphden(iatom)%trho1(jmesh,ilmp,1)*rfac
                 end do
                 do jmesh=imesh+1, mesh_size
                   rp = pawrad(itypat)%rad(jmesh)
                   rfac = (rp**2)*(rr**lp)/(rp**(lp+1))
                   ff1(jmesh) = pawsphden(iatom)%rho1(jmesh,ilmp,1)*rfac
                   fft1(jmesh) = pawsphden(iatom)%trho1(jmesh,ilmp,1)*rfac
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

   ABI_FREE(ff)
   ABI_FREE(ff1)
   ABI_FREE(fft1)

 end do ! end loop over iat
 
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
          lm1b = pack1a(adir)
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
          l1b = pack1a(adir)
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

!!****f* ABINIT/dterm_free
!! NAME
!! dterm_free
!!
!! FUNCTION
!! free dterm_type
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

subroutine dterm_free(dterm)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm

  !arrays

  !Local variables -------------------------
  !scalars
 
  !arrays
!--------------------------------------------------------------------

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if

  if(allocated(dterm%dqij)) then
    ABI_FREE(dterm%dqij)
  end if
 
end subroutine dterm_free
!!***

!!****f* ABINIT/paw_sph_den_free
!! NAME
!! paw_sph_den_free
!!
!! FUNCTION
!! free paw_sph_den_type
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

subroutine paw_sph_den_free(paw_sph_den)

  !Arguments ------------------------------------
  !scalars
  type(paw_sph_den_type),intent(inout) :: paw_sph_den

  !arrays

  !Local variables -------------------------
  !scalars
 
  !arrays
!--------------------------------------------------------------------

  if(allocated(paw_sph_den%lmselectin)) then
    ABI_FREE(paw_sph_den%lmselectin)
  end if
  
  if(allocated(paw_sph_den%lmselectout)) then
    ABI_FREE(paw_sph_den%lmselectout)
  end if
  
  if(allocated(paw_sph_den%rho1)) then
    ABI_FREE(paw_sph_den%rho1)
  end if

  if(allocated(paw_sph_den%trho1)) then
    ABI_FREE(paw_sph_den%trho1)
  end if

  if(allocated(paw_sph_den%nhat1)) then
    ABI_FREE(paw_sph_den%nhat1)
  end if
 
end subroutine paw_sph_den_free
!!***

!!****f* ABINIT/dterm_alloc
!! NAME
!! dterm_alloc
!!
!! FUNCTION
!! allocate dterm_type
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

subroutine dterm_alloc(dterm,lmn2max,natom)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,natom
  type(dterm_type),intent(inout) :: dterm

  !arrays

  !Local variables -------------------------
  !scalars
 
  !arrays
!--------------------------------------------------------------------

  dterm%lmn2max = lmn2max
  dterm%natom = natom

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if
  ABI_MALLOC(dterm%qij,(natom,lmn2max))

  if(allocated(dterm%dqij)) then
    ABI_FREE(dterm%dqij)
  end if
  ABI_MALLOC(dterm%dqij,(natom,lmn2max,3))
  
end subroutine dterm_alloc
!!***

!!****f* ABINIT/paw_sph_den_alloc
!! NAME
!! paw_sph_den_alloc
!!
!! FUNCTION
!! allocate paw_sph_den_type
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

subroutine paw_sph_den_alloc(cplex,lm_size,mesh_size,nspden,paw_sph_den)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: cplex,lm_size,mesh_size,nspden
  type(paw_sph_den_type),intent(inout) :: paw_sph_den

  !arrays

  !Local variables -------------------------
  !scalars
 
  !arrays
!--------------------------------------------------------------------

  paw_sph_den%lm_size = lm_size
  paw_sph_den%mesh_size = mesh_size
  paw_sph_den%cplex_density = cplex
  paw_sph_den%nspden = nspden

  if(allocated(paw_sph_den%lmselectin)) then
    ABI_FREE(paw_sph_den%lmselectin)
  end if
  ABI_MALLOC(paw_sph_den%lmselectin,(lm_size))
  
  if(allocated(paw_sph_den%lmselectout)) then
    ABI_FREE(paw_sph_den%lmselectout)
  end if
  ABI_MALLOC(paw_sph_den%lmselectout,(lm_size))
  
  if(allocated(paw_sph_den%rho1)) then
    ABI_FREE(paw_sph_den%rho1)
  end if
  ABI_MALLOC(paw_sph_den%rho1,(cplex*mesh_size,lm_size,nspden))

  if(allocated(paw_sph_den%trho1)) then
    ABI_FREE(paw_sph_den%trho1)
  end if
  ABI_MALLOC(paw_sph_den%trho1,(cplex*mesh_size,lm_size,nspden))

  if(allocated(paw_sph_den%nhat1)) then
    ABI_FREE(paw_sph_den%nhat1)
  end if
  ABI_MALLOC(paw_sph_den%nhat1,(cplex*mesh_size,lm_size,nspden))
 
end subroutine paw_sph_den_alloc
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
      if ( (ABS(rr) .GT. tol16) .OR. (ABS(ri) .GT. tol16) ) then
        isel = isel + 1
      end if
    end do ! end loop over klmn
    pawrhoij(iatom)%nrhoijsel = isel

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
      if ( (ABS(rr) .GT. tol16) .OR. (ABS(ri) .GT. tol16) ) then
        isel = isel + 1
        pawrhoij(iatom)%rhoijselect(isel) = klmn
        if (cprho) then
          pawrhoij(iatom)%rhoijp(2*isel-1,1) = rr
          pawrhoij(iatom)%rhoijp(2*isel,1) = ri
        else
          pawrhoij(iatom)%rhoijp(isel,1) = rr
        end if ! if cprho
      end if ! if rhoij nonzero
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

subroutine make_d(atindx,atindx1,cprj,dimlmn,dterm,dtset,gprimd,&
    & mcprj,mpi_enreg,occ,paw_an,pawang,pawrad,pawtab,psps)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcprj
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
  integer,intent(in) :: dimlmn(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,iband,isel,itypat,ilmn,jlmn,klmn
  integer :: my_lmax,ngnt,nzlmopt
  integer :: opt_compch,opt_dens,opt_l,opt_print
  real(dp) :: compch_sph
  complex(dpc) :: cpi,cpj,rhoij
  type(paw_dmft_type) :: paw_dmft
 
  !arrays
  integer,allocatable :: gntselect(:,:)
  real(dp),allocatable :: realgnt(:)
  type(pawrhoij_type),allocatable :: pawrhoij(:)
  type(paw_sph_den_type),allocatable :: pawsphden(:)
!--------------------------------------------------------------------

 !make pawrhoij
 ABI_MALLOC(pawrhoij,(dtset%natom))
 call pawrhoij_alloc(pawrhoij,dtset%pawcpxocc,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%typat,&
   & use_rhoijp=1,use_rhoij_=1,pawtab=pawtab)
 paw_dmft%use_sc_dmft=0
 paw_dmft%use_dmft=0
 call pawmkrhoij(atindx,atindx1,cprj,dimlmn,dtset%istwfk,dtset%kptopt,dtset%mband,dtset%mband,&
   & mcprj,dtset%mkmem,mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,dtset%nsppol,&
   & occ,dtset%paral_kgb,paw_dmft,pawrhoij,0,dtset%usewvl,dtset%wtk)
 call pack_pawrhoij(dtset,pawrhoij)

 ! make Gaunt integrals
 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 ! make onsite densities
 nzlmopt = -1
 opt_compch = 0
 opt_dens = 0
 opt_l = -1
 opt_print = 0
 ABI_MALLOC(pawsphden,(dtset%natom))
 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   call paw_sph_den_alloc(pawrhoij(iatom)%qphase,paw_an(iatom)%lm_size,&
     & pawtab(itypat)%mesh_size,dtset%nspden,pawsphden(iatom))

   pawsphden(iatom)%lmselectin = .TRUE. 
   call pawdensities(compch_sph,pawsphden(iatom)%cplex_density,iatom,&
     & pawsphden(iatom)%lmselectin,pawsphden(iatom)%lmselectout,&
     & pawsphden(iatom)%lm_size,pawsphden(iatom)%nhat1,dtset%nspden,&
     & nzlmopt,opt_compch,opt_dens,opt_l,opt_print,pawang,dtset%pawprtvol,&
     & pawrad(itypat),pawrhoij(iatom),pawtab(itypat),&
     & pawsphden(iatom)%rho1,pawsphden(iatom)%trho1)
 
 end do ! iat

 ! generate d terms

 call make_ddir_sij(atindx,dterm,dtset,gprimd,pawtab)

 !! term idp2 due to onsite p^2/2, corresponds to term 1 of Torrent PAW roadmap paper appendix E
 !! Comp. Mat. Sci. 42, 337-351 (2008)
 !call make_ddir_p2(dterms(idp2,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)

 !! term idpa due to onsite A.p, not part of the roadmap paper but similar to term 1
 !call make_ddir_ap(dterms(idpa,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)
 !
 !! term idvha due to v_H[n1], corresponds to term 2a of roadmap paper
 !call make_ddir_vha(atindx,dterms(idvha,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,&
 !  & pawrad,pawsphden,pawtab,realgnt)

 !! term idvhnzc due to v_H[nZ], corresponds to term 2b of roadmap paper
 !call make_ddir_vhnzc(dterms(idvhnzc,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,pawrad,pawtab,realgnt)

 !! term idvhnhat due to v_H[nhat], corresponds to term 2c of roadmap paper
 !call make_ddir_vhnhat(atindx,dterms(idvhnhat,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,&
 !  & pawrad,pawsphden,pawtab,realgnt)

 !! term idvxc due to v_xc[n1+nc] corresponds to term 3a in roadmap paper
 !call make_ddir_vxc(atindx,dterms(idvxc,:,:,:,:),dtset,gntselect,gprimd,psps%lmnmax,my_lmax,&
 !   & pawang,pawrad,pawsphden,pawtab,realgnt)

 ABI_FREE(realgnt)
 ABI_FREE(gntselect)

 do iat = 1, dtset%natom
   call paw_sph_den_free(pawsphden(iat))
 end do
 ABI_FREE(pawsphden)
 
 call pawrhoij_free(pawrhoij)
 ABI_FREE(pawrhoij)
 
end subroutine make_d
!!***

!!****f* ABINIT/make_fgh_c
!! NAME
!! make_fgh_c
!!
!! FUNCTION
!! compute all terms chern, H, E for Pc/Pc
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

subroutine make_fgh_c(atindx,cprj1_k,dimlmn,dtset,eig_k,ff_k,gg_k,gs_hamk,&
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

end subroutine make_fgh_c
!!***



end module m_orbmag
