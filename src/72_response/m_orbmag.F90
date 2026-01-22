!!*** ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2026 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! These routines implement the theory developed in Zwanziger, Torrent, Gonze
!! Phys Rev B 107, 165157 (2023). This paper will be referred to the comments as ZTG23.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! note: in a typical index over lmn2_size, think of it as row ilmn, column jlmn for element d_ij.
! In pawinit line 356, klmn is constructed such that ilmn <= jlmn. Thus we have the upper triangular
! part of the dij matrix. When looping over both ilmn and jlmn, element dij with i>j is constructed
! by symmetry from element dji.
#define MATPACK(row,col) (MAX(row,col)*(MAX(row,col)-1)/2 + MIN(row,col))

module m_orbmag

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_crystal,          only : crystal_t
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : projbd
  use m_dtfil
  use m_ebands
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : gs_hamiltonian_type, gspot_transgrid_and_pack
  use m_hdr
  use m_kg,               only : getph,mkkin,mkkpg,ph1d3d
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle,proc_distrb_nband
  use m_nctk
  use netcdf
  use m_nonlop,           only : nonlop
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get, pawcprj_put
  use m_pawdij,           only : pawv1
  use m_pawfgr,           only : pawfgr_type
  use m_pawfgrtab,        only : pawfgrtab_type
  use m_paw_ij,           only : paw_ij_type
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen,poisson
  use m_paw_sphharm,      only : setsym_ylm,slxyzs,realgaunt,make_dyadic
  use m_pawtab,           only : pawtab_type
  use m_spacepar,         only : make_vectornd
  use m_time,             only : cwtime, timab

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


  ! these parameters name the various output terms

  ! these parameters name the various output terms
  integer,parameter :: chern_nterms=3
  integer,parameter :: ibcc=1,ibvv1=2,ibvv2=3
  integer,parameter :: orbmag_nterms=6
  integer,parameter :: incc=1,invv1=2,invv2=3
  integer,parameter :: innl=4,inlr=5,inbm=6

  ! these parameters are constants used repeatedly

  ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  real(dp),parameter :: c2=one/(two_pi*two_pi)
  complex(dp),parameter :: com=-half*j_dpc  ! Orbital magnetism pre-factor
  complex(dp),parameter :: cbc=-com ! Berry curvature pre-factor

  ! local datatype for orbmag data on kpt mesh, for eventual output to netcdf
  type,private :: orbmag_mesh_type
    ! scalars

    integer :: mband, nkpt, nsppol
    ! number of bands, kpts, spin polarizations

    integer :: natom, ntypat
    ! atoms and types of atoms

    integer :: chern_nterms
    ! number of chern terms to store on the kpt mesh
    ! CC, VV1, VV2

    integer :: orbmag_nterms
    ! number of orbmag terms to store on the kpt mesh
    ! CC, VV1, VV2, NL, L_R, B.M

    real(dp),allocatable :: lambsig(:)
    ! lambsig(ntypat)

    real(dp),allocatable :: nucdipmom(:,:)
    ! nucdipmom(3,natom)

    real(dp),allocatable :: cmesh(:,:,:,:,:)
    ! 3 for the 3 directions
    ! cmesh(mband,nkpt,nsppol,3,chern_terms)

    real(dp),allocatable :: omesh(:,:,:,:,:)
    ! 3 for the 3 directions
    ! omesh(mband,nkpt,nsppol,3,orbmag_terms)

  end type orbmag_mesh_type

  ! local datatype for various onsite terms. Probably overkill, but convenient.
  type,private :: dterm_type
    ! scalars
    integer :: lmnmax
    integer :: lmn2max
    integer :: natom
    integer :: ndij
    integer :: has_aij=0
    integer :: has_qij=0
    integer :: has_LR=0
    integer :: has_BM=0

    ! sum of \Delta A_ij
    ! typically will be just paw_ij
    ! aij(natom,lmn2max,ndij)
    complex(dp),allocatable :: aij(:,:,:)

    ! <phi|phi> - <tphi|tphi>
    ! qij(natom,lmn2max,ndij)
    complex(dp),allocatable :: qij(:,:,:)

    ! onsite L_R/2
    ! <phi|L_R/2|phi> - <tphi|L_R/2|tphi>
    ! LR(natom,lmn2max,ndij,3)
    complex(dp),allocatable :: LR(:,:,:,:)

    ! onsite BM
    ! <phi|Bxr . mxr|phi> - <tphi|Bxr . mxr|tphi>
    ! BM(natom,lmn2max,ndij,3)
    complex(dp),allocatable :: BM(:,:,:,:)

  end type dterm_type

  ! Bound methods:

  public :: orbmag

  private :: orbmag_cc_k
  private :: orbmag_vv_k
  private :: orbmag_nl_k
  private :: orbmag_nl1_k
  private :: make_d
  private :: dterm_aij
  private :: dterm_qij
  private :: dterm_LR
  private :: dterm_BM
  private :: tt_me
  private :: txt_me
  private :: local_fermie

  private :: lamb_core
  private :: make_pcg1
  private :: orbmag_mesh_alloc
  private :: orbmag_mesh_free
  private :: sum_orbmag_mesh
  private :: orbmag_output
  private :: orbmag_ncwrite   ! Write orbmag_mesh contributions to netcdf file.
  private :: dterm_alloc
  private :: dterm_free

CONTAINS  !========================================================================================
!!***

!!****f* ABINIT/orbmag
!! NAME
!! orbmag
!!
!! FUNCTION
!! This routine computes the orbital magnetization and Berry curvature based on input
!! wavefunctions and DDK wavefuntions.
!!
!! INPUTS
!!  cg(2,mcg)=all ground state wavefunctions
!!  cg1(2,mcg1,3)=all DDK wavefunctions in all 3 directions
!!  cprj(dtset%natom,mcprj)<type(pawcprj_type)>=all ground state cprj
!!  crystal(crystal_t)=structured datatype holding details about unit cell
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ebands_k(ebands_t)=structured datatype holding GS eigenvalues
!!  gsqcut=large sphere cut-off
!!  hdr(hdr_type)=structured dataype with header info for eventual output
!!  kg(3,mpw*mkmem_rbz)=basis sphere of planewaves at k
!!  mcg=dimension of cg
!!  mcg1=dimension of cg1
!!  mcprj=dimension of cprj
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  mpw=max number of planewaves at k
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(18)=FFT grid size information (from pawfgr%ngfft)
!!  paw_ij(dtset%natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(dtset%ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=real space translation vectors
!!  usevxctau=1 if kinetic energy density contribution has to be included (mGGA)
!!  vtrial(nfftf,dtset%nspden)=GS potential (Hartree)
!!  vxctau(nfftf,nspden,4*usevxctau)=derivative of e_xc with respect to kinetic energy density, for mGGA
!!  ylm(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm)=all ylm's
!!  ylmgr(mpw*mkmem_rbz,3,psps%mpsang*psps%mpsang*psps%useylm)=gradients of ylm's
!!
!! OUTPUT
!!  only printing in call to orbmag_output
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! See Zwanziger, Torrent, and Gonze Phys Rev B 107, 165157 (2023), "ZTG23"
!! DDK wavefunctions are used for the derivatives.
!!
!! SOURCE

subroutine orbmag(cg,cg1,cprj,crystal,dtfil,dtset,ebands_k,gsqcut,hdr,kg,mcg,mcg1,&
    & mcprj,mkmem_rbz,mpi_enreg,mpw,nfftf,ngfftf,paw_ij,pawfgr,pawrad,&
    & pawtab,psps,usevxctau,vtrial,vxctau,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcprj,mcg,mcg1,mkmem_rbz,mpw,nfftf,usevxctau
 real(dp),intent(in) :: gsqcut
 type(crystal_t),intent(in) :: crystal
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands_k
 type(hdr_type),intent(in) :: hdr
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps

 !arrays
 integer,intent(in) :: kg(3,mpw*mkmem_rbz),ngfftf(18)
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(inout) :: vxctau(nfftf,dtset%nspden,4*usevxctau)
 real(dp),intent(in) :: ylm(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem_rbz,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 character(len=fnlen) :: fname
 integer :: adir,bdtot_index,buff_size,choice,cpopt,dimffnl,exchn2n3d
 integer :: iat,iatom,icg,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,itypat,lmn2max
 integer :: me,mcgk,mcprjk,my_lmax,my_nspinor,nband_k,nband_me,ncid,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nl1_option,nn,nkpg,npw_k,npwsp,nproc,spaceComm
 integer,parameter :: master=0
 real(dp) :: arg,ecut_eff,fermie
 logical :: has_nucdip
 type(dterm_type) :: dterm
 type(gs_hamiltonian_type) :: gs_hamk
 type(orbmag_mesh_type) :: orbmag_mesh

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: kpoint(3),omlamb(3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: chern_terms(:,:,:,:),chern_trace(:,:),cg_k(:,:),cg1_k(:,:,:),cwavef(:,:)
 real(dp),allocatable :: eig_k(:),ffnl_k(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: occ_k(:),orbmag_terms(:,:,:,:),orbmag_trace(:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: vxctaulocal(:,:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj1_k(:,:,:),cwaveprj(:,:)

 !----------------------------------------------

 ! set up basic FFT parameters
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 istwf_k = 1
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me = mpi_enreg%me_kpt
 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0; ikg1 = 0

 ! Fermi energy
 call local_fermie(dtset,ebands_k,fermie,mpi_enreg)

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
 call getph(atindx,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,crystal%xred)

 ABI_MALLOC(kg_k,(3,mpw))
 ABI_MALLOC(kinpw,(mpw))

 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,0,dimlmn)

 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 lmn2max = psps%lmnmax*(psps%lmnmax+1)/2
 ! note: in make_d, terms will be filled as iatom using atindx
 call dterm_alloc(dterm,psps%lmnmax,lmn2max,dtset%natom,paw_ij(1)%ndij)
 call make_d(atindx,dterm,dtset,crystal%gprimd,paw_ij,pawrad,pawtab,psps)

 ! number of terms to store on the kpt mesh
 ! CC, VV1, VV2, NL, L_R, B.M
 call orbmag_mesh_alloc(dtset,orbmag_mesh)
 orbmag_mesh%nucdipmom=dtset%nucdipmom
 ! if user input lambsig specifically in the input file, use it
 if ( any ( abs(dtset%lambsig).GT.tol8 ) ) then
   orbmag_mesh%lambsig=dtset%lambsig
 ! else use the value read in to pawtab structure (which might well be zero)
 else
   orbmag_mesh%lambsig=pawtab(1:dtset%ntypat)%lamb_shielding
 end if

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k
 call gs_hamk%init(psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,crystal%xred,dtset%nfft,dtset%mgfft,dtset%ngfft,crystal%rprimd,&
      & dtset%nloalg,nucdipmom=dtset%nucdipmom,paw_ij=paw_ij)

 ! iterate over spin channels
 bdtot_index=0
 icg = 0
 icprj = 0
 do isppol = 1, dtset%nsppol

   !========= construct local potential ==================
   ABI_MALLOC(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
   call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, nfftf, &
     & dtset%nspden, gs_hamk%nvloc, 1, pawfgr, mpi_enreg, vtrial, vlocal)
   call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)

   !========  compute nuclear dipole vector potential (may be zero) ==========
   has_nucdip = ANY( ABS(dtset%nucdipmom) .GT. tol8 )
   if(has_nucdip) then
     ABI_MALLOC(vectornd,(nfftf,dtset%nspden,3))
     vectornd = zero
     call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,&
       & dtset%nspden,dtset%nucdipmom,crystal%rprimd,vectornd,crystal%xred)
     ABI_MALLOC(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
     call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, nfftf, &
          & dtset%nspden, gs_hamk%nvloc, 3, pawfgr, mpi_enreg, vectornd,vectornd_pac)
     ABI_FREE(vectornd)
     call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
   end if

   !========  compute vxctaulocal if vxctau present =====================

   if (usevxctau==1) then
     ABI_MALLOC(vxctaulocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,4))
     call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, nfftf, &
       & dtset%nspden, gs_hamk%nvloc, 4, pawfgr, mpi_enreg, vxctau, vxctaulocal)
     call gs_hamk%load_spin(isppol, vxctaulocal=vxctaulocal)
   end if

   ikg = 0
   !============= BIG FAT KPT LOOP :) ===========================
   do ikpt = 1, dtset%nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     nband_me = proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,nband_k,isppol,me)

     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     kpoint(:)=ebands_k%kptns(:,ikpt)
     npw_k = ebands_k%npwarr(ikpt)
     npwsp = npw_k*dtset%nspinor

     ! retrieve kg_k at this k point
     kg_k(1:3,1:npw_k) = kg(1:3,ikg+1:ikg+npw_k)

     ! retrieve ylm at this k point
     ABI_MALLOC(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
     ABI_MALLOC(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
     do ilm=1,psps%mpsang*psps%mpsang
       ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
     end do

     ! retrieve occupation numbers at this k point
     ABI_MALLOC(occ_k,(nband_k))
     !occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
     occ_k(:)=ebands_k%occ(1:nband_k,ikpt,isppol)

     ! Compute kinetic energy at kpt
     kinpw(:) = zero
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,crystal%gmet,kg_k,kinpw,kpoint,npw_k,0,0)

     ! Compute k+G at this k point
     nkpg = 3
     ABI_MALLOC(kpg_k,(npw_k,nkpg))
     call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

     ! Make 3d phase factors
     ABI_MALLOC(phkxred,(2,dtset%natom))
     do iat = 1, dtset%natom
       iatom = atindx(iat)
       arg=two_pi*DOT_PRODUCT(kpoint,crystal%xred(:,iat))
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
       & crystal%gmet,crystal%gprimd,ider,idir,psps%indlmn,&
       & kg_k,kpg_k,kpoint,psps%lmnmax,&
       & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
       & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,crystal%rmet,&
       & psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
     !  - Load k-dependent quantities in the Hamiltonian
     call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
       & kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,&
       & compute_gbound=.TRUE.)

     ! retrieve ground state wavefunctions at this k point and isppol
     mcgk = npw_k*nband_k*dtset%nspinor
     ABI_MALLOC(cg_k,(2,mcgk))
     cg_k = cg(1:2,icg+1:icg+mcgk)

     ! retrieve first order wavefunctions at this k point and isppol
     ABI_MALLOC(cg1_k,(2,mcgk,3))
     cg1_k = cg1(1:2,icg+1:icg+mcgk,1:3)

     ! retrieve zeroth order eigenvalues at this k point and isppol
     ABI_MALLOC(eig_k,(nband_k))
     !eig_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     eig_k(:)=ebands_k%eig(1:nband_k,ikpt,isppol)

     ! retrieve cprj_k at this k point and isppol
     mcprjk = nband_k*dtset%nspinor
     ABI_MALLOC(cprj_k,(dtset%natom,mcprjk))
     call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimlmn)
     call pawcprj_get(atindx,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,nband_k,nband_k,dtset%nspinor,dtset%nsppol,0)

     ! compute P_c|cg1>
     ABI_MALLOC(pcg1_k,(2,mcgk,3))
     if (dtset%berryopt /= -2) then
       ! DDK wavefunctions computed in DFPT, must be projected onto conduction space
       call make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,ikpt,isppol,&
         & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,pcg1_k)
     else
       ! we are using berryopt PEAD DDK functions which are already projected by construction
       pcg1_k(1:2,1:mcgk,1:3) = cg1_k(1:2,1:mcgk,1:3)
     end if

     ! compute <p|Pc cg1> cprjs
     ABI_MALLOC(cprj1_k,(dtset%natom,mcprjk,3))
     do adir = 1, 3
       call pawcprj_alloc(cprj1_k(:,:,adir),0,dimlmn)
     end do
     choice = 1
     cpopt = 0
     idir = 0
     ABI_MALLOC(cwavef,(2,npwsp))
     do nn = 1, nband_k
       do adir = 1, 3
         cwavef(1:2,1:npwsp) = pcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,adir)
         call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl_k,idir,psps%indlmn,istwf_k,&
           & kg_k,kpg_k,kpoint,psps%lmnmax,dtset%mgfft,mpi_enreg,1,dtset%natom,nattyp,dtset%ngfft,&
           & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,phkxred,ph1d,ph3d,crystal%ucvol,psps%useylm)
         call pawcprj_put(atindx,cwaveprj,cprj1_k(:,:,adir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & mkmem_rbz,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)
       end do
     end do
     ABI_FREE(cwavef)

     !--------------------------------------------------------------------------------
     ! Finally ready to compute contributions to orbital magnetism and Berry curvature
     !--------------------------------------------------------------------------------

     ! ZTG23 Eq. 36 term 2 and Eq. 46 term 1
     call orbmag_cc_k(atindx,cprj1_k,dimlmn,dtset,eig_k,fermie,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,pcg1_k)

     ! ZTG23 Eq. 36 terms 3 and 4 and Eq. 46 term 2
     call orbmag_vv_k(atindx,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gs_hamk,&
      & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,orbmag_mesh,pcg1_k)

     ! ZTG23 Eq. 36 term 1
     call orbmag_nl_k(atindx,cprj_k,dimlmn,dterm,dtset,eig_k,ikpt,isppol,&
       & mcprjk,mkmem_rbz,nband_k,orbmag_mesh,pawtab)

     ! ZTG23 text after Eq. 42
     nl1_option = 1 ! LR
     call orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,mcprjk,mkmem_rbz,&
       & nband_k,nl1_option,orbmag_mesh,pawtab)

     ! ZTG23 Eq. 43
     nl1_option = 2 ! BM
     call orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,mcprjk,mkmem_rbz,&
       & nband_k,nl1_option,orbmag_mesh,pawtab)

     icg = icg + mcgk
     icprj = icprj + mcprjk
     ikg = ikg + npw_k
     bdtot_index=bdtot_index+nband_k

     ABI_FREE(cg_k)
     ABI_FREE(cg1_k)
     ABI_FREE(pcg1_k)
     ABI_FREE(eig_k)
     ABI_FREE(occ_k)
     ABI_FREE(ylm_k)
     ABI_FREE(ylmgr_k)
     ABI_FREE(kpg_k)
     ABI_FREE(ffnl_k)
     ABI_FREE(ph3d)
     ABI_FREE(phkxred)
     call pawcprj_free(cprj_k)
     ABI_FREE(cprj_k)
     do adir = 1, 3
       call pawcprj_free(cprj1_k(:,:,adir))
     end do
     ABI_FREE(cprj1_k)

   end do ! end loop over kpts

   ABI_FREE(vlocal)
   if(allocated(vectornd_pac)) then
     ABI_FREE(vectornd_pac)
   end if
   if(allocated(vxctaulocal)) then
      ABI_FREE(vxctaulocal)
   end if

 end do ! end loop over isppol

 ! prepare orbmag_terms for output to abo file
 ABI_MALLOC(orbmag_terms,(dtset%mband,dtset%nsppol,3,orbmag_nterms))
 ABI_MALLOC(chern_terms,(dtset%mband,dtset%nsppol,3,chern_nterms))
 call sum_orbmag_mesh(chern_terms,crystal,dtset,ebands_k,mpi_enreg,orbmag_mesh,orbmag_terms)

 !! collect orbmag_mesh if distributed over different processes
 if (nproc > 1) then
   buff_size=size(orbmag_mesh%omesh)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = &
     & reshape(orbmag_mesh%omesh,(/orbmag_mesh%mband*orbmag_mesh%nkpt*orbmag_mesh%nsppol*3*orbmag_nterms/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   orbmag_mesh%omesh(1:orbmag_mesh%mband,1:orbmag_mesh%nkpt,1:orbmag_mesh%nsppol,1:3,1:orbmag_nterms)=&
     & reshape(buffer2,(/orbmag_mesh%mband,orbmag_mesh%nkpt,orbmag_mesh%nsppol,3,orbmag_nterms/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

 !! collect chern_mesh if distributed over different processes
 if (nproc > 1) then
   buff_size=size(orbmag_mesh%cmesh)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = &
     & reshape(orbmag_mesh%cmesh,(/orbmag_mesh%mband*orbmag_mesh%nkpt*orbmag_mesh%nsppol*3*chern_nterms/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   orbmag_mesh%cmesh(1:orbmag_mesh%mband,1:orbmag_mesh%nkpt,1:orbmag_mesh%nsppol,1:3,1:chern_nterms)=&
     & reshape(buffer2,(/orbmag_mesh%mband,orbmag_mesh%nkpt,orbmag_mesh%nsppol,3,chern_nterms/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(3,orbmag_nterms))
 ABI_MALLOC(chern_trace,(3,chern_nterms))
 orbmag_trace = zero
 chern_trace = zero
 do isppol = 1, dtset%nsppol
   do nn = 1, nband_k
     orbmag_trace(1:3,1:orbmag_nterms) = orbmag_trace(1:3,1:orbmag_nterms) + &
       & orbmag_terms(nn,isppol,1:3,1:orbmag_nterms)
     chern_trace(1:3,1:chern_nterms) = chern_trace(1:3,1:chern_nterms) + &
       & chern_terms(nn,isppol,1:3,1:chern_nterms)
   end do ! nn
 end do ! isppol

 ! get the Lamb term
 call lamb_core(atindx,dtset,omlamb,pawtab)

 ! output raw data to netcdf file for more detailed postprocessing
 if (me == master) then
   fname = trim(dtfil%filnam_ds(4))//'_ORBMAG.nc'
   NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))
   call orbmag_ncwrite(crystal,dtset,ebands_k,hdr,ncid,orbmag_mesh)
   NCF_CHECK(nf90_close(ncid))
 end if

 ! output summary to abo file
 call orbmag_output(chern_terms,chern_trace,dtset,omlamb,orbmag_terms,orbmag_trace)

!---------------------------------------------------
! deallocate memory
!---------------------------------------------------

 call gs_hamk%free()

 ABI_FREE(kg_k)
 ABI_FREE(kinpw)
 ABI_FREE(ph1d)

 ABI_FREE(realgnt)
 ABI_FREE(gntselect)

 ABI_FREE(atindx)
 ABI_FREE(atindx1)
 ABI_FREE(nattyp)

 ABI_FREE(dimlmn)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 call dterm_free(dterm)
 call orbmag_mesh_free(orbmag_mesh)

 ABI_FREE(orbmag_terms)
 ABI_FREE(chern_terms)
 ABI_FREE(orbmag_trace)
 ABI_FREE(chern_trace)

end subroutine orbmag
!!***

!!****f*m_orbmag/sum_orbmag_mesh
!! NAME
!! sum_orbmag_mesh
!!
!! FUNCTION
!! sum terms in orbmag_mesh
!!
!! INPUTS
!!
!! OUTPUT
!!   orbmag_terms(mband,nsppol,3,nterms)
!!
!! SOURCE

subroutine sum_orbmag_mesh(chern_terms,crystal,dtset,ebands_k,mpi_enreg,orbmag_mesh,orbmag_terms)

  !Arguments ------------------------------------
  !scalars
  type(crystal_t),intent(in) :: crystal
  type(dataset_type),intent(in) :: dtset
  type(ebands_t),intent(in) :: ebands_k
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(in) :: orbmag_mesh

  !arrays
  real(dp),intent(out) :: chern_terms(dtset%mband,dtset%nsppol,3,chern_nterms)
  real(dp),intent(out) :: orbmag_terms(dtset%mband,dtset%nsppol,3,orbmag_nterms)

  !Local variables -------------------------
  !scalars
  integer :: buff_size,ierr,ikpt,isppol,iterm,nband_k
  integer :: me,nn,nproc,spaceComm
  real(dp) :: trnrm
  real(dp),allocatable :: buffer1(:),buffer2(:)
  !arrays

!--------------------------------------------------------------------

  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  chern_terms = zero
  orbmag_terms = zero
  do isppol = 1, dtset%nsppol
    do ikpt = 1, dtset%nkpt
      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      ! if the current kpt is not on the current processor, cycle
      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

      do nn = 1, nband_k
        trnrm = ebands_k%occ(nn,ikpt,isppol)*dtset%wtk(ikpt)/crystal%ucvol
        if(abs(trnrm).LT.tol8) cycle

        chern_terms(nn,isppol,1:3,ibcc:ibvv2) = chern_terms(nn,isppol,1:3,ibcc:ibvv2) + &
            & trnrm*orbmag_mesh%cmesh(nn,ikpt,isppol,1:3,ibcc:ibvv2)

        orbmag_terms(nn,isppol,1:3,incc:inbm) = orbmag_terms(nn,isppol,1:3,incc:inbm) + &
            & trnrm*orbmag_mesh%omesh(nn,ikpt,isppol,1:3,incc:inbm)

      end do ! loop on bands
    end do ! loop on kpts
  end do ! loop on nsppol

  !! collect orbmag_terms if distributed over different processes
  if (nproc > 1) then
    buff_size=size(orbmag_terms)
    ABI_MALLOC(buffer1,(buff_size))
    ABI_MALLOC(buffer2,(buff_size))
    buffer1=zero;buffer2=zero
    buffer1(1:buff_size) = reshape(orbmag_terms,(/nband_k*dtset%nsppol*3*orbmag_nterms/))
    call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
    orbmag_terms(1:nband_k,1:dtset%nsppol,1:3,1:orbmag_nterms)=reshape(buffer2,(/nband_k,dtset%nsppol,3,orbmag_nterms/))
    ABI_FREE(buffer1)
    ABI_FREE(buffer2)
  end if

  !! collect chern_terms if distributed over different processes
  if (nproc > 1) then
    buff_size=size(chern_terms)
    ABI_MALLOC(buffer1,(buff_size))
    ABI_MALLOC(buffer2,(buff_size))
    buffer1=zero;buffer2=zero
    buffer1(1:buff_size) = reshape(chern_terms,(/nband_k*dtset%nsppol*3*chern_nterms/))
    call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
    chern_terms(1:nband_k,1:dtset%nsppol,1:3,1:chern_nterms)=reshape(buffer2,(/nband_k,dtset%nsppol,3,chern_nterms/))
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
 do iterm = 1, orbmag_nterms
   do isppol = 1, dtset%nsppol
     do nn = 1, nband_k
       if((iterm.EQ.inlr).OR.(iterm.EQ.inbm)) then
         orbmag_terms(nn,isppol,1:3,iterm) = &
           & MATMUL(crystal%rprimd,orbmag_terms(nn,isppol,1:3,iterm))
       else
         orbmag_terms(nn,isppol,1:3,iterm) = &
           & crystal%ucvol*MATMUL(crystal%gprimd,orbmag_terms(nn,isppol,1:3,iterm))
       end if
     end do ! nn
   end do !isppol
 end do

 do iterm = 1, chern_nterms
   do isppol = 1, dtset%nsppol
     do nn = 1, nband_k
       chern_terms(nn,isppol,1:3,iterm) = &
         & crystal%ucvol*MATMUL(crystal%gprimd,chern_terms(nn,isppol,1:3,iterm))
     end do ! nn
   end do !isppol
 end do

 !! convert orbmag magnetization to orbital moment
 !! Berry curvature terms are ignored
 orbmag_terms(:,:,:,incc:inbm)=crystal%ucvol*orbmag_terms(:,:,:,incc:inbm)

end subroutine sum_orbmag_mesh
!!***


!!****f* ABINIT/orbmag_nl1_k
!! NAME
!! orbmag_nl1_k
!!
!! FUNCTION
!! make NL(1) term at k
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cprj_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for cg_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dterm <type(dterm_type)> data related to onsite interactions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  nband_k=bands at this kpt
!!  nl1_option=chooses which onsite term to apply
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  orbmag_mesh%omesh
!!  if nl1_option = 1, orbmag contribution of <L_R> is returned
!!  if nl1_option = 2, orbmag contribution of <A0.An> is returned
!!
!! TODO
!!
!! NOTES
!!  returns \sum_{Rij}<u|p_i>a_ij<p_j|u> for various a_ij inputs
!!
!! SOURCE

subroutine orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,mcprjk,&
    & mkmem_rbz,nband_k,nl1_option,orbmag_mesh,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcprjk,mkmem_rbz,nband_k,nl1_option
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,nn
  complex(dp) :: tt
  !arrays
  type(pawcprj_type),allocatable :: cwaveprj(:,:)
!--------------------------------------------------------------------

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do adir = 1, 3
   do nn = 1, nband_k

     call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

     select case (nl1_option)
     case(1)
       call tt_me(dterm%LR(:,:,:,adir),atindx,cwaveprj,dtset,cwaveprj,dterm%lmn2max,dterm%ndij,pawtab,tt)
       orbmag_mesh%omesh(nn,ikpt,isppol,adir,inlr) = real(tt)
     case(2)
       call tt_me(dterm%BM(:,:,:,adir),atindx,cwaveprj,dtset,cwaveprj,dterm%lmn2max,dterm%ndij,pawtab,tt)
       orbmag_mesh%omesh(nn,ikpt,isppol,adir,inbm) = real(tt)
     case default
       tt = czero
     end select

   end do !nn
 end do !adir

 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

end subroutine orbmag_nl1_k
!!***


!!****f* ABINIT/orbmag_nl_k
!! NAME
!! orbmag_nl_k
!!
!! FUNCTION
!! make NL term at k
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cprj_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for cg_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dterm <type(dterm_type)> data related to onsite interactions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(nband_k)=gs eigenvalues at this kpt
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  nband_k=bands at this kpt
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  orbmag_mesh%omesh
!!
!! TODO
!!
!! NOTES
!! computes -\frac{i}{2}\sum_{Rij}<u|d_b p_i>D^0_{ij} - E^0s^0_{ij}<d_g p_j|u>
!! This is ZTG23 Eq. 36 term 1
!!
!! SOURCE

subroutine orbmag_nl_k(atindx,cprj_k,dimlmn,dterm,dtset,eig_k,ikpt,isppol,&
    & mcprjk,mkmem_rbz,nband_k,orbmag_mesh,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcprjk,mkmem_rbz,nband_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,nn
  real(dp) :: epsabg
  complex(dp) :: m1,prefac_m,txt_d,txt_q
  !arrays
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do adir = 1, 3
   do nn = 1, nband_k
     call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

     m1 = czero
     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_m = com*c2*epsabg

         call txt_me(dterm%aij,atindx,cwaveprj,bdir,dtset,gdir,cwaveprj,&
           & dterm%lmn2max,dterm%ndij,pawtab,txt_d)

         call txt_me(dterm%qij,atindx,cwaveprj,bdir,dtset,gdir,cwaveprj,&
           & dterm%lmn2max,dterm%ndij,pawtab,txt_q)

         ! note rho^0 H^1 term has opposite sign of rho^1 H^0
         m1 = m1 - prefac_m*(txt_d - eig_k(nn)*txt_q)

       end do !gdir
     end do !bdir

     orbmag_mesh%omesh(nn,ikpt,isppol,adir,innl) = real(m1)

   end do !nn
 end do !adir

 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

end subroutine orbmag_nl_k
!!***

!!****f* ABINIT/orbmag_cc_k
!! NAME
!! orbmag_cc_k
!!
!! FUNCTION
!! computes <P_c du/dk|H + E*S|P_c du/dk> term in orbital magnetism
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cprj1_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for Pc cg1_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(nband_k)=gs eigenvalues at this kpt
!!  fermie=offset energy to use
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcgk=dimension of cg_k
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  nband_k=bands at this kpt
!!  npw_k=number of planewaves at this kpt
!!  occ_k=band occupations at this kpt
!!  pcg1_k(2,mcgk,3)=cg1_k projected on conduction space
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  orbmag_mesh%omesh, orbmag_mesh%cmesh
!!
!! TODO
!!
!! NOTES
!! ZTG23 Eq. 36 term 2 and Eq. 46 term 1
!!
!! SOURCE

subroutine orbmag_cc_k(atindx,cprj1_k,dimlmn,dtset,eig_k,fermie,gs_hamk,ikpt,isppol,&
    & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),pcg1_k(2,mcgk,3)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,mcprjk,3)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,gdir,ndat,nn,npwsp,sij_opt,tim_getghc,type_calc
  real(dp) :: doti,dotr,epsabg,lams
  complex(dp) :: b1,m1,m1_mu,prefac_b,prefac_m
  !arrays
  real(dp),allocatable :: bra(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),ket(:,:)
  type(pawcprj_type),allocatable :: cwaveprj1(:,:)
!--------------------------------------------------------------------

 npwsp = npw_k*dtset%nspinor

 ABI_MALLOC(bra,(2,npwsp))
 ABI_MALLOC(ket,(2,npwsp))
 ABI_MALLOC(ghc,(2,npwsp))
 ABI_MALLOC(gsc,(2,npwsp))
 ABI_MALLOC(gvnlxc,(2,npwsp))
 ABI_MALLOC(cwaveprj1,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj1,0,dimlmn)

 tim_getghc = 0
 lams = zero
 ndat = 1
 sij_opt = 1
 type_calc = 0

 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero
     m1_mu = czero
     b1 = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg

         cpopt = 2
         ket(1:2,1:npwsp) = pcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,gdir)

         call pawcprj_get(atindx,cwaveprj1,cprj1_k(:,:,gdir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

         ! compute H|Pc du> and S|Pc du>
         call getghc(cpopt,ket,cwaveprj1,ghc,gsc,gs_hamk,gvnlxc,lams,mpi_enreg,&
           & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)

         bra(1:2,1:npwsp) = pcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,bdir)

         dotr = DOT_PRODUCT(bra(1,:),ghc(1,:))+DOT_PRODUCT(bra(2,:),ghc(2,:))
         doti = DOT_PRODUCT(bra(1,:),ghc(2,:))-DOT_PRODUCT(bra(2,:),ghc(1,:))
         m1 = m1 + prefac_m*CMPLX(dotr,doti)

         dotr = DOT_PRODUCT(bra(1,:),gsc(1,:))+DOT_PRODUCT(bra(2,:),gsc(2,:))
         doti = DOT_PRODUCT(bra(1,:),gsc(2,:))-DOT_PRODUCT(bra(2,:),gsc(1,:))
         b1 = b1 -two*prefac_b*CMPLX(dotr,doti)
         m1 = m1 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)
         m1_mu = m1_mu - two*prefac_m*CMPLX(dotr,doti)*fermie

       end do !gdir
     end do !bdir

     orbmag_mesh%omesh(nn,ikpt,isppol,adir,incc) = real(m1 + m1_mu)
     orbmag_mesh%cmesh(nn,ikpt,isppol,adir,ibcc) = real(b1)

   end do !nn
 end do !adir

 ABI_FREE(bra)
 ABI_FREE(ket)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(gvnlxc)
 call pawcprj_free(cwaveprj1)
 ABI_FREE(cwaveprj1)

end subroutine orbmag_cc_k
!!***

!!****f* ABINIT/orbmag_vv_k
!! NAME
!! orbmag_vv_k
!!
!! FUNCTION
!! orbmag_vv_k
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg_k(2,mcgk)=ground state wavefunctions at this k point
!!  cprj_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for cg_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(nband_k)=gs eigenvalues at this kpt
!!  fermie=offset energy to use
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcgk=dimension of cg_k
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  nband_k=bands at this kpt
!!  npw_k=number of planewaves at this kpt
!!  pcg1_k(2,mcgk,3)=cg1_k projected on conduction space
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  orbmag_mesh%omesh, orbmag_mesh%cmesh
!!
!! TODO
!!
!! NOTES
!! contributions (1) <Pc d_b u|E d_gS|u> + <u|E d_bS|Pc d_g u> and
!! (2) \sum_n' <u |d_b ES|u_n'><u_n'|d_g ES|u> to orbital magnetization
!! these are ZTG23 Eq 36 terms 3 and 4, and Eq. 46 term 2
!!
!! SOURCE

subroutine orbmag_vv_k(atindx,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gs_hamk,&
    & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,orbmag_mesh,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k),pcg1_k(2,mcgk,3),occ_k(nband_k)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,gdir,ndat,nn,nnlout,np,npwsp
  integer :: paw_opt,signs,tim_getghc
  real(dp) :: doti,dotr,epsabg
  complex(dp) :: b1,bv2b,m1,m1_mu,mb,mg,mv2b,mv2b_mu,prefac_b,prefac_m
  !arrays
  real(dp) :: enlout(1),lamv(1)
  real(dp),allocatable :: bra(:,:),ket(:,:),svectoutb(:,:),svectoutg(:,:),vectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)
!--------------------------------------------------------------------

 npwsp = npw_k*dtset%nspinor

 ABI_MALLOC(bra,(2,npwsp))
 ABI_MALLOC(ket,(2,npwsp))
 ABI_MALLOC(svectoutb,(2,npwsp))
 ABI_MALLOC(svectoutg,(2,npwsp))
 ABI_MALLOC(vectout,(2,npwsp))
 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 tim_getghc = 0
 lamv = zero
 ndat = 1
 cpopt = 4 ! cprj and derivs in memory
 choice = 5 ! apply dS/dk
 paw_opt = 3 ! retain dS/dk|u>
 signs = 2
 nnlout = 1

 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero; mv2b = czero
     m1_mu = czero; mv2b_mu = czero
     b1 = czero; bv2b = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg

         ! extract |u_nk>
         ket(1:2,1:npwsp) = cg_k(1:2,(nn-1)*npwsp+1:nn*npwsp)
         call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

         ! compute dS/dk_g |u_nk>
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,gdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutg,tim_getghc,ket,vectout)
         ! compute dS/dk_b|u_nk>
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,bdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutb,tim_getghc,ket,vectout)

         ! extract |Pc du/dk_b>
         bra(1:2,1:npwsp) = pcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,bdir)
         !overlap
         dotr = DOT_PRODUCT(bra(1,:),svectoutg(1,:))+DOT_PRODUCT(bra(2,:),svectoutg(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectoutg(2,:))-DOT_PRODUCT(bra(2,:),svectoutg(1,:))
         ! add <Pc du/dk_b|dS/dk_g|u_nk>*E_nk
         b1 = b1 - prefac_b*CMPLX(dotr,doti)
         m1 = m1 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)
         m1_mu = m1_mu - prefac_m*CMPLX(dotr,doti)*fermie

         ! extract |Pc du/dk_g>
         bra(1:2,1:npwsp) = pcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,gdir)
         !overlap
         dotr = DOT_PRODUCT(bra(1,:),svectoutb(1,:))+DOT_PRODUCT(bra(2,:),svectoutb(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectoutb(2,:))-DOT_PRODUCT(bra(2,:),svectoutb(1,:))

         ! add CONJG(<Pc du/dk_b|dS/dk_g|u_nk>)*E_nk
         b1 = b1 - prefac_b*CMPLX(dotr,-doti)
         m1 = m1 + prefac_m*CMPLX(dotr,-doti)*eig_k(nn)
         m1_mu = m1_mu - prefac_m*CMPLX(dotr,-doti)*fermie

         do np = 1, nband_k

           if (occ_k(np).LT.tol8) cycle

           bra(1:2,1:npwsp) = cg_k(1:2,(np-1)*npwsp+1:np*npwsp)

           dotr = DOT_PRODUCT(bra(1,:),svectoutg(1,:))+DOT_PRODUCT(bra(2,:),svectoutg(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutg(2,:))-DOT_PRODUCT(bra(2,:),svectoutg(1,:))
           mg = CMPLX(dotr,doti)

           dotr = DOT_PRODUCT(bra(1,:),svectoutb(1,:))+DOT_PRODUCT(bra(2,:),svectoutb(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutb(2,:))-DOT_PRODUCT(bra(2,:),svectoutb(1,:))
           mb = CMPLX(dotr,doti)

           ! terms in <u|dS|u'><u'|dS|u>
           bv2b = bv2b + prefac_b*CONJG(mb)*mg
           mv2b = mv2b - prefac_m*CONJG(mb)*mg*eig_k(nn)
           mv2b_mu = mv2b_mu + prefac_m*CONJG(mb)*mg*fermie
         end do ! np

       end do !gdir
     end do !bdir

     orbmag_mesh%cmesh(nn,ikpt,isppol,adir,ibvv1) = real(b1)
     orbmag_mesh%cmesh(nn,ikpt,isppol,adir,ibvv2) = real(bv2b)

     orbmag_mesh%omesh(nn,ikpt,isppol,adir,invv1) = real(m1+m1_mu)
     orbmag_mesh%omesh(nn,ikpt,isppol,adir,invv2) = real(mv2b+mv2b_mu)

   end do !nn
 end do !adir

 ABI_FREE(bra)
 ABI_FREE(ket)
 ABI_FREE(svectoutb)
 ABI_FREE(svectoutg)
 ABI_FREE(vectout)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

end subroutine orbmag_vv_k
!!***

!!****f* ABINIT/make_pcg1
!! NAME
!! make_pcg1
!!
!! FUNCTION
!! compute Pc|cg1> from |cg1> and |cg>
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg_k(2,mcgk)=ground state wavefunctions at this k point
!!  cg1_k(2,mcgk,3)=DDK wavefunctions at this k point, all 3 directions
!!  cprj_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for cg_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcgk=dimension of cg_k
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  nband_k=bands at this kpt
!!  npw_k=number of planewaves at this kpt
!!  occ_k=band occupations at this kpt
!!
!! OUTPUT
!!  pcg1_k(2,mcgk,3)=cg1_k projected on conduction space
!!
!! NOTES
!! see Audouze et al PRB 78, 035105 (2008) Eq. 40
!!
!! SOURCE

subroutine make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,&
    & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),cg1_k(2,mcgk,3),occ_k(nband_k)
  real(dp),intent(out) :: pcg1_k(2,mcgk,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars
  integer :: adir,choice,cpopt,iband,jband
  integer :: ndat,nnlout,npwsp,paw_opt,signs,tim_nonlop
  real(dp) :: doti,dotr
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

  npwsp = npw_k*dtset%nspinor

  ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
  call pawcprj_alloc(cwaveprj,3,dimlmn)
  ABI_MALLOC(cwavef,(2,npwsp))
  ABI_MALLOC(vectout,(2,npwsp))
  ABI_MALLOC(svectout,(2,npwsp))
  ABI_MALLOC(vcg1,(2,npwsp))

  pcg1_k = zero

  do adir = 1, 3

    do iband = 1, nband_k

      cwavef(1:2,1:npwsp)=cg_k(1:2,(iband-1)*npwsp+1:iband*npwsp)

      call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,iband,0,ikpt,0,isppol,dtset%mband,&
        & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

      ! compute S^1|u_i^0> where S^1 = \partial S/\partial k_adir, the k derivative of S in
      ! direction adir
      call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,adir,lambda,mpi_enreg,ndat,&
        & nnlout,paw_opt,signs,svectout,tim_nonlop,cwavef,vectout)

      !! form vcg1 = -1/2 \sum |u_j^0><u_j^0|S^1|u_i^0>, the valence band part of cg1
      vcg1 = zero
      do jband = 1, nband_k
        if(abs(occ_k(jband)).LT.tol8) cycle
        cwavef(1:2,1:npwsp)=cg_k(1:2,(jband-1)*npwsp+1:jband*npwsp)
        dotr = DOT_PRODUCT(cwavef(1,:),svectout(1,:))+DOT_PRODUCT(cwavef(2,:),svectout(2,:))
        doti = DOT_PRODUCT(cwavef(1,:),svectout(2,:))-DOT_PRODUCT(cwavef(2,:),svectout(1,:))
        vcg1(1,:) = vcg1(1,:) - half*( dotr*cwavef(1,:) - doti*cwavef(2,:))
        vcg1(2,:) = vcg1(2,:) - half*( dotr*cwavef(2,:) + doti*cwavef(1,:))
      end do

      ! subtract vcg1 from cg1_k to obtain pcg1, the conduction band part of cg1
      pcg1_k(1:2,(iband-1)*npwsp+1:iband*npwsp,adir) =cg1_k(1:2,(iband-1)*npwsp+1:iband*npwsp,adir)-&
        &  vcg1(1:2,1:npwsp)

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
!! INPUTS
!!  atindx(dtset%natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  omlamb(2,3)=contribution of Lamb shielding to magnetic moment
!!
!! NOTES
!!  lamb shielding of core electrons contributes -m.lambsig to orbital magnetic
!!  moment
!!
!! SOURCE

subroutine lamb_core(atindx,dtset,omlamb,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: omlamb(3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,itypat
  real(dp) :: lambsig

!--------------------------------------------------------------------

  omlamb = zero
  do iat=1,dtset%natom
    iatom = atindx(iat)
    itypat = dtset%typat(iat)
    ! if user input lambsig specifically in the input file, use it
    if (abs(dtset%lambsig(itypat)).GT.tol8) then
      lambsig=dtset%lambsig(itypat)
    ! else use the value read in to pawtab structure (which might well be zero)
    else
      lambsig=pawtab(itypat)%lamb_shielding
    end if
    do adir = 1, 3
      omlamb(adir) = omlamb(adir) - lambsig*dtset%nucdipmom(adir,iat)
    end do ! end loop over adir
  end do ! end loop over iat

end subroutine lamb_core
!!***

!!****f* ABINIT/txt_me
!! NAME
!! txt_me
!!
!! FUNCTION
!! Onsite part of matrix element <u_n|dp>a_ij<dp|u_m>
!!
!! INPUTS
!!  aij(dtset%natom,lmn2max,ndij)=(complex)scalar ij couplings
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  bcp(dtset%natom,dtset%nspinor)<type(pawcprj_type)>=bra side cprj <p|bra>
!!  bdir=direction of bcp derivative to use
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gdir=direction of kcp derivative to use
!!  kcp(dtset%natom,dtset%nspinor)<type(pawcprj_type)>=ket side cprj <p|ket>
!!  lmn2max=max value of lmn2 over all psps
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ndij=spin channels in dij
!!
!! OUTPUT
!! txt=(complex) computed matrix element
!!
!! OUTPUT
!!
!! NOTES
!! computes on-site \sum_{Rij}<u|d_bdir p_i>aij<d_gdir p_j|u> for generic aij input
!!
!! SOURCE

subroutine txt_me(aij,atindx,bcp,bdir,dtset,gdir,kcp,lmn2max,ndij,pawtab,txt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,gdir,lmn2max,ndij
  complex(dp),intent(out) :: txt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dp),intent(in) :: aij(dtset%natom,lmn2max,ndij)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom,dtset%nspinor)
  type(pawcprj_type),intent(in) :: kcp(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,isp,jlmn,klmn
  complex(dp) :: dcpi,dcpj,dij
!--------------------------------------------------------------------

  txt = czero
  do iat = 1, dtset%natom
    itypat=dtset%typat(iat)
    iatom = atindx(iat)
    do isp = 1, dtset%nspinor
      do ilmn = 1, pawtab(itypat)%lmn_size
        do jlmn = 1, pawtab(itypat)%lmn_size
          klmn = MATPACK(ilmn,jlmn)
          ! in ndij = 4 case, isp 1 delivers up-up, isp 2 delivers down-down
          dij = aij(iatom,klmn,isp)
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          dcpi = CMPLX(bcp(iatom,isp)%dcp(1,bdir,ilmn),bcp(iatom,isp)%dcp(2,bdir,ilmn))
          dcpj = CMPLX(kcp(iatom,isp)%dcp(1,gdir,jlmn),kcp(iatom,isp)%dcp(2,gdir,jlmn))
          txt = txt + CONJG(dcpi)*dcpj*dij
          if (ndij == 4) then
            if (isp == 1) then
              dij = aij(iatom,klmn,3) ! up-down
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,4))
              dcpi = CMPLX(bcp(iatom,1)%dcp(1,bdir,ilmn),bcp(iatom,1)%dcp(2,bdir,ilmn))
              dcpj = CMPLX(kcp(iatom,2)%dcp(1,gdir,jlmn),kcp(iatom,2)%dcp(2,gdir,jlmn))
            else
              dij = aij(iatom,klmn,4) ! down-up
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,3))
              dcpi = CMPLX(bcp(iatom,2)%dcp(1,bdir,ilmn),bcp(iatom,2)%dcp(2,bdir,ilmn))
              dcpj = CMPLX(kcp(iatom,1)%dcp(1,gdir,jlmn),kcp(iatom,1)%dcp(2,gdir,jlmn))
            end if
            txt = txt + CONJG(dcpi)*dcpj*dij
          end if
        end do !jlmn
      end do !ilmn
    end do ! isp
  end do !iat

end subroutine txt_me
!!***


!!****f* ABINIT/tt_me
!! NAME
!! tt_me
!!
!! FUNCTION
!! Onsite part of matrix element <u_n|a_ij|u_m>
!!
!! INPUTS
!!  aij(dtset%natom,lmn2max,ndij)=(complex)scalar ij couplings
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  bcp(dtset%natom,dtset%nspinor)<type(pawcprj_type)>=bra side cprj <p|bra>
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  kcp(dtset%natom,dtset%nspinor)<type(pawcprj_type)>=ket side cprj <p|ket>
!!  lmn2max=max value of lmn2 over all psps
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ndij=spin channels in dij
!!
!! OUTPUT
!! tt=(complex) computed matrix element
!!
!! NOTES
!! computes on-site \sum_{Rij}<u|p_i>aij<p_j|u> for generic aij input
!!
!! SOURCE

subroutine tt_me(aij,atindx,bcp,dtset,kcp,lmn2max,ndij,pawtab,tt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,ndij
  complex(dp),intent(out) :: tt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dp),intent(in) :: aij(dtset%natom,lmn2max,ndij)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom,dtset%nspinor),kcp(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,isp,itypat,ilmn,jlmn,klmn
  complex(dp) :: cpi,cpj,dij
!--------------------------------------------------------------------

  tt = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat=dtset%typat(iat)
    do isp = 1, dtset%nspinor
      do ilmn = 1, pawtab(itypat)%lmn_size
        do jlmn = 1, pawtab(itypat)%lmn_size
          klmn=MATPACK(ilmn,jlmn)
          ! in ndij = 4 case, isp 1 delivers up-up, isp 2 delivers down-down
          dij = aij(iatom,klmn,isp)
          cpi =  CMPLX(bcp(iatom,isp)%cp(1,ilmn),bcp(iatom,isp)%cp(2,ilmn))
          cpj =  CMPLX(kcp(iatom,isp)%cp(1,jlmn),kcp(iatom,isp)%cp(2,jlmn))
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          ! note use of CONJG(cpi), because cpi is from the bra side cprj
          tt = tt + CONJG(cpi)*dij*cpj
          if (ndij == 4) then
            if (isp == 1) then
              dij = aij(iatom,klmn,3) ! up-down
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,4))
              cpi = CMPLX(bcp(iatom,1)%cp(1,ilmn),bcp(iatom,1)%cp(2,ilmn))
              cpj = CMPLX(kcp(iatom,2)%cp(1,jlmn),kcp(iatom,2)%cp(2,jlmn))
            else
              dij = aij(iatom,klmn,4) ! down-up
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,3))
              cpi = CMPLX(bcp(iatom,2)%cp(1,ilmn),bcp(iatom,2)%cp(2,ilmn))
              cpj = CMPLX(kcp(iatom,1)%cp(1,jlmn),kcp(iatom,1)%cp(2,jlmn))
            end if
            tt = tt + CONJG(cpi)*cpj*dij
          end if
        end do !jlmn
      end do !ilmn
    end do ! isp
  end do !iat

end subroutine tt_me
!!***

!!****f* ABINIT/dterm_qij
!! NAME
!! dterm_qij
!!
!! FUNCTION
!! Transfer pawtab%sij to dterm, as complex, solely for convenience
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dterm <type(dterm_type)> data related to onsite interactions
!!
!! NOTES
!! Transfer pawtab%sij to dterm, as complex, solely for convenience
!!
!! SOURCE

subroutine dterm_qij(atindx,dterm,dtset,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,lmn2_size

  !arrays

!--------------------------------------------------------------------

 dterm%qij = czero
 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
   lmn2_size = pawtab(itypat)%lmn2_size
   dterm%qij(iatom,1:lmn2_size,1) = &
     & CMPLX(pawtab(itypat)%sij(1:lmn2_size),zero)
   if (dterm%ndij > 1) then
     dterm%qij(iatom,1:lmn2_size,2) = &
       & CMPLX(pawtab(itypat)%sij(1:lmn2_size),zero)
   end if
 end do ! iat

 dterm%has_qij=2

end subroutine dterm_qij
!!***

!!****f* ABINIT/dterm_BM
!! NAME
!! dterm_BM
!!
!! FUNCTION
!! Compute onsite <A0.AN>
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)=nonzero gaunt integral indices
!!  gprimd(3,3)=reciprocal space lattice vectors
!!  my_lmax=augmented l_max over all psp
!!  pawrad(dtset%ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  realgnt((2*my_lmax-1)**2*(my_lmax)**4)=nonzero gaunt integral values
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dterm <type(dterm_type)> data related to onsite interactions
!!
!! NOTES
!! ZTG23 Eq. 43
!! this term is A0.An = \frac{1}{2}(B x r).\alpha^2(m x r) which can be rewritten
!! as \frac{\alpha^2}{2} [B.(1-\hat{r}\hat{r}).m]/r .
!!
!! SOURCE

subroutine dterm_BM(atindx,dterm,dtset,gntselect,gprimd,my_lmax,pawrad,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: my_lmax
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,itypat,gs1,gs2
  integer :: klmn,klm,kln,mdir,mesh_size,ngnt,pwave_size
  real(dp) :: a2,afact,intg

  !arrays
  complex(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: dyadic(:,:,:),ff(:),radint(:)
!--------------------------------------------------------------------

  dterm%BM = czero
  a2 = FineStructureConstant2

  do iat = 1, dtset%natom
    if(.NOT. ANY(ABS(dtset%nucdipmom(:,iat))>tol8) ) cycle

    iatom = atindx(iat)
    itypat = dtset%typat(iat)

    mesh_size=pawtab(itypat)%mesh_size
    pwave_size=size(pawtab(itypat)%phiphj(:,1))

    ! compute angular integrals of S_i (1-rr) S_j
    gs1=size(gntselect,1)
    gs2=size(gntselect,2)
    ngnt=size(realgnt)
    ABI_MALLOC(dyadic,(3,3,gs2))
    call make_dyadic(one,one,dyadic,gntselect,gs1,gs2,gs2,ngnt,realgnt)

    ! compute radial integrals of (ui*uj - tilde{ui}tilde{uj})/r
    ABI_MALLOC(radint,(pawtab(itypat)%ij_size))
    ABI_MALLOC(ff,(mesh_size))
    do kln=1,pawtab(itypat)%ij_size
      ff(2:pwave_size) = &
        & (pawtab(itypat)%phiphj(2:pwave_size,kln)-pawtab(itypat)%tphitphj(2:pwave_size,kln))/&
        & (pawrad(itypat)%rad(2:pwave_size))
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))
      radint(kln)=intg
    end do
    ABI_FREE(ff)

    do klmn=1, pawtab(itypat)%lmn2_size
      klm  = pawtab(itypat)%indklmn(1,klmn)
      kln  = pawtab(itypat)%indklmn(2,klmn)

      dij_cart=zero
      do adir = 1, 3 ! B field direction
        do mdir = 1, 3 ! mag dipole direction
          afact=half*a2*radint(kln)*dyadic(adir,mdir,klm)*dtset%nucdipmom(mdir,iat)
          dij_cart(adir)=dij_cart(adir)-CMPLX(afact,zero)
        end do
      end do

      dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

      dterm%BM(iatom,klmn,1,1:3) = dij_red(1:3)
      if (dterm%ndij > 1) then
        dterm%BM(iatom,klmn,2,1:3) = dij_red(1:3)
      end if

    end do ! end loop over klmn

    ABI_FREE(dyadic)
    ABI_FREE(radint)
  end do ! end loop over iatom

  dterm%has_BM = 2

end subroutine dterm_BM
!!***

!!****f* ABINIT/dterm_LR
!! NAME
!! dterm_LR
!!
!! FUNCTION
!! Compute onsite <L_R/2>
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=reciprocal space lattice vectors
!!  pawrad(dtset%ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dterm <type(dterm_type)> data related to onsite interactions
!!
!! NOTES
!! ZTG23 text after Eq 42, the on-site angular momentum
!!
!! SOURCE

subroutine dterm_LR(atindx,dterm,dtset,gprimd,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,ilmn,il,im,itypat,jlmn,jl,jm
  integer :: klmn,kln,mesh_size,pwave_size
  real(dp) :: intg
  complex(dp) :: orbl_me
  !arrays
  complex(dp) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)
!--------------------------------------------------------------------

  dterm%LR = czero

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
      ! compute <L_dir>/2
        call slxyzs(il,im,adir,jl,jm,orbl_me)
        dij_cart(adir) = -half*orbl_me*intg
      end do ! end loop over adir

      ! convert to crystal frame
      dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

      do iat=1,dtset%natom
        iatom = atindx(iat)
        if(dtset%typat(iat) .EQ. itypat) then
          dterm%LR(iatom,klmn,1,1:3) = dij_red(1:3)
          if (dterm%ndij > 1) then
            dterm%LR(iatom,klmn,2,1:3) = dij_red(1:3)
          end if
        end if
      end do
    end do ! end loop over klmn
    ABI_FREE(ff)
  end do ! end loop over itypat

  dterm%has_LR = 2

end subroutine dterm_LR
!!***

!!****f* ABINIT/orbmag_output
!! NAME
!! orbmag_output
!!
!! FUNCTION
!! Only printing. This routine outputs orbmag terms to the normal abi out file
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  orbmag_terms(dtset%mband,dtset%nsppol,3,nterms)=all computed terms per band of orb mag
!!  orbmag_trace(3,nterms)=traces of orbmag_terms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine orbmag_output(chern_terms,chern_trace,dtset,omlamb,orbmag_terms,orbmag_trace)

 !Arguments ------------------------------------
 !scalars
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: omlamb(3)
 real(dp),intent(in) :: chern_terms(dtset%mband,dtset%nsppol,3,chern_nterms),chern_trace(3,chern_nterms)
 real(dp),intent(in) :: orbmag_terms(dtset%mband,dtset%nsppol,3,orbmag_nterms),orbmag_trace(3,orbmag_nterms)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,isppol,iterms
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(dtset%mband,3),berry_total(3),orbmag_bb(dtset%mband,3),orbmag_total(3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = 1,orbmag_nterms
   orbmag_total(1:3)=orbmag_total(1:3) + orbmag_trace(1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       orbmag_bb(iband,1:3) = orbmag_bb(iband,1:3) + orbmag_terms(iband,isppol,1:3,iterms)
     end do ! iband
   end do ! isppol
 end do

 orbmag_total=orbmag_total+omlamb

 berry_bb=zero;berry_total=zero
 do iterms = 1,chern_nterms
   berry_total(1:3)=berry_total(1:3) + chern_trace(1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       berry_bb(iband,1:3) = berry_bb(iband,1:3) + chern_terms(iband,isppol,1:3,iterms)
     end do ! iband
   end do ! isppol
 end do

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a,a)')' Orbital magnetic moment computed with DFPT derivative wavefunctions ',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a)')' Orbital magnetic moment, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (orbmag_total(adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')
 write(message,'(a)')ch10
 call wrtout(ab_out,message,'COLL')
 write(message,'(a)')' Chern vector, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (berry_total(adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .EQ. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(1) CC : ',(orbmag_trace(adir,incc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(orbmag_trace(adir,invv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(orbmag_trace(adir,invv2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(0) NL : ',(orbmag_trace(adir,innl),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <L_R> terms : ',(orbmag_trace(adir,inlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <A0.An> terms : ',(orbmag_trace(adir,inbm),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    Lamb terms : ',(omlamb(adir),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Chern vector, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  Ch CC : ',(chern_trace(adir,ibcc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv1 : ',(chern_trace(adir,ibvv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv2 : ',(chern_trace(adir,ibvv2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
 end if

 if(dtset%orbmag .EQ. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Term-by-term breakdowns for each band : '
   call wrtout(ab_out,message,'COLL')
   do iband = 1, dtset%mband
     do isppol = 1, dtset%nsppol
       write(message,'(a)')ch10
       call wrtout(ab_out,message,'COLL')
       if (dtset%nsppol .EQ. 2) then
         write(message,'(a,i3,a,i3,a,i2,a,i2)') ' band ',iband,' of ',dtset%mband,'; spin polarization ',isppol,' of ',dtset%nsppol
       else
         write(message,'(a,i3,a,i3)') ' band ',iband,' of ',dtset%mband
       end if
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Orbital magnetic moment : ',(orbmag_bb(iband,adir),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '     rho(1) CC : ',(orbmag_terms(iband,isppol,adir,incc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(orbmag_terms(iband,isppol,adir,invv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(orbmag_terms(iband,isppol,adir,invv2),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '     rho(0) NL : ',(orbmag_terms(iband,isppol,adir,innl),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '   <L_R> terms : ',(orbmag_terms(iband,isppol,adir,inlr),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' <A0.An> terms : ',(orbmag_terms(iband,isppol,adir,inbm),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a)')ch10
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Chern vector : ',(berry_bb(iband,adir),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '  Ch CC : ',(chern_terms(iband,isppol,adir,ibcc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV1 : ',(chern_terms(iband,isppol,adir,ibvv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV2 : ',(chern_terms(iband,isppol,adir,ibvv2),adir=1,3)
       call wrtout(ab_out,message,'COLL')
     end do !isppol
   end do ! iband
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_output
!!***

!!****f* ABINIT/dterm_free
!! NAME
!! dterm_free
!!
!! FUNCTION
!! free space in dterm_type
!!
!! SOURCE

subroutine dterm_free(dterm)

  !Arguments ------------------------------------
  !scalars
  class(dterm_type),intent(inout) :: dterm
!--------------------------------------------------------------------

  ABI_SFREE(dterm%aij)
  dterm%has_aij=0

  ABI_SFREE(dterm%qij)
  dterm%has_qij=0

  ABI_SFREE(dterm%LR)
  dterm%has_LR=0

  ABI_SFREE(dterm%BM)
  dterm%has_BM=0

end subroutine dterm_free
!!***

!!****f* ABINIT/dterm_alloc
!! NAME
!! dterm_alloc
!!
!! FUNCTION
!! allocate space in dterm_type
!!
!! INPUTS
!!  lmnmax=max value of lmn over all psps
!!  lmn2max=max value of lmn2 over all psps
!!  natom=number of atoms in cell
!!  ndij=spin channels in dij
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dterm <type(dterm_type)> data related to onsite interactions
!!
!! SOURCE

subroutine dterm_alloc(dterm,lmnmax,lmn2max,natom,ndij)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,lmn2max,natom,ndij
  class(dterm_type),intent(inout) :: dterm
!--------------------------------------------------------------------

  dterm%lmnmax = lmnmax
  dterm%lmn2max = lmn2max
  dterm%natom = natom
  dterm%ndij = ndij

  ABI_REMALLOC(dterm%aij,(natom,lmn2max,ndij))
  dterm%has_aij=1

  ABI_REMALLOC(dterm%qij,(natom,lmn2max,ndij))
  dterm%has_qij=1

  ABI_REMALLOC(dterm%LR,(natom,lmn2max,ndij,3))
  dterm%has_LR=1

  ABI_REMALLOC(dterm%BM,(natom,lmn2max,ndij,3))
  dterm%has_BM=1

end subroutine dterm_alloc
!!***

!!****f* ABINIT/orbmag_mesh_alloc
!! NAME
!! orbmag_mesh_alloc
!!
!! FUNCTION
!! allocate space in orbmag_mesh_type
!!
!! INPUTS
!!  lmnmax=max value of lmn over all psps
!!  lmn2max=max value of lmn2 over all psps
!!  natom=number of atoms in cell
!!  ndij=spin channels in dij
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! orbmag_mesh <type(orbmag_mesh_type)> data related to orbmag terms on kpt mesh
!!
!! SOURCE

subroutine orbmag_mesh_alloc(dtset,orbmag_mesh)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),intent(in) :: dtset
  class(orbmag_mesh_type),intent(inout) :: orbmag_mesh
!--------------------------------------------------------------------

  orbmag_mesh%mband = dtset%mband
  orbmag_mesh%nkpt = dtset%nkpt
  orbmag_mesh%nsppol = dtset%nsppol
  orbmag_mesh%natom = dtset%natom
  orbmag_mesh%ntypat = dtset%ntypat
  orbmag_mesh%chern_nterms = chern_nterms
  orbmag_mesh%orbmag_nterms = orbmag_nterms

  ABI_REMALLOC(orbmag_mesh%lambsig,(orbmag_mesh%ntypat))
  ABI_REMALLOC(orbmag_mesh%nucdipmom,(3,orbmag_mesh%natom))
  ABI_REMALLOC(orbmag_mesh%cmesh,(orbmag_mesh%mband,orbmag_mesh%nkpt,orbmag_mesh%nsppol,3,chern_nterms))
  orbmag_mesh%cmesh=zero
  ABI_REMALLOC(orbmag_mesh%omesh,(orbmag_mesh%mband,orbmag_mesh%nkpt,orbmag_mesh%nsppol,3,orbmag_nterms))
  orbmag_mesh%omesh=zero

end subroutine orbmag_mesh_alloc
!!***

!!****f* ABINIT/orbmag_mesh_free
!! NAME
!! orbmag_mesh_free
!!
!! FUNCTION
!! free space in orbmag_mesh_type
!!
!! SIDE EFFECTS
!! orbmag_mesh <type(orbmag_mesh_type)> data related to orbmag terms on kpt mesh
!!
!! SOURCE

subroutine orbmag_mesh_free(orbmag_mesh)

  !Arguments ------------------------------------
  !scalars
  class(orbmag_mesh_type),intent(inout) :: orbmag_mesh
!--------------------------------------------------------------------

    ABI_SFREE(orbmag_mesh%lambsig)
    ABI_SFREE(orbmag_mesh%nucdipmom)
    ABI_SFREE(orbmag_mesh%cmesh)
    ABI_SFREE(orbmag_mesh%omesh)

end subroutine orbmag_mesh_free
!!***

!!****f* ABINIT/dterm_aij
!! NAME
!! dterm_aij
!!
!! FUNCTION
!! transfer paw_ij to dterm%aij in more convenient format
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  paw_ij(dtset%natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dterm <type(dterm_type)> data related to onsite interactions
!!
!! SOURCE

subroutine dterm_aij(atindx,dterm,dtset,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,idij,itypat,klmn
  !arrays
!--------------------------------------------------------------------

  dterm%aij = czero

  ! note that paw_ij has atoms ordered by input, while
  ! we would like to order by groups of atoms with atindx
  do iat=1,dtset%natom
    iatom=atindx(iat)
    itypat=dtset%typat(iat)
    do klmn=1,pawtab(itypat)%lmn2_size
      do idij = 1, dterm%ndij
        if (paw_ij(iat)%cplex_dij .EQ. 2) then
          dterm%aij(iatom,klmn,idij) = &
            & CMPLX(paw_ij(iat)%dij(2*klmn-1,idij),paw_ij(iat)%dij(2*klmn,idij))
        else
          dterm%aij(iatom,klmn,idij) = CMPLX(paw_ij(iat)%dij(klmn,idij),zero)
        end if
      end do ! idij
    end do ! klmn
  end do ! iat

  dterm%has_aij = 2

end subroutine dterm_aij

!!***

!!****f* ABINIT/make_d
!! NAME
!! make_d
!!
!! FUNCTION
!! this is a driver to compute different onsite terms, in convenient (complex) format
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=reciprocal space lattice vectors
!!  paw_ij(dtset%natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawrad(dtset%ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dterm <type(dterm_type)> data related to onsite interactions
!!
!! SOURCE

subroutine make_d(atindx,dterm,dtset,gprimd,paw_ij,pawrad,pawtab,psps)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type), intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: my_lmax,ngnt

  !arrays
  integer,allocatable :: gntselect(:,:)
  real(dp),allocatable :: realgnt(:)
!--------------------------------------------------------------------

 ! make Gaunt integrals
 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 ! generate d terms

 ! normal PAW sij overlap, in complex form because it's convenient
 call dterm_qij(atindx,dterm,dtset,pawtab)

 ! onsite angular momentum expectation values
 call dterm_LR(atindx,dterm,dtset,gprimd,pawrad,pawtab)

 ! onsite <A_0.A_n> interaction between magnetic field and nuclear dipole
 !call dterm_BM(atindx,dterm,dtset,gntselect,gprimd,my_lmax,pawrad,pawtab,realgnt)
 call dterm_BM(atindx,dterm,dtset,gntselect,gprimd,my_lmax,pawrad,pawtab,realgnt)

 ! transfers paw_ij to dterm%aij because it's convenient
 call dterm_aij(atindx,dterm,dtset,paw_ij,pawtab)

 ABI_FREE(realgnt)
 ABI_FREE(gntselect)

end subroutine make_d
!!***

!!****f* ABINIT/local_fermie
!! NAME
!! local_fermie
!!
!! FUNCTION
!! estimate Fermi energy as max value of all occupied bands/kpts
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)=ground state eigenvalues at each band and kpt
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  occ(dtset%mband*dtset%nkpt*dtset%nsppol)=occup number for each band (often 2) at each k point
!!
!! OUTPUT
!!  fermie=maximum energy (real(dp)) found over all occupied input bands
!!
!! CHILDREN
!!
!! SOURCE

subroutine local_fermie(dtset,ebands_k,fermie,mpi_enreg)

  !Arguments ------------------------------------
  !scalars
  real(dp),intent(out) :: fermie
  type(dataset_type),intent(in) :: dtset
  type(ebands_t) :: ebands_k
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays

  !Local variables -------------------------
  !scalars
  integer :: ierr,ikpt,isppol,me
  integer :: nband_k,nn,nproc,spaceComm
  real(dp) :: fermie_proc

  !arrays
  real(dp),allocatable :: eig_k(:),occ_k(:)

!--------------------------------------------------------------------

  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  fermie_proc = -1.0D99
  do isppol = 1, dtset%nsppol
    do ikpt = 1, dtset%nkpt

      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)

      ! if the current kpt is not on the current processor, cycle
      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

      ABI_MALLOC(occ_k,(nband_k))
      occ_k(:)=ebands_k%occ(1:nband_k,ikpt,isppol)

      ABI_MALLOC(eig_k,(nband_k))
      eig_k(:)=ebands_k%eig(1:nband_k,ikpt,isppol)

      do nn = 1, nband_k
        if ( (abs(occ_k(nn)).GT.tol8) .AND. (eig_k(nn).GT.fermie_proc) ) then
          fermie_proc = eig_k(nn)
        end if
      end do ! nn

      ABI_FREE(occ_k)
      ABI_FREE(eig_k)

    end do ! end loop over kpts
  end do ! end loop over isppol

  call xmpi_max(fermie_proc,fermie,spaceComm,ierr)

end subroutine local_fermie
!!***

!!****f* m_orbmag/orbmag_ncwrite
!! NAME
!! orbmag_ncwrite
!!
!! FUNCTION
!!  Write orbmag_mesh contributions to netcdf file.
!!
!! INPUTS
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  dtset<dtset_type>=Dataset type
!!  ebands<ebands_t>=Band structure data.
!!  hdr<hdr_t>=Abinit header
!!  orbmag_mesh<orbmag_mesh_type>=orbmag_mesh data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! SOURCE

subroutine orbmag_ncwrite(crystal,dtset,ebands,hdr,ncid,orbmag_mesh)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(crystal_t),intent(in) :: crystal
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(hdr_type),intent(in) :: hdr
 type(orbmag_mesh_type),intent(in) :: orbmag_mesh
!arrays

!Local variables-------------------------------
!scalars
 integer :: ncerr,fform
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
!arrays
!*************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 fform = fform_from_ext("ORBMAG.nc")
 ABI_CHECK(fform /= 0, "Cannot find fform associated to ORBMAG.nc")

 ! Write header, crystal structure and band energies.
 NCF_CHECK(hdr%ncwrite(ncid, fform, nc_define=.True.))
 NCF_CHECK(crystal%ncwrite(ncid))
 NCF_CHECK(ebands%ncwrite(ncid))

 !! Add orbmag-mesh-specific quantities
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("ntypat", dtset%ntypat),&
   nctkdim_t("mband", orbmag_mesh%mband),&
   nctkdim_t("nkpt", orbmag_mesh%nkpt),&
   nctkdim_t("nsppol", orbmag_mesh%nsppol),&
   nctkdim_t("chern_nterms", orbmag_mesh%chern_nterms),&
   nctkdim_t("orbmag_nterms", orbmag_mesh%orbmag_nterms),&
   nctkdim_t("ndir",3),&
   nctkdim_t("natom",dtset%natom)],defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("chern_mesh", "dp", "mband, nkpt, nsppol, ndir, chern_nterms"),&
   nctkarr_t("orbmag_mesh", "dp", "mband, nkpt, nsppol, ndir, orbmag_nterms"),&
   nctkarr_t("lambsig", "dp", "ntypat"),&
   nctkarr_t("nucdipmom", "dp", "ndir, natom")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))

 NCF_CHECK(nf90_put_var(ncid, vid("chern_mesh"), orbmag_mesh%cmesh))
 NCF_CHECK(nf90_put_var(ncid, vid("orbmag_mesh"), orbmag_mesh%omesh))
 NCF_CHECK(nf90_put_var(ncid, vid("nucdipmom"), orbmag_mesh%nucdipmom))
 NCF_CHECK(nf90_put_var(ncid, vid("lambsig"), orbmag_mesh%lambsig))

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" orbmag_ncwrite: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine orbmag_ncwrite
!!***

end module m_orbmag
