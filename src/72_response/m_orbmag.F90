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
!! Phys Rev B 107, 165157 (2023). This paper will be referred to in the comments as ZTG23.
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
  use m_cgtools,          only : cg_zdotc,cg_zdotu,projbd
  use m_dtfil
  use m_ebands
  use m_fft,              only : fourwf
  use m_getghc,           only : getghc
  use m_getgh1c
  use m_hamiltonian
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

    integer :: n4,n5,n6
    ! real space grid dimenions for rmesh

    real(dp),allocatable :: lambsig(:)
    ! lambsig(ntypat)

    real(dp),allocatable :: nucdipmom(:,:)
    ! nucdipmom(3,natom)

    real(dp),allocatable :: cmesh(:,:,:,:,:)
    ! 3 for the 3 directions
    ! cmesh(mband,nkpt,nsppol,3,chern_terms)
    
    real(dp),allocatable :: chern_terms(:,:,:,:)
    ! 3 for the 3 directions
    ! chern_terms(dtset%mband,dtset%nsppol,3,chern_nterms)

    real(dp),allocatable :: chern_trace(:,:)
    ! 3 for the 3 directions
    ! chern_trace(3,chern_nterms)

    real(dp),allocatable :: omesh(:,:,:,:,:)
    ! 3 for the 3 directions
    ! omesh(mband,nkpt,nsppol,3,orbmag_terms)
    
    real(dp),allocatable :: orbmag_terms(:,:,:,:)
    ! 3 for the 3 directions
    ! orbmag_terms(dtset%mband,dtset%nsppol,3,orbmag_nterms)
    
    real(dp),allocatable :: orbmag_trace(:,:)
    ! 3 for the 3 directions
    ! orbmag_trace(3,orbmag_nterms)
    
    real(dp),allocatable :: rmesh(:,:,:,:,:)
    ! total orbmag on real mesh
    ! 3 for the 3 directions
    ! rmesh(n4,n5,n6,3,orbmag_nterms)

    contains

      procedure :: init => orbmag_init
      procedure :: free => orbmag_free
      procedure :: accum_rmesh => orbmag_rmesh
      procedure :: mpisum => orbmag_mpisum
      procedure :: term_scale => orbmag_term_scale
      procedure :: output => orbmag_output

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

    contains

      procedure :: init => dterm_init
      procedure :: free => dterm_free

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
  private :: gauge_treatment
  private :: para_to_diag
  private :: orbmag_init
  private :: orbmag_free
  private :: orbmag_mpisum
  private :: orbmag_rmesh
  private :: orbmag_term_scale
  private :: orbmag_output
  private :: orbmag_ncwrite   ! Write orbmag_mesh contributions to netcdf file.
  private :: dterm_init
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
 integer :: me,mcgk,mcprjk,my_nspinor,nband_k,nband_me,ncid,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nl1_option,nn,nkpg,npw_k,npwsp,nproc,nucdip_dirs,spaceComm
 integer,parameter :: master=0
 real(dp) :: arg,ecut_eff,fermie
 logical :: has_nucdip
 type(dterm_type) :: dterm
 type(gs_hamiltonian_type) :: gs_hamk
 type(orbmag_mesh_type) :: orbmag_mesh

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),kg_k(:,:),nattyp(:)
 real(dp) :: kpoint(3),omlamb(3)
 real(dp),allocatable :: cg1_k(:,:,:),cwavef(:,:),dkinpw(:,:),eig_k(:)
 real(dp),allocatable :: ffnl_k(:,:,:,:),kinpw(:),kpg_k(:,:),occ_k(:)
 real(dp),allocatable,target :: cg_k(:,:),gcg1_k(:,:,:)
 real(dp),allocatable,target :: ph1d(:,:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),trnrm(:)
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
 if (abs(dtset%orbmag) .EQ. 3) then
   ABI_MALLOC(dkinpw,(mpw,3))
 end if

 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,0,dimlmn)

 lmn2max = psps%lmnmax*(psps%lmnmax+1)/2
 ! note: in make_d, terms will be filled as iatom using atindx
 call dterm%init(psps%lmnmax,lmn2max,dtset%natom,paw_ij(1)%ndij)
 call make_d(atindx,dterm,dtset,crystal%gprimd,paw_ij,pawrad,pawtab,psps)

 ! initialize orbmag_mesh datatype
 call orbmag_mesh%init(dtset)
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
     nucdip_dirs=3
     ABI_MALLOC(vectornd,(nfftf,dtset%nspden,nucdip_dirs))
     vectornd = zero
     call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,&
       & dtset%nspden,dtset%nucdipmom,crystal%rprimd,vectornd,crystal%xred)
     ABI_MALLOC(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,nucdip_dirs))
     call gspot_transgrid_and_pack(isppol, psps%usepaw, dtset%paral_kgb, dtset%nfft, dtset%ngfft, nfftf, &
          & dtset%nspden, gs_hamk%nvloc, nucdip_dirs, pawfgr, mpi_enreg, vectornd,vectornd_pac)
     ABI_FREE(vectornd)
     call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
   else
     nucdip_dirs=0
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
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,crystal%gmet,&
       & kg_k,kinpw,kpoint,npw_k,0,0)
     if (abs(dtset%orbmag).EQ.3) then
       do adir=1,3
         call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,crystal%gmet,&
           & kg_k,dkinpw(:,adir),kpoint,npw_k,adir,0)
       end do
     end if

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

     ABI_SFREE(ylm_k)
     ABI_SFREE(ylmgr_k)
     
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

     ! gauge treatment of cg1_k
     ABI_MALLOC(gcg1_k,(2,mcgk,3))
     call gauge_treatment(atindx,cg_k,cg1_k,cprj_k,dimlmn,dkinpw,dtset,eig_k,gcg1_k,gs_hamk,&
         & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,mpw,nband_k,ngfft4,ngfft5,ngfft6,npw_k,&
         & nucdip_dirs,occ_k,vectornd_pac)
     
     ! compute <p|gcg1> cprjs
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
         cwavef(1:2,1:npwsp) = gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,adir)
         call getcprj(choice,cpopt,cwavef,cwaveprj,gs_hamk%ffnl_k,idir,&
           & psps%indlmn,istwf_k,kg_k,gs_hamk%kpg_k,kpoint,psps%lmnmax,dtset%mgfft,&
           & mpi_enreg,1,dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,&
           & dtset%nspinor,dtset%ntypat,phkxred,ph1d,gs_hamk%ph3d_k,&
           & crystal%ucvol,psps%useylm)
         call pawcprj_put(atindx,cwaveprj,cprj1_k(:,:,adir),dtset%natom,&
           & nn,0,ikpt,0,isppol,dtset%mband,mkmem_rbz,dtset%natom,1,nband_k,&
           & dimlmn,dtset%nspinor,dtset%nsppol,0)
       end do
     end do
     ABI_SFREE(cwavef)
     ABI_SFREE(phkxred)

     ! set up normalization factors at this k point
     ABI_MALLOC(trnrm,(nband_k))
     trnrm(1:nband_k) = ebands_k%occ(1:nband_k,ikpt,isppol)*dtset%wtk(ikpt)/crystal%ucvol

     !--------------------------------------------------------------------------------
     ! Finally ready to compute contributions to orbital magnetism and Berry curvature
     !--------------------------------------------------------------------------------

     ! ZTG23 Eq. 36 term 2 and Eq. 46 term 1
     call orbmag_cc_k(atindx,cprj1_k,dimlmn,dtset,eig_k,fermie,gcg1_k,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,ph1d,trnrm,suppress_ormesh=.FALSE.)

     ! ZTG23 Eq. 36 terms 3 and 4 and Eq. 46 term 2
     call orbmag_vv_k(atindx,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gcg1_k,gs_hamk,&
      & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,orbmag_mesh,&
      & ph1d,trnrm,suppress_ormesh=.TRUE.)

     ! ZTG23 Eq. 36 term 1
     call orbmag_nl_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,pawtab,ph1d,trnrm,&
       & suppress_ormesh=.FALSE.)

     ! ZTG23 text after Eq. 42
     nl1_option = 1 ! LR
     call orbmag_nl1_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,nl1_option,npw_k,orbmag_mesh,pawtab,&
       & ph1d,trnrm,suppress_ormesh=.FALSE.)

     ! ZTG23 Eq. 43
     nl1_option = 2 ! A0.An
     call orbmag_nl1_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,nl1_option,npw_k,orbmag_mesh,pawtab,&
       & ph1d,trnrm,suppress_ormesh=.FALSE.)

     ! accumulate terms
     do nn = 1, nband_k
       if(abs(trnrm(nn)).LT.tol8) cycle
       orbmag_mesh%chern_terms(nn,isppol,1:3,ibcc:ibvv2) = orbmag_mesh%chern_terms(nn,isppol,1:3,ibcc:ibvv2) + &
           & trnrm(nn)*orbmag_mesh%cmesh(nn,ikpt,isppol,1:3,ibcc:ibvv2)
       orbmag_mesh%orbmag_terms(nn,isppol,1:3,incc:inbm) = orbmag_mesh%orbmag_terms(nn,isppol,1:3,incc:inbm) + &
           & trnrm(nn)*orbmag_mesh%omesh(nn,ikpt,isppol,1:3,incc:inbm)
     end do ! loop on bands

     icg = icg + mcgk
     icprj = icprj + mcprjk
     ikg = ikg + npw_k
     bdtot_index=bdtot_index+nband_k

     ABI_SFREE(ffnl_k)
     ABI_SFREE(ph3d)
     ABI_SFREE(kpg_k)
     ABI_SFREE(cg_k)
     ABI_SFREE(cg1_k)
     ABI_SFREE(gcg1_k)
     ABI_SFREE(eig_k)
     ABI_SFREE(occ_k)
     call pawcprj_free(cprj_k)
     ABI_SFREE(cprj_k)
     do adir = 1, 3
       call pawcprj_free(cprj1_k(:,:,adir))
     end do
     ABI_SFREE(cprj1_k)
     ABI_SFREE(trnrm)

   end do ! end loop over kpts

   ABI_SFREE(vlocal)
   ABI_SFREE(vectornd_pac)
   ABI_SFREE(vxctaulocal)

 end do ! end loop over isppol

 ! accumulate data over processors
 call orbmag_mesh%mpisum(nproc,spaceComm)
 
 ! prepare terms for output to abo file
 call orbmag_mesh%term_scale(crystal,dtset)

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
 call orbmag_mesh%output(dtset,omlamb)

!---------------------------------------------------
! deallocate memory
!---------------------------------------------------

 call gs_hamk%free()

 ABI_SFREE(kg_k)
 ABI_SFREE(kinpw)
 ABI_SFREE(dkinpw)
 ABI_SFREE(ph1d)

 ABI_SFREE(atindx)
 ABI_SFREE(atindx1)
 ABI_SFREE(nattyp)

 ABI_FREE(dimlmn)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 call dterm%free()
 call orbmag_mesh%free()

end subroutine orbmag
!!***

!!****f*m_orbmag/orbmag_mpisum
!! NAME
!! orbmag_mpisum
!!
!! FUNCTION
!! accumulate data in orbmag_mesh_type over processes
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine orbmag_mpisum(omag,nproc,spaceComm)
  !Arguments
  class(orbmag_mesh_type),intent(inout),target :: omag
  integer,intent(in) :: nproc,spaceComm

  !Local variables
  integer :: buff_size,ierr
  real(dp),allocatable :: buffer1(:),buffer2(:)

  if (nproc > 1) then
    if (allocated(omag%cmesh)) then
      buff_size=size(omag%cmesh)
      ABI_MALLOC(buffer1,(buff_size))
      ABI_MALLOC(buffer2,(buff_size))
      buffer1=zero;buffer2=zero
      buffer1(1:buff_size) = &
        & reshape(omag%cmesh,(/omag%mband*omag%nkpt*omag%nsppol*3*chern_nterms/))
      call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
      omag%cmesh(1:omag%mband,1:omag%nkpt,1:omag%nsppol,1:3,1:chern_nterms)=&
        & reshape(buffer2,(/omag%mband,omag%nkpt,omag%nsppol,3,chern_nterms/))
      ABI_FREE(buffer1)
      ABI_FREE(buffer2)
    end if
    if (allocated(omag%chern_terms)) then
      buff_size=size(omag%chern_terms)
      ABI_MALLOC(buffer1,(buff_size))
      ABI_MALLOC(buffer2,(buff_size))
      buffer1=zero;buffer2=zero
      buffer1(1:buff_size) = &
        & reshape(omag%chern_terms,(/omag%mband*omag%nsppol*3*chern_nterms/))
      call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
      omag%chern_terms(1:omag%mband,1:omag%nsppol,1:3,1:chern_nterms)=&
        & reshape(buffer2,(/omag%mband,omag%nsppol,3,chern_nterms/))
      ABI_FREE(buffer1)
      ABI_FREE(buffer2)
    end if
    if (allocated(omag%omesh)) then
      buff_size=size(omag%omesh)
      ABI_MALLOC(buffer1,(buff_size))
      ABI_MALLOC(buffer2,(buff_size))
      buffer1=zero;buffer2=zero
      buffer1(1:buff_size) = &
        & reshape(omag%omesh,(/omag%mband*omag%nkpt*omag%nsppol*3*orbmag_nterms/))
      call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
      omag%omesh(1:omag%mband,1:omag%nkpt,1:omag%nsppol,1:3,1:orbmag_nterms)=&
        & reshape(buffer2,(/omag%mband,omag%nkpt,omag%nsppol,3,orbmag_nterms/))
      ABI_FREE(buffer1)
      ABI_FREE(buffer2)
    end if
    if (allocated(omag%orbmag_terms)) then
      buff_size=size(omag%orbmag_terms)
      ABI_MALLOC(buffer1,(buff_size))
      ABI_MALLOC(buffer2,(buff_size))
      buffer1=zero;buffer2=zero
      buffer1(1:buff_size) = &
        & reshape(omag%orbmag_terms,(/omag%mband*omag%nsppol*3*orbmag_nterms/))
      call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
      omag%orbmag_terms(1:omag%mband,1:omag%nsppol,1:3,1:orbmag_nterms)=&
        & reshape(buffer2,(/omag%mband,omag%nsppol,3,orbmag_nterms/))
      ABI_FREE(buffer1)
      ABI_FREE(buffer2)
    end if
    if (allocated(omag%rmesh)) then
      buff_size=size(omag%rmesh)
      ABI_MALLOC(buffer1,(buff_size))
      ABI_MALLOC(buffer2,(buff_size))
      buffer1=zero;buffer2=zero
      buffer1(1:buff_size) = &
        & reshape(omag%rmesh,(/omag%n4*omag%n5*omag%n6*3*orbmag_nterms/))
      call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
      omag%rmesh(1:omag%n4,1:omag%n5,1:omag%n6,1:3,1:orbmag_nterms)=&
        & reshape(buffer2,(/omag%n4,omag%n5,omag%n6,3,orbmag_nterms/))
      ABI_FREE(buffer1)
      ABI_FREE(buffer2)
    end if
 
  end if

end subroutine orbmag_mpisum
!!***
  

!!****f*m_orbmag/orbmag_term_scale
!! NAME
!! orbmag_term_scale
!!
!! FUNCTION
!! change frames and scale terms as needed
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine orbmag_term_scale(omag,crystal,dtset)

  !Arguments ------------------------------------
  !scalars
  class(orbmag_mesh_type),intent(inout),target :: omag
  type(crystal_t),intent(in) :: crystal
  type(dataset_type),intent(in) :: dtset

  !arrays

  !Local variables -------------------------
  !scalars
  integer :: i4,i5,i6,isppol,iterm,nn
  !arrays

!--------------------------------------------------------------------

 do iterm = 1, orbmag_nterms
   do isppol = 1, dtset%nsppol
     do nn = 1, omag%mband
       if((iterm.EQ.inlr).OR.(iterm.EQ.inbm)) then
         omag%orbmag_terms(nn,isppol,1:3,iterm) = &
           & MATMUL(crystal%rprimd,omag%orbmag_terms(nn,isppol,1:3,iterm))
       else
         omag%orbmag_terms(nn,isppol,1:3,iterm) = &
           & crystal%ucvol*MATMUL(crystal%gprimd,omag%orbmag_terms(nn,isppol,1:3,iterm))
       end if
     end do ! nn
   end do !isppol
 end do

 if (dtset%orbmag .EQ. 4) then
   do iterm = 1, orbmag_nterms
     do i4=1,omag%n4
       do i5=1,omag%n5
         do i6=1,omag%n6
           if((iterm.EQ.inlr).OR.(iterm.EQ.inbm)) then
             omag%rmesh(i4,i5,i6,1:3,iterm) = &
               & MATMUL(crystal%rprimd,omag%rmesh(i4,i5,i6,1:3,iterm))
           else
             omag%rmesh(i4,i5,i6,1:3,iterm) = &
               & crystal%ucvol*MATMUL(crystal%gprimd,omag%rmesh(i4,i5,i6,1:3,iterm))
           end if
         end do
       end do
     end do
   end do
 end if


 do iterm = 1, chern_nterms
   do isppol = 1, dtset%nsppol
     do nn = 1, omag%mband
       omag%chern_terms(nn,isppol,1:3,iterm) = &
         & crystal%ucvol*MATMUL(crystal%gprimd,omag%chern_terms(nn,isppol,1:3,iterm))
     end do ! nn
   end do !isppol
 end do

 !! convert orbmag magnetization to orbital moment
 !! Berry curvature terms are ignored
 omag%orbmag_terms(:,:,:,incc:inbm)=crystal%ucvol*omag%orbmag_terms(:,:,:,incc:inbm)

 !! accumulate trace of terms 
 do isppol = 1, dtset%nsppol
   do nn = 1, omag%mband
     omag%orbmag_trace(1:3,1:orbmag_nterms) = omag%orbmag_trace(1:3,1:orbmag_nterms) + &
       & omag%orbmag_terms(nn,isppol,1:3,1:orbmag_nterms)
     omag%chern_trace(1:3,1:chern_nterms) = omag%chern_trace(1:3,1:chern_nterms) + &
       & omag%chern_terms(nn,isppol,1:3,1:chern_nterms)
   end do ! nn
 end do ! isppol

end subroutine orbmag_term_scale
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
!!  cg_k(2,mcgk) ground state wavefunctions at this k point
!!  cprj_k(dtset%natom,mcprjk)<type(pawcprj_type)>=cprj for cg_k
!!  dimlmn(dtset%natom)=cprj lmn dimensions
!!  dterm <type(dterm_type)> data related to onsite interactions
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcgk=2nd dimension of cg_k
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  nband_k=bands at this kpt
!!  nl1_option=chooses which onsite term to apply
!!  npw_k=planewaves at this k point
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

subroutine orbmag_nl1_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,gs_hamk,ikpt,isppol,&
    & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,nl1_option,npw_k,orbmag_mesh,pawtab,ph1d,trnrm,&
    & suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,nl1_option,npw_k
  logical,intent(in),optional :: suppress_ormesh
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  !real(dp),intent(in) :: ph1d(2,dtset%natom*(2*(dtset%ngfft(1)+dtset%ngfft(2)+dtset%ngfft(3))))
  real(dp),intent(in),pointer,dimension(:,:) :: ph1d
  real(dp),intent(in) :: trnrm(nband_k)
  real(dp),intent(in),target :: cg_k(2,mcgk)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,nn,npwsp
  complex(dp) :: tt
  logical :: my_suppress_ormesh
  !arrays
  real(dp),pointer :: cwavef(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)
!--------------------------------------------------------------------

 if(present(suppress_ormesh)) then
   my_suppress_ormesh=suppress_ormesh
 else
   my_suppress_ormesh=.FALSE.
 end if
 npwsp = npw_k*dtset%nspinor
 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 !ABI_MALLOC(cwavef,(2,npwsp))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do nn = 1, nband_k
   cwavef => cg_k(1:2,(nn-1)*npwsp+1:nn*npwsp)
   call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
     & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)
  
   do adir = 1, 3

     select case (nl1_option)
     case(1)
       call tt_me(adir,dterm%LR(:,:,:,adir),atindx,cwavef,dtset,gs_hamk,dterm%lmn2max,mpi_enreg,&
         & dterm%ndij,nband_k,npw_k,orbmag_mesh,inlr,pawtab,ph1d,tt,trnrm(nn),cwaveprj,&
         & suppress_ormesh=my_suppress_ormesh)
       orbmag_mesh%omesh(nn,ikpt,isppol,adir,inlr) = real(tt)
     case(2)
       call tt_me(adir,dterm%BM(:,:,:,adir),atindx,cwavef,dtset,gs_hamk,dterm%lmn2max,mpi_enreg,&
         & dterm%ndij,nband_k,npw_k,orbmag_mesh,inbm,pawtab,ph1d,tt,trnrm(nn),cwaveprj,&
         & suppress_ormesh=my_suppress_ormesh)
       orbmag_mesh%omesh(nn,ikpt,isppol,adir,inbm) = real(tt)
     case default
       tt = czero
     end select

   end do !adir
 
 end do !nn

 call pawcprj_free(cwaveprj)
 ABI_SFREE(cwaveprj)
 IF(ASSOCIATED(cwavef)) NULLIFY(cwavef)

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

subroutine orbmag_nl_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,gs_hamk,ikpt,isppol,&
    & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,pawtab,ph1d,trnrm,&
    & suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  logical,intent(in),optional :: suppress_ormesh
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in),target :: cg_k(2,mcgk)
  real(dp),intent(in) :: eig_k(nband_k),trnrm(nband_k)
  real(dp),intent(in),pointer,dimension(:,:) :: ph1d
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,nn,npwsp
  real(dp) :: epsabg
  complex(dp) :: prefac_m,txt
  logical :: my_suppress_ormesh
  !arrays
  real(dp),pointer :: cwavef(:,:)
  complex(dp) :: m1(3)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 if(present(suppress_ormesh)) then
   my_suppress_ormesh=suppress_ormesh
 else
   my_suppress_ormesh=.FALSE.
 end if
 
 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do nn = 1, nband_k

   cwavef => cg_k(1:2,(nn-1)*npwsp+1:nn*npwsp)
   call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
     & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

   m1(1:3) = czero
   do bdir = 1, 3
     do gdir = 1, 3
       do adir = 1, 3
         
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         ! note rho^0 H^1 terms have opposite sign of rho^1 H^0
         prefac_m = -com*c2*epsabg

         call txt_me(adir,atindx,bdir,cwavef,dterm,dtset,eig_k,gdir,gs_hamk,&
           & nn,mpi_enreg,nband_k,npw_k,orbmag_mesh,innl,pawtab,ph1d,prefac_m,txt,&
           & trnrm(nn),cwaveprj,&
           & suppress_ormesh=my_suppress_ormesh)
       
         m1(adir) = m1(adir) + txt

       end do ! adir
     end do !gdir
   end do !bdir

   orbmag_mesh%omesh(nn,ikpt,isppol,1:3,innl) = real(m1(1:3))

 end do !nn

 call pawcprj_free(cwaveprj)
 ABI_SFREE(cwaveprj)
 IF(ASSOCIATED(cwavef)) NULLIFY(cwavef)

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
!!  gcg1_k(2,mcgk,3)=gauge adjusted cg1_k
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

subroutine orbmag_cc_k(atindx,cprj1_k,dimlmn,dtset,eig_k,fermie,gcg1_k,gs_hamk,ikpt,isppol,&
    & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,orbmag_mesh,ph1d,trnrm,&
    & suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie
  logical,intent(in),optional :: suppress_ormesh
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),trnrm(nband_k)
  real(dp),intent(in),target :: gcg1_k(2,mcgk,3)
  real(dp),intent(in),pointer,dimension(:,:) :: ph1d
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,mcprjk,3)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,fourwf_cplex,fourwf_option,gdir,iatom,ipw,ndat
  integer :: nn,npwsp,sij_opt,t_atom,tim_fourwf,tim_getghc,type_calc
  integer :: kg1,kg2,kg3,n1,n2,n3,shift1,shift2,shift3
  real(dp) :: epsabg,lams,weight_i,weight_r
  complex(dp) :: bc,kc,ormesh_fac,ph1,ph2,ph3,prefac_b,prefac_m 
  logical :: my_suppress_ormesh,need_ormesh
  !arrays
  real(dp) bdot(2),mdot(2)
  real(dp),allocatable :: denpot(:,:,:),fofgout(:,:)
  real(dp),allocatable :: ghc(:,:),gsc(:,:),gvnlxc(:,:)
  real(dp),allocatable,target :: ghc_local(:,:)
  real(dp),pointer :: bra(:,:),ket(:,:)
  complex(dp) :: m1(3),b1(3)
  type(pawcprj_type),allocatable :: cwaveprj1(:,:)
!--------------------------------------------------------------------

 if(present(suppress_ormesh)) then
   my_suppress_ormesh=suppress_ormesh
 else
   my_suppress_ormesh=.FALSE.
 end if
 npwsp = npw_k*dtset%nspinor
 need_ormesh = ((dtset%orbmag .EQ. 4) .AND. (.NOT. my_suppress_ormesh))

 ABI_MALLOC(ghc,(2,npwsp))
 ABI_MALLOC(gsc,(2,npwsp))
 ABI_MALLOC(gvnlxc,(2,npwsp))
 ABI_MALLOC(cwaveprj1,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj1,0,dimlmn)

 tim_getghc = 0
 lams = zero
 ndat = 1
 
 if (need_ormesh) then
   ABI_MALLOC(ghc_local,(2,npwsp))
   ! need atom index with dipole for ph3d use below
   do iatom = 1, dtset%natom
     if ( ANY(ABS(dtset%nucdipmom(1:3,iatom))>tol8) ) then
       t_atom = atindx(iatom)
       exit
     end if
   end do
 end if

 do nn = 1, nband_k

   m1 = czero
   b1 = czero
     
   do gdir = 1, 3

     cpopt = 2
     ket => gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,gdir)

     call pawcprj_get(atindx,cwaveprj1,cprj1_k(:,:,gdir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

     ! compute H|Pc d_\gamma u> and S|Pc d_\gamma u>
     ! [H+E*S - 2\mu*S]|ket> is needed for orbmag
     ! -2*S|ket> needed for Chern
     type_calc = 0 ! apply local and non-local Hamiltonian
     sij_opt = 1 ! compute gsc in addition to ghc
     call getghc(cpopt,ket,cwaveprj1,ghc,gsc,gs_hamk,gvnlxc,lams,mpi_enreg,&
       & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)

     ghc(1:2,1:npwsp) = ghc(1:2,1:npwsp) + gsc(1:2,1:npwsp)*(eig_k(nn) - two*fermie)
     
     if (need_ormesh) then
       type_calc = 3 ! apply local and kinetic only
       sij_opt = 0 ! compute ghc only
       call getghc(cpopt,ket,cwaveprj1,ghc_local,gsc,gs_hamk,gvnlxc,lams,mpi_enreg,&
         & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)
     end if

     do bdir = 1, 3
       !bra(1:2,1:npwsp) = gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,bdir)
       bra => gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,bdir)

       mdot = cg_zdotc(npwsp,bra,ghc); bdot = cg_zdotc(npwsp,bra,gsc)

       ! assemble contributions alpha_dir \propto \beta_dir x \gamma_dir
       do adir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg
         m1(adir) = m1(adir) + prefac_m*CMPLX(mdot(1),mdot(2))
         b1(adir) = b1(adir) - two*prefac_b*CMPLX(bdot(1),bdot(2))
         
         if (need_ormesh) then
           ormesh_fac = trnrm(nn)*prefac_m
           call orbmag_mesh%accum_rmesh(adir,bra,dtset,gs_hamk,ghc_local,.TRUE.,mpi_enreg,&
             & npwsp,ph1d,ormesh_fac,t_atom,incc)
         end if
       
       end do ! adir
   
     end do !bdir
   end do !gdir

   orbmag_mesh%omesh(nn,ikpt,isppol,1:3,incc) = real(m1(1:3))
   orbmag_mesh%cmesh(nn,ikpt,isppol,1:3,ibcc) = real(b1(1:3))

 end do !nn

 if(ASSOCIATED(ket)) NULLIFY(ket)
 if(ASSOCIATED(bra)) NULLIFY(bra)

 ABI_SFREE(ghc)
 ABI_SFREE(gsc)
 ABI_SFREE(gvnlxc)
 call pawcprj_free(cwaveprj1)
 ABI_SFREE(cwaveprj1)
 ABI_SFREE(ghc_local)

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
!!  gcg1_k(2,mcgk,3)=gauge treated cg1_k
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  ikpt=current k pt
!!  isppol=current spin polarization
!!  mcgk=dimension of cg_k
!!  mcprjk=dimension of cprj_k
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  nband_k=bands at this kpt
!!  npw_k=number of planewaves at this kpt
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

subroutine orbmag_vv_k(atindx,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gcg1_k,gs_hamk,&
    & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,orbmag_mesh,&
    & ph1d,trnrm,suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie
  logical,intent(in),optional :: suppress_ormesh
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),occ_k(nband_k),trnrm(nband_k)
  real(dp),intent(in),target :: cg_k(2,mcgk),gcg1_k(2,mcgk,3)
  real(dp),intent(in),pointer,dimension(:,:) :: ph1d
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,fourwf_cplex,fourwf_option,gdir,iatom,ndat,nn,nnlout,np,npwsp
  integer :: paw_opt,signs,t_atom,tim_fourwf,tim_getghc
  real(dp) :: epsabg,weight_i,weight_r
  complex(dp) :: bdotc,bpdotc,gdotc,gpdotc,ormesh_fac,prefac_b,prefac_m
  logical :: my_suppress_ormesh,need_ormesh
  !arrays
  real(dp) :: bdot(2),bpdot(2),gdot(2),gpdot(2),enlout(1),lamv(1)
  real(dp),allocatable :: denpot(:,:,:)
  real(dp),allocatable :: fofgin(:,:),fofgout(:,:),vectout(:,:)
  real(dp),allocatable,target :: fofrb(:,:,:,:),fofrg(:,:,:,:),svectoutb(:,:),svectoutg(:,:)
  real(dp),pointer :: bra(:,:),brab(:,:),brag(:,:),ket(:,:)
  complex(dp) :: b1(3),bv2b(3),m1(3),mv2b(3),m1_mu(3),mv2b_mu(3)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)
!--------------------------------------------------------------------

 if(present(suppress_ormesh)) then
   my_suppress_ormesh=suppress_ormesh
 else
   my_suppress_ormesh=.FALSE.
 end if
 fourwf_cplex = 1
 fourwf_option = 0
 tim_fourwf = 1
 npwsp = npw_k*dtset%nspinor
 need_ormesh = ((dtset%orbmag .EQ. 4) .AND. (.NOT. my_suppress_ormesh))

 ABI_MALLOC(svectoutb,(2,npwsp))
 ABI_MALLOC(svectoutg,(2,npwsp))
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

 if (need_ormesh) then
   ABI_MALLOC(fofrb,(2,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6*ndat))
   ABI_MALLOC(fofrg,(2,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6*ndat))
   ! need atom index with dipole for ph3d use below
   do iatom = 1, dtset%natom
     if ( ANY(ABS(dtset%nucdipmom(1:3,iatom))>tol8) ) then
       t_atom = atindx(iatom)
       exit
     end if
   end do
 end if

 do nn = 1, nband_k
     
   m1 = czero; mv2b = czero
   m1_mu = czero; mv2b_mu = czero
   b1 = czero; bv2b = czero

   ! extract |u_nk>
   ket => cg_k(1:2,(nn-1)*npwsp+1:nn*npwsp)
   call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
     & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

   do bdir = 1, 3

     ! compute dS/dk_b|u_nk>
     call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,bdir,lamv,mpi_enreg,ndat,nnlout,&
       & paw_opt,signs,svectoutb,tim_getghc,ket,vectout)

     if (need_ormesh) then
       call fourwf(fourwf_cplex,denpot,svectoutb,fofgout,fofrb,gs_hamk%gbound_k,&
         & gs_hamk%gbound_k,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,&
         & gs_hamk%mgfft,mpi_enreg,ndat,gs_hamk%ngfft,npwsp,npwsp,&
         & gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,fourwf_option,&
         & tim_fourwf,weight_r,weight_i)
     end if
     
     do gdir = 1, 3

       ! compute dS/dk_g |u_nk>
       call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,gdir,lamv,mpi_enreg,ndat,nnlout,&
         & paw_opt,signs,svectoutg,tim_getghc,ket,vectout)

       if (need_ormesh) then
         call fourwf(fourwf_cplex,denpot,svectoutg,fofgout,fofrg,gs_hamk%gbound_k,&
           & gs_hamk%gbound_k,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,&
           & gs_hamk%mgfft,mpi_enreg,ndat,gs_hamk%ngfft,npwsp,npwsp,&
           & gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,fourwf_option,&
           & tim_fourwf,weight_r,weight_i)
       end if
   
       ! extract |Pc du/dk_b>
       brab => gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,bdir)
       gdot=cg_zdotc(npwsp,brab,svectoutg); gdotc=CMPLX(gdot(1),gdot(2))

       ! extract |Pc du/dk_g>
       brag => gcg1_k(1:2,(nn-1)*npwsp+1:nn*npwsp,gdir)
       bdot=cg_zdotc(npwsp,brag,svectoutb); bdotc=CMPLX(bdot(1),bdot(2))
        
       do adir=1,3 
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg
         
         ! add <Pc du/dk_b|dS/dk_g|u_nk>*E_nk
         b1(adir) = b1(adir) - prefac_b*gdotc
         m1(adir) = m1(adir) + prefac_m*gdotc*eig_k(nn)
         m1_mu(adir) = m1_mu(adir) - prefac_m*gdotc*fermie

         ! add CONJG(<Pc du/dk_b|dS/dk_g|u_nk>)*E_nk
         b1(adir) = b1(adir) - prefac_b*CONJG(bdotc)
         m1(adir) = m1(adir) + prefac_m*CONJG(bdotc)*eig_k(nn)
         m1_mu(adir) = m1_mu(adir) - prefac_m*CONJG(bdotc)*fermie

         !if (need_ormesh) then
         !  ormesh_fac=trnrm(nn)*prefac_m*(eig_k(nn)-fermie)
         !  call orbmag_mesh%accum_rmesh(adir,brab,dtset,fofrg,gs_hamk,npw_k,&
         !    & ph1d,ormesh_fac,t_atom,invv1,conjg_flag=.FALSE.)
         !  call orbmag_mesh%accum_rmesh(adir,brag,dtset,fofrb,gs_hamk,npw_k,&
         !    & ph1d,ormesh_fac,t_atom,invv1,conjg_flag=.TRUE.)
         !endif

       end do

       do np = 1, nband_k
         if (occ_k(np).LT.tol8) cycle
         bra => cg_k(1:2,(np-1)*npwsp+1:np*npwsp)
         gpdot=cg_zdotc(npwsp,bra,svectoutg); gpdotc=CMPLX(gpdot(1),gpdot(2))
         bpdot=cg_zdotc(npwsp,bra,svectoutb); bpdotc=CMPLX(bpdot(1),bpdot(2))

         !if (need_ormesh) then
         !  ABI_MALLOC(fofgin,(2,npwsp))
         !  fofgin(1,1:npwsp) = bra(1,1:npwsp)*svectoutg(1,1:npwsp) + bra(2,1:npwsp)*svectoutg(2,1:npwsp)
         !  fofgin(2,1:npwsp) = bra(1,1:npwsp)*svectoutg(2,1:npwsp) - bra(2,1:npwsp)*svectoutg(1,1:npwsp)
         !  !fofgin(1:2,1:npwsp)=trnrm(np)*fofgin(1:2,1:npwsp)
         !  call fourwf(fourwf_cplex,denpot,fofgin,fofgout,fofrg,gs_hamk%gbound_k,&
         !    & gs_hamk%gbound_k,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,&
         !    & gs_hamk%mgfft,mpi_enreg,ndat,gs_hamk%ngfft,npwsp,npwsp,&
         !    & gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,fourwf_option,&
         !    & tim_fourwf,weight_r,weight_i)
         !  ABI_SFREE(fofgin)
         !end if

         do adir=1,3 
           epsabg = eijk(adir,bdir,gdir)
           if (ABS(epsabg) .LT. half) cycle
           prefac_b = cbc*c2*epsabg
           prefac_m = com*c2*epsabg
           ! terms in <u|dS|u'><u'|dS|u>
           bv2b(adir) = bv2b(adir) + prefac_b*CONJG(bpdotc)*gpdotc
           mv2b(adir) = mv2b(adir) - prefac_m*CONJG(bpdotc)*gpdotc*eig_k(nn)
           mv2b_mu(adir) = mv2b_mu(adir) + prefac_m*CONJG(bpdotc)*gpdotc*fermie
         
           !if (need_ormesh) then
           !  eig_shift = CMPLX(eig_k(nn)-fermie,zero)
           !  call orbmag_mesh%accum_rmesh(adir,svectoutb,fofrg,gs_hamk%n4,gs_hamk%n5,gs_hamk%n6,&
           !    & dtset%natom,npw_k,gs_hamk%ph3d_k,-prefac_m,t_atom,invv2,&
           !    & mult_fact=eig_shift,conjg_flag=.FALSE.)
           !endif

         end do
       end do ! np

     end do !gdir
   end do !bdir

   orbmag_mesh%cmesh(nn,ikpt,isppol,1:3,ibvv1) = real(b1(1:3))
   orbmag_mesh%cmesh(nn,ikpt,isppol,1:3,ibvv2) = real(bv2b(1:3))

   orbmag_mesh%omesh(nn,ikpt,isppol,1:3,invv1) = real(m1(1:3)+m1_mu(1:3))
   orbmag_mesh%omesh(nn,ikpt,isppol,1:3,invv2) = real(mv2b(1:3)+mv2b_mu(1:3))

 end do !nn

 IF(ASSOCIATED(brab)) NULLIFY(brab)
 IF(ASSOCIATED(brag)) NULLIFY(brag)
 IF(ASSOCIATED(ket)) NULLIFY(ket)
 IF(ASSOCIATED(bra)) NULLIFY(bra)
 
 ABI_FREE(svectoutb)
 ABI_FREE(svectoutg)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 ABI_SFREE(fofrb)
 ABI_SFREE(fofrg)

end subroutine orbmag_vv_k
!!***

!!****f* ABINIT/para_to_diag
!! NAME
!! para_to_diag
!!
!! FUNCTION
!! convert cg1_k wavefunction from parallel to diagonal gauge
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
!!  gcg1_k(2,mcgk,3)=cg1_k converted to requested gauge and/or projection
!!
!! NOTES
!!
!! SOURCE

subroutine para_to_diag(atindx,cg_k,cg1_k,cprj_k,dimlmn,dkinpw,dtset,eig_k,gcg1_k,gs_hamk,&
    & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,mpw,nband_k,ngfft4,ngfft5,ngfft6,npw_k,&
    & nucdip_dirs,occ_k,vectornd_pac)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpw,nband_k,ngfft4,ngfft5,ngfft6
  integer,intent(in) :: npw_k,nucdip_dirs
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),cg1_k(2,mcgk,3),eig_k(nband_k),dkinpw(mpw,3),occ_k(nband_k)
  real(dp),intent(in) :: vectornd_pac(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,nucdip_dirs)
  real(dp),intent(out) :: gcg1_k(2,mcgk,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars
  integer :: adir,berryopt,cplex,iband,ipert,jband,ndat,npwsp
  integer :: optlocal,optnl,opt_gvnlx1,sij_opt,tim_getgh1c,usevnl
  real(dp) :: corrfac,deltae,deltapert,doti,dotr,pertr,perti,pertsize
  type(rf_hamiltonian_type) :: rf_hamk
  !arrays
  real(dp) :: hij(2),lambda(1),sij(2)
  real(dp),allocatable :: cwavef(:,:),dcg1(:,:),gh1c(:,:)
  real(dp),allocatable :: grad_berry(:,:),gs1c(:,:),gvnlx1(:,:)
  real(dp),allocatable :: vectornd_pac_idir(:,:,:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

  berryopt=0
  ndat=1
  optlocal=0
  optnl=2
  opt_gvnlx1=0
  sij_opt=1
  tim_getgh1c=0
  usevnl = 0
  lambda(1) = zero

  ipert=dtset%natom+1 ! DDK
  cplex=1 ! real space 1-order functions on FFT grid are REAL 
  call rf_hamk%init(cplex,gs_hamk,ipert)

  npwsp = npw_k*dtset%nspinor

  ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
  call pawcprj_alloc(cwaveprj,3,dimlmn)
  ABI_MALLOC(cwavef,(2,npwsp))
  ABI_MALLOC(gh1c,(2,gs_hamk%npw_kp*gs_hamk%nspinor*ndat))
  ABI_MALLOC(gs1c,(2,gs_hamk%npw_kp*gs_hamk%nspinor*ndat))
  ABI_MALLOC(gvnlx1,(2,gs_hamk%npw_kp*gs_hamk%nspinor*ndat))
  ABI_MALLOC(dcg1,(2,npwsp))

  if (nucdip_dirs .EQ. 3) then
    ABI_MALLOC(vectornd_pac_idir,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
  end if

  gcg1_k = zero

  do adir = 1, 3

    call rf_hamk%load_k(dkinpw_k=dkinpw(:,adir))

    if (nucdip_dirs .EQ. 3) then
      vectornd_pac_idir(:,:,:,:)=vectornd_pac(:,:,:,:,adir)
      call rf_hamk%load_spin(isppol, vectornd=vectornd_pac_idir)
    end if

    do iband = 1, nband_k

      cwavef(1:2,1:npwsp)=cg_k(1:2,(iband-1)*npwsp+1:iband*npwsp)

      call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,iband,0,ikpt,0,isppol,dtset%mband,&
        & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

      call getgh1c(berryopt,cwavef,cwaveprj,gh1c,grad_berry,gs1c,gs_hamk,gvnlx1,adir,ipert,&
        & lambda(1),mpi_enreg,ndat,optlocal,optnl,opt_gvnlx1,rf_hamk,sij_opt,&
        & tim_getgh1c,usevnl)

      dcg1=zero
      do jband = 1, nband_k
        if (jband .EQ. iband) cycle
        deltae = eig_k(jband) - eig_k(iband)
        ! deltae test seems to work best compared to deltapert test
        !if (abs(deltae) .LT. dtset%ggtrcut) cycle
        cwavef(1:2,1:npwsp)=cg_k(1:2,(jband-1)*npwsp+1:jband*npwsp)
        hij=cg_zdotc(npwsp,cwavef,gh1c)
        sij=cg_zdotc(npwsp,cwavef,gs1c)
        select case (dtset%orbmag)
        case ( 3 )
          lambda(1) = half*(eig_k(jband)+eig_k(iband))
          corrfac=-one
        case ( -3 )
          lambda(1) = eig_k(iband)
          corrfac=-one
        end select
        pertr = (hij(1)-lambda(1)*sij(1))/deltae
        perti = (hij(2)-lambda(1)*sij(2))/deltae
        pertsize=sqrt(pertr*pertr+perti*perti)
        if (pertsize .GT. dtset%ggtrcut) cycle
        dcg1(1,:) = dcg1(1,:) + pertr*cwavef(1,:) - perti*cwavef(2,:)
        dcg1(2,:) = dcg1(2,:) + pertr*cwavef(2,:) + perti*cwavef(1,:)
      end do
      gcg1_k(1:2,(iband-1)*npwsp+1:iband*npwsp,adir) =&
        &cg1_k(1:2,(iband-1)*npwsp+1:iband*npwsp,adir)+corrfac*dcg1(1:2,1:npwsp)
    end do
  end do

  call rf_hamk%free()
  if(allocated(vectornd_pac_idir)) then
    ABI_FREE(vectornd_pac_idir)
  end if
  ABI_FREE(cwavef)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)
  ABI_FREE(gh1c)
  ABI_FREE(gs1c)
  ABI_FREE(gvnlx1)
  ABI_FREE(dcg1)

end subroutine para_to_diag
!!***

!!****f* ABINIT/gauge_treatment
!! NAME
!! gauge_treatment
!!
!! FUNCTION
!! convert cg1_k wavefunction to requested gauge and/or projection
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
!!  gcg1_k(2,mcgk,3)=cg1_k converted to requested gauge and/or projection
!!
!! NOTES
!!
!! SOURCE

subroutine gauge_treatment(atindx,cg_k,cg1_k,cprj_k,dimlmn,dkinpw,dtset,eig_k,gcg1_k,gs_hamk,&
    & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,mpw,nband_k,ngfft4,ngfft5,ngfft6,npw_k,&
    & nucdip_dirs,occ_k,vectornd_pac)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpw,nband_k,ngfft4,ngfft5,ngfft6
  integer,intent(in) :: npw_k,nucdip_dirs
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),cg1_k(2,mcgk,3),eig_k(nband_k),dkinpw(mpw,3),occ_k(nband_k)
  real(dp),intent(in) :: vectornd_pac(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,nucdip_dirs)
  real(dp),intent(out) :: gcg1_k(2,mcgk,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars

!--------------------------------------------------------------------

  ! orbmag > 0: cg1_k contains PAW DDK in parallel gauge, which has a ground state part
  ! orbmag < 0: cg1_k contains Berry phase DDK, which is projected onto the 
  !             conduction space by construction

  select case (dtset%orbmag)
  case ( -3 )
    ! Convert Berry DDK to diagonal gauge
    call para_to_diag(atindx,cg_k,cg1_k,cprj_k,dimlmn,dkinpw,dtset,eig_k,gcg1_k,gs_hamk,&
      & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,mpw,nband_k,ngfft4,ngfft5,ngfft6,npw_k,&
      & nucdip_dirs,occ_k,vectornd_pac)
  case ( -2:-1 )
    ! Berry DDK already projected onto conduction space, leave in berry gauge
    gcg1_k(1:2,1:mcgk,1:3) = cg1_k(1:2,1:mcgk,1:3)
  case ( 3 ) 
    ! Convert PAW DDK to diagonal gauge
    call para_to_diag(atindx,cg_k,cg1_k,cprj_k,dimlmn,dkinpw,dtset,eig_k,gcg1_k,gs_hamk,&
      & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,mpw,nband_k,ngfft4,ngfft5,ngfft6,npw_k,&
      & nucdip_dirs,occ_k,vectornd_pac)
  case default
    ! project cg1_k onto conduction space by removing ground PAW part; 
    ! stay in parallel transport gauge
    call make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,&
      & ikpt,isppol,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,gcg1_k)
  end select

end subroutine gauge_treatment
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
  !arrays
  real(dp) :: dotp(2),lambda(1)
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
        dotp=cg_zdotc(npwsp,cwavef,svectout)
        vcg1(1,:) = vcg1(1,:) - half*( dotp(1)*cwavef(1,:) - dotp(2)*cwavef(2,:))
        vcg1(2,:) = vcg1(2,:) - half*( dotp(1)*cwavef(2,:) + dotp(2)*cwavef(1,:))
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
!!  bdir=direction of bcp derivative to use
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gdir=direction of kcp derivative to use
!!  kcp(dtset%natom,dtset%nspinor)<type(pawcprj_type)>=ket side cprj <p|ket>
!!  lmn2max=max value of lmn2 over all psps
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ndij=spin channels in dij
!!  ucprj(dtset%natom,dtset%nspinor)<type(pawcprj_type)> input cprj
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

subroutine txt_me(adir,atindx,bdir,cwavef,dterm,dtset,eig_k,gdir,gs_hamk,iband,&
    & mpi_enreg,nband_k,npw_k,orbmag_mesh,oterm,pawtab,ph1d,prefac_m,txt,trnrm,ucprj,&
    & suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,bdir,gdir,iband,nband_k,npw_k,oterm
  real(dp),intent(in) :: trnrm
  complex(dp),intent(in) :: prefac_m
  complex(dp),intent(out) :: txt
  logical,intent(in),optional :: suppress_ormesh
  type(dataset_type),intent(in) :: dtset
  type(dterm_type),intent(in) :: dterm
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  real(dp),intent(in),pointer :: cwavef(:,:)
  real(dp),intent(in),pointer :: ph1d(:,:)
  type(pawcprj_type),intent(in) :: ucprj(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,il,ilmn,isp,jl,jlmn,klmn,npwsp,t_atom
  complex(dp) :: dcpi,dcpj,dij,ormesh_fac
  logical :: my_suppress_ormesh,need_ormesh
  !arrays
  real(dp),allocatable,target :: bra(:,:),ket(:,:)
!--------------------------------------------------------------------

  if(present(suppress_ormesh)) then
    my_suppress_ormesh=suppress_ormesh
  else
    my_suppress_ormesh=.FALSE.
  end if
  npwsp = npw_k*dtset%nspinor
  need_ormesh = ((dtset%orbmag .EQ. 4).AND.(.NOT. my_suppress_ormesh))

  if (need_ormesh) then
    ABI_CHECK(ASSOCIATED(cwavef),"txt_me: input wavefunction needed for ormesh is not associated")
    ABI_CHECK(dtset%nspinor.EQ.1,"txt_me: orbmag_rmesh not coded for spinors yet")
    do iat = 1, dtset%natom
      if ( ANY(ABS(dtset%nucdipmom(1:3,iat))>tol8) ) then
        t_atom = atindx(iat)
        exit
      end if
    end do
    ABI_MALLOC(bra,(2,npwsp))
    ABI_MALLOC(ket,(2,npwsp))
  end if

  txt = czero
  do iat = 1, dtset%natom
    itypat=dtset%typat(iat)
    iatom = atindx(iat)
    do isp = 1, dtset%nspinor
      do jlmn = 1, pawtab(itypat)%lmn_size
  
        if (need_ormesh .AND. iatom.EQ.t_atom) then
          ket(1,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1+gdir,jlmn,itypat)*cwavef(1,1:npwsp)
          ket(2,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1+gdir,jlmn,itypat)*cwavef(2,1:npwsp)
        end if
  
        do ilmn = 1, pawtab(itypat)%lmn_size
          klmn = MATPACK(ilmn,jlmn)
          ! in ndij = 4 case, isp 1 delivers up-up, isp 2 delivers down-down
          dij = dterm%aij(iatom,klmn,isp)-eig_k(iband)*dterm%qij(iatom,klmn,isp)
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          dcpi = CMPLX(ucprj(iatom,isp)%dcp(1,bdir,ilmn),ucprj(iatom,isp)%dcp(2,bdir,ilmn))
          dcpj = CMPLX(ucprj(iatom,isp)%dcp(1,gdir,jlmn),ucprj(iatom,isp)%dcp(2,gdir,jlmn))
          txt = txt + prefac_m*CONJG(dcpi)*dij*dcpj

          if (need_ormesh .AND. iatom.EQ.t_atom) then
            jl = pawtab(itypat)%indlmn(1,jlmn)
            il = pawtab(itypat)%indlmn(1,ilmn)
            ormesh_fac = trnrm*prefac_m*dij*four_pi*CONJG(j_dpc**il)*four_pi*(j_dpc**jl)
            bra(1,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1+bdir,ilmn,itypat)*cwavef(1,1:npwsp)
            bra(2,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1+bdir,ilmn,itypat)*cwavef(2,1:npwsp)
            call orbmag_mesh%accum_rmesh(adir,bra,dtset,gs_hamk,ket,.FALSE.,mpi_enreg,&
              & npwsp,ph1d,ormesh_fac,t_atom,oterm)
          end if
          
          if (dterm%ndij == 4) then
            if (isp == 1) then
              dij = dterm%aij(iatom,klmn,3)-eig_k(iband)*dterm%qij(iatom,klmn,3) ! up-down
              ! D^ss'_ij=D^s's_ji^*
              !if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,4))
              if (ilmn .GT. jlmn) then
                dij = CONJG(dterm%aij(iatom,klmn,4)-eig_k(iband)*dterm%qij(iatom,klmn,4)) 
              end if
              dcpi = CMPLX(ucprj(iatom,1)%dcp(1,bdir,ilmn),ucprj(iatom,1)%dcp(2,bdir,ilmn))
              dcpj = CMPLX(ucprj(iatom,2)%dcp(1,gdir,jlmn),ucprj(iatom,2)%dcp(2,gdir,jlmn))
            else
              dij = dterm%aij(iatom,klmn,4)-eig_k(iband)*dterm%qij(iatom,klmn,4) ! down-up
              ! D^ss'_ij=D^s's_ji^*
              !if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,3))
              if (ilmn .GT. jlmn) then
                dij = CONJG(dterm%aij(iatom,klmn,3)-eig_k(iband)*dterm%qij(iatom,klmn,3))
              end if
              dcpi = CMPLX(ucprj(iatom,2)%dcp(1,bdir,ilmn),ucprj(iatom,2)%dcp(2,bdir,ilmn))
              dcpj = CMPLX(ucprj(iatom,1)%dcp(1,gdir,jlmn),ucprj(iatom,1)%dcp(2,gdir,jlmn))
            end if
            txt = txt + prefac_m*CONJG(dcpi)*dcpj*dij
          end if
        end do !ilmn
      end do !jlmn
    end do ! isp
  end do !iat

  ABI_SFREE(bra)
  ABI_SFREE(ket)

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
!!  cwavef(:,:),pointer=points to current wavefunction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk<type(gs_hamiltonian_type)>=ground state Hamiltonian at this k
!!  lmn2max=max value of lmn2 over all psps
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ndij=spin channels in dij
!!  ucprj(dtset%natom,dtset%nspinor)<type(pawcprj_type)> input cprj
!!
!! OUTPUT
!! tt=(complex) computed matrix element
!!
!! NOTES
!! computes on-site \sum_{Rij}<u|p_i>aij<p_j|u> for generic aij input
!!
!! SOURCE

subroutine tt_me(adir,aij,atindx,cwavef,dtset,gs_hamk,lmn2max,mpi_enreg,&
    & ndij,nband_k,npw_k,orbmag_mesh,oterm,pawtab,ph1d,tt,trnrm,ucprj,&
    & suppress_ormesh)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,lmn2max,ndij,nband_k,npw_k,oterm
  real(dp),intent(in) :: trnrm
  complex(dp),intent(out) :: tt
  logical,intent(in),optional :: suppress_ormesh
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_mesh_type),intent(inout) :: orbmag_mesh

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in),pointer :: cwavef(:,:),ph1d(:,:)
  complex(dp),intent(in) :: aij(dtset%natom,lmn2max,ndij)
  type(pawcprj_type),intent(in) :: ucprj(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,ig,isp,itypat,il,ilmn,ipw,jl,jlmn,klmn,n4,n5,n6,npwsp
  integer :: t_atom
  logical :: my_suppress_ormesh,need_ormesh
  real(dp) :: weight_i,weight_r
  complex(dp) :: cpi,cpj,dij,ormesh_fac
  !arrays
  real(dp),allocatable,target :: bra(:,:),ket(:,:)
!--------------------------------------------------------------------

  if(present(suppress_ormesh)) then
    my_suppress_ormesh=suppress_ormesh
  else
    my_suppress_ormesh=.FALSE.
  end if
  npwsp = npw_k*dtset%nspinor
  need_ormesh = ((dtset%orbmag .EQ. 4) .AND. (.NOT. my_suppress_ormesh))

  n4=dtset%ngfft(4); n5=dtset%ngfft(5); n6=dtset%ngfft(6)

  if (need_ormesh) then
    ABI_CHECK(ASSOCIATED(cwavef),"tt_me: input wavefunction needed for ormesh is not associated")
    ABI_CHECK(dtset%nspinor.EQ.1,"tt_me: orbmag_rmesh not coded for spinors yet")
    do iat = 1, dtset%natom
      if ( ANY(ABS(dtset%nucdipmom(1:3,iat))>tol8) ) then
        t_atom = atindx(iat)
        exit
      end if
    end do
    ABI_MALLOC(bra,(2,npwsp))
    ABI_MALLOC(ket,(2,npwsp))
  end if

  tt = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat=dtset%typat(iat)
    do isp = 1, dtset%nspinor
      do jlmn = 1, pawtab(itypat)%lmn_size
  
        if (need_ormesh .AND. iatom.EQ.t_atom) then
          ! FFT ket-side to fofr, real space representation
          ket(1,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1,jlmn,itypat)*cwavef(1,1:npwsp)
          ket(2,1:npwsp) = gs_hamk%ffnl_k(1:npwsp,1,jlmn,itypat)*cwavef(2,1:npwsp)
        end if
 
        do ilmn = 1, pawtab(itypat)%lmn_size
          klmn=MATPACK(ilmn,jlmn)
          ! in ndij = 4 case, isp 1 delivers up-up, isp 2 delivers down-down
          dij = aij(iatom,klmn,isp)
          cpi =  CMPLX(ucprj(iatom,isp)%cp(1,ilmn),ucprj(iatom,isp)%cp(2,ilmn))
          cpj =  CMPLX(ucprj(iatom,isp)%cp(1,jlmn),ucprj(iatom,isp)%cp(2,jlmn))
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          ! note use of CONJG(cpi), because cpi is from the bra side cprj
          tt = tt + CONJG(cpi)*dij*cpj
          
          if (need_ormesh .AND. iatom.EQ.t_atom) then
            jl = pawtab(itypat)%indlmn(1,jlmn)
            il = pawtab(itypat)%indlmn(1,ilmn)
            ormesh_fac = trnrm*dij*four_pi*CONJG(j_dpc**il)*four_pi*(j_dpc**jl)
            
            bra(1,1:npwsp)=cwavef(1,1:npwsp)*gs_hamk%ffnl_k(1:npwsp,1,ilmn,itypat)
            bra(2,1:npwsp)=cwavef(2,1:npwsp)*gs_hamk%ffnl_k(1:npwsp,1,ilmn,itypat)
            call orbmag_mesh%accum_rmesh(adir,bra,dtset,gs_hamk,ket,.FALSE.,mpi_enreg,&
              & npwsp,ph1d,ormesh_fac,t_atom,oterm)
          end if
            
          if (ndij == 4) then
            if (isp == 1) then
              dij = aij(iatom,klmn,3) ! up-down
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,4))
              cpi = CMPLX(ucprj(iatom,1)%cp(1,ilmn),ucprj(iatom,1)%cp(2,ilmn))
              cpj = CMPLX(ucprj(iatom,2)%cp(1,jlmn),ucprj(iatom,2)%cp(2,jlmn))
            else
              dij = aij(iatom,klmn,4) ! down-up
              ! D^ss'_ij=D^s's_ji^*
              if (ilmn .GT. jlmn) dij = CONJG(aij(iatom,klmn,3))
              cpi = CMPLX(ucprj(iatom,2)%cp(1,ilmn),ucprj(iatom,2)%cp(2,ilmn))
              cpj = CMPLX(ucprj(iatom,1)%cp(1,jlmn),ucprj(iatom,1)%cp(2,jlmn))
            end if
            tt = tt + CONJG(cpi)*cpj*dij
          end if
        end do !ilmn
      end do !jlmn
    end do ! isp
  end do !iat

  ABI_SFREE(bra)
  ABI_SFREE(ket)

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
!! Only printing. This routine outputs orbmag terms to the normal abinit output file
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  omlamb(3)=Lamb shielding
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine orbmag_output(omag,dtset,omlamb)

 !Arguments ------------------------------------
 !scalars
 class(orbmag_mesh_type),intent(inout),target :: omag
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: omlamb(3)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,isppol,iterms
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(dtset%mband,3),berry_total(3),orbmag_bb(dtset%mband,3),orbmag_total(3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = 1,orbmag_nterms
   orbmag_total(1:3)=orbmag_total(1:3) + omag%orbmag_trace(1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       orbmag_bb(iband,1:3) = orbmag_bb(iband,1:3) + omag%orbmag_terms(iband,isppol,1:3,iterms)
     end do ! iband
   end do ! isppol
 end do

 orbmag_total=orbmag_total+omlamb

 berry_bb=zero;berry_total=zero
 do iterms = 1,chern_nterms
   berry_total(1:3)=berry_total(1:3) + omag%chern_trace(1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       berry_bb(iband,1:3) = berry_bb(iband,1:3) + omag%chern_terms(iband,isppol,1:3,iterms)
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

 if(abs(dtset%orbmag) .GE. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(1) CC : ',(omag%orbmag_trace(adir,incc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(omag%orbmag_trace(adir,invv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(omag%orbmag_trace(adir,invv2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(0) NL : ',(omag%orbmag_trace(adir,innl),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <L_R> terms : ',(omag%orbmag_trace(adir,inlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <A0.An> terms : ',(omag%orbmag_trace(adir,inbm),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    Lamb terms : ',(omlamb(adir),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Chern vector, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  Ch CC : ',(omag%chern_trace(adir,ibcc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv1 : ',(omag%chern_trace(adir,ibvv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv2 : ',(omag%chern_trace(adir,ibvv2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
 end if

 if(abs(dtset%orbmag) .GE. 2) then
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
       write(message,'(a,3es16.8)') '     rho(1) CC : ',(omag%orbmag_terms(iband,isppol,adir,incc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(omag%orbmag_terms(iband,isppol,adir,invv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(omag%orbmag_terms(iband,isppol,adir,invv2),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '     rho(0) NL : ',(omag%orbmag_terms(iband,isppol,adir,innl),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '   <L_R> terms : ',(omag%orbmag_terms(iband,isppol,adir,inlr),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' <A0.An> terms : ',(omag%orbmag_terms(iband,isppol,adir,inbm),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a)')ch10
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Chern vector : ',(berry_bb(iband,adir),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '  Ch CC : ',(omag%chern_terms(iband,isppol,adir,ibcc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV1 : ',(omag%chern_terms(iband,isppol,adir,ibvv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV2 : ',(omag%chern_terms(iband,isppol,adir,ibvv2),adir=1,3)
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
  class(dterm_type),intent(inout),target :: dterm
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

!!****f* ABINIT/dterm_init
!! NAME
!! dterm_init
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

subroutine dterm_init(dterm,lmnmax,lmn2max,natom,ndij)

  !Arguments ------------------------------------
  !scalars
  class(dterm_type),intent(inout),target :: dterm
  integer,intent(in) :: lmnmax,lmn2max,natom,ndij
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

end subroutine dterm_init
!!***

!!****f* ABINIT/orbmag_init
!! NAME
!! orbmag_init
!!
!! FUNCTION
!! allocate space in orbmag_mesh_type
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! orbmag_mesh <type(orbmag_mesh_type)> data related to orbmag terms on kpt mesh
!!
!! SOURCE

subroutine orbmag_init(omag,dtset)

  !Arguments ------------------------------------
  !scalars
  class(orbmag_mesh_type),intent(inout),target :: omag
  type(dataset_type),intent(in) :: dtset
!--------------------------------------------------------------------

  omag%mband = dtset%mband
  omag%nkpt = dtset%nkpt
  omag%nsppol = dtset%nsppol
  omag%natom = dtset%natom
  omag%ntypat = dtset%ntypat
  omag%n4=dtset%ngfft(4)
  omag%n5=dtset%ngfft(5)
  omag%n6=dtset%ngfft(6)
  omag%chern_nterms = chern_nterms
  omag%orbmag_nterms = orbmag_nterms

  ABI_REMALLOC(omag%lambsig,(omag%ntypat))
  ABI_REMALLOC(omag%nucdipmom,(3,omag%natom))
  ABI_REMALLOC(omag%cmesh,(omag%mband,omag%nkpt,omag%nsppol,3,chern_nterms))
  omag%cmesh=zero
  ABI_REMALLOC(omag%chern_terms,(dtset%mband,dtset%nsppol,3,chern_nterms))
  omag%chern_terms=zero
  ABI_REMALLOC(omag%chern_trace,(3,chern_nterms))
  omag%chern_trace=zero
  ABI_REMALLOC(omag%omesh,(omag%mband,omag%nkpt,omag%nsppol,3,orbmag_nterms))
  omag%omesh=zero
  ABI_REMALLOC(omag%orbmag_terms,(dtset%mband,dtset%nsppol,3,orbmag_nterms))
  omag%orbmag_terms=zero
  ABI_REMALLOC(omag%orbmag_trace,(3,orbmag_nterms))
  omag%orbmag_trace=zero
  if (dtset%orbmag .EQ. 4) then
    ABI_REMALLOC(omag%rmesh,(omag%n4,omag%n5,omag%n6,3,orbmag_nterms))
    omag%omesh=zero
  end if

end subroutine orbmag_init
!!***

!!****f* ABINIT/orbmag_free
!! NAME
!! orbmag_free
!!
!! FUNCTION
!! free space in orbmag_mesh_type
!!
!! SIDE EFFECTS
!! orbmag_mesh <type(orbmag_mesh_type)> data related to orbmag terms on kpt mesh
!!
!! SOURCE

subroutine orbmag_free(omag)

  !Arguments ------------------------------------
  !scalars
  class(orbmag_mesh_type),intent(inout),target :: omag
!--------------------------------------------------------------------

    ABI_SFREE(omag%lambsig)
    ABI_SFREE(omag%nucdipmom)
    ABI_SFREE(omag%cmesh)
    ABI_SFREE(omag%chern_terms)
    ABI_SFREE(omag%chern_trace)
    ABI_SFREE(omag%omesh)
    ABI_SFREE(omag%orbmag_terms)
    ABI_SFREE(omag%orbmag_trace)
    ABI_SFREE(omag%rmesh)

end subroutine orbmag_free
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

!!****f* ABINIT/orbmag_rmesh
!! NAME
!! orbmag_rmesh
!!
!! FUNCTION
!! accumulate orbmag density on real mesh into orbmag structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! orbmag_mesh%rmesh updated
!! 
!! CHILDREN
!!
!! SOURCE

subroutine orbmag_rmesh(omag,adir,bra,dtset,gs_hamk,ket,local_term,mpi_enreg,&
    & npw_k,ph1d,scalar_factor,t_atom,term_index)

  !Arguments ------------------------------------
  !scalars
  class(orbmag_mesh_type),intent(inout),target :: omag
  integer,intent(in) :: adir,npw_k,t_atom,term_index
  complex(dp),intent(in) :: scalar_factor
  logical,intent(in) :: local_term
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg
  !arrays
  real(dp),intent(in),pointer :: bra(:,:),ket(:,:),ph1d(:,:)

  !Local variables -------------------------
  !scalars
  integer :: fourwf_cplex,fourwf_option,ig,kg1,kg2,kg3,n1,n2,n3,n4,n5,n6,ndat
  integer :: shift1,shift2,shift3,tim_fourwf
  real(dp) :: weight_i,weight_r
  complex(dp) :: cpw,ph1,ph2,ph3
  !arrays
  real(dp),allocatable :: denpot(:,:,:),fofgout(:,:),fofr(:,:,:,:),work(:,:)
  complex(dp),allocatable :: phgr(:)

!--------------------------------------------------------------------

  ndat=1
  n1=dtset%ngfft(1); n2=dtset%ngfft(2); n3=dtset%ngfft(3)
  n4=dtset%ngfft(4); n5=dtset%ngfft(5); n6=dtset%ngfft(6)
  ABI_MALLOC(fofr,(2,n4,n5,n6*ndat))

  shift1=1+n1+(t_atom-1)*(2*n1+1)
  shift2=1+n2+(t_atom-1)*(2*n2+1)+dtset%natom*(2*n1+1)
  shift3=1+n3+(t_atom-1)*(2*n3+1)+dtset%natom*(2*n1+1+2*n2+1) 
  ABI_MALLOC(phgr,(npw_k))
  do ig=1,npw_k
    kg1=gs_hamk%kg_k(1,ig)+shift1
    kg2=gs_hamk%kg_k(2,ig)+shift2
    kg3=gs_hamk%kg_k(3,ig)+shift3
    ph1=CMPLX(ph1d(1,kg1),ph1d(2,kg1))
    ph2=CMPLX(ph1d(1,kg2),ph1d(2,kg2))
    ph3=CMPLX(ph1d(1,kg3),ph1d(2,kg3))
    phgr(ig)=ph1*ph2*ph3
  end do

  if (local_term) then
    ! if local field, compute exp(-iG.R)*conjg(bra)*scalar_factor*ket and then transform with fourwf
    ABI_MALLOC(work,(2,npw_k))
    do ig = 1, npw_k
      cpw = CONJG(phgr(ig))*CONJG(CMPLX(bra(1,ig),bra(2,ig)))*scalar_factor*CMPLX(ket(1,ig),ket(2,ig))
      work(1,ig) = REAL(cpw); work(2,ig) = AIMAG(cpw)
    end do
    fourwf_cplex = 1
    fourwf_option = 0
    tim_fourwf = 1
    call fourwf(fourwf_cplex,denpot,work,fofgout,fofr,gs_hamk%gbound_k,&
      & gs_hamk%gbound_k,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,&
      & gs_hamk%mgfft,mpi_enreg,ndat,gs_hamk%ngfft,npw_k,npw_k,&
      & n4,n5,n6,fourwf_option,tim_fourwf,weight_r,weight_i)
    ABI_SFREE(work)
  else
    ! if nonlocal, transform ket with fourwf, then slow FT as scalar_factor*(sum_G exp(+iG.R)*conjg(bra))*fofr
    fourwf_cplex = 1
    fourwf_option = 0
    tim_fourwf = 1
    call fourwf(fourwf_cplex,denpot,ket,fofgout,fofr,gs_hamk%gbound_k,&
      & gs_hamk%gbound_k,gs_hamk%istwf_k,gs_hamk%kg_k,gs_hamk%kg_k,&
      & gs_hamk%mgfft,mpi_enreg,ndat,gs_hamk%ngfft,npw_k,npw_k,&
      & n4,n5,n6,fourwf_option,tim_fourwf,weight_r,weight_i)
    cpw = czero
    do ig = 1, npw_k
      cpw = cpw + phgr(ig)*CONJG(CMPLX(bra(1,ig),bra(2,ig)))
    end do
    cpw = cpw*scalar_factor
    fofr(1,:,:,:) = fofr(1,:,:,:)*REAL(cpw) - fofr(2,:,:,:)*AIMAG(cpw)
  end if

  omag%rmesh(:,:,:,adir,term_index)=omag%rmesh(:,:,:,adir,term_index)+fofr(1,:,:,:)
  
  !crr=REAL(crvec);cri=AIMAG(crvec);sfr=REAL(scalar_factor);sfi=AIMAG(scalar_factor)

  !if (term_index.EQ.incc) then
  !  omag%rmesh(:,:,:,adir,term_index) = omag%rmesh(:,:,:,adir,term_index) + &
  !    & fofr(1,:,:,:)
  !else
  !  omag%rmesh(:,:,:,adir,term_index) = omag%rmesh(:,:,:,adir,term_index) + &
  !    & fofr(1,:,:,:)*(sfr*crr+conjg_fac*sfi*cri) + &
  !    & fofr(2,:,:,:)*(-sfr*cri+conjg_fac*sfi*crr)
  !end if

  ABI_SFREE(fofr)
  ABI_SFREE(phgr)

end subroutine orbmag_rmesh
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
 logical :: has_ormesh
 character(len=500) :: msg
!arrays
!*************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 has_ormesh = (dtset%orbmag .EQ. 4)

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

 !! add orbmag_rmesh_cplex,n4,n5,n6 only if orbmag_rmesh will be output
 if (has_ormesh) then
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("n4", orbmag_mesh%n4),&
     nctkdim_t("n5", orbmag_mesh%n5),&
     nctkdim_t("n6", orbmag_mesh%n6)],defmode=.True.)
   NCF_CHECK(ncerr)
 endif

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("chern_mesh", "dp", "mband, nkpt, nsppol, ndir, chern_nterms"),&
   nctkarr_t("orbmag_mesh", "dp", "mband, nkpt, nsppol, ndir, orbmag_nterms"),&
   nctkarr_t("lambsig", "dp", "ntypat"),&
   nctkarr_t("nucdipmom", "dp", "ndir, natom")])
 NCF_CHECK(ncerr)

 !! orbmag_rmesh dimensions, only if output
 if (has_ormesh) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("orbmag_rmesh", "dp", "n4,n5,n6,ndir,orbmag_nterms")])
   NCF_CHECK(ncerr)
 endif


 NCF_CHECK(nctk_set_datamode(ncid))

 NCF_CHECK(nf90_put_var(ncid, vid("chern_mesh"), orbmag_mesh%cmesh))
 NCF_CHECK(nf90_put_var(ncid, vid("orbmag_mesh"), orbmag_mesh%omesh))
 NCF_CHECK(nf90_put_var(ncid, vid("nucdipmom"), orbmag_mesh%nucdipmom))
 NCF_CHECK(nf90_put_var(ncid, vid("lambsig"), orbmag_mesh%lambsig))
 if ( has_ormesh ) then
   NCF_CHECK(nf90_put_var(ncid, vid("orbmag_rmesh"), orbmag_mesh%rmesh))
 end if

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
