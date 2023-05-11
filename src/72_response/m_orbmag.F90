!!*** ABINIT/m_orbmag
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
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : projbd
  use m_fft,              only : fftpac
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type
  use m_kg,               only : getph,mkkin,mkkpg,ph1d3d
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle,proc_distrb_nband
  use m_nonlop,           only : nonlop
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get, pawcprj_put
  use m_pawfgr,           only : pawfgr_type
  use m_pawfgrtab,        only : pawfgrtab_type
  use m_paw_ij,           only : paw_ij_type
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen,poisson
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


  ! these parameters name the various output terms                                             
  integer,parameter :: ibcc=1,ibvv1=2,ibvv2=3
  integer,parameter :: imcc=4,imvv1=5,imvv2=6
  integer,parameter :: imnl=7,imlr=8,imbm=9
  integer,parameter :: iomlmb=10
  integer,parameter :: nterms=10

  ! these parameters are constants used repeatedly

  ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  real(dp),parameter :: c2=one/(two_pi*two_pi)
  complex(dpc),parameter :: com=-half*j_dpc  ! Orbital magnetism pre-factor
  complex(dpc),parameter :: cbc=-com ! Berry curvature pre-factor

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
    complex(dpc),allocatable :: aij(:,:,:)

    ! <phi|phi> - <tphi|tphi>
    ! qij(natom,lmn2max,ndij)
    complex(dpc),allocatable :: qij(:,:,:)

    ! onsite L_R/2
    ! <phi|L_R/2|phi> - <tphi|L_R/2|tphi>
    ! LR(natom,lmn2max,ndij,3)
    complex(dpc),allocatable :: LR(:,:,:,:)

    ! onsite BM
    ! <phi|Bxr . mxr|phi> - <tphi|Bxr . mxr|tphi>
    ! BM(natom,lmn2max,ndij,3)
    complex(dpc),allocatable :: BM(:,:,:,:)

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
  private :: orbmag_output
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
!! COPYRIGHT
!! Copyright (C) 2003-2022 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cg(2,mcg)=all ground state wavefunctions
!!  cg1(2,mcg1,3)=all DDK wavefunctions in all 3 directions
!!  cprj(dtset%natom,mcprj)<type(pawcprj_type)>=all ground state cprj
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)=all ground state eigenvalues
!!  gsqcut=large sphere cut-off
!!  kg(3,mpw*mkmem_rbz)=basis sphere of planewaves at k
!!  mcg=dimension of cg
!!  mcg1=dimension of cg1
!!  mcprj=dimension of cprj
!!  mkmem_rbz=kpts in memory
!!  mpi_enreg<type(MPI_type)>=information about MPI parallelization
!!  mpw=max number of planewaves at k
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(18)=FFT grid size information (from pawfgr%ngfft)
!!  npwarr(dtset%nkpt)=npw_k at each kpt
!!  occ(dtset%mband*dtset%nkpt*dtset%nsppol)=occup number for each band (often 2) at each k point
!!  paw_ij(dtset%natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(dtset%ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=real space translation vectors
!!  vtrial(nfftf,dtset%nspden)=GS potential (Hartree)
!!  xred(3,dtset%natom)=reduced dimensionless atomic coordinates
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
!! See Gonze and Zwanziger, PRB 84, 064446 (2011)
!! DDK wavefunctions are used for the derivatives.
!!
!! SOURCE

subroutine orbmag(cg,cg1,cprj,dtset,eigen0,gsqcut,kg,mcg,mcg1,mcprj,mkmem_rbz,mpi_enreg,mpw,&
    & nfftf,ngfftf,npwarr,occ,paw_ij,pawfgr,pawrad,&
    & pawtab,psps,rprimd,vtrial,xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcprj,mcg,mcg1,mkmem_rbz,mpw,nfftf
 real(dp),intent(in) :: gsqcut
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps

 !arrays
 integer,intent(in) :: kg(3,mpw*mkmem_rbz),ngfftf(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem_rbz,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,bdtot_index,buff_size,choice,cpopt,dimffnl,exchn2n3d
 integer :: iat,iatom,icg,icmplx,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat,lmn2max
 integer :: me,mcgk,mcprjk,my_lmax,my_nspinor,nband_k,nband_me,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nl1_option,nn,nkpg,npw_k,npwsp,nproc,spaceComm,with_vectornd
 real(dp) :: arg,ecut_eff,fermie,ucvol
 logical :: has_nucdip
 type(dterm_type) :: dterm
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: b1_k(:,:,:),b2_k(:,:,:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: eig_k(:),ffnl_k(:,:,:,:),kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: m1_k(:,:,:),m1_mu_k(:,:,:),m2_k(:,:,:),m2_mu_k(:,:,:)
 real(dp),allocatable :: occ_k(:),orbmag_terms(:,:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
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

 write(std_out,'(a,2i12)')'JWZ debug mcg mcg1 ',mcg,mcg1

 ! Fermi energy
 call local_fermie(dtset,eigen0,fermie,mpi_enreg,occ)
 !fermie = dtset%userra
 !write(std_out,'(a,es16.8)')'JWZ debug using fermi e = ',fermie

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

 ABI_MALLOC(kg_k,(3,mpw))
 ABI_MALLOC(kinpw,(mpw))

 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,0,dimlmn)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 my_lmax = psps%mpsang + 1
 ABI_MALLOC(realgnt,((2*my_lmax-1)**2*(my_lmax)**4))
 ABI_MALLOC(gntselect,((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2))
 call realgaunt(my_lmax,ngnt,gntselect,realgnt)

 lmn2max = psps%lmnmax*(psps%lmnmax+1)/2
 ! note: in make_d, terms will be filled as iatom using atindx
 call dterm_alloc(dterm,psps%lmnmax,lmn2max,dtset%natom,paw_ij(1)%ndij)
 call make_d(atindx,dterm,dtset,gprimd,paw_ij,pawrad,pawtab,psps)

 ABI_MALLOC(orbmag_terms,(2,dtset%mband,dtset%nsppol,3,nterms))
 orbmag_terms = zero

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
      & paw_ij=paw_ij)

 ! iterate over spin channels
 bdtot_index=0
 icg = 0
 icprj = 0
 do isppol = 1, dtset%nsppol

   !========= construct local potential ==================
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
   ABI_MALLOC(vectornd,(with_vectornd*nfftf,dtset%nspden,3))
   vectornd = zero
   if(has_nucdip) then
     call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,&
       & dtset%nspden,dtset%nucdipmom,rprimd,vectornd,xred)
     ABI_MALLOC(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
     ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
     do idir = 1, 3
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,:,idir))
       call fftpac(isppol,mpi_enreg,dtset%nspden,&
         & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
     end do
     ABI_FREE(cgrvtrial)
     call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
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

     kpoint(:)=dtset%kptns(:,ikpt)
     npw_k = npwarr(ikpt)
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
     occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)

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

     ! retrieve ground state wavefunctions at this k point and isppol
     mcgk = npw_k*nband_k*dtset%nspinor
     ABI_MALLOC(cg_k,(2,mcgk))
     cg_k = cg(1:2,icg+1:icg+mcgk)

     ! retrieve first order wavefunctions at this k point and isppol
     ABI_MALLOC(cg1_k,(2,mcgk,3))
     cg1_k = cg1(1:2,icg+1:icg+mcgk,1:3)

     ! retrieve zeroth order eigenvalues at this k point and isppol
     ABI_MALLOC(eig_k,(nband_k))
     eig_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)

     ! retrieve cprj_k at this k point and isppol
     mcprjk = nband_k*dtset%nspinor
     ABI_MALLOC(cprj_k,(dtset%natom,mcprjk))
     call pawcprj_alloc(cprj_k,cprj(1,1)%ncpgr,dimlmn)
     call pawcprj_get(atindx,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,nband_k,nband_k,dtset%nspinor,dtset%nsppol,0)

     ! compute P_c|cg1>
     ABI_MALLOC(pcg1_k,(2,mcgk,3))
     call make_pcg1(atindx,cg_k,cg1_k,cprj_k,dimlmn,dtset,gs_hamk,ikpt,isppol,&
       & mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,pcg1_k)

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
           & kg_k,kpg_k,kpoint,psps%lmnmax,dtset%mgfft,mpi_enreg,dtset%natom,nattyp,dtset%ngfft,&
           & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)
         call pawcprj_put(atindx,cwaveprj,cprj1_k(:,:,adir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & mkmem_rbz,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)
       end do
     end do
     ABI_FREE(cwavef)

     !--------------------------------------------------------------------------------
     ! Finally ready to compute contributions to orbital magnetism and Berry curvature
     !--------------------------------------------------------------------------------

     !! check aij against paw_ij
     !call check_eig_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,&
     !  & gs_hamk,ikpt,isppol,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab)

     ABI_MALLOC(b1_k,(2,nband_k,3))
     ABI_MALLOC(b2_k,(2,nband_k,3))
     ABI_MALLOC(m1_k,(2,nband_k,3))
     ABI_MALLOC(m1_mu_k,(2,nband_k,3))
     ABI_MALLOC(m2_k,(2,nband_k,3))
     ABI_MALLOC(m2_mu_k,(2,nband_k,3))

     call orbmag_cc_k(atindx,b1_k,cprj1_k,dimlmn,dtset,eig_k,fermie,gs_hamk,ikpt,isppol,&
       & m1_k,m1_mu_k,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,pcg1_k,ucvol)

     orbmag_terms(:,:,isppol,:,ibcc) = orbmag_terms(:,:,isppol,:,ibcc) + b1_k(:,:,:)

     orbmag_terms(:,:,isppol,:,imcc) = orbmag_terms(:,:,isppol,:,imcc) + m1_k(:,:,:)
     orbmag_terms(:,:,isppol,:,imcc) = orbmag_terms(:,:,isppol,:,imcc) + m1_mu_k(:,:,:)

     call orbmag_vv_k(atindx,b1_k,b2_k,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gs_hamk,&
      & ikpt,isppol,m1_k,m1_mu_k,m2_k,m2_mu_k,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,&
      & npw_k,occ_k,pcg1_k,ucvol)
     orbmag_terms(:,:,isppol,:,ibvv1) = orbmag_terms(:,:,isppol,:,ibvv1) + b1_k(:,:,:)
     orbmag_terms(:,:,isppol,:,ibvv2) = orbmag_terms(:,:,isppol,:,ibvv2) + b2_k(:,:,:)

     orbmag_terms(:,:,isppol,:,imvv1) = orbmag_terms(:,:,isppol,:,imvv1) + m1_k(:,:,:)
     orbmag_terms(:,:,isppol,:,imvv1) = orbmag_terms(:,:,isppol,:,imvv1) + m1_mu_k(:,:,:)

     orbmag_terms(:,:,isppol,:,imvv2) = orbmag_terms(:,:,isppol,:,imvv2) + m2_k(:,:,:)
     orbmag_terms(:,:,isppol,:,imvv2) = orbmag_terms(:,:,isppol,:,imvv2) + m2_mu_k(:,:,:)

     call orbmag_nl_k(atindx,cprj_k,dimlmn,dterm,dtset,eig_k,ikpt,isppol,&
       & m1_k,mcprjk,mkmem_rbz,nband_k,occ_k,pawtab,ucvol)
     orbmag_terms(:,:,isppol,:,imnl) = orbmag_terms(:,:,isppol,:,imnl) + m1_k(:,:,:)

     nl1_option = 1 ! LR
     call orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,m1_k,mcprjk,mkmem_rbz,&
       & nband_k,nl1_option,occ_k,pawtab,ucvol)
     orbmag_terms(:,:,isppol,:,imlr) = orbmag_terms(:,:,isppol,:,imlr) + m1_k
     nl1_option = 2 ! BM
     call orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,m1_k,mcprjk,mkmem_rbz,&
       & nband_k,nl1_option,occ_k,pawtab,ucvol)
     orbmag_terms(:,:,isppol,:,imbm) = orbmag_terms(:,:,isppol,:,imbm) + m1_k

     ABI_FREE(b1_k)
     ABI_FREE(b2_k)
     ABI_FREE(m1_k)
     ABI_FREE(m1_mu_k)
     ABI_FREE(m2_k)
     ABI_FREE(m2_mu_k)

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
   ABI_FREE(vectornd)
   if(has_nucdip) then
     ABI_FREE(vectornd_pac)
   end if

 end do ! end loop over isppol

 !! collect orbmag_terms if distributed over different processes
 if (nproc > 1) then
   buff_size=size(orbmag_terms)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = reshape(orbmag_terms,(/2*nband_k*dtset%nsppol*3*nterms/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   orbmag_terms(1:2,1:nband_k,1:dtset%nsppol,1:3,1:nterms)=reshape(buffer2,(/2,nband_k,dtset%nsppol,3,nterms/))
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
   if ((iterm.EQ.iomlmb)) cycle
   do isppol = 1, dtset%nsppol
     do nn = 1, nband_k
       do icmplx = 1, 2
         if((iterm.EQ.imlr).OR.(iterm.EQ.imbm)) then
           orbmag_terms(icmplx,nn,isppol,1:3,iterm) = &
             & MATMUL(rprimd,orbmag_terms(icmplx,nn,isppol,1:3,iterm))
         else
           orbmag_terms(icmplx,nn,isppol,1:3,iterm) = &
             & ucvol*MATMUL(gprimd,orbmag_terms(icmplx,nn,isppol,1:3,iterm))
         end if
       end do !icmplx
     end do ! nn
   end do !isppol
 end do

 !! convert orbmag magnetization to orbital moment
 !! Berry curvature terms are ignored
 !! Lamb term ignored
 do iterm = 1, nterms
   if (iterm .EQ. ibcc) cycle
   if (iterm .EQ. ibvv1) cycle
   if (iterm .EQ. ibvv2) cycle
   if (iterm .EQ. iomlmb) cycle
   orbmag_terms(:,:,:,:,iterm) = ucvol*orbmag_terms(:,:,:,:,iterm)
 end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do isppol = 1, dtset%nsppol
   do nn = 1, nband_k
     orbmag_trace = orbmag_trace + orbmag_terms(:,nn,isppol,:,:)
   end do ! nn
 end do ! isppol

! get the Lamb term
call lamb_core(atindx,dtset,omlamb,pawtab)
orbmag_trace(:,:,iomlmb) = omlamb

call orbmag_output(dtset,orbmag_terms,orbmag_trace)

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

 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

end subroutine orbmag
!!***

!!****f* ABINIT/orbmag_nl1_k
!! NAME
!! orbmag_nl1_k
!!
!! FUNCTION
!! make NL(1) term at k
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  occ_k=band occupations at this kpt
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  if nl1_option = 1, orbmag contribution of <L_R> is returned
!!  if nl1_option = 2, orbmag contribution of <A0.An> is returned
!!    m1_k(2,nband_k,3)=Orb mag contribution
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!  returns \sum_{Rij}<u|p_i>a_ij<p_j|u> for various a_ij inputs
!!
!! SOURCE

subroutine orbmag_nl1_k(atindx,cprj_k,dimlmn,dterm,dtset,ikpt,isppol,m1_k,mcprjk,&
    & mkmem_rbz,nband_k,nl1_option,occ_k,pawtab,ucvol)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcprjk,mkmem_rbz,nband_k,nl1_option
  real(dp),intent(in) :: ucvol
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: occ_k(nband_k)
  real(dp),intent(out) :: m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,nn
  real(dp) :: trnrm
  complex(dpc) :: tt
  !arrays
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 m1_k = zero

 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do adir = 1, 3
   do nn = 1, nband_k
     trnrm = occ_k(nn)*dtset%wtk(ikpt)/ucvol
     if(abs(trnrm).LT.tol8) cycle

     call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
       & mkmem_rbz,dtset%natom,1,nband_k,dtset%nspinor,dtset%nsppol,0)

     select case (nl1_option)
     case(1)
       call tt_me(dterm%LR(:,:,:,adir),atindx,cwaveprj,dtset,cwaveprj,dterm%lmn2max,dterm%ndij,pawtab,tt)
     case(2)
       call tt_me(dterm%BM(:,:,:,adir),atindx,cwaveprj,dtset,cwaveprj,dterm%lmn2max,dterm%ndij,pawtab,tt)
     case default
       tt = czero
     end select

     m1_k(1,nn,adir) = trnrm*real(tt); m1_k(2,nn,adir) = trnrm*aimag(tt)

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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  occ_k=band occupations at this kpt
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ucvol=unit cell volume
!!
!!
!! OUTPUT
!!  m1_k(2,nband_k,3)=Orb mag contribution
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computes -\frac{i}{2}\sum_{Rij}<u|d_b p_i>D^0_{ij} - E^0s^0_{ij}<d_g p_j|u>
!!
!! SOURCE

subroutine orbmag_nl_k(atindx,cprj_k,dimlmn,dterm,dtset,eig_k,ikpt,isppol,&
    & m1_k,mcprjk,mkmem_rbz,nband_k,occ_k,pawtab,ucvol)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcprjk,mkmem_rbz,nband_k
  real(dp),intent(in) :: ucvol
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),occ_k(nband_k)
  real(dp),intent(out) :: m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,nn
  real(dp) :: epsabg,trnrm
  complex(dpc) :: m1,prefac_m,txt_d,txt_q
  !arrays
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 m1_k = zero
 ABI_MALLOC(cwaveprj,(dtset%natom,dtset%nspinor))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 do adir = 1, 3
   do nn = 1, nband_k
     trnrm = occ_k(nn)*dtset%wtk(ikpt)/ucvol
     if(abs(trnrm).LT.tol8) cycle
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

     m1_k(1,nn,adir) = trnrm*real(m1); m1_k(2,nn,adir) = trnrm*aimag(m1)

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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  b1_k(2,nband_k,3)=Berry curvature contribution from <P_cu|Pc_u> terms
!!  m1_k(2,nband_k,3)=Orb mag contribution from <P_cu|Pc_u> terms
!!  m1_mu_k(2,nband_k,3)=Orb mag contribution from <P_cu|Pc_u> terms with offset
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! SOURCE

subroutine orbmag_cc_k(atindx,b1_k,cprj1_k,dimlmn,dtset,eig_k,fermie,gs_hamk,ikpt,isppol,&
    & m1_k,m1_mu_k,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,occ_k,pcg1_k,ucvol)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie,ucvol
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),occ_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: b1_k(2,nband_k,3),m1_k(2,nband_k,3),m1_mu_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,mcprjk,3)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,gdir,ndat,nn,npwsp,sij_opt,tim_getghc,type_calc
  real(dp) :: doti,dotr,epsabg,lams,trnrm
  complex(dpc) :: b1,m1,m1_mu,prefac_b,prefac_m
  !arrays
  real(dp),allocatable :: bra(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),ket(:,:)
  type(pawcprj_type),allocatable :: cwaveprj1(:,:)

!--------------------------------------------------------------------

 b1_k = zero
 m1_k = zero
 m1_mu_k = zero
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
     trnrm = occ_k(nn)*dtset%wtk(ikpt)/ucvol
     if (abs(trnrm).LT.tol8) cycle

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

     b1_k(1,nn,adir) = trnrm*real(b1); b1_k(2,nn,adir) = trnrm*aimag(b1)
     m1_k(1,nn,adir) = trnrm*real(m1); m1_k(2,nn,adir) = trnrm*aimag(m1)
     m1_mu_k(1,nn,adir) = trnrm*real(m1_mu); m1_mu_k(2,nn,adir) = trnrm*aimag(m1_mu)

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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!  occ_k=band occupations at this kpt
!!  pcg1_k(2,mcgk,3)=cg1_k projected on conduction space
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  b1_k(2,nband_k,3)=Berry curvature contribution from <u|Pc_u> terms
!!  b2_k(2,nband_k,3)=Berry curvature contribution from <u|d S|u> terms
!!  m1_k(2,nband_k,3)=Orb mag contribution from <u|Pc_u> terms
!!  m1_mu_k(2,nband_k,3)=Orb mag contribution from <u|Pc_u> terms with offset
!!  m2_k(2,nband_k,3)=Orb mag contribution from <u|dS|u> terms
!!  m2_mu_k(2,nband_k,3)=Orb mag contribution from <u|dS|u> terms with offset
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! contributions (1) <Pc d_b u|E d_gS|u> + <u|E d_bS|Pc d_g u> and
!! (2) \sum_n' <u |d_b ES|u_n'><u_n'|d_g ES|u> to orbital magnetization
!!
!! SOURCE

subroutine orbmag_vv_k(atindx,b1_k,b2_k,cg_k,cprj_k,dimlmn,dtset,eig_k,fermie,gs_hamk,&
    & ikpt,isppol,m1_k,m1_mu_k,m2_k,m2_mu_k,mcgk,mcprjk,mkmem_rbz,mpi_enreg,nband_k,npw_k,&
    & occ_k,pcg1_k,ucvol)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,mcprjk,mkmem_rbz,nband_k,npw_k
  real(dp),intent(in) :: fermie,ucvol
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k),occ_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: b1_k(2,nband_k,3),b2_k(2,nband_k,3)
  real(dp),intent(out) :: m1_k(2,nband_k,3),m1_mu_k(2,nband_k,3)
  real(dp),intent(out) :: m2_k(2,nband_k,3),m2_mu_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,mcprjk)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,gdir,ndat,nn,nnlout,np,npwsp
  integer :: paw_opt,signs,tim_getghc
  real(dp) :: doti,dotr,epsabg,trnrm
  complex(dpc) :: b1,bv2b,m1,m1_mu,mb,mg,mv2b,mv2b_mu,prefac_b,prefac_m
  !arrays
  real(dp) :: enlout(1),lamv(1)
  real(dp),allocatable :: bra(:,:),ket(:,:),svectoutb(:,:),svectoutg(:,:),vectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 m1_k = zero; m2_k = zero
 m1_mu_k = zero; m2_mu_k = zero
 b1_k = zero; b2_k = zero
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
     trnrm = occ_k(nn)*dtset%wtk(ikpt)/ucvol
     if(abs(trnrm).LT.tol8) cycle

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
           if(abs(occ_k(np)).LT.tol8) cycle

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

     b1_k(1,nn,adir) = trnrm*real(b1); b1_k(2,nn,adir) = trnrm*aimag(b1)
     b2_k(1,nn,adir) = trnrm*real(bv2b); b2_k(2,nn,adir) = trnrm*aimag(bv2b)
     m1_k(1,nn,adir) = trnrm*real(m1); m1_k(2,nn,adir) = trnrm*aimag(m1)
     m1_mu_k(1,nn,adir) = trnrm*real(m1_mu); m1_mu_k(2,nn,adir) = trnrm*aimag(m1_mu)
     m2_k(1,nn,adir) = trnrm*real(mv2b); m2_k(2,nn,adir) = trnrm*aimag(mv2b)
     m2_mu_k(1,nn,adir) = trnrm*real(mv2b_mu); m2_mu_k(2,nn,adir) = trnrm*aimag(mv2b_mu)

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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! SIDE EFFECTS
!!
!! TODO
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
  real(dp) :: doti,dotr,lambdar
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(dtset%natom)=index table for atoms (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  omlamb(2,3)=contribution of Lamb shielding to magnetic moment
!!
!! SIDE EFFECTS
!!
!! TODO
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
  real(dp),intent(out) :: omlamb(2,3)
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
      omlamb(1,adir) = omlamb(1,adir) - lambsig*dtset%nucdipmom(adir,iat)
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computes on-site \sum_{Rij}<u|d_bdir p_i>aij<d_gdir p_j|u> for generic aij input
!!
!! SOURCE

subroutine txt_me(aij,atindx,bcp,bdir,dtset,gdir,kcp,lmn2max,ndij,pawtab,txt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,gdir,lmn2max,ndij
  complex(dpc),intent(out) :: txt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max,ndij)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom,dtset%nspinor)
  type(pawcprj_type),intent(in) :: kcp(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,isp,jlmn,klmn
  complex(dpc) :: dcpi,dcpj,dij

  !arrays

!--------------------------------------------------------------------

  txt = czero
  do iat = 1, dtset%natom
    itypat=dtset%typat(iat)
    iatom = atindx(iat)
    do isp = 1, dtset%nspinor
      do ilmn = 1, pawtab(itypat)%lmn_size
        do jlmn = 1, pawtab(itypat)%lmn_size
          klmn = MATPACK(ilmn,jlmn)
          dij = aij(iatom,klmn,isp)
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          dcpi = CMPLX(bcp(iatom,isp)%dcp(1,bdir,ilmn),bcp(iatom,isp)%dcp(2,bdir,ilmn))
          dcpj = CMPLX(kcp(iatom,isp)%dcp(1,gdir,jlmn),kcp(iatom,isp)%dcp(2,gdir,jlmn))
          txt = txt + CONJG(dcpi)*dcpj*dij
          if (ndij .GT. 2) then
            if (isp == 1) then
              dij = aij(iatom,klmn,3) ! up-down
              if (ilmn .GT. jlmn) dij = CONJG(dij)
              dcpi = CMPLX(bcp(iatom,1)%dcp(1,bdir,ilmn),bcp(iatom,1)%dcp(2,bdir,ilmn))
              dcpj = CMPLX(kcp(iatom,2)%dcp(1,gdir,jlmn),kcp(iatom,2)%dcp(2,gdir,jlmn))
            else
              dij = aij(iatom,klmn,4) ! down-up
              if (ilmn .GT. jlmn) dij = CONJG(dij)
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! computes on-site \sum_{Rij}<u|p_i>aij<p_j|u> for generic aij input
!!
!! SOURCE

subroutine tt_me(aij,atindx,bcp,dtset,kcp,lmn2max,ndij,pawtab,tt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max,ndij
  complex(dpc),intent(out) :: tt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max,ndij)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom,dtset%nspinor),kcp(dtset%natom,dtset%nspinor)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,isp,itypat,ilmn,jlmn,klmn
  complex(dpc) :: cpi,cpj,dij

  !arrays

!--------------------------------------------------------------------

  tt = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat=dtset%typat(iat)
    do isp = 1, dtset%nspinor
      do ilmn = 1, pawtab(itypat)%lmn_size
        do jlmn = 1, pawtab(itypat)%lmn_size
          klmn=MATPACK(ilmn,jlmn)
          cpi =  CMPLX(bcp(iatom,isp)%cp(1,ilmn),bcp(iatom,isp)%cp(2,ilmn))
          cpj =  CMPLX(kcp(iatom,isp)%cp(1,jlmn),kcp(iatom,isp)%cp(2,jlmn))
          dij = aij(iatom,klmn,isp)
          ! see note at top of file near definition of MATPACK macro
          if (ilmn .GT. jlmn) dij = CONJG(dij)
          ! note use of CONJG(cpi), because cpi is from the bra side cprj
          tt = tt + CONJG(cpi)*dij*cpj
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! TODO
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! TODO
!!
!! NOTES
!! this term is A0.An = \frac{1}{2}(B x r).\alpha^2(m x r) which can be rewritten
!! as \frac{\alpha^2}{2} [(B.m)r^2 - B.rr.m]
!! the first term is bm1 below; the second term involves writing rr as a rank 2 cartesian
!! tensor, which leads to the rank 0 and rank 2 spherical components below. (that is, write
!! it as xx,xy,xz,... and write these terms in terms of their real spherical harmonic components)
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
  integer :: adir,gint,iat,iatom,ilm,ilmn,itypat,jlmn,jlm
  integer :: klmn,klm,kln,lpmp,mesh_size,pwave_size
  real(dp) :: a2,bm1,bm2,cfac,d00,d20,d22,dij,intg

  !arrays
  complex(dpc) :: dij_cart(3),dij_red(3)
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

  dterm%BM = czero
  a2 = FineStructureConstant2
  d00 = sqrt(4.0*pi)/3.0
  dij = sqrt(4.0*pi/15.0)
  d20 = sqrt(16.0*pi/5.0)/6.0
  d22 = sqrt(16.0*pi/15.0)/2.0

  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat = dtset%typat(iat)

    mesh_size=pawtab(itypat)%mesh_size
    pwave_size=size(pawtab(itypat)%phiphj(:,1))
    ABI_MALLOC(ff,(mesh_size))
    do klmn=1, pawtab(itypat)%lmn2_size
      klm  = pawtab(itypat)%indklmn(1,klmn)
      kln  = pawtab(itypat)%indklmn(2,klmn)
      ilmn = pawtab(itypat)%indklmn(7,klmn)
      ilm  = pawtab(itypat)%indlmn(4,ilmn)
      jlmn = pawtab(itypat)%indklmn(8,klmn)
      jlm  = pawtab(itypat)%indlmn(4,jlmn)

      ff=0
      ff(2:pwave_size) = &
        & (pawtab(itypat)%phiphj(2:pwave_size,kln)-pawtab(itypat)%tphitphj(2:pwave_size,kln))/&
        & (pawrad(itypat)%rad(2:pwave_size))
      call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
      call simp_gen(intg,ff,pawrad(itypat))

      do adir = 1, 3
        if (ilm .EQ. jlm) then
          bm1 = half*a2*intg*dtset%nucdipmom(adir,iat)
        else
          bm1 = zero
        end if

        bm2 = zero
        ! xx, yy, zz cases all have the same contribution from S00
        lpmp=1
        gint = gntselect(lpmp,klm)
        if (gint > 0) then
          bm2=bm2+half*a2*dtset%nucdipmom(adir,iat)*d00*realgnt(gint)*intg
        end if
        ! all other contributions involve Gaunt integrals of S_{2m}
        do lpmp = 5, 9
          gint = gntselect(lpmp,klm)
           if (gint > 0) then
             cfac = half*a2*realgnt(gint)*intg
             select case (lpmp)
             case (5) ! S_{2,-2} contributes to xy term
               select case (adir)
               case (1)
                 bm2=bm2+cfac*dtset%nucdipmom(2,iat)*dij
               case (2)
                 bm2=bm2+cfac*dtset%nucdipmom(1,iat)*dij
               end select
             case (6) ! S_{2,-1} contributes to yz term
               select case (adir)
               case (2)
                 bm2=bm2+cfac*dtset%nucdipmom(3,iat)*dij
               case (3)
                 bm2=bm2+cfac*dtset%nucdipmom(2,iat)*dij
               end select
             case (7) ! S_{2,0} contributes to xx, yy, and zz terms
               select case (adir)
               case (1)
                 bm2=bm2-cfac*dtset%nucdipmom(1,iat)*d20
               case (2)
                 bm2=bm2-cfac*dtset%nucdipmom(2,iat)*d20
               case (3)
                 bm2=bm2+cfac*dtset%nucdipmom(3,iat)*2.0*d20
               end select
             case (8) ! S_{2,+1} contributes to xz term
               select case (adir)
               case (1)
                  bm2=bm2+cfac*dtset%nucdipmom(3,iat)*dij
               case (3)
                  bm2=bm2+cfac*dtset%nucdipmom(1,iat)*dij
               end select
             case (9) ! S_{2,2} contributes to xx, yy terms
               select case (adir)
               case (1)
                 bm2=bm2+cfac*dtset%nucdipmom(1,iat)*d22
               case (2)
                 bm2=bm2-cfac*dtset%nucdipmom(2,iat)*d22
               end select
             end select
           end if ! end check on nonzero gaunt integral
        end do ! end loop over lp,mp

        dij_cart(adir) = -CMPLX(bm1 - bm2,zero)
      end do

      dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

      dterm%BM(iatom,klmn,1,1:3) = dij_red(1:3)
      if (dterm%ndij > 1) then
        dterm%BM(iatom,klmn,2,1:3) = dij_red(1:3)
      end if

    end do ! end loop over klmn
    ABI_FREE(ff)
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! TODO
!!
!! NOTES
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
  complex(dpc) :: orbl_me

  !arrays
  complex(dpc) :: dij_cart(3),dij_red(3)
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
!! COPYRIGHT
!! Copyright (C) 2003-2022 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  orbmag_terms(2,dtset%mband,dtset%nsppol,3,nterms)=all computed terms per band of orb mag
!!  orbmag_trace(2,3,nterms)=traces of orbmag_terms
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

subroutine orbmag_output(dtset,orbmag_terms,orbmag_trace)


 !Arguments ------------------------------------
 !scalars
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: orbmag_terms(2,dtset%mband,dtset%nsppol,3,nterms),orbmag_trace(2,3,nterms)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,isppol,iterms
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(2,dtset%mband,3),berry_total(2,3),orbmag_bb(2,dtset%mband,3),orbmag_total(2,3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = ibvv2+1, nterms
   orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,isppol,1:3,iterms)
     end do ! iband
   end do ! isppol
 end do
 berry_bb=zero;berry_total=zero
 do iterms = ibcc,ibvv2
   berry_total(1:2,1:3)=berry_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do isppol = 1, dtset%nsppol
     do iband=1, dtset%mband
       berry_bb(1:2,iband,1:3) = berry_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,isppol,1:3,iterms)
     end do ! iband
   end do ! isppol
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
 write(message,'(a)')' Chern vector, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (berry_total(1,adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .EQ. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(1) CC : ',(orbmag_trace(1,adir,imcc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(orbmag_trace(1,adir,imvv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(orbmag_trace(1,adir,imvv2),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     rho(0) NL : ',(orbmag_trace(1,adir,imnl),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <L_R> terms : ',(orbmag_trace(1,adir,imlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <A0.An> terms : ',(orbmag_trace(1,adir,imbm),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    Lamb terms : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Chern vector, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  Ch CC : ',(orbmag_trace(1,adir,ibcc),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv1 : ',(orbmag_trace(1,adir,ibvv1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' Ch vv2 : ',(orbmag_trace(1,adir,ibvv2),adir=1,3)
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
       write(message,'(a,3es16.8)') ' Orbital magnetic moment : ',(orbmag_bb(1,iband,adir),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '     rho(1) CC : ',(orbmag_terms(1,iband,isppol,adir,imcc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(orbmag_terms(1,iband,isppol,adir,imvv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(orbmag_terms(1,iband,isppol,adir,imvv2),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '     rho(0) NL : ',(orbmag_terms(1,iband,isppol,adir,imnl),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '   <L_R> terms : ',(orbmag_terms(1,iband,isppol,adir,imlr),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' <A0.An> terms : ',(orbmag_terms(1,iband,isppol,adir,imbm),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a)')ch10
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Chern vector : ',(berry_bb(1,iband,adir),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') '  Ch CC : ',(orbmag_terms(1,iband,isppol,adir,ibcc),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV1 : ',(orbmag_terms(1,iband,isppol,adir,ibvv1),adir=1,3)
       call wrtout(ab_out,message,'COLL')
       write(message,'(a,3es16.8)') ' Ch VV2 : ',(orbmag_terms(1,iband,isppol,adir,ibvv2),adir=1,3)
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
!! dterm <type(dterm_type)> data related to onsite interactions
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

  if(allocated(dterm%aij)) then
    ABI_FREE(dterm%aij)
  end if
  dterm%has_aij=0

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if
  dterm%has_qij=0

  if(allocated(dterm%LR)) then
    ABI_FREE(dterm%LR)
  end if
  dterm%has_LR=0

  if(allocated(dterm%BM)) then
    ABI_FREE(dterm%BM)
  end if
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

subroutine dterm_alloc(dterm,lmnmax,lmn2max,natom,ndij)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,lmn2max,natom,ndij
  type(dterm_type),intent(inout) :: dterm

  !arrays

  !Local variables -------------------------
  !scalars

  !arrays
!--------------------------------------------------------------------

  dterm%lmnmax = lmnmax
  dterm%lmn2max = lmn2max
  dterm%natom = natom
  dterm%ndij = ndij

  if(allocated(dterm%aij)) then
    ABI_FREE(dterm%aij)
  end if
  ABI_MALLOC(dterm%aij,(natom,lmn2max,ndij))
  dterm%has_aij=1

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if
  ABI_MALLOC(dterm%qij,(natom,lmn2max,ndij))
  dterm%has_qij=1

  if(allocated(dterm%LR)) then
    ABI_FREE(dterm%LR)
  end if
  ABI_MALLOC(dterm%LR,(natom,lmn2max,ndij,3))
  dterm%has_LR=1

  if(allocated(dterm%BM)) then
    ABI_FREE(dterm%BM)
  end if
  ABI_MALLOC(dterm%BM,(natom,lmn2max,ndij,3))
  dterm%has_BM=1

end subroutine dterm_alloc
!!***

!!****f* ABINIT/dterm_aij
!! NAME
!! dterm_aij
!!
!! FUNCTION
!! transfer paw_ij to dterm%aij in more convenient format
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

subroutine dterm_aij(atindx,dterm,dtset,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,idij,itypat,klmn
  !arrays
!--------------------------------------------------------------------

  dterm%aij = czero

  do iat=1,dtset%natom
    iatom=atindx(iat)
    itypat=dtset%typat(iat)
    do klmn=1,pawtab(itypat)%lmn2_size
      do idij = 1, dterm%ndij
        if (paw_ij(iatom)%cplex_dij .EQ. 2) then
          dterm%aij(iatom,klmn,idij) = &
            & CMPLX(paw_ij(iatom)%dij(2*klmn-1,idij),paw_ij(iatom)%dij(2*klmn,idij))
        else
          dterm%aij(iatom,klmn,idij) = CMPLX(paw_ij(iatom)%dij(klmn,idij),zero)
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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

subroutine make_d(atindx,dterm,dtset,gprimd,paw_ij,pawrad,pawtab,psps)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

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
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!  this routine is parallelized over kpts only, not bands
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!
!! SOURCE

subroutine local_fermie(dtset,eigen0,fermie,mpi_enreg,occ)

  !Arguments ------------------------------------
  !scalars
  real(dp),intent(out) :: fermie
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

  !Local variables -------------------------
  !scalars
  integer :: bdtot_index,ierr,ikpt,isppol,me
  integer :: nband_k,nn,nproc,spaceComm
  real(dp) :: fermie_proc

  !arrays
  real(dp),allocatable :: eig_k(:),occ_k(:)

!--------------------------------------------------------------------

  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  bdtot_index=0
  fermie_proc = -1.0D99
  do isppol = 1, dtset%nsppol
    do ikpt = 1, dtset%nkpt

      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)

      ! if the current kpt is not on the current processor, cycle
      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
        bdtot_index=bdtot_index+nband_k
        cycle
      end if

      ABI_MALLOC(occ_k,(nband_k))
      occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)

      ABI_MALLOC(eig_k,(nband_k))
      eig_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)

      do nn = 1, nband_k
        if ( (abs(occ_k(nn)).GT.tol8) .AND. (eig_k(nn).GT.fermie_proc) ) then
          fermie_proc = eig_k(nn)
        end if
      end do ! nn

      ABI_FREE(occ_k)
      ABI_FREE(eig_k)

      bdtot_index=bdtot_index+nband_k

    end do ! end loop over kpts
  end do ! end loop over isppol

  call xmpi_max(fermie_proc,fermie,spaceComm,ierr)

end subroutine local_fermie
!!***


end module m_orbmag
