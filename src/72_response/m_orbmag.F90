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
#define LMPACK(l,m) (l*l+l+1+m)

module m_orbmag

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_dfpt_nstwf,       only : gaugetransfo
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
  use m_pawfgrtab,        only : pawfgrtab_type
  use m_paw_ij,           only : paw_ij_type
  use m_paw_denpot,       only : pawdensities
  use m_pawdij,           only : pawdijhat,pawdijnd
  use m_paw_occupancies,  only : pawmkrhoij
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen,poisson
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
  integer,parameter :: ibcc=1,ibvv1=2,ibvv2=3,imcc=4,imvv1=5,imvv2=6
  integer,parameter :: imnl=7,imden1=8,imlr=9,imbm=10
  integer,parameter :: iomlmb=11
  integer,parameter :: nterms=11

  ! these parameters are constants used repeatedly
  
  ! accounts for exp(i k.r) in abinit derivatives rather than exp( 2pi i k.r)
  real(dp),parameter :: c2=one/(two_pi*two_pi) 
  real(dp),parameter :: cdij=sqrt(four_pi/three) ! conversion from r_a to cdij*S_{1a}*r
  complex(dpc),parameter :: cbc=j_dpc/two_pi ! Berry curvature pre-factor
  complex(dpc),parameter :: com=-half*j_dpc  ! Orbital magnetism pre-factor

  ! local datatype for d_\alpha terms
  type,private :: dterm_type
    ! "rdn" refers to term n in the road map paper (Torrent et al
    ! Comp Mat Sci 42 337 (2008) Appendix E)
    ! scalars
    integer :: lmnmax
    integer :: lmn2max
    integer :: natom
    integer :: has_aij=0
    integer :: has_aij1=0
    integer :: has_qij=0
    integer :: has_rd1=0
    integer :: has_rd1a=0
    integer :: has_rd2a=0
    integer :: has_rd2a1=0
    integer :: has_rd2b=0
    integer :: has_rd2c=0
    integer :: has_rd2d=0
    integer :: has_rd2e=0
    integer :: has_rd2f=0
    integer :: has_rd3a=0
    integer :: has_dijhat=0
    integer :: has_LR=0
    integer :: has_BM=0

    ! sum of \Delta A_ij
    ! aij(natom,lmn2max)
    complex(dpc),allocatable :: aij(:,:)
    
    ! sum of \Delta A(1)_ij
    ! aij1(natom,lmn2max,3) 
    complex(dpc),allocatable :: aij1(:,:,:)
    
    ! <phi|phi> - <tphi|tphi>
    ! qij(natom,lmn2max)
    complex(dpc),allocatable :: qij(:,:)
    
    ! roadmap term 1 
    ! <phi|p2/2|phi> - <tphi|p2/2|tphi>
    ! rd1(natom,lmn2max)
    complex(dpc),allocatable :: rd1(:,:)
    
    ! onsite A.p dipole term
    ! rd1a(natom,lmn2max)
    complex(dpc),allocatable :: rd1a(:,:)
    
    ! roadmap 2a
    ! <phi|vH[n1]|phi> - <tphi|vH[\tilde{n}1]|tphi>
    ! rd2a(natom,lmn2max)
    complex(dpc),allocatable :: rd2a(:,:)
    
    ! roadmap 2a(1)
    ! <phi|vH[n1(1)]|phi> - <tphi|vH[\tilde{n}1(1)]|tphi>
    ! rd2a1(natom,lmn2max,3)
    complex(dpc),allocatable :: rd2a1(:,:,:)
    
    ! roadmap 2b 
    ! <phi|vH[nZc]|phi> - <tphi|vH[nZc]|tphi>
    ! rd2b(natom,lmn2max)
    complex(dpc),allocatable :: rd2b(:,:)
    
    ! roadmap 2c 
    !  -<tphi|vH[\hat{n}]|tphi>
    ! rd2c(natom,lmn2max)
    complex(dpc),allocatable :: rd2c(:,:)
    
    ! roadmap 2d 
    !  -vH[\{n}1]Q^{LM}_{ij}
    ! rd2d(natom,lmn2max)
    complex(dpc),allocatable :: rd2d(:,:)
    
    ! roadmap 2e 
    ! vH[\tilde{n}Zc]Q^{LM}_{ij}
    ! rd2e(natom,lmn2max)
    complex(dpc),allocatable :: rd2e(:,:)
    
    ! roadmap 2f
    !  -vH[\hat{n}]Q^{LM}_{ij}
    ! rd2f(natom,lmn2max)
    complex(dpc),allocatable :: rd2f(:,:)
    
    ! roadmap 3a 
    ! <phi|vxc[n1+nc]|phi> - <tphi|vxc[\tilde{n}^1+\tilde{nc}|tphi>
    ! rd3a(natom,lmn2max)
    complex(dpc),allocatable :: rd3a(:,:)
    
    ! \hat{D}_ij
    ! dijhat(natom,lmn2max)
    complex(dpc),allocatable :: dijhat(:,:)
    
    ! onsite L_R/2
    ! <phi|L_R/2|phi> - <tphi|L_R/2|tphi>
    ! LR(natom,lmn2max,3)
    complex(dpc),allocatable :: LR(:,:,:)

    ! onsite BM
    ! <phi|Bxr . mxr|phi> - <tphi|Bxr . mxr|tphi>
    ! BM(natom,lmn2max,3)
    complex(dpc),allocatable :: BM(:,:,:)

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

  public :: orbmag

  private :: make_mbare_k
  private :: orbmag_cc_k
  private :: orbmag_vv_k
  private :: make_v2b_k
  private :: orbmag_nl_k
  private :: orbmag_nl1_k
  private :: make_m1_k
!  private :: make_alt_d
  private :: make_d
  private :: sum_d
  private :: check_eig_k
  private :: dterm_qij
  private :: dterm_rd1
  private :: dterm_rd1a
  private :: dterm_rd2a
  private :: dterm_rd2a1
  private :: dterm_rd2b
  private :: dterm_rd2c
  private :: dterm_rd2d
  private :: dterm_rd2e
  private :: dterm_rd2f
  private :: dterm_rd3a
  private :: dterm_dijhat
  private :: dterm_LR
  private :: dterm_BM
  private :: tt_me
  private :: txt_me

  private :: lamb_core
  private :: make_pcg1
  private :: pack_pawrhoij
  private :: orbmag_output
  private :: paw_sph_den_alloc
  private :: paw_sph_den_free
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

subroutine orbmag(cg,cg1,cprj,dtset,eigen0,gsqcut,kg,mcg,mcg1,mcprj,mpi_enreg,&
    & nfftf,ngfftf,npwarr,occ,paw_ij,paw_an,pawang,pawfgr,pawfgrtab,pawrad,&
    & pawtab,psps,rprimd,vtrial,vxc,xred,ylm,ylmgr)

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
 real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,buff_size,choice,cpopt,dimffnl,exchn2n3d,iat,iatom,icg,icmplx,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat,klmn,lmn2max
 integer :: me,mcgk,my_lmax,my_nspinor,nband_k,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,ngnt,nl1_option,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(dterm_type) :: dterm
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),gntselect(:,:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: b1_k(:,:,:),b2_k(:,:,:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:),cwavef_d(:,:)
 real(dp),allocatable :: eig_k(:),ffnl_k(:,:,:,:),kinpw(:),kpg_k(:,:),m1_k(:,:,:),m2_k(:,:,:)
 real(dp),allocatable :: occ_k(:),orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),realgnt(:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 logical,allocatable :: distrb_cycle(:)
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

 ABI_MALLOC(distrb_cycle,(nband_k))
 
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
 call dterm_alloc(dterm,psps%lmnmax,lmn2max,dtset%natom)

 call make_d(atindx,atindx1,cprj,dimlmn,dterm,dtset,gprimd,mcprj,nfftf,&
   & mpi_enreg,occ,paw_an,pawang,pawfgrtab,paw_ij,pawrad,pawtab,psps,&
   & ucvol,vtrial,vxc,xred)

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
   
   distrb_cycle = (mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) /= mpi_enreg%me_kpt)

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

   ! retrieve occupation numbers at this k point
   ABI_MALLOC(occ_k,(nband_k))
   occ_k(1:nband_k) = occ((ikpt-1)*nband_k+1:ikpt*nband_k)
   
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

   !! check aij against paw_ij 
   !call check_eig_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,&
   !  & gs_hamk,ikpt,isppol,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab)

   ABI_MALLOC(b1_k,(2,nband_k,3))
   ABI_MALLOC(b2_k,(2,nband_k,3))
   ABI_MALLOC(m1_k,(2,nband_k,3))
   ABI_MALLOC(m2_k,(2,nband_k,3))
   
   call orbmag_cc_k(atindx,b1_k,cprj1_k,dimlmn,dtset,eig_k,gs_hamk,ikpt,isppol,m1_k,&
    & mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab,pcg1_k)
   orbmag_terms(:,:,:,ibcc) = orbmag_terms(:,:,:,ibcc) + trnrm*b1_k(:,:,:)
   orbmag_terms(:,:,:,imcc) = orbmag_terms(:,:,:,imcc) + trnrm*m1_k(:,:,:)

   call orbmag_vv_k(atindx,b1_k,b2_k,cg_k,cprj_k,dimlmn,dtset,eig_k,gs_hamk,&
    & ikpt,isppol,m1_k,m2_k,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab,pcg1_k)
   orbmag_terms(:,:,:,ibvv1) = orbmag_terms(:,:,:,ibvv1) + trnrm*b1_k(:,:,:)
   orbmag_terms(:,:,:,ibvv2) = orbmag_terms(:,:,:,ibvv2) + trnrm*b2_k(:,:,:)
   orbmag_terms(:,:,:,imvv1) = orbmag_terms(:,:,:,imvv1) + trnrm*m1_k(:,:,:)
   orbmag_terms(:,:,:,imvv2) = orbmag_terms(:,:,:,imvv2) + trnrm*m2_k(:,:,:)
   
   call orbmag_nl_k(atindx,cprj_k,dterm,dtset,eig_k,m1_k,nband_k,paw_ij,pawtab)
   orbmag_terms(:,:,:,imnl) = orbmag_terms(:,:,:,imnl) + trnrm*m1_k(:,:,:)
   
   call make_m1_k(atindx,cprj_k,dterm,dtset,m1_k,nband_k,pawtab)
   orbmag_terms(:,:,:,imden1) = orbmag_terms(:,:,:,imden1) + trnrm*m1_k(:,:,:)

   nl1_option = 1 ! LR 
   call orbmag_nl1_k(atindx,cprj_k,dterm,dtset,m1_k,nband_k,nl1_option,pawtab)
   orbmag_terms(:,:,:,imlr) = orbmag_terms(:,:,:,imlr) + trnrm*m1_k
   nl1_option = 2 ! BM 
   call orbmag_nl1_k(atindx,cprj_k,dterm,dtset,m1_k,nband_k,nl1_option,pawtab)
   orbmag_terms(:,:,:,imbm) = orbmag_terms(:,:,:,imbm) + trnrm*m1_k
   
   ABI_FREE(b1_k)
   ABI_FREE(b2_k)
   ABI_FREE(m1_k)
   ABI_FREE(m2_k)

   icg = icg + npw_k*nband_k
   ikg = ikg + npw_k
   icprj = icprj + nband_k

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
       if((iterm.EQ.imlr).OR.(iterm.EQ.imbm)) then
         orbmag_terms(icmplx,nn,1:3,iterm) = MATMUL(rprimd,orbmag_terms(icmplx,nn,1:3,iterm))
       else
         orbmag_terms(icmplx,nn,1:3,iterm) = ucvol*MATMUL(gprimd,orbmag_terms(icmplx,nn,1:3,iterm))
       end if
     end do
   end do
 end do

 !! convert orbmag magnetization to orbital moment
 !! Berry curvature terms are ignored
 !! Lamb term ignored
 do iterm = 1, nterms
   if (iterm .EQ. ibcc) cycle
   if (iterm .EQ. ibvv1) cycle
   if (iterm .EQ. ibvv2) cycle
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
 
call orbmag_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)

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

 ABI_FREE(distrb_cycle)
 
 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

end subroutine orbmag
!!***

!!****f* ABINIT/make_m1_k
!! NAME
!! make_m1_k
!!
!! FUNCTION
!! make term in orb mag due to first order densities
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

subroutine make_m1_k(atindx,cprj_k,dterm,dtset,m1_k,nband_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,nn
  complex(dpc) :: m1,prefac
  !arrays

!--------------------------------------------------------------------

 prefac = -com
 m1_k = zero
 
 do adir = 1, 3
   do nn = 1, nband_k
     call tt_me(dterm%aij1(:,:,adir),atindx,cprj_k(:,nn),dtset,cprj_k(:,nn),&
           & dterm%lmn2max,pawtab,m1)
     m1_k(1,nn,adir) = real(prefac*m1); m1_k(2,nn,adir) = aimag(prefac*m1)
   end do !nn
 end do !adir

end subroutine make_m1_k
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

subroutine orbmag_nl1_k(atindx,cprj_k,dterm,dtset,m1_k,nband_k,nl1_option,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k,nl1_option
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,nn
  complex(dpc) :: tt
  !arrays
  complex(dpc) :: m1(3)

!--------------------------------------------------------------------

 m1_k = zero
 
 do adir = 1, 3
   do nn = 1, nband_k

     select case (nl1_option)
     case(1)
       call tt_me(dterm%LR(:,:,adir),atindx,cprj_k(:,nn),dtset,cprj_k(:,nn),&
         & dterm%lmn2max,pawtab,tt)
     case(2)
       call tt_me(dterm%BM(:,:,adir),atindx,cprj_k(:,nn),dtset,cprj_k(:,nn),&
         & dterm%lmn2max,pawtab,tt)
     case default
       tt = czero
     end select
         
     m1_k(1,nn,adir) = real(tt); m1_k(2,nn,adir) = aimag(tt)

   end do !nn
 end do !adir

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

subroutine orbmag_nl_k(atindx,cprj_k,dterm,dtset,eig_k,m1_k,nband_k,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k)
  real(dp),intent(out) :: m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,nn
  real(dp) :: epsabg
  complex(dpc) :: m1,prefac_m,txt_d,txt_q
  !arrays

!--------------------------------------------------------------------

 m1_k = zero
 
 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero
     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_m = -com*c2*epsabg
        
         call txt_me(dterm%aij,atindx,cprj_k(:,nn),bdir,dtset,&
           & gdir,cprj_k(:,nn),dterm%lmnmax,dterm%lmn2max,pawtab,txt_d)
         
         call txt_me(dterm%qij,atindx,cprj_k(:,nn),bdir,dtset,&
           & gdir,cprj_k(:,nn),dterm%lmnmax,dterm%lmn2max,pawtab,txt_q)

         m1 = m1 + prefac_m*(txt_d - eig_k(nn)*txt_q)

       end do !gdir
     end do !bdir

     m1_k(1,nn,adir) = real(m1); m1_k(2,nn,adir) = aimag(m1)

   end do !nn
 end do !adir

end subroutine orbmag_nl_k
!!***

!!****f* ABINIT/make_v2b_k
!! NAME
!! make_v2b_k
!!
!! FUNCTION
!! make v2b at k
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

subroutine make_v2b_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,gs_hamk,&
    & ikpt,isppol,mcgk,mpi_enreg,mv2b_k,my_nspinor,nband_k,npw_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k)
  real(dp),intent(out) :: mv2b_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,gdir,ndat,nn,nnlout,np
  integer :: paw_opt,sij_opt,signs,tim_nonlop
  real(dp) :: doti,dotr,epsabg
  complex(dpc) :: mb,mg,mv2b,prefac_m
  !arrays
  real(dp) :: enlout(1),lamv(1)
  real(dp),allocatable :: bra(:,:),ket(:,:),svectoutb(:,:)
  real(dp),allocatable :: svectoutg(:,:),vectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 mv2b_k = zero
 
 ABI_MALLOC(bra,(2,npw_k))
 ABI_MALLOC(ket,(2,npw_k))
 ABI_MALLOC(vectout,(2,npw_k))
 ABI_MALLOC(svectoutg,(2,npw_k))
 ABI_MALLOC(svectoutb,(2,npw_k))
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 tim_nonlop = 0
 lamv = zero
 ndat = 1
 choice = 5
 paw_opt = 3
 signs = 2 
 nnlout = 1

 do adir = 1, 3
   do nn = 1, nband_k

     mv2b = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_m = com*c2*epsabg
         
         cpopt = 4
         ket(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)
         call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,gdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutg,tim_nonlop,ket,vectout)
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,bdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutb,tim_nonlop,ket,vectout)

         do np = 1, nband_k

           bra(1:2,1:npw_k) = cg_k(1:2,(np-1)*npw_k+1:np*npw_k)
       
           dotr = DOT_PRODUCT(bra(1,:),svectoutg(1,:))+DOT_PRODUCT(bra(2,:),svectoutg(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutg(2,:))-DOT_PRODUCT(bra(2,:),svectoutg(1,:))
           mg = CMPLX(dotr,doti)

           dotr = DOT_PRODUCT(bra(1,:),svectoutb(1,:))+DOT_PRODUCT(bra(2,:),svectoutb(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutb(2,:))-DOT_PRODUCT(bra(2,:),svectoutb(1,:))
           mb = CMPLX(dotr,doti)

           mv2b = mv2b - prefac_m*CONJG(mb)*mg*eig_k(nn)
         end do ! np
       end do !gdir
     end do !bdir

     mv2b_k(1,nn,adir) = real(mv2b); mv2b_k(2,nn,adir) = aimag(mv2b)

   end do !nn
 end do !adir

 ABI_FREE(bra)
 ABI_FREE(ket)
 ABI_FREE(vectout)
 ABI_FREE(svectoutg)
 ABI_FREE(svectoutb)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 
end subroutine make_v2b_k
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

subroutine orbmag_cc_k(atindx,b1_k,cprj1_k,dimlmn,dtset,eig_k,gs_hamk,ikpt,isppol,m1_k,&
    & mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: eig_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: b1_k(2,nband_k,3),m1_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,cpopt,gdir,ndat,nn,sij_opt,tim_getghc,type_calc
  real(dp) :: doti,dotr,epsabg,lams
  complex(dpc) :: b1,m1,prefac_b,prefac_m
  !arrays
  real(dp),allocatable :: bra(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),ket(:,:)
  type(pawcprj_type),allocatable :: cwaveprj1(:,:)

!--------------------------------------------------------------------

 b1_k = zero
 m1_k = zero
 
 ABI_MALLOC(bra,(2,npw_k))
 ABI_MALLOC(ket,(2,npw_k))
 ABI_MALLOC(ghc,(2,npw_k))
 ABI_MALLOC(gsc,(2,npw_k))
 ABI_MALLOC(gvnlxc,(2,npw_k))
 ABI_MALLOC(cwaveprj1,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj1,0,dimlmn)

 tim_getghc = 0
 lams = zero
 ndat = 1
 sij_opt = 1
 type_calc = 0

 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero
     b1 = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg
         
         cpopt = 2
         ket(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
         call pawcprj_get(atindx,cwaveprj1,cprj1_k(:,:,gdir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
         call getghc(cpopt,ket,cwaveprj1,ghc,gsc,gs_hamk,gvnlxc,lams,mpi_enreg,&
           & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)
         
         bra(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
       
         dotr = DOT_PRODUCT(bra(1,:),ghc(1,:))+DOT_PRODUCT(bra(2,:),ghc(2,:))
         doti = DOT_PRODUCT(bra(1,:),ghc(2,:))-DOT_PRODUCT(bra(2,:),ghc(1,:))
         m1 = m1 + prefac_m*CMPLX(dotr,doti)
         
         dotr = DOT_PRODUCT(bra(1,:),gsc(1,:))+DOT_PRODUCT(bra(2,:),gsc(2,:))
         doti = DOT_PRODUCT(bra(1,:),gsc(2,:))-DOT_PRODUCT(bra(2,:),gsc(1,:))
         b1 = b1 -two*prefac_b*CMPLX(dotr,doti)
         m1 = m1 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)
         
       end do !gdir
     end do !bdir

     b1_k(1,nn,adir) = real(b1); b1_k(2,nn,adir) = aimag(b1)
     m1_k(1,nn,adir) = real(m1); m1_k(2,nn,adir) = aimag(m1)

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

subroutine orbmag_vv_k(atindx,b1_k,b2_k,cg_k,cprj_k,dimlmn,dtset,eig_k,gs_hamk,&
    & ikpt,isppol,m1_k,m2_k,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: b1_k(2,nband_k,3),b2_k(2,nband_k,3)
  real(dp),intent(out) :: m1_k(2,nband_k,3),m2_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,gdir,ndat,nn,nnlout,np
  integer :: paw_opt,signs,tim_getghc
  real(dp) :: doti,dotr,epsabg
  complex(dpc) :: b1,bv2b,m1,mb,mg,mv2b,prefac_b,prefac_m
  !arrays
  real(dp) :: enlout(1),lamv(1)
  real(dp),allocatable :: bra(:,:),ket(:,:),svectoutb(:,:),svectoutg(:,:),vectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------

 m1_k = zero; m2_k = zero
 b1_k = zero; b2_k = zero
 
 ABI_MALLOC(bra,(2,npw_k))
 ABI_MALLOC(ket,(2,npw_k))
 ABI_MALLOC(svectoutb,(2,npw_k))
 ABI_MALLOC(svectoutg,(2,npw_k))
 ABI_MALLOC(vectout,(2,npw_k))
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 tim_getghc = 0
 lamv = zero
 ndat = 1
 cpopt = 4
 choice = 5
 paw_opt = 3
 signs = 2 
 nnlout = 1

 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero; mv2b = czero
     b1 = czero; bv2b = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_b = cbc*c2*epsabg
         prefac_m = com*c2*epsabg
        
         ! extract |u_nk> 
         ket(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)
         call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
      
         ! compute dS/dk_g |u_nk> 
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,gdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutg,tim_getghc,ket,vectout)
         ! compute dS/dk_b|u_nk>
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,bdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectoutb,tim_getghc,ket,vectout)
         
         ! extract |Pc du/dk_b>         
         bra(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
         !overlap
         dotr = DOT_PRODUCT(bra(1,:),svectoutg(1,:))+DOT_PRODUCT(bra(2,:),svectoutg(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectoutg(2,:))-DOT_PRODUCT(bra(2,:),svectoutg(1,:))
         ! add <Pc du/dk_b|dS/dk_g|u_nk>*E_nk
         b1 = b1 - prefac_b*CMPLX(dotr,doti)
         m1 = m1 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)

         ! extract |Pc du/dk_g>         
         bra(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
         !overlap
         dotr = DOT_PRODUCT(bra(1,:),svectoutb(1,:))+DOT_PRODUCT(bra(2,:),svectoutb(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectoutb(2,:))-DOT_PRODUCT(bra(2,:),svectoutb(1,:))
         ! add CONJG(<Pc du/dk_b|dS/dk_g|u_nk>)*E_nk
         b1 = b1 - prefac_b*CMPLX(dotr,-doti)
         m1 = m1 + prefac_m*CMPLX(dotr,-doti)*eig_k(nn)
         
         do np = 1, nband_k

           bra(1:2,1:npw_k) = cg_k(1:2,(np-1)*npw_k+1:np*npw_k)
       
           dotr = DOT_PRODUCT(bra(1,:),svectoutg(1,:))+DOT_PRODUCT(bra(2,:),svectoutg(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutg(2,:))-DOT_PRODUCT(bra(2,:),svectoutg(1,:))
           mg = CMPLX(dotr,doti)

           dotr = DOT_PRODUCT(bra(1,:),svectoutb(1,:))+DOT_PRODUCT(bra(2,:),svectoutb(2,:))
           doti = DOT_PRODUCT(bra(1,:),svectoutb(2,:))-DOT_PRODUCT(bra(2,:),svectoutb(1,:))
           mb = CMPLX(dotr,doti)

           bv2b = bv2b + prefac_b*CONJG(mb)*mg
           mv2b = mv2b - prefac_m*CONJG(mb)*mg*eig_k(nn)
         end do ! np
 
       end do !gdir
     end do !bdir

     b1_k(1,nn,adir) = real(b1); b1_k(2,nn,adir) = aimag(b1)
     b2_k(1,nn,adir) = real(bv2b); b2_k(2,nn,adir) = aimag(bv2b)
     m1_k(1,nn,adir) = real(m1); m1_k(2,nn,adir) = aimag(m1)
     m2_k(1,nn,adir) = real(mv2b); m2_k(2,nn,adir) = aimag(mv2b)

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




!!****f* ABINIT/make_mbare_k
!! NAME
!! make_mbare_k
!!
!! FUNCTION
!! make mbare at k
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

subroutine make_mbare_k(atindx,b1_k,b2_k,cg_k,cprj_k,cprj1_k,dimlmn,dterm,dtset,eig_k,gs_hamk,&
    & ikpt,isppol,m1_k,m2_k,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab,pcg1_k)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dterm_type),intent(in) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k),pcg1_k(2,mcgk,3)
  real(dp),intent(out) :: b1_k(2,nband_k,3),b2_k(2,nband_k,3)
  real(dp),intent(out) :: m1_k(2,nband_k,3),m2_k(2,nband_k,3)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawcprj_type),intent(in) :: cprj1_k(dtset%natom,nband_k,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,choice,cpopt,gdir,ndat,nn,nnlout,np
  integer :: paw_opt,sij_opt,signs,tim_getghc,type_calc
  real(dp) :: doti,dotr,epsabg,lams
  complex(dpc) :: b1,b2,m1,m2,prefac_b,prefac_m
  !arrays
  real(dp) :: enlout(1),lamv(1)
  real(dp),allocatable :: bra(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:),ket(:,:),svectout(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:),cwaveprj1(:,:)

!--------------------------------------------------------------------

 b1_k = zero; b2_k = zero
 m1_k = zero; m2_k = zero
 
 ABI_MALLOC(bra,(2,npw_k))
 ABI_MALLOC(ket,(2,npw_k))
 ABI_MALLOC(ghc,(2,npw_k))
 ABI_MALLOC(gsc,(2,npw_k))
 ABI_MALLOC(svectout,(2,npw_k))
 ABI_MALLOC(gvnlxc,(2,npw_k))
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)
 ABI_MALLOC(cwaveprj1,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj1,0,dimlmn)

 tim_getghc = 0
 lams = zero; lamv = zero
 ndat = 1
 sij_opt = 1
 type_calc = 0
 choice = 5
 paw_opt = 3
 signs = 2 
 nnlout = 1

 do adir = 1, 3
   do nn = 1, nband_k

     m1 = czero; b1 = czero
     m2 = czero; b2 = czero

     do bdir = 1, 3
       do gdir = 1, 3
         epsabg = eijk(adir,bdir,gdir)
         if (ABS(epsabg) .LT. half) cycle
         prefac_m = com*c2*epsabg
         prefac_b = cbc*c2*epsabg
         
         cpopt = 2
         ket(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
         call pawcprj_get(atindx,cwaveprj1,cprj1_k(:,:,gdir),dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
         call getghc(cpopt,ket,cwaveprj1,ghc,gsc,gs_hamk,gvnlxc,lams,mpi_enreg,&
           & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)
         
         cpopt = 4
         ket(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)
         call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
           & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)

         bra(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
       
         dotr = DOT_PRODUCT(bra(1,:),ghc(1,:))+DOT_PRODUCT(bra(2,:),ghc(2,:))
         doti = DOT_PRODUCT(bra(1,:),ghc(2,:))-DOT_PRODUCT(bra(2,:),ghc(1,:))
         m1 = m1 + prefac_m*CMPLX(dotr,doti)
         
         dotr = DOT_PRODUCT(bra(1,:),gsc(1,:))+DOT_PRODUCT(bra(2,:),gsc(2,:))
         doti = DOT_PRODUCT(bra(1,:),gsc(2,:))-DOT_PRODUCT(bra(2,:),gsc(1,:))
         m2 = m2 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)
         b1 = b1 + prefac_b*CMPLX(dotr,doti)
         
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,gdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectout,tim_getghc,ket,gvnlxc)
         dotr = DOT_PRODUCT(bra(1,:),svectout(1,:))+DOT_PRODUCT(bra(2,:),svectout(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectout(2,:))-DOT_PRODUCT(bra(2,:),svectout(1,:))
         m2 = m2 + prefac_m*CMPLX(dotr,doti)*eig_k(nn)
         b2 = b2 + prefac_b*CMPLX(dotr,doti)

         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,bdir,lamv,mpi_enreg,ndat,nnlout,&
           & paw_opt,signs,svectout,tim_getghc,ket,gvnlxc)
         bra(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
         dotr = DOT_PRODUCT(bra(1,:),svectout(1,:))+DOT_PRODUCT(bra(2,:),svectout(2,:))
         doti = DOT_PRODUCT(bra(1,:),svectout(2,:))-DOT_PRODUCT(bra(2,:),svectout(1,:))
         m2 = m2 + prefac_m*CMPLX(dotr,-doti)*eig_k(nn)
         
       end do !gdir
     end do !bdir

     m1_k(1,nn,adir) = real(m1); m1_k(2,nn,adir) = aimag(m1)
     b1_k(1,nn,adir) = real(b1); b1_k(2,nn,adir) = aimag(b1)
     m2_k(1,nn,adir) = real(m2); m2_k(2,nn,adir) = aimag(m2)
     b2_k(1,nn,adir) = real(b2); b2_k(2,nn,adir) = aimag(b2)

   end do !nn
 end do !adir

 ABI_FREE(bra)
 ABI_FREE(ket)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(svectout)
 ABI_FREE(gvnlxc)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 call pawcprj_free(cwaveprj1)
 ABI_FREE(cwaveprj1)
 
end subroutine make_mbare_k
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

!!****f* ABINIT/txt_me
!! NAME
!! txt_me
!!
!! FUNCTION
!! Matrix element <u_n|dp>A<dp|u_m>
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
!! kcp :: ket side cprj
!! gdir :: ket side deriv direc
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

subroutine txt_me(aij,atindx,bcp,bdir,dtset,gdir,kcp,lmnmax,lmn2max,pawtab,txt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,gdir,lmnmax,lmn2max
  complex(dpc),intent(out) :: txt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom),kcp(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,jlmn,klmn
  complex(dpc) :: dcpi,dcpj,dij

  !arrays

!--------------------------------------------------------------------

  txt = czero
  do iat = 1, dtset%natom
    itypat=dtset%typat(iat)
    iatom = atindx(iat)
    do ilmn = 1, pawtab(itypat)%lmn_size
      do jlmn = 1, pawtab(itypat)%lmn_size
        klmn = MATPACK(ilmn,jlmn)
        dij = aij(iatom,klmn)
        ! see note at top of file near definition of MATPACK macro
        if (ilmn .GT. jlmn) dij = CONJG(dij)
        dcpi = CMPLX(bcp(iatom)%dcp(1,bdir,ilmn),bcp(iatom)%dcp(2,bdir,ilmn))
        dcpj = CMPLX(kcp(iatom)%dcp(1,gdir,jlmn),kcp(iatom)%dcp(2,gdir,jlmn))
        txt = txt + CONJG(dcpi)*dcpj*dij
      end do !jlmn
    end do !ilmn
  end do !iat

end subroutine txt_me
!!***


!!****f* ABINIT/tt_me
!! NAME
!! tt_me
!!
!! FUNCTION
!! Matrix element <u_n|T^\dag A T|u_m>
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

subroutine tt_me(aij,atindx,bcp,dtset,kcp,lmn2max,pawtab,tt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmn2max
  complex(dpc),intent(out) :: tt
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  complex(dpc),intent(in) :: aij(dtset%natom,lmn2max)
  type(pawcprj_type),intent(in) :: bcp(dtset%natom),kcp(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilmn,jlmn,klmn
  complex(dpc) :: cpi,cpj,dij

  !arrays

!--------------------------------------------------------------------

  tt = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat=dtset%typat(iat)
    do ilmn = 1, pawtab(itypat)%lmn_size
      do jlmn = 1, pawtab(itypat)%lmn_size
        klmn=MATPACK(ilmn,jlmn)
        cpi =  CMPLX(bcp(iatom)%cp(1,ilmn),bcp(iatom)%cp(2,ilmn))
        cpj =  CMPLX(kcp(iatom)%cp(1,jlmn),kcp(iatom)%cp(2,jlmn))
        dij = aij(iatom,klmn)
        ! see note at top of file near definition of MATPACK macro
        if (ilmn .GT. jlmn) dij = CONJG(dij)
        tt = tt + CONJG(cpi)*dij*cpj
      end do !jlmn
    end do !ilmn
  end do !iat

end subroutine tt_me
!!***

!!****f* ABINIT/dterm_qij
!! NAME
!! dterm_qij
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

subroutine dterm_qij(atindx,dterm,dtset,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,ilm,ilmn,jlm,jlmn
  integer :: klm,klmn,kln,mesh_size
  real(dp) :: qij

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

 dterm%qij = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)

   mesh_size = pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))

   do klmn = 1, pawtab(itypat)%lmn2_size
     klm  = pawtab(itypat)%indklmn(1,klmn)
     kln  = pawtab(itypat)%indklmn(2,klmn)
     ilmn = pawtab(itypat)%indklmn(7,klmn) 
     ilm  = pawtab(itypat)%indlmn(4,ilmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn) 
     jlm  = pawtab(itypat)%indlmn(4,jlmn)

     qij = zero
     if (ilm .EQ. jlm) then
       ff(2:mesh_size) = pawtab(itypat)%phiphj(2:mesh_size,kln)
       ff(2:mesh_size) = ff(2:mesh_size) - &
         & pawtab(itypat)%tphitphj(2:mesh_size,kln)
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(qij,ff,pawrad(itypat))
     end if
     dterm%qij(iatom,klmn) = CMPLX(qij,zero)

   end do ! klmn

   ABI_FREE(ff)

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

      dterm%BM(iatom,klmn,1:3) = dij_red(1:3)
      
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
          dterm%LR(iatom,klmn,1:3) = dij_red(1:3)
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
!! This routine outputs orbmag terms 
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

subroutine orbmag_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)


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
 do iterms = ibvv2+1, nterms
   orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do iband=1, nband_k
     orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
   end do
 end do
 berry_bb=zero;berry_total=zero
 do iterms = ibcc,ibvv2
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
 write(message,'(a)')' Chern vector, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (berry_total(1,adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .GE. 2) then
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
   write(message,'(a,3es16.8)') '   rho(0) n(1) : ',(orbmag_trace(1,adir,imden1),adir=1,3)
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

 if(abs(dtset%orbmag) .GE. 3) then
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
     write(message,'(a,3es16.8)') '     rho(1) CC : ',(orbmag_terms(1,iband,adir,imcc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    rho(1) VV1 : ',(orbmag_terms(1,iband,adir,imvv1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    rho(1) VV2 : ',(orbmag_terms(1,iband,adir,imvv2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '     rho(0) NL : ',(orbmag_terms(1,iband,adir,imnl),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   rho(0) n(1) : ',(orbmag_terms(1,iband,adir,imden1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '   <L_R> terms : ',(orbmag_terms(1,iband,adir,imlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' <A0.An> terms : ',(orbmag_terms(1,iband,adir,imbm),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    Lamb terms : ',(orbmag_terms(1,iband,adir,iomlmb),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Chern vector : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  Ch CC : ',(orbmag_terms(1,iband,adir,ibcc),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Ch VV1 : ',(orbmag_terms(1,iband,adir,ibvv1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') ' Ch VV2 : ',(orbmag_terms(1,iband,adir,ibvv2),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_output
!!***

!!****f* ABINIT/dterm_rd2e
!! NAME
!! dterm_rd2e
!!
!! FUNCTION
!! Compute vH[\tilde{n}_{Zc}]Q^{LM}_{ij}
!! roadmap term 2e
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

subroutine dterm_rd2e(atindx,dterm,dtset,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,ilmn,imesh,itypat,jlmn,klmn,mesh_size
  real(dp) :: intff,rr

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

 dterm%rd2e = czero

 do itypat = 1, dtset%ntypat

   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))

   do klmn = 1, pawtab(itypat)%lmn2_size
     ilmn = pawtab(itypat)%indklmn(7,klmn)
     jlmn = pawtab(itypat)%indklmn(8,klmn)

     do imesh = 2, mesh_size
       rr = pawrad(itypat)%rad(imesh)
       ff(imesh) = (rr**2)*pawtab(itypat)%vhtnzc(imesh)*pawtab(itypat)%shapefunc(imesh,1)
     end do ! imesh

     call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
     call simp_gen(intff,ff,pawrad(itypat))

     do iat = 1, dtset%natom
       iatom=atindx(iat)
       if (dtset%typat(iat) .EQ. itypat) then
         dterm%rd2e(iatom,klmn) = CMPLX(-sqrt(four_pi)*pawtab(itypat)%qijl(1,klmn)*intff,zero)
       end if
     end do !iat
   end do ! klmn

   ABI_FREE(ff)

 end do ! itypat

 dterm%has_rd2e = 2
 
end subroutine dterm_rd2e


!!****f* ABINIT/dterm_rd2b
!! NAME
!! dterm_vhnzc
!!
!! FUNCTION
!! Compute onsite vhnzc
!! roadmap term 2b
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

subroutine dterm_rd2b(atindx,dterm,dtset,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,ilm,ilmn,itypat,jlm,jlmn
  integer :: klm,klmn,kln,mesh_size,pwave_size
  real(dp) :: intff

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

 dterm%rd2b = czero

 do itypat = 1, dtset%ntypat

   mesh_size=pawtab(itypat)%mesh_size
   pwave_size=size(pawtab(itypat)%phiphj(:,1))
   ABI_MALLOC(ff,(mesh_size))

   do ilmn = 1, pawtab(itypat)%lmn_size
     do jlmn = 1, pawtab(itypat)%lmn_size
       klmn = MATPACK(ilmn,jlmn)

       klm = pawtab(itypat)%indklmn(1,klmn)
       kln = pawtab(itypat)%indklmn(2,klmn)
       ilm = pawtab(itypat)%indklmn(5,klmn)
       jlm = pawtab(itypat)%indklmn(6,klmn)

       ff(2:pwave_size) = pawtab(itypat)%phiphj(2:pwave_size,kln)*&
         & pawtab(itypat)%vhnzc(2:pwave_size)
       ff(2:pwave_size) = ff(2:pwave_size) - &
         & pawtab(itypat)%tphitphj(2:pwave_size,kln)*&
         & pawtab(itypat)%vhtnzc(2:pwave_size)

       intff = zero    
       if (ilm == jlm) then 
         call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
         call simp_gen(intff,ff,pawrad(itypat))
       end if

       do iat = 1, dtset%natom
         iatom=atindx(iat)
         if (dtset%typat(iat) .EQ. itypat) then
           dterm%rd2b(iatom,klmn) = intff
         end if
       end do !iat
     end do !jlmn
   end do ! ilmn
   ABI_FREE(ff)

 end do ! iat

 dterm%has_rd2b = 2
 
end subroutine dterm_rd2b

!!****f* ABINIT/dterm_rd2f
!! NAME
!! dterm_rd2f
!!
!! FUNCTION
!! Compute -vH[\hat{n}]Q^{LM}_ij, term 2f of roadmap paper
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

subroutine dterm_rd2f(atindx,dterm,dtset,gntselect,lmnmax,my_lmax,&
    & pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(dterm_type),intent(inout) :: dterm

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,igij,igkl,ijlm,ijlmax,ijlmin,ijlmn,ijln,itypat
  integer :: kllm,kllmax,kllmin,kllmn,klln,ll,llmm,mesh_size,meshsz
  real(dp) :: dijterm,eijkl,qij,qkl,vh1

  !arrays
  real(dp),allocatable :: ff(:),nlff(:)

!--------------------------------------------------------------------

 dterm%rd2f = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(nlff,(mesh_size))
 
   meshsz=pawrad(itypat)%int_meshsz
   if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero

   do ijlmn = 1, pawtab(itypat)%lmn2_size
     ijlm   = pawtab(itypat)%indklmn(1,ijlmn)
     ijln   = pawtab(itypat)%indklmn(2,ijlmn)
     ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
     ijlmax = pawtab(itypat)%indklmn(4,ijlmn)

     dijterm=zero
     do kllmn = 1, pawtab(itypat)%lmn2_size
       kllm   = pawtab(itypat)%indklmn(1,kllmn)
       klln   = pawtab(itypat)%indklmn(2,kllmn)
       kllmin = pawtab(itypat)%indklmn(3,kllmn)
       kllmax = pawtab(itypat)%indklmn(4,kllmn)

       eijkl=zero
       do ll = max(ijlmin,kllmin),min(ijlmax,kllmax)
         do llmm = ll**2+1,(ll+1)**2
           igij = gntselect(llmm,ijlm)
           igkl = gntselect(llmm,kllm)
           if ( (igij .EQ. 0) .OR. (igkl .EQ. 0) ) cycle

           qij =pawtab(itypat)%qijl(llmm,ijlmn)
           qkl =pawtab(itypat)%qijl(llmm,kllmn)
  
           ff(2:mesh_size)=pawtab(itypat)%shapefunc(2:mesh_size,ll+1)*&
             &pawrad(itypat)%rad(2:mesh_size)**2
           call poisson(ff,ll,pawrad(itypat),nlff)
           
           ff(1)=zero
           ff(2:meshsz)= -pawtab(itypat)%shapefunc(2:mesh_size,ll+1)*nlff(2:meshsz)&
             & *pawrad(itypat)%rad(2:meshsz)
           call simp_gen(vh1,ff,pawrad(itypat))

           eijkl = eijkl + four_pi*vh1*qij*qkl
         end do !llmm
       end do !ll

       dijterm = dijterm + pawrhoij(iatom)%rhoij_(2*kllmn-1,1)*&
         & eijkl*pawtab(itypat)%dltij(kllmn)
     end do !kllmn
     
     dterm%rd2f(iatom,ijlmn) = CMPLX(dijterm,zero)
   end do !ijlmn

   ABI_FREE(ff)
   ABI_FREE(nlff)

 end do ! end loop over iat

 dterm%has_rd2f = 2

end subroutine dterm_rd2f

!!****f* ABINIT/dterm_rd2d
!! NAME
!! dterm_rd2d
!!
!! FUNCTION
!! Compute -vH[\tilde{n}^1]Q^{LM}_ij, term 2d of roadmap paper
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

subroutine dterm_rd2d(atindx,dterm,dtset,gntselect,lmnmax,my_lmax,&
    & pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(dterm_type),intent(inout) :: dterm

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,igij,igkl,ijlm,ijlmax,ijlmin,ijlmn,ijln,itypat
  integer :: kllm,kllmax,kllmin,kllmn,klln,ll,llmm,mesh_size,meshsz
  real(dp) :: dijterm,eijkl,rgkl,qij,vh1

  !arrays
  real(dp),allocatable :: ff(:),nlff(:)

!--------------------------------------------------------------------

 dterm%rd2d = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(nlff,(mesh_size))
 
   meshsz=pawrad(itypat)%int_meshsz
   if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero

   do ijlmn = 1, pawtab(itypat)%lmn2_size
     ijlm   = pawtab(itypat)%indklmn(1,ijlmn)
     ijln   = pawtab(itypat)%indklmn(2,ijlmn)
     ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
     ijlmax = pawtab(itypat)%indklmn(4,ijlmn)

     dijterm=zero
     do kllmn = 1, pawtab(itypat)%lmn2_size
       kllm   = pawtab(itypat)%indklmn(1,kllmn)
       klln   = pawtab(itypat)%indklmn(2,kllmn)
       kllmin = pawtab(itypat)%indklmn(3,kllmn)
       kllmax = pawtab(itypat)%indklmn(4,kllmn)

       eijkl = zero
       do ll = max(ijlmin,kllmin),min(ijlmax,kllmax)
         do llmm = ll**2+1,(ll+1)**2
           igij = gntselect(llmm,ijlm)
           igkl = gntselect(llmm,kllm)
           if ( (igij .EQ. 0) .OR. (igkl .EQ. 0) ) cycle

           rgkl=realgnt(igkl)
           qij =pawtab(itypat)%qijl(llmm,ijlmn)
  
           ff(2:meshsz) = pawtab(itypat)%tphitphj(2:meshsz,klln)
           call poisson(ff,ll,pawrad(itypat),nlff)
           
           ff(1)=zero
           ff(2:mesh_size)=-pawtab(itypat)%shapefunc(2:mesh_size,ll+1)*&
             &nlff(2:mesh_size)*pawrad(itypat)%rad(2:mesh_size)
           call simp_gen(vh1,ff,pawrad(itypat))

           eijkl = eijkl + four_pi*vh1*qij*rgkl
         end do !llmm
       end do !ll

       dijterm = dijterm + pawrhoij(iatom)%rhoij_(2*kllmn-1,1)*&
         & eijkl*pawtab(itypat)%dltij(kllmn)
     end do !kllmn

     dterm%rd2d(iatom,ijlmn) = CMPLX(dijterm,zero)
   end do !ijlmn

   ABI_FREE(ff)
   ABI_FREE(nlff)

 end do ! end loop over iat

 dterm%has_rd2d = 2

end subroutine dterm_rd2d

!!****f* ABINIT/dterm_rd2c
!! NAME
!! dterm_rd2c
!!
!! FUNCTION
!! Compute -<tphi|vH[nhat]|tphj>, term 2c of roadmap paper
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

subroutine dterm_rd2c(atindx,dterm,dtset,gntselect,lmnmax,my_lmax,&
    & pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dataset_type),intent(in) :: dtset
  type(dterm_type),intent(inout) :: dterm

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,igij,igkl,ijlm,ijlmax,ijlmin,ijlmn,ijln,itypat
  integer :: kllm,kllmax,kllmin,kllmn,klln,ll,llmm,mesh_size,meshsz
  real(dp) :: dijterm,eijkl,rgij,qkl,vh1

  !arrays
  real(dp),allocatable :: ff(:),nlff(:)

!--------------------------------------------------------------------

 dterm%rd2c = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(nlff,(mesh_size))
 
   meshsz=pawrad(itypat)%int_meshsz
   if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero

   do ijlmn = 1, pawtab(itypat)%lmn2_size
     ijlm   = pawtab(itypat)%indklmn(1,ijlmn)
     ijln   = pawtab(itypat)%indklmn(2,ijlmn)
     ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
     ijlmax = pawtab(itypat)%indklmn(4,ijlmn)

     dijterm = zero
     do kllmn = 1, pawtab(itypat)%lmn2_size
       kllm   = pawtab(itypat)%indklmn(1,kllmn)
       klln   = pawtab(itypat)%indklmn(2,kllmn)
       kllmin = pawtab(itypat)%indklmn(3,kllmn)
       kllmax = pawtab(itypat)%indklmn(4,kllmn)

       eijkl = zero
       do ll = max(ijlmin,kllmin),min(ijlmax,kllmax)
         do llmm = ll**2+1,(ll+1)**2
           igij = gntselect(llmm,ijlm)
           igkl = gntselect(llmm,kllm)
           if ( (igij .EQ. 0) .OR. (igkl .EQ. 0) ) cycle

           rgij=realgnt(igij)
           qkl =pawtab(itypat)%qijl(llmm,kllmn)
  
           ff(2:mesh_size)=pawtab(itypat)%shapefunc(2:mesh_size,ll+1)*&
             &pawrad(itypat)%rad(2:mesh_size)**2
           call poisson(ff,ll,pawrad(itypat),nlff)
           
           ff(1)=zero
           ff(2:meshsz)= -pawtab(itypat)%tphitphj(2:meshsz,ijln)*nlff(2:meshsz)&
             & /pawrad(itypat)%rad(2:meshsz)
           call simp_gen(vh1,ff,pawrad(itypat))

           eijkl = eijkl + four_pi*vh1*rgij*qkl
         end do !llmm
       end do !ll
       
       dijterm = dijterm + pawrhoij(iatom)%rhoij_(2*kllmn-1,1)*&
         & eijkl*pawtab(itypat)%dltij(kllmn)
     end do !kllmn
     
     dterm%rd2c(iatom,ijlmn) = CMPLX(dijterm,zero)
   end do !ijlmn

   ABI_FREE(ff)
   ABI_FREE(nlff)

 end do ! end loop over iat

 dterm%has_rd2c = 2
 
end subroutine dterm_rd2c

!!****f* ABINIT/dterm_rd3a
!! NAME
!! dterm_rd3a
!!
!! FUNCTION
!! Compute onsite vxc
!! roadmap term 3a (no term 3b because usexcnhat 0 enforced)
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

subroutine dterm_rd3a(atindx,dterm,dtset,gntselect,gprimd,lmnmax,my_lmax,&
    & pawang,paw_ij,pawrad,pawsphden,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(paw_sph_den_type),intent(in) :: pawsphden(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,ilm,ilmn,imesh,ipt,itypat,jlm,jlmn
  integer :: klmn,kln,mesh_size,nkxc,nk3xc,usecore,usexcnhat,xc_option
  real(dp) :: dij,eexc,eexcdc,hyb_mixing,rr,xcint,xcintr
  logical :: non_magnetic_xc

  !arrays
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

 dterm%rd3a = czero

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

     dij = zero
     do ipt = 1, pawang%angl_size
       do imesh = 2, mesh_size
         rr = pawrad(itypat)%rad(imesh)
         ff(imesh) = vxc(imesh,ipt,1)*pawtab(itypat)%phiphj(imesh,kln) - &
           & tvxc(imesh,ipt,1)*pawtab(itypat)%tphitphj(imesh,kln)
       end do !imesh
       call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
       call simp_gen(xcint,ff,pawrad(itypat))

       dij = dij + &
         & four_pi*pawang%angwgth(ipt)*pawang%ylmr(ilm,ipt)*pawang%ylmr(jlm,ipt)*xcint
       
     end do ! ipt

     dterm%rd3a(iatom,klmn) = dij

   end do ! klmn
     
   ABI_FREE(vxc)
   ABI_FREE(tvxc)
   ABI_FREE(ff)
 
 end do !iat

 dterm%has_rd3a = 2

end subroutine dterm_rd3a

!!****f* ABINIT/dterm_rd1a
!! NAME
!! dterm_rd1a
!!
!! FUNCTION
!! Compute onsite terms due to A.p nuclear dipole moment
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

subroutine dterm_rd1a(atindx,dterm,dtset,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: cplex_dij,iat,iatom,ijlmn,itypat,ndij

  !arrays
  real(dp),allocatable :: dijnd(:,:)

!--------------------------------------------------------------------

 dterm%rd1a = czero

 cplex_dij = 2
 ndij = 1

 do iat=1,dtset%natom
   iatom=atindx(iat)
   itypat=dtset%typat(iat)

   ABI_MALLOC(dijnd,(cplex_dij*pawtab(itypat)%lmn2_size,ndij))

   call pawdijnd(dijnd,cplex_dij,ndij,dtset%nucdipmom(1:3,iat),pawrad(itypat),pawtab(itypat))

   do ijlmn=1,pawtab(itypat)%lmn2_size
     dterm%rd1a(iatom,ijlmn) = CMPLX(dijnd(2*ijlmn-1,1),dijnd(2*ijlmn,1))
   end do

   ABI_FREE(dijnd)

 end do !iat

 dterm%has_rd1a = 2

end subroutine dterm_rd1a


!!****f* ABINIT/dterm_dijhat
!! NAME
!! dterm_dijhat
!!
!! FUNCTION
!! Compute dijhat
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

subroutine dterm_dijhat(atindx,dterm,dtset,gprimd,nfftf,pawang,&
    & pawfgrtab,pawtab,ucvol,vtrial,vxc)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nfftf
  real(dp),intent(in) :: ucvol
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: cplex_dij,iat,iatom,ic,ilmn,ils,ilslm,isel,itypat
  integer :: jlmn,klm,klmn,lmax,lm_size,lmin
  real(dp) :: dijhat,vr

  !arrays
  real(dp),allocatable :: prod(:),v_dijhat(:,:)

!--------------------------------------------------------------------

 dterm%dijhat = czero

 cplex_dij = 2

 ABI_MALLOC(v_dijhat,(nfftf,dtset%nspden))
 ! as usexcnhat 0 is enforced, there is no vxc in interaction term in \hat{Dij}
 v_dijhat(1:nfftf,1:dtset%nspden) = vtrial(1:nfftf,1:dtset%nspden) - vxc(1:nfftf,1:dtset%nspden)

 do iat=1,dtset%natom
   iatom=atindx(iat)
   itypat=dtset%typat(iat)
   lm_size=pawtab(itypat)%lcut_size**2

   ABI_MALLOC(prod,(lm_size))

   !call pawdijhat(dijhat,cplex_dij,qphase,gprimd,iatom,dtset%natom,ndij,nfftf,nfftf,&
   !  & dtset%nspden,dtset%nsppol,pawang,pawfgrtab(iatom),pawtab(itypat),v_dijhat,qphon,ucvol,xred)

   ! compute Int[V(r).g_l(r).Y_lm(r)]
   prod=zero
   do ilslm=1,lm_size
     do ic=1,pawfgrtab(iatom)%nfgd
       vr=v_dijhat(pawfgrtab(iatom)%ifftsph(ic),1)
       prod(ilslm)=prod(ilslm)+vr*pawfgrtab(iatom)%gylm(ic,ilslm)
     end do
   end do
   !scaling factor (unit volume)
   prod=prod*ucvol/dble(nfftf)

   !compute Sum_(i,j)_LM { q_ij^L Int[V(r).g_l(r).Y_lm(r)] }
   do klmn=1,pawtab(itypat)%lmn2_size
     klm =pawtab(itypat)%indklmn(1,klmn)
     lmin=pawtab(itypat)%indklmn(3,klmn)
     lmax=pawtab(itypat)%indklmn(4,klmn)
     ilmn=pawtab(itypat)%indklmn(7,klmn)
     jlmn=pawtab(itypat)%indklmn(8,klmn)

     dijhat = zero
     do ils=lmin,lmax,2
       do ilslm=ils**2+1,(ils+1)**2
         isel=pawang%gntselect(ilslm,klm)
         if (isel>0) then 
           dijhat=dijhat+prod(ilslm)*pawtab(itypat)%qijl(ilslm,klmn)
         end if
       end do !ilslm
     end do ! ils
     dterm%dijhat(iatom,klmn) = CMPLX(dijhat,zero)
     
   end do ! klmn
 
   ABI_FREE(prod)

 end do !iat

 ABI_FREE(v_dijhat)

 dterm%has_dijhat = 2

end subroutine dterm_dijhat


!!****f* ABINIT/dterm_rd2a
!! NAME
!! dterm_rd2a
!!
!! FUNCTION
!! Compute onsite vhartree and r*vhartree
!! roadmap term 2a
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

subroutine dterm_rd2a(atindx,dterm,dtset,gntselect,lmnmax,my_lmax,&
    & pawrad,pawrhoij,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,igij,igkl,ijlm,ijlmax,ijlmin,ijlmn,ijln,itypat
  integer :: kllm,kllmax,kllmin,kllmn,klln,ll,llmm,mesh_size,meshsz
  real(dp) :: dijterm,eijkl,rgij,rgkl,vh1

  !arrays
  real(dp),allocatable :: ff(:),gg(:),hh(:)

!--------------------------------------------------------------------

 dterm%rd2a = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(gg,(mesh_size))
   ABI_MALLOC(hh,(mesh_size))
   
   meshsz=pawrad(itypat)%int_meshsz
   if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero

   do ijlmn = 1, pawtab(itypat)%lmn2_size
     ijlm   = pawtab(itypat)%indklmn(1,ijlmn)
     ijln   = pawtab(itypat)%indklmn(2,ijlmn)
     ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
     ijlmax = pawtab(itypat)%indklmn(4,ijlmn)

     dijterm = zero
     do kllmn = 1, pawtab(itypat)%lmn2_size
       kllm   = pawtab(itypat)%indklmn(1,kllmn)
       klln   = pawtab(itypat)%indklmn(2,kllmn)
       kllmin = pawtab(itypat)%indklmn(3,kllmn)
       kllmax = pawtab(itypat)%indklmn(4,kllmn)

       eijkl = zero
       do ll = max(ijlmin,kllmin),min(ijlmax,kllmax)
         do llmm = ll**2+1,(ll+1)**2
           igij = gntselect(llmm,ijlm)
           igkl = gntselect(llmm,kllm)
           if ( (igij .EQ. 0) .OR. (igkl .EQ. 0) ) cycle

           rgij=realgnt(igij)
           rgkl=realgnt(igkl)
  
           ff(1:meshsz)=pawtab(itypat)%phiphj  (1:meshsz,klln)
           call poisson(ff,ll,pawrad(itypat),gg)
           ff(1:meshsz)=pawtab(itypat)%tphitphj(1:meshsz,klln)
           call poisson(ff,ll,pawrad(itypat),hh)
           
           ff(1)=zero
           ff(2:meshsz)=(pawtab(itypat)%phiphj(2:meshsz,ijln)*gg(2:meshsz)&
             & -pawtab(itypat)%tphitphj(2:meshsz,ijln)*hh(2:meshsz))&
             & /pawrad(itypat)%rad(2:meshsz)
           call simp_gen(vh1,ff,pawrad(itypat))

           eijkl = eijkl + four_pi*vh1*rgij*rgkl

         end do !llmm
       end do !ll

       dijterm = dijterm + pawrhoij(iatom)%rhoij_(2*kllmn-1,1)*&
         & eijkl*pawtab(itypat)%dltij(kllmn)
     end do !kllmn

     dterm%rd2a(iatom,ijlmn) = CMPLX(dijterm,zero)

   end do !ijlmn

   ABI_FREE(ff)
   ABI_FREE(gg)
   ABI_FREE(hh)

 end do ! iat

 dterm%has_rd2a = 2
 
end subroutine dterm_rd2a

!!****f* ABINIT/dterm_rd2a1
!! NAME
!! dterm_rd2a1
!!
!! FUNCTION
!! Compute onsite vhartree due to first order density
!! roadmap term 2a
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

subroutine dterm_rd2a1(atindx,dterm,dtset,gntselect,gprimd,lmnmax,my_lmax,&
    & pawrad,pawrhoij1,pawtab,realgnt)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,my_lmax
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: gntselect((2*my_lmax-1)**2,my_lmax**2*(my_lmax**2+1)/2)
  real(dp),intent(in) :: gprimd(3,3),realgnt((2*my_lmax-1)**2*(my_lmax)**4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawrhoij_type),intent(in) :: pawrhoij1(dtset%natom,3)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,igij,igkl,ijlm,ijlmax,ijlmin,ijlmn,ijln,itypat
  integer :: klimn,kljmn,kllm,kllmax,kllmin,kllmn,klln,ll,llmm,mesh_size,meshsz
  real(dp) :: eijkl,rgij,rgkl,vh1
  complex(dpc) :: rhoij

  !arrays
  real(dp),allocatable :: ff(:),gg(:),hh(:)
  complex(dpc) :: dij_cart(3), dij_red(3)

!--------------------------------------------------------------------

 dterm%rd2a1 = czero

 do iat = 1, dtset%natom
   iatom = atindx(iat)
   itypat = dtset%typat(iat)
  
   mesh_size=pawtab(itypat)%mesh_size
   ABI_MALLOC(ff,(mesh_size))
   ABI_MALLOC(gg,(mesh_size))
   ABI_MALLOC(hh,(mesh_size))
   
   meshsz=pawrad(itypat)%int_meshsz
   if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero

   do ijlmn = 1, pawtab(itypat)%lmn2_size
     ijlm   = pawtab(itypat)%indklmn(1,ijlmn)
     ijln   = pawtab(itypat)%indklmn(2,ijlmn)
     ijlmin = pawtab(itypat)%indklmn(3,ijlmn)
     ijlmax = pawtab(itypat)%indklmn(4,ijlmn)

     dij_cart(:) = czero
     do kllmn = 1, pawtab(itypat)%lmn2_size
       kllm   = pawtab(itypat)%indklmn(1,kllmn)
       klln   = pawtab(itypat)%indklmn(2,kllmn)
       kllmin = pawtab(itypat)%indklmn(3,kllmn)
       kllmax = pawtab(itypat)%indklmn(4,kllmn)
       klimn = pawtab(itypat)%indklmn(7,kllmn)
       kljmn = pawtab(itypat)%indklmn(8,kllmn)

       eijkl = zero
       do ll = max(ijlmin,kllmin),min(ijlmax,kllmax)
         do llmm = ll**2+1,(ll+1)**2
           igij = gntselect(llmm,ijlm)
           igkl = gntselect(llmm,kllm)
           if ( (igij .EQ. 0) .OR. (igkl .EQ. 0) ) cycle

           rgij=realgnt(igij)
           rgkl=realgnt(igkl)
  
           ff(1:meshsz)=pawtab(itypat)%phiphj  (1:meshsz,klln)
           call poisson(ff,ll,pawrad(itypat),gg)
           ff(1:meshsz)=pawtab(itypat)%tphitphj(1:meshsz,klln)
           call poisson(ff,ll,pawrad(itypat),hh)
           
           ff(1)=zero
           ff(2:meshsz)=(pawtab(itypat)%phiphj(2:meshsz,ijln)*gg(2:meshsz)&
             & -pawtab(itypat)%tphitphj(2:meshsz,ijln)*hh(2:meshsz))&
             & /pawrad(itypat)%rad(2:meshsz)
           call simp_gen(vh1,ff,pawrad(itypat))

           eijkl = eijkl + four_pi*vh1*rgij*rgkl

         end do !llmm
       end do !ll

       do adir = 1, 3
         rhoij = CMPLX(pawrhoij1(iatom,adir)%rhoij_(2*kllmn-1,1),pawrhoij1(iatom,adir)%rhoij_(2*kllmn,1))
         if (klimn .EQ. kljmn) then
           dij_cart(adir) = dij_cart(adir) + rhoij*eijkl
         else
           dij_cart(adir) = dij_cart(adir) + two*real(rhoij)*eijkl
         end if
       end do
     end do !kllmn

     dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)
     dterm%rd2a1(iatom,ijlmn,1:3) = dij_red(1:3)

   end do !ijlmn

   ABI_FREE(ff)
   ABI_FREE(gg)
   ABI_FREE(hh)

 end do ! iat

 dterm%has_rd2a1 = 2
 
end subroutine dterm_rd2a1


!!****f* ABINIT/dterm_rd1
!! NAME
!! dterm_p2
!!
!! FUNCTION
!! Compute onsite p^2/2 
!! roadmap term 1
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

subroutine dterm_rd1(atindx,dterm,dtset,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iat,iatom,itypat,il,ilmn,ilm,iln,imesh,jl,jlm,jlmn,jln
  integer :: klm,klmn,msz,pmesh_size
  real(dp) :: angmom,avgkij,intff,rr

  !arrays
  real(dp),allocatable :: dtuj(:),d2tuj(:),duj(:),d2uj(:),ff(:),kij(:,:)
  real(dp),allocatable :: tui(:),tuj(:),ui(:),uj(:)

!--------------------------------------------------------------------

  dterm%rd1 = czero
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat = dtset%typat(iat)

    msz=pawrad(itypat)%mesh_size
    pmesh_size=pawtab(itypat)%partialwave_mesh_size
    ABI_MALLOC(ff,(msz))
    ABI_MALLOC(ui,(msz))
    ABI_MALLOC(uj,(msz))
    ABI_MALLOC(duj,(msz))
    ABI_MALLOC(d2uj,(msz))
    ABI_MALLOC(tui,(msz))
    ABI_MALLOC(tuj,(msz))
    ABI_MALLOC(dtuj,(msz))
    ABI_MALLOC(d2tuj,(msz))
    ABI_MALLOC(kij,(pawtab(itypat)%lmn_size,pawtab(itypat)%lmn_size))

    kij = zero

    do ilmn=1,pawtab(itypat)%lmn_size
      il =pawtab(itypat)%indlmn(1,ilmn)
      ilm=pawtab(itypat)%indlmn(4,ilmn)
      iln=pawtab(itypat)%indlmn(5,ilmn)
      ui(1:msz)=pawtab(itypat)%phi(1:msz,iln)
      tui(1:msz)=pawtab(itypat)%tphi(1:msz,iln)

      do jlmn=1,pawtab(itypat)%lmn_size
        jl =pawtab(itypat)%indlmn(1,jlmn)
        jlm=pawtab(itypat)%indlmn(4,jlmn)
        jln=pawtab(itypat)%indlmn(5,jlmn)
        klm=MATPACK(ilm,jlm)

        uj(1:msz) =pawtab(itypat)%phi(1:msz,jln)
        call nderiv_gen(duj,uj,pawrad(itypat),d2uj)

        tuj(1:msz)=pawtab(itypat)%tphi(1:msz,jln)
        call nderiv_gen(dtuj,tuj,pawrad(itypat),d2tuj)

        do imesh=2,msz
          rr = pawrad(itypat)%rad(imesh)
          angmom = jl*(jl+one)/(rr**2)
          ff(imesh) = ui(imesh)*(d2uj(imesh) - angmom*uj(imesh))
          ff(imesh)=ff(imesh) - &
            & tui(imesh)*(d2tuj(imesh) - angmom*tuj(imesh))
        end do

        call pawrad_deducer0(ff,msz,pawrad(itypat))
        call simp_gen(intff,ff,pawrad(itypat))
        
        ! ilm == jlm produces an rd1 term
        if (ilm .EQ. jlm) kij(ilmn,jlmn) = -half*intff

      end do !jlmn
    end do ! ilmn

    do klmn = 1, pawtab(itypat)%lmn2_size
      ilmn=pawtab(itypat)%indklmn(7,klmn)
      jlmn=pawtab(itypat)%indklmn(8,klmn)
      avgkij = half*(kij(ilmn,jlmn)+kij(jlmn,ilmn))
      dterm%rd1(iatom,klmn)=CMPLX(avgkij,zero)
    end do ! klmn
  
    ABI_FREE(ff)
    ABI_FREE(kij)
    ABI_FREE(ui)
    ABI_FREE(uj)
    ABI_FREE(duj)
    ABI_FREE(d2uj)
    ABI_FREE(tui)
    ABI_FREE(tuj)
    ABI_FREE(dtuj)
    ABI_FREE(d2tuj)

  end do !iat

  dterm%has_rd1 = 2

end subroutine dterm_rd1
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

  if(allocated(dterm%aij)) then
    ABI_FREE(dterm%aij)
  end if
  dterm%has_aij=0

  if(allocated(dterm%aij1)) then
    ABI_FREE(dterm%aij1)
  end if
  dterm%has_aij1=0

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if
  dterm%has_qij=0

  if(allocated(dterm%rd1)) then
    ABI_FREE(dterm%rd1)
  end if
  dterm%has_rd1=0
  
  if(allocated(dterm%rd1a)) then
    ABI_FREE(dterm%rd1a)
  end if
  dterm%has_rd1a=0
  
  if(allocated(dterm%rd2a)) then
    ABI_FREE(dterm%rd2a)
  end if
  dterm%has_rd2a=0
  
  if(allocated(dterm%rd2a1)) then
    ABI_FREE(dterm%rd2a1)
  end if
  dterm%has_rd2a1=0
  
  if(allocated(dterm%rd2b)) then
    ABI_FREE(dterm%rd2b)
  end if
  dterm%has_rd2b=0
  
  if(allocated(dterm%rd2c)) then
    ABI_FREE(dterm%rd2c)
  end if
  dterm%has_rd2c=0
  
  if(allocated(dterm%rd2d)) then
    ABI_FREE(dterm%rd2d)
  end if
  dterm%has_rd2d=0
  
  if(allocated(dterm%rd2e)) then
    ABI_FREE(dterm%rd2e)
  end if
  dterm%has_rd2e=0
  
  if(allocated(dterm%rd2f)) then
    ABI_FREE(dterm%rd2f)
  end if
  dterm%has_rd2f=0
  
  if(allocated(dterm%rd3a)) then
    ABI_FREE(dterm%rd3a)
  end if
  dterm%has_rd3a=0
  
  if(allocated(dterm%dijhat)) then
    ABI_FREE(dterm%dijhat)
  end if
  dterm%has_dijhat=0
  
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

subroutine dterm_alloc(dterm,lmnmax,lmn2max,natom)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: lmnmax,lmn2max,natom
  type(dterm_type),intent(inout) :: dterm

  !arrays

  !Local variables -------------------------
  !scalars
 
  !arrays
!--------------------------------------------------------------------

  dterm%lmnmax = lmnmax
  dterm%lmn2max = lmn2max
  dterm%natom = natom

  if(allocated(dterm%aij)) then
    ABI_FREE(dterm%aij)
  end if
  ABI_MALLOC(dterm%aij,(natom,lmn2max))
  dterm%has_aij=1

  if(allocated(dterm%aij1)) then
    ABI_FREE(dterm%aij1)
  end if
  ABI_MALLOC(dterm%aij1,(natom,lmn2max,3))
  dterm%has_aij1=1

  if(allocated(dterm%qij)) then
    ABI_FREE(dterm%qij)
  end if
  ABI_MALLOC(dterm%qij,(natom,lmn2max))
  dterm%has_qij=1

  if(allocated(dterm%rd1)) then
    ABI_FREE(dterm%rd1)
  end if
  ABI_MALLOC(dterm%rd1,(natom,lmn2max))
  dterm%has_rd1=1
  
  if(allocated(dterm%rd1a)) then
    ABI_FREE(dterm%rd1a)
  end if
  ABI_MALLOC(dterm%rd1a,(natom,lmn2max))
  dterm%has_rd1a=1
  
  if(allocated(dterm%rd2a)) then
    ABI_FREE(dterm%rd2a)
  end if
  ABI_MALLOC(dterm%rd2a,(natom,lmn2max))
  dterm%has_rd2a=1
  
  if(allocated(dterm%rd2a1)) then
    ABI_FREE(dterm%rd2a1)
  end if
  ABI_MALLOC(dterm%rd2a1,(natom,lmn2max,3))
  dterm%has_rd2a1=1
  
  if(allocated(dterm%rd2b)) then
    ABI_FREE(dterm%rd2b)
  end if
  ABI_MALLOC(dterm%rd2b,(natom,lmn2max))
  dterm%has_rd2b=1
  
  if(allocated(dterm%rd2c)) then
    ABI_FREE(dterm%rd2c)
  end if
  ABI_MALLOC(dterm%rd2c,(natom,lmn2max))
  dterm%has_rd2c=1
  
  if(allocated(dterm%rd2d)) then
    ABI_FREE(dterm%rd2d)
  end if
  ABI_MALLOC(dterm%rd2d,(natom,lmn2max))
  dterm%has_rd2d=1
  
  if(allocated(dterm%rd2e)) then
    ABI_FREE(dterm%rd2e)
  end if
  ABI_MALLOC(dterm%rd2e,(natom,lmn2max))
  dterm%has_rd2e=1
  
  if(allocated(dterm%rd2f)) then
    ABI_FREE(dterm%rd2f)
  end if
  ABI_MALLOC(dterm%rd2f,(natom,lmn2max))
  dterm%has_rd2f=1
  
  if(allocated(dterm%rd3a)) then
    ABI_FREE(dterm%rd3a)
  end if
  ABI_MALLOC(dterm%rd3a,(natom,lmn2max))
  dterm%has_rd3a=1
  
  if(allocated(dterm%dijhat)) then
    ABI_FREE(dterm%dijhat)
  end if
  ABI_MALLOC(dterm%dijhat,(natom,lmn2max))
  dterm%has_dijhat=1
 
  if(allocated(dterm%LR)) then
    ABI_FREE(dterm%LR)
  end if
  ABI_MALLOC(dterm%LR,(natom,lmn2max,3))
  dterm%has_LR=1
  
  if(allocated(dterm%BM)) then
    ABI_FREE(dterm%BM)
  end if
  ABI_MALLOC(dterm%BM,(natom,lmn2max,3))
  dterm%has_BM=1
  
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

!!****f* ABINIT/check_eig_k
!! NAME
!! check_eig_k
!!
!! FUNCTION
!! compute eig_k with paw_ij and dterm%aij, for comparison
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

subroutine check_eig_k(atindx,cg_k,cprj_k,dimlmn,dterm,dtset,eig_k,&
    & gs_hamk,ikpt,isppol,mcgk,mpi_enreg,my_nspinor,nband_k,npw_k,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer :: ikpt,isppol,mcgk,my_nspinor,nband_k,npw_k
  type(dataset_type),intent(in) :: dtset
  type(dterm_type),intent(in) :: dterm
  type(gs_hamiltonian_type),intent(inout) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
  real(dp),intent(in) :: cg_k(2,mcgk),eig_k(nband_k)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: cpopt,iat,iatom,itypat,klmn,ndat,nn,sij_opt,tim_getghc,type_calc
  real(dp) :: dotr_ghc,dotr_kloc,eig_aij,lambda
  complex(dpc) :: dij
  !arrays
  real(dp),allocatable :: cwave(:,:),ghc(:,:),gsc(:,:),gvnlxc(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)
!--------------------------------------------------------------------

 ABI_MALLOC(cwave,(2,npw_k))
 ABI_MALLOC(ghc,(2,npw_k))
 ABI_MALLOC(gsc,(2,npw_k))
 ABI_MALLOC(gvnlxc,(2,npw_k))
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,0,dimlmn)
 
 cpopt = 2 ! cprj available
 lambda = zero
 ndat = 1
 sij_opt = 0 ! don't need <g|S|v>
 tim_getghc = 0

  do nn = 1, nband_k
   
    cwave(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k) 
    call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
      & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
    
    type_calc = 0 ! full hamiltonian
    call getghc(cpopt,cwave,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,lambda,mpi_enreg,&
      & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)
    dotr_ghc = DOT_PRODUCT(cwave(1,1:npw_k),ghc(1,1:npw_k))+&
      &        DOT_PRODUCT(cwave(2,1:npw_k),ghc(2,1:npw_k))
    
    type_calc = 3 ! kinetic + local only
    call getghc(cpopt,cwave,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,lambda,mpi_enreg,&
      & ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)
    dotr_kloc = DOT_PRODUCT(cwave(1,1:npw_k),ghc(1,1:npw_k))+&
      &         DOT_PRODUCT(cwave(2,1:npw_k),ghc(2,1:npw_k))
    !eig_aij = dotr_kloc

    call tt_me(dterm%aij,atindx,cprj_k(:,nn),dtset,cprj_k(:,nn),&
           & dterm%lmn2max,pawtab,dij)
    eig_aij = dotr_kloc + REAL(dij)
    
    write(std_out,'(a,i4,es16.8)')'   JWZ debug nn eig_k : ',nn,eig_k(nn)
    write(std_out,'(a,i4,es16.8)')'     JWZ debug nn ghc : ',nn,dotr_ghc
    write(std_out,'(a,i4,es16.8)')'JWZ debug nn kloc+aij : ',nn,eig_aij
    write(std_out,'(2a)')'JWZ debug ',ch10

  end do ! nn

  ABI_FREE(cwave)
  ABI_FREE(ghc)
  ABI_FREE(gsc)
  ABI_FREE(gvnlxc)
  call pawcprj_free(cwaveprj)
  ABI_FREE(cwaveprj)

end subroutine check_eig_k
!!***

!!****f* ABINIT/sum_d
!! NAME
!! sum_d
!!
!! FUNCTION
!! sum various dterms into total aij and daij
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

subroutine sum_d(atindx,dterm,dtset,paw_ij,pawtab)

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
  integer :: iat,iatom,itypat,klmn
  !arrays
!--------------------------------------------------------------------

  dterm%aij = czero; dterm%aij1 = czero

  do iat=1,dtset%natom
    iatom=atindx(iat)
    itypat=dtset%typat(iat)
    do klmn=1,pawtab(itypat)%lmn2_size
      dterm%aij(iatom,klmn) = &
        & CMPLX(paw_ij(iatom)%dij(2*klmn-1,1),paw_ij(iatom)%dij(2*klmn,1))
    end do
  end do

end subroutine sum_d
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

subroutine make_d(atindx,atindx1,cprj,dimlmn,dterm,dtset,gprimd,mcprj,nfftf,&
    & mpi_enreg,occ,paw_an,pawang,pawfgrtab,paw_ij,pawrad,pawtab,psps,&
    & ucvol,vtrial,vxc,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcprj,nfftf
  real(dp),intent(in) :: ucvol
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
  integer,intent(in) :: dimlmn(dtset%natom)
  real(dp),intent(in) :: gprimd(3,3),xred(3,3)
  real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
  real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,iband,isel,itypat,ijlmn,kllmn
  integer :: my_lmax,ngnt,nzlmopt
  integer :: opt_compch,opt_dens,opt_l,opt_print
  real(dp) :: compch_sph,eijkl,hdij,my_eijkl,rrhokl
  complex(dpc) :: cdij,cpi,cpj,rhoij,rhokl
  type(paw_dmft_type) :: paw_dmft
 
  !arrays
  integer,allocatable :: gntselect(:,:)
  real(dp),allocatable :: realgnt(:)
  type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij1(:,:)
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

 ABI_MALLOC(pawrhoij1,(dtset%natom,3))
 do adir = 1, 3
   call pawrhoij_alloc(pawrhoij1(:,adir),dtset%pawcpxocc,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%typat,&
     & use_rhoijp=1,use_rhoij_=1,pawtab=pawtab)
 end do
 call make_pawrhoij1(atindx,cprj,dterm,dtset,mcprj,mpi_enreg,pawrhoij1,pawtab,ucvol)

 ! make Gaunt integrals
 !my_lmax = psps%mpsang + 1
 my_lmax = psps%mpsang + 2
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

 call dterm_qij(atindx,dterm,dtset,pawrad,pawtab)

 call dterm_LR(atindx,dterm,dtset,gprimd,pawrad,pawtab)
 
 call dterm_BM(atindx,dterm,dtset,gntselect,gprimd,my_lmax,pawrad,pawtab,realgnt)

 !!! term idp2 due to onsite p^2/2, corresponds to term 1 of Torrent PAW roadmap paper appendix E
 !!! Comp. Mat. Sci. 42, 337-351 (2008)
 !call dterm_rd1(atindx,dterm,dtset,pawrad,pawtab)

 !!! term idpa due to onsite A.p, not part of the roadmap paper but similar to term 1
 !call dterm_rd1a(atindx,dterm,dtset,pawrad,pawtab)
 !
 !! term idvha due to v_H[n1] and v_H[\tilde{n}1], corresponds to term 2a of roadmap paper
 !call dterm_rd2a(atindx,dterm,dtset,gntselect,psps%lmnmax,my_lmax,&
 !  & pawrad,pawrhoij,pawtab,realgnt)
 !!call dterm_rd2a1(atindx,dterm,dtset,gntselect,gprimd,psps%lmnmax,my_lmax,&
 !!  & pawrad,pawrhoij1,pawtab,realgnt)

 !! terms due to v_H[nZ], corresponds to term 2b of roadmap paper
 !call dterm_rd2b(atindx,dterm,dtset,pawrad,pawtab)
 !
 !! term due to v_H[nhat], corresponds to term 2c of roadmap paper
 !call dterm_rd2c(atindx,dterm,dtset,gntselect,psps%lmnmax,my_lmax,&
 !  & pawrad,pawrhoij,pawtab,realgnt)
 !
 !! term due to v_H[\tilde{n}^1]Q^{LM}_{ij}, corresponds to term 2d of roadmap paper
 !call dterm_rd2d(atindx,dterm,dtset,gntselect,psps%lmnmax,my_lmax,&
 !  & pawrad,pawrhoij,pawtab,realgnt)
 !
 !! terms due to v_H[\tilde{n}_Zc]Q^{LJM}_{ij}, corresponds to term 2e of roadmap paper
 !call dterm_rd2e(atindx,dterm,dtset,pawrad,pawtab)

 !! term due to v_H[\hat{n}]Q^{LM}_{ij}, corresponds to term 2f of roadmap paper
 !call dterm_rd2f(atindx,dterm,dtset,gntselect,psps%lmnmax,my_lmax,&
 !  & pawrad,pawrhoij,pawtab,realgnt)

 !! term dvxc due to v_xc[n1+nc] corresponds to term 3a in roadmap paper
 !call dterm_rd3a(atindx,dterm,dtset,gntselect,gprimd,psps%lmnmax,my_lmax,&
 !   & pawang,paw_ij,pawrad,pawsphden,pawtab,realgnt)

 !! term dijhat
 !call dterm_dijhat(atindx,dterm,dtset,gprimd,nfftf,pawang,pawfgrtab,pawtab,ucvol,&
 !  & vtrial,vxc)

 call sum_d(atindx,dterm,dtset,paw_ij,pawtab)

 ABI_FREE(realgnt)
 ABI_FREE(gntselect)
 
 do iat = 1, dtset%natom
   call paw_sph_den_free(pawsphden(iat))
 end do
 ABI_FREE(pawsphden)
 
 call pawrhoij_free(pawrhoij)
 ABI_FREE(pawrhoij)
 do adir = 1, 3
   call pawrhoij_free(pawrhoij1(:,adir))
 end do
 ABI_FREE(pawrhoij1)
 
end subroutine make_d
!!***

!!****f* ABINIT/make_pawrhoij1
!! NAME
!! make_pawrhoij1
!!
!! FUNCTION
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

subroutine make_pawrhoij1(atindx,cprj,dterm,dtset,mcprj,mpi_enreg,pawrhoij1,pawtab,ucvol)
  
  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcprj
  real(dp),intent(in) :: ucvol
  type(dterm_type),intent(inout) :: dterm
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrhoij_type),intent(inout) :: pawrhoij1(dtset%natom,3)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,gdir,iat,iatom,icprj,ikpt,itypat,ierr,ilmn,jlmn,klmn
  integer :: me,nband_k,nn,nproc,spaceComm
  real(dp) :: epsabg,trnrm
  complex(dpc) :: dcpi,dcpj,rhoij
 
  !arrays

!--------------------------------------------------------------------

  nband_k = dtset%mband
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt
 
  do iat = 1, dtset%natom
    iatom = atindx(iat)
    itypat = dtset%typat(iat)

    do klmn = 1, pawtab(itypat)%lmn2_size
      ilmn = pawtab(itypat)%indklmn(7,klmn)
      jlmn = pawtab(itypat)%indklmn(8,klmn)

      do adir = 1, 3

        rhoij = czero 
        do bdir = 1, 3
          do gdir = 1, 3
            epsabg = eijk(adir,bdir,gdir)
            if (ABS(epsabg) .LT. half) cycle

            icprj = 0
            do ikpt = 1, dtset%nkpt
              if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle
              
              trnrm = two*dtset%wtk(ikpt)

              do nn = 1, nband_k
                dcpi = CMPLX(cprj(iatom,icprj+nn)%dcp(1,bdir,ilmn),cprj(iatom,icprj+nn)%dcp(2,bdir,ilmn))
                dcpj = CMPLX(cprj(iatom,icprj+nn)%dcp(1,gdir,jlmn),cprj(iatom,icprj+nn)%dcp(2,gdir,jlmn))
                rhoij = rhoij + com*epsabg*trnrm*CONJG(dcpi)*dcpj
              end do ! nn
              
              icprj = icprj + nband_k
            end do ! ikpt

            if (nproc > 1) call xmpi_sum(rhoij,spaceComm,ierr)

          end do ! gdir
        end do ! bdir

        pawrhoij1(iatom,adir)%rhoij_(2*klmn-1,1) = REAL(rhoij)
        pawrhoij1(iatom,adir)%rhoij_(2*klmn,1) = AIMAG(rhoij)

        !write(std_out,'(a,2i4,2es16.8)')'JWZ debug adir klmn rhoij ',adir,klmn,real(rhoij),aimag(rhoij)

      end do ! adir

    end do !klmn
  end do ! iat
 
end subroutine make_pawrhoij1
!!***


end module m_orbmag
