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
  use m_pawang,           only : pawang_type
 use m_pawfgr,            only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_pawrad,           only : nderiv_gen,pawrad_type,pawrad_deducer0,simp_gen
  use m_paw_sphharm,      only : setsym_ylm,slxyzs,realgaunt
  use m_pawtab,           only : pawtab_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get, pawcprj_put
  use m_spacepar,         only : make_vectornd
  use m_time,             only : timab

  implicit none

  ! Bound methods:

  public :: orbmag_tt

  private :: gs_eigenvalues
  private :: nonlocal_b_k
  private :: berry_curvature_k
  private :: onsite_LR_k
  private :: onsite_a0an_k
  private :: lamb_core
  private :: orbmag_tt_output
  
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
    & nfftf,ngfftf,npwarr,occ,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,rprimd,vtrial,&
    & xred,ylm,ylmgr)

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
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat
 integer :: me,mcgk,my_nspinor,nband_k,ncpgr,ndat,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 integer,parameter :: ibcpw=1,ibcq=2,ibcd=3,iombm=4,iomlmb=5,iomlr=6,iomnlb=7
 integer,parameter :: nterms=7
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),omlamb(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bcd_k(:,:,:),bcpw_k(:,:,:),bcq_k(:,:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: Enk(:),ffnl_k(:,:,:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: ombm_k(:,:,:),oml_k(:,:,:),omnlb_k(:,:,:)
 real(dp),allocatable :: orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
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

 ABI_MALLOC(kg_k,(3,dtset%mpw))
 ABI_MALLOC(kinpw,(dtset%mpw))

 ABI_MALLOC(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
 call getph(atindx,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

 icg = 0
 ikg = 0
 icprj = 0

 ABI_MALLOC(orbmag_terms,(2,nband_k,3,nterms))
 orbmag_terms = zero 

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
   choice = 1
   cpopt = 0
   idir = 0
   do adir = 1, 3
     do nn = 1, nband_k
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
   
   !--------------------------------------------------------------------------------
   ! Berry curvature
   !--------------------------------------------------------------------------------
   ABI_MALLOC(bcpw_k,(2,nband_k,3))
   ABI_MALLOC(bcq_k,(2,nband_k,3))
   ABI_MALLOC(bcd_k,(2,nband_k,3))
   call berry_curvature_k(atindx,bcd_k,bcpw_k,bcq_k,cg1_k,cprj_k,cprj1_k,dtset,gprimd,&
     & mcgk,nband_k,npw_k,pawtab,rprimd)

   orbmag_terms(:,:,:,ibcpw) = orbmag_terms(:,:,:,ibcpw) + trnrm*bcpw_k
   orbmag_terms(:,:,:,ibcq) = orbmag_terms(:,:,:,ibcq) + trnrm*bcq_k
   orbmag_terms(:,:,:,ibcd) = orbmag_terms(:,:,:,ibcd) + trnrm*bcd_k

   ABI_FREE(bcpw_k)
   ABI_FREE(bcq_k)
   ABI_FREE(bcd_k)

   !--------------------------------------------------------------------------------
   ! nonlocal contribution arising from GIPAW expansion
   !--------------------------------------------------------------------------------
   ABI_MALLOC(omnlb_k,(2,nband_k,3))
   call nonlocal_b_k(atindx,cprj_k,Enk,dtset%natom,nband_k,dtset%ntypat,omnlb_k,paw_ij,pawtab,dtset%typat)
   orbmag_terms(:,:,:,iomnlb) = orbmag_terms(:,:,:,iomnlb) + trnrm*omnlb_k
   ABI_FREE(omnlb_k)

   !--------------------------------------------------------------------------------
   ! onsite LR magnetic moment 
   !--------------------------------------------------------------------------------
   ABI_MALLOC(oml_k,(2,nband_k,3))
   call onsite_LR_k(atindx,cprj_k,dtset,nband_k,oml_k,pawrad,pawtab)
   orbmag_terms(:,:,:,iomlr) = orbmag_terms(:,:,:,iomlr) + trnrm*oml_k
   ABI_FREE(oml_k)

   !--------------------------------------------------------------------------------
   ! onsite A0.An magnetic moment 
   !--------------------------------------------------------------------------------
   ABI_MALLOC(ombm_k,(2,nband_k,3))
   call onsite_a0an_k(atindx,cprj_k,dtset,nband_k,ombm_k,pawang,pawrad,pawtab)
   orbmag_terms(:,:,:,iombm) = orbmag_terms(:,:,:,iombm) + trnrm*ombm_k
   ABI_FREE(ombm_k)

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
 ! iombm, due to onsite A_0.A_N is already Cartesian
 ! iomlr, due to onsite LR is already Cartesian
 ! iomlmb, due to core electrons, is already Cartesian (but not added yet)
 do iterm = 1, nterms
   if ((iterm .EQ. iombm).OR.(iterm .EQ. iomlmb).OR.(iterm .EQ. iomlr)) cycle
   do nn = 1, nband_k
     do icmplx = 1, 2
       orbmag_terms(icmplx,nn,1:3,iterm) = ucvol*MATMUL(gprimd,orbmag_terms(icmplx,nn,1:3,iterm))
     end do
   end do
 end do

 ! convert orbmag magnetization to orbital moment
 ! Berry curvature terms are ignored
 ! Lamb term ignored
 do iterm = 1, nterms
   if ((iterm.EQ.ibcpw).OR.(iterm.EQ.ibcq).OR.(iterm.EQ.ibcd).OR.(iterm.EQ.iomlmb)) cycle
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

 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

end subroutine orbmag_tt
!!***

!!****f* ABINIT/berry_curvature_k
!! NAME
!! berry_curvature_k
!!
!! FUNCTION
!! Compute the planewave contribution to the berry curvature at one kpt
!! Uses the DDK wavefunctions
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
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine berry_curvature_k(atindx,bcd_k,bcpw_k,bcq_k,cg1_k,cprj_k,cprj1_k,&
    & dtset,gprimd,mcgk,nband_k,npw_k,pawtab,rprimd)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcgk,nband_k,npw_k
 type(dataset_type),intent(in) :: dtset

 !arrays
 integer,intent(in) :: atindx(dtset%natom)
 real(dp),intent(in) :: cg1_k(2,mcgk,3),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: bcd_k(2,nband_k,3),bcpw_k(2,nband_k,3),bcq_k(2,nband_k,3)
 type(pawcprj_type),intent(inout) :: cprj_k(dtset%natom,dtset%mband)
 type(pawcprj_type),intent(inout) :: cprj1_k(dtset%natom,dtset%mband,3)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,epsabg,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nbeg,nend,nn
 real(dp) :: c2,cdij,dotr,doti,qij
 complex(dpc) :: cme,cpi,cpj,cpre,dup_b,dup_g,udp_b,udp_g

 !arrays
 integer,dimension(3) :: adir_to_sij = (/4,2,3/)
 real(dp),allocatable :: bwavef(:,:),gwavef(:,:)
 complex(dpc) :: dij_cart(3),dij_red(3)

 !-----------------------------------------------------------------------

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in 
 ! the following terms, either <du|S|du> or <du|dS|u>
 c2=1.0d0/(two_pi*two_pi)

 ! Berry curvature prefactor 
 cpre = j_dpc

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)

 ! the <du/dk| x |du/dk> term
 bcpw_k = zero
 ABI_MALLOC(gwavef,(2,npw_k))
 ABI_MALLOC(bwavef,(2,npw_k))
 do nn = 1, nband_k
   nbeg = (nn-1)*npw_k+1
   nend = nn*npw_k
   do adir = 1, 3
     do epsabg = 1, -1, -2
       if (epsabg .EQ. 1) then
         bdir = modulo(adir,3)+1
         gdir = modulo(adir+1,3)+1
       else
         bdir = modulo(adir+1,3)+1
         gdir = modulo(adir,3)+1
       end if

       gwavef(1:2,1:npw_k) = cg1_k(1:2,nbeg:nend,gdir)
       bwavef(1:2,1:npw_k) = cg1_k(1:2,nbeg:nend,bdir)

       ! gwavef and bwavef are the full first-order ddk wavefunctions
       dotr = DOT_PRODUCT(bwavef(1,:),gwavef(1,:)) + DOT_PRODUCT(bwavef(2,:),gwavef(2,:))
       doti = DOT_PRODUCT(bwavef(1,:),gwavef(2,:)) - DOT_PRODUCT(bwavef(2,:),gwavef(1,:))
       cme = cmplx(dotr,doti,KIND=dpc)

       bcpw_k(1,nn,adir) = bcpw_k(1,nn,adir) + real(cpre*epsabg*c2*cme)
       bcpw_k(2,nn,adir) = bcpw_k(2,nn,adir) + aimag(cpre*epsabg*c2*cme)

     end do
   end do
 end do
 ABI_FREE(gwavef)
 ABI_FREE(bwavef)

 ! the on-site Qij terms and (r-R)_\alpha terms
 bcd_k = zero; bcq_k = zero
 do nn = 1, nband_k
  do iat=1,dtset%natom
     iatom=atindx(iat)
     itypat=dtset%typat(iat)
     do ilmn=1,pawtab(itypat)%lmn_size
       do jlmn=1,pawtab(itypat)%lmn_size
         klmn=MATPACK(jlmn,ilmn)

         dij_cart(1) = -j_dpc*cdij*pawtab(itypat)%qijl(adir_to_sij(1),klmn)
         dij_cart(2) = -j_dpc*cdij*pawtab(itypat)%qijl(adir_to_sij(2),klmn)
         dij_cart(3) = -j_dpc*cdij*pawtab(itypat)%qijl(adir_to_sij(3),klmn)
         dij_red = MATMUL(TRANSPOSE(gprimd),dij_cart)

         do adir = 1, 3
           do epsabg = 1, -1, -2
             if (epsabg .EQ. 1) then
               bdir = modulo(adir,3)+1
               gdir = modulo(adir+1,3)+1
             else
               bdir = modulo(adir+1,3)+1
               gdir = modulo(adir,3)+1
             end if
 
             udp_b = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
             udp_g = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)     

             dup_b = cmplx(cprj1_k(iatom,nn,bdir)%cp(1,ilmn),cprj1_k(iatom,nn,bdir)%cp(2,ilmn),KIND=dpc)
             dup_g = cmplx(cprj1_k(iatom,nn,gdir)%cp(1,jlmn),cprj1_k(iatom,nn,gdir)%cp(2,jlmn),KIND=dpc)     

             qij = pawtab(itypat)%sij(klmn)
             
             cme = conjg(dup_b + udp_b)*qij*(dup_g + udp_g)
             bcq_k(1,nn,adir) = bcq_k(1,nn,adir) + real(cpre*epsabg*c2*cme)
             bcq_k(2,nn,adir) = bcq_k(2,nn,adir) + aimag(cpre*epsabg*c2*cme)
           
             cpi = cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
             cpj = cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)

             cme = conjg(dup_b+udp_b)*dij_red(gdir)*cpj + conjg(cpi*dij_red(bdir))*(dup_g+udp_g)

             bcd_k(1,nn,adir)=bcd_k(1,nn,adir) + real(cpre*epsabg*c2*cme)
             bcd_k(2,nn,adir)=bcd_k(2,nn,adir) + aimag(cpre*epsabg*c2*cme)

           end do !end loop over epsabg
         end do ! end loop over adir
       end do ! end loop over jlmn
     end do ! end loop over ilmn
   end do ! end loop over atoms
 end do ! end loop over nn
 
end subroutine berry_curvature_k
!!***

!!****f* ABINIT/onsite_a0an_k
!! NAME
!! onsite_aoan_k
!!
!! FUNCTION
!! Compute A_0.A_N onsite magnetic moment for magnetic field + nuclear magnetic dipole moment
!! for kpt
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

subroutine onsite_a0an_k(atindx,cprj_k,dtset,nband_k,ombk,pawang,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(pawang_type),intent(in) :: pawang
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: ombk(2,nband_k,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,gint,iat,iatom,il,im,ilmn,itypat
  integer :: jl,jm,jlmn,klmn,klm,kln,lpmp,mesh_size,nn
  real(dp) :: bm1,bm2,d00,d20,d22,dij,intg,scale_conversion
  complex(dpc) :: cpb,cpk,cterm

  !arrays
  real(dp),allocatable :: ff(:)

  ! ***********************************************************************

  ! this term can only be non-zero if some nucdipmom is nonzero
  scale_conversion = half*FineStructureConstant2
  d00 = sqrt(4.0*pi)/3.0
  dij = sqrt(4.0*pi/15.0)
  d20 = sqrt(16.0*pi/5.0)/6.0
  d22 = sqrt(16.0*pi/15.0)/2.0

  ombk = zero

  do nn = 1, nband_k
    do adir = 1, 3
      do iat=1,dtset%natom
        iatom=atindx(iat)
        itypat=dtset%typat(iat)
        mesh_size=pawtab(itypat)%mesh_size
        ABI_MALLOC(ff,(mesh_size))
        do jlmn=1,pawtab(itypat)%lmn_size
           jl=pawtab(itypat)%indlmn(1,jlmn)
           jm=pawtab(itypat)%indlmn(2,jlmn)
           do ilmn=1,pawtab(itypat)%lmn_size
              il=pawtab(itypat)%indlmn(1,ilmn)
              im=pawtab(itypat)%indlmn(2,ilmn)
              klmn=MATPACK(jlmn,ilmn)
              kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
              klm = pawtab(itypat)%indklmn(1,klmn) ! need this for bm2 gaunt integral selection
              ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r
              ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln) - &
                   &           pawtab(itypat)%tphitphj(2:mesh_size,kln)) / &
                   &           pawrad(itypat)%rad(2:mesh_size)
              call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
              call simp_gen(intg,ff,pawrad(itypat))
              ! term B.m r^2/r^3
              bm1=zero
              if ( (jl .EQ. il) .AND. (jm .EQ. im) .AND. (abs(dtset%nucdipmom(adir,iat)) .GT. tol8) ) then
                 bm1 = scale_conversion*dtset%nucdipmom(adir,iat)*intg
              end if
              bm2 = zero
              ! xx, yy, zz cases all have the same contribution from S00
              lpmp=1
              gint = pawang%gntselect(lpmp,klm)
              if (gint > 0) then
                bm2=bm2+scale_conversion*dtset%nucdipmom(adir,iat)*d00*pawang%realgnt(gint)*intg
              end if
              ! all other contributions involve Gaunt integrals of S_{2m}
              do lpmp = 5, 9
                ! if the basis is very small (only s waves) there can't be any Gaunt coupling to d functions
                ! and asking for it would crash the code
                if (size(pawang%gntselect(:,klm)) .LT. 5) exit
                gint = pawang%gntselect(lpmp,klm)
                if (gint > 0) then
                  select case (lpmp)
                  case (5) ! S_{2,-2} contributes to xy term
                    select case (adir)
                      case (1)
                        bm2=bm2+scale_conversion*dtset%nucdipmom(2,iat)*dij*pawang%realgnt(gint)*intg
                      case (2)
                        bm2=bm2+scale_conversion*dtset%nucdipmom(1,iat)*dij*pawang%realgnt(gint)*intg
                      end select
                    case (6) ! S_{2,-1} contributes to yz term
                      select case (adir)
                      case (2)
                        bm2=bm2+scale_conversion*dtset%nucdipmom(3,iat)*dij*pawang%realgnt(gint)*intg
                      case (3)
                        bm2=bm2+scale_conversion*dtset%nucdipmom(2,iat)*dij*pawang%realgnt(gint)*intg
                      end select
                    case (7) ! S_{2,0} contributes to xx, yy, and zz terms
                      select case (adir)
                        case (1)
                          bm2=bm2-scale_conversion*dtset%nucdipmom(1,iat)*d20*pawang%realgnt(gint)*intg
                        case (2)
                          bm2=bm2-scale_conversion*dtset%nucdipmom(2,iat)*d20*pawang%realgnt(gint)*intg
                        case (3)
                          bm2=bm2+scale_conversion*dtset%nucdipmom(3,iat)*2.0*d20*pawang%realgnt(gint)*intg
                      end select
                    case (8) ! S_{2,+1} contributes to xz term
                      select case (adir)
                        case (1)
                          bm2=bm2+scale_conversion*dtset%nucdipmom(3,iat)*dij*pawang%realgnt(gint)*intg
                        case (3)
                          bm2=bm2+scale_conversion*dtset%nucdipmom(1,iat)*dij*pawang%realgnt(gint)*intg
                      end select
                    case (9) ! S_{2,2} contributes to xx, yy terms
                       select case (adir)
                         case (1)
                           bm2=bm2+scale_conversion*dtset%nucdipmom(1,iat)*d22*pawang%realgnt(gint)*intg
                         case (2)
                           bm2=bm2-scale_conversion*dtset%nucdipmom(2,iat)*d22*pawang%realgnt(gint)*intg
                       end select
                    end select
                 end if ! end check on nonzero gaunt integral
              end do ! end loop over lp,mp
              cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
              cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)
              cterm = conjg(cpb)*(bm1-bm2)*cpk
              ombk(1,nn,adir)=ombk(1,nn,adir)-real(cterm)
              ombk(2,nn,adir)=ombk(2,nn,adir)-aimag(cterm)
           end do ! end loop over ilmn
        end do ! end loop over jlmn
        ABI_FREE(ff)
      end do ! end loop over atoms
    end do ! end loop over adir
  end do ! end loop over bands

end subroutine onsite_a0an_k
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

!!****f* ABINIT/onsite_LR_k
!! NAME
!! onsite_LR_k 
!!
!! FUNCTION
!! Compute 1/2 <L_R> onsite contribution to orbital magnetization at given k point
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

subroutine onsite_LR_k(atindx,cprj_k,dtset,nband_k,omlk,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(out) :: omlk(2,nband_k,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,ilmn,il,im,inn,itypat,jlmn,jl,jm,jn,klmn,kln,mesh_size,nn
  real(dp) :: intg
  complex(dpc) :: cpb,cpk,cterm,orbl_me

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

  omlk = zero 
  do nn = 1, nband_k
    do adir = 1, 3
      do iat=1,dtset%natom
        iatom = atindx(iat)
        itypat = dtset%typat(iat)
        mesh_size=pawtab(itypat)%mesh_size
        ABI_MALLOC(ff,(mesh_size))
        do jlmn=1,pawtab(itypat)%lmn_size
          jl=pawtab(itypat)%indlmn(1,jlmn)
          jm=pawtab(itypat)%indlmn(2,jlmn)
          jn=pawtab(itypat)%indlmn(3,jlmn)
          do ilmn=1,pawtab(itypat)%lmn_size
            il=pawtab(itypat)%indlmn(1,ilmn)
            if ( il /= jl ) cycle ! <l'm'|L|lm> = 0 if l' /= l
            if ( il == 0 ) cycle ! <00|L|00> = 0
            im=pawtab(itypat)%indlmn(2,ilmn)
            klmn=MATPACK(jlmn,ilmn)
            kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
            ! compute <L_dir>
            call slxyzs(il,im,adir,jl,jm,orbl_me)
            if(abs(orbl_me).GT.tol8)then
              ff(1:mesh_size) = pawtab(itypat)%phiphj(1:mesh_size,kln)-pawtab(itypat)%tphitphj(1:mesh_size,kln)
              call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
              call simp_gen(intg,ff,pawrad(itypat))
              cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
              cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)
              cterm = half*conjg(cpb)*orbl_me*intg*cpk
              omlk(1,nn,adir)=omlk(1,nn,adir) - real(cterm)
              omlk(2,nn,adir)=omlk(2,nn,adir) - aimag(cterm)
            end if ! end check that |L_dir| > 0, otherwise ignore term
          end do ! end loop over ilmn
        end do ! end loop over jlmn
        ABI_FREE(ff)
      end do ! end loop over atoms
    end do ! end loop over adir
  end do ! end loop over nn
 
end subroutine onsite_LR_k
!!***

!!****f* ABINIT/nonlocal_b_k
!! NAME
!! nonlocal_b_k
!!
!! FUNCTION
!! Compute the orbital magnetization due to <u|dp/dk><dp/dk|u> at one k point
!! arising from the GIPAW expansion, so arising from B dependence
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

subroutine nonlocal_b_k(atindx,cprj_k,Enk,natom,nband_k,ntypat,omdp_k,paw_ij,pawtab,typat)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: natom,nband_k,ntypat

 !arrays
 integer,intent(in) :: atindx(natom),typat(natom)
 real(dp),intent(in) :: Enk(nband_k)
 real(dp),intent(out) :: omdp_k(2,nband_k,3)
 type(pawcprj_type),intent(in) :: cprj_k(natom,nband_k)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,epsabg,gdir,iat,iatom,ilmn,itypat,jlmn,klmn,nn
 real(dp) :: c2,qij
 complex(dpc) :: cterm,dij,udp_b,udp_g
 logical :: cplex_dij

 !arrays

 !-----------------------------------------------------------------------

 cplex_dij = (paw_ij(1)%cplex_dij .EQ. 2)

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in the Chern number,
 ! one for each wavefunction derivative) 
 c2=1.0d0/(two_pi*two_pi)

 omdp_k = czero
 do nn = 1, nband_k
   do adir = 1, 3

     do epsabg = 1, -1, -2

       if (epsabg .EQ. 1) then
         bdir = modulo(adir,3)+1
         gdir = modulo(adir+1,3)+1
       else
         bdir = modulo(adir+1,3)+1
         gdir = modulo(adir,3)+1
       end if

       do iat=1,natom
         iatom=atindx(iat)
         itypat=typat(iat)
         ! paw_ij(iatom) seems be sorted in input atom order, not type order
         if (paw_ij(iat)%lmn2_size .NE. pawtab(itypat)%lmn2_size ) then
           ABI_BUG('lmn2_size mismatch in orbmag_gipaw_dpdk_k')
         end if
         do ilmn=1,pawtab(itypat)%lmn_size
           do jlmn=1,pawtab(itypat)%lmn_size
             klmn=MATPACK(jlmn,ilmn)

             udp_b = cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
             udp_g = cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)     

             if (cplex_dij) then
               dij = cmplx(paw_ij(iat)%dij(2*klmn-1,1),paw_ij(iat)%dij(2*klmn,1),KIND=dpc)
               if (jlmn .GT. ilmn) dij=conjg(dij)
             else
               dij = cmplx(paw_ij(iat)%dij(klmn,1),zero,KIND=dpc)
             end if

             ! the sign on cterm (+ or -) needs to be carefully checked
             cterm = c2*half*j_dpc*epsabg*(dij-Enk(nn)*pawtab(itypat)%sij(klmn))*conjg(udp_b)*udp_g
             omdp_k(1,nn,adir) = omdp_k(1,nn,adir) + REAL(cterm)
             omdp_k(2,nn,adir) = omdp_k(2,nn,adir) + AIMAG(cterm)

           end do ! end loop over jlmn
         end do ! end loop over ilmn
       end do ! end loop over atoms
     end do ! end loop over epsabg
   
   end do ! end loop over adir
 end do ! end loop over nn

end subroutine nonlocal_b_k
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
 integer,parameter :: ibcpw=1,ibcq=2,ibcd=3,iombm=4,iomlmb=5,iomlr=6,iomnlb=7
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
   !write(message,'(a,3es16.8)') '  <du/dk|H+ES|du/dk> : ',(orbmag_trace(1,adir,iompwb),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '   <du/dk|E dS/dk|u> : ',(orbmag_trace(1,adir,iompws),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   M_NL (from GIPAW) : ',(orbmag_trace(1,adir,iomnlb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <u|p>1/2 L_R<p|u> : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     <u|p>A0.AN<p|u> : ',(orbmag_trace(1,adir,iombm),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '    <u|p>v_H(1)<p|u> : ',(orbmag_trace(1,adir,iomvh1),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '           Lamb core : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    <du/dk|du/dk> : ',(orbmag_trace(1,adir,ibcpw),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  <du|p><dp|u>Qij : ',(orbmag_trace(1,adir,ibcq),adir=1,3)
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
     !write(message,'(a,3es16.8)') '      <du/dk|H+ES|du/dk> : ',(orbmag_terms(1,iband,adir,iompwb),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '       <du/dk|E dS/dk|u> : ',(orbmag_terms(1,iband,adir,iompws),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '       M_NL (from GIPAW) : ',(orbmag_terms(1,iband,adir,iomnlb),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '       <u|p>1/2 L_R<p|u> : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         <u|p>A0.AN<p|u> : ',(orbmag_terms(1,iband,adir,iombm),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '        <u|p>v_H(1)<p|u> : ',(orbmag_terms(1,iband,adir,iomvh1),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '    <du/dk|du/dk> : ',(orbmag_terms(1,iband,adir,ibcpw),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  <du|p><dp|u>Qij : ',(orbmag_terms(1,iband,adir,ibcq),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '     <du/dk|u>dij : ',(orbmag_terms(1,iband,adir,ibcd),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_tt_output
!!***

end module m_orbmag
