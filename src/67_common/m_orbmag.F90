!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2021 ABINIT group (JWZ)
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
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,pawcprj_getdim, pawcprj_get
  use m_spacepar,         only : make_vectornd
  use m_time,             only : timab

  implicit none

  ! Bound methods:

  public :: orbmag_gipaw

  private :: orbmag_gipaw_pw_k
  private :: orbmag_gipaw_dpdk_k
  private :: orbmag_gipaw_onsite_l_k
  private :: orbmag_gipaw_onsite_bm_k
  private :: orbmag_gipaw_onsite_vh1_k
  private :: orbmag_gipaw_rhoij
  private :: orbmag_gipaw_lamb_core
  private :: orbmag_gipaw_output
  
CONTAINS  !========================================================================================
!!***

!!****f* ABINIT/orbmag_gipaw
!! NAME
!! orbmag_gipaw
!!
!! FUNCTION
!! This routine computes the orbital magnetization and Berry curvature based on input 
!! wavefunctions and DDK wavefuntions. It is based on using the modern theory of
!! orbital magnetism directly and the GIPAW expansion of the energy. 
!! It is assumed that only completely filled bands are present.
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

subroutine orbmag_gipaw(cg,cg1,cprj,dtset,gsqcut,kg,mcg,mcg1,mcprj,mpi_enreg,&
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
 integer :: adir,buff_size,dimffnl,exchn2n3d,iat,iatom,icg,icmplx,icprj,ider,idir,ierr
 integer :: ikg,ikg1,ikpt,ilm,indx,isppol,istwf_k,iterm,itypat
 integer :: me,mcgk,my_nspinor,nband_k,ncpgr,ndat,ngfft1,ngfft2,ngfft3,ngfft4
 integer :: ngfft5,ngfft6,nn,nkpg,npw_k,nproc,spaceComm,with_vectornd
 integer,parameter :: ibcpwb=1,ibcpws=2,iompwb=3,iompws=4,iomdp=5,iomlr=6,iombm=7,iomvh1=8,iomlmb=9
 integer,parameter :: nterms=9
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),oml(2,3),omlamb(2,3),omdp(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bcpw_k(:,:,:),bc_tensor(:,:,:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:)
 real(dp),allocatable :: Enk(:),ffnl_k(:,:,:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: omdp_k(:,:,:),ompw_k(:,:,:),ombm_k(:,:,:),oml_k(:,:,:),omvh1_k(:,:,:)
 real(dp),allocatable :: orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: rhoij(:,:,:,:),rhoij1(:,:,:,:,:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:)

 !----------------------------------------------

 write(std_out,'(a)')' JWZ debug: entered orbmag_gipaw '

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

 !------------------------------------------------------------
 ! make rhoij from input cprj
 ABI_MALLOC(rhoij,(2,psps%lmnmax,psps%lmnmax,dtset%natom))
 ABI_MALLOC(rhoij1,(2,psps%lmnmax,psps%lmnmax,dtset%natom,3))
 call orbmag_gipaw_rhoij(atindx,cprj,dtset,mcprj,mpi_enreg,pawtab,psps,rhoij,rhoij1,ucvol)

 ABI_MALLOC(kg_k,(3,dtset%mpw))
 ABI_MALLOC(kinpw,(dtset%mpw))

 ABI_MALLOC(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
 call getph(atindx,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

 icg = 0
 ikg = 0
 icprj = 0

 ABI_MALLOC(orbmag_terms,(2,nband_k,3,nterms))
 orbmag_terms = zero 

 if(dtset%orbmag .EQ. 4) then
   ABI_MALLOC(bc_tensor,(2,nband_k,dtset%nkpt,3)) 
   bc_tensor=zero
 end if
 
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

   ABI_MALLOC(Enk,(nband_k))
   ABI_MALLOC(ompw_k,(2,nband_k,6))
   ABI_MALLOC(bcpw_k,(2,nband_k,6))
   call orbmag_gipaw_pw_k(atindx,bcpw_k,cg_k,cg1_k,cprj_k,dimlmn,dtset,Enk,&
    & gs_hamk,ikpt,mcgk,mpi_enreg,nband_k,npw_k,ompw_k)

   orbmag_terms(1:2,1:nband_k,1:3,iompwb) = orbmag_terms(1:2,1:nband_k,1:3,iompwb) + &
     & trnrm*ompw_k(1:2,1:nband_k,1:3)
   orbmag_terms(1:2,1:nband_k,1:3,iompws) = orbmag_terms(1:2,1:nband_k,1:3,iompws) + &
     & trnrm*ompw_k(1:2,1:nband_k,4:6)
   orbmag_terms(1:2,1:nband_k,1:3,ibcpwb) = orbmag_terms(1:2,1:nband_k,1:3,ibcpwb) + &
     & trnrm*bcpw_k(1:2,1:nband_k,1:3)
   orbmag_terms(1:2,1:nband_k,1:3,ibcpws) = orbmag_terms(1:2,1:nband_k,1:3,ibcpws) + &
     & trnrm*bcpw_k(1:2,1:nband_k,4:6)

   if (dtset%orbmag .EQ. 4) then
     do nn = 1, nband_k
       do adir = 1, 3

         bc_tensor(1,nn,ikpt,adir) = bcpw_k(1,nn,adir)+bcpw_k(1,nn,adir+3)
         bc_tensor(2,nn,ikpt,adir) = bcpw_k(2,nn,adir)+bcpw_k(2,nn,adir+3)
         
       end do
       bc_tensor(1,nn,ikpt,1:3) = ucvol*MATMUL(gprimd,bc_tensor(1,nn,ikpt,1:3))
       bc_tensor(2,nn,ikpt,1:3) = ucvol*MATMUL(gprimd,bc_tensor(2,nn,ikpt,1:3))
     end do
   end if

   ABI_MALLOC(omdp_k,(2,nband_k,3))
   call orbmag_gipaw_dpdk_k(atindx,cprj_k,Enk,dtset%natom,nband_k,dtset%ntypat,&
     & omdp_k,paw_ij,pawtab,dtset%typat)
   orbmag_terms(1:2,1:nband_k,1:3,iomdp) = orbmag_terms(1:2,1:nband_k,1:3,iomdp) + &
     & trnrm*omdp_k(1:2,1:nband_k,1:3)

   ABI_MALLOC(oml_k,(2,nband_k,3))
   call orbmag_gipaw_onsite_l_k(atindx,cprj_k,dtset,nband_k,oml_k,pawrad,pawtab)
   orbmag_terms(1:2,1:nband_k,1:3,iomlr) = orbmag_terms(1:2,1:nband_k,1:3,iomlr) + &
     & trnrm*oml_k(1:2,1:nband_k,1:3)

   ABI_MALLOC(ombm_k,(2,nband_k,3))
   call orbmag_gipaw_onsite_bm_k(atindx,cprj_k,dtset,nband_k,ombm_k,pawang,pawrad,pawtab)
   orbmag_terms(1:2,1:nband_k,1:3,iombm) = orbmag_terms(1:2,1:nband_k,1:3,iombm) + &
     & trnrm*ombm_k(1:2,1:nband_k,1:3)

   ABI_MALLOC(omvh1_k,(2,nband_k,3))
   call orbmag_gipaw_onsite_vh1_k(atindx,cprj_k,dtset,nband_k,omvh1_k,pawtab,psps,rhoij1)
   orbmag_terms(1:2,1:nband_k,1:3,iomvh1) = orbmag_terms(1:2,1:nband_k,1:3,iomvh1) + &
     & trnrm*omvh1_k(1:2,1:nband_k,1:3)


   ABI_FREE(Enk)
   ABI_FREE(ompw_k)
   ABI_FREE(omdp_k)
   ABI_FREE(bcpw_k)
   ABI_FREE(oml_k)
   ABI_FREE(ombm_k)
   ABI_FREE(omvh1_k)

   icg = icg + npw_k*nband_k
   ikg = ikg + npw_k
   icprj = icprj + nband_k

   ABI_FREE(cg_k)
   ABI_FREE(cg1_k)
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

   if (dtset%orbmag .EQ. 4) then
     buff_size=size(bc_tensor)
     ABI_MALLOC(buffer1,(buff_size))
     ABI_MALLOC(buffer2,(buff_size))
     buffer1=zero;buffer2=zero
     buffer1(1:buff_size) = reshape(bc_tensor,(/2*nband_k*dtset%nkpt*3/))
     call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
     bc_tensor(1:2,1:nband_k,1:dtset%nkpt,1:3)=reshape(buffer2,(/2,nband_k,dtset%nkpt,3/))
     ABI_FREE(buffer1)
     ABI_FREE(buffer2)
   end if
 end if

 !! convert to cartesian frame from reduced triclinic
 ! iomlr, due to onsite L_R, is already Cartesian
 ! iombm, due to onsite A_0.A_N is already Cartesian
 ! iomlmb, due to core electrons, is already Cartesian (but not added yet)
 do iterm = 1, nterms
   if ((iterm .EQ. iomlr) .OR. (iterm .EQ. iombm) .OR. (iterm .EQ. iomlmb)) cycle
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
   if ((iterm.EQ.ibcpwb).OR.(iterm.EQ.ibcpws).OR.(iterm.EQ.iomlmb)) cycle
    orbmag_terms(1:2,1:nband_k,1:3,iterm) = ucvol*orbmag_terms(1:2,1:nband_k,1:3,iterm)
  end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace(1:2,1:3,1:nterms) = orbmag_trace(1:2,1:3,1:nterms) + orbmag_terms(1:2,nn,1:3,1:nterms)
 end do

 ! get the Lamb term
 call orbmag_gipaw_lamb_core(atindx,dtset,omlamb,pawtab)
 orbmag_trace(1:2,1:3,iomlmb) = omlamb(1:2,1:3)

 if (dtset%orbmag .EQ. 4) then
   call orbmag_gipaw_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace,bc_tensor)
 else
   call orbmag_gipaw_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)
 end if

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

 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)
 if (allocated(bc_tensor)) then
   ABI_FREE(bc_tensor) 
 end if

 ABI_FREE(rhoij)
 ABI_FREE(rhoij1)

 write(std_out,'(a)')' JWZ debug: exiting orbmag_gipaw '

end subroutine orbmag_gipaw
!!***

!!****f* ABINIT/orbmag_gipaw_pw_k
!! NAME
!! orbmag_gipaw_pw_k
!!
!! FUNCTION
!! Compute the planewave contribution to the orbital magnetization at one kpt
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

subroutine orbmag_gipaw_pw_k(atindx,bcpw_k,cg_k,cg1_k,cprj_k,dimlmn,dtset,Enk,&
    & gs_hamk,ikpt,mcgk,mpi_enreg,nband_k,npw_k,ompw_k)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: ikpt,mcgk,nband_k,npw_k
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type), intent(inout) :: mpi_enreg

 !arrays
 integer,intent(in) :: atindx(dtset%natom),dimlmn(dtset%natom)
 real(dp),intent(in) :: cg_k(2,mcgk),cg1_k(2,mcgk,3)
 real(dp),intent(out) :: bcpw_k(2,nband_k,6),Enk(nband_k),ompw_k(2,nband_k,6)
 type(pawcprj_type),intent(inout) :: cprj_k(dtset%natom,dtset%mband)

 !Local variables -------------------------
 !scalars
 integer :: adir,bbeg,bend,bdir,choice,cpopt,epsabg,gdir,icg,icprj,iband0,isppol,iscg,istwf_k
 integer :: my_ndat,my_nnlout,my_nspinor,my_timer,nbeg,ndat,nend,nn,nnp
 integer :: paw_opt,prtvol,scprod_io,signs,sij_opt,type_calc,useoverlap
 real(dp) :: c2,dotr,doti,lambda

 !arrays
 real(dp),allocatable :: bwavef(:,:),cwavef(:,:),dcg_k(:,:,:)
 real(dp),allocatable :: ghc(:,:),gsc(:,:),gvnlc(:,:),gwavef(:,:)
 real(dp),allocatable :: lambda_ndat(:),my_enlout(:),s1cg_k(:,:,:),scg_k(:,:),scprod(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

 !-----------------------------------------------------------------------

 ! various arguments that will be used below in various routines without change
 ! some of these could be generalized later
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 isppol = 1
 istwf_k = 1
 prtvol=0
 my_timer = 0
 my_ndat = 1
 iband0 = -1
 icg = 0
 icprj = 0
 iscg = 0
 useoverlap = 1
 scprod_io = 0

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in 
 ! the following terms, either <du|S|du> or <du|dS|u>
 c2=1.0d0/(two_pi*two_pi)

 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,cprj_k(1,1)%ncpgr,dimlmn)

 ! compute Enk = <u_nk|H|u_nk> 
 ! and s1cg_k = dS/dk|u_nk> 
 ! and dcg_k = -1/2 sum_n |u_nk><u_nk|dS/dk|u_mk>
 
 ABI_MALLOC(cwavef,(2,npw_k))
 ABI_MALLOC(ghc,(2,npw_k))
 ABI_MALLOC(gsc,(2,npw_k))
 ABI_MALLOC(gvnlc,(2,npw_k))
 ABI_MALLOC(scg_k,(2,mcgk))
 ABI_MALLOC(s1cg_k,(2,mcgk,3))
 ABI_MALLOC(dcg_k,(2,mcgk,3))
 cpopt=4;lambda=zero ! use cprj in memory, not using lambda
 sij_opt = 1 ! save gsc as well
 type_calc = 0 ! use all of Hamiltonian: kinetic, local, nonlocal

 choice = 5
 paw_opt = 3
 signs = 2
 my_nnlout = 1
 ABI_MALLOC(my_enlout,(my_nnlout))
 ABI_MALLOC(lambda_ndat,(my_ndat))
 lambda_ndat = zero
 dcg_k = zero
 
 do nn = 1, nband_k
   nbeg = (nn-1)*npw_k+1
   nend = nn*npw_k
   cwavef(1:2,1:npw_k) = cg_k(1:2,nbeg:nend)

   ! compute Enk and S|u_nk>
   call pawcprj_get(atindx,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
        &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
   call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,my_ndat,&
      &           prtvol,sij_opt,my_timer,type_calc)
   Enk(nn) = DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))
   scg_k(1:2,nbeg:nend) = gsc(1:2,1:npw_k)

   do adir = 1, 3
     ! compute dS/dk|u_nk>
     call nonlop(choice,cpopt,cwaveprj,my_enlout,gs_hamk,adir,lambda_ndat,mpi_enreg,my_ndat,my_nnlout,&
       & paw_opt,signs,gsc,my_timer,cwavef,gvnlc)

     s1cg_k(1:2,nbeg:nend,adir) = gsc(1:2,1:npw_k)

     ! now compute dcg_k = -1/2 \sum_np |u_np k><u_np k|dS/dk|u_nk> 
     ! see Audouze et al PRB 78, 035105 (2008) Eq. 40
     do nnp = 1, nband_k
       bbeg = (nnp-1)*npw_k+1
       bend = nnp*npw_k
       dotr = DOT_PRODUCT(cg_k(1,bbeg:bend),s1cg_k(1,nbeg:nend,adir))+&
         & DOT_PRODUCT(cg_k(2,bbeg:bend),s1cg_k(2,nbeg:nend,adir))
       doti = DOT_PRODUCT(cg_k(1,bbeg:bend),s1cg_k(2,nbeg:nend,adir))-&
         & DOT_PRODUCT(cg_k(2,bbeg:bend),s1cg_k(1,nbeg:nend,adir))

       dcg_k(1,nbeg:nend,adir) = dcg_k(1,nbeg:nend,adir) &
        & - half*dotr*cg_k(1,bbeg:bend) + half*doti*cg_k(2,bbeg:bend)
       dcg_k(2,nbeg:nend,adir) = dcg_k(2,nbeg:nend,adir) &
        & - half*dotr*cg_k(2,bbeg:bend) - half*doti*cg_k(1,bbeg:bend)

     end do ! end loop over nnp

   end do ! end loop over adir

 end do ! end loop over nn
 ABI_FREE(my_enlout)
 ABI_FREE(lambda_ndat)

 ABI_MALLOC(gwavef,(2,npw_k))
 ABI_MALLOC(bwavef,(2,npw_k))
 ABI_MALLOC(scprod,(2,nband_k))
 ompw_k = zero
 bcpw_k = zero
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

       ! gwavef and bwavef are the full first-order ddk wavefunctions, but
       ! we want the projection on the conduction states, Pc\psi(1). 
       ! Obtain this by subtracting dcg_k computed above, see
       ! Audouze et al PRB 78, 035105 (2008) Eq. 40
       gwavef(1:2,1:npw_k) = gwavef(1:2,1:npw_k) - dcg_k(1:2,nbeg:nend,gdir)
       bwavef(1:2,1:npw_k) = bwavef(1:2,1:npw_k) - dcg_k(1:2,nbeg:nend,bdir)

       ! compute H0(Pc gwavef), S0(Pc gwavef)
       cpopt = -1 ! need to compute but not save cprj
       call getghc(cpopt,gwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,my_ndat,&
         & prtvol,sij_opt,my_timer,type_calc)

       ! update gwavef as (H0 + En*S0)gwavef from ghc and gsc output of getghc
       gwavef(1:2,1:npw_k) = ghc(1:2,1:npw_k) + Enk(nn)*gsc(1:2,1:npw_k)

       ! store <du|P(H+ES)P|du> in terms 1:3
       ompw_k(1,nn,adir) = ompw_k(1,nn,adir) + half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),gwavef(2,:)) - DOT_PRODUCT(bwavef(2,:),gwavef(1,:)))
       ompw_k(2,nn,adir) = ompw_k(2,nn,adir) - half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),gwavef(1,:)) + DOT_PRODUCT(bwavef(2,:),gwavef(2,:)))
       ! store <du|P(E dS)|u> in terms 4:6
       ompw_k(1,nn,adir+3) = ompw_k(1,nn,adir+3) + half*epsabg*c2*Enk(nn)*&
         & (DOT_PRODUCT(bwavef(1,:),s1cg_k(2,nbeg:nend,gdir)) - DOT_PRODUCT(bwavef(2,:),s1cg_k(1,nbeg:nend,gdir)))
       ompw_k(2,nn,adir+3) = ompw_k(2,nn,adir+3) - half*epsabg*c2*Enk(nn)*&
         & (DOT_PRODUCT(bwavef(1,:),s1cg_k(1,nbeg:nend,gdir)) + DOT_PRODUCT(bwavef(2,:),s1cg_k(2,nbeg:nend,gdir)))

       ! store <du|S|du> in terms 1:3
       bcpw_k(1,nn,adir) = bcpw_k(1,nn,adir) + half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),gsc(2,:)) - DOT_PRODUCT(bwavef(2,:),gsc(1,:)))
       bcpw_k(2,nn,adir) = bcpw_k(2,nn,adir) - half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),gsc(1,:)) + DOT_PRODUCT(bwavef(2,:),gsc(2,:)))
       ! store <du|dS|u> in terms 4:6
       bcpw_k(1,nn,adir+3) = bcpw_k(1,nn,adir+3) + half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),s1cg_k(2,nbeg:nend,gdir)) - DOT_PRODUCT(bwavef(2,:),s1cg_k(1,nbeg:nend,gdir)))
       bcpw_k(2,nn,adir+3) = bcpw_k(2,nn,adir+3) - half*epsabg*c2*&
         & (DOT_PRODUCT(bwavef(1,:),s1cg_k(1,nbeg:nend,gdir)) + DOT_PRODUCT(bwavef(2,:),s1cg_k(2,nbeg:nend,gdir)))

     end do !end loop over epsabg
   end do ! end loop over adir
 end do ! end loop over nn

 ABI_FREE(gwavef)
 ABI_FREE(bwavef)
 ABI_FREE(scprod)

 ABI_FREE(cwavef)
 ABI_FREE(gvnlc)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(scg_k)
 ABI_FREE(s1cg_k)
 ABI_FREE(dcg_k)

 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

end subroutine orbmag_gipaw_pw_k
!!***

!!****f* ABINIT/orbmag_gipaw_dpdk_k
!! NAME
!! orbmag_gipaw_dpdk_k
!!
!! FUNCTION
!! Compute the orbital magnetization due to <u|dp/dk><dp/dk|u> at one k point
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

subroutine orbmag_gipaw_dpdk_k(atindx,cprj_k,Enk,natom,nband_k,&
    & ntypat,omdp_k,paw_ij,pawtab,typat)

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
             cterm = -c2*half*j_dpc*epsabg*(dij-Enk(nn)*pawtab(itypat)%sij(klmn))*conjg(udp_b)*udp_g
             omdp_k(1,nn,adir) = omdp_k(1,nn,adir) + REAL(cterm)
             omdp_k(2,nn,adir) = omdp_k(2,nn,adir) + AIMAG(cterm)

           end do ! end loop over jlmn
         end do ! end loop over ilmn
       end do ! end loop over atoms
     end do ! end loop over epsabg
   
   end do ! end loop over adir
 end do ! end loop over nn

end subroutine orbmag_gipaw_dpdk_k
!!***

!!****f* ABINIT/orbmag_gipaw_onsite_l_k
!! NAME
!! orbmag_gipaw_onsite_l_k
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

subroutine orbmag_gipaw_onsite_l_k(atindx,cprj_k,dtset,nband_k,omlk,pawrad,pawtab)

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
            if ( il /= jl ) cycle
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
 
end subroutine orbmag_gipaw_onsite_l_k
!!***

!!****f* ABINIT/orbmag_gipaw_onsite_bm_k
!! NAME
!! orbmag_gipaw_onsite_bm_k
!!
!! FUNCTION
!! Compute A_0.A_N onsite term for magnetic field + nuclear magnetic dipole moment
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

subroutine orbmag_gipaw_onsite_bm_k(atindx,cprj_k,dtset,nband_k,ombk,pawang,pawrad,pawtab)

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

end subroutine orbmag_gipaw_onsite_bm_k
!!***

!!****f* ABINIT/orbmag_gipaw_onsite_vh1_k
!! NAME
!! orbmag_gipaw_onsite_vh1_k
!!
!! FUNCTION
!! Compute first order Hartree contribution to orbital magnetization at given k point
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

subroutine orbmag_gipaw_onsite_vh1_k(atindx,cprj_k,dtset,nband_k,omvh1k,pawtab,psps,rhoij1)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dataset_type),intent(in) :: dtset
  type(pseudopotential_type), intent(inout) :: psps

  !arrays
  integer,intent(in) :: atindx(dtset%natom)
  real(dp),intent(in) :: rhoij1(2,psps%lmnmax,psps%lmnmax,dtset%natom,3)
  real(dp),intent(out) :: omvh1k(2,nband_k,3)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,iat,iatom,ijlmn,ilmn,itypat,jlmn,kllmn,klmn,llmn,nn
  complex(dpc) :: cpi,cpj,cterm,pij,pkl
  !arrays

!--------------------------------------------------------------------

  omvh1k = zero 
  do nn = 1, nband_k
    do adir = 1, 3
      do iat=1,dtset%natom
        iatom = atindx(iat)
        itypat = dtset%typat(iat)
        do ilmn=1,pawtab(itypat)%lmn_size
           do jlmn=1,pawtab(itypat)%lmn_size
              ijlmn=MATPACK(jlmn,ilmn)
              cpi = cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
              cpj = cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)     
              pij = conjg(cpi)*cpj
              do klmn = 1, pawtab(itypat)%lmn_size
                do llmn = 1, pawtab(itypat)%lmn_size
                  kllmn=MATPACK(klmn,llmn)
                  pkl = cmplx(rhoij1(1,klmn,llmn,iatom,adir),rhoij1(2,klmn,llmn,iatom,adir))
                  cterm = pij*pkl*pawtab(itypat)%eijkl(ijlmn,kllmn)
                  omvh1k(1,nn,adir) = omvh1k(1,nn,adir)-real(cterm)
                  omvh1k(2,nn,adir) = omvh1k(2,nn,adir)-aimag(cterm)
                end do ! end loop over llmn
              end do ! end loop over klmn
           end do ! end loop over jlmn
        end do ! end loop over ilmn
      end do ! end loop over atoms
    end do ! end loop over adir
  end do ! end loop over nn
 
end subroutine orbmag_gipaw_onsite_vh1_k
!!***

!****f* ABINIT/orbmag_gipaw_rhoj
!!! NAME
!! orbmag_gipaw_rhoij
!!
!! FUNCTION
!! make local rhoij and more importantly rhoij(1) due to B field
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

subroutine orbmag_gipaw_rhoij(atindx,cprj,dtset,mcprj,mpi_enreg,pawtab,psps,rhoij,rhoij1,ucvol)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcprj
 real(dp),intent(in) :: ucvol
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps

 !arrays
 integer,intent(in) :: atindx(dtset%natom)
 real(dp),intent(out) :: rhoij(2,psps%lmnmax,psps%lmnmax,dtset%natom)
 real(dp),intent(out) :: rhoij1(2,psps%lmnmax,psps%lmnmax,dtset%natom,3)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,bdir,buff_size,epsabg,gdir,iat,iatom,icprj,ierr,ikpt,ilmn,jlmn,itypat
 integer :: me,nband_k,nn,nproc,spaceComm
 real(dp) :: trnrm
 complex(dpc) :: c2,cpi,cpj,cterm

 !arrays
 real(dp),allocatable :: buffer1(:),buffer2(:)

 !----------------------------------------------

 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me = mpi_enreg%me_kpt
 nband_k = dtset%mband
 icprj = 0 
 
 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in rhoij1,
 ! one for each projector derivative)
 c2=1.0d0/(two_pi*two_pi)


 rhoij = zero
 rhoij1 = zero

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
       do ilmn=1,pawtab(itypat)%lmn_size
         do jlmn=1,pawtab(itypat)%lmn_size
           cpi = cmplx(cprj(iatom,icprj+nn)%cp(1,ilmn),cprj(iatom,icprj+nn)%cp(2,ilmn),KIND=dpc)
           cpj = cmplx(cprj(iatom,icprj+nn)%cp(1,jlmn),cprj(iatom,icprj+nn)%cp(2,jlmn),KIND=dpc)     
           cterm = trnrm*conjg(cpi)*cpj
           rhoij(1,ilmn,jlmn,iatom) = rhoij(1,ilmn,jlmn,iatom) + real(cterm)
           rhoij(2,ilmn,jlmn,iatom) = rhoij(2,ilmn,jlmn,iatom) + aimag(cterm)
           do adir = 1, 3
             do epsabg = 1, -1, -2
               if (epsabg .EQ. 1) then
                 bdir = modulo(adir,3)+1
                 gdir = modulo(adir+1,3)+1
               else
                 bdir = modulo(adir+1,3)+1
                 gdir = modulo(adir,3)+1
               end if
               cpi = cmplx(cprj(iatom,icprj+nn)%dcp(1,bdir,ilmn),cprj(iatom,icprj+nn)%dcp(2,bdir,ilmn),KIND=dpc)
               cpj = cmplx(cprj(iatom,icprj+nn)%dcp(1,gdir,jlmn),cprj(iatom,icprj+nn)%dcp(2,gdir,jlmn),KIND=dpc)     
               cterm = -trnrm*c2*half*j_dpc*epsabg*conjg(cpi)*cpj
               rhoij1(1,ilmn,jlmn,iatom,adir) = rhoij1(1,ilmn,jlmn,iatom,adir) + real(cterm)
               rhoij1(2,ilmn,jlmn,iatom,adir) = rhoij1(2,ilmn,jlmn,iatom,adir) + aimag(cterm)
             end do ! end loop over epsabg
           end do ! end loop over adir
         end do ! end loop over jlmn
       end do ! end loop over ilmn
     end do ! end loop over nn
   end do ! end loop over atoms
   
   icprj = icprj + nband_k

 end do ! end loop over k points

 if (nproc > 1) then
   buff_size=size(rhoij)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = reshape(rhoij,(/2*psps%lmnmax*psps%lmnmax*dtset%natom/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   rhoij(1:2,1:psps%lmnmax,1:psps%lmnmax,1:dtset%natom)=reshape(buffer2,(/2,psps%lmnmax,psps%lmnmax,dtset%natom/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)

   buff_size=size(rhoij1)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1=zero;buffer2=zero
   buffer1(1:buff_size) = reshape(rhoij1,(/2*psps%lmnmax*psps%lmnmax*dtset%natom*3/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   rhoij1(1:2,1:psps%lmnmax,1:psps%lmnmax,1:dtset%natom,1:3)=reshape(buffer2,(/2,psps%lmnmax,psps%lmnmax,dtset%natom,3/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

end subroutine orbmag_gipaw_rhoij

!!****f* ABINIT/orbmag_gipaw_lamb_core
!! NAME
!! orbmag_gipaw_onsite_lamb_core
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

subroutine orbmag_gipaw_lamb_core(atindx,dtset,omlamb,pawtab)

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
 
end subroutine orbmag_gipaw_lamb_core
!!***


!!****f* ABINIT/orbmag_gipaw_output
!! NAME
!! orbmag_gipaw_output
!!
!! FUNCTION
!! This routine outputs orbmag terms tailored for the gipaw ddk routine
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
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine orbmag_gipaw_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace&
    &,bc_tensor) ! bc_tensor is optional


 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nband_k,nterms
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: orbmag_terms(2,nband_k,3,nterms),orbmag_trace(2,3,nterms)
 real(dp),optional,intent(in) :: bc_tensor(2,nband_k,dtset%nkpt,3)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,ikpt,iterms
 integer,parameter :: ibcpwb=1,ibcpws=2,iompwb=3,iompws=4,iomdp=5,iomlr=6,iombm=7,iomvh1=8,iomlmb=9
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(2,nband_k,3),berry_total(2,3),orbmag_bb(2,nband_k,3),orbmag_total(2,3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 do iterms = iompwb, nterms
   orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
   do iband=1, nband_k
     if (iterms == iomlmb) cycle
     orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
   end do
 end do
 berry_bb=zero;berry_total=zero
 do iterms = ibcpwb,ibcpws
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
   write(message,'(a,3es16.8)') '  <du/dk|H+ES|du/dk> : ',(orbmag_trace(1,adir,iompwb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <du/dk|E dS/dk|u> : ',(orbmag_trace(1,adir,iompws),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') ' <u|dp/dk>D<dp/dk|u> : ',(orbmag_trace(1,adir,iomdp),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '   <u|p>1/2 L_R<p|u> : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     <u|p>A0.AN<p|u> : ',(orbmag_trace(1,adir,iombm),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '    <u|p>v_H(1)<p|u> : ',(orbmag_trace(1,adir,iomvh1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '           Lamb core : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Berry curvature, term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     <du/dk|S|du/dk> : ',(orbmag_trace(1,adir,ibcpwb),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '     <du/dk|dS/dk|u> : ',(orbmag_trace(1,adir,ibcpws),adir=1,3)
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
     write(message,'(a,3es16.8)') '      <du/dk|H+ES|du/dk> : ',(orbmag_terms(1,iband,adir,iompwb),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '       <du/dk|E dS/dk|u> : ',(orbmag_terms(1,iband,adir,iompws),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '     <u|dp/dk>D<dp/dk|u> : ',(orbmag_terms(1,iband,adir,iomdp),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '       <u|p>1/2 L_R<p|u> : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         <u|p>A0.AN<p|u> : ',(orbmag_terms(1,iband,adir,iombm),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '        <u|p>v_H(1)<p|u> : ',(orbmag_terms(1,iband,adir,iomvh1),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         Berry curvature : ',(berry_bb(1,iband,adir),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         <du/dk|S|du/dk> : ',(orbmag_terms(1,iband,adir,ibcpwb),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '         <du/dk|dS/dk|u> : ',(orbmag_terms(1,iband,adir,ibcpws),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 if(present(bc_tensor)) then
   write(std_out,'(a)') '# JWZ debug output full bc tensor '
   write(std_out,'(a)') '# iband kx ky kz Ox Oy Oz '
   do iband = 1, nband_k
     do ikpt = 1, dtset%nkpt
       write(std_out,'(i4,3f8.4,6es16.8)')iband,dtset%kptns(1,ikpt),dtset%kptns(2,ikpt),dtset%kptns(3,ikpt),&
         & bc_tensor(1,iband,ikpt,1),bc_tensor(2,iband,ikpt,1),&
         & bc_tensor(1,iband,ikpt,2),bc_tensor(2,iband,ikpt,2),&
         & bc_tensor(1,iband,ikpt,3),bc_tensor(2,iband,ikpt,3)
     end do
   end do
   write(std_out,'(a)') '# JWZ debug finished output full bc tensor '
 end if

end subroutine orbmag_gipaw_output
!!***

end module m_orbmag
