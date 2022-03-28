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

  private :: orbmag_tt_chern_k
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
 integer,parameter :: ibcpw=1,ibcq=2,ibcd=3
 integer,parameter :: nterms=3
 real(dp) :: arg,ecut_eff,trnrm,ucvol
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: atindx(:),atindx1(:),dimlmn(:),kg_k(:,:),nattyp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bcpw_k(:,:,:),bc_tensor(:,:,:,:),buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: Enk(:),ffnl_k(:,:,:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: orbmag_terms(:,:,:,:),orbmag_trace(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj1_k(:,:,:),cwaveprj(:,:)

 !----------------------------------------------

 write(std_out,'(a)')' JWZ debug: entered orbmag_tt '

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

   ABI_MALLOC(bcpw_k,(2,nband_k,9))
   call orbmag_tt_chern_k(atindx,bcpw_k,cg1_k,cprj_k,cprj1_k,dtset,gprimd,mcgk,nband_k,npw_k,pawtab,rprimd)

   orbmag_terms(1:2,1:nband_k,1:3,ibcpw) = orbmag_terms(1:2,1:nband_k,1:3,ibcpw) + &
     & trnrm*bcpw_k(1:2,1:nband_k,1:3)
   orbmag_terms(1:2,1:nband_k,1:3,ibcq) = orbmag_terms(1:2,1:nband_k,1:3,ibcq) + &
     & trnrm*bcpw_k(1:2,1:nband_k,4:6)
   orbmag_terms(1:2,1:nband_k,1:3,ibcd) = orbmag_terms(1:2,1:nband_k,1:3,ibcd) + &
     & trnrm*bcpw_k(1:2,1:nband_k,7:9)

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

   ABI_FREE(bcpw_k)

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

 ! convert orbmag magnetization to orbital moment
 ! Berry curvature terms are ignored
 ! Lamb term ignored
 !do iterm = 1, nterms
 !  if ((iterm.EQ.ibcpwb).OR.(iterm.EQ.ibcpws).OR.(iterm.EQ.iomlmb)) cycle
 !  orbmag_terms(1:2,1:nband_k,1:3,iterm) = ucvol*orbmag_terms(1:2,1:nband_k,1:3,iterm)
 !end do

 ! compute trace over filled states of each term
 ABI_MALLOC(orbmag_trace,(2,3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace(1:2,1:3,1:nterms) = orbmag_trace(1:2,1:3,1:nterms) + orbmag_terms(1:2,nn,1:3,1:nterms)
 end do

 ! get the Lamb term

 if (dtset%orbmag .EQ. 4) then
   call orbmag_tt_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace,bc_tensor)
 else
   call orbmag_tt_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace)
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
 do adir = 1, 3
   call pawcprj_free(cprj1_k(:,:,adir))
 end do
 ABI_FREE(cprj1_k)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)
 if (allocated(bc_tensor)) then
   ABI_FREE(bc_tensor) 
 end if

 write(std_out,'(a)')' JWZ debug: exiting orbmag_tt '

end subroutine orbmag_tt
!!***

!!****f* ABINIT/orbmag_tt_chern_k
!! NAME
!! orbmag_tt_chern_k
!!
!! FUNCTION
!! Compute the planewave contribution to the berry curvature at one kpt
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

subroutine orbmag_tt_chern_k(atindx,bcpw_k,cg1_k,cprj_k,cprj1_k,dtset,gprimd,mcgk,nband_k,npw_k,pawtab,rprimd)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcgk,nband_k,npw_k
 type(dataset_type),intent(in) :: dtset

 !arrays
 integer,intent(in) :: atindx(dtset%natom)
 real(dp),intent(in) :: cg1_k(2,mcgk,3),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: bcpw_k(2,nband_k,9)
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

 write(std_out,'(a)')' JWZ debug: entering orbmag_tt_chern_k '

 ! in abinit, exp(i k.r) is used not exp(i 2\pi k.r) so the following
 ! term arises to properly normalize the derivatives (there are two in 
 ! the following terms, either <du|S|du> or <du|dS|u>
 c2=1.0d0/(two_pi*two_pi)

 ! Chern prefactor = i/2\pi
 cpre = j_dpc/two_pi

 ! dij prefactor
 ! the pawtab%qijl moments do not have the sqrt(4\pi/3) factor we need here for
 ! normalization 
 cdij = sqrt(four_pi/three)

 bcpw_k = zero
 
 ! the <du/dk| x |du/dk> term
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
       ! store <du|du> in terms 1:3
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

 ! the on-site Qij terms
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
             bcpw_k(1,nn,adir+3) = bcpw_k(1,nn,adir+3) + real(cpre*epsabg*c2*cme)
             bcpw_k(2,nn,adir+3) = bcpw_k(2,nn,adir+3) + aimag(cpre*epsabg*c2*cme)
           
             cpi = cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
             cpj = cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)

             cme = conjg(dup_b+udp_b)*dij_red(gdir)*cpj + conjg(cpi*dij_red(bdir))*(dup_g+udp_g)

             bcpw_k(1,nn,adir+6)=bcpw_k(1,nn,adir+6) + real(cpre*epsabg*c2*cme)
             bcpw_k(2,nn,adir+6)=bcpw_k(2,nn,adir+6) + aimag(cpre*epsabg*c2*cme)

           end do !end loop over epsabg
         end do ! end loop over adir
       end do ! end loop over jlmn
     end do ! end loop over ilmn
   end do ! end loop over atoms
 end do ! end loop over nn
 
 write(std_out,'(a)')' JWZ debug: leaving orbmag_tt_chern_k '

end subroutine orbmag_tt_chern_k
!!***

!!****f* ABINIT/orbmag_tt_output
!! NAME
!! orbmag_tt_output
!!
!! FUNCTION
!! This routine outputs orbmag terms tailored for the tt ddk routine
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

subroutine orbmag_tt_output(dtset,nband_k,nterms,orbmag_terms,orbmag_trace&
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
 integer,parameter :: ibcpw=1,ibcq=2,ibcd=3
 character(len=500) :: message

 !arrays
 real(dp) :: berry_bb(2,nband_k,3),berry_total(2,3),orbmag_bb(2,nband_k,3),orbmag_total(2,3)

 ! ***********************************************************************

 orbmag_bb=zero;orbmag_total=zero
 !do iterms = iompwb, nterms
 !  orbmag_total(1:2,1:3)=orbmag_total(1:2,1:3) + orbmag_trace(1:2,1:3,iterms)
 !  do iband=1, nband_k
 !    if (iterms == iomlmb) cycle
 !    orbmag_bb(1:2,iband,1:3) = orbmag_bb(1:2,iband,1:3) + orbmag_terms(1:2,iband,1:3,iterms)
 !  end do
 !end do
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

 !write(message,'(a)')' Orbital magnetic moment, Cartesian directions : '
 !call wrtout(ab_out,message,'COLL')
 !write(message,'(3es16.8)') (orbmag_total(1,adir),adir=1,3)
 !call wrtout(ab_out,message,'COLL')
 !write(message,'(a)')ch10
 !call wrtout(ab_out,message,'COLL')
 write(message,'(a)')' Integral of Berry curvature, Cartesian directions : '
 call wrtout(ab_out,message,'COLL')
 write(message,'(3es16.8)') (berry_total(1,adir),adir=1,3)
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .GE. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   !write(message,'(a)')' Orbital magnetic moment, term-by-term breakdown : '
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '  <du/dk|H+ES|du/dk> : ',(orbmag_trace(1,adir,iompwb),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '   <du/dk|E dS/dk|u> : ',(orbmag_trace(1,adir,iompws),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') ' <u|dp/dk>D<dp/dk|u> : ',(orbmag_trace(1,adir,iomdp),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '   <u|p>1/2 L_R<p|u> : ',(orbmag_trace(1,adir,iomlr),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '     <u|p>A0.AN<p|u> : ',(orbmag_trace(1,adir,iombm),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '    <u|p>v_H(1)<p|u> : ',(orbmag_trace(1,adir,iomvh1),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a,3es16.8)') '           Lamb core : ',(orbmag_trace(1,adir,iomlmb),adir=1,3)
   !call wrtout(ab_out,message,'COLL')
   !write(message,'(a)')ch10
   !call wrtout(ab_out,message,'COLL')
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
     !write(message,'(a,3es16.8)') ' Orbital magnetic moment : ',(orbmag_bb(1,iband,adir),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '      <du/dk|H+ES|du/dk> : ',(orbmag_terms(1,iband,adir,iompwb),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '       <du/dk|E dS/dk|u> : ',(orbmag_terms(1,iband,adir,iompws),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '     <u|dp/dk>D<dp/dk|u> : ',(orbmag_terms(1,iband,adir,iomdp),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '       <u|p>1/2 L_R<p|u> : ',(orbmag_terms(1,iband,adir,iomlr),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '         <u|p>A0.AN<p|u> : ',(orbmag_terms(1,iband,adir,iombm),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a,3es16.8)') '        <u|p>v_H(1)<p|u> : ',(orbmag_terms(1,iband,adir,iomvh1),adir=1,3)
     !call wrtout(ab_out,message,'COLL')
     !write(message,'(a)')ch10
     !call wrtout(ab_out,message,'COLL')
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

end subroutine orbmag_tt_output
!!***

end module m_orbmag
