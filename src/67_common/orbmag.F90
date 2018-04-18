!{\src2tex{textfont=tt}}
!!****f* ABINIT/orbmag
!! NAME
!! orbmag
!!
!! FUNCTION
!! This routine computes the orbital magnetization based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtset <type(dataset_type)>=all input variables in this dataset
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! nfftf= - PAW only - number of FFT grid points for the "fine" grid (see NOTES at beginning of scfcv)
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawrad(ntypat*psps%usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!! usecprj=1 if cprj datastructure has been allocated
!! vhartr(nfftf)=Hartree potential
!! vpsp(nfftf)=array for holding local psp
!! vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!! xred(3,natom) = location of atoms in unit cell
!! ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!! ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006), and Gonze and Zwanziger, PRB 84
!! 064446 (2011). 
!! The derivative of the density operator is obtained from a discretized formula
!! $\partial_\beta \rho_k = \frac{1}{2\Delta}(\rho_{k+b} - \rho_{k-b})$ with
!! $\Delta = |b|$. When reduced to wavefunction overlaps the computation amounts to
!! multiple calls to smatrix.F90, exactly as in other Berry phase computations, with
!! the one additional complication of overlaps like $\langle u_{n,k+b}|u_{n',k+g}\rangle$.
!! At this stage mkpwind_k is invoked, which generalizes the code in initberry
!! and initorbmag necessary to index plane waves around different k points.
!! Direct questions and comments to J Zwanziger
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

subroutine orbmag(atindx1,cg,cprj,dtset,dtorbmag,kg,&
     &            mcg,mcprj,mpi_enreg,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,&
     &            pwind,pwind_alloc,rprimd,symrec,usecprj,vhartr,vpsp,vxc,xred,ylm,ylmgr)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_geometry,     only : metric
 use m_kg,               only : mkkin 

 use m_orbmag

 use m_hamiltonian,        only : init_hamiltonian,destroy_hamiltonian,&
&                                 load_spin_hamiltonian,load_k_hamiltonian,gs_hamiltonian_type

 use m_fftcore, only : kpgsph
 use m_paw_ij,           only : paw_ij_type
 use m_pawang,           only : pawang_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_free,&
      & pawcprj_get, pawcprj_getdim, pawcprj_set_zero, pawcprj_symkn

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'orbmag'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,nfftf,pwind_alloc,usecprj
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,mcg),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfftf),vpsp(nfftf),vxc(nfftf,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,bfor,bpw,bsigma,cpopt,ddkflag,dimffnl,epsabg,found,gdir,gfor,gpw,gsigma
 integer :: icg,icgb,icgg,icprj,icprjb,icprjg,ider,idir
 integer :: ikg,ikgb,ikgg,ikpt,ikptb,ikptg,ilm,ipw,isppol,istwf_k,itrs,job,jpw
 integer :: mcg1_k,my_nspinor,nband_k,ncpgr,ndat,nkpg,nn,n1,n2,n3
 integer :: ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,npw_k,npw_kb,npw_kg
 integer :: prtvol,shiftbd,sij_opt,tim_getghc,type_calc,type_calc_123
 real(dp) :: deltab,deltag,dotr,doti,htpisq,kenergy,lambda,ucvol
 complex(dpc) :: IIA,IIA1,IIA2,IIA3,IIIA,IIIA1,IIIA2,IIIA3,tprodIIA,tprodIIIA
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk,gs_hamk123
 !arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:),kg_kb(:,:),kg_kg(:,:),nattyp_dum(:),pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
 real(dp) :: dkb(3),dkg(3),dkbg(3),dtm_k(2),gmet(3,3),gprimd(3,3),kpg(3),kpoint(3),kpointb(3),orbmagvec(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bra(:,:),cg1_k(:,:),cgrvtrial(:,:),cwavef(:,:),ffnl(:,:,:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: hmat(:,:),hmat123(:,:,:,:,:),kinpw(:),kk_paw(:,:,:),kpg_k_dummy(:,:)
 real(dp),allocatable :: my_nucdipmom(:,:),ph3d(:,:,:),pwnsfac_k(:,:),smat_all(:,:,:,:,:),smat_inv(:,:,:)
 real(dp),allocatable :: smat_kk(:,:,:),vlocal(:,:,:,:),vtrial(:,:),ylm_k(:,:)
 complex(dpc),allocatable :: nucdipmom_k(:)
 logical,allocatable :: has_hmat(:),has_hmat123(:,:),has_smat(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:),cwaveprj(:,:)

 ! ***********************************************************************
 ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 
 nband_k = dtorbmag%mband_occ

 if (psps%usepaw == 1) then ! cprj allocation
    ncpgr = cprj(1,1)%ncpgr
    ABI_ALLOCATE(dimlmn,(dtset%natom))
    call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
    ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
    ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
    ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
    call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
    if (dtset%kptopt /= 3) then
       ABI_DATATYPE_ALLOCATE(cprj_ikn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
       ABI_DATATYPE_ALLOCATE(cprj_fkn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
       call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
       call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
    end if
    ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
    call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)
 else
   message = ' usepaw /= 1 but orbital magnetization calculation requires PAW '
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(pwind_kb,(dtset%mpw))
 ABI_ALLOCATE(pwind_kg,(dtset%mpw))
 ABI_ALLOCATE(pwind_bg,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

 mcg1_k = dtset%mpw*nband_k
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))
 ABI_ALLOCATE(sflag_k,(nband_k))
 ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

 ABI_ALLOCATE(kinpw,(dtset%mpw))

 ABI_ALLOCATE(my_nucdipmom,(3,dtset%natom))
 my_nucdipmom(:,:) = dtset%nucdipmom(:,:)

 ! input parameters for calls to smatrix.F90
 ddkflag = 1
 istwf_k = 1
 ! itrs = 0 means do not invoke time reversal symmetry in smatrix.F90
 itrs = 0
 job = 1
 shiftbd = 1

 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)

 ! input parameters for calls to getghc at ikpt
 cpopt = -1
 ndat = 1
 prtvol = 0
 sij_opt = 0
 tim_getghc = 0
 ! getghc: type_calc 0 means kinetic, local, nonlocal
 type_calc = 0
 lambda = zero

 htpisq = 0.5_dp*(two_pi)**2

 ABI_ALLOCATE(has_smat,(dtorbmag%fnkpt,dtorbmag%fnkpt))
 ABI_ALLOCATE(smat_all,(2,nband_k,nband_k,dtorbmag%fnkpt,dtorbmag%fnkpt))
 has_smat(:,:)=.FALSE.
 ABI_ALLOCATE(has_hmat,(dtorbmag%fnkpt))
 ABI_ALLOCATE(hmat,(nband_k,dtorbmag%fnkpt))
 has_hmat(:) = .FALSE.
 ABI_ALLOCATE(has_hmat123,(dtorbmag%fnkpt,dtorbmag%fnkpt))
 ABI_ALLOCATE(hmat123,(2,nband_k,nband_k,dtorbmag%fnkpt,dtorbmag%fnkpt))
 has_hmat123(:,:) = .FALSE.

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k, needed for computing E_nk
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=my_nucdipmom,&
      & paw_ij=paw_ij)

 !gs_hamk123 is used to compute <u_nk1|Hk2|u_mk3>. Eventually paw_ij will be the phase-twisted versions
 call init_hamiltonian(gs_hamk123,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=my_nucdipmom,&
      & paw_ij=paw_ij)

 !---------construct local potential------------------
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ! nspden=1 is essentially hard-coded in the following line
 vtrial(1:nfftf,1)=vhartr(1:nfftf)+vxc(1:nfftf,1)+vpsp(1:nfftf)
 
 ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
 call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)

 ABI_ALLOCATE(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
 call fftpac(isppol,mpi_enreg,dtset%nspden,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)

 ABI_DEALLOCATE(cgrvtrial)
 ABI_DEALLOCATE(vtrial)

 ! same vlocal in both Hamiltonians
 call load_spin_hamiltonian(gs_hamk,isppol,vlocal=vlocal,with_nonlocal=.true.)
 call load_spin_hamiltonian(gs_hamk123,isppol,vlocal=vlocal,with_nonlocal=.true.)

 !------- now local potential is attached to gs_hamk and gs_hamk123 -------------------------
 
 do adir = 1, 3

    IIA  = czero
    IIIA = czero

    do epsabg = 1, -1, -2

       if (epsabg .EQ. 1) then
          bdir = modulo(adir,3)+1
          gdir = modulo(adir+1,3)+1
       else
          bdir = modulo(adir+1,3)+1
          gdir = modulo(adir,3)+1
       end if

       ! loop over kpts, assuming for now kptopt 3, nsppol = 1, nspinor = 1
       ! and no parallelism, no symmorphic symmetry elements
       do ikpt = 1, dtorbmag%fnkpt

          icprj = dtorbmag%cprjindex(ikpt,isppol)
          
          npw_k = npwarr(ikpt)
          icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

          ikg = dtorbmag%fkgindex(ikpt)

          ! Set up remainder of normal Hamiltonian at k if necessary

          if (.NOT. has_hmat(ikpt) ) then
             kpoint(:)=dtorbmag%fkptns(:,ikpt)

             ABI_ALLOCATE(kg_k,(3,npw_k))
             ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
             kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
             if (psps%useylm==1) then
                do ilm=1,psps%mpsang*psps%mpsang
                   ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
                end do
             end if

             !      Compute (1/2) (2 Pi)**2 (k+G)**2:
!             ABI_ALLOCATE(kinpw,(npw_k))
             call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

             !  Compute (k+G) vectors (only if useylm=1)
             ! original code from vtorho.F90
             ! nkpg=3*optforces*dtset%nloalg(3)
             ! ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
             ! if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
             !    call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
             ! end if
             ! pretty sure do not need k+g vectors, use dummy kpg_k
             ! may eventually need them for nucdipmom_k hamiltonian if
             ! generalize to k,k'
             nkpg = 0
             ABI_ALLOCATE(kpg_k_dummy,(npw_k,nkpg))

             !      Compute nonlocal form factors ffnl at all (k+G):
             ider=0;idir=0;dimffnl=1
             ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
             call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                  &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k_dummy,kpoint,psps%lmnmax,&
                  &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                  &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                  &         psps%usepaw,psps%useylm,ylm_k,ylmgr)

             !     compute and load nuclear dipole Hamiltonian at current k point
             if(any(abs(gs_hamk%nucdipmom)>0.0)) then
                if(allocated(nucdipmom_k)) then
                   ABI_DEALLOCATE(nucdipmom_k)
                end if
                ABI_ALLOCATE(nucdipmom_k,(npw_k*(npw_k+1)/2))
                call mknucdipmom_k(gmet,kg_k,kpoint,dtset%natom,gs_hamk%nucdipmom,&
                     & nucdipmom_k,npw_k,rprimd,ucvol,xred)
                call load_k_hamiltonian(gs_hamk,nucdipmom_k=nucdipmom_k)
                ABI_DEALLOCATE(nucdipmom_k)
             end if
       

             !      Load k-dependent part in the Hamiltonian datastructure
             !       - Compute 3D phase factors
             !       - Prepare various tabs in case of band-FFT parallelism
             !       - Load k-dependent quantities in the Hamiltonian
             ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))

             call load_k_hamiltonian(gs_hamk,kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
                  &         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k_dummy,ffnl_k=ffnl,ph3d_k=ph3d,&
                  &         compute_ph3d=.TRUE.,compute_gbound=.TRUE.)

             ! apply gs_hamk to wavefunctions at k to compute E_nk eigenvalues
             ABI_ALLOCATE(cwavef,(2,npw_k))
             ABI_ALLOCATE(ghc,(2,npw_k))
             ABI_ALLOCATE(gsc,(2,npw_k))
             ABI_ALLOCATE(gvnlc,(2,npw_k))
             do nn = 1, nband_k
                cwavef(1,1:npw_k) = cg(1,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
                cwavef(2,1:npw_k) = cg(2,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
                call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
                     &                 prtvol,sij_opt,tim_getghc,type_calc)
                hmat(nn,ikpt)= DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) &
                     &       + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))
             end do
             has_hmat(ikpt) = .TRUE.

             ABI_DEALLOCATE(cwavef)
             ABI_DEALLOCATE(ghc)
             ABI_DEALLOCATE(gsc)
             ABI_DEALLOCATE(gvnlc)

             ABI_DEALLOCATE(kg_k)
             ABI_DEALLOCATE(ylm_k)
!             ABI_DEALLOCATE(kinpw)
             ABI_DEALLOCATE(kpg_k_dummy)
             ABI_DEALLOCATE(ffnl)
             if(any(abs(gs_hamk%nucdipmom)>0.0)) then
                if(allocated(nucdipmom_k)) then
                   ABI_DEALLOCATE(nucdipmom_k)
                end if
             end if
             ABI_DEALLOCATE(ph3d)

          end if ! end check on has_hmat
                    
          call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
               & dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

          do bfor = 1, 2
             if (bfor .EQ. 1) then
                bsigma = 1
             else
                bsigma = -1
             end if

             dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
             deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))

             ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
             icprjb = dtorbmag%cprjindex(ikptb,isppol)
             
             npw_kb = npwarr(ikptb)
             icgb = dtorbmag%cgindex(ikptb,dtset%nsppol)

             pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

             call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjb,&
                  & ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                  & my_nspinor,dtset%nsppol,0)

             if (.NOT. has_smat(ikpt,ikptb)) then
                
                call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                     & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                sflag_k=0
                call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
                     &  mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                     &  pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                do nn = 1, nband_k
                   do n1 = 1, nband_k
                      smat_all(1,nn,n1,ikpt,ikptb) =  smat_kk(1,nn,n1)
                      smat_all(2,nn,n1,ikpt,ikptb) =  smat_kk(2,nn,n1)
                      smat_all(1,n1,nn,ikptb,ikpt) =  smat_kk(1,nn,n1)
                      smat_all(2,n1,nn,ikptb,ikpt) = -smat_kk(2,nn,n1)
                   end do
                end do
                
                has_smat(ikpt,ikptb) = .TRUE.
                has_smat(ikptb,ikpt) = .TRUE.

             end if
             
             do gfor = 1, 2
                if (gfor .EQ. 1) then
                   gsigma = 1
                else
                   gsigma = -1
                end if

                dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
                deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))

                ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                icprjg = dtorbmag%cprjindex(ikptg,isppol)
             
                npw_kg = npwarr(ikptg)
                icgg = dtorbmag%cgindex(ikptg,dtset%nsppol)

                pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjg,&
                     & ikptg,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                     & my_nspinor,dtset%nsppol,0)

                if (.NOT. has_smat(ikpt,ikptg)) then
                
                   call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
                        &  mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                        &  pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   do nn = 1, nband_k
                      do n1 = 1, nband_k
                         smat_all(1,nn,n1,ikpt,ikptg) =  smat_kk(1,nn,n1)
                         smat_all(2,nn,n1,ikpt,ikptg) =  smat_kk(2,nn,n1)
                         smat_all(1,n1,nn,ikptg,ikpt) =  smat_kk(1,nn,n1)
                         smat_all(2,n1,nn,ikptg,ikpt) = -smat_kk(2,nn,n1)
                      end do
                   end do

                   has_smat(ikpt,ikptg) = .TRUE.
                   has_smat(ikptg,ikpt) = .TRUE.

                end if

                dkbg = dkg - dkb

                if (.NOT. has_smat(ikptb,ikptg)) then
                
                   call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        & dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,dtorbmag%indkk_f2ibz,ikptb,ikptg,&
                        & kg,dtorbmag%kgindex,mpi_enreg,npw_kb,pwind_bg,symrec)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
                        &  mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                        &  pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   do nn = 1, nband_k
                      do n1 = 1, nband_k
                         smat_all(1,nn,n1,ikptb,ikptg) =  smat_kk(1,nn,n1)
                         smat_all(2,nn,n1,ikptb,ikptg) =  smat_kk(2,nn,n1)
                         smat_all(1,n1,nn,ikptg,ikptb) =  smat_kk(1,nn,n1)
                         smat_all(2,n1,nn,ikptg,ikptb) = -smat_kk(2,nn,n1)
                      end do
                   end do

                   has_smat(ikptb,ikptg) = .TRUE.
                   has_smat(ikptg,ikptb) = .TRUE.

                end if

                if (.NOT. has_hmat123(ikptg,ikptb)) then

                   ! Compute (1/2) (2 Pi)**2 (k+G)**2
                   kpoint(:)=dtorbmag%fkptns(:,ikpt)
                   kpointb(:)=dtorbmag%fkptns(:,ikptb)
                   ikgb = dtorbmag%fkgindex(ikptb)
                   ABI_ALLOCATE(kg_kb,(3,npw_kb))
                   kg_kb(:,1:npw_kb)=kg(:,ikgb+1:ikgb+npw_kb)
                   ikgg = dtorbmag%fkgindex(ikptg)
                   ABI_ALLOCATE(kg_kg,(3,npw_kg))
                   kg_kg(:,1:npw_kg)=kg(:,ikgg+1:ikgg+npw_kg)
                   
                   ! ABI_ALLOCATE(kinpw,(npw_kb))
                   ! make kinetic energy at k + Gb 
                   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_kb,kinpw,kpoint,npw_kb,0,0)

                   call mkpwind_k(-dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,dtorbmag%indkk_f2ibz,ikptg,ikptb,&
                        & kg,dtorbmag%kgindex,mpi_enreg,npw_kg,pwind_bg,symrec)

                   !  Compute (k+G) vectors (only if useylm=1)
                   ! original code from vtorho.F90
                   ! nkpg=3*optforces*dtset%nloalg(3)
                   ! ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
                   ! if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
                   !    call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
                   ! end if
                   ! pretty sure do not need k+g vectors, use dummy kpg_k
                   ! may eventually need them for nucdipmom_k hamiltonian if
                   ! generalize to k,k'
                   nkpg = 0
                   ABI_ALLOCATE(kpg_k_dummy,(npw_kb,nkpg))

                   !     compute and load nuclear dipole Hamiltonian at current k point
                   if(any(abs(gs_hamk123%nucdipmom)>0.0)) then
                      if(allocated(nucdipmom_k)) then
                         ABI_DEALLOCATE(nucdipmom_k)
                      end if
                      ABI_ALLOCATE(nucdipmom_k,(npw_kb*(npw_kb+1)/2))
                      ! call mknucdipmom_k(gmet,kg_k,kpoint,dtset%natom,gs_hamk123%nucdipmom,&
                      !      & nucdipmom_k,npw_kb,rprimd,ucvol,xred)
                      nucdipmom_k = czero
                      call load_k_hamiltonian(gs_hamk123,nucdipmom_k=nucdipmom_k)
                      ABI_DEALLOCATE(nucdipmom_k)
                   end if

                   ABI_ALLOCATE(ylm_k,(npw_kb,psps%mpsang*psps%mpsang*psps%useylm))
                   if (psps%useylm==1) then
                      do ilm=1,psps%mpsang*psps%mpsang
                         ylm_k(1:npw_kb,ilm)=ylm(1+ikgb:npw_kb+ikgb,ilm)
                      end do
                   end if

                   ! !      Compute nonlocal form factors ffnl at all (kb+Gb):
                   ider=0;idir=0;dimffnl=1
                   ABI_ALLOCATE(ffnl,(npw_kb,dimffnl,psps%lmnmax,dtset%ntypat))
                   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                        &         gmet,gprimd,ider,idir,psps%indlmn,kg_kb,kpg_k_dummy,kpointb,psps%lmnmax,&
                        &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                        &         npw_kb,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                        &         psps%usepaw,psps%useylm,ylm_k,ylmgr)

                   !      Load k-dependent part in the Hamiltonian datastructure
                   !       - Compute 3D phase factors
                   !       - Prepare various tabs in case of band-FFT parallelism
                   !       - Load k-dependent quantities in the Hamiltonian
                   ABI_ALLOCATE(ph3d,(2,npw_kb,gs_hamk%matblk))

                   call load_k_hamiltonian(gs_hamk123,kpt_k=kpointb(:),istwf_k=istwf_k,npw_k=npw_kb,&
                        &         kinpw_k=kinpw,kg_k=kg_kb,kpg_k=kpg_k_dummy,ffnl_k=ffnl,ph3d_k=ph3d,&
                        &         compute_ph3d=.TRUE.,compute_gbound=.TRUE.)


                   ! apply gs_hamk123 to wavefunctions at kb 
                   ABI_ALLOCATE(cwavef,(2,npw_kb))
                   ABI_ALLOCATE(ghc,(2,npw_kb))
                   ABI_ALLOCATE(gsc,(2,npw_kb))
                   ABI_ALLOCATE(gvnlc,(2,npw_kb))

                   ABI_ALLOCATE(bra,(2,npw_kg))

                   ! getghc: type_calc_123 1 means local only, 3 means kinetic, local only
                   type_calc_123 = 0

                   do nn = 1, nband_k
                      cwavef(1,1:npw_kb) = cg(1,icgb+(nn-1)*npw_kb+1:icgb+nn*npw_kb)
                      cwavef(2,1:npw_kb) = cg(2,icgb+(nn-1)*npw_kb+1:icgb+nn*npw_kb)
                      ! apply vlocal
                      call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk123,gvnlc,lambda,mpi_enreg,ndat,&
                           &                 prtvol,sij_opt,tim_getghc,type_calc_123)
                      do n1 = 1, nband_k
                         bra(1,1:npw_kg) = cg(1,icgg+(n1-1)*npw_kg+1:icgg+n1*npw_kg)
                         bra(2,1:npw_kg) = cg(2,icgg+(n1-1)*npw_kg+1:icgg+n1*npw_kg)
                         dotr=zero;doti=zero
                         do ipw=1,npw_kg
                            jpw=pwind_bg(ipw)
                            if(jpw .GT. 0) then
                               dotr=dotr+bra(1,ipw)*ghc(1,jpw)+bra(2,ipw)*ghc(2,jpw)
                               doti=doti+bra(1,ipw)*ghc(2,jpw)-bra(2,ipw)*ghc(1,jpw)
                            end if
                         end do
                         hmat123(1,n1,nn,ikptg,ikptb) = dotr
                         hmat123(2,n1,nn,ikptg,ikptb) = doti
                         ! hmat123(1,nn,n1,ikptb,ikptg) = dotr
                         ! hmat123(2,nn,n1,ikptb,ikptg) = -doti
                      end do
                   end do
                   has_hmat123(ikptg,ikptb) = .TRUE.
                   ! has_hmat123(ikptb,ikptg) = .TRUE.

                   ABI_DEALLOCATE(cwavef)
                   ABI_DEALLOCATE(bra)
                   ABI_DEALLOCATE(ghc)
                   ABI_DEALLOCATE(gsc)
                   ABI_DEALLOCATE(gvnlc)

!                   ABI_DEALLOCATE(kinpw)
                   ABI_DEALLOCATE(kpg_k_dummy)
                   ABI_DEALLOCATE(ph3d)
                   ABI_DEALLOCATE(kg_kb)
                   ABI_DEALLOCATE(kg_kg)
                   ABI_DEALLOCATE(ylm_k)
                   ABI_DEALLOCATE(ffnl)
                   
                end if

                do nn = 1, nband_k
                   do n1 = 1, nband_k

                      IIA1  = cmplx(smat_all(1,nn,n1,ikpt,ikptg),smat_all(2,nn,n1,ikpt,ikptg))

                      IIIA1 = cmplx(smat_all(1,nn,n1,ikpt,ikptb),smat_all(2,nn,n1,ikpt,ikptb))

                      do n2 = 1, nband_k

                         IIA2 = cmplx(hmat123(1,n1,n2,ikptg,ikptb),hmat123(2,n1,n2,ikptg,ikptb))
                         IIA3 = cmplx(smat_all(1,n2,nn,ikptb,ikpt),smat_all(2,n2,nn,ikptb,ikpt))

                         IIIA2 = cmplx(smat_all(1,n1,n2,ikptb,ikptg),smat_all(2,n1,n2,ikptb,ikptg))
                         IIIA3 = conjg(cmplx(smat_all(1,nn,n2,ikpt,ikptg),smat_all(2,nn,n2,ikpt,ikptg)))

                         tprodIIA  = IIA1*IIA2*IIA3
                         IIA = IIA + epsabg*bsigma*gsigma*tprodIIA/(2.0*deltab*2.0*deltag)

                         tprodIIIA = hmat(nn,ikpt)*IIIA1*IIIA2*IIIA3
                         IIIA = IIIA - epsabg*bsigma*gsigma*tprodIIIA/(2.0*deltab*2.0*deltag)

                      end do ! end n2
                   end do ! end n1
                end do ! end nn
                
             end do ! end gfor

          end do ! end bfor

       end do ! end loop over fnkpt
    end do ! end loop over epsabg

    orbmagvec(1,adir) = real(IIA+IIIA)
    orbmagvec(2,adir) = aimag(IIA+IIIA)
    ! orbmagvec(1,adir) = real(IIIA)
    ! orbmagvec(2,adir) = aimag(IIIA)
 end do ! end loop over adir

 orbmagvec(1,1:3) = MATMUL(gprimd,orbmagvec(1,1:3))
 orbmagvec(2,1:3) = MATMUL(gprimd,orbmagvec(2,1:3))
 dtorbmag%orbmagvec(1,1:3) =  orbmagvec(2,1:3)/(two*dtorbmag%fnkpt)
 dtorbmag%orbmagvec(2,1:3) = -orbmagvec(1,1:3)/(two*dtorbmag%fnkpt)

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')
 
 write(message,'(a)')' Orbital magnetization '
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,a)')'----Orbital magnetization is a real vector, given along Cartesian directions----',ch10
 call wrtout(ab_out,message,'COLL')

 do adir = 1, 3
    write(message,'(a,i4,a,2es16.8)')' Orb Mag(',adir,') : real, imag ',&
         & dtorbmag%orbmagvec(1,adir),dtorbmag%orbmagvec(2,adir)
    call wrtout(ab_out,message,'COLL')
 end do

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')


 if (psps%usepaw == 1) then
    ABI_DEALLOCATE(dimlmn)
    call pawcprj_free(cprj_k)
    ABI_DATATYPE_DEALLOCATE(cprj_k)
    call pawcprj_free(cprj_kb)
    ABI_DATATYPE_DEALLOCATE(cprj_kb)
    call pawcprj_free(cprj_kg)
    ABI_DATATYPE_DEALLOCATE(cprj_kg)
    if (dtset%kptopt /= 3) then
       call pawcprj_free(cprj_ikn)
       call pawcprj_free(cprj_fkn)
       ABI_DATATYPE_DEALLOCATE(cprj_ikn)
       ABI_DATATYPE_DEALLOCATE(cprj_fkn)
    end if
    call pawcprj_free(cwaveprj)
    ABI_DATATYPE_DEALLOCATE(cwaveprj)
 end if

 ABI_DEALLOCATE(kk_paw)
 ABI_DEALLOCATE(cg1_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_inv)
 ABI_DEALLOCATE(smat_kk)
 ABI_DEALLOCATE(pwind_kb)
 ABI_DEALLOCATE(pwind_kg)
 ABI_DEALLOCATE(pwind_bg)
 ABI_DEALLOCATE(pwnsfac_k)

 ABI_DEALLOCATE(kinpw)
 
 ABI_DEALLOCATE(has_smat)
 ABI_DEALLOCATE(smat_all)
 ABI_DEALLOCATE(has_hmat)
 ABI_DEALLOCATE(hmat)
 ABI_DEALLOCATE(has_hmat123)
 ABI_DEALLOCATE(hmat123)

 ABI_DEALLOCATE(my_nucdipmom)

 ABI_DEALLOCATE(vlocal)
 call destroy_hamiltonian(gs_hamk)
 call destroy_hamiltonian(gs_hamk123)

end subroutine orbmag
!!***
