!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigmaph
!! NAME
!!  m_sigmaph
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sigmaph

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ebands
 use m_fstab
 use m_wfk
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth
 use m_io_tools,       only : iomode_from_fname
 use m_special_funcs,  only : dirac_delta
 use m_fftcore,        only : ngfft_seq
 use m_cgtools,        only : dotprod_g !set_istwfk
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_fftcore,        only : get_kg
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
! use m_paw_an,		     only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
! use m_paw_ij,		     only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
! use m_pawfgrtab,	   only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
! use m_pawrhoij,	     only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
! use m_pawdij,		     only : pawdij, symdij

 implicit none

 private
!!***

 public :: sigmaph_driver
!!***

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/sigmaph_driver
!! NAME
!!  sigmaph_driver
!!
!! FUNCTION
!!  Compute phonon-contribution to the electron self-energy.
!!
!! INPUTS
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigmaph_driver(wfk0_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
                          pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigmaph_driver'
 use interfaces_14_hidewrite
 use interfaces_56_recipspace
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),target,intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: dummy_npw=1,nsig=1,tim_getgh1c=1,berryopt0=0,timrev1=1,brav1=1
 integer,parameter :: useylmgr=0,useylmgr1=0,master=0,option1=1,ndat1=1,copt0=0
 integer :: my_rank,nproc,iomode,mband,my_minb,my_maxb,nsppol,nkpt,idir,ipert,iq_ibz
 integer :: cplex,db_iqpt,natom,natom3,ipc,ipc1,ipc2,nspinor,onpw
 integer :: ib1,ib2,band,num_smallw
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq,nqpt_max !timerev_q,
 integer :: spin,istwf_k,istwf_kq,istwf_kirr,npw_k,npw_kq,npw_kirr
 integer :: ii,ipw,mpw,my_mpw,ierr,cnt,ncid
 integer :: n1,n2,n3,n4,n5,n6,nspden,do_ftv1q
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,nqibz_k,nqbz
 integer :: sph_nb,sph_bstart,sph_nshiftq,ntemp
 real(dp) :: cpu,wall,gflops
 real(dp) :: ecut,eshift,eig0nk,dotr,doti,dksqmax,tstep
 logical,parameter :: have_ktimerev=.True.
 logical :: isirr_k,isirr_kq,gen_eigenpb
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 !type(phgamma_t) :: gams
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0bz_kq(3),g0_kq(3),symq(4,2,cryst%nsym),dummy_gvec(3,dummy_npw)
 integer :: qptrlatt(3,3)
 integer :: work_ngfft(18),gmax(3),my_gmax(3) !g0ibz_kq(3),
 integer :: sph_ngqpt(3)
 integer :: indkk_k(1,6),indkk_kq(1,6)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom),lf(2),rg(2),res(2)
 real(dp) :: temp_range(2)
 real(dp) :: eta,f_nk,f_mkq !,omega,wtk,gkk2,term1,term2
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom),displ_red(2,3,cryst%natom,3*cryst%natom)
 !real(dp) :: temp_tgam(2,3*cryst%natom,3*cryst%natom)
 real(dp),allocatable :: sph_shiftq(:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),tgam(:,:,:,:),gvals_qibz(:,:,:,:,:,:)
 real(dp),allocatable :: gkk_atm(:,:,:),gkk_mu(:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:),kets_k(:,:,:),h1kets_kq(:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:),work(:,:,:,:)
 real(dp),allocatable ::  gs1c(:,:) !,gvnl_direc(:,:),pcon(:),sconjgr(:,:)
 !real(dp),allocatable :: eloc0_k(:),enl0_k(:),enl1_k(:),vlocal_tmp(:,:,:),vlocal1_tmp(:,:,:), rho1wfg(:,:),rho1wfr(:,:)
 !real(dp),allocatable :: wt_kq(:,:)
 real(dp),allocatable :: qibz_k(:,:),qbz(:,:),tmesh(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)
 !real(dp),allocatable :: cwave0(:,:),gvnl1(:,:)

!************************************************************************

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkpt = ebands%nkpt; mband=ebands%mband

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 ! Build (linear) mesh of temperatures.
 ntemp = 10
 ABI_CHECK(ntemp > 1, "ntemp cannot be 1")
 temp_range = [0.6_dp, 1.2_dp]
 tstep = (temp_range(2) - temp_range(1)) / (ntemp - 1)
 ABI_MALLOC(tmesh, (ntemp))
 tmesh = arth(temp_range(1), tstep, ntemp)

 ! Broadening parameter
 if (dtset%elph2_imagden .gt. tol12) then
   eta = dtset%elph2_imagden
 else
   eta = 0.0001_dp
 end if

 ! Define q-mesh for integration.
 ! Note: The paramenters used to compute the self-energy start with sph
 ! sph stands for Sigma_phonon.
 sph_ngqpt = ifc%ngqpt; if (all(dtset%eph_ngqpt_fine /= 0)) sph_ngqpt = dtset%eph_ngqpt_fine
 sph_nshiftq = 1
 ABI_MALLOC(sph_shiftq, (3, sph_nshiftq))
 sph_shiftq = 0

 qptrlatt = 0
 qptrlatt(1,1) = sph_ngqpt(1)
 qptrlatt(2,2) = sph_ngqpt(2)
 qptrlatt(3,3) = sph_ngqpt(3)
 nqpt_max = (product(sph_ngqpt) * sph_nshiftq)

 ABI_ALLOCATE(qbz, (3,nqpt_max))
 call smpbz(brav1,std_out,qptrlatt,nqpt_max,nqbz,sph_nshiftq,option1,sph_shiftq,qbz)

 ! For the time being we do not use symmetries --> nqibz_k == nqbz
 nqibz_k = nqbz
 ABI_ALLOCATE(qibz_k, (3, nqibz_k))
 qibz_k = qbz(:,1:nqbz)
 ABI_FREE(qbz)
 ABI_FREE(sph_shiftq)
 call wrtout(std_out, sjoin("Will compute", itoa(nqibz_k), "q-points in the IBZ(k)"))

 ! Open the DVDB file
 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)

 ! Activate Fourier interpolation if irred q-points are not in the DVDB file.
 do_ftv1q = 0
 do iq_ibz=1,nqibz_k
   qpt = qibz_k(:,iq_ibz)
   if (dvdb_findq(dvdb, qpt) == -1) do_ftv1q = do_ftv1q + 1
 end do
 if (do_ftv1q /= 0) then
   write(msg, "(2(a,i0),a)")"Will use Fourier interpolation of DFPT potentials [",do_ftv1q,"/",nqibz_k,"]"
   call wrtout(std_out, msg)
   call wrtout(std_out, sjoin("From ngqpt", ltoa(ifc%ngqpt), "to", ltoa(sph_ngqpt)))
   call dvdb_ftinterp_setup(dvdb,ifc%ngqpt,1,[zero,zero,zero],nfft,ngfft,comm)
 end if

 ! Initialize the wave function descriptor.
 ! For the time being, no memory distribution, each node has the full set of states.
 my_minb = 1; my_maxb = mband

 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.True.; keep_ur=.False.

 ecut = dtset%ecut ! dtset%dilatmx
 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,dtset%paral_kgb,dummy_npw,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,ebands%istwfk,ebands%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd,header="Wavefunctions on the Fermi Surface",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)

 iomode = iomode_from_fname(wfk0_path)
 call wfd_read_wfk(wfd,wfk0_path,iomode)
 if (.False.) call wfd_test_ortho(wfd,cryst,pawtab,unit=std_out,mode_paral="PERS")

 ! k-point and bands where corrections are wanted
 ! These parameters will be passed via the input file (similar to kptgw, bdgw)
 kk = zero
 !kk = [0.5, 0.0, 0.0]
 sph_bstart = 1; sph_nb = 1
 ABI_CHECK(sph_nb <= ebands%mband, "sph_nb > ebands%mband")

 ABI_MALLOC(gkk_atm, (2, sph_nb, natom3))
 ABI_MALLOC(gkk_mu, (2, sph_nb, natom3))

 ! The k-point and the symmetries relating the BZ point to the IBZ.
 call listkk(dksqmax,cryst%gmet,indkk_k,ebands%kptns,kk,ebands%nkpt,1,cryst%nsym,&
    1,cryst%symafm,cryst%symrel,timrev1,use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(7a,es16.6,4a)' )&
    'The WFK file cannot be used to start thee present calculation ',ch10,&
    'It was asked that the wavefunctions be accurate, but',ch10,&
    'at least one of the k points could not be generated from a symmetrical one.',ch10,&
    'dksqmax=',dksqmax,ch10,&
    'Action: check your WFK file and k point input variables',ch10,&
    '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
    MSG_ERROR(msg)
 end if

 ik_ibz = indkk_k(1,1); isym_k = indkk_k(1,2)
 trev_k = indkk_k(1, 6); g0_k = indkk_k(1, 3:5)
 isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
 ABI_CHECK(isirr_k, "k-point must be in IBZ")
 kk_ibz = ebands%kptns(:,ik_ibz)
 npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

 ! ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! that will be used to symmetrize the wavefunctions in G-space.
 mpw = 0; gmax=0; cnt=0
 do iq_ibz=1,nqibz_k
   cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
   kq = kk + qibz_k(:,iq_ibz)

   ! TODO: g0 umklapp here can enter into play!
   ! fstab should contains the max of the umlapp G-vectors.
   ! gmax could not be large enough!
   call get_kg(kq,1,ecut,cryst%gmet,onpw,gtmp)
   mpw = max(mpw, onpw)
   do ipw=1,onpw
     do ii=1,3
      gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
     end do
   end do
   ABI_FREE(gtmp)
 end do

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, comm, ierr)
 call wrtout(std_out,sjoin('optimal value of mpw= ',itoa(mpw)),'COLL')

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4),work_ngfft(5),work_ngfft(6)))

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, npw_k))
 kg_k = wfd%kdata(ik_ibz)%kg_k
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm==1.
 ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

 ! TODO FOR PAW
 usecprj = 0
 ABI_DT_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1  ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2     ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnl1 = 0 ! gvnl1 is output
 ABI_MALLOC(gvnl1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nspden,natom,&
   dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
   usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 !PAW:allocate memory for non-symetrized 1st-order occupancies matrix (pawrhoij1)
 ! pawrhoij1_unsym => pawrhoij1
 ! if (psps%usepaw==1.and.iscf_mod>0) then
 !   if (paral_atom) then
 !     ABI_DATATYPE_ALLOCATE(pawrhoij1_unsym,(natom))
 !     cplex_rhoij=max(cplex,dtset%pawcpxocc);nspden_rhoij=nspden
 !     call pawrhoij_alloc(pawrhoij1_unsym,cplex_rhoij,nspden_rhoij,nspinor,&
 !       dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
 !   else
 !     pawrhoij1_unsym => pawrhoij1
 !     call pawrhoij_init_unpacked(pawrhoij1_unsym)
 !   end if
 ! end if

 ! Allocate vlocal. Note nvloc
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(tgam, (2,natom3,natom3,nsig))
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))
 !ABI_CALLOC(gvals_qibz, (2,natom3,natom3,nsig,gams%nqibz,nsppol))

 do spin=1,nsppol

   ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
   ! TODO: symmetrize them if kk is not irred
   ABI_MALLOC(kets_k, (2, npw_k*nspinor, sph_nb))
   do ib1=1,sph_nb
     band = ib1 + sph_bstart - 1
     call wfd_copy_cg(wfd, band, ik_ibz, spin, kets_k(1,1,ib1))
   end do

   ! Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal, &
     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   ! Integrations over q-points.
   do iq_ibz=1,nqibz_k
     !write(std_out,*)"in q-point loopq",iq_ibz,nqibz_k
     call cwtime(cpu,wall,gflops,"start")
     qpt = qibz_k(:,iq_ibz)
     tgam = zero

     ! Find the index of the q-point in the DVDB.
     db_iqpt = dvdb_findq(dvdb, qpt)

     ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials in the DVDB.
     if (db_iqpt /= -1) then
       if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
       ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
       ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
       call dvdb_readsym_allv1(dvdb, db_iqpt, cplex, nfftf, ngfftf, v1scf, xmpi_comm_self)
     else
       if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Could not find: ",ktoa(qpt), "in DVDB - interpolating"))
       ! Fourier interpolate of the potential
       ABI_CHECK(any(abs(qpt) > tol12), "qpt cannot be zero if Fourier interpolation is used")
       cplex = 2
       ABI_MALLOC(v1scf, (cplex,nfftf,nspden,natom3))
       call dvdb_ftinterp_qpt(dvdb, qpt, nfftf, ngfftf, v1scf, xmpi_comm_self)
     end if

     ! Examine the symmetries of the q wavevector
     !call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

     ! TODO: Make sure that symmetries in Q-space are preserved.
     ! Avoid fourq if qpt is in ddb

     ! Get phonon frequencies and displacements in reduced coordinates for this q-point
     call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

     ! Find k+q in the extended zone and extract symmetry info.
     ! Be careful here because there are two umklapp vectors to be considered:
     !
     !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
     !
     kq = kk + qpt

     call listkk(dksqmax,cryst%gmet,indkk_kq,ebands%kptns,kq,ebands%nkpt,1,cryst%nsym,&
       1,cryst%symafm,cryst%symrel,timrev1,use_symrec=.False.)

     if (dksqmax > tol12) then
       write(msg, '(7a,es16.6,4a)' )&
        'The WFK file cannot be used to start thee present calculation ',ch10,&
        'It was asked that the wavefunctions be accurate, but',ch10,&
        'at least one of the k points could not be generated from a symmetrical one.',ch10,&
        'dksqmax=',dksqmax,ch10,&
        'Action: check your WFK file and k point input variables',ch10,&
        '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
       MSG_ERROR(msg)
     end if

     ikq_ibz = indkk_kq(1,1); isym_kq = indkk_kq(1,2)
     trev_kq = indkk_kq(1, 6); g0_kq = indkk_kq(1, 3:5)
     isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
     kq_ibz = ebands%kptns(:,ikq_ibz)

     ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
     ! Be careful with time-reversal symmetry.
     if (isirr_kq) then
       ! Copy u_kq(G)
       !write(std_out,*)"before isirr"
       istwf_kq = wfd%istwfk(ikq_ibz); npw_kq = wfd%npwarr(ikq_ibz)
       ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
       kg_kq(:,1:npw_kq) = wfd%kdata(ikq_ibz)%kg_k
     else
       ! Will Reconstruct u_kq(G) from the IBZ image.
       !istwf_kq = set_istwfk(kq); if (.not. have_ktimerev) istwf_kq = 1
       !call change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat1,from_cg,to_cg)
       istwf_kq = 1
       call get_kg(kq,istwf_kq,ecut,cryst%gmet,npw_kq,gtmp)
       ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
       kg_kq(:,1:npw_kq) = gtmp(:,:npw_kq)
       ABI_FREE(gtmp)
     end if

     ! This array will store H1 |psi_nk>
     ABI_MALLOC(h1kets_kq, (2, npw_kq*nspinor, sph_nb))

     ! Allocate vlocal1 with correct cplex. Note nvloc
     ABI_STAT_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)
     ABI_CHECK(ierr==0, "oom vlocal1")

     ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
     do ipc=1,natom3
       call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
         pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))
     end do
     !write(std_out,*)"before loas spin"

     ! if PAW, one has to solve a generalized eigenproblem
     ! BE careful here because I will need sij_opt==-1
     gen_eigenpb = (psps%usepaw==1)
     sij_opt = 0; if (gen_eigenpb) sij_opt = 1
     ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

     ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
     !if (psps%useylm==1) then
     !   optder=0; if (useylmgr==1) optder=1
     !   call initylmg(cryst%gprimd,kg_k,kk,mkmem1,mpi_enreg,psps%mpsang,mpw,nband,mkmem1,&
     !     [npw_k],dtset%nsppol,optder,cryst%rprimd,ylm_k,ylmgr)
     !   call initylmg(cryst%gprimd,kg_kq,kq,mkmem1,mpi_enreg,psps%mpsang,mpw,nband,mkmem1,&
     !     [npw_kq],dtset%nsppol,optder,cryst%rprimd,ylm_kq,ylmgr_kq)
     !end if

     ! Loop over all 3*natom perturbations.
     do ipc=1,natom3
       idir = mod(ipc-1, 3) + 1
       ipert = (ipc - idir) / 3 + 1

       ! Prepare application of the NL part.
       !write(std_out,*)"before init_rf spin"
       call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=1)
       call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc), &
         comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
       !write(std_out,*)"after load_rf spin"

       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&  ! In
         cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&          ! In
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&         ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)     ! Out
       !write(std_out,*)"after getgh1c_setup"

       ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,sph_nb
         band = ib2 + sph_bstart - 1
         eig0nk = ebands%eig(band,ik_ibz,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,copt0,kets_k(:,:,ib2),cwaveprj0,h1kets_kq(:,:,ib2),&
           grad_berry,gs1c,gs_hamkq,gvnl1,idir,ipert,eshift,mpi_enreg,optlocal,&
           optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

         ! TODO: Solve Sternheimer equations non-self-consistently (??)
       end do

       call destroy_rf_hamiltonian(rf_hamkq)

       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
       if (allocated(ph3d1)) then
         ABI_FREE(ph3d1)
       end if

     end do ! ipc (loop over 3*natom atomic perturbations)

     ABI_FREE(gs1c)

     ! Sum over bands (all bands are taken for the time being)
     istwf_kirr = wfd%istwfk(ikq_ibz); npw_kirr = wfd%npwarr(ikq_ibz)
     ABI_MALLOC(bras_kq, (2, npw_kq*nspinor, 1))
     ABI_MALLOC(cgwork, (2, npw_kirr*nspinor))

     do ib1=1,wfd%nband(ikq_ibz, spin)

       ! symmetrize wavefunctions from IBZ (if needed).
       ! Be careful with time-reversal symmetry.
       if (isirr_kq) then
         ! Copy u_kq(G)
         call wfd_copy_cg(wfd, ib1, ikq_ibz, spin, bras_kq(1,1,1))
       else
         ! Reconstruct u_kq(G) from the IBZ image.
         !istwf_kq = set_istwfk(kq); if (.not. have_ktimerev) istwf_kq = 1
         !call change_istwfk(from_npw,from_kg,from_istwfk,to_npw,to_kg,to_istwfk,n1,n2,n3,ndat1,from_cg,to_cg)

         ! Use cgwork as workspace array, results stored in bras_kq
         !g0_kq =  g0ibz_kq + g0bz_kq
         call wfd_copy_cg(wfd, ib1, ikq_ibz, spin, cgwork)
         call cg_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1,&
                        npw_kirr, wfd%kdata(ikq_ibz)%kg_k,&
                        npw_kq, kg_kq, istwf_kirr, istwf_kq, cgwork, bras_kq(:,:,1), work_ngfft, work)
       end if

       do ipc=1,natom3
         ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
         !The array eig1_k contains:
         !
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
         do ib2=1,sph_nb
           call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras_kq(1,1,1),h1kets_kq(1,1,ib2),&
                mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           gkk_atm(:,ib2,ipc) = [dotr, doti]
         end do
       end do

       call gkkmu_from_atm(1, sph_nb, 1, natom, gkk_atm, phfrq, displ_red, gkk_mu, num_smallw)

#if 0
       do mu
         ! sum contribution to phonon self-energy
         do ib2=1,mband
           f_nk = ebands_k%occ(ib2,ik,spin)
           eig0nk = ebands_k%eig(ib2,ik,spin)
           do ib1=1,mband_kq
             f_mkq = ebands_kq%occ(ib1,ikq,spin)
             eig0mkq = ebands_kq%eig(ib1,ikq,spin)
             if (abs(f_mkq - f_nk) <= tol12) cycle

             gkk2 = gkk_mu(1,ib1,ib2) ** 2 + gkk_mu(2,ib1,ib2) ** 2
             term1 = (f_mkq - f_nk) * (eig0mkq - eig0nk - omega) / ((eig0mkq - eig0nk - omega) ** 2 + eta ** 2)
             term2 = (f_mkq - f_nk) * (eig0mkq - eig0nk        ) / ((eig0mkq - eig0nk        ) ** 2 + eta ** 2)
             !Pi_ph(imode) = Pi_ph(imode) + wtk * gkk2 * (term1 - term2)
           end do
         end do
       end do  ! mu
#endif

     end do ! ib1 (sum over bands)

     ABI_FREE(bras_kq)
     ABI_FREE(cgwork)

     !call xmpi_sum(tgam, comm, ierr)
     !do isig=1,nsig
     !  ! Save results for this (q-point, spin)
     !  !write(std_out,*)tgam(:,:,:,isig)
     !  gvals_qibz(:,:,:,isig,iq_ibz,spin) = tgam(:,:,:,isig)
     !end do ! isig

     ABI_FREE(h1kets_kq)
     ABI_FREE(v1scf)
     ABI_FREE(vlocal1)

     call cwtime(cpu,wall,gflops,"stop")
     write(msg,'(a,i0,2(a,f8.2))')"q-point [",iq_ibz,"] completed. cpu:",cpu,", wall:",wall
     call wrtout(std_out,msg,"COLL",do_flush=.True.)
   end do ! iq_ibz

   ABI_FREE(kets_k)
 end do ! spin

 ! Collect gvals_qibz on each node and divide by the total number of q-points in the full mesh.
 !do spin=1,nsppol
 !  gvals_qibz(:,:,:,:,:,spin) = gvals_qibz(:,:,:,:,:,spin) / nqbz
 !end do
 call wrtout(std_out, "Computation of tgamma matrices completed", "COLL", do_flush=.True.)

 ! Save results to file
 ncid = nctk_noid
#ifdef HAVE_NETCDF
 ! Open the netcdf file used to store the results of the calculation.
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_SIGMAPH.nc"), xmpi_comm_self))
   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 ! Free memory
 ABI_FREE(tmesh)
 ABI_FREE(gkk_atm)
 ABI_FREE(gkk_mu)
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(tgam)
 ABI_FREE(qibz_k)

 call destroy_hamiltonian(gs_hamkq)
 call wfd_free(wfd)

 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)

end subroutine sigmaph_driver
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/gkkmu_from_atm
!! NAME
!!  gkkmu_from_atm
!!
!! FUNCTION
!!  Transform the gkk from (atom, red_direction) basis to phonon-mode basis
!!
!! INPUTS
!!  nb1,nb2=Number of bands in gkk_atm matrix.
!!  nk=Number of k-points (usually 1)
!!  natom=Number of atoms.
!!  gkk_atm(2,nb1,nb2,3*natom)=EPH matrix elements in the atomic basis.
!!  phfrq(3*natom)=Phonon frequencies in Ha
!!  displ_red(2,3*natom,3*natom)=Phonon displacement in reduded coordinates.
!!
!! OUTPUT
!!  gkk_mu(2,nb1,nb2,3*natom)=EPH matrix elements in the phonon-mode basis.
!!  num_smallw=Number of negative/too small frequencies that have been ignored
!!    by setting the corresponding gkk_mu to zero.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gkkmu_from_atm(nb1, nb2, nk, natom, gkk_atm, phfrq, displ_red, gkk_mu, num_smallw)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gkkmu_from_atm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nb1, nb2, nk, natom
 integer,intent(out) :: num_smallw
!arrays
 real(dp),intent(in) :: phfrq(3*natom),displ_red(2,3*natom,3*natom)
 real(dp),intent(in) :: gkk_atm(2,nb1,nb2,nk,3*natom)
 real(dp),intent(out) :: gkk_mu(2,nb1,nb2,nk,3*natom)

!Local variables-------------------------
!scalars
 integer :: mu,ipc

! *************************************************************************

 gkk_mu = zero; num_smallw = 0

 ! Loop over phonon branches.
 do mu=1,3*natom
   ! Ignore negative or too small frequencies
   if (phfrq(mu) < tol6) then
     num_smallw = num_smallw + 1; cycle
   end if

   ! Transform the gkk from atom,red_dir basis to phonon mode basis
   do ipc=1,3*natom
     gkk_mu(1,:,:,:,mu) = gkk_mu(1,:,:,:,mu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(1,ipc,mu) &
       - gkk_atm(2,:,:,:,ipc) * displ_red(2,ipc,mu)
     gkk_mu(2,:,:,:,mu) = gkk_mu(2,:,:,:,mu) &
       + gkk_atm(1,:,:,:,ipc) * displ_red(2,ipc,mu) &
       + gkk_atm(2,:,:,:,ipc) * displ_red(1,ipc,mu)
   end do

   gkk_mu(:,:,:,:,mu) = gkk_mu(:,:,:,:,mu) / sqrt(two * phfrq(mu))
 end do

end subroutine gkkmu_from_atm
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/nfd
!! NAME
!!  nfd
!!
!! FUNCTION
!!   Fermi–Dirac statistic
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   mu=Chemical potential in Ha.
!!   kT=Value of K_Boltzmann x T in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function nfd(ee, mu, kT)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nfd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, mu, kT

! *************************************************************************

 !TODO: Find decent value.
 if (kT > tol16) then
   nfd = one / (exp((ee-mu) / kT) + one)
 else
   ! Heaviside
   if (ee > mu) then
     nfd = zero
   else if (ee < mu) then
     nfd = one
   else
     nfd = half
   end if
 end if

end function nfd
!!***

!----------------------------------------------------------------------

!!****f* m_sigmaph/nbe
!! NAME
!!  nbe
!!
!! FUNCTION
!!   Bose-Einstein statistic
!!
!! INPUTS
!!   ee=Single particle energy in Ha
!!   mu=Chemical potential in Ha.
!!   kT=Value of K_Boltzmann x T in Ha.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

elemental real(dp) function nbe(ee, mu, kT)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nbe'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp),intent(in) :: ee, mu, kT

! *************************************************************************

 if (kT > tol16) then
   nbe = one / (exp((ee-mu) / kT) - one)
 else
   nbe = zero
 end if

end function nbe
!!***

end module m_sigmaph
!!***
