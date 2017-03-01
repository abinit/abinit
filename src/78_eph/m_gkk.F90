!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gkk
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of electron-phonon coupling matrix elements (gkk)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2017 ABINIT group (GKA)
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


module m_gkk

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ebands
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,           only : cwtime
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_fftcore,        only : ngfft_seq
 use defs_datatypes,   only : ebands_t
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_bz_mesh,        only : findqg0

 implicit none

 private
!!***

 public :: eph_gkk


contains  !=========================================================================================================================
!!***

!!****f* m_gkk/eph_gkk
!! NAME
!!  eph_gkk
!!
!! FUNCTION
!!  Compute electron-phonon coupling matrix elements.
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
!! NOTES
!!
!! CHILDREN
!!      appdig,cwtime,destroy_hamiltonian,destroy_rf_hamiltonian,dotprod_g
!!      dvdb_ftinterp_qpt,dvdb_open_read,dvdb_readsym_allv1,findqg0,get_kg
!!      getgh1c,getgh1c_setup,getph,gkk_free,gkk_init,gkk_ncwrite
!!      init_hamiltonian,init_rf_hamiltonian,littlegroup_q
!!      load_spin_hamiltonian,load_spin_rf_hamiltonian,pawcprj_free
!!      rf_transgrid_and_pack,wfd_copy_cg,wfd_free,wfd_init,wfd_print
!!      wfd_read_wfk,wfd_test_ortho,wrtout,xmpi_split_work,xmpi_sum_master
!!
!! SOURCE

subroutine eph_gkk(wfk0_path,wfq_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands_k,ebands_kq,dvdb,ifc,&
                       pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wfk
 use m_ddb
 use m_dvdb
 use m_ifc
 use m_fft
 use m_hamiltonian
 use m_pawcprj

 use m_time,            only : sec2str
 use m_fstrings,        only : sjoin, itoa, ftoa, ktoa
 use m_io_tools,        only : iomode_from_fname
 use m_cgtools,         only : dotprod_g
 use m_fftcore,         only : get_kg, kpgsph, sphere
 use m_crystal,         only : crystal_t
 use m_wfd,             only : wfd_init, wfd_free, wfd_print, wfd_t, wfd_test_ortho, wfd_copy_cg,&
                               wfd_read_wfk, wfd_wave_free, wfd_rotate, wfd_reset_ur_cprj, wfd_get_ur
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_pawfgr,          only : pawfgr_type
 use m_eig2d,           only : gkk_t, gkk_init, gkk_ncwrite,gkk_free
! use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
! use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
! use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
! use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
! use m_pawdij,          only : pawdij, symdij
! use m_pawcprj,         only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph_gkk'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path, wfq_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands_k, ebands_kq
 type(dvdb_t),target,intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: dummy_npw=1,tim_getgh1c=1,berryopt0=0
 integer,parameter :: useylmgr1=0,master=0
 integer :: my_rank,nproc,iomode,mband,mband_kq,my_minb,my_maxb,nsppol,nkpt,nkpt_kq,idir,ipert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,onpw
 integer :: ib1,ib2,band
 integer :: ik,ikq,timerev_q
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw,my_mpw,ierr,my_kstart,my_kstop,cnt,ncid
 integer :: n1,n2,n3,n4,n5,n6,nspden
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnl1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1
 real(dp) :: cpu,wall,gflops
 real(dp) :: ecut,eshift,eig0nk,dotr,doti
 logical :: i_am_master, gen_eigenpb
 type(wfd_t) :: wfd_k, wfd_kq
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 character(len=500) :: msg
 character(len=fnlen) :: fname, gkkfilnam
!arrays
 integer :: g0_k(3),symq(4,2,cryst%nsym),dummy_gvec(3,dummy_npw)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:),nband_kq(:,:),blkflg(:,:)
 real(dp) :: kk(3),kq(3),qpt(3)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),gkk(:,:,:,:,:)
 real(dp),allocatable :: bras(:,:,:),kets(:,:,:),h1_kets(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnl1(:,:)
 real(dp),allocatable ::  gs1c(:,:)
 logical,allocatable :: bks_mask(:,:,:),bks_mask_kq(:,:,:),keep_ur(:,:,:),keep_ur_kq(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)
 type(gkk_t)     :: gkk2d

!************************************************************************

 write(msg, '(2a)') "Computation of electron-phonon coupling matrix elements (gkk)", ch10
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 i_am_master = (my_rank == master)

 ! Copy important dimensions
 natom = cryst%natom
 natom3 = 3 * natom
 nsppol = ebands_k%nsppol
 nspinor = ebands_k%nspinor
 nspden = dtset%nspden
 nkpt = ebands_k%nkpt
 mband = ebands_k%mband
 nkpt_kq = ebands_kq%nkpt
 mband_kq = ebands_kq%mband
 ecut = dtset%ecut

! GKA TODO: Make sure there is a single q-point present.
 qpt = dtset%qptn(:)

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)


 ! Open the DVDB file
 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)


 ! Initialize the wave function descriptors.
 ! For the time being, no memory distribution, each node has the full set of states.
 my_minb = 1; my_maxb = mband

 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.False.; keep_ur=.False.

 ABI_MALLOC(nband_kq, (nkpt_kq, nsppol))
 ABI_MALLOC(bks_mask_kq,(mband_kq, nkpt_kq, nsppol))
 ABI_MALLOC(keep_ur_kq,(mband_kq, nkpt_kq ,nsppol))
 nband_kq=mband_kq; bks_mask_kq=.False.; keep_ur_kq=.False.

 ! Distribute the k-points over the processors
 call xmpi_split_work(nkpt,comm,my_kstart,my_kstop,msg,ierr)
 do ik=1,nkpt
 if (.not. ((ik .ge. my_kstart) .and. (ik .le. my_kstop))) cycle

   kk = ebands_k%kptns(:,ik)
   kq = kk + qpt
   call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:),(/1,1,1/))  ! Find the index of the k+q point

   bks_mask(:,ik,:) = .True.
   bks_mask_kq(:,ikq,:) = .True.
 end do

 ! Initialize the wavefunction descriptors
 call wfd_init(wfd_k,cryst,pawtab,psps,keep_ur,dtset%paral_kgb,dummy_npw,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,ebands_k%istwfk,ebands_k%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd_k,header="Wavefunctions on the k-points grid",mode_paral='PERS')

 call wfd_init(wfd_kq,cryst,pawtab,psps,keep_ur_kq,dtset%paral_kgb,dummy_npw,mband_kq,nband_kq,nkpt_kq,nsppol,bks_mask_kq,&
   nspden,nspinor,dtset%ecutsm,dtset%dilatmx,ebands_kq%istwfk,ebands_kq%kptns,ngfft,&
   dummy_gvec,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm,opt_ecut=ecut)

 call wfd_print(wfd_kq,header="Wavefunctions on the q-shifted k-points grid",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband_kq)
 ABI_FREE(bks_mask_kq)
 ABI_FREE(keep_ur_kq)

 ! Read wafefunctions on the k-points grid and q-shifted k-points grid.
 iomode = iomode_from_fname(wfk0_path)
 call wfd_read_wfk(wfd_k,wfk0_path,iomode)
 if (.False.) call wfd_test_ortho(wfd_k,cryst,pawtab,unit=std_out,mode_paral="PERS")

 iomode = iomode_from_fname(wfq_path)
 call wfd_read_wfk(wfd_kq,wfq_path,iomode)
 if (.False.) call wfd_test_ortho(wfd_kq,cryst,pawtab,unit=std_out,mode_paral="PERS")

 ! ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! Find the appropriate value of mpw
 mpw = 0; cnt=0
 do spin=1,nsppol
   do ik=1,nkpt
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kk = ebands_k%kptns(:,ik)
     call get_kg(kk,1,ecut,cryst%gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 cnt=0
 do spin=1,nsppol
   do ikq=1,nkpt_kq
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kq = ebands_kq%kptns(:,ikq)
     call get_kg(kq,1,ecut,cryst%gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)


 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
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

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,natom,&
&  dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
&  usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda,&
&  comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)

! Allocate vlocal. Note nvloc
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))

 ! Allocate work space arrays.
 ABI_MALLOC(blkflg, (natom3,natom3))
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))


 call cwtime(cpu,wall,gflops,"start")

 ! Find the index of the q-point in the DVDB.
 db_iqpt = dvdb_findq(dvdb, qpt)

 if (db_iqpt /= -1) then
   if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
   ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
   ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
   call dvdb_readsym_allv1(dvdb, db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
 else
   if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Could not find: ",ktoa(qpt), "in DVDB - interpolating"))
   ! Fourier interpolate of the potential
   ABI_CHECK(any(abs(qpt) > tol12), "qpt cannot be zero if Fourier interpolation is used")
   cplex = 2
   call dvdb_ftinterp_setup(dvdb,ifc%ngqpt,ifc%nqshft,ifc%qshft,nfft,ngfft,comm,cryst)
   ABI_MALLOC(v1scf, (cplex,nfftf,nspden,natom3))
   call dvdb_ftinterp_qpt(dvdb, qpt, nfftf, ngfftf, v1scf, comm)
 end if

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

 ! Allocate vlocal1 with correct cplex. Note nvloc
 ABI_STAT_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)
 ABI_CHECK(ierr==0, "oom vlocal1")

 ABI_MALLOC(gkk, (2*mband*nsppol,nkpt,1,1,mband_kq))

 ! Loop over all 3*natom perturbations.
 do ipc=1,natom3
   idir = mod(ipc-1, 3) + 1
   ipert = (ipc - idir) / 3 + 1

   write(msg, '(a,2i4)')  "Treating ipert, idir = ", ipert, idir
   call wrtout(std_out, msg, "COLL", do_flush=.True.)

   gkk = zero

   do spin=1,nsppol

     ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
     !do ipc=1,natom3
     call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
               pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))
     !end do

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ABI_MALLOC(bras, (2, mpw*nspinor, mband))
     ABI_MALLOC(kets, (2, mpw*nspinor, mband))
     ABI_MALLOC(h1_kets, (2, mpw*nspinor, mband))


     ! GKA : This little block used to be right after the perturbation loop
     ! Prepare application of the NL part.
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=.true.)
               ! paw_ij1=paw_ij1,comm_atom=mpi_enreg%comm_atom,&
               !&mpi_atmtab=mpi_enreg%my_atmtab,my_spintab=mpi_enreg%my_isppoltab)
     call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

     do ik=1,nkpt

       ! Only do a subset a k-points
       if (.not. ((ik .ge. my_kstart) .and. (ik .le. my_kstop))) cycle

       kk = ebands_k%kptns(:,ik)

       kq = kk + qpt
       call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:),(/1,1,1/))  ! Find the index of the k+q point

       ! Copy u_k(G)
       istwf_k = wfd_k%istwfk(ik); npw_k = wfd_k%npwarr(ik)
       ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
       kg_k(:,1:npw_k) = wfd_k%kdata(ik)%kg_k
       do ib2=1,mband
         call wfd_copy_cg(wfd_k, ib2, ik, spin, kets(1,1,ib2))
       end do

       ! Copy u_kq(G)
       istwf_kq = wfd_kq%istwfk(ikq); npw_kq = wfd_kq%npwarr(ikq)
       ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
       kg_kq(:,1:npw_kq) = wfd_kq%kdata(ikq)%kg_k
       do ib1=1,mband_kq
         call wfd_copy_cg(wfd_kq, ib1, ikq, spin, bras(1,1,ib1))
       end do

       ! if PAW, one has to solve a generalized eigenproblem
       ! BE careful here because I will need sij_opt==-1
       gen_eigenpb = (psps%usepaw==1)
       sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

       ! GKA : Previous loop on 3*natom perturbations used to start here

       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&                   ! In
         cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&                           ! In
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&                          ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)                      ! Out

       ! Calculate dvscf * psi_k, results stored in h1_kets on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,mband
         eig0nk = ebands_k%eig(ib2,ik,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,kets(:,:,ib2),cwaveprj0,h1_kets(:,:,ib2),&
&                     grad_berry,gs1c,gs_hamkq,gvnl1,idir,ipert,eshift,mpi_enreg,optlocal,&
&                     optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
       end do

       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
       ABI_FREE(gs1c)

       if (allocated(ph3d1)) then
         ABI_FREE(ph3d1)
       end if

       ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
       !The array eig1_k contains:
       !
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
       do ib2=1,mband
         do ib1=1,mband_kq
           call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras(1,1,ib1),h1_kets(1,1,ib2),&
             mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           band = 2*ib2-1 + (spin-1) * 2 * mband
           gkk(band,ik,1,1,ib1) = dotr
           gkk(band+1,ik,1,1,ib1) = doti
         end do
       end do

     end do ! ikfs

     call destroy_rf_hamiltonian(rf_hamkq)

     ABI_FREE(bras)
     ABI_FREE(kets)
     ABI_FREE(h1_kets)

   end do ! spin

   ! Gather the k-points computed by all processes
   call xmpi_sum_master(gkk,master,comm,ierr)

   ! Init a gkk_t object
   call gkk_init(gkk,gkk2d,mband,nsppol,nkpt,1,1)

   ! Write the netCDF file.
   call appdig(ipc,dtfil%fnameabo_gkk,gkkfilnam)
   fname = strcat(gkkfilnam,".nc")
#ifdef HAVE_NETCDF
   if (i_am_master) then
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
     NCF_CHECK(crystal_ncwrite(cryst, ncid))
     NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
     call gkk_ncwrite(gkk2d,qpt,1.0_dp, ncid)
     NCF_CHECK(nf90_close(ncid))
   end if
#endif

   ! Free memory
   call gkk_free(gkk2d)

 end do ! ipc (loop over 3*natom atomic perturbations)

 call cwtime(cpu,wall,gflops,"stop")

 write(msg, '(2a)') "Computation of electron-phonon coupling matrix elements (gkk) completed", ch10
 call wrtout(ab_out, msg, "COLL", do_flush=.True.)
 call wrtout(std_out, msg, "COLL", do_flush=.True.)

 ! Free memory
 ABI_FREE(gkk)
 ABI_FREE(v1scf)
 ABI_FREE(vlocal1)
 ABI_FREE(gvnl1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(blkflg)

 call destroy_hamiltonian(gs_hamkq)
 call wfd_free(wfd_k)
 call wfd_free(wfd_kq)

 call pawcprj_free(cwaveprj0)
 ABI_DT_FREE(cwaveprj0)

end subroutine eph_gkk
!!***

!----------------------------------------------------------------------

end module m_gkk
!!***
