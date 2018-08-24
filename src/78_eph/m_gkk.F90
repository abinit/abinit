!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gkk
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of electron-phonon coupling matrix elements (gkk)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (GKA, MG)
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
 use m_abicore
 use m_xmpi
 use m_errors
 use m_ifc
 use m_ebands
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfk
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,           only : cwtime, sec2str
 use m_io_tools,       only : iomode_from_fname
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_symtk,          only : littlegroup_q
 use m_fftcore,        only : ngfft_seq, get_kg, kpgsph, sphere
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_crystal,        only : crystal_t
 use m_crystal_io,     only : crystal_ncwrite
 use m_bz_mesh,        only : findqg0
 use m_cgtools,        only : dotprod_g
 use m_kg,             only : getph
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_eig2d,          only : gkk_t, gkk_init, gkk_ncwrite, gkk_free
 use m_wfd,            only : wfd_init, wfd_free, wfd_print, wfd_t, wfd_test_ortho, wfd_copy_cg,&
                              wfd_read_wfk, wfd_wave_free, wfd_rotate, wfd_reset_ur_cprj, wfd_get_ur
 use m_getgh1c,          only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_fourier_interpol, only : transgrid
! use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
! use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
! use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
! use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
! use m_pawdij,          only : pawdij, symdij
! use m_pawcprj,         only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy

 implicit none

 private
!!***

 public :: eph_gkk
 public :: ncwrite_v1qnu          ! Compute \delta V_{q,nu)(r) and dump results to netcdf file.

contains  !===========================================================================
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
!!      get_kg
!!
!! SOURCE

subroutine eph_gkk(wfk0_path,wfq_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands_k,ebands_kq,dvdb,ifc,&
                       pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eph_gkk'
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
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor
 integer :: ib1,ib2,band
 integer :: ik,ikq,timerev_q
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: mpw,mpw_k,mpw_kq,ierr,my_kstart,my_kstop,ncid
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
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:),nband_kq(:,:),blkflg(:,:)
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
 cplex = 2
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
 call find_mpw(mpw_k, ebands_k%kptns(:,:), nsppol, nkpt, cryst%gmet,ecut,comm)
 call find_mpw(mpw_kq, ebands_kq%kptns(:,:), nsppol, nkpt_kq, cryst%gmet,ecut,comm)
 mpw = max(mpw_k, mpw_kq)

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

   ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
   call dvdb_interpolate_v1scf(dvdb,cryst,qpt,ifc%ngqpt,ifc%nqshft,ifc%qshft, &
   &                           nfft, ngfft, nfftf, ngfftf, v1scf, comm)

 end if

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

 ! Allocate vlocal1 with correct cplex. Note nvloc
 ABI_STAT_MALLOC(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)
 ABI_CHECK(ierr==0, "oom vlocal1")

 ABI_MALLOC(gkk, (2*mband*nsppol,nkpt,1,1,mband_kq))

 ! ========================================================================== !
 ! Begin loop on perturbations, spins and k-points
 ! ========================================================================== !

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
     call load_spin_rf_hamiltonian(rf_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

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

 ! ========================================================================== !
 ! Free memory
 ! ========================================================================== !

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

!!****f* m_gkk/ncwrite_v1qnu
!! NAME
!!  ncwrite_v1qnu
!!
!! FUNCTION
!!  Compute \delta V_{q,nu)(r) and dump results to netcdf file.
!!  This routine should be called by a single processor.
!!
!! INPUT
!!  dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!!  cryst(crystal_t)=Crystalline structure
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  nqlist=Number of q-points
!!  qlist(3,nqlist)=List of q-points where \delta V_{q,nu)(r) is wanted.
!!    Potentials will be Fourier interpolated if qpt is not in DVDB.
!!  prtvol=Verbosity level
!!  path=Name of the netcdf file.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ncwrite_v1qnu(dvdb, cryst, ifc, nqlist, qlist, prtvol, path)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ncwrite_v1qnu'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqlist,prtvol
 type(dvdb_t),intent(inout) :: dvdb
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 character(len=*),intent(in) :: path
!arrays
 real(dp),intent(in) :: qlist(3,nqlist)

!Local variables-------------------------------
!scalars
 integer :: ip,nu,iq,db_iqpt,do_ftv1q,cplex,ispden,nfftf,comm
#ifdef HAVE_NETCDF
 integer :: ncid,ncerr
#endif
 character(len=500) :: msg
!arrays
 integer :: ngfftf(18)
 real(dp) :: phfrq(3*cryst%natom),qpt(3)
 real(dp) :: displ_cart(2,3*cryst%natom,3*cryst%natom),displ_red(2,3*cryst%natom,3*cryst%natom)
 real(dp),allocatable :: v1scf(:,:,:,:),v1qnu(:,:,:,:)

!************************************************************************

 call wrtout(std_out, sjoin("Writing \delta V_{q,nu)(r) potentials to file:", path), do_flush=.True.)

 comm = xmpi_comm_self
 call ngfft_seq(ngfftf, dvdb%ngfft3_v1(:,1)); nfftf = product(ngfftf(1:3))

 call dvdb_open_read(dvdb, ngfftf, xmpi_comm_self)

 ! Create netcdf file.
#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_open_create(ncid, path, comm))
 NCF_CHECK(crystal_ncwrite(cryst, ncid))

 ! Add other dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("nfftf", nfftf), nctkdim_t("nspden", dvdb%nspden), &
   nctkdim_t("natom3", 3 * cryst%natom), nctkdim_t("nqlist", nqlist)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Define arrays
 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("ngfftf", "int", "three"), &
   nctkarr_t("qlist", "dp", "three, nqlist"), &
   nctkarr_t("wqlist", "dp", "natom3, nqlist"), &
   nctkarr_t("displ_cart", "dp", "two, natom3, natom3, nqlist"), &
   nctkarr_t("v1qnu", "dp", "two, nfftf, nspden, natom3, nqlist")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngfftf"), ngfftf(1:3)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qlist"), qlist))
#endif

 ! Activate Fourier interpolation if q-points are not in the DVDB file.
 ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials already available in the DVDB.
 ! without performing FT interpolation.
 do_ftv1q = 0
 do iq=1,nqlist
   if (dvdb_findq(dvdb, qlist(:,iq)) == -1) do_ftv1q = do_ftv1q + 1
 end do
 if (do_ftv1q /= 0) then
   write(msg, "(2(a,i0),a)")"Will use Fourier interpolation of DFPT potentials [",do_ftv1q,"/",nqlist,"]"
   call wrtout(std_out, msg)
   call dvdb_ftinterp_setup(dvdb, ifc%ngqpt, 1, [zero,zero,zero], nfftf, ngfftf, comm)
 end if

 ABI_MALLOC(v1qnu, (2, nfftf, dvdb%nspden, dvdb%natom3))

 do iq=1,nqlist
   qpt = qlist(:,iq)
   call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb_findq(dvdb, qpt)

   ! TODO: handle q_bz = S q_ibz case by symmetrizing the potentials already available in the DVDB.
   if (db_iqpt /= -1) then
     if (prtvol > 0) call wrtout(std_out, sjoin("Found:", ktoa(qpt), "in DVDB with index", itoa(db_iqpt)))
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
     call dvdb_readsym_allv1(dvdb, db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
   else
     if (prtvol > 0) call wrtout(std_out, sjoin("Could not find:", ktoa(qpt), "in DVDB - interpolating"))
     ! Fourier interpolation of the potential
     ABI_CHECK(any(abs(qpt) > tol12), "qpt cannot be zero if Fourier interpolation is used")
     cplex = 2
     ABI_MALLOC(v1scf, (cplex,nfftf, dvdb%nspden, dvdb%natom3))
     call dvdb_ftinterp_qpt(dvdb, qpt, nfftf, ngfftf, v1scf, comm)
   end if

   do nu=1,dvdb%natom3
     ! v1qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
     ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
     !v1qnu(2, nfftf, nspden, natom3), v1scf(cplex, nfftf, nspden, natom3)
     v1qnu(:, :, :, nu) = zero
     do ip=1,dvdb%natom3
       do ispden=1,dvdb%nspden
         if (cplex == 2) then
           v1qnu(1, :, ispden, nu) = v1qnu(1, :, ispden, nu) + &
             displ_red(1,ip,nu) * v1scf(1,:,ispden,ip) - displ_red(2,ip,nu) * v1scf(2,:,ispden,ip)
           v1qnu(2, :, ispden, nu) = v1qnu(2, :, ispden, nu) + &
             displ_red(2,ip,nu) * v1scf(1,:,ispden,ip) + displ_red(1,ip,nu) * v1scf(2,:,ispden,ip)
         else
           ! Gamma point. d(q) = d(-q)* --> d is real.
           v1qnu(1, :, ispden, nu) = v1qnu(1, :, ispden, nu) + &
             displ_red(1,ip,nu) * v1scf(1,:,ispden,ip)
         end if
       end do
     end do
   end do
   ! Save results to file.
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wqlist"), phfrq, start=[1,iq]))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "displ_cart"), displ_cart, start=[1,1,1,iq]))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1qnu"), v1qnu, start=[1,1,1,1,iq]))
#endif
   ABI_FREE(v1scf)
 end do

 ABI_FREE(v1qnu)
#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_close(ncid))
#endif
 call dvdb_close(dvdb)

 call wrtout(std_out, "dvqnu file written", do_flush=.True.)

end subroutine ncwrite_v1qnu
!!***

!----------------------------------------------------------------------

!!****f* m_gkk/find_mpw
!! NAME
!!  find_mpw
!!
!! FUNCTION
!!  Look at all k-points and spins to find the maximum
!!  number of plane waves.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gkk
!!
!! NOTES
!!
!! CHILDREN
!!      get_kg
!!
!! SOURCE

subroutine find_mpw(mpw, kpts, nsppol, nkpt, gmet, ecut, comm)

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_fftcore,         only : get_kg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_mpw'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: mpw
 integer,intent(in) :: nsppol, nkpt
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: my_rank, cnt, nproc, ierr
 integer :: ispin, ikpt
 integer :: my_mpw, onpw
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kpt(3)

 my_rank = xmpi_comm_rank(comm)
 nproc = xmpi_comm_size(comm)

 mpw = 0; cnt=0
 do ispin=1,nsppol
   do ikpt=1,nkpt
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kpt = kpts(:,ikpt)
     call get_kg(kpt,1,ecut,gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)

end subroutine find_mpw
!!***

end module m_gkk
!!***
