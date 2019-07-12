!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gkk
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of electron-phonon coupling matrix elements (gkk)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (GKA, MG)
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
 use m_fftcore,        only : ngfft_seq, get_kg
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : findqg0
 use m_cgtools,        only : dotprod_g
 use m_kg,             only : getph
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_eig2d,          only : gkk_t, gkk_init, gkk_ncwrite, gkk_free
 use m_wfd,            only : wfd_init, wfd_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup

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
 integer,parameter :: tim_getgh1c=1,berryopt0=0, useylmgr1=0,master=0
 integer :: my_rank,nproc,mband,mband_kq,my_minb,my_maxb,nsppol,nkpt,nkpt_kq,idir,ipert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor
 integer :: ib1,ib2,band,ik,ikq,timerev_q
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq, comm_rpt, method
 integer :: mpw,mpw_k,mpw_kq,ierr,my_kstart,my_kstop,ncid
 integer :: n1,n2,n3,n4,n5,n6,nspden,ncerr
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1, interpolated
 real(dp) :: cpu,wall,gflops,ecut,eshift,eig0nk,dotr,doti
 logical :: i_am_master, gen_eigenpb
 type(wfd_t) :: wfd_k, wfd_kq
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(gkk_t) :: gkk2d
 character(len=500) :: msg, what
 character(len=fnlen) :: fname, gkkfilnam, wr_path
!arrays
 integer :: g0_k(3),symq(4,2,cryst%nsym)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:),nband_kq(:,:),blkflg(:,:), wfd_istwfk(:)
 real(dp) :: kk(3),kq(3),qpt(3),phfrq(3*cryst%natom)
 real(dp),allocatable :: displ_cart(:,:,:),displ_red(:,:,:), eigens_kq(:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),gkk(:,:,:,:,:)
 real(dp),allocatable :: bras(:,:,:),kets(:,:,:),h1_kets(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnlx1(:,:), gs1c(:,:), gkq_atm(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),bks_mask_kq(:,:,:),keep_ur(:,:,:),keep_ur_kq(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)

!************************************************************************

 what = "(GKK files)"; if (dtset%eph_task == -2) what = "GKQ file"
 write(msg, '(3a)') " Computation of electron-phonon coupling matrix elements ", trim(what), ch10
 call wrtout([std_out, ab_out], msg, do_flush=.True.)

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm); i_am_master = my_rank == master

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
 !write(std_out, *)"ebands dims (b, k, s): ", ebands_k%mband, ebands_k%nkpt, ebands_k%nsppol
 !write(std_out, *)"ebands_kq dims (b, k, s): ", ebands_kq%mband, ebands_kq%nkpt, ebands_kq%nsppol

 qpt = dtset%qptn(:)

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

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
 call xmpi_split_work(nkpt,comm,my_kstart,my_kstop)
 do ik=1,nkpt
   if (.not. (ik >= my_kstart .and. ik <= my_kstop)) cycle
   kk = ebands_k%kptns(:,ik)
   kq = kk + qpt
   ! Find the index of the k+q point
   call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:), [1,1,1])
   bks_mask(:,ik,:) = .True.
   bks_mask_kq(:,ikq,:) = .True.
 end do

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 ! Initialize the wavefunction descriptors
 call wfd_init(wfd_k,cryst,pawtab,psps,keep_ur,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_k%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_k%print(header="Wavefunctions on the k-points grid",mode_paral='PERS')

 ABI_MALLOC(wfd_istwfk, (nkpt_kq))
 wfd_istwfk = 1

 call wfd_init(wfd_kq,cryst,pawtab,psps,keep_ur_kq,mband_kq,nband_kq,nkpt_kq,nsppol,bks_mask_kq,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_kq%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_kq%print(header="Wavefunctions on the q-shifted k-points grid",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband_kq)
 ABI_FREE(bks_mask_kq)
 ABI_FREE(keep_ur_kq)

 ! Read wavefunctions on the k-points grid and q-shifted k-points grid.
 call wfd_k%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))
 if (.False.) call wfd_k%test_ortho(cryst,pawtab,unit=std_out,mode_paral="PERS")

 call wfd_kq%read_wfk(wfq_path, iomode_from_fname(wfq_path))
 if (.False.) call wfd_kq%test_ortho(cryst,pawtab,unit=std_out,mode_paral="PERS")

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
 opt_gvnlx1 = 0 ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,natom,&
   dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
   usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda,&
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)

 ! Allocate vlocal. Note nvloc
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ! Allocate work space arrays.
 ABI_MALLOC(blkflg, (natom3,natom3))
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))

 call cwtime(cpu, wall, gflops, "start")

 interpolated = 0
 if (dtset%eph_use_ftinterp /= 0) then
   MSG_WARNING(sjoin("Enforcing FT interpolation for q-point", ktoa(qpt)))
   comm_rpt = xmpi_comm_self
   wr_path = ""
   method = dtset%userid
   !method = 0
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, method, wr_path, comm_rpt)
   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfftf, nspden, dvdb%my_npert))
   call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)
   interpolated = 1
 else
   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
     call dvdb%readsym_allv1(db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
   else
     MSG_WARNING(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
   end if
 end if

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

 ! Allocate vlocal1 with correct cplex. Note nvloc
 ABI_MALLOC_OR_DIE(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)

 ABI_MALLOC(displ_cart, (2,3*cryst%natom,3*cryst%natom))
 ABI_MALLOC(displ_red, (2,3*cryst%natom,3*cryst%natom))

 if (dtset%eph_task == 2) then
   ! Write GKK files (1 file for perturbation)
   ABI_MALLOC(gkk, (2*mband*nsppol,nkpt,1,1,mband_kq))

 else if (dtset%eph_task == -2) then
   ! Write GKQ file with all perturbations. gkq are given in the atom representation.
   ! TODO: mband_kq == mband
   ABI_MALLOC(gkq_atm, (2, mband_kq, mband, nkpt))
   if (i_am_master) then
     call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)
     fname = strcat(dtfil%filnam_ds(4), "_GKQ.nc")
#ifdef HAVE_NETCDF
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKQ file")
     NCF_CHECK(cryst%ncwrite(ncid))
     ! Write bands on k mesh.
     NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
     ncerr = nctk_def_dims(ncid, [nctkdim_t('number_of_phonon_modes', natom3)], defmode=.True.)
     NCF_CHECK(ncerr)
     ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
       "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"])
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"]))

     ! Define EPH arrays
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t('qpoint', "dp" , 'number_of_reduced_dimensions'), &
       nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'), &
       nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms"), &
       nctkarr_t("eigenvalues_kq", "dp", "max_number_of_states, number_of_kpoints, number_of_spins"), &
       nctkarr_t('phfreqs', "dp", 'number_of_phonon_modes'), &
       nctkarr_t('phdispl_cart', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       !nctkarr_t('phdispl_cart_qvers', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       nctkarr_t('phdispl_red', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       nctkarr_t("gkq_representation", "char", "character_string_length"), &
       nctkarr_t('gkq', "dp", &
         'complex, max_number_of_states, max_number_of_states, number_of_phonon_modes, number_of_kpoints, number_of_spins') &
     ])
     NCF_CHECK(ncerr)
     ! Write data.
     NCF_CHECK(nctk_set_datamode(ncid))
     ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
       "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"], &
       [dtset%symdynmat, dtset%symv1scf, dtset%dvdb_add_lr, interpolated])
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"], [dvdb%qdamp]))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpoint"), qpt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "emacro_cart"), dvdb%dielt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "becs_cart"), dvdb%zeff))
     ABI_MALLOC(eigens_kq, (ebands_kq%mband, nkpt, nsppol))
     do ik=1,nkpt
       kk = ebands_k%kptns(:,ik)
       kq = kk + qpt
       ! Find the index of the k+q point
       call findqg0(ikq, g0_k, kq, nkpt_kq, ebands_kq%kptns, [1,1,1])
       eigens_kq(:, ik, :) = ebands_kq%eig(:, ikq, :)
     end do
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eigenvalues_kq"), eigens_kq))
     ABI_FREE(eigens_kq)
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreqs"), phfrq))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_cart'), displ_cart))
     ! Add phonon displacement for qvers
     !call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red, nanaqdir="reduced")
     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_cart_qvers'), displ_cart))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_red'), displ_red))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gkq_representation"), "atom"))
#endif
   end if
 else
   MSG_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 ! Loop over all 3*natom perturbations.
 do ipc=1,natom3
   idir = mod(ipc-1, 3) + 1
   ipert = (ipc - idir) / 3 + 1
   write(msg, '(a,2i4)') " Treating ipert, idir = ", ipert, idir
   call wrtout(std_out, msg, do_flush=.True.)
   if (dtset%eph_task == 2) gkk = zero

   do spin=1,nsppol
     if (dtset%eph_task == -2) gkq_atm = zero

     ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
     call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
               pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))

     ! Continue to initialize the Hamiltonian
     call load_spin_hamiltonian(gs_hamkq,spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ABI_MALLOC(bras, (2, mpw*nspinor, mband))
     ABI_MALLOC(kets, (2, mpw*nspinor, mband))
     ABI_MALLOC(h1_kets, (2, mpw*nspinor, mband))

     ! GKA: This little block used to be right after the perturbation loop
     ! Prepare application of the NL part.
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=.true.)
     call load_spin_rf_hamiltonian(rf_hamkq,spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

     do ik=1,nkpt
       ! Only do a subset a k-points
       if (.not. (ik >= my_kstart .and. ik <= my_kstop)) cycle

       kk = ebands_k%kptns(:,ik)
       kq = kk + qpt
       ! Find the index of the k+q point
       call findqg0(ikq, g0_k, kq, nkpt_kq, ebands_kq%kptns, [1,1,1])

       ! Copy u_k(G)
       istwf_k = wfd_k%istwfk(ik); npw_k = wfd_k%npwarr(ik)
       ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
       kg_k(:,1:npw_k) = wfd_k%kdata(ik)%kg_k
       do ib2=1,mband
         call wfd_k%copy_cg(ib2, ik, spin, kets(1,1,ib2))
       end do

       ! Copy u_kq(G)
       istwf_kq = wfd_kq%istwfk(ikq); npw_kq = wfd_kq%npwarr(ikq)
       ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
       kg_kq(:,1:npw_kq) = wfd_kq%kdata(ikq)%kg_k
       do ib1=1,mband_kq
         call wfd_kq%copy_cg(ib1, ikq, spin, bras(1,1,ib1))
       end do

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = (psps%usepaw==1)
       sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

       ! GKA: Previous loop on 3*natom perturbations used to start here
       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&    ! In
         cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&            ! In
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&           ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)       ! Out

       ! Calculate dvscf * psi_k, results stored in h1_kets on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,mband
         eig0nk = ebands_k%eig(ib2,ik,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,kets(:,:,ib2),cwaveprj0,h1_kets(:,:,ib2),&
                      grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,mpi_enreg,optlocal,&
                      optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
       end do

       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
       ABI_FREE(gs1c)
       ABI_SFREE(ph3d1)

       ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
       ! The array eig1_k contains:
       !
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
       do ib2=1,mband
         do ib1=1,mband_kq
           call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras(1,1,ib1),h1_kets(1,1,ib2),&
             mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           band = 2*ib2-1 + (spin-1) * 2 * mband

           if (dtset%eph_task == 2) then
             gkk(band,ik,1,1,ib1) = dotr
             gkk(band+1,ik,1,1,ib1) = doti
           else
             gkq_atm(:, ib1, ib2, ik) = [dotr, doti]
           end if

         end do
       end do

     end do ! ikpt

     ABI_FREE(bras)
     ABI_FREE(kets)
     ABI_FREE(h1_kets)
     call destroy_rf_hamiltonian(rf_hamkq)

     if (dtset%eph_task == -2) then
       ! Gather the k-points computed by all processes
       call xmpi_sum_master(gkq_atm, master, comm, ierr)
       if (i_am_master) then
         ! Write the netCDF file.
#ifdef HAVE_NETCDF
         ncerr = nf90_put_var(ncid, nctk_idname(ncid, "gkq"), gkq_atm, &
           start=[1, 1, 1, ipc, 1, 1], count=[2, mband, mband, 1, nkpt, spin])
         NCF_CHECK(ncerr)
#endif
       end if
     end if
   end do ! spin

   if (dtset%eph_task == 2) then
     ! Gather the k-points computed by all processes
     call xmpi_sum_master(gkk,master,comm,ierr)
     ! Init a gkk_t object
     call gkk_init(gkk,gkk2d,mband,nsppol,nkpt,1,1)
     ! Write the netCDF file.
     call appdig(ipc,dtfil%fnameabo_gkk,gkkfilnam)
     fname = strcat(gkkfilnam, ".nc")
#ifdef HAVE_NETCDF
     if (i_am_master) then
       NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
       NCF_CHECK(cryst%ncwrite(ncid))
       NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
       call gkk_ncwrite(gkk2d, qpt, 1.0_dp,  ncid)
       NCF_CHECK(nf90_close(ncid))
     end if
#endif
     ! Free memory
     call gkk_free(gkk2d)
   end if
 end do ! ipc (loop over 3*natom atomic perturbations)

 call cwtime(cpu, wall, gflops, "stop")
 write(msg, '(2a)') " Computation of gkq matrix elements with ", trim(what)
 call wrtout([std_out, ab_out], msg, do_flush=.True.)
 call wrtout(std_out, sjoin("cpu-time:", sec2str(cpu), ",wall-time:", sec2str(wall)), do_flush=.True.)

 if (dtset%eph_task == -2 .and. i_am_master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 ! ===========
 ! Free memory
 ! ===========
 ABI_SFREE(gkk)
 ABI_SFREE(gkq_atm)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_FREE(v1scf)
 ABI_FREE(vlocal1)
 ABI_FREE(gvnlx1)
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
 call wfd_k%free()
 call wfd_kq%free()
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
 integer :: iq, db_iqpt, cplex, nfft, comm, ip, idir, ipert
 logical :: with_lr_model
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
!arrays
 integer :: ngfft(18)
 real(dp) :: phfreqs(3*cryst%natom),qpt(3)
 real(dp) :: displ_cart(2,3*cryst%natom,3*cryst%natom), displ_red(2,3*cryst%natom,3*cryst%natom)
 real(dp),allocatable :: v1scf(:,:,:,:), v1_qnu(:,:,:,:)
 real(dp),allocatable :: v1lr_atm(:,:,:,:), v1lr_qnu(:,:,:,:)

!************************************************************************

 call wrtout(std_out, sjoin("Writing Delta V_{q,nu)(r) potentials to file:", path), do_flush=.True.)

 comm = xmpi_comm_self
 call ngfft_seq(ngfft, dvdb%ngfft3_v1(:,1)); nfft = product(ngfft(1:3))

 call dvdb%open_read(ngfft, xmpi_comm_self)
 call dvdb%print()
 !call dvdb%list_perts([-1, -1, -1])
 with_lr_model = .True.

 ! Create netcdf file.
#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_open_create(ncid, path, comm))
 NCF_CHECK(cryst%ncwrite(ncid))

 ! Add other dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("nfft", nfft), nctkdim_t("nspden", dvdb%nspden), &
   nctkdim_t("natom3", 3 * cryst%natom), nctkdim_t("nqlist", nqlist)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Define arrays
 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("ngfft", "int", "three"), &
   nctkarr_t("qlist", "dp", "three, nqlist"), &
   nctkarr_t("phfreqs", "dp", "natom3, nqlist"), &
   nctkarr_t("displ_cart", "dp", "two, natom3, natom3, nqlist"), &
   nctkarr_t("v1_qnu", "dp", "two, nfft, nspden, natom3, nqlist")])
 NCF_CHECK(ncerr)

 if (with_lr_model) then
   NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t("v1lr_qnu", "dp", "two, nfft, nspden, natom3, nqlist")]))
 end if

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngfft"), ngfft(1:3)))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qlist"), qlist))
#endif

 ABI_MALLOC(v1_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 if (with_lr_model) then
   ABI_MALLOC(v1lr_atm, (2, nfft, dvdb%nspden, dvdb%natom3))
   ABI_MALLOC(v1lr_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 end if

 do iq=1,nqlist
   qpt = qlist(:,iq)
   call ifc_fourq(ifc, cryst, qpt, phfreqs, displ_cart, out_displ_red=displ_red)

   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     if (prtvol > 0) call wrtout(std_out, sjoin("Found:", ktoa(qpt), "in DVDB with index", itoa(db_iqpt)))
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfft, nspden, 3*natom))
     call dvdb%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)
   else
     MSG_ERROR(sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file"))
   end if

   ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
   ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
   ! v1_qnu(2, nfft, nspden, natom3), v1scf(cplex, nfft, nspden, natom3)
   call v1atm_to_vqnu(cplex, nfft, dvdb%nspden, dvdb%natom3, v1scf, displ_red, v1_qnu)

   if (with_lr_model) then
     ! Compute LR model in the atomic representation then compute phonon representation in v1lr_qnu.
     v1lr_atm = zero
     do idir=1,3
       do ipert=1,dvdb%natom
         ip = (ipert - 1) * 3 + idir
         call dvdb%v1r_long_range(qpt, ipert, idir, nfft, ngfft, v1lr_atm(:,:,1,ip))
         if (dvdb%nspden == 2) v1lr_atm(:,:,2,ip) = v1lr_atm(:,:,1,ip)
       end do
     end do
     call v1atm_to_vqnu(2, nfft, dvdb%nspden, dvdb%natom3, v1lr_atm, displ_red, v1lr_qnu)
   end if

#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreqs"), phfreqs, start=[1, iq]))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "displ_cart"), displ_cart, start=[1, 1, 1, iq]))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1_qnu"), v1_qnu, start=[1, 1, 1, 1, iq]))
   if (with_lr_model) then
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1lr_qnu"), v1lr_qnu, start=[1, 1, 1, 1, iq]))
   end if
#endif
   ABI_FREE(v1scf)
 end do

 ABI_FREE(v1_qnu)
 ABI_SFREE(v1lr_atm)
 ABI_SFREE(v1lr_qnu)

#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_close(ncid))
#endif
 call dvdb%close()

 call wrtout(std_out, "dvqnu file written", do_flush=.True.)

end subroutine ncwrite_v1qnu
!!***

!----------------------------------------------------------------------

!!****f* m_gkk/v1atm_to_vqnu
!! NAME
!!  v1atm_to_vqnu
!!
!! FUNCTION
!!  Receive potentials in atomic representation and return potential in phonon representation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine v1atm_to_vqnu(cplex, nfft, nspden, natom3, v1_atm, displ_red, v1_qnu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex, nfft, nspden, natom3
!arrays
 real(dp),intent(in) :: v1_atm(cplex, nfft, nspden, natom3)
 real(dp),intent(out) :: v1_qnu(2, nfft, nspden, natom3)
 real(dp),intent(in) :: displ_red(2, natom3, natom3)

!Local variables-------------------------------
!scalars
 integer :: nu, ip, ispden

!************************************************************************

 do nu=1,natom3
   ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
   ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
   ! v1_qnu(2, nfft, nspden, natom3), v1_atm(cplex, nfft, nspden, natom3)
   v1_qnu(:, :, :, nu) = zero
   do ip=1,natom3
     do ispden=1,nspden
       if (cplex == 2) then
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + &
           displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip) - displ_red(2,ip,nu) * v1_atm(2,:,ispden,ip)
         v1_qnu(2, :, ispden, nu) = v1_qnu(2, :, ispden, nu) + &
           displ_red(2,ip,nu) * v1_atm(1,:,ispden,ip) + displ_red(1,ip,nu) * v1_atm(2,:,ispden,ip)
       else
         ! Gamma point. d(q) = d(-q)* --> d is real.
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip)
       end if
     end do
   end do
 end do

end subroutine v1atm_to_vqnu
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

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

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
