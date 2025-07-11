!!****m* ABINIT/m_gkk
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of electron-phonon coupling matrix elements (gkk)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (GKA, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gkk

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_dtset
 use m_ifc
 use m_ebands
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfk
 use m_nctk
 use m_dtfil
 use netcdf

 use defs_abitypes,    only : MPI_type
 use m_time,           only : cwtime, sec2str
 use m_io_tools,       only : iomode_from_fname
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_symtk,          only : littlegroup_q
 use m_fftcore,        only : get_kg
 use defs_datatypes,   only : pseudopotential_type
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
 use m_ephtk,          only : ephtk_v1atm_to_vqnu

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
 integer,parameter :: tim_getgh1c=1, berryopt0=0, useylmgr1=0, master=0, qptopt1 = 1
 integer :: my_rank,nproc,mband,mband_kq,my_minb,my_maxb,nsppol,nkpt,nkpt_kq,idir,ipert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor
 integer :: ib1,ib2,band,ik,ikq,timerev_q
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq, comm_rpt
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
 character(len=fnlen) :: fname, gkkfilnam
!arrays
 integer :: g0_k(3),symq(4,2,cryst%nsym), units(2)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:),nband_kq(:,:),wfd_istwfk(:)
 real(dp) :: ylmgr_kq_dum(1,1,1) ! ylmgr_k_dum(1,1,1),
 real(dp) :: kk(3),kq(3),qpt(3),phfrq(3*cryst%natom),dvdb_qdamp(1)
 real(dp),allocatable :: displ_cart(:,:,:),displ_red(:,:,:), eigens_kq(:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw_kq(:),kpg_kq(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnl_k(:,:,:,:),ffnl_kq(:,:,:,:),ph3d_k(:,:,:),ph3d_kq(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),gkk(:,:,:,:,:), bras(:,:,:),kets(:,:,:),h1_kets(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:), ylm_kq(:,:),ylm_k(:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnlx1(:,:), gs1c(:,:), gkq_atm(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),bks_mask_kq(:,:,:),keep_ur(:,:,:),keep_ur_kq(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)
!************************************************************************

 units = [std_out, ab_out]

 what = "(GKK files)"; if (dtset%eph_task == -2) what = "GKQ file"
 write(msg, '(3a)') " Computation of electron-phonon coupling matrix elements ", trim(what), ch10
 call wrtout(units, msg, do_flush=.True.)

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm); i_am_master = my_rank == master

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands_k%nsppol; nspinor = ebands_k%nspinor; nspden = dtset%nspden
 nkpt = ebands_k%nkpt; mband = ebands_k%mband; nkpt_kq = ebands_kq%nkpt; mband_kq = ebands_kq%mband
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

 call wfd_k%print([std_out], header="Wavefunctions on the k-points grid")

 ABI_MALLOC(wfd_istwfk, (nkpt_kq))
 wfd_istwfk = 1

 call wfd_init(wfd_kq,cryst,pawtab,psps,keep_ur_kq,mband_kq,nband_kq,nkpt_kq,nsppol,bks_mask_kq,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_kq%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_kq%print([std_out], header="Wavefunctions on the q-shifted k-points grid")

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband_kq)
 ABI_FREE(bks_mask_kq)
 ABI_FREE(keep_ur_kq)

 ! Read wavefunctions on the k-points grid and q-shifted k-points grid.
 call wfd_k%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))
 call wfd_kq%read_wfk(wfq_path, iomode_from_fname(wfq_path))

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
 ! TODO: These arrays should be allocated with npw_k and npw_kq
 ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))

 ! TODO FOR PAW
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

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

 call gs_hamkq%init(psps,pawtab,nspinor,nsppol,nspden,natom,&
   dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
   usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,gpu_option=dtset%gpu_option,&
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)

 ! Allocate vlocal. Note nvloc
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ! Allocate work space arrays.
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))

 call cwtime(cpu, wall, gflops, "start")

 interpolated = 0
 if (dtset%eph_use_ftinterp /= 0) then
   ABI_WARNING(sjoin("Enforcing FT interpolation for q-point", ktoa(qpt)))
   comm_rpt = xmpi_comm_self
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt1, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
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
     ABI_WARNING(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
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
   ! TODO: Assuming mband_kq == mband
   ABI_MALLOC(gkq_atm, (2, mband_kq, mband, nkpt))
   if (i_am_master) then
     call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)
     fname = strcat(dtfil%filnam_ds(4), "_GKQ.nc")
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKQ file")
     NCF_CHECK(cryst%ncwrite(ncid))
     ! Write bands on k mesh.
     NCF_CHECK(ebands_k%ncwrite(ncid))
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
     dvdb_qdamp = dvdb%qdamp
     NCF_CHECK(nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"], dvdb_qdamp))
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
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_red'), displ_red))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gkq_representation"), "atom"))
   end if ! master

 else
   ABI_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 ! Loop over all 3*natom perturbations.
 do ipc=1,natom3
   idir = mod(ipc-1, 3) + 1
   ipert = (ipc - idir) / 3 + 1
   write(msg, '(a,2(i0,1x))') " Treating ipert, idir = ", ipert, idir
   call wrtout(std_out, msg, do_flush=.True.)
   if (dtset%eph_task == 2) gkk = zero

   do spin=1,nsppol
     if (dtset%eph_task == -2) gkq_atm = zero

     ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
     call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
               pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))

     ! Continue to initialize the Hamiltonian
     call gs_hamkq%load_spin(spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ABI_MALLOC(bras, (2, mpw*nspinor, mband))
     ABI_MALLOC(kets, (2, mpw*nspinor, mband))
     ABI_MALLOC(h1_kets, (2, mpw*nspinor, mband))

     ! GKA: This little block used to be right after the perturbation loop
     ! Prepare application of the NL part.
     call rf_hamkq%init(cplex,gs_hamkq,ipert,has_e1kbsc=.true.)
     call rf_hamkq%load_spin(spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

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
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq_dum,&       ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg_kq,kinpw_kq,ffnl_k,ffnl_kq,ph3d_k,ph3d_kq)       ! Out

       ! Calculate dvscf * psi_k, results stored in h1_kets on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,mband
         eig0nk = ebands_k%eig(ib2,ik,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,kets(:,:,ib2),cwaveprj0,h1_kets(:,:,ib2),&
                      grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,(/eshift/),mpi_enreg,1,optlocal,&
                      optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
       end do

       ABI_FREE(kinpw_kq)
       ABI_FREE(kpg_k)
       ABI_FREE(kpg_kq)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnl_k)
       ABI_FREE(ffnl_kq)
       ABI_FREE(gs1c)
       ABI_FREE(ph3d_k)
       ABI_SFREE(ph3d_kq)

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
         end do ! ib1
       end do ! ib2

     end do ! ikpt

     ABI_FREE(bras)
     ABI_FREE(kets)
     ABI_FREE(h1_kets)
     call rf_hamkq%free()

     if (dtset%eph_task == -2) then
       ! Gather the k-points computed by all processes
       call xmpi_sum_master(gkq_atm, master, comm, ierr)
       if (i_am_master) then
         ! Write the netCDF file.
         ncerr = nf90_put_var(ncid, nctk_idname(ncid, "gkq"), gkq_atm, &
                              start=[1, 1, 1, ipc, 1, spin], count=[2, mband, mband, 1, nkpt, 1])
         NCF_CHECK(ncerr)
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
     if (i_am_master) then
       NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
       NCF_CHECK(cryst%ncwrite(ncid))
       NCF_CHECK(ebands_k%ncwrite(ncid))
       call gkk_ncwrite(gkk2d, qpt, 1.0_dp,  ncid)
       NCF_CHECK(nf90_close(ncid))
     end if
     ! Free memory
     call gkk_free(gkk2d)
   end if
 end do ! ipc (loop over 3*natom atomic perturbations)

 call cwtime(cpu, wall, gflops, "stop")
 write(msg, '(2a)') " Computation of gkq matrix elements with ", trim(what)
 call wrtout(units, msg, do_flush=.True.)
 call wrtout(std_out, sjoin("cpu-time:", sec2str(cpu), ",wall-time:", sec2str(wall)), do_flush=.True.)

 if (dtset%eph_task == -2 .and. i_am_master) then
   NCF_CHECK(nf90_close(ncid))
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

 call gs_hamkq%free()
 call wfd_k%free()
 call wfd_kq%free()
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)

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
!!  dtset<dataset_type>= Input variables.
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  out_ncpath=Name of the netcdf file.
!!
!! OUTPUT
!!  Only writing
!!
!! SOURCE

subroutine ncwrite_v1qnu(dvdb, dtset, ifc, out_ncpath)

 use m_bz_mesh, only : kpath_t, kpath_new

!Arguments ------------------------------------
 class(dvdb_t),intent(inout) :: dvdb
 type(dataset_type),target,intent(in) :: dtset
 type(ifc_type),intent(in) :: ifc
 character(len=*),intent(in) :: out_ncpath

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, qptopt1 = 1
 integer :: db_iqpt, cplex, nfft, comm, ip, idir, ipert, my_rank, interpolated, comm_rpt, ncid, ncerr
 integer :: iq, nu, iatom, ii, jj, kk
 real(dp) :: inv_qepsq, qtau, phre, phim, rtmp
 logical :: with_lr_model
 type(kpath_t) :: qpath
!arrays
 integer :: ngfft(18), units(2)
 real(dp) :: phfreqs(dvdb%natom3),qpt(3)
 real(dp) :: displ_cart(2,3, dvdb%cryst%natom, dvdb%natom3), displ_red(2,dvdb%natom3,dvdb%natom3)
 real(dp),allocatable :: v1scf(:,:,:,:), v1_qnu(:,:,:,:), v1lr_atm(:,:,:,:), v1lr_qnu(:,:,:,:)
 real(dp) :: bounds(3,6), qpt_red(3), qpt_cart(3), glr(3), values(dvdb%natom3)

!************************************************************************

 ! +0.50000  +0.50000  +0.50000  # L
 ! +0.00000  +0.00000  +0.00000  # $\Gamma$
 ! +0.50000  +0.00000  +0.50000  # X
 ! +0.50000  +0.25000  +0.75000  # W
 ! +0.37500  +0.37500  +0.75000  # K
 ! +0.00000  +0.00000  +0.00000  # $\Gamma$
 ! +0.37500  +0.37500  +0.75000  # K

 ! +0.62500  +0.25000  +0.62500  # U
 ! +0.50000  +0.50000  +0.50000  # L
 ! +0.37500  +0.37500  +0.75000  # K
 ! +0.62500  +0.25000  +0.62500  # U
 ! +0.50000  +0.00000  +0.50000  # X

 bounds(:, 1) = tol3 * [+0.50000,  +0.50000, +0.50000] !  # L
 bounds(:, 2) = tol3 * [+0.00000,  +0.00000, +0.00000] !  # $\Gamma$
 bounds(:, 3) = tol3 * [+0.50000,  +0.00000, +0.50000] !  # X
 bounds(:, 4) = tol3 * [+0.37500,  +0.37500, +0.75000] !  # K
 bounds(:, 5) = tol3 * [+0.00000,  +0.00000, +0.00000] !  # $\Gamma$
 bounds(:, 6) = tol3 * [+0.50000,  +0.25000, +0.75000] !  # W

 qpath = kpath_new(bounds, dvdb%cryst%gprimd, dtset%ndivsm)

 units = [std_out, ab_out]

 do iq=1,qpath%npts
   qpt_red = qpath%points(:, iq)
   qpt_cart = two_pi * matmul(dvdb%cryst%gprimd, qpt_red)
   inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
   call ifc%fourq(dvdb%cryst, qpt_red, phfreqs, displ_cart)
   do nu=1, dvdb%natom3
     glr = zero
     do iatom=1, dvdb%cryst%natom
       ! Phase factor exp(-i (q+G) . tau)
       qtau = - two_pi * dot_product(qpt_red, dvdb%cryst%xred(:,iatom))
       phre = cos(qtau); phim = sin(qtau)
       do jj=1,3
         do ii=1,3
           do kk=1,3
             rtmp = dvdb%qstar(ii, jj, kk, iatom) * qpt_cart(ii) * qpt_cart(jj)
             glr(1) = glr(1) + rtmp * (displ_cart(1, kk, iatom, nu) * phre - displ_cart(2, kk, iatom, nu) * phim)
             glr(2) = glr(2) + rtmp * (displ_cart(2, kk, iatom, nu) * phre + displ_cart(1, kk, iatom, nu) * phre)
           end do
         end do
       end do
     end do
     glr = half * (glr / inv_qepsq) * (four_pi / dvdb%cryst%ucvol)
     values(nu) = (glr(1) ** 2 + glr(2) ** 2) / (two *  phfreqs(nu))
   end do ! nu
   write(std_out, "(i0, 4(f9.6), /, (es18.6, 1x))") iq, qpt_red, phfreqs(nu), (values(nu), nu=1, 3*dvdb%natom)
 end do ! iqpt

 call qpath%free()
 return

 my_rank = xmpi_comm_rank(dvdb%comm)
 comm = dvdb%comm
 qpt = dtset%qptn

 call wrtout(std_out, sjoin(" Writing Delta V_{q,nu)(r) potentials to file:", out_ncpath), do_flush=.True.)
 call wrtout(units, sjoin(ch10, "- Results stored in: ", out_ncpath))
 call wrtout(std_out, sjoin(" Using qpt:", ktoa(qpt)))
 !call wrtout(units, " Use `abiopen.py out_V1QAVG.nc -e` to visualize results")
 call dvdb%print(unit=std_out)

 ! Define FFT mesh
 ngfft = dvdb%ngfft
 nfft = product(ngfft(1:3))

 if (dtset%eph_task == -16) then
   call wrtout(units, " Assuming q-point already in the DVDB file. No interpolation.")
   interpolated = 0

 else if (dtset%eph_task == +16) then
   call wrtout(units, " Using Fourier interpolation.")
    comm_rpt = xmpi_comm_self
    call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt1, 1, dtset%ddb_shiftq, nfft, ngfft, comm_rpt)
    interpolated = 1
 else
   ABI_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 with_lr_model = .True.

 ! Create netcdf file.
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, out_ncpath, comm))
   NCF_CHECK(dvdb%cryst%ncwrite(ncid))

   ! Add other dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nfft", nfft), nctkdim_t("nspden", dvdb%nspden), &
     nctkdim_t("natom3", 3 * dvdb%cryst%natom)], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngfft", "int", "three"), &
     nctkarr_t("qpt", "dp", "three"), &
     nctkarr_t("phfreqs", "dp", "natom3"), &
     nctkarr_t("displ_cart", "dp", "two, natom3, natom3"), &
     nctkarr_t("v1_qnu", "dp", "two, nfft, nspden, natom3")])
   NCF_CHECK(ncerr)

   if (with_lr_model) then
     NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t("v1lr_qnu", "dp", "two, nfft, nspden, natom3")]))
   end if

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngfft"), ngfft(1:3)))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpt"), qpt))
 end if

 ABI_MALLOC(v1_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 if (with_lr_model) then
   ABI_MALLOC(v1lr_atm, (2, nfft, dvdb%nspden, dvdb%natom3))
   ABI_MALLOC(v1lr_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 end if

 ! Get phonon freqs and displacemented for this q-point.
 call ifc%fourq(dvdb%cryst, qpt, phfreqs, displ_cart, out_displ_red=displ_red)

 if (interpolated == 0) then
   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfft, nspden, 3*natom))
     call dvdb%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)
   else
     ABI_ERROR(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
   end if
 else

   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfft, dvdb%nspden, dvdb%my_npert))
   call dvdb%ftinterp_qpt(qpt, nfft, ngfft, v1scf, dvdb%comm_rpt)
 end if

 ! Compute scattering potential the in phonon representations instead ot atomic one.
 ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
 ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
 ! v1_qnu(2, nfft, nspden, natom3), v1scf(cplex, nfft, nspden, natom3)
 call ephtk_v1atm_to_vqnu(cplex, nfft, dvdb%nspden, dvdb%natom3, v1scf, displ_red, v1_qnu)

 if (with_lr_model) then
   ! Compute LR model in the atomic representation then compute phonon representation in v1lr_qnu.
   v1lr_atm = zero
   do idir=1,3
     do ipert=1,dvdb%natom
       ip = (ipert - 1) * 3 + idir
       call dvdb%get_v1r_long_range(qpt, idir, ipert, nfft, ngfft, v1lr_atm(:,:,1,ip))
       if (dvdb%nspden == 2) v1lr_atm(:,:,2,ip) = v1lr_atm(:,:,1,ip)
     end do
   end do
   call ephtk_v1atm_to_vqnu(2, nfft, dvdb%nspden, dvdb%natom3, v1lr_atm, displ_red, v1lr_qnu)
 end if

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreqs"), phfreqs))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "displ_cart"), displ_cart))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1_qnu"), v1_qnu))
 if (with_lr_model) then
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1lr_qnu"), v1lr_qnu))
 end if

 ABI_FREE(v1scf)
 ABI_FREE(v1_qnu)
 ABI_SFREE(v1lr_atm)
 ABI_SFREE(v1lr_qnu)

 NCF_CHECK(nf90_close(ncid))
 call dvdb%close()

 call wrtout(std_out, "dvqnu file written", do_flush=.True.)

end subroutine ncwrite_v1qnu
!!***

!----------------------------------------------------------------------

!!****f* m_gkk/find_mpw
!! NAME
!!  find_mpw
!!
!! FUNCTION
!!  Look at all k-points and spins to find the maximum number of plane waves.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine find_mpw(mpw, kpts, nsppol, nkpt, gmet, ecut, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: mpw
 integer,intent(in) :: nsppol, nkpt, comm
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: kpts(3,nkpt), gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: my_rank, cnt, nproc, ierr, ispin, ikpt, my_mpw, onpw
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kpt(3)

!************************************************************************

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
