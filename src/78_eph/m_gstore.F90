!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_gstore

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_krank
 use m_htetra
 use m_tetrahedron
 use libtetrabz
 use m_ifc
 use m_ebands
 use m_fstab
 use iso_c_binding
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_crystal
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_dtset
 use m_dtfil
 use m_wfd
 use m_ephtk
 use m_mkffnl

 use defs_abitypes,    only : mpi_type
 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : toupper, itoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, wrap2_pmhalf, simpson_int, simpson, mkherm, get_diag, isdiagmat
 use m_io_tools,       only : open_file, iomode_from_fname
 use m_symtk,          only : littlegroup_q
 use m_geometry,       only : normv
 use m_special_funcs,  only : gaussian
 use m_fftcore,        only : ngfft_seq, get_kg
 use m_cgtools,        only : cg_zdotc
 use m_kg,             only : getph, mkkpg
 use m_dynmat,         only : symdyma, ftgam_init, ftgam, asrif9
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_bz_mesh,        only : kpath_t, kpath_new
 use m_special_funcs,  only : fermi_dirac
 use m_kpts,           only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt, listkk, kpts_timrev_from_kptopt
 use defs_elphon,      only : complete_gamma !, complete_gamma_tr
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type

 !use defs_datatypes,    only : ebands_t
 !use m_time,            only : cwtime, cwtime_report
! use m_symtk,           only : matr3inv
! use m_numeric_tools,   only : arth, inrange, wrap2_pmhalf
! use m_special_funcs,   only : gaussian
 !use m_fstrings,        only : strcat, ltoa, itoa, ftoa, ktoa, sjoin
 !use m_kpts,            only : kpts_timrev_from_kptopt, listkk, kpts_ibz_from_kptrlatt
! use m_occ,             only : occ_fd, occ_be

  use m_phgamma, only : phgamma_t

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gqk_t
!! NAME
!! gqk_t
!!
!! FUNCTION
!!  Stores e-p matrix elements for a given spin as several dimensions depende on
!!  isppol when nsppol == 2.
!!
!! NOTES
!!  global and local mean that the ...
!!
!!
!! SOURCE

type, public :: gqk_t

  integer :: cplex
  ! 1 if |g|^2,
  ! 2 if complex g (global)

  integer :: nb
  ! Number of bands included in the calculation for this spin
  ! Global as this dimension is not MPI-distributed due to (m,n) pairs.

  integer :: bstart = 1
  ! The fist band starts at bstart. (global index)
  ! Used to select bands inside energy window.

  integer :: glob_nk, glob_nq
  ! Total number of k/q points in global matrix.
  ! Use kzone_type and qzone_type to interpret these values (IBZ or BZ)
  ! Note that k-points can be filtered by erange.

  integer :: my_nk, my_nq
  ! Number of k/q points treated by this MPI proc.
  ! Used to loop and allocate local arrays.

  integer :: my_kstart, my_qstart
  ! Index of the first k/q point in the global matrix treated by this MPI proc

  integer :: nkibz, nqibz
  ! Number of k/q points in the IBZ.

  real(dp),allocatable :: kibz(:,:)
  ! k-points in the IBZ
  ! (3, nkibz)

  real(dp),allocatable :: qibz(:,:)
  ! q-points in the IBZ
  ! (3, nqibz)

  !integer,allocatable :: my_k2ibz(:,:)
  ! (6, my_nk)
  ! Mapping my_kpoints --> kibz

  !integer,allocatable :: my_kpq2ibz_k(:,:,:)
  ! (6, my_nk, my_nq)
  ! Mapping k + q --> kibz for each (k, q) treated by this MPI proc.

  !integer,allocatable :: my_q2ibz(:,:)
  ! (6, my_nq)
  ! Mapping my_qpoints --> qibz

  !complex(dp),allocatable :: my_vk(:,:,:,:)
  ! group velocities v_{nm,k} for the k-points treated by this MPI proc.
  ! (3, nb, nb, my_nk)

  !real(dp) :: erange(2) = zero

  character(len=10) :: kzone_type, qzone_type
    ! Either "ibz" or "bz"

  real(dp), allocatable :: vals(:,:,:,:,:,:,:)
    ! (cplex, nb, my_nq, nb, my_nk, my_npert)

  type(xcomm_t) :: k_comm
    ! MPI communicator for parallelism over k-points

  type(xcomm_t) :: q_comm
   ! MPI communicator for parallelism over q-points

  type(xcomm_t) :: pert_comm
   ! MPI communicator for parallelism over atomic perturbations.

  !type(xcomm_t) :: gkq_comm

  !integer :: coords_pqbks(3)
  ! ! Cartesian coordinates of this processor in the Cartesian grid.

  !integer :: natom3
  ! 3 * natom

  !integer :: kptopt, qptopt
  ! Option for k-point generation.

  !integer :: timrev
  ! 1 if the use of time-reversal is allowed; 0 otherwise

  !real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  !real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states
  ! (not necessarly equal to global nband, see also bstart and bcount)

  !type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure (does not own memory)

  !type(htetra_t) :: tetra_k
  ! Used to evaluate delta(w - e_{k+q} +/- phw_q) with tetrahedron method.

 contains

     !procedure :: setup_kpoint => gstore_setup_kpoint
     ! Prepare tetrahedron method for given external k-point.

     procedure :: free => gqk_free
     ! Free memory

 end type gqk_t
!!***

!----------------------------------------------------------------------

!!****t* m_gstore/gstore_t
!! NAME
!! gstore_t
!!
!! FUNCTION
!!
!! NOTES
!!
!! SOURCE

type, public :: gstore_t

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: my_nspins
   ! Number of spins treated by this MPI rank

  integer,allocatable :: my_spins(:)
   ! my_spins(my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2 and nspinor == 1

  !kmesh, qmesh
  integer :: kptrlatt(3, 3), qptrlatt(3, 3)
   ! k-mesh and q-mesh

  real(dp) :: kshift(3, 1), qshift(3, 1)
  ! k/q-mesh shift (well, q-mesh is usually gamma-centered)

  real(dp) :: klatt(3, 3), qlatt(3, 3)
  ! Reciprocal of lattice vectors for full kpoint grid.
  ! Used by tetra routines.

  type(xcomm_t) :: spin_comm
    ! MPI communicator over spin (nsppol = 2)

  type(gqk_t), allocatable :: gqk(:)
  ! (nsppol)

  contains

    procedure :: free => gstore_free
    ! Free memory

end type gstore_t
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_build
!! NAME
!! gstore_build
!!
!! FUNCTION
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

type(gstore_t) function gstore_build(cplex, kzone_type, qzone_type, cryst, ebands, ifc, comm) result (gstore)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex, comm
 character(len=*),intent(in) :: kzone_type, qzone_type
 class(ebands_t),intent(in) :: ebands
 class(crystal_t),target,intent(in) :: cryst
 class(ifc_type),intent(in) :: ifc
!arrays

!Local variables-------------------------------
!scalars
 integer :: nprocs, my_rank, ik, ierr, out_nkibz, spin
 !integer :: kptopt, nshiftk, nkibz
 real(dp) :: cpu, wall, gflops
!arrays
 !real(dp) :: rlatt(3,3)
 !integer :: out_kptrlatt(3,3)
 !real(dp),allocatable :: out_kibz(:,:), out_wtk(:)

!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 call cwtime(cpu, wall, gflops, "start")

 ! Set basic parameters.
 gstore%nsppol = ebands%nsppol
 !if (gstore%nsppol == 1) then
 !else
 !end if

 !new%natom3 = ifc%natom * 3
 !new%bstart = bstart
 !new%kptopt = kptopt
 !new%timrev = kpts_timrev_from_kptopt(new%kptopt)
 !new%nibz = nkibz
 !new%cryst => cryst
 !call alloc_copy(kibz, new%ibz)

 !! Get full BZ (new%nbz, new%bz) and new kptrlatt for tetra.
 !call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, out_nkibz, out_kibz, out_wtk, new%nbz, new%bz, &
 !                            new_kptrlatt=out_kptrlatt)

 !new%kptrlatt = out_kptrlatt
 !rlatt = out_kptrlatt; call matr3inv(rlatt, new%klatt)

 !ABI_CHECK(size(out_kibz, dim=2) == new%nibz, "mismatch in nkibz!")
 !ABI_FREE(out_kibz)
 !ABI_FREE(out_wtk)

 !! Copy eigenvalues in IBZ. Change shape for better performance in other routines.
 !ABI_MALLOC(new%eigkbs_ibz, (new%nibz, new%nbcount, new%nsppol))
 !do ik=1,new%nibz
 !  new%eigkbs_ibz(ik, :, :) = eig_ibz(:, ik, :)
 !end do

 !! Fourier interpolate phonon frequencies on the same mesh.
 !ABI_CALLOC(new%phfrq_ibz, (new%nibz, new%natom3))

 !do ik=1,new%nibz
 !  if (mod(ik, nprocs) /= my_rank) cycle ! mpi-parallelism
 !  call ifc%fourq(cryst, new%ibz(:, ik), phfrq, displ_cart)
 !  new%phfrq_ibz(ik, :) = phfrq
 !end do

 !! Collect results on each rank
 !call xmpi_sum(new%phfrq_ibz, comm, ierr)

 !!new%max_phfrq = maxval(phfrq_ibz)

 call cwtime_report(" gstore_build: ifc_fourq", cpu, wall, gflops)

end function gstore_build
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_free
!! NAME
!! gstore_free
!!
!! FUNCTION
!!  Free dynamic memory.
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

subroutine gstore_free(gstore)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore

!Local variables-------------------------------
 integer :: is

!----------------------------------------------------------------------

 do is=1,size(gstore%gqk)
  call gstore%gqk(is)%free()
 end do

 ABI_SFREE(gstore%my_spins)
 call gstore%spin_comm%free()

end subroutine gstore_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gqk_free
!! NAME
!! gqk_free
!!
!! FUNCTION
!!  Free dynamic memory.
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

subroutine gqk_free(gqk)

!Arguments ------------------------------------
 class(gqk_t),intent(inout) :: gqk

!----------------------------------------------------------------------

 ABI_SFREE(gqk%kibz)
 ABI_SFREE(gqk%qibz)
 ABI_SFREE(gqk%vals)

 ! Communicators
 call gqk%k_comm%free()
 call gqk%q_comm%free()
 call gqk%pert_comm%free()

end subroutine gqk_free
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_compute
!! NAME
!!  gstore_compute
!!
!! FUNCTION
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
!!      m_eph_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine gstore_compute(gstore, wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, &
                          pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(inout) :: gstore
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb
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
 integer,parameter :: tim_getgh1c = 1, berryopt0 = 0, ider0 = 0, idir0 = 0
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1
 integer :: my_rank,nproc,mband,nsppol,nkibz,idir,ipert,iq_ibz, timrev
 integer :: cplex,db_iqpt,natom,natom3,ipc,ipc1,ipc2,nspinor,onpw
 integer :: bstart_k,bstart_kq,nband_k,nband_kq,band_k, band_kq, ib_k, ib_kq !ib1,ib2,
 integer :: ik_ibz,ik_bz,ikq_bz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq,timerev_q
 integer :: ik_fs, myik, mys !, myiq
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq
 integer :: ii,jj,ipw,mpw,my_mpw,mnb,ierr,cnt
 integer :: n1,n2,n3,n4,n5,n6,nspden
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,kq_count,nkpg,nkpg1
 integer :: jene, iene, comm_rpt, nesting, my_npert, imyp, imyq
 real(dp) :: cpu, wall, gflops, cpu_q, wall_q, gflops_q, cpu_k, wall_k, gflops_k, cpu_all, wall_all, gflops_all
 real(dp) :: sigma, ecut, eshift, eig0nk, dksqmax
 logical :: gen_eigenpb, need_velocities, isirr_k, isirr_kq
 type(wfd_t) :: wfd
 type(fstab_t),pointer :: fs
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(phgamma_t) :: gams
 type(ddkop_t) :: ddkop
 type(xcomm_t) :: pert_comm, qs_comm, qpt_comm, bsum_comm, kpt_comm, spin_comm, pkb_comm
 type(krank_t) :: krank
 character(len=500) :: msg
!arrays
 integer :: g0_k(3),g0bz_kq(3),g0_kq(3),symq(4,2,cryst%nsym)
 integer :: work_ngfft(18),gmax(3),my_gmax(3),gamma_ngqpt(3) !g0ibz_kq(3),
 integer :: indkk_kq(6,1)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),gtmp(:,:),nband(:,:),wfd_istwfk(:)
 integer,allocatable :: my_pinfo(:,:), pert_table(:,:) !, qibz_done(:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3), lf(2),rg(2),res(2), vk(3), vkq(3)
 real(dp) :: phfrq(3*cryst%natom)
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: displ_cart(:,:,:,:), displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:), kinpw1(:), kpg1_k(:,:), kpg_k(:,:), dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:), ph3d(:,:,:), ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:), gkk_atm(:,:,:,:) !,gkq_nu(:,:,:,:)
 real(dp),allocatable :: bras_kq(:,:,:), kets_k(:,:,:), h1kets_kq(:,:,:), cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:), ylm_k(:,:), ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:), gvnlx1(:,:), work(:,:,:,:)
 real(dp),allocatable :: gs1c(:,:), v1_work(:,:,:,:), vcar_ibz(:,:,:,:)
 real(dp),allocatable :: wt_ek(:,:), wt_ekq(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(fstab_t),target,allocatable :: fstab(:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)
#ifdef HAVE_MPI
 integer :: ndims, comm_cart, me_cart
 logical :: reorder
 integer :: coords(5)
 integer,allocatable :: dims(:)
 logical,allocatable :: periods(:), keepdim(:)
#endif

!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband
 timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! FFT meshes
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Find Fermi surface k-points
 ! TODO: support kptopt, change setup of k-points if tetra: fist tetra weights then k-points on the Fermi surface!
 ABI_MALLOC(fstab, (nsppol))
 call fstab_init(fstab, ebands, cryst, dtset, comm)

 ! Define q-mesh. eph_ngqpt_fine activates the Fourier interpolation of the DFPT potentials.
 gamma_ngqpt = ifc%ngqpt; if (all(dtset%eph_ngqpt_fine /= 0)) gamma_ngqpt = dtset%eph_ngqpt_fine

 call wrtout(std_out, sjoin("q-mesh for the phonon linewidths:", ltoa(gamma_ngqpt)))
 call wrtout(std_out, sjoin("Will compute", itoa(gams%nqibz), "q-points in the IBZ"))

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Init for sequential execution.
 my_npert = natom3

 if (any(dtset%eph_np_pqbks /= 0)) then
   ! Use parameters from input file.
   pert_comm%nproc = dtset%eph_np_pqbks(1)
   qpt_comm%nproc  = dtset%eph_np_pqbks(2)
   bsum_comm%nproc = dtset%eph_np_pqbks(3)
   kpt_comm%nproc = dtset%eph_np_pqbks(4)
   spin_comm%nproc = dtset%eph_np_pqbks(5)
   my_npert = natom3 / pert_comm%nproc
   ABI_CHECK(my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
   ABI_CHECK(mod(natom3, pert_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")
 else
   ! Automatic grid generation
   ! Not the most efficient distribution if large number of MPI procs.
   ! TODO: Add parallelism over perturbations although it's less efficient than the parallelism over k
   ! It starts to be interesting if we implement symmetries in the k-integration though.

   ! TODO: Spin
   ! Automatic grid generation over q-points and spins.
   !if (new%nsppol == 2 .and. mod(nprocs, 2) == 0) then
   !  spin_comm%nproc = 2
   !  new%qpt_comm%nproc = nprocs / 2
   !else
   !  new%qpt_comm%nproc = nprocs
   !end if

   ! Handle parallelism over perturbations first.
   ! Use MPI communicator to distribute the 3 * natom perturbations to reduce memory requirements for DFPT potentials.
   ! Ideally, perturbations are equally distributed --> total number of CPUs should be divisible by 3 * natom.
   ! or at least, divisible by one integer i for i in [2, 3 * natom - 1].

   if (nsppol == 2 .and. mod(nproc, 2) == 0) then
     spin_comm%nproc = 2
     kpt_comm%nproc = nproc / 2
   else

     ! Try to have 3 perts per proc first because the q-point parallelism is more efficient.
     ! The memory for W(R,r,ipert) will increase though.
     !do cnt=natom,2,-1
     !  if (mod(nprocs, cnt) == 0 .and. mod(natom3, cnt) == 0) then
     !    pert_comm%nproc = cnt; new%my_npert = natom3 / cnt; exit
     !  end if
     !end do

     !if (pert_comm%nproc == 1) then
     !  ! Try again with more procs.
     !  do cnt=natom3,2,-1
     !    if (mod(nprocs, cnt) == 0 .and. mod(natom3, cnt) == 0) then
     !      pert_comm%nproc = cnt; new%my_npert = natom3 / cnt; exit
     !    end if
     !  end do
     !end if

     !if (new%my_npert == natom3 .and. nprocs > 1) then
     !  ABI_WARNING("The number of MPI procs should be divisible by 3*natom to reduce memory requirements!")
     !end if

     kpt_comm%nproc = nproc
   end if
 end if

 ! Consistency check.
 if (pert_comm%nproc * qpt_comm%nproc * bsum_comm%nproc * kpt_comm%nproc * spin_comm%nproc /= nproc) then
   write(msg, "(a,i0,3a, 6(a,1x,i0))") &
     "Cannot create 5d Cartesian grid with total nproc: ", nproc, ch10, &
     "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
     "pert_nproc (", pert_comm%nproc, ") x qpt_nproc (", qpt_comm%nproc, ") x bsum_nproc (", bsum_comm%nproc, &
     ") x kcalc_nproc (", kpt_comm%nproc, ") x spin_nproc (", spin_comm%nproc, ") != ", &
     pert_comm%nproc * qpt_comm%nproc * bsum_comm%nproc * kpt_comm%nproc * spin_comm%nproc
   ABI_ERROR(msg)
 end if

 ABI_CHECK(bsum_comm%nproc == 1, "Band parallelism not implemented in m_phgamma")

#ifdef HAVE_MPI
 ! Create 5d cartesian communicator: 3*natom perturbations, q-points in IBZ, bands in FS, kpoints in FS, spins
 ndims = 5
 ABI_MALLOC(dims, (ndims))
 ABI_MALLOC(periods, (ndims))
 ABI_MALLOC(keepdim, (ndims))
 periods(:) = .False.; reorder = .False.
 dims = [pert_comm%nproc, qpt_comm%nproc, bsum_comm%nproc, kpt_comm%nproc, spin_comm%nproc]

 call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm_cart, ierr)
 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, coords, ierr)

 ! Create communicator to distribute natom3 perturbations.
 keepdim = .False.; keepdim(1) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, pert_comm%value, ierr); pert_comm%me = xmpi_comm_rank(pert_comm%value)

 ! Create communicator for qpoints in self-energy integration.
 keepdim = .False.; keepdim(2) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, qpt_comm%value, ierr); qpt_comm%me = xmpi_comm_rank(qpt_comm%value)

 ! Create communicator for bands for band summation
 keepdim = .False.; keepdim(3) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, bsum_comm%value, ierr); bsum_comm%me = xmpi_comm_rank(bsum_comm%value)

 ! Create communicator for kpoints.
 keepdim = .False.; keepdim(4) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, kpt_comm%value, ierr); kpt_comm%me = xmpi_comm_rank(kpt_comm%value)

 ! Create communicator for spins.
 keepdim = .False.; keepdim(5) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, spin_comm%value, ierr); spin_comm%me = xmpi_comm_rank(spin_comm%value)

 ! Create communicator for the (qpoint, spin) loops
 keepdim = .False.; keepdim(2) = .True.; keepdim(5) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, qs_comm%value, ierr); qs_comm%me = xmpi_comm_rank(qs_comm%value)

 ! Create communicator for the (perturbation, k-point, band_sum)
 keepdim = .False.; keepdim(1) = .True.; keepdim(3:4) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, pkb_comm%value, ierr); pkb_comm%me = xmpi_comm_rank(pkb_comm%value)

 ABI_FREE(dims)
 ABI_FREE(periods)
 ABI_FREE(keepdim)
 call xmpi_comm_free(comm_cart)
#endif

 ! Distribute spins and q-points in IBZ.
 !call xmpi_split_cyclic(nsppol, spin_comm%value, gams%my_nspins, gams%my_spins)
 !ABI_CHECK(gams%my_nspins > 0, sjoin("nsppol (", itoa(nsppol), ") < spin_comm_nproc (", itoa(spin_comm%nproc), ")"))

 !call xmpi_split_cyclic(gams%nqibz, qpt_comm%value, gams%my_nqibz, gams%my_iqibz)
 !ABI_CHECK(gams%my_nqibz > 0, sjoin("nqibz (", itoa(gams%nqibz), ") < qpt_comm_nproc (", itoa(qpt_comm%nproc), ")"))

 if (my_rank == master) then
   write(std_out, "(/,a)")" === MPI parallelism ==="
   !write(std_out, "(2(a,i0))")"P Allocating and summing bands from my_bsum_start: ", self%my_bsum_start, &
   !    " up to my_bsum_stop: ", self%my_bsum_stop
   write(std_out, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", pert_comm%nproc
   write(std_out, "(a,i0)")"P Number of perturbations treated by this CPU: ", my_npert
   write(std_out, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", qpt_comm%nproc
   !write(std_out, "(2(a,i0))")"P Number of q-points in the IBZ treated by this proc: " , &
   !    count(self%itreat_qibz == 1), " of ", self%nqibz
   write(std_out, "(a,i0)")"P Number of CPUs for parallelism over bands: ", bsum_comm%nproc
   write(std_out, "(a,i0)")"P Number of CPUs for parallelism over spins: ", spin_comm%nproc
   write(std_out, "(a,i0)")"P Number of CPUs for parallelism over k-points: ", kpt_comm%nproc
   !write(std_out, "(2(a,i0),/)")"P Number of k-point in Sigma_nk treated by this proc: ", self%my_nkcalc, " of ", self%nkcalc

   ! Master creates the netcdf file used to store the results of the calculation.
 end if

 !call xmpi_barrier(comm)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 if (pert_comm%nproc > 1) then
   ! Activate parallelism over perturbations
   ! Build table with list of perturbations treated by this CPU inside pert_comm
   call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, pert_comm%value)
   call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, pert_comm%value)
   ABI_FREE(my_pinfo)
   ABI_FREE(pert_table)
 end if

 call wrtout([std_out, ab_out], " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
 ! Prepare Fourier interpolation of DFPT potentials.
 comm_rpt = xmpi_comm_self
 !comm_rpt = bqs_comm%value
 call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions on the FS are stored in wfd.
 ! Need all k-points on the FS because of k+q, spin is not distributed for the time being.
 ! It would be possible to reduce the memory allocated per MPI-rank via OpenMP.

 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz ,nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 do mys=1,gams%my_nspins
   spin = gams%my_spins(mys)
   fs => fstab(spin)
   do ik_bz=1,fs%nkfs
     ik_ibz = fs%indkk_fs(1, ik_bz)
     bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
     bks_mask(bstart_k:bstart_k+nband_k-1, ik_ibz, spin) = .True.
   end do
 end do
 !bks_mask(1:mband,:,:) = .True. ! no memory distribution, each node has the full set of states.

 ! Impose istwfk = 1 for all k-points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 ecut = dtset%ecut
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
   nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
   dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print(header="Wavefunctions on the Fermi surface")

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ! one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2, 3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! that will be used to symmetrize the wavefunctions in G-space.
 call cwtime(cpu, wall, gflops, "start")
 call wrtout(std_out, " Computing mpw. This may take some time for dense k/q meshes...")
 !call fstab_get_mpw_gmax(nsppol, fstab, mpw, gmax, comm)
 mpw = 0; gmax = 0; cnt = 0

 ! FIXME: This is an hotspot due to the loop over nkfs.
 ! Should use a geometric approach to compute mpw gmax.
 do spin=1,nsppol
   fs => fstab(spin)
   do ik_bz=1,fs%nkfs
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm
     kk = fs%kpts(:, ik_bz)
     ! Compute G sphere, returning npw. Note istwfk == 1.
     call get_kg(kk, 1, ecut, cryst%gmet, onpw, gtmp)
     mpw = max(mpw, onpw)
     do ipw=1,onpw
       do ii=1,3
        gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
       end do
     end do
     ABI_FREE(gtmp)

     do iq_ibz=1,gams%nqibz
       qpt = gams%qibz(:,iq_ibz)
       ! TODO: g0 umklapp here can enter into play!
       ! fstab should contains the max of the umlapp G-vectors.
       ! gmax could not be large enough!
       kq = kk + qpt
       call get_kg(kq, 1, ecut, cryst%gmet, onpw, gtmp)
       mpw = max(mpw, onpw)
       do ipw=1,onpw
         do ii=1,3
          gmax(ii) = max(gmax(ii), abs(gtmp(ii,ipw)))
         end do
       end do
       ABI_FREE(gtmp)
     end do

   end do ! ik_bz
 end do ! spin

 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)
 my_gmax = gmax; call xmpi_max(my_gmax, gmax, comm, ierr)
 call wrtout(std_out, sjoin(' Optimal value of mpw: ', itoa(mpw)))
 call cwtime_report(" gmax and mpw", cpu, wall, gflops)

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

 ! TODO FOR PAW
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1    ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2       ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0  ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2, usevnl))
 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq, psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate vlocal. Note nvloc
 ! I set vlocal to huge to trigger possible bugs (DFPT routines should not access the data)
 ABI_MALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))
 vlocal = huge(one)

 ! Allocate work space arrays.
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red, (2, 3, cryst%natom, natom3))

 ABI_CALLOC(dummy_vtrial, (nfftf, nspden))
 ! TODO: Save data to netcdf file for each q in IBZ and then read data to build a2Fw once all big
 ! datastructures (wfd, dvdb) have been deallocated.
 ! As a side effect, one can also implement restart over q-points

 ! Create ddkop object to compute group velocities if needed.
 !
 !   1) precompute group velocities in the IBZ and the ihave_ikibz_spin file (common to all procs)
 !   2) Use symmetries to reconstruct v_kq from vcar_ibz
 !
 ! NB: All procs store in memory the same set of Bloch states inside the energy window.

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 call cwtime(cpu, wall, gflops, "start", msg=" Computing v_nk matrix elements for all states on the FS...")
 ii = huge(1); jj = -1
 do spin=1,nsppol
   ii = min(ii, fstab(spin)%bmin)
   jj = max(jj, fstab(spin)%bmax)
 end do
 ABI_CALLOC(vcar_ibz, (3, ii:jj, nkibz, nsppol))
 ABI_MALLOC(cgwork, (2, mpw * wfd%nspinor))

 cnt = 0
 do spin=1,nsppol
   fs => fstab(spin)
   do ik_ibz=1,ebands%nkpt
     kk = ebands%kptns(:, ik_ibz)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     ! NB: The two checks below are global --> all procs will cycle.
     if (all(bks_mask(:, ik_ibz, spin) .eqv. .False.)) cycle
     if (npw_k == 1) cycle
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do band_k=fs%bmin,fs%bmax
       if (.not. bks_mask(band_k, ik_ibz, spin)) cycle
       !if (.not. wfd%ihave_ug(band_k, ik_ibz, spin)) cycle
       call wfd%copy_cg(band_k, ik_ibz, spin, cgwork)
       eig0nk = ebands%eig(band_k, ik_ibz, spin)
       vk = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
       vcar_ibz(:, band_k, ik_ibz, spin) = vk

     end do
   end do
 end do ! spin

 call xmpi_sum(vcar_ibz, comm, ierr)
 ABI_FREE(cgwork)
 call cwtime_report(" Velocities", cpu, wall, gflops)

 ABI_FREE(bks_mask)

 ! Build krank object to find k-points
 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 ! Loop over my q-points in the IBZ.
 do imyq=1,gams%my_nqibz
   iq_ibz = gams%my_iqibz(imyq)

   call cwtime(cpu_q, wall_q, gflops_q, "start")

   qpt = gams%qibz(:, iq_ibz)
   msg = sjoin("[", itoa(iq_ibz), "/", itoa(gams%nqibz), "]")
   call wrtout(std_out, sjoin(" Computing phonon linewidths for IBZ q-point:", ktoa(qpt), msg))

   ! Use Fourier interpolation of DFPT potentials to get my_npert potentials.
   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfft, nspden, dvdb%my_npert))
   call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)
   !end if

   ! Examine the symmetries of the q wavevector.
   call littlegroup_q(cryst%nsym, qpt, symq, cryst%symrec, cryst%symafm, timerev_q, prtvol=dtset%prtvol)

   ! Get phonon frequencies and eigenvectors for this q-point.
   call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

   ! Allocate vlocal1 with correct cplex. Note nvloc and my_npert.
   ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert))

   ! Loop over my spins.
   do mys=1,gams%my_nspins
     spin = gams%my_spins(mys)
     fs => fstab(spin)

     ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
     do imyp=1,my_npert
       call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc,&
                 pawfgr, mpi_enreg, dummy_vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))
     end do

     ! Continue to initialize the GS Hamiltonian
     call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ! maxnb is the maximum number of bands crossing the FS, used to dimension arrays.
     mnb = fs%maxnb
     ABI_MALLOC(bras_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(kets_k, (2, mpw*nspinor, mnb))
     ABI_MALLOC(h1kets_kq, (2, mpw*nspinor, mnb))
     ABI_MALLOC(gkk_atm, (2, mnb, mnb, natom3))
     !ABI_MALLOC(gkq_nu, (2, mnb, natom3))

     ! =====================================
     ! Integration over the FS for this spin
     ! =====================================
     ! Compute integration weights and distribute k-points (gams%my_nfsk_q)
     !call phgamma_setup_qpoint(gams, fs, cryst, ebands, spin, ltetra, qpt, nesting, kpt_comm%value)

     do myik=1,gams%my_nfsk_q
       call cwtime(cpu_k, wall_k, gflops_k, "start")

       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       ik_fs = gams%my_ifsk_q(myik)
       kk = fs%kpts(:, ik_fs)
       ik_ibz = fs%indkk_fs(1, ik_fs); isym_k = fs%indkk_fs(2, ik_fs)
       trev_k = fs%indkk_fs(6, ik_fs); g0_k = fs%indkk_fs(3:5,ik_fs)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kk_ibz = ebands%kptns(:,ik_ibz)

       ! Number of bands crossing the Fermi level at k
       bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)

       ! Find k+q in the extended zone and extract symmetry info. cycle if k+q not in FS.
       ! Be careful here because there are two umklapp vectors to be considered:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt; ikq_bz = fs%findkg0(kq, g0bz_kq)

       ! Skip this point if kq does not belong to the FS window.
       if (ikq_bz == -1) cycle

       call krank%get_mapping(1, kq, dksqmax, cryst%gmet, indkk_kq, cryst%nsym, cryst%symafm, cryst%symrel, timrev, &
                              use_symrec=.False.)

       if (dksqmax > tol12) then
         write(msg, '(3a,es16.6,6a)' ) &
          "The WFK file cannot be used to compute phonon linewidths.",ch10, &
          "At least one of the k-points on the FS could not be generated from a symmetrical one. dksqmax: ", dksqmax,ch10, &
          "Q-mesh: ", trim(ltoa(gamma_ngqpt)), ", K-mesh (from kptrlatt) ", trim(ltoa(get_diag(ebands%kptrlatt))), &
          'Action: check your WFK file and the (k, q) point input variables.'
          ABI_ERROR(msg)
       end if

       ikq_ibz = indkk_kq(1, 1); isym_kq = indkk_kq(2, 1)
       trev_kq = indkk_kq(6, 1); g0_kq = indkk_kq(3:5, 1)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)

       ! If we have used the KERANGE trick, we may have k or k+q points with just one G component set to zero
       ! so we skip this transition immediately. This should happen only if fsewin > sigma_erange.
       if (wfd%npwarr(ik_ibz) == 1 .or. wfd%npwarr(ikq_ibz) == 1) cycle

       ! Number of bands crossing the Fermi level at k+q
       bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz)
       ABI_CHECK(nband_k <= mnb .and. nband_kq <= mnb, "wrong nband")

       ! Get npw_k, kg_k and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kk, kk_ibz, bstart_k, nband_k, spin, mpw, fs%indkk_fs(:,ik_fs), cryst, &
                          work_ngfft, work, istwf_k, npw_k, kg_k, kets_k)

       ! Get npw_kq, kg_kq and symmetrize wavefunctions from IBZ (if needed).
       call wfd%sym_ug_kg(ecut, kq, kq_ibz, bstart_kq, nband_kq, spin, mpw, indkk_kq(:,1), cryst, &
                          work_ngfft, work, istwf_kq, npw_kq, kg_kq, bras_kq)

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm == 1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1,&
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       ! Compute k+G vectors
       nkpg = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk, nkpg, npw_k)

       ! Compute nonlocal form factors ffnlk at (k+G)
       ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnlk, psps%ffspl,&
                   cryst%gmet, cryst%gprimd, ider0, idir0, psps%indlmn, kg_k, kpg_k, kk, psps%lmnmax, &
                   psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, psps%ntypat, &
                   psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum, &
                   comm=pert_comm%value)

       ! Compute k+q+G vectors
       nkpg1 = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnl1, psps%ffspl, cryst%gmet, cryst%gprimd, ider0, idir0, &
                   psps%indlmn, kg_kq, kpg1_k, kq, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg1, &
                   npw_kq, psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_kq, ylmgr_kq, &
                   comm=pert_comm%value)

       ! Loop over all my atomic perturbations and compute gkk_atm.
       gkk_atm = zero
       do imyp=1,my_npert
         idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &                    ! In
                             cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
                             npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
                             dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
                             reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

         ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do ib_k=1,nband_k
           band_k = ib_k + bstart_k - 1
           eig0nk = ebands%eig(band_k, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, kets_k(:,:,ib_k), cwaveprj0, h1kets_kq(:,:,ib_k), &
                        grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
                        optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c, usevnl)
         end do

         call rf_hamkq%free()

         ABI_FREE(kinpw1)
         ABI_FREE(dkinpw)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)

         ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
         ! No need to handle istwf_kq because it's always 1.
         ! The array eig1_k contains:
         !
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
         do ib_k=1,nband_k
           do ib_kq=1,nband_kq
             gkk_atm(:, ib_kq, ib_k, ipc) = cg_zdotc(npw_kq*nspinor, bras_kq(1,1,ib_kq), h1kets_kq(1,1,ib_k))
           end do
         end do

       end do ! imyp (loop over my_npert atomic perturbations)

       ABI_FREE(gs1c)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)

       ! Collect gkk_atm inside pert_comm so that all procs can operate on the data.
       if (pert_comm%nproc > 1) call xmpi_sum(gkk_atm, pert_comm%value, ierr)

       ! Get gkq in the phonon representation.
       !call ephtk_gkknu_from_atm(mnb, mnb, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu)

       ! Compute group velocities if we are in transport mode or adaptive gaussian or
       ! tetrahedron with libtetrabz returning nesting condition.
       need_velocities = .False.

       if (need_velocities) then
         ! Compute diagonal matrix elements of velocity operator with DFPT routines
         ! Velocities are in Cartesian coordinates.
         !
         ! If k+q is not in the IBZ, we need to recostruct the value by symmetry using v(Sq) = S v(q).
         ! Use transpose(R) because we are using the tables for the wavefunctions
         ! In this case listkk has been called with symrec and use_symrec=False
         ! so q_bz = S^T q_ibz where S is the isym_kq symmetry

         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, kg_k)
         do ib_k=1,nband_k
           band_k = ib_k + bstart_k - 1
           vk = vcar_ibz(:, band_k, ik_ibz, spin)
           if (.not. isirr_k) then
             vk = matmul(transpose(cryst%symrel_cart(:,:,isym_k)), vk)
             if (trev_k /= 0) vk = -vk
           end if
           !vk = ddkop%get_vdiag(ebands%eig(band_k, ik_ibz, spin), &
           !                     istwf_k, npw_k, wfd%nspinor, kets_k(:,:,ib_k), cwaveprj0)
           fs%vk(:,ib_k) = vk
         end do

         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kq, istwf_kq, npw_kq, kg_kq)
         do ib_kq=1,nband_kq
           band_kq = ib_kq + bstart_kq - 1
           vkq = vcar_ibz(:, band_kq, ikq_ibz, spin)
           if (.not. isirr_kq) then
             vkq = matmul(transpose(cryst%symrel_cart(:,:,isym_kq)), vkq)
             if (trev_kq /= 0) vkq = -vkq
           end if
           !vkq = ddkop%get_vdiag(ebands%eig(band_kq, ikq_ibz, spin), &
           !                                 istwf_kq, npw_kq, wfd%nspinor, bras_kq(:,:,ib_kq), cwaveprj0)
           fs%vkq(:,ib_kq) = vkq
         end do
       end if

       if (myik < 20 .or. (fs%nkfs > 100 .and. mod(myik, 200) == 0)) then
         write(msg,'(4(a,i0),a,f8.2)')" q-point [", iq_ibz, "/", gams%nqibz, "] k-point [", myik, "/", gams%my_nfsk_q, "]"
         call cwtime_report(msg, cpu_k, wall_k, gflops_k)
       end if

     end do ! ik_fs: sum over k-points on the BZ FS for this spin.

     ABI_FREE(bras_kq)
     ABI_FREE(kets_k)
     ABI_FREE(h1kets_kq)
     ABI_FREE(gkk_atm)
     !ABI_FREE(gkq_nu)

     ! Save results in gams.
     !gams%vals_qibz(:,:,:,iq_ibz, spin) = tgam
   end do ! spin

   ABI_FREE(v1scf)
   ABI_FREE(vlocal1)

   write(msg,'(2(a,i0),a)')" Computation of q-point [", iq_ibz, "/", gams%nqibz, "]"
   call cwtime_report(msg, cpu_q, wall_q, gflops_q, end_str=ch10)
 end do ! iq_ibz

 call cwtime_report(" phonon linewidths k-loop", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

 ! Free memory
 ABI_FREE(gvnlx1)
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
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_SFREE(vcar_ibz)
 call krank%free()

 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call ddkop%free()
 call gs_hamkq%free()
 call wfd%free()
 !do spin=1,ebands%nsppol
 !  call fstab(spin)%free()
 !end do
 !ABI_FREE(fstab)

 ! Deallocate MPI communicators.
 ! TODO Use gstore communicators
 !call pert_comm%free()
 !call qpt_comm%free()
 !call bsum_comm%free()
 !call qs_comm%free()
 !call kpt_comm%free()
 !call spin_comm%free()
 !call pkb_comm%free()

end subroutine gstore_compute
!!***

!----------------------------------------------------------------------

end module m_gstore
!!***
