!!****m* ABINIT/m_ddk
!! NAME
!!  m_ddk
!!
!! FUNCTION
!!  Objects and methods to extract data from DDK files.
!!  The DDK files are binary (soon also netcdf) files with Hamiltonian derivatives
!!  wrt k, and the corresponding wave functions
!!
!! COPYRIGHT
!! Copyright (C) 2016-2020 ABINIT group (MJV, HM, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ddk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr
 use m_dtset
 use m_krank
 use m_crystal
 use m_mpinfo
 use m_cgtools
 use m_hamiltonian
 use m_initylmg
 use m_ebands
 use m_pawcprj
 use m_getgh1c
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,      only : strcat, sjoin, itoa, ktoa
 use m_io_tools,      only : iomode_from_fname
 use m_time,          only : cwtime, cwtime_report
 use defs_abitypes,   only : MPI_type
 use defs_datatypes,  only : ebands_t, pseudopotential_type
 use m_vkbr,          only : vkbr_t, nc_ihr_comm, vkbr_init, vkbr_free
 use m_pawtab,        only : pawtab_type
 use m_wfk,           only : wfk_read_ebands !, wfk_read_h1mat
 use m_wfd,           only : wfd_t, wfd_init, wave_t

 implicit none

 private
!!***

 public :: ddk_red2car           ! Convert band velocities from cartesian to reduced coordinates

!!***

 type, private :: ham_targets_t
   real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:)
   real(dp),allocatable :: kpg_k(:,:), kpg1_k(:,:)
   real(dp),allocatable :: ph3d(:,:,:), ph3d1(:,:,:)
   real(dp),allocatable :: dkinpw(:), kinpw1(:)
   contains
     procedure :: free => ham_targets_free   ! Free memory.
 end type ham_targets_t


!!****t* m_ddk/ddkop_t
!! NAME
!!  ddkop_t
!!
!! FUNCTION
!!  This object provides a simplified interface to compute matrix elements of the
!!  velocity operator with the DFPT routines.
!!
!! SOURCE

 type,public :: ddkop_t

  integer :: ipert
  ! Perturbation type: natom + 1

  integer :: inclvkb
  ! Option for calculating the matrix elements of [Vnl,r].
  ! 0 to exclude commutator, 2 to include it

  integer :: usepaw
  ! 0 for NC, 1 for PAW

  integer :: mpw
  ! Maximum number of plane-waves over k-points (used to dimension arrays)

  real(dp) :: kpoint(3)
  ! K-point (set in setup_spin_kpoint)

  real(dp) :: eig0nk

  real(dp) :: dfpt_sciss = zero

  real(dp) :: rprimd(3,3)

  type(MPI_type),pointer :: mpi_enreg => null()

  type(gs_hamiltonian_type) :: gs_hamkq(3)

  type(rf_hamiltonian_type) :: rf_hamkq(3)

  type(ham_targets_t), private :: htg(3)
  ! Store arrays targetted by the hamiltonians.

  real(dp), private, allocatable :: gh1c(:,:,:)
   !gh1c, (2, mpw*nspinor, 3))

  real(dp), private, allocatable :: gs1c(:,:,:)
   ! gs1c, (2, mpw*nspinor, 3*psps%usepaw))

 contains

   procedure :: setup_spin_kpoint => ddkop_setup_spin_kpoint
    ! Prepare application of dH/dk for given spin, k-point.

   procedure :: apply => ddkop_apply
    ! Apply dH/dk to input wavefunction.

   procedure :: get_braket => ddkop_get_braket
    ! Compute matrix element (complex results) in cartesian coords.

   procedure :: get_vdiag => ddkop_get_vdiag
    ! Compute diagonal matrix element (real) in cartesian coords.

   procedure :: free => ddkop_free
    ! Free memory.

 end type ddkop_t

 public :: ddkop_new   ! Build object

 !interface ddkop_t
 !  procedure ddkop_new
 !end interface ddkop_t
!!***

!!****t* m_ddk/ddkstore_t
!! NAME
!!  ddkstore_t
!!
!! FUNCTION
!!  This object stores the matrix elements of the velocity operator computed with the DFPT routines.
!!
!! SOURCE

 type,public :: ddkstore_t

   integer :: bmin = 1, bmax = -1
   ! Min and max band index

   character(len=50) :: mode = "reduced"
    ! "cart" or "reduced"

   logical :: only_diago = .False.
   ! True if we are computing only the diagonal elements

   real(dp),allocatable :: dipoles(:,:,:,:,:,:)
    ! (3, 2, mband, mband, nkpt, nsppol))

   real(dp),allocatable :: vdiago(:,:,:,:)
    ! (3, bmin:bmax, nkpt, nsppol)

   real(dp),allocatable :: vmat(:,:,:,:,:,:)
    ! (2, 3, bmin:bmax, bmin:bmax, nkpt, nsppol))

 contains

   procedure :: compute_ddk => ddkstore_compute_ddk
    ! Calculate DDK matrix elements (diago or full b,b' matrix).
    ! Return results in datatype. Optionally, save results to disk in EVK format.

   procedure :: free => ddkstore_free
    ! Free memory.

 end type ddkstore_t
!!***

 !public :: ddkop_new   ! Build object

CONTAINS

!----------------------------------------------------------------------

!!****f* m_ddk/ddkstore_compute_ddk
!! NAME
!!  ddkstore_compute_ddk
!!
!! FUNCTION
!!  Calculate the DDK matrix elements using the commutator formulation.
!!
!! INPUTS
!!  prefix: Prefix for output EVK file. Empty if output files are not wanted
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkstore_compute_ddk(ds, wfk_path, prefix, dtset, psps, pawtab, ngfftc, comm)

!Arguments ------------------------------------
!scalars
 class(ddkstore_t),intent(inout) :: ds
 character(len=*),intent(in) :: wfk_path, prefix
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!arrays
 integer,intent(in) :: ngfftc(18)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: mband, nbcalc, nsppol, ib_v, ib_c, mpw, spin, nspinor, nkpt, nband_k, npw_k
 integer :: ii, ik, bmin, bmax, istwf_k, idir, my_rank, nproc, ierr, bstop
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
#ifdef HAVE_NETCDF
 integer :: ncerr, ncid
#endif
 character(len=500) :: msg
 character(len=fnlen) :: fname
 logical :: write_ncfile
 type(wfd_t) :: wfd
 type(vkbr_t) :: vkbr
 type(ebands_t) :: ebands
 type(crystal_t) :: cryst
 type(hdr_type) :: tmp_hdr, hdr
 type(ddkop_t) :: ddkop
 type(wave_t),pointer :: wave_v, wave_c
!arrays
 integer,allocatable :: distrib_mat(:,:,:,:), distrib_diago(:,:,:), nband(:,:), kg_k(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
 real(dp) :: kpt(3), vv(2, 3)
 real(dp),allocatable :: cg_c(:,:), cg_v(:,:)
 complex(dpc) :: vg(3), vr(3)
 complex(gwpc),allocatable :: ihrc(:,:), ug_c(:), ug_v(:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 if (my_rank == master) call wrtout([std_out, ab_out], " Computation of velocity matrix elements (DDK)", newlines=1)

 ABI_CHECK(psps%usepaw == 0, "PAW not implemented")

 ! Get ebands and hdr from WFK file.
 ebands = wfk_read_ebands(wfk_path, comm, out_hdr=hdr)
 cryst = hdr%get_crystal()

 ! Extract important dimensions from hdr%
 nkpt    = hdr%nkpt
 nsppol  = hdr%nsppol
 nspinor = hdr%nspinor
 mband   = hdr%mband

 ! Define band range
 ! TODO: Perhaps one should allocate output arrays using bmin:bmax
 ! and allow for nc output only if bmin == 1 and bmax == mband
 if (ds%bmax == -1) ds%bmax = mband

 if (ds%bmin < 1 .or. ds%bmin > mband .or. ds%bmax > mband .or. ds%bmin > ds%bmax) then
   MSG_ERROR(sjoin("Invalid value for bmin, bmax", itoa(ds%bmin), itoa(ds%bmax), "with mband:", itoa(mband)))
 end if

 bmin = ds%bmin; bmax = ds%bmax
 nbcalc  = bmax - bmin + 1
 write_ncfile = len_trim(prefix) > 0
 if (write_ncfile .and. .not. (bmin == 1 .and. bmax == mband) ) then
   write_ncfile = .False.
   MSG_WARNING("Cannot write ncfile if .not. (bmin == 1 .and. bmax == mband)")
 end if

 if (my_rank == master) then
   write(ab_out, "(a)")" Parameters extracted from the Abinit header:"
   write(ab_out, "(a, f5.1)") '    ecut:    ', hdr%ecut
   write(ab_out, "(a, i0)")   '    nkpt:    ', nkpt
   write(ab_out, "(a, i0)")   '    mband:   ', mband
   write(ab_out, "(a, i0)")   '    nsppol:  ', nsppol
   write(ab_out, "(a, i0)")   '    nspinor: ', nspinor
   write(ab_out, "(a, i0)")   '    useylm:  ', dtset%useylm
   write(ab_out, "(a, i0)")   '    inclvkb: ', dtset%inclvkb
   write(ab_out, "(2(a, i0))")'    bmin: ', bmin, ", bmax: ", bmax
   if (ds%only_diago) then
     write(ab_out, "(a)")'    Computing diagonal matrix elements only'
   else
     write(ab_out, "(a)")'    Computing diagonal and off-diagonal matrix elements'
   end if
   write(ab_out, "(2(a, i0))")'    Between band index bmin: ', bmin, ", bmax: ", bmax
   write(ab_out, "(a)")""
 end if

 ! Create distribution of the wavefunctions mask.
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkpt, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkpt, nsppol))
 keep_ur = .false.; bks_mask = .false.; nband = mband

 if (ds%only_diago) then
   ! Distribute k-points, spin and (b, b) diagonal over MPIR processors.
   ABI_MALLOC(distrib_diago, (bmin:bmax, nkpt, nsppol))
   distrib_diago = -1

   ! Create bks_mask to load the wavefunctions.
   ii = 0
   do spin=1,nsppol
     do ik=1,nkpt
       do ib_v=bmin,bmax
          ii = ii + 1; if (mod(ii, nproc) /= my_rank) cycle ! MPI parallelism.
          distrib_diago(ib_v, ik, spin) = my_rank
          bks_mask(ib_v, ik, spin) = .true.
       end do
     end do
   end do
   call wrtout(std_out, sjoin(" Rank: ", itoa(my_rank), "will treat", itoa(count(distrib_diago == my_rank))))

 else
   ! Distribute k-points, spin and (b, b') pairs over the processors
   ABI_MALLOC(distrib_mat, (bmin:bmax, bmin:bmax, nkpt, nsppol))
   call xmpi_distab(nproc, distrib_mat)

   ! Create bks_mask to load the wavefunctions
   do spin=1,nsppol
     do ik=1,nkpt
       ! Loop over v bands
       do ib_v=bmin,bmax
        ! Loop over c bands
         do ib_c=bmin,bmax
           if (distrib_mat(ib_c, ib_v, ik, spin) == my_rank) then
             bks_mask(ib_v, ik, spin) = .true.
             bks_mask(ib_c, ik, spin) = .true.
           end if
         end do
       end do
     end do
   end do

   call wrtout(std_out, sjoin(" Rank: ", itoa(my_rank), "will treat", itoa(count(distrib_mat == my_rank))))
 end if

 ! Initialize distributed wavefunctions object
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkpt, nsppol,&
   bks_mask, dtset%nspden, nspinor, hdr%ecut, dtset%ecutsm, dtset%dilatmx, ebands%istwfk, ebands%kptns,&
   ngfftc, dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband)

 call wfd%print(header="Wavefunctions on the k-points grid")

 ! Read wavefunctions from WFK file.
 call wfd%read_wfk(wfk_path, iomode_from_fname(wfk_path))

 ! Allocate workspace arrays
 mpw = maxval(wfd%npwarr)
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(ug_c, (mpw*nspinor))
 ABI_MALLOC(ug_v, (mpw*nspinor))
 if (dtset%useria /= 666) then
   ABI_MALLOC(cg_c, (2, mpw*nspinor))
   ABI_MALLOC(cg_v, (2, mpw*nspinor))
 end if

 ABI_MALLOC(cwaveprj, (0, 0))
 ABI_CALLOC(ds%dipoles, (3, 2, bmin:bmax, bmin:bmax, nkpt, nsppol))
 ABI_MALLOC(ihrc, (3, nspinor**2))

 if (ds%only_diago) then
   ABI_CALLOC(ds%vdiago, (3, bmin:bmax, nkpt, nsppol))
 else
   ABI_CALLOC(ds%vmat, (2, 3, bmin:bmax, bmin:bmax, nkpt, nsppol))
 end if

 if (dtset%useria /= 666) then
   ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
   !if (my_rank == master) call ddkop%print(ab_out)
 end if

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 do spin=1,nsppol
   do ik=1,nkpt

     ! Only do a subset a k-points
     if (ds%only_diago) then
       if (all(distrib_diago(:, ik, spin) /= my_rank)) cycle
     else
       if (all(distrib_mat(bmin:bmax, bmin:bmax, ik, spin) /= my_rank)) cycle
     end if
     call cwtime(cpu, wall, gflops, "start")

     nband_k  = wfd%nband(ik, spin)
     istwf_k  = wfd%istwfk(ik)
     npw_k    = wfd%npwarr(ik)
     kpt      = wfd%kibz(:,ik)
     kg_k(:,1:npw_k) = wfd%kdata(ik)%kg_k

     if (dtset%useria /= 666) then
       call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kpt, istwf_k, npw_k, kg_k)
     else
       ! Allocate KB form factors
       ! Prepare term i <n,k|[Vnl,r]|n"k>
       if (dtset%inclvkb /= 0) call vkbr_init(vkbr, cryst, psps, dtset%inclvkb, istwf_k, npw_k, kpt, kg_k)
     end if

     ! Loop over bands
     do ib_v=bmin,bmax
       if (ds%only_diago) then
         if (distrib_diago(ib_v,ik,spin) /= my_rank) cycle
       else
         if (all(distrib_mat(:,ib_v,ik,spin) /= my_rank)) cycle
       end if

       if (dtset%useria /= 666) then
         call wfd%copy_cg(ib_v, ik, spin, cg_v)
         call ddkop%apply(ebands%eig(ib_v, ik, spin), npw_k, wfd%nspinor, cg_v, cwaveprj)
       else
         ABI_CHECK(wfd%get_wave_ptr(ib_v, ik, spin, wave_v, msg) == 0, msg)
         ug_v(1:npw_k*nspinor) = wave_v%ug
       end if

       ! Loop over bands
       bstop = bmax; if (ds%only_diago) bstop = ib_v
       do ib_c=ib_v,bstop
         if (.not. ds%only_diago) then
           if (distrib_mat(ib_c, ib_v, ik, spin) /= my_rank) cycle
         end if

         if (dtset%useria /= 666) then
           call wfd%copy_cg(ib_c, ik, spin, cg_c)
           vv = ddkop%get_braket(ebands%eig(ib_c, ik, spin), istwf_k, npw_k, nspinor, cg_c, mode=ds%mode)
           !if (ib_v == ib_c) vv(2, :) = zero

           if (ds%only_diago) then
             ds%vdiago(:,ib_c,ik,spin) = vv(1, :)
           else
             ds%vmat(:,:,ib_c,ib_v,ik,spin) = vv
             ! Hermitian conjugate
             if (ib_v /= ib_c) then
               ds%vmat(1,:,ib_v,ib_c,ik,spin) =  vv(1, :)
               ds%vmat(2,:,ib_v,ib_c,ik,spin) = -vv(2, :)
             end if
           end if

           do idir=1,3
             ds%dipoles(idir,:,ib_c,ib_v,ik,spin) = vv(:, idir)
             ! Hermitian conjugate
             if (ib_v /= ib_c) ds%dipoles(idir,:,ib_v,ib_c,ik,spin) = [vv(1, idir), -vv(2, idir)]
           end do

         else
           ABI_CHECK(wfd%get_wave_ptr(ib_c, ik, spin, wave_c, msg) == 0, msg)
           ug_c(1:npw_k*nspinor) = wave_c%ug

           ! Calculate matrix elements of i[H,r] for NC pseudopotentials.
           ihrc = nc_ihr_comm(vkbr, cryst, psps, npw_k, nspinor, istwf_k, dtset%inclvkb, kpt, ug_c, ug_v, kg_k)

           ! HM: 24/07/2018
           ! Transform dipoles to be consistent with results from DFPT
           ! Perturbations with DFPT are along the reciprocal lattice vectors
           ! Perturbations with Commutator are along real space lattice vectors
           ! dot(A, DFPT) = X
           ! dot(B, COMM) = X
           ! B = 2 pi (A^{-1})^T =>
           ! dot(B^T B,COMM) = 2 pi DFPT
           vr = (2*pi)*(2*pi)*sum(ihrc(:,:),dim=2)
           vg(1) = dot_product(Cryst%gmet(1,:), vr)
           vg(2) = dot_product(Cryst%gmet(2,:), vr)
           vg(3) = dot_product(Cryst%gmet(3,:), vr)

           ! Save matrix elements of i*r in the IBZ
           ds%dipoles(:,1,ib_c,ib_v,ik,spin) = real(vg, kind=dp)
           ds%dipoles(:,1,ib_v,ib_c,ik,spin) = real(vg, kind=dp) ! Hermitian conjugate
           if (ib_v == ib_c) then
             ds%dipoles(:,2,ib_c,ib_v,ik,spin) = zero
             ds%dipoles(:,2,ib_v,ib_c,ik,spin) = zero
           else
             ds%dipoles(:,2,ib_c,ib_v,ik,spin) =  aimag(vg)
             ds%dipoles(:,2,ib_v,ib_c,ik,spin) = -aimag(vg) ! Hermitian conjugate
           end if
         end if

       end do
     end do

     ! Free KB form factors
     call vkbr_free(vkbr)

     if (nkpt < 1000 .or. (nkpt > 1000 .and. mod(ik, 200) == 0) .or. ik <= nproc) then
       write(msg,'(2(a,i0),a)')" k-point [", ik, "/", nkpt, "]"
       call cwtime_report(msg, cpu, wall, gflops)
     end if

   end do ! k-points
 end do ! spin

 call cwtime_report(msg, cpu_all, wall_all, gflops_all)

 ABI_FREE(ug_c)
 ABI_FREE(ug_v)
 ABI_FREE(kg_k)
 ABI_FREE(ihrc)
 ABI_FREE(cwaveprj)
 ABI_SFREE(distrib_mat)
 ABI_SFREE(distrib_diago)

 if (dtset%useria /= 666) then
   ABI_FREE(cg_c)
   ABI_FREE(cg_v)
   call ddkop%free()
 end if

 ! Gather the k-points computed by all processes
 call xmpi_sum_master(ds%dipoles, master, comm, ierr)

 if (ds%only_diago) then
   call xmpi_sum_master(ds%vdiago, master, comm, ierr)
 else
   call xmpi_sum_master(ds%vmat, master, comm, ierr)
 end if

 ! Write matrix elements to disk.
#ifdef HAVE_NETCDF

 ! Output EVK file in netcdf format.
 if (my_rank == master .and. write_ncfile) then
   ! Have to build hdr on k-grid with info about perturbation.
   call hdr_copy(hdr, tmp_hdr)
   tmp_hdr%qptn = zero

   !fname = strcat(prefix, "NEW_EVK.nc")
   !call wrtout(ab_out, sjoin("- Writing file: ", fname))
   !NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EVK.nc file")
   !tmp_hdr%pertcase = 0
   !NCF_CHECK(tmp_hdr%ncwrite(ncid, 43, nc_define=.True.))
   !NCF_CHECK(cryst%ncwrite(ncid))
   !NCF_CHECK(ebands_ncwrite(ebands, ncid))
   !if (ds%only_diago) then
   !  ncerr = nctk_def_arrays(ncid, [ &
   !    nctkarr_t('vred_diagonal', "dp", "three, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
   !else
   !  ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vred_matrix', "dp", &
   !      "two, three, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
   !end if
   !NCF_CHECK(ncerr)
   !NCF_CHECK(nctk_set_datamode(ncid))
   !if (ds%only_diago) then
   !  NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vred_diagonal"), ds%vdiago))
   !else
   !  NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vred_matrix"), ds%vmat))
   !end if
   !NCF_CHECK(nf90_close(ncid))

   do ii=1,3
     fname = strcat(prefix, '_', itoa(ii), "_EVK.nc")
     call wrtout(ab_out, sjoin("- Writing EVK file: ", fname, "for reduced direction:", itoa(ii)))
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EVK.nc file")
     tmp_hdr%pertcase = 3 * cryst%natom + ii
     NCF_CHECK(tmp_hdr%ncwrite(ncid, 43, nc_define=.True.))
     NCF_CHECK(cryst%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(ebands, ncid))
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t('h1_matrix_elements', "dp", &
        "two, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "h1_matrix_elements"), ds%dipoles(ii,:,:,:,:,:)))
     NCF_CHECK(nf90_close(ncid))
   end do
   call tmp_hdr%free()
 end if
#endif

 if (my_rank == master .and. dtset%prtvol > 0) then
   write(ab_out, "(2a)")ch10,"Writing velocity matrix elements (only diagonal terms, real part) for testing purpose:"
   do spin=1,nsppol
     do ik=1,min(nkpt, 4)
       write(ab_out, "(2(a,1x,i0),2x,2a)")"For spin: ", spin, ", ikbz: ", ik, ", kpt: ", trim(ktoa(wfd%kibz(:,ik)))
       do ib_c=bmin,min(bmin+8, bmax)
         write(ab_out, "(3(es16.6,2x))") ds%dipoles(:,1,ib_c,ib_c,ik,spin)
       end do
       write(ab_out,*)""
       !do ib_c=bmin,min(bmin+8, bmax)
       !  write(ab_out, "(a, 6(es16.6,2x))")"Sum_k: ", sum(ds%dipoles(:,:,ib_c,ib_c,:,spin), dim=3) / nkpt
       !end do
     end do
   end do
 end if

 ! Free memory
 call wfd%free()
 call ebands_free(ebands)
 call cryst%free()
 call hdr%free()

 ! Block all procs here so that we know output files are available when code returns.
 call xmpi_barrier(comm)

end subroutine ddkstore_compute_ddk
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkstore_free
!! NAME
!!  ddkstore_free
!!
!! FUNCTION
!!  Free memory
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkstore_free(self)

!Arguments ------------------------------------
!scalars
 class(ddkstore_t),intent(inout) :: self

!************************************************************************

 ABI_SFREE(self%vdiago)
 ABI_SFREE(self%vmat)
 ABI_SFREE(self%dipoles)

end subroutine ddkstore_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_red2car
!! NAME
!!  ddk_red2car
!!
!! FUNCTION
!!  Convert ddk matrix element from reduced coordinates to cartesian coordinates.
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

pure subroutine ddk_red2car(rprimd, vred, vcar)

!Arguments -------------------------------------
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: vred(2,3)
 real(dp),intent(out) :: vcar(2,3)

!Local variables -------------------------------
 real(dp) :: vtmp(2,3)

!************************************************************************
 ! vcar = vred; return

 ! Go to cartesian coordinates (same as pmat2cart routine)
 ! V_cart = 1/(2pi) * Rprimd x V_red
 ! where V_red is the derivative computed in the DFPT routines (derivative wrt reduced component).
 vtmp(1,:) = rprimd(:,1)*vred(1,1) &
            +rprimd(:,2)*vred(1,2) &
            +rprimd(:,3)*vred(1,3)
 vtmp(2,:) = rprimd(:,1)*vred(2,1) &
            +rprimd(:,2)*vred(2,2) &
            +rprimd(:,3)*vred(2,3)
 vcar = vtmp / two_pi

end subroutine ddk_red2car
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_new
!! NAME
!!  ddkop_new
!!
!! FUNCTION
!!  Build new object. Use dtset%inclvkb to determine whether non-local part should be included.
!!
!! INPUTS
!! dtset<dataset_type>=All input variables for this dataset.
!! cryst<crystal_t>=Crystal structure.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! mpi_enreg=information about MPI parallelization
!! mpw=Maximum number of plane-waves over k-points.
!! ngfft(18)=contain all needed information about 3D FFT
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ddkop_t) function ddkop_new(dtset, cryst, pawtab, psps, mpi_enreg, mpw, ngfft) result(new)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),target,intent(in) :: mpi_enreg
 integer,intent(in) :: mpw
!arrays
 integer,intent(in) :: ngfft(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1 = 1
 integer :: nfft, mgfft, idir

! *************************************************************************

 ABI_CHECK(dtset%usepaw == 0, "PAW not tested/implemented!")

 new%inclvkb = dtset%inclvkb
 new%usepaw = dtset%usepaw
 new%ipert = cryst%natom + 1
 new%dfpt_sciss = dtset%dfpt_sciss
 new%mpw = mpw

 new%rprimd = cryst%rprimd
 new%mpi_enreg => mpi_enreg

 ! Not used because vlocal1 is not applied.
 nfft = product(ngfft(1:3))
 mgfft = maxval(ngfft(1:3))

 ABI_MALLOC(new%gh1c, (2, new%mpw*dtset%nspinor, 3))
 ABI_MALLOC(new%gs1c, (2, new%mpw*dtset%nspinor, 3))

 do idir=1,3
   ! ==== Initialize most of the Hamiltonian (and derivative) ====
   ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
   ! 2) Perform the setup needed for the non-local factors:
   ! * Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
   ! * PAW: Initialize the overlap coefficients and allocate the Dij coefficients.
   call init_hamiltonian(new%gs_hamkq(idir), psps, pawtab, dtset%nspinor, dtset%nsppol, dtset%nspden, cryst%natom,&
     cryst%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg)
     !paw_ij=paw_ij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
     !usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

   ! Prepare application of the NL part.
   call init_rf_hamiltonian(cplex1, new%gs_hamkq(idir), new%ipert, new%rf_hamkq(idir), has_e1kbsc=.true.)
     !&paw_ij1=paw_ij1,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
     !&mpi_spintab=mpi_enreg%my_isppoltab)
 end do

end function ddkop_new
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_setup_spin_kpoint
!! NAME
!!  ddkop_setup_spin_kpoint
!!
!! FUNCTION
!!  Prepare internal tables that depend on k-point/spin
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t>=Crystal structure.
!!  psps<pseudopotential_type>=Variables related to pseudopotentials.
!!  spin: spin index
!!  kpoint(3): K-point in reduced coordinates.
!!  istwkf_k: defines storage of wavefunctions for this k-point
!!  npw_k: Number of planewaves.
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkop_setup_spin_kpoint(self, dtset, cryst, psps, spin, kpoint, istwf_k, npw_k, kg_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin, npw_k, istwf_k
 class(ddkop_t),intent(inout) :: self
 type(crystal_t) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt1=1, nsppol1=1
 type(mpi_type) :: mpienreg_seq
!arrays
 integer :: npwarr(nkpt1), dummy_nband(nkpt1*nsppol1)
 integer :: idir, nkpg, nkpg1, useylmgr1, optder !, nylmgr1
 real(dp),allocatable :: ylm_k(:,:),ylmgr1_k(:,:,:)

!************************************************************************

 ABI_CHECK(npw_k <= self%mpw, "npw_k > mpw!")
 self%kpoint = kpoint

 ! Set up the spherical harmonics (Ylm) at k+q if useylm = 1
 useylmgr1 = 0; optder = 0
 if (psps%useylm == 1) then
   useylmgr1 = 1; optder = 1
 end if

 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylmgr1_k, (npw_k, 3+6*(optder/2), psps%mpsang**2*psps%useylm*useylmgr1))

 if (psps%useylm == 1) then
   ! Fake MPI_type for sequential part. dummy_nband and nsppol1 are not used in sequential mode.
   call initmpi_seq(mpienreg_seq)
   dummy_nband = 0; npwarr = npw_k
   call initylmg(cryst%gprimd, kg_k, kpoint, nkpt1, mpienreg_seq, psps%mpsang, npw_k, dummy_nband, nkpt1, &
      npwarr, nsppol1, optder, cryst%rprimd, ylm_k, ylmgr1_k)
   call destroy_mpi_enreg(mpienreg_seq)
 end if

 do idir=1,3
   call self%htg(idir)%free()

   ! Continue to initialize the Hamiltonian
   call self%gs_hamkq(idir)%load_spin(spin, with_nonlocal=.true.)
   call self%rf_hamkq(idir)%load_spin(spin, with_nonlocal=.true.)

   ! We need ffnl1 and dkinpw for 3 dirs. Note that the Hamiltonian objects use pointers to keep a reference
   ! to the output results of this routine.
   ! This is the reason why we need to store the targets in self%htg
   call getgh1c_setup(self%gs_hamkq(idir), self%rf_hamkq(idir), dtset, psps, kpoint, kpoint, idir, self%ipert, & ! In
     cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, npw_k, npw_k, &            ! In
     useylmgr1, kg_k, ylm_k, kg_k, ylm_k, ylmgr1_k, &                                       ! In
     self%htg(idir)%dkinpw, nkpg, nkpg1, self%htg(idir)%kpg_k, self%htg(idir)%kpg1_k, &     ! Out
     self%htg(idir)%kinpw1, self%htg(idir)%ffnlk, self%htg(idir)%ffnl1, &                   ! Out
     self%htg(idir)%ph3d, self%htg(idir)%ph3d1)                                             ! Out
 end do

 ABI_FREE(ylm_k)
 ABI_FREE(ylmgr1_k)

end subroutine ddkop_setup_spin_kpoint
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_apply
!! NAME
!!  ddkop_apply
!!
!! FUNCTION
!!  Apply velocity operator dH/dk to wavefunction in G-space. Store results in object.
!!
!! INPUTS
!!  eig0nk: Eigenvalue associated to the wavefunction.
!!  npw_k: Number of planewaves.
!!  nspinor: Number of spinor components.
!!  cwave(2,npw_k*nspinor)=input wavefunction in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!     if not allocated or size=0, they are locally computed (and not sorted)!!
!!
!! SIDE EFFECTS
!! Stores:
!!  gh1c(2,npw1*nspinor)= <G|H^(1)|C> or  <G|H^(1)-lambda.S^(1)|C> on the k+q sphere
!!                        (only kinetic+non-local parts if optlocal=0)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkop_apply(self, eig0nk, npw_k, nspinor, cwave, cwaveprj)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(inout) :: self
 integer,intent(in) :: npw_k, nspinor
 real(dp),intent(in) :: eig0nk
!arrays
 real(dp),intent(inout) :: cwave(2,npw_k*nspinor)
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: berryopt0 = 0, optlocal0 = 0, tim_getgh1c = 1, usevnl0 = 0, opt_gvnlx1 = 0
 integer :: idir, sij_opt, ispinor, ipws, ipw, optnl
 real(dp) :: eshift
!arrays
 real(dp) :: grad_berry(2,(berryopt0/4)), gvnlx1(2,usevnl0)
 real(dp),pointer :: dkinpw(:),kinpw1(:)

!************************************************************************

 self%eig0nk = eig0nk

 if (self%inclvkb /= 0) then
   ! optlocal0 = 0: local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
   ! optnl = 2: non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
   ! opt_gvnlx1 = option controlling the use of gvnlx1 array:
   optnl = 2 !; if (self%inclvkb == 0) optnl = 0

   eshift = self%eig0nk - self%dfpt_sciss
   do idir=1,3
     sij_opt = self%gs_hamkq(idir)%usepaw
     call getgh1c(berryopt0, cwave, cwaveprj, self%gh1c(:,:,idir), &
       grad_berry, self%gs1c(:,:,idir), self%gs_hamkq(idir), gvnlx1, idir, self%ipert, eshift, self%mpi_enreg, optlocal0, &
       optnl, opt_gvnlx1, self%rf_hamkq(idir), sij_opt, tim_getgh1c, usevnl0)
   end do

 else
   ! optnl 0 with DDK does not work as expected.
   ! So I treat the kinetic term explicitly without calling getgh1c.
   do idir=1,3
     kinpw1 => self%gs_hamkq(idir)%kinpw_kp
     dkinpw => self%rf_hamkq(idir)%dkinpw_k
     do ispinor=1,nspinor
       do ipw=1,npw_k
         ipws = ipw + npw_k*(ispinor-1)
         if (kinpw1(ipw) < huge(zero)*1.d-11) then
           self%gh1c(:,ipws,idir) = dkinpw(ipw) * cwave(:,ipws)
         else
           self%gh1c(:,ipws,idir) = zero
         end if
       end do
     end do
   end do
 end if

end subroutine ddkop_apply
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_get_braket
!! NAME
!!  ddkop_get_braket
!!
!! FUNCTION
!!  Compute matrix element in Cartesian coordinates.
!!
!! INPUTS
!!  eig0mk: Eigenvalue associated to the "bra" wavefunction
!!  istwkf_k: defines storage of wavefunctions for this k-point
!!  npw_k: Number of planewaves.
!!  nspinor: Number of spinor components.
!!  brag(2,npw_k*nspinor)=input wavefunction in reciprocal space
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ddkop_get_braket(self, eig0mk, istwf_k, npw_k, nspinor, brag, mode) result(vk)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(in) :: self
 integer,intent(in) :: istwf_k, npw_k, nspinor
 real(dp),intent(in) :: eig0mk
 character(len=*),optional,intent(in) :: mode
!arrays
 real(dp),intent(in) :: brag(2*npw_k*nspinor)
 real(dp) :: vk(2,3)

!Local variables-------------------------------
!scalars
 integer :: idir
 real(dp) :: doti
!arrays
 real(dp) :: dotarr(2), vk_red(2, 3)
 character(len=50) :: my_mode

!************************************************************************

 if (self%usepaw == 0) then
   ! <u_(iband,k+q)^(0)|H_(k+q,k)^(1)|u_(jband,k)^(0)>  (NC psps)
   do idir=1,3
     dotarr = cg_zdotc(npw_k * nspinor, brag, self%gh1c(:,:,idir))
     if (istwf_k > 1) then
       !dum = two * j_dpc * AIMAG(dum); if (vkbr%istwfk==2) dum = dum - j_dpc * AIMAG(gamma_term)
       doti = two * dotarr(2)
       if (istwf_k == 2 .and. self%mpi_enreg%me_g0 == 1) then
         ! nspinor always 1
         ! TODO: Recheck this part but it should be ok.
         doti = doti - (brag(1) * self%gh1c(2,1,idir) - brag(2) * self%gh1c(1,1,idir))
       end if
       dotarr(2) = doti; dotarr(1) = zero
     end if
     vk(:, idir) = dotarr
   end do
 else
   MSG_ERROR("PAW Not Implemented")
   ! <u_(iband,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(jband,k)^(0)> (PAW)
   ! eshiftkq = half * (eig0mk - self%eig0nk)
   ABI_UNUSED(eig0mk)
 end if

 my_mode = "cart"; if (present(mode)) my_mode = mode
 select case (mode)
 case ("cart")
   vk_red = vk
   call ddk_red2car(self%rprimd, vk_red, vk)
 case ("reduced")
   continue
 case default
   MSG_ERROR(sjoin("Invalid vaue for mode:", mode))
 end select

end function ddkop_get_braket
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_get_vdiag
!! NAME
!!  ddkop_get_vdiag
!!
!! FUNCTION
!!  Simplified interface to compute the diagonal matrix element of the velocity operator in cartesian coords.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ddkop_get_vdiag(self, eig0nk, istwf_k, npw_k, nspinor, cwave, cwaveprj, mode) result(vk)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(inout) :: self
 integer,intent(in) :: istwf_k, npw_k, nspinor
 real(dp),intent(in) :: eig0nk
 character(len=*),optional,intent(in) :: mode
!arrays
 real(dp),intent(inout) :: cwave(2,npw_k*nspinor)
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)
 real(dp) :: vk(3)

!Local variables-------------------------------
 character(len=50) :: my_mode
!arrays
 real(dp) :: cvk(2, 3)

!************************************************************************

 my_mode = "cart"; if (present(mode)) my_mode = mode
 call self%apply(eig0nk, npw_k, nspinor, cwave, cwaveprj)
 cvk = self%get_braket(eig0nk, istwf_k, npw_k, nspinor, cwave, mode=my_mode)
 vk = cvk(1, :)

end function ddkop_get_vdiag
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_free
!! NAME
!!  ddkop_free
!!
!! FUNCTION
!!  Free memory
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkop_free(self)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: idir

!************************************************************************

 ABI_SFREE(self%gh1c)
 ABI_SFREE(self%gs1c)

 do idir=1,3
   call self%gs_hamkq(idir)%free()
   call self%htg(idir)%free()
   call self%rf_hamkq(idir)%free()
 end do

 self%mpi_enreg => null()

end subroutine ddkop_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ham_targets_free
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ham_targets_free(self)

!Arguments ------------------------------------
!scalars
 class(ham_targets_t),intent(inout) :: self

!************************************************************************

 ABI_SFREE(self%ffnlk)
 ABI_SFREE(self%ffnl1)
 ABI_SFREE(self%kpg_k)
 ABI_SFREE(self%kpg1_k)
 ABI_SFREE(self%dkinpw)
 ABI_SFREE(self%kinpw1)
 ABI_SFREE(self%ph3d)
 ABI_SFREE(self%ph3d1)

end subroutine ham_targets_free
!!***

!----------------------------------------------------------------------

end module m_ddk
!!***
