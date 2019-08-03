!{\src2tex{textfont=tt}}
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
!! Copyright (C) 2016-2019 ABINIT group (MJV, HM, MG)
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
 use m_krank
 use m_fstab
 use m_wfd
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

 use m_fstrings,      only : strcat, sjoin, itoa, endswith, ktoa
 use m_symtk,         only : matr3inv
 use m_io_tools,      only : iomode_from_fname
 use m_time,          only : cwtime, sec2str
 use defs_abitypes,   only : hdr_type, dataset_type, MPI_type
 use defs_datatypes,  only : ebands_t, pseudopotential_type
 use m_geometry,      only : mkradim
 use m_crystal,       only : crystal_t
 use m_vkbr,          only : vkbr_t, nc_ihr_comm, vkbr_init, vkbr_free
 use m_pawtab,        only : pawtab_type
 use m_wfk,           only : wfk_read_ebands, wfk_read_h1mat

 implicit none

 private
!!***

 integer,private,parameter :: DDK_NOMODE    = 0
 integer,private,parameter :: DDK_READMODE  = 1
 integer,private,parameter :: DDK_WRITEMODE = 2

!----------------------------------------------------------------------

!!****t* m_ddk/ddk_t
!! NAME
!!  ddk_t
!!
!! FUNCTION
!!  object containing ddk derivatives ([H,r] proportional to band velocities)
!!
!! SOURCE

 type,public :: ddk_t

  integer :: comm
  ! MPI communicator used for IO.

  !integer :: version
  ! File format version read from file.

  integer :: iomode=IO_MODE_FORTRAN
  ! Method used to access the DDK file:
  !   IO_MODE_FORTRAN for usual Fortran IO routines
  !   IO_MODE_MPI if MPI/IO routines.
  !   IO_MODE_ETSF netcdf files in etsf format

  integer :: rw_mode = DDK_NOMODE
   ! (Read|Write) mode

  integer :: nsppol
   ! Number of spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: nene
    ! Number of energy points we may need the fsavg at

  integer :: nkfs
    ! Number of k-points on the Fermi-surface (FS-BZ).

  integer :: maxnb
   ! Max number of bands on the FS.
   ! TODO: Maybe maxnbfs

  integer :: usepaw
   ! 1 if PAW calculation, 0 otherwise

  integer :: prtvol=0
   ! Verbosity level

  logical :: debug=.False.
   ! Debug flag

  character(len=fnlen) :: paths(3) = ABI_NOFILE
   ! File name

  real(dp) :: acell(3),rprim(3,3),gprim(3,3)
   ! TODO: Are these really needed?

  real(dp), allocatable :: velocity (:,:,:,:)
   ! (3,maxnb,nkfs,nsppol)
   ! velocity on the FS in cartesian coordinates.

  real(dp), allocatable :: velocity_fsavg (:,:,:)
   ! (3,nene,nsppol)
   ! velocity on the FS in cartesian coordinates.

  logical :: use_ncddk(3)
   ! True if we are readin DDK matrix elements from EVK.nc instead of WFK file

  type(crystal_t) :: cryst
   ! Crystal structure read from file

 end type ddk_t

 public :: ddk_init              ! Initialize the object.
 public :: ddk_read_fsvelocities ! Read FS velocities from file.
 public :: ddk_fs_average_veloc  ! find FS average of velocity squared
 public :: ddk_free              ! Close the file and release the memory allocated.
 public :: ddk_print             ! output values
 public :: ddk_red2car           ! Convert band velocities from cartesian to reduced coordinates
 public :: ddk_compute           ! Calculate ddk matrix elements. Save result on disk.
!!***

 type, private :: ham_targets_t
   real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:)
   real(dp),allocatable :: kpg_k(:,:),kpg1_k(:,:)
   real(dp),allocatable :: ph3d(:,:,:),ph3d1(:,:,:)
   real(dp),allocatable :: dkinpw(:),kinpw1(:)
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

  type(gs_hamiltonian_type) :: gs_hamkq(3)

  type(rf_hamiltonian_type) :: rf_hamkq(3)

  type(ham_targets_t), private :: htg(3)
  ! Store arrays that are targetted by the hamiltonians.

  real(dp), private, allocatable :: gh1c(:,:,:)
   !gh1c, (2, mpw*nspinor, 3))

  real(dp), private, allocatable :: gs1c(:,:,:)
   ! gs1c, (2, mpw*nspinor, 3*psps%usepaw))

 contains

   procedure :: setup_spin_kpoint => ddkop_setup_spin_kpoint
    ! Prepare application of dH/dk for given spin, k-point.

   procedure :: apply => ddkop_apply
    ! Apply dH/dk to input wavefunction.

   procedure :: get_velocity => ddkop_get_velocity
    ! Compute matrix element.

   procedure :: free => ddkop_free
    ! Free memory.

 end type ddkop_t

 public :: ddkop_new   ! Build object

 !interface ddkop_t
 !  procedure ddkop_new
 !end interface ddkop_t
!!***

CONTAINS

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_init
!! NAME
!!  ddk_init
!!
!! FUNCTION
!!  Initialize the object from file. This is a COLLECTIVE procedure that must be called
!!  by each process in the MPI communicator comm.
!!
!! INPUTS
!!   paths=3 Filenames (could be either DFPT WFK files of DKK.nc files.
!!   comm=MPI communicator.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_init(ddk, paths, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: paths(3)
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk

!Local variables-------------------------------
!scalars
 integer,parameter :: timrev2=2,fform2=2
 integer :: my_rank, restart, restartpaw, ii
!arrays
 integer :: fforms(3)
 type(hdr_type) :: hdrs(3)

!************************************************************************

 my_rank = xmpi_comm_rank(comm)
 ddk%paths = paths; ddk%comm = comm
 ddk%iomode = iomode_from_fname(paths(1))

 ! In this calls everything is broadcast properly to the whole comm
 do ii=1,3
   ddk%use_ncddk(ii) = endswith(paths(ii), "_EVK.nc")
   call hdr_read_from_fname(hdrs(ii), paths(ii), fforms(ii), comm)
   if (ddk%debug) call hdr_echo(hdrs(ii), fforms(ii), 4, unit=std_out)
   ! check that 2 headers are compatible
   if (ii > 1) call hdr_check(fform2, fform2, hdrs(ii-1), hdrs(ii), 'COLL', restart, restartpaw)
 end do

 ddk%nsppol = hdrs(1)%nsppol
 ddk%nspinor = hdrs(1)%nspinor
 ddk%usepaw = hdrs(1)%usepaw
 ABI_CHECK(ddk%usepaw == 0, "PAW not yet supported")

 ! Init crystal_t
 ddk%cryst = hdr_get_crystal(hdrs(1), timrev2)

 ! Compute rprim, and gprimd. Used for slow FFT q--r if multiple shifts
 call mkradim(ddk%acell,ddk%rprim,ddk%cryst%rprimd)
 call matr3inv(ddk%rprim,ddk%gprim)

 do ii=1,3
   call hdr_free(hdrs(ii))
 end do

end subroutine ddk_init
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_compute
!! NAME
!!  ddk_compute
!!
!! FUNCTION
!!  Calculate the DDK matrix elements using the commutator formulation.
!!
!! INPUTS
!!  prefix: Prefix for output EVK file.
!!
!! PARENTS
!!      wfk_analyse
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddk_compute(wfk_path, prefix, dtset, psps, pawtab, ngfftc, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk_path, prefix
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
!arrays
 integer,intent(in) :: ngfftc(18)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: mband, nbcalc, nsppol, ib_v, ib_c
 integer :: mpw, spin, nspinor, nkpt, nband_k, npw_k
 integer :: ii, jj, ik, bandmin, bandmax, istwf_k, idir, edos_intmeth
 integer :: my_rank, nproc, ierr, bstop
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
 real(dp) :: edos_step, edos_broad
#ifdef HAVE_NETCDF
 integer :: ncerr, ncid
#endif
 character(len=500) :: msg
 character(len=fnlen) :: fname
 logical :: only_diago, is_kmesh
 type(wfd_t) :: wfd
 type(vkbr_t) :: vkbr
 type(ebands_t) :: ebands, ebands_tmp
 type(edos_t)  :: edos
 !type(jdos_t)  :: jdos
 type(gaps_t)  :: gaps
 type(crystal_t) :: cryst
 type(hdr_type) :: hdr_tmp, hdr
 type(ddkop_t) :: ddkop
!arrays
 integer,parameter :: voigt2ij(2, 6) = reshape([1, 1, 2, 2, 3, 3, 2, 3, 1, 3, 1, 2], [2, 6])
 integer,allocatable :: distrib_mat(:,:,:,:), distrib_diago(:,:,:),nband(:,:), kg_k(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
 real(dp) :: kpt(3), vv(2, 3), vred(2,3), vcar(2,3)
 real(dp) :: dummy_vals(1,1,1,1),dummy_vecs(1,1,1,1,1)
 real(dp) :: eminmax_spin(2,2)
 real(dp) :: emin, emax
 real(dp),allocatable :: dipoles(:,:,:,:,:,:)
 real(dp),allocatable :: dummy_dosvals(:,:,:,:), dummy_dosvecs(:,:,:,:,:), vvdos_tens(:,:,:,:,:,:)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:)
 real(dp),allocatable :: vdiago(:,:,:,:),vmat(:,:,:,:,:,:)
 real(dp),allocatable :: cg_c(:,:), cg_v(:,:)
 real(dp),allocatable :: vvdos_mesh(:) !, vvdos_vals(:,:,:,:)
 complex(dpc) :: vg(3), vr(3)
 complex(gwpc),allocatable :: ihrc(:,:), ug_c(:), ug_v(:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 if (my_rank == master) then
   write(msg, '(2a)') "Computation of velocity matrix elements (ddk)", ch10
   call wrtout(ab_out, msg); call wrtout(std_out, msg)
 end if

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
 end if

#ifndef HAVE_NETCDF
  MSG_ERROR("The matrix elements are only written if NETCDF is activated")
#endif

 ! Get ebands and hdr from WFK file.
 ebands = wfk_read_ebands(wfk_path, comm, out_hdr=hdr)
 cryst = hdr_get_crystal(hdr, 2)

 ! Extract important dimensions from hdr%
 nkpt    = hdr%nkpt
 nsppol  = hdr%nsppol
 nspinor = hdr%nspinor
 mband   = hdr%mband

 only_diago = .False.; if (dtset%wfk_task == WFK_TASK_DDK_DIAGO) only_diago = .True.

 ! Define band range.
 !TODO: hardcoded for now but should be an argument.
 bandmin = 1; bandmax = mband
 nbcalc  = bandmax - bandmin + 1

 if (my_rank == master) then
   call wrtout(ab_out, "Parameters extracted from Abinit header:")
   write(ab_out, "(a, f5.1)") 'ecut:    ', hdr%ecut
   write(ab_out, "(a, i0)")   'nkpt:    ', nkpt
   write(ab_out, "(a, i0)")   'mband:   ', mband
   write(ab_out, "(a, i0)")   'nsppol:  ', nsppol
   write(ab_out, "(a, i0)")   'nspinor: ', nspinor
   write(ab_out, "(a, i0)")   'useylm:  ', dtset%useylm
   write(ab_out, "(a, i0)")   'inclvkb: ', dtset%inclvkb
   write(ab_out, "(a)")""
 end if

 ! Create distribution of the wavefunctions mask
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkpt, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkpt, nsppol))

 keep_ur = .false.; bks_mask = .false.
 nband = mband

 if (only_diago) then
   ! Distribute k-points, spin and (b, b) diago over the processors.
   ABI_MALLOC(distrib_diago, (bandmin:bandmax,nkpt,nsppol))
   distrib_diago = -1

   ! Create bks_mask to load the wavefunctions
   ii = 0
   do spin=1,nsppol
     do ik=1,nkpt
       do ib_v=bandmin,bandmax
          ii = ii + 1; if (mod(ii, nproc) /= my_rank) cycle ! MPI parallelism.
          distrib_diago(ib_v, ik, spin) = my_rank
          bks_mask(ib_v,ik,spin) = .true.
       end do
     end do
   end do
   call wrtout(std_out, sjoin(" Rank: ", itoa(my_rank), "will treat", itoa(count(distrib_diago == my_rank))))

 else
   ! Distribute k-points, spin and (b, b') pairs over the processors
   ABI_MALLOC(distrib_mat, (bandmin:bandmax,bandmin:bandmax,nkpt,nsppol))
   call xmpi_distab(nproc, distrib_mat)

   ! Create bks_mask to load the wavefunctions
   do spin=1,nsppol
     do ik=1,nkpt
       ! Loop over v bands
       do ib_v=bandmin,bandmax
        ! Loop over c bands
         do ib_c=bandmin,bandmax
           if (distrib_mat(ib_c,ib_v,ik,spin) == my_rank) then
             bks_mask(ib_v,ik,spin) = .true.
             bks_mask(ib_c,ik,spin) = .true.
           end if
         end do
       end do
     end do
   end do

   call wrtout(std_out, sjoin(" Rank: ", itoa(my_rank), "will treat", itoa(count(distrib_mat == my_rank))))
 end if

 ! Initialize distributed wavefunctions object
 call wfd_init(wfd,cryst,pawtab,psps,keep_ur,mband,nband,nkpt,nsppol,&
   bks_mask,dtset%nspden,nspinor,hdr%ecut,dtset%ecutsm,dtset%dilatmx,ebands%istwfk,ebands%kptns,&
   ngfftc,dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

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
   ABI_MALLOC(cg_c, (2,mpw*nspinor))
   ABI_MALLOC(cg_v, (2,mpw*nspinor))
 end if

 ABI_DT_MALLOC(cwaveprj, (0, 0))

 ABI_CALLOC(dipoles, (3,2,mband,mband,nkpt,nsppol))
 ABI_MALLOC(ihrc,    (3, nspinor**2))

 if (only_diago) then
   ABI_CALLOC(vdiago, (3,mband,nkpt,nsppol))
 else
   ABI_CALLOC(vmat, (2,3,mband,mband,nkpt,nsppol))
 end if

 if (dtset%useria /= 666) then
   ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
   !if (my_rank == master) call ddkop%print(ab_out)
 end if

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 ! Loop over spins
 do spin=1,nsppol
   ! Loop over kpoints
   do ik=1,nkpt
     ! Only do a subset a k-points
     if (only_diago) then
       if (all(distrib_diago(:,ik,spin) /= my_rank)) cycle
     else
       if (all(distrib_mat(bandmin:bandmax,bandmin:bandmax,ik,spin) /= my_rank)) cycle
     end if
     call cwtime(cpu,wall,gflops,"start")

     nband_k  = wfd%nband(ik,spin)
     istwf_k  = wfd%istwfk(ik)
     npw_k    = wfd%npwarr(ik)
     kpt      = wfd%kibz(:,ik)
     kg_k(:,1:npw_k) = wfd%kdata(ik)%kg_k

     if (dtset%useria /= 666) then
       call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kpt, istwf_k, npw_k, kg_k)
     else
       ! Allocate KB form factors
       ! Prepare term i <n,k|[Vnl,r]|n"k>
       if (dtset%inclvkb/=0) call vkbr_init(vkbr,cryst,psps,dtset%inclvkb,istwf_k,npw_k,kpt,kg_k)
     end if

     ! Loop over bands
     do ib_v=bandmin,bandmax
       if (only_diago) then
         if (distrib_diago(ib_v,ik,spin) /= my_rank) cycle
       else
         if (all(distrib_mat(:,ib_v,ik,spin) /= my_rank)) cycle
       end if

       if (dtset%useria /= 666) then
         call wfd%copy_cg(ib_v, ik, spin, cg_v)
         call ddkop%apply(ebands%eig(ib_v, ik, spin), npw_k, wfd%nspinor, cg_v, cwaveprj, wfd%mpi_enreg)
       else
         ug_v(1:npw_k*nspinor) = wfd%wave(ib_v,ik,spin)%ug
       end if

       ! Loop over bands
       bstop = bandmax; if (only_diago) bstop = ib_v
       do ib_c=ib_v,bstop
         if (.not. only_diago) then
           if (distrib_mat(ib_c,ib_v,ik,spin) /= my_rank) cycle
         end if

         if (dtset%useria /= 666) then
           call wfd%copy_cg(ib_c, ik, spin, cg_c)
           vv = ddkop%get_velocity(ebands%eig(ib_c, ik, spin), istwf_k, npw_k, nspinor, wfd%mpi_enreg%me_g0, cg_c)
           !if (ib_v == ib_c) vv(2, :) = zero

           if (only_diago) then
             vdiago(:,ib_c,ik,spin) = vv(1, :)
           else
             vmat(:,:,ib_c,ib_v,ik,spin) = vv
             ! Hermitian conjugate
             if (ib_v /= ib_c) then
               vmat(1,:,ib_v,ib_c,ik,spin) = vv(1, :)
               vmat(2,:,ib_v,ib_c,ik,spin) = -vv(2, :)
             end if
           end if

           do idir=1,3
             dipoles(idir,:,ib_c,ib_v,ik,spin) = vv(:, idir)
             ! Hermitian conjugate
             if (ib_v /= ib_c) dipoles(idir,:,ib_v,ib_c,ik,spin) = [vv(1, idir), -vv(2, idir)]
           end do

         else
           ug_c(1:npw_k*nspinor) = wfd%wave(ib_c,ik,spin)%ug

           ! Calculate matrix elements of i[H,r] for NC pseudopotentials.
           ihrc = nc_ihr_comm(vkbr,cryst,psps,npw_k,nspinor,istwf_k,dtset%inclvkb, kpt,ug_c,ug_v,kg_k)

           ! HM: 24/07/2018
           ! Transform dipoles to be consistent with results from DFPT
           ! Perturbations with DFPT are along the reciprocal lattice vectors
           ! Perturbations with Commutator are along real space lattice vectors
           ! dot(A, DFPT) = X
           ! dot(B, COMM) = X
           ! B = 2 pi (A^{-1})^T =>
           ! dot(B^T B,COMM) = 2 pi DFPT
           vr = (2*pi)*(2*pi)*sum(ihrc(:,:),dim=2)
           vg(1) = dot_product(Cryst%gmet(1,:),vr)
           vg(2) = dot_product(Cryst%gmet(2,:),vr)
           vg(3) = dot_product(Cryst%gmet(3,:),vr)

           ! Save matrix elements of i*r in the IBZ
           dipoles(:,1,ib_c,ib_v,ik,spin) = real(vg, kind=dp)
           dipoles(:,1,ib_v,ib_c,ik,spin) = real(vg, kind=dp) ! Hermitian conjugate
           if (ib_v == ib_c) then
              dipoles(:,2,ib_c,ib_v,ik,spin) = zero
              dipoles(:,2,ib_v,ib_c,ik,spin) = zero
           else
              dipoles(:,2,ib_c,ib_v,ik,spin) =  aimag(vg)
              dipoles(:,2,ib_v,ib_c,ik,spin) = -aimag(vg) ! Hermitian conjugate
           end if
         end if

       end do
     end do

     ! Free KB form factors
     call vkbr_free(vkbr)

     call cwtime(cpu,wall,gflops,"stop")
     write(msg,'(2(a,i0),2(a,f8.2))')"k-point [",ik,"/",nkpt,"] completed. cpu:",cpu,", wall:",wall
     call wrtout(std_out, msg, do_flush=.True.)

   end do ! k-points
 end do ! spin

 call cwtime(cpu_all, wall_all, gflops_all, "stop")
 call wrtout(std_out, sjoin("Calculation completed. cpu-time:", sec2str(cpu_all), ",wall-time:", &
   sec2str(wall_all)), do_flush=.True.)

 ABI_FREE(ug_c)
 ABI_FREE(ug_v)
 ABI_FREE(kg_k)
 ABI_FREE(ihrc)
 ABI_DT_FREE(cwaveprj)

 ABI_SFREE(distrib_mat)
 ABI_SFREE(distrib_diago)

 if (dtset%useria /= 666) then
   ABI_FREE(cg_c)
   ABI_FREE(cg_v)
   call ddkop%free()
 end if

 ! Gather the k-points computed by all processes
 call xmpi_sum_master(dipoles, master, comm, ierr)

 if (only_diago) then
   call xmpi_sum_master(vdiago, master, comm, ierr)
 else
   call xmpi_sum_master(vmat, master, comm, ierr)
 end if

 is_kmesh = hdr%kptopt > 0

 if (is_kmesh .and. only_diago) then
   ! Compute electron DOS with tetra.
   edos_intmeth = 2
   if (dtset%prtdos == 1) edos_intmeth = 1
   if (dtset%prtdos == -2) edos_intmeth = 3
   edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
   if (edos_step == 0) edos_step = 0.001
   !edos = ebands_get_edos(ebands, cryst, edos_intmeth, edos_step, edos_broad, comm)
   !jdos = ebands_get_jdos(ebands, cryst, edos_intmeth, edos_step, edos_broad, comm, ierr)

   ! Compute (v x v) DOS. Upper triangle in Voigt format.
   ABI_MALLOC(vv_tens, (3, 3, 1, mband, nkpt, nsppol))
   do spin=1,nsppol
     do ik=1,nkpt
       do ib_v=bandmin,bandmax
         !vr = vdiago(:,ib_v,ik,spin)
         ! Go to cartesian coordinates (same as pmat2cart routine).
         vred(1,:) = vdiago(:,ib_v,ik,spin)
         call ddk_red2car(cryst%rprimd,vred,vcar)
         vr = vcar(1,:)
         do ii=1,3
           do jj=1,3
             vv_tens(ii, jj, 1, ib_v, ik, spin) = vr(ii) * vr(jj)
           end do
         end do
       end do
     end do
   end do

   !set default erange
   eminmax_spin(:,:ebands%nsppol) = get_minmax(ebands, "eig")
   emin = minval(eminmax_spin(1,:)); emin = emin - 0.1_dp * abs(emin)
   emax = maxval(eminmax_spin(2,:)); emax = emax + 0.1_dp * abs(emax)

   ! If sigma_erange is set, get emin and emax
   ierr = get_gaps(ebands,gaps)
   if (ierr/=0.and.ebands%occopt.eq.1) then
     call ebands_copy(ebands,ebands_tmp)
     call ebands_set_scheme(ebands_tmp, ebands%occopt, ebands%tsmear, dtset%spinmagntarget, dtset%prtvol)
     call gaps%free()
     ierr = get_gaps(ebands_tmp,gaps)
     call ebands_free(ebands_tmp)
   end if
   do spin=1,ebands%nsppol
     if (dtset%sigma_erange(1) >= zero) emin = gaps%vb_max(spin) + tol2 * eV_Ha - dtset%sigma_erange(1)
     if (dtset%sigma_erange(2) >= zero) emax = gaps%cb_min(spin) - tol2 * eV_Ha + dtset%sigma_erange(2)
   end do
   call gaps%free()

   edos = ebands_get_dos_matrix_elements(ebands, cryst, &
                                         dummy_vals, 0, dummy_vecs, 0, vv_tens, 1, &
                                         edos_intmeth, edos_step, edos_broad, comm, vvdos_mesh, &
                                         dummy_dosvals, dummy_dosvecs, vvdos_tens, &
                                         emin=emin, emax=emax)
   ABI_SFREE(dummy_dosvals)
   ABI_SFREE(dummy_dosvecs)
   ABI_SFREE(vv_tens)
 end if

 ! Write the matrix elements
#ifdef HAVE_NETCDF
 ! Output EVK file in netcdf format.
 if (my_rank == master) then
   ! Have to build hdr on k-grid with info about perturbation.
   call hdr_copy(hdr, hdr_tmp)
   hdr_tmp%qptn = zero

   fname = strcat(prefix, "NEW_EVK.nc")
   call wrtout(ab_out, sjoin("- Writing file: ", fname))
   NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EVK.nc file")
   hdr_tmp%pertcase = 0
   NCF_CHECK(hdr_ncwrite(hdr_tmp, ncid, 43, nc_define=.True.))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(ebands, ncid))
   if (only_diago) then
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t('vred_diagonal', "dp", "three, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
   else
     ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vred_matrix', "dp", &
         "two, three, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
   end if
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_set_datamode(ncid))
   if (only_diago) then
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vred_diagonal"), vdiago))
   else
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vred_matrix"), vmat))
   end if
   if (is_kmesh .and. only_diago) then
     NCF_CHECK(edos%ncwrite(ncid))
     !NCF_CHECK(jdos%ncwrite(ncid))
     ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vvdos_mesh', "dp", "edos_nw")], defmode=.True.)
     ncerr = nctk_def_arrays(ncid, [ nctkarr_t('vvdos_vals', "dp", "edos_nw, nsppol_plus1, three, three")], defmode=.True.)
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_mesh"), vvdos_mesh))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvdos_vals"), vvdos_tens(:,1,:,:,:,1)))
   end if
   NCF_CHECK(nf90_close(ncid))

   do ii=1,3
     fname = strcat(prefix, '_', itoa(ii), "_EVK.nc")
     !call wrtout(ab_out, sjoin("- Writing file: ", fname))
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EVK.nc file")
     hdr_tmp%pertcase = 3 * cryst%natom + ii
     NCF_CHECK(hdr_ncwrite(hdr_tmp, ncid, 43, nc_define=.True.))
     NCF_CHECK(cryst%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(ebands, ncid))
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t('h1_matrix_elements', "dp", &
        "two, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ncid))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "h1_matrix_elements"), dipoles(ii,:,:,:,:,:)))
     NCF_CHECK(nf90_close(ncid))
   end do
   call hdr_free(hdr_tmp)
 end if
#endif

 if (my_rank == master .and. dtset%prtvol > 0) then
   write(ab_out, "(2a)")ch10,"Writing velocity matrix elements (only diagonal terms, real part) for testing purpose:"
   do spin=1,nsppol
     do ik=1,min(nkpt, 4)
       write(ab_out, "(2(a,1x,i0),2x,2a)")"For spin: ", spin, ", ikbz: ", ik, ", kpt: ", trim(ktoa(wfd%kibz(:,ik)))
       do ib_c=bandmin,min(bandmin+8, bandmax)
         write(ab_out, "(3(es16.6,2x))") dipoles(:,1,ib_c,ib_c,ik,spin)
       end do
       write(ab_out,*)""
       !do ib_c=bandmin,min(bandmin+8, bandmax)
       !  write(ab_out, "(a, 6(es16.6,2x))")"Sum_k: ", sum(dipoles(:,:,ib_c,ib_c,:,spin), dim=3) / nkpt
       !end do
     end do
   end do
 end if

 ! Free memory
 ABI_SFREE(dipoles)
 ABI_SFREE(vdiago)
 ABI_SFREE(vmat)

 ABI_SFREE(vvdos_mesh)
 ABI_SFREE(vvdos_tens)
 call edos%free()
 !call jdos%free()
 call wfd%free()
 call ebands_free(ebands)
 call cryst%free()
 call hdr_free(hdr)

 ! Block all procs here so that we know output files are available when code returns.
 call xmpi_barrier(comm)

end subroutine ddk_compute
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_read_fsvelocities
!! NAME
!!  ddk_read_fsvelocities
!!
!! FUNCTION
!!  Read FS velocities from DDK files. Returned in ddk%velocity
!!
!! INPUTS
!!   fstab(ddk%nsppol)=Tables with the correspondence between points of the Fermi surface (FS)
!!     and the k-points in the IBZ
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_read_fsvelocities(ddk, fstab, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(ddk_t),intent(inout) :: ddk
 type(fstab_t),target,intent(in) :: fstab(ddk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: idir, ikfs, isppol, ik_ibz, ii
 integer :: symrankkpt, ikpt_ddk,iband, bd2tot_index
 integer :: bstart_k, nband_k, nband_in, vdim
#ifdef HAVE_NETCDF
 integer :: ncid, varid, nc_fform
#endif
 type(hdr_type) :: hdr1
 type(krank_t) :: krank
 type(fstab_t), pointer :: fs
 character(len=500) :: msg
!arrays
 real(dp), allocatable :: eigen1(:), velocityp(:,:)

!************************************************************************

 if (ddk%rw_mode /= DDK_NOMODE) then
   MSG_ERROR("ddk should be in ddk_NOMODE before open_read is called.")
 end if
 ddk%rw_mode = DDK_READMODE

 ddk%maxnb = maxval(fstab(:)%maxnb)
 ddk%nkfs = maxval(fstab(:)%nkfs)
 ABI_MALLOC(ddk%velocity, (3,ddk%maxnb,ddk%nkfs,ddk%nsppol))

 call wrtout(std_out, sjoin('Read DDK FILES with iomode=', itoa(ddk%iomode)), 'COLL')
 do idir = 1,3
   ! Open the files. All procs in comm receive hdr1 and eigen1
   if (ddk%use_ncddk(idir)) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(ncid, ddk%paths(ii), comm))
     call hdr_ncread(hdr1, ncid, nc_fform)
     varid = nctk_idname(ncid, "h1_matrix_elements")
     ABI_MALLOC(eigen1, (2*hdr1%mband*hdr1%mband*hdr1%nkpt*ddk%nsppol))
     !NCF_CHECK(nctk_set_collective(ncid, varid))
     NCF_CHECK(nf90_get_var(ncid, varid, eigen1, count=[2, hdr1%mband, hdr1%mband, hdr1%nkpt, hdr1%nsppol]))
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Netcdf not available!")
#endif
   else
     call wfk_read_h1mat(ddk%paths(idir), eigen1, hdr1, comm)
   end if
   nband_in = maxval(hdr1%nband)

   ! need correspondence hash between the DDK and the fs k-points
   krank = krank_new(hdr1%nkpt, hdr1%kptns)
   do isppol=1,ddk%nsppol
     fs => fstab(isppol)
     do ikfs=1,fs%nkfs
       ik_ibz = fs%istg0(1,ikfs)
       symrankkpt = krank%get_rank (fs%kpts(:,ikfs))
       ikpt_ddk = krank%invrank(symrankkpt)
       if (ikpt_ddk == -1) then
         write(msg, "(3a)")&
           "Error in correspondence between ddk and fsk kpoint sets",ch10,&
           "kpt sets in fsk and ddk files must agree."
         MSG_ERROR(msg)
       end if

       bd2tot_index=2*nband_in**2*(ikpt_ddk-1) + 2*nband_in**2*hdr1%nkpt*(isppol-1)
       bstart_k = fs%bstcnt_ibz(1, ik_ibz)
       nband_k = fs%bstcnt_ibz(2, ik_ibz)
       ! first derivative eigenvalues for k-point. Diagonal of eigen1 is real -> only use that part
       do iband = bstart_k, bstart_k+nband_k-1
         ddk%velocity(idir, iband-bstart_k+1, ikfs, isppol)=eigen1(bd2tot_index + 2*nband_in*(iband-1) + iband)
       end do
     end do
   end do

   ABI_FREE(eigen1)
   call krank%free()
   call hdr_free(hdr1)
 end do ! idir

 ! process the eigenvalues(1): rotate to cartesian and divide by 2 pi
 ! use DGEMM better here on whole matrix and then reshape?
 vdim = ddk%maxnb*ddk%nkfs*ddk%nsppol
 ABI_MALLOC(velocityp, (3,vdim))
 velocityp = zero
 call dgemm('n','n',3,vdim,3,one,ddk%cryst%rprimd,3,ddk%velocity,3,zero,velocityp,3)

! do isppol = 1, ddk%nsppol
!   do ikpt = 1, ddk%nkfs
!     do iband = 1, ddk%maxnb
!       tmpveloc = ddk%velocity (:, iband, ikpt, isppol)
!       call dgemm('n','n',3,3,3,one,ddk%cryst%rprimd,3,tmpveloc,3,zero,tmpveloc2,3)
!       ddk%velocity (:, iband, ikpt, isppol) = tmpveloc2
!     end do
!   end do
! end do

 ddk%velocity = reshape (velocityp, [3,ddk%maxnb,ddk%nkfs,ddk%nsppol])

 ABI_FREE(velocityp)
 call krank%free()

end subroutine ddk_read_fsvelocities
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_fs_average_veloc
!! NAME
!!  ddk_fs_average_veloc
!!
!! FUNCTION
!!  Perform Fermi surface average of velocity squared then square rooted, print and store in ddk object
!!
!! INPUTS
!!   ddk = object with electron band velocities
!!   fstab(ddk%nsppol)=Tables with the correspondence between points of the Fermi surface (FS)
!!     and the k-points in the IBZ
!!   comm=MPI communicator
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_fs_average_veloc(ddk, ebands, fstab, sigmas)

!Arguments ------------------------------------
!scalars
!integer,intent(in) :: comm  ! could distribute this over k in the future
 real(dp),intent(in) :: sigmas(:)
 type(ebands_t),intent(in) :: ebands
 type(ddk_t),intent(inout) :: ddk
 type(fstab_t),target,intent(in) :: fstab(ddk%nsppol)

!Local variables-------------------------------
!scalars
 integer :: idir, ikfs, isppol, ik_ibz, iene
 integer :: iband, mnb, nband_k, nsig
 type(fstab_t), pointer :: fs
!arrays
 real(dp), allocatable :: wtk(:,:)

!************************************************************************

 ddk%nene = fstab(1)%nene
 ABI_MALLOC(ddk%velocity_fsavg, (3,ddk%nene,ddk%nsppol))
 ddk%velocity_fsavg = zero

 nsig = size(sigmas, dim=1)
 mnb = 1
 do isppol=1,ddk%nsppol
   fs => fstab(isppol)
   mnb = max(mnb, maxval(fs%bstcnt_ibz(2, :)))
 end do
 ABI_MALLOC(wtk, (nsig,mnb))

 do isppol=1,ddk%nsppol
   fs => fstab(isppol)
   do iene = 1, fs%nene
     do ikfs=1,fs%nkfs
       ik_ibz = fs%istg0(1,ikfs)
       nband_k = fs%bstcnt_ibz(2, ik_ibz)
       call fstab_weights_ibz(fs, ebands, ik_ibz, isppol, sigmas, wtk, iene)

       do idir = 1,3
         do iband = 1, nband_k
           ddk%velocity_fsavg(idir, iene, isppol) = ddk%velocity_fsavg(idir, iene, isppol) + &
&             wtk(1,iband) * ddk%velocity(idir, iband, ikfs, isppol)**2
!&             fs%tetra_wtk_ene(iband,ik_ibz,iene) * ddk%velocity(idir, iband, ikfs, isppol)**2
         end do
       end do ! idir
     end do ! ikfs
   end do ! iene
  ! sqrt is element wise on purpose
   ddk%velocity_fsavg(:,:,isppol) = sqrt(ddk%velocity_fsavg(:,:,isppol)) / dble(fs%nkfs)
 end do ! isppol

 ABI_DEALLOCATE(wtk)

end subroutine ddk_fs_average_veloc
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_free
!! NAME
!!  ddk_free
!!
!! FUNCTION
!! Close the file and release the memory allocated.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_free(ddk)

!Arguments ------------------------------------
!scalars
 type(ddk_t),intent(inout) :: ddk

!************************************************************************

 ! integer arrays

 ! real arrays
 ABI_SFREE(ddk%velocity)
 ABI_SFREE(ddk%velocity_fsavg)

 ! types
 call ddk%cryst%free()

end subroutine ddk_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddk_red2car
!! NAME
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

subroutine ddk_red2car(rprimd, vred, vcar)

!Arguments -------------------------------------
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: vred(2,3)
 real(dp),intent(out) :: vcar(2,3)

!Local variables -------------------------------
 real(dp) :: vtmp(2,3)

 ! Go to cartesian coordinates (same as pmat2cart routine)
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

!!****f* m_ddk/ddk_print
!! NAME
!!  ddk_print
!!
!! FUNCTION
!!  Print info on the object.
!!
!! INPUTS
!! [unit]=the unit number for output
!! [prtvol]=verbosity level
!! [mode_paral]=either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only printing.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddk_print(ddk, header, unit, prtvol, mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(ddk_t),intent(in) :: ddk

!Local variables-------------------------------
!scalars
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt =std_out; if (PRESENT(unit)) my_unt   =unit
 my_prtvol=0    ; if (PRESENT(prtvol)) my_prtvol=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the ddk% object ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(std_out,*)"Number of FS bands: ",ddk%maxnb
 write(std_out,*)"Number of FS k-points: ",ddk%nkfs
 write(std_out,*)"Number of spin channels: ",ddk%nsppol
 write(std_out,*)"Paths to files: "
 write(std_out,*)"  ", trim(ddk%paths(1))
 write(std_out,*)"  ", trim(ddk%paths(2))
 write(std_out,*)"  ", trim(ddk%paths(3))

end subroutine ddk_print
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
 type(MPI_type),intent(in) :: mpi_enreg
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
 ABI_UNUSED((/mpi_enreg%paral_kgb/))

 new%inclvkb = dtset%inclvkb
 new%usepaw = dtset%usepaw
 new%ipert = cryst%natom + 1
 new%dfpt_sciss = dtset%dfpt_sciss
 new%mpw = mpw

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
!!
!! FUNCTION
!!
!! INPUTS
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
 ABI_MALLOC(ylmgr1_k, (npw_k,3+6*(optder/2),psps%mpsang**2*psps%useylm*useylmgr1))

 if (psps%useylm == 1) then ! .and. self%inclvkb /= 0
   ! Fake MPI_type for sequential part. dummy_nband and nsppol1 are not used in sequential mode.
   call initmpi_seq(mpienreg_seq)
   dummy_nband = 0
   npwarr = npw_k
   call initylmg(cryst%gprimd,kg_k,kpoint,nkpt1,mpienreg_seq,psps%mpsang,npw_k,dummy_nband,nkpt1,&
      npwarr,nsppol1,optder,cryst%rprimd,ylm_k,ylmgr1_k)
   call destroy_mpi_enreg(mpienreg_seq)
 end if

 do idir=1,3
   call self%htg(idir)%free()

   ! Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(self%gs_hamkq(idir), spin, with_nonlocal=.true.)
   call load_spin_rf_hamiltonian(self%rf_hamkq(idir), spin, with_nonlocal=.true.)

   !if (self%inclvkb /= 0) then

   ! We need ffnl1 and dkinpw for 3 dirs. Note that the Hamiltonian objects use pointers to keep a reference
   ! to the output results of this routine.
   ! This is the reason why we need to store the targets in self%htg
   call getgh1c_setup(self%gs_hamkq(idir),self%rf_hamkq(idir),dtset,psps,kpoint,kpoint,idir,self%ipert, & ! In
     cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,npw_k,npw_k, &             ! In
     useylmgr1,kg_k,ylm_k,kg_k,ylm_k,ylmgr1_k, &                                       ! In
     self%htg(idir)%dkinpw,nkpg,nkpg1,self%htg(idir)%kpg_k,self%htg(idir)%kpg1_k, &    ! Out
     self%htg(idir)%kinpw1,self%htg(idir)%ffnlk,self%htg(idir)%ffnl1, &                ! Out
     self%htg(idir)%ph3d, self%htg(idir)%ph3d1)                                        ! Out
 end do

 ABI_FREE(ylm_k)
 ABI_FREE(ylmgr1_k)

end subroutine ddkop_setup_spin_kpoint
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_apply
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!  cwave(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!     if not allocated or size=0, they are locally computed (and not sorted)!!
!!
!! OUTPUT
!! gh1c(2,npw1*nspinor)= <G|H^(1)|C> or  <G|H^(1)-lambda.S^(1)|C> on the k+q sphere
!!                     (only kinetic+non-local parts if optlocal=0)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkop_apply(self, eig0nk, npw_k, nspinor, cwave, cwaveprj, mpi_enreg)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(inout) :: self
 integer,intent(in) :: npw_k, nspinor
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(in) :: eig0nk

!Local variables-------------------------------
!scalars
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
     call getgh1c(berryopt0,cwave,cwaveprj,self%gh1c(:,:,idir),&
       grad_berry,self%gs1c(:,:,idir),self%gs_hamkq(idir),gvnlx1,idir,self%ipert,eshift,mpi_enreg,optlocal0, &
       optnl,opt_gvnlx1,self%rf_hamkq(idir),sij_opt,tim_getgh1c,usevnl0)
   end do

 else
   ! optnl 0 with DDK does not work as expected. So I treat the kinetic term explicitly
   ! without calling getgh1c.
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

!!****f* m_ddk/ddkop_get_velocity
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

function ddkop_get_velocity(self, eig0mk, istwf_k, npw_k, nspinor, me_g0, brag) result(vv)

!Arguments ------------------------------------
!scalars
 class(ddkop_t),intent(in) :: self
 integer,intent(in) :: istwf_k, npw_k, nspinor, me_g0
 real(dp),intent(in) :: eig0mk
!arrays
 real(dp),intent(in) :: brag(npw_k*nspinor)
 real(dp) :: vv(2, 3)

!Local variables-------------------------------
!scalars
 integer :: idir
 real(dp) :: doti
!arrays
 real(dp) :: dotarr(2)

!************************************************************************

 if (self%usepaw == 0) then
   !<u_(iband,k+q)^(0)|H_(k+q,k)^(1)|u_(jband,k)^(0)>  (NC psps)
   do idir=1,3
     dotarr = cg_zdotc(npw_k * nspinor, brag, self%gh1c(:,:,idir))
     if (istwf_k > 1) then
       !dum = two * j_dpc * AIMAG(dum); if (vkbr%istwfk==2) dum = dum - j_dpc * AIMAG(gamma_term)
       doti = two * dotarr(2)
       if (istwf_k == 2 .and. me_g0 == 1) then
         ! nspinor always 1
         ! TODO: Recheck this part but it should be ok.
         doti = doti - (brag(1) * self%gh1c(2,1,idir) - brag(2) * self%gh1c(1,1,idir))
       end if
       dotarr(2) = doti; dotarr(1) = zero
     end if
     vv(:, idir) = dotarr
   end do
 else
   MSG_ERROR("PAW Not Implemented")
   ! <u_(iband,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(jband,k)^(0)> (PAW)
   ! eshiftkq = half * (eig0mk - self%eig0nk)
   ABI_UNUSED(eig0mk)
 end if

end function ddkop_get_velocity
!!***

!----------------------------------------------------------------------

!!****f* m_ddk/ddkop_free
!! NAME
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
   call destroy_hamiltonian(self%gs_hamkq(idir))
   call self%htg(idir)%free()
   call destroy_rf_hamiltonian(self%rf_hamkq(idir))
 end do

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
