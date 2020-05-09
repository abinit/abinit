!!****m* ABINIT/m_ephwg
!! NAME
!! m_ephwg
!!
!! FUNCTION
!!  Tools and objects to compute the weights used for the BZ integration of EPH quantities.
!!  More specifically the integration of quantities such as the imaginary part of the self-energy
!!  involving delta functions. Different approaches are available:
!!
!!    1.
!!    2.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG, HM)
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

module m_ephwg

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_copy
 use m_dtset
 use m_htetra
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_crystal
 use m_ifc
 use m_lgroup
 use m_ebands
 use m_eph_double_grid

 use defs_datatypes,    only : ebands_t
 use m_time,            only : cwtime, cwtime_report
 use m_symtk,           only : matr3inv
 use m_numeric_tools,   only : arth, inrange, wrap2_pmhalf
 use m_special_funcs,   only : gaussian
 use m_fstrings,        only : strcat, ltoa, itoa, ftoa, ktoa, sjoin
 use m_simtet,          only : sim0onei, SIM0TWOI
 use m_kpts,            only : kpts_timrev_from_kptopt, listkk, kpts_ibz_from_kptrlatt
 use m_occ,             only : occ_fd, occ_be

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_ephwg/ephwg_t
!! NAME
!! ephwg_t
!!
!! FUNCTION
!!  Stores electron eigevalues and phonon frequencies in the IBZ (assume same mesh for e and ph).
!!  Provides tools to compute (e_{k+q} - w{q}) in the IBZ(k)
!!  and integrate the delta functions for phonon emission/absorption with the tetrahedron method.
!!
!! SOURCE

type, public :: ephwg_t

  integer :: natom3
  ! 3 * natom

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: nbcount
  ! Number of bands treated.

  integer :: bstart
  ! The fist band (global index) starts at bstart.
  ! Used to select bands around the Fermi level.

  integer :: kptopt
  ! Option for k-point generation.

  integer :: timrev
  ! 1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nibz, nbz
  ! Number of q-points in IBZ and full BZ.

  integer :: nq_k
  ! Number of points in IBZ(k) i.e. the irreducible wedge
  ! defined by the operations of the little group of k.

  ! Integer controling whether to compute and store the electron-phonon matrix elements
  ! computed from generalized Frohlich model
  ! C. Verdi and F. Giustino, Phys. Rev. Lett. 115, 176401 (2015).
  ! 1,2: Unused (for compatibility with sigma_t)
  ! 3: Compute Frohlich electron-phonon matrix elements in the IBZ for later integration

  integer, allocatable :: kq2ibz(:)
  ! kq2ibz(nq_k)
  ! Mapping (k + q) --> initial IBZ array

  real(dp),allocatable :: ibz(:,:)
  ! ibz(3, nibz)
  ! The initial IBZ.

  real(dp),allocatable :: bz(:,:)
  ! bz(3, nbz)
  ! points in full BZ.

  real(dp) :: klatt(3, 3)
  ! Reciprocal of lattice vectors for full kpoint grid. Used by init_tetra

  integer,allocatable :: lgk2ibz(:)
  ! lgk2ibz(nq_k)
  ! Mapping Little-group IBZ_k --> initial IBZ
  ! TODO: This should be generalized to have the symmetry indices as well so
  ! that we can use it in sigmaph but then we have to implement similar algo for double grid.

  real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states
  ! (not necessarly equal to global nband, see also bstart and bcount)

  type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure (does not own memory)

  type(lgroup_t) :: lgk
  ! Little group of the k-point

  type(htetra_t) :: tetra_k
  ! Used to evaluate delta(w - e_{k+q} +/- phw_q) with tetrahedron method.

 contains

     procedure :: setup_kpoint => ephwg_setup_kpoint
     ! Prepare tetrahedron method for given external k-point.

     procedure :: double_grid_setup_kpoint => ephwg_double_grid_setup_kpoint
     ! Prepare tetrahedron method for given external k-point using double grid routines.

     procedure :: report_stats => ephwg_report_stats
     ! Report how much memory is being used by this object

     procedure :: get_deltas => ephwg_get_deltas
     ! Compute weights for $ \int \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $

     procedure :: get_deltas_wvals => ephwg_get_deltas_wvals
     ! Compute weights for $ \int \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $

     procedure :: get_deltas_qibzk => ephwg_get_deltas_qibzk
     ! Compute weights for $ \delta(\omega - \omega_{q\nu}} $ using IBZ(k) ordering

     procedure :: get_zinv_weights => ephwg_get_zinv_weights
     ! Compute weights for $ \int 1 / (\omega - \ee_{k+q, b} \pm \omega_{q\nu} $

     procedure :: free => ephwg_free
     ! Free memory
 end type ephwg_t

 public :: ephwg_new             ! Basic Constructor
 public :: ephwg_from_ebands     ! Build object from ebands_t

contains
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_new
!! NAME
!! ephwg_new
!!
!! FUNCTION
!!  Initialize the object from the electronic eigenvalues given in the IBZ.
!!
!! INPUTS
!!  cryst<cryst_t>=Crystalline structure.
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  bstart=Index of the first band to be included.
!!  nbcount=Number of bands included
!!  kptopt=Option for the k-point generation.
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk= number of shift vectors.
!!  shiftk(3,nshiftk)=shift vectors for k point generation
!!  nkibz=Number of points in the IBZ
!!  kibz(3,nkibz)=Reduced coordinates of the k-points in the IBZ.
!!  nsppol=Number of independent spin polarizations.
!!  eig_ibz(nbcount, nkibz, nsppol) = Electron eigenvalues for nbcount states
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ephwg_t) function ephwg_new( &
&  cryst, ifc, bstart, nbcount, kptopt, kptrlatt, nshiftk, shiftk, nkibz, kibz, nsppol, eig_ibz, comm) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt, nshiftk, nkibz, bstart, nbcount, nsppol, comm
 type(crystal_t),target,intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3, nshiftk), kibz(3, nkibz)
 real(dp),intent(in) :: eig_ibz(nbcount, nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: nprocs, my_rank, ik, ierr, out_nkibz
 real(dp) :: cpu, wall, gflops
!arrays
 real(dp) :: rlatt(3,3)
 integer :: out_kptrlatt(3,3)
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom), phfrq(3*cryst%natom)
 real(dp),allocatable :: out_kibz(:,:), out_wtk(:)

!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 new%natom3 = ifc%natom * 3
 new%nsppol = nsppol
 new%nbcount = nbcount
 new%bstart = bstart
 new%kptopt = kptopt
 new%timrev = kpts_timrev_from_kptopt(new%kptopt)
 new%nibz = nkibz
 new%cryst => cryst
 call alloc_copy(kibz, new%ibz)

 call cwtime(cpu, wall, gflops, "start")

 ! Get full BZ (new%nbz, new%bz) and new kptrlatt for tetra.
 call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, out_nkibz, out_kibz, out_wtk, new%nbz, new%bz, &
   new_kptrlatt=out_kptrlatt)
 call cwtime_report(" ephwg_new: kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 rlatt = out_kptrlatt; call matr3inv(rlatt, new%klatt)

 ABI_CHECK(size(out_kibz, dim=2) == new%nibz, "mismatch in nkibz!")
 ABI_FREE(out_kibz)
 ABI_FREE(out_wtk)

 ! Copy eigenvalues in IBZ. Change shape for better performance in other routines.
 ABI_MALLOC(new%eigkbs_ibz, (new%nibz, new%nbcount, new%nsppol))
 do ik=1,new%nibz
   new%eigkbs_ibz(ik, :, :) = eig_ibz(:, ik, :)
 end do

 ! Fourier interpolate phonon frequencies on the same mesh.
 ABI_CALLOC(new%phfrq_ibz, (new%nibz, new%natom3))

 do ik=1,new%nibz
   if (mod(ik, nprocs) /= my_rank) cycle ! mpi-parallelism
   call ifc%fourq(cryst, new%ibz(:, ik), phfrq, displ_cart)
   new%phfrq_ibz(ik, :) = phfrq
 end do

 call xmpi_sum(new%phfrq_ibz, comm, ierr)
 call cwtime_report(" ephwg_new: ifc_fourq", cpu, wall, gflops)

end function ephwg_new
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_from_ebands
!! NAME
!! ephwg_from_ebands
!!
!! FUNCTION
!!  Convenience constructor to initialize the object from an ebands_t object

type(ephwg_t) function ephwg_from_ebands(cryst, ifc, ebands, bstart, nbcount, comm) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::  bstart, nbcount, comm
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 real(dp),allocatable :: eig_ibz(:, :, :)

!----------------------------------------------------------------------

 if (bstart == 1 .and. nbcount == ebands%mband) then
   new = ephwg_new(cryst, ifc, bstart, nbcount, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, &
      ebands%kptns, ebands%nsppol, ebands%eig, comm)
 else
   ABI_CHECK(inrange(bstart, [1, ebands%mband]), "Wrong bstart")
   ABI_CHECK(inrange(bstart + nbcount - 1, [1, ebands%mband]), "Wrong nbcount")
   ! Copy submatrix of eigenvalues
   ABI_MALLOC(eig_ibz, (nbcount, ebands%nkpt, ebands%nsppol))
   eig_ibz = ebands%eig(bstart:bstart+nbcount-1, : , :)
   new = ephwg_new(cryst, ifc, bstart, nbcount, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, &
      ebands%kptns, ebands%nsppol, eig_ibz, comm)
   ABI_FREE(eig_ibz)
 end if

end function ephwg_from_ebands
!!***

!!****f* m_ephwg/ephwg_setup_kpoint
!! NAME
!! ephwg_setup_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute integration weights for a given k-point.
!!
!! INPUTS
!!  kpoint(3): k-point in reduced coordinates.
!!  prtvol: Verbosity level
!!  comm: MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_setup_kpoint(self, kpoint, prtvol, comm, skip_mapping )

!Arguments ------------------------------------
!scalars
 class(ephwg_t),target,intent(inout) :: self
 integer,intent(in) :: prtvol, comm
 logical,optional,intent(in) :: skip_mapping
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1
 integer :: ierr,ii
 logical :: do_mapping
 real(dp) :: dksqmax, cpu, wall, gflops
 character(len=80) :: errorstring
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: indkk(:,:)

!----------------------------------------------------------------------

 do_mapping = .true.; if (present(skip_mapping)) do_mapping = .not. skip_mapping
 cryst => self%cryst
 call cwtime(cpu, wall, gflops, "start")

 ! Get little group of the (external) kpoint.
 call self%lgk%free()
 self%lgk = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz, comm)
 if (prtvol > 0) call self%lgk%print()
 self%nq_k = self%lgk%nibz

 call cwtime_report(" lgroup_new", cpu, wall, gflops)

 if (do_mapping) then
   ! TODO: Use symrec conventions although this means that we cannot reuse these tables
   ! to symmetrize wavefunctions and potentials that require S-1 i.e. the symrel convention.

   ! Get mapping IBZ_k --> initial IBZ (self%lgk%ibz --> self%ibz)
   ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
   call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, exit_loop=.true., use_symrec=.False.)

   !call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, self%lgk%nsym_lg,&
   !   sppoldbl1, self%lgk%symafm_lg, self%lgk%symrec_lg, 0, comm, use_symrec=.True.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if
   ABI_SFREE(self%lgk2ibz)
   call alloc_copy(indkk(:, 1), self%lgk2ibz)
   ABI_FREE(indkk)
   call cwtime_report(" listkk1", cpu, wall, gflops)

   ! Get mapping (k + q) --> initial IBZ.
   do ii=1,self%nq_k
     self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) + kpoint
   end do
   ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))

   call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, exit_loop=.true., use_symrec=.False.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) + q could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if
   call cwtime_report(" listkk2", cpu, wall, gflops)

   ABI_SFREE(self%kq2ibz)
   call alloc_copy(indkk(:, 1), self%kq2ibz)
   ABI_FREE(indkk)
   do ii=1,self%nq_k
     self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) - kpoint
   end do
 end if

 ! Get mapping BZ --> IBZ_k (self%bz --> self%lgrp%ibz) required for tetrahedron method
 ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
 indkk(:, 1) = self%lgk%bz2ibz_smap(1, :)

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgk%kibz_q
 call self%tetra_k%free()
 call htetra_init(self%tetra_k, indkk(:, 1), cryst%gprimd, self%klatt, self%bz, self%nbz, &
                  self%lgk%ibz, self%nq_k, ierr, errorstring, comm)
 !call tetra_write(self%tetra_k, self%lgk%nibz, self%lgk%ibz, strcat("tetrak_", ktoa(kpoint)))
 ABI_CHECK(ierr == 0, errorstring)
 ABI_FREE(indkk)

 call cwtime_report(" init_tetra", cpu, wall, gflops)

end subroutine ephwg_setup_kpoint
!!***

!!****f* m_ephwg/ephwg_double_grid_setup_kpoint
!! NAME
!! ephwg_setup_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute integration weights for a given k-point
!!  using the double grid routines to map the different k-points.
!!  This version should be more efficient than its counterpart without the double grid.
!!
!! INPUTS
!!  kpoint(3): k-point in reduced coordinates.
!!  prtvol: Verbosity level
!!  comm: MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_double_grid_setup_kpoint(self, eph_doublegrid, kpoint, prtvol, comm)

!Arguments ------------------------------------
!scalars
 class(ephwg_t),target,intent(inout) :: self
 type(eph_double_grid_t),intent(inout) :: eph_doublegrid
 integer,intent(in) :: prtvol, comm
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1,timrev0=0
 integer :: ierr,ii,ik_idx
 !real(dp) :: dksqmax
 character(len=80) :: errorstring
 !character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: lgkibz2bz(:) !indkk(:,:),
 integer,allocatable :: bz2lgkibz(:), bz2lgkibzkq(:) !, bz2bz(:), mapping(:,:)
 !real(dp) :: kpt(3), wrap_kpt(3), shift
!----------------------------------------------------------------------

 cryst => self%cryst

 ! Get little group of the (external) kpoint.
 call self%lgk%free()
 self%lgk = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz, comm)
 if (prtvol > 0) call self%lgk%print()
 self%nq_k = self%lgk%nibz

 ! get dg%bz --> self%lgrp%ibz
 ABI_SFREE(eph_doublegrid%bz2lgkibz)
 ABI_MALLOC(eph_doublegrid%bz2lgkibz,(eph_doublegrid%dense_nbz))

 ! Old version using all crystal symmetries
 !call eph_doublegrid%bz2ibz(self%lgk%ibz, self%lgk%nibz,&
 !                            cryst%symrel, cryst%nsym, &
 !                            eph_doublegrid%bz2lgkibz, has_timrev=1)
 call eph_doublegrid%bz2ibz(self%lgk%ibz, self%lgk%nibz,&
                            self%lgk%symrec_lg, self%lgk%nsym_lg, &
                            eph_doublegrid%bz2lgkibz, timrev0, use_symrec=.true.)

 ! self%lgrp%ibz --> dg%bz
 ABI_ICALLOC(lgkibz2bz, (self%lgk%nibz))
 do ii=1,self%nbz
   ik_idx = eph_doublegrid%bz2lgkibz(ii)
   lgkibz2bz(ik_idx) = ii
 enddo

 ! get self%lgrp%ibz --> dg%bz --> self%ibz
 ABI_SFREE(self%lgk2ibz)
 ABI_MALLOC(self%lgk2ibz, (self%nq_k))
 do ii=1,self%nq_k
   ik_idx = lgkibz2bz(ii)
   self%lgk2ibz(ii) = eph_doublegrid%bz2ibz_dense(ik_idx)
 enddo

 ! calculate k+q
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) + kpoint
 end do

 ! get dg%bz --> lgk%ibz (k+q)
 ABI_MALLOC(bz2lgkibzkq, (eph_doublegrid%dense_nbz))
 ! Old version using all crystal symmetries
 !call eph_doublegrid%bz2ibz(self%lgk%ibz, self%lgk%nibz,&
 !                            cryst%symrel, cryst%nsym, &
 !                            bz2lgkibzkq, has_timrev=1)
 call eph_doublegrid%bz2ibz(self%lgk%ibz, self%lgk%nibz,&
                            self%lgk%symrec_lg, self%lgk%nsym_lg, &
                            bz2lgkibzkq, timrev0, use_symrec=.true.)

 ! self%lgrp%ibz (k+q) --> dg%bz
 do ii=1,self%nbz
   ik_idx = bz2lgkibzkq(ii)
   lgkibz2bz(ik_idx) = ii
 enddo
 ABI_FREE(bz2lgkibzkq)

 ! get self%lgrp%ibz (k+q) --> dg%bz --> self%ibz
 ABI_SFREE(self%kq2ibz)
 ABI_MALLOC(self%kq2ibz, (self%nq_k))
 do ii=1,self%nq_k
   ik_idx = lgkibz2bz(ii)
   self%kq2ibz(ii) = eph_doublegrid%bz2ibz_dense(ik_idx)
 enddo
 ABI_FREE(lgkibz2bz)

 ! revert change
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) - kpoint
 end do

 ! get self%bz --> dg%bz --> self%lgrp%ibz
 ABI_MALLOC(bz2lgkibz,(self%nbz))

 do ii=1,self%nbz
    ! get self%bz --> dg%bz
    ik_idx = eph_doublegrid%get_index(self%bz(:,ii),2)
    ! dg%bz --> self%lgrp%ibz
    bz2lgkibz(ii) = eph_doublegrid%bz2lgkibz(ik_idx)
 end do

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgk%kibz_q
 call self%tetra_k%free()
 call htetra_init(self%tetra_k, bz2lgkibz, cryst%gprimd, self%klatt, self%bz, self%nbz, &
                  self%lgk%ibz, self%nq_k, ierr, errorstring, comm)
 if (ierr /= 0) then
   MSG_ERROR(errorstring)
 end if
 ABI_FREE(bz2lgkibz)

end subroutine ephwg_double_grid_setup_kpoint
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_report_stats
!! NAME
!! ephwg_report_stats
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

subroutine ephwg_report_stats(self)

!Arguments ------------------------------------
!scalars
 class(ephwg_t),intent(in) :: self

!Variables
 real(dp) :: mem_tot
!----------------------------------------------------------------------

 ! IBZ qpoints
 mem_tot = 3 * self%nibz * dp
 ! BZ qpoints
 mem_tot = mem_tot + 3 * self%nbz * dp
 ! lgk2ibz and kq2ibz
 mem_tot = mem_tot + self%nq_k * 2 * 4
 ! phonon frequencies
 mem_tot = mem_tot + self%nibz * self%natom3 * dp
 ! eigenvalues
 mem_tot = mem_tot + self%nibz * self%nbcount * self%nsppol * dp

 write(std_out,"(a,f8.1,a)") " Memory allocated for ephwg:", mem_tot * b2Mb, " [Mb] <<< MEM"

end subroutine ephwg_report_stats
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_get_deltas
!! NAME
!! ephwg_get_deltas
!!
!! FUNCTION
!! Compute weights for $ \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
!! for a given (band, spin) and phonon mode nu.
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nene=number of energies for DOS
!! eminmax=min and  energy in delta (linear mesh)
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!! [broad]=Gaussian broadening
!!
!! OUTPUT
!!  deltaw_pm(nene, nq_k, 2)  (plus, minus) including the weights for BZ integration.
!   These arrays have the same order as the little group used in sigmaph.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_deltas(self, band, spin, nu, nene, eminmax, bcorr, deltaw_pm, comm, &
                            broad)  ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band, spin, nu, nene, bcorr, comm
 class(ephwg_t),intent(in) :: self
 real(dp),optional,intent(in) :: broad
!arrays
 real(dp),intent(in) :: eminmax(2)
 real(dp),intent(out) :: deltaw_pm(nene, self%nq_k, 2)

!Local variables-------------------------------
!scalars
 integer :: iq,iq_ibz,ikpq_ibz,ib,ie
 real(dp),parameter :: max_occ1 = one
 real(dp) :: omega_step
!arrays
 real(dp) :: thetaw(nene, self%nq_k), wme0(nene), pme_k(self%nq_k, 2)

!----------------------------------------------------------------------

 ib = band - self%bstart + 1

 ! Fill array for e_{k+q, b} +- w_{q,nu)
 do iq=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
   pme_k(iq, 1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
   pme_k(iq, 2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
 end do

 if (present(broad)) then
   omega_step = (eminmax(2) - eminmax(1)) / (nene - 1)
   ! Use thetaw as workspace array
   thetaw(:, 1) = arth(eminmax(1), omega_step, nene)
   do iq=1,self%nq_k
     do ie=1,2
       wme0 = thetaw(:, 1) - pme_k(iq, ie)
       deltaw_pm(:, iq, ie) = gaussian(wme0, broad)
     end do
   end do

   ! Multiply by weights
   do ie=1,nene
     deltaw_pm(ie, :, 1) = deltaw_pm(ie, :, 1) * self%lgk%weights
     deltaw_pm(ie, :, 2) = deltaw_pm(ie, :, 2) * self%lgk%weights
   end do

 else
   ! TODO Add routine to compute only delta
   call self%tetra_k%blochl_weights(pme_k(:,1), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, thetaw, deltaw_pm(:,:,1), comm)
   call self%tetra_k%blochl_weights(pme_k(:,2), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, thetaw, deltaw_pm(:,:,2), comm)
 end if

end subroutine ephwg_get_deltas
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_get_deltas_wvals
!! NAME
!! ephwg_get_deltas_wvals
!!
!! FUNCTION
!! Compute weights for $ \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
!! for a given (band, spin) and phonon mode nu.
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nene=number of energies for DOS
!! eminmax=min and  energy in delta (linear mesh)
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!! [broad]=Gaussian broadening
!!
!! OUTPUT
!!  deltaw_pm(nene, nq_k, 2)  (plus, minus) including the weights for BZ integration.
!   These arrays have the same order as the little group used in sigmaph.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_deltas_wvals(self, band, spin, nu, neig, eig, bcorr, deltaw_pm, comm, &
                                  broad)  ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band, spin, nu, neig, bcorr, comm
 real(dp),intent(in) :: eig(neig)
 class(ephwg_t),intent(in) :: self
 real(dp),optional,intent(in) :: broad
!arrays
 real(dp),intent(out) :: deltaw_pm(neig,self%nq_k, 2)

!Local variables-------------------------------
!scalars
 integer :: iq,iq_ibz,ikpq_ibz,ib
 integer :: nprocs, my_rank
 real(dp),parameter :: max_occ1 = one
 real(dp) :: wme0(neig)
!arrays
 real(dp) :: pme_k(self%nq_k,2)

!----------------------------------------------------------------------

 ib = band - self%bstart + 1
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 deltaw_pm = zero

 ! Fill array for e_{k+q, b} +- w_{q,nu)
 do iq=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
   pme_k(iq,1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
   pme_k(iq,2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
 end do

 ! Compute the tetrahedron or gaussian weights
 if (present(broad)) then
   do iq_ibz=1,self%nq_k
     if (mod(iq_ibz, nprocs) /= my_rank) cycle ! MPI parallelism
     wme0 = eig - pme_k(iq_ibz, 1)
     deltaw_pm(:,iq_ibz,1) = gaussian(wme0, broad) * self%lgk%weights(iq_ibz)
     wme0 = eig - pme_k(iq_ibz, 2)
     deltaw_pm(:,iq_ibz,2) = gaussian(wme0, broad) * self%lgk%weights(iq_ibz)
   end do
 else
   call self%tetra_k%wvals_weights_delta(pme_k(:,1),neig,eig,max_occ1,self%nq_k,bcorr,deltaw_pm(:,:,1),comm)
   call self%tetra_k%wvals_weights_delta(pme_k(:,2),neig,eig,max_occ1,self%nq_k,bcorr,deltaw_pm(:,:,2),comm)
 end if

end subroutine ephwg_get_deltas_wvals
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_get_deltas_qibzk
!! NAME
!! ephwg_get_deltas_qibzk
!!
!! FUNCTION
!! Compute weights for $ \delta(\omega - \omega_{q\nu} $ for given nu, using q-opints in the IBZ(k)
!!
!! INPUTS
!! nu=Phonon branch.
!! nene=number of energies for DOS
!! eminmax=min and  energy in delta (linear mesh)
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!! [with_qweights]= .False. if q-point weights should not be included in dt_weights
!!
!! OUTPUT
!!  dt_weights(nene, nq_k, 2)  weights for BZ integration (delta and theta function)
!   These arrays have the same order as the q-points in the little group of the k-point.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_deltas_qibzk(self, nu, nene, eminmax, bcorr, dt_weights, comm, with_qweights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nu, nene, bcorr, comm
 class(ephwg_t),intent(in) :: self
 logical,optional,intent(in) :: with_qweights
!arrays
 real(dp),intent(in) :: eminmax(2)
 real(dp),intent(out) :: dt_weights(nene, self%nq_k, 2)

!Local variables-------------------------------
!scalars
 integer :: iq, iq_ibz, ie, ii
 real(dp),parameter :: max_occ1 = one
!arrays
 real(dp) :: eigen_in(self%nq_k)

!----------------------------------------------------------------------

 ! Fill eigen_in
 do iq=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq)  ! IBZ_k --> IBZ
   eigen_in(iq) = self%phfrq_ibz(iq_ibz, nu)
 end do

 call self%tetra_k%blochl_weights(eigen_in, eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
   bcorr, dt_weights(:,:,2), dt_weights(:,:,1), comm)

 if (present(with_qweights)) then
  if (.not. with_qweights) then
    do ii=1,2
      do ie=1,nene
        dt_weights(ie, :, ii) = dt_weights(ie, :, ii) * self%lgk%weights
      end do
    end do
   end if
 end if

end subroutine ephwg_get_deltas_qibzk
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_get_zinv_weights
!! NAME
!! ephwg_get_zinv_weights
!!
!! FUNCTION
!! Compute weights for a given (kpoint, qpoint, spin) for all phonon modes.
!!
!! INPUTS
!! qpt(3)
!! iband_sum=band index in self-energy sum. (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nbcalc=Number of bands in self-energy matrix elements.
!! zvals
!! opt: 1 for S. Kaprzyk routines, 2 for Lambin.
!! comm=MPI communicator
!! [use_bzsum]= By default the weights are multiplied by the Nstar(q) / Nq where
!!   Nstar(q) is the number of points in the star of the q-point (using the symmetries of the little group of k)
!!   If use_bzsum is set to True, the Nstar(q) coefficient is removed so that the caller can
!!   integrate over the BZ without using symmetries.
!!
!! OUTPUT
!!  cweights(nz, 2, nbcalc, %nq_k)  (plus, minus)
!!  include weights for BZ integration.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_zinv_weights(self, nz, nbcalc, zvals, iband_sum, spin, nu, opt, cweights, comm, use_bzsum)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband_sum, spin, nu, nz, nbcalc, opt, comm
 class(ephwg_t),intent(in) :: self
 logical, optional, intent(in) :: use_bzsum
!arraye
 complex(dpc),intent(in) :: zvals(nz, nbcalc)
 complex(dpc),intent(out) :: cweights(nz, 2, nbcalc, self%nq_k)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iq_ibz,ikpq_ibz,ib,ii,iq,nprocs, my_rank, ierr
 real(dp),parameter :: max_occ1=one
 logical :: use_bzsum_
!arrays
 real(dp),allocatable :: pme_k(:,:)
 complex(dp),allocatable :: cweights_tmp(:,:)
!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 use_bzsum_ = .False.; if (present(use_bzsum)) use_bzsum_ = use_bzsum

 ! Allocate array for e_{k+q, b} +- w_{q,nu)
 ABI_MALLOC(pme_k, (self%nq_k, 2))

 ib = iband_sum - self%bstart + 1
 do iq=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
   pme_k(iq, 1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
   pme_k(iq, 2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
 end do

 cweights = zero
 ABI_MALLOC(cweights_tmp,(nz,self%nq_k))
 do ib=1,nbcalc
   do ii=1,2
     call self%tetra_k%weights_wvals_zinv(pme_k(:, ii), nz, zvals(:,ib), &
                                          max_occ1, self%nq_k, opt, cweights_tmp, comm)
     do iq=1,self%nq_k
       cweights(:,ii,ib,iq) = cweights_tmp(:,iq)
     end do
   end do
 end do
 ABI_FREE(cweights_tmp)

 ! Rescale weights so that the caller can sum over the full BZ.
 !if (use_bzsum_) cweights = cweights / ( self%lgk%weights(iqlk) * self%nbz )

 call xmpi_sum(cweights, comm, ierr)

 ! Compare results with naive one-point integration.
 !if (.False. .and. my_rank == master) then
 !  volconst_mult = one
 !  volconst_mult = self%lgk%weights(iqlk)
 !  do ib=1,nbcalc
 !    !do nu=1,self%natom3
 !      write(std_out,*)"# naive vs tetra integration for band, nu", ib - 1 + self%bstart, nu
 !      do iz=1,nz
 !        write(std_out,"(5es16.8)") &
 !          dble(zvals(iz, ib)), one / (zvals(iz, ib) - pme_k(iqlk, 1)) * volconst_mult, cweights(iz, 1, ib, nu)
 !        write(std_out,"(5es16.8)") &
 !          dble(zvals(iz, ib)), one / (zvals(iz, ib) - pme_k(iqlk, 2)) * volconst_mult, cweights(iz, 2, ib, nu)
 !      end do
 !    !end do
 !  end do
 !end if

 ABI_FREE(pme_k)

end subroutine ephwg_get_zinv_weights
!!***

!!****f* m_ephwg/ephwg_free
!! NAME
!! ephwg_free
!!
!! FUNCTION
!!  Deallocate memory
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

subroutine ephwg_free(self)

!Arguments ------------------------------------
 class(ephwg_t),intent(inout) :: self

!----------------------------------------------------------------------

 ! integer
 ABI_SFREE(self%kq2ibz)

 ! Real
 ABI_SFREE(self%ibz)
 ABI_SFREE(self%bz)
 ABI_SFREE(self%lgk2ibz)
 ABI_SFREE(self%phfrq_ibz)
 ABI_SFREE(self%eigkbs_ibz)

 ! types
 call self%tetra_k%free()
 call self%lgk%free()

 ! nullify pointers
 self%cryst => null()

end subroutine ephwg_free
!!***

end module m_ephwg
!!***
