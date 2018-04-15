!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ephwg
!! NAME
!! m_ephwg
!!
!! FUNCTION
!!  Tools and objects to compute the weights used for the BZ integration of EPH quantities.
!!  More specifically the integration of quantities such as the imaginary part of the self-energy
!!  involving delta functions.
!!  Different approaches are available:
!!
!!    1.
!!    2.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
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
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_copy
 use m_tetrahedron
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_crystal
 use m_crystal_io
 use m_ifc
 use m_lgroup
 use m_ebands

 use defs_datatypes,    only : ebands_t
 use m_time,            only : cwtime, sec2str
 use m_fstrings,        only : strcat
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
!!  Stores electron eigevalues and phonon frequencies in the IBZ.
!!  (use same mesh for e and ph)
!!  Provides tools to compute (e_{k+q} - w{q}) in the IBZ(k)
!!  and integrate the delta functions with the tetrahedron method.
!!
!! SOURCE

type, public :: ephwg_t

  integer :: natom3
  ! 3 * natom

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: nbcount
  ! Number of bands treated

  integer :: bstart
  ! The fist band starts at bstart. Used to select bands
  ! aroudh the Fermi level.

  integer :: kptopt
  ! Option for k-point generation.

  integer :: timrev
  ! 1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nibz, nbz
  ! Number of q-points in IBZ and full BZ.

  integer :: nq_k
  ! Nunber of points in IBZ(k) i.e. the irreducible wedge
  ! defined by the operations of the little group of k.

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
  ! reciprocal of lattice vectors for full kpoint grid.
  ! Used by init_tetra

  integer,allocatable :: lgk2ibz(:)
  ! lgk2ibz(nq_k)
  ! Mapping Litte-group --> initial IBZ

  real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states (not necessarly equal to nband)

  type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure

  type(lgroup_t) :: lgrp_k
  ! Little group of the k-point

  type(t_tetrahedron) :: tetra_k
  ! tetrahedron object used to evaluate delta(w - e_{k+q} - phw_q)

 end type ephwg_t

 public :: ephwg_new             ! Constructor
 public :: ephwg_set_kpoint      ! Prepare tetrahedron method for k-point
 public :: ephwg_getwk           ! Compute delta functions $\delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
 public :: ephwg_free            ! Free memory

 public :: ephwg_test

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
!!  bstart
!!  nbcount
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_new'
!End of the abilint section

 implicit none

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

 ! Get full BZ (new%nbz, new%bz) and new kptrlatt for tetra.
 call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, out_nkibz, out_kibz, out_wtk, new%nbz, new%bz, &
   new_kptrlatt=out_kptrlatt)

 rlatt = out_kptrlatt; call matr3inv(rlatt, new%klatt)

 ABI_CHECK(size(out_kibz, dim=2) == new%nibz, "mismatch in nkibz!")
 ABI_FREE(out_kibz)
 ABI_FREE(out_wtk)

 ! Copy eigenvalues in IBZ. Change shape for better performance.
 ABI_MALLOC(new%eigkbs_ibz, (new%nibz, new%nbcount, new%nsppol))
 do ik=1,new%nibz
   new%eigkbs_ibz(ik, :, :) = eig_ibz(:, ik, :)
 end do

 ! Fourier interpolate frequencies on the same mesh.
 ABI_CALLOC(new%phfrq_ibz, (new%nibz, new%natom3))
 do ik=1,new%nibz
   if (mod(ik, nprocs) /= my_rank) cycle ! mpi-parallelism
   call ifc_fourq(ifc, cryst, new%ibz(:, ik), phfrq, displ_cart)
   new%phfrq_ibz(ik, :) = phfrq
 end do
 call xmpi_sum(new%phfrq_ibz, comm, ierr)

end function ephwg_new
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_set_kpoint
!! NAME
!! ephwg_set_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute the delta functions with a given k-point.
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

subroutine ephwg_set_kpoint(self, kpoint)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_set_kpoint'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ephwg_t),target,intent(inout) :: self
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1
 integer :: ierr,ii
 real(dp) :: dksqmax
 character(len=80) :: errorstring
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: indkk(:,:)

!----------------------------------------------------------------------

 cryst => self%cryst

 ! Get little group of the kpoint.
 self%lgrp_k = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz)
 self%nq_k = self%lgrp_k%nibz

 ! Get mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgrp_k%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6)' ) &
    "At least one of the points in IBZ(k) could not be generated from a symmetrical one. dksqmax: ",dksqmax
   MSG_ERROR(msg)
 end if
 call alloc_copy(indkk(:, 1), self%lgk2ibz)
 ABI_FREE(indkk)

 ! Get mapping (k + q) --> initial IBZ.
 do ii=1,self%nq_k
   self%lgrp_k%ibz(:, ii) = self%lgrp_k%ibz(:, ii) + kpoint
 end do
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))

 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgrp_k%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6)' ) &
    "At least one of the points in IBZ(k) + k could not be generated from a symmetrical one. dksqmax: ",dksqmax
   MSG_ERROR(msg)
 end if
 call alloc_copy(indkk(:, 1), self%kq2ibz)
 ABI_FREE(indkk)
 do ii=1,self%nq_k
   self%lgrp_k%ibz(:, ii) = self%lgrp_k%ibz(:, ii) - kpoint
 end do

 ! Get mapping BZ --> IBZ_k (self%bz --> self%lgrp%ibz) required for tetrahedron method
 ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%lgrp_k%ibz, self%bz, self%nq_k, self%nbz, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgrp_k%kibz_q
 call init_tetra(indkk(:, 1), cryst%gprimd, self%klatt, self%bz, self%nbz, self%tetra_k, ierr, errorstring)
 if (ierr /= 0) then
   MSG_ERROR(errorstring)
 end if
 ABI_FREE(indkk)

end subroutine ephwg_set_kpoint
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_getwk
!! NAME
!! ephwg_getwk
!!
!! FUNCTION
!!  Compute delta functions for a given (band, spin) and phonon mode nu
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nene=number of energies for DOS
!! eminmax=min and  energy in delta (linear mesh)
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!!
!! OUTPUT
!!  wdt_plus(nene, nq_k, 2)
!!  wdt_minus(nene, nq_k, 2)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_getwk(self, band, spin, nu, nene, eminmax, bcorr, wdt_plus, wdt_minus, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_getwk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band, spin, nu, nene, bcorr, comm
 type(ephwg_t),intent(in) :: self
!arrays
! This arrays have the same order as the little group used in sigmaph.
 real(dp),intent(in) :: eminmax(2)
 real(dp),intent(out) :: wdt_plus(nene, self%nq_k, 2), wdt_minus(nene, self%nq_k, 2)

!Local variables-------------------------------
!scalars
 integer :: iq_lk,iq_ibz,ikpq_ibz,ib
 real(dp),parameter :: max_occ1 = one
!arrays
 real(dp),allocatable :: pme_k(:,:)

!----------------------------------------------------------------------

 ! Allocate array for e_{k+q, b} +- w_{q,nu)
 ABI_MALLOC(pme_k, (self%nq_k, 2))
 ib = band - self%bstart + 1

 do iq_lk=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq_lk)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq_lk)  ! k + q --> IBZ
   pme_k(iq_lk, 1) =  self%phfrq_ibz(iq_ibz, nu) - self%eigkbs_ibz(ikpq_ibz, ib, spin)
   pme_k(iq_lk, 2) = -self%phfrq_ibz(iq_ibz, nu) - self%eigkbs_ibz(ikpq_ibz, ib, spin)
 end do

 call tetra_blochl_weights(self%tetra_k, pme_k(:,1), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
   bcorr, wdt_plus(:,:,2), wdt_plus(:,:,1), comm)

 ABI_FREE(pme_k)

 !call tetra_get_onewk(self%tetra_k, ik_ibz, bcorr, nene, nkibz, eig_ibz, eminmax(1), eminmax(2), max_occ1, weights)

end subroutine ephwg_getwk
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_free
!! NAME
!! ephwg_free
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

subroutine ephwg_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ephwg_t),intent(inout) :: self

!----------------------------------------------------------------------

 ! integer
 if (allocated(self%kq2ibz)) then
   ABI_FREE(self%kq2ibz)
 end if

 ! Real
 if (allocated(self%ibz)) then
   ABI_FREE(self%ibz)
 end if
 if (allocated(self%bz)) then
   ABI_FREE(self%bz)
 end if
 if (allocated(self%lgk2ibz)) then
   ABI_FREE(self%lgk2ibz)
 end if
 if (allocated(self%phfrq_ibz)) then
   ABI_FREE(self%phfrq_ibz)
 end if
 if (allocated(self%eigkbs_ibz)) then
   ABI_FREE(self%eigkbs_ibz)
 end if

 ! types
 call destroy_tetra(self%tetra_k)
 call lgroup_free(self%lgrp_k)

end subroutine ephwg_free
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_free
!! NAME
!! ephwg_free
!!
!! FUNCTION
!!
!! INPUTS
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_test(cryst, ebands, ifc, prefix, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_free'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalarster
 integer,intent(in) :: comm
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
 character(len=*),intent(in) :: prefix

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0 = 0, intp_nshiftk = 1, master = 0
 integer :: band, spin, nu, bstart, nbcount, ik, nene, nprocs, my_rank, cnt, ierr
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 real(dp) :: omega_step, kT !, f_mkq, nqnu
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 type(ephwg_t) :: ephwg
 type(ebands_t) :: eb_dense
!arrays
 integer :: intp_kptrlatt(3, 3), band_block(2)
 real(dp):: params(4), intp_shiftk(3, intp_nshiftk), eminmax(2)
 real(dp),allocatable :: wdt_plus(:,:,:), wdt_minus(:,:,:), sdelta(:,:,:), omega_mesh(:)

!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 !call ebands_interpolate_kpath(QP_BSt, dtset, cryst, [sigp%minbdgw, sigp%maxbdgw], dtfil%filnam_ds(4), comm)
 !eb_kpath = ebands_interp_kpath(ebands, cryst, kpath, params, band_block, comm)

 ! Interpolate band energies with star-functions
 params = 0; params(1) = 1; params(2) = 5
 band_block = [1, 8]
 intp_kptrlatt = reshape([8, 0, 0, 0, 8, 0, 0, 0, 8], [3, 3])
 intp_shiftk = zero

 eb_dense = ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt, intp_nshiftk, intp_shiftk, band_block, comm)

 ! Init object for eph weights.
 bstart = band_block(1); nbcount = band_block(2) - band_block(1) + 1
 ephwg = ephwg_new(cryst, ifc, bstart, nbcount, eb_dense%kptopt, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
   eb_dense%nkpt, eb_dense%kptns, eb_dense%nsppol, eb_dense%eig, comm)

 eminmax(1) = minval(eb_dense%eig(band_block(1):band_block(2), :, :))
 eminmax(2) = maxval(eb_dense%eig(band_block(1):band_block(2), :, :))
 omega_step = 0.05_dp * eV_Ha
 nene = 1 + nint((eminmax(2) - eminmax(1)) / omega_step)
 ABI_MALLOC(omega_mesh, (nene))
 do ik = 1, nene
   omega_mesh(ik) = eminmax(1) + (ik - 1) * omega_step
 end do

 kT = zero
 ABI_CALLOC(sdelta, (nene, eb_dense%nkpt, eb_dense%nsppol))
 cnt = 0
 do ik=1,eb_dense%nkpt
   cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
   call cwtime(cpu, wall, gflops, "start")
   call ephwg_set_kpoint(ephwg, eb_dense%kptns(:, ik))
   ABI_MALLOC(wdt_plus, (nene, ephwg%nq_k, 2))
   ABI_MALLOC(wdt_minus, (nene, ephwg%nq_k, 2))
   do nu = 1, ifc%natom * 3
     !nqnu = occ_be(wqnu, kT, zero)
     do band = band_block(1), band_block(2)
       do spin = 1, eb_dense%nsppol
         call ephwg_getwk(ephwg, band, spin, nu, nene, eminmax, bcorr0, wdt_plus, wdt_minus, xmpi_comm_self)
         !f_mkq = occ_fd(ep_edens%eig(band, ikq_ibz, spin), kT, ebands%fermie)
         !wdt_plus = wdt_plus * ephwg%lgrp_k%weights
         sdelta(:, ik, spin) = sum(wdt_plus(:, :, 1), dim=2) + sum(wdt_minus(:, :, 1), dim=2)
       end do
     end do
   end do
   ABI_FREE(wdt_plus)
   ABI_FREE(wdt_minus)
   call cwtime(cpu, wall, gflops, "stop")
   write(msg,'(2(a,i0),2(a,f8.2))')"k-point [", ik, "/", eb_dense%nkpt, "] completed. cpu:",cpu,", wall:",wall
   call wrtout(std_out, msg, do_flush=.True.)
 end do

 ! Collect data, then master writes data.
 call xmpi_sum_master(sdelta, master, comm, ierr)
 if (my_rank == master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ncid, strcat(prefix, "_SDELTA.nc"), xmpi_comm_self))
   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   NCF_CHECK(ebands_ncwrite(eb_dense, ncid))
   ! Add dimensions.
   NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nomega", nene)], defmode=.True.))
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("omega_mesh", "dp", "nomega"), &
     nctkarr_t("sdelta", "dp", "nomega, number_of_kpoints, number_of_spins") &
   ])
   NCF_CHECK(ncerr)
   ! Write data
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "omega_mesh"), omega_mesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sdelta"), sdelta))
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 ABI_FREE(omega_mesh)
 ABI_FREE(sdelta)
 call ebands_free(eb_dense)
 call ephwg_free(ephwg)

end subroutine ephwg_test
!!***

end module m_ephwg
!!***
