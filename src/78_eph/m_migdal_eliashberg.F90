!!****m* ABINIT/m_migdal_eliashberg
!! NAME
!! m_migdal_eliashberg
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

module m_migdal_eliashberg

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_krank
 use m_htetra

 !use m_ebands
 !use iso_c_binding
 use netcdf
 use m_nctk
 !use m_crystal
 use m_dtset
 use m_dtfil
 !use m_ephtk

 use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : strcat !,sjoin, tolower, itoa, ftoa, ktoa, ltoa, strcat
 !use m_numeric_tools,  only : arth, get_diag !, isdiagmat
 use m_copy,           only : alloc_copy
 !use defs_datatypes,   only : ebands_t
 !use m_kpts,           only : kpts_timrev_from_kptopt
 !use m_ifc,            only : ifc_type
 use m_gstore, only : gstore_t

 implicit none

 private

 public :: migdal_eliashberg_iso
 !public :: migdal_eliashberg_aniso
!!***

!----------------------------------------------------------------------

!!****t* m_migdal_eliashberg/iso_solver_t
!! NAME
!! iso_solver_t
!!
!! FUNCTION
!!
!! NOTES
!!
!! SOURCE

type, public :: iso_solver_t

  integer :: ntemp

  integer :: max_niter = -1
  ! Maximum number of iterations.

  integer :: max_nmix = 4

  integer :: ncid

  integer :: comm = xmpi_undefined

  !integer :: niw = -1
  ! Number of Matsubara frequencies.

  !real(dp) :: kt = -one
  ! K * T in Ha.

  real(dp) :: tolerance = -one

  real(dp),allocatable :: zeta_iw(:), delta_iw(:)
  real(dp),allocatable :: delta_iw_mix(:,:)

contains

  procedure :: free => iso_solver_free
  ! Free dynamic memory

  procedure :: solve => iso_solver_solve

end type iso_solver_t

!private :: iso_solver_new

!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/iso_solver_free
!! NAME
!! iso_solver_free
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

subroutine iso_solver_free(solver)

!Arguments ------------------------------------
 class(iso_solver_t),intent(inout) :: solver
!----------------------------------------------------------------------

 ABI_SFREE(solver%zeta_iw)
 ABI_SFREE(solver%delta_iw)
 ABI_SFREE(solver%delta_iw_mix)

end subroutine iso_solver_free
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/iso_solver_solve
!! NAME
!! iso_solver_solve
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

subroutine iso_solver_solve(solver, itemp, kt, niw, imag_w, lambda_ij)

!Arguments ------------------------------------
!scalars
 class(iso_solver_t),intent(inout) :: solver
 integer,intent(in) :: itemp, niw
 real(dp),intent(in) :: kt
!arrays
 real(dp),intent(in) :: imag_w(niw), lambda_ij(2 * niw)

!Local variables-------------------------------
!scalars
 integer :: iter !, ii
!arrays
 real(dp),allocatable :: prev_vals(:)

!----------------------------------------------------------------------

 ABI_REMALLOC(solver%delta_iw_mix, (niw, solver%max_nmix))

 if (itemp == 1) then
   ! Init values from scratch
   ABI_MALLOC(solver%zeta_iw, (niw))
   ABI_MALLOC(solver%delta_iw, (niw))
 else
   ! Init values from previous temp. TODO: May use spline
   call alloc_copy(solver%zeta_iw, prev_vals)
   ABI_MOVE_ALLOC(prev_vals, solver%zeta_iw)
   call alloc_copy(solver%delta_iw, prev_vals)
   ABI_MOVE_ALLOC(prev_vals, solver%delta_iw)
 end if

 do iter=1,solver%max_niter

   ! Write SCF cycle to stdout.

   ! Check for convergence.

   ! Mix delta.

 end do ! iter

 ! Pade' to go to real axis
 ! Compute Delta F
 ! Compute QP DOS

 ! Write results to netcdff

end subroutine iso_solver_solve
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/migdal_eliashberg_iso
!! NAME
!! migdal_eliashber_iso
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

subroutine migdal_eliashberg_iso(gstore, dtset, dtfil)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(gstore_t),target,intent(inout) :: gstore

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nproc, my_rank, ierr, itemp, ntemp, niw, ncid !, my_nshiftq, nsppol !, iq_glob, ik_glob, ii ! out_nkibz,
 !integer :: spin, natom3, cnt !, band, ib, nb, my_ik, my_iq, my_is,
 !integer :: ik_ibz, ik_bz, ebands_timrev, max_nq, max_nk, gstore_cplex_
 !integer :: iq_bz, iq_ibz !, ikq_ibz, ikq_bz
 !integer :: ncid, spin_ncid, ncerr, gstore_fform
 real(dp) :: kt, cpu, wall, gflops
 character(len=5000) :: msg
 !class(crystal_t),target,intent(in) :: cryst
 !class(ebands_t),target,intent(in) :: ebands
 !class(ifc_type),target,intent(in) :: ifc
 !type(gqk_t),pointer :: gqk
 !type(krank_t) :: qrank
 type(iso_solver_t) :: solver
!arrays
 real(dp),allocatable :: ktmesh(:), lambda_ij(:), imag_w(:), imag_2w(:)

!----------------------------------------------------------------------

 call cwtime(cpu, wall, gflops, "start")

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)

 !cryst => gstore%cryst
 !ebands => gstore%ebands
 !ifc => gstore%ifc
 !kibz => gstore%kibz
 !natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "bz", "migdal_eliashberg_iso requires qzone == 'bz'", ierr)
 ABI_CHECK(ierr == 0, "Wrong GSTORE.nc file for migdal_eliashberg_iso. See messages above")

 ncid = nctk_noid
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_ISOME.nc") , xmpi_comm_self))
 end if

 call dtset%get_ktmesh(ntemp, ktmesh)
 solver = iso_solver_t(ntemp=ntemp, max_niter=10, tolerance=tol10, ncid=ncid, comm=gstore%comm)

 do itemp=1,ntemp
   ! Generate Matsubara imaginary-mesh for this T.
   kt = ktmesh(itemp)
   !call matsubara_mesh(kt, mats_wcut, niw, imag_w)
   niw = 1

   ! Compute lambda(w_i - w_j)
   ABI_MALLOC(lambda_ij, (2 * niw))
   ABI_MALLOC(imag_2w, (2 * niw))
   call gstore%get_lambda_iso_iw(2 * niw, imag_2w, lambda_ij)
   ABI_FREE(imag_2w)

   call solver%solve(itemp, kt, niw, imag_w, lambda_ij)
   ABI_FREE(lambda_ij)
   ABI_FREE(imag_w)
 end do ! itemp

 call solver%free()
 ABI_FREE(ktmesh)

 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if

 call cwtime_report(" migdal_eliashberg_iso:", cpu, wall, gflops)

end subroutine migdal_eliashberg_iso
!!**

end module m_migdal_eliashberg
!!***
