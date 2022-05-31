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
 use netcdf
 use m_nctk
 use m_crystal
 use m_dtset
 use m_dtfil

 use m_time,            only : cwtime, cwtime_report, sec2str
 use m_fstrings,        only : strcat, sjoin !, tolower, itoa, ftoa, ktoa, ltoa, strcat
 !use m_numeric_tools,  only : arth, get_diag
 use m_copy,            only : alloc_copy
 use defs_datatypes ,   only : ebands_t
 !use m_kpts,           only : kpts_timrev_from_kptopt
 use m_ebands,          only : edos_t, ebands_get_edos
 use m_gstore,          only : gstore_t

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
  real(dp),allocatable :: prev_zeta_iw(:), prev_delta_iw(:)
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

subroutine iso_solver_free(iso)

!Arguments ------------------------------------
 class(iso_solver_t),intent(inout) :: iso
!----------------------------------------------------------------------

 ABI_SFREE(iso%delta_iw)
 ABI_SFREE(iso%zeta_iw)
 ABI_SFREE(iso%prev_delta_iw)
 ABI_SFREE(iso%prev_zeta_iw)
 ABI_SFREE(iso%delta_iw_mix)

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

subroutine iso_solver_solve(iso, itemp, kt, niw, imag_w, lambda_ij)

!Arguments ------------------------------------
!scalars
 class(iso_solver_t),intent(inout) :: iso
 integer,intent(in) :: itemp, niw
 real(dp),intent(in) :: kt
!arrays
 real(dp),intent(in) :: imag_w(niw), lambda_ij(2 * niw)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nproc, my_rank, iter, ii, jj, converged
 real(dp) :: rr
!arrays
 real(dp),allocatable :: prev_vals(:)

!----------------------------------------------------------------------

 nproc = xmpi_comm_size(iso%comm); my_rank = xmpi_comm_rank(iso%comm)

 ABI_REMALLOC(iso%delta_iw_mix, (niw, iso%max_nmix))

 if (itemp == 1) then
   ! Init values from scratch
   ABI_CALLOC(iso%zeta_iw, (niw))
   ABI_CALLOC(iso%prev_zeta_iw, (niw))
   ABI_CALLOC(iso%delta_iw, (niw))
   ABI_CALLOC(iso%prev_delta_iw, (niw))
 else
   ! Init values from previous temperature. TODO: May use spline
   call alloc_copy(iso%zeta_iw, prev_vals)
   ABI_RECALLOC(iso%zeta_iw, (niw))
   ABI_MOVE_ALLOC(prev_vals, iso%prev_zeta_iw)
   call alloc_copy(iso%delta_iw, prev_vals)
   ABI_RECALLOC(iso%delta_iw, (niw))
   ABI_MOVE_ALLOC(prev_vals, iso%prev_delta_iw)
 end if

 converged = 0
iter_loop: do iter=1,iso%max_niter

   do ii=1,niw
     !if (mod(ii, nproc) /= my_rank) cycle ! MPI parallelism inside comm
     do jj=1,niw
       rr = one / sqrt(imag_w(jj) ** 2 + iso%prev_delta_iw(jj) ** 2)
       iso%zeta_iw(ii) = iso%zeta_iw(ii) + imag_w(jj) * rr  !* lambda(ii - jj)
       iso%delta_iw(ii) = iso%delta_iw(ii) + rr * iso%prev_delta_iw(jj) !* (lambda(ii - jj) - mustar)
     end do
      iso%zeta_iw(ii) = one + pi * kt / imag_w(ii) * iso%zeta_iw(ii)
      iso%delta_iw(ii) = pi * kt * iso%delta_iw(ii) / iso%zeta_iw(ii)
   end do ! ii

   if (my_rank == master) then
     ! Write SCF cycle to stdout.
     ! Check for convergence.
     converged = 0
   end if

   if (converged == 2) exit iter_loop

   ! TODO: Mixing
   iso%prev_zeta_iw = iso%zeta_iw
   iso%prev_delta_iw = iso%delta_iw

 end do iter_loop

 ! Pade' to go to real axis
 ! Compute Delta F
 ! Compute QP DOS

 ! Write results to netcdf file
 if (my_rank == master) then
 end if

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
 integer :: edos_intmeth
 !integer :: spin, natom3, cnt !, band, ib, nb, my_ik, my_iq, my_is
 !integer :: ik_ibz, ik_bz, ebands_timrev
 !integer :: iq_bz, iq_ibz !, ikq_ibz, ikq_bz
 !integer :: ncid, spin_ncid, ncerr, gstore_fform
 integer :: phmesh_size, iw
 real(dp) :: kt, wmax, cpu, wall, gflops
 real(dp) :: edos_step, edos_broad !, sigma, ecut, eshift, eig0nk
 !character(len=5000) :: msg
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 !class(ifc_type),target,intent(in) :: ifc
 !type(gqk_t),pointer :: gqk
 type(iso_solver_t) :: iso
 type(edos_t) :: edos
!arrays
 real(dp),allocatable :: ktmesh(:), lambda_ij(:), imag_w(:), imag_2w(:), phmesh(:), a2fw(:)

!----------------------------------------------------------------------

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)

 call wrtout(std_out, " Solving isotropic Migdal-Eliashberg equations on the imaginary axis", pre_newlines=2)
 call cwtime(cpu, wall, gflops, "start")

 cryst => gstore%cryst
 ebands => gstore%ebands
 !natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "bz", "qzone == 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 1, "cplex == 1 is required", ierr)
 ABI_CHECK(ierr == 0, "Wrong gstore object for migdal_eliashberg_iso. See messages above")

 ! Compute electron DOS.
 edos_intmeth = 2; if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 edos_step = 0.01 * eV_Ha; edos_broad = 0.3 * eV_Ha
 edos = ebands_get_edos(ebands, cryst, edos_intmeth, edos_step, edos_broad, gstore%comm)

 !! Store DOS per spin channel
 !n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos%print(unit=std_out)
   !call edos%print(unit=ab_out)
   !path = strcat(dtfil%filnam_ds(4), "_EDOS")
   !call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path, ch10))
   !call edos%write(path)
 end if

 ! Compute phonon frequency mesh.
 call gstore%ifc%get_phmesh(dtset%ph_wstep, phmesh_size, phmesh)

 ! Compute and store my phonon quantities
 call gstore%calc_my_phonons(store_phdispl=.False.)

 ! Compute Eliashberg function a2F(w)
 ABI_MALLOC(a2fw, (phmesh_size))
 call gstore%get_a2fw(dtset, phmesh_size, phmesh, a2fw)

 ncid = nctk_noid
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_ISOME.nc") , xmpi_comm_self))
   write(777, *)"# phmesh (meV), a2fw"
   do iw=1, phmesh_size
     write(777, *) phmesh(iw) * Ha_meV, a2fw(iw) / (edos%gef(0) / two)
   end do
 end if

 ABI_FREE(a2fw)
 ABI_FREE(phmesh)
 call edos%free()

 call dtset%get_ktmesh(ntemp, ktmesh)
 iso = iso_solver_t(ntemp=ntemp, max_niter=10, tolerance=tol10, ncid=ncid, comm=gstore%comm)

 do itemp=1,ntemp
   ! Generate Matsubara mesh for this T with cutoff wmax.
   kt = ktmesh(itemp)
   wmax = one
   call matsubara_mesh("bosons", kt, wmax, niw, imag_w)

   ! Compute lambda(w_i - w_j)
   ABI_MALLOC(lambda_ij, (2 * niw))
   ABI_MALLOC(imag_2w, (2 * niw))

   !call wrtout(std_out, " Computing lambda_iso_iw...")
   !call gstore%get_lambda_iso_iw(dtset, 2 * niw, imag_2w, lambda_ij)
   ABI_FREE(imag_2w)

   call iso%solve(itemp, kt, niw, imag_w, lambda_ij)
   ABI_FREE(lambda_ij)
   ABI_FREE(imag_w)
 end do ! itemp

 ABI_FREE(ktmesh)
 call iso%free()

 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if

 call cwtime_report(" migdal_eliashberg_iso:", cpu, wall, gflops)

end subroutine migdal_eliashberg_iso
!!!**

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/matsubara_mesh
!! NAME
!! mastubara_mesh
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

subroutine matsubara_mesh(bosons_or_fermions, kt, wmax, niw, imag_w)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: bosons_or_fermions
 real(dp),intent(in) :: kt, wmax
 integer,intent(out) :: niw
 real(dp),allocatable,intent(out) :: imag_w(:)

!Local variables-------------------------------
!scalars
 integer :: nn

!----------------------------------------------------------------------

 select case (bosons_or_fermions)
 case ("bosons")
   ! 2 n pi kT
   niw = nint(wmax / (two * pi * kt)) + 1
   ABI_MALLOC(imag_w, (niw))
   do nn=0,niw-1
     imag_w(nn) = two * nn * pi * kt
   end do

 case ("fermions")
   ! 2 (n + 1) pi kT
   niw = nint(wmax / (two * pi * kt))
   ABI_MALLOC(imag_w, (niw))
   do nn=0,niw-1
     imag_w(nn) = two * (nn  + 1) * pi * kt
   end do

 case default
   ABI_ERROR(sjoin("Wrong bosons_or_fermions:", bosons_or_fermions))
 end select

end subroutine matsubara_mesh
!!!***

end module m_migdal_eliashberg
!!***
