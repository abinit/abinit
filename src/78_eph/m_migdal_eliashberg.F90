!!****m* ABINIT/m_migdal_eliashberg
!! NAME
!! m_migdal_eliashberg
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2026 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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
 use m_fstrings,        only : strcat, sjoin, itoa, ftoa, ktoa, ltoa
 use m_copy,            only : alloc_copy
 use m_special_funcs,   only : gaussian
 use m_ebands,          only : ebands_t, edos_t
 use m_kpts,            only : kpts_timrev_from_kptopt
 use m_lgroup,          only : lgroup_t
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
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/iso_solver_free
!! NAME
!! iso_solver_free
!!
!! FUNCTION
!!  Free dynamic memory
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

 ABI_UNUSED(lambda_ij)

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
 integer :: nproc, my_rank, ierr, itemp, ntemp, niw, ncid
 integer :: edos_intmeth
 !integer :: spin, natom3, cnt !, band, ib, nb, my_ik, my_iq, my_is
 !integer :: ik_ibz, ik_bz, ebands_timrev, iq_bz, iq_ibz !, ikq_ibz, ikq_bz
 !integer :: ncid, spin_ncid, ncerr, gstore_fform
 integer :: phmesh_size, units(2) !, iw
 real(dp) :: kt, wmax, cpu, wall, gflops, edos_step, edos_broad !, sigma, ecut, eshift, eig0nk
 !character(len=5000) :: msg
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(iso_solver_t) :: iso
 type(edos_t) :: edos
!arrays
 real(dp),allocatable :: ktmesh(:), lambda_ij(:), imag_w(:), imag_2w(:), phmesh(:), a2fw(:)
!----------------------------------------------------------------------

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)
 units = [std_out, ab_out]

 call wrtout(std_out, " Solving isotropic Migdal-Eliashberg equations on the imaginary axis", pre_newlines=2)
 call cwtime(cpu, wall, gflops, "start")

 cryst => gstore%cryst; ebands => gstore%ebands
 !natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "bz", "qzone == 'bz' is required", ierr)
 ABI_CHECK(ierr == 0, "Wrong gstore object for migdal_eliashberg_iso. See messages above")

 ! Compute electron DOS.
 call dtset%get_edos_params(edos_intmeth, edos_step, edos_broad)
 edos = ebands%get_edos(cryst, edos_intmeth, edos_step, edos_broad, gstore%comm)

 !! Store DOS per spin channel
 !n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos%print(units)
   !path = strcat(dtfil%filnam_ds(4), "_EDOS")
   !call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path, ch10))
   !call edos%write(path)
 end if

 ! Compute phonon frequency mesh.
 call gstore%ifc%get_phmesh(dtset%ph_wstep, phmesh_size, phmesh)

 ! Compute Eliashberg function a2F(w)
 ABI_MALLOC(a2fw, (phmesh_size))
 call get_a2fw(gstore, dtset, phmesh_size, phmesh, a2fw)

 ncid = nctk_noid
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_ISOME.nc") , xmpi_comm_self))
   !write(777, *)"# phmesh (meV), a2fw"
   !do iw=1, phmesh_size
   !  write(777, *) phmesh(iw) * Ha_meV, a2fw(iw) / (edos%gef(0) / two)
   !end do
   !close(777)
 end if

 ABI_FREE(a2fw)
 ABI_FREE(phmesh)
 call edos%free()

 call dtset%get_ktmesh(ntemp, ktmesh)

 !NVHPC and LLVM don't like using this constructor because allocatable arrays aren't set.
#if defined FC_NVHPC || defined FC_LLVM
  iso%ntemp=ntemp
  iso%max_niter=10
  iso%tolerance=tol10
  iso%ncid=ncid
  iso%comm=gstore%comm
#else
 iso = iso_solver_t(ntemp=ntemp, max_niter=10, tolerance=tol10, ncid=ncid, comm=gstore%comm)
#endif

 do itemp=1,ntemp
   ! Generate Matsubara mesh for this T with cutoff wmax.
   kt = ktmesh(itemp)
   wmax = one
   call matsubara_mesh("bosons", kt, wmax, niw, imag_w)

   ! Compute lambda(w_i - w_j)
   ABI_MALLOC(lambda_ij, (2 * niw))
   ABI_MALLOC(imag_2w, (2 * niw))

   !call wrtout(std_out, " Computing lambda_iso_iw...")
   !call get_lambda_iso_iw(gstore, 2 * niw, imag_2w, lambda_ij)
   ABI_FREE(imag_2w)

   !call iso%solve(itemp, kt, niw, imag_w, lambda_ij)
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
!!***

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
!! SOURCE

subroutine matsubara_mesh(bosons_or_fermions, kt, wmax, niw, imag_w)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: bosons_or_fermions
 real(dp),intent(in) :: kt, wmax
 integer,intent(out) :: niw
 real(dp),allocatable,intent(out) :: imag_w(:)

!Local variables-------------------------------
 integer :: nn
!----------------------------------------------------------------------

 select case (bosons_or_fermions)
 case ("bosons")
   ! 2 n pi kT
   niw = nint(wmax / (two * pi * kt)) + 1
   ABI_MALLOC(imag_w, (niw))
   do nn=0,niw-1
     imag_w(nn + 1) = two * nn * pi * kt
   end do

 case ("fermions")
   ! 2 (n + 1) pi kT
   niw = nint(wmax / (two * pi * kt))
   ABI_MALLOC(imag_w, (niw))
   do nn=0,niw-1
     imag_w(nn + 1) = two * (nn  + 1) * pi * kt
   end do

 case default
   ABI_ERROR(sjoin("Wrong values for bosons_or_fermions:", bosons_or_fermions))
 end select

end subroutine matsubara_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/get_lambda_iso_iw
!! NAME
!! get_lambda_iso_iw
!!
!! FUNCTION
!!  Compute isotropic lambda along the imaginary axis
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine get_lambda_iso_iw(gstore, nw, imag_w, lambda)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 integer,intent(in) :: nw
 real(dp),intent(in) :: imag_w(nw)
 real(dp),intent(out) :: lambda(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, nb_k, nb_kq
 real(dp) :: g2, wqnu, weight_k, weight_q
!arrays
 real(dp) :: qpt(3)
 real(dp),allocatable :: dbl_delta_q(:,:,:), g2_pmnk(:,:,:,:)
!----------------------------------------------------------------------

 ABI_CHECK(gstore%qzone == "bz", "get_lambda_iso_iw assumes qzone == `bz`")
 !if (gstore%check_cplex_qkzone_gmode(cplex1, "bz", kzone, gmode, kfilter) result(ierr)

 lambda = zero
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbl_delta_q, (nb_kq, nb_k, gqk%my_nk))
   ABI_MALLOC(g2_pmnk, (gqk%my_npert, nb_kq, nb_k, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     ! Compute integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, gstore%dtset%eph_intmeth, gstore%dtset%eph_fsmear, qpt, weight_q, dbl_delta_q)

     ! Copy data to improve memory access in the loops below.
     g2_pmnk = gqk%my_g2(:,:,my_iq,:,:)

     do my_ik=1,gqk%my_nk
       weight_k = gqk%my_wtk(my_ik)
       do in_k=1,nb_k
         do im_kq=1,nb_kq
           do my_ip=1,gqk%my_npert
             g2 = g2_pmnk(my_ip, im_kq, in_k, my_ik)
             ! TODO: handle wqnu ~ 0
             wqnu = gqk%my_wnuq(my_ip, my_iq)
             lambda(:) = lambda(:) + &
               two * wqnu / (imag_w(:) ** 2 + wqnu ** 2) * g2 * weight_k * weight_q * dbl_delta_q(im_kq, in_k, my_ik)
           end do
         end do
       end do
     end do
   end do ! my_iq

   ABI_FREE(dbl_delta_q)
   ABI_FREE(g2_pmnk)
   end associate
 end do ! my_is

 ! Take into account collinear spin
 lambda = lambda * (two / (gstore%nsppol * gstore%dtset%nspinor))
 call xmpi_sum(lambda, gstore%comm, ierr)

end subroutine get_lambda_iso_iw
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/get_a2fw
!! NAME
!! get_a2fw
!!
!! FUNCTION
!!  Compute Eliashberg function a^2F(omega).
!!
!! INPUTS
!!  nw: Number of frequencies.
!!  wmesh: Frequency mesh.
!!
!! OUTPUT
!! a2fw(nw): Eliashberg function.
!!
!! SOURCE

subroutine get_a2fw(gstore, dtset, nw, wmesh, a2fw)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nw
 real(dp),intent(in) :: wmesh(nw)
 real(dp),intent(out) :: a2fw(nw)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, timrev_q, ii, ik_ibz, nb_k, nb_kq
 real(dp) :: g2_qnu, wqnu, weight_k, weight_q, cpu, wall, gflops
 type(lgroup_t) :: lg_myq
 character(len=500) :: msg !, kk_string !, qq_bz_string
!arrays
 integer :: units(2)
 real(dp) :: qpt(3), kk(3)
 real(dp),allocatable :: dbl_delta_q(:,:,:), g2_mnkp(:,:,:,:), deltaw_nuq(:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(units, sjoin(" Computing a^2F(w) with ph_smear:", ftoa(gstore%dtset%ph_smear * Ha_meV), "(meV)"), pre_newlines=1)

 !if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "phonon") /= 0) then
 !  ABI_ERROR("The gstore object is inconsistent with gstore_wannierize_and_write_gwan. See messages above.")
 !end if

 ABI_CHECK(gstore%qzone == "bz", "get_a2fw assumes qzone == `bz`")
 ! Check consistency of little group options.
 ABI_CHECK(gstore%check_little_group(dtset, msg) == 0, msg)

 ABI_MALLOC(deltaw_nuq, (nw))

 a2fw = zero

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbl_delta_q, (nb_kq, nb_k, gqk%my_nk))
   ABI_MALLOC(g2_mnkp, (nb_kq, nb_k, gqk%my_nk, gqk%my_npert))

   ! Loop over my q-points.
   do my_iq=1,gqk%my_nq
     ! Compute all integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, gstore%dtset%eph_intmeth, gstore%dtset%eph_fsmear, qpt, weight_q, dbl_delta_q)

     ! Copy data to improve memory access in the loops below.
     do my_ip=1,gqk%my_npert
       g2_mnkp(:,:,:,my_ip) = gqk%my_g2(my_ip,:,my_iq,:,:)
     end do

     ! Compute the little group of the q-point so that we only need to sum g(k,q) for k in the IBZ_q
     if (dtset%gstore_use_lgq /= 0) then
       timrev_q = kpts_timrev_from_kptopt(gstore%qptopt)
       call lg_myq%init(cryst, qpt, timrev_q, gstore%nkbz, gstore%kbz, gstore%nkibz, gstore%kibz, xmpi_comm_self)
     end if

     ! Loop over my phonon modes.
     do my_ip=1,gqk%my_npert
       wqnu = gqk%my_wnuq(my_ip, my_iq)
       ! delta(w - omega_qnu)
       deltaw_nuq = gaussian(wmesh - wqnu, gstore%dtset%ph_smear)

       ! Loop over my k-points.
       do my_ik=1,gqk%my_nk
         kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik); weight_k = gqk%my_wtk(my_ik)

         ! Handle little group and integration weight.
         if (dtset%gstore_use_lgq /= 0) then
           ii = lg_myq%findq_ibzk(kk); if (ii == -1) cycle; weight_k = lg_myq%weights(ii)
         end if

         ! Sum over m_kq and n_k and accumulate.
         do in_k=1,nb_k
           do im_kq=1,nb_kq
             g2_qnu = g2_mnkp(im_kq, in_k, my_ik, my_ip)
             a2fw(:) = a2fw(:) + deltaw_nuq(:) * g2_qnu * weight_k * weight_q * dbl_delta_q(im_kq, in_k, my_ik)
           end do
         end do
       end do
     end do

     call lg_myq%free()
   end do ! my_iq

   ABI_FREE(dbl_delta_q)
   ABI_FREE(g2_mnkp)
   end associate
 end do ! my_is

 ABI_FREE(deltaw_nuq)

 ! Take into account collinear spin and N(eF) TODO
 a2fw = a2fw * (two / (gstore%nsppol * gstore%dtset%nspinor))
 call xmpi_sum(a2fw, gstore%comm, ierr)

 call cwtime_report(" get_a2fw", cpu, wall, gflops)

end subroutine get_a2fw
!!***

end module m_migdal_eliashberg
!!***
