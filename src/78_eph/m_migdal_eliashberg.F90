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
 use m_numeric_tools,   only : simpson_int
 use m_special_funcs,   only : gaussian
 use m_ebands,          only : ebands_t, edos_t
 use m_bz_mesh,         only : kpath_t
 use m_geometry,        only : phdispl_cart2red_nmodes
 use m_ephtk,           only : ephtk_gkknu_from_atm, EPHTK_WTOL
 use m_gstore,          only : gstore_t

 implicit none

 private

 public :: migdal_eliashberg_iso
 !public :: migdal_eliashberg_aniso

 ! Internal policy for mode-resolved lambda on a q-path. When enabled, all
 ! modes in a degenerate phonon subspace are assigned the subspace-averaged
 ! lambda. The sum over the subspace is therefore preserved.
 logical,parameter :: average_lambda_qpath_degenerate = .True.
 real(dp),parameter :: lambda_qpath_degen_tol = tol6
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
 integer :: nproc, my_rank, ierr, itemp, ntemp, niw, ncid, iw, spin
 integer :: edos_intmeth
 !integer :: spin, natom3, cnt !, band, ib, nb, my_ik, my_iq, my_is
 !integer :: ik_ibz, ik_bz, ebands_timrev, iq_bz, iq_ibz !, ikq_ibz, ikq_bz
 !integer :: ncid, spin_ncid, ncerr, gstore_fform
 integer :: phmesh_size, units(2) !, iw
 real(dp) :: kt, wmax, cpu, wall, gflops, edos_step, edos_broad, lambda_iso, omega_log, omega_2, alpha !, sigma, ecut, eshift, eig0nk
 character(len=500) :: msg
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(iso_solver_t) :: iso
 type(edos_t) :: edos
!arrays
 real(dp),allocatable :: ktmesh(:), lambda_ij(:), imag_w(:), imag_2w(:), phmesh(:), a2fw(:), a2fw_raw(:)
 real(dp),allocatable :: a2f_1mom(:), a2f_1mom_int(:)
 real(dp),allocatable :: phfreq_qibz(:,:), phlambda_qibz(:,:,:)
 real(dp),allocatable :: qpath(:,:), phfreq_qpath(:,:), phdispl_cart_qpath(:,:,:,:), phlambda_qpath(:,:,:)
!----------------------------------------------------------------------

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)
 units = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(units, " Solving isotropic Migdal-Eliashberg equations on the imaginary axis", pre_newlines=2)

 cryst => gstore%cryst; ebands => gstore%ebands
 !natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "gstore_kzone == 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(any(gstore%qzone == ["bz ", "ibz"]), "gstore_qzone must be 'bz' or 'ibz'", ierr)
 ABI_CHECK(ierr == 0, "Wrong gstore object for migdal_eliashberg_iso. See messages above")

 ! Compute electron DOS.
 call dtset%get_edos_params(edos_intmeth, edos_step, edos_broad)
 edos = ebands%get_edos(cryst, edos_intmeth, edos_step, edos_broad, gstore%comm)

 ! A disentangled Wannier Hamiltonian generally spans only a subspace of the
 ! original bands. In this case ebands%nelect still describes the complete
 ! ab-initio manifold, so locating eF by integrating the DOS can fail (e.g. a
 ! four-band model with nelect = 8 appears completely filled). The electronic
 ! delta functions used below are centered at ebands%fermie, copied from the
 ! ab-initio bands, hence evaluate N(eF) at the same chemical potential.
 if (edos%ief == 0 .and. gstore%has_wannier) then
   iw = int((ebands%fermie - edos%mesh(1)) / edos%step) + 1
   ABI_CHECK(iw >= 1 .and. iw < edos%nw, "The ab-initio Fermi level lies outside the energy range of the Wannier-interpolated bands")
   alpha = (ebands%fermie - edos%mesh(iw)) / edos%step
   do spin=0,edos%nsppol
     edos%gef(spin) = (one - alpha) * edos%dos(iw,spin) + alpha * edos%dos(iw+1,spin)
     edos%ghf(spin) = edos%gef(spin)
   end do
   edos%ief = iw
   edos%ihf = iw
   call wrtout(units, " Using the ab-initio Fermi level to evaluate the DOS of the Wannier band subspace.")
 end if

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
 call get_a2fw(gstore, edos%gef(0), phmesh_size, phmesh, a2fw, phfreq_qibz, phlambda_qibz)

 ! Compute mode-resolved lambda on the phonon q-path when the matrix
 ! elements can be evaluated at arbitrary q with Wannier interpolation.
 call get_lambda_qpath_wan(gstore, dtset, edos%gef(0), qpath, phfreq_qpath, phdispl_cart_qpath, phlambda_qpath)

 ! Save results.
 ncid = nctk_noid
 if (my_rank == master) then
   call alloc_copy(a2fw, a2fw_raw)
   ABI_CHECK(edos%gef(0) > zero, "The electronic DOS at the Fermi level must be positive")
   a2fw = a2fw / (edos%gef(0) / two)

   ABI_MALLOC(a2f_1mom, (phmesh_size))
   ABI_MALLOC(a2f_1mom_int, (phmesh_size))
   a2f_1mom = zero
   where (phmesh > tol12) a2f_1mom = a2fw / phmesh
   call simpson_int(phmesh_size, dtset%ph_wstep, a2f_1mom, a2f_1mom_int)
   lambda_iso = two * a2f_1mom_int(phmesh_size)
   write(msg, "(a,es16.8)")" Isotropic lambda from a2F(w): ", lambda_iso
   call wrtout(units, msg)

   ABI_CHECK(lambda_iso > zero, "Cannot compute omega_log because the isotropic lambda is not positive")
   a2f_1mom = zero
   do iw=1,phmesh_size
     if (phmesh(iw) > tol12) a2f_1mom(iw) = a2fw(iw) * log(phmesh(iw)) / phmesh(iw)
   end do
   call simpson_int(phmesh_size, dtset%ph_wstep, a2f_1mom, a2f_1mom_int)
   omega_log = exp(two * a2f_1mom_int(phmesh_size) / lambda_iso)
   write(msg, "(a,es16.8,a,es16.8,a)")" Isotropic omega_log from a2F(w): ", omega_log, &
     " (Ha), ", omega_log * Ha_K, " (K)"
   call wrtout(units, msg)

   ! Allen-Dynes square-root second moment:
   ! omega_2^2 = (2 / lambda) integral dw w a2F(w).
   a2f_1mom = phmesh * a2fw
   call simpson_int(phmesh_size, dtset%ph_wstep, a2f_1mom, a2f_1mom_int)
   omega_2 = sqrt(two * a2f_1mom_int(phmesh_size) / lambda_iso)
   write(msg, "(a,es16.8,a,es16.8,a)")" Isotropic omega_2 from a2F(w): ", omega_2, &
     " (Ha), ", omega_2 * Ha_K, " (K)"
   call wrtout(units, msg)
   ABI_FREE(a2f_1mom)
   ABI_FREE(a2f_1mom_int)

   NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), "_ISOME.nc") , xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands%ncwrite(ncid))
   NCF_CHECK(edos%ncwrite(ncid))
   call isome_ncwrite_spectral(ncid, dtset, gstore, phmesh_size, phmesh, a2fw_raw, a2fw, edos%gef(0), omega_2)
   call isome_ncwrite_qibz(ncid, gstore, phfreq_qibz, phlambda_qibz)
   if (allocated(qpath)) then
     call isome_ncwrite_qpath(ncid, qpath, phfreq_qpath, phdispl_cart_qpath, phlambda_qpath)
   end if
   ABI_FREE(a2fw_raw)
 end if

 ABI_SFREE(qpath)
 ABI_SFREE(phfreq_qpath)
 ABI_SFREE(phdispl_cart_qpath)
 ABI_SFREE(phlambda_qpath)
 ABI_FREE(phfreq_qibz)
 ABI_FREE(phlambda_qibz)
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

 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if

 ABI_FREE(ktmesh)
 call iso%free()

 call cwtime_report(" migdal_eliashberg_iso:", cpu, wall, gflops)

end subroutine migdal_eliashberg_iso
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/isome_ncwrite_spectral
!! NAME
!! isome_ncwrite_spectral
!!
!! FUNCTION
!!  Write the first revision of the spectral section of the ISOME.nc file.
!!  The raw matrix-element spectral sum is stored together with the
!!  provisionally DOS-normalized Eliashberg function so that the normalization
!!  can be audited without recomputing the electron-phonon matrix elements.
!!
!! INPUTS
!!  ncid=NetCDF file identifier, opened on xmpi_comm_self by the master rank.
!!  dtset<dataset_type>=Input variables.
!!  gstore<gstore_t>=Electron-phonon matrix-element container.
!!  nomega=Number of points in the phonon-frequency mesh.
!!  omega(nomega)=Phonon-frequency mesh in Hartree.
!!  a2f_raw(nomega)=Raw matrix-element spectral sum.
!!  a2f(nomega)=DOS-normalized Eliashberg function.
!!  edos_fermie=Total electronic DOS at the Fermi level.
!!  omega_2=Allen-Dynes square-root second moment in Hartree.
!!
!! SOURCE

subroutine isome_ncwrite_spectral(ncid, dtset, gstore, nomega, omega, a2f_raw, a2f, edos_fermie, omega_2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid, nomega
 real(dp),intent(in) :: edos_fermie, omega_2
 type(dataset_type),intent(in) :: dtset
 type(gstore_t),intent(in) :: gstore
!arrays
 real(dp),intent(in) :: omega(nomega), a2f_raw(nomega), a2f(nomega)

!Local variables-------------------------------
!scalars
 integer :: ncerr
!----------------------------------------------------------------------

 ncerr = nctk_def_dims(ncid, nctkdim_t("a2f_nomega", nomega), defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   "isome_schema_version", "eph_intmeth", "ph_intmeth"])
 NCF_CHECK(ncerr)
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   "eph_fsmear", "ph_smear", "ph_wstep", "a2f_edos_fermie", "a2f_dos_normalization", "omega_2"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("gstore_ngqpt", "int", "three"), &
   nctkarr_t("eph_ngqpt_fine", "int", "three"), &
   nctkarr_t("ddb_ngqpt", "int", "three"), &
   nctkarr_t("a2f_mesh", "dp", "a2f_nomega"), &
   nctkarr_t("a2f_values_raw", "dp", "a2f_nomega"), &
   nctkarr_t("a2f_values", "dp", "a2f_nomega")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_atomic_units(ncid, "a2f_mesh"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "eph_fsmear"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "ph_smear"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "ph_wstep"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "a2f_edos_fermie"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "a2f_dos_normalization"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "omega_2"))

 ncerr = nf90_put_att(ncid, nf90_global, "isome_section", "spectral")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "omega_2"), "long_name", &
   "Allen-Dynes square-root second moment of the isotropic Eliashberg function")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nf90_global, "isome_status", "experimental")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "a2f_values_raw"), "long_name", &
   "raw matrix-element spectral sum before division by the electronic DOS")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "a2f_values"), "long_name", &
   "isotropic Eliashberg spectral function")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "a2f_values"), "normalization", &
   "a2f_values_raw divided by a2f_edos_fermie / 2")
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 ncerr = nctk_write_iscalars(ncid, &
   [character(len=nctk_slen) :: "isome_schema_version", "eph_intmeth", "ph_intmeth"], &
   [3, dtset%eph_intmeth, dtset%ph_intmeth])
 NCF_CHECK(ncerr)
 ncerr = nctk_write_dpscalars(ncid, &
   [character(len=nctk_slen) :: &
     "eph_fsmear", "ph_smear", "ph_wstep", "a2f_edos_fermie", "a2f_dos_normalization", "omega_2"], &
   [dtset%eph_fsmear, dtset%ph_smear, dtset%ph_wstep, edos_fermie, edos_fermie / two, omega_2])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gstore_ngqpt"), gstore%ngqpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_ngqpt_fine"), dtset%eph_ngqpt_fine))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ddb_ngqpt"), dtset%ddb_ngqpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a2f_mesh"), omega))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a2f_values_raw"), a2f_raw))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a2f_values"), a2f))

end subroutine isome_ncwrite_spectral
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/isome_ncwrite_qibz
!! NAME
!! isome_ncwrite_qibz
!!
!! FUNCTION
!!  Write phonon frequencies and the raw mode-resolved lambda(q,nu) on
!!  the phonon IBZ. No averaging is applied inside degenerate subspaces;
!!  their gauge-invariant contribution enters a2F through the mode sum.
!!
!! SOURCE

subroutine isome_ncwrite_qibz(ncid, gstore, phfreq, phlambda)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 type(gstore_t),intent(in) :: gstore
 real(dp),intent(in) :: phfreq(:,:), phlambda(:,:,:)

!Local variables-------------------------------
 integer :: ncerr, natom3
!----------------------------------------------------------------------

 natom3 = size(phfreq, dim=1)
 ABI_CHECK_IEQ(size(phfreq, dim=2), gstore%nqibz, "Invalid phfreq q-IBZ dimension")
 ABI_CHECK_IEQ(size(phlambda, dim=1), natom3, "Invalid phlambda mode dimension")
 ABI_CHECK_IEQ(size(phlambda, dim=2), gstore%nqibz, "Invalid phlambda q-IBZ dimension")
 ABI_CHECK_IEQ(size(phlambda, dim=3), gstore%nsppol, "Invalid phlambda spin dimension")

 ncerr = nctk_def_dims(ncid, [nctkdim_t("isome_nqibz", gstore%nqibz), &
   nctkdim_t("isome_natom3", natom3)], defmode=.True.)
 NCF_CHECK(ncerr)
 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("qibz", "dp", "number_of_reduced_dimensions, isome_nqibz"), &
   nctkarr_t("wtq", "dp", "isome_nqibz"), &
   nctkarr_t("phfreq_qibz", "dp", "isome_natom3, isome_nqibz"), &
   nctkarr_t("phlambda_qibz", "dp", "isome_natom3, isome_nqibz, number_of_spins")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_atomic_units(ncid, "phfreq_qibz"))
 NCF_CHECK(nf90_put_att(ncid, nctk_idname(ncid, "qibz"), "long_name", "phonon q-points in the irreducible Brillouin zone"))
 NCF_CHECK(nf90_put_att(ncid, nctk_idname(ncid, "wtq"), "long_name", "phonon IBZ integration weights"))
 NCF_CHECK(nf90_put_att(ncid, nctk_idname(ncid, "phlambda_qibz"), "long_name", "raw mode-resolved electron-phonon coupling lambda(q,nu)"))
 NCF_CHECK(nf90_put_att(ncid, nctk_idname(ncid, "phlambda_qibz"), "degenerate_mode_averaging", "none"))

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qibz"), gstore%qibz))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wtq"), gstore%wtq))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreq_qibz"), phfreq))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phlambda_qibz"), phlambda))

end subroutine isome_ncwrite_qibz
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/isome_ncwrite_qpath
!! NAME
!! isome_ncwrite_qpath
!!
!! FUNCTION
!!  Write phonon frequencies, displacements, and mode-resolved
!!  electron-phonon coupling lambda(q,nu) along the input phonon path.
!!
!! SOURCE

subroutine isome_ncwrite_qpath(ncid, qpath, phfreq, phdispl_cart, phlambda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
!arrays
 real(dp),intent(in) :: qpath(:,:), phfreq(:,:), phdispl_cart(:,:,:,:), phlambda(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ncerr, nqpath, natom3, nsppol
!----------------------------------------------------------------------

 nqpath = size(qpath, dim=2)
 natom3 = size(phfreq, dim=1)
 nsppol = size(phlambda, dim=3)

 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("isome_nqpath", nqpath), &
   nctkdim_t("number_of_phonon_modes", natom3), &
   nctkdim_t("isome_nsppol", nsppol)], defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "phlambda_qpath_average_degenerate"])
 NCF_CHECK(ncerr)
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "phlambda_qpath_degen_tol"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("qpath", "dp", "number_of_reduced_dimensions, isome_nqpath"), &
   nctkarr_t("phfreq_qpath", "dp", "number_of_phonon_modes, isome_nqpath"), &
   nctkarr_t("phdispl_cart_qpath", "dp", "two, number_of_phonon_modes, number_of_phonon_modes, isome_nqpath"), &
   nctkarr_t("phlambda_qpath", "dp", "number_of_phonon_modes, isome_nqpath, isome_nsppol")])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_atomic_units(ncid, "phfreq_qpath"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "phlambda_qpath_degen_tol"))
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "qpath"), "long_name", &
   "q-point path in reduced coordinates")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "phlambda_qpath"), "long_name", &
   "mode-resolved electron-phonon coupling lambda(q,nu)")
 NCF_CHECK(ncerr)
 ncerr = nf90_put_att(ncid, nctk_idname(ncid, "phlambda_qpath"), "fermi_surface_integration", &
   "Gaussian double delta with width eph_fsmear")
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nctk_write_iscalars(ncid, [character(len=nctk_slen) :: "phlambda_qpath_average_degenerate"], [merge(1, 0, average_lambda_qpath_degenerate)]))
 NCF_CHECK(nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "phlambda_qpath_degen_tol"], [lambda_qpath_degen_tol]))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpath"), qpath))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreq_qpath"), phfreq))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phdispl_cart_qpath"), phdispl_cart))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phlambda_qpath"), phlambda))

end subroutine isome_ncwrite_qpath
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/get_lambda_qpath_wan
!! NAME
!! get_lambda_qpath_wan
!!
!! FUNCTION
!!  Compute lambda(q,nu) on the path defined by ph_qpath and ph_ndivsm.
!!  The path construction is the same as in ifc_mkphbs. Arbitrary-q
!!  electronic energies and matrix elements are obtained from the Wannier
!!  Hamiltonian and GWAN interpolation, respectively.
!!
!! SOURCE

subroutine get_lambda_qpath_wan(gstore, dtset, edos_fermie, qpoints, phfreq, phdispl_cart, phlambda)

!Arguments ------------------------------------
!scalars
 class(gstore_t),intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 real(dp),intent(in) :: edos_fermie
!arrays
 real(dp),allocatable,intent(out) :: qpoints(:,:), phfreq(:,:), phdispl_cart(:,:,:,:), phlambda(:,:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: my_is, my_ik, my_ip, iq, ipc, nu, in_k, im_kq, ierr, my_rank
 integer :: natom, natom3, nwan, nqpath, spin, ik_start, ik_stop, ikb, nkb, nk_batch
 real(dp) :: weight_k, g2, fs_weight, spin_factor, cpu, wall, gflops
 logical :: has_gwan
 type(kpath_t) :: qpath
 character(len=500) :: msg
!arrays
 integer :: units(2)
 real(dp),allocatable :: displ_cart4(:,:,:,:), displ_red4(:,:,:,:), displ_red(:,:,:)
 real(dp),allocatable :: gatm_real(:,:,:,:,:), gnu_real(:,:,:,:,:)
 real(dp),allocatable :: eig_k(:,:), eig_kq(:,:)
 complex(dp),allocatable :: intp_gatm(:,:,:,:), gatm_full(:,:,:,:), g_req(:,:,:,:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]

 if (dtset%ph_nqpath <= 0 .or. dtset%ph_ndivsm <= 0) return

 has_gwan = gstore%my_nspins > 0
 if (has_gwan) has_gwan = allocated(gstore%gqk(1)%wan%grpe_wwp)
 if (.not. has_gwan) then
   call wrtout([std_out, ab_out], &
     " Skipping lambda(q,nu): GWAN interpolation is not available for this GSTORE.")
   return
 end if

 ABI_CHECK(dtset%eph_fsmear > zero, "lambda(q,nu) along a path requires a positive eph_fsmear for Gaussian Fermi-surface integration")
 ABI_CHECK(edos_fermie > zero, "The electronic DOS at the Fermi level must be positive")
 ABI_CHECK(gstore%kzone == "bz", "get_lambda_qpath_wan requires gstore_kzone = 'bz'")
 call cwtime(cpu, wall, gflops, "start")

 natom = gstore%cryst%natom; natom3 = 3 * natom
 call qpath%init(dtset%ph_qpath(:,1:dtset%ph_nqpath), gstore%cryst%gprimd, dtset%ph_ndivsm)
 nqpath = qpath%npts

 ABI_MALLOC(qpoints, (3, nqpath))
 ABI_MALLOC(phfreq, (natom3, nqpath))
 ABI_MALLOC(phdispl_cart, (2, natom3, natom3, nqpath))
 ABI_CALLOC(phlambda, (natom3, nqpath, gstore%nsppol))
 ABI_MALLOC(displ_cart4, (2, 3, natom, natom3))
 ABI_MALLOC(displ_red4, (2, 3, natom, natom3))
 ABI_MALLOC(displ_red, (2, natom3, natom3))

 qpoints = qpath%points
 do iq=1,nqpath
   call gstore%ifc%fourq(gstore%cryst, qpoints(:,iq), phfreq(:,iq), displ_cart4, out_displ_red=displ_red4)
   phdispl_cart(:,:,:,iq) = reshape(displ_cart4, [2, natom3, natom3])
 end do

 write(msg, "(a,i0,a,es12.4,a)")" Computing lambda(q,nu) at ", nqpath, &
   " path points with eph_fsmear: ", dtset%eph_fsmear * Ha_meV, " meV"
 call wrtout(units, msg, pre_newlines=1)

 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is))
   spin = gqk%spin; nwan = gqk%wan%nwan
   ABI_CHECK_IEQ(nwan, gqk%nb_k, "Wannier and GSTORE band dimensions differ")
   ABI_CHECK_IEQ(nwan, gqk%nb_kq, "Wannier and GSTORE band dimensions differ")

   nk_batch = gqk%wan%eph_kbatch_size(gqk%my_nk, natom3)
   ABI_MALLOC(intp_gatm, (nwan, nwan, gqk%my_npert, nk_batch))
   ABI_MALLOC(gatm_full, (nwan, nwan, natom3, nk_batch))
   ABI_MALLOC(g_req, (gqk%wan%nr_e, nwan, nwan, gqk%my_npert))
   ABI_MALLOC(gatm_real, (2, nwan, nwan, 1, natom3))
   ABI_MALLOC(gnu_real, (2, nwan, nwan, 1, natom3))
   ABI_MALLOC(eig_k, (nwan, nk_batch))
   ABI_MALLOC(eig_kq, (nwan, nk_batch))

   do iq=1,nqpath
     if (gqk%qpt_comm%skip(iq)) cycle
     call gstore%ifc%fourq(gstore%cryst, qpoints(:,iq), phfreq(:,iq), displ_cart4, out_displ_red=displ_red4)
     displ_red = reshape(displ_red4, [2, natom3, natom3])
     call gqk%wan%prepare_eph_q(qpoints(:,iq), g_req)

     do ik_start=1,gqk%my_nk,nk_batch
       ik_stop = min(gqk%my_nk, ik_start + nk_batch - 1); nkb = ik_stop - ik_start + 1

       ! Interpolate one bounded k block and reuse the returned electronic
       ! energies in the Fermi-surface weights.
       call gqk%wan%interp_eph_manyk_from_q(gstore%cryst, nkb, gqk%my_kpts(:,ik_start:ik_stop), &
                                            qpoints(:,iq), g_req, intp_gatm(:,:,:,1:nkb), &
                                            out_eigens_k=eig_k(:,1:nkb), out_eigens_kq=eig_kq(:,1:nkb))

       gatm_full(:,:,:,1:nkb) = czero
       do ikb=1,nkb
         do ipc=1,gqk%my_npert
           gatm_full(:,:,gqk%my_pertcases(ipc),ikb) = intp_gatm(:,:,ipc,ikb)
         end do
       end do
       if (gqk%pert_comm%nproc > 1) call xmpi_sum(gatm_full(:,:,:,1:nkb), gqk%pert_comm%value, ierr)

       do ikb=1,nkb
         my_ik = ik_start + ikb - 1
         weight_k = gqk%my_wtk(my_ik)
         gatm_real(1,:,:,1,:) = real(gatm_full(:,:,:,ikb), kind=dp)
         gatm_real(2,:,:,1,:) = aimag(gatm_full(:,:,:,ikb))
         call ephtk_gkknu_from_atm(nwan, nwan, 1, natom, gatm_real, phfreq(:,iq), displ_red, gnu_real)

         do my_ip=1,gqk%my_npert
           nu = gqk%my_pertcases(my_ip)
           if (phfreq(nu,iq) < EPHTK_WTOL) cycle
           do in_k=1,nwan
             do im_kq=1,nwan
               g2 = gnu_real(1,im_kq,in_k,1,nu)**2 + gnu_real(2,im_kq,in_k,1,nu)**2
               fs_weight = gaussian(eig_k(in_k,ikb) - gstore%ebands%fermie, dtset%eph_fsmear) * &
                           gaussian(eig_kq(im_kq,ikb) - gstore%ebands%fermie, dtset%eph_fsmear)
               phlambda(nu,iq,spin) = phlambda(nu,iq,spin) + two * g2 * weight_k * fs_weight / phfreq(nu,iq)
             end do
           end do
         end do
       end do ! ikb
     end do ! ik_start
   end do

   ABI_FREE(intp_gatm)
   ABI_FREE(gatm_full)
   ABI_FREE(g_req)
   ABI_FREE(gatm_real)
   ABI_FREE(gnu_real)
   ABI_FREE(eig_k)
   ABI_FREE(eig_kq)
   end associate
 end do

 spin_factor = two / (gstore%nsppol * dtset%nspinor)
 phlambda = phlambda * spin_factor / (edos_fermie / two)
 call xmpi_sum(phlambda, gstore%comm, ierr)
 if (average_lambda_qpath_degenerate) call average_lambda_degenerate_modes(nqpath, natom3, gstore%nsppol, phfreq, phlambda)
 my_rank = xmpi_comm_rank(gstore%comm)
 if (my_rank == master) then
   call wrtout(ab_out, " Phonon frequencies and lambda(q,nu) along the q-path:", pre_newlines=1)
   do iq=1,nqpath
     write(msg, "(a,i0,a,3es16.8)")" q-path point ", iq, ": ", qpoints(:,iq)
     call wrtout(ab_out, msg)
     select case (gstore%nsppol)
     case (1)
       call wrtout(ab_out, "   nu       omega (meV)          lambda")
       do nu=1,natom3
         write(msg, "(i5,2x,es16.8,2x,es16.8)")nu, phfreq(nu,iq) * Ha_meV, phlambda(nu,iq,1)
         call wrtout(ab_out, msg)
       end do
     case (2)
       call wrtout(ab_out, "   nu       omega (meV)         lambda_spin1        lambda_spin2")
       do nu=1,natom3
         write(msg, "(i5,2x,es16.8,2x,es16.8,2x,es16.8)")nu, phfreq(nu,iq) * Ha_meV, phlambda(nu,iq,1:2)
         call wrtout(ab_out, msg)
       end do
     case default
       ABI_ERROR("Printing lambda(q,nu) supports only nsppol = 1 or 2")
     end select
   end do
 end if

 ABI_FREE(displ_cart4)
 ABI_FREE(displ_red4)
 ABI_FREE(displ_red)
 call qpath%free()
 call cwtime_report(" Wannier interpolation of lambda(q,nu) along q-path", cpu, wall, gflops)

end subroutine get_lambda_qpath_wan
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/average_lambda_degenerate_modes
!! NAME
!! average_lambda_degenerate_modes
!!
!! FUNCTION
!!  Replace lambda(q,nu) inside each degenerate phonon subspace with its
!!  arithmetic average. Frequencies are assumed to be ordered by branch,
!!  as returned by ifc%fourq. The sum over every subspace is preserved.
!!
!! SOURCE

subroutine average_lambda_degenerate_modes(nqpath, nmode, nsppol, phfreq, phlambda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpath, nmode, nsppol
!arrays
 real(dp),intent(in) :: phfreq(nmode,nqpath)
 real(dp),intent(inout) :: phlambda(nmode,nqpath,nsppol)

!Local variables-------------------------------
!scalars
 integer :: iq, spin, first_mode, last_mode, ndeg
 real(dp) :: lambda_avg
!----------------------------------------------------------------------

 do iq=1,nqpath
   first_mode = 1
   do while (first_mode <= nmode)
     last_mode = first_mode
     do while (last_mode < nmode)
       if (abs(phfreq(last_mode + 1,iq) - phfreq(first_mode,iq)) > lambda_qpath_degen_tol) exit
       last_mode = last_mode + 1
     end do

     ndeg = last_mode - first_mode + 1
     if (ndeg > 1) then
       do spin=1,nsppol
         lambda_avg = sum(phlambda(first_mode:last_mode,iq,spin)) / ndeg
         phlambda(first_mode:last_mode,iq,spin) = lambda_avg
       end do
     end if
     first_mode = last_mode + 1
   end do
 end do

end subroutine average_lambda_degenerate_modes
!!***

!----------------------------------------------------------------------

!!****f* m_migdal_eliashberg/get_a2fw
!! NAME
!! get_a2fw
!!
!! FUNCTION
!!  Compute the raw, mode-resolved lambda(q,nu) on the phonon IBZ and
!!  construct a^2F(omega) from its weighted sum. Phonon modes are not
!!  averaged when they are degenerate.
!!
!! INPUTS
!!  gstore: Electron-phonon matrix elements. Electronic k-points must cover
!!    the full BZ; phonon q-points may cover either the BZ or the IBZ.
!!  edos_fermie: Electronic DOS at the Fermi level.
!!  nw: Number of frequencies.
!!  wmesh: Frequency mesh.
!!
!! OUTPUT
!!  a2fw(nw): Eliashberg function in the historical raw normalization.
!!  phfreq_qibz: Phonon frequencies on the q-point IBZ.
!!  phlambda_qibz: Raw mode- and spin-resolved lambda(q,nu).
!!
!! SOURCE

subroutine get_a2fw(gstore, edos_fermie, nw, wmesh, a2fw, phfreq_qibz, phlambda_qibz)

!Arguments ------------------------------------
 class(gstore_t),intent(inout) :: gstore
 integer,intent(in) :: nw
 real(dp),intent(in) :: edos_fermie
 real(dp),intent(in) :: wmesh(nw)
 real(dp),intent(out) :: a2fw(nw)
 real(dp),allocatable,intent(out) :: phfreq_qibz(:,:), phlambda_qibz(:,:,:)

!Local variables-------------------------------
 integer :: my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, iq_ibz, isym_q, trev_q, nb_k, nb_kq
 integer :: natom, natom3, g0_q(3), nk_batch, ik_start, ik_stop, ikb, nkb, ipc, nwan
 real(dp) :: g2_qnu, wqnu, weight_k, weight_q, cpu, wall, gflops, spin_factor
 logical :: isirr_q
!arrays
 integer :: units(2)
 real(dp) :: qpt(3)
 real(dp),allocatable :: dbl_delta_q(:,:,:), g2_mnkp(:,:,:,:), deltaw_nuq(:), displ_cart_dum(:,:,:,:)
 real(dp),allocatable :: displ_red_local(:,:,:)
 complex(dp),allocatable :: intp_gatm(:,:,:,:), gatm_full(:,:,:,:), g_req(:,:,:,:), gnu(:,:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")
 call wrtout(units, sjoin(" Computing a^2F(w) with ph_smear:", ftoa(gstore%dtset%ph_smear * Ha_meV), "(meV)"), pre_newlines=1)

 ABI_CHECK(gstore%kzone == "bz", "get_a2fw requires kzone == `bz`")
 ABI_CHECK(any(gstore%qzone == ["bz ", "ibz"]), "get_a2fw requires qzone == `bz` or `ibz`")
 ABI_CHECK(any(gstore%with_cplex == [0, 1]), "get_a2fw requires with_cplex=0 or 1")
 if (gstore%with_cplex == 0) then
   ABI_CHECK(gstore%has_wannier, "get_a2fw with_cplex=0 requires a Wannier interpolator")
 end if
 ABI_CHECK(edos_fermie > zero, "The electronic DOS at the Fermi level must be positive")

 ABI_MALLOC(deltaw_nuq, (nw))
 natom = gstore%cryst%natom
 natom3 = 3 * natom
 ABI_MALLOC(phfreq_qibz, (natom3, gstore%nqibz))
 ABI_CALLOC(phlambda_qibz, (natom3, gstore%nqibz, gstore%nsppol))
 ABI_MALLOC(displ_cart_dum, (2, 3, gstore%cryst%natom, natom3))
 do iq_ibz=1,gstore%nqibz
   call gstore%ifc%fourq(gstore%cryst, gstore%qibz(:,iq_ibz), phfreq_qibz(:,iq_ibz), displ_cart_dum)
 end do
 ABI_FREE(displ_cart_dum)
 a2fw = zero

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   if (gstore%with_cplex == 1) then
     ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   end if
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   ABI_CHECK_IEQ(nb_k, nb_kq, "gqk_dbldelta_qpt does not support nb_k != nb_kq")

   ! Weights for delta(e_{m k+q}) delta(e_{n k}) for my list of k-points.
   ABI_MALLOC(dbl_delta_q, (nb_kq, nb_k, gqk%my_nk))
   ABI_MALLOC(g2_mnkp, (nb_kq, nb_k, gqk%my_nk, gqk%my_npert))

   if (gstore%with_cplex == 0) then
     nwan = gqk%wan%nwan
     ABI_CHECK_IEQ(nwan, nb_k, "Wannier nwan must agree with the gstore band range")
     nk_batch = gqk%wan%eph_kbatch_size(gqk%my_nk, natom3)
     ABI_MALLOC(intp_gatm, (nwan, nwan, gqk%my_npert, nk_batch))
     ABI_MALLOC(gatm_full, (nwan, nwan, natom3, nk_batch))
     ABI_MALLOC(g_req, (gqk%wan%nr_e, nwan, nwan, gqk%my_npert))
     ABI_MALLOC(gnu, (nwan, nwan))
     ABI_MALLOC(displ_red_local, (2, natom3, gqk%my_npert))
   end if

   ! Loop over my q-points.
   do my_iq=1,gqk%my_nq
     iq_ibz = gqk%my_q2ibz(1,my_iq)
     isym_q = gqk%my_q2ibz(2,my_iq)
     trev_q = gqk%my_q2ibz(6,my_iq)
     g0_q = gqk%my_q2ibz(3:5,my_iq)
     isirr_q = isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0)
     if (gstore%qzone == "bz" .and. .not. isirr_q) cycle

     ! Compute all integration weights for the double delta.
     call gqk%dbldelta_qpt(my_iq, gstore, gstore%dtset%eph_intmeth, gstore%dtset%eph_fsmear, qpt, weight_q, dbl_delta_q)

     if (gstore%with_cplex == 1) then
       ! Copy data to improve memory access in the loops below.
       do my_ip=1,gqk%my_npert
         g2_mnkp(:,:,:,my_ip) = gqk%my_g2(my_ip,:,my_iq,:,:)
       end do
     else
       ! Convert the symmetry-consistent phonon eigenvectors already stored
       ! in gstore to reduced coordinates for the locally owned modes.
       do my_ip=1,gqk%my_npert
         call phdispl_cart2red_nmodes(natom, 1, cryst%gprimd, &
                                      gqk%my_displ_cart(:,:,:,my_ip,my_iq), &
                                      displ_red_local(:,:,my_ip:my_ip))
       end do

       ! Interpolate bounded k blocks and immediately form |g_mnnu(k,q)|^2.
       call gqk%wan%prepare_eph_q(qpt, g_req)
       do ik_start=1,gqk%my_nk,nk_batch
         ik_stop = min(gqk%my_nk, ik_start + nk_batch - 1)
         nkb = ik_stop - ik_start + 1
         call gqk%wan%interp_eph_manyk_from_q(cryst, nkb, gqk%my_kpts(:,ik_start:ik_stop), qpt, &
                                              g_req, intp_gatm(:,:,:,1:nkb))

         gatm_full(:,:,:,1:nkb) = czero
         do ikb=1,nkb
           do ipc=1,gqk%my_npert
             gatm_full(:,:,gqk%my_pertcases(ipc),ikb) = intp_gatm(:,:,ipc,ikb)
           end do
         end do
         if (gqk%pert_comm%nproc > 1) call xmpi_sum(gatm_full(:,:,:,1:nkb), gqk%pert_comm%value, ierr)

         do ikb=1,nkb
           my_ik = ik_start + ikb - 1
           do my_ip=1,gqk%my_npert
             wqnu = gqk%my_wnuq(my_ip,my_iq)
             if (wqnu < EPHTK_WTOL) then
               g2_mnkp(:,:,my_ik,my_ip) = zero
               cycle
             end if
             gnu = czero
             do ipc=1,natom3
               gnu = gnu + gatm_full(:,:,ipc,ikb) * &
                 (displ_red_local(1,ipc,my_ip) + j_dpc * displ_red_local(2,ipc,my_ip))
             end do
             gnu = gnu / sqrt(two * wqnu)
             g2_mnkp(:,:,my_ik,my_ip) = real(gnu * conjg(gnu), kind=dp)
           end do
         end do
       end do
     end if

     ! Loop over my phonon modes.
     do my_ip=1,gqk%my_npert
       wqnu = gqk%my_wnuq(my_ip, my_iq)
       if (wqnu < EPHTK_WTOL) cycle
       ! Loop over my k-points.
       do my_ik=1,gqk%my_nk
         weight_k = gqk%my_wtk(my_ik)

         ! Sum over m_kq and n_k and accumulate.
         do in_k=1,nb_k
           do im_kq=1,nb_kq
             g2_qnu = g2_mnkp(im_kq, in_k, my_ik, my_ip)
             phlambda_qibz(gqk%my_pertcases(my_ip),iq_ibz,gqk%spin) = &
               phlambda_qibz(gqk%my_pertcases(my_ip),iq_ibz,gqk%spin) + two * g2_qnu * weight_k * &
               dbl_delta_q(im_kq, in_k, my_ik) / wqnu
           end do
         end do
       end do
     end do
   end do ! my_iq

   ABI_FREE(dbl_delta_q)
   ABI_FREE(g2_mnkp)
   if (gstore%with_cplex == 0) then
     ABI_FREE(intp_gatm)
     ABI_FREE(gatm_full)
     ABI_FREE(g_req)
     ABI_FREE(gnu)
     ABI_FREE(displ_red_local)
   end if
   end associate
 end do ! my_is

 ABI_FREE(deltaw_nuq)

 ! Normalize lambda(q,nu), then construct the normalized Eliashberg function
 ! with the phonon-IBZ weights. The factor edos_fermie / 2 at the end restores
 ! the historical raw a2F convention expected by the caller.
 spin_factor = two / (gstore%nsppol * gstore%dtset%nspinor)
 phlambda_qibz = phlambda_qibz * spin_factor / (edos_fermie / two)
 call xmpi_sum(phlambda_qibz, gstore%comm, ierr)
 a2fw = zero
 do iq_ibz=1,gstore%nqibz
   do my_ip=1,natom3
     wqnu = phfreq_qibz(my_ip,iq_ibz)
     if (wqnu < EPHTK_WTOL) cycle
     deltaw_nuq = gaussian(wmesh - wqnu, gstore%dtset%ph_smear)
     a2fw = a2fw + half * gstore%wtq(iq_ibz) * wqnu * sum(phlambda_qibz(my_ip,iq_ibz,:)) * deltaw_nuq
   end do
 end do
 a2fw = a2fw * (edos_fermie / two)

 call cwtime_report(" get_a2fw", cpu, wall, gflops)

end subroutine get_a2fw
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
!! NOTES
!!  This routine currently requires with_cplex=1, i.e. squared e-ph matrix
!!  elements precomputed and stored in gqk%my_g2. The with_cplex=0 on-demand
!!  Wannier backend implemented in get_a2fw is not yet supported here.
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
 ABI_CHECK(gstore%with_cplex == 1, "get_lambda_iso_iw requires squared e-ph matrix elements in memory (with_cplex=1)")
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

end module m_migdal_eliashberg
!!***
