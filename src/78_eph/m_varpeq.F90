!!****m* ABINIT/m_varpeq
!! NAME
!!  m_varpeq
!!
!! FUNCTION
!!  Description
!!
!! COPYRIGHT
!!  Copyright (C) 2023-2024 ABINIT group (VV, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_varpeq

 use defs_basis
 use m_abicore
 use m_dtset
 use m_dtfil
 use m_crystal
 use m_ebands
 use m_errors
 use m_krank
 use netcdf
 use m_nctk
 use m_xmpi
 use m_ifc
 use m_wfd

 use defs_datatypes,    only : ebands_t, pseudopotential_type
 use m_fstrings,        only : sjoin, ktoa, ftoa, strcat, ltoa, itoa
 use m_time,            only : cwtime_report, cwtime
 use m_io_tools,        only : file_exists, iomode_from_fname, open_file
 use m_pptools,         only : write_xsf
 use m_kpts,            only : kpts_map, kpts_timrev_from_kptopt, bzlint_t, kptrlatt_from_ngkpt
 use m_fft_mesh,        only : supercell_fft, calc_ceikr
 use m_pawtab,          only : pawtab_type
 use m_gstore,          only : gstore_t, gqk_t
 use m_supercell,       only : supercell_type
 use m_paw_sphharm,     only : ylm_angular_mesh
 use m_fftcore,         only : ngfft_seq !, get_kg
 use m_ephtk,           only : ephtk_get_mpw_gmax
 use m_dynmat,          only : phdispl_from_eigvec
 use m_phonons,         only : pheigvec_rotate


 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/polstate_t
!! NAME
!!  polstate_t
!!
!! FUNCTION
!!  Description
!!
!! SOURCE

 type, public :: polstate_t

  integer :: nkbz = -1
  integer :: nqbz = -1

  real(dp) :: gradres
  real(dp) :: eps
  real(dp) :: enel
  real(dp) :: enph
  real(dp) :: enelph

  logical :: has_prev_grad = .false.

  real(dp), allocatable :: eig(:,:)

  real(dp), allocatable :: my_qpts(:,:)
  real(dp), pointer :: my_kpts(:,:) => null()

  complex(dp), allocatable :: my_a(:,:)
  complex(dp), allocatable :: a_glob(:,:)

  complex(dp), allocatable :: my_b(:,:)
  complex(dp), allocatable :: b_glob(:,:)

  complex(dp), allocatable :: my_grad_a(:,:)
  complex(dp), allocatable :: grad_a_glob(:,:)

  complex(dp), allocatable :: my_prevgrad_a(:,:)
  complex(dp), allocatable :: my_pcond(:,:)

  class(gqk_t), pointer :: gqk => null()

  type(krank_t) :: krank_kpts
  type(krank_t) :: krank_qpts

  type(crystal_t) :: cryst

  contains

    procedure :: get_linemin_param => polstate_get_linemin_param

    procedure :: localize => polstate_localize

    procedure :: get_enelph => polstate_get_enelph
    procedure :: get_enph => polstate_get_enph
    procedure :: get_enel => polstate_get_enel

    procedure :: get_grad_a => polstate_get_grad_a
    procedure :: get_conjgrad_a => polstate_get_conjgrad_a
    procedure :: update_pcond => polstate_update_pcond
    procedure :: update_a => polstate_update_a

    procedure :: free => polstate_free

    procedure :: get_b_from_a => polstate_get_b_from_a
    !procedure :: get_b_from_displ =>_get_b_from_displ

    procedure :: seed_a => polstate_seed_a
    procedure :: load_a => polstate_load_a

    procedure :: get_norm => polstate_get_norm
    procedure :: gather => polstate_gather
    procedure :: get_krank_glob => polstate_get_krank_glob
    procedure :: get_mapping => polstate_get_mapping

 end type polstate_t
!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/varpeq_t
!! NAME
!!  varpeq_t
!!
!! FUNCTION
!!
!! SOURCE

 type, public :: varpeq_t

   character(len=fnlen) :: pkind = " "
   character(len=fnlen) :: aseed = " "

   integer :: nstep = -1
   integer :: pc_nupdate = -1
   integer :: ncid = nctk_noid

   integer :: nsppol
   integer :: max_nk
   integer :: max_nq
   integer :: max_nb

   integer :: frohl_ntheta

   real(dp) :: e_frohl

   real(dp) :: pc_factor
   real(dp) :: tolgrs

   real(dp) :: gau_params(2)

   logical :: is_complete = .False.
   logical :: restart = .False.
   logical :: interpolate = .False.

   integer :: ngkpt(3)

   integer, allocatable :: nk_spin(:)
   integer, allocatable :: nq_spin(:)
   integer, allocatable :: nb_spin(:)
   integer, allocatable :: brange_spin(:,:)

   real(dp), allocatable :: kpts_spin(:,:,:)
   real(dp), allocatable :: qpts_spin(:,:,:)
   real(dp), allocatable :: a_spin(:,:,:,:)
   real(dp), allocatable :: b_spin(:,:,:,:)

   integer, allocatable :: nstep2cv_spin(:)
   real(dp), allocatable :: iter_rec_spin(:,:,:)

   class(gstore_t), pointer :: gstore => null()
   type(gaps_t) :: gaps
   type(crystal_t) :: cryst_trinv
   type(polstate_t), allocatable :: polstate(:)

 contains

    procedure :: init => varpeq_init
    procedure :: free => varpeq_free
    procedure :: solve => varpeq_solve
    procedure :: record => varpeq_record
    procedure :: collect => varpeq_collect
    procedure :: print => varpeq_print
    procedure :: ncwrite => varpeq_ncwrite
    procedure :: ncread => varpeq_ncread
    procedure :: compare => varpeq_compare
    procedure :: setup => varpeq_setup

    procedure :: avg_frohlich => varpeq_avg_frohlich

 end type varpeq_t
!!***

 public :: varpeq_run
   ! Main entry point

 public :: varpeq_plot
  ! Compute polaron wavefunctions and atomic displacements in the supercell and write results to XSF files

contains !=====================================================================


!!****f* m_varpeq/varpeq_run
!! NAME
!!  varpeq_run
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_run(gstore, dtset, dtfil)

!Arguments ------------------------------------
 class(gstore_t), intent(in) :: gstore
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type), intent(in) :: dtfil

!Local variables-------------------------------
!scalars
 type(varpeq_t) :: vpq

!----------------------------------------------------------------------

 call vpq%init(gstore, dtset)
 call vpq%setup(dtfil)
 call vpq%solve()
 call vpq%print()
 call vpq%ncwrite(dtset, dtfil)
 call vpq%free()

end subroutine varpeq_run
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_free
!! NAME
!!  varpeq_free
!!
!! FUNCTION
!!  Free allocatable arrays
!!
!! SOURCE

subroutine varpeq_free(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 integer :: my_is

!----------------------------------------------------------------------

 ! Free allocatable arrays
 ABI_SFREE(self%iter_rec_spin)
 ABI_SFREE(self%nstep2cv_spin)
 ABI_SFREE(self%kpts_spin)
 ABI_SFREE(self%qpts_spin)
 ABI_SFREE(self%a_spin)
 ABI_SFREE(self%b_spin)
 ABI_SFREE(self%nk_spin)
 ABI_SFREE(self%nq_spin)
 ABI_SFREE(self%nb_spin)
 ABI_SFREE(self%brange_spin)

 ! Free local datatypes
 call self%cryst_trinv%free()

 ! Close netcdf file
 if (self%ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(self%ncid))
   self%ncid = nctk_noid
 end if

 if (self%is_complete) then
   call self%gaps%free()
   do my_is=1,self%gstore%my_nspins
     call self%polstate(my_is)%free()
   enddo
   self%gstore => null()
 endif

end subroutine varpeq_free
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_compare
!! NAME
!!  varpeq_compare
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_compare(self, other, allow_mesh_mismatch)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(in) :: self, other
 logical, optional, intent(in) :: allow_mesh_mismatch

!Local variables-------------------------------
 character(len=500) :: msg
 integer :: ierr
 !logical :: mismatch
!----------------------------------------------------------------------

 ierr = 0

 ABI_CHECK_NOSTOP(self%pkind == other%pkind, "Difference found in pkind.", ierr)
 ABI_CHECK_NOSTOP(self%nsppol == other%nsppol, "Difference found in nsppol.", ierr)
 msg = " Comparing VARPEQ crystal with GSTORE crystal (time-reversal and inversion symmetries only) "
 ABI_CHECK_NOSTOP(self%cryst_trinv%compare(other%cryst_trinv, header=msg) == 0, "Difference found in cryst.", ierr)
 ABI_CHECK_NOSTOP(self%max_nb == other%max_nb, "Difference found in max_nb.", ierr)
 ABI_CHECK_NOSTOP(all(self%nb_spin == other%nb_spin), "Difference found in nb_spin.", ierr)

 if (present(allow_mesh_mismatch) .and. (.not. allow_mesh_mismatch)) then
   ABI_CHECK_NOSTOP(self%max_nk == other%max_nk, "Difference found in max_nk.", ierr)
   ABI_CHECK_NOSTOP(self%max_nq == other%max_nq, "Difference found in max_nq.", ierr)
   ABI_CHECK_NOSTOP(all(self%nk_spin == other%nk_spin), "Difference found in nk_spin.", ierr)
   ABI_CHECK_NOSTOP(all(self%nq_spin == other%nq_spin), "Difference found in nq_spin.", ierr)
   ABI_CHECK_NOSTOP(all(self%qpts_spin == other%qpts_spin), "Difference found in qpts_spin.", ierr)
   ABI_CHECK_NOSTOP(all(self%kpts_spin == other%kpts_spin), "Difference found in kpts_spin.", ierr)
 endif

 ABI_CHECK(ierr == 0, "Fatal error in varpeq_compare, see previous messages!")

end subroutine varpeq_compare
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_ncread
!! NAME
!!  varpeq_ncread
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_ncread(self, path, comm, keep_open)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 character(len=fnlen), intent(in) :: path
 integer, intent(in) :: comm
 logical, optional, intent(in) :: keep_open

!Local variables-------------------------------
 !character(len=500) :: msg
 integer :: ncid, nsppol, natom3
 real(dp) :: cpu, wall, gflops

!----------------------------------------------------------------------

 ABI_CHECK(file_exists(path), sjoin(" varpeq_ncread: cannot find *VARPEQ.nc file", path))

 call cwtime(cpu, wall, gflops, "start")

 NCF_CHECK(nctk_open_read(ncid, path, comm))

 ! Read crystal structure
 call self%cryst_trinv%ncread(ncid)

 ! Read varpeq dimensions
 NCF_CHECK(nctk_get_dim(ncid, "max_nk", self%max_nk))
 NCF_CHECK(nctk_get_dim(ncid, "max_nq", self%max_nq))
 NCF_CHECK(nctk_get_dim(ncid, "max_nb", self%max_nb))
 NCF_CHECK(nctk_get_dim(ncid, "nsppol", self%nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "natom3", natom3))

 ! Read data
 ! arrays
 nsppol = self%nsppol
 !print *, "self%nsppol, self%max_nk, self%max_nq, self%max_nb", self%nsppol, self%max_nk, self%max_nq, self%max_nb
 ABI_MALLOC(self%kpts_spin, (3, self%max_nk, nsppol))
 ABI_MALLOC(self%qpts_spin, (3, self%max_nq, nsppol))
 ABI_MALLOC(self%a_spin, (2, self%max_nb, self%max_nk, nsppol))
 ABI_MALLOC(self%b_spin, (2, natom3, self%max_nq, nsppol))
 ABI_MALLOC(self%nk_spin, (nsppol))
 ABI_MALLOC(self%nq_spin, (nsppol))
 ABI_MALLOC(self%nb_spin, (nsppol))
 ABI_MALLOC(self%brange_spin, (2, self%nsppol))

 NCF_CHECK(nf90_get_var(ncid, vid("kpts_spin"), self%kpts_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("qpts_spin"), self%qpts_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("a_spin"), self%a_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("b_spin"), self%b_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("varpeq_pkind"), self%pkind))
 NCF_CHECK(nf90_get_var(ncid, vid("nk_spin"), self%nk_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("nq_spin"), self%nq_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("nb_spin"), self%nb_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("brange_spin"), self%brange_spin))
 NCF_CHECK(nf90_get_var(ncid, vid("ngkpt"), self%ngkpt))

 if (present(keep_open) .and. keep_open) then
   self%ncid = ncid
 else
   NCF_CHECK(nf90_close(ncid))
   self%ncid = nctk_noid
 end if

 call cwtime_report(" varpeq_ncread", cpu, wall, gflops)

 self%is_complete = .False.

!----------------------------------------------------------------------

 contains
  integer function vid(var_name)
    character(len=*),intent(in) :: var_name
    vid = nctk_idname(ncid, var_name)
 end function vid

end subroutine varpeq_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_ncwrite
!! NAME
!!  varpeq_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_ncwrite(self, dtset, dtfil)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
 character(len=fnlen) :: path
 integer, parameter :: master = 0
 integer :: my_rank, ncid, ncerr
 real(dp) :: cpu, wall, gflops

!----------------------------------------------------------------------

 ! Shamelessly copied from (inspired by) the sigmaph%write routine

 my_rank = xmpi_comm_rank(self%gstore%comm)

 call cwtime(cpu, wall, gflops, "start")

 ! Create netcdf file (only master works, HDF5 + MPI-IO is handled afterwards by reopening the file inside ncwrite_comm)
 path = strcat(dtfil%filnam_ds(4), "_VARPEQ.nc")
 if (my_rank == master) then
   ! Master creates the netcdf file used to store the results of the calculation.
   NCF_CHECK(nctk_open_create(self%ncid, path, xmpi_comm_self))
   ncid = self%ncid

   ! Write the crystal (TR & invsersion symmetry only) & ebands dataset_type
   NCF_CHECK(self%cryst_trinv%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(self%gstore%ebands, ncid))

   ! Add varpeq dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("max_nk", self%max_nk), nctkdim_t("max_nq", self%max_nq), &
     nctkdim_t("max_nb", self%max_nb), nctkdim_t("nsppol", self%gstore%nsppol), &
     nctkdim_t("natom3", 3*self%gstore%cryst%natom), &
     nctkdim_t("nstep", self%nstep)], &
     defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define scalars
   ! integers
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "varpeq_nstep", "nkbz", "nqbz", "frohl_ntheta"])
   NCF_CHECK(ncerr)
   ! real(dp)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "varpeq_tolgrs", "e_frohl"])
   NCF_CHECK(ncerr)

   ! Define arrays with results
   ! FIXME: correct a_spin/b_spin representation based on Matteo's advice
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("varpeq_pkind", "c", "fnlen"), &
     nctkarr_t("varpeq_aseed", "c", "fnlen"), &
     nctkarr_t("varpeq_gau_params", "dp", "two"), &
     nctkarr_t("nstep2cv", "int", "nsppol"), &
     nctkarr_t("iter_rec", "dp", "six, nstep, nsppol"), &
     nctkarr_t("nk_spin", "int", "nsppol"), &
     nctkarr_t("nq_spin", "int", "nsppol"), &
     nctkarr_t("nb_spin", "int", "nsppol"), &
     nctkarr_t("brange_spin", "int", "two, nsppol"), &
     nctkarr_t("kpts_spin", "dp", "three, max_nk, nsppol"), &
     nctkarr_t("qpts_spin", "dp", "three, max_nq, nsppol"), &
     nctkarr_t("a_spin", "dp", "two, max_nb, max_nk, nsppol"), &
     nctkarr_t("b_spin", "dp", "two, natom3, max_nq, nsppol"), &
     nctkarr_t("cb_min_spin", "dp", "nsppol"), &
     nctkarr_t("vb_max_spin", "dp", "nsppol"), &
     nctkarr_t("gstore_ngqpt", "i", "three"), &
     nctkarr_t("ngkpt", "i", "three") &
   ])
   NCF_CHECK(ncerr)

   ! Write data
   NCF_CHECK(nctk_set_datamode(ncid))
   ! scalars
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task", "varpeq_nstep", "nkbz", "nqbz", "frohl_ntheta"], &
     [dtset%eph_task, self%nstep, self%gstore%nkbz, self%gstore%nqbz, self%frohl_ntheta])
   NCF_CHECK(ncerr)
   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "varpeq_tolgrs", "e_frohl"], &
     [self%tolgrs, self%e_frohl])
   NCF_CHECK(ncerr)
   ! arrays
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "varpeq_pkind"), self%pkind))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "varpeq_aseed"), self%aseed))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "varpeq_gau_params"), self%gau_params))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nstep2cv"), self%nstep2cv_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "iter_rec"), self%iter_rec_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nk_spin"), self%nk_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nq_spin"), self%nq_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nb_spin"), self%nb_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "brange_spin"), self%brange_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kpts_spin"), self%kpts_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpts_spin"), self%qpts_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "a_spin"), self%a_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "b_spin"), self%b_spin))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cb_min_spin"), self%gaps%cb_min))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vb_max_spin"), self%gaps%vb_max))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gstore_ngqpt"), self%gstore%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngkpt"), self%ngkpt))
 end if ! master

 call xmpi_barrier(self%gstore%comm)
 call cwtime_report(" varpeq: netcdf", cpu, wall, gflops)

end subroutine varpeq_ncwrite
!!***

!!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_print
!! NAME
!!  varpeq_print
!!
!! FUNCTION
!! Compute polaron wavefunctions and atomic displacements in the supercell and write results to files
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_print(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
!scalars
 character(len=5000) :: msg
 integer, parameter :: master = 0
 integer :: comm, nproc, my_rank, spin, ii
 real(dp) :: enpol, enel, enph, enelph, eps, grs
!arrays
 integer :: units(2)

!----------------------------------------------------------------------

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! FIXME: the current state of this output procedure is ugly
 units = [std_out, ab_out]
 if (my_rank == master) then

   write(msg, '(a)') " "
   call wrtout(units, msg)
   write(msg, '(a)') "  === Variational Polaron Equations ==="
   call wrtout(units, msg)
   write(msg, '(a)') repeat('-', 80)
   call wrtout(units, msg)

   ! For each spin, print the SCF cycle convergence
   do spin=1,self%gstore%nsppol

     write(msg, '(a,i1,a,i1)') "  * spin: ", spin, "/", self%gstore%nsppol
     call wrtout(units, msg)

     write(msg,'(a4,a13,2a12,a13,a13,a13)') 'Step', 'E_pol', 'E_el', 'E_ph', &
       'E_elph', 'epsilon', '||gradient||'
     call wrtout(units, msg)

     do ii=1,self%nstep2cv_spin(spin)
       enpol = self%iter_rec_spin(1, ii, spin); enel = self%iter_rec_spin(2, ii, spin)
       enph = self%iter_rec_spin(3, ii, spin); enelph = self%iter_rec_spin(4, ii, spin)
       eps = self%iter_rec_spin(5, ii, spin); grs = self%iter_rec_spin(6, ii, spin)

       write(msg,'(i4,es13.4,2es12.4,es13.4,es13.4,es13.4)') ii, enpol, enel, &
         enph, enelph, eps, grs
       call wrtout(units, msg)

     enddo
     write(msg, '(a)') repeat('-', 80)
     call wrtout(units, msg)

     ! Tell the user whether convergence is achieved
     if (grs < self%tolgrs) then
       write(msg, '(a,es8.2,a,es8.2)') "gradient norm ", grs, " < varpeq_tolgrs=", self%tolgrs
       call wrtout(units, msg)
       call wrtout(units, "Calculation is completed and convergence reached.")
     else
       call wrtout(units, "Calculation is completed but not converged.")
       write(msg, '(a,es8.2,a,es8.2,a,a)') "gradient norm ", grs, &
         " >= varpeq_tolgrs=", self%tolgrs, ch10, &
         "You may want to restart the SCF cycle with varpeq_restart=1 option"
       ABI_WARNING(msg)
     endif

     ! Final results
     write(ab_out, '(a,es13.6,a)') "polaron binding energy = ", enpol, 'a.u.'
     write(ab_out, '(a,es13.6,a)') "polaron localized state energy = ", eps, 'a.u.'
     write(ab_out, '(a,es13.6,a)') "long-range divergence correction = ", self%e_frohl, 'a.u.'

   enddo

 endif

end subroutine varpeq_print
!!***

!!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_collect
!! NAME
!!  varpeq_collect
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_collect(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 class(polstate_t), pointer :: polstate
 integer :: ierr, my_is, spin, my_ik, ik_glob, my_iq, iq_glob
 integer :: my_pert, pert_glob
 integer :: oc_iter, oc_a, oc_b, oc_k, oc_q

!----------------------------------------------------------------------

 ! FIXME
 ! Hack to mimic the summation over a non-existing spin commnicator

 ! Gather the SCF process evolution data
 call xmpi_sum(self%iter_rec_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%nstep2cv_spin, self%gstore%comm, ierr)

 ! Gather electron/phonon vectors and k/q points
 self%a_spin(:,:,:,:) = zero
 self%b_spin(:,:,:,:) = zero
 self%qpts_spin(:,:,:) = zero
 self%kpts_spin(:,:,:) = zero
 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   gqk => self%gstore%gqk(my_is)
   polstate => self%polstate(spin)

   ! electronic vector & k-points
   do my_ik=1,gqk%my_nk
     ik_glob = gqk%my_kstart + my_ik - 1
     self%a_spin(1, :, ik_glob, spin) = real(polstate%my_a(:, my_ik))  ! real part
     self%a_spin(2, :, ik_glob, spin) = aimag(polstate%my_a(:, my_ik)) ! imaginary part
     self%kpts_spin(:, ik_glob, spin) = polstate%my_kpts(:, my_ik)
   enddo

   ! phonon vector & q-points
   do my_iq=1,gqk%my_nq
     iq_glob = gqk%my_qstart + my_iq - 1
     do my_pert=1,gqk%my_npert
       pert_glob = gqk%my_pert_start + my_pert - 1
       self%b_spin(1, pert_glob, iq_glob, spin) = real(polstate%my_b(my_pert, my_iq))  ! real part
       self%b_spin(2, pert_glob, iq_glob, spin) = aimag(polstate%my_b(my_pert, my_iq)) ! imaginary part
     enddo
     self%qpts_spin(:, iq_glob, spin) = polstate%my_qpts(:, my_iq)
   enddo

 enddo
 call xmpi_sum(self%a_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%b_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%kpts_spin, self%gstore%comm, ierr)
 call xmpi_sum(self%qpts_spin, self%gstore%comm, ierr)

 ! Don't forget to divide all by the #OverCount, since we use global communicator
 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   gqk => self%gstore%gqk(my_is)
   polstate => self%polstate(spin)

   oc_iter = gqk%comm%nproc
   oc_a = gqk%qpt_pert_comm%nproc
   oc_b = gqk%kpt_comm%nproc
   oc_k = oc_a
   oc_q = oc_b * gqk%pert_comm%nproc

   self%iter_rec_spin(:,:,spin) = self%iter_rec_spin(:,:,spin) / oc_iter
   self%nstep2cv_spin(spin) = self%nstep2cv_spin(spin) / oc_iter
   self%a_spin(:,:,:,spin) = self%a_spin(:,:,:,spin) / oc_a
   self%b_spin(:,:,:,spin) = self%b_spin(:,:,:,spin) / oc_b
   self%kpts_spin(:,:,spin) = self%kpts_spin(:,:,spin) / oc_k
   self%qpts_spin(:,:,spin) = self%qpts_spin(:,:,spin) / oc_q
 enddo

end subroutine varpeq_collect
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_setup
!! NAME
!!  varpeq_setup
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_setup(self, dtfil)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(datafiles_type), intent(in) :: dtfil

!Local variables-------------------------------
!scalars
 character(len=5000) :: msg
 class(polstate_t), pointer :: polstate
 type(varpeq_t) :: vpq_loaded
 type(bzlint_t) :: bzlint
 integer, parameter :: master = 0
 integer :: my_rank, comm, ierr, my_is, spin, nk, nb, ik, ib
!arrays
 integer :: units(2), ngkpt(3)
 real(dp) :: kpt(3)
 real(dp), allocatable :: a_loaded(:,:), kpts_loaded(:,:), ank(:)

!----------------------------------------------------------------------

 ! If a compatible VARPEQ.nc file is avaliable along with special flags, use it to either
 ! restart the calculation or interpolate initial charge localization

 units = [std_out, ab_out]
 comm = self%gstore%comm; my_rank = xmpi_comm_rank(comm)

 if (my_rank == master .and. (self%interpolate .or. self%restart)) then
   call wrtout(units, " - getting charge localiztion from a previous VARPEQ.nc file")
   call vpq_loaded%ncread(dtfil%filvarpeqin, xmpi_comm_self, keep_open=.False.)

   if (self%interpolate .and. self%restart) then
     msg = " Both varpeq_interpolate and eph_restart flags are provided, &
       priority is given to the interpolation"
     ABI_WARNING(msg)
   endif

   if (self%interpolate) then
     call wrtout(units, " - interpolating previous A_nk")
     call self%compare(vpq_loaded, allow_mesh_mismatch=.True.)

     do spin=1,self%nsppol
       ! prepare for interpolation
       nk = vpq_loaded%nk_spin(spin)
       nb = vpq_loaded%nb_spin(spin)
       ngkpt(:) = vpq_loaded%ngkpt(:)
       ABI_MALLOC(a_loaded, (2*nb, nk))
       ABI_MALLOC(kpts_loaded, (3, nk))
       ABI_MALLOC(ank, (2*nb))

       do ik=1,nk
         kpts_loaded(:, ik) = vpq_loaded%kpts_spin(:, ik, spin)
         do ib=1,nb
           a_loaded(2*ib-1, ik) = vpq_loaded%a_spin(1, ib, ik, spin) ! real
           a_loaded(2*ib, ik) = vpq_loaded%a_spin(2, ib, ik, spin) ! imaginary
         enddo
       enddo

       ! Now, the actual interpolation is performed
       call bzlint%init(ngkpt, 2*nb, nk, kpts_loaded, a_loaded)
       self%a_spin(:,:,:,spin) = zero

       do ik=1,self%nk_spin(spin)
         kpt = self%kpts_spin(:, ik, spin)
         call bzlint%interp(kpt, ank)

         do ib=1,self%nb_spin(spin)
           self%a_spin(1, ib, ik, spin) = ank(2*ib-1)
           self%a_spin(2, ib, ik, spin) = ank(2*ib)

         enddo
       enddo

       call bzlint%free()
       ABI_FREE(ank)
       ABI_FREE(kpts_loaded)
       ABI_FREE(a_loaded)
     enddo

   else
     call wrtout(units, " - restarting from previous A_nk")
     call self%compare(vpq_loaded, allow_mesh_mismatch=.False.)
     self%a_spin(:,:,:,:) = vpq_loaded%a_spin(:,:,:,:)

   endif

   call vpq_loaded%free()
 endif
 call xmpi_bcast(self%a_spin, master, comm, ierr)

 ! Initialize charge localization in each polstate

 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   polstate => self%polstate(my_is)

   if (self%interpolate .or. self%restart) then
     call polstate%load_a(self%a_spin(:,:,:,spin))
   else
     call polstate%seed_a(self%aseed, gau_params=self%gau_params)
   endif

 enddo
end subroutine varpeq_setup
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_solve
!! NAME
!!  varpeq_solve
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_solve(self)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self

!Local variables-------------------------------
 !character(len=5000) :: msg
 class(polstate_t), pointer :: polstate
 integer :: my_is, spin, ii

!----------------------------------------------------------------------

 ! Variational Polaron Equations
 self%iter_rec_spin(:,:,:) = zero

 do my_is=1,self%gstore%my_nspins
   spin = self%gstore%my_spins(my_is)
   polstate => self%polstate(my_is)

   ! call polstate%seed_a(self%aseed, gau_params=self%gau_params)

   do ii=1,self%nstep
     call polstate%localize()

     ! Save the necessary data at each iterations
     call self%record(ii, my_is)

     if (polstate%gradres < self%tolgrs) exit

     ! TODO: add more control for this feature
     ! After some iterations update the preconditioner depending on the value of eps
     if (self%pc_nupdate == 1) call polstate%update_pcond(self%pc_factor)
     if (mod(ii, self%pc_nupdate) == 1) call polstate%update_pcond(self%pc_factor)

     call polstate%get_conjgrad_a()
     call polstate%update_a()
   enddo

 enddo

 ! Correction for the infrared Fr\"ohlich divergence
 call self%avg_frohlich()

 ! Collect results from each polstate
 call self%collect()

end subroutine varpeq_solve
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_record
!! NAME
!!  varpeq_record
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_record(self, iter, my_is)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 integer, intent(in) :: iter, my_is

!Local variables-------------------------------
 class(polstate_t), pointer :: polstate
 integer :: spin, psign
 real(dp) :: shift

!----------------------------------------------------------------------

 spin = self%gstore%my_spins(my_is)
 polstate => self%polstate(my_is)

 select case(self%pkind)
 case ("electron")
   psign = 1; shift = self%gaps%cb_min(spin)
 case ("hole")
   psign = -1; shift = -self%gaps%vb_max(spin)
 end select

 self%iter_rec_spin(1, iter, spin) = &
   psign*(polstate%enel + polstate%enph + polstate%enelph - shift)
 self%iter_rec_spin(2, iter, spin) = psign*(polstate%enel - shift)
 self%iter_rec_spin(3, iter, spin) = psign*polstate%enph
 self%iter_rec_spin(4, iter, spin) = psign*polstate%enelph
 self%iter_rec_spin(5, iter, spin) = psign*(polstate%eps - shift)
 self%iter_rec_spin(6, iter, spin) = polstate%gradres
 self%nstep2cv_spin(spin) = iter

end subroutine varpeq_record
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_init
!! NAME
!!  varpeq_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_init(self, gstore, dtset)

!Arguments ------------------------------------
 class(varpeq_t), target, intent(inout) :: self
 type(gstore_t), target, intent(in) :: gstore
 type(dataset_type), intent(in) :: dtset

!Local variables-------------------------------
!scalars
 character(len=5000) :: msg
 class(gqk_t), pointer :: gqk
 class(polstate_t), pointer :: polstate
 integer :: ierr, my_is, spin, bstart, my_iq
 real(dp) :: wtq

!----------------------------------------------------------------------

 ! Consistency check
 !if (gstore%check_cplex_qkzone_gmode(2, "bz", "bz", "phonon") /= 0) then
 !  ABI_ERROR("The gstore object is inconsistent with varpeq. See messages above.")
 !end if

 self%gstore => gstore

 self%nstep = dtset%varpeq_nstep
 self%tolgrs = dtset%varpeq_tolgrs
 self%pc_nupdate = dtset%varpeq_pc_nupdate
 self%pc_factor = dtset%varpeq_pc_factor
 self%aseed = dtset%varpeq_aseed
 self%pkind = dtset%varpeq_pkind
 self%gau_params = dtset%varpeq_gau_params
 self%frohl_ntheta = dtset%eph_frohl_ntheta
 self%gaps = ebands_get_gaps(gstore%ebands, ierr)
 self%cryst_trinv = gstore%cryst%new_trinv_only()

 self%nsppol = gstore%nsppol
 self%max_nk = maxval(gstore%glob_nk_spin)
 self%max_nq = maxval(gstore%glob_nq_spin)
 self%max_nb = maxval(gstore%brange_spin(2,:) - gstore%brange_spin(1,:)) + 1

 self%restart = (dtset%eph_restart == 1)
 self%interpolate = (dtset%varpeq_interpolate == 1)

 ABI_MALLOC(self%polstate, (gstore%my_nspins))
 ABI_MALLOC(self%iter_rec_spin, (6, self%nstep, gstore%nsppol))
 ABI_MALLOC(self%nstep2cv_spin, (gstore%nsppol))
 ABI_MALLOC(self%kpts_spin, (3, self%max_nk, gstore%nsppol))
 ABI_MALLOC(self%qpts_spin, (3, self%max_nq, gstore%nsppol))
 ABI_MALLOC(self%a_spin, (2, self%max_nb, self%max_nk, gstore%nsppol))
 ABI_MALLOC(self%b_spin, (2, 3*gstore%cryst%natom, self%max_nq, gstore%nsppol))
 ABI_MALLOC(self%nk_spin, (gstore%nsppol))
 ABI_MALLOC(self%nq_spin, (gstore%nsppol))
 ABI_MALLOC(self%nb_spin, (gstore%nsppol))
 ABI_MALLOC(self%brange_spin, (2, gstore%nsppol))

 self%nk_spin(:) = gstore%glob_nk_spin(:)
 self%nq_spin(:) = gstore%glob_nq_spin(:)
 self%nb_spin(:) = gstore%brange_spin(2,:) - gstore%brange_spin(1,:) + 1
 self%brange_spin(:,:) = gstore%brange_spin(:,:)
 self%ngkpt(:) = dtset%ngkpt(:)

 ! Loop over my spins and initialize polaronic states
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is)
   gqk => gstore%gqk(my_is)
   polstate => self%polstate(my_is)

   polstate%gqk => gqk
   polstate%nkbz = gstore%nkbz
   polstate%nqbz = gstore%nqbz
   polstate%cryst = gstore%cryst%new_without_symmetries()

   ! Basic dimensions
   ABI_MALLOC(polstate%my_a, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_grad_a, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_prevgrad_a, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_pcond, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(polstate%my_b, (gqk%my_npert, gqk%my_nq))

   ABI_MALLOC(polstate%a_glob, (gqk%nb, gqk%glob_nk))
   ABI_MALLOC(polstate%grad_a_glob, (gqk%nb, gqk%glob_nk))
   ABI_MALLOC(polstate%b_glob, (gqk%my_npert, gqk%glob_nq))
   ABI_MALLOC(polstate%eig, (gqk%nb, gstore%ebands%nkpt))

   ! Bands taking part in the polaron formation process
   ! TODO: shift the bands wrt vbm/cbm?
   msg = sjoin(self%gaps%errmsg_spin(spin), ". VarPEq is incompatible with metals &
     and needs CBM/VBM for electron/hole polaron calculations.")
   ABI_CHECK(self%gaps%ierr(spin) == 0, msg)

   bstart = self%gstore%brange_spin(1, spin)
   select case(self%pkind)
   case ("electron")
     polstate%eig = gstore%ebands%eig(bstart:bstart+gqk%nb-1, :, spin)
   case ("hole")
     ! here we flip the valence bands to deal with the minimization process later on
     polstate%eig = -gstore%ebands%eig(bstart:bstart+gqk%nb-1, :, spin)
   end select

   ! k-space and q-space treated by this proc
   polstate%my_kpts => gqk%my_kpts(:,:)

   ABI_MALLOC(polstate%my_qpts, (3, gqk%my_nq))
   do my_iq=1,gqk%my_nq
     call gqk%myqpt(my_iq, gstore, wtq, polstate%my_qpts(:, my_iq))
   enddo

   ! kranks (require initializaton of polstata%my_kpts/my_qpts)
   polstate%krank_kpts = polstate%get_krank_glob("k", gstore%ebands%kptrlatt)
   polstate%krank_qpts = polstate%get_krank_glob("q", gstore%ebands%kptrlatt)

 enddo

 ! FIXME: fix this hack; it's needed only to setup kpts_spn & qpts_spin prior to all calculations
 call self%collect()

 self%is_complete = .True.

end subroutine varpeq_init
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/varpeq_avg_frohlich
!! NAME
!!  varpeq_avg_frohlich
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_avg_frohlich(self)

!Arguments ------------------------------------
 class(varpeq_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 type(ifc_type), pointer :: ifc
 type(crystal_t), pointer :: cryst
 integer :: comm, my_rank, nproc, ierr
 integer :: ntheta, nphi, angl_size
 integer :: iang, iatom
 integer :: nu, natom3
 real(dp) :: inv_qepsq, wqnu
 complex(dpc) :: cnum
!arrays
 real(dp) :: qpt_cart(3)
 real(dp), allocatable :: phfreq(:), displ_cart(:, :, :, :)
 real(dp), allocatable :: qvers_cart(:, :)
 real(dp), allocatable :: angweight(:)
 complex(dpc) :: cp3(3)

!----------------------------------------------------------------------

 comm = self%gstore%comm
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ifc => self%gstore%ifc
 cryst => self%gstore%cryst

 natom3 = 3*cryst%natom
 ABI_MALLOC(phfreq, (natom3))
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))

 ntheta = self%frohl_ntheta
 nphi = 2*ntheta
 call ylm_angular_mesh(ntheta, nphi, angl_size, qvers_cart, angweight)

 self%e_frohl = zero
 do iang=1,angl_size
   if (mod(iang, nproc) /= my_rank) cycle ! MPI parallelism inside comm.

   qpt_cart = qvers_cart(:, iang)
   inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
   call ifc%fourq(cryst, qpt_cart, phfreq, displ_cart, nanaqdir="cart")

   do nu=4,natom3
     wqnu = phfreq(nu)

     cp3 = zero
     do iatom=1,cryst%natom
       cp3 = cp3 + matmul(ifc%zeff(:, :, iatom), &
         cmplx(displ_cart(1, :, iatom, nu), displ_cart(2, :, iatom, nu), kind=dpc))
     enddo
     cnum = dot_product(qpt_cart, cp3)
     if (abs(cnum) < tol12) cycle

     self%e_frohl = self%e_frohl + angweight(iang)*(abs(cnum)*inv_qepsq/wqnu)**2
   enddo
 enddo

 call xmpi_sum(self%e_frohl, comm, ierr)

 self%e_frohl = self%e_frohl * &
   eight*pi/cryst%ucvol * (three / (four_pi * cryst%ucvol * self%gstore%nqbz))**third
 ! for an electron polaron, the correction has to be negative
 if (self%pkind == "electron") self%e_frohl = -self%e_frohl

 ABI_FREE(phfreq)
 ABI_FREE(displ_cart)
 ABI_FREE(qvers_cart)
 ABI_FREE(angweight)

end subroutine varpeq_avg_frohlich
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_free
!! NAME
!!  polstate_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_free(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!----------------------------------------------------------------------

 self%has_prev_grad = .false.

 self%my_kpts => null()
 self%gqk => null()

 ABI_SFREE(self%my_qpts)

 ABI_SFREE(self%my_a)
 ABI_SFREE(self%my_b)
 ABI_SFREE(self%my_grad_a)
 ABI_SFREE(self%my_prevgrad_a)
 ABI_SFREE(self%my_pcond)

 ABI_SFREE(self%a_glob)
 ABI_SFREE(self%b_glob)
 ABI_SFREE(self%grad_a_glob)

 call self%krank_kpts%free()
 call self%krank_qpts%free()
 call self%cryst%free()

end subroutine polstate_free
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_update_a
!! NAME
!!  polstate_update_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_update_a(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
 real(dp) :: theta

!----------------------------------------------------------------------

 theta = self%get_linemin_param()
 self%my_a(:,:) = cos(theta)*self%my_a(:,:) + sin(theta)*self%my_grad_a(:,:)
 call self%gather('a')

end subroutine polstate_update_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_update_pcond
!! NAME
!!  polstate_update_pcond
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_update_pcond(self, factor)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self
 real(dp),intent(in) :: factor

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: my_ik, ik_ibz, ib

!----------------------------------------------------------------------

 gqk => self%gqk

 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     self%my_pcond(ib, my_ik) = one/(self%eig(ib, ik_ibz) - factor*self%eps)
   enddo
 enddo

end subroutine polstate_update_pcond
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_conjgrad_a
!! NAME
!!  polstate_get_conjgrad_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_get_conjgrad_a(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr
 real(dp) :: beta_den
 complex(dp) :: beta, beta_num

!----------------------------------------------------------------------

 gqk => self%gqk

 ! Apply the preconditioner and orthonormalize wrt A
 self%my_grad_a(:,:) = self%my_pcond(:,:)*self%my_grad_a(:,:)
 call orthonorm_to_a_()

 ! If previous gradient is known, get conjugate gradient using Polak-Ribi'ere formula
 if (self%has_prev_grad) then
   beta_num = &
     sum(self%my_grad_a(:,:)*conjg(self%my_grad_a(:,:) - self%my_prevgrad_a(:,:)))
   beta_den = sum(abs(self%my_prevgrad_a(:,:))**2)
   call xmpi_sum(beta_num, gqk%kpt_comm%value, ierr)
   call xmpi_sum(beta_den, gqk%kpt_comm%value, ierr)
   beta = beta_num/beta_den

   self%my_grad_a(:,:) = self%my_grad_a(:,:) + beta*self%my_prevgrad_a(:,:)
   call orthonorm_to_a_()
 endif

 ! Save the current gradient for each processor and globally
 call self%gather('grad_a')
 self%my_prevgrad_a(:,:) = self%my_grad_a(:,:)
 self%has_prev_grad = .true.

!----------------------------------------------------------------------

 contains
 subroutine orthonorm_to_a_()

  integer :: ierr
  real(dp) :: orth_factor, norm

 !----------------------------------------------------------------------

  orth_factor = real(sum(self%my_a(:,:)*conjg(self%my_grad_a(:,:))))
  call xmpi_sum(orth_factor, gqk%kpt_comm%value, ierr)
  self%my_grad_a(:,:) = &
    self%my_grad_a(:,:) - (orth_factor/self%nkbz)*self%my_a(:,:)

  norm = self%get_norm('grad_a')
  self%my_grad_a(:,:) = self%my_grad_a(:,:)*sqrt(one*self%nkbz)/norm

 end subroutine orthonorm_to_a_

end subroutine polstate_get_conjgrad_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_linemin_param
!! NAME
!!  polstate_get_linemin_param
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

real(dp) function polstate_get_linemin_param(self) result(theta)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: ierr
 integer :: my_iq, my_pert
 integer :: my_ik, ik_ibz, ik_forw, ib, jb
 real(dp) :: term_sin2, term_sincos
 complex(dp) :: a_from, a_forw, d_from, d_forw
 complex(dp) :: g_forw
 complex(dp) :: b
!arrays
 integer, allocatable :: kpq_map(:, :)
 real(dp) :: kpt(3), kpq(3)

!----------------------------------------------------------------------

 gqk => self%gqk
 ABI_MALLOC(kpq_map, (6, gqk%my_nq))

 ! Scattering-dependent part
 term_sin2 = zero
 term_sincos = zero
 do my_ik=1,gqk%my_nk
   kpt(:) = self%my_kpts(:, my_ik)

   !call self%get_mapping('kq', kpt, kpq_map)

   do my_iq=1,gqk%my_nq

     !ik_forw = kpq_map(1, my_iq)

     kpq(:) = kpt(:) + self%my_qpts(:, my_iq)
     ik_forw = self%krank_kpts%get_index(kpq)
     if (ik_forw == -1) cycle

     do ib=1,gqk%nb
       a_from = self%my_a(ib, my_ik)
       d_from = self%my_grad_a(ib, my_ik)
       do jb=1,gqk%nb
         a_forw = self%a_glob(jb, ik_forw)
         d_forw = self%grad_a_glob(jb, ik_forw)
         do my_pert=1,gqk%my_npert
           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
           b = self%my_b(my_pert, my_iq)

           term_sin2 = term_sin2 + real(d_from*conjg(b)*g_forw*conjg(d_forw))
           term_sincos = term_sincos + &
             real((a_from*conjg(d_forw) + d_from*conjg(a_forw))*conjg(b)*g_forw)
         enddo
       enddo
     enddo

   enddo
 enddo
 call xmpi_sum(term_sin2, gqk%qpt_pert_comm%value, ierr)
 call xmpi_sum(term_sincos, gqk%qpt_pert_comm%value, ierr)
 term_sin2 = -two*term_sin2/self%nqbz
 term_sincos = -two*term_sincos/self%nqbz

 ! Scattering-independent part
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     a_from = self%my_a(ib, my_ik)
     d_from = self%my_grad_a(ib, my_ik)

     term_sin2 = term_sin2 + self%eig(ib, ik_ibz)*abs(d_from)**2
     term_sincos = term_sincos + &
       self%eig(ib, ik_ibz)*real(d_from*conjg(a_from) + conjg(d_from)*a_from)
   enddo
 enddo
 call xmpi_sum(term_sin2, gqk%kpt_comm%value, ierr)
 call xmpi_sum(term_sincos, gqk%kpt_comm%value, ierr)
 term_sin2 = term_sin2/self%nkbz
 term_sincos = term_sincos/self%nkbz

 theta = half*atan2(-term_sincos, term_sin2 - self%eps)

 ABI_FREE(kpq_map)

end function polstate_get_linemin_param
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_grad_a
!! NAME
!!  polstate_get_grad_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_get_grad_a(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: ierr
 integer :: my_iq, my_pert
 integer :: my_ik, ik_ibz, ik_forw, ik_back, ib, jb
 complex(dp) :: a_forw, a_back
 complex(dp) :: g_forw, g_back
 complex(dp) :: b
!arrays
 integer, allocatable :: qpk_map(:, :), qmk_map(:, :)
 real(dp) :: qpt(3), kpq(3), kmq(3)
 complex(dp), allocatable :: gq_gathered(:, :, :, :)

!----------------------------------------------------------------------

 gqk => self%gqk

 ABI_MALLOC(qpk_map, (6, gqk%my_nk))
 ABI_MALLOC(qmk_map, (6, gqk%my_nk))
 ABI_MALLOC(gq_gathered, (gqk%my_npert, gqk%nb, gqk%nb, gqk%glob_nk))

 ! Scattering-dependent part
 self%my_grad_a(:, :) = zero
 do my_iq=1,gqk%my_nq
   qpt(:) = self%my_qpts(:, my_iq)

   !call self%get_mapping('qk', qpt, qpk_map)
   !call self%get_mapping('qk', -qpt, qmk_map)

   call gqk%gather("q", my_iq, gq_gathered)

   do my_ik=1,gqk%my_nk

     !ik_forw = qpk_map(1, my_ik)
     !ik_back = qmk_map(1, my_ik)

     kpq(:) = self%my_kpts(:, my_ik) + qpt(:)
     kmq(:) = self%my_kpts(:, my_ik) - qpt(:)
     ik_forw = self%krank_kpts%get_index(kpq)
     ik_back = self%krank_kpts%get_index(kmq)
     if ((ik_forw == -1) .and. (ik_back == -1)) cycle

     ! forward scattering
     if (ik_forw /= -1) then
       do ib=1,gqk%nb
         do jb=1,gqk%nb
           a_forw = self%a_glob(jb, ik_forw)

           do my_pert=1,gqk%my_npert
             b = self%my_b(my_pert, my_iq)
             g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)

             self%my_grad_a(ib, my_ik) = self%my_grad_a(ib, my_ik) + a_forw*b*conjg(g_forw)
           enddo
         enddo
       enddo
     endif

     ! backward scattering
     if (ik_back /= -1) then
       do ib=1,gqk%nb
         do jb=1,gqk%nb
           a_back = self%a_glob(jb, ik_back)

           do my_pert=1,gqk%my_npert
             b = self%my_b(my_pert, my_iq)
             g_back= gq_gathered(my_pert, ib, jb, ik_back)

             self%my_grad_a(ib, my_ik) = self%my_grad_a(ib, my_ik) + a_back*conjg(b)*g_back
           enddo
         enddo
       enddo
     endif

     !do ib=1,gqk%nb
     !  do jb=1,gqk%nb
     !    a_forw = self%a_glob(jb, ik_forw)
     !    a_back = self%a_glob(jb, ik_back)
     !    do my_pert=1,gqk%my_npert
     !      b = self%my_b(my_pert, my_iq)
     !      g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
     !      g_back = gq_gathered(my_pert, ib, jb, ik_back)
     !      self%my_grad_a(ib, my_ik) = self%my_grad_a(ib, my_ik) + &
     !        (a_forw*b*conjg(g_forw) + a_back*conjg(b)*g_back)
     !    enddo
     !  enddo
     !enddo

   enddo
 enddo
 call xmpi_sum(self%my_grad_a, gqk%qpt_pert_comm%value, ierr)
 self%my_grad_a(:, :) = -two*self%my_grad_a(:, :)/(self%nkbz*self%nqbz)

 ! Scattering-independent part
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     self%my_grad_a(ib, my_ik) = self%my_grad_a(ib, my_ik) + &
       (self%eig(ib, ik_ibz) - self%eps)*self%my_a(ib, my_ik)*two/self%nkbz
   enddo
 enddo

 ABI_FREE(qpk_map)
 ABI_FREE(qmk_map)
 ABI_FREE(gq_gathered)

end subroutine polstate_get_grad_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_localize
!! NAME
!!  polstate_localize
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_localize(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!----------------------------------------------------------------------

 call self%get_b_from_a()

 self%enel = self%get_enel()
 self%enph = self%get_enph()
 self%enelph = self%get_enelph()
 self%eps = self%enel + self%enelph

 ! Get gradient
 call self%get_grad_a()

 ! Save the bare gradient norm
 self%gradres = self%get_norm('grad_a')

end subroutine polstate_localize
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enelph
!! NAME
!!  polstate_get_enelph
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

real(dp) function polstate_get_enelph(self) result(enelph)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_iq, my_pert, my_ik, ik_forw, ib, jb
 complex(dp) :: a_from, a_forw, g_forw, b
!arrays
 integer, allocatable :: kpq_map(:, :)
 real(dp) :: kpt(3), kpq(3)

!----------------------------------------------------------------------

 gqk => self%gqk
 ABI_MALLOC(kpq_map, (6, gqk%my_nq))

 enelph = zero
 do my_ik=1,gqk%my_nk
   kpt(:) = self%my_kpts(:, my_ik)

   !call self%get_mapping('kq', kpt, kpq_map)

   do my_iq=1,gqk%my_nq

     !ik_forw = kpq_map(1, my_iq)

     kpq(:) = kpt(:) + self%my_qpts(:, my_iq)
     ik_forw = self%krank_kpts%get_index(kpq)
     if (ik_forw == -1) cycle

     do ib=1,gqk%nb
       a_from = self%my_a(ib, my_ik)
       do jb=1,gqk%nb
         a_forw = self%a_glob(jb, ik_forw)
         do my_pert=1,gqk%my_npert
           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)
           b = self%my_b(my_pert, my_iq)

           enelph = enelph + real(a_from*conjg(b)*g_forw*conjg(a_forw))
         enddo
       enddo
     enddo

   enddo
 enddo
 call xmpi_sum(enelph, gqk%comm%value, ierr)
 enelph = -two*enelph/(self%nkbz*self%nqbz)

 ABI_FREE(kpq_map)

 end function polstate_get_enelph
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enph
!! NAME
!!  polstate_get_enph
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

real(dp) function polstate_get_enph(self) result(enph)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_iq, my_pert
!----------------------------------------------------------------------

 gqk => self%gqk

 enph = zero
 do my_iq=1,gqk%my_nq
   do my_pert=1,gqk%my_npert
     enph = enph + gqk%my_wnuq(my_pert, my_iq)*abs(self%my_b(my_pert, my_iq))**2
   enddo
 enddo
 call xmpi_sum(enph, gqk%qpt_pert_comm%value, ierr)
 enph = enph/self%nqbz

end function polstate_get_enph
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_enel
!! NAME
!!  polstate_get_enel
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

real(dp) function polstate_get_enel(self) result(enel)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_ik, ik_ibz, ib

!----------------------------------------------------------------------

 gqk => self%gqk

 enel = zero
 do my_ik=1,gqk%my_nk
   ik_ibz = gqk%my_k2ibz(1, my_ik)
   do ib=1,gqk%nb
     enel = enel + self%eig(ib, ik_ibz)*abs(self%my_a(ib, my_ik))**2
   enddo
 enddo
 call xmpi_sum(enel, gqk%kpt_comm%value, ierr)
 enel = enel/self%nkbz

end function polstate_get_enel
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_b_from_a
!! NAME
!!  polstate_get_b_from_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_get_b_from_a(self)

!Arguments ------------------------------------
 class(polstate_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: ierr, my_iq, my_pert, my_ik, ik_forw, ib, jb
 complex(dp) :: a_from, a_forw, g_forw, b_tmp
!arrays
 integer, allocatable :: qpk_map(:, :)
 real(dp) :: qpt(3), kpq(3)

!----------------------------------------------------------------------

 gqk => self%gqk
 ABI_MALLOC(qpk_map, (6, gqk%my_nk))

 do my_iq=1,gqk%my_nq
   qpt(:) = self%my_qpts(:, my_iq)

   !call self%get_mapping('qk', qpt, qpk_map)

   do my_pert=1,gqk%my_npert

     if (gqk%my_wnuq(my_pert, my_iq) == zero) then
       self%my_b(my_pert, my_iq) = zero
       cycle
     endif

     b_tmp = zero
     do my_ik=1,self%gqk%my_nk

       !ik_forw = qpk_map(1, my_ik)

       kpq(:) = qpt(:) + self%my_kpts(:, my_ik)
       ik_forw = self%krank_kpts%get_index(kpq)
       if (ik_forw == -1) cycle

       do ib=1,gqk%nb
         a_from = self%my_a(ib, my_ik)
         do jb=1,gqk%nb
           a_forw = self%a_glob(jb, ik_forw)
           g_forw = gqk%my_g(my_pert, jb, my_iq, ib, my_ik)

           b_tmp = b_tmp + a_from*g_forw*conjg(a_forw)
         enddo
       enddo
     enddo
     call xmpi_sum(b_tmp, gqk%kpt_comm%value, ierr)

     self%my_b(my_pert, my_iq) = b_tmp/(self%nkbz*gqk%my_wnuq(my_pert, my_iq))
   enddo
 enddo

 ABI_FREE(qpk_map)

end subroutine polstate_get_b_from_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_load_a
!! NAME
!!  polstate_load_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_load_a(self, varpeq_a)

!Arguments ------------------------------------
!scalars
 class(polstate_t), intent(inout) :: self
!arrays
 real(dp), intent(in) :: varpeq_a(:,:,:)

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: my_ik, ik_glob, ib
 real(dp) :: anorm
 complex(dp) :: ank

!----------------------------------------------------------------------

 gqk => self%gqk

 do my_ik=1,gqk%my_nk
   ik_glob = gqk%my_kstart + my_ik - 1
   do ib=1,gqk%nb
     ! FIXME: find a better way to represnet vpq%a_spin so I don't have to treat Re & Im part separately
     ank = varpeq_a(1, ib, ik_glob) + j_dpc*varpeq_a(2, ib, ik_glob)
     self%my_a(ib, my_ik) = ank
   enddo
 enddo

 anorm = self%get_norm('a')
 self%my_a(:, :) = self%my_a(:, :)*sqrt(one*self%nkbz)/anorm

 call self%gather('a')

end subroutine polstate_load_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_seed_a
!! NAME
!!  polstate_seed_a
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_seed_a(self, mode, gau_params)

!Arguments ------------------------------------
!scalars
 class(polstate_t), intent(inout) :: self
 character(len=*), intent(in) :: mode
!arrays
 real(dp), optional, intent(in) :: gau_params(2)

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk
 integer :: ierr
 real(dp) :: anorm
 real(dp), allocatable :: re_rand(:,:), im_rand(:,:)

!----------------------------------------------------------------------

 gqk => self%gqk

 select case(mode)
 case ("gaussian")
   ABI_CHECK(present(gau_params), "polstate_seed_a: missing gau_params argument")
   call gaussian_(gau_params(1), gau_params(2))
 case ("random")
   ABI_MALLOC(re_rand, (gqk%nb, gqk%my_nk))
   ABI_MALLOC(im_rand, (gqk%nb, gqk%my_nk))

   call random_number(re_rand)
   call random_number(im_rand)
   self%my_a(:,:) = re_rand(:,:) + j_dpc*im_rand(:,:)

   call xmpi_sum(self%my_a, gqk%qpt_pert_comm%value, ierr)

   ABI_FREE(re_rand)
   ABI_FREE(im_rand)
 case default
   ABI_ERROR(sjoin("polstate_seed_a, unsuported mode: ", mode))
 end select

 anorm = self%get_norm('a')
 self%my_a(:, :) = self%my_a(:, :)*sqrt(one*self%nkbz)/anorm

 call self%gather('a')

!----------------------------------------------------------------------

 contains
 subroutine gaussian_(mu, sigma2)

  real(dp), intent(in) :: mu, sigma2
  character(len=5000) :: msg
  integer :: my_ik, ik_ibz, ib
  real(dp) :: eig

 !----------------------------------------------------------------------

  msg = sjoin("Seeding initial electronic vector with gaussian: &
    variance parameter sigma^2 has to be > 0, but got ", ftoa(sigma2))
  ABI_CHECK(sigma2 > 0, msg)

  do my_ik=1,gqk%my_nk
    ik_ibz = gqk%my_k2ibz(1, my_ik)
    do ib=1,gqk%nb
      eig = self%eig(ib, ik_ibz)
      self%my_a(ib, my_ik) = exp(-(eig - mu)**2/(2*sigma2))
    enddo
  enddo

 end subroutine gaussian_

end subroutine polstate_seed_a
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_norm
!! NAME
!!  polstate_get_norm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

real(dp) function polstate_get_norm(self, mode) result(norm)

!Arguments ------------------------------------
 class(polstate_t), target, intent(inout) :: self
 character(len=*), intent(in) :: mode

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk
 select case(mode)
 case ("a")
   norm = get_norm_(self%my_a, gqk%kpt_comm%value)
 case ("grad_a")
   norm = get_norm_(self%my_grad_a, gqk%kpt_comm%value)
 case ("b")
   norm = get_norm_(self%my_b, gqk%qpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_get_norm, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 real(dp) function get_norm_(my_arr, comm) result(norm)

  integer, intent(in) :: comm
  complex(dp), intent(in) :: my_arr(:, :)

  integer :: ierr
  real(dp) :: norm2_tmp

 !----------------------------------------------------------------------

  norm2_tmp = sum(abs(my_arr(:,:))**2)
  call xmpi_sum(norm2_tmp, comm, ierr)
  norm = sqrt(norm2_tmp)

 end function get_norm_

end function polstate_get_norm
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_gather
!! NAME
!!  polstate_gather
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_gather(self, mode)

!Arguments ------------------------------------
 class(polstate_t), target, intent(inout) :: self
 character(len=*),intent(in) :: mode

!Local variables-------------------------------
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk
 select case(mode)
 case ("a")
   call gather_(self%my_a, gqk%my_nk, gqk%my_kstart, self%a_glob, gqk%kpt_comm%value)
 case ("grad_a")
   call gather_(self%my_grad_a, gqk%my_nk, gqk%my_kstart, self%grad_a_glob, gqk%kpt_comm%value)
 case ("b")
   call gather_(self%my_b, gqk%my_nq, gqk%my_qstart, self%b_glob, gqk%qpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_gather, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 subroutine gather_(my_arr, my_nk, my_kstart, glob_arr, comm)

  integer, intent(in) :: comm, my_nk, my_kstart
  complex(dp), intent(in) :: my_arr(:, :)
  complex(dp), intent(out) :: glob_arr(:, :)

  integer :: ierr, my_ik, ik_glob

 !----------------------------------------------------------------------

  glob_arr(:, :) = zero
  do my_ik=1,my_nk
    ik_glob = my_ik + my_kstart - 1
    glob_arr(:, ik_glob) = my_arr(:, my_ik)
  enddo
  call xmpi_sum(glob_arr, comm, ierr)

 end subroutine gather_

end subroutine polstate_gather
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_mapping
!! NAME
!!  polstate_get_mapping
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine polstate_get_mapping(self, mode, kpt, map)

!Arguments ------------------------------------
!scalars
 class(polstate_t), target, intent(inout) :: self
 character(len=*),intent(in) :: mode
!arrays
 integer, intent(out) :: map(:, :)
 real(dp), intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk
 select case(mode)
 case ("qk")
   call get_mapping_("symrel", gqk%my_nk, self%my_kpts, self%krank_kpts)
 case ("kq")
   call get_mapping_("symrel", gqk%my_nq, self%my_qpts, self%krank_kpts)
 case default
   ABI_ERROR(sjoin("polstate_get_mapping, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 subroutine get_mapping_(symmode, from_nk, from_kpts, to_krank)

  character(len=*),intent(in) :: symmode
  type(krank_t), intent(inout) :: to_krank
  integer, intent(in) :: from_nk
  real(dp), intent(in) :: from_kpts(3, from_nk)

  integer :: ierr, timrev

 !----------------------------------------------------------------------

  timrev = 0
  ierr = kpts_map(symmode, timrev, self%cryst, to_krank, from_nk, from_kpts, &
    map, qpt=kpt)

  if (ierr /= 0) then
     ABI_ERROR(sjoin("polstate_get_mapping: cannot map ", mode, "with point ", ktoa(kpt)))
   endif

 end subroutine get_mapping_

end subroutine polstate_get_mapping
!!***

!----------------------------------------------------------------------

!!****f* m_varpeq/polstate_get_krank_glob
!! NAME
!!  polstate_get_krank_glob
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

type(krank_t) function polstate_get_krank_glob(self, mode, kptrlatt) result(krank_kpts)

!Arguments ------------------------------------
!scalars
 class(polstate_t), target, intent(in) :: self
 character(len=*), intent(in) :: mode
!arrays
 integer, intent(in) :: kptrlatt(3, 3)

!Local variables-------------------------------
!scalars
 class(gqk_t), pointer :: gqk

!----------------------------------------------------------------------

 gqk => self%gqk

 select case(mode)
 case ("k")
   krank_kpts = get_krank_glob_(self%my_kpts, gqk%my_nk, gqk%my_kstart, &
     gqk%glob_nk, gqk%kpt_comm%value)
 case ("q")
   krank_kpts = get_krank_glob_(self%my_qpts, gqk%my_nq, gqk%my_qstart, &
     gqk%glob_nq, gqk%qpt_comm%value)
 case default
   ABI_ERROR(sjoin("polstate_get_krank_glob, unsuported mode: ", mode))
 end select

!----------------------------------------------------------------------

 contains
 type(krank_t) function get_krank_glob_(my_kpts, my_nk, my_kstart, glob_nk, comm) result(krank_kpts)

  integer, intent(in) :: comm, my_nk, my_kstart, glob_nk
  real(dp), intent(in) :: my_kpts(3, my_nk)

  type(krank_t) :: krank_tmp
  integer :: ierr, my_ik, ik_glob
  real(dp) :: kpts(3, glob_nk)

 !----------------------------------------------------------------------

  kpts(:, :) = zero
  do my_ik=1,my_nk
    ik_glob = my_ik + my_kstart - 1
    kpts(:, ik_glob) = my_kpts(:, my_ik)
  enddo
  call xmpi_sum(kpts, comm, ierr)

  krank_tmp = krank_from_kptrlatt(glob_nk, kpts, kptrlatt, compute_invrank=.True.)
  krank_kpts = krank_tmp%copy()
  call krank_tmp%free()

 end function get_krank_glob_

end function polstate_get_krank_glob
!!***

!!****f* m_varpeq/varpeq_plot
!! NAME
!!  varpeq_plot
!!
!! FUNCTION
!! Compute polaron wavefunctions and atomic displacements in the supercell and write results to XSF files
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_plot(wfk0_path, ngfft, dtset, dtfil, cryst, ebands, pawtab, psps, comm)

!Arguments ------------------------------------
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: ngfft(18) !,ngfftf(18)
 type(dataset_type), intent(in) :: dtset
 type(datafiles_type), intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 integer,intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, ndat1 = 1
 integer :: my_rank, ib, irsp, ir, uc_idx, mpw, band, bstart ! nfft, nfftf, mgfft, mgfftf, n1, n2, n3, n4, n5, n6,
 integer :: natom, natom3, nsppol, nspinor, nspden, nkibz, mband, spin, ik, sc_nfft, ebands_timrev, cnt, nproc
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_kq_ibz, istwf_k_ibz, nkpg_k, ierr, nk, spinor, spad, nkbz, ncid
 integer :: nqibz, iq, qtimrev, iq_ibz, isym_q, trev_q, qptopt, uc_iat, sc_iat, nu, nqbz ! icell,
 logical :: isirr_k, isirr_q, has_scell_q
 real(dp) :: cpu_all, wall_all, gflops_all
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(varpeq_t) :: vpq
 type(wfd_t) :: wfd
 type(supercell_type), target :: scell_q, scell_k
 type(krank_t) :: krank_ibz, qrank_ibz
 complex(dp) :: a_nk, bstar_qnu, cphase, c3tmp(3)
!arrays
 integer :: sc_ngfft(18), mapl_k(6), kptrlatt_(3,3), qptrlatt_(3,3)
 integer :: units(2), work_ngfft(18), gmax(3), g0_k(3), mapl_qq(6), g0_q(3), ngqpt(3)
 integer,allocatable :: nband(:,:), wfd_istwfk(:), kg_k(:,:), sc2uc(:), gbound_k(:,:)
 real(dp),parameter :: origin0(3) = zero
 real(dp) :: kk(3), kk_ibz(3), kk_sc(3), qq(3), qq_ibz(3)
 real(dp),allocatable :: scred(:,:), kpg_k(:,:), ug_k(:,:), work(:,:,:,:), pol_rho(:), qibz(:,:)
 real(dp),allocatable :: pheigvec_qibz(:,:,:,:)
 real(dp),allocatable :: displ_cart_qbz(:,:,:,:), pheigvec_qbz(:,:,:,:)  !displ_red_qbz(:,:,:,:), displ_cart_qibz(:,:,:,:),
 real(dp),allocatable :: phfreqs_ibz(:,:), pheigvec_cart_ibz(:,:,:,:,:) !, pheigvec_cart_qbz(:,:,:,:)
 real(dp),allocatable :: sc_displ_cart_re(:,:,:), sc_displ_cart_im(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: xcart_ptr(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 !complex(dp),allocatable :: ceiqr(:)
 complex(gwpc),allocatable :: ur_k(:), pol_wf(:,:), ceikr(:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Read A_nk and B_qnu and other useful tables from file
 call vpq%ncread(dtfil%filvarpeqin, comm, keep_open=.False.)
 call wrtout(std_out, "Reading done")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 nkibz = ebands%nkpt; mband = ebands%mband

 if (dtfil%filgstorein == ABI_NOFILE) then
   call wrtout(units, "gstore_filepath is not specified in input. Cannot compute polaron-induced displacements!")
   has_scell_q = .False.

 else
   ! Start by reading ph displacements and frequencies in the IBZ from the gstore file.
   ! First compute displaced supercell then polaron wf so that we can use both when writing the XSF file.
   call wrtout(units, sjoin(" Computing polaron-induced displacements. Reading phonons from: ", dtfil%filgstorein))
   call cwtime(cpu_all, wall_all, gflops_all, "start")
   has_scell_q = .True.

   NCF_CHECK(nctk_open_read(ncid, dtfil%filgstorein, comm))
   NCF_CHECK(nctk_get_dim(ncid, "gstore_nqibz", nqibz))
   !NCF_CHECK(nctk_get_dim(ncid, "gstore_nqbz", nqbz))

   ! TODO: Wrap phstore API?
   !call gstore_read_ph_qibz(dtfil%filgstorein, ph, comm)
   !call ph%free()

   ABI_MALLOC(qibz, (3, nqibz))
   ABI_MALLOC(phfreqs_ibz, (natom3, nqibz))
   ABI_MALLOC(pheigvec_cart_ibz, (2, 3, cryst%natom, cryst%natom * 3, nqibz))
   if (nproc > 1) then
     NCF_CHECK(nctk_set_collective(ncid, vid("gstore_qibz")))
     NCF_CHECK(nctk_set_collective(ncid, vid("phfreqs_ibz")))
     NCF_CHECK(nctk_set_collective(ncid, vid("pheigvec_cart_ibz")))
   end if
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qibz"), qibz))
   NCF_CHECK(nf90_get_var(ncid, vid("phfreqs_ibz"), phfreqs_ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("pheigvec_cart_ibz"), pheigvec_cart_ibz))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_ngqpt"), ngqpt))
   NCF_CHECK(nf90_get_var(ncid, vid("gstore_qptopt"), qptopt))
   NCF_CHECK(nf90_close(ncid))

   ABI_MALLOC(pheigvec_qibz, (2, 3, cryst%natom, natom3))
   ABI_MALLOC(displ_cart_qbz, (2, 3, cryst%natom, cryst%natom * 3))
   ABI_MALLOC(pheigvec_qbz, (2, 3, cryst%natom, 3*cryst%natom))

   qtimrev = kpts_timrev_from_kptopt(qptopt)
   nqbz = product(ngqpt)
   call kptrlatt_from_ngkpt(ngqpt, qptrlatt_)
   qrank_ibz = krank_from_kptrlatt(nqibz, qibz, qptrlatt_, compute_invrank=.False.)

   call scell_q%init(cryst%natom, qptrlatt_, cryst%rprimd, cryst%typat, cryst%xcart, cryst%znucl, xyz_order="xyz")

   !ABI_MALLOC(ceiqr, (scell_q%ncells))
   ABI_CALLOC(sc_displ_cart_re, (3, scell_q%natom, nsppol))
   ABI_CALLOC(sc_displ_cart_im, (3, scell_q%natom, nsppol))

   ! TODO: Finalize the implementation.
   cnt = 0
   do spin=1,nsppol
     do iq=1,vpq%nq_spin(spin)
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm.

       !if (iq == 1) then
       !  do nu=1,3; write(std_out, *)"hello:", vpq%b_spin(1, nu, iq, spin) + j_dpc * vpq%b_spin(2, nu, iq, spin); end do
       !end if

       qq = vpq%qpts_spin(:, iq, spin)
       ! Note symrec here
       if (kpts_map("symrec", qtimrev, cryst, qrank_ibz, 1, qq, mapl_qq) /= 0) then
         ABI_ERROR("Cannot map qBZ to IBZ!")
       end if
       iq_ibz = mapl_qq(1); isym_q = mapl_qq(2)
       trev_q = mapl_qq(6); g0_q = mapl_qq(3:5)
       ! Don't test if umklapp == 0 because we use the periodic gauge:
       !
       !      phfreq(q+G) = phfreq(q) and eigvec(q) = eigvec(q+G)
       !
       isirr_q = (isym_q == 1 .and. trev_q == 0)
       qq_ibz = qibz(:, iq_ibz)
       !if (all(abs(qq_ibz) < tol6)) cycle
       pheigvec_qibz = pheigvec_cart_ibz(:,:,:,:,iq_ibz)

       if (isirr_q) then
         ! Compute phonon displacements in Cartesian coordinates
         call phdispl_from_eigvec(cryst%natom, cryst%ntypat, cryst%typat, cryst%amu, pheigvec_qibz, displ_cart_qbz)

       else
         ! Rotate phonon eigenvectors from q_ibz to q_bz.
         ! This part is needed to enforce the gauge in the ph eigenvectors, including e(-q) = e(q)^*
         call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_qibz, pheigvec_qbz, displ_cart_qbz)
       end if

       ! Compute phase e^{iq.R}
       !do icell=1,scell_q%ncells
       !  ceiqr(icell) = exp(+j_dpc * two_pi * dot_product(qq, scell_q%rvecs(:, icell)))
       !end do

       do sc_iat=1, scell_q%natom
         uc_iat = scell_q%atom_indexing(sc_iat)
         cphase = exp(+j_dpc * two_pi * dot_product(qq, scell_q%uc_indexing(:, sc_iat)))
         ! Summing over ph modes.
         do nu=1,natom3
           bstar_qnu = conjg(vpq%b_spin(1, nu, iq, spin) + j_dpc * vpq%b_spin(2, nu, iq, spin))
           c3tmp = (displ_cart_qbz(1,:,uc_iat,nu) + j_dpc * displ_cart_qbz(2,:,uc_iat,nu)) * bstar_qnu * cphase
           sc_displ_cart_re(:,sc_iat,spin) = sc_displ_cart_re(:,sc_iat,spin) + real(c3tmp)
           sc_displ_cart_im(:,sc_iat,spin) = sc_displ_cart_im(:,sc_iat,spin) + aimag(c3tmp) ! This to check if the imag part is zero.
         end do
       end do ! sc_iat

     end do ! iq
   end do ! spin

   ABI_FREE(qibz)
   ABI_FREE(phfreqs_ibz)
   ABI_FREE(pheigvec_cart_ibz)
   ABI_FREE(displ_cart_qbz)
   ABI_FREE(pheigvec_qbz)
   ABI_FREE(pheigvec_qibz)
   !ABI_FREE(ceiqr)
   call qrank_ibz%free()

   call xmpi_sum_master(sc_displ_cart_re, master, comm, ierr)
   call xmpi_sum_master(sc_displ_cart_im, master, comm, ierr)
   sc_displ_cart_re = -sqrt2 * sc_displ_cart_re / nqbz
   sc_displ_cart_im = -sqrt2 * sc_displ_cart_im / nqbz

   ! Here we displace the atoms in the supercell for this spin (only master has the correct values)
   !do sc_iat=1,scell_q%natom; write(std_out, *) sc_displ_cart_re(:,sc_iat); end do
   scell_q%xcart = scell_q%xcart_ref + sc_displ_cart_re(:,:,spin)

   ! Write polaron-induced displacements in XSF format.
   if (my_rank == master) then
     ! Handle spin-polarized case by writing two XSF files.
     do spin=1,nsppol
       !write(msg, "(a,i0,a,es16.6)")" For spin: ", spin, ": 1/N_q \sum_qnu |B_qnu|^2 = ", sum(vpq%b_spin(:,:,:, spin) ** 2) / nqbz
       !call wrtout(units, msg)
       write(msg, "(a,i0,a,es16.6)") " For spin: ", spin, ": maxval(abs(sc_displ_cart_im)): ", maxval(abs(sc_displ_cart_im(:,:,spin)))
       call wrtout(units, msg)
       path = strcat(dtfil%filnam_ds(4), "_POLARON_DISPL.xsf")
       if (nsppol == 2) path = strcat(dtfil%filnam_ds(4), "_spin_", itoa(spin), "_POLARON_DISPL.xsf")
       call scell_q%write_xsf(path)
     end do
   end if
   ABI_FREE(sc_displ_cart_re)
   ABI_FREE(sc_displ_cart_im)
   call cwtime_report(" Computation of polaron-induced displacements completed:", cpu_all, wall_all, gflops_all, &
                      pre_str=ch10, end_str=ch10)
 end if

 call wrtout(std_out, " varpeq_plot: computing polaron wavefunction in real space.", pre_newlines=1)

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 krank_ibz = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 ! Initialize the wave function descriptor.
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkibz, nsppol))
 nband = mband; bks_mask = .False.; keep_ur = .False.

 ! Here we use brange_spin to select the number of bands that should be read and stored in memory.
 ! For the time being, spin and k-points are not MPI-distributed inside comm.
 do spin=1,nsppol
   bstart = vpq%brange_spin(1, spin)
   do ik=1, vpq%nk_spin(spin)
     kk = vpq%kpts_spin(:, ik, spin)
     ! Note symrel option
     if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, 1, kk, mapl_k) /= 0) then
       write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k-point could not be generated from a symmetrical one.",trim(ltoa(kk))
       ABI_ERROR(msg)
     end if

     ik_ibz = mapl_k(1)
     do ib=1,vpq%nb_spin(spin)
       band = bstart + ib - 1
       bks_mask(band, ik_ibz, spin) = .True.
     end do
   end do
 end do

 ! Impose istwfk = 1 for all k-points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkibz, nsppol, bks_mask,&
               dtset%nspden, dtset%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft, &
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print(header="Wavefunctions for varpeq_plot")

 if (dtset%boxcutmin >= two) then
   call wrtout(std_out, " To reduce the size of the FFT mesh and the size of the XSF file, reduce boxcutmin from 2 to e.g. 1.1")
 end if

 ABI_FREE(nband)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(bks_mask)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 call wrtout(std_out, " Computing mpw and gmax needed to allocate workspace array.")
 do spin=1,nsppol
   nk = vpq%nk_spin(spin)
   do ik=1,nk
     call ephtk_get_mpw_gmax(nk, vpq%kpts_spin(:, 1:nk, spin), dtset%ecut, cryst%gmet, mpw, gmax, comm, &
                             init_with_zero=spin==1)
   end do
 end do

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 call wrtout(std_out, sjoin(" Building supercell from ngkpt:", ltoa(vpq%ngkpt)))
 ! (note xyz_order)
 call kptrlatt_from_ngkpt(vpq%ngkpt, kptrlatt_)
 call scell_k%init(cryst%natom, kptrlatt_, cryst%rprimd, cryst%typat, cryst%xcart, cryst%znucl, xyz_order="xyz")

 nkbz = product(vpq%ngkpt)
 call supercell_fft(vpq%ngkpt, ngfft, sc_nfft, sc_ngfft, sc2uc, scred)
 ABI_FREE(scred)

 call wrtout(std_out, " Computing polaron wavefunction in the real-space supercell...")

 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(gbound_k, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(ur_k, (wfd%nfft*nspinor))
 ABI_MALLOC(ceikr, (sc_nfft*nspinor))
 ABI_CALLOC(pol_wf, (sc_nfft*nspinor, nsppol)) ! Init output with zeros.

 cnt = 0
 do spin=1,nsppol
   bstart = vpq%brange_spin(1, spin)
   nk = vpq%nk_spin(spin)
   do ik=1, nk
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism inside comm.

     kk = vpq%kpts_spin(:, ik, spin)
     if (kpts_map("symrel", ebands_timrev, cryst, krank_ibz, 1, kk, mapl_k) /= 0) then
       write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k-point could not be generated from a symmetrical one.",trim(ltoa(kk))
       ABI_ERROR(msg)
     end if

     ik_ibz = mapl_k(1); isym_k = mapl_k(2); trev_k = mapl_k(6); g0_k = mapl_k(3:5)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     kk_ibz = ebands%kptns(:, ik_ibz)
     istwf_k_ibz = wfd%istwfk(ik_ibz); npw_kq_ibz = wfd%npwarr(ik_ibz)

     ! Get npw_k, kg_k for k
     call wfd%get_gvec_gbound(cryst%gmet, dtset%ecut, kk, ik_ibz, isirr_k, dtset%nloalg, &  ! in
                              istwf_k, npw_k, kg_k, nkpg_k, kpg_k, gbound_k)                ! out

     ABI_MALLOC(ug_k, (2, npw_k*nspinor))
     kk_sc = kk * vpq%ngkpt
     call calc_ceikr(kk_sc, sc_ngfft, sc_nfft, nspinor, ceikr)

     do ib=1,vpq%nb_spin(spin)
       band = bstart + ib - 1
       call wfd%rotate_cg(band, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                          cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_k)
       !print *, "int_omega dr |u(r)}^2:", sum(abs(ur_k) ** 2) / wfd%nfft

       a_nk = vpq%a_spin(1, ib, ik, spin) + j_dpc * vpq%a_spin(2, ib, ik, spin)
       do spinor=1,nspinor
         spad = (spinor - 1) * sc_nfft
         do ir=1,sc_nfft
           uc_idx = sc2uc(ir) + (spinor - 1) * wfd%nfft
           irsp = ir + spad
           pol_wf(irsp, spin) = pol_wf(irsp, spin) + a_nk * ur_k(uc_idx) * ceikr(irsp)
         end do
       end do
     end do ! ib

     ABI_FREE(ug_k)
     ABI_FREE(kpg_k)
   end do ! ik
 end do ! spin

 ABI_FREE(kg_k)
 ABI_FREE(work)
 ABI_FREE(ur_k)
 ABI_FREE(gbound_k)
 ABI_FREE(sc2uc)
 ABI_FREE(ceikr)
 call krank_ibz%free()

 ! Collect pol_wf on the master rank who's gonna write the polaron density in XSF format.
 call xmpi_sum_master(pol_wf, master, comm, ierr)

 if (my_rank == master) then
   pol_wf = pol_wf / (nkbz * sqrt(cryst%ucvol))
   ABI_MALLOC(pol_rho, (sc_nfft))

   ! Here decide if we write the polaron wavefunction with_diplaced atoms or not.
   if (all(kptrlatt_ == qptrlatt_) .and. has_scell_q) then
     call wrtout(units, " Writing the polaron wavefunction with diplaced atoms")
     xcart_ptr => scell_q%xcart
   else
     call wrtout(units, " Writing the polaron wavefunction with undiplaced atoms")
     xcart_ptr => scell_k%xcart
   end if

   if (nspinor == 1) then
     ! Handle spin-polarized case by writing two XSF files.
     do spin=1,nsppol
       write(msg, "(a,i0,a,es16.6)")&
         " For spin: ", spin, ": 1/N_k \sum_nk |A_nk|^2 = ", sum(vpq%a_spin(:,:,:,spin)**2) / nkbz
       call wrtout(units, msg)
       pol_rho = abs(pol_wf(:, spin))**2
       write(msg, "(a,i0,a,es16.6)")" Polaron density for spin: ", spin, " integrates to: ", sum(pol_rho) * cryst%ucvol/wfd%nfft
       call wrtout(units, msg)
       write(msg, "(a,es16.6)")" maxval(abs(aimag(pol_wf))): ", maxval(abs(aimag(pol_wf(:, spin))))
       call wrtout(units, msg)
       path = strcat(dtfil%filnam_ds(4), "_POLARON.xsf")
       if (nsppol == 2) path = strcat(dtfil%filnam_ds(4), "_spin_", itoa(spin), "_POLARON.xsf")
       call write_xsf(path, &
                      sc_ngfft(1), sc_ngfft(2), sc_ngfft(3), pol_rho, scell_k%rprimd, origin0, &
                      scell_k%natom, scell_k%ntypat, scell_k%typat, xcart_ptr, scell_k%znucl, 0)
     end do ! spin

   else
     write(msg, "(a,es16.6)")" 1/N_k \sum_nk |A_nk|^2 = ", sum(vpq%a_spin(:,:,:, spin)**2) / nkbz
     call wrtout(units, msg)

     pol_rho(:) = abs(pol_wf(1:sc_nfft, 1)) ** 2
     pol_rho(:) = abs(pol_wf(sc_nfft+1:, 1)) ** 2 + pol_rho(:)
     write(msg, "(a,i0,a,es16.6)")" Polaron density integrates to: ", sum(pol_rho) * cryst%ucvol/wfd%nfft
     call wrtout(units, msg)
     write(msg, "(a,es16.6)")" maxval(abs(aimag(pol_wf))): ", maxval(abs(aimag(pol_wf(:, 1))))
     call wrtout(units, msg)

     call write_xsf(strcat(dtfil%filnam_ds(4), "_POLARON.xsf"), &
                    sc_ngfft(1), sc_ngfft(2), sc_ngfft(3), pol_rho, scell_k%rprimd, origin0, &
                    scell_k%natom, scell_k%ntypat, scell_k%typat, xcart_ptr, scell_k%znucl, 0)
   end if
   ABI_FREE(pol_rho)
 end if ! master

 call cwtime_report(" Computation of polaron wavefunction completed", cpu_all, wall_all, gflops_all, pre_str=ch10, end_str=ch10)

 ABI_FREE(pol_wf)
 call wfd%free()
 call scell_q%free()
 call scell_k%free()
 call vpq%free()

contains
integer function vid(var_name)
  character(len=*),intent(in) :: var_name
  vid = nctk_idname(ncid, var_name)
end function vid

end subroutine varpeq_plot
!!***

end module m_varpeq
!!***
