!!****m* ABINIT/m_rmm_diis
!! NAME
!!  m_rmm_diis
!!
!! FUNCTION
!!  This module contains routines for the RMM-DIIS eigenvalue solver.
!!
!! COPYRIGHT
!!  Copyright (C) 2020-2020 ABINIT group (MG)
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

module m_rmm_diis

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_dtset
 use m_cgtools
 use m_hide_blas
 use m_yaml
 use m_linalg_interfaces
 use m_prep_kgb

 use defs_abitypes,   only : mpi_type
 use m_fstrings,      only : sjoin, itoa, ftoa
 use m_time,          only : timab, cwtime, cwtime_report
 use m_numeric_tools, only : pack_matrix, imin_loc, stats_t, stats_eval
 use m_hide_lapack,   only : xhegv_cplex, xhesv_cplex
 use m_pair_list,     only : pair_list
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_fftcore,       only : fftcore_set_mixprec
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 use m_nonlop,        only : nonlop
 use m_cgtk,          only : cgtk_fixphase
 use m_abi_linalg,    only : abi_zgemm_2r
 !use m_fock,         only : fock_set_ieigen, fock_set_getghc_call

 implicit none

 private
!!***

 public :: rmm_diis
!!***

 type,private :: rmm_diis_t

   integer :: accuracy_level
   ! Defines tolerances, activates/deactivates tricks.

   integer :: usepaw
   ! 1 if we are running PAW.

   integer :: istwf_k
   ! wavefunction storage mode.

   integer :: cplex
   ! 1 if matrices are real (e.g. Gamma-point), 2 for complex

   integer :: bsize
   ! (Max) block size for bands

   integer :: max_niter
   ! Maximum number of iterations

   integer :: npwsp
   ! Total number of planewaves treated by this proc
   ! npw * my_nspinor

   integer :: prtvol
   ! vervosity level

   integer :: last_iter
   ! Last RMM-DIIS iteration performed.

   real(dp) :: tol_occupied
   ! Tolerance for partial occupied states

   type(pair_list) :: stats

   real(dp),allocatable :: hist_ene(:,:)
   real(dp),allocatable :: hist_resid(:,:)
   real(dp),allocatable :: hist_enlx(:,:)
   character(len=7),allocatable :: step_type(:,:)
   ! (0:max_niter+2, bsize)
   ! 0 is the initial step, then DIIS iterations whose number may depend on the block
   ! followed by the computation of eigens after ortho.

   real(dp),allocatable :: resmat(:,:,:,:)
   ! (2, 0:max_niter, 0:max_niter, bsize))

   real(dp),allocatable :: chain_phi(:,:,:,:)
   real(dp),allocatable :: chain_sphi(:,:,:,:)
   real(dp),allocatable :: chain_resv(:,:,:,:)
   ! (2, npwsp, 0:max_niter, bsize))

 contains
   procedure :: free => rmm_diis_free                   ! Free dynamic memory
   procedure :: update_block => rmm_diis_update_block   ! DIIS uppdate of wavefuntions and residuals.
   procedure :: eval_mats => rmm_diis_eval_mats         ! Compute DIIS matrices
   procedure :: exit_iter => rmm_diis_exit_iter         ! Return True if can exit the DIIS iteration.
   procedure :: print_block => rmm_diis_print_block     ! Print energies, residuals and diffs for a given block.
   ! TODO: Fix problem with last_iter and hist
   procedure :: push_iter => rmm_diis_push_iter         ! Save results required the by DIIS algorithm

 end type rmm_diis_t

 integer,parameter, private :: level = 432
 logical,parameter, private :: timeit = .False.
 !logical,parameter, private :: timeit = .True.

contains
!!***

!!****f* ABINIT/rmm_diis
!! NAME
!! rmm_diis
!!
!! FUNCTION
!! This routine updates the wave functions at a given (k-point, spin), using the RMM-DIIS method.
!!
!! INPUTS
!!  istep,ikpt,isppol=Iteration step, k-point index, spin index (mainly for printing purposes).
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  my_nspinor=number of plane waves at this k point
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!  If usepaw==1:
!!    gsc(2,*)=<g|s|c> matrix elements (s=overlap)
!!  If usepaw==0
!!    enlx(nband)=contribution from each band to nonlocal psp + potential Fock ACE part
!!                of total energy, at this k-point
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  resid(nband)=residuals for each states. In input: previous residuals for this k-point, spin.
!!   In output: new residuals.
!!  rmm_diis_status(2): Status of the eigensolver.
!!    The first entry gives the previous accuracy.
!!    The second entry gives the number of iterations already performed with this level.
!!
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine rmm_diis(istep, ikpt, isppol, cg, dtset, eig, occ, enlx, gs_hamk, kinpw, gsc, &
                    mpi_enreg, nband, npw, my_nspinor, resid, rmm_diis_status)

!Arguments ------------------------------------
 integer,intent(in) :: istep, ikpt, isppol, nband, npw, my_nspinor
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*my_nspinor*nband)
 real(dp),target,intent(inout) :: gsc(2,npw*my_nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband), resid(nband)
 real(dp),intent(in) :: occ(nband), kinpw(npw)
 real(dp),intent(out) :: eig(nband)
 integer,intent(inout) :: rmm_diis_status(2)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer,parameter :: choice1 = 1, signs1 = 1, signs2 = 2, tim_nonlop = 0, paw_opt0 = 0, paw_opt3 = 3
 integer :: ierr, prtvol, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, idat, paral_kgb !, ortalgo
 integer :: cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, optekin, usepaw, iter, max_niter, max_niter_block
 integer :: me_g0, nb_pocc, jj, kk, accuracy_level, raise_acc, prev_mixprec, after_ortho, me_cell
 integer :: comm_bsf, prev_accuracy_level, ncalls_with_prev_accuracy, signs, paw_opt, savemem
 logical :: first_call, use_fft_mixprec, has_fock
 real(dp),parameter :: rdummy = zero
 real(dp) :: accuracy_ene,  max_res_pocc, tol_occupied, lock_tolwfr !, max_absimag
 real(dp) :: cpu, wall, gflops, cpu_all, wall_all, gflops_all
 character(len=500) :: msg
 type(yamldoc_t) :: rmm_ydoc
 type(rmm_diis_t) :: diis
 type(stats_t) :: res_stats
!arrays
 real(dp) :: tsec(2)
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp),allocatable :: lambda_bk(:), kres_bk(:,:), dots_bk(:,:), residv_bk(:,:)
 real(dp),allocatable :: umat(:,:,:), gwork(:,:), dots(:, :)
 real(dp),target,allocatable :: ghc(:,:), gvnlxc(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_bk(:,:), cg_bk(:,:), ghc_bk(:,:), gvnlxc_bk(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 ! Define useful vars.
 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; paral_kgb = mpi_enreg%paral_kgb
 me_g0 = mpi_enreg%me_g0; comm_bsf = mpi_enreg%comm_bandspinorfft
 npwsp = npw * my_nspinor; mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 me_cell = mpi_enreg%me_cell; prtvol = dtset%prtvol !; prtvol = -level
 has_fock = associated(gs_hamk%fockcommon)

 if (timeit) then
   call cwtime(cpu_all, wall_all, gflops_all, "start")
   call cwtime(cpu, wall, gflops, "start")
 end if

 ! =================
 ! Prepare DIIS loop
 ! =================
 ! accuracy_level is computed from the maxval of the previous residuals received in input.
 ! The different levels are:
 !
 !  1: Used at the beginning of the SCF cycle. Use loosy convergence criteria in order
 !     to reduce the number of H|psi> applications as much as possible so that
 !     we can start to mix densities/potentials.
 !     Allow for incosistent data in output rediduals and Vnl matrix elements.
 !     Move to the next level after Max 15 iterations.
 !
 !  2: Intermediate step. Decrease convergence criteria in order to perform more wavefuction iterations.
 !     Allow for incosistent data in output rediduals and Vnl matrix elements.
 !     Move to the next level after Max 25 iterations.
 !
 !  3: Approaching convergence. Use stricter convergence criteria.
 !     Move to the next level after Max 25 iterations.
 !
 !  4: Ultimate precision. Try to reach the same accuracy as the other eigenvalue solvers.
 !     This implies: using similar convergence criteria as in the other solvers.

 ! Note:
 !
 ! * Accuracy_level is not allowed to increase during the SCF cycle
 !
 ! * Since we operate on blocks of bands, all the states in the block will receive the same treatment.
 !   This means that one can observe different convergence behaviour depending on bsize.
 !
 if (all(rmm_diis_status == 0)) then
   ! This is the first time we call rmm_diis for this (k-point, spin)
   prev_accuracy_level = 1; ncalls_with_prev_accuracy = 0
   first_call = .True.
 else
   prev_accuracy_level = rmm_diis_status(1); ncalls_with_prev_accuracy = rmm_diis_status(2)
   first_call = .False.
 end if

 ! Decide whether we should move to the next level.
 raise_acc = 0
 if (prev_accuracy_level == 1 .and. ncalls_with_prev_accuracy >= 15) raise_acc = 2
 if (prev_accuracy_level == 2 .and. ncalls_with_prev_accuracy >= 25) raise_acc = 3
 if (prev_accuracy_level == 3 .and. ncalls_with_prev_accuracy >= 25) raise_acc = 4
 if (raise_acc > 0) then
   MSG_COMMENT("Accuracy_level is automatically increased as we reached the max number of NSCF iterations.")
 end if
 raise_acc = max(raise_acc, prev_accuracy_level)

 ! Define tolerance for occupied states on the basis of prev_accuracy_level
 ! and compute max of residuals for these bands.
 tol_occupied = zero
 if (dtset%iscf > 0) then
   tol_occupied = tol3; if (any(prev_accuracy_level == [1])) tol_occupied = tol2
 end if
 nb_pocc = count(occ >= tol_occupied)
 max_res_pocc = maxval(resid(1:nb_pocc))

 ! Define accuracy_level for this run.
 accuracy_level = 1
 if (max_res_pocc < tol8)  accuracy_level = 2
 if (max_res_pocc < tol12) accuracy_level = 3
 if (max_res_pocc < tol16) accuracy_level = 4
 !if (max_res_pocc < tol18) accuracy_level = 4
 accuracy_level = max(prev_accuracy_level, accuracy_level, raise_acc)
 if (istep == 1) accuracy_level = 2  ! FIXME: Differenciate between restart or rmm_diis - 3.
 if (first_call .and. max_res_pocc == zero) accuracy_level = 1
 !print *, "rmm_diis_status:", rmm_diis_status
 !print *, "rmm_prev_acc:", prev_accuracy_level, "rmm_raise_acc:", raise_acc
 !print *, "accuracy_level:", accuracy_level, "rmm_raise_acc:", raise_acc

 ! Update rmm_diis_status. Reset number of calls if we've just moved to a new accuracy_level.
 rmm_diis_status(1) = accuracy_level
 if (accuracy_level /= prev_accuracy_level) rmm_diis_status(2) = 0
 rmm_diis_status(2) = rmm_diis_status(2) + 1

 ! Will perform max_niter DIIS steps. Usually 3 as nline by default is 4.
 ! Note that, unlike in Vasp's recipe, here we don't end with a trial step after DIIS.
 max_niter = max(dtset%nline - 1, 1)
 if (accuracy_level >= 4) max_niter = dtset%nline
 if (dtset%iscf < 0) max_niter = dtset%nline + 1

 ! Define accuracy_ene for SCF.
 accuracy_ene = zero
 if (dtset%iscf > 0) then
   if (dtset%toldfe /= zero) then
     accuracy_ene = dtset%toldfe * ten**(-accuracy_level + 2) / nb_pocc
   else
     ! We are not using toldfe to stop the SCF cycle
     ! so we are forced to hardcode a tolerance for the absolute diff in the KS eigenvalue.
     accuracy_ene = tol8 * ten**(-accuracy_level + 2) / nb_pocc
   end if
 end if

 ! Tolerance on residuals used for band locking after subdiago.
 if (dtset%tolwfr > zero) then
   lock_tolwfr = tol2 * dtset%tolwfr
 else
   lock_tolwfr = tol14
   if (accuracy_level >= 2) lock_tolwfr = tol16
   if (accuracy_level >= 3) lock_tolwfr = tol18
   if (accuracy_level >= 4) lock_tolwfr = tol20 * tol2
 end if

 ! Use mixed precisions if requested by the user but only for low accuracy_level
 use_fft_mixprec = dtset%mixprec == 1 .and. accuracy_level < 2
 if (use_fft_mixprec) prev_mixprec = fftcore_set_mixprec(1)

 ! Select preconditioning.
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1 ! optekin = 0

 ! Will treat states in groups of bsize bands even when paral_kgb = 0
 bsize = 8; if (paral_kgb == 1) bsize = mpi_enreg%nproc_band * mpi_enreg%bandpp
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 ! Build DIIS object.
 diis = rmm_diis_new(accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol)
 diis%tol_occupied = tol_occupied
 !call wrtout(std_out, sjoin(" Using Max", itoa(max_niter), "RMM-DIIS iterations"))
 !call wrtout(std_out, sjoin( &
 !  " Max_input_resid_pocc", ftoa(max_res_pocc), "accuracy_level:", itoa(accuracy_level), &
 !  ", accuracy_ene: ", ftoa(accuracy_ene)))
 call timab(1634, 1, tsec) ! "rmm_diis:band_opt"

 rmm_ydoc = yamldoc_open("RMM-DIIS", with_iter_state=.False.)
 call rmm_ydoc%add_ints("ikpt, isppol, istep, accuracy_level", [ikpt, isppol, istep, accuracy_level])
 call rmm_ydoc%open_tabular("RESIDS_POCC") !, tag, indent, newline, comment)
 write(msg, "(1x, a12, 4(a10))")"level", "mean", "min", "max", "stdev"
 call rmm_ydoc%add_tabular_line(msg, indent=0)
 call rmm_ydoc%add_tabular_line(resids2str("input"), indent=0)

 ! =========================
 ! === Subspace rotation ===
 ! =========================
 ! Allocate big (scalable) array with <G|H|C> for all nband so that we can recompute the residuals after the rotation.
 ! This approach requires more memory but we avoid one extra call to H|Psi> per band.
 ! Alternatively, one can compute ghc and the residuals by applying H|psi>
 ! inside the loop over blocks (less memory but slower).
 savemem = dtset%rmm_diis_savemem
 !savemem = 1
 !if (savemem == 0) then
 !  ABI_MALLOC_OR_DIE(ghc, (2, npwsp*nband), ierr)
 !  ABI_MALLOC_OR_DIE(gvnlxc, (2, npwsp*nband), ierr)
 !end if

 call subspace_rotation(gs_hamk, dtset%prtvol, mpi_enreg, nband, npw, my_nspinor, savemem, &
                        enlx, eig, cg, gsc, ghc, gvnlxc)

 gsc_bk => fake_gsc_bk
 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ABI_MALLOC(lambda_bk, (bsize))
 ABI_MALLOC(dots_bk, (2, bsize))
 ABI_MALLOC(residv_bk, (2, npwsp*bsize))
 ABI_MALLOC(kres_bk, (2, npwsp*bsize))

 if (savemem == 1) then
   ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
   ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))
 end if
 !write(msg, "(a,f8.1,a)") &
 !  " Memory required: ", 2 * natom3**2 * (my_q2 - my_q1 + 1) * dp * b2Mb, " [Mb] <<< MEM"
 !call wrtout(std_out, msg)


 ! We loop over nblocks, each block contains ndat states.
 !
 ! - Convergence behaviour may depend on bsize as branches are taken according to
 !   the status of all bands in the block.
 ! TODO: Transpose only once per block and then work with already_transposed = .True.
 if (timeit) call cwtime(cpu, wall, gflops, "start")

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

   ! Reduce number of niter iterations if block contains "empty" states.
   ! This should happen only if npband is small wrt nband and nband >> nbocc.
   ! TODO: Don't reduce niter if MD
   max_niter_block = max_niter
   if (dtset%iscf > 0) then
     if (all(occ(ib_start:ib_stop) < diis%tol_occupied)) max_niter_block = max(1 + max_niter / 2, 2)
   end if

   ! Compute H |phi_0> with cg block after subdiago.
   cg_bk => cg(:,igs:ige); if (usepaw == 1) gsc_bk => gsc(1:2,igs:ige)

   ! Compute residual vectors after subspace_rotation.
   if (savemem == 0) then
     ghc_bk => ghc(:,igs:ige); gvnlxc_bk => gvnlxc(:,igs:ige)
     call cg_get_residvecs(usepaw, npwsp, ndat, eig(ib_start), cg_bk, ghc_bk, gsc_bk, residv_bk)
     call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, resid(ib_start), me_g0, comm_bsf)
   else
     call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, cg_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                         eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, normalize=.False.)
   end if

   ! Band locking.
   if (all(resid(ib_start:ib_stop) < lock_tolwfr)) then
     call diis%stats%increment("locked", ndat)
     cycle ! iblock
   end if

   ! Save <R0|R0> and <phi_0|S|phi_0>, |phi_0>, |S phi_0>. Assume input cg_bk is already S-normalized.
   call diis%push_iter(0, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), cg_bk, residv_bk, gsc_bk, "SDIAG")

   ! Line minimization with preconditioned steepest descent:
   !
   !    |phi_1> = |phi_0> + lambda |K R_0>
   !
   ! where lambda minimizes the residual (we don't try to find the stationary
   ! point of the Rayleigh quotient as in Kresse's paper).
   !
   !    lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
   !
   ! more expensive than finding the stationary point of the Rayleigh quotient as it requires
   ! an extra H application but it should be more stable and more consistent with the RMM approach.
   !
   ! Precondition |R_0>, output in kres_bk = |K R_0>
   call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
   call cg_precon_many(istwf_k, npw, my_nspinor, ndat, cg_bk, optekin, kinpw, kres_bk, me_g0, comm_bsf)

   ! Compute H |K R_0>
   if (paral_kgb == 0) then
     call getghc(cpopt, kres_bk, cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                 rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
   else
     call prep_getghc(kres_bk, gs_hamk, gvnlxc_bk, ghc_bk, gsc_bk, rdummy, ndat, &
                      mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
   end if

   ! Compute (H - e_0 S) |K R_0>
   call cg_get_residvecs(usepaw, npwsp, ndat, eig(ib_start), kres_bk, ghc_bk, gsc_bk, residv_bk)
   call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_bsf)

   ! Compute lambda
   dots_bk = zero
   do idat=1,ndat
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     call dotprod_g(dots_bk(1,idat), dots_bk(2,idat), istwf_k, npwsp, option1, &
                    diis%chain_resv(:,:,0,idat), residv_bk(:,jj), me_g0, xmpi_comm_self)
   end do
   call xmpi_sum(dots_bk, comm_bsf, ierr)

   ! Build |Psi_1> = |Phi_0> + lambda |K R_0>
   do idat=1,ndat
     lambda_bk(idat) = -dots_bk(1,idat) / lambda_bk(idat)
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     cg_bk(:,jj:kk) = diis%chain_phi(:,:,0,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
   end do

   ! ===============
   ! DIIS iterations
   ! ===============
   iter_loop: do iter=1,max_niter_block

     if (iter > 1) then
       ! Solve DIIS equations and update cg_bk and residv_bk for iter > 1
       call diis%update_block(iter, npwsp, ndat, cg_bk, residv_bk, comm_bsf)

       ! Precondition residual, output in kres_bk.
       call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
       call cg_precon_many(istwf_k, npw, my_nspinor, ndat, cg_bk, optekin, kinpw, kres_bk, me_g0, comm_bsf)

       ! Compute cg_bk with the same lambda(ndat) obtained at iteration #0
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, cg_bk)
     end if

     ! Compute H |phi_now> and evaluate new enlx for NC.
     call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, cg_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, normalize=.True.)

     ! Store new wavevefunctions and residuals.
     call diis%push_iter(iter, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), cg_bk, residv_bk, gsc_bk, "DIIS")

     ! CHECK FOR CONVERGENCE
     if (diis%exit_iter(iter, ndat, max_niter_block, occ(ib_start), accuracy_ene, dtset, comm_bsf)) exit iter_loop

     ! Compute <R_i|R_j> and <i|S|j> for j=iter
     if (iter /= max_niter_block) call diis%eval_mats(iter, ndat, me_g0, comm_bsf)
   end do iter_loop

   if (prtvol == -level) call diis%print_block(ib_start, ndat, istep, ikpt, isppol)
 end do ! iblock

 call timab(1634, 2, tsec) !"rmm_diis:band_opt"
 if (timeit) call cwtime_report(" rmm_diis:band_opt", cpu, wall, gflops)
 call rmm_ydoc%add_tabular_line(resids2str("rmm-diis"), indent=0)

 ! ===============================
 ! Orthogonalize states after DIIS
 ! ===============================
 call timab(583,1,tsec) ! "vtowfk(pw_orthon)"

 !ortalgo = 3 !; ortalgo = mpi_enreg%paral_kgb
 !call pw_orthon(0, 0, istwf_k, mcg, mgsc, npwsp, nband, ortalgo, gsc, usepaw, cg, me_g0, comm_bsf)

 ! TODO: Merge the two routines.
 if (usepaw == 1) then
   !call cgtk_fixphase(cg, gsc, 0, 0, istwf_k, mcg, mgsc, mpi_enreg, nband, npwsp, usepaw)
   !call cgpaw_normalize(npwsp, nband, cg, gsc, istwf_k, me_g0, comm_bsf)

   call cgpaw_cholesky(npwsp, nband, cg, gsc, istwf_k, me_g0, comm_bsf, umat=umat)

   !call cgtk_fixphase(cg, gsc, 0, 0, istwf_k, mcg, mgsc, mpi_enreg, nband, npwsp, usepaw)
   !call cgpaw_normalize(npwsp, nband, cg, gsc, istwf_k, me_g0, comm_bsf)
 else
   call cgnc_cholesky(npwsp, nband, cg, istwf_k, me_g0, comm_bsf, use_gemm=.False., umat=umat)
 end if

 call timab(583,2,tsec)
 if (timeit) call cwtime_report(" pw_orthon ", cpu, wall, gflops)

 ! Recompute eigenvalues, residuals, and enlx after orthogonalization.
 ! This step is important to improve the convergence of the NC total energy
 ! and it guarantees that eigenvalues and residuals are consistent with the output wavefunctions.
 ! but we try to avoid it at the beginning of the SCF cycle.
 ! NB: In principle, one can rotate Vnl(b,b') using U^-1 from the Cholesky decomposition
 ! but the full Vnl matrix should be computed before the ortho step.

 ! Select value of after_ortho:
 !
 !   0: return with inconsistent eigenvalues, residuals and enlx_bk to avoid final H |Psi>.
 !   1: recompute enlx_bx after ortho. Return inconsistent eigens and residuals (last DIIS iteration).
 !   2: fully consistent mode: execute final H|Psi> after ortho step to update enlx_bx, eigens, residuals
 !
 ! Total number of H |Psi> applications:
 !
 !   1 for subdiago.
 !   1 for preconditioned steepest descent.
 !   (nline - 1) for DIIS or nline if ultimate accuracy is reached.
 !   1 if after_ortho > 0
 !
 after_ortho = 0
 if (accuracy_level >= 2) after_ortho = 1
 if (accuracy_level >= 4) after_ortho = 2
 if (after_ortho >= 1 .and. savemem == 0) after_ortho = 1
 ! It seems that PAW is more sensitive to after_ortho. Perhaps I can avoid the final H|phi> if accuracy_level == 1
 !if (usepaw == 1) after_ortho = 1
 !if (usepaw == 1) after_ortho = 0
 if (usepaw == 1) after_ortho = 2 ! FIXME ??

 if (after_ortho == 0) then
   !if (prtvol == -level)
   call wrtout(std_out, " VERY-FAST: Won't recompute data after orthogonalization.")

 !else if (after_ortho == 1 .and. savemem == 0) then
 else if (after_ortho == 1 .and. savemem == 0 .and. usepaw == 0) then
   !if (prtvol == -level)
   call wrtout(std_out, " FAST: Recomputing data by rotating matrix elements.")

   if (usepaw == 0 .or. has_fock) then
     ! Rotate gvnlxc by solving X_new U = Y_old for X with U upper triangle.
     ! Compute enlx with rotated cg and gvnlxc.
     if (istwf_k == 1) then
       call ZTRSM('R', 'U', 'N', 'N', npwsp, nband, cone, umat, nband, gvnlxc, npwsp)
     else
       call DTRSM('R', 'U', 'N', 'N', 2*npwsp, nband, one, umat, nband, gvnlxc, 2*npwsp)
     end if
     ABI_MALLOC(dots, (2, nband))
     call cg_zdotg_zip(istwf_k, npwsp, nband, option1, cg, gvnlxc, dots, me_g0, comm_bsf)
     enlx = dots(1,:)
     ABI_FREE(dots)
   end if

   if (.False.) then
   !if (usepaw == 1) then
     ! Compute new eigenvalues, residual vectors and norms.
     ! Rotate ghc by solving X_new U = Y_old for X with U upper triangle.
     if (istwf_k == 1) then
       call ZTRSM('R', 'U', 'N', 'N', npwsp, nband, cone, umat, nband, ghc, npwsp)
     else
       call DTRSM('R', 'U', 'N', 'N', 2*npwsp, nband, one, umat, nband, ghc, 2*npwsp)
     end if
     ABI_MALLOC_OR_DIE(gwork, (2, npwsp*nband), ierr)
     call cg_get_eigens(usepaw, istwf_k, npwsp, nband, cg, ghc, gsc, eig, me_g0, comm_bsf)
     call cg_get_residvecs(usepaw, npwsp, nband, eig, cg, ghc, gsc, gwork)
     call cg_norm2g(istwf_k, npwsp, nband, gwork, resid, me_g0, comm_bsf)
     call rmm_ydoc%add_tabular_line(resids2str("ortho_rot"), indent=0)
     ABI_FREE(gwork)
   end if

 else
   !if (prtvol == -level)
   if (after_ortho == 1) call wrtout(std_out, " SLOW: Recomputing enlx gvnlx by calling nonlop.")
   if (after_ortho == 2) call wrtout(std_out, " VERY-SLOW: Recomputing eigens and residues by calling getghc.")

   do iblock=1,nblocks
     igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
     ndat = (ige - igs + 1) / npwsp
     ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)
     cg_bk => cg(:,igs:ige); if (usepaw == 1) gsc_bk => gsc(1:2,igs:ige)
     if (savemem == 0) then
       ghc_bk => ghc(:,igs:ige); gvnlxc_bk => gvnlxc(:,igs:ige)
     end if

     select case (after_ortho)
     case (1)
       ! recompute NC enlx_bx after ortho.
       ! eigens and residuals are inconsistent as they have been computed before pw_orthon.
       signs = 1; paw_opt = 0
       if (usepaw == 1) then
         signs = 2; paw_opt = 3
       end if
       if (paral_kgb == 0) then
         call nonlop(choice1, cpopt, cprj_dum, enlx(ib_start:), gs_hamk, 0, eig(ib_start), &
                     mpi_enreg, ndat, 1, paw_opt, signs, gsc_bk, tim_nonlop, cg_bk, gvnlxc_bk)
       else
         call prep_nonlop(choice1, cpopt, cprj_dum, enlx(ib_start), gs_hamk, 0, eig(ib_start), &
                          ndat, mpi_enreg, 1, paw_opt, signs, gsc_bk, tim_nonlop, &
                          cg_bk, gvnlxc_bk, already_transposed=.False.)
       end if

     case (2)
       ! Consistent mode: update enlx_bx, eigens, residuals after orthogonalizalization.
       call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, cg_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                            eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, &
                            normalize=usepaw == 1)
     case default
       MSG_BUG(sjoin("Wrong after_ortho:", itoa(after_ortho)))
     end select
   end do ! iblock

   call rmm_ydoc%add_tabular_line(resids2str("after_ortho"), indent=0)
 end if ! after_ortho > 0

 !if (usepaw == 1) then
 !  !call cgtk_fixphase(cg, gsc, 0, 0, istwf_k, mcg, mgsc, mpi_enreg, nband, npwsp, usepaw)
 !  call cg_set_imag0_to_zero(istwf_k, me_g0, npwsp, nband, cg, max_absimag)
 !  call cg_set_imag0_to_zero(istwf_k, me_g0, npwsp, nband, gsc, max_absimag)
 !  call cgpaw_normalize(npwsp, nband, cg, gsc, istwf_k, me_g0, comm_bsf)
 !end if

 if (timeit) call cwtime_report(" after_ortho ", cpu, wall, gflops)

 ! Revert mixprec to previous status before returning.
 if (use_fft_mixprec) prev_mixprec = fftcore_set_mixprec(prev_mixprec)

 !if (dtset%prtvol > 0) then
 if (diis%stats%length() > 0) call rmm_ydoc%add_dict("skip_stats", diis%stats)
 call rmm_ydoc%write_and_free(std_out)

 if (timeit) call cwtime_report(" rmm_diis total: ", cpu_all, wall_all, gflops_all)

 ! Final cleanup.
 ABI_FREE(lambda_bk)
 ABI_FREE(dots_bk)
 ABI_FREE(residv_bk)
 ABI_FREE(kres_bk)
 ABI_FREE(umat)
 if (savemem == 0) then
   ABI_FREE(ghc)
   ABI_FREE(gvnlxc)
 else
   ABI_FREE(ghc_bk)
   ABI_FREE(gvnlxc_bk)
 end if
 call diis%free()

contains

function resids2str(level) result(str)
  character(len=*),intent(in) :: level
  character(len=500) :: str
  res_stats = stats_eval(resid(1:nb_pocc))
  !res_stats = stats_eval(resid(1:nband))
  write(str, "(1x, a12, 4(es10.3))") trim(level), res_stats%mean, res_stats%min, res_stats%max, res_stats%stdev
end function resids2str

end subroutine rmm_diis
!!***

!!****f* m_rmm_diis/rmm_diis_push_iter
!! NAME
!!  rmm_diis_push_iter
!!
!! FUNCTION
!!  Save one iteration of the DIIS algorithm.
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

subroutine rmm_diis_push_iter(diis, iter, ndat, eig_bk, resid_bk, enlx_bk, cg_bk, residv_bk, gsc_bk, tag)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat
 real(dp),intent(in) :: eig_bk(ndat), resid_bk(ndat), enlx_bk(ndat)
 real(dp),intent(in) :: cg_bk(2, diis%npwsp*ndat), residv_bk(2, diis%npwsp*ndat), gsc_bk(2, diis%npwsp*ndat*diis%usepaw)
 character(len=*),intent(in) :: tag

!Local variables-------------------------------
 integer :: idat, ibk
! *************************************************************************

 diis%last_iter = iter
 diis%hist_ene(iter, 1:ndat) = eig_bk
 diis%hist_resid(iter, 1:ndat) = resid_bk
 diis%hist_enlx(iter, 1:ndat) = enlx_bk
 diis%step_type(iter, 1:ndat) = tag

 do idat=1,ndat
   if (iter == 0) then
     if (diis%cplex == 2) then
       diis%resmat(:, 0, 0, idat) = [resid_bk(idat), zero]
     else
       diis%resmat(:, 0, 0, idat) = resid_bk(idat)
     end if
   end if
   !write(std_out, *)"res0", diis%resmat(:, 0, 0, idat)
   diis%step_type(iter, idat) = tag
   ibk = 1 + (idat - 1) * diis%npwsp
   call cg_zcopy(diis%npwsp, cg_bk(:,ibk), diis%chain_phi(:,:,iter,idat))
   call cg_zcopy(diis%npwsp, residv_bk(:,ibk), diis%chain_resv(:,:,iter,idat))
   if (diis%usepaw == 1) call cg_zcopy(diis%npwsp, gsc_bk(:,ibk), diis%chain_sphi(:,:,iter,idat))
 end do

end subroutine rmm_diis_push_iter
!!***

!!****f* m_rmm_diis/rmm_diis_exit_iter
!! NAME
!!  rmm_diis_exit_iter
!!
!! FUNCTION
!!  Return true if we can exit the DIIS iteration
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

logical function rmm_diis_exit_iter(diis, iter, ndat, niter_block, occ_bk, accuracy_ene, dtset, comm) result(ans)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat, niter_block, comm
 real(dp),intent(in) :: occ_bk(ndat)
 real(dp),intent(in) :: accuracy_ene
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: idat, ierr, nok, checks(ndat) !nbocc,
 real(dp) :: resid, deltae, deold , fact
 character(len=50) :: msg_list(ndat)
! *************************************************************************

 diis%last_iter = iter !; ans = .False.; return
 if (xmpi_comm_rank(comm) /= master) goto 10

 ! Tolerances depend on accuracy_level and occupation of the state.
 checks = 0

 do idat=1,ndat
   resid = diis%hist_resid(iter, idat)
   deold = diis%hist_ene(1, idat) - diis%hist_ene(0, idat)
   deltae = diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat)

   ! Relative criterion on eigenvalue differerence.
   ! Abinit default in the CG part is 0.005 that is really low (0.3 in V).
   ! Here we increase it depending whether the state is occupied or empty
   fact = one !; if (dtset%iscf > 0 .and. abs(occ_bk(idat)) < diis%tol_occupied) fact = three
   if (diis%accuracy_level == 1) fact = fact * 18
   if (diis%accuracy_level == 2) fact = fact * 12
   if (diis%accuracy_level == 3) fact = fact * 6
   if (abs(deltae) < fact * dtset%tolrde * abs(deold)) then
     checks(idat) = 1; msg_list(idat) = "deltae < fact * tolrde * deold"; cycle
   end if

   if (dtset%iscf < 0) then
     ! This is the only condition available for NSCF run.
     if (resid < dtset%tolwfr) then
       checks(idat) = 1; msg_list(idat) = 'resid < tolwfr'; cycle
     end if

   else
     ! Conditions available for SCF run.
     if (resid < dtset%tolwfr) then
       checks(idat) = 1; msg_list(idat) = 'resid < tolwfr'; cycle
     end if

     ! Absolute criterion on eigenvalue difference. Assuming error on Etot ~ band_energy.
     fact = one; if (dtset%iscf > 0 .and. abs(occ_bk(idat)) < diis%tol_occupied) fact = ten
     if (sqrt(abs(resid)) < fact * accuracy_ene) then
       checks(idat) = 1; msg_list(idat) = 'resid < accuracy_ene'; cycle
     end if
   end if
 end do ! idat

 ! Depending on the accuracy_level either full block or a fraction of it must pass the test in order to exit.
 nok = count(checks /= 0)
 if (diis%accuracy_level == 1) ans = nok >= 0.65_dp * ndat
 if (diis%accuracy_level == 2) ans = nok >= 0.75_dp * ndat
 if (diis%accuracy_level == 3) ans = nok >= 0.90_dp * ndat
 if (diis%accuracy_level == 4) ans = nok == ndat

 if (ans) then
   ! Log exit only if this is not the last iteration.
   if (iter /= niter_block) then
     do idat=1,ndat
       if (checks(idat) /= 0) call diis%stats%increment(msg_list(idat), 1)
     end do
   end if
 end if

 ! Broadcast final decision to all ranks.
 10 call xmpi_bcast(ans, master, comm, ierr)

end function rmm_diis_exit_iter
!!***

!!****f* m_rmm_diis/rmm_diis_print_block
!! NAME
!!  rmm_diis_print_block
!!
!! FUNCTION
!!  Print energies, residuals for a block of states.
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

subroutine rmm_diis_print_block(diis, ib_start, ndat, istep, ikpt, isppol)

 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: ib_start, ndat, istep, ikpt, isppol

!Local variables-------------------------------
 integer :: iter, idat, iband
 real(dp) :: deltae, deold, dedold, absdiff
 character(len=500) :: msg
! *************************************************************************

 call wrtout(std_out, &
   sjoin("<BEGIN RMM-DIIS-BLOCK, istep:", itoa(istep), ", ikpt:", itoa(ikpt), ", spin: ", itoa(isppol), ">"), &
   pre_newlines=1)

 do idat=1,ndat
   write(msg,'(1a, 2(a5), 4(a14), 1x, a6)') &
     "#", 'iter', "band", "eigen_eV", "eigde_meV", "de/dold", "resid", "type"; call wrtout(std_out, msg)

   iband = ib_start + idat - 1
   deold = diis%hist_ene(1, idat) - diis%hist_ene(0, idat)
   do iter=0,diis%last_iter
     dedold = zero; absdiff = zero
     if (iter > 0) then
       deltae = diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat)
       dedold = deltae / deold
       absdiff = (diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat))
     end if

     write(msg,"(1x, 2(i5), 4(es14.6), 1x, a6)") &
       iter, iband, diis%hist_ene(iter, idat) * Ha_eV, absdiff * Ha_meV, dedold, &
       diis%hist_resid(iter, idat), diis%step_type(iter, idat); call wrtout(std_out, msg)
   end do
 end do

 call wrtout(std_out, "<END RMM-DIIS-BLOCK>", newlines=1)

end subroutine rmm_diis_print_block
!!***

!!****f* m_rmm_diis/getghc_eigresid
!! NAME
!!  getghc_eigresid
!!
!! FUNCTION
!!  Compute new eigenvalues, residuals, H |psi> and enlx from cg and gsc.
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

subroutine getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, cg, ghc, gsc, mpi_enreg, prtvol, &
                           eig, resid, enlx, residvecs, gvnlxc, normalize)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: npw, my_nspinor, ndat, prtvol
 real(dp),intent(inout) :: cg(2, npw*my_nspinor*ndat)
 real(dp),intent(out) :: ghc(2,npw*my_nspinor*ndat), gsc(2,npw*my_nspinor*ndat*gs_hamk%usepaw)
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),intent(out) :: eig(ndat), resid(ndat), enlx(ndat)
 real(dp),intent(out) :: residvecs(2, npw*my_nspinor*ndat), gvnlxc(2, npw*my_nspinor*ndat)
 logical,optional,intent(in) :: normalize

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: istwf_k, usepaw, cpopt, sij_opt, npwsp, me_g0, comm_bsf
 real(dp),parameter :: rdummy = zero
 !real(dp) :: cpu, wall, gflops
 logical :: normalize_, has_fock
!arrays
 real(dp) :: dots(2, ndat)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 !if (timeit) call cwtime(cpu, wall, gflops, "start")
 normalize_ = .True.; if (present(normalize)) normalize_ = normalize
 npwsp = npw * my_nspinor
 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k; me_g0 = mpi_enreg%me_g0
 comm_bsf = mpi_enreg%comm_spinorfft; if (mpi_enreg%paral_kgb == 1) comm_bsf = mpi_enreg%comm_bandspinorfft
 has_fock = associated(gs_hamk%fockcommon)

 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! NC normalization.
 if (usepaw == 0 .and. normalize_) call cgnc_normalize(npwsp, ndat, cg, istwf_k, me_g0, comm_bsf)

 ! Compute H |cg>
 !call fock_set_ieigen(gs_hamk%fockcommon, iband)
 if (mpi_enreg%paral_kgb == 0) then
   call getghc(cpopt, cg, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
               rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
 else
   call prep_getghc(cg, gs_hamk, gvnlxc, ghc, gsc, rdummy, ndat, &
                    mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
 end if

 ! PAW normalization must be done here once gsc is known.
 if (usepaw == 1 .and. normalize_) call cgpaw_normalize(npwsp, ndat, cg, gsc, istwf_k, me_g0, comm_bsf)

 ! Compute new eigenvalues, residual vectors and norms.
 call cg_get_eigens(usepaw, istwf_k, npwsp, ndat, cg, ghc, gsc, eig, me_g0, comm_bsf)
 call cg_get_residvecs(usepaw, npwsp, ndat, eig, cg, ghc, gsc, residvecs)
 call cg_norm2g(istwf_k, npwsp, ndat, residvecs, resid, me_g0, comm_bsf)

 if (usepaw == 0 .or. has_fock) then
   ! Evaluate new enlx from gvnlxc.
   call cg_zdotg_zip(istwf_k, npwsp, ndat, option1, cg, gvnlxc, dots, me_g0, comm_bsf)
   enlx = dots(1,:)
 end if

 !if (timeit) call cwtime_report(" getghc_eigresid", cpu, wall, gflops)

end subroutine getghc_eigresid
!!***

!!****f* m_rmm_diis/rmm_diis_new
!! NAME
!!  rmm_diis_new
!!
!! FUNCTION
!!  Build new rmm_diis_t instance.
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

type(rmm_diis_t) function rmm_diis_new(accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol) result(diis)

!Arguments ------------------------------------
 integer,intent(in) :: accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol

! *************************************************************************

 diis%accuracy_level = accuracy_level
 diis%usepaw = usepaw
 diis%istwf_k = istwf_k
 diis%cplex = 2; if (istwf_k == 2) diis%cplex = 1
 diis%npwsp = npwsp
 diis%max_niter = max_niter
 diis%bsize = bsize
 diis%prtvol = prtvol

 ABI_MALLOC(diis%hist_ene, (0:max_niter+2, bsize))
 ABI_MALLOC(diis%hist_resid, (0:max_niter+2, bsize))
 ABI_CALLOC(diis%hist_enlx, (0:max_niter+2, bsize))
 ABI_MALLOC(diis%step_type, (0:max_niter+2, bsize))
 ABI_MALLOC(diis%chain_phi, (2, npwsp, 0:max_niter, bsize))
 ABI_MALLOC(diis%chain_sphi, (2, npwsp*usepaw, 0:max_niter, bsize))
 ABI_MALLOC(diis%chain_resv, (2, npwsp, 0:max_niter, bsize))
 ABI_CALLOC(diis%resmat, (diis%cplex, 0:max_niter, 0:max_niter, bsize)) ! <R_i|R_j>

end function rmm_diis_new
!!***

!!****f* m_rmm_diis/rmm_diis_free
!! NAME
!!  rmm_diis_free
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

subroutine rmm_diis_free(diis)

!Arguments ------------------------------------
 class(rmm_diis_t),intent(inout) :: diis

! *************************************************************************

 ABI_SFREE(diis%hist_ene)
 ABI_SFREE(diis%hist_resid)
 ABI_SFREE(diis%hist_enlx)
 ABI_SFREE(diis%step_type)
 ABI_SFREE(diis%chain_phi)
 ABI_SFREE(diis%chain_sphi)
 ABI_SFREE(diis%chain_resv)
 ABI_SFREE(diis%resmat)

 call diis%stats%free()

end subroutine rmm_diis_free
!!***

!!****f* m_rmm_diis/rmm_diis_update_block
!! NAME
!!  rmm_diis_update_block
!!
!! FUNCTION
!!  Compute new trial wavefunctions and residuals from the DIIS chain.
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

subroutine rmm_diis_update_block(diis, iter, npwsp, ndat, cg_bk, residv_bk, comm)

!Arguments ------------------------------------
 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iter, npwsp, comm, ndat
 real(dp),intent(inout) :: cg_bk(2, npwsp, ndat), residv_bk(2, npwsp, ndat)

!local variables
 integer,parameter :: master = 0
 integer :: cplex, ierr, nprocs, my_rank, idat
 real(dp) :: noise !, cpu, wall, gflops
 real(dp),allocatable :: wmat1(:,:,:), wvec(:,:,:), alphas(:,:)
 character(len=500) :: msg
 ! *************************************************************************

 !if (timeit) call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 cplex = diis%cplex

 ! Solve system of linear equations.
 ! Only master works so that we are sure we have the same solution.
 ABI_CALLOC(wvec, (cplex, 0:iter, ndat))

 if (my_rank == master) then
   ABI_MALLOC(wmat1, (cplex, 0:iter, 0:iter))

   do idat=1,ndat
     !if (mod(idat, nprocs) /= my_rank) cycle ! MPI parallelism
     wvec(1, iter, idat) = -one
     wmat1 = zero
     wmat1(1,:,iter) = -one
     wmat1(1,iter,:) = -one
     wmat1(1,iter,iter) = zero
     wmat1(:,0:iter-1, 0:iter-1) = diis%resmat(:, 0:iter-1, 0:iter-1, idat)

     call xhesv_cplex("U", cplex, iter+1, 1, wmat1, wvec(:,:,idat), msg, ierr)
     ABI_CHECK(ierr == 0, msg)

     !if (diis%prtvol == -level) then
     !  write(std_out,*)"wvec:", wvec(:,:,idat)
     !  write(std_out,*)"sum(wvec):", sum(wvec(:, 0:iter-1, idat), dim=2)
     !end if
     if (cplex == 2) then
       ! coefficients should sum up to 1 but sometimes we get a small imaginary part. here we remove it
       noise = sum(wvec(2, 0:iter-1, idat))
       wvec(2, 0:iter-1, idat) = wvec(2, 0:iter-1, idat) - noise * iter
     end if

   end do
   ABI_FREE(wmat1)
 end if

 ! Master broadcasts data.
 if (nprocs > 1) call xmpi_bcast(wvec, master, comm, ierr)
 !if (nprocs > 1) call xmpi_sum(wvec, comm, ierr)

 ! Take linear combination of chain_phi and chain_resv.
 do idat=1,ndat
   if (cplex == 2) then
     call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), wvec(:,:,idat), cg_bk(:,:,idat))
     call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), wvec(:,:,idat), residv_bk(:,:,idat))
   else
     ! coefficients are real --> use DGEMV
     ABI_MALLOC(alphas, (1, 0:iter))
     alphas(1,:) = wvec(1,:,idat)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_phi(:,:,:,idat), 2*npwsp, alphas, 1, zero, cg_bk(:,:,idat), 1)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_resv(:,:,:,idat), 2*npwsp, alphas, 1, zero, residv_bk(:,:,idat), 1)
     ABI_FREE(alphas)
   end if
 end do

 ABI_FREE(wvec)
 !if (timeit) call cwtime_report(" update_block", cpu, wall, gflops)

end subroutine rmm_diis_update_block
!!***

!!****f* m_rmm_diis/rmm_diis_eval_mats
!! NAME
!!  rmm_diis_eval_mats
!!
!! FUNCTION
!!  Compute matrix elements required by the RMM-DIIS method.
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

subroutine rmm_diis_eval_mats(diis, iter, ndat, me_g0, comm)

!Arguments ------------------------------------
 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat, me_g0, comm

!local variables
 integer :: ii, ierr, idat, nprocs, option
 real(dp) :: dotr, doti !, cpu, wall, gflops
 !integer :: requests(ndat)
! *************************************************************************

 !if (timeit) call cwtime(cpu, wall, gflops, "start")
 nprocs = xmpi_comm_size(comm)
 option = 2; if (diis%cplex == 1) option = 1

 do idat=1,ndat

   do ii=0,iter
     ! <R_i|R_j>
     call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option, &
                    diis%chain_resv(:,:,ii,idat), diis%chain_resv(:,:,iter,idat), me_g0, xmpi_comm_self)
     if (ii == iter) doti = zero
     if (diis%cplex == 2) then
       diis%resmat(:, ii, iter, idat) = [dotr, doti]
     else
       diis%resmat(1, ii, iter, idat) = dotr
     end if
   end do ! ii

   !if (nprocs > 1) then
   !  call xmpi_sum(diis%resmat(:,0:iter,iter,idat), comm, ierr)
   !  !call xmpi_isum_ip(diis%resmat(:,0:iter,iter,idat), comm, requests(idat), ierr)
   !endif
   !if (diis%prtvol == -level) write(std_out,*)"iter, idat, resmat:", iter, idat, diis%resmat(:,0:iter,iter,idat)
 end do ! idat

 if (nprocs > 1) call xmpi_sum(diis%resmat(:,0:iter,iter,1:ndat), comm, ierr)
 !if (nprocs > 1) call xmpi_waitall(requests, ierr)
 !if (timeit) call cwtime_report(" eval_mats", cpu, wall, gflops)

end subroutine rmm_diis_eval_mats
!!***

!!****f* m_numeric_tools/my_pack_matrix
!! NAME
!! my_pack_matrix
!!
!! FUNCTION
!! Packs a matrix into hermitian format
!!
!! INPUTS
!! N: size of matrix
!! cplx: 2 if matrix is complex, 1 for real matrix.
!! mat_in(cplx, N*N)= matrix to be packed
!!
!! OUTPUT
!! mat_out(cplx*N*N+1/2)= packed matrix (upper triangle)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine my_pack_matrix(n, mat_in, mat_out)

!Arguments ------------------------------------
 integer, intent(in) :: N
 real(dp), intent(in) :: mat_in(N, N)
 real(dp), intent(out) :: mat_out(2, N*(N+1)/2)

!local variables
 integer :: isubh, i, j
! *************************************************************************

 isubh = 1
 do j=1,N
   do i=1,j
     mat_out(1,isubh) = mat_in(i, j)
     mat_out(2,isubh) = zero
     isubh = isubh + 1
   end do
 end do

end subroutine my_pack_matrix
!!***

!!****f* ABINIT/subspace_rotation
!! NAME
!! subspace_rotation
!!
!! FUNCTION
!!  This routine computes the <i|H|j> matrix elements and then performs the subspace rotation
!!  of the orbitals (rayleigh-ritz procedure)
!!  The main difference with respect to other similar routines is that this implementation does not require
!!  the <i|H|j> matrix elements as input so it can be used before starting the wavefunction optimation
!!  as required e.g. by the RMM-DIIS method.
!!  Moreover, the routine computes the new rediduals after the subspace rotation by rotating the
!!  matrix elements of the Hamiltonian in the new basis (requires more memory but client code
!!  can avoid calling getghc after subspace_rotation.
!!
!! INPUTS
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  ptrvol
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point and spin
!!  npw=number of plane waves at this k point
!!  my_nspinor=number of spinor components treated by this MPI proc.
!!
!! OUTPUT
!!  eig(nband): New eigvalues from subspace rotation.
!!  If usepaw==1:
!!    gsc(2,*)=<g|S|c> matrix elements (S=overlap)
!!  enlx(nband)=contribution from each band to nonlocal psp + potential Fock ACE part of total energy, at this k-point
!!  ghc
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  gsc(2,*)=update <G|S|C>
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine subspace_rotation(gs_hamk, prtvol, mpi_enreg, nband, npw, my_nspinor, savemem, enlx, eig, cg, gsc, ghc, gvnlxc)

!Arguments ------------------------------------
 integer,intent(in) :: prtvol, nband, npw, my_nspinor, savemem
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*my_nspinor*nband)
 real(dp),target,intent(inout) :: gsc(2,npw*my_nspinor*nband*gs_hamk%usepaw)
 !real(dp),target,intent(out) :: ghc(2,npw*my_nspinor*nband)
 !real(dp),target,intent(out) :: gvnlxc(2,npw*my_nspinor*nband)
 real(dp),target,allocatable,intent(out) :: ghc(:,:), gvnlxc(:,:)
 real(dp),intent(out) :: eig(nband), enlx(nband)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, tim_getghc = 0, use_subovl0 = 0, option1 = 1
 integer :: ig, ig0, ib, ierr, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, paral_kgb
 integer :: iband, cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, usepaw, me_g0, cplex, comm_bsf
 logical :: has_fock
 real(dp),parameter :: rdummy = zero
 real(dp) :: cpu, wall, gflops
!arrays
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp) :: subovl(use_subovl0)
 real(dp),allocatable :: subham(:), h_ij(:,:,:), evec(:,:,:), evec_re(:,:), gwork(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: ghc_bk(:,:), gvnlxc_bk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_bk(:,:)
 real(dp) :: dots(2, nband)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k
 paral_kgb = mpi_enreg%paral_kgb; me_g0 = mpi_enreg%me_g0
 comm_bsf = mpi_enreg%comm_spinorfft; if (mpi_enreg%paral_kgb == 1) comm_bsf = mpi_enreg%comm_bandspinorfft
 npwsp = npw * my_nspinor
 has_fock = associated(gs_hamk%fockcommon)

 ! =======================================
 ! Apply H to input cg to compute <i|H|j>
 ! =======================================
 gsc_bk => fake_gsc_bk
 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! Treat states in groups of bsize bands even when paral_kgb = 0
 bsize = 8; if (paral_kgb == 1) bsize = mpi_enreg%nproc_band * mpi_enreg%bandpp
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 cplex = 2; if (istwf_k /= 1) cplex = 1
 !cplex = 2; if (istwf_k == 2) cplex = 1

 ABI_CALLOC(h_ij, (cplex, nband, nband))
 ! Allocate full ghc and gvnlxc to be able to rotate residuals and Vnlx matrix elements
 ! after subdiago. More memory but we can save a call to H|psi>.
 if (savemem == 0) then
   ABI_MALLOC_OR_DIE(ghc, (2, npwsp*nband), ierr)
   ABI_MALLOC_OR_DIE(gvnlxc, (2, npwsp*nband), ierr)
 else if (savemem == 1) then
   ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
   ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))
 else
   MSG_ERROR(sjoin("Invalid savemem:", itoa(savemem)))
 end if

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)
   if (usepaw == 1) gsc_bk => gsc(1:2,igs:ige)
   if (savemem == 0) then
     ghc_bk => ghc(:, igs:ige); gvnlxc_bk => gvnlxc(:, igs:ige)
   end if

   if (paral_kgb == 0) then
     call getghc(cpopt, cg(:,igs:ige), cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                 rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
   else
     call prep_getghc(cg(:,igs:ige), gs_hamk, gvnlxc_bk, ghc_bk, gsc_bk, rdummy, ndat, &
                      mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
   end if

   ! Compute <i|H|j> for i=1,nband and all j in block
   if (cplex == 2) then
     call cg_zgemm("C", "N", npwsp, nband, ndat, cg, ghc_bk, h_ij(:,:,ib_start))
   else
     call dgemm("T", "N", nband, ndat, 2*npwsp, one, cg, 2*npwsp, ghc_bk, 2*npwsp, zero, h_ij(:,:,ib_start), nband)
   end if

   if (istwf_k /= 1) then
     do iband=ib_start, ib_stop
       h_ij(:,:,iband) = two * h_ij(:,:,iband)

       if (istwf_k == 2 .and. me_g0 == 1) then
         ! Gamma k-point and I have G=0. Remove double counting term.
         ig = 1 + (iband - ib_start) * npwsp
         do ib=1,nband
           ig0 = 1 + npwsp * (ib - 1)
           h_ij(1,ib,iband) = h_ij(1,ib,iband) - cg(1,ig0) * ghc_bk(1,ig)
         end do
       end if

       ! Force real matrix.
       if (cplex == 2) h_ij(2,:,iband) = zero
     end do
   end if

 end do ! iblock

 ! Pack <i|H|j> to prepare call to subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 if (cplex == 2) then
   do iband=1,nband
     h_ij(2,iband,iband) = zero ! Force diagonal elements to be real
   end do
 end if

 if (cplex == 2) then
   call pack_matrix(h_ij, subham, nband, 2)
 else
   call my_pack_matrix(nband, h_ij, subham)
 end if

 ABI_FREE(h_ij)
 call xmpi_sum(subham, comm_bsf, ierr)
 if (timeit) call cwtime_report(" subspace build Hij", cpu, wall, gflops)

 ! ========================
 ! Subspace diagonalization
 ! =======================
 ! Rotate cg, gsc and compute new eigevalues.
 ABI_MALLOC(evec, (2, nband, nband))
 mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 call subdiago(cg, eig, evec, gsc, 0, 0, istwf_k, mcg, mgsc, nband, npw, my_nspinor, paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)

 ABI_FREE(subham)
 if (timeit) call cwtime_report(" subspace subdiago", cpu, wall, gflops)

 if (savemem == 0) then
   ! Rotate ghc matrix in the new subspace:
   !
   !      new_{g,b} = old_{g,i} evec_{i,b}
   !
   ! cg and PAW gsc have been already rotated in subdiago
   !
   if (cplex == 1) then
     ! Eigenvectors are real.
     ABI_MALLOC(evec_re, (nband, nband))
     evec_re = evec(1,:,:)
   end if

   ABI_MALLOC_OR_DIE(gwork, (2, npwsp*nband), ierr)
   if (cplex == 1) then
     call DGEMM("N", "N", 2*npwsp, nband, nband, one, ghc, 2*npwsp, evec_re, nband, zero, gwork, 2*npwsp)
   else
     call abi_zgemm_2r("N", "N", npwsp, nband, nband, cone, ghc, npwsp, evec, nband, czero, gwork, npwsp)
   end if
   call cg_zcopy(npwsp * nband, gwork, ghc)

   ! Rotate <G|Vnlx|Psi_n> and evaluate new enlx for NC.
   if (usepaw == 0 .or. has_fock) then
     if (cplex == 1) then
       call DGEMM("N", "N", 2*npwsp, nband, nband, one, gvnlxc, 2*npwsp, evec_re, nband, zero, gwork, 2*npwsp)
     else
       call abi_zgemm_2r("N", "N", npwsp, nband, nband, cone, gvnlxc, npwsp, evec, nband, czero, gwork, npwsp)
     end if
     !call abi_xgemm('N','N', vectsize, nband, nband, cone, gvnlxc, vectsize, evec, nband, czero, gwork, vectsize, x_cplx=cplx)
     call cg_zcopy(npwsp * nband, gwork, gvnlxc)
     call cg_zdotg_zip(istwf_k, npwsp, nband, option1, cg, gvnlxc, dots, me_g0, comm_bsf)
     enlx = dots(1,:)
   end if
   ABI_FREE(gwork)
   ABI_SFREE(evec_re)
 end if

 ABI_FREE(evec)
 if (savemem == 1) then
   ABI_FREE(ghc_bk)
   ABI_FREE(gvnlxc_bk)
 end if

 if (timeit) call cwtime_report(" subspace final rotation", cpu, wall, gflops)

end subroutine subspace_rotation
!!***

end module m_rmm_diis
!!***
