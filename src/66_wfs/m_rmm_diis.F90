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
 use m_numeric_tools, only : pack_matrix, imin_loc
 use m_hide_lapack,   only : xhegv_cplex, xhesv_cplex
 use m_pair_list,     only : pair_list
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_fftcore,       only : fftcore_set_mixprec
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 use m_nonlop,        only : nonlop
 !use m_fock,         only : fock_set_ieigen, fock_set_getghc_call

 implicit none

 private
!!***

 public :: rmm_diis
!!***

 type,private :: rmm_diis_t

   integer :: accuracy_level
   ! Defines tolerances, activates/deactivates tricks

   integer :: usepaw
   ! 1 if we are running PAW.

   integer :: istwf_k
   ! wavefunction storage mode.

   integer :: cplex
   ! 1 if matrices are real (Gamma-point), 2 for complex

   integer :: bsize
   ! (Max) block size for bands

   integer :: max_niter
   ! Maximum number of iterations

   integer :: npwsp
   ! npw * my_nspinor

   integer :: prtvol
   ! vervosity level

   integer :: last_iter
   ! Last RMM-DIIS iteration performed.

   integer :: use_smat = 0

   real(dp) :: tol_occupied
   ! Tolerance for partial occupied states

   type(pair_list) :: stats

   real(dp),allocatable :: hist_ene(:,:)
   real(dp),allocatable :: hist_resid(:,:)
   real(dp),allocatable :: hist_enlx(:,:)
   character(len=7),allocatable :: step_type(:,:)
   ! (0:max_niter+2, bsize)
   ! 0 is the initial step, then DIIS iterations whose number may depend on the block
   ! followed by an optional trial step and the computation of eigens after ortho.

   real(dp),allocatable :: resmat(:,:,:,:)
   real(dp),allocatable :: smat(:,:,:,:)
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
   procedure :: rollback => rmm_diis_rollback           ! select the trial states in the chain with smaller resid.
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
!!  istep,ikpt,isppol
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  my_nspinor=number of plane waves at this k point
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!  resid(nband)=residuals for each states
!!  If usepaw==1:
!!    gsc_all(2,*)=<g|s|c> matrix elements (s=overlap)
!!  If usepaw==0
!!    enlx(nband)=contribution from each band to nonlocal psp + potential Fock ACE part
!!                of total energy, at this k-point
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  rmm_diis_status(2): Status of the eigensolver.
!!    The first entry gives the previous accuracy.
!!    The second entry gives the number of iterations already performed with this level.
!!
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine rmm_diis(istep, ikpt, isppol, cg, dtset, eig, occ, enlx, gs_hamk, kinpw, gsc_all, &
                    mpi_enreg, nband, npw, my_nspinor, resid, rmm_diis_status)

!Arguments ------------------------------------
 integer,intent(in) :: istep, ikpt, isppol, nband, npw, my_nspinor
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*my_nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*my_nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband), resid(nband)
 real(dp),intent(in) :: occ(nband), kinpw(npw)
 real(dp),intent(out) :: eig(nband)
 integer,intent(inout) :: rmm_diis_status(2)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer,parameter :: choice1 = 1, signs1 = 1, signs2 = 2, tim_nonlop = 0, paw_opt0 = 0, paw_opt3 = 3
 integer :: ierr, prtvol, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, idat, paral_kgb, ortalgo
 integer :: cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, optekin, usepaw, iter, max_niter, max_niter_block
 integer :: me_g0, nb_pocc, jj, kk, it, accuracy_level, raise_acc, prev_mixprec, after_ortho
 integer :: comm_bandspinorfft, prev_accuracy_level, ncalls_with_prev_accuracy !, nspinor
 logical :: end_with_trial_step, first_call, use_fft_mixprec
 real(dp),parameter :: rdummy = zero
 real(dp) :: accuracy_ene, cpu, wall, gflops, max_res_pocc, tol_occupied
 !character(len=500) :: msg
 !character(len=6) :: tag
 type(rmm_diis_t) :: diis
!arrays
 real(dp) :: tsec(2)
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp),allocatable :: evec(:,:,:)
 real(dp),allocatable :: ghc_bk(:,:), gvnlxc_bk(:,:), lambda_bk(:), residv_bk(:,:), kres_bk(:,:), dots_bk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_bk(:,:), phi_bk(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; paral_kgb = mpi_enreg%paral_kgb
 me_g0 = mpi_enreg%me_g0; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
 !nspinor = dtset%nspinor
 npwsp = npw * my_nspinor; mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 prtvol = dtset%prtvol !; prtvol = -level

 ! =================
 ! Prepare DIIS loop
 ! =================
 ! accuracy_level is computed from the maxval of the previous residuals that are received in input.
 !
 !  1: Used at the beginning of the SCF cycle. Use loosy convergence criteria in order
 !     to reduce the number of H|psi> applications as much as possible so that
 !     we can start to mix densities/potentials.
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

 ! NB: Accuracy_level is not allowed to increase during the SCF cycle
 !
 ! Since we operate on blocks of bands, all the states in the block will receive the same treatment.
 ! This means that one can observe different convergence behaviour depending on bsize.
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
   call wrtout(std_out, "Accuracy_level is automatically increased as we reached the max number of NSCF iterations.")
 end if
 raise_acc = max(raise_acc, prev_accuracy_level)

 ! Define tolerance for occupied states on the basis of prev_accuracy_level
 ! and compute max of residuals for this bands.
 tol_occupied = tol3; if (any(prev_accuracy_level == [1])) tol_occupied = tol2
 nb_pocc = count(occ > tol_occupied)
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
 !if (usepaw == 1) accuracy_level = 4

 ! Update rmm_diis_status. Reset number of calls if we've just moved to a new accuracy_level.
 rmm_diis_status(1) = accuracy_level
 if (accuracy_level /= prev_accuracy_level) rmm_diis_status(2) = 0
 rmm_diis_status(2) = rmm_diis_status(2) + 1

 ! Will perform max_niter DIIS steps. Usually 3 as nline by default is 4.
 ! Optionally end with trial step after DIIS.
 max_niter = max(dtset%nline - 1, 1); end_with_trial_step = .True.
 if (accuracy_level >= 4) then
   ! Prefer more DIIS iterations with rollback over the final trial step.
   max_niter = dtset%nline; end_with_trial_step = .False.
 end if
 ! FIXME: There's no real evidence that trial step is really needed
 end_with_trial_step = .False.

 ! Define accuracy_ene for SCF.
 accuracy_ene = zero
 if (dtset%iscf > 0) then
   if (dtset%toldfe /= zero) then
     accuracy_ene = dtset%toldfe * ten**(-accuracy_level + 2) !/ nb_pocc
   else
     ! We are not using toldfe to stop the SCF cycle
     ! so we are forced to hardcode a tolerance for the absolute diff in the KS eigenvalue.
     accuracy_ene = tol8 * ten**(-accuracy_level + 2) !/ nb_pocc
   end if
 end if

 ! Select value of after_ortho:
 !   0: return with inconsistent eigenvalues, residuals and enlx_bk to avoid final H |Psi>.
 !   1: recompute enlx_bx after ortho. Return inconsistent eigens and residuals (last DIIS iteration).
 !   2: fully consistent mode: execute final H|Psi> after ortho step to update enlx_bx, eigens, residuals
 !
 ! Total number of H |Psi> applications:
 !
 !   1 for subdiago.
 !   2 for preconditioned steepest descent.
 !   (nline - 1) for DIIS or nline if ultimate accuracy is reached.
 !   1 if after_ortho > 0
 !
 after_ortho = 0
 if (accuracy_level >= 2) after_ortho = 1
 if (accuracy_level >= 4) after_ortho = 2
 !after_ortho = 2
 if (usepaw == 1) after_ortho = 2 ! FIXME

 ! Use mixed precisions if requested by user but only at the beginning.
 use_fft_mixprec = dtset%mixprec == 1 .and. accuracy_level < 2
 if (use_fft_mixprec) prev_mixprec = fftcore_set_mixprec(1)

 ! Select preconditioning.
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1 ! optekin = 0

 ! Will treat states in groups of bsize bands.
 bsize = 8  ! default value for paral_kgb = 0
 if (paral_kgb == 1) bsize = mpi_enreg%nproc_band * mpi_enreg%bandpp
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 ! Build DIIS object.
 diis = rmm_diis_new(accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol)
 diis%tol_occupied = tol_occupied
 if (end_with_trial_step) then
   call wrtout(std_out, sjoin(" Using Max", itoa(max_niter), "RMM-DIIS iterations + optional final trial step."))
 else
   call wrtout(std_out, sjoin(" Using Max", itoa(max_niter), "RMM-DIIS iterations"))
 end if
 call wrtout(std_out, sjoin( &
   " Number of blocks:", itoa(nblocks), ", nb_pocc:", itoa(nb_pocc)))
 call wrtout(std_out, sjoin( &
   " Max_input_resid_pocc", ftoa(max_res_pocc), "accuracy_level:", itoa(accuracy_level), &
   ", accuracy_ene: ", ftoa(accuracy_ene)))
 call timab(1634, 1, tsec) ! "rmm_diis:band_opt"

 ! =========================
 ! === Subspace rotation ===
 ! =========================
 call subspace_rotation(gs_hamk, dtset, mpi_enreg, nband, npw, my_nspinor, eig, cg, gsc_all, evec)
 ABI_FREE(evec)

 ABI_MALLOC(lambda_bk, (bsize))
 ABI_MALLOC(dots_bk, (2, bsize))
 ABI_MALLOC(residv_bk, (2, npwsp*bsize))
 ABI_MALLOC(kres_bk, (2, npwsp*bsize))
 ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
 ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))

 gsc_bk => fake_gsc_bk
 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! We loop over nblocks, each block has ndat states.
 ! - Convergence behaviour may depend on bsize as branches are taken according to
 !   the status of all bands in the block.
 ! TODO: Transpose only once per block and then work with already_transposed = .True.

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

   ! Reduce number of niter iterations if block contains "empty" states.
   ! This should happen only if npband is small wrt nband and nband >> nbocc.
   ! TODO: Don't reduce niter if MD
   !if dtset%
   !if dtset%occopt == 2
   max_niter_block = max_niter
   if (all(occ(ib_start:ib_stop) < diis%tol_occupied)) max_niter_block = max(1 + max_niter / 2, 2)

   ! Compute H |phi_0> with cg block after subdiago.
   phi_bk => cg(:,igs:ige); if (usepaw == 1) gsc_bk => gsc_all(1:2,igs:ige)

   call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, phi_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                        eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, normalize=.False.)

   ! Save <R0|R0> and <phi_0|S|phi_0>, |phi_0>, |S phi_0>. Assume input phi_bk is already S-normalized.
   call diis%push_iter(0, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), phi_bk, residv_bk, gsc_bk, "SDIAG")
   if (timeit) call cwtime_report(" first getghc_eigresid ", cpu, wall, gflops)

   ! Line minimization with preconditioned steepest descent:
   !
   !    |phi_1> = |phi_0> + lambda |K R_0>
   !
   ! where lambda minimizes the residual (not the Rayleigh quotient as in Kresse's paper).
   !
   !    lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
   !
   ! more expensive than finding the stationary point of the Rayleigh quotient as it requires an extra H application
   ! but it should be more stable and more consistent with the RMM approach.
   !
   ! Precondition |R_0>, output in kres_bk = |K R_0>
   call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
   call cg_precon_many(istwf_k, npw, my_nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)
   if (timeit) call cwtime_report(" first cg_precon ", cpu, wall, gflops)

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
   call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_bandspinorfft)
   if (timeit) call cwtime_report(" cg_norm2g ", cpu, wall, gflops)

   ! Compute lambda
   dots_bk = zero
   do idat=1,ndat
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     call dotprod_g(dots_bk(1,idat), dots_bk(2,idat), istwf_k, npwsp, option1, &
                    diis%chain_resv(:,:,0,idat), residv_bk(:,jj), me_g0, xmpi_comm_self)
   end do
   call xmpi_sum(dots_bk, comm_bandspinorfft, ierr)

   ! Build |Psi_1> = |Phi_0> + lambda |K R_0>
   do idat=1,ndat
     lambda_bk(idat) = -dots_bk(1,idat) / lambda_bk(idat)
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     phi_bk(:,jj:kk) = diis%chain_phi(:,:,0,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
   end do
   if (timeit) call cwtime_report(" trial_step ", cpu, wall, gflops)

   ! ===============
   ! DIIS iterations
   ! ===============
   iter_loop: do iter=1,max_niter_block

     if (iter > 1) then
       ! Solve DIIS equations and update phi_bk and residv_bk for iter > 1
       call diis%update_block(iter, npwsp, ndat, lambda_bk, phi_bk, residv_bk, comm_bandspinorfft)

       ! Precondition residual, output in kres_bk.
       call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
       call cg_precon_many(istwf_k, npw, my_nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)

       ! Compute phi_bk with the lambda(ndat) obtained at iteration #1
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)
     end if

     ! Compute H |phi_now> and evaluate new enlx for NC.
     call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, phi_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, normalize=.True.)

     ! Store new residuals.
     call diis%push_iter(iter, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), phi_bk, residv_bk, gsc_bk, "DIIS")

     ! CHECK FOR CONVERGENCE
     if (diis%exit_iter(iter, ndat, max_niter_block, occ(ib_start), accuracy_ene, &
                        dtset, comm_bandspinorfft)) exit iter_loop

     ! Compute <R_i|R_j> and <i|S|j> for j=iter
     if (iter /= max_niter_block) call diis%eval_mats(iter, ndat, me_g0, comm_bandspinorfft)
   end do iter_loop
   if (timeit) call cwtime_report(" iterloop ", cpu, wall, gflops)

   call diis%rollback(npwsp, ndat, phi_bk, gsc_bk, residv_bk,  &
                      eig(ib_start), resid(ib_start), enlx(ib_start), comm_bandspinorfft)

   if (timeit) call cwtime_report(" rollback ", cpu, wall, gflops)

   if (end_with_trial_step) then
     ! The final trial step is cheap and effective provided the previous lambda is still accurate.
     ! Computing new residuals at this level is expensive, moreover we still need to orthogonalizalize
     ! the states before returning so here we just update the trial states without recomputing
     ! eigenvalues, residuals and enlx.
     call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
     call cg_precon_many(istwf_k, npw, my_nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)

     ! Build |Psi_1> = |Phi_0> + lambda |K R_0>
     ! Reuse previous lambda.
     !tag = "FIXLAM"
     !lambda_bk = limit_lambda(lambda_bk)
     call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)

     if (usepaw == 1) then
       ! Need to recompute <G|S|Psi> before calling pw_orthon
#if 0
       call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, phi_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                            eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, normalize=.True.)
#else
       !dimenl1=gs_ham%dimekb1;dimenl2=natom;tim_nonlop=0
       !choice=1;signs=2;cpopt=-1+3*gs_ham%usecprj;paw_opt=3;useylm=1

       if (paral_kgb == 0) then
         call nonlop(choice1, cpopt, cprj_dum, enlx(ib_start:), gs_hamk, 0, eig(ib_start), &
                     mpi_enreg, ndat, 1, paw_opt3, signs2, gsc_bk, tim_nonlop, phi_bk, gvnlxc_bk)
       else
         call prep_nonlop(choice1, cpopt, cprj_dum, enlx(ib_start), gs_hamk, 0, eig(ib_start), &
            ndat, mpi_enreg, 1, paw_opt3, signs2, gsc_bk, tim_nonlop, &
            phi_bk, gvnlxc_bk, already_transposed=.False.)
       end if
       !call cgpaw_normalize(npwsp, ndat, phi_bk, gsc_bk, istwf_k, me_g0, comm_bandspinorfft)
#endif
     else
        ! In principle this is not needed but as it will be done in pw_orthon
        !cgnc_normalize(npwsp, ndat, phi_bk, istwf_k, me_g0, comm_bandspinorfft)
     end if

     if (timeit) call cwtime_report(" last_trial_step", cpu, wall, gflops)
   end if ! end_with_trial_step

   if (prtvol == -level) call diis%print_block(ib_start, ndat, istep, ikpt, isppol)
 end do ! iblock

 call timab(1634, 2, tsec) !"rmm_diis:band_opt"
 if (timeit) call cwtime_report(" rmm_diis:band_opt", cpu, wall, gflops)

 ! ===============================
 ! Orthogonalize states after DIIS
 ! ===============================
 call timab(583,1,tsec) ! "vtowfk(pw_orthon)"
 ortalgo = mpi_enreg%paral_kgb !; ortalgo = 3
 call pw_orthon(0, 0, istwf_k, mcg, mgsc, npwsp, nband, ortalgo, gsc_all, usepaw, cg, me_g0, comm_bandspinorfft)
 !call fxphas(cg, gsc_all, 0, 0, istwf_k, mcg, mgsc, mpi_enreg, nband, npwsp, usepaw)
 call timab(583,2,tsec)
 if (timeit) call cwtime_report(" pw_orthon ", cpu, wall, gflops)

 if (after_ortho > 0) then
   ! Recompute eigenvalues, residuals, and enlx after orthogonalization.
   ! This step is important to improve the convergence of the NC total energy
   ! and it guarantees that eigenvalues and residuals are consistent with the output wavefunctions.
   ! but we try to avoid it at the beginning of the SCF cycle.
   ! NB: In principle, one can rotate Vnl(b,b') using the U^-1 from the Cholesky decomposition
   ! but the full Vnl matrix should be computed before the ortho step.

   !if (prtvol > 0)
   call wrtout(std_out, " Recomputing eigenvalues and residues after orthogonalization.")
   do iblock=1,nblocks
     igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
     ndat = (ige - igs + 1) / npwsp
     ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)
     phi_bk => cg(:,igs:ige); if (usepaw == 1) gsc_bk => gsc_all(1:2,igs:ige)

     select case (after_ortho)
     case (1)
       ! recompute NC enlx_bx after ortho. eigens and residuals are inconsistent
       ! as they have been computed before the ortho step
       if (usepaw == 0) then
         if (paral_kgb == 0) then
           call nonlop(choice1, cpopt, cprj_dum, enlx(ib_start:), gs_hamk, 0, eig(ib_start), &
                       mpi_enreg, ndat, 1, paw_opt0, signs1, gsc_bk, tim_nonlop, phi_bk, gvnlxc_bk)
         else
           call prep_nonlop(choice1, cpopt, cprj_dum, enlx(ib_start), gs_hamk, 0, eig(ib_start), &
              ndat, mpi_enreg, 1, paw_opt0, signs1, gsc_bk, tim_nonlop, &
              phi_bk, gvnlxc_bk, already_transposed=.False.)
         end if

       else
          ! FIXME This call seems to be needed for PAW. Strange as <G|S|psi> is already recomputed in pw_orthon!
          !call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, phi_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
          !                     eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, &
          !                     normalize=.False.)
          !                     !normalize=.True.)
       end if

     case (2)
       ! Consistent mode: update enlx_bx, eigens, residuals after orthogonalizalization.
       call getghc_eigresid(gs_hamk, npw, my_nspinor, ndat, phi_bk, ghc_bk, gsc_bk, mpi_enreg, prtvol, &
                            eig(ib_start), resid(ib_start), enlx(ib_start), residv_bk, gvnlxc_bk, &
                            normalize=.False.)
                            !normalize=.True.)
     case default
       MSG_BUG(sjoin("Wrong after_ortho:", itoa(after_ortho)))
     end select
   end do
   if (timeit) call cwtime_report(" recompute_eigens ", cpu, wall, gflops)

 else
   !if (prtvol)
   call wrtout(std_out, " Still far from convergence. Eigenvalues, residuals and NC enlx won't be recomputed.")
 end if

 ! Revert mixprec to previous status before returning.
 if (use_fft_mixprec) prev_mixprec = fftcore_set_mixprec(prev_mixprec)

 call yaml_write_dict('RMM-DIIS', "stats", diis%stats, std_out, with_iter_state=.False.)

 ! Final cleanup.
 ABI_FREE(lambda_bk)
 ABI_FREE(dots_bk)
 ABI_FREE(ghc_bk)
 ABI_FREE(residv_bk)
 ABI_FREE(kres_bk)
 ABI_FREE(gvnlxc_bk)
 call diis%free()

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

subroutine rmm_diis_push_iter(diis, iter, ndat, eig_bk, resid_bk, enlx_bk, phi_bk, residv_bk, gsc_bk, tag)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat
 real(dp),intent(in) :: eig_bk(ndat), resid_bk(ndat), enlx_bk(ndat)
 real(dp),intent(in) :: phi_bk(2, diis%npwsp*ndat), residv_bk(2, diis%npwsp*ndat), gsc_bk(2, diis%npwsp*ndat*diis%usepaw)
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
       if (diis%use_smat == 1) diis%smat(:, 0, 0, idat) = [one, zero]
     else
       diis%resmat(:, 0, 0, idat) = resid_bk(idat)
       if (diis%use_smat == 1) diis%smat(:, 0, 0, idat) = one
     end if
   end if
   !write(std_out, *)"res0", diis%resmat(:, 0, 0, idat)
   diis%step_type(iter, idat) = tag
   ibk = 1 + (idat - 1) * diis%npwsp
   call cg_zcopy(diis%npwsp, phi_bk(:,ibk), diis%chain_phi(:,:,iter,idat))
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

 diis%last_iter = iter
 ! FIXME: Temporarily disabled
 !ans = .False.; return
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
   fact = one !; if (abs(occ_bk(idat)) < diis%tol_occupied) fact = three
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
     fact = one; if (abs(occ_bk(idat)) < diis%tol_occupied) fact = ten
     if (sqrt(abs(resid)) < fact * accuracy_ene) then
       checks(idat) = 1; msg_list(idat) = 'resid < accuracy_ene'; cycle
     end if
   end if
 end do ! idat

 ! Depending on the accuracy_level either full block or a fraction of it must pass the test in order to exit.
 !nbocc = count(abs(occ_bk) < diis%tol_occupied)
 !nbocc_ok = count(checks /= 0 .and. abs(occ_bk) < diis%tol_occupied)
 !nbempty = ndat - nbocc
 !ans = nbocc_ok == nbocc

 nok = count(checks /= 0)
 if (diis%accuracy_level == 1) ans = nok >= 0.70_dp * ndat
 if (diis%accuracy_level == 2) ans = nok >= 0.80_dp * ndat
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
!!  Compute new eigenvalues, residuals, H|psi> and enlx from cg and gsc.
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
 real(dp),intent(out) :: residvecs(2, npw*my_nspinor*ndat)
 real(dp),intent(out) :: gvnlxc(2, npw*my_nspinor*ndat)
 logical,optional,intent(in) :: normalize

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: istwf_k, usepaw, cpopt, sij_opt, npwsp, me_g0, comm
 real(dp),parameter :: rdummy = zero
 real(dp) :: cpu, wall, gflops
 logical :: normalize_
!arrays
 real(dp) :: dots(2, ndat)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 normalize_ = .True.; if (present(normalize)) normalize_ = normalize
 npwsp = npw * my_nspinor
 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k; me_g0 = mpi_enreg%me_g0
 comm = mpi_enreg%comm_spinorfft; if (mpi_enreg%paral_kgb == 1) comm = mpi_enreg%comm_bandspinorfft

 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! NC normalization.
 if (usepaw == 0 .and. normalize_) call cgnc_normalize(npwsp, ndat, cg, istwf_k, me_g0, comm)

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
 if (usepaw == 1 .and. normalize_) call cgpaw_normalize(npwsp, ndat, cg, gsc, istwf_k, me_g0, comm)

 ! Compute new approximated eigenvalues, residual vectors and norms
 call cg_get_eigens(usepaw, istwf_k, npwsp, ndat, cg, ghc, gsc, eig, me_g0, comm)
 call cg_get_residvecs(usepaw, npwsp, ndat, eig, cg, ghc, gsc, residvecs)
 call cg_norm2g(istwf_k, npwsp, ndat, residvecs, resid, me_g0, comm)

 if (usepaw == 0) then
   ! Evaluate new enlx for NC.
   call cg_zdotg_zip(istwf_k, npwsp, ndat, option1, cg, gvnlxc, dots, me_g0, comm)
   enlx = dots(1,:)
 end if

 if (timeit) call cwtime_report(" getghc_eigresid", cpu, wall, gflops)

end subroutine getghc_eigresid
!!***

!!****f* m_rmm_diis/rmm_diis_new
!! NAME
!!  rmm_diis_new
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

type(rmm_diis_t) function rmm_diis_new(accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol) result(diis)

 integer,intent(in) :: accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol

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
 if (diis%use_smat == 1) then
   ABI_CALLOC(diis%smat, (diis%cplex, 0:max_niter, 0:max_niter, bsize)) ! <i|S|j>
 end if

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

 class(rmm_diis_t),intent(inout) :: diis

 ABI_SFREE(diis%hist_ene)
 ABI_SFREE(diis%hist_resid)
 ABI_SFREE(diis%hist_enlx)
 ABI_SFREE(diis%step_type)
 ABI_SFREE(diis%chain_phi)
 ABI_SFREE(diis%chain_sphi)
 ABI_SFREE(diis%chain_resv)
 ABI_SFREE(diis%resmat)
 ABI_SFREE(diis%smat)

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

subroutine rmm_diis_update_block(diis, iter, npwsp, ndat, lambda_bk, phi_bk, residv_bk, comm)

 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iter, npwsp, comm, ndat
 real(dp),intent(in) :: lambda_bk(ndat)
 real(dp),intent(inout) :: phi_bk(2, npwsp, ndat), residv_bk(2, npwsp, ndat)

!local variables
 integer,parameter :: master = 0
 integer :: cplex, ierr, nprocs, my_rank, idat !, ii !, ibk, iek
 real(dp) :: noise, cpu, wall, gflops
 !integer :: failed(ndat)
 real(dp),allocatable :: diis_eig(:), wmat1(:,:,:), wmat2(:,:,:), wvec(:,:,:), alphas(:,:)
 character(len=500) :: msg
 ! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 cplex = diis%cplex
 !failed = 0
 ABI_UNUSED(lambda_bk)

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 if (diis%use_smat == 1) then
   ! try_to_solve_eigproblem = .False. Numerically unstable.

   !do idat=1,ndat ! TODO
   ABI_MALLOC(diis_eig, (0:iter-1))
   ABI_MALLOC(wmat1, (cplex, 0:iter-1, 0:iter-1))
   ABI_MALLOC(wmat2, (cplex, 0:iter-1, 0:iter-1))
   wmat1 = diis%resmat(:, 0:iter-1, 0:iter-1, idat)
   wmat2 = diis%smat(:, 0:iter-1, 0:iter-1, idat)
   !do ii=0,iter-1; write(std_out, *) "diis_resmat:", wmat1(:,ii,:); end do
   !do ii=0,iter-1; write(std_out, *) "diis_smat:", wmat2(:,ii,:); end do
   ABI_CHECK(cplex == 2, "cplex 1 not coded")

   call xhegv_cplex(1, "V", "U", cplex, iter, wmat1, wmat2, diis_eig, msg, ierr)
   !write(std_out,*)"diis_eig:", diis_eig(0)
   !write(std_out,*)"RE diis_vec  :", wmat1(1,:,0)
   !write(std_out,*)"IMAG diis_vec:", wmat1(2,:,0)
   ABI_CHECK(ierr == 0, "xhegv returned ierr != 0")
   !if (ierr /= 0) then
   !  !call wrtout(std_out, sjoin("xhegv failed with:", msg, ch10, "at iter: ", itoa(iter), "exit iter_loop!"))
   !  ABI_FREE(diis_eig)
   !  ABI_FREE(wmat1)
   !  ABI_FREE(wmat2)
   !  !goto 10
   !end if
   !end do

   ! Take linear combination of chain_phi and chain_resv.
   call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), wmat1(:,:,0), phi_bk(:,:, idat))
   call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), wmat1(:,:,0), residv_bk(:,:,idat))

   ABI_FREE(wmat1)
   ABI_FREE(wmat2)
   ABI_FREE(diis_eig)
   !cycle

 else
   ! Solve system of linear equations.
   ! Only master works so that we are sure we have the same solution.
   ABI_CALLOC(wvec, (cplex, 0:iter, ndat))

   if (my_rank == master) then
     ABI_MALLOC(wmat1, (cplex, 0:iter, 0:iter))

     do idat=1,ndat
       wvec(1, iter, idat) = -one
       wmat1 = zero
       wmat1(1,:,iter) = -one
       wmat1(1,iter,:) = -one
       wmat1(1,iter,iter) = zero
       wmat1(:,0:iter-1, 0:iter-1) = diis%resmat(:, 0:iter-1, 0:iter-1, idat)

       !if (ierr /= 0) then
       !  write(std_out, *)"iter:", iter, "cplex:", cplex
       !  do ii=0,iter
       !    !write(std_out, sjoin("(", itoa(2*iter), ", (es14.6)")) "wmat1:", wmat1(:, ii, :)
       !    write(std_out, "(a, *(es14.6))") " wmat1:", wmat1(:, ii, :)
       !  end do
       !end if

       call xhesv_cplex("U", cplex, iter+1, 1, wmat1, wvec(:,:,idat), msg, ierr)
       ABI_CHECK(ierr == 0, msg)
       !failed(idat) = ierr

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
 end if ! use_smat

 ! Master broadcasts data.
 if (nprocs > 1) then
   call xmpi_bcast(wvec, master, comm, ierr)
   !call xmpi_bcast(failed, master, comm, ierr)
 end if

 do idat=1,ndat

   if (cplex == 2) then
     call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), wvec(:,:,idat), phi_bk(:,:,idat))
     call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), wvec(:,:,idat), residv_bk(:,:,idat))
   else
     !ABI_CALLOC(alphas, (2, 0:iter))
     !alphas(1,:) = wvec(1,:,idat)
     !call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), alphas, phi_bk(:,:,idat))
     !call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), alphas, residv_bk(:,:,idat))
     !ABI_FREE(alphas)

     ! Use DGEMV
     ABI_MALLOC(alphas, (1, 0:iter))
     alphas(1,:) = wvec(1,:,idat)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_phi(:,:,:,idat), 2*npwsp, alphas, 1, zero, phi_bk(:,:,idat), 1)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_resv(:,:,:,idat), 2*npwsp, alphas, 1, zero, residv_bk(:,:,idat), 1)
     ABI_FREE(alphas)
   end if
 end do

 ABI_FREE(wvec)

 if (timeit) call cwtime_report(" update_block", cpu, wall, gflops)

end subroutine rmm_diis_update_block
!!***

!!****f* m_rmm_diis/rmm_diis_eval_mats
!! NAME
!!  rmm_diis_eval_mats
!!
!! FUNCTION
!!  Compute matrix elements required by DIIS method.
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
 real(dp) :: dotr, doti, cpu, wall, gflops
! *************************************************************************

 nprocs = xmpi_comm_size(comm)
 option = 2; if (diis%cplex == 1) option = 1

 if (timeit) call cwtime(cpu, wall, gflops, "start")

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

     ! <i|S|j> assume normalized wavefunctions for ii == iter. Note rescaling by 1/np
     if (diis%use_smat == 1) then
       if (ii == iter) then
         dotr = one / nprocs; doti = zero
       else
         if (diis%usepaw == 0) then
           call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option, &
                          diis%chain_phi(:,:,ii,idat), diis%chain_phi(:,:,iter,idat), me_g0, xmpi_comm_self)
         else
           call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option, &
                          diis%chain_phi(:,:,ii,idat), diis%chain_sphi(:,:,iter,idat), me_g0, xmpi_comm_self)
         end if
       end if
       if (diis%cplex == 2) then
         diis%smat(:, ii, iter, idat) = [dotr, doti]
       else
         diis%smat(1, ii, iter, idat) = dotr
       end if
    end if
   end do ! ii

   if (nprocs > 1) then
     call xmpi_sum(diis%resmat(:,0:iter,iter,idat), comm, ierr)
     if (diis%use_smat == 1) call xmpi_sum(diis%smat(:,0:iter,iter, idat), comm, ierr)
   endif
   !if (diis%prtvol == -level) then
   !  write(std_out,*)"iter, idat, resmat:", iter, idat, diis%resmat(:,0:iter,iter,idat)
   !  !write(std_out,*)"iter, idat smat:", iter, idat, diis%smat(:,0:iter,iter,idat)
   !end if

 end do ! idat

 !if (nprocs > 1) then
 !  call xmpi_sum(diis%resmat(:,0:iter,iter,1:ndat), comm, ierr)
 !  call xmpi_sum(diis%smat(:,0:iter,iter, 1:ndat), comm, ierr)
 !endif

 if (timeit) call cwtime_report(" eval_mats", cpu, wall, gflops)

end subroutine rmm_diis_eval_mats
!!***

!!****f* m_rmm_diis/rmm_diis_rollback
!! NAME
!!  rmm_diis_rollback
!!
!! FUNCTION
!!  High energy states may have larger residuals at the end of the iter_loop especially if we reduce the
!!  number of iterations. Here we select the trial states in the chain with smaller resid.
!!  We are allowed to do so because we are still othogonalization-free.
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

subroutine rmm_diis_rollback(diis, npwsp, ndat, phi_bk, gsc_bk, residv_bk, eig, resid, enlx, comm)

!Arguments ------------------------------------
 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: npwsp, comm, ndat
 real(dp),intent(inout) :: phi_bk(2, npwsp, ndat), gsc_bk(2, npwsp, ndat*diis%usepaw), residv_bk(2, npwsp, ndat)
 real(dp),intent(inout) :: eig(ndat), resid(ndat), enlx(ndat)

!local variables
 integer,parameter :: master = 0
 integer :: nprocs, my_rank, idat, iter, ierr, ilast, take_iter(ndat)
 real(dp) :: cpu, wall, gflops
! *************************************************************************

 ! FIXME: Temporarily disabled
 !return
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 take_iter = -1
 ilast = diis%last_iter

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 if (my_rank == master) then
   do idat=1,ndat
     iter = imin_loc(diis%hist_resid(0:ilast, idat)) - 1
     ! Move stuff only if it's really worth it.
     if (iter /= ilast .and. &
         diis%hist_resid(iter, idat) < third * diis%hist_resid(ilast, idat)) take_iter(idat) = iter
   end do
 end if
 if (nprocs > 1) call xmpi_bcast(take_iter, master, comm, ierr)

 if (any(take_iter /= -1)) then
   do idat=1,ndat
     iter = take_iter(idat); if (iter == -1) cycle
     diis%step_type(iter, idat) = trim(diis%step_type(iter, idat)) // "*"
     eig(idat) = diis%hist_ene(iter, idat)
     resid(idat) = diis%hist_resid(iter, idat)
     enlx(idat) = diis%hist_enlx(iter, idat)
     ! Copy |psi>, S|psi> and residual.
     call cg_zcopy(npwsp, diis%chain_phi(:,:,iter,idat), phi_bk(:,:,idat))
     call cg_zcopy(npwsp, diis%chain_resv(:,:,iter,idat), residv_bk(:,:,idat))
     if (diis%usepaw == 1) call cg_zcopy(npwsp, diis%chain_sphi(:,:,iter,idat), gsc_bk(:,:,idat))
     call diis%stats%increment("rollback", 1)
   end do
 end if

 if (timeit) call cwtime_report(" diis_rollback", cpu, wall, gflops)

end subroutine rmm_diis_rollback
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

!!****f* ABINIT/rmm_diis
!! NAME
!! rmm_diis
!!
!! FUNCTION
!! This routine updates the wave functions at a given k-point, using the RMM-DIIS method.
!!
!! INPUTS
!!  istep,ikpt,isppol
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  my_nspinor=number of plane waves at this k point
!!
!! OUTPUT
!!  If usepaw==1:
!!    gsc_all(2,*)=<g|s|c> matrix elements (s=overlap)
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine subspace_rotation(gs_hamk, dtset, mpi_enreg, nband, npw, my_nspinor, eig, cg, gsc_all, evec)

!Arguments ------------------------------------
 integer,intent(in) :: nband, npw, my_nspinor
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*my_nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*my_nspinor*nband*dtset%usepaw)
 real(dp),intent(out) :: eig(nband)
 real(dp),allocatable,intent(out) :: evec(:,:,:)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, tim_getghc = 0, use_subovl0 = 0
 integer :: ig, ig0, ib, ierr, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, paral_kgb
 integer :: iband, cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, usepaw !, nspinor
 integer :: me_g0, cplex, comm_bandspinorfft
 real(dp) :: cpu, wall, gflops
 real(dp),parameter :: rdummy = zero
 !character(len=500) :: msg
!arrays
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp) :: subovl(use_subovl0)
 real(dp),allocatable :: subham(:), h_ij(:,:,:), ghc_bk(:,:), gvnlxc_bk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_bk(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; paral_kgb = mpi_enreg%paral_kgb
 me_g0 = mpi_enreg%me_g0
 comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
 !nspinor = dtset%nspinor
 npwsp = npw * my_nspinor

 ! =======================================
 ! Apply H to input cg to compute <i|H|j>
 ! =======================================
 gsc_bk => fake_gsc_bk
 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! Treat states in groups of bsize bands.
 bsize = 8  ! default value for paral_kgb = 0
 if (paral_kgb == 1) bsize = mpi_enreg%nproc_band * mpi_enreg%bandpp
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 cplex = 2; if (istwf_k == 2) cplex = 1
 ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
 ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))
 ABI_CALLOC(h_ij, (cplex, nband, nband))

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize; ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)
   if (usepaw == 1) gsc_bk => gsc_all(1:2,igs:ige)

   if (paral_kgb == 0) then
     call getghc(cpopt, cg(:,igs:ige), cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                 rdummy, mpi_enreg, ndat, dtset%prtvol, sij_opt, tim_getghc, type_calc0)
   else
     call prep_getghc(cg(:,igs:ige), gs_hamk, gvnlxc_bk, ghc_bk, gsc_bk, rdummy, ndat, &
                      mpi_enreg, dtset%prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
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
 call xmpi_sum(subham, comm_bandspinorfft, ierr)
 !do it=1,nband*(nband+1); write(std_out,*)"subham:", it, subham(it); end do

 ! ========================
 ! Subspace diagonalization
 ! ========================
 ABI_MALLOC(evec, (2, nband, nband))
 mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 call subdiago(cg, eig, evec, gsc_all, 0, 0, istwf_k, mcg, mgsc, nband, npw, my_nspinor, paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)

 ABI_FREE(subham)
 ABI_FREE(ghc_bk)
 ABI_FREE(gvnlxc_bk)

 if (timeit) call cwtime_report(" subspace rotation", cpu, wall, gflops)

end subroutine subspace_rotation
!!***

end module m_rmm_diis
!!***
