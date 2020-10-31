!!****m* ABINIT/m_rmm_diis
!! NAME
!!  m_rmm_diis
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (MG)
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
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 use m_fftcore,       only : fftcore_set_mixprec
 !use m_fock,         only : fock_set_ieigen, fock_set_getghc_call

 implicit none

 private
!!***

 public :: rmm_diis
!!***

 type,private :: rmm_diis_t

   integer :: accuracy_level
   integer :: usepaw
   integer :: istwf_k
   integer :: cplex
   integer :: bsize
   integer :: max_niter
   integer :: npwsp
   integer :: prtvol
   integer :: last_iter
   integer :: use_smat = 0

   real(dp) :: tol_occupied

   type(pair_list) :: stats

   real(dp),allocatable :: hist_ene(:,:)
   real(dp),allocatable :: hist_resid(:,:)
   real(dp),allocatable :: hist_enlx(:,:)
   character(len=7),allocatable :: step_type(:,:)
   ! (0:max_niter+2, bsize)
   ! 0 is the initial step, then DIIS iterations (whose number may depend on the block)
   ! followed by an optional trial step and the computation of eigens after ortho.

   real(dp),allocatable :: resmat(:,:,:,:)
   real(dp),allocatable :: smat(:,:,:,:)
   ! (2, 0:max_niter, 0:max_niter, bsize))

   real(dp),allocatable :: chain_phi(:,:,:,:)
   real(dp),allocatable :: chain_sphi(:,:,:,:)
   real(dp),allocatable :: chain_resv(:,:,:,:)
   ! (2, npwsp, 0:max_niter, bsize))

 contains
   procedure :: free => rmm_diis_free
   procedure :: update_block => rmm_diis_update_block
   procedure :: eval_mats => rmm_diis_eval_mats
   procedure :: exit_iter => rmm_diis_exit_iter
   procedure :: print_block => rmm_diis_print_block
   procedure :: rollback => rmm_diis_rollback
   procedure :: push_hist => rmm_diis_push_hist
   procedure :: push_iter => rmm_diis_push_iter

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
!! This routine updates the wave functions at a given k-point, using the RMM-DIIS method.
!!
!! INPUTS
!!  istep,ikpt,isppol
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  nspinor=number of plane waves at this k point
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!  resid(nband)=residuals for each states
!!  If usepaw==1:
!!    gsc_all(2,*)=<g|s|c> matrix elements (s=overlap)
!!  If usepaw==0
!!    enlx(nband)=contribution from each band to nonlocal psp + potential Fock ACE part of total energy, at this k-point
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!  rmm_diis_status(2): Status of the eigensolver.
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
                    mpi_enreg, nband, npw, nspinor, resid, rmm_diis_status)

!Arguments ------------------------------------
 integer,intent(in) :: istep, ikpt, isppol, nband, npw, nspinor
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband), resid(nband)
 real(dp),intent(in) :: occ(nband), kinpw(npw)
 real(dp),intent(out) :: eig(nband)
 integer,intent(inout) :: rmm_diis_status(2)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0, use_subovl0 = 0
 integer :: ig, ig0, ib, ierr, prtvol, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, idat, paral_kgb !, comm
 integer :: iband, cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, optekin, usepaw, iter, max_niter, max_niter_block
 integer :: me_g0, nb_pocc, jj, kk, it, ibk, iek, cplex, ortalgo, accuracy_level, raise_acc, prev_mixprec !ii, ld1, ld2,
 integer :: comm_bandspinorfft, prev_accuracy_level, ncalls_with_this_accuracy !, comm_spinorfft, comm_fft
 logical :: end_with_trial_step, recompute_lambda, recompute_eigresid_after_ortho, first_call, do_rollback
 real(dp),parameter :: rdummy = zero
 real(dp) :: accuracy_ene, cpu, wall, gflops, max_res_pocc, tol_occupied !, dotr, doti
 !character(len=500) :: msg
 character(len=6) :: tag
!arrays
 real(dp) :: tsec(2)
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp) :: subovl(use_subovl0)
 real(dp),allocatable :: evec(:,:), subham(:), h_ij(:,:,:)
 real(dp),allocatable :: ghc_bk(:,:), gsc_bk(:,:), gvnlxc_bk(:,:), lambda_bk(:), dots_bk(:,:)
 real(dp),allocatable :: residv_bk(:,:), kres_bk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: ptr_gsc_bk(:,:), phi_bk(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(rmm_diis_t) :: diis
 type(yamldoc_t) :: ydoc

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; prtvol = dtset%prtvol
 paral_kgb = mpi_enreg%paral_kgb
 me_g0 = mpi_enreg%me_g0 !; comm_fft = mpi_enreg%comm_fft; comm_spinorfft = mpi_enreg%comm_spinorfft;
 !comm = mpi_enreg%comm_spinorfft; if (mpi_enreg%paral_kgb == 1) comm = mpi_enreg%comm_bandspinorfft
 comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
 npwsp = npw * nspinor

 ! =======================================
 ! Apply H to input cg to compute <i|H|j>
 ! =======================================
 ptr_gsc_bk => fake_gsc_bk
 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! Treat states in groups of bsize bands to be able to call zgemm.
 bsize = 8  ! default for paral_kgb = 0 !if (dtset%userib /= 0) bsize = abs(dtset%userib)
 if (paral_kgb == 1) bsize = mpi_enreg%nproc_band * mpi_enreg%bandpp
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 cplex = 2; if (istwf_k == 2) cplex = 1
 ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
 ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))
 ABI_CALLOC(h_ij, (cplex, nband, nband))

 call timab(1633, 1, tsec) !"rmm_diis:build_hij
 if (timeit) call cwtime(cpu, wall, gflops, "start")

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize
   ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

   if (usepaw == 1) ptr_gsc_bk => gsc_all(1:2,igs:ige)

   if (paral_kgb == 0) then
     call getghc(cpopt, cg(:,igs:ige), cprj_dum, ghc_bk, ptr_gsc_bk, gs_hamk, gvnlxc_bk, &
                 rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
   else
     call prep_getghc(cg(:,igs:ige), gs_hamk, gvnlxc_bk, ghc_bk, ptr_gsc_bk, rdummy, ndat, &
                      mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
   end if
   !ghc_all(:,igs:ige) = ghc_bk(:,igs:ige)

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

 if (timeit) call cwtime_report(" build_hij", cpu, wall, gflops)

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
 call timab(1633, 2, tsec) !"rmm_diis:build_hij

 !do it=1,nband*(nband+1)
 !  write(std_out,*)"subham:", it, subham(it)
 !end do

 ! ========================
 ! Subspace diagonalization
 ! ========================
 call timab(585, 1, tsec) !"vtowfk(subdiago)"
 ABI_MALLOC(evec, (2*nband, nband))
 mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 call subdiago(cg, eig, evec, gsc_all, 0, 0, istwf_k, mcg, mgsc, nband, npw, nspinor, paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)

 call timab(585, 2, tsec)
 if (timeit) call cwtime_report(" subdiago", cpu, wall, gflops)

 ABI_FREE(subham)
 ABI_FREE(evec)

 ! =================
 ! Prepare DIIS loop
 ! =================
 ! accuracy_level
 !  1: Used at the beginning of the SCF cycle. Use loosy convergence criteria in order
 !     to reduce the number of H|psi> applications as much as possible so that we can mix densities/potentials.
 !     Move to the next level after Max 15 iterations.
 !  2: Intermediate step. Decrease convergence criteria in order to perform more wavefuction iterations.
 !     Allow for incosistent data in rediduals and Vnl matrix elements.
 !     Move to the next level after Max 25 iterations.
 !  3: Approaching convergence. Use stricter convergence criteria.
 !     Move to the next level after Max 25 iterations.
 !  4: Ultimate precision. Try to reach the same accuracy as the other eigenvalue solvers.
 !     This means: using similar convergence criteria as in the other solvers.

 ! We are now allowed to decrease accuracy_level during the SCF cycle this means that the routine
 ! must receive a table diis_accuracy(2, nkpt, nsppol)
 ! The first entry gives the previous accuracy.
 ! The second entry gives the number of iterations already performed with this level.
 !
 if (all(rmm_diis_status == 0)) then
   ! This is the first time we call rmm_diis.
   prev_accuracy_level = 1; ncalls_with_this_accuracy = 0
   first_call = .True.
 else
   prev_accuracy_level = rmm_diis_status(1); ncalls_with_this_accuracy = rmm_diis_status(2)
   first_call = .False.
 end if
 ! Decide whether we should move to the next level
 raise_acc = 0
 if (prev_accuracy_level == 1 .and. ncalls_with_this_accuracy >= 15) raise_acc = 2
 if (prev_accuracy_level == 2 .and. ncalls_with_this_accuracy >= 25) raise_acc = 3
 if (prev_accuracy_level == 3 .and. ncalls_with_this_accuracy >= 25) raise_acc = 4
 raise_acc = max(raise_acc, prev_accuracy_level)
 !print *, "rmm_diis_status:", rmm_diis_status
 !print *, "rmm_prev_acc:", prev_accuracy_level, "rmm_raise_acc:", raise_acc

 ! Define tolerance for occupied states on the basis of prev_accuracy_level.
 ! and compute max of residuals for this set.
 tol_occupied = tol3
 if (any(prev_accuracy_level == [1])) tol_occupied = tol2
 nb_pocc = count(occ > tol_occupied)
 max_res_pocc = maxval(resid(1:nb_pocc))

 ! Defin accuracy_level of this iteration.
 accuracy_level = 1
 if (max_res_pocc < tol6) accuracy_level = 2
 if (max_res_pocc < tol10) accuracy_level = 3
 if (max_res_pocc < tol15) accuracy_level = 4
 accuracy_level = max(prev_accuracy_level, accuracy_level, raise_acc)
 if (istep == 1) accuracy_level = 2  ! FIXME: Differenciate between restart or rmm_diis - 3.
 if (first_call) accuracy_level = 1
 if (usepaw == 1) accuracy_level = 3 ! # FIXME

 ! Update rmm_diis_status. Reset number of calls if we moved to a new accuracy_level.
 rmm_diis_status(1) = accuracy_level
 if (accuracy_level /= prev_accuracy_level) rmm_diis_status(2) = 0
 rmm_diis_status(2) = rmm_diis_status(2) + 1

 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1 ! optekin = 0
 !write(std_out,*)"optekin:", optekin

 ! nline - 1 DIIS steps (usually 3), then end with trial step (optional, see below)
 max_niter = max(dtset%nline - 1, 1)
 !if (accuracy_ene == 1) max_niter = max(dtset%nline - 2, 1)
 if (accuracy_ene >= 4) max_niter = dtset%nline

 ! Define accuracy_ene for SCF
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

 recompute_eigresid_after_ortho = diis%accuracy_level >= 2
 !recompute_eigresid_after_ortho = diis%accuracy_level >= 3
 !recompute_eigresid_after_ortho = .False.
 if (usepaw == 1) recompute_eigresid_after_ortho = .True. ! # FIXME

 !if (dtset%fftcore_mixprec == 1 .and. accuracy_level <= 2)
 !if (accuracy_level <= 2) prev_mixprec = fftcore_set_mixprec(1)

 diis = rmm_diis_new(accuracy_level, usepaw, istwf_k, npwsp, max_niter, bsize, prtvol)
 diis%tol_occupied = tol_occupied
 call wrtout(std_out, sjoin(" Using Max", itoa(max_niter), "RMM-DIIS iterations + optional final trial step."))
 call wrtout(std_out, sjoin( &
   " Number of blocks:", itoa(nblocks), ", nb_pocc:", itoa(nb_pocc)))
 call wrtout(std_out, sjoin( &
   " Max_input_resid_pocc", ftoa(max_res_pocc), "accuracy_level:", itoa(accuracy_level), &
   ", accuracy_ene: ", ftoa(accuracy_ene)))
 call timab(1634, 1, tsec) ! "rmm_diis:band_opt"

 ABI_MALLOC(lambda_bk, (bsize))
 ABI_MALLOC(dots_bk, (2, bsize))
 ABI_MALLOC(gsc_bk, (2, npwsp*bsize*usepaw))
 ABI_MALLOC(residv_bk, (2, npwsp*bsize))
 ABI_MALLOC(kres_bk, (2, npwsp*bsize))

 ! We loop over nblocks, each block has ndat states.
 ! - Convergence behaviour may depend on bsize as branches are taken according to
 !   the status of all bands in the block.
 ! - Using gemm_nonlop leads to a significant speedup when applying Vnl.
 ! TODO: Transpose only once per block and then work with already_transposed = .True.

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize
   ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

   ! Reduce number of niter iterations for empty states.
   ! TODO: Don't reduce niter if MD
   !if dtset%
   !if dtset%occopt == 2
   max_niter_block = max_niter
   if (all(occ(ib_start:ib_stop) < diis%tol_occupied)) max_niter_block = max(1 + max_niter / 2, 2)

   ! Compute H |phi_0> using cg from subdiago. Blocked call.
   ! Alternatively, one can compute the residuals before the subspace rotation and then rotate
   ! to save one call to getghc_eigresid.

   ! Point bands in the block.
   phi_bk => cg(:,igs:ige)
   if (usepaw == 1) ptr_gsc_bk => gsc_all(1:2,igs:ige)

   call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                        eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, gvnlxc_bk, normalize=.False.)
   !write(std_out,*)"eig", eig(ib_start:ib_stop), resid(ib_start:ib_stop), enlx(ib_start:ib_stop)

   ! Save <R0|R0> and <phi_0|S|phi_0>, |phi_0>, |S phi_0>. Assume input cg are already S-normalized.
   call diis%push_iter(0, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), phi_bk, residv_bk, ptr_gsc_bk, "SDIAG")

   if (timeit) call cwtime_report(" first getghc_eigresid ", cpu, wall, gflops)

   ! Precondition |R_0>, output in kres_bk = |K R_0>
   call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
   call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)
   if (timeit) call cwtime_report(" first cg_precon ", cpu, wall, gflops)

   ! Compute H |K R_0>
   if (paral_kgb == 0) then
     call getghc(cpopt, kres_bk, cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                 rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
   else
     call prep_getghc(kres_bk, gs_hamk, gvnlxc_bk, ghc_bk, ptr_gsc_bk, rdummy, ndat, &
                      mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
   end if

   ! Compute residuals: (H - e_0 S) |K R_0>
   call cg_residvecs(usepaw, npwsp, ndat, eig(ib_start:), kres_bk, ghc_bk, gsc_bk, residv_bk)
   if (timeit) call cwtime_report(" cg_residvecs ", cpu, wall, gflops)

   ! Line minimization with preconditioned steepest descent:
   !
   !    |phi_1> = |phi_0> + lambda |K R_0>
   !
   ! where lambda minimizes the residual (not the Rayleigh quotient as in Kresse's paper).
   !
   !    lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
   !
   call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_bandspinorfft)
   if (timeit) call cwtime_report(" cg_norm2g ", cpu, wall, gflops)

   dots_bk = zero
   do idat=1,ndat
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     call dotprod_g(dots_bk(1,idat), dots_bk(2,idat), istwf_k, npwsp, option1, &
                    diis%chain_resv(:,:,0,idat), residv_bk(:,jj), me_g0, xmpi_comm_self)
   end do
   call xmpi_sum(dots_bk, comm_bandspinorfft, ierr)

   do idat=1,ndat
     !write(std_out,*)"Computing lambda with:", -dots_bk(1,idat), lambda_bk(idat)
     lambda_bk(idat) = limit_lambda(-dots_bk(1,idat) / lambda_bk(idat))
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     phi_bk(:,jj:kk) = diis%chain_phi(:,:,0,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
   end do
   if (timeit) call cwtime_report(" KR0 ", cpu, wall, gflops)

   iter_loop: do iter=1,max_niter_block

     if (iter > 1) then
       ! Solve DIIS equations to get phi_bk for iter > 1
       call diis%update_block(iter, npwsp, ndat, lambda_bk, phi_bk, residv_bk, comm_bandspinorfft)

       ! Precondition residual, output in kres_bk.
       call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
       call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)

       ! Compute phi_bk with the lambda(ndat) obtained at iteration #1
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)
     end if

     ! Compute H |phi_now> and evaluate new enlx for NC.
     call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, gvnlxc_bk, normalize=.True.)

     ! Take another step with the same lambda
     !if (iter == 1) then
     !  call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
     !  call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)
     !  do idat=1,ndat
     !    jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     !    phi_bk(:,jj:kk) = diis%chain_phi(:,:,iter-1,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
     !  end do
     !end if

     ! Store new residual.
     call diis%push_iter(iter, ndat, eig(ib_start), resid(ib_start), enlx(ib_start), phi_bk, residv_bk, ptr_gsc_bk, "DIIS")

     ! === CHECK FOR CONVERGENCE ====
     if (diis%exit_iter(iter, ndat, max_niter_block, occ(ib_start:), accuracy_ene, &
                        dtset, comm_bandspinorfft)) exit iter_loop

     ! Compute <R_i|R_j> and <i|S|j> for j=iter
     if (iter /= max_niter_block) call diis%eval_mats(iter, ndat, me_g0, comm_bandspinorfft)
   end do iter_loop

   ! End with trial step but only if we performed all the iterations (i.e. no exit from diis%exit)
   ! Since we operate on blocks of bands, all the states in the block will receive the same treatment.
   ! This means that one can observe a (hopefully) slightly different convergence behaviour depending on bsize.
   !
   end_with_trial_step = .True.
   if (timeit) call cwtime_report(" iterloop ", cpu, wall, gflops)

   if (end_with_trial_step) then
     call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
     call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, kinpw, kres_bk, me_g0, comm_bandspinorfft)

     recompute_lambda = diis%accuracy_level >= 4
     ! Recomputing the preconditioned-direction and the new lambda that minimizes the residual
     ! is more expensive but more accurate.
     ! Do it only for the occupied states plus a buffer
     ! if these bands still far from convergence. In all the other cases, we reuse the previous lambda.
     !recompute_lambda = .False.
     !recompute_lambda = .True.

     if (recompute_lambda) then
       tag = "NEWLAM"
       !write(std_out, *)" Performing last trial step with new computation of lambda"
       ! Final preconditioned steepest descent with new lambda.
       ! Compute H |K R_0>

       if (paral_kgb == 0) then
         call getghc(cpopt, kres_bk, cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                     rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)
       else
         call prep_getghc(kres_bk, gs_hamk, gvnlxc_bk, ghc_bk, ptr_gsc_bk, rdummy, ndat, &
                          mpi_enreg, prtvol, sij_opt, cpopt, cprj_dum, already_transposed=.False.)
       end if

       ! Compute residual: (H - e_0 S) |K R_0>
       call cg_residvecs(usepaw, npwsp, ndat, eig(ib_start:), kres_bk, ghc_bk, gsc_bk, residv_bk)
       call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_bandspinorfft)

       it = diis%last_iter
       do idat=1,ndat
         jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
         call dotprod_g(dots_bk(1,idat), dots_bk(2,idat), istwf_k, npwsp, option1, &
                        diis%chain_resv(:,:,it,idat), residv_bk(:,jj), me_g0, xmpi_comm_self)
       end do
       call xmpi_sum(dots_bk, comm_bandspinorfft, ierr)

       do idat=1,ndat
         lambda_bk(idat) = limit_lambda(-dots_bk(1,idat) / lambda_bk(idat))
         jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
         phi_bk(:,jj:kk) = diis%chain_phi(:,:,it,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
       end do

     else
       ! Reuse previous values of lambda.
       tag = "FIXLAM"
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)
     endif

     if (.not. recompute_eigresid_after_ortho .or. usepaw == 1) then ! FIXME
       ! Finally recompute eig, resid, enlx (and ptr_gsc_bk if PAW).
       call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                            eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, gvnlxc_bk, normalize=.True.)

       ! Insert the last results in the history just for printing purposes and increment last_iter.
       call diis%push_hist(diis%last_iter + 1, ndat, eig, resid, enlx, tag)
       if (timeit) call cwtime_report(" last_trial_step ", cpu, wall, gflops)

       !TODO: For PAW one should recompute S

     end if

   end if ! end_with_trial_step

   do_rollback = .False.
   if (do_rollback) then
     call diis%rollback(npwsp, ndat, phi_bk, ptr_gsc_bk, eig(ib_start), resid(ib_start), enlx(ib_start), comm_bandspinorfft)
   end if
   if (prtvol == -level) call diis%print_block(ib_start, ndat, istep, ikpt, isppol)

   ! wavefunction block. phi_bk => cg(:,igs:ige) has been updated.
 end do ! iblock

 call timab(1634, 2, tsec) !"rmm_diis:band_opt"
 if (timeit) call cwtime_report(" rmm_diis:band_opt", cpu, wall, gflops)

 call timab(583,1,tsec) ! "vtowfk(pw_orthon)"
 ortalgo = mpi_enreg%paral_kgb
 !ortalgo = 3
 if (prtvol > 0) call wrtout(std_out, " Calling pw_orthon to orthonormalize bands.")
 call pw_orthon(0, 0, istwf_k, mcg, mgsc, npwsp, nband, ortalgo, gsc_all, gs_hamk%usepaw, cg, &
                mpi_enreg%me_g0, mpi_enreg%comm_bandspinorfft)
 call timab(583,2,tsec)
 if (timeit) call cwtime_report(" pw_orthon ", cpu, wall, gflops)

 ! Recompute eigenvalues, residuals, and enlx(ndat)  after orthogonalization.
 ! This step is important to improve the convergence of the total energy
 ! but we try to avoid it at the beginning of the SCF cycle
 ! TODO: Rotate enlx using output of Cholesky decomposition but this also means that rollback must be disabled.
 ! and full matrix Vnl_ij is needed since:
 !
 !  if (usepaw == 0) then
 !    ABI_MALLOC(totvnlx, (2*nband_k,nband_k))
 !    !call cg_hrotate_and_get_diag(istwf_k, nband_k, totvnlx, evec, enlx_k)
 !    ABI_FREE(totvnlx)
 !  end if
 !

 if (recompute_eigresid_after_ortho) then
   call wrtout(std_out, " Recomputing eigenvalues and residues after orthogonalization.")
   do iblock=1,nblocks
     igs = 1 + (iblock - 1) * npwsp * bsize
     ige = min(iblock * npwsp * bsize, npwsp * nband)
     ndat = (ige - igs + 1) / npwsp
     ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

     if (usepaw == 1) ptr_gsc_bk => gsc_all(1:2,igs:ige)

     call getghc_eigresid(gs_hamk, npw, nspinor, ndat, cg(:,igs), ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, gvnlxc_bk, normalize=.False.)

     ! Insert the last results in the history just for printing purposes and increment last_iter.
     call diis%push_hist(diis%last_iter + 1, ndat, eig, resid, enlx, "ORTH")
   end do
   if (timeit) call cwtime_report(" recompute_eigens ", cpu, wall, gflops)
 else
   call wrtout(std_out, " Still far from convergenve. Eigenvalues, residuals and NC enlx won't be recomputed.")
 end if

 !call yaml_write_dict('RMM-DIIS', "stats", dict%stats, std_out, with_iter_state=.False.)
 ydoc = yamldoc_open("RMM-DIIS", with_iter_state=.False.)
 call ydoc%add_dict("stats", diis%stats)
 call ydoc%write_and_free(std_out)
 call diis%free()

 ! Final cleanup
 ABI_FREE(lambda_bk)
 ABI_FREE(dots_bk)
 ABI_FREE(ghc_bk)
 ABI_FREE(gsc_bk)
 ABI_FREE(residv_bk)
 ABI_FREE(kres_bk)
 ABI_FREE(gvnlxc_bk)

 !if (accuracy_level <= 2) prev_mixprec = fftcore_set_mixprec(prev_mixprec)

contains

 real(dp) function limit_lambda(lambda) result(new_lambda)

  real(dp),intent(in) :: lambda
  real(dp) :: max_lambda, min_lambda

  !write(std_out, *)"input lambda:", lambda
  new_lambda = lambda
  return
  if (accuracy_level == 1) return
  !if (accuracy_level == 3) return

  ! restrict the value of abs(lambda) in [0.1, 1.0]
  max_lambda = one
  min_lambda = tol1
  if (abs(lambda) > max_lambda) then
    new_lambda = sign(max_lambda, lambda)
  else if (abs(lambda) < min_lambda) then
    !if (abs(lambda) > tol12)  then
      new_lambda = sign(min_lambda, lambda)
    !else
    !  new_lambda = tol12
    !end if
  end if
  !write(std_out, *)"new lambda:", new_lambda

end function limit_lambda

end subroutine rmm_diis
!!***

!!****f* m_rmm_diis/rmm_diis_push_hist
!! NAME
!!  rmm_diis_push_hist
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

subroutine rmm_diis_push_hist(diis, iter, ndat, eig_bk, resid_bk, enlx_bk, tag)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat
 real(dp),intent(in) :: eig_bk(ndat), resid_bk(ndat), enlx_bk(ndat)
 character(len=*),intent(in) :: tag

 diis%last_iter = iter
 diis%hist_ene(iter, 1:ndat) = eig_bk
 diis%hist_resid(iter, 1:ndat) = resid_bk
 diis%hist_enlx(iter, 1:ndat) = enlx_bk
 diis%step_type(iter, 1:ndat) = tag

end subroutine rmm_diis_push_hist
!!***

!!****f* m_rmm_diis/rmm_diis_push_iter
!! NAME
!!  rmm_diis_push_iter
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

subroutine rmm_diis_push_iter(diis, iter, ndat, eig_bk, resid_bk, enlx_bk, phi_bk, residv_bk, gsc_bk, tag)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat
 real(dp),intent(in) :: eig_bk(ndat), resid_bk(ndat), enlx_bk(ndat)
 real(dp),intent(in) :: phi_bk(2, diis%npwsp*ndat), residv_bk(2, diis%npwsp*ndat), gsc_bk(2, diis%npwsp*ndat*diis%usepaw)
 character(len=*),intent(in) :: tag

 integer :: idat, ibk

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

 integer,parameter :: master = 0
 integer :: idat, ierr, nok, checks(ndat)
 real(dp) :: resid, deltae, deold , fact
 character(len=50) :: msg_list(ndat)

 diis%last_iter = iter; if (xmpi_comm_rank(comm) /= master) goto 10

 ! Tolerances depend on accuracy_level and occupation of the state.
 checks = 0

 do idat=1,ndat
   resid = diis%hist_resid(iter, idat)
   deold = diis%hist_ene(1, idat) - diis%hist_ene(0, idat)
   deltae = diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat)

   ! Relative criterion on eigenvalue differerence.
   ! Abinit default in the CG part is 0.005 that is really to low (it's 0.3 in V).
   ! Here we increase it depending whether the state is occupied or empty
   fact = one; if (abs(occ_bk(idat)) < diis%tol_occupied) fact = five
   if (diis%accuracy_level == 1) fact = fact * 50
   if (diis%accuracy_level == 2) fact = fact * 25
   if (diis%accuracy_level == 3) fact = fact * six
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
 nok = count(checks /= 0)
 if (diis%accuracy_level == 1) ans = nok >= half * ndat
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

 integer :: iter, idat, iband
 real(dp) :: deltae, deold, dedold, absdiff
 character(len=500) :: msg

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

subroutine getghc_eigresid(gs_hamk, npw, nspinor, ndat, cg, ghc, gsc, mpi_enreg, prtvol, &
                           eig, resid, enlx, residvecs, gvnlxc, normalize)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: npw, nspinor, ndat, prtvol
 real(dp),intent(inout) :: cg(2, npw*nspinor*ndat)
 real(dp),intent(out) :: ghc(2,npw*nspinor*ndat), gsc(2,npw*nspinor*ndat*gs_hamk%usepaw)
 type(mpi_type),intent(inout) :: mpi_enreg
 real(dp),intent(out) :: eig(ndat), resid(ndat), enlx(ndat)
 real(dp),intent(out) :: residvecs(2, npw*nspinor*ndat)
 real(dp),intent(out) :: gvnlxc(2, npw*nspinor*ndat)
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
 npwsp = npw * nspinor
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

 ! PAW normalization must be done here.
 if (usepaw == 1 .and. normalize_) call cgpaw_normalize(npwsp, ndat, cg, gsc, istwf_k, me_g0, comm)

 ! Compute new approximated eigenvalues, residue vectors and residuals.
 call cg_eigens(usepaw, istwf_k, npwsp, ndat, cg, ghc, gsc, eig, me_g0, comm)
 call cg_residvecs(usepaw, npwsp, ndat, eig, cg, ghc, gsc, residvecs)
 call cg_norm2g(istwf_k, npwsp, ndat, residvecs, resid, me_g0, comm)

 if (usepaw == 0) then
   ! Evaluate new enlx for NC.
   call cg_zdotg(istwf_k, npwsp, ndat, option1, cg, gvnlxc, dots, me_g0, comm)
   enlx = dots(1,:)
 end if

 !write(std_out, *)" In getghc_eigresid:"
 !do idat=1,ndat
 !  write(std_out, *)" eig:", eig(idat), "idat:", idat
 !  write(std_out, *)" resid:", resid(idat), "idat:", idat
 !  write(std_out, *)" enlx:", enlx(idat), "idat:", idat
 !end do
 !write(std_out, *)""

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
 ABI_MALLOC(diis%hist_enlx, (0:max_niter+2, bsize))
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

 integer,parameter :: master = 0
 integer :: cplex, ierr, nprocs, my_rank, idat, ii !, ibk, iek
 real(dp) :: noise, cpu, wall, gflops
 !integer :: failed(ndat)
 real(dp),allocatable :: diis_eig(:), wmat1(:,:,:), wmat2(:,:,:), wvec(:,:,:), alphas(:,:)
 character(len=500) :: msg

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 cplex = diis%cplex
 !failed = 0

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

     ! TODO: DGEMV
     ABI_MALLOC(alphas, (1, 0:iter))
     alphas(1,:) = wvec(1,:,idat)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_phi(:,:,:,idat), 2*npwsp, alphas, 1, zero, phi_bk(:,:,idat), 1)
     call dgemv("N", 2*npwsp, iter, one, diis%chain_resv(:,:,:,idat), 2*npwsp, alphas, 1, zero, residv_bk(:,:,idat), 1)

     !call dgemv(character 	TRANS, integer 	M, integer 	N, double precision 	ALPHA,
     !           double precision, dimension(lda,*) 	A, integer 	LDA,
     !           double precision, dimension(*) 	X, integer 	INCX,
     !           double precision 	BETA, double precision, dimension(*) 	Y,integer 	INCY)

     ABI_FREE(alphas)
   end if
 end do

 ABI_FREE(wvec)

 if (timeit) call cwtime_report(" build_hij", cpu, wall, gflops)

end subroutine rmm_diis_update_block
!!***

!!****f* m_rmm_diis/rmm_diis_eval_mats
!! NAME
!!  rmm_diis_eval_mats
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

subroutine rmm_diis_eval_mats(diis, iter, ndat, me_g0, comm)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat, me_g0, comm

 integer :: ii, ierr, idat, nprocs, option
 real(dp) :: dotr, doti, cpu, wall, gflops

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

subroutine rmm_diis_rollback(diis, npwsp, ndat, phi_bk, gsc_bk, eig, resid, enlx, comm)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: npwsp, comm, ndat
 real(dp),intent(inout) :: phi_bk(2, npwsp, ndat), gsc_bk(2, npwsp, ndat*diis%usepaw)
 real(dp),intent(inout) :: eig(ndat), resid(ndat), enlx(ndat)

 integer,parameter :: master = 0
 integer :: nprocs, my_rank, idat, iter, ierr, ilast, take_iter(ndat)
 real(dp) :: cpu, wall, gflops

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 take_iter = -1
 ilast = diis%last_iter

 if (timeit) call cwtime(cpu, wall, gflops, "start")

 if (my_rank == master) then
   do idat=1,ndat
     iter = imin_loc(diis%hist_resid(0:ilast, idat)) ! back = .True.
     iter = iter - 1
     ! Move stuff only if it's really worth it.
     if (iter /= ilast .and. &
         diis%hist_resid(iter, idat) < half * diis%hist_resid(ilast, idat)) take_iter(idat) = iter
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
     call cg_zcopy(npwsp, diis%chain_phi(:,:,iter,idat), phi_bk(:,:,idat))
     if (diis%usepaw == 1) call cg_zcopy(npwsp, diis%chain_sphi(:,:,iter,idat), gsc_bk(:,:,idat))
     call diis%stats%increment("rollback", 1)
   end do
 end if

 if (timeit) call cwtime_report(" diis_rollback", cpu, wall, gflops)

end subroutine rmm_diis_rollback
!!***

!!****f* m_cgtools/cg_eigens
!! NAME
!!  cg_eigens
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

subroutine cg_eigens(usepaw, istwf_k, npwsp, ndat, cg, ghc, gsc, eig, me_g0, comm)

 integer,intent(in) :: usepaw, istwf_k, npwsp, ndat, me_g0, comm
 real(dp),intent(in) :: ghc(2*npwsp, ndat), cg(2*npwsp, ndat), gsc(2*npwsp, ndat*usepaw)
 real(dp),intent(out) :: eig(ndat)

 integer,parameter :: option1 = 1
 integer :: idat, ierr
 real(dp) :: dotr(ndat), doti

!$OMP PARALLEL DO
 do idat=1,ndat
   call dotprod_g(eig(idat), doti, istwf_k, npwsp, option1, ghc(:,idat), cg(:,idat), me_g0, xmpi_comm_self)
   if (usepaw == 1) then
     call dotprod_g(dotr(idat), doti, istwf_k, npwsp, option1, gsc(:,idat), cg(:,idat), me_g0, xmpi_comm_self)
   end if
 end do

 if (xmpi_comm_size(comm) > 1) then
   call xmpi_sum(eig, comm, ierr)
   if (usepaw == 1) call xmpi_sum(dotr, comm, ierr)
 end if

 if (usepaw == 1) eig(:) = eig(:) / dotr(:)

end subroutine cg_eigens
!!***

!!****f* m_cgtools/cg_residvecs
!! NAME
!!  cg_residvecs
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

subroutine cg_residvecs(usepaw, npwsp, ndat, eig, cg, ghc, gsc, residvecs)

 integer,intent(in) :: usepaw, npwsp, ndat
 real(dp),intent(in) :: eig(ndat)
 real(dp),intent(in) :: ghc(2*npwsp, ndat), cg(2*npwsp, ndat), gsc(2*npwsp, ndat*usepaw)
 real(dp),intent(out) :: residvecs(2*npwsp, ndat)

 integer :: idat

 if (usepaw == 1) then
   ! (H - e) |psi>
!$OMP PARALLEL DO
   do idat=1,ndat
     residvecs(:,idat) = ghc(:,idat) - eig(idat) * gsc(:,idat)
   end do
 else
   ! (H - eS) |psi>
!$OMP PARALLEL DO
   do idat=1,ndat
     residvecs(:,idat) = ghc(:,idat) - eig(idat) * cg(:,idat)
   end do
 end if

end subroutine cg_residvecs
!!***

!!****f* m_cgtools/cg_norm2g
!! NAME
!!  cg_norm2g
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

subroutine cg_norm2g(istwf_k, npwsp, ndat, cg, norms, me_g0, comm)

 integer,intent(in) :: istwf_k, npwsp, ndat, me_g0, comm
 real(dp),intent(in) :: cg(2*npwsp, ndat)
 real(dp),intent(out) :: norms(ndat)

 integer :: idat, ierr

!$OMP PARALLEL DO
 do idat=1,ndat
   call sqnorm_g(norms(idat), istwf_k, npwsp, cg(:,idat), me_g0, xmpi_comm_self)
 end do

 if (xmpi_comm_size(comm) > 1) call xmpi_sum(norms, comm, ierr)

end subroutine cg_norm2g
!!***

!!****f* m_cgtools/cg_zdotg
!! NAME
!!  cg_zdotg
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

subroutine cg_zdotg(istwf_k, npwsp, ndat, option, cg1, cg2, dots, me_g0, comm)

 integer,intent(in) :: istwf_k, npwsp, ndat, option, me_g0, comm
 real(dp),intent(in) :: cg1(2*npwsp,ndat), cg2(2*npwsp,ndat)
 real(dp),intent(out) :: dots(2,ndat)

 integer :: idat, ierr
 real(dp) :: dotr, doti, re_dots(ndat)

!$OMP PARALLEL DO PRIVATE(dotr, doti)
 do idat=1,ndat
   call dotprod_g(dotr, doti, istwf_k, npwsp, option, cg1(:,idat), cg2(:,idat), me_g0, xmpi_comm_self)
   if (istwf_k == 2) then
     re_dots(idat) = dotr
   else
     dots(:, idat) = [dotr, doti]
   end if
 end do

 if (xmpi_comm_size(comm) > 1) then
   if (istwf_k == 2) then
     call xmpi_sum(re_dots, comm, ierr)
   else
     call xmpi_sum(dots, comm, ierr)
   end if
 end if

 if (istwf_k == 2) then
   do idat=1,ndat
     dots(1,idat) = re_dots(idat)
     dots(2,idat) = zero
   end do
 end if

end subroutine cg_zdotg
!!***

!!****f* m_cgtools/cg_precon_many
!! NAME
!!  cg_precon_many
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

subroutine cg_precon_many(istwf_k, npw, nspinor, ndat, cg, optekin, kinpw, vect, me_g0, comm)

 integer,intent(in) :: istwf_k, npw, nspinor, optekin, ndat, me_g0, comm
 real(dp),intent(in) :: cg(2*npw*nspinor,ndat), kinpw(npw)
 real(dp),intent(inout) :: vect(2*npw*nspinor,ndat)

 integer :: idat
 real(dp),allocatable :: pcon(:)

 ! TODO: Optimized version for MPI with ndat > 1
 !return
 ABI_MALLOC(pcon, (npw))
 do idat=1,ndat
   call cg_precon(cg(:,idat), zero, istwf_k, kinpw, npw, nspinor, me_g0, optekin, pcon, vect(:,idat), comm)
 end do
 ABI_FREE(pcon)

 !call cg_kinene(istwf_k, npw, nspinor, ndat, cg, me_g0, comm)
 !call cg_zprecon_block(cg,eval,blocksize,iterationnumber,kinpw, npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize,comm)

end subroutine cg_precon_many
!!***

!----------------------------------------------------------------------

!!****f* m_cgtools/cg_zaxpy_many_areal
!! NAME
!!  cg_zaxpy_many_areal
!!
!! FUNCTION
!!  Computes y = alpha*x + y
!!
!! INPUTS
!!  n = Specifies the number of elements in vectors x and y.
!!  ndat
!!  alpha(ndat) = Specifies the scalar alpha.
!!  x = Array
!!
!! SIDE EFFECTS
!!  y = Array. In output, y contains the updated vector.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cg_zaxpy_many_areal(npwsp, ndat, alphas, x, y)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwsp, ndat
 real(dp),intent(in) :: alphas(ndat)
!arrays
 real(dp),intent(in) :: x(2*npwsp, ndat)
 real(dp),intent(inout) :: y(2*npwsp, ndat)

!Local variables-------------------------------
 integer :: idat

! *************************************************************************

!$OMP PARALLEL DO
 do idat=1,ndat
   call daxpy(2*npwsp, alphas(idat), x(1,idat), 1, y(1,idat), 1)
 end do

end subroutine cg_zaxpy_many_areal
!!***


!!****f* m_numeric_tools/pack_matrix
!! NAME
!! pack_matrix
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

end module m_rmm_diis
!!***
