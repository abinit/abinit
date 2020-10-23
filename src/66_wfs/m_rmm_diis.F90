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

 use defs_abitypes,   only : mpi_type
 use m_fstrings,      only : sjoin, itoa, ftoa
 use m_time,          only : timab, cwtime, cwtime_report
 use m_numeric_tools, only : pack_matrix
 use m_hide_lapack,   only : xhegv_cplex, xhesv_cplex
 use m_pair_list,     only : pair_list
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 !use m_fock,         only : fock_set_ieigen, fock_set_getghc_call

 implicit none

 private
!!***

 public :: rmm_diis
!!***

 type,private :: rmm_diis_t
   integer :: usepaw
   integer :: istwf_k
   integer :: cplex
   integer :: bsize
   integer :: max_niter
   integer :: npwsp
   type(pair_list) :: stats

   integer :: last_iter
   ! (bsize)

   real(dp),allocatable :: hist_ene(:,:)
   real(dp),allocatable :: hist_resid(:,:)
   character(len=6),allocatable :: step_type(:,:)
   ! (0:max_niter+1, bsize)
   ! 0 is the initial step, then DIIS iterations (whose number may depend on the block)
   ! followed by an optional trial step.

   real(dp),allocatable :: resmat(:,:,:,:)
   real(dp),allocatable :: smat(:,:,:,:)
   ! (2, 0:max_niter, 0:max_niter, bsize))

   real(dp),allocatable :: chain_phi(:,:,:,:)
   real(dp),allocatable :: chain_sphi(:,:,:,:)
   real(dp),allocatable :: chain_resv(:,:,:,:)
   ! (2, npwsp, 0:max_niter, bsize))

 contains
   procedure :: free => rmm_diis_free
   procedure :: solve => rmm_diis_solve
   procedure :: eval_mats => rmm_diis_eval_mats
   procedure :: exit_iter => rmm_diis_exit_iter
   procedure :: print_block => rmm_diis_print_block

 end type rmm_diis_t

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
!!
!! PARENTS
!!      m_vtowfk
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine rmm_diis(istep, ikpt, isppol, cg, dtset, eig, occ, enlx, gs_hamk, gsc_all, &
                    mpi_enreg, nband, npw, nspinor, resid)

!Arguments ------------------------------------
 integer,intent(in) :: istep, ikpt, isppol, nband, npw, nspinor
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(in) :: mpi_enreg
 real(dp),intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband)
 real(dp),intent(inout) :: resid(nband)
 real(dp),intent(in) :: occ(nband)
 real(dp),intent(out) :: eig(nband)

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2
 integer,parameter :: tim_getghc = 0, level = 432, use_subovl0 = 0
 integer :: ig, ig0, ib, ierr, prtvol, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop, idat
 integer :: iband, cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, optekin, usepaw, iter, max_niter, max_niter_block
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft, nbocc, jj, kk, it, ld1, ld2, ibk, iek !ii,
 logical :: end_with_trial_step, new_lambda
 real(dp),parameter :: rdummy = zero
 real(dp) :: accuracy_ene, cpu, wall, gflops, dotr, doti
 !character(len=500) :: msg
 character(len=6) :: tag
!arrays
 real(dp) :: tsec(2)
 real(dp),target :: fake_gsc_bk(0,0)
 real(dp) :: subovl(use_subovl0)
 real(dp),allocatable :: evec(:,:), subham(:), h_ij(:,:,:)
 real(dp),allocatable :: ghc_bk(:,:), gsc_bk(:,:), gvnlxc_bk(:,:), lambda_bk(:), dots_bk(:,:)
 real(dp),allocatable :: residv_bk(:,:), kres_bk(:,:), phi_bk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: ptr_gsc_bk(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(rmm_diis_t) :: diis
 type(yamldoc_t) :: ydoc

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; prtvol = dtset%prtvol
 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
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
 !bsize = 4
 bsize = 8
 if (dtset%userib /= 0) bsize = abs(dtset%userib)
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 ABI_MALLOC(ghc_bk, (2, npwsp*bsize))
 ABI_MALLOC(gvnlxc_bk, (2, npwsp*bsize))
 ABI_CALLOC(h_ij, (2, nband, nband))

 call timab(1633, 1, tsec) !"rmm_diis:build_hij
 call cwtime(cpu, wall, gflops, "start")

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize
   ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize
   ib_stop = min(iblock * bsize, nband)

   if (usepaw == 1) ptr_gsc_bk => gsc_all(:,igs:ige)
   call getghc(cpopt, cg(:,igs:ige), cprj_dum, ghc_bk, ptr_gsc_bk, gs_hamk, gvnlxc_bk, &
               rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

   ! Compute <i|H|j> for i=1,nband and all j in block
   call cg_zgemm("C", "N", npwsp, nband, ndat, cg, ghc_bk, h_ij(:,:,ib_start))

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
       h_ij(2,:,iband) = zero
     end do
   end if

 end do ! iblock

 call cwtime_report(" build_hij", cpu, wall, gflops)

 ! Pack <i|H|j> to prepare call to subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 do iband=1,nband
   h_ij(2,iband,iband) = zero ! Force diagonal elements to be real
 end do
 call pack_matrix(h_ij, subham, nband, 2)
 ABI_FREE(h_ij)
 call xmpi_sum(subham, comm_spinorfft, ierr)
 call timab(1633, 2, tsec) !"rmm_diis:build_hij

 ! =================
 ! Subspace rotation
 ! =================
 call timab(585, 1, tsec) !"vtowfk(subdiago)"
 ABI_MALLOC(evec, (2*nband, nband))
 mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 call subdiago(cg, eig, evec, gsc_all, 0, 0, istwf_k, mcg, mgsc, nband, npw, nspinor, dtset%paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)
 call timab(585, 2, tsec)
 call cwtime_report(" subdiago", cpu, wall, gflops)

 ABI_FREE(subham)
 ABI_FREE(evec)

 ! =================
 ! Prepare DIIS loop
 ! =================
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1
 !write(std_out,*)"optekin:", optekin
 ! nline - 1 DIIS steps, then end with trial step (optional, see below)
 max_niter = dtset%nline - 1

 ! Define accuracy_ene for SCF
 nbocc = count(occ > zero)
 accuracy_ene = zero
 if (dtset%iscf > 0) then
   if (dtset%toldfe /= 0) then
     accuracy_ene = dtset%toldfe / nbocc / four
   else
     ! We are not using toldfe to stop the SCF cycle
     ! so we are forced to hardcode a tolerance for the absolute diff in the KS eigenvalue.
     accuracy_ene = tol6 / nbocc / four
   end if
 end if

 diis = rmm_diis_new(usepaw, istwf_k, npwsp, max_niter, bsize)
 call wrtout(std_out, sjoin(" Using Max", itoa(max_niter), "RMM-DIIS iterations + final trial step."))
 call wrtout(std_out, sjoin(" Block size:", itoa(bsize), " accuracy_ene: ", ftoa(accuracy_ene)))
 call timab(1634, 1, tsec) ! "rmm_diis:band_opt"

 ABI_MALLOC(lambda_bk, (bsize))
 ABI_MALLOC(dots_bk, (2, bsize))
 ABI_MALLOC(phi_bk, (2, npwsp*bsize))
 ABI_MALLOC(gsc_bk, (2, npwsp*bsize*usepaw))
 ABI_MALLOC(residv_bk, (2, npwsp*bsize))
 ABI_MALLOC(kres_bk, (2, npwsp*bsize))

 ! We loop over nblocks, each block has ndat states.
 ! Convergence behaviour may depend on bsize as branches are taken according to
 ! the status of all bands in the block.
 ! Using gemm_nonlop leads to a significant speedup when applying Vnl.

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize
   ige = min(iblock * npwsp * bsize, npwsp * nband)
   ndat = (ige - igs + 1) / npwsp
   ib_start = 1 + (iblock - 1) * bsize
   ib_stop = min(iblock * bsize, nband)

   ! Reduce number of niter iterations for empty states.
   ! TODO: Don't reduce niter if MD
   !if dtset%
   !if dtset%occopt == 2
   max_niter_block = max_niter
   if (all(occ(ib_start:ib_stop) < tol3)) max_niter_block = max(1 + max_niter / 2, 2)

   ! Compute H |phi_0> using cg from subdiago. Blocked call.
   ! Alternatively, one can compute the residuals before the subspace rotation and then rotate
   ! to save one call to getghc_eigresid.

   ! Extract bands in the block. phi_bk = cg(:,igs:ige)
   call cg_zcopy(npwsp * ndat, cg(:,igs), phi_bk)
   if (usepaw == 1) ptr_gsc_bk => gsc_all(:,igs:ige)

   call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                        eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, normalize=.False.)

   ! Save <R0|R0> and <phi_0|S|phi_0>, |phi_0>, |S phi_0>. Assume input cg are already S-normalized.
   do idat=1,ndat
     diis%hist_ene(0, idat) = eig(ib_start+idat-1)
     diis%hist_resid(0, idat) = resid(ib_start+idat-1)
     if (diis%cplex == 2) then
       diis%resmat(:, 0, 0, idat) = [resid(ib_start+idat-1), zero]
       diis%smat(:, 0, 0, idat) = [one, zero]
     else
       diis%resmat(:, 0, 0, idat) = resid(ib_start+idat-1)
       diis%smat(:, 0, 0, idat) = one
     end if
     diis%step_type(0, idat) = "SDIAG"
     jj = 1 + (idat - 1) * npwsp
     call cg_zcopy(npwsp, residv_bk(1,jj), diis%chain_resv(:,:,0,idat))
     call cg_zcopy(npwsp, phi_bk(1,jj), diis%chain_phi(:,:,0,idat))
     if (usepaw == 1) call cg_zcopy(npwsp, ptr_gsc_bk(1:,jj), diis%chain_sphi(:,:,0,idat))
   end do

   ! Precondition |R_0>, output in kres_bk = |K R_0>
   call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
   call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, gs_hamk%kinpw_k, kres_bk, me_g0, comm_spinorfft)

   ! Compute H |K R_0>
   call getghc(cpopt, kres_bk, cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
               rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

   ! Compute residuals: (H - e_0 S) |K R_0>
   call cg_residvecs(usepaw, npwsp, ndat, eig(ib_start:), kres_bk, ghc_bk, gsc_bk, residv_bk)

   ! Line minimization with preconditioned steepest descent:
   !
   !    |phi_1> = |phi_0> + lambda |K R_0>
   !
   ! where lambda minimizes the norm of the residual (not the Rayleigh quotient as in Kresse's paper).
   !
   !    lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
   !
   call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_spinorfft)

#if 0
   ld1 = npwsp * (max_niter + 1); ld2 = npwsp
   call cg_zdotg_lds(istwf_k, npwsp, ndat, option1, ld1, diis%chain_resv, ld2, residv_bk, dots_bk, me_g0, comm_spinorfft)
   do idat=1,ndat
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     lambda_bk(idat) = -dots_bk(1,idat) / lambda_bk(idat)
     phi_bk(:,jj:kk) = diis%chain_phi(:,:,0,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
   end do

#else
   !lamdas(ii) = -dots_bk(1, ii) / lambda_bk(ii)
   do idat=1,ndat
     jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
     call dotprod_g(dotr, doti, istwf_k, npwsp, option1, &
                    diis%chain_resv(:,:,0,idat), residv_bk(:,jj), me_g0, comm_spinorfft)

     lambda_bk(idat) = -dotr / lambda_bk(idat)
     phi_bk(:,jj:kk) = diis%chain_phi(:,:,0,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
     !write(std_out, *)"lambda, dotr, rval ", lambda_bk(1), dotr, rval
   end do
#endif

   iter_loop: do iter=1,max_niter_block

     if (iter > 1) then
       ! Solve DIIS equations to get phi_bk for iter > 1
       do idat=1,ndat
         ibk = 1 + (idat - 1) * npwsp; iek = idat * npwsp
         call diis%solve(iter, npwsp, idat, phi_bk(:,ibk), residv_bk(:,ibk), comm_spinorfft)
       end do

       ! Precondition residual, output in kres_bk.
       call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
       call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, gs_hamk%kinpw_k, &
                          kres_bk, me_g0, comm_spinorfft)

       ! Compute phi_bk with the lambda(ndat) obtained at iteration #1
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)
     end if

     ! Compute H |phi_now> and evaluate new enlx for NC.
     call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, normalize=.True.)

     ! Store new residual.
     do idat=1,ndat
       diis%hist_ene(iter, idat) = eig(ib_start+idat-1)
       diis%hist_resid(iter, idat) = resid(ib_start+idat-1)
       diis%step_type(iter, idat) = "DIIS"
       ibk = 1 + (idat - 1) * npwsp; iek = idat * npwsp
       diis%chain_phi(:,:,iter, idat) = phi_bk(:,ibk:iek)
       diis%chain_resv(:, :, iter, idat) = residv_bk(:,ibk:iek)
       if (usepaw == 1) diis%chain_sphi(:,:,iter,idat) = ptr_gsc_bk(:,ibk:iek)
     end do

     ! ============== CHECK FOR CONVERGENCE ========================
     if (diis%exit_iter(iter, ndat, max_niter_block, accuracy_ene, dtset)) exit iter_loop

     ! Compute <R_i|R_j> and <i|S|j> for j=iter
     if (iter /= max_niter_block) call diis%eval_mats(iter, ndat, me_g0, comm_spinorfft)

   end do iter_loop

   ! End with trial step but only if we performed all the iterations (no exit from diis%exit)
   ! Since we operate on blocks of bands, all the states in the block will receive the same treatment.
   ! This means that one can observe a (hopefully) slightly different convergence behaviour depending on bsize.
   !
   end_with_trial_step = iter == max_niter_block + 1
   !end_with_trial_step = .False.
   !end_with_trial_step = .True.

   if (end_with_trial_step) then

     call cg_zcopy(npwsp * ndat, residv_bk, kres_bk)
     call cg_precon_many(istwf_k, npw, nspinor, ndat, phi_bk, optekin, gs_hamk%kinpw_k, kres_bk, me_g0, comm_spinorfft)

     ! Recomputing the preconditioned-direction and the new lambda that minimizes the residual
     ! is more expensive but more accurate. Do it only for the occupied states plus a buffer
     ! if these bands still far from convergence. In all the other cases, we reuse the previous lambda.

     new_lambda = .False.
     if (ib_stop <= nbocc + 4 .and. any(resid(ib_start:ib_stop) > tol10)) new_lambda = .True.
     !new_lambda = .False.
     !new_lambda = .True.

     if (new_lambda) then
       tag = "NEWLAM"
       !write(std_out, *)" Performing last trial step with new computation of lambda"
       ! Final preconditioned steepest descent with new lambda.
       ! Compute H |K R_0>
       call getghc(cpopt, kres_bk, cprj_dum, ghc_bk, gsc_bk, gs_hamk, gvnlxc_bk, &
                   rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

       ! Compute residual: (H - e_0 S) |K R_0>
       call cg_residvecs(usepaw, npwsp, ndat, eig(ib_start:), kres_bk, ghc_bk, gsc_bk, residv_bk)

       !ld1 = npwsp * (max_niter + 1); ld2 = npwsp
       !call cg_zdotg_lds(istwf_k, npwsp, ndat, option1, ld1, cg1, ld2, residv_bk, dots_bk, me_g0, comm)

       call cg_norm2g(istwf_k, npwsp, ndat, residv_bk, lambda_bk, me_g0, comm_spinorfft)

       it = diis%last_iter
       do idat=1,ndat
         jj = 1 + (idat - 1) * npwsp; kk = idat * npwsp
         call dotprod_g(dotr, doti, istwf_k, npwsp, option1, &
                        diis%chain_resv(:,:,it,idat), residv_bk(:,jj), me_g0, comm_spinorfft)

         lambda_bk(idat) = -dotr / lambda_bk(idat)
         phi_bk(:,jj:kk) = diis%chain_phi(:,:,it,idat) + lambda_bk(idat) * kres_bk(:,jj:kk)
       end do

     else
       ! Reuse previous values of lambda.
       tag = "FIXLAM"
       !phi_bk = phi_bk + 0.1 * kres_bk
       call cg_zaxpy_many_areal(npwsp, ndat, lambda_bk, kres_bk, phi_bk)
     endif

     ! Finally recompute eig, resid, enlx and ptr_gsc_bk if PAW.
     call getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_bk, ghc_bk, ptr_gsc_bk, mpi_enreg, prtvol, &
                          eig(ib_start:), resid(ib_start:), enlx(ib_start:), residv_bk, normalize=.True.)

     ! Push the last results just for printing purposes.
     it = diis%last_iter + 1
     do idat=1,ndat
       diis%hist_ene(it, idat) = eig(ib_start+idat-1)
       diis%hist_resid(it, idat) = resid(ib_start+idat-1)
       diis%step_type(it, idat) = tag
     end do

   end if ! end_with_trial_step

   ! Update wavefunction block. cg(:,igs:ige) = phi_bk
   call cg_zcopy(npwsp * ndat, phi_bk, cg(:,igs))

   if (prtvol == -level) call diis%print_block(ib_start, ndat, istep, ikpt, isppol)
 end do ! iblock

 call timab(1634, 2, tsec) !"rmm_diis:band_opt"
 call cwtime_report(" rmm_diis:band_opt", cpu, wall, gflops)

 !call yaml_write_and_free_dict('RMM-DIIS', dict, std_out)
 ydoc = yamldoc_open('RMM-DIIS', with_iter_state=.False.)
 call ydoc%add_dict("stats", diis%stats)
 call ydoc%write_and_free(std_out)
 call diis%stats%free()

 ! Final cleanup
 ABI_FREE(lambda_bk)
 ABI_FREE(dots_bk)
 ABI_FREE(ghc_bk)
 ABI_FREE(gsc_bk)
 ABI_FREE(residv_bk)
 ABI_FREE(kres_bk)
 ABI_FREE(phi_bk)
 ABI_FREE(gvnlxc_bk)

 call diis%free()

end subroutine rmm_diis
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

logical function rmm_diis_exit_iter(diis, iter, ndat, niter_block, accuracy_ene, dtset) result(ans)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat, niter_block
 real(dp),intent(in) :: accuracy_ene
 type(dataset_type),intent(in) :: dtset

 integer :: idat, checks(ndat)
 real(dp) :: resid, deltae, deold
 character(len=50) :: msg_list(ndat)

 diis%last_iter = iter
 checks = 0

 do idat=1,ndat
   resid = diis%hist_resid(iter, idat)
   deold = diis%hist_ene(1, idat) - diis%hist_ene(0, idat)
   deltae = diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat)

   ! Relative criterion on eig diff
   ! Abinit default in the CG part is 0.005 (0.3 V).
   ! Here we enlarge it a bit
   if (abs(deltae) < 6 * dtset%tolrde * abs(deold) .and. iter /= niter_block) then
   !if (abs(deltae) < 0.3 * abs(deold) .and. iter /= niter_block) then
     checks(idat) = 1; msg_list(idat) = "deltae < 6 * tolrde * deold"; cycle
   end if

   if (dtset%iscf < 0) then
     ! This is the only condition available for NSCF run.
     if (resid < dtset%tolwfr) then
       checks(idat) = 1; msg_list(idat) = 'resid < tolwfr'; cycle
     end if

   else
     ! Condition available for SCF run.
     if (resid < dtset%tolwfr) then
       checks(idat) = 1; msg_list(idat) = 'resid < tolwfr'; cycle
     end if
     if (sqrt(abs(resid)) < accuracy_ene) then
       ! Absolute criterion on eig diff
       checks(idat) = 1; msg_list(idat) = 'resid < accuracy_ene'; cycle
     end if
   end if
 end do ! idat

 ans = all(checks /= 0)
 if (ans) then
   do idat=1,ndat
     call diis%stats%increment(msg_list(idat), 1)
   end do
   !if (prtvol == -level) then
   !  write(msg, '(a,i4,a,i2,a,es12.4,a)' )&
   !   ' band: ',iband,' converged after: ',iter,' iterations with resid: ',resid(iband), ch10
   !  call wrtout(std_out, sjoin("<END RMM-DIIS, msg='", msg, "'>"))
   !end if
 end if

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
   sjoin("<BEGIN RMM-DIIS-BLOCK, istep:", itoa(istep), ", ikpt:", itoa(ikpt), ", spin: ", itoa(isppol), ">"))

 do idat=1,ndat
   write(msg,'(1a, 2(a5), 4(a14), 1x, a6)') &
     "#", 'iter', "band", "eigen_eV", "eigde_meV", "de/dold", "resid", "type"; call wrtout(std_out, msg)

   iband = ib_start + idat - 1
   deold = diis%hist_ene(1, idat) - diis%hist_ene(0, idat)
   do iter=0,diis%last_iter

     deltae = diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat)
     dedold = zero; absdiff = zero
     if (iter > 0) then
       dedold = deltae / deold
       absdiff = (diis%hist_ene(iter, idat) - diis%hist_ene(iter-1, idat))
     end if

     write(msg,"(1x, 2(i5), 2(f14.6), 2(es14.6), 1x, a6)") &
       iter, iband, diis%hist_ene(iter, idat) * Ha_eV, absdiff * Ha_meV, dedold, &
       diis%hist_resid(iter, idat), diis%step_type(iter, idat); call wrtout(std_out, msg)
   end do
 end do
 call wrtout(std_out, "<END RMM-DIIS-BLOCK>")

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
                           eig, resid, enlx, residvecs, normalize)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: npw, nspinor, ndat, prtvol
 real(dp),intent(inout) :: cg(2, npw*nspinor*ndat)
 real(dp),intent(out) :: ghc(2,npw*nspinor*ndat), gsc(2,npw*nspinor*ndat*gs_hamk%usepaw)
 type(mpi_type),intent(in) :: mpi_enreg
 real(dp),intent(out) :: eig(ndat), resid(ndat), enlx(ndat)
 real(dp),intent(out) :: residvecs(2, npw*nspinor*ndat)
 logical,optional,intent(in) :: normalize

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: istwf_k, usepaw, cpopt, sij_opt, npwsp, me_g0, comm_spinorfft, comm_fft
 real(dp),parameter :: rdummy = zero
 logical :: normalize_
!arrays
 real(dp),allocatable :: gvnlxc(:, :)
 real(dp) :: dots(2, ndat)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k
 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft
 normalize_ = .True.; if (present(normalize)) normalize_ = normalize
 npwsp = npw * nspinor

 cpopt = -1; sij_opt = 0
 if (usepaw == 1) then
   sij_opt = 1 ! matrix elements <G|S|C> have to be computed in gsc in addition to ghc
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 end if

 ! NC normalization.
 if (usepaw == 0 .and. normalize_) call cgnc_normalize(npwsp, ndat, cg, istwf_k, me_g0, comm_spinorfft)

 ABI_MALLOC(gvnlxc, (2, npwsp*ndat))

 ! Compute H |cg>
 !call fock_set_ieigen(gs_hamk%fockcommon, iband)
 call getghc(cpopt, cg, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
             rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

 ! PAW normalization must be done here.
 if (usepaw == 1 .and. normalize_) call cgpaw_normalize(npwsp, ndat, cg, gsc, istwf_k, me_g0, comm_spinorfft)

 ! Compute new approximated eigenvalues, resid-vectors and residuals.
 call cg_eigens(usepaw, istwf_k, npwsp, ndat, cg, ghc, gsc, eig, me_g0, comm_spinorfft)
 call cg_residvecs(usepaw, npwsp, ndat, eig, cg, ghc, gsc, residvecs)
 call cg_norm2g(istwf_k, npwsp, ndat, residvecs, resid, me_g0, comm_spinorfft)

 if (usepaw == 0) then
   ! Evaluate new enlx for NC.
   call cg_zdotg(istwf_k, npwsp, ndat, option1, cg, gvnlxc, dots, me_g0, comm_spinorfft)
   enlx = dots(1,:)
 end if

 ABI_FREE(gvnlxc)

end subroutine getghc_eigresid
!!!***

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

type(rmm_diis_t) function rmm_diis_new(usepaw, istwf_k, npwsp, max_niter, bsize) result(diis)

 integer,intent(in) :: usepaw, istwf_k, npwsp, max_niter, bsize

 diis%usepaw = usepaw
 diis%istwf_k = istwf_k
 diis%cplex = 2; if (istwf_k == 2) diis%cplex = 1
 diis%npwsp = npwsp
 diis%max_niter = max_niter
 diis%bsize = bsize

 ABI_MALLOC(diis%hist_ene, (0:max_niter+1, bsize))
 ABI_MALLOC(diis%hist_resid, (0:max_niter+1, bsize))
 ABI_MALLOC(diis%step_type, (0:max_niter+1, bsize))

 ABI_MALLOC(diis%chain_phi, (2, npwsp, 0:max_niter, bsize))
 ABI_MALLOC(diis%chain_sphi, (2, npwsp*usepaw, 0:max_niter, bsize))
 ABI_MALLOC(diis%chain_resv, (2, npwsp, 0:max_niter, bsize))
 ABI_CALLOC(diis%resmat, (diis%cplex, 0:max_niter, 0:max_niter, bsize)) ! <R_i|R_j>
 ABI_CALLOC(diis%smat, (diis%cplex, 0:max_niter, 0:max_niter, bsize))   ! <i|S|j>

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
 ABI_SFREE(diis%step_type)
 ABI_SFREE(diis%chain_phi)
 ABI_SFREE(diis%chain_sphi)
 ABI_SFREE(diis%chain_resv)
 ABI_SFREE(diis%resmat)
 ABI_SFREE(diis%smat)

 call diis%stats%free()

end subroutine rmm_diis_free
!!***

!!****f* m_rmm_diis/rmm_diis_solve
!! NAME
!!  rmm_diis_solve
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

subroutine rmm_diis_solve(diis, iter, npwsp, idat, new_psi, new_residvec, comm_spinorfft)

 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iter, npwsp, comm_spinorfft, idat
 real(dp),intent(out) :: new_psi(2, npwsp), new_residvec(2, npwsp)

 integer,parameter :: master = 0
 integer :: cplex, ierr, nprocs, my_rank
 real(dp),allocatable :: diis_eig(:)
 real(dp),allocatable :: wmat1(:,:,:), wmat2(:,:,:), wvec(:,:), alphas(:,:)
 character(len=500) :: msg
 logical,parameter :: try_to_solve_eigproblem = .False.

 cplex = diis%cplex
 my_rank = xmpi_comm_rank(comm_spinorfft); nprocs = xmpi_comm_size(comm_spinorfft)

 !do idat=1,ndat

 if (try_to_solve_eigproblem) then
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
   !ABI_CHECK(ierr == 0, "xhegv returned ierr != 0")
   if (ierr /= 0) then
     !call wrtout(std_out, sjoin("xhegv failed with:", msg, ch10, "at iter: ", itoa(iter), "exit iter_loop!"))
     ABI_FREE(diis_eig)
     ABI_FREE(wmat1)
     ABI_FREE(wmat2)
     goto 10
   end if

   ! Take linear combination of chain_phi and chain_resv.
   call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), wmat1(:,:,0), new_psi)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), wmat1(:,:,0), new_residvec)

   ABI_FREE(wmat1)
   ABI_FREE(wmat2)
   ABI_FREE(diis_eig)
   return
 end if

10 continue

 ! Solve system of linear equations.
 ! Only master works so that we are sure we have the same solution.
 ABI_CALLOC(wvec, (cplex, 0:iter))

 if (my_rank == master) then
   wvec(1, iter) = -one
   ABI_CALLOC(wmat1, (cplex, 0:iter, 0:iter))
   wmat1(1,:,iter) = -one
   wmat1(1,iter,:) = -one
   wmat1(1,iter,iter) = zero
   wmat1(:,0:iter-1, 0:iter-1) = diis%resmat(:, 0:iter-1, 0:iter-1, idat)

   call xhesv_cplex("U", cplex, iter+1, 1, wmat1, wvec, msg, ierr)
   ABI_CHECK(ierr == 0, msg)
   ABI_FREE(wmat1)
 end if

 ! Master broadcasts data.
 if (nprocs > 1) call xmpi_bcast(wvec, master, comm_spinorfft, ierr)

 if (cplex == 2) then
   call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), wvec(:,0), new_psi)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), wvec(:,0), new_residvec)
 else
   ABI_CALLOC(alphas, (2, 0:iter))
   alphas(1,:) = wvec(1,:)
   call cg_zgemv("N", npwsp, iter, diis%chain_phi(:,:,:,idat), alphas, new_psi)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv(:,:,:,idat), alphas, new_residvec)
   ABI_FREE(alphas)
 end if

 ABI_FREE(wvec)

end subroutine rmm_diis_solve
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

subroutine rmm_diis_eval_mats(diis, iter, ndat, me_g0, comm_spinorfft)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, ndat, me_g0, comm_spinorfft

 integer,parameter :: option2 = 2
 integer :: ii, ierr, idat
 real(dp) :: dotr, doti

 do idat=1,ndat
   do ii=0,iter
     ! <R_i|R_j>
     call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                    diis%chain_resv(:,:,ii,idat), diis%chain_resv(:,:,iter,idat), me_g0, xmpi_comm_self)
     if (ii == iter) doti = zero
     if (diis%cplex == 2) then
       diis%resmat(:, ii, iter, idat) = [dotr, doti]
     else
       diis%resmat(:, ii, iter, idat) = dotr
     end if

     ! <i|S|j> assume normalized wavefunctions for ii == iter
     if (ii == iter) then
       dotr = one; doti = zero
     else
       if (diis%usepaw == 0) then
         call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                        diis%chain_phi(:,:,ii,idat), diis%chain_phi(:,:,iter,idat), me_g0, xmpi_comm_self)
       else
         call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                        diis%chain_phi(:,:,ii,idat), diis%chain_sphi(:,:,iter,idat), me_g0, xmpi_comm_self)
       end if
     end if
     if (diis%cplex == 2) then
       diis%smat(:, ii, iter, idat) = [dotr, doti]
     else
       diis%smat(:, ii, iter, idat) = dotr
     end if
   end do

   if (xmpi_comm_size(comm_spinorfft) > 1) then
     call xmpi_sum(diis%resmat(:,0:iter,iter,idat), comm_spinorfft, ierr)
     call xmpi_sum(diis%smat(:,0:iter,iter, idat), comm_spinorfft, ierr)
   endif
 end do

end subroutine rmm_diis_eval_mats
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
 real(dp) :: dotr, doti

!$OMP PARALLEL DO PROVATE(dotr)
 do idat=1,ndat
   call dotprod_g(eig(idat), doti, istwf_k, npwsp, option1, ghc(:,idat), cg(:,idat), me_g0, xmpi_comm_self)
   if (usepaw == 1) then
     call dotprod_g(dotr, doti, istwf_k, npwsp, option1, gsc(:,idat), cg(:,idat), me_g0, xmpi_comm_self)
     eig(idat) = eig(idat) / dotr
   end if
 end do

 if (xmpi_comm_size(comm) > 1) call xmpi_sum(eig, comm, ierr)

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
!$OMP PARALLEL DO
   do idat=1,ndat
     residvecs(:,idat) = ghc(:,idat) - eig(idat) * gsc(:,idat)
   end do
 else
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

!$OMP PARALLEL DO PRIVATE(dotr, doti)
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
 real(dp) :: dotr, doti

!$OMP PARALLEL DO PRIVATE(dotr, doti)
 do idat=1,ndat
   call dotprod_g(dotr, doti, istwf_k, npwsp, option, cg1(:,idat), cg2(:,idat), me_g0, xmpi_comm_self)
   dots(:, idat) = [dotr, doti]
 end do

 if (xmpi_comm_size(comm) > 1) call xmpi_sum(dots, comm, ierr)

end subroutine cg_zdotg
!!***

!!****f* m_cgtools/cg_zdotg_lds
!! NAME
!!  cg_zdotg_lds
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

subroutine cg_zdotg_lds(istwf_k, npwsp, ndat, option, ld1, cg1, ld2, cg2, dots, me_g0, comm)

 integer,intent(in) :: istwf_k, npwsp, ndat, option, ld1, ld2, me_g0, comm
 real(dp),intent(in) :: cg1(2,ld1,ndat), cg2(2,ld2,ndat)
 real(dp),intent(out) :: dots(2,ndat)

 integer :: idat, ierr, ig1, ig2
 real(dp) :: dotr, doti

 do idat=1,ndat
   ig1 = 1 + (idat - 1) * ld1
   ig2 = 1 + (idat - 1) * ld2
   call dotprod_g(dotr, doti, istwf_k, npwsp, option, cg1(1,ig1,idat), cg2(1,ig2,idat), me_g0, xmpi_comm_self)
   dots(:, idat) = [dotr, doti]
 end do

 if (xmpi_comm_size(comm) > 1) call xmpi_sum(dots, comm, ierr)

end subroutine cg_zdotg_lds
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

 ABI_MALLOC(pcon, (npw))
 do idat=1,ndat
   call cg_precon(cg(:,idat), zero, istwf_k, kinpw, npw, nspinor, me_g0, optekin, pcon, vect(:,idat), comm)
 end do
 ABI_FREE(pcon)

 !if (xmpi_comm_size(comm) > 1) call xmpi_sum(dots, comm, ierr)

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

end module m_rmm_diis
!!***
