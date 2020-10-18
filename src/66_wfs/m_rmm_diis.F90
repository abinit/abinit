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
 use m_fstrings,      only : sjoin, itoa
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
   integer :: max_niter
   integer :: npwsp
   type(pair_list) :: stats

   real(dp),allocatable :: hist_ene(:) ! (0:max_niter+1)
   real(dp),allocatable :: hist_resid(:) !(0:dtset%nline)

   real(dp),allocatable :: resmat(:,:,:)
   real(dp),allocatable :: smat(:,:,:)

   real(dp),allocatable :: chain_phi(:,:,:)
   real(dp),allocatable :: chain_sphi(:,:,:)
   real(dp),allocatable :: chain_resv(:,:,:)

 contains
   procedure :: free => rmm_diis_free
   procedure :: solve => rmm_diis_solve
   procedure :: eval_mats => rmm_diis_eval_mats
   procedure :: exit => rmm_diis_exit

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
 integer,parameter :: ndat1 = 1, type_calc0 = 0, option1 = 1, option2 = 2
 integer,parameter :: tim_getghc = 0, level = 432, use_subovl0 = 0
 integer :: ig, ig0, ib, cplex, ierr, prtvol, bsize, nblocks, iblock, npwsp, ndat, ib_start, ib_stop
 integer :: iband, cpopt, sij_opt, igs, ige, mcg, mgsc, istwf_k, optekin, usepaw, iter, max_niter
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft, nbocc, ii
 logical :: end_with_trial_step, use_prec_steepest_descent
 real(dp),parameter :: rdummy = zero
 real(dp) :: dotr, doti, lambda, rval, accuracy_ene
 !real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
!arrays
 integer :: niter_band(nband)
 real(dp) :: tsec(2) !, hist_ene(0:dtset%nline), hist_resid(0:dtset%nline)
 real(dp),target :: fake_gsc(0,0)
 real(dp) :: subovl(use_subovl0)
 real(dp),allocatable :: ghc(:,:), gsc(:,:), gvnlxc(:,:), pcon(:)
 real(dp),allocatable :: residvec(:,:), kres(:,:), phi_now(:,:), evec(:,:), subham(:), h_ij(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_ptr(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(rmm_diis_t) :: diis
 type(yamldoc_t) :: ydoc

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; prtvol = dtset%prtvol
 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
 npwsp = npw * nspinor
 cplex = 2

 ! =======================================
 ! Apply H to input cg to compute <i|H|j>
 ! =======================================
 gsc_ptr => fake_gsc
 cpopt = -1
 if (usepaw == 1) then
   sij_opt = 1 ! Recompute S
   !cpopt = 0  ! <p_lmn|in> are computed and saved
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 else
   sij_opt = 0
 end if

 ! Treat states in groups of bsize bands to be able to use zgemm.
 bsize = 4
 if (dtset%userib /= 0) bsize = abs(dtset%userib)
 nblocks = nband / bsize; if (mod(nband, bsize) /= 0) nblocks = nblocks + 1

 ABI_MALLOC(ghc, (2, npwsp*bsize))
 ABI_MALLOC(gsc, (2, npwsp*bsize*usepaw))
 ABI_MALLOC(gvnlxc, (2, npwsp*bsize))
 ABI_CALLOC(h_ij, (2, nband, nband))

 call timab(1633, 1, tsec) !"rmm_diis:build_hij
 !call cwtime(cpu, wall, gflops, "start")

 do iblock=1,nblocks
   igs = 1 + (iblock - 1) * npwsp * bsize
   ige = igs + npwsp * bsize - 1; ndat = bsize
   if (ige > npwsp * nband) then
     ige = npwsp * nband; ndat = (ige - igs + 1) / npwsp
   end if
   ib_start = 1 + (iblock - 1) * bsize; ib_stop = min(iblock * bsize, nband)

   if (usepaw == 1) gsc_ptr => gsc_all(:,igs:ige)
   call getghc(cpopt, cg(:,igs:ige), cprj_dum, ghc, gsc_ptr, gs_hamk, gvnlxc, &
               rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

   ! Compute <i|H|j> for i=1,nband and all j in block
   call cg_zgemm("C", "N", npwsp, nband, ndat, cg, ghc, h_ij(:,:,ib_start))

   if (istwf_k /= 1) then
     do iband=ib_start, ib_stop
       h_ij(:,:,iband) = two * h_ij(:,:,iband)

       if (istwf_k == 2 .and. me_g0 == 1) then
         ! Gamma k-point and I have G=0. Remove double counting term.
         ig = 1 + (iband - ib_start) * npwsp
         do ib=1,nband
           ig0 = 1 + npwsp * (ib - 1)
           h_ij(1,ib,iband) = h_ij(1,ib,iband) - cg(1,ig0) * ghc(1,ig)
         end do
       end if
       ! Force real matrix.
       h_ij(2,:,iband) = zero
     end do
   end if

 end do

 !call cwtime_report(" build_hij", cpu, wall, gflops)

 ! Pack <i|H|j> to prepare call to subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 do iband=1,nband
   h_ij(2,iband,iband) = zero ! Force diagonal elements to be real
 end do
 call pack_matrix(h_ij, subham, nband, 2)
 ABI_FREE(h_ij)
 call xmpi_sum(subham, comm_spinorfft, ierr)
 call timab(1633, 2, tsec) !"rmm_diis:build_hij
 !call cwtime_report(" pack_hij", cpu, wall, gflops)

 ! =================
 ! Subspace rotation
 ! =================
 call timab(585, 1, tsec) !"vtowfk(subdiago)"
 ABI_MALLOC(evec, (2*nband, nband))
 mcg = npwsp * nband; mgsc = npwsp * nband * usepaw
 call subdiago(cg, eig, evec, gsc_all, 0, 0, istwf_k, mcg, mgsc, nband, npw, nspinor, dtset%paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)
 call timab(585, 2, tsec)
 !call cwtime_report(" subdiago", cpu, wall, gflops)

 ABI_FREE(subham)
 ABI_FREE(evec)

 ! =================
 ! Prepare DIIS loop
 ! =================
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1
 !write(std_out,*)"optekin:", optekin
 max_niter = dtset%nline
 !max_niter = dtset%nline - 1

 diis = rmm_diis_new(usepaw, istwf_k, npwsp, max_niter)

 ABI_REMALLOC(ghc, (2, npwsp))
 ABI_REMALLOC(gsc, (2, npwsp*usepaw))
 ABI_REMALLOC(gvnlxc, (2, npwsp))

 ! Remalloc here because blocking algorithm is not yet implemented in the DIIS part.
 ABI_MALLOC(pcon, (npw))
 ABI_MALLOC(residvec, (2, npwsp))
 ABI_MALLOC(kres, (2, npwsp))
 ABI_MALLOC(phi_now, (2, npwsp))

 ! Reduce number of niter iterations for empty states.
 niter_band = max_niter
 ! TODO: Don't reduce niter if MD
 !if dtset%
 !if dtset%occopt == 2
 do iband=1,nband
   if (occ(iband) < tol3) niter_band(iband) = max(max_niter / 2, 2)
   !if (occ(iband) < tol3) niter_band(iband) = max(1 + max_niter / 2, 2)
 end do

 ! Define accuracy_ene for SCF
 nbocc = count(occ > zero)
 accuracy_ene = zero
 if (dtset%iscf > 0) then
   if (dtset%toldfe /= 0) then
     accuracy_ene = dtset%toldfe / nbocc / four
   else
     ! User is not using toldfe to stop the SCF cycle
     ! so we are forced to hardcode a tolerance for the absolute diff in the KS eigenvalue.
     accuracy_ene = tol8 / nbocc / four
   end if
 end if

 call timab(1634, 1, tsec) ! "rmm_diis:band_opt"
 iband_loop: do iband=1,nband
   ! Compute H |phi_0> from subdiago cg.
   igs = 1 + npwsp * (iband - 1); ige = igs - 1 + npwsp
   phi_now = cg(:,igs:ige)
   if (usepaw == 1) gsc_ptr => gsc_all(:,igs:ige)

   call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                        eig(iband), resid(iband), enlx(iband), residvec, normalize=.False.)

   ! Compute <R0|R0> and <phi_0|S|phi_0>
   ! Assume input cg are already S-normalized.
   call sqnorm_g(resid(iband), istwf_k, npwsp, residvec, me_g0, comm_fft)
   diis%hist_resid(0) = resid(iband)
   diis%resmat(:, 0, 0) = [resid(iband), zero]
   diis%smat(:, 0, 0) = [one, zero]

   ! BAND LOCKING for NSCF or SCF if tolwfr > 0 is used.
   if (resid(iband) < dtset%tolwfr .or. sqrt(abs(resid(iband))) < accuracy_ene) then
      call diis%stats%increment("locked_bands", 1)
      cycle iband_loop
   end if

   diis%hist_ene(0) = eig(iband)
   diis%chain_resv(:,:,0) = residvec

   ! Store |phi_0>, |S phi_0>.
   diis%chain_phi(:,:,0) = phi_now
   if (usepaw == 1) diis%chain_sphi(:,:,0) = gsc_ptr

   ! Precondition |R_0>, output in kres = |K R_0>
   kres = residvec
   !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
   call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                  optekin, pcon, kres, comm_spinorfft)
   !write(std_out,*)"kres(1:2):", kres(:, 1:2)

   if (prtvol == -level) then
     call wrtout(std_out, &
           sjoin("<BEGIN RMM-DIIS, istep:", itoa(istep), ", ikpt:", itoa(ikpt), ", spin: ", itoa(isppol), ">"))
     write(msg,'(1a, 2(a5),2(a14),a12)')"#", 'iter', "band", "eigen_eV", "eigde_meV", "resid"
     call wrtout(std_out, msg)
     write(msg,'(1x, 2(i5),2(f14.6,1x),es12.4)')0, iband, eig(iband) * Ha_eV, zero, resid(iband)
     call wrtout(std_out, msg)
   end if

   iter_loop: do iter=1,niter_band(iband)

     if (iter == 1) then
       ! Compute H |K R_0>
       call getghc(cpopt, kres, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
                   rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

       ! Compute residual: (H - e_0 S) |K R_0>
       if (usepaw == 1) then
         residvec = ghc - eig(iband) * gsc
       else
         residvec = ghc - eig(iband) * kres
       end if

       ! Preconditioned steepest descent:
       !
       !    |phi_1> = |phi_0> + lambda |K R_0>
       !
       ! where lambda minimizes the norm of the residual
       ! (not the Rayleigh quotient as in the original Kresse's paper)
       !
       !    lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
       !
       call dotprod_g(dotr, doti, istwf_k, npwsp, option1, &
                      diis%chain_resv(:,:,0), residvec, me_g0, comm_spinorfft)
       call sqnorm_g(rval, istwf_k, npwsp, residvec, me_g0, comm_spinorfft)

       lambda = -dotr / rval
       phi_now = diis%chain_phi(:,:,0) + lambda * kres
       !write(std_out, *)"lambda, dotr, rval ", lambda, dotr, rval

     else
       ! Solve DIIS equations to get phi_now for iter > 1
       call diis%solve(iter, npwsp, phi_now, residvec, comm_spinorfft)

       ! Precondition residual, output in kres.
       kres = residvec
       !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
       call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                      optekin, pcon, kres, comm_spinorfft)

       ! Compute phi_now with the lambda obtained at iteration #1
       phi_now = phi_now + lambda * kres
     end if

     ! Compute H |phi_now>.
     if (usepaw == 1) gsc_ptr => gsc_all(:,igs:ige)

     call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                          eig(iband:), resid(iband:), enlx(iband:), residvec, normalize=.True.)

     ! Store residual R_i for i == iter and evaluate new enlx for NC.
     diis%hist_ene(iter) = eig(iband)
     diis%hist_resid(iter) = resid(iband)
     diis%chain_phi(:,:,iter) = phi_now
     diis%chain_resv(:, :, iter) = residvec
     if (usepaw == 1) diis%chain_sphi(:,:,iter) = gsc_ptr

     ! ============== CHECK FOR CONVERGENCE ========================
     if (prtvol == -level) then
       write(msg,"(1x, 2(i5), 2(f14.6),es12.4)") &
         iter, iband, diis%hist_ene(iter) * Ha_eV, (diis%hist_ene(iter) - diis%hist_ene(iter-1)) * Ha_meV, resid(iband)
       call wrtout(std_out, msg)
     end if

     if (diis%exit(iter, niter_band(iband), accuracy_ene, dtset)) then
       if (prtvol == -level) then
         write(msg, '(a,i4,a,i2,a,es12.4,a)' )&
          ' band: ',iband,' converged after: ',iter,' iterations with resid: ',resid(iband), ch10
         call wrtout(std_out, sjoin("<END RMM-DIIS, msg='", msg, "'>"))
       end if
       exit iter_loop
     end if

     ! Compute <R_i|R_j> and <i|S|j> for j = iter
     if (iter /= niter_band(iband)) call diis%eval_mats(iter, me_g0, comm_spinorfft)

   end do iter_loop

   ! End with trial step but only then if performed all the iterations (no exit from diis%exit)
   !end_with_trial_step = .True.
   end_with_trial_step = iter == niter_band(iband) + 1
   !end_with_trial_step = .False.

   if (end_with_trial_step) then
     kres = residvec
     call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                    optekin, pcon, kres, comm_spinorfft)

     use_prec_steepest_descent = .False.
     if (iband <= nbocc + 4 .and. resid(iband) > tol10) use_prec_steepest_descent = .True.
     !write(std_out, *)" Performing last trial step with use_prec_steepest_descent: ", use_prec_steepest_descent
     !use_prec_steepest_descent = .True.
     !use_prec_steepest_descent = .False.

     if (use_prec_steepest_descent) then
       ! Final Preconditioned steepest descent:
       ! Compute H |K R_0>
       call getghc(cpopt, kres, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
                   rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

       ! Compute residual: (H - e_0 S) |K R_0>
       if (usepaw == 1) then
         residvec = ghc - eig(iband) * gsc
       else
         residvec = ghc - eig(iband) * kres
       end if

       ii = niter_band(iband)
       call dotprod_g(dotr, doti, istwf_k, npwsp, option1, &
                      diis%chain_resv(:,:,ii), residvec, me_g0, comm_spinorfft)
       call sqnorm_g(rval, istwf_k, npwsp, residvec, me_g0, comm_spinorfft)

       lambda = -dotr / rval
       phi_now = diis%chain_phi(:,:,ii) + lambda * kres

     else
       !phi_now = phi_now + 0.1 * kres
       phi_now = phi_now + lambda * kres
     endif

     if (usepaw == 1) gsc_ptr => gsc_all(:,igs:ige)

     call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                          eig(iband:), resid(iband:), enlx(iband), residvec, normalize=.True.)

     if (prtvol == -level) then
       write(msg,"(1x, 2(i5), 2(f14.6),es12.4)") &
         iter+1, iband, eig(iband) * Ha_eV, (eig(iband) - diis%hist_ene(iter)) * Ha_meV, resid(iband)
       call wrtout(std_out, msg)
     end if

   end if

   ! Update wavefunction for this band
   cg(:,igs:ige) = phi_now

   if (prtvol == -level .and. iter == niter_band(iband) + 1) then
     call wrtout(std_out, "<END RMM-DIIS, msg='All iterations performed'>")
   end if

 end do iband_loop
 call timab(1634, 2, tsec) !"rmm_diis:band_opt"
 !call cwtime_report(" rmm_diis:band_opt", cpu, wall, gflops)

 !call yaml_write_and_free_dict('RMM-DIIS', dict, std_out)
 ydoc = yamldoc_open('RMM-DIIS', with_iter_state=.False.)
 call ydoc%add_dict("stats", diis%stats)
 call ydoc%write_and_free(std_out)
 call diis%stats%free()

 ! Final cleanup
 ABI_FREE(pcon)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(residvec)
 ABI_FREE(kres)
 ABI_FREE(phi_now)
 ABI_FREE(gvnlxc)
 call diis%free()

end subroutine rmm_diis
!!***

logical function rmm_diis_exit(diis, iter, niter_band, accuracy_ene, dtset) result(ans)

 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iter, niter_band
 real(dp),intent(in) :: accuracy_ene
 type(dataset_type),intent(in) :: dtset

 real(dp) :: resid, deltae, deold

 resid = diis%hist_resid(iter)
 deold = diis%hist_ene(1) - diis%hist_ene(0)
 deltae = diis%hist_ene(iter) - diis%hist_ene(iter-1)

 ans = .False.
 ! Abinit default is 0.005 (0.3 V)
 !if (abs(deltae) < dtset%tolrde * abs(deold) .and. iter /= niter_band) then
 if (abs(deltae) < 12 * dtset%tolrde * abs(deold) .and. iter /= niter_band) then

 !if (abs(deltae) < 0.3_dp * abs(deold) .and. iter /= niter_band) then
 ! This one seems much better.
 !if (abs(deltae) < 0.03_dp * abs(deold) .and. iter /= niter_band) then
   call diis%stats%increment("deltae < 12 * tolrde * deold", 1)
   ans = .True.
   return
 end if

 if (dtset%iscf < 0) then
   ! This is the only condition available for NSCF run.
   if (resid < dtset%tolwfr) then
     call diis%stats%increment("resid < tolwfr", 1)
     ans = .True.
   end if

 else
   if (sqrt(abs(resid)) < accuracy_ene) then
     call diis%stats%increment("resid < accuracy_ene", 1)
     ans = .True.  ! Similar to EDIFF/NBANDS/4
   end if
   if (resid < dtset%tolwfr) then
     call diis%stats%increment("resid < tolwfr", 1)
     ans = .True.
   end if
 end if

end function rmm_diis_exit
!!***

subroutine getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_now, ghc, gsc, mpi_enreg, prtvol, &
                           eig, resid, enlx, residvec, normalize)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: npw, nspinor, ndat, prtvol
 real(dp),intent(inout) :: phi_now(2, npw*nspinor*ndat)
 real(dp),intent(out) :: ghc(2,npw*nspinor*ndat), gsc(2,npw*nspinor*ndat*gs_hamk%usepaw)
 type(mpi_type),intent(in) :: mpi_enreg
 real(dp),intent(out) :: eig(ndat), resid(ndat), enlx(ndat)
 real(dp),intent(out) :: residvec(2, npw*nspinor*ndat)
 logical,optional,intent(in) :: normalize

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: istwf_k, usepaw, cpopt, sij_opt, idat, is, ie, npwsp
 integer :: me_g0, comm_spinorfft, comm_fft
 real(dp),parameter :: rdummy = zero
 logical :: normalize_
 real(dp) :: dotr, doti
!arrays
 real(dp),allocatable :: gvnlxc(:, :)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k
 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft
 normalize_ = .True.; if (present(normalize)) normalize_ = normalize
 npwsp = npw * nspinor

 cpopt = -1
 if (usepaw == 1) then
   sij_opt = 1 ! Recompute S
   !cpopt = 0  ! <p_lmn|in> are computed and saved
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 else
   sij_opt = 0
 end if

 if (usepaw == 0 .and. normalize_) then
   ! NC normalization.
   call cgnc_normalize(npwsp, ndat, phi_now, istwf_k, me_g0, comm_spinorfft)
 end if

 ABI_MALLOC(gvnlxc, (2, npwsp*ndat))

 ! Compute H |phi_now>
 !call fock_set_ieigen(gs_hamk%fockcommon, iband)
 call getghc(cpopt, phi_now, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
             rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

 if (usepaw == 1 .and. normalize_) then
   ! PAW normalization must be done here.
   call cgpaw_normalize(npwsp, ndat, phi_now, gsc, istwf_k, me_g0, comm_spinorfft)
 end if

 do idat=1,ndat
   is = 1 + (idat - 1) * npwsp; ie = is - 1 + npwsp
   call dotprod_g(eig(idat), doti, istwf_k, npwsp, option1, ghc(:,is), phi_now(:,is), me_g0, comm_spinorfft)

   if (usepaw == 1) then
     call dotprod_g(dotr, doti, istwf_k, npwsp, option1, gsc(:,is), phi_now(:,is), me_g0, comm_spinorfft)
     eig(idat) = eig(idat) / dotr
   end if

   ! Residual.
   if (usepaw == 1) then
     residvec(:,is:ie) = ghc(:,is:ie) - eig(idat) * gsc(:,is:ie)
   else
     residvec(:,is:ie) = ghc(:,is:ie) - eig(idat) * phi_now(:,is:ie)
   end if

   ! Store residual R_i for i == iter and evaluate new enlx for NC.
   call sqnorm_g(resid(idat), istwf_k, npwsp, residvec(1,is), me_g0, comm_spinorfft)

   if (usepaw == 0) then
     call dotprod_g(enlx(idat), doti, istwf_k, npwsp, option1, phi_now(1,is), gvnlxc(1,is), me_g0, comm_spinorfft)
   end if
 end do

 ABI_FREE(gvnlxc)

end subroutine getghc_eigresid
!!!***

type(rmm_diis_t) function rmm_diis_new(usepaw, istwf_k, npwsp, max_niter) result(diis)

 integer,intent(in) :: usepaw, istwf_k, npwsp, max_niter

 diis%usepaw = usepaw
 diis%istwf_k = istwf_k
 diis%cplex = 2; if (istwf_k == 2) diis%cplex = 1
 diis%npwsp = npwsp
 diis%max_niter = max_niter

 ABI_MALLOC(diis%hist_ene, (0:max_niter+1))
 ABI_MALLOC(diis%hist_resid, (0:max_niter+1))

 ABI_MALLOC(diis%chain_phi, (2, npwsp, 0:max_niter))
 ABI_MALLOC(diis%chain_sphi, (2, npwsp*usepaw, 0:max_niter))
 ABI_MALLOC(diis%chain_resv, (2, npwsp, 0:max_niter))
 ABI_CALLOC(diis%resmat, (2, 0:max_niter, 0:max_niter)) ! <R_i|R_j>
 ABI_CALLOC(diis%smat, (2, 0:max_niter, 0:max_niter))   ! <i|S|j>

end function rmm_diis_new

subroutine rmm_diis_free(diis)
 class(rmm_diis_t),intent(inout) :: diis

 ABI_SFREE(diis%hist_ene)
 ABI_SFREE(diis%hist_resid)
 ABI_SFREE(diis%chain_phi)
 ABI_SFREE(diis%chain_sphi)
 ABI_SFREE(diis%chain_resv)
 ABI_SFREE(diis%resmat)
 ABI_SFREE(diis%smat)
 call diis%stats%free()

end subroutine rmm_diis_free

subroutine rmm_diis_solve(diis, iter, npwsp, phi_now, residvec, comm_spinorfft)

 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iter, npwsp, comm_spinorfft
 real(dp),intent(inout) :: phi_now(2, npwsp), residvec(2, npwsp)

 integer,parameter :: master = 0
 integer :: cplex, ierr, nprocs, my_rank
 real(dp),allocatable :: diis_eig(:)
 real(dp),allocatable :: wmat1(:,:,:), wmat2(:,:,:), wvec(:,:), alphas(:,:)
 character(len=500) :: msg
 logical,parameter :: try_to_solve_eigproblem = .False.

 cplex = diis%cplex
 my_rank = xmpi_comm_rank(comm_spinorfft); nprocs = xmpi_comm_size(comm_spinorfft)

 if (try_to_solve_eigproblem) then
   ABI_MALLOC(diis_eig, (0:iter-1))
   ABI_MALLOC(wmat1, (cplex, 0:iter-1, 0:iter-1))
   ABI_MALLOC(wmat2, (cplex, 0:iter-1, 0:iter-1))
   wmat1 = diis%resmat(:, 0:iter-1, 0:iter-1)
   wmat2 = diis%smat(:, 0:iter-1, 0:iter-1)
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
   call cg_zgemv("N", npwsp, iter, diis%chain_phi, wmat1(:,:,0), phi_now)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv, wmat1(:,:,0), residvec)

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
   wmat1(:,0:iter-1, 0:iter-1) = diis%resmat(:, 0:iter-1, 0:iter-1)

   call xhesv_cplex("U", cplex, iter+1, 1, wmat1, wvec, msg, ierr)
   ABI_CHECK(ierr == 0, msg)
   ABI_FREE(wmat1)
 end if

 ! Master broadcasts data.
 if (nprocs > 1) call xmpi_bcast(wvec, master, comm_spinorfft, ierr)

 if (cplex == 2) then
   call cg_zgemv("N", npwsp, iter, diis%chain_phi, wvec(:,0), phi_now)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv, wvec(:,0), residvec)
 else
   ABI_CALLOC(alphas, (2, 0:iter))
   alphas(1,:) = wvec(1,:)
   call cg_zgemv("N", npwsp, iter, diis%chain_phi, alphas, phi_now)
   call cg_zgemv("N", npwsp, iter, diis%chain_resv, alphas, residvec)
   ABI_FREE(alphas)
 end if

 ABI_FREE(wvec)

end subroutine rmm_diis_solve

subroutine rmm_diis_eval_mats(diis, iter, me_g0, comm_spinorfft)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iter, me_g0, comm_spinorfft

 integer,parameter :: option2 = 2
 integer :: ii, ierr
 real(dp) :: dotr, doti

 do ii=0,iter
   ! <R_i|R_j>
   call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                  diis%chain_resv(:,:,ii), diis%chain_resv(:,:,iter), me_g0, xmpi_comm_self)
   if (ii == iter) doti = zero
   diis%resmat(:, ii, iter) = [dotr, doti]

   ! <i|S|j> assume normalized wavefunctions for ii == iter
   if (ii == iter) then
     dotr = one; doti = zero
   else
     if (diis%usepaw == 0) then
       call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                      diis%chain_phi(:,:,ii), diis%chain_phi(:,:,iter), me_g0, xmpi_comm_self)
     else
       call dotprod_g(dotr, doti, diis%istwf_k, diis%npwsp, option2, &
                      diis%chain_phi(:,:,ii), diis%chain_sphi(:,:,iter), me_g0, xmpi_comm_self)
     end if
   end if
   diis%smat(:, ii, iter) = [dotr, doti]
 end do

 if (xmpi_comm_size(comm_spinorfft) > 1) then
   call xmpi_sum(diis%resmat(:,0:iter,iter), comm_spinorfft, ierr)
   call xmpi_sum(diis%smat(:,0:iter,iter), comm_spinorfft, ierr)
 endif

end subroutine rmm_diis_eval_mats
!!***

end module m_rmm_diis
!!***
