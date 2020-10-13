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

 use defs_abitypes,   only : mpi_type
 use m_fstrings,      only : sjoin, itoa
 use m_time,          only : timab
 use m_numeric_tools, only : pack_matrix
 use m_hide_lapack,   only : xhegv_cplex !xheev,
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_axpby, pawcprj_copy
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 !use m_fock,          only : fock_set_ieigen,fock_set_getghc_call

 implicit none

 private
!!***

 public :: rmm_diis
!!***

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

subroutine rmm_diis(cg, dtset, eig, enlx, gs_hamk, gsc_all, mpi_enreg, nband, npw, nspinor, resid)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(in) :: mpi_enreg
 integer,intent(in) :: nband, npw, nspinor
 real(dp),intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband)
 real(dp),intent(inout) :: resid(nband)
 real(dp),intent(out) :: eig(nband)

!Local variables-------------------------------

 integer,parameter :: ndat1 = 1, type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0, level = 114
 integer,parameter :: use_subovl0 = 0
 integer :: ii, cplex, nline, ierr, prtvol
 integer :: iband, cpopt, sij_opt, is, ie, mcg, mgsc, istwf_k, optekin, usepaw, iter
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft
 real(dp),parameter :: rdummy = zero
 real(dp) :: dprod_r, dprod_i, lambda, rval
 character(len=500) :: msg
!arrays
 !integer :: max_niter_band(nband)
 real(dp) :: tsec(2), eig_save(nband)
 real(dp),allocatable :: ghc(:,:), gsc(:,:), gvnlxc(:,:), pcon(:), chain_phi(:,:,:), chain_res(:,:,:)
 real(dp),allocatable :: residvec(:,:), kres(:,:), phi_now(:,:)
 real(dp),allocatable :: evec(:,:), subham(:), h_ij(:,:,:)
 real(dp),allocatable :: diis_resmat(:,:,:), diis_smat(:,:,:), diis_eig(:), wmat1(:,:,:), wmat2(:,:,:)
 real(dp) :: subovl(use_subovl0)
 real(dp),pointer :: gsc_ptr(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 !type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; prtvol = dtset%prtvol
 mcg = npw * nspinor * nband; mgsc = npw * nspinor * nband * usepaw

 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft

 gsc_ptr => null()

 !eig_save = eig
 !eaccuracy =
 !nband_occ

 if (usepaw == 1) then
   ABI_MALLOC(cwaveprj, (gs_hamk%natom, nspinor*nband))
   call pawcprj_alloc(cwaveprj, 0, gs_hamk%dimcprj)
   sij_opt = 1 ! Recompute S
   cpopt = 0   ! <p_lmn|in> are computed and saved
   !cpopt = -1
 else
   sij_opt = 0
   cpopt = -1
 end if

 ABI_MALLOC(ghc, (2, npw*nspinor))
 ABI_MALLOC(gsc, (2, npw*nspinor*usepaw))
 ABI_MALLOC(gvnlxc, (2, npw*nspinor))

 ! ====================
 ! Apply H to input cg
 ! ====================
 ABI_MALLOC(h_ij, (2, nband, nband))

 do iband=1,nband
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)
   call getghc(cpopt, cg(:,is:ie), cwaveprj, ghc, gsc_ptr, gs_hamk, gvnlxc, &
               rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

   ! <i|H|j> and <i|Vnl + FockACE|j>
   if (istwf_k == 1) then
     call cg_zgemv("C", npw*nspinor, nband, cg, ghc, h_ij(:,:,iband))
   else
     MSG_ERROR("Only istwf_k == 1 is supported")
   end if
 end do

 ! Pack <i|H|j> to prepare call to subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 call pack_matrix(h_ij, subham, nband, 2)
 ABI_FREE(h_ij)
 !call xmpi_sum(subham, comm_spinorfft, ierr)
 !call xmpi_sum(resid, comm_fft, ierr)

 ! =================
 ! Subspace rotation
 ! =================
 call timab(585,1,tsec) !"vtowfk(subdiago)"
 ABI_MALLOC(evec, (2*nband, nband))
 call subdiago(cg, eig, evec, gsc_all, 0, 0, istwf_k, mcg, mgsc, nband, npw, nspinor, dtset%paral_kgb, &
               subham, subovl, use_subovl0, usepaw, me_g0)
 call timab(585,2,tsec)

 ABI_FREE(subham)
 ABI_FREE(evec)

 ! ==========
 ! DIIS loop
 ! ==========
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1
 !write(std_out,*)"optekin:", optekin
 nline = dtset%nline

 ABI_ALLOCATE(pcon, (npw))
 ABI_MALLOC(chain_phi, (2, npw*nspinor, 0:nline))
 ABI_MALLOC(chain_res, (2, npw*nspinor, 0:nline))

 ABI_MALLOC(diis_resmat, (2, 0:nline, 0:nline)) ! <R_i|R_j>
 ABI_MALLOC(diis_smat, (2, 0:nline, 0:nline))   ! <i|S|j>
 ABI_MALLOC(residvec, (2, npw*nspinor))
 ABI_MALLOC(kres, (2, npw*nspinor))
 ABI_MALLOC(phi_now, (2, npw*nspinor))

 eig_save = eig

#define _normalize

 do iband=1,nband
   ! Step 0:
   ! |phi_0> taken from cg
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   phi_now = cg(:,is:ie)
   chain_phi(:,:,0) = phi_now

   ! Compute H |phi_0>
   if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)
   call getghc(cpopt, phi_now, cwaveprj, ghc, gsc_ptr, gs_hamk, gvnlxc, &
               rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

   call dotprod_g(eig(iband), dprod_i, istwf_k, npw*nspinor, option1, ghc, &
                  phi_now, me_g0, comm_spinorfft)

   if (usepaw == 1) then
     call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option1, gsc_all(:,is:), &
                    phi_now, me_g0, comm_spinorfft)
     eig(iband) = eig(iband) / dprod_r
   end if

   ! Get residual R_0 = (H - e_0.S) |phi_0>
   if (usepaw == 1) then
     residvec = ghc - eig(iband) * gsc_all(:, is:ie)
   else
     residvec = ghc - eig(iband) * phi_now
   end if

   ! Store R_0 in chain array.
   chain_res(:,:,0) = residvec
   ! Evaluate <R0|R0> and <phi_0|S|phi_0>
   call sqnorm_g(resid(iband), istwf_k, npw*nspinor, residvec, me_g0, comm_fft)
   diis_resmat(:, 0, 0) = [resid(iband), zero]
   call sqnorm_g(rval, istwf_k, npw*nspinor, phi_now, me_g0, comm_spinorfft)
   diis_smat(:, 0, 0) = [rval, zero]

   kres = residvec
   ! Precondition R0, output in kres = K.R_0
   !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
   call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                  optekin, pcon, kres, comm_bandspinorfft)
   !write(std_out,*)"kres(1:2):", kres(:, 1:2)

   if (prtvol == -level) then
     write(msg,'(2(a5),a14,a12)')'iter', "band", "eigen", "resid"
     call wrtout(std_out, msg, "PERS")
     write(msg,'(2(i5),f14.6,es12.4)')0, iband, eig(iband) * Ha_eV, resid(iband)
     call wrtout(std_out, msg, "PERS")
   end if

iter_loop: do iter=1,nline

     if (iter == 1) then
       ! Compute H |K.R_0>
       call getghc(cpopt, kres, cwaveprj, ghc, gsc, gs_hamk, gvnlxc, &
                   rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

       ! Compute residual: (H - e_0.S) |K.R_0>
       if (usepaw == 1) then
         residvec = ghc - eig(iband) * gsc
       else
         residvec = ghc - eig(iband) * kres
       end if

       ! Compute lambda that minimizes the norm of the residual.
       ! lambda = - Re{<R_0|(H - e_0.S)}|K.R_0>} / |(H - e_0.S) K.R_0|**2
       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option1, &
                      chain_res(:,:,0), residvec, me_g0, comm_spinorfft)
       call sqnorm_g(rval, istwf_k, npw*nspinor, residvec, me_g0, comm_spinorfft)

       lambda = -dprod_r / rval
       !write(std_out, *)"lambda, dprod_r, rval ", lambda, dprod_r, rval
       !if (lambda < tol1 .or. lambda > one) then
       !   lambda = 0.1
       !end if
       !lambda = one

       ! phi_1 = phi_0 + lambda K.R_0
       phi_now = chain_phi(:,:,0) + lambda * kres
#ifdef _normalize
       call sqnorm_g(rval, istwf_k, npw*nspinor, phi_now, me_g0, comm_spinorfft)
       phi_now = phi_now / sqrt(rval)
#endif
       chain_phi(:,:,1) = phi_now

     else
       ! Solve DIIS to get phi_now for iter >= 2
       cplex = 2
       ABI_MALLOC(diis_eig, (0:iter-1))
       ABI_MALLOC(wmat1, (cplex, 0:iter-1, 0:iter-1))
       ABI_MALLOC(wmat2, (cplex, 0:iter-1, 0:iter-1))

       wmat1 = diis_resmat(:, 0:iter-1, 0:iter-1)
       wmat2 = diis_smat(:, 0:iter-1, 0:iter-1)
       !do ii=0,iter-1; write(std_out, *) "diis_resmat:", wmat1(:,ii,:); end do
       !do ii=0,iter-1; write(std_out, *) "diis_smat:", wmat2(:,ii,:); end do

       call xhegv_cplex(1, "V", "U", cplex, iter, wmat1, wmat2, diis_eig, msg, ierr)
       !write(std_out,*)"diis_eig:", diis_eig(0)
       !write(std_out,*)"RE diis_vec  :", wmat1(1,:,0)
       !write(std_out,*)"IMAG diis_vec:", wmat1(2,:,0)
       if (ierr /= 0) then
         call wrtout(std_out, sjoin("xhegv failed with:", msg, ch10, "at iter: ", itoa(iter), "exit iter_loop!"))
         ABI_FREE(wmat1)
         ABI_FREE(wmat2)
         ABI_FREE(diis_eig)
         exit iter_loop
       end if

       ! Take linear combination
       call cg_zgemv("N", npw*nspinor, iter, chain_phi, wmat1(:,:,0), phi_now)
       call cg_zgemv("N", npw*nspinor, iter, chain_res, wmat1(:,:,0), residvec)
#ifdef _normalize
       !call sqnorm_g(rval, istwf_k, npw*nspinor, phi_now, me_g0, comm_spinorfft)
       !phi_now = phi_now / sqrt(rval)
#endif

       ! Precondition residual, output in kres.
       kres = residvec
       !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
       call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                      optekin, pcon, kres, comm_bandspinorfft)

       ! Compute phi_now with lambda obtained in step #1
       phi_now = phi_now + lambda * kres
#ifdef _normalize
       call sqnorm_g(rval, istwf_k, npw*nspinor, phi_now, me_g0, comm_spinorfft)
       phi_now = phi_now / sqrt(rval)
#endif
       chain_phi(:,:,iter) = phi_now

       ABI_FREE(wmat1)
       ABI_FREE(wmat2)
       ABI_FREE(diis_eig)
     end if

     ! Compute H |phi_now>
     if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)
     call getghc(cpopt, phi_now, cwaveprj, ghc, gsc_ptr, gs_hamk, gvnlxc, &
                 rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

     ! New approximated eigenvalue.
!if (iter == 1) then
!if (iter <= 2) then
     call dotprod_g(eig(iband), dprod_i, istwf_k, npw*nspinor, option1, ghc, &
                    phi_now, me_g0, comm_spinorfft)

     if (usepaw == 1) then
       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option1, gsc_all(:, is:), &
                      phi_now, me_g0, comm_spinorfft)
       eig(iband) = eig(iband) / dprod_r
     end if
!end if

     !if (abs(eig(iband) - eig_save(iband))) < (eaccuracy / nband_occ / dtset%nline) then
     !  max_niter_band(iband) = 0
     !end if

     ! Residual.
     if (usepaw == 1) then
       residvec = ghc - eig(iband) * gsc_all(:, is:ie)
     else
       residvec = ghc - eig(iband) * phi_now
     end if

     ! Store residual R_iter.
     chain_res(:, :, iter) = residvec
     call sqnorm_g(resid(iband), istwf_k, npw*nspinor, residvec, me_g0, comm_fft)

     call dotprod_g(enlx(iband), dprod_i, istwf_k, npw*nspinor, option1, &
                    phi_now, gvnlxc, me_g0, comm_bandspinorfft)

     ! ============== CHECK FOR CONVERGENCE ========================
     if (prtvol == -level) then
       write(msg,"(2(i5),f14.6,es12.4)")iter, iband, eig(iband) * Ha_eV, resid(iband)
       call wrtout(std_out, msg, 'PERS')
     end if

     !if (resid(iband) < dtset%tolwfr) then
     if (resid(iband) < tol12) then
       if (prtvol == -level) then
         write(msg, '(2a,i4,a,i2,a,es12.4,a)' ) ch10, &
          ' rmm_diis: band ',iband,' converged after ',iter,' iterations with resid: ',resid(iband), ch10
         call wrtout(std_out, msg , "PERS")
       end if
       exit iter_loop
     end if

     ! <R_i|R_j> for j = iter
     ! <i|S|j>
     do ii=0,iter
       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option2, &
                      chain_res(:,:,ii), chain_res(:,:,iter), me_g0, comm_spinorfft)
       if (ii == iter) dprod_i = zero
       diis_resmat(:, ii, iter) = [dprod_r, dprod_i]

       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option2, &
                      chain_phi(:,:,ii), chain_phi(:,:,iter), me_g0, comm_spinorfft)
       if (ii == iter) dprod_i = zero
       diis_smat(:, ii, iter) = [dprod_r, dprod_i]
     end do

   end do iter_loop

   ! Update wavefunction for this band
   cg(:,is:ie) = phi_now
 end do ! iband

 ! Final cleanup
 ABI_FREE(pcon)
 ABI_FREE(chain_phi)
 ABI_FREE(chain_res)
 ABI_FREE(diis_resmat)
 ABI_FREE(diis_smat)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(residvec)
 ABI_FREE(kres)
 ABI_FREE(phi_now)
 ABI_FREE(gvnlxc)

 if (usepaw == 1) then
   call pawcprj_free(cwaveprj)
   ABI_FREE(cwaveprj)
 end if

end subroutine rmm_diis
!!***

end module m_rmm_diis
!!***
