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
 use m_time,          only : timab
 use m_numeric_tools, only : pack_matrix
 use m_hide_lapack,   only : xheev
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
!!    gsc(2,*)=<g|s|c> matrix elements (s=overlap)
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

subroutine rmm_diis(cg, dtset, eig, enlx, gs_hamk, gsc, kinpw, mpi_enreg, nband, npw, nspinor, resid)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: nband, npw, nspinor
 real(dp),target,intent(inout) :: cg(2,npw*nspinor*nband), gsc(2,npw*nspinor*nband*dtset%usepaw)
 !real(dp),intent(in) :: kinpw(npw)
 real(dp),intent(inout) :: enlx(nband)
 real(dp),intent(inout) :: resid(nband)
 real(dp),intent(out) :: eig(nband)

!Local variables-------------------------------
 integer,parameter :: ndat1 = 1, type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: ii, jj
 integer :: iband, cpopt, sij_opt, is, ie, mcg, mgsc, use_subovl, istwf_k, optekin, usepaw, iter
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft
 real(dp),parameter :: rdummy = zero
 real(dp) :: dprod_r, dprod_i, lambda, r1
!arrays
 integer :: max_niter_band(nband)
 !integer, allocatable :: nline_bands(:)
 real(dp) :: tsec(2), eig_save(nband)
 real(dp),allocatable :: ghc(:,:), gvnlxc(:,:), pcon(:), chain_ug(:,:,:), chain_res(:,:,:)
 real(dp),allocatable :: residvec(:,:), kres(:,:), phi_now(:,:)
 real(dp),allocatable :: evec(:,:),evec_loc(:,:), subham(:), subovl(:)
 real(dp),allocatable :: h_ij(:,:,:), vnlx_ij(:,:,:)
 real(dp),allocatable :: diis_resmat(:,:,:), diis_smat(:,:,:), diis_workmat(:,:,:), diis_evec(:,:,:), diis_eig(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)
 !type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k
 mcg = npw * nspinor * nband; mgsc = npw * nspinor * nband * usepaw

 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft

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
 ABI_MALLOC(gvnlxc, (2, npw*nspinor))

 ! ====================
 ! Apply H to input cg
 ! ====================
 ABI_MALLOC(h_ij, (2, nband, nband))
 ABI_MALLOC(vnlx_ij, (2, nband, nband))

 do iband=1,nband
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   call getghc(cpopt, cg(:,is:ie), cwaveprj, ghc, gsc(:,is:ie), gs_hamk, gvnlxc,&
               rdummy, mpi_enreg, ndat1, dtset%prtvol, sij_opt, tim_getghc, type_calc0)

   ! <i|H|j> and <i|Vnl + FockACE|j>
   call cg_zgemv("C", npw*nspinor, nband, cg, ghc, h_ij(:,:,iband))
   !TODO: if istwf_k ==
   call cg_zgemv("C", npw*nspinor, nband, cg, gvnlxc, vnlx_ij(:,:,iband))
 end do

 ! Pack <i|H|j> to prepare call subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 call pack_matrix(h_ij, subham, nband, 2)
 ABI_FREE(h_ij)
 !call xmpi_sum(subham, comm_spinorfft, ierr)
 !call xmpi_sum(resid, comm_fft, ierr)

 use_subovl = 0
 if (use_subovl == 1) then
   ABI_MALLOC(subovl, (nband*(nband+1)))
 else
   ABI_MALLOC(subovl, (0))
 end if

 ! =================
 ! Subspace rotation
 ! =================
 call timab(585,1,tsec) !"vtowfk(subdiago)"
 ABI_MALLOC(evec, (2*nband, nband))
 call subdiago(cg, eig, evec, gsc, 0, 0, istwf_k, mcg, mgsc, nband, npw, nspinor, dtset%paral_kgb, &
               subham, subovl, use_subovl, usepaw, me_g0)
 call timab(585,2,tsec)

 ! Compute Vnl in new cg. Use subham as Workspace.
 call pack_matrix(vnlx_ij, subham, nband, 2)
 call cg_hprotate_and_get_diag(nband, subham, evec, enlx)
 ABI_FREE(vnlx_ij)

 ABI_FREE(subham)
 ABI_FREE(subovl)
 ABI_FREE(evec)

 ! ==========
 ! DIIS loop
 ! ==========
 optekin = 0 !; if (dtset%wfoptalg>=10) optekin = 1
 max_niter_band = dtset%nline

 ABI_ALLOCATE(pcon, (npw))
 ABI_MALLOC(chain_ug, (2, npw*nspinor, 0:dtset%nline))
 ABI_MALLOC(chain_res, (2, npw*nspinor, 0:dtset%nline))

 ABI_MALLOC(diis_resmat, (2, dtset%nline, dtset%nline))          ! <R_i|R_j>
 ABI_MALLOC(diis_workmat, (2, dtset%nline, dtset%nline))
 ABI_MALLOC(diis_smat, (2, dtset%nline, dtset%nline * usepaw))   ! <i|S|j> for PAW
 ABI_MALLOC(diis_evec, (2, dtset%nline, dtset%nline))

 ABI_MALLOC(residvec, (2, npw*nspinor))
 ABI_MALLOC(kres, (2, npw*nspinor))
 ABI_MALLOC(phi_now, (2, npw*nspinor))
 !ABI_MALLOC(nline_bands, (nband))
 !ABI_FREE(nline_bands)

 eig_save = eig

 do iband=1,nband
   ! phi_0 is stored in cg
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   phi_now = cg(:,is:ie)
   chain_ug(:,:,0) = phi_now

   do iter=1,1 !, max_niter_band(iband)
   !do iter=1,2
     call getghc(cpopt, phi_now, cwaveprj, ghc, gsc(:,is:ie), gs_hamk, gvnlxc, &
                 rdummy, mpi_enreg, ndat1, dtset%prtvol, sij_opt, tim_getghc, type_calc0)

     ! New approximated eigenvalue.
     !if (iter > 1)
     call dotprod_g(eig(iband), dprod_i, istwf_k, npw*nspinor, option1, ghc, &
                    phi_now, me_g0, comm_spinorfft)

     if (usepaw == 1) then
       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option1, gsc(:, is:), &
                      phi_now, me_g0, comm_spinorfft)
       eig(iband) = eig(iband) / dprod_r
     end if

     !if (abs(eig(iband) - eig_save(iband))) < (eaccuracy / nband_occ / dtset%nline) then
     !  max_niter_band(iband) = 0
     !end if

     ! Residual.
     if (usepaw == 1) then
       residvec = ghc - eig(iband) * gsc(:, is:ie)
     else
       residvec = ghc - eig(iband) * phi_now
     end if

     ! Store residual in R_{i-1}.
     chain_res(:, :, iter-1) = residvec
     call sqnorm_g(resid(iband), istwf_k, npw*nspinor, residvec, me_g0, comm_fft)
     !resid(iband) = sum(residvec**2)

     ! ============== CHECK FOR CONVERGENCE ========================

     ! Precondition residual, output in kres.
     !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
     !               optekin, pcon, kres, comm_bandspinorfft)
     kres = residvec

     ! <R_i|R_j> for j = iter-1
     do ii=0,iter-1
       call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option2, &
                      chain_res(:,:,ii), chain_res(:,:,iter-1), me_g0, comm_spinorfft)
       diis_resmat(:, ii+1, iter) = [dprod_r, dprod_i]
     end do

     if (iter == 1) then
       ! Compute lambda = -<R0|R1> / <R1|R1>
       !call dotprod_g(dprod_r, dprod_i, istwf_k, npw*nspinor, option1, &
       !               chain_res(:,:,0), kres, me_g0, comm_spinorfft)
       !call sqnorm_g(r1, istwf_k, npw*nspinor, chain_res(:,:,0), me_g0, comm_spinorfft)
       !lambda = -dprod_r / r1 ** 2
       !write(std_out, *)"lambda: ", lambda
       lambda = 0.1_dp

       ! phi_1 = phi_0 + lambda K.R_0
       chain_ug(:,:,iter) = phi_now + lambda * kres
       phi_now = chain_ug(:,:,iter)

     else
       ! Solve DIIS
       if (usepaw == 0) then
         !diis_workmat = diis_resmat
         !call xheev("V", "U", dtset%nline, diis_workmat, diis_eig)
       else
       end if
       ! Compute new trial wavefunction
       !call cg_lincom(iter, betas, npw*nspinor, chain_ug(:,:,:), phi_now)
       !call cg_lincom(iter, betas, npw*nspinor, chain_res(:,:,:), kres)
       phi_now = phi_now + lambda * kres
     end if

   end do ! iter

   ! Update wavefunction for this band
   cg(:,is:ie) = phi_now
 end do ! iband

 ABI_FREE(pcon)
 ABI_FREE(chain_ug)
 ABI_FREE(chain_res)
 ABI_FREE(diis_resmat)
 ABI_FREE(diis_workmat)
 ABI_FREE(diis_smat)
 ABI_FREE(diis_evec)

 ! Final cleanup
 ABI_FREE(ghc)
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
