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
 use m_hide_lapack,   only : xhegv_cplex
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
   integer :: ndat
   integer :: istwf_k
   integer :: cplex
   integer :: nline
   integer :: npw
   integer :: nspinor

   real(dp),allocatable :: resmat(:,:,:)
   real(dp),allocatable :: smat(:,:,:)

   real(dp),allocatable :: chain_phi(:,:,:)
   real(dp),allocatable :: chain_sphi(:,:,:)
   real(dp),allocatable :: chain_resv(:,:,:)

 contains
   procedure :: free => rmm_diis_free
   procedure :: solve => rmm_diis_solve
   procedure :: eval_mats => rmm_diis_eval_mats

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

subroutine rmm_diis(cg, dtset, eig, occ, enlx, gs_hamk, gsc_all, mpi_enreg, nband, npw, nspinor, resid)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(in) :: mpi_enreg
 integer,intent(in) :: nband, npw, nspinor
 real(dp),intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),target,intent(inout) :: gsc_all(2,npw*nspinor*nband*dtset%usepaw)
 real(dp),intent(inout) :: enlx(nband)
 real(dp),intent(inout) :: resid(nband)
 real(dp),intent(in) :: occ(nband)
 real(dp),intent(out) :: eig(nband)

!Local variables-------------------------------
 integer,parameter :: ndat1 = 1, type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0, level = 432, use_subovl0 = 0
 integer :: ii, ig, ib, cplex, nline, ierr, prtvol
 integer :: iband, cpopt, sij_opt, is, ie, mcg, mgsc, istwf_k, optekin, usepaw, iline
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft, nbocc
 real(dp),parameter :: rdummy = zero
 real(dp) :: dotr, doti, lambda, rval, accuracy_ene
 character(len=500) :: msg
!arrays
 integer :: max_nlines_band(nband)
 real(dp) :: tsec(2), hist_ene(0:dtset%nline), hist_resid(0:dtset%nline)
 real(dp),allocatable :: ghc(:,:), gsc(:,:), gvnlxc(:,:), pcon(:)
 real(dp),allocatable :: residvec(:,:), kres(:,:), phi_now(:,:)
 real(dp),allocatable :: evec(:,:), subham(:), h_ij(:,:,:)
 real(dp) :: subovl(use_subovl0)
 real(dp), ABI_CONTIGUOUS pointer :: gsc_ptr(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(rmm_diis_t) :: diis

! *************************************************************************

 usepaw = dtset%usepaw; istwf_k = gs_hamk%istwf_k; prtvol = dtset%prtvol
 mcg = npw * nspinor * nband; mgsc = npw * nspinor * nband * usepaw

 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft

 gsc_ptr => null()

 cpopt = -1
 if (usepaw == 1) then
   sij_opt = 1 ! Recompute S
   !cpopt = 0  ! <p_lmn|in> are computed and saved
   cpopt = -1  ! <p_lmn|in> (and derivatives) are computed here (and not saved)
 else
   sij_opt = 0
 end if

 ABI_MALLOC(ghc, (2, npw*nspinor))
 ABI_MALLOC(gsc, (2, npw*nspinor*usepaw))
 ABI_MALLOC(gvnlxc, (2, npw*nspinor))

 ! =======================================
 ! Apply H to input cg to compute <i|H|j>
 ! =======================================
 ABI_CALLOC(h_ij, (2, nband, nband))

 do iband=1,nband
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)

   ! By setting ieigen to iband, Fock contrib of this iband to the energy will be calculated
   !call fock_set_ieigen(gs_hamk%fockcommon, iband)
   call getghc(cpopt, cg(:,is:ie), cprj_dum, ghc, gsc_ptr, gs_hamk, gvnlxc, &
               rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

   ! Compute <i|H|iband> for i=1,nband
   call cg_zgemv("C", npw*nspinor, nband, cg, ghc, h_ij(:,:,iband))

   if (istwf_k /= 1) then
     h_ij(:,:,iband) = two * h_ij(:,:,iband)

     if (istwf_k == 2 .and. me_g0 == 1) then
       ! Gamma k-point and I have G=0. Remove double counting term.
       do ib=1,nband
         ig = 1 + npw * nspinor * (ib - 1)
         h_ij(1,ib,iband) = h_ij(1,ib,iband) - cg(1,ig) * ghc(1,1)
       end do
     end if
     ! Force real matrix.
     h_ij(2,:,iband) = zero
   end if

 end do

 ! Pack <i|H|j> to prepare call to subdiago.
 ABI_MALLOC(subham, (nband*(nband+1)))
 call pack_matrix(h_ij, subham, nband, 2)
 ABI_FREE(h_ij)
 call xmpi_sum(subham, comm_spinorfft, ierr)

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

 ! =================
 ! Prepare DIIS loop
 ! =================
 optekin = 0; if (dtset%wfoptalg >= 10) optekin = 1
 optekin = 1
 !write(std_out,*)"optekin:", optekin
 nline = dtset%nline

 cplex = 2
 diis = rmm_diis_new(usepaw, istwf_k, npw, nspinor, ndat1, dtset%nline)

 ABI_MALLOC(pcon, (npw))
 ABI_MALLOC(residvec, (2, npw*nspinor))
 ABI_MALLOC(kres, (2, npw*nspinor))
 ABI_MALLOC(phi_now, (2, npw*nspinor))
 max_nlines_band = dtset%nline
 ! TODO: Don't reduce nline if MD
 !if dtset%
 do iband=1,nband
   if (occ(iband) < tol3) max_nlines_band(iband) = max(dtset%nline / 2, 2)
 end do
 nbocc = nint(dtset%nelect / two)
 !accuracy_ene = dtset%toldfe
 accuracy_ene = tol10 / nbocc / four

 iband_loop: do iband=1,nband
   ! Compute H |phi_0> from subdiago cg.
   is = 1 + npw * nspinor * (iband - 1); ie = is - 1 + npw * nspinor
   phi_now = cg(:,is:ie)
   if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)

   call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                        eig(iband), resid(iband), enlx(iband), residvec, normalize=.False.)
   hist_ene(0) = eig(iband)
   diis%chain_resv(:,:,0) = residvec

   ! Store |phi_0>, |S phi_0>.
   diis%chain_phi(:,:,0) = phi_now
   if (usepaw == 1) diis%chain_sphi(:,:,0) = gsc_ptr

   ! Compute <R0|R0> and <phi_0|S|phi_0>
   ! Assume input cg are already S-normalized.
   call sqnorm_g(resid(iband), istwf_k, npw*nspinor, residvec, me_g0, comm_fft)
   hist_resid(0) = resid(iband)
   diis%resmat(:, 0, 0) = [resid(iband), zero]
   diis%smat(:, 0, 0) = [one, zero]

   ! Band locking for NSCF or SCF if tolwfr > 0 is used.
   !if (resid(iband) < dtset%tolwfr) then
   if (resid(iband) < dtset%tolwfr .or. sqrt(abs(resid(iband))) < accuracy_ene) then
      cycle iband_loop
   end if

   ! Precondition |R_0>, output in kres = |K R_0>
   kres = residvec
   !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
   call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                  optekin, pcon, kres, comm_bandspinorfft)
   !write(std_out,*)"kres(1:2):", kres(:, 1:2)

   if (prtvol == -level) then
     write(msg,'(2(a5),2(a14),a12)')'iline', "band", "eigen_eV", "eigde_meV", "resid"
     call wrtout(std_out, msg, "PERS")
     write(msg,'(2(i5),2(f14.6,1x),es12.4)')0, iband, eig(iband) * Ha_eV, zero, resid(iband)
     call wrtout(std_out, msg, "PERS")
   end if

   iline_loop: do iline=1,max_nlines_band(iband)

     if (iline == 1) then
       ! Compute H |K R_0>
       call getghc(cpopt, kres, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
                   rdummy, mpi_enreg, ndat1, prtvol, sij_opt, tim_getghc, type_calc0)

       ! Compute residual: (H - e_0 S) |K R_0>
       if (usepaw == 1) then
         residvec = ghc - eig(iband) * gsc
       else
         residvec = ghc - eig(iband) * kres
       end if

       ! phi_1 = phi_0 + lambda K.R_0
       ! Compute lambda that minimizes the norm of the residual.
       ! lambda = - Re{<R_0|(H - e_0 S)} |K R_0>} / |(H - e_0 S) |K R_0>|**2
       call dotprod_g(dotr, doti, istwf_k, npw*nspinor, option1, &
                      diis%chain_resv(:,:,0), residvec, me_g0, comm_spinorfft)
       call sqnorm_g(rval, istwf_k, npw*nspinor, residvec, me_g0, comm_spinorfft)

       lambda = -dotr / rval
       !write(std_out, *)"lambda, dotr, rval ", lambda, dotr, rval
       phi_now = diis%chain_phi(:,:,0) + lambda * kres

     else
       ! Solve DIIS equations to get phi_now for iline >= 2
       call diis%solve(iline, npw*nspinor, phi_now, residvec)

       ! Precondition residual, output in kres.
       kres = residvec
       !call cg_precon(residvec, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
       call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                      optekin, pcon, kres, comm_bandspinorfft)

       ! Compute phi_now with lambda obtained in step #1
       phi_now = phi_now + lambda * kres
     end if

     ! Compute H |phi_now>.
     if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)

     call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                          eig(iband:), resid(iband:), enlx(iband), residvec, normalize=.True.)

     ! Store residual R_i for i == iline and evaluate new enlx for NC.
     hist_ene(iline) = eig(iband)
     hist_resid(iline) = resid(iband)
     diis%chain_phi(:,:,iline) = phi_now
     diis%chain_resv(:, :, iline) = residvec
     if (usepaw == 1) diis%chain_sphi(:,:,iline) = gsc_ptr

     ! ============== CHECK FOR CONVERGENCE ========================
     if (prtvol == -level) then
       write(msg,"(2(i5),2(f14.6,1x),es12.4)") &
         iline, iband, hist_ene(iline) * Ha_eV, (hist_ene(iline) - hist_ene(iline-1)) * Ha_meV, resid(iband)
       call wrtout(std_out, msg, 'PERS')
     end if

     if (exit_diis(iline, nline, hist_ene, hist_resid, accuracy_ene, dtset)) then
       if (prtvol == -level) then
         write(msg, '(2a,i4,a,i2,a,es12.4,a)' ) ch10, &
          ' band: ',iband,' converged after: ',iline,' iterations with resid: ',resid(iband), ch10
         call wrtout(std_out, msg , "PERS")
       end if
       exit iline_loop
     end if

     ! Compute <R_i|R_j> and <i|S|j> for j = iline
     if (iline /= nline) call diis%eval_mats(iline, me_g0, comm_spinorfft)

   end do iline_loop

   ! Update wavefunction for this band
   ! TODO: End with trial step but then I have to move the check for exit at the end of the loop.
   cg(:,is:ie) = phi_now

   if (.False. .and. iline == nline + 1) then
     !write(std_out, *)" Performing Last trial step"
     kres = diis%chain_resv(:,:, nline)
     call cg_precon(phi_now, zero, istwf_k, gs_hamk%kinpw_k, npw, nspinor, me_g0, &
                    optekin, pcon, kres, comm_bandspinorfft)
     !phi_now = phi_now + 0.1 * kres
     phi_now = phi_now + lambda * kres
     if (usepaw == 1) gsc_ptr => gsc_all(:,is:ie)

     call getghc_eigresid(gs_hamk, npw, nspinor, ndat1, phi_now, ghc, gsc_ptr, mpi_enreg, prtvol, &
                          eig(iband:), resid(iband:), enlx(iband), residvec, normalize=.True.)
     cg(:,is:ie) = phi_now
   end if
 end do iband_loop

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

logical function exit_diis(iline, nline, hist_ene, hist_resid, accuracy_ene, dtset) result(ans)
  integer,intent(in) :: iline, nline
  real(dp),intent(in) :: accuracy_ene
  real(dp),intent(in) :: hist_ene(0:nline), hist_resid(0:nline)
  type(dataset_type),intent(in) :: dtset

  real(dp) :: resid, deltae, deold

  resid = hist_resid(iline)
  deold = hist_ene(1) - hist_ene(0)
  deltae = hist_ene(iline) - hist_ene(iline-1)
  ans = .False.

  if (abs(deltae) < dtset%tolrde * abs(deold) .and. iline /= nline) ans = .True.

  if (dtset%iscf < 0) then
    ! This is the only condition available for NSCF run.
    if (resid < dtset%tolwfr) ans = .True.
  else
    if (sqrt(abs(resid)) < accuracy_ene) ans = .True.  ! Similar to EDIFF/NBANDS/4
    ! band locking (occupied states).
    !if (resid < 1.0d0-20) ans = .True.
    !if (resid < dtset%tolwfr) ans = .True.
  end if

end function exit_diis
!!***

subroutine getghc_eigresid(gs_hamk, npw, nspinor, ndat, phi_now, ghc, gsc, mpi_enreg, prtvol, &
                           eig, resid, enlx, residvec, normalize)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 integer,intent(in) :: npw, nspinor, ndat, prtvol
 real(dp),intent(inout) :: phi_now(2, npw*nspinor*ndat)
 real(dp),intent(out) :: ghc(2,npw*nspinor*ndat), gsc(2,npw*nspinor*ndat*gs_hamk%usepaw)
 type(mpi_type),intent(in) :: mpi_enreg
 real(dp),intent(out) :: eig(ndat), resid(ndat)
 real(dp),intent(out) :: enlx(ndat)
 real(dp),intent(out) :: residvec(2, npw*nspinor*ndat)
 logical,optional,intent(in) :: normalize

!Local variables-------------------------------
 integer,parameter :: type_calc0 = 0, option1 = 1, option2 = 2, tim_getghc = 0
 integer :: istwf_k, usepaw, cpopt, sij_opt, idat, is, ie, npws
 integer :: me_g0, comm_spinorfft, comm_bandspinorfft, comm_fft
 real(dp),parameter :: rdummy = zero
 logical :: normalize_
 real(dp) :: rval, dotr, doti
!arrays
 real(dp),allocatable :: gvnlxc(:, :)
 type(pawcprj_type) :: cprj_dum(1,1)

! *************************************************************************

 usepaw = gs_hamk%usepaw; istwf_k = gs_hamk%istwf_k !; prtvol = dtset%prtvol
 me_g0 = mpi_enreg%me_g0; comm_fft = mpi_enreg%comm_fft
 comm_spinorfft = mpi_enreg%comm_spinorfft; comm_bandspinorfft = mpi_enreg%comm_bandspinorfft
 normalize_ = .True.; if (present(normalize)) normalize_ = normalize
 npws = npw * nspinor

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
   call cgnc_normalize(npw * nspinor, ndat, phi_now, istwf_k, me_g0, comm_spinorfft)
 end if

 ABI_MALLOC(gvnlxc, (2, npw*nspinor*ndat))

 ! Compute H |phi_now>
 !call fock_set_ieigen(gs_hamk%fockcommon, iband)
 call getghc(cpopt, phi_now, cprj_dum, ghc, gsc, gs_hamk, gvnlxc, &
             rdummy, mpi_enreg, ndat, prtvol, sij_opt, tim_getghc, type_calc0)

 if (usepaw == 1 .and. normalize_) then
   ! PAW normalization must be done here.
   call cgpaw_normalize(npw * nspinor, ndat, phi_now, gsc, istwf_k, me_g0, comm_spinorfft)
 end if

 do idat=1,ndat
   is = 1 + (idat - 1) * npw * nspinor
   ie = is - 1 + npw * nspinor
   call dotprod_g(eig(idat), doti, istwf_k, npw*nspinor, option1, ghc(:,is), phi_now(:,is), me_g0, comm_spinorfft)

   if (usepaw == 1) then
     call dotprod_g(dotr, doti, istwf_k, npw*nspinor, option1, gsc(:,is), phi_now(:,is), me_g0, comm_spinorfft)
     eig(idat) = eig(idat) / dotr
   end if

   ! Residual.
   if (usepaw == 1) then
     residvec(:,is:ie) = ghc(:,is:ie) - eig(idat) * gsc(:,is:ie)
   else
     residvec(:,is:ie) = ghc(:,is:ie) - eig(idat) * phi_now(:,is:ie)
   end if

   ! Store residual R_i for i == iline and evaluate new enlx for NC.
   call sqnorm_g(resid(idat), istwf_k, npw*nspinor, residvec(1,is), me_g0, comm_fft)

   if (usepaw == 0) then
     call dotprod_g(enlx(idat), doti, istwf_k, npw*nspinor, option1, &
                    phi_now(1,is), gvnlxc(1,is), me_g0, comm_bandspinorfft)
   end if
 end do

 ABI_FREE(gvnlxc)

end subroutine getghc_eigresid
!!!***

type(rmm_diis_t) function rmm_diis_new(usepaw, istwf_k, npw, nspinor, ndat, nline) result(diis)
  integer,intent(in) :: usepaw, istwf_k, npw, nspinor, ndat, nline

 diis%usepaw = usepaw
 diis%ndat = ndat
 diis%istwf_k = istwf_k
 diis%cplex = 2
 diis%nline = nline
 diis%npw = npw
 diis%nspinor = nspinor

 ABI_MALLOC(diis%chain_phi, (2, npw*nspinor, 0:nline))
 ABI_MALLOC(diis%chain_sphi, (2, npw*nspinor*usepaw, 0:nline))
 ABI_MALLOC(diis%chain_resv, (2, npw*nspinor, 0:nline))
 ABI_MALLOC(diis%resmat, (2, 0:nline, 0:nline)) ! <R_i|R_j>
 ABI_MALLOC(diis%smat, (2, 0:nline, 0:nline))   ! <i|S|j>

end function rmm_diis_new

subroutine rmm_diis_free(diis)
  class(rmm_diis_t),intent(inout) :: diis

  ABI_FREE(diis%chain_phi)
  ABI_FREE(diis%chain_sphi)
  ABI_FREE(diis%chain_resv)
  ABI_FREE(diis%resmat)
  ABI_FREE(diis%smat)

end subroutine rmm_diis_free

subroutine rmm_diis_solve(diis, iline, npws, phi_now, residvec)
 class(rmm_diis_t),intent(in) :: diis
 integer,intent(in) :: iline, npws
 real(dp),intent(inout) :: phi_now(2, npws), residvec(2, npws)

 integer :: cplex, ierr, lwork, npw, nspinor
 integer,allocatable :: ipiv(:)
 real(dp),allocatable :: diis_eig(:)
 real(dp),allocatable :: wmat1(:,:,:), wmat2(:,:,:), wvec(:,:)
 complex(dp),allocatable :: work(:)
 character(len=500) :: msg

 cplex = diis%cplex; npw = diis%npw; nspinor = diis%nspinor

#if 0
   ABI_MALLOC(diis_eig, (0:iline-1))
   ABI_MALLOC(wmat1, (cplex, 0:iline-1, 0:iline-1))
   ABI_MALLOC(wmat2, (cplex, 0:iline-1, 0:iline-1))
   wmat1 = diis%resmat(:, 0:iline-1, 0:iline-1)
   wmat2 = diis%smat(:, 0:iline-1, 0:iline-1)
   !do ii=0,iline-1; write(std_out, *) "diis_resmat:", wmat1(:,ii,:); end do
   !do ii=0,iline-1; write(std_out, *) "diis_smat:", wmat2(:,ii,:); end do

   call xhegv_cplex(1, "V", "U", cplex, iline, wmat1, wmat2, diis_eig, msg, ierr)
   !write(std_out,*)"diis_eig:", diis_eig(0)
   !write(std_out,*)"RE diis_vec  :", wmat1(1,:,0)
   !write(std_out,*)"IMAG diis_vec:", wmat1(2,:,0)
   ABI_CHECK(ierr == 0, "xhegv returned ierr != 0")
   if (ierr /= 0) then
     call wrtout(std_out, sjoin("xhegv failed with:", msg, ch10, "at iline: ", itoa(iline), "exit iline_loop!"))
     ABI_FREE(wmat1)
     ABI_FREE(wmat2)
     ABI_FREE(diis_eig)
     !exit iline_loop
   end if

   ! Take linear combination of chain_phi and chain_resv.
   call cg_zgemv("N", npw*nspinor, iline, diis%chain_phi, wmat1(:,:,0), phi_now)
   call cg_zgemv("N", npw*nspinor, iline, diis%chain_resv, wmat1(:,:,0), residvec)
   !call sqnorm_g(rval, istwf_k, npw*nspinor, phi_now, me_g0, comm_spinorfft)
   !phi_now = phi_now / sqrt(rval)

   ABI_FREE(wmat1)
   ABI_FREE(wmat2)
   ABI_FREE(diis_eig)

#else
   ! Solve system of linear equations.
   ABI_CALLOC(wvec, (cplex, 0:iline))
   wvec(1, iline) = -one
   ABI_CALLOC(wmat1, (cplex, 0:iline, 0:iline))
   wmat1(1,:,iline) = -one
   wmat1(1,iline,:) = -one
   wmat1(1,iline,iline) = zero
   wmat1(:,0:iline-1, 0:iline-1) = diis%resmat(:, 0:iline-1, 0:iline-1)

   ABI_MALLOC(ipiv, (0:iline))
   ABI_MALLOC(work, (1))
   lwork = -1
   call zhesv("U", iline+1, 1, wmat1, iline+1, ipiv, wvec, iline+1, work, lwork, ierr)
   lwork = int(work(1))
   ABI_REMALLOC(work, (lwork))
   call zhesv("U", iline+1, 1, wmat1, iline+1, ipiv, wvec, iline+1, work, lwork, ierr)
   ABI_CHECK(ierr == 0, "zhesv returned ierr != 0")
   ABI_FREE(ipiv)
   ABI_FREE(work)
   ABI_FREE(wmat1)

   call cg_zgemv("N", npw*nspinor, iline, diis%chain_phi, wvec(:,0), phi_now)
   call cg_zgemv("N", npw*nspinor, iline, diis%chain_resv, wvec(:,0), residvec)
   !lambda = one
   ABI_FREE(wvec)
#endif

end subroutine rmm_diis_solve

subroutine rmm_diis_eval_mats(diis, iline, me_g0, comm_spinorfft)

 class(rmm_diis_t),intent(inout) :: diis
 integer,intent(in) :: iline, me_g0, comm_spinorfft

 integer,parameter :: option2 = 2
 integer :: ii, npw, nspinor, istwf_k
 real(dp) :: dotr, doti

 npw = diis%npw; nspinor = diis%nspinor; istwf_k = diis%istwf_k

 do ii=0,iline
   ! <R_i|R_j>
   call dotprod_g(dotr, doti, istwf_k, npw*nspinor, option2, &
                  diis%chain_resv(:,:,ii), diis%chain_resv(:,:,iline), me_g0, comm_spinorfft)
   if (ii == iline) doti = zero
   diis%resmat(:, ii, iline) = [dotr, doti]

   ! <i|S|j>
   if (ii == iline) then
     dotr = one; doti = zero
   else
   if (diis%usepaw == 0) then
     call dotprod_g(dotr, doti, istwf_k, npw*nspinor, option2, &
                    diis%chain_phi(:,:,ii), diis%chain_phi(:,:,iline), me_g0, comm_spinorfft)
   else
     call dotprod_g(dotr, doti, istwf_k, npw*nspinor, option2, &
                    diis%chain_phi(:,:,ii), diis%chain_sphi(:,:,iline), me_g0, comm_spinorfft)
   end if
   if (ii == iline) doti = zero
   end if
   diis%smat(:, ii, iline) = [dotr, doti]
 end do

end subroutine rmm_diis_eval_mats
!!***

end module m_rmm_diis
!!***
