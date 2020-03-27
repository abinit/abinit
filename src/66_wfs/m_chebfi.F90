!!****m* ABINIT/m_chebfi
!! NAME
!!  m_chebfi
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (AL)
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

module m_chebfi

 use defs_basis
 use m_errors
 use m_xmpi
 use m_abicore
 use m_abi_linalg
 use m_rayleigh_ritz
 use m_invovl
 use m_dtset

 use defs_abitypes, only : mpi_type
 use m_time,          only : timab
 use m_cgtools,       only : dotprod_g
 use m_bandfft_kpt,   only : bandfft_kpt, bandfft_kpt_get_ikpt
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_axpby, pawcprj_copy
 use m_hamiltonian,   only : gs_hamiltonian_type
 use m_getghc,        only : getghc
 use m_prep_kgb,      only : prep_getghc, prep_index_wavef_bandpp

 implicit none

 private
!!***

 public :: chebfi
!!***

contains
!!***


!!****f* ABINIT/chebfi
!! NAME
!! chebfi
!!
!! FUNCTION
!! this routine updates the wave functions at a given k-point,
!! using the ChebFi method (see paper by A. Levitt and M. Torrent)
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands at this k point for that spin polarization
!!  npw=number of plane waves at this k point
!!  nspinor=number of plane waves at this k point
!!  prtvol=control print volume and debugging output
!!
!! OUTPUT
!!  eig(nband)=array for holding eigenvalues (hartree)
!!  resid(nband)=residuals for each states
!!  If gs_hamk%usepaw==1:
!!    gsc(2,*)=<g|s|c> matrix elements (s=overlap)
!!  If gs_hamk%usepaw==0
!!    enlx(nband)=contribution from each band to nonlocal psp + potential Fock ACE part of total energy, at this k-point
!!
!! SIDE EFFECTS
!!  cg(2,*)=updated wavefunctions
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      apply_invovl,dotprod_g,getghc,pawcprj_alloc,pawcprj_axpby,pawcprj_copy
!!      pawcprj_free,prep_getghc,prep_index_wavef_bandpp
!!      rayleigh_ritz_distributed,rayleigh_ritz_subdiago,timab,wrtout
!!      xmpi_alltoallv,xmpi_barrier,xmpi_max,xmpi_min,xmpi_sum
!!
!! NOTES
!!  -- TODO --
!!  Normev?
!!  Ecutsm
!!  nspinor 2
!!  spinors parallelisation
!!  fock
!!  -- Performance --
!!  Improve load balancing
!!  Don't diagonalize converged eigenvectors, just orthogonalize
!!  Maybe don't diagonalize so often (once every two outer iterations?)
!!  Benchmark diagonalizations, choose np_slk
!!  How to chose npfft?
!!  Implement MINRES for invovl
!!  -- LOBPCG --
!!  Improve stability (see paper by Lehoucq Sorensen, maybe use bunch-kaufman factorizations?)
!!
!! SOURCE

subroutine chebfi(cg,dtset,eig,enlx,gs_hamk,gsc,kinpw,mpi_enreg,nband,npw,nspinor,prtvol,resid)

!Arguments ------------------------------------
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(dataset_type),intent(in) :: dtset
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: nband,npw,prtvol,nspinor
 real(dp),intent(inout), target :: cg(2,npw*nspinor*nband),gsc(2,npw*nspinor*nband)
 real(dp),intent(in) :: kinpw(npw)
 real(dp),intent(out) :: resid(nband)
 real(dp),intent(out) :: enlx(nband)
 real(dp),intent(out) :: eig(nband)

!Local variables-------------------------------
 real(dp) :: pcon(npw)
 real(dp) :: filter_low
 real(dp) :: filter_center, filter_radius
 real(dp), dimension(2, npw*nspinor*nband), target :: ghc, gvnlxc
 real(dp), allocatable, dimension(:,:) :: cg_filter_next, cg_filter_prev, gsm1hc_filter, gsc_filter_prev, gsc_filter_next
 real(dp), allocatable, dimension(:,:), target :: cg_alltoall1,gsc_alltoall1,ghc_alltoall1,gvnlxc_alltoall1
 real(dp), allocatable, dimension(:,:), target :: cg_alltoall2,gsc_alltoall2,ghc_alltoall2,gvnlxc_alltoall2
 real(dp), pointer, dimension(:,:) :: cg_filter, gsc_filter, ghc_filter, gvnlxc_filter
 real(dp) :: resid_vec(2, npw*nspinor)
 logical :: has_fock,paw
 integer :: shift, shift_cg_loadbalanced
 integer :: iband, iline, ispinor
 integer :: sij_opt, cpopt
 real(dp) :: eval, tsec(2)
 integer :: tim_getghc = 5, ierr
 integer :: i
 integer, allocatable :: index_wavef_band(:)
 real(dp) :: maxeig, mineig
 real(dp), allocatable :: resids_filter(:), residvec_filter(:,:)
 integer, allocatable :: nline_bands(:)
 integer :: iactive, nactive
 real(dp) :: ampfactor
 integer :: nline_max, nline_decrease, nline_tolwfr
 ! real(dp) :: load_imbalance
 integer :: mcg
 real(dp) :: dprod_r, dprod_i
 character(len=500) :: message
 integer :: rdisplsloc(mpi_enreg%nproc_band),recvcountsloc(mpi_enreg%nproc_band)
 integer :: sdisplsloc(mpi_enreg%nproc_band), sendcountsloc(mpi_enreg%nproc_band)
 integer :: ikpt_this_proc, npw_filter, nband_filter
 type(pawcprj_type), allocatable :: cwaveprj(:,:), cwaveprj_next(:,:), cwaveprj_prev(:,:)
 ! integer :: nline_total

 ! timers
 integer, parameter :: timer_chebfi = 1600, timer_alltoall = 1601, timer_apply_inv_ovl = 1602, timer_rotation = 1603
 integer, parameter :: timer_subdiago = 1604, timer_subham = 1605, timer_ortho = 1606, timer_getghc = 1607
 integer, parameter :: timer_residuals = 1608, timer_update_eigen = 1609, timer_sync = 1610

! *************************************************************************

 !======================================================================================================
 ! Initialize, transpose input cg if paral_kgb
 !======================================================================================================
 call timab(timer_chebfi,1,tsec)

 !Initializations
 paw = gs_hamk%usepaw == 1
 has_fock=(associated(gs_hamk%fockcommon))
 mcg = npw*nspinor*nband

 ! Init pcon
 pcon = (27+kinpw*(18+kinpw*(12+8*kinpw))) / (27+kinpw*(18+kinpw*(12+8*kinpw)) + 16*kinpw**4)

 ghc=zero; gvnlxc=zero

 ! Initialize the _filter pointers. Depending on paral_kgb, they might point to the actual arrays or to _alltoall variables
 if (dtset%paral_kgb == 1) then
   ikpt_this_proc = bandfft_kpt_get_ikpt()
   npw_filter = bandfft_kpt(ikpt_this_proc)%ndatarecv
   nband_filter = mpi_enreg%bandpp

   ABI_ALLOCATE(cg_alltoall1, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(gsc_alltoall1, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(ghc_alltoall1, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(gvnlxc_alltoall1, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(cg_alltoall2, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(gsc_alltoall2, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(ghc_alltoall2, (2, npw_filter*nspinor*nband_filter))
   ABI_ALLOCATE(gvnlxc_alltoall2, (2, npw_filter*nspinor*nband_filter))

   ! Init tranpose variables
   recvcountsloc=bandfft_kpt(ikpt_this_proc)%recvcounts*2*nspinor*mpi_enreg%bandpp
   rdisplsloc=bandfft_kpt(ikpt_this_proc)%rdispls*2*nspinor*mpi_enreg%bandpp
   sendcountsloc=bandfft_kpt(ikpt_this_proc)%sendcounts*2*nspinor
   sdisplsloc=bandfft_kpt(ikpt_this_proc)%sdispls*2*nspinor

   ! Load balancing, so that each processor has approximately the same number of converged and non-converged bands
   ! for two procs, rearrange 1 2 3 4 5 6 as 1 4 2 5 3 6
   !
   ! trick to save memory: ghc has the necessary size, and will be overwritten afterwards anyway
#define cg_loadbalanced ghc
   shift = 0
   do i=1, mpi_enreg%nproc_band
     do iband=1, mpi_enreg%bandpp
       shift_cg_loadbalanced = (i-1 + (iband-1)*mpi_enreg%nproc_band)*npw*nspinor
       cg_loadbalanced(:, shift+1:shift+npw*nspinor) = cg(:, shift_cg_loadbalanced+1:shift_cg_loadbalanced+npw*nspinor)
       shift = shift + npw*nspinor
     end do
   end do

   ! Transpose input cg into cg_alloall1. cg_alltoall1 is now (npw_filter, nband_filter)
   call timab(timer_alltoall, 1, tsec)
   call xmpi_alltoallv(cg_loadbalanced,sendcountsloc,sdisplsloc,cg_alltoall1,&
&   recvcountsloc,rdisplsloc,mpi_enreg%comm_band,ierr)
   call timab(timer_alltoall, 2, tsec)
#undef cg_loadbalanced

   ! sort according to bandpp (from lobpcg, I don't fully understand what's going on but it works and it's fast)
   call prep_index_wavef_bandpp(mpi_enreg%nproc_band,mpi_enreg%bandpp,&
&   nspinor,bandfft_kpt(ikpt_this_proc)%ndatarecv,&
&   bandfft_kpt(ikpt_this_proc)%recvcounts,bandfft_kpt(ikpt_this_proc)%rdispls,&
&   index_wavef_band)
   cg_alltoall2(:,:) = cg_alltoall1(:,index_wavef_band)

   cg_filter => cg_alltoall2
   gsc_filter => gsc_alltoall2
   ghc_filter => ghc_alltoall2
   gvnlxc_filter => gvnlxc_alltoall2
 else
   npw_filter = npw
   nband_filter = nband

   cg_filter => cg
   gsc_filter => gsc
   ghc_filter => ghc
   gvnlxc_filter => gvnlxc
 end if
 ! from here to the next alltoall, all computation is done on _filter variables, agnostic
 ! to whether it's nband x npw (paral_kgb == 0) or ndatarecv*bandpp (paral_kgb = 1)

 ! Allocate filter variables for the application of the Chebyshev polynomial
 ABI_ALLOCATE(cg_filter_next, (2, npw_filter*nspinor*nband_filter))
 ABI_ALLOCATE(cg_filter_prev, (2, npw_filter*nspinor*nband_filter))
 ABI_ALLOCATE(gsc_filter_prev, (2, npw_filter*nspinor*nband_filter))
 ABI_ALLOCATE(gsc_filter_next, (2, npw_filter*nspinor*nband_filter))
 ABI_ALLOCATE(gsm1hc_filter, (2, npw_filter*nspinor*nband_filter))

 ! PAW init
 if(paw) then
   ABI_DATATYPE_ALLOCATE(cwaveprj, (gs_hamk%natom,nspinor*nband_filter))
   ABI_DATATYPE_ALLOCATE(cwaveprj_next, (gs_hamk%natom,nspinor*nband_filter))
   ABI_DATATYPE_ALLOCATE(cwaveprj_prev, (gs_hamk%natom,nspinor*nband_filter))
   call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
   call pawcprj_alloc(cwaveprj_next,0,gs_hamk%dimcprj)
   call pawcprj_alloc(cwaveprj_prev,0,gs_hamk%dimcprj)

   sij_opt = 1 ! recompute S
   cpopt = 0 ! save cprojs
 else
   sij_opt = 0
   cpopt = -1
 end if



 !======================================================================================================
 ! Data in npfft x npband distribution. First getghc, update eigenvalues and residuals
 !======================================================================================================
 write(message, *) 'First getghc'
 call wrtout(std_out,message,'COLL')

 ! get_ghc on cg
 call timab(timer_getghc, 1, tsec)
 if (dtset%paral_kgb == 0) then
   call getghc(cpopt,cg_filter,cwaveprj,ghc_filter,gsc_filter,gs_hamk,gvnlxc_filter,&
&   eval,mpi_enreg,nband,prtvol,sij_opt,tim_getghc,0)
 else
   call prep_getghc(cg_filter,gs_hamk,gvnlxc_filter,ghc_filter,gsc_filter,eval,nband,mpi_enreg,&
&   prtvol,sij_opt,cpopt,cwaveprj,already_transposed=.true.)
 end if
 call timab(timer_getghc, 2, tsec)

 ! Debug barrier: should be invisible
 call timab(timer_sync, 1, tsec)
 call xmpi_barrier(mpi_enreg%comm_band)
 call timab(timer_sync, 2, tsec)

 write(message, *) 'Computing residuals'
 call wrtout(std_out,message,'COLL')
 ! update eigenvalues and residuals
 call timab(timer_update_eigen, 1, tsec)
 ABI_ALLOCATE(resids_filter, (nband_filter))
 ABI_ALLOCATE(residvec_filter, (2, npw_filter*nspinor))
 ABI_ALLOCATE(nline_bands, (nband_filter))
 do iband=1, nband_filter
   shift = npw_filter*nspinor*(iband-1)
   call dotprod_g(eig(iband),dprod_i,gs_hamk%istwf_k,npw_filter*nspinor,1,ghc_filter(:, shift+1:shift+npw_filter*nspinor),&
&   cg_filter(:, shift+1:shift+npw_filter*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   if(paw) then
     call dotprod_g(dprod_r,dprod_i,gs_hamk%istwf_k,npw_filter*nspinor,1,gsc_filter(:, shift+1:shift+npw_filter*nspinor),&
&     cg_filter(:, shift+1:shift+npw_filter*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     eig(iband) = eig(iband)/dprod_r
   end if

   if(paw) then
     residvec_filter = ghc_filter(:, shift+1 : shift+npw_filter*nspinor) &
&     - eig(iband)*gsc_filter(:, shift+1 : shift+npw_filter*nspinor)
   else
     residvec_filter = ghc_filter(:, shift+1 : shift+npw_filter*nspinor) &
&     - eig(iband)*cg_filter(:, shift+1 : shift+npw_filter*nspinor)
   end if
   resids_filter(iband) = SUM(residvec_filter**2)
 end do
 call xmpi_sum(resids_filter,mpi_enreg%comm_fft,ierr)
 call xmpi_max(MAXVAL(eig(1:nband_filter)),maxeig,mpi_enreg%comm_band,ierr)
 call xmpi_min(MINVAL(eig(1:nband_filter)),mineig,mpi_enreg%comm_band,ierr)
 filter_low = maxeig
 call timab(timer_update_eigen, 2, tsec)

 ! Decide how many iterations per band are needed
 ! don't go above this, or face bad conditioning of the Gram matrix.
 nline_max = cheb_oracle(mineig, filter_low, dtset%ecut, 1e-16_dp, 40)
 ! if(mpi_enreg%me == 0) write(0, *) nline_max
 do iband=1, nband_filter
   ! nline necessary to converge to tolwfr
   nline_tolwfr = cheb_oracle(eig(iband), filter_low, dtset%ecut, dtset%tolwfr / resids_filter(iband), dtset%nline)
   ! nline necessary to decrease residual by a constant factor
   nline_decrease = cheb_oracle(eig(iband), filter_low, dtset%ecut, 0.1_dp, dtset%nline)

   nline_bands(iband) = MAX(MIN(nline_tolwfr, nline_decrease, nline_max, dtset%nline), 1)
   nline_bands(iband) = dtset%nline ! fiddle with this to use locking
 end do


 !!!!! Uncomment for diagnostics
 ! nline_total = SUM(nline_bands)
 ! call xmpi_sum(nline_total, mpi_enreg%comm_band, ierr)
 ! load_imbalance = (SUM(nline_bands) - REAL(nline_total)/REAL(mpi_enreg%nproc_band)) / &
 ! &                (REAL(nline_total)/REAL(mpi_enreg%nproc_band))
 ! call xmax_mpi(load_imbalance, mpi_enreg%comm_band, ierr)

 ! write(message, *) 'Mean nline', REAL(nline_total)/REAL(nband), 'max imbalance (%)', load_imbalance*100
 ! call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(resids_filter)
 ABI_DEALLOCATE(residvec_filter)

 !======================================================================================================
 ! Chebyshev polynomial application
 !======================================================================================================
 ! Filter by a chebyshev polynomial of degree nline
 do iline=1,dtset%nline
   ! Filter only on [iactive, iactive+nactive-1]
   iactive = nband_filter
   do iband = 1, nband_filter
     ! does iband need an iteration?
     if (nline_bands(iband) >= iline) then
       iactive = iband
       exit
     end if
   end do
   nactive = nband_filter - iactive + 1
   shift = npw_filter*nspinor*(iactive-1) + 1
   ! trick the legacy prep_getghc
   mpi_enreg%bandpp = nactive

   ! Define the filter position
   filter_center = (dtset%ecut+filter_low)/2
   filter_radius = (dtset%ecut-filter_low)/2

   ! write(message, *) 'Applying invovl, iteration', iline
   ! call wrtout(std_out,message,'COLL')

   ! If paw, have to apply S^-1
   if(paw) then
     call timab(timer_apply_inv_ovl, 1, tsec)
     call apply_invovl(gs_hamk, ghc_filter(:,shift:), gsm1hc_filter(:,shift:), cwaveprj_next(:,iactive:), &
&     npw_filter, nactive, mpi_enreg, nspinor)
     call timab(timer_apply_inv_ovl, 2, tsec)
   else
     gsm1hc_filter(:,shift:) = ghc_filter(:,shift:)
   end if

   ! Chebyshev iteration: UPDATE cg
   if(iline == 1) then
     cg_filter_next(:,shift:) = one/filter_radius * (gsm1hc_filter(:,shift:) - filter_center*cg_filter(:,shift:))
   else
     cg_filter_next(:,shift:) = two/filter_radius * (gsm1hc_filter(:,shift:) - filter_center*cg_filter(:,shift:)) &
&     - cg_filter_prev(:,shift:)
   end if
   ! Update gsc and cwaveprj
   if(paw) then
     if(iline == 1) then
       gsc_filter_next(:,shift:) = one/filter_radius * (ghc_filter(:,shift:) - filter_center*gsc_filter(:,shift:))
       !cwaveprj_next = one/filter_radius * (cwaveprj_next - filter_center*cwaveprj)
       call pawcprj_axpby(-filter_center/filter_radius, one/filter_radius,cwaveprj(:,iactive:),cwaveprj_next(:,iactive:))
     else
       gsc_filter_next(:,shift:) = two/filter_radius * (ghc_filter(:,shift:) - filter_center*gsc_filter(:,shift:))&
&       - gsc_filter_prev(:,shift:)
       !cwaveprj_next = two/filter_radius * (cwaveprj_next - filter_center*cwaveprj) - cwaveprj_prev
       call pawcprj_axpby(-two*filter_center/filter_radius, two/filter_radius,cwaveprj(:,iactive:),cwaveprj_next(:,iactive:))
       call pawcprj_axpby(-one, one,cwaveprj_prev(:,iactive:),cwaveprj_next(:,iactive:))
     end if
   end if

   ! Bookkeeping of the _prev variables
   cg_filter_prev(:,shift:) = cg_filter(:,shift:)
   cg_filter(:,shift:) = cg_filter_next(:,shift:)
   if(paw) then
     gsc_filter_prev(:,shift:) = gsc_filter(:,shift:)
     gsc_filter(:,shift:) = gsc_filter_next(:,shift:)

     call pawcprj_copy(cwaveprj(:,iactive:),cwaveprj_prev(:,iactive:))
     call pawcprj_copy(cwaveprj_next(:,iactive:),cwaveprj(:,iactive:))
   end if

   ! Update ghc
   if(paw) then
     !! DEBUG use this to remove the optimization and recompute gsc/cprojs
     ! sij_opt = 1
     ! cpopt = 0

     sij_opt = 0 ! gsc is already computed
     cpopt = 2 ! reuse cprojs
   else
     sij_opt = 0
     cpopt = -1
   end if

   write(message, *) 'Getghc, iteration', iline
   call wrtout(std_out,message,'COLL')

   call timab(timer_getghc, 1, tsec)
   if (dtset%paral_kgb == 0) then
     call getghc(cpopt,cg_filter(:,shift:),cwaveprj(:,iactive:),ghc_filter(:,shift:),&
&     gsc_filter(:,shift:),gs_hamk,gvnlxc_filter(:,shift:),eval,mpi_enreg,&
&     nband,prtvol,sij_opt,tim_getghc,0)
   else
     call prep_getghc(cg_filter(:,shift:),gs_hamk,gvnlxc_filter(:,shift:),ghc_filter(:,shift:),&
&     gsc_filter(:,shift:),eval,nband,mpi_enreg,prtvol,sij_opt,cpopt,&
&     cwaveprj(:,iactive:),already_transposed=.true.)
   end if

   ! end of the trick
   mpi_enreg%bandpp = nband_filter

   call timab(timer_getghc, 2, tsec)
 end do ! end loop on nline

 ! normalize according to the previously computed rayleigh quotients (inaccurate, but cheap)
 do iband = 1, nband_filter
   ampfactor = cheb_poly(eig(iband), nline_bands(iband), filter_low, dtset%ecut)
   if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 ! just in case, avoid amplifying too much
   shift = npw_filter*nspinor*(iband-1)
   cg_filter(:, shift+1:shift+npw_filter*nspinor) = cg_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
   ghc_filter(:, shift+1:shift+npw_filter*nspinor) = ghc_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
   if(paw) then
     gsc_filter(:, shift+1:shift+npw_filter*nspinor) = gsc_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
   endif
   if(.not.paw .or. has_fock)then
     gvnlxc_filter(:, shift+1:shift+npw_filter*nspinor) = gvnlxc_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
   end if
 end do

 ! Cleanup
 if(paw) then
   call pawcprj_free(cwaveprj)
   call pawcprj_free(cwaveprj_next)
   call pawcprj_free(cwaveprj_prev)
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
   ABI_DATATYPE_DEALLOCATE(cwaveprj_next)
   ABI_DATATYPE_DEALLOCATE(cwaveprj_prev)
 end if
 ABI_DEALLOCATE(nline_bands)
 ABI_DEALLOCATE(cg_filter_next)
 ABI_DEALLOCATE(cg_filter_prev)
 ABI_DEALLOCATE(gsc_filter_prev)
 ABI_DEALLOCATE(gsc_filter_next)
 ABI_DEALLOCATE(gsm1hc_filter)

 !======================================================================================================
 ! Filtering done, tranpose back
 !======================================================================================================

 write(message, *) 'Filtering done, transposing back'
 call wrtout(std_out,message,'COLL')

 ! transpose back
 if(dtset%paral_kgb == 1) then
   cg_alltoall1(:,index_wavef_band) = cg_alltoall2(:,:)
   ghc_alltoall1(:,index_wavef_band) = ghc_alltoall2(:,:)
   if(paw) then
     gsc_alltoall1(:,index_wavef_band) = gsc_alltoall2(:,:)
   else
     gvnlxc_alltoall1(:,index_wavef_band) = gvnlxc_alltoall2(:,:)
   end if

   ABI_DEALLOCATE(index_wavef_band)

   call timab(timer_sync, 1, tsec)
   call xmpi_barrier(mpi_enreg%comm_band)
   call timab(timer_sync, 2, tsec)

   call timab(timer_alltoall, 1, tsec)

  ! Do we pack the arrays in the alltoall, saving latency, or do we do it separately, saving memory and copies?
   call xmpi_alltoallv(cg_alltoall1,recvcountsloc,rdisplsloc,cg,&
&   sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   call xmpi_alltoallv(ghc_alltoall1,recvcountsloc,rdisplsloc,ghc,&
&   sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   if(paw) then
     call xmpi_alltoallv(gsc_alltoall1,recvcountsloc,rdisplsloc,gsc,&
&     sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   else
     call xmpi_alltoallv(gvnlxc_alltoall1,recvcountsloc,rdisplsloc,gvnlxc,&
&     sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   end if
   call timab(timer_alltoall, 2, tsec)

   if(mpi_enreg%paral_kgb == 1) then
     ABI_DEALLOCATE(cg_alltoall1)
     ABI_DEALLOCATE(gsc_alltoall1)
     ABI_DEALLOCATE(ghc_alltoall1)
     ABI_DEALLOCATE(gvnlxc_alltoall1)
     ABI_DEALLOCATE(cg_alltoall2)
     ABI_DEALLOCATE(gsc_alltoall2)
     ABI_DEALLOCATE(ghc_alltoall2)
     ABI_DEALLOCATE(gvnlxc_alltoall2)
   end if
 else
   ! nothing to do, the _filter variables already point to the right ones
 end if



 !======================================================================================================
 ! Data in (npfft*npband) x 1 distribution. Rayleigh-Ritz step
 !======================================================================================================

 ! _subdiago might use less memory when using only one proc, should maybe call it, or just remove it
 ! and always call _distributed
#if defined HAVE_LINALG_SCALAPACK
 call rayleigh_ritz_distributed(cg,ghc,gsc,gvnlxc,eig,has_fock,gs_hamk%istwf_k,mpi_enreg,nband,npw,nspinor,gs_hamk%usepaw)
#else
 call rayleigh_ritz_subdiago(cg,ghc,gsc,gvnlxc,eig,has_fock,gs_hamk%istwf_k,mpi_enreg,nband,npw,nspinor,gs_hamk%usepaw)
#endif

 ! Build residuals
 call timab(timer_residuals, 1, tsec)
 do iband=1,nband
   shift = npw*nspinor*(iband-1)
   if(paw) then
     resid_vec = ghc(:, shift+1 : shift+npw*nspinor) - eig(iband)*gsc(:, shift+1 : shift+npw*nspinor)
   else
     resid_vec = ghc(:, shift+1 : shift+npw*nspinor) - eig(iband)*cg (:, shift+1 : shift+npw*nspinor)
   end if

   ! precondition resid_vec
   do ispinor = 1,nspinor
     resid_vec(1, npw*(ispinor-1)+1:npw*ispinor) = resid_vec(1, npw*(ispinor-1)+1:npw*ispinor) * pcon
     resid_vec(2, npw*(ispinor-1)+1:npw*ispinor) = resid_vec(2, npw*(ispinor-1)+1:npw*ispinor) * pcon
   end do

   call dotprod_g(resid(iband),dprod_i,gs_hamk%istwf_k,npw*nspinor,1,resid_vec,&
&   resid_vec,mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)

   if(.not. paw .or. has_fock) then
     call dotprod_g(enlx(iband),dprod_i,gs_hamk%istwf_k,npw*nspinor,1,cg(:, shift+1:shift+npw*nspinor),&
&     gvnlxc(:, shift+1:shift+npw_filter*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
   end if
 end do
 call timab(timer_residuals, 2, tsec)

 ! write(message, '(a,4e10.2)') 'Resids (1, N, min, max) ', resid(1), resid(nband), MINVAL(resid), MAXVAL(resid)
 ! call wrtout(std_out,message,'COLL')

 ! write(message,*)'Eigens(1,nocc,nband) ',eig(1), eig(ilastocc),eig(nband)
 ! call wrtout(std_out,message,'COLL')
 ! write(message,*)'Resids(1,nocc,nband) ',resid(1), resid(ilastocc),resid(nband)
 ! call wrtout(std_out,message,'COLL')

 call timab(timer_chebfi,2,tsec)

end subroutine chebfi
!!***

!!****f* ABINIT/cheb_poly
!! NAME
!! cheb_poly
!!
!! FUNCTION
!! Computes the value of the Chebyshev polynomial of degree n on the interval [a,b] at x
!!
!! INPUTS
!! x= input variable
!! n= degree
!! a= left bound of the interval
!! b= right bound of the interval
!!
!! OUTPUT
!! y= Tn(x)
!!
!! PARENTS
!!      chebfi
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

function cheb_poly(x, n, a, b) result(y)

 integer, intent(in) :: n
 integer :: i
 real(dp), intent(in) :: x, a, b
 real(dp) :: y, xred, temp
 real(dp) :: yim1

! *************************************************************************

 xred = (x-(a+b)/2)/(b-a)*2
 y = xred
 yim1 = one
 do i=2, n
   temp = y
   y = 2*xred*y - yim1
   yim1 = temp
 end do

end function cheb_poly
!!***

!!****f* ABINIT/cheb_oracle
!! NAME
!! cheb_oracle
!!
!! FUNCTION
!! Returns the number of necessary iterations to decrease residual by at least tol
!! Here as in the rest of the code, the convention is that residuals are squared (||Ax-lx||^2)
!!
!! INPUTS
!! x= input variable
!! a= left bound of the interval
!! b= right bound of the interval
!! tol= needed precision
!! nmax= max number of iterations
!!
!! OUTPUT
!! n= number of iterations needed to decrease residual by tol
!!
!! PARENTS
!!      chebfi
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

function cheb_oracle(x, a, b, tol, nmax) result(n)

 real(dp) :: tol

 integer :: nmax
 integer :: n, i
 real(dp), intent(in) :: x, a, b
 real(dp) :: y, xred, temp
 real(dp) :: yim1

! *************************************************************************

 xred = (x-(a+b)/2)/(b-a)*2
 y = xred
 yim1 = one

 n = nmax
 if(1/(y**2) < tol) then
   n = 1
 else
   do i=2, nmax-1
     temp = y
     y = 2*xred*y - yim1
     yim1 = temp
     if(1/(y**2) < tol) then
       n = i
       exit
     end if
   end do
 end if

end function cheb_oracle
!!***

end module m_chebfi
!!***
