!!****m* ABINIT/m_rttddft_exponential
!! NAME
!!  m_rttddft_exponential
!!
!! FUNCTION
!!  Contains subroutines to compute the exponential
!!  part of the propagator using various approximations
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB)
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

module m_rttddft_exponential

 use defs_basis
 use defs_abitypes,        only: MPI_type

 use m_bandfft_kpt,        only: bandfft_kpt, &
                               & bandfft_kpt_get_ikpt
 use m_dtset,              only: dataset_type
 use m_getghc,             only: multithreaded_getghc
 use m_hamiltonian,        only: gs_hamiltonian_type
 use m_invovl,             only: apply_invovl
 use m_pawcprj,            only: pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_prep_kgb,           only: prep_getghc, prep_index_wavef_bandpp
 use m_profiling_abi,      only: abimem_record
 use m_rttddft_properties, only: rttddft_calc_eig, rttddft_calc_enl
 use m_xmpi,               only: xmpi_alltoallv

 implicit none

 private
!!***

 public :: rttddft_exp_taylor
!!***

contains
!!***

!!****f* m_rttddft_exponential/rttddft_exp_taylor
!!
!! NAME
!!  rttddft_propagator_exp
!!
!! FUNCTION
!!  Applies the propagator exp(-i*dt*H(*S^-1)) on the cg coeffs
!!  approximating the exponential using a Taylor exapansion
!!
!! INPUTS
!!  cg <real(npw*nspinor*nband)> = the wavefunction coefficients
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ham_k <type(gs_hamiltonian_type)> = the Hamiltonian H
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nband_k <integer> = number of bands
!!  npw_k <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!
!! OUTPUT
!!  cg <real(npw*nspinor*nband)> = the new cg after application
!!  of the exponential propagator
!!  enl <real(bandpp)> = non local contribution to the energy in the
!!                       NC case - optional
!!  eig <real(bandpp)> = eigenvalues - optional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine rttddft_exp_taylor(cg,dtset,ham_k,mpi_enreg,nband_k,npw_k,nspinor,enl,eig)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)              :: nband_k
 integer,                   intent(in)              :: npw_k
 integer,                   intent(in)              :: nspinor
 type(dataset_type),        intent(in)              :: dtset
 type(gs_hamiltonian_type), intent(inout)           :: ham_k
 type(MPI_type),            intent(inout)           :: mpi_enreg
 !arrays
 real(dp), target,          intent(inout)           :: cg(2,npw_k*nband_k*nspinor)
 real(dp),                  intent(out),   optional :: enl(:)
 real(dp), pointer,         intent(inout), optional :: eig(:)

 !Local variables-------------------------------
 !scalars
 integer                         :: cpopt
 integer                         :: iorder
 integer                         :: nband_t
 integer                         :: npw_t
 integer                         :: sij_opt
 integer                         :: tim_getghc = 5
 integer                         :: nfact
 logical                         :: l_paw, l_eig, l_enl
 real(dp)                        :: dt
 !arrays
 integer,            allocatable :: index_wavef_band(:)
 type(pawcprj_type), allocatable :: cwaveprj(:,:)
 real(dp), pointer               :: cg_t(:,:) => null()
 real(dp), target,   allocatable :: cg_trans(:,:)
 real(dp),           allocatable :: ghc(:,:)
 real(dp),           allocatable :: gvnlxc_dummy(:,:)
 real(dp),           allocatable :: gsm1hc(:,:)
 real(dp),           allocatable :: gsc(:,:)
 real(dp),           allocatable :: tmp(:,:)

! ***********************************************************************

 !Check what properties we need to compute
 l_eig = .false.; l_enl = .false.
 if (present(eig)) then
    if (associated(eig)) l_eig = .true.
 end if
 if (present(enl)) l_enl = (size(enl)>0)

 !Transpose if paral_kgb
 if (dtset%paral_kgb == 1 .and. mpi_enreg%nproc_band>1) then
   call paral_kgb_transpose(cg,cg_trans,mpi_enreg,nband_t,npw_t,nspinor,1,index_wavef_band)
   cg_t => cg_trans
 else
   npw_t = npw_k
   nband_t = nband_k
   cg_t => cg
 end if

 l_paw = (ham_k%usepaw == 1)
 if(l_paw) then
   ABI_MALLOC(cwaveprj, (ham_k%natom,nspinor*nband_t))
   call pawcprj_alloc(cwaveprj,0,ham_k%dimcprj)
 else
   ABI_MALLOC(cwaveprj,(0,0))
 end if
 cpopt = 0

 dt = dtset%dtele !electronic timestep

 ABI_MALLOC(ghc,   (2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gsc,   (2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gsm1hc,(2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gvnlxc_dummy, (0, 0))

 !*** Taylor expansion ***
 ABI_MALLOC(tmp,(2, npw_t*nspinor*nband_t))
 tmp(:,:) = cg_t(:,:)
 nfact = 1

 do iorder = 1, dtset%td_exp_order

   nfact = nfact*iorder

   !** Apply Hamiltonian
   sij_opt = 0
   if (iorder == 1 .and. l_eig .and. l_paw) sij_opt = 1
   if (dtset%paral_kgb /= 1) then
      call multithreaded_getghc(cpopt,tmp,cwaveprj,ghc,gsc,ham_k,gvnlxc_dummy,1.0_dp, &
                              & mpi_enreg,nband_t,dtset%prtvol,sij_opt,tim_getghc,0)
   else
      call prep_getghc(tmp,ham_k,gvnlxc_dummy,ghc,gsc,1.0_dp,nband_t,mpi_enreg, &
                     & dtset%prtvol,sij_opt,cpopt,cwaveprj,already_transposed=.true.)
   end if

   !** Also apply S^-1 in PAW case
   if (l_paw) then
      call apply_invovl(ham_k,ghc,gsm1hc,cwaveprj,npw_t,nband_t,mpi_enreg,nspinor,dtset%diago_apply_block_sliced)
      tmp(1,:) =  dt*gsm1hc(2,:)/real(nfact,dp)
      tmp(2,:) = -dt*gsm1hc(1,:)/real(nfact,dp)
   else
      tmp(1,:) =  dt*ghc(2,:)/real(nfact,dp)
      tmp(2,:) = -dt*ghc(1,:)/real(nfact,dp)
   end if

   !** Compute properties if requested
   if (iorder == 1) then
      !eigenvalues
      if (l_eig) then
         if (l_paw) then
            call rttddft_calc_eig(cg_t,eig,ghc,ham_k%istwf_k,nband_t,npw_t,nspinor, &
                                & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft,gsc=gsc)
         else
            call rttddft_calc_eig(cg_t,eig,ghc,ham_k%istwf_k,nband_t,npw_t,nspinor, &
                                & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         end if
      end if

      !non local PSP energy contribution in NC case
      if (l_enl) then
         call rttddft_calc_enl(cg_t,enl,ham_k,nband_t,npw_t,nspinor,mpi_enreg)
      end if
   end if

   !** Update cg
   cg_t(:,:) = cg_t(:,:) + tmp(:,:)

 end do

 ABI_FREE(tmp)
 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(gsm1hc)
 ABI_FREE(gvnlxc_dummy)

 if(l_paw) call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 !Transpose back if paral_kgb
 if (dtset%paral_kgb == 1 .and. mpi_enreg%nproc_band > 1) then
    call paral_kgb_transpose(cg,cg_trans,mpi_enreg,nband_t,npw_t,nspinor,-1,index_wavef_band)
 end if

 nullify(cg_t)

 end subroutine rttddft_exp_taylor
!!***

!!****f* m_rttddft_exponential/paral_kgb_transpose
!!
!! NAME
!!  paral_kgb_transpose
!!
!! FUNCTION
!!  if option = 1: Forward transpose
!!    Transpose cg_1 in linalg ((npw/npband),nband) distribution
!!    into cg_2 in fft (npw,bandpp) distribution
!!
!!  if option = -1: Backward transpose
!!    Transpose back cg_2 in fft (npw,bandpp) distribution
!!    into cg_1 linakg (npw/npband),nband) distribution
!!
!! INPUTS
!!  if option = 1:
!!    cg_1 <real((npw/nband)*nspinor*nband)>
!!  if option = -1:
!!    cg_2 <real(npw*nspinor*bandpp)>
!!    nband_t <integer> = number of bands after forward transpose (bandpp)
!!    npw_t <integer> = number of pw after forward transpose (npw_k)
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nspinor <integer> = "dimension of spinors"
!!  option <integer> = option for forward or backward transpose
!!  index_wavef_band <integer> = order of the bands after transpose
!!
!! OUTPUT
!!  if option = 1 :
!!    cg_2 <real(npw*nspinor*bandpp)>
!!    nband_t <integer> = number of bands after forward transpose (bandpp)
!!    npw_t <integer> = number of pw after forward transpose (npw_k)
!!  if option = -1:
!!    cg_1 <real((npw/nband)*nspinor*nband)>

!! SIDE EFFECTS
!!  if option = 1 :
!!    cg_2 has been allocated
!!    index_wavef_band has been allocated
!!  if option = -1:
!!    cg_2 has been deallocated
!!    index_wavef_band has been allocated
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine paral_kgb_transpose(cg_1,cg_2,mpi_enreg,nband_t,npw_t,nspinor,option,index_wavef_band)

implicit none

!Arguments ------------------------------------
!scalars
integer,               intent(inout) :: nband_t
integer,               intent(inout) :: npw_t
integer,               intent(in)    :: nspinor
integer,               intent(in)    :: option
type(MPI_type),        intent(inout) :: mpi_enreg
!arrays
integer,  allocatable, intent(inout) :: index_wavef_band(:)
real(dp),              intent(inout) :: cg_1(:,:)
real(dp), allocatable, intent(inout) :: cg_2(:,:)

!Local variables-------------------------------
!scalars
integer               :: bandpp
integer               :: ierr
integer               :: ikpt_this_proc
!arrays
integer               :: recvcountsloc(mpi_enreg%nproc_band)
integer               :: rdisplsloc(mpi_enreg%nproc_band)
integer               :: sendcountsloc(mpi_enreg%nproc_band)
integer               :: sdisplsloc(mpi_enreg%nproc_band)
real(dp), allocatable :: cg_work(:,:)

! ***********************************************************************

 !Init useful MPI variables
 bandpp = mpi_enreg%bandpp
 ikpt_this_proc = bandfft_kpt_get_ikpt()
 recvcountsloc = bandfft_kpt(ikpt_this_proc)%recvcounts*2*nspinor*bandpp
 rdisplsloc = bandfft_kpt(ikpt_this_proc)%rdispls*2*nspinor*bandpp
 sendcountsloc = bandfft_kpt(ikpt_this_proc)%sendcounts*2*nspinor
 sdisplsloc = bandfft_kpt(ikpt_this_proc)%sdispls*2*nspinor

 !Forward transpose: cg_1 -> cg_2
 if (option == 1) then
   nband_t = bandpp
   npw_t = bandfft_kpt(ikpt_this_proc)%ndatarecv
   ABI_MALLOC(cg_2,  (2, npw_t*nspinor*nband_t))
   ABI_MALLOC(cg_work, (2, npw_t*nspinor*nband_t))
   !Transpose input cg_1 into cg_work
   call xmpi_alltoallv(cg_1,sendcountsloc,sdisplsloc,cg_work, &
                     & recvcountsloc,rdisplsloc,mpi_enreg%comm_band,ierr)
   !properly sort array according to bandd after alltoall
   call prep_index_wavef_bandpp(mpi_enreg%nproc_band,mpi_enreg%bandpp,          &
                              & nspinor,bandfft_kpt(ikpt_this_proc)%ndatarecv , &
                              & bandfft_kpt(ikpt_this_proc)%recvcounts,         &
                              & bandfft_kpt(ikpt_this_proc)%rdispls, index_wavef_band)
   cg_2(:,:) = cg_work(:,index_wavef_band)
   ABI_FREE(cg_work)
 end if

 !Transpose back: cg_2 -> cg_1
 if (option == -1) then
   ABI_MALLOC(cg_work, (2, npw_t*nspinor*nband_t))
   cg_work(:,index_wavef_band) = cg_2(:,:)
   !Transpose cg_work to input cg_1
   call xmpi_alltoallv(cg_work,recvcountsloc,rdisplsloc,cg_1, &
                     & sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   !Free memory
   ABI_FREE(cg_work)
   ABI_FREE(cg_2)
   ABI_FREE(index_wavef_band)
 end if

end subroutine paral_kgb_transpose
!!***

end module m_rttddft_exponential
!!***
