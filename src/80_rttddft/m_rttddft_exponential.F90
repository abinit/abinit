!!****m* ABINIT/m_rttddft_exponential
!! NAME
!!  m_rttddft_exponential
!!
!! FUNCTION
!!  Contains subroutines to compute the exponential 
!!  part of the propagator using various approximations
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  m_rttddft_propagators
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
 use defs_abitypes,   only: MPI_type

 use m_bandfft_kpt,   only: bandfft_kpt,          &
                          & bandfft_kpt_set_ikpt, &
                          & bandfft_kpt_get_ikpt
 use m_cgtools,       only: dotprod_g
 use m_dtset,         only: dataset_type
 use m_energies,      only: energies_type
 use m_getghc,        only: multithreaded_getghc
 use m_hamiltonian,   only: gs_hamiltonian_type
 use m_invovl,        only: apply_invovl
 use m_pawcprj,       only: pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_prep_kgb,      only: prep_getghc
 use m_xmpi,          only: xmpi_alltoallv, xmpi_comm_rank

 implicit none

 private
!!***

 public :: rttddft_exp_taylor
!!***

contains


!!****f* m_rttddft_exponential/rttddft_exp_taylor
!!
!! NAME
!!  rttddft_propagator_exp
!!
!! FUNCTION
!!  Applies the propagator exp(-i*dt*H(*S^-1)) on the cg coeffs
!!  approximating the exponential using different methods
!!
!! INPUTS
!!  cg <real(npw*nspinor*nband)> = the wavefunction coefficients
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = the Hamiltonian H
!!  method <integer> = method use to approximate the exponential
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!
!! OUTPUT
!!  cg <real(npw*nspinor*nband)> = the exponential propagator 
!!   applied to the WF coeeficients
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_propagators/rttddft_propagator_er
!!
!! CHILDREN
!!
!! SOURCE
 subroutine rttddft_exp_taylor(cg,dt,dtset,gs_hamk,method,mpi_enreg,nband_k,npw_k,nspinor,eig)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)            :: method
 integer,                   intent(in)            :: nband_k
 integer,                   intent(in)            :: npw_k
 integer,                   intent(in)            :: nspinor
 type(dataset_type),        intent(in)            :: dtset
 real(dp),                  intent(in)            :: dt
 type(gs_hamiltonian_type), intent(inout)         :: gs_hamk
 type(MPI_type),            intent(inout)         :: mpi_enreg
 !arrays
 real(dp), target,          intent(inout)         :: cg(2,npw_k*nband_k*nspinor)
 real(dp),                  intent(out), optional :: eig(nband_k)
 
 !Local variables-------------------------------
 !scalars
 integer                         :: iband
 integer                         :: nband_t
 integer                         :: npw_t
 integer                         :: iorder
 integer                         :: sij_opt,cpopt
 integer                         :: shift
 integer                         :: tim_getghc = 5
 integer                         :: nfact
 logical                         :: paw
 real(dp)                        :: dprod_r, dprod_i
 !arrays
 type(pawcprj_type), allocatable :: cwaveprj(:,:)
 real(dp), pointer               :: cg_t(:,:)
 real(dp),           allocatable :: cg_work(:,:)
 real(dp),           allocatable :: ghc(:,:)
 real(dp),           allocatable :: gvnlxc_dummy(:,:)
 real(dp),           allocatable :: gsm1hc(:,:)
 real(dp),           allocatable :: gsc(:,:)
 real(dp),           allocatable :: tmp(:,:)
 
! ***********************************************************************

 !Transpose if paral_kgb
 if (dtset%paral_kgb == 1) then
   call paral_kgb_transpose(cg,cg_t,cg_work,mpi_enreg,nband_k,nband_t,npw_k,npw_t,nspinor,1)
   npw_t = npw_k
   nband_t = nband_k
   cg_t => cg
 else
   npw_t = npw_k
   nband_t = nband_k
   cg_t => cg
 end if

 paw = (gs_hamk%usepaw == 1)
 if(paw) then
   !FB: Needs cwaveprj for invovl
   ABI_MALLOC(cwaveprj, (gs_hamk%natom,nspinor*nband_t))
   call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
 else
   ABI_MALLOC(cwaveprj,(0,0))
 end if
 cpopt = 0

 ABI_MALLOC(ghc,    (2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gsc,    (2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gsm1hc, (2, npw_t*nspinor*nband_t))
 ABI_MALLOC(gvnlxc_dummy, (0, 0))
      
 !*** Taylor expansion ***
 if (method == 1) then 

   ABI_MALLOC(tmp, (2, npw_t*nspinor*nband_t))
   tmp(:,:) = cg_t(:,:)
   nfact = 1

   do iorder = 1, dtset%td_exp_order

      nfact = nfact*iorder

      !** Apply Hamiltonian 
      sij_opt = 0
      if (iorder == 1 .and. paw .and. present(eig)) sij_opt = 1
      if (dtset%paral_kgb /= 1) then
         call multithreaded_getghc(cpopt,tmp,cwaveprj,ghc,gsc,gs_hamk,gvnlxc_dummy,1.0_dp, &
                                 & mpi_enreg,nband_t,dtset%prtvol,sij_opt,tim_getghc,0)
      else
         call prep_getghc(tmp,gs_hamk,gvnlxc_dummy,ghc,gsc,1.0_dp,nband_t,mpi_enreg, &
                        & dtset%prtvol,sij_opt,cpopt,cwaveprj,already_transposed=.true.)
      end if

      !** Also apply S^-1 in PAW case
      if (paw) then 
         call apply_invovl(gs_hamk, ghc, gsm1hc, cwaveprj, npw_t, nband_t, mpi_enreg, nspinor)
         tmp(1,:) =  dt*gsm1hc(2,:)/real(nfact,dp)
         tmp(2,:) = -dt*gsm1hc(1,:)/real(nfact,dp)
      else
         tmp(1,:) =  dt*ghc(2,:)/real(nfact,dp)
         tmp(2,:) = -dt*ghc(1,:)/real(nfact,dp)
      end if

      !Update eigenvalues if required
      if (present(eig) .and. iorder == 1) then 
         do iband=1, nband_t
            shift = npw_t*nspinor*(iband-1)
            !Compute eigenvalues
            call dotprod_g(eig(iband),dprod_i,gs_hamk%istwf_k,                   &
                         & npw_t*nspinor,1,ghc(:, shift+1:shift+npw_t*nspinor),  &
                         & cg_t(:, shift+1:shift+npw_t*nspinor),mpi_enreg%me_g0, &
                         & mpi_enreg%comm_spinorfft)
            if(paw) then
               call dotprod_g(dprod_r,dprod_i,gs_hamk%istwf_k,                      &
                            & npw_t*nspinor,1,gsc(:, shift+1:shift+npw_t*nspinor),  &
                            & cg_t(:, shift+1:shift+npw_t*nspinor),mpi_enreg%me_g0, &
                            & mpi_enreg%comm_spinorfft)
               eig(iband) = eig(iband)/dprod_r
            end if
         end do
      end if

      cg_t(:,:) = cg_t(:,:) + tmp(:,:)

   end do

   ABI_FREE(tmp)

 end if

 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(gsm1hc)
 ABI_FREE(gvnlxc_dummy)

 if(paw) call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 !Transpose back if paral_kgb
 !if (dtset%paral_kgb == 1) call paral_kgb_transpose(cg,cg_t,cg_work,mpi_enreg,nband_k,nband_t, &
 !                                                & npw_k,npw_t,nspinor,-1)
 ABI_FREE(cg_work)

 end subroutine rttddft_exp_taylor

!!****f* m_rttddft_exponential/paral_kgb_transpose
!!
!! NAME
!!  paral_kgb_transpose
!!
!! FUNCTION
!!  option = 1: 
!!    Transpose the cg from ((npw/npband),nband) distribution
!!    to (npw,bandpp) distribution 
!!    if paral_kgb: cg_in -> cg_out => cg_work 
!!    else: cg_in -> cg_out => cg_in (no transposition needed if no paral_kgb)
!!
!!  option = -1:
!!    Transpose back from (npw,bandpp) to (npw/npband),nband)
!!    if paral_kgb: cg_work -> cg_in
!!    else: nothing to be done
!!
!! INPUTS
!!  cg <real(npw*nspinor*nband)> = the wavefunction coefficients
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = the Hamiltonian H
!!  method <integer> = method use to approximate the exponential
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = "dimension of spinors"
!!
!! OUTPUT
!!  exp_op <real(npw*nspinor*nband)> = the exponential propagator 
!!                                     applied on the WF (exp(-i*dt*H)|psi>)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine paral_kgb_transpose(cg_in,cg_out,cg_work,mpi_enreg,nband_in,nband_out,npw_in,npw_out,nspinor,option)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,            intent(in)                 :: nband_in
 integer,            intent(out)                :: nband_out
 integer,            intent(in)                 :: npw_in
 integer,            intent(out)                :: npw_out
 integer,            intent(in)                 :: nspinor
 integer,            intent(in)                 :: option
 type(MPI_type),     intent(inout)              :: mpi_enreg
 !arrays
 real(dp), target,   intent(inout)              :: cg_in(2,npw_in*nband_in*nspinor)
 real(dp), pointer,  intent(inout)              :: cg_out(:,:)
 real(dp), target,   intent(inout), allocatable :: cg_work(:,:)
 
 !Local variables-------------------------------
 !scalars
 integer :: bandpp
 integer :: ierr
 integer :: ikpt_this_proc
 !arrays
 integer :: recvcountsloc(mpi_enreg%nproc_band)
 integer :: rdisplsloc(mpi_enreg%nproc_band)
 integer :: sendcountsloc(mpi_enreg%nproc_band)
 integer :: sdisplsloc(mpi_enreg%nproc_band)
 
! ***********************************************************************

 ! Init useful MPI variables
 bandpp = nband_in / mpi_enreg%nproc_band
 !bandpp = mpi_enreg%bandpp
 ikpt_this_proc = bandfft_kpt_get_ikpt()
 recvcountsloc = bandfft_kpt(ikpt_this_proc)%recvcounts*2*nspinor*bandpp
 sendcountsloc = bandfft_kpt(ikpt_this_proc)%sendcounts*2*nspinor*bandpp
 rdisplsloc = bandfft_kpt(ikpt_this_proc)%rdispls*2*nspinor*bandpp
 sdisplsloc = bandfft_kpt(ikpt_this_proc)%sdispls*2*nspinor
 print*, 'FB-test: ikpt_this_proc =', ikpt_this_proc
 print*, 'FB-test: recvcountsloc =', recvcountsloc
 print*, 'FB-test: rdisplsloc =', rdisplsloc
 print*, 'FB-test: sendcountsloc =', sendcountsloc
 print*, 'FB-test: sdisplsloc =', sdisplsloc

 !Transpose: cg_in -> cg_out => cg_work
 if (option == 1) then 

   write(400+xmpi_comm_rank(mpi_enreg%comm_cell),*) 'Transposing..'
   nband_out = bandpp
   !nband_out = nband_in/mpi_enreg%nproc_band
   npw_out = bandfft_kpt(ikpt_this_proc)%ndatarecv
   print*, 'FB-test: nband_out, npw_out, nproc_band =', nband_out, npw_out, mpi_enreg%nproc_band
   
   ABI_MALLOC(cg_work, (2, npw_out*nspinor*nband_out))
   
   print*, 'FB-test: dim cg_work, cg_in=', 2*npw_out*nspinor*nband_out, 2*npw_in*nspinor*nband_in
   
   ! Transpose input cg_in into cg_work
   call xmpi_alltoallv(cg_in,sendcountsloc,sdisplsloc,cg_work, &
                     & recvcountsloc,rdisplsloc,mpi_enreg%comm_band,ierr)
   
   cg_out => cg_work

   write(400+xmpi_comm_rank(mpi_enreg%comm_cell),*) 'Transposed'
 end if
 
 !Transpose back: cg_work -> cg_in
 if (option == -1) then
   write(400+xmpi_comm_rank(mpi_enreg%comm_cell),*) 'Transposing back..'

   !Transpose cg_work to input cg_in
   call xmpi_alltoallv(cg_work,recvcountsloc,rdisplsloc,cg_in, &
                     & sendcountsloc,sdisplsloc,mpi_enreg%comm_band,ierr)
   !Free memory
   ABI_FREE(cg_work)
   write(400+xmpi_comm_rank(mpi_enreg%comm_cell),*) 'Transposed back'
 end if

end subroutine paral_kgb_transpose

end module m_rttddft_exponential
!!***
