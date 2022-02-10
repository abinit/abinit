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

 use m_cgtools,       only: dotprod_g
 use m_dtset,         only: dataset_type
 use m_energies,      only: energies_type
 use m_getghc,        only: getghc
 use m_hamiltonian,   only: gs_hamiltonian_type
 use m_invovl,        only: apply_invovl
 use m_pawcprj,       only: pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_prep_kgb,      only: prep_getghc
 use m_spacepar,      only: meanvalue_g

 implicit none

 private
!!***

 public :: rttddft_exp_taylor
!!***

contains


!!****f* m_rttddft/rttddft_exp_taylor
!!
!! NAME
!!  rttddft_propagator_exp
!!
!! FUNCTION
!!  Applies the propagator exp(-i*dt*H) on the cg coeffs
!!  approximating the exponential using different methods
!!
!! INPUTS
!!  cg <real(npw*nspinor*nband)> = the wavefunction coefficients
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = the Hamiltonian H
!!  ikpt <integer> = indice of the k-point considered here
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
!!  m_rttddft_propagators/rttddft_propagator_er
!!
!! CHILDREN
!!
!! SOURCE
 subroutine rttddft_exp_taylor(cg,dt,dtset,eig,gs_hamk,ikpt,method,mpi_enreg,nband_k,npw_k,nspinor,occ_k,energies)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)              :: ikpt
 integer,                   intent(in)              :: method
 integer,                   intent(in)              :: nband_k
 integer,                   intent(in)              :: npw_k
 integer,                   intent(in)              :: nspinor
 type(dataset_type),        intent(in)              :: dtset
 real(dp),                  intent(in)              :: dt
 type(gs_hamiltonian_type), intent(inout)           :: gs_hamk
 type(MPI_type),            intent(inout)           :: mpi_enreg
 type(energies_type),       intent(inout), optional :: energies
 !arrays
 real(dp),                  intent(inout)           :: cg(2,npw_k*nband_k*nspinor)
 real(dp),                  intent(out)             :: eig(nband_k)
 real(dp),                  intent(in)              :: occ_k(nband_k)
 
 !Local variables-------------------------------
 !scalars
 integer                         :: iband
 integer                         :: iorder
 integer                         :: sij_opt,cpopt
 integer                         :: shift
 integer                         :: tim_getghc = 5
 integer                         :: nfact
 logical                         :: paw
 real(dp)                        :: ar
 real(dp)                        :: dprod_r, dprod_i
 real(dp)                        :: enlx
 !arrays
 real(dp),           allocatable :: ghc(:,:)
 real(dp),           allocatable :: gsm1hc(:,:)
 real(dp),           allocatable :: gsc(:,:)
 real(dp),           allocatable :: gvnlxc(:,:)
 real(dp),           allocatable :: tmp(:,:)
 type(pawcprj_type), allocatable :: cwaveprj(:,:)
 
! ***********************************************************************

 !FB: cwaveprj probably just need to point to the right part of the tdks%cprj array
 !and we should use cptopt=2 option in getghc to avoid recomputing cprj all the time
 paw = (gs_hamk%usepaw == 1)
 if(paw) then
   ABI_MALLOC(cwaveprj, (gs_hamk%natom,nspinor*nband_k))
   call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
   !FB: Not sure about the values of sij_opt and cpopt
   sij_opt = 1 ! recompute S
   cpopt = 0   ! save cprojs
 else
   ABI_MALLOC(cwaveprj,(0,0))
   sij_opt = 0
   cpopt = -1
 end if

 !FB: Probably not needed to create so many arrrays ?
 ABI_MALLOC(ghc,    (2, npw_k*nspinor*nband_k))
 ABI_MALLOC(gsc,    (2, npw_k*nspinor*nband_k))
 ABI_MALLOC(gsm1hc, (2, npw_k*nspinor*nband_k))
 ABI_MALLOC(gvnlxc, (2, npw_k*nspinor*nband_k)) ! @MT: gvnlxc not needed according to Lucas?
      
 !*** Taylor expansion ***
 if (method == 1) then 
   ABI_MALLOC(tmp, (2, npw_k*nspinor*nband_k))
   tmp(:,:) = cg(:,:)
   nfact = 1
   do iorder = 1, dtset%td_exp_order
      !** Apply Hamiltonian 
      ! Using gsm1hc only as a temporary array here
      if (dtset%paral_kgb /= 1) then
         call getghc(cpopt,tmp,cwaveprj,ghc,gsc,gs_hamk,gvnlxc,1.0_dp,mpi_enreg,nband_k, &
                   & dtset%prtvol,sij_opt,tim_getghc,0)
      else
         !FB: already_transposed put to false since we have not done any all_to_all before contrarily to chebfi
         call prep_getghc(tmp,gs_hamk,gvnlxc,ghc,gsc,1.0_dp,nband_k,mpi_enreg,dtset%prtvol, &
                        & sij_opt,cpopt,cwaveprj,already_transposed=.false.)
      end if

      !** Also apply S^-1 in PAW case
      if(paw) then 
         call apply_invovl(gs_hamk, ghc, gsm1hc, cwaveprj, npw_k, nband_k, mpi_enreg, nspinor)
         tmp(1,:) =  dt*gsm1hc(2,:)/real(nfact,dp)
         tmp(2,:) = -dt*gsm1hc(1,:)/real(nfact,dp)
      else
         tmp(1,:) =  dt*ghc(2,:)/real(nfact,dp)
         tmp(2,:) = -dt*ghc(1,:)/real(nfact,dp)
      end if

      if (present(energies) .and. iorder == 1) then 
         do iband=1,nband_k
            if (abs(occ_k(iband))>tol8) then 
               shift = npw_k*nspinor*(iband-1)
               !Compute kinetic energy contribution 
               call meanvalue_g(ar,gs_hamk%kinpw_k,0,gs_hamk%istwf_k,mpi_enreg,npw_k,nspinor,&
                              & cg(:,shift+1:shift+npw_k*nspinor),cg(:,shift+1:shift+npw_k*nspinor),0)
               energies%e_kinetic = energies%e_kinetic + dtset%wtk(ikpt)*occ_k(iband)*ar
               !Compute non local psp energy contribution
               if (gs_hamk%usepaw /= 1) then
                  call dotprod_g(enlx,dprod_i,gs_hamk%istwf_k,npw_k*nspinor,1, &
                               & cg(:,shift+1:shift+npw_k*nspinor),            &
                               & gvnlxc(:,shift+1:shift+npw_k*nspinor),        &
                               & mpi_enreg%me_g0,mpi_enreg%comm_bandspinorfft)
                  energies%e_nlpsp_vfock=energies%e_nlpsp_vfock+dtset%wtk(ikpt)*occ_k(iband)*enlx                    
               end if
               !Compute eigenvalues
               call dotprod_g(eig(iband),dprod_i,gs_hamk%istwf_k,npw_k*nspinor,1,ghc(:, shift+1:shift+npw_k*nspinor),&
                            & cg(:, shift+1:shift+npw_k*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               if(gs_hamk%usepaw == 1) then
                  call dotprod_g(dprod_r,dprod_i,gs_hamk%istwf_k,npw_k*nspinor,1,gsc(:, shift+1:shift+npw_k*nspinor),&
                               & cg(:, shift+1:shift+npw_k*nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                  eig(iband) = eig(iband)/dprod_r
               end if
            end if
         end do
      end if

      cg(:,:) = cg(:,:) + tmp(:,:)
      nfact = nfact*(iorder+1)
   end do
   ABI_FREE(tmp)
 end if

 ABI_FREE(ghc)
 ABI_FREE(gsc)
 ABI_FREE(gsm1hc)
 ABI_FREE(gvnlxc)

 if(paw) call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

 end subroutine rttddft_exp_taylor

end module m_rttddft_exponential
!!***
