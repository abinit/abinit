!!****m* ABINIT/m_gwls_hamiltonian
!! NAME
!! m_gwls_hamiltonian
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (JLJ, BR, MC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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


module m_gwls_hamiltonian

! local modules
use m_gwls_utility
use m_gwls_wf
use m_dtset

! abinit modules
use m_bandfft_kpt
use m_cgtools
use defs_basis
use m_abicore
use m_xmpi
use m_pawang
use m_errors
use m_ab7_mixing
use m_mpinfo
use m_crystal

use defs_abitypes,      only : MPI_type
use m_io_tools,         only : get_unit
use m_hamiltonian,      only : gs_hamiltonian_type, copy_hamiltonian
use m_pawcprj,          only : pawcprj_type
use m_vcoul_dt
use m_vcoul,            only : vcoul_init, vcoul_free
use m_gsphere,          only : gsphere_t, gsph_init, gsph_free, print_gsphere
use m_bz_mesh,          only : kmesh_t, kmesh_init, kmesh_free, kmesh_print, find_qmesh
use m_fft,              only : fftpac, fourwf
use m_getghc,           only : getghc
use m_io_kss,           only : make_gvec_kss

implicit none
save
private
!!***

!Data passed in argument to build_H and build_vxc that we need to copy for the GW calculation.
!build_H
type(dataset_type)            :: dtset                  !Public
type(MPI_type)                :: mpi_enreg              !Public
integer                       :: cpopt
real(dp), allocatable         :: cg(:,:)                !Public
integer                       :: dimffnl
type(gs_hamiltonian_type)     :: gs_hamk
real(dp), allocatable, target :: ffnl(:,:,:,:)
integer,  allocatable, target :: kg_k(:,:)
real(dp), allocatable, target :: kinpw(:)
real(dp), allocatable, target :: ph3d(:,:,:)
real(dp), allocatable, target :: vlocal(:,:,:,:)
!build_vxc
real(dp), allocatable         :: vxc_dg(:,:)

!Other variables initialized in build_H
integer :: ndat
integer :: sij_opt
integer :: tim_getghc
integer :: type_calc
integer :: nfft                                      !Public
integer :: nline                                     !Public
integer :: ngfft(18)
integer :: e                                         !Public
integer :: nspinor                                   !Public
integer :: n1, n2, n3
integer :: n4, n5, n6                                !Public
integer :: npw_k                                     !Public
integer :: ckpt                                      !Public
integer :: mgfft
integer :: tim_fourwf
integer :: i
integer :: nband                                     !Public
integer :: ispden                                    !Public
integer :: v
integer :: nbandv                                    !Public
integer :: mcg
integer :: ktot                                      !Public
!integer :: tmp2i(2)
integer, parameter   :: iovar=6 !137                 !Public
integer, allocatable :: gbound(:,:)
integer, allocatable :: istwfk(:)                    !Public
real(dp) :: eshift
real(dp) :: tolwfr                                   !Public
real(dp) :: ucvol                                    !Public
real(dp) :: weight
real(dp) :: tmpc(2)
real(dp), allocatable :: vxc(:,:,:,:)                !Public (Just the transcript of vxc_dg on the "small" (wfk) real grid)
real(dp), allocatable :: psik1(:,:)                  !Working wf in k space for module subtourines.
real(dp), allocatable :: psik2(:,:)                  !Working wf in k space for module subtourines.
real(dp), allocatable :: psik3(:,:)                  !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psik4(:,:)                  !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psikb1(:,:)                 !Working wf in k space for module subtourines.
real(dp), allocatable :: psikb2(:,:)                 !Working wf in k space for module subtourines.
real(dp), allocatable :: psikb3(:,:)                 !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psikb4(:,:)                 !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psig1(:,:)                  !Working wf in k space for module subtourines.
real(dp), allocatable :: psig2(:,:)                  !Working wf in k space for module subtourines.
real(dp), allocatable :: psig3(:,:)                  !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psig4(:,:)                  !Working wf in k space for module subtourines. RESERVED FOR Hpsik.
real(dp), allocatable :: psir1(:,:,:,:)              !Working wf in real space (sg) for module subtourines.
real(dp), allocatable :: psir2(:,:,:,:)              !Working wf in real space (sg) for module subtourines.
real(dp), allocatable :: psir3(:,:,:,:)              !Working wf in real space (sg) for module subtourines.
real(dp), allocatable :: psidg(:,:)                  !Working wf in real space (dg) for module subtourines.

real(dp), allocatable :: denpot(:,:,:)               !Working array for real space product of wf with fourwf.
real(dp), allocatable :: pcon(:)
real(dp), allocatable :: eig(:)                      !Public
real(dp), allocatable :: dummy2(:,:), dummy3(:,:,:)
real(dp), allocatable :: scprod2(:,:)
type(pawcprj_type), allocatable     :: conjgrprj(:,:)   !Not allocated (currently).

!To enable the use of ABINIT GW tools for constructing the square root of the coulomb operator.
integer                      :: timrev
real(dp)                     :: ecut_eff
integer, pointer             :: gvec(:,:)
type(vcoul_t)                :: Vcp
type(crystal_t)              :: Cryst
type(gsphere_t)          :: Gsphere
type(kmesh_t)           :: Kmesh, Qmesh
character(len=132),pointer   :: title(:)   !SET2NULL
complex(dpc), allocatable    :: vc_sqrt(:)

!MPI over bands requires :
integer :: blocksize                                 !Public
integer :: nbdblock                                  !Public
integer :: npw_g                                     !Public
integer :: ikpt_this_proc                            !Public
integer :: npw_kb
integer :: iblock
integer :: iband
integer :: n
integer :: npw_serial
integer,  pointer :: kg_k_gather(:,:)
real(dp), pointer :: kinpw_gather(:)
real(dp), pointer :: ph3d_gather(:,:,:)
real(dp), pointer :: ffnl_gather(:,:,:,:)
!!***


public :: dtset, mpi_enreg, cg, kinpw
public :: nfft, n4, n5, n6, nline, tolwfr, npw_k, nbandv, eig, nspinor, ckpt
public :: istwfk, nband, iovar, e, ucvol, vxc, ispden, ktot, kg_k
public :: Hpsik, Hpsikc, g_to_r, gr_to_g, kbkb_to_kb
public :: build_H, destroy_H, set_precondition, build_vxc
public :: unset_precondition, precondition, exchange, dft_xc_energy, sqrt_vc_k, precondition_cplx
public :: pcon, wf_block_distribute
public :: blocksize, nbdblock, ikpt_this_proc, npw_g, npw_kb, kg_k_gather, kinpw_gather

public :: dummy2, dummy3,  mgfft, ngfft, weight, gbound, psir1, psir3, tim_fourwf

real(dp), public, allocatable :: kernel_wavefunctions_FFT(:,:,:)
real(dp), public, allocatable :: valence_wavefunctions_FFT(:,:,:)

public :: pc_k_valence_kernel
! public :: CleanupSQMRKernel
!!***

contains


!!****f* m_gwls_hamiltonian/DistributeValenceWavefunctions
!! NAME
!!  DistributeValenceWavefunctions
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_hamiltonian
!!
!! CHILDREN
!!
!! SOURCE

subroutine DistributeValenceWavefunctions()
!--------------------------------------------------------------------------------
!
! This subroutine distributes, once and for all, the valence wavefunctions to
! the FFT configuration. They will thus be ready to be used within the
! susceptibility operator.
!
!--------------------------------------------------------------------------------
integer  :: iblk, mb, v

real(dp), allocatable :: psik_v(:,:)             !wavefunctions in LA format
real(dp), allocatable :: psik_v_alltoall(:,:)    !wavefunctions in FFT format

! *************************************************************************



!===================================================================
! Allocate the global array which will contain the valence states
!===================================================================
ABI_ALLOCATE( valence_wavefunctions_FFT, (2,npw_g,nbdblock))

valence_wavefunctions_FFT = zero


!====================================================
! Allocate working arrays
!====================================================

ABI_ALLOCATE(psik_v,         (2,npw_kb))
ABI_ALLOCATE(psik_v_alltoall,(2,npw_g))



! loop on all blocks of states,
do iblk = 1, nbdblock

! loop on valence states for this block; if the state is conduction, fill with zeros
do mb = 1, blocksize

v = (iblk-1)*blocksize+mb

if (v <= nbandv) then
  psik_v(:,(mb-1)*npw_k+1:mb*npw_k)   = cg(:,(v-1)*npw_k+1:v*npw_k)
else
  psik_v(:,(mb-1)*npw_k+1:mb*npw_k)   = zero
end if

end do

! change configuration of the data
call wf_block_distribute(psik_v,  psik_v_alltoall,1) ! LA -> FFT

! copy data in global array

valence_wavefunctions_FFT(:,:,iblk) = psik_v_alltoall(:,:)

end do

!====================================================
! cleanup
!====================================================

ABI_DEALLOCATE(psik_v)
ABI_DEALLOCATE(psik_v_alltoall)


end subroutine DistributeValenceWavefunctions
!!***

!!****f* m_gwls_hamiltonian/DistributeValenceKernel
!! NAME
!!  DistributeValenceKernel
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gwls_hamiltonian
!!
!! CHILDREN
!!
!! SOURCE

subroutine DistributeValenceKernel()
!--------------------------------------------------------------------------------
!
! This subroutine distributes, once and for all, the kernel of the static
! SQMR operator.
!
! In this first, quick and dirty implementation, we simply distribute
! ALL the valence bands on ALL the FFT groups. This is not efficient, but
! kernel projections are not a bottleneck of the computation.
!
! A better (forthcoming) algorithm would only distribute the actual kernel,
! not all valence bands.
!--------------------------------------------------------------------------------

integer  :: mb, n

real(dp), allocatable :: psik_n(:,:)             !wavefunctions in LA format
real(dp), allocatable :: psik_n_alltoall(:,:)    !wavefunctions in FFT format

! *************************************************************************


!===================================================================
! Allocate the global array which will contain the valence states
!===================================================================
ABI_ALLOCATE( kernel_wavefunctions_FFT, (2,npw_g, nband))

kernel_wavefunctions_FFT = zero

!====================================================
! Allocate working arrays
!====================================================

ABI_ALLOCATE(psik_n,         (2,npw_kb))
ABI_ALLOCATE(psik_n_alltoall,(2,npw_g))


! loop on all valence states
do n = 1, nband

! Copy multiple instances of this valence state in the array;
! this way, each FFT group will have ALL the valence states!
do mb = 1, blocksize
psik_n(:,(mb-1)*npw_k+1:mb*npw_k)   = cg(:,(n-1)*npw_k+1:n*npw_k)
end do

! change configuration of the data
call wf_block_distribute(psik_n,  psik_n_alltoall,1) ! LA -> FFT

! copy data in global array
kernel_wavefunctions_FFT(:,:,n) = psik_n_alltoall(:,:)

end do

!====================================================
! cleanup
!====================================================

ABI_DEALLOCATE(psik_n)
ABI_DEALLOCATE(psik_n_alltoall)

end subroutine DistributeValenceKernel
!!***

!!****f* m_gwls_hamiltonian/pc_k_valence_kernel
!! NAME
!!  pc_k_valence_kernel
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_Projected_AT,gwls_hamiltonian
!!      gwls_lineqsolver,gwls_model_polarisability,gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine pc_k_valence_kernel(psi_inout,n)

!================================================================================
! This routine projects out of the valence kernel. It assumes the input/output
! array is distributed in the FFT configuration, and that the global
! array containing the kernel (defined in this module) is already prepared
! and ready to be used.
!================================================================================

real(dp), intent(inout) :: psi_inout(2,npw_g)
integer , intent(in), optional :: n


real(dp), allocatable :: psi_projected(:,:)

real(dp) :: tmpc(2)

integer :: v, n_max

integer :: mpi_communicator, ierr

! *************************************************************************


ABI_ALLOCATE( psi_projected, (2,npw_g))

mpi_communicator =  mpi_enreg%comm_fft

if (present(n)) then
  n_max = n
else
  n_max = nbandv
end if

!====================================================
! Compute projection
!====================================================

psi_projected(:,:) = zero

do v = 1, n_max

! compute overlap of kernel member with function
tmpc  = cg_zdotc(npw_g ,kernel_wavefunctions_FFT(:,:,v),psi_inout(:,:))

! Communicate results
call xmpi_sum(tmpc, mpi_communicator, ierr) ! sum on all processors working on FFT!

! add overlap with array to project
psi_projected(1,:) = psi_projected(1,:) + (tmpc(1)*kernel_wavefunctions_FFT(1,:,v)-tmpc(2)*kernel_wavefunctions_FFT(2,:,v))
psi_projected(2,:) = psi_projected(2,:) + (tmpc(1)*kernel_wavefunctions_FFT(2,:,v)+tmpc(2)*kernel_wavefunctions_FFT(1,:,v))

!psi_inout(1,:) = psi_inout(1,:) -( tmpc(1)*kernel_wavefunctions_FFT(1,:,v)-tmpc(2)*kernel_wavefunctions_FFT(2,:,v) )
!psi_inout(2,:) = psi_inout(2,:) -( tmpc(1)*kernel_wavefunctions_FFT(2,:,v)+tmpc(2)*kernel_wavefunctions_FFT(1,:,v) )


end do

if (present(n)) then
  psi_inout(:,:) = psi_projected(:,:)
else
  psi_inout(:,:) = psi_inout(:,:) - psi_projected(:,:)
end if



ABI_DEALLOCATE( psi_projected)

end subroutine pc_k_valence_kernel
!!***

!!****f* m_gwls_hamiltonian/wf_block_distribute
!! NAME
!!  wf_block_distribute
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_GWlanczos,gwls_LanczosBasis,gwls_Projected_AT
!!      gwls_Projected_BT,gwls_hamiltonian,gwls_model_polarisability
!!      gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine wf_block_distribute(psik, psik_alltoall, direction)

!================================================================================
!
! This subroutine distributes a block of wavefunctions in the "linear algebra"
! configuration, to a configuration appropriate to perform FFT (or apply Hamiltonian).
!
! The code below is inspired by the subroutine prep_getghc, which rearranges
! data over MPI processors prior to applying the Hamiltonian.
!
! input:
!             npw_kb :   dimension of wavefunctions in the "linear algebra" configuration,
!                               which means the G-vectors are distributed over all
!                               processors, and every processor has information about all the bands.
!
!
!             npw_g       :   dimension of wavefunctions in the "FFT" configuration,
!                               which means the G-vectors are distributed along rows of the
!                               MPI topology, and a given row only has one band.
!
!              direction    :   Are we going from LA to FFT, or from FFT to LA?
!                              1 : LA  -> FFT
!                              2 : LA <-  FFT
!
! input/output:
!              psik          : block of wavefunctions, in LA configuration
!              psik_alltoall : wavefunction, in FFT configuration
!================================================================================

integer , intent(in)    :: direction                   ! flag which determines the direction of the transfer
real(dp), intent(inout) :: psik(2,npw_kb)       ! block of wavefunctions, in "linear algebra" configuration
real(dp), intent(inout) :: psik_alltoall(2,npw_g)    ! wavefunction in "FFT" configuration; a single band, but more G-vectors


integer  :: ier, spaceComm

integer,allocatable :: rdisplsloc(:)
integer,allocatable :: recvcountsloc(:)
integer,allocatable :: sdisplsloc(:)
integer,allocatable :: sendcountsloc(:)

integer :: nproc_band,  bandpp

integer,pointer :: rdispls(:)
integer,pointer :: recvcounts(:)
integer,pointer :: sdispls(:)
integer,pointer :: sendcounts(:)

! *************************************************************************


! extract information in order to perform MPI communication.
! This code comes from prep_getghc
nproc_band = mpi_enreg%nproc_band
bandpp     = mpi_enreg%bandpp

if(mpi_enreg%nproc_band*mpi_enreg%bandpp > 1) then

  spaceComm=mpi_enreg%comm_fft
  if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band

  ABI_ALLOCATE(sendcountsloc,(nproc_band))
  ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
  ABI_ALLOCATE(recvcountsloc,(nproc_band))
  ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

  recvcounts   =>bandfft_kpt(ikpt_this_proc)%recvcounts(:)
  sendcounts   =>bandfft_kpt(ikpt_this_proc)%sendcounts(:)
  rdispls      =>bandfft_kpt(ikpt_this_proc)%rdispls   (:)
  sdispls      =>bandfft_kpt(ikpt_this_proc)%sdispls   (:)

  recvcountsloc(:)= recvcounts(:)*2*nspinor*bandpp
  rdisplsloc(:)   = rdispls(:)*2*nspinor*bandpp
  sendcountsloc(:)= sendcounts(:)*2*nspinor
  sdisplsloc(:)   = sdispls(:)*2*nspinor

  ! use MPI to communicate information!
  if (direction == 1) then
    ! LA -> FFT
    call xmpi_alltoallv(psik, sendcountsloc, sdisplsloc, psik_alltoall, recvcountsloc,rdisplsloc, spaceComm, ier)

  else if (direction == 2) then
    ! FFT -> LA

    call xmpi_alltoallv(psik_alltoall,recvcountsloc,rdisplsloc, psik, sendcountsloc,sdisplsloc,spaceComm,ier)

  end if

  ABI_DEALLOCATE(sendcountsloc)
  ABI_DEALLOCATE(sdisplsloc   )
  ABI_DEALLOCATE(recvcountsloc)
  ABI_DEALLOCATE(rdisplsloc   )

else

  if(direction == 1) psik_alltoall = psik

  if(direction == 2) psik = psik_alltoall

end if

end subroutine wf_block_distribute
!!***



!!****f* m_gwls_hamiltonian/exchange
!! NAME
!!  exchange
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function exchange(e, Lbasis_lanczos)

!use m_bandfft_kpt
use m_cgtools
!================================================================================
! This subroutine computes the exchange energy in band+FFT parallel
!
!================================================================================
real(dp) :: exchange

integer, intent(in) :: e

! If these arguments are provided, the exchange energy is to be projected on this subspace
complex(dpc), optional, intent(in) :: Lbasis_lanczos(:,:)  ! complex array which contains the Lanczos basis

real(dp), allocatable :: psik_e(:,:)             !Working array to store the wavefunction

real(dp), allocatable :: psik_v(:,:)             !Working array to store the wavefunction

real(dp), allocatable :: psik_out(:,:)           !Working array to store the wavefunction


integer :: iblk, mb

integer :: l, lmax

real(dp) :: tmpc(2)

logical  :: project

! *************************************************************************

!--------------------------------------------------------------------------------
! Determine if the exhcange energy must be projected on the Lanczos basis
! a truncated Coulomb potential
!--------------------------------------------------------------------------------

project = .false.
if (present(Lbasis_lanczos)) then
  project = .true.
  lmax = size(Lbasis_lanczos, 2)
end if

!--------------------------------------------------------------------------------
! The goal of this routine is to compute the exchange energy using
! a truncated Coulomb potential
!--------------------------------------------------------------------------------

exchange = 0.0_dp

!====================================================
! build the block of wavefunctions which will
! contain copies the e-state
!====================================================

ABI_ALLOCATE(psik_e,(2,npw_kb))

! fill psik_e with as many copies of the e-state as there are band processors;
! that way, upon LA -> FFT, each row of fft processors will have the e-state
do v = 1, blocksize
psik_e(:,(v-1)*npw_k+1:v*npw_k)   = cg(:,(e-1)*npw_k+1:e*npw_k)
end do

!====================================================
! Allocate the block of wavefunctions which will
! contain the valence states
!====================================================

ABI_ALLOCATE(psik_v,         (2,npw_kb))
ABI_ALLOCATE(psik_out,         (2,npw_kb))

! loop on all blocks of states,
do iblk = 1, nbdblock

! loop on valence states for this block; if the state is conduction, fill with zeros
do mb = 1, blocksize

v = (iblk-1)*blocksize+mb

if (v <= nbandv) then
  psik_v(:,(mb-1)*npw_k+1:mb*npw_k)   = cg(:,(v-1)*npw_k+1:v*npw_k)
else
  psik_v(:,(mb-1)*npw_k+1:mb*npw_k)   = zero
end if

end do

call kbkb_to_kb(psik_out,psik_v,psik_e)

! apply Coulomb potential, and take norm: cumulate the exchange energy
do mb = 1, blocksize

psik1 = psik_out(:,(mb-1)*npw_k+1:mb*npw_k)
call sqrt_vc_k(psik1)

if (project) then
  ! project on the Lanczos basis
  do l = 1, lmax
  psik2(1,:) = dble(Lbasis_lanczos(:,l))
  psik2(2,:) = dimag(Lbasis_lanczos(:,l))
  tmpc = scprod_k(psik2,psik1)
  exchange = exchange - (tmpc(1)**2+tmpc(2)**2)
  end do
else
  ! compute the "exact" exchange energy
  exchange = exchange - norm_k(psik1)**2
end if

end do

end do

ABI_DEALLOCATE(psik_e)

ABI_DEALLOCATE(psik_v)

ABI_DEALLOCATE(psik_out)

end function exchange
!!***


!!****f* m_gwls_hamiltonian/dft_xc_energy
!! NAME
!!  dft_xc_energy
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

function dft_xc_energy(e)

real(dp) :: dft_xc_energy
integer, intent(in) :: e


integer :: cplex, option, ierr
real(dp), allocatable :: psik_e(:,:)             !Working array to store the wavefunction
real(dp), allocatable :: psik_e_alltoall(:,:)    !Working array to store the wavefunction

real(dp), allocatable :: psik_out(:,:)  !Working array to store the wavefunction

real(dp) :: tmpc(2)
integer  :: mpi_communicator
real(dp) :: dft_xc_energy_tmp

! *************************************************************************

!--------------------------------------------------------------------------------
! The goal of this routine is to compute the xc energy using
! the DFT Vxc array.
!--------------------------------------------------------------------------------

! Parallel-MPI code. This is inspired from code within the getghc routine, which
! applies a complex potential to a wavefunction.


cplex  = 1 ! real potential
option = 2 ! multiply wavefunction by potential
nspinor= 1
! Allocate wavefunction arrays which contains the coefficients of the e-state
ABI_ALLOCATE(psik_e,         (2,npw_kb))
ABI_ALLOCATE(psik_e_alltoall,(2,npw_g))
ABI_ALLOCATE(psik_out,       (3,npw_g))


! Only fetch the e-state, setting other states in block to zero!
psik_e(:,:)       = zero
psik_e(:,1:npw_k) = cg(:,(e-1)*npw_k+1:e*npw_k)
psik_out(:,:)     = zero


! change the configuration of the data
call wf_block_distribute(psik_e,  psik_e_alltoall,1) ! LA -> FFT


! Call fourwf to generate the product, in k-space
! Computation:
!                       psik_e_alltoall(k) -> psi(r)
!                       res(r)   =  (vxc(r) x psi(r))
!                       psik_out(k) <- res(r)
!  psir3 is a dummy, not used here.

call fourwf(cplex,vxc(:,:,:,ispden),psik_e_alltoall,psik_out,psir3,gbound,gbound,istwfk(ckpt),kg_k_gather,kg_k_gather,mgfft,&
&             mpi_enreg,1,ngfft,npw_g,npw_g,n4,n5,n6,option,tim_fourwf,weight,weight)


tmpc = cg_zdotc(npw_g, psik_e_alltoall,psik_out)

mpi_communicator =  mpi_enreg%comm_fft
call xmpi_sum(tmpc,mpi_communicator, ierr) ! sum on all processors working on FFT!

dft_xc_energy_tmp = tmpc(1)

mpi_communicator =  mpi_enreg%comm_band
call xmpi_sum(dft_xc_energy_tmp,mpi_communicator, ierr) ! sum on all processors working on FFT!

dft_xc_energy = dft_xc_energy_tmp

ABI_DEALLOCATE(psik_e)
ABI_DEALLOCATE(psik_e_alltoall)
ABI_DEALLOCATE(psik_out)

end function dft_xc_energy
!!***

!!****f* m_gwls_hamiltonian/set_precondition
!! NAME
!!  set_precondition
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_LanczosResolvents,gwls_lineqsolver
!!
!! CHILDREN
!!
!! SOURCE

subroutine set_precondition(lambda,omega)
!--------------------------------------------------------------------------------
! This subroutine preconditions the problem
!
!                                 A x = b,
!
! with:
!
!        omega         lambda              Operator
!        ------------------------------------------------------
!     absent        absent         A =   (H - lambda_0)  (value of lambda_0 not important)
!     present       present        A =   (H - lambda)^2 + omega^2
!                         other cases not implemented
!
!
! In the above, b = psik.
!
!--------------------------------------------------------------------------------

! TODO :
! - eliminate the 2 "if(kinpw(i) < huge(0.0_dp)*1.0d-11)"
!   since ecutsm = 0.0 always (check if that's true in this gw_sternheimer subroutine).

real(dp), intent(in), optional :: lambda, omega

real(dp) :: poly, x
logical :: omega_imaginary


!integer,save :: counter = 0
!integer      :: io_unit
!character(18):: filename
!logical      :: file_exists

! *************************************************************************


if (present(lambda) .and. present(omega))  then
  omega_imaginary        =.true.
else
  omega_imaginary        =.false.
end if

!io_unit  = get_unit()
!filename = "PRECONDITIONER.log"
!inquire(file=filename,exist=file_exists)

!if (file_exists) then
!      open(io_unit,file=filename,position='append',status=files_status_old)
!else
!      open(io_unit,file=filename,status=files_status_new)
!        write(io_unit,10) "#======================================================================================="
!        write(io_unit,10) "#                                                                                       "
!        write(io_unit,10) "#   This file contains information regarding the preconditioning scheme for SQMR.       "
!        write(io_unit,10) "#                                                                                       "
!        write(io_unit,10) "#======================================================================================="
!end if

!counter = counter + 1
!write(io_unit,10) "#                                                                                       "
!write(io_unit,15) "#   Call #:", counter
!if (present(lambda) .and. present(omega))  then
!        write(io_unit,10) "#    lambda and omega are present: imaginary frequency case"
!else
!        write(io_unit,10) "#    lambda and omega are absent: omega = 0 case"
!end if

!do i=1,npw_k

do i = 1, npw_g

if(omega_imaginary) then
  x = (kinpw_gather(i)-lambda)**2 + omega**2
else
  x = kinpw_gather(i)
end if

if(x < huge(0.0_dp)*1.0d-11) then
  poly    = 27.0 + x*(18.0 + x*(12.0 + 8.0*x))
  pcon(i) = poly/(poly + 16.0*(x**4))
  !pcon(i) = 1.0/(1.0+x) !I don't know why, it gives better results for Silane than the above polynomial.
else
  pcon(i) = zero
end if
end do

!write(io_unit,30) "                prec(       1:   10) =  ",pcon(1:10)
!write(io_unit,30) "                prec(npw_k-10:npw_k) =  ",pcon(npw_k-10:npw_k)

!close(io_unit)

!10 format(A)
!15 format(A,I8)
!30 format(A,1000(ES24.12))

end subroutine set_precondition
!!***

!!****f* m_gwls_hamiltonian/unset_precondition
!! NAME
!!  unset_precondition
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_lineqsolver
!!
!! CHILDREN
!!
!! SOURCE

subroutine unset_precondition()

! *************************************************************************

pcon = one
end subroutine unset_precondition
!!***

!!****f* m_gwls_hamiltonian/precondition
!! NAME
!!  precondition
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_lineqsolver
!!
!! CHILDREN
!!
!! SOURCE

subroutine precondition(psi_out,psi_in)

real(dp), intent(out) :: psi_out(2,npw_g)
real(dp), intent(in)  :: psi_in(2,npw_g)

! *************************************************************************

do i=1,npw_g
psi_out(:,i) = psi_in(:,i)*pcon(i)
end do

end subroutine precondition
!!***

!!****f* m_gwls_hamiltonian/precondition_cplx
!! NAME
!!  precondition_cplx
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_lineqsolver
!!
!! CHILDREN
!!
!! SOURCE

subroutine precondition_cplx(psi_out,psi_in)

complex(dpc), intent(out) :: psi_out(npw_g)
complex(dpc), intent(in)  :: psi_in(npw_g)

! *************************************************************************

do i=1,npw_g
psi_out(i) = psi_in(i)*pcon(i)
end do
end subroutine precondition_cplx
!!***

!!****f* m_gwls_hamiltonian/sqrt_vc_k
!! NAME
!!  sqrt_vc_k
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_GWlanczos,gwls_LanczosBasis,gwls_hamiltonian
!!      gwls_model_polarisability,gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine sqrt_vc_k(psi_inout)

!External variables
real(dp), intent(inout) :: psi_inout(2,npw_k)

!Internal variable
complex(dpc) :: c

! *************************************************************************

do i=1,npw_k
c = vc_sqrt(i) * cmplx(psi_inout(1,i),psi_inout(2,i),dpc)
psi_inout(1,i) = dble (c)
psi_inout(2,i) = dimag(c)
end do
end subroutine sqrt_vc_k
!!***

!!****f* m_gwls_hamiltonian/Hpsik
!! NAME
!!  Hpsik
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_LanczosResolvents,gwls_hamiltonian
!!      gwls_lineqsolver,gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine Hpsik(psi_out,psi_in,cte)

!External variables
real(dp), intent(inout) :: psi_out(2,npw_g)
real(dp), intent(inout), optional :: psi_in(2,npw_g)
real(dp), intent(in), optional :: cte

! *************************************************************************

!psikb1 is allocated in init_hamiltonian and destroyed in destroy_hamiltonian, to avoid doing it at each application of
!the hamiltonian...

!We are keeping track of how the module treat each argument of the hamiltonian (subroutine getghc) here.
!Legend : T = transmitted from 79_seqpar_mpi/vtowfk.F90 through the call gwls_sternheimer
!         I = input argument (allocated and filled)
!         O = output argument (allocated)
!         D = dummy argument (declared but not allocated. Read/Write attempts should trigger a segfault.)
!call getghc(cpopt,             T
!psi,              I
!conjgrprj,        D
!Hpsik,            O
!gsc,              D              (output for PAW : <G|S|C>)
!gs_hamk,          T
!gvnlxc,            D              (<G|Vnonlocal+VFockACE|C>)
!eshift,           I              (<G|H-eshift.S|C>)
!mpi_enreg,        T
!ndat,             Fixed to 1     (# of FFTs to do in //)
!dtset%prtvol,     T
!sij_opt,          Fixed to 0     (PAW dependant : 0-><G|H|C> ; 1-><G|H|C> & <G|S|C> ; -1-><G|H-cte.S|C>)
!tim_getghc,       Fixed to 0     (identity of the timer of the calling subroutine. 1:cgwf 5:lobpcg)
!type_calc)        Fixed to 0     (0:whole H N:some part only)

! It is necessary to pass this optional argument to getghc if paral_kgb=1, or else the code exits cleanly with a "BUG" message...
! It will be important to properly understand what this means when we start using kpoints...

if(present(psi_in)) then

  call getghc(cpopt,psi_in,conjgrprj,psi_out,dummy2,gs_hamk,psig4,eshift,&
&             mpi_enreg,ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)

else
  psig3 = psi_out

  call getghc(cpopt,psig3,conjgrprj,psi_out,dummy2,gs_hamk,psig4,eshift,&
&             mpi_enreg,ndat,dtset%prtvol,sij_opt,tim_getghc,type_calc)

end if
if(present(cte)) then
  if(present(psi_in)) then
    psi_out = psi_out - cte*psi_in
  else
    psi_out = psi_out - cte*psig3
  end if
end if

end subroutine Hpsik
!!***

!!****f* m_gwls_hamiltonian/Hpsikc
!! NAME
!!  Hpsikc
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_lineqsolver
!!
!! CHILDREN
!!
!! SOURCE

subroutine Hpsikc(psi_out,psi_in,cte)

!External variables
complex(dpc), intent(out) :: psi_out(npw_g)
complex(dpc), intent(in)  :: psi_in(npw_g)
complex(dpc), intent(in), optional :: cte

! *************************************************************************


psig1(1,:) = dble(psi_in)
psig1(2,:) = dimag(psi_in)

call Hpsik(psig1)

psi_out = dcmplx(psig1(1,:),psig1(2,:))

if(present(cte)) psi_out = psi_out - cte*psi_in
end subroutine Hpsikc
!!***

!!!****f* m_gwls_hamiltonian/pc_k
!!! NAME
!!!  pc_k
!!!
!!! FUNCTION
!!!  .
!!!
!!! INPUTS
!!!
!!! OUTPUT
!!!
!!! PARENTS
!!!      gwls_DielectricArray,gwls_GWDebugAlgorithms_3,gwls_ProjectedEpsilon
!!!      gwls_Projected_AT,gwls_lineqsolver,gwls_model_polarisability
!!!
!!! CHILDREN
!!!
!!! SOURCE
!
!subroutine pc_k(psi_inout,n,eig_e,above)
!
!real(dp), intent(inout) :: psi_inout(2,npw_kb)
!integer , intent(in), optional :: n
!real(dp), intent(in), optional :: eig_e
!logical, intent(in), optional :: above      !Has an effect only if n is also given in argument
!
!!Local variables
!real(dp),parameter :: degeneracy_tolerance = 2.0e-8
!real(dp),parameter :: projection_tolerance = 1.0e-16
!real(dp) :: z(2)
!
!integer  :: mpi_communicator
!
!! *************************************************************************
!
!
!mpi_communicator =  mpi_enreg%comm_bandfft
!
!if( present(eig_e) ) then
!  do i = 1, nband
!  if (abs(eig_e-eig(i)) < degeneracy_tolerance) then
!    z(:) = scprod_k(cg(:,(i-1)*npw_k+1:i*npw_k),psi_inout)
!    ! project it out!
!    if ( sqrt(z(1)**2+z(2)**2) > projection_tolerance) then
!      ! only project if z is large enough; if we project when z is
!      ! very small, we introduce roundoff error!
!      psi_inout(1,1:npw_k)  = psi_inout(1,:) - ( z(1)*cg(1,(i-1)*npw_k+1:i*npw_k)-z(2)*cg(2,(i-1)*npw_k+1:i*npw_k) )
!      psi_inout(2,1:npw_k)  = psi_inout(2,:) - ( z(1)*cg(2,(i-1)*npw_k+1:i*npw_k)+z(2)*cg(1,(i-1)*npw_k+1:i*npw_k) )
!    end if
!  end if
!  end do
!  elseif( present(n) ) then
!  !If there is a state "n" given in argument, then project on subspace <=n.
!
!  if ( n == 0) then
!    ! if n = 0, there are no states with a lower energy! Thus, the projection vanishes.
!    if(.not. (present(above) .and. above)) psi_inout = zero
!  else
!    if(present(above) .and. above) then
!      !call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,mpi_enreg,0,n,npw_k,nspinor,dtset%ortalg,1,dummy2,scprod2,0,0,0)
!      call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,0,n,npw_k,nspinor,dummy2,scprod2,0,0,0,mpi_enreg%me_g0,mpi_communicator)
!    else
!      psikb4 = psi_inout
!      !call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,mpi_enreg,0,n,npw_k,nspinor,dtset%ortalg,1,dummy2,scprod2,0,0,0)
!      call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,0,n,npw_k,nspinor,dummy2,scprod2,0,0,0,mpi_enreg%me_g0,mpi_communicator)
!      psi_inout = psikb4 - psi_inout
!    end if
!  end if
!else
!  !If there is no "n" given in argument, then project on conduction states.
!  !call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,mpi_enreg,0,nbandv,npw_k,nspinor,dtset%ortalg,1,dummy2,scprod2,0,0,0)
!  call projbd(cg,psi_inout,-1,0,0,istwfk(ckpt),mcg,0,nbandv,npw_k,nspinor,dummy2,scprod2,0,0,0,mpi_enreg%me_g0,mpi_communicator)
!end if
!end subroutine pc_k
!!!***

!!****f* m_gwls_hamiltonian/g_to_r
!! NAME
!!  g_to_r
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_GWlanczos,gwls_GenerateEpsilon
!!      gwls_LanczosBasis,gwls_hamiltonian,gwls_model_polarisability
!!      gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine g_to_r(psi_out,psi_in)

real(dp), intent(out) :: psi_out(2,n4,n5,n6)
real(dp), intent(in)  :: psi_in(2,npw_g)
integer :: option, cplex
integer :: nproc_fft, me_fft, nd3 !, ierr

! *************************************************************************

option = 0 ! fft wavefunction to real space
cplex  = 2 ! complex potential

psig4 = psi_in
psi_out = zero
call fourwf(cplex,dummy3,psig4,dummy2,psi_out,gbound,gbound,istwfk(ckpt),kg_k_gather,kg_k_gather,mgfft,mpi_enreg, &
1,ngfft,npw_g,npw_g,n4,n5,n6,option,tim_fourwf,weight,weight)
psi_out = psi_out/sqrt(ucvol)

!! This comes from prep_fourwf
!nproc_fft=mpi_enreg%nproc_fft
!if (nproc_fft>1) then
!  me_fft=mpi_enreg%me_fft
!  if (me_fft>0) then
!    nd3=(ngfft(3)-1)/nproc_fft+1
!    psi_out(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=psi_out(:,:,:,1:nd3)
!    psi_out(:,:,:,1:nd3)=zero
!  end if
!  call xmpi_sum(psi_out,mpi_enreg%comm_fft,ierr)
!end if

!Instead of communications the whole real space vector on all comm_fft CPUs,
!it is possible to let each CPU keep it's part only and communicate only the sums.
!Then, however, we need to clean the trash in the part of the real space
!vector that's not attributed to the given comm_fft CPU.
nproc_fft=mpi_enreg%nproc_fft
if (nproc_fft>1) then
  me_fft=mpi_enreg%me_fft
  if (me_fft>0) then
    nd3=(ngfft(3)-1)/nproc_fft+1
    psi_out(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=psi_out(:,:,:,1:nd3)
    psi_out(:,:,:,1:nd3)=zero
  end if
  !!!  call xmpi_sum(psi_out,mpi_enreg%comm_fft,ierr)
end if

end subroutine g_to_r
!!***

!!****f* m_gwls_hamiltonian/gr_to_g
!! NAME
!!  gr_to_g
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_DielectricArray,gwls_GWlanczos,gwls_LanczosBasis,gwls_hamiltonian
!!      gwls_model_polarisability,gwls_polarisability
!!
!! CHILDREN
!!
!! SOURCE

subroutine gr_to_g(psig_out,psir_in,psig_in)

real(dp), intent(in)  :: psir_in(2,n4,n5,n6)
real(dp), intent(in), optional :: psig_in(2,npw_g)
real(dp), intent(out) :: psig_out(2,npw_g)

integer :: i1, i2, i3
integer :: cplex, option

! *************************************************************************

cplex = 2 ! complex potential
option= 2 ! multiply wavefunction by potential

if(.not. present(psig_in)) then
  psig4(:,:) = psig_out(:,:)
else
  psig4(:,:) = psig_in(:,:)
end if

do i3=1,n6
do i2=1,n5
do i1=1,n4
denpot(2*i1-1,i2,i3)= psir_in(1,i1,i2,i3)
denpot(2*i1  ,i2,i3)= psir_in(2,i1,i2,i3)
end do
end do
end do

call fourwf(          cplex, & ! complex potential
denpot, & ! real space wavefunction, in denpot format
psig4, & ! fourier space wavefunction
psig_out, & ! result, in FFT configuration
psir3,gbound,gbound,istwfk(ckpt),kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1, & ! Various other arguments
ngfft,npw_g,npw_g,n4,n5,n6,option,tim_fourwf,weight,weight)

end subroutine gr_to_g
!!***

!!****f* m_gwls_hamiltonian/kbkb_to_kb
!! NAME
!!  kbkb_to_kb
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_hamiltonian
!!
!! CHILDREN
!!
!! SOURCE

subroutine kbkb_to_kb(psik_out,psik_in_1,psik_in_2)
!----------------------------------------------------------------------------------------------------
! This function computes the direct product of two wavefunctions in real space,
!         psi_out(r) = psi_in_1^*(r)*psi_in_2(r)  (without complex conjugating any of the 2 wavefunctions)
! and returns the result in fourier space, psi_out(k).
!
! The two input wavefunctions are in k space.
!
!
!----------------------------------------------------------------------------------------------------
real(dp), intent(out) :: psik_out(2,npw_kb)
real(dp), intent(inout)  :: psik_in_1(2,npw_kb), psik_in_2(2,npw_kb)

! *************************************************************************

! change configuration of the data : LA -> FFT
call wf_block_distribute(psik_in_1, psig1,1)
call wf_block_distribute(psik_in_2, psig2,1)

! Fourier transform the first input
call g_to_r(psir1, psig1)

! Complex conjugate in real space the first input
psir1(2,:,:,:) = -psir1(2,:,:,:)

! Fourrier transform the second input, multiply component-wise
call gr_to_g(psig2,psir1)

! change configuration of the output : FFT -> LA
call wf_block_distribute(psik_out, psig2,2)

end subroutine kbkb_to_kb
!!***

!!****f* m_gwls_hamiltonian/build_vxc
!! NAME
!!  build_vxc
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine build_vxc(vxc2,nfft2,nspden2)
!Only transcribe the argument vxc2 in the module; the change from dg to sg is done in build_H (and stored in vxc), since the
!arguments of fftpac are built in build_H.

!We need the dimensions of vxc since they don't exist yet in the module; build_vxc being called before build_H.
integer, intent(in) :: nfft2, nspden2
real(dp), intent(in) :: vxc2(nfft2,nspden2)

! *************************************************************************

ABI_ALLOCATE(vxc_dg,(nfft2,nspden2))
vxc_dg = vxc2

end subroutine build_vxc
!!***

!!****f* m_gwls_hamiltonian/destroy_H
!! NAME
!!  destroy_H
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gwls_sternheimer
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_H

call dtset%free()
call gs_hamk%free()

call cryst%free()
call kmesh_free(Kmesh)
call kmesh_free(Qmesh)
call gsph_free(Gsphere)
call vcoul_free(Vcp)

call bandfft_kpt_destroy_array(bandfft_kpt,mpi_enreg)
call destroy_mpi_enreg(mpi_enreg)

!NOTE : the syntax if(allocated(a)) ABI_DEALLOCATE(a) result in an error if "a" is not allocated; since the macro replace
!ABI_ALLOCATE by more than one line of text, the second lines and up get outside the if... if() then syntax is equired.
if(allocated(cg)) then
  ABI_DEALLOCATE(cg)
end if
if(allocated(gbound)) then
  ABI_DEALLOCATE(gbound)
end if
if(allocated(kg_k)) then
  ABI_DEALLOCATE(kg_k)
end if
if(allocated(ffnl)) then
  ABI_DEALLOCATE(ffnl)
end if
if(allocated(ph3d)) then
  ABI_DEALLOCATE(ph3d)
end if
if(allocated(kinpw)) then
  ABI_DEALLOCATE(kinpw)
end if
if(allocated(vxc)) then
  ABI_DEALLOCATE(vxc)
end if
if(allocated(vlocal)) then
  ABI_DEALLOCATE(vlocal)
end if
if(allocated(conjgrprj)) then
  ABI_DATATYPE_DEALLOCATE(conjgrprj)
end if
if(allocated(istwfk)) then
  ABI_DEALLOCATE(istwfk)
end if
if(allocated(dummy2)) then
  ABI_DEALLOCATE(dummy2)
end if
if(allocated(dummy3)) then
  ABI_DEALLOCATE(dummy3)
end if
if(allocated(eig)) then
  ABI_DEALLOCATE(eig)
end if
if(allocated(scprod2)) then
  ABI_DEALLOCATE(scprod2)
end if
if(allocated(pcon)) then
  ABI_DEALLOCATE(pcon)
end if
if(allocated(psik1)) then
  ABI_DEALLOCATE(psik1)
end if
if(allocated(psik2)) then
  ABI_DEALLOCATE(psik2)
end if
if(allocated(psik3)) then
  ABI_DEALLOCATE(psik3)
end if
if(allocated(psik4)) then
  ABI_DEALLOCATE(psik4)
end if
if(allocated(psikb1)) then
  ABI_DEALLOCATE(psikb1)
end if
if(allocated(psikb2)) then
  ABI_DEALLOCATE(psikb2)
end if
if(allocated(psikb3)) then
  ABI_DEALLOCATE(psikb3)
end if
if(allocated(psikb4)) then
  ABI_DEALLOCATE(psikb4)
end if
if(allocated(psig1)) then
  ABI_DEALLOCATE(psig1)
end if
if(allocated(psig2)) then
  ABI_DEALLOCATE(psig2)
end if
if(allocated(psig3)) then
  ABI_DEALLOCATE(psig3)
end if
if(allocated(psig4)) then
  ABI_DEALLOCATE(psig4)
end if
if(allocated(psir1)) then
  ABI_DEALLOCATE(psir1)
end if
if(allocated(psir2)) then
  ABI_DEALLOCATE(psir2)
end if
if(allocated(psir3)) then
  ABI_DEALLOCATE(psir3)
end if
if(allocated(psidg)) then
  ABI_DEALLOCATE(psidg)
end if
if(allocated(vxc_dg)) then
  ABI_DEALLOCATE(vxc_dg)
end if
if(allocated(denpot)) then
  ABI_DEALLOCATE(denpot)
end if
if(allocated(kernel_wavefunctions_FFT)) then
  ABI_DEALLOCATE(kernel_wavefunctions_FFT)
end if
if(allocated(valence_wavefunctions_FFT)) then
  ABI_DEALLOCATE(valence_wavefunctions_FFT)
end if
if(associated(gvec)) then
  ABI_DEALLOCATE(gvec)
end if
if(allocated(vc_sqrt)) then
  ABI_DEALLOCATE(vc_sqrt)
end if

end subroutine destroy_H
!!***

!!****f* m_gwls_hamiltonian/build_H
!! NAME
!!  build_H
!!
!! FUNCTION
!!  .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine build_H(dtset2,mpi_enreg2,cpopt2,cg2,gs_hamk2,kg_k2,kinpw2)

!use m_bandfft_kpt
use m_cgtools
use m_wfutils

!Arguments of gw_sternheimer, reveived as argument by build_H-------------------------
type(dataset_type),  intent(in) :: dtset2
type(MPI_type),   intent(inout) :: mpi_enreg2
type(gs_hamiltonian_type), intent(inout) :: gs_hamk2

integer, intent(in) :: cpopt2
integer, intent(in) :: kg_k2(3,gs_hamk2%npw_k)
!Since mcg = npw_k*nband when there is only gamma, the next line works. If there is not only Gamma, mcg is the proper size of cg.
real(dp), intent(in) :: cg2(2,dtset2%mpw*dtset2%nspinor*dtset2%mband*dtset2%mkmem*dtset2%nsppol)
real(dp), intent(in) :: kinpw2(gs_hamk2%npw_k)
!Local variables : most of them are in the module for now.
real(dp) :: commutation_error
integer  :: io_unit
integer  :: io_unit_debug
integer  :: ierr
integer  :: cplx
integer  :: j, k
integer  :: dimph3d
integer  :: mb

integer  :: mpi_communicator

character(128):: filename_debug


real(dp), allocatable :: wfk_tmp1(:,:) ,wfk_tmp2(:,:)

! *************************************************************************

! Hartree-Fock cannot be used with GWLS.
if(dtset2%usefock==1 .and. associated(gs_hamk2%fockcommon)) then
  MSG_ERROR(' build_H : Hartree-Fock option can not be used with optdriver==66 (GWLS calculations).')
end if

!First we copy the data structure types
dtset = dtset2%copy()

call copy_mpi_enreg(mpi_enreg2,mpi_enreg)

call copy_hamiltonian(gs_hamk,gs_hamk2)

!Then we copy the standard types
cpopt   = cpopt2
npw_k = gs_hamk2%npw_k

ABI_ALLOCATE(cg,(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
cg = cg2

ABI_ALLOCATE(vlocal,(gs_hamk2%n4,gs_hamk2%n5,gs_hamk2%n6,gs_hamk2%nvloc))
vlocal = gs_hamk2%vlocal
gs_hamk%vlocal => vlocal

ABI_ALLOCATE(kg_k,(3,npw_k))
kg_k = kg_k2
ABI_ALLOCATE(kinpw,(npw_k))
kinpw = kinpw2

dimffnl=0; if (blocksize<=1) dimffnl=size(gs_hamk2%ffnl_k,2)
ABI_ALLOCATE(ffnl,(npw_k,dimffnl,gs_hamk%lmnmax,gs_hamk%ntypat))
dimph3d=0; if (blocksize<=1) dimph3d=gs_hamk2%matblk
ABI_ALLOCATE(ph3d,(2,npw_k,dimph3d))

!Initializing variables from dataset
nfft  =  dtset%nfft
nline = dtset%nline
tolwfr = dtset%tolwfr
ngfft  = dtset%ngfft
mcg = dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
nspinor=dtset%nspinor
n1=ngfft(1)
n2=ngfft(2)
n3=ngfft(3)
n4=ngfft(4)
n5=ngfft(5)
n6=ngfft(6)
mgfft=maxval(ngfft(1:3))
nbandv = int(dtset%nelect)/2
ABI_ALLOCATE(istwfk,(dtset%nkpt))
istwfk(:)=dtset%istwfk

!Initializing variables from gs_hamk
ucvol = gs_hamk%ucvol
ABI_ALLOCATE(gbound,(2*mgfft+8,2))
!gbound = gs_hamk%gbound_k !Must be done later for bandft paralelism

!Parameters which need to be set by hand for now...
weight           = 1           ! The weight of the k-pts, which sum to 1.
tim_fourwf       = 0
ckpt             = 1           ! The kpt of the wf to FFT. We only consider gamma point for now.

! ndat = 1 used to be hard coded, which works fine for FFT only parallelism.
ndat              = 1                                                ! # of FFT to do in parallel in fourwf

eshift           = 0.0         ! For PAW, opt_sij=-1 gives <G|H-lambda.S|C> in ghc instead of <G|H|C>
sij_opt          = 0           ! Option in PAW to tell getghc what to compute...
tim_getghc       = 0           ! timing code of the calling subroutine(can be set to 0 if not attributed)
type_calc        = 0           ! option governing which part of Hamitonian is to be applied: 0: whole Hamiltonian
nband            = dtset%mband
ispden           = 1           !When required, the spin index to be used. We don't support spin polarized calculations yet...



!========================================================================================================================
! Trying to implement band parallelism (Bruno, 14/10/2013)
!------------------------------------------------------------------------------------------------------------------------
!
! Band parallelism is described in the article on the LOBPCG method by F. Bottin et al,
! "Large scale ab initio calculations based on three levels of parallelisation", Computational Materials Science 42 (2008) 329-336
!
! The article describes the ABINIT implementation specifically, suggesting that the information there is directly related
! to the code.
!
! The main points in order to understand band + FFT parallelism are:
!
!     1) The processors are arranged in a 2D Cartesian MPI topology, or dimensions M x P, with
!                              M  = nproc_bands
!                              P  = nproc_fft
!        There are of course M x P processors
!
!     2) the wavefunction coefficients (namely the array cg) are distributed on these processors. Unfortunately,
!        two different distribution schemes are necessary: one for linear algebra (matrix products, etc), and one
!        consistent with the Goedecker FFTs. The code must thus change the data distribution back and forth.
!
!     3) the "linear algebra" distribution describes the array cg directly. The G vectors are distributed among  all
!        M x P processors, such that npw = npw_tot / (nproc_bands x nproc_fft). Every processor has information about
!        all bands. Schematically:
!
!
!                                        <- P = nproc_fft ->
!                             ----------------------------------
!                             |  n    |  n     |  n    |  n    |       ( for all band indices n)
!                     ^       |  G_11 | ...    | ...   |  G_P1 |
!                     |       ----------------------------------
!             M = nproc_bands |  n    |  n     |  n    |  n    |
!                     |       |  G_1. | ...    | ...   |  G_P. |
!                     v       ----------------------------------
!                             |  n    |  n     |  n    |  n    |
!                             |  G_1M | ...    | ...   |  G_PM |
!                             ----------------------------------
!
!
!        the LOBPCG algorithm acts on blocks of bands. Thus, although each processor has information about all bands,
!        we must conceptually view the bands grouped in blocks of size M, as {c_{1}(G),...,c_{M}(G)},
!        {c_{M+1}(G)...,c_{2M}(G)}, ...
!
!        With this conceptual grouping, for a given block each processor has M x NG/(M x P) coefficients, where NG is the
!        total number of G-vectors.
!
!     4) The distribution of information above is not appropriate for parallel FFTs. The data for one block is thus
!        transposed and takes the form
!
!                                        <- P = nproc_fft ->
!                             ----------------------------------
!                             |  1    |  1     | ...   |  1    |
!                     ^       |  G_1  |  G_2   | ...   |  G_P  |
!                     |       ----------------------------------
!             M = nproc_bands |  ..   | ...    | ...   |  ...  |
!                     |       |  G_1  |  G_2   | ...   |  G_P  |
!                     v       ----------------------------------
!                             |  M    |  M     | ...   |  M    |
!                             |  G_1  |  G_2   | ...   |  G_P  |
!                             ----------------------------------
!
!     where the set of G vectors G_i = (G_{i1}, G_{i2}, ..., G_{iM}). Thus, a given processor has information about a single
!     band, but more G vectors. Each processor has
!                             NG/(M x P) x M = NG/P coefficients, just like before.
!     Also, it is clear that the information is only communicated in the columns of the diagram above, avoiding
!     "all processors to all processors" communication.
!
!     The data is now distributed properly to do parallel FFTs! Each row in the diagram above corresponds to FFTs done
!     in parallel over nproc_fft processors, on a given band. There are M rows running thus in parallel!
!
!     The underlying ABINIT routines are equiped to handle FFT parallelism, not band distributed parallelism.
!     prep_getghc.F90 and lobpgcwf.F90 show how to rearange information to be able to apply basic ABINIT routines.
!
!     The code below is inspired / guessed from lobpcgwf and prep_getghc. WE ASSUME THERE IS NO SPINORS!  CODE SHOULD
!     BE HEAVILY REVIEWED for k-points and spinors, if and when they come to be useful.
!
!========================================================================================================================

! variables for band parallelism, following the language in lobpcgwf
nspinor          = 1                                      ! how many spinors are present. Hardcoded to 1 for now.
blocksize        = mpi_enreg%nproc_band*mpi_enreg%bandpp  ! how many bands treated in a block; # of FFT done in parallel in getghc
nbdblock         = nband/blocksize                        ! how many blocks of bands are there
ikpt_this_proc   = 1                                      ! Assuming only 1 kpoint, for molecules.
! This will have to be reviewed for crystals!

npw_kb           = npw_k*blocksize
npw_g            = gs_hamk2%npw_fft_k

if (blocksize>1) then
  kg_k_gather  => bandfft_kpt(ikpt_this_proc)%kg_k_gather
  ffnl_gather  => bandfft_kpt(ikpt_this_proc)%ffnl_gather
  ph3d_gather  => bandfft_kpt(ikpt_this_proc)%ph3d_gather
  kinpw_gather => bandfft_kpt(ikpt_this_proc)%kinpw_gather
else
  ffnl  = gs_hamk2%ffnl_k
  ph3d  = gs_hamk2%ph3d_k
  kg_k_gather  => kg_k
  ffnl_gather  => ffnl
  ph3d_gather  => ph3d
  kinpw_gather => kinpw
endif
call gs_hamk%load_k(kinpw_k=kinpw_gather,kg_k=kg_k_gather,ffnl_k=ffnl_gather,ph3d_k=ph3d_gather)

gbound = gs_hamk%gbound_k

!Set up wfk routines
call set_wf(ucvol/nfft,n4,n5,n6,npw_k,blocksize,npw_g,mpi_enreg%comm_bandfft,mpi_enreg%comm_fft,mpi_enreg%comm_band)

! prepare the valence wavefunctions and projection operator to work in band+FFT parallel
call DistributeValenceWavefunctions()
call DistributeValenceKernel()

cplx             = 2                                      ! wavefunctions have complex coefficients

! Initialize the dummy variables
ABI_DATATYPE_ALLOCATE(conjgrprj,(0,0))
ABI_ALLOCATE(dummy2,(0,0))
ABI_ALLOCATE(dummy3,(0,0,0))

!Initialisation of the total counter for iterations of SQMR :
ktot = 0

!Allocations of working arrays for the module
!- for pc_k function
ABI_ALLOCATE(scprod2,(2,nband))

!- for precondition function (and set_precondition subroutine)
ABI_ALLOCATE(pcon,(npw_g))
pcon = one

!- for (private) working wf
ABI_ALLOCATE(psik1,(2,npw_k))
ABI_ALLOCATE(psik2,(2,npw_k))
ABI_ALLOCATE(psik3,(2,npw_k))
ABI_ALLOCATE(psik4,(2,npw_k))
ABI_ALLOCATE(psikb1,(2,npw_kb))
ABI_ALLOCATE(psikb2,(2,npw_kb))
ABI_ALLOCATE(psikb3,(2,npw_kb))
ABI_ALLOCATE(psikb4,(2,npw_kb))
ABI_ALLOCATE(psig1,(2,npw_g))
ABI_ALLOCATE(psig2,(2,npw_g))
ABI_ALLOCATE(psig3,(2,npw_g))
ABI_ALLOCATE(psig4,(2,npw_g))
ABI_ALLOCATE(psir1,(2,n4,n5,n6))
ABI_ALLOCATE(psir2,(2,n4,n5,n6))
ABI_ALLOCATE(psir3,(2,n4,n5,n6))
ABI_ALLOCATE(psidg,(2,nfft))
ABI_ALLOCATE(denpot,(2*n4,n5,n6))

psir1 = zero
psir2 = zero
psir3 = zero
denpot = zero

!Construct the vector of eigenvalues and write then to std output
ABI_ALLOCATE(eig,(nband))
eig=zero

write(std_out,*) ch10,"Eigenvalues computation check, routine build_H:",ch10
io_unit_debug = get_unit()
write(filename_debug,'(A,I0.4,A)') "DEBUG_PROC=",mpi_enreg%me,".dat"



mpi_communicator =  mpi_enreg%comm_bandfft

open(file=filename_debug,status=files_status_old,unit=io_unit_debug)

write(io_unit_debug,'(A)') " Parameters:"
write(io_unit_debug,'(A,I5)') "                 nband = ", nband
write(io_unit_debug,'(A,I5)') "             blocksize = ", blocksize
write(io_unit_debug,'(A,I5)') "                 npw_k = ", npw_k
write(io_unit_debug,'(A,I5)') "              nbdblock = ", nbdblock

! temporary wavefunction array, for data in the "linear algebra" distribution
ABI_ALLOCATE( wfk_tmp1, (2,npw_k))
ABI_ALLOCATE( wfk_tmp2, (2,npw_k))

do n=1, nband
! Extract i^t/h wavefunction
wfk_tmp1(:,1:npw_k) = cg(:,(n-1)*npw_k+1:n*npw_k)

! Are the wavefunctions normalized?
!tmpc = cg_zdotc(npw_k,wfk_tmp1,wfk_tmp1)
!call xmpi_sum(tmpc,mpi_enreg%comm_bandfft,ierr) ! sum on all processors

! DEBUGGING CODE
write(std_out,'(A,I5,A,2F24.16)') "band ", n, ", <psi_n | psi_n > =", norm_k(wfk_tmp1)
write(io_unit_debug,'(A,I5,A,2F24.16)') "band ", n, ", <psi_n | psi_n > =", norm_k(wfk_tmp1)
flush(io_unit_debug)
end do


! loop on blocks of bands. All bands in a given block will be done in parallel!

! loop on block of bands
do iblock = 1, nbdblock

! loop on states for this block
do iband = 1, blocksize

n = (iblock-1)*blocksize+iband

psikb1(:,(iband-1)*npw_k+1:iband*npw_k)   = cg(:,(n-1)*npw_k+1:n*npw_k)

end do

! change configuration of the data
call wf_block_distribute(psikb1,  psig1,1) ! LA -> FFT

! Apply hamiltonian on wavefunction. In bandFFT parallel, Hpsik calls prep_getghc, which knows how to transform the
! distribution of the data from "linear algebra"-like to "FFT"-like. The output is *probably* in the "linear algebra"-like
! distribution.
call Hpsik(psig2,psig1)

! change configuration of the data
call wf_block_distribute(psikb2,  psig2,2) ! LA -> FFT

! extract the wavefunctions band by band
do iband=1, blocksize

n = blocksize*(iblock-1) + iband ! band index

wfk_tmp1(:,1:npw_k) = psikb1(:,(iband-1)*npw_k+1:iband*npw_k)
wfk_tmp2(:,1:npw_k) = psikb2(:,(iband-1)*npw_k+1:iband*npw_k)

tmpc   = cg_zdotc(npw_k,wfk_tmp1,wfk_tmp2)

call xmpi_sum(tmpc,mpi_enreg%comm_bandfft,ierr) ! sum on all processors
eig(n) = tmpc(1)

write(std_out,'(A,I5,A,F24.16,A)') "(build_H) band ", n, ", eig =", eig(n), " Ha."

! DEBUGGING CODE
write(io_unit_debug,'(A,I5,A,F24.16,A)') "band ", n, ", eig =", eig(n), " Ha."
flush(io_unit_debug)
flush(io_unit_debug)
end do

end do

ABI_DEALLOCATE( wfk_tmp1 )
ABI_DEALLOCATE( wfk_tmp2 )
close(io_unit_debug)



! Finishing the construction of vxc (transcribing it from the double real grid (for the density)
! to the single real grid (for the wfs).
! Assummes we only need one spin component; one transcription per spin being needed.
if(allocated(vxc_dg)) then
  ABI_ALLOCATE(vxc,(n4,n5,n6,dtset%nspden))
  vxc = zero
  call fftpac(ispden,mpi_enreg,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vxc_dg(:,ispden),vxc(:,:,:,ispden),2)
end if

timrev = 1 !Assumes istwfk *1. 1:time reversal symmetry NOT present | 2:" " " present
ABI_ALLOCATE(title,(dtset%ntypat))
do i=1,dtset%ntypat
title(i) = "Bloup" ! The clean way would be to get the psps structure in this module
! (build_vxc is called from a place in GS calculations where it is available;
! should be the easiest way). For now, this allows the code to run.
end do
call crystal_init(dtset%amu_orig(:,1),Cryst,dtset%spgroup,dtset%natom,dtset%npsp,&
&                 dtset%ntypat,dtset%nsym,dtset%rprimd_orig(:,:,1),dtset%typat,&
&                 dtset%xred_orig(:,:,1),dtset%ziontypat,dtset%znucl,timrev,.false.,.false.,title,&
&                 dtset%symrel,dtset%tnons,dtset%symafm)
ABI_DEALLOCATE(title)
call Cryst%print()

!TODO : Should be put in a separate build_vc constructor, and should be called right after build_H in the context of optdriver 66.
if(dtset%optdriver==66) then

  !Set up of the k-points and tables in the whole BZ
  call kmesh_init(Kmesh,Cryst,dtset%nkpt,dtset%kptns,dtset%kptopt,wrap_1zone=.false.)
  call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out)
  call find_qmesh(Qmesh,Cryst,Kmesh)
  call kmesh_print(Qmesh,"Q-mesh for the screening function",std_out)


  !------------------------------
  !Building the vc_sqrt structure
  !------------------------------
  ! The sqrt_vc code relies on a different parallelism scheme than Vanilla ABINIT.
  ! The Gsphere object must be built on a single processor having all the G-vectors
  ! available.

  !gsph_init need gvec built according to the KSS convention, so that kg_k is not suitable (not properly sorted).
  npw_serial=npw_k
  call xmpi_sum(npw_serial,mpi_enreg%comm_bandfft,ierr)

  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  call make_gvec_kss(dtset%nkpt,dtset%kptns,ecut_eff,dtset%symmorphi,dtset%nsym,dtset%symrel,dtset%tnons,Cryst%gprimd,&
  &                      dtset%prtvol,npw_serial,gvec,ierr)

  call gsph_init(Gsphere,Cryst,npw_serial,gvec=gvec)

  call print_gsphere(Gsphere)

  call vcoul_init(Vcp,Gsphere,Cryst,Qmesh,Kmesh,dtset%rcut,dtset%icutcoul,dtset%vcutgeo,dtset%ecutsigx,npw_serial,&
  &               dtset%nkpt,dtset%kptns,dtset%ngfft,mpi_enreg%comm_world)

  ! Since Vcp%vc_sqrt is sorted according to the KSS convention for G vectors
  ! BUT will be applied to GS wavefunctions (where G vectors are sorted otherwise)
  ! It is necessary to construct a vc_sqrt vector with the GS sorting.
  ! Moreover, if we are in parallel (over band / FFTs), Vcp%vc_sqrt is the serial version
  ! and need to be distributed according to the GS scheme (Linear Algebra configuration).
  ABI_ALLOCATE(vc_sqrt,(npw_k))
  vc_sqrt=zero
  k=0
  do i=1,npw_k
  do j=1,npw_serial
  if(all(kg_k(:,i)==gvec(:,j))) k=j
  end do
  vc_sqrt(i)=Vcp%vc_sqrt(k,1)
  end do

end if

!--------------------------------------------------------------------------------
! Security check : The eigenstates of the projector and the hamiltonian need to
! agree down to the precision requested (tolwfr). Otherwise, SQMR is doomed.
! Now that the density is read and the eigenstates are calculated in the GW run,
! it is less useful. A test on tolwfr could be sufficient.
!
! This remains a good tool to debug band+fft parallelism
!--------------------------------------------------------------------------------

! only write on the head node!
if (mpi_enreg%me == 0) then
  io_unit = get_unit()
  open(file='build_H.log',status=files_status_old,unit=io_unit)
  write(io_unit,10) '#----------------------------------------------------------------------------'
  write(io_unit,10) '#'
  write(io_unit,10) '#               This file presents the results of a small test to check      '
  write(io_unit,10) '#               how well the Hamiltonian commutes with the projection        '
  write(io_unit,10) '#               operator.                                                    '
  write(io_unit,10) '#'
  write(io_unit,10) '# Definitions:'
  write(io_unit,10) '#'
  write(io_unit,10) '#                    P     : projections on conduction states'
  write(io_unit,10) '#                    H     : Hamiltonian operator'
  write(io_unit,10) '#                  | psi > : eigenstate'
  write(io_unit,10) '#'
  flush(io_unit)
end if

psikb1 = zero

! sum all valence states, on copy in every band block.
do i=1,nbandv

do mb = 1, blocksize

psikb1(:,(mb-1)*npw_k+1:mb*npw_k) = psikb1(:,(mb-1)*npw_k+1:mb*npw_k)  + cg(:,(i-1)*npw_k+1:i*npw_k)

end do

end do

! normalize to one! Only sum the fist band block
tmpc = cg_zdotc(npw_k,psikb1,psikb1)
call xmpi_sum(tmpc,mpi_communicator ,ierr) ! sum on all processors
psikb1 = psikb1/sqrt(tmpc(1))


! change data distribution
call wf_block_distribute(psikb1,  psig1,1) ! LA -> FFT

! Apply P.H operator
call Hpsik(psig2 ,psig1)
call pc_k_valence_kernel(psig2)

! Apply H.P operator
call pc_k_valence_kernel(psig1)
call Hpsik(psig1)

! compute error
psig3 = psig1 - psig2 ! explicitly write the difference in an array
tmpc   = cg_zdotc(npw_g,psig3,psig3)
commutation_error = tmpc(1)

call xmpi_sum(commutation_error , mpi_enreg%comm_fft, ierr) ! sum on all processors working on FFT!

! only write on the head node!
if (mpi_enreg%me == 0) then
  write(io_unit,20) '   || (PH -HP) |b> ||^2 =  ',commutation_error
  write(io_unit,20) '           tolwfr       =  ',tolwfr
end if

if(commutation_error > tolwfr) then
  !write(io_unit,10) '# || (PH -HP) |b> ||^2 > tolwfr ==> Decision taken exit!'

  if (mpi_enreg%me == 0) write(io_unit,10) '# || (PH -HP) |b> ||^2 > tolwfr ==> This must be fixed!'

  write(std_out,20) "WARNING-build_H: The tolerance tolwfr=",tolwfr
  write(std_out,20) "                 is smaller than the commutation error ||(PH-HP)|b>||^2=",commutation_error
  write(std_out,10) "                 Either set tolwfr to a less stringent value in this calculation"
  write(std_out,10) "                 or to a more stringent value in the wavefunction calculation."
end if

if (mpi_enreg%me == 0) close(io_unit)

10 format(A)
20 format(A,ES12.3)
end subroutine build_H
!!***

end module m_gwls_hamiltonian
!!***
