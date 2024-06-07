!!****f* ABINIT/m_chebfiwf_cprj
!! NAME
!! m_chebfiwf_cprj
!!
!! FUNCTION
!! This module contains a routine updating the whole wave functions at a given k-point,
!! using the Chebyshev filtering method (2021 implementation using xG abstraction layer)
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2023-2023 ABINIT group (LB)
!! This file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_chebfiwf_cprj

 use defs_abitypes
 use defs_basis
 use m_abicore
 use m_errors
 use m_fstrings
 use m_time
 use m_xg
 use m_xg_nonlop
 use m_chebfi2_cprj
 use m_dtset

 use m_chebfi
 use m_invovl

 use m_dtset,       only : dataset_type

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_getghc,      only : multithreaded_getghc

 use m_xg
 use m_xgTransposer

 use iso_c_binding, only: c_associated,c_loc,c_ptr,c_f_pointer

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 implicit none

 private

 integer, parameter :: l_tim_getghc=5

 integer, parameter :: CPRJ_ALLOC=1
 integer, parameter :: CPRJ_FREE=2

! For use in getghc_gsc1
 integer, save :: l_prtvol
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: chebfiwf2_cprj

 contains

subroutine chebfiwf2_cprj(cg,cprj_cwavef_bands,dtset,eig,enl_out,gs_hamk,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid,xg_nonlop)

!!****f* m_chebfiwf/chebfiwf2
!! NAME
!! chebfiwf2
!!
!! FUNCTION
!! This routine updates the whole wave functions set at a given k-point,
!! using the Chebfi method (2021 version using xG abstraction layer)
!!
!! INPUTS
!!  dtset= input variables for this dataset
!!  kinpw(npw)= kinetic energy for each plane wave (Hartree)
!!  mpi_enreg= MPI-parallelisation information
!!  nband= number of bands at this k point
!!  npw= number of plane waves at this k point
!!  nspinor= number of spinorial components of the wavefunctions
!!  prtvol= control print volume and debugging
!!
!! OUTPUT
!!  eig(nband)= eigenvalues (hartree) for all bands
!!  enl_out(nband)= contribution of each band to the nl part of energy
!!  resid(nband)= residuals for each band
!!
!! SIDE EFFECTS
!!  cg(2,npw*nspinor*nband)= planewave coefficients of wavefunctions
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!
!! SOURCE

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk
 type(dataset_type)              ,intent(in   ) :: dtset
 type(mpi_type)           ,target,intent(in)    :: mpi_enreg
 type(pawcprj_type)       ,target,intent(inout) :: cprj_cwavef_bands(:,:)
 real(dp)                 ,target,intent(inout) :: cg(2,nspinor*nband*npw)
 real(dp)                        ,intent(in   ) :: kinpw(npw)
 real(dp)                 ,target,intent(  out) :: resid(nband)
 real(dp)                        ,intent(  out) :: enl_out(nband)
 real(dp)                 ,target,intent(  out) :: eig(nband)
 type(xg_nonlop_t)        ,       intent(in   ) :: xg_nonlop

!Local variables-------------------------------

 type(xgBlock_t) :: xgx0
 type(xgBlock_t) :: cprj_xgx0
 type(xgBlock_t) :: xgeigen
 type(xgBlock_t) :: xgresidu
 type(xgBlock_t) :: xg_precond,xg_kin
 type(chebfi_t) :: chebfi

 logical :: paw

 integer :: space, space_cprj, blockdim, cprjdim, nband_cprj
 integer :: me_g0,me_g0_fft

 integer, parameter :: tim_chebfiwf2 = 2060
 double precision :: tsec(2)

 ! Important things for NC
 real(dp), allocatable :: pcon(:),kin(:)
 real(dp), allocatable :: cprj_contiguous(:,:)

! *********************************************************************

 call timab(tim_chebfiwf2,1,tsec)

 ! Set module variables
 !LTEST
 if (cprj_cwavef_bands(1,1)%ncpgr==3) then
   ABI_ERROR('chebfi with cprj not implemented with cprj%ncpgr==3')
 end if
 !LTEST

 paw = gs_hamk%usepaw==1

 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk

 cprjdim = xg_nonlop%cprjdim

!Variables
 blockdim=mpi_enreg%nproc_band*mpi_enreg%bandpp
 nband_cprj=nband/mpi_enreg%nproc_band

!Depends on istwfk
 if ( gs_hamk%istwf_k > 1 ) then ! Real only
   ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
   ! and im*im so that a vector of complex is consider as a long vector of real
   ! therefore the number of data is (2*npw*nspinor)*nband
   ! This space is completely equivalent to SPACE_R but will correctly set and
   ! get the array data into the xgBlock
   space      = SPACE_CR
   space_cprj = SPACE_R
 else ! complex
   space      = SPACE_C
   space_cprj = SPACE_C
 end if

 !For kinetic part of the Hamiltonian
 ABI_MALLOC(kin,(gs_hamk%npw_fft_k))
 call build_kin(kin,gs_hamk%kinpw_k,gs_hamk%npw_fft_k)
 call xgBlock_map_1d(xg_kin,kin,SPACE_R,gs_hamk%npw_fft_k)

 !For preconditionning
 ABI_MALLOC(pcon,(npw))
 call build_pcon(pcon,kinpw,npw)

 ! Local variables for chebfi
 me_g0 = -1
 me_g0_fft = -1
 if (space==SPACE_CR) then
   me_g0 = 0
   me_g0_fft = 0
   if (gs_hamk%istwf_k == 2) then
     if (l_mpi_enreg%me_g0 == 1) me_g0 = 1
     if (l_mpi_enreg%me_g0_fft == 1) me_g0_fft = 1
   end if
 end if
 call xgBlock_map(xgx0,cg,space,npw*nspinor,nband,l_mpi_enreg%comm_band,me_g0=me_g0)

 call xgBlock_map_1d(xg_precond,pcon,SPACE_R,npw)

 call xgBlock_map_1d(xgeigen,eig,SPACE_R,nband)

 call xgBlock_map_1d(xgresidu,resid,SPACE_R,nband)

 call xg_cprj_copy(cprj_cwavef_bands,cprj_contiguous,space_cprj,nband_cprj,cprj_xgx0,&
   & xg_nonlop,l_mpi_enreg%comm_band,CPRJ_ALLOC)

 call chebfi_init(chebfi,nband,npw*nspinor,cprjdim,dtset%tolwfr_diago,dtset%ecut, &
&                 mpi_enreg%bandpp, &
&                 dtset%nline, space,space_cprj,1, &
&                 l_mpi_enreg%comm_band,me_g0,paw,&
&                 xg_nonlop,me_g0_fft)


 ! Run chebfi
 call chebfi_run_cprj(chebfi,xgx0,cprj_xgx0,xg_getghc,xg_kin,xg_precond,xgeigen,xgresidu,nspinor)

 ! Free preconditionning since not needed anymore
 ABI_FREE(pcon)
 ABI_FREE(kin)

 call xg_cprj_copy(cprj_cwavef_bands,cprj_contiguous,space_cprj,nband_cprj,cprj_xgx0,&
   & xg_nonlop,l_mpi_enreg%comm_band,CPRJ_FREE)

 ! Free chebfi
 call chebfi_free(chebfi)

 call timab(tim_chebfiwf2,2,tsec)

 DBG_EXIT("COLL")

end subroutine chebfiwf2_cprj

!!****f* m_chebfi/xg_getghc
!! NAME
!! xg_getghc
!!
!! FUNCTION
!! This routine computes H|C> and possibly S|C> for a given wave function C.
!!  It acts as a driver for getghc, taken into account parallelism, multithreading, etc.
!!
!! SIDE EFFECTS
!!  X  <type(xgBlock_t)>= memory block containing |C>
!!  AX <type(xgBlock_t)>= memory block containing H|C>
!!
!! PARENTS
!!
!! CHILDREN
!!      xgBlock_getSize,xgBlock_reverseMap,xgBlock_scale,xgBlock_copy
!!      multithreaded_getghc
!!
!! SOURCE
!
subroutine xg_getghc(X,AX)

 use iso_c_binding
 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: AX

!Local variables-------------------------------
!scalars
 integer         :: blockdim
 integer         :: spacedim
 integer,parameter :: sij_opt=0,cpopt=-1,type_calc=1 ! Compute local part only
! integer :: iatom,iband,ispinor,cprj_index,cprj_rows,cprj_cols,ncpgr,nlmn
 real(dp) :: eval
 type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,1)
!arrays
 real(dp), pointer :: cg(:,:)
 real(dp), pointer :: ghc(:,:)
 real(dp) :: gsc(1,1),gvnlxc(1,1)

! *********************************************************************

 call xgBlock_getSize(X,spacedim,blockdim)
 call xgBlock_check(X,AX)

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)

 ! Apply only local part of the Hamiltonian
 call multithreaded_getghc(cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,sij_opt,l_tim_getghc,type_calc)

end subroutine xg_getghc
!!***

subroutine build_pcon(pcon,kinpw,npw)

  implicit none

  integer,intent(in) :: npw
  real(dp),intent(in) :: kinpw(:)
  real(dp),intent(out) :: pcon(:)

  integer :: ipw

  !$omp parallel do schedule(static), shared(pcon,kinpw)
  do ipw=1,npw
    if(kinpw(ipw)>huge(0.0_dp)*1.d-11) then
      pcon(ipw)=0.d0
    else
      pcon(ipw) = (27+kinpw(ipw)*(18+kinpw(ipw)*(12+8*kinpw(ipw)))) &
&     / (27+kinpw(ipw)*(18+kinpw(ipw)*(12+8*kinpw(ipw))) + 16*kinpw(ipw)**4)
    end if
  end do

end subroutine build_pcon

subroutine build_kin(kin,kinpw,npw)

  implicit none

  integer,intent(in) :: npw
  real(dp),intent(in) :: kinpw(:)
  real(dp),intent(out) :: kin(:)

  integer :: ipw

  !$omp parallel do schedule(static), shared(kin,kinpw)
  do ipw=1,npw
    if(kinpw(ipw)>huge(0.0_dp)*1.d-11) then
      kin(ipw)=0.d0
    else
      kin(ipw) = kinpw(ipw)
    end if
  end do

end subroutine build_kin

subroutine xg_cprj_copy(cprj_cwavef_bands,cprj_contiguous,space_cprj,ncprj,xg_cprj,xg_nonlop,comm,option)

  implicit none

  integer, intent(in) :: comm,space_cprj,ncprj,option
  type(pawcprj_type),intent(inout)   :: cprj_cwavef_bands(:,:)
  real(dp),allocatable,intent(inout) :: cprj_contiguous(:,:)
  type(xgBlock_t), intent(inout) :: xg_cprj
  type(xg_nonlop_t), intent(in)  :: xg_nonlop

  integer :: cplex,cprj_index,iatom,iband,ispinor,nlmn,nspinor

  if (option/=CPRJ_ALLOC.and.option/=CPRJ_FREE) then
    ABI_ERROR('Bad option')
  end if

  cplex=2;if (space_cprj==SPACE_R) cplex=1

  nspinor = xg_nonlop%nspinor

  if (option==CPRJ_ALLOC) then
    ABI_MALLOC(cprj_contiguous,(cplex,xg_nonlop%cprjdim*ncprj*nspinor))
    call xgBlock_map(xg_cprj,cprj_contiguous,space_cprj,xg_nonlop%cprjdim,ncprj*nspinor,comm)
  end if

  cprj_index=1
  do iband=1,ncprj
    do iatom=1,xg_nonlop%natom
      nlmn=xg_nonlop%nlmn_natom(iatom)
      do ispinor=1,nspinor
        if (option==CPRJ_ALLOC) then
          cprj_contiguous(:,cprj_index:cprj_index+nlmn-1) = &
            cprj_cwavef_bands(iatom,(iband-1)*nspinor+ispinor)%cp(1:cplex,1:nlmn)
        else
          cprj_cwavef_bands(iatom,(iband-1)*nspinor+ispinor)%cp(1:cplex,1:nlmn) = &
            & cprj_contiguous(:,cprj_index:cprj_index+nlmn-1)
        end if
        cprj_index=cprj_index+nlmn
      end do
    end do
  end do

  if (option==CPRJ_FREE) then
    ABI_FREE(cprj_contiguous)
  end if

end subroutine xg_cprj_copy

end module m_chebfiwf_cprj
!!***
