!!****f* ABINIT/m_lobpcgwf_cprj
!! NAME
!! m_lobpcgwf_cprj
!!
!! FUNCTION
!! This routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2026 ABINIT group (LB)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lobpcgwf_cprj

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_lobpcg
 use m_errors
 use m_time
 use m_xomp
 use m_fstrings
 use m_xg
 use m_xg_nonlop
 use m_xgTransposer
 use m_lobpcg2_cprj
 use m_dtset

 use defs_abitypes, only : mpi_type
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_getghc,      only : multithreaded_getghc
 use m_pawcprj,     only : pawcprj_type

 implicit none

 private

 integer, parameter :: l_tim_getghc=5

! For use in getghc_gsc1
 integer, save :: l_prtvol
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: lobpcgwf2_cprj

 contains
!!***

subroutine lobpcgwf2_cprj(cg,dtset,eig,occ,enl_out,gs_hamk,isppol,ikpt,inonsc,istep,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid,nbdbuf,xg_nonlop)


 use m_cgtools, only : dotprod_g

!Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 integer,intent(in) :: isppol,ikpt,inonsc,istep,nbdbuf
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk
 type(dataset_type)              ,intent(in   ) :: dtset
 type(mpi_type)           ,target,intent(in)    :: mpi_enreg
 real(dp)                 ,target,intent(inout) :: cg(2,nspinor*nband*npw)
 real(dp)                        ,intent(in   ) :: kinpw(npw)
 real(dp)                 ,target,intent(  out) :: resid(nband)
 real(dp)                        ,intent(  out) :: enl_out(nband)
 real(dp)                 ,target,intent(  out) :: eig(nband)
 real(dp)                 ,target,intent(in   ) :: occ(nband)
 type(xg_nonlop_t)        ,target,intent(in   ) :: xg_nonlop

!Local variables-------------------------------

 type(xgBlock_t) :: xgx0
 !type(xgBlock_t) :: cprj_xgx0
 type(xg_t) :: cprj_xgx0
 type(xgBlock_t) :: xgeigen
 type(xgBlock_t) :: xgresidu
 type(xgBlock_t) :: xgocc
 type(xgBlock_t) :: xgenl
 type(xgBlock_t) :: xg_precond,xg_kin
 type(lobpcg_t) :: lobpcg

 integer :: space, space_cprj, blockdim, cprjdim, nband_cprj
 integer :: me_g0, me_g0_fft

 integer, parameter :: tim_lobpcgwf2 = 2030
 real(dp) :: tsec(2)

 real(dp), allocatable :: occ_tmp(:)
 ! Important things for NC
 real(dp), allocatable :: pcon(:),kin(:)
! real(dp), allocatable :: cprj_contiguous(:,:)

! *********************************************************************

 call timab(tim_lobpcgwf2,1,tsec)

 ! Set module variables
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk

 cprjdim = xg_nonlop%cprjdim

!Variables
 blockdim=nband/dtset%nblock_lobpcg
 if (blockdim/=mpi_enreg%nproc_band*mpi_enreg%bandpp) then
   ABI_ERROR('blockdim is not consistent with nproc_band and bandpp')
 end if
 nband_cprj=nband/mpi_enreg%nproc_band

!Depends on istwfk
 if ( gs_hamk%istwf_k > 1 ) then ! Real only
   ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
   ! and im*im so that a vector of complex is consider as a long vector of real
   ! therefore the number of data is (2*npw*nspinor)*nband
   ! This space is completely equivalent to SPACE_R but will correctly set and
   ! get the array data into the xgBlock
   space      = SPACE_CR
 else ! complex
   space      = SPACE_C
 end if
 space_cprj = xg_nonlop%space_cprj

 !For kinetic part of the Hamiltonian
 ABI_MALLOC(kin,(npw))
 call build_kin(kin,kinpw,npw)
 call xgBlock_map_1d(xg_kin,kin,SPACE_R,npw)

 !For preconditionning
 ABI_MALLOC(pcon,(npw))
 call build_pcon(pcon,kinpw,npw)

 ! Local variables for lobpcg
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

 call xgBlock_map_1d(xgenl,enl_out,SPACE_R,nband)

 call xg_init(cprj_xgx0,space_cprj,xg_nonlop%cprjdim,nband_cprj*nspinor,comm=l_mpi_enreg%comm_band)

 ! Occupancies in LOBPCG are used for convergence criteria only
 if (dtset%nbdbuf==-101.and.nspinor==1.and.dtset%nsppol==1) then
   ABI_MALLOC(occ_tmp,(nband))
   occ_tmp(:) = half*occ(:)
   call xgBlock_map_1d(xgocc,occ_tmp,SPACE_R,nband,gpu_option=dtset%gpu_option)
 else
   call xgBlock_map_1d(xgocc,occ,SPACE_R,nband,gpu_option=dtset%gpu_option)
 end if

 call lobpcg_init(lobpcg,mpi_enreg%bandpp,nband,npw*nspinor,cprjdim,blockdim,dtset%tolwfr_diago,dtset%nline,&
   space,space_cprj,l_mpi_enreg%comm_band,dtset%paral_kgb,xg_nonlop,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
   me_g0,me_g0_fft)

 ! Run lobpcg
 call lobpcg_run_cprj(lobpcg,xgx0,cprj_xgx0%self,xg_getghc,xg_kin,xg_precond,xgeigen,xgocc,xgresidu,xgenl,&
   prtvol,nspinor,isppol,ikpt,inonsc,istep,nbdbuf)

 if (allocated(occ_tmp)) then
   ABI_FREE(occ_tmp)
 end if
 ! Free preconditionning since not needed anymore
 ABI_FREE(pcon)
 ABI_FREE(kin)

 call xg_free(cprj_xgx0)

 ! Free lobpcg
 call lobpcg_free(lobpcg)

 call timab(tim_lobpcgwf2,2,tsec)

 DBG_EXIT("COLL")

end subroutine lobpcgwf2_cprj
!!***

!!****f* m_lobpcg/xg_getghc
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
!! SOURCE
!
subroutine xg_getghc(X,AX)

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

end module m_lobpcgwf_cprj
!!***
