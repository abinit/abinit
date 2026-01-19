!!****f* ABINIT/m_lobpcgwf
!! NAME
!! m_lobpcgwf
!!
!! FUNCTION
!! This routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2025 ABINIT group (JB)
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

! nvtx related macro definition
#include "nvtx_macros.h"

module m_lobpcgwf

 use defs_basis
 use m_abicore
 use m_lobpcg
 use m_xmpi
 use m_errors
 use m_time
 use m_xomp
 use m_fstrings
 use m_xg
 use m_xgTransposer
 use m_lobpcg2
 use m_dtset

 use defs_abitypes, only : mpi_type
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_getghc,prep_nonlop
 use m_getghc,      only : multithreaded_getghc

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

 use, intrinsic :: iso_c_binding

 implicit none
 private

 integer, parameter :: l_tim_getghc=5
 double precision, parameter :: inv_sqrt2 = 1/sqrt2

 ! For use in getghc_gsc1
 integer,save  :: l_cpopt
 logical,save  :: l_paw
 integer,save  :: l_prtvol
 integer,save  :: l_sij_opt
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: lobpcgwf2

 contains

subroutine lobpcgwf2(cg,dtset,eig,occ,enl_out,gs_hamk,isppol,ikpt,inonsc,istep,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid,nbdbuf)



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

!Local variables-------------------------------

 type(xgBlock_t) :: xgx0
 type(xgBlock_t) :: xgeigen
 type(xgBlock_t) :: xgresidu
 type(xgBlock_t) :: xgocc
 type(xgBlock_t) :: xg_precond
 type(lobpcg_t) :: lobpcg

 integer :: space, blockdim

 integer, parameter :: tim_lobpcgwf2 = 1640
 integer, parameter :: tim_enl = 1657
 double precision :: tsec(2)

 ! Important things for NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 type(pawcprj_type) :: cprj_dum(1,1)
 integer :: iblock, shift, me_g0, me_g0_fft
 real(dp) :: gsc_dummy(0,0)
 real(dp), allocatable :: gvnlxc(:,:)
 real(dp), allocatable :: pcon(:),occ_tmp(:)

! *********************************************************************

 call timab(tim_lobpcgwf2,1,tsec)

 ! Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk

!Variables
 blockdim=nband/dtset%nblock_lobpcg
 if (blockdim/=mpi_enreg%nproc_band*mpi_enreg%bandpp) then ! without this check computation of enl_out can be wrong
   ABI_ERROR('blockdim is not consistent with nproc_band and bandpp')
 end if

!Depends on istwfk
 if ( gs_hamk%istwf_k > 1 ) then ! Real only
   ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
   ! and im*im so that a vector of complex is consider as a long vector of real
   ! therefore the number of data is (2*npw*nspinor)*nband
   ! This space is completely equivalent to SPACE_R but will correctly set and
   ! get the array data into the xgBlock
   space = SPACE_CR
 else ! complex
   space = SPACE_C
 end if

 !For preconditionning
 ABI_MALLOC(pcon,(npw))
 call build_pcon(pcon,kinpw,npw)

#ifdef HAVE_OPENMP_OFFLOAD
 if(gs_hamk%gpu_option==ABI_GPU_OPENMP) then
   !$OMP TARGET ENTER DATA MAP(to:cg,eig,resid,occ,pcon)
 end if
#endif

 ! Local variables for lobpcg
 me_g0 = -1
 me_g0_fft = -1
 if (space==SPACE_CR) then
   me_g0 = 0
   me_g0_fft = 0
   if (gs_hamk%istwf_k == 2) then
     if (mpi_enreg%me_g0 == 1) me_g0 = 1
     if (mpi_enreg%me_g0_fft == 1) me_g0_fft = 1
   end if
 end if
 call xgBlock_map(xgx0,cg,space,npw*nspinor,nband,comm=mpi_enreg%comm_bandspinorfft,me_g0=me_g0,&
   & gpu_option=dtset%gpu_option)

 call xgBlock_map_1d(xg_precond,pcon,SPACE_R,npw,gpu_option=dtset%gpu_option)

 call xgBlock_map_1d(xgeigen,eig,SPACE_R,nband,gpu_option=dtset%gpu_option)

 call xgBlock_map_1d(xgresidu,resid,SPACE_R,nband,gpu_option=dtset%gpu_option)

 ! Occupancies in LOBPCG are used for convergence criteria only
 if (dtset%nbdbuf==-101.and.nspinor==1.and.dtset%nsppol==1) then
   ABI_MALLOC(occ_tmp,(nband))
   occ_tmp(:) = half*occ(:)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(to:occ_tmp) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
   call xgBlock_map_1d(xgocc,occ_tmp,SPACE_R,nband,gpu_option=dtset%gpu_option)
 else
   call xgBlock_map_1d(xgocc,occ,SPACE_R,nband,gpu_option=dtset%gpu_option)
 end if

 call lobpcg_init(lobpcg,nband,npw*nspinor,blockdim,dtset%tolwfr_diago,dtset%nline,&
   space,mpi_enreg%comm_bandspinorfft,dtset%paral_kgb,mpi_enreg%comm_spinorfft,mpi_enreg%comm_band,&
   me_g0,me_g0_fft,gs_hamk%gpu_option,gpu_thread_limit=dtset%gpu_thread_limit)

 ! Run lobpcg
 call lobpcg_run(lobpcg,xgx0,getghc_gsc1,xg_precond,xgeigen,xgocc,xgresidu,prtvol,nspinor,isppol,ikpt,inonsc,istep,nbdbuf)

 if (allocated(occ_tmp)) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET EXIT DATA MAP(delete:occ_tmp) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
   ABI_FREE(occ_tmp)
 end if
 ! Free preconditionning since not needed anymore
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET EXIT DATA MAP(delete:pcon) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
 ABI_FREE(pcon)

 if ( .not. l_paw ) then
   call timab(tim_enl,1,tsec)
#ifdef FC_CRAY
   ABI_MALLOC(gvnlxc,(1,1))
#else
   ABI_MALLOC(gvnlxc,(0,0))
#endif

   !Call nonlop
   if (dtset%paral_kgb==0) then

     call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,gvnlxc)

   else
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
     do iblock=1,nband/blockdim
       shift = (iblock-1)*blockdim*npw*nspinor
       call prep_nonlop(choice,l_cpopt,cprj_dum, &
&        enl_out((iblock-1)*blockdim+1:iblock*blockdim),gs_hamk,0,&
&        eig((iblock-1)*blockdim+1:iblock*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&        gsc_dummy,l_tim_getghc,cg(:,shift+1:shift+blockdim*npw*nspinor),gvnlxc(:,:),&
&        already_transposed=.false.)
     end do
   end if
   ABI_FREE(gvnlxc)
   call timab(tim_enl,2,tsec)
 end if

 ! Free lobpcg
 call lobpcg_free(lobpcg)

#ifdef HAVE_OPENMP_OFFLOAD
 if(gs_hamk%gpu_option==ABI_GPU_OPENMP) then
   !$OMP TARGET EXIT DATA MAP(from:cg,eig,resid,occ,pcon)
 end if
#endif

 call timab(tim_lobpcgwf2,2,tsec)

 DBG_EXIT("COLL")

end subroutine lobpcgwf2

subroutine getghc_gsc1(X,AX,BX)

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: AX
 type(xgBlock_t), intent(inout) :: BX

!Local variables-------------------------------
!scalars
 integer         :: blockdim
 integer         :: spacedim
 real(dp) :: eval,dum
 type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,1)
!arrays
 real(dp), pointer :: cg(:,:)
 real(dp), pointer :: ghc(:,:)
 real(dp), pointer :: gsc(:,:)
 real(dp), allocatable :: gvnlxc(:,:)

! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_GETGHC)

 call xgBlock_getSize(X,spacedim,blockdim)
 call xgBlock_check(X,AX)
 call xgBlock_check(X,BX)

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,rows=1,cols=spacedim*blockdim)

#ifdef FC_CRAY
 ABI_MALLOC(gvnlxc,(1,1))
#else
 ABI_MALLOC(gvnlxc,(0,0))
#endif

 if (l_mpi_enreg%nproc_fft==1.or.l_gs_hamk%istwf_k==1) then
   call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
     l_gs_hamk,gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)
 else if (l_gs_hamk%istwf_k==2) then ! nproc_fft>1 and istwfk==2
    call prep_getghc(cg(:,1:blockdim*spacedim),l_gs_hamk,gvnlxc,ghc,gsc(:,1:blockdim*spacedim),dum,blockdim,&
      l_mpi_enreg,l_prtvol,l_sij_opt,l_cpopt,cprj_dum,already_transposed=.true.)
 else ! nproc_fft>1 and istwfk>2
   ABI_ERROR('getghc in lobpcg not implemented for npfft>1 and istwfk>2')
 end if

 ABI_FREE(gvnlxc)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 call gpu_device_synchronize()
#endif

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

 ABI_NVTX_END_RANGE()

end subroutine getghc_gsc1

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

end module m_lobpcgwf
!!***
