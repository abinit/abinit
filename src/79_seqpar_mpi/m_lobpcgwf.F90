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
!! Copyright (C) 1998-2022 ABINIT group (JB)
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
 use m_prep_kgb,    only : prep_getghc, prep_nonlop
 use m_getghc,      only : multithreaded_getghc
 use m_cgtools,     only : dotprod_g

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

 use, intrinsic :: iso_c_binding

 private

 integer, parameter :: l_tim_getghc=5
 double precision, parameter :: inv_sqrt2 = 1/sqrt2

 ! For use in getghc_gsc1
 integer,save  :: l_cpopt
 integer,save  :: l_icplx
 integer,save  :: l_istwf
 integer,save  :: l_npw
 integer,save  :: l_nspinor
 logical,save  :: l_paw
 integer, save :: l_paral_kgb
 integer,save  :: l_prtvol
 integer,save  :: l_sij_opt
 real(dp), allocatable,save ::  l_pcon(:)
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 public :: lobpcgwf2

 contains

subroutine lobpcgwf2(cg,dtset,eig,occ,enl_out,gs_hamk,isppol,ikpt,inonsc,istep,kinpw,mpi_enreg,&
&                   nband,npw,nspinor,prtvol,resid,nbdbuf)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 integer,intent(in) :: isppol,ikpt,inonsc,istep,nbdbuf
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk
 type(dataset_type)              ,intent(in   ) :: dtset
 type(mpi_type)           ,target,intent(in)    :: mpi_enreg
 real(dp)                 ,target,intent(inout) :: cg(2,nspinor*nband*npw)!,gsc(2,nspinor*nband*npw)
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
 type(lobpcg_t) :: lobpcg

 integer :: space, blockdim,  nline
 integer :: ipw

 integer, parameter :: tim_lobpcgwf2 = 1640
 double precision :: cputime, walltime
 double precision :: tsec(2)

 type(c_ptr) :: cptr
 real(dp), pointer :: eig_ptr(:,:) => NULL()
 real(dp), pointer :: resid_ptr(:,:) => NULL()
 real(dp), pointer :: occ_ptr(:,:) => NULL()

 ! Important things for NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 type(pawcprj_type) :: cprj_dum(1,1)
 integer :: iband, shift, me_g0, me_g0_fft
 real(dp) :: gsc_dummy(0,0)
 real(dp), allocatable :: l_gvnlxc(:,:)

! *********************************************************************


!###########################################################################
!################ INITIALISATION  ##########################################
!###########################################################################

 call timab(tim_lobpcgwf2,1,tsec)
 cputime = abi_cpu_time()
 walltime = abi_wtime()

 ! Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_istwf=gs_hamk%istwf_k
 l_npw = npw
 l_nspinor = nspinor
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk
 l_paral_kgb = dtset%paral_kgb

!Variables
 nline=dtset%nline
 blockdim=nband/dtset%nblock_lobpcg !=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp

!Depends on istwfk
 if ( l_istwf > 1 ) then ! Real only
   ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
   ! and im*im so that a vector of complex is consider as a long vector of real
   ! therefore the number of data is (2*npw*nspinor)*nband
   ! This space is completely equivalent to SPACE_R but will correctly set and
   ! get the array data into the xgBlock
   space = SPACE_CR
   l_icplx = 2
 else ! complex
   space = SPACE_C
   l_icplx = 1
 end if

 !For preconditionning
 ABI_MALLOC(l_pcon,(1:l_icplx*npw))
 !$omp parallel do schedule(static), shared(l_pcon,kinpw)
 do ipw=1-1,l_icplx*npw-1
   if(kinpw(ipw/l_icplx+1)>huge(zero)*1.d-11) then
     l_pcon(ipw+1)=0.d0
   else
     l_pcon(ipw+1) = (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1)))) &
&    / (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1))) + 16*kinpw(ipw/l_icplx+1)**4)
   end if
 end do

#ifdef HAVE_OPENMP_OFFLOAD
 if(gs_hamk%gpu_option==ABI_GPU_OPENMP) then
   !$OMP TARGET ENTER DATA MAP(to:cg,eig,resid)
 end if
#endif

 ! Local variables for lobpcg
 me_g0 = -1
 me_g0_fft = -1
 if (space==SPACE_CR) then
   me_g0 = 0
   me_g0_fft = 0
   if (l_istwf == 2) then
     if (l_mpi_enreg%me_g0 == 1) me_g0 = 1
     if (l_mpi_enreg%me_g0_fft == 1) me_g0_fft = 1
   end if
 end if
 call xgBlock_map(xgx0,cg,space,l_npw*l_nspinor,nband,comm=l_mpi_enreg%comm_bandspinorfft,me_g0=me_g0,&
   & gpu_option=dtset%gpu_option)

 ! Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(eig)
 call c_f_pointer(cptr,eig_ptr,(/ nband,1 /))
 call xgBlock_map(xgeigen,eig_ptr,SPACE_R,nband,1,gpu_option=dtset%gpu_option)

 !call xg_init(xgresidu,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 ! Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(resid)
 call c_f_pointer(cptr,resid_ptr,(/ nband,1 /))
 call xgBlock_map(xgresidu,resid_ptr,SPACE_R,nband,1,gpu_option=dtset%gpu_option)

 !call xg_init(xgeigen,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
 ! Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(occ)
 call c_f_pointer(cptr,occ_ptr,(/ nband,1 /))
 call xgBlock_map(xgocc,occ_ptr,SPACE_R,nband,1,gpu_option=dtset%gpu_option)

 call lobpcg_init(lobpcg,nband,l_npw*l_nspinor,blockdim,dtset%tolwfr,nline,&
   space,l_mpi_enreg%comm_bandspinorfft,l_paral_kgb,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
   me_g0,me_g0_fft,gs_hamk%gpu_option)

!###########################################################################
!################    RUUUUUUUN    ##########################################
!###########################################################################

 ! Run lobpcg
 call lobpcg_run(lobpcg,xgx0,getghc_gsc1,precond,xgeigen,xgocc,xgresidu,prtvol,nspinor,isppol,ikpt,inonsc,istep,nbdbuf)

 ! Free preconditionning since not needed anymore
 ABI_FREE(l_pcon)

 if ( .not. l_paw ) then
#ifdef FC_CRAY
   ABI_MALLOC(l_gvnlxc,(1,1))
#else
   ABI_MALLOC(l_gvnlxc,(0,0))
#endif

   !Call nonlop
   if (l_paral_kgb==0) then

     call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,l_gvnlxc)

   else
     do iband=1,nband/blockdim
       shift = (iband-1)*blockdim*l_npw*l_nspinor
       call prep_nonlop(choice,l_cpopt,cprj_dum, &
&        enl_out((iband-1)*blockdim+1:iband*blockdim),l_gs_hamk,0,&
&        eig((iband-1)*blockdim+1:iband*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&        gsc_dummy,l_tim_getghc, &
&        cg(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
!&        l_gvnlxc(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
&        l_gvnlxc(:,:),&
&        already_transposed=.false.)
     end do
   end if
   ABI_FREE(l_gvnlxc)
 end if

 ! Free lobpcg
 call lobpcg_free(lobpcg)

#ifdef HAVE_OPENMP_OFFLOAD
 if(gs_hamk%gpu_option==ABI_GPU_OPENMP) then
   !$OMP TARGET EXIT DATA MAP(from:cg,eig,resid)
 end if
#endif
!###########################################################################
!################    SORRY IT'S ALREADY FINISHED : )  ######################
!###########################################################################


 call timab(tim_lobpcgwf2,2,tsec)

 DBG_EXIT("COLL")

end subroutine lobpcgwf2

subroutine getghc_gsc1(X,AX,BX)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: AX
 type(xgBlock_t), intent(inout) :: BX
 integer         :: blockdim
 integer         :: spacedim
 type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,1)

!Local variables-------------------------------
!scalars
 real(dp) :: eval
!arrays
 real(dp), pointer :: cg(:,:)
 real(dp), pointer :: ghc(:,:)
 real(dp), pointer :: gsc(:,:)
 real(dp)          :: l_gvnlxc(1,1)

! *********************************************************************

 call xgBlock_getSize(X,spacedim,blockdim)
 call xgBlock_check(X,AX)
 call xgBlock_check(X,BX)

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,rows=1,cols=spacedim*blockdim)

 call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,l_gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 call gpu_device_synchronize()
#endif

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

end subroutine getghc_gsc1

 subroutine precond(W)
   use m_xg, only : xg_t, xgBlock_colwiseMul
   type(xgBlock_t), intent(inout) :: W
   integer :: ispinor

   ! precondition resid_vec
   do ispinor = 1,l_nspinor
     call xgBlock_colwiseMul(W,l_pcon,l_npw*(ispinor-1))
   end do

 end subroutine precond

 end module m_lobpcgwf
!!***
