!!****f* ABINIT/m_chebfiwf
!! NAME
!! m_chebfiwf
!!
!! FUNCTION
!! This module contains a routine updating the whole wave functions at a given k-point,
!! using the Chebyshev filtering method (2021 implementation using xG abstraction layer)
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2022 ABINIT group (BS)
!! This file is distributed under the terms of the
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

module m_chebfiwf

 use defs_abitypes
 use defs_basis
 use m_abicore
 use m_errors
 use m_fstrings
 use m_time

 use m_chebfi
 use m_chebfi2
 use m_invovl

 use m_cgtools,     only : dotprod_g
 use m_dtset,       only : dataset_type

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_getghc, prep_nonlop
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_getghc,      only : multithreaded_getghc

 use m_xg
 use m_xgTransposer

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
 use m_nvtx_data
#endif

 use, intrinsic :: iso_c_binding, only: c_associated,c_loc,c_ptr,c_f_pointer

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 implicit none

 private

 integer, parameter :: l_tim_getghc=7
 real(dp), parameter :: inv_sqrt2 = 1/sqrt2

! For use in getghc_gsc1
 integer, save :: l_cpopt
 integer, save :: l_icplx
 integer, save :: l_istwf
 integer, save :: l_npw
 integer, save :: l_nband_filter
 integer, save :: l_nspinor
 logical, save :: l_paw
 integer, save :: l_prtvol
 integer, save :: l_sij_opt
 integer, save :: l_paral_kgb
 integer, save :: l_useria
 integer, save :: l_block_sliced
 real(dp), allocatable,save ::  l_pcon(:)
 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 integer, parameter :: DEBUG_ROWS = 5
 integer, parameter :: DEBUG_COLUMNS = 5

 public :: chebfiwf2

 CONTAINS  !========================================================================================
!!***

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

subroutine chebfiwf2(cg,dtset,eig,enl_out,gs_hamk,kinpw,mpi_enreg,&
&                    nband,npw,nspinor,prtvol,resid)

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 type(mpi_type),target,intent(inout) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),intent(in) :: kinpw(npw)
 real(dp),target,intent(out) :: resid(nband)
 real(dp),intent(out) :: enl_out(nband)
 real(dp),target,intent(out) :: eig(nband)
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk

!Local variables-------------------------------
!scalars
 integer, parameter :: tim_chebfiwf2 = 1750
 integer :: ipw,space,blockdim,nline,total_spacedim,ierr,nthreads
 real(dp) :: cputime,walltime,localmem
 type(c_ptr) :: cptr
 type(chebfi_t) :: chebfi
 type(xgBlock_t) :: xgx0,xgeigen,xgresidu
!arrays
 real(dp) :: tsec(2),chebfiMem(2)
 real(dp),pointer :: eig_ptr(:,:) => NULL()
 real(dp),pointer :: resid_ptr(:,:) => NULL()
 real(dp), allocatable :: l_gvnlxc(:,:)

 !Stupid things for NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 real(dp) :: gsc_dummy(0,0)
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,1)

! *********************************************************************

!################ INITIALIZATION  #####################################
!######################################################################

  call timab(tim_chebfiwf2,1,tsec)
  cputime = abi_cpu_time()
  walltime = abi_wtime()

!Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_istwf=gs_hamk%istwf_k
 l_npw = npw
 l_nspinor = nspinor
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk
 l_nband_filter = nband
 l_paral_kgb = dtset%paral_kgb
 l_block_sliced = dtset%diago_apply_block_sliced

!Variables
 nline=dtset%nline
 blockdim=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp
 !for debug
 l_useria=dtset%useria

!Depends on istwfk
 if ( l_istwf == 2 ) then ! Real only
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

!Memory info
 if ( prtvol >= 3 ) then
   if (l_mpi_enreg%paral_kgb == 1) then
     total_spacedim = l_icplx*l_npw*l_nspinor
     call xmpi_sum(total_spacedim,l_mpi_enreg%comm_bandspinorfft,ierr)
   else
     total_spacedim = 0
   end if
   chebfiMem = chebfi_memInfo(nband,l_icplx*l_npw*l_nspinor,space,l_mpi_enreg%paral_kgb, &
&                             total_spacedim,l_mpi_enreg%bandpp) !blockdim
   localMem = (l_npw+2*l_npw*l_nspinor+2*nband)*kind(1.d0) !blockdim
   write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling chebfi should need around ", &
   (localMem+sum(chebfiMem))/1e9,"GB of peak memory as follows :"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in chebfiwf : ",(localMem)/1e9,"GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_chebfi : ",(chebfiMem(1))/1e9,"GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_chebfi : ",(chebfiMem(2))/1e9,"GB"
 end if

!For preconditionning
 ABI_MALLOC(l_pcon,(1:l_icplx*npw))

!$omp parallel do schedule(static), shared(l_pcon,kinpw)
 do ipw=1-1,l_icplx*npw-1
   if(kinpw(ipw/l_icplx+1)>huge(0.0_dp)*1.d-11) then
     l_pcon(ipw+1)=0.d0
   else
     l_pcon(ipw+1) = (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1)))) &
&    / (27+kinpw(ipw/l_icplx+1)*(18+kinpw(ipw/l_icplx+1)*(12+8*kinpw(ipw/l_icplx+1))) + 16*kinpw(ipw/l_icplx+1)**4)
   end if
 end do

 call xgBlock_map(xgx0,cg,space,l_icplx*l_npw*l_nspinor,nband,l_mpi_enreg%comm_bandspinorfft)

 if ( l_istwf == 2 ) then ! Real only
   ! Scale cg
   call xgBlock_scale(xgx0,sqrt2,1)  !ALL MPI processes do this

   ! This is possible since the memory in cg and xgx0 is the same
   ! Don't know yet how to deal with this with xgBlock
   !MPI HANDLES THIS AUTOMATICALLY (only proc 0 is me_g0)
   if(l_mpi_enreg%me_g0 == 1) cg(:, 1:npw*nspinor*nband:npw) = cg(:, 1:npw*nspinor*nband:npw) * inv_sqrt2
 end if

!Trick with C is to change rank of arrays (:) to (:,:)
 cptr = c_loc(eig)
 call c_f_pointer(cptr,eig_ptr,(/ nband,1 /))
 call xgBlock_map(xgeigen,eig_ptr,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)
!Trick the with C to change rank of arrays (:) to (:,:)
 cptr = c_loc(resid)
 call c_f_pointer(cptr,resid_ptr,(/ nband,1 /))
 call xgBlock_map(xgresidu,resid_ptr,SPACE_R,nband,1,l_mpi_enreg%comm_bandspinorfft)

! ABI_MALLOC(l_gvnlxc,(2,l_npw*l_nspinor*l_nband_filter))
 call timab(tim_chebfiwf2,2,tsec)

 cputime = abi_cpu_time() - cputime
 walltime = abi_wtime() - walltime

 nthreads = xomp_get_num_threads(open_parallel = .true.)

 if ( cputime/walltime/dble(nthreads) < 0.75 .and. (int(cputime/walltime)+1) /= nthreads) then
   if ( prtvol >= 3 ) then
     write(std_out,'(a)',advance='no') sjoin(" Chebfi took", sec2str(cputime), "of cpu time")
     write(std_out,*) sjoin("for a wall time of", sec2str(walltime))
     write(std_out,'(a,f6.2)') " -> Ratio of ", cputime/walltime
   end if
   ABI_COMMENT(sjoin("You should set the number of threads to something close to",itoa(int(cputime/walltime)+1)))
 end if

 call chebfi_init(chebfi,nband,l_icplx*l_npw*l_nspinor,dtset%tolwfr,dtset%ecut, &
&                 dtset%paral_kgb,l_mpi_enreg%nproc_band,l_mpi_enreg%bandpp, &
&                 l_mpi_enreg%nproc_fft,nline, space,1,l_gs_hamk%istwf_k, &
&                 l_mpi_enreg%comm_bandspinorfft,l_mpi_enreg%me_g0,l_paw)

!################    RUUUUUUUN    #####################################
!######################################################################

 call chebfi_run(chebfi,xgx0,getghc_gsc1,getBm1X,precond1,xgeigen,xgresidu,l_mpi_enreg)

!Free preconditionning since not needed anymore
 ABI_FREE(l_pcon)

!Compute enlout (nonlocal energy for each band if necessary) This is the best
!  quick and dirty trick to compute this part in NC. gvnlc cannot be part of
!  chebfi algorithm
 if ( .not. l_paw ) then
   !Check l_gvnlc size
   !if ( size(l_gvnlxc) < 2*nband*l_npw*l_nspinor ) then
   !if ( size(l_gvnlxc) /= 0 ) then
   !  ABI_FREE(l_gvnlxc)
   ABI_MALLOC(l_gvnlxc,(0,0))
   !end if

   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NONLOP)
   !Call nonlop
   call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
        &            signs,gsc_dummy,l_tim_getghc,cg,l_gvnlxc)
   ABI_NVTX_END_RANGE()
   ABI_FREE(l_gvnlxc)
 end if

!Free chebfi
 call chebfi_free(chebfi)

!################    SORRY IT'S ALREADY FINISHED : )  #################
!######################################################################

 call timab(tim_chebfiwf2,2,tsec)

 cputime = abi_cpu_time() - cputime
 walltime = abi_wtime() - walltime

 nthreads = xomp_get_num_threads(open_parallel = .true.)

 if ( cputime/walltime/dble(nthreads) < 0.75 .and. (int(cputime/walltime)+1) /= nthreads) then
   if ( prtvol >= 3 ) then
     write(std_out,'(a)',advance='no') sjoin(" Chebfi took", sec2str(cputime), "of cpu time")
     write(std_out,*) sjoin("for a wall time of", sec2str(walltime))
     write(std_out,'(a,f6.2)') " -> Ratio of ", cputime/walltime
   end if
   ABI_COMMENT(sjoin("You should set the number of threads to something close to",itoa(int(cputime/walltime)+1)))
 end if

 DBG_EXIT("COLL")

end subroutine chebfiwf2
!!***

!----------------------------------------------------------------------

!!****f* m_chebfiwf/getghc_gsc1
!! NAME
!! getghc_gsc1
!!
!! FUNCTION
!! This routine computes H|C> and possibly S|C> for a given wave function C.
!!  It acts as a driver for getghc, taken into account parallelism, multithreading, etc.
!!
!! SIDE EFFECTS
!!  X  <type(xgBlock_t)>= memory block containing |C>
!!  AX <type(xgBlock_t)>= memory block containing H|C>
!!  BX <type(xgBlock_t)>= memory block containing S|C>
!!  transposer <type(xgTransposer_t)>= data used for array transpositions
!!
!! SOURCE

subroutine getghc_gsc1(X,AX,BX,transposer)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: AX
 type(xgBlock_t), intent(inout) :: BX
 type(xgTransposer_t), optional, intent(inout) :: transposer
 integer         :: blockdim
 integer         :: spacedim
 type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,1)

!Local variables-------------------------------
!scalars
 integer :: cpuRow
 real(dp) :: eval
!arrays
 real(dp), pointer :: cg(:,:)
 real(dp), pointer :: ghc(:,:)
 real(dp), pointer :: gsc(:,:)
 real(dp), allocatable :: l_gvnlxc(:,:)

! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_GETGHC)

 call xgBlock_getSize(X,spacedim,blockdim)

 spacedim = spacedim/l_icplx

 call xgBlock_reverseMap(X,cg,l_icplx,spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,l_icplx,spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,l_icplx,spacedim*blockdim)

!Scale back cg
 if(l_istwf == 2) then
   call xgBlock_scale(X,inv_sqrt2,1)

   if (l_paral_kgb == 0) then
     if(l_mpi_enreg%me_g0 == 1) cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * sqrt2
   else
     cpuRow = xgTransposer_getRank(transposer, 2)
     if (cpuRow == 0) then
       cg(:, 1:spacedim*blockdim:spacedim) = cg(:, 1:spacedim*blockdim:spacedim) * sqrt2
     end if
   end if
 end if

 !if ( size(l_gvnlxc) < 2*blockdim*spacedim ) then
 !  ABI_FREE(l_gvnlxc)
 !  ABI_MALLOC(l_gvnlxc,(2,blockdim*spacedim))
 !end if
 ABI_MALLOC(l_gvnlxc,(0,0))

 call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,l_gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

 ABI_FREE(l_gvnlxc)

!Scale cg, ghc, gsc
 if ( l_istwf == 2 ) then
   call xgBlock_scale(X,sqrt2,1)
   call xgBlock_scale(AX,sqrt2,1)

   if (l_paral_kgb == 0) then
     if(l_mpi_enreg%me_g0 == 1) then
       cg(:, 1:spacedim*blockdim:l_npw) = cg(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
       ghc(:, 1:spacedim*blockdim:l_npw) = ghc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
     endif
   else
     if (cpuRow == 0) then
       cg(:, 1:spacedim*blockdim:spacedim) = cg(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
       ghc(:, 1:spacedim*blockdim:spacedim) = ghc(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
     end if
   end if
   if(l_paw) then
     call xgBlock_scale(BX,sqrt2,1)
     if (l_paral_kgb == 0) then
       if(l_mpi_enreg%me_g0 == 1) gsc(:, 1:spacedim*blockdim:l_npw) = gsc(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
     else
       if (cpuRow == 0) then
         gsc(:, 1:spacedim*blockdim:spacedim) = gsc(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
       end if
     end if
   end if
 end if

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

 ABI_NVTX_END_RANGE()

end subroutine getghc_gsc1
!!***

!----------------------------------------------------------------------

!!****f* m_chebfiwf/getBm1X
!! NAME
!! getBm1X
!!
!! FUNCTION
!! This routine computes S^-1|C> for a given wave function C.
!!  It acts as a driver for apply_invovl.
!!
!! SIDE EFFECTS
!!  X  <type(xgBlock_t)>= memory block containing |C>
!!  Bm1X <type(xgBlock_t)>= memory block containing S^-1|C>
!!  transposer <type(xgTransposer_t)>= data used for array transpositions
!!
!! SOURCE

subroutine getBm1X(X,Bm1X,transposer)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: Bm1X
 type(xgTransposer_t), optional, intent(inout) :: transposer

!Local variables-------------------------------
!scalars
 integer :: blockdim
 integer :: spacedim
 integer :: cpuRow
!arrays
 real(dp), pointer :: ghc_filter(:,:)
 real(dp), pointer :: gsm1hc_filter(:,:)
 type(pawcprj_type), allocatable :: cwaveprj_next(:,:) !dummy

! *********************************************************************

 call xgBlock_getSize(X,spacedim,blockdim)

 spacedim = spacedim/l_icplx

 call xgBlock_reverseMap(X,ghc_filter,l_icplx,spacedim*blockdim)

 call xgBlock_reverseMap(Bm1X,gsm1hc_filter,l_icplx,spacedim*blockdim)

 !scale back cg
 if(l_istwf == 2) then
   call xgBlock_scale(X,inv_sqrt2,1)
   if (l_paral_kgb == 0) then
     if(l_mpi_enreg%me_g0 == 1) ghc_filter(:, 1:spacedim*blockdim:l_npw) = ghc_filter(:, 1:spacedim*blockdim:l_npw) * sqrt2
   else
     cpuRow = xgTransposer_getRank(transposer, 2)
     if (cpuRow == 0) then
       ghc_filter(:, 1:spacedim*blockdim:spacedim) = ghc_filter(:, 1:spacedim*blockdim:spacedim) * sqrt2
     end if
   end if
   if(l_paw) then
     call xgBlock_scale(Bm1X,inv_sqrt2,1)
     if (l_paral_kgb == 0) then
       if(l_mpi_enreg%me_g0 == 1) &
&        gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) = gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) * sqrt2
     else
       if (cpuRow == 0) then
         gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) = gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) * sqrt2
       end if
     end if
   end if
 end if

 if(l_paw) then
   ABI_MALLOC(cwaveprj_next, (l_gs_hamk%natom,l_nspinor*blockdim))
   call pawcprj_alloc(cwaveprj_next,0,l_gs_hamk%dimcprj)
   call apply_invovl(l_gs_hamk, ghc_filter(:,:), gsm1hc_filter(:,:), cwaveprj_next(:,:), &
&       spacedim, blockdim, l_mpi_enreg, l_nspinor, l_block_sliced)
 else
   gsm1hc_filter(:,:) = ghc_filter(:,:)
 end if

!Scale cg, ghc, gsc
 if ( l_istwf == 2 ) then
   call xgBlock_scale(X,sqrt2,1)
   if (l_paral_kgb == 0) then
     if(l_mpi_enreg%me_g0 == 1) then
       ghc_filter(:, 1:spacedim*blockdim:l_npw) = ghc_filter(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
     endif
   else
     if (cpuRow == 0) then
       ghc_filter(:, 1:spacedim*blockdim:spacedim) = ghc_filter(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
     end if
   end if
   if(l_paw) then
     call xgBlock_scale(Bm1X,sqrt2,1)
     if (l_paral_kgb == 0) then
       if(l_mpi_enreg%me_g0 == 1) &
&        gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) = gsm1hc_filter(:, 1:spacedim*blockdim:l_npw) * inv_sqrt2
     else
       if (cpuRow == 0) then
         gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) = gsm1hc_filter(:, 1:spacedim*blockdim:spacedim) * inv_sqrt2
       end if
     end if
   end if
 end if

 if (l_paw) then
   if (l_useria /= 121212) then
     call pawcprj_free(cwaveprj_next)
     ABI_FREE(cwaveprj_next)
   end if
 end if

end subroutine getBm1X
!!***

!----------------------------------------------------------------------

!!****f* m_chebfiwf/precond1
!! NAME
!! precond1
!!
!! FUNCTION
!! This routine applies a preconditionning to a block of memory
!!
!! SIDE EFFECTS
!!  W <type(xgBlock_t)>= memory block
!!
!! SOURCE

subroutine precond1(W)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: W

!Local variables-------------------------------
!scalars
 integer :: ispinor

! *********************************************************************

!Precondition resid_vec
 do ispinor = 1,l_nspinor
   call xgBlock_colwiseMul(W,l_pcon,l_icplx*l_npw*(ispinor-1))
 end do

end subroutine precond1
!!***

!----------------------------------------------------------------------

end module m_chebfiwf
!!***
