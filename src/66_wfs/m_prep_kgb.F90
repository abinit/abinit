!!****m* ABINIT/m_prep_kgb
!! NAME
!!  m_prep_kgb
!!
!! FUNCTION
!!  This module provides wrappers that used to apply the full Hamiltonian or just the Vnl part
!!  or to perform the FFT of the wavefunctions when the orbitals are distributed in linalg mode (paral_kgb = 1).
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (FBottin,MT,GZ,MD,FDahm)
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

module m_prep_kgb

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_bandfft_kpt, only : bandfft_kpt, bandfft_kpt_get_ikpt, bandfft_kpt_type
 use m_pawcprj,     only : pawcprj_type
 use m_hamiltonian, only : gs_hamiltonian_type
 use m_nonlop,      only : nonlop
 use m_getghc,      only : multithreaded_getghc
 use m_fft,         only : fourwf

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

 implicit none

 private
!!***

 public :: prep_getghc
 public :: prep_nonlop
 public :: prep_fourwf
 public :: prep_wavef_sym_do
 public :: prep_wavef_sym_undo
 public :: prep_index_wavef_bandpp
 public :: prep_sort_wavef_spin
!!***

contains
!!***

!!****f* ABINIT/prep_getghc
!! NAME
!! prep_getghc
!!
!! FUNCTION
!! this routine prepares the data to the call of getghc.
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cprjin%cp(:)=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*my_nspinor*blocksize)=planewave coefficients of wavefunction.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  gvnlxc=matrix elements <G|Vnonlocal+VFockACE|C>
!!  lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!         Typically lambda is the eigenvalue (or its guess)
!!  mpi_enreg=information about mpi parallelization
!!  prtvol=control print volume and debugging output
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                       if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!!
!! OUTPUT
!!  gwavef=(2,npw*my_nspinor*blocksize)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                  or <G|H-lambda.S|C> (if sij_opt=-1).
!!  swavef=(2,npw*my_nspinor*blocksize)=matrix elements <G|S|C>.
!!
!! SIDE EFFECTS
!!  ====== if gs_hamk%usepaw==1
!!  cwaveprj(natom,my_nspinor*bandpp)= wave functions at k projected with nl projectors
!!
!! PARENTS
!!      m_chebfi,m_dft_energy,m_lobpcgwf,m_lobpcgwf_old
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_getghc(cwavef, gs_hamk, gvnlxc, gwavef, swavef, lambda, blocksize, &
                       mpi_enreg, prtvol, sij_opt, cpopt, cwaveprj, &
                       already_transposed) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,cpopt,prtvol,sij_opt
 logical, intent(in),optional :: already_transposed
 real(dp),intent(in) :: lambda
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cwavef(:,:)
 real(dp),intent(inout) :: gvnlxc (:,:),gwavef(:,:),swavef(:,:)
 type(pawcprj_type), intent(inout) :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=6
 integer :: bandpp,bandpp_sym,idatarecv0,ier,ikpt_this_proc,iscalc,mcg,my_nspinor
 integer :: nbval,ndatarecv,ndatarecv_tot,ndatasend_sym,nproc_band,nproc_fft
 integer :: old_me_g0,spaceComm
 logical :: flag_inv_sym, do_transpose
 !character(len=500) :: msg
!arrays
 integer,allocatable :: index_wavef_band(:),index_wavef_send(:),index_wavef_spband(:)
 integer,allocatable :: rdisplsloc(:),recvcountsloc(:),sdisplsloc(:),sendcountsloc(:)
 integer,ABI_CONTIGUOUS pointer :: kg_k_gather_sym(:,:)
 integer,ABI_CONTIGUOUS pointer :: rdispls(:),rdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: recvcounts(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,ABI_CONTIGUOUS pointer :: sdispls(:),sdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: sendcounts(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,ABI_CONTIGUOUS pointer :: tab_proc(:)
 real(dp) :: tsec(2)
 real(dp),allocatable,target :: cwavef_alltoall1(:,:),gvnlxc_alltoall1(:,:)
 real(dp),allocatable,target :: gwavef_alltoall1(:,:),swavef_alltoall1(:,:)
 real(dp),allocatable,target :: cwavef_alltoall2(:,:),gvnlxc_alltoall2(:,:)
 real(dp),allocatable,target :: gwavef_alltoall2(:,:),swavef_alltoall2(:,:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),gvnlxc_alltoall_sym(:,:)
 real(dp),pointer :: gwavef_alltoall_sym(:,:)
 real(dp),pointer :: swavef_alltoall_sym(:,:)

! *************************************************************************

 call timab(630,1,tsec)
 call timab(631,3,tsec)

!Some inits
 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp
 my_nspinor = max(1,gs_hamk%nspinor/mpi_enreg%nproc_spinor)

 do_transpose = .true.
 if(present(already_transposed)) then
   if(already_transposed) do_transpose = .false.
 end if

 flag_inv_sym = (gs_hamk%istwf_k==2 .and. any(gs_hamk%ngfft(7) == [401,402,312,512]))
 if (flag_inv_sym) then
   gs_hamk%istwf_k = 1
   if (modulo(bandpp,2)==0) bandpp_sym = bandpp/2
   if (modulo(bandpp,2)/=0) bandpp_sym = bandpp
 end if

!Check sizes
 mcg=2*gs_hamk%npw_fft_k*my_nspinor*bandpp
 if (do_transpose) mcg=2*gs_hamk%npw_k*my_nspinor*blocksize
 if (size(cwavef)<mcg) then
   ABI_BUG('wrong size for cwavef!')
 end if
 if (size(gwavef)<mcg) then
   ABI_BUG('wrong size for gwavef!')
 end if
 if (size(gvnlxc)<mcg) then
   ABI_BUG('wrong size for gvnlxc!')
 end if
 if (sij_opt==1) then
   if (size(swavef)<mcg) then
     ABI_BUG('wrong size for swavef!')
   end if
 end if
 if (gs_hamk%usepaw==1.and.cpopt>=0) then
   if (size(cwaveprj)<gs_hamk%natom*my_nspinor*bandpp) then
     ABI_BUG('wrong size for cwaveprj!')
   end if
 end if

!====================================================================================

 spaceComm=mpi_enreg%comm_fft
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band

 ikpt_this_proc=bandfft_kpt_get_ikpt()

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   =>bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   =>bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      =>bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      =>bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    = bandfft_kpt(ikpt_this_proc)%ndatarecv

 if (flag_inv_sym ) then
   idatarecv0           = bandfft_kpt(ikpt_this_proc)%idatarecv0
   ndatarecv_tot        = bandfft_kpt(ikpt_this_proc)%ndatarecv_tot
   ndatasend_sym        = bandfft_kpt(ikpt_this_proc)%ndatasend_sym
   kg_k_gather_sym      =>bandfft_kpt(ikpt_this_proc)%kg_k_gather_sym(:,:)
   rdispls_sym          =>bandfft_kpt(ikpt_this_proc)%rdispls_sym(:)
   recvcounts_sym       =>bandfft_kpt(ikpt_this_proc)%recvcounts_sym(:)
   recvcounts_sym_tot   =>bandfft_kpt(ikpt_this_proc)%recvcounts_sym_tot(:)
   sdispls_sym          =>bandfft_kpt(ikpt_this_proc)%sdispls_sym(:)
   sendcounts_sym       =>bandfft_kpt(ikpt_this_proc)%sendcounts_sym(:)
   sendcounts_sym_all   =>bandfft_kpt(ikpt_this_proc)%sendcounts_sym_all(:)
   tab_proc             =>bandfft_kpt(ikpt_this_proc)%tab_proc(:)
 end if
 iscalc=(sij_opt+1)/2  ! 0 if S not calculated, 1 otherwise
 nbval=(ndatarecv*my_nspinor*bandpp)*iscalc

 if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. my_nspinor==2 ).or. &
& ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   ABI_ALLOCATE(gwavef_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   ABI_ALLOCATE(swavef_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   ABI_ALLOCATE(gvnlxc_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   swavef_alltoall1(:,:)=zero
   gvnlxc_alltoall1(:,:)=zero
   cwavef_alltoall1(:,:)=zero
   gwavef_alltoall1(:,:)=zero
 end if
 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 ABI_ALLOCATE(gwavef_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 ABI_ALLOCATE(swavef_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 ABI_ALLOCATE(gvnlxc_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 swavef_alltoall2(:,:)=zero
 gvnlxc_alltoall2(:,:)=zero
 cwavef_alltoall2(:,:)=zero
 gwavef_alltoall2(:,:)=zero

 recvcountsloc(:)=recvcounts(:)*2*my_nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*my_nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*my_nspinor
 sdisplsloc(:)=sdispls(:)*2*my_nspinor
 call timab(631,2,tsec)

 if(do_transpose) then
   call timab(545,3,tsec)
   if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. my_nspinor==2 ).or. &
&   ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall1,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   else
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   end if
   call timab(545,2,tsec)
 else
   ! Here, we cheat, and use DCOPY to bypass some compiler's overzealous bound-checking
   ! (ndatarecv*my_nspinor*bandpp might be greater than the declared size of cwavef)
   call DCOPY(2*ndatarecv*my_nspinor*bandpp, cwavef, 1, cwavef_alltoall2, 1)
 end if

 if(gs_hamk%istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if

!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then
   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.my_nspinor==2)then
     call timab(632,3,tsec)
!    Sort to have all ispinor=1 first, then all ispinor=2
     call prep_sort_wavef_spin(nproc_band,my_nspinor,ndatarecv,recvcounts,rdispls,index_wavef_spband)
     cwavef_alltoall2(:,:)=cwavef_alltoall1(:,index_wavef_spband)
     call timab(632,2,tsec)
   end if

   call timab(635,3,tsec)
   call multithreaded_getghc(cpopt,cwavef_alltoall2,cwaveprj,gwavef_alltoall2,swavef_alltoall2(:,1:nbval),&
&   gs_hamk,gvnlxc_alltoall2,lambda,mpi_enreg,1,prtvol,sij_opt,tim_getghc,0)
   call timab(635,2,tsec)

   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.my_nspinor==2)then
     call timab(634,3,tsec)
     gwavef_alltoall1(:,index_wavef_spband)=gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_spband)=swavef_alltoall2(:,:)
     gvnlxc_alltoall1(:,index_wavef_spband)=gvnlxc_alltoall2(:,:)
     ABI_DEALLOCATE(index_wavef_spband)
     call timab(634,2,tsec)
   end if

 else if ((.not.(flag_inv_sym)) .and. (bandpp>1)) then
!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------

   if(do_transpose) then
     call timab(632,3,tsec)
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     my_nspinor,ndatarecv, recvcounts,rdispls, index_wavef_band)
!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
     call timab(632,2,tsec)
   end if

!  ----------------------
!  Fourier transformation
!  ----------------------
   call timab(636,3,tsec)
   call multithreaded_getghc(cpopt,cwavef_alltoall2,cwaveprj,gwavef_alltoall2,swavef_alltoall2,gs_hamk,&
&   gvnlxc_alltoall2,lambda,mpi_enreg,bandpp,prtvol,sij_opt,tim_getghc,0)
   call timab(636,2,tsec)

!  -----------------------------------------------------
!  Sorting of waves functions below the processors
!  -----------------------------------------------------
   if(do_transpose) then
     call timab(634,3,tsec)
     gwavef_alltoall1(:,index_wavef_band) = gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_band) = swavef_alltoall2(:,:)
     gvnlxc_alltoall1(:,index_wavef_band)  = gvnlxc_alltoall2(:,:)
     ABI_DEALLOCATE(index_wavef_band)
     call timab(634,2,tsec)
   end if


 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   if(do_transpose) then
     call timab(632,3,tsec)
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     my_nspinor,ndatarecv,&
&     recvcounts,rdispls,&
&     index_wavef_band)

!  -------------------------------------------------------
!  Sorting the wave functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
   end if

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,my_nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,tab_proc,&
&   cwavef_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

!  ------------------------------------------------------------
!  Allocation
!  ------------------------------------------------------------
   ABI_ALLOCATE(gwavef_alltoall_sym,(2,ndatarecv_tot*bandpp_sym))
   ABI_ALLOCATE(swavef_alltoall_sym,(2,(ndatarecv_tot*bandpp_sym)*iscalc))
   ABI_ALLOCATE(gvnlxc_alltoall_sym ,(2,ndatarecv_tot*bandpp_sym))

   gwavef_alltoall_sym(:,:)=zero
   swavef_alltoall_sym(:,:)=zero
   gvnlxc_alltoall_sym(:,:)=zero

   call timab(632,2,tsec)

!  ------------------------------------------------------------
!  Fourier calculation
!  ------------------------------------------------------------
   call timab(637,3,tsec)
   call multithreaded_getghc(cpopt,ewavef_alltoall_sym,cwaveprj,gwavef_alltoall_sym,swavef_alltoall_sym,gs_hamk,&
&   gvnlxc_alltoall_sym,lambda,mpi_enreg,bandpp_sym,prtvol,sij_opt,tim_getghc,1,&
&   kg_fft_k=kg_k_gather_sym)
   call timab(637,2,tsec)

   call timab(633,3,tsec)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,my_nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   gwavef_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   gwavef_alltoall_sym,&
&   index_wavef_send)
   if (sij_opt==1)then
     call prep_wavef_sym_undo(mpi_enreg,bandpp,my_nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,idatarecv0,&
&     swavef_alltoall2,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     swavef_alltoall_sym,&
&     index_wavef_send)
   end if
   call prep_wavef_sym_undo(mpi_enreg,bandpp,my_nspinor,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   gvnlxc_alltoall2,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   gvnlxc_alltoall_sym,&
&   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)
   ABI_DEALLOCATE(gwavef_alltoall_sym)
   ABI_DEALLOCATE(swavef_alltoall_sym)
   ABI_DEALLOCATE(gvnlxc_alltoall_sym)

!  -------------------------------------------
!  We call getghc to calculate the nl matrix elements.
!  --------------------------------------------
   gs_hamk%istwf_k=2
   !!write(std_out,*)"Setting iswfk_k to 2"

   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
   call timab(633,2,tsec)

   call timab(638,3,tsec)
   call multithreaded_getghc(cpopt,cwavef_alltoall2,cwaveprj,gwavef_alltoall2,swavef_alltoall2,gs_hamk,&
&   gvnlxc_alltoall2,lambda,mpi_enreg,bandpp,prtvol,sij_opt,tim_getghc,2)
   call timab(638,2,tsec)

   call timab(634,3,tsec)
   mpi_enreg%me_g0=old_me_g0

   gs_hamk%istwf_k=1

!  -------------------------------------------------------
!  Sorting the wave functions below the processors
!  -------------------------------------------------------
   if(do_transpose) then
!    cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)   ! NOT NEEDED
     gwavef_alltoall1(:,index_wavef_band) = gwavef_alltoall2(:,:)
     if (sij_opt==1) swavef_alltoall1(:,index_wavef_band) = swavef_alltoall2(:,:)
     gvnlxc_alltoall1(:,index_wavef_band)  = gvnlxc_alltoall2(:,:)
     ABI_DEALLOCATE(index_wavef_band)
     call timab(634,2,tsec)
   end if

 end if
!====================================================================

 if (gs_hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0

 if(do_transpose) then

   call timab(545,3,tsec)
   if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. my_nspinor==2 ).or. &
&   ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
     if (sij_opt==1) then
       call xmpi_alltoallv(swavef_alltoall1,recvcountsloc,rdisplsloc,swavef,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     call xmpi_alltoallv(gvnlxc_alltoall1,recvcountsloc,rdisplsloc,gvnlxc,&
&     sendcountsloc,sdisplsloc,spaceComm,ier)
     call xmpi_alltoallv(gwavef_alltoall1,recvcountsloc,rdisplsloc,gwavef,&
&     sendcountsloc,sdisplsloc,spaceComm,ier)
   else
     if (sij_opt==1) then
       call xmpi_alltoallv(swavef_alltoall2,recvcountsloc,rdisplsloc,swavef,&
&       sendcountsloc,sdisplsloc,spaceComm,ier)
     end if
     call xmpi_alltoallv(gvnlxc_alltoall2,recvcountsloc,rdisplsloc,gvnlxc,&
&     sendcountsloc,sdisplsloc,spaceComm,ier)
     call xmpi_alltoallv(gwavef_alltoall2,recvcountsloc,rdisplsloc,gwavef,&
&     sendcountsloc,sdisplsloc,spaceComm,ier)
   end if

   call timab(545,2,tsec)
 else
   if(sij_opt == 1) then
     call DCOPY(2*ndatarecv*my_nspinor*bandpp, swavef_alltoall2, 1, swavef, 1)
   end if
   call DCOPY(2*ndatarecv*my_nspinor*bandpp, gvnlxc_alltoall2, 1, gvnlxc, 1)
   call DCOPY(2*ndatarecv*my_nspinor*bandpp, gwavef_alltoall2, 1, gwavef, 1)
 end if

!====================================================================
 if (flag_inv_sym) then
   gs_hamk%istwf_k = 2
 end if
!====================================================================
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 ABI_DEALLOCATE(gwavef_alltoall2)
 ABI_DEALLOCATE(gvnlxc_alltoall2)
 ABI_DEALLOCATE(swavef_alltoall2)

 if ( ((.not.flag_inv_sym) .and. bandpp==1 .and. mpi_enreg%paral_spinor==0 .and. my_nspinor==2 ).or. &
& ((.not.flag_inv_sym) .and. bandpp>1) .or.  flag_inv_sym  ) then
   ABI_DEALLOCATE(cwavef_alltoall1)
   ABI_DEALLOCATE(gwavef_alltoall1)
   ABI_DEALLOCATE(gvnlxc_alltoall1)
   ABI_DEALLOCATE(swavef_alltoall1)
 end if

 call timab(630,2,tsec)

end subroutine prep_getghc
!!***

!!****f* abinit/prep_nonlop
!! NAME
!! prep_nonlop
!!
!! FUNCTION
!! this routine prepares the data to the call of nonlop.
!!
!! INPUTS
!!  choice: chooses possible output:
!!    choice=1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a 2nd derivative with respect to atomic pos.
!!          =24=> a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!          =7 => no operator, just projections
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cwaveprj=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*my_nspinor*blocksize)=planewave coefficients of wavefunction.
!!  gvnlxc=matrix elements <G|Vnonlocal+VFockACE|C>
!!  hamk <type(gs_hamiltonian_type)>=data defining the Hamiltonian at a given k (NL part involved here)
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  lambdablock(blocksize)=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  mpi_enreg=information about mpi parallelization
!!  nnlout=dimension of enlout (when signs=1):
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!
!! OUTPUT
!!  ==== if (signs==1) ====
!!  enlout_block(nnlout)=
!!    if paw_opt==0, 1 or 2: contribution of this block of states to the nl part of various properties
!!    if paw_opt==3:         contribution of this block of states to <c|S|c>  (where S=overlap when PAW)
!!  ==== if (signs==2) ====
!!    if paw_opt==0, 1, 2 or 4:
!!       gvnlc(2,my_nspinor*npw)=result of the application of the nl operator
!!                        or one of its derivative to the input vect.
!!    if paw_opt==3 or 4:
!!       gsc(2,my_nspinor*npw*(paw_opt/3))=result of the aplication of (I+S)
!!                        to the input vect. (where S=overlap when PAW)
!!
!! SIDE EFFECTS
!!  ==== ONLY IF useylm=1
!!  cwaveprj(natom,my_nspinor) <type(pawcprj_type)>=projected input wave function |c> on non-local projector
!!                                  =<p_lmn|c> and derivatives
!!                                  Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  have to be computed (and not saved)
!!                     if cpopt= 0, <p_lmn|in> have to be computed and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives have to be computed and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  only derivatives are computed here and not saved
!! (if useylm=0, should have cpopt=-1)
!!
!! PARENTS
!!      m_dft_energy,m_forstr,m_invovl,m_lobpcgwf,m_vtowfk
!!
!! CHILDREN
!!
!! NOTES
!!  cprj (as well as cg) is distributed over band processors.
!!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projected WFs are stored on each proc.
!!
!! SOURCE

subroutine prep_nonlop(choice,cpopt,cwaveprj,enlout_block,hamk,idir,lambdablock,&
                       blocksize,mpi_enreg,nnlout,paw_opt,signs,gsc, tim_nonlop,cwavef,gvnlc, &
                       already_transposed) ! optional

!Arguments ------------------------------------
 integer,intent(in) :: blocksize,choice,cpopt,idir,signs,nnlout,paw_opt
 logical,optional,intent(in) :: already_transposed
 real(dp),intent(in) :: lambdablock(blocksize)
 real(dp),intent(out) :: enlout_block(nnlout*blocksize),gvnlc(:,:),gsc(:,:)
 real(dp),intent(inout) :: cwavef(:,:)
 type(gs_hamiltonian_type),intent(in) :: hamk
 type(mpi_type),intent(inout) :: mpi_enreg
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp,ier,ikpt_this_proc,my_nspinor,ndatarecv,nproc_band,npw,nspinortot
 integer :: old_me_g0,spaceComm=0,tim_nonlop
 logical :: do_transpose
 !character(len=500) :: msg
!arrays
 integer,allocatable :: index_wavef_band(:)
 integer,  allocatable :: rdisplsloc(:),recvcountsloc(:),sdisplsloc(:),sendcountsloc(:)
 integer,ABI_CONTIGUOUS  pointer :: rdispls(:),recvcounts(:),sdispls(:),sendcounts(:)
 real(dp) :: lambda_nonlop(mpi_enreg%bandpp),tsec(2)
 real(dp), allocatable :: cwavef_alltoall2(:,:),gvnlc_alltoall2(:,:),gsc_alltoall2(:,:)
 real(dp), allocatable :: cwavef_alltoall1(:,:),gvnlc_alltoall1(:,:),gsc_alltoall1(:,:)
 real(dp),allocatable :: enlout(:)

! *************************************************************************

 DBG_ENTER('COLL')

 call timab(570,1,tsec)

 do_transpose = .true.
 if(present(already_transposed)) then
   if(already_transposed) do_transpose = .false.
 end if

 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp
 spaceComm=mpi_enreg%comm_fft
 if(mpi_enreg%paral_kgb==1) spaceComm=mpi_enreg%comm_band
 my_nspinor=max(1,hamk%nspinor/mpi_enreg%nproc_spinor)
 nspinortot=hamk%nspinor

!Check sizes
 npw=hamk%npw_k;if (.not.do_transpose) npw=hamk%npw_fft_k
 if (size(cwavef)/=2*npw*my_nspinor*blocksize) then
   ABI_BUG('Incorrect size for cwavef!')
 end if
 if(choice/=0.and.signs==2) then
   if (paw_opt/=3) then
     if (size(gvnlc)/=2*npw*my_nspinor*blocksize) then
       ABI_BUG('Incorrect size for gvnlc!')
     end if
   end if
   if(paw_opt>=3) then
     if (size(gsc)/=2*npw*my_nspinor*blocksize) then
       ABI_BUG('Incorrect size for gsc!')
     end if
   end if
 end if
 if(cpopt>=0) then
   if (size(cwaveprj)/=hamk%natom*my_nspinor*mpi_enreg%bandpp) then
     ABI_BUG('Incorrect size for cwaveprj!')
   end if
 end if

 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 ikpt_this_proc=bandfft_kpt_get_ikpt()

 recvcounts   => bandfft_kpt(ikpt_this_proc)%recvcounts(:)
 sendcounts   => bandfft_kpt(ikpt_this_proc)%sendcounts(:)
 rdispls      => bandfft_kpt(ikpt_this_proc)%rdispls   (:)
 sdispls      => bandfft_kpt(ikpt_this_proc)%sdispls   (:)
 ndatarecv    =  bandfft_kpt(ikpt_this_proc)%ndatarecv

 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*my_nspinor*bandpp))
 ABI_ALLOCATE(gsc_alltoall2,(2,ndatarecv*my_nspinor*(paw_opt/3)*bandpp))
 ABI_ALLOCATE(gvnlc_alltoall2,(2,ndatarecv*my_nspinor*bandpp))

 if(do_transpose .and. (bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2)))then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
   if(signs==2)then
     if (paw_opt/=3) then
       ABI_ALLOCATE(gvnlc_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
     end if
     if (paw_opt==3.or.paw_opt==4) then
       ABI_ALLOCATE(gsc_alltoall1,(2,ndatarecv*my_nspinor*bandpp))
     end if
   end if
 end if

 ABI_ALLOCATE(enlout,(nnlout*bandpp))
 enlout = zero

 recvcountsloc(:)=recvcounts(:)*2*my_nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*my_nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*my_nspinor
 sdisplsloc(:)=sdispls(:)*2*my_nspinor

 if(do_transpose) then
   call timab(581,1,tsec)
   if(bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2))then
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall1,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   else
     call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
&     recvcountsloc,rdisplsloc,spaceComm,ier)
   end if
   call timab(581,2,tsec)
 else
   ! Here, we cheat, and use DCOPY to bypass some compiler's overzealous bound-checking
   ! (ndatarecv*my_nspinor*bandpp might be greater than the declared size of cwavef)
   call DCOPY(2*ndatarecv*my_nspinor*bandpp,cwavef,1,cwavef_alltoall2,1)
 end if

 if(hamk%istwf_k==2) then
   old_me_g0=mpi_enreg%me_g0
   if (mpi_enreg%me_fft==0) then
     mpi_enreg%me_g0=1
   else
     mpi_enreg%me_g0=0
   end if
 end if

!=====================================================================
 if (bandpp==1) then


   if (do_transpose .and. mpi_enreg%paral_spinor==0.and.nspinortot==2) then !Sort WF by spin
     call prep_sort_wavef_spin(nproc_band,my_nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)
     cwavef_alltoall2(:, :)=cwavef_alltoall1(:,index_wavef_band)
   end if

   if (paw_opt==2) lambda_nonlop(1)=lambdablock(mpi_enreg%me_band+1)
   call nonlop(choice,cpopt,cwaveprj,enlout,hamk,idir,lambda_nonlop,mpi_enreg,1,nnlout,paw_opt,&
&   signs,gsc_alltoall2,tim_nonlop,cwavef_alltoall2,gvnlc_alltoall2)

   if (do_transpose .and. mpi_enreg%paral_spinor == 0 .and. nspinortot==2.and.signs==2) then
     if (paw_opt/=3) gvnlc_alltoall1(:,index_wavef_band)=gvnlc_alltoall2(:,:)
     if (paw_opt==3.or.paw_opt==4) gsc_alltoall1(:,index_wavef_band)=gsc_alltoall2(:,:)
   end if

 else   ! bandpp/=1

!  -------------------------------------------------------------
!  Computation of the index used to sort the waves functions below bandpp
!  -------------------------------------------------------------
   if(do_transpose) then
     call prep_index_wavef_bandpp(nproc_band,bandpp,&
&     my_nspinor,ndatarecv,recvcounts,rdispls,index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the waves functions below bandpp
!  -------------------------------------------------------
     cwavef_alltoall2(:,:) = cwavef_alltoall1(:,index_wavef_band)
   end if

!  -------------------------------------------------------
!  Call nonlop
!  -------------------------------------------------------
   if(paw_opt == 2) lambda_nonlop(1:bandpp) = lambdablock((mpi_enreg%me_band*bandpp)+1:((mpi_enreg%me_band+1)*bandpp))
   call nonlop(choice,cpopt,cwaveprj,enlout,hamk,idir,lambda_nonlop,mpi_enreg,bandpp,nnlout,paw_opt,&
&   signs,gsc_alltoall2,tim_nonlop,cwavef_alltoall2,gvnlc_alltoall2)

!  -----------------------------------------------------
!  Sorting of waves functions below the processors
!  -----------------------------------------------------
   if(do_transpose.and.signs==2) then
     if (paw_opt/=3) gvnlc_alltoall1(:,index_wavef_band)=gvnlc_alltoall2(:,:)
     if (paw_opt==3.or.paw_opt==4) gsc_alltoall1(:,index_wavef_band)=gsc_alltoall2(:,:)
   end if

 end if

!=====================================================================
!  -------------------------------------------------------
!  Deallocation
!  -------------------------------------------------------
 if (allocated(index_wavef_band)) then
   ABI_DEALLOCATE(index_wavef_band)
 end if

!Transpose the gsc_alltoall or gvlnc_alltoall tabs
!according to the paw_opt and signs values
 if(do_transpose) then
   if (signs==2) then
     call timab(581,1,tsec)
     if(bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2))then
       if (paw_opt/=3) then
         call xmpi_alltoallv(gvnlc_alltoall1,recvcountsloc,rdisplsloc,gvnlc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
       if (paw_opt==3.or.paw_opt==4) then
         call xmpi_alltoallv(gsc_alltoall1,recvcountsloc,rdisplsloc,gsc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
     else
       if (paw_opt/=3) then
         call xmpi_alltoallv(gvnlc_alltoall2,recvcountsloc,rdisplsloc,gvnlc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
       if (paw_opt==3.or.paw_opt==4) then
         call xmpi_alltoallv(gsc_alltoall2,recvcountsloc,rdisplsloc,gsc,&
&         sendcountsloc,sdisplsloc,spaceComm,ier)
       end if
     end if
     call timab(581,2,tsec)
   end if
 else
   ! TODO check other usages, maybe
   call DCOPY(2*ndatarecv*my_nspinor*bandpp, gsc_alltoall2, 1, gsc, 1)
 end if
 if (hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0

 if (nnlout>0) then
   call xmpi_allgather(enlout,nnlout*bandpp,enlout_block,spaceComm,ier)
 end if
 ABI_DEALLOCATE(enlout)
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 ABI_DEALLOCATE(gvnlc_alltoall2)
 ABI_DEALLOCATE(gsc_alltoall2)
 if(do_transpose .and. (bandpp/=1 .or. (bandpp==1 .and. mpi_enreg%paral_spinor==0.and.nspinortot==2)))then
   ABI_DEALLOCATE(cwavef_alltoall1)
   if(signs==2)then
     if (paw_opt/=3) then
       ABI_DEALLOCATE(gvnlc_alltoall1)
     end if
     if (paw_opt==3.or.paw_opt==4) then
       ABI_DEALLOCATE(gsc_alltoall1)
     end if
   end if
 end if

 call timab(570,2,tsec)

 DBG_EXIT('COLL')

end subroutine prep_nonlop
!!***

!!****f* ABINIT/prep_fourwf
!! NAME
!! prep_fourwf
!!
!! FUNCTION
!! this routine prepares the data to the call of fourwf.
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cwavef(2,npw*ndat)=planewave coefficients of wavefunction (one spinorial component?).
!!  dtfil <type(datafiles_type)>=variables related to files
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=information about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndat=number of FFT to do in //
!!  ngfft(18)= contain all needed information about 3D FFT
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6 used for dimensionning of vlocal
!!  option_fourwf=option for fourwf (see fourwf.F90)
!!  prtvol=control print volume and debugging output
!!  ucvol=unit cell volume
!!  [bandfft_kpt_tab]= (optional) if present, contains tabs used to implement
!!                     the "band-fft" parallelism
!!                      if not present, the bandfft_kpt global variable is used
!!  [use_gpu_cuda]= (optional) 1 if Cuda (GPU) is on
!!
!! OUTPUT
!!  gwavef=(2,npw*ndat)=matrix elements <G|H|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_mkrho,m_positron,m_vtowfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,iblock,istwf_k,mgfft,&
&          mpi_enreg,nband_k,ndat,ngfft,npw_k,n4,n5,n6,occ_k,option_fourwf,ucvol,wtk,&
&          bandfft_kpt_tab,use_gpu_cuda) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,iblock,istwf_k,mgfft,n4,n5,n6,nband_k,ndat,npw_k
 integer,intent(in) :: option_fourwf
 integer,intent(in),optional :: use_gpu_cuda
 real(dp),intent(in) :: ucvol,wtk
 type(bandfft_kpt_type),optional,target,intent(in) :: bandfft_kpt_tab
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(out) :: rhoaug(n4,n5,n6)
 real(dp),intent(in) :: cwavef(2,npw_k*blocksize)
 real(dp),target,intent(inout) :: wfraug(2,n4,n5,n6*ndat)

!Local variables-------------------------------
!scalars
 integer :: bandpp,bandpp_sym,ier,iibandpp,ikpt_this_proc,ind_occ,ind_occ1,ind_occ2,ipw
 integer :: istwf_k_,jjbandpp,me_fft,nd3,nproc_band,nproc_fft,npw_fft
 integer :: old_me_g0=0,spaceComm=0,tim_fourwf,use_gpu_cuda_
 integer,pointer :: idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 logical :: flag_inv_sym,have_to_reequilibrate
 real(dp) :: weight,weight1,weight2
 type(bandfft_kpt_type),pointer :: bandfft_kpt_ptr
!arrays
 integer,ABI_CONTIGUOUS pointer :: indices_pw_fft(:),kg_k_fft(:,:),kg_k_gather(:,:),kg_k_gather_sym(:,:)
 integer,ABI_CONTIGUOUS pointer :: rdispls(:),rdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: recvcounts(:),recvcount_fft(:),recvcounts_sym(:),recvcounts_sym_tot(:)
 integer,ABI_CONTIGUOUS pointer :: recvdisp_fft(:),sdispls(:),sdispls_sym(:)
 integer,ABI_CONTIGUOUS pointer :: sendcounts(:),sendcount_fft(:),sendcounts_sym(:),sendcounts_sym_all(:)
 integer,ABI_CONTIGUOUS pointer :: senddisp_fft(:),tab_proc(:)
 integer,allocatable :: rdisplsloc(:)
 integer,allocatable :: recvcountsloc(:),sdisplsloc(:)
 integer,allocatable :: sendcountsloc(:)
 integer,allocatable :: index_wavef_band(:),index_wavef_send(:)
 integer,pointer :: gbound_(:,:)
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: buff_wf(:,:),cwavef_alltoall1(:,:),cwavef_alltoall2(:,:)
 real(dp),allocatable :: cwavef_fft(:,:), cwavef_fft_tr(:,:)
 real(dp),allocatable :: weight_t(:),weight1_t(:),weight2_t(:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),wfraug_ptr(:,:,:,:)

! *************************************************************************

 ABI_CHECK((option_fourwf/=3),'Option=3 (FFT r->g) not implemented')
 ABI_CHECK((mpi_enreg%bandpp==ndat),'BUG: bandpp/=ndat')

 spaceComm=mpi_enreg%comm_band
 nproc_band = mpi_enreg%nproc_band
 nproc_fft  = mpi_enreg%nproc_fft
 bandpp     = mpi_enreg%bandpp
 me_fft     = mpi_enreg%me_fft

 use_gpu_cuda_=0;if (present(use_gpu_cuda)) use_gpu_cuda_=use_gpu_cuda

 if (present(bandfft_kpt_tab)) then
   bandfft_kpt_ptr => bandfft_kpt_tab
 else
   ikpt_this_proc=bandfft_kpt_get_ikpt()
   bandfft_kpt_ptr => bandfft_kpt(ikpt_this_proc)
 end if

 istwf_k_=istwf_k
 flag_inv_sym = (istwf_k_==2 .and. any(ngfft(7) == [401,402,312,512]))
 if (option_fourwf==0) flag_inv_sym=((flag_inv_sym).and.(use_gpu_cuda_==0))

 if (flag_inv_sym) then
   istwf_k_       = 1
   if (modulo(bandpp,2)==0) then
     bandpp_sym   = bandpp/2
   else
     bandpp_sym   = bandpp
   end if
 end if

!====================================================================================
 ABI_ALLOCATE(sendcountsloc,(nproc_band))
 ABI_ALLOCATE(sdisplsloc   ,(nproc_band))
 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc   ,(nproc_band))

 recvcounts   =>bandfft_kpt_ptr%recvcounts(:)
 sendcounts   =>bandfft_kpt_ptr%sendcounts(:)
 rdispls      =>bandfft_kpt_ptr%rdispls   (:)
 sdispls      =>bandfft_kpt_ptr%sdispls   (:)
 ndatarecv    =>bandfft_kpt_ptr%ndatarecv

 kg_k_gather  =>bandfft_kpt_ptr%kg_k_gather(:,:)
 gbound_      =>bandfft_kpt_ptr%gbound(:,:)

 if (flag_inv_sym ) then
   idatarecv0           =>bandfft_kpt_ptr%idatarecv0
   ndatarecv_tot        =>bandfft_kpt_ptr%ndatarecv_tot
   ndatasend_sym        =>bandfft_kpt_ptr%ndatasend_sym
   kg_k_gather_sym      =>bandfft_kpt_ptr%kg_k_gather_sym(:,:)
   rdispls_sym          =>bandfft_kpt_ptr%rdispls_sym(:)
   recvcounts_sym       =>bandfft_kpt_ptr%recvcounts_sym(:)
   recvcounts_sym_tot   =>bandfft_kpt_ptr%recvcounts_sym_tot(:)
   sdispls_sym          =>bandfft_kpt_ptr%sdispls_sym(:)
   sendcounts_sym       =>bandfft_kpt_ptr%sendcounts_sym(:)
   sendcounts_sym_all   =>bandfft_kpt_ptr%sendcounts_sym_all(:)
   tab_proc             =>bandfft_kpt_ptr%tab_proc(:)
 end if

 ABI_ALLOCATE(cwavef_alltoall2,(2,ndatarecv*bandpp))
 if ( ((.not.flag_inv_sym) .and. (bandpp>1) ) .or. flag_inv_sym )then
   ABI_ALLOCATE(cwavef_alltoall1,(2,ndatarecv*bandpp))
 end if

 recvcountsloc(:)=recvcounts(:)*2*bandpp
 rdisplsloc(:)=rdispls(:)*2*bandpp
 sendcountsloc(:)=sendcounts(:)*2
 sdisplsloc(:)=sdispls(:)*2

 call timab(547,1,tsec)
 call xmpi_alltoallv(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall2,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(547,2,tsec)

!If me_fft==0, I have the G=0 vector, but keep for the record the old value
 if (me_fft==0) then
   old_me_g0=mpi_enreg%me_g0
   mpi_enreg%me_g0=1
 end if

 tim_fourwf=16

!Eventually adjust load balancing for FFT (by changing FFT distrib)
 have_to_reequilibrate = bandfft_kpt_ptr%have_to_reequilibrate
 if(have_to_reequilibrate) then
   npw_fft =  bandfft_kpt_ptr%npw_fft
   sendcount_fft  => bandfft_kpt_ptr%sendcount_fft(:)
   recvcount_fft  => bandfft_kpt_ptr%recvcount_fft(:)
   senddisp_fft   => bandfft_kpt_ptr%senddisp_fft(:)
   recvdisp_fft   => bandfft_kpt_ptr%recvdisp_fft(:)
   indices_pw_fft => bandfft_kpt_ptr%indices_pw_fft(:)
   kg_k_fft       => bandfft_kpt_ptr%kg_k_fft(:,:)
   ABI_ALLOCATE( buff_wf, (2,ndatarecv*bandpp) ) ! for sorting cgwavef
   ABI_ALLOCATE( cwavef_fft, (2,npw_fft*bandpp) )
   if(bandpp>1) then
     ABI_ALLOCATE( cwavef_fft_tr, (2,npw_fft*bandpp))
   end if
 end if

 if (option_fourwf==0) wfraug(:,:,:,:)=zero

!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

!  Compute the index of the band
   ind_occ = (iblock-1)*blocksize + mpi_enreg%me_band + 1

   if(abs(occ_k(ind_occ))>=tol8.or.option_fourwf==0) then

!    Compute the weight of the band
     weight=occ_k(ind_occ)*wtk/ucvol

     if(have_to_reequilibrate) then
!      filling of sorted send buffers before exchange
       do ipw = 1 ,ndatarecv
         buff_wf(1:2, indices_pw_fft(ipw) ) = cwavef_alltoall2(1:2,ipw)
       end do
       call xmpi_alltoallv(buff_wf,2*sendcount_fft,2*senddisp_fft,  &
&       cwavef_fft,2*recvcount_fft, 2*recvdisp_fft, mpi_enreg%comm_fft,ier)
       call fourwf(1,rhoaug,cwavef_fft,dummy,wfraug,gbound_,gbound_,&
&       istwf_k_,kg_k_fft,kg_k_fft,mgfft,mpi_enreg,1,&
&       ngfft,npw_fft,1,n4,n5,n6,option_fourwf,tim_fourwf,weight,weight,&
&       use_gpu_cuda=use_gpu_cuda_)
     else
       call fourwf(1,rhoaug,cwavef_alltoall2,dummy,wfraug,gbound_,gbound_,&
&       istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&       ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,tim_fourwf,weight,weight,&
&       use_gpu_cuda=use_gpu_cuda_)
     end if
     if (option_fourwf==0.and.nproc_fft>1) then
       if (me_fft>0) then
         nd3=(ngfft(3)-1)/nproc_fft+1
         wfraug(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug(:,:,:,1:nd3)
         wfraug(:,:,:,1:nd3)=zero
       end if
       call xmpi_sum(wfraug,mpi_enreg%comm_fft,ier)
     end if
   end if

!====================================================================
 else if ((.not.(flag_inv_sym)) .and. (bandpp>1) ) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting of the wave functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall1(:,:) = cwavef_alltoall2(:,index_wavef_band)

   if(have_to_reequilibrate) then
!    filling of sorted send buffers before exchange
     do iibandpp=1,bandpp
       do ipw = 1 ,ndatarecv
         buff_wf(1:2, iibandpp + bandpp*(indices_pw_fft(ipw)-1)) = &
&         cwavef_alltoall1(1:2,ipw + ndatarecv*(iibandpp-1))
       end do
     end do
     call xmpi_alltoallv(buff_wf,2*bandpp*sendcount_fft,2*bandpp*senddisp_fft,  &
&     cwavef_fft_tr,2*bandpp*recvcount_fft, 2*bandpp*recvdisp_fft, mpi_enreg%comm_fft,ier)
     do iibandpp=1,bandpp
       do ipw = 1 ,npw_fft
         cwavef_fft(1:2,  ipw + npw_fft*(iibandpp -1)) = cwavef_fft_tr(1:2,  iibandpp + bandpp*(ipw-1))
       end do
     end do
   end if

!  -------------------
!  Fourier calculation
!  -------------------
!  Cuda version
   if(use_gpu_cuda_==1) then
     ABI_ALLOCATE(weight_t,(bandpp))
     do iibandpp=1,bandpp
!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight_t(iibandpp)=occ_k(ind_occ)*wtk/ucvol
       if(abs(occ_k(ind_occ)) < tol8) weight_t(iibandpp) = zero
     end do
!    Accumulate time because it is not done in gpu_fourwf
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     cwavef_alltoall1,&
&     dummy,wfraug,gbound_,gbound_,&
&     istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,bandpp,&
&     ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&     tim_fourwf,weight_t,weight_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight_t)

!  Standard version
   else
     do iibandpp=1,bandpp
!      Compute the index of the band
       ind_occ = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + iibandpp
!      Compute the weight of the band
       weight=occ_k(ind_occ)*wtk/ucvol
       if (option_fourwf==0) then
         wfraug_ptr => wfraug(:,:,:,(iibandpp-1)*n6+1:iibandpp*n6)
       else
         wfraug_ptr => wfraug
       end if
       if (abs(occ_k(ind_occ)) >=tol8.or.option_fourwf==0) then
         if(have_to_reequilibrate) then
           call fourwf(1,rhoaug, &
&           cwavef_fft(:,(npw_fft*(iibandpp-1))+1:(npw_fft*iibandpp)), &
&           dummy,wfraug_ptr,gbound_,gbound_,&
&           istwf_k_,kg_k_fft,kg_k_fft,mgfft,mpi_enreg,1,&
&           ngfft,npw_fft,1,n4,n5,n6,option_fourwf,tim_fourwf,weight,weight,&
&           use_gpu_cuda=use_gpu_cuda_)
         else
           call fourwf(1,rhoaug,&
&           cwavef_alltoall1(:,(ndatarecv*(iibandpp-1))+1:(ndatarecv*iibandpp)),&
&           dummy,wfraug_ptr,gbound_,gbound_,&
&           istwf_k_,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&           ngfft,ndatarecv,1,n4,n5,n6,option_fourwf,&
&           tim_fourwf,weight,weight)
         end if
         if (option_fourwf==0.and.nproc_fft>1) then
           if (me_fft>0) then
             nd3=(ngfft(3)-1)/nproc_fft+1
             wfraug_ptr(:,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug_ptr(:,:,:,1:nd3)
             wfraug_ptr(:,:,:,1:nd3)=zero
           end if
           call xmpi_sum(wfraug_ptr,mpi_enreg%comm_fft,ier)
         end if
       end if
     end do
   end if ! (use_gpu_cuda==1)

!  -----------------------------------------------------
!  Sorting waves functions below the processors
!  -----------------------------------------------------
!  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)   ! NOT NEEDED
   ABI_DEALLOCATE(index_wavef_band)

!====================================================================
 else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Computation of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Sorting the wave functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall1(:,:) = cwavef_alltoall2(:,index_wavef_band)

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,bandpp,1,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,tab_proc,&
&   cwavef_alltoall1,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

!  ------------------------------------------------------------
!  Fourier calculation
!  ------------------------------------------------------------
!  Cuda version
   if (use_gpu_cuda_==1) then
     ABI_ALLOCATE(weight1_t,(bandpp_sym))
     ABI_ALLOCATE(weight2_t,(bandpp_sym))
     do iibandpp=1,bandpp_sym
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if
       weight1_t(iibandpp) = occ_k(ind_occ1)*wtk/ucvol
       weight2_t(iibandpp) = occ_k(ind_occ2)*wtk/ucvol
     end do
     call timab(240+tim_fourwf,1,tsec)
#if defined HAVE_GPU_CUDA
     call gpu_fourwf(1,rhoaug,&
&     ewavef_alltoall_sym,&
&     dummy,wfraug,gbound_,gbound_,&
&     istwf_k_,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,bandpp_sym,&
&     ngfft,ndatarecv_tot,1,n4,n5,n6,option_fourwf,mpi_enreg%paral_kgb,&
&     tim_fourwf,weight1_t,weight2_t)
#endif
     call timab(240+tim_fourwf,2,tsec)
     ABI_DEALLOCATE(weight1_t)
     ABI_DEALLOCATE(weight2_t)

!  Standard version
   else
     if (option_fourwf==0.and.bandpp>1) then
       ABI_ALLOCATE(wfraug_ptr,(2,n4,n5,n6))
     else
       wfraug_ptr => wfraug
     end if
     do iibandpp=1,bandpp_sym
       if (bandpp/=1) then
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp-1)
         ind_occ2 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + (2*iibandpp  )
       else
         ind_occ1 = (iblock-1)*blocksize + (mpi_enreg%me_band * bandpp) + 1
         ind_occ2 = ind_occ1
       end if
       weight1 = occ_k(ind_occ1)*wtk/ucvol
       weight2 = occ_k(ind_occ2)*wtk/ucvol
       call fourwf(1,rhoaug,&
&       ewavef_alltoall_sym(:,(ndatarecv_tot*(iibandpp-1))+1:(ndatarecv_tot*iibandpp)),&
&       dummy,wfraug_ptr,gbound_,gbound_,&
&       istwf_k_,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,1,&
&       ngfft,ndatarecv_tot,1,n4,n5,n6,option_fourwf,&
&       tim_fourwf,weight1,weight2)
       if (option_fourwf==0) then
         if (modulo(bandpp,2)==0) then
           jjbandpp=2*iibandpp-1
           wfraug(1,:,:,(jjbandpp-1)*n6+1:jjbandpp*n6)=wfraug_ptr(1,:,:,1:n6)
           wfraug(1,:,:,(jjbandpp)*n6+1:(jjbandpp+1)*n6)=wfraug_ptr(2,:,:,1:n6)
         else if (bandpp>1) then
           wfraug(1,:,:,(iibandpp-1)*n6+1:iibandpp*n6)=wfraug_ptr(1,:,:,1:n6)
         end if
         if (nproc_fft>1) then
           if (me_fft>0) then
             nd3=(ngfft(3)-1)/nproc_fft+1
             wfraug(1,:,:,me_fft*nd3+1:me_fft*nd3+nd3)=wfraug(1,:,:,1:nd3)
             wfraug(1,:,:,1:nd3)=zero
           end if
           call xmpi_sum(wfraug,mpi_enreg%comm_fft,ier)
         end if
       end if
     end do
     if (option_fourwf==0.and.bandpp>1) then
       ABI_DEALLOCATE(wfraug_ptr)
     end if
   end if ! (use_gpu_cuda==1)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,bandpp,1,&
&   ndatarecv,&
&   ndatarecv_tot,ndatasend_sym,idatarecv0,&
&   cwavef_alltoall1,&
&   sendcounts_sym,sdispls_sym,&
&   recvcounts_sym,rdispls_sym,&
&   ewavef_alltoall_sym,&
&   index_wavef_send)

   ABI_DEALLOCATE(ewavef_alltoall_sym)
   ABI_DEALLOCATE(index_wavef_send)

!  -------------------------------------------------------
!  Sorting waves functions below the processors
!  -------------------------------------------------------
!  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:) ! NOT NEEDED

   ABI_DEALLOCATE(index_wavef_band)

 end if

!====================================================================
 if (me_fft==0) mpi_enreg%me_g0=old_me_g0
 if(have_to_reequilibrate) then
   ABI_DEALLOCATE(buff_wf)
   ABI_DEALLOCATE(cwavef_fft)
   if(bandpp > 1) then
     ABI_DEALLOCATE(cwavef_fft_tr)
   end if
 end if
 ABI_DEALLOCATE(sendcountsloc)
 ABI_DEALLOCATE(sdisplsloc)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)
 ABI_DEALLOCATE(cwavef_alltoall2)
 if ( ((.not.flag_inv_sym) .and. (bandpp>1) ) .or. flag_inv_sym ) then
   ABI_DEALLOCATE(cwavef_alltoall1)
 end if

end subroutine prep_fourwf
!!***

!!****f* ABINIT/prep_wavef_sym_do
!! NAME
!! prep_wavef_sym_do
!!
!! FUNCTION
!! this routine associates waves functions by two as following
!!      E(G)  = C(G) + D(G)
!!      E(-G) = C*(G) + iD*(G)
!! the values are distributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! INPUTS
!!  mpi_enreg          = information about mpi parallelization
!!  bandpp             = number of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!!  cwavef_alltoall    = planewave coefficients of wavefunction
!!                      ( initial of the processor + sended by other processors band)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor from each processor fft
!!  rdispls_sym        = postions of the received values by the processor from each processor fft
!!
!! OUTPUT
!!  ewavef_alltoall_sym = planewave coefficients of wavefunction
!!                        initial of the processor +
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and compisited if bandpp >1
!!  index_wavef_send    = index to send the values in blocks to the other processor fft
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_prep_kgb
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,tab_proc,&
&     cwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     ewavef_alltoall_sym,&
&     index_wavef_send)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,allocatable,intent(out) :: index_wavef_send(:)
 integer,intent(in) :: rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 integer,intent(in) :: tab_proc(:)
 real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),pointer :: ewavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,idatarecv,ideb_loc,idebc,idebd,idebe
 integer :: ier,ifin_loc,ifinc,ifind,ifine,iproc,jbandpp,jsendloc
 integer :: kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: ewavef_alltoall_loc(:,:),ewavef_alltoall_send(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_do : enter '
!ENDDEBUG

!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 newspacecomm = mpi_enreg%comm_fft

 if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
   flag_compose = .TRUE.
 else
   bandpp_sym   = bandpp
   flag_compose = .FALSE.
 end if

!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(ewavef_alltoall_sym     ,(2,ndatarecv_tot*bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_loc     ,(2,ndatarecv    *bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_send    ,(2,ndatasend_sym*bandpp_sym))
 ABI_ALLOCATE(index_wavef_send        ,(  ndatasend_sym*bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))


!Initialisation
!--------------
 ewavef_alltoall_sym(:,:) =0.
 ewavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0

 index_wavef_send(:)   =0


!-------------------------------------------------
!We are bandpp blocks which we want to :
!associate by two      (band_sym==bandpp/2)
!or not associate by two  (band_sym==bandpp)
!
!So We'll have got bandpp_sym blocks
!So we loop on the bandpp_sym blocks
!--------------------------------------------------

 do kbandpp=1,bandpp_sym

!  position of the two blocks
!  --------------------------
   ibandpp = (kbandpp-1) * 2
   jbandpp =  ibandpp    + 1

   idebe = (kbandpp-1) * ndatarecv_tot + 1
   ifine = idebe       + ndatarecv     - 1

   idebc = ibandpp * ndatarecv     + 1
   ifinc = idebc   + ndatarecv     - 1

   idebd = jbandpp * ndatarecv     + 1
   ifind = idebd   + ndatarecv     - 1

   ideb_loc = (kbandpp-1) * ndatarecv  + 1
   ifin_loc = ideb_loc    + ndatarecv  - 1


   if (flag_compose) then


!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine) =    &
&     cwavef_alltoall(1,idebc:ifinc) &
&     - cwavef_alltoall(2,idebd:ifind)

     ewavef_alltoall_sym(2,idebe:ifine) =    &
&     cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)

!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =  &
&     cwavef_alltoall(1,idebc:ifinc) &
&     + cwavef_alltoall(2,idebd:ifind)

     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) =  &
&     - cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)
   else

!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine)   = cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_sym(2,idebe:ifine)   = cwavef_alltoall(2,idebc:ifinc)

!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =   cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) = - cwavef_alltoall(2,idebc:ifinc)

   end if

 end do



!------------------------------------------------------------------------
!Creation of datas blocks for each processor fft from ewavef_alltoall_loc
!to send datas by blocks with a alltoall...
!------------------------------------------------------------------------

!Position of the blocks
!----------------------
 jsendloc=0
 do ibandpp=1,bandpp_sym
   do iproc=1,nproc_fft
     do idatarecv=1,ndatarecv
       if (tab_proc(idatarecv)==(iproc-1)) then
         jsendloc=jsendloc+1
         index_wavef_send(jsendloc)  = idatarecv + ndatarecv * (ibandpp-1)
       end if
     end do
   end do
 end do

!Classment
!----------
 ewavef_alltoall_send(:,:)=ewavef_alltoall_loc(:,index_wavef_send)


!-------------------------------------------------
!Calcul of the number of received and sended datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2

!------------------------------------------
!Exchange of the datas ewavef_allto_all_loc
!------------------------------------------
 do ibandpp=1,bandpp_sym

!  ------------------------------------------------
!  Deplacment of the sended datas because of bandpp
!  ------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2

!  --------------------------------------------------
!  Deplacment of the received datas because of bandpp
!  --------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2


   call xmpi_alltoallv(&
&   ewavef_alltoall_send(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   ewavef_alltoall_sym(:,:)  ,recvcounts_sym_loc,rdispls_sym_loc,&
&   newspacecomm,ier)

 end do

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)

 ABI_DEALLOCATE(ewavef_alltoall_loc)
 ABI_DEALLOCATE(ewavef_alltoall_send)

end subroutine prep_wavef_sym_do
!!***

!!****f* ABINIT/prep_wavef_sym_undo
!! NAME
!! prep_wavef_sym_undo
!!
!! FUNCTION
!! this routine dissociates each wave function in two waves functions as following
!!      C(G) =   ( E*(-G) + E(G))/2
!!      D(G) = i*( E*(-G) - E(G))/2
!! the values are redistributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! INPUTS
!!  mpi_enreg          = information about mpi parallelization
!!  bandpp             = number of groups of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor to each processor fft
!!!  rdispls_sym        = postions of the received values by the processor to each processor fft
!!
!!  gwavef_alltoall_sym = planewave coefficients of wavefunction
!!                        initial of the processor +
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and composited if bandpp >1
!!  index_wavef_send    = index to send the values by block to the other processor fft
!!
!! OUTPUT
!!  gwavef_alltoall     = planewave coefficients of wavefunction
!!                        ( for of the processor + to send to other processors band)
!!
!! PARENTS
!!      m_prep_kgb
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,idatarecv0,&
&     gwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     gwavef_alltoall_sym,&
&     index_wavef_send)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: index_wavef_send(:),rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),intent(inout) :: gwavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,ideb_loc,idebc,idebd
 integer :: idebe,ier,ifin_loc,ifinc,ifind,ifine
 integer :: jbandpp,kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: gwavef_alltoall_loc(:,:),gwavef_alltoall_rcv(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_undo : enter '
!ENDDEBUG


!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 newspacecomm = mpi_enreg%comm_fft

 if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
   flag_compose = .TRUE.
 else
   bandpp_sym   = bandpp
   flag_compose = .FALSE.
 end if

!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(gwavef_alltoall_loc     ,(2,ndatarecv     *bandpp_sym))
 ABI_ALLOCATE(gwavef_alltoall_rcv     ,(2,ndatasend_sym *bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))


!---------------------------------------------
!Initialisation
!---------------------------------------------
 gwavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0


!-------------------------------------------------
!Calcul of number of the sended and received datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2

!----------------------------------------------------
!Sending of the values
!----------------------------------------------------
 do ibandpp = 1,bandpp_sym

!  -------------------------------------------------
!  Deplacment of the sended values because of bandpp
!  -------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2

!  ---------------------------------------------------
!  Deplacment of the received values because of bandpp
!  ---------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2


   call xmpi_alltoallv(&
&   gwavef_alltoall_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,&
&   gwavef_alltoall_rcv(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   newspacecomm,ier)

 end do


!----------------------
!Dispatching the blocks
!----------------------
 gwavef_alltoall_loc(:,index_wavef_send(:)) = gwavef_alltoall_rcv(:,:)

!----------------------
!Case  -kg = [0 0 0]
!----------------------
 if (idatarecv0/=-1) then
   do kbandpp=1,bandpp_sym
     gwavef_alltoall_loc(:,(kbandpp-1)*ndatarecv     + idatarecv0)= &
     gwavef_alltoall_sym(:,(kbandpp-1)*ndatarecv_tot + idatarecv0)
   end do
 end if

!---------------------------------------------------
!Build of hwavef_alltoall
!
!We have got :
!bandpp_sym blocks to dissociate
!or   bandpp_sym blokcs to not dissociate
!--------------------------------------------------
 do kbandpp=1,bandpp_sym

!  position of the 2 blocks
!  ----------------------------------
   ibandpp = (kbandpp-1) * 2
   jbandpp =  ibandpp    + 1

   idebe = (kbandpp-1) * ndatarecv_tot + 1
   ifine = idebe       + ndatarecv     - 1

   idebc = ibandpp * ndatarecv     + 1
   ifinc = idebc   + ndatarecv     - 1

   idebd = jbandpp * ndatarecv     + 1
   ifind = idebd   + ndatarecv     - 1

   ideb_loc = (kbandpp-1) * ndatarecv  + 1
   ifin_loc = ideb_loc    + ndatarecv  - 1


   if (flag_compose) then

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)

!    calcul dwavef(G)
!    ------------------
     gwavef_alltoall(1,idebd:ifind) =   gwavef_alltoall_sym(2,idebe:ifine) &
&     + gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebd:ifind) = - gwavef_alltoall_sym(1,idebe:ifine) &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
   else

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
   end if

 end do

!We divise by two
 gwavef_alltoall(:,:)    = gwavef_alltoall(:,:)/2

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)

 ABI_DEALLOCATE(gwavef_alltoall_loc)
 ABI_DEALLOCATE(gwavef_alltoall_rcv)

end subroutine prep_wavef_sym_undo
!!***

!!****f* ABINIT/prep_index_wavef_bandpp
!! NAME
!! prep_index_wavef_bandpp
!!
!! FUNCTION
!! this routine sorts the waves functions by bandpp and by processors
!! after the alltoall
!!
!! INPUTS
!!  nproc_band = number of processors below the band
!!  bandpp     = number of groups of couple of waves functions
!!  nspinor    = number of spin
!!  ndatarecv  = total number of values received by the processor and sended
!!               by the other processors band
!!  recvcounts = number of values sended by each processor band and received
!!               by the processor
!!  rdispls    = positions of the values received by the processor and
!!               sended by each processor band
!!
!! OUTPUT
!!  index_wavef_band = position of the sorted values
!!
!! PARENTS
!!      m_chebfi,m_prep_kgb
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&
                             nspinor,ndatarecv,&
                             recvcounts,rdispls,&
                             index_wavef_band)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef_band(:)

!Local variables-------------------------------
!scalars
 integer :: delta,idebc,idebe,ifinc,ifine,iindex,iproc,kbandpp,nb

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_index_wavef_banpp : enter '
!write(std_out,*) 'ndatarecv = ', ndatarecv
!write(std_out,*) 'rdispls(:) = ', rdispls(:)
!write(std_out,*) 'recvcounts(:) = ', recvcounts(:)
!ENDDEBUG


!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(index_wavef_band ,(bandpp*nspinor*ndatarecv))
 index_wavef_band(:)   =0

!---------------------------------------------
!Calcul : loops on bandpp and processors band
!---------------------------------------------
 nb = sum(recvcounts(1:nproc_band))
 do kbandpp=1,bandpp

   do iproc=1,nproc_band

     idebe = (rdispls(iproc) + 1)  + (kbandpp-1) * ndatarecv*nspinor
     ifine = idebe + recvcounts(iproc) -1

     if (iproc==1) then
       idebc =   (kbandpp-1)* recvcounts(iproc)*nspinor + 1
     else
       idebc = (bandpp)  * sum(recvcounts(1:iproc-1))*nspinor &
       + (kbandpp-1)* recvcounts(iproc)*nspinor &
       + 1
     end if
     ifinc = idebc + recvcounts(iproc) -1
     index_wavef_band(idebe:ifine) = (/( iindex,iindex=idebc,ifinc)/)
     delta=ifine-idebe
     if (nspinor==2) then
       index_wavef_band(idebe+nb :idebe+nb +delta)=(/( iindex,iindex=ifinc+1,ifinc+1+delta)/)
     end if
   end do
 end do

end subroutine prep_index_wavef_bandpp
!!***

!!****f* ABINIT/prep_sort_wavef_spin
!! NAME
!! prep_sort_wavef_spin
!!
!! FUNCTION
!! Compute index used to sort a spinorial wave-function by spin
!! Sort to have all nspinor=1 fisrt, then all nspinor=2
!!
!! INPUTS
!!  nproc_band=size of "band" communicator
!!  nspinor=number of spinorial components of the wavefunction
!!  ndatarecv=total number of values on all processors
!!  recvcounts(nproc_band)= number of received values by the processor
!!  rdispls(nproc_band)= offsets of the received values by the processor
!!
!! OUTPUT
!!  index_wavef(:)=array containing the sorted indexes (pointer, allocated in this routine)
!!
!! PARENTS
!!      m_prep_kgb
!!
!! CHILDREN
!!
!! SOURCE

subroutine prep_sort_wavef_spin(nproc_band,nspinor,ndatarecv,recvcounts,rdispls,index_wavef)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef(:)

!Local variables-------------------------------
!scalars
 integer :: isft,isft1,iproc,iindex
!arrays
 integer,allocatable :: recvcountsloc(:),rdisplsloc(:)

! *********************************************************************

 ABI_ALLOCATE(index_wavef,(ndatarecv*nspinor))

 ABI_ALLOCATE(recvcountsloc,(nproc_band))
 ABI_ALLOCATE(rdisplsloc,(nproc_band))
 recvcountsloc(:)=recvcounts(:)*2*nspinor
 rdisplsloc(:)=rdispls(:)*2*nspinor

!---------------------------------------------
!Loops on bandpp and processors band
!---------------------------------------------
 isft=0
 do iproc=1,nproc_band

!  ===== Spin up
   if (iproc==1) then
     isft= 0
   else
     isft= sum(recvcounts(1: (iproc-1)))
   end if
   isft1 = 0.5*rdisplsloc(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+recvcounts(iproc))/)

!  =====Spin down
   if (iproc==1)then
     isft=sum(recvcounts(1:nproc_band))
   else
     isft=sum(recvcounts(1:nproc_band)) &
&     +sum(recvcounts(1:iproc-1))
   end if
   isft1 = 0.5 * rdisplsloc(iproc) + recvcounts(iproc)

   index_wavef(1+isft:isft+recvcounts(iproc))= &
&   (/(iindex,iindex=isft1+1,isft1+ recvcounts(iproc))/)

 end do

 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

end subroutine prep_sort_wavef_spin
!!***

end module m_prep_kgb
!!***
