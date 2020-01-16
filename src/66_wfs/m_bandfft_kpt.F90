!!****m* ABINIT/m_bandfft_kpt
!! NAME
!!  m_bandfft_kpt
!!
!! FUNCTION
!!  This module provides the definition of the bandfft_kpt_type
!!  used for kgb parallelization.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (FJ, FB, MT)
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

MODULE m_bandfft_kpt

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes, only : MPI_type
 use m_time,      only : timab
 use m_kg,        only : mkkpg
 use m_fftcore,   only : sphereboundary
 use m_mpinfo,    only : proc_distrb_cycle
 use m_hamiltonian, only : gs_hamiltonian_type

 implicit none

 private
 public :: bandfft_kpt_init1
 public :: bandfft_kpt_init2
 public :: bandfft_kpt_reset
 public :: bandfft_kpt_destroy
 public :: bandfft_kpt_destroy_array
 public :: bandfft_kpt_copy
 public :: bandfft_kpt_mpi_send
 public :: bandfft_kpt_mpi_recv
 public :: bandfft_kpt_savetabs
 public :: bandfft_kpt_restoretabs
 public :: bandfft_kpt_set_ikpt
 public :: bandfft_kpt_get_ikpt
 public :: prep_bandfft_tabs
!!***

!!****t* m_bandfft_kpt/bandfft_kpt_type
!! NAME
!! bandfft_kpt_type
!!
!! FUNCTION
!! The bandfft_kpt_type structured datatype gather different information
!! about the triple band-fft-kpt parallelisation :
!! tabs which are distributed over all the three dimensions and stored during
!! the calculation, dimensions of messages exchange during the calculations...
!! i.e.: all the information which were spread over the entire code before and
!! recomputed at each iline, istep or itime STEP with a large probability to
!! make a mistake.
!!
!! SOURCE

 type, public :: bandfft_kpt_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: flag1_is_allocated             ! determine if the following data are allocated or not
  integer :: npw_tot                        ! array holding the total number of plane waves for each k point
  integer :: ndatarecv                      ! total number of values received by the processor and sent
                                            ! by the other processors band
  integer, allocatable :: kg_k_gather(:,:)  ! planewave coordinates
                                            ! (of the processor + sent by other processors band)
  integer, allocatable :: recvcounts(:)     ! number of values received by the  processor from each processor band
  integer, allocatable :: sendcounts(:)     ! number of values sent   by the  processor to   each processor band
  integer, allocatable :: rdispls   (:)     ! positions of values received by the processor from each processor band
  integer, allocatable :: sdispls   (:)     ! postions of values sent by the processor to each processor band
  integer, allocatable :: gbound(:,:)       ! sphere boundary info: gbound(2*mgfft+8,2)

  integer :: flag2_is_allocated                 ! determine if the following data are allocated or not
  real(dp), allocatable :: ffnl_gather(:,:,:,:) ! ffnl tab (of the processor + sent by other processors band)
  real(dp), allocatable :: kinpw_gather(:)      ! kinpw tab (of the processor + sent by other processors band)
  real(dp), allocatable :: ph3d_gather(:,:,:)   ! ph3d tab (of the processor + sent by other processors band)
  real(dp), allocatable :: kpg_k_gather(:,:)    ! kpg_k tab (of the processor + sent by other processors band)


  integer :: flag3_is_allocated             ! determine if the following data are allocated or not
  integer :: istwf_k                        ! input option parameter that describes the storage of wfs
  integer :: idatarecv0                     ! position of the planewave coordinates (0,0,0)
  integer :: ndatarecv_tot                  ! total number of received values by the processor
                                            ! (ndatarecv   + number of received opposited planewave coordinates)
  integer :: ndatasend_sym                  ! number of sent values to the processors fft to create opposited
                                            ! planewave coordinates
  integer, allocatable :: kg_k_gather_sym(:,:)  ! planewave coordinates
                                            ! (kg_k_gather + opposited planewave coordinates sent by the processors fft)
  integer, allocatable :: rdispls_sym(:)        ! positions of values received by the processor from each processor fft
  integer, allocatable :: recvcounts_sym(:)     ! number of values received by the  processor from each processor fft
  integer, allocatable :: recvcounts_sym_tot(:) ! number of values received by each processor from the  other processors fft
  integer, allocatable :: sdispls_sym(:)        ! postions of values sent by the processor to each processor fft
  integer, allocatable :: sendcounts_sym(:)     ! number of values sent   by the  processor to each processor fft
  integer, allocatable :: sendcounts_sym_all(:) ! number of values sent   by each processor to the other processors fft
  integer, allocatable :: tab_proc(:)           ! positions of opposited planewave coordinates in the list of the processors fft

  logical              :: have_to_reequilibrate ! indicates weather we will have to reequilibrate and allocate all these stuff
  integer              :: npw_fft               ! Number of plane waves during fft step
  integer, allocatable :: indices_pw_fft(:)     !  Indices for sorting pw like pw
  integer, allocatable :: sendcount_fft (:)     ! Number of pw to send to others proc fft
  integer, allocatable :: senddisp_fft(:)       ! Positions for sending
  integer, allocatable :: recvcount_fft(:)      ! Number of pw to receive from others proc fft
  integer, allocatable :: recvdisp_fft(:)       ! Positions for receiving
  integer, allocatable :: kg_k_fft(:,:)         ! planewaves coordinates

 end type bandfft_kpt_type
!!***

 type(bandfft_kpt_type),save,public,pointer :: bandfft_kpt(:) => null()
    ! Contains all the information related to the band/FFT parallelism
    ! which depends on kpt exists only if mpi_enreg%paral_kgb==1

 integer,save,private :: bandfft_kpt_current_ikpt=-1
    ! Index of current k point processed by current proc,
    ! i.e. current index of bandfft_kpt loaded in memory
    !   => bandfft_kpt_current_ikpt=my_kpttab(ikpt)
    !   Please, use bandfft_kpt_set_ikpt method to change the value
    !    and bandfft_kpt_get_ikpt to get it.

CONTAINS

!===========================================================
!!***

!!****f* m_bandfft_kpt/bandfft_kpt_init1
!! NAME
!!  bandfft_kpt_init1
!!
!! FUNCTION
!!  Init all (or part of) scalars and pointers in a bandfft_kpt datastructure
!!
!! INPUTS
!!  istwfk(nkpt)       = input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)    = dimensionless coords of G vecs in basis sphere at k point
!!  mgfft              = maximum single fft dimension (IN)
!!  mkmem              = number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg          = information about MPI parallelization
!!  mpw                = maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol) = number of bands at each k point
!!  nkpt               = number of k points
!!  npwarr(nkpt)       = array holding npw for each k point, taking into account
!!                       the effect of istwfk, and the spreading over processors
!!  nsppol             = 1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!
!!---------------------------------------------------------------------
!!  within the bandfft_kpt data_type : Initialize and compute
!!---------------------------------------------------------------------
!!  gbound             = sphere boundary info
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  istwf_k            = input option parameter that describes the storage of wfs
!!  kg_k_gather        = planewave coordinates
!!                       (of the processor + sended by other processors band)
!!  kg_k_gather_sym    = planewave coordinates
!!                       (kg_k_gather + opposited planewave coordinates sended by the processors fft)
!!  ndatarecv          = total number of values received by the processor and sended
!!                       by the other processors band
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited
!!                       planewave coordinates
!!  ndatarecv_tot      = total number of received values by the processor
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  recvcounts         = number of values received by the  processor from each processor band
!!  recvcounts_sym     = number of values received by the  processor from each processor fft
!!  recvcounts_sym_tot = number of values received by each processor from the  other processors fft
!!  rdispls            = positions of values received by the processor from each processor band
!!  rdispls_sym        = positions of values received by the processor from each processor fft
!!  sendcounts         = number of values sended   by the  processor to   each processor band
!!  sendcounts_sym     = number of values sended   by the  processor to   each processor fft
!!  sendcounts_sym_all = number of values sended   by each processor to the other processors fft
!!  sdispls            = postions of values sended by the processor to each processor band
!!  sdispls_sym        = postions of values sended by the processor to each processor fft
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!! SIDE EFFECTS
!!  bandfft_kpt_in=<type(bandfft_kpt)>=bandfft_kpt datastructure
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_init1(bandfft_kpt_in,istwfk,kg,mgfft,mkmem,mpi_enreg,mpw,nband,nkpt,npwarr,nsppol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mkmem,mpw,nkpt,nsppol
 type(bandfft_kpt_type),pointer :: bandfft_kpt_in(:)
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)

!Local variables-------------------------------
!scalars
 integer :: comm_band,comm_fft,idatarecv,idatarecv0,ierr,ikg,ikpt,ikpt_this_proc,iproc,isppol,istwf_k,itest,jsendloc
 integer :: me_fft,me_kpt,n2,nband_k,ndatarecv,ndatarecv_tot,ndatasend_sym,nproc_band,nproc_fft,npw_fft,npw_k,npw_tot
 logical :: reequilibrate_allocated
 character(len=500)  :: message
!arrays
 integer,allocatable :: buff_kg(:,:),gbound(:,:),indices_pw_fft(:),kg_k(:,:),kg_k_fft(:,:),kg_k_gather(:,:)
 integer,allocatable :: kg_k_gather_all(:,:),kg_k_gather_send(:,:),kg_k_gather_sym(:,:)
 integer,allocatable :: npw_per_proc(:)
 integer,allocatable :: rdispls(:),rdispls_all(:),rdispls_sym(:),rdispls_sym_loc(:)
 integer,allocatable :: recvcounts(:),recvcounts_sym(:),recvcounts_fft(:),recvcounts_sym_loc(:),recvcounts_sym_tot(:)
 integer,allocatable :: recvdisp_fft(:),sdispls(:),sdispls_sym_loc(:),sdispls_sym(:)
 integer,allocatable :: sendcounts(:),sendcounts_fft(:),sendcounts_sym(:),sendcounts_sym_all(:),sendcounts_sym_loc(:)
 integer,allocatable :: senddisp_fft(:),sum_kg(:),tab_proc(:)

! *********************************************************************

 If(mpi_enreg%paral_kgb/=1) then
   nullify(bandfft_kpt_in)
   return
 end if

!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft
 nproc_band   = mpi_enreg%nproc_band

 me_fft       = mpi_enreg%me_fft
 me_kpt       = mpi_enreg%me_kpt

 comm_band    = mpi_enreg%comm_band
 comm_fft     = mpi_enreg%comm_fft

!=============================================================================
!Compute and store various tabs in bandfft_kpt(ikpt) data_struc
!These ones will be used in following subroutines:
!vtorho, mkrho, prep_nonlop, prep_fourwf, prep_getghc...
!=============================================================================

 ABI_DATATYPE_ALLOCATE(bandfft_kpt_in,(mkmem))

 ABI_ALLOCATE(sdispls       ,(nproc_band))
 ABI_ALLOCATE(sendcounts    ,(nproc_band))
 ABI_ALLOCATE(rdispls       ,(nproc_band))
 ABI_ALLOCATE(recvcounts    ,(nproc_band))
 ABI_ALLOCATE(sendcounts_fft,(nproc_fft))
 ABI_ALLOCATE(senddisp_fft  ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_fft,(nproc_fft))
 ABI_ALLOCATE(recvdisp_fft  ,(nproc_fft))

 bandfft_kpt_in(:)%flag1_is_allocated=0
 bandfft_kpt_in(:)%flag2_is_allocated=0
 bandfft_kpt_in(:)%flag3_is_allocated=0

 do isppol=1,nsppol
   ikg=0
   do ikpt=1,nkpt
     npw_k=npwarr(ikpt)
     istwf_k=istwfk(ikpt)
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt))cycle
     ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
     if ((ikpt_this_proc > mkmem).or.(ikpt_this_proc==0)) then
       message = ' this bandfft tab is not allocated !'
       MSG_ERROR(message)
     end if
     if (bandfft_kpt_in(ikpt_this_proc)%flag1_is_allocated==0) then
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%gbound    ,(2*mgfft+8,2))
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%recvcounts,(nproc_band))
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sendcounts,(nproc_band))
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%rdispls   ,(nproc_band))
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sdispls   ,(nproc_band))
       bandfft_kpt_in(ikpt_this_proc)%flag1_is_allocated=1
     end if

!    Initialize various quantities for the reequilibration step
     reequilibrate_allocated=(isppol==2.and.bandfft_kpt_in(ikpt_this_proc)%have_to_reequilibrate)
     bandfft_kpt_in(ikpt_this_proc)%have_to_reequilibrate = .false.
     bandfft_kpt_in(ikpt_this_proc)%npw_fft = 0

     call xmpi_allgather(npw_k,recvcounts,comm_band,ierr)
     rdispls(1)=0
     do iproc=2,nproc_band
       rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
     end do
     ndatarecv=rdispls(nproc_band)+recvcounts(nproc_band)

     ABI_ALLOCATE(kg_k_gather,(3,ndatarecv))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     call xmpi_allgatherv(kg_k,3*npw_k,kg_k_gather,3*recvcounts(:),3*rdispls(:),comm_band,ierr)

     sendcounts(:)=npw_k*mpi_enreg%bandpp
     do iproc=1,nproc_band
       sdispls(iproc)=(iproc-1)*npw_k*mpi_enreg%bandpp
     end do

!    ============================================================================
!    Here we compute gbound, as well for istwf_k=1 as for istwf_k=2 and store it
!    ============================================================================
     !MG: Why don't we compute gbound locally by just computing the full G-sphere from ecut
     ABI_ALLOCATE(npw_per_proc,(nproc_fft))
     ABI_ALLOCATE(rdispls_all,(nproc_fft))
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     if (mgfft>0) gbound(:,:)=0
     if (istwf_k==1) then
       call xmpi_allgather(ndatarecv,npw_per_proc,mpi_enreg%comm_fft,ierr)
       rdispls_all(1)=0
       do iproc=2,nproc_fft
         rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
       end do
       npw_tot=rdispls_all(nproc_fft)+npw_per_proc(nproc_fft)
       ABI_ALLOCATE(kg_k_gather_all,(3,npw_tot))
       call xmpi_allgatherv(kg_k_gather,&
&       3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,mpi_enreg%comm_fft,ierr)
       if (mgfft>0) then
         call sphereboundary(gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
       end if

     else if (istwf_k==2) then

!      ============================================================================
!      In this case, we have to add the opposite values in the kg_k_gather tab
!      before computing gbound
!      ============================================================================

!      Allocation
       ABI_ALLOCATE(tab_proc          ,(ndatarecv))
       ABI_ALLOCATE(sendcounts_sym    ,(nproc_fft))
       ABI_ALLOCATE(sendcounts_sym_all,(nproc_fft*nproc_fft))
       ABI_ALLOCATE(sdispls_sym       ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym    ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym_tot,(nproc_fft))
       ABI_ALLOCATE(rdispls_sym       ,(nproc_fft))
       ABI_ALLOCATE(sendcounts_sym_loc,(nproc_fft))
       ABI_ALLOCATE(sdispls_sym_loc   ,(nproc_fft))
       ABI_ALLOCATE(recvcounts_sym_loc,(nproc_fft))
       ABI_ALLOCATE(rdispls_sym_loc   ,(nproc_fft))

!      Initialisation
       tab_proc(:)            = 0
       sendcounts_sym(:)      = 0
       sendcounts_sym_all(:)  = 0
       sdispls_sym(:)         = 0
       recvcounts_sym(:)      = 0
       recvcounts_sym_tot(:)  = 0

!      Localisation of kg_k==[0 0 0]
       ABI_ALLOCATE(sum_kg,(ndatarecv))
       idatarecv0    = -1
       ndatasend_sym = ndatarecv
       sum_kg=sum(abs(kg_k_gather),1)
       if (count(sum_kg==0)/=0) then
         do idatarecv=1,ndatarecv
           if (sum_kg(idatarecv)==0) idatarecv0=idatarecv
         end do
         ndatasend_sym = ndatarecv-1
       end if

!      Localisation of the processor where the vector -k2 is
       n2 = mpi_enreg%distribfft%n2_coarse
       do idatarecv=1,ndatarecv
         if (idatarecv/=idatarecv0) then
           tab_proc(idatarecv)   = mpi_enreg%distribfft%tab_fftwf2_distrib(modulo(-kg_k_gather(2,idatarecv),n2) + 1)
         else
           tab_proc(idatarecv) = -1
         end if
       end do

!      Number of values send by the processor to the others
       do iproc=1,nproc_fft
         sendcounts_sym(iproc) = count(tab_proc(:)==(iproc-1))
       end do

!      Save sendcounts_sym for each processor in sendcounts_sym_all
!      knowed by all processors of comm_fft
       rdispls_sym(1)=0
       do iproc=2,nproc_fft
         rdispls_sym(iproc)= nproc_fft*(iproc-1)
       end do
       recvcounts_sym(:)=nproc_fft
       call xmpi_allgatherv(sendcounts_sym(:),nproc_fft,&
&        sendcounts_sym_all(:),recvcounts_sym,rdispls_sym,comm_fft,ierr)

!      Calculation of the dimension of kg_k_gather_sym for each processor
!      recvcounts_sym_tot is knowed by all processors of comm_fft
       call xmpi_sum(sendcounts_sym,recvcounts_sym_tot,nproc_fft,comm_fft,ierr)

!      Dimension of kg_k_gather_sym
       ndatarecv_tot = ndatarecv+recvcounts_sym_tot(me_fft+1)

!      Intialize kg_k_gather_sym
       ABI_ALLOCATE(kg_k_gather_sym,(3,ndatarecv_tot))
       kg_k_gather_sym(:,:)=0
       kg_k_gather_sym(:,1:ndatarecv) = kg_k_gather(:,:)

!      Allocation and initialisation
       ABI_ALLOCATE(kg_k_gather_send,(3,ndatasend_sym))
       kg_k_gather_send(:,:)=0

!      The values are sorted in blocks
       jsendloc=0
       do iproc=1,nproc_fft

!        Position of the beginning of the block
         sdispls_sym(iproc)=jsendloc

!        Creation of the blocks
         do idatarecv=1,ndatarecv
           if (tab_proc(idatarecv)==(iproc-1)) then
             jsendloc=jsendloc+1
             kg_k_gather_send(:,jsendloc)  = -kg_k_gather(:,idatarecv)
           end if
         end do
       end do

!      Position of received data
       rdispls_sym(1)= ndatarecv
       recvcounts_sym(1)= sendcounts_sym_all((me_fft+1))
       do iproc=2,nproc_fft
         rdispls_sym(iproc)    = rdispls_sym(iproc-1) + &
         sendcounts_sym_all((me_fft+1)+(iproc-2)*nproc_fft)
         recvcounts_sym(iproc) = sendcounts_sym_all((me_fft+1)+(iproc-1)*nproc_fft)
       end do

!      Exchange of kg_k
       sendcounts_sym_loc = sendcounts_sym*3
       sdispls_sym_loc    = sdispls_sym   *3
       recvcounts_sym_loc = recvcounts_sym*3
       rdispls_sym_loc    = rdispls_sym   *3
       call xmpi_alltoallv(kg_k_gather_send(:,:),sendcounts_sym_loc,sdispls_sym_loc,&
&       kg_k_gather_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,comm_fft,ierr)

!      Store the following data in the bandfft_kpt_in data_struc
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       if (bandfft_kpt_in(ikpt_this_proc)%flag3_is_allocated==0) then
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kg_k_gather_sym,(3,ndatarecv_tot))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%rdispls_sym,(nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%recvcounts_sym,(nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%recvcounts_sym_tot,(nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sdispls_sym,(nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sendcounts_sym,(nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sendcounts_sym_all,(nproc_fft*nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%tab_proc,(ndatarecv))
         bandfft_kpt_in(ikpt_this_proc)%flag3_is_allocated=1
       end if

       bandfft_kpt_in(ikpt_this_proc)%idatarecv0           =idatarecv0
       bandfft_kpt_in(ikpt_this_proc)%ndatarecv_tot        =ndatarecv_tot
       bandfft_kpt_in(ikpt_this_proc)%ndatasend_sym        =ndatasend_sym
       bandfft_kpt_in(ikpt_this_proc)%kg_k_gather_sym(:,:) =kg_k_gather_sym(:,:)
       bandfft_kpt_in(ikpt_this_proc)%rdispls_sym(:)       =rdispls_sym(:)
       bandfft_kpt_in(ikpt_this_proc)%recvcounts_sym(:)    =recvcounts_sym(:)
       bandfft_kpt_in(ikpt_this_proc)%recvcounts_sym_tot(:)=recvcounts_sym_tot(:)
       bandfft_kpt_in(ikpt_this_proc)%sdispls_sym(:)       =sdispls_sym(:)
       bandfft_kpt_in(ikpt_this_proc)%sendcounts_sym(:)    =sendcounts_sym(:)
       bandfft_kpt_in(ikpt_this_proc)%sendcounts_sym_all(:)=sendcounts_sym_all(:)
       bandfft_kpt_in(ikpt_this_proc)%tab_proc(:)          =tab_proc(:)

       ABI_DEALLOCATE(tab_proc)
       ABI_DEALLOCATE(sendcounts_sym)
       ABI_DEALLOCATE(sendcounts_sym_all)
       ABI_DEALLOCATE(sdispls_sym)
       ABI_DEALLOCATE(recvcounts_sym)
       ABI_DEALLOCATE(recvcounts_sym_tot)
       ABI_DEALLOCATE(rdispls_sym)
       ABI_DEALLOCATE(kg_k_gather_sym)
       ABI_DEALLOCATE(sendcounts_sym_loc)
       ABI_DEALLOCATE(recvcounts_sym_loc)
       ABI_DEALLOCATE(sdispls_sym_loc)
       ABI_DEALLOCATE(rdispls_sym_loc)
       ABI_DEALLOCATE(kg_k_gather_send)
       ABI_DEALLOCATE(sum_kg)

!      Then compute gbound
       call xmpi_allgather(ndatarecv_tot,npw_per_proc,mpi_enreg%comm_fft,ierr)
       rdispls_all(1)=0
       do iproc=2,nproc_fft
         rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
       end do
       npw_tot=rdispls_all(nproc_fft)+npw_per_proc(nproc_fft)
       ABI_ALLOCATE(kg_k_gather_all,(3,npw_tot))
       call xmpi_allgatherv(bandfft_kpt_in(ikpt_this_proc)%kg_k_gather_sym,&
&       3*ndatarecv_tot,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,mpi_enreg%comm_fft,ierr)
       if (mgfft>0) then
         call sphereboundary(gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
       end if

!      Only calculations with istwfk=1 or 2
     else
       write(message, '(a,i0,a)' )' the value istwfk=',istwf_k,' is not allowed in case of bandfft parallelization!'
       MSG_BUG(message)
     end if
     ABI_DEALLOCATE(kg_k_gather_all)
     ABI_DEALLOCATE(npw_per_proc)
     ABI_DEALLOCATE(rdispls_all)
!    ============================================================================
!    End of gbound
!    ============================================================================

!    Check if there is k_G vectors have been redistributed (after unbalancing dectecion)
!    Note that FFT load balancing is directly related to unbalancing detection
!    made in kpgsph routine.
     itest=0 ; n2 = mpi_enreg%distribfft%n2_coarse
     do idatarecv=1,ndatarecv
      iproc=mpi_enreg%distribfft%tab_fftwf2_distrib(modulo(kg_k_gather(2,idatarecv),n2)+1)
      if (iproc/=me_fft) itest=itest+1
     end do
     call xmpi_sum(itest,mpi_enreg%comm_fft,ierr)
     if (itest>0) then
       write(message, '(a,i4,3a)' ) &
&        'There is a load unbalancing for the FFT parallelization (kpt',ikpt,').',ch10,&
&        'Plane-wave components will be redistributed before each FFT!'
       MSG_COMMENT(message)
     end if
     bandfft_kpt_in(ikpt_this_proc)%have_to_reequilibrate=(itest>0)

!    If yes, store relevant data
     if(bandfft_kpt_in(ikpt_this_proc)%have_to_reequilibrate) then
       n2 = mpi_enreg%distribfft%n2_coarse
       sendcounts_fft(:) = 0
       recvcounts_fft(:) = 0
       senddisp_fft(:) = 0
       recvdisp_fft(:) = 0
       do idatarecv = 1 ,ndatarecv
          iproc = mpi_enreg%distribfft%tab_fftwf2_distrib(modulo(kg_k_gather(2,idatarecv),n2) + 1)
          sendcounts_fft(iproc + 1) =  sendcounts_fft(iproc + 1)  + 1
       end do
       call xmpi_alltoall(sendcounts_fft,1,recvcounts_fft,1,mpi_enreg%comm_fft,ierr)
       do iproc =2, nproc_fft
          senddisp_fft(iproc) =  senddisp_fft(iproc - 1)  + sendcounts_fft(iproc-1)
          recvdisp_fft(iproc) =  recvdisp_fft(iproc - 1)  + recvcounts_fft(iproc-1)
       end do
       npw_fft = recvdisp_fft(nproc_fft) + recvcounts_fft(nproc_fft) ! nb plane wave for fourwf call
       ABI_ALLOCATE(buff_kg,(3,ndatarecv)) ! for sorting kg_k
       ABI_ALLOCATE(kg_k_fft,(3,npw_fft))
       ABI_ALLOCATE(indices_pw_fft,(ndatarecv))
       !filling of sorted send buffers
       sendcounts_fft(:) = 0
       do idatarecv = 1 ,ndatarecv
          iproc = mpi_enreg%distribfft%tab_fftwf2_distrib(modulo(kg_k_gather(2,idatarecv),n2) + 1)
          sendcounts_fft(iproc + 1) =  sendcounts_fft(iproc + 1)  + 1
          indices_pw_fft(idatarecv) =  senddisp_fft(iproc+1) + sendcounts_fft(iproc+1)
          buff_kg(1:3,indices_pw_fft(idatarecv))  =  kg_k_gather(1:3,idatarecv)
       end do

       call xmpi_alltoallv(buff_kg, 3*sendcounts_fft, 3*senddisp_fft, &
&            kg_k_fft,3*recvcounts_fft, 3*recvdisp_fft, mpi_enreg%comm_fft,ierr)

       if (.not.reequilibrate_allocated) then
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kg_k_fft, (3,npw_fft) )
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%indices_pw_fft, (ndatarecv))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%sendcount_fft, (nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%senddisp_fft, (nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%recvcount_fft, (nproc_fft))
         ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%recvdisp_fft, (nproc_fft))
       end if
       bandfft_kpt_in(ikpt_this_proc)%npw_fft           = npw_fft
       bandfft_kpt_in(ikpt_this_proc)%indices_pw_fft(:) = indices_pw_fft(:)
       bandfft_kpt_in(ikpt_this_proc)%sendcount_fft (:) = sendcounts_fft (:)
       bandfft_kpt_in(ikpt_this_proc)%senddisp_fft(:)   = senddisp_fft(:)
       bandfft_kpt_in(ikpt_this_proc)%recvcount_fft(:)  = recvcounts_fft(:)
       bandfft_kpt_in(ikpt_this_proc)%recvdisp_fft(:)   = recvdisp_fft(:)
       bandfft_kpt_in(ikpt_this_proc)%kg_k_fft(:,:)     = kg_k_fft(:,:)
       ABI_DEALLOCATE(buff_kg )
       ABI_DEALLOCATE(kg_k_fft)
       ABI_DEALLOCATE(indices_pw_fft)
     end if

!    Tabs which are common to istwf_k=1 and 2
     if (.not. allocated(bandfft_kpt_in(ikpt_this_proc)%kg_k_gather)) then
       ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kg_k_gather,(3,ndatarecv))
     end if
     bandfft_kpt_in(ikpt_this_proc)%recvcounts(:)   =recvcounts(:)
     bandfft_kpt_in(ikpt_this_proc)%sendcounts(:)   =sendcounts(:)
     bandfft_kpt_in(ikpt_this_proc)%rdispls(:)      =rdispls(:)
     bandfft_kpt_in(ikpt_this_proc)%sdispls(:)      =sdispls(:)
     bandfft_kpt_in(ikpt_this_proc)%gbound(:,:)     =gbound(:,:)
     bandfft_kpt_in(ikpt_this_proc)%kg_k_gather(:,:)=kg_k_gather(:,:)
     bandfft_kpt_in(ikpt_this_proc)%ndatarecv       =ndatarecv
     bandfft_kpt_in(ikpt_this_proc)%istwf_k         =istwf_k
     bandfft_kpt_in(ikpt_this_proc)%npw_tot         =npw_tot
     ABI_DEALLOCATE(kg_k_gather)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(gbound)

     ikg=ikg+npw_k
   end do
 end do
 ABI_DEALLOCATE(recvcounts)
 ABI_DEALLOCATE(sendcounts)
 ABI_DEALLOCATE(rdispls)
 ABI_DEALLOCATE(sdispls)
 ABI_DEALLOCATE(sendcounts_fft)
 ABI_DEALLOCATE(senddisp_fft)
 ABI_DEALLOCATE(recvcounts_fft)
 ABI_DEALLOCATE(recvdisp_fft)

!=============================================================================
!End of computation and storage of the bandfft_kpt(ikpt) data_struc
!=============================================================================

end subroutine bandfft_kpt_init1
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_init2
!! NAME
!!  bandfft_kpt_init2
!!
!! FUNCTION
!!  Init part of scalars and pointers in a bandfft_kpt datastructure
!!
!! INPUTS
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  kinpw_gather(:)=(modified) kinetic energy for each plane wave (Hartree)
!!  kpg_k_gather(ndatarecv,nkpg)=k+G vector for a given k point
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  ndatarecv= dimension of the arrays
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  ntypat=number of types of atoms in unit cell.
!!  ph3d_gather(2,ndatarecv,matblk)=3-dim structure factors, for each atom and plane wave.

!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_bandfft_tabs
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_init2(bandfft_kpt_in,dimffnl,ffnl_gather,ikpt_this_proc,kinpw_gather,&
&                            kpg_k_gather,lmnmax,matblk,mkmem,ndatarecv,nkpg,ntypat,ph3d_gather)

!Arguments -------------------------------
 integer, intent(in) :: dimffnl,ikpt_this_proc,lmnmax,matblk,mkmem,ndatarecv,nkpg,ntypat
!Local variables-------------------------------
 real(dp),intent(in) :: ffnl_gather(:,:,:,:),kinpw_gather(:)
 real(dp),intent(in) :: kpg_k_gather(:,:),ph3d_gather(:,:,:)
 type(bandfft_kpt_type), intent(inout) :: bandfft_kpt_in (mkmem)

! *********************************************************************

 if (allocated(bandfft_kpt_in(ikpt_this_proc)%ffnl_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ffnl_gather)
 end if
 if (size(ffnl_gather)>0) then
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ffnl_gather,(ndatarecv,dimffnl,lmnmax,ntypat))
   bandfft_kpt_in(ikpt_this_proc)%ffnl_gather(:,:,:,:)=ffnl_gather(:,:,:,:)
 else
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ffnl_gather,(0,0,0,0))
 end if

 if (allocated(bandfft_kpt_in(ikpt_this_proc)%ph3d_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ph3d_gather)
 end if
 if (size(ph3d_gather)>0) then
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ph3d_gather,(2,ndatarecv,matblk))
   bandfft_kpt_in(ikpt_this_proc)%ph3d_gather(:,:,:)  =ph3d_gather(:,:,:)
 else
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%ph3d_gather,(0,0,0))
 end if

 if (allocated(bandfft_kpt_in(ikpt_this_proc)%kpg_k_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kpg_k_gather)
 end if
 if (size(kpg_k_gather)>0) then
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kpg_k_gather,(ndatarecv,nkpg))
   bandfft_kpt_in(ikpt_this_proc)%kpg_k_gather(:,:)   =kpg_k_gather(:,:)
 else
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kpg_k_gather,(0,0))
 end if

 if (allocated(bandfft_kpt_in(ikpt_this_proc)%kinpw_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kinpw_gather)
 end if
 if (size(kinpw_gather)>0) then
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kinpw_gather,(ndatarecv))
   bandfft_kpt_in(ikpt_this_proc)%kinpw_gather(:)     =kinpw_gather(:)
 else
   ABI_ALLOCATE(bandfft_kpt_in(ikpt_this_proc)%kinpw_gather,(0))
 end if

 bandfft_kpt_in(ikpt_this_proc)%flag2_is_allocated=1

end subroutine bandfft_kpt_init2
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_reset
!! NAME
!!  bandfft_kpt_reset
!!
!! FUNCTION
!!  Reset flags a bandfft_kpt datastructure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  bandfft_kpt_in=the datastructure to nullify
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_reset(bandfft_kpt_in)

!Arguments ------------------------------------
 type(bandfft_kpt_type) :: bandfft_kpt_in
!Local variables-------------------------------

! ***********************************************************************

 bandfft_kpt_in%flag1_is_allocated=0
 bandfft_kpt_in%flag2_is_allocated=0
 bandfft_kpt_in%flag3_is_allocated=0
 bandfft_kpt_in%have_to_reequilibrate=.false.

end subroutine bandfft_kpt_reset
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_destroy
!! NAME
!!  bandfft_kpt_destroy
!!
!! FUNCTION
!!  Destroy a bandfft_kpt datastructure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  bandfft_kpt_in=the datastructure to destroy
!!
!! PARENTS
!!      m_bandfft_kpt,posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_destroy(bandfft_kpt_in)

!Arguments ------------------------------------
 type(bandfft_kpt_type) :: bandfft_kpt_in
!Local variables-------------------------------

! ***********************************************************************

 bandfft_kpt_in%flag1_is_allocated=0
 bandfft_kpt_in%flag2_is_allocated=0
 bandfft_kpt_in%flag3_is_allocated=0
 bandfft_kpt_in%have_to_reequilibrate=.false.

 if (allocated(bandfft_kpt_in%kg_k_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in%kg_k_gather)
 end if
 if (allocated(bandfft_kpt_in%gbound)) then
   ABI_DEALLOCATE(bandfft_kpt_in%gbound)
 end if
 if (allocated(bandfft_kpt_in%recvcounts)) then
   ABI_DEALLOCATE(bandfft_kpt_in%recvcounts)
 end if
 if (allocated(bandfft_kpt_in%sendcounts)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sendcounts)
 end if
 if (allocated(bandfft_kpt_in%rdispls)) then
   ABI_DEALLOCATE(bandfft_kpt_in%rdispls)
 end if
 if (allocated(bandfft_kpt_in%sdispls)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sdispls)
 end if
 if (allocated(bandfft_kpt_in%ffnl_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in%ffnl_gather)
 end if
 if (allocated(bandfft_kpt_in%kinpw_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in%kinpw_gather)
 end if
 if (allocated(bandfft_kpt_in%kpg_k_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in%kpg_k_gather)
 end if
 if (allocated(bandfft_kpt_in%ph3d_gather)) then
   ABI_DEALLOCATE(bandfft_kpt_in%ph3d_gather)
 end if
 if (allocated(bandfft_kpt_in%kg_k_gather_sym)) then
   ABI_DEALLOCATE(bandfft_kpt_in%kg_k_gather_sym)
 end if
 if (allocated(bandfft_kpt_in%rdispls_sym)) then
   ABI_DEALLOCATE(bandfft_kpt_in%rdispls_sym)
 end if
 if (allocated(bandfft_kpt_in%recvcounts_sym)) then
   ABI_DEALLOCATE(bandfft_kpt_in%recvcounts_sym)
 end if
 if (allocated(bandfft_kpt_in%recvcounts_sym_tot)) then
   ABI_DEALLOCATE(bandfft_kpt_in%recvcounts_sym_tot)
 end if
 if (allocated(bandfft_kpt_in%sdispls_sym)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sdispls_sym)
 end if
 if (allocated(bandfft_kpt_in%sendcounts_sym)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sendcounts_sym)
 end if
 if (allocated(bandfft_kpt_in%sendcounts_sym_all)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sendcounts_sym_all)
 end if
 if (allocated(bandfft_kpt_in%tab_proc)) then
   ABI_DEALLOCATE(bandfft_kpt_in%tab_proc)
 end if
 if (allocated(bandfft_kpt_in%indices_pw_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%indices_pw_fft)
 end if
 if (allocated(bandfft_kpt_in%sendcount_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%sendcount_fft)
 end if
 if (allocated(bandfft_kpt_in%senddisp_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%senddisp_fft)
 end if
 if (allocated(bandfft_kpt_in%recvcount_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%recvcount_fft)
 end if
 if (allocated(bandfft_kpt_in%recvdisp_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%recvdisp_fft)
 end if
 if (allocated(bandfft_kpt_in%kg_k_fft)) then
   ABI_DEALLOCATE(bandfft_kpt_in%kg_k_fft)
 end if

end subroutine bandfft_kpt_destroy
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_destroy_array
!! NAME
!!  bandfft_kpt_destroy_array
!!
!! FUNCTION
!!  Clean and destroy an array of bandfft_kpt datastructures
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  bandfft_kpt_in(:)=the array of datastructure to destroy
!!
!!
!! PARENTS
!!      gstate,gwls_hamiltonian
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_destroy_array(bandfft_kpt_in,mpi_enreg)

!Arguments ------------------------------------
 type(bandfft_kpt_type),pointer :: bandfft_kpt_in(:)
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: ikpt_this_proc,isppol,ikpt,mkmem,nband,nkpt,nsppol
 character(len=500) :: msg

! ***********************************************************************

 if (xmpi_paral==0) return

 if (associated(bandfft_kpt_in)) then
   mkmem =size(bandfft_kpt_in)
   nkpt=size(mpi_enreg%proc_distrb,1)
   nband=size(mpi_enreg%proc_distrb,2)
   nsppol=size(mpi_enreg%proc_distrb,3)
   if (nsppol==0.or.nkpt==0) then
     msg=' mpi_enreg%proc_distrb should be allocated !'
     MSG_BUG(msg)
   end if
   nkpt=size(mpi_enreg%my_kpttab)
   if (nkpt==0) then
     msg=' mpi_enreg%my_kpttab should be allocated !'
     MSG_BUG(msg)
   end if
   do isppol=1,nsppol
     do ikpt=1,nkpt
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband,isppol,mpi_enreg%me_kpt)) cycle
       ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
       if ((ikpt_this_proc>mkmem) .or.(ikpt_this_proc<=0)) then
         msg=' The bandfft tab cannot be deallocated !'
         MSG_BUG(msg)
       end if
       call bandfft_kpt_destroy(bandfft_kpt_in(ikpt_this_proc))
     end do
   end do
   ABI_DATATYPE_DEALLOCATE(bandfft_kpt_in)
   nullify(bandfft_kpt_in)
 end if

end subroutine bandfft_kpt_destroy_array
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_copy
!! NAME
!!  bandfft_kpt_copy
!!
!! FUNCTION
!!  Copy a bandfft_kpt datastructure into another
!!
!! INPUTS
!!  bandfft_kpt_in=<type(bandfft_kpt_type)>=input bandfft_kpt datastructure
!!
!! OUTPUT
!!  bandfft_kpt_out=<type(bandfft_kpt_type)>=output bandfft_kpt datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_copy(bandfft_kpt_in,bandfft_kpt_out,mpi_enreg1,opt_bandfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: opt_bandfft
 type(bandfft_kpt_type),pointer :: bandfft_kpt_in(:)
 type(bandfft_kpt_type),pointer :: bandfft_kpt_out(:)
 type(MPI_type),intent(inout) :: mpi_enreg1

!Local variables-------------------------------
!scalars
 integer :: ikpt,isppol,jkpt,sz1,sz2,sz3,sz4

! *********************************************************************

!Optional pointers
 if (opt_bandfft==0) then
   nullify(bandfft_kpt_out)
 else if (opt_bandfft==1) then
   if (associated(bandfft_kpt_in)) then
     ABI_DATATYPE_ALLOCATE(bandfft_kpt_out,(size(bandfft_kpt_in)))
     do isppol=1,size(mpi_enreg1%proc_distrb,3)
       do ikpt=1,size(mpi_enreg1%proc_distrb,1)
         sz1=size(mpi_enreg1%proc_distrb,2)
         if(proc_distrb_cycle(mpi_enreg1%proc_distrb,ikpt,1,sz1,isppol,mpi_enreg1%me_kpt)) then
           cycle
         end if
         jkpt=mpi_enreg1%my_kpttab(ikpt)
!          if (allocated(bandfft_kpt_in(jkpt)%ind_kg_mpi_to_seq)) then
!            sz1=size(bandfft_kpt_in(jkpt)%ind_kg_mpi_to_seq)
!            ABI_ALLOCATE(bandfft_kpt_out(jkpt)%ind_kg_mpi_to_seq,(sz1))
!            bandfft_kpt_out(jkpt)%ind_kg_mpi_to_seq= &
! &           bandfft_kpt_in(jkpt)%ind_kg_mpi_to_seq
!          end if
         if (allocated(bandfft_kpt_in(jkpt)%kg_k_gather)) then
           sz1=size(bandfft_kpt_in(jkpt)%kg_k_gather,1)
           sz2=size(bandfft_kpt_in(jkpt)%kg_k_gather,2)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%kg_k_gather,(sz1,sz2))
           bandfft_kpt_out(jkpt)%kg_k_gather= &
&           bandfft_kpt_in(jkpt)%kg_k_gather
         end if
         bandfft_kpt_out(jkpt)%flag1_is_allocated=bandfft_kpt_in(jkpt)%flag1_is_allocated
         if (allocated(bandfft_kpt_in(jkpt)%gbound)) then
           sz1=size(bandfft_kpt_in(jkpt)%gbound,1)
           sz2=size(bandfft_kpt_in(jkpt)%gbound,2)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%gbound,(sz1,sz2))
           bandfft_kpt_out(jkpt)%gbound= &
&           bandfft_kpt_in(jkpt)%gbound
         end if
         if (allocated(bandfft_kpt_in(jkpt)%recvcounts)) then
           sz1=size(bandfft_kpt_in(jkpt)%recvcounts)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%recvcounts,(sz1))
           bandfft_kpt_out(jkpt)%recvcounts= &
&           bandfft_kpt_in(jkpt)%recvcounts
         end if
         if (allocated(bandfft_kpt_in(jkpt)%sendcounts)) then
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sendcounts,(size(bandfft_kpt_in(jkpt)%sendcounts)))
           bandfft_kpt_out(jkpt)%sendcounts= &
&           bandfft_kpt_in(jkpt)%sendcounts
         end if
         if (allocated(bandfft_kpt_in(jkpt)%rdispls)) then
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%rdispls,(size(bandfft_kpt_in(jkpt)%rdispls)))
           bandfft_kpt_out(jkpt)%rdispls= &
&           bandfft_kpt_in(jkpt)%rdispls
         end if
         if (allocated(bandfft_kpt_in(jkpt)%sdispls)) then
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sdispls,(size(bandfft_kpt_in(jkpt)%sdispls)))
           bandfft_kpt_out(jkpt)%sdispls= &
&           bandfft_kpt_in(jkpt)%sdispls
         end if
         bandfft_kpt_out(jkpt)%flag2_is_allocated=bandfft_kpt_in(jkpt)%flag2_is_allocated
         if (allocated(bandfft_kpt_in(jkpt)%ffnl_gather)) then
           sz1=size(bandfft_kpt_in(jkpt)%ffnl_gather,1)
           sz2=size(bandfft_kpt_in(jkpt)%ffnl_gather,2)
           sz3=size(bandfft_kpt_in(jkpt)%ffnl_gather,3)
           sz4=size(bandfft_kpt_in(jkpt)%ffnl_gather,4)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%ffnl_gather,(sz1,sz2,sz3,sz4))
           bandfft_kpt_out(jkpt)%ffnl_gather= &
&           bandfft_kpt_in(jkpt)%ffnl_gather
         end if
         if (allocated(bandfft_kpt_in(jkpt)%kinpw_gather)) then
           sz1=size(bandfft_kpt_in(jkpt)%kinpw_gather)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%kinpw_gather,(sz1))
           bandfft_kpt_out(jkpt)%kinpw_gather= &
&           bandfft_kpt_in(jkpt)%kinpw_gather
         end if
         if (allocated(bandfft_kpt_in(jkpt)%ph3d_gather)) then
           sz1=size(bandfft_kpt_in(jkpt)%ph3d_gather,1)
           sz2=size(bandfft_kpt_in(jkpt)%ph3d_gather,2)
           sz3=size(bandfft_kpt_in(jkpt)%ph3d_gather,3)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%ph3d_gather,(sz1,sz2,sz3))
           bandfft_kpt_out(jkpt)%ph3d_gather= &
&           bandfft_kpt_in(jkpt)%ph3d_gather
         end if
         if (allocated(bandfft_kpt_in(jkpt)%kpg_k_gather)) then
           sz1=size(bandfft_kpt_in(jkpt)%kpg_k_gather,1)
           sz2=size(bandfft_kpt_in(jkpt)%kpg_k_gather,2)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%kpg_k_gather,(sz1,sz2))
           bandfft_kpt_out(jkpt)%kpg_k_gather= &
&           bandfft_kpt_in(jkpt)%kpg_k_gather
         end if
         bandfft_kpt_out(jkpt)%flag3_is_allocated=bandfft_kpt_in(jkpt)%flag3_is_allocated
         if (allocated(bandfft_kpt_in(jkpt)%kg_k_gather_sym)) then
           sz1=size(bandfft_kpt_in(jkpt)%kg_k_gather_sym,1)
           sz2=size(bandfft_kpt_in(jkpt)%kg_k_gather_sym,2)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%kg_k_gather_sym,(sz1,sz2))
           bandfft_kpt_out(jkpt)%kg_k_gather_sym= &
&           bandfft_kpt_in(jkpt)%kg_k_gather_sym
         end if
         if (allocated(bandfft_kpt_in(jkpt)%rdispls_sym)) then
           sz1=size(bandfft_kpt_in(jkpt)%rdispls_sym)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%rdispls_sym,(sz1))
           bandfft_kpt_out(jkpt)%rdispls_sym= &
&           bandfft_kpt_in(jkpt)%rdispls_sym
         end if
         if (allocated(bandfft_kpt_in(jkpt)%recvcounts_sym)) then
           sz1=size(bandfft_kpt_in(jkpt)%recvcounts_sym)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%recvcounts_sym,(sz1))
           bandfft_kpt_out(jkpt)%recvcounts_sym= &
&           bandfft_kpt_in(jkpt)%recvcounts_sym
         end if
         if (allocated(bandfft_kpt_in(jkpt)%recvcounts_sym_tot)) then
           sz1=size(bandfft_kpt_in(jkpt)%recvcounts_sym_tot)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%recvcounts_sym_tot,(sz1))
           bandfft_kpt_out(jkpt)%recvcounts_sym_tot= &
&           bandfft_kpt_in(jkpt)%recvcounts_sym_tot
         end if
         if (allocated(bandfft_kpt_in(jkpt)%sdispls_sym)) then
           sz1=size(bandfft_kpt_in(jkpt)%sdispls_sym)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sdispls_sym,(sz1))
           bandfft_kpt_out(jkpt)%sdispls_sym= &
&           bandfft_kpt_in(jkpt)%sdispls_sym
         end if
         if (allocated(bandfft_kpt_in(jkpt)%sendcounts_sym)) then
           sz1=size(bandfft_kpt_in(jkpt)%sendcounts_sym)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sendcounts_sym,(sz1))
           bandfft_kpt_out(jkpt)%sendcounts_sym= &
&           bandfft_kpt_in(jkpt)%sendcounts_sym
         end if
         if (allocated(bandfft_kpt_in(jkpt)%sendcounts_sym_all)) then
           sz1=size(bandfft_kpt_in(jkpt)%sendcounts_sym_all)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sendcounts_sym_all,(sz1))
           bandfft_kpt_out(jkpt)%sendcounts_sym_all= &
&           bandfft_kpt_in(jkpt)%sendcounts_sym_all
         end if
         if (allocated(bandfft_kpt_in(jkpt)%tab_proc)) then
           sz1=size(bandfft_kpt_in(jkpt)%tab_proc)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%tab_proc,(sz1))
           bandfft_kpt_out(jkpt)%tab_proc= &
           bandfft_kpt_in(jkpt)%tab_proc
         end if
         if (bandfft_kpt_in(jkpt)%have_to_reequilibrate) then
           bandfft_kpt_out(jkpt)%have_to_reequilibrate = .true.
           bandfft_kpt_out(jkpt)%npw_fft = bandfft_kpt_in(jkpt)%npw_fft
           sz1=size(bandfft_kpt_in(jkpt)%indices_pw_fft)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%indices_pw_fft,(sz1))
           bandfft_kpt_out(jkpt)%indices_pw_fft=bandfft_kpt_in(jkpt)%indices_pw_fft
           sz1=size(bandfft_kpt_in(jkpt)%sendcount_fft)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%sendcount_fft,(sz1))
           bandfft_kpt_out(jkpt)%sendcount_fft=bandfft_kpt_in(jkpt)%sendcount_fft
           sz1=size(bandfft_kpt_in(jkpt)%senddisp_fft)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%senddisp_fft,(sz1))
           bandfft_kpt_out(jkpt)%senddisp_fft=bandfft_kpt_in(jkpt)%senddisp_fft
           sz1=size(bandfft_kpt_in(jkpt)%recvcount_fft)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%recvcount_fft,(sz1))
           bandfft_kpt_out(jkpt)%recvcount_fft=bandfft_kpt_in(jkpt)%recvcount_fft
           sz1=size(bandfft_kpt_in(jkpt)%recvdisp_fft)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%recvdisp_fft,(sz1))
           bandfft_kpt_out(jkpt)%recvdisp_fft=bandfft_kpt_in(jkpt)%recvdisp_fft
           sz1=size(bandfft_kpt_in(jkpt)%kg_k_fft,1)
           sz2=size(bandfft_kpt_in(jkpt)%kg_k_fft,2)
           ABI_ALLOCATE(bandfft_kpt_out(jkpt)%kg_k_fft,(sz1,sz2))
           bandfft_kpt_out(jkpt)%kg_k_fft=bandfft_kpt_in(jkpt)%kg_k_fft
         end if
       end do
     end do
   end if
 end if

end subroutine bandfft_kpt_copy
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/bandfft_kpt_mpi_send
!! NAME
!! bandfft_kpt_mpi_send
!!
!! FUNCTION
!! Send a bandfft_kpt_type inside a MPI communicator
!!
!! INPUTS
!!  input=The datatype to be transmitted
!!  receiver=ID of the receiver in spaceComm
!!  tag=message tag
!!  spaceComm=MPI Communicator
!!  [profile]=(character, optional, default=full)
!!            Define the part of the datastructure to be sent; possible values:
!!            'full'= the entire datastructure
!!            'fourwf'= only arrays needed by the fourwf routine
!!
!! OUTPUT
!!  ierr=Error status
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_mpi_send(input,receiver,tag,spaceComm,ierr,profile)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: receiver,spaceComm,tag
 integer,intent(out) :: ierr
 character(len=*),optional,intent(in) :: profile
!arrays
 type(bandfft_kpt_type),intent(in) :: input

!Local variables-------------------------------
!scalars
 integer :: ipck,nsize,size_dp,size_int
 integer :: size1_kg_k_gather=0,size2_kg_k_gather=0,size_recvcounts=0
 integer :: size_sendcounts=0,size_rdispls=0,size_sdispls=0,size1_gbound=0,size2_gbound
 integer :: size1_kg_k_gather_sym=0,size2_kg_k_gather_sym=0,size_rdispls_sym=0
 integer :: size_sdispls_sym=0,size_recvcounts_sym=0,size_recvcounts_sym_tot=0
 integer :: size_sendcounts_sym=0,size_sendcounts_sym_all=0,size_tab_proc=0
 integer :: size_indices_pw_fft=0,size_sendcount_fft=0,size_senddisp_fft=0
 integer :: size_recvcount_fft=0,size_recvdisp_fft=0,size1_kg_k_fft=0,size2_kg_k_fft=0
 integer :: size1_ffnl_gather=0,size2_ffnl_gather=0,size3_ffnl_gather=0,size4_ffnl_gather=0
 integer :: size_kinpw_gather=0,size1_ph3d_gather=0,size2_ph3d_gather=0,size3_ph3d_gather=0
 integer :: size1_kpg_k_gather=0,size2_kpg_k_gather=0
 logical :: fourwf,full
!arrays
 integer,allocatable :: buffer_int(:)
 real(dp),allocatable :: buffer_dp(:)

! *************************************************************************

 if (xmpi_comm_size(spaceComm)<=1) return

 fourwf=.false.;if (present(profile)) fourwf=(trim(profile)=='fourwf')
 full=(.not.fourwf)

!=== Store sizes ====
 if (fourwf.or.full) then
   if (allocated(input%kg_k_gather)) size1_kg_k_gather=size(input%kg_k_gather,1)
   if (allocated(input%kg_k_gather)) size2_kg_k_gather=size(input%kg_k_gather,2)
   if (input%flag1_is_allocated==1) then
     if (allocated(input%recvcounts)) size_recvcounts=size(input%recvcounts)
     if (allocated(input%sendcounts)) size_sendcounts=size(input%sendcounts)
     if (allocated(input%rdispls)) size_rdispls=size(input%rdispls)
     if (allocated(input%sdispls)) size_sdispls=size(input%sdispls)
     if (allocated(input%gbound)) size1_gbound=size(input%gbound,1)
     if (allocated(input%gbound)) size2_gbound=size(input%gbound,2)
   end if
   if (input%flag3_is_allocated==1.and.input%istwf_k>1) then
     if (allocated(input%kg_k_gather_sym)) size1_kg_k_gather_sym=size(input%kg_k_gather_sym,1)
     if (allocated(input%kg_k_gather_sym)) size2_kg_k_gather_sym=size(input%kg_k_gather_sym,2)
     if (allocated(input%rdispls_sym)) size_rdispls_sym=size(input%rdispls_sym)
     if (allocated(input%sdispls_sym)) size_sdispls_sym=size(input%sdispls_sym)
     if (allocated(input%recvcounts_sym)) size_recvcounts_sym=size(input%recvcounts_sym)
     if (allocated(input%recvcounts_sym_tot)) size_recvcounts_sym_tot=size(input%recvcounts_sym_tot)
     if (allocated(input%sendcounts_sym)) size_sendcounts_sym=size(input%sendcounts_sym)
     if (allocated(input%sendcounts_sym_all)) size_sendcounts_sym_all=size(input%sendcounts_sym_all)
     if (allocated(input%tab_proc)) size_tab_proc=size(input%tab_proc)
   end if
   if (input%have_to_reequilibrate) then
     if (allocated(input%indices_pw_fft)) size_indices_pw_fft=size(input%indices_pw_fft)
     if (allocated(input%sendcount_fft)) size_sendcount_fft=size(input%sendcount_fft)
     if (allocated(input%senddisp_fft)) size_senddisp_fft=size(input%senddisp_fft)
     if (allocated(input%recvcount_fft)) size_recvcount_fft=size(input%recvcount_fft)
     if (allocated(input%recvdisp_fft)) size_recvdisp_fft=size(input%recvdisp_fft)
     if (allocated(input%kg_k_fft)) size1_kg_k_fft=size(input%kg_k_fft,1)
     if (allocated(input%kg_k_fft)) size2_kg_k_fft=size(input%kg_k_fft,2)
   end if
 end if
 if (input%flag2_is_allocated==1.and.full) then
   if (allocated(input%ffnl_gather)) size1_ffnl_gather=size(input%ffnl_gather,1)
   if (allocated(input%ffnl_gather)) size2_ffnl_gather=size(input%ffnl_gather,2)
   if (allocated(input%ffnl_gather)) size3_ffnl_gather=size(input%ffnl_gather,3)
   if (allocated(input%ffnl_gather)) size4_ffnl_gather=size(input%ffnl_gather,4)
   if (allocated(input%kinpw_gather)) size_kinpw_gather=size(input%kinpw_gather)
   if (allocated(input%ph3d_gather)) size1_ph3d_gather=size(input%ph3d_gather,1)
   if (allocated(input%ph3d_gather)) size2_ph3d_gather=size(input%ph3d_gather,2)
   if (allocated(input%ph3d_gather)) size3_ph3d_gather=size(input%ph3d_gather,3)
   if (allocated(input%kpg_k_gather)) size1_kpg_k_gather=size(input%kpg_k_gather,1)
   if (allocated(input%kpg_k_gather)) size2_kpg_k_gather=size(input%kpg_k_gather,2)
 end if

!=== Compute amount of data to transmit ====
 size_int=size1_kg_k_gather*size2_kg_k_gather &
&        +size_recvcounts+size_sendcounts+size_rdispls+size_sdispls &
&        +size1_gbound*size2_gbound &
&        +size1_kg_k_gather_sym*size2_kg_k_gather_sym &
&        +size_rdispls_sym+size_sdispls_sym &
&        +size_recvcounts_sym+size_recvcounts_sym_tot &
&        +size_sendcounts_sym+size_sendcounts_sym_all+size_tab_proc &
&        +size_indices_pw_fft+size_sendcount_fft+size_senddisp_fft &
&        +size_recvcount_fft+size_recvdisp_fft &
&        +size1_kg_k_fft*size2_kg_k_fft
 size_dp=size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather &
&       +size_kinpw_gather &
&       +size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather &
&       +size1_kpg_k_gather*size2_kpg_k_gather

!=== Send flags and array sizes ====
 nsize=45
 ABI_ALLOCATE(buffer_int,(nsize))
 buffer_int(1)=input%flag1_is_allocated
 buffer_int(2)=input%flag2_is_allocated
 buffer_int(3)=input%flag3_is_allocated
 buffer_int(4)=0;if (input%have_to_reequilibrate) buffer_int(4)=1
 buffer_int(5)=input%istwf_k
 buffer_int(6)=input%npw_tot
 buffer_int(7)=input%npw_fft
 buffer_int(8)=input%ndatarecv
 buffer_int(9)=input%idatarecv0
 buffer_int(10)=input%ndatarecv_tot
 buffer_int(11)=input%ndatasend_sym
 buffer_int(12)=size_recvcounts
 buffer_int(13)=size_sendcounts
 buffer_int(14)=size_rdispls
 buffer_int(15)=size_sdispls
 buffer_int(16)=size1_gbound
 buffer_int(17)=size2_gbound
 buffer_int(18)=size1_ffnl_gather
 buffer_int(19)=size2_ffnl_gather
 buffer_int(20)=size3_ffnl_gather
 buffer_int(21)=size4_ffnl_gather
 buffer_int(22)=size_kinpw_gather
 buffer_int(23)=size1_ph3d_gather
 buffer_int(24)=size2_ph3d_gather
 buffer_int(25)=size3_ph3d_gather
 buffer_int(26)=size1_kpg_k_gather
 buffer_int(27)=size2_kpg_k_gather
 buffer_int(28)=size_indices_pw_fft
 buffer_int(29)=size_sendcount_fft
 buffer_int(30)=size_senddisp_fft
 buffer_int(31)=size_recvcount_fft
 buffer_int(32)=size_recvdisp_fft
 buffer_int(33)=size1_kg_k_fft
 buffer_int(34)=size2_kg_k_fft
 buffer_int(35)=size1_kg_k_gather
 buffer_int(36)=size2_kg_k_gather
 buffer_int(37)=size1_kg_k_gather_sym
 buffer_int(38)=size2_kg_k_gather_sym
 buffer_int(39)=size_rdispls_sym
 buffer_int(40)=size_sdispls_sym
 buffer_int(41)=size_recvcounts_sym
 buffer_int(42)=size_recvcounts_sym_tot
 buffer_int(43)=size_sendcounts_sym
 buffer_int(44)=size_sendcounts_sym_all
 buffer_int(45)=size_tab_proc
 call xmpi_send(buffer_int,receiver,3*tag-2,spaceComm,ierr)
 ABI_DEALLOCATE(buffer_int)

!=== Pack and send integer data ====
 if (size_int>0) then
   ipck=0
   ABI_ALLOCATE(buffer_int,(size_int))
   if (size1_kg_k_gather*size2_kg_k_gather>0) then
     nsize=size1_kg_k_gather*size2_kg_k_gather
     buffer_int(ipck+1:ipck+nsize)=reshape(input%kg_k_gather(:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   if (size_recvcounts>0) then
     buffer_int(ipck+1:ipck+size_recvcounts)=input%recvcounts(:)
     ipck=ipck+size_recvcounts
   end if
   if (size_sendcounts>0) then
     buffer_int(ipck+1:ipck+size_sendcounts)=input%sendcounts(:)
     ipck=ipck+size_sendcounts
   end if
   if (size_rdispls>0) then
     buffer_int(ipck+1:ipck+size_rdispls)=input%rdispls(:)
     ipck=ipck+size_rdispls
   end if
   if (size_sdispls>0) then
     buffer_int(ipck+1:ipck+size_sdispls)=input%sdispls(:)
     ipck=ipck+size_sdispls
   end if
   if (size1_gbound*size2_gbound>0) then
     nsize=size1_gbound*size2_gbound
     buffer_int(ipck+1:ipck+nsize)=reshape(input%gbound(:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   if (size1_kg_k_gather_sym*size2_kg_k_gather_sym>0) then
     nsize=size1_kg_k_gather_sym*size2_kg_k_gather_sym
     buffer_int(ipck+1:ipck+nsize)=reshape(input%kg_k_gather_sym(:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   if (size_rdispls_sym>0) then
     buffer_int(ipck+1:ipck+size_rdispls_sym)=input%rdispls_sym(:)
     ipck=ipck+size_rdispls_sym
   end if
   if (size_sdispls_sym>0) then
     buffer_int(ipck+1:ipck+size_sdispls_sym)=input%sdispls_sym(:)
     ipck=ipck+size_sdispls_sym
   end if
   if (size_recvcounts_sym>0) then
     buffer_int(ipck+1:ipck+size_recvcounts_sym)=input%recvcounts_sym(:)
     ipck=ipck+size_recvcounts_sym
   end if
   if (size_recvcounts_sym_tot>0) then
     buffer_int(ipck+1:ipck+size_recvcounts_sym_tot)=input%recvcounts_sym_tot(:)
     ipck=ipck+size_recvcounts_sym_tot
   end if
   if (size_sendcounts_sym>0) then
     buffer_int(ipck+1:ipck+size_sendcounts_sym)=input%sendcounts_sym(:)
     ipck=ipck+size_sendcounts_sym
   end if
   if (size_sendcounts_sym_all>0) then
     buffer_int(ipck+1:ipck+size_sendcounts_sym_all)=input%sendcounts_sym_all(:)
     ipck=ipck+size_sendcounts_sym_all
   end if
   if (size_tab_proc>0) then
     buffer_int(ipck+1:ipck+size_tab_proc)=input%tab_proc(:)
     ipck=ipck+size_tab_proc
   end if
   if (size_indices_pw_fft>0) then
     buffer_int(ipck+1:ipck+size_indices_pw_fft)=input%indices_pw_fft(:)
     ipck=ipck+size_indices_pw_fft
   end if
   if (size_sendcount_fft>0) then
     buffer_int(ipck+1:ipck+size_sendcount_fft)=input%sendcount_fft(:)
     ipck=ipck+size_sendcount_fft
   end if
   if (size_senddisp_fft>0) then
     buffer_int(ipck+1:ipck+size_senddisp_fft)=input%senddisp_fft(:)
     ipck=ipck+size_senddisp_fft
   end if
   if (size_recvcount_fft>0) then
     buffer_int(ipck+1:ipck+size_recvcount_fft)=input%recvcount_fft(:)
     ipck=ipck+size_recvcount_fft
   end if
   if (size_recvdisp_fft>0) then
     buffer_int(ipck+1:ipck+size_recvdisp_fft)=input%recvdisp_fft(:)
     ipck=ipck+size_recvdisp_fft
   end if
   if (size1_kg_k_fft*size2_kg_k_fft>0) then
     nsize=size1_kg_k_fft*size2_kg_k_fft
     buffer_int(ipck+1:ipck+nsize)=reshape(input%kg_k_fft(:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   call xmpi_send(buffer_int,receiver,3*tag-1,spaceComm,ierr)
   ABI_DEALLOCATE(buffer_int)
 end if

!=== Pack and send real data ====
 if (size_dp>0) then
   ipck=0
   ABI_ALLOCATE(buffer_dp,(size_dp))
   if (size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather>0) then
     nsize=size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather
     buffer_dp(ipck+1:ipck+nsize)=reshape(input%ffnl_gather(:,:,:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   if (size_kinpw_gather>0) then
     buffer_dp(ipck+1:ipck+size_kinpw_gather)=input%kinpw_gather(:)
     ipck=ipck+size_kinpw_gather
   end if
   if (size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather>0) then
     nsize=size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather
     buffer_dp(ipck+1:ipck+nsize)=reshape(input%ph3d_gather(:,:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   if (size1_kpg_k_gather*size2_kpg_k_gather>0) then
     nsize=size1_kpg_k_gather*size2_kpg_k_gather
     buffer_dp(ipck+1:ipck+nsize)=reshape(input%kpg_k_gather(:,:),(/nsize/))
     ipck=ipck+nsize
   end if
   call xmpi_send(buffer_dp,receiver,3*tag,spaceComm,ierr)
   ABI_DEALLOCATE(buffer_dp)
 end if

end subroutine bandfft_kpt_mpi_send
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/bandfft_kpt_mpi_recv
!! NAME
!! bandfft_kpt_mpi_recv
!!
!! FUNCTION
!! Receive a bandfft_kpt_type inside a MPI communicator.
!!
!! INPUTS
!!  sender=ID of the sender in spaceComm
!!  tag=message tag
!!  spaceComm=MPI Communicator
!!
!! OUTPUT
!!  ierr=Error status
!!  output=# of on proc. sender
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_mpi_recv(output,sender,tag,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sender,spaceComm,tag
 integer,intent(out) :: ierr
!arrays
 type(bandfft_kpt_type),intent(out) :: output

!Local variables-------------------------------
!scalars
 integer :: ipck,nsize,size_dp,size_int
 integer :: size1_kg_k_gather=0,size2_kg_k_gather=0,size_recvcounts=0
 integer :: size_sendcounts=0,size_rdispls=0,size_sdispls=0,size1_gbound=0,size2_gbound
 integer :: size1_kg_k_gather_sym=0,size2_kg_k_gather_sym=0,size_rdispls_sym=0
 integer :: size_sdispls_sym=0,size_recvcounts_sym=0,size_recvcounts_sym_tot=0
 integer :: size_sendcounts_sym=0,size_sendcounts_sym_all=0,size_tab_proc=0
 integer :: size_indices_pw_fft=0,size_sendcount_fft=0,size_senddisp_fft=0
 integer :: size_recvcount_fft=0,size_recvdisp_fft=0,size1_kg_k_fft=0,size2_kg_k_fft=0
 integer :: size1_ffnl_gather=0,size2_ffnl_gather=0,size3_ffnl_gather=0,size4_ffnl_gather=0
 integer :: size_kinpw_gather=0,size1_ph3d_gather=0,size2_ph3d_gather=0,size3_ph3d_gather=0
 integer :: size1_kpg_k_gather=0,size2_kpg_k_gather=0,sz1,sz2,sz3,sz4
!arrays
 integer,allocatable :: buffer_int(:)
 real(dp),allocatable :: buffer_dp(:)

! *************************************************************************

 if (xmpi_comm_size(spaceComm)<=1) return

!=== Receive flags and array sizes ====
 nsize=45
 ABI_ALLOCATE(buffer_int,(nsize))
 call xmpi_recv(buffer_int,sender,3*tag-2,spaceComm,ierr)
 output%flag1_is_allocated=buffer_int(1)
 output%flag2_is_allocated=buffer_int(2)
 output%flag3_is_allocated=buffer_int(3)
 output%have_to_reequilibrate=(buffer_int(4)/=0)
 output%istwf_k=buffer_int(5)
 output%npw_tot=buffer_int(6)
 output%npw_fft=buffer_int(7)
 output%ndatarecv=buffer_int(8)
 output%idatarecv0=buffer_int(9)
 output%ndatarecv_tot=buffer_int(10)
 output%ndatasend_sym=buffer_int(11)
 size_recvcounts=buffer_int(12)
 size_sendcounts=buffer_int(13)
 size_rdispls=buffer_int(14)
 size_sdispls=buffer_int(15)
 size1_gbound=buffer_int(16)
 size2_gbound=buffer_int(17)
 size1_ffnl_gather=buffer_int(18)
 size2_ffnl_gather=buffer_int(19)
 size3_ffnl_gather=buffer_int(20)
 size4_ffnl_gather=buffer_int(21)
 size_kinpw_gather=buffer_int(22)
 size1_ph3d_gather=buffer_int(23)
 size2_ph3d_gather=buffer_int(24)
 size3_ph3d_gather=buffer_int(25)
 size1_kpg_k_gather=buffer_int(26)
 size2_kpg_k_gather=buffer_int(27)
 size_indices_pw_fft=buffer_int(28)
 size_sendcount_fft=buffer_int(29)
 size_senddisp_fft=buffer_int(30)
 size_recvcount_fft=buffer_int(31)
 size_recvdisp_fft=buffer_int(32)
 size1_kg_k_fft=buffer_int(33)
 size2_kg_k_fft=buffer_int(34)
 size1_kg_k_gather=buffer_int(35)
 size2_kg_k_gather=buffer_int(36)
 size1_kg_k_gather_sym=buffer_int(37)
 size2_kg_k_gather_sym=buffer_int(38)
 size_rdispls_sym=buffer_int(39)
 size_sdispls_sym=buffer_int(40)
 size_recvcounts_sym=buffer_int(41)
 size_recvcounts_sym_tot=buffer_int(42)
 size_sendcounts_sym=buffer_int(43)
 size_sendcounts_sym_all=buffer_int(44)
 size_tab_proc=buffer_int(45)
 ABI_DEALLOCATE(buffer_int)

!=== Compute amount of transmitted data ====
 size_int=size1_kg_k_gather*size2_kg_k_gather &
&        +size_recvcounts+size_sendcounts+size_rdispls+size_sdispls &
&        +size1_gbound*size2_gbound &
&        +size1_kg_k_gather_sym*size2_kg_k_gather_sym &
&        +size_rdispls_sym+size_sdispls_sym &
&        +size_recvcounts_sym+size_recvcounts_sym_tot &
&        +size_sendcounts_sym+size_sendcounts_sym_all+size_tab_proc &
&        +size_indices_pw_fft+size_sendcount_fft+size_senddisp_fft &
&        +size_recvcount_fft+size_recvdisp_fft &
&        +size1_kg_k_fft*size2_kg_k_fft
 size_dp=size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather &
&       +size_kinpw_gather &
&       +size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather &
&       +size1_kpg_k_gather*size2_kpg_k_gather

!=== Receive and unpack integer data ====
 if (size_int>0) then
   ipck=0
   ABI_ALLOCATE(buffer_int,(size_int))
   call xmpi_recv(buffer_int,sender,3*tag-1,spaceComm,ierr)
   if (allocated(output%kg_k_gather)) then
     ABI_DEALLOCATE(output%kg_k_gather)
   end if
   if (size1_kg_k_gather*size2_kg_k_gather>0) then
     nsize=size1_kg_k_gather*size2_kg_k_gather
     sz1=size1_kg_k_gather;sz2=size2_kg_k_gather
     ABI_ALLOCATE(output%kg_k_gather,(sz1,sz2))
     output%kg_k_gather(:,:)=reshape(buffer_int(ipck+1:ipck+nsize),(/sz1,sz2/))
     ipck=ipck+nsize
   end if
   if (allocated(output%recvcounts)) then
     ABI_DEALLOCATE(output%recvcounts)
   end if
   if (size_recvcounts>0) then
     ABI_ALLOCATE(output%recvcounts,(size_recvcounts))
     output%recvcounts(:)=buffer_int(ipck+1:ipck+size_recvcounts)
     ipck=ipck+size_recvcounts
   end if
   if (allocated(output%sendcounts)) then
     ABI_DEALLOCATE(output%sendcounts)
   end if
   if (size_sendcounts>0) then
     ABI_ALLOCATE(output%sendcounts,(size_sendcounts))
     output%sendcounts(:)=buffer_int(ipck+1:ipck+size_sendcounts)
     ipck=ipck+size_sendcounts
   end if
   if (allocated(output%rdispls)) then
     ABI_DEALLOCATE(output%rdispls)
   end if
   if (size_rdispls>0) then
     ABI_ALLOCATE(output%rdispls,(size_rdispls))
     output%rdispls(:)=buffer_int(ipck+1:ipck+size_rdispls)
     ipck=ipck+size_rdispls
   end if
   if (allocated(output%sdispls)) then
     ABI_DEALLOCATE(output%sdispls)
   end if
   if (size_sdispls>0) then
     ABI_ALLOCATE(output%sdispls,(size_sdispls))
     output%sdispls(:)=buffer_int(ipck+1:ipck+size_sdispls)
     ipck=ipck+size_sdispls
   end if
   if (allocated(output%gbound)) then
     ABI_DEALLOCATE(output%gbound)
   end if
   if (size1_gbound*size2_gbound>0) then
     nsize=size1_gbound*size2_gbound
     sz1=size1_gbound;sz2=size2_gbound
     ABI_ALLOCATE(output%gbound,(sz1,sz2))
     output%gbound(:,:)=reshape(buffer_int(ipck+1:ipck+nsize),(/sz1,sz2/))
     ipck=ipck+nsize
   end if
   if (allocated(output%kg_k_gather_sym)) then
     ABI_DEALLOCATE(output%kg_k_gather_sym)
   end if
   if (size1_kg_k_gather_sym*size2_kg_k_gather_sym>0) then
     nsize=size1_kg_k_gather_sym*size2_kg_k_gather_sym
     sz1=size1_kg_k_gather_sym;sz2=size2_kg_k_gather_sym
     ABI_ALLOCATE(output%kg_k_gather_sym,(sz1,sz2))
     output%kg_k_gather_sym(:,:)=reshape(buffer_int(ipck+1:ipck+nsize),(/sz1,sz2/))
     ipck=ipck+nsize
   end if
   if (allocated(output%rdispls_sym)) then
     ABI_DEALLOCATE(output%rdispls_sym)
   end if
   if (size_rdispls_sym>0) then
     ABI_ALLOCATE(output%rdispls_sym,(size_rdispls_sym))
     output%rdispls_sym(:)=buffer_int(ipck+1:ipck+size_rdispls_sym)
     ipck=ipck+size_rdispls_sym
   end if
   if (allocated(output%sdispls_sym)) then
     ABI_DEALLOCATE(output%sdispls_sym)
   end if
   if (size_sdispls_sym>0) then
     ABI_ALLOCATE(output%sdispls_sym,(size_sdispls_sym))
     output%sdispls_sym(:)=buffer_int(ipck+1:ipck+size_sdispls_sym)
     ipck=ipck+size_sdispls_sym
   end if
   if (allocated(output%recvcounts_sym)) then
     ABI_DEALLOCATE(output%recvcounts_sym)
   end if
   if (size_recvcounts_sym>0) then
     ABI_ALLOCATE(output%recvcounts_sym,(size_recvcounts_sym))
     output%recvcounts_sym(:)=buffer_int(ipck+1:ipck+size_recvcounts_sym)
     ipck=ipck+size_recvcounts_sym
   end if
   if (allocated(output%recvcounts_sym_tot)) then
     ABI_DEALLOCATE(output%recvcounts_sym_tot)
   end if
   if (size_recvcounts_sym_tot>0) then
     ABI_ALLOCATE(output%recvcounts_sym_tot,(size_recvcounts_sym_tot))
     output%recvcounts_sym_tot(:)=buffer_int(ipck+1:ipck+size_recvcounts_sym_tot)
     ipck=ipck+size_recvcounts_sym_tot
   end if
   if (allocated(output%sendcounts_sym)) then
     ABI_DEALLOCATE(output%sendcounts_sym)
   end if
   if (size_sendcounts_sym>0) then
     ABI_ALLOCATE(output%sendcounts_sym,(size_sendcounts_sym))
     output%sendcounts_sym(:)=buffer_int(ipck+1:ipck+size_sendcounts_sym)
     ipck=ipck+size_sendcounts_sym
   end if
   if (allocated(output%sendcounts_sym_all)) then
     ABI_DEALLOCATE(output%sendcounts_sym_all)
   end if
   if (size_sendcounts_sym_all>0) then
     ABI_ALLOCATE(output%sendcounts_sym_all,(size_sendcounts_sym_all))
     output%sendcounts_sym_all(:)=buffer_int(ipck+1:ipck+size_sendcounts_sym_all)
     ipck=ipck+size_sendcounts_sym_all
   end if
   if (allocated(output%tab_proc)) then
     ABI_DEALLOCATE(output%tab_proc)
   end if
   if (size_tab_proc>0) then
     ABI_ALLOCATE(output%tab_proc,(size_tab_proc))
     output%tab_proc(:)=buffer_int(ipck+1:ipck+size_tab_proc)
     ipck=ipck+size_tab_proc
   end if
   if (allocated(output%indices_pw_fft)) then
     ABI_DEALLOCATE(output%indices_pw_fft)
   end if
   if (size_indices_pw_fft>0) then
     ABI_ALLOCATE(output%indices_pw_fft,(size_indices_pw_fft))
     output%indices_pw_fft(:)=buffer_int(ipck+1:ipck+size_indices_pw_fft)
     ipck=ipck+size_indices_pw_fft
   end if
   if (allocated(output%sendcount_fft)) then
     ABI_DEALLOCATE(output%sendcount_fft)
   end if
   if (size_sendcount_fft>0) then
     ABI_ALLOCATE(output%sendcount_fft,(size_sendcount_fft))
     output%sendcount_fft(:)=buffer_int(ipck+1:ipck+size_sendcount_fft)
     ipck=ipck+size_sendcount_fft
   end if
   if (allocated(output%senddisp_fft)) then
     ABI_DEALLOCATE(output%senddisp_fft)
   end if
   if (size_senddisp_fft>0) then
     ABI_ALLOCATE(output%senddisp_fft,(size_senddisp_fft))
     output%senddisp_fft(:)=buffer_int(ipck+1:ipck+size_senddisp_fft)
     ipck=ipck+size_senddisp_fft
   end if
   if (allocated(output%recvcount_fft)) then
     ABI_DEALLOCATE(output%recvcount_fft)
   end if
   if (size_recvcount_fft>0) then
     ABI_ALLOCATE(output%recvcount_fft,(size_recvcount_fft))
     output%recvcount_fft(:)=buffer_int(ipck+1:ipck+size_recvcount_fft)
     ipck=ipck+size_recvcount_fft
   end if
   if (allocated(output%recvdisp_fft)) then
     ABI_DEALLOCATE(output%recvdisp_fft)
   end if
   if (size_recvdisp_fft>0) then
     ABI_ALLOCATE(output%recvdisp_fft,(size_recvdisp_fft))
     output%recvdisp_fft(:)=buffer_int(ipck+1:ipck+size_recvdisp_fft)
     ipck=ipck+size_recvdisp_fft
   end if
   if (allocated(output%kg_k_fft)) then
     ABI_DEALLOCATE(output%kg_k_fft)
   end if
   if (size1_kg_k_fft*size2_kg_k_fft>0) then
     nsize=size1_kg_k_fft*size2_kg_k_fft
     sz1=size1_kg_k_fft;sz2=size2_kg_k_fft
     ABI_ALLOCATE(output%kg_k_fft,(sz1,sz2))
     output%kg_k_fft(:,:)=reshape(buffer_int(ipck+1:ipck+nsize),(/sz1,sz2/))
     ipck=ipck+nsize
   end if
   ABI_DEALLOCATE(buffer_int)
 end if

!=== Receive and unpack real data ====
 if (size_dp>0) then
   ipck=0
   ABI_ALLOCATE(buffer_dp,(size_dp))
   call xmpi_recv(buffer_dp,sender,3*tag,spaceComm,ierr)
   if (allocated(output%ffnl_gather)) then
     ABI_DEALLOCATE(output%ffnl_gather)
   end if
   if (size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather>0) then
     nsize=size1_ffnl_gather*size2_ffnl_gather*size3_ffnl_gather*size4_ffnl_gather
     sz1=size1_ffnl_gather;sz2=size2_ffnl_gather;sz3=size3_ffnl_gather;sz4=size4_ffnl_gather
     ABI_ALLOCATE(output%ffnl_gather,(sz1,sz2,sz3,sz4))
     output%ffnl_gather(:,:,:,:)=reshape(buffer_dp(ipck+1:ipck+nsize),(/sz1,sz2,sz3,sz4/))
     ipck=ipck+nsize
   end if
   if (allocated(output%kinpw_gather)) then
     ABI_DEALLOCATE(output%kinpw_gather)
   end if
   if (size_kinpw_gather>0) then
     ABI_ALLOCATE(output%kinpw_gather,(size_kinpw_gather))
     output%kinpw_gather(:)=buffer_dp(ipck+1:ipck+size_kinpw_gather)
     ipck=ipck+size_kinpw_gather
   end if
   if (allocated(output%ph3d_gather)) then
     ABI_DEALLOCATE(output%ph3d_gather)
   end if
   if (size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather>0) then
     nsize=size1_ph3d_gather*size2_ph3d_gather*size3_ph3d_gather
     sz1=size1_ph3d_gather;sz2=size2_ph3d_gather;sz3=size3_ph3d_gather
     ABI_ALLOCATE(output%ph3d_gather,(sz1,sz2,sz3))
     output%ph3d_gather(:,:,:)=reshape(buffer_dp(ipck+1:ipck+nsize),(/sz1,sz2,sz3/))
     ipck=ipck+nsize
   end if
   if (allocated(output%kpg_k_gather)) then
     ABI_DEALLOCATE(output%kpg_k_gather)
   end if
   if (size1_kpg_k_gather*size2_kpg_k_gather>0) then
     nsize=size1_kpg_k_gather*size2_kpg_k_gather
     sz1=size1_kpg_k_gather;sz2=size2_kpg_k_gather
     ABI_ALLOCATE(output%kpg_k_gather,(sz1,sz2))
     output%kpg_k_gather(:,:)=reshape(buffer_dp(ipck+1:ipck+nsize),(/sz1,sz2/))
     ipck=ipck+nsize
   end if
   ABI_DEALLOCATE(buffer_dp)
 end if

end subroutine bandfft_kpt_mpi_recv
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_savetabs
!! NAME
!!  bandfft_kpt_savetabs
!!
!! FUNCTION
!!  Save some (kpt-dependent) tabs from a bandfft_kpt datastructure
!!
!! INPUTS
!!  bandfft_kpt_in=<type(bandfft_kpt)>=bandfft_kpt datastructure
!!
!! OUTPUT
!!  ffnl(:,:,:,:)=nonlocal form factors on basis sphere
!!  ph3d(:,:,:)=3-dim structure factors, for each atom and plane wave
!!  kpg(:,:)=k+G vector for a given k point
!!  kinpw(:)=kinetic energy for each plane wave (Hartree)
!!
!! PARENTS
!!      energy,fock2ACE,forstrnps
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_savetabs(bandfft_kpt_in,ffnl,ph3d,kpg,kinpw)

!Arguments -------------------------------
 type(bandfft_kpt_type),intent(inout) :: bandfft_kpt_in
 real(dp),intent(inout),allocatable,optional :: ffnl(:,:,:,:),ph3d(:,:,:),kpg(:,:),kinpw(:)
!Local variables-------------------------------
 integer :: is1,is2,is3,is4

! *********************************************************************

 if (present(ffnl)) then
   if (allocated(ffnl)) then
     ABI_DEALLOCATE(ffnl)
   end if
   if (allocated(bandfft_kpt_in%ffnl_gather)) then
     is1=size(bandfft_kpt_in%ffnl_gather,1)
     is2=size(bandfft_kpt_in%ffnl_gather,2)
     is3=size(bandfft_kpt_in%ffnl_gather,3)
     is4=size(bandfft_kpt_in%ffnl_gather,4)
     ABI_ALLOCATE(ffnl,(is1,is2,is3,is4))
     ffnl(:,:,:,:)=bandfft_kpt_in%ffnl_gather(:,:,:,:)
   end if
 end if
 if (present(ph3d)) then
   if (allocated(ph3d)) then
     ABI_DEALLOCATE(ph3d)
   end if
   if (allocated(bandfft_kpt_in%ph3d_gather)) then
     is1=size(bandfft_kpt_in%ph3d_gather,1)
     is2=size(bandfft_kpt_in%ph3d_gather,2)
     is3=size(bandfft_kpt_in%ph3d_gather,3)
     ABI_ALLOCATE(ph3d,(is1,is2,is3))
     ph3d(:,:,:)=bandfft_kpt_in%ph3d_gather(:,:,:)
   end if
 end if
 if (present(kpg)) then
   if (allocated(kpg)) then
     ABI_DEALLOCATE(kpg)
   end if
   if (allocated(bandfft_kpt_in%kpg_k_gather)) then
     is1=size(bandfft_kpt_in%kpg_k_gather,1)
     is2=size(bandfft_kpt_in%kpg_k_gather,2)
     ABI_ALLOCATE(kpg,(is1,is2))
     kpg(:,:)=bandfft_kpt_in%kpg_k_gather(:,:)
   end if
 end if
 if (present(kinpw)) then
   if (allocated(kinpw)) then
     ABI_DEALLOCATE(kinpw)
   end if
   if (allocated(bandfft_kpt_in%kinpw_gather)) then
     is1=size(bandfft_kpt_in%kinpw_gather,1)
     ABI_ALLOCATE(kinpw,(is1))
     kinpw(:)=bandfft_kpt_in%kinpw_gather(:)
   end if
 end if

end subroutine bandfft_kpt_savetabs
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_restoretabs
!! NAME
!!  bandfft_kpt_restoretabs
!!
!! FUNCTION
!!  Restore some (kpt-dependent) tabs into a bandfft_kpt datastructure
!!
!! INPUT
!!
!! OUTPUT
!!  === Arrays to be eventually restored =====
!!  ffnl(:,:,:,:)=nonlocal form factors on basis sphere
!!  ph3d(:,:,:)=3-dim structure factors, for each atom and plane wave
!!  kpg(:,:)=k+G vector for a given k point
!!  kinpw(:)=kinetic energy for each plane wave (Hartree)
!!
!! SIDE EFFECTS
!!  bandfft_kpt_out=<type(bandfft_kpt)>=bandfft_kpt datastructure
!!
!! PARENTS
!!      energy,fock2ACE,forstrnps
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_restoretabs(bandfft_kpt_out,ffnl,ph3d,kpg,kinpw)

!Arguments -------------------------------
 type(bandfft_kpt_type),intent(inout) :: bandfft_kpt_out
 real(dp),intent(inout),allocatable,optional :: ffnl(:,:,:,:),ph3d(:,:,:),kpg(:,:),kinpw(:)
!Local variables-------------------------------
 integer :: is1,is2,is3,is4

! *********************************************************************

 if (present(ffnl)) then
   if (allocated(bandfft_kpt_out%ffnl_gather)) then
     ABI_DEALLOCATE(bandfft_kpt_out%ffnl_gather)
   end if
   if (allocated(ffnl)) then
     is1=size(ffnl,1)
     is2=size(ffnl,2)
     is3=size(ffnl,3)
     is4=size(ffnl,4)
     ABI_ALLOCATE(bandfft_kpt_out%ffnl_gather,(is1,is2,is3,is4))
     bandfft_kpt_out%ffnl_gather(:,:,:,:)=ffnl(:,:,:,:)
     ABI_DEALLOCATE(ffnl)
     bandfft_kpt_out%flag2_is_allocated=1
   end if
 end if
 if (present(ph3d)) then
   if (allocated(bandfft_kpt_out%ph3d_gather)) then
     ABI_DEALLOCATE(bandfft_kpt_out%ph3d_gather)
   end if
   if (allocated(ph3d)) then
     is1=size(ph3d,1)
     is2=size(ph3d,2)
     is3=size(ph3d,3)
     ABI_ALLOCATE(bandfft_kpt_out%ph3d_gather,(is1,is2,is3))
     bandfft_kpt_out%ph3d_gather(:,:,:)=ph3d(:,:,:)
     ABI_DEALLOCATE(ph3d)
     bandfft_kpt_out%flag2_is_allocated=1
   end if
 end if
 if (present(kpg)) then
   if (allocated(bandfft_kpt_out%kpg_k_gather)) then
     ABI_DEALLOCATE(bandfft_kpt_out%kpg_k_gather)
   end if
   if (allocated(kpg)) then
     is1=size(kpg,1)
     is2=size(kpg,2)
     ABI_ALLOCATE(bandfft_kpt_out%kpg_k_gather,(is1,is2))
     bandfft_kpt_out%kpg_k_gather(:,:)=kpg(:,:)
     ABI_DEALLOCATE(kpg)
     bandfft_kpt_out%flag2_is_allocated=1
   end if
 end if
 if (present(kinpw)) then
   if (allocated(bandfft_kpt_out%kinpw_gather)) then
     ABI_DEALLOCATE(bandfft_kpt_out%kinpw_gather)
   end if
   if (allocated(kinpw)) then
     is1=size(kinpw,1)
     ABI_ALLOCATE(bandfft_kpt_out%kinpw_gather,(is1))
     bandfft_kpt_out%kinpw_gather(:)=kinpw(:)
     ABI_DEALLOCATE(kinpw)
   end if
 end if

end subroutine bandfft_kpt_restoretabs
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_set_ikpt
!! NAME
!!  bandfft_kpt_set_ikpt
!!
!! FUNCTION
!!  Set the value of bandfft_kpt_current_ikpt variable,
!!    i.e. current index of bandfft_kpt loaded in memory
!!
!! INPUT
!!  ikpt=index of k-point (in the global array dtset%kpt)
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  bandfft_kpt_current_ikpt value changed
!!
!! PARENTS
!!      mkrho,prep_bandfft_tabs,vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

!Arguments -------------------------------
 integer,intent(in) :: ikpt
 type(MPI_type),intent(inout) :: mpi_enreg
!Local variables-------------------------------

! *********************************************************************

 if (ikpt>0) then
   bandfft_kpt_current_ikpt=mpi_enreg%my_kpttab(ikpt)
 else
   bandfft_kpt_current_ikpt=-1
 end if

end subroutine bandfft_kpt_set_ikpt
!!***

!----------------------------------------------------------------------

!!****f* m_bandfft_kpt/bandfft_kpt_get_ikpt
!! NAME
!!  bandfft_kpt_get_ikpt
!!
!! FUNCTION
!!  Get the value of bandfft_kpt_current_ikpt variable,
!!    i.e. current index of bandfft_kpt loaded in memory
!!
!! INPUT
!!
!! OUTPUT
!!  bandfft_kpt_get_ikpt= current index of bandfft_kpt
!!
!! PARENTS
!!      apply_invovl,chebfi,getghc,make_invovl,prep_fourwf
!!      prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

function bandfft_kpt_get_ikpt()

!Arguments -------------------------------
 integer :: bandfft_kpt_get_ikpt
!Local variables-------------------------------

! *********************************************************************

 bandfft_kpt_get_ikpt=bandfft_kpt_current_ikpt

end function bandfft_kpt_get_ikpt
!!***

!!****f* ABINIT/prep_bandfft_tabs
!! NAME
!! prep_bandfft_tabs
!!
!! FUNCTION
!! This routine transpose various tabs needed in bandfft parallelization
!!
!! INPUTS
!!  ikpt=index of the k-point
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!
!! SIDE EFFECTS
!!  bandfft_kpt tabs (defined in m_bandfft_kpt module)
!!
!! PARENTS
!!      energy,fock2ACE,forstrnps,vtorho
!!
!! CHILDREN
!!      bandfft_kpt_init2,bandfft_kpt_set_ikpt,mkkpg,timab,xmpi_allgatherv
!!
!! SOURCE


subroutine prep_bandfft_tabs(gs_hamk,ikpt,mkmem,mpi_enreg)

!Arguments -------------------------------
 integer,intent(in) :: ikpt,mkmem
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: dimffnl,ierr,ikpt_this_proc,ipw,lmnmax,matblk,ndatarecv,nkpg,npw_k,ntypat,spaceComm
 logical :: tabs_allocated
 real(dp) :: tsec(2)
 character(len=500)   :: message
 integer, allocatable :: recvcounts(:),rdispls(:)
 integer, allocatable :: recvcountsloc(:),rdisplsloc(:)
 real(dp),allocatable :: ffnl_gather(:,:,:,:),ffnl_little(:,:,:,:),ffnl_little_gather(:,:,:,:)
 real(dp),allocatable :: kinpw_gather(:),kpg_k_gather(:,:)
 real(dp),allocatable :: ph3d_gather(:,:,:),ph3d_little(:,:,:),ph3d_little_gather(:,:,:)

! *********************************************************************

 call timab(575,1,tsec)

 ikpt_this_proc=mpi_enreg%my_kpttab(ikpt)
 tabs_allocated=((bandfft_kpt(ikpt_this_proc)%flag1_is_allocated==1).and.&
& (ikpt_this_proc <= mkmem).and.(ikpt_this_proc/=0))

 if (.not.tabs_allocated) then
   message = ' the bandfft tabs are not allocated !'
   MSG_BUG(message)
 end if

 ntypat=gs_hamk%ntypat
 lmnmax=gs_hamk%lmnmax
 matblk=gs_hamk%matblk
 npw_k=gs_hamk%npw_k
 dimffnl=size(gs_hamk%ffnl_k,2)
 nkpg=size(gs_hamk%kpg_k,2)

 ABI_ALLOCATE(rdispls,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(recvcounts,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(rdisplsloc,(mpi_enreg%nproc_band))
 ABI_ALLOCATE(recvcountsloc,(mpi_enreg%nproc_band))

 spaceComm    =mpi_enreg%comm_band
 ndatarecv    =bandfft_kpt(ikpt_this_proc)%ndatarecv
 rdispls(:)   =bandfft_kpt(ikpt_this_proc)%rdispls(:)
 recvcounts(:)=bandfft_kpt(ikpt_this_proc)%recvcounts(:)

!---- Process FFNL
 if (associated(gs_hamk%ffnl_k)) then
   ABI_ALLOCATE(ffnl_gather,(ndatarecv,dimffnl,lmnmax,ntypat))
   ABI_ALLOCATE(ffnl_little,(dimffnl,lmnmax,ntypat,npw_k))
   ABI_ALLOCATE(ffnl_little_gather,(dimffnl,lmnmax,ntypat,ndatarecv))
   do ipw=1,npw_k
     ffnl_little(:,:,:,ipw)=gs_hamk%ffnl_k(ipw,:,:,:)
   end do
   recvcountsloc(:)=recvcounts(:)*dimffnl*lmnmax*ntypat
   rdisplsloc(:)=rdispls(:)*dimffnl*lmnmax*ntypat
   call xmpi_allgatherv(ffnl_little,npw_k*dimffnl*lmnmax*ntypat,ffnl_little_gather,&
&   recvcountsloc(:),rdisplsloc,spaceComm,ierr)
   do ipw=1,ndatarecv
     ffnl_gather(ipw,:,:,:)=ffnl_little_gather(:,:,:,ipw)
   end do
   ABI_DEALLOCATE(ffnl_little)
   ABI_DEALLOCATE(ffnl_little_gather)
 else
   ABI_ALLOCATE(ffnl_gather,(0,0,0,0))
 end if

!---- Process PH3D
 if (associated(gs_hamk%ph3d_k)) then
   ABI_ALLOCATE(ph3d_gather,(2,ndatarecv,matblk))
   ABI_ALLOCATE(ph3d_little,(2,matblk,npw_k))
   ABI_ALLOCATE(ph3d_little_gather,(2,matblk,ndatarecv))
   recvcountsloc(:)=recvcounts(:)*2*matblk
   rdisplsloc(:)=rdispls(:)*2*matblk
   do ipw=1,npw_k
     ph3d_little(:,:,ipw)=gs_hamk%ph3d_k(:,ipw,:)
   end do
   call xmpi_allgatherv(ph3d_little,npw_k*2*matblk,ph3d_little_gather,&
&   recvcountsloc(:),rdisplsloc,spaceComm,ierr)
   do ipw=1,ndatarecv
     ph3d_gather(:,ipw,:)=ph3d_little_gather(:,:,ipw)
   end do
   ABI_DEALLOCATE(ph3d_little_gather)
   ABI_DEALLOCATE(ph3d_little)
 else
   ABI_ALLOCATE(ph3d_gather,(0,0,0))
 end if

!---- Process KPG_K
 if (associated(gs_hamk%kpg_k)) then
   ABI_ALLOCATE(kpg_k_gather,(ndatarecv,nkpg))
   if (nkpg>0) then
     call mkkpg(bandfft_kpt(ikpt_this_proc)%kg_k_gather,kpg_k_gather,gs_hamk%kpt_k,nkpg,ndatarecv)
   end if
 else
   ABI_ALLOCATE(kpg_k_gather,(0,0))
 end if

!---- Process KINPW
 if (associated(gs_hamk%kinpw_k)) then
   ABI_ALLOCATE(kinpw_gather,(ndatarecv))
   recvcountsloc(:)=recvcounts(:)
   rdisplsloc(:)=rdispls(:)
   call xmpi_allgatherv(gs_hamk%kinpw_k,npw_k,kinpw_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
 else
   ABI_ALLOCATE(kinpw_gather,(0))
 end if

 ABI_DEALLOCATE(recvcounts)
 ABI_DEALLOCATE(rdispls)
 ABI_DEALLOCATE(recvcountsloc)
 ABI_DEALLOCATE(rdisplsloc)

 call bandfft_kpt_init2(bandfft_kpt,dimffnl,ffnl_gather,ikpt_this_proc,kinpw_gather,kpg_k_gather,&
& lmnmax,matblk,mkmem,ndatarecv,nkpg,ntypat,ph3d_gather)

 ABI_DEALLOCATE(ffnl_gather)
 ABI_DEALLOCATE(ph3d_gather)
 ABI_DEALLOCATE(kinpw_gather)
 ABI_DEALLOCATE(kpg_k_gather)

!---- Store current kpt index
 call bandfft_kpt_set_ikpt(ikpt,mpi_enreg)

 call timab(575,2,tsec)

end subroutine prep_bandfft_tabs
!!***

END MODULE m_bandfft_kpt
!!***
