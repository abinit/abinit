!!****m* ABINIT/defs_abitypes
!! NAME
!! defs_abitypes
!!
!! FUNCTION
!! This module contains definitions of high-level structured datatypes for the ABINIT package.
!!
!! If you are sure a new high-level structured datatype is needed,
!! write it here, and DOCUMENT it properly (not all datastructure here are well documented, it is a shame ...).
!! Do not forget: you will likely be the major winner if you document properly.
!!
!! Proper documentation of a structured datatype means:
!!  (1) Mention it in the list just below
!!  (2) Describe it in the NOTES section
!!  (3) Put it in alphabetical order in the the main section of this module
!!  (4) Document each of its records, except if they are described elsewhere
!!      (this exception is typically the case of the dataset associated with
!!      input variables, for which there is a help file)
!!  (5) Declare variables on separated lines in order to reduce the occurence of git conflicts.
!!
!! List of datatypes:
!! * MPI_type: the data related to MPI parallelization
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_abitypes

 use defs_basis
 use m_abicore
 use m_distribfft

 implicit none
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/MPI_type
!! NAME
!! MPI_type
!!
!! FUNCTION
!! The MPI_type structured datatype gather different information
!! about the MPI parallelisation: number of processors,
!! the index of my processor, the different groups of processors, etc ...
!!
!! SOURCE

 type MPI_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.
! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.

! *****************************************************************************************
! Please make sure that initmpi_seq is changed so that any variable or any flag in MPI_type
! is initialized with the value used for sequential executions.
! In particular any MPI communicator should be set to MPI_COMM_SELF
! ************************************************************************************

  ! Set of variables for parallelism, that do NOT depend on input variables.
  ! These are defined for each dataset

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Main variables for parallelisation
  integer :: comm_world
  ! number of the world communicator MPI COMM WORLD

  integer :: me
  ! number of my processor in the group of all processors

  integer :: nproc
  ! number of processors

  integer :: me_g0
  ! if set to 1, means that the current processor is taking care of the G(0 0 0) planewave.

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over atoms (PAW)
   integer :: comm_atom
   ! Communicator over atoms

   integer :: nproc_atom
   ! Size of the communicator over atoms

   integer :: my_natom
   ! Number of atoms treated by current proc

   integer,pointer :: my_atmtab(:) => null()
   ! Indexes of the atoms treated by current processor

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over perturbations
   integer :: paral_pert
   ! to activate parallelisation over perturbations for linear response

   integer :: comm_pert
   ! communicator for calculating perturbations

   integer :: comm_cell_pert
   ! general communicator over all processors treating the same cell

   integer :: me_pert
   ! number of my processor in my group of perturbations

   integer :: nproc_pert
   ! number of processors in my group of perturbations

   integer, allocatable :: distrb_pert(:)
   ! distrb_pert(1:npert)
   ! index of processor treating each perturbation

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over images
   integer :: paral_img
   ! Flag activated if parallelization over image is on

   integer :: my_nimage
   ! Number of images of the cell treated by current proc (i.e. local nimage)

   integer :: comm_img
   ! Communicator over all images

   integer :: me_img
   ! Index of my processor in the comm. over all images

   integer :: nproc_img
   ! Size of the communicator over all images

   integer,allocatable :: distrb_img(:)
   ! distrb_img(1:dtset%nimage)
   ! index of processor treating each image (in comm_img communicator)

   integer,allocatable :: my_imgtab(:)
   ! index_img(1:my_nimage) indexes of images treated by current proc

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over the cell
   integer :: comm_cell
   ! local Communicator over all processors treating the same cell

   integer :: me_cell
   ! Index of my processor in the comm. over one cell

   integer :: nproc_cell
   ! Size of the communicator over one cell

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over fft
   integer :: comm_fft
   ! Communicator over fft

   integer :: me_fft
   ! Rank of my processor in my group of FFT

   integer :: nproc_fft
   ! number of processors in my group of FFT

   type(distribfft_type),pointer :: distribfft  => null()
   ! Contains all the information related to the FFT distribution

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over bands
   integer :: paralbd
    ! paralbd=0 : (no //ization on bands)
    ! paralbd=1 : (//ization on bands)

   integer :: comm_band
   ! Communicator over bands

   integer :: me_band
   ! Rank of my proc in my group of bands

   integer :: nproc_band
   ! Number of procs on which we distribute bands

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the spinor parallelisation
   integer :: paral_spinor
   ! Flag: activation of parallelization over spinors

   integer :: comm_spinor
   ! Communicator over spinors

   integer :: me_spinor
   ! Rank of my proc in the communicator over spinors
   ! Note: me_spinor is related to the index treated by current proc
   ! (nspinor_index= mpi_enreg%me_spinor + 1)

   integer :: nproc_spinor
   ! Number of procs on which we distribute spinors

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the kpt/nsppol parallelisation
   integer :: comm_kpt
   ! Communicator over kpt

   integer :: me_kpt
   ! Rank of my proc in the communicator over kpt

   integer :: nproc_kpt
   ! Number of procs on which we distribute kpt

   integer, allocatable :: proc_distrb(:,:,:)
    ! proc_distrb(nkpt,mband,nsppol)
    ! number of the processor that will treat
    ! each band in each k point.

   integer :: my_isppoltab(2)
    ! my_isppoltab(2) contains the flags telling which value of isppol is treated by current proc
    ! in sequential, its value is (1,0) when nsppol=1 and (1,1) when nsppol=2
    ! in parallel,   its value is (1,0) when nsppol=1 and (1,0) when nsppol=2 and up-spin is treated
    !                                                  or (0,1) when nsppol=2 and dn-spin is treated
    !                                                  or (1,1) when nsppol=2 and both spins are treated

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the band-FFT-kpt-spinor parallelisation
   integer :: paral_kgb
   ! Flag: activation of parallelization over kpt/band/fft

   integer :: bandpp
   ! # of Bands Per Processor

   integer :: comm_bandspinorfft
   ! Cartesian communicator over band-fft-spinor

   integer :: comm_bandfft
   ! Cartesian communicator over the band-fft

   integer :: comm_kptband
   ! Communicator over kpt-band subspace

   integer :: comm_spinorfft
   ! Communicator over fft-spinors subspace

   integer :: comm_bandspinor
   ! Communicator over band-spinors subspace

   integer, allocatable :: my_kgtab(:,:)
    ! (mpw, mkmem)
    ! Indexes of kg treated by current proc
    ! i.e. mapping betwee the G-vector stored by this proc and the list of G-vectors
    ! one would have in the sequential version. See kpgsph in m_fftcore.

   integer, allocatable :: my_kpttab(:)
    ! Indicates the correspondence between the ikpt and ikpt_this_proc

   real(dp) :: pw_unbal_thresh
    !Threshold (in %) activating the plane-wave load balancing process (see kpgsph routine)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over kpt/nsppol in the Berry Phase case
   integer, allocatable :: kptdstrb(:,:,:)
    ! kptdstrb(me,ineigh,ikptloc)
    ! tab of processors required for dfptnl_mv.f and berryphase_new.f

   integer, allocatable :: kpt_loc2fbz_sp(:,:,:)
    ! kpt_loc2fbz_sp(nproc, dtefield%fmkmem_max, 2)
    ! K-PoinT LOCal TO Full Brilloin Zone and Spin Polarization
    ! given a processor and the local number of the k-point for this proc,
    ! give the number of the k-point in the FBZ and the isppol;
    ! necessary for synchronisation in berryphase_new
    ! kpt_loc2fbz(iproc, ikpt_loc,1) = ikpt
    ! kpt_loc2fbz(iproc, ikpt_loc,2) = isppol

   integer, allocatable :: kpt_loc2ibz_sp(:,:,:)

   ! TODO: Is is still used?
   integer, allocatable :: mkmem(:)

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for Hartree-Fock's parallelisation
   integer :: paral_hf
   ! Flag: activation of parallelization for Hartree-Fock

   integer :: comm_hf
   ! Communicator over the k-points and bands of occupied states for Hartree-Fock

   integer :: me_hf
   ! Rank of my proc in the communicator for Hartree-Fock

   integer :: nproc_hf
   ! Number of procs on which we distribute the occupied states for Hartree-Fock

   integer, allocatable :: distrb_hf(:,:,:)
    ! distrb_hf(nkpthf,nbandhf,1)
    ! index of processor treating each occupied states for Hartree Fock.
    ! No spin-dependence because only the correct spin is treated (in parallel) or both spins are considered (sequential)
    ! but we keep the third dimension (always equal to one) to be able to use the same routines as the one for proc_distrb

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the wavelet parallelisation
   integer :: comm_wvl
   ! Communicator over real space grids for WVLs

   integer :: me_wvl
   ! Rank of my proc for WVLs

   integer :: nproc_wvl
   ! Number of procs for WVLs
   ! Array to store the description of the scaterring in real space of
   ! the potentials and density. It is allocated to (0:nproc-1,4).
   ! The four values are:
   ! - the density size in z direction ( = ngfft(3)) ;
   ! - the potential size in z direction ( <= ngfft(3)) ;
   ! - the position of the first value in the complete array ;
   ! - the shift for the potential in the array.
   ! This array is a pointer on a BigDFT handled one.

   integer, pointer :: nscatterarr(:,:) => null()
   ! Array to store the total size (of this proc) of the potentails arrays when
   ! the memory is distributed following nscatterarr.
   ! This array is a pointer on a BigDFT handled one.

   integer, pointer :: ngatherarr(:,:) => null()
   ! Store the ionic potential size in z direction.
   ! This array is a pointer on a BigDFT handled one.

   integer :: ngfft3_ionic
   ! End wavelet additions

 end type MPI_type
!!***

end module defs_abitypes
!!***
