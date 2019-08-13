!{\src2tex{textfont=tt}}
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
!! * datafiles_type: gather all the variables related to files
!! * dataset_type: the "dataset" for the main abinit code
!! * MPI_type: the data related to MPI parallelization
!! * macro_uj_type: TO BE COMPLETED
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

!----------------------------------------------------------------------

!!****t* defs_datatypes/datafiles_type
!! NAME
!! datafiles_type
!!
!! FUNCTION
!! The datafiles_type structures datatype gather all the variables
!! related to files, such as filename, and file units.
!! For one dataset, it is initialized in 95_drive/dtfil_init1.F90,
!! and will not change at all during the treatment of the dataset.
!!
!! SOURCE

 type datafiles_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! These keywords are only used in algorithms using images of the cell
  integer :: getwfk_from_image
   ! index of image from which read WFK file (0 if standard WFK)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getden_from_image
   ! index of image from which read DEN file (0 if standard DEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getpawden_from_image
   ! index of image from which read PAWDEN file (0 if standard PAWDEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: ireadddb
   ! ireadddb non-zero  if the ddb file must be read

  integer :: ireadden
   ! ireadden non-zero  if the den file must be read

  integer :: ireadkden
   ! ireadkden non-zero  if the kden file must be read

  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffk file must be read
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1 file must be read
   !   (if ird1wf non-zero or get1wf non-zero)

  integer :: unchi0  ! unit number for chi0 files
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: unddk   ! unit number for ddk 1WF file
  integer :: undkdk  ! unit number for 2WF file (dkdk)
  integer :: undkde  ! unit number for 2WF file (dkde)
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unkss   ! unit number for KSS file
  integer :: unqps   ! unit number for QPS file
  integer :: unscr   ! unit number for SCR file
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwff3  ! unit number for wavefunctions, number three
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwft3  ! unit number for wavefunctions, temporary three
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: unpaw   ! unit number for temporary PAW data (for ex. rhoij_nk) (Paw only)
  integer :: unpaw1  ! unit number for temporary PAW first-order cprj1=<c1_k,q|p>(1) data
  integer :: unpawq  ! unit number for temporary PAW cprjq=<c+_k+q|p> at k+qdata
  integer :: unpos   ! unit number for restart molecular dynamics

  ! TODO: All this strings should be initialized with ABI_NOFILE
  ! so that we can easily test for path /= ABI_NOFILE instead of getwfk /= 0 or irdwfk /= 0

  character(len=fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input:
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character(len=fnlen) :: filddbsin
   ! if no dataset mode             : abi//'DDB'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DDB'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DDB'

  character(len=fnlen) :: fildensin
   ! if no dataset mode             : abi//'DEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character(len=fnlen) :: fildvdbin
   ! if no dataset mode             : abi//'DVDB'
   ! if dataset mode, and getdvdb==0 : abi//'_DS'//trim(jdtset)//'DVDB'
   ! if dataset mode, and getdvdb/=0 : abo//'_DS'//trim(jgetden)//'DVDB'

  character(len=fnlen) :: filpotin
   ! Filename used to read POT file.
   ! Initialize via getpot_path

  character(len=fnlen) :: filkdensin
   ! if no dataset mode             : abi//'KDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'KDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'KDEN'

  character(len=fnlen) :: filpawdensin
   ! if no dataset mode             : abi//'PAWDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'PAWDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'PAWDEN'

! character(len=fnlen) :: filpsp(ntypat)
   ! the filenames of the pseudopotential files, from the standard input.

  character(len=fnlen) :: filstat
   ! tmp//'_STATUS'

  character(len=fnlen) :: fnamewffk
   ! the name of the ground-state wavefunction file to be read (see driver.F90)

  character(len=fnlen) :: fnamewffq
   ! the name of the k+q ground-state wavefunction file to be read (see driver.F90)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdelfd
   ! the generic name of the electric field response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkdk
   ! the generic name of the 2nd order dkdk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkde
   ! the generic name of the 2nd order dkde response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewff1
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fildens1in   ! to be described by MVeithen
  character(len=fnlen) :: fname_tdwf
  character(len=fnlen) :: fname_w90

  character(len=fnlen) :: fnametmp_wf1
  character(len=fnlen) :: fnametmp_wf2
  character(len=fnlen) :: fnametmp_1wf1
  character(len=fnlen) :: fnametmp_1wf2
  character(len=fnlen) :: fnametmp_wfgs
  character(len=fnlen) :: fnametmp_wfkq
   ! Set of filenames formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _WF1, _WF2 ...
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_kg
  character(len=fnlen) :: fnametmp_kgq
  character(len=fnlen) :: fnametmp_kg1
  character(len=fnlen) :: fnametmp_dum
  character(len=fnlen) :: fnametmp_ylm
  character(len=fnlen) :: fnametmp_ylm1
  character(len=fnlen) :: fnametmp_paw
  character(len=fnlen) :: fnametmp_paw1
  character(len=fnlen) :: fnametmp_pawq
   ! Set of filenames formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _KG, _DUM, followed
   ! by the index of the processor.
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_cg
  character(len=fnlen) :: fnametmp_cprj
  character(len=fnlen) :: fnametmp_eig
  character(len=fnlen) :: fnametmp_1wf1_eig
  character(len=fnlen) :: fnametmp_fft
  character(len=fnlen) :: fnametmp_kgs
  character(len=fnlen) :: fnametmp_sustr
  character(len=fnlen) :: fnametmp_tdexcit
  character(len=fnlen) :: fnametmp_tdwf

!@Bethe-Salpeter
! New files introduced for the Bethe-Salpeter part.

   character(len=fnlen) :: fnameabi_bsham_reso
    ! if no dataset mode             : abi//'BSR'
    ! if dataset mode, and getbsreso==0 : abi//'_DS'//trim(jdtset)//'BSR'
    ! if dataset mode, and getbsreso/=0 : abo//'_DS'//trim(jget_reso_bsham)//'BSR'

   character(len=fnlen) :: fnameabi_bsham_coup
    ! if no dataset mode             : abi//'BSC'
    ! if dataset mode, and getbscoup==0 : abi//'_DS'//trim(jdtset)//'BSC'
    ! if dataset mode, and getbscoup/=0 : abo//'_DS'//trim(jget_coup_bsham)//'BSC'

  character(len=fnlen) :: fnameabi_bseig
   ! The name of the file containing the eigenstates and eigenvalues of the Bethe-Salpeter Hamiltonian
   ! if no dataset mode             : abi//'BS_EIG'
   ! if dataset mode, and getbseig==0 : abi//'_DS'//trim(jdtset)//'BS_EIG'
   ! if dataset mode, and getbseig/=0 : abo//'_DS'//trim(jget_bseig)//'BS_EIG'

   character(len=fnlen) :: fnameabi_haydock
   ! The prefix used to construct the names of the files containing the coefficients of the
   ! continued fractions produced by the Haydock iterative algorithm.
   ! if no dataset mode             : abi//'HAYDOCK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'HAYDOCK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'HAYDOCK'

   character(len=fnlen) :: fnameabi_wfkfine
   ! The name of the file containing the wavefunctions on a fine grid
   ! if no dataset mode             : abi//'WFK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'WFK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'WFK'

!END @BEthe-Salpeter

!The following filenames do not depend on itimimage, iimage and itime loops.
!Note the following convention:
!  fnameabo_* are filenames used for ouput results (writing)
!  fnameabi_* are filenames used for data that should be read by the code.
!  fnametmp_* are filenames used for temporary files that should be erased at the end of each dataset.
!
!If a file does not have the corresponding "abi" or the corresponding "abo" name, that means that
!that particular file is only used for writing or for reading results, respectively.

  character(len=fnlen) :: fnameabi_efmas
  character(len=fnlen) :: fnameabi_hes
  character(len=fnlen) :: fnameabi_phfrq
  character(len=fnlen) :: fnameabi_phvec
  character(len=fnlen) :: fnameabi_qps
  character(len=fnlen) :: fnameabi_scr            ! SCReening file (symmetrized inverse dielectric matrix)
  character(len=fnlen) :: fnameabi_sus            ! KS independent-particle polarizability file
  character(len=fnlen) :: fnameabo_ddb
  character(len=fnlen) :: fnameabo_den
  character(len=fnlen) :: fnameabo_dos
  character(len=fnlen) :: fnameabo_dvdb
  character(len=fnlen) :: fnameabo_eelf
  character(len=fnlen) :: fnameabo_eig
  character(len=fnlen) :: fnameabo_eigi2d
  character(len=fnlen) :: fnameabo_eigr2d
  character(len=fnlen) :: fnameabo_em1
  character(len=fnlen) :: fnameabo_em1_lf
  character(len=fnlen) :: fnameabo_em1_nlf
  character(len=fnlen) :: fnameabo_fan
  character(len=fnlen) :: fnameabo_gkk
  character(len=fnlen) :: fnameabo_gw
  character(len=fnlen) :: fnameabo_gw_nlf_mdf
  character(len=fnlen) :: fnameabo_kss
  character(len=fnlen) :: fnameabo_moldyn
  character(len=fnlen) :: fnameabo_pot
  character(len=fnlen) :: fnameabo_qps            ! Quasi-Particle band structure file.
  character(len=fnlen) :: fnameabo_qp_den
  character(len=fnlen) :: fnameabo_qp_pawden      ! Full QP density
  character(len=fnlen) :: fnameabo_qp_dos
  character(len=fnlen) :: fnameabo_qp_eig
  character(len=fnlen) :: fnameabo_rpa
  character(len=fnlen) :: fnameabo_rpa_nlf_mdf
  character(len=fnlen) :: fnameabo_scr
  character(len=fnlen) :: fnameabo_sgm
  character(len=fnlen) :: fnameabo_sgr
  character(len=fnlen) :: fnameabo_sig
  character(len=fnlen) :: fnameabo_spcur
  character(len=fnlen) :: fnameabo_sus
  character(len=fnlen) :: fnameabo_vha
  character(len=fnlen) :: fnameabo_vpsp
  character(len=fnlen) :: fnameabo_vso
  character(len=fnlen) :: fnameabo_vxc
  character(len=fnlen) :: fnameabo_wan
  character(len=fnlen) :: fnameabo_wfk
  character(len=fnlen) :: fnameabo_wfq
  character(len=fnlen) :: fnameabo_w90
  character(len=fnlen) :: fnameabo_1wf
  character(len=fnlen) :: fnameabo_gwdiag
  character(len=fnlen) :: fnameabo_nlcc_derivs
  character(len=fnlen) :: fnameabo_pspdata

!The following filenames are initialized only iniside itimimage, iimage and itime loops,
!and are appended with the adequate specifier 'app'.

  character(len=fnlen) :: fnameabo_app
  character(len=fnlen) :: fnameabo_app_atmden_core
  character(len=fnlen) :: fnameabo_app_atmden_full
  character(len=fnlen) :: fnameabo_app_atmden_val
  character(len=fnlen) :: fnameabo_app_n_tilde
  character(len=fnlen) :: fnameabo_app_n_one
  character(len=fnlen) :: fnameabo_app_nt_one
  character(len=fnlen) :: fnameabo_app_bxsf
  character(len=fnlen) :: fnameabo_app_cif
  character(len=fnlen) :: fnameabo_app_den
  character(len=fnlen) :: fnameabo_app_dos
  character(len=fnlen) :: fnameabo_app_elf
  character(len=fnlen) :: fnameabo_app_elf_down
  character(len=fnlen) :: fnameabo_app_elf_up
  character(len=fnlen) :: fnameabo_app_eig
  character(len=fnlen) :: fnameabo_app_fatbands
  character(len=fnlen) :: fnameabo_app_gden1
  character(len=fnlen) :: fnameabo_app_gden2
  character(len=fnlen) :: fnameabo_app_gden3
  character(len=fnlen) :: fnameabo_app_geo
  character(len=fnlen) :: fnameabo_app_kden
  character(len=fnlen) :: fnameabo_app_lden
  character(len=fnlen) :: fnameabo_app_nesting
  character(len=fnlen) :: fnameabo_app_pawden
  character(len=fnlen) :: fnameabo_app_pot
  character(len=fnlen) :: fnameabo_app_opt
  character(len=fnlen) :: fnameabo_app_opt2
  character(len=fnlen) :: fnameabo_app_stm
  character(len=fnlen) :: fnameabo_app_vclmb
  character(len=fnlen) :: fnameabo_app_vha
  character(len=fnlen) :: fnameabo_app_vhxc
  character(len=fnlen) :: fnameabo_app_vhpsp
  character(len=fnlen) :: fnameabo_app_vpsp
  character(len=fnlen) :: fnameabo_app_vxc
  character(len=fnlen) :: fnameabo_app_wfk
  character(len=fnlen) :: fnameabo_app_1dm
  character(len=fnlen) :: fnameabo_app_vha_1dm
  character(len=fnlen) :: fnameabo_app_vclmb_1dm
  character(len=fnlen) :: fnametmp_app_den
  character(len=fnlen) :: fnametmp_app_kden

 end type datafiles_type
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/ab_dimensions
!! NAME
!! ab_dimensions
!!
!! FUNCTION
!! One record for each dimension of arrays used in ABINIT.
!! Will be used to e.g.:
!! - contain the maximum size attained over all datasets (mxvals)
!! - indicate whether this dimension is the same for all datasets or not (multivals).
!! Used for example inside outvars
!!
!! SOURCE

 type ab_dimensions

    integer :: ga_n_rules   ! maximal value of input ga_n_rules for all the datasets
    integer :: gw_nqlwl     ! maximal value of input gw_nqlwl for all the datasets
    integer :: lpawu        ! maximal value of input lpawu for all the datasets
    integer :: mband
    integer :: mband_upper ! maximal value of input nband for all the datasets
                           ! Maybe this one could be removed
    integer :: natom
    integer :: natpawu     ! maximal value of number of atoms on which +U is applied for all the datasets
    integer :: natsph      ! maximal value of input natsph for all the datasets
    integer :: natsph_extra  ! maximal value of input natsph_extra for all the datasets
    integer :: natvshift   ! maximal value of input natvshift for all the datasets
    integer :: nberry = 20 ! This is presently a fixed value. Should be changed.
    integer :: nbandhf
    integer :: nconeq      ! maximal value of input nconeq for all the datasets
    integer :: n_efmas_dirs
    integer :: nfreqsp
    integer :: n_projection_frequencies
    integer :: nimage
    integer :: nimfrqs
    integer :: nkpt       ! maximal value of input nkpt for all the datasets
    integer :: nkptgw     ! maximal value of input nkptgw for all the datasets
    integer :: nkpthf     ! maximal value of input nkpthf for all the datasets
    integer :: nnos       ! maximal value of input nnos for all the datasets
    integer :: nqptdm     ! maximal value of input nqptdm for all the datasets
    integer :: nshiftk
    integer :: nsp
    integer :: nspinor    ! maximal value of input nspinor for all the datasets
    integer :: nsppol     ! maximal value of input nsppol for all the datasets
    integer :: nsym       ! maximum number of symmetries
    integer :: ntypalch
    integer :: ntypat     ! maximum number of types of atoms
    integer :: nzchempot  ! maximal value of input nzchempot for all the datasets

 end type ab_dimensions
!!***

!----------------------------------------------------------------------

!!****t* defs_abitypes/macro_uj_type
!! NAME
!! dtmacro_uj
!!
!! FUNCTION
!! This data type contains the potential shifts and the occupations
!! for the determination of U and J for the DFT+U calculations.
!! iuj=1,2: non-selfconsistent calculations. iuj=3,4 selfconsistent calculations.
!! iuj=2,4  => pawujsh<0 ; iuj=1,3 => pawujsh >0,
!!
!! SOURCE

 type macro_uj_type

! Integer
  integer :: iuj        ! dataset treated
  integer :: nat        ! number of atoms U (J) is determined on
  integer :: ndtset     ! total number of datasets
  integer :: nspden     ! number of densities treated
  integer :: macro_uj   ! which mode the determination runs in
  integer :: pawujat    ! which atom U (J) is determined on
  integer :: pawprtvol  ! controlling amount of output
  integer :: option     ! controls the determination of U (1 with compensating charge bath, 2 without)
  integer :: dmatpuopt  ! controls the renormalisation of the PAW projectors

! Real
  real(dp) :: diemix    ! mixing parameter
  real(dp) :: mdist     ! maximal distance of ions
  real(dp) :: pawujga   ! gamma for inversion of singular matrices
  real(dp) :: ph0phiint ! integral of phi(:,1)*phi(:,1)
  real(dp) :: pawujrad  ! radius to which U should be extrapolated.
  real(dp) :: pawrad    ! radius of the paw atomic data

! Integer arrays
  integer , allocatable  :: scdim(:)
  ! size of supercell

! Real arrays
  real(dp) , allocatable :: occ(:,:)
  ! occupancies after a potential shift: occ(ispden,nat)

  real(dp) , allocatable :: rprimd(:,:)
  ! unit cell for symmetrization

  real(dp) , allocatable :: vsh(:,:)
  ! potential shifts on atoms, dimensions: nspden,nat

  real(dp) , allocatable :: xred(:,:)
  ! atomic position for symmetrization

  real(dp) , allocatable :: wfchr(:)
  ! wfchr(1:3): zion, n and l of atom on which projection is done
  ! wfchr(4:6): coefficients ai of a0+a1*r+a2*r^2, fit to the wfc for r< r_paw

  real(dp), allocatable :: zioneff(:)
  ! zioneff(ij_proj), "Effective charge"*n "seen" at r_paw, deduced from Phi at r_paw, n:
  ! pricipal quantum number; good approximation to model wave function outside PAW-sphere through

 end type macro_uj_type
!!***

end module defs_abitypes
!!***
