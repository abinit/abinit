!!****m* ABINIT/m_mpinfo
!! NAME
!! m_mpinfo
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MT, GG, XG, FJ, AR, MB, CMartins)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! TODO
!!  Change the name of the datatype: (MPI_|mpi_) is a reserved keyword
!!  and should not be used in client code!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_mpinfo

 use defs_basis
 use m_errors
 use m_abicore
#if defined HAVE_MPI2
 use mpi
#endif
 use m_xmpi
 use m_sort
 use m_distribfft
 use m_dtset

 use defs_abitypes,   only : MPI_type
 use m_fstrings,      only : sjoin, ltoa
 use m_io_tools,      only : file_exists, open_file
 use m_libpaw_tools,  only : libpaw_write_comm_set
 use m_paral_atom,    only : get_my_natom, get_my_atmtab

 implicit none

 private

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 public :: init_mpi_enreg        ! Initialise a mpi_enreg structure with dataset independent values.
 public :: nullify_mpi_enreg     ! nullify a mpi_enreg datastructure
 public :: destroy_mpi_enreg     ! Free memory
 public :: copy_mpi_enreg        ! Copy a mpi_enreg datastructure into another.
 public :: set_mpi_enreg_fft     ! Set the content of a MPI datastructure in order to call fourwf/fourdp
 public :: unset_mpi_enreg_fft   ! Unset the content of a MPI datastructure used to call fourwf/fourdp
 public :: ptabs_fourdp          ! Return *pointers* to the internal MPI-FFT tables used in fourdp
 public :: ptabs_fourwf          ! Return *pointers* to the internal MPI-FFT tables used in fourwf
 public :: mpi_distrib_is_ok     ! Check if a MPI datastructure contains number of processors
                                 ! compatible (in terms of efficiency) with the number of spins/kpts/bands
 public :: proc_distrb_cycle     ! Test a condition to cycle

 public :: initmpi_seq           ! Initializes the MPI information for sequential use.
 public :: initmpi_world         ! %comm_world is redifined for the number of processors on which ABINIT is launched

 public :: initmpi_atom          ! Initializes the mpi information for parallelism over atoms (PAW).
 public :: clnmpi_atom           ! Cleans-up the mpi information for the parallelism over atoms (PAW).

 public :: initmpi_grid          ! Initializes the MPI information for the grid
 public :: clnmpi_grid           ! Cleans-up the mpi information for parallelism over grid (kpt/band/fft).

 public :: initmpi_img           ! Initializes the mpi information for parallelism over images of the cell (npimage>1).
 public :: clnmpi_img            ! Cleans-up the mpi information for parallelism over images of the cell (npimage>1).

 public :: initmpi_pert          ! Creates group for Parallelization over Perturbations.
 public :: clnmpi_pert           ! Cleans-up the mpi information for parallelization over perturbations.

 public :: initmpi_band          !  Initializes the mpi information for band parallelism (paralbd=1).

! Helper functions.
 public :: pre_gather
 public :: pre_scatter
 public :: iwrite_fftdatar  ! Select the subset of processors that will write density/potential files.

 public :: distrb2          ! Creates the tabs of repartition of processors for sharing the jobs on k-points, spins and bands.
 public :: distrb2_hf       ! Ceate the tabs of repartition for Hartree-Fock calculations.
!!***

CONTAINS  !========================================================================================
!!***

!!****f* m_mpinfo/init_mpi_enreg
!! NAME
!! init_mpi_enreg
!!
!! FUNCTION
!!  Initialise a mpi_enreg structure with dataset independent values.
!!  Other values of mpi_enreg are dataset dependent, and should NOT be initialized
!!  inside abinit.F90 .
!!  XG 071118 : At present several other values are
!!  initialized temporarily inside invars1.F90, FROM THE DTSET
!!  VALUES. In order to releave the present constraint of having mpi_enreg
!!  equal for all datasets, they should be reinitialized from the dtset values
!!  inside invars2m.F90 (where there is a loop over datasets, and finally,
!!  reinitialized from the dataset values inside each big routine called by driver,
!!  according to the kind of parallelisation that is needed there.
!!  One should have one init_mpi_dtset routine (or another name) per big routine (well, there is also
!!  the problem of TDDFT ...). Also, one should have a clean_mpi_dtset called at the end
!!  of each big routine, as well as invars1.F90 or invars2m.F90 .
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=All pointer set to null().
!!
!! PARENTS
!!      lapackprof,m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_mpi_enreg(mpi_enreg)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

!Default for sequential use
 call initmpi_seq(mpi_enreg)
!Initialize MPI
#if defined HAVE_MPI
 mpi_enreg%comm_world=xmpi_world
 mpi_enreg%me = xmpi_comm_rank(xmpi_world)
 mpi_enreg%nproc = xmpi_comm_size(xmpi_world)
#endif

end subroutine init_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/nullify_mpi_enreg
!! NAME
!! nullify_mpi_enreg
!!
!! FUNCTION
!!  nullify a mpi_enreg datastructure
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=All pointer set to null().
!!
!! PARENTS
!!      m_fft_prof,m_mpinfo
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_mpi_enreg(MPI_enreg)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 nullify(mpi_enreg%nscatterarr)
 nullify(mpi_enreg%ngatherarr)
 nullify(mpi_enreg%my_atmtab)
 nullify(mpi_enreg%distribfft)

 end subroutine nullify_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/destroy_mpi_enreg
!! NAME
!! destroy_mpi_enreg
!!
!! FUNCTION
!!  Destroy a mpi_enreg datastructure
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism.
!!
!! PARENTS
!!      abinit,conducti,cut3d,fftprof,lapackprof,m_bethe_salpeter,m_cut3d,m_ddk
!!      m_dfpt_nstwf,m_dvdb,m_eph_driver,m_epjdos,m_fft,m_fft_prof,m_fftcore
!!      m_fock,m_gsphere,m_gwls_hamiltonian,m_inwffil,m_ioarr,m_ksdiago,m_kxc
!!      m_mlwfovlp_qp,m_mover_effpot,m_paw_optics,m_pawpwij,m_positron
!!      m_ppmodel,m_prcref,m_scfcv_core,m_screening,m_screening_driver
!!      m_sigma_driver,m_suscep_stat,m_vhxc_me,m_wfd,m_wfk,m_wfk_analyze,mrgscr
!!      ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_mpi_enreg(MPI_enreg)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 if (associated(mpi_enreg%distribfft)) then
   call destroy_distribfft(mpi_enreg%distribfft)
   ABI_FREE(mpi_enreg%distribfft)
   nullify(mpi_enreg%distribfft)
 end if

 ABI_SFREE(mpi_enreg%proc_distrb)
 ABI_SFREE(mpi_enreg%kptdstrb)
 ABI_SFREE(mpi_enreg%kpt_loc2fbz_sp)
 ABI_SFREE(mpi_enreg%kpt_loc2ibz_sp)
 ABI_SFREE(mpi_enreg%mkmem)
 ABI_SFREE(mpi_enreg%my_kpttab)
 if (associated(mpi_enreg%my_atmtab)) then
   ABI_FREE(mpi_enreg%my_atmtab)
   nullify(mpi_enreg%my_atmtab)
 end if
 ABI_SFREE(mpi_enreg%distrb_pert)
 ABI_SFREE(mpi_enreg%distrb_img)
 ABI_SFREE(mpi_enreg%my_imgtab)
 ABI_SFREE(mpi_enreg%my_kgtab)
 ABI_SFREE(mpi_enreg%distrb_hf)

!Do not deallocate wavelet denspot distribution arrays, they are handled by BigDFT.

end subroutine destroy_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/copy_mpi_enreg
!! NAME
!! copy_mpi_enreg
!!
!! FUNCTION
!!  Copy a mpi_enreg datastructure into another
!!
!! INPUTS
!!  MPI_enreg1<MPI_type>=input mpi_enreg datastructure
!!
!! OUTPUT
!!  MPI_enreg2<MPI_type>=output mpi_enreg datastructure
!!
!! PARENTS
!!      m_fft_prof,m_fock,m_gwls_hamiltonian,m_inwffil,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_mpi_enreg(MPI_enreg1,MPI_enreg2)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg1
 type(MPI_type),intent(out) :: MPI_enreg2

!Local variables-------------------------------
!scalars
 integer :: sz1,sz2,sz3

! *********************************************************************

!scalars
 mpi_enreg2%comm_world=mpi_enreg1%comm_world
 mpi_enreg2%me=mpi_enreg1%me
 mpi_enreg2%nproc=mpi_enreg1%nproc
 mpi_enreg2%paral_spinor=mpi_enreg1%paral_spinor
 mpi_enreg2%paralbd=mpi_enreg1%paralbd
 mpi_enreg2%me_fft=mpi_enreg1%me_fft
 mpi_enreg2%me_band=mpi_enreg1%me_band
 mpi_enreg2%nproc_fft=mpi_enreg1%nproc_fft
 mpi_enreg2%paral_kgb=mpi_enreg1%paral_kgb
 mpi_enreg2%me_g0=mpi_enreg1%me_g0
 mpi_enreg2%paral_pert=mpi_enreg1%paral_pert
 mpi_enreg2%me_pert=mpi_enreg1%me_pert
 mpi_enreg2%nproc_pert=mpi_enreg1%nproc_pert
 mpi_enreg2%comm_pert=mpi_enreg1%comm_pert
 mpi_enreg2%comm_bandfft=mpi_enreg1%comm_bandfft
 mpi_enreg2%comm_band=mpi_enreg1%comm_band
 mpi_enreg2%comm_fft=mpi_enreg1%comm_fft
 mpi_enreg2%nproc_band=mpi_enreg1%nproc_band
 mpi_enreg2%comm_bandspinorfft=mpi_enreg1%comm_bandspinorfft
 mpi_enreg2%comm_kpt=mpi_enreg1%comm_kpt
 mpi_enreg2%me_kpt=mpi_enreg1%me_kpt
 mpi_enreg2%nproc_kpt=mpi_enreg1%nproc_kpt
 mpi_enreg2%my_isppoltab=mpi_enreg1%my_isppoltab
 mpi_enreg2%my_natom=mpi_enreg1%my_natom
 mpi_enreg2%comm_atom=mpi_enreg1%comm_atom
 mpi_enreg2%nproc_atom=mpi_enreg1%nproc_atom
 mpi_enreg2%comm_kptband=mpi_enreg1%comm_kptband
 mpi_enreg2%bandpp=mpi_enreg1%bandpp
 mpi_enreg2%paral_img=mpi_enreg1%paral_img
 mpi_enreg2%comm_img=mpi_enreg1%comm_img
 mpi_enreg2%me_img=mpi_enreg1%me_img
 mpi_enreg2%nproc_img=mpi_enreg1%nproc_img
 mpi_enreg2%comm_cell=mpi_enreg1%comm_cell
 mpi_enreg2%comm_cell_pert=mpi_enreg1%comm_cell_pert
 mpi_enreg2%me_cell=mpi_enreg1%me_cell
 mpi_enreg2%nproc_cell=mpi_enreg1%nproc_cell
 mpi_enreg2%nproc_spinor=mpi_enreg1%nproc_spinor
 mpi_enreg2%me_spinor=mpi_enreg1%me_spinor
 mpi_enreg2%comm_spinorfft=mpi_enreg1%comm_spinorfft
 mpi_enreg2%me_wvl      =mpi_enreg1%me_wvl
 mpi_enreg2%nproc_wvl   =mpi_enreg1%nproc_wvl
 mpi_enreg2%comm_wvl    =mpi_enreg1%comm_wvl
 mpi_enreg2%me_hf      =mpi_enreg1%me_hf
 mpi_enreg2%nproc_hf   =mpi_enreg1%nproc_hf
 mpi_enreg2%comm_hf    =mpi_enreg1%comm_hf
 mpi_enreg2%paral_hf=mpi_enreg1%paral_hf
 mpi_enreg2%pw_unbal_thresh=mpi_enreg1%pw_unbal_thresh

!pointers
 if (associated(mpi_enreg1%distribfft)) then
   if (.not.associated(mpi_enreg2%distribfft)) then
     ABI_MALLOC(mpi_enreg2%distribfft,)
   end if
   call copy_distribfft(mpi_enreg1%distribfft,mpi_enreg2%distribfft)
 end if

 if (allocated(mpi_enreg1%proc_distrb)) then
   sz1=size(mpi_enreg1%proc_distrb,1)
   sz2=size(mpi_enreg1%proc_distrb,2)
   sz3=size(mpi_enreg1%proc_distrb,3)
   ABI_MALLOC(mpi_enreg2%proc_distrb,(sz1,sz2,sz3))
   mpi_enreg2%proc_distrb=mpi_enreg1%proc_distrb
 end if
 if (allocated(mpi_enreg1%kptdstrb)) then
   sz1=size(mpi_enreg1%kptdstrb,1)
   sz2=size(mpi_enreg1%kptdstrb,2)
   sz3=size(mpi_enreg1%kptdstrb,3)
   ABI_MALLOC(mpi_enreg2%kptdstrb,(sz1,sz2,sz3))
   mpi_enreg2%kptdstrb=mpi_enreg1%kptdstrb
 end if
 if (allocated(mpi_enreg1%kpt_loc2fbz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2fbz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2fbz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2fbz_sp,3)
   ABI_MALLOC(mpi_enreg2%kpt_loc2fbz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2fbz_sp=mpi_enreg1%kpt_loc2fbz_sp
 end if
 if (allocated(mpi_enreg1%kpt_loc2ibz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2ibz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2ibz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2ibz_sp,3)
   ABI_MALLOC(mpi_enreg2%kpt_loc2ibz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2ibz_sp=mpi_enreg1%kpt_loc2ibz_sp
 end if
 if (allocated(mpi_enreg1%mkmem)) then
   ABI_MALLOC(mpi_enreg2%mkmem,(0:size(mpi_enreg1%mkmem,1)-1))
   mpi_enreg2%mkmem=mpi_enreg1%mkmem
 end if
 if (associated(mpi_enreg1%my_atmtab)) then
   ABI_MALLOC(mpi_enreg2%my_atmtab,(size(mpi_enreg1%my_atmtab)))
   mpi_enreg2%my_atmtab=mpi_enreg1%my_atmtab
 else
   nullify(mpi_enreg2%my_atmtab)
 end if
 if (allocated(mpi_enreg1%my_kgtab)) then
   sz1=size(mpi_enreg1%my_kgtab,1)
   sz2=size(mpi_enreg1%my_kgtab,2)
   ABI_MALLOC(mpi_enreg2%my_kgtab,(sz1,sz2))
   mpi_enreg2%my_kgtab=mpi_enreg1%my_kgtab
 end if
 if (allocated(mpi_enreg1%distrb_pert)) then
   ABI_MALLOC(mpi_enreg2%distrb_pert,(size(mpi_enreg1%distrb_pert)))
   mpi_enreg2%distrb_pert=mpi_enreg1%distrb_pert
 end if
 if (allocated(mpi_enreg1%distrb_img)) then
   ABI_MALLOC(mpi_enreg2%distrb_img,(size(mpi_enreg1%distrb_img)))
   mpi_enreg2%distrb_img=mpi_enreg1%distrb_img
 end if
 if (allocated(mpi_enreg1%my_imgtab)) then
   ABI_MALLOC(mpi_enreg2%my_imgtab,(size(mpi_enreg1%my_imgtab)))
   mpi_enreg2%my_imgtab=mpi_enreg1%my_imgtab
 end if
 if (allocated(mpi_enreg1%distrb_hf)) then
   sz1=size(mpi_enreg1%distrb_hf,1)
   sz2=size(mpi_enreg1%distrb_hf,2)
   sz3=size(mpi_enreg1%distrb_hf,3)
   ABI_MALLOC(mpi_enreg2%distrb_hf,(sz1,sz2,sz3))
   mpi_enreg2%distrb_hf=mpi_enreg1%distrb_hf
 end if

!Optional pointers
 if (allocated(mpi_enreg1%my_kpttab)) then
   ABI_MALLOC(mpi_enreg2%my_kpttab,(size(mpi_enreg1%my_kpttab)))
   mpi_enreg2%my_kpttab=mpi_enreg1%my_kpttab
 end if

!Do not copy wavelet pointers, just associate.
 mpi_enreg2%nscatterarr => mpi_enreg1%nscatterarr
 mpi_enreg2%ngatherarr => mpi_enreg1%ngatherarr

end subroutine copy_mpi_enreg
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/set_mpi_enreg_fft
!! NAME
!! set_mpi_enreg_fft
!!
!! FUNCTION
!!  Set the content of a MPI datastructure in order to call fourwf/fourdp
!!  (in view of a wrapper for these routines)
!!
!! INPUTS
!!  me_g0=1 if the current process treat the g=0 plane-wave
!!  comm_fft= MPI communicator over FFT components
!!  paral_kgb= flag used to activate "band-FFT" parallelism
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=FFT pointer/flags intialized
!!
!! PARENTS
!!      m_atm2fft,m_paw_nhat,m_positron
!!
!! CHILDREN
!!
!! SOURCE

subroutine set_mpi_enreg_fft(MPI_enreg,comm_fft,distribfft,me_g0,paral_kgb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: me_g0,comm_fft,paral_kgb
 type(distribfft_type),intent(in),target :: distribfft
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 mpi_enreg%comm_fft=comm_fft
 mpi_enreg%paral_kgb=paral_kgb
 mpi_enreg%me_g0=me_g0
 mpi_enreg%nproc_fft=xmpi_comm_size(comm_fft)
 mpi_enreg%me_fft=xmpi_comm_rank(comm_fft)
 if (associated(mpi_enreg%distribfft)) then
   call destroy_distribfft(mpi_enreg%distribfft)
   ABI_FREE(mpi_enreg%distribfft)
 end if
 mpi_enreg%distribfft => distribfft

end subroutine set_mpi_enreg_fft
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/unset_mpi_enreg_fft
!! NAME
!! unset_mpi_enreg_fft
!!
!! FUNCTION
!!  Unset the content of a MPI datastructure used to call fourwf/fourdp
!!  (in view of a wrapper for these routines)
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  MPI_enreg<MPI_type>=FFT pointer/flags intialized
!!
!! PARENTS
!!      m_atm2fft,m_paw_nhat,m_positron
!!
!! CHILDREN
!!
!! SOURCE

subroutine unset_mpi_enreg_fft(MPI_enreg)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 mpi_enreg%me_g0=1
 mpi_enreg%comm_fft=xmpi_comm_self
 mpi_enreg%nproc_fft=1
 mpi_enreg%me_fft=0
 mpi_enreg%paral_kgb=0
 nullify(mpi_enreg%distribfft)

end subroutine unset_mpi_enreg_fft
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/ptabs_fourdp
!! NAME
!!  ptabs_fourdp
!!
!! FUNCTION
!!  Returns pointers to the tables used for the MPI FFT of densities and potentials (fourdp routine).
!!
!! NOTES
!!   1) These pointers are references to the internal tables stored in MPI_enreg hence
!!      *** DO NOT DEALLOCATE THE POINTERS YOU HAVE RECEIVED! ***
!!
!!   2) Client code should declare the pointers with the attribute ABI_CONTIGUOUS
!!      (this macro expands to F2008 CONTIGUOUS if the compiler supports it)
!!
!! INPUTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism.
!!  n2,n3=Number of FFT divisions along y and z
!!
!! OUTPUT
!!  fftn2_distrib(:)=  rank of the processor which own fft planes in 2nd dimension for fourdp
!!  ffti2_local(:) = local i2 indices in fourdp
!!  fftn3_distrib(:) = rank of the processor which own fft planes in 3rd dimension for fourdp
!!  ffti3_local(:) = local i3 indices in fourdp
!!
!! PARENTS
!!      m_dens,m_dfpt_elt,m_fft,m_fock,m_ioarr,m_mkcore,m_mklocl
!!      m_mklocl_realspace,m_mkrho,m_multipoles,m_nucprop,m_positron,m_prcref
!!      m_spacepar,m_stress,m_xctk
!!
!! CHILDREN
!!
!! SOURCE

subroutine ptabs_fourdp(MPI_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n2,n3
 type(MPI_type),intent(in) :: MPI_enreg
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

!Local variables-------------------------------
!scalars
 logical :: grid_found

! *********************************************************************

 grid_found=.false.

 ! Get the distrib associated with this fft_grid => for i2 and i3 planes
 grid_found=.false.
 if (n2== mpi_enreg%distribfft%n2_coarse) then
   if( n3 == size(mpi_enreg%distribfft%tab_fftdp3_distrib) )then
     fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2_distrib
     ffti2_local   => mpi_enreg%distribfft%tab_fftdp2_local
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3_distrib
     ffti3_local   => mpi_enreg%distribfft%tab_fftdp3_local
     grid_found=.true.
   end if
 end if

 if((n2 == mpi_enreg%distribfft%n2_fine).and.(.not.(grid_found))) then
   if( n3 == size(mpi_enreg%distribfft%tab_fftdp3dg_distrib) )then
     fftn2_distrib => mpi_enreg%distribfft%tab_fftdp2dg_distrib
     ffti2_local   => mpi_enreg%distribfft%tab_fftdp2dg_local
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3dg_distrib
     ffti3_local   => mpi_enreg%distribfft%tab_fftdp3dg_local
     grid_found=.true.
   end if
 end if

 if (.not.grid_found) then
   ABI_BUG(sjoin("Unable to find an allocated distrib for this fft grid with n2, n3 = ", ltoa([n2, n3])))
 end if

end subroutine ptabs_fourdp
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/ptabs_fourwf
!! NAME
!!  ptabs_fourwf
!!
!! FUNCTION
!!  Returns pointers to the tables used for the MPI FFT of the wavefunctions (fourwf routine).
!!
!! NOTES
!!   1) These pointers are references to the internal tables stored in MPI_enreg hence
!!      *** DO NOT DEALLOCATE THE POINTERS YOU HAVE RECEIVED! ***
!!
!!   2) Client code should declare the pointers with the attribute ABI_CONTIGUOUS
!!      (this macro expands to F2008 CONTIGUOUS if the compiler supports it)
!!
!! INPUTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism.
!!  n2,n3=Number of FFT divisions along y and z
!!
!! OUTPUT
!!  fftn2_distrib(:)=  rank of the processors which own fft planes in 2nd dimension for fourwf
!!  ffti2_local(:) = local i2 indices in fourwf
!!  fftn3_distrib(:) = rank of the processors which own fft planes in 3rd dimension for fourwf
!!  ffti3_local(:) = local i3 indices in fourwf
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!
!! SOURCE

subroutine ptabs_fourwf(MPI_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n2,n3
 type(MPI_type),intent(in) :: MPI_enreg
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

!Local variables-------------------------------
!scalars
 logical :: grid_found

! *********************************************************************

 grid_found=.false.

 ! Get the distrib associated with this fft_grid => for i2 and i3 planes
 if (n2 == mpi_enreg%distribfft%n2_coarse) then
   if (n3 == size(mpi_enreg%distribfft%tab_fftdp3_distrib))then
     fftn2_distrib => mpi_enreg%distribfft%tab_fftwf2_distrib
     ffti2_local   => mpi_enreg%distribfft%tab_fftwf2_local
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3_distrib
     ffti3_local   => mpi_enreg%distribfft%tab_fftdp3_local
     grid_found=.true.
   end if
 end if

 if((n2 == mpi_enreg%distribfft%n2_fine).and.(.not.(grid_found))) then
   if (n3 == size(mpi_enreg%distribfft%tab_fftdp3dg_distrib) )then
     fftn2_distrib => mpi_enreg%distribfft%tab_fftwf2dg_distrib
     ffti2_local   => mpi_enreg%distribfft%tab_fftwf2dg_local
     fftn3_distrib => mpi_enreg%distribfft%tab_fftdp3dg_distrib
     ffti3_local   => mpi_enreg%distribfft%tab_fftdp3dg_local
     grid_found=.true.
   end if
 end if

 if(.not. grid_found) then
   ABI_BUG(sjoin("Unable to find an allocated distrib for this fft grid", ltoa([n2, n3])))
 end if

end subroutine ptabs_fourwf
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/mpi_distrib_is_ok
!! NAME
!!  mpi_distrib_is_ok
!!
!! FUNCTION
!!  Check if a MPI datastructure contains number of processors
!!  compatible (in terms of efficiency) with the number of spins/k-points/bands
!!
!! INPUTS
!!  MPI_enreg<MPI_type>=Datatype gathering information on the parallelism
!!  nband=number of bands
!!  nkpt=number of k-points
!!  nptk_current_proc=number of k-points handled by current MPI process
!!  nsppol= number of spins (1 or 2)
!!
!! OUTPUT
!!  mpi_distrib_is_ok (current function)=TRUE if the current MPI distribution is optimal
!!                                       FALSE otherwise
!!  [msg]= -optional- warning message to be printed out
!!
!! PARENTS
!!  driver,mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

logical function mpi_distrib_is_ok(MPI_enreg,nband,nkpt,nkpt_current_proc,nsppol,msg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nkpt_current_proc,nsppol
 type(MPI_type),intent(in) :: MPI_enreg
 character(len=*),optional,intent(out) :: msg

! *********************************************************************

 mpi_distrib_is_ok=.true.

 if (MPI_enreg%paralbd==0) then
   if (MPI_enreg%nproc_kpt-floor(nsppol*nkpt*one/nkpt_current_proc)>=nkpt_current_proc) then
     mpi_distrib_is_ok=.false.
     if (present(msg)) then
       write(msg,'(a,i0,4a,i0,3a)') &
        'Your number of spins*k-points (=',nsppol*nkpt,') ',&
        'will not distribute correctly',ch10, &
        'with the current number of processors (=',MPI_enreg%nproc_kpt,').',ch10,&
        'You will leave some empty.'
     end if
   end if
 else
   if (mod(nband,max(1,MPI_enreg%nproc_kpt/(nsppol*nkpt)))/=0) then
     mpi_distrib_is_ok=.false.
     if (present(msg)) then
       write(msg,'(a,i0,2a,i0,4a,i0,7a)')&
        'Your number of spins*k-points (=',nsppol*nkpt,') ',&
         'and bands (=',nband,') ',&
         'will not distribute correctly',ch10,&
         'with the current number of processors (=',MPI_enreg%nproc_kpt,').',ch10,&
         'You will leave some empty.'
     end if
   end if
 end if

end function mpi_distrib_is_ok
!!***

!!****f* ABINIT/proc_distrb_cycle
!! NAME
!!  proc_distrb_cycle
!!
!! FUNCTION
!!  test a condition to cycle
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function proc_distrb_cycle(distrb,ikpt,iband1,iband2,isppol,me)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,iband1,iband2,isppol,me
 integer,allocatable,intent(in) :: distrb(:,:,:)
 logical :: proc_distrb_cycle

! *************************************************************************

 proc_distrb_cycle=.false.
 if (allocated(distrb)) then
   if (isppol==-1) then
     proc_distrb_cycle=(minval(abs(distrb(ikpt,iband1:iband2,:)-me))/=0)
   else
     proc_distrb_cycle=(minval(abs(distrb(ikpt,iband1:iband2,isppol)-me))/=0)
   end if
 end if

end function proc_distrb_cycle
!!***

!!****f* ABINIT/initmpi_world
!! NAME
!!  initmpi_world
!!
!! FUNCTION
!!  %comm_world is redifined for the number of processors on which ABINIT is launched
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_world(mpi_enreg,nproc)

!Arguments ------------------------------------
 integer, intent(in)::nproc
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: ii
!arrays
 integer,allocatable :: ranks(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 if(nproc==mpi_enreg%nproc) return

 ABI_MALLOC(ranks,(0:nproc-1))
 ranks(0:nproc-1)=(/((ii),ii=0,nproc-1)/)
 mpi_enreg%comm_world=xmpi_subcomm(xmpi_world,nproc,ranks)
 ABI_FREE(ranks)

 if(mpi_enreg%me<nproc)  then
   mpi_enreg%me=xmpi_comm_rank(mpi_enreg%comm_world)
   mpi_enreg%nproc=xmpi_comm_size(mpi_enreg%comm_world)
   call abi_io_redirect(new_io_comm=mpi_enreg%comm_world)
   call libpaw_write_comm_set(mpi_enreg%comm_world)
 else
   mpi_enreg%me=-1
 end if

 DBG_EXIT("COLL")

end subroutine initmpi_world
!!***

!!****f* ABINIT/initmpi_seq
!! NAME
!!  initmpi_seq
!!
!! FUNCTION
!!  Initializes the MPI information for sequential use.
!!
!! INPUTS
!!
!! OUTPUT
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      cut3d,fftprof,m_atm2fft,m_bethe_salpeter,m_cut3d,m_ddk,m_dfpt_nstwf
!!      m_dvdb,m_eph_driver,m_epjdos,m_fft,m_fft_prof,m_fftcore,m_gsphere
!!      m_ioarr,m_ksdiago,m_kxc,m_mlwfovlp_qp,m_mpinfo,m_paw_nhat,m_paw_optics
!!      m_pawpwij,m_positron,m_ppmodel,m_prcref,m_scfcv_core,m_screening
!!      m_screening_driver,m_sigma_driver,m_suscep_stat,m_vhxc_me,m_wfd,m_wfk
!!      m_wfk_analyze,mrgscr,ujdet
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_seq(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type),intent(out) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

!Set default seq values for scalars
 mpi_enreg%bandpp=1
 mpi_enreg%me=0
 mpi_enreg%me_band=0
 mpi_enreg%me_cell=0
 mpi_enreg%me_fft=0
 mpi_enreg%me_g0=1
 mpi_enreg%me_img=0
 mpi_enreg%me_hf=0
 mpi_enreg%me_kpt=0
 mpi_enreg%me_pert=0
 mpi_enreg%me_spinor=0
 mpi_enreg%me_wvl=0
 mpi_enreg%my_natom=0       ! Should be natom
 mpi_enreg%my_isppoltab=0   ! Should be (1,0) if nsppol=1 or (1,1) if nsppol=2
 mpi_enreg%ngfft3_ionic=1
 mpi_enreg%my_nimage=1
 mpi_enreg%nproc=1
 mpi_enreg%nproc_atom=1
 mpi_enreg%nproc_band=1
 mpi_enreg%nproc_cell=1
 mpi_enreg%nproc_fft=1
 mpi_enreg%nproc_img=1
 mpi_enreg%nproc_hf=1
 mpi_enreg%nproc_kpt=1
 mpi_enreg%nproc_pert=1
 mpi_enreg%nproc_spinor=1
 mpi_enreg%nproc_wvl=1
 mpi_enreg%paralbd=0
 mpi_enreg%paral_img=0
 mpi_enreg%paral_hf=0
 mpi_enreg%paral_kgb=0
 mpi_enreg%paral_pert=0
 mpi_enreg%paral_spinor=0
 mpi_enreg%pw_unbal_thresh=-1._dp

!Set default seq values for communicators
 mpi_enreg%comm_world          = xmpi_world
 mpi_enreg%comm_atom           = xmpi_comm_self
 mpi_enreg%comm_band           = xmpi_comm_self
 mpi_enreg%comm_bandspinor     = xmpi_comm_self
 mpi_enreg%comm_bandfft        = xmpi_comm_self
 mpi_enreg%comm_bandspinorfft  = xmpi_comm_self
 mpi_enreg%comm_cell           = xmpi_comm_self
 mpi_enreg%comm_cell_pert      = xmpi_comm_self
 mpi_enreg%comm_fft            = xmpi_comm_self
 mpi_enreg%comm_hf             = xmpi_comm_self
 mpi_enreg%comm_img            = xmpi_comm_self
 mpi_enreg%comm_kpt            = xmpi_comm_self
 mpi_enreg%comm_kptband        = xmpi_comm_self
 mpi_enreg%comm_pert           = xmpi_comm_self
 mpi_enreg%comm_spinor         = xmpi_comm_self
 mpi_enreg%comm_spinorfft      = xmpi_comm_self
 mpi_enreg%comm_wvl            = xmpi_comm_self

!Nullify all pointers
 call nullify_mpi_enreg(mpi_enreg)

!Allocate and nullify distribfft datastructure
! This is not good since distribfft is not initialized here (even with 0s).
! It can be dangerous if use with no care (Valgrind might complain)
 ABI_MALLOC(mpi_enreg%distribfft,)

 DBG_EXIT("COLL")

end subroutine initmpi_seq
!!***

!!****f* ABINIT/initmpi_atom
!! NAME
!!  initmpi_atom
!!
!! FUNCTION
!!  Initializes the mpi information for parallelism over atoms (PAW).
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  mpi_enreg= information about MPI parallelization
!!
!! OUTPUT
!!  mpi_enreg= information about MPI parallelization
!!    comm_atom                 =communicator over atoms
!!    nproc_atom                =size of the communicator over atoms
!!    my_natom                  =number of atoms treated by current proc
!!    my_atmtab(mpi_enreg%natom)=indexes of the atoms treated by current processor
!!
!! PARENTS
!!      m_mpi_setup,m_paral_pert
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_atom(dtset,mpi_enreg)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: msg
 integer :: iatom

! ***********************************************************************

 DBG_ENTER("COLL")

 mpi_enreg%nproc_atom=1
 mpi_enreg%comm_atom=xmpi_comm_self
 mpi_enreg%my_natom=dtset%natom
 if (associated(mpi_enreg%my_atmtab))then
   ABI_FREE(mpi_enreg%my_atmtab)
 end if
 nullify(mpi_enreg%my_atmtab)

 if (xmpi_paral==0) then
   mpi_enreg%nproc_atom=0
   ABI_MALLOC(mpi_enreg%my_atmtab,(0))
   return
 end if

!Check compatibility
 if (dtset%paral_atom>0) then
   msg=''
   if (dtset%usepaw==0)  msg= 'Parallelisation over atoms not compatible with usepaw=0 !'
   if (dtset%usedmft==1) msg=' Parallelisation over atoms not compatible with usedmft=1 !'
   if (dtset%usewvl==1)  msg= 'Parallelisation over atoms not compatible with usewvl=1 !'
   if (dtset%prtden>1.and.dtset%paral_kgb<=0) &
&   msg= 'Parallelisation over atoms not compatible with prtden>1 (PAW AE densities) !'
   if (dtset%optdriver/=RUNL_GSTATE.and.dtset%optdriver/=RUNL_RESPFN) &
&   msg=' Parallelisation over atoms only compatible with GS or RF !'
   if (dtset%macro_uj/=0)msg=' Parallelisation over atoms not compatible with macro_uj!=0 !'
   if (msg/='') then
     ABI_ERROR(msg)
   end if
 end if

 if (mpi_enreg%comm_atom==xmpi_comm_null) then
   mpi_enreg%nproc_atom=0;mpi_enreg%my_natom=0
   ABI_MALLOC(mpi_enreg%my_atmtab,(0))
   return
 end if

 if (dtset%paral_atom>0) then

!  Build correct atom communicator
   if (dtset%optdriver==RUNL_GSTATE.and.dtset%paral_kgb==1) then
     mpi_enreg%comm_atom=mpi_enreg%comm_kptband
   else
     mpi_enreg%comm_atom=mpi_enreg%comm_cell
   end if

!  Get number of processors sharing the atomic data distribution
   mpi_enreg%nproc_atom=xmpi_comm_size(mpi_enreg%comm_atom)

!  Get local number of atoms
   call get_my_natom(mpi_enreg%comm_atom,mpi_enreg%my_natom,dtset%natom)
   paral_atom=(mpi_enreg%my_natom/=dtset%natom)

!  Build atom table
   if (mpi_enreg%my_natom>0.and.paral_atom) then
     my_atmtab_allocated=.false.
     call get_my_atmtab(mpi_enreg%comm_atom,mpi_enreg%my_atmtab,my_atmtab_allocated, &
&     paral_atom,dtset%natom)
   else if (.not.paral_atom) then
     ABI_MALLOC(mpi_enreg%my_atmtab,(dtset%natom))
     mpi_enreg%my_atmtab(1:dtset%natom)=(/(iatom, iatom=1,dtset%natom)/)
   else if (mpi_enreg%my_natom==0) then
     ABI_MALLOC(mpi_enreg%my_atmtab,(0))
   end if

 end if

 DBG_EXIT("COLL")

end subroutine initmpi_atom
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_atom
!! NAME
!!  clnmpi_atom
!!
!! FUNCTION
!!  Cleans-up the mpi information for the parallelism over atoms (PAW).
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_atom(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

 if (xmpi_paral==0) return

 if (mpi_enreg%comm_atom/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_atom)
   mpi_enreg%comm_atom=xmpi_comm_null
 end if

 if(associated(mpi_enreg%my_atmtab)) then
   ABI_FREE(mpi_enreg%my_atmtab)
 end if

 mpi_enreg%nproc_atom=1
 mpi_enreg%my_natom=0 ! should be natom

 DBG_EXIT("COLL")

end subroutine clnmpi_atom
!!***

!!****f* ABINIT/initmpi_grid
!! NAME
!!  initmpi_grid
!!
!! FUNCTION
!!  Initializes the MPI information for the grid:
!!    * 2D if parallization KPT/FFT (paral_kgb == 0 & MPI)
!!    * 3D if parallization KPT/FFT/BAND (paral_kgb == 1 & MPI)
!!    * 2D in case of an Hartree-Fock calculation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_grid(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: nproc,nproc_eff,spacecomm
 character(len=500) :: msg
#if defined HAVE_MPI
 integer :: commcart_4d,dimcart,ierr,me_cart_4d
 integer :: commcart_2d,me_cart_2d
 logical :: reorder
!arrays
 integer,allocatable :: coords(:),sizecart(:)
 logical,allocatable :: periode(:), keepdim(:)
#endif

! *********************************************************************

 DBG_ENTER("COLL")

!Select the correct "world" communicator"
 nproc=mpi_enreg%nproc_cell
 if(mpi_enreg%paral_pert==1) nproc=mpi_enreg%nproc_cell
 spacecomm=mpi_enreg%comm_cell

!Fake values for null communicator
 if (nproc==0) then
   mpi_enreg%nproc_fft    = 0
   mpi_enreg%nproc_band   = 0
   mpi_enreg%nproc_hf    = 0
   mpi_enreg%nproc_kpt    = 0
   mpi_enreg%nproc_spinor   = 0
   mpi_enreg%comm_fft            = xmpi_comm_null
   mpi_enreg%comm_band           = xmpi_comm_null
   mpi_enreg%comm_hf             = xmpi_comm_null
   mpi_enreg%comm_kpt            = xmpi_comm_null
   mpi_enreg%comm_kptband        = xmpi_comm_null
   mpi_enreg%comm_spinor         = xmpi_comm_null
   mpi_enreg%comm_bandspinor     = xmpi_comm_null
   mpi_enreg%comm_spinorfft      = xmpi_comm_null
   mpi_enreg%comm_bandfft        = xmpi_comm_null
   mpi_enreg%comm_bandspinorfft  = xmpi_comm_null
   mpi_enreg%bandpp       = 1
   return
 end if

#if defined HAVE_MPI
 if (mpi_enreg%paral_hf==0) then
   ! either the option Fock exchange is not active or there is no parallelization on Fock exchange calculation.

   if (mpi_enreg%nproc_spinor>1) mpi_enreg%paral_spinor=1

    !Effective number of processors used for the grid
   nproc_eff=mpi_enreg%nproc_fft*mpi_enreg%nproc_band *mpi_enreg%nproc_kpt*mpi_enreg%nproc_spinor
   if(nproc_eff/=nproc) then
     write(msg,'(4a,5(a,i0))') &
      '  The number of band*FFT*kpt*spinor processors, npband*npfft*npkpt*npspinor should be',ch10,&
      '  equal to the total number of processors, nproc.',ch10,&
      '  However, npband   =',mpi_enreg%nproc_band,&
      '           npfft    =',mpi_enreg%nproc_fft,&
      '           npkpt    =',mpi_enreg%nproc_kpt,&
      '           npspinor =',mpi_enreg%nproc_spinor,&
      '       and nproc    =',nproc
     ABI_WARNING(msg)
   end if

   !Nothing to do if only 1 proc
   if (nproc_eff==1) return

   ! Initialize the communicator for Hartree-Fock to xmpi_comm_self
   mpi_enreg%me_hf =0
   mpi_enreg%comm_hf=xmpi_comm_self

   if(mpi_enreg%paral_kgb==0) then
     mpi_enreg%me_fft =0
     mpi_enreg%me_band=0
     mpi_enreg%me_kpt =mpi_enreg%me_cell
     mpi_enreg%me_spinor=0
     mpi_enreg%comm_fft=xmpi_comm_self
     mpi_enreg%comm_band=xmpi_comm_self
     mpi_enreg%comm_kpt=mpi_enreg%comm_cell
     mpi_enreg%comm_spinor=xmpi_comm_self
     mpi_enreg%comm_bandspinor=xmpi_comm_self
     mpi_enreg%comm_kptband=mpi_enreg%comm_cell
     mpi_enreg%comm_spinorfft=xmpi_comm_self
     mpi_enreg%comm_bandfft=xmpi_comm_self
     mpi_enreg%comm_bandspinorfft=xmpi_comm_self
   else
     !  CREATE THE 4D GRID
     !  ==================================================

     !  Create the global cartesian 4D- communicator
     !  valgrind claims this is not deallocated in test v5/72
     !  Can someone knowledgable check?
     dimcart=4
     ABI_MALLOC(sizecart,(dimcart))
     ABI_MALLOC(periode,(dimcart))
!    MT 2012-june: change the order of the indexes; not sure this is efficient
!    (not efficient on TGCC-Curie).
     sizecart(1)=mpi_enreg%nproc_kpt  ! mpi_enreg%nproc_kpt
     sizecart(2)=mpi_enreg%nproc_band ! mpi_enreg%nproc_band
     sizecart(3)=mpi_enreg%nproc_spinor ! mpi_enreg%nproc_spinor
     sizecart(4)=mpi_enreg%nproc_fft  ! mpi_enreg%nproc_fft
     periode(:)=.false.;reorder=.false.
     call MPI_CART_CREATE(spacecomm,dimcart,sizecart,periode,reorder,commcart_4d,ierr)
     ABI_FREE(periode)
     ABI_FREE(sizecart)

!    Find the index and coordinates of the current processor
     call MPI_COMM_RANK(commcart_4d, me_cart_4d, ierr)
     ABI_MALLOC(coords,(dimcart))
     call MPI_CART_COORDS(commcart_4d, me_cart_4d,dimcart,coords,ierr)
     mpi_enreg%me_kpt =coords(1)
     mpi_enreg%me_band=coords(2)
     mpi_enreg%me_spinor=coords(3)
     mpi_enreg%me_fft =coords(4)
     ABI_FREE(coords)

     ABI_MALLOC(keepdim,(dimcart))

!    Create the communicator for fft distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.false.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_fft,ierr)

!    Create the communicator for band distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_band,ierr)

!    Create the communicator for kpt distribution
     keepdim(1)=.true.
     keepdim(2)=.false.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_kpt,ierr)

!    Create the communicator for spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.true.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_spinor,ierr)

!    Create the communicator for band-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.true.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandspinor,ierr)
     if (ierr /= MPI_SUCCESS ) then
       call xmpi_abort(mpi_enreg%comm_world,ierr)
     end if

!    Create the communicator for kpt-band distribution
     keepdim(1)=.true.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_kptband,ierr)

!    Create the communicator for fft-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.true.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_spinorfft,ierr)

!    Create the communicator for fft-band distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandfft,ierr)

!    Create the communicator for fft-band-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.true.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandspinorfft,ierr)

     ABI_FREE(keepdim)
     call xmpi_comm_free(commcart_4d)
   end if

   !Write some data
   !write(msg,'(a,4i5)') 'npfft, npband, npspinor and npkpt: ',&
   !mpi_enreg%nproc_fft,mpi_enreg%nproc_band, mpi_enreg%nproc_spinor,mpi_enreg%nproc_kpt
   !call wrtout(std_out,msg,'COLL')
   !write(msg,'(a,4i5)') 'me_fft, me_band, me_spinor , me_kpt: ',&
   !mpi_enreg%me_fft,mpi_enreg%me_band,mpi_enreg%me_spinor, mpi_enreg%me_kpt
   !call wrtout(std_out,msg,'COLL')

 else ! paral_hf==1
!* Option Hartree-Fock is active and more than 1 processor is dedicated to the parallelization over occupied states.

!* Initialize the values of fft, band and spinor communicators, as in the case paral_kgb==0.
   mpi_enreg%me_fft =0
   mpi_enreg%me_band=0
   mpi_enreg%me_spinor=0
   mpi_enreg%comm_fft=xmpi_comm_self
   mpi_enreg%comm_band=xmpi_comm_self
   mpi_enreg%comm_spinor=xmpi_comm_self
   mpi_enreg%comm_bandspinor=xmpi_comm_self
   mpi_enreg%comm_kptband=mpi_enreg%comm_cell
   mpi_enreg%comm_spinorfft=xmpi_comm_self
   mpi_enreg%comm_bandfft=xmpi_comm_self
   mpi_enreg%comm_bandspinorfft=xmpi_comm_self

!* Create the global cartesian 2D- communicator
   dimcart=2
   ABI_MALLOC(sizecart,(dimcart))
   ABI_MALLOC(periode,(dimcart))
   sizecart(1)=mpi_enreg%nproc_kpt  ! mpi_enreg%nproc_kpt
   sizecart(2)=mpi_enreg%nproc_hf   ! mpi_enreg%nproc_hf
   periode(:)=.false.;reorder=.false.
   call MPI_CART_CREATE(spacecomm,dimcart,sizecart,periode,reorder,commcart_2d,ierr)
   ABI_FREE(periode)
   ABI_FREE(sizecart)

!* Find the index and coordinates of the current processor
   call MPI_COMM_RANK(commcart_2d, me_cart_2d, ierr)
   ABI_MALLOC(coords,(dimcart))
   call MPI_CART_COORDS(commcart_2d, me_cart_2d,dimcart,coords,ierr)
   mpi_enreg%me_kpt =coords(1)
   mpi_enreg%me_hf=coords(2)
   ABI_FREE(coords)

   ABI_MALLOC(keepdim,(dimcart))

!* Create the communicator for kpt distribution
   keepdim(1)=.true.
   keepdim(2)=.false.
   call MPI_CART_SUB(commcart_2d, keepdim, mpi_enreg%comm_kpt,ierr)

!* Create the communicator for hf distribution
   keepdim(1)=.false.
   keepdim(2)=.true.
   call MPI_CART_SUB(commcart_2d, keepdim, mpi_enreg%comm_hf,ierr)

   ABI_FREE(keepdim)
   call xmpi_comm_free(commcart_2d)

!* Write some data
   write(msg,'(a,2(1x,i0))') 'nphf and npkpt: ',mpi_enreg%nproc_hf, mpi_enreg%nproc_kpt
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,2(1x,i0))') 'me_hf, me_kpt: ',mpi_enreg%me_hf, mpi_enreg%me_kpt
   call wrtout(std_out,msg,'COLL')
 end if
#endif

 DBG_EXIT("COLL")

end subroutine initmpi_grid
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_grid
!! NAME
!!  clnmpi_grid
!!
!! FUNCTION
!!  Cleans-up the mpi information for parallelism over grid (kpt/band/fft).
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_grid(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

 if (xmpi_paral==0) return

 if (mpi_enreg%comm_bandspinorfft/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_bandspinorfft)
   mpi_enreg%comm_bandspinorfft=xmpi_comm_null
 end if

 if (mpi_enreg%comm_bandfft/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_bandfft)
   mpi_enreg%comm_bandfft=xmpi_comm_null
 end if

 if (mpi_enreg%comm_spinorfft/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_spinorfft)
   mpi_enreg%comm_spinorfft=xmpi_comm_null
 end if

 if (mpi_enreg%comm_bandspinor/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_bandspinor)
   mpi_enreg%comm_bandspinor=xmpi_comm_null
 end if

 if (mpi_enreg%comm_kptband/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_kptband)
   mpi_enreg%comm_kptband=xmpi_comm_null
 end if

 if (mpi_enreg%comm_fft/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_fft)
   mpi_enreg%comm_fft=xmpi_comm_null
 end if

 if (mpi_enreg%comm_band/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_band)
   mpi_enreg%comm_band=xmpi_comm_null
 end if

 if (mpi_enreg%comm_spinor/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_spinor)
   mpi_enreg%comm_spinor=xmpi_comm_null
 end if

 if (mpi_enreg%comm_kpt/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_kpt)
   mpi_enreg%comm_kpt=xmpi_comm_null
 end if

 DBG_EXIT("COLL")

end subroutine clnmpi_grid
!!***

!!****f* ABINIT/initmpi_img
!! NAME
!!  initmpi_img
!!
!! FUNCTION
!!  Initializes the mpi information for parallelism over images of the cell (npimage>1).
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg= information about MPI parallelization
!!  option= see below
!!
!! OUTPUT
!!  mpi_enreg%my_nimage= number of images of the cell treated by current proc
!!  ===== if option==1 or option==-1
!!    mpi_enreg%my_imgtab= indexes of images of the cell treated by current proc
!!  ===== if option==2 or option==3 or option==-1
!!    mpi_enreg%comm_cell=Communicator over all processors treating the same image
!!    mpi_enreg%nproc_cell=size of comm_cell
!!    mpi_enreg%me_cell=my rank in comm_cell
!!  ===== if option==3 or option==-1
!!    mpi_enreg%comm_img=Communicator over all images
!!    mpi_enreg%nproc_img=size of comm_img
!!    mpi_enreg%me_img=my rank in comm_img
!!    mpi_enreg%distrb_img(:)=index of processor treating each image (in comm_img communicator)
!!
!! PARENTS
!!      m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_img(dtset,mpi_enreg,option)

!Arguments ------------------------------------
 integer,intent(in) :: option
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: imod,irank,iprocmax,iprocmin,jrank
 integer :: ndynimage_eff,nimage_eff,nproc_per_image,nrank
 logical,parameter :: debug=.false.
 integer,allocatable :: ranks(:)
 character(len=500) :: msg

!integer :: group_cell,ierr

! ***********************************************************************

 DBG_ENTER("COLL")
 if (option/=0) then
   mpi_enreg%comm_img=xmpi_comm_self
   mpi_enreg%comm_cell=mpi_enreg%comm_world
 end if

 if (xmpi_paral==1.and.dtset%npimage>1.and.dtset%npimage<=mpi_enreg%nproc.and. dtset%optdriver==RUNL_GSTATE) then

!  Activate flag for parallelization over images
   mpi_enreg%paral_img=1

   ndynimage_eff=dtset%ndynimage;if (dtset%ntimimage<=1) ndynimage_eff=0

!  Print several warnings
   if (option==0) then
     nimage_eff=max(ndynimage_eff,dtset%nimage-ndynimage_eff)
     if (dtset%npimage>nimage_eff) then
       write(unit=msg,fmt='(3a,i4,a,i4,4a)') &
        'The number of processors used for the parallelization',ch10,&
        ' over images (npimage=',dtset%npimage,&
        ') is greater than the number of dynamic (or static) images (',nimage_eff,') !',ch10,&
        ' This is inefficient.',ch10
       ABI_WARNING(msg)
     end if
     if (dtset%npimage>mpi_enreg%nproc) then
       write(unit=msg,fmt='(3a,i6,a,i4,4a)') &
        'The number of processors used for the parallelization',ch10,&
        ' over images (nproc=',mpi_enreg%nproc,&
        ') is smaller than npimage in input file (',dtset%npimage,&
        ')!',ch10,' This is unconsistent.',ch10
       ABI_ERROR(msg)
     end if
     if (mod(nimage_eff,dtset%npimage)/=0) then
       write(unit=msg,fmt='(3a,i4,a,i4,4a)') &
        'The number of processors used for the parallelization',ch10,&
        ' over images (npimage=',dtset%npimage,&
        ') does not divide the number of dynamic images (',nimage_eff,&
        ') !',ch10,' This is inefficient (charge unbalancing).',ch10
       ABI_WARNING(msg)
     end if
   end if

!  # of images treated by current proc
   nproc_per_image=mpi_enreg%nproc/dtset%npimage
   iprocmax=nproc_per_image*dtset%npimage-1
   if (mpi_enreg%me<=iprocmax) then
     mpi_enreg%my_nimage=(ndynimage_eff/dtset%npimage)+((dtset%nimage-ndynimage_eff)/dtset%npimage)
     imod=mod(ndynimage_eff,dtset%npimage)-1
     if (mpi_enreg%me/nproc_per_image<=imod) mpi_enreg%my_nimage=mpi_enreg%my_nimage+1
     imod=mod((dtset%nimage-ndynimage_eff),dtset%npimage)-1
     if (mpi_enreg%me/nproc_per_image<=imod) mpi_enreg%my_nimage=mpi_enreg%my_nimage+1
   else
     mpi_enreg%my_nimage=0
   end if
   if (option==1.or.option==-1) then
!    Indexes of images treated by current proc
     if (mpi_enreg%me<=iprocmax) then
       ABI_MALLOC(mpi_enreg%my_imgtab,(mpi_enreg%my_nimage))
       nrank=0
       imod=mpi_enreg%me/nproc_per_image+1;imod=mod(imod,dtset%npimage)
!      Dynamic images
       irank=0
       do jrank=1,dtset%nimage
         if (dtset%dynimage(jrank)/=0.and.dtset%ntimimage>1) then
           irank=irank+1
           if (mod(irank,dtset%npimage)==imod) then
             nrank=nrank+1
             mpi_enreg%my_imgtab(nrank)=jrank
           end if
         end if
       end do
!      Static images
       irank=0
       do jrank=1,dtset%nimage
         if (dtset%dynimage(jrank)==0.or.dtset%ntimimage<=1) then
           irank=irank+1
           if (mod(irank,dtset%npimage)==imod) then
             nrank=nrank+1
             mpi_enreg%my_imgtab(nrank)=jrank
           end if
         end if
       end do
       if (nrank/=mpi_enreg%my_nimage) then
         ABI_BUG('Error on nrank !')
       end if
!      Sort images by increasing index (this step is MANDATORY !!)
       ABI_MALLOC(ranks,(nrank))
       call sort_int(nrank,mpi_enreg%my_imgtab,ranks)
       ABI_FREE(ranks)
     else
       ABI_MALLOC(mpi_enreg%my_imgtab,(0))
     end if
   end if
   if (option==2.or.option==3.or.option==-1) then
!    Communicator over one image
     if (mpi_enreg%me<=iprocmax) then
       ABI_MALLOC(ranks,(nproc_per_image))
       iprocmin=(mpi_enreg%me/nproc_per_image)*nproc_per_image
       ranks=(/((iprocmin+irank-1),irank=1,nproc_per_image)/)
       mpi_enreg%comm_cell=xmpi_subcomm(mpi_enreg%comm_world,nproc_per_image,ranks)
       ABI_FREE(ranks)
       mpi_enreg%me_cell=xmpi_comm_rank(mpi_enreg%comm_cell)
       mpi_enreg%nproc_cell=nproc_per_image
       if (mpi_enreg%me_cell==0.and.mod(mpi_enreg%me,nproc_per_image)/=0) then
         ABI_BUG('Error on me_cell !')
       end if
     else
       mpi_enreg%comm_img=xmpi_comm_null
       mpi_enreg%nproc_cell=0
       mpi_enreg%me_cell=-1
     end if
   end if
   if (option==3.or.option==-1) then
!    Communicator over all images
     if (mpi_enreg%me<=iprocmax) then
       ABI_MALLOC(ranks,(dtset%npimage))
       iprocmin=mod(mpi_enreg%me,nproc_per_image)
       ranks=(/((iprocmin+(irank-1)*nproc_per_image),irank=1,dtset%npimage)/)
       mpi_enreg%comm_img=xmpi_subcomm(mpi_enreg%comm_world,dtset%npimage,ranks)
       ABI_FREE(ranks)
       mpi_enreg%me_img=xmpi_comm_rank(mpi_enreg%comm_img)
       mpi_enreg%nproc_img=dtset%npimage
       if (iprocmin==0.and.mpi_enreg%me_img==0.and.mpi_enreg%me/=0) then
         ABI_BUG('Error on me_img!')
       end if
       ABI_MALLOC(mpi_enreg%distrb_img,(dtset%nimage))
!      Dynamic images
       nrank=0
       do irank=1,dtset%nimage
         if (dtset%dynimage(irank)/=0.and.dtset%ntimimage>1) then
           nrank=nrank+1
           mpi_enreg%distrb_img(irank)=mod(nrank,dtset%npimage)-1
           if (mpi_enreg%distrb_img(irank)==-1) mpi_enreg%distrb_img(irank)=dtset%npimage-1
         end if
       end do
!      Static images
       nrank=0
       do irank=1,dtset%nimage
         if (dtset%dynimage(irank)==0.or.dtset%ntimimage<=1) then
           nrank=nrank+1
           mpi_enreg%distrb_img(irank)=mod(nrank,dtset%npimage)-1
           if (mpi_enreg%distrb_img(irank)==-1) mpi_enreg%distrb_img(irank)=dtset%npimage-1
         end if
       end do
     else
       mpi_enreg%comm_img=xmpi_comm_null
       mpi_enreg%nproc_img=0
       mpi_enreg%me_img=-1
       ABI_MALLOC(mpi_enreg%distrb_img,(0))
     end if
   end if

!  if (debug) then
!  write(200+mpi_enreg%me,*) "=================================="
!  write(200+mpi_enreg%me,*) "DEBUGGING STATMENTS IN INITMPI_IMG"
!  write(200+mpi_enreg%me,*) "=================================="
!  write(200+mpi_enreg%me,*) "option         =",option
!  write(200+mpi_enreg%me,*) "MPI_UNDEFINED  =",MPI_UNDEFINED
!  write(200+mpi_enreg%me,*) "MPI_IDENT      =",MPI_IDENT
!  write(200+mpi_enreg%me,*) "MPI_CONGRUENT  =",MPI_CONGRUENT
!  write(200+mpi_enreg%me,*) "MPI_SIMILAR    =",MPI_SIMILAR
!  write(200+mpi_enreg%me,*) "MPI_UNEQUAL    =",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "null_comm      =",MPI_COMM_NULL
!  write(200+mpi_enreg%me,*) "self_comm      =",xmpi_comm_self
!  write(200+mpi_enreg%me,*) "world_comm     =",mpi_enreg%comm_world
!  write(200+mpi_enreg%me,*) "empty_group    =",MPI_GROUP_EMPTY
!  write(200+mpi_enreg%me,*) "nimage         =",mpi_enreg%my_nimage
!  write(200+mpi_enreg%me,*) "nproc_per_image=",nproc_per_image
!  call MPI_COMM_SIZE(mpi_enreg%comm_world,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of world_comm    =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_world,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in world_comm =",irank
!  if (option==1.or.option==-1) then
!  write(200+mpi_enreg%me,*) "index_img=",mpi_enreg%my_imgtab(:)
!  end if
!  if (option==2.or.option==3.or.option==-1) then
!  write(200+mpi_enreg%me,*) "nproc_cell  =",mpi_enreg%nproc_cell
!  write(200+mpi_enreg%me,*) "me_cell     =",mpi_enreg%me_cell
!  call xmpi_comm_group(mpi_enreg%comm_cell,group_cell,ierr)
!  write(200+mpi_enreg%me,*) "group_cell  =",group_cell
!  write(200+mpi_enreg%me,*) "comm_cell   =",mpi_enreg%comm_cell
!  if (group_cell/=MPI_GROUP_EMPTY) then
!  call MPI_GROUP_SIZE(group_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of group_cell   =",irank
!  call MPI_GROUP_RANK(group_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in group_cell=",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of group_cell   =",0
!  write(200+mpi_enreg%me,*) "My rank in group_cell=",-1
!  end if
!  if (mpi_enreg%comm_cell/=MPI_COMM_NULL) then
!  call MPI_COMM_SIZE(mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of comm_cell   =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in comm_cell=",irank
!  call MPI_COMM_COMPARE(mpi_enreg%comm_world,mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison world_comm/comm_cell=",irank
!  call MPI_COMM_COMPARE(xmpi_comm_self,mpi_enreg%comm_cell,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison self_comm/comm_cell =",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of comm_cell   =",0
!  write(200+mpi_enreg%me,*) "My rank in comm_cell=",-1
!  write(200+mpi_enreg%me,*) "Comparison world_comm/comm_cell=",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "Comparison self_comm/comm_cell =",MPI_UNEQUAL
!  end if
!  end if
!  if (option==3.or.option==-1) then
!  write(200+mpi_enreg%me,*) "nproc_img  =",mpi_enreg%nproc_img
!  write(200+mpi_enreg%me,*) "me_img     =",mpi_enreg%me_img
!  write(200+mpi_enreg%me,*) "img_comm   =",mpi_enreg%comm_img
!  if (mpi_enreg%comm_img/=MPI_COMM_NULL) then
!  call MPI_COMM_SIZE(mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Size of img_comm   =",irank
!  call MPI_COMM_RANK(mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "My rank in img_comm=",irank
!  call MPI_COMM_COMPARE(mpi_enreg%comm_world,mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison world_comm/img_comm=",irank
!  call MPI_COMM_COMPARE(xmpi_comm_self,mpi_enreg%comm_img,irank,ierr)
!  write(200+mpi_enreg%me,*) "Comparison self_comm/img_comm =",irank
!  else
!  write(200+mpi_enreg%me,*) "Size of img_comm   =",0
!  write(200+mpi_enreg%me,*) "My rank in img_comm=",-1
!  write(200+mpi_enreg%me,*) "Comparison world_comm/img_comm=",MPI_UNEQUAL
!  write(200+mpi_enreg%me,*) "Comparison self_comm/img_comm =",MPI_UNEQUAL
!  end if
!  write(200+mpi_enreg%me,*) "distrb_img=",mpi_enreg%distrb_img(:)
!  end if
!  write(200+mpi_enreg%me,*)
!  call flush_unit(200+mpi_enreg%me)
!  if (option==-1) stop
!  end if

 else

!  Do not activate flag for parallelization over images
   mpi_enreg%paral_img=0
!  # of images treated by current proc
   if (dtset%optdriver==RUNL_GSTATE) then
     mpi_enreg%my_nimage=dtset%nimage
   else
     mpi_enreg%my_nimage=1
   end if
!  Indexes of images treated by current proc
   if (option==1.or.option==-1) then
     ABI_MALLOC(mpi_enreg%my_imgtab,(mpi_enreg%my_nimage))
     mpi_enreg%my_imgtab=(/(irank,irank=1,mpi_enreg%my_nimage)/)
   end if
!  Communicator over all images
   if (option==2.or.option==3.or.option==-1) then
!    Communicator for one image
     mpi_enreg%nproc_cell=mpi_enreg%nproc
     mpi_enreg%me_cell=mpi_enreg%me
   end if
   if (option==3.or.option==-1) then
!    Communicator over all images
     mpi_enreg%nproc_img=1
     mpi_enreg%comm_img=xmpi_comm_self
     mpi_enreg%me_img=0
     ABI_MALLOC(mpi_enreg%distrb_img,(dtset%nimage))
     mpi_enreg%distrb_img(:)=0
   end if
 end if

 DBG_EXIT("COLL")

end subroutine initmpi_img
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_img
!! NAME
!!  clnmpi_img
!!
!! FUNCTION
!!  Cleans-up the mpi information for parallelism over images of the cell (npimage>1).
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_img(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type), intent(inout) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

 if (xmpi_paral==0) return

 if (mpi_enreg%comm_cell/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_cell)
   mpi_enreg%comm_cell=xmpi_comm_null
 end if

 if (mpi_enreg%comm_img/=mpi_enreg%comm_world) then
   call xmpi_comm_free(mpi_enreg%comm_img)
   mpi_enreg%comm_img=xmpi_comm_null
 end if

 if (allocated(mpi_enreg%my_imgtab))  then
   ABI_FREE(mpi_enreg%my_imgtab)
 end if
 if (allocated(mpi_enreg%distrb_img))  then
   ABI_FREE(mpi_enreg%distrb_img)
 end if

 mpi_enreg%paral_img=0
 mpi_enreg%my_nimage=1
 mpi_enreg%me_img=0
 mpi_enreg%me_cell=0
 mpi_enreg%nproc_img=1
 mpi_enreg%nproc_cell=1

 DBG_EXIT("COLL")

end subroutine clnmpi_img
!!***

!!****f* ABINIT/initmpi_pert
!! NAME
!!  initmpi_pert
!!
!! FUNCTION
!!  Creates group for Parallelization over Perturbations.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_pert(dtset,mpi_enreg)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
 integer:: iprocmin,irank,npert,nproc_per_cell,nrank,numproc
 integer,allocatable :: ranks(:)
 !character(len=500) :: msg
!arrays
 integer,pointer :: nkpt_rbz(:)
 real(dp),pointer :: nband_rbz(:,:)

! ***********************************************************************

 if (mpi_enreg%me_pert<0) then
   ABI_ERROR('Error in MPI distribution! Change your proc(s) distribution or use autoparal>0.')
 end if

 call dtset%get_npert_rbz(nband_rbz, nkpt_rbz, npert)

 if (dtset%nppert>=1) then
   if (mpi_enreg%comm_cell/=mpi_enreg%comm_world) then
     call xmpi_comm_free(mpi_enreg%comm_cell)
   end if
   mpi_enreg%comm_cell=mpi_enreg%comm_world

   ! These values will be properly set in set_pert_comm
   mpi_enreg%me_cell=mpi_enreg%me
   mpi_enreg%nproc_cell=mpi_enreg%nproc

   if (mpi_enreg%me>=0) then
     nproc_per_cell=mpi_enreg%nproc/dtset%nppert
     ABI_MALLOC(ranks,(dtset%nppert))
     iprocmin=mod(mpi_enreg%me,nproc_per_cell)
     ranks=(/((iprocmin+(irank-1)*nproc_per_cell),irank=1,dtset%nppert)/)
     mpi_enreg%comm_pert=xmpi_subcomm(mpi_enreg%comm_world,dtset%nppert,ranks)
     ABI_FREE(ranks)
     mpi_enreg%me_pert=xmpi_comm_rank(mpi_enreg%comm_pert)
     mpi_enreg%nproc_pert=dtset%nppert
     if (iprocmin==0.and.mpi_enreg%me_pert==0.and.mpi_enreg%me/=0) then
       ABI_BUG('Error on me_pert!')
     end if
!    Define mpi_enreg%distrb_pert
     ABI_MALLOC(mpi_enreg%distrb_pert,(npert))
     nrank=0
     do irank=1,npert
       nrank=nrank+1
       mpi_enreg%distrb_pert(irank)=mod(nrank,dtset%nppert)-1
       if (mpi_enreg%distrb_pert(irank)==-1) mpi_enreg%distrb_pert(irank)=dtset%nppert-1
     end do
     ! Make sure that subrank 0 is working on the last perturbation
     ! Swap the ranks if necessary
     numproc=mpi_enreg%distrb_pert(npert)
     if(numproc/=0) then
       do irank=1,npert
         if (mpi_enreg%distrb_pert(irank)==numproc) mpi_enreg%distrb_pert(irank)=-2
         if (mpi_enreg%distrb_pert(irank)==0) mpi_enreg%distrb_pert(irank)=-3
       end do
       do irank=1,npert
         if (mpi_enreg%distrb_pert(irank)==-2) mpi_enreg%distrb_pert(irank)=0
         if (mpi_enreg%distrb_pert(irank)==-3) mpi_enreg%distrb_pert(irank)=numproc
       end do
     end if
!    Communicator over one cell
     ABI_MALLOC(ranks,(nproc_per_cell))
     iprocmin=(mpi_enreg%me/nproc_per_cell)*nproc_per_cell
     ranks=(/((iprocmin+irank-1),irank=1,nproc_per_cell)/)
     mpi_enreg%comm_cell_pert=xmpi_subcomm(mpi_enreg%comm_world,nproc_per_cell,ranks)
     ABI_FREE(ranks)
   end if

 else  !nppert<=1
   mpi_enreg%nproc_pert=1
   mpi_enreg%comm_pert=xmpi_comm_self
   mpi_enreg%me_pert=0
   ABI_MALLOC(mpi_enreg%distrb_pert,(npert))
   mpi_enreg%distrb_pert(:)=0
 end if

 ABI_FREE(nband_rbz)
 ABI_FREE(nkpt_rbz)

end subroutine initmpi_pert
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_pert
!! NAME
!!  clnmpi_pert
!!
!! FUNCTION
!!  Cleans-up the mpi information for parallelization over perturbations.
!!
!! INPUTS
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_pert(mpi_enreg)

!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg

! ***********************************************************************

 DBG_ENTER("COLL")

 if (xmpi_paral==0) return

 if(mpi_enreg%paral_pert == 1) then

   !  Reset communicators
   if (mpi_enreg%comm_pert/=mpi_enreg%comm_world) then
     call xmpi_comm_free(mpi_enreg%comm_pert)
     mpi_enreg%comm_pert=xmpi_comm_null
   end if

   if (allocated(mpi_enreg%distrb_pert))  then
     ABI_FREE(mpi_enreg%distrb_pert)
   end if

   mpi_enreg%me_pert=0
   mpi_enreg%me_cell=0
   mpi_enreg%nproc_pert=1
   mpi_enreg%nproc_cell=1
 end if

 DBG_EXIT("COLL")

end subroutine clnmpi_pert
!!***

!!****f* ABINIT/initmpi_band
!! NAME
!!  initmpi_band
!!
!! FUNCTION
!!  Initializes the mpi information for band parallelism (paralbd=1).
!!
!! INPUTS
!!  mpi_enreg= information about MPI parallelization
!!  nband(nkpt*nsppol)= number of bands per k point, for each spin
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!  mpi_enreg=information about MPI parallelization
!!  mpi_enreg%comm_band=communicator of BAND set
!!
!! PARENTS
!!      m_dfpt_looppert,m_dfpt_lw
!!
!! CHILDREN
!!
!! SOURCE

subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,iproc_min,iproc_max,irank,isppol
 integer :: me,nband_k,nproc,nbsteps,nrank,nstates,spacecomm
 character(len=500) :: msg
!arrays
 integer,allocatable :: ranks(:)

! ***********************************************************************

 mpi_enreg%comm_band=xmpi_comm_self

 if (mpi_enreg%paralbd==1.and.xmpi_paral==1) then

!  Comm_kpt is supposed to treat spins, k-points and bands
   spacecomm=mpi_enreg%comm_kpt
   nproc=mpi_enreg%nproc_kpt
   me=mpi_enreg%me_kpt

   nstates=sum(nband(1:nkpt*nsppol))
   nbsteps=nstates/nproc
   if (mod(nstates,nproc)/=0) nbsteps=nbsteps+1

   if (nbsteps<maxval(nband(1:nkpt*nsppol))) then

     nrank=0
     do isppol=1,nsppol
       do ikpt=1,nkpt
         ii=ikpt+(isppol-1)*nkpt
         nband_k=nband(ii)
         if (nbsteps<nband_k) then
           iproc_min=minval(mpi_enreg%proc_distrb(ikpt,:,isppol))
           iproc_max=maxval(mpi_enreg%proc_distrb(ikpt,:,isppol))
           if ((me>=iproc_min).and.(me<=iproc_max)) then
             nrank=iproc_max-iproc_min+1
             if (.not.allocated(ranks)) then
               ABI_MALLOC(ranks,(nrank))
               if (nrank>0) ranks=(/((iproc_min+irank-1),irank=1,nrank)/)
             else if (nrank/=size(ranks)) then
               msg='Number of bands per proc should be the same for all k-points!'
               ABI_BUG(msg)
             end if
           end if
         end if
       end do
     end do
     if (.not.allocated(ranks)) then
       ABI_MALLOC(ranks,(0))
     end if

     mpi_enreg%comm_band=xmpi_subcomm(spacecomm,nrank,ranks)

     ABI_FREE(ranks)
   end if
 end if

end subroutine initmpi_band
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/pre_gather
!!
!! NAME
!!  pre_gather
!!
!! FUNCTION
!!  Gathers data from FFT processors.
!!
!! INPUTS
!!  n1,n2,n3= FFT grid dimensions
!!  n4= n3/mpi_enreg%nproc_fft
!!  array= data to gather among procs
!!
!! OUTPUT
!!  None
!!
!! SIDE EFFECTS
!!  array_allgather= gathered data
!!
!! PARENTS
!!      m_forces
!!
!! CHILDREN
!!
!! SOURCE


subroutine pre_gather(array,array_allgather,n1,n2,n3,n4,mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,n3,n4
 real(dp),intent(in) :: array(n1,n2,n4,1)
 real(dp),intent(inout) :: array_allgather(n1,n2,n3,1)
 type(mpi_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
 integer :: ier

! *********************************************************************

!Gather the array on all procs
 call xmpi_allgather(array,n1*n2*n3/mpi_enreg%nproc_fft,array_allgather,mpi_enreg%comm_fft,ier)

end subroutine pre_gather
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/pre_scatter
!!
!! NAME
!!  pre_scatter
!!
!! FUNCTION
!!  Scatters data among FFT processors.
!!
!! INPUTS
!!  n1,n2,n3= FFT grid dimensions
!!  n4= n3/mpi_enreg%nproc_fft
!!  array_allgather= data to scatter among FFT procs
!!
!! OUTPUT
!!  array= scattered data
!!
!! PARENTS
!!      m_forces
!!
!! CHILDREN
!!
!! SOURCE

subroutine pre_scatter(array,array_allgather,n1,n2,n3,n4,mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: n1,n2,n3,n4
 real(dp),intent(out) :: array(n1,n2,n4,1)
 real(dp),intent(in) :: array_allgather(n1,n2,n3,1)
 type(mpi_type),intent(in) :: mpi_enreg

! *********************************************************************

!Perform the reverse operation
 array(:,:,:,:) = &
&  array_allgather(:,:,n3/mpi_enreg%nproc_fft*mpi_enreg%me_fft+1:n3/mpi_enreg%nproc_fft*(mpi_enreg%me_fft+1),:)

end subroutine pre_scatter
!!***

!!****f* m_mpinfo/iwrite_fftdatar
!! NAME
!!  iwrite_fftdatar
!!
!! FUNCTION
!!  This function selects the subset of processors that should write density/potential
!!  Return True if the processors should do IO.
!!
!! INPUTS
!!  mpi_enreg<MPI_type>=Datatype gathering information on the parallelism
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

logical function iwrite_fftdatar(mpi_enreg) result(ans)

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg

! *********************************************************************

 ans = (xmpi_paral==0 .or. &                                  ! No MPI
  (mpi_enreg%paral_kgb==0 .and. mpi_enreg%me==0) .or. &       ! paral_kgb=0 does not use MPI-FFT and cartesian communicators.
  (mpi_enreg%paral_kgb==1 .and. mpi_enreg%me_band==0 .and. &  ! select procs in one FFT communicator.
  mpi_enreg%me_kpt==0 .and. mpi_enreg%me_spinor==0) .or.   &
  (mpi_enreg%paral_pert==1 .and. mpi_enreg%me_cell==0) & ! Group master in perturbation communicator.
  )

end function iwrite_fftdatar
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/distrb2
!! NAME
!!  distrb2
!!
!! FUNCTION
!!  Creates the tabs of repartition of processors for sharing the jobs on k-points, spins and bands.
!!
!! INPUTS
!!  mband = maximum number of bands
!!  nband(nkpt*nsppol) = number of bands per k point, for each spin
!!  nkpt = number of k-points
!!  nproc= number of processors available for this distribution
!!  nsppol = 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!   mpi_enreg%proc_distrb(nkpt,mband,nsppol)=number of the processor
!!       that will treat each band in each k point.
!!   mpi_enreg%nproc_kpt is set
!!
!! NOTES
!!  For the time being, the band parallelisation works only
!!  when the number of bands is identical for spin up and spin down
!!  at the same k point. The problem is the most clearly seen
!!  in the kpgio routine, where a different parallel repartition
!!  of k points for spin up and spin down would conflict with the
!!  present computation of k+G sphere, independent of the spin.
!!
!! PARENTS
!!      m_dfpt_looppert,m_dfpt_lw,m_eig2d,m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine distrb2(mband,nband,nkpt,nproc,nsppol,mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: mband,nkpt,nproc,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: inb,inb1,ind,ind0,nband_k,proc_max,proc_min
 integer :: iiband,iikpt,iisppol,ikpt_this_proc,nbsteps,nproc_kpt,temp_unit
 integer :: kpt_distrb(nkpt)
 logical,save :: first=.true.,has_file
 character(len=500) :: msg

!******************************************************************

 nproc_kpt=mpi_enreg%nproc_kpt
 if (mpi_enreg%paral_pert==1) nproc_kpt=nproc

!Initialization of proc_distrb
 do iisppol=1,nsppol
   do iiband=1,mband
     do iikpt=1,nkpt
       mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=nproc_kpt-1
     end do
   end do
 end do
!That s all for an empty communication space
 if (nproc==0) return

!Some checks
 !if (mpi_enreg%paralbd==0 .and. any(dtset%optdriver == [RUNL_GSTATE, RUNL_RESPFN])) then
 if (mpi_enreg%paralbd==0) then
!  Check if nkpt and nproc_kpt match
   if(nproc_kpt>nkpt*nsppol) then
!    Too many proc. with respect to nkpt
     write(msg,'(a,i0,a,i0,a,i0,2a)')&
      'nproc_kpt= ',nproc_kpt,' >= nkpt= ',nkpt,'* nsppol= ',nsppol,ch10,&
      'The number of processors is larger than nkpt*nsppol. This is a waste. (Ignore this warning if this is not a GS run)'
     ABI_WARNING(msg)
   else if (mod(nkpt*nsppol,nproc_kpt) /= 0) then
!    nkpt not a multiple of nproc_kpt
     write(msg,'(a,i0,a,i0,3a)')&
      'nkpt*nsppol (', nkpt*nsppol, ') is not a multiple of nproc_kpt (',nproc_kpt, ')', ch10,&
      'The k-point parallelisation is inefficient. (Ignore this warning if this is not a GS run)'
     ABI_WARNING(msg)
   end if
 end if

!Inquire whether there exist a file containing the processor distribution
 if (first) then
!  Case first time: test file to do
!  Open the file containing the k-point distribution
   first=.false.; has_file = file_exists("kpt_distrb")
 end if

!Initialize the processor distribution, either from a file, or from an algorithm
 if (has_file) then
   if (open_file('kpt_distrb',msg,newunit=temp_unit,form='formatted',status='old') /= 0) then
     ABI_ERROR(msg)
   end if
   rewind(unit=temp_unit)
   if (mpi_enreg%paralbd == 1) then
!    -> read bands distribution
     read(temp_unit,*) mpi_enreg%proc_distrb
   else
     read(temp_unit,*) kpt_distrb
   end if
   close(temp_unit)
   proc_max=0
   proc_min=nproc_kpt
!  -> determine the range of proc. requested
   if (mpi_enreg%paralbd == 1) then
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         proc_max=maxval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
         proc_min=minval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
       end do
     end do
   else
     proc_max=maxval(kpt_distrb(1:nkpt))
     proc_min=minval(kpt_distrb(1:nkpt))
!    -> fill the tab proc_distrb with kpt_distrb
     do iisppol=1,nsppol
       do iikpt=1,nkpt
         nband_k = nband(iikpt+(iisppol-1)*nkpt)
         do iiband=1,nband_k
           mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=kpt_distrb(iikpt)
         end do
       end do
     end do
   end if ! mpi_enreg%paralbd

   if(proc_max>(nproc_kpt-1)) then
!    Too much proc. requested
     write(msg, '(a,a,a,i0,a,a,a)' )&
      'The number of processors mentioned in the kpt_distrb file',ch10,&
      'must be lower or equal to the actual number of processors =',nproc_kpt-1,ch10,&
      'Action: change the kpt_distrb file, or increase the','  number of processors.'
     ABI_ERROR(msg)
   end if

   if(proc_max/=(nproc_kpt-1)) then
!    Too few proc. used
     write(msg, '(a,i0,a,a,a,i0,a,a,a)' )&
      'Only ',proc_max+1,' processors are used (from kpt_distrb file),',ch10,&
      'when',nproc_kpt,' processors are available.',ch10,&
      'Action: adjust number of processors and kpt_distrb file.'
     ABI_ERROR(msg)
   end if

   if(proc_min<0) then
     write(msg, '(a,a,a)' )&
      'The number of processors must be bigger than 0 in kpt_distrb file.',ch10,&
      'Action: modify kpt_distrb file.'
     ABI_ERROR(msg)
   end if

 else
!  'kpt_distrb' file does not exist

   if (mpi_enreg%paralbd==1) then

!    No possible band parallelization
     if (nproc<(nkpt*nsppol)) then

!      Does not allow a processor to treat different spins
       ind0=0
       inb1=(nkpt*nsppol)/nproc;if (mod((nkpt*nsppol),nproc)/=0) inb1=inb1+1
       do iikpt=1,nkpt
         nband_k=nband(iikpt)
         ind=ind0/inb1
         do iiband=1,nband_k
           mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
           if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc-ind-1
         end do
         ind0=ind0+1
       end do

!      MT130831 : OLD CODING
!      do iisppol=1,nsppol;do iikpt=1,nkpt
!      ind=(iikpt+(iisppol-1)*nkpt-1)/inb1
!      nband_k=nband(iikpt+(iisppol-1)*nkpt)
!      do iiband=1,nband_k
!      mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind
!      end do;end do;end do
!      MT130831 : END OF OLD CODING

!    Possible band parallelization
     else
!      Does not allow a processor to treat different spins
       ind0=0;inb=nproc/(nkpt*nsppol)
       do iikpt=1,nkpt
         nband_k=nband(iikpt)
         inb1=nband_k/inb;if (mod(nband_k,inb)/=0) inb1=inb1+1
         do iiband=1,nband_k
           ind=(iiband-1)/inb1+ind0
           mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
           if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc-ind-1
         end do
         ind0=ind+1
       end do

!      MT130831 : OLD CODING
!      ind0=0;inb=nproc/(nkpt*nsppol)
!      do iisppol=1,nsppol;do iikpt=1,nkpt
!      nband_k=nband(iikpt+(iisppol-1)*nkpt)
!      inb1=nband_k/inb;if (mod(nband_k,inb)/=0) inb1=inb1+1
!      do iiband=1,nband_k
!      ind=(iiband-1)/inb1+ind0
!      mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind
!      end do
!      ind0=ind+1
!      end do;end do
!      MT130831 : END OF OLD CODING

     end if

!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol;do iikpt=1,nkpt
!    nband_k=nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    ind=ind+1
!    end do;end do;end do
!    XG060807 : END OF OLD CODING

   elseif (mpi_enreg%paralbd==0) then

!    Does not allow a processor to treat different spins
     ind0=0
     nbsteps=(nsppol*nkpt)/nproc_kpt
     if (mod((nsppol*nkpt),nproc_kpt)/=0) nbsteps=nbsteps+1
     do iikpt=1,nkpt
       nband_k=nband(iikpt)
       ind=ind0/nbsteps
       do iiband=1,nband_k
         mpi_enreg%proc_distrb(iikpt,iiband,1)=ind
         if (nsppol==2) mpi_enreg%proc_distrb(iikpt,iiband,2)=nproc_kpt-ind-1
       end do
       ind0=ind0+1
     end do

!    XG060807 : OLD CODING
!    ind=0
!    do iisppol=1,nsppol;do iikpt=1,nkpt
!    nband_k = nband(iikpt+(iisppol-1)*nkpt)
!    do iiband=1,nband_k
!    Distribute k-points homogeneously
!    proc_distrb(iikpt,iiband,iisppol)=mod(iikpt-1,nproc_kpt)
!    mpi_enreg%proc_distrb(iikpt,iiband,iisppol)=ind/nbsteps
!    end do
!    ind=ind + 1
!    end do;end do
!    XG060807 : END OF OLD CODING

   end if ! mpi_enreg%paralbd

 end if ! has_file

 mpi_enreg%my_kpttab(:)=0
 mpi_enreg%my_isppoltab(:)=0
 do iisppol=1,nsppol
   ikpt_this_proc=0
   do iikpt=1,nkpt
     nband_k=nband(iikpt+(iisppol-1)*nkpt)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,iikpt,1,nband_k,iisppol,mpi_enreg%me_kpt)) cycle
     ikpt_this_proc=ikpt_this_proc+1
     mpi_enreg%my_kpttab(iikpt)=ikpt_this_proc
     mpi_enreg%my_isppoltab(iisppol)=1
   end do
 end do

end subroutine distrb2
!!***

!!****f* ABINIT/distrb2_hf
!! NAME
!!  distrb2_hf
!!
!! FUNCTION
!!  Ceate the tabs of repartition of processors for sharing the jobs
!!  on occupied states (labeled by k-points, bands and spin indices) for Hartree-Fock calculation.
!!
!! INPUTS
!!  nbandhf = maximum number of occupied bands
!!  nkpthf = number of k-points in full BZ
!!  nproc= number of processors available for this distribution
!!  nsppol = 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!   mpi_enreg%proc_distrb(nkpthf,nbandhf,nsppol)=number of the processor
!!       that will treat each band in each k point.
!!   mpi_enreg%nproc_kpt is set
!!
!! NOTES
!!  For the time being, the band parallelisation works only
!!  when the number of bands is identical for spin up and spin down
!!  at the same k point. The problem is the most clearly seen
!!  in the kpgio routine, where a different parallel repartition
!!  of k points for spin up and spin down would conflict with the
!!  present computation of k+G sphere, independent of the spin.
!!
!! PARENTS
!!      m_mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine distrb2_hf(nbandhf,nkpthf, nproc, nsppol, mpi_enreg)

!Arguments ------------------------------------
 integer,intent(in) :: nbandhf,nkpthf,nproc,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: ind,iiband,iikpt,iistep,nproc_hf
 character(len=500) :: msg

!******************************************************************

 nproc_hf=mpi_enreg%nproc_hf

!* Initialize distrb_hf (the array exists necessarily)
 do iiband=1,nbandhf
   do iikpt=1,nkpthf
     mpi_enreg%distrb_hf(iikpt,iiband,1)=nproc_hf-1
   end do
 end do

!* End of the routine for an empty communication space
 if (nproc==0) return

!*** Testing section ***

 if (nsppol==2) then
!* Check that the distribution over (spin,k point) allow to consider spin up and spin dn independently.
   if (mpi_enreg%nproc_kpt/=1.and.mod(mpi_enreg%nproc_kpt,2)/=0) then
     ABI_ERROR( 'The variable nproc_kpt is not even but nssppol= 2')
!* In this case, one processor will carry both spin. (this will create pbms for the calculation)
   end if
!* Check that the number of band is the same for each spin, at each k-point. (It should be)
!*   do iikpt=1,nkpthf
!*     if (nband(iikpt)/=nband(iikpt+nkpthf)) then
!*     msg = ' WARNING - the number of bands is different for spin up or spin down. '
!*     ABI_ERROR(msg)
!*     end if
!*    end do
!* If one of this test is not good for one proc then other procs fall in deadlock, according to distrb2.
!* What does it mean ???
 end if


!* Check if nkpthf and nproc_hf match
 if (nproc_hf>nkpthf*nbandhf) then
!* There are too many processors with respect to nkpthf*nbandhf
   write(msg, '(a,a,i4,a,i4,a,i4,a,a)' ) ch10,&
&   'nproc_hf=',nproc_hf,' >= nkpthf=',nkpthf,'* nbandhf=',nbandhf,ch10,&
&   'The number of processors is larger than nkpthf*nbandhf. This is a waste.'
   ABI_WARNING(msg)

 else if(mod(nkpthf*nbandhf,nproc_hf)/=0) then
!* nkpthf*nbandhf is not a multiple of nproc_hf
   write(msg, '(2a,i5,a,i5,3a)' ) ch10,&
&   'nkpthf*nbandhf (', nkpthf*nbandhf, ') is not a multiple of nproc_hf (',nproc_hf, ')', ch10,&
&   'The parallelisation may not be efficient.'
   ABI_WARNING(msg)
 end if

!*** End of testing section ***

!*** Initialize the processor distribution from a simple algorithm ***

 if (nproc_hf<nkpthf) then
!* In this case, a parallelization over kpts only.
   iistep=nkpthf/nproc_hf
   if (mod(nkpthf,nproc_hf) /=0) iistep=iistep+1
   ind=0
   do iikpt=1,nkpthf
!*** Only the first "nbandhf" bands are considered (they are assumed to be the only occupied ones)
     do iiband=1,nbandhf
       mpi_enreg%distrb_hf(iikpt,iiband,1)=ind/iistep
     end do
     ind=ind+1
   end do

 else
!* In this case, a parallelization over all the occupied states is possible.
   if (nproc_hf < nbandhf*nkpthf) then
     iistep=(nbandhf*nkpthf)/nproc_hf;
     if (mod((nbandhf*nkpthf),nproc_hf) /=0) iistep=iistep+1
   else
     iistep=1
   end if
   ind=0
   do iikpt=1,nkpthf
!*** Only the first "nbandhf" bands are considered (they are assumed to be the only occupied ones)
     do iiband=1,nbandhf
       mpi_enreg%distrb_hf(iikpt,iiband,1)=ind/iistep
       ind=ind+1
     end do
   end do
 end if

!*** Initialization of processor distribution from a file (simple copy from distrb2, not yet implemented) ***

! !* Inquire whether there exist a file containing the processor distribution
!  if (first) then
! !  Case first time : test file to do
! !  Open the file containing the (k-points,bands) distribution
!    open(unit=temp_unit,file='kpt_distrb_hf',form='formatted',status='old',iostat=ios)
!    if(ios==0) then
! !    'kpt_distrb_hf' file exists
!      file_exist=1
!      close(temp_unit)
!    else
!      file_exist=0
!    end if
!    first=.false.
!  end if
!
! !* Initialize the processor distribution, either from a file, or from an algorithm
!  if (file_exist == 1) then
! !* Read (k-points,bands) distribution out of the file
!    open(unit=temp_unit,file='kpt_distrb_hf',form='formatted',status='old',iostat=ios)
!    rewind(unit=temp_unit)
!    read(temp_unit,*) mpi_enreg%distrb_hf
!    close(temp_unit)
! !* Determine the range of processors requested
!    proc_max=0
!    proc_min=nproc_hf
!    do iikpt=1,nkpthf
!      mband_occ_k = mband_occ(iikpt+(iisppol-1)*nkpthf)
!      proc_max=maxval(mpi_enreg%distrb_hf(iikpt,1:mband_occ_k,1))
!      proc_min=minval(mpi_enreg%distrb_hf(iikpt,1:mband_occ_k,1))
!    end do
!
!    if(proc_max>(nproc_hf-1)) then
! !*    Too much proc. requested
!      write(msg, '(a,a,a,i4,a,a,a)' )&
! &     '  The number of processors mentioned in the kpt_distrb file',ch10,&
! &     '  must be lower or equal to the actual number of processors =',&
! &     nproc_hf-1,ch10,&
! &     '  Action: change the kpt_distrb file, or increase the',&
! &     '  number of processors.'
!      ABI_ERROR(msg)
!    end if
!
!    if(proc_max/=(nproc_hf-1)) then
! !*    Too few proc. used
!      write(msg, '(a,i4,a,a,a,i4,a,a,a)' )&
! &     '  Only ',proc_max+1,' processors are used (from kpt_distrb file),',ch10,&
! &     '  when',nproc_hf,' processors are available.',ch10,&
! &     '  Action: adjust number of processors and kpt_distrb file.'
!      ABI_ERROR(msg)
!    end if
!
!    if(proc_min<0) then
!      write(msg, '(a,a,a)' )&
! &     '  The number of processors must be bigger than 0 in kpt_distrb file.',ch10,&
! &     ' Action: modify kpt_distrb file.'
!      ABI_ERROR(msg)
!    end if
!  else
! !* The file does not exist...
!  end if ! file_exist

!DEBUG
!write(std_out,*)' distrb2_hf: exit '
!ENDDEBUG

end subroutine distrb2_hf
!!***

end module m_mpinfo
!!***
