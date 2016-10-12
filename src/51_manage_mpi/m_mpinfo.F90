!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mpinfo
!! NAME
!! m_mpinfo
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MT, GG)
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
 use m_profiling_abi
 use m_xmpi
 use m_distribfft

 use defs_abitypes,   only : MPI_type

 implicit none

 private 

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


! Destructor methods
 public :: clnmpi_atom
 public :: clnmpi_grid
 public :: clnmpi_img
 public :: clnmpi_pert

! Helper functions.
 public :: pre_gather
 public :: pre_scatter
 public :: iwrite_fftdatar     ! Select the subset of processors that will write density/potential files.
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
!!      lapackprof,mpi_setup
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_mpi_enreg(mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_mpi_enreg'
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

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
!!      initmpi_seq,m_fft_prof,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine nullify_mpi_enreg(MPI_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nullify_mpi_enreg'
!End of the abilint section

 implicit none

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
!!      abinit,bethe_salpeter,bsepostproc,calc_vhxc_me,conducti,cut3d
!!      debug_tools,dfpt_nstpaw,dieltcel,eph,fftprof,gwls_hamiltonian,inwffil
!!      ks_ddiago,lapackprof,linear_optics_paw,m_cut3d,m_ddk,m_dvdb,m_fft
!!      m_fft_prof,m_fftcore,m_gsphere,m_hamiltonian,m_io_kss,m_ioarr,m_kxc
!!      m_pawpwij,m_ppmodel,m_screening,m_wfd,m_wfk,mlwfovlp_qp,mrgddb,mrggkk
!!      mrgscr,posdoppler,scfcv,screening,sigma,suscep_stat,susk,suskmm,ujdet
!!      vdw_kernelgen,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_mpi_enreg(MPI_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_mpi_enreg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: MPI_enreg

! *********************************************************************

 if (associated(mpi_enreg%distribfft)) then
   call destroy_distribfft(mpi_enreg%distribfft)
   ABI_DATATYPE_DEALLOCATE(mpi_enreg%distribfft)
   nullify(mpi_enreg%distribfft)
 end if

 if (allocated(mpi_enreg%proc_distrb)) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if
 if (allocated(mpi_enreg%kptdstrb)) then
   ABI_DEALLOCATE(mpi_enreg%kptdstrb)
 end if
 if (allocated(mpi_enreg%kpt_loc2fbz_sp)) then
   ABI_DEALLOCATE(mpi_enreg%kpt_loc2fbz_sp)
 end if
 if (allocated(mpi_enreg%kpt_loc2ibz_sp)) then
   ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
 end if
 if (allocated(mpi_enreg%mkmem)) then
   ABI_DEALLOCATE(mpi_enreg%mkmem)
 end if
 if (allocated(mpi_enreg%my_kpttab)) then
   ABI_DEALLOCATE(mpi_enreg%my_kpttab)
 end if
 if (associated(mpi_enreg%my_atmtab)) then
   ABI_DEALLOCATE(mpi_enreg%my_atmtab)
   nullify(mpi_enreg%my_atmtab)
 end if
 if (allocated(mpi_enreg%distrb_pert)) then
   ABI_DEALLOCATE(mpi_enreg%distrb_pert)
 end if
 if (allocated(mpi_enreg%distrb_img)) then
   ABI_DEALLOCATE(mpi_enreg%distrb_img)
 end if
 if (allocated(mpi_enreg%my_imgtab)) then
   ABI_DEALLOCATE(mpi_enreg%my_imgtab)
 end if
 if (allocated(mpi_enreg%my_kgtab)) then
   ABI_DEALLOCATE(mpi_enreg%my_kgtab)
 end if
 if (allocated(mpi_enreg%distrb_hf)) then
   ABI_DEALLOCATE(mpi_enreg%distrb_hf)
 end if

 if(allocated(mpi_enreg%my_cells))then
   ABI_DEALLOCATE(mpi_enreg%my_cells)
 end if

 if(allocated(mpi_enreg%my_index_cells))then
   ABI_DEALLOCATE(mpi_enreg%my_index_cells)
 end if


!Do not deallocate wavelet denspot distribution arrays,
!they are handled by BigDFT.

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
!!      gwls_hamiltonian,inwffil,m_fft_prof,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_mpi_enreg(MPI_enreg1,MPI_enreg2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'copy_mpi_enreg'
!End of the abilint section

 implicit none

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
!mpi_enreg2%flag_ind_kg_mpi_to_seq=mpi_enreg1%flag_ind_kg_mpi_to_seq
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

!pointers
 if (associated(mpi_enreg1%distribfft)) then
   if (.not.associated(mpi_enreg2%distribfft)) then
     ABI_DATATYPE_ALLOCATE(mpi_enreg2%distribfft,)
   end if
   call copy_distribfft(mpi_enreg1%distribfft,mpi_enreg2%distribfft)
 end if

 if (allocated(mpi_enreg1%proc_distrb)) then
   sz1=size(mpi_enreg1%proc_distrb,1)
   sz2=size(mpi_enreg1%proc_distrb,2)
   sz3=size(mpi_enreg1%proc_distrb,3)
   ABI_ALLOCATE(mpi_enreg2%proc_distrb,(sz1,sz2,sz3))
   mpi_enreg2%proc_distrb=mpi_enreg1%proc_distrb
 end if
 if (allocated(mpi_enreg1%kptdstrb)) then
   sz1=size(mpi_enreg1%kptdstrb,1)
   sz2=size(mpi_enreg1%kptdstrb,2)
   sz3=size(mpi_enreg1%kptdstrb,3)
   ABI_ALLOCATE(mpi_enreg2%kptdstrb,(sz1,sz2,sz3))
   mpi_enreg2%kptdstrb=mpi_enreg1%kptdstrb
 end if
 if (allocated(mpi_enreg1%kpt_loc2fbz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2fbz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2fbz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2fbz_sp,3)
   ABI_ALLOCATE(mpi_enreg2%kpt_loc2fbz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2fbz_sp=mpi_enreg1%kpt_loc2fbz_sp
 end if
 if (allocated(mpi_enreg1%kpt_loc2ibz_sp)) then
   sz1=size(mpi_enreg1%kpt_loc2ibz_sp,1)-1
   sz2=size(mpi_enreg1%kpt_loc2ibz_sp,2)
   sz3=size(mpi_enreg1%kpt_loc2ibz_sp,3)
   ABI_ALLOCATE(mpi_enreg2%kpt_loc2ibz_sp,(0:sz1,1:sz2,1:sz3))
   mpi_enreg2%kpt_loc2ibz_sp=mpi_enreg1%kpt_loc2ibz_sp
 end if
 if (allocated(mpi_enreg1%mkmem)) then
   ABI_ALLOCATE(mpi_enreg2%mkmem,(0:size(mpi_enreg1%mkmem,1)-1))
   mpi_enreg2%mkmem=mpi_enreg1%mkmem
 end if
 if (associated(mpi_enreg1%my_atmtab)) then
   ABI_ALLOCATE(mpi_enreg2%my_atmtab,(size(mpi_enreg1%my_atmtab)))
   mpi_enreg2%my_atmtab=mpi_enreg1%my_atmtab
 else
   nullify(mpi_enreg2%my_atmtab)
 end if
 if (allocated(mpi_enreg1%my_kgtab)) then
   sz1=size(mpi_enreg1%my_kgtab,1)
   sz2=size(mpi_enreg1%my_kgtab,2)
   ABI_ALLOCATE(mpi_enreg2%my_kgtab,(sz1,sz2))
   mpi_enreg2%my_kgtab=mpi_enreg1%my_kgtab
 end if
 if (allocated(mpi_enreg1%distrb_pert)) then
   ABI_ALLOCATE(mpi_enreg2%distrb_pert,(size(mpi_enreg1%distrb_pert)))
   mpi_enreg2%distrb_pert=mpi_enreg1%distrb_pert
 end if
 if (allocated(mpi_enreg1%distrb_img)) then
   ABI_ALLOCATE(mpi_enreg2%distrb_img,(size(mpi_enreg1%distrb_img)))
   mpi_enreg2%distrb_img=mpi_enreg1%distrb_img
 end if
 if (allocated(mpi_enreg1%my_imgtab)) then
   ABI_ALLOCATE(mpi_enreg2%my_imgtab,(size(mpi_enreg1%my_imgtab)))
   mpi_enreg2%my_imgtab=mpi_enreg1%my_imgtab
 end if
 if (allocated(mpi_enreg1%distrb_hf)) then
   sz1=size(mpi_enreg1%distrb_hf,1)
   sz2=size(mpi_enreg1%distrb_hf,2)
   sz3=size(mpi_enreg1%distrb_hf,3)
   ABI_ALLOCATE(mpi_enreg2%distrb_hf,(sz1,sz2,sz3))
   mpi_enreg2%distrb_hf=mpi_enreg1%distrb_hf
 end if

!Optional pointers
 if (allocated(mpi_enreg1%my_kpttab)) then
   ABI_ALLOCATE(mpi_enreg2%my_kpttab,(size(mpi_enreg1%my_kpttab)))
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
!!      atm2fft,dfpt_atm2fft,pawmknhat,pawmknhat_psipsi,pawsushat,posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine set_mpi_enreg_fft(MPI_enreg,comm_fft,distribfft,me_g0,paral_kgb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'set_mpi_enreg_fft'
!End of the abilint section

 implicit none

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
   ABI_DATATYPE_DEALLOCATE(mpi_enreg%distribfft)
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
!!      atm2fft,dfpt_atm2fft,pawmknhat,pawmknhat_psipsi,pawsushat,posdoppler
!!
!! CHILDREN
!!
!! SOURCE

subroutine unset_mpi_enreg_fft(MPI_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unset_mpi_enreg_fft'
!End of the abilint section

 implicit none

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
!!      dfpt_eltfrhar,dfpt_eltfrloc,dfpt_vlocal,fftpac,fourdp,hartre,hartrestr
!!      indirect_parallel_Fourier,initro,laplacian,m_fock,m_ioarr,mag_constr
!!      make_efg_el,mkcore,mkcore_paw,mklocl_realspace,mklocl_recipspace
!!      moddiel,out1dm,posdoppler,prcrskerker2,strhar,symrhg,vlocalstr,xcden
!!      xcpot
!!
!! CHILDREN
!!
!! SOURCE

subroutine ptabs_fourdp(MPI_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ptabs_fourdp'
!End of the abilint section

 implicit none

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

 if(.not.(grid_found)) then
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
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
!!      fourwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine ptabs_fourwf(MPI_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ptabs_fourwf'
!End of the abilint section

 implicit none

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

 if(.not.(grid_found)) then
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpi_distrib_is_ok'
!End of the abilint section

 implicit none

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
&        'Your number of spins*k-points (=',nsppol*nkpt,') ',&
&        'will not distribute correctly',ch10, &
&        'with the current number of processors (=',MPI_enreg%nproc_kpt,').',ch10,&
&        'You will leave some empty.'
     end if
   end if
 else
   if (mod(nband,max(1,MPI_enreg%nproc_kpt/(nsppol*nkpt)))/=0) then
     mpi_distrib_is_ok=.false.
     if (present(msg)) then
       write(msg,'(a,i0,2a,i0,4a,i0,7a)')&
&        'Your number of spins*k-points (=',nsppol*nkpt,') ',&
&         'and bands (=',nband,') ',&
&         'will not distribute correctly',ch10,&
&         'with the current number of processors (=',MPI_enreg%nproc_kpt,').',ch10,&
&         'You will leave some empty.'
     end if
   end if
 end if

end function mpi_distrib_is_ok
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_atom
!! NAME
!!  clnmpi_atom
!!
!! FUNCTION
!!  Cleans-up the mpi informations for the parallelism over atoms (PAW).
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_atom(mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_atom'
!End of the abilint section

 implicit none

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
   ABI_DEALLOCATE(mpi_enreg%my_atmtab)
 end if

 mpi_enreg%nproc_atom=1
 mpi_enreg%my_natom=0 ! should be natom

 DBG_EXIT("COLL")

end subroutine clnmpi_atom
!!***

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_grid
!! NAME
!!  clnmpi_grid
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelism over grid (kpt/band/fft).
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_grid(mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_grid'
!End of the abilint section

 implicit none

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

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_img
!! NAME
!!  clnmpi_img
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelism over images of the cell (npimage>1).
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine clnmpi_img(mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_img'
!End of the abilint section

 implicit none

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
   ABI_DEALLOCATE(mpi_enreg%my_imgtab)
 end if
 if (allocated(mpi_enreg%distrb_img))  then
   ABI_DEALLOCATE(mpi_enreg%distrb_img)
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

!----------------------------------------------------------------------

!!****f* m_mpinfo/clnmpi_pert
!! NAME
!!  clnmpi_pert
!!
!! FUNCTION
!!  Cleans-up the mpi informations for parallelization over perturbations.
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnmpi_pert'
!End of the abilint section

 implicit none

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
     ABI_DEALLOCATE(mpi_enreg%distrb_pert)
   end if

   mpi_enreg%me_pert=0
   mpi_enreg%me_cell=0
   mpi_enreg%nproc_pert=1
   mpi_enreg%nproc_cell=1
 end if

 DBG_EXIT("COLL")

end subroutine clnmpi_pert
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
!!      fresid
!!
!! CHILDREN
!!
!! SOURCE


subroutine pre_gather(array,array_allgather,n1,n2,n3,n4,mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pre_gather'
!End of the abilint section

 implicit none

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
!!      fresid
!!
!! CHILDREN
!!
!! SOURCE

subroutine pre_scatter(array,array_allgather,n1,n2,n3,n4,mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pre_scatter'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iwrite_fftdatar'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg

! *********************************************************************

 ans = (xmpi_paral==0 .or. &                                  ! No MPI
  (mpi_enreg%paral_kgb==0 .and. mpi_enreg%me==0) .or. &       ! paral_kgb=0 does not use MPI-FFT and cartesian communicators. 
  (mpi_enreg%paral_kgb==1 .and. mpi_enreg%me_band==0 .and. &  ! select procs in one FFT communicator.
  mpi_enreg%me_kpt==0 .and. mpi_enreg%me_spinor==0))

end function iwrite_fftdatar
!!***

!----------------------------------------------------------------------

end module m_mpinfo
