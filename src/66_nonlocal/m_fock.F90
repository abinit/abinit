!!****m* ABINIT/m_fock
!! NAME
!!  m_fock
!!
!! FUNCTION
!!  This module provides the definition of
!!  the fock_type used to store data for the calculation of Fock exact exchange term
!!  and the procedures to perform this calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (CMartins,FJ,FA,MT)
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

module m_fock

 use m_vcoul 

 use defs_basis
 use m_abicore
 use m_errors
 use m_mpinfo
 use m_xmpi
 use libxc_functionals
 use m_pawang
 use m_pawtab
 use m_pawfgr
 use m_pawfgrtab
 use m_pawcprj
 use m_cgtools
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_time,            only : timab
 use m_fstrings,        only : itoa, ftoa, sjoin
 use m_symtk,           only : mati3inv, matr3inv
 use m_fftcore,         only : sphereboundary
 use m_fft,             only : zerosym, fourwf
 use m_kg,              only : ph1d3d, getph
 use m_kpts,            only : listkk

 implicit none

 private
!!***

!!****t* m_fock/fock_type
!! NAME
!!  fock_type
!!
!! FUNCTION
!!   This object stores the occupied wavefunctions and other quantities
!!   needed to calculate Fock exact exchange
!!
!! SOURCE

 type, public :: fock_type
  type(fock_common_type), pointer :: fock_common=> null()
  type(fock_BZ_type), pointer :: fock_BZ=> null()
  type(fock_ACE_type), pointer :: fockACE(:,:)=> null()
 end type fock_type


 type, public :: fock_common_type

! Integer scalars
  !integer :: mcgocc_bz,mkg_bz,mocc
  !integer :: natom,ntypat

  integer :: usepaw
    ! 0 if norm-conserving psps, 1 for PAW (not implemented)

  integer :: ikpt,isppol,ieigen,iband
    ! data relative to the current states.

  integer :: mband
    ! maximum number of bands

  integer :: my_nsppol
   ! my_nsppol=1 when nsppol=1 or nsppol=2 and only one spin is treated by the processor.
   ! my_nsppol=2 when nsppol=2 and no parallelization over kpt (both spins are treated by the processor).

  integer :: natom
   ! Number of atoms, input variable

  integer :: nsppol
   ! Number of independent spin polarizations, input variable
   ! Note that this value does not take into account the MPI distribution of the wavefunctions.

  integer :: ntypat
   ! Number of type of atoms

  integer :: nnsclo_hf
    ! Number of iterations with fixed occupied states when calculating the exact exchange contribution.

  integer :: ixc
    ! XC option (abinit input variable)

  integer :: use_ACE
    ! option to use the ACE method of Lin Lin
    !==0 if the normal Fock operator is to be created and/or used
    !==1 if the ACE operator is to be created and/or used

  integer ABI_PRIVATE :: getghc_call_ = 1
  ! 1 if fock_getghc should be called in getghc, 0 otherwise

! Logical
  logical :: optfor
    ! option to calculate forces

  logical :: optstr
    ! option to calculate stresses

  logical :: fock_converged
    ! .false. if the Fock cycle (with changing Fock/ACE operator) is not converged
    ! .true. if the Fock cycle (with changing Fock/ACE operator) has converged

  logical :: scf_converged
    ! .false. if the SCF cycle (with fixed Fock/ACE operator) is not converged
    ! .true. if the SCF cycle (with fixed Fock/ACE operator) has converged

! Real(dp) scalars

  real(dp) :: gsqcut
    !  cutoff value on G**2 for sphere inside fft box.
    !   (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2)). Used in hartre

  real(dp) :: hyb_mixing
    ! hybrid mixing coefficient for the Fock contribution

  real(dp) :: hyb_mixing_sr
    ! hybrid mixing coefficient for the short-range Fock contribution

  real(dp) :: hyb_range_dft
    ! hybrid range for separation, used in the DFT functional
    ! (should be equal to hyb_range_fock, but this is not true for HSE03)

  real(dp) :: hyb_range_fock
    ! hybrid range for separation, used in the fock contribution

  real(dp) :: e_fock0
    ! contribution of the Fock term to energy (computed and stored here in case of ACE)

  integer, allocatable :: atindx(:)
    !  atindx(natom)=index table for atoms (see gstate.f)

  integer, allocatable  :: nband(:)
   ! nband(nkpt)
   ! Number of bands for each k point

  integer, allocatable :: symrec(:,:,:)

  integer,allocatable :: typat(:)
   ! typat(natom)
   ! type of each atom

! Real(dp) arrays
  real(dp) :: stress(6)
    ! stress(6)
    ! contribution of the fock term to stresses

  real(dp), allocatable  :: stress_ikpt(:,:)
    ! stress(6,nband)
    ! contribution of the fock term to stresses for the current band

  real(dp), allocatable :: forces_ikpt(:,:,:)
    ! forces(3,natom,nband))
    ! contribution of the fock term to forces for the current band

  real(dp), allocatable :: forces(:,:)
    ! forces(3,natom))
    ! contribution of the fock term to forces

  real(dp), allocatable :: eigen_ikpt(:)
    ! eigen_ikpt,(nband))
    !  Will contain the band index of the current state
    !  if the value is 0, the Fock contribution to the eigenvalue is not calculated.

! Pointers to PAW-types (associated only if usepaw==1)
! Note that these are references to already existing objects.

  type(pawtab_type), pointer :: pawtab(:)
  type(pawfgr_type),pointer :: pawfgr
  type(pawfgrtab_type),allocatable :: pawfgrtab(:)

 end type fock_common_type

 type, public :: fock_BZ_type

  integer :: mcprj
    ! dimension of cwaveocc_cprj

  integer :: mkpt
    ! maximum number of k-points for Fock treated by this node

  integer :: mkptband
    ! size of occupied states stored by this node.

  integer :: nkpt_bz
    ! Number of k-points in the BZ for Fock operator

  integer, allocatable   :: gbound_bz(:,:,:)
    ! gbound_bz(2*mgfft+8,2,mkpt)
    ! Tables for zero-padded FFT of wavefunctions.

  integer, allocatable :: kg_bz(:,:)
    ! kg_bz(3,mpw*mkpt)
    ! G-vectors for each k-point in the BZ treate by this node

  integer, allocatable :: nbandocc_bz(:,:)
    ! nbandocc_bz,(mkpt,my_nsppol))
    ! nb of bands at each k point

  integer, allocatable :: npwarr(:)
    ! npwarr(mkpt)

  integer, allocatable :: istwfk_bz(:)
    ! istwfk_bz,(mkpt))
    ! storage mode of the wavefunction at each k-point

  integer, allocatable :: calc_phase(:)
    ! calc_phase,(mkpt))
    ! 1 if a phase factor must be considered (0 otherwise) at each k point

  integer, allocatable :: tab_symkpt(:)
    ! tab_symkpt,(mkpt))
    ! indices of symmetry operation to apply to get jkpt in full BZ from ikpt in IBZ

  integer, allocatable :: timerev(:)
    ! timerev,(mkpt))
    ! 1 if time reversal symmetry must be used (0 otherwise) at each k point

  integer, allocatable :: tab_ibg(:,:)
    ! tab_ibg,(mkpt,my_nsppol))
    ! indices of cprj(ikpt)/occ(ikpt) in the arrays cprj/occ for each k-point jkpt

  integer, allocatable :: tab_icg(:,:)
    ! tab_icg,(mkpt,my_nsppol))
    ! indices of cg(ikpt) in the arrays cg for each k-point jkpt

  integer, allocatable :: tab_icp(:,:)
    ! tab_icg,(mkpt,my_nsppol))
    ! indices of cprj(ikpt) in the arrays cprj for each k-point jkpt


  integer, allocatable :: tab_ikpt(:)
    ! tab_ikpt,(mkpt))
    ! indices of k-point ikpt in IBZ which corresponds to each k-point jkpt in full BZ

  real(dp), allocatable :: cgocc(:,:,:)
    ! cgocc(2,npw*mkptband,my_nsppol)
    ! wavefunction in the G-space

  real(dp), allocatable :: cwaveocc_bz(:,:,:,:,:,:)
    ! (2,n4,n5,n6,mkptband,my_nsppol))
    ! occupied states of each bands at each k point (used to construct Fock operator), in the real space
  real(dp), allocatable :: occ_bz(:,:)
    ! occ_bz(mkptband,my_nsppol))
    ! occupancy of each bands at each k point

  real(dp), allocatable :: wtk_bz(:)
    ! wtk_bz,(mkpt))
    ! weights assigned to each k point in the BZ
    ! Caution, the definition takes into account "ucvol" !

  real(dp), allocatable :: kptns_bz(:,:)
    ! kptns_bz(3,mkpt)
    ! k-points in full BZ

  real(dp), allocatable :: phase(:,:)
    ! phase(2,mpw*mkpt))
    ! phase factor the cg array will be multiplied with at each k point

  type(MPI_type) :: mpi_enreg
  type(pawang_type),pointer :: pawang
  type(pawcprj_type), allocatable :: cwaveocc_prj(:,:)

 end type fock_BZ_type
!----------------------------------------------------------------------

 type,public :: fock_ACE_type

! ===== Real pointers
  real(dp), allocatable :: xi(:,:,:)

 end type fock_ACE_type
!----------------------------------------------------------------------

 public :: fock_init                  ! Initialize the object.
 public :: fock_set_ieigen            ! Set the value of ieigen to the value given in argument.
 public :: fock_updateikpt            ! Update the value of energies%e_xc and energies%e_xcdc with Fock contribution.
 public :: fock_destroy               ! Free memory.
 public :: fock_ACE_destroy           ! Free memory.
 public :: fock_common_destroy        ! Free memory.
 public :: fock_bz_destroy            ! Free memory.
 public :: fock_calc_ene              ! Calculate the Fock contribution to the total energy.
 public :: fock_update_exc            ! Update the value of energies%e_xc and energies%e_xcdc with Fock contribution.
 public :: fock_updatecwaveocc        ! Update in the fock datastructure the fields relative to the occupied states.
 public :: fock_set_getghc_call       ! Enable/disable the call to fock_getghc in getghc.
 public :: fock_get_getghc_call       ! Return the value of the flag used to enable/disable the call to fock_getghc in getghc.
 public :: fock_print                 ! Print info on the object.
!!***

 ! Help functions
 public :: bare_vqg
 public :: strfock

contains
!!***

!!****f* m_fock/fockbz_create
!! NAME
!!  fockbz_create
!!
!! FUNCTION
!!  Create a fock_type structure.
!!
!! INPUTS
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange are allocated
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!!  The current version is restricted to the case nsym=1, nspinor=1 and mkmem/=0.
!!
!! PARENTS
!!      m_fock
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fockbz_create(fockbz,mgfft,mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,use_ACE)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mgfft,mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,use_ACE
 type(fock_BZ_type) , intent(inout) :: fockbz

!Local variables-------------------------------
!scalars
!arrays
! character(len=500) :: message                   ! to be uncommented, if needed

! *************************************************************************

 !write (std_out,*) ' fockbz_create : enter'

!* Create the array %kptns_bz = the k points in full BZ
 ABI_ALLOCATE(fockbz%kptns_bz,(3,mkpt))
 fockbz%kptns_bz=zero
!* Create the array %jstwfk = how is stored the wavefunction at each k-point
!* By default, the table is initialized to 1 (do NOT take advantage of the time-reversal symmetry)
 ABI_ALLOCATE(fockbz%istwfk_bz,(mkpt))
 fockbz%istwfk_bz=1
!* Create the array %wtk_bz = weight assigned to each k point.
 ABI_ALLOCATE(fockbz%wtk_bz,(mkpt))
 fockbz%wtk_bz=zero
!* Create the array %npwarr_bz = number of planewaves in basis at each k-point
!    ABI_ALLOCATE(fockbz%npwarr_bz,(mkpt))
!    fockbz%npwarr_bz=0
!* Create the array %kg_bz = reduced planewave coordinates at each k-point
 ABI_ALLOCATE(fockbz%kg_bz,(3,mpw*mkpt))
 fockbz%kg_bz=0
!* Create the array %gbound_bz = boundary of the basis sphere of G vectors at each k-point
 ABI_ALLOCATE(fockbz%gbound_bz,(2*mgfft+8,2,mkpt))
 fockbz%gbound_bz=0

!* Create the array %tab_ikpt = indices of k-point ikpt in IBZ which corresponds to each k-point jkpt in full BZ
 ABI_ALLOCATE(fockbz%tab_ikpt,(mkpt))
 fockbz%tab_ikpt=0
!* Create the array %tab_symkpt =indices of symmetry operation to apply to get jkpt in full BZ from ikpt in IBZ
 ABI_ALLOCATE(fockbz%tab_symkpt,(mkpt))
 fockbz%tab_symkpt=0
!* Create the array %tab_ibg = indices of occ(ikpt) in the arrays cprj/occ for each k-point jkpt
 ABI_ALLOCATE(fockbz%tab_ibg,(mkpt,my_nsppol))
 fockbz%tab_ibg=0
!* Create the array %tab_icp = indices of cprj(ikpt) in the arrays cprj/occ for each k-point jkpt
 ABI_ALLOCATE(fockbz%tab_icp,(mkpt,my_nsppol))
 fockbz%tab_icp=0
!* Create the array %tab_icg = indices of cg(ikpt) in the arrays cg for each k-point jkpt
 ABI_ALLOCATE(fockbz%tab_icg,(mkpt,my_nsppol))
 fockbz%tab_icg=0

!* Create the array %calc_phase = 1 if a phase factor must be considered (0 otherwise) at each k point
 ABI_ALLOCATE(fockbz%calc_phase,(mkpt))
 fockbz%calc_phase=0
!* Create the array %phase = phase factor the cg array will be multiplied with at each k point
 ABI_ALLOCATE(fockbz%phase,(2,mpw*mkpt))
 fockbz%phase=zero

!* Create the array %timerev i= 1 if time reversal symmetry must be used (0 otherwise) at each k point
 ABI_ALLOCATE(fockbz%timerev,(mkpt))
 fockbz%timerev=0

!* Create the array %cwaveocc_bz = wavefunctions of each bands at each k point

 if (use_ACE==1) then
   ABI_ALLOCATE(fockbz%cgocc,(2,mpw*mkptband,my_nsppol))
   fockbz%cgocc=zero
 else
   ABI_ALLOCATE(fockbz%cwaveocc_bz,(2,n4,n5,n6,mkptband,my_nsppol))
   fockbz%cwaveocc_bz=zero
 end if
!* Create the array %occ_bz = occupancy of each bands at each k point => will be limited to only the occupied states
 ABI_ALLOCATE(fockbz%occ_bz,(mkptband,my_nsppol))
 fockbz%occ_bz=zero
!* Create the array %nbandocc_bz = nb of bands at each k point
 ABI_ALLOCATE(fockbz%nbandocc_bz,(mkpt,my_nsppol))
 fockbz%nbandocc_bz=0

 ABI_ALLOCATE(fockbz%npwarr,(mkpt))
 fockbz%npwarr=0

end subroutine fockbz_create
!!***

!!****f* m_fock/fock_init
!! NAME
!!  fock_init
!!
!! FUNCTION
!!  Init all scalars, arrays and pointers in the structure.
!!
!! INPUTS
!!  cg(2,mcg)= wavefunctions
!!  dtset <type(dataset_type)>= all input variables for this dataset
!!  gsqcut= Fourier cutoff on G^2 used to calculate charge density
!!  kg(3,mpw*mkmem)= reduced planewave coordinates.
!!  mcg= size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  npwarr_bz(nkpt)= number of planewaves in basis at this k point
!!  occ(mband*nkpt*nsppol)= occupation number for each band (often 2) at each k point
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange are initialized
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!!  The current version is restricted to the case nsym=1, nspinor=1 and mkmem/=0.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_init(atindx,cplex,dtset,fock,gsqcut,kg,mpi_enreg,nattyp,npwarr,pawang,pawfgr,pawtab,rprimd)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: cplex
 real(dp),intent(in) :: gsqcut
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(fock_type),intent(inout),pointer :: fock
 type(pawfgr_type),intent(in),target :: pawfgr
 type(pawang_type),intent(in),target :: pawang
!arrays
 integer, intent(in) :: atindx(dtset%natom),nattyp(dtset%ntypat), npwarr(dtset%nkpt)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 real(dp), intent(in) :: rprimd(3,3)
 type(pawtab_type), intent(in),target :: pawtab(dtset%ntypat*dtset%usepaw)
!Local variables-------------------------------
!scalars
 integer :: iatom,ibg,icg,icp,ier,ik,ikg,ikpt,isppol,isym,itypat,jkpt,jpw,jsym,mband,mgfft,mkpt,mkptband
 integer :: n1,n2,n3,n4,n5,n6,nband,ncpgr,nkpt_bz,nproc_hf,npwj,timrev,use_ACE,v1,v2,v3
 integer :: my_jkpt,jkg_this_proc,my_nsppol,my_nspinor
 real(dp) :: dksqmax,arg
 character(len=500) :: msg
!arrays
 integer :: indx(1),l_size_atm(dtset%natom),shiftg(3),symm(3,3),ident(3,3),symrec(3,3,dtset%nsym)
 real(dp) :: gmet(3,3),gprimd(3,3),tau_nons(3),phktnons(2,1),tsec(2),Rtnons(3,dtset%nsym)
 integer,allocatable :: dimcprj(:),indkk(:,:),kg_tmp(:),my_ikgtab(:),my_ibgtab(:,:),my_icgtab(:,:),my_icptab(:,:),invsym(:)
 real(dp),allocatable :: kptns_hf(:,:), phase1d(:,:)
 type(fock_common_type),pointer :: fockcommon
 type(fock_BZ_type),pointer :: fockbz

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(1501,1,tsec)

 if (dtset%nspinor/=1) then
   MSG_ERROR('Hartree-Fock option can be used only with option nspinor=1.')
 end if


! =====================================
! === Define useful local variables ===
! =====================================

 nkpt_bz=dtset%nkpthf
 nproc_hf=mpi_enreg%nproc_hf
 mband=dtset%nbandhf

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

!* Allocations
 ABI_ALLOCATE(kptns_hf,(3,nkpt_bz))
 kptns_hf=zero
 ABI_ALLOCATE(indkk,(nkpt_bz,6))
 indkk=0
 ABI_ALLOCATE(phase1d,(2,(2*n1+1)*(2*n2+1)*(2*n3+1)))
 phase1d=zero
 ABI_ALLOCATE(kg_tmp,(3*dtset%mpw))

!* Initialize the array my_ikgtab = shifts in arrays kg(ikg) associated to ikpt
 ABI_ALLOCATE(my_ikgtab,(dtset%nkpt))
 ikg=0
 do ikpt=1,dtset%nkpt
   nband=dtset%nband(ikpt)
   if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband,-1,mpi_enreg%me_kpt))) then
!* The point ikpt is treated on this processor.
     my_ikgtab(ikpt)=ikg
!* The array kg is distributed, the shift ikg is incremented only on this proc.
     ikg=ikg+npwarr(ikpt)
   else
     my_ikgtab(ikpt)=-1
!* Default value is -1.
   end if
 end do

!* Initialize the array my_ibgtab = shifts in arrays occ(ibg) associated to ikpt
!* Initialize the array my_icgtab = shifts in arrays cg(icg) associated to ikpt
 ABI_ALLOCATE(my_ibgtab,(dtset%nkpt,dtset%nsppol))
 ABI_ALLOCATE(my_icgtab,(dtset%nkpt,dtset%nsppol))
 ABI_ALLOCATE(my_icptab,(dtset%nkpt,dtset%nsppol))
 ibg=0; icg=0 ;icp=0
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     nband=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband,isppol,mpi_enreg%me_kpt))) then
!* The states with (ikpt,isppol) are stored on this processor.
       my_icgtab(ikpt,isppol)=icg
       my_icptab(ikpt,isppol)=icp
!* The array cg is distributed, the shift icg is incremented only on this proc.
       icg=icg+npwarr(ikpt)*nband
       icp=icp+nband
     else
       my_icgtab(ikpt,isppol)=-1
       my_icptab(ikpt,isppol)=-1
!* Otherwise, the states with (ikpt,isspol) are not stored on this processor and default value is -1.
     end if
!* The array occ is shared among the proc, the shift ibg is always incremented.
     my_ibgtab(ikpt,isppol)=ibg
     ibg=ibg+nband
   end do
 end do

 if (.not.(associated(fock))) then

! =================================
! === Create the fock structure ===
! =================================
   ABI_DATATYPE_ALLOCATE(fock,)
   ABI_DATATYPE_ALLOCATE(fock%fock_common,)
   ABI_DATATYPE_ALLOCATE(fock%fock_BZ,)
! ========================================================
! === Set all the other state-dependent fields to zero ===
! ========================================================
   fockcommon=>fock%fock_common
   fockbz=> fock%fock_BZ
   ABI_ALLOCATE(fockcommon%nband,(dtset%nkpt*dtset%nsppol))
   do ikpt=1,dtset%nkpt*dtset%nsppol
     fockcommon%nband(ikpt)=dtset%nband(ikpt)
   end do

   nband=dtset%mband
   fockcommon%ikpt= 0
!* Will contain the k-point ikpt of the current state
   fockcommon%isppol= 0
!* Will contain the spin isppol of the current state
   fockcommon%ieigen=0
!* Will contain the band index of the current state
!* if the value is 0, the Fock contribution to the eigenvalue is not calculated.
   ABI_ALLOCATE(fockcommon%eigen_ikpt,(nband))
   fockcommon%eigen_ikpt=0.d0
!* Will contain the Fock contributions to the eigenvalue of the current state

!* Compute the dimension of arrays in "spin" w.r.t parallelism
   my_nsppol=dtset%nsppol
   if (mpi_enreg%nproc_kpt>1) my_nsppol=1
!* my_nsppol=1 when nsppol=1 or nsppol=2 and only one spin is treated by the processor.
!* my_nsppol=2 when nsppol=2 and no parallelization over kpt (both spins are treated by the processor).

!* Compute mkpt the size of arrays/pointers for k points w.r.t. parallelism
!* Compute mkptband the size of arrays/pointers for occupied states w.r.t. parallelism
   if (nproc_hf<nkpt_bz) then
!* Parallelization over kpts only
     mkpt=nkpt_bz/nproc_hf
     if (mod(nkpt_bz,nproc_hf) /=0) mkpt=mkpt+1
     mkptband=mkpt*mband
   else
!* Parallelization over occupied states
     if (nproc_hf<nkpt_bz*mband) then
       mkptband=(nkpt_bz*mband)/nproc_hf
       if (mod((nkpt_bz*mband),nproc_hf) /=0) mkptband=mkptband+1
       mkpt=1
       if (mod(nproc_hf,nkpt_bz) /=0) mkpt=2
     else
       mkptband=1
       mkpt=1
     end if
   end if

! mpi_enreg settings
   call copy_mpi_enreg(mpi_enreg,fockbz%mpi_enreg)
   fockbz%mpi_enreg%me_kpt=mpi_enreg%me_hf
   fockbz%mpi_enreg%comm_kpt=mpi_enreg%comm_hf
   fockbz%mpi_enreg%nproc_kpt=mpi_enreg%nproc_hf
   if (allocated(fockbz%mpi_enreg%proc_distrb)) then
     ABI_DEALLOCATE(fockbz%mpi_enreg%proc_distrb)
   end if
   ABI_ALLOCATE(fockbz%mpi_enreg%proc_distrb,(nkpt_bz,mband,1))
   do jkpt=1,nkpt_bz
     fockbz%mpi_enreg%proc_distrb(jkpt,:,1)=fockbz%mpi_enreg%me_kpt
   end do

   mgfft=dtset%mgfft
   fockcommon%usepaw=dtset%usepaw
   if (fockcommon%usepaw==1)then
     mgfft=dtset%mgfftdg
     n4=dtset%ngfftdg(4) ; n5=dtset%ngfftdg(5) ; n6=dtset%ngfftdg(6)
   end if
   fockcommon%optfor=.false.; fockcommon%optstr=.false.
   if(dtset%optforces==1) fockcommon%optfor=.true.
   if (fockcommon%optfor) then
     ABI_ALLOCATE(fockcommon%forces_ikpt,(3,dtset%natom,nband))
     ABI_ALLOCATE(fockcommon%forces,(3,dtset%natom))
     fockcommon%forces=zero
   endif
   use_ACE=1 ! Default. Normal users do not have access to this variable, although the next line allows experts to make tests.
   if(dtset%userie==1729)use_ACE=0 ! Hidden possibility to disable ACE

   fockcommon%use_ACE=use_ACE
   call fockbz_create(fockbz,mgfft,dtset%mpw,mkpt,mkptband,my_nsppol,n4,n5,n6,use_ACE)

!* Initialize %mband, %mkpt, %mkptband = size of arrays
   fockcommon%mband=mband
   fockbz%mkpt=mkpt
   fockbz%mkptband=mkptband
   fockcommon%my_nsppol = my_nsppol
   fockcommon%nsppol = dtset%nsppol
   if (fockcommon%use_ACE/=0) then
     ABI_DATATYPE_ALLOCATE(fock%fockACE,(dtset%nkpt,dtset%nsppol))
     my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
     do isppol=1,dtset%nsppol
       do ikpt=1,dtset%nkpt
         nband=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
         ABI_ALLOCATE(fock%fockACE(ikpt,isppol)%xi,(2,npwarr(ikpt)*my_nspinor,nband))
       end do
     end do
   end if
!========Initialze PAW data========
   fockcommon%ntypat=dtset%ntypat
   fockcommon%natom=dtset%natom
   if (fockcommon%usepaw==1) then
     fockbz%mcprj=mkptband*my_nsppol
     fockcommon%pawfgr => pawfgr
     fockbz%pawang => pawang
     fockcommon%pawtab => pawtab
     ABI_DATATYPE_ALLOCATE(fockcommon%pawfgrtab,(dtset%natom))
     do iatom = 1, dtset%natom
       itypat=dtset%typat(iatom)
       l_size_atm(iatom) = pawtab(itypat)%lcut_size
     end do
     call pawfgrtab_init(fockcommon%pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat)
     call pawfgrtab_nullify(fockcommon%pawfgrtab)
     ABI_DATATYPE_ALLOCATE(fockbz%cwaveocc_prj,(dtset%natom,fockbz%mcprj))
     ABI_ALLOCATE(dimcprj,(dtset%natom))
     call pawcprj_getdim(dimcprj,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
     ncpgr = 0
     if (dtset%optforces== 1) ncpgr = 3
!     if (dtset%optstress /= 0) ncpgr = 6
!     ncpgr=3*dtset%optforces+6*dtset%optstress
     call pawcprj_alloc(fockbz%cwaveocc_prj,ncpgr,dimcprj)
     ABI_DEALLOCATE(dimcprj)
     ABI_ALLOCATE(fockcommon%atindx,(dtset%natom))
     fockcommon%atindx=atindx
     ABI_ALLOCATE(fockcommon%typat,(dtset%natom))
     fockcommon%typat=dtset%typat
   end if
! ==========================================
! === Initialize the convergence options ===
! ==========================================
   write(msg,'(2a)') ch10,'Fock_init: initialization of Fock operator parameters:'
   call wrtout(std_out,msg,'COLL')

   fockcommon%fock_converged=.false.
   fockcommon%scf_converged=.false.

!* Number of iterations with fixed occupied states when calculating the exact exchange contribution.
   if (dtset%nnsclohf<0) then
     MSG_ERROR('The parameter nnsclohf must be a non-negative integer.')
   end if
   if (dtset%nnsclohf==0) then
     fockcommon%nnsclo_hf=1
     msg=' - The parameter nnsclohf is set to its default value 1.'
     call wrtout(std_out,msg,'COLL')
!* Default value is set to 1 (updating cgocc at each step)
!* May be useful to put default to 3
   else
     fockcommon%nnsclo_hf=dtset%nnsclohf
     write(msg,'(a,i3)') ' - The parameter nnsclohf is set to the value:', dtset%nnsclohf
     call wrtout(std_out,msg,'COLL')
!* value chosen by the user
   end if

! =========================================
! === Initialize the hybrid coefficient ===
! =========================================
   fockcommon%ixc = dtset%ixc
!  By convention, positive values are the default values for the ixc,
!  while negative values have been set by the user (and stored as negative numbers)
   fockcommon%hyb_mixing=abs(dtset%hyb_mixing)
   fockcommon%hyb_mixing_sr=abs(dtset%hyb_mixing_sr)
   fockcommon%hyb_range_dft=abs(dtset%hyb_range_dft)
   fockcommon%hyb_range_fock=abs(dtset%hyb_range_fock)

!  Set the hybrid parameters if functional from libxc for which parameters can be changed, or if the user asked to do so.
!  Usually, these parameters were obtained from libxc,
!  but the user might have possibly modified them. By the way, must define them here for the usual changeable fonctionals,
!  since otherwise might inherit them from the previous dataset !
   if(dtset%ixc<0)then
     if (dtset%ixc==-406.or.dtset%ixc==-427.or.dtset%ixc==-428 .or. &
&      min(dtset%hyb_mixing,dtset%hyb_mixing_sr,dtset%hyb_range_dft,dtset%hyb_range_fock)<-tol8)then
       call libxc_functionals_set_hybridparams(hyb_mixing=fockcommon%hyb_mixing,&
&                                              hyb_mixing_sr=fockcommon%hyb_mixing_sr,&
&                                              hyb_range=fockcommon%hyb_range_dft)
     end if
   end if


! ======================================================
! === Initialize the data relative to Poisson solver ===
! ======================================================

!* gsqcut = cutoff value on G^2 for sphere inside the fft box (input for vhartre).
   fockcommon%gsqcut= gsqcut

! =======================================================
! === Initialize the properties of the k-points in BZ ===
! =======================================================
!* Initialize %nkpt_bz = nb of k point in BZ for the calculation of exchange
   fockbz%nkpt_bz=nkpt_bz
!* Initialize the array %wtk_bz = weight assigned to each k point.
   fockbz%wtk_bz=1.0_dp/dble(nkpt_bz)


   if (dtset%kptopt>=1 .and. dtset%kptopt<=4) then
! ============================================
! === Initialize the set of k-points in BZ ===
! ============================================
     kptns_hf(:,1:nkpt_bz)=dtset%kptns_hf(:,1:nkpt_bz)
!* kptns_hf contains the special k points obtained by the Monkhorst & Pack method, in reduced coordinates. (output)

! =======================================================
! === Compute the transformation to go from IBZ to BZ ===
! =======================================================
!* Compute the reciprocal space metric.
     call matr3inv(rprimd,gprimd)
     gmet = MATMUL(TRANSPOSE(gprimd),gprimd)

!* Calculate the array indkk which describes how to get IBZ from BZ
!* dksqmax=maximal value of the norm**2 of the difference between a kpt2 vector and the closest k-point found from the kptns1 set, using symmetries. (output)
!* sppoldbl=1, no spin-polarisation doubling is required.
     timrev=1 ; if (dtset%kptopt==3 .or. dtset%kptopt==4) timrev=0
!* timrev=1 if the use of time-reversal is allowed ; 0 otherwise
     if (dtset%kptopt==2 .or. dtset%kptopt==3) then
!* No space symmetry is used, if kptopt==2 time reversal symmetry is used.
       symm=0 ; symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1
       call listkk(dksqmax,gmet,indkk(1:nkpt_bz,:),dtset%kptns,kptns_hf,dtset%nkpt, &
&          nkpt_bz,1,1,indx,symm,timrev,xmpi_comm_self)
     else
!* As in getkgrid, no use of antiferromagnetic symmetries thans to the option sppoldbl=1
       call listkk(dksqmax,gmet,indkk(1:nkpt_bz,:),dtset%kptns,kptns_hf,dtset%nkpt, &
&          nkpt_bz,dtset%nsym,1,dtset%symafm,dtset%symrel,timrev, xmpi_comm_self)
     end if
!* indkk(nkpt_bz,6) describes the k point of IBZ that generates each k point of BZ
!*    indkk(:,1)   = k point of IBZ, kpt_ibz
!*    indkk(:,2)   = symmetry operation to apply to kpt_ibz to give the k point of BZ
!*                   (if 0, means no symmetry operation, equivalent to identity )
!*    indkk(:,3:5) = Umklapp vectors to apply to remain in BZ
!*    indkk(:,6)   = 1 if time-reversal was used to generate the k point of BZ, 0 otherwise
!* No use of symafm to generate spin down wfs from spin up wfs for the moment

   else
     if (dtset%kptopt==0) then
!* kptopt =0 : read directly nkpt, kpt, kptnrm and wtk in the input file
!*              => this case is not allowed for the moment
       MSG_ERROR('Hartree-Fock option can not be used with option kptopt=0.')
     else
!* kptopt <0 : rely on kptbounds, and ndivk to set up a band structure calculation
!*              => a band structure calculation is not yet allowed.
       MSG_ERROR('Hartree-Fock option can not be used with option kptopt<0.')
     end if
   end if

!! =======================================================
!! === Initialize the properties of the k-points in BZ ===
!! =======================================================
!       jkg=0
!!* Initialize the arrays %npwarr_bz, %kg_j, %phase_j, %gbound_j
!       do jkpt=1,nkpt_bz
!         ikpt=indkk(jkpt,1)
!!* ikpt = the point of IBZ that jkpt is an image of in BZ
!         npwj=npwarr(ikpt)
!!* npwj = number of planewaves in basis at point jkpt = at point ikpt
!         jsym=indkk(jkpt,2)
!!* jsym = symmetry operation to apply to get jkpt from ikpt
!         shiftg(:)=indkk(jkpt,3:5)
!!* shiftg = Bravais vector G0 to add to remain in BZ
!         if (jsym/=0) then
!           symm(:,:)=dtset%symrel(:,:,jsym)
!           tau_nons(:)=dtset%tnons(:,jsym)
!!* The symmetry operation in k-space (symm) and the non-symorphic translation (tau_nons) are now defined.
!           if(sum(tau_nons(:)**2)>tol8) then
!!* Initialize %calc_phase(jkpt) to 1
!             fock%calc_phase(jkpt)=1
!!* Compute the phase factor exp(i*2*pi*G.tau) for all G.
!             indx(1)=1
!             phase1d=zero
!             call getph(indx,1,n1,n2,n3,phase1d,tau_nons)
!!* Although the routine getph is orignally written for atomic phase factors, it does precisely what we want
!             arg=two_pi*(dtset%kptns(1,ikpt)*tau_nons(1) + dtset%kptns(2,ikpt)*tau_nons(2) &
!&                + dtset%kptns(3,ikpt)*tau_nons(3))
!             phktnons(1,1)=cos(arg)
!             phktnons(2,1)=sin(arg)
!!              phktnons(1,1)=one
!!              phktnons(2,1)=zero
!!* Convert 1D phase factors to 3D phase factors exp(i*2*pi*(k+G).tau) and store it in %phase_j
!             call ph1d3d(1,1,kg(:,1+tab_indikpt(1,ikpt):npwj+tab_indikpt(1,ikpt)),1,1,npwj,n1, &
!&              n2,n3,phktnons,phase1d,fock%phase(:,1+jkg:npwj+jkg))
!           end if
!         else
!           symm=0 ; symm(1,1)=1 ; symm(2,2)=1 ; symm(3,3)=1
!           tau_nons(:)=zero
!           shiftg(:)=0
!         end if
!!* Apply time-reversal symmetry if required
!         if(indkk(jkpt,6)/=0) then
!!* Initialize %timerev(jkpt) to 1
!           fock%timerev(jkpt)=1
!           symm(:,:)=-symm(:,:)
!         end if

!!* Initialize %istwfk_bz(jkpt) to
!         fock%istwfk_bz(jkpt)=dtset%istwfk(ikpt)

!!* Initialize %tab_ikpt and %tab_ibgcg
!         fock%tab_ikpt(jkpt)=ikpt
!         fock%tab_ibgcg(1:dtset%nsppol,jkpt)=tab_indikpt(2:1+dtset%nsppol,ikpt)
!         fock%tab_ibgcg(1+dtset%nsppol:2*dtset%nsppol,jkpt)= &
!&          tab_indikpt(2+dtset%nsppol:2*dtset%nsppol+1,ikpt)

!!* Initialize %npwarr_bz
!         fock%npwarr_bz(jkpt)=npwj

!!* Initialize %kg_bz
!         do jpw=1,npwj
!           v1=kg(1,jpw+tab_indikpt(1,ikpt)) ; v2=kg(2,jpw+tab_indikpt(1,ikpt)) ; v3=kg(3,jpw+tab_indikpt(1,ikpt))
!           fock%kg_bz(1,jpw+jkg)=-shiftg(1)+symm(1,1)*v1+symm(2,1)*v2+symm(3,1)*v3
!           fock%kg_bz(2,jpw+jkg)=-shiftg(2)+symm(1,2)*v1+symm(2,2)*v2+symm(3,2)*v3
!           fock%kg_bz(3,jpw+jkg)=-shiftg(3)+symm(1,3)*v1+symm(2,3)*v2+symm(3,3)*v3
!!* The symmetry operation symm must be transposed when used. (cf. docs about wfconv)
!         end do

!!* Initialize %gbound_bz
!         call sphereboundary(fock%gbound_bz(:,:,jkpt),fock%istwfk_bz(jkpt), &
!&          fock%kg_bz(:,1+jkg:npwj+jkg),dtset%mgfft,npwj)

!!* Update of the shift to be applied
!         jkg=jkg+npwj
!       end do

! ==========================================================
! === Initialize the k-points in BZ and their properties ===
! ==========================================================
!   jkg=0;

   do isym=1,dtset%nsym
         call mati3inv(dtset%symrel(:,:,isym),symrec(:,:,isym))
         Rtnons (:,isym)= MATMUL(TRANSPOSE(symrec(:,:,isym)),dtset%tnons(:,isym))
   end do
   ABI_ALLOCATE(fockcommon%symrec,(3,3,dtset%nsym))
   fockcommon%symrec=symrec

   ABI_ALLOCATE(invsym,(dtset%nsym))
   invsym=0
   ident(1,:3)=(/1,0,0/)
   ident(2,:3)=(/0,1,0/)
   ident(3,:3)=(/0,0,1/)
   do isym=1,dtset%nsym
     symm(:,:)=MATMUL(dtset%symrel(:,:,isym),dtset%symrel(:,:,isym))
     if (all(symm(:,:)==ident(:,:))) then
       invsym(isym)=isym
     else
       do jsym=1,dtset%nsym
         symm(:,:)=MATMUL(dtset%symrel(:,:,isym),dtset%symrel(:,:,jsym))
         if (all(symm(:,:)==ident(:,:))) then
            invsym(isym)=jsym
            cycle
         end if
       end do
     end if
     if(invsym(isym)==0) then
       MSG_ERROR('No inverse has been found for isym')
     end if
   end do

   jkg_this_proc=0;my_jkpt=0
!indkk(1:nkpt_bz,2)=(/1,1,3,1,11,7,9,1/)
   do jkpt=1,nkpt_bz

!* If this processor does not calculate exchange with the k point jkpt, skip the rest of the k-point loop.
     if (proc_distrb_cycle(mpi_enreg%distrb_hf,jkpt,1,mband,1,mpi_enreg%me_hf)) cycle
!       if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,jkpt,1,dtset%nbandhf,1,mpi_enreg%me_kpt))) then
!* The processor does own a copy of the array kg of ikpt ; increment the shift.
!         jkg=jkg+npwj
!        end if
! Skip the rest of the k-point loop
!       cycle
!     end if
     my_jkpt=my_jkpt+1

     ikpt=indkk(jkpt,1)
!* ikpt = the point of IBZ that jkpt is an image of in BZ
     npwj=npwarr(ikpt)
     fockbz%npwarr(my_jkpt)=npwarr(ikpt)
!* npwj = number of planewaves in basis at point jkpt = at point ikpt
     jsym=indkk(jkpt,2)
!* jsym = symmetry operation to apply to get jkpt from ikpt
     fockbz%tab_symkpt(my_jkpt)=invsym(jsym)
     shiftg(:)=indkk(jkpt,3:5)
!* shiftg = Bravais vector G0 to add to remain in BZ

!* Initialize the array %kptns_bz = the k points in full BZ
     fockbz%kptns_bz(:,my_jkpt)=kptns_hf(:,jkpt)

!* Initialize the array %jstwfk = how is stored the wavefunction at each k point
     if (dtset%istwfk(ikpt)/=1) then
       fockbz%istwfk_bz(my_jkpt)=set_istwfk(kptns_hf(:,jkpt))
     end if

!* One can take advantage of the time-reversal symmetry in this case.
!* Initialize the array %wtk_bz = weight assigned to each k point.
!     fock%wtk_bz(my_jkpt)=dtset%wtk(jkpt)/ucvol
!* Caution, the definition takes into account "ucvol" !

!* Initialize the array %npwarr_bz = number of planewaves in basis at each k point
!     fock%npwarr_bz(my_jkpt)=npwj

!!* Initialize the array %tab_ikpt = indices of k-point in IBZ ikpt for each k point jkpt in BZ (here,ikpt=jkpt)
     fockbz%tab_ikpt(my_jkpt)=ikpt


!!* Initialize the array %tab_ibgcg = indices of cprj(ikpt)/occ(ikpt) and cg(ikpt) for each k point jkpt
!     if (my_nsppol==2) then
!!* In this case, my_nsppol=dtset%nsppol=2
!       fock%tab_ibgcg(1:2,my_jkpt)=tab_indikpt(2:3,ikpt)
!       fock%tab_ibgcg(3:4,my_jkpt)=tab_indikpt(4:5,ikpt)
!     else
!       if(mpi_enreg%my_isppoltab(1)==1) then
!!* In this case, my_nsppol=1 and the up spin is treated (dtset%nsppol= 1 or 2)
!         fock%tab_ibgcg(1,my_jkpt)=tab_indikpt(2,ikpt)
!         fock%tab_ibgcg(2,my_jkpt)=tab_indikpt(2+dtset%nsppol,ikpt)
!       else
!!* In this case, my_nsppol=1 and the dn spin is treated (so dtset%nsppol=2)
!         fock%tab_ibgcg(1,my_jkpt)=tab_indikpt(3,ikpt)
!         fock%tab_ibgcg(2,my_jkpt)=tab_indikpt(5,ikpt)
!       end if
!     end if

!* Initialize the array %kg_bz = reduced planewave coordinates at each k point
     if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me_kpt))) then
!* We perform the test with isppol=-1 (both spins) and the occupied band (dtset%nbandhf).
!* We assume that paral_kgb==0 (a k-point may not be present on several proc.)
!* The array kg for ikpt is stored on this processor and copied in kg_tmp.
       ikg=my_ikgtab(ikpt)
!* ikg = the shift in kg to get the G-vectors associated to ikpt
       do ik=1,3
!         kg_tmp(1+(ik-1)*npwj:ik*npwj)=kg(ik,1+tab_indikpt(1,ikpt):npwj+tab_indikpt(1,ikpt))
         kg_tmp(1+(ik-1)*npwj:ik*npwj)=kg(ik,1+ikg:npwj+ikg)
       end do
!       jkg=jkg+npwj
     end if
!* Broadcast the array kg_tmp to all the processors of comm_kpt.
!* Since paral_kgb==0, all the bands of a k-point are treated on the same proc.
     call xmpi_bcast(kg_tmp,mpi_enreg%proc_distrb(ikpt,1,1),mpi_enreg%comm_kpt,ier)
     do ik=1,3
       fockbz%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)=kg_tmp(1+(ik-1)*npwj:ik*npwj)
     end do

!* Apply a symmetry operation on kg_bz if necessary
     if (jsym/=0) then
       symm(:,:)=dtset%symrel(:,:,jsym)
!      tau_nons(:)=dtset%tnons(:,jsym)
       tau_nons(:)=-Rtnons(:,invsym(jsym))
!* The symmetry operation in k-space (symm) and the non-symorphic translation (tau_nons) are now defined.
       if(sum(tau_nons(:)**2)>tol8) then
!* Initialize %calc_phase(jkpt) to 1
         fockbz%calc_phase(my_jkpt)=1
!* Compute the phase factor exp(i*2*pi*G.tau) for all G.
         indx(1)=1
         phase1d=zero
         call getph(indx,1,n1,n2,n3,phase1d,tau_nons)
!* Although the routine getph is orignally written for atomic phase factors, it does precisely what we want
         arg=two_pi*(dtset%kptns(1,ikpt)*tau_nons(1) + dtset%kptns(2,ikpt)*tau_nons(2) &
&            + dtset%kptns(3,ikpt)*tau_nons(3))
         phktnons(1,1)=cos(arg)
         phktnons(2,1)=sin(arg)
!          phktnons(1,1)=one
!          phktnons(2,1)=zero
!* Convert 1D phase factors to 3D phase factors exp(i*2*pi*(k+G).tau) and store it in %phase_j
         call ph1d3d(1,1,fockbz%kg_bz(:,1+jkg_this_proc:npwj+jkg_this_proc),1,1,npwj,n1,n2,n3, &
&          phktnons,phase1d,fockbz%phase(:,1+jkg_this_proc:npwj+jkg_this_proc))
       end if
!* Apply time-reversal symmetry if required
       if(indkk(jkpt,6)/=0) then
!* Initialize %timerev(jkpt) to 1
         fockbz%timerev(my_jkpt)=1
         symm(:,:)=-symm(:,:)
       end if
!* Initialize %kg_bz
       do jpw=1,npwj
         v1=fockbz%kg_bz(1,jpw+jkg_this_proc) ; v2=fockbz%kg_bz(2,jpw+jkg_this_proc) ; v3=fockbz%kg_bz(3,jpw+jkg_this_proc)
         fockbz%kg_bz(1,jpw+jkg_this_proc)=-shiftg(1)+symm(1,1)*v1+symm(2,1)*v2+symm(3,1)*v3
         fockbz%kg_bz(2,jpw+jkg_this_proc)=-shiftg(2)+symm(1,2)*v1+symm(2,2)*v2+symm(3,2)*v3
         fockbz%kg_bz(3,jpw+jkg_this_proc)=-shiftg(3)+symm(1,3)*v1+symm(2,3)*v2+symm(3,3)*v3
!* The symmetry operation symm must be transposed when used. (cf. docs about wfconv)
       end do
     else
!* Ths symmetry operation is the identity.
!* Apply time-reversal symmetry if required
       if(indkk(jkpt,6)/=0) then
!* Initialize %timerev(jkpt) to 1
         fockbz%timerev(my_jkpt)=1
         fockbz%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)=-fockbz%kg_bz(ik,1+jkg_this_proc:npwj+jkg_this_proc)
       end if
     end if

!* Initialize the array %gbound_bz = boundary of the basis sphere of G vectors at each k point
     call sphereboundary(fockbz%gbound_bz(:,:,my_jkpt),fockbz%istwfk_bz(my_jkpt),&
&      fockbz%kg_bz(:,1+jkg_this_proc:npwj+jkg_this_proc),mgfft,npwj)

     jkg_this_proc=jkg_this_proc+npwj

!* Initialize the arrays %tab_ibg = shifts in arrays cprj and occ (ibg) for each k point jkpt
!* Initialize the arrays %tab_icg = shifts in arrays cg(icg) for each k point jkpt
     if (my_nsppol==1) then
         fockbz%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,1+mpi_enreg%my_isppoltab(2))
         fockbz%tab_icg(my_jkpt,1)=my_icgtab(ikpt,1+mpi_enreg%my_isppoltab(2))
         fockbz%tab_icp(my_jkpt,1)=my_icptab(ikpt,1+mpi_enreg%my_isppoltab(2))
!* if mpy_isppoltab(2)=0, the up spin is treated (dtset%nsppol= 1 or 2)
!* if mpy_isppoltab(2)=1, the dn spin is treated (so dtset%nsppol=2)

!       if(mpi_enreg%my_isppoltab(2)==1) then
!* In this case, my_nsppol=1 and the up spin is treated (dtset%nsppol= 1 or 2)
!         fock%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,1)
!         fock%tab_icg(my_jkpt,1)=my_icgtab(ikpt,1)
!       else
!* In this case, my_nsppol=1 and the dn spin is treated (so dtset%nsppol=2)
!         fock%tab_ibg(my_jkpt,1)=my_ibgtab(ikpt,2)
!         fock%tab_icg(my_jkpt,1)=my_icgtab(ikpt,2)
!       end if
     else
!* In this case, my_nsppol=dtset%nsppol=2
       fockbz%tab_ibg(my_jkpt,:)=my_ibgtab(ikpt,:)
       fockbz%tab_icg(my_jkpt,:)=my_icgtab(ikpt,:)
       fockbz%tab_icp(my_jkpt,:)=my_icptab(ikpt,:)
     end if

   enddo

!* Deallocation
   ABI_DEALLOCATE(invsym)

 end if
 ABI_DEALLOCATE(indkk)
 ABI_DEALLOCATE(kg_tmp)
 ABI_DEALLOCATE(kptns_hf)
 ABI_DEALLOCATE(my_ibgtab)
 ABI_DEALLOCATE(my_icgtab)
 ABI_DEALLOCATE(my_icptab)
 ABI_DEALLOCATE(my_ikgtab)
 ABI_DEALLOCATE(phase1d)
 call fock_print(fockcommon,fockbz,unit=std_out)

 call timab(1501,2,tsec)

 DBG_EXIT("COLL")

end subroutine fock_init
!!***

!!****f* m_fock/fock_updateikpt
!! NAME
!!  fock_updateikpt
!!
!! FUNCTION
!!  Update the value of ikpt,isppol for the next exact exchange calculation.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  gs_ham <type(gs_hamiltonian_type)>= all data for the Hamiltonian to be applied
!!  ikpt= reduced planewave coordinates.
!!  isppol= number of planewaves in basis at this k point
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!   The field fock%eigen_ikpt is also set to 0.d0.
!!
!! NOTES
!!  May be improved to calculate the star of ikpt. => I think NO finally
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_updateikpt(fock,ikpt,isppol)

!Arguments ------------------------------------
 integer, intent(in) :: ikpt,isppol
 type(fock_common_type),pointer :: fock

! *************************************************************************

 !write (std_out,*) ' fock_updateikpt : enter'

! ======================================================
! === Update the data relative to the current states ===
! ======================================================
!* Copy of the value ikpt in the field ikpt
   fock%ikpt=ikpt
!* Copy of the value isppol in the field isppol
   fock%isppol=isppol
!* Set all the Fock contributions to the eigenvalues to 0.d0.
   fock%eigen_ikpt=zero
!* Set all the Fock contributions to the forces to 0.d0.
   if ((fock%optfor).and.(fock%use_ACE==0)) then
     fock%forces_ikpt=zero
   endif

end subroutine fock_updateikpt
!!***

!!****f* m_fock/fock_set_ieigen
!! NAME
!!  fock_set_ieigen
!!
!! FUNCTION
!!  Set the value of ieigen to the value given in argument.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  iband= index of the band iband
!!
!! OUTPUT
!!  none
!!
!! PARENTS
!!      cgwf
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_set_ieigen(fock,iband)

!Arguments ------------------------------------
 integer, intent(in) :: iband
 type(fock_common_type),pointer :: fock

! *************************************************************************

!Nothing to do if fock pointer is not associated...

! ======================================================
! === Update the data relative to the current states ===
! ======================================================

!* Copy of the value iband in the field ieigen
 if (associated(fock)) then
   fock%ieigen=iband
   fock%iband=iband
 end if

end subroutine fock_set_ieigen
!!***

!!****f* m_fock/fock_destroy
!! NAME
!!  fock_destroy
!!
!! FUNCTION
!!  Clean and destroy fock datastructure.
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE
subroutine fock_destroy(fock)

!Arguments ------------------------------------
 type(fock_type),pointer :: fock

! *************************************************************************

 DBG_ENTER("COLL")
 if (fock%fock_common%use_ACE/=0) then
   ABI_DATATYPE_DEALLOCATE(fock%fockACE)
 end if
 ABI_DATATYPE_DEALLOCATE(fock%fock_common)
 ABI_DATATYPE_DEALLOCATE(fock%fock_BZ)
 ABI_DATATYPE_DEALLOCATE(fock)

 DBG_EXIT("COLL")
end subroutine fock_destroy

subroutine fock_common_destroy(fock)

!Arguments ------------------------------------
 type(fock_common_type),pointer :: fock

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_SFREE(fock%atindx)
 ABI_SFREE(fock%typat)
 ! real arrays
 ABI_SFREE(fock%forces)
 ABI_SFREE(fock%nband)
 ABI_SFREE(fock%forces_ikpt)
 ABI_SFREE(fock%stress_ikpt)
 ABI_SFREE(fock%eigen_ikpt)
 ! Deallocate datatypes
 if (allocated(fock%pawfgrtab)) then
    call pawfgrtab_free(fock%pawfgrtab)
    ABI_FREE(fock%pawfgrtab)
 endif

 ! Put the integer to 0
 fock%ieigen=0
 fock%ikpt=0
 fock%isppol=0

 if (allocated(fock%symrec)) then
    ABI_DEALLOCATE(fock%symrec)
 endif

!* [description of divergence in |q+G|=0]
!* Put the real (dp) to 0
 fock%gsqcut=zero
 fock%hyb_mixing=zero
 fock%hyb_mixing_sr=zero
 fock%hyb_range_dft=zero
 fock%hyb_range_fock=zero

 DBG_EXIT("COLL")
end subroutine fock_common_destroy


subroutine fock_BZ_destroy(fock)

!Arguments ------------------------------------
 type(fock_BZ_type),pointer :: fock

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_SFREE(fock%cwaveocc_bz)
 ABI_SFREE(fock%cgocc)
 ABI_SFREE(fock%npwarr)
 ABI_SFREE(fock%occ_bz)
 if (allocated(fock%cwaveocc_prj)) then
   call pawcprj_free(fock%cwaveocc_prj)
   ABI_FREE(fock%cwaveocc_prj)
 endif
 ! Deallocate integer arrays

 ABI_SFREE(fock%kg_bz)
 ABI_SFREE(fock%nbandocc_bz)
 ABI_SFREE(fock%istwfk_bz)
 ABI_SFREE(fock%calc_phase)
 ABI_SFREE(fock%timerev)
 ABI_SFREE(fock%tab_ibg)
 ABI_SFREE(fock%tab_icg)
 ABI_SFREE(fock%tab_icp)
 ABI_SFREE(fock%tab_ikpt)
 ABI_SFREE(fock%tab_symkpt)

!* [description of IBZ and BZ]
!* Deallocate real arrays
 ABI_SFREE(fock%wtk_bz)
 ABI_SFREE(fock%kptns_bz)
 ABI_SFREE(fock%phase)
!* Put the integer to 0
 fock%nkpt_bz=0

!* Deallocate real arrays

!* Deallocate integer arrays
 ABI_SFREE(fock%gbound_bz)

!* [description of size of arrays/pointers]
!* Put the integer to 0
 fock%mkpt=0
 fock%mkptband=0
 call destroy_mpi_enreg(fock%mpi_enreg)

 DBG_EXIT("COLL")

end subroutine fock_BZ_destroy
!!***
!!****f* m_fock/fock_ACE_destroy
!! NAME
!!  fock_ACE_destroy
!!
!! FUNCTION
!!  Clean and destroy fock datastructure.
!!
!! INPUTS
!!  fockACE <type(fock_ACE_type)>= all the quantities to calculate Fock exact exchange in the ACE context
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_ACE_destroy(fockACE)

!Arguments ------------------------------------
 type(fock_ACE_type),pointer :: fockACE(:,:)
!Local variables-------------------------------
 integer :: dim1,dim2,ii,jj
! *************************************************************************
 DBG_ENTER("COLL")

 dim1=size(fockACE,1)
 dim2=size(fockACE,2)
 do jj=1,dim2
   do ii=1,dim1
     if (allocated(fockACE(ii,jj)%xi)) then
       ABI_DEALLOCATE(fockACE(ii,jj)%xi)
     end if
   end do
 end do
 DBG_EXIT("COLL")

end subroutine fock_ACE_destroy
!!***


!!****f* m_fock/fock_calc_ene
!! NAME
!!  fock_calc_ene
!!
!! FUNCTION
!!  Calculate the Fock contribution to the total energy
!!
!! INPUTS
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  ikpt= reduced planewave coordinates.
!!
!! OUTPUT
!!  none
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_exactX = Fock contribution to the total energy (Hartree)
!!
!! NOTES
!! If the cgocc_bz are not updated at each iteration, be careful to calculate Fock energy at the same frequency.
!! TO CHECK == CHANGE IN SOME DEFINTIONS
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_calc_ene(dtset,fock,fock_energy,ikpt,nband,occ)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,nband
 real(dp),intent(inout) :: fock_energy
 type(dataset_type),intent(in) :: dtset
 type(fock_common_type),pointer :: fock
!arrays
 real(dp),intent(in) :: occ(nband)

!Local variables-------------------------------
 integer :: iband

! *************************************************************************

 ABI_UNUSED(fock_energy)

 do iband=1,nband

   ! Select only the occupied states (such that fock%occ_bz > 10^-8)
   if (abs(occ(iband))>tol8) then
!     fock_energy=fock_energy + half*fock%eigen_ikpt(iband)*occ(iband)*dtset%wtk(ikpt)
     !* Sum the contribution of each occupied states at point k_i
     !* No need to multiply %wtk by ucvol since there is no factor 1/ucvol in the definition of %wtk

!* accumulate Fock contributions to the forces.
!     if (fock%optfor) then
       fock%forces(:,:)=fock%forces(:,:)+occ(iband)*dtset%wtk(ikpt)*fock%forces_ikpt(:,:,iband)
!     endif
   end if
 end do

end subroutine fock_calc_ene
!!***

!!****f* m_fock/fock_update_exc
!! NAME
!!  fock_update_exc
!!
!! FUNCTION
!!  Update the value of energies%e_xc and energies%e_xcdc with Fock contribution
!!
!! INPUTS
!!
!! OUTPUT
!!  none
!!
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_fock= Fock contribution to the total energy (Hartree)
!!
!! NOTES
!!   If the cgocc_bz are not updated at each iteration, be careful to calculate Fock energy at the same frequency.
!!
!! PARENTS
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_update_exc(fock_energy,xc_energy,xcdc_energy)

!Arguments ------------------------------------
 real(dp),intent(in) :: fock_energy
 real(dp),intent(inout) :: xc_energy,xcdc_energy

! *************************************************************************

!xc_energy = fock%hyb_mixing*fock_energy
!xcdc_energy = two*fock%hyb_mixing*fock_energy
 xc_energy =  fock_energy
 xcdc_energy = two*fock_energy
!CMartins : For an atom, ewald should be set to zero (at the beginning of the loop) and
!the contribution in !|q+G|=0 should be an approximation to the missing component of Vloc in G=0
!energies%e_ewald=energies%e_ewald-half*fock%divgq0*fock%wtk_bz(1)*piinv

end subroutine fock_update_exc
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_updatecwaveocc
!! NAME
!!  fock_updatecwaveocc
!!
!! FUNCTION
!!  Update in the fock datastructure the fields relative to the occupied states.
!!
!! INPUTS
!!  cg(2,mcg)= Input wavefunctions
!!  cprj(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  fock <type(fock_type)>= all the quantities to calculate Fock exact exchange
!!  indsym(4,nsym,natom) :: 1:3 shift, and 4 final atom, of symmetry isym operating on iatom
!!                            (S^{-1}(R - t) = r0 + L, see symatm.F90
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  occ(mband*nkpt*nsppol)= occupation number for each band (often 2) at each k point
!!  ucvol= unit cell volume ($\textrm{bohr}^{3}$)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   The field fock%cgocc_bz contains the table cg at the end.
!!   The fields kg_bz, occ_bz and fock%cwaveocc_prj are simultaneously updated.
!!
!! NOTES
!!
!!  ############################
!!  ### Not fully tested yet ###
!!  ############################
!!
!! May be improved by selecting only the occupied states with the same spin isppol.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_updatecwaveocc(cg,cprj,dtset,fock,indsym,mcg,mcprj,&
&                              mpi_enreg,nattyp,npwarr,occ,ucvol)

!scalars
 integer, intent(in) :: mcg,mcprj
 real(dp), intent(in) :: ucvol
 type(dataset_type),intent(in) :: dtset
 type(fock_type),intent(inout),pointer :: fock
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom),nattyp(dtset%ntypat),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0
 integer :: iatm,iatom,iband,iband0,iband_cprj,ibg,icg,icp,ier,ikpt,ilmn,isize,ispinor,isppol,itypat,jbg,jcg,jkg,jkpt,jpw,jstwfk!,ii1,ii2
 integer :: lmnmax,mband,mband0,mgfft,mkpt,mpw,my_jsppol,my_jband,my_jkpt
 integer :: nband,ncpgr,n4,n5,n6,nkpt_bz,npwj,nsppol,nspinor
 real(dp),parameter :: weight1=one
 real(dp) :: cgre,cgim,invucvol
 character(len=500) :: message
! arrays
 integer :: ngfft(18)
 integer, ABI_CONTIGUOUS pointer :: gbound_k(:,:),kg_k(:,:)
 integer,allocatable :: dimlmn(:),indlmn(:,:,:),indsym_(:,:,:),typat_srt(:)
 real(dp) :: tsec(2),tsec2(2),dcp(3)
 real(dp),allocatable :: cgocc_tmp(:),cgocc(:,:),dummytab2(:,:),dummytab3(:,:,:),phase_jkpt(:,:)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)
 type(fock_common_type),pointer :: fockcommon
 type(fock_BZ_type),pointer :: fockbz
! *************************************************************************

 call timab(1502,1,tsec)

 ABI_CHECK(associated(fock),"fock must be associated")

 if (associated(fock)) then

   fockcommon=>fock%fock_common
   fockbz=> fock%fock_BZ

     invucvol=1.d0/sqrt(ucvol)
! Local variables = useful dimensions
     mband=fockcommon%mband
     mkpt=fockbz%mkpt
     mpw=dtset%mpw
     mgfft=dtset%mgfft
     ngfft=dtset%ngfft
     nkpt_bz=fockbz%nkpt_bz
     nsppol=dtset%nsppol
     nspinor=1
     ncpgr=0

! Local variables : useful arrays
     ABI_ALLOCATE(cgocc,(2,mpw))
     cgocc=zero
     ABI_ALLOCATE(cgocc_tmp,(2*mpw+1))
     cgocc_tmp=zero
     if (fockcommon%usepaw==1) then
       mgfft=dtset%mgfftdg
       ngfft=dtset%ngfftdg
       ABI_DATATYPE_ALLOCATE(cprj_tmp,(dtset%natom,nspinor))
       ABI_ALLOCATE(dimlmn,(dtset%natom))
       call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,fockcommon%pawtab,"O")
       ncpgr = 0
       if (dtset%optforces== 1) ncpgr = 3
!       if (dtset%optstress /= 0) ncpgr = 6
       call pawcprj_alloc(cprj_tmp,ncpgr,dimlmn)

       lmnmax=maxval(fockcommon%pawtab(:)%lmn_size)
       ABI_ALLOCATE(indlmn,(6,lmnmax,dtset%ntypat))
       do itypat=1,dtset%ntypat
         isize=size(fockcommon%pawtab(itypat)%indlmn,2)
         indlmn(:,1:isize,itypat)=fockcommon%pawtab(itypat)%indlmn(:,1:isize)
       end do
       ABI_ALLOCATE(indsym_,(4,dtset%nsym,dtset%natom))
       ABI_ALLOCATE(typat_srt,(dtset%natom))

       if (dtset%nsym==1) then
         indsym_=0
         do iatom=1,dtset%natom
           iatm=fockcommon%atindx(iatom)
           typat_srt(iatm)=dtset%typat(iatom)
           indsym_(4,:,iatom)=iatom
         end do
       else
         do iatom=1,dtset%natom
           iatm=fockcommon%atindx(iatom)
           typat_srt(iatm)=dtset%typat(iatom)
           indsym_(1:3,:,iatm)=indsym(1:3,:,iatom)
           indsym_(4,:,iatm)=fockcommon%atindx(indsym(4,:,iatom))
       end do
       end if
     end if

! Local variables to perform FFT
     n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
     ABI_ALLOCATE(dummytab3,(n4,n5,n6))


     if(ANY(fockbz%calc_phase(:)/=0)) then
       ABI_ALLOCATE(phase_jkpt,(2,mpw))
       phase_jkpt=zero
     end if


! =======================================================
! === Update the data relative to the occupied states ===
! =======================================================
!* The arrays cgocc_bz, kg_bz, occ_bz and npwarr_bz are already allocated with the maximal size.
!     if ((dtset%kptopt>=1).and.(dtset%kptopt<=4)) then
!       if (dtset%kptopt/=3) then

     do isppol=1,nsppol
       jbg=0 ; jcg=0 ; jkg=0 ; icp=0
       my_jsppol=isppol
       if ((isppol==2).and.(mpi_enreg%nproc_kpt/=1)) my_jsppol=1
!* Both spins are treated on the same proc., only in the case where nproc_kpt=1;
!* otherwise each proc. treats only one spin.

       ! MG: This loop is not effient!
       ! Ok the number of k-points in the BZ is usually `small` when hybrids are used
       ! but what happens if we use a 12x12x12.
       ! One should loop over the IBZ, broadcast and reconstruct the star of the k-point.
       my_jkpt=0
       do jkpt=1,nkpt_bz

         if (proc_distrb_cycle(mpi_enreg%distrb_hf,jkpt,1,mband,1,mpi_enreg%me_hf)) cycle

!* In this case, the processor does not calculate the exchange with any occupied state on jkpt.

!               if (.NOT.(proc_distrb_cycle(mpi_enreg%proc_distrb,jkpt,1,dtset%nbandhf,jsppol,mpi_enreg%me_kpt))) then
!* The state (jkpt,jband,jsppol) is stored in the array cg of this processor and copied in cgocc_tmp.
!                 icg=icg+dtset%nband(jkpt)*npwj
!               end if
!               ibg=ibg+dtset%nband(jkpt)
!               cycle
!             end if
         my_jkpt=my_jkpt+1

         ikpt=fockbz%tab_ikpt(my_jkpt)

!* ikpt = the point of IBZ that jkpt is an image of in BZ
         npwj=npwarr(ikpt)
!* npwj= number of plane wave in basis for the wavefunction
         jstwfk=fockbz%istwfk_bz(my_jkpt)
!* jstwfk= how is stored the wavefunction
         ibg=fockbz%tab_ibg(my_jkpt,my_jsppol)
!* ibg = shift to be applied on the location of data in the array occ
         icg=fockbz%tab_icg(my_jkpt,my_jsppol)
!* icg = shift to be applied on the location of data in the array cg
         icp=fockbz%tab_icp(my_jkpt,my_jsppol)
!* icp = shift to be applied on the location of data in the array cprj
         gbound_k => fockbz%gbound_bz(:,:,my_jkpt)
!* boundary of the basis sphere of G vectors
         kg_k => fockbz%kg_bz(:,1+jkg:npwj+jkg)
!* reduced plean wave coordinates
         if (fockbz%calc_phase(my_jkpt)==1) then
           phase_jkpt(:,1:npwj)=fockbz%phase(:,1+jkg:npwj+jkg)
         end if
!* phase factor at k-point j

!* Initialize the band counter
         my_jband=0
         do iband=1,dtset%nband(ikpt)


           cgocc_tmp=zero
           if (fockcommon%usepaw==1) then
             call pawcprj_set_zero(cprj_tmp)
           end if

           if(ABS(occ(iband+ibg))>tol8) then
!* If the band is occupied

!* To avoid segmentation fault, my_jband should not be greater than nbandhf
             if ((my_jband+1)>mband) then
               write(message,*) 'The number of occupied band',my_jband+1,' at k-point',&
&                 ikpt,' is greater than the value of nbandhf ', mband
               MSG_ERROR(message)
             end if

!* If the processor does not calculate the exchange with the occupied state (jkpt,my_jband), cycle
!             if (mpi_enreg%distrb_hf(jkpt,(my_jband+1),1)/=mpi_enreg%me_hf) cycle
             if (mpi_enreg%distrb_hf(jkpt,iband,1)/=mpi_enreg%me_hf) cycle
!                   if (mpi_enreg%proc_distrb(jkpt,jband,jsppol)==mpi_enreg%me_kpt) then
!* The state (jkpt,jband,jsppol) is stored in the array cg of this processor ; shift are incremented.
!                     icg=icg+npwj
!                   end if
!                   ibg=ibg+1
!* Skip the end of the loop
!                   cycle
!                 end if

!* increment the number of occupied bands treated on this processor
             my_jband = my_jband+1

!* In this case, the processor calculates the exchange with the occupied state (jkpt,my_jband).
             if (mpi_enreg%proc_distrb(ikpt,iband,isppol)==mpi_enreg%me_kpt) then
!* The state (ikpt,iband,isppol) is stored in the array cg of this processor and copied in cgocc_tmp.
               if(icg==-1) then
                 write(100,*) 'icg=-1',mpi_enreg%me,isppol,my_jsppol,jkpt,my_jkpt,ikpt,iband
               end if
             ! MG: Why packing re and im part?
               cgocc_tmp(1)=occ(iband+ibg)
               cgocc_tmp(2:npwj+1)=cg(1,1+(iband-1)*npwj+icg:iband*npwj+icg)
               cgocc_tmp(npwj+2:2*npwj+1)=cg(2,1+(iband-1)*npwj+icg:iband*npwj+icg)
               if (fockcommon%usepaw==1) then
                 call pawcprj_copy(cprj(:,icp+iband:icp+iband+nspinor-1),cprj_tmp)
               end if
             end if

!* Broadcast the state (ikpt,iband,isppol) to all the processors of comm_kpt for cgocc
             call timab(1503,1,tsec2)
             call xmpi_bcast(cgocc_tmp,mpi_enreg%proc_distrb(ikpt,iband,isppol),mpi_enreg%comm_kpt,ier)

!* Broadcast the state (ikpt,iband,isppol) to all the processors of comm_kpt for cprj
             if (fockcommon%usepaw==1) then
               call pawcprj_bcast(cprj_tmp,dtset%natom,nspinor,dimlmn,ncpgr,mpi_enreg%proc_distrb(ikpt,iband,isppol),&
&               mpi_enreg%comm_kpt,ier)
             end if
             call timab(1503,2,tsec2)
!* Keep the processors in %comm_kpt which needs the values in cgocc_tmp to build their own %cwaveocc and %occ_bz.
             if ((mpi_enreg%nproc_kpt/=1).and.(nsppol==2)) then
               if (fockbz%timerev(my_jkpt)==mpi_enreg%my_isppoltab(isppol)) cycle
!* In the case of a parallel spin-polarized calculation
!* when time reversal symmetry is applied at this k-point (timrev==1), only the processors with the opposite spin (my_isppoltab==0) are kept.
!* when time reversal symmetry is not applied at this k-point (timrev==0), only the processors with the same spin (my_isppoltab==1) are kept.

!             if (fock%timerev(my_jkpt)==1)) then
!               if (mpi_enreg%my_isppoltab(isppol)==1) cycle
!* In the case of a parallel spin-polarized calculation and when time reversal symmetry is applied at this k-point,
!* only the processors with the opposite spin are kept.
!             else
!               if (mpi_enreg%my_isppoltab(isppol)==0) cycle
!* only the processors with isppol are kept.
!             end if
             end if

!* Copy the values of cgocc_tmp in the arrays cgocc and %occ_bz
             fockbz%occ_bz(my_jband+jbg,my_jsppol) = cgocc_tmp(1)
             cgocc(1,1:npwj) = cgocc_tmp(2:npwj+1)
             cgocc(2,1:npwj) = cgocc_tmp(npwj+2:2*npwj+1)

!* calculate cg and store it in cgocc_bz
             if (fockbz%calc_phase(my_jkpt)==1) then
               do jpw=1,npwj
                 cgre=cgocc(1,jpw) ; cgim=cgocc(2,jpw)
                 cgocc(1,jpw) = phase_jkpt(1,jpw)*cgre - phase_jkpt(2,jpw)*cgim
                 cgocc(2,jpw) = phase_jkpt(1,jpw)*cgim + phase_jkpt(2,jpw)*cgre
               end do
             end if ! phase

!* apply time reversal symmetry if necessary
             if (fockbz%timerev(my_jkpt)==1) then
               cgocc(2,:) = - cgocc(2,:)
               if((mpi_enreg%nproc_kpt==1).and.(nsppol==2)) my_jsppol=mod(my_jsppol,2)+1
!* exchange spin (1 ->2 ; 2-> 1) in the sequential case.
             end if

!* apply FFT to get cwaveocc in real space

             if (allocated(fockbz%cwaveocc_bz)) then

               ABI_ALLOCATE(dummytab2,(2,npwj))
               call fourwf(1,dummytab3,cgocc(:,1:npwj),dummytab2,fockbz%cwaveocc_bz(:,:,:,:,my_jband+jbg,my_jsppol), &
&               gbound_k,gbound_k,jstwfk,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,&
&               npwj,npwj,n4,n5,n6,tim_fourwf0,0,weight1,weight1,use_gpu_cuda=dtset%use_gpu_cuda)
               ABI_DEALLOCATE(dummytab2)

             else
               fockbz%cgocc(:,jcg+1+(my_jband-1)*npwj:jcg+my_jband*npwj,my_jsppol)=cgocc(:,1:npwj)
             end if

!* calculate cprj and store it in cwaveocc_prj
             if (fockcommon%usepaw==1) then
               iband_cprj=(my_jsppol-1)*fockbz%mkptband+jbg+my_jband
               nband=1;mband0=1;iband0=1
               call pawcprj_symkn(fockbz%cwaveocc_prj(:,iband_cprj:iband_cprj+nspinor-1),cprj_tmp(:,1:nspinor),&
&               indsym_,dimlmn,iband0,indlmn,&
&               fockbz%tab_symkpt(my_jkpt),fockbz%timerev(my_jkpt),dtset%kptns(:,ikpt),fockbz%pawang%l_max-1,lmnmax,&
&               mband0,dtset%natom,nband,nspinor,dtset%nsym,dtset%ntypat,typat_srt,fockbz%pawang%zarot)

               if(dtset%optforces==1) then
                 do iatom=1,dtset%natom
                   iatm=fockcommon%atindx(iatom)
                   do ispinor=iband_cprj,iband_cprj+nspinor-1
                     do ilmn=1,fockcommon%pawtab(dtset%typat(iatom))%lmn_size
                       dcp(:)= MATMUL(TRANSPOSE(fockcommon%symrec(:,:,fockbz%tab_symkpt(my_jkpt))),&
&                                               fockbz%cwaveocc_prj(iatm,ispinor)%dcp(1,:,ilmn))
                       fockbz%cwaveocc_prj(iatm,ispinor)%dcp(1,:,ilmn)=dcp(:)
                       dcp(:)= MATMUL(TRANSPOSE(fockcommon%symrec(:,:,fockbz%tab_symkpt(my_jkpt))),&
&                                               fockbz%cwaveocc_prj(iatm,ispinor)%dcp(2,:,ilmn))
                       fockbz%cwaveocc_prj(iatm,ispinor)%dcp(2,:,ilmn)=dcp(:)
                     end do
                   end do
                 end do
               end if

             end if

           end if ! band occupied

!* update the shift to apply to occ in all case because this array is not distributed among the proc.
!               ibg=ibg+1

         end do ! iband

!* Save the true number of occupied bands in the array %nbandocc_bz
         fockbz%nbandocc_bz(my_jkpt,my_jsppol) = my_jband

!* update the shifts to apply
         jbg=jbg+my_jband
         jcg=jcg+npwj*my_jband
         jkg=jkg+npwj
       end do ! ikpt
     end do ! isppol
     if (allocated(fockbz%cwaveocc_bz)) then
       fockbz%cwaveocc_bz=fockbz%cwaveocc_bz*invucvol
     end if

     ABI_DEALLOCATE(cgocc_tmp)
     ABI_DEALLOCATE(cgocc)
     if (fockcommon%usepaw==1) then
       ABI_DEALLOCATE(indlmn)
       ABI_DEALLOCATE(indsym_)
       ABI_DEALLOCATE(typat_srt)
       ABI_DEALLOCATE(dimlmn)
       call pawcprj_free(cprj_tmp)
       ABI_DATATYPE_DEALLOCATE(cprj_tmp)
     end if
     if(allocated(phase_jkpt)) then
       ABI_DEALLOCATE(phase_jkpt)
     end if
     ABI_DEALLOCATE(dummytab3)


! Restricted or unrestricted HF
     if (nsppol==1) then
!* Update the array %occ_bz => May be limited to the occupied states only
       fockbz%occ_bz(:,:)=half*fockbz%occ_bz(:,:)

! If nsppol=1, this is a restricted Hartree-Fock calculation.
! If nsppol=2, this is an unrestricted Hartree-Fock calculation.
     end if

 end if

 call timab(1502,2,tsec)

end subroutine fock_updatecwaveocc
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_set_getghc_call
!! NAME
!!  fock_set_getghc_call
!!
!! FUNCTION
!!  Set the value of fock%getghc_call, Returns the old value
!!
!! PARENTS
!!
!! SOURCE

integer function fock_set_getghc_call(fock, new) result(old)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: new
 type(fock_common_type),intent(inout) :: fock

! *************************************************************************

 old = fock%getghc_call_
 fock%getghc_call_ = new

end function fock_set_getghc_call
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_get_getghc_call
!! NAME
!!  fock_get_getghc_call
!!
!! FUNCTION
!!  Returns the value of fock%getghc_call_
!!
!! PARENTS
!!
!! SOURCE

pure integer function fock_get_getghc_call(fock)

!Arguments ------------------------------------
!scalars
 type(fock_common_type),intent(in) :: fock

! *************************************************************************

 fock_get_getghc_call = fock%getghc_call_

end function fock_get_getghc_call
!!***

!----------------------------------------------------------------------

!!****f* m_fock/fock_print
!! NAME
!!  fock_print
!!
!! FUNCTION
!!  Print info on the fock_type data type
!!
!! INPUTS
!!  fock<crystal_t>=The object
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS"
!!  [header]=String to be printed as header for additional info.
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!      m_fock
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine fock_print(fockcommon,fockbz,header,unit,mode_paral,prtvol)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: unit,prtvol
 character(len=4),optional,intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: header
 type(fock_common_type),intent(in) :: fockcommon
 type(fock_BZ_type),intent(in) :: fockbz
!Local variables-------------------------------
 integer :: my_unt,my_prtvol
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt=std_out; if (PRESENT(unit)) my_unt=unit
 my_prtvol=0 ; if (PRESENT(prtvol)) my_prtvol=prtvol
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 msg=' ==== Info on fock_type ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 ! Important dimensions
 call wrtout(my_unt,sjoin(" my_nsppol ...",itoa(fockcommon%my_nsppol)),my_mode)
 call wrtout(my_unt,sjoin(" nkpt_bz .....",itoa(fockbz%nkpt_bz)),my_mode)

 ! Options
 call wrtout(my_unt,sjoin(" nnsclo_hf .......",itoa(fockcommon%nnsclo_hf)),my_mode)
 call wrtout(my_unt,sjoin(" ixc .............",itoa(fockcommon%ixc)),my_mode)
 call wrtout(my_unt,sjoin(" hybrid mixing....",ftoa(fockcommon%hyb_mixing)),my_mode)
 call wrtout(my_unt,sjoin(" hybrid SR mixing ",ftoa(fockcommon%hyb_mixing_sr)),my_mode)
 call wrtout(my_unt,sjoin(" hybrid range DFT ",ftoa(fockcommon%hyb_range_dft)),my_mode)
 call wrtout(my_unt,sjoin(" hybrid range Fock",ftoa(fockcommon%hyb_range_fock)),my_mode)

! write(msg,"(a,f12.1,a)")" Memory required for HF u(r) states: ",product(shape(fockbz%cwaveocc_bz)) * dp * b2Mb, " [Mb]"
! call wrtout(my_unt,msg,my_mode)

 ! Extra info.
 if (my_prtvol > 0) then
   call wrtout(my_unt,"Extra info not available",my_mode)
 end if

end subroutine fock_print
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/bare_vqg
!! NAME
!! bare_vqg
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2
!!
!! INPUTS
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  icutcoul=Option for the Coulomb potential cutoff technique
!!  divgq0= value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq. Used if q = Gamma
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero
!!  hyb_mixing=hybrid mixing coefficient for the Fock contribution
!!  hyb_mixing_sr=hybrid mixing coefficient for the short-range Fock contribution
!!  hyb_range_fock=hybrid range for separation
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  vqg(nfft)=4pi/(G+q)**2, G=0 component is set to divgq0/pi if q = Gamma.
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! PARENTS
!!      fock_getghc
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine bare_vqg(qphon,gsqcut,icutcoul,gmet,izero,hyb_mixing,hyb_mixing_sr,hyb_range_fock,nfft,nkpt_bz,ngfft,ucvol,vqg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,nfft,nkpt_bz,icutcoul
 real(dp),intent(in) :: gsqcut,hyb_mixing,hyb_mixing_sr,hyb_range_fock,ucvol
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3)
 real(dp),intent(out) ::  vqg(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1=1
 integer :: i1,i2,i23,i3,id1,id2,id3
 integer :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05
 real(dp),parameter :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,rcut,divgq0
 character(len=100) :: msg,method
!arrays
 integer :: id(3)
 real(dp),allocatable :: gq(:,:)

! *************************************************************************

 if (abs(hyb_mixing_sr)>tol8.and.abs(hyb_range_fock)<tol8) then
   msg='SR mixing<>0 while range separation=0!'
   MSG_BUG(msg)
 end if

! Re-use variable defined initially in m_vcoul
if (icutcoul == 0) method = 'SPHERE' ! Default value for the moment
if (icutcoul /= 0) method = 'unknown' ! Default value for the moment

!Treatment of the divergence at q+g=zero
!For the time being, only Spencer-Alavi scheme...
 rcut= (three*nkpt_bz*ucvol/four_pi)**(one/three)
 divgq0= two_pi*rcut**two
!divgq0=zero
!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff=gsqcut*tolfix
 vqg=zero

!Some peculiar values of q
 qeq0=0; if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qphon(1))-half)<tol12.or.abs(abs(qphon(2))-half)<tol12.or. &
&   abs(abs(qphon(3))-half)<tol12) qeq05=1
 end if

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     i23=n1*(i2-1 +(n2)*(i3-1))
     ! Do the test that eliminates the Gamma point outside of the inner loop
     ii1=1
     if (i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0) then
       ii1=2
       ! value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
       vqg(1+i23)=hyb_mixing*divgq0

!      Note the combination of Spencer-Alavi and Erfc screening
       if (abs(hyb_range_fock)>tol8)then
         vqg(1+i23)=vqg(1+i23)+hyb_mixing_sr*(pi/hyb_range_fock**2)
!        This would give a combination of Spencer-Alavi and Erfc screening,
!        unfortunately, it modifies also the tests for pure HSE06, so was not retained.
!        vqg(1+i23)=vqg(1+i23)+hyb_mixing_sr*min(divgq0,pi/(hyb_range_fock**2))
       endif

     end if

     ! Final inner loop on the first dimension (note the lower limit)
     do i1=ii1,n1
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ii=i1+i23

       if(gs<=cutoff)then
         ! Identify min/max indexes (to cancel unbalanced contributions later)
         ! Count (q+g)-vectors with similar norm
         if ((qeq05==1).and.(izero==1)) then
           ig1=i1-(i1/id1)*n1-1
           ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
           ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
           ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
         end if

         den=piinv/gs

         ! Treat the Coulomb potential cut-off by selected method
         if (abs(hyb_mixing)>tol8)then
!           SELECT CASE ( trim(method) )
!           CASE ('SPHERE')
             vqg(ii)=vqg(ii)+hyb_mixing*den*(one-cos(rcut*sqrt(four_pi/den)))
             !& vqg(ii)=vqg(ii)+hyb_mixing*den
!           CASE DEFAULT
!             msg = sjoin('Cut-off method: ',method)
!             MSG_ERROR(msg)
!           END SELECT  
         endif
!        Erfc screening
         if (abs(hyb_mixing_sr)>tol8) then
           vqg(ii)=vqg(ii)+hyb_mixing_sr*den*(one-exp(-pi/(den*hyb_range_fock**2)))
!          This other possibility combines Erfc and Spencer-Alavi screening in case rcut is too small or hyb_range_fock too large
!          if(divgq0<pi/(hyb_range_fock**2))then
!            vqg(ii)=vqg(ii)+hyb_mixing_sr*den*&
!&             (one-exp(-pi/(den*hyb_range_fock**2)))*(one-cos(rcut*sqrt(four_pi/den)))
!          endif
         endif

       end if ! Cut-off
     end do ! End loop on i1
   end do ! End loop on i2
 end do ! End loop on i3

 if (izero==1) then
   ! Set contribution of unbalanced components to zero
   if (qeq0==1) then !q=0
     call zerosym(vqg,cplex1,n1,n2,n3)
   else if (qeq05==1) then
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qphon(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qphon(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qphon(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(vqg,cplex1,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
   end if
 end if

 ABI_DEALLOCATE(gq)

end subroutine bare_vqg
!!***

!!****f* ABINIT/strfock
!!
!! NAME
!! strfock
!!
!! FUNCTION
!! Compute Fock energy contribution to stress tensor (Cartesian coordinates).
!!
!! INPUTS
!!  gsqcut=cutoff value on $G^2$ for (large) sphere inside fft box.
!!  $gsqcut=(boxcut^2)*ecut/(2._dp*(\pi^2))$
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  hyb_mixing=hybrid mixing coefficient for the Fock contribution
!!  hyb_mixing_sr=hybrid mixing coefficient for the short-range Fock contribution
!!  hyb_range_fock=hybrid range for separation
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt_bz= number of k points in the BZ
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhog2(2,nfft)= optional argument: Fourier transform of a second charge density (bohr^-3)
!!  ucvol=unit cell volume (bohr^3)
!!  vqg(nfft)=4pi/(G+q)**2
!!
!! OUTPUT
!!  fockstr(6)=components of Fock part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/bohr^3
!!   Definition of symmetric tensor storage: store 6 unique components
!!   in the order 11, 22, 33, 32, 31, 21 (suggested by Xavier Gonze).
!!
!! PARENTS
!!      fock_getghc
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine strfock(gprimd,gsqcut,fockstr,hyb_mixing,hyb_mixing_sr,hyb_range_fock,mpi_enreg,nfft,ngfft,&
&                  nkpt_bz,rhog,ucvol,qphon,&
&                 rhog2) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nkpt_bz
 real(dp),intent(in) :: gsqcut,hyb_mixing,hyb_mixing_sr,hyb_range_fock,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),rhog(2,nfft),qphon(3)
 real(dp),intent(in),optional :: rhog2(2,nfft)
 real(dp),intent(out) :: fockstr(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,id1,id2,id3,ierr,ig1,ig2,ig3,ii,irho2,me_fft,n1,n2,n3,nproc_fft
 real(dp) :: arg,cutoff,gsquar,rcut,rhogsq,tolfix=1.000000001_dp,tot,tot1,divgq0
 character(len=100) :: msg
!arrays
 real(dp) :: gcart(3),tsec(2)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 call timab(568,1,tsec)

 if (abs(hyb_mixing_sr)>tol8.and.abs(hyb_range_fock)<tol8) then
   msg='strfock: SR mixing<>0 while range separation=0!'
   MSG_BUG(msg)
 end if

 fockstr(:)=zero
 rcut= (three*nkpt_bz*ucvol/four_pi)**(one/three)
 irho2=0;if (present(rhog2)) irho2=1
 divgq0=two_pi/three*rcut**2

!Conduct looping over all fft grid points to find G vecs inside gsqcut
!Include G**2 on surface of cutoff sphere as well as inside:
 cutoff=gsqcut*tolfix
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2


 ii=0
 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2)==me_fft) then
       do i1=1,n1
         tot=zero; tot1=zero
         ig1=i1-(i1/id1)*n1-1
         ii=i1+n1*(ffti2_local(i2)-1+(n2/nproc_fft)*(i3-1))

!        Compute cartesian components of G
         gcart(1)=gprimd(1,1)*(dble(ig1)+qphon(1))+gprimd(1,2)*(dble(ig2)+qphon(2))+gprimd(1,3)*(dble(ig3)+qphon(3))
         gcart(2)=gprimd(2,1)*(dble(ig1)+qphon(1))+gprimd(2,2)*(dble(ig2)+qphon(2))+gprimd(2,3)*(dble(ig3)+qphon(3))
         gcart(3)=gprimd(3,1)*(dble(ig1)+qphon(1))+gprimd(3,2)*(dble(ig2)+qphon(2))+gprimd(3,3)*(dble(ig3)+qphon(3))
!        Compute |G+q|^2
         gsquar=gcart(1)**2+gcart(2)**2+gcart(3)**2
!        take |rho(G)|^2 for complex rhog
         if (irho2==0) then
           rhogsq=rhog(re,ii)**2+rhog(im,ii)**2
         else
           rhogsq=rhog(re,ii)*rhog2(re,ii)+rhog(im,ii)*rhog2(im,ii)
         end if
!        Case G=0:
         if(gsquar<tol10) then
           if (abs(hyb_mixing_sr)>tol8) cycle
           if (abs(hyb_mixing)>tol8) then
             fockstr(1)=fockstr(1)+hyb_mixing*divgq0*rhogsq
             fockstr(2)=fockstr(2)+hyb_mixing*divgq0*rhogsq
             fockstr(3)=fockstr(3)+hyb_mixing*divgq0*rhogsq
             cycle
           end if
         end if

!        Spencer-Alavi screening
         if (abs(hyb_mixing)>tol8) then
           arg=two_pi*rcut*sqrt(gsquar)
           tot=hyb_mixing*rhogsq*piinv/(gsquar**2)*(1-cos(arg)-arg*sin(arg)/two)
           tot1=hyb_mixing*rhogsq/three*rcut*sin(arg)/sqrt(gsquar)
         end if

!        Erfc screening
         if (abs(hyb_mixing_sr)>tol8) then
           arg=-gsquar*pi**2/(hyb_range_fock**2)
           tot=tot+hyb_mixing_sr*rhogsq*piinv/(gsquar**2)*(1.d0-exp(arg)*(1-arg))
         end if
         fockstr(1)=fockstr(1)+tot*gcart(1)*gcart(1)+tot1
         fockstr(2)=fockstr(2)+tot*gcart(2)*gcart(2)+tot1
         fockstr(3)=fockstr(3)+tot*gcart(3)*gcart(3)+tot1
         fockstr(4)=fockstr(4)+tot*gcart(3)*gcart(2)
         fockstr(5)=fockstr(5)+tot*gcart(3)*gcart(1)
         fockstr(6)=fockstr(6)+tot*gcart(2)*gcart(1)
       end do
     end if
   end do
 end do

!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
#ifdef FC_IBM
 write(std_out,*)' strfock : before mpi_comm, fockstr=',fockstr
#endif

!Init mpi_comm
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(fockstr,mpi_enreg%comm_fft ,ierr)
   call timab(48,2,tsec)
 end if

#ifdef FC_IBM
!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
 write(std_out,*)' strfock : after mpi_comm, fockstr=',fockstr
#endif

!Normalize and add term -efock/ucvol on diagonal
!efock has been set to zero because it is not yet known. It will be added later.
 fockstr(1)=-fockstr(1)
 fockstr(2)=-fockstr(2)
 fockstr(3)=-fockstr(3)
 fockstr(4)=-fockstr(4)
 fockstr(5)=-fockstr(5)
 fockstr(6)=-fockstr(6)

 call timab(568,2,tsec)

end subroutine strfock
!!***

end module m_fock
!!***
