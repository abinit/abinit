!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2021 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
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

module m_orbmag

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : projbd
  use m_fft,              only : fftpac
  use m_fftcore,          only : kpgsph
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type
  use m_kg,               only : getph,mkkin,mkkpg,ph1d3d
  use m_kpts,             only : listkk, smpbz
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle
  use m_nonlop,           only : nonlop
  use m_pawang,           only : pawang_type
  use m_pawfgr,           only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_pawrad,           only : pawrad_type,pawrad_deducer0,simp_gen
  use m_paw_sphharm,      only : setsym_ylm,slxyzs
  use m_pawtab,           only : pawtab_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free,&
       &                         pawcprj_put, pawcprj_getdim 
  use m_spacepar,         only : make_vectornd
  use m_symtk,            only : symatm
  use m_time,             only : timab

  implicit none

  private
!!***

!!****t* m_orbmag/orbmag_type
!! NAME
!! orbmag_type
!!
!! FUNCTION
!! variables used in orbital magnetism calculation
!!
!! SOURCE

  type, public :: orbmag_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
     integer :: orbmag              ! value of orbmag input variable in use
     integer :: fmkmem              ! number of k-points in the FBZ per cpu
     integer :: fmkmem_max          ! max of fmkmem
     integer :: fnkpt               ! number of k-points in the FBZ
     integer :: lmax
     integer :: lmnmax
     integer :: lmn2max
     integer :: mkmem_max           ! max of mkmem
     integer :: natom               ! number of atoms in unit cell
     integer :: my_natom            ! number of atoms treated by current proc
     integer :: mband_occ           ! max number of occupied bands (over spin)
     ! this number must be the same for every k
     integer :: nspinor             ! nspinor input from data set
     integer :: nsym
     integer :: usepaw              ! 1 if a PAW calculation, 0 else

     ! Real(dp) scalars
     real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1

     ! Real(dp) arrays
     real(dp) :: chern(2,3)           ! result of chern number calculation

     real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-point and its nearest neighbour along idir

     real(dp) :: orbmagvec(2,3)     ! result of orbital magnetization calculation

     ! Integer pointers
     integer, allocatable :: atom_indsym(:,:,:) ! atom_indsym(4,nsym,natom)
     ! this is data on how the symmetries map the atoms in the cell
     ! see symatm.F90 for full description
     integer, allocatable :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
     ! for each k-point, stores the location
     ! of the WF in the cg array
     integer, allocatable :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
     ! for each k-point, stores the location
     ! of the cprj in the cprj array (used only
     ! for PAW calculations)
     integer, allocatable :: fkgindex(:)     ! same as kgindex, but defined
     ! for the FBZ and intended to use
     ! with pwindf
     integer, allocatable :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
     ! ikpt_dp(ikpt,ii,idir) = index of the
     ! k-point at k+dk (ii=1) and k-dk (ii=2)
     integer, allocatable :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtorbmag%fnkpt,1:6)
     ! information needed to fold a
     ! k-point in the FBZ into the IBZ;
     ! the second index (1:6)
     ! is as described in listkk
     integer, allocatable :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
     ! k-points in the FBZ k-point list

     integer, allocatable :: kg(:,:) ! reduced (integer) coordinates of G vecs in basis sphere

     integer, allocatable :: kgindex(:)      ! kgind(nkpt) on current processor
     ! kgindex(ikpt) = ikg

     integer, allocatable :: lmn_size(:)        ! lmn_size(ntypat)
     integer, allocatable :: lmn2_size(:)       ! lmn2_size(ntypat)

     integer, allocatable :: nband_occ(:)       ! nband_occ(nsppol) = actual number of occupied bands
     !  can be different for spin up and down!!!
     ! Real(dp) allocatables

     real(dp), allocatable :: fkptns(:,:)       ! fkptns(3,1:dtorbmag%fnkpt) k-points in FBZ

     real(dp), allocatable :: zarot(:,:,:,:)
     !  zarot(l_size_max,l_size_max,l_max,nsym)
     !  Coeffs of the transformation of real spherical
     !  harmonics under the symmetry operations. These are needed when the
     ! cprj's need to be computed in the full BZ, that is,
     ! in the PAW case with kptopt /= 3.

     ! complex(dpc) allocatable

  end type orbmag_type


  ! Bound methods:
  public :: destroy_orbmag
  public :: initorbmag
  public :: orbmag_ddk

  private :: make_onsite_l_k_n
  private :: make_onsite_bm_k_n
  private :: make_rhorij1_k_n
  private :: make_S1trace_k_n
  private :: orbmag_ddk_output
  
CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_orbmag/destroy_orbmag
!! NAME
!!
!! FUNCTION
!!   deallocate fields in orbmag structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_orbmag(dtorbmag)

  !Arguments ------------------------------------
  !array
  type(orbmag_type),intent(inout) :: dtorbmag

  ! ************************************************************************

  ! Integer pointers
  if(allocated(dtorbmag%atom_indsym))  then
     ABI_FREE(dtorbmag%atom_indsym)
  end if
  if(allocated(dtorbmag%cgindex))  then
     ABI_FREE(dtorbmag%cgindex)
  end if
  if(allocated(dtorbmag%cprjindex))  then
     ABI_FREE(dtorbmag%cprjindex)
  end if
  if(allocated(dtorbmag%fkgindex))  then
     ABI_FREE(dtorbmag%fkgindex)
  end if
  if(allocated(dtorbmag%ikpt_dk))  then
     ABI_FREE(dtorbmag%ikpt_dk)
  end if
  if(allocated(dtorbmag%indkk_f2ibz))  then
     ABI_FREE(dtorbmag%indkk_f2ibz)
  end if
  if(allocated(dtorbmag%i2fbz))  then
     ABI_FREE(dtorbmag%i2fbz)
  end if
  if(allocated(dtorbmag%kg)) then
     ABI_FREE(dtorbmag%kg)
  end if
  if(allocated(dtorbmag%kgindex))  then
     ABI_FREE(dtorbmag%kgindex)
  end if
  if(allocated(dtorbmag%lmn_size))  then
     ABI_FREE(dtorbmag%lmn_size)
  end if
  if(allocated(dtorbmag%lmn2_size))  then
     ABI_FREE(dtorbmag%lmn2_size)
  end if
  if(allocated(dtorbmag%nband_occ))  then
     ABI_FREE(dtorbmag%nband_occ)
  end if
  ! Real(dp) pointers

  if(allocated(dtorbmag%fkptns))  then
     ABI_FREE(dtorbmag%fkptns)
  end if
  if(allocated(dtorbmag%zarot))  then
     ABI_FREE(dtorbmag%zarot)
  end if

end subroutine destroy_orbmag
!!***

!!****f* ABINIT/initorbmag
!! NAME
!! initorbmag
!!
!! FUNCTION
!! Initialization of orbital magnetization calculation; similar to initberry
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group.
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! SIDE EFFECTS
!!  mpi_enreg = information about MPI parallelization
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      kpgsph,listkk,setsym_ylm,smpbz,symatm,timab,wrtout,xmpi_max,xmpi_sum
!!
!! SOURCE

subroutine initorbmag(dtorbmag,dtset,gmet,gprimd,kg,mpi_enreg,npwarr,occ,&
     &                     pawtab,psps,pwind,pwind_alloc,pwnsfac,&
     &                     rprimd,symrec,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(out) :: pwind_alloc
  type(MPI_type),intent(inout) :: mpi_enreg
  type(dataset_type),intent(inout) :: dtset
  type(orbmag_type),intent(out) :: dtorbmag
  type(pseudopotential_type),intent(in) :: psps
  !arrays
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  integer,pointer :: pwind(:,:,:)
  real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
  real(dp),pointer :: pwnsfac(:,:)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables-------------------------------
  !scalars
  integer :: brav,exchn2n3d,fnkpt_computed
  integer :: iband,icg,icprj,idir,idum,idum1,ierr,ifor,ikg,ikg1
  integer :: ikpt,ikpt_loc,ikpti,ikpt1,ikpt1f,ikpt1i
  integer :: index,ipw,ipwnsfac,isign,isppol,istwf_k,isym,isym1,itrs,itypat
  integer :: jpw,lmax,lmn2_size_max
  integer :: mband_occ_k,me,me_g0,mkmem_,mkpt,my_nspinor,nband_k,nkptlatt,nproc,npw_k,npw_k1
  integer :: option,spaceComm
  real(dp) :: diffk1,diffk2,diffk3,ecut_eff
  real(dp) :: kpt_shifted1,kpt_shifted2,kpt_shifted3,rdum
  character(len=500) :: message
  !arrays
  integer :: iadum(3),iadum1(3),dg(3)
  integer,allocatable :: kg1_k(:,:)
  real(dp) :: diffk(3),dk(3),dum33(3,3),kpt1(3),tsec(2)
  real(dp),allocatable :: spkpt(:,:)

  ! *************************************************************************

  DBG_ENTER("COLL")

  call timab(1001,1,tsec)
  call timab(1002,1,tsec)

  !save the current value of nspinor
  dtorbmag%nspinor = dtset%nspinor

  !----------------------------------------------------------------------------
  !-------------------- Obtain k-point grid in the full BZ --------------------
  !----------------------------------------------------------------------------

  if(dtset%kptopt==1 .or. dtset%kptopt==2 .or. dtset%kptopt==4)then
     !  Compute the number of k points in the G-space unit cell
     nkptlatt=dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
          &   +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
          &   +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
          &   -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
          &   -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
          &   -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)

     !  Call smpbz to obtain the list of k-point in the full BZ - without symmetry reduction
     option = 0
     brav = 1
     mkpt=nkptlatt*dtset%nshiftk
     ABI_MALLOC(spkpt,(3,mkpt))
     call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
     dtorbmag%fnkpt = fnkpt_computed
     ABI_MALLOC(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
     dtorbmag%fkptns(:,:)=spkpt(:,1:dtorbmag%fnkpt)
     ABI_FREE(spkpt)
  else if(dtset%kptopt==3.or.dtset%kptopt==0)then
     dtorbmag%fnkpt=dtset%nkpt
     ABI_MALLOC(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
     dtorbmag%fkptns(1:3,1:dtorbmag%fnkpt)=dtset%kpt(1:3,1:dtorbmag%fnkpt)
     if(dtset%kptopt==0)then
        write(message,'(10a)') ch10,&
             &     ' initorbmag : WARNING -',ch10,&
             &     '  you have defined manually the k-point grid with kptopt = 0',ch10,&
             &     '  the orbital magnetization calculation works only with a regular k-points grid,',ch10,&
             &     '  abinit doesn''t check if your grid is regular...'
        call wrtout(std_out,message,'PERS')
     end if
  end if

  !call listkk to get mapping from FBZ to IBZ
  rdum=1.0d-5  ! cutoff distance to decide when two k points match
  ABI_MALLOC(dtorbmag%indkk_f2ibz,(dtorbmag%fnkpt,6))

  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  spaceComm=mpi_enreg%comm_cell

  !JWZ: The following may need modification in the future
  !**** no spin-polarization doubling ; do not allow use of time reversal symmetry ****

  call timab(1002,2,tsec)
  call timab(1003,1,tsec)

  call listkk(rdum,gmet,dtorbmag%indkk_f2ibz,dtset%kptns,dtorbmag%fkptns,dtset%nkpt,&
       & dtorbmag%fnkpt,dtset%nsym,1,dtset%symafm,symrec,0,spaceComm,use_symrec=.True.)

  call timab(1003,2,tsec)
  call timab(1004,1,tsec)

  !Construct i2fbz and f2ibz
  ABI_MALLOC(dtorbmag%i2fbz,(dtset%nkpt))
  idum=0
  do ikpt=1,dtorbmag%fnkpt
     if (dtorbmag%indkk_f2ibz(ikpt,2)==1 .and. &
          &   dtorbmag%indkk_f2ibz(ikpt,6) == 0 .and. &
          &   maxval(abs(dtorbmag%indkk_f2ibz(ikpt,3:5))) == 0 ) then
        dtorbmag%i2fbz(dtorbmag%indkk_f2ibz(ikpt,1))=ikpt
        idum=idum+1
     end if
  end do
  if (idum/=dtset%nkpt)then
     message = ' Found wrong number of k-points in IBZ'
     ABI_ERROR(message)
  end if

  !----------------------------------------------------------------------------
  !------------- Allocate PAW space as necessary ------------------------------
  !----------------------------------------------------------------------------

  dtorbmag%usepaw   = psps%usepaw
  dtorbmag%natom    = dtset%natom
  dtorbmag%my_natom = mpi_enreg%my_natom

  ABI_MALLOC(dtorbmag%lmn_size,(dtset%ntypat))
  ABI_MALLOC(dtorbmag%lmn2_size,(dtset%ntypat))
  do itypat = 1, dtset%ntypat
     dtorbmag%lmn_size(itypat) = pawtab(itypat)%lmn_size
     dtorbmag%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
  end do

  lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
  dtorbmag%lmn2max = lmn2_size_max

  ABI_MALLOC(dtorbmag%cprjindex,(dtset%nkpt,dtset%nsppol))
  dtorbmag%cprjindex(:,:) = 0

  if (dtset%kptopt /= 3) then
     ABI_MALLOC(dtorbmag%atom_indsym,(4,dtset%nsym,dtorbmag%natom))
     call symatm(dtorbmag%atom_indsym,dtorbmag%natom,dtset%nsym,symrec,dtset%tnons,tol8,dtset%typat,xred)
     lmax = psps%mpsang - 1
     ABI_MALLOC(dtorbmag%zarot,(2*lmax+1,2*lmax+1,lmax+1,dtset%nsym))
     call setsym_ylm(gprimd,lmax,dtset%nsym,1,rprimd,symrec,dtorbmag%zarot)
     dtorbmag%nsym = dtset%nsym
     dtorbmag%lmax = lmax
     dtorbmag%lmnmax = psps%lmnmax
  end if

  ! !------------------------------------------------------------------------------
  ! !------------------- Compute variables related to MPI // ----------------------
  ! !------------------------------------------------------------------------------
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me=xmpi_comm_rank(spaceComm)

  if (nproc==1) then
     dtorbmag%fmkmem = dtorbmag%fnkpt
     dtorbmag%fmkmem_max = dtorbmag%fnkpt
     dtorbmag%mkmem_max = dtset%nkpt
  else
     dtorbmag%fmkmem = 0
     do ikpt = 1, dtorbmag%fnkpt
        ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
        nband_k = dtset%nband(ikpti)
        if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,-1,me))) &
             &     dtorbmag%fmkmem = dtorbmag%fmkmem + 1
     end do
     !  Maximum value of mkmem and fmkmem
     call xmpi_max(dtorbmag%fmkmem,dtorbmag%fmkmem_max,spaceComm,ierr)
     !  I have to use the dummy variable mkmem_ because
     !  mkmem is declared as intent(in) while the first
     !  argument of xmpi_max must be intent(inout)
     mkmem_ = dtset%mkmem
     call xmpi_max(mkmem_,dtorbmag%mkmem_max,spaceComm,ierr)
  end if

  ABI_MALLOC(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtorbmag%fmkmem_max*dtset%nsppol, 1:2))
  ABI_MALLOC(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtorbmag%mkmem_max*dtset%nsppol, 1:2))
  ABI_MALLOC(mpi_enreg%kptdstrb,(nproc,6,dtorbmag%fmkmem_max*dtset%nsppol*2))
  ABI_MALLOC(mpi_enreg%mkmem,(0:nproc-1))
  mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
  mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
  mpi_enreg%kptdstrb(:,:,:)       = 0
  mpi_enreg%mkmem(:)              = 0

  pwind_alloc = dtset%mpw*dtorbmag%fmkmem_max

  ABI_MALLOC(pwind,(pwind_alloc,2,3))
  ABI_MALLOC(pwnsfac,(2,pwind_alloc))

  ! !------------------------------------------------------------------------------
  ! !---------------------- Compute orbmag_type variables -------------------------
  ! !------------------------------------------------------------------------------

  !Initialization of orbmag_type variables
  dtorbmag%dkvecs(:,:) = zero
  ABI_MALLOC(dtorbmag%ikpt_dk,(dtorbmag%fnkpt,2,3))
  ABI_MALLOC(dtorbmag%cgindex,(dtset%nkpt,dtset%nsppol))
  ABI_MALLOC(dtorbmag%kgindex,(dtset%nkpt))
  ABI_MALLOC(dtorbmag%fkgindex,(dtorbmag%fnkpt))
  dtorbmag%ikpt_dk(:,:,:) = 0
  dtorbmag%cgindex(:,:) = 0
  dtorbmag%mband_occ = 0
  ABI_MALLOC(dtorbmag%nband_occ,(dtset%nsppol))
  dtorbmag%kgindex(:) = 0
  dtorbmag%fkgindex(:) = 0
  ABI_MALLOC(dtorbmag%kg,(3,dtset%mpw*dtset%mkmem))
  dtorbmag%kg(:,:) = kg(:,:)

  !Compute spin degeneracy
  if (dtset%nsppol == 1 .and. dtset%nspinor == 1) then
     dtorbmag%sdeg = two
  else if (dtset%nsppol == 2 .or. my_nspinor == 2) then
     dtorbmag%sdeg = one
  end if

  !Compute the number of occupied bands and check that
  !it is the same for each k-point

  index = 0
  do isppol = 1, dtset%nsppol
     dtorbmag%nband_occ(isppol) = 0
     do ikpt = 1, dtset%nkpt

        mband_occ_k = 0
        nband_k = dtset%nband(ikpt + (isppol - 1)*dtset%nkpt)

        do iband = 1, nband_k
           index = index + 1
           if (abs(occ(index) - dtorbmag%sdeg) < tol8) mband_occ_k = mband_occ_k + 1
        end do

        if (ikpt > 1) then
           if (dtorbmag%nband_occ(isppol) /= mband_occ_k) then
              message = "The number of valence bands is not the same for every k-point of present spin channel"
              ABI_ERROR(message)
           end if
        else
           dtorbmag%mband_occ         = max(dtorbmag%mband_occ, mband_occ_k)
           dtorbmag%nband_occ(isppol) = mband_occ_k
        end if

     end do                ! close loop over ikpt
  end do                ! close loop over isppol

  !Compute the location of each wavefunction

  icg = 0
  icprj = 0
  !ikg = 0
  do isppol = 1, dtset%nsppol
     do ikpt = 1, dtset%nkpt

        nband_k = dtset%nband(ikpt + (isppol-1)*dtset%nkpt)

        if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) cycle

        dtorbmag%cgindex(ikpt,isppol) = icg
        npw_k = npwarr(ikpt)
        icg = icg + npw_k*dtorbmag%nspinor*nband_k

        if (psps%usepaw == 1) then
           dtorbmag%cprjindex(ikpt,isppol) = icprj
           icprj = icprj + dtorbmag%nspinor*nband_k
        end if

     end do
  end do

  ikg = 0
  do ikpt = 1, dtset%nkpt
     if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,1,me)).and.&
          &   (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,dtset%nsppol,me))) cycle

     npw_k = npwarr(ikpt)
     dtorbmag%kgindex(ikpt) = ikg
     ikg = ikg + npw_k
  end do

  call timab(1004,2,tsec)

  !------------------------------------------------------------------------------
  !---------------------- Compute dk --------------------------------------------
  !------------------------------------------------------------------------------

  call timab(1005,1,tsec)

  do idir = 1, 3

     !    Compute dk(:), the vector between a k-point and its nearest
     !    neighbour along the direction idir

     dk(:) = zero
     dk(idir) = 1._dp   ! 1 mean there is no other k-point un the direction idir
     do ikpt = 2, dtorbmag%fnkpt
        diffk(:) = abs(dtorbmag%fkptns(:,ikpt) - dtorbmag%fkptns(:,1))
        if ((diffk(1) < dk(1)+tol8).and.(diffk(2) < dk(2)+tol8).and.&
             &     (diffk(3) < dk(3)+tol8)) dk(:) = diffk(:)
     end do
     dtorbmag%dkvecs(:,idir) = dk(:)
     !    DEBUG
     !    write(std_out,*)' initorbmag : idir, dk', idir, dk
     !    ENDDEBUG

     !    For each k point, find k_prim such that k_prim= k + dk mod(G)
     !    where G is a vector of the reciprocal lattice

     do ikpt = 1, dtorbmag%fnkpt

        !      First k+dk, then k-dk
        do isign=-1,1,2
           kpt_shifted1=dtorbmag%fkptns(1,ikpt)- isign*dk(1)
           kpt_shifted2=dtorbmag%fkptns(2,ikpt)- isign*dk(2)
           kpt_shifted3=dtorbmag%fkptns(3,ikpt)- isign*dk(3)
           !        Note that this is still a order fnkpt**2 algorithm.
           !        It is possible to implement a order fnkpt algorithm, see listkk.F90.
           do ikpt1 = 1, dtorbmag%fnkpt
              diffk1=dtorbmag%fkptns(1,ikpt1) - kpt_shifted1
              if(abs(diffk1-nint(diffk1))>tol8)cycle
              diffk2=dtorbmag%fkptns(2,ikpt1) - kpt_shifted2
              if(abs(diffk2-nint(diffk2))>tol8)cycle
              diffk3=dtorbmag%fkptns(3,ikpt1) - kpt_shifted3
              if(abs(diffk3-nint(diffk3))>tol8)cycle
              dtorbmag%ikpt_dk(ikpt,(isign+3)/2,idir) = ikpt1
              exit
           end do   ! ikpt1
        end do     ! isign

     end do     ! ikpt

  end do     ! close loop over idir

  call timab(1005,2,tsec)
  call timab(1006,1,tsec)

  !------------------------------------------------------------------------------
  !------------ Build the array pwind that is needed to compute the -------------
  !------------ overlap matrices at k +- dk                         -------------
  !------------------------------------------------------------------------------

  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
  pwind(:,:,:) = 0
  pwnsfac(1,:) = 1.0_dp
  pwnsfac(2,:) = 0.0_dp
  ABI_MALLOC(kg1_k,(3,dtset%mpw))

  ipwnsfac = 0

  do idir = 1, 3

     dk(:) = dtorbmag%dkvecs(:,idir)

     do ifor = 1, 2

        if (ifor == 2) dk(:) = -1._dp*dk(:)

        !      Build pwind and kgindex
        !      NOTE: The array kgindex is important for parallel execution.
        !      In case nsppol = 2, it may happen that a particular processor
        !      treats k-points at different spin polarizations.
        !      In this case, it is not possible to address the elements of
        !      pwind correctly without making use of the kgindex array.

        ikg = 0 ; ikpt_loc = 0 ; isppol = 1
        do ikpt = 1, dtorbmag%fnkpt

           ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
           nband_k = dtset%nband(ikpti)
           ikpt1f = dtorbmag%ikpt_dk(ikpt,ifor,idir)
           ikpt1i = dtorbmag%indkk_f2ibz(ikpt1f,1)

           if ((proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,1,me)).and.&
                &       (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,dtset%nsppol,me))) cycle

           ikpt_loc = ikpt_loc + 1

           !        Build basis sphere of plane waves for the nearest neighbour of
           !        the k-point (important for MPI //)

           kg1_k(:,:) = 0
           kpt1(:) = dtset%kptns(:,ikpt1i)
           call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,&
                &       1,mpi_enreg,dtset%mpw,npw_k1)
           me_g0=mpi_enreg%me_g0


           !        ji: fkgindex is defined here !
           dtorbmag%fkgindex(ikpt) = ikg

           !
           !        Deal with symmetry transformations
           !

           !        bra k-point k(b) and IBZ k-point kIBZ(b) related by
           !        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
           !        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
           !
           !        For the ket k-point:
           !        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
           !        where GBZ(k) takes k(k) to the BZ
           !

           isym  = dtorbmag%indkk_f2ibz(ikpt,2)
           isym1 = dtorbmag%indkk_f2ibz(ikpt1f,2)

           !        Construct transformed G vector that enters the matching condition:
           !        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

           dg(:) = -dtorbmag%indkk_f2ibz(ikpt,3:5) &
                &       -nint(-dtorbmag%fkptns(:,ikpt) - dk(:) - tol10 &
                &       +dtorbmag%fkptns(:,ikpt1f)) &
                &       +dtorbmag%indkk_f2ibz(ikpt1f,3:5)

           iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

           dg(:) = iadum(:)

           if ( dtorbmag%indkk_f2ibz(ikpt1f,6) == 1 ) dg(:) = -dg(:)

           !        Construct S(k)^{t,-1} S(b)^{t}

           dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

           !        Construct alpha(k) alpha(b)

           if (dtorbmag%indkk_f2ibz(ikpt,6) == dtorbmag%indkk_f2ibz(ikpt1f,6)) then
              itrs=0
           else
              itrs=1
           end if


           npw_k  = npwarr(ikpti)
           !        npw_k1 = npwarr(ikpt1i)

           !        loop over bra G vectors
           do ipw = 1, npw_k

              !          NOTE: the bra G vector is taken for the sym-related IBZ k point,
              !          not for the FBZ k point
              iadum(:) = kg(:,dtorbmag%kgindex(ikpti) + ipw)

              !          Store non-symmorphic operation phase factor exp[i2\pi \alpha G \cdot t]

              if ( ipwnsfac == 0 ) then
                 rdum=0.0_dp
                 do idum=1,3
                    rdum=rdum+dble(iadum(idum))*dtset%tnons(idum,isym)
                 end do
                 rdum=two_pi*rdum
                 if ( dtorbmag%indkk_f2ibz(ikpt,6) == 1 ) rdum=-rdum
                 pwnsfac(1,ikg+ipw) = cos(rdum)
                 pwnsfac(2,ikg+ipw) = sin(rdum)
              end if

              !          to determine r.l.v. matchings, we transformed the bra vector
              !          Rotation
              iadum1(:)=0
              do idum1=1,3
                 iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
              end do
              iadum(:)=iadum1(:)
              !          Time reversal
              if (itrs==1) iadum(:)=-iadum(:)
              !          Translation
              iadum(:) = iadum(:) + dg(:)

              do jpw = 1, npw_k1
                 iadum1(1:3) = kg1_k(1:3,jpw)
                 if ( (iadum(1) == iadum1(1)).and. &
                      &           (iadum(2) == iadum1(2)).and. &
                      &           (iadum(3) == iadum1(3)) ) then
                    pwind(ikg + ipw,ifor,idir) = jpw
                    !              write(std_out,'(a,2x,3i4,2x,i4)') 'Found !:',iadum1(:),jpw
                    exit
                 end if
              end do
           end do

           ikg  = ikg + npw_k

        end do    ! close loop over ikpt

        ipwnsfac = 1

     end do    ! close loop over ifor

  end do        ! close loop over idir


  call timab(1008,2,tsec)
  call timab(1009,1,tsec)

  !Build mpi_enreg%kptdstrb
  !array required to communicate the WFs between cpus
  !(MPI // over k-points)
  if (nproc>1) then
     do idir = 1, 3
        do ifor = 1, 2

           ikpt_loc = 0
           do isppol = 1, dtset%nsppol

              do ikpt = 1, dtorbmag%fnkpt

                 ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
                 nband_k = dtset%nband(ikpti)
                 ikpt1f = dtorbmag%ikpt_dk(ikpt,ifor,idir)
                 ikpt1i = dtorbmag%indkk_f2ibz(ikpt1f,1)

                 if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle

                 ikpt_loc = ikpt_loc + 1
                 mpi_enreg%kptdstrb(me + 1,ifor+2*(idir-1),ikpt_loc) = &
                      &           ikpt1i + (isppol - 1)*dtset%nkpt

                 mpi_enreg%kptdstrb(me+1,ifor+2*(idir-1),ikpt_loc+dtorbmag%fmkmem_max*dtset%nsppol) = &
                      &           ikpt1f + (isppol - 1)*dtorbmag%fnkpt

              end do   ! ikpt
           end do     ! isppol
        end do       ! ifor
     end do           ! idir
  end if             ! nproc>1

  !build mpi_enreg%kpt_loc2fbz_sp
  ikpt_loc = 0
  do isppol = 1, dtset%nsppol
     do ikpt = 1, dtorbmag%fnkpt

        ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
        nband_k = dtset%nband(ikpti)

        if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpti,1,nband_k,isppol,me)) cycle

        ikpt_loc = ikpt_loc + 1

        mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 1) = ikpt
        mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc, 2) = isppol

     end do
  end do

  !should be temporary
  !unassigned mpi_enreg%kpt_loc2fbz_sp are empty ; inform other cpu (there are better ways...)
  mpi_enreg%mkmem(me) = dtset%mkmem
  !do ii=ikpt_loc+1,dtefield%fmkmem_max
  !mpi_enreg%kpt_loc2fbz_sp(me, ii, 1) = -1
  !end do

  call xmpi_sum(mpi_enreg%kptdstrb,spaceComm,ierr)
  call xmpi_sum(mpi_enreg%kpt_loc2fbz_sp,spaceComm,ierr)

  ABI_FREE(kg1_k)


  call timab(1009,2,tsec)
  call timab(1001,2,tsec)

  DBG_EXIT("COLL")

end subroutine initorbmag
!!***


!----------------------------------------------------------------------

!!****f* ABINIT/make_onsite_l_k_n
!! NAME
!! make_onsite_l_k_n
!!
!! FUNCTION
!! Compute 1/2 <L_R> onsite contribution to orbital magnetization at given k point, band, and idir
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine make_onsite_l_k_n(cprj_k,dtset,iband,idir,nband_k,onsite_l_k_n,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: iband,idir,nband_k
  complex(dpc),intent(out) :: onsite_l_k_n
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iatom,ilmn,il,im,itypat,jlmn,jl,jm,klmn,kln,mesh_size
  real(dp) :: intg
  complex(dpc) :: cpb,cpk,orbl_me

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------

  onsite_l_k_n = czero
  do iatom=1,dtset%natom
    itypat=dtset%typat(iatom)
    mesh_size=pawtab(itypat)%mesh_size
    ABI_MALLOC(ff,(mesh_size))
    do jlmn=1,pawtab(itypat)%lmn_size
       jl=pawtab(itypat)%indlmn(1,jlmn)
       jm=pawtab(itypat)%indlmn(2,jlmn)
       do ilmn=1,pawtab(itypat)%lmn_size
          il=pawtab(itypat)%indlmn(1,ilmn)
          im=pawtab(itypat)%indlmn(2,ilmn)
          klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
          kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
          ! compute <L_dir>
          call slxyzs(il,im,idir,jl,jm,orbl_me)
          ! compute integral of phi_i*phi_j - tphi_i*tphi_j
          if (abs(orbl_me) > tol8) then
             ff(1:mesh_size)=pawtab(itypat)%phiphj(1:mesh_size,kln) - pawtab(itypat)%tphitphj(1:mesh_size,kln)
             call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
             call simp_gen(intg,ff,pawrad(itypat))
             cpb=cmplx(cprj_k(iatom,iband)%cp(1,ilmn),cprj_k(iatom,iband)%cp(2,ilmn),KIND=dpc)
             cpk=cmplx(cprj_k(iatom,iband)%cp(1,jlmn),cprj_k(iatom,iband)%cp(2,jlmn),KIND=dpc)
             onsite_l_k_n=onsite_l_k_n+conjg(cpb)*half*orbl_me*intg*cpk
          end if ! end check that |L_dir| > 0, otherwise ignore term
       end do ! end loop over ilmn
    end do ! end loop over jlmn
    ABI_FREE(ff)
  end do ! end loop over atoms

end subroutine make_onsite_l_k_n
!!***

!!****f* ABINIT/make_onsite_bm_k_n
!! NAME
!! make_onsite_bm_k_n
!!
!! FUNCTION
!! Compute A_0.A_N onsite term for magnetic field + nuclear magnetic dipole moment
!! for k pt and 1 band
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine make_onsite_bm_k_n(cprj_k,dtset,iband,idir,nband_k,onsite_bm_k_n,&
     & pawang,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,iband,nband_k
  complex(dpc),intent(out) :: onsite_bm_k_n
  type(pawang_type),intent(in) :: pawang
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: gint,iatom,il,im,ilmn,itypat
  integer :: jl,jm,jlmn,klmn,klm,kln,lpmp,mesh_size
  real(dp) :: bm1,bm2,d00,d20,d22,dij,intg,scale_conversion
  complex(dpc) :: cpb,cpk

  !arrays
  real(dp),allocatable :: ff(:)

  ! ***********************************************************************

  ! this term can only be non-zero if some nucdipmom is nonzero
  scale_conversion = half*FineStructureConstant2
  d00 = sqrt(4.0*pi)/3.0
  dij = sqrt(4.0*pi/15.0)
  d20 = sqrt(16.0*pi/5.0)/6.0
  d22 = sqrt(16.0*pi/15.0)/2.0
  onsite_bm_k_n = czero

  do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     mesh_size=pawtab(itypat)%mesh_size
     ABI_MALLOC(ff,(mesh_size))
     do jlmn=1,pawtab(itypat)%lmn_size
        jl=pawtab(itypat)%indlmn(1,jlmn)
        jm=pawtab(itypat)%indlmn(2,jlmn)
        do ilmn=1,pawtab(itypat)%lmn_size
           il=pawtab(itypat)%indlmn(1,ilmn)
           im=pawtab(itypat)%indlmn(2,ilmn)
           klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
           kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
           klm = pawtab(itypat)%indklmn(1,klmn) ! need this for bm2 gaunt integral selection
           ! compute integral of (phi_i*phi_j - tphi_i*tphi_j)/r
           ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln) - &
                &           pawtab(itypat)%tphitphj(2:mesh_size,kln)) / &
                &           pawrad(itypat)%rad(2:mesh_size)
           call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
           call simp_gen(intg,ff,pawrad(itypat))
           ! term B.m r^2/r^3
           bm1=zero
           if ( (jl .EQ. il) .AND. (jm .EQ. im) .AND. (abs(dtset%nucdipmom(idir,iatom)) .GT. tol8) ) then
              bm1 = scale_conversion*dtset%nucdipmom(idir,iatom)*intg
           end if
           bm2 = zero
           ! xx, yy, zz cases all have the same contribution from S00
           lpmp=1
           gint = pawang%gntselect(lpmp,klm)
           if (gint > 0) then
              bm2=bm2+scale_conversion*dtset%nucdipmom(idir,iatom)*d00*pawang%realgnt(gint)*intg
           end if
           ! all other contributions involve Gaunt integrals of S_{2m}
           do lpmp = 5, 9
              gint = pawang%gntselect(lpmp,klm)
              if (gint > 0) then
                 select case (lpmp)
                 case (5) ! S_{2,-2} contributes to xy term
                    select case (idir)
                    case (1)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(2,iatom)*dij*pawang%realgnt(gint)*intg
                    case (2)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(1,iatom)*dij*pawang%realgnt(gint)*intg
                    end select
                 case (6) ! S_{2,-1} contributes to yz term
                    select case (idir)
                    case (2)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(3,iatom)*dij*pawang%realgnt(gint)*intg
                    case (3)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(2,iatom)*dij*pawang%realgnt(gint)*intg
                    end select
                 case (7) ! S_{2,0} contributes to xx, yy, and zz terms
                    select case (idir)
                       case (1)
                          bm2=bm2-scale_conversion*dtset%nucdipmom(1,iatom)*d20*pawang%realgnt(gint)*intg
                       case (2)
                          bm2=bm2-scale_conversion*dtset%nucdipmom(2,iatom)*d20*pawang%realgnt(gint)*intg
                       case (3)
                          bm2=bm2+scale_conversion*dtset%nucdipmom(3,iatom)*2.0*d20*pawang%realgnt(gint)*intg
                       end select
                 case (8) ! S_{2,+1} contributes to xz term
                    select case (idir)
                    case (1)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(3,iatom)*dij*pawang%realgnt(gint)*intg
                    case (3)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(1,iatom)*dij*pawang%realgnt(gint)*intg
                    end select
                 case (9) ! S_{2,2} contributes to xx, yy terms
                    select case (idir)
                    case (1)
                       bm2=bm2+scale_conversion*dtset%nucdipmom(1,iatom)*d22*pawang%realgnt(gint)*intg
                    case (2)
                       bm2=bm2-scale_conversion*dtset%nucdipmom(2,iatom)*d22*pawang%realgnt(gint)*intg
                    end select
                 end select
              end if ! end check on nonzero gaunt integral
           end do ! end loop over lp,mp
           cpb=cmplx(cprj_k(iatom,iband)%cp(1,ilmn),cprj_k(iatom,iband)%cp(2,ilmn),KIND=dpc)
           cpk=cmplx(cprj_k(iatom,iband)%cp(1,jlmn),cprj_k(iatom,iband)%cp(2,jlmn),KIND=dpc)
           onsite_bm_k_n=onsite_bm_k_n+conjg(cpb)*(bm1-bm2)*cpk
        end do ! end loop over ilmn
     end do ! end loop over jlmn
     ABI_FREE(ff)
  end do ! end loop over atoms

end subroutine make_onsite_bm_k_n
!!***

!!****f* ABINIT/make_S1trace_k_n
!! NAME
!! make_S1trace_k_n
!!
!! FUNCTION
!! Compute single band contribution to Trace[\rho^0_k S_k^{(1)} ] 
!! in orbital magnetism context
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine make_S1trace_k_n(adir,cprj_k,dtset,ENK,iband,nband_occ,pawtab,S1trace)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,iband,nband_occ
  real(dp),intent(in) :: ENK
  complex(dpc),intent(out) :: S1trace
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_occ)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,epsabg,gdir,iatom,ilmn,itypat,jlmn,klmn
  complex(dpc) :: cpb,cpk

!----------------------------------------------------------------

  S1trace = czero

  do epsabg = 1, -1, -2

    if (epsabg .EQ. 1) then
       bdir = modulo(adir,3)+1
       gdir = modulo(adir+1,3)+1
    else
       bdir = modulo(adir+1,3)+1
       gdir = modulo(adir,3)+1
    end if

    do iatom=1,dtset%natom
      itypat=dtset%typat(iatom)
      do ilmn=1,pawtab(itypat)%lmn_size
        do jlmn=1,pawtab(itypat)%lmn_size
          klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
          cpb=cmplx(cprj_k(iatom,iband)%dcp(1,bdir,ilmn),cprj_k(iatom,iband)%dcp(2,bdir,ilmn),KIND=dpc)
          cpk=cmplx(cprj_k(iatom,iband)%dcp(1,gdir,jlmn),cprj_k(iatom,iband)%dcp(2,gdir,jlmn),KIND=dpc)
          S1trace=S1trace-half*j_dpc*epsabg*ENK*conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
        end do ! end loop over jlmn
      end do ! end loop over ilmn
    end do ! end loop over atoms
  end do ! end loop over epsabg

end subroutine make_S1trace_k_n

!!***

!!****f* ABINIT/make_rhorij1_k_n
!! NAME
!! make_rhorij1_k_n
!!
!! FUNCTION
!! Compute Trace[\rho^0_k \rho_Rij(1)_k ] in orbital magnetism context
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine make_rhorij1_k_n(adir,cprj_k,dtset,iband,nband_occ,&
    & paw_ij,pawtab,rhorij1)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,iband,nband_occ
  complex(dpc),intent(out) :: rhorij1
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_occ)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,epsabg,gdir,iatom,ilmn,itypat,jlmn,klmn
  complex(dpc) :: cpb,cdij,cpk

!----------------------------------------------------------------

  rhorij1 = czero

  do epsabg = 1, -1, -2

    if (epsabg .EQ. 1) then
       bdir = modulo(adir,3)+1
       gdir = modulo(adir+1,3)+1
    else
       bdir = modulo(adir+1,3)+1
       gdir = modulo(adir,3)+1
    end if

    do iatom=1,dtset%natom
      itypat=dtset%typat(iatom)
      do ilmn=1,pawtab(itypat)%lmn_size
        do jlmn=1,pawtab(itypat)%lmn_size
          klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
          cpb=cmplx(cprj_k(iatom,iband)%dcp(1,bdir,ilmn),cprj_k(iatom,iband)%dcp(2,bdir,ilmn),KIND=dpc)
          cpk=cmplx(cprj_k(iatom,iband)%dcp(1,gdir,jlmn),cprj_k(iatom,iband)%dcp(2,gdir,jlmn),KIND=dpc)
          if (paw_ij(iatom)%cplex_dij .EQ. 2) then
             cdij=cmplx(paw_ij(iatom)%dij(2*klmn-1,1),paw_ij(iatom)%dij(2*klmn,1),KIND=dpc)
             if (jlmn .GT. ilmn) cdij=conjg(cdij)
          else
             cdij=cmplx(paw_ij(iatom)%dij(klmn,1),zero,KIND=dpc)
          end if
          rhorij1=rhorij1-half*j_dpc*epsabg*conjg(cpb)*cdij*cpk
        end do ! end loop over jlmn
      end do ! end loop over ilmn
    end do ! end loop over atoms
  end do ! end loop over epsabg

end subroutine make_rhorij1_k_n
!!***

!!****f* ABINIT/orbmag_ddk
!! NAME
!! orbmag_ddk
!!
!! FUNCTION
!! This routine computes the orbital magnetization and Berry curvature based on input 
!! wavefunctions and DDK wavefuntions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006) [[cite:Ceresoli2006]],
!! and Gonze and Zwanziger, PRB 84, 064445 (2011) [[cite:Gonze2011a]].
!! DDK wavefunctions are used for the derivatives.
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!
!! SOURCE

subroutine orbmag_ddk(atindx1,cg,cg1,dtset,gsqcut,kg,mcg,mcg1,mpi_enreg,&
    & nattyp,nfftf,ngfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,rprimd,vtrial,&
    & xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcg1,nfftf
 real(dp),intent(in) :: gsqcut
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
 integer,intent(in) :: nattyp(dtset%natom),ngfftf(18),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),cg1(2,mcg1,3),rprimd(3,3),xred(3,dtset%natom)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

 !Local
 !scalars
 integer :: adir,bdir,buff_size,dimffnl,exchn2n3d
 integer :: getcprj_choice,getcprj_cpopt
 integer :: getghc_cpopt,getghc_prtvol,getghc_sij_opt,getghc_tim,getghc_type_calc
 integer :: gdir,iatom,icg,ider,idir,ierr,ikg,ikg1,ikpt,ilm,isppol,istwf_k
 integer :: me,mcgk,my_nspinor
 integer :: nband_k,ncpgr,ndat,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,nn,nnp
 integer :: nkpg,npw_k
 integer :: nonlop_choice,nonlop_cpopt,nonlop_nnlout,nonlop_pawopt,nonlop_signs,nonlop_tim
 integer :: nproc,nterms,projbd_scprod_io,projbd_tim,projbd_useoverlap,spaceComm
 integer :: with_vectornd
 integer,parameter :: cci=1,vvii=2,vvia=3,vvib=4,rho0h1=5,rho0s1=6,lrr3=7,a0an=8,berrycurve=9
 real(dp) :: arg,dbi,dbr,dgi,dgr,doti,dub_dsg_i,dug_dsb_i
 real(dp) :: ecut_eff,Enk,lambda,local_fermie,trnrm,ucvol
 complex(dpc) :: dbc,dgc,onsite_bm_k_n,onsite_l_k_n,rhorij1,S1trace
 logical :: has_nucdip
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:),nattyp_dum(:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),lambda_ndat(1),nonlop_enlout(1),rhodum(1),rmet(3,3)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 real(dp),allocatable :: cg_k(:,:),cg1_k(:,:,:),cgrvtrial(:,:),cwaveb1(:,:),cwavef(:,:),cwavefp(:,:),cwaveg1(:,:)
 real(dp),allocatable :: cwavedsdb(:,:),cwavedsdg(:,:)
 real(dp),allocatable :: dscg_k(:,:,:),ffnl_k(:,:,:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:),orbmag_terms(:,:,:),orbmag_trace(:,:)
 real(dp),allocatable :: pcg1_k(:,:,:),ph1d(:,:),ph3d(:,:,:),phkxred(:,:),scg_k(:,:),scg1_k(:,:,:),scprod(:,:)
 real(dp),allocatable :: vectornd(:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cwaveprj(:,:)

 !----------------------------------------------

 ! set up basic FFT parameters
 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 nband_k = dtset%mband
 istwf_k = 1
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 me = mpi_enreg%me_kpt
 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0; ikg1 = 0

 ! input parameters for calls to nonlop
 ! nonlop_choice will be changed from call to call
 nonlop_cpopt = -1  ! cprj computed and not saved
 nonlop_pawopt = 3 ! apply only S
 nonlop_signs = 2  ! get <G|Op|C> vector
 nonlop_nnlout = 1
 nonlop_tim = 0

 ! input parameters to projbd
 projbd_scprod_io = 0
 projbd_useoverlap = 1
 projbd_tim = 0 

 ! input parameters for calls to getghc at ikpt
 getghc_cpopt = -1 ! cprj computed and not saved
 getghc_sij_opt = 1 ! compute both H|C> and S|C>
 ndat = 1           ! number of fft's in parallel
 getghc_prtvol = 0
 getghc_type_calc = 0 ! type_calc 0 means kinetic, local, nonlocal
 getghc_tim = 0
 lambda = zero 
 lambda_ndat = zero 

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ncpgr = 3
 ABI_MALLOC(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
 ABI_MALLOC(cprj_k,(dtset%natom,dtset%mband))
 call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
 ABI_MALLOC(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
      & paw_ij=paw_ij)

 !========= construct local potential ==================
 ! nspden=1 is essentially hard-coded in the following line
 ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
 call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
 ABI_MALLOC(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
 call fftpac(isppol,mpi_enreg,dtset%nspden,&
      & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)
 ABI_FREE(cgrvtrial)
 call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)
 
 !========  compute nuclear dipole vector potential (may be zero) ==========
 with_vectornd=0
 has_nucdip = ANY( ABS(dtset%nucdipmom) .GT. tol8 )
 if (has_nucdip) with_vectornd=1
 ABI_MALLOC(vectornd,(with_vectornd*nfftf,3))
 vectornd = zero
 if(has_nucdip) then
   call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,&
     & dtset%nucdipmom,rprimd,vectornd,xred)
   ABI_MALLOC(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
   ABI_MALLOC(cgrvtrial,(dtset%nfft,dtset%nspden))
   do idir = 1, 3
     call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
     call fftpac(isppol,mpi_enreg,dtset%nspden,&
       & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
   end do
   ABI_FREE(cgrvtrial)
   call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
 end if

 ABI_MALLOC(kg_k,(3,dtset%mpw))
 ABI_MALLOC(kinpw,(dtset%mpw))

 ABI_MALLOC(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
 call getph(atindx1,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

 icg = 0
 ikg = 0
 nterms = 9 ! various contributing terms in orbmag and berrycurve
 ! 1 orbmag CC
 ! 2 orbmag VV II
 ! 3 orbmag VV I+III part a
 ! 4 orbmag VV I+III part b 
 ! 5 orbmag Tr[\rho^0 H^1] with D^0_ij part
 ! 6 orbmag -Tr[\rho^0 S^1] part
 ! 7 orbmag onsite L_R/r^3
 ! 8 orbmag onsite A0.An
 ! 9 berrycurve
 ABI_MALLOC(orbmag_terms,(3,nterms,nband_k))
 orbmag_terms = zero
 local_fermie = -1.0D10
 
 !============= BIG FAT KPT LOOP :) ===========================
 do ikpt = 1, dtset%nkpt

   ! if the current kpt is not on the current processor, cycle
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

   ! trace norm: assume occupation of two for each band and weight by kpts
   trnrm = two*dtset%wtk(ikpt)

   kpoint(:)=dtset%kptns(:,ikpt)
   npw_k = npwarr(ikpt)

   ! retrieve kg_k at this k point
   kg_k(1:3,1:npw_k) = kg(1:3,ikg+1:ikg+npw_k)

   ! retrieve ylm at this k point
   ABI_MALLOC(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
   ABI_MALLOC(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
   do ilm=1,psps%mpsang*psps%mpsang
     ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
   end do

   ! Compute kinetic energy at kpt
   kinpw(:) = zero
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

   ! Compute k+G at this k point
   nkpg = 3
   ABI_MALLOC(kpg_k,(npw_k,nkpg))
   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

   ! Make 3d phase factors
   ABI_MALLOC(phkxred,(2,dtset%natom))
   do iatom=1, dtset%natom
     arg=two_pi*DOT_PRODUCT(kpoint,xred(:,iatom))
     phkxred(1,iatom)=cos(arg);phkxred(2,iatom)=sin(arg)
   end do
   ABI_MALLOC(ph3d,(2,npw_k,dtset%natom))
   call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,&
     & npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)

   ! Compute nonlocal form factors ffnl at all (k+G):
   ider=1 ! ffnl and 1st derivatives
   idir=4 ! ignored when ider = 0; idir=0 means d ffnl/ dk in reduced units referenced 
          ! to reciprocal translations
          ! idir=4 meand d ffnl / dk in reduced units referenced to real space
          ! translations. rfddk = 1 wavefunctions are computed using this convention.
   dimffnl=4 ! 1 + number of derivatives
   ABI_MALLOC(ffnl_k,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,&
     & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
     & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
     & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
     & psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

   !  - Load k-dependent quantities in the Hamiltonian
   call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
     & kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,&
     & compute_gbound=.TRUE.)

   ! retrieve ground state wavefunctions at this k point
   mcgk = npw_k*nband_k
   ABI_MALLOC(cg_k,(2,mcgk))
   cg_k = cg(1:2,icg+1:icg+mcgk)
   ABI_MALLOC(scg_k,(2,mcgk))
   ABI_MALLOC(scg1_k,(2,mcgk,3))

   ! retrieve first order wavefunctions at this k point
   ABI_MALLOC(cg1_k,(2,mcgk,3))
   cg1_k = cg1(1:2,icg+1:icg+mcgk,1:3)

   ! compute cprj_k, S|u_nk>, and S|du/dk>
   ABI_MALLOC(cwavef,(2,npw_k))
   ABI_MALLOC(gsc,(2,npw_k))
   ABI_MALLOC(gvnlc,(2,npw_k))
   ! input parameters for calls to nonlop
   nonlop_choice =  1! apply (I+S)
   ! input parameters for calls to getcprj
   getcprj_choice = 5 ! cprj and d cprj/dk
   getcprj_cpopt = 0 ! compute both cprj and d cprj/dk
   do nn = 1, nband_k
     cwavef = cg_k(:,(nn-1)*npw_k+1:nn*npw_k)
     call nonlop(nonlop_choice,nonlop_cpopt,cwaveprj,nonlop_enlout,gs_hamk,0,&
       & lambda_ndat,mpi_enreg,ndat,nonlop_nnlout,nonlop_pawopt,nonlop_signs,gsc,&
       & nonlop_tim,cwavef,gvnlc)
     scg_k(1:2,(nn-1)*npw_k+1:nn*npw_k) = gsc(1:2,1:npw_k)

     do adir = 1, 3
       call getcprj(getcprj_choice,getcprj_cpopt,cwavef,cwaveprj,ffnl_k,&
         & adir,psps%indlmn,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,&
         & dtset%mgfft,mpi_enreg,dtset%natom,nattyp,dtset%ngfft,&
         & dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
         & phkxred,ph1d,ph3d,ucvol,psps%useylm)

       call pawcprj_put(atindx1,cwaveprj,cprj_k,dtset%natom,&
         & nn,0,ikpt,0,isppol,nband_k,dtset%mkmem,&
         & dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0,&
         & mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)

     end do

     do adir = 1, 3
       cwavef = cg1_k(:,(nn-1)*npw_k+1:nn*npw_k,adir)
       call nonlop(nonlop_choice,nonlop_cpopt,cwaveprj,nonlop_enlout,gs_hamk,0,&
         & lambda_ndat,mpi_enreg,ndat,nonlop_nnlout,nonlop_pawopt,nonlop_signs,gsc,&
         & nonlop_tim,cwavef,gvnlc)
       scg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir) = gsc(1:2,1:npw_k)
     end do
   end do ! end loop over nn

   ! compute \partial S/\partial k |u_nk>
   ABI_MALLOC(dscg_k,(2,mcgk,3))
   ! input parameters for calls to nonlop
   nonlop_choice =  5! apply dS/dk
   do adir = 1, 3
     do nn = 1, nband_k
       cwavef = cg_k(:,(nn-1)*npw_k+1:nn*npw_k)
       call nonlop(nonlop_choice,nonlop_cpopt,cwaveprj,nonlop_enlout,gs_hamk,adir,&
         & lambda_ndat,mpi_enreg,ndat,nonlop_nnlout,nonlop_pawopt,nonlop_signs,gsc,&
         & nonlop_tim,cwavef,gvnlc)
       dscg_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir) = gsc(1:2,1:npw_k)
     end do ! end loop over nn
   end do

   ! compute projection of cg1_k on conduction space
   ABI_MALLOC(pcg1_k,(2,nband_k*npw_k,3))
   ABI_MALLOC(scprod,(2,nband_k))
   do adir = 1, 3
     do nn = 1, nband_k

       cwavef = cg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir)
       call projbd(cg_k,cwavef,-1,0,0,istwf_k,mcgk,mcgk,nband_k,npw_k,dtset%nspinor,&
         & scg_k,scprod,projbd_scprod_io,projbd_tim,projbd_useoverlap,&
         & mpi_enreg%me_g0,mpi_enreg%comm_fft)
       pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,adir) = cwavef

     end do ! end loop over nn
   end do ! end loop over adir

   ABI_MALLOC(ghc,(2,npw_k))
   ABI_MALLOC(cwaveb1,(2,npw_k))
   ABI_MALLOC(cwaveg1,(2,npw_k))
   ABI_MALLOC(cwavedsdb,(2,npw_k))
   ABI_MALLOC(cwavedsdg,(2,npw_k))
   ABI_MALLOC(cwavefp,(2,npw_k))
   do nn = 1, nband_k

     ! compute H^0|u_nk> and <u_nk|H^0|u_nk>
     cwavef(1:2,1:npw_k) = cg_k(1:2,(nn-1)*npw_k+1:nn*npw_k)
     call getghc(getghc_cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
       & getghc_prtvol,getghc_sij_opt,getghc_tim,getghc_type_calc)
     Enk = DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) &
       & + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))

     if (Enk .GT. local_fermie) local_fermie = Enk

     do adir =1, 3
       bdir = modulo(adir,3)+1
       gdir = modulo(adir+1,3)+1

       ! 1 orbmag CC
       ! -i/2 eps_abg <du/dg|P_c H0 P_c|du/db> =
       ! -i/2 (<du/dg|P_c H0 P_c|du/db> - <du/db|P_c H0 P_c|du/dg>) =
       ! -i/2 (2 i Im<du/dg|P_c H0 P_c|du/db>) =
       ! Im<du/dg|P_c H0 P_c|du/db>

       cwaveb1(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
       call getghc(getghc_cpopt,cwaveb1,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
         & getghc_prtvol,getghc_sij_opt,getghc_tim,getghc_type_calc)

       cwaveg1(1:2,1:npw_k) = pcg1_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)
       doti=-DOT_PRODUCT(cwaveg1(2,:),ghc(1,:))+DOT_PRODUCT(cwaveg1(1,:),ghc(2,:))

       orbmag_terms(adir,cci,nn) = orbmag_terms(adir,cci,nn) + doti*trnrm
       
       ! 2 orbmag VV II
       ! vv needs (+i/2)*eps_abg*<du/db|P_c S P_c|du/dg>Enk =
       ! +i/2 (<du/db|P_c S P_c|du/dg> - <du/dg|P_c S P_c|du/db>)Enk =
       ! -i/2 (<du/dg|P_c S P_c|du/db> - <du/db|P_c S P_c|du/dg>)Enk =
       ! Im<du/dg|P_c S P_c|du/db>Enk
       doti=-DOT_PRODUCT(cwaveg1(2,:),gsc(1,:))+DOT_PRODUCT(cwaveg1(1,:),gsc(2,:))
       orbmag_terms(adir,vvii,nn) = orbmag_terms(adir,vvii,nn) + doti*Enk*trnrm


       !VV Ib term gives (i/2)eps_abg <du/db|P_c dS/dg|u>Enk
       !VV IIIb term gives (i/2)eps_abg <du|dS/db P_c|du/dg>Enk
       ! combined with eps_abg contraction they contribute
       ! -Im(VVI)*Enk
       ! 4 orbmag VV I+III part b
       cwavedsdb(1:2,1:npw_k) = dscg_k(1:2,(nn-1)*npw_k+1:nn*npw_k,bdir)
       cwavedsdg(1:2,1:npw_k) = dscg_k(1:2,(nn-1)*npw_k+1:nn*npw_k,gdir)

       dug_dsb_i = -DOT_PRODUCT(cwaveg1(2,:),cwavedsdb(1,:)) + DOT_PRODUCT(cwaveg1(1,:),cwavedsdb(2,:))
       dub_dsg_i = -DOT_PRODUCT(cwaveb1(2,:),cwavedsdg(1,:)) + DOT_PRODUCT(cwaveb1(1,:),cwavedsdg(2,:))
       orbmag_terms(adir,vvib,nn)= orbmag_terms(adir,vvib,nn) - (dub_dsg_i-dug_dsb_i)*Enk*trnrm

       ! VVIa term gives (i/2)eps_abg sum_n' (-)<u_n|dS/db|u_n'><u_n'|dS/dg|u_n>Enk
       ! = + sum_n' Im{<u_n|dS/db|u_n'><u_n'|dS/dg|u_n>Enk}
       ! VVIIIa is identical. VVIIa is the negative of VVIa. The total contribution of all
       ! three terms is thus the same as VVIa itself.
       ! term 3 
       do nnp=1,nband_k
         cwavefp(1:2,1:npw_k) = cg_k(1:2,(nnp-1)*npw_k+1:nnp*npw_k)
         dbr= DOT_PRODUCT(cwavefp(1,:),cwavedsdb(1,:))+DOT_PRODUCT(cwavefp(2,:),cwavedsdb(2,:))
         dbi=-DOT_PRODUCT(cwavefp(2,:),cwavedsdb(1,:))+DOT_PRODUCT(cwavefp(1,:),cwavedsdb(2,:))
         dbc=cmplx(dbr,dbi,kind=dpc)
         dgr= DOT_PRODUCT(cwavefp(1,:),cwavedsdg(1,:))+DOT_PRODUCT(cwavefp(2,:),cwavedsdg(2,:))
         dgi=-DOT_PRODUCT(cwavefp(2,:),cwavedsdg(1,:))+DOT_PRODUCT(cwavefp(1,:),cwavedsdg(2,:))
         dgc=cmplx(dgr,dgi,kind=dpc)
         orbmag_terms(adir,vvia,nn) = orbmag_terms(adir,vvia,nn) + AIMAG(CONJG(dbc)*dgc*Enk)*trnrm
       end do
      
       ! 5 Tr[-\rho^0 S^1 \rho^0 H^0] contribution 
       call make_S1trace_k_n(adir,cprj_k,dtset,Enk,nn,nband_k,pawtab,S1trace)
       orbmag_terms(adir,rho0s1,nn) = orbmag_terms(adir,rho0s1,nn) - real(S1trace)*trnrm
       
       ! 6 Tr[\rho^0 H^1] contribution:
       ! -i/2 eps_abg <u|dp/db>D_ij^0<dp/dg|u>
       call make_rhorij1_k_n(adir,cprj_k,dtset,nn,nband_k,paw_ij,pawtab,rhorij1)
       orbmag_terms(adir,rho0h1,nn) = orbmag_terms(adir,rho0h1,nn) + real(rhorij1)*trnrm

       ! 7 onsite L_R/r^3 contribution
       call make_onsite_l_k_n(cprj_k,dtset,nn,adir,nband_k,onsite_l_k_n,pawrad,pawtab)
       orbmag_terms(adir,lrr3,nn) = orbmag_terms(adir,lrr3,nn) + real(onsite_l_k_n)*trnrm

       ! 8 onsite A0.An contiribution
       call make_onsite_bm_k_n(cprj_k,dtset,nn,adir,nband_k,onsite_bm_k_n,&
         & pawang,pawrad,pawtab)
       orbmag_terms(adir,a0an,nn) = orbmag_terms(adir,a0an,nn) + real(onsite_bm_k_n)*trnrm

       ! berrycurve needs i*eps_abg*<du/db|S|du/dg> 
       ! N.B. the Berry curvature does not involve H0 so no projection onto conduction
       ! and valence bands, the "S" here is really I+S from PAW
       ! i eps_abg <du/db|S|du/dg> = -2*Im<du/db|S|du/dg> 
       ! 9 berrycurve
       doti = -DOT_PRODUCT(cg1_k(2,(nn-1)*npw_k+1:nn*npw_k,bdir),scg1_k(1,(nn-1)*npw_k+1:nn*npw_k,gdir)) + &
             & DOT_PRODUCT(cg1_k(1,(nn-1)*npw_k+1:nn*npw_k,bdir),scg1_k(2,(nn-1)*npw_k+1:nn*npw_k,gdir))
       orbmag_terms(adir,berrycurve,nn) = orbmag_terms(adir,berrycurve,nn) - two*doti*trnrm

     end do

   end do

   ABI_FREE(cwavef)
   ABI_FREE(cwavefp)
   ABI_FREE(cwaveb1)
   ABI_FREE(cwaveg1)
   ABI_FREE(cwavedsdb)
   ABI_FREE(cwavedsdg)
   ABI_FREE(ghc)
   ABI_FREE(gsc)
   ABI_FREE(gvnlc)
   ABI_FREE(cg_k)
   ABI_FREE(scg_k)
   ABI_FREE(scg1_k)
   ABI_FREE(dscg_k)
   ABI_FREE(cg1_k)
   ABI_FREE(pcg1_k)
   ABI_FREE(scprod)

   icg = icg + npw_k*nband_k
   ikg = ikg + npw_k

   ABI_FREE(ylm_k)
   ABI_FREE(ylmgr_k)
   ABI_FREE(kpg_k)
   ABI_FREE(ffnl_k)
   ABI_FREE(ph3d)
   ABI_FREE(phkxred)

 end do ! end loop over kpts

 if (nproc > 1) then
   buff_size=size(orbmag_terms)
   ABI_MALLOC(buffer1,(buff_size))
   ABI_MALLOC(buffer2,(buff_size))
   buffer1(1:buff_size) = reshape(orbmag_terms,(/3*nterms*nband_k/))
   call xmpi_sum(buffer1,buffer2,buff_size,spaceComm,ierr)
   orbmag_terms(1:3,1:nterms,1:nband_k)=reshape(buffer2,(/3,nterms,nband_k/))
   ABI_FREE(buffer1)
   ABI_FREE(buffer2)
 end if

 ! convert to cartesian frame, supply 1/(2\pi)^2 factor
 ! but not to lrr3 and a0an terms, they are already cartesian and don't require 
 ! additional normalization
 do nn = 1, nband_k
   orbmag_terms(1:3,cci,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,cci,nn))
   orbmag_terms(1:3,vvii,nn) = (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,vvii,nn))
   orbmag_terms(1:3,vvib,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,vvib,nn))
   orbmag_terms(1:3,vvia,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,vvia,nn))
   orbmag_terms(1:3,rho0h1,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,rho0h1,nn))
   orbmag_terms(1:3,rho0s1,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,rho0s1,nn))
   orbmag_terms(1:3,berrycurve,nn) =  (ucvol/(two_pi*two_pi))*MATMUL(gprimd,orbmag_terms(1:3,berrycurve,nn))
 end do

 ! compute trace of each term
 ABI_MALLOC(orbmag_trace,(3,nterms))
 orbmag_trace = zero
 do nn = 1, nband_k
   orbmag_trace(1:3,1:nterms) = orbmag_trace(1:3,1:nterms) + orbmag_terms(1:3,1:nterms,nn)
 end do

 call orbmag_ddk_output(dtset,local_fermie,nband_k,nterms,orbmag_terms,orbmag_trace)

!---------------------------------------------------
! deallocate memory
!---------------------------------------------------
 call gs_hamk%free()

 ABI_FREE(vlocal)
 ABI_FREE(vectornd)
 if(has_nucdip) then
   ABI_FREE(vectornd_pac)
 end if
 ABI_FREE(kg_k)
 ABI_FREE(kinpw)
 ABI_FREE(ph1d)
 ABI_FREE(orbmag_terms)
 ABI_FREE(orbmag_trace)

 ABI_FREE(dimlmn)
 call pawcprj_free(cprj_k)
 ABI_FREE(cprj_k)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)

end subroutine orbmag_ddk
!!***

!!****f* ABINIT/orbmag_ddk_output
!! NAME
!! orbmag_ddk_output
!!
!! FUNCTION
!! This routine outputs orbmag terms tailored for ddk routine
!!
!! COPYRIGHT
!! Copyright (C) 2003-2021 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_mpi_recv
!!      pawcprj_mpi_send,xmpi_sum
!!
!! SOURCE

subroutine orbmag_ddk_output(dtset,fermie,nband_k,nterms,orbmag_terms,orbmag_trace)


 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: nband_k,nterms
 real(dp),intent(in) :: fermie
 type(dataset_type),intent(in) :: dtset

 !arrays
 real(dp),intent(in) :: orbmag_terms(3,nterms,nband_k),orbmag_trace(3,nterms)

 !Local variables -------------------------
 !scalars
 integer :: adir,iband,iterms
 integer,parameter :: cci=1,vvii=2,vvia=3,vvib=4,rho0h1=5,rho0s1=6,lrr3=7,a0an=8,berrycurve=9
 character(len=500) :: message

 !arrays
 real(dp) :: berrycurve_total(3),orbmag_total(3)

 ! ***********************************************************************

 orbmag_total=zero;berrycurve_total=zero
 do iterms = 1, nterms-1
   orbmag_total(1:3)=orbmag_total(1:3) + orbmag_trace(1:3,iterms)
 end do
 berrycurve_total(1:3)=orbmag_trace(1:3,berrycurve)

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 if(dtset%orbmag .GE. 1) then
   write(message,'(a)')' Orbital magnetic moment, Cartesian directions : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(3es16.8)') (orbmag_total(adir),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Integral of Berry curvature, Cartesian directions : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(3es16.8)') (berrycurve_total(adir),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,es16.8)')' Fermie energy : ',fermie
   call wrtout(ab_out,message,'COLL')
 end if

 if(dtset%orbmag .GE. 2) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, Term-by-term breakdown : '
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '           Conduction space : ',(orbmag_trace(adir,cci),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '          Valence space IIb : ',(orbmag_trace(adir,vvii),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '  Valence space Ia+IIa+IIIa : ',(orbmag_trace(adir,vvia),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '          Valence space Ib  : ',(orbmag_trace(adir,vvib),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '           S(1) PAW overlap : ',(orbmag_trace(adir,rho0s1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '                  H(1) cprj : ',(orbmag_trace(adir,rho0h1),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '           H(1) on-site L_R : ',(orbmag_trace(adir,lrr3),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '         H(1) on-site A0.An : ',(orbmag_trace(adir,a0an),adir=1,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3es16.8)') '            Berry curvature : ',(orbmag_trace(adir,berrycurve),adir=1,3)
   call wrtout(ab_out,message,'COLL')
 end if

 if(dtset%orbmag .EQ. 3) then
   write(message,'(a)')ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)')' Orbital magnetic moment, Term-by-term breakdown for each band : '
   call wrtout(ab_out,message,'COLL')
   do iband = 1, nband_k
     write(message,'(a)')ch10
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,i2,a,i2)') ' band ',iband,' of ',nband_k
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '            Conduction space : ',(orbmag_terms(adir,cci,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '           Valence space IIb : ',(orbmag_terms(adir,vvii,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '  Valence space Ia+IIa+IIIa  : ',(orbmag_terms(adir,vvia,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '           Valence space Ib  : ',(orbmag_terms(adir,vvib,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '            S(1) PAW overlap : ',(orbmag_terms(adir,rho0s1,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '                   H(1) cprj : ',(orbmag_terms(adir,rho0h1,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '            H(1) on-site L_R : ',(orbmag_terms(adir,lrr3,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '          H(1) on-site A0.An : ',(orbmag_terms(adir,a0an,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,3es16.8)') '             Berry curvature : ',(orbmag_terms(adir,berrycurve,iband),adir=1,3)
     call wrtout(ab_out,message,'COLL')
   end do
 end if

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine orbmag_ddk_output
!!***

end module m_orbmag
