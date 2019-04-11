!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (JWZ)
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

  use defs_abitypes
  use defs_basis
  use defs_datatypes
  use m_errors
  use m_abicore
  use m_xmpi

  use m_berrytk,          only : smatrix
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : overlap_g
  use m_fft,              only : fftpac
  use m_fftcore,          only : kpgsph
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian,destroy_hamiltonian,&
       &                         load_spin_hamiltonian,load_k_hamiltonian,gs_hamiltonian_type
  use m_initylmg,         only : initylmg
  use m_kg,               only : getph,mkkin,mkkpg,mkpwind_k,mknucdipmom_k,ph1d3d
  use m_kpts,             only : listkk, smpbz
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle
  use m_nonlop,           only : nonlop
  use m_pawang,           only : pawang_type
  use m_pawfgr,           only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_paw_overlap,      only : overlap_k1k2_paw
  use m_pawrad,           only : pawrad_type,pawrad_deducer0,simp_gen
  use m_paw_sphharm,      only : initylmr,setsym_ylm,slxyzs
  use m_pawtab,           only : pawtab_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_free,&
       &                         pawcprj_get, pawcprj_put, pawcprj_getdim, pawcprj_set_zero,&
       &                         pawcprj_symkn,pawcprj_mpi_recv,pawcprj_mpi_send
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
  public :: rho_norm_check
  public :: chern_number
  public :: mpi_chern_number
  public :: make_dsdk
  public :: make_dsdk_FD
  public :: make_onsite_l_k
  public :: make_S1trace_k
  public :: make_smat
  public :: make_CCIV_k
  public :: make_CCIV_k_FD
  public :: orbmag
  public :: ctocprjb

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

    implicit none

    !Arguments ------------------------------------
    !array
    type(orbmag_type),intent(inout) :: dtorbmag

    ! ************************************************************************

    ! Integer pointers
    if(allocated(dtorbmag%atom_indsym))  then
       ABI_DEALLOCATE(dtorbmag%atom_indsym)
    end if
    if(allocated(dtorbmag%cgindex))  then
       ABI_DEALLOCATE(dtorbmag%cgindex)
    end if
    if(allocated(dtorbmag%cprjindex))  then
       ABI_DEALLOCATE(dtorbmag%cprjindex)
    end if
    if(allocated(dtorbmag%fkgindex))  then
       ABI_DEALLOCATE(dtorbmag%fkgindex)
    end if
    if(allocated(dtorbmag%ikpt_dk))  then
       ABI_DEALLOCATE(dtorbmag%ikpt_dk)
    end if
    if(allocated(dtorbmag%indkk_f2ibz))  then
       ABI_DEALLOCATE(dtorbmag%indkk_f2ibz)
    end if
    if(allocated(dtorbmag%i2fbz))  then
       ABI_DEALLOCATE(dtorbmag%i2fbz)
    end if
    if(allocated(dtorbmag%kg)) then
       ABI_DEALLOCATE(dtorbmag%kg)
    end if
    if(allocated(dtorbmag%kgindex))  then
       ABI_DEALLOCATE(dtorbmag%kgindex)
    end if
    if(allocated(dtorbmag%lmn_size))  then
       ABI_DEALLOCATE(dtorbmag%lmn_size)
    end if
    if(allocated(dtorbmag%lmn2_size))  then
       ABI_DEALLOCATE(dtorbmag%lmn2_size)
    end if
    if(allocated(dtorbmag%nband_occ))  then
       ABI_DEALLOCATE(dtorbmag%nband_occ)
    end if
    ! Real(dp) pointers

    if(allocated(dtorbmag%fkptns))  then
       ABI_DEALLOCATE(dtorbmag%fkptns)
    end if
    if(allocated(dtorbmag%zarot))  then
       ABI_DEALLOCATE(dtorbmag%zarot)
    end if

  end subroutine destroy_orbmag
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/initorbmag
!! NAME
!! initorbmag
!!
!! FUNCTION
!! Initialization of orbital magnetization calculation; similar to initberry
!!
!! COPYRIGHT
!! Copyright (C) 2004-2019 ABINIT group.
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
!!      m_gstate
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine initorbmag(dtorbmag,dtset,gmet,gprimd,kg,mpi_enreg,npwarr,occ,&
     &                     pawtab,psps,pwind,pwind_alloc,pwnsfac,&
     &                     rprimd,symrec,xred)

  implicit none

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
     ABI_ALLOCATE(spkpt,(3,mkpt))
     call smpbz(1,ab_out,dtset%kptrlatt,mkpt,fnkpt_computed,dtset%nshiftk,option,dtset%shiftk,spkpt)
     dtorbmag%fnkpt = fnkpt_computed
     ABI_ALLOCATE(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
     dtorbmag%fkptns(:,:)=spkpt(:,1:dtorbmag%fnkpt)
     ABI_DEALLOCATE(spkpt)
  else if(dtset%kptopt==3.or.dtset%kptopt==0)then
     dtorbmag%fnkpt=dtset%nkpt
     ABI_ALLOCATE(dtorbmag%fkptns,(3,dtorbmag%fnkpt))
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
  ABI_ALLOCATE(dtorbmag%indkk_f2ibz,(dtorbmag%fnkpt,6))

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
  ABI_ALLOCATE(dtorbmag%i2fbz,(dtset%nkpt))
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
     MSG_ERROR(message)
  end if

  !----------------------------------------------------------------------------
  !------------- Allocate PAW space as necessary ------------------------------
  !----------------------------------------------------------------------------

  dtorbmag%usepaw   = psps%usepaw
  dtorbmag%natom    = dtset%natom
  dtorbmag%my_natom = mpi_enreg%my_natom

  ABI_ALLOCATE(dtorbmag%lmn_size,(dtset%ntypat))
  ABI_ALLOCATE(dtorbmag%lmn2_size,(dtset%ntypat))
  do itypat = 1, dtset%ntypat
     dtorbmag%lmn_size(itypat) = pawtab(itypat)%lmn_size
     dtorbmag%lmn2_size(itypat) = pawtab(itypat)%lmn2_size
  end do

  lmn2_size_max = psps%lmnmax*(psps%lmnmax+1)/2
  dtorbmag%lmn2max = lmn2_size_max

  ABI_ALLOCATE(dtorbmag%cprjindex,(dtset%nkpt,dtset%nsppol))
  dtorbmag%cprjindex(:,:) = 0

  if (dtset%kptopt /= 3) then
     ABI_ALLOCATE(dtorbmag%atom_indsym,(4,dtset%nsym,dtorbmag%natom))
     call symatm(dtorbmag%atom_indsym,dtorbmag%natom,dtset%nsym,symrec,dtset%tnons,tol8,dtset%typat,xred)
     lmax = psps%mpsang - 1
     ABI_ALLOCATE(dtorbmag%zarot,(2*lmax+1,2*lmax+1,lmax+1,dtset%nsym))
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

  ABI_ALLOCATE(mpi_enreg%kpt_loc2fbz_sp,(0:nproc-1,1:dtorbmag%fmkmem_max*dtset%nsppol, 1:2))
  ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:nproc-1,1:dtorbmag%mkmem_max*dtset%nsppol, 1:2))
  ABI_ALLOCATE(mpi_enreg%kptdstrb,(nproc,6,dtorbmag%fmkmem_max*dtset%nsppol*2))
  ABI_ALLOCATE(mpi_enreg%mkmem,(0:nproc-1))
  mpi_enreg%kpt_loc2fbz_sp(:,:,:) = 0
  mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
  mpi_enreg%kptdstrb(:,:,:)       = 0
  mpi_enreg%mkmem(:)              = 0

  pwind_alloc = dtset%mpw*dtorbmag%fmkmem_max

  ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
  ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))

  ! !------------------------------------------------------------------------------
  ! !---------------------- Compute orbmag_type variables -------------------------
  ! !------------------------------------------------------------------------------

  !Initialization of orbmag_type variables
  dtorbmag%dkvecs(:,:) = zero
  ABI_ALLOCATE(dtorbmag%ikpt_dk,(dtorbmag%fnkpt,2,3))
  ABI_ALLOCATE(dtorbmag%cgindex,(dtset%nkpt,dtset%nsppol))
  ABI_ALLOCATE(dtorbmag%kgindex,(dtset%nkpt))
  ABI_ALLOCATE(dtorbmag%fkgindex,(dtorbmag%fnkpt))
  dtorbmag%ikpt_dk(:,:,:) = 0
  dtorbmag%cgindex(:,:) = 0
  dtorbmag%mband_occ = 0
  ABI_ALLOCATE(dtorbmag%nband_occ,(dtset%nsppol))
  dtorbmag%kgindex(:) = 0
  dtorbmag%fkgindex(:) = 0
  ABI_ALLOCATE(dtorbmag%kg,(3,dtset%mpw*dtset%mkmem))
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
              MSG_ERROR(message)
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
  ABI_ALLOCATE(kg1_k,(3,dtset%mpw))

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

  ABI_DEALLOCATE(kg1_k)


  call timab(1009,2,tsec)
  call timab(1001,2,tsec)

  DBG_EXIT("COLL")

end subroutine initorbmag
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/rho_norm_check
!! NAME
!! rho_norm_check
!!
!! FUNCTION
!! Routine to play with and check density operator normalization.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!! only printing
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE


subroutine rho_norm_check(atindx1,cg,cprj,dtorbmag,dtset,mpi_enreg,mcg,mcprj,&
     & npwarr,pawtab,usecprj,usepaw)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,usecprj,usepaw
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(in) :: dtorbmag

  !arrays
  integer,intent(in) :: atindx1(dtset%natom)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: cg(2,mcg)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

  !Local variables -------------------------
  !scalars
  integer :: iatom,iband,icg,icprj,ikpt,ilmn,ipw,isppol,itypat,jlmn,klmn
  integer :: my_nspinor,nband_k,ncpgr,npw_k
  real(dp) :: cr,ci,dotr,trace
  complex(dpc) :: cpb,cpk,consite
  !arrays
  integer,allocatable :: dimlmn(:),nattyp_dum(:)
  type(pawcprj_type),allocatable :: cprj_k(:,:)

  !----------------------------------------------------
  
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  
  nband_k = dtorbmag%mband_occ

  trace=zero
  do ikpt = 1, dtorbmag%fnkpt
           
     icprj = dtorbmag%cprjindex(ikpt,isppol)

     npw_k = npwarr(ikpt)
     icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

     do iband = 1, nband_k
        dotr=zero
        do ipw=1, npw_k
           cr=cg(1,icg+npw_k*(iband-1)+ipw);ci=cg(2,icg+npw_k*(iband-1)+ipw)
           dotr=dotr+cr*cr+ci*ci
        end do
        do iatom=1,dtset%natom
           itypat=dtset%typat(iatom)
           do klmn=1,pawtab(itypat)%lmn2_size
              ilmn=pawtab(itypat)%indklmn(7,klmn)
              jlmn=pawtab(itypat)%indklmn(8,klmn)
              cpb=cmplx(cprj_k(iatom,iband)%cp(1,ilmn),cprj_k(iatom,iband)%cp(2,ilmn))
              cpk=cmplx(cprj_k(iatom,iband)%cp(1,jlmn),cprj_k(iatom,iband)%cp(2,jlmn))
              consite=conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk*pawtab(itypat)%dltij(klmn)
              dotr=dotr+real(consite)
           end do
        end do
        trace=trace+dotr
     end do

  end do ! end loop over ikpt

  ! final trace: factor of two assumes two electrons per band (normal occupance for an insulator)
  trace = trace*two/dtorbmag%fnkpt
  
  write(std_out,'(a,2i4,es16.8)')'JWZ debug nkpt nband_k trace ',dtorbmag%fnkpt,nband_k,trace
  
  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)

end subroutine rho_norm_check
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/chern_number
!! NAME
!! chern_number
!!
!! FUNCTION
!! This routine computes the Chern number based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtset <type(dataset_type)>=all input variables in this dataset
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! rprimd(3,3) = real space translation vectors
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw=1 if PAW calculation
!! xred(3,natom) = location of atoms in unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006) [[cite:Ceresoli2006]],
!! and Gonze and Zwanziger, PRB 84 064445 (2011) [[cite:Gonze2011a]].
!! This routine computes the Chern number as
!! $C_\alpha = \frac{i}{2\pi}\int_{\mathrm{BZ}} dk \epsilon_{\alpha\beta\gamma}
!! \mathrm{Tr}[\rho_k \partial_\beta \rho_k (1 - \rho_k) \partial_gamma\rho_k] $
!! The derivative of the density operator is obtained from a discretized formula
!! $\partial_\beta \rho_k = \frac{1}{2\Delta}(\rho_{k+b} - \rho_{k-b})$ with
!! $\Delta = |b|$. When reduced to wavefunction overlaps the computation amounts to
!! multiple calls to smatrix.F90, exactly as in other Berry phase computations, with
!! the one additional complication of overlaps like $\langle u_{n,k+b}|u_{n',k+g}\rangle$.
!! At this stage mkpwind_k is invoked, which generalizes the code in initberry
!! and initorbmag necessary to index plane waves around different k points.
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine chern_number(atindx1,cg,cprj,dtset,dtorbmag,kg,&
     &            mcg,mcprj,mpi_enreg,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,&
     &            rprimd,symrec,usecprj,usepaw,xred)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,pwind_alloc,usecprj,usepaw
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: cg(2,mcg),rprimd(3,3),xred(3,dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,bdx,bdxc,bfor,bsigma,ddkflag,epsabg,gdir,gdx,gdxc,gfor,gsigma
  integer :: icg,icgb,icgg,icprj,icprji,icprjb,icprjbi,icprjg,icprjgi
  integer :: ikg,ikpt,ikpti,ikptb,ikptbi,ikptg,ikptgi,isppol,itrs,job
  integer :: mcg1_k,my_nspinor,nband_k,ncpgr,nn,n1,n2,n3,npw_k,npw_kb,npw_kg,shiftbd
  real(dp) :: deltab,deltag,ucvol
  complex(dpc) :: IA,IB,t1A,t2A,t3A,t1B,t2B,t3B,t4B
  character(len=500) :: message
  !arrays
  integer,allocatable :: dimlmn(:),nattyp_dum(:),pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
  real(dp) :: cnum(2,3),dkb(3),dkg(3),dkbg(3),dtm_k(2),gmet(3,3),gprimd(3,3),rmet(3,3)
  real(dp),allocatable :: cg1_k(:,:),kk_paw(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_all_indx(:,:,:,:,:,:),smat_inv(:,:,:),smat_kk(:,:,:)
  logical,allocatable :: has_smat_indx(:,:,:)
  type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)
  type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

  ! ***********************************************************************
  ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

  ! ! check density operator norm
  ! call rho_norm_check(atindx1,cg,cprj,dtorbmag,dtset,mpi_enreg,mcg,mcprj,&
  !    & npwarr,pawtab,usecprj,usepaw)

  ! TODO: generalize to nsppol > 1
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  nband_k = dtorbmag%mband_occ

  if (usepaw == 1) then ! cprj allocation
     ncpgr = cprj(1,1)%ncpgr
     ABI_ALLOCATE(dimlmn,(dtset%natom))
     call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
     ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
     ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
     ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
     if (dtset%kptopt /= 3) then
        ABI_DATATYPE_ALLOCATE(cprj_ikn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
        ABI_DATATYPE_ALLOCATE(cprj_fkn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
        call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
        call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
     end if
  else
     message = ' usepaw /= 1 but Chern number calculation requires PAW '
     MSG_ERROR(message)
  end if

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))
  ABI_ALLOCATE(pwind_kg,(dtset%mpw))
  ABI_ALLOCATE(pwind_bg,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = 1.0_dp ! bra real
  pwnsfac_k(2,:) = 0.0_dp ! bra imag
  pwnsfac_k(3,:) = 1.0_dp ! ket real
  pwnsfac_k(4,:) = 0.0_dp ! ket imag

  mcg1_k = dtset%mpw*nband_k
  ABI_ALLOCATE(cg1_k,(2,mcg1_k))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

  ddkflag = 1

  !itrs = 0 means do not invoke time reversal symmetry in smatrix.F90
  itrs = 0

  job = 1
  shiftbd = 1

  ABI_ALLOCATE(has_smat_indx,(dtorbmag%fnkpt,0:6,0:6))
  ABI_ALLOCATE(smat_all_indx,(2,nband_k,nband_k,dtorbmag%fnkpt,0:6,0:6))
  has_smat_indx(:,:,:)=.FALSE.

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

  ! loop over kpts, assuming for now kptopt 3 or 4, nsppol = 1, nspinor = 1
  ! and no parallelism, no symmorphic symmetry elements

  cnum(:,:) = zero
  do ikpt = 1, dtorbmag%fnkpt
     
     ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
     
     icprji = dtorbmag%cprjindex(ikpti,isppol)
     
     npw_k = npwarr(ikpti)
     icg = dtorbmag%cgindex(ikpti,dtset%nsppol)
     
     ikg = dtorbmag%fkgindex(ikpt)
     
     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
     if ( ikpti /= ikpt ) then
        call pawcprj_copy(cprj_k,cprj_ikn)
        call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
             & dtorbmag%indkk_f2ibz(ikpt,2),dtorbmag%indkk_f2ibz(ikpt,6),&
             & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikpti)),&
             & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
             & dtorbmag%mband_occ,my_nspinor,&
             & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
        call pawcprj_copy(cprj_fkn,cprj_k)
     end if

     do adir = 1, 3

        do epsabg = 1, -1, -2

           if (epsabg .EQ. 1) then
              bdir = modulo(adir,3)+1
              gdir = modulo(adir+1,3)+1
           else
              bdir = modulo(adir+1,3)+1
              gdir = modulo(adir,3)+1
           end if

           do bfor = 1, 2
              if (bfor .EQ. 1) then
                 bsigma = 1
              else
                 bsigma = -1
              end if
              ! index of neighbor 1..6
              bdx = 2*bdir-2+bfor
              ! index of ikpt viewed from neighbor
              bdxc = 2*bdir-2+bfor+bsigma

              dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
              ! deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))
              deltab = sqrt(DOT_PRODUCT(dkb,dkb))

              ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
              ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)

              icprjbi = dtorbmag%cprjindex(ikptbi,isppol)

              npw_kb = npwarr(ikptbi)
              icgb = dtorbmag%cgindex(ikptbi,dtset%nsppol)

              pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

              call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjbi,&
                   &         ikptbi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                   &         my_nspinor,dtset%nsppol,0)
              if ( ikptbi /= ikptb ) then
                 call pawcprj_copy(cprj_kb,cprj_ikn)
                 call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                      & dtorbmag%indkk_f2ibz(ikptb,2),dtorbmag%indkk_f2ibz(ikptb,6),&
                      & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptbi)),&
                      & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                      & dtorbmag%mband_occ,my_nspinor,&
                      & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                 call pawcprj_copy(cprj_fkn,cprj_kb)
              end if
              
              if (.NOT. has_smat_indx(ikpt,bdx,0)) then

                 call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                      &           dtorbmag%lmn_size,dtset%mband,&
                      &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                 sflag_k=0
                 call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
                      &           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                      &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                 smat_all_indx(:,:,:,ikpt,bdx,0) = smat_kk(:,:,:)
                 smat_all_indx(1,:,:,ikptb,bdxc,0) = TRANSPOSE(smat_kk(1,:,:))
                 smat_all_indx(2,:,:,ikptb,bdxc,0) = -TRANSPOSE(smat_kk(2,:,:))

                 has_smat_indx(ikpt,bdx,0) = .TRUE.
                 has_smat_indx(ikptb,bdxc,0) = .TRUE.

              end if

              do gfor = 1, 2
                 if (gfor .EQ. 1) then
                    gsigma = 1
                 else
                    gsigma = -1
                 end if
                 ! index of neighbor 1..6
                 gdx = 2*gdir-2+gfor
                 ! index of ikpt viewed from neighbor
                 gdxc = 2*gdir-2+gfor+gsigma

                 dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
                 ! deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))
                 deltag = sqrt(DOT_PRODUCT(dkg,dkg))

                 ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                 ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)

                 icprjgi = dtorbmag%cprjindex(ikptgi,isppol)

                 npw_kg = npwarr(ikptgi)
                 icgg = dtorbmag%cgindex(ikptgi,dtset%nsppol)

                 pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                 call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjgi,&
                      &           ikptgi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                      &           my_nspinor,dtset%nsppol,0)
                 if ( ikptgi /= ikptg ) then
                    call pawcprj_copy(cprj_kg,cprj_ikn)
                    call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                         & dtorbmag%indkk_f2ibz(ikptg,2),dtorbmag%indkk_f2ibz(ikptg,6),&
                         & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptgi)),&
                         & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                         & dtorbmag%mband_occ,my_nspinor,&
                         & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                    call pawcprj_copy(cprj_fkn,cprj_kg)
                 end if

                 if (.NOT. has_smat_indx(ikpt,gdx,0)) then

                    call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &             dtorbmag%lmn_size,dtset%mband,&
                         &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                    sflag_k=0
                    call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
                         &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                         &             pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                    smat_all_indx(:,:,:,ikpt,gdx,0) = smat_kk(:,:,:)
                    smat_all_indx(1,:,:,ikptg,gdxc,0) = TRANSPOSE(smat_kk(1,:,:))
                    smat_all_indx(2,:,:,ikptg,gdxc,0) = -TRANSPOSE(smat_kk(2,:,:))

                    has_smat_indx(ikpt,gdx,0) = .TRUE.
                    has_smat_indx(ikptg,gdxc,0) = .TRUE.

                 end if

                 dkbg = dkg - dkb

                 if (.NOT. has_smat_indx(ikpt,bdx,gdx)) then

                    call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &             dtorbmag%lmn_size,dtset%mband,&
                         &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                    call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                         &             dtorbmag%indkk_f2ibz,ikptb,ikptg,&
                         &             mpi_enreg,npwarr,pwind_bg,symrec)

                    sflag_k=0
                    call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
                         &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                         &             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                    smat_all_indx(:,:,:,ikpt,bdx,gdx) = smat_kk(:,:,:)
                    smat_all_indx(1,:,:,ikpt,gdx,bdx) = TRANSPOSE(smat_kk(1,:,:))
                    smat_all_indx(2,:,:,ikpt,gdx,bdx) = -TRANSPOSE(smat_kk(2,:,:))

                    has_smat_indx(ikpt,bdx,gdx) = .TRUE.
                    has_smat_indx(ikpt,gdx,bdx) = .TRUE.

                 end if

                 IA=czero
                 IB=czero
                 do nn = 1, nband_k
                    do n1 = 1, nband_k

                       t1A = cmplx(smat_all_indx(1,nn,n1,ikpt,bdx,0),smat_all_indx(2,nn,n1,ikpt,bdx,0))
                       t1B = t1A

                       do n2 = 1, nband_k

                          t2A = cmplx(smat_all_indx(1,n1,n2,ikpt,bdx,gdx),smat_all_indx(2,n1,n2,ikpt,bdx,gdx))
                          t3A = conjg(cmplx(smat_all_indx(1,nn,n2,ikpt,gdx,0),smat_all_indx(2,nn,n2,ikpt,gdx,0)))

                          t2B = conjg(cmplx(smat_all_indx(1,n2,n1,ikpt,bdx,0),smat_all_indx(2,n2,n1,ikpt,bdx,0)))

                          do n3 = 1, nband_k

                             t3B = cmplx(smat_all_indx(1,n2,n3,ikpt,gdx,0),smat_all_indx(2,n2,n3,ikpt,gdx,0))
                             t4B=conjg(cmplx(smat_all_indx(1,nn,n3,ikpt,gdx,0),smat_all_indx(2,nn,n3,ikpt,gdx,0)))

                             IB = IB + t1B*t2B*t3B*t4B
                             ! IB = IB - epsabg*bsigma*gsigma*tprodB/(2.0*deltab*2.0*deltag)
                          end do ! end loop over n3

                          IA = IA + t1A*t2A*t3A
                          ! IA = IA + epsabg*bsigma*gsigma*tprodA/(2.0*deltab*2.0*deltag) 

                       end do ! end loop over n2
                    end do ! end loop over n1
                 end do ! end loop over nn

                 cnum(1,adir) = cnum(1,adir) + epsabg*bsigma*gsigma*real(IA-IB)/(2.0*deltab*2.0*deltag) 
                 cnum(2,adir) = cnum(2,adir) + epsabg*bsigma*gsigma*aimag(IA-IB)/(2.0*deltab*2.0*deltag) 

              end do ! end loop over gfor

           end do ! end loop over bfor

        end do ! end loop over epsabg

     end do ! end loop over adir

  end do ! end loop over kpts

  cnum(1,1:3) = ucvol*MATMUL(gprimd,cnum(1,1:3))
  cnum(2,1:3) = ucvol*MATMUL(gprimd,cnum(2,1:3))

  ! factor of 2 in the numerator is the occupation number--each band contains
  ! two electrons, by assumption. Necessary such that trace over density operator
  ! gives number of electrons as expected.
  dtorbmag%chern(1,1:3) = -cnum(2,1:3)*two/(two_pi*ucvol*dtorbmag%fnkpt)
  dtorbmag%chern(2,1:3) =  cnum(1,1:3)*two/(two_pi*ucvol*dtorbmag%fnkpt)
  
  write(message,'(a,a,a)')ch10,'====================================================',ch10
  call wrtout(ab_out,message,'COLL')

  write(message,'(a)')' Chern number C from orbital magnetization '
  call wrtout(ab_out,message,'COLL')
  write(message,'(a,a)')'----C is a real vector, given along Cartesian directions----',ch10
  call wrtout(ab_out,message,'COLL')

  do adir = 1, 3
     write(message,'(a,i4,a,2es16.8)')' C(',adir,') : real, imag ',&
          &   dtorbmag%chern(1,adir),dtorbmag%chern(2,adir)
     call wrtout(ab_out,message,'COLL')
  end do

  write(message,'(a,a,a)')ch10,'====================================================',ch10
  call wrtout(ab_out,message,'COLL')

  if (usepaw == 1) then
     ABI_DEALLOCATE(dimlmn)
     call pawcprj_free(cprj_k)
     ABI_DATATYPE_DEALLOCATE(cprj_k)
     call pawcprj_free(cprj_kb)
     ABI_DATATYPE_DEALLOCATE(cprj_kb)
     call pawcprj_free(cprj_kg)
     ABI_DATATYPE_DEALLOCATE(cprj_kg)
     if (dtset%kptopt /= 3) then
        call pawcprj_free(cprj_ikn)
        call pawcprj_free(cprj_fkn)
        ABI_DATATYPE_DEALLOCATE(cprj_ikn)
        ABI_DATATYPE_DEALLOCATE(cprj_fkn)
     end if
  end if

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(cg1_k)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(pwind_kb)
  ABI_DEALLOCATE(pwind_kg)
  ABI_DEALLOCATE(pwind_bg)
  ABI_DEALLOCATE(pwnsfac_k)

  ABI_DEALLOCATE(has_smat_indx)
  ABI_DEALLOCATE(smat_all_indx)

end subroutine chern_number
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/mpi_chern_number
!! NAME
!! mpi_chern_number
!!
!! FUNCTION
!! This routine computes the Chern number based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtset <type(dataset_type)>=all input variables in this dataset
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! rprimd(3,3) = real space translation vectors
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!! usecprj=1 if cprj datastructure has been allocated
!! usepaw=1 if PAW calculation
!! xred(3,natom) = location of atoms in unit cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006) [[cite:Ceresoli2006]],
!! and Gonze and Zwanziger, PRB 84 064445 (2011) [[cite:Gonze2011a]].
!! This routine computes the Chern number as
!! $C_\alpha = \frac{i}{2\pi}\int_{\mathrm{BZ}} dk \epsilon_{\alpha\beta\gamma}
!! \mathrm{Tr}[\rho_k \partial_\beta \rho_k (1 - \rho_k) \partial_gamma\rho_k] $
!! The derivative of the density operator is obtained from a discretized formula
!! $\partial_\beta \rho_k = \frac{1}{2\Delta}(\rho_{k+b} - \rho_{k-b})$ with
!! $\Delta = |b|$. When reduced to wavefunction overlaps the computation amounts to
!! multiple calls to smatrix.F90, exactly as in other Berry phase computations, with
!! the one additional complication of overlaps like $\langle u_{n,k+b}|u_{n',k+g}\rangle$.
!! At this stage mkpwind_k is invoked, which generalizes the code in initberry
!! and initorbmag necessary to index plane waves around different k points.
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!      m_afterscfloop
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine mpi_chern_number(atindx1,cg,cprj,dtset,dtorbmag,kg,&
     &            mcg,mcprj,mpi_enreg,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,&
     &            rprimd,symrec,usecprj,usepaw,xred)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,pwind_alloc,usecprj,usepaw
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: cg(2,mcg),rprimd(3,3),xred(3,dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,bdx,bdxc,bfor,bsigma,epsabg,gdir,gdx,gdxstor,gfor,gsigma
  integer :: ikpt,isppol
  integer :: my_nspinor,nband_k,nn,n1,n2,n3
  real(dp) :: deltab,deltag,ucvol
  complex(dpc) :: IA,IB,t1A,t2A,t3A,t1B,t2B,t3B,t4B
  character(len=500) :: message
  !arrays
  real(dp) :: cnum(2,3),dkb(3),dkg(3),gmet(3,3),gprimd(3,3),rmet(3,3)
  real(dp),allocatable :: smat_all_indx(:,:,:,:,:,:)

  ! ***********************************************************************
  ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

  ! ! check density operator norm
  ! call rho_norm_check(atindx1,cg,cprj,dtorbmag,dtset,mpi_enreg,mcg,mcprj,&
  !    & npwarr,pawtab,usecprj,usepaw)

  ! TODO: generalize to nsppol > 1
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  nband_k = dtorbmag%mband_occ

  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

  ! the smat_all_indx structure holds the <u_nk|S|u_n'k'> overlap matrix
  ! elements. k ranges over the k pts in the FBZ.
  ! k' can be k + bsigma*dkb, so k +/- an increment in the b direction,
  ! where b ranges over bdir = 1 .. 3. In these cases (equivalent to berryphase_new.F90)
  ! storage is in smat_all_indx(1:2,n,n',ikpt,bdx,0) where bdx maps bdir, bfor onto
  ! the range 1..6. In addition, we also need twist matrix elements of the form
  ! <u_nk+bsigma*dkb|S|u_n'k+gsigma*dkg>, that is, overlap between two functions
  ! that both are neighbors to ikpt, so depend on both bdir and gdir. But we never need
  ! the case bdir // gdir, so we only need four additional slots, not 6. Thus if
  ! bdir = 1, gdir = 2 or 3 and gdx = 3,4,5,6; if bdir = 2, gdir = 3 or 1 and gdx = 5,6,1,2;
  ! if bdir = 3, gdir = 1 or 2, gdx = 1,2,3,4.
  ! This storage is mapped as gdxstor = mod(gdx+6-2*bdir,6)
  ABI_ALLOCATE(smat_all_indx,(2,nband_k,nband_k,dtorbmag%fnkpt,1:6,0:4))
  call make_smat(atindx1,cg,cprj,dtorbmag,dtset,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
     & nband_k,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,smat_all_indx,symrec,xred)

  cnum(:,:) = zero
  do adir = 1, 3
     do epsabg = 1, -1, -2
        if (epsabg .EQ. 1) then
           bdir = modulo(adir,3)+1
           gdir = modulo(adir+1,3)+1
        else
           bdir = modulo(adir+1,3)+1
           gdir = modulo(adir,3)+1
        end if
        do bfor = 1, 2
           if (bfor .EQ. 1) then
              bsigma = 1
           else
              bsigma = -1
           end if
           ! index of neighbor 1..6
           bdx = 2*bdir-2+bfor
           ! index of ikpt viewed from neighbor
           bdxc = 2*bdir-2+bfor+bsigma
           dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
           deltab = sqrt(DOT_PRODUCT(dkb,dkb))
           do gfor = 1, 2
              if (gfor .EQ. 1) then
                 gsigma = 1
              else
                 gsigma = -1
              end if
              ! index of neighbor 1..6
              gdx = 2*gdir-2+gfor
              gdxstor=mod(gdx+6-2*bdir,6)

              dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
              deltag = sqrt(DOT_PRODUCT(dkg,dkg))
              do ikpt = 1, dtorbmag%fnkpt
                 IA=czero
                 IB=czero
                 do nn = 1, nband_k
                    do n1 = 1, nband_k
                       t1A = cmplx(smat_all_indx(1,nn,n1,ikpt,bdx,0),smat_all_indx(2,nn,n1,ikpt,bdx,0))
                       t1B = t1A
                       do n2 = 1, nband_k
                          t2A = cmplx(smat_all_indx(1,n1,n2,ikpt,bdx,gdxstor),smat_all_indx(2,n1,n2,ikpt,bdx,gdxstor))
                          t3A = conjg(cmplx(smat_all_indx(1,nn,n2,ikpt,gdx,0),smat_all_indx(2,nn,n2,ikpt,gdx,0)))
                          t2B = conjg(cmplx(smat_all_indx(1,n2,n1,ikpt,bdx,0),smat_all_indx(2,n2,n1,ikpt,bdx,0)))
                          do n3 = 1, nband_k
                             t3B = cmplx(smat_all_indx(1,n2,n3,ikpt,gdx,0),smat_all_indx(2,n2,n3,ikpt,gdx,0))
                             t4B=conjg(cmplx(smat_all_indx(1,nn,n3,ikpt,gdx,0),smat_all_indx(2,nn,n3,ikpt,gdx,0)))
                             IB = IB + t1B*t2B*t3B*t4B
                          end do ! end loop over n3
                          IA = IA + t1A*t2A*t3A
                       end do ! end loop over n2
                    end do ! end loop over n1
                 end do ! end loop over nn

                 cnum(1,adir) = cnum(1,adir) + epsabg*bsigma*gsigma*real(IA-IB)/(2.0*deltab*2.0*deltag) 
                 cnum(2,adir) = cnum(2,adir) + epsabg*bsigma*gsigma*aimag(IA-IB)/(2.0*deltab*2.0*deltag) 

              end do ! end loop over kpts
           end do ! end loop over gfor
        end do ! end loop over bfor
     end do ! end loop over epsabg
  end do ! end loop over adir

  cnum(1,1:3) = ucvol*MATMUL(gprimd,cnum(1,1:3))
  cnum(2,1:3) = ucvol*MATMUL(gprimd,cnum(2,1:3))

  ! factor of 2 in the numerator is the occupation number--each band contains
  ! two electrons, by assumption. Necessary such that trace over density operator
  ! gives number of electrons as expected.
  dtorbmag%chern(1,1:3) = -cnum(2,1:3)*two/(two_pi*ucvol*dtorbmag%fnkpt)
  dtorbmag%chern(2,1:3) =  cnum(1,1:3)*two/(two_pi*ucvol*dtorbmag%fnkpt)
  
  write(message,'(a,a,a)')ch10,'====================================================',ch10
  call wrtout(ab_out,message,'COLL')

  write(message,'(a)')' Chern number C from orbital magnetization '
  call wrtout(ab_out,message,'COLL')
  write(message,'(a,a)')'----C is a real vector, given along Cartesian directions----',ch10
  call wrtout(ab_out,message,'COLL')

  do adir = 1, 3
     write(message,'(a,i4,a,2es16.8)')' C(',adir,') : real, imag ',&
          &   dtorbmag%chern(1,adir),dtorbmag%chern(2,adir)
     call wrtout(ab_out,message,'COLL')
  end do

  write(message,'(a,a,a)')ch10,'====================================================',ch10
  call wrtout(ab_out,message,'COLL')

  ABI_DEALLOCATE(smat_all_indx)

end subroutine mpi_chern_number
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_smat
!! NAME
!! make_smat
!!
!! FUNCTION
!! This routine computes overlaps <u_k'n'|S|u_kn> as needed by orbital
!! magnetization and Chern numbers. It is assumed that only completely
!! filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group (JWZ)
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_smat(atindx1,cg,cprj,dtorbmag,dtset,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
     & nband_k,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,smat_all,symrec,xred)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,nband_k,pwind_alloc
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
  real(dp),intent(out) :: smat_all(2,nband_k,nband_k,dtorbmag%fnkpt,1:6,0:4)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,bdx,bdxstor,bdxc,bfor,bsigma,countb,countjb,countg,countjg,ddkflag,dest,gdir,gdx,gdxstor,gfor,gg,gsigma
  integer :: icg,icgb,icgg,icprjbi,icprjgi,icprji,ierr,ikg,ikpt,ikpt_loc
  integer :: ikptb,ikptbi,ikptg,ikptgi,ikpti,isppol,itrs
  integer :: jcgb,jcgg,jcprjbi,jcprjgi,jkpt,jkptb,jkptbi,jkptg,jkptgi,jsppol
  integer :: job,mcg1_k,me,my_nspinor,n2dim,ncpgr,nproc,nproc0,npw_k,npw_kb,npw_kg,ntotcp,nn,n1
  integer :: shiftbd,sourceb,sourceg,spaceComm,tagb,tagg,usepaw

  !arrays
  integer,allocatable :: dimlmn(:),nattyp_dum(:),pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
  real(dp) :: dkb(3),dkg(3),dkbg(3),dtm_k(2)
  real(dp),allocatable :: buffer1(:),buffer2(:),buffer(:,:),cg1_k(:,:),cgqb(:,:),cgqg(:,:),kk_paw(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_inv(:,:,:),smat_kk(:,:,:)
  logical,allocatable :: has_smat(:,:,:)
  type(pawcprj_type),allocatable :: cprj_buf(:,:),cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)
  type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

  ! ***********************************************************************

  !Init MPI
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  nproc0=nproc-1
  me=mpi_enreg%me_kpt

  ! TODO: generalize to nsppol > 1
  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
  if (dtset%kptopt /= 3) then
     ABI_DATATYPE_ALLOCATE(cprj_ikn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     ABI_DATATYPE_ALLOCATE(cprj_fkn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
     call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
     call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
  end if

  n2dim = dtorbmag%nspinor*nband_k
  ntotcp = n2dim*SUM(dimlmn(:))
  if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(dtset%natom,n2dim))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
  end if

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))
  ABI_ALLOCATE(pwind_kg,(dtset%mpw))
  ABI_ALLOCATE(pwind_bg,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = 1.0_dp ! bra real
  pwnsfac_k(2,:) = 0.0_dp ! bra imag
  pwnsfac_k(3,:) = 1.0_dp ! ket real
  pwnsfac_k(4,:) = 0.0_dp ! ket imag

  mcg1_k = dtset%mpw*nband_k
  ABI_ALLOCATE(cg1_k,(2,mcg1_k))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

  ABI_ALLOCATE(has_smat,(dtorbmag%fnkpt,1:6,0:6))

  ! eventually must generalize to nsppol, nspinor .NE. 1
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  ddkflag = 1
  !itrs = 0 means do not invoke time reversal symmetry in smatrix.F90
  itrs = 0
  job = 1
  shiftbd = 1
  usepaw = 1
  
  smat_all(:,:,:,:,:,:) = zero
  has_smat(:,:,:) = .FALSE.

  do bdir = 1, 3
     do gg = bdir,bdir+1
        gdir=mod(gg,3)+1

        do bfor = 1, 2
           bsigma = -2*bfor+3
           ! index of neighbor 1..6
           bdx = 2*bdir-2+bfor
           ! index of ikpt viewed from neighbor
           bdxc = 2*bdir-2+bfor+bsigma
           dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)        

           do gfor = 1, 2
              gsigma=-2*gfor+3
              ! index of neighbor 1..6
              gdx = 2*gdir-2+gfor
              dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
              dkbg = dkg - dkb

              !    Loop on the values of ikpt_loc and ikpt1 :
              !    ikpt1 is incremented one by one, and number the k points in the FBZ
              !    ikpt1i refer to the k point numbering in the IBZ
              !    ikpt_loc differs from ikpt1 only in the parallel case, and gives
              !    the index of the k point in the FBZ, in the set treated by the present processor
              !    NOTE : in order to allow synchronisation, ikpt_loc contain information about
              !    ikpt AND ISPPOL !
              !    It means that the following loop is equivalent to a double loop :
              !    do isppol = 1, nsppol
              !    do ikpt1 =  1, dtefield%fmkmem
              !
              ! do ikpt_loc = 1, dtorbmag%fmkmem_max*nsppol
              do ikpt_loc = 1, dtorbmag%fmkmem_max
  
                 ikpt=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
                 ! isppol=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,2)

                 ! if this k and spin are for me do it
                 ! if (ikpt1 > 0 .and. isppol > 0) then
                 if (ikpt > 0) then

                    ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
                    icprji = dtorbmag%cprjindex(ikpti,isppol)
                    npw_k = npwarr(ikpti)
                    icg = dtorbmag%cgindex(ikpti,dtset%nsppol)
                    ikg = dtorbmag%fkgindex(ikpt)

                    ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
                    ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)
                    npw_kb = npwarr(ikptbi)
                    pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

                    ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                    ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)
                    npw_kg = npwarr(ikptgi)
                    pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                    call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
                         &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
                    if ( ikpti /= ikpt ) then
                       call pawcprj_copy(cprj_k,cprj_ikn)
                       call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                            & dtorbmag%indkk_f2ibz(ikpt,2),dtorbmag%indkk_f2ibz(ikpt,6),&
                            & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikpti)),&
                            & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                            & dtorbmag%mband_occ,my_nspinor,&
                            & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                       call pawcprj_copy(cprj_fkn,cprj_k)
                    end if

                 end if ! end check that ikpt > 0

                 !      --------------------------------------------------------------------------------
                 !      Communication
                 !      --------------------------------------------------------------------------------
                 if (ikpt > 0 .and. isppol > 0) then ! I currently have a true kpt to use
                    countb = npw_kb*my_nspinor*nband_k
                    if(allocated(cgqb)) then
                       ABI_DEALLOCATE(cgqb)
                    endif
                    ABI_ALLOCATE(cgqb,(2,countb))
                    cgqb = zero
                    sourceb = me
                    if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikptbi,1,nband_k,isppol,me)) then
                       ! I need the datas from someone else
                       sourceb = mpi_enreg%proc_distrb(ikptbi,1,isppol)
                    end if
                    countg = npw_kg*my_nspinor*nband_k
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    end if
                    ABI_ALLOCATE(cgqg,(2,countg))
                    cgqg = zero
                    sourceg = me
                    if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikptgi,1,nband_k,isppol,me)) then
                       ! I need the datas from someone else
                       sourceg = mpi_enreg%proc_distrb(ikptgi,1,isppol)
                    end if
                 else
                    sourceb = -1 ! I do not have a kpt to use
                    sourceg = -1
                 end if

                 do dest=0,nproc0
                    if ((dest.EQ.me) .AND. (ikpt.GT.0) .AND. (isppol.GT.0)) then
                       ! I am destination and I have something to do
                       if(sourceb.EQ.me) then
                          ! I am destination and source for kptb
                          icprjbi = dtorbmag%cprjindex(ikptbi,isppol)
                          icgb = dtorbmag%cgindex(ikptbi,dtset%nsppol)
                          call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjbi,&
                               &         ikptbi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                               &         my_nspinor,dtset%nsppol,0)
                          cgqb(1:2,1:countb) = cg(1:2,icgb+1:icgb+countb)
                       else ! sourceb .NE. me
                          ! receive cgqb (and cprj_kb)
                          tagb = ikptbi + (isppol - 1)*dtset%nkpt
                          call xmpi_recv(cgqb,sourceb,tagb,spaceComm,ierr)
                          call pawcprj_mpi_recv(dtset%natom,n2dim,dimlmn,ncpgr,cprj_kb,sourceb,spaceComm,ierr)
                       end if
                       if(sourceg.EQ.me) then
                          ! I am destination and source for kptg
                          icprjgi = dtorbmag%cprjindex(ikptgi,isppol)
                          icgg = dtorbmag%cgindex(ikptgi,dtset%nsppol)
                          call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjgi,&
                               &           ikptgi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                               &           my_nspinor,dtset%nsppol,0)
                          cgqg(1:2,1:countg) = cg(1:2,icgg+1:icgg+countg)
                       else ! sourceg .NE. me
                          ! receive cgqg (and cprj_kg)
                          tagg = ikptgi + (isppol - 1)*dtset%nkpt
                          call xmpi_recv(cgqg,sourceg,tagg,spaceComm,ierr)
                          call pawcprj_mpi_recv(dtset%natom,n2dim,dimlmn,ncpgr,cprj_kg,sourceg,spaceComm,ierr)
                       end if
                    else if (dest.NE.me) then
                       ! jkpt is the kpt which is being treated by dest
                       ! jsppol is his isppol
                       jkpt = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,1)
                       jsppol = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,2)
                       if (jkpt > 0 .and. jsppol > 0) then ! dest is treating a true kpt

                          jkptb = dtorbmag%ikpt_dk(jkpt,bfor,bdir)
                          jkptbi = dtorbmag%indkk_f2ibz(jkptb,1)
                          jkptg = dtorbmag%ikpt_dk(jkpt,gfor,gdir)
                          jkptgi = dtorbmag%indkk_f2ibz(jkptg,1)

                          if((mpi_enreg%proc_distrb(jkptbi,1,jsppol) == me))  then
                             jcgb = dtorbmag%cgindex(jkptbi,jsppol)
                             jcprjbi=dtorbmag%cprjindex(jkptbi,jsppol)
                             call pawcprj_get(atindx1,cprj_buf,cprj,dtset%natom,1,jcprjbi,jkptbi,0,jsppol,&
                                  & dtset%mband,dtset%mkmem,dtset%natom,dtorbmag%mband_occ,dtorbmag%mband_occ,&
                                  & my_nspinor,dtset%nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
                                  & proc_distrb=mpi_enreg%proc_distrb)
                             tagb = jkptbi + (jsppol - 1)*dtset%nkpt
                             countjb = npwarr(jkptbi)*my_nspinor*nband_k
                             ABI_ALLOCATE(buffer,(2,countjb))
                             buffer(:,1:countjb)  = cg(:,jcgb+1:jcgb+countjb)
                             call xmpi_send(buffer,dest,tagb,spaceComm,ierr)
                             ABI_DEALLOCATE(buffer)
                             call pawcprj_mpi_send(dtset%natom,n2dim,dimlmn,ncpgr,cprj_buf,dest,spaceComm,ierr)
                          end if ! end check that I am his source
                          if((mpi_enreg%proc_distrb(jkptgi,1,jsppol) == me))  then
                             jcgg = dtorbmag%cgindex(jkptgi,jsppol)
                             jcprjgi=dtorbmag%cprjindex(jkptgi,jsppol)
                             call pawcprj_get(atindx1,cprj_buf,cprj,dtset%natom,1,jcprjgi,jkptgi,0,jsppol,&
                                  & dtset%mband,dtset%mkmem,dtset%natom,dtorbmag%mband_occ,dtorbmag%mband_occ,&
                                  & my_nspinor,dtset%nsppol,0,mpicomm=mpi_enreg%comm_kpt,&
                                  & proc_distrb=mpi_enreg%proc_distrb)
                             tagg = jkptgi + (jsppol - 1)*dtset%nkpt
                             countjg = npwarr(jkptgi)*my_nspinor*nband_k
                             ABI_ALLOCATE(buffer,(2,countjg))
                             buffer(:,1:countjg)  = cg(:,jcgg+1:jcgg+countjg)
                             call xmpi_send(buffer,dest,tagg,spaceComm,ierr)
                             ABI_DEALLOCATE(buffer)
                             call pawcprj_mpi_send(dtset%natom,n2dim,dimlmn,ncpgr,cprj_buf,dest,spaceComm,ierr)
                          end if ! end check that I am his source

                       end if ! end check that jkpt > 0 and jsppol > 0

                    end if ! test dest .EQ. me and ikpt .GT. 0

                 end do ! end loop over dest

                 if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps
                    
                    if ( ikptbi /= ikptb ) then
                       call pawcprj_copy(cprj_kb,cprj_ikn)
                       call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                            & dtorbmag%indkk_f2ibz(ikptb,2),dtorbmag%indkk_f2ibz(ikptb,6),&
                            & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptbi)),&
                            & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                            & dtorbmag%mband_occ,my_nspinor,&
                            & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                       call pawcprj_copy(cprj_fkn,cprj_kb)
                    end if
                 
                    if ( ikptgi /= ikptg ) then
                       call pawcprj_copy(cprj_kg,cprj_ikn)
                       call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                            & dtorbmag%indkk_f2ibz(ikptg,2),dtorbmag%indkk_f2ibz(ikptg,6),&
                            & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptgi)),&
                            & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                            & dtorbmag%mband_occ,my_nspinor,&
                            & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                       call pawcprj_copy(cprj_fkn,cprj_kg)
                    end if
                 
                    if (.NOT. has_smat(ikpt,bdx,0)) then
                       call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                            &           dtorbmag%lmn_size,dtset%mband,&
                            &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                       
                       sflag_k=0
                       call smatrix(cg,cgqb,cg1_k,ddkflag,dtm_k,icg,0,itrs,job,nband_k,&
                            &           mcg,countb,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                            &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)
                       
                       smat_all(:,:,:,ikpt,bdx,0) = smat_kk(:,:,:)
                       has_smat(ikpt,bdx,0) = .TRUE.
                       if(sourceb.EQ.me) then
                          smat_all(1,:,:,ikptb,bdxc,0) = TRANSPOSE(smat_kk(1,:,:))
                          smat_all(2,:,:,ikptb,bdxc,0) = -TRANSPOSE(smat_kk(2,:,:))
                          has_smat(ikptb,bdxc,0) = .TRUE.
                       end if
                    end if
                    
                    if (.NOT. has_smat(ikpt,bdx,gdx) .AND. .NOT. has_smat(ikpt,gdx,bdx) ) then
                       
                       call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,&
                            &             dtorbmag%lmn_size,dtset%mband,&
                            &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                       
                       call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                            &             dtorbmag%indkk_f2ibz,ikptb,ikptg,&
                            &             mpi_enreg,npwarr,pwind_bg,symrec)
                       
                       sflag_k=0
                       call smatrix(cgqb,cgqg,cg1_k,ddkflag,dtm_k,0,0,itrs,job,nband_k,&
                            &             countb,countg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                            &             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)
                       
                       gdxstor = mod(gdx+6-2*bdir,6)
                       smat_all(:,:,:,ikpt,bdx,gdxstor) = smat_kk(:,:,:)
                    
                       bdxstor = mod(bdx+6-2*gdir,6)
                       smat_all(1,:,:,ikpt,gdx,bdxstor) = TRANSPOSE(smat_kk(1,:,:))
                       smat_all(2,:,:,ikpt,gdx,bdxstor) = -TRANSPOSE(smat_kk(2,:,:))
                       
                       has_smat(ikpt,bdx,gdx) = .TRUE.
                       has_smat(ikpt,gdx,bdx) = .TRUE.
                       
                    end if

                    if(allocated(cgqb)) then
                       ABI_DEALLOCATE(cgqb)
                    end if
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    end if

                 end if ! end check on ikpt > 0
                 
              end do ! end loop over ikpt_loc
              
           end do ! end loop over gfor

        end do ! end loop over bfor
        
     end do ! end loop over gg

  end do ! end loop over bdir

  !  MPI communicate stuff between everyone
  if (nproc>1) then
     countb = 2*nband_k*nband_k*dtorbmag%fnkpt*6*5
     ABI_ALLOCATE(buffer1,(countb))
     ABI_ALLOCATE(buffer2,(countb))
     buffer1(1:countb) = reshape(smat_all(1:2,1:nband_k,1:nband_k,1:dtorbmag%fnkpt,1:6,0:4),(/countb/))
     call xmpi_sum(buffer1,buffer2,countb,spaceComm,ierr)
     smat_all(1:2,1:nband_k,1:nband_k,1:dtorbmag%fnkpt,1:6,0:4) = &
&     reshape(buffer2(1:countb),(/2,nband_k,nband_k,dtorbmag%fnkpt,6,5/))
     ABI_DEALLOCATE(buffer1)
     ABI_DEALLOCATE(buffer2)
  end if

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cprj_kb)
  ABI_DATATYPE_DEALLOCATE(cprj_kb)
  call pawcprj_free(cprj_kg)
  ABI_DATATYPE_DEALLOCATE(cprj_kg)
  if (dtset%kptopt /= 3) then
     call pawcprj_free(cprj_ikn)
     call pawcprj_free(cprj_fkn)
     ABI_DATATYPE_DEALLOCATE(cprj_ikn)
     ABI_DATATYPE_DEALLOCATE(cprj_fkn)
  end if
  if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
  end if

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(cg1_k)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(pwind_kb)
  ABI_DEALLOCATE(pwind_kg)
  ABI_DEALLOCATE(pwind_bg)
  ABI_DEALLOCATE(pwnsfac_k)

  ABI_DEALLOCATE(has_smat)

end subroutine make_smat
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_onsite_l_k
!! NAME
!! make_onsite_l_k
!!
!! FUNCTION
!! Compute 1/2 <L_R> onsite contribution to orbital magnetization at given k point
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_onsite_l_k(cprj_k,dtset,idir,nband_k,onsite_l_k,pawrad,pawtab)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,nband_k
  complex(dpc),intent(out) :: onsite_l_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iatom,ilmn,il,im,itypat,jlmn,jl,jm,klmn,kln,mesh_size,nn
  real(dp) :: intg
  complex(dpc) :: cpb,cpk,orbl_me

  !arrays
  real(dp),allocatable :: ff(:)

!--------------------------------------------------------------------  

  onsite_l_k = czero
  do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     mesh_size=pawtab(itypat)%mesh_size
     ABI_ALLOCATE(ff,(mesh_size))
     do jlmn=1,pawtab(itypat)%lmn_size
        jl=pawtab(itypat)%indlmn(1,jlmn)
        jm=pawtab(itypat)%indlmn(2,jlmn)
        do ilmn=1,pawtab(itypat)%lmn_size
           il=pawtab(itypat)%indlmn(1,ilmn)
           im=pawtab(itypat)%indlmn(2,ilmn)
           klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
           kln = pawtab(itypat)%indklmn(2,klmn) ! need this for mesh selection below
           ! compute <L_dir>
           call slxyzs(jl,jm,idir,il,im,orbl_me)
           ! compute integral of phi_i*phi_j - tphi_i*tphi_j
           if (abs(orbl_me) > tol8) then
              ff(1:mesh_size)=pawtab(itypat)%phiphj(1:mesh_size,kln) - pawtab(itypat)%tphitphj(1:mesh_size,kln)
              call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
              call simp_gen(intg,ff,pawrad(itypat))
              do nn = 1, nband_k
                 cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn))
                 cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn))
                 onsite_l_k=onsite_l_k+conjg(cpb)*half*orbl_me*intg*cpk
              end do ! end loop over nn
           end if ! end check that |L_dir| > 0, otherwise ignore term
        end do ! end loop over ilmn
     end do ! end loop over jlmn
     ABI_DEALLOCATE(ff)
  end do ! end loop over atoms

end subroutine make_onsite_l_k
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_S1trace_k
!! NAME
!! make_S1trace_k
!!
!! FUNCTION
!! Compute Trace[\rho_0 S^{(1)} \rho_0] in orbital magnetism context
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_S1trace_k(adir,cprj_k,dtset,eeig,nband_k,pawrad,pawtab,S1trace_k)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,nband_k
  complex(dpc),intent(out) :: S1trace_k
  type(dataset_type),intent(in) :: dtset

  !arrays
  real(dp),intent(in) :: eeig(nband_k,nband_k)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,nband_k)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,epsabg,gdir,iatom,ilmn,itypat,jlmn,klmn,nn
  real(dp) :: ENK
  complex(dpc) :: cpb,cpk

  !arrays

!----------------------------------------------------------------


  S1trace_k = czero

  do epsabg = 1, -1, -2

     if (epsabg .EQ. 1) then
        bdir = modulo(adir,3)+1
        gdir = modulo(adir+1,3)+1
     else
        bdir = modulo(adir+1,3)+1
        gdir = modulo(adir,3)+1
     end if

     do nn = 1, nband_k
        ENK = eeig(nn,nn)
        do iatom=1,dtset%natom
           itypat=dtset%typat(iatom)
           do ilmn=1,pawtab(itypat)%lmn_size
              do jlmn=1,pawtab(itypat)%lmn_size
                 klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
                 cpb=cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn))
                 cpk=cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn))
                 S1trace_k=S1trace_k-half*j_dpc*epsabg*ENK*conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
              end do ! end loop over jlmn
           end do ! end loop over ilmn
        end do ! end loop over atoms
     end do ! end loop over bands

  end do ! end loop over epsabg
  
     
end subroutine make_S1trace_k
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_CCIV_k
!! NAME
!! make_CCIV_k
!!
!! FUNCTION
!! Compute Trace[dS_k/db * dS_k/dg * H_k] arising in orbital magnetism context
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_CCIV_k(adir,CCIV_k,dsdk,eeig,nband_k)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,nband_k
  complex(dpc),intent(out) :: CCIV_k

  !arrays
  real(dp),intent(in) :: dsdk(2,nband_k,nband_k,3),eeig(nband_k,nband_k)

  !Local variables -------------------------
  !scalars
  integer :: bdir,epsabg,gdir,nn,n1
  real(dp) :: ENK
  complex(dpc) :: c1,c2

  !arrays

!----------------------------------------------------------------


  CCIV_k = czero

  do epsabg = 1, -1, -2

     if (epsabg .EQ. 1) then
        bdir = modulo(adir,3)+1
        gdir = modulo(adir+1,3)+1
     else
        bdir = modulo(adir+1,3)+1
        gdir = modulo(adir,3)+1
     end if

     do nn = 1, nband_k
        ENK = eeig(nn,nn)

        do n1 = 1, nband_k
           c1 = cmplx(dsdk(1,nn,n1,bdir),dsdk(2,nn,n1,bdir))
           c2 = cmplx(dsdk(1,n1,nn,gdir),dsdk(2,n1,nn,gdir))
           CCIV_k=CCIV_k-half*j_dpc*epsabg*ENK*c1*c2
        end do ! end loop over n1

     end do ! end loop over nn

  end do ! end loop over epsabg
     
end subroutine make_CCIV_k
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_CCIV_k_FD
!! NAME
!! make_CCIV_k_FD
!!
!! FUNCTION
!! Compute Trace[dS_k/db * dS_k/dg * H_k] arising in orbital magnetism context using
!! finite difference approximation for the derivatives
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
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
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_CCIV_k_FD(adir,dtorbmag,dtset,CCIV_k_FD,cprj_kb_k,eeig,gmet,ikpt,nband_k,pawtab)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,ikpt,nband_k
  complex(dpc),intent(out) :: CCIV_k_FD
  type(dataset_type),intent(in) :: dtset
  type(orbmag_type), intent(inout) :: dtorbmag

  !arrays
  real(dp),intent(in) :: eeig(nband_k,nband_k),gmet(3,3)
  type(pawcprj_type),intent(in) :: cprj_kb_k(dtorbmag%fnkpt,6,0:6,dtset%natom,dtorbmag%nspinor*dtset%mband)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,bdx,bfor,bsigma,epsabg,gdir,gdx,gfor,gsigma,iatom,itypat,ilmn,jlmn,klmn,nn,n1
  real(dp) :: deltab,deltag,ENK,sij
  complex(dpc) :: dsdb,dsdb_1,dsdb_2,dsdg,dsdg_1,dsdg_2

  !arrays
  real(dp) :: dkb(3),dkg(3)
!----------------------------------------------------------------


  CCIV_k_FD = czero

  do epsabg = 1, -1, -2

     if (epsabg .EQ. 1) then
        bdir = modulo(adir,3)+1
        gdir = modulo(adir+1,3)+1
     else
        bdir = modulo(adir+1,3)+1
        gdir = modulo(adir,3)+1
     end if

     do bfor=1, 2
        if (bfor .EQ. 1) then
           bsigma = 1
        else
           bsigma = -1
        end if
        bdx = 2*bdir-2+bfor
        dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
        deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))
        do gfor=1, 2
           if (gfor .EQ. 1) then
              gsigma = 1
           else
              gsigma = -1
           end if
           gdx = 2*gdir-2+gfor
           dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
           deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))

           do nn = 1, nband_k
              ENK = eeig(nn,nn)

              do n1 = 1, nband_k

                 dsdb = czero; dsdg = czero
                 do iatom = 1, dtset%natom
                    itypat = dtset%typat(iatom)
                    do ilmn = 1, pawtab(itypat)%lmn_size
                       dsdb_1=cmplx(cprj_kb_k(ikpt,bdx,0,iatom,nn)%cp(1,ilmn),&
                                  & cprj_kb_k(ikpt,bdx,0,iatom,nn)%cp(2,ilmn))
                       dsdg_1=cmplx(cprj_kb_k(ikpt,gdx,0,iatom,n1)%cp(1,ilmn),&
                                  & cprj_kb_k(ikpt,gdx,0,iatom,n1)%cp(2,ilmn))
                       do jlmn = 1, pawtab(itypat)%lmn_size
                          dsdb_2=cmplx(cprj_kb_k(ikpt,bdx,0,iatom,n1)%cp(1,jlmn),&
                                     & cprj_kb_k(ikpt,bdx,0,iatom,n1)%cp(2,jlmn))
                          dsdg_2=cmplx(cprj_kb_k(ikpt,gdx,0,iatom,nn)%cp(1,jlmn),&
                                     & cprj_kb_k(ikpt,gdx,0,iatom,nn)%cp(2,jlmn))
                          klmn = max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
                          sij=pawtab(itypat)%sij(klmn)
                          dsdb=dsdb+conjg(dsdb_1)*dsdb_2*sij
                          dsdg=dsdg+conjg(dsdg_1)*dsdg_2*sij
                       end do ! end loop over jlmn
                    end do ! end loop over ilmn
                 end do ! end loop over iatom
                 CCIV_k_FD=CCIV_K_FD-half*j_dpc*epsabg*bsigma*gsigma*ENK*dsdb*dsdg/(2.0*deltab*2.0*deltag)
              end do ! end loop over n1
           end do ! end loop over nn
        end do ! end loop over gfor
     end do ! end loop over bfor
  end do ! end loop over epsabg
     
end subroutine make_CCIV_k_FD
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/ctocprjb
!! NAME
!! ctocprjb
!!
!! FUNCTION
!! Compute <p_k+b|u_k> cprj's as needed by orbital magnetization,
!! at all k points and all bands
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine ctocprjb(atindx1,cg,cprj_kb_k,dtorbmag,dtset,gmet,gprimd,&
     & istwf_k,kg,mcg,mpi_enreg,nattyp,ncpgr,npwarr,pawtab,psps,rmet,rprimd,ucvol,xred)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: istwf_k,mcg,ncpgr
  real(dp),intent(in) :: ucvol
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(dtset%ntypat),npwarr(dtset%nkpt)
  real(dp),intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),xred(3,dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  type(pawcprj_type),intent(inout) :: cprj_kb_k(dtorbmag%fnkpt,6,0:6,dtset%natom,dtorbmag%nspinor*dtset%mband)

  !Locals------------------------------------
  !scalars
  integer :: bdir,bdx,bfor,bsigma,choice,cpopt,dimffnl,dimph1d,gdir,gdx,gfor,gsigma
  integer :: ia,iband,icg,ider,idir,ikg,ikpt
  integer :: n1,n2,n3,nband_k,nkpg,npw_k,optder
  real(dp) :: arg

  !arrays
  integer :: nband_dum(1),npwarr_dum(1)
  integer,allocatable :: dimlmn(:),kg_k(:,:)
  real(dp) :: dkb(3),kpoint(3),kpointb(3),kptns(3,1)
  real(dp),allocatable :: cwavef(:,:),ffnl(:,:,:,:),kpg_k(:,:)
  real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:)
  real(dp),allocatable :: ylm_k(:,:),ylm_k_gr(:,:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

  ! ***********************************************************************

  nband_k = dtorbmag%mband_occ

  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')

  ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

  ABI_ALLOCATE(phkxred,(2,dtset%natom))

  n1=dtset%ngfft(1); n2=dtset%ngfft(2); n3=dtset%ngfft(3)
  dimph1d=dtset%natom*(2*(n1+n2+n3)+3)
  ABI_ALLOCATE(ph1d,(2,dimph1d))
  call getph(atindx1,dtset%natom,n1,n2,n3,ph1d,xred)

  do ikpt=1,dtorbmag%fnkpt

     kpoint(:)=dtorbmag%fkptns(:,ikpt)
     npw_k = npwarr(ikpt)
     ABI_ALLOCATE(cwavef,(2,npw_k))

     dimffnl=1 ! 1 + number of derivatives
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))

     icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

     ikg = dtorbmag%fkgindex(ikpt)
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     nkpg = 3
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))

     ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))

     ! data for initylmg call below
     optder=0 ! do not need gradients
     nband_dum(1) = nband_k
     npwarr_dum(1) = npw_k
     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
     ABI_ALLOCATE(ylm_k_gr,(npw_k,3+6*(optder/2),psps%mpsang*psps%mpsang))

     do bdir=1, 3
        do bfor=1, 2

           if (bfor .EQ. 1) then
              bsigma = 1
           else
              bsigma = -1
           end if

           bdx = 2*bdir-2+bfor

           dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)

           kpointb(:) = kpoint(:) + dkb(:)

           do ia=1, dtset%natom
              arg=two_pi*(kpointb(1)*xred(1,ia)+kpointb(2)*xred(2,ia)+kpointb(3)*xred(3,ia))
              phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
           end do

           call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

           call mkkpg(kg_k,kpg_k,kpointb,nkpg,npw_k)

           kptns(:,1) = kpointb(:)
           call initylmg(gprimd,kg_k,kptns,1,mpi_enreg,psps%mpsang,npw_k,&
                & nband_dum,1,npwarr_dum,dtset%nsppol,optder,rprimd,ylm_k,ylm_k_gr)

           !      Compute nonlocal form factors ffnl at all (k+b+G_k):
           ider=0 ! no derivatives
           idir=0 ! not applicable when ider = 0
           call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpointb,psps%lmnmax,&
                & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                & psps%usepaw,psps%useylm,ylm_k,ylm_k_gr)

           choice=1 ! cprj only, no derivatives
           cpopt=0
           idir=0 ! not applicable for choice 1
           do iband = 1, nband_k
              cwavef(1,1:npw_k) = cg(1,icg+(iband-1)*npw_k+1:icg+iband*npw_k)
              cwavef(2,1:npw_k) = cg(2,icg+(iband-1)*npw_k+1:icg+iband*npw_k)

              call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
                   & idir,psps%indlmn,istwf_k,kg_k,kpg_k,kpointb,psps%lmnmax,&
                   & dtset%mgfft,mpi_enreg,&
                   & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                   & phkxred,ph1d,ph3d,ucvol,psps%useylm)
              
              call pawcprj_put(atindx1,cwaveprj,cprj_kb_k(ikpt,bdx,0,:,:),dtset%natom,&
                   & iband,0,ikpt,0,1,nband_k,1,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)

           end do ! end loop over bands

           do gdir=1, 3
              if (gdir .EQ. bdir) cycle
              do gfor=1, 2

                 if (gfor .EQ. 1) then
                    gsigma = 1
                 else
                    gsigma = -1
                 end if

                 gdx = 2*gdir-2+gfor

                 dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir) + gsigma*dtorbmag%dkvecs(1:3,gdir)

                 kpointb(:) = kpoint(:) + dkb(:)

                 do ia=1, dtset%natom
                    arg=two_pi*(kpointb(1)*xred(1,ia)+kpointb(2)*xred(2,ia)+kpointb(3)*xred(3,ia))
                    phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
                 end do

                 call mkkpg(kg_k,kpg_k,kpointb,nkpg,npw_k)

                 call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

                 kptns(:,1) = kpointb(:)
                 call initylmg(gprimd,kg_k,kptns,1,mpi_enreg,psps%mpsang,npw_k,&
                      & nband_dum,1,npwarr_dum,dtset%nsppol,optder,rprimd,ylm_k,ylm_k_gr)
           
                 !      Compute nonlocal form factors ffnl at all (k+b+G_k):
                 ider=0 ! no derivatives
                 idir=0 ! not applicable
                 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                      & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpointb,psps%lmnmax,&
                      & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                      & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                      & psps%usepaw,psps%useylm,ylm_k,ylm_k_gr)
           
                 choice=1 ! no derivs
                 cpopt=0
                 idir=0 ! not applicable for choice 1
                 do iband = 1, nband_k
                    cwavef(1,1:npw_k) = cg(1,icg+(iband-1)*npw_k+1:icg+iband*npw_k)
                    cwavef(2,1:npw_k) = cg(2,icg+(iband-1)*npw_k+1:icg+iband*npw_k)

                    call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
                         & idir,psps%indlmn,istwf_k,kg_k,kpg_k,kpointb,psps%lmnmax,&
                         & dtset%mgfft,mpi_enreg,&
                         & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                         & phkxred,ph1d,ph3d,ucvol,psps%useylm)
              
                    call pawcprj_put(atindx1,cwaveprj,cprj_kb_k(ikpt,bdx,gdx,:,:),dtset%natom,&
                         & iband,0,ikpt,0,1,nband_k,1,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)

                 end do ! end loop over bands

              end do ! end loop over gfor
           end do ! end loop over gdir

        end do ! end loop over bfor
     end do ! end loop over bdir

     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm_k_gr)

  end do ! end loop over nkpt

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cwaveprj)
  ABI_DATATYPE_DEALLOCATE(cwaveprj)

  ABI_DEALLOCATE(phkxred)
  ABI_DEALLOCATE(ph1d)

end subroutine ctocprjb
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_dsdk_FD
!! NAME
!! make_dsdk_FD
!!
!! FUNCTION
!! Compute <u_n'k|dS/dk_idir|u_nk> by finite differences
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! cprj_k
!! dtset
!! idir
!! nband_k
!!
!! OUTPUT
!! dsdk(2,nband_k,nband_k)
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
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_dsdk_FD(bdir,cprj_k,cprj_kb_k,dsdk_FD,dtorbmag,dtset,gmet,ikpt,nband_k,pawtab)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,ikpt,nband_k
  type(dataset_type),intent(in) :: dtset
  type(orbmag_type), intent(inout) :: dtorbmag

  !arrays
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out) :: dsdk_FD(2,nband_k,nband_k)
  type(pawcprj_type),intent(in) :: cprj_k(dtset%natom,dtorbmag%nspinor*dtset%mband)
  type(pawcprj_type),intent(in) :: cprj_kb_k(dtorbmag%fnkpt,6,0:6,dtset%natom,dtorbmag%nspinor*dtset%mband)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdx,bfor,bsigma,iatom,ilmn,jlmn,klmn,itypat,nn, n1
  real(dp) :: deltab,sij
  complex(dpc) :: cp,cp1,cp2,cp3,cp4
  !arrays
  real(dp) :: dkb(3)

!--------------------------------------------------------------------  

  dsdk_FD(1:2,1:nband_k,1:nband_k) = zero
  do nn = 1, nband_k
     do n1 = 1, nband_k

        cp=czero
        do iatom = 1, dtset%natom
           itypat = dtset%typat(iatom)
           do ilmn = 1, pawtab(itypat)%lmn_size
              do jlmn = 1, pawtab(itypat)%lmn_size
                 klmn = max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
                 ! cp1=cmplx(cprj_k(iatom,n1)%dcp(1,bdir,ilmn),cprj_k(iatom,n1)%dcp(2,bdir,ilmn))
                 ! cp2=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn))
                 ! cp3=cmplx(cprj_k(iatom,n1)%cp(1,ilmn),cprj_k(iatom,n1)%cp(2,ilmn))
                 ! cp4=cmplx(cprj_k(iatom,nn)%dcp(1,bdir,jlmn),cprj_k(iatom,nn)%dcp(2,bdir,jlmn))
                 ! sij=pawtab(itypat)%sij(klmn)
                 ! cp = cp + sij*(conjg(cp1)*cp2 + conjg(cp3)*cp4)
                 do bfor = 1, 2
                    bsigma = 1
                    if (bfor .EQ. 2) bsigma = -1
                    bdx = 2*bdir-2+bfor
                    dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
                    ! deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))
                    deltab = sqrt(DOT_PRODUCT(dkb,dkb))

                    cp1=cmplx(cprj_kb_k(ikpt,bdx,0,iatom,n1)%cp(1,ilmn),&
                         &    cprj_kb_k(ikpt,bdx,0,iatom,n1)%cp(2,ilmn))
                    cp2=cmplx(cprj_kb_k(ikpt,bdx,0,iatom,nn)%cp(1,jlmn),&
                         &    cprj_kb_k(ikpt,bdx,0,iatom,nn)%cp(2,jlmn))
                    sij=pawtab(itypat)%sij(klmn)
                    cp=cp+conjg(cp1)*cp2*sij*bsigma/(2.0*deltab)

                 end do ! end loop over bfor
              end do ! end loop over jlmn
           end do ! end loop over ilmn
        end do ! end loop over iatom
        dsdk_FD(1,n1,nn) = real(cp)
        dsdk_FD(2,n1,nn) = aimag(cp)
     end do ! end loop over n1
  end do ! end loop over nn

end subroutine make_dsdk_FD
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_dsdk
!! NAME
!! make_dsdk
!!
!! FUNCTION
!! Compute <u_n'k'|dS/dk_idir|u_nk>
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! cprj_k
!! dtset
!! idir
!! nband_k
!!
!! OUTPUT
!! dsdk(2,nband_k,nband_k,0:6) = <u_n'k'|dS/dk_idir|u_nk>
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!! dsdk(2,nband_k,nband_k,0:6) = <u_n'k'|dS/dk_idir|u_nk>
!! where 0 entry is for k' = k and 1-6 entries are for k' = k + dk, with indexing as in berryphase_new
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine make_dsdk(atindx1,cg,cprj_k,dimlmn,dsdk,dtorbmag,dtset,gs_hamk,&
     & idir,ikpt,isppol,mcg,mpi_enreg,my_nspinor,nband_k,ncpgr,npwarr,&
     & pwind,pwind_alloc)

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,ikpt,isppol,mcg,my_nspinor,nband_k,ncpgr,pwind_alloc
  type(orbmag_type), intent(inout) :: dtorbmag
  type(dataset_type),intent(in) :: dtset
  type(gs_hamiltonian_type),intent(in) :: gs_hamk
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),dimlmn(dtset%natom)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: dsdk(2,nband_k,nband_k,0:6)
  type(pawcprj_type),intent(in) ::  cprj_k(dtset%natom,nband_k)

  !Local variables -------------------------
  !scalars
  integer :: bdir,bdx,bfor,choice,cpopt,icg,icgb,ikg,ikptb,ikptbi
  integer :: n1,ndat,nn,nnlout,npw_k,npw_kb,paw_opt,signs,tim_nonlop
  real(dp) :: doti,dotr

  !arrays
  integer,allocatable :: pwind_kb(:)
  real(dp) :: enlout(1),lambda(1),vectout(2,1)
  real(dp),allocatable :: cwaveb(:,:),cwavef(:,:),swavef(:,:)
  type(pawcprj_type),allocatable :: cwaveprj(:,:)

!--------------------------------------------------------------------  

  ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

  npw_k = npwarr(ikpt)
  icg = dtorbmag%cgindex(ikpt,dtset%nsppol)
  ikg = dtorbmag%fkgindex(ikpt)

  ABI_ALLOCATE(cwavef,(2,0:dtset%mpw))
  ABI_ALLOCATE(cwaveb,(2,0:dtset%mpw))
  ABI_ALLOCATE(swavef,(2,0:dtset%mpw))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))

  choice = 5 ! 1st derivative(s) with respect to k wavevector
  cpopt = 4 ! <p_lmn|in> and first derivatives are already in memory
  nnlout = 1 ! dimension of enlout, not used in nonlop call
  lambda = zero ! not used in nonlop call
  ndat = 1 ! number of wavefunctions on which to apply nonlop
  paw_opt = 3 ! paw_opt=3 : PAW overlap matrix (Sij)
  signs = 2 ! if 2, applies the non-local operator to a function in reciprocal space
  tim_nonlop = 0 ! not using timing

  dsdk(1:2,1:nband_k,1:nband_k,0:6) = zero
  
  do nn = 1, nband_k

     call pawcprj_get(atindx1,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
          &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)

     cwavef(1:2,1:npw_k) = cg(1:2,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
     cwavef(1:2,0) = zero

     call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,lambda,mpi_enreg,ndat,nnlout, &
          &                 paw_opt,signs,swavef(1:2,1:npw_k),tim_nonlop,cwavef(1:2,1:npw_k),vectout)
     swavef(1:2,0) = zero
     
     do n1 = 1, nband_k
        
        cwavef(1:2,1:npw_k) = cg(1:2,icg+(n1-1)*npw_k+1:icg+n1*npw_k)

        dsdk(1,n1,nn,0) = DOT_PRODUCT(cwavef(1,1:npw_k),swavef(1,1:npw_k))+DOT_PRODUCT(cwavef(2,1:npw_k),swavef(2,1:npw_k))
        dsdk(2,n1,nn,0) = DOT_PRODUCT(cwavef(1,1:npw_k),swavef(2,1:npw_k))-DOT_PRODUCT(cwavef(2,1:npw_k),swavef(1,1:npw_k))

        do bdir = 1, 3
           ! never need bdir // idir so ignore
           if (bdir .EQ. idir) cycle
           do bfor = 1, 2
              ! index of neighbor 1..6
              bdx = 2*bdir-2+bfor
              ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
              ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)
              npw_kb = npwarr(ikptbi)
              icgb = dtorbmag%cgindex(ikptbi,dtset%nsppol)
              pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

              cwaveb(1:2,1:npw_kb) = cg(1:2,icgb+(n1-1)*npw_kb+1:icgb+n1*npw_kb)
              cwaveb(1:2,0) = zero

              ! overlap_g is computing <swavef|cwaveb>, that is <u_nk|dS/dk|u_n'k'> which is the
              ! conjugate of of the term we want
              call overlap_g(doti,dotr,dtset%mpw,npw_k,npw_kb,dtset%nspinor,pwind_kb,swavef,cwaveb)
              dsdk(1,n1,nn,bdx) = dotr; dsdk(2,n1,nn,bdx) = -doti

           end do ! end loop over bfor
        end do ! end loop over bdir

     end do ! end loop over n1
        
  end do ! end loop over ket bands

  call pawcprj_free(cwaveprj)
  ABI_DATATYPE_DEALLOCATE(cwaveprj)
  ABI_DEALLOCATE(cwavef)
  ABI_DEALLOCATE(cwaveb)
  ABI_DEALLOCATE(swavef)
  ABI_DEALLOCATE(pwind_kb)

end subroutine make_dsdk
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/orbmag
!! NAME
!! orbmag
!!
!! FUNCTION
!! This routine computes the orbital magnetization based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2019 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!! cg(2,mcg)=planewave coefficients of wavefunctions
!! cprj(natom,mcprj*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!! dtset <type(dataset_type)>=all input variables in this dataset
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! nfftf= - PAW only - number of FFT grid points for the "fine" grid (see NOTES at beginning of scfcv)
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawrad(ntypat*psps%usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!! usecprj=1 if cprj datastructure has been allocated
!! vhartr(nfftf)=Hartree potential
!! vpsp(nfftf)=array for holding local psp
!! vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!! xred(3,natom) = location of atoms in unit cell
!! ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!! ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtorbmag <type(orbmag_type)> = variables related to orbital magnetization
!!
!! TODO
!!
!! NOTES
!! See Ceresoli et al, PRB 74, 024408 (2006) [[cite:Ceresoli2006]],
!! and Gonze and Zwanziger, PRB 84, 064445 (2011) [[cite:Gonze2011a]].
!! The derivative of the density operator is obtained from a discretized formula
!! $\partial_\beta \rho_k = \frac{1}{2\Delta}(\rho_{k+b} - \rho_{k-b})$ with
!! $\Delta = |b|$. When reduced to wavefunction overlaps the computation amounts to
!! multiple calls to smatrix.F90, exactly as in other Berry phase computations, with
!! the one additional complication of overlaps like $\langle u_{n,k+b}|u_{n',k+g}\rangle$.
!! At this stage mkpwind_k is invoked, which generalizes the code in initberry
!! and initorbmag necessary to index plane waves around different k points.
!! Direct questions and comments to J Zwanziger
!!
!! PARENTS
!!      m_afterscfloop
!!
!! CHILDREN
!!      ctocprjb,destroy_hamiltonian,fftpac,getghc,init_hamiltonian
!!      load_k_hamiltonian,load_spin_hamiltonian,make_cciv_k,make_dsdk
!!      make_dsdk_fd,make_onsite_l_k,make_s1trace_k,metric,mkffnl,mkkin,mkkpg
!!      mknucdipmom_k,mkpwind_k,overlap_k1k2_paw,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawcprj_get,pawcprj_getdim,pawcprj_symkn,smatrix,transgrid
!!      wrtout
!!
!! SOURCE

subroutine orbmag(atindx1,cg,cprj,dtset,dtorbmag,kg,&
     &            mcg,mcprj,mpi_enreg,nattyp,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,&
     &            pwind,pwind_alloc,rprimd,symrec,usecprj,vhartr,vpsp,vxc,xred,ylm,ylmgr)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,nfftf,pwind_alloc,usecprj
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(orbmag_type), intent(inout) :: dtorbmag
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%ntypat)
 integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,mcg),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfftf),vpsp(nfftf),vxc(nfftf,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdir,bdx,bdxc,bfor,bsigma,cpopt,ddkflag,dimffnl,epsabg
 integer :: gdir,gdx,gdxc,gfor,gsigma
 integer :: iatom,icg,icgb,icgg,icprj,icprji,icprjb,icprjbi,icprjg,icprjgi,ider,idir
 integer :: ikg,ikgb,ikgg,ikpt,ikpti,ikptb,ikptbi,ikptg,ikptgi
 integer :: il,im,ilm,ilmn,ipw,isppol,istwf_k,itrs,itypat
 integer :: jl,jm,jlmn,job,jpw
 integer :: klmn,kln,mesh_size,mcg1_k,my_cpopt,my_nspinor,nband_k,ncpgr,ndat,nkpg,nn,n1,n2,n3
 integer :: ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,npw_k,npw_kb,npw_kg
 integer :: prtvol,shiftbd,sij_opt,tim_getghc,type_calc,type_calc_123
 real(dp) :: deltab,deltag,dkg2,dotr,doti,ENK,EN2K,htpisq,intg,keg,lambda,ucvol
 complex(dpc) :: cdij,cgdijcb,cpb,cpg,cpk
 complex(dpc) :: CCI,CCI_1,CCI_2,CCI_3,CCIV_k
 complex(dpc) :: CCII,CCII_1,CCII_2,CCII_3,CCII_4,CCVV_k
 complex(dpc) :: onsite_l_k,orbl_me
 complex(dpc) :: S1trace_k,VVI,VVI_1,VVI_2
 complex(dpc) :: VVIII,VVIII_1,VVIII_2
 complex(dpc) :: VVII,VVII_1,VVII_2,VVII_3
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk,gs_hamk123
 !arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:),kg_kb(:,:),kg_kg(:,:)
 integer,allocatable :: pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
 real(dp) :: CCIV(2,3),CCVV(2,3),dkb(3),dkg(3),dkbg(3),dtm_k(2),gmet(3,3),gprimd(3,3)
 real(dp) :: kpoint(3),kpointb(3),kpointg(3)
 real(dp) :: onsite_l(2,3),orbmagvec(2,3),rhodum(1),rmet(3,3),S1trace(2,3)
 real(dp),allocatable :: bra(:,:),cg1_k(:,:),cgrvtrial(:,:),cwavef(:,:),ff(:),ffnl(:,:,:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: hmat(:,:,:,:,:,:),kinpw(:),kk_paw(:,:,:),kpg_k(:,:),kpg_k_dummy(:,:)
 real(dp),allocatable :: my_nucdipmom(:,:),ph3d(:,:,:),pwnsfac_k(:,:),smat_all(:,:,:,:,:,:),smat_inv(:,:,:)
 real(dp),allocatable :: dsmatdk_all(:,:,:,:,:,:),dsdk_FD(:,:,:,:,:)
 real(dp),allocatable :: smat_kk(:,:,:),vlocal(:,:,:,:),vtrial(:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 complex(dpc),allocatable :: nucdipmom_k(:)
 logical,allocatable :: has_dsmatdk(:,:,:),has_hmat(:,:,:),has_smat(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kb_k(:,:,:,:,:)
 type(pawcprj_type),allocatable :: cprj_kg(:,:),cwaveprj(:,:)
 type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

 ! ***********************************************************************
 ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 nband_k = dtorbmag%mband_occ

 if (psps%usepaw == 1) then ! cprj allocation
   ncpgr = cprj(1,1)%ncpgr
   ABI_ALLOCATE(dimlmn,(dtset%natom))
   call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')

   ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
   call pawcprj_alloc(cprj_k,ncpgr,dimlmn)

   ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
   call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)

   ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
   call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)

   ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
   call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

   ABI_DATATYPE_ALLOCATE(cprj_kb_k,(dtorbmag%fnkpt,6,0:6,dtset%natom,dtorbmag%nspinor*dtset%mband))
   do ikpt=1,dtorbmag%fnkpt
      do bdx=1, 6
         do gdx = 0, 6
            call pawcprj_alloc(cprj_kb_k(ikpt,bdx,gdx,:,:),ncpgr,dimlmn)
         end do
      end do
   end do

   if (dtset%kptopt /= 3) then
      ABI_DATATYPE_ALLOCATE(cprj_ikn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
      ABI_DATATYPE_ALLOCATE(cprj_fkn,(dtset%natom,dtorbmag%nspinor*dtset%mband))
      call pawcprj_alloc(cprj_ikn,ncpgr,dimlmn)
      call pawcprj_alloc(cprj_fkn,ncpgr,dimlmn)
   end if

 else
   message = ' usepaw /= 1 but orbital magnetization calculation requires PAW '
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(pwind_kb,(dtset%mpw))
 ABI_ALLOCATE(pwind_kg,(dtset%mpw))
 ABI_ALLOCATE(pwind_bg,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 pwnsfac_k(1,:) = 1.0_dp ! bra real
 pwnsfac_k(2,:) = 0.0_dp ! bra imag
 pwnsfac_k(3,:) = 1.0_dp ! ket real
 pwnsfac_k(4,:) = 0.0_dp ! ket imag

 mcg1_k = dtset%mpw*nband_k
 ABI_ALLOCATE(cg1_k,(2,mcg1_k))
 ABI_ALLOCATE(sflag_k,(nband_k))
 ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
 ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

 ABI_ALLOCATE(kinpw,(dtset%mpw))

 ABI_ALLOCATE(my_nucdipmom,(3,dtset%natom))
 my_nucdipmom(:,:) = dtset%nucdipmom(:,:)

 ! input parameters for calls to smatrix.F90
 ddkflag = 1
 istwf_k = 1
 ! itrs = 0 means do not invoke time reversal symmetry in smatrix.F90
 itrs = 0
 job = 1
 shiftbd = 1

 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)

 ! input parameters for calls to getghc at ikpt
 cpopt = -1
 ndat = 1
 prtvol = 0
 sij_opt = 0
 tim_getghc = 0
 ! getghc: type_calc 0 means kinetic, local, nonlocal
 type_calc = 0
 lambda = zero

 htpisq = 0.5_dp*(two_pi)**2

 ABI_ALLOCATE(has_smat,(dtorbmag%fnkpt,0:6,0:6))
 ABI_ALLOCATE(smat_all,(2,nband_k,nband_k,dtorbmag%fnkpt,0:6,0:6))
 has_smat(:,:,:)=.FALSE.
 ABI_ALLOCATE(has_dsmatdk,(dtorbmag%fnkpt,3,0:6))
 ABI_ALLOCATE(dsmatdk_all,(2,nband_k,nband_k,dtorbmag%fnkpt,3,0:6))
 ABI_ALLOCATE(dsdk_FD,(2,nband_k,nband_k,dtorbmag%fnkpt,3))
 has_dsmatdk(:,:,:) = .FALSE.
 ABI_ALLOCATE(has_hmat,(dtorbmag%fnkpt,0:6,0:6))
 ABI_ALLOCATE(hmat,(2,nband_k,nband_k,dtorbmag%fnkpt,0:6,0:6))
 has_hmat(:,:,:) = .FALSE.

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k, needed for computing E_nk
!  call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
! & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
! & paw_ij=paw_ij)
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=my_nucdipmom,&
      & paw_ij=paw_ij)

 !gs_hamk123 is used to apply vlocal in <u_nk1|Hk2|u_mk3>
 ! my_nucdipmom can be used to override the input nuclear dipoles
 ! call init_hamiltonian(gs_hamk123,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
 !      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom)
 call init_hamiltonian(gs_hamk123,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=my_nucdipmom)

 !---------construct local potential------------------
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ! nspden=1 is essentially hard-coded in the following line
 vtrial(1:nfftf,1)=vhartr(1:nfftf)+vxc(1:nfftf,1)+vpsp(1:nfftf)

 ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
 call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)

 ABI_ALLOCATE(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
 call fftpac(isppol,mpi_enreg,dtset%nspden,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)

 ABI_DEALLOCATE(cgrvtrial)
 ABI_DEALLOCATE(vtrial)

 ! same vlocal in both Hamiltonians
 call load_spin_hamiltonian(gs_hamk,isppol,vlocal=vlocal,with_nonlocal=.true.)
 call load_spin_hamiltonian(gs_hamk123,isppol,vlocal=vlocal,with_nonlocal=.false.)

 !------- now local potential is attached to gs_hamk and gs_hamk123 -------------------------

 ! compute the shifted cprj's <p_k+b|u_k>
 call ctocprjb(atindx1,cg,cprj_kb_k,dtorbmag,dtset,gmet,gprimd,&
      & istwf_k,kg,mcg,mpi_enreg,nattyp,ncpgr,npwarr,pawtab,psps,rmet,rprimd,ucvol,xred)

 ! loop over kpts, assuming for now kptopt 3 or 4, nsppol = 1, nspinor = 1
 ! and no parallelism, no symmorphic symmetry elements
 orbmagvec(:,:) = zero
 S1trace(:,:) = zero
 onsite_l(:,:) = zero
 CCIV(:,:) = zero
 CCVV(:,:) = zero
 do ikpt = 1, dtorbmag%fnkpt

    kpoint(:)=dtorbmag%fkptns(:,ikpt)

    ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
    icprji = dtorbmag%cprjindex(ikpti,isppol)

    npw_k = npwarr(ikpti)
    icg = dtorbmag%cgindex(ikpti,dtset%nsppol)

    ikg = dtorbmag%fkgindex(ikpt)
    ABI_ALLOCATE(kg_k,(3,npw_k))
    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

    call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
         &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
    if ( ikpti /= ikpt ) then
       call pawcprj_copy(cprj_k,cprj_ikn)
       call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
            & dtorbmag%indkk_f2ibz(ikpt,2),dtorbmag%indkk_f2ibz(ikpt,6),&
            & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikpti)),&
            & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
            & dtorbmag%mband_occ,my_nspinor,&
            & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
       call pawcprj_copy(cprj_fkn,cprj_k)
    end if

    ! kpoint(:)=dtorbmag%fkptns(:,ikpt)

    ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
    do ilm=1,psps%mpsang*psps%mpsang
       ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
    end do

    ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
    do ilm=1,psps%mpsang*psps%mpsang
       ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
    end do

    !      Compute (1/2) (2 Pi)**2 (k+G)**2:
    kinpw(:) = zero
    call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

    !  Compute (k+G) vectors (only if useylm=1)
    ! original code from vtorho.F90
    ! nkpg=3*optforces*dtset%nloalg(3)
    ! ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
    ! if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
    !    call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
    ! end if
    nkpg = 3
    ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
    call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
    
    !      Compute nonlocal form factors ffnl at all (k+G):
    ider=1 ! want ffnl and 1st derivative
    idir=4 ! d ffnl/ dk_red in all 3 directions
    dimffnl=4 ! 1 + number of derivatives
    ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
    call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
         &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
         &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
         &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
         &         psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
    
    !     compute and load nuclear dipole Hamiltonian at current k point
    if(any(abs(gs_hamk%nucdipmom)>0.0)) then
       if(allocated(nucdipmom_k)) then
          ABI_DEALLOCATE(nucdipmom_k)
       end if
       ABI_ALLOCATE(nucdipmom_k,(npw_k*(npw_k+1)/2))
       call mknucdipmom_k(gmet,kg_k,kpoint,dtset%natom,gs_hamk%nucdipmom,&
            &           nucdipmom_k,npw_k,rprimd,ucvol,xred)
       if(allocated(gs_hamk%nucdipmom_k)) then
          ABI_DEALLOCATE(gs_hamk%nucdipmom_k)
       end if
       ABI_ALLOCATE(gs_hamk%nucdipmom_k,(npw_k*(npw_k+1)/2))
       call load_k_hamiltonian(gs_hamk,nucdipmom_k=nucdipmom_k)
       ABI_DEALLOCATE(nucdipmom_k)
    end if
    
    !      Load k-dependent part in the Hamiltonian datastructure
    !       - Compute 3D phase factors
    !       - Prepare various tabs in case of band-FFT parallelism
    !       - Load k-dependent quantities in the Hamiltonian
    ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
    
    call load_k_hamiltonian(gs_hamk,kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
         &         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k_dummy,ffnl_k=ffnl,ph3d_k=ph3d,&
         &         compute_ph3d=.TRUE.,compute_gbound=.TRUE.)
    
    if (.NOT. has_hmat(ikpt,0,0) ) then
       ! apply gs_hamk to wavefunctions at k to compute E_nk eigenvalues
       ABI_ALLOCATE(cwavef,(2,npw_k))
       ABI_ALLOCATE(ghc,(2,npw_k))
       ABI_ALLOCATE(gsc,(2,npw_k))
       ABI_ALLOCATE(gvnlc,(2,npw_k))
       hmat(:,:,:,ikpt,0,0) = zero
       do nn = 1, nband_k
          cwavef(1,1:npw_k) = cg(1,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
          cwavef(2,1:npw_k) = cg(2,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
          call pawcprj_get(atindx1,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
               &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
          my_cpopt=2
          call getghc(my_cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
               &           prtvol,sij_opt,tim_getghc,type_calc)
          hmat(1,nn,nn,ikpt,0,0)= DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) &
               &           + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))
       end do
       has_hmat(ikpt,0,0) = .TRUE.
       
       ABI_DEALLOCATE(cwavef)
       ABI_DEALLOCATE(ghc)
       ABI_DEALLOCATE(gsc)
       ABI_DEALLOCATE(gvnlc)
    end if ! end check on has_hmat
       
    ABI_DEALLOCATE(ylm_k)
    ABI_DEALLOCATE(ylmgr_k)
    ABI_DEALLOCATE(kpg_k)
    ABI_DEALLOCATE(ffnl)
    if(any(abs(gs_hamk%nucdipmom)>0.0)) then
       if(allocated(nucdipmom_k)) then
          ABI_DEALLOCATE(nucdipmom_k)
       end if
    end if
    ABI_DEALLOCATE(ph3d)

    ! make the <u_n'k|dS/dk|u_nk> terms
    do adir = 1, 3
       if (.NOT. has_dsmatdk(ikpt,adir,0)) then
          call make_dsdk(atindx1,cg,cprj_k,dimlmn,dsmatdk_all(1:2,1:nband_k,1:nband_k,ikpt,adir,0:6),&
               & dtorbmag,dtset,gs_hamk,adir,ikpt,isppol,mcg,mpi_enreg,my_nspinor,nband_k,ncpgr,npwarr,&
               & pwind,pwind_alloc)
          call make_dsdk_FD(adir,cprj_k,cprj_kb_k,dsdk_FD(1:2,1:nband_k,1:nband_k,ikpt,adir),&
               & dtorbmag,dtset,gmet,ikpt,nband_k,pawtab)
       end if
    end do

    do adir = 1, 3

       call make_onsite_l_k(cprj_k,dtset,adir,nband_k,onsite_l_k,pawrad,pawtab)

       call make_S1trace_k(adir,cprj_k,dtset,hmat(1,1:nband_k,1:nband_k,ikpt,0,0),nband_k,pawrad,pawtab,S1trace_k)

       call make_CCIV_k(adir,CCIV_k,dsmatdk_all(1:2,1:nband_k,1:nband_k,ikpt,1:3,0),&
            & hmat(1,1:nband_k,1:nband_k,ikpt,0,0),nband_k)
       ! call make_CCIV_k(adir,CCIV_k,dsdk_FD(1:2,1:nband_k,1:nband_k,ikpt,1:3),&
       !      & hmat(1,1:nband_k,1:nband_k,ikpt,0,0),nband_k)

       CCVV_k = czero
       do epsabg = 1, -1, -2

          if (epsabg .EQ. 1) then
             bdir = modulo(adir,3)+1
             gdir = modulo(adir+1,3)+1
          else
             bdir = modulo(adir+1,3)+1
             gdir = modulo(adir,3)+1
          end if

          do bfor = 1, 2
             if (bfor .EQ. 1) then
                bsigma = 1
             else
                bsigma = -1
             end if
             ! index of neighbor 1..6
             bdx = 2*bdir-2+bfor
             ! index of ikpt viewed from neighbor
             bdxc = 2*bdir-2+bfor+bsigma

             dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
             ! deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))
             deltab = sqrt(DOT_PRODUCT(dkb,dkb))

             ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
             kpointb(:)=dtorbmag%fkptns(:,ikptb)

             ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)
             icprjbi = dtorbmag%cprjindex(ikptbi,isppol)

             ikgb = dtorbmag%fkgindex(ikptb)
             npw_kb = npwarr(ikptbi)

             ABI_ALLOCATE(kg_kb,(3,npw_kb))
             kg_kb(:,1:npw_kb)=kg(:,ikgb+1:ikgb+npw_kb)

             icgb = dtorbmag%cgindex(ikptbi,dtset%nsppol)

             pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

             call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjbi,&
                  &         ikptbi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                  &         my_nspinor,dtset%nsppol,0)
             if ( ikptbi /= ikptb ) then
                call pawcprj_copy(cprj_kb,cprj_ikn)
                call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                     & dtorbmag%indkk_f2ibz(ikptb,2),dtorbmag%indkk_f2ibz(ikptb,6),&
                     & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptbi)),&
                     & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                     & dtorbmag%mband_occ,my_nspinor,&
                     & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                call pawcprj_copy(cprj_fkn,cprj_kb)
             end if

             if (.NOT. has_smat(ikpt,bdx,0)) then

                call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                     &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                sflag_k=0
                call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
                     &           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                     &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                smat_all(:,:,:,ikpt,bdx,0) = smat_kk(:,:,:)
                smat_all(1,:,:,ikptb,bdxc,0) = TRANSPOSE(smat_kk(1,:,:))
                smat_all(2,:,:,ikptb,bdxc,0) = -TRANSPOSE(smat_kk(2,:,:))

                has_smat(ikpt,bdx,0) = .TRUE.
                has_smat(ikptb,bdxc,0) = .TRUE.

             end if

             do gfor = 1, 2
                if (gfor .EQ. 1) then
                   gsigma = 1
                else
                   gsigma = -1
                end if
                ! index of neighbor 1..6
                gdx = 2*gdir-2+gfor
                ! index of ikpt viewed from neighbor
                gdxc = 2*gdir-2+gfor+gsigma

                dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
                ! deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))
                deltag = sqrt(DOT_PRODUCT(dkg,dkg))

                ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)

                icprjgi = dtorbmag%cprjindex(ikptgi,isppol)

                kpointg(:)=dtorbmag%fkptns(:,ikptg)

                npw_kg = npwarr(ikptgi)
                ABI_ALLOCATE(kg_kg,(3,npw_kg))
                ikgg = dtorbmag%fkgindex(ikptg)
                kg_kg(:,1:npw_kg)=kg(:,ikgg+1:ikgg+npw_kg)

                icgg = dtorbmag%cgindex(ikptgi,dtset%nsppol)

                pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjgi,&
                     &           ikptgi,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                     &           my_nspinor,dtset%nsppol,0)
                if ( ikptgi /= ikptg ) then
                   call pawcprj_copy(cprj_kg,cprj_ikn)
                   call pawcprj_symkn(cprj_fkn,cprj_ikn,dtorbmag%atom_indsym,dimlmn,-1,psps%indlmn,&
                        & dtorbmag%indkk_f2ibz(ikptg,2),dtorbmag%indkk_f2ibz(ikptg,6),&
                        & dtorbmag%fkptns(:,dtorbmag%i2fbz(ikptgi)),&
                        & dtorbmag%lmax,dtorbmag%lmnmax,dtset%mband,dtset%natom,&
                        & dtorbmag%mband_occ,my_nspinor,&
                        & dtorbmag%nsym,dtset%ntypat,dtset%typat,dtorbmag%zarot)
                   call pawcprj_copy(cprj_fkn,cprj_kg)
                end if

                if (.NOT. has_smat(ikpt,gdx,0)) then

                   call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
                        &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                        &             pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   smat_all(:,:,:,ikpt,gdx,0) = smat_kk(:,:,:)
                   smat_all(1,:,:,ikptg,gdxc,0) = TRANSPOSE(smat_kk(1,:,:))
                   smat_all(2,:,:,ikptg,gdxc,0) = -TRANSPOSE(smat_kk(2,:,:))

                   has_smat(ikpt,gdx,0) = .TRUE.
                   has_smat(ikptg,gdxc,0) = .TRUE.

                end if

                dkbg = dkg - dkb

                if (.NOT. has_smat(ikpt,bdx,gdx)) then

                   call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                        &             dtorbmag%indkk_f2ibz,ikptb,ikptg,&
                        &             mpi_enreg,npwarr,pwind_bg,symrec)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
                        &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                        &             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   smat_all(:,:,:,ikpt,bdx,gdx) = smat_kk(:,:,:)
                   smat_all(1,:,:,ikpt,gdx,bdx) = TRANSPOSE(smat_kk(1,:,:))
                   smat_all(2,:,:,ikpt,gdx,bdx) = -TRANSPOSE(smat_kk(2,:,:))

                   has_smat(ikpt,bdx,gdx) = .TRUE.
                   has_smat(ikpt,gdx,bdx) = .TRUE.

                end if

                if (.NOT. has_hmat(ikpt,gdx,bdx)) then

                   call mkpwind_k(-dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                        &             dtorbmag%indkk_f2ibz,ikptg,ikptb,&
                        &             mpi_enreg,npwarr,pwind_bg,symrec)

                   nkpg = 0
                   ABI_ALLOCATE(kpg_k_dummy,(npw_kb,nkpg))

                   !     compute and load nuclear dipole Hamiltonian at current k point
                   ! this may need to be modified to take into account "twist"
                   if(any(abs(gs_hamk123%nucdipmom)>0.0)) then
                      if(allocated(nucdipmom_k)) then
                         ABI_DEALLOCATE(nucdipmom_k)
                      end if
                      ABI_ALLOCATE(nucdipmom_k,(npw_kb*(npw_kb+1)/2))
                      call mknucdipmom_k(gmet,kg_kb,kpointb,dtset%natom,gs_hamk123%nucdipmom,&
                           &           nucdipmom_k,npw_kb,rprimd,ucvol,xred)
                      if(allocated(gs_hamk123%nucdipmom_k)) then
                         ABI_DEALLOCATE(gs_hamk123%nucdipmom_k)
                      end if
                      ABI_ALLOCATE(gs_hamk123%nucdipmom_k,(npw_kb*(npw_kb+1)/2))
                      call load_k_hamiltonian(gs_hamk123,nucdipmom_k=nucdipmom_k)
                      ABI_DEALLOCATE(nucdipmom_k)
                   end if

                   ! this is minimal Hamiltonian information, to apply vlocal (and only vlocal) to |u_kb>
                   call load_k_hamiltonian(gs_hamk123,kpt_k=kpointb(:),istwf_k=istwf_k,npw_k=npw_kb,&
                        &             kg_k=kg_kb,kpg_k=kpg_k_dummy,compute_gbound=.TRUE.)


                   ! apply gs_hamk123 to wavefunctions at kb
                   ABI_ALLOCATE(cwavef,(2,npw_kb))
                   ABI_ALLOCATE(ghc,(2,npw_kb))
                   ABI_ALLOCATE(gsc,(2,npw_kb))
                   ABI_ALLOCATE(gvnlc,(2,npw_kb))

                   ABI_ALLOCATE(bra,(2,npw_kg))

                   ! getghc: type_calc_123 1 means local only, 3 means kinetic, local only
                   type_calc_123 = 1
                   dkg2=DOT_PRODUCT(dkg(:),MATMUL(gmet(:,:),dkg(:)))

                   do nn = 1, nband_k
                      cwavef(1,1:npw_kb) = cg(1,icgb+(nn-1)*npw_kb+1:icgb+nn*npw_kb)
                      cwavef(2,1:npw_kb) = cg(2,icgb+(nn-1)*npw_kb+1:icgb+nn*npw_kb)
                      ! apply only vlocal
                      call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk123,gvnlc,lambda,mpi_enreg,ndat,&
                           &               prtvol,sij_opt,tim_getghc,type_calc_123)
                      do n1 = 1, nband_k
                         bra(1,1:npw_kg) = cg(1,icgg+(n1-1)*npw_kg+1:icgg+n1*npw_kg)
                         bra(2,1:npw_kg) = cg(2,icgg+(n1-1)*npw_kg+1:icgg+n1*npw_kg)
                         dotr=zero;doti=zero

                         ! apply local potential through ghc
                         do ipw=1,npw_kg
                            jpw=pwind_bg(ipw)
                            if(jpw .GT. 0) then
                               dotr=dotr+bra(1,ipw)*ghc(1,jpw)+bra(2,ipw)*ghc(2,jpw)
                               doti=doti+bra(1,ipw)*ghc(2,jpw)-bra(2,ipw)*ghc(1,jpw)

                               ! twisted kinetic energy: here we are computing
                               ! -\frac{1}{2}<u_kg|e^{-i.k.r}\nabla^2 e^{i.k.r}|u_kb>,
                               ! that is, kinetic energy at k between wavefunctions at kg and kb. The correct
                               ! formula is htpisq*(ikpt + G_right)^2\delta(G_left,G_right) but it's hard to apply
                               ! because the G's have wrap-around shifts (output of mkpwind_k)
                               ! for the kpts near and at the edge of the IBZ.
                               ! The following approach is based on the bra <u_kg| expansion, because
                               ! these G vectors are unshifted (indexed by ipw, not jpw). So we are using
                               ! k+G_left = (k-kg) + (kg+G_left) = -dkg + (kg+G_left). When squared we obtain
                               ! |kg+G_left|^2 - 2*dkg*(kg+G_left) + |dkg|^2. In this way we only use the G_left
                               ! expansion vectors, with no shift, for each k point.

                               ! normal kinetic energy for bra
                               keg=htpisq*dot_product((kpointg(:)+kg_kg(:,ipw)),MATMUL(gmet,(kpointg(:)+kg_kg(:,ipw))))

                               ! addition of |dkg|^2
                               keg=keg+htpisq*dkg2

                               ! addition of -2*dkg*(kg+G_left)
                               keg=keg-2.0*htpisq*DOT_PRODUCT(dkg(:),MATMUL(gmet,(kpointg(:)+kg_kg(:,ipw))))

                               ! application of ecut filter and wavefunction
                               if (keg < dtset%ecut) then
                                  dotr=dotr+bra(1,ipw)*keg*cwavef(1,jpw)+bra(2,ipw)*keg*cwavef(2,jpw)
                                  doti=doti+bra(1,ipw)*keg*cwavef(2,jpw)-bra(2,ipw)*keg*cwavef(1,jpw)
                               end if ! end keg filter
                            end if ! end check on jpw > 0
                         end do ! end loop over ipw

                         ! apply onsite terms through cprjk+b
                         cgdijcb = czero
                         do iatom = 1, dtset%natom
                            itypat = dtset%typat(iatom)
                            do ilmn = 1, pawtab(itypat)%lmn_size
                               cpg=cmplx(cprj_kb_k(ikptg,gdxc,0,iatom,n1)%cp(1,ilmn),&
                                    & cprj_kb_k(ikptg,gdxc,0,iatom,n1)%cp(2,ilmn))
                               do jlmn = 1, pawtab(itypat)%lmn_size
                                  cpb=cmplx(cprj_kb_k(ikptb,bdxc,0,iatom,nn)%cp(1,jlmn),&
                                       & cprj_kb_k(ikptb,bdxc,0,iatom,nn)%cp(2,jlmn))
                                  if (jlmn .LE. ilmn) then
                                     klmn = (ilmn-1)*ilmn/2 + jlmn
                                  else
                                     klmn = (jlmn-1)*jlmn/2 + ilmn
                                  end if
                                  if (paw_ij(iatom)%cplex_dij .EQ. 2) then
                                     cdij=cmplx(paw_ij(iatom)%dij(2*klmn-1,1),paw_ij(iatom)%dij(2*klmn,1))
                                     if (jlmn .GT. ilmn) cdij=conjg(cdij)
                                  else
                                     cdij=cmplx(paw_ij(iatom)%dij(klmn,1),zero)
                                  end if
                                  cgdijcb = cgdijcb + conjg(cpg)*cdij*cpb
                               end do
                            end do
                         end do
                         hmat(1,n1,nn,ikpt,gdx,bdx) = dotr + real(cgdijcb)
                         hmat(2,n1,nn,ikpt,gdx,bdx) = doti + aimag(cgdijcb)
                         ! hmat(1,n1,nn,ikpt,gdx,bdx) = dotr
                         ! hmat(2,n1,nn,ikpt,gdx,bdx) = doti
                         ! hmat(1,nn,n1,ikpt,bdx,gdx) = dotr
                         ! hmat(2,nn,n1,ikpt,bdx,gdx) = -doti
                      end do ! end loop over n1
                   end do ! end loop over nn
                   has_hmat(ikpt,gdx,bdx) = .TRUE.
                   ! has_hmat(ikpt,bdx,gdx) = .TRUE.

                   ABI_DEALLOCATE(cwavef)
                   ABI_DEALLOCATE(bra)
                   ABI_DEALLOCATE(ghc)
                   ABI_DEALLOCATE(gsc)
                   ABI_DEALLOCATE(gvnlc)

                   ABI_DEALLOCATE(kpg_k_dummy)

                end if

                VVI = czero; VVII = czero; VVIII=czero
                CCI = czero; CCII = czero
                do nn = 1, nband_k
                   ENK = hmat(1,nn,nn,ikpt,0,0)

                   do n1 = 1, nband_k

                      VVI_1 = cmplx(smat_all(1,nn,n1,ikpt,bdx,0),smat_all(2,nn,n1,ikpt,bdx,0))
                      VVI_2 = cmplx(dsmatdk_all(1,n1,nn,ikpt,gdir,bdx),dsmatdk_all(2,n1,nn,ikpt,gdir,bdx))

                      VVIII_1 = cmplx(dsmatdk_all(1,n1,nn,ikpt,bdir,gdx),dsmatdk_all(2,n1,nn,ikpt,bdir,gdx))
                      VVIII_2 = cmplx(smat_all(1,nn,n1,ikpt,gdx,0),smat_all(2,nn,n1,ikpt,gdx,0))

                      VVII_1 = cmplx(smat_all(1,nn,n1,ikpt,bdx,0),smat_all(2,nn,n1,ikpt,bdx,0))
                      
                      CCI_1 = cmplx(smat_all(1,nn,n1,ikpt,gdx,0),smat_all(2,nn,n1,ikpt,gdx,0))

                      ! CCII_1 = cmplx(smat_all(1,nn,n1,ikpt,bdx,0),smat_all(2,nn,n1,ikpt,bdx,0))

                      do n2 = 1, nband_k
                         VVII_2 = cmplx(smat_all(1,n1,n2,ikpt,bdx,gdx),smat_all(2,n1,n2,ikpt,bdx,gdx))
                         VVII_3 = cmplx(smat_all(1,n2,nn,ikptg,gdxc,0),smat_all(2,n2,nn,ikptg,gdxc,0))
                         
                         CCI_2 = cmplx(hmat(1,n1,n2,ikpt,gdx,bdx),hmat(2,n1,n2,ikpt,gdx,bdx))
                         CCI_3 = cmplx(smat_all(1,n2,nn,ikptb,bdxc,0),smat_all(2,n2,nn,ikptb,bdxc,0))

                         ! CCII_2 = cmplx(smat_all(1,n1,n2,ikptb,bdxc,0),smat_all(2,n1,n2,ikptb,bdxc,0))

                         ! do n3 = 1, nband_k

                         !    CCII_3 = cmplx(smat_all(1,n2,n3,ikpt,gdx,0),smat_all(2,n2,n3,ikpt,gdx,0))
                         !    CCII_4 = cmplx(smat_all(1,n3,nn,ikptg,gdxc,0),smat_all(2,n3,nn,ikptg,gdxc,0))
                         !    CCII = CCII - ENK*CCII_1*CCII_2*CCII_3*CCII_4

                         ! end do ! end n3

                         CCI = CCI + CCI_1*CCI_2*CCI_3

                         VVII = VVII + ENK*VVII_1*VVII_2*VVII_3

                      end do ! end n2

                      VVI = VVI + ENK*VVI_1*VVI_2

                      VVIII = VVIII + ENK*conjg(VVIII_1)*conjg(VVIII_2)

                   end do ! end n1

                end do ! end nn

                CCVV_k = CCVV_k - half*j_dpc*epsabg*bsigma*(-half*VVI)/(2.0*deltab) 
                CCVV_k = CCVV_k - half*j_dpc*epsabg*gsigma*(-half*VVIII)/(2.0*deltag) ! VVI and VVIII are not good
                CCVV_k = CCVV_k - half*j_dpc*epsabg*bsigma*gsigma*(CCI-VVII)/(2.0*deltab*2.0*deltag) ! 

                ABI_DEALLOCATE(kg_kg)

             end do ! end gfor

             ABI_DEALLOCATE(kg_kb)

          end do ! end bfor

       end do ! end loop over epsabg

       ! orbmagvec(1,adir) = orbmagvec(1,adir) + real(CCVV_k)
       ! orbmagvec(2,adir) = orbmagvec(2,adir) + aimag(CCVV_k)
       CCVV(1,adir) = CCVV(1,adir) + real(CCVV_k)
       CCVV(2,adir) = CCVV(2,adir) + aimag(CCVV_k)

       CCIV(1,adir) = CCIV(1,adir) + real(CCIV_k)
       CCIV(2,adir) = CCIV(2,adir) + aimag(CCIV_k)

       S1trace(1,adir) = S1trace(1,adir) - real(S1trace_k)
       S1trace(2,adir) = S1trace(2,adir) - aimag(S1trace_k)

       onsite_l(1,adir) = onsite_l(1,adir) + real(onsite_l_k)
       onsite_l(2,adir) = onsite_l(2,adir) + aimag(onsite_l_k)

    end do ! end loop over adir
    
    ABI_DEALLOCATE(kg_k)

 end do ! end loop over fnkpt

 ! convert terms to cartesian coordinates as needed
 ! note that terms like <dv/dk| x |dw/dk> computed in reduced coords,
 ! become ucvol*gprimd*<dv/dk| x |dw/dk> when expressed in cartesian coords

 CCVV(1,1:3) = ucvol*MATMUL(gprimd,CCVV(1,1:3))
 CCVV(2,1:3) = ucvol*MATMUL(gprimd,CCVV(2,1:3))
 
 S1trace(1,1:3) = ucvol*MATMUL(gprimd,S1trace(1,1:3))
 S1trace(2,1:3) = ucvol*MATMUL(gprimd,S1trace(2,1:3))

 CCIV(1,1:3) = ucvol*MATMUL(gprimd,CCIV(1,1:3))
 CCIV(2,1:3) = ucvol*MATMUL(gprimd,CCIV(2,1:3))

 ! onsite_l is already cartesian
 
 ! accumulate in orbmagvec
 ! terms are: CCI to CCIV, VVI to VVIII, S1trace, and onsite_l.
 ! Signs appear as CCI - CCII - CCIII + CCIV - VVI - VVII - VVIII
 ! also, |CCII| = |CCIII| = |CCIV|
 ! CCVV includes CCI and all VV terms

 orbmagvec(1:2,1:3) = CCVV(1:2,1:3) - CCIV(1:2,1:3) +  S1trace(1:2,1:3) + onsite_l(1:2,1:3)

 ! orbmagvec(1,1:3) = CCVV(1,1:3)
 ! orbmagvec(2,1:3) = CCVV(2,1:3)
 ! orbmagvec(1,1:3) = S1trace(1,1:3)
 ! orbmagvec(2,1:3) = S1trace(2,1:3)
 ! orbmagvec(1,1:3) = onsite_l(1,1:3)
 ! orbmagvec(2,1:3) = onsite_l(2,1:3)
 ! orbmagvec(1,1:3) = CCIV(1,1:3)
 ! orbmagvec(2,1:3) = CCIV(2,1:3)
 
 ! pre factor is occ/ucvol*N_k
 ! factor of 2 in numerator is the band occupation (two electrons in normal insulator)
 ! converting integral over k space to a sum gives a factor of Omega_BZ/N_k or 1/ucvol*N_k

 dtorbmag%orbmagvec(1:2,1:3) = two*orbmagvec(1:2,1:3)/(ucvol*dtorbmag%fnkpt)

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 write(message,'(a)')' Orbital magnetization '
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,a)')'----Orbital magnetization is a real vector, given along Cartesian directions----',ch10
 call wrtout(ab_out,message,'COLL')

 do adir = 1, 3
   write(message,'(a,i4,a,2es16.8)')' Orb Mag(',adir,') : real, imag ',&
&   dtorbmag%orbmagvec(1,adir),dtorbmag%orbmagvec(2,adir)
   call wrtout(ab_out,message,'COLL')
 end do

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')


 if (psps%usepaw == 1) then
   ABI_DEALLOCATE(dimlmn)
   call pawcprj_free(cprj_k)
   ABI_DATATYPE_DEALLOCATE(cprj_k)
   call pawcprj_free(cprj_kb)
   ABI_DATATYPE_DEALLOCATE(cprj_kb)
   call pawcprj_free(cprj_kg)
   ABI_DATATYPE_DEALLOCATE(cprj_kg)
   call pawcprj_free(cwaveprj)
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
   do ikpt=1,dtorbmag%fnkpt
      do bdx = 1, 6
         do gdx = 0, 6
            call pawcprj_free(cprj_kb_k(ikpt,bdx,gdx,:,:))
         end do
      end do
   end do
   ABI_DATATYPE_DEALLOCATE(cprj_kb_k)
 end if

 ABI_DEALLOCATE(kk_paw)
 ABI_DEALLOCATE(cg1_k)
 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_inv)
 ABI_DEALLOCATE(smat_kk)
 ABI_DEALLOCATE(pwind_kb)
 ABI_DEALLOCATE(pwind_kg)
 ABI_DEALLOCATE(pwind_bg)
 ABI_DEALLOCATE(pwnsfac_k)

 ABI_DEALLOCATE(kinpw)

 ABI_DEALLOCATE(has_smat)
 ABI_DEALLOCATE(smat_all)
 ABI_DEALLOCATE(has_hmat)
 ABI_DEALLOCATE(hmat)
 ABI_DEALLOCATE(has_dsmatdk)
 ABI_DEALLOCATE(dsmatdk_all)

 ABI_DEALLOCATE(my_nucdipmom)

 ABI_DEALLOCATE(vlocal)
 call destroy_hamiltonian(gs_hamk)
 call destroy_hamiltonian(gs_hamk123)

end subroutine orbmag
!!***


end module m_orbmag
!!***
