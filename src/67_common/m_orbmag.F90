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
!! Copyright (C) 2011-2017 ABINIT group (JWZ)
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
  use m_fft,              only : fftpac
  use m_fftcore,          only : kpgsph
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian,destroy_hamiltonian,&
       &                         load_spin_hamiltonian,load_k_hamiltonian,gs_hamiltonian_type
  use m_initylmg,         only : initylmg
  use m_kg,               only : getph,mkkin,mkpwind_k,mknucdipmom_k,ph1d3d
  use m_kpts,             only : listkk, smpbz
  use m_mkffnl,           only : mkffnl
  use m_mpinfo,           only : proc_distrb_cycle
  use m_pawang,           only : pawang_type
  use m_pawfgr,           only : pawfgr_type
  use m_paw_ij,           only : paw_ij_type
  use m_paw_overlap,      only : overlap_k1k2_paw
  use m_pawrad,           only : pawrad_type
  use m_paw_sphharm,      only : initylmr,setsym_ylm
  use m_pawtab,           only : pawtab_type
  use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_free,&
       &                         pawcprj_get, pawcprj_put, pawcprj_getdim, pawcprj_set_zero
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

     integer, allocatable :: kgindex(:)      ! kgind(nkpt)
     ! kgind(ikpt) = ikg

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
  public :: chern_number
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_orbmag'
!End of the abilint section

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
!! Copyright (C) 2004-2017 ABINIT group.
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initorbmag'
!End of the abilint section

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

  !JWZ: The following may need modification in the future
  !**** no spin-polarization doubling ; do not allow use of time reversal symmetry ****

  call timab(1002,2,tsec)
  call timab(1003,1,tsec)

  call listkk(rdum,gmet,dtorbmag%indkk_f2ibz,dtset%kptns,dtorbmag%fkptns,dtset%nkpt,&
       & dtorbmag%fnkpt,dtset%nsym,1,dtset%symafm,symrec,0,use_symrec=.True.)

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
!!****f* ABINIT/chern_number
!! NAME
!! chern_number
!!
!! FUNCTION
!! This routine computes the Chern number based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!! gmet(3,3)=metric in reciprocal space
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!! mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!! mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!! mpi_enreg=information about MPI parallelization
!! npwarr(nkpt)=number of planewaves in basis at this k point
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!! pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!! pwind_alloc = first dimension of pwind
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
!!
!! SOURCE

subroutine chern_number(atindx1,cg,cprj,dtset,dtorbmag,gmet,gprimd,kg,&
     &            mcg,mcprj,mpi_enreg,npwarr,pawang,pawrad,pawtab,pwind,pwind_alloc,&
     &            symrec,usecprj,usepaw,xred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chern_number'
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,pwind_alloc,usecprj,usepaw
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,bdx,bdxc,bfor,bsigma,ddkflag,epsabg,gdir,gdx,gdxc,gfor,gsigma
  integer :: icg,icgb,icgg,icprj,icprjb,icprjg
  integer :: ikg,ikpt,ikptb,ikptg,isppol,itrs,job
  integer :: mcg1_k,my_nspinor,nband_k,ncpgr,nn,n1,n2,n3,npw_k,npw_kb,npw_kg,shiftbd
  real(dp) :: deltab,deltag
  complex(dpc) :: IA,IB,t1A,t2A,t3A,t1B,t2B,t3B,t4B,tprodA,tprodB
  character(len=500) :: message
  !arrays
  integer,allocatable :: dimlmn(:),nattyp_dum(:),pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
  real(dp) :: cnum(2,3),dkb(3),dkg(3),dkbg(3),dtm_k(2)
  real(dp),allocatable :: cg1_k(:,:),kk_paw(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_all_indx(:,:,:,:,:,:),smat_inv(:,:,:),smat_kk(:,:,:)
  logical,allocatable :: has_smat_indx(:,:,:)
  type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kg(:,:)
  type(pawcprj_type),allocatable :: cprj_fkn(:,:),cprj_ikn(:,:)

  ! ***********************************************************************
  ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

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

  do adir = 1, 3

     IA = czero
     IB = czero

     do epsabg = 1, -1, -2

        if (epsabg .EQ. 1) then
           bdir = modulo(adir,3)+1
           gdir = modulo(adir+1,3)+1
        else
           bdir = modulo(adir+1,3)+1
           gdir = modulo(adir,3)+1
        end if

        ! loop over kpts, assuming for now kptopt 3, nsppol = 1, nspinor = 1
        ! and no parallelism, no symmorphic symmetry elements


        do ikpt = 1, dtorbmag%fnkpt

           icprj = dtorbmag%cprjindex(ikpt,isppol)

           npw_k = npwarr(ikpt)
           icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

           ikg = dtorbmag%fkgindex(ikpt)

           call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
                &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

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
              deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))

              ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
              icprjb = dtorbmag%cprjindex(ikptb,isppol)

              npw_kb = npwarr(ikptb)
              icgb = dtorbmag%cgindex(ikptb,dtset%nsppol)

              pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

              call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjb,&
                   &         ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                   &         my_nspinor,dtset%nsppol,0)

              if (.NOT. has_smat_indx(ikpt,bdx,0)) then

                 call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                      &           dtorbmag%lmn_size,dtset%mband,&
                      &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                 sflag_k=0
                 call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
                      &           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                      &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                 do nn = 1, nband_k
                    do n1 = 1, nband_k
                       smat_all_indx(1,nn,n1,ikpt,bdx,0) =  smat_kk(1,nn,n1)
                       smat_all_indx(2,nn,n1,ikpt,bdx,0) =  smat_kk(2,nn,n1)
                       smat_all_indx(1,n1,nn,ikptb,bdxc,0) =  smat_kk(1,nn,n1)
                       smat_all_indx(2,n1,nn,ikptb,bdxc,0) = -smat_kk(2,nn,n1)
                    end do
                 end do

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
                 deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))

                 ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                 icprjg = dtorbmag%cprjindex(ikptg,isppol)

                 npw_kg = npwarr(ikptg)
                 icgg = dtorbmag%cgindex(ikptg,dtset%nsppol)

                 pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                 call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjg,&
                      &           ikptg,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                      &           my_nspinor,dtset%nsppol,0)

                 if (.NOT. has_smat_indx(ikpt,gdx,0)) then

                    call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &             dtorbmag%lmn_size,dtset%mband,&
                         &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                    sflag_k=0
                    call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
                         &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                         &             pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                    do nn = 1, nband_k
                       do n1 = 1, nband_k
                          smat_all_indx(1,nn,n1,ikpt,gdx,0) =  smat_kk(1,nn,n1)
                          smat_all_indx(2,nn,n1,ikpt,gdx,0) =  smat_kk(2,nn,n1)
                          smat_all_indx(1,n1,nn,ikptg,gdxc,0) =  smat_kk(1,nn,n1)
                          smat_all_indx(2,n1,nn,ikptg,gdxc,0) = -smat_kk(2,nn,n1)
                       end do
                    end do

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
                         &             kg,dtorbmag%kgindex,mpi_enreg,npw_kb,pwind_bg,symrec)

                    sflag_k=0
                    call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
                         &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                         &             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,usepaw)

                    do nn = 1, nband_k
                       do n1 = 1, nband_k
                          smat_all_indx(1,nn,n1,ikpt,bdx,gdx) =  smat_kk(1,nn,n1)
                          smat_all_indx(2,nn,n1,ikpt,bdx,gdx) =  smat_kk(2,nn,n1)
                          smat_all_indx(1,n1,nn,ikpt,gdx,bdx) =  smat_kk(1,nn,n1)
                          smat_all_indx(2,n1,nn,ikpt,gdx,bdx) = -smat_kk(2,nn,n1)
                       end do
                    end do

                    has_smat_indx(ikpt,bdx,gdx) = .TRUE.
                    has_smat_indx(ikpt,gdx,bdx) = .TRUE.

                 end if

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

                             tprodB = t1B*t2B*t3B*t4B
                             IB = IB - epsabg*bsigma*gsigma*tprodB/(2.0*deltab*2.0*deltag)
                          end do ! end loop over n3

                          tprodA = t1A*t2A*t3A
                          IA = IA + epsabg*bsigma*gsigma*tprodA/(2.0*deltab*2.0*deltag)


                       end do ! end loop over n2
                    end do ! end loop over n1
                 end do ! end loop over nn

              end do ! end loop over gfor

           end do ! end loop over bfor

        end do ! end loop over fnkpt

     end do ! end loop over epsabg

     cnum(1,adir) = real(IA+IB)
     cnum(2,adir) = aimag(IA+IB)

  end do ! end loop over adir

  cnum(1,1:3) = MATMUL(gprimd,cnum(1,1:3))
  cnum(2,1:3) = MATMUL(gprimd,cnum(2,1:3))
  dtorbmag%chern(1,1:3) = -cnum(2,1:3)/(two_pi*dtorbmag%fnkpt)
  dtorbmag%chern(2,1:3) =  cnum(1,1:3)/(two_pi*dtorbmag%fnkpt)

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
!!****f* ABINIT/orbmag
!! NAME
!! orbmag
!!
!! FUNCTION
!! This routine computes the orbital magnetization based on input wavefunctions.
!! It is assumed that only completely filled bands are present.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine orbmag(atindx1,cg,cprj,dtset,dtorbmag,kg,&
     &            mcg,mcprj,mpi_enreg,nattyp,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,&
     &            pwind,pwind_alloc,rprimd,symrec,usecprj,vhartr,vpsp,vxc,xred,ylm,ylmgr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'orbmag'
!End of the abilint section

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
 integer :: iatom,icg,icgb,icgg,icprj,icprjb,icprjg,ider,idir
 integer :: ikg,ikgb,ikgg,ikpt,ikptb,ikptg,ilm,ilmn,ipw,isppol,istwf_k,itrs,itypat
 integer :: jlmn,job,jpw
 integer :: klmn,mcg1_k,my_cpopt,my_nspinor,nband_k,ncpgr,ndat,nkpg,nn,n1,n2
 integer :: ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,npw_k,npw_kb,npw_kg
 integer :: prtvol,shiftbd,sij_opt,tim_getghc,type_calc,type_calc_123
 real(dp) :: deltab,deltag,dkg2,dotr,doti,htpisq,keg,lambda,ucvol
 complex(dpc) :: cdij,cgdijcb,cpb,cpg
 complex(dpc) :: IIA,IIA1,IIA2,IIA3,IIIA,IIIA1,IIIA2,IIIA3,tprodIIA,tprodIIIA
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk,gs_hamk123
 !arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:),kg_kb(:,:),kg_kg(:,:)
 integer,allocatable :: pwind_kb(:),pwind_kg(:),pwind_bg(:),sflag_k(:)
 real(dp) :: dkb(3),dkg(3),dkbg(3),dtm_k(2),gmet(3,3),gprimd(3,3)
 real(dp) :: kpoint(3),kpointb(3),kpointg(3)
 real(dp) :: orbmagvec(2,3),rhodum(1),rmet(3,3)
 real(dp),allocatable :: bra(:,:),cg1_k(:,:),cgrvtrial(:,:),cwavef(:,:),ffnl(:,:,:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: hmat(:,:,:,:,:,:),kinpw(:),kk_paw(:,:,:),kpg_k_dummy(:,:)
 real(dp),allocatable :: my_nucdipmom(:,:),ph3d(:,:,:),pwnsfac_k(:,:),smat_all(:,:,:,:,:,:),smat_inv(:,:,:)
 real(dp),allocatable :: smat_kk(:,:,:),vlocal(:,:,:,:),vtrial(:,:),ylm_k(:,:)
 complex(dpc),allocatable :: nucdipmom_k(:)
 logical,allocatable :: has_hmat(:,:,:),has_smat(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kb_k(:,:,:,:)
 type(pawcprj_type),allocatable :: cprj_kg(:,:),cwaveprj(:,:)

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

   ABI_DATATYPE_ALLOCATE(cprj_kb_k,(dtorbmag%fnkpt,6,dtset%natom,dtorbmag%nspinor*dtset%mband))
   do ikpt=1,dtorbmag%fnkpt
      do bdx=1, 6
         call pawcprj_alloc(cprj_kb_k(ikpt,bdx,:,:),ncpgr,dimlmn)
      end do
   end do
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

 do adir = 1, 3

    IIA  = czero
    IIIA = czero

    do epsabg = 1, -1, -2

       if (epsabg .EQ. 1) then
          bdir = modulo(adir,3)+1
          gdir = modulo(adir+1,3)+1
       else
          bdir = modulo(adir+1,3)+1
          gdir = modulo(adir,3)+1
       end if

       ! loop over kpts, assuming for now kptopt 3, nsppol = 1, nspinor = 1
       ! and no parallelism, no symmorphic symmetry elements
       do ikpt = 1, dtorbmag%fnkpt

          kpoint(:)=dtorbmag%fkptns(:,ikpt)
          icprj = dtorbmag%cprjindex(ikpt,isppol)

          npw_k = npwarr(ikpt)
          icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

          ikg = dtorbmag%fkgindex(ikpt)
          ABI_ALLOCATE(kg_k,(3,npw_k))
          kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

          call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
               &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

          ! Set up remainder of normal Hamiltonian at k if necessary

          if (.NOT. has_hmat(ikpt,0,0) ) then
             ! kpoint(:)=dtorbmag%fkptns(:,ikpt)

             ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
             if (psps%useylm==1) then
                do ilm=1,psps%mpsang*psps%mpsang
                   ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
                end do
             end if

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
             ! pretty sure do not need k+g vectors, use dummy kpg_k
             ! may eventually need them for nucdipmom_k hamiltonian if
             ! generalize to k,k'
             nkpg = 0
             ABI_ALLOCATE(kpg_k_dummy,(npw_k,nkpg))

             !      Compute nonlocal form factors ffnl at all (k+G):
             ider=0;idir=0;dimffnl=1
             ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
             call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                  &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k_dummy,kpoint,psps%lmnmax,&
                  &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                  &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                  &         psps%usepaw,psps%useylm,ylm_k,ylmgr)

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

             ABI_DEALLOCATE(ylm_k)
             ABI_DEALLOCATE(kpg_k_dummy)
             ABI_DEALLOCATE(ffnl)
             if(any(abs(gs_hamk%nucdipmom)>0.0)) then
                if(allocated(nucdipmom_k)) then
                   ABI_DEALLOCATE(nucdipmom_k)
                end if
             end if
             ABI_DEALLOCATE(ph3d)

          end if ! end check on has_hmat

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
             deltab = sqrt(DOT_PRODUCT(dkb,MATMUL(gmet,dkb)))

             ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
             kpointb(:)=dtorbmag%fkptns(:,ikptb)

             icprjb = dtorbmag%cprjindex(ikptb,isppol)

             ikgb = dtorbmag%fkgindex(ikptb)
             npw_kb = npwarr(ikptb)
             ABI_ALLOCATE(kg_kb,(3,npw_kb))
             kg_kb(:,1:npw_kb)=kg(:,ikgb+1:ikgb+npw_kb)

             icgb = dtorbmag%cgindex(ikptb,dtset%nsppol)

             pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

             call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjb,&
                  &         ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                  &         my_nspinor,dtset%nsppol,0)

             if (.NOT. has_smat(ikpt,bdx,0)) then

                call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                     &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                sflag_k=0
                call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgb,itrs,job,nband_k,&
                     &           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                     &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                do nn = 1, nband_k
                   do n1 = 1, nband_k
                      smat_all(1,nn,n1,ikpt,bdx,0) =  smat_kk(1,nn,n1)
                      smat_all(2,nn,n1,ikpt,bdx,0) =  smat_kk(2,nn,n1)
                      smat_all(1,n1,nn,ikptb,bdxc,0) =  smat_kk(1,nn,n1)
                      smat_all(2,n1,nn,ikptb,bdxc,0) = -smat_kk(2,nn,n1)
                   end do
                end do

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
                deltag = sqrt(DOT_PRODUCT(dkg,MATMUL(gmet,dkg)))

                ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                kpointg(:)=dtorbmag%fkptns(:,ikptg)

                icprjg = dtorbmag%cprjindex(ikptg,isppol)

                npw_kg = npwarr(ikptg)
                ABI_ALLOCATE(kg_kg,(3,npw_kg))
                ikgg = dtorbmag%fkgindex(ikptg)
                kg_kg(:,1:npw_kg)=kg(:,ikgg+1:ikgg+npw_kg)

                icgg = dtorbmag%cgindex(ikptg,dtset%nsppol)

                pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                call pawcprj_get(atindx1,cprj_kg,cprj,dtset%natom,1,icprjg,&
                     &           ikptg,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                     &           my_nspinor,dtset%nsppol,0)

                if (.NOT. has_smat(ikpt,gdx,0)) then

                   call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icg,icgg,itrs,job,nband_k,&
                        &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                        &             pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   do nn = 1, nband_k
                      do n1 = 1, nband_k
                         smat_all(1,nn,n1,ikpt,gdx,0) =  smat_kk(1,nn,n1)
                         smat_all(2,nn,n1,ikpt,gdx,0) =  smat_kk(2,nn,n1)
                         smat_all(1,n1,nn,ikptg,gdxc,0) =  smat_kk(1,nn,n1)
                         smat_all(2,n1,nn,ikptg,gdxc,0) = -smat_kk(2,nn,n1)
                      end do
                   end do

                   has_smat(ikpt,gdx,0) = .TRUE.
                   has_smat(ikptg,gdxc,0) = .TRUE.

                end if

                dkbg = dkg - dkb

                if (.NOT. has_smat(ikpt,bdx,gdx)) then

                   call overlap_k1k2_paw(cprj_kb,cprj_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,dtorbmag%lmn_size,dtset%mband,&
                        &             dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

                   call mkpwind_k(dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                        &             dtorbmag%indkk_f2ibz,ikptb,ikptg,&
                        &             kg,dtorbmag%kgindex,mpi_enreg,npw_kb,pwind_bg,symrec)

                   sflag_k=0
                   call smatrix(cg,cg,cg1_k,ddkflag,dtm_k,icgb,icgg,itrs,job,nband_k,&
                        &             mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_kb,npw_kg,my_nspinor,&
                        &             pwind_bg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                   do nn = 1, nband_k
                      do n1 = 1, nband_k
                         smat_all(1,nn,n1,ikpt,bdx,gdx) =  smat_kk(1,nn,n1)
                         smat_all(2,nn,n1,ikpt,bdx,gdx) =  smat_kk(2,nn,n1)
                         smat_all(1,n1,nn,ikpt,gdx,bdx) =  smat_kk(1,nn,n1)
                         smat_all(2,n1,nn,ikpt,gdx,bdx) = -smat_kk(2,nn,n1)
                      end do
                   end do

                   has_smat(ikpt,bdx,gdx) = .TRUE.
                   has_smat(ikpt,gdx,bdx) = .TRUE.

                end if

                if (.NOT. has_hmat(ikpt,gdx,bdx)) then

                   call mkpwind_k(-dkbg,dtset,dtorbmag%fnkpt,dtorbmag%fkptns,gmet,&
                        &             dtorbmag%indkk_f2ibz,ikptg,ikptb,&
                        &             kg,dtorbmag%kgindex,mpi_enreg,npw_kg,pwind_bg,symrec)

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
                               cpg=cmplx(cprj_kb_k(ikptg,gdxc,iatom,n1)%cp(1,ilmn),&
                                    & cprj_kb_k(ikptg,gdxc,iatom,n1)%cp(2,ilmn))
                               do jlmn = 1, pawtab(itypat)%lmn_size
                                  cpb=cmplx(cprj_kb_k(ikptb,bdxc,iatom,nn)%cp(1,jlmn),&
                                       & cprj_kb_k(ikptb,bdxc,iatom,nn)%cp(2,jlmn))
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

                do nn = 1, nband_k
                   do n1 = 1, nband_k

                      IIA1  = cmplx(smat_all(1,nn,n1,ikpt,gdx,0),smat_all(2,nn,n1,ikpt,gdx,0))

                      IIIA1 = cmplx(smat_all(1,nn,n1,ikpt,bdx,0),smat_all(2,nn,n1,ikpt,bdx,0))

                      do n2 = 1, nband_k

                         IIA2 = cmplx(hmat(1,n1,n2,ikpt,gdx,bdx),hmat(2,n1,n2,ikpt,gdx,bdx))
                         IIA3 = cmplx(smat_all(1,n2,nn,ikptb,bdxc,0),smat_all(2,n2,nn,ikptb,bdxc,0))

                         IIIA2 = cmplx(smat_all(1,n1,n2,ikpt,bdx,gdx),smat_all(2,n1,n2,ikpt,bdx,gdx))
                         IIIA3 = conjg(cmplx(smat_all(1,nn,n2,ikpt,gdx,0),smat_all(2,nn,n2,ikpt,gdx,0)))

                         tprodIIA  = IIA1*IIA2*IIA3
                         IIA = IIA + epsabg*bsigma*gsigma*tprodIIA/(2.0*deltab*2.0*deltag)

                         tprodIIIA = hmat(1,nn,nn,ikpt,0,0)*IIIA1*IIIA2*IIIA3
                         IIIA = IIIA - epsabg*bsigma*gsigma*tprodIIIA/(2.0*deltab*2.0*deltag)

                      end do ! end n2
                   end do ! end n1
                end do ! end nn

                ABI_DEALLOCATE(kg_kg)

             end do ! end gfor

             ABI_DEALLOCATE(kg_kb)

          end do ! end bfor

          ABI_DEALLOCATE(kg_k)

       end do ! end loop over fnkpt

    end do ! end loop over epsabg

    orbmagvec(1,adir) = real(IIA+IIIA)
    orbmagvec(2,adir) = aimag(IIA+IIIA)
    ! orbmagvec(1,adir) = real(IIIA)
    ! orbmagvec(2,adir) = aimag(IIIA)
 end do ! end loop over adir

 orbmagvec(1,1:3) = MATMUL(gprimd,orbmagvec(1,1:3))
 orbmagvec(2,1:3) = MATMUL(gprimd,orbmagvec(2,1:3))
 dtorbmag%orbmagvec(1,1:3) =  orbmagvec(2,1:3)/(two*dtorbmag%fnkpt)
 dtorbmag%orbmagvec(2,1:3) = -orbmagvec(1,1:3)/(two*dtorbmag%fnkpt)

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
         call pawcprj_free(cprj_kb_k(ikpt,bdx,:,:))
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

 ABI_DEALLOCATE(my_nucdipmom)

 ABI_DEALLOCATE(vlocal)
 call destroy_hamiltonian(gs_hamk)
 call destroy_hamiltonian(gs_hamk123)

end subroutine orbmag
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
!! Copyright (C) 2003-2017 ABINIT  group
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
!!
!! SOURCE

subroutine ctocprjb(atindx1,cg,cprj_kb_k,dtorbmag,dtset,gmet,gprimd,&
     & istwf_k,kg,mcg,mpi_enreg,nattyp,ncpgr,npwarr,pawtab,psps,rmet,rprimd,ucvol,xred)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ctocprjb'
!End of the abilint section

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
  type(pawcprj_type),intent(inout) :: cprj_kb_k(dtorbmag%fnkpt,6,dtset%natom,dtorbmag%nspinor*dtset%mband)

  !Locals------------------------------------
  !scalars
  integer :: bdir,bdx,bfor,bsigma,choice,cpopt,dimffnl,dimph1d,ia,iband,icg,ider,idir,ikg,ikpt
  integer :: n1,n2,n3,nband_k,nkpg,npw_k,optder
  real(dp) :: arg

  !arrays
  integer :: nband_dum(1),npwarr_dum(1)
  integer,allocatable :: dimlmn(:),kg_k(:,:)
  real(dp) :: dkb(3),kpoint(3),kpointb(3),kptns(3,1)
  real(dp),allocatable :: cwavef(:,:),ffnl(:,:,:,:),kpg_k_dummy(:,:)
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

     dimffnl=1
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))

     icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

     ikg = dtorbmag%fkgindex(ikpt)
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     nkpg = 0
     ABI_ALLOCATE(kpg_k_dummy,(npw_k,nkpg))

     ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))

     ! data for initylmg call below
     optder=0
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

           kptns(:,1) = kpointb(:)
           call initylmg(gprimd,kg_k,kptns,1,mpi_enreg,psps%mpsang,npw_k,&
                & nband_dum,1,npwarr_dum,dtset%nsppol,optder,rprimd,ylm_k,ylm_k_gr)

           !      Compute nonlocal form factors ffnl at all (k+b+G_k):
           ider=0;idir=0
           call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
                & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k_dummy,kpointb,psps%lmnmax,&
                & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                & psps%usepaw,psps%useylm,ylm_k,ylm_k_gr)

           choice=1;cpopt=0;idir=0
           do iband = 1, nband_k
              cwavef(1,1:npw_k) = cg(1,icg+(iband-1)*npw_k+1:icg+iband*npw_k)
              cwavef(2,1:npw_k) = cg(2,icg+(iband-1)*npw_k+1:icg+iband*npw_k)

              call getcprj(choice,cpopt,cwavef,cwaveprj,ffnl,&
                   & idir,psps%indlmn,istwf_k,kg_k,kpg_k_dummy,kpointb,psps%lmnmax,&
                   & dtset%mgfft,mpi_enreg,&
                   & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                   & phkxred,ph1d,ph3d,ucvol,psps%useylm)

              call pawcprj_put(atindx1,cwaveprj,cprj_kb_k(ikpt,bdx,:,:),dtset%natom,&
                   & iband,0,ikpt,0,1,nband_k,1,dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0)

           end do ! end loop over bands

        end do ! end loop over bfor
     end do ! end loop over bdir

     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(kpg_k_dummy)
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


end module m_orbmag
!!***
