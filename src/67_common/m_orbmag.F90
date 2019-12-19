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

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_dtset

  use defs_datatypes,     only : pseudopotential_type
  use defs_abitypes,      only : MPI_type
  use m_berrytk,          only : smatrix
  use m_cgprj,            only : getcprj
  use m_cgtools,          only : overlap_g
  use m_fft,              only : fftpac,fourwf
  use m_fftcore,          only : kpgsph,sphereboundary
  use m_fourier_interpol, only : transgrid
  use m_geometry,         only : metric
  use m_getghc,           only : getghc
  use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type
  use m_initylmg,         only : initylmg
  use m_kg,               only : getph,mkkin,mkkpg,mkpwind_k,ph1d3d
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
  public :: orbmag

  private :: covar_test
  private :: covar_cprj
  private :: rho_norm_check
  private :: duqdu
  private :: mpicomm_helper
  private :: duqhqdu
  private :: cpg_dij_cpb
  private :: udsqdu
  private :: udsdsu
  private :: make_eeig
  private :: make_onsite_l
  private :: make_onsite_l_k
  private :: make_onsite_bm
  private :: make_S1trace
  private :: make_rhorij1
  private :: output_orbmag
  
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
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE


subroutine rho_norm_check(atindx1,cg,cprj,dtorbmag,dtset,mpi_enreg,mcg,mcprj,&
     & npwarr,pawtab,usecprj,usepaw)

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
              cpb=cmplx(cprj_k(iatom,iband)%cp(1,ilmn),cprj_k(iatom,iband)%cp(2,ilmn),KIND=dpc)
              cpk=cmplx(cprj_k(iatom,iband)%cp(1,jlmn),cprj_k(iatom,iband)%cp(2,jlmn),KIND=dpc)
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
!!****f* ABINIT/covar_cprj
!! NAME
!! covar_cprj
!!
!! FUNCTION
!! Generate cprj multiplied by S^{-1}, similarly to the wavefunctions in smatrix
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE


subroutine covar_cprj(cprj_kb,cprj_kb_covar,dtset,nband_k,pawtab,smat_inv)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nband_k
  type(dataset_type),intent(in) :: dtset
  type(pawcprj_type),intent(in) ::  cprj_kb(dtset%natom,nband_k)
  type(pawcprj_type),intent(inout) ::  cprj_kb_covar(dtset%natom,nband_k)

  !arrays
  real(dp),intent(in) :: smat_inv(2,nband_k,nband_k)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iatom,iband,ilmn,itypat,jband
  complex(dpc) :: cpb,smi

  !arrays

  !----------------------------------------------------

  ! make covar cprj same as covar |unk>
  do iatom=1,dtset%natom
     itypat = dtset%typat(iatom)
     do ilmn=1,pawtab(itypat)%lmn_size
        do jband = 1, nband_k
           cprj_kb_covar(iatom,jband)%cp(:,ilmn) = zero
           do iband = 1, nband_k
              cpb=cmplx(cprj_kb(iatom,iband)%cp(1,ilmn),cprj_kb(iatom,iband)%cp(2,ilmn),KIND=dpc)
              smi=cmplx(smat_inv(1,iband,jband),smat_inv(2,iband,jband),KIND=dpc)
              cprj_kb_covar(iatom,jband)%cp(1,ilmn) = cprj_kb_covar(iatom,jband)%cp(1,ilmn) + &
                   & real(cpb*smi)
              cprj_kb_covar(iatom,jband)%cp(2,ilmn) = cprj_kb_covar(iatom,jband)%cp(2,ilmn) + &
                   & aimag(cpb*smi)
           end do ! end loop over iband
        end do ! end loop over jband
     end do ! end loop over ilmn
  end do ! end loop over iatom

end subroutine covar_cprj
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/covar_test
!! NAME
!! covar_test
!!
!! FUNCTION
!! Routine to play with and check covariant derivative Q|\partial_k u_{nk}\rangle$
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE


subroutine covar_test(atindx1,cg,cprj,dtorbmag,dtset,gmet,gprimd,mcg,mcprj,mpi_enreg,&
      & nband_k,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,xred)

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
  integer,intent(in) :: atindx1(dtset%natom)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,bfor,bsigma,iatom,iband,icg,icgb,icprj,icprjb
  integer :: ikg,ikpt,ikptb,ilmn,ipw,isppol,itrs,itypat
  integer :: jband,mcg1_k,my_nspinor,ncpgr,npw_k,npw_kb
  integer :: shiftbd,smatrix_ddkflag,smatrix_job
  real(dp) :: doti,dotr
  complex(dpc) :: cpb,smi

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:),pwind_kb(:),sflag_k(:)
  real(dp) :: dkb(3),dtm_k(2)
  real(dp),allocatable :: cg_k(:,:),cg1_k(:,:),kk_paw(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_inv(:,:,:),smat_kk(:,:,:)
  complex(dpc),allocatable :: smat_invc(:,:),smat_kkc(:,:),smat(:,:)
  type(pawcprj_type),allocatable :: cprj_k(:,:),cprj_kb(:,:),cprj_kb_covar(:,:)

  !----------------------------------------------------

  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kb_covar,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kb_covar,ncpgr,dimlmn)

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = one; pwnsfac_k(2,:) = zero
  pwnsfac_k(3,:) = one; pwnsfac_k(4,:) = zero

  mcg1_k = dtset%mpw*dtset%nsppol*my_nspinor*nband_k
  ABI_ALLOCATE(cg1_k,(2,mcg1_k))
  ABI_ALLOCATE(cg_k,(2,mcg1_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_invc,(nband_k,nband_k))
  ABI_ALLOCATE(smat_kkc,(nband_k,nband_k))
  ABI_ALLOCATE(smat,(nband_k,nband_k))

  smatrix_ddkflag = 1
  itrs = 0
  smatrix_job = 1
  shiftbd = 1

  do ikpt = 1, dtorbmag%fnkpt

     icprj = dtorbmag%cprjindex(ikpt,isppol)

     npw_k = npwarr(ikpt)
     ikg = dtorbmag%fkgindex(ikpt)
     icg = dtorbmag%cgindex(ikpt,dtset%nsppol)

     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

     do bdir = 1, 3
        do bfor = 1, 2
           bsigma = -2*bfor+3
           ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
           npw_kb = npwarr(ikptb)
           icprjb = dtorbmag%cprjindex(ikptb,isppol)
           icgb = dtorbmag%cgindex(ikptb,dtset%nsppol)
           dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
           pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)
           call pawcprj_get(atindx1,cprj_kb,cprj,dtset%natom,1,icprjb,&
                &         ikptb,0,isppol,dtset%mband,dtset%mkmem,dtset%natom,nband_k,nband_k,&
                &         my_nspinor,dtset%nsppol,0)
           call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                &           dtorbmag%lmn_size,dtset%mband,&
                &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

           sflag_k=0
           cg1_k(:,:) = zero
           call smatrix(cg,cg,cg1_k,smatrix_ddkflag,dtm_k,icg,icgb,itrs,smatrix_job,nband_k,&
                &           mcg,mcg,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

           cg_k(1:2,1:mcg1_k) = cg(1:2,icg+1:icg+mcg1_k)

           call covar_cprj(cprj_kb,cprj_kb_covar,dtset,nband_k,pawtab,smat_inv)

           call overlap_k1k2_paw(cprj_k,cprj_kb_covar,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                &           dtorbmag%lmn_size,dtset%mband,&
                &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)

           do iband = 1, nband_k
              do jband = 1, nband_k
                 dotr=DOT_PRODUCT(cg_k(1,(iband-1)*npw_k+1:iband*npw_k),cg1_k(1,(jband-1)*npw_k+1:jband*npw_k)) + &
                      & DOT_PRODUCT(cg_k(2,(iband-1)*npw_k+1:iband*npw_k),cg1_k(2,(jband-1)*npw_k+1:jband*npw_k))
                 doti=DOT_PRODUCT(cg_k(1,(iband-1)*npw_k+1:iband*npw_k),cg1_k(2,(jband-1)*npw_k+1:jband*npw_k)) - &
                      & DOT_PRODUCT(cg_k(2,(iband-1)*npw_k+1:iband*npw_k),cg1_k(1,(jband-1)*npw_k+1:jband*npw_k))
                 write(std_out,'(a,4i4,2es16.8)')'JWZ debug ikpt ikptb iband jband ',ikpt,ikptb,iband,jband,&
                      & dotr+kk_paw(1,iband,jband),doti+kk_paw(2,iband,jband)
              end do
           end do
           
        end do ! end loop over bfor
     end do ! end loop over bdir

  end do ! end loop over ikpt

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cprj_kb)
  ABI_DATATYPE_DEALLOCATE(cprj_kb)
  call pawcprj_free(cprj_kb_covar)
  ABI_DATATYPE_DEALLOCATE(cprj_kb_covar)

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(smat_invc)
  ABI_DEALLOCATE(smat_kkc)
  ABI_DEALLOCATE(smat)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(pwind_kb)
  ABI_DEALLOCATE(cg1_k)
  ABI_DEALLOCATE(cg_k)
  ABI_DEALLOCATE(pwnsfac_k)

end subroutine covar_test
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/duqdu
!! NAME
!! duqdu
!!
!! FUNCTION
!! Return i*epsabg\sum_n E_nk <\partial_b u_kn|Q|\partial_g u_kn> where
!! Q projects onto the conduction space.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!! energies is an optional argument. When not used, E_nk will be set to 1 and this
!! is appropriate for a Chern number call. When energies are set to the Bloch eigenvalues,
!! this is appropriate for a magnetization call.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


subroutine duqdu(atindx1,cg,cprj,dtorbmag,dtset,duqduchern,duqdumag,energies,&
     & gmet,gprimd,mcg,mcprj,mpi_enreg,nband_k,npwarr,pawang,pawrad,pawtab,&
     & psps,pwind,pwind_alloc,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,nband_k,pwind_alloc
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom)
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),xred(3,dtset%natom)
  real(dp),intent(in) :: energies(nband_k,dtset%nkpt)
  real(dp), intent(out) :: duqduchern(2,3),duqdumag(2,3)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,bfor,bsigma,countb,countg,countjb,countjg,countk,dest
  integer :: epsabg,gdir,gfor,gsigma,iatom,iband
  integer :: icg,icgb,icgg,icprji,icprjbi,icprjgi,ierr
  integer :: ikg,ikpt,ikpt_loc,ikpti,ikptb,ikptbi,ikptg,ikptgi,ish1,ish2,isppol,itrs,itypat
  integer :: jcgb,jcgg,jcprjbi,jcprjgi,jkpt,jkptb,jkptbi,jkptg,jkptgi,jsppol
  integer :: me,mcg1_k,my_nspinor,n2dim,ncpgr,npw_k,npw_kb,npw_kg,nproc,ntotcp
  integer :: shiftbd,smatrix_ddkflag,smatrix_job,sourceb,sourceg,spaceComm
  integer :: tagb,tagg
  real(dp) :: deltab,deltag,doti,dotr,ENK
  complex(dpc) :: cprefac,duqduchern_term,duqdumag_term

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:),pwind_kb(:),pwind_kg(:),sflag_k(:)
  real(dp) :: dkb(3),dkbg(3),dkg(3),dtm_k(2)
  real(dp),allocatable :: buffer(:,:)
  real(dp),allocatable :: cg_k(:,:),cg1_kb(:,:),cg1_kg(:,:),cgqb(:,:),cgqg(:,:)
  real(dp),allocatable :: kk_paw(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_inv(:,:,:),smat_kk(:,:,:)
  type(pawcprj_type),allocatable :: cprj_buf(:,:),cprj_k(:,:),cprj_kb(:,:),cprj1_kb(:,:)
  type(pawcprj_type),allocatable :: cprj_kg(:,:),cprj1_kg(:,:)

  !----------------------------------------------------

  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me=mpi_enreg%me_kpt

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj1_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj1_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj1_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj1_kg,ncpgr,dimlmn)
  n2dim = dtorbmag%nspinor*nband_k
  ntotcp = n2dim*SUM(dimlmn(:))
  if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(dtset%natom,n2dim))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
  end if

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))
  ABI_ALLOCATE(pwind_kg,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = one; pwnsfac_k(2,:) = zero
  pwnsfac_k(3,:) = one; pwnsfac_k(4,:) = zero

  mcg1_k = dtset%mpw*dtset%nsppol*my_nspinor*nband_k
  ABI_ALLOCATE(cg_k,(2,mcg1_k))
  ABI_ALLOCATE(cg1_kb,(2,mcg1_k))
  ABI_ALLOCATE(cg1_kg,(2,mcg1_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

  smatrix_ddkflag = 1
  itrs = 0
  smatrix_job = 1
  shiftbd = 1

  duqduchern(:,:) = zero
  duqdumag(:,:) = zero

  do ikpt_loc = 1,dtorbmag%fmkmem_max

     ikpt=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
     ! if this k and spin are for me do it
     ! if (ikpt1 > 0 .and. isppol > 0) then
     if (ikpt > 0) then

        ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
        icprji = dtorbmag%cprjindex(ikpti,isppol)
        npw_k = npwarr(ikpti)
        icg = dtorbmag%cgindex(ikpti,dtset%nsppol)
        ikg = dtorbmag%fkgindex(ikpt)
        countk = npw_k*my_nspinor*nband_k
        cg_k(1:2,1:countk) = cg(1:2,icg+1:icg+countk)
                    
        call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
             &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
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
              bsigma = 3-2*bfor
              dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
              deltab = sqrt(DOT_PRODUCT(dkb,dkb))

              if (ikpt > 0) then
                 ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
                 ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)
                 npw_kb = npwarr(ikptbi)
                 pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)
              end if

              if (ikpt > 0 .AND. isppol > 0) then
                 countb = npw_kb*my_nspinor*nband_k
                 if(allocated(cgqb)) then
                    ABI_DEALLOCATE(cgqb)
                 endif
                 ABI_ALLOCATE(cgqb,(2,countb))
                 call mpicomm_helper(atindx1,bdir,bfor,cg,cgqb,cprj,cprj_kb,dimlmn,dtorbmag,dtset,&
                      & ikpt,ikpt_loc,ikptbi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
                      & nproc,npw_kb,npwarr,spaceComm)
              end if
                
              if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps
                 
                 ! get covariant |u_{n,k+b}> and associated cprj
                 call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                      &           dtorbmag%lmn_size,dtset%mband,&
                      &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                 sflag_k=0
                 cg1_kb(:,:) = zero
                 ! cg1_kb will hold |\tilde{u}_{n,k+b}>
                 call smatrix(cg_k,cgqb,cg1_kb,smatrix_ddkflag,dtm_k,0,0,itrs,smatrix_job,nband_k,&
                      &           mcg1_k,mcg1_k,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                      &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)
                 ! cprj1_kb will hold cprj for cg1_kb
                 call covar_cprj(cprj_kb,cprj1_kb,dtset,nband_k,pawtab,smat_inv)

                 if(allocated(cgqb)) then
                    ABI_DEALLOCATE(cgqb)
                 end if

              end if

              do gfor = 1, 2
                 gsigma=3-2*gfor
                 dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
                 deltag = sqrt(DOT_PRODUCT(dkg,dkg))

                 cprefac = j_dpc*epsabg*bsigma*gsigma/(two*deltab*two*deltag)

                 if (ikpt > 0) then
                    ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                    ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)
                    npw_kg = npwarr(ikptgi)
                    pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                 end if
                 

                 if (ikpt > 0 .AND. isppol > 0) then
                    countg = npw_kg*my_nspinor*nband_k
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    endif
                    ABI_ALLOCATE(cgqg,(2,countg))
                    call mpicomm_helper(atindx1,gdir,gfor,cg,cgqg,cprj,cprj_kg,dimlmn,dtorbmag,dtset,&
                         & ikpt,ikpt_loc,ikptgi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
                         & nproc,npw_kg,npwarr,spaceComm)
                 end if
                    
                 if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps

                    ! get covariant |u_{n,k+g}> and associated cprj
                    call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &           dtorbmag%lmn_size,dtset%mband,&
                         &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                    sflag_k=0
                    cg1_kg(:,:) = zero
                    ! cg1_kg will hold |\tilde{u}_{n,k+g}>
                    call smatrix(cg_k,cgqg,cg1_kg,smatrix_ddkflag,dtm_k,0,0,itrs,smatrix_job,nband_k,&
                         &           mcg1_k,mcg1_k,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                         &           pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)
                    ! cprj1_kg will hold cprj for cg1_kg
                    call covar_cprj(cprj_kg,cprj1_kg,dtset,nband_k,pawtab,smat_inv)

                    dkbg = dkg - dkb
                    ! overlap of covariant cprj at kb and kg
                    call overlap_k1k2_paw(cprj1_kb,cprj1_kg,dkbg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &           dtorbmag%lmn_size,dtset%mband,&
                         &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                    do iband = 1, nband_k
                          
                       ish1 = (iband-1)*npw_k+1
                       ish2 = iband*npw_k
                       ENK = energies(iband,ikpt)
                       
                       dotr=DOT_PRODUCT(cg1_kb(1,ish1:ish2),cg1_kg(1,ish1:ish2)) + &
                            & DOT_PRODUCT(cg1_kb(2,ish1:ish2),cg1_kg(2,ish1:ish2))
                       doti=DOT_PRODUCT(cg1_kb(1,ish1:ish2),cg1_kg(2,ish1:ish2)) - &
                            & DOT_PRODUCT(cg1_kb(2,ish1:ish2),cg1_kg(1,ish1:ish2))
                       
                       ! accumulate i*epsabg*ENK*\sum_occ [<d_bdir u|Q|d_gdir u>]
                       duqduchern_term = cprefac*cmplx((dotr+kk_paw(1,iband,iband)),(doti+kk_paw(2,iband,iband)))
                       duqdumag_term = duqduchern_term*ENK
                       
                       duqduchern(1,adir) = duqduchern(1,adir) + real(duqduchern_term)
                       duqduchern(2,adir) = duqduchern(2,adir) + aimag(duqduchern_term)
                       duqdumag(1,adir) = duqdumag(1,adir) + real(duqdumag_term)
                       duqdumag(2,adir) = duqdumag(2,adir) + aimag(duqdumag_term)
                       
                    end do ! end loop over iband
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    end if
                 end if ! end check on ikpt > 0
                 
              end do ! end loop over gfor
           end do ! end loop over bfor
        end do ! end loop over epsabg
     end do ! end loop over adir
  end do ! end loop over ikpt_loc

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(duqduchern,spaceComm,ierr)
     call xmpi_sum(duqdumag,spaceComm,ierr)
  end if
  
  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cprj_kb)
  ABI_DATATYPE_DEALLOCATE(cprj_kb)
  call pawcprj_free(cprj1_kb)
  ABI_DATATYPE_DEALLOCATE(cprj1_kb)
  call pawcprj_free(cprj_kg)
  ABI_DATATYPE_DEALLOCATE(cprj_kg)
  call pawcprj_free(cprj1_kg)
  ABI_DATATYPE_DEALLOCATE(cprj1_kg)
  if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
  end if

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(pwind_kb)
  ABI_DEALLOCATE(pwind_kg)
  ABI_DEALLOCATE(cg1_kb)
  ABI_DEALLOCATE(cg1_kg)
  ABI_DEALLOCATE(cg_k)
  ABI_DEALLOCATE(pwnsfac_k)

end subroutine duqdu
!!***n

!{\src2tex{textfont=tt}}
!!****f* ABINIT/mpicomm_helper
!! NAME
!! mpicomm_helper
!!
!! FUNCTION
!! get wavefunction and cprj in mpi communication loop
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE

subroutine mpicomm_helper(atindx1,bdir,bfor,cg,cgqb,cprj,cprj_kb,dimlmn,dtorbmag,dtset,&
     & ikpt,ikpt_loc,ikptbi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
     & nproc,npw_kb,npwarr,spaceComm)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: bdir,bfor,ikpt,ikpt_loc,ikptbi,isppol,mcg,mcprj,me,my_nspinor
  integer,intent(in) :: nband_k,nproc,npw_kb,spaceComm
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),dimlmn(dtset%natom),npwarr(dtset%nkpt)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(out) :: cgqb(2,npw_kb*my_nspinor*nband_k)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawcprj_type),intent(inout) ::  cprj_kb(dtset%natom,nband_k)

  !Local variables -------------------------
  !scalars
  integer :: countb,countjb,dest,icgb,icprjbi,ierr
  integer :: jcgb,jcprjbi,jkpt,jkptb,jkptbi,jsppol,n2dim,ncpgr,sourceb,tagb

  !arrays
  real(dp),allocatable :: buffer(:,:)
  type(pawcprj_type),allocatable :: cprj_buf(:,:)

  !----------------------------------------------------

  n2dim = dtorbmag%nspinor*nband_k
  ncpgr = cprj(1,1)%ncpgr
  if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(dtset%natom,n2dim))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
  end if

  if (ikpt > 0 .and. isppol > 0) then ! I currently have a true kpt to use
     countb = npw_kb*my_nspinor*nband_k
     cgqb = zero
     sourceb = me
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikptbi,1,nband_k,isppol,me)) then
        ! I need the data from someone else
        sourceb = mpi_enreg%proc_distrb(ikptbi,1,isppol)
     end if
  else
     sourceb = -1 ! I do not have a kpt to use
  end if

  do dest=0,nproc-1
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
           ! receive cgqb and cprj_kb
           tagb = ikptbi + (isppol - 1)*dtset%nkpt
           call xmpi_recv(cgqb,sourceb,tagb,spaceComm,ierr)
           call pawcprj_mpi_recv(dtset%natom,n2dim,dimlmn,ncpgr,cprj_kb,sourceb,spaceComm,ierr)
        end if
     else if (dest.NE.me) then
        ! jkpt is the kpt which is being treated by dest
        ! jsppol is his isppol
        jkpt = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,1)
        jsppol = mpi_enreg%kpt_loc2fbz_sp(dest, ikpt_loc,2)
        if (jkpt > 0 .and. jsppol > 0) then ! dest is treating a true kpt
           
           jkptb = dtorbmag%ikpt_dk(jkpt,bfor,bdir)
           jkptbi = dtorbmag%indkk_f2ibz(jkptb,1)
                       
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
           
        end if ! end check that jkpt > 0 and jsppol > 0

     end if ! test dest .EQ. me and ikpt .GT. 0

  end do ! end loop over dest

  if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
  end if

end subroutine mpicomm_helper
            
!{\src2tex{textfont=tt}}
!!****f* ABINIT/duqhqdu
!! NAME
!! duqhqdu
!!
!! FUNCTION
!! Return i*epsabg\sum_n <\partial_b u_kn|QH_kQ|\partial_g u_kn> where
!! Q projects onto the conduction space.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE


subroutine duqhqdu(atindx1,cg,cnum_duqhqdu,cprj,dtorbmag,dtset,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
     & nattyp,nband_k,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,pwind,pwind_alloc,&
     & rmet,rprimd,ucvol,vectornd,vhartr,vpsp,vxc,with_vectornd,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,nband_k,nfftf,pwind_alloc,with_vectornd
  real(dp),intent(in) :: ucvol
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer, intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(dtset%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
  real(dp), intent(in) :: vhartr(nfftf),vpsp(nfftf),vxc(nfftf,dtset%nspden),xred(3,dtset%natom)
  real(dp), intent(inout) :: vectornd(with_vectornd*nfftf,3)
  real(dp), intent(out) :: cnum_duqhqdu(2,3)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgr_type),intent(in) :: pawfgr
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,bfor,bsigma,countb,countg,countjb,countjg,countk,cpopt,dest,dimffnl
  integer :: getcprj_choice,getcprj_cpopt,getcprj_idir
  integer :: epsabg,exchn2n3d,gdir,gfor,gsigma,ia,iatom,iband,ibs1,ibs2
  integer :: icg,icgb,icgg,icprji,icprjbi,icprjgi,ider,idir,ierr
  integer :: ikg,ikg1,ikpt,ikpt_ctr,ikpt_loc,ikpti,ikptb,ikptbi,ikptg,ikptgi
  integer :: ilm,isppol,istwf_k,itrs,itypat
  integer :: jcgb,jcgg,jcprjbi,jcprjgi,jkpt,jkptb,jkptbi,jkptg,jkptgi,jsppol
  integer :: me,mcg1_k,my_nspinor,n2dim,ncpgr,ndat
  integer :: ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,nkpg,npw_k,npw_k_,npw_kb,npw_kg,nproc,ntotcp
  integer :: optder,prtvol,shiftbd,sij_opt,smatrix_ddkflag,smatrix_job,sourceb,sourceg,spaceComm
  integer :: tagb,tagg,tim_getghc,type_calc
  real(dp) :: arg,deltab,deltag,doti,dotr,ecut_eff,lambda
  complex(dpc) :: cgdijcb,cprefac,duqhqdu_term
  logical :: has_vectornd
  type(gs_hamiltonian_type) :: gs_hamk

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:),kg_k(:,:),nband_dum(:),pwind_kb(:),pwind_kg(:),sflag_k(:)
  real(dp) :: dkb(3),dkbg(3),dkg(3),dtm_k(2),kpoint(3),kpointb(3),kpointg(3),lambdarr(1),rhodum(1)
  real(dp),allocatable :: buffer(:,:)
  real(dp),allocatable :: cg_k(:,:),cg1_kb(:,:),cg1_kg(:,:),cgrvtrial(:,:)
  real(dp),allocatable :: cgqb(:,:),cgqg(:,:),cwavef(:,:),ffnlb(:,:,:,:),ffnlg(:,:,:,:)
  real(dp),allocatable :: ghc(:,:),ghcall(:,:),gsc(:,:),gvnlc(:,:)
  real(dp),allocatable :: kinpw(:),kk_paw(:,:,:),kpg_k(:,:),kpg_k_dummy(:,:),kptnsb(:,:),kptnsball(:,:,:,:),kptnsg(:,:)
  real(dp),allocatable :: ph1d(:,:),ph3d(:,:,:),phkxred(:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_inv(:,:,:),smat_kk(:,:,:)
  real(dp),allocatable :: vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:),vtrial(:,:)
  real(dp),allocatable :: ylm_k(:,:),ylmb(:,:),ylmball(:,:,:,:),ylmgr_k(:,:,:)
  real(dp),allocatable :: ylmgrb(:,:,:),ylmgrball(:,:,:,:,:)
  type(pawcprj_type),allocatable :: cprj_buf(:,:),cprj_k(:,:),cprj_kb(:,:),cprj1_kb(:,:)
  type(pawcprj_type),allocatable :: cprj_kg(:,:),cprj1_kg(:,:),cwaveprj(:,:)

  !----------------------------------------------------

  isppol = 1
  ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
  ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me=mpi_enreg%me_kpt

  ! input parameters for calls to getghc at ikpt
  cpopt = -1
  ndat = 1
  prtvol = 0
  sij_opt = 0
  tim_getghc = 0
  ! getghc: type_calc 3 means kinetic, local only
  ! type_calc 1 means local only
  type_calc = 3
  lambda = zero; lambdarr(1) = zero

  ! input parameters for calls to getcprj
  getcprj_choice = 1 ! just cprj no gradients
  getcprj_cpopt = 0 ! no cprj in memory already
  getcprj_idir = 0 ! gradient directions irrelevant here

  !==== Initialize most of the Hamiltonian ====
  !Allocate all arrays and initialize quantities that do not depend on k and spin.
  !gs_hamk is the normal hamiltonian at k
  call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
       & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
       & paw_ij=paw_ij)

  !---------construct local potential------------------
  ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
  ! nspden=1 is essentially hard-coded in the following line
  vtrial(1:nfftf,1)=vhartr(1:nfftf)+vxc(1:nfftf,1)+vpsp(1:nfftf)

  ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
  call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)

  ABI_ALLOCATE(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
  call fftpac(isppol,mpi_enreg,dtset%nspden,&
       & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)
  
  ABI_DEALLOCATE(cgrvtrial)
  ABI_DEALLOCATE(vtrial)

  ! if vectornd is present, set it up for addition to gs_hamk similarly to how it's done for
  ! vtrial. Note that it must be done for the three directions. Also, the following
  ! code assumes explicitly and implicitly that nvloc = 1. This should eventually be generalized.
   has_vectornd = (with_vectornd .EQ. 1)
   if(has_vectornd) then
     ABI_ALLOCATE(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
     ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
     do adir = 1, 3
        call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,adir))
        call fftpac(isppol,mpi_enreg,dtset%nspden,&
             & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,adir),2)
     end do
     ABI_DEALLOCATE(cgrvtrial)
  end if

  ! add vlocal
  call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)

  ! add vectornd if available
  if(has_vectornd) then
     call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
  end if

  ABI_ALLOCATE(kg_k,(3,dtset%mpw))
  ABI_ALLOCATE(kinpw,(dtset%mpw))

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj1_kb,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj1_kb,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj1_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj1_kg,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)
  n2dim = dtorbmag%nspinor*nband_k
  ntotcp = n2dim*SUM(dimlmn(:))
  if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(dtset%natom,n2dim))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
  end if

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(pwind_kb,(dtset%mpw))
  ABI_ALLOCATE(pwind_kg,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = one; pwnsfac_k(2,:) = zero
  pwnsfac_k(3,:) = one; pwnsfac_k(4,:) = zero

  mcg1_k = dtset%mpw*dtset%nsppol*my_nspinor*nband_k
  ABI_ALLOCATE(cg_k,(2,mcg1_k))
  ABI_ALLOCATE(cg1_kb,(2,mcg1_k))
  ABI_ALLOCATE(cg1_kg,(2,mcg1_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

  optder=0 ! ylm only, no d ylm/dk
  ABI_ALLOCATE(kptnsb,(3,dtset%nkpt))
  ABI_ALLOCATE(kptnsball,(3,dtset%nkpt,3,2))
  ABI_ALLOCATE(ylmball,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang,3,2))
  ABI_ALLOCATE(ylmgrball,(dtset%mpw*dtset%mkmem,3+6*(optder/2),psps%mpsang*psps%mpsang,3,2))
  ABI_ALLOCATE(kptnsg,(3,dtset%nkpt))
  ABI_ALLOCATE(nband_dum,(dtset%nkpt*dtset%nsppol))
  nband_dum(:)=nband_k

  ! set up all ylm's at k+dkb. Notice that this kpt might differ by a periodic wrap-around
  ! from the one determined from dtorbmag%ikpt_dk(ikpt,bfor,bdir)
  ABI_ALLOCATE(ylmb,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang))
  ABI_ALLOCATE(ylmgrb,(dtset%mpw*dtset%mkmem,3+6*(optder/2),psps%mpsang*psps%mpsang))
  do bdir = 1, 3
     do bfor = 1, 2
        bsigma = 3-2*bfor
        dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
        do ikpt_ctr=1, dtset%nkpt
           kptnsball(:,ikpt_ctr,bdir,bfor) = dtset%kptns(:,ikpt_ctr)+dkb(:)
        end do
        kptnsb(:,:) = kptnsball(:,:,bdir,bfor)
        call initylmg(gprimd,kg,kptnsb,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
             & nband_dum,dtset%nkpt,npwarr,dtset%nsppol,optder,rprimd,ylmb,ylmgrb)
        ylmball(:,:,bdir,bfor) = ylmb(:,:)
        ylmgrball(:,:,:,bdir,bfor) = ylmgrb(:,:,:)
     end do
  end do
  ABI_DEALLOCATE(ylmb)
  ABI_DEALLOCATE(ylmgrb)

  smatrix_ddkflag = 1
  itrs = 0
  smatrix_job = 1
  shiftbd = 1

  cnum_duqhqdu(:,:) = zero

  do ikpt_loc = 1,dtorbmag%fmkmem_max

     ikpt=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
     ! if this k and spin are for me do it
     ! if (ikpt1 > 0 .and. isppol > 0) then
     if (ikpt > 0) then

        ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
        kpoint(:)=dtorbmag%fkptns(:,ikpt)
        icprji = dtorbmag%cprjindex(ikpti,isppol)
        npw_k = npwarr(ikpti)
        icg = dtorbmag%cgindex(ikpti,dtset%nsppol)
        ikg = dtorbmag%fkgindex(ikpt)

        countk = npw_k*my_nspinor*nband_k
        cg_k(1:2,1:countk) = cg(1:2,icg+1:icg+countk)
                    
        call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
             &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

        ! this is minimal Hamiltonian information, to apply vlocal and kinetic to |u_kn>
        ! Build basis sphere of plane waves for the k-point
        kg_k(:,:) = 0
        call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg_k,kpoint,1,mpi_enreg,dtset%mpw,npw_k_)
        if (npw_k .NE. npw_k_) then
           write(std_out,'(a)')'JWZ debug duqhqdu npw_k inconsistency'
        end if
        ! Compute kinetic energy at kpt
        kinpw(:) = zero
        call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)
        nkpg = 0
        ABI_ALLOCATE(kpg_k_dummy,(npw_kg,nkpg))
        call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
             &             kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k_dummy,compute_gbound=.TRUE.)
        ABI_DEALLOCATE(kpg_k_dummy)


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
              bsigma = 3-2*bfor
              dkb(1:3) = bsigma*dtorbmag%dkvecs(1:3,bdir)
              deltab = sqrt(DOT_PRODUCT(dkb,dkb))

              if(ikpt > 0) then
                 ikptb = dtorbmag%ikpt_dk(ikpt,bfor,bdir)
                 ikptbi = dtorbmag%indkk_f2ibz(ikptb,1)
                 npw_kb = npwarr(ikptbi)
                 pwind_kb(1:npw_k) = pwind(ikg+1:ikg+npw_k,bfor,bdir)

                 ! for use in mkffnl, differs by wraparounds from ikpt_dk
                 kpointb(:) = kptnsball(:,ikpt,bdir,bfor)

                 ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
                 do ilm=1,psps%mpsang*psps%mpsang
                    ylm_k(1:npw_k,ilm)=ylmball(1+ikg:npw_k+ikg,ilm,bdir,bfor)
                 end do

                 ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
                 do ilm=1,psps%mpsang*psps%mpsang
                    ylmgr_k(1:npw_k,1:3,ilm)=ylmgrball(1+ikg:npw_k+ikg,1:3,ilm,bdir,bfor)
                 end do

                 ABI_ALLOCATE(phkxred,(2,dtset%natom))
                 do ia=1, dtset%natom
                    arg=two_pi*(kpointb(1)*xred(1,ia)+kpointb(2)*xred(2,ia)+kpointb(3)*xred(3,ia))
                    phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
                 end do

                 ABI_ALLOCATE(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
                 call getph(atindx1,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

                 ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))
                 call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,&
                      & npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)

                 nkpg = 3
                 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
                 call mkkpg(kg_k,kpg_k,kpointb,nkpg,npw_k)

                 !      Compute nonlocal form factors ffnl at all (k+b+G_k):
                 dimffnl=1 ! 1 + number of derivatives
                 ABI_ALLOCATE(ffnlb,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
                 ider=0 ! no derivs
                 idir=4 ! derivs in all directions of cart (not used if ider = 0)
                 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlb,psps%ffspl,&
                      & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpointb,psps%lmnmax,&
                      & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                      & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                      & psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

                 ABI_DEALLOCATE(ylm_k)
                 ABI_DEALLOCATE(ylmgr_k)

              end if

              if (ikpt > 0 .AND. isppol > 0) then
                 countb = npw_kb*my_nspinor*nband_k
                 if(allocated(cgqb)) then
                    ABI_DEALLOCATE(cgqb)
                 endif
                 ABI_ALLOCATE(cgqb,(2,countb))
                 call mpicomm_helper(atindx1,bdir,bfor,cg,cgqb,cprj,cprj_kb,dimlmn,dtorbmag,dtset,&
                      & ikpt,ikpt_loc,ikptbi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
                      & nproc,npw_kb,npwarr,spaceComm)
              end if

              if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps
                 ! get covariant |u_{n,k+b}> and associated cprj
                 call overlap_k1k2_paw(cprj_k,cprj_kb,dkb,gprimd,kk_paw,dtorbmag%lmn2max,&
                      &           dtorbmag%lmn_size,dtset%mband,&
                      &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                 sflag_k=0
                 cg1_kb(:,:) = zero
                 ! cg1_kb will hold |\tilde{u}_{n,k+b}>
                 call smatrix(cg_k,cgqb,cg1_kb,smatrix_ddkflag,dtm_k,0,0,itrs,smatrix_job,nband_k,&
                      &           mcg1_k,mcg1_k,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kb,my_nspinor,&
                      &           pwind_kb,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)
                 ! cprj1_kb will hold cprj for cg1_kb
                 call covar_cprj(cprj_kb,cprj1_kb,dtset,nband_k,pawtab,smat_inv)

                 ABI_DEALLOCATE(cgqb)
                 
                 ABI_ALLOCATE(cwavef,(2,npw_k))
                 ABI_ALLOCATE(ghc,(2,npw_k))
                 ABI_ALLOCATE(gsc,(2,npw_k))
                 ABI_ALLOCATE(gvnlc,(2,npw_k))
                 if(allocated(ghcall)) then
                    ABI_DEALLOCATE(ghcall)
                 end if
                 ABI_ALLOCATE(ghcall,(2,npw_k*nband_k))
                 do iband = 1, nband_k
                    ibs1 = (iband-1)*npw_k+1
                    ibs2 = iband*npw_k
                    cwavef(1:2,1:npw_k) = cg1_kb(1:2,ibs1:ibs2)
                    call pawcprj_get(atindx1,cwaveprj,cprj1_kb,dtset%natom,iband,0,ikptb,0,isppol,dtset%mband,&
                         & dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
                    call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
                         & prtvol,sij_opt,tim_getghc,type_calc)
                    ghcall(1:2,ibs1:ibs2) = ghc(1:2,1:npw_k)
                 end do

                 call pawcprj_set_zero(cprj1_kb)
                 do iband = 1, nband_k
                    ibs1 = (iband-1)*npw_k+1
                    ibs2 = iband*npw_k
                    cwavef(1:2,1:npw_k) = cg1_kb(1:2,ibs1:ibs2)
                    call getcprj(getcprj_choice,getcprj_cpopt,cwavef,cwaveprj,ffnlb,&
                         & getcprj_idir,psps%indlmn,istwf_k,kg_k,kpg_k,kpointb,psps%lmnmax,&
                         & dtset%mgfft,mpi_enreg,&
                         & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                         & phkxred,ph1d,ph3d,ucvol,psps%useylm)

                    call pawcprj_put(atindx1,cwaveprj,cprj1_kb,dtset%natom,&
                         & iband,0,ikpt,0,isppol,nband_k,dtset%mkmem,&
                         & dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0,&
                         & mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
                 end do

                 ABI_DEALLOCATE(cwavef)
                 ABI_DEALLOCATE(ghc)
                 ABI_DEALLOCATE(gsc)
                 ABI_DEALLOCATE(gvnlc)

                 ABI_DEALLOCATE(kpg_k)
                 ABI_DEALLOCATE(phkxred)
                 ABI_DEALLOCATE(ph1d)
                 ABI_DEALLOCATE(ph3d)
                 ABI_DEALLOCATE(ffnlb)

              end if
              
              do gfor = 1, 2
                 gsigma=3-2*gfor
                 dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
                 deltag = sqrt(DOT_PRODUCT(dkg,dkg))

                 cprefac = j_dpc*epsabg*bsigma*gsigma/(two*deltab*two*deltag)

                 if (ikpt > 0) then
                    ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
                    ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)
                    npw_kg = npwarr(ikptgi)
                    pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

                    ! for use in mkffnl, differs by wraparounds from ikpt_dk
                    kpointg(:) = kptnsball(:,ikpt,gdir,gfor)

                    ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
                    do ilm=1,psps%mpsang*psps%mpsang
                       ylm_k(1:npw_k,ilm)=ylmball(1+ikg:npw_k+ikg,ilm,gdir,gfor)
                    end do

                    ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
                    do ilm=1,psps%mpsang*psps%mpsang
                       ylmgr_k(1:npw_k,1:3,ilm)=ylmgrball(1+ikg:npw_k+ikg,1:3,ilm,gdir,gfor)
                    end do

                    ABI_ALLOCATE(phkxred,(2,dtset%natom))
                    do ia=1, dtset%natom
                       arg=two_pi*(kpointg(1)*xred(1,ia)+kpointg(2)*xred(2,ia)+kpointg(3)*xred(3,ia))
                       phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
                    end do

                    ABI_ALLOCATE(ph1d,(2,dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)))
                    call getph(atindx1,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)

                    ABI_ALLOCATE(ph3d,(2,npw_k,dtset%natom))
                    call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,&
                         & npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)

                    nkpg = 3
                    ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
                    call mkkpg(kg_k,kpg_k,kpointg,nkpg,npw_k)

                    !      Compute nonlocal form factors ffnl at all (k+b+G_k):
                    dimffnl=1 ! 1 + number of derivatives
                    ABI_ALLOCATE(ffnlg,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
                    ider=0 ! no derivs
                    idir=4 ! derivs in all directions of cart (not used if ider = 0)
                    call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlg,psps%ffspl,&
                         & gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpointg,psps%lmnmax,&
                         & psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
                         & npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
                         & psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

                    ABI_DEALLOCATE(ylm_k)
                    ABI_DEALLOCATE(ylmgr_k)
                 end if
                    
                 if (ikpt > 0 .AND. isppol > 0) then
                    countg = npw_kg*my_nspinor*nband_k
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    endif
                    ABI_ALLOCATE(cgqg,(2,countg))
                    call mpicomm_helper(atindx1,gdir,gfor,cg,cgqg,cprj,cprj_kg,dimlmn,dtorbmag,dtset,&
                         & ikpt,ikpt_loc,ikptgi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
                         & nproc,npw_kg,npwarr,spaceComm)
                 end if

                 if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps

                    ! get covariant |u_{n,k+g}> and associated cprj
                    call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,&
                         &           dtorbmag%lmn_size,dtset%mband,&
                         &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                    sflag_k=0
                    cg1_kg(:,:) = zero
                    ! cg1_kg will hold |\tilde{u}_{n,k+g}>
                    call smatrix(cg_k,cgqg,cg1_kg,smatrix_ddkflag,dtm_k,0,0,itrs,smatrix_job,nband_k,&
                         &           mcg1_k,mcg1_k,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                         &           pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)
                    if(allocated(cgqg)) then
                       ABI_DEALLOCATE(cgqg)
                    end if

                    ABI_ALLOCATE(cwavef,(2,npw_k))
                    call pawcprj_set_zero(cprj1_kg)
                    do iband = 1, nband_k
                       ibs1 = (iband-1)*npw_k+1
                       ibs2 = iband*npw_k
                       cwavef(1:2,1:npw_k) = cg1_kg(1:2,ibs1:ibs2)
                       call getcprj(getcprj_choice,getcprj_cpopt,cwavef,cwaveprj,ffnlg,&
                            & getcprj_idir,psps%indlmn,istwf_k,kg_k,kpg_k,kpointg,psps%lmnmax,&
                            & dtset%mgfft,mpi_enreg,&
                            & dtset%natom,nattyp,dtset%ngfft,dtset%nloalg,npw_k,dtset%nspinor,dtset%ntypat,&
                            & phkxred,ph1d,ph3d,ucvol,psps%useylm)

                       call pawcprj_put(atindx1,cwaveprj,cprj1_kg,dtset%natom,&
                            & iband,0,ikpt,0,isppol,nband_k,dtset%mkmem,&
                            & dtset%natom,1,nband_k,dimlmn,dtset%nspinor,dtset%nsppol,0,&
                            & mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
                    end do
                    ABI_DEALLOCATE(cwavef)
                    
                    do iband = 1, nband_k
                       ibs1=(iband-1)*npw_k+1
                       ibs2=iband*npw_k
                       dotr=DOT_PRODUCT(cg1_kg(1,ibs1:ibs2),ghcall(1,ibs1:ibs2)) + &
                            & DOT_PRODUCT(cg1_kg(2,ibs1:ibs2),ghcall(2,ibs1:ibs2))
                       doti=DOT_PRODUCT(cg1_kg(1,ibs1:ibs2),ghcall(2,ibs1:ibs2)) - &
                            & DOT_PRODUCT(cg1_kg(2,ibs1:ibs2),ghcall(1,ibs1:ibs2))

                       ! compute onsite contribution due to paw_ij%dij
                         call cpg_dij_cpb(cgdijcb,cprj1_kb,cprj1_kg,dtset,iband,iband,dtorbmag%nspinor,paw_ij,pawtab)
        
                       ! accumulate i*epsabg*\sum_occ [<d_bdir u|QHQ|d_gdir u>]
                       duqhqdu_term = cprefac*(cmplx(dotr,doti)+cgdijcb)
                       ! duqhqdu_term = cprefac*cmplx(dotr,doti)
                       ! duqhqdu_term = cprefac*cgdijcb
                       cnum_duqhqdu(1,adir) = cnum_duqhqdu(1,adir) + real(duqhqdu_term) 
                       cnum_duqhqdu(2,adir) = cnum_duqhqdu(2,adir) + aimag(duqhqdu_term)

                    end do ! end loop over iband

                    ABI_DEALLOCATE(kpg_k)
                    ABI_DEALLOCATE(phkxred)
                    ABI_DEALLOCATE(ph1d)
                    ABI_DEALLOCATE(ph3d)
                    ABI_DEALLOCATE(ffnlg)

                 end if ! end check on ikpt > 0
                 
              end do ! end loop over gfor
              if(allocated(ghcall)) then
                 ABI_DEALLOCATE(ghcall)
              end if
           end do ! end loop over bfor
        end do ! end loop over epsabg
     end do ! end loop over adir

  end do ! end loop over ikpt_loc

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(cnum_duqhqdu,spaceComm,ierr)
  end if

  call gs_hamk%free()
  ABI_DEALLOCATE(vlocal)
  if(has_vectornd) then
     ABI_DEALLOCATE(vectornd_pac)
  end if

  ABI_DEALLOCATE(kg_k)
  ABI_DEALLOCATE(kinpw)
  
  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cprj_kb)
  ABI_DATATYPE_DEALLOCATE(cprj_kb)
  call pawcprj_free(cprj1_kb)
  ABI_DATATYPE_DEALLOCATE(cprj1_kb)
  call pawcprj_free(cprj_kg)
  ABI_DATATYPE_DEALLOCATE(cprj_kg)
  call pawcprj_free(cprj1_kg)
  ABI_DATATYPE_DEALLOCATE(cprj1_kg)
  call pawcprj_free(cwaveprj)
  ABI_DATATYPE_DEALLOCATE(cwaveprj)
  if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
  end if

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(pwind_kb)
  ABI_DEALLOCATE(pwind_kg)
  ABI_DEALLOCATE(cg1_kb)
  ABI_DEALLOCATE(cg1_kg)
  ABI_DEALLOCATE(cg_k)
  ABI_DEALLOCATE(pwnsfac_k)

  ABI_DEALLOCATE(kptnsb)
  ABI_DEALLOCATE(kptnsball)
  ABI_DEALLOCATE(ylmball)
  ABI_DEALLOCATE(ylmgrball)
  ABI_DEALLOCATE(kptnsg)
  ABI_DEALLOCATE(nband_dum)

end subroutine duqhqdu
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/udsqdu
!! NAME
!! udsqdu
!!
!! FUNCTION
!! Return i*epsabg\sum_n<u_kn|\partial_b S Q |\partial_g u_kn> where
!! Q projects onto the conduction space. This term contributes to
!! the Chern number (integral over Berry curvature).
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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
!!
!! SOURCE


subroutine udsqdu(atindx1,cg,cprj,dtorbmag,dtset,energies,gmet,gprimd,kg,&
     & mcg,mcprj,mpi_enreg,nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,psps,&
     pwind,pwind_alloc,rmet,rprimd,udsqduchern,udsqdumag,xred,ylm,ylmgr)

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
  integer,intent(in) :: npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),xred(3,dtset%natom)
  real(dp), intent(out) :: udsqduchern(2,3),udsqdumag(2,3)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: energies(nband_k,dtset%nkpt)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,countg,countjg,countk,dest,dimph1d,dimffnl,dir_tmp,exchn2n3d,epsabg
  integer :: gdir,gfor,gsigma,ia,iatom,iband,ibs1,ibs2
  integer :: icg,icgg,icprji,icprjg,icprjgi,ider,idir,ierr
  integer :: ikg,ikg1,ikpt,ikpt_loc,ikpti,ikptg,ikptgi,ilm,isppol,istwf_k,itrs,itypat
  integer :: jcgg,jcprjgi,jkpt,jkptg,jkptgi,jsppol
  integer :: mcg1_k,me,my_nspinor,n2dim,ncpgr,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6
  integer :: nkpg,nproc,npw_k,npw_k_,npw_kg
  integer :: nonlop_choice,nonlop_cpopt,nonlop_nnlout,nonlop_ndat
  integer :: nonlop_paw_opt,nonlop_signs,ntotcp
  integer :: shiftbd,smatrix_ddkflag,smatrix_job,sourceg,spaceComm,tagg,tim_nonlop
  real(dp) :: arg,deltab,deltag,doti,dotr,ecut_eff,ENK
  complex(dpc) :: cprefac,udsqduchern_term,udsqdumag_term
  type(gs_hamiltonian_type) :: gs_hamk

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:),kg_k(:,:),pwind_kg(:),sflag_k(:)
  real(dp) :: dkg(3),dtm_k(2),kpoint(3),nonlop_lambda(1)
  real(dp),allocatable :: buffer(:,:),cg_k(:,:),cg1_kg(:,:),cgqg(:,:),cwavef(:,:),ffnl_k(:,:,:,:)
  real(dp),allocatable :: kk_paw(:,:,:),kpg_k(:,:),nonlop_enlout(:)
  real(dp),allocatable :: phkxred(:,:),ph1d(:,:),ph3d(:,:,:),pwnsfac_k(:,:)
  real(dp),allocatable :: smat_inv(:,:,:),smat_kk(:,:,:),svect(:,:,:),svectout(:,:),vectout(:,:)
  real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
  type(pawcprj_type),allocatable :: cprj_buf(:,:),cprj_k(:,:),cprj_kg(:,:),cwaveprj(:,:)

  !----------------------------------------------------

  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cprj_kg,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_kg,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)
  n2dim = dtorbmag%nspinor*nband_k
  ntotcp = n2dim*SUM(dimlmn(:))
  if (nproc>1) then
     ABI_DATATYPE_ALLOCATE(cprj_buf,(dtset%natom,n2dim))
     call pawcprj_alloc(cprj_buf,ncpgr,dimlmn)
  end if

  ABI_ALLOCATE(kk_paw,(2,dtset%mband,dtset%mband))
  ABI_ALLOCATE(sflag_k,(nband_k))
  ABI_ALLOCATE(pwind_kg,(dtset%mpw))
  ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
  pwnsfac_k(1,:) = one; pwnsfac_k(2,:) = zero
  pwnsfac_k(3,:) = one; pwnsfac_k(4,:) = zero

  ABI_ALLOCATE(kg_k,(3,dtset%mpw))

  mcg1_k = dtset%mpw*dtset%nsppol*my_nspinor*nband_k
  ABI_ALLOCATE(cg1_kg,(2,mcg1_k))
  ABI_ALLOCATE(cg_k,(2,mcg1_k))
  ABI_ALLOCATE(smat_inv,(2,nband_k,nband_k))
  ABI_ALLOCATE(smat_kk,(2,nband_k,nband_k))

  smatrix_ddkflag = 1
  itrs = 0
  smatrix_job = 1
  shiftbd = 1

  !==== Initialize most of the Hamiltonian ====
  !Allocate all arrays and initialize quantities that do not depend on k and spin.
  !gs_hamk is the normal hamiltonian at k
  ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
  ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
  istwf_k = 1
  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; ikg1 = 0

  ABI_ALLOCATE(phkxred,(2,dtset%natom))
  dimph1d=dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)
  ABI_ALLOCATE(ph1d,(2,dimph1d))
  call getph(atindx1,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)
  call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
       & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
       & paw_ij=paw_ij)
  call gs_hamk%load_spin(isppol,with_nonlocal=.true.)

  ! nonlop parameters
  nonlop_choice = 5 ! first derivative wrt k
  nonlop_cpopt = -1
  nonlop_nnlout = 1 ! size of enlout, not used in call
  ABI_ALLOCATE(nonlop_enlout,(nonlop_nnlout))
  nonlop_lambda(1) = 0.0 ! shift for eigenvalues, not used
  nonlop_ndat = 1 ! number of wavefunctions to apply nonlop
  nonlop_paw_opt = 3 ! use Sij matrix
  nonlop_signs = 2 ! apply to function in recip space
  tim_nonlop = 0 ! timing not used

  udsqduchern(:,:) = zero
  udsqdumag(:,:) = zero

  do ikpt_loc = 1,dtorbmag%fmkmem_max

     ikpt=mpi_enreg%kpt_loc2fbz_sp(me, ikpt_loc,1)
     ! if this k and spin are for me do it
     if (ikpt > 0) then

        ikpti = dtorbmag%indkk_f2ibz(ikpt,1)
        icprji = dtorbmag%cprjindex(ikpti,isppol)
        npw_k = npwarr(ikpti)
        icg = dtorbmag%cgindex(ikpti,dtset%nsppol)
        ikg = dtorbmag%fkgindex(ikpt)

        countk = npw_k*my_nspinor*nband_k
        cg_k(1:2,1:countk) = cg(1:2,icg+1:icg+countk)
                    
        call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprji,ikpti,0,isppol,dtset%mband,&
             &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
        
        kpoint(:)=dtset%kptns(:,ikpt)
        ! Build basis sphere of plane waves for the k-point
        kg_k(:,:) = 0
        call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg_k,kpoint,1,mpi_enreg,dtset%mpw,npw_k_)
        if (npw_k .NE. npw_k_) then
           write(std_out,'(a)')'JWZ debug udsqdu npw_k inconsistency'
        end if

        ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
        ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
        do ilm=1,psps%mpsang*psps%mpsang
           ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
           ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
        end do

        nkpg = 3
        ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
        call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

        ! Compute nonlocal form factors ffnl at all (k+G):
        ider=1 ! want ffnl and 1st derivative
        idir=4 ! d ffnl/ dk 
        dimffnl=4 ! 1 + number of derivatives
        ABI_ALLOCATE(ffnl_k,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
        call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,&
             &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
             &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
             &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
             &         psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
        
        ! Load k-dependent part in the Hamiltonian datastructure
        !  - Compute 3D phase factors
        !  - Prepare various tabs in case of band-FFT parallelism
        !  - Load k-dependent quantities in the Hamiltonian
        ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
        do ia=1, dtset%natom
           arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
           phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
        end do

        call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)
        
        call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
             &         kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,compute_gbound=.TRUE.)

        ABI_ALLOCATE(cwavef,(2,npw_k))
        ABI_ALLOCATE(svectout,(2,npw_k))
        ABI_ALLOCATE(vectout,(2,npw_k))
        ABI_ALLOCATE(svect,(2,3,npw_k*nband_k))
        
        do bdir = 1, 3
           do iband = 1, nband_k

              ibs1=(iband-1)*npw_k+1
              ibs2=iband*npw_k
                    
              cwavef(1:2,1:npw_k) = cg_k(1:2,ibs1:ibs2)
              call pawcprj_get(atindx1,cwaveprj,cprj_k,dtset%natom,iband,0,ikpt,0,isppol,dtset%mband,&
                   &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
                    
              ! compute dS/dk_bdir|u_nk>
              call nonlop(nonlop_choice,nonlop_cpopt,cwaveprj,nonlop_enlout,gs_hamk,bdir,nonlop_lambda,&
                   & mpi_enreg,nonlop_ndat,nonlop_nnlout,nonlop_paw_opt,nonlop_signs,svectout,&
                   & tim_nonlop,cwavef,vectout)

              svect(1:2,bdir,ibs1:ibs2) = svectout(1:2,1:npw_k)
              
           end do ! end loop over iband
        end do ! end loop over bdir
        ABI_DEALLOCATE(cwavef)
        ABI_DEALLOCATE(svectout)
        ABI_DEALLOCATE(vectout)

     end if ! end check that ikpt > 0

     do adir = 1, 3
        do epsabg = 1, -1, -2
           if (epsabg .EQ. 1) then
              bdir = modulo(adir,3)+1
              gdir = modulo(adir+1,3)+1
           else
              bdir = modulo(adir+1,3)+1
              gdir = modulo(adir,3)+1
           end if
        
           do gfor = 1, 2
              gsigma = 3-2*gfor
              dkg(1:3) = gsigma*dtorbmag%dkvecs(1:3,gdir)
              deltag = sqrt(DOT_PRODUCT(dkg,dkg))

              ikptg = dtorbmag%ikpt_dk(ikpt,gfor,gdir)
              ikptgi = dtorbmag%indkk_f2ibz(ikptg,1)
              npw_kg = npwarr(ikptgi)
              pwind_kg(1:npw_k) = pwind(ikg+1:ikg+npw_k,gfor,gdir)

              cprefac = j_dpc*epsabg*gsigma/(two*deltag)

              if (ikpt > 0 .AND. isppol > 0) then
                 countg = npw_kg*my_nspinor*nband_k
                 if(allocated(cgqg)) then
                    ABI_DEALLOCATE(cgqg)
                 endif
                 ABI_ALLOCATE(cgqg,(2,countg))
                 call mpicomm_helper(atindx1,gdir,gfor,cg,cgqg,cprj,cprj_kg,dimlmn,dtorbmag,dtset,&
                      & ikpt,ikpt_loc,ikptgi,isppol,mcg,mcprj,me,mpi_enreg,my_nspinor,nband_k,&
                      & nproc,npw_kg,npwarr,spaceComm)
              end if

              if (ikpt > 0 .and. isppol > 0) then ! if I am treating a kpt, compute the overlaps

                 ! get covariant |u_{n,k+g}> and associated cprj
                 call overlap_k1k2_paw(cprj_k,cprj_kg,dkg,gprimd,kk_paw,dtorbmag%lmn2max,&
                      &           dtorbmag%lmn_size,dtset%mband,&
                      &           dtset%natom,my_nspinor,dtset%ntypat,pawang,pawrad,pawtab,dtset%typat,xred)
                 sflag_k=0
                 cg1_kg(:,:) = zero
                 ! cg1_kg will hold |\tilde{u}_{n,k+g}>
                 call smatrix(cg_k,cgqg,cg1_kg,smatrix_ddkflag,dtm_k,0,0,itrs,smatrix_job,nband_k,&
                      &           mcg1_k,mcg1_k,mcg1_k,1,dtset%mpw,nband_k,nband_k,npw_k,npw_kg,my_nspinor,&
                      &           pwind_kg,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_kk,kk_paw,psps%usepaw)

                 do iband = 1, nband_k
                    
                    ibs1=(iband-1)*npw_k+1
                    ibs2=iband*npw_k
                    
                    dotr=DOT_PRODUCT(svect(1,bdir,ibs1:ibs2),cg1_kg(1,ibs1:ibs2))+&
                         & DOT_PRODUCT(svect(2,bdir,ibs1:ibs2),cg1_kg(2,ibs1:ibs2))
                    doti=DOT_PRODUCT(svect(1,bdir,ibs1:ibs2),cg1_kg(2,ibs1:ibs2))-&
                         & DOT_PRODUCT(svect(2,bdir,ibs1:ibs2),cg1_kg(1,ibs1:ibs2))
                    
                    ! accumulate i*epsabg*ENK\sum_occ [<u|dS_bdir Q|d_gdir u>]
                    ENK = energies(iband,ikpt)
                    udsqduchern_term = cprefac*cmplx(dotr,doti)
                    udsqdumag_term = udsqduchern_term*ENK

                    udsqduchern(1,adir) = udsqduchern(1,adir) + real(udsqduchern_term)
                    udsqduchern(2,adir) = udsqduchern(2,adir) + aimag(udsqduchern_term)
                    udsqdumag(1,adir) = udsqdumag(1,adir) + real(udsqdumag_term)
                    udsqdumag(2,adir) = udsqdumag(2,adir) + aimag(udsqdumag_term)

                 end do ! end loop over iband

                 if(allocated(cgqg)) then
                    ABI_DEALLOCATE(cgqg)
                 end if

              end if ! end check that ikpt > 0
           end do ! end loop for gfor
        end do ! end loop over epsabg
     end do ! end loop over adir
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ffnl_k)
     ABI_DEALLOCATE(ph3d)
     ABI_DEALLOCATE(svect)
  end do ! end loop over ikpt_loc

  ! accumulate result from all processes
  if(nproc > 1) then
     call xmpi_sum(udsqduchern,spaceComm,ierr)
     call xmpi_sum(udsqdumag,spaceComm,ierr)
  end if

  call gs_hamk%free()

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cprj_kg)
  ABI_DATATYPE_DEALLOCATE(cprj_kg)
  call pawcprj_free(cwaveprj)
  ABI_DATATYPE_DEALLOCATE(cwaveprj)
  if (nproc>1) then
     call pawcprj_free(cprj_buf)
     ABI_DATATYPE_DEALLOCATE(cprj_buf)
  end if

  ABI_DEALLOCATE(kk_paw)
  ABI_DEALLOCATE(sflag_k)
  ABI_DEALLOCATE(pwind_kg)
  ABI_DEALLOCATE(pwnsfac_k)
  ABI_DEALLOCATE(kg_k)
  ABI_DEALLOCATE(cg1_kg)
  ABI_DEALLOCATE(cg_k)
  ABI_DEALLOCATE(smat_inv)
  ABI_DEALLOCATE(smat_kk)
  ABI_DEALLOCATE(phkxred)
  ABI_DEALLOCATE(ph1d)
  ABI_DEALLOCATE(nonlop_enlout)

end subroutine udsqdu
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/udsdsu
!! NAME
!! udsdsu
!!
!! FUNCTION
!! Return i*epsabg\sum_{n,n}<u_kn|\partial_b S|u_kn'><u_kn|\partial_g S|u_kn>E_nk
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group (JWZ)
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


subroutine udsdsu(atindx1,cg,cnum_udsdsu,cprj,dtorbmag,dtset,energies,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
     & nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,psps,rmet,rprimd,xred,ylm,ylmgr)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: mcg,mcprj,nband_k
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(orbmag_type), intent(inout) :: dtorbmag
  type(pawang_type),intent(in) :: pawang
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pseudopotential_type),intent(in) :: psps

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(in) :: cg(2,mcg),energies(nband_k,dtset%nkpt),gmet(3,3),gprimd(3,3)
  real(dp), intent(in) :: rmet(3,3),rprimd(3,3),xred(3,dtset%natom)
  real(dp), intent(out) :: cnum_udsdsu(2,3)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: adir,bdir,dimph1d,dimffnl,dir_tmp,exchn2n3d,epsabg
  integer :: gdir,ia,iatom,iband,ibs1,ibs2
  integer :: icg,icprj,ider,idir,ierr
  integer :: ikg,ikg1,ikpt,ikpt_loc,ilm,isppol,istwf_k,itypat,jband
  integer :: me,my_nspinor,ncpgr,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6
  integer :: nkpg,npw_k,npw_k_
  integer :: nonlop_choice,nonlop_cpopt,nonlop_nnlout,nonlop_ndat,nonlop_paw_opt,nonlop_signs
  integer :: nproc,shiftbd,smatrix_ddkflag,smatrix_job,spaceComm,tim_nonlop
  real(dp) :: arg,deltab,deltag,doti,dotr,ecut_eff,ENK
  complex(dpc) :: udsdsu_term,ujdsbu,ujdsgu
  type(gs_hamiltonian_type) :: gs_hamk

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:),kg_k(:,:)
  real(dp) :: kpoint(3),nonlop_lambda(1)
  real(dp),allocatable :: cwavef(:,:),cwavefp(:,:),ffnl_k(:,:,:,:)
  real(dp),allocatable :: kpg_k(:,:),nonlop_enlout(:)
  real(dp),allocatable :: phkxred(:,:),ph1d(:,:),ph3d(:,:,:)
  real(dp),allocatable :: svect(:,:,:),svectout(:,:),svectoutp(:,:),vectout(:,:)
  real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
  type(pawcprj_type),allocatable :: cprj_k(:,:),cwaveprj(:,:)

  !----------------------------------------------------

  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,dtorbmag%nspinor*dtset%mband))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
  ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
  call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

  !==== Initialize most of the Hamiltonian ====
  !Allocate all arrays and initialize quantities that do not depend on k and spin.
  !gs_hamk is the normal hamiltonian at k
  ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
  ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
  istwf_k = 1
  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; ikg1 = 0

  ABI_ALLOCATE(phkxred,(2,dtset%natom))
  dimph1d=dtset%natom*(2*(ngfft1+ngfft2+ngfft3)+3)
  ABI_ALLOCATE(ph1d,(2,dimph1d))
  call getph(atindx1,dtset%natom,ngfft1,ngfft2,ngfft3,ph1d,xred)
  call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
       & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
       & paw_ij=paw_ij)
  call gs_hamk%load_spin(isppol,with_nonlocal=.true.)

  ABI_ALLOCATE(kg_k,(3,dtset%mpw))

  ! nonlop parameters
  nonlop_choice = 5 ! first derivative wrt k
  nonlop_cpopt = -1
  nonlop_nnlout = 1 ! size of enlout, not used in call
  ABI_ALLOCATE(nonlop_enlout,(nonlop_nnlout))
  nonlop_lambda(1) = 0.0 ! shift for eigenvalues, not used
  nonlop_ndat = 1 ! number of wavefunctions to apply nonlop
  nonlop_paw_opt = 3 ! use Sij matrix
  nonlop_signs = 2 ! apply to function in recip space
  tim_nonlop = 0 ! timing not used

  cnum_udsdsu(:,:) = zero
  icg = 0
  ikg = 0
  icprj = 0
  do ikpt = 1, dtset%nkpt
     
     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

     kpoint(:)=dtset%kptns(:,ikpt)
     npw_k = npwarr(ikpt)

     ! Build basis sphere of plane waves for the k-point
     kg_k(:,:) = 0
     call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg_k,kpoint,1,mpi_enreg,dtset%mpw,npw_k_)
     if (npw_k .NE. npw_k_) then
        write(std_out,'(a)')'JWZ debug udsdsu npw_k inconsistency'
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
     ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
     do ilm=1,psps%mpsang*psps%mpsang
        ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
        ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
     end do

     nkpg = 3
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

     ! Compute nonlocal form factors ffnl at all (k+G):
     ider=1 ! want ffnl and 1st derivative
     idir=4 ! d ffnl/ dk 
     dimffnl=4 ! 1 + number of derivatives
     ABI_ALLOCATE(ffnl_k,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,&
          &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
          &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
          &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
          &         psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

     ! Load k-dependent part in the Hamiltonian datastructure
     !  - Compute 3D phase factors
     !  - Prepare various tabs in case of band-FFT parallelism
     !  - Load k-dependent quantities in the Hamiltonian
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
     do ia=1, dtset%natom
        arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
        phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
     end do

     call ph1d3d(1,dtset%natom,kg_k,dtset%natom,dtset%natom,npw_k,ngfft1,ngfft2,ngfft3,phkxred,ph1d,ph3d)
     
     call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
          &         kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,compute_gbound=.TRUE.)
     
     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)
     

     ABI_ALLOCATE(cwavef,(2,npw_k))
     ABI_ALLOCATE(svectout,(2,npw_k))
     ABI_ALLOCATE(vectout,(2,npw_k))
     ABI_ALLOCATE(svect,(2,3,npw_k*nband_k))
     do adir = 1, 3
        do iband = 1, nband_k
           ibs1=(iband-1)*npw_k+1
           ibs2=iband*npw_k
           cwavef(1:2,1:npw_k) = cg(1:2,icg+ibs1:icg+ibs2)
           call pawcprj_get(atindx1,cwaveprj,cprj_k,dtset%natom,iband,0,ikpt,0,isppol,dtset%mband,&
                &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
           ! compute dS/dk_adir|u_nk>, store in svectout
           call nonlop(nonlop_choice,nonlop_cpopt,cwaveprj,nonlop_enlout,gs_hamk,adir,nonlop_lambda,&
                & mpi_enreg,nonlop_ndat,nonlop_nnlout,nonlop_paw_opt,nonlop_signs,svectout,&
                & tim_nonlop,cwavef,vectout)
           svect(1:2,adir,ibs1:ibs2) = svectout(1:2,1:npw_k)
        end do
     end do
     ABI_DEALLOCATE(cwavef)
     ABI_DEALLOCATE(svectout)
     ABI_DEALLOCATE(vectout)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ffnl_k)
     ABI_DEALLOCATE(ph3d)

     do adir = 1, 3
        
        do epsabg = 1, -1, -2
           if (epsabg .EQ. 1) then
              bdir = modulo(adir,3)+1
              gdir = modulo(adir+1,3)+1
           else
              bdir = modulo(adir+1,3)+1
              gdir = modulo(adir,3)+1
           end if

           ABI_ALLOCATE(svectout,(2,npw_k))
           ABI_ALLOCATE(svectoutp,(2,npw_k))
           ABI_ALLOCATE(cwavefp,(2,npw_k))
           do iband = 1, nband_k

              ENK = energies(iband,ikpt)
              svectout(1:2,1:npw_k) = svect(1:2,gdir,(iband-1)*npw_k+1:iband*npw_k)
              svectoutp(1:2,1:npw_k) = svect(1:2,bdir,(iband-1)*npw_k+1:iband*npw_k)

              do jband = 1, nband_k

                 cwavefp(1:2,1:npw_k) = cg(1:2,icg+(jband-1)*npw_k+1:icg+jband*npw_k)

                 dotr=DOT_PRODUCT(cwavefp(1,:),svectout(1,:))+DOT_PRODUCT(cwavefp(2,:),svectout(2,:))
                 doti=DOT_PRODUCT(cwavefp(1,:),svectout(2,:))-DOT_PRODUCT(cwavefp(2,:),svectout(1,:))

                 ujdsgu = cmplx(dotr,doti)

                 dotr=DOT_PRODUCT(cwavefp(1,:),svectoutp(1,:))+DOT_PRODUCT(cwavefp(2,:),svectoutp(2,:))
                 doti=DOT_PRODUCT(cwavefp(1,:),svectoutp(2,:))-DOT_PRODUCT(cwavefp(2,:),svectoutp(1,:))

                 ujdsbu = cmplx(dotr,doti)

                 ! accumulate i*epsabg*ENK\sum_occ [<u_nk|dS_bdir|u_n'k><u_n'k|dS_gdir|u_nk>]
                 udsdsu_term = j_dpc*epsabg*ENK*CONJG(ujdsbu)*ujdsgu
                 cnum_udsdsu(1,adir) = cnum_udsdsu(1,adir) + real(udsdsu_term)
                 cnum_udsdsu(2,adir) = cnum_udsdsu(2,adir) + aimag(udsdsu_term)

              end do !end loop over jband
           end do ! end loop over iband
           ABI_DEALLOCATE(svectout)
           ABI_DEALLOCATE(svectoutp)
           ABI_DEALLOCATE(cwavefp)
           
        end do ! end loop over epsabg
     end do ! end loop over adir

     icg = icg + npw_k*nband_k
     ikg = ikg + npw_k
     icprj = icprj + nband_k

     ABI_DEALLOCATE(svect)

  end do ! end loop over ikpt

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(cnum_udsdsu,spaceComm,ierr)
  end if

  call gs_hamk%free()

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)
  call pawcprj_free(cwaveprj)
  ABI_DATATYPE_DEALLOCATE(cwaveprj)

  ABI_DEALLOCATE(phkxred)
  ABI_DEALLOCATE(ph1d)
  ABI_DEALLOCATE(kg_k)
  ABI_DEALLOCATE(nonlop_enlout)

end subroutine udsdsu
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_onsite_l_k
!! NAME
!! make_onsite_l_k
!!
!! FUNCTION
!! Compute 1/2 <L_R> onsite contribution to orbital magnetization at given k point and idir
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

subroutine make_onsite_l_k(cprj_k,dtset,idir,nband_k,onsite_l_k,pawrad,pawtab)

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
                 cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
                 cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)
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
!!****f* ABINIT/make_onsite_l
!! NAME
!! make_onsite_l
!!
!! FUNCTION
!! Compute 1/2 <L_R> onsite contribution to orbital magnetization in direction idir
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

subroutine make_onsite_l(atindx1,cprj,dtset,idir,mcprj,mpi_enreg,nband_k,onsite_l,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,mcprj,nband_k
  complex(dpc),intent(out) :: onsite_l
  type(MPI_type), intent(inout) :: mpi_enreg
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx1(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: icprj,ierr,ikpt,ikpt_loc,isppol,me,my_nspinor,ncpgr,nproc,spaceComm
  complex(dpc) :: onsite_l_k

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:)
  type(pawcprj_type),allocatable :: cprj_k(:,:)

  ! ***********************************************************************

  ! TODO: generalize to nsppol > 1
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  !Init MPI
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,nband_k))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)

  ! loop over kpts on each processor
  onsite_l = czero
  ikpt_loc = 0
  ! loop over all the kpts
  do ikpt = 1, dtset%nkpt

     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

     ikpt_loc = ikpt_loc + 1
     icprj= (ikpt_loc - 1)*nband_k
     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt_loc,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

     call make_onsite_l_k(cprj_k,dtset,idir,nband_k,onsite_l_k,pawrad,pawtab)
     onsite_l = onsite_l + onsite_l_k

  end do

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(onsite_l,spaceComm,ierr)
  end if

  !---------clean up memory-------------------

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)

end subroutine make_onsite_l
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_onsite_bm
!! NAME
!! make_onsite_bm
!!
!! FUNCTION
!! Compute A_0.A_N onsite term for magnetic field + nuclear magnetic dipole moment
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

subroutine make_onsite_bm(atindx1,cprj,dtset,idir,mcprj,mpi_enreg,nband_k,onsite_bm,&
     & pawang,pawrad,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: idir,mcprj,nband_k
  complex(dpc),intent(out) :: onsite_bm
  type(MPI_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx1(dtset%natom)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: gint,iatom,icprj,ierr,ikpt,ikpt_loc,il,im,ilmn,isppol,itypat
  integer :: jl,jm,jlmn,klmn,klm,kln,lpmp,me,mesh_size,my_nspinor,ncpgr,nn,nproc,spaceComm
  real(dp) :: bm1,bm2,d00,d20,d22,dij,intg,scale_conversion
  complex(dpc) :: cpb,cpk

  !arrays
  integer :: nattyp_dum(dtset%ntypat)
  integer,allocatable :: dimlmn(:)
  real(dp),allocatable :: ff(:)
  type(pawcprj_type),allocatable :: cprj_k(:,:)

  ! ***********************************************************************

  ! this term can only be non-zero if some nucdipmom is nonzero
  scale_conversion = half*FineStructureConstant2
  d00 = sqrt(4.0*pi)/3.0
  dij = sqrt(4.0*pi/15.0)
  d20 = sqrt(16.0*pi/5.0)/6.0
  d22 = sqrt(16.0*pi/15.0)/2.0
  onsite_bm = czero

  ! TODO: generalize to nsppol > 1
  isppol = 1
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

  !Init MPI
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  me = mpi_enreg%me_kpt

  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp_dum,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,nband_k))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)

  ! loop over kpts on each processor

  ikpt_loc = 0
  ! loop over all the kpts
  do ikpt = 1, dtset%nkpt

     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

     ikpt_loc = ikpt_loc + 1
     icprj= (ikpt_loc - 1)*nband_k
     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt_loc,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

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
              do nn = 1, nband_k
                 cpb=cmplx(cprj_k(iatom,nn)%cp(1,ilmn),cprj_k(iatom,nn)%cp(2,ilmn),KIND=dpc)
                 cpk=cmplx(cprj_k(iatom,nn)%cp(1,jlmn),cprj_k(iatom,nn)%cp(2,jlmn),KIND=dpc)
                 onsite_bm=onsite_bm+conjg(cpb)*(bm1-bm2)*cpk
              end do ! end loop over nn
           end do ! end loop over ilmn
        end do ! end loop over jlmn
        ABI_DEALLOCATE(ff)
     end do ! end loop over atoms
  end do ! end loop over local k points

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(onsite_bm,spaceComm,ierr)
  end if

  !---------clean up memory-------------------

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)

end subroutine make_onsite_bm
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_eeig
!! NAME
!! make_eeig
!!
!! FUNCTION
!! Compute the energy eigenvalues at each k point
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
!! eeig
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

subroutine make_eeig(atindx1,cg,cprj,dtset,eeig,gmet,gprimd,mcg,mcprj,mpi_enreg,&
     & nattyp,nband_k,nfftf,npwarr,&
     & paw_ij,pawfgr,pawtab,psps,rmet,rprimd,&
     & vectornd,vhartr,vpsp,vxc,with_vectornd,xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,nband_k,nfftf,with_vectornd
 type(dataset_type),intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

 !arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfftf),vpsp(nfftf),vxc(nfftf,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: vectornd(with_vectornd*nfftf,3)
 real(dp),intent(out) :: eeig(nband_k,dtset%nkpt)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

 !Local variables -------------------------
 !scalars
 integer :: cpopt,dimffnl,eeig_size,exchn2n3d
 integer :: ierr,icg,icprj,ider,idir,ikg,ikg1,ikpt,ilm,isppol,istwf_k
 integer :: me,my_nspinor,ncpgr,ndat,ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,nkpg,nn
 integer :: nproc,npw_k,npw_k_,prtvol,sij_opt,spaceComm,tim_getghc,type_calc
 logical :: has_vectornd
 real(dp) :: ecut_eff,lambda
 type(gs_hamiltonian_type) :: gs_hamk

 !arrays
 integer,allocatable :: dimlmn(:),kg_k(:,:)
 real(dp) :: kpoint(3),lambdarr(1),rhodum(1)
 real(dp),allocatable :: buffer1(:),buffer2(:),cgrvtrial(:,:),cwavef(:,:)
 real(dp),allocatable :: ffnl_k(:,:,:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp),allocatable :: kinpw(:),kpg_k(:,:)
 real(dp),allocatable :: ph3d(:,:,:),vectornd_pac(:,:,:,:,:),vlocal(:,:,:,:),vtrial(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:),cwaveprj(:,:)


 !-----------------------------------------------------------------------

 !Init MPI
 spaceComm=mpi_enreg%comm_cell
 nproc=xmpi_comm_size(spaceComm)
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 me = mpi_enreg%me_kpt

 ! TODO: generalize to nsppol > 1
 isppol = 1
 ngfft1=dtset%ngfft(1) ; ngfft2=dtset%ngfft(2) ; ngfft3=dtset%ngfft(3)
 ngfft4=dtset%ngfft(4) ; ngfft5=dtset%ngfft(5) ; ngfft6=dtset%ngfft(6)
 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0

 has_vectornd = (with_vectornd .EQ. 1)

 ! input parameters for calls to getghc at ikpt
 cpopt = 4 ! was 4
 ndat = 1
 prtvol = 0
 sij_opt = 0
 tim_getghc = 0
 ! getghc: type_calc 0 means kinetic, local, nonlocal
 type_calc = 0
 lambda = zero; lambdarr(1) = zero

 !==== Initialize most of the Hamiltonian ====
 !Allocate all arrays and initialize quantities that do not depend on k and spin.
 !gs_hamk is the normal hamiltonian at k, needed for computing E_nk
 call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
      & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,nucdipmom=dtset%nucdipmom,&
      & paw_ij=paw_ij)

 !---------construct local potential------------------
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 ! nspden=1 is essentially hard-coded in the following line
 vtrial(1:nfftf,1)=vhartr(1:nfftf)+vxc(1:nfftf,1)+vpsp(1:nfftf)

 ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
 call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)

 ABI_ALLOCATE(vlocal,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc))
 call fftpac(isppol,mpi_enreg,dtset%nspden,&
      & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vlocal,2)

 ABI_DEALLOCATE(cgrvtrial)
 ABI_DEALLOCATE(vtrial)

 ! if vectornd is present, set it up for addition to gs_hamk similarly to how it's done for
 ! vtrial. Note that it must be done for the three directions. Also, the following
 ! code assumes explicitly and implicitly that nvloc = 1. This should eventually be generalized.
 if(has_vectornd) then
    ABI_ALLOCATE(vectornd_pac,(ngfft4,ngfft5,ngfft6,gs_hamk%nvloc,3))
    ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
    do idir = 1, 3
       call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vectornd(:,idir))
       call fftpac(isppol,mpi_enreg,dtset%nspden,&
            & ngfft1,ngfft2,ngfft3,ngfft4,ngfft5,ngfft6,dtset%ngfft,cgrvtrial,vectornd_pac(:,:,:,1,idir),2)
    end do
    ABI_DEALLOCATE(cgrvtrial)
 end if

 ! add vlocal
 call gs_hamk%load_spin(isppol,vlocal=vlocal,with_nonlocal=.true.)

 ! add vectornd if available
 if(has_vectornd) then
    call gs_hamk%load_spin(isppol,vectornd=vectornd_pac)
 end if

 ncpgr = cprj(1,1)%ncpgr
 ABI_ALLOCATE(dimlmn,(dtset%natom))
 call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')

 ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,nband_k))
 call pawcprj_alloc(cprj_k,ncpgr,dimlmn)
 ABI_DATATYPE_ALLOCATE(cwaveprj,(dtset%natom,1))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

 ABI_ALLOCATE(kg_k,(3,dtset%mpw))
 ABI_ALLOCATE(kinpw,(dtset%mpw))

 eeig(:,:) = zero
 icg = 0
 ikg = 0
 icprj = 0
 do ikpt = 1, dtset%nkpt

    ! if the current kpt is not on the current processor, cycle
    if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

    kpoint(:)=dtset%kptns(:,ikpt)
    npw_k = npwarr(ikpt)

    ! Build basis sphere of plane waves for the k-point
    kg_k(:,:) = 0
    call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg_k,kpoint,1,mpi_enreg,dtset%mpw,npw_k_)
    if (npw_k .NE. npw_k_) then
       write(std_out,'(a)')'JWZ debug mpi_eeig npw_k inconsistency'
    end if

    ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang))
    ABI_ALLOCATE(ylmgr_k,(npw_k,3,psps%mpsang*psps%mpsang*psps%useylm))
    do ilm=1,psps%mpsang*psps%mpsang
       ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       ylmgr_k(1:npw_k,1:3,ilm)=ylmgr(1+ikg:npw_k+ikg,1:3,ilm)
    end do

    ! Compute kinetic energy at kpt
    kinpw(:) = zero
    call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,kinpw,kpoint,npw_k,0,0)

    nkpg = 3
    ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
    call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

    ! Compute nonlocal form factors ffnl at all (k+G):
    ider=0 ! want ffnl and 1st derivative
    idir=4 ! d ffnl/ dk_red in all 3 directions
    dimffnl=1 ! 1 + number of derivatives
    ABI_ALLOCATE(ffnl_k,(npw_k,dimffnl,psps%lmnmax,dtset%ntypat))
    call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_k,psps%ffspl,&
         &         gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
         &         psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
         &         npw_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
         &         psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

    !  - Compute 3D phase factors
    !  - Prepare various tabs in case of band-FFT parallelism
    !  - Load k-dependent quantities in the Hamiltonian
    ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
    call gs_hamk%load_k(kpt_k=kpoint(:),istwf_k=istwf_k,npw_k=npw_k,&
         &         kinpw_k=kinpw,kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl_k,ph3d_k=ph3d,&
         &         compute_ph3d=.TRUE.,compute_gbound=.TRUE.)


    call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
         &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

    ! apply gs_hamk to wavefunctions at k to compute E_nk eigenvalues
    ABI_ALLOCATE(cwavef,(2,npw_k))
    ABI_ALLOCATE(ghc,(2,npw_k))
    ABI_ALLOCATE(gsc,(2,npw_k))
    ABI_ALLOCATE(gvnlc,(2,npw_k))
    do nn = 1, nband_k
       cwavef(1:2,1:npw_k) = cg(1:2,icg+(nn-1)*npw_k+1:icg+nn*npw_k)
       call pawcprj_get(atindx1,cwaveprj,cprj_k,dtset%natom,nn,0,ikpt,0,isppol,dtset%mband,&
            &           dtset%mkmem,dtset%natom,1,nband_k,my_nspinor,dtset%nsppol,0)
       call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,lambda,mpi_enreg,ndat,&
            &           prtvol,sij_opt,tim_getghc,type_calc)
       eeig(nn,ikpt) = DOT_PRODUCT(cwavef(1,1:npw_k),ghc(1,1:npw_k)) &
            &           + DOT_PRODUCT(cwavef(2,1:npw_k),ghc(2,1:npw_k))
    end do

    icg = icg + npw_k*nband_k
    ikg = ikg + npw_k
    icprj = icprj + nband_k

    ABI_DEALLOCATE(ylm_k)
    ABI_DEALLOCATE(ylmgr_k)
    ABI_DEALLOCATE(kpg_k)
    ABI_DEALLOCATE(ffnl_k)
    ABI_DEALLOCATE(ph3d)
    ABI_DEALLOCATE(cwavef)
    ABI_DEALLOCATE(ghc)
    ABI_DEALLOCATE(gsc)
    ABI_DEALLOCATE(gvnlc)

 end do ! end loop over kpts on current processor

 !  MPI communicate stuff between everyone
 if (nproc>1) then
    eeig_size = size(eeig)
    ABI_ALLOCATE(buffer1,(eeig_size))
    ABI_ALLOCATE(buffer2,(eeig_size))
    buffer1(1:eeig_size) = reshape(eeig,(/eeig_size/))
    call xmpi_sum(buffer1,buffer2,eeig_size,spaceComm,ierr)
    eeig(1:nband_k,1:dtset%nkpt)=reshape(buffer2,(/nband_k,dtset%nkpt/))
    ABI_DEALLOCATE(buffer1)
    ABI_DEALLOCATE(buffer2)
 end if

 call gs_hamk%free()
 ABI_DEALLOCATE(vlocal)
 if(has_vectornd) then
    ABI_DEALLOCATE(vectornd_pac)
 end if

 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kinpw)

 ABI_DEALLOCATE(dimlmn)
 call pawcprj_free(cprj_k)
 ABI_DATATYPE_DEALLOCATE(cprj_k)
 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

end subroutine make_eeig
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_S1trace
!! NAME
!! make_S1trace
!!
!! FUNCTION
!! Compute Trace[\rho_0 S^{(1)} \rho_0] in orbital magnetism context
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

subroutine make_S1trace(adir,atindx1,cprj,dtset,eeig,&
     & mcprj,mpi_enreg,nattyp,nband_k,pawtab,S1trace)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,mcprj,nband_k
  complex(dpc),intent(out) :: S1trace
  type(MPI_type), intent(inout) :: mpi_enreg
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
  real(dp),intent(in) :: eeig(nband_k,dtset%nkpt)
  type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,icprj,epsabg,gdir,iatom,ierr,ikpt,ilmn,isppol,itypat
  integer :: jlmn,klmn,me,my_nspinor,ncpgr,nn,nproc,spaceComm
  real(dp) :: ENK
  complex(dpc) :: cpb,cpk

  !arrays
  integer,allocatable :: dimlmn(:)
  type(pawcprj_type),allocatable :: cprj_k(:,:)

!----------------------------------------------------------------

  !Init MPI
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  me = mpi_enreg%me_kpt

  isppol = 1
  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,nband_k))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)

  S1trace = czero
  icprj = 0
  do ikpt = 1, dtset%nkpt

     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

     do epsabg = 1, -1, -2

        if (epsabg .EQ. 1) then
           bdir = modulo(adir,3)+1
           gdir = modulo(adir+1,3)+1
        else
           bdir = modulo(adir+1,3)+1
           gdir = modulo(adir,3)+1
        end if

        do nn = 1, nband_k
           ENK = eeig(nn,ikpt)
           do iatom=1,dtset%natom
              itypat=dtset%typat(iatom)
              do ilmn=1,pawtab(itypat)%lmn_size
                 do jlmn=1,pawtab(itypat)%lmn_size
                    klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
                    cpb=cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
                    cpk=cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)
                    S1trace=S1trace+half*j_dpc*epsabg*ENK*conjg(cpb)*pawtab(itypat)%sij(klmn)*cpk
                 end do ! end loop over jlmn
              end do ! end loop over ilmn
           end do ! end loop over atoms
        end do ! end loop over bands
     end do ! end loop over epsabg

     icprj = icprj + nband_k

  end do ! end loop over kpt

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(S1trace,spaceComm,ierr)
  end if

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)

end subroutine make_S1trace
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_rhorij1
!! NAME
!! make_rhorij1
!!
!! FUNCTION
!! Compute Trace[\rho_0 \rho_Rij(1) ] in orbital magnetism context
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

subroutine make_rhorij1(adir,atindx1,cprj,dtset,mcprj,mpi_enreg,&
     & nattyp,nband_k,paw_ij,pawtab,rhorij1)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: adir,mcprj,nband_k
  complex(dpc),intent(out) :: rhorij1
  type(MPI_type), intent(inout) :: mpi_enreg
  type(dataset_type),intent(in) :: dtset

  !arrays
  integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat)
  type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj)
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: bdir,icprj,epsabg,gdir,iatom,ierr,ikpt,ilmn,isppol,itypat
  integer :: jlmn,klmn,me,my_nspinor,ncpgr,nn,nproc,spaceComm
  complex(dpc) :: cpb,cdij,cpk

  !arrays
  integer,allocatable :: dimlmn(:)
  type(pawcprj_type),allocatable :: cprj_k(:,:)

!----------------------------------------------------------------

  !Init MPI
  spaceComm=mpi_enreg%comm_cell
  nproc=xmpi_comm_size(spaceComm)
  my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  me = mpi_enreg%me_kpt

  isppol = 1
  ncpgr = cprj(1,1)%ncpgr
  ABI_ALLOCATE(dimlmn,(dtset%natom))
  call pawcprj_getdim(dimlmn,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'R')
  ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natom,nband_k))
  call pawcprj_alloc(cprj_k,ncpgr,dimlmn)

  rhorij1 = czero
  icprj = 0
  do ikpt = 1, dtset%nkpt

     ! if the current kpt is not on the current processor, cycle
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle

     call pawcprj_get(atindx1,cprj_k,cprj,dtset%natom,1,icprj,ikpt,0,isppol,dtset%mband,&
          &       dtset%mkmem,dtset%natom,nband_k,nband_k,my_nspinor,dtset%nsppol,0)

     do epsabg = 1, -1, -2

        if (epsabg .EQ. 1) then
           bdir = modulo(adir,3)+1
           gdir = modulo(adir+1,3)+1
        else
           bdir = modulo(adir+1,3)+1
           gdir = modulo(adir,3)+1
        end if

        do nn = 1, nband_k
           do iatom=1,dtset%natom
              itypat=dtset%typat(iatom)
              do ilmn=1,pawtab(itypat)%lmn_size
                 do jlmn=1,pawtab(itypat)%lmn_size
                    klmn=max(jlmn,ilmn)*(max(jlmn,ilmn)-1)/2 + min(jlmn,ilmn)
                    cpb=cmplx(cprj_k(iatom,nn)%dcp(1,bdir,ilmn),cprj_k(iatom,nn)%dcp(2,bdir,ilmn),KIND=dpc)
                    cpk=cmplx(cprj_k(iatom,nn)%dcp(1,gdir,jlmn),cprj_k(iatom,nn)%dcp(2,gdir,jlmn),KIND=dpc)
                    if (paw_ij(iatom)%cplex_dij .EQ. 2) then
                       cdij=cmplx(paw_ij(iatom)%dij(2*klmn-1,1),paw_ij(iatom)%dij(2*klmn,1),KIND=dpc)
                       if (jlmn .GT. ilmn) cdij=conjg(cdij)
                    else
                       cdij=cmplx(paw_ij(iatom)%dij(klmn,1),zero,KIND=dpc)
                    end if
                    rhorij1=rhorij1+half*j_dpc*epsabg*conjg(cpb)*cdij*cpk
                 end do ! end loop over jlmn
              end do ! end loop over ilmn
           end do ! end loop over atoms
        end do ! end loop over bands
     end do ! end loop over epsabg

     icprj = icprj + nband_k

  end do ! end loop over kpt

  ! ---- parallel communication
  if(nproc > 1) then
     call xmpi_sum(rhorij1,spaceComm,ierr)
  end if

  ABI_DEALLOCATE(dimlmn)
  call pawcprj_free(cprj_k)
  ABI_DATATYPE_DEALLOCATE(cprj_k)

end subroutine make_rhorij1
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/cpg_dij_cpb
!! NAME
!! cpg_dij_cpb
!!
!! FUNCTION
!! Compute <u_kg|p_i>dij<p_j|u_kb> energy contribution
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
!! cgdijcb
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

subroutine cpg_dij_cpb(cgdijcb,cprj_kb,cprj_kg,dtset,nb,ng,nspinor,paw_ij,pawtab)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: nb,ng,nspinor
  complex(dpc),intent(out) :: cgdijcb
  type(dataset_type),intent(in) :: dtset

  !arrays
  type(pawcprj_type),intent(in) :: cprj_kb(dtset%natom,nspinor*dtset%mband)
  type(pawcprj_type),intent(in) :: cprj_kg(dtset%natom,nspinor*dtset%mband)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

  !Local variables -------------------------
  !scalars
  integer :: iatom,ilmn,jlmn,klmn,itypat
  complex(dpc) :: cdij,cpg,cpb

!-----------------------------------------------------------------------

  cgdijcb = czero
  do iatom = 1, dtset%natom
     itypat = dtset%typat(iatom)
     do ilmn = 1, pawtab(itypat)%lmn_size
        cpg=cmplx(cprj_kg(iatom,ng)%cp(1,ilmn),cprj_kg(iatom,ng)%cp(2,ilmn),KIND=dpc)
        do jlmn = 1, pawtab(itypat)%lmn_size
           cpb=cmplx(cprj_kb(iatom,nb)%cp(1,jlmn),cprj_kb(iatom,nb)%cp(2,jlmn),KIND=dpc)
           if (jlmn .LE. ilmn) then
              klmn = (ilmn-1)*ilmn/2 + jlmn
           else
              klmn = (jlmn-1)*jlmn/2 + ilmn
           end if
           if (paw_ij(iatom)%cplex_dij .EQ. 2) then
              cdij=cmplx(paw_ij(iatom)%dij(2*klmn-1,1),paw_ij(iatom)%dij(2*klmn,1),KIND=dpc)
              if (jlmn .GT. ilmn) cdij=conjg(cdij)
           else
              cdij=cmplx(paw_ij(iatom)%dij(klmn,1),zero,KIND=dpc)
           end if
           cgdijcb = cgdijcb + conjg(cpg)*cdij*cpb
        end do
     end do
  end do

end subroutine cpg_dij_cpb
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
     & mcg,mcprj,mpi_enreg,nattyp,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,&
     & pwind,pwind_alloc,rprimd,symrec,usecprj,vectornd,&
     & vhartr,vpsp,vxc,with_vectornd,xred,ylm,ylmgr)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: mcg,mcprj,nfftf,pwind_alloc,usecprj,with_vectornd
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
 real(dp),intent(inout) :: vectornd(with_vectornd*nfftf,3)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawcprj_type),intent(in) ::  cprj(dtset%natom,mcprj*usecprj)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

 !Local variables -------------------------
 !scalars
 integer :: adir,bdx,gdxstor
 integer :: isppol,istwf_k,my_nspinor
 integer :: nband_k,ncpgr,ncpgrb
 real(dp) :: chernnorm,magnorm,ucvol,finish_time,start_time
 complex(dpc) :: CCI_dir,VVI_dir,VVII_dir,VVIII_dir
 complex(dpc) :: CCIV_dir,dpds_dir,onsite_bm_dir,onsite_l_dir,rhorij1_dir,s1trace_dir
 character(len=500) :: message

 !arrays
 integer,allocatable :: dimlmn(:)
 real(dp) :: CCI(2,3),CCIV(2,3),cnum_dpdp(2,3),cnum_dpds(2,3),gmet(3,3),gprimd(3,3)
 real(dp) :: duqduchern(2,3),duqdumag(2,3),udsqduchern(2,3),udsqdumag(2,3)
 real(dp) :: onsite_bm(2,3),onsite_l(2,3),orbmagvec(2,3),rhorij1(2,3)
 real(dp) :: rmet(3,3),s1trace(2,3),VVI(2,3),VVI_dsdu(2,3)
 real(dp) :: VVII(2,3),VVII_dusdu(2,3),VVII_udsdsu(2,3),VVIII(2,3),VVIII_dsdu(2,3)
 real(dp),allocatable :: eeig(:,:)

 ! ***********************************************************************
 ! my_nspinor=max(1,dtorbmag%nspinor/mpi_enreg%nproc_spinor)

 ncpgr = cprj(1,1)%ncpgr

 ! TODO: generalize to nsppol > 1
 isppol = 1
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 nband_k = dtorbmag%mband_occ
 istwf_k = 1

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ABI_ALLOCATE(eeig,(nband_k,dtset%nkpt))
 eeig(:,:) = zero
 if (dtset%orbmag .GT. 1) then
    call make_eeig(atindx1,cg,cprj,dtset,eeig,gmet,gprimd,mcg,mcprj,mpi_enreg,nattyp,nband_k,nfftf,npwarr,&
         & paw_ij,pawfgr,pawtab,psps,rmet,rprimd,&
         & vectornd,vhartr,vpsp,vxc,with_vectornd,xred,ylm,ylmgr)
 end if

 ! compute i*\epsilon_{abg}\sum_n <du|Q|du> with and without E_nk weights (needed respectively
 ! by Chern number and by magnetization)
 call duqdu(atindx1,cg,cprj,dtorbmag,dtset,duqduchern,duqdumag,eeig,gmet,gprimd,mcg,mcprj,mpi_enreg,&
      & nband_k,npwarr,pawang,pawrad,pawtab,psps,pwind,pwind_alloc,xred)

 ! compute i*\epsilon_{abg}\sum_n <u|dS Q|du> with and without E_nk weights (needed respectively
 ! by Chern number and by magnetization)
 call udsqdu(atindx1,cg,cprj,dtorbmag,dtset,eeig,gmet,gprimd,kg,&
     & mcg,mcprj,mpi_enreg,nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,psps,&
     pwind,pwind_alloc,rmet,rprimd,udsqduchern,udsqdumag,xred,ylm,ylmgr)

 ! ! call chern number routines if necessary
 if ( (dtset%orbmag .EQ. 1) .OR. (dtset%orbmag .EQ. 3) ) then

    ! convert from crystallographic to cartesian
    duqduchern(1,1:3) = ucvol*MATMUL(gprimd,duqduchern(1,1:3))
    duqduchern(2,1:3) = ucvol*MATMUL(gprimd,duqduchern(2,1:3))

    ! factor of two in numerator is band occupations
    ! factor of two_pi in denominator is in definition of Chern number
    ! factor of ucvol*fnkpt is discrete representation of integral over BZ
    chernnorm = two/(two_pi*ucvol*dtorbmag%fnkpt)
    
    duqduchern(1:2,1:3) = duqduchern(1:2,1:3)*chernnorm

    udsqduchern(1,1:3) = ucvol*MATMUL(gprimd,udsqduchern(1,1:3))
    udsqduchern(2,1:3) = ucvol*MATMUL(gprimd,udsqduchern(2,1:3))
    udsqduchern(1:2,1:3) = udsqduchern(1:2,1:3)*chernnorm

    dtorbmag%chern(1:2,1:3) = duqduchern(1:2,1:3)+udsqduchern(1:2,1:3)
    
    call output_orbmag(2,dtorbmag%chern)
    
 end if

 ! continue with computation of orbital magnetization if necessary
 if ( dtset%orbmag .GT. 1 ) then

    do adir = 1, 3

       call make_onsite_l(atindx1,cprj,dtset,adir,mcprj,mpi_enreg,nband_k,onsite_l_dir,pawrad,pawtab)
       onsite_l(1,adir) = real(onsite_l_dir)
       onsite_l(2,adir) = aimag(onsite_l_dir)

       call make_S1trace(adir,atindx1,cprj,dtset,eeig,mcprj,mpi_enreg,nattyp,nband_k,pawtab,s1trace_dir)
       s1trace(1,adir) = real(s1trace_dir)
       s1trace(2,adir) = aimag(s1trace_dir)

       call make_rhorij1(adir,atindx1,cprj,dtset,mcprj,mpi_enreg,nattyp,nband_k,paw_ij,pawtab,rhorij1_dir)
       rhorij1(1,adir) = real(rhorij1_dir)
       rhorij1(2,adir) = aimag(rhorij1_dir)

       if (any(abs(dtset%nucdipmom)>tol8)) then
          call make_onsite_bm(atindx1,cprj,dtset,adir,mcprj,mpi_enreg,nband_k,onsite_bm_dir,&
               & pawang,pawrad,pawtab)
          onsite_bm(1,adir) = real(onsite_bm_dir)
          onsite_bm(2,adir) = aimag(onsite_bm_dir)
       else
          onsite_bm(:,adir) = zero
       end if

    end do ! end loop over adir

    CCI=zero
    call duqhqdu(atindx1,cg,CCI,cprj,dtorbmag,dtset,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
         & nattyp,nband_k,nfftf,npwarr,paw_ij,pawang,pawfgr,pawrad,pawtab,psps,pwind,pwind_alloc,&
         & rmet,rprimd,ucvol,vectornd,vhartr,vpsp,vxc,with_vectornd,xred)
    CCI=-half*CCI

    ! call udsqdu(atindx1,cg,VVI_dsdu,cprj,dtorbmag,dtset,gmet,gprimd,&
    !      & kg,mcg,mcprj,mpi_enreg,&
    !      & nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,&
    !      & psps,pwind,pwind_alloc,&
    !      & rmet,rprimd,xred,ylm,ylmgr,energies=eeig,swap_option=.TRUE.)
    VVI=half*udsqdumag

    call udsdsu(atindx1,cg,VVII_udsdsu,cprj,dtorbmag,dtset,eeig,gmet,gprimd,kg,mcg,mcprj,mpi_enreg,&
         & nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,psps,rmet,rprimd,xred,ylm,ylmgr)

    VVII=half*(duqdumag - VVII_udsdsu) ! check this

    ! call udsqdu(atindx1,cg,VVIII_dsdu,cprj,dtorbmag,dtset,gmet,gprimd,&
    !      & kg,mcg,mcprj,mpi_enreg,&
    !      & nband_k,npwarr,paw_ij,pawang,pawrad,pawtab,&
    !      & psps,pwind,pwind_alloc,&
    !      & rmet,rprimd,xred,ylm,ylmgr,energies=eeig)
    VVIII=half*udsqdumag
    
    ! convert terms to cartesian coordinates as needed
    ! note that terms like <dv/dk| x |dw/dk> computed in reduced coords,
    ! become ucvol*gprimd*<dv/dk| x |dw/dk> when expressed in cartesian coords
    ! onsite_l and onsite_bm are already cartesian

    s1trace(1,1:3) = ucvol*MATMUL(gprimd,s1trace(1,1:3))
    s1trace(2,1:3) = ucvol*MATMUL(gprimd,s1trace(2,1:3))

    rhorij1(1,1:3) = ucvol*MATMUL(gprimd,rhorij1(1,1:3))
    rhorij1(2,1:3) = ucvol*MATMUL(gprimd,rhorij1(2,1:3))
    
    CCI(1,1:3) = ucvol*MATMUL(gprimd,CCI(1,1:3))
    CCI(2,1:3) = ucvol*MATMUL(gprimd,CCI(2,1:3))
    
    VVII(1,1:3) = ucvol*MATMUL(gprimd,VVII(1,1:3))
    VVII(2,1:3) = ucvol*MATMUL(gprimd,VVII(2,1:3))
    
    VVI(1,1:3) = ucvol*MATMUL(gprimd,VVI(1,1:3))
    VVI(2,1:3) = ucvol*MATMUL(gprimd,VVI(2,1:3))

    VVIII(1,1:3) = ucvol*MATMUL(gprimd,VVIII(1,1:3))
    VVIII(2,1:3) = ucvol*MATMUL(gprimd,VVIII(2,1:3))

    ! scale for integration over Brillouin zone
    ! pre factor is occ/ucvol*N_k
    ! factor of 2 in numerator is the band occupation (two electrons in normal insulator)
    ! converting integral over k space to a sum gives a factor of Omega_BZ/N_k or 1/ucvol*N_k
    magnorm = two/(ucvol*dtorbmag%fnkpt)
    onsite_l(1:2,1:3) = onsite_l(1:2,1:3)*magnorm
    onsite_bm(1:2,1:3) = onsite_bm(1:2,1:3)*magnorm
    s1trace(1:2,1:3) = s1trace(1:2,1:3)*magnorm
    rhorij1(1:2,1:3) = rhorij1(1:2,1:3)*magnorm
    CCI(1:2,1:3) = CCI(1:2,1:3)*magnorm
    VVII(1:2,1:3) = VVII(1:2,1:3)*magnorm
    VVI(1:2,1:3) = VVI(1:2,1:3)*magnorm
    VVIII(1:2,1:3) = VVIII(1:2,1:3)*magnorm

    write(std_out,'(a,3es16.8)')' JWZ debug onsite_l ',onsite_l(1,1),onsite_l(1,2),onsite_l(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug onsite_bm ',onsite_bm(1,1),onsite_bm(1,2),onsite_bm(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug s1trace ',s1trace(1,1),s1trace(1,2),s1trace(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug rhorij1 ',rhorij1(1,1),rhorij1(1,2),rhorij1(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug CCI ',CCI(1,1),CCI(1,2),CCI(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug VVII ',VVII(1,1),VVII(1,2),VVII(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug VVI ',VVI(1,1),VVI(1,2),VVI(1,3)
    write(std_out,'(a,3es16.8)')' JWZ debug VVIII ',VVIII(1,1),VVIII(1,2),VVIII(1,3)

    ! accumulate in orbmagvec

    orbmagvec(1:2,1:3) = onsite_l(1:2,1:3)  &
         & + onsite_bm(1:2,1:3) &
         & - s1trace(1:2,1:3) &
         & + rhorij1(1:2,1:3) &
         & + VVII(1:2,1:3) &
         & + VVI(1:2,1:3) &
         & + VVIII(1:2,1:3) &
         & + CCI(1:2,1:3)       

    dtorbmag%orbmagvec(1:2,1:3) = two*orbmagvec(1:2,1:3)/(ucvol*dtorbmag%fnkpt) 

    ! orbmagvec(1:2,1:3) = CCI(1:2,1:3)       
    ! dtorbmag%orbmagvec(1:2,1:3) = orbmagvec(1:2,1:3)

    call output_orbmag(1,dtorbmag%orbmagvec)

 end if ! end computation and output of orbital magnetization

 if(allocated(eeig)) then
    ABI_DEALLOCATE(eeig)
 end if

end subroutine orbmag
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/output_orbmag
!! NAME
!! output_orbmag
!!
!! FUNCTION
!! This routine outputs orbmag and chern number
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! output_option 1 or 2 for orbital magnetization, chern number respectively
!! output_vector real(2,3) to be printed
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

subroutine output_orbmag(output_option,output_vector)

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: output_option

 !arrays
 real(dp),intent(in) :: output_vector(2,3)

 !Local variables -------------------------
 !scalars
 integer :: adir
 character(len=500) :: message

 ! ***********************************************************************

 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

 select case (output_option)
    ! case 1 is orbital magnetization
 case(1)
    write(message,'(a)')' Orbital magnetization '
    call wrtout(ab_out,message,'COLL')
    write(message,'(a,a)')'----Orbital magnetization is a real vector, given along Cartesian directions----',ch10
    call wrtout(ab_out,message,'COLL')
    do adir = 1, 3
       write(message,'(a,i4,a,2es16.8)')' Orb Mag(',adir,') : real, imag ',&
            &   output_vector(1,adir),output_vector(2,adir)
       call wrtout(ab_out,message,'COLL')
    end do
    ! case 2 is chern number
 case(2)
    write(message,'(a)')' Chern number C from orbital magnetization '
    call wrtout(ab_out,message,'COLL')
    write(message,'(a,a)')'----C is a real vector, given along Cartesian directions----',ch10
    call wrtout(ab_out,message,'COLL')
    do adir = 1, 3
       write(message,'(a,i4,a,2es16.8)')' C(',adir,') : real, imag ',&
            &   output_vector(1,adir),output_vector(2,adir)
       call wrtout(ab_out,message,'COLL')
    end do
 end select
    
 write(message,'(a,a,a)')ch10,'====================================================',ch10
 call wrtout(ab_out,message,'COLL')

end subroutine output_orbmag
!!***


end module m_orbmag
