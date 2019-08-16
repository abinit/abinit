!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_efield
!! NAME
!!  m_efield
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle electric fields
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (MJV)
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

module m_efield

 use defs_basis
 use m_abicore
 use m_errors

 use m_pawcprj, only : pawcprj_type, pawcprj_free

 implicit none

 private
!!***


!!****t* m_efield/efield_type
!! NAME
!! efield_type
!!
!! FUNCTION
!! First-principles calculations in a finite electric field
!!
!! SOURCE

 type, public :: efield_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
  integer :: berryopt            ! value of berryopt in use
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: has_epawf3          ! 2 if epawf3 computed, 1 if only allocated, zero else
  integer :: has_epaws3          ! 2 if epaws3 computed, 1 if only allocated, zero else
  integer :: has_expibi          ! 2 if expibi computed, 1 if only allocated, zero else
  integer :: has_rij             ! 2 if paw_rij computed, 1 if only allocated, zero else
  integer :: has_qijb            ! 2 if paw_qijb computed, 1 if only allocated, zero else
  integer :: lmax
  integer :: lmnmax
  integer :: lmn2max
  integer :: maxnstr             ! max number of strings along idir=1,2,3
  integer :: maxnkstr            ! max number of k-points per string
  integer :: mkmem_max           ! max of mkmem
  integer :: natom               ! number of atoms in unit cell
  integer :: my_natom            ! number of atoms treated by current proc
  integer :: mband_occ           ! max number of occupied bands (over spin)
                                 ! this number must be the same for every k
  integer :: nspinor             ! nspinor input from data set
  integer :: nsym
  integer :: usecprj             ! 1 if efield%cprj allocated (see below), 0 else
  integer :: usepaw              ! 1 if a PAW calculation, 0 else

! Integer arrays

  integer :: nstr(3)             ! nstr(idir) = number of strings along idir
  integer :: nkstr(3)            ! nkstr(idir) = number of k-points per string

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1
                                 !                         1 if nsppol = 2

! Real(dp) arrays
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
!! if berryopt=4,6,7, dtset%red_efieldbar and dtset%red_dfield will be initialized in 67_common/initberry.F90
  real(dp) :: efield_dot(3)      ! reciprocal lattice coordinates of the
                                 ! electric field
  real(dp) :: red_ptot1(3)       ! reduced total polarization
  real(dp) :: efield2(3)         ! unreduced electric field, only used when berryopt == 14 in order to save real electric field for print out.

  real(dp) :: gmet_str(2,2,3)    ! gmet_str(:,:,idir) is the metric of the metric of
                                 ! the space of strings of direction idir
! Integer pointers
  integer, allocatable :: atom_indsym(:,:,:) ! atom_indsym(4,nsym,natom)
                                         ! this is data on how the symmetries map the atoms in the cell
                                         ! see symatm.F90 for full description
  integer, allocatable :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, allocatable :: cgqindex(:,:,:) ! cgqindex(3,6,nkpt*nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cgq and pwnsfacq
                                      ! arrays
                                      ! (see vtorho.f and initberry.f)
  integer, allocatable :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the cprj in the cprj array (used only
                                      ! for PAW calculations)
  integer, allocatable :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, allocatable :: idxkstr(:,:,:)  ! idxkstr(maxnkstr,maxnstr,3)
                                      ! idxkstr(ikstr,istr,idir) index (ikpt) of
                                      ! k-point ikstr on string istr along idir
  integer, allocatable :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, allocatable :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtefield%fnkpt,1:6)
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
  integer, allocatable :: nneigh(:)          ! nneigh(nkpt)
                                         ! for each k-point, nneigh stores
                                         ! the number of its nearest neighbours
                                         ! that are not related by symmetry
  integer, allocatable :: sflag(:,:,:,:)  ! sflag(mband_occ,nkpt*nsppol,2,3)
                                      ! sflag = 0 : compute the whole row of
                                      !             smat
                                      ! sflag = 1 : the row is up to date

  integer, allocatable :: str_neigh(:,:,:)
  integer, allocatable :: strg_neigh(:,:,:,:)
! str_neigh(ineigh, istr, idir) is the index ineigh-th neighbour of the istr-th string in
! the direction idir
! str_neigh(ineigh, istr, :, idir) is a 2-dimensional vector which coordinates are 0 or 1,
! useful only if the k-point mesh isn't a full mesh - if it's a single point, a line or a plane.


! Real(dp) allocatables

! the coordinates of the ineigh-th neighbour of the istr-th string in the direction idir are :
! coord_str(:,str_neigh(ineigh,istr,idir),idir) + real(str_neigh(ineigh, istr, :, idir),dp)
  real(dp),allocatable :: coord_str(:,:,:)
! coord_str(1:2,istr,idir) are the coordinate of the istr-th string in the direction idir.

  real(dp),allocatable :: epawf3(:,:,:)
! epawf3(natom,3,3) ! F3-type force term (derivatives of projectors with respect to ion posiion)
! that arises in force for finite electric field with PAW
! epawf3(iatom,idir,fdir) is derivative of polarization component idir with respect to iatom
! displaced in direction fdir
! see equation 32 of Torrent et al. CMS 42, 337 (2008)

  real(dp),allocatable :: epaws3(:,:,:)
! epaws3(natom,3,6) ! F3-type stress term (derivatives of projectors with respect to strain)
! that arises in stress for finite electric field with PAW
! epaws3(iatom,idir,strain) is derivative of polarization component idir with respect to strain
! component for atom iatom (note that these are on-site terms)
! see equation D.7 of Torrent et al. CMS 42, 337 (2008)

  real(dp), allocatable :: expibi(:,:,:)
! expibi(2,my_natom,3)
! used for PAW field calculations (distributed over atomic sites)
! stores the on-site phase factors arising from
! $\langle\phi_{i,k}|\phi_{j,k+\sigma_k k_k}\rangle$
! where $\sigma = \pm 1$. These overlaps arise in various Berry
! phase calculations of electric and magnetic polarization. The on-site
! phase factor is $\exp[-i\sigma_k k_k)\cdot I]$ where
! $I$ is the nuclear position. Only the following
! are computed and saved, in the given order:
! 1)    -k_1
! 2)    -k_2
! 3)    -k_3

  real(dp), allocatable :: fkptns(:,:)       ! fkptns(3,1:dtefield%fnkpt)
                                         ! k-points in FBZ

  real(dp), allocatable :: qijb_kk(:,:,:,:)
! qijb_kk(2,lmn2max,natom,3)
! on-site part of <u_nk|u_mk+b> matrix elements, relevant for PAW only
! vector b described by idir (1,2,3), forward direction; value for
! reverse direction (ifor = 2 in berryphase_new and cgwf) obtained by
! complex conjugation

! pointer to on-site dipole moment
  real(dp),allocatable :: rij(:,:,:) ! rij(lmn2_size_max,ntypat,3)
 ! gives <r-R> at each atom in each of 3 directions
 ! these are used only in the PAW case with electric field

  real(dp), allocatable :: smat(:,:,:,:,:,:)
! smat(2,mband_occ,mband_occ,nkpt*nsppol,2,3)
! Overlap matrix for every k-point. In an electric field calculation,
! smat is updated at every iteration.

  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations. These are needed when the
   ! cprj's need to be computed in the full BZ, that is,
   ! in the PAW case with kptopt /= 3.

! pointer to cprj
   type(pawcprj_type),allocatable :: cprj(:,:)
! used with finite efield and PAW

 end type efield_type

 ! Bound methods:
 public :: destroy_efield
!!***

contains

!!****f* m_efield/destroy_efield
!! NAME
!!
!! FUNCTION
!!   deallocate fields in efield structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_efield(dtefield)

!Arguments ------------------------------------
!array
 type(efield_type),intent(inout) :: dtefield !vz_i

! ************************************************************************

! Integer pointers
  if(allocated(dtefield%atom_indsym))  then
    ABI_DEALLOCATE(dtefield%atom_indsym)
  end if
  if(allocated(dtefield%cgindex))  then
    ABI_DEALLOCATE(dtefield%cgindex)
  end if
  if(allocated(dtefield%cgqindex))  then
    ABI_DEALLOCATE(dtefield%cgqindex)
  end if
  if(allocated(dtefield%cprjindex))  then
    ABI_DEALLOCATE(dtefield%cprjindex)
  end if
  if(allocated(dtefield%fkgindex))  then
    ABI_DEALLOCATE(dtefield%fkgindex)
  end if
  if(allocated(dtefield%idxkstr))  then
    ABI_DEALLOCATE(dtefield%idxkstr)
  end if
  if(allocated(dtefield%ikpt_dk))  then
    ABI_DEALLOCATE(dtefield%ikpt_dk)
  end if
  if(allocated(dtefield%indkk_f2ibz))  then
    ABI_DEALLOCATE(dtefield%indkk_f2ibz)
  end if
  if(allocated(dtefield%i2fbz))  then
    ABI_DEALLOCATE(dtefield%i2fbz)
  end if
  if(allocated(dtefield%kgindex))  then
    ABI_DEALLOCATE(dtefield%kgindex)
  end if
  if(allocated(dtefield%lmn_size))  then
    ABI_DEALLOCATE(dtefield%lmn_size)
  end if
  if(allocated(dtefield%lmn2_size))  then
    ABI_DEALLOCATE(dtefield%lmn2_size)
  end if
  if(allocated(dtefield%nband_occ))  then
    ABI_DEALLOCATE(dtefield%nband_occ)
  end if
  if(allocated(dtefield%nneigh))  then
    ABI_DEALLOCATE(dtefield%nneigh)
  end if
  if(allocated(dtefield%sflag))  then
    ABI_DEALLOCATE(dtefield%sflag)
  end if
  if(allocated(dtefield%str_neigh))  then
    ABI_DEALLOCATE(dtefield%str_neigh)
  end if
  if(allocated(dtefield%strg_neigh))  then
    ABI_DEALLOCATE(dtefield%strg_neigh)
  end if

! Real(dp) pointers

  if(allocated(dtefield%coord_str))  then
    ABI_DEALLOCATE(dtefield%coord_str)
  end if
  if(allocated(dtefield%epawf3))  then
    ABI_DEALLOCATE(dtefield%epawf3)
  end if
  if(allocated(dtefield%epaws3))  then
    ABI_DEALLOCATE(dtefield%epaws3)
  end if
  if(allocated(dtefield%expibi))  then
    ABI_DEALLOCATE(dtefield%expibi)
  end if
  if(allocated(dtefield%fkptns))  then
    ABI_DEALLOCATE(dtefield%fkptns)
  end if
  if(allocated(dtefield%qijb_kk))  then
    ABI_DEALLOCATE(dtefield%qijb_kk)
  end if
  if(allocated(dtefield%rij))  then
    ABI_DEALLOCATE(dtefield%rij)
  end if
  if(allocated(dtefield%smat))  then
    ABI_DEALLOCATE(dtefield%smat)
  end if
  if(allocated(dtefield%zarot))  then
    ABI_DEALLOCATE(dtefield%zarot)
  end if

! pointer to cprj
  if(allocated(dtefield%cprj)) then
    call pawcprj_free(dtefield%cprj)
    ABI_DATATYPE_DEALLOCATE(dtefield%cprj)
  end if

end subroutine destroy_efield
!!***

end module m_efield
!!***
