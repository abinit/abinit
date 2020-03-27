!!****m* ABINIT/m_ab7_kpoints
!! NAME
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ab7_kpoints

  use defs_basis
  use m_ab7_symmetry
  use m_abicore

  use m_kpts,      only : getkgrid, testkgrid
  use m_spacepar,  only : irrzg

  implicit none

  private

  logical, private, parameter :: AB_DBG = .false.

  public :: kpoints_get_irreductible_zone

  public :: kpoints_get_mp_k_grid
  public :: kpoints_get_auto_k_grid

  public :: kpoints_binding_mp_k_1
  public :: kpoints_binding_mp_k_2
  public :: kpoints_binding_auto_k_1
  public :: kpoints_binding_auto_k_2

contains
!!***


!!****f* m_ab7_kpoints/kpoints_get_irreductible_zone
!! NAME
!!  kpoints_get_irreductible_zone
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE

  subroutine kpoints_get_irreductible_zone(irrzon, phnons, &
       & n1, n2, n3, nsppol, nspden, symid, errno)

    integer, intent(in)   :: symid
    integer, intent(in)   :: n1, n2, n3, nsppol, nspden
    integer, intent(out)  :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
    real(dp), intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))
    integer, intent(out)  :: errno

    type(symmetry_type), pointer  :: sym

    if (AB_DBG) write(std_err,*) "AB kpoints: call get irreductible zone."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    if (sym%withSpin /= nspden) then
       errno = AB7_ERROR_ARG
       return
    end if

    call irrzg(irrzon, nspden, nsppol, sym%nSym, n1, n2, n3, phnons, &
         & sym%symAfm, sym%sym, sym%transNon)
  end subroutine kpoints_get_irreductible_zone
!!***


!!****f* m_ab7_kpoints/kpoints_binding_mp_k_1
!! NAME
!!  kpoints_binding_mp_k_1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_ab7_kpoints
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE

  subroutine kpoints_binding_mp_k_1(symid, nkpt, ngkpt, &
       & kptrlatt, kptrlen, nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3,MAX_NSHIFTK)
    real(dp), intent(out) :: kptrlen
    integer, intent(out) :: kptrlatt(3,3)
    integer, intent(out) :: nkpt

    type(symmetry_type), pointer  :: sym
    real(dp), allocatable :: kpt(:,:), wkpt(:)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid1."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! First, compute the number of kpoints
    kptrlatt(:,:) = 0
    kptrlatt(1,1) = ngkpt(1)
    kptrlatt(2,2) = ngkpt(2)
    kptrlatt(3,3) = ngkpt(3)
    kptrlen = 20.

    call getkgrid(0, 0, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB7_MAX_SYMMETRIES, 0, nkpt, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_mp_k_1
!!***


!!****f* m_ab7_kpoints/kpoints_binding_mp_k_2
!! NAME
!!  kpoints_binding_mp_k_2
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_ab7_kpoints
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE
  subroutine kpoints_binding_mp_k_2(symid, nkpt, kpt, wkpt, &
       & kptrlatt, kptrlen, nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3,MAX_NSHIFTK)
    integer, intent(in) :: nkpt
    real(dp), intent(out) :: kpt(3,nkpt), wkpt(nkpt)
    real(dp), intent(inout) :: kptrlen
    integer, intent(inout) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    integer :: nkpt_

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid2."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! Then, we call it again to get the actual values for the k points.
    call getkgrid(0, 0, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB7_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_mp_k_2
!!***


!!****f* m_ab7_kpoints/kpoints_get_mp_k_grid
!! NAME
!! kpoints_get_mp_k_grid
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE
  subroutine kpoints_get_mp_k_grid(symid, nkpt, kpt, wkpt, &
       & ngkpt, nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: ngkpt(3)
    integer, intent(in) :: nshiftk
    real(dp), intent(in) :: shiftk(3, nshiftk)
    integer, intent(out) :: nkpt
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: kptrlen
    integer :: kptrlatt(3,3)
    integer :: nshiftk_
    real(dp) :: shiftk_(3,MAX_NSHIFTK)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get k grid."

    nshiftk_ = nshiftk
    shiftk_(:,1:nshiftk_) = shiftk(:,:)

    call kpoints_binding_mp_k_1(symid, nkpt, ngkpt, kptrlatt, kptrlen, &
         & nshiftk_, shiftk_, errno)
    if (errno /= AB7_NO_ERROR) return
    ABI_ALLOCATE(kpt,(3, nkpt))
    ABI_ALLOCATE(wkpt,(nkpt))
    call kpoints_binding_mp_k_2(symid, nkpt, kpt, wkpt, &
       & kptrlatt, kptrlen, nshiftk_, shiftk_, errno)
  end subroutine kpoints_get_mp_k_grid
!!***


!!****f* m_ab7_kpoints/kpoints_binding_auto_k_1
!! NAME
!!  kpoints_binding_auto_k_1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_ab7_kpoints
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE
  subroutine kpoints_binding_auto_k_1(symid, nkpt, kptrlatt, kptrlen, &
       & nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(out) :: nkpt
    real(dp), intent(inout) :: kptrlen
    integer, intent(out) :: nshiftk
    real(dp), intent(out) :: shiftk(3,MAX_NSHIFTK)
    integer, intent(out) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    real(dp), allocatable :: kpt(:,:), wkpt(:)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid1."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    !  The parameters of the k lattice are not known, compute
    !  kptrlatt, nshiftk, shiftk.
    call testkgrid(sym%bravais,6,kptrlatt,kptrlen,&
         & AB7_MAX_SYMMETRIES,nshiftk,sym%nSym,0,sym%rprimd,&
         & shiftk,sym%symAfm,sym%sym,sym%vacuum)
    if (AB_DBG) write(std_err,*) "AB symmetry: testkgrid -> kptrlatt=", kptrlatt

    nkpt=0
    ABI_ALLOCATE(kpt,(3, nkpt))
    ABI_ALLOCATE(wkpt,(nkpt))

    call getkgrid(0, 0, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB7_MAX_SYMMETRIES, 0, nkpt, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
    if (AB_DBG) write(std_err,*) "AB symmetry: getkgrid -> nkpt=", nkpt

    ABI_DEALLOCATE(kpt)
    ABI_DEALLOCATE(wkpt)

  end subroutine kpoints_binding_auto_k_1
!!***

!!****f* m_ab7_kpoints/kpoints_binding_auto_k_2
!! NAME
!!  kpoints_binding_auto_k_2
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_ab7_kpoints
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE

  subroutine kpoints_binding_auto_k_2(symid, nkpt, kpt, wkpt, kptrlatt, kptrlen, &
       & nshiftk, shiftk, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(in) :: nkpt
    real(dp), intent(out) :: kpt(3,nkpt), wkpt(nkpt)
    real(dp), intent(inout) :: kptrlen
    integer, intent(inout) :: nshiftk
    real(dp), intent(inout) :: shiftk(3,MAX_NSHIFTK)
    integer, intent(inout) :: kptrlatt(3,3)

    type(symmetry_type), pointer  :: sym
    integer :: nkpt_

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid2."

    errno = AB7_NO_ERROR
    call symmetry_get_from_id(sym, symid, errno)
    if (errno /= AB7_NO_ERROR) return

    ! Then, we call it again to get the actual values for the k points.
    call getkgrid(0, 0, 1, kpt, 1, kptrlatt, kptrlen, &
         & AB7_MAX_SYMMETRIES, nkpt, nkpt_, nshiftk, sym%nSym, &
         & sym%rprimd, shiftk, sym%symAfm, sym%sym, &
         & sym%vacuum, wkpt)
  end subroutine kpoints_binding_auto_k_2
!!***


!!****f* m_ab7_kpoints/kpoints_get_auto_k_grid
!! NAME
!!  kpoints_get_auto_k_grid
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      kpoints_binding_auto_k_1,kpoints_binding_auto_k_2
!!
!! SOURCE
  subroutine kpoints_get_auto_k_grid(symid, nkpt, kpt, wkpt, &
       & kptrlen, errno)

    integer, intent(in)  :: symid
    integer, intent(out) :: errno
    integer, intent(out) :: nkpt
    real(dp), intent(in) :: kptrlen
    real(dp), pointer :: kpt(:,:), wkpt(:)

    real(dp) :: kptrlen_
    integer :: kptrlatt(3,3)
    integer :: nshiftk
    real(dp) :: shiftk(3,MAX_NSHIFTK)

    if (AB_DBG) write(std_err,*) "AB symmetry: call get auto k grid."

    kptrlen_ = kptrlen
    call kpoints_binding_auto_k_1(symid, nkpt, kptrlatt, kptrlen_, &
       & nshiftk, shiftk, errno)
    if (errno /= AB7_NO_ERROR) return
    ABI_ALLOCATE(kpt,(3, nkpt))
    ABI_ALLOCATE(wkpt,(nkpt))
    call kpoints_binding_auto_k_2(symid, nkpt, kpt, wkpt, kptrlatt, kptrlen_, &
       & nshiftk, shiftk, errno)
  end subroutine kpoints_get_auto_k_grid
!!***

end module m_ab7_kpoints
!!***
