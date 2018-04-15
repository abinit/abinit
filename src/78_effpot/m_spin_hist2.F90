
!{\src2tex{textfont=tt}}
  !!****m* ABINIT/m_abihist
  !! NAME
  !! m_abihist
  !!
  !! FUNCTION
  !! This module contains definition the type abihist_spin
  !! and its related routines
  !!
  !! Datatypes:
  !!
  !! * abihist: Historical record of spin orientations and amplitudes
  !!
  !! Subroutines:
  !!
  !! * abispinhist_init
  !! * abispinhist_free
  !! * abispinhist_bcast
  !! * abispinhist_compare
  !! * spinhist2var
  !! * var2spinhist
  !! * dSdt2spinhist
  !!
  !! COPYRIGHT
  !! Copyright (C) 2001-2017 ABINIT group (XG, SE)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !!
  !! SOURCE

! TODO hexu:
! sync ihist_latt when lattice dynamics
! add average , variance, etc (should they be here?)

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abihist

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_nctk
#if defined HAVE_NETCDF
 use netcdf
#endif

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_abihist/abihist
!! NAME
!! abihist
!!
!! FUNCTION
!! This type has several vectors, and index scalars to store
!! a proper history of previous evaluations of forces and
!! stresses,velocities,positions and energies
!!
!! It contains:
!! * mxhist                  : Maximum size of history
!! * ihist                   : index of history

 !! natom : number of atoms
 !! nmatoms: number of magnetic atoms
!! * acell(3)         : Acell (acell , rprimd, xred: only initial value kept if there is!!  no lattice dynamics. Other wise for each step, the corresponding lattice step number is kept)
!! * rprimd(3,3)      : Rprimd
!! * xred(3,natom)    : Xred
!! * index_spin     : the index of atom in spin model, -1 if it is not in the spin model 
!! * heff(3,nmatom,mxhist)   : effective magnetic field (cartesian)
!! * snorm(nmatom, mxhist) : magnetitude of spin.
!! * S(3,nmatom,mxhist)   : spin orientation of atoms (cartesian)
!! * dSdt(3, nmatom, mxhist) : dS/dt (cartesian)
!! * etot(mxhist)            : Electronic total Energy
!! * entropy(mxhist)         : Entropy
!! * time(mxhist)            : Time (or iteration number for GO)
!!
!! * has_latt (whether lattice dynamics is also present)
!! * ihist_latt(mxhist): the corresponding lattice step. 0 if none.
!! SOURCE

 type, public :: abispinhist
    ! scalars
    ! Index of the last element on all records
    integer :: ihist = 0
    ! Maximun size of the historical records
    integer :: mxhist = 0

    integer :: natom
    integer :: nmatom
    ! whether lattice dynamics is also present
    logical :: has_latt

    ! arrays
    ! structure
    real(dp) :: acell(3)
    real(dp) :: rprimd(3,3)
    ! xred(3, natom)
    real(dp), allocatable :: xred(3, :)

    ! spin
    !heff(3, nmatom, mxhist)
    real(dp), allocatable :: heff(3, :, :)
    !snorm(nmatom, mxhist)
    real(dp), allocatable :: snorm(:, :)

    !S(3, nmatom, mxhist)
    real(dp), allocatable :: S(3, nmatom, mxhist)
    !dSdt(3, nmatom, mxhist)
    real(dp), allocatable :: dSdt(3, nmatom, mxhist)

    ! etot(mxhist)
    real(dp), allocatable :: etot(:)
    real(dp) :: entropy(:)
    real(dp) :: time(:)

 end type abispinhist

 public :: abispinhist_init
 public :: abispinhist_free
 public :: spinhist2var
 public :: var2spinhist
 public :: abispinhist_findIndex
 public :: write_sd_hist
 public :: read_md_hist
 public :: get_dims_spinhist

contains
  subroutine abispinhist_init(hist, natom, nmatom, mxhist, has_latt)
    implicit none
    class(abispinhist), intent(inout) :: hist
    integer, intent(in) :: natom, nmatom, mxhist
    logical, intent(in) :: has_latt

    hist%ihist=1
    hist%mxhist=mxhist

    hist%has_latt=has_latt

    ABI_ALLOCATE(hist%xred, (3, natoms))
    ABI_ALLOCATE(hist%heff, (3, nmatoms, mxhist))
    ABI_ALLOCATE(hist%snorm, (nmatoms, mxhist))
    ABI_ALLOCATE(hist%S, (3, nmatoms, mxhist))
    ABI_ALLOCATE(hist%dSdt, (3, nmatoms, mxhist))

    ABI_ALLOCATE(hist%etot, (mxhist))
    ABI_ALLOCATE(hist%entropy, (mxhist))
    ABI_ALLOCATE(hist%time, (mxhist))
    hist%etot(1) =zero
    hist%entropy(1) =zero
    hist%time(1) =zero

    hist%acell(:)=zero
    hist%rprimd(:, :)=zero
    hist%xred(:,:) =zero
    hist%heff(:,:,1)=zero
    hist%S(:,:,1)=zero
    hist%dSdt(:,:,1)=zero
    hist%snorm(:,:,1)=zero
  end subroutine abispinhist_init

  subroutine abispinhist_free(hist)
    implicit none
    class(abispinhist) , intent(inout) :: hist

    if (allocated(hist%xred)) then
       ABI_DEALLOCATE(hist%xred)
    end if
    if (allocated(hist%heff)) then
       ABI_DEALLOCATE(hist%heff)
    end if
    if (allocated(hist%snorm)) then
       ABI_DEALLOCATE(hist%snorm)
    end if
    if (allocated(hist%S)) then
       ABI_DEALLOCATE(hist%S)
    end if
    if (allocated(hist%dSdt)) then
       ABI_DEALLOCATE(hist%dSdt)
    end if
    if (allocated(hist%etot)) then
       ABI_DEALLOCATE(hist%etot)
    end if
    if (allocated(hist%entropy)) then
       ABI_DEALLOCATE(hist%entropy)
    end if
    if (allocated(hist%time)) then
       ABI_DEALLOCATE(hist%time)
    end if
  end subroutine abispinhist_free


  subroutine abispinhist_get_S(hist, S, ihist)
    class (abispinhist), intent(in) :: hist
    real(dp), intent(out) :: S(3, hist%nmatoms)
    integer, intent(in), optional:: ihist
    if (.not. present(ihist)) then
       ihist=hist%ihist
    end if
    S(:,:)=hist%S(:,:,ihist)
  end subroutine abispinhist_get_S

  subroutine abispinhist_inc(hist)
    class (abispinhist), intent(inout) :: hist
    hist%ihist=abispinhist_findIndex(1)
  end subroutine abispinhist_inc

  function abispinhist_findIndex(hist, step) result(index)
    type(abispinhist), intent(inout) :: hist
    integer , intent(in) :: step
    integer :: index
    !Local variables-------------------------------
    !scalars
    integer :: ii,mxhist
    !arrays
    character(len=500) :: msg
    ! *************************************************************

    mxhist = hist%mxhist
    if ((mxhist ==1.and.step/=+1).or.&
         &    (mxhist /=1.and.abs(step) >=mxhist)) then
       write(msg,'(a,I0,2a)')' The requested step must be lass than ',mxhist,ch10,&
            &                     'Action: increase the number of history store in the hist'
      MSG_BUG(msg)
    end if
    index= mod(hist%ihist+step, hist%mxhist)
  end function abispinhist_findIndex

  subroutine abispinhist_set_vars(hist, S, Snorm, dSdt, Heff, inc)
    class (abispinhist), intent(inout) :: hist
    real(dp), intent(in), optional :: S(3, hist%nmatoms), Snorm(hist%nmatom), dSdt(3, hist%nmatoms), Heff(3, hist%nmatoms)
    logical, intent(in), optional :: inc
    if(present(inc) .and. inc) then
       call abispinhist_inc(hist)
    end if
    if(present(S)) then
       hist%S(:, :, ihist)=S(:,:)
    end if
    if(present(Snorm)) then
       hist%Snorm(:,  ihist)=Snorm(:)
    endif
    if(present(dSdt)) then
       hist%dSdt(:, :, ihist)=dSdt(:,:)
    end if
    if(present(Heff)) then
       hist%Heff(:, :, ihist)=Heff(:,:)
    end if
  end subroutine abispinhist_set_vars

  subroutine write_sd_hist(hist, filename, ifirst, itime, natom, nmatom, ntypat, & typat, amu, znucl, dtspin, sdtemp, index_spin)
    integer,intent(in) :: ifirst,itime,natom,ntypat
    real(dp),intent(in) :: dtion
    character(len=*),intent(in) :: filename
    !arrays
    integer,intent(in) :: typat(natom), index_spin(natom)
    real(dp),intent(in) :: amu(ntypat),znucl(:),mdtemp(2)
    type(abispinhist),intent(inout),target :: hist
  end subroutine write_sd_hist


