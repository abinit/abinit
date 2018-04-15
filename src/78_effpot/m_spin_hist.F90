
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_hist
!! NAME
!! m_spin_hist
!!
!! FUNCTION
!! This module contains definition the type spin_hist
!! and its related routines
!!
!! Datatypes:
!!
!! * spin_hist: Historical record of spin orientations and amplitudes
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

module m_spin_hist

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

  !!****t* m_spin_hist/spin_hist_t
  !! NAME
  !! spin_hist_t
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

  type, public :: spin_hist_t
     ! scalars
     ! Index of the last element on all records
     integer :: ihist = 0
     integer :: ihist_prev = -1
     ! Maximun size of the historical records
     integer :: mxhist = 0

     integer :: natom
     integer :: nmatom
     ! whether lattice dynamics is also present
     interger, allocatable :: ihist_latt(:)
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
     ! TODO hexu: is it useful?
     real(dp), allocatable :: dSdt(3, nmatom, mxhist)

     ! etot(mxhist)
     real(dp), allocatable :: etot(:)
     real(dp) :: entropy(:)
     real(dp) :: time(:)

  end type spin_hist_t

  public :: spin_hist_t_init
  public :: spin_hist_t_free
  public :: spinhist2var
  public :: var2spinhist
  public :: spin_hist_t_findIndex
  public :: write_sd_hist
  public :: read_md_hist
  public :: get_dims_spinhist

contains
  subroutine spin_hist_t_init(hist, natom, nmatom, mxhist, has_latt)
    implicit none
    class(spin_hist_t), intent(inout) :: self
    integer, intent(in) :: natom, nmatom, mxhist
    logical, intent(in) :: has_latt

    hist%ihist=1
    hist%ihist_prev=0
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
  end subroutine spin_hist_t_init

  subroutine spin_hist_t_free(hist)
    implicit none
    class(spin_hist_t) , intent(inout) :: hist

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
  end subroutine spin_hist_t_free


  subroutine spin_hist_t_get_S(hist, S, step)
    class(spin_hist_t), intent(in) :: hist
    real(dp), intent(out) :: S(3, hist%nmatoms)
    integer, intent(in), optional:: step
    integer :: i
    if (.not. present(step)) then
       step=0
    end if
    i=spin_hist_t_findIndex(step=step)
    S(:,:)=hist%S(:,:,i)
  end subroutine spin_hist_t_get_S

  subroutine spin_hist_t_inc(hist)
    class(spin_hist_t), intent(inout) :: hist
    hist%ihist_prev=hist%ihist
    hist%ihist=spin_hist_t_findIndex(1)
  end subroutine spin_hist_t_inc

  function spin_hist_t_findIndex(hist, step) result(index)
    type(spin_hist_t), intent(inout) :: hist
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
  end function spin_hist_t_findIndex

  subroutine spin_hist_t_set_vars(hist, S, Snorm, dSdt, Heff, etot, entropy, time, ihist_latt, inc)
    class(spin_hist_t), intent(inout) :: hist
    real(dp), intent(in), optional :: S(3, hist%nmatoms), Snorm(hist%nmatom), &
         & dSdt(3, hist%nmatoms), Heff(3, hist%nmatoms), etot, entropy, time
    integer, optional :: ihist_latt
    logical, intent(in), optional :: inc
    if(present(inc) .and. inc) then
       call spin_hist_t_inc(hist)
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
    if(present(etot)) then
       hist%etot( ihist)=etot
    end if
    if(present(entropy)) then
       hist%entropy(ihist)=entropy
    end if
    if(present(time)) then
       hist%time( ihist)=time
    end if
    if(present(ihist_latt)) then
       hist%ihist_latt(ihist)=ihist_latt
    endif
  end subroutine spin_hist_t_set_vars

  subroutine write_spin_hist(hist, filename, ifirst, itime, natom, nmatom, ntypat, & typat, amu, znucl, dtspin, sdtemp, index_spin)
    integer,intent(in) :: ifirst,itime,natom,ntypat
    real(dp),intent(in) :: dtion
    character(len=*),intent(in) :: filename
    !arrays
    integer,intent(in) :: typat(natom), index_spin(natom)
    real(dp),intent(in) :: amu(ntypat),znucl(:),mdtemp(2)
    type(spin_hist_t),intent(inout),target :: hist
  end subroutine write_spin_hist

end module m_spin_hist
