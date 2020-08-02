!!****m* ABINIT/m_lwf_hist
!! NAME
!! m_lwf_hist
!!
!! FUNCTION
!! This module contains definition the type lwf_hist_t
!! and its related routines
!! The observables are also calculated. 
!!
!! Datatypes:
!!
!! * lwf_hist_t: history record of lwf orientations and amplitudes
!!
!! Subroutines:
!!
!! * lwf_hist_t
!! * set_params
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
module m_lwf_hist
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_lwf_observables
  !use m_lwf_ncfile
  implicit none

  private
!!***

  type, public :: lwf_hist_t
     integer :: mxhist = 0
     integer :: nlwf = 0
     integer :: ihist=0
     real(dp), allocatable :: hist(:,:)
     real(dp), allocatable :: energy(:)
     real(dp), pointer :: current_lwf(:)
     real(dp), pointer :: current_energy
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_hist
     procedure :: get_hist
  end type lwf_hist_t

contains

  subroutine initialize(self, nlwf, mxhist)
    class(lwf_hist_t), intent(inout) :: self
    integer, intent(in) :: nlwf, mxhist
    self%nlwf=nlwf
    self%mxhist=mxhist
    ABI_MALLOC(self%hist, (nlwf, mxhist))
    ABI_MALLOC(self%energy, (mxhist))
    self%hist(:, :)=0.0_dp
  end subroutine initialize

  subroutine finalize(self)
    class(lwf_hist_t), intent(inout) :: self
    self%mxhist=0
    self%nlwf=0
    self%ihist=0
    nullify(self%current_lwf)
    nullify(self%current_energy)
    ABI_SFREE(self%hist)
    ABI_SFREE(self%energy)
  end subroutine finalize

  subroutine set_hist(self, lwf, energy)
    class(lwf_hist_t), target, intent(inout) :: self
    real(dp), intent(in) :: lwf(:), energy
    self%ihist=modulo(self%ihist+1, self%mxhist)+1
    self%hist(:,self%ihist)=lwf(:)
    self%energy(self%ihist)=energy
    self%current_lwf => self%hist(:,self%ihist)
    self%current_energy => self%energy(self%ihist)
  end subroutine set_hist


  function get_hist(self, rel_ihist) result(lwf_out)
    class(lwf_hist_t), target, intent(inout) :: self
    integer, optional, intent(in) :: rel_ihist
    real(dp), pointer :: lwf_out(:)
    integer ::i
    if (present(rel_ihist)) then
       if (rel_ihist>0 .or. abs(rel_ihist)>self%mxhist) then
          MSG_BUG("Asking for lwf hist which is beyond mxhist.")
       end if
       i=modulo(self%ihist+rel_ihist, self%mxhist)+1
    else
       i=self%ihist
    end if
    lwf_out => self%hist(:, i)
  end function get_hist

end module m_lwf_hist
