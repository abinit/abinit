!!****m* ABINIT/m_lwf_ncfile
!! NAME
!! m_lwf_ncfile
!!
!! FUNCTION
!! This module contains the subroutines for output netcdf file
!!
!!
!! Datatypes:
!! lwf_ncfile_t: store data to calculate lwf_ncfile
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

module m_lwf_ncfile

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_nctk
  use m_lwf_hist , only: lwf_hist_t
  use m_lwf_primitive_potential, only: lwf_primitive_potential_t
  use m_lwf_potential , only: lwf_potential_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_lwf_observables, only : lwf_observables_t
!#if defined HAVE_NETCDF
  use netcdf
!#endif
  implicit none

  private
  !!***
  type, public :: lwf_ncfile_t
     logical :: isopen=.False.  ! if the file is open
     ! dimensions
     integer :: three, nlwf
     ! three: 3

     ! file id
     integer :: ncerr, ncid
     ! variable id
     integer :: heff_id, time_id, itime_id
     integer :: itime
     ! itime: time index
     integer :: write_traj=0
     !whether to write the trajectory 
     character(len=fnlen) :: filename
     ! netcdf filename
   contains
     ! initialize
     procedure :: initialize    
     procedure :: finalize
  end type lwf_ncfile_t

contains
  subroutine initialize(self, filename, write_traj)
    class(lwf_ncfile_t) :: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: write_traj
    integer :: ncerr
    self%itime=0
    self%write_traj=write_traj
    self%filename=trim(filename)
    self%isopen=.False.
!#if defined HAVE_NETCDF
    write(std_out,*) "Write iteration in lwf history file "//trim(self%filename)//"."
    !  Create netCDF file
    ncerr = nf90_create(path=trim(filename), cmode=NF90_CLOBBER, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating netcdf history file")
    self%isopen=.True.
    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in spin netcdf history file")
!#endif
  end subroutine initialize

  subroutine finalize(self)
    class(lwf_ncfile_t), intent(inout) :: self
!#if defined HAVE_NETCDF
    integer :: ncerr
    if (self%isopen) then
       write(std_out, *) "Closing lwf history file "//trim(self%filename)//"."
       ncerr=nf90_close(self%ncid)
       NCF_CHECK_MSG(ncerr, "close netcdf lwf history file"//trim(self%filename)//".")
    end if
!#endif
  end subroutine finalize

end module m_lwf_ncfile
