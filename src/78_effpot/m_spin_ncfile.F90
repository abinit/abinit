!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_ncfile
!! NAME
!! m_spin_ncfile
!!
!! FUNCTION
!! This module contains the wrapper for writting spin hist netcdf file.
!! Unlike the m_spin_terms, inside netcdf, there should be not only the
!! data of magnetic atoms, but also the whole lattice (which do not move).
!!
!! Datatypes:
!!  spin_ncfile_t
!!
!! Subroutines:
!!  * spin_ncfile_t_init
!!  * spin_ncfile_t_write_parameters (write parameters)
!!  * spin_ncfile_t_def_sd (define spin dynamics related dimensions and ids)
!!  * spin_ncfile_t_write_primitive_cell (write primitive cell information)
!!  * spin_ncfile_t_write_supercell (write supercell information)
!!  * spin_ncfile_t_write_one_step (write one step of spin dynamics)
!!  * spin_ncfile_t_close (close and save netcdf file)
!!
!! TODO hexu: should consider carefully what to write.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_spin_ncfile
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_nctk
  use m_spin_hist , only: spin_hist_t
  use m_spin_model_primitive, only: spin_model_primitive_t
  use m_spin_terms , only: spin_terms_t
  use m_multibinit_dataset, only: multibinit_dtset_type
#if defined HAVE_NETCDF
  use netcdf
#endif
  implicit none

!!***

  type spin_ncfile_t
     ! dimensions
     integer :: xyz_dimid, nmatom_dimid, natom_dimid, ntime_dimid, ntypat_dimid
     integer :: itime
     ! file id
     integer :: ncerr,ncid
     ! variable id
     integer :: xred_id,  typat_id, znucl_id,  label_id, spin_index_id
     integer ::  acell_id, rprimd_id
     ! variable ids for spin dynamics
     integer :: entropy_id, etotal_id, S_id, snorm_id, dsdt_id, heff_id, time_id, itime_id
     ! variable ids for spin/lattice coupling
     integer :: ihist_latt_id
  end type spin_ncfile_t

contains

  subroutine spin_ncfile_t_init(self, filename)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_init'
!End of the abilint section

    class(spin_ncfile_t), intent(inout):: self
    character(len=*),intent(in) :: filename
    integer :: ncerr
    self%itime=0

#if defined HAVE_NETCDF
    write(std_out,*) 'Write iteration in spin HIST netCDF file'
    !  Create netCDF file
    ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER,ncid=self%ncid)
    !NCF_CHECK_MSG(ncerr, "create netcdf history file")
#endif
  end subroutine spin_ncfile_t_init

  subroutine spin_ncfile_t_def_sd(self, hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_def_sd'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_hist_t),intent(in) :: hist
    integer :: ncerr
    ! define dimensions
#if defined HAVE_NETCDF
    write(std_out,*) "Defining variables in spin hist netcdf file."
    ncerr=nf90_def_dim(self%ncid, "xyz_dimid", 3, self%xyz_dimid)
    ncerr=nf90_def_dim(self%ncid, "nmatom", hist%nmatom, self%nmatom_dimid )
    ncerr=nf90_def_dim(self%ncid, "ntime_dimid", nf90_unlimited, self%ntime_dimid)

    !ncerr=nf90_def_var(self%ncid, "label", NF90_INT, (/ self%nmatom_dimid/), self%label_id)

    ncerr=nf90_def_var(self%ncid, "S", NF90_DOUBLE, (/ self%xyz_dimid, self%nmatom_dimid, self%ntime_dimid /), self%S_id)
    ncerr=nf90_def_var(self%ncid, "snorm", NF90_DOUBLE, (/ self%nmatom_dimid, self%ntime_dimid /), self%snorm_id)
    ncerr=nf90_def_var(self%ncid, "dsdt", NF90_DOUBLE, (/ self%xyz_dimid, self%nmatom_dimid, self%ntime_dimid /), self%dsdt_id)
    ncerr=nf90_def_var(self%ncid, "Heff", NF90_DOUBLE, (/ self%xyz_dimid, self%nmatom_dimid, self%ntime_dimid /), self%heff_id)
    ncerr=nf90_def_var(self%ncid, "time", NF90_DOUBLE, (/ self%ntime_dimid /), self%time_id)
    ncerr=nf90_def_var(self%ncid, "entropy", NF90_DOUBLE, (/ self%ntime_dimid /), self%entropy_id)
    ncerr=nf90_def_var(self%ncid, "etotal", NF90_DOUBLE, (/ self%ntime_dimid /), self%etotal_id)

    ncerr=nf90_def_var(self%ncid, "itime", NF90_INT, (/ self%ntime_dimid /), self%itime_id)
    ncerr=nf90_def_var(self%ncid, "ihist_latt", NF90_INT, (/ self%ntime_dimid /), self%ihist_latt_id)
    ncerr=nf90_enddef(self%ncid)
#endif
  end subroutine spin_ncfile_t_def_sd

  subroutine spin_ncfile_t_write_one_step(self, hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_write_one_step'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_hist_t), intent(in) :: hist
    integer :: ncerr, itime
    itime=self%itime+1
#if defined HAVE_NETCDF
    !write(std_out, *) "writing spin dynamics step into spin hist netcdf file: itime: ", itime
    ncerr=nf90_put_var(self%ncid, self%S_id, hist%S(:,:,hist%ihist_prev), start=[1, 1, itime], count=[3, hist%nmatom, 1])
    ncerr=nf90_put_var(self%ncid, self%dsdt_id, hist%dsdt(:,:,hist%ihist_prev), start=[1, 1, itime], count=[3, hist%nmatom, 1])
    ncerr=nf90_put_var(self%ncid, self%heff_id, hist%heff(:,:,hist%ihist_prev), start=[1, 1, itime], count=[3, hist%nmatom, 1])
    ncerr=nf90_put_var(self%ncid, self%snorm_id, hist%snorm(:,hist%ihist_prev), start=[1, itime], count=[hist%nmatom, 1])
    ncerr=nf90_put_var(self%ncid, self%ihist_latt_id, [hist%ihist_latt(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%itime_id, [hist%itime(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%time_id, [hist%time(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%entropy_id, [hist%entropy(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%etotal_id, [hist%etot(hist%ihist_prev)], start=[itime], count=[1])
    self%itime=itime
#endif
  end subroutine spin_ncfile_t_write_one_step


  subroutine spin_ncfile_t_write_primitive_cell(self, prim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_write_primitive_cell'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_model_primitive_t) :: prim
    integer :: ncerr

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)
    ncerr=nf90_def_dim(self%ncid, "prim_natom", 1, self%natom_dimid )
    ncerr=nf90_def_dim(self%ncid, "ntime_dimid", nf90_unlimited, self%ntime_dimid)
    !ncerr=nf90_def_dim(self%ncid, "ntypat", 1, self%ntypat_dimid)

    ncerr=nf90_def_var(self%ncid, "prim_acell", NF90_DOUBLE, (/ self%xyz_dimid /), self%acell_id)
    ncerr=nf90_def_var(self%ncid, "prim_rprimd", NF90_DOUBLE, (/self%xyz_dimid, self%xyz_dimid /), self%rprimd_id)
    ncerr=nf90_def_var(self%ncid, "prim_xred", NF90_DOUBLE, (/self%xyz_dimid, self%natom_dimid /), self%xred_id)
    !ncerr=nf90_def_var(self%ncid, "prim_typat", NF90_INT, [self%natom_dimid],  self%typat_id)
    !ncerr=nf90_def_var(self%ncid, "prim_znucl", NF90_DOUBLE, [self%ntypat_dimid],  self%znucl_id)
    ncerr=nf90_def_var(self%ncid, "prim_spin_index", NF90_DOUBLE, [self%natom_dimid], self%spin_index_id)

    ncerr=nf90_enddef(self%ncid)

    !ncerr=nf90_put_var(self%ncid, self%acell_id, prim%unitcell)
    ncerr=nf90_put_var(self%ncid, self%rprimd_id, prim%unitcell)
    ncerr=nf90_put_var(self%ncid, self%xred_id, prim%positions)
    !ncerr=nf90_put_var(self%ncid, typat_id, hist%typat)
    !ncerr=nf90_put_var(self%ncid, znucl_id, hist%znucl)
    ncerr=nf90_put_var(self%ncid, self%spin_index_id, prim%index_spin)
#endif
  end subroutine spin_ncfile_t_write_primitive_cell

  subroutine spin_ncfile_t_write_supercell(self, scell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_write_supercell'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_terms_t), intent(in) :: scell
    integer :: unitcell_id, pos_id, ispin_prim_id, rvec_id, ncerr
    ! sc_matric

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)
    ncerr=nf90_def_var(self%ncid, "unitcell", NF90_DOUBLE, [self%xyz_dimid, self%xyz_dimid], unitcell_id)
    ncerr=nf90_def_var(self%ncid, "pos", NF90_DOUBLE, [self%xyz_dimid, self%nmatom_dimid], pos_id)
    ncerr=nf90_def_var(self%ncid, "ispin_prim", NF90_INT, [self%nmatom_dimid], ispin_prim_id)
    ncerr=nf90_def_var(self%ncid, "rvec", NF90_INT, [self%xyz_dimid,self%nmatom_dimid], rvec_id)
    ncerr=nf90_enddef(self%ncid)

    ncerr=nf90_put_var(self%ncid, unitcell_id, scell%cell)
    ncerr=nf90_put_var(self%ncid, pos_id, scell%pos)
    ncerr=nf90_put_var(self%ncid, ispin_prim_id, scell%ispin_prim)
    ncerr=nf90_put_var(self%ncid, rvec_id, scell%rvec)

#endif
  end subroutine spin_ncfile_t_write_supercell

  subroutine spin_ncfile_t_write_parameters(self, params)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_write_parameters'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    integer :: qpoint_id, temperature_id, dt_id, mfield_id, ncell_id
    integer :: ncerr

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)

    ! dims 
    ! vars
    ncerr=nf90_def_var(self%ncid, "spin_qpoint", NF90_DOUBLE, [self%xyz_dimid], qpoint_id)
    ncerr=nf90_def_var(self%ncid, "ncell", NF90_INT, [self%xyz_dimid], ncell_id)
    ncerr=nf90_def_var(self%ncid, "spin_temperature", NF90_DOUBLE,  temperature_id)
    ncerr=nf90_def_var(self%ncid, "spin_dt", NF90_DOUBLE,  dt_id)
    ncerr=nf90_def_var(self%ncid, "spin_mag_field", NF90_DOUBLE, [self%xyz_dimid],  mfield_id)

    ncerr=nf90_enddef(self%ncid)
    ! put vars
    ncerr=nf90_put_var(self%ncid, qpoint_id, params%spin_qpoint)
    ncerr=nf90_put_var(self%ncid, ncell_id, params%ncell)
    ncerr=nf90_put_var(self%ncid, temperature_id, params%spin_temperature)
    ncerr=nf90_put_var(self%ncid, dt_id, params%spin_dt)
    ncerr=nf90_put_var(self%ncid, mfield_id, params%spin_mag_field)
#endif
  end subroutine spin_ncfile_t_write_parameters

  subroutine spin_ncfile_t_close(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ncfile_t_close'
!End of the abilint section

    class(spin_ncfile_t), intent(inout) :: self
    integer :: ncerr
#if defined HAVE_NETCDF
    write(std_out, *) "Closing spin hist file"
    ncerr=nf90_close(self%ncid)
    !NCF_CHECK_MSG(ncerr, "close netcdf spin history file")
#endif
  end subroutine spin_ncfile_t_close

end module m_spin_ncfile
