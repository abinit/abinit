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
  use m_spin_observables, only : spin_observable_t
#if defined HAVE_NETCDF
  use netcdf
#endif
  implicit none

  !!***

  type spin_ncfile_t
     ! dimensions
     integer :: three, nspins, natoms, ntime, ntypat, nsublatt
     ! file id
     integer :: ncerr, ncid
     ! variable id
     integer :: xred_id,  typat_id, znucl_id,  label_id, spin_index_id
     integer ::  acell_id, rprimd_id
     ! variable ids for spin dynamics
     integer :: entropy_id, etotal_id, S_id, snorm_id, dsdt_id, heff_id, time_id, itime_id
     integer :: Mst_sub_id, Mst_sub_norm_id, Mst_norm_total_id, Snorm_total_id

     ! thermo obs
     integer :: chi_id, binderU4_id, Cv_id
     ! variable ids for spin/lattice coupling
     integer :: ihist_g_id

     integer :: itime
     integer :: write_traj
     character(len=fnlen) :: filename
  end type spin_ncfile_t

contains

  subroutine spin_ncfile_t_init(self, filename, write_traj)

    class(spin_ncfile_t), intent(inout):: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: write_traj
    integer :: ncerr

    self%itime=0
    self%write_traj=write_traj
    self%filename=trim(filename)

#if defined HAVE_NETCDF
    write(std_out,*) 'Write iteration in spin HIST netCDF file'
    !  Create netCDF file
    ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER,ncid=self%ncid)
    !NCF_CHECK_MSG(ncerr, "create netcdf history file")
#endif
  end subroutine spin_ncfile_t_init

  subroutine spin_ncfile_t_def_sd(self, hist)

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_hist_t),intent(in) :: hist
    integer :: ncerr
    ! define dimensions
#if defined HAVE_NETCDF
    write(std_out,*) "Defining variables in spin hist netcdf file."
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%three)
    ncerr=nf90_def_dim(self%ncid, "nspins", hist%nspins, self%nspins )
    ncerr=nf90_def_dim(self%ncid, "ntime", nf90_unlimited, self%ntime)
    !call ab_define_var(ncid,dim1,typat_id,NF90_DOUBLE,&
    !     &  "typat","types of atoms","dimensionless" )

    if(self%write_traj==1) then
       call ab_define_var(self%ncid, (/ self%three, self%nspins, self%ntime /), &
            &         self%S_id, NF90_DOUBLE, "S", "Spin orientations", "dimensionless")
       call ab_define_var(self%ncid, (/ self%nspins, self%ntime /), &
            &         self%snorm_id, NF90_DOUBLE, "snorm", "Spin norm2", "Mu_B")
       call ab_define_var(self%ncid, (/ self%three, self%nspins, self%ntime /), &
            &         self%dsdt_id, NF90_DOUBLE, "dsdt", "Spin orientations derivative to time", "1/s")
       call ab_define_var(self%ncid, (/ self%three, self%nspins, self%ntime /), &
            &         self%Heff_id, NF90_DOUBLE, "Heff", "Effective spin torque", "Tesla")
    endif

    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%time_id, NF90_DOUBLE, "time", "time", "s")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%etotal_id, NF90_DOUBLE, "etotal", "TOTAL energy", "Joule")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%entropy_id, NF90_DOUBLE, "entropy", "Entropy", "Joule/K")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%itime_id, NF90_INT, "itime", "index of time in spin timeline", "1")
    ncerr=nf90_enddef(self%ncid)
#endif
  end subroutine spin_ncfile_t_def_sd

  subroutine spin_ncfile_t_def_ob(self, ob)
    class(spin_ncfile_t), intent(inout) :: self
    type(spin_observable_t), intent(in) :: ob
    integer ncerr

#if defined HAVE_NETCDF
    ncerr = nf90_redef(self%ncid)
    ncerr = nf90_def_dim(self%ncid, "nsublatt", ob%nsublatt, self%nsublatt)

    call ab_define_var(self%ncid, (/self%three, self%nsublatt, self%ntime/),& 
           & self%Mst_sub_id, NF90_DOUBLE, "Mst_sub", "Sublattice staggered M", "Bohr magneton")
    call ab_define_var(self%ncid, (/ self%nsublatt, self%ntime/), & 
          &  self%Mst_sub_norm_id, NF90_DOUBLE, "Mst_sub_norm", &
          &  "Norm of sublattice staggered M", "Bohr magneton")
    call ab_define_var(self%ncid, (/self%ntime/), self%Mst_norm_total_id, &
           & NF90_DOUBLE, "Mst_norm_total", "total Norm of sublattice M", "Bohr magneton")
    call ab_define_var(self%ncid, (/self%ntime/), self%Snorm_total_id, &
           & NF90_DOUBLE, "Snorm_sub", "Snorm of sublattice", "Bohr magneton")

    if(ob%calc_thermo_obs)then
       call ab_define_var(self%ncid, (/self%ntime/), self%binderU4_id, & 
               & NF90_DOUBLE, "BinderU4", "Binder U4", "1")
       call ab_define_var(self%ncid, (/self%ntime/), self%Cv_id, &
               & NF90_DOUBLE, "Cv", "Specific heat", "Joule/K")
       call ab_define_var(self%ncid, (/self%ntime/), self%chi_id, &
               &NF90_DOUBLE, "chi", "magnetic susceptibility", "1")
    endif

    if(ob%calc_traj_obs)then
     !TODO define traj obs here
  endif
  if(ob%calc_correlation_obs)then
     !TODO define correlation obs here
  endif

  ncerr=nf90_enddef(self%ncid)

#endif
end subroutine spin_ncfile_t_def_ob

  subroutine spin_ncfile_t_write_one_step(self, hist, ob)

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_hist_t), intent(in) :: hist
    type(spin_observable_t), optional, intent(in) :: ob
    integer :: ncerr, itime
    itime=self%itime+1
#if defined HAVE_NETCDF
    !write(std_out, *) "writing spin dynamics step into spin hist netcdf file: itime: ", itime
    if(self%write_traj ==1) then
       ncerr=nf90_put_var(self%ncid, self%S_id, hist%S(:,:,hist%ihist_prev), &
            &      start=[1, 1, itime], count=[3, hist%nspins, 1])

       ncerr=nf90_put_var(self%ncid, self%dsdt_id, &
            &      hist%dsdt(:,:,hist%ihist_prev), start=[1, 1, itime], &
            &      count=[3, hist%nspins, 1])

       ncerr=nf90_put_var(self%ncid, self%heff_id, &
            &      hist%heff(:,:,hist%ihist_prev), start=[1, 1, itime], count=[3, hist%nspins, 1])

       ncerr=nf90_put_var(self%ncid, self%snorm_id, &
            &      hist%snorm(:,hist%ihist_prev), start=[1, itime], count=[hist%nspins, 1])
    end if
    !ncerr=nf90_put_var(self%ncid, self%ihist_g_id, [hist%ihist_latt(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%itime_id, &
         &       [hist%itime(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%time_id, &
         & [hist%time(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%entropy_id, &
         & [hist%entropy(hist%ihist_prev)], start=[itime], count=[1])
    ncerr=nf90_put_var(self%ncid, self%etotal_id,  &
         & [hist%etot(hist%ihist_prev)], start=[itime], count=[1])
    self%itime=itime

    if(present(ob)) then
       ncerr=nf90_put_var(self%ncid, self%Mst_sub_id, ob%Mst_sub, &
            &      start=[1, 1, itime], count=[3, ob%nsublatt, 1])
       ncerr=nf90_put_var(self%ncid, self%Mst_sub_norm_id, ob%Mst_sub_norm, &
            &      start=[ 1, itime], count=[ob%nsublatt, 1])
       ncerr=nf90_put_var(self%ncid, self%Mst_norm_total_id, [ob%Mst_norm_total], &
            &      start=[itime], count=[1])

       ncerr=nf90_put_var(self%ncid, self%Snorm_total_id, [ob%Snorm_total], &
            &      start=[itime], count=[1])

       if(ob%calc_traj_obs)then
       endif
       if(ob%calc_thermo_obs)then
       endif
       if(ob%calc_correlation_obs)then
       endif
    end if

#endif
  end subroutine spin_ncfile_t_write_one_step

  subroutine spin_ncfile_t_write_primitive_cell(self, prim)

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_model_primitive_t) :: prim
    integer :: ncerr

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)
    ncerr=nf90_def_dim(self%ncid, "prim_natoms", 1, self%natoms )
    ncerr=nf90_def_dim(self%ncid, "ntime", nf90_unlimited, self%ntime)
    !ncerr=nf90_def_dim(self%ncid, "ntypat", 1, self%ntypat)

    ncerr=nf90_def_var(self%ncid, "prim_acell", NF90_DOUBLE, &
         & (/ self%three /), self%acell_id)
    ncerr=nf90_def_var(self%ncid, "prim_rprimd", NF90_DOUBLE, &
         & (/self%three, self%three /), self%rprimd_id)
    ncerr=nf90_def_var(self%ncid, "prim_xred", NF90_DOUBLE, &
         & (/self%three, self%natoms /), self%xred_id)
    !ncerr=nf90_def_var(self%ncid, "prim_typat", NF90_INT, [self%natoms],  self%typat_id)
    !ncerr=nf90_def_var(self%ncid, "prim_znucl", NF90_DOUBLE, [self%ntypat],  self%znucl_id)
    ncerr=nf90_def_var(self%ncid, "prim_spin_index", NF90_DOUBLE, &
         & [self%natoms], self%spin_index_id)

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

    class(spin_ncfile_t), intent(inout) :: self
    type(spin_terms_t), intent(in) :: scell
    integer :: rprimd_id, pos_id, ispin_prim_id, rvec_id, iatoms_id, ncerr
    ! sc_matric

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)

    call ab_define_var(self%ncid, (/self%three, self%three /), rprimd_id, &
         & NF90_DOUBLE, "rprimd", "primitive cell vectors in real space with units", "bohr")
    call ab_define_var(self%ncid, (/self%three, self%nspins/), pos_id, &
         & NF90_DOUBLE, "xcart_spin","position of spin in cartesian coordinates", "bohr")
    call ab_define_var(self%ncid, (/ self%nspins/), ispin_prim_id, &
         & NF90_INT, "ispin_prim", "index of spin in primitive cell", "dimensionless")
    call ab_define_var(self%ncid, (/self%three, self%nspins/), rvec_id, &
         & NF90_INT, "Rvec", "R vector for spin in supercell", "dimensionless")
    call ab_define_var(self%ncid, (/ self%nspins/), iatoms_id, &
         & NF90_INT, "iatoms", "indices of atoms with spin", "dimensionless")

    !ncerr=nf90_def_var(self%ncid, "unitcell", NF90_DOUBLE, [self%three, self%three], unitcell_id)
    !ncerr=nf90_def_var(self%ncid, "xred", NF90_DOUBLE, [self%three, self%nspins], pos_id)
    !ncerr=nf90_def_var(self%ncid, "ispin_prim", NF90_INT, [self%nspins], ispin_prim_id)
    !ncerr=nf90_def_var(self%ncid, "rvec", NF90_INT, [self%three,self%nspins], rvec_id)
    ncerr=nf90_enddef(self%ncid)

    ncerr=nf90_put_var(self%ncid, rprimd_id, scell%cell)
    ncerr=nf90_put_var(self%ncid, pos_id, scell%pos)
    ncerr=nf90_put_var(self%ncid, ispin_prim_id, scell%ispin_prim)
    ncerr=nf90_put_var(self%ncid, rvec_id, scell%rvec)
    ncerr=nf90_put_var(self%ncid, iatoms_id, scell%iatoms)
#endif
  end subroutine spin_ncfile_t_write_supercell

  subroutine spin_ncfile_t_write_parameters(self, params)

    class(spin_ncfile_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    integer :: qpoint_id, temperature_id, dt_id, mfield_id, ncell_id
    integer :: dim0(0)
    integer :: ncerr

#if defined HAVE_NETCDF
    ncerr=nf90_redef(self%ncid)

    ! dims 
    ! vars
    call ab_define_var(self%ncid, (/self%three/), qpoint_id, &
         & NF90_DOUBLE, "spin_qpoint", "spin QPOINT", "dimensionless")
    ! TODO should change ncell to 3*3 matrix
    call ab_define_var(self%ncid, (/self%three/), ncell_id, &
         & NF90_INT, "ncell", "supercell matrix (only diagonal)", "dimensionless")
    call ab_define_var(self%ncid, dim0, temperature_id, &
         & NF90_DOUBLE, "spin_temperature", "Spin temperature", "Kelvin")
    call ab_define_var(self%ncid, dim0, dt_id, &
         & NF90_DOUBLE, "spin_dt", "Spin time step", "second")
    call ab_define_var(self%ncid, (/self%three/), mfield_id, &
         & NF90_DOUBLE, "spin_mag_field", "magnetic field for spin dynamics", "Tesla")
    !ncerr=nf90_def_var(self%ncid, "spin_qpoint", NF90_DOUBLE, [self%three], qpoint_id)
    !ncerr=nf90_def_var(self%ncid, "ncell", NF90_INT, [self%three], ncell_id)
    !ncerr=nf90_def_var(self%ncid, "spin_temperature", NF90_DOUBLE,  temperature_id)
    !ncerr=nf90_def_var(self%ncid, "spin_dt", NF90_DOUBLE,  dt_id)
    !ncerr=nf90_def_var(self%ncid, "spin_mag_field", NF90_DOUBLE, [self%three],  mfield_id)

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

    class(spin_ncfile_t), intent(inout) :: self
    integer :: ncerr
#if defined HAVE_NETCDF
    write(std_out, *) "Closing spin hist file"//trim(self%filename)//"."
    ncerr=nf90_close(self%ncid)
    !NCF_CHECK_MSG(ncerr, "close netcdf spin history file"//trim(self%filename)//".")
#endif
  end subroutine spin_ncfile_t_close


end module m_spin_ncfile
