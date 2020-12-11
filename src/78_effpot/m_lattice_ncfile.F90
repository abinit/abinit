!!****m* ABINIT/m_lattice_ncfile
!! NAME
!! m_lattice_ncfile
!!
!! FUNCTION
!! This module contains the subroutines for output netcdf file
!!
!!
!! Datatypes:
!! lattice_ncfile_t: store data to calculate lattice_ncfile
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

#define HAVE_NETCDF 1

module m_lattice_ncfile

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_nctk
  use m_lattice_harmonic_primitive_potential, only: lattice_harmonic_primitive_potential_t
  use m_lattice_harmonic_potential , only: lattice_harmonic_potential_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
#if defined HAVE_NETCDF
  use netcdf
#endif
  implicit none

  private
  !!***
  type, public :: lattice_ncfile_t
     logical :: isopen=.False.  ! if the file is open
     ! dimensions
     integer :: three, natom
     ! three: 3

     ! file id
     integer :: ncerr, ncid
     ! variable id
     integer :: vcart_id, xcart_id, cell_id, time_id, itime_id, etotal_id, ekin_id, entropy_id
     integer :: ilatt_prim_id, rvec_id, zion_id, masses_id, ref_xcart_id, ref_cell_id
     integer :: itime, ntime,ncell
     ! itime: time index
     integer :: write_traj=1
     !whether to write the trajectory
     character(len=fnlen) :: filename
     ! netcdf filename
   contains
     ! initialize
     procedure :: initialize
     procedure :: finalize
     procedure :: def_lattice_var
     procedure :: write_cell
     procedure :: write_one_step
  end type lattice_ncfile_t

contains
  subroutine initialize(self, filename, write_traj)
    class(lattice_ncfile_t) :: self
    character(len=*),intent(in) :: filename
    integer, intent(in) :: write_traj
    integer :: ncerr
    self%itime=0
    self%write_traj=write_traj
    self%filename=trim(filename)
#if defined HAVE_NETCDF
    write(std_out,*) "Write iteration in lattice history file "//trim(self%filename)//"."
    !  Create netCDF file
    ncerr = nf90_create(path=trim(filename), cmode=NF90_CLOBBER, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating netcdf history file")
    self%isopen=.True.
    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in lattice netcdf history file")
#endif
  end subroutine initialize

  subroutine finalize(self)
    class(lattice_ncfile_t), intent(inout) :: self
#if defined HAVE_NETCDF
    integer :: ncerr
    if (self%isopen) then
       write(std_out, *) "Closing lattice history file "//trim(self%filename)//"."
       ncerr=nf90_close(self%ncid)
       NCF_CHECK_MSG(ncerr, "close netcdf lattice history file"//trim(self%filename)//".")
    end if
#endif
  end subroutine finalize

  subroutine write_cell(self, supercell)
    class(lattice_ncfile_t), intent(inout) :: self
    type(mbsupercell_t), intent(in) :: supercell
    integer :: ncerr
#if defined HAVE_NETCDF
    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when starting defining trajectory variables in lattice history file.")
    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%three)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension three in history file.")

    ncerr=nf90_def_dim(self%ncid, "natom", supercell%lattice%natom, self%natom)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension natom in history file.")

    call ab_define_var(self%ncid, [self%three, self%natom], &
         &         self%rvec_id, NF90_INT, "lattice_rvec", "R-vectors for LATTICE ", "dimensionless")

    call ab_define_var(self%ncid, [self%natom], &
         &         self%ilatt_prim_id, NF90_INT, "ilatt_prim", &
         & "index of lattice in primitive cell", "dimensionless")

    
    call ab_define_var(self%ncid, [self%three,self%three], &
         &         self%ref_cell_id, NF90_DOUBLE, "ref_cell", &
         & "REFerence CELL", "bohr")



    call ab_define_var(self%ncid, [self%three,self%natom], &
         &         self%ref_xcart_id, NF90_DOUBLE, "ref_xcart", &
         & "REFerence XCART", "bohr")

    call ab_define_var(self%ncid, [self%natom], &
         &         self%zion_id, NF90_INT, "zion", &
         & "ZION", "dimensionless")

    call ab_define_var(self%ncid, [self%natom], &
         &         self%masses_id, NF90_INT, "masses", &
         & "MASSES", "dimensionless")


    Ncerr=nf90_enddef(self%ncid)

    ! ncerr=nf90_put_var(self%ncid, self%ilatt_prim_id, [supercell%lattice%ilatt_prim], &
    !      &      start=[1], count=[supercell%lattice%natom])
    ! NCF_CHECK_MSG(ncerr, "Error when writting ilatt_prim in lattice history file.")

    ! ncerr=nf90_put_var(self%ncid, self%rvec_id, [supercell%lattice%rvec], &
    !      &      start=[1,1], count=[3, supercell%lattice%natom])
    ! NCF_CHECK_MSG(ncerr, "Error when writting ilatt_prim in lattice history file.")

    
     ncerr=nf90_put_var(self%ncid, self%zion_id, [supercell%lattice%zion], &
          &      start=[1], count=[supercell%lattice%natom])
     NCF_CHECK_MSG(ncerr, "Error when writting zion in lattice history file.")

     ncerr=nf90_put_var(self%ncid, self%masses_id, [supercell%lattice%masses], &
          &      start=[1], count=[supercell%lattice%natom])
     NCF_CHECK_MSG(ncerr, "Error when writting masses in lattice history file.")

     ncerr=nf90_put_var(self%ncid, self%ref_xcart_id, [supercell%lattice%xcart], &
          &      start=[1,1], count=[3,supercell%lattice%natom])
     NCF_CHECK_MSG(ncerr, "Error when writting ref_xcart in lattice history file.")

     ncerr=nf90_put_var(self%ncid, self%ref_cell_id, [supercell%lattice%cell], &
          &      start=[1,1], count=[3,3])
     NCF_CHECK_MSG(ncerr, "Error when writting ref_cell in lattice history file.")

#endif

  end subroutine write_cell

  subroutine def_lattice_var(self)
    class(lattice_ncfile_t), intent(inout) :: self
    !type(lattice_hist_t),intent(in) :: hist
    integer :: ncerr

#if defined HAVE_NETCDF
    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when starting defining trajectory variables in lattice history file.")
    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "ntime", nf90_unlimited, self%ntime)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension ntime in lattice history file.")
    !call ab_define_var(ncid,dim1,typat_id,NF90_DOUBLE,&
    !     &  "typat","types of atoms","dimensionless" )

    !if(self%write_traj==1) then
       call ab_define_var(self%ncid, (/self%three, self%natom, self%ntime /), &
            &         self%xcart_id, NF90_DOUBLE, "xcart", "lattice amplitude", "bohr")
       call ab_define_var(self%ncid, (/self%three, self%natom, self%ntime /), &
            &         self%vcart_id, NF90_DOUBLE, "vcart", "cartesian velocity", "bohr/a.u.")
       call ab_define_var(self%ncid, (/self%three, self%three, self%ntime /), &
            &         self%cell_id, NF90_DOUBLE, "cell", "cell parameters", "bohr")


    !endif

    !call ab_define_var(self%ncid, (/ self%ntime /), &
    !     &         self%time_id, NF90_DOUBLE, "time", "time", "s")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%etotal_id, NF90_DOUBLE, "etotal", "TOTAL energy", "Ha")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%ekin_id, NF90_DOUBLE, "ekin", "Kinetic energy", "Ha")

    !call ab_define_var(self%ncid, (/ self%ntime /), &
    !     &         self%entropy_id, NF90_DOUBLE, "entropy", "Entropy", "Joule/K")

    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%itime_id, NF90_INT, "itime", "index of time in timeline", "1")

    ncerr=nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when finishing defining variables in lattice history file.")
#endif


  end subroutine def_lattice_var



  subroutine write_one_step(self, xcart, vcart, Etotal, Ek)
    class(lattice_ncfile_t), intent(inout) :: self
    !type(lattice_hist_t), intent(in) :: hist
    real(dp) , intent(in) :: Etotal, Ek, xcart(:,:), vcart(:,:)
    integer :: ncerr, itime, natom
    natom=size(xcart, 2)
    self%itime=self%itime+1
    itime = self%itime

#if defined HAVE_NETCDF
    !if(self%write_traj ==1) then
    ncerr=nf90_put_var(self%ncid, self%xcart_id, xcart, &
            &      start=[1,1, itime], count=[3,natom, 1])

    ncerr=nf90_put_var(self%ncid, self%vcart_id, vcart, &
         &      start=[1,1, itime], count=[3,natom, 1])

    NCF_CHECK_MSG(ncerr, "Error when writting lattice amplitudes in lattice history file.")
    !end if

    ncerr=nf90_put_var(self%ncid, self%etotal_id, [Etotal], &
         &      start=[itime], count=[1])
    NCF_CHECK_MSG(ncerr, "Error when writting total energy in lattice history file.")

    ncerr=nf90_put_var(self%ncid, self%ekin_id, [Ek], &
         &      start=[itime], count=[1])
    NCF_CHECK_MSG(ncerr, "Error when writting kinetic energy in lattice history file.")

    ncerr=nf90_put_var(self%ncid, self%itime_id, [self%itime], &
         &      start=[itime], count=[1])
    NCF_CHECK_MSG(ncerr, "Error when writting itime in lattice history file.")

#endif

  end subroutine write_one_step

  subroutine close(self)
    class(lattice_ncfile_t) :: self

  end subroutine close


end module m_lattice_ncfile
