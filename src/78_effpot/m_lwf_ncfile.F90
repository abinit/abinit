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
!! Copyright (C) 2001-2025 ABINIT group (hexu)
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
#if defined HAVE_NETCDF
  use netcdf
#endif
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
     integer :: lwf_id, time_id, itime_id, etotal_id, entropy_id, displacement_id
     integer :: ilwf_prim_id, rvec_id, lwf_masses_id, natom_id
     integer :: itime, ntime,ncell
     ! itime: time index
     integer :: write_traj=1
     !whether to write the trajectory
     character(len=fnlen) :: filename

     ! netcdf filename
     logical :: has_constrain=.False.
     integer :: n_fixed_lwf
     integer, allocatable :: fixed_lwf_ids(:)
     real(dp), allocatable :: fixed_lwf_values(:)
   contains
     ! initialize
     procedure :: initialize
     procedure :: finalize
     procedure :: def_lwf_var
     procedure :: write_cell
     procedure :: write_one_step
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
#if defined HAVE_NETCDF
    write(std_out,*) "Write iteration in lwf history file "//trim(self%filename)//"."
    !  Create netCDF file
    ncerr = nf90_create(path=trim(self%filename), cmode=NF90_NETCDF4, ncid=self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when creating netcdf history file")
    self%isopen=.True.
    ncerr =nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when ending def mode in lwf netcdf history file")
#endif
  end subroutine initialize

  subroutine finalize(self)
    class(lwf_ncfile_t), intent(inout) :: self
#if defined HAVE_NETCDF
    integer :: ncerr
    if (self%isopen) then
       write(std_out, *) "Closing lwf history file "//trim(self%filename)//"."
       ncerr=nf90_close(self%ncid)
       NCF_CHECK_MSG(ncerr, "close netcdf lwf history file"//trim(self%filename)//".")
    end if
#endif
  end subroutine finalize

  subroutine write_cell(self, supercell)
    class(lwf_ncfile_t), intent(inout) :: self
    type(mbsupercell_t), intent(in) :: supercell
    integer :: ncerr, natom3, nnz, id_nnz,  id_natom3, id_map_ilwf, id_map_idisp, id_map_val
    integer :: latt_rvec_id,  ilatt_prim_id, ref_cell_id, ref_xcart_id, zion_id, masses_id
#if defined HAVE_NETCDF
    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when starting defining trajectory variables in lwf history file.")
    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "three", 3, self%three)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension three in history file.")

    ncerr=nf90_def_dim(self%ncid, "nlwf", supercell%lwf%nlwf, self%nlwf)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nlwf in history file.")


    ncerr=nf90_def_dim(self%ncid, "natom", supercell%lattice%natom, self%natom_id)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension natom in history file.")


    natom3=supercell%lwf%lwf_latt_coeffs%coeffs%mshape(1)
    ncerr=nf90_def_dim(self%ncid, "natom3", natom3, id_natom3)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension natom3 in history file.")

    nnz=supercell%lwf%lwf_latt_coeffs%coeffs%nnz
    ncerr=nf90_def_dim(self%ncid, "lwf_latt_map_nnz",nnz, id_nnz)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension nlwf in history file.")

    !ncerr=nf90_def_dim(self%ncid, "nR", 3, self%nlwf)
    !NCF_CHECK_MSG(ncerr, "Error when defining dimension three in history file.")

    call ab_define_var(self%ncid, [self%three, self%nlwf], &
         &         self%rvec_id, NF90_INT, "lwf_rvec", "R-vectors for LWF ", "dimensionless")

    call ab_define_var(self%ncid, [self%nlwf], &
         &         self%ilwf_prim_id, NF90_INT, "ilwf_prim", &
         & "index of lwf in primitive cell", "dimensionless")

    call ab_define_var(self%ncid, [self%nlwf], &
         &         self%lwf_masses_id, NF90_DOUBLE, "lwf_masses", "LWF MASSES", "dimensionless")


    ! Lattice
    call ab_define_var(self%ncid, [self%three, self%natom_id], &
         &         latt_rvec_id, NF90_INT, "lattice_rvec", "R-vectors for LATTICE ", "dimensionless")

    call ab_define_var(self%ncid, [self%natom_id], &
         &         ilatt_prim_id, NF90_INT, "ilatt_prim", &
         & "index of lattice in primitive cell", "dimensionless")


    call ab_define_var(self%ncid, [self%three,self%three], &
         &         ref_cell_id, NF90_DOUBLE, "ref_cell", &
         & "REFerence CELL", "bohr")

    call ab_define_var(self%ncid, [self%three,self%natom_id], &
         &         ref_xcart_id, NF90_DOUBLE, "ref_xcart", &
         & "REFerence XCART", "bohr")

    call ab_define_var(self%ncid, [self%natom_id], &
         &         zion_id, NF90_INT, "zion", &
         & "ZION", "dimensionless")

    call ab_define_var(self%ncid, [self%natom_id], &
         &         masses_id, NF90_DOUBLE, "masses", &
         & "MASSES", "dimensionless")



    ! define vars for lwf lattice displacement mapping in the format of a COO matrix.
    call ab_define_var(self%ncid, [id_nnz], id_map_idisp, &
         & NF90_INT, "lwf_latt_map_id_displacement", &
         & "LWF lattice mapping coefficient COO matrix displacement id",  "dimensionless")

    ncerr=nf90_def_var_deflate(self%ncid, id_map_idisp, shuffle=1, deflate=1, deflate_level=2)
    NCF_CHECK_MSG(ncerr, "Error when defining deflating for variable id_map_idisp")



    call ab_define_var(self%ncid, [id_nnz], id_map_ilwf,&
         & NF90_INT, "lwf_latt_map_id_lwf", &
         & "LWF lattice mapping coefficient COO matrix LWF id","dimensionless")

    ncerr=nf90_def_var_deflate(self%ncid, id_map_ilwf, shuffle=1, deflate=1, deflate_level=2)
    NCF_CHECK_MSG(ncerr, "Error when defining delfating for variable id_map_ilwf")


    call ab_define_var(self%ncid, [id_nnz], id_map_val, &
         & NF90_DOUBLE, "lwf_latt_map_values", &
         & "LWF lattice mapping coefficient COO matrix values","dimensionless")

    ncerr=nf90_def_var_deflate(self%ncid, id_map_val, shuffle=1, deflate=1, deflate_level=2)
    NCF_CHECK_MSG(ncerr, "Error when defining delfating for variable id_map_val")

    ncerr=nf90_enddef(self%ncid)



    ncerr=nf90_put_var(self%ncid, zion_id, [supercell%lattice%zion], &
         &      start=[1], count=[supercell%lattice%natom])
    NCF_CHECK_MSG(ncerr, "Error when writting zion in lattice history file.")

    ncerr=nf90_put_var(self%ncid, masses_id, [supercell%lattice%masses], &
         &      start=[1], count=[supercell%lattice%natom])
    NCF_CHECK_MSG(ncerr, "Error when writting masses in lattice history file.")

    ncerr=nf90_put_var(self%ncid, ref_xcart_id, [supercell%lattice%xcart], &
         &      start=[1,1], count=[3,supercell%lattice%natom])
    NCF_CHECK_MSG(ncerr, "Error when writting ref_xcart in lattice history file.")

    ncerr=nf90_put_var(self%ncid, ref_cell_id, [supercell%lattice%cell], &
         &      start=[1,1], count=[3,3])
    NCF_CHECK_MSG(ncerr, "Error when writting ref_cell in lattice history file.")



    ncerr=nf90_put_var(self%ncid, self%ilwf_prim_id, [supercell%lwf%ilwf_prim], &
         &      start=[1], count=[supercell%lwf%nlwf])
    NCF_CHECK_MSG(ncerr, "Error when writting ilwf_prim in lwf history file.")

    ncerr=nf90_put_var(self%ncid, self%rvec_id, [supercell%lwf%rvec], &
         &      start=[1,1], count=[3, supercell%lwf%nlwf])
    NCF_CHECK_MSG(ncerr, "Error when writting lwf rvec in lwf history file.")

    ncerr=nf90_put_var(self%ncid, self%lwf_masses_id, [supercell%lwf%lwf_masses], &
         &      start=[1], count=[supercell%lwf%nlwf])
    NCF_CHECK_MSG(ncerr, "Error when writting lwf_masses in lwf history file.")

    ncerr=nf90_put_var(self%ncid, id_map_idisp,  &
         &         supercell%lwf%lwf_latt_coeffs%coeffs%ind%data(1, 1:nnz), &
         &      start=[1], count=[nnz])
    NCF_CHECK_MSG(ncerr, "Error when writting id_map_idisp in lwf history file.")


    ncerr=nf90_put_var(self%ncid, id_map_ilwf,  &
         &         supercell%lwf%lwf_latt_coeffs%coeffs%ind%data(2, 1:nnz), &
         &      start=[1], count=[nnz])
    NCF_CHECK_MSG(ncerr, "Error when writting id_map_ilwf in lwf history file.")


    ncerr=nf90_put_var(self%ncid, id_map_val,  &
         &         supercell%lwf%lwf_latt_coeffs%coeffs%val%data(1:nnz), &
         &      start=[1], count=[nnz])
    NCF_CHECK_MSG(ncerr, "Error when writting id_map_ilwf in lwf history file.")


#endif

  end subroutine write_cell

  subroutine def_lwf_var(self, hist)
    class(lwf_ncfile_t), intent(inout) :: self
    type(lwf_hist_t),intent(in) :: hist
    integer :: ncerr
    ABI_UNUSED_A(hist)

#if defined HAVE_NETCDF
    ncerr = nf90_redef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when starting defining trajectory variables in lwf history file.")
    ! define dimensions
    ncerr=nf90_def_dim(self%ncid, "ntime", nf90_unlimited, self%ntime)
    NCF_CHECK_MSG(ncerr, "Error when defining dimension ntime in lwf history file.")
    !call ab_define_var(ncid,dim1,typat_id,NF90_DOUBLE,&
    !     &  "typat","types of atoms","dimensionless" )

    !if(self%write_traj==1) then
       call ab_define_var(self%ncid, (/ self%nlwf, self%ntime /), &
            &         self%lwf_id, NF90_DOUBLE, "lwf", "lwf amplitude", "dimensionless")

        ncerr=nf90_def_var_deflate(self%ncid, self%lwf_id, shuffle=1, deflate=1, deflate_level=2)
        NCF_CHECK_MSG(ncerr, "Error when defining deflating for variable lwf")

    !endif

    !call ab_define_var(self%ncid, (/ self%ntime /), &
    !     &         self%time_id, NF90_DOUBLE, "time", "time", "s")
    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%etotal_id, NF90_DOUBLE, "etotal", "TOTAL energy", "Ha")
    !call ab_define_var(self%ncid, (/ self%ntime /), &
    !     &         self%entropy_id, NF90_DOUBLE, "entropy", "Entropy", "Joule/K")

    call ab_define_var(self%ncid, (/ self%ntime /), &
         &         self%itime_id, NF90_INT, "itime", "index of time in timeline", "1")

    ncerr=nf90_enddef(self%ncid)
    NCF_CHECK_MSG(ncerr, "Error when finishing defining variables in lwf history file.")
#endif

  end subroutine def_lwf_var



  subroutine write_one_step(self, hist)
    class(lwf_ncfile_t), intent(inout) :: self
    type(lwf_hist_t), intent(in) :: hist
    integer :: ncerr, itime
    self%itime=self%itime+1
    itime = self%itime

#if defined HAVE_NETCDF
    !if(self%write_traj ==1) then
    !print *, "Write one step of netcdf."
    ncerr=nf90_put_var(self%ncid, self%lwf_id, hist%current_lwf, &
            &      start=[1, itime], count=[hist%nlwf, 1])
    NCF_CHECK_MSG(ncerr, "Error when writting lwf amplitudes in lwf history file.")
    !end if

    ncerr=nf90_put_var(self%ncid, self%etotal_id, [hist%current_energy], &
         &      start=[itime], count=[1])
    NCF_CHECK_MSG(ncerr, "Error when writting energy in lwf history file.")

    ncerr=nf90_put_var(self%ncid, self%itime_id, [self%itime], &
         &      start=[itime], count=[1])
    NCF_CHECK_MSG(ncerr, "Error when writting energy in lwf history file.")

    !ncerr=nf90_put_var(self%ncid, self%displacement_id, [displacement], &
    !     &      start=[1,itime], count=[3, hist%natom, 1])
    !NCF_CHECK_MSG(ncerr, "Error when writting energy in lwf history file.")

#endif

  end subroutine write_one_step


end module m_lwf_ncfile
