!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sspin_primitive_potential
!! NAME
!! m_slc_primitive_potential
!!
!! FUNCTION
!! This module contains the atomic structures and the spin-lattice coupling terms inside the primitive cell
!! which can be directly mapped to the netcdf file. It is not for the calculation, but for constructing
!! the hamiltonian in supercell.
!!
!! Datatypes:
!!  slc_primitive_potential_t
!!
!! Subroutines:
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu,nehelbig)
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

module m_slc_primitive_potential

  use m_nctk
  !#if defined HAVE_NETCDF
  use netcdf
  !#endif
  use m_abicore
  use defs_basis
  use m_errors

  use m_dynamic_array, only: int2d_array_type
  use m_multibinit_cell, only: mbcell_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_primitive_potential, only: primitive_potential_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_xmpi


  type, public, extends(primitive_potential_t) :: slc_primitive_potential_t
     integer :: natom, nspin            ! on every mpi node
     !two possible coupling terms, quadratic-linear and biquadratic
     !all the following are only on master node 
     type(ndcoo_mat_t) :: oiju          ! parameter values
     type(ndcoo_mat_t) :: tijuv         ! parameter values
     !R vector of i index is always zero
     type(int2d_array_type) :: oRjlist  ! R vector for j index of oiju 
     type(int2d_array_type) :: oRulist  ! R vector for u index of oiju
     type(int2d_array_type) :: tRjlist  ! R vector for j index of tijuv
     type(int2d_array_type) :: tRulist  ! R vector for u index of tijuv
     type(int2d_array_type) :: tRvlist  ! R vector for v index of tijuv

   contains
     procedure:: initialize
     procedure:: finalize
     !procedure :: set_slc_primcell
     procedure :: set_oiju
     procedure :: set_oiju_1term
     !procedure:: set_tijuv
     procedure :: load_from_files
     procedure :: read_netcdf
     !procedure:: fill_supercell
  end type slc_primitive_potential_t

contains

  subroutine initialize(self, primcell)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(mbcell_t), target, intent(inout) :: primcell
    !integer, intent(in) :: nspin
    self%primcell=>primcell
    self%label="SLC_primitive_potential"
    self%has_spin=.True.
    self%has_displacement=.True.
    self%has_strain=.False.
    self%has_lwf=.False.
  end subroutine initialize


  subroutine finalize(self)
    class(slc_primitive_potential_t), intent(inout) :: self
    call self%oiju%finalize()
    call self%oRjlist%finalize()
    call self%oRulist%finalize
    !call self%tijuv%finalize()
    nullify(self%primcell)
    self%nspin=0
    self%natom=0
    self%label="Destroyed SLC_primitive_potential"
    call self%primitive_potential_t%finalize()
  end subroutine finalize

  subroutine load_from_files(self, params, fnames)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    character(len=fnlen) :: ncdf_fname
    character(len=500) :: message
    integer :: ii

    if (xmpi_comm_rank(xmpi_world)==0) then
       ncdf_fname=fnames(4)
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
            &     'reading spin-lattice coupling terms.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
    endif
    call self%read_netcdf(trim(ncdf_fname)//char(0), params%slc_coupling)
  end subroutine load_from_files
 

  subroutine read_netcdf(self, ncdf_fname, slc_coupling)
    class(slc_primitive_potential_t), intent(inout) :: self
    character (len = *), intent(in) :: ncdf_fname
    integer, intent(in) :: slc_coupling
    integer:: ncid, ncerr
    integer:: dimid, varid
    integer:: oiju_ndata, tijuv_ndata
    
    !O_iju: quadratic in spin, linear in lattice variables
    integer, allocatable :: oiju_ilist(:)
    integer, allocatable :: oiju_jlist(:)
    integer, allocatable :: oiju_ulist(:)
    integer, allocatable :: oiju_Rjlist(:,:)
    integer, allocatable :: oiju_Rulist(:,:)
    real(dp), allocatable :: oiju_vallist(:)

    !T_ijuv: quadratic in spin and lattice variables
    integer, allocatable :: tijuv_ilist(:)
    integer, allocatable :: tijuv_jlist(:)
    integer, allocatable :: tijuv_ulist(:)
    integer, allocatable :: tijuv_vlist(:)
    integer, allocatable :: tijuv_Rjlist(:,:)
    integer, allocatable :: tijuv_Rulist(:,:)
    integer, allocatable :: tijuv_Rvlist(:,:)
    real(dp), allocatable :: tijuv_vallist(:)

!#if defined HAVE_NETCDF

    ncerr = nf90_open(ncdf_fname,NF90_NOWRITE,ncid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A24)') 'Could not open netcdf file'
    else 
      write(std_out,'(A55)') 'Reading spin-lattice coupling parameters from netcdf file'
    endif
    
    ! read primcell info
    ncerr=nctk_get_dim(ncid, "natom", self%natom)
    ncerr=nctk_get_dim(ncid, "nspin", self%nspin)

    ! read O_iju terms    
    ncerr = nf90_inq_dimid(ncid, "Oiju_ndata", dimid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A20)') 'Spin-lattice coupling does not contain O_iju term'
    else
      ncerr = nctk_get_dim(ncid, "Oiju_ndata", oiju_ndata)
      ABI_ALLOCATE(oiju_ilist, (oiju_ndata))
      ABI_ALLOCATE(oiju_jlist, (oiju_ndata))
      ABI_ALLOCATE(oiju_ulist, (oiju_ndata))
      ABI_ALLOCATE(oiju_Rjlist, (3, oiju_ndata))
      ABI_ALLOCATE(oiju_Rulist, (3, oiju_ndata))
      ABI_ALLOCATE(oiju_vallist, (oiju_ndata))

      ierr =nf90_inq_varid(ncid, "Oiju_ilist", varid)
      call nc_handle_err(ierr, "Oiju_ilist")
      ierr = nf90_get_var(ncid, varid, oiju_ilist)
      call nc_handle_err(ierr, "Oiju_ilist")

      ierr =nf90_inq_varid(ncid, "Oiju_jlist", varid)
      call nc_handle_err(ierr, "Oiju_jlist")
      ierr = nf90_get_var(ncid, varid, oiju_jlist)
      call nc_handle_err(ierr, "Oiju_jlist")

      ierr =nf90_inq_varid(ncid, "Oiju_ulist", varid)
      call nc_handle_err(ierr, "Oiju_uist")
      ierr = nf90_get_var(ncid, varid, oiju_ulist)
      call nc_handle_err(ierr, "Oiju_ulist")

      ierr =nf90_inq_varid(ncid, "Oiju_Rjlist", varid)
      call nc_handle_err(ierr, "Oiju_Rjlist")
      ierr = nf90_get_var(ncid, varid, oiju_Rjlist)
      call nc_handle_err(ierr, "Oiju_Rjlist")

      ierr =nf90_inq_varid(ncid, "Oiju_Rulist", varid)
      call nc_handle_err(ierr, "Oiju_Rulist")
      ierr = nf90_get_var(ncid, varid, oiju_Rulist)
      call nc_handle_err(ierr, "Oiju_Rulist")

      ierr =nf90_inq_varid(ncid, "Oiju_Rulist", varid)
      call nc_handle_err(ierr, "Oiju_Rulist")
      ierr = nf90_get_var(ncid, varid, oiju_Rulist)
      call nc_handle_err(ierr, "Oiju_Rulist")

      ierr =nf90_inq_varid(ncid, "Oiju_vallist", varid)
      call nc_handle_err(ierr, "Oiju_vallist")
      ierr = nf90_get_var(ncid, varid, oiju_vallist)
      call nc_handle_err(ierr, "Oiju_vallist")

      write(std_out,'(A7,I10,A11)') 'O_iju:', oiju_ndata, 'terms read'  

      !change units from eV to Ha
      oiju_vallist(:) = oiju_vallist(:)*eV_Ha
  
      !fill the sparse matrix for oiju parameters
      call self%set_oiju(oiju_ndata, oiju_ilist, oiju_jlist, oiju_ulist, oiju_Rjlist, oiju_Rulist, oiju_vallist)

      ABI_SFREE(oiju_ilist)
      ABI_SFREE(oiju_jlist)
      ABI_SFREE(oiju_ulist)
      ABI_SFREE(oiju_Rjlist)
      ABI_SFREE(oiju_Rulist)
      ABI_SFREE(oiju_vallist)

    endif


    !if(slc_coupling>1) then
    !  ncerr = nf90_inq_varid(ncid, "Tiju_vallist", tiju_varid)
    !  if(ncerr /= NF90_NOERR) then
    !    write(std_out,'(A30)') 'Could not find values for Tijuv'
    !  else
    !    ABI_ALLOCATE(tijuv, (ndata))
    !    ncerr = nf90_get_var(ncid, tijuv_varid, tijuv)
    !    write(std_out,'(A7,I10,A11)') 'T_ijuv:', ndata, 'terms read'  
    !  endif
    !endif
    
    ncerr = nf90_close(ncid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A25)') 'Could not close netcdf file'
    endif


  end subroutine read_netcdf

  subroutine set_oiju(self, nn, ilist, jlist, ulist, Rjlist, Rulist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: jlist(nn)
    integer,                          intent(in)    :: ulist(nn)
    integer,                          intent(in)    :: Rjlist(3,nn)
    integer,                          intent(in)    :: Rulist(3,nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx

    call self%oiju%initialize(mshape=[-1, -1, self%nspin*3, self%nspin*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_oiju_1term(ilist(idx), jlist(idx), ulist(idx), Rjlist(:,idx), Rulist(:,idx), vallist(idx))
      end do
    endif
  end subroutine set_oiju

  subroutine set_oiju_1term(self, ii, jj, uu, Rj, Ru, val)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ii    
    integer,                          intent(in)    :: jj
    integer,                          intent(in)    :: uu
    integer,                          intent(in)    :: Rj(3)
    integer,                          intent(in)    :: Ru(3)
    real(dp),                         intent(in)    :: val

    integer :: indRj, indRu
    
    call self%oRjlist%push_unique(Rj, position=indRj)
    call self%oRulist%push_unique(Ru, position=indRu)
    
    call self%oiju%add_entry(ind=[indRj, indRu, ii, jj, uu ], val=val)

  end subroutine set_oiju_1term

  subroutine nc_handle_err(ierr, name)
    integer, intent ( in) ::ierr
    character(*), optional, intent(in) :: name
    if(ierr/= nf90_noerr) then
       if (present(name)) then
          write(std_out, *)  trim(nf90_strerror(ierr)), "when trying to read ", name
       else
          write(std_out, *)  trim(nf90_strerror(ierr))
       end if
       stop "Stopped"
    end if
  end subroutine nc_handle_err


end module m_slc_primitive_potential
