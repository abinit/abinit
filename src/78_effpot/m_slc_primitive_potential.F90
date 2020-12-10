!!****m* ABINIT/m_slc_primitive_potential
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
!! Copyright (C) 2001-2020 ABINIT group (hexu,nehelbig)
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
#if defined HAVE_NETCDF
  use netcdf
#endif
  use m_abicore
  use defs_basis
  use m_errors

  use m_slc_potential
  use m_spmat_ndcoo

  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_cell, only: mbcell_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_primitive_potential, only: primitive_potential_t

  use m_supercell_maker, only: supercell_maker_t
  use m_xmpi

  implicit none

  type, public, extends(primitive_potential_t) :: slc_primitive_potential_t
     integer :: natom, nspin        ! on every mpi node
     logical :: has_bilin=.False.   ! bilinear coupling term, i.e. liu        
     logical :: has_linquad=.False. ! spin first then lattice, i.e. niuv
     logical :: has_quadlin=.False. ! spin first then lattice, i.e. oiju
     logical :: has_biquad=.False.  ! biquadratic coupling term, i.e. tijuv

     !four possible coupling terms, bilinear, quadratic-linear, linear-quadratic, and biquadratic
     !all the following are only on master node 
     type(ndcoo_mat_t) :: liu           ! parameter values
     type(ndcoo_mat_t) :: niuv          ! parameter values
     type(ndcoo_mat_t) :: oiju          ! parameter values
     type(ndcoo_mat_t) :: tijuv         ! parameter values
     !R vector of i index is always zero
     type(int2d_array_type) :: lRulist  ! R vector for u index of liu 
     type(int2d_array_type) :: nRulist  ! R vector for u index of niuv 
     type(int2d_array_type) :: nRvlist  ! R vector for v index of niuv 
     type(int2d_array_type) :: oRjlist  ! R vector for j index of oiju 
     type(int2d_array_type) :: oRulist  ! R vector for u index of oiju
     type(int2d_array_type) :: tRjlist  ! R vector for j index of tijuv
     type(int2d_array_type) :: tRulist  ! R vector for u index of tijuv
     type(int2d_array_type) :: tRvlist  ! R vector for v index of tijuv

   contains
     procedure:: initialize
     procedure:: finalize
     procedure :: set_liu
     procedure :: set_liu_1term
     procedure :: set_niuv
     procedure :: set_niuv_1term
     procedure :: set_oiju
     procedure :: set_oiju_1term
     procedure :: set_tijuv
     procedure :: set_tijuv_1term
     procedure :: load_from_files
     procedure :: read_netcdf
     procedure :: read_liu
     procedure :: read_niuv
     procedure :: read_oiju
     procedure :: read_tijuv
     procedure :: fill_supercell
     procedure :: fill_liu
     procedure :: fill_niuv
     procedure :: fill_oiju
     procedure :: fill_tijuv
     procedure :: set_liu_sc
     procedure :: set_niuv_sc
     procedure :: set_oiju_sc
     procedure :: set_tijuv_sc
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
    call self%oRulist%finalize()
    call self%tijuv%finalize()
    call self%tRjlist%finalize()
    call self%tRulist%finalize()
    call self%tRvlist%finalize()
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

    ABI_UNUSED_A(params)

    if (xmpi_comm_rank(xmpi_world)==0) then
       ncdf_fname=fnames(3)
       write(message,'(a,(81a, 80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
            &     '- Reading spin-lattice coupling terms from ', trim(ncdf_fname)
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
    endif
    call self%read_netcdf(trim(ncdf_fname)//char(0))
    
    ABI_UNUSED_A(params)
  end subroutine load_from_files

  !-----------------------------------
  ! reading parameters from ncdf file
  !-----------------------------------
  subroutine read_netcdf(self, ncdf_fname)
    class(slc_primitive_potential_t), intent(inout) :: self
    character (len = *),              intent(in)    :: ncdf_fname

    integer:: ncid, ncerr, comm

#if defined HAVE_NETCDF

    comm = xmpi_world

    ncerr = nctk_open_read(ncid, ncdf_fname, comm)

    ! read primcell info
    ncerr=nctk_get_dim(ncid, "natom", self%natom)
    ncerr=nctk_get_dim(ncid, "nspin", self%nspin)

    ! read different coupling terms if they are present in the file
    call self%read_liu(ncid)
    call self%read_niuv(ncid)
    call self%read_oiju(ncid)
    call self%read_tijuv(ncid)
    
    ncerr = nf90_close(ncid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A25)') 'Could not close netcdf file'
    endif
#else
    ABI_ERROR("Multibint should be installed with netcdf enabled to run this.")

#endif

  end subroutine read_netcdf

  !---------------------------------------
  ! reading liu parameters from ncdf file
  ! TODO: test
  !---------------------------------------
  subroutine read_liu(self, ncid)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ncid    

    integer :: ncerr, dimid, ndata, varid
    integer, allocatable :: ilist(:), ulist(:,:)
    real(dp), allocatable :: vallist(:)
#if defined HAVE_NETCDF

    ncerr = nf90_inq_dimid(ncid, "spin_lattice_Liu_number_of_entries", dimid)
    if(ncerr /= NF90_NOERR) then
      ndata = 0
    else
      self%has_bilin=.True.
      ncerr = nctk_get_dim(ncid, "spin_lattice_Liu_number_of_entries", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(ulist, (4,ndata))
      ABI_ALLOCATE(vallist, (ndata))

      varid = nctk_idname(ncid, "spin_lattice_Liu_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call netcdf_check(ncerr, "when reading Liu_ilist") 

      varid = nctk_idname(ncid, "spin_lattice_Liu_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call netcdf_check(ncerr, "when reading Liu_ulist")

      varid = nctk_idname(ncid, "spin_lattice_Liu_valuelist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call netcdf_check(ncerr, "when reading Liu_valuelist")

      !change units from eV to Ha and Ang to Bohr
      vallist(:) = vallist(:)*eV_Ha*Bohr_Ang
  
      !fill the sparse matrix for liu parameters
      call self%set_liu(ndata, ilist, ulist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(ulist)
      ABI_SFREE(vallist)
    endif

    write(std_out,'(A8,I10,A11)') 'L_iu:  ', ndata, 'terms read'  
#else
    ABI_ERROR("Multibinit should be installed with netcdf.")
#endif

  end subroutine read_liu

  !----------------------------------------
  ! reading niuv parameters from ncdf file
  ! TODO: test
  !----------------------------------------
  subroutine read_niuv(self, ncid)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ncid    

    integer :: ncerr, dimid, ndata, varid
    integer, allocatable :: ilist(:), ulist(:,:), vlist(:,:)
    real(dp), allocatable :: vallist(:)
#if defined HAVE_NETCDF
    ncerr = nf90_inq_dimid(ncid, "spin_lattice_Niuv_number_of_entries", dimid)
    if(ncerr /= NF90_NOERR) then
      ndata = 0
    else
      self%has_linquad=.True.
      ncerr = nctk_get_dim(ncid, "spin_lattice_Niuv_number_of_entries", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(ulist, (4,ndata))
      ABI_ALLOCATE(vlist, (4,ndata))
      ABI_ALLOCATE(vallist, (ndata))

      varid = nctk_idname(ncid, "spin_lattice_Niuv_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call netcdf_check(ncerr, "when reading Niuv_ilist") 

      varid = nctk_idname(ncid, "spin_lattice_Niuv_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call netcdf_check(ncerr, "when reading Niuv_ulist")

      varid = nctk_idname(ncid, "spin_lattice_Niuv_vlist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call netcdf_check(ncerr, "when reading Niuv_vlist")

      varid = nctk_idname(ncid, "spin_lattice_Niuv_valuelist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call netcdf_check(ncerr, "when reading Niuv_valuelist")

      !change units from eV to Ha and Ang to Bohr
      vallist(:) = vallist(:)*eV_Ha*(Bohr_Ang*Bohr_Ang)
  
      !fill the sparse matrix for liu parameters
      call self%set_niuv(ndata, ilist, ulist, vlist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(ulist)
      ABI_SFREE(vlist)
      ABI_SFREE(vallist)
    endif

    write(std_out,'(A8,I10,A11)') 'N_iuv: ', ndata, 'terms read'  
#else
    ABI_ERROR("Multibinit should be installed with netcdf") 
#endif

  end subroutine read_niuv

  !---------------------------------------
  ! reading oiju parameters from ncdf file
  !---------------------------------------
  subroutine read_oiju(self, ncid)  
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ncid    

    integer :: ncerr, dimid, ndata, varid
    integer, allocatable :: ilist(:), jlist(:,:), ulist(:,:)
    real(dp), allocatable :: vallist(:)

#if defined HAVE_NETCDF
    ncerr = nf90_inq_dimid(ncid, "spin_lattice_Oiju_number_of_entries", dimid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A20)') 'No O_iju term found'
    else
      self%has_quadlin=.True.
      ncerr = nctk_get_dim(ncid, "spin_lattice_Oiju_number_of_entries", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(jlist, (4,ndata))
      ABI_ALLOCATE(ulist, (4,ndata))
      ABI_ALLOCATE(vallist, (ndata))

      varid = nctk_idname(ncid, "spin_lattice_Oiju_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call netcdf_check(ncerr, "when reading Oiju_ilist") 

      varid = nctk_idname(ncid, "spin_lattice_Oiju_jlist")
      ncerr = nf90_get_var(ncid, varid, jlist)
      call netcdf_check(ncerr, "when reading Oiju_jlist")

      varid = nctk_idname(ncid, "spin_lattice_Oiju_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call netcdf_check(ncerr, "when reading Oiju_ulist")

      varid = nctk_idname(ncid, "spin_lattice_Oiju_valuelist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call netcdf_check(ncerr, "when reading spin_lattice_Oiju_valuelist")

      write(std_out,'(A8,I10,A11)') 'O_iju: ', ndata, 'terms read'  

      !change units from eV to Ha and Ang to Bohr
      vallist(:) = vallist(:)*eV_Ha*Bohr_Ang

      !fill the sparse matrix for oiju parameters
      call self%set_oiju(ndata, ilist, jlist, ulist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(jlist)
      ABI_SFREE(ulist)
      ABI_SFREE(vallist)
    endif
#else
    ABI_ERROR('Multibinit should be install with netcdf to run this.')
#endif

  end subroutine read_oiju

  !-----------------------------------------
  ! reading tijuv parameters from ncdf file
  !-----------------------------------------
  subroutine read_tijuv(self, ncid)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ncid    

    integer :: ncerr, dimid, ndata, varid
    integer, allocatable :: ilist(:), jlist(:,:), ulist(:,:), vlist(:,:)
    real(dp), allocatable :: vallist(:)

#if defined HAVE_NETCDF
    ncerr = nf90_inq_dimid(ncid, "spin_lattice_Tijuv_number_of_entries", dimid)
    if(ncerr /= NF90_NOERR) then
      ndata = 0
    else
      self%has_biquad=.True.
      ncerr = nctk_get_dim(ncid, "spin_lattice_Tijuv_number_of_entries", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(jlist, (4, ndata))
      ABI_ALLOCATE(ulist, (4, ndata))
      ABI_ALLOCATE(vlist, (4, ndata))
      ABI_ALLOCATE(vallist, (ndata))

      varid = nctk_idname(ncid, "spin_lattice_Tijuv_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call netcdf_check(ncerr, "when reading Tijuv_ilist") 

      varid = nctk_idname(ncid, "spin_lattice_Tijuv_jlist")
      ncerr = nf90_get_var(ncid, varid, jlist)
      call netcdf_check(ncerr, "when reading Tijuv_jlist")

      varid = nctk_idname(ncid, "spin_lattice_Tijuv_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call netcdf_check(ncerr, "when reading Tijuv_ulist")

      varid = nctk_idname(ncid, "spin_lattice_Tijuv_vlist")
      ncerr = nf90_get_var(ncid, varid, vlist)
      call netcdf_check(ncerr, "when reading Tijuv_vlist")

      varid = nctk_idname(ncid, "spin_lattice_Tijuv_valuelist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call netcdf_check(ncerr, "when reading Tijuv_valuelist")

      write(std_out,'(A8,I10,A11)') 'T_ijuv:', ndata, 'terms read'  

      !change units from eV to Ha and Ang to Bohr
      vallist(:) = vallist(:)*eV_Ha*(Bohr_Ang*Bohr_Ang)
  
      !fill the sparse matrix for tijuv parameters
      call self%set_tijuv(ndata, ilist, jlist, ulist, vlist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(jlist)
      ABI_SFREE(ulist)
      ABI_SFREE(vlist)
      ABI_SFREE(vallist)
    endif

#else
    ABI_ERROR('Multibinit should be install with netcdf to run this.')
#endif


  end subroutine read_tijuv

  !---------------------------------------
  ! store liu parameters in sparse matrix
  ! TODO: test
  !---------------------------------------
  subroutine set_liu(self, nn, ilist, ulist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: ulist(4,nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx

    call self%liu%initialize(mshape=[-1, self%nspin*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_liu_1term(ilist(idx), ulist(1,idx), ulist(2:4,idx), vallist(idx))
      end do
    endif
  end subroutine set_liu

  !------------------------------------
  ! add one entry to sparse liu matrix
  ! TODO: test
  !------------------------------------
  subroutine set_liu_1term(self, ii, uu, Ru, val)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ii    
    integer,                          intent(in)    :: uu
    integer,                          intent(in)    :: Ru(3)
    real(dp),                         intent(in)    :: val

    integer :: indRu
    
    call self%lRulist%push_unique(Ru, position=indRu)
    
    call self%liu%add_entry(ind=[indRu, ii, uu], val=val)

  end subroutine set_liu_1term

  !----------------------------------------
  ! store niuv parameters in sparse matrix
  ! TODO: test
  !----------------------------------------
  subroutine set_niuv(self, nn, ilist, ulist, vlist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: ulist(4,nn)
    integer,                          intent(in)    :: vlist(4,nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx
   
    call self%niuv%initialize(mshape=[-1, -1, self%nspin*3, self%natom*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_niuv_1term(ilist(idx), ulist(1,idx), vlist(1,idx), ulist(2:4,idx), vlist(2:4,idx), vallist(idx))
      end do
    endif

  end subroutine set_niuv

  !-------------------------------------
  ! add one entry to sparse niuv matrix
  !-------------------------------------
  subroutine set_niuv_1term(self, ii, uu, vv, Ru, Rv, val)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ii    
    integer,                          intent(in)    :: uu
    integer,                          intent(in)    :: vv
    integer,                          intent(in)    :: Ru(3)
    integer,                          intent(in)    :: Rv(3)
    real(dp),                         intent(in)    :: val

    integer :: indRu, indRv
    
    call self%nRulist%push_unique(Ru, position=indRu)
    call self%nRvlist%push_unique(Rv, position=indRv)
    
    call self%niuv%add_entry(ind=[indRu, indRv, ii, uu, vv], val=val)

  end subroutine set_niuv_1term


  !----------------------------------------
  ! store oiju parameters in sparse matrix
  ! TODO: test
  !----------------------------------------
  subroutine set_oiju(self, nn, ilist, jlist, ulist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: jlist(4,nn)
    integer,                          intent(in)    :: ulist(4,nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx
   
    call self%oiju%initialize(mshape=[-1, -1, self%nspin*3, self%nspin*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_oiju_1term(ilist(idx), jlist(1,idx), ulist(1,idx), jlist(2:4,idx), ulist(2:4,idx), vallist(idx))
      end do
    endif

    call self%oiju%sum_duplicates()
  end subroutine set_oiju

  !-------------------------------------
  ! add one entry to sparse oiju matrix
  ! TODO: test
  !-------------------------------------
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
    
    call self%oiju%add_entry(ind=[indRj, indRu, ii, jj, uu], val=val)

  end subroutine set_oiju_1term

  !----------------------------------------
  ! store tijuv parameters in sparse matrix
  !----------------------------------------
  subroutine set_tijuv(self, nn, ilist, jlist, ulist, vlist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: jlist(4, nn)
    integer,                          intent(in)    :: ulist(4, nn)
    integer,                          intent(in)    :: vlist(4, nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx
   
    call self%tijuv%initialize(mshape=[-1, -1, -1, self%nspin*3, self%nspin*3, self%natom*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_tijuv_1term(ilist(idx), jlist(1, idx), ulist(1, idx), vlist(1, idx), & 
                                & jlist(2:4,idx), ulist(2:4,idx), vlist(2:4,idx), vallist(idx))
      end do
    endif

    call self%tijuv%sum_duplicates()

  end subroutine set_tijuv

  !-------------------------------------
  ! add one entry to sparse tijuv matrix
  !-------------------------------------
  subroutine set_tijuv_1term(self, ii, jj, uu, vv, Rj, Ru, Rv, val)
    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(in)    :: ii    
    integer,                          intent(in)    :: jj
    integer,                          intent(in)    :: uu
    integer,                          intent(in)    :: vv
    integer,                          intent(in)    :: Rj(3)
    integer,                          intent(in)    :: Ru(3)
    integer,                          intent(in)    :: Rv(3)
    real(dp),                         intent(in)    :: val

    integer :: indRj, indRu, indRv
    
    call self%tRjlist%push_unique(Rj, position=indRj)
    call self%tRulist%push_unique(Ru, position=indRu)
    call self%tRvlist%push_unique(Rv, position=indRv)
    
    call self%tijuv%add_entry(ind=[indRj, indRu, indRv, ii, jj, uu, vv], val=val)

  end subroutine set_tijuv_1term


  !-----------------------------------------------------------------
  ! transfer parameter information from primitive cell to supercell
  ! TODO: test
  !-----------------------------------------------------------------
  subroutine fill_supercell(self, scmaker, params, scpot)
    class(slc_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    type(multibinit_dtset_type),       intent(inout) :: params
    class(abstract_potential_t), pointer, intent(inout) :: scpot

    integer :: nspin, sc_nspin, natom, sc_natom
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    nspin=self%nspin
    natom=self%natom
    sc_nspin=nspin * scmaker%ncells
    sc_natom=natom * scmaker%ncells
    
    call xmpi_bcast(sc_nspin, master, comm, ierr)
    call xmpi_bcast(sc_natom, master, comm, ierr)
    ABI_DATATYPE_ALLOCATE_SCALAR(slc_potential_t, scpot)

    select type(scpot) ! use select type because properties only defined for slc_potential are used
    type is (slc_potential_t) 
      call scpot%initialize(sc_nspin, sc_natom)
      call scpot%set_params(params)
      ! fill different coupling terms
      if (iam_master) then
        call self%fill_liu(scpot, scmaker)
        call self%fill_oiju(scpot, scmaker)
        call self%fill_niuv(scpot, scmaker)
        call self%fill_tijuv(scpot, scmaker)

        !Write information which terms are used
        write(std_out,'(A55)') 'Using the following terms for the spin-lattice coupling'
        if(scpot%has_bilin)   write(std_out,'(A19)') 'Bilinear term: Liu'
        if(scpot%has_quadlin) write(std_out,'(A28)') 'Quadratic-linear term: Oiju'
        if(scpot%has_linquad) write(std_out,'(A28)') 'Linear-quadratic term: Niuv'
        if(scpot%has_biquad)  write(std_out,'(A24)') 'Biquadratic term: Tijuv'
      endif
    end select
  end subroutine fill_supercell

  !--------------------------------------------------------
  ! Check for each term if it is needed in the calculation
  ! and present in the ncdf input then put into supercell
  !--------------------------------------------------------
  subroutine fill_liu(self, scpot, scmaker)
    class(slc_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    class(slc_potential_t),            intent(inout) :: scpot

    if(scpot%has_bilin) then
      if(self%has_bilin) then !did we actually find this in the netcdf file
        call self%liu%sum_duplicates()
        call self%set_liu_sc(scpot, scmaker)
      else
        ABI_ERROR("No parameters for bilinear coupling available. Check your input and parameter files.")
        scpot%has_bilin = .False.
      endif
    endif
  end subroutine fill_liu
  
  subroutine fill_niuv(self, scpot, scmaker)
    class(slc_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    class(slc_potential_t),            intent(inout) :: scpot

    if(scpot%has_linquad) then
      if(self%has_linquad) then
        call self%niuv%sum_duplicates()
        call self%set_niuv_sc(scpot, scmaker)
      else
        ABI_ERROR("No parameters for linear-quadratic coupling available. Check your input and parameter files.")
        scpot%has_linquad = .False.
      endif
    endif
  end subroutine fill_niuv

  subroutine fill_oiju(self, scpot, scmaker)
    class(slc_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    class(slc_potential_t),            intent(inout) :: scpot

    if(scpot%has_quadlin) then
      if(self%has_quadlin) then
        call self%oiju%sum_duplicates()
        call self%set_oiju_sc(scpot, scmaker)
       else
         ABI_ERROR("No parameters for quadratic-linear coupling available. Check your input and parameter files.")
         scpot%has_quadlin = .False.
       endif
     endif
  end subroutine fill_oiju

  subroutine fill_tijuv(self, scpot, scmaker)
    class(slc_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t),           intent(inout) :: scmaker
    class(slc_potential_t),            intent(inout) :: scpot

    if(scpot%has_biquad) then
      if(self%has_biquad) then
        call self%tijuv%sum_duplicates()
        call self%set_tijuv_sc(scpot, scmaker)
      else
        ABI_ERROR("No parameters for biquadratic coupling available. Check your input and parameter files.")
        scpot%has_biquad = .False.
      endif
    endif
  end subroutine fill_tijuv

  !------------------------------
  ! fill liu terms in supercell
  ! TODO: test
  !------------------------------
  subroutine set_liu_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Ru_prim(3), liu_ind(3), iRu, i_prim, u_prim, inz
    integer :: ngroup
    integer, allocatable :: i_sc(:), u_sc(:), Ru_sc(:,:)
    integer, allocatable :: i1list(:), ise(:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%liu_sc%initialize(mshape=[self%nspin*3, self%natom*3])
    endif

    do inz=1, self%liu%nnz
      liu_ind=self%liu%get_ind_inz(inz)
      iRu=liu_ind(1)
      i_prim=liu_ind(2)
      u_prim=liu_ind(3)
      Ru_prim=self%lRulist%data(:,iRu)
      call scmaker%trans_i(nbasis=self%nspin*3, i=i_prim, i_sc=i_sc) 
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=u_prim, Rj=Ru_prim, j_sc=u_sc, Rj_sc=Ru_sc)
      val_sc(:)= self%liu%val%data(inz)
      do icell=1, scmaker%ncells
        call scpot%add_liu_term(i_sc(icell), u_sc(icell), val_sc(icell))
      end do
      ABI_SFREE(i_sc)
      ABI_SFREE(u_sc)
      ABI_SFREE(Ru_sc)
    end do
    
    call scpot%liu_sc%group_by_1dim(ngroup, i1list, ise)
    ABI_SFREE(i1list)
    ABI_SFREE(ise)

  end subroutine set_liu_sc

  !------------------------------
  ! fill niuv terms in supercell
  ! TODO: test
  !------------------------------
  subroutine set_niuv_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Ru_prim(3), Rv_prim(3), niuv_ind(5), iRu, iRv, i_prim, u_prim, v_prim, inz
    integer :: ngroup
    integer, allocatable :: i_sc(:), u_sc(:), v_sc(:), Ru_sc(:,:), Rv_sc(:,:)
    integer, allocatable :: i1list(:), ise(:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%niuv_sc%initialize(mshape=[self%nspin*3, self%natom*3, self%natom*3])
    endif

    do inz=1, self%niuv%nnz
      niuv_ind=self%niuv%get_ind_inz(inz)
      iRu=niuv_ind(1)
      iRv=niuv_ind(2)
      i_prim=niuv_ind(3)
      u_prim=niuv_ind(4)
      v_prim=niuv_ind(5)
      Ru_prim=self%nRulist%data(:,iRu)
      Rv_prim=self%nRvlist%data(:,iRv)
      call scmaker%trans_i(nbasis=self%nspin*3, i=i_prim, i_sc=i_sc) 
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=u_prim, Rj=Ru_prim, j_sc=u_sc, Rj_sc=Ru_sc)
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=v_prim, Rj=Rv_prim, j_sc=v_sc, Rj_sc=Rv_sc)
      val_sc(:)= self%niuv%val%data(inz)
      do icell=1, scmaker%ncells
        !call scpot%add_niuv_term(i_sc(icell), u_sc(icell), v_sc(icell), val_sc(icell))
      end do
      ABI_SFREE(i_sc)
      ABI_SFREE(u_sc)
      ABI_SFREE(v_sc)
      ABI_SFREE(Ru_sc)
      ABI_SFREE(Rv_sc)
    end do
    
    call scpot%niuv_sc%group_by_1dim(ngroup, i1list, ise)
    ABI_SFREE(i1list)
    ABI_SFREE(ise)

  end subroutine set_niuv_sc


  !------------------------------
  ! fill oiju terms in supercell
  ! TODO: test
  !------------------------------
  subroutine set_oiju_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    class(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Rj_prim(3), Ru_prim(3), oiju_ind(5), iRj, iRu, i_prim, j_prim, u_prim, inz
    integer :: ngroup
    integer, allocatable :: i_sc(:), j_sc(:), u_sc(:), Rj_sc(:, :), Ru_sc(:,:)
    integer, allocatable :: i1list(:), ise(:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%oiju_sc%initialize(mshape=[self%nspin*3, self%nspin*3, self%natom*3])
    endif

    do inz=1, self%oiju%nnz
      oiju_ind=self%oiju%get_ind_inz(inz)
      iRj=oiju_ind(1)
      iRu=oiju_ind(2)
      i_prim=oiju_ind(3)
      j_prim=oiju_ind(4)
      u_prim=oiju_ind(5)
      Rj_prim=self%oRjlist%data(:,iRj)
      Ru_prim=self%oRulist%data(:,iRu)
      call scmaker%trans_i(nbasis=self%nspin*3, i=i_prim, i_sc=i_sc) 
      call scmaker%trans_j_and_Rj(nbasis=self%nspin*3, j=j_prim, Rj=Rj_prim, j_sc=j_sc, Rj_sc=Rj_sc)
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=u_prim, Rj=Ru_prim, j_sc=u_sc, Rj_sc=Ru_sc)
      val_sc(:)= self%oiju%val%data(inz)
      do icell=1, scmaker%ncells
        call scpot%add_oiju_term(i_sc(icell), j_sc(icell), u_sc(icell), val_sc(icell))
      end do
      ABI_SFREE(i_sc)
      ABI_SFREE(j_sc)
      ABI_SFREE(u_sc)
      ABI_SFREE(Rj_sc)
      ABI_SFREE(Ru_sc)
    end do
    
    call scpot%oiju_sc%group_by_1dim(ngroup, i1list, ise)
    ABI_SFREE(i1list)
    ABI_SFREE(ise)

  end subroutine set_oiju_sc

  subroutine set_tijuv_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Rj(3), Ru(3), Rv(3), tijuv_ind(7), iRj, iRu, iRv, ii, ij, iu, iv, inz
    integer, allocatable :: i_sc(:), j_sc(:), u_sc(:), v_sc(:), Rj_sc(:, :), Ru_sc(:,:), Rv_sc(:,:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%tijuv_sc%initialize(mshape=[self%nspin*3, self%nspin*3, self%natom*3, self%natom*3])
      call scpot%tuvij_sc%initialize(mshape=[self%natom*3, self%natom*3, self%nspin*3, self%nspin*3])
    endif

    do inz=1, self%tijuv%nnz
      tijuv_ind=self%tijuv%get_ind_inz(inz)
      iRj=tijuv_ind(1)
      iRu=tijuv_ind(2)
      iRv=tijuv_ind(3)
      ii=tijuv_ind(4)
      ij=tijuv_ind(5)
      iu=tijuv_ind(6)
      iv=tijuv_ind(7)
      Rj=self%tRjlist%data(:,iRj)
      Ru=self%tRulist%data(:,iRu)
      Rv=self%tRvlist%data(:,iRv)
      call scmaker%trans_i(nbasis=self%nspin*3, i=ii, i_sc=i_sc)
      call scmaker%trans_j_and_Rj(nbasis=self%nspin*3, j=ij, Rj=Rj, j_sc=j_sc, Rj_sc=Rj_sc)
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=iu, Rj=Ru, j_sc=u_sc, Rj_sc=Ru_sc)
      call scmaker%trans_j_and_Rj(nbasis=self%natom*3, j=iv, Rj=Rv, j_sc=v_sc, Rj_sc=Rv_sc)
      val_sc(:)= self%tijuv%val%data(inz)
      do icell=1, scmaker%ncells
        call scpot%add_tijuv_term(i_sc(icell), j_sc(icell), u_sc(icell), v_sc(icell), val_sc(icell))
      end do
      ABI_SFREE(i_sc)
      ABI_SFREE(j_sc)
      ABI_SFREE(u_sc)
      ABI_SFREE(v_sc)
      ABI_SFREE(Rj_sc)
      ABI_SFREE(Ru_sc)
      ABI_SFREE(Rv_sc)
    end do
    
    call scpot%tijuv_sc%group_by_pair()
    call scpot%tuvij_sc%group_by_pair()

  end subroutine set_tijuv_sc

  !!***
end module m_slc_primitive_potential
