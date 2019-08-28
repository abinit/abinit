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

  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_cell, only: mbcell_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_primitive_potential, only: primitive_potential_t
  use m_slc_potential, only: slc_potential_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_supercell_maker, only: supercell_maker_t
  use m_xmpi


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
     !procedure :: set_slc_primcell
     procedure :: set_liu
     procedure :: set_liu_1term
     procedure :: set_oiju
     procedure :: set_oiju_1term
     !procedure:: set_tijuv
     procedure :: load_from_files
     procedure :: read_netcdf
     procedure :: fill_supercell
     !procedure :: set_liu_sc
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
       ncdf_fname=fnames(3)
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
            &     'Reading spin-lattice coupling terms from'
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
    integer:: ndata
    
    integer, allocatable :: ilist(:), jlist(:), ulist(:), vlist(:)
    integer, allocatable :: Rjlist(:,:), Rulist(:,:), Rvlist(:,:)
    real(dp), allocatable :: vallist(:)

!#if defined HAVE_NETCDF

    ncerr = nf90_open(ncdf_fname,NF90_NOWRITE,ncid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A24)') 'Could not open netcdf file'
    else 
      write(std_out,'(A30)') ncdf_fname
    endif
    
    ! read primcell info
    ncerr=nctk_get_dim(ncid, "natom", self%natom)
    ncerr=nctk_get_dim(ncid, "nspin", self%nspin)

    ! read L_iu terms
    ncerr = nf90_inq_dimid(ncid, "Liu_ndata", dimid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A20)') 'Spin-lattice coupling does not contain L_iu term'
    else
      self%has_bilin=.True.
      ncerr = nctk_get_dim(ncid, "Liu_ndata", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(ulist, (ndata))
      ABI_ALLOCATE(Rulist, (3, ndata))
      ABI_ALLOCATE(vallist, (ndata))

      ncerr =nf90_inq_varid(ncid, "Liu_ilist", varid)
      call nc_handle_err(ncerr, "Liu_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call nc_handle_err(ncerr, "Liu_ilist")

      ncerr =nf90_inq_varid(ncid, "Liu_ulist", varid)
      call nc_handle_err(ncerr, "Liu_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call nc_handle_err(ncerr, "Liu_ulist")

      ncerr =nf90_inq_varid(ncid, "Liu_Rulist", varid)
      call nc_handle_err(ncerr, "Liu_Rulist")
      ncerr = nf90_get_var(ncid, varid, Rulist)
      call nc_handle_err(ncerr, "Liu_Rulist")

      ncerr =nf90_inq_varid(ncid, "Liu_vallist", varid)
      call nc_handle_err(ncerr, "Liu_vallist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call nc_handle_err(ncerr, "Liu_vallist")

      write(std_out,'(A7,I10,A11)') 'L_iu:', ndata, 'terms read'  

      !change units from eV to Ha
      vallist(:) = vallist(:)*eV_Ha
  
      !fill the sparse matrix for oiju parameters
      call self%set_liu(ndata, ilist, ulist, Rulist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(ulist)
      ABI_SFREE(Rulist)
      ABI_SFREE(vallist)

    endif


    ! read O_iju terms    
    ncerr = nf90_inq_dimid(ncid, "Oiju_ndata", dimid)
    if(ncerr /= NF90_NOERR) then
      write(std_out,'(A20)') 'Spin-lattice coupling does not contain O_iju term'
    else
      self%has_quadlin=.True.
      ncerr = nctk_get_dim(ncid, "Oiju_ndata", ndata)
      ABI_ALLOCATE(ilist, (ndata))
      ABI_ALLOCATE(jlist, (ndata))
      ABI_ALLOCATE(ulist, (ndata))
      ABI_ALLOCATE(Rjlist, (3, ndata))
      ABI_ALLOCATE(Rulist, (3, ndata))
      ABI_ALLOCATE(vallist, (ndata))

      ncerr =nf90_inq_varid(ncid, "Oiju_ilist", varid)
      call nc_handle_err(ncerr, "Oiju_ilist")
      ncerr = nf90_get_var(ncid, varid, ilist)
      call nc_handle_err(ncerr, "Oiju_ilist")

      ncerr =nf90_inq_varid(ncid, "Oiju_jlist", varid)
      call nc_handle_err(ncerr, "Oiju_jlist")
      ncerr = nf90_get_var(ncid, varid, jlist)
      call nc_handle_err(ncerr, "Oiju_jlist")

      ncerr =nf90_inq_varid(ncid, "Oiju_ulist", varid)
      call nc_handle_err(ncerr, "Oiju_ulist")
      ncerr = nf90_get_var(ncid, varid, ulist)
      call nc_handle_err(ncerr, "Oiju_ulist")

      ncerr =nf90_inq_varid(ncid, "Oiju_Rjlist", varid)
      call nc_handle_err(ncerr, "Oiju_Rjlist")
      ncerr = nf90_get_var(ncid, varid, Rjlist)
      call nc_handle_err(ncerr, "Oiju_Rjlist")

      ncerr =nf90_inq_varid(ncid, "Oiju_Rulist", varid)
      call nc_handle_err(ncerr, "Oiju_Rulist")
      ncerr = nf90_get_var(ncid, varid, Rulist)
      call nc_handle_err(ncerr, "Oiju_Rulist")

      ncerr =nf90_inq_varid(ncid, "Oiju_vallist", varid)
      call nc_handle_err(ncerr, "Oiju_vallist")
      ncerr = nf90_get_var(ncid, varid, vallist)
      call nc_handle_err(ncerr, "Oiju_vallist")

      write(std_out,'(A7,I10,A11)') 'O_iju:', ndata, 'terms read'  

      !change units from eV to Ha
      vallist(:) = vallist(:)*eV_Ha
  
      !fill the sparse matrix for oiju parameters
      call self%set_oiju(ndata, ilist, jlist, ulist, Rjlist, Rulist, vallist)

      ABI_SFREE(ilist)
      ABI_SFREE(jlist)
      ABI_SFREE(ulist)
      ABI_SFREE(Rjlist)
      ABI_SFREE(Rulist)
      ABI_SFREE(vallist)

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

  subroutine set_liu(self, nn, ilist, ulist, Rulist, vallist)

    class(slc_primitive_potential_t), intent(inout) :: self
    integer,                          intent(inout) :: nn
    integer,                          intent(in)    :: ilist(nn)
    integer,                          intent(in)    :: ulist(nn)
    integer,                          intent(in)    :: Rulist(3,nn)
    real(dp),                         intent(in)    :: vallist(nn)

    integer :: idx

    call self%liu%initialize(mshape=[-1, self%nspin*3, self%natom*3])
    
    if (xmpi_comm_rank(xmpi_world)==0) then
      do idx=1, nn
        call self%set_liu_1term(ilist(idx), ulist(idx), Rulist(:,idx), vallist(idx))
      end do
    endif
  end subroutine set_liu

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
    select type(scpot) ! use select type because properties only defined for slc_potential is used.
    type is (slc_potential_t) 
      call scpot%initialize(sc_nspin, sc_natom)
      call scpot%set_params(params)
      if (iam_master) then
        if(scpot%has_bilin) then
          if(self%has_bilin) then !did we actually find this in the netcdf file
            call self%liu%sum_duplicates()
            !call self%set_liu_sc(scpot, scmaker)
          else
            write(std_out,'(A47)') 'No parameters for bilinear coupling available'
            scpot%has_bilin = .False.
          endif
        endif
        if(scpot%has_quadlin) then
          if(self%has_quadlin) then
            call self%oiju%sum_duplicates()
            call self%set_oiju_sc(scpot, scmaker)
          else
            write(std_out,'(A55)') 'No parameters for quadratic-linear coupling available'
            scpot%has_quadlin = .False.
          endif
        endif
        if(scpot%has_linquad) then
          if(self%has_linquad) then
            call self%niuv%sum_duplicates()
            !call self%set_niuv_sc(scpot, scmaker)
          else
            write(std_out,'(A55)') 'No parameters for linear-quadratic coupling available'
            scpot%has_linquad = .False.
          endif
        endif
        if(scpot%has_biquad) then
          if(self%has_biquad) then
            call self%tijuv%sum_duplicates()
            !call self%set_tijuv_sc(scpot, scmaker)
          else
            write(std_out,'(A50)') 'No parameters for biquadratic coupling available'
            scpot%has_biquad = .False.
          endif
        endif
        !Write information which terms are used
        write(std_out,'(A55)') 'Using the following terms for the spin-lattice coupling'
        if(scpot%has_bilin)   write(std_out,'(A19)') 'Bilinear term: Liu'
        if(scpot%has_quadlin) write(std_out,'(A28)') 'Quadratic-linear term: Oiju'
        if(scpot%has_linquad) write(std_out,'(A28)') 'Linear-quadratic term: Niuv'
        if(scpot%has_biquad)  write(std_out,'(A24)') 'Biquadratic term: Tijuv'
      endif
    end select
  end subroutine fill_supercell

  subroutine set_oiju_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Rj(3), Ru(3), oiju_ind(5), iRj, iRu, ii, ij, iu, inz
    integer :: ngroup
    integer, allocatable :: i_sc(:), j_sc(:), u_sc(:), Rj_sc(:, :), Ru_sc(:,:)
    integer, allocatable :: i1list(:), ise(:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%oiju_sc%initialize(mshape=[-1, -1, self%nspin*3, self%nspin*3, self%natom*3])
    endif

    do inz=1, self%oiju%nnz
      oiju_ind=self%oiju%get_ind_inz(inz)
      iRj=oiju_ind(1)
      iRu=oiju_ind(2)
      ii=oiju_ind(3)
      ij=oiju_ind(4)
      iu=oiju_ind(5)
      Rj=self%oRjlist%data(:,iRj)
      Ru=self%oRulist%data(:,iRu)
      call scmaker%trans_i(nbasis=nspin*3, i=ii, i_sc=i_sc)
      call scmaker%trans_j_and_Rj(nbasis=nspin*3, j=ij, Rj=Rj, j_sc=j_sc, Rj_sc=Rj_sc)
      call scmaker%trans_j_and_Rj(nbasis=natom*3, j=iu, Rj=Ru, j_sc=u_sc, Rj_sc=Ru_sc)
      val_sc(:)= self%oiju%val%data(inz)
      do icell=1, scmaker%ncells
        call scpot%add_oiju_term(i_sc(icell), j_sc(icell), u_sc(icell), val_sc(icell))
      end do
      if(allocated(i_sc)) ABI_DEALLOCATE(i_sc)
      if(allocated(j_sc)) ABI_DEALLOCATE(j_sc)
      if(allocated(u_sc)) ABI_DEALLOCATE(u_sc)
      if(allocated(Rj_sc)) ABI_DEALLOCATE(Rj_sc)
      if(allocated(Ru_sc)) ABI_DEALLOCATE(Ru_sc)
    end do
    
    call scpot%oiju_sc%group_by_1dim(ngroup, i1list, ise)

  end subroutine set_oiju_sc

  subroutine set_tijuv_sc(self, scpot, scmaker)
    class(slc_primitive_potential_t), intent(inout) :: self
    type(slc_potential_t),            intent(inout) :: scpot
    type(supercell_maker_t),          intent(inout) :: scmaker
        
    integer :: icell, Rj(3), Ru(3), Rv(3), tijuv_ind(7), iRj, iRu, iRv, ii, ij, iu, iv, inz
    integer :: ngroup
    integer, allocatable :: i_sc(:), j_sc(:), u_sc(:), v_sc(:), Rj_sc(:, :), Ru_sc(:,:), Rv_sc(:,:)
    integer, allocatable :: i1list(:), ise(:)
    real(dp) :: val_sc(scmaker%ncells)

    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master

    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if(iam_master) then
      call scpot%tijuv_sc%initialize(mshape=[-1, -1, -1, self%nspin*3, self%nspin*3, self%natom*3, self%natom*3])
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
      call scmaker%trans_i(nbasis=nspin*3, i=ii, i_sc=i_sc)
      call scmaker%trans_j_and_Rj(nbasis=nspin*3, j=ij, Rj=Rj, j_sc=j_sc, Rj_sc=Rj_sc)
      call scmaker%trans_j_and_Rj(nbasis=natom*3, j=iu, Rj=Ru, j_sc=u_sc, Rj_sc=Ru_sc)
      call scmaker%trans_j_and_Rj(nbasis=natom*3, j=iv, Rj=Rv, j_sc=v_sc, Rj_sc=Rv_sc)
      val_sc(:)= self%tijuv%val%data(inz)
      do icell=1, scmaker%ncells
        call scpot%add_tijuv_term(i_sc(icell), j_sc(icell), u_sc(icell), v_sc(icell), val_sc(icell))
      end do
      if(allocated(i_sc)) ABI_DEALLOCATE(i_sc)
      if(allocated(j_sc)) ABI_DEALLOCATE(j_sc)
      if(allocated(u_sc)) ABI_DEALLOCATE(u_sc)
      if(allocated(v_sc)) ABI_DEALLOCATE(v_sc)
      if(allocated(Rj_sc)) ABI_DEALLOCATE(Rj_sc)
      if(allocated(Ru_sc)) ABI_DEALLOCATE(Ru_sc)
      if(allocated(Rv_sc)) ABI_DEALLOCATE(Rv_sc)
    end do
    
    call scpot%tijuv_sc%group_by_1dim(ngroup, i1list, ise)

  end subroutine set_tijuv_sc


end module m_slc_primitive_potential
