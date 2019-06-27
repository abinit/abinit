!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_primitive_potential
!! NAME
!! m_spin_primitive_potential
!!
!! FUNCTION
!! This module contains the atomic structures and the spin hamiltonian inside the primitive cell
!! which can be directly mapped to the xml file. It is not for the calculation, but for constructing
!! the hamiltonian in supercell. It is also used as input for the magnon band structure calculation.
!!
!! Datatypes:
!!  spin_primitive_potential_t
!!
!! Subroutines:
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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

module m_spin_primitive_potential
  use iso_c_binding
  use m_dynamic_array, only: int_array_type, real_array_type, int2d_array_type
  use m_mathfuncs
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_io_xml, only: xml_read_spin, xml_free_spin
  use m_multibinit_cell, only: mbcell_t
  use m_primitive_potential, only: primitive_potential_t
  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_supercell_maker, only: supercell_maker_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_spin_potential, only: spin_potential_t
  implicit none
  private
  !!*** 
  type, public, extends(primitive_potential_t) :: spin_primitive_potential_t
     integer :: natoms, nspin    ! on every mpi node
     type(ndcoo_mat_t) :: coeff  ! on only master node
     type(int2d_array_type) :: Rlist  ! on only master node
   contains
     procedure:: initialize
     procedure:: finalize
     procedure :: set_spin_primcell
     procedure :: set_bilinear_1term
     procedure:: set_bilinear
     procedure:: set_exchange
     procedure:: set_dmi
     procedure:: set_sia
     procedure :: add_input_sia
     procedure :: load_from_files
     procedure:: read_xml
     procedure:: fill_supercell
  end type spin_primitive_potential_t

contains

  subroutine initialize(self, primcell)
    class(spin_primitive_potential_t), intent(inout) :: self
    type(mbcell_t), target, intent(inout) :: primcell
    !integer, intent(in) :: nspin
    self%primcell=>primcell
    self%label="Spin_primitive_potential"
    self%has_spin=.True.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%has_lwf=.False.
  end subroutine initialize


  subroutine finalize(self)
    class(spin_primitive_potential_t), intent(inout) :: self
    call self%coeff%finalize()
    call self%Rlist%finalize()
    nullify(self%primcell)
    self%nspin=0
    self%natoms=0
    self%label="Destroyed Spin_primitive_potential"
    call self%primitive_potential_t%finalize()
  end subroutine finalize


  subroutine set_spin_primcell(self, natoms, unitcell, positions, &
       nspin, index_spin, spinat, gyroratios, damping_factors )
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(inout):: natoms, nspin, index_spin(:)
    real(dp), intent(inout):: unitcell(3, 3),  positions(3,natoms), &
         spinat(3,natoms), gyroratios(nspin), damping_factors(nspin)
    integer :: iatom, ispin
    real(dp) :: ms(nspin), spin_positions(3, nspin)
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 
    
    ABI_UNUSED_A(unitcell)

    self%nspin=nspin
    call xmpi_bcast(self%nspin, master, comm, ierr)
    if (iam_master) then
       call self%coeff%initialize(mshape=[-1, self%nspin*3, self%nspin*3])
       do iatom=1, natoms
          ispin=index_spin(iatom)
          if(ispin>0) then
             spin_positions(:,ispin)= positions(:, iatom)
             ms(ispin) = sqrt(sum(spinat(:,iatom)**2, dim=1))* mu_B
          end if
       end do
    endif
    call self%primcell%set_spin(nspin, ms, spin_positions, gyroratios, damping_factors)
  end subroutine set_spin_primcell


  subroutine load_from_files(self, params, fnames)
    class(spin_primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    character(len=fnlen) :: xml_fname
    character(len=500) :: message
    integer :: ii
    logical:: use_sia, use_exchange, use_dmi, use_bi

    if (xmpi_comm_rank(xmpi_world)==0) then
       xml_fname=fnames(3)
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
            &     'reading spin terms.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
    endif
    use_exchange=.True.
    use_sia=.True.
    use_dmi=.True.
    use_bi=.True.
    ! Do not use sia term in xml if spin_sia_add is set to 1.
    if(params%spin_sia_add == 1) use_sia=.False.
    call self%read_xml( trim(xml_fname)//char(0), &
         & use_exchange=use_exchange,  use_sia=use_sia, use_dmi=use_dmi, use_bi=use_bi)
    if (params%spin_sia_add /= 0 ) then
       call self%add_input_sia(params%spin_sia_k1amp, &
            & params%spin_sia_k1dir)
    end if
  end subroutine load_from_files


  subroutine set_bilinear_1term(self, i, j, R, val)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: i, j, R(3)
    real(dp), intent(in) :: val(3,3)
    real(dp) :: v
    integer :: indR, iv, jv



    if (xmpi_comm_rank(xmpi_world)==0) then
       call self%Rlist%push_unique(R, position=indR)
       do jv=1,3
          do iv=1,3
             v=val(iv,jv)
             call self%coeff%add_entry(ind=[indR, (i-1)*3+iv, (j-1)*3+jv ], val=v)
          end do
       end do
    endif
  end subroutine set_bilinear_1term

  subroutine  set_bilinear(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(n), jlist(n), Rlist(3,n)
    real(dp), intent(in) :: vallist(3, 3,n)
    integer :: idx

    if (xmpi_comm_rank(xmpi_world)==0) then
       do idx = 1, n
          call self%set_bilinear_1term(ilist(idx), jlist(idx), Rlist(:,idx), vallist(:,:, idx))
       end do
    endif
  end subroutine set_bilinear

  subroutine set_exchange(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: idx
    real(dp) :: bivallist(3,3, n)

    if (xmpi_comm_rank(xmpi_world)==0) then
       bivallist(:,:,:)=0.0d0
       do idx = 1, n, 1
          bivallist(1,1,idx)=vallist(1, idx)
          bivallist(2,2,idx)=vallist(2, idx)
          bivallist(3,3,idx)=vallist(3, idx)
       end do
       call self%set_bilinear(n,ilist,jlist,Rlist,bivallist)
    endif
  end subroutine set_exchange


  subroutine set_dmi(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: idx
    real(dp) :: bivallist(3,3, n), D(3)

    if (xmpi_comm_rank(xmpi_world)==0) then
       bivallist(:,:,:)=0.0d0
       do idx=1,n, 1
          D(:)=vallist(:, idx)
          ! 0 Dz -Dy
          ! -Dz 0 Dx
          ! Dy -Dx 0
          bivallist(:,:, idx)=reshape ( (/0.0d0, -D(3), D(2),  &
               D(3), 0.0d0, -D(1),  &
               -D(2), D(1), 0.0d0 /),(/3,3/) )
       end do
       call self%set_bilinear(n,ilist,jlist,Rlist,bivallist)
    endif
  end subroutine set_dmi

  subroutine set_sia(self, n, ilist, k1list, k1dirlist)

    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:)
    real(dp), intent(in) :: k1list(:), k1dirlist(:, :)
    integer :: idx, Rlist(3, n)
    real(dp) :: bivallist(3,3, n)
    if (xmpi_comm_rank(xmpi_world)==0) then
       bivallist(:,:,:)=0.0d0
       Rlist(:, :)=0.0d0
       do idx=1,n, 1
          bivallist(:,:, idx)= (- k1list(idx))*  &
               outer_product(k1dirlist(:,idx), k1dirlist(:, idx))
       end do
       call self%set_bilinear(n,ilist,ilist,Rlist,bivallist)
    endif
  end subroutine set_sia

  ! add user inputed SIA.
  subroutine add_input_sia(self,  input_sia_k1amp, input_sia_k1dir)
    class(spin_primitive_potential_t), intent(inout) :: self
    real(dp), intent(in):: input_sia_k1amp, input_sia_k1dir(3)
    integer :: in_sia_ind(self%nspin)
    real(dp)::  in_sia_k1amp(self%nspin), in_sia_k1dir(3, self%nspin)

    integer :: i

    if (xmpi_comm_rank(xmpi_world)==0) then
       write(std_out,'(A28)') "Adding SIA terms from input"
       do i =1, self%nspin
          in_sia_ind(i)=i
          in_sia_k1amp(i)=input_sia_k1amp
          in_sia_k1dir(:,i)=input_sia_k1dir
       end do
       call self%set_sia(self%nspin, in_sia_ind, in_sia_k1amp, in_sia_k1dir )
    endif
  end subroutine add_input_sia


  subroutine read_xml(self, xml_fname, use_exchange, use_dmi, use_sia, use_bi)
    class(spin_primitive_potential_t), intent(inout) :: self
    character(kind=C_CHAR) :: xml_fname(*)
    integer :: natoms, nspin, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
    logical, optional, intent(in) :: use_exchange, use_dmi, use_sia, use_bi
    logical :: uexc, udmi, usia, ubi
    real(dp) :: ref_energy

    type(c_ptr) ::  p_unitcell,         &
         p_masses,  p_index_spin, p_gyroratios, p_damping_factors, p_positions, p_spinat, &
         p_exc_ilist, p_exc_jlist, p_exc_Rlist, p_exc_vallist, &
         p_dmi_ilist, p_dmi_jlist, p_dmi_Rlist, p_dmi_vallist, &
         p_uni_ilist, p_uni_amplitude_list, p_uni_direction_list, &
         p_bi_ilist, p_bi_jlist, p_bi_Rlist, p_bi_vallist

    integer(c_int),pointer :: index_spin(:)=>null() ,&
         exc_ilist(:)=>null(), exc_jlist(:)=>null(),  exc_Rlist(:)=>null(), &
         dmi_ilist(:)=>null(), dmi_jlist(:)=>null(),  dmi_Rlist(:)=>null(), &
         bi_ilist(:)=>null(), bi_jlist(:)=>null(), bi_Rlist(:)=>null(), &
         uni_ilist(:)=>null()

    real(c_double), pointer:: unitcell(:)=>null(), masses(:)=>null(),  &
         gyroratios(:)=>null(), damping_factors(:)=>null(), &
         positions(:)=>null(), spinat(:)=>null(), &
         exc_vallist(:)=>null(), dmi_vallist(:)=>null(), &
         uni_amplitude_list(:)=>null(), uni_direction_list(:)=>null(), &
         bi_vallist(:)=>null()

    real(dp) :: uc(3,3)

    integer :: master, my_rank, comm, nproc
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    if (iam_master) then
       write(std_out,'(A58)') "Reading parameters from xml file and setting up spin model"
       write(std_out,'(A80)') " "
       call xml_read_spin(xml_fname, ref_energy, p_unitcell,                 &
            natoms, p_masses, nspin, p_index_spin, p_gyroratios, p_damping_factors, p_positions, p_spinat, &
            exc_nnz, p_exc_ilist, p_exc_jlist, p_exc_Rlist, p_exc_vallist, &
            dmi_nnz, p_dmi_ilist, p_dmi_jlist, p_dmi_Rlist, p_dmi_vallist, &
            uni_nnz, p_uni_ilist, p_uni_amplitude_list, p_uni_direction_list, &
            bi_nnz, p_bi_ilist, p_bi_jlist, p_bi_Rlist, p_bi_vallist)
       call c_f_pointer(p_unitcell, unitcell, [9])
       call c_f_pointer(p_masses, masses, [natoms])
       call c_f_pointer(p_index_spin, index_spin, [natoms])
       call c_f_pointer(p_gyroratios, gyroratios, [nspin])
       call c_f_pointer(p_damping_factors, damping_factors, [nspin])
       call c_f_pointer(p_positions, positions, [natoms*3])
       call c_f_pointer(p_spinat, spinat, [natoms*3])
       call c_f_pointer(p_exc_ilist, exc_ilist, [exc_nnz])
       call c_f_pointer(p_exc_jlist, exc_jlist, [exc_nnz])
       call c_f_pointer(p_exc_Rlist, exc_Rlist, [exc_nnz*3])
       call c_f_pointer(p_exc_vallist, exc_vallist, [exc_nnz*3])
       call c_f_pointer(p_dmi_ilist, dmi_ilist, [dmi_nnz])
       call c_f_pointer(p_dmi_jlist, dmi_jlist, [dmi_nnz])
       call c_f_pointer(p_dmi_Rlist, dmi_Rlist, [dmi_nnz*3])
       call c_f_pointer(p_dmi_vallist, dmi_vallist, [dmi_nnz*3])
       call c_f_pointer(p_uni_ilist, uni_ilist, [uni_nnz])
       call c_f_pointer(p_uni_amplitude_list, uni_amplitude_list, [uni_nnz])
       call c_f_pointer(p_uni_direction_list, uni_direction_list, [uni_nnz*3])
       call c_f_pointer(p_bi_ilist, bi_ilist, [bi_nnz])
       call c_f_pointer(p_bi_jlist, bi_jlist, [bi_nnz])
       call c_f_pointer(p_bi_Rlist, bi_Rlist, [bi_nnz*3])
       call c_f_pointer(p_bi_vallist, bi_vallist, [bi_nnz*9])
       
       ! change of units to a.u.
       
       ! unitcell already Bohr
       
       !gyroratios already in a.u. (unit=1)
       
       ! masses already in a.u.
       
       ! J, DMI, k1, bi are in eV
       exc_vallist(:) =exc_vallist(:) * eV_Ha
       dmi_vallist(:) = dmi_vallist(:) * eV_Ha
       uni_amplitude_list(:) = uni_amplitude_list(:) * eV_Ha
       bi_vallist(:) = bi_vallist(:) * eV_Ha
       
       write(std_out,'(A80)') " "
       write(std_out,'(A21)') "Setting up spin model"
       write(std_out,'(A15)') "Setting system"
       uc(:,:)=transpose(reshape(unitcell, [3,3]))
       !call set_atoms(self,)
    endif

    ! (MPI) Only this runs on non-master node
    call self%set_spin_primcell(natoms, uc, positions, &
         nspin, index_spin, spinat, gyroratios, damping_factors )

    if (iam_master) then
       if(.not. present(use_exchange))  then
          uexc=.True.
       else
          uexc=use_exchange
       end if

       if(uexc .and. exc_nnz>0) then
          write(std_out,'(A23)') "Setting exchange terms"
          call self%set_exchange(exc_nnz,exc_ilist,exc_jlist,&
               reshape(exc_Rlist, (/3, exc_nnz /)), &
               reshape(exc_vallist, (/3, exc_nnz/)))
       else
        if (.not. uexc)  write(std_out, '(A38)') " Exchange term from xml file not used."
       endif

       if(.not. present(use_dmi))  then
          udmi=.True.
       else
          udmi=use_dmi
       end if
       if (udmi .and. dmi_nnz>0) then
          write(std_out,'(A19)') "Setting DMI terms."
          call self%set_dmi( n=dmi_nnz, ilist=dmi_ilist, jlist=dmi_jlist, &
               Rlist=reshape(dmi_Rlist, (/3, dmi_nnz /)), &
               vallist = reshape(dmi_vallist, (/3, dmi_nnz/)))
       else
          if (.not. udmi) write(std_out, '(A35)') " DMI term from xml file not used."
       end if

       if(.not. present(use_sia)) then
          usia=.True.
       else
          usia=use_sia
       end if
       if (usia .and. uni_nnz>0) then
          write(std_out,'(A18)') "Setting SIA terms"
          call self%set_sia(uni_nnz, uni_ilist, uni_amplitude_list, &
               reshape(uni_direction_list, [3, uni_nnz]) )
       else
         if(.not. usia) write(std_out,'(A34)') " SIA term in xml file not used."
       end if
       if(.not. present(use_bi)) then
          ubi=.True.
       else
          ubi=use_bi
       endif
       if (ubi .and. bi_nnz>0) then
          write(std_out,'(A23)') "Setting bilinear terms."
          call self%set_bilinear(bi_nnz, bi_ilist, bi_jlist,  &
               Rlist=reshape(bi_Rlist, (/3, bi_nnz /)), &
               vallist = reshape(bi_vallist, (/3,3, bi_nnz/)))
       else
          if(.not. ubi) write(std_out, '(A38)') " Bilinear term in xml file not used."
       endif
    endif

    if (iam_master) then
       call xml_free_spin(xml_fname, ref_energy, p_unitcell,                 &
            natoms, p_masses, nspin, p_index_spin, p_gyroratios, p_damping_factors, p_positions, p_spinat, &
            exc_nnz, p_exc_ilist, p_exc_jlist, p_exc_Rlist, p_exc_vallist, &
            dmi_nnz, p_dmi_ilist, p_dmi_jlist, p_dmi_Rlist, p_dmi_vallist, &
            uni_nnz, p_uni_ilist, p_uni_amplitude_list, p_uni_direction_list, &
            bi_nnz, p_bi_ilist, p_bi_jlist, p_bi_Rlist, p_bi_vallist)
    endif

  end subroutine read_xml


  subroutine fill_supercell(self, scmaker, scpot)
    class(spin_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t), intent(inout):: scmaker
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    integer :: nspin, sc_nspin, i, R(3), ind_Rij(3), iR, ii, ij, inz
    integer, allocatable :: i_sc(:), j_sc(:), Rj_sc(:, :)
    real(dp) :: val_sc(scmaker%ncells)
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    nspin=self%nspin
    sc_nspin= nspin * scmaker%ncells
    call xmpi_bcast(sc_nspin, master, comm, ierr)
    !ABI_MALLOC_SCALAR(spin_potential_t::scpot)
    ABI_DATATYPE_ALLOCATE_SCALAR(spin_potential_t, scpot)
    select type(scpot)
    type is (spin_potential_t)
       call scpot%initialize(sc_nspin)
       if (iam_master) then
          call self%coeff%sum_duplicates()
          do inz=1, self%coeff%nnz
             ind_Rij=self%coeff%get_ind_inz(inz)
             iR=ind_Rij(1)
             ii=ind_Rij(2)
             ij=ind_Rij(3)
             R=self%Rlist%data(:,iR)
             call scmaker%trans_i(nbasis=nspin*3, i=ii, i_sc=i_sc)
             call scmaker%trans_j_and_Rj(nbasis=nspin*3, j=ij, Rj=R, j_sc=j_sc, Rj_sc=Rj_sc)
             val_sc(:)= self%coeff%val%data(inz)
             do i=1, scmaker%ncells
                call scpot%add_bilinear_term(i_sc(i), j_sc(i), val_sc(i))
             end do
             if(allocated(i_sc)) ABI_DEALLOCATE(i_sc)
             if(allocated(j_sc)) ABI_DEALLOCATE(j_sc)
             if(allocated(Rj_sc)) ABI_DEALLOCATE(Rj_sc)
          end do
       endif
    end select
  end subroutine fill_supercell

end module m_spin_primitive_potential
