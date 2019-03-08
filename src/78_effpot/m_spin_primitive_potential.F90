!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sspin_primitive_potential
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
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_io_xml, only: xml_read_spin, xml_free_spin
  use m_unitcell, only: unitcell_t
  use m_primitive_potential, only: primitive_potential_t
  use m_abstract_potential, only: abstract_potential_t
  use m_dynamic_array, only: int2d_array_type
  use m_supercell_maker, only: supercell_maker_t
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_multibinit_global
  use m_spin_potential, only: spin_potential_t
  implicit none
  private
  !!*** 
  type, public, extends(primitive_potential_t) :: spin_primitive_potential_t
     integer :: natoms, nspin
     ! TOTAL
     !integer :: total_nnz=0
     !type(int_array_type):: total_ilist, total_jlist, total_Rlist(3)
     !type(real_array_type) :: total_val_list(3,3)
     type(ndcoo_mat_t) :: coeff
     type(int2d_array_type) :: Rlist
     !integer, allocatable :: index_spin(:)
   contains
     procedure:: initialize
     procedure:: finalize
     procedure:: set_atoms
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
     procedure :: print_terms
  end type spin_primitive_potential_t

contains

  subroutine initialize(self)
    class(spin_primitive_potential_t), intent(inout) :: self
    !integer, intent(in) :: nspin
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

  subroutine set_atoms(self, primcell)
    class(spin_primitive_potential_t), intent(inout) :: self
    class(unitcell_t), target, intent(inout) :: primcell
    self%primcell=>primcell
  end subroutine set_atoms


  subroutine set_spin_primcell(self, natoms, unitcell, positions, &
       nspin, index_spin, spinat, gyroratios, damping_factors )
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in):: natoms, nspin, index_spin(:)
    real(dp), intent(in):: unitcell(3, 3),  positions(3,natoms), &
         spinat(3,natoms), gyroratios(nspin), damping_factors(nspin)
    integer :: iatom, ispin
    real(dp) :: ms(nspin), spin_positions(3, nspin)
    self%nspin=nspin
    call self%coeff%initialize(shape=[-1, nspin*3, nspin*3])
    do iatom=1, natoms
       ispin=index_spin(iatom)
       if(iatom>0) then
          spin_positions(:,ispin)= positions(:, iatom)
          ms(ispin) = sqrt(sum(spinat(:,iatom)**2, dim=1))* mu_B
       end if
       call self%primcell%set_spin(nspin, ms, spin_positions, gyroratios, damping_factors)
    end do
  end subroutine set_spin_primcell


  subroutine load_from_files(self, params, fnames)
    class(spin_primitive_potential_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(:)
    character(len=fnlen) :: xml_fname
    character(len=500) :: message
    integer :: ii
    logical:: use_sia, use_exchange, use_dmi, use_bi

    xml_fname=fnames(3)
    if (iam_master) then
       write(message,'(a,(80a),3a)') ch10,('=',ii=1,80),ch10,ch10,&
            &     'reading spin terms.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
    end if

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
    call self%Rlist%push_unique(R, position=indR)
    do jv=1,3
       do iv=1,3
          v=val(iv,jv)
          call self%coeff%add_entry(ind=[indR, (i-1)*3+iv, (j-1)*3+jv ], val=v)
       end do
    end do
  end subroutine set_bilinear_1term

  subroutine  set_bilinear(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(n), jlist(n), Rlist(3,n)
    real(dp), intent(in) :: vallist(3, 3,n)
    integer :: idx
    do idx = 1, n
       call self%set_bilinear_1term(ilist(idx), jlist(idx), Rlist(:,idx), vallist(:,:, idx))
    end do
  end subroutine set_bilinear

  subroutine set_exchange(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: idx
    real(dp) :: bivallist(3,3, n)
    bivallist(:,:,:)=0.0d0
    do idx = 1, n, 1
       bivallist(1,1,idx)=vallist(1, idx)
       bivallist(2,2,idx)=vallist(2, idx)
       bivallist(3,3,idx)=vallist(3, idx)
    end do
    call self%set_bilinear(n,ilist,jlist,Rlist,bivallist)
  end subroutine set_exchange


  subroutine set_dmi(self, n, ilist, jlist, Rlist, vallist)
    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: idx
    real(dp) :: bivallist(3,3, n), D(3)
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
  end subroutine set_dmi

  subroutine set_sia(self, n, ilist, k1list, k1dirlist)

    class(spin_primitive_potential_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:)
    real(dp), intent(in) :: k1list(:), k1dirlist(:, :)
    integer :: idx, Rlist(3, n)
    real(dp) :: bivallist(3,3, n)
    bivallist(:,:,:)=0.0d0
    Rlist(:, :)=0.0d0
    do idx=1,n, 1
       bivallist(:,:, idx)= (- k1list(idx))*  &
            outer_product(k1dirlist(:,idx), k1dirlist(:, idx))
    end do
    call self%set_bilinear(n,ilist,ilist,Rlist,bivallist)
  end subroutine set_sia

  ! add user inputed SIA.
  subroutine add_input_sia(self,  input_sia_k1amp, input_sia_k1dir)
    class(spin_primitive_potential_t), intent(inout) :: self
    real(dp), intent(in):: input_sia_k1amp, input_sia_k1dir(3)
    integer :: in_sia_ind(self%nspin)
    real(dp)::  in_sia_k1amp(self%nspin), in_sia_k1dir(3, self%nspin)

    integer :: i

    write(*,'(A28)') "Adding SIA terms from input"
    do i =1, self%nspin
       in_sia_ind=i
       in_sia_k1amp(i)=input_sia_k1amp
       in_sia_k1dir(:,i)=input_sia_k1dir
    end do
    call self%set_sia(self%nspin, in_sia_ind, in_sia_k1amp, in_sia_k1dir )
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

    write(*,'(A58)') "Reading parameters from xml file and setting up spin model"
    write(*,'(A80)') " "

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

    write(*,'(A80)') " "
    write(*,'(A21)') "Setting up spin model"
    write(*,'(A15)') "Setting system"
    uc(:,:)=transpose(reshape(unitcell, [3,3]))
    !call set_atoms(self,)
    
    call self%set_spin_primcell(natoms, unitcell, positions, &
         nspin, index_spin, spinat, gyroratios, damping_factors )


    if(.not. present(use_exchange))  then
       uexc=.True.
    else
       uexc=use_exchange
    end if

    if(uexc) then
       write(*,'(A23)') "Setting exchange terms"
       call self%set_exchange(exc_nnz,exc_ilist,exc_jlist,&
            reshape(exc_Rlist, (/3, exc_nnz /)), &
            reshape(exc_vallist, (/3, exc_nnz/)))
    else
       write(*, '(A38)') " Exchange term from xml file not used."
    endif

    if(.not. present(use_dmi))  then
       udmi=.True.
    else
       udmi=use_dmi
    end if
    if (udmi) then
       write(*,'(A18)') "Setting DMI terms"
       call self%set_dmi( n=dmi_nnz, ilist=dmi_ilist, jlist=dmi_jlist, &
            Rlist=reshape(dmi_Rlist, (/3, dmi_nnz /)), &
            vallist = reshape(dmi_vallist, (/3, dmi_nnz/)))
       write(*,'(A27)') " Setting uniaxial SIA terms"
    else
       print *, " DMI term from xml file not used"
    end if

    if(.not. present(use_sia)) then
       usia=.True.
    else
       usia=use_sia
    end if
    if (usia) then
       write(*,'(A18)') "Setting SIA terms"
       call self%set_sia(uni_nnz, uni_ilist, uni_amplitude_list, &
            reshape(uni_direction_list, [3, uni_nnz]) )
    else
       print *, " SIA term in xml file not used"
    end if

    if(.not. present(use_bi)) then
       ubi=.True.
    else
       ubi=use_bi
    endif

    if (ubi) then
       write(*,'(A23)') "Setting bilinear terms"
       call self%set_bilinear(bi_nnz, bi_ilist, bi_jlist,  &
            Rlist=reshape(bi_Rlist, (/3, bi_nnz /)), &
            vallist = reshape(bi_vallist, (/3,3, bi_nnz/)))
    else
       print *, " Bilinear term in xml file not used"
    endif

    call xml_free_spin(xml_fname, ref_energy, p_unitcell,                 &
         natoms, p_masses, nspin, p_index_spin, p_gyroratios, p_damping_factors, p_positions, p_spinat, &
         exc_nnz, p_exc_ilist, p_exc_jlist, p_exc_Rlist, p_exc_vallist, &
         dmi_nnz, p_dmi_ilist, p_dmi_jlist, p_dmi_Rlist, p_dmi_vallist, &
         uni_nnz, p_uni_ilist, p_uni_amplitude_list, p_uni_direction_list, &
         bi_nnz, p_bi_ilist, p_bi_jlist, p_bi_Rlist, p_bi_vallist)

  end subroutine read_xml


  subroutine fill_supercell(self, scmaker, scpot)
    class(spin_primitive_potential_t) , intent(inout) :: self
    type(supercell_maker_t), intent(inout):: scmaker
    class(abstract_potential_t), pointer, intent(inout) :: scpot
    type(spin_potential_t), pointer:: tpot
    integer :: nspin, sc_nspin, i, R(3), ind_Rij(3), iR, ii, ij, inz
    integer, allocatable :: i_sc(:), j_sc(:), R_sc(:, :)
    real(dp) :: val_sc(scmaker%ncells)

    nspin=self%nspin
    sc_nspin= nspin * scmaker%ncells
    allocate(spin_potential_t::tpot)
    call tpot%initialize(sc_nspin)
    call self%coeff%sum_duplicates()
    do inz=1, self%coeff%nnz
       ind_Rij=self%coeff%get_ind_inz(inz)
       iR=ind_Rij(1)
       ii=ind_Rij(2)
       ij=ind_Rij(3)
       R=self%Rlist%data(:,iR)
       call scmaker%trans_i(nbasis=nspin*3, i=ii, i_sc=i_sc)
       call scmaker%trans_j_and_Rj(nbasis=nspin*3, j=ij, Rj=R, j_sc=j_sc, Rj_sc=R_sc)
       val_sc(:)= self%coeff%val%data(inz)
       do i=1, scmaker%ncells
          call tpot%add_bilinear_term(i_sc(i), j_sc(i), val_sc(i))
       end do
    end do
    scpot=>tpot
  end subroutine fill_supercell



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine make_supercell(self, sc_matrix, sc_ham)                                                       !
  !   class(spin_primitive_potential_t) , intent(inout) :: self                                              !
  !   type(spin_potential_t) , intent(inout) :: sc_ham                                                           !
  !   integer :: sc_matrix(3,3), iatoms(self%nspin)                                                         !
  !                                                                                                          !
  !   integer ::  typat_primcell(self%natoms), sc_nspin,  i, counter, icell                                 !
  !   real(dp) :: znucl(self%natoms), tmp(3,3)                                                               !
  !   integer, allocatable ::sc_index_spin(:), sc_znucl(:), sc_iatoms(:)                                     !
  !   integer, allocatable :: sc_ispin_prim(:), sc_rvec(:, :)                                                !
  !   real(dp), allocatable ::sc_spinat(:,:), sc_gyroratios(:), sc_damping_factors(:), sc_spinpos(:,:)       !
  !   integer :: ii, jj, icol, irow, rr(3), R_sc(3), iatom                                                   !
  !                                                                                                          !
  !   if(iam_master) then                                                                                    !
  !      typat_primcell(:)=1                                                                                 !
  !      ! TODO hexu:should use supercell generate with lattice model                                        !
  !      !init_supercell(natom_primcell, sc_matrix, rprimd_primcell,                                         !
  !      !             typat_primcell, xcart_primcell, znucl, scell)                                         !
  !                                                                                                          !
  !      znucl(:)=0                                                                                          !
  !                                                                                                          !
  !      counter=0                                                                                           !
  !      do i=1, self%natoms                                                                                 !
  !         if(self%index_spin(i)>0) then                                                                    !
  !            counter=counter+1                                                                             !
  !            iatoms(counter)=i                                                                             !
  !         endif                                                                                            !
  !      enddo                                                                                               !
  !                                                                                                          !
  !      call init_supercell(self%natoms, sc_matrix, self%unitcell, typat_primcell, &                        !
  !           self%positions, znucl, scell)                                                                  !
  !      ! cell, positions                                                                                   !
  !      !scell%rprimd                                                                                       !
  !      !scell%xcart                                                                                        !
  !      !nspin                                                                                             !
  !      sc_nspin=scell%ncells*self%nspin                                                                  !
  !      ABI_ALLOCATE(sc_index_spin, (self%natoms*scell%ncells))                                             !
  !      ABI_ALLOCATE(sc_ispin_prim, (sc_nspin) )                                                           !
  !      ABI_ALLOCATE(sc_rvec, (3, sc_nspin) )                                                              !
  !      ABI_ALLOCATE(sc_spinat, (3, sc_nspin) )                                                            !
  !      ABI_ALLOCATE(sc_spinpos, (3, sc_nspin) )                                                           !
  !      ABI_ALLOCATE(sc_gyroratios, (sc_nspin))                                                            !
  !      ABI_ALLOCATE(sc_damping_factors, (sc_nspin))                                                       !
  !      !ABI_ALLOCATE(sc_znucl, (sc_nspin))                                                                !
  !      ABI_ALLOCATE(sc_iatoms, (sc_nspin))                                                                !
  !      ! sc_index_spin                                                                                     !
  !      counter=0                                                                                           !
  !      do i = 1, scell%natom, 1                                                                            !
  !         iatom=scell%atom_indexing(i)                                                                     !
  !         !if(scell%atom_indexing(i)>0) then                                                               !
  !         if(self%index_spin(iatom)>0) then                                                                !
  !            counter=counter+1                                                                             !
  !            !map from spin to atom and atom to spin.                                                      !
  !            sc_index_spin(i)=counter                                                                      !
  !            ! variables which every atom have in primitive cell                                           !
  !            sc_iatoms(counter)=iatom                                                                      !
  !            sc_spinpos(:, counter)=scell%xcart(:,i)                                                       !
  !            sc_spinat(:, counter)=self%spinat(:,iatom)                                                    !
  !                                                                                                          !
  !            sc_ispin_prim(counter) = self%index_spin(iatom)                                               !
  !                                                                                                          !
  !            ! variables only atoms with spin have.                                                        !
  !            sc_gyroratios( counter)=self%gyroratios(self%index_spin(iatom))                               !
  !            sc_damping_factors( counter)=self%damping_factors(self%index_spin(iatom))                     !
  !                                                                                                          !
  !            ! Rvec                                                                                        !
  !            sc_rvec(:,counter)=scell%uc_indexing(:,i)                                                     !
  !         else                                                                                             !
  !            sc_index_spin(i)=-1                                                                           !
  !         endif                                                                                            !
  !      end do                                                                                              !
  !                                                                                                          !
  !   endif                                                                                                  !
  !                                                                                                          !
  !                                                                                                          !
  !   call spin_potential_t_initialize(sc_ham, cell=scell%rprimd, pos=sc_spinpos, &                              !
  !        & spinat=sc_spinat, iatoms=sc_iatoms, ispin_prim=sc_ispin_prim, &                                 !
  !        & rvec=sc_rvec, gyro_ratio=sc_gyroratios, damping=sc_damping_factors)                             !
  !                                                                                                          !
  !   call  xmpi_bcast(self%total_nnz, master, comm, ierr)                                                   !
  !   call  xmpi_bcast(scell%ncells, master, comm, ierr)                                                     !
  !   do i =1, self%total_nnz, 1                                                                             !
  !      do icell=1, scell%ncells, 1                                                                         !
  !         if(iam_master) then                                                                              !
  !            ! Note i0 and j0 are in spin index, while find_supercell_ijR work in atom index               !
  !            call find_supercell_ijR(scell=scell, i0=iatoms(self%total_ilist%data(i)), &                   !
  !                 j0=iatoms(self%total_jlist%data(i)), &                                                   !
  !                 R0=[self%total_Rlist(1)%data(i),  &                                                      !
  !                 self%total_Rlist(2)%data(i),  &                                                          !
  !                 self%total_Rlist(3)%data(i)], &                                                          !
  !                 R=scell%rvecs(:,icell), &                                                                !
  !                 i1=ii,j1=jj,R1=rr,R_sc=R_sc)                                                             !
  !            do irow = 1, 3                                                                                !
  !               do icol=1, 3                                                                               !
  !                  tmp(icol, irow)=self%total_val_list(icol,irow)%data(i)                                  !
  !               end do                                                                                     !
  !            end do                                                                                        !
  !            call spin_potential_t_set_bilinear_term_single(sc_ham, sc_index_spin(ii), sc_index_spin(jj), tmp) !
  !         endif                                                                                            !
  !      enddo                                                                                               !
  !   enddo                                                                                                  !
  !   if (allocated(sc_ispin_prim)) then                                                                     !
  !      ABI_DEALLOCATE(sc_ispin_prim)                                                                       !
  !   endif                                                                                                  !
  !   if (allocated(sc_iatoms)) then                                                                         !
  !      ABI_DEALLOCATE(sc_iatoms)                                                                           !
  !   endif                                                                                                  !
  !   if (allocated(sc_spinpos)) then                                                                        !
  !      ABI_DEALLOCATE(sc_spinpos)                                                                          !
  !   endif                                                                                                  !
  !   if (allocated(sc_rvec)) then                                                                           !
  !      ABI_DEALLOCATE(sc_rvec)                                                                             !
  !   endif                                                                                                  !
  !   if (allocated(sc_index_spin)) then                                                                     !
  !      ABI_DEALLOCATE(sc_index_spin)                                                                       !
  !   endif                                                                                                  !
  !                                                                                                          !
  !   if (allocated(sc_gyroratios)) then                                                                     !
  !      ABI_DEALLOCATE(sc_gyroratios)                                                                       !
  !   endif                                                                                                  !
  !   if (allocated(sc_damping_factors)) then                                                                !
  !      ABI_DEALLOCATE(sc_damping_factors)                                                                  !
  !  end if                                                                                                  !
  !   if (allocated(sc_spinat)) then                                                                         !
  !      ABI_DEALLOCATE(sc_spinat)                                                                           !
  !   endif                                                                                                  !
  ! end subroutine make_supercell                                                                            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_terms(self)
    class(spin_primitive_potential_t) :: self
    integer :: i, ii, jj,  R(3)
    character(len=80) :: msg
    real(dp) :: tmp(3, 3)
  !   msg=repeat("=", 34)
  !   write(msg, '(A34, 1X, A10, 1x, A34)') msg, 'spin terms', msg
  !   call wrtout(std_out,msg,'COLL')

  !   write(msg, '(1X, A16, 2X, I5.3)') 'Number of terms:',  self%coeff
  !   call wrtout(std_out,msg,'COLL')

  !   do i=1, self%total_nnz
  !      do ii=1,3
  !         R(ii)=self%total_Rlist(ii)%data(i)
  !         do jj=1, 3
  !            tmp(jj, ii)=self%total_val_list(jj, ii)%data(i)
  !         enddo
  !      enddo

  !      msg=repeat("-", 64)
  !      write(msg, '(1X, A64)') msg
  !      call wrtout(std_out,msg,'COLL')
  !      call wrtout(ab_out, msg, 'COLL')

  !      write(msg, "(2X, A3, I5.4, 2X, A3, I5.4, 5X, A3, I4.2, I5.2, I5.2)")  & 
  !           'i =', self%total_ilist%data(i), 'j =', self%total_jlist%data(i), 'R =', R
  !      call wrtout(std_out,msg,'COLL')
  !      call wrtout(ab_out, msg, 'COLL')
  !      write(msg, "(2X, A16)") 'Matrix Form (Ha)'
  !      call wrtout(std_out,msg,'COLL')
  !      call wrtout(ab_out, msg, 'COLL')
  !      do ii=1, 3
  !         write(msg, "(3ES21.12)") tmp(ii, :)
  !         call wrtout(std_out,msg,'COLL')
  !         call wrtout(ab_out, msg, 'COLL')
  !      end do
  !   end do
  !   msg=repeat("=", 80)
  !   call wrtout(std_out,msg,'COLL')
  !   call wrtout(ab_out, msg, 'COLL')
  end subroutine print_terms

end module m_spin_primitive_potential
