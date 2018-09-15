!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_model_primitive
!! NAME
!! m_spin_model_primitive
!!
!! FUNCTION
!! This module contains the atomic structures and the spin hamiltonian inside the primitive cell
!! which can be directly mapped to the xml file. It is not for the calculation, but for constructing
!! the hamiltonian in supercell. It is also used as input for the magnon band structure calculation.
!!
!! Datatypes:
!!  spin_model_primitive_t
!!
!! Subroutines:
!! 
!!  * spin_model_primitive_t_initialize
!!  * spin_model_primitive_t_read_xml
!!  * spin_model_primitive_t_make_supercell
!!  * spin_model_primitive_t_finalize
!!
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
#include "abi_common.h"

module m_spin_model_primitive
  use iso_c_binding
  use m_dynmaic_array
  use m_mathfuncs
  use defs_basis
  use m_abicore
  use m_errors
  use m_supercell
  use m_spin_terms
  implicit none

!!*** 
  interface
     ! C function:
     ! void xml_read_spin(char *fname, double *ref_energy, double *unitcell[9],
     ! int *natoms, double *masses[], int *nmatoms,
     ! int *index_spin[], double *gyroratios[], double *damping_factors[],
     ! double *positions[], double *spinat[],
     ! // exchange
     ! int *exc_nnz, int *exc_ilist[],
     ! int *exc_jlist[], int *exc_Rlist[],
     ! double *exc_vallist[],
     ! //dmi
     ! int *dmi_nnz, int *dmi_ilist[],
     ! int *dmi_jlist[], int *dmi_Rlist[],
     ! double *dmi_vallist[],
     ! //uniaxial SIA
     ! int *uni_nnz, int *uni_ilist[],
     ! double *uni_amplitude_list[],
     ! double *uni_direction_list[],
     ! //bilinear
     ! int *bi_nnz, int *bi_ilist[],
     ! int *bi_jlist[], int *bi_Rlist[],
     ! double *bi_vallist[])

     subroutine xml_read_spin(xml_fname, ref_energy, unitcell,                 &
          natoms, masses, nmatoms, index_spin, gyroratios, damping_factors, positions, spinat, &
          exc_nnz, exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
          dmi_nnz, dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
          uni_nnz, uni_ilist, uni_amplitude_list, uni_direction_list, &
          bi_nnz, bi_ilist, bi_jilst, bi_Rlist, bi_vallist) bind(C, name="xml_read_spin")
       import
       character(len=1), dimension(*), intent(in) :: xml_fname
       integer (c_int), intent(out):: natoms, nmatoms, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
       real  (c_double), intent(out) :: ref_energy
       type(c_ptr)::  unitcell,  &
            masses,  index_spin, gyroratios, damping_factors, positions, spinat, &
            exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
            dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
            uni_ilist, uni_amplitude_list, uni_direction_list, &
            bi_ilist, bi_jilst, bi_Rlist, bi_vallist
     end subroutine xml_read_spin
  end interface

  type spin_model_primitive_t
     integer :: natoms, nmatoms, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
     real (dp) :: ref_energy, unitcell(3,3)
     ! integer, allocatable :: masses,  index_spin, gyroratios, damping_factors, positions, spinat, &
     !exc_ilist, exc_jlist, exc_Rlist, exc_vallist, &
     !dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist, &
     !uni_ilist, uni_amplitude_list, uni_direction_list, &
     !bi_ilist, bi_jilst, bi_Rlist, bi_vallist
     integer, allocatable :: index_spin(:), &
          exc_ilist(:), exc_jlist(:), exc_Rlist(:,:), &
          dmi_ilist(:), dmi_jlist(:), dmi_Rlist(:,:), &
          uni_ilist(:), &
          bi_ilist(:) ,bi_jlist(:), bi_Rlist(:,:)
     real(dp), allocatable ::gyroratios(:), damping_factors(:), &
          positions(:,:), spinat(:,:), &
          exc_vallist(:,:), dmi_vallist(:,:), &
          uni_amplitude_list(:), uni_direction_list(:,:), &
          bi_vallist(:,:,:)
     ! TOTAL
     integer :: total_nnz=0
     type(int_array_type):: total_ilist, total_jlist, total_Rlist(3)
     type(real_array_type) :: total_val_list(3,3)

     !  contains
     !    procedure:: initialize=> spin_model_primitive_t_initialize
     !    procedure:: finalize=> spin_model_primitive_t_finalize
     !    procedure:: set_atoms => spin_model_primitive_t_set_atoms
     !    procedure:: set_bilinear => spin_model_primitive_t_set_bilinear
     !    procedure:: set_exchange=> spin_model_primitive_t_set_exchange
     !    procedure:: set_dmi => spin_model_primitive_t_set_dmi
     !    procedure:: set_uni => spin_model_primitive_t_set_uni
     !    procedure:: read_xml => spin_model_primitive_t_read_xml
     !    procedure:: make_supercell => spin_model_primitive_t_make_supercell
     !    procedure :: print_terms => spin_model_primitive_t_print_terms
     !    procedure:: get_total_terms => spin_model_primitive_t_get_total_term
  end type spin_model_primitive_t

contains

  subroutine spin_model_primitive_t_initialize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_initialize'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
   !TODO should something  be done here?
  end subroutine spin_model_primitive_t_initialize

  subroutine spin_model_primitive_t_set_atoms(self, natoms, unitcell, positions, &
       nmatoms, index_spin, spinat, gyroratios, damping_factors )


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_set_atoms'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    integer, intent(in):: natoms, nmatoms, index_spin(:)
    real(dp), intent(in):: unitcell(3, 3),  positions(3,natoms), &
         spinat(3,nmatoms), gyroratios(nmatoms), damping_factors(nmatoms)

    !print *, "natoms",natoms
    !print *, "nmatoms", nmatoms
    !print *, "positions", positions
    !print *, "spinat", spinat
    !print *, "gyroratios", gyroratios
    !print *, "damping_factors", damping_factors
    ABI_ALLOCATE(self%positions, (3, natoms))
    ABI_ALLOCATE(self%index_spin, (natoms))
    ABI_ALLOCATE(self%spinat, (3, nmatoms))
    ABI_ALLOCATE(self%gyroratios, (nmatoms))
    ABI_ALLOCATE(self%damping_factors, (nmatoms))

    self%natoms=natoms
    self%unitcell(:,:)=unitcell(:,:)
    self%positions(:,:)=positions(:,:)
    self%nmatoms=nmatoms
    self%index_spin(:)=index_spin(:)
    self%spinat(:,:)=spinat(:,:)
    self%gyroratios(:)=gyroratios(:)
    self%damping_factors(:)=damping_factors(:)

  end subroutine spin_model_primitive_t_set_atoms


  subroutine  spin_model_primitive_t_set_bilinear(self, n, ilist, jlist, Rlist, vallist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_set_bilinear'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:, :,:)
    integer :: idx, ii, ix, iy
    self%total_nnz=self%total_nnz+n
    do idx = 1, n, 1
       !call self%total_ilist%push(ilist(idx))
       !call self%total_jlist%push(jlist(idx))
       call int_array_type_push(self%total_ilist, ilist(idx))
       call int_array_type_push(self%total_jlist, jlist(idx))
       do ii = 1, 3, 1
          !call self%total_Rlist(ii)%push(Rlist(ii, idx))
          call int_array_type_push(self%total_Rlist(ii), Rlist(ii, idx))
       end do
       do iy = 1, 3, 1
          do ix = 1, 3, 1
             !call self%total_val_list(ix, iy)%push(vallist(ix, iy, idx))
             call real_array_type_push( self%total_val_list(ix, iy), vallist(ix, iy, idx))
          end do
       end do
    end do
  end subroutine spin_model_primitive_t_set_bilinear

  subroutine spin_model_primitive_t_set_exchange(self, n, ilist, jlist, Rlist, vallist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_set_exchange'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
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
    !call self%set_bilinear(n,ilist,jlist,Rlist,bivallist)
    call  spin_model_primitive_t_set_bilinear(self,n,ilist,jlist,Rlist,bivallist)
  end subroutine spin_model_primitive_t_set_exchange


  subroutine spin_model_primitive_t_set_dmi(self, n, ilist, jlist, Rlist, vallist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_set_dmi'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: idx
    real(dp) :: bivallist(3,3, n), D(3)
    bivallist(:,:,:)=0.0d0
    ! -2.0 * self.Ku[i]) * np.outer(self.e[i], self.e[i]
    do idx=1,n, 1
       D(:)=vallist(:, idx)
       ! 0 Dz -Dy
       ! -Dz 0 Dx
       ! Dy -Dx 0
       bivallist(:,:, idx)=reshape ( (/0.0d0, -D(3), D(2),  &
            D(3), 0.0d0, -D(1),  &
            -D(2), D(1), 0.0d0 /),(/3,3/) )
    end do
    !call self%set_bilinear(n,ilist,jlist,Rlist,bivallist)
    call  spin_model_primitive_t_set_bilinear(self,n,ilist,jlist,Rlist,bivallist)
  end subroutine spin_model_primitive_t_set_dmi

  subroutine spin_model_primitive_t_set_uni(self, n, ilist, k1list, k1dirlist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_set_uni'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    integer, intent(in) :: n, ilist(:)
    real(dp), intent(in) :: k1list(:), k1dirlist(:, :)
    integer :: idx, Rlist(3, n)
    real(dp) :: bivallist(3,3, n)
    bivallist(:,:,:)=0.0d0
    Rlist(:, :)=0.0d0
    ! at
    ! -2.0 * self.Ku[i]) * np.outer(self.e[i], self.e[i]
    do idx=1,n, 1
       bivallist(:,:, idx)= (-2.0d0* k1list(idx))*  &
            outer_product(k1dirlist(:,idx), k1dirlist(:, idx))
    end do
    !call self%set_bilinear(n,ilist,ilist,Rlist,bivallist)
    call  spin_model_primitive_t_set_bilinear(self,n,ilist,ilist,Rlist,bivallist)
  end subroutine spin_model_primitive_t_set_uni

  subroutine spin_model_primitive_t_read_xml(self, xml_fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_read_xml'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    character(len=*), intent(in):: xml_fname
    integer :: natoms, nmatoms, exc_nnz, dmi_nnz, uni_nnz, bi_nnz
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

    call xml_read_spin(trim(xml_fname)//C_NULL_CHAR, ref_energy, p_unitcell,                 &
         natoms, p_masses, nmatoms, p_index_spin, p_gyroratios, p_damping_factors, p_positions, p_spinat, &
         exc_nnz, p_exc_ilist, p_exc_jlist, p_exc_Rlist, p_exc_vallist, &
         dmi_nnz, p_dmi_ilist, p_dmi_jlist, p_dmi_Rlist, p_dmi_vallist, &
         uni_nnz, p_uni_ilist, p_uni_amplitude_list, p_uni_direction_list, &
         bi_nnz, p_bi_ilist, p_bi_jlist, p_bi_Rlist, p_bi_vallist)

    call c_f_pointer(p_unitcell, unitcell, [9])
    call c_f_pointer(p_masses, masses, [natoms])
    call c_f_pointer(p_index_spin, index_spin, [natoms])
    call c_f_pointer(p_gyroratios, gyroratios, [nmatoms])
    call c_f_pointer(p_damping_factors, damping_factors, [nmatoms])
    call c_f_pointer(p_positions, positions, [natoms*3])
    call c_f_pointer(p_spinat, spinat, [nmatoms*3])
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

    print *, "Spin model: setting structure."
    !print *, "natoms", natoms
    !print *, "nmatoms", nmatoms
    !print *, "positions", positions
    !print *, "unitcell: ", unitcell
    !print *, "spinat", spinat
    !print *, "index_spin", index_spin
    call spin_model_primitive_t_set_atoms(self,natoms,reshape(unitcell, [3,3]), & 
            & reshape(positions, [3, natoms]), &
            & nmatoms, &
            & index_spin, &
            & reshape(spinat, [3, nmatoms]), &
            & gyroratios,damping_factors)

    print *, "Spin model: setting exchange terms."
    ! call self%set_exchange(exc_nnz,exc_ilist,exc_jlist,&
    !      reshape(exc_Rlist, (/3, exc_nnz /)), &
    !      reshape(exc_vallist, (/3, exc_nnz/)) )
    call spin_model_primitive_t_set_exchange(self, exc_nnz,exc_ilist,exc_jlist,&
         reshape(exc_Rlist, (/3, exc_nnz /)), &
         reshape(exc_vallist, (/3, exc_nnz/)))

    print *, "Spin model: setting dmi terms."
    ! call self%set_dmi( n=dmi_nnz, ilist=dmi_ilist, jlist=dmi_jlist, &
    !    Rlist=reshape(dmi_Rlist, (/3, dmi_nnz /)), &
    !     vallist = reshape(dmi_vallist, (/3, dmi_nnz/)) )
    call spin_model_primitive_t_set_dmi(self, n=dmi_nnz, ilist=dmi_ilist, jlist=dmi_jlist, &
         Rlist=reshape(dmi_Rlist, (/3, dmi_nnz /)), &
         vallist = reshape(dmi_vallist, (/3, dmi_nnz/)))
    print *, "Spin model: setting uniaxial SIA terms."

    !call self%set_uni(uni_nnz, uni_ilist, uni_amplitude_list, &
    !     reshape(uni_direction_list, [3, uni_nnz]) )
    call spin_model_primitive_t_set_uni(self, uni_nnz, uni_ilist, uni_amplitude_list, &
         reshape(uni_direction_list, [3, uni_nnz]) )

    print *, "Spin model: setting bilinear terms."
    !call self%set_bilinear(bi_nnz, bi_ilist, bi_jlist,  &
    !     Rlist=reshape(bi_Rlist, (/3, bi_nnz /)), &
    !     vallist = reshape(bi_vallist, (/3,3, bi_nnz/)) )
    call spin_model_primitive_t_set_bilinear(self, bi_nnz, bi_ilist, bi_jlist,  &
         Rlist=reshape(bi_Rlist, (/3, bi_nnz /)), &
         vallist = reshape(bi_vallist, (/3,3, bi_nnz/)))


    ! TODO hexu: should use free in C code, not here.
    if(associated(unitcell))  then
       ABI_DEALLOCATE(unitcell)
    end if
    nullify(unitcell)

    if(associated(masses))  then
       ABI_DEALLOCATE(masses)
    end if
    nullify(masses)

    if(associated(index_spin))  then
       ABI_DEALLOCATE(index_spin)
    end if
    nullify(index_spin)

    if(associated(gyroratios))  then
       ABI_DEALLOCATE(gyroratios)
    end if
    nullify(gyroratios)

    if(associated(damping_factors))  then
       ABI_DEALLOCATE(damping_factors)
    end if
    nullify(damping_factors)
    if(associated(positions))  then
       ABI_DEALLOCATE(positions)
    end if
    nullify(positions)

    if(associated(spinat))  then
       ABI_DEALLOCATE(spinat)
    end if
    nullify(spinat)

    if(exc_nnz /=0) then
       if(associated(exc_ilist))  then
          ABI_DEALLOCATE(exc_ilist)
       end if
       nullify(exc_ilist)


       if(associated(exc_jlist))  then
          ABI_DEALLOCATE(exc_jlist)
       end if
       nullify(exc_jlist)

       if(associated(exc_Rlist))  then
          ABI_DEALLOCATE(exc_Rlist)
       end if
       nullify(exc_Rlist)

       if(associated(exc_vallist))  then
          ABI_DEALLOCATE(exc_vallist)
       end if
       nullify(exc_vallist)
    endif

    if(dmi_nnz/=0) then
       if(associated(dmi_ilist))  then
          ABI_DEALLOCATE(dmi_ilist)
       end if
       nullify(dmi_ilist)

       if(associated(dmi_jlist))  then
          ABI_DEALLOCATE(dmi_jlist)
       end if
       nullify(dmi_jlist)

       if(associated(dmi_Rlist))  then
          ABI_DEALLOCATE(dmi_Rlist)
       end if
       nullify(dmi_Rlist)

       if(associated(dmi_vallist))  then
          ABI_DEALLOCATE(dmi_vallist)
       end if
       nullify(dmi_vallist)
    end if

    if(uni_nnz/=0) then
       if(associated(uni_ilist))  then
          ABI_DEALLOCATE(uni_ilist)
       end if
       nullify(uni_ilist)

       if(associated(uni_amplitude_list))  then
          ABI_DEALLOCATE(uni_amplitude_list)
       end if
       nullify(uni_amplitude_list)

       if(associated(uni_direction_list))  then
          ABI_DEALLOCATE(uni_direction_list)
       end if
       nullify(uni_direction_list)
    end if

    if(bi_nnz/=0) then
       if(associated(bi_ilist))  then
          ABI_DEALLOCATE(bi_ilist)
       end if
       nullify(bi_ilist)

       if(associated(bi_jlist))  then
          ABI_DEALLOCATE(bi_jlist)
       end if
       nullify(bi_jlist)

       if(associated(bi_Rlist))  then
          ABI_DEALLOCATE(bi_Rlist)
       end if
       nullify(bi_Rlist)

       if(associated(bi_vallist))  then
          ABI_DEALLOCATE(bi_vallist)
       end if
       nullify(bi_vallist)
    endif
  end subroutine spin_model_primitive_t_read_xml

  subroutine spin_model_primitive_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_finalize'
!End of the abilint section

    class(spin_model_primitive_t), intent(inout) :: self
    integer :: i, j

    !if(allocated(self%masses))  then
    !    ABI_DEALLOCATE(self%masses)
    !  end if
    if(allocated(self%index_spin))  then
       ABI_DEALLOCATE(self%index_spin)
    end if

    if(allocated(self%gyroratios))  then
       ABI_DEALLOCATE(self%gyroratios)
    end if

    if(allocated(self%damping_factors))  then
       ABI_DEALLOCATE(self%damping_factors)
    end if

    if(allocated(self%positions))  then
       ABI_DEALLOCATE(self%positions)
    end if

    if(allocated(self%spinat))  then
       ABI_DEALLOCATE(self%spinat)
    end if

    if(allocated(self%exc_ilist))  then
       ABI_DEALLOCATE(self%exc_ilist)
    end if

    if(allocated(self%exc_jlist))  then
       ABI_DEALLOCATE(self%exc_jlist)
    end if

    if(allocated(self%exc_Rlist))  then
       ABI_DEALLOCATE(self%exc_Rlist)
    end if

    if(allocated(self%exc_vallist))  then
       ABI_DEALLOCATE(self%exc_vallist)
    end if

    if(allocated(self%dmi_ilist))  then
       ABI_DEALLOCATE(self%dmi_ilist)
    end if

    if(allocated(self%dmi_jlist))  then
       ABI_DEALLOCATE(self%dmi_jlist)
    end if

    if(allocated(self%dmi_Rlist))  then
       ABI_DEALLOCATE(self%dmi_Rlist)
    end if

    if(allocated(self%dmi_vallist))  then
       ABI_DEALLOCATE(self%dmi_vallist)
    end if

    if(allocated(self%uni_ilist))  then
       ABI_DEALLOCATE(self%uni_ilist)
    end if

    if(allocated(self%uni_amplitude_list))  then
       ABI_DEALLOCATE(self%uni_amplitude_list)
    end if

    if(allocated(self%uni_direction_list))  then
       ABI_DEALLOCATE(self%uni_direction_list)
    end if

    if(allocated(self%bi_ilist))  then
       ABI_DEALLOCATE(self%bi_ilist)
    end if

    if(allocated(self%bi_jlist))  then
       ABI_DEALLOCATE(self%bi_jlist)
    end if

    if(allocated(self%bi_Rlist))  then
       ABI_DEALLOCATE(self%bi_Rlist)
    end if

    if(allocated(self%bi_vallist))  then
       ABI_DEALLOCATE(self%bi_vallist)
    end if

    !call self%total_ilist%finalize()
    !call self%total_jlist%finalize()
    call int_array_type_finalize(self%total_ilist)
    call int_array_type_finalize(self%total_jlist)

    do i=1, 3, 1
       !call self%total_Rlist(i)%finalize()
       call int_array_type_finalize(self%total_Rlist(i))
    enddo

    do i=1, 3, 1
       do j=1, 3, 1
          !call self%total_val_list(i, j)%finalize()
          call real_array_type_finalize(self%total_val_list(j, i))
       end do
    end do

  end subroutine spin_model_primitive_t_finalize

  subroutine spin_ham_set_exchange(self, nnz,  ilist, jlist, Rlist, vallist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_ham_set_exchange'
!End of the abilint section

    type(spin_model_primitive_t) , intent(inout) :: self
    integer, intent(in) :: nnz,  ilist(:), jlist(:), Rlist(:,:)
    real(dp), intent(in) :: vallist(:,:)
    integer :: err
    allocate(self%exc_ilist(nnz), stat=err)
    allocate(self%exc_jlist(nnz), stat=err)
    allocate(self%exc_Rlist(3,nnz), stat=err)
    allocate(self%exc_vallist(3,nnz), stat=err)
    self%exc_nnz=nnz
    self%exc_ilist=ilist
    self%exc_jlist=jlist
    self%exc_Rlist=Rlist
    self%exc_vallist=vallist
  end subroutine spin_ham_set_exchange

  ! R (in term of primitive cell) to R_sc(in term of supercell) + R_prim
  subroutine find_R_PBC(scell, R, R_sc, R_prim)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_R_PBC'
!End of the abilint section

    type(supercell_type) , intent(in):: scell
    integer, intent(in):: R(3)
    integer, intent(out):: R_sc(3), R_prim(3)
    real(dp) :: R_sc_d(3), sc_mat(3,3)

    integer:: ipriv(3), info
    !call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
    sc_mat(:,:)=scell%rlatt
    R_sc_d(:)=R(:)
    !print *, sc_mat
    !print *, R_sc_d
    call dgesv(3, 1, sc_mat, 3, ipriv, R_sc_d, 3, info)
    if ( info/=0 ) then
       print *, "Failed to find R_sc"
    end if

    ! if only diagonal of rlatt works.
    !R_sc_d(1)= real(R(1))/real(scell%rlatt(1,1))
    !R_sc_d(2)= real(R(2))/real(scell%rlatt(2,2))
    !R_sc_d(3)= real(R(3))/real(scell%rlatt(3,3))
    ! TODO hexu: R_prim should be non-negative, which is assumed in m_supercell.
    ! But should we make it more general?
    R_sc(1)=floor(R_sc_d(1))
    R_sc(2)=floor(R_sc_d(2))
    R_sc(3)=floor(R_sc_d(3))
    R_prim(1)=(R(1)-R_sc(1)*scell%rlatt(1,1))
    R_prim(2)=(R(2)-R_sc(2)*scell%rlatt(2,2))
    R_prim(3)=(R(3)-R_sc(3)*scell%rlatt(3,3))
  end subroutine find_R_PBC

  ! TODO hexu: move this to m_supercell?
  ! find the spercelll atom index from index of atom in primitive cell and R vector
  function find_supercell_index(scell, iatom_prim, rvec) result(iatom_supercell)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_supercell_index'
!End of the abilint section

    type(supercell_type) , intent(in):: scell
    integer, intent(in) :: iatom_prim, rvec(3)
    integer  :: iatom_supercell
    integer :: i
    iatom_supercell=-1
    do i=1, scell%natom, 1
       if ( scell%atom_indexing(i) == iatom_prim .and. &
            all(scell%uc_indexing(:,i)==rvec) ) then
          iatom_supercell=i
          return
       end if
    end do
    print *, "cannot find iatom_prim, rvec pair", iatom_prim, rvec, "in supercell"
  end function find_supercell_index

  !i0, j0+R0 shifted by R to i1=i0+0+R->periodic, j1=j0+R0+R->periodic
  subroutine find_supercell_ijR(scell, i0, j0, R0, R, i1, j1, R1, R_sc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_supercell_ijR'
!End of the abilint section

    type(supercell_type) , intent(in):: scell
    integer, intent(in) :: i0, j0, R0(3), R(3)
    integer, intent(out) :: i1, j1, R1(3), R_sc(3)
    i1=find_supercell_index(scell,i0,R)
    call find_R_PBC(scell,R0+R,R_sc,R1)
    j1=find_supercell_index(scell, j0, R1)
  end subroutine find_supercell_ijR

  subroutine spin_model_primitive_t_make_supercell(self, sc_matrix, sc_ham)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_make_supercell'
!End of the abilint section

    class(spin_model_primitive_t) , intent(in) :: self
    type(spin_terms_t) , intent(inout) :: sc_ham
    integer :: sc_matrix(3,3), atom_index(self%nmatoms)

    integer ::  typat_primcell(self%natoms), sc_nmatoms,  i, counter, icell
    real(dp) :: znucl(self%natoms), tmp(3,3)
    type(supercell_type) :: scell
    integer, allocatable ::sc_index_spin(:), sc_znucl(:)
    integer, allocatable :: sc_ispin_prim(:), sc_rvec(:, :)
    real(dp), allocatable ::sc_spinat(:,:), sc_gyroratios(:), sc_damping_factors(:), sc_spinpos(:,:)
    integer :: ii, jj, icol, irow, rr(3), R_sc(3), iatom

    typat_primcell(:)=1
    ! TODO hexu:should use supercell generate with lattice model
    !init_supercell(natom_primcell, sc_matrix, rprimd_primcell,
    !             typat_primcell, xcart_primcell, znucl, scell)

    znucl(:)=0

    counter=0
    do i=1, self%natoms
      if(self%index_spin(i)>0) then
          counter=counter+1
          atom_index(counter)=i
      endif
    enddo 

    call init_supercell(self%natoms, sc_matrix, self%unitcell, typat_primcell, &
         self%positions, znucl, scell)
    ! cell, positions
    !scell%rprimd
    !scell%xcart
    !nmatoms
    sc_nmatoms=scell%ncells*self%nmatoms
    ABI_ALLOCATE(sc_index_spin, (self%natoms*scell%ncells))
    ABI_ALLOCATE(sc_ispin_prim, (sc_nmatoms) )
    ABI_ALLOCATE(sc_rvec, (3, sc_nmatoms) )
    ABI_ALLOCATE(sc_spinat, (3, sc_nmatoms) )
    ABI_ALLOCATE(sc_spinpos, (3, sc_nmatoms) )
    ABI_ALLOCATE(sc_gyroratios, (sc_nmatoms))
    ABI_ALLOCATE(sc_damping_factors, (sc_nmatoms))
    ABI_ALLOCATE(sc_znucl, (sc_nmatoms))
    ! sc_index_spin
    counter=0
    do i = 1, scell%natom, 1
       iatom=scell%atom_indexing(i)
       !if(scell%atom_indexing(i)>0) then
       if(self%index_spin(iatom)>0) then
          counter=counter+1
          sc_index_spin(i)=counter
          sc_spinpos(:, counter)=scell%xcart(:,counter)
          sc_spinat(:, counter)=self%spinat(:,self%index_spin(iatom))
          sc_gyroratios( counter)=self%gyroratios(self%index_spin(iatom))
          sc_damping_factors( counter)=self%damping_factors(self%index_spin(iatom))
          sc_ispin_prim(counter) = self%index_spin(iatom)
          sc_rvec(:,counter)=scell%uc_indexing(:,i)
       else
          sc_index_spin(i)=-1
       endif
    end do
    !print *, sc_index_spin

    do i=1, sc_nmatoms
       sc_znucl(i)=1
    enddo

    !call spin_terms_t_initialize(sc_ham,)
    !call sc_ham%initialize( cell=scell%rprimd, pos=sc_spinpos, &
    !     spinat=sc_spinat, zion=sc_znucl)
    call spin_terms_t_initialize(sc_ham, cell=scell%rprimd, pos=sc_spinpos, &
         spinat=sc_spinat, zion=sc_znucl, spin_index=sc_index_spin, ispin_prim=sc_ispin_prim, rvec=sc_rvec)
    sc_ham%gyro_ratio(:)=sc_gyroratios(:)
    sc_ham%gilbert_damping(:)=sc_damping_factors(:)

    !print *, "The total number of terms in primitive cell: ", self%total_nnz
    do i =1, self%total_nnz, 1
       do icell=1, scell%ncells, 1
          ! Note i0 and j0 are in spin index, while find_supercell_ijR work in atom index
          call find_supercell_ijR(scell=scell, i0=atom_index(self%total_ilist%data(i)), &
               j0=atom_index(self%total_jlist%data(i)), &
               R0=[self%total_Rlist(1)%data(i),  &
               self%total_Rlist(2)%data(i),  &
               self%total_Rlist(3)%data(i)], &
               R=scell%rvecs(:,icell), &
               i1=ii,j1=jj,R1=rr,R_sc=R_sc)
          ! ii is the index in atom supercell not in spin supercell 
          !scell%rvecs(icell)
          do irow = 1, 3
             do icol=1, 3
                tmp(icol, irow)=self%total_val_list(icol,irow)%data(i)
             end do
          end do
          !print *, "i0", self%total_ilist%data(i)
          !print *, "ii", ii, sc_index_spin(ii)
          !print *, "jj", jj, sc_index_spin(jj)
          call spin_terms_t_set_bilinear_term_single(sc_ham, sc_index_spin(ii), sc_index_spin(jj), tmp)
       enddo
    enddo

    if (allocated(sc_ispin_prim)) then
       ABI_DEALLOCATE(sc_ispin_prim)
    endif
    if (allocated(sc_rvec)) then
       ABI_DEALLOCATE(sc_rvec)
    endif
    if (allocated(sc_index_spin)) then
       ABI_DEALLOCATE(sc_index_spin)
    endif

    if (allocated(sc_gyroratios)) then
       ABI_DEALLOCATE(sc_gyroratios)
    endif
    if (allocated(sc_damping_factors)) then
       ABI_DEALLOCATE(sc_damping_factors)
    end if
    if (allocated(sc_spinat)) then
       ABI_DEALLOCATE(sc_spinat)
    endif
    ! TODO hexu: should it be destroyed here? or do all magnetic moment should be stored there?
    call destroy_supercell(scell)
  end subroutine spin_model_primitive_t_make_supercell

  subroutine spin_model_primitive_t_print_terms(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_model_primitive_t_print_terms'
!End of the abilint section

    class(spin_model_primitive_t) :: self
    integer :: i, ii, jj,  R(3)
    real(dp) :: tmp(3, 3)
    print *, "==========spin terms==============="
    print *, "Number of terms: ",  self%total_nnz
    do i=1, self%total_nnz
       do ii=1,3
          R(ii)=self%total_Rlist(ii)%data(i)
          do jj=1, 3
             tmp(jj, ii)=self%total_val_list(jj, ii)%data(i)
          enddo
       enddo
       print *, 'i=',  self%total_ilist%data(i), '  j=', self%total_jlist%data(i), &
            '  R=', R , ' Value= '
       do ii=1, 3
          write(*, *) tmp(ii, :)
       end do
    end do
    print *, "==================================="
  end subroutine spin_model_primitive_t_print_terms

end module m_spin_model_primitive
