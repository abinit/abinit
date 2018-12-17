!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_terms
!! NAME
!! m_spin_terms
!!
!! FUNCTION
!! This module contains the spin hamiltonian, and the methods for
!! calculating effective magnetic field (torque), dS/dt, and total_energy
!!
!!
!! Datatypes:
!!
!! * spin_terms_t
!!
!! Subroutines:
!!
!! * spin_terms_t_initialize
!! * spin_terms_t_finalize
!! * spin_terms_t_total_Heff : calculate total Heff (no Langevin term)
!! * spin_terms_t_Heff_to_dSdt: 
!!  * spin_terms_t_get_dSdt : dSdt, Langevin term is an input.
!!  * spin_terms_t_get_Langevin_Heff

!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (TO, hexu)
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
module  m_spin_terms
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_multibinit_global
  use m_spmat_csr, only : CSR_mat_t
  use m_spmat_lil, only : LIL_mat_t
  use m_spmat_convert, only : LIL_to_CSR
  use m_effpot_api, only : effpot_t
  implicit none
  !!***

  type, extends(effpot_t) :: spin_terms_t
     real(dp) :: etot
     ! ispin_prim: index in the spin model in primitive cell, which is used as the index of sublattice.
     ! rvec: supercell R vector.
     integer, allocatable :: ispin_prim(:), rvec(:,:)
     real(dp), allocatable ::  pos(:,:),  spinat(:,:), S(:,:)
     real(dp):: cell(3,3)
     ! index of atoms in the spin hamiltonian. -1 if the atom is not in the spin_hamiltonian.
     integer, allocatable :: iatoms(:)
     !integer :: seed
     ! random seed
     logical :: has_external_hfield, has_uniaxial_anistropy, has_exchange, &
          has_DMI, has_dipdip, has_bilinear


     ! Array or scalar?
     real(dp), allocatable :: external_hfield(:,:)

     ! Exchange/DMI/dipdip stored like COO sparse matrix form.

     ! nint: total number of i,j pairs. size nint
     integer :: exchange_nint
     integer, allocatable:: exchange_i(:)
     integer, allocatable:: exchange_j(:)
     real(dp), allocatable:: exchange_val(:)


     ! nint: total number of i,j pairs. size nint
     integer :: DMI_nint
     integer, allocatable:: DMI_i(:)
     integer, allocatable:: DMI_j(:)
     ! DMI item should be vector. (3*nint)
     real(dp), allocatable:: DMI_val(:,:)

     ! dipole-dipole. hexu: Do we use self?
     integer :: dipdip_nint
     integer, allocatable:: dipdip_i(:)
     integer, allocatable:: dipdip_j(:)
     real(dp), allocatable:: dipdip_val(:,:,:)

     ! Matrix form of bilinear interactions.
     !integer :: bilinear_nint
     !integer, allocatable:: bilinear_i(:)
     !integer, allocatable:: bilinear_j(:)
     !real(dp), allocatable:: bilinear_val(:,:,:)
     type(LIL_mat_t) :: bilinear_lil_mat
     logical :: csr_mat_ready= .False.
     type(CSR_mat_t) :: bilinear_csr_mat
     ! 3, 3, ninit

     real(dp), allocatable :: gyro_ratio(:)
     !     Gyromagnetic ration can be on-site
     !     Size should be nspins

     real(dp), allocatable :: k1(:)
     !     On-site uniaxial anisotropy
     !     Size should be nspins

     real(dp), allocatable :: k1dir(:,:)
     !     Direction of On-site uniaxial anisotropy
     !     Size should be 3, nspins
     !   Should be normalized when being set.

     real(dp), allocatable :: gilbert_damping(:)
     ! Heat bath
     !     Coupling to thermal bath can be on-site
     !     Size should be nspins?


     ! vector for calculating effective field
     real(dp), allocatable :: Htmp(:, :)
   CONTAINS
     procedure, non_overridable :: initialize => spin_terms_t_initialize
     procedure, non_overridable :: finalize => spin_terms_t_finalize
     procedure, non_overridable :: get_Heff => spin_terms_t_total_Heff
     procedure, non_overridable :: calculate => spin_terms_t_calculate
     procedure, non_overridable :: get_energy => spin_terms_t_get_energy
     procedure, non_overridable :: get_delta_E => spin_terms_t_get_delta_E
     procedure, non_overridable :: set_bilinear_term => spin_terms_t_set_bilinear_term
  end type spin_terms_t

contains
  subroutine spin_terms_t_initialize(self, cell, pos, spinat,  iatoms, &
       & ispin_prim, rvec, gyro_ratio, damping)

    implicit none
    !Arguments ------------------------------------
    !scalars
    class(spin_terms_t), intent(out) :: self
    integer, intent(in) :: ispin_prim(:), rvec(:,:)
    integer, intent(in) ::  iatoms(:)
    real(dp), intent(in) :: cell(3,3), pos(:,:), spinat(:,:), gyro_ratio(:), damping(:)
    !Local variables-------------------------------
    integer :: nspins,  i

    self%has_spin=.True.
    self%has_displacement=.False.
    self%has_strain=.False.
    self%is_null=.False.

    nspins=size(pos, 2)
    self%nspins=nspins
    call xmpi_bcast(nspins, master, comm, ierr)
    call xmpi_bcast(self%nspins, master, comm, ierr)

    ABI_ALLOCATE(self%iatoms, (nspins))
    self%iatoms(:)=iatoms(:)
    call xmpi_bcast(self%iatoms, master, comm, ierr)
    self%cell(:,:)=cell(:,:)
    call xmpi_bcast(self%cell, master, comm, ierr)

    ABI_ALLOCATE( self%pos, (3,nspins) )
    self%pos(:,:)=pos(:,:)
    call xmpi_bcast(self%pos, master, comm, ierr)

    ABI_ALLOCATE( self%spinat, (3, nspins))
    self%spinat(:,:)=spinat(:,:)
    call xmpi_bcast(self%spinat, master, comm, ierr)

    ABI_ALLOCATE( self%ms, (nspins) )
    self%ms(:)= sqrt(sum(spinat(:,:)**2, dim=1))* mu_B
    call xmpi_bcast(self%ms, master, comm, ierr)

    ABI_ALLOCATE( self%ispin_prim, (nspins))
    ABI_ALLOCATE(self%rvec, (3, nspins))
    
    self%ispin_prim(:)=ispin_prim(:)

    call xmpi_bcast(self%ispin_prim, master, comm, ierr)
    self%rvec(:,:)=rvec(:,:)
    call xmpi_bcast(self%rvec, master, comm, ierr)

    ABI_ALLOCATE( self%S, (3, nspins))
    do i=1,nspins
       self%S(:,i)=self%spinat(:,i)/self%ms(i)
    end do
    call xmpi_bcast(self%S, master, comm, ierr)
    
    self%has_external_hfield=.False.
    !self%has_uniaxial_anistropy=.False.
    !self%has_exchange=.False.
    !self%has_DMI=.False.
    self%has_dipdip=.False.
    self%has_bilinear=.False.
    call xmpi_bcast(self%has_external_hfield, master, comm, ierr)
    call xmpi_bcast(self%has_dipdip, master, comm, ierr)
    call xmpi_bcast(self%has_bilinear, master, comm, ierr)

    ABI_ALLOCATE( self%gyro_ratio, (nspins))
    ABI_ALLOCATE( self%gilbert_damping, (nspins) )
    ! Defautl gyro_ratio
    self%gyro_ratio(:)=gyro_ratio !gyromagnetic_ratio
    call xmpi_bcast(self%gyro_ratio, master, comm, ierr)
    self%gilbert_damping(:)=damping
    call xmpi_bcast(self%gilbert_damping, master, comm, ierr)
    call self%bilinear_lil_mat%initialize(self%nspins*3,self%nspins*3)

    ABI_ALLOCATE( self%Htmp, (3, nspins))

    ! set dt default value so x/dt would not crash

  end subroutine spin_terms_t_initialize

  subroutine spin_terms_t_set_terms(self, &
       &     external_hfield, &
                                !&     exchange_i, exchange_j, exchange_val, &
                                !&     DMI_i, DMI_j, DMI_val, &
                                !&   k1,k1dir, &
       & bilinear_i, bilinear_j, bilinear_val)

    implicit none
    !Arguments ------------------------------------
    !scalars
    !arrays
    type(spin_terms_t), intent(inout) :: self

    ! Terms.
    real(dp), optional, intent(in) :: external_hfield(:,:)
    integer, optional, intent(in) :: bilinear_i(:), bilinear_j(:)
    real(dp), optional,intent(in) :: bilinear_val(:,:,:)

    !Local variables-------------------------------
    ! *************************************************************************


    if(present(external_hfield)) then
       call spin_terms_t_set_external_hfield(self, external_hfield)
    end if

    if ( present(bilinear_i) .and. present( bilinear_j) .and. present(bilinear_val) ) then
       !call spin_terms_t_set_bilinear_term(bilinear_i, bilinear_j, bilinear_val)
       call self%set_bilinear_term(bilinear_i, bilinear_j, bilinear_val)
    endif

  end subroutine spin_terms_t_set_terms

  subroutine spin_terms_t_set_external_hfield(self, external_hfield)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: external_hfield(:,:)
    ABI_ALLOCATE(self%external_hfield, (3,self%nspins))
    self%has_external_hfield = .true.
    self%external_hfield = external_hfield
  end subroutine spin_terms_t_set_external_hfield

  subroutine spin_terms_t_calc_external_Heff(self, Heff)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(out) :: Heff(:,:)
    Heff(:,:)=self%external_hfield(:,:)
  end subroutine spin_terms_t_calc_external_Heff

  subroutine spin_terms_t_set_bilinear_term_single(self, i, j, val)

    class(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val(:,:)
    integer :: ia, ib
    self%has_bilinear=.True.
    do ia = 1, 3, 1
       do ib=1, 3, 1
          call self%bilinear_lil_mat%insert(irow=(i-1)*3+ia, &
               icol=(j-1)*3+ib,val=val(ia,ib),mode=1)
       end do
    end do
    call xmpi_bcast(self%has_bilinear, master, comm, ierr)
  end subroutine spin_terms_t_set_bilinear_term_single

  subroutine spin_terms_t_set_bilinear_term(self, idx_i, idx_j, val)

    class(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: idx_i(:), idx_j(:)
    real(dp), intent(in) :: val(:,:,:)
    integer :: i, j, ia, ib, nnz
    self%has_bilinear=.True.
    nnz=size(idx_i)
    do i = 1, nnz, 1
       do ia = 1, 3, 1
          do ib=1, 3, 1
             call self%bilinear_lil_mat%insert(irow=(idx_i(i)-1)*3+ia, &
                  icol=(idx_j(i)-1)*3+ib,val=val(ia,ib,i),mode=1)
          end do
       end do
    end do
    call xmpi_bcast(self%has_bilinear, master, comm, ierr)
  end subroutine spin_terms_t_set_bilinear_term

  subroutine spin_terms_t_calc_bilinear_term_Heff(self, S, Heff)

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(inout) :: S(:,:)
    real(dp), intent(out) :: Heff(3,self%nspins)
    integer :: i, iatom, jatom
    if (.not. self%csr_mat_ready) then
       if(iam_master) then
          call LIL_to_CSR(self%bilinear_lil_mat, self%bilinear_csr_mat)
       endif
       call self%bilinear_csr_mat%sync()
       self%csr_mat_ready=.True.
    endif
    call self%bilinear_csr_mat%mv_mpi(S ,Heff)
    if(iam_master) then
       do i =1, self%nspins
          Heff(:, i)=Heff(:,i)/self%ms(i)*2.0_dp
       end do
    endif
  end subroutine spin_terms_t_calc_bilinear_term_Heff

  subroutine spin_terms_t_calculate(self, displacement, strain, spin, force, stress, bfield, energy)
    class(spin_terms_t), intent(inout) :: self  
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
    ! if present in input
    ! calculate if required
    if (present(displacement) .or. present(strain) .or. present(force) .or. present(stress)) then
       write(std_err, *) "spin terms cannot accept atomic input and output"
    end if
    if (present(bfield)) then
       call self%get_Heff(spin, bfield, energy)
    else
       call self%get_energy(spin, energy)
    end if
  end subroutine spin_terms_t_calculate


  subroutine spin_terms_t_total_Heff(self,S, Heff, energy)
    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(inout):: S(3,self%nspins)
    real(dp), intent(inout):: Heff(3,self%nspins)
    real(dp), intent(inout) :: energy
    integer :: i, j

    Heff(:,:) =0.0_dp
    energy=0.0_dp

    if(self%has_bilinear) then
       call spin_terms_t_calc_bilinear_term_Heff(self,S,self%Htmp)
       Heff=Heff+self%Htmp
    endif

    if (self%has_dipdip) then
       continue
       ! TODO implement dipdip and add it
    endif

    ! calculate energy from bilinear terms (all the above ones)
    do i=1, self%nspins
       do j=1, 3
          energy=energy-(Heff(j, i)*S(j,i)*self%ms(i))*0.5_dp
       end do
    end do

    if (self%has_external_hfield) then
       call spin_terms_t_calc_external_Heff(self,self%Htmp)
       Heff = Heff+self%Htmp
       energy= energy- self%Htmp(j, i)*S(j, i)*self%ms(i)
    endif

  end subroutine spin_terms_t_total_Heff

  subroutine spin_terms_t_get_energy(self, S, energy)
    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(inout):: S(3,self%nspins)
    real(dp), intent(inout) :: energy
    integer :: i, j
    energy=0.0_dp

    if(self%has_bilinear) then
       if (.not. self%csr_mat_ready) then
          call LIL_to_CSR(self%bilinear_lil_mat, self%bilinear_csr_mat)
          call self%bilinear_csr_mat%sync()
          self%csr_mat_ready=.True.
       endif
       call self%bilinear_csr_mat%mv_mpi(S ,self%Htmp)
       if(iam_master) then
          energy=energy - sum(sum(self%Htmp* S, dim=1))
       endif
    end if

    if(iam_master) then
       if (self%has_external_hfield) then
          call spin_terms_t_calc_external_Heff(self,self%Htmp)
          do i=1, self%nspins
             do j=1, 3
                energy= energy- self%Htmp(j, i)*S(j, i)*self%ms(i)
             end do
          end do
       endif
    end if
  end subroutine spin_terms_t_get_energy


  subroutine spin_terms_t_get_delta_E(self, S, ispin, Snew, deltaE)
    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(inout):: S(:,:), Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(out) ::deltaE
    !real(dp) ::  Eold, Enew
    real(dp) :: tmp(3), dS(3)
    ! naive implementation, just for test
    !call self%get_Heff(S, self%Htmp, Eold)
    !call self%bilinear_csr_mat%mv(S, self%Htmp)
    !Eold=-sum(sum(self%Htmp(:,:)*S(:,:), dim=1))

    !Stmp(:,:)=S(:,:)
    !Stmp(:, ispin)= Snew(:)
    !call self%get_Heff(Stmp, self%Htmp, Enew)
    !deltaE=Enew-Eold

    ! test
    deltaE=0.0_dp
    dS(:)=Snew(:)-S(:, ispin)
    S(:, ispin)= S(:, ispin)+ dS

    if(self%has_bilinear) then
       if (.not. self%csr_mat_ready) then
          call LIL_to_CSR(self%bilinear_lil_mat, self%bilinear_csr_mat)
          self%csr_mat_ready=.True.
       endif
       call self%bilinear_csr_mat%mv_select_row(3, [3*ispin-2, 3*ispin-1, 3*ispin], S, tmp)

       deltaE=deltaE-dot_product(tmp, dS ) *2.0
    end if

    if(self%has_external_hfield) then
       deltaE=deltaE - dot_product(self%external_hfield(:, ispin), dS)*self%ms(ispin)
    end if
    S(:, ispin)=S(:,ispin)-dS
  end subroutine spin_terms_t_get_delta_E



  subroutine spin_terms_t_finalize(self)

    class(spin_terms_t), intent(inout):: self

    self%is_null=.True.
    if (allocated(self%S)) then
       ABI_DEALLOCATE(self%S)
    endif

    if (allocated(self%Htmp)) then
       ABI_DEALLOCATE(self%Htmp)
    endif


    if (allocated(self%ms))  then
       ABI_DEALLOCATE(self%ms)
    endif

    if (allocated(self%pos))  then
       ABI_DEALLOCATE(self%pos)
    endif


    if (allocated(self%spinat))  then
       ABI_DEALLOCATE(self%spinat)
    endif


    if (allocated(self%iatoms)) then
       ABI_DEALLOCATE(self%iatoms)
    endif
    if (allocated(self%ispin_prim)) then
       ABI_DEALLOCATE(self%ispin_prim)
    endif
    if (allocated(self%rvec)) then
       ABI_DEALLOCATE(self%rvec)
    endif


    if (allocated(self%gyro_ratio)) then
       ABI_DEALLOCATE(self%gyro_ratio)
    endif


    if (allocated(self%external_hfield)) then 
       ABI_DEALLOCATE(self%external_hfield)
    endif


    if (allocated(self%k1))  then
       ABI_DEALLOCATE(self%k1)
    endif


    if (allocated(self%k1dir)) then
       ABI_DEALLOCATE(self%k1dir)
    endif


    if (allocated(self%gilbert_damping))  then
       ABI_DEALLOCATE(self%gilbert_damping)
    endif

    self%has_exchange=.False.
    if (allocated(self%exchange_i))  then
       ABI_DEALLOCATE(self%exchange_i)
    endif


    if (allocated(self%exchange_j))  then
       ABI_DEALLOCATE(self%exchange_j)
    endif


    if (allocated(self%exchange_val))  then
       ABI_DEALLOCATE(self%exchange_val)
    endif



    self%has_DMI=.False.
    if (allocated(self%DMI_i))  then
       ABI_DEALLOCATE(self%DMI_i)
    endif


    if (allocated(self%DMI_j)) then
       ABI_DEALLOCATE(self%DMI_j)
    endif

    if (allocated(self%DMI_val))  then
       ABI_DEALLOCATE(self%DMI_val)
    endif

    self%has_dipdip=.False.
    self%has_bilinear=.False.
    ! destroy LIL an CSR
    call self%bilinear_csr_mat%finalize()
    call self%bilinear_lil_mat%finalize()

  end subroutine spin_terms_t_finalize

end module m_spin_terms
