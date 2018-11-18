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
  use m_spin_terms_funcs
  use m_sparse_matrix
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
     type(LIL_mat) :: bilinear_lil_mat
     logical :: csr_mat_ready= .False.
     type(CSR_mat) :: bilinear_csr_mat
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
     real(dp), allocatable :: Htmp(:)
       CONTAINS
         procedure :: initialize => spin_terms_t_initialize
         procedure :: finalize => spin_terms_t_finalize
         procedure :: get_Heff => spin_terms_t_total_Heff
         procedure :: calculate => spin_terms_t_calculate
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
    ABI_ALLOCATE(self%iatoms, (nspins))
    self%iatoms(:)=iatoms(:)
    self%cell(:,:)=cell(:,:)
    ABI_ALLOCATE( self%pos, (3,nspins) )
    self%pos(:,:)=pos(:,:)
    ABI_ALLOCATE( self%spinat, (3, nspins))
    self%spinat(:,:)=spinat(:,:)

    ABI_ALLOCATE( self%ms, (nspins) )
    self%ms(:)= sqrt(sum(spinat(:,:)**2, dim=1))* mu_B

    ABI_ALLOCATE( self%ispin_prim, (nspins))
    ABI_ALLOCATE(self%rvec, (3, nspins))

    self%ispin_prim(:)=ispin_prim(:)
    self%rvec(:,:)=rvec(:,:)

    ABI_ALLOCATE( self%S, (3, nspins))
    do i=1,nspins
       self%S(:,i)=self%spinat(:,i)/self%ms(i)
    end do

    self%has_external_hfield=.False.
    self%has_uniaxial_anistropy=.False.
    self%has_exchange=.False.
    self%has_DMI=.False.
    self%has_dipdip=.False.
    self%has_bilinear=.False.

    ABI_ALLOCATE( self%gyro_ratio, (nspins))
    ABI_ALLOCATE( self%gilbert_damping, (nspins) )
    ! Defautl gyro_ratio
    self%gyro_ratio(:)=gyro_ratio !gyromagnetic_ratio
    self%gilbert_damping(:)=damping

    call LIL_mat_initialize(self%bilinear_lil_mat,self%nspins*3,self%nspins*3)
    ABI_ALLOCATE( self%Htmp, (nspins*3))

    ! set dt default value so x/dt would not crash
    
  end subroutine spin_terms_t_initialize

  subroutine spin_terms_t_set_terms(self, &
       &     external_hfield, &
       &     exchange_i, exchange_j, exchange_val, &
       &     DMI_i, DMI_j, DMI_val, &
       &   k1,k1dir, &
       & bilinear_i, bilinear_j, bilinear_val)

    implicit none
    !Arguments ------------------------------------
    !scalars
    !arrays
    type(spin_terms_t), intent(inout) :: self

    ! Terms. 
    real(dp), optional, intent(in) :: external_hfield(:,:)
    integer, optional,intent(in) :: exchange_i(:), exchange_j(:), &
         & dmi_i(:), dmi_j(:), bilinear_i(:), bilinear_j(:)
    real(dp), optional,intent(in) :: exchange_val(:), dmi_val(:,:), bilinear_val(:,:,:)
    real(dp),optional,intent(in) :: k1(self%nspins)
    real(dp),optional,intent(in) :: k1dir(3,self%nspins)

    !Local variables-------------------------------

    ! *************************************************************************


    if(present(external_hfield)) then
       call spin_terms_t_set_external_hfield(self, external_hfield)
    end if

    if (present(bilinear_i) .and. present(bilinear_j) &
         & .and. present(bilinear_val)) then
       call spin_terms_t_set_bilinear_term(self, bilinear_i, bilinear_j, bilinear_val)
    endif

    if (present(exchange_i) .and. present(exchange_j) .and. present(exchange_val) ) then
       call spin_terms_t_set_exchange(self, exchange_i, exchange_j, exchange_val)
    end if


    if ( present(DMI_i) .and. present(DMI_j) .and. present(DMI_val) ) then
       call spin_terms_t_set_DMI(self, DMI_i, DMI_j, DMI_val)
    end if

    if ( present(k1) .and. present( k1dir) ) then
       call spin_terms_t_set_uniaxial_MCA(self, k1, k1dir)
    endif

  end subroutine spin_terms_t_set_terms



  subroutine spin_terms_t_add_SIA(self, mode, k1, k1dir)

    class(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: mode
    real(dp), intent(in) :: k1, k1dir(3)
    real(dp) :: k1_tmp(self%nspins), k1dir_tmp(3, self%nspins)
    integer :: i
    ! Add
    if(mode==1 .or. mode==2) then
       ! Note: mode 1: add
       ! mode 2, override. If spin_sia_add==2, then the sia is not set from xml file.
       k1_tmp(:)=k1
       do i=1, self%nspins
          k1dir_tmp(:,i)=k1dir(:)
       end do
       call spin_terms_t_set_uniaxial_MCA(self, k1_tmp, k1dir_tmp)
    endif
  end subroutine spin_terms_t_add_SIA

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

  subroutine spin_terms_t_set_uniaxial_MCA(self, k1, k1dir)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: k1(:), k1dir(:,:)
    integer :: i
    real(dp) :: norm

    ABI_ALLOCATE(self%k1, (self%nspins))
    ABI_ALLOCATE(self%k1dir, (3,self%nspins))
    self%has_uniaxial_anistropy= .true.
    self%k1=k1
    do i = 1, self%nspins
       norm = sqrt(sum( k1dir(:,i)**2))
       self%k1dir(:,i)=k1dir(:,i)/norm
    end do
  end subroutine spin_terms_t_set_uniaxial_MCA

  subroutine spin_terms_t_calc_uniaxial_MCA_Heff(self, S, Heff)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    call uniaxial_MCA_Heff(self%nspins,self%k1,self%k1dir,self%ms,S,Heff)
  end subroutine spin_terms_t_calc_uniaxial_MCA_Heff

  subroutine spin_terms_t_set_bilinear_term_single(self, i, j, val)

    type(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: val(:,:)
    integer :: ia, ib
    self%has_bilinear=.True.
    do ia = 1, 3, 1
       do ib=1, 3, 1
          call LIL_mat_insert(self%bilinear_lil_mat,irow=(i-1)*3+ia, &
               icol=(j-1)*3+ib,val=val(ia,ib),mode=1)
       end do
    end do
  end subroutine spin_terms_t_set_bilinear_term_single

  subroutine spin_terms_t_set_bilinear_term(self, idx_i, idx_j, val)

    type(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: idx_i(:), idx_j(:)
    real(dp), intent(in) :: val(:,:,:)
    integer :: i, j, ia, ib, nnz
    !self%bilinear_nint = size(idx_i)
    !ABI_ALLOCATE(self%bilinear_i, (self%bilinear_nint))
    !ABI_ALLOCATE(self%bilinear_j, (self%bilinear_nint))
    !ABI_ALLOCATE(self%bilinear_val, (3,3,self%bilinear_nint))
    self%has_bilinear=.True.
    nnz=size(idx_i)
    do i = 1, nnz, 1
       do ia = 1, 3, 1
          do ib=1, 3, 1
             call LIL_mat_insert(self%bilinear_lil_mat,irow=(idx_i(i)-1)*3+ia, &
                  icol=(idx_j(i)-1)*3+ib,val=val(ia,ib,i),mode=1)
          end do
       end do
    end do
  end subroutine spin_terms_t_set_bilinear_term

  subroutine spin_terms_t_calc_bilinear_term_Heff(self, S, Heff)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(3,self%nspins)
    integer :: i, iatom, jatom
    !real(dp) ::  H(3,1)
    !Svec=reshape(S, (/self%nspins*3/))
    !do i = 1, self%bilinear_nint, 1
    !    iatom=self%bilinear_i(i)
    !    jatom=self%bilinear_j(i)
    ! Heff(:,i)=Heff(:,i) + matmul(S(:,i))
    !call dgemv('n', 3, 3, 1.0_dp/self%ms(iatom), self%bilinear_val(:,:,i), 3, S(:,iatom), 1, 1.0_dp , H, 1 )
    !  H=matmul(self%bilinear_val(:,:,i), reshape(S(:,iatom), (/3,1/)))
    !  Heff(:,i)=Heff(:,i)+ H(:,1)
    !end do
    !print *, "calculate bilinear Heff"
    !print *, "Heff: ", Heff
    if (.not. self%csr_mat_ready) then
       call LIL_mat_to_CSR(self%bilinear_lil_mat, self%bilinear_csr_mat)
       !call CSR_mat_print(self%bilinear_csr_mat)
       self%csr_mat_ready=.True.
    endif
    !call CSR_mat_mv(self%bilinear_csr_mat, reshape(S, [3*self%nspins]),self%Htmp)
    !Heff(:,:) = reshape (self%Htmp, [3, self%nspins])
    call CSR_mat_mv(self%bilinear_csr_mat, S ,Heff)
    !print *, "ms", self%ms
    !$OMP PARALLEL DO private(i)
    do i =1, self%nspins
       Heff(:, i)=Heff(:,i)/self%ms(i)*2.0_dp
    end do
    !$OMP END PARALLEL DO
  end subroutine spin_terms_t_calc_bilinear_term_Heff

  subroutine spin_terms_t_set_exchange(self, exchange_i, exchange_j, exchange_val)

    type(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: exchange_i(:), exchange_j(:)
    real(dp), intent(in) :: exchange_val(:)

    self%has_exchange=.true.
    self%exchange_nint=size(exchange_i)

    ABI_ALLOCATE(self%exchange_i, (self%exchange_nint))
    ABI_ALLOCATE(self%exchange_j, (self%exchange_nint))
    ABI_ALLOCATE(self%exchange_val, (self%exchange_nint))

    self%exchange_i(:)=exchange_i(:)
    self%exchange_j(:)=exchange_j(:)
    self%exchange_val(:)=exchange_val(:)
  end subroutine spin_terms_t_set_exchange

  subroutine spin_terms_t_calc_exchange_Heff(self, S, Heff)

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    call exchange_Heff(self%exchange_nint,self%nspins,self%exchange_i, &
         self%exchange_j,self%exchange_val,S,self%ms,Heff)
  end subroutine spin_terms_t_calc_exchange_Heff

  subroutine spin_terms_t_set_DMI(self, DMI_i, DMI_j, DMI_val)

    type(spin_terms_t), intent(inout) :: self
    integer, intent(in) :: DMI_i(:), DMI_j(:)
    real(dp), intent(in) :: DMI_val(:,:)

    self%has_DMI=.true.
    self%DMI_nint=size(DMI_i)
    ABI_ALLOCATE(self%DMI_i, (self%DMI_nint))
    ABI_ALLOCATE(self%DMI_j, (self%DMI_nint))
    ABI_ALLOCATE(self%DMI_val, (3,self%DMI_nint))
    self%DMI_i=DMI_i
    self%DMI_j=DMI_j
    self%DMI_val(:,:)=DMI_val(:,:)
  end subroutine spin_terms_t_set_DMI

  subroutine spin_terms_t_calc_DMI_Heff(self, S, Heff)

    type(spin_terms_t), intent(in) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    integer :: nij
    nij=size(self%DMI_i)
    call DMI_Heff(nij,self%nspins,self%DMI_i,self%DMI_j,self%DMI_val,S,self%ms,Heff)
  end subroutine spin_terms_t_calc_DMI_Heff


  subroutine spin_terms_t_calculate(self, displacement, strain, spin, force, stress, bfield, energy)
    class(spin_terms_t), intent(inout) :: self  
    real(dp), optional, intent(in) :: displacement(:,:), strain(:,:), spin(:,:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), energy
    ! if present in input
    ! calculate if required
    if (present(displacement) .or. present(strain) .or. present(force) .or. present(stress)) then
       write(std_err, *) "spin terms cannot accept atomic input and output"
    end if
    call self%get_Heff(spin, bfield, energy)
  end subroutine spin_terms_t_calculate



  subroutine spin_terms_t_total_Heff(self,S, Heff, energy)

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(in):: S(3,self%nspins)
    real(dp), intent(inout):: Heff(3,self%nspins)
    real(dp), intent(inout) :: energy
    real(dp) :: Htmp(3,self%nspins)
    integer :: i, j

    Heff(:,:) =0.0_dp
    energy=0.0_dp

    if(self%has_bilinear) then
       call spin_terms_t_calc_bilinear_term_Heff(self,S,Htmp)
       Heff=Heff+Htmp
    endif

    if ( self%has_exchange) then
       call spin_terms_t_calc_exchange_Heff(self,S,Htmp)
       Heff = Heff + Htmp
    end if

    if (self%has_DMI) then
       call spin_terms_t_calc_DMI_Heff(self,S,Htmp)
       Heff = Heff + Htmp
    endif

    if ( self%has_uniaxial_anistropy ) then
       call spin_terms_t_calc_uniaxial_MCA_Heff(self,S,Htmp)
       Heff = Heff+Htmp
    end if

    if (self%has_dipdip) then
       continue
       ! TODO implement dipdip and add it
    endif

    ! calculate energy from bilinear terms (all the above ones)
    do i=1, self%nspins
       do j=1, 3
          energy=energy-Heff(j, i)*S(j,i)*self%ms(i)*0.5_dp
       end do
    end do

    if (self%has_external_hfield) then
       call spin_terms_t_calc_external_Heff(self,Htmp)
       Heff = Heff+Htmp
       energy= energy- Htmp(j, i)*S(j, i)*self%ms(i)
    endif

  end subroutine spin_terms_t_total_Heff


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
    !if (allocated(self%dipdip_i))  then
    !ABI_DEALLOCATE(self%dipdip_i)
    !endif


    !if (allocated(self%dipdip_j))  then
    !ABI_DEALLOCATE(self%dipdip_j)
    !endif


    !if (allocated(self%dipdip_val))  then
    !ABI_DEALLOCATE(self%dipdip_val)
    !endif



    !if (allocated(self%bilinear_i)) deallocate(self%bilinear_i, stat=err)
    !if (allocated(self%bilinear_j)) deallocate(self%bilinear_j, stat=err)
    !if (allocated(self%bilinear_val)) deallocate(self%bilinear_val, stat=err)
    self%has_bilinear=.False.
    ! destroy LIL an CSR
    call CSR_mat_finalize(self%bilinear_csr_mat)
    call LIL_mat_finalize(self%bilinear_lil_mat)

  end subroutine spin_terms_t_finalize

end module m_spin_terms
