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
  use m_mathfuncs
  use m_spin_terms_funcs
  use m_sparse_matrix
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
  implicit none
!!***

  ! TODO move parameters to somewhere (where?)
  real(dp), parameter :: bohr_mag=9.27400995e-24_dp, gyromagnetic_ratio = 1.76e11_dp
  type spin_terms_t
     integer :: nmatoms, natoms
     real(dp) :: etot
     ! ispin_prim: index in the spin model in primitive cell, which is used as the index of sublattice.
     ! rvec: supercell R vector.
     integer, allocatable :: ispin_prim(:), rvec(:,:)
     real(dp), allocatable :: ms(:), pos(:,:),  spinat(:,:), S(:,:)
     real(dp):: cell(3,3)
     integer, allocatable :: zion(:)
     ! index of atoms in the spin hamiltonian. -1 if the atom is not in the spin_hamiltonian.
     integer, allocatable :: spin_index(:)
     !integer :: seed
     type(rng_t) :: rng 
     ! random seed
     ! has_xxx : which terms are included in hamiltonian
     ! use_bilinear_mat_form: if false, for each ij pair,
     ! J (exchange) is a scalar, D (DMI) is a vector,
     ! which may be more efficient than a matrix. But a matrix form
     ! is more general.
     logical :: has_external_hfield, has_uniaxial_anistropy, has_exchange, &
          has_DMI, has_dipdip, has_bilinear

     !   real(dp), allocatable :: bilinear_exch(:,:,:,:)
     !     Bilinear exchange, 3rd dimension unknown until all neighbours found
     !     Size should be 3, 3, max_num_int, nmatoms (max_num_int is the number of interactions. Each spin can have a different number of neighbours.)

     !  comment hexu . removed. instead three list i, j, val is used for exchange.
     ! For exchange val is a scalar for each pair i, j. For DMI, val is a 3-vector.
     ! For dipole dipole, val is a 3*3 matrix.

     ! current in python version, i,j,R,val_ijR is stored in a class derived from self one.

     ! Array or scalar?
     real(dp), allocatable :: external_hfield(:,:)

     ! Exchange/DMI/dipdip stored like COO sparse matrix form.

     ! nint: total number of i,j pairs. size nint
     integer :: exchange_nint
     integer, allocatable:: exchange_i(:)
     integer, allocatable:: exchange_j(:)
     ! J_ij is scalar
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
     !     Size should be nmatoms

     real(dp), allocatable :: k1(:)
     !     On-site uniaxial anisotropy
     !     Size should be nmatoms

     real(dp), allocatable :: k1dir(:,:)
     !     Direction of On-site uniaxial anisotropy
     !     Size should be 3, nmatoms
     !   Should be normalized when being set.

     real(dp), allocatable :: gilbert_damping(:)
     ! Heat bath
     ! HexuComm: use gilbert_damping instead of lambda to avoid conflict with python
     ! (lambda is a python keyword)
     !     Coupling to thermal bath can be on-site
     !     Size should be nmatoms

     logical :: gamma_l_calculated = .False.
     real(dp), allocatable :: gamma_l(:)
     ! gamma_l= gyro_ratio/(1+gilbert_damping)**2


     !   integer, allocatable :: bilin_num_int(:)
     !     The number of interactions of spin i (the array value) with its neighbours
     !     Size should be nmatoms. No value should be larger than nmatoms

     ! TOMCOM - NOT SURE IF THE TOTAL NUMBER OF INTERACTIONS SHOULD BE ifc_type
     !   type(ifc_type) :: max_num_int
     !     total number of interactions
     ! TOMCOM - UNTIL HERE

     ! vector for calculating effective field
     real(dp), allocatable :: Htmp(:)
  !  CONTAINS
  !    procedure :: initialize => spin_terms_t_initialize
  !    procedure :: finalize => spin_terms_t_finalize
  !    procedure :: get_Heff => total_Heff
  !    procedure :: Heff_to_dSdt => Heff_to_dSdt
  !    procedure :: get_dSdt => get_dSdt
  !    procedure :: get_Langevin_Heff => get_Langevin_Heff
   end type spin_terms_t
contains
  subroutine spin_terms_t_initialize(self, cell, pos, spinat, zion, spin_index, ispin_prim, rvec)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_initialize'
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    !scalars
    class(spin_terms_t), intent(out) :: self
    integer, intent(in) :: ispin_prim(:), rvec(:,:)
    integer, intent(in) :: zion(:), spin_index(:)
    real(dp), intent(in) :: cell(3,3), pos(:,:), spinat(:,:)
    !Local variables-------------------------------
    integer :: nmatoms, natoms, i
    ! TODO: should be sth else, not zion
    nmatoms=size(zion)
    natoms=size(spin_index)
    self%nmatoms=nmatoms
    self%natoms=natoms
    ABI_ALLOCATE(self%zion, (nmatoms))
    ABI_ALLOCATE(self%spin_index, (natoms))
    self%zion(:)=zion(:)
    self%spin_index(:)=spin_index(:)
    self%cell(:,:)=cell(:,:)
    ABI_ALLOCATE( self%pos, (3,nmatoms) )
    self%pos(:,:)=pos(:,:)
    ABI_ALLOCATE( self%spinat, (3, nmatoms))
    self%spinat(:,:)=spinat(:,:)

    ABI_ALLOCATE( self%ms, (nmatoms) )
    self%ms(:)= sqrt(sum(spinat(:,:)**2, dim=1))*bohr_mag

    ABI_ALLOCATE( self%ispin_prim, (nmatoms))
    ABI_ALLOCATE(self%rvec, (3, nmatoms))

    self%ispin_prim(:)=ispin_prim(:)
    self%rvec(:,:)=rvec(:,:)

    ABI_ALLOCATE( self%S, (3, nmatoms))
    do i=1,nmatoms
       self%S(:,i)=self%spinat(:,i)/self%ms(i)*bohr_mag
    end do

    self%has_external_hfield=.False.
    self%has_uniaxial_anistropy=.False.
    self%has_exchange=.False.
    self%has_DMI=.False.
    self%has_dipdip=.False.
    self%has_bilinear=.False.

    ABI_ALLOCATE( self%gyro_ratio, (nmatoms))
    ! Defautl gyro_ratio
    ! TODO remove this, use gyroration from xml
    self%gyro_ratio(:)=gyromagnetic_ratio

    ABI_ALLOCATE( self%gilbert_damping, (nmatoms) )
    ABI_ALLOCATE( self%gamma_l, (nmatoms))

    call LIL_mat_initialize(self%bilinear_lil_mat,self%nmatoms*3,self%nmatoms*3)
    ABI_ALLOCATE( self%Htmp, (nmatoms*3))

    call set_seed(self%rng, [111111_dp, 2_dp])
  end subroutine spin_terms_t_initialize

  subroutine spin_terms_t_initialize_all(self,nmatoms, ms, &
       &     external_hfield, &
       &     exchange_i, exchange_j, exchange_val, &
       &     DMI_i, DMI_j, DMI_val, &
       &    gyro_ratio,k1,k1dir,gilbert_damping, &
       & bilinear_i, bilinear_j, bilinear_val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_initialize_all'
!End of the abilint section

    implicit none
    !Arguments ------------------------------------
    !scalars
    integer, intent(in) :: nmatoms
    !arrays
    type(spin_terms_t), intent(out) :: self
    real(dp), optional, intent(in) :: external_hfield(:,:)
    !integer, optional,intent(in) :: exchange_nint, DMI_nint
    integer, optional,intent(in) :: exchange_i(:), exchange_j(:), &
         & dmi_i(:), dmi_j(:), bilinear_i(:), bilinear_j(:)
    real(dp), intent(in) :: ms(:)
    real(dp), optional,intent(in) :: exchange_val(:), dmi_val(:,:), bilinear_val(:,:,:)
    real(dp),optional,intent(in) :: gyro_ratio(nmatoms), gilbert_damping(nmatoms)
    real(dp),optional,intent(in) :: k1(nmatoms)
    real(dp),optional,intent(in) :: k1dir(3,nmatoms)
    !Local variables-------------------------------

    ! *************************************************************************

    ABI_ALLOCATE( self%ms, (nmatoms))
    self%ms = ms

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

    if ( present(gyro_ratio) ) then
       ABI_ALLOCATE(self%gyro_ratio, (nmatoms))
    endif
    self%gyro_ratio=gyro_ratio

    if ( present(k1) .and. present( k1dir) ) then
       call spin_terms_t_set_uniaxial_MCA(self, k1, k1dir)
    endif

    if (present(gilbert_damping)) then
       ABI_ALLOCATE(self%gilbert_damping, (self%nmatoms))
       self%gilbert_damping=gilbert_damping
    endif

    call set_seed(self%rng, [111111_dp, 2_dp])
  end subroutine spin_terms_t_initialize_all

  subroutine spin_terms_t_get_gamma_l(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_get_gamma_l'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    self%gamma_l(:)= self%gyro_ratio(:)/(1.0_dp+ self%gilbert_damping(:)**2)
    self%gamma_l_calculated=.True.
  end subroutine spin_terms_t_get_gamma_l

  subroutine spin_terms_t_set_external_hfield(self, external_hfield)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_external_hfield'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: external_hfield(:,:)
    ABI_ALLOCATE(self%external_hfield, (3,self%nmatoms))
    self%has_external_hfield = .true.
    self%external_hfield = external_hfield
  end subroutine spin_terms_t_set_external_hfield

  subroutine spin_terms_t_calc_external_Heff(self, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_calc_external_Heff'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(out) :: Heff(:,:)
    Heff(:,:)=self%external_hfield(:,:)
  end subroutine spin_terms_t_calc_external_Heff

  subroutine spin_terms_t_set_uniaxial_MCA(self, k1, k1dir)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_uniaxial_MCA'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: k1(:), k1dir(:,:)
    integer :: i
    real(dp) :: norm

    ABI_ALLOCATE(self%k1, (self%nmatoms))
    ABI_ALLOCATE(self%k1dir, (3,self%nmatoms))
    self%has_uniaxial_anistropy= .true.
    self%k1=k1
    do i = 1, self%nmatoms
       norm = sqrt(sum( k1dir(:,i)**2))
       self%k1dir(:,i)=k1dir(:,i)/norm
    end do
  end subroutine spin_terms_t_set_uniaxial_MCA

  subroutine spin_terms_t_calc_uniaxial_MCA_Heff(self, S, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_calc_uniaxial_MCA_Heff'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    call uniaxial_MCA_Heff(self%nmatoms,self%k1,self%k1dir,self%ms,S,Heff)
  end subroutine spin_terms_t_calc_uniaxial_MCA_Heff

  subroutine spin_terms_t_set_bilinear_term_single(self, i, j, val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_bilinear_term_single'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_bilinear_term'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_calc_bilinear_term_Heff'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(3,self%nmatoms)
    integer :: i, iatom, jatom
    real(dp) ::  H(3,1)
    !Svec=reshape(S, (/self%nmatoms*3/))
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
    !call CSR_mat_mv(self%bilinear_csr_mat, reshape(S, [3*self%nmatoms]),self%Htmp)
    !Heff(:,:) = reshape (self%Htmp, [3, self%nmatoms])
    call CSR_mat_mv(self%bilinear_csr_mat, S ,Heff)
    !print *, "ms", self%ms
    !$OMP PARALLEL DO
    do i =1, self%nmatoms
       Heff(:, i)=Heff(:,i)/self%ms(i)
    end do
    !$OMP END PARALLEL DO
    !print *, "Heff", Heff
  end subroutine spin_terms_t_calc_bilinear_term_Heff

  subroutine spin_terms_t_set_exchange(self, exchange_i, exchange_j, exchange_val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_exchange'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_calc_exchange_Heff'
!End of the abilint section

    type(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    call exchange_Heff(self%exchange_nint,self%nmatoms,self%exchange_i, &
         self%exchange_j,self%exchange_val,S,self%ms,Heff)
  end subroutine spin_terms_t_calc_exchange_Heff

  subroutine spin_terms_t_set_DMI(self, DMI_i, DMI_j, DMI_val)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_set_DMI'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_calc_DMI_Heff'
!End of the abilint section

    type(spin_terms_t), intent(in) :: self
    real(dp), intent(in) :: S(:,:)
    real(dp), intent(out) :: Heff(:,:)
    integer :: nij
    nij=size(self%DMI_i)
    call DMI_Heff(nij,self%nmatoms,self%DMI_i,self%DMI_j,self%DMI_val,S,self%ms,Heff)
  end subroutine spin_terms_t_calc_DMI_Heff


  subroutine spin_terms_t_total_Heff(self,S, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_total_Heff'
!End of the abilint section

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(in):: S(3,self%nmatoms)
    real(dp), intent(out):: Heff(3,self%nmatoms)
    real(dp) :: Htmp(3,self%nmatoms)
    Heff(:,:) =0.0_dp

    if(.not. self%gamma_l_calculated) then
       call spin_terms_t_get_gamma_l(self)
    end if

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
       !write (*,*) "Htmp", Htmp
       !write (*,*) "Heff", Heff
    endif

    if (self%has_dipdip) then
       continue
       ! TODO implement dipdip and add it
    endif

    if (self%has_external_hfield) then
       call spin_terms_t_calc_external_Heff(self,Htmp)
       Heff = Heff+Htmp
    endif

    if ( self%has_uniaxial_anistropy ) then
       call spin_terms_t_calc_uniaxial_MCA_Heff(self,S,Htmp)
       Heff = Heff+Htmp
    end if

  end subroutine spin_terms_t_total_Heff
  
  ! A effective torque from Langevin heat bath
  subroutine spin_terms_t_get_Langevin_Heff(self, dt, temperature, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_get_Langevin_Heff'
!End of the abilint section

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: dt, temperature
    real(dp), intent(out):: Heff(3,self%nmatoms)
    real(dp) :: x(3, self%nmatoms), C
    integer :: i, j
    if ( temperature .gt. 1d-7) then
       !call rand_normal_ziggurat(x)
       call rand_normal_array(self%rng, x, 3*self%nmatoms)
       do i = 1, self%nmatoms
          C=sqrt(2.0*self%gilbert_damping(i)*boltzmann* temperature &
               &  /(self%gyro_ratio(i)* dt *self%ms(i)))
          do j = 1, 3
             Heff(j, i)= x(j,i)*C
          end do
       end do
    else
       Heff(:,:)=0.0_dp
    end if
  end subroutine spin_terms_t_get_Langevin_Heff

  ! ds/dt = f(Heff, S)
  subroutine spin_terms_t_Heff_to_dsdt(self, Heff, S, dSdt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_Heff_to_dsdt'
!End of the abilint section

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: Heff(3,self%nmatoms), S(3,self%nmatoms)
    real(dp), intent(out) :: dSdt(3, self%nmatoms)
    integer :: i
    real(dp) :: Ri(3)
    !!$OMP PARALLEL DO private(Ri)
    do i=1,self%nmatoms
       Ri = cross(S(:,i),Heff(:,i))
       dSdt(:,i) = -self%gamma_L(i)*(Ri+self%gilbert_damping(i)* cross(S(:,i), Ri))
    end do
    !!$OMP END PARALLEL DO
  end subroutine spin_terms_t_Heff_to_dsdt

  subroutine spin_terms_t_get_etot(self, S, Heff, etot)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_get_etot'
!End of the abilint section

    class(spin_terms_t), intent(inout) :: self
    real(dp), intent(in) :: S(3, self%nmatoms), Heff(3, self%nmatoms)
    real(dp), intent(out) :: etot
    integer :: i, j
    etot=0.0_dp
    do i=1, self%nmatoms
       do j=1, 3
          etot=etot-Heff(j, i)*S(j,i)*self%ms(i)
       end do
    end do
  end subroutine spin_terms_t_get_etot

  subroutine spin_terms_t_get_dSdt(self,  S, H_lang, dSdt )


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_get_dSdt'
!End of the abilint section

     class(spin_terms_t) ,intent(inout) :: self
     real(dp), intent(in):: S(3, self%nmatoms), H_lang(3, self%nmatoms)
     real(dp), intent(out):: dSdt(3, self%nmatoms)
     real(dp):: Heff(3, self%nmatoms)
     !call self%get_Heff(S=S, Heff=Heff)
     call spin_terms_t_total_Heff(self=self, S=S, Heff=Heff)
     call spin_terms_t_get_etot(self=self, S=S, Heff=Heff, etot=self%etot)
     Heff(:,:)=Heff(:,:)+H_lang(:,:)
     !call self%Heff_to_dsdt(Heff,S, dSdt)
     call spin_terms_t_Heff_to_dsdt(self, Heff, S, dSdt)

   end subroutine spin_terms_t_get_dSdt

   subroutine spin_terms_t_finalize(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spin_terms_t_finalize'
!End of the abilint section

    class(spin_terms_t), intent(inout):: self
    integer::  err
    if (allocated(self%ms))  then
    ABI_DEALLOCATE(self%ms)
    endif

    if (allocated(self%pos))  then
    ABI_DEALLOCATE(self%pos)
    endif


    if (allocated(self%spinat))  then
     ABI_DEALLOCATE(self%spinat)
    endif

   
    if (allocated(self%zion))  then
    ABI_DEALLOCATE(self%zion)
    endif


    if (allocated(self%spin_index)) then
       ABI_DEALLOCATE(self%spin_index)
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


    if (allocated(self%gamma_l))  then
    ABI_DEALLOCATE(self%gamma_l)
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
