! Functions used in spin_terms.
! Including the calculating of effective H field.
! ds/dt= f(Heff)
! Langevin

!TODO hexu: merge this file with m_spin_terms.F90
! this file exists only for historical reasons.

module m_spin_terms_funcs
  use defs_basis
  use m_mathfuncs, only: cross, outer_product, rand_normal, rand_normal_ziggurat
  implicit none
  real(dp), parameter :: boltzmann=1.38064852d-23

CONTAINS

  ! External H field. (Too simple to be called?)
  subroutine Zeeman_Heff(N,Hext, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Zeeman_Heff'
!End of the abilint section

    integer, intent(in) :: N
    real(dp), intent(in) :: Hext(3,N)
    real(dp), intent(out) :: Heff(3,N)
    Heff(:,:)=Hext(:,:)
  end subroutine Zeeman_Heff

  ! Homogeneous uniaxial single ion anistropy (not used, to be removed?)
  subroutine homo_uniaxial_MCA_Heff(N, k1 ,k1dir , ms, S, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'homo_uniaxial_MCA_Heff'
!End of the abilint section

    integer, intent(in) :: N
    real(dp), intent(in) :: k1(N), k1dir(3,N), ms(N), S(3,N)
    real(dp), intent(out) :: Heff(3,N)
    integer :: i
    do i = 1, N, 1
       Heff(:,i)=2.0*k1(i)* (dot_product(S(:,i),k1dir(:,i))/ms(i)*k1dir(:,i))
    end do
  end subroutine homo_uniaxial_MCA_Heff

  ! Uniaxial single ion anistropy
  subroutine uniaxial_MCA_Heff(N, k1 ,k1dir , ms, S, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'uniaxial_MCA_Heff'
!End of the abilint section

    integer, intent(in) :: N
    real(dp), intent(in) :: k1(N), k1dir(3,N), ms(N), S(3,N)
    real(dp), intent(out) :: Heff(3,N)
    real(dp) :: sk
    integer :: i

    do i = 1, N, 1
       sk= dot_product(S(:,i),k1dir(:,i))
       Heff(:,i)=2.0*k1(i)* sk/ms(i)*k1dir(:,i)
    end do
  end subroutine uniaxial_MCA_Heff

  ! Exchange
  subroutine exchange_Heff(Nij, N, ilist, jlist, vallist, S, ms, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'exchange_Heff'
!End of the abilint section

    integer, intent(in):: Nij, N
    integer, intent(in)::  ilist(Nij), jlist(Nij)
    real(dp), intent(in) :: vallist(Nij), S(3,N), ms(N)
    real(dp), intent(out) :: Heff(3, N)
    integer :: i, iatom, jatom, rowS=3
    Heff(:,:)=0.0

    ! MKL version , which worked for S(N,3) , Heff(N,3).
    ! TODO Modify it so it works for S(3, N). How to do it without needing to transpose?
    !  real(dp) :: alpha, beta
    !  alpha=1.0
    !  beta=0.0
    ! ! MKL version
    ! call  mkl_dcoomm('n', N, rowS, N, alpha, 'GLNFOO', vallist, ilist, jlist, nij, S, N, beta, Heff, N)
    ! do iatom =1, N
    !      Heff(iatom,:)=Heff(iatom,:)/ms(iatom)
    ! end do

    !    No MKL version
    do i = 1, Nij, 1
       iatom=ilist(i)
       jatom=jlist(i)
       Heff(:,iatom) = Heff(:,iatom)+vallist(i)*S(:,jatom)/ms(iatom)
    end do
  end subroutine exchange_Heff

  ! DM interaction
  subroutine DMI_Heff(Nint, N, ilist, jlist, vallist, S, ms, Heff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'DMI_Heff'
!End of the abilint section

    integer, intent(in):: Nint, N, ilist(:), jlist(:)
    real(dp), intent(in) :: vallist(3, Nint), S(3,N), ms(N)
    real(dp), intent(out) :: Heff(3, N)
    integer :: i, iatom, jatom
    Heff(:,:)=0.0
    do i = 1, Nint, 1
       iatom=ilist(i)
       jatom=jlist(i)
       !write(*,*) "jatom", jatom
       !write(*,*) "iatom", iatom
       !write(*,*) "=================="
       !write(*,*) "DMI value", vallist(:,i)
       !write(*,*) "=================="
       !write(*,*) "v:", v
       !write(*,*) "=================="
       !write(*,*) "S:", S(:, jatom)
       !write(*,*) "=================="
       Heff(:,iatom) = Heff(:,iatom)+ cross(vallist(:,i),S(:,jatom))/ms(iatom)
    end do
  end subroutine DMI_Heff

  ! Langevin term heat bath
  subroutine langevin_term( n, gilbert_damping, T, gyro_ratio, ms, dt, Heff,seed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'langevin_term'
!End of the abilint section

    integer, intent(in) :: n
    real(dp), intent(in) :: gilbert_damping(:), T, gyro_ratio(:), ms(:), dt
    real(dp), intent(out):: Heff(3,n)
    real(dp) :: x(3,n), C
    integer, intent(inout):: seed
    integer :: i, j
    !print*, "dtspin: ", dt
    !print *, "damping:", gilbert_damping(1)
    !print *, "gyro:", gyro_ratio(1)
    !print *, "T", T
    if (T .gt. 1e-15_dp) then
       call rand_normal(x)
       !x(:,:)=1.0d0
       !call r8vec_normal_01(3*n, seed, x)
       do i = 1,n
          C=sqrt(2.0*gilbert_damping(1)*boltzmann*T/(gyro_ratio(1)*dt*ms(1)))
          do j = 1, 3
             Heff(j, i)= x(j,i)*C  !sqrt(2.0*gilbert_damping(i)*boltzmann*T/(gyro_ratio(i)*dt*ms(i)))
          end do
       end do
    else
       Heff(:,:)=0.0_dp
    end if
  end subroutine langevin_term

end module m_spin_terms_funcs
