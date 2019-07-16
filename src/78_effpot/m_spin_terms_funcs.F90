!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spin_terms_funcs
!! NAME
!! m_spin_terms_funcs
!!
!! FUNCTION
!! This module contains the subroutines to calcuate Heff from spin Hamiltonian terms.
!!
!!
!! Datatypes:
!!
!!
!! Subroutines:
!!
!!
!!TODO hexu: merge this file with m_spin_terms.F90
!! this file exists only for historical reasons.
!!
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

module m_spin_terms_funcs
  use defs_basis
   use m_errors
    use m_abicore
  use m_mathfuncs, only: cross, outer_product 
  implicit none
  real(dp), parameter :: boltzmann=1.38064852d-23 ! TODO where is it in abinit.

CONTAINS

  !!***

  ! External H field. (Too simple to be called?)
  subroutine Zeeman_Heff(N,Hext, Heff)

    integer, intent(in) :: N
    real(dp), intent(in) :: Hext(3,N)
    real(dp), intent(out) :: Heff(3,N)
    Heff(:,:)=Hext(:,:)
  end subroutine Zeeman_Heff

  ! Homogeneous uniaxial single ion anistropy (not used, to be removed?)
  subroutine homo_uniaxial_MCA_Heff(N, k1 ,k1dir , ms, S, Heff)

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

    integer, intent(in):: Nij, N
    integer, intent(in)::  ilist(Nij), jlist(Nij)
    real(dp), intent(in) :: vallist(Nij), S(3,N), ms(N)
    real(dp), intent(out) :: Heff(3, N)
    integer :: i, iatom, jatom !, rowS=3
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

    integer, intent(in):: Nint, N, ilist(:), jlist(:)
    real(dp), intent(in) :: vallist(3, Nint), S(3,N), ms(N)
    real(dp), intent(out) :: Heff(3, N)
    integer :: i, iatom, jatom
    Heff(:,:)=0.0
    do i = 1, Nint, 1
       iatom=ilist(i)
       jatom=jlist(i)
       Heff(:,iatom) = Heff(:,iatom)+ cross(vallist(:,i),S(:,jatom))/ms(iatom)
    end do
  end subroutine DMI_Heff

end module m_spin_terms_funcs
