!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_supercell_maker
!! NAME
!! m_supercell_maker
!!
!! FUNCTION
!! This module define the supercell_maker file, which provide functions to help build
!! potentials in supercell.
!!
!! Datatypes:
!!  supercell_maker_t
!!
!! Subroutines:
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
#include "abi_common.h"

module m_supercell_maker
  use defs_basis
  use m_abicore
  use m_errors
  use m_symtk,    only : matr3inv
  use m_mathfuncs , only: mat33det

  implicit none
  private

  type ,public :: supercell_maker_t
     integer :: scmat(3,3)
     real(dp) :: inv_scmat(3,3)
     integer :: ncells
     integer, allocatable :: rvecs(:,:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: to_red_sc
     procedure :: build_rvec
     procedure :: sc_cell
     procedure :: R_to_sc

     ! Translations: from one primitive cell element to ncell supercell elements
     ! with values unchanged (therefore just repeat ncells times)
     generic :: repeat => repeat_int1d, repeat_real1d
     procedure :: repeat_int1d
     procedure :: repeat_real1d

     ! Translations: from one primitive cell element to ncell supercell elements
     ! with values changed (different from repeat)
     procedure :: trans_xred
     procedure :: trans_xcart
     procedure :: trans_ind
     procedure :: trans_Rj_and_jlist
     procedure :: trans_jR
  end type supercell_maker_t

  public :: scmaker_unittest
contains

  subroutine initialize(self, sc_matrix)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: sc_matrix(3, 3)
    real(dp) :: tmp(3,3)
    self%scmat(:,:)=sc_matrix
    self%ncells=mat33det(self%scmat)
    ABI_ALLOCATE(self%rvecs, (3, self%ncells))
    ! for inv
    tmp(:,:)=self%scmat(:,:)
    call matr3inv(transpose(tmp), self%inv_scmat)
    call self%build_rvec()
  end subroutine initialize

  subroutine finalize(self)
    class(supercell_maker_t), intent(inout) :: self
    self%scmat=0
    self%ncells=0
    if (allocated(self%rvecs)) then
       ABI_DEALLOCATE(self%rvecs)
    end if
  end subroutine finalize


  function to_red_sc(self, pos) result(ret)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: pos(3)
    real(dp) :: ret(3)
    ret(:)=matmul(self%inv_scmat, pos)
  end function to_red_sc

  subroutine build_rvec(self)
    class(supercell_maker_t), intent(inout) :: self
    real(dp):: scorners_newcell(8,3), corners(8,3), x(3), tmp(3)
    integer :: rep(3), ix, iy, iz, counter=0
    real(dp) :: eps=1d-8
    scorners_newcell=reshape([0., 0., 0., 0., 0., 1., 0., 1., 0., &
         0., 1., 1., 1., 0., 0., 1., 0., 1.,  &
         1., 1., 0., 1., 1., 1.], [8,3] )
    corners = matmul(scorners_newcell, self%scmat)
    rep=ceiling(maxval(corners, dim=1) - minval(corners, dim=1))
    do ix = 1, rep(1)
       do iy = 1, rep(2)
          do iz = 1, rep(3)
             x(:)=[ix, iy, iz]
             tmp=self%to_red_sc(x)
             if ( (.not. any(tmp<=-1.0*eps)) &
                  .and. (.not. any(tmp>1.0-eps)) ) then
                counter=counter+1
                self%rvecs(:,counter)=nint(tmp)
             end if
          end do
       end do
    end do
    if (counter /= self%ncells ) then
       print *, "Wrong number of supercell found"
    end if
  end subroutine build_rvec

  ! cell parameters from unitcell to supercell
  function sc_cell(self, cell) result(sccell)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: cell(3,3)
    real(dp) :: sccell(3,3)
    sccell=matmul(self%scmat, cell)
  end function sc_cell

  ! xred in list of xred in supercell by translation
  subroutine trans_xred(self, xred, scxred)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: xred(:,:)
    real(dp), intent(inout) :: scxred(3,size(xred, dim=2)*self%ncells)
    integer :: npos, icell, ipos, counter=0
    npos=size(xred, dim=2)
    do icell = 1, self%ncells
       do ipos=1 , npos
          counter=counter+1
          scxred(:, counter) = self%to_red_sc(xred(:, ipos)+ self%rvecs(:, icell))
       end do
    end do
  end subroutine trans_xred

  subroutine trans_xcart(self, primcell, xcart, scxcart)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: xcart(:,:), primcell(3,3)
    real(dp), intent(inout) :: scxcart(3,size(xcart, dim=2)*self%ncells)
    integer :: npos, icell, ipos, counter=0
    do icell = 1, self%ncells
       do ipos=1 , npos
          counter=counter+1
          scxcart(:, counter) = matmul(primcell, self%rvecs(:, icell)) + xcart(:,ipos)
       end do
    end do
  end subroutine trans_xcart

  ! R: R index using primitive cell parameter
  ! R_sc: R index using supercell parameter
  ! ind_sc: index of cell INSIDE supercell in primitive cell.
  subroutine R_to_sc(self, R, R_sc, ind_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: R(3)
    integer, intent(inout) :: ind_sc, R_sc(3)
    integer :: rprim(3)
    integer :: i
    R_sc=floor(self%to_red_sc(R*1.0d0))
    rprim(:)= R-matmul(R_sc, self%scmat)
    do i =0, self%ncells
       if (all(self%rvecs(:,i) == rprim)) then
          ind_sc=i
          exit
       end if
    end do
  end subroutine R_to_sc

  ! ind: index in primitive cell
  ! ind_sc: indices (PLURAL) in supercell
  subroutine trans_ind(self, nbasis, ind, ind_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: ind, nbasis
    integer, intent(out) :: ind_sc(self%ncells)
    integer :: i
    do i =1, self%ncells
       ind_sc(i)=nbasis*(i-1)+ind
    end do
  end subroutine trans_ind

  subroutine trans_jR(self, nbasis, j, Rj, ind_sc, R_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: j, Rj(3), nbasis
    integer, intent(out) :: ind_sc(self%ncells), R_sc(3, self%ncells)
    integer :: i,jj
    do i =1, self%ncells
       call self%R_to_sc(Rj + self%rvecs(:,i), R_sc(:,i), jj)
       ind_sc(i)=nbasis*(jj-1)+j
    end do
  end subroutine trans_jR


  subroutine trans_Rj_and_jlist(self, nbasis, jlist, Rj, ind_sc, R_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: jlist(:), Rj(3), nbasis
    integer, intent(out) :: ind_sc(self%ncells*size(jlist)), R_sc(3, self%ncells)
    integer :: i,jj, counter, indj
    counter=0
    do i =1, self%ncells
       call self%R_to_sc(Rj + self%rvecs(:,i), R_sc(:,i), jj)
       do indj=1, size(jlist)
          counter=counter+1
          ind_sc(counter)=nbasis*(jj-1)+jlist(indj)
       end do
    end do
  end subroutine trans_Rj_and_jlist


  subroutine repeat_int1d(self, a, ret)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: a(:)
    integer :: ret(size(a)*self%ncells)
    integer :: n, i
    n=size(a)
    do i =1, self%ncells
       ret((i-1)*n+1: i*n) = a(:)
    end do
  end subroutine repeat_int1d

  subroutine repeat_real1d(self, a, ret)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: a(:)
    real(dp):: ret(size(a)*self%ncells)
    integer :: n, i
    n=size(a)
    do i =1, self%ncells
       ret((i-1)*n+1: i*n) = a(:)
    end do
  end subroutine repeat_real1d

  subroutine test1()
    type(supercell_maker_t) :: maker
    integer :: scmat(3,3) 
    scmat=reshape([0,1, 1, 1,0, 1, 1,1,0], [3,3])
    call maker%initialize(scmat)
  end subroutine test1

  subroutine scmaker_unittest()
    call test1()
  end subroutine scmaker_unittest

end module m_supercell_maker
