!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_supercell_maker
!! NAME
!! m_supercell_maker
!!
!! FUNCTION
!! This module define the supercell_maker file, which provide functions to help build
!! potentials in supercell.
!! Note: This module works at lower level than m_supercell. The purpose is not to get
!! an supercell of crystal structure, but it provide the functions to help building such thing (but not limited to).
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
  use m_mathfuncs , only: mat33det, binsearch_left_integerlist

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
     procedure :: trans_i
     procedure :: trans_j_and_Rj
     procedure :: trans_jlist_and_Rj
  end type supercell_maker_t

  public :: scmaker_unittest
contains

  subroutine initialize(self, sc_matrix)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: sc_matrix(3, 3)
    real(dp) :: tmp(3,3)
    self%scmat(:,:)=sc_matrix
    self%ncells=abs(mat33det(self%scmat))
    ABI_ALLOCATE(self%rvecs, (3, self%ncells))
    ! for inv
    tmp(:,:)=self%scmat(:,:)
    ! Why transpose?
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
    ret(:)=matmul(pos, self%inv_scmat)
  end function to_red_sc

  ! build R vectors for supercells.
  subroutine build_rvec(self)
    class(supercell_maker_t), intent(inout) :: self
    real(dp):: scorners_newcell(8,3), corners(8,3), x(3), tmp(3)
    integer :: rep(3), ix, iy, iz, counter
    real(dp) :: eps=1d-8
    counter=0
    scorners_newcell=transpose(reshape([0., 0., 0., 0., 0., 1., 0., 1., 0., &
         0., 1., 1., 1., 0., 0., 1., 0., 1.,  &
         1., 1., 0., 1., 1., 1.],[3,8]))
    corners = matmul(scorners_newcell, self%scmat)
    rep=ceiling(maxval(corners, dim=1) - minval(corners, dim=1))
    ! NOTE: DO NOT CHANGE THE ORDER. It is used in the binary search. 
    do ix = 0, rep(1)
       do iy = 0, rep(2)
          do iz = 0, rep(3)
             x(:)=[ix, iy, iz]
             tmp=self%to_red_sc(x)
             if ( (.not. any(tmp<=-1.0*eps)) &
                  .and. (.not. any(tmp>1.0-eps)) ) then
                counter=counter+1
                self%rvecs(:,counter)=nint(x)
             end if
          end do
       end do
    end do
    if (counter /= self%ncells ) then
       print *, "Wrong number of supercell found. Should have ", &
            self%ncells,"But found ", counter, "."
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
  ! TODO: The speed can be improved by using a bisect search instead.
  ! TODO: The speed can be improved by caching at least last n result, since in many cases,
  !       the next R is sill the same. 
  subroutine R_to_sc(self, R, R_sc, ind_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: R(3)
    integer, intent(inout) :: ind_sc, R_sc(3)
    integer :: rprim(3)
    R_sc=floor(self%to_red_sc(R*1.0d0))
    rprim(:)= R-matmul(R_sc, self%scmat)
    ind_sc=binsearch_left_integerlist(self%rvecs, rprim)
    !do i =0, self%ncells
    !  if (all(self%rvecs(:,i) == rprim)) then
    !      ind_sc=i
    !      exit
    !   end if
    !end do
  end subroutine R_to_sc

  ! ind: index in primitive cell
  ! ind_sc: indices (PLURAL) in supercell
  subroutine trans_i(self, nbasis, i, ind_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: i, nbasis
    integer, intent(out) :: ind_sc(self%ncells)
    integer :: icell
    do icell =1, self%ncells
       ind_sc(icell)=nbasis*(icell-1)+i
    end do
  end subroutine trans_i

  subroutine trans_ilist(self, nbasis, ilist, ilist_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: ilist(:), nbasis
    integer, intent(out) :: ilist_sc(size(ilist)*self%ncells)
    integer :: i
    do i =1, size(ilist)
       call self%trans_i(nbasis, ilist(i), ilist_sc(self%ncells*(i-1)+1:self%ncells*i ))
    end do
  end subroutine trans_ilist


  subroutine trans_j_and_Rj(self, nbasis, j, Rj, ind_sc, R_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: j, Rj(3), nbasis
    integer, intent(out) :: ind_sc(self%ncells), R_sc(3, self%ncells)
    integer :: i,jj
    do i =1, self%ncells
       call self%R_to_sc(Rj + self%rvecs(:,i), R_sc(:,i), jj)
       ind_sc(i)=nbasis*(jj-1)+j
    end do
  end subroutine trans_j_and_Rj


  subroutine trans_jlist_and_Rj(self, nbasis, jlist, Rj, ind_sc, R_sc)
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
  end subroutine trans_jlist_and_Rj


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

  function test1() result(err)
    type(supercell_maker_t) :: maker
    integer :: err
    integer :: scmat(3,3)
    integer :: ind_sc(8)
    integer :: R(3), R_sc(3), ind_sc2
    integer ::  j_sc(8), R_sc3(3, 8)
    scmat=reshape([2,0, 0, 0,2, 0, 0,0,2], [3,3])
    call maker%initialize(scmat)
    err=0
    call maker%trans_i(nbasis=3, i=1, ind_sc=ind_sc)

    R=[0,0,3]
    call maker%R_to_sc(R, R_sc, ind_sc2)
    if (.not.( (all(R_sc==[0,0,1])) .and. (ind_sc2==2 )))  then
       print *, "R_to_sc is wrong"
       err=1
    end if

    call maker%trans_j_and_Rj(nbasis=3, j=1, Rj=R, ind_sc=j_sc, R_sc=R_sc3  )
    if (.not. (  all(j_sc(1:3)==[4,1,10]) &
         .and. all(R_sc3(:,1)==[0,0,1]) &
         .and. all(R_sc3(:,3)==[0,0,1]) ) ) then
       print *, "Wrong trans_j_and_Rj"
       err=1
    end if
  end function test1

  function test2() result(err)
    type(supercell_maker_t) :: maker
    integer :: err
    integer :: scmat(3,3)
    real(dp) :: sccell(3,3)
    real(dp) :: scxred(3, 2*9)
    integer :: rep1(18)
    real(dp) :: rep2(18)
    scmat=transpose(reshape([1,2, 3, 4,5,6,7,8,6], [3,3]))
    !scmat=reshape([1,2, 3, 4,5,6,7,8,6], [3,3])
    call maker%initialize(scmat)
    !print *, "inv_scmat:",  maker%inv_scmat
    !print*, "Rvecs: ", maker%rvecs
 
    err=0

    ! test build_rvecs
    if (.not. all(maker%rvecs==reshape([0,  0,  0,  2,  3,  4,  3,  4,  4, &
         3,  4,  5,  4,  5,  5,  5,  6, 5, &
         5,  6,  6,  6,  7, 6,  8, 10, 10], [3, 8]))) then
       print *, "Wrong Rvecs found !"
       err=1
    end if

    ! test sc_cell
    sccell(:,:)=maker%sc_cell(transpose(reshape([1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,3.d0,2.d0,3.d0], [3,3])))
    if (.not. all(abs(sccell-transpose(reshape([18,18, 24, 42, 45, 60, 57, 66, 87], [3,3])))<1e-4)) then
       print *, "Wrong cell paramter found !"
       err=1
    end if

    ! test to_red_sc
    ! print*, maker%to_red_sc([1.0d0, 2.0d0,3.0d0])

    ! test trans_xred
    call maker%trans_xred(reshape([0.1d0, 0.2d0, 0.3d0, 0.1d0, 0.4d0, 0.6d0], [3, 2]), scxred)
    if (.not. all(abs(scxred(:,18)-[1.0666666666666684,0.53333333333333277,0.70000000000000040])<1e-4)) then
       print *, "Wrong sc_xred !"
       err=1
    end if

    ! test repeat
    call maker%repeat([1,2], rep1)
    print *, rep1

    call maker%repeat([1.0d0,2.0d0], rep2)
    print *, rep2

  end function test2




  subroutine scmaker_unittest()
    if (test1()/=0) print *, "Supercell maker: Test1 Failed"
    if (test2()/=0) print *, "Supercell maker: Test2 Failed"
  end subroutine scmaker_unittest

end module m_supercell_maker
