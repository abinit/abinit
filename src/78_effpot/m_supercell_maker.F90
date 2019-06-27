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
!! TODO: update m_supercell based on this so that it can use non-diagonal supercell matrix.
!! Then this need to be moved to level same as m_supercell.
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
  use m_xmpi
  use m_symtk,    only : matr3inv
  use m_mathfuncs , only: mat33det, binsearch_left_integerlist
  use m_mpi_scheduler, only: init_mpi_info
  use m_supercell
  implicit none
  private
!!***

  type ,public :: supercell_maker_t
     integer :: scmat(3,3)  ! supercell matrix
     real(dp) :: inv_scmat(3,3) ! inverse of supercell matrix
     integer :: ncells         ! number of cells in supercell
     integer, allocatable :: rvecs(:,:) ! R vectors for cells in supercell. dim:(3, ncells)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: to_red_sc  
     procedure :: build_rvec
     procedure :: sc_cell
     procedure :: R_to_sc

     ! Translations: from one primitive cell element to ncell supercell elements
     ! with values unchanged (therefore just repeat ncells times)
     generic :: repeat => repeat_int1d, repeat_real1d, repeat_real2d, repeat_realmat
     procedure :: repeat_int1d
     procedure :: repeat_real1d
     procedure :: repeat_real2d
     procedure :: repeat_realmat

     ! Translations: from one primitive cell element to ncell supercell elements
     ! with values changed (different from repeat)
     procedure :: trans_xred
     procedure :: trans_xcart
     procedure :: trans_i
     procedure :: trans_ilist
     procedure :: trans_j_and_Rj
     procedure :: trans_jlist_and_Rj
     procedure :: rvec_for_each
  end type supercell_maker_t

  public :: scmaker_unittest
contains

  subroutine initialize(self, sc_matrix)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: sc_matrix(3, 3)
    real(dp) :: tmp(3,3)
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 

    self%scmat(:,:)=sc_matrix
    call xmpi_bcast(self%scmat, master, comm, ierr)
    self%ncells=abs(mat33det(self%scmat))
    ! call xmpi_bcast(self%ncells, master, comm, ierr)
    ABI_ALLOCATE(self%rvecs, (3, self%ncells))
    ! Why transpose?
    tmp(:,:)=transpose(self%scmat)
    call matr3inv(tmp, self%inv_scmat)
    call self%build_rvec()
    !call xmpi_bcast(self%inv_scmat, master, comm, ierr)
    !call xmpi_bcast(self%ncells, master, comm, ierr)
    !call xmpi_bcast(self%rvecs, master, comm, ierr)
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
    ! cornors
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
       MSG_ERROR("Bug found. supercell_maker: Wrong number of supercell found in build_rvec. ")
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
    real(dp), allocatable, intent(inout) :: scxred(:, :)
    integer :: npos, icell, ipos, counter=0
    npos=size(xred, dim=2)
    if (.not. allocated(scxred)) then
       ABI_ALLOCATE(scxred, (3,size(xred, dim=2)*self%ncells))
    end if
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
    R_sc=floor(self%to_red_sc(R*1.0d0))
    rprim(:)= R-matmul(R_sc, self%scmat)
    ind_sc=binsearch_left_integerlist(self%rvecs, rprim)
    if (ind_sc==0) then
            MSG_ERROR("Bug found. supercell_maker%R_to_sc: Cannot find rprim")
    end if
  end subroutine R_to_sc

  ! ind: index in primitive cell
  ! ind_sc: indices (PLURAL) in supercell
  subroutine trans_i(self, nbasis, i, i_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: i, nbasis
    integer, allocatable, intent(inout) :: i_sc(:)
    integer :: icell
    if(.not. allocated(i_sc)) then
       ABI_ALLOCATE(i_sc, (self%ncells))
    end if
    do icell =1, self%ncells
       i_sc(icell)=nbasis*(icell-1)+i
    end do
  end subroutine trans_i

  subroutine trans_i_noalloc(self, nbasis, i, i_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: i, nbasis
    integer, intent(inout) :: i_sc(:)
    integer :: icell
    do icell =1, self%ncells
       i_sc(icell)=nbasis*(icell-1)+i
    end do
  end subroutine trans_i_noalloc


  subroutine trans_ilist(self, nbasis, ilist, ilist_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: ilist(:), nbasis
    integer, allocatable,  intent(inout) :: ilist_sc(:)
    integer :: i
    if(.not. allocated(ilist_sc)) then
       ABI_ALLOCATE(ilist_sc, (size(ilist)*self%ncells))
    end if
    do i =1, size(ilist)
       call trans_i_noalloc(self,nbasis, ilist(i), ilist_sc(self%ncells*(i-1)+1:self%ncells*i ))
    end do
  end subroutine trans_ilist


  subroutine trans_j_and_Rj(self, nbasis, j, Rj, j_sc, Rj_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: j, Rj(3), nbasis
    integer, allocatable , intent(inout) :: j_sc(:), Rj_sc(:, :)
    integer :: i,jj
    if(.not. allocated(j_sc)) then
       ABI_ALLOCATE(j_sc, (self%ncells))
    endif
    if(.not. allocated(Rj_sc)) then
       ABI_ALLOCATE(Rj_sc, (3, self%ncells))
    endif
    do i =1, self%ncells
       call self%R_to_sc(Rj + self%rvecs(:,i), Rj_sc(:,i), jj)
       j_sc(i)=nbasis*(jj-1)+j
    end do
  end subroutine trans_j_and_Rj

  subroutine trans_j_and_Rj_noalloc(self, nbasis, j, Rj, j_sc, Rj_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: j, Rj(3), nbasis
    integer, intent(inout) :: j_sc(:), Rj_sc(:, :)
    integer :: i,jj
    do i =1, self%ncells
       call self%R_to_sc(Rj + self%rvecs(:,i), Rj_sc(:,i), jj)
       j_sc(i)=nbasis*(jj-1)+j
    end do
  end subroutine trans_j_and_Rj_noalloc

  subroutine trans_jlist_and_Rj(self, nbasis, jlist, Rj, ind_sc, R_sc)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: jlist(:), Rj(3), nbasis
    integer, allocatable, intent(inout) :: ind_sc(:), R_sc(:,:)
    integer :: i,jj, counter, indj
    if (.not. allocated(ind_sc)) then
       ABI_ALLOCATE(ind_sc, (self%ncells*size(jlist)) )
    endif
    if (.not. allocated(R_sc)) then
       ABI_ALLOCATE(R_sc, (3, self%ncells))
    endif
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
    integer, allocatable :: ret(:)
    integer :: n, i
    n=size(a)
    if (.not. allocated(ret)) then
       ABI_ALLOCATE(ret, (n*self%ncells))
    end if
    do i =1, self%ncells
       ret((i-1)*n+1: i*n) = a(:)
    end do
  end subroutine repeat_int1d

  subroutine repeat_real1d(self, a, ret)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: a(:)
    real(dp), allocatable:: ret(:)
    integer :: n, i
    n=size(a)
    if (.not. allocated(ret)) then
       ABI_ALLOCATE(ret, (n*self%ncells))
    end if
    do i =1, self%ncells
       ret((i-1)*n+1: i*n) = a(:)
    end do
  end subroutine repeat_real1d

  subroutine repeat_real2d(self, a, ret)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable:: ret(:, :)
    integer :: n1, n2, i
    n1=size(a, dim=1)
    n2=size(a,dim=2)
    if (.not. allocated(ret)) then
       ABI_ALLOCATE(ret, (n1, n2*self%ncells))
    end if
    do i =1, self%ncells
       ret(:,(i-1)*n2+1: i*n2) = a(:,:)
    end do
  end subroutine repeat_real2d

  subroutine repeat_realmat(self, a, ret)
    class(supercell_maker_t), intent(inout) :: self
    real(dp), intent(in) :: a(:,:,:)
    real(dp), allocatable, intent(inout) :: ret(:,:,:)
    integer :: n,  i
    n=size(a, 3)
    if (.not. allocated(ret)) then
       ABI_ALLOCATE(ret, (size(a, dim=1), size(a, dim=2), n*self%ncells))
    end if
    do i=1, self%ncells
       ret(:,:,(i-1)*n+1: i*n)=a(:, :, :)
    end do
  end subroutine repeat_realmat

  subroutine rvec_for_each(self, nbasis, ret)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: nbasis
    integer, allocatable, intent(inout) :: ret(:,:)
    integer :: icell, ibasis
    if (.not. allocated(ret)) then
       ABI_ALLOCATE(ret, (3, nbasis*self%ncells))
    end if
    do icell=1, self%ncells
       do ibasis=1, nbasis
          ret(:, (icell-1)*nbasis+ibasis)= self%rvecs(:, icell)
       end do
    end do
  end subroutine rvec_for_each

  !============================Unit test=====================================

  function test1() result(err)
    type(supercell_maker_t) :: maker
    integer :: err
    integer :: scmat(3,3)
    integer, allocatable :: ind_sc(:)
    integer :: R(3), R_sc(3), ind_sc2
    integer, allocatable ::  j_sc(:), R_sc3(:, :)
    scmat=reshape([2,0, 0, 0,2, 0, 0,0,2], [3,3])
    call maker%initialize(scmat)
    err=0
    call maker%trans_i(nbasis=3, i=1, i_sc=ind_sc)
    R=[0,0,0]
    call maker%R_to_sc(R, R_sc, ind_sc2)
    if (.not.( (all(R_sc==[0,0,0])) .and. (ind_sc2==1 )))  then
       MSG_ERROR("R_to_sc is wrong")
    end if

    R=[0,0,3]
    call maker%trans_j_and_Rj(nbasis=3, j=1, Rj=R, j_sc=j_sc, Rj_sc=R_sc3  )
    if (.not. (  all(j_sc(1:3)==[4,1,10]) &
         .and. all(R_sc3(:,1)==[0,0,1]) &
         .and. all(R_sc3(:,3)==[0,0,1]) ) ) then
       MSG_ERROR("Wrong trans_j_and_Rj")
       err=1
    end if
    ABI_DEALLOCATE(ind_sc)
    ABI_DEALLOCATE(j_sc)
    ABI_DEALLOCATE(R_sc3)
  end function test1

  function test2() result(err)
    type(supercell_maker_t) :: maker
    integer :: err
    integer :: scmat(3,3)
    real(dp) :: sccell(3,3)
    real(dp), allocatable :: scxred(:, :)
    integer, allocatable :: rep1(:)
    real(dp), allocatable :: rep2(:)
    scmat=transpose(reshape([1,2, 3, 4,5,6,7,8,6], [3,3]))
    !scmat=reshape([1,2, 3, 4,5,6,7,8,6], [3,3])
    call maker%initialize(scmat)
 
    err=0

    ! test build_rvecs
    if (.not. all(maker%rvecs==reshape([0,  0,  0,  2,  3,  4,  3,  4,  4, &
         3,  4,  5,  4,  5,  5,  5,  6, 5, &
         5,  6,  6,  6,  7, 6,  8, 10, 10], [3, 8]))) then
       MSG_ERROR("Wrong Rvecs found !")
       err=1
    end if

    ! test sc_cell
    sccell(:,:)=maker%sc_cell(transpose(reshape([1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,3.d0,2.d0,3.d0], [3,3])))
    if (.not. all(abs(sccell-transpose(reshape([18,18, 24, 42, 45, 60, 57, 66, 87], [3,3])))<1e-4)) then
       MSG_ERROR("Wrong cell paramter found !")
       err=1
    end if


    ! test trans_xred
    call maker%trans_xred(reshape([0.1d0, 0.2d0, 0.3d0, 0.1d0, 0.4d0, 0.6d0], [3, 2]), scxred)
    if (.not. all(abs(scxred(:,18)-[1.0666666666666684,0.53333333333333277,0.70000000000000040])<1e-4)) then
       MSG_ERROR("Wrong sc_xred !")
       err=1
    end if

    ! test repeat
    call maker%repeat([1,2], rep1)

    call maker%repeat([1.0d0,2.0d0], rep2)

  end function test2

  subroutine scmaker_unittest()
    if (test1()/=0) MSG_ERROR("Supercell maker: Test1 Failed")
    if (test2()/=0) MSG_ERROR("Supercell maker: Test2 Failed")
  end subroutine scmaker_unittest


!===================================================================================
! The functions below are deprecated!!!
! They are used with the m_supercell module. 
  ! R (in term of primitive cell) to R_sc(in term of supercell) + R_prim
  subroutine find_R_PBC(scell, R, R_sc, R_prim)
    type(supercell_type) , intent(in):: scell
    integer, intent(in):: R(3)
    integer, intent(out):: R_sc(3), R_prim(3)
    real(dp) :: R_sc_d(3), sc_mat(3,3)

    integer:: ipriv(3), info
    !call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
    sc_mat(:,:)=scell%rlatt
    R_sc_d(:)=R(:)
    call dgesv(3, 1, sc_mat, 3, ipriv, R_sc_d, 3, info)
    if ( info/=0 ) then
       MSG_ERROR("Failed to find R_sc")
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
    MSG_ERROR("BUG found. supercell_maker%find_supercell_index cannot find iatom_prim, rvec pair in supercell")
  end function find_supercell_index

  !i0, j0+R0 shifted by R to i1=i0+0+R->periodic, j1=j0+R0+R->periodic
  subroutine find_supercell_ijR(scell, i0, j0, R0, R, i1, j1, R1, R_sc)

    type(supercell_type) , intent(in):: scell
    integer, intent(in) :: i0, j0, R0(3), R(3)
    integer, intent(out) :: i1, j1, R1(3), R_sc(3)
    i1=find_supercell_index(scell,i0,R)
    call find_R_PBC(scell,R0+R,R_sc,R1)
    j1=find_supercell_index(scell, j0, R1)
  end subroutine find_supercell_ijR



end module m_supercell_maker
