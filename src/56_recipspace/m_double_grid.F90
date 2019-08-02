!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_double_grid
!! NAME
!!  m_double_grid
!!
!! FUNCTION
!! This module defines the double grid object. This object contains the coarse mesh
!! and the dense mesh used for the interpolation of the BSE Hamiltonian,
!! and contains the mapping between the two meshes.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (YG, SP, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_double_grid

 use defs_basis
 use m_errors
 use m_abicore
 use m_hide_blas
 use m_bz_mesh
 use m_kptrank

 use m_numeric_tools,  only : wrap2_zero_one, interpol3d_indices
 use m_symtk,          only : matr3inv

 implicit none

 private
!!***

!!****t* m_double_grid/double_grid_t
!! NAME
!! double_grid_t
!!
!! FUNCTION
!! The double grid contains a coarse mesh and a dense mesh
!! It also contains the mapping between the two meshes
!!
!! SOURCE

 type,public :: double_grid_t

  integer :: kmult(3)
  ! Number of subdivisions in the coarse box in each direction

  integer :: ndiv
  ! Total number of small sub-boxes in the large box associated to the coarse k-mesh.
  ! i.e. product(kmult)

  integer :: maxcomp_coarse(3)
  ! Dimensions of the box containing the coarse points in integer coord

  integer :: nbz_coarse
  ! Number of k-points in the coarse BZ (open mesh)

  integer :: nbz_closedcoarse
  ! Number of k-points inside the coarse BZ (closed mesh)
  ! = PROD(maxcomp_coarse+1)

  integer :: nbz_dense
  ! Number of k-point inside the dense BZ (open mesh)
  ! = PROD(maxcomp_coarse.*kmult)

  integer, allocatable :: inttoik_coarse(:)
  ! inttoik_coarse(nbz_closedcoarse)
  ! Index of the kpoint in the coarse BZ.

  integer, allocatable :: iktoint_coarse(:)
  ! iktoint_coarse(nbz_coarse)
  ! Index of int by kpoint

  integer, allocatable :: inttoik_dense(:)

  integer, allocatable :: iktoint_dense(:)

  integer,allocatable :: indices_coarse(:,:)
  ! indices_coarse(3,nbz_closedcoarse)
  ! Indices (i1,i2,i3) for each point, ordinated by the integer coord

  integer,allocatable :: indices_dense(:,:)
  ! indices_dense(6,nbz_dense)
  ! Indices (i1,i2,i3);(j1,j2,j3) for each point, ordinated by integer coord

  integer :: kptrlatt_dense(3,3)
  ! kptrlatt of the dense mesh

  real(dp) :: klatt_dense(3,3)

  integer :: nshiftk_dense
  ! Number of shifts in the dense mesh.

  real(dp),allocatable :: shiftk_dense(:,:)
  ! shift_dense(3,nshiftk_dense)
  ! Shifts of the dense mesh.

  ! Dense lattice
  integer :: kptrlatt_coarse(3,3)
  ! kptrlatt of the coarse mesh

  real(dp) :: klatt_coarse(3,3)

  integer :: nshiftk_coarse
  ! Number of shifts in the coarse mesh.

  real(dp),allocatable :: shiftk_coarse(:,:)
  ! shift_coarse(3,nshiftk_coarse)
  ! Shifts of the coarse mesh.

  ! Coarse lattice
  integer, allocatable :: g0_coarse(:,:)
  integer, allocatable :: g0_dense(:,:)
  ! g0_dense/coarse(3,nkpt_closedcoarse/dense)
  ! G0 vector between the kpt obtained with indices
  ! and the kpt obtained insize bz

  integer, allocatable :: dense_to_coarse(:)
  ! dense_to_coarse(nbz_dense)
  ! Give the ibz_coarse corresponding to the dense mesh (the (0,0,0) point)

  integer, allocatable :: coarse_to_dense(:,:)
  ! coarse_to_dense(nbz_coarse,ndiv)
  ! Give all the ibz_dense corresponding to the (0,0,0) coarse point

 end type double_grid_t

 public :: double_grid_init             ! Initializes the double grid with coarse mesh and dense mesh read from file
 public :: double_grid_free             ! Deallocate all memory
 public :: get_kpt_from_indices_coarse  ! Returns the k-point index and g0 vector associated to the set of indices
 public :: compute_corresp              ! Compute correspondance data between k-dense and k-coarse
 !public :: get_kpt_from_indices_dense

 public :: kptfine_av                   ! Find the k-points of a fine grid that are around a k-point of a coarse mesh.
 public :: k_neighbors                  ! Find 8 neighbors of given k-point on a coarse grid, and return
!!***

!----------------------------------------------------------------------

CONTAINS  !=============================================================================
!!***

!!****f* m_double_grid/double_grid_init
!! NAME
!! double_grid_init
!!
!! FUNCTION
!! Initialize the double_grid datatype "grid" from coarse and dense mesh
!!
!! INPUTS
!!  Kmesh_coarse = descriptor of the coarse BZ sampling
!!  Kmesh_dense = descriptor of the dense BZ sampling
!!  kptrlatt_coarse(3,3) = vectors in R space that defines the reciprocal cell
!!  kmult(3) = multiplication factors from coarse to dense
!!
!! OUTPUT
!!  grid = double_grid to be created
!!
!! PARENTS
!!      setup_bse_interp
!!
!! CHILDREN
!!
!! SOURCE

subroutine double_grid_init(Kmesh_coarse,Kmesh_dense,kptrlatt_coarse,kmult,grid)

!Argument ------------------------------------
!scalars
 type(kmesh_t),intent(in) :: Kmesh_coarse,Kmesh_dense
 type(double_grid_t),intent(out) :: grid
!arrays
 integer,intent(in) :: kptrlatt_coarse(3,3),kmult(3)

!Local variables -----------------------------
!scalars
 integer :: ii, info
!arrays
 integer :: ipiv(3)
 real(dp) :: rlatt_coarse(3,3),klatt_coarse(3,3),curmat(3,3)

!*********************************************

 ABI_CHECK(Kmesh_coarse%nshift == 1, "Coarse mesh works only with nshiftk=1")
 ABI_CHECK(Kmesh_dense%nshift == 1, "Dense mesh : Works only with nshiftk=1")

 grid%nshiftk_coarse = Kmesh_coarse%nshift
 grid%nshiftk_dense = Kmesh_dense%nshift

 grid%nbz_coarse = Kmesh_coarse%nbz

 ABI_ALLOCATE(grid%shiftk_coarse,(3,grid%nshiftk_coarse))
 ABI_ALLOCATE(grid%shiftk_dense,(3,grid%nshiftk_dense))

 grid%shiftk_coarse(:,:) = Kmesh_coarse%shift(:,:)
 grid%shiftk_dense(:,:) = Kmesh_dense%shift(:,:)

 grid%kptrlatt_coarse(:,:) = kptrlatt_coarse(:,:)
 rlatt_coarse(:,:) = kptrlatt_coarse(:,:)
 call matr3inv(rlatt_coarse,klatt_coarse)
 grid%klatt_coarse(:,:) = klatt_coarse(:,:)

 grid%nbz_dense = Kmesh_dense%nbz
 grid%nbz_coarse = Kmesh_coarse%nbz

 grid%kmult(:) = kmult(:)
 grid%ndiv = kmult(1)*kmult(2)*kmult(3)

 ABI_ALLOCATE(grid%indices_dense,(6,Kmesh_dense%nbz))
 ABI_ALLOCATE(grid%g0_dense,(3,Kmesh_dense%nbz))
 ABI_ALLOCATE(grid%iktoint_dense,(Kmesh_dense%nbz))
 ABI_ALLOCATE(grid%inttoik_dense,(Kmesh_dense%nbz))

 grid%maxcomp_coarse(:) = -1

 curmat(:,:) = grid%kptrlatt_coarse(:,:)

 ! Gaussian elimination
 call dgetrf(3,3,curmat,3,ipiv,info)

 grid%nbz_closedcoarse = 1

 do ii = 1,3
   grid%maxcomp_coarse(ii) = ABS(NINT(curmat(ipiv(ii),ipiv(ii))))
   grid%nbz_closedcoarse = grid%nbz_closedcoarse*(grid%maxcomp_coarse(ii)+1)
 end do

 ABI_ALLOCATE(grid%indices_coarse,(3,grid%nbz_closedcoarse))
 ABI_ALLOCATE(grid%g0_coarse,(3,grid%nbz_closedcoarse))
 ABI_ALLOCATE(grid%iktoint_coarse,(Kmesh_coarse%nbz))
 ABI_ALLOCATE(grid%inttoik_coarse,(grid%nbz_closedcoarse))

 ! We should pass 'grid' at this stage !

 call create_indices_coarse(Kmesh_coarse%bz, Kmesh_coarse%nbz, grid%klatt_coarse, &
&    grid%nshiftk_coarse, grid%shiftk_coarse, grid%maxcomp_coarse, grid%nbz_closedcoarse, grid%indices_coarse, &
&    grid%g0_coarse,grid%iktoint_coarse,grid%inttoik_coarse)

 call create_indices_dense(grid%klatt_coarse, grid%maxcomp_coarse, Kmesh_dense%bz, Kmesh_dense%nbz, &
&    grid%nshiftk_dense, grid%shiftk_dense, grid%kmult, grid%indices_dense, grid%g0_dense, grid%iktoint_dense, &
&    grid%inttoik_dense)

 ABI_ALLOCATE(grid%dense_to_coarse,(Kmesh_dense%nbz))
 ABI_ALLOCATE(grid%coarse_to_dense,(Kmesh_coarse%nbz,grid%ndiv))

 call compute_neighbours(grid%nbz_dense, grid%iktoint_dense, grid%indices_dense, &
& grid%maxcomp_coarse, grid%inttoik_coarse, grid%g0_coarse, grid%nbz_closedcoarse, grid%nbz_coarse,&
& grid%ndiv, grid%dense_to_coarse, grid%coarse_to_dense)

end subroutine double_grid_init
!!***

!----------------------------------------------------------------------------

!!****f* m_double_grid/create_indices_coarse
!! NAME
!! create_indices_coarse
!!
!! FUNCTION
!!  Create mapping between kpoints and integer indexing
!!
!! INPUTS
!!  bz(3,nbz) = k-points in the Brillouin Zone
!!  nbz = number of k-points
!!  klatt(3,3) = reciprocal space vectors defining the reciprocal cell
!!  nshiftk = Number of shifts
!!  shiftk(3,nshiftk) = Shiftks of the Brillouin Zone
!!  maxcomp(3) = Maximum int along each direction
!!  nbz_closed = Number of k-points inside the closed Brillouin Zone (adding periodic images)
!!
!! OUTPUT
!!  indices(3,nbz_closed) = indices for each k-point in the closed BZ
!!  g0(3,nbz_closed) = g vectors between k-point inside bz and k-point given by indices
!!  iktoint(nbz) = mapping between k-points in the bz and int indices
!!  inttoik(nbz_closed) = mapping between int indices and k-points in the bz
!!
!! PARENTS
!!      m_double_grid
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_indices_coarse(bz, nbz, klatt, nshiftk, shiftk, maxcomp, nbz_closed, indices, g0, iktoint, inttoik)

!Argument ------------------------------------
!scalars
 integer,intent(in) :: nbz,nshiftk,nbz_closed
!arrays
 integer,intent(in) :: maxcomp(3)
 integer,intent(out) :: indices(3,nbz_closed)
 integer,intent(out) :: g0(3,nbz_closed)
 integer,intent(out) :: iktoint(nbz), inttoik(nbz_closed)
 real(dp),intent(in) :: bz(3,nbz),klatt(3,3),shiftk(3,nshiftk)

!Local variables -----------------------------
!scalars
 integer :: ik,ii,i1,i2, i3
 logical :: found
!arrays
 integer :: curg0(3)
 real(dp) :: curk1(3),ktoget(3)

!*********************************************

 ABI_CHECK(nshiftk==1,"nshiftk != 1 not supported")

 do i1 = 0,maxcomp(1)
   do i2 = 0,maxcomp(2)
     do i3 = 0,maxcomp(3)
       ii = (i1*(maxcomp(2)+1)+i2)*(maxcomp(3)+1)+i3+1
       ktoget(:) = shiftk(:,1)+(/i1,i2,i3/)
       curk1(:) = MATMUL(klatt(:,:),ktoget(:))
       found = .FALSE.
       do ik = 1,nbz
         if(isamek(curk1(:),bz(:,ik),curg0)) then
           indices(:,ii) = (/i1,i2,i3/)
           g0(:,ii) = curg0
           if (i1 /= maxcomp(1) .and. i2 /= maxcomp(2) .and. i3 /= maxcomp(3)) then
              iktoint(ik) = ii
           end if
           inttoik(ii) = ik
           found = .TRUE.
           exit
         end if
       end do
       if (.not. found) then
         write(std_out,*) "curk1 = ",curk1
         write(std_out,*) bz
         MSG_ERROR("A k-point generated from kptrlatt cannot be found in the BZ")
       end if
     end do
   end do
 end do

end subroutine create_indices_coarse
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/get_kpt_from_indices_coarse
!! NAME
!! get_kpt_from_indices_coarse
!!
!! FUNCTION
!!  Returns the k-point index and g0 vector associated to the set of indices
!!
!! INPUTS
!!  indices(3) = index of the searched k-point
!!  maxcomp(3) = Maximum int along each direction
!!  inttoik(nkpt) = mapping between int indices and k-points in the bz
!!  allg0(3,nkpt) = g vectors between k-point inside bz and k-point given by indices
!!  nkpt = number of k-points
!!
!! OUTPUT
!!  ikpt = index of k-point we search
!!  g0(3) = g-vector obtained
!!
!! PARENTS
!!      m_bseinterp,m_double_grid
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_kpt_from_indices_coarse(indices,maxcomp,inttoik,allg0,nkpt,ikpt,g0)

!Argument ------------------------------------
!scalars
 integer,intent(in) :: nkpt
 integer,intent(out) :: ikpt
!arrays
 integer,intent(in) :: indices(3),maxcomp(3)
 integer,intent(in) :: inttoik(nkpt),allg0(3,nkpt)
 integer,intent(out) :: g0(3)

!Local variables -----------------------------
!scalars
 integer :: curicoord

!*********************************************

 curicoord = (indices(1)*(maxcomp(2)+1)+indices(2))*(maxcomp(3)+1)+indices(3)+1
 ikpt = inttoik(curicoord)
 g0 = allg0(:,curicoord)

end subroutine get_kpt_from_indices_coarse
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/create_indices_dense
!! NAME
!! create_indices_dense
!!
!! FUNCTION
!!  Create mapping between kpoints and integer indexing
!!
!! INPUTS
!!  klatt_coarse(3,3) = reciprocal space vectors defining the reciprocal cell of coarse BZ
!!  maxcomp(3) = Maximum int along each direction
!!  bz_dense(3,nbz_dense) = k-points in the dense BZ
!!  nbz_dense = number of k-points in the dense BZ
!!  nshiftk = Number of shifts
!!  shiftk(3,nshiftk) = Shiftks of the Brillouin Zone
!!  kmult(3) = multiplication factors
!!  nbz_coarse = number of k-points in the coarse BZ
!!  kptrlatt_coarse(3,3) = real space vectors defining the reciprocal cell of coarse BZ
!!
!! OUTPUT
!!  indices(6,nbz_dense) = indices for each k-point in the closed BZ
!!  g0(3,nbz_dense) = g vectors between k-point inside bz and k-point given by indices
!!  iktoint(nbz_dense) = mapping between k-points in the bz and int indices
!!  inttoik(nbz_dense) = mapping between int indices and k-points in the bz
!!
!! PARENTS
!!      m_double_grid
!!
!! CHILDREN
!!
!! SOURCE

subroutine create_indices_dense(klatt_coarse, maxcomp, &
& bz_dense, nbz_dense, nshiftk, shiftk, kmult, indices, g0, inttoik, iktoint)

!Argument ------------------------------------
!scalars
 integer,intent(in) :: nbz_dense, nshiftk
!arrays
 integer,intent(in) :: kmult(3),maxcomp(3)
 integer,intent(out) :: indices(6,nbz_dense),g0(3,nbz_dense)
 integer,intent(out) :: inttoik(nbz_dense),iktoint(nbz_dense)
 real(dp),intent(in) :: bz_dense(3,nbz_dense),klatt_coarse(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)

!Local variables -----------------------------
 integer :: ik,ii,ii_coarse
 integer :: i1,i2,i3,j1,j2,j3
 logical :: found
!arrays
 integer :: curg0(3)
 real(dp) :: curk1(3),ktoget(3)

!*********************************************

 call wrtout(std_out, "Create Indices Dense", "COLL")

 ABI_CHECK(nshiftk==1,"nshiftk != 1 not supported")

 do i1 = 0,maxcomp(1)-1
   do i2 = 0,maxcomp(2)-1
     do i3 = 0,maxcomp(3)-1
       ii_coarse = (i1*(maxcomp(2)+1)+i2)*(maxcomp(3)+1)+i3+1

       do j1 = 0,kmult(1)-1
         do j2 = 0,kmult(2)-1
           do j3 = 0,kmult(3)-1
             ii = ((i1*kmult(1)+j1)*(maxcomp(2)*kmult(2)) +&
&                  (i2*kmult(2)+j2))*(maxcomp(3)*kmult(3))+&
&                  (i3*kmult(3)+j3)+1

             ktoget(1) = i1+((REAL(j1)+shiftk(1,1))/kmult(1))
             ktoget(2) = i2+((REAL(j2)+shiftk(2,1))/kmult(2))
             ktoget(3) = i3+((REAL(j3)+shiftk(3,1))/kmult(3))

             curk1(:) = MATMUL(klatt_coarse(:,:),ktoget(:))
             found = .FALSE.
             do ik = 1,nbz_dense
               if(isamek(curk1(:),bz_dense(:,ik),curg0)) then
                 indices(:,ii) = (/i1,i2,i3,j1,j2,j3/)
                 g0(:,ii) = curg0
                 inttoik(ii) = ik
                 iktoint(ik) = ii
                 found = .TRUE.
                 exit
               end if
             end do
             if(.not. found) then
               write(std_out,*) "curk1 = ",curk1
               write(std_out,*) bz_dense
               MSG_ERROR("Problem when creating indices")
             end if
           end do
         end do
       end do

     end do
   end do
 end do

end subroutine create_indices_dense
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/get_kpt_from_indices_dense
!! NAME
!! get_kpt_from_indices_coarse
!!
!! FUNCTION
!!  Returns the k-point index and g0 vector associated to the set of indices
!!
!! INPUTS
!!  indices(6) = index of the searched k-point
!!  maxcomp(3) = Maximum int along each direction
!!  kmult(3) = multiplication factors
!!  inttoik(nkpt) = mapping between int indices and k-points in the bz
!!  allg0(3,nkpt) = g vectors between k-point inside bz and k-point given by indices
!!  nkpt = number of k-points
!!
!! OUTPUT
!!  ikpt = index of k-point we search
!!  g0(3) = g-vector obtained
!!
!! PARENTS
!!
!! SOURCE

subroutine get_kpt_from_indices_dense(indices,maxcomp,kmult,inttoik,allg0,nkpt,ikpt,g0)

!Argument ------------------------------------
!scalars
 integer, intent(in) :: nkpt
 integer, intent(out) :: ikpt
!arrays
 integer, intent(in) :: indices(6),maxcomp(3),inttoik(nkpt)
 integer, intent(in) :: allg0(3,nkpt),kmult(3)
 integer, intent(out) :: g0(3)

!Local variables -----------------------------
!scalars
 integer :: curicoord

!*********************************************

 curicoord = ((indices(1)*kmult(1)+indices(4))*(maxcomp(2)*kmult(2))+&
&             (indices(2)*kmult(2)+indices(5)))*(maxcomp(3)*kmult(3))+&
&             (indices(3)*kmult(3)+indices(6))+1

 ikpt = inttoik(curicoord)
 g0 = allg0(:,curicoord)

end subroutine get_kpt_from_indices_dense
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/compute_neighbours
!! NAME
!! compute_neighbours
!!
!! FUNCTION
!! Compute correspondance between points in the dense BZ and in the coarse BZ
!!
!! INPUTS
!!   nbz_dense, nbz_closedcoarse, nbz_coarse, ndiv
!!   iktoint_dense(nbz_dense)
!!   indices_dense(6,nbz_dense)
!!   maxcomp_coarse(3)
!!   inttoik_coarse(nbz_closedcoarse)
!!   g0_coarse(3,nbz_closedcoarse)
!!
!! OUTPUT
!!  dense_to_coarse(nbz_dense)
!!  coarse_to_dense(nbz_coarse,ndiv)
!!
!! PARENTS
!!      m_double_grid
!!
!! CHILDREN
!!
!! SOURCE

subroutine compute_neighbours(nbz_dense, iktoint_dense, indices_dense, maxcomp_coarse, &
&  inttoik_coarse, g0_coarse, nbz_closedcoarse, nbz_coarse, ndiv, dense_to_coarse, coarse_to_dense)

!Argument ------------------------------------
!scalars
 integer,intent(in) :: nbz_dense, nbz_closedcoarse, nbz_coarse, ndiv
!arrays
 integer,intent(in) :: iktoint_dense(nbz_dense)
 integer,intent(in) :: indices_dense(6,nbz_dense)
 integer,intent(in) :: maxcomp_coarse(3)
 integer,intent(in) :: inttoik_coarse(nbz_closedcoarse)
 integer,intent(in) :: g0_coarse(3,nbz_closedcoarse)
 integer,intent(out) :: dense_to_coarse(nbz_dense)
 integer,intent(out) :: coarse_to_dense(nbz_coarse,ndiv)

!Local variables -----------------------------
!scalars
 integer :: ik_dense, iorder, ik_coarse
!arrays
 integer :: curindex(nbz_coarse)
 integer :: curindices_dense(6), curindices_coarse(3)
 integer :: g0(3)

!*********************************************

 DBG_ENTER("COLL")

 coarse_to_dense = 1
 dense_to_coarse = 1

 curindex = 1
 do ik_dense = 1, nbz_dense
  ! From ik_ibz in the dense mesh -> indices_dense
  iorder = iktoint_dense(ik_dense)

  ! From indices_dense -> indices_coarse
  curindices_dense = indices_dense(:,iorder)
  curindices_coarse = curindices_dense(1:3)
  ! From indices_coarse -> ik_ibz in the coarse mesh
  call get_kpt_from_indices_coarse(curindices_coarse,maxcomp_coarse,&
&   inttoik_coarse,g0_coarse,nbz_closedcoarse,ik_coarse,g0)

  dense_to_coarse(ik_dense) = ik_coarse
  coarse_to_dense(ik_coarse, curindex(ik_coarse)) = ik_dense

  curindex(ik_coarse) = curindex(ik_coarse) + 1
 end do

 DBG_EXIT("COLL")

end subroutine compute_neighbours
!!***

!---------------------------------------------------------------------

!!****f* m_double_grid/compute_corresp
!! NAME
!! compute_corresp
!!
!! FUNCTION
!! Pre-process tables with mapping between divisions and coarse k-points
!!
!! INPUTS
!! double_grid
!!
!! OUTPUT
!! div2kdense(double_grid%nbz_coarse,double_grid%ndiv)
!! (k_coarse,idiv) -> k_dense
!! kdense2div(double_grid%nbz_dense)
!! k_dense -> idiv
!!
!! PARENTS
!!      m_hexc
!!
!! CHILDREN
!!
!! SOURCE

subroutine compute_corresp(double_grid, div2kdense, kdense2div)

!Argument ------------------------------------
!scalars
 type(double_grid_t),intent(in) :: double_grid
!arrays
 integer,intent(out) :: div2kdense(double_grid%nbz_coarse,double_grid%ndiv)
 integer,intent(out) :: kdense2div(double_grid%nbz_dense)

!Local variables -----------------------------
!scalars
 integer :: iorder,ik_dense,ik_coarse
!arrays
 integer :: curindices_dense(6)
 integer,allocatable :: curindex(:)

!*********************************************
 ABI_MALLOC(curindex,(double_grid%nbz_coarse))
 curindex = 1

 do ik_dense = 1,double_grid%nbz_dense

   ! From ik_ibz in the dense mesh -> indices_dense
   iorder = double_grid%iktoint_dense(ik_dense)
   !g01 = double_grid%g0_dense(:,iorder)

   ! From indices_dense -> indices_coarse
   curindices_dense = double_grid%indices_dense(:,iorder)

   ik_coarse = double_grid%dense_to_coarse(ik_dense)
   div2kdense(ik_coarse,curindex(ik_coarse)) = ik_dense
   kdense2div(ik_dense) = curindex(ik_coarse)

   curindex(ik_coarse) = curindex(ik_coarse) + 1

 end do

 ABI_FREE(curindex)

end subroutine compute_corresp
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/double_grid_free
!! NAME
!! double_grid_free
!!
!! FUNCTION
!! Deallocate all dynamics entities present in a double_grid structure.
!!
!! INPUTS
!! grid<double_grid>=The datatype to be freed.
!!
!! SIDE EFFECTS
!! All allocated memory is released.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!
!! SOURCE

subroutine double_grid_free(grid)

!Arguments ------------------------------------
 type(double_grid_t),intent(inout) :: grid

! *********************************************************************

!integer
 ABI_SFREE(grid%inttoik_coarse)
 ABI_SFREE(grid%inttoik_dense)
 ABI_SFREE(grid%iktoint_coarse)
 ABI_SFREE(grid%iktoint_dense)
 ABI_SFREE(grid%indices_coarse)
 ABI_SFREE(grid%indices_dense)
 ABI_SFREE(grid%g0_coarse)
 ABI_SFREE(grid%g0_dense)
 ABI_SFREE(grid%dense_to_coarse)
 ABI_SFREE(grid%coarse_to_dense)

!real
 ABI_SFREE(grid%shiftk_dense)
 ABI_SFREE(grid%shiftk_coarse)

end subroutine double_grid_free
!!***

!----------------------------------------------------------------------

!!****f* m_double_grid/kptfine_av
!! NAME
!! kptfine_av
!!
!! FUNCTION
!! Find the k-points of a fine grid that are around a k-point of a coarse mesh.
!!
!! INPUTS
!!  center(3) = the point of the coarse mesh around which you want know which
!!              k-points of the fine mesh belong to.
!!  qptrlatt(3,3) = qptrlatt of the considered calculation (this is obtained
!!              from the input variable ngqpt and shiftq.
!!  kpt_fine(3,nkpt_fine) = this table contain all the k-points of the fine grid
!!              in the full BZ (no sym op. allowed) and is read from the header
!!              of the dense WF file.
!!  nkpt_fine = number of k-points of the fine grid read from the header of the
!!              dense WF file.
!!
!! OUTPUT
!!  kpt_fine_sub(nkpt_sub) = k-points of the fine grid that are around center(3)
!!  nkpt_sub = number of k-points of the fine grid that are around center(3)
!!  wgt_sub(nkpt_sub) = weight of the k-points of the fine grid that are around center(3).
!!
!! PARENTS
!!      eig2stern,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

subroutine kptfine_av(center,qptrlatt,kpt_fine,nkpt_fine,kpt_fine_sub,nkpt_sub,wgt_sub)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: nkpt_fine
 integer,intent(out)  :: nkpt_sub
!arrays
 integer,intent(in)   :: qptrlatt(3,3)
 real(dp),intent(in)  :: kpt_fine(3,nkpt_fine)
 real(dp),intent(in)  :: center(3)
 integer,pointer      :: kpt_fine_sub(:)
 real(dp),pointer     :: wgt_sub(:)

!Local variables-------------------------------
!scalars
 integer :: ikpt,aa,bb,cc
 integer :: ii,jj
!arrays
 real(dp) :: center_ref(3)
 real(dp) :: kpt_fine_ref(3)
 real(dp) :: kpt_tmp(3),kpt_tmp2(3)
 integer,allocatable  :: kpt_fine_sub_tmp(:)
 real(dp),allocatable :: wgt_sub_tmp(:)
 logical :: found(3)

! *************************************************************************

 ABI_ALLOCATE(kpt_fine_sub_tmp,(nkpt_fine))
 ABI_ALLOCATE(wgt_sub_tmp,(nkpt_fine))

!It is easier to work in real space using the qptrlatt matrices because in this
!referential any k-points sampling will be cast into an orthorhombic shape.
!In that space we can simply take all k-points of the fine grid that between
!center_ref-0.5 and center_ref+0.5

 center_ref = MATMUL(qptrlatt,center)

!When considering points center(3) that lying close or on a BZ edge we need to
!take the k-points of the fine grid taking into account unklamp vectors. This
!is done with the aa, bb and cc loops.

 ii = 1
 do ikpt=1,nkpt_fine
   kpt_tmp = kpt_fine(:,ikpt)
   do aa=-1,1
     kpt_tmp2(1) = kpt_tmp(1)+aa
     do bb=-1,1
       kpt_tmp2(2) = kpt_tmp(2)+bb
       do cc=-1,1
         kpt_tmp2(3) = kpt_tmp(3)+cc
         kpt_fine_ref = MATMUL(qptrlatt,kpt_tmp2)
         if((kpt_fine_ref(1)>=center_ref(1)-0.5-tol8).and.&
&         (kpt_fine_ref(1)<=center_ref(1)+0.5+tol8)) then
           if((kpt_fine_ref(2)>=center_ref(2)-0.5-tol8).and.&
&           (kpt_fine_ref(2)<=center_ref(2)+0.5+tol8)) then
             if((kpt_fine_ref(3)>=center_ref(3)-0.5-tol8).and.&
&             (kpt_fine_ref(3)<=center_ref(3)+0.5+tol8)) then
               kpt_fine_sub_tmp(ii) = ikpt
               ii = ii +1
             end if
           end if
         end if
       end do
     end do
   end do
 end do

 nkpt_sub = ii-1
 ABI_ALLOCATE(kpt_fine_sub,(nkpt_sub))
 ABI_ALLOCATE(wgt_sub,(nkpt_sub))

 do jj=1,nkpt_sub
   kpt_fine_sub(jj) = kpt_fine_sub_tmp(jj)
 end do

!We then compute a weight function. This weight function is simply a
!rectangular weight function that take the value 1 for k-points of the fine
!grid inside the cube, 0.5 for k-points that are lying on one face of the cube,
!0.25 for k-points that are lying on an edge of the cube and 0.125 for k-points
!that are lying on a peak of the cube.

 wgt_sub(:) = 1.0

 do ikpt=1,nkpt_sub
   found(:) = .True.
   kpt_tmp = kpt_fine(:,kpt_fine_sub(ikpt))
   do aa=-1,1
     kpt_tmp2(1) = kpt_tmp(1)+aa
     do bb=-1,1
       kpt_tmp2(2) = kpt_tmp(2)+bb
       do cc=-1,1
         kpt_tmp2(3) = kpt_tmp(3)+cc
         kpt_fine_ref = MATMUL(qptrlatt,kpt_tmp2)
         if((ABS(kpt_fine_ref(1)-center_ref(1)-0.5)< tol8) .or.&
&         (ABS(kpt_fine_ref(1)-center_ref(1)+0.5) < tol8)) then
           if(found(1)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(1) = .False.
           end if
         end if
         if((ABS(kpt_fine_ref(2)-center_ref(2)-0.5) < tol8) .or.&
&         (ABS(kpt_fine_ref(2)-center_ref(2)+0.5) < tol8)) then
           if(found(2)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(2) = .False.
           end if
         end if
         if((ABS(kpt_fine_ref(3)-center_ref(3)-0.5)< tol8) .or.&
&         (ABS(kpt_fine_ref(3)-center_ref(3)+0.5) < tol8)) then
           if(found(3)) then
             wgt_sub(ikpt) = wgt_sub(ikpt)*0.5
             found(3) = .False.
           end if
         end if
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(kpt_fine_sub_tmp)
 ABI_DEALLOCATE(wgt_sub_tmp)

end subroutine kptfine_av
!!***

!!****f* m_double_grid/k_neighbors
!!
!! NAME
!!   k_neighbors
!!
!! FUNCTION
!!   find 8 neighbors of given k-point on a coarse grid, and return
!!   them along with relative k-shift within coarse grid cell
!!
!! INPUTS
!!   kpt        = k-point to be interpolated to, in full BZ
!!   kptrlatt   = lattice vectors for coarse k-grid
!!   invrankkpt = rank list to find k-points
!!
!! OUTPUT
!!   rel_kpt = k-point coordinates renormalized to coarse grid cell
!!   kpt_phon_indices = indices of k-points on corners of cell
!!
!! TODO
!!  This routine is not used anymore. Deprecate or Remove?
!!
!! PARENTS
!!
!! CHILDREN
!!      get_rank_1kpt,interpol3d_indices,wrap2_zero_one
!!
!! SOURCE

subroutine k_neighbors (kpt, kptrlatt,krank, rel_kpt, kpt_phon_indices)

! inputs
 real(dp), intent(in) :: kpt(3)
 integer, intent(in) :: kptrlatt(3,3)
 type(krank_t), intent(in) :: krank

! outputs
 real(dp), intent(out) :: rel_kpt(3)
 integer, intent(out) :: kpt_phon_indices(8)

! local vars
 integer :: symrankkpt
 integer :: ir1,ir2,ir3, pr1,pr2,pr3
 real(dp) :: redkpt(3), cornerkpt(3), res

! *************************************************************************

!wrap fine kpt to [0,1]
 call wrap2_zero_one(kpt(1),redkpt(1),res)
 call wrap2_zero_one(kpt(2),redkpt(2),res)
 call wrap2_zero_one(kpt(3),redkpt(3),res)
!find 8 indices of points neighboring ikpt_phon, for interpolation
 call interpol3d_indices (redkpt,kptrlatt(1,1),kptrlatt(2,2),kptrlatt(3,3), &
& ir1,ir2,ir3, pr1,pr2,pr3)

!transpose ir pr to ikpt_phon indices
!order of kpt_phons:
!ir1 ir2 ir3
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt(cornerkpt)
 kpt_phon_indices(1) = krank%invrank(symrankkpt)
!pr1 ir2 ir3
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(2) = krank%invrank(symrankkpt)
!ir1 pr2 ir3
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(3) = krank%invrank(symrankkpt)
!pr1 pr2 ir3
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(4) = krank%invrank(symrankkpt)
!ir1 ir2 pr3
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(5) = krank%invrank(symrankkpt)
!pr1 ir2 pr3
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(6) = krank%invrank(symrankkpt)
!ir1 pr2 pr3
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(7) = krank%invrank(symrankkpt)
!pr1 pr2 pr3
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 symrankkpt = krank%get_rank_1kpt (cornerkpt)
 kpt_phon_indices(8) = krank%invrank(symrankkpt)

!retrieve the gkq matrix for all q, at the neighbor k vectors
 rel_kpt(1) = redkpt(1)*kptrlatt(1,1)-real(ir1-1)
 rel_kpt(2) = redkpt(2)*kptrlatt(2,2)-real(ir2-1)
 rel_kpt(3) = redkpt(3)*kptrlatt(3,3)-real(ir3-1)

end subroutine k_neighbors
!!***

END MODULE m_double_grid
!!***

