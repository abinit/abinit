!{\src2tex{textfont=tt}}
!!****f* ABINIT/atomden
!! NAME
!! atomden
!!
!! FUNCTION
!! Construct atomic proto-bulk density (i.e. the superposed density 
!! from neutral, isolated atoms at the bulk atomic positions).
!! This is useful if one wants to construct the bonding density:
!! 
!! rho^{bnd} = rho^{bulk}(r) 
!!                 - \sum_{\alpha}\rho^{atm}_{\alpha}(r-R_{\alpha})
!!
!! Where rho^{bulk} is the bulk density, rho^{atm} the atomic density
!! and the index \alpha sums over all atoms. the R_{\alpha} are the
!! atomic positions in the bulk. This routine calculates the sum over
!! rho^{atm}_{\alpha}(r-R_{\alpha}) on a grid.
!!
!! Units are atomic.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! calctype : type of calculation
!!          'replace' zero the input/output density array
!!          'add'     add to the input/output density array
!! natom : number of atoms in cell
!! ntypat : number of different types of atoms in cell
!! typat(natom) : type of each atom
!! ngrid : number of gridpoints
!! r_vec_grid(3,ngrid) : real (non-reduced) coordinates for grid points
!! rho(ngrid) : input/output density array
!! a(3),b(3),c(3) : real-space basis vectors
!! atom_pos(3,natom) : reduced coordinates for atomic positions
!! natomgr(ntypat) : number of gridpoints for each atomic density grid
!! natomgrmax : max(natomgr(ntypat))
!! atomrgrid(natomgrmax,ntypat)
!! density(natomgrmax,ntypat)
!!
!! OUTPUT
!! rho(ngrid) : input/output density array
!!
!! SIDE EFFECTS
!!
!! NOTES
!! There are two ways to compile the proto density in real space
!! for a solid. One alternative is that the density is calculated
!! for an extended grid encompassing the sphere of points around
!! one atom, and the results are folded back into the unit cell.
!! On the other hand one can, around each grid point, identify the
!! number of atoms in a sphere equivalent to the length of the radial
!! grid for each type of atom.
!! The second approach, with some modification, is taken here. The
!! numer of atoms in a supercell cell are listed such that the supercell
!! encompasses the atoms which could contribute to any point in the grid.
!! That list is kept and cycled through, to avoid recalculating it at
!! each point.
!! Note that the density calculated from the atom is the spherical
!! average, since there is no preferred direction without any
!! external field (and it's simpler)
!!
!!
!! PARENTS
!!      read_atomden
!!
!! CHILDREN
!!      sort_dp,spline,splint,wrtout,xmpi_barrier,xmpi_sum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,atom_pos, &
&                  natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_splines

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atomden'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat,ngrid,natomgrmax,prtvol
 character(len=7),intent(in) :: calctype
!arrays
 type(MPI_type),intent(in) :: MPI_enreg
 integer,intent(in) :: typat(natom),natomgr(ntypat)
 real(dp),intent(in) :: r_vec_grid(3,ngrid),a(3),b(3),c(3)
 real(dp),intent(in) :: atom_pos(3,natom),atomrgrid(natomgrmax,ntypat)
 real(dp),intent(in) :: density(natomgrmax,ntypat)
 real(dp),intent(inout) :: rho(ngrid)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 integer :: cnt,delta,i,l,m,n,iatom,itypat,igrid,n_cells,n_grid_p
 integer :: ierr,spaceComm,nprocs,master,rank,remainder
 real(dp) :: a_norm,b_norm,c_norm
 real(dp) :: r_max,R_sphere_max,dp_dummy,ybcbeg,ybcend
!arrays
 integer :: n_equiv_atoms(ntypat),grid_index(ngrid)
 integer :: my_start_equiv_atoms(ntypat)
 integer :: my_end_equiv_atoms(ntypat)
 integer :: l_min(ntypat),m_min(ntypat),n_min(ntypat)
 integer :: l_max(ntypat),m_max(ntypat),n_max(ntypat)
 real(dp) :: center(3),dp_vec_dummy(3),delta_a(3),delta_b(3),delta_c(3)
 real(dp) :: r_atom(3),grid_distances(ngrid)
 integer, allocatable :: new_index(:),i_1d_dummy(:)
 real(dp),allocatable :: equiv_atom_dist(:,:),equiv_atom_pos(:,:,:),rho_temp(:,:)
 real(dp),allocatable :: dp_1d_dummy(:),dp_2d_dummy(:,:),ypp(:)
 real(dp),allocatable :: x_fit(:),y_fit(:)


! ************************************************************************

!initialise and check parallel execution
 spaceComm=MPI_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 rank=MPI_enreg%me_kpt

 master=0
 
!initialise variables and vectors
 a_norm = sqrt(dot_product(a,a))
 b_norm = sqrt(dot_product(b,b))
 c_norm = sqrt(dot_product(c,c))
 center = (a+b+c)*half
 dp_dummy = dot_product(a,b)/(b_norm*b_norm)
 dp_vec_dummy = dp_dummy*b
 delta_a = a - dp_vec_dummy
 dp_dummy = dot_product(b,a)/(a_norm*a_norm)
 dp_vec_dummy = dp_dummy*a
 delta_b = b - dp_vec_dummy
 dp_dummy = dot_product(c,(a+b))/(dot_product((a+b),(a+b)))
 dp_vec_dummy = dp_dummy*(a+b)
 delta_c = c - dp_vec_dummy
 ABI_ALLOCATE(rho_temp,(ngrid,ntypat))
 rho_temp = zero

!write(std_out,*) '*** --- In atomden --- ***'
!write(std_out,*) ' a_norm:',a_norm,' b_norm:',b_norm,' c_norm:',c_norm
!write(std_out,*) 'delta_a:',delta_a,'delta_b:',delta_b,'delta_c:',delta_c
!write(std_out,*) ' center:',center

!Find supercell which will contain all possible contributions
!for all atoms, and enumerate positions for all atoms
!TODO list of atoms can be "pruned", i.e identify all atoms
!that can't possibly contribute and remove from list.
!Should be most important for very oblique cells
 do itypat=1,ntypat
   R_sphere_max = atomrgrid(natomgr(itypat),itypat)
   l_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_a,delta_a)))
   l_max(itypat) = -l_min(itypat)
   m_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_b,delta_b)))
   m_max(itypat) = -m_min(itypat) 
   n_min(itypat) = -ceiling(R_sphere_max/sqrt(dot_product(delta_c,delta_c)))
   n_max(itypat) = -n_min(itypat)
   n_cells = (l_max(itypat)-l_min(itypat)+1) &
&   *(m_max(itypat)-m_min(itypat)+1) &
&   *(n_max(itypat)-n_min(itypat)+1)
   n_equiv_atoms(itypat) = 0
   do iatom=1,natom
     if (typat(iatom)==itypat) then
       n_equiv_atoms(itypat) = n_equiv_atoms(itypat) + n_cells
     end if ! if type=itypat
   end do ! number of atoms per cell
   if ((rank==master).and.(prtvol>9)) then
     write(message,'(a)') '*** --- In atomden --- find box ***'
     call wrtout(std_out,message,'COLL')
     write(message,'(a,I4)') ' itypat:',itypat
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' l_min:',l_min(itypat),' l_max:',l_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' m_min:',m_min(itypat),' m_max:',m_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' n_min:',n_min(itypat),' n_max:',n_max(itypat)
     call wrtout(std_out,message,'COLL')
     write(message,'(2(a,I4))') ' n_equiv_atoms:',n_equiv_atoms(itypat)
     call wrtout(std_out,message,'COLL')
   end if
 end do !atom type

!allocate arrays
 n = maxval(n_equiv_atoms) 
 ABI_ALLOCATE(equiv_atom_pos,(3,n,ntypat))
 ABI_ALLOCATE(equiv_atom_dist,(n,ntypat))
 equiv_atom_pos = zero
 equiv_atom_dist = zero

!Find positions and distance of atoms from center of cell
 do itypat=1,ntypat
   i = 1
   do l=l_min(itypat),l_max(itypat)
     do m=m_min(itypat),m_max(itypat)
       do n=n_min(itypat),n_max(itypat)
         do iatom=1,natom
           if (typat(iatom)==itypat) then
             if (i>n_equiv_atoms(itypat)) then 
               MSG_ERROR('atomden: i>n_equiv_atoms')
             end if
             equiv_atom_pos(:,i,itypat) = (atom_pos(1,iatom)+dble(l))*a &
&             + (atom_pos(2,iatom)+dble(m))*b &
&             + (atom_pos(3,iatom)+dble(n))*c
             dp_vec_dummy = equiv_atom_pos(:,i,itypat)-center
             equiv_atom_dist(i,itypat) = &
&             sqrt(dot_product(dp_vec_dummy,dp_vec_dummy))
             i = i + 1
           end if
         end do
       end do !n
     end do !m
   end do !l
!  write(std_out,*) '*** --- In atomden --- find equiv ***'
!  write(std_out,*) ' itypat:',itypat
!  write(std_out,*) ' equiv_atom_pos:'
!  write(std_out,*) equiv_atom_pos(:,:,itypat)
!  write(std_out,*) ' equiv_atom_dist:',equiv_atom_dist(:,itypat)
 end do !atom type

!Sort the atoms after distance so that the density from the ones
!furthest away can be added first. This is to prevent truncation error.
 do itypat=1,ntypat
   n = n_equiv_atoms(itypat)
   ABI_ALLOCATE(dp_1d_dummy,(n))
   ABI_ALLOCATE(new_index,(n))
   ABI_ALLOCATE(dp_2d_dummy,(3,n))
   dp_1d_dummy = equiv_atom_dist(1:n,itypat)
   dp_2d_dummy = equiv_atom_pos(1:3,1:n,itypat)
   do i=1,n
     new_index(i) = i
   end do
   call sort_dp(n,dp_1d_dummy,new_index,tol14)
   do i=1,n
!    write(std_out,*) i,' -> ',new_index(i)
     equiv_atom_pos(1:3,n+1-i,itypat) = dp_2d_dummy(1:3,new_index(i))
     equiv_atom_dist(1:n,itypat) = dp_1d_dummy 
   end do
   ABI_DEALLOCATE(dp_1d_dummy)
   ABI_DEALLOCATE(new_index)
   ABI_DEALLOCATE(dp_2d_dummy)
!  write(std_out,*) '*** --- In atomden ---  sorting atoms ***'
!  write(std_out,*) ' itypat:',itypat
!  write(std_out,*) ' equiv_atom_pos:'
!  write(std_out,*) equiv_atom_pos(:,:,itypat) 
!  write(std_out,*) ' equiv_atom_dist:',equiv_atom_dist(:,itypat) 
 end do ! atom type

!Divide the work in case of parallel execution
 if (nprocs==1) then ! Make sure everything runs with one proc
   if (prtvol>9) then
     write(message,'(a)') '  In atomden - number of processors:     1'
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of proto-atomic density done in serial'
     call wrtout(std_out,message,'COLL')
   end if
   do itypat=1,ntypat
     if (prtvol>9) then
       write(message,'(a,I6)') '  Number of equivalent atoms:',n_equiv_atoms(itypat)
       call wrtout(std_out,message,'COLL')
     end if
     my_start_equiv_atoms(itypat) = 1
     my_end_equiv_atoms(itypat) = n_equiv_atoms(itypat)
   end do
 else
   if (rank==master.and.prtvol>9) then
     write(message,'(a,I5)') '  In atomden - number of processors:',nprocs
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of proto-atomic density done in parallel'
     call wrtout(std_out,message,'COLL')
   end if
   do itypat=1,ntypat
     if (rank==master.and.prtvol>9) then
       write(message,'(a,I6)') '  Number of equivalent atoms:',n_equiv_atoms(itypat)
       call wrtout(std_out,message,'COLL')
     end if
!    Divide the atoms among the processors by shuffling indices
     delta = int(floor(real(n_equiv_atoms(itypat))/real(nprocs)))
     remainder = n_equiv_atoms(itypat)-nprocs*delta
     my_start_equiv_atoms(itypat) = 1+rank*delta
     my_end_equiv_atoms(itypat) = (rank+1)*delta
!    Divide the remainder points among the processors
!    by shuffling indices
     if ((rank+1)>remainder) then
       my_start_equiv_atoms(itypat) = my_start_equiv_atoms(itypat) + remainder
       my_end_equiv_atoms(itypat) = my_end_equiv_atoms(itypat) + remainder
     else
       my_start_equiv_atoms(itypat) = my_start_equiv_atoms(itypat) + rank
       my_end_equiv_atoms(itypat) = my_end_equiv_atoms(itypat) + rank + 1
     end if
     if (prtvol>9) then
       write(message,'(a,I3)') '          For atom type: ',itypat
       call wrtout(std_out,message,'PERS')
!      write(message,'(a,I6)') '  I''ll take atoms from: ',my_start_equiv_atoms(itypat)
!      call wrtout(std_out,message,'PERS')
!      write(message,'(a,I6)') '           total for me: ',my_end_equiv_atoms(itypat)
!      call wrtout(std_out,message,'PERS')
       write(message,'(a,I6)') '            total for me: ', &
&       my_end_equiv_atoms(itypat)+1-my_start_equiv_atoms(itypat)
       call wrtout(std_out,message,'PERS')
     end if
   end do
 end if

!Loop over types of atoms and equivalent atoms and
!interpolate density onto grid
 do itypat=1,ntypat
!  do iatom=my_start_equiv_atoms(itypat),my_end_equiv_atoms(itypat)
   
   cnt = 0
   iatom = rank+1 - nprocs
!  Parallel execution of loop
   do
     cnt = cnt + 1
     iatom = iatom + nprocs
     if (iatom>n_equiv_atoms(itypat)) exit ! Exit if index is too large

     if (mod(cnt,100)==0.and.prtvol>0) then
       write(message,'(2(a,I6))') ' atoms so far',cnt,' of: ',n_equiv_atoms(itypat)/nprocs
       call wrtout(std_out,message,'PERS')
     end if
     
     r_max = atomrgrid(natomgr(itypat),itypat) 
     r_atom = equiv_atom_pos(:,iatom,itypat)

!    Set up an array with the gridpoint distances
     i = 1
     grid_distances = zero
     grid_index = 0
     do igrid=1,ngrid
       dp_vec_dummy(:) = r_vec_grid(:,igrid) - r_atom(:)
       dp_dummy = sqrt(dot_product(dp_vec_dummy,dp_vec_dummy))
       if (dp_dummy <= r_max) then
         grid_distances(i) = dp_dummy
         grid_index(i) = igrid      
         i = i + 1
       else
         cycle ! cycle if point is too far away
       end if
     end do
     n_grid_p = i - 1

     if (n_grid_p==0) cycle ! Cycle if no point needs
!    to be interpolated
     
!    Sort points to be interpolated in ascending order
     ABI_ALLOCATE(dp_1d_dummy,(n_grid_p))
     ABI_ALLOCATE(new_index,(n_grid_p))
     ABI_ALLOCATE(i_1d_dummy,(n_grid_p))
     dp_1d_dummy = grid_distances(1:n_grid_p)
     do i=1,n_grid_p
       new_index(i) = i 
     end do
     call sort_dp(n_grid_p,dp_1d_dummy,new_index,tol16)
     grid_distances(1:n_grid_p) = dp_1d_dummy
     i_1d_dummy = grid_index(1:n_grid_p)
     do i=1,n_grid_p
!      write(std_out,*) i_1d_dummy(i),' -> ',i_1d_dummy(new_index(i))
       grid_index(i) = i_1d_dummy(new_index(i))
     end do
     ABI_DEALLOCATE(dp_1d_dummy)
     ABI_DEALLOCATE(new_index)
     ABI_DEALLOCATE(i_1d_dummy)

!    Interpolate density onto all grid points
     ABI_ALLOCATE(ypp,(natomgr(itypat)))
     ABI_ALLOCATE(x_fit,(n_grid_p))
     ABI_ALLOCATE(y_fit,(n_grid_p))
     ypp = zero; y_fit = zero
     ybcbeg = zero; ybcend = zero
     x_fit = grid_distances(1:n_grid_p)
     call spline(atomrgrid(1:natomgr(itypat),itypat), &
&     density(1:natomgr(itypat),itypat), &
&     natomgr(itypat),ybcbeg,ybcend,ypp)
     call splint(natomgr(itypat),atomrgrid(1:natomgr(itypat),itypat), &
&     density(1:natomgr(itypat),itypat),ypp,n_grid_p, &
&     x_fit,y_fit)
     
!    Save the interpolated points to grid
     do i=1,n_grid_p
       rho_temp(grid_index(i),itypat) = rho_temp(grid_index(i),itypat) + y_fit(i)
     end do
     ABI_DEALLOCATE(ypp)
     ABI_DEALLOCATE(x_fit)
     ABI_DEALLOCATE(y_fit)

   end do ! n equiv atoms
 end do ! type of atom

!Collect all contributions to rho_temp if
!we are running in parallel
 if (nprocs>1) then
   call xmpi_barrier(spaceComm)
   call xmpi_sum_master(rho_temp,master,spaceComm,ierr)
   call xmpi_barrier(spaceComm)
   if (prtvol>9) then
     write(message,'(a)') '  In atomden - contributions to rho_temp collected'
     call wrtout(std_out,message,'PERS')
   end if
 end if

!Now rho_temp contains the atomic protodensity for each atom.
!Check whether this is to replace or be added to the input/output array
!and sum up contributions
 if (trim(calctype)=='replace') rho = zero
 do itypat=1,ntypat
   rho(:) = rho(:) + rho_temp(:,itypat)
 end do

!deallocations
 if (allocated(rho_temp))  then
   ABI_DEALLOCATE(rho_temp)
 end if
 if (allocated(equiv_atom_pos))  then
   ABI_DEALLOCATE(equiv_atom_pos)
 end if
 if (allocated(equiv_atom_dist))  then
   ABI_DEALLOCATE(equiv_atom_dist)
 end if
!if (allocated()) deallocate()

 return

 end subroutine atomden
!!***
