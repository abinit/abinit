!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_atomden
!! NAME
!! read_atomden
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! natom : number of atoms in cell
!! nspden : number of spin densities 
!! ntypat : number of types of atoms in the cell
!! pawfgr <type(pawfgr_type)> : data about the fine grid
!! typat(natom) : list of atom types
!!
!! OUTPUT
!! rhor_paw(pawfgr%nfft,nspden) : full electron density on the fine grid
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomden
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine read_atomden(MPI_enreg,natom,nspden,ntypat,pawfgr, &
&                       rhor_paw,typat,rprimd,xred,prtvol,file_prefix)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_io_tools, only : open_file
 use m_pawfgr,   only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_atomden'
 use interfaces_65_paw, except_this_one => read_atomden
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,ntypat,prtvol
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 type(MPI_type),intent(in) :: MPI_enreg
 integer,intent(in) :: typat(natom)
 real(dp), intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(inout) :: rhor_paw(pawfgr%nfft,nspden)
 character(len=7), intent(in) :: file_prefix

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 character(len=120) :: filename
 character(len=7) :: calctype='replace'
 integer :: igrid,i,i1,i2,i3,io_err,itypat,unt
 integer :: natomgrmax,nlines,ngrid,n1,n2,n3
 real(dp) :: difx,dify,difz,ucvol!,norm
!arrays
 integer :: natomgr(ntypat)
 real(dp) :: a(3),b(3),c(3)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: atomrgrid(:,:),r_vec_grid(:,:),density(:,:)
 real(dp),allocatable :: rho(:)

! ************************************************************************

!Initialise various variables
 ngrid = pawfgr%nfft
 a(:) = rprimd(:,1)
 b(:) = rprimd(:,2)
 c(:) = rprimd(:,3)
 ABI_ALLOCATE(rho,(ngrid))
 if (nspden/=1) then 
   MSG_ERROR('read_atomden: Only nspden=1 allowed.')
 end if
 rho = rhor_paw(:,1)
 gmet=zero;gprimd=zero;rmet=zero;ucvol=zero


!Calculate the r vector (reduced coord.) of the fine gridpoints
 ABI_ALLOCATE(r_vec_grid,(3,ngrid))
 igrid = 0
 n1 = pawfgr%ngfft(1)
 n2 = pawfgr%ngfft(2)
 n3 = pawfgr%ngfft(3)
 do i3=0,n3-1
   difz=dble(i3)/dble(n3)
   do i2=0,n2-1
     dify=dble(i2)/dble(n2)
     do i1=0,n1-1
       difx=dble(i1)/dble(n1)
       igrid = igrid + 1
       r_vec_grid(1,igrid)=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
       r_vec_grid(2,igrid)=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
       r_vec_grid(3,igrid)=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
     end do
   end do
 end do
 if (igrid/=ngrid) then 
   MSG_ERROR('read_atomden: igrid not equal to ngrid')
 end if

!Read in atomic density data for each atom type
!first check how many datapoints are in each file
 do itypat=1,ntypat
   filename='';io_err=0;
   if (itypat>0)  write(filename,'(a,a,i1,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'  
   if (itypat>10) write(filename,'(a,a,i2,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (open_file(filename, message, newunit=unt, status='old',action='read') /= 0) then
     write(std_out,*) 'ERROR in read_atomden: Could not open file: ',filename
     write(std_out,*) ' Current implementation requires this file to be present'
     write(std_out,*) ' for each type of atom.'
     write(std_out,*)trim(message)
     MSG_ERROR("Cannot continue")
   end if
!  Check number of lines in file
   nlines = 1;io_err=0;
   do
     read(unt,*,iostat=io_err)
     if (io_err<0) exit
     nlines = nlines + 1   
   end do
   close(unt)
   natomgr(itypat) = nlines - 2
 end do ! Atom type
!Allocate arrays and read in data
 natomgrmax = maxval(natomgr)
 ABI_ALLOCATE(atomrgrid,(natomgrmax,ntypat))
 ABI_ALLOCATE(density,(natomgrmax,ntypat))
 atomrgrid = zero ; density = zero
 do itypat=1,ntypat
   filename='';io_err=0;
   if (itypat>0)  write(filename,'(a,a,i1,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (itypat>10) write(filename,'(a,a,i2,a)') trim(file_prefix), &
&   '_density_atom_type',itypat,'.dat'
   if (open_file(filename,message,newunit=unt,status='old',action='read') /= 0) then
     MSG_ERROR(message)
   end if
   read(unt,*) ! Skip comment line
   do i=1,natomgr(itypat)
     read(unt,*) atomrgrid(i,itypat),density(i,itypat)
   end do
   close(unt)
   if (atomrgrid(1,itypat)/=zero) then
     write(std_out,*) 'ERROR in read_atomden, in file: ',filename
     write(std_out,*) ' First gridpoint has to be the origin.'
     MSG_ERROR("Cannot continue")
   end if
 end do ! Atom type

!write(std_out,*) '*** --- In read_atomden before call--- ***'
!write(std_out,*) '  calctype:',calctype,' natom:',natom
!write(std_out,*) '    ntypat:',ntypat,' typat:',typat
!write(std_out,*) '     ngrid:',ngrid
!write(std_out,*) '         a:',a
!write(std_out,*) '         b:',b
!write(std_out,*) '         c:',c
!write(std_out,*) '      xred:',xred
!write(std_out,*) '   natomgr:',natomgr
!write(std_out,*) 'natomgrmax:',natomgrmax
!write(std_out,*) ' atomrgrid:',atomrgrid
!write(std_out,*) '   density:',density
!write(std_out,*) 'r_vec_grid:'
!write(std_out,*) r_vec_grid

!Call atomden
 call atomden(MPI_enreg,natom,ntypat,typat,ngrid,r_vec_grid,rho,a,b,c,xred, &
& natomgr,natomgrmax,atomrgrid,density,prtvol,calctype)

!if (prtvol>9) then ! calculate norm
!call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
!norm = SUM(rho(:))*ucvol/dble(n1*n2*n3)
!write(message,'(a,F8.4)') '  In read_atomden - NORM OF DENSITY: ',norm
!call wrtout(std_out,message,'COLL')
!end if

 rhor_paw(:,1) = rho

 if (allocated(atomrgrid))  then
   ABI_DEALLOCATE(atomrgrid)
 end if
 if (allocated(density))  then
   ABI_DEALLOCATE(density)
 end if
 if (allocated(r_vec_grid))  then
   ABI_DEALLOCATE(r_vec_grid)
 end if
 if (allocated(rho))  then
   ABI_DEALLOCATE(rho)
 end if

 return

end subroutine read_atomden
!!***
