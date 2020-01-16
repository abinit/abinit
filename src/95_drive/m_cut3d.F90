!!****m* ABINIT/m_cut3d
!! NAME
!!  m_cut3d
!!
!! FUNCTION
!!  This module the predures used by cut3d
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (XG,MVerstraete,GMR,RC,LSI,JFB,MCote,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_cut3d

 use defs_basis
 use m_abicore
 use m_errors
 use m_splines
 use m_hdr
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_wfk
 use m_xmpi
 use m_sort
 use m_distribfft

 use defs_abitypes,      only : MPI_type
 use m_io_tools,         only : get_unit, iomode_from_fname, open_file, file_exists, read_string
 use m_numeric_tools,    only : interpol3d
 use m_symtk,            only : matr3inv
 use m_fstrings,         only : int2char10, sjoin, itoa
 use m_geometry,         only : xcart2xred, metric
 use m_special_funcs,    only : jlspline_t, jlspline_new, jlspline_free, jlspline_integral
 use m_pptools,          only : print_fofr_ri, print_fofr_xyzri , print_fofr_cube
 use m_mpinfo,           only : destroy_mpi_enreg, initmpi_seq
 use m_cgtools,          only : cg_getspin
 use m_gsphere,          only : getkpgnorm
 use m_epjdos,           only : recip_ylm, dens_in_sph
 use m_dens,             only : dens_hirsh
 use m_kg,               only : kpgio, ph1d3d, getph
 use m_fftcore,          only : sphereboundary
 use m_initylmg,         only : initylmg
 use m_fft,              only : fourwf

 implicit none

 private

 public :: cut3d_hirsh
 public :: cut3d_rrho
 public :: cut3d_volumeint
 public :: cut3d_planeint
 public :: cut3d_lineint
 public :: cut3d_pointint
 public :: cut3d_wffile

CONTAINS  !===========================================================
!!***

!!****f* m_cut3d/cut3d_hirsh
!! NAME
!! cut3d_hirsh
!!
!! FUNCTION
!! Compute the Hirshfeld charges
!!
!! INPUTS
!!  grid_den(nrx,nry,nrz)= density on the grid
!!  natom = number of atoms in the unit cell
!!  nrx,nry,nrz= number of points in the grid for the three directions
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  xcart(3,natom) = different positions of the atoms in the unit cell
!!  zion=(ntypat)gives the ionic charge for each type of atom
!!  znucl(ntypat)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  write the Hirshfeld charge decomposition
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_hirsh(grid_den,natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrx,nry,nrz,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: grid_den(nrx,nry,nrz),rprimd(3,3),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(in) :: xcart(3,natom)

!Local variables -------------------------
!scalars
 integer,parameter :: prtcharge1=1
 integer :: ierr,ipoint,itypat,mpoint,temp_unit
 real(dp) :: minimal_den
 real(dp) :: param1,param2,xx,yy
 character(len=fnlen) :: file_allelectron
 character(len=500) :: msg
!arrays
 integer,allocatable :: npoint(:)
 real(dp),allocatable :: aeden(:,:),hcharge(:),hden(:),hweight(:),radii(:,:)

! *********************************************************************

!1. Read the 1D all-electron atomic files
!Store the radii in radii(:,itypat), and the all-electron
!densities in aeden(:,itypat). The number of the last
!point with significant density is stored in npoint(itypat)

 minimal_den=tol6
 mpoint=4000
 ABI_ALLOCATE(npoint,(ntypat))
 ABI_ALLOCATE(radii,(4000,ntypat))
 ABI_ALLOCATE(aeden,(4000,ntypat))
 do itypat=1,ntypat
   write(std_out,'(a)' )' Please, give the filename of the all-electron density file'
   write(std_out,'(a,es16.6)' )' for the first type of atom, with atomic number=',znucl(itypat)
   if (read_string(file_allelectron, unit=std_in) /= 0) then
     MSG_ERROR("Fatal error!")
   end if
   write(std_out,*)' The name you entered is : ',trim(file_allelectron),ch10
   ierr = open_file(file_allelectron,msg,newunit=temp_unit,form='formatted',status='old')
   if (ierr/=0) then
     MSG_ERROR(msg)
   else
     read(temp_unit, *) param1, param2
     do ipoint=1,mpoint
!      Either the file is finished
       read(temp_unit, *, end=888) xx,yy
       radii(ipoint,itypat)=xx
       aeden(ipoint,itypat)=yy
!      Or the density is lower than the minimal significant value
       if(yy<minimal_den)exit
     end do
     888 continue
     npoint(itypat)=ipoint-1
     if(ipoint==mpoint)then
       write(std_out,*)' hirsh : mpoint is too low, increase its value to match ipoint.'
     end if
   end if
   close(temp_unit)
 end do

 ABI_MALLOC(hden,(natom))
 ABI_MALLOC(hcharge,(natom))
 ABI_MALLOC(hweight,(natom))

 call dens_hirsh(mpoint,radii,aeden,npoint,minimal_den,grid_den, &
  natom,nrx,nry,nrz,ntypat,rprimd,xcart,typat,zion,prtcharge1,hcharge,hden,hweight)

 ABI_FREE(hweight)
 ABI_FREE(aeden)
 ABI_FREE(hcharge)
 ABI_FREE(hden)
 ABI_FREE(npoint)
 ABI_FREE(radii)

end subroutine cut3d_hirsh
!!***

!!****f* m_cut3d/cut3d_lineint
!! NAME
!! cut3d_lineint
!!
!! FUNCTION
!! Computes the values along a line defined by two points
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

 subroutine cut3d_lineint(gridtt,gridux,griddy,gridmz,nr1,nr2,nr3,nspden,rprimd)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: griddy(nr1,nr2,nr3)
 real(dp),intent(in) :: gridmz(nr1,nr2,nr3),gridtt(nr1,nr2,nr3)
 real(dp),intent(in) :: gridux(nr1,nr2,nr3),rprimd(3,3)

!Local variables--------------------------------------------------------
!scalars
 integer :: inpopt,inpopt2,k2,nresol,okline,unt
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux,dx,dy,dz,length
 character(len=fnlen) :: filnam
 character(len=500) :: msg
!arrays
 real(dp) :: cent(3),r1(3),r2(3),rcart(3),rr(3),x1(3),x2(3)

! *********************************************************************

 okline=0
 do while (okline==0)
   write(std_out,*) ' Type 1) for a line between two cartesian-defined points'
   write(std_out,*) '   or 2) for a line between two crystallographic-defined points '
   write(std_out,*) '   or 3) for a line defined by its direction in cartesion coordinates'
   write(std_out,*) '   or 4) for a line defined by its direction in crystallographic coordinates'
   read(std_in,*) inpopt
   write(std_out,*) ' You typed ',inpopt,ch10
   if (inpopt==1 .or. inpopt ==2 .or. inpopt==3 .or. inpopt==4) okline=1
 end do

!In the case of a line defined by its two extreme points
 if (inpopt==1) then
   write(std_out,*) ' Type the first point coordinates (Bohrs):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(std_in,*) x1
   write(std_out,'(a,3es16.6,a)') ' You typed ',x1,ch10
   call reduce(r1,x1,rprimd)

   write(std_out,*) ' Type the second point coordinates (Bohrs):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(std_in,*) x2
   write(std_out,'(a,3es16.6,a)') ' You typed ',x2,ch10
   call reduce(r2,x2,rprimd)
 end if

 if (inpopt==2) then
   write(std_out,*) ' Type the first point coordinates (fractional):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(std_in,*) r1
   write(std_out,'(a,3es16.6,a)') ' You typed ',r1,ch10

   write(std_out,*) ' Type the second point coordinates (fractional):'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(std_in,*) r2
   write(std_out,'(a,3es16.6,a)') ' You typed ',r2,ch10
 end if

 if(inpopt==3 .or. inpopt==4 )then

   write(std_out,*) 'Please enter now the line direction:'
   write(std_out,*) '    -> X-dir   Y-dir   Z-dir:'
   read(std_in,*) x2
   write(std_out,'(a,3es16.6,a)') 'The line direction is:',x2(1),x2(2),x2(3),ch10

   if (inpopt == 4) then
     rcart=matmul(x2,rprimd)
     x2(:)=rcart(:)
     write(std_out,'(a,3es16.6,a)') 'Expressed in cartesian coordinates: ',x2(1),x2(2),x2(3),ch10
   end if

   call normalize(x2)

   write(std_out,*) 'Enter now the central point of line:'
   write(std_out,*) 'Type 1) for cartesian coordinates'
   write(std_out,*) '  or 2) for crystallographic coordinates'
   read(std_in,*) inpopt2
   if (inpopt2==1 .or. inpopt2==2) then
     write(std_out,*) 'Type the point coordinates:'
     write(std_out,*) '    -> X-Coord   Y-Coord   Z-Coord:'
     read(std_in,*) cent
     write(std_out,'(a,3es16.6,a)') 'Central point coordinates:', cent(1),cent(2),cent(3),ch10
     if (inpopt2==2) then
       rcart=matmul(cent,rprimd)
       cent(:)=rcart(:)
       write(std_out,'(a,3es16.6,a)') 'Expressed in cartesian coordinates:',cent(1),cent(2),cent(3),ch10
     end if
     write(std_out,*) 'Enter line length (in cartesian coordinates, in Bohr):'
     read(std_in,*) length

!    Compute the extremal points in cartesian coordinates
     x1(:)=cent(:)-length*x2(:)*half
     x2(:)=cent(:)+length*x2(:)*half

!    Transfer to crystallographic coordinates
     call reduce(r1,x1,rprimd)
     call reduce(r2,x2,rprimd)

   end if

 end if ! inpopt

 write(std_out,*)
 write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates of the first point  :',r1
 write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates of the second point :',r2
 write(std_out,*)

 write(std_out,*) '  Enter line resolution:   (integer, number of points on the line)'
 read(std_in,*) nresol
 write(std_out,*) ' You typed',nresol,ch10

!At this moment the code knows everything about the geometric input, the data and
!the line direction. It will further calculate the values along this line using
!an interpolation

 write(std_out,*) ch10,'  Enter the name of an output file:'
 if (read_string(filnam, unit=std_in) /= 0) then
   MSG_ERROR("Fatal error!")
 end if
 write(std_out,*) '  The name of your file is : ',trim(filnam),ch10

 if (open_file(filnam,msg,newunit=unt,status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if

 dx=(r2(1)-r1(1))/nresol
 dy=(r2(2)-r1(2))/nresol
 dz=(r2(3)-r1(3))/nresol

!DEBUG
!write(std_out,*)' nspden=',nspden
!ENDDEBUG

 if(nspden==1)then
   write(std_out,*)' Index of point   value '
 else if (nspden==2)then
   write(std_out,*)' Index of point   non-spin-polarized   spin up       spin down     difference '
 else if (nspden==4)then
   write(std_out,*)' Index of point   non-spin-polarized      x              y              z '
 end if

 do k2=0,nresol

   rr(1)=r1(1)+k2*dx
   rr(2)=r1(2)+k2*dy
   rr(3)=r1(3)+k2*dz

   rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
   rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
   rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

   denvaltt = interpol3d(rr,nr1,nr2,nr3,gridtt)
   if(nspden==1)then
     write(unt, '(i13,es22.12)' ) k2,denvaltt
     write(std_out,'(i13,es22.12)' ) k2,denvaltt

   else if(nspden==2 .or. nspden==4)then
     denvalux = interpol3d(rr,nr1,nr2,nr3,gridux)
     denvaldy = interpol3d(rr,nr1,nr2,nr3,griddy)
     denvalmz = interpol3d(rr,nr1,nr2,nr3,gridmz)
     write(unt, '(i13,4(es22.12))' ) k2,denvaltt,denvalux,denvaldy,denvalmz
     write(std_out,'(i13,4es22.12)' ) k2,denvaltt,denvalux,denvaldy,denvalmz
   end if
 end do

 close(unt)

end subroutine cut3d_lineint
!!***

!!****f* m_cut3d/normalize
!! NAME
!! normalize
!!
!! FUNCTION
!! Normalizes the value of v
!!
!! INPUTS
!!  v = on entry, vector to be normalized
!!
!! OUTPUT
!!  v = on exit, vector normalized

!! SIDE EFFECTS
!!   v=value to be normalized
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine normalize(v)

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(inout) :: v(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: idir
 real(dp) :: norm

! *************************************************************************

 norm=0.0
 do idir=1,3
   norm=norm+v(idir)**2
 end do
 norm=sqrt(norm)

 do idir=1,3
   v(idir)=v(idir)/norm
 end do

end subroutine normalize
!!***

!!****f* m_cut3d/cut3d_planeint
!! NAME
!! cut3d_planeint
!!
!! FUNCTION
!! Computes the values within a plane
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! natom=integer number of atoms
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!! tau(3,nat)=atomic positions in 3D cartesian space (from XMOL format)
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_planeint(gridtt,gridux,griddy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,tau)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: griddy(nr1,nr2,nr3)
 real(dp),intent(in) :: gridmz(nr1,nr2,nr3),gridtt(nr1,nr2,nr3)
 real(dp),intent(in) :: gridux(nr1,nr2,nr3),rprimd(3,3),tau(3,natom)

!Local variables -------------------------
!scalars
 integer :: iat,idir,ii,inpopt,itypat,k2,k3,mu,nresoll,nresolw,okhkl,okinp
 integer :: okparam,oksure,unt
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux,length
 real(dp) :: width,xcoord,ycoord
 character(len=fnlen) :: filnam
 character(len=500) :: msg
!arrays
 integer :: hkl(3)
 real(dp) :: cent(3),mminv(3,3),r1(3),r2(3),r3(3),rcart(3),rr(3),x1(3)
 real(dp) :: x2(3),x3(3),xcart(3)

! *********************************************************************

!Several lines to compute the transformation matrix from crystallographic to cartesian

 call matr3inv(rprimd,mminv)

!Start of the real input of the plane orientation

 okinp=0
 do while (okinp==0)
   write(std_out,*)
   write(std_out,*) '  Type 1) for a plane passing through 3 atoms'
   write(std_out,*) '    or 2) for a plane passing through 3 cartesian points'
   write(std_out,*) '    or 3) for a plane passing through 3 crystallographic points'
   write(std_out,*) '    or 4) for a plane parallel to a crystallographic plane'
   write(std_out,*) '    or 5) for a plane orthogonal to a cartesian direction'
   write(std_out,*) '    or 6) for a plane orthogonal to a crystallographic direction'
   write(std_out,*) '    or 0) to stop'
   read(std_in,*) itypat
   select case (itypat)

   case (0)
     stop

!      A plane passing through 3 atoms
   case (1)
     write(std_out,*) '  The X axis will be through atms: 1,2 '
     write(std_out,*) '  Define each atom by its species and its number:'
     write(std_out,*) '    -> atom 1 (iat):'
     read(std_in,*) iat
     x1(1)=tau(1,iat)
     x1(2)=tau(2,iat)
     x1(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x1
     write(std_out,*)
     write(std_out,*) '    -> atom 2 (iat):'
     read(std_in,*) iat
     x2(1)=tau(1,iat)
     x2(2)=tau(2,iat)
     x2(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x2
     write(std_out,*)
     write(std_out,*) '    -> atom 3 (iat):'
     read(std_in,*) iat
     x3(1)=tau(1,iat)
     x3(2)=tau(2,iat)
     x3(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x3
     write(std_out,*)

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

!      A plane passing through 3 cartesian points
   case (2)
     write(std_out,*) '  The X axis will be through points: 1,2 '
     write(std_out,*) '  Define each :point coordinates'
     write(std_out,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x1(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x1
     write(std_out,*)
     write(std_out,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x2(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x2
     write(std_out,*)
     write(std_out,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x3(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x3
     write(std_out,*)

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

!      A plane passing through 3 crystallographic points
   case (3)
     write(std_out,*) '  The X axis will be through points: 1,2 '
     write(std_out,*) '  Define each :point coordinates'
     write(std_out,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r1
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r1
     write(std_out,*)
     write(std_out,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r2
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r2
     write(std_out,*)
     write(std_out,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r3
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r3
     write(std_out,*)

!      Transforms the points coordinates into cartesian
     do mu=1,3
       x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
       x2(mu)=rprimd(mu,1)*r2(1)+rprimd(mu,2)*r2(2)+rprimd(mu,3)*r2(3)
       x3(mu)=rprimd(mu,1)*r3(1)+rprimd(mu,2)*r3(2)+rprimd(mu,3)*r3(3)
     end do

     write(std_out,*) ' Cartesian positions:'
     write(std_out,*) x1
     write(std_out,*) x2
     write(std_out,*) x3

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

!      A plane parallel to a crystallographic plane
   case (4)
     okhkl=0
     do while (okhkl==0)
       write(std_out,*) '  Enter plane coordinates:'
       write(std_out,*) '    -> H  K  L '
       read(std_in,*) hkl
       if (.not. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3)==0)) okhkl=1
     end do
     write(std_out,*) ' Miller indices are:',hkl

     do ii=1,3
       x1(ii)=mminv(ii,1)*hkl(1) + mminv(ii,2)*hkl(2) + mminv(ii,3)*hkl(3)
     end do
     write(std_out,*) ' Orthogonal vector to the plane',x1

     call normalize(x1)
     if((x1(1).ne.0).or.(x1(2).ne.0)) then
       x2(1)=-x1(2)
       x2(2)=x1(1)
       x2(3)=0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

!      A plane orthogonal to a cartesian direction
   case (5)
     write(std_out,*) '  Enter the cartesian coordinates of the vector orthogonal to plane:'
     write(std_out,*) '    -> X-dir   Y-dir   Z-dir (Angstroms or Bohrs):'
     read(std_in,*) x1
     call normalize(x1)
     if((x1(1).ne.0).or.(x1(2).ne.0)) then
       x2(1)=-x1(2)
       x2(2)=x1(1)
       x2(3)=0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

!      A plane orthogonal to a crystallographic direction
   case (6)
     write(std_out,*) '  Enter crystallographic vector orthogonal to plane:'
     write(std_out,*) '    -> X-dir   Y-dir   Z-dir (Fractional coordinates):'
     read(std_in,*) r1
     okinp=1
     do mu=1,3
       x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
     end do
     call normalize(x1)
     if((x1(1).ne.0).or.(x1(2).ne.0)) then
       x2(1)=-x1(2)
       x2(2)=x1(1)
       x2(3)=0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x2,x1,x3)
     call normalize(x3)
     okinp=1

   case default
     okinp=0
     write(std_out,*) 'Input option do not correspond to the available options'
     write(std_out,*) 'Please try again'
   end select

 end do

!At this moment the family of planes was defined
!The code knows also some of the geometric input
!It will proceed to the anchorage of the plane onto a point and then
!to the effective calculation

 write(std_out,*) '  Vectors: (orthogonal & normalized)   '
 write(std_out,'(11x,a,3f10.6)') ' X-dir in the plot         ',x2
 write(std_out,'(11x,a,3f10.6)') ' Y-dir in the plot         ',x3
 write(std_out,'(11x,a,3f10.6)') ' Z-dir (orth. to the plot) ',x1

 write(std_out,*)
 write(std_out,*) '  Enter central point of plane (Bohrs):'
 write(std_out,*) '  Type 1) for Cartesian coordinates.'
 write(std_out,*) '    or 2) for Crystallographic coordinates.'
 read(std_in,*) inpopt
 write(std_out,*) '    -> X-Coord   Y-Coord   Z-Coord:'
 read(std_in,*) cent

 if (inpopt==2) then

   do mu=1,3
     rcart(mu)=rprimd(mu,1)*cent(1)+rprimd(mu,2)*cent(2)+rprimd(mu,3)*cent(3)
   end do

   cent(:)=rcart(:)
   write(std_out,'(a,3f16.6)' ) ' Expressed in cartesian coordinates: ',cent(1),cent(2),cent(3)

 end if

 okparam=0
 do while(okparam==0)
   write(std_out,*)
   write(std_out,*) '  Enter plane width:'
   read(std_in,*) width
   write(std_out,*) '  Enter plane length:'
   read(std_in,*) length
   write(std_out,*)
   write(std_out,*) '  Enter plane resolution in width:'
   read(std_in,*) nresolw
   write(std_out,*) '  Enter plane resolution in length:'
   read(std_in,*) nresoll
   write(std_out,*) ch10,'  Enter the name of an output file:'
   if (read_string(filnam, unit=std_in) /= 0) then
     MSG_ERROR("Fatal error!")
   end if
   write(std_out,*) '  The name of your file is : ',trim(filnam)
   write(std_out,*)
   write(std_out,*) '  You asked for a plane of ',length,' x ',width
   write(std_out,*) '  With a resolution of ',nresoll,' x ',nresolw
   write(std_out,*) '  The result will be redirected to the file:  ',trim(filnam)
   write(std_out,*) '  These parameters may still be changed.'
   write(std_out,*) '  Are you sure you want to keep them? (1=default=yes,2=no) '
   read(std_in,*) oksure
   if (oksure/=2) okparam=1
 end do

 if (open_file(filnam,msg,newunit=unt,status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if

 do k2=-nresoll/2,nresoll/2
   do k3=-nresolw/2,nresolw/2
     rcart(1)=cent(1) + k2*x2(1)*length/nresoll + k3*x3(1)*width/nresolw
     rcart(2)=cent(2) + k2*x2(2)*length/nresoll + k3*x3(2)*width/nresolw
     rcart(3)=cent(3) + k2*x2(3)*length/nresoll + k3*x3(3)*width/nresolw
     xcoord=k2*length/nresoll
     ycoord=k3*width/nresolw
     call reduce(rr,rcart,rprimd)
     rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
     rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
     rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)
     denvaltt = interpol3d(rr,nr1,nr2,nr3,gridtt)
     if(nspden==2 .or. nspden==4)then
       denvalux = interpol3d(rr,nr1,nr2,nr3,gridux)
       denvaldy = interpol3d(rr,nr1,nr2,nr3,griddy)
       denvalmz = interpol3d(rr,nr1,nr2,nr3,gridmz)
     end if
     if(nspden==1)then
       write(unt, '(3e16.8)' ) xcoord,ycoord,denvaltt
     else
       write(unt, '(3e16.8)' ) xcoord,ycoord,denvaltt,denvalux,denvaldy,denvalmz
     end if
   end do
 end do

 close(unt)

 end subroutine cut3d_planeint
!!***

!!****f* m_cut3d/cut3d_pointint
!! NAME
!! cut3d_pointint
!!
!! FUNCTION
!! Computes the values at any point rr (this point is input from keyboard)
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!!   only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_pointint(gridt,gridu,gridd,gridm,nr1,nr2,nr3,nspden,rprimd)

!Arguments--------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: gridd(nr1,nr2,nr3),gridm(nr1,nr2,nr3)
 real(dp),intent(in) :: gridt(nr1,nr2,nr3),gridu(nr1,nr2,nr3),rprimd(3,3)

!Local variables--------------------------------------------------------
!scalars
 integer :: inpopt,mu,okinp
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux
!arrays
 real(dp) :: rcart(3),rr(3)

! *************************************************************************

 okinp=0
 do while (okinp==0)
   write(std_out,*) ' Select the coordinate system:'
   write(std_out,*) ' Type 1) for cartesian coordinates'
   write(std_out,*) '  or 2) for crystallographic coordinates'
   read(std_in,*) inpopt
   if (inpopt==1 .or. inpopt==2) okinp=1
 end do

 if (inpopt==1) then

   write(std_out,*) ' Input point Cartesian Coord:  X  Y  Z'
   read(std_in,*) rcart(1),rcart(2),rcart(3)
   call reduce(rr,rcart,rprimd)
   write(std_out,'(a,3es16.6)' ) ' Crystallographic coordinates: ',rr(1:3)

 else

   write(std_out,*) ' Input point Crystallographic Coord:  X  Y  Z'
   read(std_in,*) rr(1),rr(2),rr(3)

   do mu=1,3
     rcart(mu)=rprimd(mu,1)*rr(1)+rprimd(mu,2)*rr(2)+rprimd(mu,3)*rr(3)
   end do

   write(std_out,*) ' Cartesian coordinates : '
   write(std_out,'(3es16.6)' ) rcart(1),rcart(2),rcart(3)

 end if

!At this moment the code knows everything needed about the geometric input
!It will further proceed to calculate the interpolation at the demanded point

 rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
 rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
 rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

 write(std_out,'(a,es16.6)' ) ' X coordinate, r1 is:',rr(1)
 write(std_out,'(a,es16.6)' ) ' Y coordinate, r2 is:',rr(2)
 write(std_out,'(a,es16.6)' ) ' Z coordinate, r3 is:',rr(3)

!devalt = total density value
!devalu = spin-up density value
!devald = spin-down density value
!devalm = magnetization density value
 denvaltt = interpol3d(rr,nr1,nr2,nr3,gridt)
 if(nspden==2 .or. nspden==4)then
   denvalux = interpol3d(rr,nr1,nr2,nr3,gridu)
   denvaldy = interpol3d(rr,nr1,nr2,nr3,gridd)
   denvalmz = interpol3d(rr,nr1,nr2,nr3,gridm)
 end if
 write(std_out,*)
 write(std_out,*)'---------------------------------------------'
 write(std_out,'(a,es16.6)') ' Non-spin-polarized value= ',denvaltt
 if(nspden==2)then
   write(std_out,'(a,es16.6)')' Spin-up value           = ',denvalux
   write(std_out,'(a,es16.6)')' Spin-down value         = ',denvaldy
   write(std_out,'(a,es16.6)')' Spin difference value   = ',denvalmz
 else if(nspden==4)then
   write(std_out,'(a,es16.6)')' x component             = ',denvalux
   write(std_out,'(a,es16.6)')' y component             = ',denvaldy
   write(std_out,'(a,es16.6)')' z component             = ',denvalmz
 end if
 write(std_out,*)'---------------------------------------------'

end subroutine cut3d_pointint
!!***

!!****f* m_cut3d/reduce
!! NAME
!! reduce
!!
!! FUNCTION
!! Transforms coordinates of an input point
!! from cartesian to crystallographic
!!
!! INPUTS
!! rcart(3)=position vector in crystallographic coordinates
!! rprimd(3,3)=orientation of the unit cell in 3D
!!
!! OUTPUT
!! r(3)=position vector in cartesian coordinates
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine reduce(r,rcart,rprimd)

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(in) :: rcart(3),rprimd(3,3)
 real(dp),intent(out) :: r(3)

!Local variables--------------------------------------------------------
!scalars
!arrays
 real(dp) :: mminv(3,3)

! *************************************************************************

 call matr3inv(rprimd,mminv)
 r(1)=rcart(1)*mminv(1,1)+rcart(2)*mminv(2,1)+rcart(3)*mminv(3,1)
 r(2)=rcart(1)*mminv(1,2)+rcart(2)*mminv(2,2)+rcart(3)*mminv(3,2)
 r(3)=rcart(1)*mminv(1,3)+rcart(2)*mminv(2,3)+rcart(3)*mminv(3,3)

end subroutine reduce
!!***

!!****f* m_cut3d/cut3d_rrho
!! NAME
!! cut3d_rrho
!!
!! FUNCTION
!! Reads in the charge in mkdens3D format
!! The file was opened in the calling program, unit number 19.
!! The header was already read in the case of the unformatted file
!!
!! INPUTS
!! path=File name
!! varname=Name of the netcdf variable to be read.
!! iomode=flag specifying the IO library.
!! nr1=grid_full size along x
!! nr2=grid_full size along y
!! nr3=grid_full size along z
!! nspden=number of spin polartized densities (1 for non-spin polarized, 2 for spin-polarized)
!!
!! OUTPUT
!! grid_full(nr1,nr2,nr3)=grid_full matrix
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_rrho(path,varname,iomode,grid_full,nr1,nr2,nr3,nspden)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: iomode,nr1,nr2,nr3,nspden
 character(len=*),intent(in) :: path,varname
!arrays
 real(dp),intent(out),target :: grid_full(nr1,nr2,nr3,nspden)

!Local variables--------------------------------------------------------
!scalars
 integer :: ispden,unt,fform
#ifdef HAVE_NETCDF
 integer :: varid
#endif
 character(len=500) :: msg
 type(hdr_type) :: hdr

! *************************************************************************

 select case (iomode)
 case (IO_MODE_FORTRAN)
   !Unformatted, on one record
   if (open_file(path, msg, newunit=unt, form='unformatted', status='old', action="read") /= 0) then
     MSG_ERROR(msg)
   end if
   call hdr_fort_read(hdr, unt, fform)
   ABI_CHECK(fform /= 0, sjoin("Error while reading:", path))
   call hdr%free()

   do ispden=1,nspden
     read(unit=unt) grid_full(1:nr1,1:nr2,1:nr3,ispden)
   end do

   close(unt)

 case (IO_MODE_ETSF)
   ! ETSF case
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_read(unt, path, xmpi_comm_self))
   NCF_CHECK(nf90_inq_varid(unt, varname, varid))
   ! [cplex, n1, n2, n3, nspden]
   NCF_CHECK(nf90_get_var(unt, varid, grid_full, start=[1,1,1,1,1], count=[1, nr1,nr2,nr3,nspden]))
   NCF_CHECK(nf90_close(unt))
#else
   MSG_ERROR('netcdf support is not compiled. Reconfigure with --enable-netcdf.')
#endif

 case default
   MSG_BUG(sjoin("invalid iomode:", itoa(iomode)))
 end select

end subroutine cut3d_rrho
!!***

!!****f* m_cut3d/vdot
!! NAME
!! vdot
!!
!! FUNCTION
!! Computes the cross product of two vectors
!!
!! INPUTS
!! x1(3)=first vector
!! x2(3)=second vector
!!
!! OUTPUT
!! x3(3)=cross product of x1 * x2
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine vdot(x1,x2,x3)

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(in) :: x1(3),x2(3)
 real(dp),intent(out) :: x3(3)

!Local variables-------------------------------

! *************************************************************************

 x3(1)=x1(2)*x2(3)-x2(2)*x1(3)
 x3(2)=x1(3)*x2(1)-x2(3)*x1(1)
 x3(3)=x1(1)*x2(2)-x2(1)*x1(2)

end subroutine vdot
!!***

!!****f* m_cut3d/cut3d_volumeint
!! NAME
!! cut3d_volumeint
!!
!! FUNCTION
!! Computes the values within a volume
!!
!! INPUTS
!! gridtt(nr1,nr2,nr3)=Total density
!! gridux(nr1,nr2,nr3)=spin-Up density, or magnetization density in X direction
!! griddy(nr1,nr2,nr3)=spin-Down density, or magnetization density in Y direction
!! gridmz(nr1,nr2,nr3)=spin-polarization density or magnetization density in Z direction
!! natom=integer number of atoms
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!! tau(3,natom)=list of atoms in cartesian coordinates
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_volumeint(gridtt,gridux,griddy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,tau)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: griddy(nr1,nr2,nr3)
 real(dp),intent(in) :: gridmz(nr1,nr2,nr3),gridtt(nr1,nr2,nr3)
 real(dp),intent(in) :: gridux(nr1,nr2,nr3),rprimd(3,3),tau(3,natom)

!Local variables -------------------------
!scalars
 integer :: fileformattype,iat,idir,ii,inpopt,itypat,k1,k2,k3,mu,nresolh
 integer :: nresoll,nresolw,okhkl,okparam,planetype
 integer :: referenceposition,unt
 real(dp) :: denvaldy,denvalmz,denvaltt,denvalux,height
 real(dp) :: length,width
 real(dp) :: xm,xp,ym,yp,zm,zp
 character(len=fnlen) :: filnam
 character(len=500) :: msg
!arrays
 integer :: hkl(3)
 real(dp) :: cent(3),centpl(3),mminv(3,3),r1(3),r2(3),r3(3),rcart(3)
 real(dp) :: rr(3),x1(3),x2(3),x3(3),xcart(3)
 real(dp),allocatable :: rhomacudy(:,:),rhomacumz(:,:),rhomacutt(:,:)
 real(dp),allocatable :: rhomacuux(:,:)

! *********************************************************************

 call matr3inv(rprimd,mminv)
!Start of the real input of the volume orientation

 write(std_out,*)
 write(std_out,*) ' The volume is an orthogonal prism, that is defined by: '
 write(std_out,*) ' the basal plane and'
 write(std_out,*) ' the height perpendicular to the basal plane'
 write(std_out,*)
 write(std_out,*) ' First you will define the basal plane '
 write(std_out,*) ' second you will define the height'
 write(std_out,*) ' and third you will define the basal plane position '
 write(std_out,*) ' along the height vector'

 do
   write(std_out,*)
   write(std_out,*) '  Type 1) for a plane passing through 3 atoms'
   write(std_out,*) '    or 2) for a plane passing through 3 cartesian points'
   write(std_out,*) '    or 3) for a plane passing through 3 crystallographic points'
   write(std_out,*) '    or 4) for a plane parallel to a crystallographic plane'
   write(std_out,*) '    or 5) for a plane orthogonal to a cartesian direction'
   write(std_out,*) '    or 6) for a plane orthogonal to a crystallographic direction'
   write(std_out,*) '    or 0) to stop'
   read(std_in,*) itypat

   select case (itypat)

   case (0)
     stop

!      A plane passing through 3 atoms
   case (1)
     write(std_out,*) '  The X axis will be through atoms: 1,2 '
     write(std_out,*) '  Define each atom by its species and its number:'
     write(std_out,*) '    -> atom 1 (iat):'
     read(std_in,*) iat
     x1(1)=tau(1,iat)
     x1(2)=tau(2,iat)
     x1(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x1
     write(std_out,*)
     write(std_out,*) '    -> atom 2 (iat):'
     read(std_in,*) iat
     x2(1)=tau(1,iat)
     x2(2)=tau(2,iat)
     x2(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x2
     write(std_out,*)
     write(std_out,*) '    -> atom 3 (iat):'
     read(std_in,*) iat
     x3(1)=tau(1,iat)
     x3(2)=tau(2,iat)
     x3(3)=tau(3,iat)
     write(std_out,'(a,3f10.6)') '        position: ',x3
     write(std_out,*)

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     exit

!      A plane passing through 3 cartesian points
   case (2)
     write(std_out,*) '  The X axis will be through points: 1,2 '
     write(std_out,*) '  Define each :point coordinates'
     write(std_out,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x1(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x1
     write(std_out,*)
     write(std_out,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x2(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x2
     write(std_out,*)
     write(std_out,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) xcart
     x3(:)=xcart(:)
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',x3
     write(std_out,*)

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     exit

!      A plane passing through 3 crystallographic points
   case (3)
     write(std_out,*) '  The X axis will be through points: 1,2 '
     write(std_out,*) '  Define each :point coordinates'
     write(std_out,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r1
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r1
     write(std_out,*)
     write(std_out,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r2
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r2
     write(std_out,*)
     write(std_out,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
     read(std_in,*) r3
     write(std_out,'(a,3f10.6)') ' crystallographic position: ',r3
     write(std_out,*)

!      Transforms the points coordinates into cartesian
     do mu=1,3
       x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
       x2(mu)=rprimd(mu,1)*r2(1)+rprimd(mu,2)*r2(2)+rprimd(mu,3)*r2(3)
       x3(mu)=rprimd(mu,1)*r3(1)+rprimd(mu,2)*r3(2)+rprimd(mu,3)*r3(3)
     end do

     write(std_out,*) ' Cartesian positions:'
     write(std_out,*) x1
     write(std_out,*) x2
     write(std_out,*) x3

!      Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
     do idir=1,3
       x2(idir)=x2(idir)-x1(idir)
       x3(idir)=x3(idir)-x1(idir)
     end do
     call normalize(x2)
     call vdot(x3,x2,x1)
     call normalize(x1)
     call vdot(x2,x1,x3)
     call normalize(x3)
     exit

!      A plane parallel to a crystallographic plane
   case (4)
     okhkl=0
     do while (okhkl==0)
       write(std_out,*) '  Enter plane coordinates:'
       write(std_out,*) '    ->H  K  L '
       read(std_in,*) hkl
       if (.not. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3)==0)) okhkl=1
     end do
     write(std_out,*) ' Miller indices are:',hkl

     do ii=1,3
       x1(ii)=mminv(ii,1)*hkl(1) + mminv(ii,2)*hkl(2) + mminv(ii,3)*hkl(3)
     end do
     write(std_out,*) ' Orthogonal vector to the plane',x1

     call normalize(x1)
     if((x1(1).ne.0).or.(x1(2).ne.0)) then
       x2(1)=-x1(2)
       x2(2)=x1(1)
       x2(3)=0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x2,x1,x3)
     call normalize(x3)
     exit

!      A plane orthogonal to a cartesian direction
   case (5)
     write(std_out,*) '  Enter cartesian vector orthogonal to plane:'
     write(std_out,*) '    -> X-dir   Y-dir   Z-dir (Angstroms or Bohrs):'
     read(std_in,*) x1
     call normalize(x1)
     if((x1(1).ne.0).or.(x1(2).ne.0)) then
       x2(1)=-x1(2)
       x2(2)=x1(1)
       x2(3)=0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x1,x2,x3)
     call normalize(x3)
     exit

!      A plane orthogonal to a crystallographic direction
   case (6)
     write(std_out,*) '  Enter crystallographic vector orthogonal to plane:'
     write(std_out,*) '    -> X-dir   Y-dir   Z-dir (Fractional coordinates):'
     read(std_in,*) r1
     do mu=1,3
       x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
     end do
     call normalize(x1)
     if(abs(x1(1))<tol10 .or. abs(x1(2)) < tol10) then
       x2(1)=-x1(2)
       x2(2)= x1(1)
       x2(3)= 0
       call normalize(x2)
     else
       x2(1)=1
       x2(2)=0
       x2(3)=0
     end if
     call vdot(x1,x2,x3)
     call normalize(x3)
     exit

   case default
     write(std_out,*) ' Input option does not correspond to one available option'
     write(std_out,*) ' Please try again'
     cycle

   end select
 end do

!At this moment the family of planes was defined
!The code knows also some of the geometric input
!It will proceed to the anchorage of the plane onto a point and then
!to the effective calculation

 write(std_out,*) '  Vectors: (orthogonal & normalized)   '
 write(std_out,'(11x,a,3f10.6)') '  X-dir in the plot         ',x2
 write(std_out,'(11x,a,3f10.6)') '  Y-dir in the plot         ',x3
 write(std_out,'(11x,a,3f10.6)') '  Z-dir (orth. to the plot) ',x1

 do
   write(std_out,*)
   write(std_out,*) '  Enter reference point of plane (Bohr):'
   write(std_out,*) '  Type 1) for Cartesian coordinates.'
   write(std_out,*) '    or 2) for Crystallographic coordinates.'
   read(std_in,*) inpopt

   select case (inpopt)

   case(1)
     write(std_out,*) '    -> X-Coord   Y-Coord   Z-Coord:'
     read(std_in,*) cent
     exit
   case(2)
     write(std_out,*) '    -> X-Coord   Y-Coord   Z-Coord:'
     read(std_in,*) cent
     do mu=1,3
       rcart(mu)=rprimd(mu,1)*cent(1)+rprimd(mu,2)*cent(2)+rprimd(mu,3)*cent(3)
     end do
     cent(:)=rcart(:)
     write(std_out,'(a,3es16.6)' ) ' Expressed in cartesian coordinates: ',cent(1:3)
     exit
   case (3)
     cycle

   end select
 end do

!End of basal plane orientation

!Input box dimensions now

 write(std_out,*)
 write(std_out,*) ' It is now time to input the 3D box dimensions.'
 write(std_out,*) ' and the position of the basal plane in the box.'

 do
   write(std_out,*)
   write(std_out,*) '  Enter in-plane width:'
   read(std_in,*) width
   write(std_out,*) '  Enter in-plane length:'
   read(std_in,*) length
   write(std_out,*) '  Enter box height:'
   read(std_in,*) height
   write(std_out,*)
   write(std_out,*) ' Enter the position of the basal plane in the box:'
   do
     write(std_out,*)
     write(std_out,*) ' Type 1) for DOWN'
     write(std_out,*) ' Type 2) for MIDDLE'
     write(std_out,*) ' Type 3) for UP'
     read(std_in,*) planetype

     select case(planetype)

     case (1)
       exit
     case (2)
       exit
     case (3)
       exit
     case default
       cycle

     end select
   end do

   write(std_out,*) ' Enter the position of the reference point in the basal plane '

   do
     write(std_out,*)
     write(std_out,*) ' Type 1) for CENTRAL position '
     write(std_out,*) ' Type 2) for CORNER(0,0) position '
     read(std_in,*) referenceposition

     select case(referenceposition)

     case (1)
       exit
     case (2)
       exit
     case default
       cycle

     end select
   end do

   write(std_out,*)
   write(std_out,*) ' Enter the box grid values:'
   write(std_out,*) '  Enter plane resolution in width:'
   read(std_in,*) nresolw
   write(std_out,*) '  Enter plane resolution in lenth:'
   read(std_in,*) nresoll
   write(std_out,*) '  Enter height resolution:'
   read(std_in,*) nresolh
   write(std_out,*)
   write(std_out,*) ch10,'  Enter the name of an output file:'
   if (read_string(filnam, unit=std_in) /= 0) then
     MSG_ERROR("Fatal error!")
   end if
   write(std_out,*) '  The name of your file is : ',trim(filnam)

   do
     write(std_out,*) '  Enter the format of the output file:'
     write(std_out,*) '   Type 1=> ASCII formatted'
     write(std_out,*) '   Type 2=> 3D index + data, ASCII formatted'
     write(std_out,*) '   Type 3=> Molekel formatted'
     read(std_in,*) fileformattype
     if (fileformattype>=1 .and. fileformattype<=3) then
       exit
     else
       cycle
     end if
   end do

   write(std_out,*) ' You asked for a 3d box of:'
   write(std_out,*) length,' x ',width,' x ',height
   write(std_out,*) ' With a resolution of ;'
   write(std_out,*) nresoll,' x ',nresolw,' x ',nresolh
   write(std_out,*) ' The result will be redirected to the file:  ',trim(filnam)
   if      (fileformattype==1) then
     write(std_out,*) ' ASCII formatted'
   else if (fileformattype==2) then
     write(std_out,*) ' 3d index + data, ASCII formatted'
   else if (fileformattype==3) then
     write(std_out,*) ' Molekel formatted'
   end if
   write(std_out,*) ' These parameters may still be changed.'
   write(std_out,*) ' Are you sure you want to keep them? (1=default=yes,2=no) '
   read(std_in,*) okparam
   if (okparam==2) then
     cycle
   else
     exit
   end if
 end do

!Write the header of the Molekel input file

 if (fileformattype==1 .or. fileformattype==2) then
   if (open_file(filnam,msg,newunit=unt,status='unknown') /= 0) then
     MSG_ERROR(msg)
   end if
 else if (fileformattype==3) then
   if (open_file(filnam,msg,newunit=unt,form='unformatted') /= 0) then
     MSG_ERROR(msg)
   end if

   xm=0
   xp=length
   ym=0
   yp=width
   zm=0
   zp=height

   write(std_out,'(a)' )&
&   ' Extremas of the cube in which the system is placed (x,y,z), in Angs.:'
   write(std_out,'(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
   write(std_out,'(a,a,3i5)' ) ch10,&
&   ' Number of points per side:   ',nresolw+1,nresoll+1,nresolh+1
   write(std_out,'(a,a,i10,a,a)' ) ch10,&
&   ' Total number of points:  ',(nresolw+1)*(nresoll+1)*(nresolh+1),&
&   ch10,ch10

   write(unt) xm,xp,ym,yp,zm,zp,nresolw+1,nresoll+1,nresolh+1

 end if

!Allocate rhomacu in case of molekel output format
 ABI_ALLOCATE(rhomacutt,(nresoll+1,nresolw+1))
 ABI_ALLOCATE(rhomacuux,(nresoll+1,nresolw+1))
 ABI_ALLOCATE(rhomacudy,(nresoll+1,nresolw+1))
 ABI_ALLOCATE(rhomacumz,(nresoll+1,nresolw+1))

 do k1=0,nresolh

   select case (planetype)

!    Basal plane at the bottom
   case (1)
     centpl(1)=cent(1)+k1*x1(1)*height/nresolh
     centpl(2)=cent(2)+k1*x1(2)*height/nresolh
     centpl(3)=cent(3)+k1*x1(3)*height/nresolh

!      Basal plane in the middle
   case (2)
     centpl(1)=cent(1)+(k1-nresolh/2)*x1(1)*height/nresolh
     centpl(2)=cent(2)+(k1-nresolh/2)*x1(2)*height/nresolh
     centpl(3)=cent(3)+(k1-nresolh/2)*x1(3)*height/nresolh

!      Basal plane on the top
   case (3)
     centpl(1)=cent(1)+(k1-nresolh)*x1(1)*height/nresolh
     centpl(2)=cent(2)+(k1-nresolh)*x1(2)*height/nresolh
     centpl(3)=cent(3)+(k1-nresolh)*x1(3)*height/nresolh

   end select

   do k3=0,nresolw
     do k2=0,nresoll

       select case(referenceposition)

!        Reference point in the middle of the basal plane
       case (1)
         rcart(1)=centpl(1) + (k2-nresoll/2)*x2(1)*length/nresoll + (k3-nresolw/2)*x3(1)*width/nresolw
         rcart(2)=centpl(2) + (k2-nresoll/2)*x2(2)*length/nresoll + (k3-nresolw/2)*x3(2)*width/nresolw
         rcart(3)=centpl(3) + (k2-nresoll/2)*x2(3)*length/nresoll + (k3-nresolw/2)*x3(3)*width/nresolw

!          Reference point in the corner of the basal plane
       case (2)
         rcart(1)=centpl(1) + k2*x2(1)*length/nresoll + k3*x3(1)*width/nresolw
         rcart(2)=centpl(2) + k2*x2(2)*length/nresoll + k3*x3(2)*width/nresolw
         rcart(3)=centpl(3) + k2*x2(3)*length/nresoll + k3*x3(3)*width/nresolw

       end select

       call reduce(rr,rcart,rprimd)
       rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
       rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
       rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

       denvaltt = interpol3d(rr,nr1,nr2,nr3,gridtt)
       if(nspden==2 .or. nspden==4)then
         denvalux = interpol3d(rr,nr1,nr2,nr3,gridux)
         denvaldy = interpol3d(rr,nr1,nr2,nr3,griddy)
         denvalmz = interpol3d(rr,nr1,nr2,nr3,gridmz)
       end if

       if (fileformattype==1) then
         if(nspden==1)then
           write(unt, '(es22.12)' ) denvaltt
         else if(nspden==2 .or. nspden==4)then
           write(unt, '(4(es22.12))' ) denvaltt,denvalux,denvaldy,denvalmz
         end if

       else if (fileformattype==2) then
         if(nspden==1)then
           write(unt, '(4es22.12)' ) rcart, denvaltt
         else if(nspden==2 .or. nspden==4)then
           write(unt, '(3(e22.12), 4(es22.12))' ) rcart, denvaltt,denvalux,denvaldy,denvalmz
         end if

       else if (fileformattype==3) then
         rhomacutt(k2+1,k3+1)=denvaltt
         if(nspden==2 .or. nspden==4)then
           rhomacuux(k2+1,k3+1)=denvalux
           rhomacudy(k2+1,k3+1)=denvaldy
           rhomacumz(k2+1,k3+1)=denvalmz
         end if
       end if

     end do ! resoll
     write(unt, * )
   end do ! resolw

   if (fileformattype==3) then
     write(unt) rhomacutt(:,:)
     if(nspden==2 .or. nspden==4)then
       write(unt) rhomacuux(:,:)
       write(unt) rhomacudy(:,:)
       write(unt) rhomacumz(:,:)
     end if
   end if

 end do

 close(unt)

 ABI_DEALLOCATE(rhomacutt)
 ABI_DEALLOCATE(rhomacuux)
 ABI_DEALLOCATE(rhomacudy)
 ABI_DEALLOCATE(rhomacumz)

end subroutine cut3d_volumeint
!!***

!!****f* m_cut3d/cut3d_wffile
!! NAME
!! cut3d_wffile
!!
!! FUNCTION
!! Part of cut3d that gives the wavefunction for one kpt,one band
!! and one spin polarisation in real space.  The output depends on
!! the chosen option.
!!
!! INPUTS
!! wfk_fname=Name of the WFK file.
!! Needs an unformatted wave function from abinit.
!! ecut= effective ecut (ecut*dilatmx**2)
!! exchn2n3d= if 1, n2 and n3 are exchanged
!! istwfk= input variable indicating the storage option of each k-point
!! natom = number of atoms in the unit cell
!! nband= size of e_kpt
!! nkpt= number of k-points
!! npwarr= array holding npw for each k point
!! nr1,nr2,nr3 = grid size (nr1 x nr2 x nr3 = filrho dimension)
!! nspinor= number of spinorial components of the wavefunctions
!! nsppol= number of spin polarization
!! ntypat = number of atom type
!! rprim = orientation of the unit cell axes
!! xcart = cartesian coordinates
!! typat= input variable typat(natom)
!! znucl= znucltypat(ntypat) from alchemy
!!
!! OUTPUT
!! Depends on the option chosen.
!! It is the wave function for the k point, band and spin polarisation
!! chosen.  It can be written in different ways. The option are describe
!! with the option list.  It is possible to output a Data Explorer file.
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!
!! SOURCE

subroutine cut3d_wffile(wfk_fname,ecut,exchn2n3d,istwfk,kpt,natom,nband,nkpt,npwarr,&
&  nr1,nr2,nr3,nspinor,nsppol,ntypat,rprimd,xcart,typat,znucl)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,natom,nkpt,nr1,nr2,nr3,nspinor,nsppol
 integer,intent(in) :: ntypat
 real(dp),intent(in) :: ecut
 character(len=*),intent(in) :: wfk_fname
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt),typat(natom)
 real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3),znucl(ntypat)
 real(dp),intent(in) :: xcart(3,natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf0=0,tim_rwwf0=0,ndat1=1,formeig0=0
 integer :: cband,cgshift,ckpt,cplex,cspinor,csppol,gridshift1
 integer :: gridshift2,gridshift3,ia,iatom,iband,ichoice,ifile,iomode
 integer :: ii1,ii2,ii3,ikpt,ilang,ioffkg,iout,iprompt,ipw
 integer :: ir1,ir2,ir3,ivect,ixint,mband,mbess,mcg,mgfft
 integer :: mkmem,mlang,mpw,n4,n5,n6,nfit,npw_k
 integer :: nradintmax,oldcband,oldckpt,oldcspinor,oldcsppol
 integer :: prtsphere,select_exit,unout,iunt,rc_ylm
 integer :: ikpt_qps,nkpt_qps,nband_qps,iscf_qps
 real(dp) :: arg,bessargmax,bessint_delta,kpgmax,ratsph,tmpi,tmpr,ucvol,weight,eig_k_qps
 character(len=*), parameter :: INPUTfile='cut.in'
 character(len=1) :: outputchar
 character(len=10) :: string
 character(len=4) :: mode_paral
 character(len=500) :: msg
 character(len=fnlen) :: output,output1
 type(MPI_type) :: mpi_enreg
 type(wfk_t) :: Wfk
 type(jlspline_t) :: jlspl
!arrays
 integer :: atindx(natom),iatsph(natom),ngfft(18),nradint(natom),mlang_type(ntypat)
 integer,allocatable :: gbound(:,:),iindex(:),kg(:,:),kg_dum(:,:),kg_k(:,:)
 integer,allocatable :: npwarr1(:),npwarrk1(:),npwtot1(:)
 real(dp) :: cmax(natom),gmet(3,3),gprimd(3,3)
 real(dp) :: phkxred(2,natom),ratsph_arr(natom),rmet(3,3),shift_tau(3)
 real(dp) :: tau2(3,natom),xred(3,natom),kpt_qps(3)
 real(dp) :: znucl_atom(natom)
 integer  :: znucl_atom_int(natom)
 real(dp),allocatable :: bess_fit(:,:,:)
 real(dp),allocatable :: cg_k(:,:),cgcband(:,:),denpot(:,:,:),eig_k(:)
 real(dp),allocatable :: fofgout(:,:),fofr(:,:,:,:),k1(:,:)
 real(dp),allocatable :: kpgnorm(:),occ_k(:),ph1d(:,:),ph3d(:,:,:),rint(:)
 real(dp),allocatable :: sum_1ll_1atom(:,:,:),sum_1lm_1atom(:,:,:)
 real(dp),allocatable :: cplx_1lm_1atom(:,:,:,:)
 real(dp),allocatable :: xfit(:),yfit(:),ylm_k(:,:)
 real(dp),allocatable :: ylmgr_dum(:,:,:)
 character(len=fnlen) :: fileqps
 character(len=fnlen),allocatable :: filename(:)
 complex(dp),allocatable :: ccoeff(:,:),wfg(:,:),wfg_qps(:)
 real(dp) :: spinvec(3)

! ***********************************************************************

 call initmpi_seq(mpi_enreg)
 mband=maxval(nband)
 ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt,mband,nsppol))
 mpi_enreg%proc_distrb=0
 mpi_enreg%me_g0 = 1
 oldckpt=0
 oldcband=0
 oldcsppol=0
 oldcspinor=0

 iout=-1
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!get xred
 call xcart2xred(natom,rprimd,xcart,xred)

!znucl indexed by atoms
 znucl_atom     =     znucl(typat(1:natom))
 znucl_atom_int = INT(znucl(typat(1:natom)))

 do iatom=1,natom
   iatsph(iatom) = iatom
   atindx(iatom) = iatom
 end do

!max ang mom + 1
 mlang = 5

 ABI_ALLOCATE(kg_dum,(3,0))

 ABI_ALLOCATE(ph1d,(2,(2*nr1+1+2*nr2+1+2*nr3+1)*natom))
 call getph(atindx,natom,nr1,nr2,nr3,ph1d,xred)

 do
!  Get k-point, band and spin polarisation for the output
   if(nkpt/=1)then
     write(std_out,*)
     write(std_out,'(a,i4,a)') ' For which k-points? (1 to ',nkpt,')'
     read(std_in,*)ckpt
!    Check if kpt exist
     if(ckpt<1 .or. ckpt>nkpt) then
       write(msg,'(a,i0)') 'Invalid k-point ',ckpt
       MSG_ERROR(msg)
     end if
   else
     ckpt=nkpt
   end if
   write(std_out,*) ' => Your k-point is : ',ckpt
   write(std_out,*)

   if(nband(ckpt)/=1)then
     write(std_out,*)
     write(std_out,'(a,i5,a)') ' For which band ? (1 to ',nband(ckpt),')'
     read(std_in,*)cband
!    Check if band number exist

     if(cband<1 .or. cband>nband(ckpt)) then
       write(msg,'(a,i0)')'Invalid band number',cband
       MSG_ERROR(msg)
     end if
   else
     cband=nband(ckpt)
   end if
   write(std_out,*) ' => Your band number is : ',(cband)
   write(std_out,*)

   if(nsppol/=1)then
     write(std_out,*)
     write(std_out,*) ' For which spin polarisation ?'
     read(std_in,*)csppol
!    Check if spin polarisation exist
     if(csppol<1 .or. csppol>nsppol) then
       write(msg,'(a,i0)')'Invalid spin polarisation ',csppol
       MSG_ERROR(msg)
     end if
   else
     csppol=1
   end if

   write(std_out,*) ' => Your spin polarisation number is : ',(csppol)
   write(std_out,*)

   if(nspinor/=1) then
     write(std_out,*) ' nspinor = ', nspinor
     write(std_out,*)
     write(std_out,*) ' For which spinor component ?'
     read(std_in,*) cspinor
!    Check if spin polarisation exist
     if(cspinor<1 .or. cspinor>nspinor) then
       write(msg,'(a,i0)')'Invalid spinor index ',cspinor
       MSG_ERROR(msg)
     end if
     write(std_out,*) ' => Your spinor component is : ',(cspinor)
     write(std_out,*)
   else
     cspinor=1
   end if

!  Reading of the data if the value of ckpt and csppol are different from oldckpt and oldcsppol
!  formeig=0 gstate calculation
!  formeig=1 for response function calculation
   if(csppol/=oldcsppol .or. ckpt/=oldckpt)then
     mband=maxval(nband)
     mpw=maxval(npwarr)
     mcg=mpw*nspinor*mband
     if (allocated (cg_k))   then
       ABI_DEALLOCATE(cg_k)
       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(occ_k)
     end if
     ABI_ALLOCATE(cg_k,(2,mcg))
     ABI_ALLOCATE(eig_k,((2*mband)**formeig0*mband))
     ABI_ALLOCATE(occ_k,(mband))

!    FIXME
!    nband depends on (kpt,spin)
     iomode = iomode_from_fname(wfk_fname)
     call wfk_open_read(Wfk,wfk_fname,formeig0,iomode,get_unit(),xmpi_comm_self)
     call wfk%read_band_block([1,nband(ckpt)],ckpt,csppol,xmpio_single,cg_k=cg_k,eig_k=eig_k,occ_k=occ_k)
     call wfk%close()
   end if

   if (csppol/=oldcsppol .or. ckpt/=oldckpt .or. cband/=oldcband .or. cspinor/=oldcspinor ) then
!    The data of ckpt,cnsspol are in cg_k. Now we have to do the Fourier Transform of the datas

     ngfft(1)=nr1
     ngfft(2)=nr2
     ngfft(3)=nr3
!    ngfft(4) and ngfft(5) can not be even (see getng.f)
     if (mod(nr1,2)==0)then
       ngfft(4)=nr1+1
     else
       ngfft(4)=nr1
     end if
     if (mod(nr2,2)==0)then
       ngfft(5)=nr2+1
     else
       ngfft(5)=nr2
     end if
     ngfft(6)=nr3
!    XG 020829 : 112 does not work yet for all istwfk values
     ngfft(7)=111
     ngfft(8)=256
     ngfft(9)=0
     ngfft(10)=1
     ngfft(11)=0
     ngfft(12)=ngfft(2)
     ngfft(13)=ngfft(3)
     ngfft(14)=0

!    if iout<0, the output of metric will not be print
     mode_paral='PERS'
     mkmem=nkpt
     mgfft=maxval(ngfft(1:3))
     ABI_ALLOCATE(npwarr1,(nkpt))
     ABI_ALLOCATE(kg,(3,mpw*mkmem))
     ABI_ALLOCATE(npwtot1,(nkpt))
     call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfft(2),ngfft(3),'all')

!    Create positions index for pw
     call kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kpt,mkmem,nband,nkpt,&
&     mode_paral,mpi_enreg,mpw,npwarr1,npwtot1,nsppol)

     ioffkg=0
     do ikpt=1,ckpt-1
       ioffkg=ioffkg+npwarr1(ikpt)
     end do
     npw_k=npwarr(ckpt)
     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     ABI_ALLOCATE(kg_k,(3,npw_k))
     kg_k(:,1:npw_k)=kg(:,1+ioffkg:npw_k+ioffkg)

     ABI_ALLOCATE(ylm_k,(npw_k,mlang*mlang))
     ABI_ALLOCATE(ylmgr_dum,(npw_k,3,mlang*mlang))

!    call for only the kpoint we are interested in !
     ABI_ALLOCATE(k1,(3,1))
     k1(:,1)=kpt(:,ckpt)
     ABI_ALLOCATE(npwarrk1,(1))
     npwarrk1 = (/npw_k/)
     call initylmg(gprimd,kg_k,k1,1,mpi_enreg,mlang,npw_k,nband,1,&
&     npwarrk1,nsppol,0,rprimd,ylm_k,ylmgr_dum)
     ABI_DEALLOCATE(ylmgr_dum)
     ABI_DEALLOCATE(k1)
     ABI_DEALLOCATE(npwarrk1)

!    Compute the norms of the k+G vectors
     ABI_ALLOCATE(kpgnorm,(npw_k))
     call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)

     call sphereboundary(gbound,istwfk(ckpt),kg_k,mgfft,npw_k)
!    Do the Fourier Transform
     n4=ngfft(4)
     n5=ngfft(5)
     n6=ngfft(6)
!    cplex=0
     cplex=1
!    Complex can be set to 0 with this option(0) of fourwf

!    Read the QPS file if GW wavefunctions are to be analysed
     write(std_out,*) 'Do you want to analyze a GW wavefunction? (1=yes,0=no)'
     read(std_in,*) ii1
     write(std_out,*) '=> Your choice is :',ii1
     write(std_out,*)

     if(ii1==1) then
       write(std_out,*) 'What is the name of the QPS file?'
       if (read_string(fileqps, unit=std_in) /= 0) then
         MSG_ERROR("Fatal error!")
       end if
!      Checking the existence of data file
       if (.not. file_exists(fileqps)) then
         MSG_ERROR(sjoin('Missing data file:', fileqps))
       end if

       if (open_file(fileqps, msg, newunit=iunt, status='old',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if

       read(iunt,*) iscf_qps
       read(iunt,*) nkpt_qps
       read(iunt,*) nband_qps
       read(iunt,*) ikpt_qps

       ABI_ALLOCATE(ccoeff,(nband_qps,nband_qps))
       do ikpt=1,ckpt ! nkpt_qps
         read(iunt,*) kpt_qps(:)
         do iband=1,nband_qps
           read(iunt,*) eig_k_qps
           read(iunt,*) ccoeff(:,iband)
         end do
       end do
       close(iunt)

       ABI_ALLOCATE(wfg,(npw_k,nband_qps))
       ABI_ALLOCATE(wfg_qps,(npw_k))
       do iband=1,nband_qps
         cgshift=(iband-1)*npw_k*nspinor + (cspinor-1)*npw_k
         wfg(:,iband) = dcmplx( cg_k(1,cgshift+1:cgshift+npw_k),cg_k(2,cgshift+1:cgshift+npw_k) )
       end do

       wfg_qps = matmul( wfg(:,:) , ccoeff(:,cband) )

!      write(std_out,*) 'norm',SUM( abs(wfg(:,cband))**2 )
!      write(std_out,*) 'norm',SUM( abs(wfg_qps(:))**2 )
       ABI_DEALLOCATE(ccoeff)
       ABI_DEALLOCATE(wfg)
       ABI_ALLOCATE(cgcband,(2,npw_k*nspinor))
       cgcband = zero
       cgcband(1,:)= real(wfg_qps(:))
       cgcband(2,:)= aimag(wfg_qps(:))
       ABI_DEALLOCATE(wfg_qps)

     else ! not a GW wavefunction

! get spin vector for present state
       cgshift=(cband-1)*npw_k*nspinor
       ABI_ALLOCATE(cgcband,(2,npw_k*nspinor))
       cgcband(:,1:npw_k*nspinor)=cg_k(:,cgshift+1:cgshift+nspinor*npw_k)
     end if ! test QPS wavefunction from GW

     if (nspinor == 2) then
       call cg_getspin(cgcband, npw_k, spinvec)
       write(std_out,'(a,6E20.10)' ) ' spin vector for this state = ', (spinvec)
     end if

!    Fix the phase of cgcband, for portability reasons
!    call fxphas(cgcband,cgcband,0,npw_k,1,npw_k,0)

     ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
     ABI_ALLOCATE(fofgout,(2,npw_k))
     ABI_ALLOCATE(fofr,(2,n4,n5,n6))

     call fourwf(cplex,denpot,cgcband(:,(cspinor-1)*npw_k+1:cspinor*npw_k),fofgout,fofr,gbound,gbound,&
&     istwfk(ckpt),kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,&
&     npw_k,n4,n5,n6,0,tim_fourwf0,weight,weight)

!    TODO
!    call fft_ug_dp(npw_k,nfft,nspinor,ndat1,mgfft,ngfft,istwf_k(ckpt),kg_k,gbound,cgcband,fofr)

!    Analyse wavefunction inside atomic sphere

     write(std_out,'(a)' ) ' Do you want the atomic analysis for this state : '
     write(std_out,'(a,2i5,a)' ) ' (kpt,band)= (',ckpt,cband,')? '
     write(std_out,'(a)' ) ' If yes, enter the radius of the atomic spheres, in bohr '
     write(std_out,'(a)' ) ' If no, enter 0 '
     read (std_in,*) ratsph
     write(std_out,'(a,f16.8,a)' ) ' You entered ratsph=',ratsph,' Bohr '

     if (ratsph >= tol10) then

       write(std_out,'(3a)' ) ch10,' Atomic sphere analysis ',ch10

!      Init bessel function integral for recip_ylm: max ang mom + 1
       mlang = 5
       bessint_delta = 0.1_dp
       kpgmax = sqrt(ecut)
       bessargmax = ratsph*two_pi*kpgmax
       mbess = int (bessargmax / bessint_delta) + 1
       bessargmax = bessint_delta*mbess

!      Intervals in radial integration
       nradintmax = mbess
       nradint(1:natom)=nradintmax

       write(std_out,'(a,2es16.6,i6)')' wffile : kpgmax, bessargmax, nradint = ', kpgmax, bessargmax,nradintmax

!      Initialize general Bessel function array on uniform grid xx, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
       ABI_ALLOCATE(rint,(nradintmax))

       jlspl = jlspline_new(mbess, bessint_delta, mlang)

       ABI_ALLOCATE(bess_fit,(mpw,nradintmax,mlang))
       ABI_ALLOCATE(xfit,(npw_k))
       ABI_ALLOCATE(yfit,(npw_k))
       ABI_ALLOCATE(iindex,(npw_k))
       nfit = npw_k

       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*ratsph / (nradintmax-1)

         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw) * rint(ixint)
           iindex(ipw) = ipw
         end do
         call sort_dp (npw_k,xfit,iindex,tol14)
         do ilang=1,mlang
           call splint(mbess,jlspl%xx,jlspl%bess_spl(:,ilang),jlspl%bess_spl_der(:,ilang),nfit,xfit,yfit)
!          Re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do ! ipw
       end do ! ixint

!      Construct phases ph3d for all G vectors in present sphere make phkred for all atoms
       do ia=1,natom
         iatom=atindx(ia)
         arg=two_pi*( kpt(1,ckpt)*xred(1,ia) + kpt(2,ckpt)*xred(2,ia) + kpt(3,ckpt)*xred(3,ia))
         phkxred(1,iatom)=cos(arg)
         phkxred(2,iatom)=sin(arg)
       end do

       ABI_ALLOCATE(ph3d,(2,npw_k,natom))
!      Get full phases exp (2 pi i (k+G).x_tau) in ph3d
       call ph1d3d(1,natom,kg_k,natom,natom,npw_k,nr1,nr2,nr3,phkxred,ph1d,ph3d)

       ABI_ALLOCATE(sum_1ll_1atom,(nspinor**2,mlang,natom))
       ABI_ALLOCATE(sum_1lm_1atom,(nspinor**2,mlang**2,natom))
       ABI_ALLOCATE(cplx_1lm_1atom,(2,nspinor,mlang**2,natom))
       prtsphere=1
       ratsph_arr(:)=ratsph

       rc_ylm = 1 ! Real or Complex spherical harmonics.
       mlang_type = 5

       call recip_ylm (bess_fit,cgcband,istwfk(ckpt),mpi_enreg,&
&       nradint,nradintmax,mlang,mpw,natom,typat,mlang_type,npw_k,nspinor,ph3d,prtsphere,rint,&
&       ratsph_arr,rc_ylm,sum_1ll_1atom,sum_1lm_1atom,cplx_1lm_1atom,ucvol,ylm_k,znucl_atom)

       call dens_in_sph(cmax,cgcband(:,(cspinor-1)*npw_k+1:cspinor*npw_k),gmet,istwfk(ckpt),&
&       kg_k,natom,ngfft,mpi_enreg,npw_k,ph1d,ratsph_arr,ucvol)

       write(std_out,'(a)' )' Charge in the sphere around each atom '
       do iatom=1,natom
         write(std_out,'(a,i4,a,f14.8)' ) ' Atom number ',iatom,' :  charge =',cmax(iatom)
       end do

       ABI_DEALLOCATE(sum_1ll_1atom)
       ABI_DEALLOCATE(sum_1lm_1atom)
       ABI_DEALLOCATE(cplx_1lm_1atom)
       ABI_DEALLOCATE(ph3d)
       ABI_DEALLOCATE(iindex)
       ABI_DEALLOCATE(yfit)
       ABI_DEALLOCATE(xfit)
       ABI_DEALLOCATE(bess_fit)
       call jlspline_free(jlspl)
       ABI_DEALLOCATE(rint)
     end if ! ratsph < 0     = end if for atomic sphere analysis

     ABI_DEALLOCATE(cgcband)
     ABI_DEALLOCATE(fofgout)
     ABI_DEALLOCATE(denpot)
     ABI_DEALLOCATE(gbound)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(npwarr1)
     ABI_DEALLOCATE(kg)
     ABI_DEALLOCATE(npwtot1)
     ABI_DEALLOCATE(kpgnorm)
     ABI_DEALLOCATE(ylm_k)
     call destroy_distribfft(mpi_enreg%distribfft)
   end if

   write(std_out,*)
   write(std_out,*) ' 3D wave function was read. ','Ready for further treatment.'
   write(std_out,*)
   write(std_out,*) '============================','==============================='
   write(std_out,*)

!  ------------------------------------------------------------------------

!  At this moment all the input is done
!  The code knows the geometry of the system, and the data file.


   select_exit = 0
   do while (select_exit == 0)
     write(std_out,*) ' What is your choice ? Type:'
     write(std_out,*) '  0 => exit to k-point / band / spin-pol loop'
     write(std_out,*) '  1 => 3D formatted real and imaginary data'
     write(std_out,*) '       (output the bare 3D data - two column,R,I)'
     write(std_out,*) '  2 => 3D formatted real data'
     write(std_out,*) '       (output the bare 3D data - one column)'
     write(std_out,*) '  3 => 3D formatted imaginary data'
     write(std_out,*) '       (output the bare 3D data - one column)'
     write(std_out,*) '  4 => 3D indexed real and imaginary data'
     write(std_out,*) '       (3D data, preceeded by 3D index)'
     write(std_out,*) '  5 => 3D indexed real data'
     write(std_out,*) '       (bare 3D data, preceeded by 3D index)'
     write(std_out,*) '  6 => 3D indexed imaginary data'
     write(std_out,*) '       (bare 3D data, preceeded by 3D index)'
     write(std_out,*) '  7 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Real file and Imaginary file)'
     write(std_out,*) '  8 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Only the Real file)'
     write(std_out,*) '  9 => 3D Data Explorer formatted data '
     write(std_out,*) '       (Only the Imaginary file)'
     write(std_out,*) ' 10 => 3D Data Explorer formatted data and position files'
     write(std_out,*) ' 11 => XCrysden formatted data (norm of wf) and position files'
     write(std_out,*) ' 12 => NetCDF data and position file'
     write(std_out,*) ' 13 => XCrysden/VENUS wavefunction (real part of data)'
     write(std_out,*) ' 14 => Gaussian/cube wavefunction module'
     read(std_in,*) ichoice
     write(std_out,'(a,a,i2,a)' ) ch10,' Your choice is ',ichoice,char(10)

     if (ichoice>0 .and. ichoice<15)then
       write(std_out,*) ch10,'  Enter the root of an output file:'
       if (read_string(output1, unit=std_in) /= 0) then
         MSG_ERROR("Fatal error!")
       end if
       write(std_out,*) '  The root of your file is : ',trim(output1)
       output=trim(output1)
       call int2char10(ckpt,string)
       output=trim(output)//'_k'//trim(string)
       call int2char10(cband,string)
       output=trim(output)//'_b'//trim(string)
       if (nsppol > 1) then
         call int2char10(csppol,string)
         output=trim(output)//'_sppol'//trim(string)
       end if
       if (nspinor > 1) then
         call int2char10(cspinor,string)
         output=trim(output)//'_spinor'//trim(string)
       end if

       write(std_out,*) '  The corresponding filename is : ',trim(output)
     end if

     select case (ichoice)

     case (1) ! data R,I
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted real and imaginary data'
       write(std_out,*) 'The first column is the real data'
       write(std_out,*) 'The second column is the imaginary data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_ri("RI",nr1,nr2,nr3,n4,n5,n6,fofr,unit=unout)
       close(unout)
       exit

     case (2) ! data R
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted real data'
       write(std_out,*) 'The only column is the real data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_ri("R",nr1,nr2,nr3,n4,n5,n6,fofr,unit=unout)
       close(unout)
       exit

     case (3) ! data I
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted real data'
       write(std_out,*) 'The only column is the imaginary data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_ri("I",nr1,nr2,nr3,n4,n5,n6,fofr,unit=unout)
       close(unout)
       exit

     case (4) ! coord(x,y,z) data R,I
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted data'
       write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
       write(std_out,*) 'The fourth column is the real data'
       write(std_out,*) 'The fifth column is the imaginary data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_xyzri("RI",nr1,nr2,nr3,n4,n5,n6,fofr,rprimd,conv_fact=Bohr_Ang,unit=unout)
       close(unout)
       exit

     case (5) ! coord(x,y,z) data R
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted data'
       write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
       write(std_out,*) 'The fourth column is the real data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_xyzri("R",nr1,nr2,nr3,n4,n5,n6,fofr,rprimd,conv_fact=Bohr_Ang,unit=unout)
       close(unout)
       exit

     case(6) ! coord(x,y,z) data I
       write(std_out,*)
       write(std_out,*) 'Give 1 file of 3D formatted data'
       write(std_out,*) 'The first three columns are the x,y,z positions(Angstrom)'
       write(std_out,*) 'The fourth column is the imaginary data'
       write(std_out,*)
       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_xyzri("I",nr1,nr2,nr3,n4,n5,n6,fofr,rprimd,conv_fact=Bohr_Ang,unit=unout)
       close(unout)
       exit

     case(7) !OpenDX format, data R and data I
       write(std_out,*)
       write(std_out,*) 'Give 2 files of 3D formatted data'
       write(std_out,*) 'The file is ready to be use with OpenDX'
       write(std_out,*) 'The eig_kvalues and occupations numbers are in comments'
       write(std_out,*)
       ABI_ALLOCATE(filename,(2))
       filename(1)=trim(output)//'Real.dx'
       filename(2)=trim(output)//'Imag.dx'
       write(std_out,*) '  The name of your files is : '
       write(std_out,*) trim(filename(1)),'  for the real part,'
       write(std_out,*) trim(filename(2)),'  for the imaginary part.'
       write(std_out,*)

       do ifile=1,2
         if (open_file(filename(ifile), msg, newunit=unout, status='replace',form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if
         rewind(unout)
         write(unout,*)'# band,  eig_kvalues   and   occupations'
         do iband=1,nband(ckpt)
           write(unout,'(a,i4,2f20.16)')'#',iband,eig_k(iband),occ_k(iband)
         end do
         write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',nr1*nr2*nr3,' data follows'
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)')fofr(ifile,ir1,ir2,ir3)
             end do
           end do
         end do

         write(unout,'(a)')'# this is the object defining the grid connections'
         write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',nr3,nr2,nr1
         write(unout,*)
         write(unout,*)
         write(unout,'(a)')'# this is the object defining the grid'
         write(unout,'(a,3i5)')'object "positions" class gridpositions counts',nr3,nr2,nr1

         write(unout,'(a)') 'origin 0 0 0'
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

         write(unout,'(a)')'# this is the collective object, one for each grid '
         write(unout,'(a)')'object "densite" class field '
         write(unout,'(a)')'component "positions"   value "positions"'
         write(unout,'(a)')'component "connections" value "gridconnections" '
         write(unout,'(a)')'component "data"        value "donnees"'

         close(unit=unout)
       end do
       ABI_DEALLOCATE(filename)
       exit

     case(8) ! OpenDX format, data R and data I
       write(std_out,*)
       write(std_out,*) 'Give 2 files of 3D formatted data'
       write(std_out,*) 'The file is ready to be use with OpenDX'
       write(std_out,*) 'The eig_kvalues and occupations numbers are in comments'
       write(std_out,*)
       ABI_ALLOCATE(filename,(1))
       filename(1)=trim(output)//'Real.dx'
       write(std_out,*) '  The name of your file is : '
       write(std_out,*) trim(filename(1)),'  for the real part,'
       write(std_out,*)

       if (open_file(filename(1), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       rewind(unout)
       write(unout,*)'# band,  eig_kvalues   and   occupations'
       do iband=1,nband(ckpt)
         write(unout,'(a,i4,2f20.16)')'#',iband,eig_k(iband),occ_k(iband)
       end do
       write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',nr1*nr2*nr3,' data follows'
       do ir3=1,nr3
         do ir2=1,nr2
           do ir1=1,nr1
             write(unout,'(f20.16)')fofr(1,ir1,ir2,ir3)
           end do
         end do
       end do

       write(unout,'(a)')'# this is the object defining the grid connections'
       write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',nr3,nr2,nr1
       write(unout,*)
       write(unout,*)
       write(unout,'(a)')'# this is the object defining the grid'
       write(unout,'(a,3i5)')'object "positions" class gridpositions counts',nr3,nr2,nr1

       write(unout,'(a)') 'origin 0 0 0'
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

       write(unout,'(a)')'# this is the collective object, one for each grid '
       write(unout,'(a)')'object "densite" class field '
       write(unout,'(a)')'component "positions"   value "positions"'
       write(unout,'(a)')'component "connections" value "gridconnections" '
       write(unout,'(a)')'component "data"        value "donnees"'

       close(unit=unout)
       ABI_DEALLOCATE(filename)
       exit

     case(9) !OpenDX format, data R and data I
       write(std_out,*)
       write(std_out,*) 'Give 2 files of 3D formatted data'
       write(std_out,*) 'The file is ready to be use with OpenDX'
       write(std_out,*) 'The eig_kvalues and occupations numbers are in comments'
       write(std_out,*)
       ABI_ALLOCATE(filename,(1))
       filename(1)=trim(output)//'Imag.dx'
       write(std_out,*) '  The name of your file is : '
       write(std_out,*) trim(filename(1)),'  for the imaginary part.'
       write(std_out,*)

       if (open_file(filename(1), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       rewind(unout)
       write(unout,*)'# band,  eig_kvalues   and   occupations'
       do iband=1,nband(ckpt)
         write(unout,'(a,i4,2f20.16)')'#',iband,eig_k(iband),occ_k(iband)
       end do
       write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',nr1*nr2*nr3,' data follows'
       do ir3=1,nr3
         do ir2=1,nr2
           do ir1=1,nr1
             write(unout,'(f20.16)')fofr(2,ir1,ir2,ir3)
           end do
         end do
       end do

       write(unout,'(a)')'# this is the object defining the grid connections'
       write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',nr3,nr2,nr1
       write(unout,*)
       write(unout,*)
       write(unout,'(a)')'# this is the object defining the grid'
       write(unout,'(a,3i5)')'object "positions" class gridpositions counts',nr3,nr2,nr1

       write(unout,'(a)') 'origin 0 0 0'
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
       write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

       write(unout,'(a)')'# this is the collective object, one for each grid '
       write(unout,'(a)')'object "densite" class field '
       write(unout,'(a)')'component "positions"   value "positions"'
       write(unout,'(a)')'component "connections" value "gridconnections" '
       write(unout,'(a)')'component "data"        value "donnees"'

       close(unit=unout)
       ABI_DEALLOCATE(filename)
       exit

     case(10)           !OpenDX format, data R and data I, atoms positions, lattice and cell
       write(std_out,*)
       write(std_out,*) 'Give 5 files of formatted data'
       write(std_out,*) 'The files are ready to be use with Data Explorer'
       write(std_out,*) 'The eig_kvalues and occupations numbers are in comments'
       write(std_out,*) 'of the two data files'
       write(std_out,*)
       ABI_ALLOCATE(filename,(2))
       filename(1)=trim(output)//'Real.dx'
       filename(2)=trim(output)//'Imag.dx'
       write(std_out,*) '  The name of your data files is : '
       write(std_out,*) trim(filename(1)),'  for the real part,'
       write(std_out,*) trim(filename(2)),'  for the imaginary part.'
       write(std_out,*)

       do ifile=1,2
         if (open_file(filename(ifile), msg, newunit=unout, status='replace',form='formatted') /= 0) then
           MSG_ERROR(msg)
         end if
         rewind(unout)
         do iband=1,nband(ckpt)
           write(unout,'(a,2f20.16)')'#', eig_k(iband),occ_k(iband)
         end do
         write(unout,'(a,i10,a)')'object "donnees" class array type float rank 0 items',nr1*nr2*nr3,' data follows'
         do ir3=1,nr3
           do ir2=1,nr2
             do ir1=1,nr1
               write(unout,'(f20.16)')fofr(ifile,ir1,ir2,ir3)
             end do
           end do
         end do

         write(unout,'(a)')'# this is the object defining the grid connections'
         write(unout,'(a,3i5)')'object "gridconnections" class gridconnections counts',nr3,nr2,nr1
         write(unout,*)
         write(unout,*)
         write(unout,'(a)')'# this is the object defining the grid'
         write(unout,'(a,3i5)')'object "positions" class gridpositions counts',nr3,nr2,nr1

         write(unout,'(a)') 'origin 0 0 0'
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,3)/nr3, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,2)/nr2, ii1=1,3)
         write(unout,'(a,3f16.10)')'delta ',(Bohr_Ang*rprimd(ii1,1)/nr1, ii1=1,3)

         write(unout,'(a)')'# this is the collective object, one for each grid '
         write(unout,'(a)')'object "densite" class field '
         write(unout,'(a)')'component "positions"   value "positions"'
         write(unout,'(a)')'component "connections" value "gridconnections" '
         write(unout,'(a)')'component "data"        value "donnees"'

         close(unit=unout)
       end do
       ABI_DEALLOCATE(filename)
!
!        write LATTICE_VEC.dx file
!
       ABI_ALLOCATE(filename,(3))
       filename(1)=trim(output1)//'_LATTICE_VEC.dx'
       filename(2)=trim(output1)//'_ATOM_POS.dx'
       filename(3)=trim(output1)//'_UCELL_FRAME.dx'
       write(std_out,*)
       write(std_out,*)'Give the lattice file, ', trim(filename(1))
       if (open_file(filename(1), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if

       write(unout,'("#",/,"#",/,"#    LATTICE VECTOR INFO:",/,"#",/,"#")')
       write(unout,'(a)') 'object "lattices" class array type float rank 1 shape 3 items 3 data follows'
       do ivect=1,3
         write(unout,'(3f16.10)')  Bohr_Ang*rprimd(1,ivect),Bohr_Ang*rprimd(2,ivect),Bohr_Ang*rprimd(3,ivect)
       end do
       write(unout,'(a,a)') 'object "lattices_location" class array type float ','rank 1 shape 3 items 3 data follows'
       do ivect=1,3
         write(unout,'(3f16.10)')  0_dp,0_dp,0_dp
       end do
       write(unout,'("object   3 class field")')
       write(unout,'(a)') 'component "data" value "lattices"'
       write(unout,'(a)') 'component "positions" value "lattices_location"'
       close(unout)
!
!        write ATOM_POS.dx file
!
       write(std_out,*)'Give the atoms positions file, ', trim(filename(2))

       if (open_file(filename(2), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if

       write(unout,'("#",/,"#",/,"#    BALL AND STICK INFO:",/,"#",/,"#")')
       write(unout,'(a,i5,a)') 'object "atomcoord" array type float rank 1 shape 3 items ',natom,' data follows'
       do iatom=1,natom
         write(unout,'(3f16.10)')  Bohr_Ang*xcart(1:3,iatom)
       end do
!        write(unout,'(a,i5,a)') 'object "data" array type string rank 0 shape 2 items ',natom,' data follows'
       write(unout,'(a,i5,a)') 'object "colorcode" array type float rank 0 items ',natom,' data follows'
       do iatom=1,natom
         write(unout,'(f10.4)') znucl(typat(iatom))
       end do
       write(unout,'(a)') 'object "molecule" field'
       write(unout,'(a)') 'component "positions" value "atomcoord"'
       write(unout,'(a)') 'component "data" value "colorcode"'
       close(unout)

!
!        write UCELL_FRAME.dx file
!
       write(std_out,*)'Give the enveloppe of the cell file, ',trim(filename(3))
       if (open_file(filename(3), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if

       write(unout,'("#",/,"#",/,"#    UNIT CELL FRAME INFO:",/,"#",/,"#")')
       write(unout,'(a)')'object 3 class array type int rank 1 shape 2 items 12 data follows'
       write(unout,'(" 0  1",/," 0  2",/," 0  3",/," 1  4",/," 1  5",/," 3  5")')
       write(unout,'(" 3  6",/," 2  6",/," 2  4",/," 7  5",/," 7  6",/," 7  4")')
       write(unout,'(a)') 'attribute "element type" string "lines"'
       write(unout,'("object  4 class array type float rank 1 shape 3 items    8 data follows")')
       write(unout,'("      .00000000      .00000000      .00000000")')
       write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,1)
       write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,2)
       write(unout,'(3f20.10)') Bohr_Ang*rprimd(:,3)
       write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,2))
       write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,3))
       write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,2)+rprimd(:,3))
       write(unout,'(3f20.10)') Bohr_Ang*(rprimd(:,1)+rprimd(:,2)+rprimd(:,3))
       write(unout,'("object 5 array type float rank 0 items 12 data follows")')
       do ivect=1,12
         write(unout,'("1.0")')
       end do
       write(unout,'(a)') 'attribute "dep" string "connections"'
       write(unout,'("object 6 class field")')
       write(unout,'(a)') 'component "data" value 5'
       write(unout,'(a)') 'component "positions" value 4'
       write(unout,'(a)') 'component "connections" value 3'
       close(unout)
       ABI_DEALLOCATE(filename)

       write(std_out,*)
       exit

     case(11)
       write(std_out,*)
       write(std_out,*) 'Give 1 files of formatted data'
       write(std_out,*) 'The files are ready to be used with XCrysDen'
       write(std_out,*)
       gridshift1 = 0
       gridshift2 = 0
       gridshift3 = 0
       write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
       write(std_out,*)
       shift_tau(:) = 0.0
       if (read_string(outputchar, unit=std_in) /= 0) then
         MSG_ERROR("Fatal error!")
       end if
       if (outputchar == 'y' .or. outputchar == 'Y') then
         MSG_ERROR("Shift is buggy, don't use it")
         write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
         write(std_out,*)
         read (std_in,*) gridshift1, gridshift2, gridshift3
         shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
       end if

       ABI_ALLOCATE(filename,(1))
       filename(1)=trim(output)
       write(std_out,*) '  The name of your data files is : '
       write(std_out,*) trim(filename(1)),'  for the density (norm of the wfk),'
       write(std_out,*)

       if (open_file(filename(1), msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       rewind(unout)
       do iband=1,nband(ckpt)
         write(unout,'(a,2f20.16)')'#', eig_k(iband),occ_k(iband)
       end do

       write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
       write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
       write(std_out,*) ' znucl = ', znucl, ' typat = ', typat, ' ntypat = ', ntypat

       write(unout,'(1X,A)')  'DIM-GROUP'
       write(unout,*) '3  1'
       write(unout,'(1X,A)') 'PRIMVEC'
       do ir1 = 1,3
         write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
       end do
       write(unout,'(1X,A)') 'PRIMCOORD'
       write(unout,*) natom, ' 1'
!
!        generate translated coordinates to match density shift
!
       do iatom = 1,natom
         tau2 (:,iatom) = xcart(:,iatom) - shift_tau(:)
       end do

       do iatom = 1,natom
         write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&         Bohr_Ang*tau2(1,iatom), &
&         Bohr_Ang*tau2(2,iatom), &
&         Bohr_Ang*tau2(3,iatom)
       end do
       write(unout,'(1X,A)') 'ATOMS'
       do iatom = 1,natom
         write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&         Bohr_Ang*tau2(1,iatom), &
&         Bohr_Ang*tau2(2,iatom), &
&         Bohr_Ang*tau2(3,iatom)
       end do

!        write(unout,'(1X,A)') 'FRAMES'
       write(unout,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
       write(unout,*) 'datagrids'
       write(unout,'(1X,A)') 'DATAGRID_3D_DENSITY'
       write(unout,*) nr1+1,nr2+1,nr3+1
       write(unout,*) '0.0 0.0 0.0 '
       do ir1 = 1,3
         write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
       end do

       do ir3=gridshift3+1,nr3+1
         ii3=mod(ir3-1,nr3) + 1
         do ir2=gridshift2+1,nr2+1
           ii2=mod(ir2-1,nr2) + 1
           do ir1=gridshift1+1,nr1+1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
           do ir1=1,gridshift1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
         end do
         do ir2=1,gridshift2
           ii2=mod(ir2-1,nr2) + 1
           do ir1=gridshift1+1,nr1+1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
           do ir1=1,gridshift1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
         end do
       end do
       do ir3=1,gridshift3
         ii3=mod(ir3-1,nr3) + 1
         do ir2=gridshift2+1,nr2+1
           ii2=mod(ir2-1,nr2) + 1
           do ir1=gridshift1+1,nr1+1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
           do ir1=1,gridshift1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
         end do
         do ir2=1,gridshift2
           ii2=mod(ir2-1,nr2) + 1
           do ir1=gridshift1+1,nr1+1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
           do ir1=1,gridshift1
             ii1=mod(ir1-1,nr1) + 1
             tmpr=fofr(1,ii1,ii2,ii3)
             tmpi=fofr(2,ii1,ii2,ii3)
             write(unout,'(e12.5)') tmpr*tmpr + tmpi*tmpi
           end do
         end do
       end do


       write(unout,*)
       write(unout,'(1X,A)') 'END_DATAGRID_3D'
       write(unout,'(1X,A)') 'END_BLOCK_DATAGRID3D'
       close(unout)

       ABI_DEALLOCATE(filename)

       write(std_out,*)
       exit


     case(12)
       write(std_out,*)"NetCDF output is not available anymore"
       exit

!        ************************************************************

     case(13)
       write(std_out,*)
       write(std_out,*) 'Give 1 files of formatted data'
       write(std_out,*) 'The files are ready to be used with XCrysDen'
       write(std_out,*)
       gridshift1 = 0
       gridshift2 = 0
       gridshift3 = 0
       write(std_out,*) 'Do you want to shift the grid along the x,y or z axis (y/n)?'
       write(std_out,*)
       shift_tau(:) = 0.0
       if (read_string(outputchar, unit=std_in) /= 0) then
         MSG_ERROR("Fatal error!")
       end if
       if (outputchar == 'y' .or. outputchar == 'Y') then
         MSG_ERROR("Shift is buggy, don't use it")
         write(std_out,*) 'Give the three shifts (x,y,z < ',nr1,nr2,nr3,') :'
         write(std_out,*)
         read (std_in,*) gridshift1, gridshift2, gridshift3
         shift_tau(:) = gridshift1*rprimd(:,1)/(nr1+1) + gridshift2*rprimd(:,2)/(nr2+1) + gridshift3*rprimd(:,3)/(nr3+1)
       end if

       ABI_ALLOCATE(filename,(1))
       filename(1)=trim(output)
       write(std_out,*) '  The name of your data files is : '
       write(std_out,*) trim(filename(1)),'  for the density (norm of the wfk),'
       write(std_out,*)

       if (open_file(filename(1), msg, newunit=unout, status='unknown',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       rewind(unout)

       do iband=1,nband(ckpt)
         write(unout,'(a,2f20.16)')'#', eig_k(iband),occ_k(iband)
       end do

       write(std_out,'(/,a,2x,3i5)' )' Number of points per side: ',nr1+1,nr2+1,nr3+1
       write(std_out,'(/,a,2x,i10,//)' )' Total number of points:', (nr1+1)*(nr2+1)*(nr3+1)
       write(std_out,*) ' znucl = ', znucl, ' typat = ', typat, ' ntypat = ', ntypat

       write(unout,'(1X,A)')  'DIM-GROUP'
       write(unout,*) '3  1'
       write(unout,'(1X,A)') 'PRIMVEC'
       do ir1 = 1,3
         write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
       end do
       write(unout,'(1X,A)') 'PRIMCOORD'
       write(unout,*) natom, ' 1'
!
!        generate translated coordinates to match density shift
!
       do iatom = 1,natom
         tau2 (:,iatom) = xcart(:,iatom) - shift_tau(:)
       end do

       do iatom = 1,natom
         write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&         Bohr_Ang*tau2(1,iatom), &
&         Bohr_Ang*tau2(2,iatom), &
&         Bohr_Ang*tau2(3,iatom)
       end do
       write(unout,'(1X,A)') 'ATOMS'
       do iatom = 1,natom
         write(unout,'(i9,3(3X,ES17.10))') int(znucl(typat(iatom))), &
&         Bohr_Ang*tau2(1,iatom), &
&         Bohr_Ang*tau2(2,iatom), &
&         Bohr_Ang*tau2(3,iatom)
       end do

!        write(unout,'(1X,A)') 'FRAMES'
       write(unout,'(1X,A)') 'BEGIN_BLOCK_DATAGRID3D'
       write(unout,*) 'datagrids'
       write(unout,'(1X,A)') 'DATAGRID_3D_DENSITY'
       write(unout,*) nr1+1,nr2+1,nr3+1
       write(unout,*) '0.0 0.0 0.0 '
       do ir1 = 1,3
         write(unout,'(3(ES17.10,2X))') (Bohr_Ang*rprimd(ir2,ir1), ir2=1,3)
       end do

       do ir3=1,nr3+1
         ii3=mod(ir3-1+gridshift3, nr3) + 1
         do ir2=1,nr2+1
           ii2=mod(ir2-1+gridshift2, nr2) + 1
           do ir1=1,nr1+1
             ii1=mod(ir1-1+gridshift1, nr1) + 1
             write(unout,'(ES17.10)') fofr(1,ii1,ii2,ii3)
           end do
         end do
       end do
       write(unout,*)
       write(unout,'(1X,A)') 'END_DATAGRID_3D'
       write(unout,'(1X,A)') 'END_BLOCK_DATAGRID3D'
       close(unout)

       ABI_DEALLOCATE(filename)

       write(std_out,*)
       exit

     case(14) ! CUBE file format from GAUSSIAN

       write(std_out,*)
       write(std_out,*) 'Output a cube file of 3D volumetric data'
       write(std_out,*)

       if (open_file(output, msg, newunit=unout, status='replace',form='formatted') /= 0) then
         MSG_ERROR(msg)
       end if
       call print_fofr_cube(nr1,nr2,nr3,n4,n5,n6,fofr,rprimd,natom,znucl_atom_int,xcart,unit=unout)
       close(unout)
       exit

     case(0)
       write(std_out,*)' Exit inner loop'
       select_exit = 1

     case default
       write(std_out,*) ' This choice is not valid.'
       write(std_out,*)
       cycle

     end select

   end do

   ckpt=oldckpt
   cband=oldcband
   csppol=oldcsppol
   cspinor=oldcspinor
!  deallocate the datas
   ABI_DEALLOCATE(fofr)

   write(std_out,*) ' Task ',ichoice,' has been done !'
   write(std_out,*)
   write(std_out,*) ' Run interpolation again? (1=default=yes,0=no)'
   read(std_in,*) iprompt
   if(iprompt==0) then
     exit
   else
     cycle
   end if
 end do

!Deallocate the datas
 ABI_DEALLOCATE(cg_k)
 ABI_DEALLOCATE(eig_k)
 ABI_DEALLOCATE(kg_dum)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(occ_k)

 call destroy_mpi_enreg(mpi_enreg)

end subroutine cut3d_wffile
!!***

!----------------------------------------------------------------------

end module  m_cut3d
