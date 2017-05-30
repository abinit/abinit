!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pptools
!! NAME
!! m_pptools
!!
!! FUNCTION
!!  Helper functions used for simple post-processing.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (MG,ZL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_pptools

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 private

 public :: prmat                ! print real(dp) matrices in an attractive format.
 public :: printxsf             ! Write a generic array in the XSF format (XCrysden format)
 public :: print_fofr_ri        ! Print the [real, imaginary] part of an array
 public :: print_fofr_xyzri     ! Print the Cartesian coordinates and the [real,imaginary] part of an array
 public :: print_fofr_cube      ! Print ||fofr|| in the CUBE format.


CONTAINS  !===========================================================
!!***

!!****f* m_pptools/prmat
!! NAME
!! prmat
!!
!! FUNCTION
!! This subroutine prints real*8 matrices in an attractive format.
!!
!! INPUTS
!!  mat(mi,nj)= matrix to be printed
!!  mi        = no rows of mat
!!  ni        = no rows to print
!!  nj        = no colums of mat
!!  unitm     = unit to print to, if not provided std_out is chosen
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      chiscwrt,ioniondist,newkpt,pawuj_det,pawuj_red,pawuj_utils,shellstruct
!!
!! CHILDREN
!!
!! SOURCE

subroutine prmat (mat, ni, nj, mi, unitm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prmat'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)           :: mi,ni,nj
 integer,intent(in), optional :: unitm
!arrays
 real(dp),intent(in)          :: mat(mi,nj)

!Local variables-------------------------------
!scalars
 character(len=1000)    :: message
 integer,parameter      :: nline=10
 integer                :: ii,jj,jstart,jstop,unitn

! *************************************************************************

 if (present(unitm)) then ! standard printing to std_out
   unitn=unitm
 else
   unitn=std_out
 end if

 do  jstart = 1, nj, nline
   jstop = min(nj, jstart+nline-1)
   write(message, '(3x,10(i4,8x))' ) (jj,jj=jstart,jstop)
   call wrtout(unitn,message,'COLL') 
 end do

 do ii = 1,ni
   do jstart= 1, nj, nline
     jstop = min(nj, jstart+nline-1)
     if (jstart==1) then
       write(message, '(i3,1p,10e12.4)' ) ii, (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     else
       write(message, '(3x,1p,10e12.4)' )    (mat(ii,jj),jj=jstart,jstop)
       call wrtout(unitn,message,'COLL')
     end if
   end do
 end do

end subroutine prmat
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/printxsf
!! NAME
!! printxsf
!!
!! FUNCTION
!! Write a generic array in the XSF format (XCrysden format)
!!
!! INPUTS
!! basis(3,3) = basis vectors of the direct real lattice or of the reciprocal lattice (fortran convention)
!!              (Bohr units if realrecip=0, Bohr^-1 if realrecip=1, see below)
!! realrecip = 0  for a plot in real space
!!             1  for a plot in reciprocal space
!! nunit   = unit number of the output file (already open by the caller, not closed here!)
!! n1=grid size along x
!! n2=grid size along y
!! n3=grid size along z
!! origin(3) = origin of the grid
!! datagrid(n1*n2*n3) = datagrid values stored using the fortran convention
!!
!! OUTPUT
!! Only write
!!
!! PARENTS
!!      denfgr,exc_plot,m_kxc,m_nesting,m_wfd,spin_current
!!
!! CHILDREN
!!
!! SOURCE

subroutine printxsf(n1,n2,n3,datagrid,basis,origin,natom,ntypat,typat,xcart,znucl,nunit,realrecip)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'printxsf'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nunit,realrecip
 integer,intent(in) :: natom,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: basis(3,3),datagrid(n1*n2*n3),origin(3)
 real(dp),intent(in) :: xcart(3,natom), znucl(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,nslice,nsym,iatom
 real(dp) :: fact
 character(len=500) :: msg
!arrays
 real(dp) :: tau(3,natom)

! *************************************************************************

 DBG_ENTER("COLL")

 if ( ALL(realrecip/=(/0,1/)) )then
   write(msg,'(a,i6)')' The argument realrecip should be 0 or 1; however, realrecip= ',realrecip
   MSG_BUG(msg)
 end if

!conversion between ABINIT default units and XCrysden units
 fact=Bohr_Ang; if (realrecip ==1) fact=one/fact  !since we are in reciprocal space

!TODO insert crystalline structure and dummy atoms in case of reciprocal space
!need to convert basis too

 write(nunit,'(1X,A)')  'DIM-GROUP'
 write(nunit,*) '3  1'
 write(nunit,'(1X,A)') 'PRIMVEC'
 do iy = 1,3
   write(nunit,'(3(ES17.10,2X))') (Bohr_Ang*basis(ix,iy), ix=1,3)
 end do
!
!generate translated coordinates to fit origin shift
!
 do iatom = 1,natom
   tau (:,iatom) = xcart(:,iatom) - origin(:)
 end do

 write(nunit,'(1X,A)') 'PRIMCOORD'
 write(nunit,*) natom, ' 1'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), &  ! WARNING alchemy not supported by XCrysden
&  Bohr_Ang*tau(1,iatom), &
&  Bohr_Ang*tau(2,iatom), &
&  Bohr_Ang*tau(3,iatom)
 end do
 write(nunit,'(1X,A)') 'ATOMS'
 do iatom = 1,natom
   write(nunit,'(i9,3(3X,ES17.10))') NINT(znucl(typat(iatom))), & ! WARNING alchemy not supported by XCrysden
&  Bohr_Ang*tau(1,iatom), &
&  Bohr_Ang*tau(2,iatom), &
&  Bohr_Ang*tau(3,iatom)
 end do

 write(nunit,'(a)')' BEGIN_BLOCK_DATAGRID3D'
 write(nunit,'(a)')' datagrid'
 write(nunit,'(a)')' DATAGRID_3D_DENSITY'
!NOTE: XCrysden uses aperiodical data grid
 write(nunit,*)n1+1,n2+1,n3+1
 write(nunit,*)origin
 write(nunit,*)basis(:,1)*fact
 write(nunit,*)basis(:,2)*fact
 write(nunit,*)basis(:,3)*fact

 nslice=1
 do iz=1,n3
   do iy=1,n2
     write(nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
     nslice=nslice+1
   end do
   nsym=nslice-n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
 end do

!Now write upper plane
 nslice=1
 do iy=1,n2
   write (nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
   nslice=nslice+1
 end do

 nsym=nslice-n2
 write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))

 write (nunit,'(a)')' END_DATAGRID_3D'
 write (nunit,'(a)')' END_BLOCK_DATAGRID3D'

 DBG_EXIT("COLL")

end subroutine printxsf
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_ri
!! NAME
!!  print_fofr_ri
!!
!! FUNCTION
!!  Print the [real,imaginary] part of fofr on unit unit
!!
!! INPUTS
!!  ri_mode = 
!!    "RI" if both real and imag part are wanted
!!    "R"  for real part
!!    "I"  for imaginary part
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE


subroutine print_fofr_ri(ri_mode,nx,ny,nz,ldx,ldy,ldz,fofr,unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_fofr_ri'
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz
 integer,optional,intent(in) :: unit
 character(len=*),intent(in) :: ri_mode
!arrays
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer :: ount,ix,iy,iz
!arrays

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 SELECT CASE (ri_mode)
 CASE ("RI","ri")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(2f20.16)') fofr(:,ix,iy,iz)
       end do
     end do
   end do

 CASE ("R","r")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(f20.16)') fofr(1,ix,iy,iz)
       end do
     end do
   end do

 CASE ("I","i")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         write(ount,'(f20.16)') fofr(2,ix,iy,iz)
       end do
     end do
   end do

 CASE DEFAULT
   MSG_ERROR("Wrong ri_mode")
 END SELECT

end subroutine print_fofr_ri
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_xyzri
!! NAME
!!  print_fofr_xyzri
!!
!! FUNCTION
!!  Print the Cartesian coordinates and the [real,imaginary] part of fofr on unit unit
!!
!! INPUTS
!!  ri_mode = 
!!    "RI" if both real and imag part are wanted
!!    "R"  for real part
!!    "I"  for imaginary part
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  rprimd(3,3)=Lattive vectors in Bohr
!!  [conv_fact] = Conversion factor for rprimd (rprimd is multiplied by conv_fact). Default is one
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE


subroutine print_fofr_xyzri(ri_mode,nx,ny,nz,ldx,ldy,ldz,fofr,rprimd,conv_fact,unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_fofr_xyzri'
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz
 integer,optional,intent(in) :: unit
 real(dp),optional,intent(in) :: conv_fact
 character(len=*),intent(in) :: ri_mode
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer :: ount,ix,iy,iz
 real(dp) :: xnow,ynow,znow,my_cfact
!arrays

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit
 my_cfact = one; if (PRESENT(conv_fact)) my_cfact = conv_fact

 SELECT CASE (ri_mode)
 CASE ("RI","ri")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,2f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(:,ix,iy,iz)
       end do
     end do
   end do

 CASE ("R","r")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(1,ix,iy,iz)
       end do
     end do
   end do

 CASE ("I","i")
   do iz=1,nz
     do iy=1,ny
       do ix=1,nx
         xnow = rprimd(1,1)*(ix-1)/nx + rprimd(1,2)*(iy-1)/ny + rprimd(1,3)*(iz-1)/nz
         ynow = rprimd(2,1)*(ix-1)/nx + rprimd(2,2)*(iy-1)/ny + rprimd(2,3)*(iz-1)/nz
         znow = rprimd(3,1)*(ix-1)/nx + rprimd(3,2)*(iy-1)/ny + rprimd(3,3)*(iz-1)/nz
         write(ount,'(3f16.10,f20.16)') my_cfact*xnow, my_cfact*ynow, my_cfact*znow,fofr(2,ix,iy,iz)
       end do
     end do
   end do

 CASE DEFAULT
   MSG_ERROR("Wrong ri_mode")
 END SELECT

end subroutine print_fofr_xyzri
!!***

!----------------------------------------------------------------------

!!****f* m_pptools/print_fofr_cube
!! NAME
!!  print_fofr_cube
!!
!! FUNCTION
!!  Print array fofr in the cube file format 
!!
!! INPUTS
!!  nx,ny,nz,ldx,ldy,ldz = Logical and physical dimensions of the array.
!!  fofr(2,ldx,ldy,ldz) = Input data
!!  rprimd(3,3)=Lattive vectors in Bohr
!!  [unit] = Fortran unit number. Default: std_out
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!
!! SOURCE


subroutine print_fofr_cube(nx,ny,nz,ldx,ldy,ldz,fofr,rprimd,natom,znucl_atom,xcart,unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'print_fofr_cube'
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,natom
 integer,optional,intent(in) :: unit
!arrays
 integer,intent(in) :: znucl_atom(natom)
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(in) :: fofr(2,ldx,ldy,ldz)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplx=2
 integer :: ount,ix,iy,iz,iatom
!arrays

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

! EXAMPLE FROM THE WEB
! CPMD CUBE FILE.
! OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
! 3    0.000000    0.000000    0.000000
! 40    0.283459    0.000000    0.000000
! 40    0.000000    0.283459    0.000000
! 40    0.000000    0.000000    0.283459
! 8    0.000000    5.570575    5.669178    5.593517
! 1    0.000000    5.562867    5.669178    7.428055
! 1    0.000000    7.340606    5.669178    5.111259
! -0.25568E-04  0.59213E-05  0.81068E-05  0.10868E-04  0.11313E-04  0.35999E-05

 write(ount,'(a)') 'ABINIT generated cube file'
 write(ount,'(a)') 'from cut3d tool'

 write(ount,'(i9,3(1x,f12.6))') natom,0.,0.,0.
 write(ount,'(i9,3(1x,f12.6))') nx,(rprimd(iy,1)/nx, iy=1,3)
 write(ount,'(i9,3(1x,f12.6))') ny,(rprimd(iy,2)/ny, iy=1,3)
 write(ount,'(i9,3(1x,f12.6))') nz,(rprimd(iy,3)/nz, iy=1,3)

 do iatom=1,natom
   write(ount,'(i9,4(3X,ES17.10))') znucl_atom(iatom),0.d0,xcart(1:3,iatom)
 end do

! Note C ordering of the indexes 
 if (cplx==2) then
   do ix=1,nx
     do iy=1,ny
       do iz=1,nz
         write(ount,'(6(f12.6,2x))') sqrt(fofr(1,ix,iy,iz)**2 + fofr(2,ix,iy,iz)**2 )
       end do
     end do
   end do
 else
   do ix=1,nx
     do iy=1,ny
       do iz=1,nz
         write(ount,'(6(f12.6,2x))') fofr(1,ix,iy,iz)
       end do
     end do
   end do
 end if

end subroutine print_fofr_cube
!!***

!----------------------------------------------------------------------

END MODULE m_pptools
