!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtlammps
!! NAME
!!  prtlammps
!!
!! FUNCTION
!!  Output LAMMPS MD style dump format for use with visualizers and post-processing. Assumes fixed dump atom format as follows:
!!  Atom-ID Atom-Type x(uw/us) y(uw/us) z(uw/us) fx fy fz
!!  (uw/us) - unwrapped and unscaled coordinates in periodic cell.
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (Stefan Bringuier)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! AUTHOR(S): Stefan Bringuier
!! CONTACT(S): stefanbringuier@gmail.com
!!
!! INPUTS
!!  fnameradix = prefix used for the output file (e.g. name of simulation)
!!  natom = number of atoms
!!  ntypat = number of types of atoms
!!  rprimd = lattice vectors for the primitive cell
!!  typat = type for each of the natom atoms
!!  istep = current ion timestep
!!  xcart = cartesian positions of the atoms
!!  fcart = forces on atoms in cartesian coordinates
!!  
!!
!! OUTPUT
!!  No return value, only writes output files.
!!
!! SIDE EFFECTS
!!
!! NOTES
!! The dump format here is fixed unlike the LAMMPS dump command which is customizable
!!
!! TODO
!! If needed it it is straighforward to add velocities.
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


subroutine prtlammps(fnameradix, natom, ntypat, rprimd, typat, istep, xred, fcart)

 use defs_basis
 use m_profiling_abi
 use m_atomdata
 use m_errors

 use m_io_tools,     only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtlammps'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, ntypat
 integer, intent(in) :: istep
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: fcart(3,natom)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,natom)
 character(len=fnlen), intent(in) :: fnameradix

!Local variables-------------------------------
!scalars
 integer :: iatom, itypat, iout
 type(atomdata_t) :: atom
 !Unit cell/box variables
 real(dp) :: magr1,magr2,magr3
 real(dp) :: cosalpha,cosbeta,cosgamma
 real(dp) :: ly,lz,xy,xz,yz
 real(dp) :: xlo_bound,xhi_bound
 real(dp) :: ylo_bound,yhi_bound
 real(dp) :: zlo_bound,zhi_bound
! arrays
 character(len=500) :: msg
 character(len=fnlen) :: fname
 real(dp) :: xcart(3,natom)
 character(len=4) :: sistep
!************************************************************************

 !!Open output file
 write(sistep,'(I4.4)') istep !! Convert istep to string
 fname = trim(fnameradix)//".LAMMPS.dump."//trim(sistep)
 if (open_file(fname,msg,newunit=iout) /= 0) then
   MSG_ERROR(msg)
 end if


 !!Convert cell to LAMMPS cell dump format. For details see http://lammps.sandia.gov/doc/Section_howto.html#howto-12
 !!Norm of cell vectors
 magr1=rprimd(1,1)*rprimd(1,1)+rprimd(1,2)*rprimd(1,2)+rprimd(1,3)*rprimd(1,3)
 magr2=rprimd(2,1)*rprimd(2,1)+rprimd(2,2)*rprimd(2,2)+rprimd(2,3)*rprimd(2,3)
 magr3=rprimd(3,1)*rprimd(3,1)+rprimd(3,2)*rprimd(3,2)+rprimd(3,3)*rprimd(3,3)
 
 !!Cosine angle terms
 cosgamma=(rprimd(1,1)*rprimd(1,2)+rprimd(2,1)*rprimd(2,2)+rprimd(3,1)*rprimd(3,2))/(magr1*magr2)
 cosbeta=(rprimd(1,2)*rprimd(1,3)+rprimd(2,2)*rprimd(2,3)+rprimd(3,2)*rprimd(3,3))/(magr2*magr3)
 cosalpha=(rprimd(1,1)*rprimd(1,3)+rprimd(2,1)*rprimd(2,3)+rprimd(3,1)*rprimd(3,3))/(magr1*magr3)
 
 !!Cell box lengths and skew
 xy=rprimd(2,2)*cosgamma
 xz=rprimd(3,3)*cosbeta
 ly=sqrt(rprimd(2,2)*rprimd(2,2)-xy*xy)
 yz=(rprimd(2,2)*rprimd(3,3)*cosalpha-xy*xz)/ly
 lz=sqrt(rprimd(3,3)*rprimd(3,3)-xz*xz-yz*yz)

 !!Set the box bounds for LAMMPS dump format
 xlo_bound=MIN(0.0,xy,xz,xy+xz)
 xhi_bound=rprimd(1,1)+MAX(0.0,xy,xz,xy+xz)
 ylo_bound=MIN(00.0,yz)
 yhi_bound=ly+MAX(0.0,xy,xz,xy+xz)
 zlo_bound=0.00
 zhi_bound=lz

 !!LAMMPS dumpfile header content
 write (iout,'(a)') "ITEM: TIMESTEP"
 write (iout,'(I7)') istep
 write (iout,'(a)') "ITEM: NUMBER OF ATOMS"
 write (iout,'(I7)') natom
 write (iout,'(a)') "ITEM: BOX BOUNDS pp pp pp xy xz yz xx yy zz"
 write (iout,'(3E24.14,1x)') Bohr_Ang*xlo_bound,Bohr_Ang*xhi_bound,xy
 write (iout,'(3E24.14,1x)') Bohr_Ang*ylo_bound,Bohr_Ang*yhi_bound,xz
 write (iout,'(3E24.14,1x)') Bohr_Ang*zlo_bound,Bohr_Ang*zhi_bound,yz
 write (iout,'(a)') "ITEM: ATOMS id type xu yu zu fx fy fz"

 !!convert xred to xcart
 call xred2xcart(natom,rprimd,xcart,xred)
 !!Output atomic data in Angstroms and eV (metal units in LAMMPS). Assumes unwrapped and unscaled coordinates.
 do iatom=1,natom
    write(iout,'(I7,I7,6(E24.14,1x))') iatom,typat(iatom),Bohr_Ang*xcart(:,iatom), Ha_eV/Bohr_Ang*fcart(:,iatom) 
 end do
 close(iout)


end subroutine prtlammps
!!***
