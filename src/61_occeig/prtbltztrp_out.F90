! {\src2tex{textfont=tt}}
!!****f* ABINIT/prtbltztrp_out
!! NAME
!! prtbltztrp_out
!!
!! FUNCTION
!!   output files for BoltzTraP code, which integrates Boltzmann transport quantities
!!   over the Fermi surface for different T and chemical potentials. Abinit provides
!!   all necessary input files: struct, energy, input file, and def file for the unit
!!   definitions of fortran files in BT.
!!   See http://www.icams.de/content/departments/ams/madsen/boltztrap.html
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol) = array for holding eigenvalues (hartree)
!!  fermie = Fermi level
!!  fname_radix = radix of file names for output
!!  natom = number of atoms in cell.
!!  nband = number of bands
!!  nkpt = number of k points.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetries in space group
!!  rprimd(3,3) = dimensional primitive translations for real space (bohr)
!!  symrel = symmetry operations in reduced coordinates, real space
!!  to be used in future  xred(3,natom) = reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! PARENTS
!!      eph,outscfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtbltztrp_out (eigen, fermie, fname_radix, kpt, &
&       natom, nband, nelec, nkpt, nspinor, nsppol, nsym, &
&       rprimd, symrel, tau_k)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools,     only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtbltztrp_out'
!End of the abilint section

 implicit none

! arguments
 integer, intent(in) :: natom, nsym, nband, nkpt, nsppol, nspinor
 real(dp), intent(in) :: fermie
 real(dp), intent(in) :: nelec(nsppol)

 integer, intent(in) :: symrel(3,3,nsym)
 real(dp), intent(in) :: kpt(3,nkpt)
 real(dp), intent(in) :: eigen(nband, nkpt, nsppol)
 character(len=fnlen), intent(in) :: fname_radix
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in), optional :: tau_k(nsppol,nkpt,nband)

!local
 integer :: iout, isym, iband, isppol, ikpt
 character(len=fnlen) :: filename
 character(len=3) :: spinsuffix(nsppol)
 character(len=2) :: so_suffix
 character(len=500) :: msg
 
 real(dp) :: ha2ryd

! *************************************************************************

 ha2ryd = two

!bogus line using natom - we may include it in the future
 iout = natom

 so_suffix=""
 if (nsppol > 1 .or. nspinor > 1) so_suffix="so"

 if (nsppol == 1) then
   spinsuffix(1) = "ns_"
 else 
   spinsuffix = (/"up_", "dn_"/)
 end if

 do isppol = 1, nsppol
!input file for boltztrap: general info, Ef, Nelec, etc...

   filename= trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.intrans"
   if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(a)') "GENE                      # Format of input: generic format, with Symmetries"
   write (iout, '(a)') "0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap"
   write (iout, '(E15.5,a,F10.4,a)') fermie*two, " 0.0005 0.4  ", nelec(isppol), &
&   "  # Fermilevel (Ry), energy grid spacing, energy span around Fermilevel, number of electrons for this spin"
   write (iout, '(a)') "CALC                      # CALC (calculate expansion coeff), NOCALC read from file"
   write (iout, '(a)') "3                         # lpfac, number of latt-points per k-point"
   write (iout, '(a)') "BOLTZ                     # run mode (only BOLTZ is supported)"
   write (iout, '(a)') ".15                       # (efcut) energy range of chemical potential"
   write (iout, '(a)') "300. 10.                  # Tmax, temperature grid spacing"
   write (iout, '(2a)') "-1                        # energyrange of bands given ",&
&   "individual DOS output sig_xxx and dos_xxx (xxx is band number)"
   write (iout, '(a)') "HISTO                     # DOS calculation method. Other possibility is TETRA"
   write (iout, '(a)') "No                        # not using model for relaxation time"
   write (iout, '(a)') "3                         # Number of doping levels coefficients will be output for"
   write (iout, '(a)') "-1.e16 0.0d0 1.e16        # Values of doping levels (in carriers / cm^3" 
   
   close(iout)

!files file, with association of all units for Boltztrap
   filename= trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.def"
   if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(3a)') "5, '", trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.intrans',      'old',    'formatted',0"
   write (iout, '(3a)') "6, '", trim(fname_radix)//"_BLZTRP", ".outputtrans',      'unknown',    'formatted',0"
   write (iout, '(3a)') "20,'", trim(fname_radix)//"_BLZTRP", ".struct',         'old',    'formatted',0"
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP."//trim(spinsuffix(isppol))//"energy"//trim(so_suffix),&
&   "',         'old',    'formatted',0"
   if (present (tau_k)) then
     write (iout, '(3a)') "11,'", trim(fname_radix)//"_BLZTRP", ".tau_k',         'old',    'formatted',0"
   end if
   write (iout, '(3a)') "48,'", trim(fname_radix)//"_BLZTRP", ".engre',         'unknown',    'unformatted',0"
   write (iout, '(3a)') "49,'", trim(fname_radix)//"_BLZTRP", ".transdos',        'unknown',    'formatted',0"
   write (iout, '(3a)') "50,'", trim(fname_radix)//"_BLZTRP", ".sigxx',        'unknown',    'formatted',0"
   write (iout, '(3a)') "51,'", trim(fname_radix)//"_BLZTRP", ".sigxxx',        'unknown',    'formatted',0"
   write (iout, '(3a)') "21,'", trim(fname_radix)//"_BLZTRP", ".trace',           'unknown',    'formatted',0"
   write (iout, '(3a)') "22,'", trim(fname_radix)//"_BLZTRP", ".condtens',           'unknown',    'formatted',0"
   write (iout, '(3a)') "24,'", trim(fname_radix)//"_BLZTRP", ".halltens',           'unknown',    'formatted',0"
   write (iout, '(3a)') "25,'", trim(fname_radix)//"_BLZTRP", ".trace_fixdoping',     'unknown',    'formatted',0"
   write (iout, '(3a)') "26,'", trim(fname_radix)//"_BLZTRP", ".condtens_fixdoping',           'unknown',    'formatted',0"
   write (iout, '(3a)') "27,'", trim(fname_radix)//"_BLZTRP", ".halltens_fixdoping',           'unknown',    'formatted',0"
   write (iout, '(3a)') "30,'", trim(fname_radix)//"_BLZTRP", "_BZ.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "31,'", trim(fname_radix)//"_BLZTRP", "_fermi.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "32,'", trim(fname_radix)//"_BLZTRP", "_sigxx.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "33,'", trim(fname_radix)//"_BLZTRP", "_sigyy.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "34,'", trim(fname_radix)//"_BLZTRP", "_sigzz.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "35,'", trim(fname_radix)//"_BLZTRP", "_band.dat',           'unknown',    'formatted',0"
   write (iout, '(3a)') "36,'", trim(fname_radix)//"_BLZTRP", "_band.gpl',           'unknown',    'formatted',0"
   write (iout, '(3a)') "37,'", trim(fname_radix)//"_BLZTRP", "_deriv.dat',           'unknown',    'formatted',0"
   write (iout, '(3a)') "38,'", trim(fname_radix)//"_BLZTRP", "_mass.dat',           'unknown',    'formatted',0"
   
   close(iout)
 end do !isppol

!file is for geometry symmetries etc
 filename= trim(fname_radix)//"_BLZTRP.struct"
 if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 write (iout, '(a)') "BoltzTraP geometry file generated by ABINIT."

!here we need to print out the unit cell vectors 
 write (iout, '(3E20.10)') rprimd(:,1)
 write (iout, '(3E20.10)') rprimd(:,2)
 write (iout, '(3E20.10)') rprimd(:,3)

 write (iout, '(I7)') nsym

 do isym=1, nsym
   write (iout,'(3(3I5,2x), a, I5)') &
&   symrel(1,:,isym), &
&   symrel(2,:,isym), &
&   symrel(3,:,isym), &
&   ' ! symmetry rotation matrix isym = ', isym
 end do

 close (iout)

!second file is for eigenvalues

 do isppol = 1, nsppol 
! two file names for each spin, if necessary
   filename=trim(fname_radix)//"_BLZTRP."//spinsuffix(isppol)//"energy"//trim(so_suffix)
   
   if (open_file (filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(a,I5)') "BoltzTraP eigen-energies file generated by ABINIT. ispin = ", isppol
   write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie, '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'
   do ikpt = 1, nkpt
!    these need to be in reduced coordinates
     write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
     do iband = 1, nband
!      output in eV
       write (iout, '(E20.10)') ha2ryd*eigen(iband, ikpt, isppol)
     end do
   end do

   close (iout)
 end do

!this file is for tau_k
 if (present (tau_k)) then
   do isppol = 1, nsppol 
     filename= trim(fname_radix)//"_"//spinsuffix(isppol)//"BLZTRP.tau_k"
     if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if

     write (iout, '(a)') "BoltzTraP tau_k file generated by ANADDB."
     write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie, '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'

     do ikpt = 1, nkpt
!      these need to be in reduced coordinates
       write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
       do iband = 1, nband
!        output in eV
         write (iout, '(E20.10)') tau_k(isppol,ikpt,iband)
       end do
     end do

     close (iout)
   end do

 end if

end subroutine prtbltztrp_out
!!***
