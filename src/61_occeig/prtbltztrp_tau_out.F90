! {\src2tex{textfont=tt}}
!!****f* ABINIT/prtbltztrp_tau_out
!! NAME
!! prtbltztrp_tau_out
!!
!! FUNCTION
!!   output files for BoltzTraP code, which integrates Boltzmann transport quantities
!!   over the Fermi surface for different T and chemical potentials. Abinit provides
!!   all necessary input files: struct, energy, input file, and def file for the unit
!!   definitions of fortran files in BT.
!!   See http://www.icams.de/content/departments/ams/madsen/boltztrap.html
!!   Output T-depedent tau_k, modified from prtbltztrp_out
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (BXu)
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
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      get_tau_k
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtbltztrp_tau_out (eigen, tempermin, temperinc, ntemper, fermie, fname_radix, kpt, &
&       natom, nband, nelec, nkpt, nspinor, nsppol, nsym, &
&       rprimd, symrel, tau_k)

 use defs_basis
 use m_profiling_abi
 use m_io_tools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtbltztrp_tau_out'
!End of the abilint section

 implicit none

! arguments
 integer, intent(in) :: natom, nsym, nband, nkpt, nsppol, nspinor, ntemper
 real(dp), intent(in) :: tempermin, temperinc
 real(dp), intent(in) :: fermie(ntemper)
 real(dp), intent(in) :: nelec

 integer, intent(in) :: symrel(3,3,nsym)
 real(dp), intent(in) :: kpt(3,nkpt)
 real(dp), intent(in) :: eigen(nband, nkpt, nsppol)
 character(len=fnlen), intent(in) :: fname_radix
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: tau_k(ntemper,nsppol,nkpt,nband)

!local
 integer :: iout, isym, iband, isppol, ikpt, itemp
 character(len=fnlen) :: filename
 character(len=fnlen) :: appendix
 real(dp) :: ha2ryd, Temp

! *************************************************************************

 ha2ryd = two

!bogus line using natom - we may include it in the future
 iout = natom

 iout = get_unit()

!input file for boltztrap: general info, Ef, Nelec, etc...

 do itemp = 1, ntemper
   write(appendix,"(i0)") itemp
   filename= trim(fname_radix)//"_BLZTRP.intrans_"//trim(appendix)
   open (iout, file=filename, form='formatted')
   
   write (iout, '(a)') "GENE                      # Format of input: generic format, with Symmetries"
   write (iout, '(a)') "0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap"
   write (iout, '(E15.5,a,F10.4,a)') fermie(itemp)*two, " 0.0005 0.4  ", nelec, &
&   "  # Fermilevel (Ry), energy grid spacing, energy span around Fermilevel, number of electrons"
   write (iout, '(a)') "CALC                      # CALC (calculate expansion coeff), NOCALC read from file"
   write (iout, '(a)') "3                         # lpfac, number of latt-points per k-point"
   write (iout, '(a)') "BOLTZ                     # run mode (only BOLTZ is supported)"
   write (iout, '(a)') ".15                       # (efcut) energy range of chemical potential"
   write (iout, '(2f8.2,a)')&
&   tempermin+temperinc*dble(itemp),tempermin+temperinc*dble(itemp), "                  # Tmax, temperature grid spacing"
   write (iout, '(2a)') "-1                        # energyrange of bands given ",&
&   "individual DOS output sig_xxx and dos_xxx (xxx is band number)"
   write (iout, '(a)') "TETRA                     # DOS calculation method. Other possibility is TETRA"
   write (iout, '(a)') "No                        # not using model for relaxation time"
   write (iout, '(a)') "3                         # Number of doping levels coefficients will be output for"
   write (iout, '(a)') "-1.e16 0.0d0 1.e16        # Values of doping levels (in carriers / cm^3" 
   
   close(iout)
 end do

!files file, with association of all units for Boltztrap
 filename= trim(fname_radix)//"_BLZTRP.def"
 open (iout, file=filename, form='formatted')

 write (iout, '(3a)') "5, '", trim(fname_radix)//"_BLZTRP", ".intrans',      'old',    'formatted',0"
 write (iout, '(3a)') "6, '", trim(fname_radix)//"_BLZTRP", ".outputtrans',      'unknown',    'formatted',0"
 write (iout, '(3a)') "20,'", trim(fname_radix)//"_BLZTRP", ".struct',         'old',    'formatted',0"
 if (nspinor == 1) then
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energy',         'old',    'formatted',0"
 else if (nspinor == 2) then
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energyso',         'old',    'formatted',0"
 end if
 write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energy',         'old',    'formatted',0"
 write (iout, '(3a)') "11,'", trim(fname_radix)//"_BLZTRP", ".tau_k',         'old',    'formatted',0"
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

!file is for geometry symmetries etc
 filename= trim(fname_radix)//"_BLZTRP.struct"
 open (iout, file=filename, form='formatted')

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
 if (nspinor == 1) then
   filename= trim(fname_radix)//"_BLZTRP.energy"
 else if (nspinor == 2) then
   filename= trim(fname_radix)//"_BLZTRP.energyso"
 end if
 open (iout, file=filename, form='formatted')

 write (iout, '(a)') "BoltzTraP eigen-energies file generated by ABINIT."
 write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie(1), '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'
 do isppol = 1, nsppol 
   do ikpt = 1, nkpt
!    these need to be in reduced coordinates
     write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
     do iband = 1, nband
!      output in eV
       write (iout, '(E20.10)') ha2ryd*eigen(iband, ikpt, isppol)
     end do
   end do
 end do

 close (iout)

!this file is for tau_k
 do itemp = 1, ntemper
   Temp=tempermin+temperinc*dble(itemp)

   write(appendix,"(i0)") itemp
   filename= trim(fname_radix)//"_BLZTRP.tau_k_"//trim(appendix)
   open (iout, file=filename, form='formatted')

   write (iout, '(a,f12.6)') "BoltzTraP tau_k file generated by ANADDB for T= ", Temp
   write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie(itemp), &
   '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'
   do isppol = 1, nsppol 
     do ikpt = 1, nkpt
!      these need to be in reduced coordinates
       write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
       do iband = 1, nband
!        output in sec
         write (iout, '(E20.10)') tau_k(itemp,isppol,ikpt,iband)
       end do
     end do
   end do

   close (iout)

 end do

end subroutine prtbltztrp_tau_out
!!***


