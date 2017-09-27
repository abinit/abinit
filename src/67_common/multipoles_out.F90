!{\src2tex{textfont=tt}}
!!****f* ABINIT/multipoles_out
!! NAME
!! multipoles_out
!!
!! FUNCTION
!!  Output multipole moments of input array on FFT grid, calculated with multipoles_fftr
!!  Namely, the electrical dipole, quadrupole, etc... of the electron density
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (MJV,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rhor(nfft,nspden)=electronic density
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nfft=number of FFT points stored by one proc
!!  ngfft =number of subdivisions along each lattice vector
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  rprimd(3,3)=dimensionful lattice vectors
!!  typat(ntypat)=types of atoms
!!  ucvol=unit cell volume
!!  unit_out=file unit to print out
!!  ziontypat(ntypat)=ionic charge of each type of atom
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      multipoles_fftr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine multipoles_out(rhor,mpi_enreg,natom,nfft,ngfft,nspden,&
&                         ntypat,rprimd,typat,ucvol,unit_out,xred,ziontypat)

 use m_profiling_abi
 use defs_basis
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multipoles_out'
 use interfaces_14_hidewrite
 use interfaces_53_spacepar
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,unit_out
 real(dp), intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer, intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),xred(3,natom),ziontypat(ntypat)
!Local variables ------------------------------
!scalars
 integer :: iatom,nspden_updn
 real(dp) :: ziontotal
 character(len=500) :: message
!arrays
 real(dp) :: center_of_charge(3),dipole_el(3,2),dipole_ions_cart(3)
 real(dp) :: dipole_ions_red(3),dipole_tot(3),tmp(3)

! *************************************************************************

!Separate spins only for nspden=2
 nspden_updn=merge(1,2,nspden/=2)

!Title

!Get nuclear part of dipole
 dipole_ions_red(:) = zero ; ziontotal = zero
 do iatom = 1, natom
   dipole_ions_red(:) = dipole_ions_red(:) + xred(:,iatom)*ziontypat(typat(iatom))
   ziontotal = ziontotal + ziontypat(typat(iatom))
 end do
 dipole_ions_cart(:) = matmul(rprimd,dipole_ions_red(:))

!Find coordinates of center of charge on FFT grid
 center_of_charge(1:3) = dipole_ions_red(1:3)/ziontotal

!Get electronic part of dipole with respect to center of charge of ions (in cart. coord.)
 dipole_el = zero
 call multipoles_fftr(rhor(:,1:nspden_updn),dipole_el(:,1:nspden_updn),nfft,ngfft,nspden_updn,&
& rprimd,center_of_charge,distribfft=mpi_enreg%distribfft,mpi_comm_grid=mpi_enreg%comm_fft)
 dipole_el(1:3,1:nspden_updn)=dipole_el(1:3,1:nspden_updn)*ucvol

!Take into account storage of rhor (up+dn,up)
 if (nspden==2) then
   tmp(1:3)=dipole_el(1:3,1)
   dipole_el(1:3,1)=dipole_el(1:3,2)
   dipole_el(1:3,2)=tmp(1:3)-dipole_el(1:3,2)
 end if

!Compute total dipole
!NOTE: wrt center of charge, dipole_ions is 0
 dipole_tot(1) = - sum(dipole_el(1,1:nspden_updn))
 dipole_tot(2) = - sum(dipole_el(2,1:nspden_updn))
 dipole_tot(3) = - sum(dipole_el(3,1:nspden_updn))

!Output
 write (message, '(2a)') ch10,' ----- Electric nuclear dipole wrt the center of ionic charge ----- '
 call wrtout(unit_out, message, 'COLL')
 write (message, '(a,3(1x,ES12.5))') &
& ' Center of charge for ionic distribution (red. coord.): ',center_of_charge(1:3)
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Ionic dipole (cart. coord.)     = ',dipole_ions_cart, ' (a.u.)',ch10, &
& '                                 = ',dipole_ions_cart/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 if (nspden/=2) then
   !This is compatible with nspden=4
   write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' (D)'
 else
   write (message, '(3a,3(1x,E16.6),a,3(2a,3(1x,E16.6),a))') ' -----',ch10,&
&   ' Electronic dipole (cart. coord.)= ',dipole_el(:,1),' up (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,2),' dn (a.u.)',ch10,&
&   '                                 = ',dipole_el(:,1)/dipole_moment_debye,' up (D)',ch10,&
&   '                                 = ',dipole_el(:,2)/dipole_moment_debye,' dn (D)'
 end if
 call wrtout(unit_out, message, 'COLL')
 write (message, '(3a,3(1x,E16.6),3a,3(1x,E16.6),a)') ' -----',ch10,&
& ' Total dipole (cart. coord.)     = ',dipole_tot,' (a.u.)',ch10,&
& '                                 = ',dipole_tot/dipole_moment_debye,' (D)'
 call wrtout(unit_out, message, 'COLL')
 call wrtout(unit_out, ' ', 'COLL')

end subroutine multipoles_out
!!***
