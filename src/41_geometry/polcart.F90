!{\src2tex{textfont=tt}}
!!****f* ABINIT/polcart
!! NAME
!! polcart
!!
!! FUNCTION
!! Transform polarization from reduced to cartesian coordinates,
!! divide by ucvol and write the result to an output file
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT  group (MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! pel(3)   = reduced coordinates of the electronic polarization
!! pion(3)  = reduced coordinates of the ionic polarization
!! polunit  = units used for the output of the polarization
!!             1 : use atomic units
!!             2 : use MKS units
!!             3 : use both atomic and MKS units
!! rprimd(3,3) = dimensional primitive translations (bohr)
!! ucvol    = volume of the primitive unit cell
!! unit_out = unit for output of the results
!! usepaw = 1 for PAW computation, zero else
!!
!! OUTPUT
!! pel_cart(3) = cartesian coords of the electronic polarization
!!               in atomic units
!! pelev(3)= expectation value polarization term (PAW only) in cartesian coordinates
!! pion_cart(3)= cartesian coords of the ionic polarization
!!               in atomic units
!! ptot_cart(3)= cartesian coords of the total polarization
!!               in atomic units
!!
!! NOTES
!! The sum of the electronic and ionic Berry phase is folded into
!! [-1,1] before it is transformed to cartesian coordinates.
!! This means that in some cases, ptot_cart /= pel_cart + pion_cart
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!! - pel and pion do not take into account the factor 1/ucvol.
!!   At the opposite, this factor is taken into account in
!!   pel_cart and pion_cart
!! - unit_out = 0 is allowed, in this case, there will be no
!!   output of the results
!!
!! PARENTS
!!      berryphase_new,relaxpol
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine polcart(red_ptot,pel,pel_cart,pelev,pion,pion_cart,polunit,&
&  ptot_cart,rprimd,ucvol,unit_out,usepaw)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polcart'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: polunit,unit_out,usepaw
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: red_ptot(3) !!REC
 real(dp),intent(in) :: pel(3),pelev(3),pion(3),rprimd(3,3)
 real(dp),intent(out) :: pel_cart(3),pion_cart(3),ptot_cart(3)

!Local variables -------------------------
!scalars
 integer :: idir
 character(len=500) :: message
!arrays
 real(dp) :: pel_mks(3),pelev_mks(3),pion_mks(3),ptot(3),ptot_mks(3)

! ***********************************************************************
!!REC Note ptot has already been folded and kept onto same branch
!unless ptot=0d0, in which case ptot has not been computed yet
 if( sum(abs(red_ptot(:))) < tol8 )then
   ptot(:) = pel(:) + pion(:)
!  Fold ptot into [-1, 1]
   do idir = 1, 3
     ptot(idir) = ptot(idir) - 2_dp*nint(ptot(idir)/2_dp)
   end do
 else !!REC
   ptot=red_ptot !!REC
 end if !!REC

!Transform pel, pion and ptot to cartesian coordinates
 pel_cart(:) = zero ; pion_cart(:) = zero ; ptot_cart(:) = zero
 do idir = 1, 3
   pel_cart(idir) =  rprimd(idir,1)*pel(1) + rprimd(idir,2)*pel(2) + &
&   rprimd(idir,3)*pel(3)
   pion_cart(idir) = rprimd(idir,1)*pion(1) + rprimd(idir,2)*pion(2) + &
&   rprimd(idir,3)*pion(3)
   ptot_cart(idir) = rprimd(idir,1)*ptot(1) + rprimd(idir,2)*ptot(2) + &
&   rprimd(idir,3)*ptot(3)
 end do

!Divide by the unit cell volume
 pel_cart(:)  = pel_cart(:)/ucvol
 pion_cart(:) = pion_cart(:)/ucvol
 ptot_cart(:)  = ptot_cart(:)/ucvol
!pelev is either zero (in NCPP case), or possibly non-zero (in PAW case)
!however, in the PAW case, it is already implicitly included in the computation
!of pel and should not be added in here. Below in the PAW case we report its
!value to the user, which is helpful for comparison to the USPP theory.
!13 June 2012 J Zwanziger
!note that pelev was AlREADY in cartesian frame.
!ptot_cart(:)  = (ptot_cart(:)+pelev(:))/ucvol

!Write the results to unit_out (if /= 0)
!Use the coordinates specified by the value polunit

!Atomic units

 if (((polunit == 1).or.(polunit == 3)).and.(unit_out /= 0)) then

!  write(message,'(7(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
!  &   ' Polarization in cartesian coordinates (a.u.):',ch10,&
!  &   ' (the sum of the electronic and ionic Berry phase',&
!  &   ' has been fold into [-1, 1])',ch10,&
!  &   '     Electronic berry phase:       ', (pel_cart(idir), idir = 1, 3), ch10,&
!  &   '     Expectation value (PAW only): ', (pelev(idir)/ucvol, idir = 1, 3), ch10,&
!  &   '     Ionic:                        ', (pion_cart(idir), idir = 1, 3), ch10, &
!  &   '     Total:                        ', (ptot_cart(idir), idir = 1, 3)
!  call wrtout(unit_out,message,'COLL')
   write(message,'(7(a),3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&   ' Polarization in cartesian coordinates (a.u.):',ch10,&
&   ' (the sum of the electronic and ionic Berry phase',&
&   ' has been folded into [-1, 1])',ch10,&
&   '     Electronic berry phase:       ', (pel_cart(idir), idir = 1, 3)
   call wrtout(unit_out,message,'COLL')
   if(usepaw==1) then
     write(message,'(a,3(e16.9,2x))')&
&     '     ...includes PAW on-site term: ', (pelev(idir)/ucvol, idir = 1, 3)
     call wrtout(unit_out,message,'COLL')
   end if
   write(message,'(a,3(e16.9,2x),a,a,3(e16.9,2x))')&
&   '     Ionic:                        ', (pion_cart(idir), idir = 1, 3), ch10, &
&   '     Total:                        ', (ptot_cart(idir), idir = 1, 3)
   call wrtout(unit_out,message,'COLL')

 end if

!MKS units

 if (((polunit == 2).or.(polunit == 3)).and.(unit_out /= 0)) then

   pel_mks(:)  = pel_cart(:)*(e_Cb)/(Bohr_Ang*1d-10)**2
   pion_mks(:) = pion_cart(:)*(e_Cb)/(Bohr_Ang*1d-10)**2
   pelev_mks(:) = pelev(:)/ucvol*(e_Cb)/(Bohr_Ang*1d-10)**2
   ptot_mks(:) = (ptot_cart(:))*(e_Cb)/(Bohr_Ang*1d-10)**2

   write(message,'(7(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&   ' Polarization in cartesian coordinates (C/m^2):',ch10,&
&   ' (the sum of the electronic and ionic Berry phase',&
&   ' has been folded into [-1, 1])',ch10,&
&   '     Electronic berry phase:       ', (pel_mks(idir), idir = 1, 3)
   call wrtout(unit_out,message,'COLL')
   if(usepaw==1) then
     write(message,'(a,3(e16.9,2x))')&
&     '     ...includes PAW on-site term: ', (pelev_mks(idir), idir = 1, 3)
     call wrtout(unit_out,message,'COLL')
   end if
   write(message,'(a,3(e16.9,2x),a,a,3(e16.9,2x))')&
&   '     Ionic:                        ', (pion_mks(idir), idir = 1, 3), ch10, &
&   '     Total:                        ', (ptot_mks(idir), idir = 1, 3)
   call wrtout(unit_out,message,'COLL')

 end if

 end subroutine polcart
!!***
