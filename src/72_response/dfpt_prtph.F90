!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_prtph
!!
!! NAME
!! dfpt_prtph
!!
!! FUNCTION
!! Print the phonon frequencies, on unit 6 as well as the printing
!! unit (except if the associated number -iout- is negative),
!! and for the latter, in Hartree, meV, Thz, Kelvin or cm-1.
!! If eivec==1,2, also print the eigenmodes : displacements in cartesian coordinates.
!! If eivec==4, generate output files for band2eps (drawing tool for the phonon band structure
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  displ(2,3*natom,3*natom)= contains the displacements of atoms in cartesian coordinates.
!!  The first index means either the real or the imaginary part,
!!  The second index runs on the direction and the atoms displaced
!!  The third index runs on the modes.
!!  eivec=(if eivec==0, the eigendisplacements are not printed,
!!    if eivec==1,2, the eigendisplacements are printed,
!!    if eivec==4, files for band2eps
!!  enunit=units for output of the phonon frequencies :
!!    0=> Hartree and cm-1, 1=> eV and Thz, other=> Ha,Thz,eV,cm-1 and K
!!  iout= unit for long print (if negative, the routine only print on unit 6, and in Hartree only).
!!  natom= number of atom
!!  phfreq(3*natom)= phonon frequencies in Hartree
!!  qphnrm=phonon wavevector normalisation factor
!!  qphon(3)=phonon wavevector
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!! called by one processor only
!!
!! PARENTS
!!      anaddb,m_ifc,m_phonons,respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_prtph(displ,eivec,enunit,iout,natom,phfrq,qphnrm,qphon)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_prtph'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: eivec,enunit,iout,natom
 real(dp),intent(in) :: qphnrm
!arrays
 real(dp),intent(in) :: displ(2,3*natom,3*natom),phfrq(3*natom),qphon(3)

!Local variables -------------------------
!scalars
 integer :: i,idir,ii,imode,jj
 real(dp) :: tolerance
 logical :: t_degenerate
 character(len=500) :: message
!arrays
 real(dp) :: vecti(3),vectr(3)
 character(len=1),allocatable :: metacharacter(:)

! *********************************************************************

!Check the value of eivec
 if (all(eivec /= [0,1,2,4])) then
   write(message, '(a,i0,a,a)' )&
&   'In the calling subroutine, eivec is',eivec,ch10,&
&   'but allowed values are between 0 and 4.'
   MSG_BUG(message)
 end if

!write the phonon frequencies on unit std_out
 write(message,'(4a)' )' ',ch10,&
& ' phonon wavelength (reduced coordinates) , ','norm, and energies in hartree'
 call wrtout(std_out,message,'COLL')

!The next format should be rewritten
 write(message,'(a,4f5.2)' )' ',(qphon(i),i=1,3),qphnrm
 call wrtout(std_out,message,'COLL')
 do jj=1,3*natom,5
   if (3*natom-jj<5) then
     write(message,'(5es17.9)') (phfrq(ii),ii=jj,3*natom)
   else
     write(message,'(5es17.9)') (phfrq(ii),ii=jj,jj+4)
   end if
   call wrtout(std_out,message,'COLL')
 end do
 write(message,'(a,es17.9)')' Zero Point Motion energy (sum of freqs/2)=',sum(phfrq(1:3*natom))/2
 call wrtout(std_out,message,'COLL')

!Put the wavevector in nice format
 if(iout>=0)then
   write(iout, '(a)' )' '
   if(qphnrm/=0.0_dp)then
     write(message, '(a,3f9.5)' )&
&     '  Phonon wavevector (reduced coordinates) :',(qphon(i)/qphnrm+tol10,i=1,3)
   else
     write(message, '(a,/,a,3f9.5)' )&
&     '  Phonon at Gamma, with non-analyticity in the',&
&     '  direction (cartesian coordinates)',qphon(1:3)+tol10
   end if
   call wrtout(iout,message,'COLL')

!  Write it, in different units.
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Hartree :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(iout, '(1x,5es14.6)') (phfrq(ii),ii=jj,3*natom)
       else
         write(message, '(1x,5es14.6)') (phfrq(ii),ii=jj,jj+4)
       end if
       call wrtout(iout,message,'COLL')
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon energies in meV     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,3*natom)
       else
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_eV*1.0d3,ii=jj,jj+4)
       end if
       call wrtout(iout,message,'COLL')
     end do
   end if
   if(enunit/=1)then
     write(iout, '(a)' )' Phonon frequencies in cm-1    :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,3*natom)
       else
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_cmm1,ii=jj,jj+4)
       end if
       call wrtout(iout,message,'COLL')
     end do
   end if
   if(enunit/=0)then
     write(iout, '(a)' )' Phonon frequencies in Thz     :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,3*natom)
       else
         write(message, '("-",5es14.6)') (phfrq(ii)*Ha_THz,ii=jj,jj+4)
       end if
       call wrtout(iout,message,'COLL')
     end do
   end if
   if(enunit/=0.and.enunit/=1)then
     write(iout, '(a)' )' Phonon energies in Kelvin  :'
     do jj=1,3*natom,5
       if (3*natom-jj<5) then
         write(message, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,3*natom)
       else
         write(message, '("-",5es14.6)') (phfrq(ii)/kb_HaK,ii=jj,jj+4)
       end if
       call wrtout(iout,message,'COLL')
     end do
   end if
 end if

!Take care of the eigendisplacements
 if(eivec==1 .or. eivec==2)then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' Eigendisplacements ',ch10,&
&   ' (will be given, for each mode : in cartesian coordinates',ch10,&
&   '   for each atom the real part of the displacement vector,',ch10,&
&   '   then the imaginary part of the displacement vector - absolute values smaller than 1.0d-7 are set to zero)'
   call wrtout(std_out,message,'COLL')
   if(iout>=0) then
     call wrtout(iout,message,'COLL')
   end if

!  Examine the degeneracy of each mode. The portability of the echo of the eigendisplacements
!  is very hard to obtain, and has not been attempted.
   ABI_ALLOCATE(metacharacter,(3*natom))
   do imode=1,3*natom
!    The degenerate modes are not portable
     t_degenerate=.false.
     if(imode>1)then
       if(phfrq(imode)-phfrq(imode-1)<tol6)t_degenerate=.true.
     end if
     if(imode<3*natom)then
       if(phfrq(imode+1)-phfrq(imode)<tol6)t_degenerate=.true.
     end if
     metacharacter(imode)=';'
     if(t_degenerate)metacharacter(imode)='-'
   end do

   do imode=1,3*natom
     write(message,'(a,i4,a,es16.6)' )'  Mode number ',imode,'   Energy',phfrq(imode)
     call wrtout(std_out,message,'COLL')
     if(iout>=0)then
       write(message, '(a,i4,a,es16.6)' )'  Mode number ',imode,'   Energy',phfrq(imode)
       call wrtout(iout,message,'COLL')
     end if
     tolerance=1.0d-7
     if(abs(phfrq(imode))<1.0d-5)tolerance=2.0d-7
     if(phfrq(imode)<1.0d-5)then
       write(message,'(3a)' )' Attention : low frequency mode.',ch10,&
&                           '   (Could be unstable or acoustic mode)'
       call wrtout(std_out,message,'COLL')
       if(iout>=0)then
         write(iout, '(3a)' )' Attention : low frequency mode.',ch10,&
&                            '   (Could be unstable or acoustic mode)'
       end if
     end if
     do ii=1,natom
       do idir=1,3
         vectr(idir)=displ(1,idir+(ii-1)*3,imode)
         if(abs(vectr(idir))<tolerance)vectr(idir)=0.0_dp
         vecti(idir)=displ(2,idir+(ii-1)*3,imode)
         if(abs(vecti(idir))<tolerance)vecti(idir)=0.0_dp
       end do
       write(message,'(i4,3es16.8,a,4x,3es16.8)' ) ii,vectr(:),ch10,vecti(:)
       call wrtout(std_out,message,'COLL')
       if(iout>=0)then
         write(message,'(a,i3,3es16.8,2a,3x,3es16.8)') metacharacter(imode),ii,vectr(:),ch10,&
&                                                      metacharacter(imode),   vecti(:)
         call wrtout(iout,message,'COLL')
       end if
     end do
   end do

   ABI_DEALLOCATE(metacharacter)
 end if

end subroutine dfpt_prtph
!!***
