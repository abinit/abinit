!!****m* ABINIT/m_raman
!! NAME
!!  m_raman
!!
!! FUNCTION
!!  Raman susceptibilities of zone-center phonons and electroo tensor.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (MVeithen)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_raman

 use defs_basis
 use m_errors
 use m_abicore

 use m_symtk,           only : matr3inv

 implicit none

 private
!!***

 public :: ramansus        ! Raman susceptibilities of zone-center phonons.
 public :: electrooptic    ! Electrooptic tensor and the raman tensors of zone-center phonons.
!!***

contains
!!***

!!****f* m_raman/ramansus
!!
!! NAME
!! ramansus
!!
!! FUNCTION
!! Compute the raman susceptibilities of zone-center phonons
!!
!! INPUTS
!!  d2cart = second order derivatives of the energy wrt all perturbations
!!  dchide(3,3,3) = non-linear optical coefficients from dtchi
!!  dchidt(natom,3,3,3) = first-order change of the electronic dielectric
!!     tensor induced by an individual atomic displacement
!!  displ = phonon mode atomic displacements
!!  mpert = maximum number of perturbations
!!  natom = number of atoms
!!  phfrq = phonon frequencies
!!  qphnrm=(described below)
!!  ucvol = unit cell volume
!!
!! OUTPUT
!!  qphon(3)= to be divided by qphnrm, give the phonon wavevector;
!!     if qphnrm==0.0_dp, then the wavevector is zero (Gamma point)
!!     and qphon gives the direction of
!!     the induced electric field; in the latter case, if qphon is
!!     zero, no non-analytical contribution is included.
!!  rsus
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ramansus(d2cart,dchide,dchidt,displ,mpert,natom,phfrq,qphon,qphnrm,rsus,ucvol)

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: qphnrm,ucvol
!arrays
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert),dchide(3,3,3)
 real(dp),intent(in) :: dchidt(natom,3,3,3),displ(2,3*natom,3*natom)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(inout) :: qphon(3)
 real(dp),intent(out) :: rsus(3*natom,3,3)

!Local variables-------------------------------
!scalars
 integer :: analyt,i1,i1dir,i1pert,i2dir,iatom,idir,imode
 real(dp) :: epsq,fac,g0,g1,g2,qphon2
 logical :: t_degenerate,iwrite
 character(len=500) :: message
!arrays
 real(dp) :: dijk_q(3,3)
 real(dp),allocatable :: zeff(:,:)
 character(len=1),allocatable :: metacharacter(:)

! *********************************************************************

 iwrite = ab_out > 0

 ABI_ALLOCATE(zeff,(3,natom))

 rsus(:,:,:) = zero
 epsq        = zero
 zeff(:,:)   = zero
 dijk_q(:,:) = zero

!Determine the analyticity of the matrix.
 analyt=1
 if(abs(qphnrm)<tol8)analyt=0
 if(abs(qphon(1))<tol8.and.abs(qphon(2))<tol8.and.abs(qphon(3))<tol8) analyt=1

!In the case the non-analyticity is required :
 if(analyt == 0) then

!  Normalize the limiting direction
   qphon2=qphon(1)**2+qphon(2)**2+qphon(3)**2
   qphon(1)=qphon(1)/sqrt(qphon2)
   qphon(2)=qphon(2)/sqrt(qphon2)
   qphon(3)=qphon(3)/sqrt(qphon2)

!  Get the dielectric constant for the limiting direction
   epsq= 0._dp
   do i1dir=1,3
     do i2dir=1,3
       epsq=epsq+qphon(i1dir)*qphon(i2dir)*d2cart(1,i1dir,natom+2,i2dir,natom+2)
     end do
   end do

!  Check if epsq > 0
   if (epsq < tol8) then
     write(message,'(a,es14.6)')'  The value of epsq must be > 0 while it is found to be',epsq
     MSG_BUG(message)
   end if

!  Get the effective charges for the limiting direction
   do i1dir=1,3
     do i1pert=1,natom
       zeff(i1dir,i1pert)=zero
       do i2dir=1,3
         zeff(i1dir,i1pert)=zeff(i1dir,i1pert)+qphon(i2dir)*&
&         d2cart(1,i1dir,i1pert,i2dir,natom+2)
       end do
     end do
   end do

!  Get the NLO tensor for the limiting direction !$\sum_{k} d_{ijk} \cdot q_k$

   dijk_q(:,:) = zero
   do i1dir = 1, 3
     do i2dir = 1, 3
       do idir = 1, 3
         dijk_q(i1dir,i2dir) = dijk_q(i1dir,i2dir) + dchide(i1dir,i2dir,idir)*qphon(idir)
       end do
     end do
   end do

   fac = 16._dp*pi/(ucvol*epsq)
   do imode = 1, 3*natom
     do iatom = 1, natom
       do idir = 1, 3
         i1=idir + (iatom - 1)*3
         rsus(imode,:,:) = rsus(imode,:,:) + &
&         (dchidt(iatom,idir,:,:) - fac*zeff(idir,iatom)*dijk_q(:,:))* &
&         displ(1,i1,imode)
       end do  ! disp
     end do  ! iatom
   end do  ! imode
   rsus(:,:,:) = rsus(:,:,:)*sqrt(ucvol)

 else
   do imode = 1, 3*natom
     do iatom = 1, natom
       do idir = 1, 3
         i1=idir + (iatom - 1)*3
         rsus(imode,:,:) = rsus(imode,:,:) + dchidt(iatom,idir,:,:)*displ(1,i1,imode)
       end do  ! disp
     end do  ! iatom
   end do  ! imode
   rsus(:,:,:) = rsus(:,:,:)*sqrt(ucvol)
 end if      ! analyt == 0

 if (analyt == 0) then
   if (iwrite) then
     write(ab_out,*) ch10
     write(ab_out, '(a,/,a,3f9.5)' )&
&     ' Raman susceptibility of zone-center phonons, with non-analyticity in the',&
&     '  direction (cartesian coordinates)',qphon(1:3)+tol10
     write(ab_out,'(a)')&
&     ' -----------------------------------------------------------------------'
     write(ab_out,*) ch10
   end if

 else
   if (iwrite) then
     write(ab_out,*) ch10
     write(ab_out,*)' Raman susceptibilities of transverse zone-center phonon modes'
     write(ab_out,*)' -------------------------------------------------------------'
     write(ab_out,*) ch10
   end if
 end if

!Examine the degeneracy of each mode. The portability of the echo of the Raman susceptibility
!for each degenerate mode is very hard to guarantee. On the contrary,
!the scalar reductions of these quantities are OK.
 ABI_ALLOCATE(metacharacter,(3*natom))
 do imode=1,3*natom
!  The degenerate modes are not portable
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

 do imode = 1, 3*natom
   if (iwrite) then
     write(ab_out,'(a4,i3,2x,a2,f7.2,a6)')' Mode',imode,' (',phfrq(imode)*Ha_cmm1,' cm-1)'
     do idir = 1,3
       write(ab_out,'(a,4x,3(f16.9,2x))')metacharacter(imode),rsus(imode,idir,:)
     end do
   end if
!  See R. Caracas and X. Gonze, Thermodynamic Properties of Solids : experiment and modeling, Wiley-VCH,
!  Ed. S. Chaplot and R. Mittal and N. Choudhury , chap. 8, pp 291-312.
   g0=(rsus(imode,1,1)+rsus(imode,2,2)+rsus(imode,3,3))**2*third
   g1=((rsus(imode,1,2)-rsus(imode,2,1))**2+&
&   (rsus(imode,1,3)-rsus(imode,3,1))**2+&
&   (rsus(imode,2,3)-rsus(imode,3,2))**2)*half
   g2=((rsus(imode,1,2)+rsus(imode,2,1))**2+&
&   (rsus(imode,1,3)+rsus(imode,3,1))**2+&
&   (rsus(imode,2,3)+rsus(imode,3,2))**2)*half +&
&   ((rsus(imode,1,1)-rsus(imode,2,2))**2+&
&   (rsus(imode,2,2)-rsus(imode,3,3))**2+&
&   (rsus(imode,3,3)-rsus(imode,1,1))**2)*third
   if (iwrite) then
     write(ab_out,'(3(a,f16.9))')' Spherical averages : G0=',g0,'    G1=',g1,'    G2=',g2
     write(ab_out,*)
   end if
 end do

 ABI_DEALLOCATE(metacharacter)
 ABI_DEALLOCATE(zeff)

end subroutine ramansus
!!***

!!****f* m_raman/electrooptic
!! NAME
!! electrooptic
!!
!! FUNCTION
!! Compute the electrooptic tensor and the raman tensors of zone-center phonons
!!
!! INPUTS
!! dchide(3,3,3) = non-linear optical coefficients
!! dieflag= dielectric tensor flag. 0=> no dielectric tensor,
!!  1=> frequency-dependent dielectric tensor,
!!  2=> only the electronic dielectric tensor.
!! epsinf=electronic dielectric tensor
!! fact_oscstr(2,3,3*natom)=factors of the oscillator strengths for the different eigenmodes,
!!  for different direction of the electric field
!! natom=number of atoms in unit cell
!! phfrq(3*natom)=phonon frequencies (square root of the dynamical
!!  matrix eigenvalues, except if these are negative, and in this
!!  case, give minus the square root of the absolute value
!!  of the matrix eigenvalues). Hartree units.
!! prtmbm= if equal to 1 write out the mode by mode decomposition of the EO tensor
!! rsus = Raman susceptibilities
!! ucvol=unit cell volume
!!
!! OUTPUT
!!  (to be completed ?)
!!
!! NOTES
!! 1. The phonon frequencies phfrq should correspond to the
!! wavevector at Gamma, without any non-analyticities.
!! 2. Should clean for no imaginary part ...
!! This routine should be used only by one processor.
!! 3. frdiel(3,3,nfreq)= frequency-dependent dielectric tensor
!! mode effective charges for the different eigenmodes,
!! for different direction of the electric field
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

subroutine electrooptic(dchide,dieflag,epsinf,fact_oscstr,natom,phfrq,prtmbm,rsus,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: dieflag,natom,prtmbm
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: dchide(3,3,3),epsinf(3,3),fact_oscstr(2,3,3*natom)
 real(dp),intent(in) :: phfrq(3*natom),rsus(3*natom,3,3)

!Local variables -------------------------
!scalars
 integer :: flag,i1,i2,ii,imode,jj,kk
 real(dp) :: dtm,fac
 logical :: iwrite
 character(len=500) :: message
!arrays
 integer :: voigtindex(6,2)
 real(dp) :: eta(3,3),rvoigt(6,3),work(3,3,3)
 real(dp),allocatable :: rijk(:,:,:,:),rijk_tot(:,:,:)

! *********************************************************************

!rijk(1:3*natom,:,:,:) = mode by mode decomposition of the electrooptic tensor
!rijk(3*natom+1,:,:,:) = electronic contribution
 iwrite = ab_out > 0

 voigtindex(1,1) = 1 ; voigtindex(1,2) = 1
 voigtindex(2,1) = 2 ; voigtindex(2,2) = 2
 voigtindex(3,1) = 3 ; voigtindex(3,2) = 3
 voigtindex(4,1) = 2 ; voigtindex(4,2) = 3
 voigtindex(5,1) = 1 ; voigtindex(5,2) = 3
 voigtindex(6,1) = 1 ; voigtindex(6,2) = 2

 ABI_ALLOCATE(rijk,(3*natom+1,3,3,3))
 ABI_ALLOCATE(rijk_tot,(3,3,3))
 rijk(:,:,:,:) = 0._dp
 rijk_tot(:,:,:) = 0._dp


!In case there is no mode with truly negative frequency
!and the electronic dielectric tensor is available
!compute the electro-optic tensor

 flag = 1

 if (abs(phfrq(1)) > abs(phfrq(4))) then
   flag = 0
   write(message,'(6a)')&
&   'The lowest mode appears to be a "true" negative mode,',ch10,&
&   'and not an acoustic mode. This precludes the computation',ch10,&
&   'of the EO tensor.',ch10
   MSG_WARNING(message)
 end if

 dtm = epsinf(1,1)*epsinf(2,2)*epsinf(3,3) + &
& epsinf(1,2)*epsinf(2,3)*epsinf(3,1) + &
& epsinf(1,3)*epsinf(2,1)*epsinf(3,2) - &
& epsinf(3,1)*epsinf(2,2)*epsinf(1,3) - &
& epsinf(3,2)*epsinf(2,3)*epsinf(1,1) - &
& epsinf(3,3)*epsinf(2,1)*epsinf(1,2)

 if (abs(dtm) < tol6) then
   flag = 0
   write(message,'(a,a,a,a,a,a,a,a)')&
&   'The determinant of the electronic dielectric tensor is zero.',ch10,&
&   'This preludes the computation fo the EO tensor since',ch10,&
&   'this quantity requires the inverse of epsilon.',ch10,&
&   'Action : check you database and the value of dieflag in the input file.',ch10
   MSG_WARNING(message)
 end if

!dieflag is required to be one since the EO tensor
!requires the oscillator strengths

 if ((flag == 1).and.(dieflag==1)) then

!  Factor to convert atomic units to MKS units

   fac = -16._dp*pi*pi*eps0*(Bohr_Ang**2)*1.0d-8/(e_Cb*sqrt(ucvol))

!  Compute inverse of dielectric tensor
!  needed to convert the nonlinear optical susceptibility tensor
!  to the electrooptic tensor

   call matr3inv(epsinf,eta)

   if (iwrite) then
     write(ab_out,*)ch10
     write(ab_out,*)'Output of the EO tensor (pm/V) in Voigt notations'
     write(ab_out,*)'================================================='
     write(ab_out,*)
     if (prtmbm == 1) then
       write(ab_out,*)'Mode by mode decomposition'
       write(ab_out,*)
     end if
   end if

!  Compute the ionic contribution to the EO tensor

   do imode = 4, 3*natom

     if (prtmbm == 1 .and. iwrite) then
       write(ab_out,*)
       write(ab_out,'(a4,i3,2x,a2,f7.2,a6)')'Mode',imode,' (',phfrq(imode)*Ha_cmm1,' cm-1)'
     end if

     do ii = 1, 3
       do jj = 1, 3
         do kk = 1, 3
           rijk(imode,ii,jj,kk) = rsus(imode,ii,jj)*fact_oscstr(1,kk,imode)/(phfrq(imode)**2)
         end do
       end do
     end do

     work(:,:,:) = 0._dp
     do ii = 1,3
       do jj = 1, 3
         do kk = 1, 3

           do i1 = 1, 3
             do i2 = 1, 3
               work(ii,jj,kk) = work(ii,jj,kk) + eta(ii,i1)*rijk(imode,i1,i2,kk)*eta(i2,jj)
             end do  ! i2
           end do   ! i1

           rijk(imode,ii,jj,kk) = fac*work(ii,jj,kk)
           rijk_tot(ii,jj,kk) = rijk_tot(ii,jj,kk) + rijk(imode,ii,jj,kk)
         end do

       end do
     end do

     if (prtmbm == 1) then
       rvoigt(:,:) = 0._dp
       do i1 = 1, 6
         ii = voigtindex(i1,1)
         jj = voigtindex(i1,2)
         do kk = 1, 3
           rvoigt(i1,kk) = (rijk(imode,ii,jj,kk) + rijk(imode,jj,ii,kk))/2._dp
         end do
         if (iwrite) write(ab_out,'(5x,3(2x,f16.9))')rvoigt(i1,:)
       end do
     end if

   end do     ! imode

!  Compute the electronic contribution to the EO tensor

   if (prtmbm == 1 .and. iwrite) then
     write(ab_out,*)
     write(ab_out,*)'Electronic contribution to the EO tensor'
   end if

   fac = 16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

   do ii = 1,3
     do jj = 1, 3
       do kk = 1, 3

         do i1 = 1, 3
           do i2 = 1, 3
             rijk(3*natom+1,ii,jj,kk) = rijk(3*natom+1,ii,jj,kk) + &
&             eta(ii,i1)*dchide(i1,i2,kk)*eta(i2,jj)
           end do  ! i2
         end do   ! i1

         rijk(3*natom+1,ii,jj,kk) = -4._dp*rijk(3*natom+1,ii,jj,kk)*fac
         rijk_tot(ii,jj,kk) = rijk_tot(ii,jj,kk) + rijk(3*natom+1,ii,jj,kk)

       end do

     end do
   end do

   if (prtmbm == 1) then
     rvoigt(:,:) = 0._dp
     do i1 = 1, 6
       ii = voigtindex(i1,1)
       jj = voigtindex(i1,2)
       do kk = 1, 3
         rvoigt(i1,kk) = (rijk(3*natom+1,ii,jj,kk) + rijk(3*natom+1,jj,ii,kk))/2._dp
       end do
       if (iwrite) write(ab_out,'(5x,3(2x,f16.9))')rvoigt(i1,:)
     end do
     if (iwrite) write(ab_out,*)ch10
   end if

   if (iwrite) write(ab_out,*)'Total EO tensor (pm/V) in Voigt notations'
   rvoigt(:,:) = 0._dp
   do i1 = 1, 6
     ii = voigtindex(i1,1)
     jj = voigtindex(i1,2)
     do kk = 1, 3
       rvoigt(i1,kk) = (rijk_tot(ii,jj,kk) + rijk_tot(jj,ii,kk))/2._dp
     end do
     if (iwrite) write(ab_out,'(5x,3(2x,f16.9))')rvoigt(i1,:)
   end do

 end if  ! flag

 ABI_DEALLOCATE(rijk)
 ABI_DEALLOCATE(rijk_tot)

end subroutine electrooptic
!!***

end module m_raman
!!***
