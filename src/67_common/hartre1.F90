!{\src2tex{textfont=tt}}
!!****f* ABINIT/hartre1
!! NAME
!! hartre1
!!
!! FUNCTION
!! Compute the Hartree energy from the density in reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, MF, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! The Hartree energy is given for the Coulomb interaction with a
!! real space cutoff, $\theta(r-R)/r$, as
!!  ehvalues(3)=$\displaystyle \frac{1}{2\pi}\sum_{\vec G}
!!          \frac{|n(\vec G)|^2}{G^2}\left( 1-cos(GR) \right)$,
!! and for the true Coulomb interaction, $1/r$, as
!!  ehvalues(2)=$\displaystyle \frac{1}{2\pi}\sum_{\vec G \neq 0}\frac{|n(\vec G)|^2}{G^2}$.
!! Also give
!!  ehvalues(1)=$\displaystyle \frac{1}{2\pi}\sum_{\vec G \neq 0}
!!          \frac{|n(\vec G)|^2}{G^2}\left( 1-cos(GR) \right)$.
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!     (gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2))
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rcut=cutoff radius for Coulomb interaction in bohr
!!  rhog(2,nfft)=electron density in G space
!!  ucvol=unit cell volume
!!
!! OUTPUT
!!  ehvalues(1)=Hartree energy with cutoff interaction without $G=0$ term
!!  ehvalues(2)=Hartree energy with full interaction without $G=0$ term
!!  ehvalues(2)=Hartree energy with cutoff interaction including $G=0$ term
!!  vhart(cplex*nfft)=Hartree potential in real space, either REAL or COMPLEX
!!
!! WARNINGS
!! Case cplex=2 not implemented.
!! Hartree potential is not computed.
!!
!! TODO
!! Implement case cplex=2, clean up and calculate Hartree potential.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hartre1(cplex,gmet,gsqcut,nfft,ngfft,paral_kgb,qphon,rhog,vhartr,ehvalues,rcut,ucvol)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hartre1'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,paral_kgb
 real(dp),intent(in) :: gsqcut,rcut,ucvol
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3),rhog(2,nfft)
 real(dp),intent(out) :: ehvalues(3),vhartr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,ig,ii,ii1,ing,qeq0
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg
!arrays
 integer :: id(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************

!Keep track of total time spent in hartre
 call timab(10,1,tsec)

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message, '(a,i0,a,a)' )&
&   'From the calling routine, cplex=',cplex,ch10,&
&   'but the only value allowed are 1 and 2.'
   MSG_BUG(message)
 end if

!Initialize a few quantities
 cutoff=gsqcut*tolfix
!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message,'(a,3e12.4,a,a)')&
&   'cplex=1 but qphon=',qphon,ch10,&
&   'qphon should be 0 0 0.'
   MSG_BUG(message)
 end if

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(ngfft(1),ngfft(2),ngfft(3))))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do

 ABI_ALLOCATE(work1,(2,nfft))

!Triple loop on each dimension
 ehvalues(:)=zero
 do i3=1,ngfft(3)

!  Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,ngfft(2)
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     i23=ngfft(1)*((i2-1)+ngfft(2)*(i3-1))
!    Do the test that eliminates the Gamma point outside
!    of the inner loop
     ii1=1
     if(i23==0 .and. qeq0==1)then
       ii1=2
       work1(re,1+i23)=0.0_dp
       work1(im,1+i23)=0.0_dp
     end if

!    Final inner loop on the first dimension
!    (note the lower limit)
     do i1=ii1,ngfft(1)
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ii=i1+i23
       if(gs<=cutoff)then
         den=piinv/gs
         work1(re,ii)=rhog(re,ii)*den
         work1(im,ii)=rhog(im,ii)*den
!        MF evaluate reciprocal space hartree energy
         den=(rhog(re,ii)**2+rhog(im,ii)**2)/gs
!        MF use cutoff Coulomb interaction
         ehvalues(re)=ehvalues(re)+den*(1._dp-cos(2._dp*pi*sqrt(gs)*rcut))
!        MF use full Coulomb interaction
         ehvalues(im)=ehvalues(im)+den
       else
         work1(re,ii)=0.0_dp
         work1(im,ii)=0.0_dp
       end if
!      End loop on i1
     end do

!    End loop on i2
   end do

!  End loop on i3
 end do
 ehvalues(3)=ehvalues(re)+(rhog(re,1)**2+rhog(im,1)**2)*0.5_dp*(2*pi*rcut)**2

 ABI_DEALLOCATE(gq)

!MF
!these are the correct pi factors,
!numerator:   $1/2$ [double counting] $\times  4\pi$ [Poisson eq] $= 2\pi$
!denominator: $2\pi$ [reciprocal lattice vectors] squared         $= (2\pi)^2$
!gives the same result as real space evaluation via
!$\frac{1}{2}\int n(\vec r)*V_{H}(\vec r) d^3r$
 ehvalues(:)=ehvalues(:)*ucvol*2*pi/(2*pi)**2
!MF

!DEBUG
!write(std_out,*)' hartre : before fourdp'
!write(std_out,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(std_out,*)' maxval work1=',maxval(abs(work1))
!ENDDEBUG

!Fourier Transform Vhartree.
!Vh in reciprocal space was stored in work1
!MF no need to calculate potential
 if(.false.)then
   mpi_enreg%me_fft=0
   mpi_enreg%nproc_fft=1
   call fourdp(cplex,work1,vhartr,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
 end if

 ABI_DEALLOCATE(work1)

 call timab(10,2,tsec)

end subroutine hartre1
!!***
