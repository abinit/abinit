!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkcore
!! NAME
!! dfpt_mkcore
!!
!! FUNCTION
!! Compute the derivative of the core electron density
!! with respect to one specific atom displacement
!! In case of derivative with respect to k or
!! electric field perturbation, the 1st-order core electron density
!! vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis) or cartesian coordinate
!!    pair for strain perturbation
!!  ipert=number of the atom being displaced or natom+3,4 for strain
!!    perturbation
!!  natom=number of atoms in cell.
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid dimensions.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  xccc3d1(cplex*n1*n2*n3)=3D core electron density for XC core correction,
!!    bohr^-3
!!
!! NOTES
!! Note that this routine is tightly connected to the mkcore.f routine
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_eltfrxc,dfpt_looppert,dfpt_nselt,dfpt_nstdy,dfpt_nstpaw
!!      dfptnl_loop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_mkcore(cplex,idir,ipert,natom,ntypat,n1,n1xccc,&
& n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xccc3d1,xred)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkcore'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n1,n1xccc,n2,n3,natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: qphon(3),rprimd(3,3),xccc1d(n1xccc,6,ntypat)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(out) :: xccc3d1(cplex*n1*n2*n3)

!Local variables-------------------------------
!scalars
 integer,parameter :: mshift=401
 integer :: i1,i2,i3,iatom,ifft,ishift,ishift1,ishift2,ishift3,istr
 integer :: ixp,jj,ka,kb,mrange,mu,nu
 real(dp) :: aa,bb,cc,dd,delta,delta2div6,deltam1,diff,difmag
 real(dp) :: difmag2,difmag2_fact,difmag2_part,func,phase,phi,phr,prod
 real(dp) :: range,range2,rangem1,rdiff1,rdiff2,rdiff3,term
 real(dp) :: yy
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: igrid(3),irange(3),ngfft(3)
 integer,allocatable :: ii(:,:)
 real(dp) :: drmetds(3,3),lencp(3),rmet(3,3),scale(3),tau(3)
 real(dp),allocatable :: rrdiff(:,:)

! *************************************************************************

! if( ipert<1 .or. ipert> natom+7) then
!   write(message,'(a,i0,a,a,a,i0,a)')&
!&   ' The argument ipert must be between 1 and natom+7=',natom+7,',',ch10,&
!&   ' while it is ipert=',ipert,'.'
!   MSG_BUG(message)
! end if

 if( (ipert==natom+3 .or. ipert==natom+4) .and. cplex/=1) then
   write(message,'(3a,i4,a)')&
&   'The argument cplex must be 1 for strain perturbationh',ch10,&
&   'while it is cplex=',cplex,'.'
   MSG_BUG(message)
 end if

!Zero out array
 xccc3d1(:)=0.0_dp

!For a non-linear XC core correction, the perturbation must be phonon-type or strain type
 if(ipert<=natom .or. ipert==natom+3 .or. ipert==natom+4) then

   if( idir<1 .or. idir> 3) then
     write(message,'(a,a,a,i4,a)')&
&     'The argument idir must be between 1 and 3,',ch10,&
&     'while it is idir=',idir,'.'
     MSG_BUG(message)
   end if

!  Compute lengths of cross products for pairs of primitive
!  translation vectors (used in setting index search range below)
   lencp(1)=cross_mk(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
&   rprimd(1,3),rprimd(2,3),rprimd(3,3))
   lencp(2)=cross_mk(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
&   rprimd(1,1),rprimd(2,1),rprimd(3,1))
   lencp(3)=cross_mk(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
&   rprimd(1,2),rprimd(2,2),rprimd(3,2))

!  Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!  (recall ucvol=R1.(R2xR3))
   scale(:)=ucvol/lencp(:)

!  Compute metric tensor in real space rmet
   do nu=1,3
     rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+&
&     rprimd(3,:)*rprimd(3,nu)
   end do

!  Section to be executed only for strain perturbation
!  Compute derivative of metric tensor wrt strain component istr
   if(ipert==natom+3 .or. ipert==natom+4) then
     istr=idir + 3*(ipert-natom-3)

     ka=idx(2*istr-1);kb=idx(2*istr)
     do jj = 1,3
       drmetds(:,jj)=(rprimd(ka,:)*rprimd(kb,jj)+rprimd(kb,:)*rprimd(ka,jj))
     end do
!    For historical reasons:
     drmetds(:,:)=0.5_dp*drmetds(:,:)

!    end of strain perturbation section
   end if

   ngfft(1)=n1
   ngfft(2)=n2
   ngfft(3)=n3

   delta=1.0_dp/(n1xccc-1)
   deltam1=n1xccc-1
   delta2div6=delta**2/6.0_dp

!  Loop over atoms in unit cell
!  Note that we cycle immediately for all except the displaced atom
!  for such a perturbation.  The loop is executed over all the
!  atoms for a strain peturbation.
   do iatom=1,natom
     if(ipert<=natom .and. iatom/=ipert) cycle
!    Set search range (density cuts off perfectly beyond range)
!    Cycle if no range.
     range=0.0_dp
     range=xcccrc(typat(iatom))
     if(range<1.d-16) cycle

     range2=range**2
     rangem1=1.0_dp/range

!    compute mrange for ii(:,3), rrdiff(:,3), inserted by MM (2005/12/06)
!    Consider each component in turn : compute range
     do mu=1,3

!      Convert reduced coord of given atom to [0,1)
       tau(mu)=mod(xred(mu,iatom)+1._dp-aint(xred(mu,iatom)),1._dp)

!      Use tau to find nearest grid point along R(mu)
!      (igrid=0 is the origin; shift by 1 to agree with usual index)
       igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!      Use range to compute an index range along R(mu)
!      (add 1 to make sure it covers full range)
       irange(mu)=1+nint((range/scale(mu))*dble(ngfft(mu)))

     end do

!    Allocate arrays that depends on the range
     mrange=maxval(irange(1:3))
     ABI_ALLOCATE(ii,(2*mrange+1,3))
     ABI_ALLOCATE(rrdiff,(2*mrange+1,3))

!    Consider each component in turn
     do mu=1,3

!      temporarily suppressed by MM (2005/12/02)
!      Convert reduced coord of given atom to [0,1)
!      tau(mu)=mod(xred(mu,iatom)+1._dp-aint(xred(mu,iatom)),1._dp)

!      Use tau to find nearest grid point along R(mu)
!      (igrid=0 is the origin; shift by 1 to agree with usual index)
!      igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!      Use range to compute an index range along R(mu)
!      (add 1 to make sure it covers full range)
!      irange(mu)=1+nint((range/scale(mu))*dble(ngfft(mu)))

!      Check that the largest range is smallest than the maximum
!      allowed one
!      if(2*irange(mu)+1 > mshift)then
!      write(message, '(a,a,a,a,i6,a)' ) ch10,&
!      &    ' dfpt_mkcore : BUG -',ch10,&
!      &    '  The range around atom',iatom,' is too large.'
!      MSG_BUG(message)
!      end if

!      Set up a counter that explore the relevant range
!      of points around the atom
       ishift=0
       do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
         ishift=ishift+1
         ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
         rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
       end do

!      End loop on mu
     end do

!    Conduct triple loop over restricted range of grid points for iatom

     do ishift3=1,1+2*irange(3)
!      map back to [1,ngfft(3)] for usual fortran index in unit cell
       i3=ii(ishift3,3)
!      find vector from atom location to grid point (reduced)
       rdiff3=rrdiff(ishift3,3)

       do ishift2=1,1+2*irange(2)
         i2=ii(ishift2,2)
         rdiff2=rrdiff(ishift2,2)
!        Prepare the computation of difmag2
         difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&         +2.0_dp*rmet(3,2)*rdiff3*rdiff2
         difmag2_fact=2.0_dp*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

         do ishift1=1,1+2*irange(1)
           rdiff1=rrdiff(ishift1,1)

!          Compute (rgrid-tau-Rprim)**2
           difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!          Only accept contribution inside defined range
           if (difmag2<range2) then

!            Prepare computation of core charge function and derivative,
!            using splines
             difmag=sqrt(difmag2)
             if (difmag>=1.0d-10) then
               i1=ii(ishift1,1)
               yy=difmag*rangem1

!              Compute index of yy over 1 to n1xccc scale
               jj=1+int(yy*(n1xccc-1))
               diff=yy-(jj-1)*delta

!              Will evaluate spline fit (p. 86 Numerical Recipes, Press et al;
!              NOTE error in book for sign of "aa" term in derivative;
!              also see splfit routine).
               bb = diff*deltam1
               aa = 1.0_dp-bb
               cc = aa*(aa**2-1.0_dp)*delta2div6
               dd = bb*(bb**2-1.0_dp)*delta2div6
               
!              Evaluate spline fit of 1st der of core charge density
!              from xccc1d(:,2,:) and (:,4,:)
               func=aa*xccc1d(jj,2,typat(iatom))+bb*xccc1d(jj+1,2,typat(iatom)) +&
&               cc*xccc1d(jj,4,typat(iatom))+dd* xccc1d(jj+1,4,typat(iatom))

               if(ipert<=natom) then
                 phase=2*pi*(qphon(1)*(rdiff1+xred(1,iatom))  &
&                 +qphon(2)*(rdiff2+xred(2,iatom))  &
&                 +qphon(3)*(rdiff3+xred(3,iatom)))
                 prod=rmet(idir,1)*rdiff1+rmet(idir,2)*rdiff2+rmet(idir,3)*rdiff3

                 term=-func*rangem1/difmag*prod
                 ifft=i1+n1*(i2-1+n2*(i3-1))
                 phr=cos(phase)
                 if(cplex==1)then
                   xccc3d1(ifft)=xccc3d1(ifft)+term*phr
                 else
                   phi=sin(phase)
                   xccc3d1(2*ifft-1)=xccc3d1(2*ifft-1)+term*phr
                   xccc3d1(2*ifft  )=xccc3d1(2*ifft  )-term*phi
                 end if
               else
                 prod=&
&                 (rdiff1*(drmetds(1,1)*rdiff1+drmetds(1,2)*rdiff2+drmetds(1,3)*rdiff3)&
&                 +rdiff2*(drmetds(2,1)*rdiff1+drmetds(2,2)*rdiff2+drmetds(2,3)*rdiff3)&
&                 +rdiff3*(drmetds(3,1)*rdiff1+drmetds(3,2)*rdiff2+drmetds(3,3)*rdiff3))
                 term=prod*func*rangem1/difmag
                 
                 ifft=i1+n1*(i2-1+n2*(i3-1))
                 xccc3d1(ifft)=xccc3d1(ifft)+term

               end if

!              End of the condition for the distance not to vanish
             end if

!            End of condition to be inside the range
           end if

!          End loop on ishift1
         end do

!        End loop on ishift2
       end do

!      End loop on ishift3
     end do

     ABI_DEALLOCATE(ii)
     ABI_DEALLOCATE(rrdiff)
!    End loop on atoms
   end do

!  End of the condition ipert corresponds to a phonon type perturbation
!  or strain type perturbation
 end if

 contains

   function cross_mk(xx,yy,zz,aa,bb,cc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cross_mk'
!End of the abilint section

   real(dp) :: cross_mk
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_mk=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_mk

end subroutine dfpt_mkcore
!!***
