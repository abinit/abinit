!{\src2tex{textfont=tt}}
!!****f* ABINIT/eltxccore
!! NAME
!! eltxccore
!!
!! FUNCTION
!! Compute the core charge contributions to the 2nd derivatives
!! of the exchange-correlation energy with respect to all pairs of
!! strain or strain and atomic displacement for the frozen wavefunction
!! contribution to the elastic tensor. 1st-order potentials representing
!! the perturbation by one strain are supplied, and the routine loops
!! over the second strain and over all atomic displacements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid dimensions.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  vxc_core(nfft)=spin-averaged xc potential
!!  vxc10_core(nfft)=spin-averaged 1st-order xc potential for elastic tensor
!!  vxc1is_core(nfft)=spin-averaged 1st-order xc potential for internal strain
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic tensor
!!
!! SIDE EFFECTS
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic and internal-strain tensor.  One column is incremented
!!   by the core contribution.
!!
!! NOTES
!! Note that this routine is related to the mkcore.f routine
!!
!! PARENTS
!!      dfpt_eltfrxc
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine eltxccore(eltfrxc,is2_in,my_natom,natom,nfft,ntypat,&
& n1,n1xccc,n2,n3,rprimd,typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&
& xcccrc,xccc1d,xred, &
& mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_time,       only : timab
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_xmpi,       only : xmpi_comm_self,xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eltxccore'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: is2_in,n1,n1xccc,n2,n3,my_natom,natom,nfft,ntypat
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rprimd(3,3),vxc10_core(nfft),vxc1is_core(nfft)
 real(dp),intent(in) :: vxc_core(nfft),xccc1d(n1xccc,6,ntypat)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: eltfrxc(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: mshift=401
 integer :: i1,i2,i3,iat,iatom,ierr,ifft,is1,is2,ishift,ishift1,ishift2
 integer :: ishift3,ixp,jj,js,ka,kb,kd,kg,mu,my_comm_atom,nu
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: aa,bb,cc,d2rss,dd,delta,delta2div6,deltam1,diff
 real(dp) :: difmag,difmag2,difmag2_fact,difmag2_part,drss1,drss2,func1
 real(dp) :: func2,range,range2,rangem1,rdiff1,rdiff2,rdiff3
 real(dp) :: term1,term2,yy
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: igrid(3),ii(mshift,3),irange(3),ngfft(3)
 integer,pointer :: my_atmtab(:)
 real(dp) :: drm(3,3,6),eltfrxc_core(6+3*natom,6),lencp(3),rmet(3,3),rrdiff(mshift,3)
 real(dp) :: scale(3),tau(3),ts2(3),tsec(2),tt(3)
 real(dp),allocatable :: d2rm(:,:,:,:)

! *************************************************************************

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_elt(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_elt(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_elt(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+&
&   rprimd(3,:)*rprimd(3,nu)
 end do

!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.

 ABI_ALLOCATE(d2rm,(3,3,6,6))

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     drm(:,jj,is2)=rprimd(kg,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kg,jj)
   end do

!  Loop over 1st strain index
   do is1=1,6

     ka=idx(2*is1-1);kb=idx(2*is1)
     d2rm(:,:,is1,is2)=0._dp
     do jj = 1,3
       if(ka==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kb,jj)
       if(ka==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(kb,jj)
       if(kb==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(ka,jj)
       if(kb==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(ka,jj)
     end do
   end do !is1
 end do !is2

 ngfft(1)=n1
 ngfft(2)=n2
 ngfft(3)=n3
 delta=1.0_dp/(n1xccc-1)
 deltam1=n1xccc-1
 delta2div6=delta**2/6.0_dp

!Loop over atoms in unit cell
 eltfrxc_core(:,:)=zero

 do iat=1,my_natom
   iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
   js=7+3*(iatom-1)
!  Set search range (density cuts off perfectly beyond range)
!  Cycle if no range.
   range=0.0_dp
   range=xcccrc(typat(iatom))
   if(range<1.d-16) cycle

   range2=range**2
   rangem1=1.0_dp/range

!  Consider each component in turn
   do mu=1,3
     tau(mu)=mod(xred(mu,iatom)+1._dp-aint(xred(mu,iatom)),1._dp)

!    Use tau to find nearest grid point along R(mu)
!    (igrid=0 is the origin; shift by 1 to agree with usual index)
     igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!    Use range to compute an index range along R(mu)
!    (add 1 to make sure it covers full range)
     irange(mu)=1+nint((range/scale(mu))*dble(ngfft(mu)))

!    Check that the largest range is smallest than the maximum
!    allowed one
     if(2*irange(mu)+1 > mshift)then
       write(message, '(a,i0,a)' )' The range around atom',iatom,' is too large.'
       MSG_BUG(message)
     end if

!    Set up a counter that explore the relevant range
!    of points around the atom
     ishift=0
     do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
       ishift=ishift+1
       ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
       rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
     end do

!    End loop on mu
   end do

!  Conduct triple loop over restricted range of grid points for iatom

   do ishift3=1,1+2*irange(3)
!    map back to [1,ngfft(3)] for usual fortran index in unit cell
     i3=ii(ishift3,3)
!    find vector from atom location to grid point (reduced)
     rdiff3=rrdiff(ishift3,3)

     do ishift2=1,1+2*irange(2)
       i2=ii(ishift2,2)
       rdiff2=rrdiff(ishift2,2)
!      Prepare the computation of difmag2
       difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&       +2.0_dp*rmet(3,2)*rdiff3*rdiff2
       difmag2_fact=2.0_dp*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

       do ishift1=1,1+2*irange(1)
         rdiff1=rrdiff(ishift1,1)
         
!        Compute (rgrid-tau-Rprim)**2
         difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!        Only accept contribution inside defined range
         if (difmag2<range2) then

!          Prepare computation of core charge function and derivatives,
!          using splines
           difmag=sqrt(difmag2)
           if (difmag>=1.0d-10) then
             i1=ii(ishift1,1)
             yy=difmag*rangem1

!            Compute index of yy over 1 to n1xccc scale
             jj=1+int(yy*(n1xccc-1))
             diff=yy-(jj-1)*delta

!            Will evaluate spline fit (p. 86 Numerical Recipes, Press et al;
!            NOTE error in book for sign of "aa" term in derivative;
!            also see splfit routine).
             bb = diff*deltam1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*delta2div6
             dd = bb*(bb**2-1.0_dp)*delta2div6
             
!            Evaluate spline fit of 1st der of core charge density
!            from xccc1d(:,2,:) and (:,4,:)
             func1=aa*xccc1d(jj,2,typat(iatom))+bb*xccc1d(jj+1,2,typat(iatom)) +&
&             cc*xccc1d(jj,4,typat(iatom))+dd*xccc1d(jj+1,4,typat(iatom))
             term1=func1*rangem1
!            Evaluate spline fit of 2nd der of core charge density
!            from xccc1d(:,3,:) and (:,5,:)
             func2=aa*xccc1d(jj,3,typat(iatom))+bb*xccc1d(jj+1,3,typat(iatom)) +&
&             cc*xccc1d(jj,5,typat(iatom))+dd*xccc1d(jj+1,5,typat(iatom))
             term2=func2*rangem1**2
             
             ifft=i1+n1*(i2-1+n2*(i3-1))
             tt(:)=rmet(:,1)*rdiff1+rmet(:,2)*rdiff2+rmet(:,3)*rdiff3
             
!            Add contributions to 2nd derivative tensor
             drss2=&
&             (rdiff1*(drm(1,1,is2_in)*rdiff1+drm(1,2,is2_in)*rdiff2&
&             +drm(1,3,is2_in)*rdiff3)&
&             +rdiff2*(drm(2,1,is2_in)*rdiff1+drm(2,2,is2_in)*rdiff2&
&             +drm(2,3,is2_in)*rdiff3)&
&             +rdiff3*(drm(3,1,is2_in)*rdiff1+drm(3,2,is2_in)*rdiff2&
&             +drm(3,3,is2_in)*rdiff3))

!            Loop over 1st strain index
             do is1=1,6
               
               drss1=&
&               (rdiff1*(drm(1,1,is1)*rdiff1+drm(1,2,is1)*rdiff2&
&               +drm(1,3,is1)*rdiff3)&
&               +rdiff2*(drm(2,1,is1)*rdiff1+drm(2,2,is1)*rdiff2&
&               +drm(2,3,is1)*rdiff3)&
&               +rdiff3*(drm(3,1,is1)*rdiff1+drm(3,2,is1)*rdiff2&
&               +drm(3,3,is1)*rdiff3))
               
               d2rss=&
&               (rdiff1*(d2rm(1,1,is1,is2_in)*rdiff1+d2rm(1,2,is1,is2_in)*rdiff2&
&               +d2rm(1,3,is1,is2_in)*rdiff3)&
&               +rdiff2*(d2rm(2,1,is1,is2_in)*rdiff1+d2rm(2,2,is1,is2_in)*rdiff2&
&               +d2rm(2,3,is1,is2_in)*rdiff3)&
&               +rdiff3*(d2rm(3,1,is1,is2_in)*rdiff1+d2rm(3,2,is1,is2_in)*rdiff2&
&               +d2rm(3,3,is1,is2_in)*rdiff3))

!              Vall(0) X Rhocore(2) term
               eltfrxc_core(is1,is2_in)=eltfrxc_core(is1,is2_in)+0.25_dp*&
&               (vxc_core(ifft)*(term1*(d2rss/difmag&
&               -drss1*drss2/difmag**3)&
&               +term2*drss1*drss2/difmag**2))

!              Vall(1) X Rhocore(1) term
               eltfrxc_core(is1,is2_in)=eltfrxc_core(is1,is2_in)+0.25_dp*&
&               vxc10_core(ifft)*drss1*term1/difmag
               eltfrxc_core(is2_in,is1)=eltfrxc_core(is2_in,is1)+0.25_dp*&
&               vxc10_core(ifft)*drss1*term1/difmag
               
!              End loop in is1
             end do
!            Internal strain contributions
             ts2(:)=drm(:,1,is2_in)*rdiff1+drm(:,2,is2_in)*rdiff2&
&             +drm(:,3,is2_in)*rdiff3

             eltfrxc_core(js:js+2,is2_in)=eltfrxc_core(js:js+2,is2_in)&
&             -(vxc1is_core(ifft)*term1/difmag&
&             +0.5_dp*vxc_core(ifft)*(term2-term1/difmag)*drss2/difmag**2)*tt(:)&
&             -(vxc_core(ifft)*term1/difmag)*ts2(:)

!            End of the condition for the distance not to vanish
           end if

!          End of condition to be inside the range
         end if
         
!        End loop on ishift1
       end do

!      End loop on ishift2
     end do

!    End loop on ishift3
   end do

!  End loop on atoms
 end do

!In case of parallelism over atoms: communicate
 if (paral_atom) then
   call timab(48,1,tsec)
   call xmpi_sum(eltfrxc_core,my_comm_atom,ierr)
   call timab(48,2,tsec)
 end if 

!Add core contribution to XC elastic tensor
 eltfrxc(:,:)=eltfrxc(:,:)+eltfrxc_core(:,:)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 ABI_DEALLOCATE(d2rm)

 contains

   function cross_elt(xx,yy,zz,aa,bb,cc)
!Define magnitude of cross product of two vectors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cross_elt'
!End of the abilint section

   real(dp) :: cross_elt
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_elt=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_elt

end subroutine eltxccore
!!***
