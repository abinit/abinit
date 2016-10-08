!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkcore
!! NAME
!! mkcore
!!
!! FUNCTION
!! Optionally compute
!!  (1) core electron density throughout unit cell, or
!!  (2) contribution to Exc gradient wrt xred, or
!!  (3) contribution to stress tensor, or (response function code)
!!  (4) contribution to frozen-wavefunction part of
!!      the dynamical matrix (part 2)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid dimensions.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  option: 1 for computing xccc3d (core charge density),
!!   2 for computing core charge contribution to $d(E_{xc})/d(tau)$,
!!   3 for computing core charge contribution to stress tensor corstr,
!!   4 for contribution to frozen-wavefunction part of dynamical matrix
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real
!!   space--only used when option=2,3, or 4,  else ignored
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  corstr(6)=core charge contribution to stress tensor, only if option=3
!!  dyfrx2(3,3,natom)=non-linear xc core correction part of the
!!    frozen-wavefunction part of the dynamical matrix, only for option=4
!!  grxc(3,natom)=d(Exc)/d(xred), hartree (only computed when option=2, else
!!   ignored)
!!
!! SIDE EFFECTS
!!  xccc3d(n1*n2*n3)=3D core electron density for XC core correction, bohr^-3
!!   (computed and returned when option=1, needed as input when option=3)
!!
!! NOTES
!! Note that this routine is tightly connected to the dfpt_mkcore.f routine
!!
!! PARENTS
!!      dfpt_dyfro,forces,nonlinear,prcref,prcref_PMA,respfn,setvtr,stress
!!
!! CHILDREN
!!      ptabs_fourdp,strconv,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkcore(corstr,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,ntypat,n1,n1xccc,&
& n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_linalg_interfaces

 use m_mpinfo,     only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkcore'
 use interfaces_18_timing
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n1xccc,n2,n3,natom,nfft,nspden,ntypat,option
 real(dp),intent(in) :: ucvol
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rprimd(3,3),vxc(nfft,nspden),xccc1d(n1xccc,6,ntypat)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: xccc3d(nfft)
 real(dp),intent(out) :: corstr(6),dyfrx2(3,3,natom)  !vz_i
 real(dp),intent(inout) :: grxc(3,natom) !vz_i

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,iatom,ier,ifft,ishift,ishift1,ishift2
 integer :: ishift3,itypat,ixp,jj,me_fft,mrange,mu,nfftot,nu
 real(dp) :: dd,delta,delta2div6,deltam1,diff,difmag,difmag2
 real(dp) :: difmag2_fact,difmag2_part,fact,func,grxc1,grxc2,grxc3,range,range2
 real(dp) :: rangem1,rdiff1,rdiff2,rdiff3,strdia,t1,t2,t3,term,term1,term2
 character(len=500) :: message
!arrays
 integer :: igrid(3),irange(3),ngfft(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 integer,allocatable :: ii(:,:)
 real(dp) :: yy,aa,bb,cc
 real(dp) :: corfra(3,3),lencp(3),rmet(3,3),scale(3),tau(3),tsec(2),tt(3)
 real(dp),allocatable :: rrdiff(:,:),work(:,:,:)

!************************************************************************

 call timab(12,1,tsec)

!Make sure option is acceptable
 if (option<0.or.option>4) then
   write(message, '(a,i12,a,a,a)' )&
&   'option=',option,' is not allowed.',ch10,&
&   'Must be 1, 2, 3 or 4.'
   MSG_BUG(message)
 end if

!Zero out only the appropriate array according to option:
!others are dummies with no storage

 if (option==1) then
!  Zero out array to permit accumulation over atom types below:
   xccc3d(:)=zero
 else if (option==2) then
!  Zero out gradient of Exc array
   grxc(:,:)=zero
 else if (option==3) then
!  Zero out locally defined stress array
   corfra(:,:)=zero
   strdia=zero
 else if (option==4) then
!  Zero out fr-wf part of the dynamical matrix
   dyfrx2(:,:,:)=zero
 else
   MSG_BUG(" Can't be here! (bad option)")
 end if

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_mkcore(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_mkcore(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_mkcore(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+rprimd(3,:)*rprimd(3,nu)
 end do

 ngfft(1)=n1
 ngfft(2)=n2
 ngfft(3)=n3
 nfftot=n1*n2*n3
 me_fft = mpi_enreg%me_fft

 ! Get the distrib associated with this fft_grid 
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 delta=one/(n1xccc-1)
 deltam1=n1xccc-1
 delta2div6=delta**2/6.0d0

 if (option>=2) then
   ABI_ALLOCATE(work,(n1,n2,n3))
!  For spin-polarization, replace vxc by (1/2)*(vxc(up)+vxc(down))
!  For non-collinear magnetism, replace vxc by (1/2)*(vxc^{11}+vxc^{22})
   if (nspden>=2) then
     ifft=1
     do i3=1,n3
       if(me_fft==fftn3_distrib(i3)) then 
         do i2=1,n2
           do i1=1,n1
             work(i1,i2,i3)=half*(vxc(ifft,1)+vxc(ifft,2))
             ifft=ifft+1
           end do
         end do
       end if
     end do
   else
     ifft=1
     do i3=1,n3
       if(me_fft==fftn3_distrib(i3)) then 
         do i2=1,n2
           do i1=1,n1
             work(i1,i2,i3)=vxc(ifft,1)
             ifft=ifft+1
           end do
         end do
       end if
     end do
!    call DCOPY(nfft,vxc,1,work,1)
   end if
 end if

!Loop over atoms in unit cell
 do iatom=1,natom

   if(option==2)then
     grxc1=zero
     grxc2=zero
     grxc3=zero
   end if

!  Set search range (density cuts off perfectly beyond range)
   itypat=typat(iatom)
   range=xcccrc(itypat)

!  Skip loop if this atom has no core charge
   if (abs(range)<1.d-16) cycle

   range2=range**2
   rangem1=one/range

!  Consider each component in turn : compute range
   do mu=1,3

!    Convert reduced coord of given atom to [0,1)
     tau(mu)=mod(xred(mu,iatom)+one-aint(xred(mu,iatom)),one)

!    Use tau to find nearest grid point along R(mu)
!    (igrid=0 is the origin; shift by 1 to agree with usual index)
     igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!    Use range to compute an index range along R(mu)
!    (add 1 to make sure it covers full range)
     irange(mu)=1+nint((range/scale(mu))*dble(ngfft(mu)))

   end do

!  Allocate arrays that depends on the range
   mrange=maxval(irange(1:3))
   ABI_ALLOCATE(ii,(2*mrange+1,3))
   ABI_ALLOCATE(rrdiff,(2*mrange+1,3))

!  Set up counters that explore the relevant range
!  of points around the atom
   do mu=1,3
     ishift=0
     do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
       ishift=ishift+1
       ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
       rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
     end do
   end do

!  Conduct triple loop over restricted range of grid points for iatom
   do ishift3=1,1+2*irange(3)
!    map back to [1,ngfft(3)] for usual fortran index in unit cell
     i3=ii(ishift3,3)
     if(fftn3_distrib(i3)/=mpi_enreg%me_fft) cycle
!    find vector from atom location to grid point (reduced)
     rdiff3=rrdiff(ishift3,3)

     do ishift2=1,1+2*irange(2)
       i2=ii(ishift2,2)
       rdiff2=rrdiff(ishift2,2)
!      Prepare the computation of difmag2
       difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&       +2.0d0*rmet(3,2)*rdiff3*rdiff2
       difmag2_fact=2.0d0*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

       do ishift1=1,1+2*irange(1)
         rdiff1=rrdiff(ishift1,1)

!        Compute (rgrid-tau-Rprim)**2
         difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!        Only accept contribution inside defined range
         if (difmag2<range2-tol12) then

!          Prepare computation of core charge function and derivative,
!          using splines
           i1=ii(ishift1,1)
           difmag=sqrt(difmag2)
           yy=difmag*rangem1

!          Compute index of yy over 1 to n1xccc scale
           jj=1+int(yy*(n1xccc-1))
           diff=yy-(jj-1)*delta

!          Will evaluate spline fit (p. 86 Numerical Recipes, Press et al;
!          NOTE error in book for sign of "aa" term in derivative;
!          also see splfit routine).
           bb = diff*deltam1
           aa = one-bb
           cc = aa*(aa**2-one)*delta2div6
           dd = bb*(bb**2-one)*delta2div6


!          Test first for option 2, the most frequently used
           if (option==2) then

!            Accumulate contributions to Exc gradients

             if (difmag>1.0d-10) then

!              Evaluate spline fit of 1st der of core charge density
!              from xccc1d(:,2,:) and (:,4,:)
               func=aa*xccc1d(jj,2,itypat)+bb*xccc1d(jj+1,2,itypat) +&
&               cc*xccc1d(jj,4,itypat)+dd*xccc1d(jj+1,4,itypat)
               term=work(i1,i2,i3)*func/difmag
               grxc1=grxc1+rdiff1*term
               grxc2=grxc2+rdiff2*term
               grxc3=grxc3+rdiff3*term
             end if

           else if (option==1) then

!            Evaluate spline fit of core charge density
!            from xccc1d(:,1,:) and (:,3,:)
             func=aa*xccc1d(jj,1,itypat)+bb*xccc1d(jj+1,1,itypat) +&
&             cc*xccc1d(jj,3,itypat)+dd*xccc1d(jj+1,3,itypat)

!            Accumulate contributions to core electron density
!            throughout unit cell
             ifft=i1+n1*(i2-1+n2*(ffti3_local(i3)-1))
             xccc3d(ifft)=xccc3d(ifft)+func

           else if (option==3) then

!            Accumulate contributions to stress tensor
!            in reduced coordinates

             if (difmag>1.0d-10) then

!              Evaluate spline fit of 1st der of core charge density
!              from xccc1d(:,2,:) and (:,4,:)
               func=aa*xccc1d(jj,2,itypat)+bb*xccc1d(jj+1,2,itypat) +&
&               cc*xccc1d(jj,4,itypat)+dd*xccc1d(jj+1,4,itypat)
               term=work(i1,i2,i3)*func*rangem1/difmag/dble(n1*n2*n3)
!              Write out the 6 symmetric components
               corfra(1,1)=corfra(1,1)+term*rdiff1**2
               corfra(2,2)=corfra(2,2)+term*rdiff2**2
               corfra(3,3)=corfra(3,3)+term*rdiff3**2
               corfra(3,2)=corfra(3,2)+term*rdiff3*rdiff2
               corfra(3,1)=corfra(3,1)+term*rdiff3*rdiff1
               corfra(2,1)=corfra(2,1)+term*rdiff2*rdiff1
!              (the above still needs to be transformed to cartesian coords)

             end if

!            Also compute a diagonal term
!            Evaluate spline fit of core charge density
!            from xccc1d(:,1,:) and (:,3,:)
             func=aa*xccc1d(jj,1,itypat)+bb*xccc1d(jj+1,1,itypat) +&
&             cc*xccc1d(jj,3,itypat)+dd*xccc1d(jj+1,3,itypat)
             strdia=strdia+work(i1,i2,i3)*func

           else if (option==4) then

!            Compute frozen-wf contribution to Dynamical matrix

             tt(1)=rmet(1,1)*rdiff1+rmet(1,2)*rdiff2+rmet(1,3)*rdiff3
             tt(2)=rmet(2,1)*rdiff1+rmet(2,2)*rdiff2+rmet(2,3)*rdiff3
             tt(3)=rmet(3,1)*rdiff1+rmet(3,2)*rdiff2+rmet(3,3)*rdiff3

             if (difmag>1.d-10) then

!              Accumulate contributions to dynamical matrix
               term=(ucvol/dble(nfftot))*work(i1,i2,i3)*rangem1/difmag
!              Evaluate spline fit of 1st der of core charge density
!              from xccc1d(:,2,:) and (:,4,:)
               func=aa*xccc1d(jj,2,itypat)+bb*xccc1d(jj+1,2,itypat) +&
&               cc*xccc1d(jj,4,itypat)+dd*xccc1d(jj+1,4,itypat)
               term1=term*func
!              Evaluate spline fit of 2nd der of core charge density
!              from xccc1d(:,3,:) and (:,5,:)
               func=aa*xccc1d(jj,3,itypat)+bb*xccc1d(jj+1,3,itypat) +&
&               cc*xccc1d(jj,5,itypat)+dd*xccc1d(jj+1,5,itypat)
               term2=term*func*rangem1/difmag
               do mu=1,3
                 do nu=1,3
                   dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)&
&                   +(term2-term1/difmag**2)*tt(mu)*tt(nu)&
&                   +term1*rmet(mu,nu)
                 end do
               end do

             else

!              There is a contribution from difmag=zero !
!              Evaluate spline fit of 2nd der of core charge density
!              from xccc1d(:,3,:) and (:,5,:)
               func=aa*xccc1d(jj,3,itypat)+bb*xccc1d(jj+1,3,itypat) +&
&               cc*xccc1d(jj,5,itypat)+dd*xccc1d(jj+1,5,itypat)
               term=(ucvol/dble(nfftot))*work(i1,i2,i3)*func*rangem1**2
               do mu=1,3
                 do nu=1,3
                   dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)+term*rmet(mu,nu)
                 end do
               end do

!              End of condition not to be precisely on the point (difmag=zero)
             end if

!            If option is not 1, 2, 3, or 4.
           else
             MSG_BUG("Can't be here in mkcore")
!            End of choice of option
           end if

!          End of condition on the range
         end if

!        End loop on ishift1
       end do

!      End loop on ishift2
     end do

!    End loop on ishift3
   end do

   ABI_DEALLOCATE(ii)
   ABI_DEALLOCATE(rrdiff)

   if(option==2)then
     fact=-(ucvol/dble(nfftot))/range
     grxc(1,iatom)=grxc1*fact
     grxc(2,iatom)=grxc2*fact
     grxc(3,iatom)=grxc3*fact

     if(mpi_enreg%nproc_fft == 1) then
       call timab(539,1,tsec)
       call xmpi_sum(grxc1,mpi_enreg%comm_fft,ier)
       call xmpi_sum(grxc2,mpi_enreg%comm_fft,ier)
       call xmpi_sum(grxc3,mpi_enreg%comm_fft,ier)
       call timab(539,2,tsec)
     end if
   end if

!  End big loop on atoms
 end do

 if (option==2) then

!  Apply rmet as needed to get reduced coordinate gradients
   do iatom=1,natom
     t1=grxc(1,iatom)
     t2=grxc(2,iatom)
     t3=grxc(3,iatom)
     grxc(:,iatom)=rmet(:,1)*t1+rmet(:,2)*t2+rmet(:,3)*t3

   end do
 end if

 if (option==3) then

!  Transform stress tensor from full storage mode to symmetric storage mode
   corstr(1)=corfra(1,1)
   corstr(2)=corfra(2,2)
   corstr(3)=corfra(3,3)
   corstr(4)=corfra(3,2)
   corstr(5)=corfra(3,1)
   corstr(6)=corfra(2,1)

!  Transform stress tensor from reduced coordinates to cartesian coordinates
   call strconv(corstr,rprimd,corstr)

!  Compute diagonal contribution to stress tensor (need input xccc3d)
!  strdia = (1/N) Sum(r) [mu_xc_avg(r) * rho_core(r)]
   ifft=0 ; strdia=zero
   do i3=1,n3
     if(me_fft==fftn3_distrib(i3)) then 
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
           strdia=strdia+work(i1,i2,i3)*xccc3d(ifft)
         end do
       end do
     end if
   end do
   strdia=strdia/dble(nfftot)
!  strdia=DDOT(nfft,work,1,xccc3d,1)/dble(nfftot)

!  Add diagonal term to stress tensor
   corstr(1)=corstr(1)+strdia
   corstr(2)=corstr(2)+strdia
   corstr(3)=corstr(3)+strdia
 end if

 if(option>=2)  then
   ABI_DEALLOCATE(work)
 end if

 if(mpi_enreg%nproc_fft > 1) then
   call timab(539,1,tsec)
   if(option==3) then
     call xmpi_sum(corstr,mpi_enreg%comm_fft,ier)
   end if
   if(option==2) then
     call xmpi_sum(grxc,mpi_enreg%comm_fft,ier)
   end if
   call timab(539,2,tsec)
 end if

 call timab(12,2,tsec)

 contains

   real(dp) pure function cross_mkcore(xx,yy,zz,aa,bb,cc)
!    Define magnitude of cross product of two vectors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cross_mkcore'
!End of the abilint section

   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_mkcore=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_mkcore

end subroutine mkcore
!!***
