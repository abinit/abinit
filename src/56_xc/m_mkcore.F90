!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mkcore
!! NAME
!!  m_mkcore
!!
!! FUNCTION
!!  Routines related to non-linear core correction.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, TRangel, MT)
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

module m_mkcore

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_linalg_interfaces

 use defs_abitypes, only : mpi_type
 use m_geometry,   only : strconv
 use m_time,       only : timab
 use m_mpinfo,     only : ptabs_fourdp
 use m_sort,        only : sort_dp
 use m_pawrad,      only : pawrad_type, pawrad_init, pawrad_free
 use m_pawtab,      only : pawtab_type
 use m_paw_numeric, only : paw_splint

 implicit none

 private
!!***

 public :: mkcore
 public :: mkcore_alt
 public :: dfpt_mkcore       ! Derivative of the core electron density with respect to one specific atom displacement
!!***

contains
!!***

!!****f* ABINIT/mkcore
!! NAME
!! mkcore
!!
!! FUNCTION
!! Optionally compute:
!!  (1) pseudo core electron density throughout unit cell
!!  (2) pseudo-core contribution to forces
!!  (3) pseudo-core contribution to stress tensor
!!  (4) pseudo-core contrib. to frozen-wf part the dynamical matrix (part 2)
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
!!      xchybrid_ncpp_cc
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkcore(corstr,dyfrx2,grxc,mpi_enreg,natom,nfft,nspden,ntypat,n1,n1xccc,&
& n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)

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
 real(dp),intent(out) :: corstr(6),dyfrx2(3,3,natom)
 real(dp),intent(inout) :: grxc(3,natom)

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
    'option=',option,' is not allowed.',ch10,&
    'Must be 1, 2, 3 or 4.'
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

         end if ! End of condition on the range
       end do ! End loop on ishift1
     end do ! End loop on ishift2
   end do ! End loop on ishift3

   ABI_DEALLOCATE(ii)
   ABI_DEALLOCATE(rrdiff)

   if(option==2)then
     fact=-(ucvol/dble(nfftot))/range
     grxc(1,iatom)=grxc1*fact
     grxc(2,iatom)=grxc2*fact
     grxc(3,iatom)=grxc3*fact
   end if

 end do !  End big loop on atoms

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
   if(option==2) then
     call xmpi_sum(grxc,mpi_enreg%comm_fft,ier)
   end if
   if(option==3) then
     call xmpi_sum(corstr,mpi_enreg%comm_fft,ier)
   end if
   if(option==4) then
     call xmpi_sum(dyfrx2,mpi_enreg%comm_fft,ier)
   end if
   call timab(539,2,tsec)
 end if

 call timab(12,2,tsec)

 contains
!!***

!!****f* ABINIT/cross_mkcore
!! NAME
!!  cross_mkcore
!!
!! FUNCTION
!!  Define magnitude of cross product of two vectors
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

   function cross_mkcore(xx,yy,zz,aa,bb,cc)

   real(dp) :: cross_mkcore
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
! *************************************************************************
   cross_mkcore=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_mkcore

end subroutine mkcore
!!***

!--------------------------------------------------------------------------------------

!!****f* ABINIT/mkcore_alt
!! NAME
!! mkcore_alt
!!
!! FUNCTION
!! Optionally compute:
!!  (1) pseudo core electron density throughout unit cell
!!  (2) pseudo-core contribution to forces
!!  (3) pseudo-core contribution to stress tensor
!!  (4) pseudo-core contrib. to frozen-wf part the dynamical matrix (part 2)
!! This routine is an alternative to mkcore, to be used for PAW and/or WVL.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  icoulomb= periodic treatment of Hartree potential: 0=periodic, 1=free BC, 2=surface BC
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell
!!  n1,n2,n3=fft grid dimensions.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  option: 1 for computing core charge density
!!          2 for computing core charge contribution to forces
!!          3 for computing core charge contribution to stress tensor
!!          4 for computing contribution to frozen-wf part of dynamical matrix
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  ucvol=unit cell volume (bohr**3)
!!  usepaw=flag for PAW method
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and 5 derivatives for each atom type
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  === if option==1 ===
!!  xccc3d(n1*n2*n3)=3D core electron density for XC core correction (bohr^-3)
!!  === if option==2 ===
!!  grxc(3,natom)=core charge contribution to forces
!!  === if option==3 ===
!!  corstr(6)=core charge contribution to stress tensor
!!  === if option==4 ===
!!  dyfrx2(3,3,natom)=non-linear xc core correction part of the
!!    frozen-wavefunction part of the dynamical matrix
!!
!! SIDE EFFECTS
!!  xccc3d(n1*n2*n3)=3D core electron density for XC core correction (bohr^-3)
!!   (computed and returned when option=1, needed as input when option=3)
!!
!! NOTES
!!  Based on mkcore.F90
!!
!! PARENTS
!!      forces,setvtr,stress
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkcore_alt(atindx1,corstr,dyfrx2,grxc,icoulomb,mpi_enreg,natom,nfft,nspden,&
&          nattyp,ntypat,n1,n1xccc,n2,n3,option,rprimd,ucvol,vxc,xcccrc,xccc1d,&
&          xccc3d,xred,pawrad,pawtab,usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icoulomb,n1,n1xccc,n2,n3,natom,nfft,nspden,ntypat,option,usepaw
 real(dp),intent(in) :: ucvol
 type(mpi_type),intent(in) :: mpi_enreg
 type(pawrad_type),intent(in) :: pawrad(:)
 type(pawtab_type),intent(in) :: pawtab(:)
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(in),target :: vxc(nfft,nspden)
 real(dp),intent(out) :: corstr(6),grxc(3,natom),dyfrx2(3,3,natom)
 real(dp),intent(inout) :: xccc3d(nfft)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,iat,iatm,iatom,ier,ipts
 integer :: ishift,ishift1,ishift2,ishift3
 integer :: itypat,ixp,jj,jpts,me_fft,mrange,mu
 integer :: nfftot,npts,npts12,nu
 logical :: letsgo
 real(dp) :: aa,bb,cc,dd,delta,delta2div6,deltam1
 real(dp) :: diff,difmag,fact,range,range2
 real(dp) :: rangem1,rdiff1,rdiff2,rdiff3
 real(dp) :: rnorm2,rnorm2_fact,rnorm2_part
 real(dp) :: strdia,t1,t2,t3,term,term1,term2,yy
 character(len=1) :: geocode
 character(len=500) :: message
 type(pawrad_type) :: core_mesh
!arrays
 integer :: igrid(3),irange(3),ishiftmax(3),ngfft(3)
 integer,allocatable :: ii(:,:),iindex(:),indx1(:),indx2(:)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 logical :: per(3)
 real(dp) :: corfra(3,3),corgr(3),lencp(3),rmet(3,3)
 real(dp) :: scale(3),tau(3),tsec(2),tt(3)
 real(dp),allocatable :: dtcore(:),d2tcore(:),rnorm(:)
 real(dp),allocatable :: rrdiff(:,:),tcore(:)
 real(dp), ABI_CONTIGUOUS pointer :: vxc_eff(:)

!************************************************************************

 call timab(12,1,tsec)

!Make sure option is acceptable
 if (option<0.or.option>4) then
   write(message, '(a,i12,a,a,a)' )&
    'option=',option,' is not allowed.',ch10,&
    'Must be 1, 2, 3 or 4.'
   MSG_BUG(message)
 end if

!Zero out only the appropriate array according to option:
 if (option==1) then
   xccc3d(:)=zero
 else if (option==2) then
   grxc(:,:)=zero
 else if (option==3) then
   corfra(:,:)=zero
   strdia=zero
 else if (option==4) then
   dyfrx2(:,:,:)=zero
 end if

!Conditions for periodicity in the three directions
 geocode='P'
 if (icoulomb==1) geocode='F'
 if (icoulomb==2) geocode='S'
 per(1)=(geocode /= 'F')
 per(2)=(geocode == 'P')
 per(3)=(geocode /= 'F')

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_mkcore_alt(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_mkcore_alt(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_mkcore_alt(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+rprimd(3,:)*rprimd(3,nu)
 end do

!Get the distrib associated with this fft_grid
 ngfft(1)=n1;ngfft(2)=n2;ngfft(3)=n3
 nfftot=n1*n2*n3 ; me_fft=mpi_enreg%me_fft
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 delta=one/(n1xccc-1)
 deltam1=n1xccc-1
 delta2div6=delta**2/6.0_dp

 if (option>=2) then
!  For spin-polarization, replace vxc by (1/2)*(vxc(up)+vxc(down))
!  For non-collinear magnetism, replace vxc by (1/2)*(vxc^{11}+vxc^{22})
   if (nspden>=2) then
     ABI_ALLOCATE(vxc_eff,(nfft))
     do jj=1,nfft
       vxc_eff(jj)=half*(vxc(jj,1)+vxc(jj,2))
     end do
   else
     vxc_eff => vxc(1:nfft,1)
   end if
 end if

!Loop over atom types
 iatm=0
 do itypat=1,ntypat

!  Set search range (density cuts off perfectly beyond range)
   range=xcccrc(itypat);if (usepaw==1) range=pawtab(itypat)%rcore
   range2=range**2 ; rangem1=one/range

!  Skip loop if this type has no core charge
   if (abs(range)<1.d-16) cycle

!  PAW: create mesh for core density
   if (usepaw==1) then
     call pawrad_init(core_mesh,mesh_size=pawtab(itypat)%core_mesh_size,&
&     mesh_type=pawrad(itypat)%mesh_type,&
&     rstep=pawrad(itypat)%rstep,lstep=pawrad(itypat)%lstep)
   end if

!  Loop over atoms of the type
   do iat=1,nattyp(itypat)
     iatm=iatm+1;iatom=atindx1(iatm)

     if(option==2) corgr(:)=zero

!    Consider each component in turn : compute range
     do mu=1,3
!      Convert reduced coord of given atom to [0,1)
       tau(mu)=mod(xred(mu,iatom)+one-aint(xred(mu,iatom)),one)
!      Use tau to find nearest grid point along R(mu)
!      (igrid=0 is the origin; shift by 1 to agree with usual index)
       igrid(mu)=nint(tau(mu)*real(ngfft(mu),dp))
!      Use range to compute an index range along R(mu)
!      (add 1 to make sure it covers full range)
       irange(mu)=1+nint((range/scale(mu))*real(ngfft(mu),dp))
     end do

!    Allocate arrays that depends on the range
     mrange=maxval(irange(1:3))
     ABI_ALLOCATE(ii,(2*mrange+1,3))
     ABI_ALLOCATE(rrdiff,(2*mrange+1,3))

!    Set up counters that explore the relevant range of points around the atom
     if (geocode=='P') then
!      Fully periodic version
       do mu=1,3
         ishift=0
         do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
           ishift=ishift+1
           ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
           rrdiff(ishift,mu)=real(ixp,dp)/real(ngfft(mu),dp)-tau(mu)
         end do
         ishiftmax(mu)=ishift
       end do
     else
!      Free or surface conditions
       do mu=1,3
         ishift=0
         do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
           call indpos_mkcore_alt(per(mu),ixp,ngfft(mu),jj,letsgo)
           if (letsgo) then
             ishift=ishift+1;ii(ishift,mu)=1+jj
             rrdiff(ishift,mu)=real(ixp,dp)/real(ngfft(mu),dp)-tau(mu)
           end if
         end do
         ishiftmax(mu)=ishift
       end do
     end if
     npts12=ishiftmax(1)*ishiftmax(2)
     ABI_ALLOCATE(indx1,(npts12))
     ABI_ALLOCATE(indx2,(npts12))
     ABI_ALLOCATE(iindex,(npts12))
     ABI_ALLOCATE(rnorm,(npts12))
     if (option==1.or.option==3) then
       ABI_ALLOCATE(tcore,(npts12))
     end if
     if (option>=2) then
       ABI_ALLOCATE(dtcore,(npts12))
     end if
     if (option==4) then
       ABI_ALLOCATE(d2tcore,(npts12))
     end if

!    Conduct loop over restricted range of grid points for iatom
     do ishift3=1,ishiftmax(3)
       i3=ii(ishift3,3)
       rdiff3=rrdiff(ishift3,3)

       if(fftn3_distrib(i3)/=mpi_enreg%me_fft) cycle

!      Select the vectors located around the current atom
!        TR: all of the following  could be done inside or
!        outside the loops (i2,i1,i3).
!        Outside: the memory consumption increases.
!        Inside: the time of calculation increases.
!        Here, I choose to do it here, somewhere in the middle.
       npts=0
       do ishift2=1,ishiftmax(2)
         i2=ii(ishift2,2) ; rdiff2=rrdiff(ishift2,2)
         rnorm2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2 &
&         +2.0d0*rmet(3,2)*rdiff3*rdiff2
         rnorm2_fact=2.0d0*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)
         do ishift1=1,ishiftmax(1)
           i1=ii(ishift1,1) ; rdiff1=rrdiff(ishift1,1)
           rnorm2=rnorm2_part+rdiff1*(rnorm2_fact+rmet(1,1)*rdiff1)
!          Only accept contributions inside defined range
           if (rnorm2<range2-tol12) then
             npts=npts+1 ; iindex(npts)=npts
             indx1(npts)=ishift1;indx2(npts)=ishift2
             rnorm(npts)=sqrt(rnorm2)
           end if
         end do
       end do
       if (npts==0) cycle
       if (npts>npts12) then
         message='npts>npts12!'
         MSG_BUG(message)
       end if

!      Evaluate core density (and derivatives) on the set of selected points
       if (usepaw==1) then
!        PAW: use splint routine
         call sort_dp(npts,rnorm(1:npts),iindex(1:npts),tol16)
         if (option==1.or.option==3) then
!          Evaluate fit of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,1), &
&           pawtab(itypat)%tcoredens(:,3),&
&           npts,rnorm(1:npts),tcore(1:npts))
         end if
         if (option>=2) then
!          Evaluate fit of 1-der of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,2), &
&           pawtab(itypat)%tcoredens(:,4),&
&           npts,rnorm(1:npts),dtcore(1:npts))
         end if
         if (option==4) then
!          Evaluate fit of 2nd-der of core density
           call paw_splint(core_mesh%mesh_size,core_mesh%rad, &
&           pawtab(itypat)%tcoredens(:,3), &
&           pawtab(itypat)%tcoredens(:,5),&
&           npts,rnorm(1:npts),d2tcore(1:npts))
         end if
       else
!        Norm-conserving PP:
!          Evaluate spline fit with method from Numerical Recipes
!          (p. 86 Numerical Recipes, Press et al;
!          NOTE error in book for sign of "aa" term in derivative)
         do ipts=1,npts
           yy=rnorm(ipts)*rangem1
           jj=1+int(yy*(n1xccc-1))
           diff=yy-(jj-1)*delta
           bb = diff*deltam1 ; aa = one-bb
           cc = aa*(aa**2-one)*delta2div6
           dd = bb*(bb**2-one)*delta2div6
           if (option==1.or.option==3) then
             tcore(ipts)=aa*xccc1d(jj,1,itypat)+bb*xccc1d(jj+1,1,itypat) +&
&             cc*xccc1d(jj,3,itypat)+dd*xccc1d(jj+1,3,itypat)
           end if
           if (option>=2) then
             dtcore(ipts)=aa*xccc1d(jj,2,itypat)+bb*xccc1d(jj+1,2,itypat) +&
&             cc*xccc1d(jj,4,itypat)+dd*xccc1d(jj+1,4,itypat)
           end if
           if (option==4) then
             d2tcore(ipts)=aa*xccc1d(jj,3,itypat)+bb*xccc1d(jj+1,3,itypat) +&
&             cc*xccc1d(jj,5,itypat)+dd*xccc1d(jj+1,5,itypat)
           end if
         end do
       end if

!      Now, perform the loop over selected grid points
       do ipts=1,npts
         ishift1=indx1(iindex(ipts))
         ishift2=indx2(iindex(ipts))
         difmag=rnorm(ipts)

         rdiff1=rrdiff(ishift1,1);rdiff2=rrdiff(ishift2,2)
         jpts=ii(ishift1,1)+n1*(ii(ishift2,2)-1+n2*(ffti3_local(i3)-1))

!        === Evaluate charge density
         if (option==1) then
           xccc3d(jpts)=xccc3d(jpts)+tcore(ipts)

!        === Accumulate contributions to forces
         else if (option==2) then
           if (difmag>tol10) then
             term=vxc_eff(jpts)*dtcore(ipts)/difmag
             corgr(1)=corgr(1)+rdiff1*term
             corgr(2)=corgr(2)+rdiff2*term
             corgr(3)=corgr(3)+rdiff3*term
           end if

!        === Accumulate contributions to stress tensor (in red. coordinates)
         else if (option==3) then
           if (difmag>tol10) then
             term=vxc_eff(jpts)*dtcore(ipts)*rangem1/difmag/real(nfftot,dp)
!            Write out the 6 symmetric components
             corfra(1,1)=corfra(1,1)+term*rdiff1*rdiff1
             corfra(2,2)=corfra(2,2)+term*rdiff2*rdiff2
             corfra(3,3)=corfra(3,3)+term*rdiff3*rdiff3
             corfra(3,2)=corfra(3,2)+term*rdiff3*rdiff2
             corfra(3,1)=corfra(3,1)+term*rdiff3*rdiff1
             corfra(2,1)=corfra(2,1)+term*rdiff2*rdiff1
!            (the above still needs to be transformed to cartesian coords)
           end if
!          Also compute a diagonal term
           strdia=strdia+vxc_eff(jpts)*tcore(ipts)

!        === Compute frozen-wf contribution to Dynamical matrix
         else if (option==4) then
           tt(1)=rmet(1,1)*rdiff1+rmet(1,2)*rdiff2+rmet(1,3)*rdiff3
           tt(2)=rmet(2,1)*rdiff1+rmet(2,2)*rdiff2+rmet(2,3)*rdiff3
           tt(3)=rmet(3,1)*rdiff1+rmet(3,2)*rdiff2+rmet(3,3)*rdiff3
           if (difmag>tol10) then
             term=(ucvol/real(nfftot,dp))*vxc_eff(jpts)*rangem1/difmag
             term1=term*tcore(ipts)
             term2=term*d2tcore(ipts)*rangem1/difmag
             do mu=1,3
               do nu=1,3
                 dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)&
&                 +(term2-term1/difmag**2)*tt(mu)*tt(nu)&
&                 +term1*rmet(mu,nu)
               end do
             end do
           else
!            There is a contribution from difmag=zero !
             term=(ucvol/real(nfftot,dp))*vxc_eff(jpts)*d2tcore(ipts)*rangem1**2
             do mu=1,3
               do nu=1,3
                 dyfrx2(mu,nu,iatom)=dyfrx2(mu,nu,iatom)+term*rmet(mu,nu)
               end do
             end do
           end if
         end if ! Choice of option

       end do ! Loop on ipts (ishift1, ishift2)

     end do ! Loop on ishift3

     ABI_DEALLOCATE(ii)
     ABI_DEALLOCATE(rrdiff)
     ABI_DEALLOCATE(indx1)
     ABI_DEALLOCATE(indx2)
     ABI_DEALLOCATE(iindex)
     ABI_DEALLOCATE(rnorm)
     if (allocated(tcore)) then
       ABI_DEALLOCATE(tcore)
     end if
     if (allocated(dtcore)) then
       ABI_DEALLOCATE(dtcore)
     end if
     if (allocated(d2tcore)) then
       ABI_DEALLOCATE(d2tcore)
     end if

     if (option==2) then
       fact=-(ucvol/real(nfftot,dp))/range
       grxc(:,iatom)=corgr(:)*fact
     end if

!  End loop on atoms
   end do

   if (usepaw==1) then
     call pawrad_free(core_mesh)
   end if

!End loop over atom types
 end do

 if(option>=2.and.nspden>=2)  then
   ABI_DEALLOCATE(vxc_eff)
 end if

!Forces: translate into reduced coordinates
 if (option==2) then
   do iatom=1,natom
     t1=grxc(1,iatom);t2=grxc(2,iatom);t3=grxc(3,iatom)
     grxc(:,iatom)=rmet(:,1)*t1+rmet(:,2)*t2+rmet(:,3)*t3
   end do
 end if

!Stress tensor: symmetrize, translate into cartesian coord., add diagonal part
 if (option==3) then
   corstr(1)=corfra(1,1) ; corstr(2)=corfra(2,2)
   corstr(3)=corfra(3,3) ; corstr(4)=corfra(3,2)
   corstr(5)=corfra(3,1) ; corstr(6)=corfra(2,1)
   call strconv(corstr,rprimd,corstr)
   corstr(1)=corstr(1)+strdia/real(nfftot,dp)
   corstr(2)=corstr(2)+strdia/real(nfftot,dp)
   corstr(3)=corstr(3)+strdia/real(nfftot,dp)
 end if

!If needed sum over MPI processes
 if(mpi_enreg%nproc_fft>1) then
   call timab(539,1,tsec)
   if (option==2) then
     call xmpi_sum(grxc,mpi_enreg%comm_fft,ier)
   end if
   if (option==3) then
     call xmpi_sum(corstr,mpi_enreg%comm_fft,ier)
   end if
   if (option==4) then
     call xmpi_sum(dyfrx2,mpi_enreg%comm_fft,ier)
   end if
   call timab(539,2,tsec)
 end if

 call timab(12,2,tsec)

 contains
!!***

!--------------------------------------------------------------

!!****f* ABINIT/cross_mkcore_alt
!! NAME
!!  cross_mkcore_alt
!!
!! FUNCTION
!!  Define magnitude of cross product of two vectors
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

   function cross_mkcore_alt(xx,yy,zz,aa,bb,cc)

    real(dp) :: cross_mkcore_alt
    real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
! *************************************************************************
   cross_mkcore_alt=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_mkcore_alt
!!***

!--------------------------------------------------------------

!!****f* ABINIT/indpos_mkcore_alt
!! NAME
!!  indpos_mkcore_alt
!!
!! FUNCTION
!!  Find the grid index of a given position in the cell according to the BC
!!  Determine also whether the index is inside or outside the box for free BC
!!
!! PARENTS
!!      mkcore
!!
!! CHILDREN
!!
!! SOURCE

   subroutine indpos_mkcore_alt(periodic,ii,nn,jj,inside)
!    Find the grid index of a given position in the cell according to the BC
!    Determine also whether the index is inside or outside the box for free BC
    integer, intent(in) :: ii,nn
    integer, intent(out) :: jj
    logical, intent(in) :: periodic
    logical, intent(out) :: inside
! *************************************************************************
   if (periodic) then
     inside=.true. ; jj=modulo(ii-1,nn)+1
   else
     jj=ii ; inside=(ii>=1.and.ii<=nn)
   end if
 end subroutine indpos_mkcore_alt

end subroutine mkcore_alt
!!***

!!****f* ABINIT/dfpt_mkcore
!! NAME
!! dfpt_mkcore
!!
!! FUNCTION
!! Compute the derivative of the core electron density
!! with respect to one specific atom displacement
!! In case of derivative with respect to k or
!! electric (magnetic) field perturbation, the 1st-order core electron density
!! vanishes.
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

subroutine dfpt_mkcore(cplex,idir,ipert,natom,ntypat,n1,n1xccc,&
& n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xccc3d1,xred)

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

   real(dp) :: cross_mk
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_mk=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_mk

end subroutine dfpt_mkcore
!!***

end module m_mkcore
!!***
