!{\src2tex{textfont=tt}}
!!****f* ABINIT/vlocalstr
!! NAME
!! vlocalstr
!!
!! FUNCTION
!! Compute strain derivatives of local ionic potential
!!                second derivative of E wrt xred
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric ($\textrm{Bohr}^{-2}$).
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on $|G|^2$: see setup1 for definition (doubled sphere).
!!  istr=1,...6 specifies cartesian strain component 11,22,33,32,31,21
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  ucvol=unit cell volume ($\textrm{Bohr}^3$).
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!
!! OUTPUT
!!  vpsp1(nfft)=first-order local crystal pseudopotential in real space.
!!
!! SIDE EFFECTS
!!
!! NOTES
!! * Note that the present routine is tightly connected to the dfpt_vlocal.f routine,
!! that compute the derivative of the local ionic potential
!! with respect to one atomic displacement. The argument list
!! and the internal loops to be considered were sufficiently different
!! as to make the two routines different.
!! * The routine was adapted from mklocl.F90
!!
!! PARENTS
!!      dfpt_looppert,dfpt_nselt,dfpt_nstpaw
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine vlocalstr(gmet,gprimd,gsqcut,istr,mgfft,mpi_enreg,&
&  mqgrid,natom,nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,qgrid,&
&  ucvol,vlspl,vpsp1)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vlocalstr'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istr,mgfft,mqgrid,natom,nfft,ntypat,paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: vpsp1(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ig1,ig2,ig3,ii,itypat,jj
 integer :: ka,kb,n1,n2,n3
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,bb,cc,cutoff,dd,dgsquards,diff
 real(dp) :: dq,dq2div6,dqdiv6,dqm1,ee,ff,gmag,gsquar
 real(dp) :: sfi,sfr,term,vion1,vion2
 real(dp) :: xnorm
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgmetds(3,3)
 real(dp),allocatable :: work1(:,:)

! *************************************************************************

!Define G^2 based on G space metric gmet.
! gsq_vl(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
!& dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
!& dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

!Define dG^2/ds based on G space metric derivative dgmetds.
! dgsqds_vl(i1,i2,i3)=dble(i1*i1)*dgmetds(1,1)+dble(i2*i2)*dgmetds(2,2)+&
!& dble(i3*i3)*dgmetds(3,3)+&
!& dble(i1*i2)*(dgmetds(1,2)+dgmetds(2,1))+&
!& dble(i1*i3)*(dgmetds(1,3)+dgmetds(3,1))+&
!& dble(i2*i3)*(dgmetds(2,3)+dgmetds(3,2))

!Real and imaginary parts of phase--statment functions:
! phr_vl(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
! phi_vl(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
! ph1_vl(nri,i1,ia)=ph1d(nri,i1+1+n1+(ia-1)*(2*n1+1))
! ph2_vl(nri,i2,ia)=ph1d(nri,i2+1+n2+(ia-1)*(2*n2+1)+&
!& natom*(2*n1+1))
! ph3_vl(nri,i3,ia)=ph1d(nri,i3+1+n3+(ia-1)*(2*n3+1)+&
!& natom*(2*n1+1+2*n2+1))
! phre_vl(i1,i2,i3,ia)=phr_vl(ph1_vl(re,i1,ia),ph1_vl(im,i1,ia),ph2_vl(re,i2,ia),&
!& ph2_vl(im,i2,ia),ph3_vl(re,i3,ia),ph3_vl(im,i3,ia))
! phimag_vl(i1,i2,i3,ia)=phi_vl(ph1_vl(re,i1,ia),ph1_vl(im,i1,ia),ph2_vl(re,i2,ia),&
!& ph2_vl(im,i2,ia),ph3_vl(re,i3,ia),ph3_vl(im,i3,ia))

!-----
!Compute derivative of metric tensor wrt strain component istr
 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   ' Input istr=',istr,' not allowed.',ch10,&
&   ' Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Zero out array to permit accumulation over atom types below:
 ABI_ALLOCATE(work1,(2,nfft))
 work1(:,:)=0.0_dp
!
 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1

   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           ii=ii+1
!          ***     GET RID OF THIS THESE IF STATEMENTS (if they slow code)
!          Skip G=0:
!          if (ii==1) cycle
           if (ig1==0 .and. ig2==0 .and. ig3==0) cycle
           gsquar=gsq_vl(ig1,ig2,ig3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)
             dgsquards=dgsqds_vl(ig1,ig2,ig3)
!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al;
!            NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine).

             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6

             vion1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) +&
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) &
&             / gsquar

!            Also get (dV(q)/dq)/q:
!            (note correction of Numerical Recipes sign error
!            before (3._dp*aa**2-1._dp)
             ee= vlspl(jj+1,1,itypat)-vlspl(jj,1,itypat)
             ff=  (3._dp*bb**2-1._dp)*vlspl(jj+1,2,itypat) &
&             - (3._dp*aa**2-1._dp)*vlspl(jj,2,itypat)
             vion2 = ( ( ee*dqm1 + ff*dqdiv6 )/gmag&
&             - 2.0_dp*vion1                 ) / gsquar


!            Assemble structure factor over all atoms of given type:
             sfr=0.0_dp
             sfi=0.0_dp
             do ia=ia1,ia2
               sfr=sfr+phre_vl(ig1,ig2,ig3,ia)
               sfi=sfi-phimag_vl(ig1,ig2,ig3,ia)
             end do

             term=dgsquards*vion2
!            Add potential for diagonal strain components
             if(istr <=3) then
               term=term-vion1
             end if

!            Multiply structure factor times vion derivatives:
             work1(re,ii)=work1(re,ii)+sfr*term
             work1(im,ii)=work1(im,ii)+sfi*term

!            End skip G**2 outside cutoff:
           end if
!          End loop on n1, n2, n3. There is a "cycle" inside the loop
         end do
       end if
     end do
   end do

   ia1=ia2+1

!  End loop on type of atoms
 end do


!Set Vloc(G=0)=0:
 work1(re,1)=0.0_dp
 work1(im,1)=0.0_dp


!Transform back to real space
 call fourdp(1,work1,vpsp1,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!Divide by unit cell volume
 xnorm=1.0_dp/ucvol
 vpsp1(:)=vpsp1(:)*xnorm

 ABI_DEALLOCATE(work1)

 contains

!Real and imaginary parts of phase.
   function phr_vl(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phr_vl'
!End of the abilint section

   real(dp) :: phr_vl,x1,x2,x3,y1,y2,y3
   phr_vl=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl

   function phi_vl(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phi_vl'
!End of the abilint section

   real(dp):: phi_vl,x1,x2,x3,y1,y2,y3
   phi_vl=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl

   function ph1_vl(nri,ig1,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1_vl'
!End of the abilint section

   real(dp):: ph1_vl 
   integer :: nri,ig1,ia
   ph1_vl=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_vl

   function ph2_vl(nri,ig2,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph2_vl'
!End of the abilint section

   real(dp):: ph2_vl 
   integer :: nri,ig2,ia
   ph2_vl=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl

   function ph3_vl(nri,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph3_vl'
!End of the abilint section

   real(dp):: ph3_vl
   integer :: nri,ig3,ia
   ph3_vl=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl

   function phre_vl(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phre_vl'
!End of the abilint section

   real(dp):: phre_vl
   integer :: ig1,ig2,ig3,ia
   phre_vl=phr_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phre_vl

   function phimag_vl(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phimag_vl'
!End of the abilint section

   real(dp) :: phimag_vl
   integer :: ig1,ig2,ig3,ia
   phimag_vl=phi_vl(ph1_vl(re,ig1,ia),ph1_vl(im,ig1,ia),&
&   ph2_vl(re,ig2,ia),ph2_vl(im,ig2,ia),ph3_vl(re,ig3,ia),ph3_vl(im,ig3,ia))
 end function phimag_vl

   function gsq_vl(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_vl'
!End of the abilint section

   real(dp) :: gsq_vl
   integer :: i1,i2,i3
!Define G^2 based on G space metric gmet.
   gsq_vl=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_vl

   function dgsqds_vl(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dgsqds_vl'
!End of the abilint section

   real(dp) :: dgsqds_vl 
   integer :: i1,i2,i3
!Define dG^2/ds based on G space metric derivative dgmetds.
   dgsqds_vl=dble(i1*i1)*dgmetds(1,1)+dble(i2*i2)*dgmetds(2,2)+&
&   dble(i3*i3)*dgmetds(3,3)+&
&   dble(i1*i2)*(dgmetds(1,2)+dgmetds(2,1))+&
&   dble(i1*i3)*(dgmetds(1,3)+dgmetds(3,1))+&
&   dble(i2*i3)*(dgmetds(2,3)+dgmetds(3,2))
 end function dgsqds_vl

end subroutine vlocalstr
!!***
