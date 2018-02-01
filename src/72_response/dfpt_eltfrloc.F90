!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_eltfrloc
!! NAME
!! dfpt_eltfrloc
!!
!! FUNCTION
!! Compute the frozen-wavefunction local pseudopotential contribution
!! to the elastic tensor and the internal strain (derivative wrt one
!! cartesian strain component and one reduced-coordinate atomic displacement).
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!!  gsqcut=cutoff on G^2 based on ecut
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  rhog(2,nfft)=electron density in G space
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!
!! OUTPUT
!!  eltfrloc(6+3*natom,6)=non-symmetrized local pseudopotenial contribution
!!   to the elastic tensor and internal strain.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_eltfrloc(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfft,&
&  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,ph1d,qgrid,rhog,vlspl)


 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_mpinfo,  only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_eltfrloc'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),rhog(2,nfft),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: eltfrloc(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ielt,ieltx,ierr,ig1,ig2,ig3,ii
 integer :: is1,is2,itypat,jj,ka,kb,kd,kg,me_fft,n1,n2,n3,nproc_fft
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,bb,cc,cutoff,d2g 
 real(dp) :: dd,dg1,dg2,diff,dq
 real(dp) :: dq2div6,dqdiv6,dqm1,ee,ff,gmag,gsquar
 real(dp) :: sfi,sfr,term,term1
!real(dp) :: ph1_elt,ph2_elt,ph3_elt,phi_elt,phr_elt
 real(dp) :: term2,term3,term4,term5,vion1,vion2,vion3
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgm(3,3,6),tsec(2)
 real(dp),allocatable :: d2gm(:,:,:,:),elt_work(:,:)

! *************************************************************************

!Define G^2 based on G space metric gmet.
! gsq_elt(i1,i2,i3)=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
!& dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
!& dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)

!Define dG^2/ds based on G space metric derivative 
! dgsqds_elt(i1,i2,i3,is)=dble(i1*i1)*dgm(1,1,is)+dble(i2*i2)*dgm(2,2,is)+&
!& dble(i3*i3)*dgm(3,3,is)+&
!& dble(i1*i2)*(dgm(1,2,is)+dgm(2,1,is))+&
!& dble(i1*i3)*(dgm(1,3,is)+dgm(3,1,is))+&
!& dble(i2*i3)*(dgm(2,3,is)+dgm(3,2,is))

!Define 2dG^2/ds1ds2  based on G space metric derivative 
! d2gsqds_elt(i1,i2,i3,is1,is2)=dble(i1*i1)*d2gm(1,1,is1,is2)+&
!& dble(i2*i2)*d2gm(2,2,is1,is2)+dble(i3*i3)*d2gm(3,3,is1,is2)+&
!& dble(i1*i2)*(d2gm(1,2,is1,is2)+d2gm(2,1,is1,is2))+&
!& dble(i1*i3)*(d2gm(1,3,is1,is2)+d2gm(3,1,is1,is2))+&
!& dble(i2*i3)*(d2gm(2,3,is1,is2)+d2gm(3,2,is1,is2))

!Real and imaginary parts of phase--statment functions:
! phr_elt(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
! phi_elt(x1,y1,x2,y2,x3,y3)=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
! ph1_elt(nri,i1,ia)=ph1d(nri,i1+1+n1+(ia-1)*(2*n1+1))
! ph2_elt(nri,i2,ia)=ph1d(nri,i2+1+n2+(ia-1)*(2*n2+1)+&
!& natom*(2*n1+1))
! ph3_elt(nri,i3,ia)=ph1d(nri,i3+1+n3+(ia-1)*(2*n3+1)+&
!& natom*(2*n1+1+2*n2+1))
! phre_elt(i1,i2,i3,ia)=phr_elt(ph1_elt(re,i1,ia),ph1_elt(im,i1,ia),ph2_elt(re,i2,ia),&
!& ph2_elt(im,i2,ia),ph3_elt(re,i3,ia),ph3_elt(im,i3,ia))
! phimag_elt(i1,i2,i3,ia)=phi_elt(ph1_elt(re,i1,ia),ph1_elt(im,i1,ia),ph2_elt(re,i2,ia),&
!& ph2_elt(im,i2,ia),ph3_elt(re,i3,ia),ph3_elt(im,i3,ia))

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11) ; nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!-----
!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.
 ABI_ALLOCATE(d2gm,(3,3,6,6))

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
   end do
!  Loop over 1st strain index, upper triangle only
   do is1=1,is2
     ka=idx(2*is1-1);kb=idx(2*is1)
     d2gm(:,:,is1,is2)=0._dp
     do jj = 1,3
       if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
       if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
       if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
       if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
     end do
   end do !is1
 end do !is2

!Zero out array to permit accumulation over atom types below:
 eltfrloc(:,:)=0.0_dp

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
       if (fftn2_distrib(i2)==me_fft) then
         ig2=i2-(i2/id2)*n2-1
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1

           ii=ii+1
!          Skip G=0:
           if (ig1==0 .and. ig2==0 .and. ig3==0) cycle

!          Skip G**2 outside cutoff:
           gsquar=gsq_elt(ig1,ig2,ig3)
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al; NOTE error in book for sign
!             of "aa" term in derivative; also see splfit routine).
             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             term1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) +&
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat))
             vion1=term1 / gsquar

!            Also get dV(q)/dq:
!            (note correction of Numerical Recipes sign error
!             before (3._dp*aa**2-1._dp)
             ee= vlspl(jj+1,1,itypat)-vlspl(jj,1,itypat)
             ff=  (3._dp*bb**2-1._dp)*vlspl(jj+1,2,itypat) &
&             - (3._dp*aa**2-1._dp)*vlspl(jj,2,itypat)
             term2 = ee*dqm1 + ff*dqdiv6
             vion2 = term2/gsquar - 2._dp*term1/(gsquar*gmag)

!            Also get V''(q)
             term3=aa*vlspl(jj,2,itypat)+bb*vlspl(jj+1,2,itypat)
             vion3 = (term3 - 4.0_dp*term2/gmag + 6._dp*term1/gsquar)/gsquar

!            Assemble structure factor over all atoms of given type:
             sfr=zero;sfi=zero
             do ia=ia1,ia2
               sfr=sfr+phre_elt(ig1,ig2,ig3,ia)
               sfi=sfi-phimag_elt(ig1,ig2,ig3,ia)
             end do
             term=(rhog(re,ii)*sfr+rhog(im,ii)*sfi)

!            Loop over 2nd strain index
             do is2=1,6
               dg2=0.5_dp*dgsqds_elt(ig1,ig2,ig3,is2)/gmag
!              Loop over 1st strain index, upper triangle only
               do is1=1,is2
                 dg1=0.5_dp*dgsqds_elt(ig1,ig2,ig3,is1)/gmag
                 d2g=(0.25_dp*d2gsqds_elt(ig1,ig2,ig3,is1,is2)-dg1*dg2)/gmag
                 eltfrloc(is1,is2)=eltfrloc(is1,is2)+&
&                 term*(vion3*dg1*dg2+vion2*d2g) 
                 if(is2<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)-term*vion2*dg1
                 if(is1<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)-term*vion2*dg2
                 if(is1<=3 .and. is2<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)+term*vion1
               end do !is1

!              Internal strain section - loop over current atoms
               do ia=ia1,ia2
                 if(is2 <=3) then
                   term4=vion2*dg2-vion1
                 else
                   term4=vion2*dg2
                 end if
                 term5=-two_pi*(rhog(re,ii)*phimag_elt(ig1,ig2,ig3,ia)&
&                 +rhog(im,ii)*phre_elt(ig1,ig2,ig3,ia))*term4
                 eltfrloc(7+3*(ia-1),is2)=eltfrloc(7+3*(ia-1),is2)+term5*dble(ig1)
                 eltfrloc(8+3*(ia-1),is2)=eltfrloc(8+3*(ia-1),is2)+term5*dble(ig2)
                 eltfrloc(9+3*(ia-1),is2)=eltfrloc(9+3*(ia-1),is2)+term5*dble(ig3)
               end do

             end do !is2

!            End skip G**2 outside cutoff:
           end if

!          End loop on n1, n2, n3. There is a "cycle" inside the loop
         end do
       end if
     end do
   end do

!  End loop on type of atoms
   ia1=ia2+1
 end do
!Init mpi_comm
 call timab(48,1,tsec)
 call xmpi_sum(eltfrloc,mpi_enreg%comm_fft,ierr)
 call timab(48,2,tsec)

!Fill in lower triangle
 do is2=2,6
   do is1=1,is2-1
     eltfrloc(is2,is1)=eltfrloc(is1,is2)
   end do
 end do

!The indexing array atindx is used to reestablish the correct
!order of atoms
 ABI_ALLOCATE(elt_work,(6+3*natom,6))
 elt_work(1:6,1:6)=eltfrloc(1:6,1:6)
 do ia=1,natom
   ielt=7+3*(ia-1)
   ieltx=7+3*(atindx(ia)-1)
   elt_work(ielt:ielt+2,1:6)=eltfrloc(ieltx:ieltx+2,1:6)
 end do
 eltfrloc(:,:)=elt_work(:,:)

 ABI_DEALLOCATE(d2gm)
 ABI_DEALLOCATE(elt_work)

 contains

!Real and imaginary parts of phase.
   function phr_elt(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phr_elt'
!End of the abilint section

   real(dp) :: phr_elt,x1,x2,x3,y1,y2,y3
   phr_elt=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_elt

   function phi_elt(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phi_elt'
!End of the abilint section

   real(dp):: phi_elt,x1,x2,x3,y1,y2,y3
   phi_elt=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_elt

   function ph1_elt(nri,ig1,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1_elt'
!End of the abilint section

   real(dp):: ph1_elt
   integer :: nri,ig1,ia
   ph1_elt=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_elt

   function ph2_elt(nri,ig2,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph2_elt'
!End of the abilint section

   real(dp):: ph2_elt
   integer :: nri,ig2,ia
   ph2_elt=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_elt

   function ph3_elt(nri,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph3_elt'
!End of the abilint section

   real(dp):: ph3_elt
   integer :: nri,ig3,ia
   ph3_elt=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_elt

   function phre_elt(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phre_elt'
!End of the abilint section

   real(dp):: phre_elt
   integer :: ig1,ig2,ig3,ia
   phre_elt=phr_elt(ph1_elt(re,ig1,ia),ph1_elt(im,ig1,ia),&
&   ph2_elt(re,ig2,ia),ph2_elt(im,ig2,ia),ph3_elt(re,ig3,ia),ph3_elt(im,ig3,ia))
 end function phre_elt

   function phimag_elt(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phimag_elt'
!End of the abilint section

   real(dp) :: phimag_elt
   integer :: ig1,ig2,ig3,ia
   phimag_elt=phi_elt(ph1_elt(re,ig1,ia),ph1_elt(im,ig1,ia),&
&   ph2_elt(re,ig2,ia),ph2_elt(im,ig2,ia),ph3_elt(re,ig3,ia),ph3_elt(im,ig3,ia))
 end function phimag_elt

   function gsq_elt(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_elt'
!End of the abilint section

   real(dp) :: gsq_elt
   integer :: i1,i2,i3
!Define G^2 based on G space metric gmet.
   gsq_elt=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_elt
 
   function dgsqds_elt(i1,i2,i3,is)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dgsqds_elt'
!End of the abilint section

   real(dp) :: dgsqds_elt
   integer :: i1,i2,i3,is
!Define dG^2/ds based on G space metric derivative 
   dgsqds_elt=dble(i1*i1)*dgm(1,1,is)+dble(i2*i2)*dgm(2,2,is)+&
&   dble(i3*i3)*dgm(3,3,is)+&
&   dble(i1*i2)*(dgm(1,2,is)+dgm(2,1,is))+&
&   dble(i1*i3)*(dgm(1,3,is)+dgm(3,1,is))+&
&   dble(i2*i3)*(dgm(2,3,is)+dgm(3,2,is))
 end function dgsqds_elt

   function d2gsqds_elt(i1,i2,i3,is1,is2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2gsqds_elt'
!End of the abilint section

   real(dp) :: d2gsqds_elt
   integer :: i1,i2,i3,is1,is2
!Define 2dG^2/ds1ds2  based on G space metric derivative 
   d2gsqds_elt=dble(i1*i1)*d2gm(1,1,is1,is2)+&
&   dble(i2*i2)*d2gm(2,2,is1,is2)+dble(i3*i3)*d2gm(3,3,is1,is2)+&
&   dble(i1*i2)*(d2gm(1,2,is1,is2)+d2gm(2,1,is1,is2))+&
&   dble(i1*i3)*(d2gm(1,3,is1,is2)+d2gm(3,1,is1,is2))+&
&   dble(i2*i3)*(d2gm(2,3,is1,is2)+d2gm(3,2,is1,is2))
 end function d2gsqds_elt

end subroutine dfpt_eltfrloc
!!***
