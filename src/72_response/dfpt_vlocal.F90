!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_vlocal
!! NAME
!! dfpt_vlocal
!!
!! FUNCTION
!! Compute local part of 1st-order potential from the appropriate
!! atomic pseudopotential with structure and derivative factor.
!! In case of derivative with respect to k or
!! electric field perturbation, the 1st-order local potential vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG,MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cplex: if 1, real space 1-order functions on FFT grid
!!    are REAL, if 2, COMPLEX
!!  gmet(3,3)=reciprocal space metric (Bohr**-2)
!!  gsqcut=cutoff G**2 for included G s in fft box.
!!  idir=direction of atomic displacement (=1,2 or 3 : displacement of
!!    atom ipert along the 1st, 2nd or 3rd axis).
!!  ipert=number of the atom being displaced in the frozen-phonon
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimension of q grid for pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information.
!!  qgrid(mqgrid)=grid of q points from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  ucvol=unit cell volume (Bohr**3).
!!  vlspl(mqgrid,2,ntypat)=spline fit of q^2 V(q) for each type of atom.
!!  xred(3,natom)=reduced atomic coordinates
!!
!! OUTPUT
!!  vpsp1(cplex*nfft)=first-order local crystal pseudopotential in real space
!!    (including the minus sign, forgotten in the paper non-linear..
!!
!! PARENTS
!!      dfpt_looppert,dfpt_nstdy,dfpt_nstpaw,dfptnl_loop
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir,ipert,&
& mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
& ntypat,n1,n2,n3,paral_kgb,ph1d,qgrid,qphon,ucvol,vlspl,vpsp1,xred)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_vlocal'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mqgrid,n1,n2,n3,natom,nfft,ntypat
 integer,intent(in) :: paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),vlspl(mqgrid,2,ntypat)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vpsp1(cplex*nfft)

!Local variables -------------------------
!scalars
 integer :: i1,i2,i3,ia1,iatom,id1,id2,id3,ig1,ig2,ig3,ii,ii1,im=2
 integer :: itypat,jj,re=1
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: aa,bb,cc,cutoff,dd,diff,dq,dq2div6,dqdiv6,dqm1,gmag,gq1
 real(dp) :: gq2,gq3,gsquar,phqim,phqre
 real(dp) :: qxred2pi,sfi,sfr,vion1,xnorm
 logical :: qeq0
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gq(3)
 real(dp),allocatable :: work1(:,:)

! *********************************************************************

 iatom=ipert

 if(iatom==natom+1 .or. iatom==natom+2 .or. iatom==natom+10  .or. iatom==natom+11)then

!  (In case of d/dk or an electric field)
   vpsp1(1:cplex*nfft)=zero

 else

!  (In case of a phonon perturbation)
   ABI_ALLOCATE(work1,(2,nfft))
   work1(1:2,1:nfft)=0.0_dp

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=1.0_dp/dq
   dqdiv6=dq/6.0_dp
   dq2div6=dq**2/6.0_dp
   cutoff=gsqcut*tolfix
   id1=n1/2+2
   id2=n2/2+2
   id3=n3/2+2

   ! Get the distrib associated with this fft_grid
   call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!  This is to allow q=0
   qeq0=.false.
   if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)qeq0=.true.

!  Determination of the atom type
   ia1=0
   itypat=0
   do ii=1,ntypat
     ia1=ia1+nattyp(ii)
     if(atindx(iatom)<=ia1.and.itypat==0)itypat=ii
   end do

!  Determination of phase qxred*
   qxred2pi=2.0_dp*pi*(qphon(1)*xred(1,iatom)+ &
&   qphon(2)*xred(2,iatom)+ &
&   qphon(3)*xred(3,iatom) )
   phqre=cos(qxred2pi)
   phqim=sin(qxred2pi)
   ii=0

   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gq3=dble(ig3)+qphon(3)
     gq(3)=gq3
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-(i2/id2)*n2-1
         gq2=dble(ig2)+qphon(2)
         gq(2)=gq2

!        Note the lower limit of the next loop
         ii1=1
         if(i3==1 .and. i2==1 .and. qeq0 .and. ig2==0 .and. ig3==0)then
           ii1=2
           ii=ii+1
         end if
         do i1=ii1,n1
           ig1=i1-(i1/id1)*n1-1
           gq1=dble(ig1)+qphon(1)
           gq(1)=gq1
           ii=ii+1
           gsquar=gsq_vl3(gq1,gq2,gq3)
!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al; NOTE error in book for sign
!            of "aa" term in derivative; also see splfit routine.
!            This bug fixed here 27 Jan 1992.)

             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             vion1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) + &
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat) ) &
&             / gsquar

!            Phase   G*xred  (complex conjugate) * -i *2pi*(g+q)*vion
             sfr=-phimag_vl3(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1
             sfi=-phre_vl3(ig1,ig2,ig3,iatom)*2.0_dp*pi*gq(idir)*vion1

!            Phase   q*xred  (complex conjugate)
             work1(re,ii)=sfr*phqre+sfi*phqim
             work1(im,ii)=-sfr*phqim+sfi*phqre
           end if

         end do
       end if
     end do
   end do

!  Transform back to real space
   call fourdp(cplex,work1,vpsp1,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

   xnorm=1.0_dp/ucvol
   vpsp1(1:cplex*nfft)=vpsp1(1:cplex*nfft)*xnorm

   ABI_DEALLOCATE(work1)

!  End the condition of non-electric-field
 end if

 contains

!Real and imaginary parts of phase.
   function phr_vl3(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phr_vl3'
!End of the abilint section

   real(dp) :: phr_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phr_vl3=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_vl3

   function phi_vl3(x1,y1,x2,y2,x3,y3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phi_vl3'
!End of the abilint section

   real(dp) :: phi_vl3
   real(dp),intent(in) :: x1,x2,x3,y1,y2,y3
   phi_vl3=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph1_vl3(nri,ig1,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph1_vl3'
!End of the abilint section

   real(dp) :: ph1_vl3
   integer,intent(in) :: nri,ig1,ia
   ph1_vl3=ph1d(nri,ig1+1+n1+(atindx(ia)-1)*(2*n1+1))
 end function ph1_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph2_vl3(nri,ig2,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph2_vl3'
!End of the abilint section

   real(dp) :: ph2_vl3
   integer,intent(in) :: nri,ig2,ia
   ph2_vl3=ph1d(nri,ig2+1+n2+(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_vl3

!  Warning : this function differ from similar ones for ground-state calculations : note the atindx !!
   function ph3_vl3(nri,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ph3_vl3'
!End of the abilint section

   real(dp) :: ph3_vl3
   integer,intent(in) :: nri,ig3,ia
   ph3_vl3=ph1d(nri,ig3+1+n3+(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_vl3

   function phre_vl3(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phre_vl3'
!End of the abilint section

   real(dp) :: phre_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phre_vl3=phr_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phre_vl3

   function phimag_vl3(ig1,ig2,ig3,ia)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'phimag_vl3'
!End of the abilint section

   real(dp) :: phimag_vl3
   integer,intent(in) :: ig1,ig2,ig3,ia
   phimag_vl3=phi_vl3(ph1_vl3(re,ig1,ia),ph1_vl3(im,ig1,ia),&
&   ph2_vl3(re,ig2,ia),ph2_vl3(im,ig2,ia),ph3_vl3(re,ig3,ia),ph3_vl3(im,ig3,ia))
 end function phimag_vl3

   function gsq_vl3(g1,g2,g3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_vl3'
!End of the abilint section

   real(dp) :: gsq_vl3
   real(dp),intent(in) :: g1,g2,g3 ! Note that they are real, unlike in other similar function definitions
!Define G^2 based on G space metric gmet.
   gsq_vl3=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+&
&   g3*g3*gmet(3,3)+2.0_dp*g1*g2*gmet(1,2)+&
&   2.0_dp*g2*g3*gmet(2,3)+2.0_dp*g3*g1*gmet(3,1)
 end function gsq_vl3

end subroutine dfpt_vlocal
!!***
