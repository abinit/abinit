!{\src2tex{textfont=tt}}
!!****f* ABINIT/prcrskerker2
!! NAME
!! prcrskerker2
!!
!! FUNCTION
!! preconditionning by a real-space conjugate gradient on residual
!! using a model dielectric function in real space
!! differing from prcrskerker1 by the
!! use of a linear response approach
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (PMA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  nfft=number of fft grid points
!!  nspden=number of spin-density components
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  gprimd(3,3)=dimensional primitive translations in fourier space (bohr**-1)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(nfft,nspden)=residual potential
!!
!! OUTPUT
!!  vrespc(nfft,nspden)=preconditioned residual of the potential
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      cgpr,dotprod_vn,frskerker2__end,frskerker2__init,laplacian,ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prcrskerker2(dtset,nfft,nspden,ngfft,dielar,gprimd,rprimd,vresid,vrespc,natom,xred,mpi_enreg,ucvol)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use frskerker2

 use m_mpinfo,  only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prcrskerker2'
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_62_cg_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dielar(7),gprimd(3,3),rprimd(3,3),vresid(nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vrespc(nfft,nspden)

!Local variables-------------------------------
  !logical,save ::ok=.FALSE.
!scalars
 integer :: cplex,i1,i2,i3,iatom,iatom27,ifft,ispden,n1,n2,n3,natom27,nfftotf
 integer :: option
 real(dp),save :: lastp1=one,lastp2=one
 real(dp) :: C1,C2,DE,core,dielng,diemac,diemix,diemixmag,doti,dr,l1,l2,l3,l4,r
 real(dp) :: rdummy1,rdummy2,rmin,xr,y,yr,zr
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: V1(nfft,nspden),V2(nfft,nspden),buffer(nfft,nspden)
 real(dp) :: deltaW(nfft,nspden)
 real(dp) :: mat(nfft,nspden)
 real(dp) :: rdielng(nfft),rdiemac(nfft),xcart(3,natom)
 real(dp) :: xcart27(3,natom*27)

! *************************************************************************

 dielng=dielar(2)
 diemac=dielar(3)
 diemix=dielar(4)
 diemixmag=dielar(7)
!******************************************************************
!compute the diemac(r)                                          **
!******************************************************************
!this task will be devoted to a general function later
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 nfftotf=n1*n2*n3
!if(.not.ok) then
 xcart(1,:)=xred(1,:)*rprimd(1,1)+xred(2,:)*rprimd(1,2)+xred(3,:)*rprimd(1,3)
 xcart(2,:)=xred(1,:)*rprimd(2,1)+xred(2,:)*rprimd(2,2)+xred(3,:)*rprimd(2,3)
 xcart(3,:)=xred(1,:)*rprimd(3,1)+xred(2,:)*rprimd(3,2)+xred(3,:)*rprimd(3,3)

 iatom27=0
 do i1=-1,1
   do i2=-1,1
     do i3=-1,1
       do iatom=1,natom
         iatom27=iatom27+1
         xcart27(:,iatom27)=xcart(:,iatom)+rprimd(:,1)*i1+rprimd(:,2)*i2+rprimd(:,3)*i3
       end do
     end do
   end do
 end do

!stop
 natom27=27*natom

 l1=0.34580850339844665
!l2=0.5123510203906797 !0.41242551019533985
!l3=0.8001489796093203 !0.90007448980466009

 l2=0.41242551019533985
 l3=0.90007448980466009
 l4=0.9666914966015534


 l1=0.31387233559896449
 l2=0.35828367346355994
 l3=0.9333829932031068
 l4=0.9777943310677023

 l1=3.5
 l2=11.5
 l3=2.5
 l4=6.5
!l1=30. !cellules pleines

 rdielng=zero
 core=1. !(value of Er at the core of atoms)
 dr=2.65 ! radius of atoms=2.65165
 y=1. ! value of Er in the empty region

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do i3=1,n3
   ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
   do i2=1,n2
     if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
       do i1=1,n1
         ifft=ifft+1
!        !!!!!!!!!!!!!!!!!!!!!!!!!
!        ! calculation of the simplest part void/metal
!        !!              x=real(real(i3,dp)/real(n3,dp),dp)
!        !!              !x=i3/n3
!        !!              if(x < l1) then
!        !!                 rdiemac(ifft)=diemac
!        !!                 rdielng(ifft)=dielng
!        !!              else if(x < l2) then
!        !!                 xp=(l2-x)/(l2-l1)
!        !!                 rdiemac(ifft)=y+(diemac-y)&
!        !!                      & * (1.-(1.-xp)**4)**4
!        !!                 rdielng(ifft)=dielng*(1.-(1.-xp)**4)**4
!        !!              else if(x < l3) then
!        !!                 rdiemac(ifft)=y
!        !!                 rdielng(ifft)=zero
!        !!              else if(x < l4) then
!        !!                 xp=(l3-x)/(l3-l4)
!        !!                 rdiemac(ifft)=y+(diemac-y)&
!        !!                      & * (1.-(1.-xp)**4)**4
!        !!                 rdielng(ifft)=dielng*(1.-(1.-xp)**4)**4
!        !!              else
!        !!                 rdiemac(ifft)=diemac
!        !!                 rdielng(ifft)=dielng
!        !!              end if
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!!! calculation of atomic core dielectric
!        !!              rmin=1e16
!        !!              xr=real(real((i1-1),dp)/n1,dp)*rprimd(1,1)+real(real((i2-1),dp)/n2,dp)*rprimd(1,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(1,3)
!        !!              yr=real(real((i1-1),dp)/n1,dp)*rprimd(2,1)+real(real((i2-1),dp)/n2,dp)*rprimd(2,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(2,3)
!        !!              zr=real(real((i1-1),dp)/n1,dp)*rprimd(3,1)+real(real((i2-1),dp)/n2,dp)*rprimd(3,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(3,3)
!        !!              do iatom=1,natom27
!        !!                 r=(xr-xcart27(1,iatom))**2+(yr-xcart27(2,iatom))**2+(zr-xcart27(3,iatom))**2
!        !!                 if (r<rmin) then
!        !!                    rmin=r
!        !!                 end if
!        !!              end do
!        !!              if(rmin < dr**2) then
!        !!                 rdiemac(ifft)=min(rdiemac(ifft),core+(diemac-core)*(1.-(1.-sqrt(rmin)/dr)**2)**2)
!        !!                 rdielng(ifft)=dielng-dielng*(1.-(1.-sqrt(rmin)/dr)**4)**4
!        !!              else
!        !!                 rdiemac(ifft)=min(rdiemac(ifft),diemac)
!        !!              end if
         rmin=1e16
         xr=real(real((i1-1),dp)/n1,dp)*rprimd(1,1)+real(real((i2-1),dp)/n2,dp)*rprimd(1,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(1,3)
         yr=real(real((i1-1),dp)/n1,dp)*rprimd(2,1)+real(real((i2-1),dp)/n2,dp)*rprimd(2,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(2,3)
         zr=real(real((i1-1),dp)/n1,dp)*rprimd(3,1)+real(real((i2-1),dp)/n2,dp)*rprimd(3,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(3,3)

         rdiemac(ifft)=y
         rdielng(ifft)=zero
         do iatom=1,natom27
           r=(xr-xcart27(1,iatom))**2+(yr-xcart27(2,iatom))**2+(zr-xcart27(3,iatom))**2
           if (r<rmin) then
             rmin=r
           end if
           if(r < l1) then
             rdiemac(ifft)= rdiemac(ifft) +  0.7_dp * (diemac-y)
           else if(r < l2) then
             rdiemac(ifft)= rdiemac(ifft) + 0.7_dp * (diemac-y)*(one-((sqrt(r)-l1)/(l2-l1))**2)**2
           else
             rdiemac(ifft)=rdiemac(ifft)
           end if
           if(r < l3) then
             rdielng(ifft)= rdielng(ifft) +  0.5_dp * (dielng)
           else if(r < l4) then
             rdielng(ifft)= rdielng(ifft) + 0.5_dp * (dielng)  *(one-((sqrt(r)-l3)/(l4-l3))**2)**2
           end if
         end do

         rdielng(ifft)=min(rdielng(ifft),dielng)
!        rdielng(ifft)=dielng
         rdiemac(ifft)=min(rdiemac(ifft),diemac)
         rdiemac(ifft)=diemac
       end do
     end if
   end do
 end do
!rdielng(:)=dielng

!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!******************************************************************
!compute V1
!******************************************************************
 V1=vresid
 call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,dtset%paral_kgb,rdfuncr=V1,laplacerdfuncr=deltaW)
 deltaW(:,1)= (((one/rdiemac(:))*V1(:,1))-(((rdielng(:))**2)*deltaW(:,1)))
!deltaW(:,1)= -diemix*(((rdielng(:))**2)*deltaW(:,ispden))
 if (nspden>1) then
   do ispden=2,nspden
     deltaW(:,ispden)= (((one/rdiemac(:))*V1(:,ispden))-(((rdielng(:))**2)*deltaW(:,ispden)))
!    deltaW(:,ispden)= -abs(diemixmag)*(((rdielng(:))**2)*deltaW(:,ispden))
   end do
 end if
 call frskerker2__init(dtset,mpi_enreg,nfft,ngfft,nspden,rdielng,deltaW,gprimd,mat)
 call cgpr(nfft,nspden,frskerker2__pf,frskerker2__dpf,&
& frskerker2__newvres2,lastp1*real(1e-6 ,dp),700,V1,rdummy1,rdummy2)
 lastp1=min(abs(rdummy1),1e-6_dp)
 call frskerker2__end()

!******************************************************************
!compute V2
!******************************************************************
 V2=vresid
 do ispden=1,nspden
   deltaW(:,ispden)= (rdielng(:)**2)
 end do
 call frskerker2__init(dtset,mpi_enreg,nfft,ngfft,nspden,rdielng,deltaW,gprimd,mat)
 call cgpr(nfft,nspden,frskerker2__pf,frskerker2__dpf,&
& frskerker2__newvres2,lastp2*real(1e-6,dp),700,V2,rdummy1,rdummy2)
 lastp2=min(abs(rdummy1),1e-6_dp)
 call frskerker2__end()


!******************************************************************
!compute C1, C2 & DE
!******************************************************************
 cplex=1;
 option=1;
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&DE,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&rdielng,&          !the potential
&ucvol,&          !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 do ispden=1,nspden
   buffer(:,ispden)=rdielng(:)*V1(:,ispden)
 end do
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&C1,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&buffer,&          !the potential
&ucvol,&         !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 do ispden=1,nspden
   buffer(:,ispden)=rdielng(:)*V2(:,ispden)
 end do
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&C2,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&buffer,&          !the potential
&ucvol,&         !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 C1=C1/DE
 C2=C2/DE
 DE=C1/(one-C2)

!******************************************************************
!compute the new preconditionned residuals
!******************************************************************
 vrespc(:,1)=diemix*(V1(:,1)+DE*V2(:,1))
 if (nspden>1) vrespc(:,2:nspden)=abs(diemixmag)*(V1(:,2:nspden)+DE*V2(:,2:nspden))

end subroutine prcrskerker2
!!***
