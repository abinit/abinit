!{\src2tex{textfont=tt}}
!!****f* ABINIT/prcrskerker1
!! NAME
!! prcrskerker1
!!
!! FUNCTION
!! preconditionning by a real-space conjugate gradient on residual
!! using a model dielectric function in real space
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (PMA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  nfft=number of fft grid points
!!  nspden=number of spin-density components
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  gprimd(3,3)=dimensional primitive translations in fourier space (bohr**-1)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(nfft,nspden)=residual potential
!!  base(nfft) = real space function used as a basis to guess a fine dielectric funtion
!!  see the calling routine to know the content
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
!!  needs severe cleaning and this is abuse of modules as common blocks...
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      cgpr,frskerker1__end,frskerker1__init,laplacian,prc_mem_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prcrskerker1(dtset,mpi_enreg,nfft,nspden,ngfft,dielar,etotal,gprimd,vresid,vrespc,base)


 use defs_basis
 use defs_abitypes
 use frskerker1
 use mod_prc_memory
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prcrskerker1'
 use interfaces_56_recipspace
 use interfaces_62_cg_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden
 real(dp) :: etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: base(nfft),dielar(7),gprimd(3,3)
 real(dp),intent(in) :: vresid(nfft,nspden)
 real(dp),intent(out) :: vrespc(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer ::  ifft,ispden,n1,n2,n3
 real(dp) :: base_delta,base_max,base_min,dielng,diemac,diemix
 real(dp) :: diemixmag
 real(dp) :: rdummy1,rdummy2
 logical ::  new_prc_func
!arrays
 real(dp) :: deltaW(nfft,nspden)
 real(dp) :: g2cart(nfft)
 real(dp) :: mat(nfft,nspden)

! *************************************************************************

!DEBUG
!write(std_out,*)' prckerker1 : enter '
!ENDDEBUG
!if(cycle==0) then
 call prc_mem_init(nfft)

 if(cycle==0) then
   new_prc_func=.TRUE.
   energy_min=etotal
 else if(etotal < energy_min) then
   new_prc_func=.TRUE.
   energy_min=etotal
 else
   new_prc_func=.FALSE.
 end if


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
!base_cp=base
 if(new_prc_func) then
   base_min=base(1)
   base_max=base(1)
   do ifft=1,nfft
     base_min = min(base_min,base(ifft))
     base_max = max(base_max,base(ifft))
   end do
   base_delta = base_max - base_min
!  if(cycle.lt.2) then
   rdiemac(:) = (((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  else
!  rdiemac(:) = rdiemac(:)*0.5_dp+0.5_dp*(((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  end if
!  if(cycle==0) rdiemac(:) = (((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  rdiemac(:) = exp(((base(:)-base_min) / (base_delta) *log(diemac)))
 end if
 cycle=cycle+1
!if(cycle==5) cycle=0
!end if
!******************************************************************
!compute deltaW                                                 **
!******************************************************************
 vrespc=vresid !starting point
 call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,dtset%paral_kgb,rdfuncr=vrespc,laplacerdfuncr=deltaW,g2cart_out=g2cart) ! put the laplacian of the residuals into deltaW
!call laplacian(vrespc,buffer,ngfft,gprimd) ! put the laplacian of the residuals into deltaW
!do ifft=1,nfft
!if (buffer(ifft,1)/=deltaW(ifft,1)) then
!stop
!end if
!end do
 deltaW(:,1)= diemix*(((one/rdiemac(:))*vresid(:,1))-(((dielng)**2)*deltaW(:,1)))
 if (nspden>1.and.(diemixmag>=zero)) then
   do ispden=2,nspden
     deltaW(:,ispden)= abs(diemixmag)*(((one/rdiemac(:))*vresid(:,ispden))-(((dielng)**2)*deltaW(:,ispden)))
   end do
 end if
!call random_number(deltaW)
!call random_number(vrespc)
!******************************************************************
!Finding the preconditionned residuals which minimizes          **
!half*(vrespc*(1-dielng2/4pi2 nabla2) vrespc) - vrespc * deltaW **
!***********************************************************************
 vrespc(:,1)=diemix*vrespc(:,1)
 if (nspden>1) vrespc(:,2:nspden)=abs(diemixmag)*vrespc(:,2:nspden)
!buffer=vrespc


!==============================================================================
!==============================================================================
!! Original loop
!==============================================================================
!==============================================================================

 call frskerker1__init(dtset,mpi_enreg,nfft,ngfft,nspden,dielng,deltaW,gprimd,mat,g2cart)

!call cgpr(pf_rscgres,dpf_rscgres,newvres,real(1e-40,dp),700,vrespc,rdummy1,rdummy2)
!rdummy1 = pf_rscgres(nfft,nspden,vrespc)
 call cgpr(nfft,nspden,frskerker1__pf,frskerker1__dpf,frskerker1__newvres,&
& real(1e-10,dp),700,vrespc,rdummy1,rdummy2)
 call frskerker1__end()

!==============================================================================
!==============================================================================
!! Original loop end
!==============================================================================
!==============================================================================


!cplex=1
!qphon(:)=zero
!call moddiel(cplex,dielar,nfft,ngfft,nspden,1,0,qphon,rprimd,vresid,buffer)
!c1=0
!do ifft=1,nfft,1
!if((abs(buffer(ifft,1)-vrespc(ifft,1))/(abs(buffer(ifft,1)+vrespc(ifft,1))*half)) > 5e-3) then
!c1=c1+1
!end if
!end do
!call laplacian(vrespc,buffer,ngfft,gprimd)
!buffer=vrespc(:,:)-buffer(:,:)*dielng**2
!c2=0
!do ifft=1,nfft,1
!if((abs(buffer(ifft,1)-deltaW(ifft,1))/(abs(buffer(ifft,1)+deltaW(ifft,1))*half)) > 5e-3) then
!c2=c2+1
!end if
!end do
!!!  !stop
!call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,&
!& g2cart_out=g2cart)

!vrespc=vresid
!do ispden=1,nspden
!call fourdp(1, gvrespc(:,:,ispden), vrespc(:,ispden),-1,mpi_enreg,nfft,ngfft,0)
!end do
!filtering
!do ispden=1,nspden
!do ifft=1,nfft
!!    gvrespc(:,ifft,ispden)=(one-exp(-g2cart(ifft)*15.0_dp))*gvrespc(:,ifft,ispden)
!!      gvrespc(:,ifft,ispden)=(exp(-g2cart(ifft)*10.0_dp))*gvrespc(:,ifft,ispden)
!!      gvrespc(:,ifft,ispden)=(one-one/(exp(-0.002_dp/g2cart(ifft)**2)+one))*gvrespc(:,ifft,ispden)
!gvrespc(:,ifft,ispden)=(two-2_dp/(exp(-0.008_dp/(g2cart(ifft)+0.0012_dp))+one))*gvrespc(:,ifft,ispden)
!gvrespc(:,ifft,ispden)=min(one,(sqrt(g2cart(ifft)/0.006_dp))**(one))*gvrespc(:,ifft,ispden)
!end do
!end do
!change resulting potential to real space
!do ispden=1,nspden
!call fourdp(1,gvrespc(:,:,ispden),vrespc(:,ispden),1,mpi_enreg,nfft,ngfft,0)
!end do
!vrespc=vrespc*diemix
!maxg2=g2cart(1)
!ming2=g2cart(5)
!do ifft=1,nfft
!maxg2=max(g2cart(ifft),maxg2)
!if(g2cart(ifft) .gt. zero) ming2=min(g2cart(ifft),ming2)
!end do
!stop

!DEBUG
!write(std_out,*)' prckerker1 : exit '
!ENDDEBUG

end subroutine prcrskerker1
!!***
