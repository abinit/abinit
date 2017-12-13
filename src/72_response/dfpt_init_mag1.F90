!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_init_mag1
!! NAME
!!  dfpt_init_mag1
!!
!! FUNCTION
!!  Initial guess of the first order magnetization/density for magnetic field perturbation.
!!  The first order magnetization is set so as to zero out the first order XC magnetic field, which
!!  should minimize the second order XC energy (without taking self-consistency into account).
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ipert = perturbtation type (works only for ipert==natom+5)
!!  idir  = direction of the applied magnetic field
!!  cplex = complex or real first order density and magnetization
!!  nfft  = dimension of the fft grid
!!  nspden= number of density matrix components
!!  nkxc  = number of kxc components
!!  vxc0(nfft,nspden)  = GS XC potential 
!!  kxc0(nfft,nspden)  = GS XC derivatives 
!!  rhor0(nfft,nspden) = GS density matrix 
!!
!! OUTPUT
!!  rhor1(cplex*nfft) = first order density magnetization guess 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!!  dfpt_looppert.F90
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_init_mag1(ipert,idir,rhor1,rhor0,cplex,nfft,nspden,vxc0,kxc0,nkxc)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_init_mag1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)    :: ipert,idir,cplex,nfft,nspden,nkxc
 real(dp), intent(in)    :: vxc0(nfft,nspden),rhor0(nfft,nspden)
 real(dp), intent(in)    :: kxc0(nfft,nkxc)
 real(dp), intent(out)   :: rhor1(cplex*nfft,nspden)                        

!Local variables-------------------------------
 integer  :: ipt                                     
 real(dp) :: bxc0,phixc0,m_dot_m1,bxc1                  
 real(dp) :: mdir(3),fdir(3)               
 real(dp) :: m1(3),m0(3),m1_norm,m0_norm
 real(dp) :: f_dot_m,f_perp(3)
 
! *************************************************************************

 if (nspden==2) then
   if(cplex==1) then
     do ipt=1,nfft
       bxc1=half*(half*(kxc0(ipt,1)+kxc0(ipt,3))-kxc0(ipt,2)) ! d/dm Bxc
       !this overestimates the first order magnetization because of n1 not taken into account
       m1_norm=-half*(1/bxc1)
       rhor1(ipt,1)=zero             ! rho_up+rho_dwn    => charge density
       rhor1(ipt,2)=half*m1_norm     ! rho_up=1/2(rho+m) => half*m
     enddo
   else
     do ipt=1,cplex*nfft
       rhor1(ipt,:)=zero
     enddo
   endif
 else if(nspden==4) then
   if(cplex==1) then

     fdir=zero
     fdir(idir)= 1.0d0
     do ipt=1,nfft  
       m0_norm=sqrt(rhor0(ipt,2)**2+rhor0(ipt,3)**2+rhor0(ipt,4)**2)
       mdir(1)=rhor0(ipt,2)/m0_norm
       mdir(2)=rhor0(ipt,3)/m0_norm
       mdir(3)=rhor0(ipt,4)/m0_norm
       f_dot_m=fdir(1)*mdir(1)+fdir(2)*mdir(2)+fdir(3)*mdir(3) ! projection of the field direction on m0

       bxc1=half*(half*(kxc0(ipt,1)+kxc0(ipt,3))-kxc0(ipt,2))  ! d/dm Bxc
       m1_norm=(-half/bxc1)*f_dot_m                            ! get an estimate of the norm of m1

       bxc0=-sqrt((half*(vxc0(ipt,1)-vxc0(ipt,2)))**2+vxc0(ipt,3)**2+vxc0(ipt,4)**2)       

       rhor1(ipt,1)=zero       ! rho_up+rho_dwn    => charge density
       rhor1(ipt,2)=m1_norm*mdir(1)-half*m0_norm/bxc0*(fdir(1)-f_dot_m*mdir(1))   ! m1x
       rhor1(ipt,3)=m1_norm*mdir(2)-half*m0_norm/bxc0*(fdir(2)-f_dot_m*mdir(2))   ! m1x
       rhor1(ipt,4)=m1_norm*mdir(3)-half*m0_norm/bxc0*(fdir(3)-f_dot_m*mdir(3))   ! m1x

       rhor1(ipt,:)=zero
       !write(*,*) ipt,rhor0(ipt,2),rhor0(ipt,3),rhor0(ipt,4),bxc1,bxc0,half*(vxc0(ipt,1)-vxc0(ipt,2))/mdir(3)
       !write(*,*) ipt,mdir(3),m1_norm,rhor1(ipt,2),rhor1(ipt,3),rhor1(ipt,4)
       !write(*,*) ipt,bxc1,bxc0
     enddo
   else
     do ipt=1,cplex*nfft
       rhor1(ipt,:)=zero
     enddo
   endif
 endif


end subroutine dfpt_init_mag1
!!***
