!{\src2tex{textfont=tt}}
!!****f* ABINIT/density_rec
!! NAME
!! density_rec
!! 
!! FUNCTION
!! This routine computes the density using  the coefficients corresponding to
!! continued fraction at a point from a fixed potential. 
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (SLeroux,MMancini).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  coordx, coordy, coordz=coordonnees of the computed point
!!  an, bn2 : coefficient given by density_rec. Input if get_rec_coef=0, output else
!!  nrec=order of density_rec
!!  fermie=fermi energy (Hartree)
!!  tsmear=temperature (Hartree)
!!  rtrotter=real trotter parameter
!!  tol=tolerance criteria for stopping density_rec
!!  inf_ucvol=infinitesimal unit cell volume
!!  dim_trott = max(0,2*trotter-1)
!! 
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!! 
!! SIDE EFFECTS
!! 
!! PARENTS
!!      fermisolverec
!!
!! CHILDREN
!!      timab,trottersum
!!
!! NOTES
!!  at this time :
!!       - exppot should be replaced by ?
!!       - coord should be replaced by ?
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine density_rec(an,bn2,rho_out,nrec, &
&                     fermie,tsmear,rtrotter, &
&                     dim_trott,tol,inf_ucvol)

 use defs_basis
 use m_profiling_abi

 use m_time,     only : timab
 use m_rec_tools,only : trottersum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'density_rec'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nrec
 integer,intent(in) :: dim_trott
 real(dp),intent(in) :: fermie,tol,tsmear,inf_ucvol,rtrotter
 real(dp), intent(out) :: rho_out
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)
!Local variables-------------------------------
!not used, debugging purpose only
!for debugging purpose, detailled printing only once for density and ekin
!scalars
 integer, parameter :: minrec = 3
 integer  :: irec
 real(dp) :: beta,mult,prod_b2,error,errold
 real(dp) :: pi_on_rtrotter,twortrotter,exp1,exp2
 complex(dpc) :: cinv2rtrotter,coeef_mu,facrec0
! character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 complex(dpc) :: acc_rho(0:nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)
!**************************************************************************

 call timab(605,1,tsec)
 
!##############################################################
!--Initialisation of metrics 
 mult = two/inf_ucvol   !non-spined system
 beta = one/tsmear

!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*fermie)/(rtrotter))
 exp2 = exp(beta*fermie/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero
!--Initialisation of accumulated density 
 acc_rho = czero  
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero
 

!##############################################################
!--Main loop
 maindo : do irec = 0, nrec
   
!  ######################################################
!  --Density computation
!  !--using the property that: sum_i(bi*c)^2|(z-ai*c)=1/c*sum_i(bi)^2|(z/c-ai)
!  !and for c =exp(-beta*fermie/(two*rtrotter)

   call trottersum(dim_trott,error,prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)
   
   if(irec/=nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit maindo
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)
 
 call timab(605,2,tsec)
 
 end subroutine density_rec
!!***
 
