!{\src2tex{textfont=tt}}
!!****f* ABINIT/recursion
!! NAME
!! recursion
!! 
!! FUNCTION
!! This routine computes the recursion coefficients and the corresponding 
!! continued fraction to get the density at a point from a fixed potential. 
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (SLeroux,MMancini).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  coordx, coordy, coordz=coordonnees of the computed point
!!  nrec=order of recursion
!!  fermie=fermi energy (Hartree)
!!  tsmear=temperature (Hartree)
!!  dim_trott=dimension of the partial fraction decomposition
!!  rtrotter=trotter parameter (real)
!!  ZT_p=fourier transform of the Green krenel (computed only once in vtorhorec)
!!  typat(:)=type of psp associated to any atom
!!  tol=tolerance criteria for stopping recursion
!!  debug=debugging variable
!!  mpi_enreg=information about MPI paralelisation
!!  nfft=number of points in FFT grid
!!  ngfft=information about FFT
!!  metrec<type(metricrec_type)>=information concerning the infinitesimal metrics
!!  inf_ucvol=infinitesimal unit cell volume
!!  tim_fourdp=time counter for fourdp
!!  natom=number of atoms
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing 
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!  tim= 0 if the time spent in the routine is not taken into account,1 otherwise. For example 
!!  when measuring time for loading  balancing, we don't want to add the time spent in this to the
!!  total time calculation
!! 
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!!  an, bn2 : coefficient given by recursion. 
!!
!! SIDE EFFECTS
!! 
!! PARENTS
!!      first_rec,vtorhorec
!!
!! CHILDREN
!!      fourdp,timab,trottersum,vn_nl_rec
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


subroutine recursion(exppot,coordx,coordy,coordz,an,bn2,rho_out, &
&                    nrec,fermie,tsmear,rtrotter,dim_trott, &
&                    ZT_p, tol,typat, &
&                    nlrec,mpi_enreg,&
&                    nfft,ngfft,metrec,&
&                    tim_fourdp,natom,projec,tim)


 use defs_basis
 use defs_abitypes
 use defs_rectypes
 use m_profiling_abi
 use m_linalg_interfaces

 use m_time,       only : timab
 use m_rec_tools,  only : trottersum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recursion'
 use interfaces_53_ffts
 use interfaces_68_recursion, except_this_one => recursion
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: coordx,coordy,coordz,nfft,nrec,tim
 integer,intent(in) :: tim_fourdp,natom,dim_trott
 real(dp),intent(in) :: fermie,tol,tsmear,rtrotter
 real(dp), intent(out) :: rho_out
 type(MPI_type),intent(in) :: mpi_enreg
 type(nlpsprec_type),intent(in) :: nlrec
 type(metricrec_type),intent(in) :: metrec
!arrays
 integer, intent(in) :: ngfft(18)
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: ZT_p(1:2, 0:nfft-1)
 real(dp), intent(in) :: exppot(0:nfft-1)
 real(dp), intent(in) :: projec(0:,0:,0:,1:,1:)
 real(dp), intent(out) :: an(0:nrec),bn2(0:nrec)
!Local variables-------------------------------
!not used, debugging purpose only
!for debugging purpose, detailled printing only once for density and ekin
!scalars
 integer, parameter :: level = 7, minrec = 3
 integer  :: irec,isign,timab_id,ii
 real(dp) :: switchimu,switchu
 real(dp) :: bb,beta,mult,prod_b2,error,errold
 real(dp) :: inf_ucvol,pi_on_rtrotter,twortrotter,exp1,exp2
 complex(dpc) :: cinv2rtrotter,coeef_mu,facrec0
! character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp) :: inf_tr(3)
 real(dp) :: Zvtempo(1:2, 0:nfft-1)
 real(dp) :: unold(0:nfft-1),vn(0:nfft-1),un(0:nfft-1)
 complex(dpc) :: acc_rho(0:nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)

! *************************************************************************

!--If count time or not
 timab_id = 616; if(tim/=0) timab_id = 606;

 call timab(timab_id,1,tsec)
 
!##############################################################
!--Initialisation of metrics
 inf_ucvol = metrec%ucvol
 inf_tr = metrec%tr
 mult = two/inf_ucvol    !non-spined system

 beta = one/tsmear
!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*fermie)/(rtrotter))
 exp2 = exp(beta*fermie/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

!--Initialisation of  an,bn,un....
 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero

 an = zero; bn2 = zero;  bn2(0) = one
 bb = zero; vn  = zero;  unold  = zero
!--u0 is a Dirac function
 un = zero
 un(coordx+ngfft(1)*(coordy+ngfft(2)*coordz)) = one/sqrt(inf_ucvol)

!--Initialisation of accumulated density 
 acc_rho = czero  
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero

!##############################################################
!--Main loop
 maindo : do irec = 0, nrec
   
!  --Get an and bn2 coef by the lanczos method
   
!  --Computation of exp(-beta*V/8*p)*un or exp(-beta*V/4*p)*un
!  depending on if nl part has to be calculated or not.
   vn = exppot * un   
   
!  --First Non-local psp contribution: (Id+sum_atom int dr1(E(r,r1))vn(r1))
!  --Computation of exp(-beta*V_NL/4*p)*vn
   if(nlrec%nlpsp) then
     call timab(timab_id,2,tsec)
     call vn_nl_rec(vn,natom,typat,ngfft(:3),inf_ucvol,nlrec,projec)
     call timab(timab_id,1,tsec)

!    --Computation of exp(-beta*V/8*p)*vn in nonlocal case
     vn = exppot * vn
   end if !--End if on nlrec%nlpsp
   
!  --Convolution with the Green kernel
!  --FFT of vn
   isign = -1
   call fourdp(1,Zvtempo,vn,isign,mpi_enreg,nfft,ngfft,1,tim_fourdp)

!  --F(T)F(vn)
   do ii = 0,nfft-1
     switchu   = Zvtempo(1,ii)
     switchimu = Zvtempo(2,ii)
     Zvtempo(1,ii) = switchu*ZT_p(1,ii) - switchimu*ZT_p(2,ii)
     Zvtempo(2,ii) = switchu*ZT_p(2,ii) + switchimu*ZT_p(1,ii)
   end do
   
!  --F^-1(F(T)F(vn))
   isign = 1
   call fourdp(1,Zvtempo,vn,isign,mpi_enreg,nfft,ngfft,1,tim_fourdp)

!  --Computation of exp(-beta*V/8*p)*un or exp(-beta*V/4*p)*un
!  depending on if nl part has to be calculated or not.

   vn = inf_ucvol * exppot * vn

!  --Second Non-local psp contribution: (Id+sum_atom E(r,r1))vn   
   if(nlrec%nlpsp) then
     call timab(timab_id,2,tsec)
     call vn_nl_rec(vn,natom,typat,ngfft(:3),inf_ucvol,nlrec,projec)
     call timab(timab_id,1,tsec)

!    --Computation of exp(-beta*V/8*p)*vn in nonlocal case
     vn = exppot * vn
   end if !--End if on nlrec%nlpsp

!  --Multiplication of a and b2 coef by exp(beta*fermie/(two*rtrotter)) must be done in the continued fraction computation  
!  --Computation of a and b2
   an(irec) = inf_ucvol*ddot(nfft,vn,1,un,1)

!  --an must be positive real       
!  --We must compute bn2 and prepare for the next iteration
   if(irec<nrec)then     
     do ii = 0,nfft-1
       switchu = un(ii)
       un(ii) = vn(ii)-an(irec)*un(ii)-bb*unold(ii)
       unold(ii) = switchu
       bn2(irec+1) = bn2(irec+1)+inf_ucvol*un(ii)*un(ii)
     end do
     bb = sqrt(bn2(irec+1))
     un = (one/bb)*un
   end if

!  ######################################################
!  --Density computation
!  density computation is done inside the main looping, juste after the calculus of a and b2, in order to make 
!  it possible to stop the recursion at the needed accuracy, without doing more recursion loop than needed - 
!  further developpement

!  !--using the property that: sum_i(bi*c)^2|(z-ai*c)=1/c*sum_i(bi)^2|(z/c-ai)
!  !and for c =exp(-beta*fermie/(two*rtrotter)

   
   call trottersum(dim_trott,error,prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)
   

   if(irec/=nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)

 
 call timab(timab_id,2,tsec)

 end subroutine recursion
!!***
