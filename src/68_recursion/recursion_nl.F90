!{\src2tex{textfont=tt}}
!!****f* ABINIT/recursion_nl
!! NAME
!! recursion_nl
!! 
!! FUNCTION
!! Given a $|un>$ vector on the real-space grid this routine calculates
!! the density in  $|un>$ by recursion method.
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  trotter=trotter parameter
!!  dim_trott=dimension of the partial fraction decomposition
!!  tsmear=temperature (Hartree)
!!  tol=tolerance criteria for stopping recursion_nl
!!  ngfft=information about FFT(dtset%ngfft a priori different from ngfftrec)
!!  rset<recursion_type> contains all parameter of recursion 
!!  typat(natom)=type of pseudo potential associated to any atom
!!  natom=number of atoms
!!  projec(ngfftrec(1),ngfftrec(2),ngfftrec(3),lmnmax,natom) is the  vector, on the ngfftrec grid containing 
!!  the non-lacal projector $Y_{lm}(r-R_A)f_{lk}(r-R_A)
!!
!! OUTPUT
!!  rho_out=result of the continued fraction multiplied by a multiplicator
!! 
!! SIDE EFFECTS
!!  un(:,:,:)=initial vector on the grid. it is changed in output
!! 
!! PARENTS
!!      nlenergyrec
!!
!! CHILDREN
!!      fourdp,timab,trottersum,vn_nl_rec,wrtout
!!
!! NOTES
!!  at this time :
!!       - need a rectangular box (rmet diagonal matrix)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine recursion_nl(exppot,un,rho_out,rset,ngfft, &
  &                     tsmear,trotter,dim_trott,tol,typat,&
  &                     natom,projec)

 use m_profiling_abi
 
 use defs_basis
 use defs_abitypes
 use defs_rectypes
 use m_rec_tools,only :         trottersum
 use m_linalg_interfaces

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recursion_nl'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_53_ffts
 use interfaces_68_recursion, except_this_one => recursion_nl
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: trotter,natom,dim_trott
 real(dp),intent(in) :: tol,tsmear
 type(recursion_type),intent(in) :: rset
 real(dp), intent(out) :: rho_out
!arrays
 integer,intent(in) ::  typat(natom),ngfft(18)
 real(dp),intent(in) :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
 real(dp),intent(inout) :: un(0:rset%nfftrec-1)
 real(dp),pointer :: projec(:,:,:,:,:)
!Local variables-------------------------------
!scalars
 integer, parameter ::  minrec = 3
 integer  :: irec,isign,ii
 real(dp) :: bb,beta,mult,prod_b2,rtrotter
 real(dp) :: inf_ucvol,pi_on_rtrotter,twortrotter,exp1
 real(dp) :: exp2,error,errold
 real(dp) :: switchu,switchimu
 complex(dpc) :: facrec0,cinv2rtrotter,coeef_mu
 character(len=500) :: msg
 type(mpi_type),pointer:: mpi_loc
!arrays
 real(dp):: tsec(2)
 real(dp):: inf_tr(3)
 real(dp):: an(0:rset%min_nrec),bn2(0:rset%min_nrec)
 real(dp):: vn(0:rset%nfftrec-1)
 real(dp):: unold(0:rset%nfftrec-1)
 real(dp):: Zvtempo(1:2,0:rset%nfftrec-1)
 complex(dpc) :: acc_rho(0:rset%min_nrec)
 complex(dpc) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc) :: N(0:dim_trott),Nold(0:dim_trott)

! *************************************************************************

 call timab(608,1,tsec) !--start time-counter: recursion_nl
 if(rset%debug)then
   msg=' '
   call wrtout(std_out,msg,'COLL')
 end if
 
!##############################################################
 beta = one/tsmear

!--Rewriting the trotter parameter
 rtrotter  = max(half,real(trotter,dp))

!--Initialisation of mpi
 mpi_loc => rset%mpi

!--Initialisation of metrics
 inf_ucvol = rset%inf%ucvol
 inf_tr = rset%inf%tr
 mult = one   !--In the case of the calculus of the NL-energy

!--Initialisation of  an,bn,un....
 N = czero;  D = cone
 facrec0 = cone
 Nold = czero; Dold = czero

 an = zero; bn2 = zero;  bn2(0) = one
 bb = zero; vn  = zero;  unold  = zero

!--Variables for optimisation
 pi_on_rtrotter = pi/rtrotter
 twortrotter = two*rtrotter
 exp1 = exp((beta*rset%efermi)/(rtrotter))
 exp2 = exp(beta*rset%efermi/(twortrotter))
 cinv2rtrotter = cmplx(one/twortrotter,zero,dp)
 coeef_mu = cmplx(one/exp2,zero,dp)

!--Initialisation of accumulated density 
 acc_rho = czero  
!--Initialisation of estimated error
 prod_b2 = twortrotter/exp1
 errold = zero

!##############################################################
!--Main loop
 maindo : do irec = 0, rset%min_nrec
!  --Get an and bn2 coef by the lanczos method
   
!  --Computation of exp(-beta*V/8*p)*un
   vn = exppot * un   
   
!  --First Non-local psp contribution: (Id+sum_atom E(r,r1))vn   
   call timab(608,2,tsec)
   call vn_nl_rec(vn,natom,typat,rset%ngfftrec(:3),inf_ucvol,rset%nl,projec)
   call timab(608,1,tsec)

!  --Computation of exp(-beta*V/8*p)*un
   vn = exppot * vn

!  --Convolution with the Green kernel
!  --FFT of vn
   isign = -1
   call fourdp(1,Zvtempo,vn,isign,mpi_loc,rset%nfftrec,rset%ngfftrec,1,6)

!  --F(T)F(vn)
   do ii = 0,rset%nfftrec-1
     switchu   = Zvtempo(1,ii)
     switchimu = Zvtempo(2,ii)
     Zvtempo(1,ii) = switchu*rset%ZT_p(1,ii) - switchimu*rset%ZT_p(2,ii)
     Zvtempo(2,ii) = switchu*rset%ZT_p(2,ii) + switchimu*rset%ZT_p(1,ii)
   end do

!  --F^-1(F(T)F(vn))
   isign = 1
   call fourdp(1,Zvtempo,vn,isign,mpi_loc,rset%nfftrec,rset%ngfftrec,1,6)

!  --Computation of exp(-beta*V/2*p)*vn
   vn = inf_ucvol * exppot * vn

!  --Second Non-local psp contribution: (Id+sum_atom E(r,r1))vn   
   call timab(608,2,tsec)
   call vn_nl_rec(vn,natom,typat,rset%ngfftrec(:3),inf_ucvol,rset%nl,projec)
   call timab(608,1,tsec)

!  --Computation of exp(-beta*V/8*p)*vn
   vn = exppot * vn


!  --Multiplication of a and b2 coef by exp(beta*fermie/(2.d0*rtrotter)) must be done in the continued fraction computation  
!  --Computation of a and b2
   an(irec) = inf_ucvol*ddot(rset%nfftrec,vn,1,un,1)      !--an must be positive real 

!  --We must compute bn2 and prepare for the next iteration
   if(irec<rset%min_nrec)then    
     do ii = 0,rset%nfftrec-1
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
!  in order to make it possible to stop the recursion_nl at the
!  needed accuracy, without doing more recursion_nl loop than needed further developpement

   call trottersum(dim_trott,error,&
&   prod_b2,pi_on_rtrotter,&
&   facrec0,coeef_mu,exp1,&
&   an(irec),bn2(irec),&
&   N,D,Nold,Dold)


   if(irec/=rset%min_nrec .and. irec>=minrec)then
     if((bn2(irec+1)<tol14).or.(mult*error<tol.and.errold<tol)) exit
   end if
   errold = mult*error
 end do maindo
!--Accumulated density
 rho_out = mult*real(cone-sum(N/D)*cinv2rtrotter,dp)

 call timab(608,2,tsec) !--stop time-counter: recursion_nl

end subroutine recursion_nl
!!***
