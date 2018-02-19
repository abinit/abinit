!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkdenpos
!! NAME
!! mkdenpos
!!
!! FUNCTION
!! Make a ground-state density positive everywhere :
!! when the density (or spin-density) is smaller than xc_denpos,
!! set it to the value of xc_denpos
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components (max. 2)
!!  option=0 if density rhonow is stored as (up,dn)
!!         1 if density rhonow is stored as (up+dn,up)
!!         Active only when nspden=2
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  Input/output
!!  iwarn=At input: iwarn=0 a warning will be printed when rho is negative
!!                  iwarn>0 no warning will be printed out
!!        At output: iwarn is increased by 1
!!  rhonow(nfft,nspden)=electron (spin)-density in real space,
!!     either on the unshifted grid (if ishift==0,
!!     then equal to rhor),or on the shifted grid
!!
!! NOTES
!!  At this stage, rhonow(:,1:nspden) contains the density in real space,
!!  on the unshifted or shifted grid. Now test for negative densities
!!  Note that, ignoring model core charge, as long as boxcut>=2
!!  the shifted density is derivable from the square of a Fourier
!!  interpolated charge density => CANNOT go < 0.
!!  However, actually can go < 0 to within machine precision;
!!  do not print useless warnings in this case, just fix it.
!!  Fourier interpolated core charge can go < 0 due to Gibbs
!!  oscillations; could avoid this by recomputing the model core
!!  charge at the new real space grid points (future work).
!!
!! PARENTS
!!      bethe_salpeter,m_pawxc,mkcore_wvl,posdoppler,poslifetime,posratecore
!!      psolver_rhohxc,rhohxcpositron,rhotoxc,wvl_initro
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkdenpos'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden,option
 integer,intent(inout) :: iwarn
 real(dp),intent(in) :: xc_denpos
!arrays
 real(dp),intent(inout) :: rhonow(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,numneg
 real(dp) :: rhotmp,worst
 character(len=500) :: message
!arrays
 real(dp) :: rho(2)

! *************************************************************************

 numneg=0
 worst=zero

 if(nspden==1)then

!  Non spin-polarized
!$OMP PARALLEL DO PRIVATE(ifft,rhotmp) REDUCTION(MIN:worst) REDUCTION(+:numneg) SHARED(nfft,rhonow)
   do ifft=1,nfft
     rhotmp=rhonow(ifft,1)
     if(rhotmp<xc_denpos)then
       if(rhotmp<-xc_denpos)then
!        This case is probably beyond machine precision considerations
         worst=min(worst,rhotmp)
         numneg=numneg+1
       end if
       rhonow(ifft,1)=xc_denpos
     end if
   end do
 else if (nspden==2) then

!  Spin-polarized

!  rhonow is stored as (up,dn)
   if (option==0) then

!$OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) REDUCTION(MIN:worst) REDUCTION(+:numneg) &
!$OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,1)
       rho(2)=rhonow(ifft,2)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rhonow(ifft,ispden)=xc_denpos
         end if
       end do
     end do

!    rhonow is stored as (up+dn,up)
   else if (option==1) then

!$OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) &
!$OMP&REDUCTION(MIN:worst) REDUCTION(+:numneg) &
!$OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,2)
       rho(2)=rhonow(ifft,1)-rho(1)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rho(ispden)=xc_denpos
           rhonow(ifft,1)=rho(1)+rho(2)
           rhonow(ifft,2)=rho(1)
         end if
       end do
     end do

   end if  ! option

 else
   MSG_BUG('nspden>2 not allowed !')
 end if ! End choice between non-spin polarized and spin-polarized.

 if (numneg>0) then
   if (iwarn==0) then
     write(message,'(a,i0,a,a,a,es10.2,a,e10.2,a,a,a,a)')&
&     'Density went too small (lower than xc_denpos) at ',numneg,' points',ch10,&
&     'and was set to xc_denpos = ',xc_denpos,'. Lowest was ',worst,'.',ch10,&
&     'Likely due to too low boxcut or too low ecut for',' pseudopotential core charge.'
     MSG_WARNING(message)
   end if
   iwarn=iwarn+1
 end if

end subroutine mkdenpos
!!***
