!{\src2tex{textfont=tt}}
!!****f* ABINIT/moddiel_csrb
!! NAME
!! moddiel_csrb
!!
!! FUNCTION
!! - Compute a model real-space dielectric functions from the input density
!! and dielar. Based on formula (1) of
!! PRB 47, nb 15, p. 9892 (by C, S, R and B)
!! - formulas for kf and qtf are from
!! J. phys.: condens. matter 2 (1990) 7597-7611
!! - The plasmon frequency/pulsation is from 
!! Theory of the inhomogeneous electron gas (ed. Lundqvist and march) P327
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!
!! OUTPUT
!! rdiemac : real space, diagonal, dielectric permittivity
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp,laplacian
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine moddiel_csrb(dielar,dtset,gprimd,mpi_enreg,rdiemac,rhor_in)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'moddiel_csrb'
 use interfaces_53_ffts
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: dielar(7),gprimd(3,3),rhor_in(dtset%nfft,dtset%nspden)
 real(dp),intent(out) :: rdiemac(dtset%nfft,dtset%nspden)

!Local variables -------------------------
  !real(dp) :: invqtf2(2)
  !real(dp) :: kf(2),wp2(2)  
!scalars
 integer :: ifft
 real(dp) :: alpha
 complex(dpc) :: invqtf2,kf,wp2
! character(len=fnlen) :: filapp
!arrays
 real(dp) :: g2cart(dtset%nfft),qdiemac(2,dtset%nfft,dtset%nspden)
 real(dp) :: rhog(2,dtset%nfft,dtset%nspden),rhor(dtset%nfft,dtset%nspden)
 real(dp) :: xred2(size(dtset%xred_orig,1),size(dtset%xred_orig,2))
 complex(dpc) :: cqdiemac(dtset%nfft,dtset%nspden)
 complex(dpc) :: crhog(dtset%nfft,dtset%nspden)

! *********************************************************************

 xred2=dtset%xred_orig(:,:,1)
 alpha=0.0d0
!presently works only with nspden=1
 rhor=rhor_in
 call laplacian(gprimd,mpi_enreg,dtset%nfft,dtset%nspden,dtset%ngfft,dtset%paral_kgb,g2cart_out=g2cart)
 


 call fourdp(1, rhog(:,:,1), rhor(:,1),-1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0) !

 crhog(:,:)=cmplx(rhog(1,:,:),rhog(2,:,:))

 cqdiemac=zero
 do ifft=1,dtset%nfft
!  wp2(:)=((4.0d0*pi*rhog(:,ifft,1)))
!  kf(:)=((3.0d0*pi*pi*rhog(:,ifft,1))**2)**(1.0d0/6.0d0)
!  invqtf2(:)=(4.0d0*kf(:)/pi)**(-1)
!  qdiemac(:,ifft,1)=one &
!  & + ((dielar(3)-1)**(-1) &
!  & + alpha*g2cart(ifft)*invqtf2(:) &
!  & + g2cart(ifft)*g2cart(ifft)*4.0d0*pi*pi/(4.0d0*wp2(:)))**(-1)
   if(crhog(ifft,1)/=(zero,zero)) then
     wp2=((4.0d0*pi*crhog(ifft,1)))
     kf=((3.0d0*pi*pi*crhog(ifft,1)))**(1.0d0/3.0d0)
     invqtf2=(4.0d0*kf/pi)**(-1)
     cqdiemac(ifft,1)=one &
&     + ((dielar(3)-1)**(-1) &
&     + alpha*g2cart(ifft)*invqtf2 &
&     + g2cart(ifft)*g2cart(ifft)*4.0d0*pi*pi/(4.0d0*wp2))**(-1)
     write(212,*) cqdiemac(ifft,1),wp2,kf,invqtf2,crhog(ifft,1)
   end if
 end do
 qdiemac(1,:,1)=real(cqdiemac(:,1),dp)
 qdiemac(2,:,1)=real(aimag(cqdiemac(:,1)),dp)
 call fourdp(1, qdiemac(:,:,1), rdiemac(:,1),1,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%paral_kgb,0)

end subroutine moddiel_csrb
!!***
