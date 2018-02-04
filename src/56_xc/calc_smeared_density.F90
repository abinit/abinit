!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_smeared_density
!! NAME
!! calc_smeared_density
!!
!! FUNCTION
!! Calculate the smeared charge density on the an FFT grid in real space.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  rhor           - the density on the real space FFT grid
!!  kappa_strategy - integer indicating strategy for the smearing parameter
!!                   1 - > fixed kappa
!!                   2 - > take kappa as n-dependent
!!  nfftf - contains all needed information about the 3D FFT grid
!!  ngfftf(18) - contains all needed information about the fine 3D FFT grid
!!  npw      - number of plane waves
!!  mpi_enreg - Datatype with parallellisation info
!!  paral_kgb - band parallelism variable
!!  kappa_in  - optional variable setting a fixed value of kappa
!!
!! OUTPUT
!!  rhotilder - the smeared density on the real space FFT grid
!!
!! NOTES
!!   This routine transforms a given density rhor(R) -> rhog(G) and calculates
!!   smeared density. This is given by the formula:
!!   
!!   \tilde{n}(\mathbf{G}) = n(\mathbf{G})\frac{{\kappa}^2}{{\kappa}^2+|\mathbf{G}|^2}
!!
!!   which corresponds to a convolution of the density with a Yukawa function
!!   in real space:
!!  
!!   \tilde{n}(\mathbf{r})=\int \mathrm{d}r'\,\frac{{\kappa}^2}{4\pi} \\
!!                         \frac{e^{-{\kappa}|\mathbf{r}-\mathbf{r}'|}}{|\mathbf{r}-\mathbf{r'}|} \\
!!                         n(\mathbf{r}')
!!
!!   When the convolution has been performed in reciprocal space, the function
!!   is transformed back to real space, giving the smeared density in real space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine calc_smeared_density(rhor,kappa_strategy,rhotilder,nfftf,ngfftf,npw,&
&              gvec,gprimd,ucvol,mpi_enreg,paral_kgb,kappa_in)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_fft_mesh,       only : g2ifft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_smeared_density'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: kappa_strategy,nfftf,npw,paral_kgb
!integer, optional, intent(in) :: prtvol
 real(dp), intent(in) :: ucvol
 real(dp), intent(in), optional :: kappa_in
!arrays
 integer,intent(in) :: ngfftf(18),gvec(3,npw)
 real(dp), intent(in) :: gprimd(3,3)
 real(dp), intent(inout) :: rhor(nfftf)
 real(dp), intent(out) :: rhotilder(nfftf)
! types
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables ------------------------------
!scalars
 integer :: ig,igp,ii,gmgp_idx,ierr
 real(dp), parameter :: two_pi_sq = two_pi*two_pi
 real(dp) :: inv_kappasq,abs_gmgp_sq,yukawa_factor,yukawa_denom
!arrays
 real(dp) :: gmet(3,3)
 integer :: gmgp(3)
 real(dp), allocatable :: rhog_dp(:,:),rhotildeg_dp(:,:)
 character(len=500) :: msg

!*************************************************************************

 DBG_ENTER("COLL")

 ABI_ALLOCATE(rhog_dp,(2,nfftf))
 ABI_ALLOCATE(rhotildeg_dp,(2,nfftf))
 rhog_dp=zero; rhotildeg_dp=zero

!write(std_out,*) ' calc. sm. den., npw=',npw
!write(std_out,*) ' calc. sm. den., gprimd=',gprimd
!write(std_out,*) ' calc. sm. den., ngfftf(1:3)=',ngfftf(1:3)
!write(std_out,*) ' calc. sm. den., nfftf=',nfftf

 if (present(kappa_in)) write(std_out,*) ' calc. sm. den., kappa_in=',kappa_in

 call fourdp(1,rhog_dp,rhor,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)

!Compute reciprocal space metrics
 do ii=1,3
   gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&   gprimd(2,ii)*gprimd(2,:)+&
&   gprimd(3,ii)*gprimd(3,:)
 end do

!Check which strategy for determining kappa is used
 select case(kappa_strategy)

 case(1) ! kappa is a constant, kappa_in is provided in the routine call

   if (present(kappa_in)) then
     inv_kappasq = 1.0_dp/(kappa_in*kappa_in)
   else
     MSG_ERROR('   Argument kappa_in must be set with kappa_strategy=1')
   end if

!    Perform the multiplication
   ierr=ierr+1
   do ig=1,npw
     do igp=1,npw
!        Calculate |G-G'|
       gmgp(:) = gvec(:,ig)-gvec(:,igp)
!        Make sure the vector is not outside the FFT box
       if (ANY (gmgp(:)>ngfftf(1:3)/2 .or. gmgp(:)<-(ngfftf(1:3)-1)/2) ) then
         ierr=ierr+1
         cycle
       end if
       abs_gmgp_sq = two_pi_sq*dot_product(gmgp,MATMUL(gmet,gmgp))
!        Calculate Yukawa factor kappa^2/(kappa^2+|G-G'|^2)
!        = 1/(1+(|G-G'|/kappa)^2)
       yukawa_denom = one + abs_gmgp_sq*inv_kappasq
       yukawa_factor = one/yukawa_denom
!        Calculate smeared density in G-space i.e. n(G-G')*kappa^2/(kappa^2+|G-G'|^2)
       gmgp_idx = g2ifft(gmgp,ngfftf)
!        rhotildeg_dp(:,igfft(ig,igp)) = rhog_dp(:,igfft(ig,igp))*yukawa_factor
       rhotildeg_dp(:,gmgp_idx) = rhog_dp(:,gmgp_idx)*yukawa_factor
     end do
   end do

   if (ierr/=0) then 
     write(msg,'(a,i4,3a)')&
&     'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     'Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if

 case(2) ! kappa is a function of the density
   MSG_ERROR('kappa_strategy=2 not coded yet')

 case default
   MSG_ERROR('error in kappa_strategy')
 end select

!Transform the smeared density back to real space

 call fourdp(1,rhotildeg_dp,rhotilder,1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)

!Debug
 abs_gmgp_sq = SUM(rhotilder(:))*ucvol/nfftf
 write(msg,'(a,ES29.10E3)') '   Integral of rhotilde(r):',abs_gmgp_sq
 MSG_WARNING(msg)
 abs_gmgp_sq = SUM(ABS(rhor(:)))*ucvol/nfftf
 write(msg,'(a,ES29.10E3)') '       Integral of  rho(r):',abs_gmgp_sq
 MSG_WARNING(msg)
 abs_gmgp_sq = SUM(rhotilder(:)-rhor(:))*ucvol/nfftf
 write(msg,'(a,ES29.10E3)') '   Integral of difference rhotilde(r)-rho(r):',abs_gmgp_sq
 MSG_WARNING(msg)
 abs_gmgp_sq = SUM(ABS(rhotilder(:)-rhor(:)))*ucvol/nfftf
 write(msg,'(a,ES29.10E3)') '   Integral of abs. difference |rhotilde(r)-rho(r)|:', abs_gmgp_sq
 MSG_WARNING(msg)

 ABI_DEALLOCATE(rhog_dp)
 ABI_DEALLOCATE(rhotildeg_dp)

 DBG_EXIT("COLL")

end subroutine calc_smeared_density
!!***

