!{\src2tex{textfont=tt}}
!!****f* ABINIT/recip_ylm
!! NAME
!! recip_ylm
!!
!! FUNCTION
!! Project input wavefunctions (real space) on to Ylm
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bess_fit(mpw,nradintmax,ll) = bessel functions for ll, splined
!!   with arguments $2 \pi |k+G| \Delta r$, for all G vectors in sphere
!!   and all points on radial grid.
!!  cgcband(2,npw_k)=wavefunction in recip space
!!  istwfk= storage mode of cgcband
!!  mlang=maximum angular momentum
!!  mpi_enreg=information about MPI parallelization
!!  natsph=number of atoms around which ang mom projection has to be done 
!!  npw_k=number of plane waves for kpt
!!  nradint(natsph)=number of points on radial real-space grid for a given atom
!!  nradintmax=dimension of rint array
!!  ph3d(2,npw_k,natsph)=3-dim structure factors, for each atom and plane wave.
!!  prtsphere= if 1, print a complete analysis of the angular momenta in atomic spheres
!!  rint(nradintmax) = points on radial real-space grid for integration
!!  rmax(natsph)=maximum radius for real space integration sphere
!!  ucvol=unit cell volume in bohr**3.
!!  ylm(mpw,mlang*mlang)=real spherical harmonics for each G and k point
!!  znucl_sph(natsph)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  sum_1ll_1atom(mlang,natsph)= projected scalars for each atom and ang. mom.
!!
!! NOTES
!!  ph3d atoms are ordered with atindx -> by typat
!!  The atoms have to be reordered !
!!  factor monetom = (-1)**ll has been removed - was it missing in one of the equations?
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      atomdata_from_znucl,dotprod_g
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine recip_ylm (bess_fit,cgcband,istwfk,&
& nradint,nradintmax,mlang,mpi_enreg,mpw,natsph,&
& npw_k,ph3d,prtsphere,rint,rmax,sum_1ll_1atom,sum_1lm_1atom,ucvol,ylm,znucl_sph)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_cgtools
 use m_atomdata

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recip_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,mlang,mpw,natsph,npw_k,nradintmax
 integer,intent(in) :: prtsphere
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nradint(natsph)
 real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang),cgcband(2,npw_k)
 real(dp),intent(in) :: ph3d(2,npw_k,natsph),rint(nradintmax)
 real(dp),intent(in) :: rmax(natsph),ylm(mpw,mlang*mlang)
 real(dp),intent(in) :: znucl_sph(natsph)
 real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
 real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)

!Local variables-------------------------------
!scalars
 integer :: illmm, iat, ipw, ixint, ll, mm, option
 real(dp) :: doti, dotr
 real(dp) :: sum_all, integ
 real(dp) :: dr, llsign1, llsign2
 type(atomdata_t) :: atom
!arrays
 real(dp) :: sum_1atom(natsph),sum_1ll(mlang),sum_1lm(mlang*mlang)
 real(dp) :: tmppsia(2,npw_k),tmppsim(2,npw_k)

 integer, allocatable :: ilang(:)
 integer,  allocatable :: reylmind(:)
 integer,  allocatable :: imylmind(:)

 real(dp), allocatable :: coef1(:)
 real(dp), allocatable :: coef2(:)
 real(dp), allocatable :: vect(:,:)
 real(dp), allocatable :: mmsign(:)
 real(dp), allocatable :: rint2(:)

! *************************************************************************

!metric has been called in calling routine.

 sum_1lm_1atom(:,:) = zero
 sum_1ll_1atom(:,:) = zero

!option for dotprod_g
 option=2

 ABI_ALLOCATE(rint2, (nradintmax))
 do ixint = 1, nradintmax
   rint2(ixint) = rint(ixint)**2
 end do

 ABI_ALLOCATE(vect,(2,npw_k))
 vect = zero

 ABI_ALLOCATE(ilang,    (mlang**2))
 ABI_ALLOCATE(reylmind, (mlang**2))
 ABI_ALLOCATE(imylmind, (mlang**2))
 ABI_ALLOCATE(mmsign,   (mlang**2))
 ABI_ALLOCATE(coef1,   (mlang**2))
 ABI_ALLOCATE(coef2,   (mlang**2))
 do ll=0,mlang-1
   llsign1 = one
   llsign2 = zero
   if (mod(ll,2) == 1) then
     llsign1 = zero
     llsign2 = one
   end if
   do mm = -ll, -1
     illmm = (ll+1)**2-ll+mm
     ilang(illmm) = ll+1
     coef1(illmm) = llsign1 ! 1 for even ll channels
     coef2(illmm) = llsign2 ! 1 for odd ll channels
     reylmind(illmm) = (ll+1)**2-ll-mm
     imylmind(illmm) = (ll+1)**2-ll+mm
     mmsign(illmm) = -one
   end do

   do mm = 0, ll
     illmm = (ll+1)**2-ll+mm
     ilang(illmm) = ll+1
     coef1(illmm) = llsign1 ! 1 for even ll channels
     coef2(illmm) = llsign2 ! 1 for odd ll channels
     reylmind(illmm) = (ll+1)**2-ll+mm
     imylmind(illmm) = (ll+1)**2-ll-mm
     mmsign(illmm) = one
   end do
 end do
 if (istwfk == 1) then
   coef1 = one
   coef2 = one
 end if

!Big loop  on all atoms
 do iat=1,natsph
   dr = rmax(iat)/(nradint(iat)-1)

!  Temporary arrays for part of psi which depends only on iat
   do ipw=1,npw_k
     tmppsia(1,ipw) = cgcband(1,ipw) * ph3d(1,ipw,iat) &
&     - cgcband(2,ipw) * ph3d(2,ipw,iat)
     tmppsia(2,ipw) = cgcband(1,ipw) * ph3d(2,ipw,iat) &
&     + cgcband(2,ipw) * ph3d(1,ipw,iat)
   end do

   do illmm=1, mlang*mlang

!    tmppsim = temporary arrays for part of psi which doesnt depend on ixint
!    Take into account the fact that ylm are REAL spherical harmonics, see initylmg.f
!    
!    For time-reversal states, detailed treatment show that only the real or imaginary
!    part of tmppsia is needed here, depending on l being even or odd: only one of the coef is 1, the other 0
     do ipw=1,npw_k
!      to get PDOS for real spherical harmonics, may be sufficient to multiply here by ylm instead of linear combination
!      tmppsim(1,ipw) = tmppsia(1,ipw)*ylm(ipw,illmm)
!      tmppsim(2,ipw) = tmppsia(2,ipw)*ylm(ipw,illmm)
!      TODO: check the prefactor part for istwfk /= 1! Could also be incorrect if we go to real spherical harmonics
       tmppsim(1, ipw) =  coef1(illmm) * tmppsia(1, ipw) * ylm(ipw, reylmind(illmm)) &
&       + mmsign(illmm) * coef2(illmm) * tmppsia(2, ipw) * ylm(ipw, imylmind(illmm))

       tmppsim(2, ipw) =  coef2(illmm) * tmppsia(2, ipw) * ylm(ipw, reylmind(illmm)) &
&       - mmsign(illmm) * coef1(illmm) * tmppsia(1, ipw) * ylm(ipw, imylmind(illmm))
     end do

     integ = zero
     do ixint = 1, nradint(iat)
       vect(1, 1:npw_k) = bess_fit(1:npw_k, ixint, ilang(illmm))
       call dotprod_g(dotr, doti, istwfk, npw_k, option, vect, tmppsim, mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)

!      Multiply by r**2 and take norm, integrate
!      MJV 5.5.2012 removed call to simpson_int - just need last value, no need to allocate full space for primitive and integrand

       if (ixint ==1 .or. ixint == nradint(iat)) then
         integ = integ + 0.375_dp * rint2(ixint) * (dotr**2 + doti**2)
       else if (ixint ==2 .or. ixint == nradint(iat)-1) then
         integ = integ + 1.166666666666666666666666667_dp * rint2(ixint) * (dotr**2 + doti**2)
       else if (ixint ==3 .or. ixint == nradint(iat)-2) then
         integ = integ + 0.958333333333333333333333333_dp * rint2(ixint) * (dotr**2 + doti**2)
       else
         integ = integ + rint2(ixint) * (dotr**2 + doti**2)
       end if
     end do ! ixint
     integ = integ * dr

!    NOTE : could exploit full r dependency of integ
!    which is calculated in the call to simpson

     sum_1lm_1atom(illmm, iat)        = sum_1lm_1atom(illmm, iat)        + integ
     sum_1ll_1atom(ilang(illmm), iat) = sum_1ll_1atom(ilang(illmm), iat) + integ

   end do !  illmm
 end do ! iat

 ABI_DEALLOCATE(coef1)
 ABI_DEALLOCATE(coef2)
 ABI_DEALLOCATE(reylmind)
 ABI_DEALLOCATE(imylmind)
 ABI_DEALLOCATE(mmsign)
 ABI_DEALLOCATE(ilang)

 ABI_DEALLOCATE(vect)
 ABI_DEALLOCATE(rint2)

!MJV 5.5.2012: removed 1/sqrt(2) above in tmppsim and (4 pi)^2 in integrand - just multiply by 8 pi^2 here
!Normalize with unit cell volume
 sum_1lm_1atom(:,:) = eight * pi**2 * sum_1lm_1atom(:,:) / ucvol
 sum_1ll_1atom(:,:) = eight * pi**2 * sum_1ll_1atom(:,:) / ucvol

!Output
 if(prtsphere==1)then
   sum_1ll = zero
   sum_1lm = zero
   sum_1atom = zero
   do iat=1,natsph
     sum_1atom(iat) = sum(sum_1lm_1atom(:,iat))
     sum_1ll(:)=sum_1ll(:)+sum_1ll_1atom(:,iat)
     sum_1lm(:)=sum_1lm(:)+sum_1lm_1atom(:,iat)
   end do
   sum_all = sum(sum_1atom)

   write(std_out,'(a)' ) ' Angular analysis '
   do iat=1,natsph
     call atomdata_from_znucl(atom, znucl_sph(iat) )
     write(std_out,'(a)' ) ' '
     write(std_out,'(a,i3,a,a,a,f10.6)' )&
&     ' Atom # ',iat, ' is  ',  atom%symbol,', in-sphere charge =',sum_1atom(iat)
     do ll=0,mlang-1
       write(std_out,'(a,i1,a,f9.6,a,9f6.3)' )&
&       ' l=',ll,', charge=',sum_1ll_1atom(ll+1,iat),&
&       ', m=-l,l splitting:',sum_1lm_1atom(1+ll**2:(ll+1)**2,iat)
     end do ! ll
   end do ! iat
   write(std_out,'(a,a)') ch10,' Sum of angular contributions for all atomic spheres '
   do ll=0,mlang-1
     write(std_out,'(a,i1,a,f9.6,a,f9.6)' )&
&     ' l=',ll,', charge =',sum_1ll(ll+1),' proportion =',sum_1ll(ll+1)/sum_all
   end do
   write(std_out,'(a,a,f10.6)' ) ch10,' Total over all atoms and l=0 to 4 :',sum_all
   write(std_out,'(a)' ) ' '
 end if

end subroutine recip_ylm
!!***
