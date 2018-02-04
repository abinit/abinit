!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourier_interpol
!! NAME
!! fourier_interpol
!!
!! FUNCTION
!!  Perform a Fourier interpolation. Just a wrapper for transgrid, the table giving the correspondence 
!!  between the coarse and the mesh FFT grids are constructed inside the routine. This allows to specify an 
!!  arbitrary FFT mesh to be used for the interpolation. Besides the routine works also in 
!!  case of NC calculations since it does not require Pawfgr.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! cplex=1 if rhor[f] is real, 2 if rhor[f] is complex
!! MPI_enreg<MPI_type>=Information about MPI parallelization
!! nspden=number of spin-density components
!! nfft_in =number of points in the input FFT box (WARNING no FFT parallelism)
!! nfft_out=number of points in the output FFT box
!! ngfft_in(18)=all needed information about 3D FFT, for the input grid
!! ngfft_out(18) =all needed information about 3D FFT, for the output grid
!! optin= 0: input density/potential is taken from rhor_in(:,nspden)
!!        1: input density/potential is taken from rhog_in(:)     (ispden=1)
!!                                             and rhor_in(:,2:4) (ispden=2,3,4)
!! optout= 0: output density/potential is given in r space in rhor_out(:,nspden)
!!         1: output density/potential is given in r space in rhor_out(:,nspden)
!!                                          and in g space in rhog_out(:)
!! ngfft_asked(18)=All info on the required FFT mesh.
!! paral_kgb=Flag related to the kpoint-band-fft parallelism (TODO not implemented)
!!  
!! OUTPUT
!!  rhor_out(cplex*nfft_out,nspden)=output density/potential in r space on the required FFT mesh.
!!  if optout=1:
!!   rhog_out(2,nfftc)=Fourier transform of output density/potential on the coarse grid
!!
!! PARENTS
!!      m_dvdb,m_ioarr,m_qparticles
!!
!! CHILDREN
!!      indgrid,pawfgr_destroy,transgrid
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fourier_interpol(cplex,nspden,optin,optout,nfft_in,ngfft_in,nfft_out,ngfft_out,&
& paral_kgb,MPI_enreg,rhor_in,rhor_out,rhog_in,rhog_out)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawfgr, only : pawfgr_type, pawfgr_destroy, indgrid

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourier_interpol'
 use interfaces_65_paw, except_this_one => fourier_interpol
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nspden,optin,optout
 integer,intent(in) :: nfft_in,nfft_out,paral_kgb
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfft_in(18),ngfft_out(18)
 real(dp),intent(inout) :: rhor_in(cplex*nfft_in,nspden)
 real(dp),intent(inout) :: rhog_in(2,nfft_in)
 real(dp),intent(out) :: rhor_out(cplex*nfft_out,nspden)
 real(dp),intent(out) :: rhog_out(2,nfft_out)

!Local variables ---------------------------------------
!scalars
 integer :: nfftf,nfftc,nfftc_tot,nfftf_tot,optgrid
 logical :: ltest
 type(Pawfgr_type) :: Pawfgr
!arrays
 integer :: ngfftc(18),ngfftf(18)

! *************************************************************************

!=== FFT parallelism not implemented ===
 ABI_CHECK(paral_kgb==0,'paral_kgb/=0 not implemented')

 ltest= ALL( ngfft_in(7:) == ngfft_out(7:) )
 ABI_CHECK(ltest,'ngfftf_in(7:18)/=ngfftf_out(7:18)')

!================================
!=== Which one is the coarse? ===
!================================
 if (nfft_out>=nfft_in) then 
   ! From coarse to fine grid. If meshes are equivalent, call transgrid anyway because of optout, optin.
   nfftf    =nfft_out 
   ngfftf(:)=ngfft_out(:)
   nfftf_tot =PRODUCT(ngfft_out(1:3))

   nfftc    =nfft_in 
   ngfftc(:)=ngfft_in(:)
   nfftc_tot =PRODUCT(ngfft_in (1:3))

   Pawfgr%usefinegrid=1 
   if (ALL(ngfft_in(1:3)==ngfft_out(1:3))) Pawfgr%usefinegrid=0 
   optgrid=1

 else 
   ! From fine towards coarse.
   nfftf    =nfft_in 
   ngfftf(:)=ngfft_in(:)
   nfftf_tot =PRODUCT(ngfft_in(1:3))

   nfftc    =nfft_out 
   ngfftc(:)=ngfft_out(:)
   nfftc_tot =PRODUCT(ngfft_out (1:3))
   Pawfgr%usefinegrid=1 
   optgrid=-1
 end if

 ABI_MALLOC(Pawfgr%coatofin,(nfftc_tot))
 ABI_MALLOC(Pawfgr%fintocoa,(nfftf_tot))

 call indgrid(Pawfgr%coatofin,Pawfgr%fintocoa,nfftc_tot,nfftf_tot,ngfftc,ngfftf)

 Pawfgr%mgfft =MAXVAL (ngfftf(1:3))
 Pawfgr%nfft  =PRODUCT(ngfftf(1:3)) ! no FFT parallelism!
 Pawfgr%ngfft =ngfftf

 Pawfgr%mgfftc =MAXVAL (ngfftc(1:3))
 Pawfgr%nfftc  =PRODUCT(ngfftc(1:3)) ! no FFT parallelism!
 Pawfgr%ngfftc =ngfftc

 if (optgrid==1) then 
   call transgrid(cplex,MPI_enreg,nspden,optgrid,optin,optout,paral_kgb,Pawfgr,rhog_in ,rhog_out,rhor_in ,rhor_out)
 else
   call transgrid(cplex,MPI_enreg,nspden,optgrid,optin,optout,paral_kgb,Pawfgr,rhog_out,rhog_in ,rhor_out,rhor_in)
 end if

 call pawfgr_destroy(Pawfgr)
 
end subroutine fourier_interpol
!!***
