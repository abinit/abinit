!{\src2tex{textfont=tt}}
!!****f* ABINIT/gammapositron_fft
!! NAME
!! gammapositron_fft
!!
!! FUNCTION
!! Compute positron electron-positron enhancement factor on a real space FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 2013-2018 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  igamma=type of enhancement factor:
!!    -1:  gamma=one (test)
!!     1:  Boronski and Nieminen [1]
!!     2:  Boronski and Nieminen, RPA limit [1]
!!     3:  Sterne and Kaiser [2]
!!     4:  Puska, Seitsonen and Nieminen [3]
!!  mpi_enreg=information about MPI parallelization
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nfft=number of FFT grid points
!!  ngfft(18)=contain all needed information about 3D FFT
!!  rhor_e(nfft)=real space density electronic density (total density)
!!  rhor_p(nfft)=real space density positronic density
!!  xccc3d(n3xccc)=3D core electron density for XC core correction
!!
!! OUTPUT
!!  gamma(nfft,2)=electron-positron enhancement factor,
!!                gamma(:,1): using total   electronic density
!!                gamma(:,2): using valence electronic density
!!
!! NOTES
!!  The input densities (rhor_e and rhor_p) should be positive everywhere
!!  (call mkdenpos routine before entering this one)
!!
!! PARENTS
!!      posdoppler,poslifetime
!!
!! CHILDREN
!!      gammapositron,xcden
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gammapositron_fft(electronpositron,gamma,gprimd,igamma,mpi_enreg,&
&                            n3xccc,nfft,ngfft,rhor_e,rhor_p,xccc3d)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 
 use m_electronpositron

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gammapositron_fft'
 use interfaces_56_xc, except_this_one => gammapositron_fft
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: igamma,n3xccc,nfft
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),rhor_e(nfft),rhor_p(nfft),xccc3d(n3xccc)
 real(dp),intent(out) :: gamma(nfft,2)

!Local variables-------------------------------
!scalars
 integer :: cplex,ishift,ngr,ngrad,nspden_ep,usecore
!arrays
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/)
 real(dp),allocatable :: grhocore2(:),grhoe2(:),rhoc(:,:,:),rhoe(:,:,:)

! *************************************************************************

!Several useful constants
 usecore=n3xccc/nfft
 cplex=1;ishift=0;ngrad=1;nspden_ep=1
 if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 ngr=0;if (ngrad==2) ngr=nfft
 
!Allocate several arrays
 ABI_ALLOCATE(rhoe,(nfft,nspden_ep,ngrad**2))
 ABI_ALLOCATE(grhoe2,(ngr))
 ABI_ALLOCATE(grhocore2,(ngr*usecore))

!Store electronic density and its gradients
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,&
& mpi_enreg%paral_kgb,qphon,rhor_e,rhoe)

!Compute squared gradient of the electronic density
 if (ngrad==2) then
   grhoe2(:)=rhoe(:,1,2)**2+rhoe(:,1,3)**2+rhoe(:,1,4)**2
   if (usecore>0) then
     ABI_ALLOCATE(rhoc,(nfft,1,ngrad**2))
     call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,&
&     mpi_enreg%paral_kgb,qphon,xccc3d,rhoc)
     grhocore2(:)=rhoc(:,1,2)**2+rhoc(:,1,3)**2+rhoc(:,1,4)**2
     ABI_DEALLOCATE(rhoc)
   end if
 end if

!Compute enhancement factor on FFT grid
 call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,nfft,xccc3d,&
& rhoe(:,1,1),rhor_p,usecore)

!Release temporary memory
 ABI_DEALLOCATE(rhoe)
 ABI_DEALLOCATE(grhoe2)
 ABI_DEALLOCATE(grhocore2)

end subroutine gammapositron_fft
!!***
