!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhohxcpositron
!! NAME
!! rhohxcpositron
!!
!! FUNCTION
!! Calculate the electrons/positron correlation term for the positron
!!
!! NOTE
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (GJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  nspden=number of spin density components
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  paral_kgb=flag for (k,band,FFT) parallelism
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  ucvol = unit cell volume (Bohr**3)
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  usepaw=flag for PAW
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  electronpositron%e_xc=electron-positron XC energy
!!  electronpositron%e_xcdc=Double-counting electron-positron XC energy
!!  strsxc(6)= contribution of xc to stress tensor (hartree/bohr^3),
!!  vhartr(nfft)=Hartree potential (returned if option/=0 and option/=10)
!!  vxcapn=XC electron-positron XC potential for the positron
!!  vxcavg=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!  kxcapn(nfft,nkxc)=electron-positron XC kernel (returned only if nkxc/=0)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      energy,rhotov,setvtr
!!
!! CHILDREN
!!      mean_fftr,mkdenpos,xcden,xcpositron,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine rhohxcpositron(electronpositron,gprimd,kxcapn,mpi_enreg,nfft,ngfft,nhat,nkxc,nspden,n3xccc,&
&                         paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxcapn,vxcavg,xccc3d,xc_denpos)

 use defs_basis
 use defs_abitypes
 use m_cgtools
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rhohxcpositron'
 use interfaces_41_xc_lowlevel
 use interfaces_56_xc, except_this_one => rhohxcpositron
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nkxc,nspden,n3xccc,paral_kgb,usexcnhat,usepaw
 real(dp),intent(in) :: ucvol,xc_denpos
 real(dp),intent(out) :: vxcavg
 type(electronpositron_type),pointer :: electronpositron
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: nhat(nfft,nspden*usepaw),rhor(nfft,nspden),xccc3d(n3xccc)
 real(dp),intent(out) :: kxcapn(nfft,nkxc),strsxc(6),vhartr(nfft),vxcapn(nfft,nspden)
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: cplex,ierr,ifft,ishift,iwarn,iwarnp,nfftot,ngr,ngrad,nspden_ep
 real(dp) :: exc,excdc,strdiag
 character(len=500) :: message
!arrays
 real(dp),parameter :: qphon(3)=(/0._dp,0._dp,0._dp/)
 real(dp) :: vxcavg_tmp(1)
 real(dp),allocatable :: fxcapn(:),grho2apn(:),rhoe(:,:,:),rhop(:,:),rhotote(:),vxc_ep(:),vxcgr_ep(:)

! *************************************************************************

 if (electronpositron_calctype(electronpositron)/=1) then
   message = 'Only electronpositron%calctype=1 allowed !'
   MSG_BUG(message)
 end if

 if (nkxc>3) then
   message = 'nkxc>3 (Kxc for GGA) not yet implemented !'
   MSG_ERROR(message)
 end if

!Hartree potential of the positron is zero
 vhartr=zero

!Some allocations/inits
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 ngr=0;if (ngrad==2) ngr=nfft
 ABI_ALLOCATE(fxcapn,(nfft))
 ABI_ALLOCATE(grho2apn,(ngr))
 nspden_ep=1;cplex=1;ishift=0
 iwarn=0;iwarnp=1

!Compute total electronic density
 ABI_ALLOCATE(rhotote,(nfft))
 rhotote(:)=electronpositron%rhor_ep(:,1)
 if (n3xccc>0) rhotote(:)=rhotote(:)+xccc3d(:)
 if (usepaw==1.and.usexcnhat==0) rhotote(:)=rhotote(:)-electronpositron%nhat_ep(:,1)

!Extra total electron/positron densities; compute gradients for GGA
 ABI_ALLOCATE(rhoe,(nfft,nspden_ep,ngrad**2))
 ABI_ALLOCATE(rhop,(nfft,nspden_ep))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,paral_kgb,qphon,rhotote,rhoe)
 if (ngrad==2) grho2apn(:)=rhoe(:,1,2)**2+rhoe(:,1,3)**2+rhoe(:,1,4)**2
 rhop(:,1)=rhor(:,1);if (usepaw==1.and.usexcnhat==0) rhop(:,1)=rhop(:,1)-nhat(:,1)
 ABI_DEALLOCATE(rhotote)

!Make the densities positive
 call mkdenpos(iwarn ,nfft,nspden_ep,1,rhoe(:,1,1),xc_denpos)
 if (.not.electronpositron%posdensity0_limit) then
   call mkdenpos(iwarnp,nfft,nspden_ep,1,rhop,xc_denpos)
 end if

!Compute electron-positron Vxc_pos, Vxc_el, Fxc, Kxc, ...
 ABI_ALLOCATE(vxc_ep,(nfft))
 ABI_ALLOCATE(vxcgr_ep,(ngr))
 if (nkxc==0) then
   call xcpositron(fxcapn,grho2apn,electronpositron%ixcpositron,ngr,nfft,electronpositron%posdensity0_limit,&
&   rhoe(:,1,1),rhop(:,1),vxc_ep,vxcgr_ep,vxcapn)
 else
   call xcpositron(fxcapn,grho2apn,electronpositron%ixcpositron,ngr,nfft,electronpositron%posdensity0_limit,&
&   rhoe(:,1,1),rhop(:,1),vxc_ep,vxcgr_ep,vxcapn,dvxce=kxcapn)
 end if
 ABI_DEALLOCATE(rhoe)
 ABI_DEALLOCATE(vxc_ep)
 ABI_DEALLOCATE(vxcgr_ep)
 ABI_DEALLOCATE(grho2apn)

!Store Vxc and Kxc according to spin components
 if (nspden>=2) vxcapn(:,2)=vxcapn(:,1)
 if (nspden==4) vxcapn(:,3:4)=zero
 if (nkxc==3) then
   kxcapn(:,1)=two*kxcapn(:,1)
   kxcapn(:,2)=kxcapn(:,1)
   kxcapn(:,3)=kxcapn(:,1)
 end if

!Compute XC energies and contribution to stress tensor
 electronpositron%e_xc  =zero
 electronpositron%e_xcdc=zero
 strdiag=zero
 nfftot=PRODUCT(ngfft(1:3))
 do ifft=1,nfft
   electronpositron%e_xc  =electronpositron%e_xc  +fxcapn(ifft)
   electronpositron%e_xcdc=electronpositron%e_xcdc+vxcapn(ifft,1)*rhor(ifft,1)
!  strdiag=strdiag+fxcapn(ifft)   ! Already stored in rhotoxc !
   strdiag=strdiag-vxcapn(ifft,1)*rhop(ifft,1)
 end do
 if (usepaw==1.and.usexcnhat==0) then
   do ifft=1,nfft
     electronpositron%e_xcdc=electronpositron%e_xcdc-vxcapn(ifft,1)*nhat(ifft,1)
   end do
 end if
 electronpositron%e_xc  =electronpositron%e_xc  *ucvol/dble(nfftot)
 electronpositron%e_xcdc=electronpositron%e_xcdc*ucvol/dble(nfftot)
 strdiag=strdiag/dble(nfftot)
 ABI_DEALLOCATE(fxcapn)
 ABI_DEALLOCATE(rhop)

!Reduction in case of parallelism
 if(mpi_enreg%paral_kgb==1)then
   if(paral_kgb/=0)then
     exc=electronpositron%e_xc;excdc=electronpositron%e_xcdc
     call xmpi_sum(exc  ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(excdc,mpi_enreg%comm_fft,ierr)
     electronpositron%e_xc=exc;electronpositron%e_xcdc=excdc
     call xmpi_sum(strsxc,mpi_enreg%comm_fft,ierr)
   end if
 end if

!Store stress tensor
 strsxc(1:3)=strdiag
 strsxc(4:6)=zero

!Compute vxcavg
 call mean_fftr(vxcapn(:,1),vxcavg_tmp,nfft,nfftot,1)
 vxcavg=vxcavg_tmp(1)

end subroutine rhohxcpositron
!!***
