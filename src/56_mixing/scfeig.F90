!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfeig
!!
!! NAME
!! scfeig
!!
!! FUNCTION
!! Compute the largest eigenvalue and eigenvector of the SCF cycle.
!! A brute force algorithm is presently used.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  istep= number of the step in the SCF cycle
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial0(nfft,nspden)= contains vtrial at istep == 1
!!  vtrial(nfft,nspden)= at input, it is the trial potential that gave vresid .
!!       at output, it is an updated trial potential
!!  vrespc(nfft,nspden)=the input preconditioned residual potential
!!  work(nfft,nspden,2)=work space
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scfeig(istep,nfft,nspden,vrespc,vtrial,vtrial0,work,errid,errmess)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfeig'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,nfft,nspden
 integer,intent(out) :: errid
 character(len = 500), intent(out) :: errmess
!arrays
 real(dp),intent(inout) :: vtrial0(nfft,nspden),work(nfft,nspden,2)
 real(dp),intent(inout) :: vrespc(nfft,nspden)
 real(dp), intent(inout) :: vtrial(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,isp
 real(dp) :: eigen_scf,factor,fix_resid,resid_new,resid_old
 character(len=500) :: message

! *************************************************************************
 
 errid = AB7_NO_ERROR

 if(nspden==4)then
   errid = AB7_ERROR_MIXING_ARG
   write(errmess, *) ' scfeig : does not work yet for nspden=4'
   return
 end if

!Set a fixed residual square for normalization of eigenvectors
 fix_resid=1.0d-4

!A few initialisations for the first istep
 if(istep==1)then

   write(message, '(a,es12.4,a,a,a,a,a,a,a)' )&
&   ' scfeig : fixed PC_residual square =',fix_resid,ch10,&
&   '    Note that fixed resid should always be much larger',ch10,&
&   '    than initial PC resid square, still sufficiently',ch10,&
&   '    small to reduce anharmonic effects ',ch10
   call wrtout(std_out,message,'COLL')

!  Compute the preconditioned residual
   resid_old=0.0_dp
   do isp=1,nspden
     do ifft=1,nfft
       resid_old=resid_old+vrespc(ifft,isp)**2
     end do
   end do
   write(message, '(a,es12.4)' )&
&   ' scfeig : initial PC_residual square =',resid_old
   call wrtout(std_out,message,'COLL')
   if(resid_old>1.0d-8)then
     errid = AB7_ERROR_MIXING_ARG
     write(errmess,'(a,a,a,a,a,a,a,a,a,a)') ch10,&
&     ' scfeig : ERROR -',ch10,&
&     '  This value is not good enough to allow',ch10,&
&     '  the computation of the eigenvectors of the SCF cycle.',ch10,&
&     '  It should be better than 1.0d-8 .',ch10,&
&     '  Action : improve the accuracy of your starting wavefunctions.'
     return
   end if

!  Also transfer vtrial in vtrial_old
   vtrial0(:,:)=vtrial(:,:)

!  In order to start the search for eigenvectors,
!  use the tiny residual vector, renormalized
   factor=sqrt(fix_resid/resid_old)
   work(:,:,1)=vrespc(:,:)*factor
   vtrial(:,:)=vtrial0(:,:)+work(:,:,1)

!  If istep is not equal to 1
 else if(istep>=2)then
!  
!  Compute the corresponding operator expectation value
!  And put the residual vector minus the difference
!  between vtrial and vtrial_old
!  (this is actually the action of the operator !) in vect(*,2)
   eigen_scf=0.0_dp
   do isp=1,nspden
     do ifft=1,nfft
       eigen_scf=eigen_scf+&
&       work(ifft,isp,1) * vrespc(ifft,isp)
     end do
   end do

   do isp=1,nspden
     do ifft=1,nfft
       vrespc(ifft,isp)=vrespc(ifft,isp)&
&       +vtrial(ifft,isp)-vtrial0(ifft,isp)
       work(ifft,isp,2)=vrespc(ifft,isp)
     end do
   end do
   eigen_scf=eigen_scf/fix_resid
   write(message, '(a,es12.4,a)' ) &
&   ' scfeig : Operator expectation value ',eigen_scf,' (extremal eigenvalue * diemix)'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
!  
!  Compute residual of vect(*,2)
   resid_new=zero
   do isp=1,min(nspden,2)
     do ifft=1,nfft
       resid_new=resid_new+ work(ifft,isp,2) ** 2
     end do
   end do
   if (nspden==4) then
     do ifft=1,nfft
       resid_new=resid_new+two*(work(ifft,3,2)**2+work(ifft,4,2)**2)
     end do
   end if
   factor=sqrt(fix_resid/resid_new)
   if(eigen_scf<zero) then
     factor=-factor ! the new vector MAY be oposite to the old one
!    if(factor<-one) factor=-factor ! the new vector is not opposed to the old one
   end if
   write(message, '(a,es12.4)' ) &
&   ' scfeig : Inverse of renormalization factor ',one/factor
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(message, '(a,es12.4)' ) &
&   ' scfeig : Convergence criterion value (->0 at convergency) ',one/factor-eigen_scf-one
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   
   work(:,:,1)=work(:,:,2)*factor
   vtrial(:,:)=vtrial0(:,:)+work(:,:,1)
!  End the different istep cases
 end if

end subroutine scfeig
!!***
