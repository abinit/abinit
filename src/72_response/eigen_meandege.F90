!{\src2tex{textfont=tt}}
!!****f* ABINIT/eigen_meandege
!! NAME
!! eigen_meandege
!!
!! FUNCTION
!! This routine takes the mean values of the responses
!! for the eigenstates that are degenerate in energy.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)= input eigenresp
!!       eigenrep(2*mband**2*nkpt*nsppol) for first-order derivatives of the eigenvalues
!!       eigenrep(mband*nkpt*nsppol) for Fan or Debye-Waller second-order derivatives of the eigenvalues
!!  mband= maximum number of bands
!!  natom= number of atoms in the unit cell
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for spin-polarized
!!  option= 1 for eigen(1), 2 for eigen(2) - Fan or Debye-Waller
!!
!! OUTPUT
!!  eigenresp_mean(mband*nkpt*nsppol)= eigenresp, averaged over degenerate states
!!
!! PARENTS
!!      dfpt_looppert,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine eigen_meandege(eigen0,eigenresp,eigenresp_mean,mband,nband,nkpt,nsppol,option)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'eigen_meandege'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,option
 integer,intent(in) :: nband(nkpt*nsppol)

!arrays
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)
 real(dp),intent(out) :: eigenresp_mean(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bd2tot_index,iband,ii,ikpt,isppol,nband_k
 real(dp) :: eig0,mean
 character(len=500) :: message
!arrays

! *********************************************************************

!DEBUG
!write(std_out,*)' eigen_meandege : enter '
!ENDDEBUG

 if(option/=1 .and. option/=2)then
   write(message, '(a,i0)' )' The argument option should be 1 or 2, while it is found that option=',option
   MSG_BUG(message)
 end if

 bdtot_index=0 ; bd2tot_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     if(option==1)then
       do iband=1,nband_k
         eigenresp_mean(iband+bdtot_index)=&
&         eigenresp(2*iband-1 + (iband-1)*2*nband_k + bd2tot_index)
       end do
     else if(option==2)then
       do iband=1,nband_k
         eigenresp_mean(iband+bdtot_index)=eigenresp(iband+bdtot_index)
       end do
     end if

!    Treat the case of degeneracies : take the mean of degenerate states
     if(nband_k>1)then
       eig0=eigen0(1+bdtot_index)
       ii=1
       do iband=2,nband_k
         if(eigen0(iband+bdtot_index)-eig0<tol8)then
           ii=ii+1
         else
           mean=sum(eigenresp_mean(iband-ii+bdtot_index:iband-1+bdtot_index))/ii
           eigenresp_mean(iband-ii+bdtot_index:iband-1+bdtot_index)=mean
           ii=1
         end if
         eig0=eigen0(iband+bdtot_index)
         if(iband==nband_k)then
           mean=sum(eigenresp_mean(iband-ii+1+bdtot_index:iband+bdtot_index))/ii
           eigenresp_mean(iband-ii+1+bdtot_index:iband+bdtot_index)=mean
         end if
       end do
     end if

     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2
   end do
 end do
 
!DEBUG
!write(std_out,*)' eigen_meandege : exit'
!ENDDEBUG

end subroutine eigen_meandege
!!***

