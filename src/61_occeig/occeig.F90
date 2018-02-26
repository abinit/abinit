!{\src2tex{textfont=tt}}
!!****f* ABINIT/occeig
!! NAME
!! occeig
!!
!! FUNCTION
!! For each pair of active bands (m,n), generates ratios
!! that depend on the difference between occupation numbers
!! and eigenvalues.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  doccde_k(nband_k)=derivative of occ_k wrt the energy
!!  doccde_kq(nband_k)=derivative of occ_kq wrt the energy
!!  eig0_k(nband_k)=GS eigenvalues at k
!!  eig0_kq(nband_k)=GS eigenvalues at k+q
!!  nband_k=number of bands
!!  occopt=option for occupancies
!!  occ_k(nband_k)=occupation number for each band at k
!!  occ_kq(nband_k)=occupation number for each band at k+q
!!
!! OUTPUT
!!  rocceig(nband_k,nband_k)$= (occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$,
!!   if this ratio has been attributed to the band n, 0.0_dp otherwise
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Supposing the occupations numbers differ :
!! if $abs(occ_{k,q}(m)) < abs(occ_k(n))$
!!  $rocceig(m,n)=(occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n)) $
!! if $abs(occ_{k,q}(m))>abs(occ_k(n))$
!!  rocceig(m,n)=0.0_dp
!!
!! If the occupation numbers are close enough, then
!! if the eigenvalues are also close, take the derivative
!!  $ rocceig(m,n)=\frac{1}{2}*docc/deig0 $
!! otherwise,
!!  $ rocceig(m,n)=\frac{1}{2}*(occ_{k,q}(m)-occ_k(n))/(eig0_{k,q}(m)-eig0_k(n))$
!!
!! PARENTS
!!      dfpt_nstpaw,dfpt_rhofermi,dfpt_vtorho
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,occopt,occ_k,occ_kq,rocceig)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'occeig'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband_k,occopt
!arrays
 real(dp),intent(in) :: doccde_k(nband_k),doccde_kq(nband_k),eig0_k(nband_k)
 real(dp),intent(in) :: eig0_kq(nband_k),occ_k(nband_k),occ_kq(nband_k)
 real(dp),intent(out) :: rocceig(nband_k,nband_k)

!Local variables-------------------------------
!scalars
 integer :: ibandk,ibandkq
 real(dp) :: diffabsocc,diffeig,diffocc,ratio,sumabsocc
 character(len=500) :: message

! *************************************************************************

!The parameter tol5 defines the treshhold for degeneracy, and the width of the step function

 rocceig(:,:)=0.0_dp

 do ibandk=1,nband_k
   do ibandkq=1,nband_k
     diffeig=eig0_kq(ibandkq)-eig0_k(ibandk)
     diffocc=occ_kq(ibandkq)-occ_k(ibandk)

     if( abs(diffeig) > tol5 ) then
       ratio=diffocc/diffeig
     else
       if(occopt<3)then
!        In a non-metallic case, if the eigenvalues are degenerate,
!        the occupation numbers must also be degenerate, in which
!        case there is no contribution from this pair of bands
         if( abs(diffocc) > tol5 ) then
           write(message,'(a,a,a,a,a,a,a,2(a,i4,a,es16.6,a,es16.6,a,a),a)' ) &
&           'In a non-metallic case (occopt<3), for a RF calculation,',ch10,&
&           'if the eigenvalues are degenerate,',&
&           ' the occupation numbers must also be degenerate.',ch10,&
&           'However, the following pair of states gave :',ch10,&
&           'k -state, band number',ibandk,', occ=',occ_k(ibandk),&
&           'eigenvalue=',eig0_k(ibandk),',',ch10,&
&           ' kq-state, band number',ibandkq,', occ=',occ_kq(ibandkq),&
&           ', eigenvalue=',eig0_kq(ibandkq),'.',ch10,&
&           'Action : change occopt, consistently, in GS and RF calculations.'
           MSG_ERROR(message)
         end if
         ratio=0.0_dp
       else
!        In the metallic case, one can compute a better approximation of the
!        ratio by using derivatives doccde
         ratio=0.5_dp*(doccde_kq(ibandkq)+doccde_k(ibandk))
!        DEBUG
!        write(std_out,*)' occeig : ibandkq,doccde_kq(ibandkq)',ibandkq,doccde_kq(ibandkq)
!        write(std_out,*)'          ibandk ,doccde_k (ibandk )',ibandk,doccde_k(ibandk)
!        ENDDEBUG
       end if
     end if

!    Here, must pay attention to the smallness of some coefficient
     diffabsocc=abs(occ_k(ibandk))-abs(occ_kq(ibandkq))
     sumabsocc=abs(occ_k(ibandk))+abs(occ_kq(ibandkq))
     if(sumabsocc>tol8)then
       if( diffabsocc > sumabsocc*tol5 ) then
         rocceig(ibandkq,ibandk)=ratio
       else if ( diffabsocc >= -sumabsocc*tol5 ) then
         rocceig(ibandkq,ibandk)=0.5_dp*ratio
       else
         rocceig(ibandkq,ibandk)=0.0_dp
       end if
     end if

   end do ! ibandkq
 end do ! ibandk

end subroutine occeig
!!***
