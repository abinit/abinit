!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuenergy
!! NAME
!! pawuenergy
!!
!! FUNCTION
!! Compute contributions to energy for PAW+U calculations
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  iatom=index of current atom (note: this is the absolute index, not the index on current proc)
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!     %lpawu=l used for lda+u
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!  paw_ij <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!     %noccmmp(2*lpawu+1,2*lpawu+1,nspden)=density matrix in the PAW augm. region
!!     %nocctot(nspden)=number of electrons in the correlated subspace
!!  dmft_dc,e_ee,e_dc,e_dcdc,u_dmft,j_dmft= optional arguments for DMFT
!!
!! OUTPUT
!!  eldaumdc= PAW+U contribution to total energy
!!  eldaumdcdc= PAW+U contribution to double-counting total energy
!!
!! PARENTS
!!      m_energy,pawdenpot
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab,paw_ij,&
 &                     dmft_dc,e_ee,e_dc,e_dcdc,u_dmft,j_dmft) ! optional arguments (DMFT)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawtab, only : pawtab_type
 use m_paw_ij, only : paw_ij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawuenergy'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iatom,pawprtvol
 integer,optional,intent(in) :: dmft_dc
 real(dp),intent(inout) :: eldaumdc,eldaumdcdc
 real(dp),optional,intent(inout) :: e_ee,e_dc,e_dcdc
 real(dp),optional,intent(in) :: j_dmft,u_dmft
 type(paw_ij_type),intent(in) :: paw_ij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-U/4.N.(N-2)
!           2: E_int=-U/2.(Nup.(Nup-1)+Ndn.(Ndn-1))
 integer,parameter :: option_interaction=1

 integer :: cplex_dij,dmftdc,ispden,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41
 real(dp) :: eks_opt3,edcdc_opt3,edcdctemp,edctemp,eldautemp,jpawu,jpawu_dc,mnorm,mx,my,mz
 real(dp) :: n_sig,n_sigs,n_msig,n_msigs,n_dndn,n_tot,n_upup
 real(dp) :: n12_ud_im,n12_du_im
 real(dp) :: n12_ud_re,n12_du_re
 real(dp) :: n34_ud_im,n34_du_im
 real(dp) :: n34_ud_re,n34_du_re
 real(dp) :: upawu
 real(dp),allocatable :: n12_sig(:),n34_msig(:),n34_sig(:)
 character(len=500) :: message

! *****************************************************

 if(present(dmft_dc))  then
   dmftdc=dmft_dc
   if(pawtab%usepawu/=10) then
     write(message,'(4x,2a,i5)') "Error, usepawu should be equal to 10 if", &
&     " dmft_dc is an argument of pawuenergy",pawtab%usepawu
     call wrtout(std_out,message,'COLL')
   end if
 else
   dmftdc=pawtab%usepawu
 end if

 DBG_ENTER("COLL")

 lpawu=pawtab%lpawu
 upawu=pawtab%upawu;if(present(u_dmft)) upawu=u_dmft
 jpawu=pawtab%jpawu;if(present(j_dmft)) jpawu=j_dmft

!======================================================
!Compute LDA+U Energy
!-----------------------------------------------------

 eldautemp=zero
 edcdc_opt3=zero
 eks_opt3=zero
 cplex_dij=paw_ij%cplex_dij
 ABI_ALLOCATE(n12_sig,(cplex_dij))
 ABI_ALLOCATE(n34_msig,(cplex_dij))
 ABI_ALLOCATE(n34_sig,(cplex_dij))
 do ispden=1,min(paw_ij%ndij,2)
   jspden=min(paw_ij%ndij,2)-ispden+1
!  compute n_sigs and n_msigs for pawtab%usepawu=3
   if (paw_ij%nspden<=2) then
     n_sig =paw_ij%nocctot(ispden)
     n_msig=paw_ij%nocctot(jspden)
     n_tot=n_sig+n_msig
   else
     n_tot=paw_ij%nocctot(1)
     mx=paw_ij%nocctot(2)
     my=paw_ij%nocctot(3)
     mz=paw_ij%nocctot(4)
     mnorm=sqrt(mx*mx+my*my+mz*mz)
     if (ispden==1) then
!      n_sig =half*(n_tot+mnorm)
!      n_msig=half*(n_tot-mnorm)
       n_sig =half*(n_tot+sign(mnorm,mz))
       n_msig=half*(n_tot-sign(mnorm,mz))
     else
!      n_sig =half*(n_tot-mnorm)
!      n_msig=half*(n_tot+mnorm)
       n_sig =half*(n_tot-sign(mnorm,mz))
       n_msig=half*(n_tot+sign(mnorm,mz))
     end if
   end if
   n_sigs =n_sig/(float(2*lpawu+1))
   n_msigs =n_msig/(float(2*lpawu+1))
!  if(pawtab%usepawu==3) then
!  write(message,fmt=12) "noccmmp11 ",ispden,paw_ij%noccmmp(1,1,1,ispden)
!  call wrtout(std_out,message,'COLL')
!  write(message,fmt=12) "noccmmp11 ",jspden,paw_ij%noccmmp(1,1,1,jspden)
!  call wrtout(std_out,message,'COLL')
!  write(message,fmt=12) "n_sig      ",ispden,n_sig
!  call wrtout(std_out,message,'COLL')
!  write(message,fmt=12) "n_msig     ",jspden,n_msig
!  call wrtout(std_out,message,'COLL')
!  write(message,fmt=12) "n_sigs     ",ispden,n_sigs
!  call wrtout(std_out,message,'COLL')
!  write(message,fmt=12) "n_msigs    ",jspden,n_msigs
!  call wrtout(std_out,message,'COLL')
!  endif
!  12 format(a,i4,e20.10)
!  compute interaction energy E_{ee}
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       n12_sig(:)=paw_ij%noccmmp(:,m11,m21,ispden)
       if(m21==m11.and.(pawtab%usepawu==3.or.dmftdc==3)) n12_sig(1)=n12_sig(1)-n_sigs
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
           n34_sig(:) =paw_ij%noccmmp(:,m31,m41,ispden)
           n34_msig(:)=paw_ij%noccmmp(:,m31,m41,jspden)
           if(m31==m41.and.(pawtab%usepawu==3.or.dmftdc==3)) then
             n34_sig(1)= n34_sig(1) - n_sigs
             n34_msig(1)= n34_msig(1) - n_msigs
           end if
           eldautemp=eldautemp &
&           + n12_sig(1)*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&           + n12_sig(1)*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
           if(cplex_dij==2) then
             eldautemp=eldautemp &
&             - n12_sig(2)*n34_msig(2)*pawtab%vee(m11,m31,m21,m41) &
&             - n12_sig(2)*n34_sig(2) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
           end if
           if (pawtab%usepawu==3.or.dmftdc==3) then
             edcdc_opt3=edcdc_opt3 &
&             + n_sigs*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&             + n_sigs*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             eks_opt3=eks_opt3 &
&             + paw_ij%noccmmp(1,m11,m21,ispden)*n34_msig(1)*pawtab%vee(m11,m31,m21,m41) &
&             + paw_ij%noccmmp(1,m11,m21,ispden)*n34_sig(1) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             if(cplex_dij==2) then
               eks_opt3=eks_opt3 &
&               - paw_ij%noccmmp(2,m11,m21,ispden)*n34_msig(2)*pawtab%vee(m11,m31,m21,m41) &
&               - paw_ij%noccmmp(2,m11,m21,ispden)*n34_sig(2) *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
             end if
           end if
         end do ! m4
       end do ! m3
     end do ! m2
   end do ! m1
 end do ! ispden
 if (paw_ij%ndij==1) eldautemp=two*eldautemp ! Non-magn. system: sum up and dn energies
 ABI_DEALLOCATE(n12_sig)
 ABI_DEALLOCATE(n34_msig)
 ABI_DEALLOCATE(n34_sig)

!Non-collinear magnetism: add non-diagonal term; see (Eq 3) in PRB 72, 024458 (2005)
 if (paw_ij%ndij==4) then
   do m1=-lpawu,lpawu
     m11=m1+lpawu+1
     do m2=-lpawu,lpawu
       m21=m2+lpawu+1
       n12_ud_re=paw_ij%noccmmp(1,m11,m21,3) ! updn
       n12_ud_im=paw_ij%noccmmp(2,m11,m21,3) ! updn
       n12_du_re=paw_ij%noccmmp(1,m11,m21,4) ! dnup
       n12_du_im=paw_ij%noccmmp(2,m11,m21,4) ! dnup
       do m3=-lpawu,lpawu
         m31=m3+lpawu+1
         do m4=-lpawu,lpawu
           m41=m4+lpawu+1
           n34_ud_re=paw_ij%noccmmp(1,m31,m41,3)  ! updn
           n34_ud_im=paw_ij%noccmmp(2,m31,m41,3)  ! updn
           n34_du_re=paw_ij%noccmmp(1,m31,m41,4)  ! dnup
           n34_du_im=paw_ij%noccmmp(2,m31,m41,4)  ! dnup
           eldautemp=eldautemp-pawtab%vee(m11,m31,m41,m21) &
&           *(n12_ud_re*n34_du_re-n12_ud_im*n34_du_im &
&           +n12_du_re*n34_ud_re-n12_du_im*n34_ud_im)
           if (pawtab%usepawu==3.or.dmftdc==3) then
             eks_opt3=eks_opt3-pawtab%vee(m11,m31,m41,m21) &
&             *(n12_ud_re*n34_du_re-n12_ud_im*n34_du_im &
&             +n12_du_re*n34_ud_re-n12_du_im*n34_ud_im)
           end if
         end do ! m4
       end do ! m3
     end do ! m2
   end do ! m1
 end if

!Divide eldautemp by 2; see (Eq 1) in PRB 77, 155104 (2008)
 eldautemp=half*eldautemp

!if (paw_ij%ndij==1) then
!n_tot=two*paw_ij%nocctot(1)
!n_upup=paw_ij%nocctot(1)
!n_dndn=paw_ij%nocctot(1)
!else if (paw_ij%ndij==2) then
!n_tot=paw_ij%nocctot(1)+paw_ij%nocctot(2)
!n_upup=paw_ij%nocctot(1)
!n_dndn=paw_ij%nocctot(2)
!else if (paw_ij%ndij==4) then
!n_tot=paw_ij%nocctot(1)
!mx=paw_ij%nocctot(2)
!my=paw_ij%nocctot(3)
!mz=paw_ij%nocctot(4)
!mnorm=sqrt(mx*mx+my*my+mz*mz)
!n_upup=half*(n_tot+mnorm)
!n_dndn=half*(n_tot-mnorm)
!end if
 n_upup=n_sig
 n_dndn=n_msig

 edcdctemp=zero;edctemp=zero

!Full localized limit
 if(pawtab%usepawu==1.or.(dmftdc==1.or.dmftdc==4.or.dmftdc==5)) then
   jpawu_dc=jpawu
   if(dmftdc==4)  then
     jpawu_dc=zero
   end if
   edcdctemp=edcdctemp-half*upawu*n_tot**2
   edctemp  =edctemp  +half*upawu*(n_tot*(n_tot-one))
   if (paw_ij%ndij/=4.or.option_interaction==2) then
     if(dmftdc/=5) then
       edcdctemp=edcdctemp+half*jpawu_dc*(n_upup**2+n_dndn**2)
       edctemp  =edctemp  -half*jpawu_dc*(n_upup*(n_upup-one)+n_dndn*(n_dndn-one))
     else if(dmftdc==5)  then
       edcdctemp=edcdctemp+quarter*jpawu_dc*n_tot**2
       edctemp  =edctemp  -quarter*jpawu_dc*(n_tot*(n_tot-two))
     end if
   else if (paw_ij%ndij==4.and.option_interaction==1) then
!    write(message,'(a)') "  warning: option_interaction==1 for test         "
!    call wrtout(std_out,message,'COLL')
     edcdctemp=edcdctemp+quarter*jpawu_dc*n_tot**2
     edctemp  =edctemp  -quarter*jpawu_dc*(n_tot*(n_tot-two))
   else if (paw_ij%ndij==4.and.option_interaction==3) then
!    edcdctemp= \frac{J}/{4}[ N(N) + \vect{m}.\vect{m}]
     edcdctemp=edcdctemp+quarter*jpawu_dc*(n_tot**2 + &
&     mx**2+my**2+mz**2)  ! +\frac{J}/{4}\vect{m}.\vect{m}
!    edctemp= -\frac{J}/{4}[ N(N-2) + \vect{m}.\vect{m}]
     edctemp  =edctemp  -quarter*jpawu_dc*(  &
&     (n_tot*(n_tot-two)) +   &
&     mx**2+my**2+mz**2)  ! -\frac{J}/{4}\vect{m}.\vect{m}
   end if

!  Around mean field
 else if(pawtab%usepawu==2.or.dmftdc==2) then
   edctemp=edctemp+upawu*(n_upup*n_dndn)&
&   +half*(upawu-jpawu)*(n_upup**2+n_dndn**2) &
&   *(dble(2*lpawu)/dble(2*lpawu+1))
   edcdctemp=-edctemp
 else if(pawtab%usepawu==3.or.dmftdc==3) then
   edcdctemp=edcdc_opt3
   if(abs(pawprtvol)>=3) then
     write(message,fmt=11) "edcdc_opt3          ",edcdc_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "eks_opt3            ",eks_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "eks+edcdc_opt3      ",eks_opt3+edcdc_opt3
     call wrtout(std_out,message,'COLL')
     write(message,fmt=11) "(eks+edcdc_opt3)/2  ",(eks_opt3+edcdc_opt3)/2.d0
     call wrtout(std_out,message,'COLL')
   end if
 end if


 eldaumdc  =eldaumdc  +eldautemp-edctemp
 eldaumdcdc=eldaumdcdc-eldautemp-edcdctemp

!if(pawtab%usepawu/=10.or.pawprtvol>=3) then
 if(abs(pawprtvol)>=3) then
   if(pawtab%usepawu<10) then
     write(message, '(5a,i4)')ch10,'======= LDA+U Energy terms (in Hartree) ====',ch10,&
&     ch10,' For Atom ',iatom
   else if (pawtab%usepawu >= 10) then
     write(message, '(5a,i4)')ch10,'  ===   LDA+U Energy terms for the DMFT occupation matrix ==',ch10,&
&     ch10,' For Atom ',iatom
   end if

   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )"   Contributions to the direct expression of energy:"
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Double counting  correction   =",edctemp
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Interaction energy            =",eldautemp
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     Total LDA+U Contribution      =",eldautemp-edctemp
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' )' '
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' )"   For the ""Double-counting"" decomposition:"
   call wrtout(std_out,  message,'COLL')
   write(message,fmt=11) "     LDA+U Contribution            =",-eldautemp-edcdctemp
   call wrtout(std_out,  message,'COLL')
   11 format(a,e20.10)
   if(abs(pawprtvol)>=2) then
     write(message,fmt=11)"     edcdctemp                     =",edcdctemp
     call wrtout(std_out,  message,'COLL')
     write(message,fmt=11)"     eldaumdcdc for current atom   =",-eldautemp-edcdctemp
     call wrtout(std_out,  message,'COLL')
     write(message, '(a)' )' '
     call wrtout(std_out,  message,'COLL')
     write(message,fmt=11)"   pawuenergy: -VUKS pred          =",eldaumdcdc-eldaumdc
     call wrtout(std_out,  message,'COLL')
   end if
   write(message, '(a)' )' '
   call wrtout(std_out,  message,'COLL')
 end if

!For DMFT calculation
 if(present(e_ee))   e_ee=e_ee+eldautemp
 if(present(e_dc))   e_dc=e_dc+edctemp
 if(present(e_dcdc)) e_dcdc=e_dcdc+edcdctemp

 DBG_EXIT("COLL")

 end subroutine pawuenergy
!!***
