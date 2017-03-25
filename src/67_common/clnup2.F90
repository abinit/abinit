!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnup2
!! NAME
!! clnup2
!!
!!
!! FUNCTION
!! Perform more "cleanup" after completion of iterations.
!! This subroutine prints out more breakdown of force
!! information, shifts of atomic positions, and stresses.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fred(3,natom)=d(E_total)/d(xred) derivatives (hartree)
!!  grchempottn(3,natom)=d(E_chempot)/d(xred) derivatives (hartree)
!!  grewtn(3,natom)=d(E_Ewald)/d(xred) derivatives (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D2 dispersion (hartree)
!!  grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!  iscf=parameter controlling scf or non-scf iterations
!!  natom=number of atoms in unit cell
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  prtfor= >0 if forces have to be printed (0 otherwise)
!!  prtstr= >0 if stresses have to be printed (0 otherwise)
!!  prtvol=control print volume and debugging output
!!  start(3,natom)=starting coordinates in terms of real space
!!   primitive translations
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  synlgr(3,natom)=d(E_nlpsp)/d(xred) derivatives (hartree)
!!  xred(3,natom)=final coordinates in terms of primitive translations
!!
!! OUTPUT
!!  (only print)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine clnup2(n1xccc,fred,grchempottn,gresid,grewtn,grvdw,grxc,iscf,natom,ngrvdw,&
&                 prtfor,prtstr,prtvol,start,strten,synlgr,xred)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clnup2'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,n1xccc,natom,ngrvdw,prtfor,prtstr,prtvol
!arrays
 real(dp),intent(in) :: fred(3,natom),grchempottn(3,natom),gresid(3,natom)
 real(dp),intent(in) :: grewtn(3,natom),grvdw(3,ngrvdw)
 real(dp),intent(in) :: grxc(3,natom),start(3,natom),strten(6),synlgr(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
 character(len=*), parameter :: format01020 ="(i5,1x,3f20.12)"
!scalars
 integer :: iatom,mu
 real(dp) :: devsqr,grchempot2
 character(len=500) :: message

! *************************************************************************
!
!DEBUG
!write(std_out,*)' clnup2 : enter '
!ENDDEBUG

!Only print additional info for scf calculations
 if (iscf>=0) then

   if((prtvol>=10).and.(prtfor>0))then

     write(message, '(a,10x,a)' ) ch10,&
&     '===> extra information on forces <==='
     call wrtout(ab_out,message,'COLL')

     write(message, '(a)' ) ' ewald contribution to reduced grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(grewtn(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

     grchempot2=sum(grchempottn(:,:)**2)
     if(grchempot2>tol16)then
       write(message, '(a)' ) ' chemical potential contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grchempottn(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     endif

     write(message, '(a)' ) ' nonlocal contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(synlgr(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

     write(message, '(a)' ) ' local psp contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     if (n1xccc/=0) then
       do iatom=1,natom
         write(message,format01020) iatom,fred(:,iatom)-&
&         (grewtn(:,iatom)+grchempottn(:,iatom)+synlgr(:,iatom)+grxc(:,iatom)+gresid(:,iatom))
         call wrtout(ab_out,message,'COLL')
       end do
     else
       do iatom=1,natom
         write(message,format01020) iatom,fred(:,iatom)-&
&         (grewtn(:,iatom)+grchempottn(:,iatom)+synlgr(:,iatom)+gresid(:,iatom))
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     if (n1xccc/=0) then
       write(message, '(a)' ) ' core charge xc contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grxc(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     if (ngrvdw==natom) then
       write(message, '(a)' ) ' Van der Waals DFT-D contribution to reduced grads'
       call wrtout(ab_out,message,'COLL')
       do iatom=1,natom
         write(message,format01020) iatom,(grvdw(mu,iatom),mu=1,3)
         call wrtout(ab_out,message,'COLL')
       end do
     end if

     write(message, '(a)' ) ' residual contribution to red. grads'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message,format01020) iatom,(gresid(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do

   end if

!  Compute mean squared deviation from starting coords
   devsqr=0.0_dp
   do iatom=1,natom
     do mu=1,3
       devsqr=devsqr+(xred(mu,iatom)-start(mu,iatom))**2
     end do
   end do

!  When shift is nonnegligible then print values
   if (devsqr>1.d-14) then
     write(message, '(a,1p,e12.4,3x,a)' ) &
&     ' rms coord change=',sqrt(devsqr/dble(3*natom)),&
&     'atom, delta coord (reduced):'
     call wrtout(ab_out,message,'COLL')
     do iatom=1,natom
       write(message, '(1x,i5,2x,3f20.12)' ) iatom,&
&       (xred(mu,iatom)-start(mu,iatom),mu=1,3)
       call wrtout(ab_out,message,'COLL')
     end do
   end if

!  Write out stress results
   if (prtstr>0) then
     write(message, '(a,a)' ) ch10,&
&     ' Cartesian components of stress tensor (hartree/bohr^3)'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

!    Also output the pressure (minus one third the trace of the stress
!    tensor.
     write(message, '(a,a,es12.4,a)' ) ch10,&
&     '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&     -(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa]'

     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(1 1)=',strten(1)*HaBohr3_GPa,&
&     '  sigma(3 2)=',strten(4)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(2 2)=',strten(2)*HaBohr3_GPa,&
&     '  sigma(3 1)=',strten(5)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&     '- sigma(3 3)=',strten(3)*HaBohr3_GPa,&
&     '  sigma(2 1)=',strten(6)*HaBohr3_GPa
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

!  Last end if above refers to iscf > 0
 end if

!DEBUG
!write(std_out,*)' clnup2 : exit '
!ENDDEBUG

end subroutine clnup2
!!***
