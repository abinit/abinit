!{\src2tex{textfont=tt}}
!!****f* ABINIT/dielmt
!! NAME
!! dielmt
!!
!!
!! FUNCTION
!! Compute dielectric matrix from susceptibility matrix
!! Diagonalize it, then invert it.

!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  npwdiel=size of the dielinv and susmat arrays.
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  prtvol=control print volume and debugging output
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! OUTPUT
!!  dielinv(2,npwdiel,(nspden+4)/3,npwdiel,(nspden+4)/3)=inverse of the (non-hermitian)
!!      TC dielectric matrix in reciprocal space.
!!
!! NOTES
!! Warning : will not work in the spin-polarized, metallic case.
!! Output (not cleaned)
!! !!! Spin behaviour is not obvious !!!
!!
!! TODO
!! Write equation below (hermitian matrix)
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      timab,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dielmt(dielinv,gmet,kg_diel,&
&  npwdiel,nspden,occopt,prtvol,susmat)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dielmt'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwdiel,nspden,occopt,prtvol
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)

!Local variables-------------------------------
!scalars
 integer :: ieig,ier,ii,index,ipw,ipw1,ipw2,isp,jj,npwsp
 real(dp) :: ai1,ai2,ar1,ar2,eiginv,gfact,gfactinv,gred1,gred2,gred3,gsquar
 real(dp) :: tpisq
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dielh(:),dielmat(:,:,:,:,:),dielvec(:,:,:)
 real(dp),allocatable :: eig_diel(:),zhpev1(:,:),zhpev2(:)
!no_abirules
!integer :: ipw3
!real(dp) :: elementi,elementr

! *************************************************************************

!DEBUG
!write(std_out,*)' dielmt : enter '
!if(.true.)stop
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 call timab(90,1,tsec)

!-Compute now the hermitian dielectric matrix------------------------------
!Following remarks are only valid within RPA approximation (Kxc=0):

!for the spin-unpolarized case, 1 - 4pi (1/G) chi0(G,Gp) (1/Gp)

!for the spin-polarized case,
!( 1  0 ) - 4pi ( 1/G  1/G )   ( chi0 upup  chi0 updn )   ( 1/Gp 1/Gp )
!( 0  1 )       ( 1/G  1/G )   ( chi0 dnup  chi0 dndn )   ( 1/Gp 1/Gp )
!which is equal to
!( 1  0 ) - 4pi (1/G  0 ) (chi0 upup+dndn+updn+dnup  chi0 upup+dndn+updn+dnup) (1/Gp 0  )
!( 0  1 )       ( 0  1/G) (chi0 upup+dndn+updn+dnup  chi0 upup+dndn+updn+dnup) ( 0  1/Gp)
!So, if spin-polarized, sum all spin contributions
!Note: chi0 updn = chi0 dnup = zero for non-metallic systems

!In the case of non-collinear magnetism, within RPA, this is the same because:
!chi0_(s1,s2),(s3,s4) = delta_s1,s2 * delta_s3,s4 * chi0_(s1,s1),(s3,s3)
!Only chi_upup,upup, chi_dndn,dndn, chi_upup,dndn and chi_dndn,upup
!have to be taken into account (stored, susmat(:,ipw1,1:2,ipw2,1:2)

 ABI_ALLOCATE(dielmat,(2,npwdiel,min(nspden,2),npwdiel,min(nspden,2)))

 if(nspden/=1)then
   if (occopt<3) then
     do ipw2=1,npwdiel
       do ipw1=1,npwdiel
         dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)+susmat(1,ipw1,2,ipw2,2)
         dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw1,2,ipw2,2)
       end do
     end do
   else
     do ipw2=1,npwdiel
       do ipw1=1,npwdiel
         dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)+susmat(1,ipw1,2,ipw2,2)+susmat(1,ipw1,1,ipw2,2)+susmat(1,ipw1,2,ipw2,1)
         dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw1,2,ipw2,2)+susmat(2,ipw1,1,ipw2,2)+susmat(2,ipw1,2,ipw2,1)
       end do
     end do
   end if
 else
   do ipw2=1,npwdiel
     do ipw1=1,npwdiel
       dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)
     end do
   end do
 end if
!Compute 1/G factors and include them in the dielectric matrix
 do ipw1=1,npwdiel
   gred1=dble(kg_diel(1,ipw1))
   gred2=dble(kg_diel(2,ipw1))
   gred3=dble(kg_diel(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +two*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3)                        )
!  Distinguish G=0 from other elements
   if(gsquar>tol12)then
!    !$ gfact=\sqrt (4.0_dp \pi/gsquar/dble(nspden))$
     gfact=sqrt(four_pi/gsquar)
     do ipw2=1,npwdiel
!      Must multiply both rows and columns, and also changes the sign
       dielmat(1,ipw2,1,ipw1,1)=-dielmat(1,ipw2,1,ipw1,1)*gfact
       dielmat(2,ipw2,1,ipw1,1)=-dielmat(2,ipw2,1,ipw1,1)*gfact
       dielmat(1,ipw1,1,ipw2,1)= dielmat(1,ipw1,1,ipw2,1)*gfact
       dielmat(2,ipw1,1,ipw2,1)= dielmat(2,ipw1,1,ipw2,1)*gfact
     end do
   else
!    Zero the G=0 elements, head and wings
     do ipw2=1,npwdiel
       dielmat(1,ipw2,1,ipw1,1)=zero
       dielmat(2,ipw2,1,ipw1,1)=zero
       dielmat(1,ipw1,1,ipw2,1)=zero
       dielmat(2,ipw1,1,ipw2,1)=zero
     end do
   end if
 end do

!Complete the matrix in the spin-polarized case
!should this be nspden==2??
 if(nspden/=1)then
   do ipw1=1,npwdiel
     do ipw2=1,npwdiel
       dielmat(1,ipw1,1,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,1,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
       dielmat(1,ipw1,2,ipw2,1)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,2,ipw2,1)=dielmat(2,ipw1,1,ipw2,1)
       dielmat(1,ipw1,2,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,2,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
     end do
   end do
 end if

!DEBUG
!write(std_out,*)' dielmt : make dielmat equal to identity matrix '
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!dielmat(1,ipw1,1,ipw2,1)=0.0_dp
!dielmat(2,ipw1,1,ipw2,1)=0.0_dp
!end do
!end do
!ENDDEBUG

!Add the diagonal part
 do isp=1,min(nspden,2)
   do ipw=1,npwdiel
     dielmat(1,ipw,isp,ipw,isp)=one+dielmat(1,ipw,isp,ipw,isp)
   end do
 end do

!-The hermitian dielectric matrix is computed ------------------------------
!-Now, diagonalize it ------------------------------------------------------

!In RPA, everything is projected on the spin-symmetrized
!space. This was coded here (for the time being).

!Diagonalize the hermitian dielectric matrix

!npwsp=npwdiel*nspden
 npwsp=npwdiel

 ABI_ALLOCATE(dielh,(npwsp*(npwsp+1)))
 ABI_ALLOCATE(dielvec,(2,npwsp,npwsp))
 ABI_ALLOCATE(eig_diel,(npwsp))
 ABI_ALLOCATE(zhpev1,(2,2*npwsp-1))
 ABI_ALLOCATE(zhpev2,(3*npwsp-2))
 ier=0
!Store the dielectric matrix in proper mode before calling zhpev
 index=1
 do ii=1,npwdiel
   do jj=1,ii
     dielh(index  )=dielmat(1,jj,1,ii,1)
     dielh(index+1)=dielmat(2,jj,1,ii,1)
     index=index+2
   end do
 end do
!If spin-polarized and non RPA, need to store other parts of the matrix
!if(nspden/=1)then
!do ii=1,npwdiel
!Here, spin-flip contribution
!do jj=1,npwdiel
!dielh(index  )=dielmat(1,jj,1,ii,2)
!dielh(index+1)=dielmat(2,jj,1,ii,2)
!index=index+2
!end do
!Here spin down-spin down upper matrix
!do jj=1,ii
!dielh(index  )=dielmat(1,jj,2,ii,2)
!dielh(index+1)=dielmat(2,jj,2,ii,2)
!index=index+2
!end do
!end do
!end if

 call ZHPEV ('V','U',npwsp,dielh,eig_diel,dielvec,npwdiel,zhpev1,&
& zhpev2,ier)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

 if(prtvol>=10)then
   write(message, '(a,a,a,5es12.4)' )ch10,&
&   ' Five largest eigenvalues of the hermitian RPA dielectric matrix:',&
&   ch10,eig_diel(npwdiel:npwdiel-4:-1)
   call wrtout(ab_out,message,'COLL')
 end if

 write(message, '(a,a)' )ch10,&
& ' dielmt : 15 largest eigenvalues of the hermitian RPA dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(npwdiel:npwdiel-4:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  6-10 :',eig_diel(npwdiel-5:npwdiel-9:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  11-15:',eig_diel(npwdiel-10:npwdiel-14:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,a)' )ch10,&
& ' dielmt : 5 smallest eigenvalues of the hermitian RPA dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(1:5)
 call wrtout(std_out,message,'COLL')

!Invert the hermitian dielectric matrix,
!Should use a BLAS call !
 do ipw2=1,npwdiel
   do ipw1=ipw2,npwdiel
     dielinv(1,ipw1,1,ipw2,1)=zero  
     dielinv(2,ipw1,1,ipw2,1)=zero  
   end do
 end do
 do ieig=1,npwdiel
   eiginv=one/eig_diel(ieig)
   do ipw2=1,npwdiel
     do ipw1=ipw2,npwdiel
       ar1=dielvec(1,ipw1,ieig)
       ai1=dielvec(2,ipw1,ieig)
       ar2=dielvec(1,ipw2,ieig)
       ai2=dielvec(2,ipw2,ieig)
       dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)+&
&       (ar1*ar2+ai1*ai2)*eiginv
       dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)+&
&       (ai1*ar2-ar1*ai2)*eiginv
     end do
   end do
 end do
 do ipw2=1,npwdiel-1
   do ipw1=ipw2+1,npwdiel
     dielinv(1,ipw2,1,ipw1,1)= dielinv(1,ipw1,1,ipw2,1)
     dielinv(2,ipw2,1,ipw1,1)=-dielinv(2,ipw1,1,ipw2,1)
   end do
 end do

 ABI_DEALLOCATE(dielh)
 ABI_DEALLOCATE(dielvec)
 ABI_DEALLOCATE(eig_diel)

!DEBUG
!Checks whether the inverse of the hermitian dielectric matrix
!has been correctly generated
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!elementr=0.0_dp
!elementi=0.0_dp
!do ipw3=1,npwdiel
!elementr=elementr+dielinv(1,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)&
!&                    -dielinv(2,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)
!elementi=elementi+dielinv(1,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)&
!&                    +dielinv(2,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)
!end do
!if(elementr**2+elementi**2 > 1.0d-12)then
!if( ipw1 /= ipw2 .or. &
!&        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!write(std_out,*)' dielmt : the inversion procedure is not correct '
!write(std_out,*)' ipw1, ipw2 =',ipw1,ipw2
!write(std_out,*)' elementr,elementi=',elementr,elementi
!stop
!end if
!end if
!end do
!end do
!write(std_out,*)'dielmt : matrix has been inverted successfully '
!stop
!ENDDEBUG

!Then get the inverse of the asymmetric
!dielectric matrix, as required for the preconditioning.

!Inverse of the dielectric matrix : ( 1 - 4pi (1/G^2) chi0(G,Gp) )^(-1)
!In dielinv there is now (1 - 4pi (1/G) chi0(G,Gp) (1/Gp) )^(-1)
!So, evaluate dielinv_after(G,Gp) =
!(4pi/G^2)^(1/2) dielinv_before(G,Gp) (4pi/Gp^2)^(-1/2)
!In RPA, can focus on the spin-averaged quantities
 do ipw1=1,npwdiel
   gred1=dble(kg_diel(1,ipw1))
   gred2=dble(kg_diel(2,ipw1))
   gred3=dble(kg_diel(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +two*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3)                        )
!  Distinguish G=0 from other elements
   if(gsquar>tol12)then
     gfact=sqrt(four_pi/gsquar)
     gfactinv=one/gfact
     do ipw2=1,npwdiel
!      Must multiply both rows and columns
       dielinv(1,ipw2,1,ipw1,1)=dielinv(1,ipw2,1,ipw1,1)*gfactinv
       dielinv(2,ipw2,1,ipw1,1)=dielinv(2,ipw2,1,ipw1,1)*gfactinv
       dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)*gfact
       dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)*gfact
     end do
   else
!    Zero the G=0 elements, head
     do ipw2=1,npwdiel
       if (ipw2/=ipw1) dielinv(1:2,ipw1,1,ipw2,1)=zero
     end do
   end if
 end do

 ABI_DEALLOCATE(dielmat)

 call timab(90,2,tsec)

end subroutine dielmt
!!***
