!{\src2tex{textfont=tt}}
!!****f* ABINIT/elt_ewald
!!
!! NAME
!! elt_ewald
!!
!! FUNCTION
!! Compute 2nd derivatives of Ewald energy wrt strain for frozen wavefunction
!! contributions to elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr^-1)
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! elteew(6+3*natom,6)=2nd derivatives of Ewald energy wrt strain
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine elt_ewald(elteew,gmet,gprimd,my_natom,natom,ntypat,rmet,rprimd,&
&                 typat,ucvol,xred,zion, &
&                 mpi_atmtab,comm_atom) ! optional arguments (parallelism)

 use defs_basis
 use m_profiling_abi
 use m_xmpi

 use m_special_funcs,  only : abi_derfc
 use m_paral_atom,     only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elt_ewald'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom),zion(ntypat)
 real(dp),intent(out) :: elteew(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer :: ia,ia0,ib,ierr,ig1,ig2,ig3,ir1,ir2,ir3,is1,is2,jj,js,ka,kb,kd,kg,my_comm_atom,newg,newr,ng,nr
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: arg,ch,chsq,cos_term,d2derfc,d2gss,d2r,d2rs,dderfc,derfc_arg
 real(dp) :: dgss1,dgss2,direct,dr1,dr2,drs1,drs2,eew,eta,fac,fraca1,fraca2
 real(dp) :: fraca3,fracb1,fracb2,fracb3,gsq,gsum,r1,r2,r3,recip,reta
 real(dp) :: rmagn,rsq,sin_term,sumg,summi,summr,sumr,t1,term
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer,pointer :: my_atmtab(:)
 real(dp) :: d2gm(3,3,6,6),d2ris(3),d2rm(3,3,6,6),dgm(3,3,6),dris(3),drm(3,3,6)
 real(dp) :: t2(3),ts2(3),tsec(2),tt(3)
 real(dp),allocatable :: d2sumg(:,:),d2sumr(:,:),drhoisi(:,:),drhoisr(:,:)
 real(dp),allocatable :: mpibuf(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' elt_ewald : enter '
!stop
!ENDDEBUG

!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     drm(:,jj,is2)=rprimd(kg,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kg,jj)
     dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
   end do

!  Loop over 1st strain index, upper triangle only
   do is1=1,is2

     ka=idx(2*is1-1);kb=idx(2*is1)
     d2rm(:,:,is1,is2)=zero
     d2gm(:,:,is1,is2)=zero
     do jj = 1,3
       if(ka==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kb,jj)
       if(ka==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(kb,jj)
       if(kb==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(ka,jj)
       if(kb==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(ka,jj)

       if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
       if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
       if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
       if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
     end do
   end do !is1
 end do !is2

!Add up total charge and sum of $charge^2$ in cell
 chsq=zero
 ch=zero
 do ia=1,natom
   ch=ch+zion(typat(ia))
   chsq=chsq+zion(typat(ia))**2
 end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!Here, a bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
 eta=pi*200._dp/33.0_dp*sqrt(1.69_dp*recip/direct)

!Conduct reciprocal space summations
 fac=pi**2/eta ; gsum=zero
 ABI_ALLOCATE(d2sumg,(6+3*natom,6))
 ABI_ALLOCATE(drhoisr,(3,natom))
 ABI_ALLOCATE(drhoisi,(3,natom))
 d2sumg(:,:)=zero

!Sum over G space, done shell after shell until all
!contributions are too small.
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          gsq is G dot G = |G|^2
           gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
&           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
&           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))

!          Skip g=0:
           if (gsq>1.0d-20) then
             arg=fac*gsq

!            Larger arg gives 0 contribution because of exp(-arg)
             if (arg <= 80._dp) then
!              When any term contributes then include next shell
               newg=1
               term=exp(-arg)/gsq
               summr = zero
               summi = zero
!              Note that if reduced atomic coordinates xred drift outside
!              of unit cell (outside [0,1)) it is irrelevant in the following
!              term, which only computes a phase.
               do ia=1,natom
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 cos_term=cos(arg)
                 sin_term=sin(arg)
                 summr=summr+zion(typat(ia))*cos_term
                 summi=summi+zion(typat(ia))*sin_term
                 drhoisr(1,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig1)
                 drhoisi(1,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig1)
                 drhoisr(2,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig2)
                 drhoisi(2,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig2)
                 drhoisr(3,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig3)
                 drhoisi(3,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig3)
               end do

!              The following two checks avoid an annoying
!              underflow error message
               if (abs(summr)<1.d-16) summr=zero
               if (abs(summi)<1.d-16) summi=zero

!              The product of term and summr**2 or summi**2 below
!              can underflow if not for checks above
               t1=term*(summr*summr+summi*summi)
               gsum=gsum+t1
!              Loop over 2nd strain index
               do is2=1,6
                 dgss2=dgm(1,1,is2)*dble(ig1*ig1)+dgm(2,2,is2)*dble(ig2*ig2)+&
&                 dgm(3,3,is2)*dble(ig3*ig3)+2._dp*(dgm(2,1,is2)*dble(ig1*ig2)+&
&                 dgm(3,1,is2)*dble(ig1*ig3)+dgm(3,2,is2)*dble(ig3*ig2))
!                Loop over 1st strain index, upper triangle only
                 do is1=1,is2
                   dgss1=dgm(1,1,is1)*dble(ig1*ig1)+dgm(2,2,is1)*dble(ig2*ig2)+&
&                   dgm(3,3,is1)*dble(ig3*ig3)+2._dp*(dgm(2,1,is1)*dble(ig1*ig2)+&
&                   dgm(3,1,is1)*dble(ig1*ig3)+dgm(3,2,is1)*dble(ig3*ig2))

                   d2gss=d2gm(1,1,is1,is2)*dble(ig1*ig1)+&
&                   d2gm(2,2,is1,is2)*dble(ig2*ig2)+&
&                   d2gm(3,3,is1,is2)*dble(ig3*ig3)+2._dp*&
&                   (d2gm(2,1,is1,is2)*dble(ig1*ig2)+&
&                   d2gm(3,1,is1,is2)*dble(ig1*ig3)+&
&                   d2gm(3,2,is1,is2)*dble(ig3*ig2))

                   d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*((fac**2 + 2.0_dp*fac/gsq + 2.0_dp/(gsq**2))*dgss1*dgss2 -&
&                   0.5_dp*(fac + 1.0_dp/gsq)*d2gss)
                   if(is1<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*(fac + 1.0_dp/gsq)*dgss2
                   if(is2<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*(fac + 1.0_dp/gsq)*dgss1
                   if(is1<=3 .and. is2<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+t1

                 end do !is1

!                Internal strain contributions
                 do ia=1,natom
                   js=7+3*(ia-1)
                   t2(:)=2.0_dp*term*(summr*drhoisr(:,ia)+summi*drhoisi(:,ia))
                   d2sumg(js:js+2,is2)=d2sumg(js:js+2,is2)-&
&                   (fac + 1.0_dp/gsq)*dgss2*t2(:)
                   if(is2<=3) d2sumg(js:js+2,is2)=d2sumg(js:js+2,is2)-t2(:)
                 end do
               end do !is2

             end if ! End condition of not larger than 80.0
           end if ! End skip g=0
         end if ! End triple loop over G s and associated new shell condition
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit
 end do !  End the loop on ng (new shells). Note that there is one exit from this loop.

 sumg=gsum/(two_pi*ucvol)
 d2sumg(:,:)=d2sumg(:,:)/(two_pi*ucvol)

 ABI_DEALLOCATE(drhoisr)
 ABI_DEALLOCATE(drhoisi)
!Stress tensor is now computed elsewhere (ewald2) hence do not need
!length scale gradients (used to compute them here).

!Conduct real space summations
 reta=sqrt(eta)
 fac=2._dp*sqrt(eta/pi)
 ABI_ALLOCATE(d2sumr,(6+3*natom,6))
 sumr=zero;d2sumr(:,:)=zero

!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0
!  
   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then

           do ia0=1,my_natom
             ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
             js=7+3*(ia-1)
!            Map reduced coordinate xred(mu,ia) into [0,1)
             fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
             fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
             fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=dble(ir1)+fracb1-fraca1
               r2=dble(ir2)+fracb2-fraca2
               r3=dble(ir3)+fracb3-fraca3
               rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
&               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Avoid zero denominators in 'term':
               if (rsq>=1.0d-24) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 term=zero
                 if (eta*rsq<64.0_dp) then
                   newr=1
                   rmagn=sqrt(rsq)
                   arg=reta*rmagn
!                  derfc computes the complementary error function
!                  dderfc is the derivative of the complementary error function
!                  d2derfc is the 2nd derivative of the complementary error function
                   dderfc=-fac*exp(-eta*rsq)
                   d2derfc=-2._dp*eta*rmagn*dderfc
                   derfc_arg = abi_derfc(arg)
                   term=derfc_arg/rmagn
                   sumr=sumr+zion(typat(ia))*zion(typat(ib))*term
                   tt(:)=rmet(:,1)*r1+rmet(:,2)*r2+rmet(:,3)*r3
                   dris(:)=tt(:)/rmagn
!                  Loop over 2nd strain index
                   do is2=1,6
                     drs2=drm(1,1,is2)*r1*r1+drm(2,2,is2)*r2*r2+&
&                     drm(3,3,is2)*r3*r3+&
&                     2.0_dp*(drm(2,1,is2)*r2*r1+drm(3,2,is2)*r3*r2+&
&                     drm(3,1,is2)*r1*r3)
                     dr2=0.5_dp*drs2/rmagn
!                    Loop over 1st strain index, upper triangle only
                     do is1=1,is2
                       drs1=drm(1,1,is1)*r1*r1+drm(2,2,is1)*r2*r2+&
&                       drm(3,3,is1)*r3*r3+&
&                       2.0_dp*(drm(2,1,is1)*r2*r1+drm(3,2,is1)*r3*r2+&
&                       drm(3,1,is1)*r1*r3)
                       dr1=0.5_dp*drs1/rmagn
                       d2rs=d2rm(1,1,is1,is2)*r1*r1+d2rm(2,2,is1,is2)*r2*r2+&
&                       d2rm(3,3,is1,is2)*r3*r3+&
&                       2.0_dp*(d2rm(2,1,is1,is2)*r2*r1+d2rm(3,2,is1,is2)*r3*r2+&
&                       d2rm(3,1,is1,is2)*r1*r3)
                       d2r=(0.25_dp*d2rs-dr1*dr2)/rmagn
                       d2sumr(is1,is2)=d2sumr(is1,is2)+&
&                       zion(typat(ia))*zion(typat(ib))*&
&                       ((d2derfc-2.0_dp*dderfc/rmagn+2.0_dp*derfc_arg/rsq)*dr1*dr2+&
&                       (dderfc-derfc_arg/rmagn)*d2r)/rmagn
                     end do !is1
!                    Internal strain contribution
                     ts2(:)=drm(:,1,is2)*r1+drm(:,2,is2)*r2+drm(:,3,is2)*r3
                     d2ris(:)=ts2(:)/rmagn-0.5_dp*drs2*tt(:)/(rsq*rmagn)

                     d2sumr(js:js+2,is2)=d2sumr(js:js+2,is2)-&
&                     2.0_dp*zion(typat(ia))*zion(typat(ib))*&
&                     ((d2derfc-2.0_dp*dderfc/rmagn+2.0_dp*derfc_arg/rsq)*dr2*dris(:)+&
&                     (dderfc-derfc_arg/rmagn)*d2ris(:))/rmagn
                   end do !is2
                 end if
               end if ! End avoid zero denominators in'term'

             end do ! end loop over ib:
           end do ! end loop over ia:
         end if ! end triple loop over real space points and associated condition of new shell
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newr==0) exit
 end do !  End loop on nr (new shells). Note that there is an exit within the loop

!In case of parallelism over atoms: communicate
 if (paral_atom) then
   call timab(48,1,tsec)
   ABI_ALLOCATE(mpibuf,((6+3*natom)*6+1)) 
   mpibuf(1:(6+3*natom)*6)=reshape(d2sumr(:,:),shape=(/((6+3*natom)*6)/))
   mpibuf((6+3*natom)*6+1)=sumr
   call xmpi_sum(mpibuf,my_comm_atom,ierr)
   sumr=mpibuf((6+3*natom)*6+1)
   d2sumr(:,:)=reshape(mpibuf(1:(6+3*natom)*6),shape=(/(6+3*natom),6/))
   ABI_DEALLOCATE(mpibuf) 
   call timab(48,2,tsec)
 end if

 sumr=0.5_dp*sumr
 d2sumr(:,:)=0.5_dp*d2sumr(:,:)
 fac=pi*ch**2/(2.0_dp*eta*ucvol)

!Finally assemble Ewald energy, eew
 eew=sumg+sumr-chsq*reta/sqrt(pi)-fac

 elteew(:,:)=d2sumg(:,:)+d2sumr(:,:)

!Additional term for all strains diagonal (from "fac" term in eew)
 elteew(1:3,1:3)=elteew(1:3,1:3)-fac

!Fill in lower triangle
 do is2=2,6
   do is1=1,is2-1
     elteew(is2,is1)=elteew(is1,is2)
   end do
 end do

 ABI_DEALLOCATE(d2sumg)
 ABI_DEALLOCATE(d2sumr)

!Output the final values of ng and nr
 write(message, '(a,i4,a,i4)' )' elt_ewald : nr and ng are ',nr,' and ',ng
 call wrtout(std_out,message,'COLL')

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine elt_ewald
!!***
