!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernl4b
!! NAME
!! opernl4b
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! from projected quantities to reciprocal space
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  if(choice==2 .or choice==4 .or. choice==5)
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(choice==3)
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  ------ Taken away in beautification because unused MS -------
!!   d2gxds2((2,mlang6,mincat,mproj) dummy argument here, not used
!!  -------------------------------------------------------------
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gxa(2,mlang3,mincat,mproj)= projected scalars
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the perturbation (needed if choice==2 or 5, and ndgxdt=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  ------ Taken away in beautification because unused MS -------
!!  istwf_k=option parameter that describes the storage of wfs
!!  -------------------------------------------------------------
!!  itypat = type of atom, needed for ffnl
!!  jproj(nlang)=number of projectors for each angular momentum
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  matblk=dimension of the array ph3d
!!  mincat= maximum increment of atoms
!!  mlang3 = one of the dimensions of the array gxa
!!  mlang4 = dimension for dgxds
!!  ------ Taken away in beautification because unused MS -------
!!  mlang6 = dimension for d2gxds2
!!  -------------------------------------------------------------
!!  mproj=maximum dimension for number of projection operators for each
!!    angular momentum for nonlocal pseudopotential
!!  ndgxdt=second dimension of dgxdt
!!  nincat = number of atoms in the subset here treated
!!  nkpg=second size of array kpg_k
!!  nffnl=third dimension of ffnl
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = number of plane waves in reciprocal space
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!  gxa(2,mlang3,nincat,mproj)=modified projected scalars;
!!  NOTE that metric contractions have already been performed
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!
!! OUTPUT
!!  vect(2*npw)=final vector in reciprocal space <G|V_nonlocal|vect_start>.
!!
!! NOTES
!! Operate with the non-local part of the hamiltonian for one type of
!! atom, and within this given type of atom, for a subset of
!! at most nincat atoms.
!!
!! This routine basically replaces getgla (gxa here is the former gla),
!! except for the calculation of <G|dVnl/dk|C> or strain gradients.
!!
!! Present version decomposed according to iffkg
!! opernl4a.f is from reciprocal space to projected quantities.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      dfpt_mkffkg
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine opernl4b(choice,dgxds,dgxdt,ffnl,gmet,gxa,&
&  ia3,idir,indlmn,ispinor,itypat,jproj,kg_k,kpg_k,kpt,&
&  lmnmax,matblk,mincat,mlang3,mlang4,mproj,ndgxdt,nffnl,nincat,&
&  nkpg,nlang,nloalg,npw,ntypat,ph3d,vect)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'opernl4b'
 use interfaces_66_nonlocal, except_this_one => opernl4b
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,ia3,idir,ispinor,itypat,lmnmax,matblk !,istwf_k
 integer,intent(in) :: mincat,mlang3,mlang4,mproj,ndgxdt,nffnl,nincat !,mlang6
 integer,intent(in) :: nkpg,nlang,npw,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),jproj(nlang),kg_k(3,npw)
 integer,intent(in) :: nloalg(3)
!real(dp),intent(in) :: d2gxds2(2,mlang6,mincat,mproj)
 real(dp),intent(in) :: dgxds(2,mlang4,mincat,mproj)
 real(dp),intent(in) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
 real(dp),intent(in) :: ffnl(1,npw,nffnl,lmnmax,ntypat),gmet(3,3)
 real(dp),intent(in) :: gxa(2,mlang3,mincat,mproj),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(in) :: ph3d(2,npw,matblk)
 real(dp),intent(inout) :: vect(2,npw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: chunk,ia,iaph3d,iffkg,iffkgk,iffkgs,ig,ii,ilang,ilang2,ilang3
 integer :: iproj,ipw,ipw1,ipw2,jffkg,jj,jump,mblkpw,nffkg
 integer :: nffkgd,nffkge,nffkgk,nffkgs,nincpw,nproj,ntens,start
 real(dp) :: ai,ar,sci1,sci2,sci3,sci4,sci5,sci6,sci7,sci8
 real(dp) :: scr1,scr2,scr3,scr4,scr5,scr6,scr7,scr8
 character(len=500) :: message
!arrays
 integer,allocatable :: parity(:)
 real(dp),allocatable :: ffkg(:,:),kpgx(:,:),scalars(:,:),teffv(:,:)

! *************************************************************************

!mblkpw sets the size of blocks of planewaves
 mblkpw=NLO_MBLKPW

!jump governs, in fine, the use of registers in the most cpu
!time consuming part of the routine. Until now, jump=8 is the maximal value.
!The optimal value will be machine-dependent !
 jump=4

!Initialisation before blocking on the plane waves

!Set up dimension of kpgx and allocate
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
 ntens=1
 if(nlang>=2 .or. choice==2 .or. choice==4 .or. choice==5)ntens=4
 if(nlang>=3 .or. choice==3)ntens=10
 if(nlang>=4 .or. (choice==3 .and. nlang>=2) )ntens=20
 if(choice==3 .and. nlang>=3)ntens=35
 if(choice==3 .and. nlang==4)ntens=56

!Set up second dimension of ffkg array, and allocate
 nffkg=0; nffkge=0; nffkgd=0; nffkgk=0; nffkgs=0
 do ilang=1,nlang
!  Get the number of projectors for that angular momentum
   nproj=jproj(ilang)
!  If there is a non-local part, accumulate the number of vectors needed
   if(nproj>0)then
     ilang2=(ilang*(ilang+1))/2
     nffkge=nffkge+nproj*ilang2
     if(choice==5)nffkgk=nffkgk+nproj*(2*ilang2-ilang)
     if(choice==2 .or. choice==4)nffkgd=nffkgd+ndgxdt*nproj*ilang2
     if(choice==3)then
       ilang3=((ilang+2)*(ilang+3))/2
       nffkgs=nffkgs+nproj*ilang3
     end if
   end if
 end do
 nffkg=nffkge+nffkgd+nffkgs+nffkgk

!!$OMP PARALLEL DEFAULT(PRIVATE) &
!!$OMP SHARED(choice,dgxds,dgxdt,ffnl,gmet,gxa,ia3,idir,indlmn,ispinor) &
!!$OMP SHARED(itypat,jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang3,mlang4,mproj) &
!!$OMP SHARED(ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,ntypat,ph3d,vect) &
!!$OMP SHARED(jump,nffkgd,nffkgk,nffkgs,mblkpw,nffkg,nffkge,ntens)

 ABI_ALLOCATE(ffkg,(nffkg,mblkpw))
 ABI_ALLOCATE(parity,(nffkg))
 ABI_ALLOCATE(kpgx,(mblkpw,ntens))
 ABI_ALLOCATE(scalars,(2,nffkg))
 ABI_ALLOCATE(teffv,(2,mblkpw))

!Loop on subsets of plane waves (blocking)
!!$OMP DO 
 do ipw1=1,npw,mblkpw

   ipw2=min(npw,ipw1+mblkpw-1)
   nincpw=ipw2-ipw1+1

!  Initialize kpgx array related to tensors defined below
   call dfpt_mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&
&   kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&   npw,ntens,ntypat,parity)

   if (choice==1 .or. choice==2 .or. choice==3 .or. choice==5) then
!    Application of non-local part from projected scalars
!    back to reciprocal space ...
!    [this section merely computes terms which add to <G|Vnl|C>;
!    nothing here is needed when various gradients are being computed]

!    Loop on atoms
     do ia=1,nincat

!      Compute the shift eventually needed to get the phases in ph3d
       iaph3d=ia
       if(nloalg(2)>0)iaph3d=ia+ia3-1

!      Transfer gxa (and eventually dgxdt) in scalars with different indexing
       iffkg=0
       iffkgk=nffkge*2
       iffkgs=nffkge
       do ilang=1,nlang
         nproj=jproj(ilang)
         if (nproj>0) then
           ilang2=(ilang*(ilang+1))/2
           ilang3=((ilang+2)*(ilang+3))/2
           do iproj=1,nproj
             do ii=1,ilang2
               jj=ii+((ilang-1)*ilang*(ilang+1))/6
               iffkg=iffkg+1
               if(choice==1 .or. choice==3)then
                 scalars(1,iffkg)=gxa(1,jj,ia,iproj)
                 scalars(2,iffkg)=gxa(2,jj,ia,iproj)
               else if (choice==2 .and. ndgxdt==1) then
                 jffkg=nffkge+iffkg
!                Pay attention to the use of reals and imaginary parts here ...
!                Also, the gxa and dgxdt arrays are switched, in order
!                to give the correct combination when multiplying ffkg,
!                see Eq.(53) of PRB55,10337(1997)
                 scalars(1,jffkg)= two_pi*gxa(2,jj,ia,iproj)
                 scalars(2,jffkg)=-two_pi*gxa(1,jj,ia,iproj)
                 scalars(1,iffkg)= dgxdt(1,1,jj,ia,iproj)
                 scalars(2,iffkg)= dgxdt(2,1,jj,ia,iproj)
               else if (choice==5) then
                 jffkg=nffkge+iffkg
!                The gxa and dgxdt arrays are switched, in order
!                to give the correct combination when multiplying ffkg,
                 scalars(1,jffkg)= gxa(1,jj,ia,iproj)
                 scalars(2,jffkg)= gxa(2,jj,ia,iproj)
                 scalars(1,iffkg)= dgxdt(1,1,jj,ia,iproj)
                 scalars(2,iffkg)= dgxdt(2,1,jj,ia,iproj)
               end if
             end do
             if(choice==3) then
               do ii=1,ilang3
                 iffkgs=iffkgs+1
                 jj=ii+((ilang+1)*(ilang+2)*(ilang+3))/6
                 scalars(1,iffkgs)=dgxds(1,jj-4,ia,iproj)
                 scalars(2,iffkgs)=dgxds(2,jj-4,ia,iproj)
               end do
             end if
             if(ilang>=2 .and. choice==5)then
               do ii=1,((ilang-1)*ilang)/2
                 jj=ii+((ilang-2)*(ilang-1)*ilang)/6
                 iffkgk=iffkgk+1
                 scalars(1,iffkgk)= dgxdt(1,2,jj,ia,iproj)
                 scalars(2,iffkgk)= dgxdt(2,2,jj,ia,iproj)
               end do
             end if
           end do
         end if
       end do

!      DEBUG
!      write(std_out,*)' opernl4b, write scalars '
!      do iffkg=1,nffkg
!      write(std_out,*)iffkg,scalars(1:2,iffkg)
!      end do
!      ENDDEBUG

!      ******* Entering the second critical part ****************************

!      First, treat small nffkg; send treat the loop needed for big nffkg;
!      finally treat the end of the loop needed for big nffkg

!      For the time being, the maximal jump allowed is 8.

!      1) Here, treat small nffkg
       if(nffkg<=jump)then

         select case(nffkg)

         case(1)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           ig=ipw1
           do ipw=1,nincpw
             ar=ffkg(1,ipw)*scr1 ; ai=ffkg(1,ipw)*sci1
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(2)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(3)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(4)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           scr4=scalars(1,4) ; sci4=scalars(2,4)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             ar=ar+ffkg(4,ipw)*scr4 ; ai=ai+ffkg(4,ipw)*sci4
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(5)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           scr4=scalars(1,4) ; sci4=scalars(2,4)
           scr5=scalars(1,5) ; sci5=scalars(2,5)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             ar=ar+ffkg(4,ipw)*scr4 ; ai=ai+ffkg(4,ipw)*sci4
             ar=ar+ffkg(5,ipw)*scr5 ; ai=ai+ffkg(5,ipw)*sci5
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(6)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           scr4=scalars(1,4) ; sci4=scalars(2,4)
           scr5=scalars(1,5) ; sci5=scalars(2,5)
           scr6=scalars(1,6) ; sci6=scalars(2,6)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             ar=ar+ffkg(4,ipw)*scr4 ; ai=ai+ffkg(4,ipw)*sci4
             ar=ar+ffkg(5,ipw)*scr5 ; ai=ai+ffkg(5,ipw)*sci5
             ar=ar+ffkg(6,ipw)*scr6 ; ai=ai+ffkg(6,ipw)*sci6
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(7)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           scr4=scalars(1,4) ; sci4=scalars(2,4)
           scr5=scalars(1,5) ; sci5=scalars(2,5)
           scr6=scalars(1,6) ; sci6=scalars(2,6)
           scr7=scalars(1,7) ; sci7=scalars(2,7)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             ar=ar+ffkg(4,ipw)*scr4 ; ai=ai+ffkg(4,ipw)*sci4
             ar=ar+ffkg(5,ipw)*scr5 ; ai=ai+ffkg(5,ipw)*sci5
             ar=ar+ffkg(6,ipw)*scr6 ; ai=ai+ffkg(6,ipw)*sci6
             ar=ar+ffkg(7,ipw)*scr7 ; ai=ai+ffkg(7,ipw)*sci7
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(8)

           scr1=scalars(1,1) ; sci1=scalars(2,1)
           scr2=scalars(1,2) ; sci2=scalars(2,2)
           scr3=scalars(1,3) ; sci3=scalars(2,3)
           scr4=scalars(1,4) ; sci4=scalars(2,4)
           scr5=scalars(1,5) ; sci5=scalars(2,5)
           scr6=scalars(1,6) ; sci6=scalars(2,6)
           scr7=scalars(1,7) ; sci7=scalars(2,7)
           scr8=scalars(1,8) ; sci8=scalars(2,8)
           ig=ipw1
           do ipw=1,nincpw
             ar=   ffkg(1,ipw)*scr1 ; ai=   ffkg(1,ipw)*sci1
             ar=ar+ffkg(2,ipw)*scr2 ; ai=ai+ffkg(2,ipw)*sci2
             ar=ar+ffkg(3,ipw)*scr3 ; ai=ai+ffkg(3,ipw)*sci3
             ar=ar+ffkg(4,ipw)*scr4 ; ai=ai+ffkg(4,ipw)*sci4
             ar=ar+ffkg(5,ipw)*scr5 ; ai=ai+ffkg(5,ipw)*sci5
             ar=ar+ffkg(6,ipw)*scr6 ; ai=ai+ffkg(6,ipw)*sci6
             ar=ar+ffkg(7,ipw)*scr7 ; ai=ai+ffkg(7,ipw)*sci7
             ar=ar+ffkg(8,ipw)*scr8 ; ai=ai+ffkg(8,ipw)*sci8
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do


         end select

       else

         do ipw=1,nincpw
           teffv(1,ipw)=0.0d0 ; teffv(2,ipw)=0.0d0
         end do

!        2) Here treart the loop for big nffkg
         do start=1,nffkg-jump,jump
           chunk=min(jump,nffkg-jump-start+1)

           select case(chunk)

           case(1)

             scr1=scalars(1,start) ; sci1=scalars(2,start)
             do ipw=1,nincpw
               ar=teffv(1,ipw)            ; ai=teffv(2,ipw)
               ar=ar+ffkg(start,ipw)*scr1 ; ai=ai+ffkg(start,ipw)*sci1
               teffv(1,ipw)=ar            ; teffv(2,ipw)=ai
             end do

           case(2)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(3)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(4)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(5)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
             scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
               ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(6)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
             scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
             scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
               ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
               ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(7)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
             scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
             scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
             scr7=scalars(1,start+6) ; sci7=scalars(2,start+6)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
               ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
               ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
               ar=ar+ffkg(start+6,ipw)*scr7 ; ai=ai+ffkg(start+6,ipw)*sci7
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           case(8)

             scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
             scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
             scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
             scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
             scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
             scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
             scr7=scalars(1,start+6) ; sci7=scalars(2,start+6)
             scr8=scalars(1,start+7) ; sci8=scalars(2,start+7)
             do ipw=1,nincpw
               ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
               ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
               ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
               ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
               ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
               ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
               ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
               ar=ar+ffkg(start+6,ipw)*scr7 ; ai=ai+ffkg(start+6,ipw)*sci7
               ar=ar+ffkg(start+7,ipw)*scr8 ; ai=ai+ffkg(start+7,ipw)*sci8
               teffv(1,ipw)=ar              ; teffv(2,ipw)=ai
             end do

           end select

         end do

!        3) Treat the end of the loops

         start=nffkg-jump+1

         select case(jump)

         case(1)

           scr1=scalars(1,start) ; sci1=scalars(2,start)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)            ; ai=teffv(2,ipw)
             ar=ar+ffkg(start,ipw)*scr1 ; ai=ai+ffkg(start,ipw)*sci1
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(2)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(3)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(4)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(5)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
           scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
             ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(6)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
           scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
           scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
             ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
             ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(7)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
           scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
           scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
           scr7=scalars(1,start+6) ; sci7=scalars(2,start+6)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
             ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
             ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
             ar=ar+ffkg(start+6,ipw)*scr7 ; ai=ai+ffkg(start+6,ipw)*sci7
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         case(8)

           scr1=scalars(1,start  ) ; sci1=scalars(2,start  )
           scr2=scalars(1,start+1) ; sci2=scalars(2,start+1)
           scr3=scalars(1,start+2) ; sci3=scalars(2,start+2)
           scr4=scalars(1,start+3) ; sci4=scalars(2,start+3)
           scr5=scalars(1,start+4) ; sci5=scalars(2,start+4)
           scr6=scalars(1,start+5) ; sci6=scalars(2,start+5)
           scr7=scalars(1,start+6) ; sci7=scalars(2,start+6)
           scr8=scalars(1,start+7) ; sci8=scalars(2,start+7)
           ig=ipw1
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             ar=ar+ffkg(start  ,ipw)*scr1 ; ai=ai+ffkg(start  ,ipw)*sci1
             ar=ar+ffkg(start+1,ipw)*scr2 ; ai=ai+ffkg(start+1,ipw)*sci2
             ar=ar+ffkg(start+2,ipw)*scr3 ; ai=ai+ffkg(start+2,ipw)*sci3
             ar=ar+ffkg(start+3,ipw)*scr4 ; ai=ai+ffkg(start+3,ipw)*sci4
             ar=ar+ffkg(start+4,ipw)*scr5 ; ai=ai+ffkg(start+4,ipw)*sci5
             ar=ar+ffkg(start+5,ipw)*scr6 ; ai=ai+ffkg(start+5,ipw)*sci6
             ar=ar+ffkg(start+6,ipw)*scr7 ; ai=ai+ffkg(start+6,ipw)*sci7
             ar=ar+ffkg(start+7,ipw)*scr8 ; ai=ai+ffkg(start+7,ipw)*sci8
             vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
             vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
             ig=ig+1
           end do

         end select

!        End if statement for small or big nffkg
       end if

!      ******* Leaving the critical part *********************************

!      End loop on atoms
     end do

!    End choice==1 or choice==2 or choice==3
   else
!    Problem: choice does not make sense
     write(message,'(a,i0,a)' )' Input choice=',choice,' not allowed. '
     MSG_BUG(message)
   end if

!  End loop on blocks of planewaves
 end do
!!$OMP END DO

 ABI_DEALLOCATE(ffkg)
 ABI_DEALLOCATE(kpgx)
 ABI_DEALLOCATE(parity)
 ABI_DEALLOCATE(scalars)
 ABI_DEALLOCATE(teffv)
!!$OMP END PARALLEL

end subroutine opernl4b
!!***
