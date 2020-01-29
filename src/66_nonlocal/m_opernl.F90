!!****m* ABINIT/m_opernl
!! NAME
!!  m_opernl
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernl

 use defs_basis
 use m_errors
 use m_abicore

 use m_mkffkg, only : mkffkg, dfpt_mkffkg

 implicit none

 private
!!***

 public :: opernl2
 public :: opernl3
 public :: opernl4a
 public :: opernl4b
!!***

contains
!!***

!!****f* ABINIT/opernl2
!! NAME
!! opernl2
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! either from reciprocal space to projected quantities (sign=1),
!! or from projected quantities to reciprocal space (sign=-1)
!!
!! INPUTS
!!  if(sign==-1 .and. (choice==2 .or choice==4 .or. choice==5))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==-1 .and. choice==3)
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the perturbation (needed if choice==2 or 5, and ndgxdt=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  istwf_k=option parameter that describes the storage of wfs
!!  itypat = type of atom, needed for ffnl
!!  jproj(nlang)=number of projectors for each angular momentum
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  matblk=dimension of the array ph3d
!!  mincat= maximum increment of atoms
!!  mlang1 = dimensions for dgxdis1
!!  mlang3 = one of the dimensions of the array gxa
!!  mlang4 = dimension for dgxds
!!  mlang5 = dimensions for dgxdis2
!!  mlang6 = dimension for d2gxds2
!!  mproj=maximum dimension for number of projection operators for each
!!    angular momentum for nonlocal pseudopotential
!!  ndgxdt=second dimension of dgxdt
!!  nffnl=second dimension of ffnl
!!  nincat = number of atoms in the subset here treated
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = number of plane waves in reciprocal space
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!  sign : if  1, go from reciprocal space to projected scalars,
!!         if -1, go from projected scalars to reciprocal space.
!!  if(sign==1), vect(2*npw)=starting vector in reciprocal space
!!  if(sign==-1) gxa(2,mlang3,nincat,mproj)=modified projected scalars;
!!   NOTE that metric contractions have already been performed on the
!!   arrays gxa if sign=-1
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!
!! OUTPUT
!!  if(sign==1)
!!   gxa(2,mlang3,mincat,mproj)= projected scalars
!!  if(sign==1 .and. (choice==2 .or choice==4 .or. choice==5 .or. choice==23))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==1 .and. (choice==3 .or. choice==23))
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  if(sign==1 .and. choice==6)
!!   dgxdis((2,mlang1,mincat,mproj) = derivatives of projected scalars
!!    wrt coord. indexed for internal strain
!!   d2gxdis((2,mlang5,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strain and coord
!!   d2gxds2((2,mlang6,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strains
!!  if(sign==-1)
!!   vect(2*npw)=final vector in reciprocal space <G|V_nonlocal|vect_start>.
!!
!! NOTES
!! Operate with the non-local part of the hamiltonian for one type of
!! atom, and within this given type of atom, for a subset
!! of at most nincat atoms.
!!
!! This routine basically replaces getgla (gxa here is the former gla),
!! except for the calculation of <G|dVnl/dk|C> or strain gradients.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      mkffkg
!!
!! SOURCE

subroutine opernl2(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&
&  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&
&  ntypat,ph3d,sign,vect)


!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,ia3,idir,ispinor,istwf_k,itypat,lmnmax,matblk
 integer,intent(in) :: mincat,mlang1,mlang3,mlang4,mlang5,mlang6,mproj,ndgxdt
 integer,intent(in) :: nffnl,nincat,nkpg,nlang,npw,ntypat,sign
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),jproj(nlang),kg_k(3,npw)
 integer,intent(in) :: nloalg(3)
 real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat),gmet(3,3),kpg_k(npw,nkpg)
 real(dp),intent(in) :: kpt(3),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj)
 real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
 real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj),vect(:,:)
 real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
 real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
 real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)

!Local variables-------------------------------
!scalars
 integer :: ia,iaph3d,iffkg,iffkgk,iffkgs,iffkgs2,ig,ii,ilang,ilang2,ilang3
 integer :: ilang4,ilang5,ilang6,ilangx,iproj,ipw,ipw1,ipw2,jffkg,jj,jjs,mblkpw
 integer :: mmproj,mu,nffkg,nffkgd,nffkge,nffkgk,nffkgs,nffkgs2,nincpw,nproj,ntens
 real(dp),parameter :: two_pi2=two_pi**2
 character(len=500) :: message
!arrays
 integer,allocatable :: parity(:)
!real(dp) :: tsec(2)
 real(dp),allocatable :: ffkg(:,:),kpgx(:,:),scalars(:,:),teffv(:,:)

! *************************************************************************

!call wrtout(std_out,"in opernl2","COLL")

!mblkpw sets the size of blocks of planewaves
 mblkpw=NLO_MBLKPW

!Get the actual maximum number of projectors
 mmproj=maxval(indlmn(3,:,itypat))

!Initialisation before blocking on the plane waves

 if (sign==1) then
!  Put projected scalars to zero
   gxa(:,:,:,1:mmproj)=0.0d0
   if (choice==2 .or. choice==4 .or. choice==5 .or. choice==23) dgxdt(:,:,:,:,1:mmproj)=0.0d0
   if (choice==3 .or. choice==6 .or. choice==23) dgxds(:,:,:,1:mmproj)=0.0d0
   if (choice==6) then
     dgxdis(:,:,:,1:mmproj)=0.0d0
     d2gxdis(:,:,:,1:mmproj)=0.0d0
     d2gxds2(:,:,:,1:mmproj)=0.0d0
   end if
 end if

!Set up dimension of kpgx and allocate
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
 ntens=1
 if(nlang>=2 .or. choice==2 .or. choice==4 .or. choice==5 .or. choice==23) ntens=4
 if(nlang>=3 .or. (choice==3.or.choice==23))ntens=10
 if(nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) )ntens=20
 if(((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6)ntens=35
 if(((choice==3.or.choice==23) .and. nlang==4) .or. (choice==6 .and. nlang>=2))ntens=56
 if(choice==6 .and. nlang>=3)ntens=84
 if(choice==6 .and. nlang==4)ntens=120

!Set up second dimension of ffkg array, and allocate
 nffkg=0 ; nffkge=0 ; nffkgd=0 ; nffkgk=0 ; nffkgs=0 ; nffkgs2=0
 do ilang=1,nlang
!  Get the number of projectors for that angular momentum
   nproj=jproj(ilang)
!  If there is a non-local part, accumulate the number of vectors needed
!  The variables ilang below are the number of independent tensors of
!  various ranks, the variable names being more historical than logical.
!  ilang2=number of rank ilang-1
!  ilang3=number of rank ilang+1
!  ilang4=number of rank ilang
!  ilang5=number of rank ilang+2
!  ilang6=number of rank ilang+3
   if(nproj>0)then
     ilang2=(ilang*(ilang+1))/2
     nffkge=nffkge+nproj*ilang2
     if(choice==5)nffkgk=nffkgk+nproj*(2*ilang2-ilang)
     if(choice==2 .or. choice==4 .or. choice==23)nffkgd=nffkgd+ndgxdt*nproj*ilang2
     if(choice==3 .or. choice==6 .or. choice==23)then
       ilang3=((ilang+2)*(ilang+3))/2
       nffkgs=nffkgs+nproj*ilang3
     end if
     if(choice==6)then
       ilang4=((ilang+1)*(ilang+2))/2
       ilang5=((ilang+3)*(ilang+4))/2
       ilang6=((ilang+4)*(ilang+5))/2
       nffkgs2=nffkgs2+nproj*(ilang4+ilang5+ilang6)
     end if
   end if
 end do
 nffkg=nffkge+nffkgd+nffkgs+nffkgs2+nffkgk

!Loop on subsets of plane waves (blocking)
!Disabled by MG on Dec  6 2011, omp sections have to be tested, this coding causes a sigfault with nthreads==1
!Feb 16 2012: The code does not crash anymore but it's not efficient.
!
!!$OMP PARALLEL DEFAULT(PRIVATE) &
!!$OMP SHARED(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt) &
!!$OMP SHARED(ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat) &
!!$OMP SHARED(jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4) &
!!$OMP SHARED(mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw) &
!!$OMP SHARED(ntypat,ph3d,sign,vect) &
!!$OMP SHARED(mblkpw,nffkg,nffkgd,nffkge,nffkgs,ntens,mmproj)

 ABI_ALLOCATE(ffkg,(mblkpw,nffkg))
 ABI_ALLOCATE(parity,(nffkg))
 ABI_ALLOCATE(kpgx,(mblkpw,ntens))
 ABI_ALLOCATE(scalars,(2,nffkg))
 ABI_ALLOCATE(teffv,(2,mblkpw))
!!$OMP DO
 do ipw1=1,npw,mblkpw

   ipw2=min(npw,ipw1+mblkpw-1)
   nincpw=ipw2-ipw1+1

!  call timab(74+choice,1,tsec)

!  Initialize kpgx array related to tensors defined below
   call mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,kg_k,&
&   kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&   npw,ntens,ntypat,parity)

!  call timab(74+choice,2,tsec)

!  Now treat the different signs
   if (sign==1) then

     do ia=1,nincat

!      Compute the shift eventually needed to get the phases in ph3d
       iaph3d=ia
       if(nloalg(2)>0)iaph3d=ia+ia3-1

!      ******* Entering the first time-consuming part of the routine *******

!      Multiply by the phase factor
!      This allows to be left with only real operations,
!      that are performed in the most inner loops
       ig=ipw1
       do ipw=1,nincpw
         teffv(1,ipw)=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
         teffv(2,ipw)=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
         ig=ig+1
       end do

       do iffkg=1,nffkg
         scalars(1,iffkg)=0.0d0
         scalars(2,iffkg)=0.0d0
         do ipw=1,nincpw
           scalars(1,iffkg)=scalars(1,iffkg)+teffv(1,ipw)*ffkg(ipw,iffkg)
           scalars(2,iffkg)=scalars(2,iffkg)+teffv(2,ipw)*ffkg(ipw,iffkg)
         end do
       end do

!      ******* Leaving the critical part *********************************

!      DEBUG
!      write(std_out,*)' opernl2 : scalars'
!      do iffkg=1,10
!      write(std_out,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!      end do
!      stop
!      ENDDEBUG

       if(istwf_k>=2)then
!        Impose parity of resulting scalar (this operation could be
!        replaced by direct saving of CPU time in the preceeding section)
         do iffkg=1,nffkg
           scalars(parity(iffkg),iffkg)=0.0d0
         end do
       end if

       iffkg=0
       iffkgs=nffkge+nffkgd
       iffkgs2=nffkge+nffkgs
       iffkgk=nffkge*2
       do ilang=1,nlang
         nproj=jproj(ilang)
         if(nproj>0)then
!          ilang2 is the number of independent tensor components
!          for symmetric tensor of rank ilang-1
           ilang2=(ilang*(ilang+1))/2

!          Loop over projectors
           do iproj=1,nproj
!            Multiply by the k+G factors (tensors of various rank)
             do ii=1,ilang2
!              Get the starting address for the relevant tensor
               jj=ii+((ilang-1)*ilang*(ilang+1))/6
               iffkg=iffkg+1
!              !$OMP CRITICAL (OPERNL2_1)
               gxa(1,jj,ia,iproj)=gxa(1,jj,ia,iproj)+scalars(1,iffkg)
               gxa(2,jj,ia,iproj)=gxa(2,jj,ia,iproj)+scalars(2,iffkg)
!              !$OMP END CRITICAL (OPERNL2_1)
!              Now, compute gradients, if needed.
               if ((choice==2.or.choice==23) .and. ndgxdt==3) then
                 do mu=1,3
                   jffkg=nffkge+(iffkg-1)*3+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  !$OMP CRITICAL (OPERNL2_2)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                  !$OMP END CRITICAL (OPERNL2_2)
                 end do
               end if
               if (choice==2 .and. ndgxdt==1) then
                 jffkg=nffkge+iffkg
!                Pay attention to the use of reals and imaginary parts here ...
!                !$OMP CRITICAL (OPERNL2_3)
                 dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                 dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                !$OMP END CRITICAL (OPERNL2_3)
               end if
               if (choice==4) then
                 do mu=1,3
                   jffkg=nffkge+(iffkg-1)*9+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  !$OMP CRITICAL (OPERNL2_4)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                  !$OMP END CRITICAL (OPERNL2_4)
                 end do
                 do mu=4,9
                   jffkg=nffkge+(iffkg-1)*9+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  Also, note the multiplication by (2 pi)**2
!                  !$OMP CRITICAL (OPERNL2_5)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi2*scalars(1,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)-two_pi2*scalars(2,jffkg)
!                  !$OMP END CRITICAL (OPERNL2_5)
                 end do
               end if
!              End loop on ii=1,ilang2
             end do

             if (choice==3 .or. choice==6 .or. choice==23) then
!              Compute additional tensors related to strain gradients
!              ilang3 is number of unique tensor components of rank ilang+1
               ilang3=((ilang+2)*(ilang+3))/2
               jjs=((ilang+1)*(ilang+2)*(ilang+3))/6
!              Compute strain gradient tensor components
               do ii=1,ilang3
!                Note that iffkgs is also used by ddk and 2nd derivative parts
                 iffkgs=iffkgs+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL2_6)
                 dgxds(1,jj-4,ia,iproj)=dgxds(1,jj-4,ia,iproj)+scalars(1,iffkgs)
                 dgxds(2,jj-4,ia,iproj)=dgxds(2,jj-4,ia,iproj)+scalars(2,iffkgs)
!                !$OMP END CRITICAL (OPERNL2_6)
               end do
             end if

             if (choice==6) then
!              Compute additional tensors related to strain 2nd derivatives
!              and internal strain derivatives
!              ilang6 is number of unique tensor components of rank ilang+3
               ilang6=((ilang+4)*(ilang+5))/2
               jjs=((ilang+3)*(ilang+4)*(ilang+5))/6
!              Compute strain gradient tensor components
               do ii=1,ilang6
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL2_6)
                 d2gxds2(1,jj-20,ia,iproj)=d2gxds2(1,jj-20,ia,iproj)+scalars(1,iffkgs2)
                 d2gxds2(2,jj-20,ia,iproj)=d2gxds2(2,jj-20,ia,iproj)+scalars(2,iffkgs2)
!                !$OMP END CRITICAL (OPERNL2_6)
               end do

!              ilang4 is number of unique tensor components of rank ilang
               ilang4=((ilang+1)*(ilang+2))/2
               jjs=((ilang)*(ilang+1)*(ilang+2))/6
!              Compute internal strain gradient tensor components
               do ii=1,ilang4
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL2_6)
!                Pay attention to the use of reals and imaginary parts here ...
                 dgxdis(1,jj-1,ia,iproj)=dgxdis(1,jj-1,ia,iproj)-two_pi*scalars(2,iffkgs2)
                 dgxdis(2,jj-1,ia,iproj)=dgxdis(2,jj-1,ia,iproj)+two_pi*scalars(1,iffkgs2)
!                !$OMP END CRITICAL (OPERNL2_6)
               end do

!              ilang5 is number of unique tensor components of rank ilang+2
               ilang5=((ilang+3)*(ilang+4))/2
               jjs=((ilang+2)*(ilang+3)*(ilang+4))/6
!              Compute internal strain gradient tensor components
               do ii=1,ilang5
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL2_6)
!                Pay attention to the use of reals and imaginary parts here ...
                 d2gxdis(1,jj-10,ia,iproj)=d2gxdis(1,jj-10,ia,iproj)-two_pi*scalars(2,iffkgs2)
                 d2gxdis(2,jj-10,ia,iproj)=d2gxdis(2,jj-10,ia,iproj)+two_pi*scalars(1,iffkgs2)
!                !$OMP END CRITICAL (OPERNL2_6)
               end do
             end if ! choice==6

             if (choice==5) then
!              Compute additional tensors related to ddk with ffnl(:,2,...)
               ilangx=(ilang*(ilang+1))/2
               jjs=((ilang-1)*ilang*(ilang+1))/6
               do ii=1,ilangx
!                Note that iffkgs is also used by strain part
                 iffkgs=iffkgs+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL2_7)
                 dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)+scalars(1,iffkgs)
                 dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+scalars(2,iffkgs)
!                !$OMP END CRITICAL (OPERNL2_7)
               end do
!              Compute additional tensors related to ddk with ffnl(:,1,...)
               if(ilang>=2)then
                 ilangx=((ilang-1)*ilang)/2
                 jjs=((ilang-2)*(ilang-1)*ilang)/6
                 do ii=1,ilangx
                   iffkgk=iffkgk+1
                   jj=ii+jjs
!                  !$OMP CRITICAL (OPERNL2_8)
                   dgxdt(1,2,jj,ia,iproj)=dgxdt(1,2,jj,ia,iproj)+scalars(1,iffkgk)
                   dgxdt(2,2,jj,ia,iproj)=dgxdt(2,2,jj,ia,iproj)+scalars(2,iffkgk)
!                  !$OMP END CRITICAL (OPERNL2_8)
                 end do
               end if
             end if

!            End projector loop
           end do

!          End condition of non-zero projectors
         end if

!        End angular momentum loop
       end do

!      End loop on atoms
     end do

   else if (sign==-1) then
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
!                see Eq.(53) of PRB55,10337(1997) [[cite:Gonze1997]]
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
!      if(choice==5)then
!      write(std_out,*)' opernl2 : write gxa(:,...) array '
!      do ii=1,10
!      write(std_out,'(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!      end do
!      write(std_out,*)' opernl2 : write dgxdt(:,1,...) array '
!      do ii=1,10
!      write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!      end do
!      write(std_out,*)' opernl2 : write dgxdt(:,2,...) array '
!      do ii=1,4
!      write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!      end do
!      end if

!      do iffkg=1,nffkg
!      write(std_out,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!      end do
!      stop
!      ENDDEBUG

!      ******* Entering the critical part *********************************

       do ipw=1,nincpw
         teffv(1,ipw)=0.0d0
         teffv(2,ipw)=0.0d0
       end do
       do iffkg=1,nffkg
         do ipw=1,nincpw
           teffv(1,ipw)=teffv(1,ipw)+ffkg(ipw,iffkg)*scalars(1,iffkg)
           teffv(2,ipw)=teffv(2,ipw)+ffkg(ipw,iffkg)*scalars(2,iffkg)
         end do
       end do
!      Multiplication by the complex conjugate of the phase
       ig=ipw1
       do ipw=1,nincpw
         vect(1,ig)=vect(1,ig)+teffv(1,ipw)*ph3d(1,ig,iaph3d)+teffv(2,ipw)*ph3d(2,ig,iaph3d)
         vect(2,ig)=vect(2,ig)+teffv(2,ipw)*ph3d(1,ig,iaph3d)-teffv(1,ipw)*ph3d(2,ig,iaph3d)
         ig=ig+1
       end do

!      ******* Leaving the critical part *********************************

!      End loop on atoms
     end do

!    End sign==-1
   else

!    Problem: sign and/or choice do not make sense
     write(message,'(a,2i10,a,a)')&
&     ' Input sign, choice=',sign,choice,ch10,&
&     ' Are not compatible or allowed. '
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


!DEBUG
!if(choice==5)then
!write(std_out,*)'opernl2 : write vect(2*npw)'
!do ipw=1,2
!write(std_out,*)ipw,vect(1:2,ipw)
!end do
!write(std_out,*)'opernl2 : write ph3d'
!do ipw=1,npw
!write(std_out,*)ipw,ph3d(1:2,ipw,1)
!end do
!write(std_out,*)' opernl2 : write gxa array '
!write(std_out,*)' ang mom ,ia '
!do iproj=1,mproj
!do ia=1,1
!do ii=1,3
!do ii=1,10
!write(std_out,'(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!end do
!end do
!end do
!end if
!if(choice==5)then
!write(std_out,*)' opernl2 : write dgxdt(:,1,...) array '
!do ii=1,10
!write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!end do
!write(std_out,*)' opernl2 : write dgxdt(:,2,...) array '
!do ii=1,4
!write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!end do
!stop
!end if
!ENDDEBUG

end subroutine opernl2
!!***

!!****f* ABINIT/opernl3
!! NAME
!! opernl3
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! either from reciprocal space to projected quantities (sign=1),
!! or from projected quantities to reciprocal space (sign=-1)
!!
!! INPUTS
!!  if(sign==-1 .and. (choice==2 .or choice==4 .or. choice==5))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==-1 .and. choice==3)
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the perturbation (needed if choice==2 or 5, and ndgxdt=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  istwf_k=option parameter that describes the storage of wfs
!!  itypat = type of atom, needed for ffnl
!!  jproj(nlang)=number of projectors for each angular momentum
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  matblk=dimension of the array ph3d
!!  mincat= maximum increment of atoms
!!  mlang1 = dimensions for dgxdis1
!!  mlang3 = one of the dimensions of the array gxa
!!  mlang4 = dimension for dgxds
!!  mlang5 = dimensions for dgxdis2
!!  mlang6 = dimension for d2gxds2
!!  mproj=maximum dimension for number of projection operators for each
!!    angular momentum for nonlocal pseudopotential
!!  ndgxdt=second dimension of dgxdt
!!  nffnl=third dimension of ffnl
!!  nincat = number of atoms in the subset here treated
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = number of plane waves in reciprocal space
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!  sign : if  1, go from reciprocal space to projected scalars,
!!         if -1, go from projected scalars to reciprocal space.
!!  if(sign==1),
!!   vect(2*npw)=starting vector in reciprocal space
!!  if(sign==-1)
!!   gxa(2,mlang3,nincat,mproj)=modified projected scalars;
!!   NOTE that metric contractions have already been performed on the
!!   arrays gxa if sign=-1
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!
!! OUTPUT
!!  if(sign==1)
!!   gxa(2,mlang3,mincat,mproj)= projected scalars
!!  if(sign==1 .and. (choice==2 .or choice==4 .or. choice==5 .or. choice==23))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==1 .and. (choice==3 .or. choice==23))
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  if(sign==1 .and. choice==6)
!!   dgxdis((2,mlang1,mincat,mproj) = derivatives of projected scalars
!!    wrt coord. indexed for internal strain
!!   d2gxdis((2,mlang5,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strain and coord
!!   d2gxds2((2,mlang6,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strains
!!  if(sign==-1)
!!   vect(2*npw)=final vector in reciprocal space <G|V_nonlocal|vect_start>.
!!
!! NOTES
!! Operate with the non-local part of the hamiltonian for one type of
!! atom, and within this given type of atom, for a subset
!! of at most nincat atoms.
!!
!! This routine basically replaces getgla (gxa here is the former gla),
!! except for the calculation of <G|dVnl/dk|C> or strain gradients.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      dfpt_mkffkg
!!
!! SOURCE

subroutine opernl3(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&
&  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&
&  ntypat,ph3d,sign,vect)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,ia3,idir,ispinor,istwf_k,itypat,lmnmax,matblk
 integer,intent(in) :: mincat,mlang1,mlang3,mlang4,mlang5,mlang6,mproj,ndgxdt
 integer,intent(in) :: nffnl,nincat,nkpg,nlang,npw,ntypat,sign
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),jproj(nlang),kg_k(3,npw)
 integer,intent(in) :: nloalg(3)
 real(dp),intent(in) :: ffnl(1,npw,nffnl,lmnmax,ntypat),gmet(3,3)
 real(dp),intent(in) :: kpg_k(npw,nkpg),kpt(3),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj)
 real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
 real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj),vect(:,:)
 real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
 real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
 real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)

!Local variables-------------------------------
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
!to handle strain 2DER, need up to rank 7
!scalars
 integer :: ia,iaph3d,iffkg,iffkgk,iffkgs,iffkgs2,ig,ii,ilang,ilang2,ilang3
 integer :: ilang4,ilang5,ilang6,ilangx,iproj,ipw,ipw1,ipw2,jffkg,jj,jjs,mblkpw
 integer :: mmproj,mu,nffkg,nffkgd,nffkge,nffkgk,nffkgs,nffkgs2,nincpw,nproj
 integer :: ntens
 real(dp) :: ai,ar
 real(dp),parameter :: two_pi2=two_pi**2
 character(len=500) :: message
!arrays
 integer,allocatable :: parity(:)
! real(dp) :: tsec(2)
 real(dp),allocatable :: ffkg(:,:),kpgx(:,:),scalars(:,:),teffv(:,:)

! *************************************************************************

!call wrtout(std_out,"in opernl3","COLL")

!mblkpw sets the size of blocks of planewaves
 mblkpw=NLO_MBLKPW

!Get the actual maximum number of projectors
 mmproj=maxval(indlmn(3,:,itypat))

!Initialisation before blocking on the plane waves

 if (sign==1) then
!  Put projected scalars to zero
   gxa(:,:,:,1:mmproj)=0.0d0
   if (choice==2 .or. choice==4 .or. choice==5 .or. choice==23) dgxdt(:,:,:,:,1:mmproj)=0.0d0
   if (choice==3 .or. choice==6 .or. choice==23) dgxds(:,:,:,1:mmproj)=0.0d0
   if (choice==6) then
     dgxdis(:,:,:,1:mmproj)=0.0d0
     d2gxdis(:,:,:,1:mmproj)=0.0d0
     d2gxds2(:,:,:,1:mmproj)=0.0d0
   end if
 end if

!Set up dimension of kpgx and allocate
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
 ntens=1
 if(nlang>=2 .or. choice==2 .or. choice==4 .or. choice==5 .or. choice==23) ntens=4
 if(nlang>=3 .or. (choice==3.or.choice==23))ntens=10
 if(nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) )ntens=20
 if(((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6)ntens=35
 if(((choice==3.or.choice==23) .and. nlang==4) .or. (choice==6 .and. nlang>=2))ntens=56
 if(choice==6 .and. nlang>=3)ntens=84
 if(choice==6 .and. nlang==4)ntens=120

!Set up second dimension of ffkg array, and allocate
 nffkg=0 ; nffkge=0 ; nffkgd=0 ; nffkgk=0 ; nffkgs=0 ; nffkgs2=0
 do ilang=1,nlang
!  Get the number of projectors for that angular momentum
   nproj=jproj(ilang)
!  If there is a non-local part, accumulate the number of vectors needed
!  The variables ilang below are the number of independent tensors of
!  various ranks, the variable names being more historical than logical.
!  ilang2=number of rank ilang-1
!  ilang3=number of rank ilang+1
!  ilang4=number of rank ilang
!  ilang5=number of rank ilang+2
!  ilang6=number of rank ilang+3
   if(nproj>0)then
     ilang2=(ilang*(ilang+1))/2
     nffkge=nffkge+nproj*ilang2
     if(choice==5)nffkgk=nffkgk+nproj*(2*ilang2-ilang)
     if(choice==2 .or. choice==4 .or. choice==23)nffkgd=nffkgd+ndgxdt*nproj*ilang2
     if(choice==3 .or. choice==6 .or. choice==23)then
       ilang3=((ilang+2)*(ilang+3))/2
       nffkgs=nffkgs+nproj*ilang3
     end if
     if(choice==6)then
       ilang4=((ilang+1)*(ilang+2))/2
       ilang5=((ilang+3)*(ilang+4))/2
       ilang6=((ilang+4)*(ilang+5))/2
       nffkgs2=nffkgs2+nproj*(ilang4+ilang5+ilang6)
     end if
   end if
 end do
 nffkg=nffkge+nffkgd+nffkgs+nffkgs2+nffkgk

!Disabled by MG on Dec  6 2011, omp sections have to be tested, this coding causes a
!sigfault with nthreads==1
!
!Loop on subsets of plane waves (blocking)
!!$OMP PARALLEL DEFAULT(PRIVATE) &
!!$OMP SHARED(choice,dgxds,dgxdt,ffnl,gmet,gxa,ia3,idir,indlmn,ispinor) &
!!$OMP SHARED(istwf_k,itypat,jproj,kg_k,kpg_k,kpt,lmnmax,mblkpw,mproj) &
!!$OMP SHARED(ndgxdt,nffkg,nffkgd,nffkge,nffkgs,nincat,nkpg,nlang) &
!!$OMP SHARED(nloalg,ph3d,npw,ntens,ntypat,sign,vect)

 ABI_ALLOCATE(ffkg,(nffkg,mblkpw))
 ABI_ALLOCATE(parity,(nffkg))
 ABI_ALLOCATE(kpgx,(mblkpw,ntens))
 ABI_ALLOCATE(scalars,(2,nffkg))
 ABI_ALLOCATE(teffv,(2,mblkpw))
!!$OMP DO
 do ipw1=1,npw,mblkpw

   ipw2=min(npw,ipw1+mblkpw-1)
   nincpw=ipw2-ipw1+1

!  Initialize kpgx array related to tensors defined below
   call dfpt_mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&
&   kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&   npw,ntens,ntypat,parity)

!  call timab(74+choice,2,tsec)

!  Now treat the different signs
   if (sign==1) then

     do ia=1,nincat

!      Compute the shift eventually needed to get the phases in ph3d
       iaph3d=ia
       if(nloalg(2)>0)iaph3d=ia+ia3-1

       do iffkg=1,nffkg
         scalars(1,iffkg)=0.0d0
         scalars(2,iffkg)=0.0d0
       end do

!      ******* Entering the first time-consuming part of the routine *******

!      Note : first multiply by the phase factor
!      This allows to be left with only real operations afterwards.
       ig=ipw1

       do ipw=1,nincpw
         ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
         ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
         do iffkg=1,nffkg
           scalars(1,iffkg)=scalars(1,iffkg)+ar*ffkg(iffkg,ipw)
           scalars(2,iffkg)=scalars(2,iffkg)+ai*ffkg(iffkg,ipw)
         end do
         ig=ig+1
       end do

!      ******* Leaving the critical part *********************************

!      DEBUG
!      write(std_out,*)' opernl2 : scalars'
!      do iffkg=1,10
!      write(std_out,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!      end do
!      stop
!      ENDDEBUG

       if(istwf_k>=2)then
!        Impose parity of resulting scalar (this operation could be
!        replaced by direct saving of CPU time in the preceeding section)
         do iffkg=1,nffkg
           scalars(parity(iffkg),iffkg)=0.0d0
         end do
       end if

       iffkg=0
       iffkgs=nffkge+nffkgd
       iffkgs2=nffkge+nffkgs
       iffkgk=nffkge*2
       do ilang=1,nlang
         nproj=jproj(ilang)
         if(nproj>0)then
!          ilang2 is the number of independent tensor components
!          for symmetric tensor of rank ilang-1
           ilang2=(ilang*(ilang+1))/2

!          Loop over projectors
           do iproj=1,nproj
!            Multiply by the k+G factors (tensors of various rank)
             do ii=1,ilang2
!              Get the starting address for the relevant tensor
               jj=ii+((ilang-1)*ilang*(ilang+1))/6
               iffkg=iffkg+1
!              !$OMP CRITICAL (OPERNL3_1)
               gxa(1,jj,ia,iproj)=gxa(1,jj,ia,iproj)+scalars(1,iffkg)
               gxa(2,jj,ia,iproj)=gxa(2,jj,ia,iproj)+scalars(2,iffkg)
!              !$OMP END CRITICAL (OPERNL3_1)
!              Now, compute gradients, if needed.
               if ((choice==2.or.choice==23) .and. ndgxdt==3) then
                 do mu=1,3
                   jffkg=nffkge+(iffkg-1)*3+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  !$OMP CRITICAL (OPERNL3_2)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                  !$OMP END CRITICAL (OPERNL3_2)
                 end do
               end if
               if (choice==2 .and. ndgxdt==1) then
                 jffkg=nffkge+iffkg
!                Pay attention to the use of reals and imaginary parts here ...
!                !$OMP CRITICAL (OPERNL3_3)
                 dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                 dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                !$OMP END CRITICAL (OPERNL3_3)
               end if
               if (choice==4) then
                 do mu=1,3
                   jffkg=nffkge+(iffkg-1)*9+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  !$OMP CRITICAL (OPERNL3_4)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                  !$OMP END CRITICAL (OPERNL3_4)
                 end do
                 do mu=4,9
                   jffkg=nffkge+(iffkg-1)*9+mu
!                  Pay attention to the use of reals and imaginary parts here ...
!                  Also, note the multiplication by (2 pi)**2
!                  !$OMP CRITICAL (OPERNL3_5)
                   dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi2*scalars(1,jffkg)
                   dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)-two_pi2*scalars(2,jffkg)
!                  !$OMP END CRITICAL (OPERNL3_5)
                 end do
               end if
!              End loop on ii=1,ilang2
             end do

             if (choice==3 .or. choice==6 .or. choice==23) then
!              Compute additional tensors related to strain gradients
!              ilang3 is number of unique tensor components of rank ilang+1
               ilang3=((ilang+2)*(ilang+3))/2
               jjs=((ilang+1)*(ilang+2)*(ilang+3))/6
!              Compute strain gradient tensor components
               do ii=1,ilang3
!                Note that iffkgs is also used by ddk and 2nd derivative parts
                 iffkgs=iffkgs+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL3_6)
                 dgxds(1,jj-4,ia,iproj)=dgxds(1,jj-4,ia,iproj)+scalars(1,iffkgs)
                 dgxds(2,jj-4,ia,iproj)=dgxds(2,jj-4,ia,iproj)+scalars(2,iffkgs)
!                !$OMP END CRITICAL (OPERNL3_6)
               end do
             end if

             if (choice==6) then
!              Compute additional tensors related to strain 2nd derivatives
!              and internal strain derivatives
!              ilang6 is number of unique tensor components of rank ilang+3
               ilang6=((ilang+4)*(ilang+5))/2
               jjs=((ilang+3)*(ilang+4)*(ilang+5))/6
!              Compute strain gradient tensor components
               do ii=1,ilang6
!                Note that iffkgs is also used by ddk part
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL3_6)
                 d2gxds2(1,jj-20,ia,iproj)=d2gxds2(1,jj-20,ia,iproj)+scalars(1,iffkgs2)
                 d2gxds2(2,jj-20,ia,iproj)=d2gxds2(2,jj-20,ia,iproj)+scalars(2,iffkgs2)
!                !$OMP END CRITICAL (OPERNL3_6)
               end do

!              ilang4 is number of unique tensor components of rank ilang
               ilang4=((ilang+1)*(ilang+2))/2
               jjs=((ilang)*(ilang+1)*(ilang+2))/6
!              Compute internal strain gradient tensor components
               do ii=1,ilang4
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL3_6)
!                Pay attention to the use of reals and imaginary parts here ...
                 dgxdis(1,jj-1,ia,iproj)=dgxdis(1,jj-1,ia,iproj)-two_pi*scalars(2,iffkgs2)
                 dgxdis(2,jj-1,ia,iproj)=dgxdis(2,jj-1,ia,iproj)+two_pi*scalars(1,iffkgs2)
!                !$OMP END CRITICAL (OPERNL3_6)
               end do

!              ilang5 is number of unique tensor components of rank ilang+2
               ilang5=((ilang+3)*(ilang+4))/2
               jjs=((ilang+2)*(ilang+3)*(ilang+4))/6
!              Compute internal strain gradient tensor components
               do ii=1,ilang5
                 iffkgs2=iffkgs2+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL3_6)
!                Pay attention to the use of reals and imaginary parts here ...
                 d2gxdis(1,jj-10,ia,iproj)=d2gxdis(1,jj-10,ia,iproj)-two_pi*scalars(2,iffkgs2)
                 d2gxdis(2,jj-10,ia,iproj)=d2gxdis(2,jj-10,ia,iproj)+two_pi*scalars(1,iffkgs2)
!                !$OMP END CRITICAL (OPERNL3_6)
               end do
             end if ! choice==6

             if (choice==5) then
!              Compute additional tensors related to ddk with ffnl(:,2,...)
               ilangx=(ilang*(ilang+1))/2
               jjs=((ilang-1)*ilang*(ilang+1))/6
               do ii=1,ilangx
!                Note that iffkgs is also used by strain part
                 iffkgs=iffkgs+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL3_7)
                 dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)+scalars(1,iffkgs)
                 dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+scalars(2,iffkgs)
!                !$OMP END CRITICAL (OPERNL3_7)
               end do
!              Compute additional tensors related to ddk with ffnl(:,1,...)
               if(ilang>=2)then
                 ilangx=((ilang-1)*ilang)/2
                 jjs=((ilang-2)*(ilang-1)*ilang)/6
                 do ii=1,ilangx
                   iffkgk=iffkgk+1
                   jj=ii+jjs
!                  !$OMP CRITICAL (OPERNL3_8)
                   dgxdt(1,2,jj,ia,iproj)=dgxdt(1,2,jj,ia,iproj)+scalars(1,iffkgk)
                   dgxdt(2,2,jj,ia,iproj)=dgxdt(2,2,jj,ia,iproj)+scalars(2,iffkgk)
!                  !$OMP END CRITICAL (OPERNL3_8)
                 end do
               end if
             end if

!            End projector loop
           end do

!          End condition of non-zero projectors
         end if

!        End angular momentum loop
       end do

!      End loop on atoms
     end do

   else if (sign==-1) then
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
!                see Eq.(53) of PRB55,10337(1997) [[cite:Gonze1997]]
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
!      if(choice==5)then
!      write(std_out,*)' opernl2 : write gxa(:,...) array '
!      do ii=1,10
!      write(std_out,'(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!      end do
!      write(std_out,*)' opernl2 : write dgxdt(:,1,...) array '
!      do ii=1,10
!      write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!      end do
!      write(std_out,*)' opernl2 : write dgxdt(:,2,...) array '
!      do ii=1,4
!      write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!      end do
!      end if

!      do iffkg=1,nffkg
!      write(std_out,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!      end do
!      stop
!      ENDDEBUG

       ig=ipw1

!      ******* Entering the critical part *********************************

       do ipw=1,nincpw
         ar=0.0d0
         ai=0.0d0
         do iffkg=1,nffkg
           ar=ar+ffkg(iffkg,ipw)*scalars(1,iffkg)
           ai=ai+ffkg(iffkg,ipw)*scalars(2,iffkg)
         end do
         vect(1,ig)=vect(1,ig)+ar*ph3d(1,ig,iaph3d)+ai*ph3d(2,ig,iaph3d)
         vect(2,ig)=vect(2,ig)+ai*ph3d(1,ig,iaph3d)-ar*ph3d(2,ig,iaph3d)
         ig=ig+1
       end do

!      ******* Leaving the critical part *********************************

!      End loop on atoms
     end do

!    End sign==-1
   else

!    Problem: sign and/or choice do not make sense
     write(message, '(a,2i10,a,a)' )&
&     '  Input sign, choice=',sign,choice,ch10,&
&     '  Are not compatible or allowed. '
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


!DEBUG
!if(choice==5)then
!write(std_out,*)'opernl2 : write vect(2*npw)'
!do ipw=1,2
!write(std_out,*)ipw,vect(1:2,ipw)
!end do
!write(std_out,*)'opernl2 : write ph3d'
!do ipw=1,npw
!write(std_out,*)ipw,ph3d(1:2,ipw,1)
!end do
!write(std_out,*)' opernl2 : write gxa array '
!write(std_out,*)' ang mom ,ia '
!do iproj=1,mproj
!do ia=1,1
!do ii=1,3
!do ii=1,10
!write(std_out,'(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!end do
!end do
!end do
!end if
!if(choice==5)then
!write(std_out,*)' opernl2 : write dgxdt(:,1,...) array '
!do ii=1,10
!write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!end do
!write(std_out,*)' opernl2 : write dgxdt(:,2,...) array '
!do ii=1,4
!write(std_out,'(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!end do
!stop
!end if
!ENDDEBUG

end subroutine opernl3
!!***


!!****f* ABINIT/opernl4a
!! NAME
!! opernl4a
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! from reciprocal space to projected quantities
!! (oprnl4b is from projected quantities to reciprocal space)
!!
!! INPUTS
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the perturbation (needed if choice==2 or 5, and ndgxdt=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  istwf_k=option parameter that describes the storage of wfs
!!  itypat = type of atom, needed for ffnl
!!  jproj(nlang)=number of projectors for each angular momentum
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  matblk=dimension of the array ph3d
!!  mincat= maximum increment of atoms
!!  mlang1 = dimensions for dgxdis1
!!  mlang3 = one of the dimensions of the array gxa
!!  mlang4 = dimension for dgxds
!!  mlang5 = dimensions for dgxdis2
!!  mlang6 = dimension for d2gxds2
!!  mproj=maximum dimension for number of projection operators for each
!!    angular momentum for nonlocal pseudopotential
!!  ndgxdt=second dimension of dgxdt
!!  nffnl=third dimension of ffnl
!!  nincat = number of atoms in the subset here treated
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = number of plane waves in reciprocal space
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!  vect(2*npw)=starting vector in reciprocal space
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!
!! OUTPUT
!!  gxa(2,mlang3,mincat,mproj)= projected scalars
!!  if(choice==2 .or choice==4 .or. choice==5 .or. choice==23)
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= gradients of projected scalars wrt coords
!!    or with respect to ddk
!!  if(choice==3 .or. choice==23)
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  if(choice==6)
!!   dgxdis((2,mlang1,mincat,mproj) = derivatives of projected scalars
!!    wrt coord. indexed for internal strain
!!   d2gxdis((2,mlang5,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strain and coord
!!   d2gxds2((2,mlang6,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strains
!!
!! NOTES
!! Operate with the non-local part of the hamiltonian for one type of
!! atom, and within this given type of atom, for a subset
!! of at most nincat atoms.
!! This routine basically replaces getgla (gxa here is the former gla),
!! except for the calculation of <G|dVnl/dk|C> or strain gradients.
!!
!! Present version decomposed according to iffkg
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      dfpt_mkffkg
!!
!! SOURCE

subroutine opernl4a(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&
&  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&
&  ntypat,ph3d,vect)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,ia3,idir,ispinor,istwf_k,itypat,lmnmax,matblk
 integer,intent(in) :: mincat,mlang1,mlang3,mlang4,mlang5,mlang6,mproj,ndgxdt
 integer,intent(in) :: nffnl,nincat,nkpg,nlang,npw,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),jproj(nlang),kg_k(3,npw)
 integer,intent(in) :: nloalg(3)
 real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat),gmet(3,3),kpg_k(npw,nkpg)
 real(dp),intent(in) :: kpt(3),ph3d(2,npw,matblk),vect(:,:)
 real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
 real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
 real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)
 real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj) !vz_i
 real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj) !vz_i
 real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj) !vz_i

!Local variables-------------------------------
!scalars
 integer :: chunk,ia,iaph3d,iffkg,iffkgk,iffkgs,iffkgs2,ig,ii,ilang,ilang2
 integer :: ilang3,ilang4,ilang5,ilang6,ilangx,iproj,ipw,ipw1,ipw2,jffkg,jj,jjs
 integer :: jump,mblkpw,mmproj,mu,nffkg,nffkgd,nffkge,nffkgk,nffkgs,nffkgs2
 integer :: nincpw,nproj,ntens,start
 real(dp) :: ai,ar,sci1,sci2,sci3,sci4,sci5,sci6,sci7,sci8
 real(dp) :: scr1,scr2,scr3,scr4,scr5,scr6,scr7,scr8
 real(dp),parameter :: two_pi2=two_pi*two_pi
!arrays
 integer,allocatable :: parity(:)
 real(dp),allocatable :: ffkg(:,:),kpgx(:,:),scalars(:,:),teffv(:,:)

! *************************************************************************

!call wrtout(std_out,"in opernl4a","COLL")

!mblkpw sets the size of blocks of planewaves
 mblkpw=NLO_MBLKPW

!jump governs, in fine, the use of registers in the most cpu
!time consuming part of the routine. Until now, jump=8 is the maximal value.
!The optimal value will be machine-dependent !
 jump=4

!Get the actual maximum number of projectors
 mmproj=maxval(indlmn(3,:,itypat))

!Initialisation before blocking on the plane waves

!Put projected scalars to zero
 gxa(:,:,:,1:mmproj)=0.0d0
 if (choice==2 .or. choice==4 .or. choice==5 .or. choice==23) dgxdt(:,:,:,:,1:mmproj)=0.0d0
 if (choice==3 .or. choice==6 .or. choice==23) dgxds(:,:,:,1:mmproj)=0.0d0
 if (choice==6) then
   dgxdis(:,:,:,1:mmproj)=0.0d0
   d2gxdis(:,:,:,1:mmproj)=0.0d0
   d2gxds2(:,:,:,1:mmproj)=0.0d0
 end if

!Set up dimension of kpgx and allocate
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
 ntens=1
 if(nlang>=2 .or. choice==2 .or. choice==4 .or. choice==5 .or. choice==23)ntens=4
 if(nlang>=3 .or. (choice==3.or.choice==23))ntens=10
 if(nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) )ntens=20
 if(((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6)ntens=35
 if(((choice==3.or.choice==23) .and. nlang==4) .or. (choice==6 .and. nlang>=2))ntens=56
 if(choice==6 .and. nlang>=3)ntens=84
 if(choice==6 .and. nlang==4)ntens=120

!Set up second dimension of ffkg array, and allocate
 nffkg=0 ; nffkge=0 ; nffkgd=0 ; nffkgk=0 ; nffkgs=0 ; nffkgs2=0
 do ilang=1,nlang
!  Get the number of projectors for that angular momentum
   nproj=jproj(ilang)
!  If there is a non-local part, accumulate the number of vectors needed
!  The variables ilang below are the number of independent tensors of
!  various ranks, the variable names being more historical than logical.
!  ilang2=number of rank ilang-1
!  ilang3=number of rank ilang+1
!  ilang4=number of rank ilang
!  ilang5=number of rank ilang+2
!  ilang6=number of rank ilang+3
   if(nproj>0)then
     ilang2=(ilang*(ilang+1))/2
     nffkge=nffkge+nproj*ilang2
     if(choice==5)nffkgk=nffkgk+nproj*(2*ilang2-ilang)
     if(choice==2 .or. choice==4 .or. choice==23)nffkgd=nffkgd+ndgxdt*nproj*ilang2
     if(choice==3 .or. choice==6 .or. choice==23)then
       ilang3=((ilang+2)*(ilang+3))/2
       nffkgs=nffkgs+nproj*ilang3
     end if
     if(choice==6)then
       ilang4=((ilang+1)*(ilang+2))/2
       ilang5=((ilang+3)*(ilang+4))/2
       ilang6=((ilang+4)*(ilang+5))/2
       nffkgs2=nffkgs2+nproj*(ilang4+ilang5+ilang6)
     end if
   end if
 end do
 nffkg=nffkge+nffkgd+nffkgs+nffkgs2+nffkgk

!DEBUG
!write(std_out,*)' jproj(1:nlang)',jproj(1:nlang)
!write(std_out,*)' nffkg,nffkge,nffkgd,nffkgs,nffkgk',nffkg,nffkge,nffkgd,nffkgs,nffkgk
!ENDDEBUG

!Loop on subsets of plane waves (blocking)

!Disabled by MG on Dec  6 2011, omp sections have to be tested, this coding causes a
!sigfault with nthreads==1
!Feb 16 2012: The code does not crash anymore but it's not efficient.
!
!!$OMP PARALLEL DEFAULT(PRIVATE) &
!!$OMP SHARED (choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt) &
!!$OMP SHARED (ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat) &
!!$OMP SHARED (jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4) &
!!$OMP SHARED (mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw) &
!!$OMP SHARED (ntypat,ph3d,vect) &
!!$OMP SHARED (mblkpw,jump,nffkgd,nffkg,nffkge,nffkgs,ntens)

 ABI_ALLOCATE(ffkg,(nffkg,mblkpw))
 ABI_ALLOCATE(parity,(nffkg))
 ABI_ALLOCATE(kpgx,(mblkpw,ntens))
 ABI_ALLOCATE(scalars,(2,nffkg))
 ABI_ALLOCATE(teffv,(2,mblkpw))

!!$OMP DO
 do ipw1=1,npw,mblkpw

   ipw2=min(npw,ipw1+mblkpw-1)
   nincpw=ipw2-ipw1+1

!  Initialize kpgx array related to tensors defined below
   call dfpt_mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&
&   kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&   npw,ntens,ntypat,parity)

   do ia=1,nincat

!    Compute the shift eventually needed to get the phases in ph3d
     iaph3d=ia
     if(nloalg(2)>0)iaph3d=ia+ia3-1

     do iffkg=1,nffkg
       scalars(1,iffkg)=0.0d0 ; scalars(2,iffkg)=0.0d0
     end do

!    DEBUG
!    write(std_out,*)'opernl4, before first time-consuming'
!    write(std_out,*)'opernl4 : nffkg,nincpw=',nffkg,nincpw
!    write(std_out,*)'ig,ipw,ffkg(1:4),vect(1:2)'
!    ig=ipw1
!    do ipw=1,nincpw
!    write(std_out,'(2i4,13es11.3)' )ig,ipw,ffkg(1:min(9,nffkg),ipw),vect(1:2,ipw),ph3d(1:2,ipw,iaph3d)
!    ig=ig+1
!    end do
!    stop
!    ENDDEBUG

!    ******* Entering the first time-consuming part of the routine *******


!    First, treat small nffkg; send treat the initial phase of big
!    nffkg; finally treat the loop needed for big nffkg

!    In the loops, first multiply by the phase factor.
!    This allows to be left with only real operations afterwards.

!    For the time being, the maximal jump allowed is 8.

!    1) Here, treat small nffkg
     if(nffkg<=jump)then

       select case(nffkg)

       case(1)

         scr1=0.0d0 ; sci1=0.0d0
         ig=ipw1
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1

       case(2)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2

       case(3)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3

       case(4)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4

       case(5)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5

       case(6)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6

       case(7)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         scr7=0.0d0 ; sci7=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           scr7=scr7+ar*ffkg(7,ipw) ; sci7=sci7+ai*ffkg(7,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6
         scalars(1,7)=scr7 ; scalars(2,7)=sci7

       case(8)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         scr7=0.0d0 ; sci7=0.0d0
         scr8=0.0d0 ; sci8=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           scr7=scr7+ar*ffkg(7,ipw) ; sci7=sci7+ai*ffkg(7,ipw)
           scr8=scr8+ar*ffkg(8,ipw) ; sci8=sci8+ai*ffkg(8,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6
         scalars(1,7)=scr7 ; scalars(2,7)=sci7
         scalars(1,8)=scr8 ; scalars(2,8)=sci8

       end select

     else
!      Now treat big nffkg

!      2) Here, initialize big nffkg. The only difference with the
!      preceeding case is that the intermediate results are stored.

       select case(jump)

       case(1)

         scr1=0.0d0 ; sci1=0.0d0
         ig=ipw1
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1

       case(2)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2

       case(3)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3

       case(4)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4

       case(5)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5

       case(6)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6

       case(7)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         scr7=0.0d0 ; sci7=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           scr7=scr7+ar*ffkg(7,ipw) ; sci7=sci7+ai*ffkg(7,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6
         scalars(1,7)=scr7 ; scalars(2,7)=sci7

       case(8)

         ig=ipw1
         scr1=0.0d0 ; sci1=0.0d0
         scr2=0.0d0 ; sci2=0.0d0
         scr3=0.0d0 ; sci3=0.0d0
         scr4=0.0d0 ; sci4=0.0d0
         scr5=0.0d0 ; sci5=0.0d0
         scr6=0.0d0 ; sci6=0.0d0
         scr7=0.0d0 ; sci7=0.0d0
         scr8=0.0d0 ; sci8=0.0d0
         do ipw=1,nincpw
           ar=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
           ai=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
           teffv(1,ipw)=ar          ; teffv(2,ipw)=ai
           scr1=scr1+ar*ffkg(1,ipw) ; sci1=sci1+ai*ffkg(1,ipw)
           scr2=scr2+ar*ffkg(2,ipw) ; sci2=sci2+ai*ffkg(2,ipw)
           scr3=scr3+ar*ffkg(3,ipw) ; sci3=sci3+ai*ffkg(3,ipw)
           scr4=scr4+ar*ffkg(4,ipw) ; sci4=sci4+ai*ffkg(4,ipw)
           scr5=scr5+ar*ffkg(5,ipw) ; sci5=sci5+ai*ffkg(5,ipw)
           scr6=scr6+ar*ffkg(6,ipw) ; sci6=sci6+ai*ffkg(6,ipw)
           scr7=scr7+ar*ffkg(7,ipw) ; sci7=sci7+ai*ffkg(7,ipw)
           scr8=scr8+ar*ffkg(8,ipw) ; sci8=sci8+ai*ffkg(8,ipw)
           ig=ig+1
         end do
         scalars(1,1)=scr1 ; scalars(2,1)=sci1
         scalars(1,2)=scr2 ; scalars(2,2)=sci2
         scalars(1,3)=scr3 ; scalars(2,3)=sci3
         scalars(1,4)=scr4 ; scalars(2,4)=sci4
         scalars(1,5)=scr5 ; scalars(2,5)=sci5
         scalars(1,6)=scr6 ; scalars(2,6)=sci6
         scalars(1,7)=scr7 ; scalars(2,7)=sci7
         scalars(1,8)=scr8 ; scalars(2,8)=sci8

       end select

!      3) Here, do-loop for big nffkg.

       do start=1+jump,nffkg,jump
         chunk=min(jump,nffkg-start+1)

         select case(chunk)

         case(1)

           scr1=0.0d0 ; sci1=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)              ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start,ipw) ; sci1=sci1+ai*ffkg(start,ipw)
           end do
           scalars(1,start)=scr1 ; scalars(2,start)=sci1

         case(2)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2

         case(3)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3

         case(4)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           scr4=0.0d0 ; sci4=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
             scr4=scr4+ar*ffkg(start+3,ipw) ; sci4=sci4+ai*ffkg(start+3,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3
           scalars(1,start+3)=scr4 ; scalars(2,start+3)=sci4

         case(5)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           scr4=0.0d0 ; sci4=0.0d0
           scr5=0.0d0 ; sci5=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
             scr4=scr4+ar*ffkg(start+3,ipw) ; sci4=sci4+ai*ffkg(start+3,ipw)
             scr5=scr5+ar*ffkg(start+4,ipw) ; sci5=sci5+ai*ffkg(start+4,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3
           scalars(1,start+3)=scr4 ; scalars(2,start+3)=sci4
           scalars(1,start+4)=scr5 ; scalars(2,start+4)=sci5

         case(6)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           scr4=0.0d0 ; sci4=0.0d0
           scr5=0.0d0 ; sci5=0.0d0
           scr6=0.0d0 ; sci6=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
             scr4=scr4+ar*ffkg(start+3,ipw) ; sci4=sci4+ai*ffkg(start+3,ipw)
             scr5=scr5+ar*ffkg(start+4,ipw) ; sci5=sci5+ai*ffkg(start+4,ipw)
             scr6=scr6+ar*ffkg(start+5,ipw) ; sci6=sci6+ai*ffkg(start+5,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3
           scalars(1,start+3)=scr4 ; scalars(2,start+3)=sci4
           scalars(1,start+4)=scr5 ; scalars(2,start+4)=sci5
           scalars(1,start+5)=scr6 ; scalars(2,start+5)=sci6

         case(7)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           scr4=0.0d0 ; sci4=0.0d0
           scr5=0.0d0 ; sci5=0.0d0
           scr6=0.0d0 ; sci6=0.0d0
           scr7=0.0d0 ; sci7=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
             scr4=scr4+ar*ffkg(start+3,ipw) ; sci4=sci4+ai*ffkg(start+3,ipw)
             scr5=scr5+ar*ffkg(start+4,ipw) ; sci5=sci5+ai*ffkg(start+4,ipw)
             scr6=scr6+ar*ffkg(start+5,ipw) ; sci6=sci6+ai*ffkg(start+5,ipw)
             scr7=scr7+ar*ffkg(start+6,ipw) ; sci7=sci7+ai*ffkg(start+6,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3
           scalars(1,start+3)=scr4 ; scalars(2,start+3)=sci4
           scalars(1,start+4)=scr5 ; scalars(2,start+4)=sci5
           scalars(1,start+5)=scr6 ; scalars(2,start+5)=sci6
           scalars(1,start+6)=scr7 ; scalars(2,start+6)=sci7

         case(8)

           scr1=0.0d0 ; sci1=0.0d0
           scr2=0.0d0 ; sci2=0.0d0
           scr3=0.0d0 ; sci3=0.0d0
           scr4=0.0d0 ; sci4=0.0d0
           scr5=0.0d0 ; sci5=0.0d0
           scr6=0.0d0 ; sci6=0.0d0
           scr7=0.0d0 ; sci7=0.0d0
           scr8=0.0d0 ; sci8=0.0d0
           do ipw=1,nincpw
             ar=teffv(1,ipw)                ; ai=teffv(2,ipw)
             scr1=scr1+ar*ffkg(start  ,ipw) ; sci1=sci1+ai*ffkg(start  ,ipw)
             scr2=scr2+ar*ffkg(start+1,ipw) ; sci2=sci2+ai*ffkg(start+1,ipw)
             scr3=scr3+ar*ffkg(start+2,ipw) ; sci3=sci3+ai*ffkg(start+2,ipw)
             scr4=scr4+ar*ffkg(start+3,ipw) ; sci4=sci4+ai*ffkg(start+3,ipw)
             scr5=scr5+ar*ffkg(start+4,ipw) ; sci5=sci5+ai*ffkg(start+4,ipw)
             scr6=scr6+ar*ffkg(start+5,ipw) ; sci6=sci6+ai*ffkg(start+5,ipw)
             scr7=scr7+ar*ffkg(start+6,ipw) ; sci7=sci7+ai*ffkg(start+6,ipw)
             scr8=scr8+ar*ffkg(start+7,ipw) ; sci8=sci8+ai*ffkg(start+7,ipw)
           end do
           scalars(1,start  )=scr1 ; scalars(2,start  )=sci1
           scalars(1,start+1)=scr2 ; scalars(2,start+1)=sci2
           scalars(1,start+2)=scr3 ; scalars(2,start+2)=sci3
           scalars(1,start+3)=scr4 ; scalars(2,start+3)=sci4
           scalars(1,start+4)=scr5 ; scalars(2,start+4)=sci5
           scalars(1,start+5)=scr6 ; scalars(2,start+5)=sci6
           scalars(1,start+6)=scr7 ; scalars(2,start+6)=sci7
           scalars(1,start+7)=scr8 ; scalars(2,start+7)=sci8

         end select

!        End loop on start
       end do

!      End if statement for small or big nffkg
     end if

!    ******* Leaving the critical part *********************************

!    DEBUG
!    write(std_out,*)' opernl4a, write scalars '
!    do iffkg=1,nffkg
!    write(std_out,*)iffkg,scalars(1:2,iffkg)
!    end do
!    ENDDEBUG

     if(istwf_k>=2)then
!      Impose parity of resulting scalar (this operation could be
!      replaced by direct saving of CPU time in the preceeding section)
       do iffkg=1,nffkg
         scalars(parity(iffkg),iffkg)=0.0d0
       end do
     end if

     iffkg=0 ; iffkgs=nffkge+nffkgd ; iffkgk=nffkge*2
     iffkgs2=nffkge+nffkgs
     do ilang=1,nlang
       nproj=jproj(ilang)
       if(nproj>0)then
!        ilang2 is the number of independent tensor components
!        for symmetric tensor of rank ilang-1
         ilang2=(ilang*(ilang+1))/2

!        Loop over projectors
         do iproj=1,nproj
!          Multiply by the k+G factors (tensors of various rank)
           do ii=1,ilang2
!            Get the starting address for the relevant tensor
             jj=ii+((ilang-1)*ilang*(ilang+1))/6
             iffkg=iffkg+1
!            !$OMP CRITICAL (OPERNL4a_1)
             gxa(1,jj,ia,iproj)=gxa(1,jj,ia,iproj)+scalars(1,iffkg)
             gxa(2,jj,ia,iproj)=gxa(2,jj,ia,iproj)+scalars(2,iffkg)
!            !$OMP END CRITICAL (OPERNL4a_1)
!            Now, compute gradients, if needed.
             if ((choice==2.or.choice==23) .and. ndgxdt==3) then
               do mu=1,3
                 jffkg=nffkge+(iffkg-1)*3+mu
!                Pay attention to the use of reals and imaginary parts here ...
!                !$OMP CRITICAL (OPERNL4a_2)
                 dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                 dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                !$OMP END CRITICAL (OPERNL4a_2)
               end do
             end if
             if (choice==2 .and. ndgxdt==1) then
               jffkg=nffkge+iffkg
!              Pay attention to the use of reals and imaginary parts here ...
!              !$OMP CRITICAL (OPERNL4a_3)
               dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)-two_pi*scalars(2,jffkg)
               dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!              !$OMP END CRITICAL (OPERNL4a_3)
             end if
             if (choice==4) then
               do mu=1,3
                 jffkg=nffkge+(iffkg-1)*9+mu
!                Pay attention to the use of reals and imaginary parts here ...
!                !$OMP CRITICAL (OPERNL4a_4)
                 dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
                 dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!                !$OMP END CRITICAL (OPERNL4a_4)
               end do
               do mu=4,9
                 jffkg=nffkge+(iffkg-1)*9+mu
!                Pay attention to the use of reals and imaginary parts here ...
!                Also, note the multiplication by (2 pi)**2
!                !$OMP CRITICAL (OPERNL4a_5)
                 dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi2*scalars(1,jffkg)
                 dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)-two_pi2*scalars(2,jffkg)
!                !$OMP END CRITICAL (OPERNL4a_5)
               end do
             end if
!            End loop on ii=1,ilang2
           end do

           if ((choice==3.or.choice==23) .or. choice==6) then
!            Compute additional tensors related to strain gradients
!            ilang3 is number of unique tensor components of rank ilang+1
             ilang3=((ilang+2)*(ilang+3))/2
             jjs=((ilang+1)*(ilang+2)*(ilang+3))/6
!            Compute strain gradient tensor components
             do ii=1,ilang3
!              Note that iffkgs is also used by ddk and 2nd derivative parts
               iffkgs=iffkgs+1
               jj=ii+jjs
!              !$OMP CRITICAL (OPERNL4a_6)
               dgxds(1,jj-4,ia,iproj)=dgxds(1,jj-4,ia,iproj)+scalars(1,iffkgs)
               dgxds(2,jj-4,ia,iproj)=dgxds(2,jj-4,ia,iproj)+scalars(2,iffkgs)
!              !$OMP END CRITICAL (OPERNL4a_6)
             end do
           end if

           if (choice==6) then
!            Compute additional tensors related to strain 2nd derivatives
!            and internal strain derivatives
!            ilang6 is number of unique tensor components of rank ilang+3
             ilang6=((ilang+4)*(ilang+5))/2
             jjs=((ilang+3)*(ilang+4)*(ilang+5))/6
!            Compute strain gradient tensor components
             do ii=1,ilang6
!              Note that iffkgs is also used by ddk part
               iffkgs2=iffkgs2+1
               jj=ii+jjs
!              !$OMP CRITICAL (OPERNL4a_7)
               d2gxds2(1,jj-20,ia,iproj)=d2gxds2(1,jj-20,ia,iproj)+scalars(1,iffkgs2)
               d2gxds2(2,jj-20,ia,iproj)=d2gxds2(2,jj-20,ia,iproj)+scalars(2,iffkgs2)
!              !$OMP END CRITICAL (OPERNL4a_7)
             end do

!            ilang4 is number of unique tensor components of rank ilang
             ilang4=((ilang+1)*(ilang+2))/2
             jjs=((ilang)*(ilang+1)*(ilang+2))/6
!            Compute internal strain gradient tensor components
             do ii=1,ilang4
               iffkgs2=iffkgs2+1
               jj=ii+jjs
!              !$OMP CRITICAL (OPERNL4a_8)
!              Pay attention to the use of reals and imaginary parts here ...
               dgxdis(1,jj-1,ia,iproj)=dgxdis(1,jj-1,ia,iproj)-two_pi*scalars(2,iffkgs2)
               dgxdis(2,jj-1,ia,iproj)=dgxdis(2,jj-1,ia,iproj)+two_pi*scalars(1,iffkgs2)
!              !$OMP END CRITICAL (OPERNL4a_8)
             end do

!            ilang5 is number of unique tensor components of rank ilang+2
             ilang5=((ilang+3)*(ilang+4))/2
             jjs=((ilang+2)*(ilang+3)*(ilang+4))/6
!            Compute internal strain gradient tensor components
             do ii=1,ilang5
               iffkgs2=iffkgs2+1
               jj=ii+jjs
!              !$OMP CRITICAL (OPERNL4a_9)
!              Pay attention to the use of reals and imaginary parts here ...
               d2gxdis(1,jj-10,ia,iproj)=d2gxdis(1,jj-10,ia,iproj)-two_pi*scalars(2,iffkgs2)
               d2gxdis(2,jj-10,ia,iproj)=d2gxdis(2,jj-10,ia,iproj)+two_pi*scalars(1,iffkgs2)
!              !$OMP END CRITICAL (OPERNL4a_9)
             end do
           end if ! choice==6

           if (choice==5) then
!            Compute additional tensors related to ddk with ffnl(:,2,...)
             ilangx=(ilang*(ilang+1))/2
             jjs=((ilang-1)*ilang*(ilang+1))/6
             do ii=1,ilangx
!              Note that iffkgs is also used by strain part
               iffkgs=iffkgs+1
               jj=ii+jjs
!              !$OMP CRITICAL (OPERNL4a_10)
               dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)+scalars(1,iffkgs)
               dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+scalars(2,iffkgs)
!              !$OMP END CRITICAL (OPERNL4a_10)
             end do
!            Compute additional tensors related to ddk with ffnl(:,1,...)
             if(ilang>=2)then
               ilangx=((ilang-1)*ilang)/2
               jjs=((ilang-2)*(ilang-1)*ilang)/6
               do ii=1,ilangx
                 iffkgk=iffkgk+1
                 jj=ii+jjs
!                !$OMP CRITICAL (OPERNL4a_11)
                 dgxdt(1,2,jj,ia,iproj)=dgxdt(1,2,jj,ia,iproj)+scalars(1,iffkgk)
                 dgxdt(2,2,jj,ia,iproj)=dgxdt(2,2,jj,ia,iproj)+scalars(2,iffkgk)
!                !$OMP END CRITICAL (OPERNL4a_11)
               end do
             end if
           end if

!          End projector loop
         end do

!        End condition of non-zero projectors
       end if

!      End angular momentum loop
     end do

!    End loop on atoms
   end do

!  End loop on blocks of planewaves
 end do
!!$OMP END DO

 ABI_DEALLOCATE(ffkg)
 ABI_DEALLOCATE(kpgx)
 ABI_DEALLOCATE(parity)
 ABI_DEALLOCATE(scalars)
 ABI_DEALLOCATE(teffv)
!!$OMP END PARALLEL

!DEBUG
!write(std_out,*)' opernl4a : exit'
!ENDDEBUG

end subroutine opernl4a
!!***

!!****f* ABINIT/opernl4b
!! NAME
!! opernl4b
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! from projected quantities to reciprocal space
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

subroutine opernl4b(choice,dgxds,dgxdt,ffnl,gmet,gxa,&
&  ia3,idir,indlmn,ispinor,itypat,jproj,kg_k,kpg_k,kpt,&
&  lmnmax,matblk,mincat,mlang3,mlang4,mproj,ndgxdt,nffnl,nincat,&
&  nkpg,nlang,nloalg,npw,ntypat,ph3d,vect)

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
 real(dp),intent(inout) :: vect(:,:) !vz_i

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
!                see Eq.(53) of PRB55,10337(1997) [[cite:Gonze1997]]
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

end module m_opernl
!!***
