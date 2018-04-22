!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spacepar
!! NAME
!! m_spacepar
!!
!! FUNCTION
!!  Relatively Low-level procedures operating on arrays defined on the FFT box (G- or R- space)
!!  Unlike the procedures in m_cgtools, the routines declared in this module can use mpi_type.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (XG, BA, MT)
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

module m_spacepar

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use defs_abitypes,     only : MPI_type
 use m_mpinfo,          only : ptabs_fourdp

 implicit none

 private
!!***

public :: meanvalue_g       ! Compute <wf|op|wf> where op is real and diagonal in G-space.
public :: laplacian         ! Compute the laplacian of a function defined in real space
!!***

contains
!!***

!!****f* ABINIT/meanvalue_g
!! NAME
!! meanvalue_g
!!
!! FUNCTION
!!  Compute the mean value of one wavefunction, in reciprocal space,
!!  for an operator that is real, diagonal in G-space: <wf|op|wf>
!!  For the time being, only spin-independent operators are treated.
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent!)
!!  filter= if 1, need to filter on the value of diag, that must be less than huge(0.0d0)*1.d-11
!!          otherwise, should be 0
!!  istwf_k=storage mode of the vectors
!!  npw=number of planewaves of the vector
!!  nspinor=number of spinor components
!!  vect(2,npw*nspinor)=vector
!!  vect1(2,npw*nspinor*use_ndo)=vector1 (=vector in most of the cases)
!!  use_ndo = says if vect=/vect1
!!
!! OUTPUT
!!  ar=mean value
!!
!! PARENTS
!!      dfpt_vtowfk,energy,forstrnps,vtowfk
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo,ar_im)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'meanvalue_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: filter,istwf_k,npw,nspinor,use_ndo
 real(dp),intent(out) :: ar
 real(dp),intent(out),optional :: ar_im
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: diag(npw),vect(2,npw*nspinor)
 real(dp),intent(in) :: vect1(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw,jpw,me_g0
 character(len=500) :: message
!arrays

! *************************************************************************
 me_g0 = mpi_enreg%me_g0

 DBG_CHECK(ANY(filter==(/0,1/)),"Wrong filter")
 DBG_CHECK(ANY(nspinor==(/1,2/)),"Wrong nspinor")
 DBG_CHECK(ANY(istwf_k==(/(ipw,ipw=1,9)/)),"Wrong istwf_k")

 if(nspinor==2 .and. istwf_k/=1)then
   write(message,'(a,a,a,i6,a,i6)')&
&   'When istwf_k/=1, nspinor must be 1,',ch10,&
&   'however, nspinor=',nspinor,', and istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(use_ndo==1 .and. (istwf_k==2 .and.me_g0==1)) then
   MSG_BUG('use_ndo==1, not tested')
 end if

 ar=zero
 if(present(ar_im)) ar_im=zero

!Normal storage mode
 if(istwf_k==1)then

!  No filter
   if(filter==0)then
!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end do
     end if
     if(use_ndo==1 .and. nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
       end do
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
       end do
     end if

!    !$OMP PARALLEL DO REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    end if
   else ! will filter

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         if(diag(jpw)<huge(0.0d0)*1.d-11)then
           ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
         end if
       end do
     end if
     if(use_ndo==1 .and. nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         if(diag(ipw)<huge(0.0d0)*1.d-11)then
           ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end if
       end do
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         if(diag(jpw)<huge(0.0d0)*1.d-11)then
           ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end if
       end do
     end if


!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    if(diag(ipw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    if(diag(ipw-npw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    end if ! nspinor==2

   end if ! filter==0

 else if(istwf_k>=2)then

   if(filter==0)then
     i1=1
     if(istwf_k==2 .and. me_g0==1)then ! MPIWF need to know which proc has G=0
       ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do

   else ! filter/=0
     i1=1
     if(istwf_k==2 .and. me_g0==1)then
       if(diag(1)<huge(0.0d0)*1.d-11)then
         ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
       end if
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
   end if ! filter==0

   ar=two*ar

 end if ! istwf_k

!MPIWF need to make reduction on ar and ai .
 if(mpi_enreg%paral_kgb==1)then
   call xmpi_sum(ar,mpi_enreg%comm_bandspinorfft ,ierr)
   if(present(ar_im))then
     call xmpi_sum(ar_im,mpi_enreg%comm_bandspinorfft,ierr)
   end if
 end if

end subroutine meanvalue_g
!!***

!!****f* m_spacepar/laplacian
!! NAME
!! laplacian
!!
!! FUNCTION
!! compute the laplacian of a function defined in real space
!! the code is written in the way of /3xc/xcden.F90
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=number of points of the fft grid
!!  nfunc=number of functions on the grid for which the laplacian is to be calculated
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  (optional) rdfuncr(nfft,nfunc)=real(dp) discretized functions in real space
!!  rdfuncg_in TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_in TO BE DESCRIBED SB 090901
!!  (optional) g2cart_in(nfft) = G**2 on the grid
!!
!! OUTPUT
!! (optional) laplacerdfuncr = laplacian in real space of the functions in rdfuncr
!! (optional) rdfuncg = real(dp) discretized functions in fourier space
!! (optional) laplacerdfuncg = real(dp) discretized laplacian of the functions in fourier space
!! (optional) g2cart_out(nfft) = G**2 on the grid
!!  rdfuncg_out TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_out TO BE DESCRIBED SB 090901
!!
!! PARENTS
!!      frskerker1,frskerker2,moddiel_csrb,prcrskerker1,prcrskerker2
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,paral_kgb,rdfuncr,&
&  laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'laplacian'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfunc,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
 real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 real(dp),intent(in),optional,target :: g2cart_in(nfft) !vz_i
 real(dp),intent(out),optional,target :: g2cart_out(nfft)  !vz_i
 real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
 real(dp),intent(in),optional,target :: rdfuncg_in(2,nfft,nfunc) !vz_i
 real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)

!Local variables-------------------------------
!scalars
 integer :: count,i1,i2,i3,id1,id2,id3,ifft,ifunc,ig1,ig2,ig3,ii1,n1,n2
 integer :: n3
 real(dp) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp),ABI_CONTIGUOUS pointer :: g2cart(:),laplacerdfuncg(:,:,:),rdfuncg(:,:,:)

! *************************************************************************

!Keep local copy of fft dimensions
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)

 if(present(laplacerdfuncg_out)) then
   laplacerdfuncg => laplacerdfuncg_out
 else
   ABI_ALLOCATE(laplacerdfuncg,(2,nfft,nfunc))
 end if

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!change the real density rdfuncr on real space on the real density
!rdfuncg in reciprocal space
 if(.not.present(rdfuncg_in)) then
   if(present(rdfuncg_out)) then
     rdfuncg => rdfuncg_out
   else
     ABI_ALLOCATE(rdfuncg,(2,nfft,nfunc))
   end if
   if(present(rdfuncr)) then
     do ifunc=1,nfunc
       call fourdp(1,rdfuncg(:,:,ifunc),rdfuncr(:,ifunc),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     end do
   end if
 else
   rdfuncg => rdfuncg_in
 end if

!apply the laplacian on laplacerdfuncr
!code from /3xc/xcden.F90
!see also src/5common/hatre.F90 and src/5common/moddiel.F90
!Keep local copy of fft dimensions
!Initialize computation of G^2 in cartesian coordinates
 if(.not.present(g2cart_in)) then
   if(present(g2cart_out)) then
     g2cart => g2cart_out
   else
     ABI_ALLOCATE(g2cart,(nfft))
   end if
   id1=int(n1/2)+2
   id2=int(n2/2)+2
   id3=int(n3/2)+2
   count=0
   do i3=1,n3
     ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
     ig3=i3-int(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-int(i2/id2)*n2-1

         ii1=1
         do i1=ii1,n1
           ig1=i1-int(i1/id1)*n1-1
           ifft=ifft+1

           b11=gprimd(1,1)*real(ig1,dp)
           b21=gprimd(2,1)*real(ig1,dp)
           b31=gprimd(3,1)*real(ig1,dp)
           b12=gprimd(1,2)*real(ig2,dp)
           b22=gprimd(2,2)*real(ig2,dp)
           b32=gprimd(3,2)*real(ig2,dp)
           b13=gprimd(1,3)*real(ig3,dp)
           b23=gprimd(2,3)*real(ig3,dp)
           b33=gprimd(3,3)*real(ig3,dp)

           g2cart(ifft)=( &
&           (b11+b12+b13)**2&
&           +(b21+b22+b23)**2&
&           +(b31+b32+b33)**2&
&           )
           do ifunc=1,nfunc
!            compute the laplacian in Fourier space that is * (i x 2pi x G)**2
             laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
             laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
           end do
         end do
       end if
     end do
   end do
   if(.not.present(g2cart_out))  then
     ABI_DEALLOCATE(g2cart)
   end if
 else
   g2cart => g2cart_in
   do ifunc=1,nfunc
     do ifft=1,nfft
!      compute the laplacian in Fourier space that is * (i x 2pi x G)**2
       laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
       laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
     end do
   end do
 end if

!get the result back into real space
 if(present(laplacerdfuncr)) then
   do ifunc=1,nfunc
     call fourdp(1,laplacerdfuncg(:,:,ifunc),laplacerdfuncr(:,ifunc),1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   end do
 end if

!deallocate pointers
 if((.not.present(rdfuncg_in)).and.(.not.present(rdfuncg_in)))  then
   ABI_DEALLOCATE(rdfuncg)
 end if
 if(.not.present(laplacerdfuncg_out))  then
   ABI_DEALLOCATE(laplacerdfuncg)
 end if

end subroutine laplacian
!!***
end module m_spacepar
!!***
