!{\src2tex{textfont=tt}}
!!****f* ABINIT/meanvalue_g
!! NAME
!! meanvalue_g
!!
!! FUNCTION
!!  Compute the mean value of one wavefunction, in reciprocal space,
!!  for an operator that is real, diagonal in reciprocal space:
!!  <wf|op|wf>
!!  For the time being, only spin-independent operators are treated.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2017 ABINIT group (XG,BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent !)
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo,ar_im)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

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
&   '  When istwf_k/=1, nspinor must be 1,',ch10,&
&   '  however, nspinor=',nspinor,', and istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(use_ndo==1 .and. (istwf_k==2 .and.me_g0==1)) then
   write(message,'(a,a)') ch10,' use_ndo==1, not tested'
   MSG_BUG(message)
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
