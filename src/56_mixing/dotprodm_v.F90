!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprodm_v
!! NAME
!! dotprodm_v
!!
!!
!! FUNCTION
!! For two sets of potentials,
!! compute dot product of each pair of two potentials (integral over FFT grid), to obtain
!! a series of square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden),
!! and sum over them.
!! Need the index of the first pair of potentials to be treated, in each array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in
!! a big array, such as to avoid copying a potential from a big array
!! to a temporary place.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex
!!  index1=index of the first potential to be treated in the potarr1 array
!!  index2=index of the first potential to be treated in the potarr2 array
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  mult1=number of potentials to be treated in the first set
!!  mult2=number of potentials to be treated in the second set
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot1= third dimension of the potarr1 array
!!  npot2= third dimension of the potarr2 array
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potentials are stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potentials are stored as V, B_x, B_y, Bz  (B=magn. field)
!!  potarr1(cplex*nfft,nspden,npot)=first array of real space potentials on FFT grid
!!    (if cplex=2 and cpldot=2, potarr1 is the array that will be complex conjugated)
!!  potarr2(cplex*nfft,nspden,npot)=second array of real space potentials on FFT grid
!!
!! OUTPUT
!!  dot(cpldot,mult1,mult2)= series of values of the dot product
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Concerning storage when nspden=4:
!!   cplex=1:
!!     opt_storage=0: V are stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     opt_storage=1: V are stored as : V, B_x, B_y, B_z               (real)
!!   cplex=2:
!!     opt_storage=0: V are stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     opt_storage=1: V are stored as : V, B_x, B_y, B_z         (complex)
!!
!! PARENTS
!!      scfopt
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dotprodm_v(cplex,cpldot,dot,index1,index2,mpicomm,mpi_summarize,&
&   mult1,mult2,nfft,npot1,npot2,nspden,opt_storage,potarr1,potarr2)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprodm_v'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpldot,cplex,index1,index2,mult1,mult2,nfft,npot1,npot2
 integer,intent(in) :: nspden,opt_storage,mpicomm
 logical, intent(in) :: mpi_summarize
!arrays
 real(dp),intent(in) :: potarr1(cplex*nfft,nspden,npot1)
 real(dp),intent(in) :: potarr2(cplex*nfft,nspden,npot2)
 real(dp),intent(out) :: dot(cpldot,mult1,mult2)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ifft,ispden
 real(dp) :: ai,ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")

!Real or complex outputs are coded
 DBG_CHECK(ANY(cpldot==(/1,2/)),"Wrong cpldot")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")
 DBG_CHECK( npot1-index1-mult1 >= -1,"npot1-index1-mult1")
 DBG_CHECK( npot2-index2-mult2 >= -1,"npot2-index2-mult2")

 if(cplex==1 .or. cpldot==1)then

   do i1=1,mult1
     do i2=1,mult2
       ar=zero
       do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
         do ifft=1,cplex*nfft
           ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
         end do
       end do
       dot(1,i1,i2)=ar
       if (nspden==4) then
         ar=zero
         do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
           do ifft=1,cplex*nfft
             ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
           end do
         end do
         if (opt_storage==0) then
           if (cplex==1) then
             dot(1,i1,i2)=dot(1,i1,i2)+two*ar
           else
             dot(1,i1,i2)=dot(1,i1,i2)+ar
           end if
         else
           dot(1,i1,i2)=half*(dot(1,i1,i2)+ar)
         end if
       end if
     end do
   end do

 else ! if (cplex==2 .and. cpldot==2)

   do i1=1,mult1
     do i2=1,mult2
       ar=zero ; ai=zero
       do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar,ai)
         do ifft=1,nfft
           ar=ar + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1) &
&           + potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1)
           ai=ai + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1) &
&           - potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1)
         end do
       end do
       dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       if (nspden==4) then
         ar=zero
         do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar,ai)
           do ifft=1,nfft
             ar=ar + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1) &
&             + potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1)
             ai=ai + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1) &
&             - potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1)
           end do
         end do
         if (opt_storage==0) then
           dot(1,i1,i2)=dot(1,i1,i2)+ar
           dot(2,i1,i2)=dot(2,i1,i2)+ai
         else
           dot(1,i1,i2)=half*(dot(1,i1,i2)+ar)
           dot(2,i1,i2)=half*(dot(2,i1,i2)+ai)
         end if
       end if
     end do
   end do
 end if

!XG030513 : MPIWF reduction (addition) on dot is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(dot,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_v
!!***
