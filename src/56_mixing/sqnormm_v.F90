!{\src2tex{textfont=tt}}
!!****f* ABINIT/sqnormm_v
!! NAME
!! sqnormm_v
!!
!! FUNCTION
!! For a series of potentials,
!! compute square of the norm (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden), and sum over them.
!! Need the index of the first potential to be treated, in the provided array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in a big array, such as to avoid 
!! copying a potential from a big array to a temporary place.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  index=index of the first potential to be treated
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  mult=number of potentials to be treated
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot= third dimension of the potarr array
!!  nspden=number of spin-density components
!!  opt_storage: 0, if potential is stored as V^up-up, V^dn-dn, Re[V^up-dn], Im[V^up-dn]
!!               1, if potential is stored as V, B_x, B_y, Bz  (B=magn. field)
!!  potarr(cplex*nfft,nspden,npot)=array of real space potentials on FFT grid
!!
!! OUTPUT
!!  norm2(mult)= value of the square of the norm of the different potentials
!!
!! PARENTS
!!      scfcge,scfopt
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sqnormm_v(cplex,index,mpicomm, mpi_summarize,mult,nfft,norm2,npot,nspden,opt_storage,potarr)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnormm_v'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,index,mult,nfft,npot,nspden,opt_storage,mpicomm
 logical, intent(in) :: mpi_summarize
!arrays
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: norm2(mult)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ii,ispden
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")

 DBG_CHECK(index>=1,"wrong index")
 DBG_CHECK(mult>=1,"wrong mult")
 DBG_CHECK(npot>=1,"wrong npot")

 DBG_CHECK(npot-index-mult>=-1,'npot-index-mult')

 do ii=1,mult
   ar=zero
   do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
       ar=ar + potarr(ifft,ispden,index+ii-1)**2
     end do
   end do
   norm2(ii)=ar
   if (nspden==4) then
     ar=zero
     do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
       do ifft=1,cplex*nfft
         ar=ar + potarr(ifft,ispden,index+ii-1)**2
       end do
     end do
     if (opt_storage==0) then
       if (cplex==1) then
         norm2(ii)=norm2(ii)+two*ar
       else
         norm2(ii)=norm2(ii)+ar
       end if
     else
       norm2(ii)=half*(norm2(ii)+ar)
     end if
   end if
 end do

!XG030513 : MPIWF reduction (addition) on norm2 is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(norm2,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

end subroutine sqnormm_v
!!***
