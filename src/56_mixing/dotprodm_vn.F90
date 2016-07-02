!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprodm_vn
!! NAME
!! dotprodm_vn
!!
!!
!! FUNCTION
!! For a set of densities and a set of potentials,
!! compute the dot product (integral over FFT grid) of each pair, to obtain
!! a series of energy-like quantity (so the usual dotproduct is divided
!! by the number of FFT points, and multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum correctly over them. Note that the storage of densities and
!! potentials is different : for potential, one stores the matrix components,
!! while for the density, one stores the trace, and then, either the
!! spin-polarisation (if nspden=2), or the magnetization vector (if nspden=4).
!! Need the index of the first density/potential pair to be treated, in each array,
!! and the number of pairs to be treated.
!! Might be used to compute just one dot product, in
!! a big array, such as to avoid copying the density and potential from a big array
!! to a temporary place.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex (not coded yet for nspden=4)
!!  denarr(cplex*nfft,nspden,nden)=real space density on FFT grid
!!  id=index of the first density to be treated in the denarr array
!!  ip=index of the first potential to be treated in the potarr array
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  multd=number of densities to be treated
!!  multp=number of potentials to be treated
!!  nden=third dimension of the denarr array
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  npot=third dimension of the potarr array
!!  nspden=number of spin-density components
!!  potarr(cplex*nfft,nspden,npot)=real space potential on FFT grid
!!                 (will be complex conjugated if cplex=2 and cpldot=2)
!!  ucvol=unit cell volume (Bohr**3)
!!
!! OUTPUT
!!  dot(cpldot,multp,multd)= series of values of the dot product potential/density
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Concerning storage when nspden=4:
!!   cplex=1:
!!     V are stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     N are stored as : n, m_x, m_y, m_z               (real)
!!   cplex=2:
!!     V are stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     N are stored as : n, m_x, m_y, mZ          (complex)
!!
!! PARENTS
!!      aprxdr
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dotprodm_vn(cplex,cpldot,denarr,dot,id,ip,mpicomm, mpi_summarize,multd,multp,&
& nden,nfft,nfftot,npot,nspden,potarr,ucvol)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprodm_vn'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpldot,cplex,id,ip,multd,multp,nden,nfft,nfftot,npot
 integer,intent(in) :: nspden,mpicomm
 logical, intent(in) :: mpi_summarize
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: denarr(cplex*nfft,nspden,nden)
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: dot(cpldot,multp,multd)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ir,jr
 real(dp) :: ai,ar,dim11,dim12,dim21,dim22,dim_dn,dim_up,dre11,dre12,dre21
 real(dp) :: dre22,dre_dn,dre_up,factor,pim11,pim12,pim21,pim22,pim_dn,pim_up
 real(dp) :: pre11,pre12,pre21,pre22,pre_dn,pre_up
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")

!Real or complex outputs are coded
 DBG_CHECK(ANY(cpldot==(/1,2/)),"Wrong cpldot")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")

 DBG_CHECK(id >= 1,'Wrong id')
 DBG_CHECK(ip >= 1,'Wrong id')

 DBG_CHECK(multd >= 1,"wrong multd")
 DBG_CHECK(multp >= 1,"wrong multp")

 DBG_CHECK(nden-id-multd >=-1,'nden-id-multd')
 DBG_CHECK(npot-ip-multp >=-1,'npot-ip-multp')

 if(nspden==1)then

   if(cpldot==1 .or. cplex==1 )then

     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar + potarr(ir,1,ip+i1-1)*denarr(ir,1,id+i2-1)
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2

   else  ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

     do i2=1,multd
       do i1=1,multp
         ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr) SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
         do ir=1,nfft
           jr=2*ir
           ar=ar + potarr(jr-1,1,ip+i1-1)*denarr(jr-1,1,id+i2-1) &
&           + potarr(jr  ,1,ip+i1-1)*denarr(jr  ,1,id+i2-1)
           ai=ai + potarr(jr-1,1,ip+i1-1)*denarr(jr  ,1,id+i2-1) &
&           - potarr(jr  ,1,ip+i1-1)*denarr(jr-1,1,id+i2-1)
         end do
         dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       end do ! i1
     end do ! i2

   end if

 else if(nspden==2)then

   if(cpldot==1 .or. cplex==1 )then

     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar + potarr(ir,1,ip+i1-1)* denarr(ir,2,id+i2-1)               &       ! This is the spin up contribution
&          + potarr(ir,2,ip+i1-1)*(denarr(ir,1,id+i2-1)-denarr(ir,2,id+i2-1)) ! This is the spin down contribution
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2

   else ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

     do i2=1,multd
       do i1=1,multp
         ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre_up,dim_up,dre_dn,dim_dn,pre_up,pim_up,pre_dn,pim_dn) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
         do ir=1,nfft
           jr=2*ir

           dre_up=denarr(jr-1,2,id+i2-1)
           dim_up=denarr(jr  ,2,id+i2-1)
           dre_dn=denarr(jr-1,1,id+i2-1)-dre_up
           dim_dn=denarr(jr  ,1,id+i2-1)-dim_up

           pre_up=potarr(jr-1,1,ip+i1-1)
           pim_up=potarr(jr  ,1,ip+i1-1)
           pre_dn=potarr(jr-1,2,ip+i1-1)
           pim_dn=potarr(jr  ,2,ip+i1-1)

           ar=ar + pre_up * dre_up &
&           + pim_up * dim_up &
&           + pre_dn * dre_dn &
&           + pim_dn * dim_dn
           ai=ai + pre_up * dim_up &
&           - pim_up * dre_up &
&           + pre_dn * dim_dn &
&           - pim_dn * dre_dn

         end do
         dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
       end do ! i1
     end do ! i2

   end if

 else if(nspden==4)then
!  \rho{\alpha,\beta} V^{\alpha,\beta} =
!  rho*(V^{11}+V^{22})/2$
!  + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2
   if (cplex==1) then
     do i2=1,multd
       do i1=1,multp
         ar=zero
!$OMP PARALLEL DO PRIVATE(ir) SHARED(id,i1,i2,ip,cplex,nfft,denarr,potarr) REDUCTION(+:ar)
         do ir=1,cplex*nfft
           ar=ar+(potarr(ir,1,ip+i1-1)+potarr(ir,2,ip+i1-1))*half*denarr(ir,1,id+i2-1)& ! This is the density contrib
&          + potarr(ir,3,ip+i1-1)                                *denarr(ir,2,id+i2-1)& ! This is the m_x contrib
&          - potarr(ir,4,ip+i1-1)                                *denarr(ir,3,id+i2-1)& ! This is the m_y contrib
&          +(potarr(ir,1,ip+i1-1)-potarr(ir,2,ip+i1-1))*half*denarr(ir,4,id+i2-1)       ! This is the m_z contrib
         end do
         dot(1,i1,i2)=ar
       end do ! i1
     end do ! i2
   else ! cplex=2
!    Note concerning storage when cplex=2:
!    V are stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!    N are stored as : n, m_x, m_y, mZ          (each are complex)
     if (cpldot==1) then
       do i2=1,multd
         do i1=1,multp
           ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre11,dim11,dre22,dim22,dre12,dim12,pre11,pim11,pre22,pim22,pre12,pim12) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar)
           do ir=1,nfft
             jr=2*ir
             dre11=half*(denarr(jr-1,1,id+i2)+denarr(jr-1,4,id+i2))
             dim11=half*(denarr(jr  ,1,id+i2)+denarr(jr-1,4,id+i2))
             dre22=half*(denarr(jr-1,1,id+i2)-denarr(jr-1,4,id+i2))
             dim22=half*(denarr(jr  ,1,id+i2)-denarr(jr-1,4,id+i2))
             dre12=half*(denarr(jr-1,2,id+i2)+denarr(jr  ,3,id+i2))
             dim12=half*(denarr(jr  ,2,id+i2)-denarr(jr-1,3,id+i2))
             dre21=half*(denarr(jr-1,2,id+i2)-denarr(jr  ,3,id+i2))
             dim21=half*(denarr(jr  ,2,id+i2)+denarr(jr-1,3,id+i2))
             pre11= potarr(jr-1,1,ip+i1)
             pim11= potarr(jr  ,1,ip+i1)
             pre22= potarr(jr-1,2,ip+i1)
             pim22= potarr(jr  ,2,ip+i1)
             pre12= potarr(jr-1,3,ip+i1)
             pim12= potarr(jr  ,3,ip+i1)
             pre21= potarr(jr  ,4,ip+i1)
             pim21=-potarr(jr-1,4,ip+i1)
             ar=ar + pre11 * dre11 &
&             + pim11 * dim11 &
&             + pre22 * dre22 &
&             + pim22 * dim22 &
&             + pre12 * dre12 &
&             + pim12 * dim12 &
&             + pre21 * dre21 &
&             + pim21 * dim21
           end do
           dot(1,i1,i2)=ar
         end do ! i1
       end do ! i2
     else !cpldot=2
       do i2=1,multd
         do i1=1,multp
           ar=zero ; ai=zero
!$OMP PARALLEL DO PRIVATE(ir,jr,dre11,dim11,dre22,dim22,dre12,dim12,pre11,pim11,pre12,pim12,pre22,pim22) &
!$OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
           do ir=1,nfft
             jr=2*ir
             dre11=half*(denarr(jr-1,1,id+i2)+denarr(jr-1,4,id+i2))
             dim11=half*(denarr(jr  ,1,id+i2)+denarr(jr-1,4,id+i2))
             dre22=half*(denarr(jr-1,1,id+i2)-denarr(jr-1,4,id+i2))
             dim22=half*(denarr(jr  ,1,id+i2)-denarr(jr-1,4,id+i2))
             dre12=half*(denarr(jr-1,2,id+i2)+denarr(jr  ,3,id+i2))
             dim12=half*(denarr(jr  ,2,id+i2)-denarr(jr-1,3,id+i2))
             dre21=half*(denarr(jr-1,2,id+i2)-denarr(jr  ,3,id+i2))
             dim21=half*(denarr(jr  ,2,id+i2)+denarr(jr-1,3,id+i2))
             pre11= potarr(jr-1,1,ip+i1)
             pim11= potarr(jr  ,1,ip+i1)
             pre22= potarr(jr-1,2,ip+i1)
             pim22= potarr(jr  ,2,ip+i1)
             pre12= potarr(jr-1,3,ip+i1)
             pim12= potarr(jr  ,3,ip+i1)
             pre21= potarr(jr  ,4,ip+i1)
             pim21=-potarr(jr-1,4,ip+i1)
             ar=ar + pre11 * dre11 &
&             + pim11 * dim11 &
&             + pre22 * dre22 &
&             + pim22 * dim22 &
&             + pre12 * dre12 &
&             + pim12 * dim12 &
&             + pre21 * dre21 &
&             + pim21 * dim21
             ai=ai + pre11 * dim11 &
&             - pim11 * dre11 &
&             + pre22 * dim22 &
&             - pim22 * dre22 &
&             + pre12 * dim12 &
&             - pim12 * dre12 &
&             + pre21 * dim21 &
&             - pim21 * dre21
           end do
           dot(1,i1,i2)=ar
           dot(2,i1,i2)=ai
         end do ! i1
       end do ! i2
     end if ! cpldot
   end if ! cplex
 end if ! nspden

 factor=ucvol/dble(nfftot)
 dot(:,:,:)=factor*dot(:,:,:)

!XG030513 : MPIWF reduction (addition) on dot is needed here
 if (mpi_summarize) then
   call timab(48,1,tsec)
   call xmpi_sum(dot,mpicomm ,ierr)
   call timab(48,2,tsec)
 end if

 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_vn
!!***
