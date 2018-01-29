!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_vn
!! NAME
!! dotprod_vn
!!
!!
!! FUNCTION
!! Compute dot product of potential and density (integral over FFT grid), to obtain
!! an energy-like quantity (so the usual dotproduct is divided
!! by the number of FFT points, and multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum correctly over them. Note that the storage of densities and
!! potentials is different : for potential, one stores the matrix components,
!! while for the density, one stores the trace, and then, either the
!! up component (if nspden=2), or the magnetization vector (if nspden=4).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  dens(cplex*nfft,nspden)=real space density on FFT grid
!!  mpi_enreg=informations about MPI parallelization
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  nspden=number of spin-density components
!!  option= if 1, only the real part is computed
!!          if 2, both real and imaginary parts are computed  (not yet coded)
!!  pot(cplex*nfft,nspden)=real space potential on FFT grid
!!                 (will be complex conjugated if cplex=2 and option=2)
!!  ucvol=unit cell volume (Bohr**3)
!!
!! OUTPUT
!!  doti= imaginary part of the dot product, output only if option=2 (and cplex=2).
!!  dotr= real part
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Concerning storage when nspden=4:
!!   cplex=1:
!!     V is stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     N is stored as : n, m_x, m_y, m_z               (real)
!!   cplex=2:
!!     V is stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     N is stored as : n, m_x, m_y, mz          (complex)
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_eltfrxc,dfpt_nselt,dfpt_nstdy,dfpt_nstpaw,dfpt_rhotov
!!      dfptnl_loop,energy,newfermie1,odamix,prcrskerker2,rhotov,rhotoxc,setvtr
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine dotprod_vn(cplex,dens,dotr,doti,nfft,nfftot,nspden,option,pot,&
& ucvol,mpi_comm_sphgrid)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_xmpi, only: xmpi_comm_size,xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotprod_vn'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nfftot,nspden,option
 integer,intent(in),optional :: mpi_comm_sphgrid
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: doti,dotr
!arrays
 real(dp),intent(in) :: dens(cplex*nfft,nspden),pot(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer  :: ierr,ifft,jfft
 real(dp) :: dim11,dim12,dim21,dim22,dim_dn,dim_up,dre11,dre12,dre21,dre22
 real(dp) :: dre_dn,dre_up,factor,nproc_sphgrid,pim11,pim12,pim21,pim22,pim_dn,pim_up,pre11
 real(dp) :: pre12,pre21,pre22,pre_dn,pre_up
 real(dp) :: bx_re,bx_im,by_re,by_im,bz_re,bz_im,v0_re,v0_im
!arrays
 real(dp) :: buffer2(2)

! *************************************************************************

!Real or complex inputs are coded
 DBG_CHECK(ANY(cplex==(/1,2/)),"Wrong cplex")
 DBG_CHECK(ANY(nspden==(/1,2,4/)),"Wrong nspden")

!Real or complex output are coded
 DBG_CHECK(ANY(option==(/1,2/)),"Wrong option")

 dotr=zero; doti=zero

 if(nspden==1)then

   if(option==1 .or. cplex==1 )then
!$OMP PARALLEL DO REDUCTION(+:dotr)
     do ifft=1,cplex*nfft
       dotr=dotr + pot(ifft,1)*dens(ifft,1)
     end do
!    dotr = ddot(cplex*nfft,pot,1,dens,1)

   else  ! option==2 and cplex==2 : one builds the imaginary part, from complex den/pot

!$OMP PARALLEL DO PRIVATE(jfft) REDUCTION(+:dotr,doti)
     do ifft=1,nfft
       jfft=2*ifft
       dotr=dotr + pot(jfft-1,1)*dens(jfft-1,1) + pot(jfft,1)*dens(jfft  ,1)
       doti=doti + pot(jfft-1,1)*dens(jfft  ,1) - pot(jfft,1)*dens(jfft-1,1)
     end do

   end if

 else if(nspden==2)then

   if(option==1 .or. cplex==1 )then
!$OMP PARALLEL DO REDUCTION(+:dotr)
     do ifft=1,cplex*nfft
       dotr=dotr + pot(ifft,1)* dens(ifft,2)     &    ! This is the spin up contribution
&      + pot(ifft,2)*(dens(ifft,1)-dens(ifft,2))      ! This is the spin down contribution
     end do

   else ! option==2 and cplex==2 : one builds the imaginary part, from complex den/pot

!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nfft,dens,pot) REDUCTION(+:dotr,doti)
     do ifft=1,nfft

       jfft=2*ifft
       dre_up=dens(jfft-1,2)
       dim_up=dens(jfft  ,2)
       dre_dn=dens(jfft-1,1)-dre_up
       dim_dn=dens(jfft  ,1)-dim_up
       pre_up=pot(jfft-1,1)
       pim_up=pot(jfft  ,1)
       pre_dn=pot(jfft-1,2)
       pim_dn=pot(jfft  ,2)

       dotr=dotr + pre_up * dre_up &
&       + pim_up * dim_up &
&       + pre_dn * dre_dn &
&       + pim_dn * dim_dn
       doti=doti + pre_up * dim_up &
&       - pim_up * dre_up &
&       + pre_dn * dim_dn &
&       - pim_dn * dre_dn

     end do
   end if

 else if(nspden==4)then
!  \rho{\alpha,\beta} V^{\alpha,\beta} =
!  rho*(V^{11}+V^{22})/2$
!  + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2
   if (cplex==1) then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(nfft,dens,pot) REDUCTION(+:dotr)
     do ifft=1,nfft
       dotr=dotr + &
&       (pot(ifft,1) + pot(ifft,2))*half*dens(ifft,1) &   ! This is the density contrib
&      + pot(ifft,3)      *dens(ifft,2) &   ! This is the m_x contrib
&      - pot(ifft,4)      *dens(ifft,3) &   ! This is the m_y contrib
&      +(pot(ifft,1) - pot(ifft,2))*half*dens(ifft,4)     ! This is the m_z contrib
     end do
   else ! cplex=2
!    Note concerning storage when cplex=2:
!    V is stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!    N is stored as : n, m_x, m_y, mZ          (each are complex)
     if (option==1) then
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nfft,dens,pot) REDUCTION(+:dotr)
       do ifft=1,nfft
         jfft=2*ifft
         dre11=half*(dens(jfft-1,1)+dens(jfft-1,4))
         dim11=half*(dens(jfft  ,1)+dens(jfft-1,4))
         dre22=half*(dens(jfft-1,1)-dens(jfft-1,4))
         dim22=half*(dens(jfft  ,1)-dens(jfft-1,4))
         dre12=half*(dens(jfft-1,2)+dens(jfft  ,3))
         dim12=half*(dens(jfft  ,2)-dens(jfft-1,3))
         dre21=half*(dens(jfft-1,2)-dens(jfft  ,3))
         dim21=half*(dens(jfft  ,2)+dens(jfft-1,3))
         pre11= pot(jfft-1,1)
         pim11= pot(jfft  ,1)
         pre22= pot(jfft-1,2)
         pim22= pot(jfft  ,2)
         pre12= pot(jfft-1,3)
         pim12= pot(jfft  ,3)
         pre21= pot(jfft  ,4)
         pim21=-pot(jfft-1,4)

         v0_re=half*(pre11+pre22)
         v0_im=half*(pim11+pim22)
         bx_re=half*(pre12+pre21)
         bx_im=half*(pim12+pim21)
         by_re=half*(-pim12+pim21)
         by_im=half*(pre12-pre21)
         bz_re=half*(pre11-pre22)
         bz_im=half*(pim11-pim22)
         dotr=dotr+v0_re * dens(jfft-1,1)&
&         + v0_im * dens(jfft  ,1) &
&         + bx_re * dens(jfft-1,2) &
&         + bx_im * dens(jfft  ,2) &
&         + by_re * dens(jfft-1,3) &
&         + by_im * dens(jfft  ,3) &
&         + bz_re * dens(jfft-1,4) &
&         + bz_im * dens(jfft  ,4)
!         dotr=dotr + pre11 * dre11 &
!&         + pim11 * dim11 &
!&         + pre22 * dre22 &
!&         + pim22 * dim22 &
!&         + pre12 * dre12 &
!&         + pim12 * dim12 &
!&         + pre21 * dre21 &
!&         + pim21 * dim21
       end do
     else ! option=2
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nfft,dens,pot) REDUCTION(+:dotr,doti)
       do ifft=1,nfft
         jfft=2*ifft
         dre11=half*(dens(jfft-1,1)+dens(jfft-1,4))
         dim11=half*(dens(jfft  ,1)+dens(jfft-1,4))
         dre22=half*(dens(jfft-1,1)-dens(jfft-1,4))
         dim22=half*(dens(jfft  ,1)-dens(jfft-1,4))
         dre12=half*(dens(jfft-1,2)+dens(jfft  ,3))
         dim12=half*(dens(jfft  ,2)-dens(jfft-1,3))
         dre21=half*(dens(jfft-1,2)-dens(jfft  ,3))
         dim21=half*(dens(jfft  ,2)+dens(jfft-1,3))
         pre11= pot(jfft-1,1)
         pim11= pot(jfft  ,1)
         pre22= pot(jfft-1,2)
         pim22= pot(jfft  ,2)
         pre12= pot(jfft-1,3)
         pim12= pot(jfft  ,3)
         pre21= pot(jfft  ,4)
         pim21=-pot(jfft-1,4)
         dotr=dotr + pre11 * dre11 &
&         + pim11 * dim11 &
&         + pre22 * dre22 &
&         + pim22 * dim22 &
&         + pre12 * dre12 &
&         + pim12 * dim12 &
&         + pre21 * dre21 &
&         + pim21 * dim21
         doti=doti + pre11 * dim11 &
&         - pim11 * dre11 &
&         + pre22 * dim22 &
&         - pim22 * dre22 &
&         + pre12 * dim12 &
&         - pim12 * dre12 &
&         + pre21 * dim21 &
&         - pim21 * dre21
       end do
     end if ! option
   end if ! cplex
 end if ! nspden

 factor=ucvol/dble(nfftot)
 dotr=factor*dotr
 doti=factor*doti

!MPIWF reduction (addition) on dotr, doti is needed here
 if(present(mpi_comm_sphgrid)) then
   nproc_sphgrid=xmpi_comm_size(mpi_comm_sphgrid)
   if(nproc_sphgrid>1) then
     buffer2(1)=dotr
     buffer2(2)=doti
     call xmpi_sum(buffer2,mpi_comm_sphgrid,ierr)
     dotr=buffer2(1)
     doti=buffer2(2)
   end if
 end if

end subroutine dotprod_vn
!!***
