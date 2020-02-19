!!****m* ABINIT/m_sigc
!! NAME
!!  m_gwrdm
!!
!! FUNCTION
!!  Compute density matrix correction for k-point 'k' using GW (imaginary freqs. are used in Sigma_c)
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

module m_gwrdm

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_time
 use m_wfd
 use m_dtset

 use defs_datatypes,   only : ebands_t
 use m_sigma,         only : sigma_t

 implicit none

 private
!!***

 public :: calc_rdm,calc_rdmc,natoccs,printdm1
!!***

contains
!!***

!!****f* ABINIT/calc_rdm
!! NAME
!! calc_rdm
!!
!! FUNCTION
!! Calculate density matrix corrections within k-point 'k' for vxc and Sigma_x
!!
!! INPUTS
!! b1gw, b2gw=min and max band indeces within ikpt for GW correction. 
!! dm1=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigmadriver.F90 for more details). 
!! pot=potential, Self-energy, or Self-energy-Potential difference, matrix size (i,j), where i and j belong to k.
!! b1=lowest band where the density matrix correction is applied

!! OUTPUT
!! 
!! Updated dm1 matrix with: Go.vxc.Go, Go.Sigma_x.Go, or Go.Sigma_x-Vxc.Go
!! PARENTS
!! CHILDREN
!! SOURCE

subroutine calc_rdm(ib1,ib2,kpoint,isgn,iinfo,pot,dm1,BSt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,kpoint,isgn,iinfo
 type(ebands_t),target,intent(in) :: BSt
 !arrays
 complex(dpc),intent(in) :: pot(:,:)
 complex(dpc),intent(inout) :: dm1(:,:)
!Local variables ------------------------------
!scalars
 character(len=500) :: msg,msg2
 integer :: ib1dm,ib2dm
 real(dp) :: sgn_spin,tol
!arrays

!************************************************************************

 DBG_ENTER("COLL")

 tol=1.0d-8
 if(isgn==0) then !Vxc
   sgn_spin=-2.0d0
   msg2='Vxc    '
 else if(isgn==1) then !Sigma_x
   sgn_spin=2.0d0       
   msg2='Sx     '
 else ! Both: Sigma_x - Vxc
   sgn_spin=2.0d0       
   msg2='Sx-Vxc '
 endif

 if(iinfo==0) then 
   write(msg,'(a37,a7,a14,3f10.5)')'Computing the 1-RDM correction for  ',msg2,' and k-point: ',BSt%kptns(1:,kpoint)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 endif

 do ib1dm=ib1,ib2-1  
   do ib2dm=ib1dm+1,ib2
     if((BSt%occ(ib1dm,kpoint,1)>tol) .and. (BSt%occ(ib2dm,kpoint,1)<tol)) then
       dm1(ib1dm,ib2dm)=sgn_spin*pot(ib1dm,ib2dm)/(BSt%eig(ib1dm,kpoint,1)-BSt%eig(ib2dm,kpoint,1))
       dm1(ib2dm,ib1dm)=dm1(ib1dm,ib2dm)
     endif
   enddo
 enddo

 DBG_EXIT("COLL")

end subroutine calc_rdm
!!***


subroutine calc_rdmc(ib1,ib2,nomega_sigc,kpoint,iinfo,Sr,weights,sigcme_k,BSt,dm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,nomega_sigc,kpoint,iinfo
 type(ebands_t),target,intent(in) :: BSt
 type(sigma_t) :: Sr
!arrays
 real(dp),intent(in) :: weights(:)
 complex(dpc),intent(inout) :: dm1(:,:)
 complex(dpc),intent(in) :: sigcme_k(:,:,:,:)
!Local variables ------------------------------
!scalars
 real(dp), parameter :: pi=3.141592653589793238462643383279502884197
 character(len=500) :: msg,msg2
 integer :: ib1dm,ib2dm,iquad 
 real(dp) :: tol
 complex(dpc) :: denominator,fact,division,dm1_mel
!arrays

!************************************************************************

 DBG_ENTER("COLL")

 tol=1.0d-8
 msg2='Sc     '
 fact=cmplx(0.0d0,-1.0d0/pi)

 if(iinfo==0) then
   write(msg,'(a37,a7,a14,3f10.5)')'Computing the 1-RDM correction for  ',msg2,' and k-point: ',BSt%kptns(1:,kpoint)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 endif

 dm1=0.0d0
 do ib1dm=ib1,ib2  
   do ib2dm=ib1dm,ib2 
     dm1_mel=0.0d0
     do iquad=1,nomega_sigc
       denominator=(Sr%omega_i(iquad)-BSt%eig(ib1dm,kpoint,1))*(Sr%omega_i(iquad)-BSt%eig(ib2dm,kpoint,1)) ! As in FHI-aims for RPA
       if(abs(denominator)>tol) then 
         if(abs(sigcme_k(iquad,ib1dm,ib2dm,1))>tol) then
           division=sigcme_k(iquad,ib1dm,ib2dm,1)/denominator 
           dm1_mel=dm1_mel+weights(iquad)*division
         endif 
         if(abs(sigcme_k(iquad,ib2dm,ib1dm,1))>tol) then  
           division=sigcme_k(iquad,ib2dm,ib1dm,1)/denominator 
           dm1_mel=dm1_mel+weights(iquad)*conjg(division)
         endif 
       endif
     enddo
     dm1_mel=fact*dm1_mel
     dm1(ib1dm,ib2dm)=dm1_mel
     dm1(ib2dm,ib1dm)=dm1(ib1dm,ib2dm)
   enddo  
 enddo

 DBG_EXIT("COLL")

end subroutine calc_rdmc
!!***

subroutine natoccs(ib1,ib2,dm1,BSt,kpoint,iinfo)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,kpoint,iinfo
 type(ebands_t),target,intent(in) :: BSt
!arrays
 complex(dpc),intent(inout) :: dm1(:,:)
!Local variables ------------------------------
!scalars
 integer::ndim,ib1dm,ib2dm,lwork,info
 character(len=500) :: msg,msg2
 real(dp) :: toccs_k
!arrays
 real(dp),allocatable :: occs(:),rwork(:)
 complex(dpc),allocatable :: work(:),dm1_tmp(:,:)
!************************************************************************

 DBG_ENTER("COLL")
 
 ndim=ib2-ib1+1
 lwork=2*ndim-1
 ABI_MALLOC(occs,(ndim))
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(dm1_tmp,(ndim,ndim))
 ABI_MALLOC(rwork,(3*ndim-2))

 dm1_tmp=0.0d0
 do ib2dm=1,ndim
   do ib1dm=ib2dm,ndim
     dm1_tmp(ib1dm,ib2dm)=dm1(ib1+(ib1dm-1),ib1+(ib2dm-1))
     dm1_tmp(ib2dm,ib1dm)=dm1_tmp(ib1dm,ib2dm)
   enddo
 enddo

 work=0.0d0
 occs=0.0d0
 info=0
 call zheev('v','u',ndim,dm1_tmp,ndim,occs,work,lwork,rwork,info)

 if(info==0) then
   if(iinfo==0) then       
     write(msg,'(a51,3f10.5)') 'Occs. after updating with the exchange at k-point:',BSt%kptns(1:,kpoint)
   else
     write(msg,'(a51,3f10.5)') 'Occs. after updating with the  ex+cor  at k-point:',BSt%kptns(1:,kpoint)
   endif 
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a2,*(f10.5))') '  ',occs(1:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 else
   write(msg,'(a39,3f10.5)') 'Error computing the occs. for k-point: ',BSt%kptns(1:,kpoint)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 endif
 
 toccs_k=0.0d0
 do ib1dm=1,ndim
   toccs_k=toccs_k+occs(ib1dm)
 enddo

 write(msg,'(a21,i5,a3,i5,a21,f10.5)') 'Total occ. from band ',ib1,' to', ib2,' at current k-point: ',toccs_k
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a5)') ' '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ABI_FREE(rwork)
 ABI_FREE(occs)
 ABI_FREE(work)
 ABI_FREE(dm1_tmp)

 DBG_EXIT("COLL")

end subroutine natoccs
!!***

subroutine printdm1(ib1,ib2,dm1) ! Basically used for debug
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2
!arrays
 complex(dpc),intent(in) :: dm1(:,:)
!Local variables ------------------------------
!scalars
 integer::ib1dm
 character(len=500) :: msg
!arrays
!************************************************************************
 do ib1dm=ib1,ib2
   write(msg,'(a2,*(f10.5))') '  ',dm1(ib1dm,ib1:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 enddo

end subroutine printdm1
!!***

end module m_gwrdm
!!***
