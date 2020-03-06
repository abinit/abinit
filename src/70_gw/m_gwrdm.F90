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
 use m_hdr
 use m_dtset

 use defs_datatypes,   only : ebands_t
 use m_sigma,          only : sigma_t
 
 implicit none

 private
!!***

 public :: calc_rdm,calc_rdmc,natoccs,printdm1,update_wfk_gw_rdm
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
 if(isgn==0) then ! Sigma_x - Vxc
   sgn_spin=2.0d0
   msg2='Sx-Vxc '
 else ! Vxc
   sgn_spin=-2.0d0       
   msg2='Vxc    '
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

subroutine natoccs(ib1,ib2,dm1,nateigv,occs_ks,BSt,kpoint,iinfo)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,kpoint,iinfo
 type(ebands_t),target,intent(in) :: BSt
!arrays
 real(dp),intent(inout) :: occs_ks(:,:)
 complex(dpc),intent(inout) :: dm1(:,:),nateigv(:,:,:)
!Local variables ------------------------------
!scalars
 integer::ndim,ib1dm,ib2dm,lwork,info
 character(len=500) :: msg,msg2
 real(dp) :: toccs_k
!arrays
 real(dp),allocatable :: occs(:),rwork(:),occs_tmp(:)
 complex(dpc),allocatable :: work(:),dm1_tmp(:,:),eigenvect(:,:)
!************************************************************************

 DBG_ENTER("COLL")
 
 ndim=ib2-ib1+1
 lwork=2*ndim-1
 ABI_MALLOC(occs,(ndim))
 ABI_MALLOC(occs_tmp,(ndim))
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(dm1_tmp,(ndim,ndim))
 ABI_MALLOC(eigenvect,(ndim,ndim))
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

 !Order from highes occ to lowest occ
 do ib1dm=1,ndim
  occs_tmp(ib1dm)=occs(ndim-(ib1dm-1))
  do ib2dm=1,ndim
   eigenvect(ib1dm,ib2dm)=dm1_tmp(ndim-(ib1dm-1),ib2dm)
  enddo
 enddo

 if(info==0) then
   if(iinfo==0) then       
     write(msg,'(a51,3f10.5)') 'Occs. after updating with the exchange at k-point:',BSt%kptns(1:,kpoint)
   else
     write(msg,'(a51,3f10.5)') 'Occs. after updating with the  ex+cor  at k-point:',BSt%kptns(1:,kpoint)
   endif 
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ib1dm=ndim-(ndim/10)*10
   do ib2dm=1,(ndim/10)*10,10
     write(msg,'(f11.5,9f10.5)') occs_tmp(ib2dm:ib2dm+9)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   enddo  
   ib1dm=(ndim/10)*10+1
   write(msg,'(f11.5,*(f10.5))') occs_tmp(ib1dm:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 else
   write(msg,'(a36,3f10.5)') 'Error computing occs. for k-point: ',BSt%kptns(1:,kpoint)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 endif

 ! Store natural orbital eigenvectors matrix and occs. Also compute total number of electrons for this k value
 toccs_k=0.0d0
 do ib1dm=1,ndim
   do ib2dm=1,ndim
     nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),kpoint)=eigenvect(ib2dm,ib1dm)
   enddo
   occs_ks(ib1+(ib1dm-1),kpoint)=occs_tmp(ib1dm)
   toccs_k=toccs_k+occs_tmp(ib1dm)
 enddo

 write(msg,'(a22,i5,a3,i5,a21,f10.5)') ' Total occ. from band ',ib1,' to', ib2,' at current k-point: ',toccs_k
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a5)') ' '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ABI_FREE(rwork)
 ABI_FREE(occs)
 ABI_FREE(work)
 ABI_FREE(dm1_tmp)
 ABI_FREE(eigenvect)
 ABI_FREE(occs_tmp)

 DBG_EXIT("COLL")

end subroutine natoccs
!!***

subroutine update_wfk_gw_rdm(wfd_i,wfd_f,nateigv,occs,b1gw,b2gw,BSt,Hdr,Hdr2)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: b1gw,b2gw
 type(wfd_t),intent(inout) :: wfd_f
 type(wfd_t),intent(in) :: wfd_i
 type(ebands_t),target,intent(inout) :: BSt
 type(Hdr_type),intent(inout) :: Hdr,Hdr2
!arrays
 real(dp),intent(in) :: occs(:,:)
 complex(dpc),intent(in) :: nateigv(:,:,:)
!Local variables ------------------------------
!scalars
 character(len=500) :: msg
 integer :: ib1dm,ib2dm,dim_bands,ikpoint,irecip_v
!arrays
!************************************************************************
 DBG_ENTER("COLL")

 !call xmpi_barrier(Wfd_i%comm)

 ! Only master
 !if(xmpi_comm_rank(Wfd_i%comm)==0) then 
 !Wfd%Wave(6,2,1)%ug(3) ! BAND 6 , k-POINT 2, SPIN 1, UG="MO Coef"= 3
   do ikpoint=1,BSt%nkpt
     write(msg,'(a31,i5,a9,i5,a1,i5,a7,i5)') ' Calc. nat. orbs coefs k-point: ',ikpoint,', bands: ',b1gw,'-',b2gw,', npw: ',wfd_i%Kdata(ikpoint)%npw ! Print for debug
     call wrtout(std_out,msg,'COLL')
     do ib1dm=b1gw,b2gw
       do irecip_v=1,wfd_i%Kdata(ikpoint)%npw ! No spin used, setting nspinor=1 
        wfd_f%Wave(ib1dm,ikpoint,1)%ug(irecip_v)=0.0d0 
        do ib2dm=b1gw,b2gw
          wfd_f%Wave(ib1dm,ikpoint,1)%ug(irecip_v)=wfd_f%Wave(ib1dm,ikpoint,1)%ug(irecip_v) &
                  +nateigv(ib2dm,ib1dm,ikpoint)*wfd_i%Wave(ib2dm,ikpoint,1)%ug(irecip_v)
        enddo
       enddo
       !write(*,'(*(f10.5))') real(nateigv(ib1dm,:,ikpoint)) !Print for debug
     enddo
     !write(*,*) ' ' !Space for debug
   enddo
 !endif
 
 !call xmpi_barrier(Wfd_i%comm)

 ! MRM BSt occ are changed and never recoverd
 MSG_COMMENT("QP_BSt occupations were updated with nat. orb. ones")
 do ikpoint=1,BSt%nkpt
   BSt%occ(b1gw:b2gw,ikpoint,1) = occs(b1gw:b2gw,ikpoint) ! No spin used
 enddo
 if((size(Hdr%occ(:))/BSt%nkpt) < (b2gw-b1gw+1)) then
   !Actually, we should never reach this point as the code should crash during Wfd initialization in m_sigma_driver.F90
   MSG_ERROR("Impossible to use the existing read WFK to build a new one!")
 endif
 ! MRM change occs in Header Hdr
 ! Update occ in Hdr before printing
 MSG_COMMENT("Hdr and Hdr_sigma occupations were updated with nat. orb. ones")
 ib1dm=1
 do ikpoint=1,BSt%nkpt
   dim_bands=size(BSt%occ(:,ikpoint,1))
   do ib2dm=1,dim_bands   
     Hdr2%occ(ib1dm)=BSt%occ(ib2dm,ikpoint,1)
     Hdr%occ(ib1dm)=BSt%occ(ib2dm,ikpoint,1)
     ib1dm=ib1dm+1
   enddo
 enddo

 !call xmpi_barrier(Wfd_i%comm)

end subroutine update_wfk_gw_rdm        
!!***

subroutine printdm1(ib1,ib2,dm1) ! Only used for debug of this file, do not use it with large arrays
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
