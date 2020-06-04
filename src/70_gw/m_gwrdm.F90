!!****m* ABINIT/m_gwrdm
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

 public :: calc_rdmx,calc_rdmc,natoccs,printdm1,update_hdr_bst,rotate_exchange ! ,update_wfd_bst ! -> The commented one only works on serial mode
!!***

contains
!!***

!!****f* ABINIT/calc_rdmx
!! NAME
!! calc_rdmx
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
!! EX_WIFE?
!! SOURCE

subroutine calc_rdmx(ib1,ib2,kpoint,isgn,iinfo,pot,dm1,BSt)

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
 real(dp) :: sgn_spin
!arrays

!************************************************************************

 DBG_ENTER("COLL")

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
     if((BSt%occ(ib1dm,kpoint,1)>tol8) .and. (BSt%occ(ib2dm,kpoint,1)<tol8)) then
       dm1(ib1dm,ib2dm)=sgn_spin*pot(ib1dm,ib2dm)/(BSt%eig(ib1dm,kpoint,1)-BSt%eig(ib2dm,kpoint,1))
       dm1(ib2dm,ib1dm)=dm1(ib1dm,ib2dm)
     endif
   enddo
 enddo

 DBG_EXIT("COLL")

end subroutine calc_rdmx
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
 complex(dpc) :: denominator,fact,division,dm1_mel
!arrays

!************************************************************************

 DBG_ENTER("COLL")

 msg2='Sc     '
 fact=cmplx(1.0d0/pi,0.0d0)

 if(iinfo==0) then
   write(msg,'(a37,a7,a14,3f10.5)')'Computing the 1-RDM correction for  ',msg2,' and k-point: ',BSt%kptns(1:,kpoint)
   call wrtout(std_out,msg,'COLL')
   !call wrtout(ab_out,msg,'COLL')
   write(msg,'(a11,i5,a8,i5,a5,i5,a7)')'from band ',ib1,' to band',ib2,' with',nomega_sigc,' omegas'
   call wrtout(std_out,msg,'COLL')
   !call wrtout(ab_out,msg,'COLL')
   write(msg,'(a11,5x,*(f17.5))')'First omega',Sr%omega_i(1)
   call wrtout(std_out,msg,'COLL')
   !call wrtout(ab_out,msg,'COLL')
   write(msg,'(a11,5x,*(f17.5))')'Last  omega',Sr%omega_i(nomega_sigc)
   call wrtout(std_out,msg,'COLL')
   !call wrtout(ab_out,msg,'COLL')
 endif

 dm1=0.0d0
 do ib1dm=ib1,ib2  
   do ib2dm=ib1dm,ib2 
     dm1_mel=0.0d0
     do iquad=1,nomega_sigc
       denominator=(Sr%omega_i(iquad)-BSt%eig(ib1dm,kpoint,1))*(Sr%omega_i(iquad)-BSt%eig(ib2dm,kpoint,1)) ! As in FHI-aims for RPA
       if(abs(denominator)>tol8) then 
         ! Sigma_pq/[(denominator)]
           division=sigcme_k(iquad,ib1dm,ib2dm,1)/denominator 
           dm1_mel=dm1_mel+weights(iquad)*division        ! +
         ! [Sigma_qp/[(denominator)]]^*
           division=sigcme_k(iquad,ib2dm,ib1dm,1)/denominator 
           dm1_mel=dm1_mel+weights(iquad)*conjg(division) ! +
       endif
     enddo
     dm1(ib1dm,ib2dm)=fact*dm1_mel
     ! Dji = Dij^*
     dm1(ib2dm,ib1dm)=conjg(dm1(ib1dm,ib2dm))
   enddo  
 enddo

 DBG_EXIT("COLL")

end subroutine calc_rdmc
!!***

subroutine natoccs(ib1,ib2,dm1,nateigv,occs_ks,BSt,ikpoint,iinfo)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ikpoint,iinfo
 type(ebands_t),target,intent(in) :: BSt
!arrays
 real(dp),intent(inout) :: occs_ks(:,:)
 complex(dpc),intent(inout) :: dm1(:,:),nateigv(:,:,:,:)
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

 !call printdm1(1,ndim,dm1_tmp) ! Uncomment for debug 
 !eigenvect=dm1_tmp
 !Order from highest occ to lowest occ
 occs_tmp=occs
 do ib1dm=1,ndim
  occs_tmp(ib1dm)=occs(ndim-(ib1dm-1))
  do ib2dm=1,ndim
   eigenvect(ib2dm,ib1dm)=dm1_tmp(ib2dm,(ndim-(ib1dm-1)))
  enddo
 enddo

 if(info==0) then
   if(iinfo==0) then       
     write(msg,'(a51,3f10.5)') 'Occs. after updating with Sx-Vxc corr. at k-point:',BSt%kptns(1:,ikpoint)
   else
     write(msg,'(a51,3f10.5)') 'Occs. after updating with S_c correct. at k-point:',BSt%kptns(1:,ikpoint)
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
   write(msg,'(a36,3f10.5)') 'Error computing occs. for k-point: ',BSt%kptns(1:,ikpoint)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 endif

 ! Store natural orbital eigenvectors matrix and occs. Also compute total number of electrons for this k value
 toccs_k=0.0d0
 do ib1dm=1,ndim
   do ib2dm=1,ndim
     nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),ikpoint,1)=eigenvect(ib1dm,ib2dm)
   enddo
   occs_ks(ib1+(ib1dm-1),ikpoint)=occs_tmp(ib1dm)
   toccs_k=toccs_k+occs_tmp(ib1dm)
 enddo
 !call printdm1(1,ndim,eigenvect) ! Uncomment for debug 

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

subroutine update_hdr_bst(Wfd,occs,b1gw,b2gw,BSt,Hdr,ngfft_in)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: b1gw,b2gw
 integer,intent(in),dimension(3) :: ngfft_in
 type(ebands_t),target,intent(inout) :: BSt
 type(Hdr_type),intent(inout) :: Hdr
 type(wfd_t),intent(in) :: Wfd
 !!type(),
!arrays
 real(dp),intent(in) :: occs(:,:)
!Local variables ------------------------------
!scalars
 character(len=500) :: msg
 integer :: ib1dm,ib2dm,dim_bands,ikpoint
!arrays
!************************************************************************
 DBG_ENTER("COLL")

 ! MRM BSt occ are changed and never recoverd
 do ikpoint=1,BSt%nkpt
   BSt%occ(b1gw:b2gw,ikpoint,1) = occs(b1gw:b2gw,ikpoint) ! No spin used
 enddo
 MSG_COMMENT("QP_BSt occupations correctly updated with nat. orb. ones")
 if((size(Hdr%occ(:))/BSt%nkpt) < (b2gw-b1gw+1)) then
   !Actually, we should never reach this point because the code should stop during Wfd initialization in m_sigma_driver.F90
   MSG_ERROR("Impossible to use the existing read WFK to build a new one!")
 endif
 ! MRM change occs in Header Hdr
 ! Update occ in Hdr before printing
 ib1dm=1
 do ikpoint=1,BSt%nkpt
   dim_bands=size(BSt%occ(:,ikpoint,1))
   do ib2dm=1,dim_bands   
     Hdr%occ(ib1dm)=BSt%occ(ib2dm,ikpoint,1)
     ib1dm=ib1dm+1
   enddo
 enddo
 MSG_COMMENT("Hdr_sigma occupations correctly updated with nat. orb. ones")

 Hdr%npwarr(:)=Wfd%npwarr(:)                                   ! Use the npw = ones used in GW calc
 Hdr%ngfft(1:3)=ngfft_in(1:3)
 MSG_COMMENT("Hdr_sigma npw and ngfft correctly updated")

end subroutine update_hdr_bst
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

subroutine rotate_exchange(ikpoint,ib1,ib2,Sr,nateigv) ! Only used for debug of this file, do not use it with large arrays
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ikpoint 
 type(sigma_t) :: Sr
!arrays
 complex(dpc),intent(in) :: nateigv(:,:,:,:)
!Local variables ------------------------------
!scalars
 integer::ib1dm,ib2dm,ndim
 character(len=500) :: msg
!arrays
 complex(dpc),allocatable :: res(:,:),Umat(:,:),Kex_tmp(:,:)
!************************************************************************
 ndim=ib2-ib1+1
 ABI_MALLOC(res,(ndim,ndim))
 ABI_MALLOC(Umat,(ndim,ndim))
 ABI_MALLOC(Kex_tmp,(ndim,ndim))
 res=czero

 do ib1dm=1,ndim
   do ib2dm=1,ndim
     Umat(ib1dm,ib2dm)=nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),ikpoint,1)
     Kex_tmp(ib1dm,ib2dm)=Sr%x_mat(ib1+(ib1dm-1),ib1+(ib2dm-1),ikpoint,1)
   enddo
 enddo

 ! Print for debug
 !do ib1dm=ib1,ib2
 !  write(msg,'(a2,*(f10.5))') '  ',REAL(Sr%x_mat(ib1dm,ib1dm,ikpoint,1))
 !  call wrtout(std_out,msg,'COLL')
 !enddo

 write(msg,'(a3)') 'MAU'
 call wrtout(std_out,msg,'COLL')

 ! <KS|K[NO]KS> = U <NO|K[NO]|NO> (U^t)*
 res=matmul(Umat,Kex_tmp)
 Kex_tmp=matmul(res,conjg(transpose(Umat)))

 do ib1dm=1,ndim
   do ib2dm=1,ndim
     Sr%x_mat(ib1+(ib1dm-1),ib1+(ib2dm-1),ikpoint,1)=Kex_tmp(ib1dm,ib2dm)
   enddo
 enddo

 ! Print for debug
 do ib1dm=ib1,ib2
   write(msg,'(a2,*(f10.5))') '  ',REAL(Sr%x_mat(ib1dm,ib1dm,ikpoint,1))
   call wrtout(std_out,msg,'COLL')
 enddo

 ABI_FREE(res)
 ABI_FREE(Umat)
 ABI_FREE(Kex_tmp)
end subroutine rotate_exchange
!!***

end module m_gwrdm
!!***

!subroutine update_wfd_bst(wfd_i,wfd_f,nateigv,occs,b1gw,b2gw,BSt,Hdr,ngfft_in)
!!Arguments ------------------------------------
!!scalars
! integer,intent(in) :: b1gw,b2gw
! integer,intent(in),dimension(3) :: ngfft_in
! type(wfd_t),intent(inout) :: wfd_f
! type(wfd_t),intent(in) :: wfd_i
! type(ebands_t),target,intent(inout) :: BSt
! type(Hdr_type),intent(inout) :: Hdr
!!arrays
! real(dp),intent(in) :: occs(:,:)
! complex(dpc),intent(in) :: nateigv(:,:,:,:)
!Local variables ------------------------------
!scalars
! character(len=500) :: msg
! integer :: ib1dm,ib2dm,dim_bands,ikpoint
!!arrays
!!************************************************************************
! DBG_ENTER("COLL")
!
! !Wfd%Wave(6,2,1)%ug(3) ! BAND 6 , k-POINT 2, SPIN 1, UG="MO Coef"= 3
! do ikpoint=1,BSt%nkpt
!   write(msg,'(a31,i5,a9,i5,a1,i5,a7,i5)') ' Calc. nat. orbs coefs k-point: ',ikpoint,', bands: ',b1gw,'-',b2gw,', npw: ',wfd_i%Kdata(ikpoint)%npw ! Print for debug
!   call wrtout(std_out,msg,'COLL')
!   do ib1dm=b1gw,b2gw
!     wfd_f%Wave(ib1dm,ikpoint,1)%ug(:)=0.0d0 
!     do ib2dm=b1gw,b2gw
!       Wfd_f%Wave(ib1dm,ikpoint,1)%ug(:)=Wfd_f%Wave(ib1dm,ikpoint,1)%ug(:) &
!               +nateigv(ib2dm,ib1dm,ikpoint,1)*Wfd_i%Wave(ib2dm,ikpoint,1)%ug(:)
!     enddo
!     !write(*,'(*(f10.5))') real(nateigv(ib1dm,:,ikpoint,1)) !Print for debug
!   enddo
!   !write(*,*) ' ' !Space for debug
! enddo
!
! ! MRM BSt occ are changed and never recoverd
! MSG_COMMENT("QP_BSt occupations correctly updated with nat. orb. ones")
! do ikpoint=1,BSt%nkpt
!   BSt%occ(b1gw:b2gw,ikpoint,1) = occs(b1gw:b2gw,ikpoint) ! No spin used
! enddo
! if((size(Hdr%occ(:))/BSt%nkpt) < (b2gw-b1gw+1)) then
!   !Actually, we should never reach this point because the code should stop during Wfd initialization in m_sigma_driver.F90
!   MSG_ERROR("Impossible to use the existing read WFK to build a new one!")
! endif
! ! MRM change occs in Header Hdr
! ! Update occ in Hdr before printing
! MSG_COMMENT("Hdr_sigma occupations correctly updated with nat. orb. ones")
! ib1dm=1
! do ikpoint=1,BSt%nkpt
!   dim_bands=size(BSt%occ(:,ikpoint,1))
!   do ib2dm=1,dim_bands   
!     Hdr%occ(ib1dm)=BSt%occ(ib2dm,ikpoint,1)
!     ib1dm=ib1dm+1
!   enddo
! enddo
!
! MSG_COMMENT("Hdr_sigma npw and ngfft correctly updated")
! Hdr%npwarr(:)=Wfd_f%npwarr(:)                                   ! Use the npw = ones used in GW calc
! Hdr%ngfft(1:3)=ngfft_in(1:3)
!
!!end subroutine update_wfd_bst
!!***
