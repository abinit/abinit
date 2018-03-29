!{\src2tex{textfont=tt}}
!!****f* ABINIT/psichi_renormalization
!! NAME
!! psichi_renormalization
!!
!! FUNCTION
!! Renormalize psichi.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>= crystal structure data.
!!  paw_dmft =  data for LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!
!! OUTPUT
!!  paw_dmft%psichi(nsppol,nkpt,mband,nspinor,dtset%natom,(2*maxlpawu+1))): projections <Psi|chi> are orthonormalized.
!!
!! NOTES
!!
!! PARENTS
!!      datafordmft,dmft_solve
!!
!! CHILDREN
!!      add_matlu,destroy_oper,gather_matlu,identity_oper,init_oper
!!      invsqrt_matrix,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine psichi_renormalization(cryst_struc,paw_dmft,pawang,opt)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_pawang, only : pawang_type
 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_oper, only : init_oper,oper_type,identity_oper,loc_oper,destroy_oper,diff_oper
 use m_matlu, only : matlu_type,sym_matlu,print_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psichi_renormalization'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
! type(green_type),intent(inout) :: greenlda
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 integer, optional, intent(in) :: opt

!Local variables ------------------------------
!scalars
 integer :: iatom,ib,ikpt,im,ispinor,isppol,jkpt
 integer :: natom,mbandc,ndim,nkpt,nspinor,nsppol,option
 real(dp), pointer ::  temp_wtk(:) => null()
 real(dp) ::  pawprtvol
 type(oper_type) :: norm
 type(oper_type) :: oper_temp
 character(len=500) :: message
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************

 DBG_ENTER("COLL")

 option=3
 if(present(opt)) then
   if(opt==2) option=2
   if(opt==3) option=3
 end if
 pawprtvol=2

 nsppol  = paw_dmft%nsppol
 nkpt    = paw_dmft%nkpt
 mbandc  = paw_dmft%mbandc
 natom   = cryst_struc%natom
 nspinor = paw_dmft%nspinor


!== Normalize psichi
 if(option==1) then
!  ====================================
!  == simply renormalize psichi =======
!  ====================================
   write(message,'(2a)') ch10," Psichi are renormalized  "
   call wrtout(std_out,  message,'COLL')
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do ib=1,mbandc
         do iatom=1,natom
           if(paw_dmft%lpawu(iatom).ne.-1) then
             ndim=2*paw_dmft%lpawu(iatom)+1
             do im=1,ndim
               do ispinor=1,nspinor
!                write(std_out,*) "psichi1",paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)
                 paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)=     &
&                 paw_dmft%psichi(isppol,ikpt,ib,ispinor,iatom,im)/    &
&                 sqrt(real(norm%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)))
               end do ! ispinor
             end do ! im
           end if
         end do ! iatom
       end do ! ib
     end do ! ikpt
   end do ! isppol
!  todo_ab introduce correct orthonormalization in the general case.

 else if(option==2) then ! option==2
!  ====================================
!  == renormalize k-point after k-point
!  ====================================
   ABI_ALLOCATE(temp_wtk,(1))

   write(message,'(2a,i5)') ch10,' Nb of K-point',nkpt
   call wrtout(std_out,message,'COLL')
   do jkpt=1,nkpt  ! jkpt
     write(message,'(2a,i5)') ch10,'  == Renormalization for K-point',jkpt
     call wrtout(std_out,message,'COLL')
     temp_wtk(1)=one

     call normalizepsichi(cryst_struc,1,paw_dmft,pawang,temp_wtk,jkpt)
   end do ! jkpt
   ABI_DEALLOCATE(temp_wtk)

 else if(option==3)  then  ! option==3
!  ====================================
!  == renormalize the sum over k-points
!  ====================================
!  allocate(temp_wtk(nkpt)) 
   temp_wtk=>paw_dmft%wtk
   write(message,'(6a)') ch10,'  ====================================== '&
&   ,ch10,'  == Renormalization for all K-points == '&
&   ,ch10,'  ======================================='
   call wrtout(std_out,message,'COLL')
   call normalizepsichi(cryst_struc,nkpt,paw_dmft,pawang,temp_wtk)
!  deallocate(temp_wtk)

 end if 
 
!== Change back repr for norm 

!===============================================
!==  Compute norm with new psichi
!===============================================
 call init_oper(paw_dmft,norm,nkpt=paw_dmft%nkpt,wtk=paw_dmft%wtk)
!== Build identity for norm%ks (option=1)
 call identity_oper(norm,1)

!== Deduce norm%matlu from norm%ks with new psichi
 call loc_oper(norm,paw_dmft,1) 

!== Print norm%matlu unsymetrized with new psichi
 if(pawprtvol>2) then
   write(message,'(4a,2a)') &
&   ch10,"  == Check: Overlap with renormalized psichi without symetrization is == "
   call wrtout(std_out,message,'COLL')
   call print_matlu(norm%matlu,natom,prtopt=1)
 end if


!== Symetrise norm%matlu with new psichi
 call sym_matlu(cryst_struc,norm%matlu,pawang)

!== Print norm%matlu symetrized with new psichi
 if(pawprtvol>2) then
   write(message,'(4a,2a)') &
&   ch10,"  == Check: Overlap with renormalized psichi and symetrization is =="
   call wrtout(std_out,message,'COLL')
   call print_matlu(norm%matlu,natom,prtopt=1,opt_diag=-1)
 end if

!== Check that norm is now the identity
 call init_oper(paw_dmft,oper_temp)
 call identity_oper(oper_temp,2)
 call diff_oper('Overlap after renormalization','Identity',&
& norm,oper_temp,1,tol6)
 call destroy_oper(oper_temp)

 call destroy_oper(norm)

 paw_dmft%lpsichiortho=1

 DBG_EXIT("COLL")

 CONTAINS
!===========================================================
!!***

!!****f* psichi_renormalization/normalizepsichi
!! NAME
!!  normalizepsichi
!!
!! FUNCTION
!!  normalizepsichi psichi from norm
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  change psichi: normalizepsichi it
!!
!! PARENTS
!!      psichi_renormalization
!!
!! CHILDREN
!!      add_matlu,destroy_oper,gather_matlu,identity_oper,init_oper
!!      invsqrt_matrix,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

subroutine normalizepsichi(cryst_struc,nkpt,paw_dmft,pawang,temp_wtk,jkpt)

 use m_profiling_abi

 use defs_basis
 use m_errors

 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_oper, only : init_oper,oper_type,identity_oper,loc_oper,destroy_oper
 use m_matlu, only : gather_matlu,sym_matlu,print_matlu,add_matlu
 use m_matrix, only : invsqrt_matrix

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normalizepsichi'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nkpt
 integer,optional,intent(in) :: jkpt
 real(dp),pointer :: temp_wtk(:)
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 type(pawang_type), intent(in) :: pawang

!Local variables ------------------------------
!scalars
 integer :: diag,iatom,ib,ikpt1,im,im1,ispinor,ispinor1,isppol,isppol1,jc,jc1
 integer :: tndim
 integer :: natom,mbandc,ndim,nspinor,nsppol
 real(dp) :: pawprtvol
 type(oper_type) :: norm1,norm2,norm3
 character(len=500) :: message
 complex(dpc),allocatable :: wan(:,:,:),sqrtmatinv(:,:)
 type(coeff2c_type), allocatable :: overlap(:)
!arrays
! real(dp),allocatable :: e0pde(:,:,:),omegame0i(:)
!************************************************************************
   nsppol  = paw_dmft%nsppol
   mbandc  = paw_dmft%mbandc
   natom   = cryst_struc%natom
   nspinor = paw_dmft%nspinor
   pawprtvol=3
   diag=0

   if(nkpt/=1.and.present(jkpt)) then
     message = 'BUG in psichi_normalization'
     MSG_ERROR(message)
   end if

!  *********************************************************************
   call init_oper(paw_dmft,norm1,nkpt=nkpt,wtk=temp_wtk)
   
!  == Build identity for norm1%ks (option=1)
   call identity_oper(norm1,1)
   
   if(nkpt==1.and.present(jkpt)) then
     call loc_oper(norm1,paw_dmft,1,jkpt=jkpt) 
   end if
   if(.not.present(jkpt)) then
     call loc_oper(norm1,paw_dmft,1)
   end if
   if(nkpt>1) then
     call sym_matlu(cryst_struc,norm1%matlu,pawang)
   end if

   if(pawprtvol>2) then
     write(message,'(2a)') ch10,'  - Print norm with current psichi '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu,natom,prtopt=1,opt_exp=1)
   end if
!  ==-------------------------------------
!  == Start loop on atoms
   ABI_DATATYPE_ALLOCATE(overlap,(natom))
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom).ne.-1) then
       ndim=2*paw_dmft%lpawu(iatom)+1
       tndim=nsppol*nspinor*ndim
       ABI_ALLOCATE(overlap(iatom)%value,(tndim,tndim))
       overlap(iatom)%value=czero
     end if
   end do
!  ==-------------------------------------
   
!  built large overlap matrix
   write(message,'(2a)') ch10,'  - Overlap (before orthonormalization) -'
   call wrtout(std_out,message,'COLL')
   call gather_matlu(norm1%matlu,overlap,cryst_struc%natom,option=1,prtopt=1)
   call destroy_oper(norm1)


   
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom).ne.-1) then
       ndim=2*paw_dmft%lpawu(iatom)+1
       tndim=nsppol*nspinor*ndim
       ABI_ALLOCATE(sqrtmatinv,(tndim,tndim))
       
!      == Compute Inverse Square root of overlap : O^{-0.5}
       !!write(message,'(a,1x,a,e21.14,a,e21.14,a)') "overlap", &
       !!"(",real(overlap(1)%value),",",aimag(overlap(1)%value),")"
       !!call wrtout(std_out,message,'COLL')
       if(diag==0) then
         call invsqrt_matrix(overlap(iatom)%value,tndim)
         sqrtmatinv=overlap(iatom)%value
       else
         sqrtmatinv(:,:)=czero
         do ib=1,tndim
           sqrtmatinv(ib,ib)=cone/(sqrt(overlap(iatom)%value(ib,ib)))
         end do
       end if
       
!      == Apply O^{-0.5} on psichi 
       ABI_ALLOCATE(wan,(nsppol,nspinor,ndim))
!      write(std_out,*) mbandc,nsppol,nspinor,ndim
!      write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)
       do ikpt=1,nkpt
         do ib=1,mbandc
           if(present(jkpt)) then 
             ikpt1=jkpt
           else
             ikpt1=ikpt
           end if
           jc=0
           wan=czero
           do isppol=1,nsppol
             do ispinor=1,nspinor
               do im=1,ndim
!                write(std_out,*) "psichi", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
                 jc=jc+1
                 jc1=0
                 do isppol1=1,nsppol
                   do ispinor1=1,nspinor
                     do im1=1,ndim
                       jc1=jc1+1
                       wan(isppol,ispinor,im)= wan(isppol,ispinor,im) &
&                       + paw_dmft%psichi(isppol1,ikpt1,ib,ispinor1,iatom,im1)*sqrtmatinv(jc,jc1)
                     end do ! ispinor1
                   end do ! isppol1
                 end do ! im1
               end do ! im
             end do ! ispinor
           end do !  isppol
           do isppol=1,nsppol
             do ispinor=1,nspinor
               do im=1,ndim
                 paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)=wan(isppol,ispinor,im)
!                write(std_out,*) "psichi2", paw_dmft%psichi(isppol,ikpt1,ib,ispinor,iatom,im)
               end do ! ispinor
             end do ! isppol
           end do ! im
         end do ! ib
       end do ! ikpt
       ABI_DEALLOCATE(wan)
       ABI_DEALLOCATE(sqrtmatinv)
!      write(std_out,*)  paw_dmft%psichi(1,1,1,1,1,1)
       
!      ==-------------------------------------
     end if ! lpawu.ne.-1
   end do ! iatom
!  == End loop on atoms
!  ==-------------------------------------
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom).ne.-1) then
       ABI_DEALLOCATE(overlap(iatom)%value)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(overlap)

!  ======================================================================
!  == Check norm with new psichi.
!  ======================================================================

   call init_oper(paw_dmft,norm1,nkpt=nkpt,wtk=temp_wtk)

   call identity_oper(norm1,1)

   if(nkpt==1.and.present(jkpt)) then
     call loc_oper(norm1,paw_dmft,1,jkpt=jkpt) 
   end if
   if(.not.present(jkpt)) then
     call loc_oper(norm1,paw_dmft,1) 
   end if

   if (nkpt>1) then
     call sym_matlu(cryst_struc,norm1%matlu,pawang)
   end if

   if(pawprtvol>2) then
     write(message,'(2a)') ch10,'  - Print norm with new psichi '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm1%matlu,natom,prtopt=1)
   end if

!  ======================================================================
!  == Check that norm-identity is zero
!  ======================================================================
   call init_oper(paw_dmft,norm2,nkpt=nkpt,wtk=temp_wtk)
   call init_oper(paw_dmft,norm3,nkpt=nkpt,wtk=temp_wtk)
   call identity_oper(norm2,2)
   call add_matlu(norm1%matlu,norm2%matlu,norm3%matlu,natom,-1)
   call destroy_oper(norm2)
   if(pawprtvol>2) then
     write(message,'(2a)') ch10,'  - Print norm with new psichi minus Identity '
     call wrtout(std_out,message,'COLL')
     call print_matlu(norm3%matlu,natom,prtopt=1,opt_exp=1)
   end if
   call destroy_oper(norm3)

   call destroy_oper(norm1)
!  call flush(std_out)           ! debug debug  debug   debug 
!  MSG_ERROR("Stop for debugging")

 end subroutine normalizepsichi

end subroutine psichi_renormalization
!!***
