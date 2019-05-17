!!****m* ABINIT/m_ldau_self
!! NAME
!!  m_ldau_self
!!
!! FUNCTION
!! Compute DFT+U self energy for DMFT
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_ldau_self

 use defs_basis

 implicit none

 private

 public :: ldau_self
!!***

contains
!{\src2tex{textfont=tt}}
!!****f* m_ldau_self/ldau_self
!! NAME
!! ldau_self
!!
!! FUNCTION
!! Use LDA+U to compute self-energy
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for LDA.
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      impurity_solve,spectral_function
!!
!! CHILDREN
!!      pawpupot,wrtout
!!
!! SOURCE

subroutine ldau_self(cryst_struc,green,paw_dmft,pawtab,self,opt_ldau,prtopt)

 use m_abicore

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_self, only : self_type
 use m_oper, only : print_oper
 use m_energy, only : energy_type,compute_energy

 use m_pawtab, only : pawtab_type
 use m_pawdij, only : pawpupot
 use m_paw_dmft, only : paw_dmft_type
 use m_paw_correlations, only : setnoccmmp

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(pawtab_type),target,intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(inout) :: self !vz_i
 integer, intent(in) :: opt_ldau,prtopt

!Local variables ------------------------------
!scalars
 character(len=500) :: message
 integer :: iatom,idijeff,isppol,ispinor,ispinor1,lpawu
 integer :: ifreq,im,im1,ldim,natom,nocc,nsppol,nspinor,nsploop
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp),allocatable :: noccmmp(:,:,:,:),nocctot(:)
 real(dp),allocatable :: vpawu(:,:,:,:)
 type(pawtab_type),pointer :: pawtab_

!************************************************************************

 ABI_UNUSED(opt_ldau)

 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 nsploop=max(nsppol,nspinor**2)
 nocc=nsploop
 isppol=0
 ispinor=0
 ispinor1=0

!===============================
!Impose density matrix
!===============================
!todo_ab build pawrhoij or suppress it in setnoccmmp
!dimdmat=2*maxval(pawtab(:)%lpawu)+1
!call setnoccmmp(0,dimdmat,paw_dmft%dmatpawu(1:dimdmat,1:dimdmat,1:nsppol*nspinor,1:natpawu),&
!&   1,1,cryst_struc%indsym,cryst_struc%natom,cryst_struc%natom,natpawu,&
!&   paw_dmft%nspinor,paw_dmft%nsppol,cryst_struc%nsym,cryst_struc%ntypat,paw_ij,pawang,3,&
!&   pawrhoij_dum,pawtab,cryst_struc%spinat,cryst_struc%symafm,struct%typat,0,10)

 do iatom=1,cryst_struc%natom
   lpawu=paw_dmft%lpawu(iatom)
   pawtab_ => pawtab(cryst_struc%typat(iatom))
   if(lpawu.ne.-1) then
     ldim=2*lpawu+1
     ABI_ALLOCATE(vpawu,(2,ldim,ldim,nocc))

     ABI_ALLOCATE(noccmmp,(2,2*pawtab_%lpawu+1,2*pawtab_%lpawu+1,nocc))
     ABI_ALLOCATE(nocctot,(nocc))
     noccmmp(:,:,:,:)=zero ; nocctot(:)=zero ! contains nmmp in the n m representation

!    ===============================
!    Begin loop over spin/spinors to initialize noccmmp
     do idijeff=1,nsploop
!      ===============================
       if(nsploop==2) then
         isppol=spinor_idxs(1,idijeff)
         ispinor=1
         ispinor1=1
       else if(nsploop==4) then
         isppol=1
         ispinor=spinor_idxs(1,idijeff)
         ispinor1=spinor_idxs(2,idijeff)
       else if(nsploop==1) then
         isppol=1
         ispinor=1
         ispinor1=1
       else
         write(message,'(2a)') " BUG in ldau_self: nsploop should be equal to 1, 2 or 4"
         call wrtout(std_out,message,'COLL')
       end if
!      ===============================
!      Initialize noccmmp
!      ===============================
       do im1 = 1 , ldim
         do im = 1 ,  ldim
!          noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!          noccmmp(2,im,im1,idijeff)=imag(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
           noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
           noccmmp(2,im,im1,idijeff)=aimag(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
         end do
       end do

!      ===============================
!      Compute nocctot for pawpupot =
!      ===============================
       do im1=1,ldim
         if(nsploop==4) then
           if(idijeff<=2) then
             nocctot(1)=nocctot(1)+noccmmp(1,im1,im1,idijeff)
           end if
         else
           nocctot(idijeff)=nocctot(idijeff)+noccmmp(1,im1,im1,idijeff)
         end if
       end do
!      write(message,'(2a)') ch10," == The noccmmp matrix is"
!      call wrtout(std_out,message,'COLL')
!      do im=1,ldim
!      write(message,'(12(1x,9(1x,f10.5)))') (noccmmp(1,im,im1,idijeff),im1=1,ldim)
!      call wrtout(std_out,message,'COLL')
!      write(message,'(12(1x,9(1x,f10.5)))') (noccmmp(2,im,im1,idijeff),im1=1,ldim)
!      call wrtout(std_out,message,'COLL')
!      end do
!      !     write(message,'(2a)') ch10," == The nocctot are"
!      call wrtout(std_out,message,'COLL')
!      write(std_out,*) nocctot(idijeff)
     end do

!    warning  dmft works here if nspden=nsppol (this is checked elsewhere)

!    ===============================
!    Compute LDA+U vpawu from noccmmp
!    ===============================
     call pawpupot(2,nocc,noccmmp,nocctot,2,pawtab_,vpawu)
!    do idijeff=1,size(vpawu,4)
!    write(message,'(2a)') ch10," == The vpawu matrix is"
!    call wrtout(std_out,message,'COLL')
!    do im=1,ldim
!    write(message,'(12(1x,9(1x,f10.5)))') (vpawu(1,im,im1,idijeff),im1=1,ldim)
!    call wrtout(std_out,message,'COLL')
!    write(message,'(12(1x,9(1x,f10.5)))') (vpawu(2,im,im1,idijeff),im1=1,ldim)
!    call wrtout(std_out,message,'COLL')
!    end do
!    end do

!    ===============================
!    Begin loop over spin/spinors to compute self%oper
     do idijeff=1,nsploop
!      ===============================
       if(nsploop==2) then
         isppol=spinor_idxs(1,idijeff)
         ispinor=1
         ispinor1=1
       else if(nsploop==4) then
         isppol=1
         ispinor=spinor_idxs(1,idijeff)
         ispinor1=spinor_idxs(2,idijeff)
       else if(nsploop==1) then
         isppol=1
         ispinor=1
         ispinor1=1
       else
         write(message,'(2a)') " BUG in ldau_self: nsploop should be equal to 1, 2 or 4"
         call wrtout(std_out,message,'COLL')
       end if

!      ===============================
!      vpawu -> self%oper
!      ===============================
       do im1 = 1 , ldim
         do im = 1 ,  ldim
           do ifreq=1,self%nw
             self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
!            &             cmplx(vpawu(1,im1,im),vpawu(2,im1,im),kind=dp)
!            One take the transpose in orbital index to be coherent with the
!            current DFT+U implementation in Abinit.
&             cmplx(vpawu(1,im,im1,idijeff),vpawu(2,im,im1,idijeff),kind=dp)
           end do
         end do
       end do

     end do ! idijeff
!    write(std_out,*) "check im,im1 in vpawu",iatom
     ABI_DEALLOCATE(vpawu)

!    ===============================
!    Compute energy
!    ===============================

     ABI_DEALLOCATE(noccmmp)
     ABI_DEALLOCATE(nocctot)
   end if
 end do

 if(abs(prtopt)>0) then
 end if

end subroutine ldau_self
!!***

END MODULE m_ldau_self
!!***
