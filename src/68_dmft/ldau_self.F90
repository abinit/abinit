!{\src2tex{textfont=tt}}
!!****f* ABINIT/ldau_self
!! NAME
!! ldau_self
!!
!! FUNCTION
!! Use LDA+U to compute self-energy
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (BAmadon)
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
!!      paw_ij_free,paw_ij_init,paw_ij_nullify,pawpupot,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine ldau_self(cryst_struc,green,paw_dmft,pawtab,self,opt_ldau,prtopt)

 use m_profiling_abi

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_crystal, only : crystal_t
 use m_green, only : green_type
 use m_self, only : self_type
 use m_oper, only : print_oper
 use m_energy, only : energy_type,compute_energy

 use m_pawtab, only : pawtab_type
 use m_paw_ij, only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawdij, only : pawpupot
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ldau_self'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(inout) :: self !vz_i
 integer, intent(in) :: opt_ldau,prtopt

!Local variables ------------------------------
!scalars
 character(len=500) :: message
 integer :: iatom,idijeff,isppol,ispinor,ispinor1,itypat,lpawu
 integer :: ifreq,im,im1,ldim,natom,nspden,nsppol,nspinor,nsploop
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp),allocatable :: vpawu(:,:,:,:)
 type(paw_ij_type), allocatable :: paw_ij(:)

!************************************************************************

 natom=cryst_struc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 nspden=paw_dmft%nspden
 nsploop=max(paw_dmft%nsppol,paw_dmft%nspinor**2)
 isppol=0
 ispinor=0
 ispinor1=0
!write(std_out,*) "nsploop ",nsploop
 if(opt_ldau==1) then
 end if
 ABI_DATATYPE_ALLOCATE(paw_ij,(cryst_struc%natom))
 call paw_ij_nullify(paw_ij)
 call paw_ij_init(paw_ij,2,paw_dmft%nspinor,paw_dmft%nsppol,paw_dmft%nspden, &
& 1,paw_dmft%natom,cryst_struc%ntypat,cryst_struc%typat,pawtab,has_pawu_occ=1)

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
   itypat=cryst_struc%typat(iatom)
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu.ne.-1) then
     ldim=2*lpawu+1
     ABI_ALLOCATE(vpawu,(paw_ij(iatom)%cplex_dij,ldim,ldim,nspden))

     paw_ij(iatom)%nocctot(:)=zero ! contains nmmp in the n m representation
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
!          paw_ij(iatom)%noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!          paw_ij(iatom)%noccmmp(2,im,im1,idijeff)=imag(green%occup%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
           paw_ij(iatom)%noccmmp(1,im,im1,idijeff)=real(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
           paw_ij(iatom)%noccmmp(2,im,im1,idijeff)=aimag(green%occup%matlu(iatom)%mat(im1,im,isppol,ispinor,ispinor1))
         end do
       end do

!      paw_ij(1)%noccmmp(1,3,3,1)=0.d0
!      paw_ij(1)%noccmmp(1,5,5,1)=0.d0
!      paw_ij(2)%noccmmp(1,3,3,2)=0.d0
!      paw_ij(2)%noccmmp(1,5,5,2)=0.d0
!      ===============================
!      Compute nocctot for pawpupot =
!      ===============================
       do im1=1,ldim
         if(nsploop==4) then
           if(idijeff<=2) then
             paw_ij(iatom)%nocctot(1)=paw_ij(iatom)%nocctot(1)+paw_ij(iatom)%noccmmp(1,im1,im1,idijeff)
           end if
         else
           paw_ij(iatom)%nocctot(idijeff)=paw_ij(iatom)%nocctot(idijeff)+paw_ij(iatom)%noccmmp(1,im1,im1,idijeff)
         end if
       end do
!      write(message,'(2a)') ch10," == The noccmmp matrix is"
!      call wrtout(std_out,message,'COLL')
!      do im=1,ldim
!      write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(1,im,im1,idijeff),im1=1,ldim)
!      call wrtout(std_out,message,'COLL')
!      write(message,'(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(2,im,im1,idijeff),im1=1,ldim)
!      call wrtout(std_out,message,'COLL')
!      end do
!      !     write(message,'(2a)') ch10," == The nocctot are"
!      call wrtout(std_out,message,'COLL')
!      write(std_out,*) paw_ij(iatom)%nocctot(idijeff)
     end do
     paw_ij(iatom)%has_pawu_occ=2

!    warning  dmft works here if nspden=nsppol (this is checked elsewhere)

!    ===============================
!    Compute LDA+U vpawu from noccmmp
!    ===============================
     call pawpupot(paw_ij(iatom)%cplex_dij,paw_ij(iatom)%ndij,&
&     paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot,&
&     nspden,2,pawtab(itypat),vpawu)
!    do idijeff=1,nspden
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

   end if
 end do

 call paw_ij_free(paw_ij)
 ABI_DATATYPE_DEALLOCATE(paw_ij)

 if(abs(prtopt)>0) then
 end if

end subroutine ldau_self
!!***
