!!****m* ABINIT/m_matlu
!! NAME
!!  m_matlu
!!
!! FUNCTION
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
!! NOTES
!!  subroutines in this module must never call
!!   a subroutine of m_oper, m_green, m_self
!!   in order to avoid circular dependancies
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_matlu

 use defs_basis
 use m_errors
 use m_abicore

 use m_hide_lapack,  only : xginv
 use m_geometry, only : symredcart

 implicit none

 private

 public :: init_matlu
 public :: inverse_matlu
 public :: destroy_matlu
 public :: diff_matlu
 public :: add_matlu
 public :: print_matlu
 public :: sym_matlu
 public :: copy_matlu
 public :: gather_matlu
 public :: zero_matlu
 public :: trace_matlu
 public :: diag_matlu
 public :: rotate_matlu
 public :: shift_matlu
 public :: checkdiag_matlu
 public :: checkreal_matlu
 public :: prod_matlu
 public :: conjg_matlu
 public :: ln_matlu
 public :: slm2ylm_matlu
 public :: fac_matlu
 public :: printplot_matlu
 public :: identity_matlu
!!***


!!****t* m_matlu/matlu_type
!! NAME
!!  matlu_type
!!
!! FUNCTION
!!  This structured datatype contains an matrix for the correlated subspace
!!
!! SOURCE

 type, public :: matlu_type ! for each atom

  integer :: lpawu
  ! value of the angular momentum for correlated electrons

!  integer :: natom
   ! number of atoms (given for each atom, not useful..could be changed)
!
!  integer :: mband
!  ! Number of bands
!
!  integer :: mbandc
!  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT
!
!  integer :: natpawu         ! Number of correlated atoms
!
!  integer :: nkpt
!  ! Number of k-point in the IBZ.
  character(len=12) :: whichmatlu
  ! describe the type of local matrix computed (greenLDA, etc..)
!
  integer :: nspinor
  ! Number of spinorial component
!
  integer :: nsppol
  ! Number of polarization

  complex(dpc), allocatable :: mat(:,:,:,:,:)
  ! local quantity

 end type matlu_type

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_matlu/init_matlu
!! NAME
!! init_matlu
!!
!! FUNCTION
!!  Allocate variables used in type matlu_type.
!!
!! INPUTS
!!  natom = number of atoms
!!  nspinor = number of spinorial components
!!  nsppol = number of polarisation components
!!  lpawu_natom(natom) = value of lpawu for all atoms
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! OUTPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! PARENTS
!!      m_datafordmft,hubbard_one,m_green,m_matlu,m_oper,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_matlu(natom,nspinor,nsppol,lpawu_natom,matlu)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer :: natom,nspinor,nsppol
!type
 integer, intent(in) :: lpawu_natom(natom)
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables ------------------------------------
 integer :: iatom,lpawu

!************************************************************************

! matlu%mband       = mband
! matlu%dmftbandf   = dmftbandf
! matlu%dmftbandi   = dmftbandi
! matlu%nkpt        = nkpt
! matlu%mbandc  = 0
 matlu%nsppol      = nsppol
 matlu%nspinor     = nspinor
 do iatom=1,natom
  lpawu=lpawu_natom(iatom)
  matlu(iatom)%lpawu=lpawu
  if(lpawu.ne.-1) then
   ABI_ALLOCATE(matlu(iatom)%mat,(2*lpawu+1,2*lpawu+1,nsppol,nspinor,nspinor))
   matlu(iatom)%mat=czero
  else
   ABI_ALLOCATE(matlu(iatom)%mat,(0,0,nsppol,nspinor,nspinor))
  endif
 enddo

end subroutine init_matlu
!!***

!!****f* m_matlu/zero_matlu
!! NAME
!! zero_matlu
!!
!! FUNCTION
!!  zero_matlu
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!
!! PARENTS
!!      m_green,m_matlu
!!
!! CHILDREN
!!
!! SOURCE

subroutine zero_matlu(matlu,natom,onlynondiag)

 use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(matlu_type),intent(inout) :: matlu(natom)
 integer, optional, intent(in) :: onlynondiag
!Local variables-------------------------------
 integer :: iatom,im,im1,ispinor,ispinor1,isppol,ndim

!*********************************************************************

 if(present(onlynondiag)) then
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       do ispinor=1,matlu(iatom)%nspinor
         ndim=(2*matlu(iatom)%lpawu+1)
         do im=1,ndim
           do im1=1,ndim
             do ispinor1=1,matlu(iatom)%nspinor
               if(im/=im1.or.ispinor/=ispinor1) then
                 do isppol=1,matlu(iatom)%nsppol
                   matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=czero
                 enddo
               end if
             end do
           end do
         end do
       end do
     endif
   enddo
 else
   do iatom=1,natom
    matlu(iatom)%mat=czero
   enddo
 endif


end subroutine zero_matlu
!!***

!!****f* m_matlu/destroy_matlu
!! NAME
!! destroy_matlu
!!
!! FUNCTION
!!  deallocate matlu
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,hubbard_one,m_green,m_matlu,m_oper,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_matlu(matlu,natom)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 type(matlu_type),intent(inout) :: matlu(natom)

!Local variables-------------------------------
 integer :: iatom

! *********************************************************************

 do iatom=1,natom
  if ( allocated(matlu(iatom)%mat) )   then
    ABI_DEALLOCATE(matlu(iatom)%mat)
  end if
 enddo

end subroutine destroy_matlu
!!***

!!****f* m_matlu/copy_matlu
!! NAME
!! copy_matlu
!!
!! FUNCTION
!!  Copy matlu1 into matlu2
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!
!! PARENTS
!!      m_datafordmft,hubbard_one,m_green,m_matlu,m_oper,m_self,qmc_prep_ctqmc
!!      spectral_function
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_matlu(nmat1,nmat2,natom,opt_diag,opt_non_diag,opt_re)

 use defs_basis
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom
 type(matlu_type),intent(in) :: nmat1(natom)
 type(matlu_type),intent(inout) :: nmat2(natom) !vz_i
 integer, optional, intent(in) :: opt_diag,opt_non_diag,opt_re

!Local variables-------------------------------
 integer :: iatom,isppol,im1,im2,ispinor,ispinor1,tndim
! *********************************************************************

 do isppol=1,nmat1(1)%nsppol
   do iatom=1,natom
     if(nmat1(iatom)%lpawu.ne.-1) then
       tndim=(2*nmat1(iatom)%lpawu+1)
       do im1=1,tndim
         do im2=1,tndim
           do ispinor=1,nmat1(1)%nspinor
             do ispinor1=1,nmat1(1)%nspinor
               if(present(opt_diag)) then
                 nmat2(iatom)%mat(im1,im1,isppol,ispinor,ispinor)=nmat1(iatom)%mat(im1,im1,isppol,ispinor,ispinor)
               else if(present(opt_non_diag)) then
                 if(im1/=im2.or.ispinor/=ispinor1) then
                   nmat2(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)=nmat1(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)
                 endif
               else
                 nmat2(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)=nmat1(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)
                 if(present(opt_re)) nmat2(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)=&
&                         cmplx(real(nmat1(iatom)%mat(im1,im2,isppol,ispinor,ispinor1)),0.d0,kind=dp)
               endif
             enddo
           enddo
         enddo
       enddo
     endif ! lpawu=-1
   enddo ! iatom
 enddo ! isppol
!   do iatom=1,natom
!    lpawu=nmat1(iatom)%lpawu
!    if(lpawu.ne.-1) then
!     nmat2(iatom)%mat=nmat1(iatom)%mat
!    endif
!   enddo


end subroutine copy_matlu
!!***

!!****f* m_matlu/print_matlu
!! NAME
!! print_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  maltu <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  natom= number of atoms
!!  prtopt= option for printing
!!  opt_diag=   0   print non diagonal matrix (real or complex according to nspinor)
!!             -1   print non diagonal complex matrix
!!            >=1   print diagonal matrix (real or complex according to nspinor)
!!  opt_ab_out=  0  print matrix on std_out
!!             /=0  print matrix on ab_out
!!
!! OUTPUT
!!
!! PARENTS
!!      compute_levels,m_datafordmft,hubbard_one,impurity_solve,m_green,m_matlu
!!      m_oper,m_self,psichi_renormalization,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE

subroutine print_matlu(matlu,natom,prtopt,opt_diag,opt_ab_out,opt_exp,argout,compl)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: natom
 type(matlu_type),intent(in) :: matlu(natom)
 integer, intent(in) :: prtopt
 integer, optional, intent(in) :: opt_diag,opt_ab_out,opt_exp,argout,compl
!Local variables-------------------------------
 integer :: iatom,ispinor,ispinor1,isppol,m1,m,lpawu,nspinor,nsppol,optdiag,optab_out,arg_out
 character(len=500) :: message
 character(len=4) :: mode_paral
 type(matlu_type), allocatable :: matlu_tmp(:)
 character(len=9),parameter :: dspinm(2,2)= RESHAPE((/"n        ","mx       ","my       ","mz       "/),(/2,2/))
 logical :: testcmplx
 real(dp) :: noccspin
! *********************************************************************
 mode_paral='COLL'
 if(present(opt_diag)) then
   optdiag=opt_diag
 else
   optdiag=0
 endif
 if(present(opt_ab_out)) then
   optab_out=opt_ab_out
 else
   optab_out=0
 endif
 if(optab_out==0) then
   arg_out=std_out
 else
   arg_out=ab_out
 endif
 if(present(argout)) then
  arg_out=argout
  mode_paral='PERS'
 endif
 nspinor=matlu(1)%nspinor
 nsppol=matlu(1)%nsppol
 testcmplx=(nspinor==2)
 if(present(compl)) testcmplx=(nspinor==2).or.(compl==1)

 do iatom = 1 , natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu/=-1) then
     write(message,'(2a,i4)')  ch10,'   -------> For Correlated Atom', iatom
     call wrtout(arg_out,  message,mode_paral)
     do isppol = 1 , nsppol
       if(present(opt_ab_out).and.nsppol==2) then
         noccspin=zero
         do m1=1, 2*lpawu +1
           noccspin=noccspin+REAL(matlu(iatom)%mat(m1,m1,isppol,1,1))
         enddo
         !write(message,fmt='(7x,a,i3,a,f10.5)') ". Occ. for lpawu and for spin",isppol," =",noccspin
         !call wrtout(arg_out, message,mode_paral)
       endif
     enddo

     do isppol = 1 , nsppol
       if(nspinor == 1) then
         write(message,'(a,10x,a,i3,i3)')  ch10,'-- polarization spin component',isppol
         call wrtout(arg_out,  message,mode_paral)
       endif
       do ispinor = 1 , nspinor
         do ispinor1 = 1, nspinor
           if(nspinor == 2) then
             write(message,'(a,10x,a,i3,i3)')  ch10,'-- spin components',ispinor,ispinor1
             call wrtout(arg_out,  message,mode_paral)
           endif
           if(optdiag<=0) then
             do m1=1, 2*lpawu+1
               if(optdiag==0) then
                 if(nspinor==1.and.abs(prtopt)>0)  then
                   if(present(opt_exp)) then
                     write(message,'(5x,20e24.14)') (REAL(matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
!                     call wrtout(arg_out,  message,mode_paral)
!                     write(message,'(5x,20e20.14)') (REAL(sqrt(matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1))),m=1,2*lpawu+1)
!                     call wrtout(arg_out,  message,mode_paral)
!                     write(message,'(5x,20e20.14)') (REAL(1.d0/sqrt(matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1))),m=1,2*lpawu+1)
                   else
                     write(message,'(5x,20f10.5)') (REAL(matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                   endif
                 endif
                 if(testcmplx.and.abs(prtopt)>0) then
                   if(present(opt_exp)) then
                     if(opt_exp==2) then
                       write(message,'(5x,14(2e18.10,1x))') ((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                     else
                       write(message,'(5x,14(2e14.4,2x))') ((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                     endif
                   else
                     write(message,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                   endif
!&                  write(message,'(5x,14(2f15.11,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
                 endif
               else if(optdiag==-1) then
                 write(message,'(5x,14(2f10.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
               endif
               call wrtout(arg_out,  message,mode_paral)
             enddo
           elseif (optdiag>=1) then
             if(nspinor==1.and.abs(prtopt)>0) &
&             write(message,'(5x,20f10.5)') (REAL(matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
             if(testcmplx.and.abs(prtopt)>0) &
&             write(message,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
!            write(std_out,'(5x,14(2f9.5,2x))')((matlu(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
             call wrtout(arg_out,  message,mode_paral)
           endif
         end do ! ispinor1
       end do ! ispinor
       if(nspinor==2.and.prtopt>=5) then
         ABI_DATATYPE_ALLOCATE(matlu_tmp,(0:natom))
         call init_matlu(natom,nspinor,nsppol,matlu(1:natom)%lpawu,matlu_tmp)
         matlu_tmp(iatom)%mat(m1,m,isppol,1,1)= matlu(iatom)%mat(m1,m,isppol,1,1)+matlu(iatom)%mat(m1,m,isppol,2,2)
         matlu_tmp(iatom)%mat(m1,m,isppol,2,2)= matlu(iatom)%mat(m1,m,isppol,1,1)+matlu(iatom)%mat(m1,m,isppol,2,2)
         matlu_tmp(iatom)%mat(m1,m,isppol,1,2)= matlu(iatom)%mat(m1,m,isppol,1,2)+matlu(iatom)%mat(m1,m,isppol,2,1)
         matlu_tmp(iatom)%mat(m1,m,isppol,2,1)= &
&           (matlu(iatom)%mat(m1,m,isppol,1,2)+matlu(iatom)%mat(m1,m,isppol,2,1))*cmplx(0_dp,1_dp,kind=dp)
         do ispinor = 1 , nspinor
           do ispinor1 = 1, nspinor
             write(message,'(a,10x,a,a)')  ch10,'-- spin components',dspinm(ispinor,ispinor1)
             call wrtout(arg_out,  message,mode_paral)
             write(message,'(5x,14(2f9.5,2x))')((matlu_tmp(iatom)%mat(m1,m,isppol,ispinor,ispinor1)),m=1,2*lpawu+1)
             call wrtout(arg_out,  message,mode_paral)
           end do ! ispinor1
         end do ! ispinor
         call destroy_matlu(matlu_tmp,natom)
         ABI_DATATYPE_DEALLOCATE(matlu_tmp)
       endif
     enddo ! isppol
!     if(nsppol==1.and.nspinor==1) then
!       write(message,'(a,10x,a,i3,a)')  ch10,'-- polarization spin component',isppol+1,' is identical'
!       call wrtout(arg_out,  message,mode_paral)
!     endif
   endif ! lpawu/=1
 enddo ! natom


end subroutine print_matlu
!!***

!!****f* m_matlu/sym_matlu
!! NAME
!! sym_matlu
!!
!! FUNCTION
!! Symetrise local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  gloc(natom) <type(matlu_type)>= density matrix in the local orbital basis and related variables
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>= density matrix symetrized in the local orbital basis and related variables
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      compute_levels,m_datafordmft,hybridization_asymptotic_coefficient,m_green
!!      psichi_renormalization,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine sym_matlu(cryst_struc,gloc,pawang,paw_dmft)

 use defs_basis
! use defs_wvltypes
 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_paw_dmft, only: paw_dmft_type

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(pawang_type),intent(in) :: pawang
 type(paw_dmft_type), intent(in) :: paw_dmft
!arrays
 type(matlu_type),intent(inout) :: gloc(cryst_struc%natom)
!scalars
!Local variables-------------------------------
!scalars
 integer :: at_indx,iatom,irot,ispinor,ispinor1,isppol,lpawu,m1,m2,m3,m4,mu
 integer :: natom,ndim,nsppol,nspinor,nu,t2g,m1s,m2s,m3s,m4s,lpawu_zarot,x2my2d
 complex(dpc) :: sumrho,summag(3),rotmag(3),ci
 real(dp) :: zarot2
!arrays
! complex(dpc),allocatable :: glocnm(:,:,:,:,:)
 type(matlu_type),allocatable :: glocnm(:)
! complex(dpc),allocatable :: glocnms(:,:,:,:,:)
 type(matlu_type),allocatable :: glocnms(:)
 type(matlu_type),allocatable :: glocsym(:)
 real(dp),allocatable :: symrec_cart(:,:,:)
 integer :: mt2g(3),mx2my2d
 mt2g(1)=1
 mt2g(2)=2
 mt2g(3)=4
 mx2my2d=5
 t2g=paw_dmft%dmftqmc_t2g
 x2my2d=paw_dmft%dmftqmc_x2my2d

! DBG_ENTER("COLL")

 ci=cone
 nspinor=gloc(1)%nspinor
 nsppol=gloc(1)%nsppol
 natom=cryst_struc%natom

 ABI_DATATYPE_ALLOCATE(glocnm,(natom))
 ABI_DATATYPE_ALLOCATE(glocnms,(natom))
 ABI_DATATYPE_ALLOCATE(glocsym,(natom))
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocnm)
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocnms)
 call init_matlu(natom,nspinor,nsppol,gloc(1:natom)%lpawu,glocsym)


!=========  Case nspinor ==1 ========================

 if (nspinor==1) then
  ispinor=1
  ispinor1=1
  do iatom=1,cryst_struc%natom
   do isppol=1,nsppol
    if(gloc(iatom)%lpawu/=-1) then
     lpawu=gloc(iatom)%lpawu
     do m1=1, 2*lpawu+1
      do m2=1, 2*lpawu+1
       do irot=1,cryst_struc%nsym
        at_indx=cryst_struc%indsym(4,irot,iatom)
        do m3=1, 2*lpawu+1
         do m4=1, 2*lpawu+1
          if(t2g==1) then
           m1s=mt2g(m1)
           m2s=mt2g(m2)
           m3s=mt2g(m3)
           m4s=mt2g(m4)
           lpawu_zarot=2
          else if (x2my2d==1) then
           m1s=mx2my2d
           m2s=mx2my2d
           m3s=mx2my2d
           m4s=mx2my2d
           lpawu_zarot=2
          else
           m1s=m1
           m2s=m2
           m3s=m3
           m4s=m4
           lpawu_zarot=lpawu
          endif
          zarot2=pawang%zarot(m3s,m1s,lpawu_zarot+1,irot)*pawang%zarot(m4s,m2s,lpawu_zarot+1,irot)
          glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
&          glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)&
&          +gloc(at_indx)%mat(m3,m4,isppol,ispinor,ispinor1)*zarot2
         end do  ! m3
        end do  ! m4
       end do  ! irot
       glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)=&
&       glocsym(iatom)%mat(m1,m2,isppol,ispinor,ispinor1)/real(cryst_struc%nsym,kind=dp)
      end do ! m2
     end do ! m1
    endif ! lpawu/=-1
   end do ! isppol
  end do ! iatom
!==  Put glocsym into gloc
  do iatom=1,cryst_struc%natom
    if(gloc(iatom)%lpawu/=-1) then
      gloc(iatom)%mat=glocsym(iatom)%mat
!      gloc(iatom)%mat(:,:,1,:,:)=(glocsym(iatom)%mat(:,:,1,:,:) &
!&      + glocsym(iatom)%mat(:,:,2,:,:))/two
!      gloc(iatom)%mat(:,:,2,:,:)= gloc(iatom)%mat(:,:,1,:,:)
!      write(std_out,*) "WARNING: SYM non mag"
!      write(ab_out,*) "WARNING: SYM non mag"
    endif
  end do ! iatom

!=========  Case nspinor ==2 ========================

 else if (nspinor==2) then

!== Allocate temporary arrays
  do iatom=1,cryst_struc%natom
   if(gloc(iatom)%lpawu/=-1) then
    ndim=2*gloc(iatom)%lpawu+1
    ABI_DEALLOCATE(glocnm(iatom)%mat)
    ABI_DEALLOCATE(glocnms(iatom)%mat)
    ABI_DEALLOCATE(glocsym(iatom)%mat)
    ABI_ALLOCATE(glocnm(iatom)%mat,(ndim,ndim,nsppol,4,1))
    ABI_ALLOCATE(glocnms(iatom)%mat,(ndim,ndim,nsppol,4,1))
    ABI_ALLOCATE(glocsym(iatom)%mat,(ndim,ndim,nsppol,2,2))
   endif
  enddo
  ABI_ALLOCATE(symrec_cart,(3,3,cryst_struc%nsym))

!==  Compute symrec_cart
  do irot=1,cryst_struc%nsym
   call symredcart(cryst_struc%gprimd,cryst_struc%rprimd,symrec_cart(:,:,irot),cryst_struc%symrec(:,:,irot))
  end do

!==  Compute density matrix in density and magnetization representation
  call chg_repr_matlu(gloc,glocnm,cryst_struc%natom,option=1,prtopt=1)

!==  Do the sum over symetrized density matrix (in n,m repr)
  isppol=1
  do iatom=1,cryst_struc%natom
   if(gloc(iatom)%lpawu/=-1) then
    lpawu=gloc(iatom)%lpawu
    ndim=2*gloc(iatom)%lpawu+1
    do m1=1, 2*lpawu+1
     do m2=1, 2*lpawu+1
      sumrho=czero
      rotmag=czero
      do irot=1,cryst_struc%nsym
       summag=czero
       at_indx=cryst_struc%indsym(4,irot,iatom)
       do m3=1, 2*lpawu+1
        do m4=1, 2*lpawu+1
          if(t2g==1) then
           m1s=mt2g(m1)
           m2s=mt2g(m2)
           m3s=mt2g(m3)
           m4s=mt2g(m4)
           lpawu_zarot=2
          else if (x2my2d==1) then
           m1s=mx2my2d
           m2s=mx2my2d
           m3s=mx2my2d
           m4s=mx2my2d
           lpawu_zarot=2
          else
           m1s=m1
           m2s=m2
           m3s=m3
           m4s=m4
           lpawu_zarot=lpawu
          endif
         zarot2=pawang%zarot(m3s,m2s,lpawu_zarot+1,irot)*pawang%zarot(m4s,m1s,lpawu_zarot+1,irot)
         sumrho=sumrho +  glocnm(at_indx)%mat(m4,m3,isppol,1,1)  * zarot2
         do mu=1,3
          summag(mu)=summag(mu) + glocnm(at_indx)%mat(m4,m3,isppol,mu+1,1) * zarot2
         enddo
        end do ! m3
       end do !m4

!       ==  special case of magnetization
       do nu=1,3
        do mu=1,3
         rotmag(mu)=rotmag(mu)+symrec_cart(mu,nu,irot)*summag(nu)
        end do
       end do
!      write(std_out,'(a,3i4,2x,3(2f10.5,2x))') "rotmag",irot,m1,m2,(rotmag(mu),mu=1,3)
      end do ! irot

!       ==  Normalizes sum
      sumrho=sumrho/real(cryst_struc%nsym,kind=dp)
!        sumrho=glocnm(isppol,1,iatom,m1,m2) ! test without sym
      glocnms(iatom)%mat(m1,m2,isppol,1,1)=sumrho
      do mu=1,3
       rotmag(mu)=rotmag(mu)/real(cryst_struc%nsym,kind=dp)
!          rotmag(mu)=glocnm(isppol,mu+1,iatom,m1,m2) ! test without sym
       glocnms(iatom)%mat(m1,m2,isppol,mu+1,1)=rotmag(mu)
      enddo
     end do  ! m2
    end do ! m1
   endif ! lpawu/=-1
  end do ! iatom

!==  Compute back density matrix in upup dndn updn dnup representation
  call chg_repr_matlu(glocsym,glocnms,cryst_struc%natom,option=-1,prtopt=1)

!==  Put glocsym into gloc
  do iatom=1,cryst_struc%natom
    if(gloc(iatom)%lpawu/=-1) then
      gloc(iatom)%mat=glocsym(iatom)%mat
!      gloc(iatom)%mat(:,:,1,:,:)=(glocsym(iatom)%mat(:,:,1,:,:) &
!&      + glocsym(iatom)%mat(:,:,2,:,:))/two
!      gloc(iatom)%mat(:,:,2,:,:)= gloc(iatom)%mat(:,:,1,:,:)
!      write(std_out,*) "WARNING: SYM non mag"
!      write(ab_out,*) "WARNING: SYM non mag"
    endif
  end do ! iatom

  ABI_DEALLOCATE(symrec_cart)
 endif

 call destroy_matlu(glocnm,cryst_struc%natom)
 call destroy_matlu(glocnms,cryst_struc%natom)
 call destroy_matlu(glocsym,cryst_struc%natom)
 ABI_DATATYPE_DEALLOCATE(glocnm)
 ABI_DATATYPE_DEALLOCATE(glocnms)
 ABI_DATATYPE_DEALLOCATE(glocsym)
!==============end of nspinor==2 case ===========


! DBG_EXIT("COLL")

 end subroutine sym_matlu
!!***

!!****f* m_matlu/inverse_matlu
!! NAME
!! inverse_matlu
!!
!! FUNCTION
!! Inverse local quantity.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to inverse
!!  natom=number of atoms in cell.
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: inverse of input matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_oper
!!
!! CHILDREN
!!
!! SOURCE
 subroutine inverse_matlu(matlu,natom,prtopt)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
 integer :: iatom,tndim
 integer :: nsppol,nspinor
!scalars
 type(coeff2c_type),allocatable :: gathermatlu(:)
 !************************************************************************


 nspinor=matlu(1)%nspinor
 nsppol=matlu(1)%nsppol
 if(prtopt>0) then
 endif
 ABI_DATATYPE_ALLOCATE(gathermatlu,(natom))
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
     ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
     gathermatlu(iatom)%value=czero
   endif
 enddo

 call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     tndim=nsppol*nspinor*(2*matlu(iatom)%lpawu+1)
     !call matcginv_dpc(gathermatlu(iatom)%value,tndim,tndim)
     call xginv(gathermatlu(iatom)%value,tndim)
   endif
 enddo
 call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)

 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ABI_DEALLOCATE(gathermatlu(iatom)%value)
   endif
 enddo
 ABI_DATATYPE_DEALLOCATE(gathermatlu)
 end subroutine inverse_matlu
!!***

!!****f* m_matlu/diff_matlu
!! NAME
!! diff_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  char1 = character describing matlu1
!!  char2 = character describing matlu2
!!  matlu1(natom) <type(matlu_type)>= density matrix 1 in the local orbital basis and related variables
!!  matlu2(natom) <type(matlu_type)>= density matrix 2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  option =1      if diff> toldiff , stop
!!          0      print diff and toldiff
!!          else   do not test and do not print
!!  toldiff = maximum value for the difference between matlu1 and matlu2
!!
!! OUTPUT
!!
!! PARENTS
!!      m_datafordmft,m_green,m_oper,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
subroutine diff_matlu(char1,char2,matlu1,matlu2,natom,option,toldiff,ierr,zero_or_one)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_t
 use m_io_tools,           only : flush_unit
 implicit none

!Arguments ------------------------------------
!type
 integer,intent(in) :: natom,option
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 character(len=*), intent(in) :: char1,char2
 real(dp),intent(in) :: toldiff
 integer,intent(out), optional :: ierr
 integer,intent(in), optional :: zero_or_one

!local variables-------------------------------
 integer :: iatom,idiff,ispinor,ispinor1,isppol,m1,m,lpawu,nspinor,nsppol
 real(dp) :: matludiff
 character(len=500) :: message
! *********************************************************************
 nsppol=matlu1(1)%nsppol
 nspinor=matlu1(1)%nspinor

 matludiff=zero
 idiff=0
 do iatom = 1 , natom
  lpawu=matlu1(iatom)%lpawu
  if(lpawu/=-1) then
   do isppol = 1 , nsppol
    do ispinor = 1 , nspinor
     do ispinor1 = 1, nspinor
      do m1 = 1 , 2*lpawu+1
       do m = 1 ,  2*lpawu+1
        idiff=idiff+1
        matludiff=matludiff+ &
&        sqrt( real(matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1)      &
&            -      matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1))**2  &
&            +aimag(matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1)      &
&            -      matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1))**2  )
!       write(std_out,*) m,m1,matlu1(iatom)%mat(m1,m,isppol,ispinor,ispinor1),matlu2(iatom)%mat(m1,m,isppol,ispinor,ispinor1),matludiff
       enddo
      enddo
     end do ! ispinor1
    end do ! ispinor
   enddo ! isppol
  endif ! lpawu/=1
 enddo ! natom
 if(.not.present(zero_or_one)) matludiff=matludiff/float(idiff)

 if(option==1.or.option==0) then
  if( matludiff < toldiff ) then
   write(message,'(5a,6x,3a,4x,e12.4,a,e12.4)') ch10,&
&   '   ** Differences between ',trim(char1),' and ',ch10,trim(char2),' are small enough:',&
&   ch10,matludiff,' is lower than',toldiff
   call wrtout(std_out,message,'COLL')
   if(present(ierr)) ierr=0
  else 
   write(message,'(5a,3x,3a,3x,e12.4,a,e12.4)') ch10,&
&   'Differences between ',trim(char1),' and ',ch10,trim(char2),' is too large:',&
&   ch10,matludiff,' is larger than',toldiff
   MSG_WARNING(message)
!   write(message,'(8a,4x,e12.4,a,e12.4)') ch10,"  Matrix for ",trim(char1)
   write(message,'(a,3x,a)') ch10,trim(char1)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu1,natom,prtopt=1,opt_diag=-1)
   write(message,'(a,3x,a)') ch10,trim(char2)
   call wrtout(std_out,message,'COLL')
   call print_matlu(matlu2,natom,prtopt=1,opt_diag=-1)
   if (present(zero_or_one).and.(mod(matludiff,1.d0)< toldiff)) then
     write(message,'(a,3x,a)') ch10," The norm is not identity for this k-point but&
    & is compatible with a high symmetry point"
     call wrtout(std_out,message,'COLL')
   else if(present(zero_or_one)) then
     write(message,'(a,3x,a)') ch10," The norm is not identity for this k-point but&
    & might be compatible with a high symmetry point: it should be checked"
     call wrtout(std_out,message,'COLL')
   else 
     if(option==1) then
       call flush_unit(std_out)
       MSG_ERROR("option==1, aborting now!")
     end if
   end if
   if(present(ierr)) ierr=-1
  endif
 endif

end subroutine diff_matlu
!!***

!!****f* m_matlu/add_matlu
!! NAME
!! add_matlu
!!
!! FUNCTION
!!
!! INPUTS
!!  maltu1 <type(matlu_type)>= density matrix matlu1 in the local orbital basis and related variables
!!  maltu2 <type(matlu_type)>= density matrix matlu2 in the local orbital basis and related variables
!!  natom = number of atoms
!!  sign_matlu2= 1 add matlu1 and matlu2
!!              -1 substract matlu2 to matlu1
!!
!! OUTPUT
!!  maltu3 <type(matlu_type)>= density matrix matlu3, sum/substract matlu1 and matlu2
!!
!! PARENTS
!!      dyson,hybridization_asymptotic_coefficient,m_green
!!      psichi_renormalization,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
subroutine add_matlu(matlu1,matlu2,matlu3,natom,sign_matlu2)

 use defs_basis
 use m_paw_dmft, only : paw_dmft_type
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer,intent(in) :: natom,sign_matlu2
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 type(matlu_type), intent(inout) :: matlu3(natom) !vz_i

!local variables-------------------------------
 integer :: iatom
! *********************************************************************

 do iatom = 1 , natom
   matlu3(iatom)%mat=matlu1(iatom)%mat+float(sign_matlu2)*matlu2(iatom)%mat
 enddo ! natom

end subroutine add_matlu
!!***

!!****f* m_matlu/chg_repr_matlu
!! NAME
!! chg_repr_matlu
!!
!! FUNCTION
!! Change representation of density matrix (useful for nspinor=2)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  glocspsp(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the spin spin representation
!!  glocnm(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the magnetization representation
!!  natom=number of atoms in cell.
!!  option= 1 glocspsp is input, glocnm is computed
!!  option= -1 glocspsp is computed, glocnm is input
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  glocspsp(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the spin spin representation
!!  glocnm(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: density matrix in the magnetization representation
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_matlu
!!
!! CHILDREN
!!
!! SOURCE
 subroutine chg_repr_matlu(glocspsp,glocnm,natom,option,prtopt)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,prtopt
!arrays
 type(matlu_type),intent(inout) :: glocspsp(natom)
 type(matlu_type),intent(inout) :: glocnm(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,isppol,lpawu,m1,m2,ndim,nsppol,mu
 complex(dpc) :: ci
 character(len=500) :: message

! DBG_ENTER("COLL")

 ci=j_dpc

!==  Compute density matrix in density magnetization representation
 if (option==1) then
  nsppol=glocspsp(1)%nsppol
  do isppol=1,nsppol
   do iatom=1,natom
    if(glocspsp(iatom)%lpawu/=-1) then
     ndim=2*glocspsp(iatom)%lpawu+1
     do m1=1,ndim
      do m2=1,ndim
       glocnm(iatom)%mat(m1,m2,isppol,1,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,1) &
&                                         +glocspsp(iatom)%mat(m1,m2,isppol,2,2)
       glocnm(iatom)%mat(m1,m2,isppol,4,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,1) &
&                                         -glocspsp(iatom)%mat(m1,m2,isppol,2,2)
       glocnm(iatom)%mat(m1,m2,isppol,2,1)=glocspsp(iatom)%mat(m1,m2,isppol,1,2) &
&                                         +glocspsp(iatom)%mat(m1,m2,isppol,2,1)
       glocnm(iatom)%mat(m1,m2,isppol,3,1)= &
&                               cmplx((aimag(glocspsp(iatom)%mat(m1,m2,isppol,2,1))   &
&                                     -aimag(glocspsp(iatom)%mat(m1,m2,isppol,1,2))), &
&                                    (-real(glocspsp(iatom)%mat(m1,m2,isppol,2,1))+  &
&                                      real(glocspsp(iatom)%mat(m1,m2,isppol,1,2))),kind=dp)
      enddo  ! m2
     enddo ! m1
     if(abs(prtopt)>=3) then
      write(message,'(a)') "        -- in n, m repr "
      call wrtout(std_out,  message,'COLL')
      do mu=1,4
       do m1=1,ndim
        write(message,'(8x,(14(2f9.5,2x)))')(glocnm(iatom)%mat(m1,m2,isppol,mu,1),m2=1,ndim)
        call wrtout(std_out,  message,'COLL')
       enddo ! m1
       write(message,'(a)') ch10
       call wrtout(std_out,  message,'COLL')
      enddo ! mu
     endif ! prtopt >3
    endif ! lpawu/=-1
   enddo
  enddo

!==  Compute back density matrix in upup dndn updn dnup representation
 else  if (option==-1) then
  isppol=1
  do iatom=1,natom
   if(glocnm(iatom)%lpawu/=-1) then
    lpawu=glocnm(iatom)%lpawu
    ndim=2*glocnm(iatom)%lpawu+1
    do m1=1, 2*lpawu+1
     do m2=1, 2*lpawu+1
      glocspsp(iatom)%mat(m1,m2,isppol,1,1)=half*(glocnm(iatom)%mat(m1,m2,isppol,1,1)+glocnm(iatom)%mat(m1,m2,isppol,4,1))
      glocspsp(iatom)%mat(m1,m2,isppol,2,2)=half*(glocnm(iatom)%mat(m1,m2,isppol,1,1)-glocnm(iatom)%mat(m1,m2,isppol,4,1))
      glocspsp(iatom)%mat(m1,m2,isppol,1,2)=half*(glocnm(iatom)%mat(m1,m2,isppol,2,1)-ci*glocnm(iatom)%mat(m1,m2,isppol,3,1))
      glocspsp(iatom)%mat(m1,m2,isppol,2,1)=half*(glocnm(iatom)%mat(m1,m2,isppol,2,1)+ci*glocnm(iatom)%mat(m1,m2,isppol,3,1))
     end do  ! m2
    end do ! m1
    if(abs(prtopt)>6) then
     write(message,'(a)') "        -- in spin spin repr "
     call wrtout(std_out,  message,'COLL')
     do mu=1,4
      do m1=1,ndim
       write(message,'(8x,14(2f9.5,2x))')(glocspsp(iatom)%mat(m1,m2,isppol,mu,1),m2=1,ndim)
       call wrtout(std_out,  message,'COLL')
      enddo
      write(message,'(a)') ch10
      call wrtout(std_out,  message,'COLL')
     enddo
    endif !prtopt>3
   endif ! lpawu/=-1
  end do ! iatom
 else
  message = "stop in chg_repr_matlu"
  MSG_ERROR(message)
 endif


! DBG_EXIT("COLL")

 end subroutine chg_repr_matlu
!!***

!!****f* m_matlu/trace_matlu
!! NAME
!! trace_matlu
!!
!! FUNCTION
!!  Compute the trace of the matlu matrix
!!
!! INPUTS
!!  maltu(natom) <type(matlu_type)>= density matrix in the
!!               local orbital basis and related variables
!!  natom = number of atoms
!!
!! OUTPUT
!!  trace_loc(natom,nsppol+1)= trace for each atoms and each polarization,
!!                             trace_loc(iatom,nsppol+1) is
!!                             the full trace over polarization also.
!!
!! PARENTS
!!      m_green,m_oper
!!
!! CHILDREN
!!
!! SOURCE
 subroutine trace_matlu(matlu,natom,trace_loc,itau)

 use defs_basis
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in) :: natom
 type(matlu_type), intent(in) :: matlu(natom)
 real(dp),intent(inout),target,optional :: trace_loc(natom,matlu(1)%nsppol+1)
 integer, intent(in),optional :: itau

!local variables-------------------------------
 integer :: iatom,isppol,ispinor,m,lpawu
 integer :: nsppol,nspinor
 real(dp), pointer :: traceloc(:,:)=>null()
 character(len=500) :: message
! *********************************************************************
 nsppol=matlu(1)%nsppol
 nspinor=matlu(1)%nspinor
 if(present(trace_loc)) then
   traceloc=>trace_loc
 else
   ABI_ALLOCATE(traceloc,(natom,matlu(1)%nsppol+1))
 endif

 traceloc=zero
 do iatom = 1 , natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu/=-1) then
      write(message,'(2a,i4)')  ch10,'   -------> For Correlated Atom', iatom
     if(.not.present(itau)) then
       call wrtout(std_out,  message,'COLL')
     end if
     if(present(itau)) then
       if (itau>0) then
         call wrtout(std_out,  message,'COLL')
       end if
     endif
     do isppol = 1 , nsppol
       do ispinor = 1 , nspinor
         do m = 1 ,  2*lpawu+1
           traceloc(iatom,isppol)=traceloc(iatom,isppol)+&
&           matlu(iatom)%mat(m,m,isppol,ispinor,ispinor)
         enddo
       enddo
      traceloc(iatom,nsppol+1)=traceloc(iatom,nsppol+1)+traceloc(iatom,isppol)
     enddo
     if(nsppol==1.and.nspinor==1)  traceloc(iatom,nsppol+1)=traceloc(iatom,nsppol+1)*two
     if(.not.present(itau)) then
       write(message,'(8x,a,f12.6)')   'Nb of Corr. elec. from G(w) is:'&
&       ,traceloc(iatom,nsppol+1)
       call wrtout(std_out,  message,'COLL')
     endif
     if(present(itau)) then
       if(itau==1) then
         write(message,'(8x,a,f12.6)')   'Nb of Corr. elec. from G(tau) is:'&
&         ,traceloc(iatom,nsppol+1)
         call wrtout(std_out,  message,'COLL')
       else if(itau==-1) then
         write(message,'(8x,a,f12.6)')   'Nb: Sum of the values of G0(tau=0-) is:'&
&       ,traceloc(iatom,nsppol+1)
         call wrtout(std_out,  message,'COLL')
       else if(itau==4) then
         write(message,'(8x,a,f12.6)')   'Trace of matlu matrix is:'&
&         ,traceloc(iatom,nsppol+1)
         call wrtout(std_out,  message,'COLL')
       endif
     endif
   endif
 enddo
 if(.not.present(trace_loc)) then
  ABI_DEALLOCATE(traceloc)
  traceloc => null()
 endif

 end subroutine trace_matlu
!!***

!!****f* m_matlu/gather_matlu
!! NAME
!! gather_matlu
!!
!! FUNCTION
!! Create new array from matlu
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!  natom=number of atoms in cell.
!!  option= 1 go from gloc to gathergloc
!!  option= -1 go from gathergloc to gloc
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  gloc(natom) <type(matlu_type)>        = density matrix in the spin spin representation
!!  gatherloc(natom) <type(coeff2c_type)> = density matrix where spin and angular momentum are gathered in the same index
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_matlu,psichi_renormalization
!!
!! CHILDREN
!!
!! SOURCE
 subroutine gather_matlu(gloc,gathergloc,natom,option,prtopt)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

! type  matlus_type
!  SEQUENCE
!  complex(dpc), pointer :: mat(:,:)
! end type matlus_type

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,option,prtopt
 type(coeff2c_type), intent(inout) :: gathergloc(natom)
 type(matlu_type),intent(inout) :: gloc(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,ispinor,ispinor1,isppol,isppol1
 integer :: jc1,jc2,ml1,ml2,ndim,nspinor,nsppol,tndim
 character(len=500) :: message

! DBG_ENTER("COLL")
 nsppol=gloc(1)%nsppol
 nspinor=gloc(1)%nspinor

 do iatom=1,natom
   if(gloc(iatom)%lpawu.ne.-1) then
!==-------------------------------------

     ndim=2*gloc(iatom)%lpawu+1
     tndim=nsppol*nspinor*ndim

!== Put norm into array "gathergloc"
     jc1=0
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ml1=1,ndim
           jc1=jc1+1
           jc2=0
           do isppol1=1,nsppol
             do ispinor1=1,nspinor
               do ml2=1,ndim
                 jc2=jc2+1
                 if(option==1) then
                   if(isppol==isppol1) then
                     gathergloc(iatom)%value(jc1,jc2)=gloc(iatom)%mat(ml1,ml2,isppol,ispinor,ispinor1)
                   endif
                 else if(option==-1) then
                   if(isppol==isppol1) then
                     gloc(iatom)%mat(ml1,ml2,isppol,ispinor,ispinor1)=gathergloc(iatom)%value(jc1,jc2)
                   endif
                 endif
               enddo
             enddo ! ispinor1
           enddo ! isppol1
         enddo
       enddo !ispinor
     enddo ! isppol
   endif
 enddo ! iatom
 if(option==1.and.prtopt==3) then
   do iatom=1,natom
     if(gloc(iatom)%lpawu.ne.-1) then
       tndim=nsppol*nspinor*(2*gloc(iatom)%lpawu+1)
       write(message,'(2a,i5)') ch10,' (gathermatlu:) For atom', iatom
       call wrtout(std_out,message,'COLL')
       do im1=1,tndim
         write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&         (gathergloc(iatom)%value(im1,im2),im2=1,tndim)
         call wrtout(std_out,message,'COLL')
       end do
     endif
   enddo ! iatom
 else if(option==-1.and.prtopt==3) then
   call print_matlu(gloc,natom,prtopt)
 endif



! DBG_EXIT("COLL")

 end subroutine gather_matlu
!!***

!!****f* m_matlu/diag_matlu
!! NAME
!! diag_matlu
!!
!! FUNCTION
!! Diagonalize matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to diagonalize
!!  natom=number of atoms
!!  prtopt= option to define level of printing
!!  nsppol_imp= if 1, one can diagonalize with the same matrix the Up
!!   and Dn matlu matrix. It is convenient because one can thus have the
!! same interaction matrix for up and dn spins.
!!
!! OUTPUT
!!  matlu_diag(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: diagonalized density matrix
!!  eigvectmatlu(natom) <type(coeff2c_type)> = Eigenvectors corresponding to the diagonalization
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      hubbard_one,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine diag_matlu(matlu,matlu_diag,natom,prtopt,eigvectmatlu,nsppol_imp,checkstop,optreal,test)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 use m_matrix,         only : blockdiago_fordsyev,blockdiago_forzheev
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 type(matlu_type),intent(inout) :: matlu_diag(natom) !vz_i
 type(coeff2c_type),optional,intent(inout) :: eigvectmatlu(natom,matlu(1)%nsppol) !vz_i
 integer,optional,intent(in) :: nsppol_imp
 logical,optional,intent(in) :: checkstop
 integer,optional,intent(in) :: optreal
 integer,optional,intent(in) :: test
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,im3,imc,imc1,info,ispinor,ispinor1,isppol,lwork,tndim,lworkr
 integer :: nsppol,nsppolimp,nspinor
 logical :: checkstop_in,blockdiag
 character(len=500) :: message
!arrays
 type(coeff2c_type),allocatable :: gathermatlu(:)
 real(dp),allocatable :: eig(:),rwork(:),work(:),valuer(:,:)!,valuer2(:,:)
 !real(dp),allocatable :: valuer3(:,:),valuer4(:,:)
! real(dp),allocatable :: eigvec(:,:)
 complex(dpc),allocatable :: zwork(:)
 logical :: donotdiag,print_temp_mat2
 complex(dpc),allocatable :: temp_mat(:,:)
 complex(dpc),allocatable :: temp_mat2(:,:)
!debug complex(dpc),allocatable :: temp_mat3(:,:)
!************************************************************************

 call zero_matlu(matlu_diag,natom)
 nsppol=matlu(1)%nsppol
 if(present(nsppol_imp)) then
   nsppolimp=nsppol_imp
 else
   nsppolimp=nsppol
 endif
 if(present(checkstop)) then
  checkstop_in=checkstop
 else
  checkstop_in=.true.
 endif
 if(present(test)) then
  blockdiag=(test==8)
 else
  blockdiag=.false.
 endif
 nspinor=matlu(1)%nspinor

 donotdiag=.true.
 donotdiag=.false.
! ===========================
! Check is diagonalization is necessary and how
! ===========================
 do isppol=1,matlu(1)%nsppol
   do iatom=1,natom
      if(matlu(iatom)%lpawu.ne.-1) then
       tndim=(2*matlu(iatom)%lpawu+1)
        do im1=1,tndim
          do im2=1,tndim
            do ispinor=1,nspinor
              do ispinor1=1,nspinor
                if(abs(matlu(iatom)%mat(im1,im2,isppol,ispinor,ispinor1))>tol8.and.&
&                 (im1/=im2.or.ispinor/=ispinor1)) then
!                 if matrix is diagonal: do not diagonalize
                  donotdiag=.false.
                  exit
                endif
              enddo
            enddo
          enddo
        enddo
      endif
   enddo
 enddo

 if(donotdiag) then
   do isppol=1,matlu(1)%nsppol
     do iatom=1,natom
        if(matlu(iatom)%lpawu.ne.-1) then
          tndim=(2*matlu(iatom)%lpawu+1)
          eigvectmatlu(iatom,isppol)%value(:,:)=czero
          do im1=1,tndim
            eigvectmatlu(iatom,isppol)%value(im1,im1)=cone
          enddo
        endif
     enddo
   enddo
   call copy_matlu(matlu,matlu_diag,natom)
   write(message,'(a)')  "   Diagonalisation of matlu will not be performed"
   call wrtout(std_out,message,'COLL')
   return
 endif


! Below, gathermatlu is used to store eigenvectors (after diagonalization)
! For nsppol=2, and if nsppolimp=1, these eigenvectors computed for isppol=1, and applied to
! rotate matrix for isppol=2. It is the reason why the sum below is only
! from 1 to nsppolimp !



 do isppol=1,nsppolimp !
! ===========================
! Define gathermatlu
! ===========================
   ABI_DATATYPE_ALLOCATE(gathermatlu,(natom))
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
       gathermatlu(iatom)%value=czero
     endif
   enddo
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
           enddo
         enddo
       endif
     enddo
   endif
! ===========================
! Diagonalize
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       if(isppol<=nsppolimp) then
         tndim=nspinor*(2*matlu(iatom)%lpawu+1)
!debug       allocate(temp_mat2(tndim,tndim))
!debug       temp_mat2=zero
         lwork=2*tndim-1
         ABI_ALLOCATE(rwork,(3*tndim-2))
         rwork = zero

         lworkr=tndim*(tndim+2)*2
         ABI_ALLOCATE(work,(lworkr))
         work = zero
         ABI_ALLOCATE(valuer,(tndim,tndim))
!         ABI_ALLOCATE(valuer2,(tndim,tndim))
!         ABI_ALLOCATE(valuer3,(tndim,tndim))
!         ABI_ALLOCATE(valuer4,(tndim,tndim))
         ABI_ALLOCATE(zwork,(lwork))
!         valuer2=zero
!         valuer3=zero
!         valuer4=zero
         zwork = czero
         ABI_ALLOCATE(eig,(tndim))
         eig = zero
         info = 0
         if(prtopt>=4) then
           write(message,'(a,i4,a,i4)')  "       BEFORE DIAGONALIZATION for atom",iatom,"  and isppol",isppol
           call wrtout(std_out,message,'COLL')
           do im1=1,tndim
             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&             (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
             call wrtout(std_out,message,'COLL')
           end do
         endif
!debug       temp_mat2(:,:)=gathermatlu(iatom)%value(:,:)
!           write(std_out,*)"diag"
         if(present(optreal).and.maxval(abs(aimag(gathermatlu(iatom)%value(:,:))))<tol6) then
           write(message,'(a,2x,a,e9.3,a)') ch10,"Imaginary part of Local Hamiltonian is lower than ",&
&                   tol8, ": the real matrix is used"
           call wrtout(std_out,message,'COLL')
           valuer=real(gathermatlu(iatom)%value,kind=dp)
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,"(",f20.15,",",f20.15,")")))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           do im1=1,tndim
!             valuer(im1,im1)=real(im1,kind=dp)*0.00000000001_dp+valuer(im1,im1)
!           enddo
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,f20.15,f20.15)))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
           !call dsyev('v','u',tndim,valuer,tndim,eig,work,lworkr,info)
           if (blockdiag) then
             call blockdiago_fordsyev(valuer,tndim,eig)
           else
             call dsyev('v','u',tndim,valuer,tndim,eig,work,lworkr,info)
           endif
!!       For reproductibility
!           ! valuer2: eigenvector for the perturb matrix
!           valuer2=real(gathermatlu(iatom)%value,kind=dp)
!           do im1=1,tndim
!             valuer2(im1,im1)=float(im1)*0.00000000001+valuer2(im1,im1)
!           enddo
!           call dsyev('v','u',tndim,valuer2,tndim,eig,work,lworkr,info)
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "       valuer2 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer2(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dgemm('n','n',tndim,tndim,tndim,cone,valuer,tndim,&
!&            valuer2,tndim,czero,valuer3                ,tndim)
!           call dgemm('c','n',tndim,tndim,tndim,cone,valuer2,tndim,&
!&            valuer3                   ,tndim,czero,valuer4,tndim)
!           ! valuer4: compute unpert matrix in the basis of the
!           ! perturb basis
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "BEFORE valuer4 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer4(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dsyev('v','u',tndim,valuer4,tndim,eig,work,lworkr,info)
!           ! valuer4: Diago valuer4 (nearly diag)
!           write(message,'(a)') ch10
!           call wrtout(std_out,message,'COLL')
!           write(message,'(a,i4,a,i4)')  "AFTER  valuer4 for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (valuer4(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
!           call dgemm('n','n',tndim,tndim,tndim,cone,valuer2,tndim,&
!&            valuer4,tndim,czero,valuer                ,tndim)
           !write(6,*) "INFO",info
           gathermatlu(iatom)%value=cmplx(valuer,0.d0,kind=dp)
!           write(message,'(a,i4,a,i4)')  "AFTER valuer for atom",iatom,"  and isppol",isppol
!           call wrtout(std_out,message,'COLL')
!           do im1=1,tndim
!             write(message,'(2(1x,18(1x,"(",f20.15,",",f20.15,")")))')&
!&             (valuer(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
         else
           if(present(optreal).and.maxval(abs(aimag(gathermatlu(iatom)%value(:,:))))>tol8) then
             write(message,'(a)') " Local hamiltonian in correlated basis is complex"
             MSG_COMMENT(message)
           endif
           call zheev('v','u',tndim,gathermatlu(iatom)%value,tndim,eig,zwork,lwork,rwork,info)
           !call blockdiago_forzheev(gathermatlu(iatom)%value,tndim,eig)
         endif
         if(prtopt>=3) then
           write(message,'(a)') ch10
           call wrtout(std_out,message,'COLL')
           write(message,'(a,i4,a,i4)')  "       EIGENVECTORS for atom",iatom,"  and isppol",isppol
           call wrtout(std_out,message,'COLL')
           do im1=1,tndim
             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&             (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
             call wrtout(std_out,message,'COLL')
           end do
          ! do im1=1,tndim
          !   xcheck=czero
          !   do im3=1,tndim
          !     do im2=1,tndim
          !       xcheck=xcheck+gathermatlu(iatom)%value(im1,im2)*conjg(gathermatlu(iatom)%value(im2,im3))
          !     end do
          !   end do
          !   write(6,*) "check",im3,im1,xcheck
          ! end do
         endif
!       write(std_out,*) "eig",eig
! ===========================
! Put eigenvalue in matlu_diag
! ===========================
         imc=0
         do ispinor=1,nspinor
           do im1=1,2*matlu(iatom)%lpawu+1
             imc=imc+1
             matlu_diag(iatom)%mat(im1,im1,isppol,ispinor,ispinor)=eig(imc)
!debug             temp_mat2(imc,imc)=eig(imc)
           enddo
         enddo
         if(prtopt>=2) then
!            write(message,'(a,12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&            ch10,(eig(im1),im1=1,tndim)
!             call wrtout(std_out,message,'COLL')
           !call wrtout(std_out,message,'COLL')
            !write(std_out,*) "EIG", eig
         endif
         ABI_DEALLOCATE(zwork)
         ABI_DEALLOCATE(rwork)
         ABI_DEALLOCATE(work)
         ABI_DEALLOCATE(valuer)
!         ABI_DEALLOCATE(valuer2)
!         ABI_DEALLOCATE(valuer3)
!         ABI_DEALLOCATE(valuer4)
         ABI_DEALLOCATE(eig)
!     endif
!   enddo
! ===========================
! Keep eigenvectors gathermatlu
! ===========================
         if (present(eigvectmatlu)) then
           tndim=nspinor*(2*matlu(iatom)%lpawu+1)
           eigvectmatlu(iatom,isppol)%value(:,:)=gathermatlu(iatom)%value(:,:)
!           write(std_out,*) "eigvect in diag_matlu"
!           do im1=1,tndim
!             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!&             (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
!             call wrtout(std_out,message,'COLL')
!           end do
         endif


! ==================================================================
       endif
       if(nsppolimp==1.and.matlu(1)%nsppol==2) then
! If necessary rotate levels for this other spin, assuming the same
! rotation matrix: it have to be checked afterwards that the matrix is
! diagonal
! ===================================================================
         ABI_ALLOCATE(temp_mat,(tndim,tndim))
         ABI_ALLOCATE(temp_mat2,(tndim,tndim))
         temp_mat(:,:)=czero
!        input matrix: gathermatlu
!        rotation matrix: eigvectmatlu
!        intermediate matrix: temp_mat
!        result matrix: temp_mat2
         do im1=1,tndim
           do im2=1,tndim
             gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,2,1,1)
           enddo
         enddo
         if(prtopt>=3) then
           write(message,'(a,i4,a,i4)')  "       GATHERMATLU for atom",iatom," inside if nsppolimp==1"
           call wrtout(std_out,message,'COLL')
           do im1=1,tndim
             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&             (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
             call wrtout(std_out,message,'COLL')
           end do
         endif
         call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,&
&          eigvectmatlu(iatom,1)%value,tndim,czero,temp_mat                ,tndim)
         call zgemm('c','n',tndim,tndim,tndim,cone,eigvectmatlu(iatom,1)%value,tndim,&
&          temp_mat                   ,tndim,czero,temp_mat2,tndim)
         eigvectmatlu(iatom,2)%value=eigvectmatlu(iatom,1)%value
         imc=0
         print_temp_mat2=.false.
         do ispinor=1,nspinor
           do im1=1,2*matlu(iatom)%lpawu+1
             imc=imc+1
             imc1=0
             do ispinor1=1,nspinor
               do im2=1,2*matlu(iatom)%lpawu+1
                 imc1=imc1+1
                 matlu_diag(iatom)%mat(im1,im2,2,ispinor,ispinor1)=temp_mat2(imc,imc1)
                 if (imc/=imc1.and.(abs(temp_mat2(imc,imc1))>tol5)) then
                   write(message,'(3a,i4,2f16.4)') ch10,'diag_matlu= Matrix for spin number 2 obtained with', &
&                   ' eigenvectors from diagonalization for spin nb 1 is non diagonal for atom:',iatom,&
&                    abs(temp_mat2(imc,imc1)),tol5
                   call wrtout(std_out,message,'COLL')
                   if(.not.checkstop_in.and.(abs(temp_mat2(imc,imc1))>0.1)) print_temp_mat2=.true.
                   if(checkstop_in) print_temp_mat2=.true.
                 endif
!                 write(std_out,*) temp_mat2(imc,imc1)
!debug             temp_mat2(imc,imc)=eig(imc)
               enddo
             enddo
           enddo
         enddo
         if(print_temp_mat2.and.prtopt>=3) then
           write(message,'(a)')  "       temp_mat2"
           call wrtout(std_out,message,'COLL')
           do im1=1,tndim
             write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
&             (temp_mat2(im1,im3),im3=1,tndim)
             call wrtout(std_out,message,'COLL')
           end do
           if(iatom==2) then
             MSG_ERROR("iatom==2")
           end if
         endif
         ABI_DEALLOCATE(temp_mat)
         ABI_DEALLOCATE(temp_mat2)
       endif


     endif ! end lpawu/=-1

!!  for check only
!debug     if(matlu(iatom)%lpawu.ne.-1) then
!debug       allocate(temp_mat(tndim,tndim))
!debug       allocate(temp_mat3(tndim,tndim))
!debug           do im1=1,tndim
!debug             do im2=1,tndim
!debug!               rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom,isppol)%value(im1,im2)
!debug               temp_mat3(im1,im2)=conjg(gathermatlu(iatom)%value(im2,im1))
!debug             enddo
!debug           enddo
!debug       temp_mat(:,:)=czero
!debug!      input matrix: temp_mat2
!debug!      rotation matrix: gathermatlu
!debug!      intermediate matrix: temp_mat
!debug!      result matrix: temp_mat2
!debug       call zgemm('n','c',tndim,tndim,tndim,cone,temp_mat2   ,tndim,&
!debug&        temp_mat3,tndim,czero,temp_mat                ,tndim)
!debug       call zgemm('n','n',tndim,tndim,tndim,cone,temp_mat3,tndim,&
!debug&        temp_mat                   ,tndim,czero,temp_mat2,tndim)
!debug!       call zgemm('n','c',tndim,tndim,tndim,cone,temp_mat2   ,tndim,&
!debug!&        gathermatlu(iatom)%value,tndim,czero,temp_mat                ,tndim)
!debug!       call zgemm('n','n',tndim,tndim,tndim,cone,gathermatlu(iatom)%value,tndim,&
!debug!&        temp_mat                   ,tndim,czero,temp_mat2,tndim)
!debug         write(std_out,*) "result"
!debug         do im1=1,tndim
!debug           write(message,'(12(1x,18(1x,"(",f9.3,",",f9.3,")")))')&
!debug&           (temp_mat2(im1,im2),im2=1,tndim)
!debug           call wrtout(std_out,message,'COLL')
!debug         end do
!debug       deallocate(temp_mat)
!debug       deallocate(temp_mat3)
!debug     endif ! lpawu

   enddo  ! iatom
! End loop over atoms
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
!debug       deallocate(temp_mat2)
       ABI_DEALLOCATE(gathermatlu(iatom)%value)
     endif
   enddo
   ABI_DATATYPE_DEALLOCATE(gathermatlu)
 enddo ! isppol

 end subroutine diag_matlu
!!***

!!****f* m_matlu/rotate_matlu
!! NAME
!! rotate_matlu
!!
!! FUNCTION
!! Rotate matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  rot_mat(natom) <type(coeff2c_type)> = Rotation matrix (usually from diag_matlu)
!!  natom=number of atoms in cell.
!!  prtopt= option to define level of printing
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: Rotated matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      hubbard_one,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine rotate_matlu(matlu,rot_mat,natom,prtopt,inverse)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
 integer, intent(in) :: prtopt
 integer, intent(in) :: inverse
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 type(coeff2c_type),optional,intent(inout) :: rot_mat(natom,matlu(1)%nsppol)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,isppol
 integer :: nsppol,nspinor,tndim
!arrays
 type(coeff2c_type),allocatable :: gathermatlu(:)
! type(coeff2c_type),allocatable :: rot_mat_orig(:,:)
 type(coeff2c_type),allocatable :: rot_mat_orig(:)
 complex(dpc),allocatable :: temp_mat(:,:)
!************************************************************************
 if(prtopt==1) then
 endif
 nsppol=matlu(1)%nsppol
 nspinor=matlu(1)%nspinor
! ABI_DATATYPE_ALLOCATE(rot_mat_orig,(natom,matlu(1)%nsppol))

 do isppol=1,nsppol

! ===========================
! Define gathermatlu and rot_mat_orig and allocate
! ===========================
   ABI_DATATYPE_ALLOCATE(rot_mat_orig,(natom))
   ABI_DATATYPE_ALLOCATE(gathermatlu,(natom))
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       ABI_ALLOCATE(gathermatlu(iatom)%value,(tndim,tndim))
       gathermatlu(iatom)%value=czero
!       ABI_ALLOCATE(rot_mat_orig(iatom,isppol)%value,(tndim,tndim))
!       rot_mat_orig(iatom,isppol)%value(:,:)=rot_mat(iatom,isppol)%value(:,:)
       ABI_ALLOCATE(rot_mat_orig(iatom)%value,(tndim,tndim))
       rot_mat_orig(iatom)%value(:,:)=rot_mat(iatom,isppol)%value(:,:)
     endif
   enddo
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=1,prtopt=1)
   else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             gathermatlu(iatom)%value(im1,im2)=matlu(iatom)%mat(im1,im2,isppol,1,1)
           enddo
         enddo
       endif
     enddo
   endif
        ! write(std_out,*) "gathermatlu in rotate matlu"
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e17.10,",",e17.10,")")))')&
        !    (gathermatlu(1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do

! ===========================
! If necessary, invert rot_mat
! ===========================
   if(inverse==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         tndim=nspinor*(2*matlu(iatom)%lpawu+1)
           do im1=1,tndim
             do im2=1,tndim
!               rot_mat(iatom,isppol)%value(im1,im2)=conjg(rot_mat_orig(iatom,isppol)%value(im2,im1))
               rot_mat(iatom,isppol)%value(im1,im2)=conjg(rot_mat_orig(iatom)%value(im2,im1))
             enddo
           enddo
       endif ! lpawu
     enddo ! iatom
   endif
        ! write(std_out,*) "rot_mat_orig "
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e18.10,",",e18.10,")")))')&
        ! &   (rot_mat_orig(1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do
        ! write(std_out,*) "rot_mat "
        ! do im1=1,tndim
        !   write(message,'(12(1x,18(1x,"(",e18.10,",",e18.10,")")))')&
        ! &   (rot_mat(1,1)%value(im1,im2),im2=1,tndim)
        !   call wrtout(std_out,message,'COLL')
        ! end do

! ===========================
! Rotate
! ===========================
   ABI_ALLOCATE(temp_mat,(tndim,tndim))
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       tndim=nspinor*(2*matlu(iatom)%lpawu+1)
       temp_mat(:,:)=czero
!      input matrix: gathermatlu
!      rotation matrix: rot_mat
!      intermediate matrix: temp_mat
!      result matrix: gathermatlu
       ! temp_mat = gathermatlu * conjg(rot_mat)
       call zgemm('n','c',tndim,tndim,tndim,cone,gathermatlu(iatom)%value   ,tndim,&
&        rot_mat(iatom,isppol)%value,tndim,czero,temp_mat                ,tndim)
       ! gathermatlu = rot_mat * temp_mat = rot_mat * gathermatlu * conjg(rot_mat)
       call zgemm('n','n',tndim,tndim,tndim,cone,rot_mat(iatom,isppol)%value,tndim,&
&        temp_mat                   ,tndim,czero,gathermatlu(iatom)%value,tndim)
     endif ! lpawu
   enddo ! iatom
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    write(std_out,*) "temp_mat in rotate_matlu 2"
   !    do im1=1,tndim
   !      write(message,'(12(1x,18(1x,"(",f17.10,",",f17.10,")")))')&
   !&       (temp_mat(im1,im2),im2=1,tndim)
   !      call wrtout(std_out,message,'COLL')
   !    end do
   !  endif
   !enddo
   !do iatom=1,natom
   !  if(matlu(iatom)%lpawu.ne.-1) then
   !    write(std_out,*) "gathermatlu in rotate_matlu 2"
   !    do im1=1,tndim
   !      write(message,'(12(1x,18(1x,"(",f17.10,",",f17.10,")")))')&
   !&       (gathermatlu(iatom)%value(im1,im2),im2=1,tndim)
   !      call wrtout(std_out,message,'COLL')
   !    end do
   !  endif
   !enddo
   ABI_DEALLOCATE(temp_mat)
     !MSG_ERROR("Aborting now")

! Choose inverse rotation: reconstruct correct rot_mat from rot_mat_orig
! ========================================================================
   if(inverse==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         tndim=nspinor*(2*matlu(iatom)%lpawu+1)
           do im1=1,tndim
             do im2=1,tndim
!               rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom,isppol)%value(im1,im2)
               rot_mat(iatom,isppol)%value(im1,im2)=rot_mat_orig(iatom)%value(im1,im2)
             enddo
           enddo
       endif ! lpawu
     enddo ! iatom
   endif

! ===========================
! Put data into matlu(iatom)
! ===========================
   if(nsppol==1.and.nspinor==2) then
     call gather_matlu(matlu,gathermatlu,natom,option=-1,prtopt=1)
   else if((nsppol==2.or.nsppol==1).and.nspinor==1) then
     do iatom=1,natom
       if(matlu(iatom)%lpawu.ne.-1) then
         do im1=1,tndim
           do im2=1,tndim
             matlu(iatom)%mat(im1,im2,isppol,1,1)= gathermatlu(iatom)%value(im1,im2)
           enddo
         enddo
       endif
     enddo
   endif ! test nsppol/nspinor
! ===========================
! Deallocations
! ===========================
   do iatom=1,natom
     if(matlu(iatom)%lpawu.ne.-1) then
       ABI_DEALLOCATE(gathermatlu(iatom)%value)
!       ABI_DEALLOCATE(rot_mat_orig(iatom,isppol)%value)
       ABI_DEALLOCATE(rot_mat_orig(iatom)%value)
     endif
   enddo
   ABI_DATATYPE_DEALLOCATE(gathermatlu)
   ABI_DATATYPE_DEALLOCATE(rot_mat_orig)
 enddo ! isppol

 end subroutine rotate_matlu
!!***

!!****f* m_matlu/shift_matlu
!! NAME
!! shift_matlu
!!
!! FUNCTION
!! shift matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!  shift= shift of the diagonal part.
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: shifted matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_green,m_self,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine shift_matlu(matlu,natom,shift,signe)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 complex(dpc),intent(in) :: shift(natom)
 integer, optional,intent(in) :: signe
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol
 integer :: lpawu,signe_used
! character(len=500) :: message
!arrays
!************************************************************************
 signe_used=1
 if(present(signe)) then
   if(signe==-1) signe_used=-1
 endif
 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     do im=1,2*lpawu+1
       do isppol=1,matlu(1)%nsppol
         do ispinor=1,matlu(1)%nspinor
           matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=&
&           matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)-signe_used*shift(iatom)
         enddo
       enddo
     enddo
   endif
 enddo

 end subroutine shift_matlu
!!***

!!****f* m_matlu/checkreal_matlu
!! NAME
!! checkreal_matlu
!!
!! FUNCTION
!! Check that matlu is real in the orbital index with given precision
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!  tol : precision
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine checkreal_matlu(matlu,natom,tol)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 integer, intent(in) :: natom
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,im1,ispinor,ispinor1,isppol
 integer :: lpawu
 character(len=500) :: message
 real(dp) :: maximag,maxoffdiag,maximagdiag
!arrays
!************************************************************************
 maximag=zero
 maximagdiag=zero
 maxoffdiag=zero
 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     do im=1,2*lpawu+1
       do im1=1,2*lpawu+1
         do isppol=1,matlu(1)%nsppol
           do ispinor=1,matlu(1)%nspinor
             do ispinor1=1,matlu(1)%nspinor
               if(abs(aimag(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)))>maximag) then
                 maximag=abs(aimag(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)))
                 if(im==im1.and.ispinor==ispinor1) maximagdiag=abs(aimag(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)))
               endif
               if((im/=im1.or.ispinor/=ispinor1).and.abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))>maxoffdiag) then
                 maxoffdiag=abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               endif
             enddo
           enddo ! ispinor
         enddo ! isppol
       enddo ! im1
     enddo ! im
   endif ! lpawu
 enddo ! iatom
 if (maximagdiag>tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
&   ' Diagonal part of the occupation matrix is complex: the imaginary part ',&
&     maximagdiag,' is larger than',tol,ch10  &
&    , "The calculation cannot handle it : check that your calculation is meaningfull"
   MSG_ERROR(message)
 endif
 if (maximag>tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
&   ' Off diag occupation matrix is complex: the imaginary part ',maximag,' is larger than',tol,ch10&
    , "Check that your calculation is meaningfull"
   MSG_WARNING(message)
 else
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
&   ' Occupation matrix is real: the imaginary part ',maximag,' is lower than',tol
   MSG_COMMENT(message)
 endif
 if (maxoffdiag>tol) then
   write(message,'(3x,2a,e12.4,a,e12.4,6a)') ch10,&
&   ' Occupation matrix is non diagonal : the maximum off-diag part ',maxoffdiag,' is larger than',tol,ch10&
&    , "The corresponding non diagonal elements will be neglected in the Weiss/Hybridization functions",ch10&
&    , "(Except if dmft_solv=8,9 where these elements are taken into accounts)",ch10&
&    , "This is an approximation"
   MSG_WARNING(message)
 else
   write(message,'(3x,2a,e12.4,a,e12.4,2a)') ch10,&
&   ' Occupation matrix is diagonal : the off-diag part ',maxoffdiag,' is lower than',tol
   MSG_COMMENT(message)
 endif

 end subroutine checkreal_matlu
!!***

!!****f* m_matlu/checkdiag_matlu
!! NAME
!! checkdiag_matlu
!!
!! FUNCTION
!! Check that matlu is diagonal in the orbital index with given precision
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!  tol : precision
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      compute_levels
!!
!! CHILDREN
!!
!! SOURCE
 subroutine checkdiag_matlu(matlu,natom,tol,nondiag)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: tol
 integer, intent(in) :: natom
 logical, intent(out) :: nondiag
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,im1,ispinor,ispinor1,isppol
 integer :: lpawu
!arrays
!************************************************************************
 nondiag=.false.
 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     do im=1,2*lpawu+1
       do im1=1,2*lpawu+1
         do isppol=1,matlu(1)%nsppol
           do ispinor=1,matlu(1)%nspinor
!              if(im/=im1) write(std_out,*) "im,im1",im,im1,matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor)
 !            if(present(nondiag).eqv..false.) then
 !              if(im/=im1.and.(abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor))>tol))  then
 !                write(message,'(5i5)') im,im1,isppol,ispinor,ispinor
 !                call wrtout(std_out,message,'COLL')
 !                write(message,'(a,3e16.5)')" checkdiag_matlu: Warning ",matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor),tol
 !                call wrtout(std_out,message,'COLL')
 !                if(.not.present(opt)) MSG_ERROR("not present(opt)")
 !                if(matlu(1)%nspinor==1) MSG_ERROR("matlu%nspinor==1")
 !              endif
!             endif
             do ispinor1=1,matlu(1)%nspinor
 !              if(present(nondiag)) then
                 if((im/=im1.or.ispinor/=ispinor1)&
&                         .and.abs(real(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)))>tol) then
                   nondiag=.true.
                  ! write(6,*) "NONDIAG", matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                 endif
               !if(ispinor/=ispinor1.and.(abs(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))>tol))  then
               !  write(message,'(a,3e16.5)')" checkdiag_matlu :i Warning ",matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1),tol
               !  call wrtout(std_out,message,'COLL')
               !  write(message,'(5i5)') im,im1,isppol,ispinor,ispinor
               !  call wrtout(std_out,message,'COLL')
               !  if(matlu(1)%nspinor==1) MSG_ERROR("matlu%nspinor==1")
               !endif
             enddo
           enddo ! ispinor
         enddo ! isppol
       enddo ! im1
     enddo ! im
   endif ! lpawu
 enddo ! iatom

 end subroutine checkdiag_matlu
!!***

!!****f* m_matlu/prod_matlu
!! NAME
!! prod_matlu
!!
!! FUNCTION
!! Do the matrix product of two matlus
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity
!!  matlu2(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity
!!
!! OUTPUT
!!  matlu3(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: output quantity
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_oper,qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine prod_matlu(matlu1,matlu2,matlu3,natom)
 use defs_basis
 use defs_wvltypes
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type), intent(in) :: matlu1(natom),matlu2(natom)
 type(matlu_type), intent(inout) :: matlu3(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,im3,ispinor1,ispinor2,ispinor3,isppol
 integer :: lpawu
!arrays
!************************************************************************
 call zero_matlu(matlu3,natom)
 do iatom=1,natom
   lpawu=matlu1(iatom)%lpawu
   if(lpawu.ne.-1) then
     do isppol=1,matlu1(1)%nsppol
       do ispinor1=1,matlu1(1)%nspinor
         do ispinor2=1,matlu1(1)%nspinor
           do ispinor3=1,matlu1(1)%nspinor
             do im1=1,2*lpawu+1
               do im2=1,2*lpawu+1
                 do im3=1,2*lpawu+1
                   matlu3(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)= &
&                    matlu3(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)+ &
&                    matlu1(iatom)%mat(im1,im3,isppol,ispinor1,ispinor3)*&
&                    matlu2(iatom)%mat(im3,im2,isppol,ispinor3,ispinor2)
                 enddo ! im3
               enddo ! im2
             enddo ! im1
           enddo ! ispinor3
         enddo ! ispinor2
       enddo ! ispinor1
     enddo ! isppol
   endif ! lpawu
 enddo ! iatom

 end subroutine prod_matlu
!!***

!!****f* m_matlu/conjg_matlu
!! NAME
!! conjg_matlu
!!
!! FUNCTION
!! conjugate of input matlu
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine conjg_matlu(matlu1,natom)
 use defs_basis
 use defs_wvltypes
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type), intent(inout) :: matlu1(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im1,im2,ispinor2,ispinor1,isppol
 integer :: lpawu
!arrays
!************************************************************************
 do iatom=1,natom
   lpawu=matlu1(iatom)%lpawu
   if(lpawu.ne.-1) then
     do isppol=1,matlu1(1)%nsppol
       do ispinor1=1,matlu1(1)%nspinor
         do ispinor2=1,matlu1(1)%nspinor
           do im1=1,2*lpawu+1
             do im2=1,2*lpawu+1
               matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2)= &
&               conjg(matlu1(iatom)%mat(im1,im2,isppol,ispinor1,ispinor2))
             enddo ! im2
           enddo ! im1
         enddo ! ispinor2
       enddo ! ispinor1
     enddo ! isppol
   endif ! lpawu
 enddo ! iatom

 end subroutine conjg_matlu
!!***

!!****f* m_matlu/ln_matlu
!! NAME
!! ln_matlu
!!
!! FUNCTION
!! Compute the logarithm of matlu (only if diagonal for the moment)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity

!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
 subroutine ln_matlu(matlu1,natom)
 use defs_basis
 use defs_wvltypes
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type), intent(inout) :: matlu1(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol
 integer :: lpawu
 character(len=500) :: message
!arrays
!************************************************************************
 !call checkdiag_matlu(matlu1,natom,tol8)
 do iatom=1,natom
   lpawu=matlu1(iatom)%lpawu
   if(lpawu.ne.-1) then
     do isppol=1,matlu1(1)%nsppol
       do ispinor=1,matlu1(1)%nspinor
         do im=1,2*lpawu+1
           if( real(matlu1(iatom)%mat(im,im,isppol,ispinor,ispinor))<zero) then
             write(message,'(2a,2es13.5,a)') ch10," ln_matlu: PROBLEM " &
&             , matlu1(iatom)%mat(im,im,isppol,ispinor,ispinor)
             MSG_ERROR(message)
           endif
           matlu1(iatom)%mat(im,im,isppol,ispinor,ispinor)= &
&           log(matlu1(iatom)%mat(im,im,isppol,ispinor,ispinor))
         enddo ! im
       enddo ! ispinor
     enddo ! isppol
   endif ! lpawu
 enddo ! iatom

 end subroutine ln_matlu
!!***

!!****f* m_matlu/slm2ylm_matlu
!! NAME
!! slm2ylm_matlu
!!
!! FUNCTION
!! Transform mat from Slm to Ylm basis or vice versa
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu1(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity
!!  natom :: number of atoms
!!  option=1 go from Slm to Ylm basis
!!  option=2 go from Ylm to Slm basis
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine slm2ylm_matlu(matlu,natom,option,optprt)
 use defs_basis
 use defs_wvltypes
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,option,optprt
!arrays
 type(matlu_type), intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol,ispinor2
 integer :: lpawu,ll,mm,jm,ii,jj,im1,im2
 character(len=500) :: message
 real(dp) :: onem
 complex(dpc),allocatable :: slm2ylm(:,:)
 complex(dpc),allocatable :: mat_inp_c(:,:)
 complex(dpc),allocatable :: mat_out_c(:,:)
 complex(dpc) :: tmp2
 real(dp),parameter :: invsqrt2=one/sqrt2
!arrays
!************************************************************************

 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     ll=lpawu
     ABI_ALLOCATE(slm2ylm,(2*ll+1,2*ll+1))
     slm2ylm=czero
     do im=1,2*ll+1
       mm=im-ll-1;jm=-mm+ll+1
       onem=dble((-1)**mm)
       if (mm> 0) then
         slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
         slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
       end if
       if (mm==0) then
         slm2ylm(im,im)=cone
       end if
       if (mm< 0) then
         slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
         slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
       end if
     end do
     if(optprt>2) then
       write(message,'(2a)') ch10,"SLM2YLM matrix"
       call wrtout(std_out,message,'COLL')
       do im1=1,ll*2+1
         write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&         (slm2ylm(im1,im2),im2=1,ll*2+1)
         call wrtout(std_out,message,'COLL')
       end do
     endif
     do isppol=1,matlu(1)%nsppol
       do ispinor=1,matlu(1)%nspinor
         do ispinor2=1,matlu(1)%nspinor
           ABI_ALLOCATE(mat_out_c,(2*ll+1,2*ll+1))
           ABI_ALLOCATE(mat_inp_c,(2*ll+1,2*ll+1))
           mat_inp_c(:,:) = matlu(iatom)%mat(:,:,isppol,ispinor,ispinor2)
           mat_out_c=czero

           if(optprt>2) then
             write(message,'(2a, i2, a, i2, a, i2)') ch10,"SLM input matrix, isppol=", isppol, ", ispinor=", ispinor,& 
&             ", ispinor2=", ispinor2
             call wrtout(std_out,message,'COLL')
             do im1=1,ll*2+1
               write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&               (mat_inp_c(im1,im2),im2=1,ll*2+1)
               call wrtout(std_out,message,'COLL')
             end do
           endif

           do jm=1,2*ll+1
             do im=1,2*ll+1
               tmp2=czero
               do ii=1,2*ll+1
                 do jj=1,2*ll+1
                   if(option==1) then
                     tmp2=tmp2+mat_inp_c(ii,jj)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))
                   else if(option==2) then
                     tmp2=tmp2+mat_inp_c(ii,jj)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))
                   end if
                 end do
               end do
               mat_out_c(im,jm)=tmp2
             end do
           end do

           if(optprt>2) then
             write(message,'(2a, i2, a, i2, a, i2)') ch10,"YLM output matrix, isppol=", isppol, ", ispinor=", ispinor,&
&             ", ispinor2=", ispinor2
             call wrtout(std_out,message,'COLL')
             do im1=1,ll*2+1
               write(message,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
      &         (mat_out_c(im1,im2),im2=1,ll*2+1)
               call wrtout(std_out,message,'COLL')
             end do
           endif

           matlu(iatom)%mat(:,:,isppol,ispinor,ispinor2)=mat_out_c(:,:)
           ABI_DEALLOCATE(mat_out_c)
           ABI_DEALLOCATE(mat_inp_c)
         enddo ! im
       enddo ! ispinor
     enddo ! isppol
     ABI_DEALLOCATE(slm2ylm)
   endif ! lpawu
 enddo ! iatom


 end subroutine slm2ylm_matlu
!!***

!!****f* m_matlu/fac_matlu
!! NAME
!! fac_matlu
!!
!! FUNCTION
!! shift matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!  shift= shift of the diagonal part.
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: shifted matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine fac_matlu(matlu,natom,fac)
 use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 complex(dpc),intent(in)      :: fac
!Local variables-------------------------------
!scalars
 integer :: iatom,im,im1,ispinor,ispinor1,isppol
 integer :: lpawu
! character(len=500) :: message
!arrays
!************************************************************************
 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     do im=1,2*lpawu+1
       do im1=1,2*lpawu+1
         do isppol=1,matlu(1)%nsppol
           do ispinor=1,matlu(1)%nspinor
             do ispinor1=1,matlu(1)%nspinor
               matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=&
&                matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)*fac
             enddo
           enddo
         enddo
       enddo
     enddo
   endif
 enddo

 end subroutine fac_matlu
!!***

!!****f* m_matlu/printplot_matlu
!! NAME
!! printplot_matlu
!!
!! FUNCTION
!! shift matlu matrix
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!  shift= shift of the diagonal part.
!!
!! OUTPUT
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: shifted matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine printplot_matlu(matlu,natom,freq,char1,units,imre)
 use defs_basis
 use m_fstrings,       only : int2char4
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,units
 integer, optional, intent(in) :: imre
 real(dp), intent(in) :: freq
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
 character(len=*), intent(in) :: char1
!Local variables-------------------------------
!scalars
 integer :: iatom,im,im1,ispinor,ispinor1,isppol
 integer :: lpawu,unitnb
 character(len=4) :: tag_at
 character(len=fnlen) :: tmpfil,tmpfilre,tmpfilim
! character(len=500) :: message
!arrays
!************************************************************************
 ! not yet tested and used
 do iatom=1,natom
   lpawu=matlu(iatom)%lpawu
   if(lpawu.ne.-1) then
     unitnb=units+iatom
     call int2char4(iatom,tag_at)
     if(present(imre)) then
       tmpfilre = trim(char1)//tag_at//"re"
       tmpfilim = trim(char1)//tag_at//"im"
       open (unit=unitnb+10,file=trim(tmpfilre),status='unknown',form='formatted')
       open (unit=unitnb+20,file=trim(tmpfilim),status='unknown',form='formatted')
       write(unitnb+10,'(400e26.16)') freq,(((((real(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
       im=1,2*lpawu+1),ispinor=1,matlu(1)%nspinor),im1=1,2*lpawu+1),ispinor1=1,matlu(1)%nspinor),isppol=1,matlu(1)%nsppol)
       write(unitnb+20,'(400e26.16)') freq,(((((aimag(matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
       im=1,2*lpawu+1),ispinor=1,matlu(1)%nspinor),im1=1,2*lpawu+1),ispinor1=1,matlu(1)%nspinor),isppol=1,matlu(1)%nsppol)
     else
       tmpfil = trim(char1)//tag_at
       open (unit=unitnb,file=trim(tmpfil),status='unknown',form='formatted')
       write(unitnb,'(400e26.16)') freq,(((((matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1),&
       im=1,2*lpawu+1),ispinor=1,matlu(1)%nspinor),im1=1,2*lpawu+1),ispinor1=1,matlu(1)%nspinor),isppol=1,matlu(1)%nsppol)
     endif
   endif
 enddo
 end subroutine printplot_matlu
!!***

!!****f* m_matlu/identity_matlu
!! NAME
!! identity_matlu
!!
!! FUNCTION
!! Make the matlu the identity
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matlu(natom)%(nsppol,nspinor,nspinor,ndim,ndim) :: input quantity to rotate
!!  natom=number of atoms in cell.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!
!! SOURCE
 subroutine identity_matlu(matlu,natom)
 use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom
!arrays
 type(matlu_type),intent(inout) :: matlu(natom)
!Local variables-------------------------------
!scalars
 integer :: iatom,im,ispinor,isppol
 integer :: ndim
! character(len=500) :: message
!arrays
!************************************************************************
 ! not yet tested and used
 do iatom=1,natom
   if(matlu(iatom)%lpawu.ne.-1) then
     ndim=2*matlu(iatom)%lpawu+1
     do isppol=1,matlu(iatom)%nsppol
       do im=1,ndim
         do ispinor=1,matlu(iatom)%nspinor
           matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)= cone
         enddo ! ispinor
       enddo ! im
     enddo
   endif ! lpawu
 enddo ! iatom
 end subroutine identity_matlu
!!***

END MODULE m_matlu
!!***
