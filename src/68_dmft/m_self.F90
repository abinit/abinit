!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_self
!! NAME
!!  m_self
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_self

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore

 use m_fstrings, only : int2char4
 use m_crystal,  only : crystal_t
 use m_hu, only : hu_type
 use m_io_tools, only : get_unit
 use m_pawang, only : pawang_type
 use m_paw_dmft, only : paw_dmft_type
 use m_matlu, only : matlu_type, copy_matlu,shift_matlu,diag_matlu,rotate_matlu,init_matlu,destroy_matlu,print_matlu,zero_matlu
 use m_oper, only : oper_type,init_oper,destroy_oper, loc_oper, print_oper
 use m_datafordmft, only : compute_levels

 implicit none

 private

 public :: alloc_self
 public :: initialize_self
 public :: destroy_self
 public :: print_self
 public :: rw_self
 public :: dc_self
 public :: new_self
 public :: make_qmcshift_self
 public :: selfreal2imag_self

!!***

!!****t* m_self/self_type
!! NAME
!!  self_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! SOURCE

 type, public :: self_type ! for each atom

  integer :: dmft_nwlo
  ! dmft frequencies

  character(len=4) :: w_type
  ! type of frequencies used

  integer :: nw
  ! dmft frequencies

  integer :: iself_cv
  ! integer for self convergence

  integer :: dmft_nwli
  ! dmft frequencies

  real(dp), pointer :: omega(:) => null()
  ! value of frequencies

  real(dp), allocatable :: qmc_shift(:)
  ! value of frequencies

  real(dp), allocatable :: qmc_xmu(:)
  ! value of frequencies

  type(oper_type), allocatable :: oper(:)
  ! self function  in different basis

  type(oper_type):: hdc
  ! self function  in different basis

 end type self_type
!!***

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_self/alloc_self
!! NAME
!! alloc_self
!!
!! FUNCTION
!!  Allocate variables used in type self_type.
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent LDA+DMFT calculations.
!!  opt_oper = 1  Allocate only quantities in the KS basis.
!!             2  Allocate only quantities in the local basis.
!!             3  Allocate quantities in both the KS and local basis.
!!  wtype = "real" Self energy will be computed for real frequencies
!!        = "imag" Self energy will be computed for imaginary frequencies
!!
!! OUTPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!! PARENTS
!!      m_self
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine alloc_self(self,paw_dmft,opt_oper,wtype)

!Arguments ------------------------------------
!scalars
!type
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, optional, intent(in) :: opt_oper
 character(len=4), optional :: wtype
!local variables ------------------------------------
 integer :: ifreq,optoper

!************************************************************************
 if(present(opt_oper)) then
   optoper=opt_oper
 else
   optoper=2
 endif
 if(present(wtype)) then
   self%w_type=wtype
 else
   self%w_type="imag"
 endif
 if(self%w_type=="imag") then
   self%nw=paw_dmft%dmft_nwlo
   self%omega=>paw_dmft%omega_lo
 else if(self%w_type=="real") then
   self%nw=size(paw_dmft%omega_r)
   self%omega=>paw_dmft%omega_r
 endif

 self%dmft_nwlo=paw_dmft%dmft_nwlo
 self%dmft_nwli=paw_dmft%dmft_nwli
 self%iself_cv=0

 call init_oper(paw_dmft,self%hdc,opt_ksloc=optoper)
 ABI_DATATYPE_ALLOCATE(self%oper,(self%nw))
 do ifreq=1,self%nw
  call init_oper(paw_dmft,self%oper(ifreq),opt_ksloc=optoper)
 enddo

 if(paw_dmft%dmft_solv==4) then
   ABI_ALLOCATE(self%qmc_shift,(paw_dmft%natom))
   ABI_ALLOCATE(self%qmc_xmu,(paw_dmft%natom))
   self%qmc_shift(:)=zero
   self%qmc_xmu(:)=zero
 endif

end subroutine alloc_self
!!***

!!****f* m_self/initialize_self
!! NAME
!! initialize_self
!!
!! FUNCTION
!!  Initialize self-energy.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=variables related to crystal structure
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent LDA+DMFT calculations.
!!  opt_read =  not used for the moment
!!  wtype = "real" Self energy will be computed for real frequencies
!!        = "imag" Self energy will be computed for imaginary frequencies
!!
!! OUTPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!!
!! PARENTS
!!      m_dmft,spectral_function
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine initialize_self(self,paw_dmft,wtype)

!Arguments ------------------------------------
!scalars
!type
 type(self_type), intent(inout) :: self
 type(paw_dmft_type), intent(inout) :: paw_dmft
 character(len=4), optional, intent(in) :: wtype
!local variables ------------------------------------
! character(len=500) :: message
 integer :: iatom,ifreq
 character(len=4) :: wtype2
!************************************************************************
 if(present(wtype)) then
   wtype2=wtype
 else
   wtype2="imag"
 endif


 call alloc_self(self,paw_dmft,opt_oper=2,wtype=wtype2) !  opt_oper=1 is not useful and not implemented
 do ifreq=1,self%nw
   do iatom=1,paw_dmft%natom
     self%oper(ifreq)%matlu(iatom)%mat=czero
   enddo
 enddo
! if(paw_dmft%dmft_rslf==1.and.opt_read==1) then
!   call rw_self(cryst_struc,self,mpi_enreg,paw_dmft,pawtab,pawprtvol=2,opt_rw=1)
! endif
! write(message,'(a,a)') ch10,"   Self-energy for large frequency is"
! call wrtout(std_out,message,'COLL')
! call print_matlu(self%oper(paw_dmft%dmft_nwlo)%matlu,  &
!&                 paw_dmft%natom,3)

end subroutine initialize_self
!!***

!!****f* m_self/destroy_self
!! NAME
!! destroy_self
!!
!! FUNCTION
!!  deallocate self
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,spectral_function
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine destroy_self(self)

!Arguments ------------------------------------
!scalars
 type(self_type),intent(inout) :: self

!local variables-------------------------------
 integer :: ifreq

! *********************************************************************
 if ( allocated(self%oper))       then
   do ifreq=1,self%nw
    call destroy_oper(self%oper(ifreq))
   enddo
   ABI_DATATYPE_DEALLOCATE(self%oper)
 end if

 call destroy_oper(self%hdc)
 if (allocated(self%qmc_shift)) then
   ABI_DEALLOCATE(self%qmc_shift)
 end if
 if (allocated(self%qmc_xmu))  then
   ABI_DEALLOCATE(self%qmc_xmu)
 end if
 self%omega => null()

end subroutine destroy_self
!!***

!!****f* m_self/print_self
!! NAME
!! print_self
!!
!! FUNCTION
!!  print occupations
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  option = 1 Do not print double counting.
!!           2 Print double counting
!!  paw_dmft <type(paw_dmft_type)> =  variables related to self-consistent LDA+DMFT calculations.
!!  prtopt = integer which precises the amount of printing in the subroutine called
!!
!! OUTPUT
!!  self <type(self_type)>= variables related to self-energy
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine print_self(self,prtdc,paw_dmft,prtopt)

!Arguments ------------------------------------
!type
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(self_type),intent(in) :: self
 character(len=*), intent(in) :: prtdc
 integer,intent(in) :: prtopt


!local variables-------------------------------
 character(len=500) :: message
! *********************************************************************

 write(message,'(2a)') ch10,"  == The self-energy for smallest frequency is   == "
 call wrtout(std_out,message,'COLL')
 call print_oper(self%oper(1),1,paw_dmft,prtopt)
! write(message,'(2a)') ch10,"  == The self-energy for small (3) frequency is   == "
! call wrtout(std_out,message,'COLL')
! call print_oper(self%oper(3),1,paw_dmft,prtopt)
 write(message,'(2a)') ch10,"  == The self-energy for large frequency is   == "
 call wrtout(std_out,message,'COLL')
 call print_oper(self%oper(self%nw),1,paw_dmft,prtopt)
 if(prtdc=="print_dc") then
   write(message,'(2a)') ch10,"  == The double counting hamiltonian is (diagonal part)  == "
   call wrtout(std_out,message,'COLL')
   call print_matlu(self%hdc%matlu,paw_dmft%natom,prtopt,opt_diag=1)
 endif

end subroutine print_self
!!***

!!****f* m_self/dc_self
!! NAME
!! dc_self
!!
!! FUNCTION
!!
!! INPUTS
!!  charge_loc(cryst_struc%natom,paw_dmft%nsppol+1)= total charge for correlated electrons on a given atom, and for spin
!!  cryst_struc <type(crystal_t)>=variables related to crystal structure
!!  hu <type(hu_type)>= variables related to the interaction between electrons
!!  self <type(self_type)>= variables related to self-energy
!!  dmft_dc = 1 Full localized Limit double counting.
!!           2 Around Mean Field (without SO)
!!           0 not double counting
!!  prtopt = integer which precises the amount of printing (not used here)
!!
!! OUTPUT
!!  self <type(self_type)>= variables related to self-energy
!!  hu <type(hu_type)>= variables related to the interaction between electrons
!!
!! PARENTS
!!      m_dmft,spectral_function
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine dc_self(charge_loc,cryst_struc,hu,self,dmft_dc,prtopt)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(self_type),intent(inout) :: self
 real(dp), intent(in) :: charge_loc(cryst_struc%natom,self%hdc%nsppol+1)
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 integer, intent(in) :: prtopt,dmft_dc

!Local variables-------------------------------
 integer :: iatom,isppol,ispinor,lpawu,m1,nspinor,nsppol
 real(dp) :: ntot,jpawu_dc
 character(len=500) :: message
! *********************************************************************
 nspinor = self%hdc%nspinor
 nsppol  = self%hdc%nsppol

 do iatom=1,self%hdc%natom
   lpawu=self%hdc%matlu(iatom)%lpawu
   ntot=charge_loc(iatom,nsppol+1)
!     nup=charge_loc(iatom,1)
!     ndn=charge_loc(iatom,nsppol)
   if(lpawu/=-1) then
     self%hdc%matlu(iatom)%mat=czero
     do isppol = 1 , nsppol
       do ispinor = 1 , nspinor
         do m1=1, 2*lpawu+1
           if(dmft_dc==1.or.dmft_dc==4.or.dmft_dc==5) then ! FLL
             if(dmft_dc==1.or.dmft_dc==5) then
               jpawu_dc=hu(cryst_struc%typat(iatom))%jpawu
             else if (dmft_dc==4) then ! FLL
               jpawu_dc=zero
             endif
             if(nspinor==2.or.dmft_dc==5) then
               self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor) =  &
&                hu(cryst_struc%typat(iatom))%upawu * ( ntot - one / two )     &
&                - one/two * jpawu_dc * ( ntot - one       )
             else
               self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor)=  &
         &       hu(cryst_struc%typat(iatom))%upawu * ( ntot - one / two ) &
!&     -  hu(cryst_struc%typat(iatom))%jpawu * ( charge_loc(iatom,2-nsppol+1) - one)
&                 -  jpawu_dc * ( charge_loc(iatom,isppol) - one / two )
             endif
           else if(dmft_dc==2) then  ! AMF
             if(nspinor==2) then
               !write(message,'(a,i4,i4,2x,e20.10)') " AMF Double counting not implemented for SO"
               !MSG_ERROR(message)
               self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor)= &
               hu(cryst_struc%typat(iatom))%upawu * ntot/two &
               + ( hu(cryst_struc%typat(iatom))%upawu - hu(cryst_struc%typat(iatom))%jpawu )&
               *ntot/two*(float(2*lpawu))/(float(2*lpawu+1))
               write(message,'(a,i4,i4,2x,e20.10)') " AMF Double counting is under test for SOC"
               MSG_WARNING(message)
             else
               if(nsppol==2) then
                 self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor)=  &
&                  hu(cryst_struc%typat(iatom))%upawu * charge_loc(iatom,2-isppol+1) &
&                  +  (hu(cryst_struc%typat(iatom))%upawu - hu(cryst_struc%typat(iatom))%jpawu )&
&                  *charge_loc(iatom,isppol)*(float(2*lpawu))/(float(2*lpawu+1))
                else  if(nsppol==1) then
                  self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor)=  &
&                   hu(cryst_struc%typat(iatom))%upawu * charge_loc(iatom,isppol) &
&                    +  (hu(cryst_struc%typat(iatom))%upawu - hu(cryst_struc%typat(iatom))%jpawu )&
&                   *charge_loc(iatom,isppol)*(float(2*lpawu))/(float(2*lpawu+1))
                endif
!                 write(std_out,*) "AMF",  charge_loc(iatom,2-isppol+1)
!                 write(std_out,*) "AMF",  charge_loc(iatom,isppol+1)
!                 write(std_out,*) "AMF",  lpawu
!                 write(std_out,*) "AMF",  hu(cryst_struc%typat(iatom))%upawu
!                 write(std_out,*) "AMF", self%hdc%matlu(iatom)%mat(m1,m1,isppol,ispinor,ispinor)
             endif
           else
             MSG_ERROR("not implemented")
           endif
         enddo  ! m1
       enddo  ! ispinor
     enddo ! isppol
   endif
 enddo ! iatom

 if(prtopt>0) then
 endif


end subroutine dc_self
!!***

!!****f* m_self/rw_self
!! NAME
!! rw_self
!!
!! FUNCTION
!!
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  prtopt = integer which precises the amount of printing
!!  opt_rw = 1  Read Self-Energy.
!!           2  Write Self-Energy.
!!           3  Impose Self-Energy.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft,spectral_function
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine rw_self(self,paw_dmft,prtopt,opt_rw,istep_iter,opt_char,opt_imagonly,opt_selflimit,opt_hdc,opt_stop,pawang,cryst_struc)

!Arguments ------------------------------------
!type
 type(self_type),intent(inout) :: self
 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer,intent(in) :: prtopt
 integer,intent(in),optional :: opt_rw,istep_iter,opt_imagonly
 character(len=4), optional :: opt_char
 integer, intent(in), optional :: opt_stop
 type(matlu_type), optional, intent(inout) :: opt_selflimit(paw_dmft%natom)
 type(matlu_type), optional, intent(in) :: opt_hdc(paw_dmft%natom)
 type(pawang_type), optional, intent(in) :: pawang
 type(crystal_t), optional, intent(in) :: cryst_struc

!local variables-------------------------------
 type(coeff2c_type), allocatable :: eigvectmatlu(:,:)
 type(matlu_type), allocatable :: level_diag(:)
 type(matlu_type), allocatable :: selfrotmatlu(:)
 type(oper_type)  :: energy_level
 logical :: lexist
 complex(dpc), allocatable :: buffer(:)
 integer :: iall,iatom,iatu,ier,iexist2,ifreq,im,im1,ioerr,ispinor,ispinor1,isppol,istepiter,istep,istep_imp
 integer :: icount,iexit,iter,iter_imp,master,mbandc,myproc,natom,ncount,ndim,nkpt,nproc,nrecl,nspinor,nsppol,spacecomm
 integer :: natom_read,nsppol_read,nspinor_read,ndim_read,nw_read,optrw,readimagonly,tndim,iflavor,unitrot
 logical :: nondiaglevels
 character(len=30000) :: message ! Big buffer to avoid buffer overflow.
 integer,allocatable :: unitselffunc_arr(:)
 integer,allocatable :: unitselfrot(:,:,:,:)
 character(len=fnlen) :: tmpfil,tmpfilrot,tmpmatrot
 character(len=1) :: tag_is
 character(len=10) :: tag_iflavor
 character(len=10) :: tag_at
 character(len=4) :: chtemp
 real(dp):: xtemp,fermie_read,x_r,x_i
 real(dp), allocatable:: s_r(:,:,:,:),s_i(:,:,:,:),fermie_read2(:)
! *********************************************************************

! Initialise spaceComm, myproc, and nproc
 istep=0 ; iter=0 ; istep_imp=0 ; iter_imp=0
 if(present(opt_rw)) then
   optrw=opt_rw
 else
   optrw=0
 endif
 readimagonly=0
 if(present(opt_imagonly)) then
   if(opt_imagonly==1.and.paw_dmft%dmft_solv>=5) then
     readimagonly=opt_imagonly
     write(message,*)
     write(message,'(4x,2a)') "About to read imaginary part of Self energy"
     call wrtout(std_out,message,'COLL')
   else
     readimagonly=0
     write(message,'(4x,2a)') "About to read both real and imaginary part of Self energy"
     call wrtout(std_out,message,'COLL')
   endif
 else
   readimagonly=0
 endif
 if(present(istep_iter)) then
   istepiter=istep_iter
 else
   istepiter=0
 endif
 if(paw_dmft%use_fixed_self>0) then
   istep=istepiter/1000
   iter=istepiter-(istepiter/1000)*1000
   istep_imp=paw_dmft%use_fixed_self/1000
   iter_imp=paw_dmft%use_fixed_self-(paw_dmft%use_fixed_self/1000)*1000
 endif

 if(paw_dmft%dmft_rslf<=0.and.optrw==1) optrw=0
 iexit=0
 ioerr=0
 iexist2=1
 lexist=.true.
 spacecomm=paw_dmft%spacecomm
 myproc=paw_dmft%myproc
 nproc=paw_dmft%nproc
 master=0

! write(std_out,*) "myproc,master",myproc,master
 if(prtopt>200) then
 endif
 natom=self%oper(1)%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt

!   - For the Tentative rotation of the self-energy file (begin init)
 if(present(pawang)) then
   ABI_ALLOCATE(unitselfrot,(natom,nsppol,nspinor,7)) ! 7 is the max ndim possible
   if(optrw==2) then
     write(message,'(a,2x,a,f13.5)') ch10,&
&     " == About to print self-energy for MAXENT code "
   else if (optrw==1)  then
     write(message,'(a,2x,a,f13.5)') ch10,&
&     " == About to read self-energy from MAXENT code "
   endif
   call wrtout(std_out,message,'COLL')

   ABI_DATATYPE_ALLOCATE(eigvectmatlu,(natom,nsppol))
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom)/=-1) then
       tndim=nspinor*(2*paw_dmft%lpawu(iatom)+1)
       do isppol=1,nsppol
         ABI_ALLOCATE(eigvectmatlu(iatom,isppol)%value,(tndim,tndim))
       end do
     end if
   end do
   ABI_DATATYPE_ALLOCATE(level_diag,(natom))
   ABI_DATATYPE_ALLOCATE(selfrotmatlu,(natom))
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,level_diag)
   call init_matlu(natom,nspinor,nsppol,paw_dmft%lpawu,selfrotmatlu)
   call init_oper(paw_dmft,energy_level,opt_ksloc=3)
 endif
!   - For the Tentative rotation of the self-energy file (end init)

!   - For the Tentative rotation of the self-energy file (begin diag)
 if(optrw==2.and.present(pawang)) then
   call compute_levels(cryst_struc,energy_level,self%hdc,pawang,paw_dmft,nondiag=nondiaglevels)
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print not Diagonalized Self Energy for Fermi Level=",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
   call print_matlu(self%oper(2)%matlu,natom,1,compl=1,opt_exp=1)
   call diag_matlu(energy_level%matlu,level_diag,natom,&
&   prtopt=prtopt,eigvectmatlu=eigvectmatlu,&
&   test=paw_dmft%dmft_solv)
   write(message,'(a,2x,a,f13.5)') ch10,&
&   " == Print Diagonalized levels for Fermi Level=",paw_dmft%fermie
   call wrtout(std_out,message,'COLL')
   call print_matlu(level_diag,natom,1,compl=1,opt_exp=1)
   !  Create file for rotation
 endif
 if(present(pawang)) then
   do iatom=1,natom
     if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
       call int2char4(iatom,tag_at)
       ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
       if(optrw==2) then
         tmpmatrot = trim(paw_dmft%filapp)//'.UnitaryMatrix_for_DiagLevel_iatom'//trim(tag_at)
       else if (optrw==1) then
         tmpmatrot = trim(paw_dmft%filnamei)//'.UnitaryMatrix_for_DiagLevel_iatom'//trim(tag_at)
       endif
       unitrot=3189+iatom
#ifdef FC_NAG
       open (unit=unitrot,file=trim(tmpmatrot),status='unknown',form='formatted',recl=ABI_RECL)
#else
       open (unit=unitrot,file=trim(tmpmatrot),status='unknown',form='formatted')
#endif
       write(std_out,*) "     Open file  ",trim(tmpmatrot)
       rewind(unitrot)
       ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
       do isppol=1,nsppol
         if(optrw==2) then
           do im=1,ndim
             do im1=1,ndim
               write(message,*) real(eigvectmatlu(iatom,isppol)%value(im,im1)),aimag(eigvectmatlu(iatom,isppol)%value(im,im1))
               call wrtout(unitrot,message,'COLL')
             enddo
           enddo
         else if (optrw==1) then
           do im=1,ndim
             do im1=1,ndim
               read(unitrot,*) x_r,x_i
               eigvectmatlu(iatom,isppol)%value(im,im1)=cmplx(x_r,x_i)
             enddo
           enddo
         endif
       enddo
       close(unitrot)
     endif
   enddo
   if(optrw==1) then
     write(message,'(a,2x,a,i4)') ch10,&
&      " == Print non rotated Self Limit read from Matsubara space=",ifreq
     call wrtout(std_out,message,'COLL')
     call print_matlu(opt_selflimit,natom,1,compl=1)
     call rotate_matlu(opt_selflimit,eigvectmatlu,natom,3,1)
     write(message,'(a,2x,a,i4)') ch10,&
&      " == Print rotated Self Limit read from Matsubara space=",ifreq
     call wrtout(std_out,message,'COLL')
     call print_matlu(opt_selflimit,natom,1,compl=1)
   endif
 endif
!   - For the Tentative rotation of the self-energy file (end diag)

 if((optrw==2.or.optrw==1).and.myproc==master)then
   ABI_ALLOCATE(unitselffunc_arr,(natom*nsppol*nspinor))
   iall=0
   do iatom=1,natom
     if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
       ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
       ABI_ALLOCATE(s_r,(ndim,ndim,nspinor,nspinor))
       ABI_ALLOCATE(s_i,(ndim,ndim,nspinor,nspinor))
!       write(std_out,*) "print_self",ndim
       call int2char4(iatom,tag_at)
       ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
       do isppol=1,nsppol
         write(tag_is,'(i1)')isppol
!         do ispinor=1,nspinor
           iall=iall+1
!           write(tag_is2,'(i1)')ispinor

!          ===========================
!           == Create name for file
!          ===========================
           if(self%w_type=="real") then
             if (optrw==1) then
               tmpfil = trim(paw_dmft%filnamei)//'Self_ra-omega_iatom'//trim(tag_at)//'_isppol'//tag_is
             else
               tmpfil = trim(paw_dmft%filapp)//'Self_ra-omega_iatom'//trim(tag_at)//'_isppol'//tag_is
             endif
           else
             tmpfil = trim(paw_dmft%filapp)//'Self-omega_iatom'//trim(tag_at)//'_isppol'//tag_is
             if(present(opt_char)) then
               tmpfil = trim(paw_dmft%filapp)//'Self_ra-omega_iatom'//trim(tag_at)//'_isppol'//tag_is//opt_char
             endif
           endif
           if(optrw==1) write(message,'(3a)') ch10,"  == Read  self function and Fermi Level on file ",trim(tmpfil)
           if(optrw==2) write(message,'(3a)') ch10,"  == Write self function and Fermi Level on file ",trim(tmpfil)
           call wrtout(std_out,message,'COLL')
           !unitselffunc_arr(iall)=300+iall-1
           unitselffunc_arr(iall) = get_unit()
           ABI_CHECK(unitselffunc_arr(iall) > 0, "Cannot find free IO unit!")

           !- For the Tentative rotation of the self-energy file (create file)
           if(optrw==2.and.present(pawang)) then
             iflavor=0
             do ispinor=1,nspinor
              do im=1,ndim
                iflavor=iflavor+1
                call int2char4(iflavor,tag_iflavor)
                unitselfrot(iatom,isppol,ispinor,im)=3000+iflavor
                ABI_CHECK(unitselfrot(iatom,isppol,ispinor,im) > 0, "Cannot find free IO unit for unitselfrot!")
                tmpfilrot = trim(paw_dmft%filapp)//'Selfrotformaxent'//&
                & trim(tag_at)//'_isppol'//tag_is//'_iflavor'//trim(tag_iflavor)
                write(std_out,*) "Create file  ",trim(tmpfilrot)," unit ",unitselfrot(iatom,isppol,ispinor,im)," for flavor",iflavor
#ifdef FC_NAG
                open (unit=unitselfrot(iatom,isppol,ispinor,im),file=trim(tmpfilrot),&
                & status='unknown',form='formatted',recl=ABI_RECL)
#else
                open (unit=unitselfrot(iatom,isppol,ispinor,im),file=trim(tmpfilrot),status='unknown',form='formatted')
#endif
                rewind(unitselfrot(iatom,isppol,ispinor,im))
              enddo
             enddo
           endif
           !- For the Tentative rotation of the self-energy file (create file)

!           write(std_out,*) "1"

!          ===========================
!           == Read: check that the file exists
!          ===========================
           if(optrw==1) then
!           write(std_out,*) "3"
             inquire(file=trim(tmpfil),exist=lexist,recl=nrecl)
             if((.not.lexist)) then
!           write(std_out,*) "4"
               iexist2=0
               write(message,'(4x,a,i5,3a)') "File number",unitselffunc_arr(iall),&
&               " called ",trim(tmpfil)," does not exist"
!               write(std_out,*) lexist,nrecl
               call wrtout(std_out,message,'COLL')
             endif
           endif
           !write(std_out,*) "2"

!          ===========================
!           == Open file
!          ===========================
           if(optrw==2.or.(optrw==1.and.iexist2==1)) then
             !write(std_out,*) "5"
#ifdef FC_NAG
             open (unit=unitselffunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted',recl=ABI_RECL)
#else
             open (unit=unitselffunc_arr(iall),file=trim(tmpfil),status='unknown',form='formatted')
#endif
             rewind(unitselffunc_arr(iall))
             !write(std_out,*) "61",nrecl
             if(prtopt>=3) then
               write(message,'(a,a,a,i4)') '    opened file : ', trim(tmpfil), ' unit', unitselffunc_arr(iall)
               call wrtout(std_out,message,'COLL')
             endif
           endif
           !write(std_out,*) "6",nrecl

!          ===========================
!           == Check Header
!          ===========================
           if(optrw==2) then
             write(message,'(3a,5i5,2x,e25.17)') "# natom,nsppol,nspinor,ndim,nw,fermilevel",ch10&
&             ,"####",natom,nsppol,nspinor,ndim,self%nw,paw_dmft%fermie
             call wrtout(unitselffunc_arr(iall),message,'COLL')
           else if(optrw==1.and.iexist2==1.and.readimagonly==0) then
             read(unitselffunc_arr(iall),*)
             read(unitselffunc_arr(iall),*,iostat=ioerr)&
&              chtemp,natom_read,nsppol_read,nspinor_read,ndim_read,nw_read,fermie_read
             if(ioerr<0) then
!              write(std_out,*)" HEADER IOERR"
!              write(std_out,'(a4,2x,31(e15.8,2x))') chtemp,natom_read,nsppol_read,nspinor_read,ndim_read,nw_read,fermie_read
             endif
             if(ioerr==0) then
               write(message,'(a,3x,3a,i12,2a,i11,2a,i10,2a,i13,2a,i11,2a,e25.8)') ch10,"Data in Self file corresponds to",&
&               ch10,"     natom",natom_read,&
&               ch10,"     nsppol",nsppol_read,&
&               ch10,"     nspinor",nspinor_read,&
&               ch10,"     ndim",ndim_read, &
&               ch10,"     nw",nw_read, &
&               ch10,"     Fermi level",fermie_read
               call wrtout(std_out,message,'COLL')
               if((natom/=natom_read).or.(nsppol_read/=nsppol).or.&
&                (nspinor/=nspinor_read).or.(nw_read/=self%nw)) then
                 write(message,'(a,3x,3a,i12,2a,i11,2a,i10,2a,i13,2a,i11,2a,e25.8)') ch10,"Data required is ",&
&                 ch10,"     natom",natom,&
&                 ch10,"     nsppol",nsppol,&
&                 ch10,"     nspinor",nspinor,&
&                 ch10,"     ndim",ndim, &
&                 ch10,"     nw",self%nw, &
&                 ch10,"     Fermi level",paw_dmft%fermie
                 call wrtout(std_out,message,'COLL')
                 message = "Dimensions in self are not correct"
                 if(readimagonly==1.or.present(opt_stop)) then
                   MSG_ERROR(message)
                 else
                   MSG_WARNING(message)
                 endif
                 iexist2=2
               endif
             else
               MSG_WARNING("Self file is empty")
             endif
           endif
           !write(std_out,*) "7"

!          ===========================
!           == Write/Read self in the file
!          ===========================

           rewind(111)
           do ifreq=1,self%nw
             if(optrw==2) then
!               write(std_out,'(a,2x,31(e15.8,2x))') &
!&              "SETEST",self%omega(ifreq),&
!&              (self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)&
!&               ,im=1,ndim)
!               write(std_out,*) self%omega(ifreq),&
!&              ((self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor)&
!&               ,im=1,ndim),im1=1,ndim)
!               write(message,'(2x,393(e25.17,2x))')  self%omega(ifreq),&
!&              ((((self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
!&              ,im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)

               !MGNAG Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output
               !Is it possible to rewrite the code below to avoid such a long message
               !What about Netcdf binary files ?
               if(nspinor==1) then
                 write(message,'(2x,393(e25.17,2x))')  self%omega(ifreq),&
&                ((((real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
&                aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
&                im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
               else
                 write(message,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
&                ((((real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
&                aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),&
&                im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
               endif
               call wrtout(unitselffunc_arr(iall),message,'COLL')

               !- For the Tentative rotation of the self-energy file (begin rot)
               !----------------------------------------------------------------
               if(optrw==2.and.present(pawang)) then
                 call copy_matlu(self%oper(ifreq)%matlu,selfrotmatlu,natom)
                 if(ifreq<3) then
                   write(message,'(a,2x,a,i4)') ch10,&
&                    " == Print non Rotated Self Energy for freq=",ifreq
                   call wrtout(std_out,message,'COLL')
                   call print_matlu(selfrotmatlu,natom,1,compl=1)
                 endif
                 call rotate_matlu(selfrotmatlu,eigvectmatlu,natom,3,1)
                 if(ifreq<3) then
                   write(message,'(a,2x,a,i4)') ch10,&
&                    " == Print Rotated Self Energy for freq=",ifreq
                   call wrtout(std_out,message,'COLL')
                   call print_matlu(selfrotmatlu,natom,1,compl=1)
                 else if(ifreq==3) then
                   write(message,'(a,2x,a,i4)') ch10,&
&                    "  (Other frequencies not printed)"
                   call wrtout(std_out,message,'COLL')
                 endif
                 !write(message,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
                !  ((real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
                !  aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
                !  im=1,ndim),ispinor=1,nspinor)
                 iflavor=0
                 do ispinor=1,nspinor
                   do im=1,ndim
                     iflavor=iflavor+1
                    ! if(ifreq<5) then
                    !   write(std_out,*) "Write in file unit",unitselfrot(iatom,isppol,ispinor,im),"for flavor",iflavor
                    ! endif
                     write(message,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
&                      real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
&                      aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor))
!                     write(6,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
!&                      real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
!&                      aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor))
                    ! if(iflavor==1) then
                    ! write(1024,*) iatom,isppol,ispinor,im,unitselfrot(iatom,isppol,ispinor,im)
                    ! write(1024,'(2x,393(e18.10,2x))')  self%omega(ifreq),&
                    ! &  real(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor)),&
                    ! &  aimag(selfrotmatlu(iatom)%mat(im,im,isppol,ispinor,ispinor))
                    ! endif
                     call wrtout(unitselfrot(iatom,isppol,ispinor,im),message,'COLL')
                   enddo
                 enddo
               endif
               !- For the Tentative rotation of the self-energy file (end rot)

!               write(std_out,*) unitselffunc_arr(iall)
             else if(optrw==1.and.iexist2==1.and.ioerr==0) then
           !write(std_out,*) "8"
!               read(unitselffunc_arr(iall),'(2x,31(e15.8,2x))',iostat=ioerr) &
!&              xtemp,(s_r(im),s_i(im),im=1,ndim)
               if(readimagonly==0) then
                 read(unitselffunc_arr(iall),*,iostat=ioerr) xtemp,&
&   ((((s_r(im,im1,ispinor,ispinor1),s_i(im,im1,ispinor,ispinor1),im=1,ndim),im1=1,ndim),ispinor=1,nspinor),ispinor1=1,nspinor)
!               if(ioerr<0) then
!                write(std_out,*)" SELF IOERR<"
!               else if(ioerr>0) then
!                write(std_out,*)" SELF IOERR>"
!                write(std_out,'(a4,2x,31(e15.8,2x))') xtemp,(s_r(im),s_i(im),im=1,ndim)
!               endif
                 do im=1,ndim
                   do im1=1,ndim
                     do ispinor=1,nspinor
                       do ispinor1=1,nspinor
                          self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
&                         =cmplx(s_r(im,im1,ispinor,ispinor1),s_i(im,im1,ispinor,ispinor1),kind=dp)
                 !if((im1==im1).and.(ispinor==ispinor1)) then
                 !  write(6,*) "Self read",ifreq, self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                 !endif
                       enddo
                     enddo
                   enddo
                 enddo
               endif
             endif
           enddo ! ifreq

           !- For the Tentative rotation of the self-energy file (begin close file)
           if(optrw==2.and.present(pawang)) then
             do ispinor=1,nspinor
               do im=1,ndim
                 close(unitselfrot(iatom,isppol,ispinor,im))
                 write(std_out,*) "Close file unit",unitselfrot(iatom,isppol,ispinor,im)
               enddo
             enddo
           endif
           !- For the Tentative rotation of the self-energy file (end close file)


           if(optrw==1.and.iexist2==1.and.ioerr==0) then
             if(readimagonly==1) then ! read from OmegaMaxent
               s_r=zero

               ! Read self energy from Maxent (imag part) on the real axis
               !----------------------------------------------------------
               do ifreq=1,self%nw
                 call zero_matlu(self%oper(ifreq)%matlu,natom)
               enddo
               do ispinor=1,nspinor
                 do im=1,ndim
                   do ifreq=1,self%nw
                     read(unitselffunc_arr(iall),*,iostat=ioerr) xtemp,s_i(im,im,ispinor,ispinor)
                      ! minus sign because - Im Sigma is the output of OmegaMaxent
                     self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)&
&                       =cmplx(s_r(im,im,ispinor,ispinor),-s_i(im,im,ispinor,ispinor),kind=dp)
                   enddo ! ifreq
                 enddo
               enddo

               ! Kramers Kronig
               !-------------------
               write(message,'(4x,2a)') " Read only diagonal self energy from Maxent"
               call wrtout(std_out,message,'COLL')
               !write(6,*) "opt_hdc",opt_hdc(1)%mat(1,1,1,1,1)
               call kramerskronig_self(self,opt_selflimit,opt_hdc,paw_dmft%filapp)

               ! Rotate back rotate_matlu
               !-----------------------------
               if(present(pawang)) then
                 write(message,'(4x,2a)') " Rotate Back self in the original basis"
                 call wrtout(std_out,message,'COLL')
                 do ifreq=1,self%nw
                   if(ifreq<20) then
                     write(message,'(a,2x,a,i4)') ch10,&
&                      " == Print Rotated real axis Self Energy for freq=",ifreq
                     call wrtout(std_out,message,'COLL')
                     call print_matlu(self%oper(ifreq)%matlu,natom,1,compl=1)
                   endif
                     call rotate_matlu(self%oper(ifreq)%matlu,eigvectmatlu,natom,3,-1)
                   if(ifreq<20) then
                     write(message,'(a,2x,a,i4)') ch10,&
&                      " == Print Rotated back real axis Self Energy for freq=",ifreq
                     call wrtout(std_out,message,'COLL')
                     call print_matlu(self%oper(ifreq)%matlu,natom,1,compl=1)
                   endif
                 enddo
               endif

             endif
           endif

!          ===========================
!           == Write/Read hdc in the file
!          ===========================
           if(optrw==2) then
!             write(std_out,'(a,2x,31(e15.8,2x))') &
!&            "SETEST #dc ",(self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim)
             write(message,'(a,2x,31(e25.17,2x))') &
&            "#dc ",((self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor),im=1,ndim),ispinor=1,nspinor)
             call wrtout(unitselffunc_arr(iall),message,'COLL')
           else if(optrw==1.and.iexist2==1.and.ioerr==0.and.readimagonly==0) then
         !write(std_out,*) "8"
             read(unitselffunc_arr(iall),*,iostat=ioerr) &
&             chtemp,((s_r(im,1,ispinor,1),s_i(im,1,ispinor,1),im=1,ndim),ispinor=1,nspinor)
             if(ioerr<0) then
!              write(std_out,*)" HDC IOERR<",ioerr
             else if(ioerr>0) then
!              write(std_out,*)" HDC IOERR>",ioerr
             endif
             do ispinor=1,nspinor
               do im=1,ndim
                 self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)&
&                 =cmplx(s_r(im,1,ispinor,1),s_i(im,1,ispinor,1),kind=dp)
               enddo
             enddo
             !write(6,*) "read selfhdc",self%hdc%matlu(1)%mat(1,1,1,1,1)
           else if(readimagonly==1.and..not.present(opt_hdc)) then
             do ispinor=1,nspinor
               do im=1,ndim
                 self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=czero
               enddo
             enddo
           else
            write(std_out,*) "     self%hdc fixed in kramerskronig_self"
           endif
           close(unitselffunc_arr(iall))
!         enddo ! ispinor
       enddo ! isppol
       ABI_DEALLOCATE(s_r)
       ABI_DEALLOCATE(s_i)
     endif ! lpawu=/-1
   enddo ! iatom
   ABI_DEALLOCATE(unitselffunc_arr)
 endif ! optrw==2.or.myproc==master
! call xmpi_barrier(spacecomm)
           !write(std_out,*) "9"
!   - For the Tentative rotation of the self-energy file (begin destroy)
 if(present(pawang)) then
   call destroy_oper(energy_level)
   call destroy_matlu(level_diag,natom)
   call destroy_matlu(selfrotmatlu,natom)
   ABI_DATATYPE_DEALLOCATE(level_diag)
   ABI_DATATYPE_DEALLOCATE(selfrotmatlu)
   do iatom=1,natom
     if(paw_dmft%lpawu(iatom)/=-1) then
       do isppol=1,nsppol
         ABI_DEALLOCATE(eigvectmatlu(iatom,isppol)%value)
       end do
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(eigvectmatlu)
   ABI_DEALLOCATE(unitselfrot) ! 7 is the max ndim possible
 endif
!   - For the Tentative rotation of the self-energy file (end destroy)

!  ===========================
!  == Error messages
!  ===========================
 if(optrw==1) then
!   call xmpi_barrier(spacecomm)
   ncount=natom*nsppol*(nspinor**2)*(self%nw+1)*(maxval(self%oper(1)%matlu(:)%lpawu)*2+1)**2
   !write(std_out,*) ncount,maxval(pawtab(:)%lpawu)*2+1
   call xmpi_bcast(iexist2,master,spacecomm ,ier)
   call xmpi_bcast(ioerr,master,spacecomm ,ier)
   if(iexist2==0.or.ioerr<0.or.ioerr>0) then
     message = "Self file does not exist or is incomplete"
     if(readimagonly==1.or.present(opt_stop)) then
      if(readimagonly==1) then
        message = "Self file does not exist or is incomplete: check the number of self data in file"
      endif
       MSG_ERROR(message)
     else
       MSG_WARNING(message)
     endif
     if(iexist2==0) then
       write(message,'(4x,2a)') "File does not exist"
       call wrtout(std_out,message,'COLL')
     endif
     if(ioerr<0) then
       write(message,'(4x,2a)') "End of file reached"
       call wrtout(std_out,message,'COLL')
     endif
     if(ioerr>0) then
       write(message,'(4x,2a)') "Error during read statement"
       call wrtout(std_out,message,'COLL')
     endif
     if(paw_dmft%dmft_solv/=4) then
       write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to double counting term"
     else if(paw_dmft%dmft_solv==4) then
       write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to dc term - shift"
       call wrtout(std_out,message,'COLL')
       write(message,'(4x,a,a,5i5,2x,e14.7)') " No self energy is given, change dmft_rslf"
       MSG_ERROR(message)
     endif
     call wrtout(std_out,message,'COLL')
     do ifreq=1,self%nw
!       write(std_out,*) "before",self%oper(1)%matlu(1)%mat(1,1,1,1,1)
!       write(std_out,*) "before",self%hdc%matlu(1)%mat(1,1,1,1,1)
       call copy_matlu(self%hdc%matlu,self%oper(ifreq)%matlu,natom)
!       write(std_out,*) "after",self%oper(1)%matlu(1)%mat(1,1,1,1,1)
!       write(std_out,*) "before",self%hdc%matlu(1)%mat(1,1,1,1,1)
       if(paw_dmft%dmft_solv==4) then
!         if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
         call shift_matlu(self%oper(ifreq)%matlu,natom,cmplx(self%qmc_shift,0.d0,kind=dp),1)
!         if(ifreq==1) write(std_out,*) "self after dc and shift",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
!         if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       endif
     enddo
   else ! test read successfull
!   call xmpi_barrier(spacecomm)
!! lignes 924-928 semblent inutiles puisque la valeur de paw_dmft%fermie creee
!! en ligne 927 est ecrasee en ligne 992. BA+jmb
!!     ABI_ALLOCATE(fermie_read2,(1))
!!     fermie_read2(1)=fermie_read
!!     call xmpi_sum(fermie_read2,spacecomm ,ier)
!!     paw_dmft%fermie=fermie_read2(1)
!!     ABI_DEALLOCATE(fermie_read2)

!  ===========================
!   bcast to other proc
!  ===========================
     ABI_ALLOCATE(buffer,(ncount))
     ABI_ALLOCATE(fermie_read2,(1))
     buffer(:)=czero
!! BA+jmb
     fermie_read2=zero
   !write(std_out,*) self%nw
     if(myproc==master) then

!               == Send read data to all process
       icount=0
       fermie_read2(1)=fermie_read
       do iatom=1,natom
         if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
           ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
           do isppol=1,nsppol
             do ispinor=1,nspinor
!               Self energy-----------
               do ifreq=1,self%nw
                 do ispinor1=1,nspinor
                   do im=1,ndim
                     do im1=1,ndim
                       icount=icount+1
                       if(icount.gt.ncount) then
                         write(message,'(2a,2i5)') ch10,"Error buffer",icount,ncount
                         iexit=1
                         MSG_ERROR(message)
                       endif
                       buffer(icount)=self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                     enddo
                   enddo
                 enddo  ! ispinor1
               enddo
!               Double counting-------
               do im=1,ndim
                 icount=icount+1
                 if(icount.gt.ncount) then
                   write(message,'(2a,2i5)') ch10,"Error buffer",icount,ncount
                   iexit=1
                   MSG_ERROR(message)
                 endif
                 buffer(icount)=self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)
               enddo
             enddo  ! ispinor
           enddo ! isppol
         endif ! lpawu=/-1
       enddo ! iatom
     endif
!    call xmpi_bcast(buffer,master,spacecomm ,ier)
!    call xmpi_sum(iexit,spacecomm ,ier)
!!JB call xmpi_barrier(spacecomm)
     call xmpi_sum(buffer,spacecomm ,ier)
!!JB call xmpi_barrier(spacecomm)

! bcast fermi level
   call xmpi_sum(fermie_read2,spacecomm ,ier)

     if(ier/=0) then
       message =  "error in xmpi_sum in rw_self"
       MSG_ERROR(message)
     endif
     paw_dmft%fermie=fermie_read2(1)
!     write(std_out,*) "Fermi level",paw_dmft%fermie
     icount=0
     do iatom=1,natom
       if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
         ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
         do isppol=1,nsppol
           do ispinor=1,nspinor
!             self ---------------
             do ifreq=1,self%nw
               do ispinor1=1,nspinor
                 do im=1,ndim
                   do im1=1,ndim
                     icount=icount+1
                     self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=buffer(icount)
                     !write(6,*)'self procs', ifreq, self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                   enddo
                 enddo
               enddo
             enddo
!             hdc  ---------------
             do im=1,ndim
               icount=icount+1
               self%hdc%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)=buffer(icount)
             enddo
           enddo
         enddo ! isppol
       endif ! lpawu=/-1
     enddo ! iatom
     ABI_DEALLOCATE(fermie_read2)
     ABI_DEALLOCATE(buffer)
   endif  ! test read successful
 endif  ! optrw==1
!   call flush_unit(std_out)
!   MSG_ERROR("Aboring now")
 if(optrw==0) then
   if(paw_dmft%dmft_rslf==0) then
     if(paw_dmft%dmft_solv/=4) then
       write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to double counting term"
     else if(paw_dmft%dmft_solv==4) then
       write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to dc term - shift"
     endif
   else if (paw_dmft%dmft_rslf==-1) then
     write(message,'(4x,a,a,5i5,2x,e14.7)') "-> Put Self-Energy Equal to zero"
   endif
   call wrtout(std_out,message,'COLL')
   do ifreq=1,self%nw
     if(paw_dmft%dmft_rslf==0) then
       call copy_matlu(self%hdc%matlu,self%oper(ifreq)%matlu,natom)
      ! if(ifreq==1) write(std_out,*) "self after dc",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
       if(paw_dmft%dmft_solv==4) then
      !   if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
         call shift_matlu(self%oper(ifreq)%matlu,natom,cmplx(self%qmc_shift,0.d0,kind=dp),1)
      !   if(ifreq==1) write(std_out,*) "self after dc and shift",self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)
      !   if(ifreq==1) write(std_out,*) "shift",self%qmc_shift(1)
       endif
     else if (paw_dmft%dmft_rslf==-1) then
       do iatom=1,natom
         if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
           ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
           do isppol=1,nsppol
             do ispinor=1,nspinor
               do im=1,ndim
                 self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)= czero
               enddo
             enddo
           enddo
         endif
       enddo
     endif
   enddo
 endif

! write(std_out,*) "optrw,use_fixed_self,istep,iter,istep_imp,iter_imp"
! write(std_out,*) optrw,paw_dmft%use_fixed_self,istep,iter,istep_imp,iter_imp
 if((optrw==1.or.optrw==3).and.paw_dmft%use_fixed_self>0.and.istep<=istep_imp.and.iter<=iter_imp) then
   write(message,'(4x,a)') "-> Put Self-Energy Equal to imposed self-energy"
   call wrtout(std_out,message,'COLL')
   iatu=0
   do iatom=1,natom
     if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
       iatu=iatu+1
       do ifreq=1,self%nw
         ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
         do isppol=1,nsppol
           do ispinor=1,nspinor
             do im=1,ndim
               if(nspinor==1) then
                 self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)= paw_dmft%fixed_self(im,im,isppol,iatu)
!            write(std_out,*) paw_dmft%fixed_self(im,im,isppol,iatu)
               else
                 self%oper(ifreq)%matlu(iatom)%mat(im,im,isppol,ispinor,ispinor)= paw_dmft%fixed_self(im,im,ispinor,iatu)
!                 write(message,'(a,i4,i4,2x,e20.10)') " Fixed self not implemented for nspinor==2"
!                 call wrtout(std_out,  message,'COLL')
!                 MSG_ERROR("Aboring now")
               endif
             enddo
           enddo
         enddo
       enddo
     endif
   enddo

 endif




end subroutine rw_self
!!***

!!****f* m_self/new_self
!! NAME
!! new_self
!!
!! FUNCTION
!!
!!  Mix Old and New self_energy with the mixing coefficient dmft_mxsf
!!
!! INPUTS
!!  self <type(self_type)>= variables related to self-energy
!!  self_new <type(self_type)>= variables related to the new self-energy
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  opt_mix not used
!!
!! OUTPUT
!!  self <type(self_type)>= variables related to mixed self-energy
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      shift_matlu,wrtout
!!
!! SOURCE

subroutine new_self(self,self_new,paw_dmft,opt_mix)

!Arguments ------------------------------------
!type
! type(crystal_t),intent(in) :: cryst_struc
 type(self_type),intent(inout) :: self
 type(self_type),intent(in) :: self_new
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer,intent(in) :: opt_mix

!local variables-------------------------------
 integer :: iatom,icount,ifreq,im,im1,ispinor,isppol,ispinor1
 integer :: natom,ndim,nspinor,nsppol
 real(dp) :: alpha,diff_self,sum_self
 complex(dpc) :: s1,s2
 character(len=500) :: message
! *********************************************************************
 natom=self%hdc%natom
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 alpha=paw_dmft%dmft_mxsf

 if(opt_mix==1) then
 endif
 sum_self=zero
 diff_self=zero
 icount=0
 do iatom=1,natom
   if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
     ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           do im=1,ndim
             do im1=1,ndim
               do ifreq=1,self%nw
!  warning: self_new is the recent self-energy, which is mixed with self
!  to give self= mixed self energy. self_new is deallocated just after.
                 self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=     &
&                  (one-alpha)*self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1) +    &
&                  (alpha)*self_new%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                     s1=self%hdc%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                 s2=self_new%hdc%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
                 if((ispinor==ispinor1).and.(im==im1)) then
                   diff_self=diff_self+dsqrt(real(s1-s2)**2+aimag(s1-s2)**2)
                   sum_self=sum_self+dsqrt(real(s1)**2+aimag(s1)**2)
                   icount=icount+1
                 endif
               enddo
               self%hdc%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)=             &
&               (one-alpha)*self%hdc%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)   +          &
&               (alpha)*self_new%hdc%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)
             enddo
           enddo
         enddo
       enddo
     enddo ! isppol
   endif ! lpawu=/-1
 enddo ! iatom
 diff_self=diff_self/icount

 write(message,'(8x,a,e12.5)') "DMFT Loop: Precision on self-energy is",diff_self
 call wrtout(std_out,message,'COLL')
 if(diff_self<paw_dmft%dmft_fermi_prec.and.sum_self>tol6.and.paw_dmft%idmftloop>=2) then
    write(message,'(a,8x,a,e9.2,a,8x,a)') ch10, "Change of self =<", paw_dmft%dmft_fermi_prec,&
&    ch10,"DMFT Loop: Self Energy is converged"
    call wrtout(std_out,message,'COLL')
    self%iself_cv=1
 else
    write(message,'(a,8x,a)') ch10,"DMFT Loop: Self Energy is not converged"
    call wrtout(std_out,message,'COLL')
    self%iself_cv=0
 endif


end subroutine new_self
!!***

!!****f* m_self/make_qmcshift_self
!! NAME
!! make_qmcshift_hu
!!
!! FUNCTION
!!
!! INPUTS
!!  hu <type(hu_type)> = U interaction
!!  paw_dmft  <type(paw_dmft_type)> = paw+dmft related data
!!
!! OUTPUT
!!  self%qmc_shift in self <type(self_type)> = Self-energy
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine make_qmcshift_self(cryst_struc,hu,self,apply)

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(hu_type),intent(in) :: hu(cryst_struc%ntypat)
 type(self_type),intent(inout) :: self
 logical, optional :: apply

!Local variables-------------------------------
 integer :: im,iatom,ifreq,itypat,lpawu,tndim
 real(dp) :: hu_shift2
 character(len=500) :: message
! *********************************************************************

 do iatom = 1 , cryst_struc%natom
   lpawu=self%hdc%matlu(iatom)%lpawu
   tndim=2*lpawu+1
   itypat=cryst_struc%typat(iatom)
   if(lpawu/=-1) then
     self%qmc_shift(iatom) = zero
     do im =1, 2*tndim-1
!       write(std_out,*)"make before",self%qmc_shift(iatom)
       self%qmc_shift(iatom) = self%qmc_shift(iatom) + hu(itypat)%uqmc(im)
!       write(std_out,*)"make after",self%qmc_shift(iatom)
     enddo
     self%qmc_shift(iatom) = self%qmc_shift(iatom) / two
     hu_shift2 = hu(itypat)%uqmc(1)

     do im = 2*tndim, 2*tndim + 2*tndim -3
       hu_shift2 = hu_shift2 + hu(itypat)%uqmc(im)
     enddo

     hu_shift2 = hu_shift2 / two
     write(message,'(2a,i4)')  ch10,'  -------> For Correlated atom',iatom
     call wrtout(std_out,  message,'COLL')

     if(abs(self%qmc_shift(iatom)-hu_shift2)>tol6) then
       write(message,'(2a,2f16.7)')  "  Shift for QMC is not correctly"&
&      ," computed",self%qmc_shift(iatom),hu_shift2
       MSG_ERROR(message)
     endif ! shifts not equals

     write(message,'(4x,a,f16.7)')  &
&     "  Shift for QMC (used to compute G(w)) is (in Ha) :",&
&     self%qmc_shift(iatom)
     call wrtout(std_out,  message,'COLL')

     self%qmc_xmu(iatom)=-self%qmc_shift(iatom)
     self%qmc_xmu(iatom)=zero
     write(message,'(4x,a,f16.7)')  &
&     "Artificial Shift used in QMC AND to compute G is (in Ha) :",self%qmc_xmu(iatom)
     MSG_WARNING(message)

   endif ! lpawu/=1
 enddo ! natom

 if(present(apply)) then
   if (apply) then
     write(message,'(5x,a,f16.7,a)')  " Shifts applied to self"
     call wrtout(std_out,  message,'COLL')
     do ifreq=1,self%nw
       call shift_matlu(self%oper(ifreq)%matlu,cryst_struc%natom,cmplx(self%qmc_shift,0.d0,kind=dp),1)
       call shift_matlu(self%oper(ifreq)%matlu,cryst_struc%natom,cmplx(self%qmc_xmu,0.d0,kind=dp),1)
     enddo
   endif
 endif


end subroutine make_qmcshift_self
!!***

!!****f* m_self/kramerskronig_self
!! NAME
!! kramerskronig_self
!!
!! FUNCTION
!!
!! INPUTS
!!  hu <type(hu_type)> = U interaction
!!  paw_dmft  <type(paw_dmft_type)> = paw+dmft related data
!!
!! OUTPUT
!!  self%qmc_shift in self <type(self_type)> = Self-energy
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine kramerskronig_self(self,selflimit,selfhdc,filapp)

!Arguments ------------------------------------
!type
 type(self_type),intent(inout) :: self
 type(matlu_type),intent(in) :: selflimit(self%hdc%natom)
 type(matlu_type),intent(in) :: selfhdc(self%hdc%natom)
 character(len=fnlen), intent(in) :: filapp

!Local variables-------------------------------
 integer :: ifreq,jfreq,isppol,ispinor,ispinor1,im,im1,iatom
 real(dp), allocatable :: selftemp_re(:)
 real(dp), allocatable :: selftemp_imag(:)
 integer :: natom,ndim,nsppol,nspinor
 real(dp) :: delta
 character(len=500) :: message
! *********************************************************************
 delta=0.0000000
 ABI_ALLOCATE(selftemp_re,(self%nw))
 ABI_ALLOCATE(selftemp_imag,(self%nw))
 natom=self%hdc%natom
 nsppol  = self%hdc%nsppol
 nspinor=self%hdc%nspinor
 write(message,'(2a,i4)')  ch10,'  ------ Limit of real part of Self'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(selflimit,natom,3)

!print norms
 write(message,'(2a,i4)')  ch10,'  ------ Double counting'
 call wrtout(std_out,  message,'COLL')

 call print_matlu(selfhdc,natom,3)
 open(unit=67,file=trim(filapp)//"_DFTDMFT_Self_realaxis_from_maxent_and_kramerskronig.dat", status='unknown',form='formatted')
 rewind(67)
!  Compute limit of Real Part and put in double counting energy.
! call copy_matlu(selfhdc,self%hdc%matlu,natom)
     !!write(6,*) "selfhdc   kramerskronig",selfhdc(1)%mat(1,1,1,1,1)
     !write(6,*) "selfr%hdc kramerskronig",self%hdc%matlu(1)%mat(1,1,1,1,1)
 do iatom=1,natom
   if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
     ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           do im=1,ndim
             do im1=1,ndim
               !write(6,*) "realpart",real(selflimit(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               do ifreq=1,self%nw
                 selftemp_re(ifreq)=zero
                 selftemp_imag(ifreq)=aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
                 do jfreq=1,self%nw-1
                   if(jfreq==ifreq) cycle
!                   selftemp_re(ifreq)=selftemp_re(ifreq) -   &
! &                   aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
! &                   *(self%omega(ifreq)-self%omega(jfreq)) &
! &                   /((self%omega(ifreq)-self%omega(jfreq))**2+delta**2)&
! &                   *(self%omega(jfreq+1)-self%omega(jfreq))
                   selftemp_re(ifreq)=selftemp_re(ifreq) -   &
 &                   aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
 &                   /(self%omega(ifreq)-self%omega(jfreq)) * (self%omega(jfreq+1)-self%omega(jfreq))
                 enddo
                 selftemp_re(ifreq)=selftemp_re(ifreq)/pi
                 !write(671,*)  self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
               enddo
!                 TEST*************************
!               do ifreq=1,self%nw
!                 selftemp_imag(ifreq)=zero
!                 do jfreq=1,self%nw-1
!                   if(jfreq==ifreq) cycle
!!                   selftemp_re(ifreq)=selftemp_re(ifreq) -   &
!! &                   aimag(self%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
!! &                   *(self%omega(ifreq)-self%omega(jfreq)) &
!! &                   /((self%omega(ifreq)-self%omega(jfreq))**2+delta**2)&
!! &                   *(self%omega(jfreq+1)-self%omega(jfreq))
!                   selftemp_imag(ifreq)=selftemp_imag(ifreq) +    &
! &                   selftemp_re(jfreq)  &
! &                   /(self%omega(ifreq)-self%omega(jfreq)) * (self%omega(jfreq+1)-self%omega(jfreq))
!                 enddo
!                 selftemp_imag(ifreq)=selftemp_imag(ifreq)/pi
!                 write(672,*)  self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
!               enddo
!                 TEST*************************
              ! write(6,*) "TWO FACTOR IS PUT BECAUSE OF MAXENT CODE ??"
               do ifreq=1,self%nw
!                 write(68,*)  self%omega(ifreq),selftemp_re(ifreq),selftemp_imag(ifreq)
                 selftemp_re(ifreq)=selftemp_re(ifreq)+ &
 &                 real(selflimit(iatom)%mat(im,im1,isppol,ispinor,ispinor1)- &
 &                 selfhdc(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
                 self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
  &                       =cmplx(selftemp_re(ifreq),selftemp_imag(ifreq),kind=dp)/two
  !&                       =cmplx(0.d0,selftemp_imag(ifreq),kind=dp)/two
!  &                       =cmplx(selftemp_re(ifreq),0.d0,kind=dp)/two
  !               self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
  !&                       =cmplx(selftemp_re(ifreq),0.d0,kind=dp)/two
  !&                       =cmplx(0.d0,0.d0,kind=dp)/two
!    The factor two is here to compensate for the factor two in OmegaMaxent..
!  &                       =cmplx(selftemp_re(ifreq),0.0,kind=dp)
                 write(67,*)  self%omega(ifreq),real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))&
                 ,aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               enddo
               write(67,*)
               !write(68,*)
               !!!!!!!!!! Z renormalization
!               i0=389
!               slope=(selftemp_re(i0+1)-selftemp_re(i0))/&
!                     (self%omega(i0+1)-self%omega(i0))
!               y0= selftemp_re(i0)
!               do ifreq=1,self%nw
!                 selftemp_re(ifreq)=slope * (self%omega(ifreq)-self%omega(i0)) + y0
!                 selftemp_imag(ifreq)=zero
!                 self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)&
!  &                       =cmplx(selftemp_re(ifreq),selftemp_imag(ifreq),kind=dp)/two
!                 write(6777,*)  self%omega(ifreq),real(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1)),aimag(self%oper(ifreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
!               enddo
               !!!!!!!!!!
             enddo
           enddo
         enddo
       enddo
     enddo ! isppol
   endif ! lpawu=/-1
 enddo ! iatom
 close(67)
     !write(6,*) "self1",aimag(self%oper(489)%matlu(1)%mat(1,1,1,1,1))

 ABI_DEALLOCATE(selftemp_re)
 ABI_DEALLOCATE(selftemp_imag)


end subroutine kramerskronig_self
!!***

!!****f* m_self/selfreal2imag_self
!! NAME
!! selfreal2imag_self
!!
!! FUNCTION
!!
!! INPUTS
!!  hu <type(hu_type)> = U interaction
!!  paw_dmft  <type(paw_dmft_type)> = paw+dmft related data
!!
!! OUTPUT
!!  self%qmc_shift in self <type(self_type)> = Self-energy
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine selfreal2imag_self(selfr,self,filapp)

!Arguments ------------------------------------
!type
 type(self_type),intent(inout) :: selfr
 type(self_type),intent(inout) :: self
 character(len=fnlen), intent(in) :: filapp

!Local variables-------------------------------
 integer :: ifreq,jfreq,isppol,ispinor,ispinor1,im,im1,iatom
 complex(dpc), allocatable :: selftempmatsub(:)
 integer :: natom,ndim,nsppol,nspinor
 real(dp) :: delta
! *********************************************************************
 delta=0.0000000
 ABI_ALLOCATE(selftempmatsub,(self%nw))
 natom=self%hdc%natom
 nsppol  = self%hdc%nsppol
 nspinor=self%hdc%nspinor
!  Compute limit of Real Part and put in double counting energy.
! call copy_matlu(selfhdc,self%hdc%matlu,natom)
     !write(6,*) "self3",aimag(selfr%oper(489)%matlu(1)%mat(1,1,1,1,1))

 open(unit=672,file=trim(filapp)//"_DFTDMFT_Self_forcheck_imagaxis_from_realaxis.dat", status='unknown',form='formatted')
 do iatom=1,natom
   if(self%oper(1)%matlu(iatom)%lpawu.ne.-1) then
     ndim=2*self%oper(1)%matlu(iatom)%lpawu+1
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do ispinor1=1,nspinor
           do im=1,ndim
             do im1=1,ndim
               do jfreq=1,selfr%nw-1
               !  write(6700,*)  selfr%omega(jfreq),aimag(selfr%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))
               enddo
               !  write(6700,*)
               do ifreq=1,self%nw
                 selftempmatsub(ifreq)=czero
                 do jfreq=1,selfr%nw-1
                   selftempmatsub(ifreq)=selftempmatsub(ifreq) -   &
 &                   aimag(selfr%oper(jfreq)%matlu(iatom)%mat(im,im1,isppol,ispinor,ispinor1))  &
 &                   /(cmplx(zero,self%omega(ifreq),kind=dp)-selfr%omega(jfreq))   &
 &                 * (selfr%omega(jfreq+1)-selfr%omega(jfreq))
                 enddo
                 selftempmatsub(ifreq)=selftempmatsub(ifreq)/pi
                 write(672,*)  self%omega(ifreq),real(selftempmatsub(ifreq)),aimag(selftempmatsub(ifreq))
               enddo
                 write(672,*)
             enddo
           enddo
         enddo
       enddo
     enddo ! isppol
   endif ! lpawu=/-1
 enddo ! iatom
 close(672)

 ABI_DEALLOCATE(selftempmatsub)


end subroutine selfreal2imag_self
!!***

END MODULE m_self
!!***
