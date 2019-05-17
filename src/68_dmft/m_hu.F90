!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hu
!! NAME
!!  m_hu
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_hu

 use m_abicore

 use defs_basis
 use m_pawtab, only : pawtab_type

 implicit none

 private

 public :: init_hu
 public :: destroy_hu
! public :: qmc_hu
 public :: print_hu
 public :: vee2udens_hu
 public :: rotatevee_hu
 public :: printvee_hu
 public :: vee2udensatom_hu
 public :: vee_slm2ylm_hu
 public :: vee_ndim2tndim_hu
 public :: vee_ndim2tndim_hu_r
 public :: udens_slatercondon_hu
 public :: udens_inglis_hu


!!***

!!****t* m_hu/hu_type
!! NAME
!!  hu_type
!!
!! FUNCTION
!!  This structured datatype contains interaction matrices for the correlated subspace
!!
!! SOURCE

 type, public :: hu_type ! for each typat

  integer :: lpawu

  logical :: jmjbasis

  real(dp) :: upawu    ! => upaw

  real(dp) :: jpawu    ! => jpaw

  real(dp) :: f2_sla    ! => f2_sla

  real(dp) :: f4of2_sla    ! => f4of2_sla

  real(dp) :: f6of2_sla    ! => f6of2_sla

  logical :: jpawu_zero  ! true if all jpawu are zero
                         ! false if one of the jpaw is not zero

  real(dp), pointer :: vee(:,:,:,:) => null() ! => vee

  real(dp), allocatable :: uqmc(:)

  real(dp), allocatable :: udens(:,:)

 end type hu_type

!----------------------------------------------------------------------


CONTAINS  !========================================================================================
!!***

!!****f* m_hu/init_hu
!! NAME
!! init_hu
!!
!! FUNCTION
!!  Allocate variables used in type hu_type.
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  pawtab <type(pawtab)>=paw related data
!!
!! OUTPUTS
!!  hu <type(hu_type)>= U interaction
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_hu(cryst_struc,pawtab,hu,t2g,x2my2d)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 type(pawtab_type), target, intent(in)  :: pawtab(cryst_struc%ntypat)
 type(hu_type), intent(inout) :: hu(cryst_struc%ntypat)
 integer, intent(in) :: t2g,x2my2d
!Local variables ------------------------------------
 integer :: itypat,i,ij,ij1,ij2,j,lpawu,ms,ms1,m,m1,ndim
 integer :: ns,ns1,n,n1
 integer, allocatable :: xij(:,:)
 real(dp) :: xtemp
 character(len=500) :: message
!************************************************************************
 write(message,'(2a)') ch10,"  == Compute Interactions for DMFT"
 call wrtout(std_out,message,'COLL')

 xtemp=zero

! ====================================
!  Compute hu(iatom)%uqmc from vee
! ====================================
 hu(1)%jpawu_zero=.true.
 do itypat=1,cryst_struc%ntypat
   hu(itypat)%lpawu=pawtab(itypat)%lpawu
   hu(itypat)%jmjbasis=.false.
   if(t2g==1.and.hu(itypat)%lpawu==2) hu(itypat)%lpawu=1
   if(x2my2d==1.and.hu(itypat)%lpawu==2) hu(itypat)%lpawu=0
   lpawu=hu(itypat)%lpawu
   if(lpawu.ne.-1) then
     hu(itypat)%upawu=pawtab(itypat)%upawu
     hu(itypat)%jpawu=pawtab(itypat)%jpawu

   ! The if below are necessary for an unknown reason with the NAG  compiler.
    if(pawtab(itypat)%f4of2_sla>=0) then
      hu(itypat)%f4of2_sla = pawtab(itypat)%f4of2_sla
    else if(pawtab(itypat)%f4of2_sla==0) then
      hu(itypat)%f4of2_sla = 0_dp
    else if(pawtab(itypat)%f4of2_sla<0) then
      hu(itypat)%f4of2_sla = pawtab(itypat)%f4of2_sla
    endif
    if(pawtab(itypat)%f6of2_sla>=0) then
      hu(itypat)%f6of2_sla = pawtab(itypat)%f6of2_sla
    else if(pawtab(itypat)%f6of2_sla==0) then
      hu(itypat)%f6of2_sla = 0_dp
    else if(pawtab(itypat)%f6of2_sla<0) then
      hu(itypat)%f6of2_sla = pawtab(itypat)%f6of2_sla
    endif

!    This is copied from pawpuxinit: it would be better not to duplicate
!    these lines.
     if(lpawu==0) then
       hu(itypat)%f2_sla=zero
     else if(lpawu==1) then
       hu(itypat)%f2_sla=pawtab(itypat)%jpawu*5._dp
     else if(lpawu==2) then
       hu(itypat)%f2_sla=pawtab(itypat)%jpawu*14._dp/(One+hu(itypat)%f4of2_sla)
     else if(lpawu==3) then
       hu(itypat)%f2_sla=pawtab(itypat)%jpawu*6435._dp/(286._dp+&
&       195._dp*hu(itypat)%f4of2_sla+250._dp*hu(itypat)%f6of2_sla)
     endif

     if(hu(itypat)%jpawu>tol4) hu(1)%jpawu_zero=.false.
     ndim=2*lpawu+1
!     ndim1=2*hu(itypat)%lpawu+1
     write(message,'(2a,i4)')  ch10,'  -------> For Correlated Species', itypat
     call wrtout(std_out,  message,'COLL')
!     allocate(hu(itypat)%vee(ndim,ndim,ndim,ndim))
     ABI_ALLOCATE(hu(itypat)%uqmc,(ndim*(2*ndim-1)))
     ABI_ALLOCATE(hu(itypat)%udens,(2*ndim,2*ndim))
     ABI_ALLOCATE(xij,(2*ndim,2*ndim))
     if(t2g==0.and.x2my2d==0) then
       hu(itypat)%vee => pawtab(itypat)%vee
!   t2g case begin
     else if(t2g==1.and.hu(itypat)%lpawu==1) then
       ABI_ALLOCATE(hu(itypat)%vee,(ndim,ndim,ndim,ndim))
       n=0
       do m=1,5
         if((m/=1.and.m/=2.and.m/=4)) cycle
         n=n+1
         ns=0
         do ms=1,5
           if((ms/=1.and.ms/=2.and.ms/=4)) cycle
           ns=ns+1
           n1=0
           do m1=1,5
             if((m1/=1.and.m1/=2.and.m1/=4)) cycle
             n1=n1+1
             ns1=0
             do ms1=1,5
               if((ms1/=1.and.ms1/=2.and.ms1/=4)) cycle
               ns1=ns1+1
               hu(itypat)%vee(n,ns,n1,ns1)=pawtab(itypat)%vee(m,ms,m1,ms1)
             enddo
           enddo
         enddo
       enddo
!   t2g case end
!   x2my2d case begin
     else if(x2my2d==1.and.hu(itypat)%lpawu==0) then
       ABI_ALLOCATE(hu(itypat)%vee,(ndim,ndim,ndim,ndim))
       hu(itypat)%vee(1,1,1,1)=pawtab(itypat)%upawu
     endif
!   x2my2d case end

     hu(itypat)%udens=zero
     ij=0
     do ms=1,2*ndim-1
         xij(ms,ms)=0
       do ms1=ms+1,2*ndim
         ij=ij+1
         xij(ms,ms1)=ij
         xij(ms1,ms)=ij
         if(ms<=ndim.and.ms1>ndim) then
           m1 = ms1 - ndim
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else if(ms<=ndim.and.ms1<=ndim) then
           m1 = ms1
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else
           m1 = ms1 - ndim
           m  = ms  - ndim
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         endif
       enddo
     enddo
     xij(2*ndim,2*ndim)=0
     write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation "
     call wrtout(std_out,  message,'COLL')
     write(message,'(1x,14(2x,i5))') (m,m=1,2*ndim)
     call wrtout(std_out,  message,'COLL')
!     xtemp1b=0.d0
! ====================================
!  Print hu(iatom)%uqmc
! ====================================
     ij1=-10
     ij2=-10
     ij=0
     do i=1,2*ndim
       do j=i+1,2*ndim
         ij=ij+1
         if(j==i+1) ij1=ij
         if(j==2*ndim) ij2=ij
       enddo
!       write(std_out,*) itypat
!       do m=1,i
!        write(std_out,*) i,m
!        write(std_out,*) xij(i,m)
!        write(std_out,*) ij1,ij2
!       enddo
       if(i==1)               write(message,'(i3,14f7.3)') &
&                              i,xtemp, (hu(itypat)%uqmc(m),m=ij1,ij2)
       if(i/=2*ndim.and.i/=1) write(message,'(i3,14f7.3)') i, &
&        (hu(itypat)%uqmc(xij(i,m)), m=1,i-1),xtemp, (hu(itypat)%uqmc(m),m=ij1,ij2)
       if(i==2*ndim)          write(message,'(i3,14f7.3)') i, &
&                  (hu(itypat)%uqmc(xij(i,m)), m=1,i-1),xtemp
       call wrtout(std_out,  message,'COLL')
     enddo
       write(message,'(5x,a)') "--------------------------------------------------------"
       call wrtout(std_out,  message,'COLL')
     ABI_DEALLOCATE(xij)
   else
     hu(itypat)%upawu=zero
     hu(itypat)%jpawu=zero
!     allocate(hu(itypat)%vee(0,0,0,0))
   endif
 enddo ! itypat

end subroutine init_hu
!!***

!!****f* m_hu/destroy_hu
!! NAME
!! destroy_hu
!!
!! FUNCTION
!!  deallocate hu
!!
!! INPUTS
!!  ntypat = number of species
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_hu(hu,ntypat,t2g,x2my2d)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ntypat
 type(hu_type),intent(inout) :: hu(ntypat)
 integer, intent(in) :: t2g,x2my2d

!Local variables-------------------------------
 integer :: itypat

! *********************************************************************

 do itypat=1,ntypat
!  if ( allocated(hu(itypat)%vee) )  deallocate(hu(itypat)%vee)
  if ( allocated(hu(itypat)%uqmc) )   then
    ABI_DEALLOCATE(hu(itypat)%uqmc)
  end if
  if ( allocated(hu(itypat)%udens) )   then
    ABI_DEALLOCATE(hu(itypat)%udens)
  end if
  if(t2g==1.and.hu(itypat)%lpawu==1) then
    ABI_DEALLOCATE(hu(itypat)%vee)
  else if(x2my2d==1.and.hu(itypat)%lpawu==0) then
    ABI_DEALLOCATE(hu(itypat)%vee)
  else
    hu(itypat)%vee => null()
  endif
 enddo

end subroutine destroy_hu
!!***

!!****f* m_hu/print_hu
!! NAME
!! print_hu
!!
!! FUNCTION
!!  print density density interaction (used for DFT+DMFT)
!!
!! INPUTS
!!  ntypat = number of species
!!  prtopt = option for printing
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_hu(hu,ntypat,prtopt)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: ntypat
 type(hu_type),intent(in) :: hu(ntypat)
 integer :: prtopt

!Local variables-------------------------------
 integer :: itypat
 integer :: lpawu,ms,ms1,m,ndim
 character(len=500) :: message
! *********************************************************************

 do itypat = 1 , ntypat
   lpawu=hu(itypat)%lpawu
   if(lpawu/=-1) then
     ndim=2*lpawu+1
     write(message,'(2a,i4)')  ch10,'  -------> For Correlated species'
     call wrtout(std_out,  message,'COLL')
     if(prtopt==0) then
       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation in Slm basis "
     else if(prtopt==1) then
       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation in diagonal basis"
     else if(prtopt==2) then
       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation in Ylm basis"
     else if(prtopt==3) then
       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation in JMJ basis"
     endif
     call wrtout(std_out,  message,'COLL')
     write(message,'(1x,14(2x,i5))') (m,m=1,2*ndim)
     call wrtout(std_out,  message,'COLL')
       do ms=1,2*ndim
          write(message,'(i3,14f7.3)') &
&          ms, (hu(itypat)%udens(ms,ms1),ms1=1,2*ndim)
          call wrtout(std_out,  message,'COLL')
       enddo
       write(message,'(5x,a)') "--------------------------------------------------------"
       call wrtout(std_out,  message,'COLL')
   endif ! lpawu/=1
 enddo ! ntypat


end subroutine print_hu
!!***

!!****f* m_hu/vee2udens_hu
!! NAME
!! print_hu
!!
!! FUNCTION
!!  interaction udens in recomputed from new vee.
!!
!! INPUTS
!!  ntypat = number of species
!!  prtopt = option for printing
!!
!! OUTPUT
!!
!! SIDE EFFECT
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee2udens_hu(hu,ntypat,prtopt)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: ntypat
 type(hu_type),intent(inout) :: hu(ntypat)
 integer :: prtopt

!Local variables-------------------------------
 integer :: ij,itypat
 integer :: lpawu,m1,ms,ms1,m,ndim
 character(len=500) :: message
! *********************************************************************
 do itypat=1,ntypat
   lpawu=hu(itypat)%lpawu
   if(lpawu.ne.-1) then
     ndim=2*lpawu+1
     write(message,'(2a,i4)')  ch10,'  -------> For Correlated Species', itypat
     call wrtout(std_out,  message,'COLL')

     hu(itypat)%udens=zero
     ij=0
     do ms=1,2*ndim-1
!         xij(ms,ms)=0
       do ms1=ms+1,2*ndim
         ij=ij+1
!         xij(ms,ms1)=ij
!         xij(ms1,ms)=ij
         if(ms<=ndim.and.ms1>ndim) then
           m1 = ms1 - ndim
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%vee(m,m1,m,m1)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else if(ms<=ndim.and.ms1<=ndim) then
           m1 = ms1
           m  = ms
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         else
           m1 = ms1 - ndim
           m  = ms  - ndim
           hu(itypat)%uqmc(ij)=hu(itypat)%vee(m,m1,m,m1)-hu(itypat)%vee(m,m1,m1,m)
           hu(itypat)%udens(ms,ms1)= hu(itypat)%uqmc(ij)
           hu(itypat)%udens(ms1,ms)= hu(itypat)%udens(ms,ms1)
         endif
       enddo
     enddo
!     xij(2*ndim,2*ndim)=0
!     write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation "
!     call wrtout(std_out,  message,'COLL')
!     write(message,'(1x,14(2x,i5))') (m,m=1,2*ndim)
!     call wrtout(std_out,  message,'COLL')
   endif
 enddo ! itypat
 call print_hu(hu,ntypat,prtopt)


end subroutine vee2udens_hu
!!***

!!****f* m_hu/rotatevee_hu
!! NAME
!! rotatevee_hu
!!
!! FUNCTION
!!  interaction udens in recomputed from new vee.
!!
!! INPUTS
!!  ntypat = number of species
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!
!! OUTPUT
!!
!! SIDE EFFECT
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! PARENTS
!!      hubbard_one,qmc_prep_ctqmc
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine rotatevee_hu(cryst_struc,hu,nspinor,nsppol,pawprtvol,rot_mat,udens_atoms,rot_type)

 use defs_basis
 use m_errors

 use m_crystal,          only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 type(crystal_t),intent(in) :: cryst_struc
 integer, intent(in):: nsppol,nspinor,pawprtvol,rot_type
 type(hu_type),intent(in) :: hu(cryst_struc%ntypat)
 type(coeff2c_type),optional,intent(in) :: rot_mat(cryst_struc%natom,nsppol)
 type(coeff2_type),intent(inout) :: udens_atoms(cryst_struc%natom)

!Local variables-------------------------------
 integer :: dim_vee,iatom,itypat
 integer :: lpawu,m1,m2,m3,m4,mi,mj,mk,ml,natom,ndim,tndim
 integer :: ms,ms1,nat_correl
 character(len=500) :: message
 character(len=30) :: basis_vee
 real(dp) :: xsum,xsum2,xsumnew,xsum2new,uaver
 complex(dpc),allocatable :: temp_mat(:,:)
 complex(dpc),allocatable :: temp_mat2(:,:)
 real(dp),allocatable :: veetemp(:,:,:,:),fk(:)
 complex(dpc),allocatable :: veeylm(:,:,:,:),veeslm(:,:,:,:)
 complex(dpc),allocatable :: veetmp(:,:,:),veetmp_slm(:,:,:),veetmp_ylm(:,:,:)
 complex(dpc),allocatable :: veejmj(:,:,:,:),veerotated(:,:,:,:),veenew(:,:,:,:)
 complex(dpc),allocatable :: veeylm2(:,:,:,:),veeslm2(:,:,:,:)
! *********************************************************************
 natom=cryst_struc%natom

!================================================
!  NSPINOR = 2 and J=0
!================================================
 if(hu(1)%jpawu_zero.and.nspinor==2) then
! if(3==4) then

!   call vee2udens_hu(hu,cryst_struc%ntypat,2)
   do iatom=1,natom
     itypat=cryst_struc%typat(iatom)
     lpawu=hu(itypat)%lpawu
     if(lpawu.ne.-1) then
       ndim=2*lpawu+1
       write(message,'(2a,i4)')  ch10,'  -------> For Correlated Species', itypat
       call wrtout(std_out,  message,'COLL')
    !   write(std_out,*)"ndim",ndim
    !   write(std_out,'(a,20e10.4)')"udens before vee2udensaomt", udens_atoms(1)%value
!       write(std_out,*)"vee",hu(cryst_struc%typat(iatom))%vee
       call vee2udensatom_hu(ndim,nspinor,udens_atoms(iatom)%value,hu(cryst_struc%typat(iatom))%vee,"Slm")
    !   write(std_out,*)"udens after vee2udensatom", udens_atoms(1)%value
     endif
   enddo
  ! write(std_out,*)"udensafter after", udens_atoms(1)%value


!================================================
!  NSPINOR = 2 and J/=0
!================================================
 else if (nspinor==2) then
   write(std_out,*)
   write(std_out,*)" == Rotate interaction in the level basis "


   !rot_type=0 ! keep original Slm basis
   !rot_type=2 ! rotation to the Ylm basis ! for tests
   !rot_type=4 ! use a rotation from Ylm basis to a new basis
   !rot_type=3 ! rotation to the JmJ Basis ! for tests
   !rot_type=1 ! use a rotation for diago of dmat,green, levels..

   do iatom=1,natom
     itypat=cryst_struc%typat(iatom)
     lpawu=hu(itypat)%lpawu
     if(lpawu.ne.-1) then
       ndim=2*lpawu+1
       if(pawprtvol>=3) then
!         write(message,'(2a)')  ch10," VEE INPUT AVANT TRANSFORMATION"
!         call wrtout(std_out,  message,'COLL')
!         call printvee_hu(ndim,hu(itypat)%vee,1,'Slm')
       endif
       write(message,'(2a,i4)')  ch10,'  -------> For Correlated atom', iatom
       call wrtout(std_out,  message,'COLL')
       tndim=nspinor*ndim
       ABI_ALLOCATE(veejmj,(tndim,tndim,tndim,tndim))
       ABI_ALLOCATE(veeslm,(ndim,ndim,ndim,ndim))
       ABI_ALLOCATE(veetmp_slm,(ndim,ndim,4))
       ABI_ALLOCATE(veetmp_ylm,(ndim,ndim,4))
       ABI_ALLOCATE(veetmp,(ndim,ndim,4))
       ABI_ALLOCATE(veeslm2,(tndim,tndim,tndim,tndim))
       ABI_ALLOCATE(veeylm,(ndim,ndim,ndim,ndim))
       ABI_ALLOCATE(veeylm2,(tndim,tndim,tndim,tndim))
       ABI_ALLOCATE(veerotated,(tndim,tndim,tndim,tndim))
       ABI_ALLOCATE(fk,(0:lpawu))
       fk(0)=hu(itypat)%upawu
       if (lpawu==1) then
         fk(1)=hu(itypat)%jpawu*5
       else if (lpawu==2) then
         fk(1)=hu(itypat)%jpawu*14._dp/(One+hu(itypat)%f4of2_sla)
         fk(2)=fk(1)*hu(itypat)%f4of2_sla
       else if (lpawu==3) then
         fk(1)=hu(itypat)%jpawu*6435._dp/(286._dp+195._dp*hu(itypat)%f4of2_sla+250._dp*hu(itypat)%f6of2_sla)
         fk(2)=fk(1)*hu(itypat)%f4of2_sla
         fk(3)=fk(1)*hu(itypat)%f6of2_sla
       endif




!      ==================================
!      First define veeslm
!      ==================================


!!     veeslm2(s1m1,s2m2,s3m3,s4m4)=vee(m1,m2,m3,m4)*delta_s1s3*delta_s2s4
       call vee_ndim2tndim_hu(lpawu,hu(itypat)%vee,veeslm2,1)

!!     veeslm(m1,m2,m3,m4)=cmplx(vee(m1,m2,m3,m4),zero)
     !  ABI_DEALLOCATE(veeslm)
       veeslm(:,:,:,:)=cmplx(hu(itypat)%vee(:,:,:,:),zero)
       !ABI_DEALLOCATE(veeslm)! ici cela plante

!!     build udens in the Slm basis and print it
       call vee2udensatom_hu(ndim,nspinor,udens_atoms(iatom)%value,real(veeslm),"slm")
       !ABI_DEALLOCATE(veeslm) ! ici cela plante

       dim_vee=ndim
       basis_vee='Slm'
!     first print veeslm
       !call printvee_hu(ndim,real(veeslm),1,basis_vee)
       call printvee_hu(tndim,real(veeslm2),1,basis_vee)

      ! ABI_DEALLOCATE(veeslm)

!      ==================================
!      Then compute veerotated
!      ==================================

!      In the basis where levels/densitymatrix/green function is  diagonalized
!      ================================================================================
       if (rot_type==1) then
!      ---------------------

         veerotated=czero
         do m1=1,tndim
           do m2=1,tndim
             do m3=1,tndim
               do m4=1,tndim
                 do mi=1,tndim
                   do mj=1,tndim
                     do mk=1,tndim
                       do ml=1,tndim
                          veerotated(m1,m2,m3,m4)= veerotated(m1,m2,m3,m4) + &
&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
&                                rot_mat(iatom,1)%value(mk,m3)* &
&                                rot_mat(iatom,1)%value(ml,m4)* &
&                              veeslm2(mi,mj,mk,ml)
                       enddo
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         enddo

         dim_vee=tndim
         basis_vee='Diagonal basis from Slm basis'
!      In the Ylm basis
!      ================================================================================
       else if (rot_type==2.or.rot_type==3.or.rot_type==4) then
!      ---------------------------


!      Change basis from slm to ylm basis
         call vee_slm2ylm_hu(lpawu,veeslm,veeylm,1,2)

         basis_vee='Ylm'
         dim_vee=ndim

!        print interaction matrix in the ylm basis
         call printvee_hu(ndim,real(veeylm),1,basis_vee,hu(itypat)%upawu)

!        print interaction matrix in the ylm basis from Slater tables
         if(pawprtvol>=3) call udens_slatercondon_hu(fk,lpawu)

         if (rot_type==3.or.rot_type==4) then
!        ---------------------------
!

!          built large matrix
           call vee_ndim2tndim_hu(lpawu,real(veeylm),veeylm2,1)


           call printvee_hu(tndim,real(veeylm2),1,basis_vee)

         endif


!      In the JmJ basis
!      ================================================================================
         if (rot_type==3) then

!          apply change of basis
           call vee_ylm2jmj_hu(lpawu,veeylm2,veejmj,1)

!        print interaction matrix in the JMJ basis from Inglis and Julien tables
           if(pawprtvol>=3) call udens_inglis_hu(fk,lpawu)

!          new dimension



           dim_vee=tndim
           basis_vee='JmJ'

         else if (rot_type==4) then

           call udens_inglis_hu(fk,lpawu)
           veerotated=czero

           do m1=1,tndim
             do m2=1,tndim
               do m3=1,tndim
                 do m4=1,tndim
                   do mi=1,tndim
                     do mj=1,tndim
                       do mk=1,tndim
                         do ml=1,tndim
                            veerotated(m1,m2,m3,m4)= veerotated(m1,m2,m3,m4) + &
&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
&                                rot_mat(iatom,1)%value(mk,m3)* &
&                                rot_mat(iatom,1)%value(ml,m4)* &
&                                veeylm2(mi,mj,mk,ml)
                         enddo
                       enddo
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo

           dim_vee=tndim
           basis_vee='Diagonal basis from Ylm'

         endif
       endif
       ABI_ALLOCATE(veenew,(dim_vee,dim_vee,dim_vee,dim_vee))
       if(rot_type==0) then   ; veenew=veeslm
       else if(rot_type==1) then  ; veenew=veerotated
       else if(rot_type==2) then  ; veenew=veeylm
       else if(rot_type==3) then  ; veenew=veejmj
       else if(rot_type==4) then  ; veenew=veerotated
       endif

       call printvee_hu(dim_vee,real(veenew),1,basis_vee,hu(itypat)%upawu,hu(itypat)%f2_sla)
!       call printvee_hu(dim_vee,real(veeylm),1,hu(itypat)%upawu)

       uaver=zero
       do ms=1,2*ndim
         do ms1=1,2*ndim
           udens_atoms(iatom)%value(ms,ms1)=veenew(ms,ms1,ms,ms1)-veenew(ms,ms1,ms1,ms)
           uaver=uaver+udens_atoms(iatom)%value(ms,ms1)
         enddo
       enddo

       call vee2udensatom_hu(ndim,nspinor,udens_atoms(iatom)%value,real(veenew),"basis_vee",prtonly=1)


       ABI_DEALLOCATE(veenew)
       ABI_DEALLOCATE(veeslm)
       ABI_DEALLOCATE(veeslm2)
       ABI_DEALLOCATE(veeylm)
       ABI_DEALLOCATE(veeylm2)
       ABI_DEALLOCATE(veejmj)
       ABI_DEALLOCATE(veerotated)
       ABI_DEALLOCATE(veetmp_slm)
       ABI_DEALLOCATE(veetmp)
       ABI_DEALLOCATE(veetmp_ylm)
       ABI_DEALLOCATE(fk)
     endif
   enddo
   !MSG_ERROR("Aborting now!")

!================================================
!  NSPINOR = 1
!================================================
 else if (nspinor==1) then

   nat_correl=0
   do iatom=1,natom
     itypat=cryst_struc%typat(iatom)
     lpawu=hu(itypat)%lpawu
     if(lpawu.ne.-1) then
       nat_correl=nat_correl+1
       if(nat_correl>1.and.(hu(itypat)%jpawu>tol4)) then
          write(message,'(3a)')  ch10,'  -------> Warning: several atoms: '&
&         ,' not extensively tested '
          MSG_WARNING(message)
       endif

!  ! ================================================================
!  !  If rotation for spin 2 and rotation for spin 1 are not equal print
!  !  then print a warning
!  !  useful only for magnetic case
!  ! ================================================================
       ndim=2*lpawu+1
       tndim=nspinor*ndim
       do m1=1,tndim
         do m2=1,tndim
           if(nsppol==2) then
             if(abs(rot_mat(iatom,1)%value(m1,m2)-rot_mat(iatom,2)%value(m1,m2))>tol4.and.pawprtvol>=3) then
               write(message,'(2a,i4)')  ch10,' rot_mat differs for value of isppol but value for isppol=2 not used'
               call wrtout(std_out,  message,'COLL')
               write(message,'(a,4e16.8)')  ch10,rot_mat(iatom,1)%value(m1,m2),rot_mat(iatom,2)%value(m1,m2)
               call wrtout(std_out,  message,'COLL', do_flush=.True.)
             endif
           endif
         end do
       end do
       ABI_ALLOCATE(temp_mat,(ndim,ndim))
       ABI_ALLOCATE(temp_mat2,(ndim,ndim))
       ABI_ALLOCATE(veetemp,(ndim,ndim,ndim,ndim))
       temp_mat(:,:)=czero
       temp_mat2(:,:)=czero

!  ! =================================================
!  !    See if rotation is complex or real
!  ! =================================================
       do mi=1,ndim
         do m1=1,ndim
         if(abs(aimag(rot_mat(iatom,1)%value(mi,m1)))>tol8.and.pawprtvol>=3) then
              write(message,'(2a,2i6,2e14.3)')  ch10,"rot_mat is complex for", &
&              mi,m1,rot_mat(iatom,1)%value(mi,m1)
              call wrtout(std_out,  message,'COLL')
           endif
         enddo
       enddo



!  !    write vee for information with a classification.
       if(pawprtvol>=3) then
         call printvee_hu(ndim,hu(itypat)%vee,2,'original')
       endif

!  !    Compute rotated vee.
       veetemp=zero
       do m1=1,ndim
         do m2=1,ndim
           do m3=1,ndim
             do m4=1,ndim
               do mi=1,ndim
                 do mj=1,ndim
                   do mk=1,ndim
                     do ml=1,ndim
!                        if((mi==mk.and.mj==ml).or.(mi==ml.and.mj==mk)) then
                        veetemp(m1,m2,m3,m4)= veetemp(m1,m2,m3,m4) + &
&                      real(   &
&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
&                                rot_mat(iatom,1)%value(mk,m3)* &
&                                rot_mat(iatom,1)%value(ml,m4)* &
&                            hu(itypat)%vee(mi,mj,mk,ml)&
                           )
!                        endif
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           enddo
         enddo
       enddo
       ABI_DEALLOCATE(temp_mat)
       ABI_DEALLOCATE(temp_mat2)
       xsum=zero
       xsum2=zero
       xsumnew=zero
       xsum2new=zero
       do m1=1,ndim
         do m2=1,ndim
           xsum=xsum+hu(itypat)%vee(m1,m2,m1,m2)
           xsum2=xsum2+hu(itypat)%vee(m1,m2,m2,m1)
           xsumnew=xsumnew+veetemp(m1,m2,m1,m2)
           xsum2new=xsum2new+veetemp(m1,m2,m2,m1)
         enddo
       enddo
       if(abs(xsum-xsumnew)>tol5.or.abs(xsum2-xsum2new)>tol5) then
         write(message,'(2a)')  ch10," BUG: New interaction after rotation do not respect sum rules"
         call wrtout(std_out,  message,'COLL')
         write(message,'(2a,2f14.3)')  ch10,' Comparison of \sum_{m1,m3} vee(m1,m3,m1,m3) before and after rotation is',&
&         xsum,xsumnew
         call wrtout(std_out,  message,'COLL')
         write(message,'(2a,2f14.3)')  ch10,' Comparison of \sum_{m1,m3} vee(m1,m3,m3,m1) before and after rotation is',&
&         xsum2,xsum2new
         call wrtout(std_out,  message,'COLL')
       endif
       if(pawprtvol>=3) then
         write(message,'(2a)')  ch10," VEE ROTATED"
         call wrtout(std_out,  message,'COLL')
         call printvee_hu(tndim,veetemp,2,'diag1')
         write(message,'(a)') ch10
         call wrtout(std_out,  message,'COLL')
       endif

       write(message,'(2a,i4)')  ch10,'  -------> For Correlated Species', itypat
       call wrtout(std_out,  message,'COLL')

       call vee2udensatom_hu(ndim,nspinor,udens_atoms(iatom)%value,veetemp,"diag")

!       udens_atoms(iatom)%value=zero
!       ij=0
!       do ms=1,2*ndim-1
!         do ms1=ms+1,2*ndim
!           ij=ij+1
!           if(ms<=ndim.and.ms1>ndim) then
!             m1 = ms1 - ndim
!             m  = ms
!             hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)
!             udens_atoms(iatom)%value(ms,ms1)= veetemp(m,m1,m,m1)
!             udens_atoms(iatom)%value(ms1,ms)= udens_atoms(iatom)%value(ms,ms1)
!           else if(ms<=ndim.and.ms1<=ndim) then
!             m1 = ms1
!             m  = ms
!             hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
!             udens_atoms(iatom)%value(ms,ms1)= hu(itypat)%uqmc(ij)
!             udens_atoms(iatom)%value(ms1,ms)= udens_atoms(iatom)%value(ms,ms1)
!           else
!             m1 = ms1 - ndim
!             m  = ms  - ndim
!             hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
!             udens_atoms(iatom)%value(ms,ms1)= hu(itypat)%uqmc(ij)
!             udens_atoms(iatom)%value(ms1,ms)= udens_atoms(iatom)%value(ms,ms1)
!           endif
!         enddo
!       enddo
!       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation "
!       call wrtout(std_out,  message,'COLL')
!       write(message,'(1x,14(2x,i5))') (m,m=1,2*ndim)
!       call wrtout(std_out,  message,'COLL')
!       do ms=1,2*ndim
!          write(message,'(i3,14f7.3)') &
!  &        ms, (udens_atoms(iatom)%value(ms,ms1),ms1=1,2*ndim)
!          call wrtout(std_out,  message,'COLL')
!       enddo
!       write(message,'(5x,a)') "--------------------------------------------------------"
!       call wrtout(std_out,  message,'COLL')
       ABI_DEALLOCATE(veetemp)
     endif ! lpawu/=1
!   call print_hu(hu,cryst_struc%ntypat,1)

   enddo ! iatom
!   call print_hu(hu,cryst_struc%ntypat,1)
!   call vee2udens_hu(hu,cryst_struc%ntypat,2)
 endif ! nspinor==1


end subroutine rotatevee_hu
!!***

!!****f* m_hu/printvee_hu
!! NAME
!! printvee_hu
!!
!! FUNCTION
!!  print vee
!!
!! INPUTS
!!  vee = number of species
!!  lpawu = value of l
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine printvee_hu(ndim,vee,prtopt,basis,upawu,f2)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer,intent(in) :: ndim,prtopt
 real(dp), intent(in) :: vee(ndim,ndim,ndim,ndim)
 real(dp), optional, intent(in) :: upawu,f2
 character(len=*), intent(in) :: basis

!Local variables-------------------------------
 integer :: abcomp,m1,m2,mi,mj,mk,ml
 real(dp),allocatable :: b0(:,:)
 real(dp),allocatable :: a2pp(:,:)
 real(dp),allocatable :: b2pp(:,:)
 real(dp):: aver52,aver72,avernondiag,averall,averu,averj,averurestricted
 character(len=2000) :: message
! *********************************************************************
  write(message,'(5a)') ch10,&
&  '  Coulomb interaction in the ', trim(basis),' basis'
  call wrtout(std_out,message,'COLL')

 if(prtopt>=2) then

   write(message,'(2a)')  ch10," <mi,mi|vee|mi mi> : U1"
   call wrtout(std_out,  message,'COLL')
   do mi=1,ndim
     write(message,'(4i4,3x,e10.3)')   mi,mi,mi,mi,vee(mi,mi,mi,mi)
     call wrtout(std_out,  message,'COLL')
   enddo

   write(message,'(2a)')  ch10," <mi,mj|vee|mi mj> : U2"
   call wrtout(std_out,  message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)')   mi,mj,mi,mj,vee(mi,mj,mi,mj)
       call wrtout(std_out,  message,'COLL')
     enddo
   enddo

   write(message,'(2a)')  ch10," <mi,mj|vee|mj mi> : J"
   call wrtout(std_out,  message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)')   mi,mj,mj,mi,vee(mi,mj,mj,mi)
       call wrtout(std_out,  message,'COLL')
     enddo
   enddo

   write(message,'(2a)')  ch10," <mi,mi|vee|mj mj> : J"
   call wrtout(std_out,  message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)')   mi,mi,mj,mj,vee(mi,mi,mj,mj)
       call wrtout(std_out,  message,'COLL')
     enddo
   enddo

   write(message,'(2a)')  ch10," vee is non zero also for"
   call wrtout(std_out,  message,'COLL')
   do mi=1,ndim
     do mj=1,ndim
       do mk=1,ndim
         do ml=1,ndim
            if((.not.(mi==mk.and.mj==ml)).and.(.not.(mi==ml.and.mj==mk)).and.&
&               .not.(mi==mj.and.mk==ml)) then
              if(vee(mi,mj,mk,ml)>tol8) then
                write(message,'(4i4,3x,e10.3)')   mi,mj,mk,ml,vee(mi,mj,mk,ml)
                call wrtout(std_out,  message,'COLL')
              endif
            endif
         enddo
       enddo
     enddo
   enddo
   write(message,'(a)')  ch10
   call wrtout(std_out,  message,'COLL')

 endif

 if(prtopt>=1) then

   write(message,'(2x,a,3x,14f10.4)') "Um1m2=Vee(m1,m2,m1,m2)"
   call wrtout(std_out,  message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,  message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(vee(m1,m2,m1,m2),m2=1,ndim)
     call wrtout(std_out,  message,'COLL')
   enddo
   write(message,'(a)') ch10;  call wrtout(std_out,  message,'COLL')

   write(message,'(2x,a,3x,14f10.4)') "Jm1m2=Vee(m1,m2,m2,m1)"
   call wrtout(std_out,  message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,  message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(vee(m1,m2,m2,m1),m2=1,ndim)
     call wrtout(std_out,  message,'COLL')
   enddo
   write(message,'(a)') ch10;  call wrtout(std_out,  message,'COLL')

   write(message,'(2x,a,3x,14f10.4)') "Udens(m1,m2)"
   call wrtout(std_out,  message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,  message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1),m2=1,ndim)
     call wrtout(std_out,  message,'COLL')
   enddo
   write(message,'(a)') ch10;  call wrtout(std_out,  message,'COLL')

   if(prtopt>=3) then
     if(ndim==7) then
       averu=zero
       do m1=1,7
         do m2=1,7
           averu=vee(m1,m2,m1,m2)+averu
         enddo
       enddo
       averu=averu/49.d0
       averurestricted=zero
       do m1=1,7
         do m2=1,7
           if(m1/=m2) averurestricted=vee(m1,m2,m1,m2)+averurestricted
         enddo
       enddo
       averurestricted=averurestricted/42.d0
       averj=zero
       do m1=1,7
         do m2=1,7
           if(m1/=m2) averj=vee(m1,m2,m2,m1)+averj
         enddo
       enddo
       averj=averj/42
       averall=zero
       do m1=1,7
         do m2=1,7
           if(m1/=m2) averall=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+averall
         enddo
       enddo
       averall=averall/42
         !write(6,*) "averages Ylm U, U restricted,J, U-J",averu,averurestricted,averj,averall

     endif

     if(ndim==14.and.trim(basis)=='jmj') then
       aver52=zero
       do m1=1,6
         do m2=1,6
           aver52=vee(m1,m2,m1,m2)+aver52
         enddo
       enddo
       aver52=aver52/36.d0
       aver72=zero
       do m1=7,14
         do m2=7,14
            aver72=vee(m1,m2,m1,m2)+aver72
         enddo
       enddo
       aver72=aver72/64.d0
       avernondiag=zero
       do m1=1,6
         do m2=7,14
           avernondiag=vee(m1,m2,m1,m2)+avernondiag
         enddo
       enddo
       avernondiag=avernondiag/48.d0
       averall=zero
       do m1=1,14
         do m2=1,14
           averall=vee(m1,m2,m1,m2)+averall
         enddo
       enddo
       averall=averall/196.d0
         !write(6,*) "U averages",aver52,aver72,avernondiag,averall



       aver52=zero
       do m1=1,6
         do m2=1,6
           if(m1/=m2) aver52=vee(m1,m2,m2,m1)+aver52
         enddo
       enddo
       aver52=aver52/30.d0
       aver72=zero
       do m1=7,14
         do m2=7,14
           if(m1/=m2) aver72=vee(m1,m2,m2,m1)+aver72
         enddo
       enddo
       aver72=aver72/56.d0
       avernondiag=zero
       do m1=1,6
         do m2=7,14
           avernondiag=vee(m1,m2,m2,m1)+avernondiag
         enddo
       enddo
       avernondiag=avernondiag/48.d0
       averall=zero
       do m1=1,14
         do m2=1,14
           if(m1/=m2) averall=vee(m1,m2,m2,m1)+averall
         enddo
       enddo
       averall=averall/182.d0
         !write(6,*) "J averages",aver52,aver72,avernondiag,averall




       aver52=zero
       do m1=1,6
         do m2=1,6
           if(m1/=m2) aver52=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+aver52
         enddo
       enddo
       aver52=aver52/30.d0
       aver72=zero
       do m1=7,14
         do m2=7,14
           if(m1/=m2) aver72=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+aver72
         enddo
       enddo
       aver72=aver72/56.d0
       avernondiag=zero
       do m1=1,6
         do m2=7,14
           avernondiag=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+avernondiag
         enddo
       enddo
       avernondiag=avernondiag/48.d0
       averall=zero
       do m1=1,14
         do m2=1,14
           if(m1/=m2) averall=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+averall
         enddo
       enddo
       averall=averall/182.d0
       !write(6,*) "U-J averages",aver52,aver72,avernondiag,averall

     endif
   endif

   if (present(upawu)) then
     ABI_ALLOCATE(a2pp,(ndim,ndim))
     ABI_ALLOCATE(b2pp,(ndim,ndim))
     ABI_ALLOCATE(b0,(ndim,ndim))

!     write(message,'(2x,a,3x,14f10.4)') "For check with respect to Slater's paper"
!     call wrtout(std_out,  message,'COLL')
!     ABI_ALLOCATE(f0,(ndim,ndim,ndim,ndim))
!     write(message,'(2x,a,3x,14f10.4)') "Vee(m1,m2,m1,m2)-F0*ao(m1,m2)"
!     call wrtout(std_out,  message,'COLL')
!     write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
!     call wrtout(std_out,  message,'COLL')
!
!     do m1=1,ndim
!       write(message,'(2x,i4,3x,14f10.4)') m1,(vee(m1,m2,m1,m2)-upawu,m2=1,ndim)
!       call wrtout(std_out,  message,'COLL')
!     enddo
!
!     f0=zero
!     do m1=1,ndim
!     f0(m1,m1,m1,m1)=upawu
!     enddo
!     write(message,'(a)') ch10
!     call wrtout(std_out,  message,'COLL')
!     write(message,'(2x,a,3x,14f10.4)') "Vee(m1,m2,m2,m1)-F0*b0(m1,m2)"
!     call wrtout(std_out,  message,'COLL')
!     write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
!     call wrtout(std_out,  message,'COLL')
!     do m1=1,ndim
!       write(message,'(2x,i4,3x,14f10.4)') m1,(vee(m1,m2,m2,m1)-f0(m1,m2,m2,m1),m2=1,ndim)
!       call wrtout(std_out,  message,'COLL')
!     enddo
!     ABI_DEALLOCATE(f0)


     b0=zero
     do m1=1,ndim
     b0(m1,m1)=upawu
     enddo
     abcomp=0
     if(ndim==3.and.present(f2).and.(trim(basis)=='slm')) then
       a2pp(:,:) = RESHAPE((/1,-2,1,-2,4,-2,-1,-2,-1/),(/3,3/))
       a2pp=a2pp/25*f2+upawu
       b2pp(:,:) = RESHAPE((/1,3,6,3,4,3,6,3,1/),(/3,3/))
       b2pp=b2pp/25*f2+b0
       abcomp=1
     else if(ndim==6.and.present(f2).and.(trim(basis)=='jmj')) then
       ABI_ALLOCATE(a2pp,(6,6))
       ABI_ALLOCATE(b2pp,(6,6))
       a2pp(:,:)=RESHAPE((/0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,&
&       0,-1,1,1,-1,0,0,-1,1,1,-1,0,0,1,-1,-1,1/),(/6,6/))
       a2pp=a2pp/25*f2+upawu
       b2pp(:,:)=RESHAPE((/0,0,1,2,3,4,0,0,4,3,2,1,1,4,1,2,2,0,2,3,&
&       2,1,0,2,3,2,2,0,1,2,4,1,0,2,2,1/),(/6,6/))
       b2pp=b2pp/25*f2+b0
       abcomp=1
     endif
     if(mod(ndim,3)==0.and.present(f2).and.abcomp==1) then
       write(message,'(2x,a)') "Exact result for Umm is"
       call wrtout(std_out,message,'COLL')
       do m1=1,ndim
         write(message,'(2x,i4,3x,14f10.4)') m1,(a2pp(m1,m2),m2=1,ndim)
         call wrtout(std_out,  message,'COLL')
       enddo
       write(message,'(a)') ch10;  call wrtout(std_out,  message,'COLL')
       write(message,'(2x,a,3x,14f10.4)') "Exact result for Jmm is"
       call wrtout(std_out,message,'COLL')
       do m1=1,ndim
         write(message,'(2x,i4,3x,14f10.4)') m1,(b2pp(m1,m2),m2=1,ndim)
         call wrtout(std_out,message,'COLL')
       enddo
       write(message,'(a)') ch10;  call wrtout(std_out,  message,'COLL')
     endif
     ABI_DEALLOCATE(a2pp)
     ABI_DEALLOCATE(b2pp)
     ABI_DEALLOCATE(b0)
   endif

 endif

end subroutine printvee_hu
!!***

!!****f* m_hu/vee2udensatom_hu
!! NAME
!! vee2udensatom_hu
!!
!! FUNCTION
!!  compute density density interaction (used for DFT+DMFT)
!!
!! INPUTS
!!  ntypat = number of species
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee2udensatom_hu(ndim,nspinor,udens_atoms,veetemp,basis,prtonly)

 use defs_basis
 use m_crystal, only : crystal_t
 implicit none

!Arguments ------------------------------------
!type
 integer, intent(in):: ndim,nspinor
 real(dp),intent(out) :: udens_atoms(2*ndim,2*ndim)
 !real(dp), intent(in) :: veetemp(nspinor*ndim,nspinor*ndim,nspinor*ndim,nspinor*ndim)
 real(dp), intent(in) :: veetemp(ndim,ndim,ndim,ndim)
 character(len=*), intent(in) :: basis
 integer, intent(in), optional:: prtonly

!Local variables-------------------------------
 integer :: ij,ms,ms1,m,m1,prt
 character(len=1000) :: message
! *********************************************************************
 prt=1
 if(present(prtonly)) prt=2
 if (nspinor==1.or.nspinor==2.and.prt<=1) then
   udens_atoms=zero
   ij=0
   do ms=1,2*ndim-1
     do ms1=ms+1,2*ndim
       ij=ij+1
       if(ms<=ndim.and.ms1>ndim) then
         m1 = ms1 - ndim
         m  = ms
!         hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)
         udens_atoms(ms,ms1)= veetemp(m,m1,m,m1)
         udens_atoms(ms1,ms)= udens_atoms(ms,ms1)
        ! write(6,*)"A", ms,ms1,udens_atoms(ms,ms1)
       else if(ms<=ndim.and.ms1<=ndim) then
         m1 = ms1
         m  = ms
!         hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
         udens_atoms(ms,ms1)= veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
         udens_atoms(ms1,ms)= udens_atoms(ms,ms1)
        ! write(6,*)"B", ms,ms1,udens_atoms(ms,ms1)
       else
         m1 = ms1 - ndim
         m  = ms  - ndim
!         hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
         udens_atoms(ms,ms1)= veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
         udens_atoms(ms1,ms)= udens_atoms(ms,ms1)
        ! write(6,*)"C", ms,ms1,udens_atoms(ms,ms1)
       endif
     enddo
   enddo

! else if(nspinor==2) then
!
!   do ms=1,2*ndim
!     do ms1=1,2*ndim
!       udens_atoms(ms,ms1)=veetemp(ms,ms1,ms,ms1)-veetemp(ms,ms1,ms1,ms)
!     enddo
!   enddo

 endif


 write(message,'(4a)') ch10,"-------- Interactions in the ",trim(basis)," basis "
 call wrtout(std_out,  message,'COLL')
 write(message,'(1x,14(2x,i5))') (m,m=1,2*ndim)
 call wrtout(std_out,  message,'COLL')
 do ms=1,2*ndim
    write(message,'(i3,14f7.3)') &
&    ms, (udens_atoms(ms,ms1),ms1=1,2*ndim)
    call wrtout(std_out,  message,'COLL')
 enddo
 write(message,'(5x,a)') "--------------------------------------------------------"
 call wrtout(std_out,  message,'COLL')

end subroutine vee2udensatom_hu
!!***

!!****f* m_hu/reddd
!! NAME
!! reddd
!!
!! FUNCTION
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

function reddd(mi,ndim)

 use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mi,ndim
 integer :: reddd
! *************************************************************************

 if(mi<ndim+1)  reddd=mi
 if(mi>=ndim+1) reddd=mi-ndim

end function reddd
!!***

!!****f* m_hu/vee_slm2ylm_hu
!! NAME
!! vee_slm2ylm_hu
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix  of interaction of dimension (2*lcor+1)
!! from the Slm to the Ylm basis if option==1 or from Ylm to Slm if !option==2
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= Input matrix
!!  option= 1  Change matrix from Slm to Ylm basis
!!          2  Change matrix from Ylm to Slm basis
!!  optspin=  1  Spin up are first
!!            2  Spin dn are first
!!  prtvol=printing volume
!!
!! OUTPUT
!!  mat_inp_c= Output matrix in Ylm or Slm basis according to option
!!
!! NOTES
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee_slm2ylm_hu(lcor,mat_inp_c,mat_out_c,option,prtvol)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor,option,prtvol
!arrays
 complex(dpc), intent(in) :: mat_inp_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)
 complex(dpc), intent(out) :: mat_out_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)

!Local variables ---------------------------------------
!scalars
 integer :: gm,hm,jm,gg,hh,ii,jj,ll,mm,im
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
 complex(dpc) :: tmp2
 character(len=500) :: message
!arrays
 complex(dpc),allocatable :: slm2ylm(:,:)

! *********************************************************************


 if (option/=1.and.option/=2) then
   message=' option=/1 or 2 !'
   MSG_BUG(message)
 end if

 if(abs(prtvol)>2) then
   write(message,'(3a)') ch10, "   vee_slm2ylm_hu"
   call wrtout(std_out,message,'COLL')
 end if

 if(abs(prtvol)>2) then
   if(option==1) then
     write(message,'(3a)') ch10,"matrix in Slm basis is changed into Ylm basis"
     call wrtout(std_out,message,'COLL')
   else if(option==2) then
     write(message,'(3a)') ch10,"matrix in Ylm basis is changed into Slm basis"
     call wrtout(std_out,message,'COLL')
   end if
 end if

 ll=lcor
 ABI_ALLOCATE(slm2ylm,(2*ll+1,2*ll+1))
 slm2ylm=czero
 mat_out_c=czero

!  ===== Definitions of slm2ylm
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1   ! mmj=-mm
!! im is in {1,....2*ll+1}
!! mm is in {-ll,....+ll}
!! jm is in {2*ll+1,....,1}
   onem=dble((-1)**mm)
   if (mm> 0) then ! im in {ll+1,2ll+1} and jm in {ll+1,1}
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
! do im=1,2*ll+1
!   write(message,'(7(2f14.5))') (slm2ylm(im,jm),jm=1,2*ll+1)
!   call wrtout(std_out,message,'COLL')
! end do

!  ===== Definitions of slm2ylm
!!!!  pawtab(itypat)%vee(m11,m31,m21,m41)= <m11 m31| vee| m21 m41 >
!!!!  pawtab(itypat)%vee(m11,m21,m31,m41)= <m11 m21| vee| m31 m41 >

 do jm=1,2*ll+1
   do im=1,2*ll+1
     do hm=1,2*ll+1
       do gm=1,2*ll+1
         tmp2=czero
         do gg=1,2*ll+1
           do hh=1,2*ll+1
             do ii=1,2*ll+1
               do jj=1,2*ll+1
                 if(option==1) then
                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))&
&                                                   *(slm2ylm(gm,gg))*CONJG(slm2ylm(hm,hh))
!                   if(gm==1.and.hm==1.and.im==1.and.jm==1) then
!                      write(6,'(4i4,2f10.5,2f10.5)') gg,hh,ii,jj,tmp2,mat_inp_c(gg,hh,ii,jj)
!                      write(6,*) "i1"
!                   endif
                 else if(option==2) then
                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))&
&                                                   *CONJG(slm2ylm(gg,gm))*(slm2ylm(hh,hm))
                 end if
               end do
             end do
           end do
         end do
!         mat_out_c(gm,hm,im,jm)=tmp2
         mat_out_c(gm,im,hm,jm)=tmp2
       end do
     end do
   end do
 end do

 ABI_DEALLOCATE(slm2ylm)

end subroutine vee_slm2ylm_hu
!!***

!!****f* m_hu/vee_ndim2tndim_hu_r
!! NAME
!! vee_ndim2tndim_hu_r
!!
!! FUNCTION
!! Change a matrix  of interaction of dimension [(2*lcor+1)]**2
!! into a full spin and orbital interaction matrix of dimension [2*(2l+1)]**4
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= real input matrix
!!  prtvol=printing volume
!!  option= 1 : Vout_(s1m1,s2m2,s3m3,s4m4)=Vinp_(m1,m2,m3,m4)*delta_s1s3*delta_s2s4
!!
!!
!! OUTPUT
!!  mat_out_c= real output matrix
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee_ndim2tndim_hu_r(lcor,mat_inp_c,mat_out_c,option)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor,option
!arrays
 real(dp), intent(in)       :: mat_inp_c(:,:,:,:) !real(dp), intent(inout) :: mat_inp_c(:,:,:,:) !(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)  real(dp), intent(in), pointer :: mat_inp_c(:,:,:,:)
 real(dp), intent(out)       :: mat_out_c(:,:,:,:) !(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: m1,m2,m3,m4,is1,is2,is3,is4,ndim,s1,s2,s3,s4

! *********************************************************************
  ndim=2*lcor+1
  mat_out_c=czero

  if(option==1) then
    do m1=1,ndim
      do m2=1,ndim
        do m3=1,ndim
          do m4=1,ndim
            do is1=1,2
              do is2=1,2

                is3=is1 ; is4=is2

                s1=(is1-1)*ndim ; s2=(is2-1)*ndim ; s3=(is3-1)*ndim ; s4=(is4-1)*ndim

                mat_out_c(m1+s1,m2+s2,m3+s3,m4+s4)=  mat_inp_c(m1,m2,m3,m4)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif


end subroutine vee_ndim2tndim_hu_r
!!***

!!****f* m_hu/vee_ndim2tndim_hu
!! NAME
!! vee_ndim2tndim_hu
!!
!! FUNCTION
!! Change a matrix  of interaction of dimension [(2*lcor+1)]**2
!! into a full spin and orbital interaction matrix of dimension [2*(2l+1)]**4
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= real input matrix
!!  prtvol=printing volume
!!  option= 1 : Vout_(s1m1,s2m2,s3m3,s4m4)=Vinp_(m1,m2,m3,m4)*delta_s1s3*delta_s2s4
!!
!!
!! OUTPUT
!!  mat_out_c= Complex output matrix
!!
!! NOTES
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee_ndim2tndim_hu(lcor,mat_inp_c,mat_out_c,option)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor,option
!arrays
 real(dp), intent(in) :: mat_inp_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)
 complex(dpc), intent(out) :: mat_out_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: m1,m2,m3,m4,is1,is2,is3,is4,ndim,s1,s2,s3,s4

! *********************************************************************
  ndim=2*lcor+1
  mat_out_c=czero
  if(option==1) then
    do m1=1,ndim
      do m2=1,ndim
        do m3=1,ndim
          do m4=1,ndim
            do is1=1,2
              do is2=1,2

                is3=is1 ; is4=is2

                s1=(is1-1)*ndim ; s2=(is2-1)*ndim ; s3=(is3-1)*ndim ; s4=(is4-1)*ndim

                mat_out_c(m1+s1,m2+s2,m3+s3,m4+s4)=  cmplx(mat_inp_c(m1,m2,m3,m4),zero)

              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif


end subroutine vee_ndim2tndim_hu
!!***

!!****f* ABINIT/vee_ylm2jmj_hu
!! NAME
!! vee_ylm2jmj_hu
!!
!! FUNCTION
!! For a given angular momentum lcor, change a matrix  of dimension [2(2*lcor+1)]**4
!! from the Ylm basis to the J,M_J basis if option==1
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!  option=  1 matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis
!!           2 matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis
!!
!! SIDE EFFECTS
!!  mat_mlms= Input/Ouput matrix in the Ylm basis, size of the matrix is (2*lcor+1,2*lcor+1,ndij)
!!  mat_jmj= Input/Output matrix in the J,M_J basis, size is 2*(2*lcor+1),2*(2*lcor+1)
!!
!! NOTES
!!  usefull only in ndij==4
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine vee_ylm2jmj_hu(lcor,mat_inp_c,mat_out_c,option)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor,option
!arrays
 complex(dpc),intent(in) :: mat_inp_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))
 complex(dpc),intent(out) :: mat_out_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))

!Local variables ---------------------------------------
!scalars
 integer :: ii,im,jc1,jj,jm,ll,ml1,ms1,gg,hh,gm,hm
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: invsqrt2lp1,xj,xmj
 complex(dpc) :: tmp2
 character(len=500) :: message
!arrays
 integer, allocatable :: ind_msml(:,:)
 complex(dpc),allocatable :: mlms2jmj(:,:)

!*********************************************************************

 if (option/=1.and.option/=2) then
   message=' option=/1 and =/2 !'
   MSG_BUG(message)
 end if

 if(option==1) then
   write(message,'(3a)') ch10,&
&   "matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis"
   call wrtout(std_out,message)
 else if(option==2) then
   write(message,'(3a)') ch10,&
&   "matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis"
   call wrtout(std_out,message)
 end if

!--------------- Built indices + allocations
 ll=lcor
 ABI_ALLOCATE(mlms2jmj,(2*(2*ll+1),2*(2*ll+1)))
 mlms2jmj=czero
 ABI_ALLOCATE(ind_msml,(2,-ll:ll))
 mlms2jmj=czero
 jc1=0
 do ms1=1,2
   do ml1=-ll,ll
     jc1=jc1+1
     ind_msml(ms1,ml1)=jc1
   end do
 end do

!--------------- built mlms2jmj
!do jj=ll,ll+1    ! the physical value of j are ll-0.5,ll+0.5
!xj(jj)=jj-0.5
 if(ll==0)then
   message=' ll should not be equal to zero !'
   MSG_BUG(message)
 end if
 jc1=0
 invsqrt2lp1=one/sqrt(float(2*lcor+1))
 do jj=ll,ll+1
   xj=float(jj)-half !  xj is in {ll-0.5, ll+0.5}
   do jm=-jj,jj-1
     xmj=float(jm)+half  ! xmj is in {-xj,xj}
     jc1=jc1+1           ! Global index for JMJ
     if(nint(xj+0.5)==ll+1) then  ! if xj=ll+0.5
       if(nint(xmj+0.5)==ll+1)  then
         mlms2jmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
       else if(nint(xmj-0.5)==-ll-1) then
         mlms2jmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
       else
         mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
         mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
       end if
     end if
     if(nint(xj+0.5)==ll) then  ! if xj=ll-0.5
       mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
       mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
     end if
   end do
 end do
 write(message,'(3a)') ch10,"Matrix to go from |M_L,M_S> to |J,M_J>"
 call wrtout(std_out,message,"COLL")
 do im=1,2*(ll*2+1)
   write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mlms2jmj(im,jm),jm=1,2*(ll*2+1))
   call wrtout(std_out,message,"COLL")
 end do

!--------------- compute change of basis
 do jm=1,2*(2*ll+1)
   do im=1,2*(2*ll+1)
     do hm=1,2*(2*ll+1)
       do gm=1,2*(2*ll+1)
         tmp2=czero
         do gg=1,2*(2*ll+1)
           do hh=1,2*(2*ll+1)
             do ii=1,2*(2*ll+1)
               do jj=1,2*(2*ll+1)
                 if(option==1) then
                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))&
&                                                  *CONJG(mlms2jmj(gg,gm))*(mlms2jmj(hh,hm))
                 else if(option==2) then
!                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))& ! inv=t*
!&                                                  *CONJG(mlms2jmj(gg,gm))*(mlms2jmj(hh,hm)) ! inv=t*
                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*(mlms2jmj(im,ii))*CONJG(mlms2jmj(jm,jj))& ! inv=t*
&                                                  *(mlms2jmj(gm,gg))*CONJG(mlms2jmj(hm,hh)) ! inv=t*
                 end if
               end do
             end do
           end do
         end do
         mat_out_c(gm,im,hm,jm)=tmp2
       end do
     end do
   end do
 end do
 ABI_DEALLOCATE(mlms2jmj)
 ABI_DEALLOCATE(ind_msml)

 end subroutine vee_ylm2jmj_hu
!!***

!!****f* ABINIT/udens_slatercondon_hu
!! NAME
!! udens_slatercondon_hu
!!
!! FUNCTION
!! For a given angular momentum l and Slater integrals, give the
!! density density interactions U(m,m') and J(m,m') from Slater and
!! Condon tables
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!  fk(lcor+1)= Slater integrals
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine udens_slatercondon_hu(fk,lcor)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
 real(dp), intent(in) :: fk(0:lcor)

!Local variables ---------------------------------------
!scalars
 character(len=500) :: message
 integer :: m1,m2,kk
!arrays
 real(dp), allocatable :: aklmlmp(:,:,:)
 real(dp), allocatable :: bklmlmp(:,:,:)
 real(dp), allocatable :: udens(:,:)
 real(dp), allocatable :: jdens(:,:)

!*********************************************************************
 ABI_ALLOCATE(aklmlmp,(0:lcor,-lcor:lcor,-lcor:lcor)) ! k,m,m'
 ABI_ALLOCATE(bklmlmp,(0:lcor,-lcor:lcor,-lcor:lcor)) ! k,m,m'
 ABI_ALLOCATE(udens,(-lcor:lcor,-lcor:lcor)) ! m,m'
 ABI_ALLOCATE(jdens,(-lcor:lcor,-lcor:lcor)) ! m,m'
! k=2*(lcor)
 aklmlmp=zero
 bklmlmp=zero
 udens=zero
 jdens=zero
 if(lcor==0) then
   aklmlmp(0,0,0)=1
!
   bklmlmp(0,0,0)=0
 else if(lcor==1) then
   aklmlmp(0, :, :)=1
   aklmlmp(1, 1, 1)= one/25._dp
   aklmlmp(1,-1,-1)= one/25._dp
   aklmlmp(1,-1, 1)= one/25._dp
   aklmlmp(1, 1,-1)= one/25._dp
   aklmlmp(1, 1, 0)=-two/25._dp
   aklmlmp(1,-1, 0)=-two/25._dp
   aklmlmp(1, 0,-1)=-two/25._dp
   aklmlmp(1, 0, 1)=-two/25._dp
   aklmlmp(1, 0, 0)= four/25._dp
!
   bklmlmp(0, 1, 1)= one
   bklmlmp(0,-1,-1)= one
   bklmlmp(0, 0, 0)= one
   bklmlmp(1, 1, 1)= one/25._dp
   bklmlmp(1,-1,-1)= one/25._dp
   bklmlmp(1,-1,+1)= six/25._dp
   bklmlmp(1,+1,-1)= six/25._dp
   bklmlmp(1,+1, 0)= three/25._dp
   bklmlmp(1,-1, 0)= three/25._dp
   bklmlmp(1, 0,-1)= three/25._dp
   bklmlmp(1, 0, 1)= three/25._dp
   bklmlmp(1, 0, 0)= four/25._dp
 else if(lcor==2) then
   aklmlmp(0, :, :)=1
   aklmlmp(1, 2, 2)=  four / 49._dp
   aklmlmp(1, 2, 1)=  -two / 49._dp
   aklmlmp(1, 2, 0)= -four / 49._dp
   aklmlmp(1, 2,-1)=  -two / 49._dp
   aklmlmp(1, 2,-2)=  four / 49._dp
   aklmlmp(1, 1, 2)=  -two / 49._dp
   aklmlmp(1, 1, 1)=   one / 49._dp
   aklmlmp(1, 1, 0)=   two / 49._dp
   aklmlmp(1, 1,-1)=   one / 49._dp
   aklmlmp(1, 1,-2)=  -two / 49._dp
   aklmlmp(1, 0, 2)= -four / 49._dp
   aklmlmp(1, 0, 1)=   two / 49._dp
   aklmlmp(1, 0, 0)=  four / 49._dp
   aklmlmp(1, 0,-1)=   two / 49._dp
   aklmlmp(1, 0,-2)= -four / 49._dp
   aklmlmp(1,-1, 2)=  -two / 49._dp
   aklmlmp(1,-1, 1)=   one / 49._dp
   aklmlmp(1,-1, 0)=   two / 49._dp
   aklmlmp(1,-1,-1)=   one / 49._dp
   aklmlmp(1,-1,-2)=  -two / 49._dp
   aklmlmp(1,-2, 2)=  four / 49._dp
   aklmlmp(1,-2, 1)=  -two / 49._dp
   aklmlmp(1,-2, 0)= -four / 49._dp
   aklmlmp(1,-2,-1)=  -two / 49._dp
   aklmlmp(1,-2,-2)=  four / 49._dp

   aklmlmp(2, 2, 2)=    one  / 441._dp
   aklmlmp(2, 2, 1)=  -four  / 441._dp
   aklmlmp(2, 2, 0)=    six  / 441._dp
   aklmlmp(2, 2,-1)=  -four  / 441._dp
   aklmlmp(2, 2,-2)=    one  / 441._dp
   aklmlmp(2, 1, 2)=  -four  / 441._dp
   aklmlmp(2, 1, 1)=  16._dp / 441._dp
   aklmlmp(2, 1, 0)= -24._dp / 441._dp
   aklmlmp(2, 1,-1)=  16._dp / 441._dp
   aklmlmp(2, 1,-2)=  -four  / 441._dp
   aklmlmp(2, 0, 2)=    six  / 441._dp
   aklmlmp(2, 0, 1)= -24._dp / 441._dp
   aklmlmp(2, 0, 0)=  36._dp / 441._dp
   aklmlmp(2, 0,-1)= -24._dp / 441._dp
   aklmlmp(2, 0,-2)=    six  / 441._dp
   aklmlmp(2,-1, 2)=  -four  / 441._dp
   aklmlmp(2,-1, 1)=  16._dp / 441._dp
   aklmlmp(2,-1, 0)= -24._dp / 441._dp
   aklmlmp(2,-1,-1)=  16._dp / 441._dp
   aklmlmp(2,-1,-2)=  -four  / 441._dp
   aklmlmp(2,-2, 2)=    one  / 441._dp
   aklmlmp(2,-2, 1)=  -four  / 441._dp
   aklmlmp(2,-2, 0)=    six  / 441._dp
   aklmlmp(2,-2,-1)=  -four  / 441._dp
   aklmlmp(2,-2,-2)=    one  / 441._dp
   !do m1=lcor,-lcor,-1
   ! do m2=lcor,-lcor,-1
   !   write(6,*) m1,m2,aklmlmp(2,m1,m2)
   ! enddo
   !enddo

   !write(message,'(2x,a,3x,14f10.4)') " Slater aklmlmp(2,:,:)"
   !call wrtout(std_out,  message,'COLL')
   !write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
   !call wrtout(std_out,  message,'COLL')
   !do m1=-lcor,lcor,1
   !  write(message,'(2x,i4,3x,14f10.4)') m1,(aklmlmp(2,m1,m2),m2=-lcor,lcor,1)
   !  call wrtout(std_out,  message,'COLL')
   !enddo

   bklmlmp(0, 2, 2)=1
   bklmlmp(0, 1, 1)=1
   bklmlmp(0, 0, 0)=1
   bklmlmp(0,-1,-1)=1
   bklmlmp(0,-2,-2)=1
   bklmlmp(1, 2, 2)= four / 49._dp
   bklmlmp(1, 2, 1)=  six / 49._dp
   bklmlmp(1, 2, 0)= four / 49._dp
   bklmlmp(1, 2,-1)=  zero
   bklmlmp(1, 2,-2)=  zero
   bklmlmp(1, 1, 2)=  six / 49._dp
   bklmlmp(1, 1, 1)=  one / 49._dp
   bklmlmp(1, 1, 0)=  one / 49._dp
   bklmlmp(1, 1,-1)=  six / 49._dp
   bklmlmp(1, 1,-2)=  zero
   bklmlmp(1, 0, 2)= four / 49._dp
   bklmlmp(1, 0, 1)=  one / 49._dp
   bklmlmp(1, 0, 0)= four / 49._dp
   bklmlmp(1, 0,-1)=  one / 49._dp
   bklmlmp(1, 0,-2)= four / 49._dp
   bklmlmp(1,-1, 2)=  zero
   bklmlmp(1,-1, 1)=  six / 49._dp
   bklmlmp(1,-1, 0)=  one / 49._dp
   bklmlmp(1,-1,-1)=  one / 49._dp
   bklmlmp(1,-1,-2)=  six / 49._dp
   bklmlmp(1,-2, 2)=  zero
   bklmlmp(1,-2, 1)=  zero
   bklmlmp(1,-2, 0)= four / 49._dp
   bklmlmp(1,-2,-1)=  six / 49._dp
   bklmlmp(1,-2,-2)= four / 49._dp

   bklmlmp(2, 2, 2)=   one  / 441._dp
   bklmlmp(2, 2, 1)=  five  / 441._dp
   bklmlmp(2, 2, 0)= 15._dp / 441._dp
   bklmlmp(2, 2,-1)= 35._dp / 441._dp
   bklmlmp(2, 2,-2)= 70._dp / 441._dp
   bklmlmp(2, 1, 2)=  five  / 441._dp
   bklmlmp(2, 1, 1)= 16._dp / 441._dp
   bklmlmp(2, 1, 0)= 30._dp / 441._dp
   bklmlmp(2, 1,-1)= 40._dp / 441._dp
   bklmlmp(2, 1,-2)= 35._dp / 441._dp
   bklmlmp(2, 0, 2)= 15._dp / 441._dp
   bklmlmp(2, 0, 1)= 30._dp / 441._dp
   bklmlmp(2, 0, 0)= 36._dp / 441._dp
   bklmlmp(2, 0,-1)= 30._dp / 441._dp
   bklmlmp(2, 0,-2)= 15._dp / 441._dp
   bklmlmp(2,-1, 2)= 35._dp / 441._dp
   bklmlmp(2,-1, 1)= 40._dp / 441._dp
   bklmlmp(2,-1, 0)= 30._dp / 441._dp
   bklmlmp(2,-1,-1)= 16._dp / 441._dp
   bklmlmp(2,-1,-2)=  five  / 441._dp
   bklmlmp(2,-2, 2)= 70._dp / 441._dp
   bklmlmp(2,-2, 1)= 35._dp / 441._dp
   bklmlmp(2,-2, 0)= 15._dp / 441._dp
   bklmlmp(2,-2,-1)=  five  / 441._dp
   bklmlmp(2,-2,-2)=   one  / 441._dp

   !write(message,'(2x,a,3x,14f10.4)') " Slater bklmlmp(2,:,:)"
   !call wrtout(std_out,  message,'COLL')
   !write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
   !call wrtout(std_out,  message,'COLL')
   !do m1=-lcor,lcor,1
   !  write(message,'(2x,i4,3x,14f10.4)') m1,(bklmlmp(2,m1,m2),m2=-lcor,lcor,1)
   !  call wrtout(std_out,  message,'COLL')
   !enddo
         ! if f4of2_sla: these data agree with the explicit calculation
 else if(lcor==3) then
 endif

 do kk=0,lcor
   do m1=-lcor,lcor,1
     do m2=-lcor,lcor,1
       !write(6,*) kk,m1,m2
       !write(6,*) "--",fk(kk),aklmlmp(kk,m1,m2)
       !write(6,*) "--",fk(kk),bklmlmp(kk,m1,m2)
       udens(m1,m2)=udens(m1,m2)+fk(kk)*aklmlmp(kk,m1,m2)
       jdens(m1,m2)=jdens(m1,m2)+fk(kk)*bklmlmp(kk,m1,m2)
     enddo
   enddo
 enddo
 write(message,'(2x,a,3x,14f10.4)') " Direct Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,  message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,  message,'COLL')
 enddo

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Exchange Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,  message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(jdens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,  message,'COLL')
 enddo

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Density density Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,  message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2)-jdens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,  message,'COLL')
 enddo

 ABI_DEALLOCATE(jdens)
 ABI_DEALLOCATE(udens)
 ABI_DEALLOCATE(aklmlmp)
 ABI_DEALLOCATE(bklmlmp)

 end subroutine udens_slatercondon_hu
!!***

!!****f* ABINIT/udens_inglis_hu
!! NAME
!! udens_inglis_hu
!!
!! FUNCTION
!! For a given angular momentum l and Slater integrals, give the
!! density density interactions U(m,m') and J(m,m') from Inglis tables
!! in JMJ Basis
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!  fk(lcor+1)= Slater integrals
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_hu
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine udens_inglis_hu(fk,lcor)

 use defs_basis
 use m_errors
 use m_abicore
 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
 real(dp), intent(in) :: fk(0:lcor)

!Local variables ---------------------------------------
!scalars
 character(len=500) :: message
 integer :: m1,m2,kk,tndim
!arrays
 real(dp), allocatable :: udens(:,:)
 real(dp), allocatable :: jdens(:,:)
 real(dp),allocatable :: app(:,:,:)
 real(dp),allocatable :: bpp(:,:,:)
 real(dp),allocatable :: a2pp(:,:)
 real(dp),allocatable :: b2pp(:,:)

!*********************************************************************
 tndim=2*(2*lcor+1)
 ABI_ALLOCATE(app,(0:lcor,tndim,tndim))
 ABI_ALLOCATE(bpp,(0:lcor,tndim,tndim))
 ABI_ALLOCATE(a2pp,(tndim,tndim))
 ABI_ALLOCATE(b2pp,(tndim,tndim))

 ABI_ALLOCATE(udens,(tndim,tndim))
 ABI_ALLOCATE(jdens,(tndim,tndim))
 udens=zero
 jdens=zero
 a2pp=zero
 b2pp=zero
 app=zero
 bpp=zero
 if(lcor==1) then
   app(0,:,:)=one
   a2pp(:,:)=RESHAPE((/0._dp,0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
&                      0._dp,0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
&                      0._dp,0._dp, 1._dp,-1._dp,-1._dp, 1._dp,&
&                      0._dp,0._dp,-1._dp, 1._dp, 1._dp,-1._dp,&
&                      0._dp,0._dp,-1._dp, 1._dp, 1._dp,-1._dp,&
&                      0._dp,0._dp, 1._dp,-1._dp,-1._dp, 1._dp/),(/6,6/))
   app(1,:,:)=a2pp(:,:)
   app(1,:,:)=app(1,:,:)/25._dp

   b2pp(:,:)=RESHAPE((/1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,1._dp /),(/6,6/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/0._dp,0._dp,1._dp,2._dp,3._dp,4._dp,&
&                      0._dp,0._dp,4._dp,3._dp,2._dp,1._dp,&
&                      1._dp,4._dp,1._dp,2._dp,2._dp,0._dp,&
&                      2._dp,3._dp,2._dp,1._dp,0._dp,2._dp,&
&                      3._dp,2._dp,2._dp,0._dp,1._dp,2._dp,&
&                      4._dp,1._dp,0._dp,2._dp,2._dp,1._dp /),(/6,6/))
   bpp(1,:,:)=b2pp(:,:)
   bpp(1,:,:)=bpp(1,:,:)/25._dp
 endif
 if(lcor==2) then
   app(0,:,:)=one
   a2pp(:,:)=RESHAPE((/ 49._dp,-49._dp,-49._dp, 49._dp, 70._dp,-14._dp,-56._dp,-56._dp,-14._dp, 70._dp,&
&                      -49._dp, 49._dp, 49._dp,-49._dp,-70._dp, 14._dp, 56._dp, 56._dp, 14._dp,-70._dp,&
&                      -49._dp, 49._dp, 49._dp,-49._dp,-70._dp, 14._dp, 56._dp, 56._dp, 14._dp,-70._dp,&
&                       49._dp,-49._dp,-49._dp, 49._dp, 70._dp,-14._dp,-56._dp,-56._dp,-14._dp, 70._dp,&
&                       70._dp,-70._dp,-70._dp, 70._dp,100._dp,-20._dp,-80._dp,-80._dp,-20._dp,100._dp,&
&                      -14._dp, 14._dp, 14._dp,-14._dp,-20._dp,  4._dp, 16._dp, 16._dp,  4._dp,-20._dp,&
&                      -56._dp, 56._dp, 56._dp,-56._dp,-80._dp, 16._dp, 64._dp, 64._dp, 16._dp,-80._dp,&
&                      -56._dp, 56._dp, 56._dp,-56._dp,-80._dp, 16._dp, 64._dp, 64._dp, 16._dp,-80._dp,&
&                      -14._dp, 14._dp, 14._dp,-14._dp,-20._dp,  4._dp, 16._dp, 16._dp,  4._dp,-20._dp,&
&                       70._dp,-70._dp,-70._dp, 70._dp,100._dp,-20._dp,-80._dp,-80._dp,-20._dp,100._dp/),(/10,10/))
   app(1,:,:)=a2pp(:,:)/1225._dp
   app(2,:,:)=zero

   b2pp(:,:)=RESHAPE((/1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,&
&                      0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp/),(/10,10/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/ 49._dp, 98._dp, 98._dp,  0._dp, 30._dp, 36._dp, 27._dp, 12._dp,  0._dp,  0._dp,&
&                       98._dp, 49._dp,  0._dp, 98._dp, 40._dp,  2._dp,  6._dp, 25._dp, 32._dp,  0._dp,&
&                       98._dp,  0._dp, 49._dp, 98._dp,  0._dp, 32._dp, 25._dp,  6._dp,  2._dp, 40._dp,&
&                        0._dp, 98._dp, 98._dp, 49._dp,  0._dp,  0._dp, 12._dp, 27._dp, 36._dp, 30._dp,&
&                       30._dp, 40._dp,  0._dp,  0._dp,100._dp,120._dp, 60._dp,  0._dp,  0._dp,  0._dp,&
&                       36._dp,  2._dp, 32._dp,  0._dp,120._dp,  4._dp, 48._dp,108._dp,  0._dp,  0._dp,&
&                       27._dp,  6._dp, 25._dp, 12._dp, 60._dp, 48._dp, 64._dp,  0._dp,108._dp,  0._dp,&
&                       12._dp, 25._dp,  6._dp, 27._dp,  0._dp,108._dp,  0._dp, 64._dp, 48._dp, 60._dp,&
&                        0._dp, 32._dp,  2._dp, 36._dp,  0._dp,  0._dp,108._dp, 48._dp,  4._dp,120._dp,&
&                        0._dp,  0._dp, 40._dp, 30._dp,  0._dp,  0._dp,  0._dp, 60._dp,120._dp,100._dp/),(/10,10/))
   bpp(1,:,:)=b2pp(:,:)/1225._dp
   write(std_out,*) "warning: this test is only valid if f4of2_sla=0"
 endif
 if(lcor==3) then
!   app(0,:,:)=one
!   a2pp(:,:)=RESHAPE((/ 100.0/1225.0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 100.0/1225.0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0,&
!&                        0 ,0, 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0,&/),(/14,14/))
   b2pp(:,:)=RESHAPE((/  1.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   1.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   1.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   1.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   1.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0/),(/14,14/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/100.0, 120.0,  60.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                      120.0,   4.0,  48.0, 108.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                       60.0,  48.0,  64.0,   0.0, 108.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0, 108.0,   0.0,  64.0,  48.0,  60.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0, 108.0,  48.0,   4.0, 120.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,  60.0, 120.0, 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
&                        0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/),(/14,14/))
   bpp(1,:,:)=b2pp(:,:)/1225._dp
   write(std_out,*) "warning: only 5/2 5/2 elements are given and only for exchange and F2 !"
!   app(1,:,:)=a2pp(:,:)
!   app(1,:,:)=app(1,:,:)/25
!
!   bpp(0,:,:)=one
!   b2pp(:,:)=RESHAPE((/0,0,1,2,3,4,&
!&                      0,0,4,3,2,1,&
!&                      1,4,1,2,2,0,&
!&                      2,3,2,1,0,2,&
!&                      3,2,2,0,1,2,&
!&                      4,1,0,2,2,1 /),(/6,6/))
!   bpp(1,:,:)=b2pp(:,:)
!   bpp(1,:,:)=bpp(1,:,:)/25
 endif

 do kk=0,lcor
   do m1=1,tndim
     do m2=1,tndim
       !write(6,*) kk,m1,m2
       !write(6,*) "--",fk(kk),aklmlmp(kk,m1,m2)
       !write(6,*) "--",fk(kk),bklmlmp(kk,m1,m2)
       udens(m1,m2)=udens(m1,m2)+fk(kk)*app(kk,m1,m2)
       jdens(m1,m2)=jdens(m1,m2)+fk(kk)*bpp(kk,m1,m2)
     enddo
   enddo
 enddo
 write(message,'(2x,a,3x,14f10.4)') " Direct Interaction Matrix from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,  message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,  message,'COLL')
 enddo

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Exchange Interaction Matrix from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,  message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(jdens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,  message,'COLL')
 enddo

 write(message,'(a,2x,a,3x,14f10.4)') ch10, " Density Density interactions from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,  message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,  message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2)-jdens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,  message,'COLL')
 enddo


 ABI_DEALLOCATE(jdens)
 ABI_DEALLOCATE(udens)
 ABI_DEALLOCATE(app)
 ABI_DEALLOCATE(bpp)
 ABI_DEALLOCATE(a2pp)
 ABI_DEALLOCATE(b2pp)

 end subroutine udens_inglis_hu
!!***

END MODULE m_hu
!!***
