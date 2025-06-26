!!****m* ABINIT/m_hu
!! NAME
!!  m_hu
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2025 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

MODULE m_hu

 use defs_basis
 use m_abicore
 use m_errors

 use m_abi_linalg, only : abi_xgemm
 use m_matlu, only : matlu_type
 use m_paw_dmft, only : paw_dmft_type
 use m_pawtab, only : pawtab_type
 use m_crystal, only : crystal_t

 implicit none

 private

 public :: init_vee
 public :: destroy_vee
 public :: init_hu
 public :: copy_hu
 public :: destroy_hu
 !public :: qmc_hu
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

!!****t* m_hu/vee_type
!! NAME
!!  vee_type
!!
!! FUNCTION
!!  Structured datatype to store the U tensor
!!  for each atom.
!!
!! SOURCE

 type, public :: vee_type ! for each atom

   complex(dpc), allocatable :: mat(:,:,:,:)

 end type vee_type
!!***

!----------------------------------------------------------------------


!!****t* m_hu/hu_type
!! NAME
!!  hu_type
!!
!! FUNCTION
!!  This structured datatype contains interaction matrices for the correlated subspace
!!
!! SOURCE

 type, public :: hu_type ! for each typat

  integer  :: lpawu

  !logical  :: jmjbasis

  real(dp) :: upawu      ! => upaw

  real(dp) :: jpawu      ! => jpaw

  logical  :: jpawu_zero  ! true if all jpawu are zero
                          ! false if one of the jpaw is not zero

  real(dp), allocatable :: fk(:)

  complex(dpc), allocatable :: udens(:,:)

  complex(dpc), allocatable :: uqmc(:)

  complex(dpc), allocatable :: vee(:,:,:,:)

  complex(dpc), allocatable :: veeslm2(:,:,:,:)

 end type hu_type

!----------------------------------------------------------------------

CONTAINS  !========================================================================================
!!***

!!****f* m_hu/init_vee
!! NAME
!! init_vee
!!
!! FUNCTION
!!  Allocate variables used in type vee_type.
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  vee = tensor for the interactions
!!
!! OUTPUTS
!!
!! SOURCE

subroutine init_vee(paw_dmft,vee)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(vee_type), intent(inout) :: vee(paw_dmft%natom)
!Local variables ------------------------------------
 integer :: iatom,lpawu,ndim
!************************************************************************

 do iatom=1,paw_dmft%natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ndim = 2 * (2*lpawu+1)
   ABI_MALLOC(vee(iatom)%mat,(ndim,ndim,ndim,ndim))
   vee(iatom)%mat(:,:,:,:) = czero
 end do ! iatom

end subroutine init_vee
!!***

!!****f* m_hu/destroy_vee
!! NAME
!! destroy_vee
!!
!! FUNCTION
!!  Destroy vee
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  vee = tensor for the interactions
!!
!! OUTPUTS
!!
!! SOURCE

subroutine destroy_vee(paw_dmft,vee)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(vee_type), intent(inout) :: vee(paw_dmft%natom)
!Local variables ------------------------------------
 integer :: iatom,lpawu
!************************************************************************

 do iatom=1,paw_dmft%natom
   lpawu = paw_dmft%lpawu(iatom)
   if (lpawu == -1) cycle
   ABI_FREE(vee(iatom)%mat)
 end do ! iatom

end subroutine destroy_vee
!!***

!!****f* m_hu/init_hu
!! NAME
!! init_hu
!!
!! FUNCTION
!!  Allocate variables used in type hu_type.
!!
!! INPUTS
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawtab <type(pawtab)>=paw related data
!!
!! OUTPUTS
!!  hu <type(hu_type)>= U interaction
!!
!! SOURCE

subroutine init_hu(hu,paw_dmft,pawtab)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawtab_type), intent(in) :: pawtab(paw_dmft%ntypat)
 type(hu_type), intent(inout) :: hu(paw_dmft%ntypat)
!Local variables ------------------------------------
 integer  :: dmft_optim,i,ij,ij1,ij2,itypat,lpawu,m
 integer  :: m1,ms,ms1,ndim,ntypat,tndim
 logical  :: t2g,x2my2d
 real(dp) :: jpawu,upawu,xtemp
 integer, parameter   :: mt2g(3) = (/1,2,4/)
 integer, allocatable :: xij(:,:)
 character(len=4) :: tag
 character(len=500) :: message
!************************************************************************

 ntypat = paw_dmft%ntypat
 t2g    = (paw_dmft%dmft_t2g == 1)
 x2my2d = (paw_dmft%dmft_x2my2d == 1)

 dmft_optim = paw_dmft%dmft_optim

 write(message,'(2a)') ch10,"  == Compute Interactions for DMFT"
 call wrtout(std_out,message,'COLL')

 xtemp = zero

! ====================================
!  Compute hu(iatom)%uqmc from vee
! ====================================
 hu(1)%jpawu_zero = .true.
 do itypat=1,ntypat
   lpawu = pawtab(itypat)%lpawu
   hu(itypat)%upawu = zero
   hu(itypat)%jpawu = zero
   !hu(itypat)%jmjbasis = .false.
   if (t2g .and. lpawu == 2) lpawu = 1
   if (x2my2d .and. lpawu == 2) lpawu = 0
   hu(itypat)%lpawu = lpawu
   if (lpawu == -1) cycle
   ndim  = 2*lpawu + 1
   tndim = 2 * ndim
   hu(itypat)%upawu = pawtab(itypat)%upawu
   hu(itypat)%jpawu = pawtab(itypat)%jpawu

   if (hu(itypat)%jpawu > tol4) hu(1)%jpawu_zero = .false.
!     ndim1=2*hu(itypat)%lpawu+1

!     allocate(hu(itypat)%vee(ndim,ndim,ndim,ndim))

   ABI_MALLOC(hu(itypat)%vee,(ndim,ndim,ndim,ndim))

!  t2g case begin
   if (t2g) then
     do ms1=1,ndim
       do m1=1,ndim
         do ms=1,ndim
           do m=1,ndim
             hu(itypat)%vee(m,ms,m1,ms1) = &
               & cmplx(pawtab(itypat)%vee(mt2g(m),mt2g(ms),mt2g(m1),mt2g(ms1)),zero,kind=dp)
           end do ! m
         end do ! ms
       end do ! m1
     end do ! ms1
!   t2g case end
!   x2my2d case begin
   else if (x2my2d) then
     hu(itypat)%vee(1,1,1,1) = cmplx(pawtab(itypat)%upawu,zero,kind=dp)
!   x2my2d case end
   else
     hu(itypat)%vee(:,:,:,:) = cmplx(pawtab(itypat)%vee(:,:,:,:),zero,kind=dp)
   end if ! t2g or xymy2d
!   x2my2d case end

   ABI_MALLOC(hu(itypat)%veeslm2,(tndim,tndim,tndim,tndim))
   call vee_ndim2tndim_hu(lpawu,hu(itypat)%vee(:,:,:,:),hu(itypat)%veeslm2(:,:,:,:))

!  This is copied from pawpuxinit: it would be better not to duplicate
!  these lines.
   ABI_MALLOC(hu(itypat)%fk,(0:lpawu)) ! not used in the t2g and x2my2d cases
   hu(itypat)%fk(0) = hu(itypat)%upawu
   if (lpawu == 1) then
     hu(itypat)%fk(1) = hu(itypat)%jpawu * dble(5)
   else if (lpawu == 2) then
     hu(itypat)%fk(1) = hu(itypat)%jpawu * dble(14) / (one+pawtab(itypat)%f4of2_sla)
     hu(itypat)%fk(2) = hu(itypat)%fk(1) * pawtab(itypat)%f4of2_sla
   else if (lpawu == 3) then
     hu(itypat)%fk(1) = hu(itypat)%jpawu * dble(6435) / (dble(286)+&
           & dble(195)*pawtab(itypat)%f4of2_sla+dble(250)*pawtab(itypat)%f6of2_sla)
     hu(itypat)%fk(2) = hu(itypat)%fk(1) * pawtab(itypat)%f4of2_sla
     hu(itypat)%fk(3) = hu(itypat)%fk(1) * pawtab(itypat)%f6of2_sla
   end if ! lpawu

   write(tag,'(i4)') itypat
   write(message,'(3a)') ch10,'  -------> For Correlated Species ',adjustl(tag)
   call wrtout(std_out,message,'COLL')

   ABI_MALLOC(hu(itypat)%uqmc,(ndim*(tndim-1)))
   ABI_MALLOC(hu(itypat)%udens,(tndim,tndim))
   ABI_MALLOC(xij,(tndim,tndim))

   hu(itypat)%udens(:,:) = czero
   ij = 0
   do ms=1,tndim-1
     xij(ms,ms) = 0
     m = mod(ms-1,ndim) + 1
     do ms1=ms+1,tndim
       ij = ij + 1
       xij(ms,ms1) = ij
       xij(ms1,ms) = ij
       m1 = mod(ms1-1,ndim) + 1
       if ((ms <= ndim) .and. (ms1 > ndim)) then
         hu(itypat)%uqmc(ij) = hu(itypat)%vee(m,m1,m,m1)
       else
         hu(itypat)%uqmc(ij) = hu(itypat)%vee(m,m1,m,m1) - hu(itypat)%vee(m,m1,m1,m)
       end if
       hu(itypat)%udens(ms,ms1) = hu(itypat)%uqmc(ij)
       hu(itypat)%udens(ms1,ms) = hu(itypat)%udens(ms,ms1)
     end do ! ms1
   end do ! ms

   if (t2g .and. dmft_optim == 1) then
     upawu = zero
     jpawu = zero
     do ms1=1,ndim
       do ms=1,ndim
         upawu = upawu + dble(hu(itypat)%vee(ms,ms1,ms,ms1))
         jpawu = jpawu + dble(hu(itypat)%vee(ms,ms1,ms,ms1)-hu(itypat)%vee(ms,ms1,ms1,ms))
       end do ! ms
     end do ! ms1
     upawu = upawu / dble(ndim**2)
     jpawu = upawu - jpawu/dble(2*lpawu*ndim)
     hu(itypat)%upawu = upawu
     hu(itypat)%jpawu = jpawu
   end if ! t2g and dmft_optim=1

   xij(tndim,tndim) = 0
   write(message,'(a,5x,a)') ch10,"-------- Interactions in the density-density representation, in the cubic basis "
   call wrtout(std_out,message,'COLL')
   write(message,'(1x,14(2x,i5))') (m,m=1,tndim)
   call wrtout(std_out,message,'COLL')
!     xtemp1b=0.d0
! ====================================
!  Print hu(iatom)%uqmc
! ====================================
   ij2 = 0
   do i=1,tndim
     if (i < tndim) then
       ij1 = ij2 + 1
       ij2 = ij2 + tndim - i
     end if
!       write(std_out,*) itypat
!       do m=1,i
!        write(std_out,*) i,m
!        write(std_out,*) xij(i,m)
!        write(std_out,*) ij1,ij2
!       enddo
     if (i == 1) write(message,'(i3,14f7.3)') i,xtemp,(dble(hu(itypat)%uqmc(m)),m=ij1,ij2)
     if (i /= tndim .and. i /= 1) write(message,'(i3,14f7.3)') i, &
         & (dble(hu(itypat)%uqmc(xij(i,m))),m=1,i-1),xtemp,(dble(hu(itypat)%uqmc(m)),m=ij1,ij2)
     if (i == tndim) write(message,'(i3,14f7.3)') i,(dble(hu(itypat)%uqmc(xij(i,m))),m=1,i-1),xtemp
     call wrtout(std_out,message,'COLL')
   end do ! i
   write(message,'(5x,a)') "--------------------------------------------------------"
   call wrtout(std_out,message,'COLL')
   ABI_FREE(xij)

 end do ! itypat

end subroutine init_hu
!!***

!!****f* m_hu/copy_hu
!! NAME
!! copy_hu
!!
!! FUNCTION
!!  Copy hu into hu_new
!!
!! INPUTS
!!  hu <type(hu_type)>= U interaction
!!
!! OUTPUTS
!!  hu_new <type(hu_type)>= U interaction
!!
!! SOURCE

subroutine copy_hu(ntypat,hu,hu_new)

!Arguments ------------------------------------
!type
 integer, intent(in) :: ntypat
 type(hu_type), intent(in) :: hu(ntypat)
 type(hu_type), intent(inout) :: hu_new(ntypat)
!Local variables ------------------------------------
 integer :: itypat,ndim
!************************************************************************

 do itypat=1,ntypat
   hu_new(itypat)%lpawu      = hu(itypat)%lpawu
   !hu_new(itypat)%jmjbasis  = hu(itypat)%jmjbasis
   hu_new(itypat)%upawu      = hu(itypat)%upawu
   hu_new(itypat)%jpawu      = hu(itypat)%jpawu
   hu_new(itypat)%jpawu_zero = hu(itypat)%jpawu_zero
   ndim=2*hu_new(itypat)%lpawu+1
   ABI_MALLOC(hu_new(itypat)%uqmc,(ndim*(2*ndim-1)))
   ABI_MALLOC(hu_new(itypat)%udens,(2*ndim,2*ndim))
   ABI_MALLOC(hu_new(itypat)%vee,(ndim,ndim,ndim,ndim))
   ABI_MALLOC(hu_new(itypat)%fk,(0:hu_new(itypat)%lpawu))
   hu_new(itypat)%vee        = hu(itypat)%vee
   hu_new(itypat)%udens      = hu(itypat)%udens
   hu_new(itypat)%uqmc       = hu(itypat)%uqmc
   hu_new(itypat)%fk         = hu(itypat)%fk
 end do ! itypat

end subroutine copy_hu
!!***

!!****f* m_hu/destroy_hu
!! NAME
!! destroy_hu
!!
!! FUNCTION
!!  Deallocate hu
!!
!! INPUTS
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!  ntypat = number of species
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_hu(hu,ntypat)

!Arguments ------------------------------------
 integer, intent(in) :: ntypat
 type(hu_type), intent(inout) :: hu(ntypat)
!Local variables-------------------------------
 integer :: itypat
! *********************************************************************

 do itypat=1,ntypat
   ABI_SFREE(hu(itypat)%uqmc)
   ABI_SFREE(hu(itypat)%udens)
   ABI_SFREE(hu(itypat)%fk)
   ABI_SFREE(hu(itypat)%vee)
   ABI_SFREE(hu(itypat)%veeslm2)
 end do ! itypat

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
!! SOURCE

subroutine print_hu(hu,ntypat,prtopt)

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
       write(message,'(a,5x,a)') ch10,"-------- Interactions in the density matrix representation in cubic basis "
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
&          ms, (dble(hu(itypat)%udens(ms,ms1)),ms1=1,2*ndim)
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
!! SOURCE

subroutine vee2udens_hu(hu,ntypat,prtopt)

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
!!  Rotate U matrix to new basis
!!
!! INPUTS
!!  hu <type(hu_type)>= U interaction
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawprtvol >= 3 : print different quantities
!!  rot_mat = rotation matrix
!!  rot_type = 0 ! keep original Slm basis
!!           = 1 ! use the rotation matrix rot_mat from diago of dmat, green, levels..
!!           = 2 ! rotation to the Ylm basis
!!           = 3 ! rotation to the JmJ Basis
!!           = 4 ! same as 1 but rot_mat is applied from the Ylm basis instead of the Slm basis
!!
!! OUTPUT
!!  udens_atoms = rotated udens
!!  vee_rotated = rotated vee
!!
!! SIDE EFFECT
!!  hu <type(hu_type)> = data for the interaction in DMFT.
!!
!! SOURCE

subroutine rotatevee_hu(hu,paw_dmft,pawprtvol,rot_mat,rot_type,udens_atoms,vee_rotated)

!Arguments ------------------------------------
 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: pawprtvol,rot_type
 type(hu_type), target, intent(inout) :: hu(paw_dmft%ntypat)
 type(matlu_type), intent(in) :: rot_mat(paw_dmft%natom)
 type(matlu_type), intent(inout) :: udens_atoms(paw_dmft%natom)
 type(vee_type), target, intent(inout) :: vee_rotated(paw_dmft%natom)
!Local variables-------------------------------
 integer  :: dmft_optim,iatom,itypat,lpawu,m,m1,m2,mi,ms,ms1,nat_correl
 integer  :: natom,ndim,nflavor,nspinor,nsppol,nsppol_,prtonly,tndim
 logical  :: triqs
 real(dp) :: f2,jpawu,xsum,xsum2,xsum2new,xsumnew
 character(len=4) :: tag_at
 character(len=30) :: basis_vee
 character(len=500) :: message
 complex(dpc), target, allocatable :: veeylm(:,:,:,:)
 complex(dpc), pointer :: veeslm(:,:,:,:) => null(),veetemp(:,:,:,:) => null()
 complex(dpc), pointer :: veetemp2(:,:,:,:) => null(),veetemp3(:,:,:,:) => null()
 complex(dpc), pointer :: veeylm2(:,:,:,:) => null()
! *********************************************************************

 natom   = paw_dmft%natom
 nspinor = paw_dmft%nspinor
 nsppol  = paw_dmft%nsppol

 triqs = (paw_dmft%dmft_solv == 6 .or. paw_dmft%dmft_solv == 7)

 dmft_optim = paw_dmft%dmft_optim

 write(message,'(a,3x,a)') ch10,"== Rotate interaction to the CTQMC basis"
 call wrtout(std_out,message,"COLL")

!================================================
!  NSPINOR = 2
!================================================

 if (nspinor == 2) then

   do iatom=1,natom
     lpawu  = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     itypat = paw_dmft%typat(iatom)
     jpawu  = hu(itypat)%jpawu
     ndim   = 2*lpawu + 1
     tndim  = nspinor * ndim
      !if(pawprtvol>=3) then
!         write(message,'(2a)')  ch10," VEE INPUT AVANT TRANSFORMATION"
!         call wrtout(std_out,  message,'COLL')
!         call printvee_hu(ndim,hu(itypat)%vee,1,'Slm')
       !endif

     write(tag_at,'(i4)') iatom
     write(message,'(3a)') ch10,'  -------> For Correlated atom ',adjustl(tag_at)
     call wrtout(std_out,message,'COLL')

!    ==================================
!    First print veeslm
!    ==================================

     ! Print udens in the Slm basis
     call vee2udensatom_hu(ndim,hu(itypat)%udens(:,:),hu(itypat)%vee(:,:,:,:),"cubic",prtonly=1)

     basis_vee = 'cubic'
!    First print veeslm
       !call printvee_hu(ndim,real(veeslm),1,basis_vee)
     call printvee_hu(tndim,hu(itypat)%veeslm2(:,:,:,:),1,basis_vee)

!    ==================================
!    Then compute veerotated
!    ==================================

!    In the basis where levels/density matrix/green function is diagonal
!    ================================================================================

     ! When J=0, vee(:,i,:,j) and vee(i,:,j,:) are homotheties for all i,j
     ! when using Slater parametrization ; so no need to rotate them
     if (rot_type == 1 .and. jpawu > tol10) then
!      ---------------------

         !veerotated=czero
         !do m1=1,tndim
         !  do m2=1,tndim
         !    do m3=1,tndim
         !      do m4=1,tndim
         !        do mi=1,tndim
         !          do mj=1,tndim
         !            do mk=1,tndim
         !              do ml=1,tndim
         !                 veerotated(m1,m2,m3,m4)= veerotated(m1,m2,m3,m4) + &
!&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
!&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
!&                                rot_mat(iatom,1)%value(mk,m3)* &
!&                                rot_mat(iatom,1)%value(ml,m4)* &
!&                              veeslm2(mi,mj,mk,ml)
        !               enddo
        !             enddo
        !           enddo
        !         enddo
        !       enddo
        !     enddo
        !   enddo
        ! enddo

       call rotate_hu(rot_mat(iatom)%mat(:,:,:),1,tndim,hu(itypat)%veeslm2(:,:,:,:),vee_rotated(iatom)%mat(:,:,:,:))
       basis_vee = 'CTQMC basis from cubic'

!    In the Ylm basis
!    ================================================================================
     else if ((rot_type == 2 .or. rot_type == 3 .or. rot_type == 4) .and. jpawu > tol10) then
!    ---------------------------

       ABI_MALLOC(veeylm,(ndim,ndim,ndim,ndim))
       if (rot_type == 2) then
         veeylm2 => vee_rotated(iatom)%mat(:,:,:,:)
       else
         ABI_MALLOC(veeylm2,(tndim,tndim,tndim,tndim))
       end if
!      Change basis from slm to ylm basis
       if (dmft_optim == 1) then
         veeslm => hu(itypat)%vee(:,:,:,:)
       else
         ABI_MALLOC(veeslm,(ndim,ndim,ndim,ndim))
         veeslm(:,:,:,:) = cmplx(real(hu(itypat)%vee(:,:,:,:)),zero,kind=sp)
       end if

       call vee_slm2ylm_hu(lpawu,veeslm(:,:,:,:),veeylm(:,:,:,:),paw_dmft,1,2)

       if (dmft_optim == 0) then
         ABI_FREE(veeslm)
       end if
       veeslm => null()

       ! The line below is not really useful
       if (.not. triqs) veeylm(:,:,:,:) = cmplx(dble(veeylm(:,:,:,:)),zero,kind=dp)

       basis_vee = 'Ylm'

!      Print interaction matrix in the ylm basis
       call printvee_hu(ndim,veeylm(:,:,:,:),1,basis_vee,hu(itypat)%upawu)

!      Print interaction matrix in the ylm basis from Slater tables
       if (pawprtvol >= 3) then
         call udens_slatercondon_hu(hu(itypat)%fk(:),lpawu)
       end if

!      Build large matrix
       call vee_ndim2tndim_hu(lpawu,veeylm(:,:,:,:),veeylm2(:,:,:,:))

       if (rot_type == 3 .or. rot_type == 4) then
         call printvee_hu(tndim,veeylm2(:,:,:,:),1,basis_vee)
       end if

!      ---------------------------
!
!      In the JmJ basis
!      ================================================================================
       if (rot_type == 3) then

!        apply change of basis
         call vee_ylm2jmj_hu(lpawu,veeylm2(:,:,:,:),vee_rotated(iatom)%mat(:,:,:,:),1,paw_dmft)

!        print interaction matrix in the JMJ basis from Inglis and Julien tables
         if (pawprtvol >= 3) then
           call udens_inglis_hu(hu(itypat)%fk(:),lpawu)
         end if

!        new dimension

         basis_vee = 'JmJ'

       else if (rot_type == 4) then

           !veerotated=czero

           !do m1=1,tndim
           !  do m2=1,tndim
           !    do m3=1,tndim
           !      do m4=1,tndim
           !        do mi=1,tndim
           !          do mj=1,tndim
           !            do mk=1,tndim
           !              do ml=1,tndim
           !                 veerotated(m1,m2,m3,m4)= veerotated(m1,m2,m3,m4) +
           !                 &
!&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
!&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
!&                                rot_mat(iatom,1)%value(mk,m3)* &
!&                                rot_mat(iatom,1)%value(ml,m4)* &
!&                                veeylm2(mi,mj,mk,ml)
           !              enddo
           !            enddo
           !          enddo
           !        enddo
           !      enddo
           !    enddo
           !  enddo
           !enddo

         call udens_inglis_hu(hu(itypat)%fk(:),lpawu)
         call rotate_hu(rot_mat(iatom)%mat(:,:,:),1,tndim,veeylm2(:,:,:,:),vee_rotated(iatom)%mat(:,:,:,:))
         basis_vee = 'CTQMC basis from Ylm'

       end if ! rot_type

     end if ! rot_type

     if (rot_type == 0 .or. (jpawu <= tol10)) then
       vee_rotated(iatom)%mat(:,:,:,:) = hu(itypat)%veeslm2(:,:,:,:)
       udens_atoms(iatom)%mat(:,:,1)   = hu(itypat)%udens(:,:)
     end if ! rot_type

     ABI_SFREE(veeylm)
     if (rot_type >= 3 .and. jpawu > tol10) then
       ABI_FREE(veeylm2)
     end if
     veeylm2 => null()

     f2 = zero
     if (lpawu /= 0) f2 = hu(itypat)%fk(1)

     call printvee_hu(tndim,vee_rotated(iatom)%mat(:,:,:,:),1,basis_vee,hu(itypat)%upawu,f2)
!       call printvee_hu(dim_vee,real(veeylm),1,hu(itypat)%upawu)

       !uaver=zero
     if (rot_type /= 0 .and. jpawu > tol10) then
       ! Careful, since vee now has spin off-diagonal elements, this is wrong to
       ! use vee2udensatom_hu in order to compute udens
       do ms1=1,tndim
         do ms=1,tndim
           udens_atoms(iatom)%mat(ms,ms1,1) = vee_rotated(iatom)%mat(ms,ms1,ms,ms1) - vee_rotated(iatom)%mat(ms,ms1,ms1,ms)
           !uaver=uaver+udens_atoms(iatom)%value(ms,ms1)
         end do ! ms
       end do ! ms1
     end if ! rot_type /= 0 and jpawu /= zero

     call vee2udensatom_hu(ndim,udens_atoms(iatom)%mat(:,:,1),vee_rotated(iatom)%mat(:,:,:,:),basis_vee,prtonly=1)

   end do ! iatom
   !ABI_ERROR("Aborting now!")

!================================================
!  NSPINOR = 1
!================================================

 else if (nspinor == 1) then

   nat_correl = 0
   do iatom=1,natom
     lpawu = paw_dmft%lpawu(iatom)
     if (lpawu == -1) cycle
     itypat = paw_dmft%typat(iatom)
     jpawu  = hu(itypat)%jpawu
     nat_correl = nat_correl + 1
     if (nat_correl > 1 .and. (hu(itypat)%jpawu > tol4)) then
       write(message,'(3a)') ch10,'  -------> Warning: several atoms: ',' not extensively tested '
       ABI_WARNING(message)
     end if

     write(tag_at,'(i4)') iatom
     write(message,'(3a)') ch10,'  -------> For Correlated atom ',adjustl(tag_at)
     call wrtout(std_out,message,'COLL')

!  ! ================================================================
!  !  If rotation for spin 2 and rotation for spin 1 are not equal
!  !  then print a warning
!  !  useful only for magnetic case
!  ! ================================================================
     ndim  = 2*lpawu + 1
     tndim = nspinor * ndim
     nflavor = 2 * ndim
     if (nsppol == 2 .and. pawprtvol >= 3 .and. (.not. triqs)) then
       do m2=1,tndim
         do m1=1,tndim
           if (abs(rot_mat(iatom)%mat(m1,m2,1)-rot_mat(iatom)%mat(m1,m2,2)) > tol4) then
             write(message,'(2a,i4)') ch10,' rot_mat differs for value of isppol but value for isppol=2 not used'
             call wrtout(std_out,message,'COLL')
             write(message,'(a,4e16.8)') ch10,rot_mat(iatom)%mat(m1,m2,1),rot_mat(iatom)%mat(m1,m2,2)
             call wrtout(std_out,message,'COLL',do_flush=.True.)
           end if
         end do ! m1
       end do ! m2
     end if ! nsppol=2 and pawprtvol>=3

!  ! =================================================
!  !    See if rotation is complex or real
!  ! =================================================
     if (pawprtvol >= 3 .and. (.not. triqs)) then
       do m1=1,ndim
         do mi=1,ndim
           if (abs(aimag(rot_mat(iatom)%mat(mi,m1,1))) > tol8) then
             write(message,'(2a,2i6,2e14.3)') ch10,"rot_mat is complex for", &
               & mi,m1,rot_mat(iatom)%mat(mi,m1,1)
             call wrtout(std_out,message,'COLL')
           end if
         end do ! mi
       end do ! m1
     end if ! pawprtvol

!    write vee for information with a classification.
     if (pawprtvol >= 3) then
       call printvee_hu(ndim,hu(itypat)%vee(:,:,:,:),2,'cubic')
     end if

     basis_vee = 'cubic'
     prtonly = 1
     if (jpawu > tol10 .and. rot_type /= 0) then

!  !    Compute rotated vee.
       !veetemp=zero
       !do m1=1,ndim
       !  do m2=1,ndim
       !    do m3=1,ndim
       !      do m4=1,ndim
       !        do mi=1,ndim
       !          do mj=1,ndim
       !            do mk=1,ndim
       !              do ml=1,ndim
!      !                  if((mi==mk.and.mj==ml).or.(mi==ml.and.mj==mk)) then
       !                 veetemp(m1,m2,m3,m4)= veetemp(m1,m2,m3,m4) + &
!&                      real(   &
!&                          conjg(rot_mat(iatom,1)%value(mi,m1))* &
!&                          conjg(rot_mat(iatom,1)%value(mj,m2))* &
!&                                rot_mat(iatom,1)%value(mk,m3)* &
!&                                rot_mat(iatom,1)%value(ml,m4)* &
!&                            hu(itypat)%vee(mi,mj,mk,ml)&
                     !      )
!                        endif
    !                 enddo
    !               enddo
    !             enddo
    !           enddo
    !         enddo
    !       enddo
    !     enddo
    !   enddo

       nsppol_ = 1
       prtonly = 0

       ! Use different rotation matrices for each spin with TRIQS
       if (triqs .and. rot_type /= 2) nsppol_ = nsppol

       if (rot_type == 2 .or. rot_type == 4) then
         basis_vee = 'Ylm'
         ABI_MALLOC(veeylm,(ndim,ndim,ndim,ndim))
         call vee_slm2ylm_hu(lpawu,hu(itypat)%vee(:,:,:,:),veeylm(:,:,:,:),paw_dmft,1,2)
         if (rot_type == 2 .or. nsppol_ == 2) then
           if (rot_type == 2) then
             veeylm2 => vee_rotated(iatom)%mat(:,:,:,:)
           else if (nsppol_ == 2) then
             ABI_MALLOC(veeylm2,(nflavor,nflavor,nflavor,nflavor))
           end if
           call vee_ndim2tndim_hu(lpawu,veeylm(:,:,:,:),veeylm2(:,:,:,:))
         end if ! rot_type
         veetemp3 => veeylm(:,:,:,:)
       end if ! rot_type=2 or 4

       if (rot_type == 1 .or. rot_type == 4) then
         basis_vee = 'CTQMC basis from cubic'
         if (rot_type == 4) basis_vee = 'CTQMC basis from Ylm'
         if (nsppol_ == 2) then
           if (rot_type == 1) then
             veetemp => hu(itypat)%veeslm2(:,:,:,:)
           else
             veetemp => veeylm2
           end if ! rot_type
           veetemp2 => vee_rotated(iatom)%mat(:,:,:,:)
         else
           if (rot_type == 1) then
             veetemp => hu(itypat)%vee(:,:,:,:)
           else
             veetemp => veeylm(:,:,:,:)
           end if ! rot_type
           ABI_MALLOC(veetemp2,(ndim,ndim,ndim,ndim))
         end if ! nsppol_
         call rotate_hu(rot_mat(iatom)%mat(:,:,:),nsppol_,ndim,veetemp(:,:,:,:),veetemp2(:,:,:,:))
         if (.not. triqs) veetemp2(:,:,:,:) = cmplx(dble(veetemp2(:,:,:,:)),zero,kind=dp) ! neglect imaginary part in Abinit
         if (nsppol_ == 1) then
           call vee_ndim2tndim_hu(lpawu,veetemp2(:,:,:,:),vee_rotated(iatom)%mat(:,:,:,:))
           veetemp3 => veetemp2(:,:,:,:)
         end if
       end if ! rot_type=1 or 4

       if (nsppol_ == 2) then
         prtonly = 1
         ! It is wrong to use vee2udensatom_hu to build udens_atoms here
         do m=1,nflavor
           do m1=1,nflavor
             udens_atoms(iatom)%mat(m,m1,1) = vee_rotated(iatom)%mat(m,m1,m,m1) - vee_rotated(iatom)%mat(m,m1,m1,m)
           end do ! m1
         end do ! m
       end if ! nsppol_=2
     else
       prtonly = 1
       udens_atoms(iatom)%mat(:,:,1)   = hu(itypat)%udens(:,:)
       vee_rotated(iatom)%mat(:,:,:,:) = hu(itypat)%veeslm2(:,:,:,:)
     end if ! jpawu=zero

     xsum     = zero
     xsum2    = zero
     xsumnew  = zero
     xsum2new = zero
     do m1=1,ndim
       do m2=1,ndim
         xsum     = xsum + dble(hu(itypat)%vee(m1,m2,m1,m2))
         xsum2    = xsum2 + dble(hu(itypat)%vee(m1,m2,m2,m1))
         xsumnew  = xsumnew + dble(vee_rotated(iatom)%mat(m1,m2,m1,m2))
         xsum2new = xsum2new + dble(vee_rotated(iatom)%mat(m1,m2,m2,m1))
       end do ! m2
     end do ! m1
     if (abs(xsum-xsumnew) > tol5 .or. abs(xsum2-xsum2new) > tol5) then
       write(message,'(2a)') ch10," BUG: New interaction after rotation do not respect sum rules"
       call wrtout(std_out,message,'COLL')
       write(message,'(2a,2f14.3)') ch10,' Comparison of \sum_{m1,m3} vee(m1,m3,m1,m3) before and after rotation is',&
         & xsum,xsumnew
       call wrtout(std_out,message,'COLL')
       write(message,'(2a,2f14.3)') ch10,' Comparison of \sum_{m1,m3} vee(m1,m3,m3,m1) before and after rotation is',&
         & xsum2,xsum2new
       call wrtout(std_out,message,'COLL')
     end if ! abs(xsum-xsumnew)>tol5
     if (pawprtvol >= 3) then
       write(message,'(2a)') ch10," VEE ROTATED"
       call wrtout(std_out,message,'COLL')
       call printvee_hu(2*ndim,vee_rotated(iatom)%mat(:,:,:,:),2,'CTQMC')
       write(message,'(a)') ch10
       call wrtout(std_out,message,'COLL')
     end if ! pawprtvol>=3

     call vee2udensatom_hu(ndim,udens_atoms(iatom)%mat(:,:,1),veetemp3(:,:,:,:),basis_vee,prtonly=prtonly)

     veetemp  => null()
     veetemp3 => null()

     ABI_SFREE(veeylm)

     if (jpawu > tol10 .and. rot_type /= 0) then
       if (rot_type == 4 .and. nsppol_ == 2) then
         ABI_FREE(veeylm2)
       end if
       if (rot_type /= 2 .and. nsppol_ == 1) then
         ABI_FREE(veetemp2)
       end if
     end if ! jpawu>tol10 and rot_type/=0

     veetemp2 => null()
     veeylm2 => null()

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
       !ABI_FREE(veetemp)
    ! endif ! lpawu/=1
!   call print_hu(hu,cryst_struc%ntypat,1)

   end do ! iatom
!   call print_hu(hu,cryst_struc%ntypat,1)
!   call vee2udens_hu(hu,cryst_struc%ntypat,2)
 end if ! nspinor

end subroutine rotatevee_hu
!!***

!!****f* m_hu/rotate_hu
!! NAME
!! rotate_hu
!!
!! FUNCTION
!!  Rotate an interaction tensor
!!
!! INPUTS
!!  rot_mat = rotation matrix
!!  tndim = dimension of tensor
!!  vee = input tensor
!!
!! OUTPUT
!!  vee_rotated = rotated tensor
!!
!! SOURCE

subroutine rotate_hu(rot_mat,nsppol,tndim,vee,vee_rotated)

!Arguments ------------------------------------
 integer, intent(in) :: nsppol,tndim
 complex(dpc), intent(in) :: rot_mat(tndim,tndim,nsppol)
 complex(dpc), intent(in) :: vee(tndim*nsppol,tndim*nsppol,tndim*nsppol,tndim*nsppol)
 complex(dpc), intent(inout) :: vee_rotated(tndim*nsppol,tndim*nsppol,tndim*nsppol,tndim*nsppol)
!Local variables-------------------------------
 integer :: is1,is2,loop,m1,m2,ms1,ms2
 complex(dpc), allocatable :: mat_tmp(:,:),vee_tmp(:,:)
! *********************************************************************

 ABI_MALLOC(mat_tmp,(tndim,tndim))
 ABI_MALLOC(vee_tmp,(tndim,tndim))

 do loop=1,2
   do is2=1,nsppol
     do m2=1,tndim
       ms2 = m2 + (is2-1)*tndim
       do m1=1,tndim
         ms1 = m1 + (is2-1)*tndim
         do is1=1,nsppol

           ! Make copy here to prevent creation of temporary when calling zgemm
           if (loop == 1) then
             vee_tmp(:,:) = vee(1+(is1-1)*tndim:is1*tndim,ms1,1+(is1-1)*tndim:is1*tndim,ms2)
           else
             vee_tmp(:,:) = vee_rotated(ms1,1+(is1-1)*tndim:is1*tndim,ms2,1+(is1-1)*tndim:is1*tndim)
           end if ! loop

           call abi_xgemm("c","n",tndim,tndim,tndim,cone,rot_mat(:,:,is1),tndim, &
                        & vee_tmp(:,:),tndim,czero,mat_tmp(:,:),tndim)
           call abi_xgemm("n","n",tndim,tndim,tndim,cone,mat_tmp(:,:),tndim, &
                        & rot_mat(:,:,is1),tndim,czero,vee_tmp(:,:),tndim)

           if (loop == 1) then
             vee_rotated(1+(is1-1)*tndim:is1*tndim,ms1,1+(is1-1)*tndim:is1*tndim,ms2) = vee_tmp(:,:)
           else
             vee_rotated(ms1,1+(is1-1)*tndim:is1*tndim,ms2,1+(is1-1)*tndim:is1*tndim) = vee_tmp(:,:)
           end if ! loop

         end do ! is1
       end do ! m1
     end do ! ms2
   end do ! is2
 end do ! loop

 ABI_FREE(mat_tmp)
 ABI_FREE(vee_tmp)

end subroutine rotate_hu
!!***

!!****f* m_hu/printvee_hu
!! NAME
!! printvee_hu
!!
!! FUNCTION
!!  Print vee
!!
!! INPUTS
!!  vee = tensor for Coulomb interactions
!!
!! OUTPUT
!!
!! SOURCE

subroutine printvee_hu(ndim,vee,prtopt,basis,upawu,f2)

!Arguments ------------------------------------
!type
 integer, intent(in) :: ndim,prtopt
 complex(dpc), intent(in) :: vee(ndim,ndim,ndim,ndim)
 real(dp), optional, intent(in) :: f2,upawu
 character(len=*), intent(in) :: basis
!Local variables-------------------------------
 integer :: abcomp,m1,m2,mi,mj,mk,ml
 real(dp), allocatable :: a2pp(:,:),b0(:,:),b2pp(:,:)
 character(len=2000) :: message
! *********************************************************************

 write(message,'(5a)') ch10,&
     & '  Coulomb interaction in the ',trim(basis),' basis'
 call wrtout(std_out,message,'COLL')

 if (prtopt >= 2) then

   write(message,'(2a)') ch10," <mi,mi|vee|mi mi> : U1"
   call wrtout(std_out,message,'COLL')
   do mi=1,ndim
     write(message,'(4i4,3x,e10.3)') mi,mi,mi,mi,dble(vee(mi,mi,mi,mi))
     call wrtout(std_out,message,'COLL')
   end do ! mi

   write(message,'(2a)') ch10," <mi,mj|vee|mi mj> : U2"
   call wrtout(std_out,message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)') mi,mj,mi,mj,dble(vee(mi,mj,mi,mj))
       call wrtout(std_out,message,'COLL')
     end do ! mj
   end do ! mi

   write(message,'(2a)') ch10," <mi,mj|vee|mj mi> : J"
   call wrtout(std_out,message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)') mi,mj,mj,mi,dble(vee(mi,mj,mj,mi))
       call wrtout(std_out,message,'COLL')
     end do ! mj
   end do ! mi

   write(message,'(2a)') ch10," <mi,mi|vee|mj mj> : J"
   call wrtout(std_out,message,'COLL')
   do mi=1,ndim
     do mj=mi+1,ndim
       write(message,'(4i4,3x,e10.3)') mi,mi,mj,mj,dble(vee(mi,mi,mj,mj))
       call wrtout(std_out,message,'COLL')
     end do ! mj
   end do ! mi

   write(message,'(2a)') ch10," vee is non zero also for"
   call wrtout(std_out,message,'COLL')
   do mi=1,ndim
     do mj=1,ndim
       do mk=1,ndim
         do ml=1,ndim
           if ((.not. (mi == mk .and. mj == ml)) .and. (.not. (mi == ml .and. mj == mk)) .and.&
               & .not. (mi == mj .and. mk == ml)) then
             if (dble(vee(mi,mj,mk,ml)) > tol8) then
               write(message,'(4i4,3x,e10.3)') mi,mj,mk,ml,dble(vee(mi,mj,mk,ml))
               call wrtout(std_out,message,'COLL')
             end if
           end if
         end do ! ml
       end do ! mk
     end do ! mj
   end do ! mi
   write(message,'(a)') ch10
   call wrtout(std_out,message,'COLL')

 end if ! prtopt>=2

 if (prtopt >= 1) then

   write(message,'(2x,a,3x,14f10.4)') "Um1m2=Vee(m1,m2,m1,m2)"
   call wrtout(std_out,message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(dble(vee(m1,m2,m1,m2)),m2=1,ndim)
     call wrtout(std_out,message,'COLL')
   end do ! m1
   write(message,'(a)') ch10
   call wrtout(std_out,message,'COLL')

   write(message,'(2x,a,3x,14f10.4)') "Jm1m2=Vee(m1,m2,m2,m1)"
   call wrtout(std_out,message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(dble(vee(m1,m2,m2,m1)),m2=1,ndim)
     call wrtout(std_out,message,'COLL')
   end do ! m1
   write(message,'(a)') ch10
   call wrtout(std_out,message,'COLL')

   write(message,'(2x,a,3x,14f10.4)') "Udens(m1,m2)"
   call wrtout(std_out,message,'COLL')
   write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,ndim)
   call wrtout(std_out,message,'COLL')
   do m1=1,ndim
     write(message,'(2x,i4,3x,14f10.4)') m1,(dble(vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)),m2=1,ndim)
     call wrtout(std_out,message,'COLL')
   end do ! m1
   write(message,'(a)') ch10
   call wrtout(std_out,message,'COLL')

   !if (prtopt >= 3) then
   !  if(ndim==7) then
   !    averu=zero
   !    do m1=1,7
   !      do m2=1,7
   !        averu=vee(m1,m2,m1,m2)+averu
   !      enddo
   !    enddo
   !    averu=averu/49.d0
   !    averurestricted=zero
   !    do m1=1,7
   !      do m2=1,7
   !        if(m1/=m2) averurestricted=vee(m1,m2,m1,m2)+averurestricted
   !      enddo
   !    enddo
   !    averurestricted=averurestricted/42.d0
   !    averj=zero
   !    do m1=1,7
   !      do m2=1,7
   !        if(m1/=m2) averj=vee(m1,m2,m2,m1)+averj
   !      enddo
   !    enddo
   !    averj=averj/42
   !    averall=zero
   !    do m1=1,7
   !      do m2=1,7
   !        if(m1/=m2) averall=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+averall
   !      enddo
   !    enddo
   !    averall=averall/42
         !write(6,*) "averages Ylm U, U restricted,J, U-J",averu,averurestricted,averj,averall

   !  endif

   !  if(ndim==14.and.trim(basis)=='jmj') then
   !    aver52=zero
   !    do m1=1,6
   !      do m2=1,6
   !        aver52=vee(m1,m2,m1,m2)+aver52
   !      enddo
   !    enddo
    !   aver52=aver52/36.d0
    !   aver72=zero
    !   do m1=7,14
    !     do m2=7,14
    !        aver72=vee(m1,m2,m1,m2)+aver72
    !     enddo
    !   enddo
    !   aver72=aver72/64.d0
    !   avernondiag=zero
    !   do m1=1,6
    !     do m2=7,14
    !       avernondiag=vee(m1,m2,m1,m2)+avernondiag
    !     enddo
    !   enddo
    !   avernondiag=avernondiag/48.d0
    !   averall=zero
    !   do m1=1,14
    !     do m2=1,14
    !       averall=vee(m1,m2,m1,m2)+averall
    !     enddo
   !    enddo
    !   averall=averall/196.d0
         !write(6,*) "U averages",aver52,aver72,avernondiag,averall



    !   aver52=zero
    !   do m1=1,6
    !     do m2=1,6
    !       if(m1/=m2) aver52=vee(m1,m2,m2,m1)+aver52
    !     enddo
    !   enddo
    !   aver52=aver52/30.d0
    !   aver72=zero
    !   do m1=7,14
    !     do m2=7,14
    !       if(m1/=m2) aver72=vee(m1,m2,m2,m1)+aver72
    !     enddo
    !   enddo
    !   aver72=aver72/56.d0
    !   avernondiag=zero
    !   do m1=1,6
    !     do m2=7,14
    !       avernondiag=vee(m1,m2,m2,m1)+avernondiag
    !     enddo
    !   enddo
    !   avernondiag=avernondiag/48.d0
    !   averall=zero
    !   do m1=1,14
    !     do m2=1,14
    !       if(m1/=m2) averall=vee(m1,m2,m2,m1)+averall
    !     enddo
    !   enddo
    !   averall=averall/182.d0
         !write(6,*) "J averages",aver52,aver72,avernondiag,averall




    !   aver52=zero
    !   do m1=1,6
    !     do m2=1,6
    !       if(m1/=m2) aver52=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+aver52
    !     enddo
    !   enddo
   !    aver52=aver52/30.d0
    !   aver72=zero
   !    do m1=7,14
   !      do m2=7,14
   !        if(m1/=m2) aver72=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+aver72
   !      enddo
   !    enddo
    !   aver72=aver72/56.d0
    !   avernondiag=zero
    !   do m1=1,6
    !     do m2=7,14
    !       avernondiag=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+avernondiag
    !     enddo
    !   enddo
    !   avernondiag=avernondiag/48.d0
    !   averall=zero
    !   do m1=1,14
    !     do m2=1,14
    !       if(m1/=m2) averall=vee(m1,m2,m1,m2)-vee(m1,m2,m2,m1)+averall
    !     enddo
    !   enddo
    !   averall=averall/182.d0
       !write(6,*) "U-J averages",aver52,aver72,avernondiag,averall

   !  endif
   !endif

   if (present(upawu)) then

     ABI_MALLOC(a2pp,(ndim,ndim))
     ABI_MALLOC(b2pp,(ndim,ndim))
     ABI_MALLOC(b0,(ndim,ndim))

!     write(message,'(2x,a,3x,14f10.4)') "For check with respect to Slater's paper"
!     call wrtout(std_out,  message,'COLL')
!     ABI_MALLOC(f0,(ndim,ndim,ndim,ndim))
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
!     ABI_FREE(f0)


     b0(:,:) = zero
     do m1=1,ndim
       b0(m1,m1) = upawu
     end do
     abcomp = 0
     if (ndim == 3 .and. present(f2) .and. (trim(basis)=='slm')) then
       a2pp(:,:) = RESHAPE((/1,-2,1,-2,4,-2,-1,-2,-1/),(/3,3/))
       a2pp(:,:) = a2pp(:,:)/dble(25)*f2 + upawu
       b2pp(:,:) = RESHAPE((/1,3,6,3,4,3,6,3,1/),(/3,3/))
       b2pp(:,:) = b2pp(:,:)/dble(25)*f2 + b0(:,:)
       abcomp = 1
     else if (ndim == 6 .and. present(f2) .and. (trim(basis)=='jmj')) then
       ABI_MALLOC(a2pp,(6,6))
       ABI_MALLOC(b2pp,(6,6))
       a2pp(:,:) = RESHAPE((/0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,&
                  & 0,-1,1,1,-1,0,0,-1,1,1,-1,0,0,1,-1,-1,1/),(/6,6/))
       a2pp(:,:) = a2pp(:,:)/dble(25)*f2 + upawu
       b2pp(:,:) = RESHAPE((/0,0,1,2,3,4,0,0,4,3,2,1,1,4,1,2,2,0,2,3,&
                  & 2,1,0,2,3,2,2,0,1,2,4,1,0,2,2,1/),(/6,6/))
       b2pp(:,:) = b2pp(:,:)/dble(25)*f2 + b0(:,:)
       abcomp = 1
     end if
     if (mod(ndim,3) == 0 .and. present(f2) .and. abcomp == 1) then
       write(message,'(2x,a)') "Exact result for Umm is"
       call wrtout(std_out,message,'COLL')
       do m1=1,ndim
         write(message,'(2x,i4,3x,14f10.4)') m1,(a2pp(m1,m2),m2=1,ndim)
         call wrtout(std_out,message,'COLL')
       end do ! m1
       write(message,'(a)') ch10
       call wrtout(std_out,message,'COLL')
       write(message,'(2x,a,3x,14f10.4)') "Exact result for Jmm is"
       call wrtout(std_out,message,'COLL')
       do m1=1,ndim
         write(message,'(2x,i4,3x,14f10.4)') m1,(b2pp(m1,m2),m2=1,ndim)
         call wrtout(std_out,message,'COLL')
       end do ! m1
       write(message,'(a)') ch10
       call wrtout(std_out,message,'COLL')
     end if
     ABI_FREE(a2pp)
     ABI_FREE(b2pp)
     ABI_FREE(b0)
   end if ! present(upawu)

 end if ! prtopt>=1

end subroutine printvee_hu
!!***

!!****f* m_hu/vee2udensatom_hu
!! NAME
!! vee2udensatom_hu
!!
!! FUNCTION
!!  Compute and print density density interaction from full tensor (used for DFT+DMFT)
!!
!! INPUTS
!!  ndim = number of orbitals (without counting the spin)
!!  udens_atoms = density-density interactions
!!  veetemp = full interaction tensor
!!  basis = basis of the interaction tensor
!!  prtonly = 0 (default) : compute and print udens_atoms
!!          = 1 : only print udens_atoms
!!
!! OUTPUT
!!
!! SOURCE

subroutine vee2udensatom_hu(ndim,udens_atoms,veetemp,basis,prtonly)

!Arguments ------------------------------------
 integer, intent(in) :: ndim
 complex(dpc), intent(inout) :: udens_atoms(2*ndim,2*ndim)
 !real(dp), intent(in) :: veetemp(nspinor*ndim,nspinor*ndim,nspinor*ndim,nspinor*ndim)
 complex(dpc), intent(in) :: veetemp(ndim,ndim,ndim,ndim)
 character(len=*), intent(in) :: basis
 integer, intent(in), optional :: prtonly
!Local variables-------------------------------
 integer :: m,m1,ms,ms1,prt_only,tndim
 character(len=1000) :: message
! *********************************************************************

 tndim = 2 * ndim
 prt_only = 0
 if (present(prtonly)) prt_only = prtonly
 if (prt_only == 0) then
   udens_atoms(:,:) = czero
   !ij = 0
   do ms=1,tndim-1
     m = mod(ms-1,ndim) + 1
     do ms1=ms+1,tndim
       m1 = mod(ms1-1,ndim) + 1
       !ij = ij + 1
       if (ms <= ndim .and. ms1 > ndim) then
!         hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)
         udens_atoms(ms,ms1) = veetemp(m,m1,m,m1)
        ! write(6,*)"A", ms,ms1,udens_atoms(ms,ms1)
       else
!         hu(itypat)%uqmc(ij)=veetemp(m,m1,m,m1)-veetemp(m,m1,m1,m)
         udens_atoms(ms,ms1) = veetemp(m,m1,m,m1) - veetemp(m,m1,m1,m)
        ! write(6,*)"B", ms,ms1,udens_atoms(ms,ms1)
       end if
       udens_atoms(ms1,ms) = udens_atoms(ms,ms1)
     end do ! ms1
   end do ! ms

! else if(nspinor==2) then
!
!   do ms=1,2*ndim
!     do ms1=1,2*ndim
!       udens_atoms(ms,ms1)=veetemp(ms,ms1,ms,ms1)-veetemp(ms,ms1,ms1,ms)
!     enddo
!   enddo

 end if ! prt_only=0


 write(message,'(4a)') ch10,"-------- Interactions in the ",trim(basis)," basis "
 call wrtout(std_out,message,'COLL')
 write(message,'(1x,14(2x,i5))') (m,m=1,tndim)
 call wrtout(std_out,message,'COLL')
 do ms=1,tndim
   write(message,'(i3,14f7.3)') ms,(dble(udens_atoms(ms,ms1)),ms1=1,tndim)
   call wrtout(std_out,message,'COLL')
 end do ! ms
 write(message,'(5x,a)') "--------------------------------------------------------"
 call wrtout(std_out,message,'COLL')

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
!! SOURCE

!function reddd(mi,ndim)

! use defs_basis

!Arguments ------------------------------------
!scalars
! integer,intent(in) :: mi,ndim
! integer :: reddd
! *************************************************************************

! if(mi<ndim+1)  reddd=mi
! if(mi>=ndim+1) reddd=mi-ndim

!end function reddd
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
!! Copyright (C) 1998-2025 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2*lcor+1
!!  mat_inp_c= Input matrix
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  option= 1  Change matrix from Slm to Ylm basis
!!          2  Change matrix from Ylm to Slm basis
!!  prtvol=printing volume
!!
!! OUTPUT
!!  mat_out_c= Output matrix in Ylm or Slm basis according to option
!!
!! NOTES
!!
!! SOURCE

subroutine vee_slm2ylm_hu(lcor,mat_inp_c,mat_out_c,paw_dmft,option,prtvol)

!Arguments ---------------------------------------------
 integer,intent(in) :: lcor,option,prtvol
 complex(dpc), intent(in) :: mat_inp_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)
 complex(dpc), intent(out) :: mat_out_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)
 type(paw_dmft_type) , target, intent(in) :: paw_dmft
!Local variables ---------------------------------------
 integer :: ndim
 complex(dpc), allocatable :: slm2ylm(:,:)
 character(len=500) :: message
! *********************************************************************

 if (option /= 1 .and. option /= 2) then
   message = ' option=/1 or 2 !'
   ABI_BUG(message)
 end if

 if (abs(prtvol) > 2) then
   write(message,'(3a)') ch10,"   vee_slm2ylm_hu"
   call wrtout(std_out,message,'COLL')
 end if

 if (abs(prtvol) > 2) then
   if (option == 1) then
     write(message,'(3a)') ch10,"matrix in cubic basis is changed into Ylm basis"
   else if (option == 2) then
     write(message,'(3a)') ch10,"matrix in Ylm basis is changed into cubic basis"
   end if
   call wrtout(std_out,message,'COLL')
 end if ! prtvol>2

 ndim = 2*lcor + 1

 ABI_MALLOC(slm2ylm,(ndim,ndim))

 if (option == 1) then
   slm2ylm(:,:) = conjg(transpose(paw_dmft%slm2ylm(1:ndim,1:ndim,lcor+1)))
 else if (option == 2) then
   ! Make copy here to prevent creation of temporary in case ndim /= ndim_max
   slm2ylm(:,:) = paw_dmft%slm2ylm(1:ndim,1:ndim,lcor+1)
 end if
 call rotate_hu(slm2ylm(:,:),1,ndim,mat_inp_c(:,:,:,:),mat_out_c(:,:,:,:))

 ABI_FREE(slm2ylm)

 !ll=lcor
 !ABI_MALLOC(slm2ylm,(2*ll+1,2*ll+1))
 !slm2ylm=czero
 !mat_out_c=czero

!  ===== Definitions of slm2ylm
 !do im=1,2*ll+1
 !  mm=im-ll-1;jm=-mm+ll+1   ! mmj=-mm
!! im is in {1,....2*ll+1}
!! mm is in {-ll,....+ll}
!! jm is in {2*ll+1,....,1}
  ! onem=dble((-1)**mm)
  ! if (mm> 0) then ! im in {ll+1,2ll+1} and jm in {ll+1,1}
  !   slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
  !   slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
  ! end if
  ! if (mm==0) then
  !   slm2ylm(im,im)=cone
  ! end if
  ! if (mm< 0) then
  !   slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
  !   slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
  ! end if
! end do
! do im=1,2*ll+1
!   write(message,'(7(2f14.5))') (slm2ylm(im,jm),jm=1,2*ll+1)
!   call wrtout(std_out,message,'COLL')
! end do

!  ===== Definitions of slm2ylm
!!!!  pawtab(itypat)%vee(m11,m31,m21,m41)= <m11 m31| vee| m21 m41 >
!!!!  pawtab(itypat)%vee(m11,m21,m31,m41)= <m11 m21| vee| m31 m41 >

 !do jm=1,2*ll+1
 !  do im=1,2*ll+1
 !    do hm=1,2*ll+1
 !      do gm=1,2*ll+1
 !        tmp2=czero
 !        do gg=1,2*ll+1
 !          do hh=1,2*ll+1
 !            do ii=1,2*ll+1
 !              do jj=1,2*ll+1
 !                if(option==1) then
 !                  tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*(slm2ylm(im,ii))*CONJG(slm2ylm(jm,jj))&
!&                                                   *(slm2ylm(gm,gg))*CONJG(slm2ylm(hm,hh))
!                   if(gm==1.and.hm==1.and.im==1.and.jm==1) then
!                      write(6,'(4i4,2f10.5,2f10.5)') gg,hh,ii,jj,tmp2,mat_inp_c(gg,hh,ii,jj)
!                      write(6,*) "i1"
!                   endif
 !                else if(option==2) then
 !                  tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(slm2ylm(ii,im))*(slm2ylm(jj,jm))&
!&                                                   *CONJG(slm2ylm(gg,gm))*(slm2ylm(hh,hm))
 !                end if
 !              end do
 !            end do
 !          end do
 !        end do
!         mat_out_c(gm,hm,im,jm)=tmp2
 !        mat_out_c(gm,im,hm,jm)=tmp2
 !      end do
 !    end do
 !  end do
 !end do

 !ABI_FREE(slm2ylm)

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
!! Copyright (C) 1998-2025 ABINIT group (BA)
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
!! SOURCE

subroutine vee_ndim2tndim_hu_r(lcor,mat_inp_c,mat_out_c,option)

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
!! Change a matrix  of interaction of dimension [(2*lcor+1)]**4
!! into a full spin and orbital interaction matrix of dimension [2*(2l+1)]**4
!!
!! COPYRIGHT
!! Copyright (C) 1998-2025 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!  mat_inp_c= input matrix
!!
!! OUTPUT
!!  mat_out_c= Complex output matrix
!!
!! NOTES
!!
!! SOURCE

subroutine vee_ndim2tndim_hu(lcor,mat_inp_c,mat_out_c)

!Arguments ---------------------------------------------
!scalars
 integer, intent(in) :: lcor
!arrays
 complex(dpc), intent(in) :: mat_inp_c(2*lcor+1,2*lcor+1,2*lcor+1,2*lcor+1)
 complex(dpc), intent(inout) :: mat_out_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))
!Local variables ---------------------------------------
!scalars
 integer :: is1,is2,m1,m2,m3,m4,ndim,s1,s2
! *********************************************************************

 ndim = 2*lcor + 1
 mat_out_c(:,:,:,:) = czero

 do is2=1,2
   s2 = (is2-1) * ndim
   do m4=1,ndim
     do is1=1,2
       s1 = (is1-1) * ndim
       do m3=1,ndim
         do m2=1,ndim
           do m1=1,ndim
             mat_out_c(m1+s1,m2+s2,m3+s1,m4+s2) = mat_inp_c(m1,m2,m3,m4)
           end do ! m1
         end do ! m2
       end do ! m3
     end do ! is1
   end do ! m4
 end do ! is2

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
!! Copyright (C) 1998-2025 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum
!!  mat_inp_c = input tensor
!!  mat_out_c = output tensor
!!  option=  1 matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis
!!           2 matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  useful only in ndij==4
!!
!! SOURCE

subroutine vee_ylm2jmj_hu(lcor,mat_inp_c,mat_out_c,option,paw_dmft)

!Arguments ---------------------------------------------
 integer, intent(in) :: lcor,option
 complex(dpc), intent(in) :: mat_inp_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))
 complex(dpc), intent(inout) :: mat_out_c(2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1),2*(2*lcor+1))
 type(paw_dmft_type), intent(in) :: paw_dmft
!Local variables ---------------------------------------
 integer :: im,jm,tndim
 complex(dpc), allocatable :: jmj2ylm(:,:)
 character(len=500) :: message
!*********************************************************************

 if (option /= 1 .and. option /= 2) then
   message = ' option=/1 and =/2 !'
   ABI_BUG(message)
 end if

 if (lcor == 0) ABI_BUG("l should not be equal to 0")

 if (option == 1) then
   write(message,'(3a)') ch10,"matrix in |l,s,m_l,m_s> basis is changed into |l,s,j,m_j> basis"
 else if (option == 2) then
   write(message,'(3a)') ch10,"matrix in |l,s,j,m_j> basis is changed into |l,s,m_l,m_s> basis"
 end if
 call wrtout(std_out,message,"COLL")

 tndim = 2 * (2*lcor+1)

 ! Make copy to prevent creation of temporary in the case ndim /= ndim_max
 ABI_MALLOC(jmj2ylm,(tndim,tndim))

 if (option == 1) then
   jmj2ylm(:,:) = paw_dmft%jmj2ylm(1:tndim,1:tndim,lcor+1)
 else if (option == 2) then
   jmj2ylm(:,:) = conjg(transpose(paw_dmft%jmj2ylm(1:tndim,1:tndim,lcor+1)))
 end if

 write(message,'(3a)') ch10,"Matrix to go from |J,M_J> to |M_L,M_S>"
 call wrtout(std_out,message,"COLL")
 do im=1,tndim
   write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (jmj2ylm(im,jm),jm=1,tndim)
   call wrtout(std_out,message,"COLL")
 end do

 call rotate_hu(jmj2ylm(:,:),1,tndim,mat_inp_c(:,:,:,:),mat_out_c(:,:,:,:))

 ABI_FREE(jmj2ylm)

!--------------- Built indices + allocations
 !ll=lcor
 !ABI_MALLOC(mlms2jmj,(2*(2*ll+1),2*(2*ll+1)))
 !mlms2jmj=czero
 !ABI_MALLOC(ind_msml,(2,-ll:ll))
 !mlms2jmj=czero
 !jc1=0
 !do ms1=1,2
 !  do ml1=-ll,ll
 !    jc1=jc1+1
 !    ind_msml(ms1,ml1)=jc1
 !  end do
 !end do

!--------------- built mlms2jmj
!do jj=ll,ll+1    ! the physical value of j are ll-0.5,ll+0.5
!xj(jj)=jj-0.5
 !if(ll==0)then
 !  message=' ll should not be equal to zero !'
 !  ABI_BUG(message)
 !end if
 !jc1=0
 !invsqrt2lp1=one/sqrt(float(2*lcor+1))
 !do jj=ll,ll+1
 !  xj=float(jj)-half !  xj is in {ll-0.5, ll+0.5}
 !  do jm=-jj,jj-1
 !    xmj=float(jm)+half  ! xmj is in {-xj,xj}
 !    jc1=jc1+1           ! Global index for JMJ
 !    if(nint(xj+0.5)==ll+1) then  ! if xj=ll+0.5
 !      if(nint(xmj+0.5)==ll+1)  then
 !        mlms2jmj(ind_msml(2,ll),jc1)=1.0   !  J=L+0.5 and m_J=L+0.5
 !      else if(nint(xmj-0.5)==-ll-1) then
 !        mlms2jmj(ind_msml(1,-ll),jc1)=1.0   !  J=L+0.5 and m_J=-L-0.5
 !      else
 !        mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
 !        mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
 !      end if
 !    end if
 !    if(nint(xj+0.5)==ll) then  ! if xj=ll-0.5
 !      mlms2jmj(ind_msml(1,nint(xmj+0.5)),jc1)=invsqrt2lp1*(sqrt(float(ll)+xmj+0.5))
 !      mlms2jmj(ind_msml(2,nint(xmj-0.5)),jc1)=-invsqrt2lp1*(sqrt(float(ll)-xmj+0.5))
 !    end if
 !  end do
 !end do
 !write(message,'(3a)') ch10,"Matrix to go from |M_L,M_S> to |J,M_J>"
 !call wrtout(std_out,message,"COLL")
 !do im=1,2*(ll*2+1)
 !  write(message,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') (mlms2jmj(im,jm),jm=1,2*(ll*2+1))
 !  call wrtout(std_out,message,"COLL")
 !end do

!--------------- compute change of basis
! do jm=1,2*(2*ll+1)
!   do im=1,2*(2*ll+1)
!     do hm=1,2*(2*ll+1)
!       do gm=1,2*(2*ll+1)
!         tmp2=czero
!         do gg=1,2*(2*ll+1)
!           do hh=1,2*(2*ll+1)
!             do ii=1,2*(2*ll+1)
!               do jj=1,2*(2*ll+1)
!                 if(option==1) then
!                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))&
!&                                                  *CONJG(mlms2jmj(gg,gm))*(mlms2jmj(hh,hm))
!                 else if(option==2) then
!                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*CONJG(mlms2jmj(ii,im))*(mlms2jmj(jj,jm))& ! inv=t*
!&                                                  *CONJG(mlms2jmj(gg,gm))*(mlms2jmj(hh,hm)) ! inv=t*
!                   tmp2=tmp2+mat_inp_c(gg,ii,hh,jj)*(mlms2jmj(im,ii))*CONJG(mlms2jmj(jm,jj))& ! inv=t*
!&                                                  *(mlms2jmj(gm,gg))*CONJG(mlms2jmj(hm,hh)) ! inv=t*
 !                end if
 !              end do
!             end do
!           end do
!         end do
!         mat_out_c(gm,im,hm,jm)=tmp2
!       end do
!     end do
!   end do
! end do
! ABI_FREE(mlms2jmj)
! ABI_FREE(ind_msml)

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
!! Copyright (C) 1998-2025 ABINIT group (BA)
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
!! SOURCE

subroutine udens_slatercondon_hu(fk,lcor)

!Arguments ---------------------------------------------
!scalars
 integer,  intent(in) :: lcor
 real(dp), intent(in) :: fk(0:lcor)
!Local variables ---------------------------------------
!scalars
 character(len=500) :: message
 integer :: m1,m2
!arrays
 real(dp), allocatable :: aklmlmp(:,:,:),bklmlmp(:,:,:),jdens(:,:),udens(:,:)
!*********************************************************************

 ABI_MALLOC(aklmlmp,(0:lcor,-lcor:lcor,-lcor:lcor)) ! k,m,m'
 ABI_MALLOC(bklmlmp,(0:lcor,-lcor:lcor,-lcor:lcor)) ! k,m,m'
 ABI_MALLOC(udens,(-lcor:lcor,-lcor:lcor)) ! m,m'
 ABI_MALLOC(jdens,(-lcor:lcor,-lcor:lcor)) ! m,m'
! k=2*(lcor)
 aklmlmp(:,:,:) = zero
 bklmlmp(:,:,:) = zero
 udens(:,:) = zero
 jdens(:,:) = zero
 if (lcor == 0) then
   aklmlmp(0,0,0)=1
!
   bklmlmp(0,0,0)=0
 else if (lcor == 1) then
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
 else if (lcor == 2) then
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
 !else if(lcor==3) then
 end if ! lcor

 do m2=-lcor,lcor,1
   do m1=-lcor,lcor,1
     udens(m1,m2) = sum(fk(:)*aklmlmp(:,m1,m2))
     jdens(m1,m2) = sum(fk(:)*bklmlmp(:,m1,m2))
       !write(6,*) kk,m1,m2
       !write(6,*) "--",fk(kk),aklmlmp(kk,m1,m2)
       !write(6,*) "--",fk(kk),bklmlmp(kk,m1,m2)
       !udens(m1,m2)=udens(m1,m2)+fk(kk)*aklmlmp(kk,m1,m2)
       !jdens(m1,m2)=jdens(m1,m2)+fk(kk)*bklmlmp(kk,m1,m2)
   end do ! m1
 end do ! m2
 write(message,'(2x,a,3x,14f10.4)') " Direct Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Exchange Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(jdens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Density density Interaction Matrix from Slater tables (in the Ylm basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=-lcor,lcor,1)
 call wrtout(std_out,message,'COLL')
 do m1=-lcor,lcor,1
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2)-jdens(m1,m2),m2=-lcor,lcor,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1

 ABI_FREE(jdens)
 ABI_FREE(udens)
 ABI_FREE(aklmlmp)
 ABI_FREE(bklmlmp)

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
!! Copyright (C) 1998-2025 ABINIT group (BA)
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
!! SOURCE

subroutine udens_inglis_hu(fk,lcor)

!Arguments ---------------------------------------------
!scalars
 integer,  intent(in) :: lcor
 real(dp), intent(in) :: fk(0:lcor)
!Local variables ---------------------------------------
!scalars
 character(len=500) :: message
 integer :: m1,m2,tndim
!arrays
 real(dp), allocatable :: a2pp(:,:),app(:,:,:),b2pp(:,:),bpp(:,:,:),jdens(:,:),udens(:,:)
!*********************************************************************

 tndim = 2 * (2*lcor+1)
 ABI_MALLOC(app,(0:lcor,tndim,tndim))
 ABI_MALLOC(bpp,(0:lcor,tndim,tndim))
 ABI_MALLOC(a2pp,(tndim,tndim))
 ABI_MALLOC(b2pp,(tndim,tndim))

 ABI_MALLOC(udens,(tndim,tndim))
 ABI_MALLOC(jdens,(tndim,tndim))

 udens(:,:) = zero
 jdens(:,:) = zero
 a2pp(:,:)  = zero
 b2pp(:,:)  = zero
 app(:,:,:) = zero
 bpp(:,:,:) = zero
 if (lcor == 1) then
   app(0,:,:)=one
   a2pp(:,:)=RESHAPE((/0._dp,0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
                    &  0._dp,0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
                    &  0._dp,0._dp, 1._dp,-1._dp,-1._dp, 1._dp,&
                    &  0._dp,0._dp,-1._dp, 1._dp, 1._dp,-1._dp,&
                    &  0._dp,0._dp,-1._dp, 1._dp, 1._dp,-1._dp,&
                    &  0._dp,0._dp, 1._dp,-1._dp,-1._dp, 1._dp/),(/6,6/))
   app(1,:,:)=a2pp(:,:)
   app(1,:,:)=app(1,:,:)/25._dp

   b2pp(:,:)=RESHAPE((/1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,1._dp /),(/6,6/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/0._dp,0._dp,1._dp,2._dp,3._dp,4._dp,&
                     & 0._dp,0._dp,4._dp,3._dp,2._dp,1._dp,&
                     & 1._dp,4._dp,1._dp,2._dp,2._dp,0._dp,&
                     & 2._dp,3._dp,2._dp,1._dp,0._dp,2._dp,&
                     & 3._dp,2._dp,2._dp,0._dp,1._dp,2._dp,&
                     & 4._dp,1._dp,0._dp,2._dp,2._dp,1._dp /),(/6,6/))
   bpp(1,:,:)=b2pp(:,:)
   bpp(1,:,:)=bpp(1,:,:)/25._dp
 end if ! lcor=1
 if (lcor == 2) then
   app(0,:,:)=one
   a2pp(:,:)=RESHAPE((/ 49._dp,-49._dp,-49._dp, 49._dp, 70._dp,-14._dp,-56._dp,-56._dp,-14._dp, 70._dp,&
                     & -49._dp, 49._dp, 49._dp,-49._dp,-70._dp, 14._dp, 56._dp, 56._dp, 14._dp,-70._dp,&
                     & -49._dp, 49._dp, 49._dp,-49._dp,-70._dp, 14._dp, 56._dp, 56._dp, 14._dp,-70._dp,&
                     &  49._dp,-49._dp,-49._dp, 49._dp, 70._dp,-14._dp,-56._dp,-56._dp,-14._dp, 70._dp,&
                     &  70._dp,-70._dp,-70._dp, 70._dp,100._dp,-20._dp,-80._dp,-80._dp,-20._dp,100._dp,&
                     & -14._dp, 14._dp, 14._dp,-14._dp,-20._dp,  4._dp, 16._dp, 16._dp,  4._dp,-20._dp,&
                     & -56._dp, 56._dp, 56._dp,-56._dp,-80._dp, 16._dp, 64._dp, 64._dp, 16._dp,-80._dp,&
                     & -56._dp, 56._dp, 56._dp,-56._dp,-80._dp, 16._dp, 64._dp, 64._dp, 16._dp,-80._dp,&
                     & -14._dp, 14._dp, 14._dp,-14._dp,-20._dp,  4._dp, 16._dp, 16._dp,  4._dp,-20._dp,&
                     &  70._dp,-70._dp,-70._dp, 70._dp,100._dp,-20._dp,-80._dp,-80._dp,-20._dp,100._dp/),(/10,10/))
   app(1,:,:)=a2pp(:,:)/1225._dp
   app(2,:,:)=zero

   b2pp(:,:)=RESHAPE((/1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp,0._dp,&
                    &  0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp,1._dp/),(/10,10/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/ 49._dp, 98._dp, 98._dp,  0._dp, 30._dp, 36._dp, 27._dp, 12._dp,  0._dp,  0._dp,&
                      & 98._dp, 49._dp,  0._dp, 98._dp, 40._dp,  2._dp,  6._dp, 25._dp, 32._dp,  0._dp,&
                      & 98._dp,  0._dp, 49._dp, 98._dp,  0._dp, 32._dp, 25._dp,  6._dp,  2._dp, 40._dp,&
                      &  0._dp, 98._dp, 98._dp, 49._dp,  0._dp,  0._dp, 12._dp, 27._dp, 36._dp, 30._dp,&
                      & 30._dp, 40._dp,  0._dp,  0._dp,100._dp,120._dp, 60._dp,  0._dp,  0._dp,  0._dp,&
                      & 36._dp,  2._dp, 32._dp,  0._dp,120._dp,  4._dp, 48._dp,108._dp,  0._dp,  0._dp,&
                      & 27._dp,  6._dp, 25._dp, 12._dp, 60._dp, 48._dp, 64._dp,  0._dp,108._dp,  0._dp,&
                      & 12._dp, 25._dp,  6._dp, 27._dp,  0._dp,108._dp,  0._dp, 64._dp, 48._dp, 60._dp,&
                      &  0._dp, 32._dp,  2._dp, 36._dp,  0._dp,  0._dp,108._dp, 48._dp,  4._dp,120._dp,&
                      &  0._dp,  0._dp, 40._dp, 30._dp,  0._dp,  0._dp,  0._dp, 60._dp,120._dp,100._dp/),(/10,10/))
   bpp(1,:,:)=b2pp(:,:)/1225._dp
   write(std_out,*) "warning: this test is only valid if f4of2_sla=0"
 end if ! lcor=2
 if (lcor == 3) then
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
                       & 0.0,   1.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   1.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   1.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   1.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,&
                       & 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0/),(/14,14/))
   bpp(0,:,:)=b2pp(:,:)
   b2pp(:,:)=RESHAPE((/100.0, 120.0,  60.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     & 120.0,   4.0,  48.0, 108.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &  60.0,  48.0,  64.0,   0.0, 108.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0, 108.0,   0.0,  64.0,  48.0,  60.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0, 108.0,  48.0,   4.0, 120.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,  60.0, 120.0, 100.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,&
                     &   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/),(/14,14/))
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
 end if ! lcor=3

 do m2=1,tndim
   do m1=1,tndim
     udens(m1,m2) = sum(fk(:)*app(:,m1,m2))
     jdens(m1,m2) = sum(fk(:)*bpp(:,m1,m2))
       !write(6,*) kk,m1,m2
       !write(6,*) "--",fk(kk),aklmlmp(kk,m1,m2)
       !write(6,*) "--",fk(kk),bklmlmp(kk,m1,m2)
       !udens(m1,m2)=udens(m1,m2)+fk(kk)*app(kk,m1,m2)
       !jdens(m1,m2)=jdens(m1,m2)+fk(kk)*bpp(kk,m1,m2)
   end do ! m1
 end do ! m2
 write(message,'(2x,a,3x,14f10.4)') " Direct Interaction Matrix from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1

 write(message,'(a,2x,a,3x,14f10.4)') ch10," Exchange Interaction Matrix from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(jdens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1

 write(message,'(a,2x,a,3x,14f10.4)') ch10, " Density Density interactions from Inglis tables (in the JMJ basis) "
 call wrtout(std_out,message,'COLL')
 write(message,'(2x,4x,14(2x,i8))') (m1,m1=1,tndim,1)
 call wrtout(std_out,message,'COLL')
 do m1=1,tndim
   write(message,'(2x,i4,3x,14f10.4)') m1,(udens(m1,m2)-jdens(m1,m2),m2=1,tndim,1)
   call wrtout(std_out,message,'COLL')
 end do ! m1


 ABI_FREE(jdens)
 ABI_FREE(udens)
 ABI_FREE(app)
 ABI_FREE(bpp)
 ABI_FREE(a2pp)
 ABI_FREE(b2pp)

 end subroutine udens_inglis_hu
!!***

END MODULE m_hu
!!***
