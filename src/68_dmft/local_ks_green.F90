!{\src2tex{textfont=tt}}
!!****f* ABINIT/local_ks_green
!! NAME
!! local_ks_green
!!
!! FUNCTION
!! Compute the sum over k-point of ks green function.
!! do the fourier transformation and print it
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc
!!  istep    =  step of iteration for LDA.
!!  lda_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      fourier_fct,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine local_ks_green(green,paw_dmft,prtopt)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_crystal, only : crystal_t
 use m_green, only : green_type,fourier_fct
 use m_paw_dmft, only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'local_ks_green'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type), intent(in) :: green
 type(paw_dmft_type), intent(in)  :: paw_dmft
 integer, intent(in) :: prtopt

!Local variables ------------------------------
 character(len=500) :: message
 integer :: iband,ifreq,ikpt,isppol,itau,lsub,ltau,mbandc,nkpt,nsppol
 character(len=1) :: tag_is
 character(len=fnlen) :: tmpfil
 integer,allocatable :: unitgreenlocks_arr(:)
 real(dp) :: beta
 real(dp), allocatable :: tau(:)
 complex(dpc), allocatable :: loc_ks(:,:,:)
 complex(dpc), allocatable :: loc_ks_tau(:,:,:),fw(:),ft(:)
!scalars
!************************************************************************
 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 ltau=128
 ABI_ALLOCATE(tau,(ltau))
 do itau=1,ltau
   tau(itau)=float(itau-1)/float(ltau)/paw_dmft%temp
 end do
 beta=one/paw_dmft%temp

!Only imaginary frequencies here
 if(green%w_type=="real") then
   message = ' compute_energy not implemented for real frequency'
   MSG_BUG(message)
 end if

!=========================================
!Compute local band ks green function
! should be computed in compute_green: it would be less costly in memory.
!=========================================
 ABI_ALLOCATE(loc_ks,(nsppol,mbandc,paw_dmft%dmft_nwlo))
 if(green%oper(1)%has_operks==1) then
   loc_ks(:,:,:)=czero
   do isppol=1,nsppol
     do iband=1,mbandc
       do ifreq=1,paw_dmft%dmft_nwlo
         do ikpt=1,nkpt
           loc_ks(isppol,iband,ifreq)=loc_ks(isppol,iband,ifreq)+  &
&           green%oper(ifreq)%ks(isppol,ikpt,iband,iband)*paw_dmft%wtk(ikpt)
         end do
       end do
     end do
   end do
 else
   message = ' green fct is not computed in ks space'
   MSG_BUG(message)
 end if

!=========================================
!Compute fourier transformation 
!=========================================

 ABI_ALLOCATE(loc_ks_tau,(nsppol,mbandc,ltau))
 ABI_ALLOCATE(fw,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(ft,(ltau))
 loc_ks_tau(:,:,:)=czero
 do isppol=1,nsppol
   do iband=1,mbandc
     do ifreq=1,paw_dmft%dmft_nwlo
       fw(ifreq)=loc_ks(isppol,iband,ifreq)
     end do
     call fourier_fct(fw,ft,.true.,ltau,-1,paw_dmft) ! inverse fourier
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=ft(itau)
     end do
   end do
 end do
 ABI_DEALLOCATE(fw)
 ABI_DEALLOCATE(ft)
 do isppol=1,nsppol
   do iband=1,mbandc
     do itau=1,ltau
       loc_ks_tau(isppol,iband,itau)=(loc_ks_tau(isppol,iband,itau)+conjg(loc_ks_tau(isppol,iband,itau)))/two
     end do
   end do
 end do

!=========================================
!Print out ksloc green function
!=========================================
 if(abs(prtopt)==1) then
   ABI_ALLOCATE(unitgreenlocks_arr,(nsppol))
   do isppol=1,nsppol
     write(tag_is,'(i1)')isppol
     tmpfil = trim(paw_dmft%filapp)//'Gtau_locks_isppol'//tag_is
     write(message,'(3a)') ch10," == Print green function on file ",tmpfil
     call wrtout(std_out,message,'COLL')
     unitgreenlocks_arr(isppol)=500+isppol-1
     open (unit=unitgreenlocks_arr(isppol),file=trim(tmpfil),status='unknown',form='formatted')
     rewind(unitgreenlocks_arr(isppol))
     write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitgreenlocks_arr(isppol)
     write(message,'(a,a)') ch10,"# New record : First 40 bands"
     call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
     do lsub=1,mbandc/40+1
       do itau=1, ltau
         write(message,'(2x,50(e10.3,2x))') tau(itau), &
&         (real(loc_ks_tau(isppol,iband,itau)),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end do
       write(message,'(2x,50(e10.3,2x))') beta, &
&       ((-one-real(loc_ks_tau(isppol,iband,1))),iband=40*(lsub-1)+1,min(40*lsub,mbandc))
       call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       if(40*lsub<mbandc) then
         write(message,'(a,a,i5,a,i5)')    &
&         ch10,"# Same record, Following bands : From ",    &
&         40*(lsub),"  to ",min(40*(lsub+1),mbandc)
         call wrtout(unitgreenlocks_arr(isppol),message,'COLL')
       end if
     end do
!    call flush(unitgreenlocks_arr(isppol))
   end do
   ABI_DEALLOCATE(unitgreenlocks_arr)
 end if

!Deallocations
 ABI_DEALLOCATE(loc_ks)
 ABI_DEALLOCATE(loc_ks_tau)
 ABI_DEALLOCATE(tau)

end subroutine local_ks_green
!!***
