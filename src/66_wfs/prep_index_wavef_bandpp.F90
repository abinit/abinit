!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_index_wavef_bandpp
!! NAME
!! prep_index_wavef_bandpp
!!
!! FUNCTION
!! this routine sorts the waves functions by bandpp and by processors 
!! after the alltoall
!!
!! COPYRIGHT
!!
!! INPUTS
!!  nproc_band = number of processors below the band
!!  bandpp     = number of groups of couple of waves functions
!!  nspinor    = number of spin
!!  ndatarecv  = total number of values received by the processor and sended 
!!               by the other processors band
!!  recvcounts = number of values sended by each processor band and received 
!!               by the processor
!!  rdispls    = positions of the values received by the processor and  
!!               sended by each processor band
!!
!! OUTPUT
!!  index_wavef_band = position of the sorted values
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      chebfi,prep_fourwf,prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&
                             nspinor,ndatarecv,&
                             recvcounts,rdispls,&
                             index_wavef_band)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_index_wavef_bandpp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,allocatable,intent(out) :: index_wavef_band(:)

!Local variables-------------------------------
!scalars
 integer :: delta,idebc,idebe,ifinc,ifine,iindex,iproc,kbandpp,nb

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_index_wavef_banpp : enter '
!write(std_out,*) 'ndatarecv = ', ndatarecv
!write(std_out,*) 'rdispls(:) = ', rdispls(:) 
!write(std_out,*) 'recvcounts(:) = ', recvcounts(:) 
!ENDDEBUG


!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(index_wavef_band ,(bandpp*nspinor*ndatarecv))
 index_wavef_band(:)   =0
 
!---------------------------------------------
!Calcul : loops on bandpp and processors band
!---------------------------------------------
 nb = sum(recvcounts(1:nproc_band))
 do kbandpp=1,bandpp
   
   do iproc=1,nproc_band
     
     idebe = (rdispls(iproc) + 1)  + (kbandpp-1) * ndatarecv*nspinor
     ifine = idebe + recvcounts(iproc) -1

     if (iproc==1) then
       idebc =   (kbandpp-1)* recvcounts(iproc)*nspinor + 1  
     else
       idebc = (bandpp)  * sum(recvcounts(1:iproc-1))*nspinor &
       + (kbandpp-1)* recvcounts(iproc)*nspinor &
       + 1 
     end if
     ifinc = idebc + recvcounts(iproc) -1
     index_wavef_band(idebe:ifine) = (/( iindex,iindex=idebc,ifinc)/)
     delta=ifine-idebe
     if (nspinor==2) then
       index_wavef_band(idebe+nb :idebe+nb +delta)=(/( iindex,iindex=ifinc+1,ifinc+1+delta)/)
     end if
   end do
 end do

end subroutine prep_index_wavef_bandpp
!!***
