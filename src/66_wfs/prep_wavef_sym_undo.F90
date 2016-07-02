!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_wavef_sym_undo
!! NAME
!! prep_wavef_sym_undo
!!
!! FUNCTION
!! this routine dissociates each wave function in two waves functions as following
!!      C(G) =   ( E*(-G) + E(G))/2
!!      D(G) = i*( E*(-G) - E(G))/2
!! the values are redistributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg          = informations about mpi parallelization
!!  bandpp             = number of groups of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended 
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values 
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited 
!!                       planewave coordinates
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor to each processor fft
!!!  rdispls_sym        = postions of the received values by the processor to each processor fft
!!
!!  gwavef_alltoall_sym = planewave coefficients of wavefunction 
!!                        initial of the processor + 
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and composited if bandpp >1                                       
!!  index_wavef_send    = index to send the values by block to the other processor fft
!!
!! OUTPUT
!!  gwavef_alltoall     = planewave coefficients of wavefunction 
!!                        ( for of the processor + to send to other processors band)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_fourwf,prep_getghc
!!
!! CHILDREN
!!      xmpi_alltoallv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prep_wavef_sym_undo(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,idatarecv0,&
&     gwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     gwavef_alltoall_sym,&
&     index_wavef_send)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_wavef_sym_undo'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,idatarecv0,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: index_wavef_send(:),rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 real(dp),intent(inout) :: gwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),intent(inout) :: gwavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,ideb_loc,idebc,idebd
 integer :: idebe,ier,ifin_loc,ifinc,ifind,ifine
 integer :: jbandpp,kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: gwavef_alltoall_loc(:,:),gwavef_alltoall_rcv(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_undo : enter '
!ENDDEBUG


!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 newspacecomm = mpi_enreg%comm_fft

 if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
   flag_compose = .TRUE.
 else
   bandpp_sym   = bandpp
   flag_compose = .FALSE.
 end if

!---------------------------------------------
!Allocation
!---------------------------------------------
 ABI_ALLOCATE(gwavef_alltoall_loc     ,(2,ndatarecv     *bandpp_sym))
 ABI_ALLOCATE(gwavef_alltoall_rcv     ,(2,ndatasend_sym *bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))


!---------------------------------------------
!Initialisation
!---------------------------------------------
 gwavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0


!-------------------------------------------------
!Calcul of number of the sended and received datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2

!----------------------------------------------------
!Sending of the values
!----------------------------------------------------
 do ibandpp = 1,bandpp_sym

!  -------------------------------------------------  
!  Deplacment of the sended values because of bandpp
!  -------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2
   
!  ---------------------------------------------------
!  Deplacment of the received values because of bandpp
!  ---------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2

   
   call xmpi_alltoallv(&
&   gwavef_alltoall_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,&
&   gwavef_alltoall_rcv(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   newspacecomm,ier)

 end do


!----------------------
!Dispatching the blocks
!----------------------
 gwavef_alltoall_loc(:,index_wavef_send(:)) = gwavef_alltoall_rcv(:,:)

!----------------------
!Case  -kg = [0 0 0]
!----------------------
 if (idatarecv0/=-1) then
   do kbandpp=1,bandpp_sym
     gwavef_alltoall_loc(:,(kbandpp-1)*ndatarecv     + idatarecv0)= &
     gwavef_alltoall_sym(:,(kbandpp-1)*ndatarecv_tot + idatarecv0)
   end do
 end if

!---------------------------------------------------
!Build of hwavef_alltoall
!
!We have got :
!bandpp_sym blocks to dissociate 
!or   bandpp_sym blokcs to not dissociate
!--------------------------------------------------
 do kbandpp=1,bandpp_sym
   
!  position of the 2 blocks  
!  ----------------------------------
   ibandpp = (kbandpp-1) * 2
   jbandpp =  ibandpp    + 1
   
   idebe = (kbandpp-1) * ndatarecv_tot + 1
   ifine = idebe       + ndatarecv     - 1

   idebc = ibandpp * ndatarecv     + 1
   ifinc = idebc   + ndatarecv     - 1
   
   idebd = jbandpp * ndatarecv     + 1
   ifind = idebd   + ndatarecv     - 1
   
   ideb_loc = (kbandpp-1) * ndatarecv  + 1
   ifin_loc = ideb_loc    + ndatarecv  - 1


   if (flag_compose) then

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc) 
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)

!    calcul dwavef(G)
!    ------------------
     gwavef_alltoall(1,idebd:ifind) =   gwavef_alltoall_sym(2,idebe:ifine) &
&     + gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
     gwavef_alltoall(2,idebd:ifind) = - gwavef_alltoall_sym(1,idebe:ifine) &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc)
   else 

!    calcul cwavef(G)
!    ----------------
     gwavef_alltoall(1,idebc:ifinc) =   gwavef_alltoall_sym(1,idebe:ifine)  &
&     + gwavef_alltoall_loc(1,ideb_loc:ifin_loc) 
     gwavef_alltoall(2,idebc:ifinc) =   gwavef_alltoall_sym(2,idebe:ifine)  &
&     - gwavef_alltoall_loc(2,ideb_loc:ifin_loc)
   end if

 end do

!We divise by two
 gwavef_alltoall(:,:)    = gwavef_alltoall(:,:)/2

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)
 
 ABI_DEALLOCATE(gwavef_alltoall_loc)
 ABI_DEALLOCATE(gwavef_alltoall_rcv)

end subroutine prep_wavef_sym_undo
!!***
