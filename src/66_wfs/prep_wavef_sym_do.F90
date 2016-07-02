!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_wavef_sym_do
!! NAME
!! prep_wavef_sym_do
!!
!! FUNCTION
!! this routine associates waves functions by two as following
!!      E(G)  = C(G) + D(G)
!!      E(-G) = C*(G) + iD*(G)
!! the values are distributed on the processors in function of
!! the value of mpi_enreg%distribfft%tab_fftwf2_distrib( (-kg_k_gather(2,i) )
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg          = informations about mpi parallelization
!!  bandpp             = number of couple of waves functions
!!  nspinor            = number of spin
!!  ndatarecv          = number of values received by the processor and sended 
!!                       by the other processors band
!!  ndatarecv_tot      = total number of received values 
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited 
!!                       planewave coordinates
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!!  cwavef_alltoall    = planewave coefficients of wavefunction 
!!                      ( initial of the processor + sended by other processors band)
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor from each processor fft
!!  rdispls_sym        = postions of the received values by the processor from each processor fft
!!
!! OUTPUT
!!  ewavef_alltoall_sym = planewave coefficients of wavefunction 
!!                        initial of the processor + 
!!                        sended by other processors band +
!!                        sended by other processors fft  +
!!                        and compisited if bandpp >1                                       
!!  index_wavef_send    = index to send the values in blocks to the other processor fft
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

subroutine prep_wavef_sym_do(mpi_enreg,bandpp,nspinor,&
&     ndatarecv,&
&     ndatarecv_tot,ndatasend_sym,tab_proc,&
&     cwavef_alltoall,&
&     sendcounts_sym,sdispls_sym,&
&     recvcounts_sym,rdispls_sym,&
&     ewavef_alltoall_sym,&
&     index_wavef_send)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prep_wavef_sym_do'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,ndatarecv_tot,ndatasend_sym
 integer,intent(in) :: nspinor
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,allocatable,intent(out) :: index_wavef_send(:) 
 integer,intent(in) :: rdispls_sym(:),recvcounts_sym(:)
 integer,intent(in) :: sdispls_sym(:),sendcounts_sym(:)
 integer,intent(in) :: tab_proc(:)
 real(dp),intent(inout) :: cwavef_alltoall(2,ndatarecv*nspinor*bandpp)
 real(dp),pointer :: ewavef_alltoall_sym(:,:)

!Local variables-------------------------------
!scalars
 integer :: bandpp_sym,ibandpp,idatarecv,ideb_loc,idebc,idebd,idebe
 integer :: ier,ifin_loc,ifinc,ifind,ifine,iproc,jbandpp,jsendloc
 integer :: kbandpp,newspacecomm,nproc_fft
 logical :: flag_compose
!arrays
 integer,allocatable :: rdispls_sym_loc(:),recvcounts_sym_loc(:)
 integer,allocatable :: sdispls_sym_loc(:),sendcounts_sym_loc(:)
 real(dp),allocatable :: ewavef_alltoall_loc(:,:),ewavef_alltoall_send(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' prep_wavef_sym_do : enter '
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
 ABI_ALLOCATE(ewavef_alltoall_sym     ,(2,ndatarecv_tot*bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_loc     ,(2,ndatarecv    *bandpp_sym))
 ABI_ALLOCATE(ewavef_alltoall_send    ,(2,ndatasend_sym*bandpp_sym))
 ABI_ALLOCATE(index_wavef_send        ,(  ndatasend_sym*bandpp_sym))

 ABI_ALLOCATE(sendcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(sdispls_sym_loc       ,(nproc_fft))
 ABI_ALLOCATE(recvcounts_sym_loc    ,(nproc_fft))
 ABI_ALLOCATE(rdispls_sym_loc       ,(nproc_fft))
 

!Initialisation
!--------------
 ewavef_alltoall_sym(:,:) =0.
 ewavef_alltoall_loc(:,:) =0.

 sendcounts_sym_loc(:) =0
 sdispls_sym_loc(:)    =0
 recvcounts_sym_loc(:) =0
 rdispls_sym_loc(:)    =0

 index_wavef_send(:)   =0

 
!-------------------------------------------------
!We are bandpp blocks which we want to :
!associate by two      (band_sym==bandpp/2)
!or not associate by two  (band_sym==bandpp)
!
!So We'll have got bandpp_sym blocks
!So we loop on the bandpp_sym blocks
!--------------------------------------------------

 do kbandpp=1,bandpp_sym
   
!  position of the two blocks  
!  --------------------------
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


!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine) =    &
&     cwavef_alltoall(1,idebc:ifinc) &
&     - cwavef_alltoall(2,idebd:ifind)
     
     ewavef_alltoall_sym(2,idebe:ifine) =    &
&     cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)
     
!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =  &
&     cwavef_alltoall(1,idebc:ifinc) &
&     + cwavef_alltoall(2,idebd:ifind)
     
     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) =  &
&     - cwavef_alltoall(2,idebc:ifinc) &
&     + cwavef_alltoall(1,idebd:ifind)
   else

!    calcul ewavef(G)
!    ----------------
     ewavef_alltoall_sym(1,idebe:ifine)   = cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_sym(2,idebe:ifine)   = cwavef_alltoall(2,idebc:ifinc) 

!    calcul ewavef_loc(-G)
!    ---------------------
     ewavef_alltoall_loc(1,ideb_loc:ifin_loc) =   cwavef_alltoall(1,idebc:ifinc)
     ewavef_alltoall_loc(2,ideb_loc:ifin_loc) = - cwavef_alltoall(2,idebc:ifinc)
     
   end if

 end do



!------------------------------------------------------------------------
!Creation of datas blocks for each processor fft from ewavef_alltoall_loc
!to send datas by blocks with a alltoall...
!------------------------------------------------------------------------
 
!Position of the blocks 
!----------------------
 jsendloc=0
 do ibandpp=1,bandpp_sym
   do iproc=1,nproc_fft
     do idatarecv=1,ndatarecv
       if (tab_proc(idatarecv)==(iproc-1)) then
         jsendloc=jsendloc+1
         index_wavef_send(jsendloc)  = idatarecv + ndatarecv * (ibandpp-1)         
       end if
     end do
   end do
 end do

!Classment
!----------
 ewavef_alltoall_send(:,:)=ewavef_alltoall_loc(:,index_wavef_send)


!-------------------------------------------------
!Calcul of the number of received and sended datas
!-------------------------------------------------
 sendcounts_sym_loc = sendcounts_sym*2
 recvcounts_sym_loc = recvcounts_sym*2
 
!------------------------------------------
!Exchange of the datas ewavef_allto_all_loc
!------------------------------------------
 do ibandpp=1,bandpp_sym

!  ------------------------------------------------  
!  Deplacment of the sended datas because of bandpp
!  ------------------------------------------------
   sdispls_sym_loc(:) = sdispls_sym(:) + ndatasend_sym * (ibandpp-1)
   sdispls_sym_loc    = sdispls_sym_loc   *2
   
!  --------------------------------------------------
!  Deplacment of the received datas because of bandpp
!  --------------------------------------------------
   rdispls_sym_loc(:) = rdispls_sym(:) + ndatarecv_tot * (ibandpp-1)
   rdispls_sym_loc    = rdispls_sym_loc   *2

   
   call xmpi_alltoallv(&
&   ewavef_alltoall_send(:,:) ,sendcounts_sym_loc,sdispls_sym_loc,&
&   ewavef_alltoall_sym(:,:)  ,recvcounts_sym_loc,rdispls_sym_loc,&
&   newspacecomm,ier)
   
 end do

!-----------------------
!Desallocation
!-----------------------

 ABI_DEALLOCATE(sendcounts_sym_loc)
 ABI_DEALLOCATE(recvcounts_sym_loc)
 ABI_DEALLOCATE(sdispls_sym_loc)
 ABI_DEALLOCATE(rdispls_sym_loc)
 
 ABI_DEALLOCATE(ewavef_alltoall_loc)
 ABI_DEALLOCATE(ewavef_alltoall_send)

end subroutine prep_wavef_sym_do
!!***
