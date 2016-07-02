!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rec
!! NAME
!!  m_rec
!!
!! FUNCTION
!!  This module  provides some functions applied to the
!!  recursion structured datatype recursion_type.
!!  It includes also some function used to change some variables
!!  of recursion_type
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
!!
!! * Routines tagged with "@type_name" are strongly connected to the definition of the data type. 
!!   Strongly connected means that the proper functioning of the implementation relies on the 
!!   assumption that the tagged procedure is consistent with the type declaration.
!!   Every time a developer changes the structure "type_name" adding new entries, he/she has to make sure 
!!   that all the strongly connected routines are changed accordingly to accomodate the modification of the data type. 
!!   Typical examples of strongly connected routines are creation, destruction or reset methods.
!!
!! PARENTS
!!     
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rec

 use defs_basis
 use defs_rectypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 implicit none

 private ::           &
   find_maxmin_proc,  &  !--To calculate max and min pt for any cpu  
   H_D_distrib

 public ::            &
   InitRec,           &  !--Main creation method.
   Init_MetricRec,    &  !--To Initalize the inf. metric in recursion
   Init_nlpspRec,     &  !--Main creation method for non-local part.
   CleanRec,          &  !--deallocate all pointers.
   Calcnrec,          &  !--calculates the new min_nrec
   cpu_distribution      !--Regulates the work load on cpu-gpu
CONTAINS  !===========================================================
!!***

!!****f* m_rec/H_D_distrib
!! NAME
!! H_D_distrib
!!
!! FUNCTION
!! Calculate the number of point,GPU,for any proc
!!
!! INPUTS   
!!  rset<recursion_type>= recursion variables
!!  cpu (-1 if there are not gpu)
!!  nfft=nuber of point of the fine grid
!!  ngfftrec=nuber of point of one edge of the coarse grid
!!  gratio=recgratio ratio between the fine and coarse grid
!!  beta_coeff=estimated time ratio between CPU_time and GPU_time
!!  
!! OUTPUT
!!  proc_pt_dev(2,0:nproc-1) which device and how many points 
!!  that proc has to compute: proc_pt_dev(1,iproc) which device
!!  associated to proc i (-1 if none), proc_pt_dev(2,iproc) how 
!!  many points 
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE

subroutine H_D_distrib(rset,nfft,gratio,proc_pt_dev,&
     &                    beta_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'H_D_distrib'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nfft,gratio
 real(dp),intent(in) :: beta_coeff
 integer,pointer :: proc_pt_dev(:,:)  
 type(recursion_type),intent(inout) :: rset
!Local ---------------------------
 integer :: me,icpu,resto,ntot,ngpu
 integer :: n_per_cpu,n_per_gpu
 character(500) :: msg
#ifdef HAVE_GPU_CUDA
 integer,pointer :: ndev(:)
#else
 integer :: ndev(0:rset%mpi%nproc-1)
#endif
! *********************************************************************


#ifdef HAVE_GPU_CUDA
 ndev => rset%GPU%map
#else
 ndev = -1
#endif

 me = rset%mpi%me
 ntot = nfft/(gratio*gratio*gratio)
 ngpu = rset%ngpu

 !--If sequential code all points are computed by the proc 0
 if(rset%mpi%nproc ==1) then
   proc_pt_dev(1,0) = ndev(0)
   proc_pt_dev(2,0) = ntot
   return
 end if

 !--Number of points for any cpu
 n_per_cpu = int(int(ntot/(rset%mpi%nproc+ngpu*(beta_coeff-1.d0))))
 n_per_gpu = int(n_per_cpu*beta_coeff)
 !write(std_out,*)'n_per_cpu',n_per_cpu
 !write(std_out,*)'rset%GPU%map',rset%GPU%map
 do icpu=0,rset%mpi%nproc-1
   proc_pt_dev(1,icpu) = ndev(icpu)
   proc_pt_dev(2,icpu) = n_per_cpu
   if(ndev(icpu)>-1) proc_pt_dev(2,icpu) = n_per_gpu
 end do

 !--Distribute the rest 
 resto = ntot-sum(proc_pt_dev(2,:))
 icpu = 0
 !write(std_out,*)'rest',resto,ngpu
 if(resto>0) then
   if(ngpu/=0) then
     !--distribute rest only on GPU
     do while(resto/=0)
       if(proc_pt_dev(1,icpu)>-1) then
         proc_pt_dev(2,icpu) = proc_pt_dev(2,icpu)+1
         resto = resto-1
       endif
       icpu = mod(icpu+1,rset%mpi%nproc)
     enddo
   else
     !--distribute rest on all CPU
     do while(resto/=0)
       proc_pt_dev(2,icpu) = proc_pt_dev(2,icpu)+1
       resto = resto-1
       icpu = mod(icpu+1,rset%mpi%nproc)
     enddo
     return 
   endif
 endif

 !--Printing GPU and load distribution on procs
 write(msg,'(3a)')&
      & ' -Load on procs------------',ch10,& 
      & '   me  device        points'
 call wrtout(std_out,msg,'COLL')
 do icpu=0,rset%mpi%nproc-1
   write(msg,'(i5,i8,i14)') icpu,proc_pt_dev(:,icpu);
   call wrtout(std_out,msg,'COLL')
 end do

end subroutine H_D_distrib
!!***
 


!!****f* m_rec/find_maxmin_proc
!! NAME
!! find_maxmin_proc
!!
!! FUNCTION
!! To calculate max and min pt for any cpu, it is useful for 
!! recgratio!=1
!!
!! INPUTS
!! nproc = number of procs
!! me = identity of the proc
!! ngfft(3) = fine grid (corresponds to dtset%ngfft(1:3))
!! proc_pt_dev(2,0:nproc-1) which device and how many points
!! recpar%npt = number of points computed by the proc me (see side effects)
!!
!! OUTPUT
!! recpar%pt0<type(vec_int)>=Intial point for this proc in x,y,z
!! recpar%pt1<type(vec_int)>=Final point for this proc in x,y,z
!! recpar%min_pt=Intial point for this proc
!! recpar%max_pt=Final point for this proc 
!!
!! SIDE EFFECTS
!! recpar%ntranche=number of pts computed by the proc me on the fine grid. 
!!
!! 
!! So when  recgratio!=1, ntranche will  not correspond to the npt!
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE
subroutine find_maxmin_proc(recpar,nproc,me,gratio,&
  &                         ngfft,proc_pt_dev)

  
  use m_rec_tools,only  :get_pt0_pt1

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'find_maxmin_proc'
!End of the abilint section

 implicit none  

!Arguments ------------------------------------
 integer,intent(in)   :: nproc,me,gratio
 integer,intent(in)   :: ngfft(3)
 type(recparall_type),intent(inout) :: recpar
 integer,pointer :: proc_pt_dev(:,:)
!Local ---------------------------
 integer :: pointoncpu
 integer :: nfft,ntot,ii
 integer :: inf,sup
 integer :: proc_limit(0:nproc-1)
! *********************************************************************
 !  write(std_out,*)'start find_maxmin_proc'
 recpar%npt = 0
 nfft = product(ngfft)
 ntot = nfft/(gratio*gratio*gratio)
 pointoncpu = ntot/nproc

 proc_limit = (/(sum(proc_pt_dev(2,:ii)),ii=0,nproc-1)/) 

 if(gratio==1)then
   recpar%ntranche = proc_limit(me)
   if(me/=0) recpar%ntranche = recpar%ntranche-proc_limit(me-1)
 endif

 inf=0 
 if(me/=0) inf = proc_limit(me-1)
 sup = proc_limit(me)
 
 
 call get_pt0_pt1(ngfft,gratio,inf,sup,recpar)
 
 recpar%npt = sup-inf

 !write(std_out,*)'exit find_maxmin_proc'
end subroutine find_maxmin_proc
!!***


!!****f* m_rec/cpu_distribution
!! NAME
!! cpu_distribution
!!
!! FUNCTION
!! Calculate the number of point,GPU,for any proc
!!
!!
!! INPUTS
!!  ngfft(3)=nuber of point of the grid
!!  gratio=recgratio ratio between the fine and coarse grid
!!  beta_coeff=estimated time ratio between CPU_time and GPU_time
!!  calc_type=if 0 takes the possible max for nptrec (to test the
!!  completly full graphic card). 1 after test to calculate the min
!!  possible value for nptrec
!!
!! OUTPUT
!!
!! PARENTS
!!      first_rec,m_rec
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE

 subroutine cpu_distribution(gratio,rset,ngfft,beta_coeff,calc_type)

#ifdef HAVE_GPU_CUDA
 use m_hidecudarec,only   :InitRecGPU
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cpu_distribution'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------ 
 integer,intent(in)  :: gratio,calc_type
 real(dp),intent(in) :: beta_coeff
 integer,intent(in)  :: ngfft(3)
 type(recursion_type),intent(inout),target :: rset  
!Local ---------------------------
 integer :: ii,nfft,ierr
 integer,pointer :: proc_pt_dev(:,:)
 type(recparall_type),pointer :: recpar
 character(500) :: msg
! *********************************************************************
 
 ! write(std_out,*)'start cpu_distribution'

 nullify(proc_pt_dev)
 ABI_ALLOCATE(proc_pt_dev,(2,0:rset%mpi%nproc-1))

 nfft = product(ngfft) 
 call H_D_distrib(rset,nfft,gratio,proc_pt_dev,beta_coeff)

 nullify(recpar)
 if(rset%load == 0)then
   ABI_ALLOCATE(rset%par%displs,(0:rset%mpi%nproc-1))
   ABI_ALLOCATE(rset%par%vcount,(0:rset%mpi%nproc-1))
   recpar => rset%par
#if defined HAVE_GPU_CUDA
 else
   if(rset%tp==4)then
     ABI_ALLOCATE(rset%GPU%par%displs,(0:rset%mpi%nproc-1))
     ABI_ALLOCATE(rset%GPU%par%vcount,(0:rset%mpi%nproc-1))
   endif
   recpar => rset%GPU%par
#endif
 endif
 
 recpar%ntranche = nfft/(rset%mpi%nproc)!equipartitioned point

 call find_maxmin_proc(recpar,rset%mpi%nproc,&
&                      rset%mpi%me,gratio,ngfft,proc_pt_dev)

 recpar%vcount = 0
 if(rset%load==0)then
   recpar%vcount(rset%mpi%me) = recpar%ntranche
 else
   recpar%vcount(rset%mpi%me) = recpar%npt
 endif

 call xmpi_sum(recpar%vcount,rset%mpi%comm_bandfft,ierr)
   
 recpar%displs = 0
 if(rset%mpi%nproc>1) recpar%displs(1:) = (/(sum(recpar%vcount(:ii)),ii=0,rset%mpi%nproc-1)/) 

 !--INITALIZATION OF CUDA FOR RECURSION
#if defined HAVE_GPU_CUDA
 if(rset%load == 0)   rset%GPU%par = rset%par
 call InitRecGPU(rset,nfft,gratio,rset%GPU%map(rset%mpi%me),calc_type)
#else
  ierr = calc_type !only of abirule when there is not HAVE_GPU_CUDA
#endif

 
! if(rset%debug ) then
 write(msg,'(a,i7,2(2a,3i7),8(2a,i7),2(2a,3i7),(2a,e14.6))')&
   & ' me                 ',  rset%mpi%me,ch10,&
   & ' ngfft              ',  ngfft(1:3),ch10,&
   & ' ngfftrec           ',  rset%ngfftrec(1:3),ch10,&
   & ' load               ',  rset%load,ch10,&
   & ' ntranche           ',  recpar%ntranche,ch10,&
   & ' min_pt             ',  recpar%min_pt,ch10,&
   & ' max_pt             ',  recpar%max_pt,ch10,&
   & ' rset%mpi%nproc     ',  rset%mpi%nproc,ch10,&
   & ' rset%mpi%nproc_fft ',  rset%mpi%nproc_fft,ch10,&
   & ' dtset%ngfft(10)    ',  rset%ngfftrec(10),ch10,&
   & ' recpar%npt         ',  recpar%npt,ch10,&
   & ' recpar%pt0         ',  recpar%pt0%x,recpar%pt0%y,recpar%pt0%z,ch10,&
   & ' recpar%pt1         ',  recpar%pt1%x,recpar%pt1%y,recpar%pt1%z,ch10,&
   & ' grid step          ',  rset%inf%tr(1)
 call wrtout(std_out,msg,'PERS')
#if defined HAVE_GPU_CUDA
 write(msg,'(a,i7,2(2a,i7),a)')&
   & ' rset%ngp           ',  rset%ngpu,ch10,&
   & ' gpudevice          ',  rset%gpudevice,ch10,&
   & ' nptrec             ',  rset%GPU%nptrec,ch10
 call wrtout(std_out,msg,'PERS')
#endif
!  write(std_out,*)'display',recpar%displs
!  write(std_out,*)'vcount',recpar%vcount
! end if


 nullify(recpar)
 if(associated(proc_pt_dev))  then
   ABI_DEALLOCATE(proc_pt_dev)
 end if

! write(std_out,*)'exit from cpu_distribution'
end subroutine cpu_distribution
!!***


!!****f* m_rec/InitRec
!! NAME
!! InitRec
!!
!! FUNCTION
!! Initialise the rset<recursion_type>=Data type concerning recursion.
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!! mpi_ab <type(mpi_type)=MPI-parallelisation information
!! mproj=0 if psp is only local
!!
!! SIDE EFFECTS
!! All pointers set to null().
!!
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE

subroutine InitRec(dtset,mpi_ab,rset,rmet,mproj)
  
 use defs_abitypes
 use m_errors
 use m_pawfgr, only : pawfgr_nullify, indgrid

#ifdef HAVE_GPU_CUDA
 use m_gpu_detect,only    :get_topo,find_set_gpu
 use m_hidecudarec,only   :InitRecGPU_0
#include "cuda_common.h"
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'InitRec'
 use interfaces_68_recursion
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
 integer,intent(in) :: mproj
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in),target :: mpi_ab
 type(recursion_type),intent(inout) :: rset
 real(dp),intent(in) :: rmet(3,3)
! arrays
!Local ---------------------------
 integer :: ii
 real(dp) :: beta,rtrotter
#if defined HAVE_GPU_CUDA
 character(500) :: msg
#endif
! *********************************************************************
 ! @recursion_type
 ! @pawfgr_type

 !--Initialisation
 beta = one/dtset%tsmear  !--Inverse of temperature
 !--Rewriting the trotter parameter
 rtrotter  = max(half,real(dtset%recptrott,dp))

 rset%debug= (dtset%prtvol==-7)
 rset%quitrec   = 0
 rset%min_nrec  = dtset%recnrec
 rset%efermi    = dtset%recefermi !initial guess for fermie
 
 rset%nfftrec   = 0
 rset%ngfftrec  = 0        
 
 rset%tronc = .False.

 rset%mpi => mpi_ab

 !--Are all pseudo-potentials local?
 rset%nl%nlpsp = (mproj/=0) 
 !--Some initialisation concerning the metrics
 !  If non-local psps then it allocates the atoms positions
 !   on the grid
 if(rset%nl%nlpsp)  then 
  ABI_ALLOCATE(rset%inf%gcart,(3,dtset%natom))
 else 
  ABI_ALLOCATE(rset%inf%gcart,(0,0))
 end if
 rset%inf%gcart = 0

 !----------------------------------------------------------
 !--TRONCATION OF THE BOX 
 !! determines new dimensions the method is similar to the one used 
 !! in getng  (except that ecut and xboxcutmin give no constraint, 
 !!  and symmetries are not handled)   
   
 call getngrec(dtset%ngfft,rmet,rset%ngfftrec,rset%nfftrec,dtset%recrcut,0.25d0*sqrt(beta/rtrotter),rset%tronc) 
 !  1/4*sqrt(beta/trotter) for guess - should be modified     

 !------------------------------------------------------------
 !--DETERMINING WHICH POINT WILL COMPUTE THAT PROC
 !----------------------------------------------------------
 !--Paralelism using the band communicator (not used in the recursion)
 !--Distribution on procs with cuda


 rset%ngpu = 0       !--Initial guess no GPU at all
 rset%gpudevice = -1 !--Inital guess no GPU associated
 rset%load = 0       !--Inital homogeneous work load 
 rset%tp = 0         !--Initial guess 1 cpu, 0 gpu

  
#ifdef HAVE_GPU_CUDA
 !--Initialise GPU variables for recursion
 call InitRecGPU_0(rset%GPU,mpi_ab)

 !--Get the distribution of GPUs on CPUs
 call find_set_gpu(mpi_ab%nproc,mpi_ab%comm_bandfft,rset%GPU%map,rset%ngpu)

 !--Get the topology of the machine
 call get_topo(rset%mpi%nproc,rset%ngpu,rset%tp) 
 if(rset%tp>4)then
   msg = 'm_rec: number of gpu>number of cpu is not implemented'
   MSG_ERROR(msg)
 endif
!  rset%tp = 0;if(rset%mpi%nproc>1)rset%tp = 1
!  rset%ngpu = 0; rset%GPU%map=-1
 !--For the moment cuda doesnt take into account non-local psp
 if(rset%nl%nlpsp) then
   rset%tp = 0;if(rset%mpi%nproc>1)rset%tp = 1
   rset%GPU%map = -1 
 endif
#else
 if(rset%mpi%nproc>1)rset%tp = 1
#endif
 
 !--Basic initialization for recursion metric (only needed for printing)
 do ii=1,3
   rset%inf%rmet(ii,:) = rmet(ii,:)/(real(dtset%ngfft(1:3)*dtset%ngfft(ii),dp))
 end do
 rset%inf%tr(:) = sqrt((/(rset%inf%rmet(ii,ii),ii=1,3)/)) !grid step
 
 !--Compute the work loqd distribution on devices (gpu,cpu)
 call cpu_distribution(dtset%recgratio,rset,dtset%ngfft(:3),1.d0,0)

 !------------------------------------------------------------
 !--DEFINITION VARIABLE COARSE-FINE GRID  TO USE TRANSGRID--INGRID FUNCTIONS
 call pawfgr_nullify(rset%pawfgr)
 !if coarse grid is used
 if (dtset%recgratio>1) then
   !fine grid--
   rset%pawfgr%mgfft = 0
   rset%pawfgr%nfft = product(dtset%ngfft(1:3))
   rset%pawfgr%ngfft(:) = dtset%ngfft(:)
   rset%pawfgr%ngfft(9:11)=(/0,1,0/)
   rset%pawfgr%ngfft(12:13)= dtset%ngfft(2:3)
   !coarse grid--
   rset%pawfgr%mgfftc = 0
   rset%pawfgr%ngfftc(:) = rset%pawfgr%ngfft(:)
   rset%pawfgr%ngfftc(:3) = floor(real(dtset%ngfft(:3)+1,dp)/real(dtset%recgratio,dp))
   rset%pawfgr%nfftc = product(rset%pawfgr%ngfftc(1:3))

   rset%pawfgr%usefinegrid = 1
   ABI_ALLOCATE(rset%pawfgr%fintocoa,(rset%pawfgr%nfft))
   ABI_ALLOCATE(rset%pawfgr%coatofin,(rset%pawfgr%nfftc))
   call indgrid(rset%pawfgr%coatofin,rset%pawfgr%fintocoa,&
     rset%pawfgr%nfftc,rset%pawfgr%nfft,&
     rset%pawfgr%ngfftc,rset%pawfgr%ngfft)

 else
   rset%pawfgr%mgfft = 0
   rset%pawfgr%ngfft = 0
   rset%pawfgr%mgfftc = 0

   rset%pawfgr%usefinegrid = 0
 end if


end subroutine InitRec
!!***

!!****f* m_rec/Init_MetricRec
!! NAME
!! Init_MetricRec
!!
!! FUNCTION
!! Initialise the rset<recursion_type>=Data type concerning recursion.
!! In particular, the information on the infinitesimal metric. 
!! Also other variable are initialized  
!!
!! INPUTS
!! rmet: metrics
!! ucvol=unit cell volume in bohr**3.
!! ngfft(1:3)=fine grid used in recursion
!! rprimd=Real space PRIMitive translations, Dimensional 
!! xred=vectors (X) of atom positions in reduced coordinates 
!! natom=number of atoms
!! debug=debug variable
!!
!! OUTPUT
!! metrec <type(metricrec_type)>= infinitesimal metrics used in recursion
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE
subroutine Init_MetricRec(metrec,nlpsp,rmet,ucvol,rprimd,xred,ngfft,&
&                        natom,debug)
 
 use m_per_cond,only     :  per_cond

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Init_MetricRec'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer, intent(in) ::natom
 real(dp), intent(in) :: ucvol
 logical,intent(in) ::nlpsp,debug
 type(metricrec_type),intent(inout) :: metrec
!arrays
 integer,intent(in) :: ngfft(3)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),xred(3,natom)

!Local ---------------------------
 integer :: ii
 real(dp) :: xcart(3,natom)
 character(500) :: msg
! *********************************************************************

 !--Intialisation of variables concerning the infinitesimal metric  
 do ii=1,3
   metrec%rmet(ii,:) = rmet(ii,:)/(real(ngfft(1:3)*ngfft(ii),dp))
 end do
 metrec%ucvol = ucvol/real(product(ngfft(1:3)),dp)
 metrec%tr(:) = sqrt((/(metrec%rmet(ii,ii),ii=1,3)/)) !grid step

 !--Initialisation of others variables
 !--In non-loc-psp case: calculate the position of ions and conversion factor
 if(nlpsp) then
   do ii = 1,natom
     xcart(:,ii) = matmul(rprimd(:,:),xred(:,ii))
   end do
   metrec%gcart(:,:) = per_cond(natom,xcart,ngfft(1:3),metrec%tr(:))
   if(debug) then
     do ii=1,natom
       write (msg,'(a,3f8.2)')'xcart=',xcart(:,ii)
       call wrtout(std_out,msg,'COLL')
       write (msg,'(a,3i4)')'gcart=',metrec%gcart(:,ii)
       call wrtout(std_out,msg,'COLL')
     end do
   end if
 end if



end subroutine Init_MetricRec
!!***          


!!****f* m_rec/Init_nlpspRec
!! NAME
!! Init_nlpspRec
!!
!! FUNCTION
!! Initialise the rset<recursion_type>=Data type concerning recursion.
!! In particular, the non-local part of pseudo-potential. 
!!
!! INPUTS
!! tempe=temperature
!! psps <type(pseudopotential_type)>=variables related to pseudo-potentials
!! metrec <type(metricrec_type)>=infinitesimal metrics used in recursion
!! ngfftrec(18)=Number of Grid points for Fast Fourier Transform for
!!   Recursion (truncated box, if different from ngfft)
!! debug=debug variable
!!
!! SIDE EFFECTS
!! nlrec <type(nlpsprec_type)>=pseudo-potentials informations for recursion
!!
!! PARENTS
!!      first_rec
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE
subroutine Init_nlpspRec(tempe,psps,nlrec,metrec,ngfftrec,debug)

 use defs_datatypes
 use defs_abitypes

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Init_nlpspRec'
 use interfaces_14_hidewrite
 use interfaces_68_recursion
!End of the abilint section

 implicit none
!Arguments ------------------------------------
! scalars
 logical,intent(in) :: debug
 real(dp), intent(in) :: tempe 
 type(pseudopotential_type),intent(in) ::psps
 type(metricrec_type),intent(inout) :: metrec
 type(nlpsprec_type),intent(inout) :: nlrec
! arrays
 integer,intent(in) :: ngfftrec(18)
!Local ---------------------------
 integer :: ii,jj
 character(500) :: msg
! *********************************************************************
 !!--Routine for the calcul of the non-local pseudo 
 if(any(psps%pspcod/=3) .and. nlrec%nlpsp ) then
   msg = "The non-local part of psp is used in Recursion only for hgh-psp"
   MSG_WARNING(msg)
   nlrec%nlpsp = .False.
   if (allocated(metrec%gcart))  then
     ABI_DEALLOCATE(metrec%gcart)
   end if
 end if

 if(any(psps%pspcod==3) .and.  nlrec%nlpsp) then
 
  nlrec%nlpsp = .True.
  nlrec%npsp  = psps%npsp
  nlrec%lmnmax = count(psps%indlmn(3,:,psps%npsp)/=0)
  ABI_ALLOCATE(nlrec%mat_exp_psp_nl,(3,3,psps%mpsang,psps%npsp))
  ABI_ALLOCATE(nlrec%eival,(3,psps%mpsang,psps%npsp))
  ABI_ALLOCATE(nlrec%eivec,(3,3,psps%mpsang,psps%npsp))
  ABI_ALLOCATE(nlrec%pspinfo,(psps%mpsang,psps%npsp))
  ABI_ALLOCATE(nlrec%radii,(psps%mpsang,psps%npsp))
  ABI_ALLOCATE(nlrec%indlmn,(6,nlrec%lmnmax,psps%npsp))
  nlrec%indlmn(:,:,:) = psps%indlmn(:,:nlrec%lmnmax,:)
  nlrec%mat_exp_psp_nl(:,:,:,:) = zero
  nlrec%eivec(:,:,:,:) = zero
  nlrec%eival(:,:,:) = zero
  nlrec%radii(:,:) = zero
  nlrec%pspinfo(:,:) = 0

  !--Get the exponential of the strength times the projectors overlap 
  !  of the non-local part of psp(hgh):
  !  h_ij=strength; g_ij=ovelap => (exp(-h.g/temp/4p)-Identity).g^(-1)
  !  And the diagonalisation of the projectors and associated eigenvectors
  call pspnl_hgh_rec(psps,tempe,nlrec,debug)

  if(debug)then
   do jj=1,psps%npsp
    write(msg,'(a)')' Exponential matrices:'
    call wrtout(std_out,msg,'COLL')
    do ii=1,psps%mpsang
     write(msg,'(a,i2,a,3f15.10,a,3f15.10,a,3f15.10)')&
       &   'angular moment',ii-1,ch10,&
       &                    nlrec%mat_exp_psp_nl(1,:,ii,jj),ch10,&
       &                    nlrec%mat_exp_psp_nl(2,:,ii,jj),ch10,&
       &                    nlrec%mat_exp_psp_nl(3,:,ii,jj)
     call wrtout(std_out,msg,'COLL')
    end do
   end do
  end if
   
  !--Now it calculates the matrix of the exp(V_NL)
  call pspnl_operat_rec(nlrec,metrec,ngfftrec,debug)

 else !--Only local pseudo potentials
  nlrec%nlpsp = .False.
  nlrec%npsp  = psps%npsp
  ABI_ALLOCATE(nlrec%mat_exp_psp_nl,(0,0,0,0))
  ABI_ALLOCATE(nlrec%pspinfo,(0,0))
  ABI_ALLOCATE(nlrec%radii,(0,0))
  ABI_ALLOCATE(nlrec%indlmn,(0,0,0))
  ABI_ALLOCATE(nlrec%projec,(0,0,0))
 endif
end subroutine Init_nlpspRec
!!***



!!****f* m_rec/CleanRec
!! NAME
!! CleanRec
!!
!! FUNCTION
!! Deallocate the pointers of rset<recursion_type>=Data type concerning recursion.
!!
!! INPUTS
!! rset<recursion_type>=Data type concerning recursion
!!
!! SIDE EFFECTS
!! All pointers are deallocated.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE
subroutine CleanRec(rset)

 use m_pawfgr, only : pawfgr_destroy
#ifdef HAVE_GPU_CUDA
 use m_hidecudarec,only     : CleanRecGPU
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'CleanRec'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
 type(recursion_type),intent(inout) :: rset
! arrays
! *********************************************************************
 
 ! @recursion_type

 if (allocated(rset%ZT_p))  then
   ABI_DEALLOCATE(rset%ZT_p)
 end if
 if (allocated(rset%par%displs))  then
   ABI_DEALLOCATE(rset%par%displs)
 end if
 if (allocated(rset%par%vcount))  then
   ABI_DEALLOCATE(rset%par%vcount)
 end if
 if (allocated(rset%nl%mat_exp_psp_nl))  then
   ABI_DEALLOCATE(rset%nl%mat_exp_psp_nl)
 end if
 if (allocated(rset%nl%eival))  then
   ABI_DEALLOCATE(rset%nl%eival)
 end if
 if (allocated(rset%nl%eivec))  then
   ABI_DEALLOCATE(rset%nl%eivec)
 end if
 if (allocated(rset%nl%pspinfo))  then
   ABI_DEALLOCATE(rset%nl%pspinfo)
 end if
 if (allocated(rset%nl%radii))  then
   ABI_DEALLOCATE(rset%nl%radii)
 end if
 if (allocated(rset%nl%indlmn))  then
   ABI_DEALLOCATE(rset%nl%indlmn)
 end if
 if (allocated(rset%nl%projec))  then
   ABI_DEALLOCATE(rset%nl%projec)
 end if
 if (allocated(rset%inf%gcart))  then
   ABI_DEALLOCATE(rset%inf%gcart)
 end if

 call pawfgr_destroy(rset%pawfgr)

 ! No is needed deallocate rset%mpi: it is a copy of mpi_enreg which
 ! pointers are deallocated in gstate 

#ifdef HAVE_GPU_CUDA
  call CleanRecGPU(rset%GPU,rset%load)
#endif

end subroutine CleanRec
!!***

!!****f* m_rec/Calcnrec
!! NAME
!! Calcnrec
!!
!! FUNCTION
!! Calculate the new min_nrec.
!!
!! INPUTS
!! rset<recursion_type>=Data type concerning recursion
!! b2(:,:) recursion coefficients
!!
!! OUTPUT
!! rset%min_nrec is changed
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      wrtout,xmpi_max
!!
!! SOURCE

subroutine Calcnrec(rset,b2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'Calcnrec'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

 !Arguments ------------------------------------
 ! scalars
 type(recursion_type),intent(inout) :: rset
 ! arrays
 real(dp), intent(in):: b2(0:rset%min_nrec,1:rset%par%ntranche)
 !Local ----------------------------------------
 ! scalars
 integer :: kk,ii,jj,ierr,loc_nrec
 character(len=500) :: msg
 ! *********************************************************************
 ! @recursion_type
 ! @pawfgr_type
 kk = 1
 loc_nrec = rset%min_nrec
 do ii=1,rset%par%ntranche
  !--Use to lbound because b2 passed as argument
  !  doesn't have the same bounds as in the calling 
  !  subroutine, the +1 because b2(lbound,ii)=1.
  jj = lbound(b2,dim=1)+1
  do while (b2(jj,ii)>tol10 .and.  jj<=rset%min_nrec-1)
   jj = jj+1
   kk = max(jj,kk)   
  end do
 enddo
 call xmpi_max(kk,rset%min_nrec,rset%mpi%comm_bandfft,ierr)
 rset%min_nrec = rset%min_nrec+1-lbound(b2,dim=1)

 write(msg,'(a,i2.2,a9,i2.2)') ' -- nrec adjustement   nrec=',loc_nrec,' => nrec=',rset%min_nrec
 call wrtout(std_out,msg,'COLL')
 write(msg,'(51a)')' ',('-',ii=1,50)
 call wrtout(std_out,msg,'COLL')
end subroutine Calcnrec
!!***


end module m_rec
!!***
