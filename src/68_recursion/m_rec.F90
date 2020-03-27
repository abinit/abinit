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
!! Copyright (C) 2002-2020 ABINIT group (MMancini)
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
 use m_abicore
 use m_errors
 use m_xmpi
 use m_sort
 use m_dtset

 use defs_datatypes,    only : pseudopotential_type
 use defs_abitypes,     only : mpi_type
 use m_exp_mat,         only : exp_mat
 use m_numeric_tools,   only : set2unit
 use m_special_funcs,   only : gamma_function
 use m_pawfgr,          only : pawfgr_nullify, indgrid, pawfgr_destroy
 use m_paw_sphharm,     only : initylmr
 use m_time,            only : timab
 use m_rec_tools,       only : get_pt0_pt1
 use m_per_cond,        only : per_cond
#ifdef HAVE_GPU_CUDA
 use m_hidecudarec,     only : InitRecGPU, CleanRecGPU
#endif


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

subroutine H_D_distrib(rset,nfft,gratio,proc_pt_dev,beta_coeff)

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

subroutine find_maxmin_proc(recpar,nproc,me,gratio,ngfft,proc_pt_dev)

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
 if(rset%mpi%nproc>1) recpar%displs(1:) = (/(sum(recpar%vcount(:ii)),ii=0,rset%mpi%nproc-2)/)

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

#ifdef HAVE_GPU_CUDA
 use m_gpu_detect,only    :get_topo,find_set_gpu
 use m_hidecudarec,only   :InitRecGPU_0
#include "cuda_common.h"
#endif
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

subroutine Init_MetricRec(metrec,nlpsp,rmet,ucvol,rprimd,xred,ngfft,natom,debug)

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

!!****f* ABINIT/getngrec
!! NAME
!! getngrec
!!
!! FUNCTION
!! This routine computes the fft box for the recursion method, accordingly to the troncation radius.
!! It is quite similar to getng, but :
!!     - there is no xboxcut and ecut consistency
!!     - ngfft (the initial fft box) is the maximum fft box
!!
!! INPUTS
!!  ngfft(18)=non truncated fft box
!!  mgfft=maximum of ngfft(1:3)
!!  inf_rmet=define the infinitesimal metric : rprimd*(transpose(rprimd)), divided by the number of discretisation point
!!  recrcut=truncating
!!  delta=to obtain radius of truncation
!!
!! OUTPUT
!!  ngfftrec=truncated fft box
!!  nfftrec= truncated nfft
!!  tronc=True if truncation is made
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      sort_int,timab
!!
!! SOURCE

subroutine getngrec(ngfft,rmet,ngfftrec,nfftrec,recrcut,delta,tronc)

implicit none

!Arguments -------------------------------
!scalars
real(dp),intent(in) :: recrcut,delta
integer,intent(out) :: nfftrec
logical,intent(out) :: tronc
!arrays
integer,intent(in)  :: ngfft(18)
real(dp),intent(in) :: rmet(3,3)
integer,intent(out) :: ngfftrec(18)

!Local variables-------------------------------
!scalars
integer :: ii,iimin,index,jj,jjmin,kk,kkmin,largest_ngfftrec,maxpow11,maxpow2
integer :: maxpow3,maxpow5,maxpow7,mmsrch,plane
real(dp) :: dsm,dsp,dsqmin,rtroncat
!arrays
integer :: get_ngfftrec(3),maxsrch(3),minsrch(3)
integer,allocatable :: iperm(:),srch(:)
real(dp) :: tsec(2)
real(dp) :: inf_rmet(3,3)
! *************************************************************************

 call timab(602,1,tsec)

 ngfftrec(:) = ngfft(:)

 if (recrcut>tol14) then  !default value dtset%recrcut = zero means no troncation
   rtroncat = recrcut+delta
   get_ngfftrec(:)=1
   plane = 1

   do ii=1,3
     inf_rmet(ii,:) = rmet(ii,:)/(real(ngfft(1:3)*ngfft(ii),dp))
   end do


!  minimum value of ngfftrec
   do ii = 1,3
     ngfftrec(ii)=floor(2*rtroncat/sqrt(inf_rmet(ii,ii)))+1  !minimum value
     if(ngfftrec(ii)>=ngfft(ii))then
       ngfftrec(ii)=ngfft(ii)
       get_ngfftrec(ii)=0
     end if
   end do


   if(sum(get_ngfftrec)/=0)then
     largest_ngfftrec=maxval(ngfft(1:3))
     maxpow2=int(log(largest_ngfftrec+0.5d0)/log(two))
     maxpow3=int(log(largest_ngfftrec+0.5d0)/log(three))
     maxpow5=int(log(largest_ngfftrec+0.5d0)/log(five))
     maxpow7=0
     maxpow11=0
     mmsrch=(maxpow2+1)*(maxpow3+1)*(maxpow5+1)*(maxpow7+1)*(maxpow11+1)
     ABI_ALLOCATE(srch,(mmsrch))
     ABI_ALLOCATE(iperm,(mmsrch))
!    Factors of 2
     srch(1)=1
     do ii=1,maxpow2
       srch(ii+1)=srch(ii)*2
     end do
!    Factors of 3
     index=maxpow2+1
     if(maxpow3>0)then
       do ii=1,maxpow3
         srch(1+ii*index:(ii+1)*index)=3*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 5
     index=(maxpow3+1)*index
     if(maxpow5>0)then
       do ii=1,maxpow5
         srch(1+ii*index:(ii+1)*index)=5*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 7
     index=(maxpow5+1)*index
     if(maxpow7>0)then
       do ii=1,maxpow7
         srch(1+ii*index:(ii+1)*index)=7*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    Factors of 11
     index=(maxpow7+1)*index
     if(maxpow11>0)then
       do ii=1,maxpow11
         srch(1+ii*index:(ii+1)*index)=11*srch(1+(ii-1)*index:ii*index)
       end do
     end if
!    srch is the set of allowed ngfftrec values

     call sort_int(mmsrch,srch,iperm)
     ABI_DEALLOCATE(iperm)

     do ii=1,3
       if(get_ngfftrec(ii)==1)then
         do jj=1,mmsrch
           if(srch(jj)>=ngfftrec(ii))then
             minsrch(ii)=jj
             ngfftrec(ii)=srch(jj)
             exit
           end if
         end do
         do jj=minsrch(ii),mmsrch
           if(srch(jj)>ngfft(ii))then
!            since ngfftrec(ii)<ngfft(ii) for get_ngfftrec(ii)==1,
!            and srch(mmsrch)maxval(ngfft(1:3)),
!            that will appens in the range minsrch(ii),mmsrch
             maxsrch(ii)=jj-1
             exit
           end if
         end do
       end if
!      since ngfft(ii) is in srch, we have here srch(maxsrch(ii))=ngfft(ii)
!      minsrch(ii), maxsrch(ii) is the range of index of srch in which we can
!      search ngfftrec(ii)

       if(ngfftrec(ii)>=ngfft(ii))then
         ngfftrec(ii)=ngfft(ii)
         get_ngfftrec(ii)=0
       end if
     end do
   end if

!  verify that the entiere truncation sphere is in the fft box ;
!  but only in the dimension in which we do not consider the entiere fft box
   do while(sum(get_ngfftrec)/=0)  !again...

!    determining the minimum distance between 0 and the boundary
!    of the fft box
!    quite similar to the subroutine "bound", but only over the plane which
!    are not the whole fft box
     dsqmin=dsq_rec(ngfftrec(1)/2,-ngfftrec(2)/2,-ngfftrec(3)/2,inf_rmet)+0.01d0

     if(get_ngfftrec(1)/=0)then
!      look at +/- g1 planes:
       do jj=-ngfftrec(2)/2,ngfftrec(2)/2
         do kk=-ngfftrec(3)/2,ngfftrec(3)/2
           dsp = dsq_rec(ngfftrec(1)/2, jj, kk,inf_rmet)
           dsm = dsq_rec( - ngfftrec(1)/2, jj, kk,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ngfftrec(1)/2
             jjmin = jj
             kkmin = kk
             plane=1
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin =  - ngfftrec(1)/2
             jjmin = jj
             kkmin = kk
             plane=1
           end if
         end do
       end do
     end if

     if(get_ngfftrec(2)/=0)then
!      +/- g2 planes:
       do ii=-ngfftrec(1)/2,ngfftrec(1)/2
         do kk=-ngfftrec(3)/2,ngfftrec(3)/2
           dsp = dsq_rec(ii,ngfftrec(2)/2,kk,inf_rmet)
           dsm = dsq_rec(ii,-ngfftrec(2)/2,kk,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ii
             jjmin = ngfftrec(2)/2
             kkmin = kk
             plane=2
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin = ii
             jjmin =  - ngfftrec(2)/2
             kkmin = kk
             plane=2
           end if
         end do
       end do
     end if

     if(get_ngfftrec(3)/=0)then
!      +/- g3 planes:
       do ii=-ngfftrec(1)/2,ngfftrec(1)/2
         do jj=-ngfftrec(2)/2,ngfftrec(2)/2
           dsp = dsq_rec(ii,jj,ngfftrec(3)/2,inf_rmet)
           dsm = dsq_rec(ii,jj,-ngfftrec(3)/2,inf_rmet)
           if (dsp<dsqmin) then
             dsqmin = dsp
             iimin = ii
             jjmin = jj
             kkmin = ngfftrec(3)/2
             plane=3
           end if
           if (dsm<dsqmin) then
             dsqmin = dsm
             iimin = ii
             jjmin = jj
             kkmin =  - ngfftrec(3)/2
             plane=3
           end if
         end do
       end do
     end if

     if(dsqmin>=rtroncat)then
       get_ngfftrec=0
       exit
     end if

!    Fix nearest boundary
     do ii=minsrch(plane),maxsrch(plane)
       if (srch(ii)>=ngfftrec(plane)) then
!        redefine ngfft(plane) to next higher choice
         ngfftrec(plane)=srch(ii+1)
!        verify if we cover the whole box
         if(ngfftrec(plane)>=ngfft(plane))then
           ngfftrec(plane)=ngfft(plane)
           get_ngfftrec(plane)=0
         end if
!        Exit the loop over ii
         exit
       end if
     end do

   end do

   if (allocated(srch)) then
     ABI_DEALLOCATE(srch)
   end if

!  if(mod(ngfftrec(1),16)/=0) then
!  ngfftrec(1) = ngfftrec(1)+(16-mod(ngfftrec(1),16))
!  ngfftrec(2:3) = ngfftrec(1)
!  endif

   ngfftrec(4)=2*(ngfftrec(1)/2)+1
   ngfftrec(5)=2*(ngfftrec(2)/2)+1
   ngfftrec(6)=ngfftrec(3)

!  --algorithm
   ngfftrec(7)=ngfft(7)   ! to be improved for a better non-parallel algorithm - here it is automatically 401
   ngfftrec(8)=ngfft(8)

 end if

!--For now, recursion method doesn't use paralelism on FFT - which would require a great number of processors
 nfftrec = product(ngfftrec(1:3))
 ngfftrec(9:11) = (/0,1,0/)   !--(/ paral, nproc, %me \)
 ngfftrec(12:13) = ngfftrec(2:3)   ! n2proc ! n3proc

 tronc  = all(ngfftrec(:3)/=ngfft(:3))
 call timab(602,2,tsec)

 contains

   function dsq_rec(ii,jj,kk,inf_rmet)

   real(dp) :: dsq_rec
   integer,intent(in) :: ii,jj,kk
   real(dp),intent(in) :: inf_rmet(3,3)
   dsq_rec=sqrt(&
&   inf_rmet(1,1)*dble(ii**2)&
&   +inf_rmet(2,2)*dble(jj**2)&
&   +inf_rmet(3,3)*dble(kk**2)&
&   +two*(inf_rmet(1,2)*dble(ii*jj)&
&   +inf_rmet(2,3)*dble(jj*kk)&
&   +inf_rmet(3,1)*dble(kk*ii)))
 end function dsq_rec


end subroutine getngrec
!!***

!!****f* ABINIT/pspnl_operat_rec
!! NAME
!! pspnl_operat_rec
!!
!! FUNCTION
!! It calculates the non-local projectors used in recursion for any psp non-local:
!! The nl interaction in recursion is $$exp{-V_{NL}/beta}=\sum_A\sum_{lm}
!! \sum{ij}Y_{lm}(\hat{r-R_A}')f^l_i(r-R_A)D^l_{i,j}Y_{lm}(\hat{r-R_A})f^l_j{r-R_A}$$
!! where $D^_{i,j}$ is a matrix  previously (see pspnl_operat_rec).
!! In this routine  the projectors $Y_{lm}(\hat{r-R_A}')f^l_i(r-R_A)$
!! are calculated. So an array of dimensions
!! rset%nl%projec(nfftrec,lmnmax,nlrec%npsp)
!!
!! INPUTS
!! metrec<metricrec_type>=contains information concerning metric in
!!         recursion: grid_step, metric, infinitesimal volume
!! ngfftrec(18)=is the ngfft grid (truncated if different from ngfft) of recursion
!! debug=debug variable
!!
!!
!! OUTPUT
!! nlrec<nlpsprec_type>%projec= array containig the projectors on the real grid
!! nlrec<nlpsprec_type>%intlen= integer linear size of the non-local grid
!!
!! SIDE EFFECTS
!! nlrec<nlpsprec_type> data set of non-local pseudo for recursion
!! The better Interaction length (Intlen) is also calculated and printed but
!! recursion use intlen=ngfftrec/2
!!
!! NOTES
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      gamma_function,initylmr,wrtout
!!
!! SOURCE


subroutine pspnl_operat_rec(nlrec,metrec,ngfftrec,debug)

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: debug
 type(metricrec_type),intent(in) ::metrec
 type(nlpsprec_type),intent(inout) :: nlrec
!arrays
 integer,intent(in) :: ngfftrec(18)
!Local variables-------------------------------
 integer :: ii,intlen
 integer :: iangol,ipsp,iproj
 integer :: mpsang,jj,kk,rr
 integer :: nfftrec
 integer :: ilmn,il,ilm,in,lmnmax
 real(dp) :: raggio,rloc,denom,step
 real(dp) :: delta_out,partial,err
 character(len=500) :: msg
 real(dp) :: part_sum(3)
 real(dp) :: ylmr_gr_dm(0,0,0)
 real(dp),allocatable :: ylmr(:,:),proj_arr(:,:,:)
 real(dp),allocatable :: radloc(:,:),nrm(:)

! *************************************************************************

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_operat_rec : enter '
   call wrtout(std_out,msg,'PERS')
 end if

!#####################################################################
!--CALCULATE THE (SEMI-)MAXIMUM INTERVAL WHERE ALL THE PROJECTORS ARE
!DIFFERENT TO ZERO.
!--For any pseudo potential:
 delta_out = zero
 step = metrec%tr(1)*half !--Scanning step= grid step/2


 do ipsp = 1, nlrec%npsp !--Loop on the pseudos

!  --For any angular moment:
   do iangol = 0,maxval(nlrec%indlmn(1,:,ipsp)) !--Loop on the angular moment
     rloc = nlrec%radii(iangol+1,ipsp) !--Local radius

!    --For any projector
     do iproj = 1,nlrec%pspinfo(iangol+1,ipsp)
!      --Starting point to searching when the projector goes to zero.
!      this correspond to twice the radius wher the projector has its maximum
       raggio = two*sqrt(real(-2+2*iproj+iangol,dp))*rloc
!      --Caculate the gamma func at the denominator
       call gamma_function(real(iangol+2*iproj,dp)-half,denom)
!      --Find the zero
!      --The following while cycle should be replaced by a bisection
!      --method. Bucause this is calculated only 1 time it is not very
!      important.
       err = one
       ii=0
!      tolloc = 1.d0*abs(minval(nlrec%mat_exp_psp_nl(:nlrec%pspinfo(iangol+1,ipsp),:nlrec%pspinfo(iangol+1,ipsp),iangol+1,ipsp)))
       do while(abs(err)>1.d-2)
         raggio = raggio + step
         err = project_prec(raggio,iproj,iangol,rloc)/sqrt(denom)
         ii = ii+1
       end do
       write(std_out,*)'local delta',raggio,ii
       delta_out=maxval((/ delta_out,raggio /))
     end do !end loop on projectors

   end do !enddo on angular moment
 end do !enddo on pseudos

!--CALCULATE how many grid steps correspond to delta_out
 intlen = int(delta_out/metrec%tr(1))
!--I want that intlen is odd
 if(mod(intlen,2)==0) intlen = intlen+1

 write(msg,'(2a,i3,a)') ch10,' Interac. length of non-local psp(grid steps)=',intlen,ch10
 call wrtout(std_out,msg,'COLL')
!#####################################################################

!--Initialisation
 nfftrec = product(ngfftrec(1:3))
 lmnmax = nlrec%lmnmax
 intlen = ngfftrec(1)/2
 nlrec%intlen = intlen !--Setted in recursion variables

!#####################################################################
!--CALCULATE E(q,q')
!--Cration of the exponential*projectors*ylm matrix

!--Initialisation
 ABI_ALLOCATE(nlrec%projec,(nfftrec,lmnmax,nlrec%npsp))
 nlrec%projec = zero
 ABI_ALLOCATE(radloc,(3,nfftrec))
 radloc = zero
 ABI_ALLOCATE(nrm,(nfftrec))
 nrm = zero

!--Loop on pseudo types
 pseudodo: do ipsp = 1, nlrec%npsp
!  --Control if the psp is non-local, else continue
   if(all(nlrec%pspinfo(:,ipsp)==0)) cycle
!  --Vector which stores localy the upper part of symmetrical
!  matrix of the exponential of the non-local operator
   mpsang = maxval(nlrec%indlmn(1,:,ipsp))+1
   ABI_ALLOCATE(proj_arr,(nfftrec,maxval(nlrec%pspinfo(:,ipsp)),mpsang))
   ABI_ALLOCATE(ylmr,(mpsang*mpsang,nfftrec))
   proj_arr = zero
   ylmr = zero

!  !debug
!  write(std_out,*)'mpsang,proj num',mpsang,maxval(nlrec%pspinfo(:,ipsp))
!  !enddebug

!  --Calculate the projctors
   do iangol = 0,mpsang-1
     rloc = nlrec%radii(iangol+1,ipsp)
     do iproj = 1,nlrec%pspinfo(iangol+1,ipsp)
       call gamma_function(real(iangol+2*iproj,dp)-half,denom)
       denom = one/sqrt(denom)
       do ii = 0,ngfftrec(1)-1 !--3-loop on coordinates
         do jj = 0,ngfftrec(2)-1
           do kk = 0,ngfftrec(3)-1
!            --Calculate the radii
             part_sum(:) = real((/ ii,jj,kk /)-intlen,dp)*(metrec%tr)
             rr = 1+ii+(jj+kk*ngfftrec(2))*ngfftrec(3)
             radloc(:,rr) = part_sum
             nrm(rr) = sqrt(sum(part_sum(:)**two))
             partial = project_prec(nrm(rr),iproj,iangol,rloc)*denom
             if(abs(partial)>tol12 ) proj_arr(rr,iproj,iangol+1) = partial
           end do
         end do
       end do  !--End 3-loop on coordinates
     end do
   end do


!  -------------------------------------------------------------
!  --Calculate the spherical harmonics (Verified: it works well)
   call initylmr(mpsang,1,nfftrec,nrm(:),1,radloc(:,:),ylmr(:,:),ylmr_gr_dm)
!  -------------------------------------------------------------


   do ilmn = 1,lmnmax
     ilm = nlrec%indlmn(4,ilmn,ipsp)
     il = nlrec%indlmn(1,ilmn,ipsp)+1
     in = nlrec%indlmn(3,ilmn,ipsp)
     write(msg,'(2a,i3,2i2)')ch10,'lm,l,n',ilm,il,in
     call wrtout(std_out,msg,'COLL')

     nlrec%projec(:,ilmn,ipsp) = ylmr(ilm,:)*proj_arr(:,in,il)
   end do

   ABI_DEALLOCATE(ylmr)
   ABI_DEALLOCATE(proj_arr)
 end do pseudodo !--end loop on pseudo types


 ABI_DEALLOCATE(radloc)
 ABI_DEALLOCATE(nrm)

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_operat_rec : exit '
   call wrtout(std_out,msg,'PERS')
 end if

 contains

   function project_prec(raggio,iproj,iangol,rloc)
!--Analytical expression of the projectors in hgh-pspeudopotential
!--The gamma function at denominator is missing
   real(dp) :: project_prec
   integer,intent(in) :: iproj,iangol
   real(dp),intent(in) :: raggio,rloc

   project_prec=sqrt2*(raggio/rloc)**real((iangol+2*(iproj-1)),dp)*&
&   exp(-((raggio/rloc)**two)*half)/rloc**onehalf
 end function project_prec

end subroutine pspnl_operat_rec
!!***

!!****f* ABINIT/pspnl_hgh_rec
!! NAME
!! pspnl_hgh_rec
!!
!! FUNCTION
!! This routine computes the matrices S_kk'^{l,A} of the projectors
!! (it is the exponential of the overlap matrix). It coorresponds to the matrix:
!!   $$\left[(U_l)^{-1}*Exp(-temperature D_l )*U_l* (g_l)^{-1} -Identity\right]_kk'
!!   where (U_l)^-1* D_l* U_l = h^lg_l.
!! $g_l = <f^l_k|f^l_{k'}>$ is the overlap matrix between projectors
!! and $h^l_{kk'}$ is the strength matrix of the projectors.
!! It calulates also the strength eigenvalues and eigenvectors of $h$,
!! used in the calculus of non-local energy
!!
!! INPUTS
!!  temperature=4*rtrotter/beta=4*rtrotter*tsmear: the effective temp. in  recursion
!!  psps <type(pseudopotential_type)>=variables related to pseudo-potentials
!!  debug=debug variable
!!
!! OUTPUT
!!
!!  nlrec%mat_exp_psp_nl=the matrix of the exponential of the projectors:
!!   for any psp, for any angular moment:
!!   h_ij=strength; g_ij=ovelap => exp(-h.g/temp/4p).g^(-1)
!!  nlrec%radii=Local radii of nl psp
!!  nlrec%pspinfo(:,:) for any typat:  (momang,typat)=number of projectors
!!  nlrec%eival(:,:,:) for any psp, any momang, the eigenvalues of the
!!    strength matrix H: eival(:,mang,psp)
!!  nlrec%eivec(:,:,:,:)for any psp, any momang, the eigenvectors of the
!!    strength matrix H: eivec(:,:,mang,psp)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      dgetrf,dgetri,dsyev,exp_mat,gamma_function,set2unit,wrtout
!!
!! SOURCE

subroutine pspnl_hgh_rec(psps,temperature,nlrec,debug)

 use m_linalg_interfaces
 implicit none

!Arguments -----------------------------------
!scalars
 real(dp),intent(in) :: temperature
 logical,intent(in) :: debug
 type(pseudopotential_type),intent(in) :: psps
 type(nlpsprec_type),intent(inout) :: nlrec
!arrays
!Local variables-------------------------------
!scalars
 integer,parameter :: maxsize=3
 integer,parameter :: lwork=(1+32)*maxsize
 integer :: iangol,ipseudo,info,nproj
 integer :: g_mat_size,ii,nproj2
 real(dp) :: denom_1,denom_2,tot_proj
 character(len=500) :: msg
!arrays
 integer :: ipvt(1)
 !real(dp) :: rwork(2*maxsize)
 real(dp) :: h_mat_init(3,3), rework(lwork)
 real(dp), allocatable :: g_mat(:,:),h_mat(:,:),eig_val_h(:)
 real(dp), allocatable :: identity(:,:),inv_g_mat(:,:),u_mat(:,:)

 complex(dpc),allocatable :: hg_mat(:,:)
! *************************************************************************

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_hgh_rec : enter '
   call wrtout(std_out,msg,'COLL')
 end if

!--For any pseudo potential:
 do ipseudo = 1, psps%npsp !--Loop on the pseudos
   write(msg,'(a,80a)')' pseudo file',('-',ii=1,10)
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a)') psps%filpsp(ipseudo)
   call wrtout(std_out,msg,'COLL')

!  --For any angular moment:
   do iangol = 0,psps%mpsang-1 !--Loop on the angular moment

!    --Local radius
     nlrec%radii(iangol+1,ipseudo) = psps%gth_params%psppar(iangol+1,0,ipseudo)

!    --Strenghts of non-local projectors  (matrix h)
!    --Diagonal part:
     h_mat_init = zero
     h_mat_init(1,1) = psps%gth_params%psppar(iangol+1,1,ipseudo)
     h_mat_init(2,2) = psps%gth_params%psppar(iangol+1,2,ipseudo)
     h_mat_init(3,3) = psps%gth_params%psppar(iangol+1,3,ipseudo)
!    --Out-diagonal part
!    --Depending on angular moment the projectors
!    strength is calculated differently
     select case(iangol)
     case(0)
       h_mat_init(1,2) = -half*sqrt(3.d0/5.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  half*sqrt(5.d0/21.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -half*sqrt(100.d0/63.d0)*h_mat_init(3,3)
     case(1)
       h_mat_init(1,2) = -half*sqrt(5.d0/7.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  sixth*sqrt(35.d0/11.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -14.d0/six/sqrt(11.d0) *h_mat_init(3,3)
     case(2)
       h_mat_init(1,2) = -half*sqrt(7.d0/9.d0)*h_mat_init(2,2)
       h_mat_init(1,3) =  half*sqrt(63.d0/143.d0)*h_mat_init(3,3)
       h_mat_init(2,3) = -nine/sqrt(143.d0)*h_mat_init(3,3)
     case(3)
       h_mat_init(1,2) = zero;  h_mat_init(1,3) = zero;  h_mat_init(2,3) = zero;
     case default
       write(msg,'(a)')' error angular: moment component'
       call wrtout(std_out,msg,'COLL')
     end select



!    --Real dimensions of projectors.
     g_mat_size = count(abs((/ (h_mat_init(ii,ii),ii=1,3) /))>1.d-8)
     nlrec%pspinfo(iangol+1,ipseudo) = g_mat_size
     write(msg,'(a,i2,a,i2)')' ang. moment=',iangol,', N projectors=',g_mat_size
     call wrtout(std_out,msg,'COLL')
     if (g_mat_size>0) then
!      --Identity matrix
       ABI_ALLOCATE(identity,(g_mat_size,g_mat_size))
       call set2unit(identity)
!      identity = zero
!      identity(:,1) = one
!      identity(:,:) = cshift(array=identity,shift=(/ (-ii,ii=0,g_mat_size) /), dim=2 )


!      ############## CALCULOUS OF THE EIGEN_SPACE OF THE PROJECTORS STRENGTHS ##################
!      --Inverse of the matrix h
       ABI_ALLOCATE(eig_val_h,(g_mat_size))
       ABI_ALLOCATE(u_mat,(g_mat_size,g_mat_size))
!      --u-mat will contain the eigenvectors of h_mat_init
       u_mat = h_mat_init(:g_mat_size,:g_mat_size)

!      write(std_out,*)'hmat_init'
!      do ii=1,g_mat_size
!      write(std_out,*)h_mat_init(ii,:)
!      end do
       call DSYEV('v','u',g_mat_size,u_mat,g_mat_size,eig_val_h,rework,lwork,info)

!      --THE DIAGONAL MATRIX IS GIVEN BY D=U^t.H.U
!      (eival=transpose(eivec).h_mat_init.eivec)
       write(msg,'(a,3d10.3)')'  eigenvalues=',eig_val_h
       call wrtout(std_out,msg,'COLL')
!      write(std_out,*)'autovec';write(std_out,*)u_mat

       nlrec%eival(:g_mat_size,1+iangol,ipseudo) = eig_val_h
       nlrec%eivec(:g_mat_size,:g_mat_size,1+iangol,ipseudo) = u_mat
       ABI_DEALLOCATE(eig_val_h)
       ABI_DEALLOCATE(u_mat)

!      ##########END  CALCULOUS OF THE EIGEN_SPACE OF THE PROJECTORS STRENGTHS ##################

       ABI_ALLOCATE(g_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(inv_g_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(h_mat,(g_mat_size,g_mat_size))
       ABI_ALLOCATE(hg_mat,(g_mat_size,g_mat_size))

       g_mat(:,:) = one
       h_mat(:,:) = zero
       h_mat(:,:) = h_mat_init(:g_mat_size,:g_mat_size)

!      -------------------------------------------------------
!      --Matrix  of the overlap between projetors (matrix g)
!      and the h matrix of strength
       do nproj = 1,g_mat_size-1
         do nproj2 = 1+nproj,g_mat_size
           tot_proj = zero
!          --Analytic value of overlap
!          g_ij=Gamma[-1/2+i+j+l]/Sqrt(Gamma[-1/2+i+iangol]*Gamma[-1/2+j+iangol])
           call gamma_function(-half+real(nproj+nproj2+iangol,dp),tot_proj)
           call gamma_function(-half+real(iangol+2*nproj,dp),denom_1)
           call gamma_function(-half+real(iangol+2*nproj2,dp),denom_2)

           g_mat(nproj,nproj2) = tot_proj/sqrt(denom_1*denom_2)
           g_mat(nproj2,nproj) = g_mat(nproj,nproj2)

           h_mat(nproj,nproj2) = h_mat_init(nproj,nproj2)
           h_mat(nproj2,nproj) = h_mat_init(nproj,nproj2)
         end do
       end do

!      --Inverse of the overlap matrix g
       inv_g_mat = g_mat
       call DGETRF(g_mat_size,g_mat_size,inv_g_mat,g_mat_size,ipvt,info)
       call DGETRI(g_mat_size,inv_g_mat,g_mat_size,ipvt,rework,lwork,info)


!      -----------------------------------------------------------
!      --Now it calculates the exponential of the matrix h.g
       hg_mat = matmul(h_mat,g_mat)

       call exp_mat(hg_mat,g_mat_size,-one/temperature)

!      --(exp(h.g)-Identity).(g^-1)
       hg_mat = hg_mat-identity(:,:)


!      --results on output
       nlrec%mat_exp_psp_nl(:g_mat_size,:g_mat_size,1+iangol,ipseudo) = matmul(real(hg_mat,dp),inv_g_mat)

!      write(std_out,*) nlrec%mat_exp_psp_nl(:g_mat_size,:g_mat_size,1+iangol,ipseudo)

       ABI_DEALLOCATE(g_mat)
       ABI_DEALLOCATE(hg_mat)
       ABI_DEALLOCATE(h_mat)
       ABI_DEALLOCATE(inv_g_mat)
       ABI_DEALLOCATE(identity)
     end if

   end do !enddo on angular moment
 end do !enddo on pseudos


 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_hgh_rec : exit '
   call wrtout(std_out,msg,'COLL')
 end if

end subroutine pspnl_hgh_rec
!!***

end module m_rec
!!***
