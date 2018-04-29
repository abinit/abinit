!{\src2tex{textfont=tt}}
!!****f* ABINIT/first_rec
!! NAME
!! first_rec
!!
!! FUNCTION
!! When recursion method is used, in the first step this routine
!! compute some quantities which are used in the rest of the calculation.
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2018 ABINIT group (MMancini)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset:
!!   | recgratio =fine/coarse grid ratio
!!   | recptrott =trotter parameter
!!   | tsmear    =temperature
!!   | recrcut   =tut radius in recursion (range of iteration)
!!   | ngfft(18) =FFT grid used as real (fine) grid in recursion
!!  psps <type(pseudopotential_type)>=variables related to pseudo-potentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  rset <type(recursion_type)>=variables related to recursion method
!!   | debug<logical> = T if debugging is used
!!   | inf <type(metricrec_type)>=information concerning the infinitesimal metrics
!!   | ngfftrec(18) =truncated (or not, if not ngfftrec=ngfft)FFT grid used as real grid in recursion.
!!   | nfftrec =product(ngfftrec(1:3))
!!   | tronc<logical> = T if truncation is effectively used
!!   | ZT_p = fourier transform of the green_kernel calculated on the fine grid
!!
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cpu_distribution,cudarec,get_pt0_pt1,green_kernel,init_nlpsprec
!!      random_number,recursion,reshape_pot,timab,timein,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#if defined HAVE_GPU_CUDA
#include "cuda_common.h"
#endif

#include "abi_common.h"

subroutine first_rec(dtset,psps,rset)

 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use defs_rectypes
 use m_profiling_abi
 use m_errors

 use m_time,       only : timein, timab
 use m_rec,        only : init_nlpsprec,cpu_distribution
 use m_rec_tools,  only : get_pt0_pt1,reshape_pot
 use m_xmpi,       only:  xmpi_sum
#if defined HAVE_GPU_CUDA
 use m_initcuda,only    : cudap
 use m_hidecudarec,only : cudarec,CleanRecGPU
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'first_rec'
 use interfaces_14_hidewrite
 use interfaces_68_recursion, except_this_one => first_rec
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rset
!Local variables-------------------------------
!scalars
 integer  :: nfftrec,trotter,ii,dim_trott
 real(dp) :: tsmear,beta,rtrotter
 character(len=500) :: msg
!arrays
 integer  :: ngfftrec(18)
 real(dp) :: tsec(2)
#ifdef HAVE_GPU_CUDA
 integer  :: max_rec,ierr,testpts,swt_tm
 real(dp) :: rho,tm_ratio
 real(dp) :: time_cu,time_f
 type(recursion_type) :: rset_test
 type(recparall_type) :: parold
 integer :: trasl(3)
 real(dp) :: tsec2(2),tsec3(2)
 real(dp) :: aloc(0,1),b2loc(0,1)
 real(dp) :: dm_projec(0,0,0,1,1)
 real(dp) :: exppot(0:dtset%nfft-1)
 real(dp),allocatable :: exppotloc(:)
 real(cudap),allocatable :: aloc_cu(:),b2loc_cu(:)
#endif

! *************************************************************************

 call timab(601,1,tsec)  !!--Start time-counter: initialisation

 MSG_WARNING("RECURSION")
 if(dtset%recgratio>1) then
   write(msg,'(a)')'COARSE GRID IS USED'
   call wrtout(std_out,msg,'COLL')
 end if

!--Initialisation
 trotter = dtset%recptrott  !--Trotter parameter
 tsmear  = dtset%tsmear     !--Temperature
 beta    = one/tsmear       !--Inverse of temperature

!--Rewriting the trotter parameter
 dim_trott = max(0,2*trotter-1)
 rtrotter  = max(half,real(trotter,dp))

 write (msg,'(2a)')ch10,'==== FIRST CYCLE RECURSION ========================='
 call wrtout(std_out,msg,'COLL')


 ngfftrec = rset%ngfftrec
 nfftrec = rset%nfftrec
!------------------------------------------------
!--TRONCATION OF THE BOX: determines new dimensions
!--Now in InitRec
!--------------------------------------------------------
!--DEFINITION PAW VARIABLES COARSE-FINE GRID  TO USE TRANSGRID--INGRID FUNCTIONS
!--Now these variables are defined into gstate by InitRec

!--------------------------------------------------------
!--COMPUTATION OF THE FOURIER TRANSFORM OF THE GREEN KERNEL (only once)
 write (msg,'(a)')' - green kernel calculation -----------------------'
 call wrtout(std_out,msg,'COLL')
 ABI_ALLOCATE(rset%ZT_p,(1:2,0: nfftrec-1))
 call timab(601,2,tsec)
 call green_kernel(rset%ZT_p,rset%inf%rmet,rset%inf%ucvol,rtrotter/beta,rset%mpi,ngfftrec,nfftrec)
 call timab(601,1,tsec)
 write(msg,'(a,50a)')' ',('-',ii=1,50)
 call wrtout(std_out,msg,'COLL')
!!--end computation of the fourier transform of the Green kernel

!!-----------------------------------
!!--ROUTINE FOR THE CALCULATION OF THE NON-LOCAL PSEUDO
!--Now these variables here by  Init_nlpspRec
 call Init_nlpspRec(four*tsmear*rtrotter,psps,rset%nl,rset%inf,rset%ngfftrec,rset%debug)

!!-----------------------------------
!--Load distribution on procs when GPU are present
#if defined HAVE_GPU_CUDA

!--Test timing only if exists GPU and they are not equal to the cpus
 if(rset%tp == 4) then
   parold = rset%par
   ii = 0
   time_f = zero
   time_cu = zero
   call random_number(exppot)  !   exppot = one

   if(rset%gpudevice == -1) then
!    --Test CPUS
     swt_tm = 0
     testpts = min(rset%par%npt, 20)
     call timein(tsec2(1),tsec2(2))
     if(rset%tronc) then
       ABI_ALLOCATE(exppotloc,(0:nfftrec-1))
       do while(ii< testpts)
         trasl = -(/1,2,3/)+ngfftrec(:3)/2
         call reshape_pot(trasl,dtset%nfft,nfftrec,dtset%ngfft(:3),ngfftrec(:3),&
&         exppot,exppotloc)
         call recursion(exppotloc,0,0,0, &
&         aloc, &
&         b2loc, &
&         rho,&
&         0, rset%efermi,tsmear,rtrotter,dim_trott, &
&         rset%ZT_p, &
&         dtset%rectolden,dtset%typat, &
&         rset%nl,&
&         rset%mpi,nfftrec,ngfftrec,rset%inf,&
&         6,dtset%natom,dm_projec,0)
         ii=ii+1
       end do
       ABI_DEALLOCATE(exppotloc)
     else
       do while(ii< testpts)
         call recursion(exppot,0,0,0, &
&         aloc, &
&         b2loc, &
&         rho,&
&         0, rset%efermi,tsmear,rtrotter,dim_trott, &
&         rset%ZT_p, &
&         dtset%rectolden,dtset%typat, &
&         rset%nl,&
&         rset%mpi,nfftrec,ngfftrec,rset%inf,&
&         6,dtset%natom,dm_projec,0)
         ii=ii+1
       end do
     end if
     call timein(tsec3(1),tsec3(2))
     time_f = (tsec3(1)-tsec2(1))/real(testpts,dp)
     time_f = time_f*time_f
   else
!    --Test GPUS
     swt_tm = 1
     rset_test = rset
     rset_test%GPU%par%npt = max(rset%GPU%nptrec,100)
     rset_test%min_nrec = 0
     call get_pt0_pt1(dtset%ngfft(:3),dtset%recgratio,0,&
&     rset_test%GPU%par%npt,rset_test%GPU%par)


     ABI_ALLOCATE(aloc_cu,(rset_test%GPU%par%npt))
     ABI_ALLOCATE(b2loc_cu,(rset_test%GPU%par%npt))
     call timein(tsec2(1),tsec2(2))
     call cudarec(rset_test, exppot,aloc_cu,b2loc_cu,&
&     beta,trotter,dtset%rectolden,dtset%recgratio,dtset%ngfft,max_rec)
     call timein(tsec3(1),tsec3(2))
     ABI_DEALLOCATE(aloc_cu)
     ABI_DEALLOCATE(b2loc_cu)

     time_cu = (tsec3(1)-tsec2(1))/real(rset_test%GPU%par%npt,dp)
     time_cu = time_cu*time_cu
   end if


!  --Get Total Times
   call xmpi_sum(time_f,rset%mpi%comm_bandfft,ierr)
   call xmpi_sum(time_cu,rset%mpi%comm_bandfft,ierr)

!  --Average Total Times
   time_f   = sqrt(time_f/real(rset%mpi%nproc-rset%ngpu,dp))
   time_cu  = sqrt(time_cu/real(rset%ngpu,dp))
   tm_ratio = time_f/time_cu


   write(msg,'(3(a25,f10.5,a))')&
&   ' Time for cpu recursion ',time_f,ch10,&
&   ' Time for gpu recursion ',time_cu,ch10,&
&   ' Time ratio             ',tm_ratio,ch10
   call wrtout(std_out,msg,'COLL')


!  tm_ratio =1.20d2! 0.d0! 1.21d0
   rset%par = parold
!  --Compute the work-load distribution on devices (gpu,cpu)
   if(tm_ratio>1.5d0 .and. time_cu>zero)then
     rset%load = 1
     call cpu_distribution(dtset%recgratio,rset,dtset%ngfft(:3),tm_ratio,1)
   else
     rset%gpudevice = -1
   end if
 end if

#endif


!------------------------------------------------------------
!--DETERMINING WHICH POINT WILL COMPUTE THAT PROC
!--Now these variables are defined into gstate by Init_rec

 write (msg,'(2a)')ch10,'==== END FIRST CYCLE RECURSION ====================='
 call wrtout(std_out,msg,'COLL')
 call timab(601,2,tsec) !!--stop time-counter: initialisation

end subroutine first_rec
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/green_kernel
!! NAME
!! green_kernel
!!
!! FUNCTION
!! this routine compute the fourrier transform of the Green kernel for the
!! recursion method
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  inf_rmet=define the  infinitesimal metric : rprimd*(transpose(rprimd)) divided by the number of discretisation point
!!  inf_ucvol=volume of infinitesimal cell
!!  mult=variance of the Gaussian (=rtrotter/beta)
!!  mpi_enreg=information about MPI parallelization
!!  ngfft=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nfft=total number of fft grid points
!!  debug_rec=debugging variable
!!
!! OUTPUT
!!  ZT_p=fourier transforme of the Green kernel
!!
!! PARENTS
!!      first_rec
!!
!! CHILDREN
!!      fourdp,timab,wrtout
!!
!! NOTES
!!  at this time :
!!       - need a rectangular box
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

 use m_time,        only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'green_kernel'
 use interfaces_14_hidewrite
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: inf_ucvol,mult
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: inf_rmet(3,3)
 real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)

!Local variables-------------------------------
!scalars
 integer,parameter :: n_green_max=5
 integer :: ii,isign,jj,kk,n_green,xx,yy,zz
 real(dp) :: acc, norme
 character(len=500) :: msg
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: T_p(:)

! *************************************************************************

 call timab(603,1,tsec)

 norme = (mult/pi)**(onehalf)

 ABI_ALLOCATE(T_p,(0:nfft-1))

!n_green should be better chosen for non rectangular cell
 do xx=1, n_green_max
   n_green = xx
   if(exp(-mult*dsq_green(xx*ngfft(1),0,0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,xx*ngfft(2),0,inf_rmet))<tol14 &
&   .and. exp(-mult*dsq_green(0,0,xx*ngfft(3),inf_rmet))<tol14 ) exit
 end do

 acc = zero
 T_p = zero
 do kk = 0,ngfft(3)-1
   do jj = 0,ngfft(2)-1
     do ii = 0,ngfft(1)-1

       do xx=-n_green,n_green-1
         do yy=-n_green,n_green-1
           do zz=-n_green,n_green-1

             T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)+ &
&             exp(-mult*dsq_green(ii+xx*ngfft(1),jj+yy*ngfft(2),kk+zz*ngfft(3),inf_rmet))

           end do
         end do
       end do

       T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = norme*T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
       acc = acc + inf_ucvol* T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)

     end do
   end do
 end do

 T_p(:)= (one/acc)*T_p(:)

!if(debug_rec)then
 write(msg,'(a,d12.3,2(2a,i8),2(2a,3d12.3),2a,d16.6)')&
& ' on the boundary    ', exp(-mult*dsq_green(ngfft(1),0,0,inf_rmet)),ch10, &
& ' no zero            ', count(T_p>tol14),ch10, &
& ' n_green            ', n_green,ch10, &
& ' erreur_n_green     ', exp(-mult*dsq_green(n_green*ngfft(1),0,0,inf_rmet)), &
& exp(-mult*dsq_green(0,n_green*ngfft(2),0,inf_rmet)), &
& exp(-mult*dsq_green(0,0,n_green*ngfft(3),inf_rmet)),ch10,&
& ' erreur_troncat     ', T_p(ngfft(1)/2),  &
& T_p(ngfft(1)*(ngfft(2)/2)), &
& T_P(ngfft(1)*ngfft(2)*(ngfft(3)/2)),ch10, &
& ' erreurT_p          ',abs(acc-1.d0)
 call wrtout(std_out,msg,'COLL')
!endif


 isign = -1
 call fourdp(1,ZT_p,T_p,isign,mpi_enreg,nfft,ngfft,1,0)

 ABI_DEALLOCATE(T_p)

 ZT_p(:,:) = real(nfft,dp)*ZT_p


 call timab(603,2,tsec)

 contains

   function dsq_green(ii,jj,kk,inf_rmet)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dsq_green'
!End of the abilint section

   real(dp) :: dsq_green
   integer,intent(in) :: ii,jj,kk
   real(dp),intent(in) :: inf_rmet(3,3)
   dsq_green= inf_rmet(1,1)*dble(ii**2)&
&   +inf_rmet(2,2)*dble(jj**2)&
&   +inf_rmet(3,3)*dble(kk**2)&
&   +two*(inf_rmet(1,2)*dble(ii*jj)&
&   +inf_rmet(2,3)*dble(jj*kk)&
&   +inf_rmet(3,1)*dble(kk*ii))
 end function dsq_green

end subroutine green_kernel
!!***
