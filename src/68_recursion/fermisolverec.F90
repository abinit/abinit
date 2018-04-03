!{\src2tex{textfont=tt}}
!!****f* ABINIT/fermisolverec
!! NAME
!! fermisolverec
!!
!! FUNCTION
!! This routine computes the fermi energy in order to have a given number of
!! valence electrons in the recursion method, using a Ridder s Method
!! 
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  debug_rec=debugging variable
!!  nb_rec=order of recursion
!!  nb_point=number of discretization point in one dimension (=n1=n2=n3)
!!  temperature=temperature (Hartree)
!!  trotter=trotter parameter
!!  nelect=number of valence electrons (dtset%nelect)
!!  acc=accuracy for the fermi energy
!!  max_it=maximum number of iteration for the Ridder's Method
!!  long_tranche=number of point computed by thi proc
!!  mpi_enreg=information about MPI parallelization
!!  inf_ucvol=infinitesimal unit cell volume
!!  gputopo=true if topology gpu-cpu= 2 or 3
!! 
!! OUTPUT
!! 
!! SIDE EFFECTS
!!  fermie=fermi energy
!!  rho=density, recomputed for the new fermi energy
!!  a, b2 : coefficient given by recursion recomputed for the new fermi energy
!! 
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      alloc_dens_cuda,dealloc_dens_cuda,density_cuda,density_rec,timab,wrtout
!!      xmpi_barrier,xmpi_sum
!!
!! NOTES
!!  at this time :
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fermisolverec(fermie,rho,a,b2,debug_rec,nb_rec, &
  &                      temperature,trotter,nelect, &
  &                      acc, max_it, &
  &                      long_tranche,mpi_enreg,&
  &                      inf_ucvol,gputopo)

 use defs_basis
 use defs_abitypes
 use defs_rectypes
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_time,         only : timab
#ifdef HAVE_GPU_CUDA
 use m_initcuda,only    : cudap
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fermisolverec'
 use interfaces_14_hidewrite
 use interfaces_68_recursion, except_this_one => fermisolverec
!End of the abilint section

 implicit none

!Arguments -------------------------------
 !scalars
 integer,intent(in) :: long_tranche,max_it,nb_rec,trotter
 logical,intent(in) :: debug_rec,gputopo
 real(dp),intent(in) :: acc,inf_ucvol,nelect,temperature
 real(dp), intent(inout) :: fermie
 type(MPI_type),intent(in) :: mpi_enreg
 !arrays 
 real(dp), intent(inout) :: a(0:nb_rec,long_tranche), b2(0:nb_rec,long_tranche)
 real(dp), intent(inout) :: rho(long_tranche)

!Local variables-------------------------------
 !scalars  
 integer  ::  ierr,ii,ipointlocal,nn,dim_trott
 real(dp) :: beta,fermieh,fermiel,fermiem,fermienew,nelecth,nelectl,nelectm
 real(dp) :: nelectnew,res_nelecth,res_nelectl,res_nelectm,res_nelectnew
 real(dp) :: rtrotter,ss,fermitol
 character(len=500) :: msg
 !arrays
 real(dp) :: tsec(2)
 real(dp) :: rhotry(long_tranche)
 !no_abirules
#ifdef HAVE_GPU_CUDA
 integer :: swt_tm,npitch
 real(cudap) :: rhocu(long_tranche)
 real(dp) :: tsec2(2)
#endif

! *************************************************************************

#ifdef HAVE_GPU_CUDA
 swt_tm = 0
#endif

 call timab(609,1,tsec)

 beta = one/temperature
 rtrotter  = max(half,real(trotter,dp))
 dim_trott = max(0,2*trotter-1)

 write(msg,'(a)')' -- fermisolverec ---------------------------------'
 call wrtout(std_out,msg,'COLL')
 if(debug_rec) then 
   write (msg,'(a,d10.3)')' nelect= ',nelect 
   call wrtout(std_out,msg,'COLL') 
 end if
!initialisation of fermiel
 fermiel = fermie
 call timab(609,2,tsec)

!initialisation fermitol
 fermitol = acc
#ifdef HAVE_GPU_CUDA_SP
 if(gputopo)  fermitol = 1.d-3
#endif

 if(gputopo) then
#ifdef HAVE_GPU_CUDA
   swt_tm = 1
!  allocate array an and bn2 on gpu for computation of trotter formula
   call alloc_dens_cuda(long_tranche,nb_rec,dim_trott,npitch,&
&   real(a,cudap),real(b2,cudap))

   call timab(617,1,tsec)
   call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&   real(fermiel,cudap),real(temperature,cudap),&
&   real(rtrotter,cudap),real(inf_ucvol,cudap),&
&   real(tol14,cudap),&
&   rhocu)
   rhotry = real(rhocu,dp)  
   call timab(617,2,tsec)
#endif
 else
   do ipointlocal = 1,long_tranche
     call density_rec(a(:,ipointlocal),& 
&     b2(:,ipointlocal),& 
&     rhotry(ipointlocal),&
&     nb_rec,fermiel,temperature,rtrotter,dim_trott, &
&     tol14,inf_ucvol)
   end do
 end if

 call timab(609,1,tsec) 
 nelectl = sum(rhotry)
 call xmpi_sum( nelectl,mpi_enreg%comm_bandfft,ierr)
 res_nelectl = inf_ucvol*nelectl - nelect

 if (res_nelectl /= zero) then 
!  initialisation of fermih
!  excess of electrons -> smaller fermi
   res_nelecth = zero
   ii = 1
   fermieh = fermie - ten*sign(one,res_nelectl)*temperature  
   do while(ii<6 .and. res_nelecth*res_nelectl>=0)
     fermieh = fermieh - ten*sign(one,res_nelectl)*temperature     
     call timab(609,2,tsec)

     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermieh,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  & 
&         b2(:,ipointlocal), & 
&         rhotry(ipointlocal), &
&         nb_rec,fermieh,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)
       end do
     end if
     call timab(609,1,tsec)
     nelecth = sum(rhotry)
     call xmpi_sum( nelecth,mpi_enreg%comm_bandfft ,ierr);
     res_nelecth = inf_ucvol*nelecth - nelect

     if(debug_rec) then
       write (msg,'(a,es11.4e2,a,es11.4e2)') ' Fermi energy interval',fermieh,' ',fermiel
       call wrtout(std_out,msg,'COLL') 
     end if
     ii = ii +1
   end do

   if (res_nelecth*res_nelectl>0) then
     write (msg,'(4a)')' fermisolverec : ERROR- ',ch10,&
&     ' initial guess for fermi energy doesnt permit to  find solutions in solver',ch10
     MSG_ERROR(msg)
   end if

!  MAIN LOOP   ------------------------------------------------------
   main : do nn=1,max_it     
!    fermiem computation
     fermiem = 0.5d0*(fermiel+fermieh) 

!    nelectm = zero
     call timab(609,2,tsec)

     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermiem,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  & 
&         b2(:,ipointlocal), & 
&         rhotry(ipointlocal), &
&         nb_rec,fermiem,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)       
       end do
     end if

     call timab(609,1,tsec)
     nelectm = sum(rhotry)
     call xmpi_sum( nelectm,mpi_enreg%comm_bandfft,ierr)
     res_nelectm = inf_ucvol*nelectm - nelect

!    new guess
     ss = sqrt(res_nelectm**two-res_nelectl*res_nelecth)
     fermienew = fermiem + (fermiem-fermiel)*sign(one, res_nelectl-res_nelecth)*res_nelectm/ss

     call timab(609,2,tsec)
     if(gputopo) then
#ifdef HAVE_GPU_CUDA
       call timab(617,1,tsec)
       call density_cuda(npitch,long_tranche,nb_rec,dim_trott,&
&       real(fermienew,cudap),real(temperature,cudap),&
&       real(rtrotter,cudap),real(inf_ucvol,cudap),&
&       real(tol14,cudap),&
&       rhocu)
       rhotry = real(rhocu,dp)
       call timab(617,2,tsec)
#endif
     else
       do ipointlocal = 1,long_tranche
         call density_rec(a(:,ipointlocal),  & 
&         b2(:,ipointlocal), & 
&         rhotry(ipointlocal), &
&         nb_rec,fermienew,temperature,rtrotter,dim_trott, &
&         tol14,inf_ucvol)       
       end do
     end if

     call timab(609,1,tsec)
     nelectnew = sum(rhotry)
     call xmpi_sum( nelectnew,mpi_enreg%comm_bandfft ,ierr); 
     res_nelectnew = inf_ucvol*nelectnew - nelect

!    fermiel et fermieh for new iteration
     if (sign(res_nelectm,res_nelectnew) /= res_nelectm) then
       fermiel = fermiem
       res_nelectl = res_nelectm
       fermieh = fermienew
       res_nelecth = res_nelectnew
     else if (sign(res_nelectl,res_nelectnew) /= res_nelectl) then
       fermieh = fermienew
       res_nelecth = res_nelectnew
     else if (sign(res_nelecth,res_nelectnew) /= res_nelecth) then
       fermiel = fermienew
       res_nelectl = res_nelectnew
     end if

!    are we within the tolerance ?
     if ((abs(res_nelectnew) < fermitol).or.(nn == max_it)) then
       fermie = fermienew
       rho = rhotry
       if(debug_rec) then
         write (msg,'(a,es11.4e2,a,i4)')' err, num_iter ', res_nelectnew, ' ',nn
         call wrtout(std_out,msg,'COLL')
         write(msg,'(a,50a)')' ',('-',ii=1,50)
         call wrtout(std_out,msg,'COLL')
       end if
       exit main
     end if

   end do main

 end if

#ifdef HAVE_GPU_CUDA
!deallocate array on GPU
 if(gputopo) then
   call dealloc_dens_cuda()
 end if
 call timab(613+swt_tm,1,tsec2)  !!--start time-counter: sync gpu-cpu
 call xmpi_barrier(mpi_enreg%comm_bandfft)
 call timab(613+swt_tm,2,tsec2)  !!--stop time-counter: sync gpu-cpu
#endif


 call timab(609,2,tsec)
end subroutine fermisolverec
!!***
