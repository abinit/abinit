!!****m* ABINIT/m_paw_gaussfit
!! NAME
!!  m_paw_gaussfit
!!
!! FUNCTION
!!  Module to fit PAW related data to sums of gaussians
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_paw_gaussfit

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_paw_numeric, only : paw_splint, paw_spline
 use m_pawrad,      only : pawrad_type, pawrad_init, pawrad_deducer0, pawrad_free, pawrad_ifromr

 implicit none

 private

 public:: gaussfit_projector !fit non-local projectors to a sum of gaussians
 public:: gaussfit_main !main routine to fit
!Routines related to MPI:
 private:: gaussfit_mpi_set_weight
 private:: gaussfit_mpi_remove_item
 private:: gaussfit_mpi_add_item
 private:: gaussfit_mpi_assign
 private:: gaussfit_mpi_main
 private:: gaussfit_mpi_calc_deviation
 private:: gaussfit_mpi_swap
!Routines related to fitting:
 private:: gaussfit_fit !fit a function to gaussians
 private:: gaussfit_calc_deriv_r  !calc. derivatives for real gauss.
 private:: gaussfit_calc_deriv_c  !calc. derivatives for cplex. gauss. of form 1
 private:: gaussfit_calc_deriv_c2 !calc. derivatives for cplex. gauss. of form 2
 private:: gaussfit_calc_deriv_c3 !calc. derivatives for cplex. gauss. of form 3
 private:: gaussfit_calc_deriv_c4 !calc. derivatives for cplex. gauss. of form 4
 private:: gaussfit_rlsf !retreined least squares fit
 private:: gaussfit_chisq_alpha_beta
!set parameters for LSF (for 5 different forms of gauss. sums):
 private:: gaussfit_set_param1
 private:: gaussfit_set_param2
 private:: gaussfit_set_param3
 private:: gaussfit_set_param4
 private:: gaussfit_set_param5
 private:: gaussfit_constrains_init !initialize constrains
 private:: gaussfit_apply_constrains !apply cons. to get new params.

!!***
 integer,private,parameter:: positive=2
 integer,private,parameter:: restricted=3
 integer,private,parameter:: restricted_and_positive=4

CONTAINS
!===========================================================
!!***

!!****f* m_paw_gaussfit/gaussfit_main
!! NAME
!!  gaussfit_main
!!
!! FUNCTION
!!  Fits a given input function f(r) to a sum of gaussians
!!
!! INPUTS
!!  nterm_bounds= sets the minimum and maximum number of terms to be taken into account.
!!  mparam= maximum number of parameters
!!         if(option==1)mparam=nterm_bounds(2)*4
!!         if(option==2)mparam=nterm_bounds(2)*6
!!         if(option==3)mparam=nterm_bounds(2)*2
!!         if(option==4)mparam=nterm_bounds(2)*4
!!  nr= number of real space points
!!  pawrad= pawrad type
!!  option=1 fit to a1 cos(a2 x^2)+ a3 sin( a4 x^2)
!!         2 fit to a1 exp(-a2 x^2)*(a3 cos (a4 x^2) + a5 sin (a6 x^2) )
!!         3 fit to a1 cos (k x^2) + a2 sin (k x^2)
!!         4 fit to a1 exp(-a2 x^2)* (a3 cos(k x^2)+ a4 sin (k x^2))
!!  Given these definitions, the number of complex gaussians are:
!!    ngauss=4*nterm (for option=1,2)
!!    ngauss=2*nterm (for option=3,4)
!!  outfile= filename to write out fitted functions.
!!  rpaw=paw radius
!!  y(nr)= function to fit
!!
!! OUTPUT
!! nparam_out= number of parameters found.
!! param_out(nparam_out)= parameters (coefficients and factors of complex gaussians).
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_main(mparam,nparam_out,nterm_bounds,nr,&
&           param_out,pawrad,option,outfile,rpaw,y,comm_mpi)

!Arguments ------------------------------------
 integer,intent(in)::mparam,nr,nterm_bounds(2)
 integer,intent(in)::option
 integer,intent(in),optional:: comm_mpi
 real(dp),intent(in)::rpaw
 real(dp),intent(inout)::y(nr)
 character(80),intent(in)::outfile
 type(pawrad_type),intent(in) :: pawrad
 integer,intent(out)::nparam_out
 real(dp),intent(out)::param_out(mparam)

!Local variables-------------------------------
!scalars
 logical,parameter::modify_y=.false. !used only for plotting purposes
 integer :: ichisq,ierr,ii,jj
 integer :: master,me
 integer :: maxiter,minterm,my_chisq_size,counts_all
 integer :: ngauss,nparam,nproc,nterm
 integer :: verbosity
 real(dp) :: chisq,chisq_min
!real(dp) :: T1,T2 !uncomment for timming
 !arrays
 integer::constrains(mparam)
 integer::proc_dist(nterm_bounds(1):nterm_bounds(2))
 integer,allocatable::counts(:),disp(:),map_nterm(:)
 real(dp)::limit(mparam),param_tmp(mparam),weight(mparam)
 real(dp),allocatable::chisq_array(:),recv_buf(:),send_buf(:)
 real(dp),allocatable::y_out(:)
 character(len=500) :: msg

! *************************************************************************

!initialize variables
 maxiter=200
!
!initialize mpi quantities:
 master=0; me=0; nproc=1; proc_dist=1;
 if(present(comm_mpi)) then
   me=xmpi_comm_rank(comm_mpi)
   nproc=xmpi_comm_size(comm_mpi)
 end if
 if(nproc>1) then
!  Find distribution (master)
   if(me==master) then
     call gaussfit_mpi_main(nproc,nterm_bounds,proc_dist)
   end if
!  send distribution to all processors
   call xmpi_bcast(proc_dist(nterm_bounds(1):nterm_bounds(2)),&
&    master,comm_mpi,ierr)
 end if
!Set size of chisq treated by each proc
 my_chisq_size=nterm_bounds(2)-nterm_bounds(1)+1 !all terms
 if(nproc>1 .and. .not. me==master) then
   do nterm=nterm_bounds(1),nterm_bounds(2)
     if(.not. proc_dist(nterm)==me+1) cycle
     my_chisq_size=my_chisq_size+1
   end do
 end if

!
!Allocate objects
!
 LIBPAW_ALLOCATE(y_out,(nr))
 LIBPAW_ALLOCATE(chisq_array,(my_chisq_size))
 if(master==me ) then
   LIBPAW_BOUND1_ALLOCATE(map_nterm,BOUNDS(nterm_bounds(1),nterm_bounds(2)))
   jj=1
   do ii=1,nproc
     do nterm=nterm_bounds(1),nterm_bounds(2)
       if(proc_dist(nterm)==ii) then
         map_nterm(nterm)=jj
         jj=jj+1
       end if
     end do
   end do
 end if
!
!fill with zeros
!
 chisq_array=zero
 nparam_out=0
 y_out=0.d0
 verbosity=1 !print the minimum at the terminal
!
 ichisq=0
 do nterm=nterm_bounds(1),nterm_bounds(2)
!  mpi distribution
   if(.not. proc_dist(nterm)==me+1) cycle
   ichisq=ichisq+1

!   call CPU_TIME(T1)

   if(option==1) nparam=nterm*4
   if(option==2) nparam=nterm*6
   if(option==3) nparam=nterm*2
   if(option==4) nparam=nterm*4
!  set initial guess
!   if(option==1) then
!     call gaussfit_set_param3(nterm,nparam,param_tmp(1:nparam),sep(ii))
!   elseif(option==2) then
!     call gaussfit_set_param1(nterm,nparam,nr,&
!&     param_tmp(1:nparam),sep(ii),pawrad%rad(1:nr),y)
   if(option==3) then
     call gaussfit_set_param4(nparam,param_tmp(1:nparam))
   elseif(option==4) then
     call gaussfit_set_param5(nterm,nparam,nr,&
&     param_tmp(1:nparam),rpaw,y)
   end if
!
!
   call gaussfit_constrains_init(weight(1:nparam),constrains(1:nparam),&
&   limit(1:nparam),nparam,nterm,nr,option,rpaw,y)
!
   call gaussfit_fit(chisq,constrains(1:nparam),&
&   limit(1:nparam),maxiter,nparam,nterm,nr,option,outfile,param_tmp,&
&   verbosity,weight(1:nparam),pawrad%rad(1:nr),y,y_out)

!  if there was an error, set chisq to very high
!  if there is an NaN, for instance:
   if(abs(chisq+1.d0)<tol8) chisq=99999
   if(abs(chisq)==chisq*chisq) chisq=99999
   if(chisq .ne. chisq) chisq=99999
!
   chisq_array(ichisq)=chisq

!   call CPU_TIME(T2)
!   print *, 'Time for fit ', T2-T1, 'seconds.'

 end do !nterm

!mpicast results
!send distribution to master
 if(nproc>1) then
!  Prepare communications:
   LIBPAW_ALLOCATE(counts,(nproc))
   counts(:)=0
   do nterm=nterm_bounds(1),nterm_bounds(2)
     counts(proc_dist(nterm))=counts(proc_dist(nterm))+1
   end do
   counts_all=sum(counts)
   LIBPAW_ALLOCATE(send_buf,(counts(me+1)))
   send_buf(:)=chisq_array(1:counts(me+1))
   if(me==master) then
     LIBPAW_ALLOCATE(recv_buf,(counts_all))
   else
     LIBPAW_ALLOCATE(recv_buf,(1))
   end if
   LIBPAW_ALLOCATE(disp,(nproc))
   disp(1)=0
   do ii=2,nproc
     disp(ii)=disp(ii-1)+counts(ii-1)
   end do
!  communicate all info to master
   call xmpi_gatherv(send_buf,counts(me+1),recv_buf,&
&    counts,disp,master,comm_mpi,ierr)
!  fill in chisq_array with all received info:
   if(master==me) then
     do ii=1,counts_all
       chisq_array(ii)=recv_buf(ii)
     end do
   end if
!  Deallocate MPI arrays:
   LIBPAW_DEALLOCATE(recv_buf)
   LIBPAW_DEALLOCATE(counts)
   LIBPAW_DEALLOCATE(disp)
   LIBPAW_DEALLOCATE(send_buf)
 end if

!Print out info:
 if(me==master) then
   write(msg,'(3a)')'Preliminary results (with only 200 iter.):',ch10,'   ngauss    chisq'
   call wrtout(std_out,msg,'COLL')
   do nterm=nterm_bounds(1),nterm_bounds(2)
     if(option==1) ngauss=nterm*4
     if(option==2) ngauss=nterm*4
     if(option==3) ngauss=nterm*2
     if(option==4) ngauss=nterm*2
     write(msg,'(i4,2x,e13.6,1x)')ngauss,&
&     chisq_array(map_nterm(nterm))
     call wrtout(std_out,msg,'COLL')
   end do
 end if

!get minterm for best accuracy:
 if(me==master) then
   chisq_min=999999999
   do nterm=nterm_bounds(1),nterm_bounds(2)
     if(chisq_array(map_nterm(nterm))<chisq_min) then
       chisq_min=chisq_array(map_nterm(nterm))
       minterm=nterm
     end if
   end do

!run again with the best parameters
   nterm=minterm
   if(option==1)nparam=4*nterm
   if(option==2)nparam=6*nterm
   if(option==3)nparam=2*nterm
   if(option==4)nparam=4*nterm

!  set initial guess
   if (option==3) then
     call gaussfit_set_param4(nparam,param_tmp(1:nparam))
   elseif (option==4) then
     call gaussfit_set_param5(nterm,nparam,nr,&
&     param_tmp(1:nparam),rpaw,y)
   end if
!
   call gaussfit_constrains_init(weight(1:nparam),constrains(1:nparam),&
&   limit(1:nparam),nparam,nterm,nr,option,rpaw,y)
!
   verbosity=1;
   maxiter=1000 !this time we do more iterations
   call gaussfit_fit(chisq,constrains(1:nparam),&
&   limit(1:nparam),maxiter,nparam,nterm,nr,option,outfile,param_tmp,&
&   verbosity,weight(1:nparam),pawrad%rad(1:nr),y,y_out)

!  Write out best solution
   write(msg,'(3a)')"Best solution (with more iterations):",ch10,"   ngauss    chisq"
   call wrtout(std_out,msg,'COLL')
   if(option==1) ngauss=nterm*4
   if(option==2) ngauss=nterm*4
   if(option==3) ngauss=nterm*2
   if(option==4) ngauss=nterm*2
   write(msg,'(i4,2x,e13.6,1x)')ngauss,chisq
   call wrtout(std_out,msg,'COLL')
 end if

!Fill output variables
!and communicate results to all procs:
 if(option==1) then
   nparam_out=minterm*4
 elseif (option==2) then
   nparam_out=minterm*6
 elseif (option==3) then
   nparam_out=minterm*2
 elseif (option==4) then
   nparam_out=minterm*4
 end if

!communicate
 if(nproc>1) then
   call xmpi_bcast(nparam_out,master,comm_mpi,ierr)
   call xmpi_bcast(param_tmp(1:nparam_out),master,comm_mpi,ierr)
 end if !nproc>1

 param_out(:)=param_tmp

 if(modify_y) then
!   call xmpi_scatterv(y_out,nr,mpi_displ,y_out,nr,0,comm_mpi,ierr)
   y=y_out !at output modify y for the fitted y
 end if

 LIBPAW_DEALLOCATE(y_out)
 LIBPAW_DEALLOCATE(chisq_array)
 if(me==master) then
   LIBPAW_DEALLOCATE(map_nterm)
 end if

end subroutine gaussfit_main
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_set_weight
!! NAME
!!  gaussfit_mpi_set_weight
!!
!! FUNCTION
!! It sets a weight to the number of
!! Gaussians used.
!! This was calculated by measuring the time
!! it takes to fit a projector with different
!! number of gaussians.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_set_weight(f,x)

!Arguments ------------------------------------
 integer,intent(in)::x
 integer,intent(out)::f

!Local variables ------------------------------
 real(dp)::a,b,c,d,ff,xx

!************************************************************************

 !The following parameters were obtained
 !from the time (in seconds)
 !it takes to fit a given projector
 !using from 1 to 100 gaussians.
 a               = -0.374137d0      !  +/- 2.013        (538.1%)
 b               = 0.207854d0       !  +/- 0.3385       (162.8%)
 c               = 0.0266371d0      !  +/- 0.01534      (57.59%)
 d               = 0.000152476d0    !  +/- 0.0001978    (129.7%)

 xx=real(x,dp)
 ff=a+b*xx+c*xx**2+d*xx**3
 f=max(1,ceiling(ff))

 end subroutine gaussfit_mpi_set_weight
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_remove_item
!! NAME
!!  gaussfit_mpi_remove_item
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_remove_item(iterm,pload)

!Arguments ------------------------------------
 integer,intent(in)::iterm
 integer,intent(inout)::pload

!Local variables ------------------------------
 integer:: f_i

!***********************************************************************

 call gaussfit_mpi_set_weight(f_i,iterm)
 pload=pload-f_i

 end subroutine gaussfit_mpi_remove_item
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_add_item
!! NAME
!!  gaussfit_mpi_add_item
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_add_item(iterm,pload)

!Arguments ------------------------------------
 integer,intent(in)::iterm
 integer,intent(inout)::pload

!Local variables ------------------------------
 integer:: f_i

!************************************************************************

 call gaussfit_mpi_set_weight(f_i,iterm)
 pload=pload+f_i

 end subroutine gaussfit_mpi_add_item
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_calc_deviation
!! NAME
!!  gaussfit_mpi_calc_deviation
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_calc_deviation(deviation,nproc,proc_load)

!Arguments ------------------------------------
 integer,intent(in)::nproc
 integer,intent(in)::proc_load(nproc)
 integer,intent(out)::deviation

!Local variables ------------------------------
 integer:: jproc,kproc,kload,jload

!************************************************************************

! deviation=0
! do jproc=1,nproc
!   deviation=deviation+abs(proc_load(jproc)-ideal)
! end do

 deviation=0
 do jproc=1,nproc
   jload=proc_load(jproc)
   do kproc=1,jproc
     kload=proc_load(kproc)
     deviation=deviation+abs(kload-jload)
   end do
 end do

 end subroutine gaussfit_mpi_calc_deviation
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_swap
!! NAME
!!  gaussfit_mpi_swap
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_swap(iterm,jterm,&
&           nproc,nterm_bounds,proc_dist,proc_load)

!Arguments ------------------------------------
 integer,intent(in)::iterm,jterm,nproc,nterm_bounds(2)
 integer,intent(inout)::proc_dist(nterm_bounds(1):nterm_bounds(2)),proc_load(nproc)

!Local variables ------------------------------
 integer:: deviation1,deviation2
 integer:: iproc,jproc

!************************************************************************

!Calculate initial state
 call gaussfit_mpi_calc_deviation(deviation1,nproc,proc_load)
 iproc=proc_dist(iterm)
 jproc=proc_dist(jterm)
!Swap terms:
 call gaussfit_mpi_add_item(jterm,proc_load(iproc))
 call gaussfit_mpi_remove_item(iterm,proc_load(iproc))
 call gaussfit_mpi_add_item(iterm,proc_load(jproc))
 call gaussfit_mpi_remove_item(jterm,proc_load(jproc))
!Calculate final state
 call gaussfit_mpi_calc_deviation(deviation2,nproc,proc_load)
!Swap them only if final state is better than the initial one
 if(deviation2<deviation1) then
   proc_dist(iterm)=jproc
   proc_dist(jterm)=iproc
 else
!  Return work load to initial state
   call gaussfit_mpi_add_item(iterm,proc_load(iproc))
   call gaussfit_mpi_remove_item(jterm,proc_load(iproc))
   call gaussfit_mpi_add_item(jterm,proc_load(jproc))
   call gaussfit_mpi_remove_item(iterm,proc_load(jproc))
!   write(*,*)'Back initial state'
!   write(*,*)'proc_load',proc_load(:)
 end if

 end subroutine gaussfit_mpi_swap
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_assign
!! NAME
!!  gaussfit_mpi_assign
!!
!! FUNCTION
!!  Set task to a processor
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_assign(iterm,nproc,nterm_bounds,&
& proc_dist,proc_load)

!Arguments ------------------------------------
 integer,intent(in)::iterm,nproc,nterm_bounds(2)
 integer,intent(inout)::proc_dist(nterm_bounds(1):nterm_bounds(2)),proc_load(nproc)

!Local variables ------------------------------
 integer:: iproc,jproc,dev
 integer:: deviation(nproc),mindev
 character(len=100) :: msg

!************************************************************************

 do iproc=1,nproc
  !add this term to iproc
  call gaussfit_mpi_add_item(iterm,proc_load(iproc))
  !calculate the deviation for this configuration
  call gaussfit_mpi_calc_deviation(dev,nproc,proc_load)
  deviation(iproc)=dev
  !remove this term from iproc
  call gaussfit_mpi_remove_item(iterm,proc_load(iproc))
 end do

!assign to jproc, the proc with minimal deviation above:
 jproc=-1; mindev=999999999
 do iproc=1,nproc
   if(deviation(iproc)<mindev) then
     mindev=deviation(iproc)
     jproc=iproc
   end if
 end do
 if(jproc==-1) then
!  One should not get here!
   msg = 'error in accomodate_mpi'
   MSG_BUG(msg)
 end if

!assign this term for jproc
 proc_dist(iterm)=jproc
 call gaussfit_mpi_add_item(iterm,proc_load(jproc))

 end subroutine gaussfit_mpi_assign
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_mpi_main
!! NAME
!!  gaussfit_mpi_main
!!
!! FUNCTION
!!  Set charge for each processor
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_mpi_main(nproc,nterm_bounds,proc_dist)

!Arguments ------------------------------------
 integer,intent(in)::nproc,nterm_bounds(2)
 integer,intent(out)::proc_dist(nterm_bounds(1):nterm_bounds(2))

!Local variables ------------------------------
 integer:: dev1,dev2,ii
 integer:: iproc,iterm,jterm,ngauss,weight
 integer:: proc_load(nproc)
 character(len=500) :: msg

!************************************************************************

 proc_load=0; proc_dist=0 !initializations

!1) get a first-trial distribution:
 do iterm=nterm_bounds(2),nterm_bounds(1),-1
     call gaussfit_mpi_assign(iterm,nproc,nterm_bounds,proc_dist,proc_load)
 end do

!Do the following 20 times
 do ii=1,20
!  Calculate initial state
   call gaussfit_mpi_calc_deviation(dev1,nproc,proc_load)
!  Try to swap tasks between two processors:
   do iterm=nterm_bounds(2),nterm_bounds(1),-1
     do jterm=nterm_bounds(2),nterm_bounds(1),-1
       call gaussfit_mpi_swap(iterm,jterm,&
&       nproc,nterm_bounds,proc_dist,proc_load)
     end do
   end do
!!  Try to reassign tasks to different processors
   do iterm=nterm_bounds(2),nterm_bounds(1),-1
     iproc=proc_dist(iterm)
!    Remove this job from this node
     call gaussfit_mpi_remove_item(iterm,proc_load(iproc))
!    Accomodate this again:
     call gaussfit_mpi_assign(iterm,nproc,nterm_bounds,proc_dist,proc_load)
   end do
!  Calculate final state
   call gaussfit_mpi_calc_deviation(dev2,nproc,proc_load)
!  If final state equals initial state, exit:
   if(dev2 == dev1) exit
!   write(*,'(a)')'Deviation: ',dev2
 end do

!Write down distribution:
 write(msg,'(3a)') 'MPI distribution',ch10,'N. gauss, iproc, weight '
 call wrtout(std_out,msg,'COLL')
 do iterm=nterm_bounds(2),nterm_bounds(1),-1
     ngauss=iterm*2
     call gaussfit_mpi_set_weight(weight,iterm)
     write(msg,'(3(i4,1x))') ngauss,proc_dist(iterm),weight
     call wrtout(std_out,msg,'COLL')
 end do
 write(msg,'(a)') 'Load per processor: '
 call wrtout(std_out,msg,'COLL')
 do iproc=1,nproc
   write(msg,'(i5,1x,i10)') iproc,proc_load(iproc)
   call wrtout(std_out,msg,'COLL')
 end do

 end subroutine gaussfit_mpi_main
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_fit
!! NAME
!!  gaussfit_fit
!!
!! FUNCTION
!!  Fits a given input function f(r) to a sum of gaussians
!!
!! INPUTS
!!  chisq= It measures how good is the fitting.
!!   It is defined here as sum_i^N_i f(x_i)-y(x_i)/N_i
!!  constrains= constraints for Gaussians
!!  limit(nparam)= it limits the widths of Gaussians
!!  maxiter=maximum number of iterations
!!  nparam= number of parameters (a constant times nterm)
!!         if(option==1)nparam=nterm*4
!!         if(option==2)nparam=nterm*6
!!         if(option==3)nparam=nterm*2
!!         if(option==4)nparam=nterm*4
!!  nterm= number of Gaussians
!!  nx=number of points along the x axis
!!  option=1 fit to a1 cos(a2 x^2)+ a3 sin( a4 x^2)
!!         2 fit to a1 exp(-a2 x^2)*(a3 cos (a4 x^2) + a5 sin (a6 x^2) )
!!         3 fit to a1 cos (k x^2) + a2 sin (k x^2)
!!         4 fit to a1 exp(-a2 x^2)* (a3 cos(k x^2)+ a4 sin (k x^2))
!!  outfile= filename for output (only written if verbosity>1)
!!  verbosity= controls output volume
!!  weight(nparam)= weights for the fitting procedure
!!  x(nx)= points along the x axis
!!  y(nx)= function to be fitted
!!  rpaw ,optional= paw radius
!!
!! OUTPUT
!! y_out(nx)= fitted function
!!
!! SIDE EFFECTS
!! if(verbosity>1) output files are written with y(x) and y_out(x)
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_fit(chisq,constrains,&
& limit,maxiter,nparam,nterm,nx,option,outfile,param,&
& verbosity,weight,x,y,y_out)

!Arguments ------------------------------------
 integer, intent(in) :: maxiter,nparam
 integer,intent(in)  :: nterm,nx,option,verbosity
 !real(dp),optional,intent(in)::rpaw
 !arrays
 integer,intent(in)::constrains(nparam)
 real(dp),intent(in)::limit(nparam),weight(nparam)
 real(dp),intent(in)::x(nx),y(nx)
 real(dp),intent(inout)::param(nparam)
 real(dp),intent(out)::chisq,y_out(nx)
 character(80),intent(in)::outfile

!Local variables ------------------------------
 integer, parameter :: wfn_unit=1007
 integer::ix
 real(dp)::rerror
 real(dp),allocatable::sy(:)

! *************************************************************************

 LIBPAW_ALLOCATE(sy,(nx))

 sy(:)=1.0d0

 call gaussfit_rlsf(&
& chisq,constrains,limit,maxiter,&
& nterm,nparam,nx,option,param(1:nparam),&
& verbosity,weight,x,y)
!
 if(verbosity>=1) then
   if(option==1) then
     call gaussfit_calc_deriv_c2(nparam,nterm,nx,1,param,x,y_out)
   elseif(option==2) then
     call gaussfit_calc_deriv_c(nparam,nterm,nx,1,param,x,y_out)
   elseif(option==3) then
     call gaussfit_calc_deriv_c3(nparam,nterm,nx,1,param,x,y_out)
   elseif(option==4) then
     call gaussfit_calc_deriv_c4(nparam,nterm,nx,1,param,x,y_out)
   end if
!
!
   open(wfn_unit,file=outfile,form='formatted',status='unknown')
!  per_error=0.d0
   do ix=1, nx
     rerror=abs(y(ix)-y_out(ix))
     write(wfn_unit,'(6(e20.12,1x))')x(ix),y(ix),y_out(ix),rerror
   end do
   close(wfn_unit)

 end if

 LIBPAW_DEALLOCATE(sy)

end subroutine gaussfit_fit
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_calc_deriv_r
!! NAME
!!  gaussfit_calc_deriv_r
!!
!! FUNCTION
!!  Calculate derivatives for Gaussians
!!   Only relevant for fitting Gaussians algorithm.
!!   The Gaussians expressions are defined in the comments of "gaussfit_main"
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_calc_deriv_r(nterm,nparam,nx,opt,param,x,y_out,&
& deriv) ! optional

!Arguments -------------------------------
 integer,intent(in)::nx     !number of point in the x grid
 integer,intent(in)::nparam !number of parameters
 integer,intent(in)::nterm !number of gaussian expressions
 integer,intent(in)::opt    !option:
                            !1) calculate only f(x)
                            !2) calculate f(x) and its derivatives
 real(dp),intent(in)::param(nparam) !parameters
 real(dp),intent(in)::x(nx) !xgrid
 real(dp),intent(out)::y_out(nx) !f(x)
 real(dp),optional,intent(out)::deriv(nx,nparam) !derivatives

!Local variables-------------------------------
 integer::iexp,ii
 real(dp)::alpha1(nterm),alpha2(nterm),alpha3(nterm)
 real(dp)::term1(nx,nterm)
 real(dp)::aux1(nx)
 !real(dp)::step

! *********************************************************************

!
!Initialize
!
 y_out(:)=0.d0
!
!Get parameters from parameters array:
!
 alpha1(:)=param(1:nterm)
 alpha2(:)=param(nterm+1:2*nterm)
 alpha3(:)=param(2*nterm+1:3*nterm)
!
!alpha3
!set to constant values of x
!step=rpaw/real(nterm,dp)
!do ii=1,nterm
!raux=step*real(ii,dp)
!alpha3(ii)=raux
!end do
!
!
!
!calculate useful quantities
!
 do iexp=1,nterm
   aux1(:)=-alpha2(iexp)*(x(:)-alpha3(iexp))**2
   term1(:,iexp)=alpha1(iexp)*exp(aux1(:))
 end do
!
 do iexp=1,nterm
   y_out(:)=y_out(:)+term1(:,iexp)
 end do
!
!Calculate derivatives:
!
 if(opt==2) then
!
!  alpha1
!
   do iexp=1,nterm
     aux1(:)=term1(:,iexp)/alpha1(iexp)
     deriv(:,iexp)=aux1(:)
   end do
!
!  alpha2
!
   do iexp=1,nterm
     ii=nterm+iexp
     aux1(:)=-term1(:,iexp)*(x(:)-alpha3(iexp))
     deriv(:,ii)=aux1(:)
     deriv(:,ii)=0.1d0
   end do
!
!  alpha3
!
   do iexp=1,nterm
     ii=2*nterm+iexp
     aux1(:)=term1(:,iexp)*2.d0*alpha2(iexp)
     aux1(:)=aux1(:)*(x(:)-alpha3(iexp))
     deriv(:,ii)=aux1(:)
   end do
 end if

end subroutine gaussfit_calc_deriv_r
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_calc_deriv_c3
!! NAME
!!  gaussfit_calc_deriv_c3
!!
!! FUNCTION
!!  Calculate expressions and derivatives for Gaussians fitting.
!!   The Gaussians expressions are defined in the comments of "gaussfit_main"
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_calc_deriv_c3(nparam,nterm,nx,opt,param,x,y_out,&
& deriv) ! optional

!Arguments -------------------------------
 integer,intent(in)::nparam !number of parameters
 integer,intent(in)::nterm !number of gaussian expressions
 integer,intent(in)::nx     !number of point in the x grid
 integer,intent(in)::opt    !option:
                            !1) calculate only f(x)
                            !2) calculate f(x) and its derivatives
 real(dp),intent(in)::param(nparam) !parameters
 real(dp),intent(in)::x(nx) !xgrid
 real(dp),intent(out)::y_out(nx) !f(x)
 real(dp),optional,intent(out)::deriv(nx,nparam) !derivatives

!Local variables-------------------------------
 integer::iexp,ii
 real(dp)::sep
 real(dp)::alpha1(nterm),alpha2(nterm),alpha3(nterm)
 real(dp)::term1(nx,nterm)
 real(dp)::sin1(nx,nterm),cos1(nx,nterm)
 real(dp)::aux1(nx),aux2(nx)

! *********************************************************************

!
!Initialize
!
 y_out(:)=0.d0
!
 sep=1.2d0
!
!Get param from param array:
!
 alpha1(:)=param(1:nterm)
 alpha2(:)=param(nterm+1:2*nterm)
!
 do ii=1,nterm
   alpha3(ii)=sep**(ii)
 end do
!
!calculate useful quantities
!
 do iexp=1,nterm
   aux1(:)=alpha3(iexp)*x(:)**2
!
   sin1(:,iexp)=sin(aux1(:))
   cos1(:,iexp)=cos(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha1(iexp)*sin1(:,iexp)
   aux2(:)=alpha2(iexp)*cos1(:,iexp)
   term1(:,iexp)=aux1(:)+aux2(:)
   y_out(:)=y_out(:)+term1(:,iexp)
 end do
!
!Calculate derivatives:
!
 if(opt==2) then
!
!  alpha1
!
   do iexp=1,nterm
     deriv(:,iexp)=sin1(:,iexp)
   end do
!
!  alpha2
!
   do iexp=1,nterm
     ii=nterm+iexp
     deriv(:,ii)=cos1(:,iexp)
   end do
 end if

end subroutine gaussfit_calc_deriv_c3
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_calc_deriv_c2
!! NAME
!!  gaussfit_calc_deriv_c2
!!
!! FUNCTION
!!  Calculate expressions and derivatives for Gaussians fitting.
!!   The Gaussians expressions are defined in the comments of "gaussfit_main"
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_calc_deriv_c2(nparam,nterm,nx,opt,param,x,y_out,&
& deriv) ! optional

!Arguments -------------------------------
 integer,intent(in)::nparam !number of param
 integer,intent(in)::nterm  !number of gaussian expressions
 integer,intent(in)::nx     !number of point in the x grid
 integer,intent(in)::opt    !option:
                            !1) calculate only f(x)
                            !2) calculate f(x) and its derivatives
 real(dp),intent(in)::param(nparam) !parameters
 real(dp),intent(in)::x(nx) !xgrid
 real(dp),intent(out)::y_out(nx) !f(x)
 real(dp),optional,intent(out)::deriv(nx,nparam) !derivatives

!Local variables-------------------------------
 integer::iexp,ii
 real(dp)::alpha1(nterm),alpha2(nterm),alpha3(nterm)
 real(dp)::alpha4(nterm)
 real(dp)::term1(nx,nterm)
 real(dp)::sin1(nx,nterm),sin2(nx,nterm),cos1(nx,nterm),cos2(nx,nterm)
 real(dp)::aux1(nx),aux2(nx)

! *********************************************************************

!
!Initialize
!
 y_out(:)=0.d0
!
!Get param from param array:
!
 alpha1(:)=param(1:nterm)
 alpha2(:)=param(nterm+1:2*nterm)
 alpha3(:)=param(2*nterm+1:3*nterm)
 alpha4(:)=param(3*nterm+1:4*nterm)
!
!calculate useful quantities
!
!
 do iexp=1,nterm
   aux1(:)=alpha2(iexp)*x(:)**2
   sin1(:,iexp)=sin(aux1(:))
!
   aux1(:)=alpha4(iexp)*x(:)**2
   sin2(:,iexp)=sin(aux1(:))
!
   aux1(:)=alpha2(iexp)*x(:)**2
   cos1(:,iexp)=cos(aux1(:))
!
   aux1(:)=alpha4(iexp)*x(:)**2
   cos2(:,iexp)=cos(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha1(iexp)*sin1(:,iexp)
   aux2(:)=alpha3(iexp)*cos2(:,iexp)
   term1(:,iexp)=aux1(:)+aux2(:)
   y_out(:)=y_out(:)+term1(:,iexp)
 end do
!
!Calculate derivatives:
!
 if(opt==2) then
!
!  alpha1
!
   do iexp=1,nterm
     deriv(:,iexp)=sin1(:,iexp)
   end do
!
!  alpha2
!
   do iexp=1,nterm
     ii=nterm+iexp
     aux1(:)=alpha1(iexp)*cos1(:,iexp)*x(:)**2
     deriv(:,ii)=aux1(:)
   end do
!
!  alpha3
!
   do iexp=1,nterm
     ii=2*nterm+iexp
     deriv(:,ii)=cos2(:,iexp)
   end do
!
!  alpha4
!
   do iexp=1,nterm
     ii=3*nterm+iexp
     aux1(:)=-alpha3(iexp)*sin2(:,iexp)*x(:)**2
     deriv(:,ii)=aux1(:)
   end do
 end if

end subroutine gaussfit_calc_deriv_c2
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_calc_deriv_c
!! NAME
!!  gaussfit_calc_deriv_c
!!
!! FUNCTION
!!  Calculate expressions and derivatives for Gaussians fitting.
!!   The Gaussians expressions are defined in the comments of "gaussfit_main"
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_calc_deriv_c(nparam,nterm,nx,opt,param,x,y_out,&
& deriv) ! optional

!Arguments -------------------------------
 integer,intent(in)::nparam !number of parameters
 integer,intent(in)::nterm !number of gaussian expressions
 integer,intent(in)::nx     !number of point in the x grid
 integer,intent(in)::opt    !option:
                            !1) calculate only f(x)
                            !2) calculate f(x) and its derivatives
 real(dp),intent(in)::param(nparam) !parameters
 real(dp),intent(in)::x(nx) !xgrid
 real(dp),intent(out)::y_out(nx) !f(x)
 real(dp),optional,intent(out)::deriv(nx,nparam) !derivatives

!Local variables-------------------------------
 integer::iexp,ii
 real(dp)::alpha1(nterm),alpha2(nterm),alpha3(nterm)
 real(dp)::alpha4(nterm),alpha5(nterm),alpha6(nterm)
 real(dp)::aux1(nx),aux2(nx)
 real(dp)::cos1(nx,nterm),cos2(nx,nterm),sin1(nx,nterm),sin2(nx,nterm)
 real(dp)::term1(nx,nterm),term2(nx,nterm)

! *********************************************************************

!
!Initialize
!
 y_out(:)=0.d0
!
!Get parameters from param array:
!
 alpha1(:)=param(1:nterm)
 alpha2(:)=param(nterm+1:2*nterm)
 alpha3(:)=param(2*nterm+1:3*nterm)
 alpha4(:)=param(3*nterm+1:4*nterm)
 alpha5(:)=param(4*nterm+1:5*nterm)
 alpha6(:)=param(5*nterm+1:6*nterm)
!
!calculate useful quantities
!
 do iexp=1,nterm
   aux1(:)=-alpha2(iexp)*x(:)**2
   term1(:,iexp)=alpha1(iexp)*exp(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha4(iexp)*x(:)**2
   sin1(:,iexp)=sin(aux1(:))
!
   aux1(:)=alpha6(iexp)*x(:)**2
   sin2(:,iexp)=sin(aux1(:))
!
   aux1(:)=alpha4(iexp)*x(:)**2
   cos1(:,iexp)=cos(aux1(:))
!
   aux1(:)=alpha6(iexp)*x(:)**2
   cos2(:,iexp)=cos(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha3(iexp)*sin1(:,iexp)
   aux2(:)=alpha5(iexp)*cos2(:,iexp)
   term2(:,iexp)=aux1(:)+aux2(:)
   y_out(:)=y_out(:)+term1(:,iexp)*term2(:,iexp)
 end do
!
!Calculate derivatives:
!
 if(opt==2) then
!
!  alpha1
!
   do iexp=1,nterm
     aux1(:)=term1(:,iexp)/alpha1(iexp)
     aux2(:)=aux1(:)*term2(:,iexp)
     deriv(:,iexp)=aux2(:)
   end do
!
!  alpha2
!
   do iexp=1,nterm
     ii=nterm+iexp
     aux1(:)=-term1(:,iexp)*term2(:,iexp)
     aux2(:)=aux1(:)*x(:)**2
     deriv(:,ii)=aux2(:)
   end do
!
!  alpha3
!
   do iexp=1,nterm
     ii=2*nterm+iexp
     aux1(:)=term1(:,iexp)*sin1(:,iexp)
     deriv(:,ii)=aux1(:)
   end do
!
!  alpha4
!
   do iexp=1,nterm
     ii=3*nterm+iexp
     aux1(:)=term1(:,iexp)*alpha3(iexp)
     aux2(:)=cos1(:,iexp)*x(:)**2
     deriv(:,ii)=aux2(:)*aux1(:)
   end do
!
!  alpha5
!
   do iexp=1,nterm
     ii=4*nterm+iexp
     aux1(:)=term1(:,iexp)*cos2(:,iexp)
     deriv(:,ii)=aux1(:)
   end do
!
!  alpha6
!
   do iexp=1,nterm
     ii=5*nterm+iexp
     aux1(:)=-term1(:,iexp)*alpha5(iexp)
     aux2(:)=sin2(:,iexp)*x(:)**2
     deriv(:,ii)=aux1(:)*aux2(:)
   end do
 end if

end subroutine gaussfit_calc_deriv_c
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_calc_deriv_c4
!! NAME
!!  gaussfit_calc_deriv_c4
!!
!! FUNCTION
!!  Calculate expressions and derivatives for Gaussians fitting.
!!   The Gaussians expressions are defined in the comments of "gaussfit_main"
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_calc_deriv_c4(nparam,nterm,nx,opt,param,x,y_out,&
& deriv) ! optional

!Arguments -------------------------------
 integer,intent(in)::nparam !number of parameters
 integer,intent(in)::nterm !number of gaussian expressions
 integer,intent(in)::nx     !number of point in the x grid
 integer,intent(in)::opt    !option:
                            !1) calculate only f(x)
                            !2) calculate f(x) and its derivatives
 real(dp),intent(in)::param(nparam) !parameters
 real(dp),intent(in)::x(nx) !xgrid
 real(dp),intent(out)::y_out(nx) !f(x)
 real(dp),optional,intent(out)::deriv(nx,nparam) !derivatives

!Local variables-------------------------------
 integer::iexp,ii
 real(dp)::raux,sep
 real(dp)::alpha1(nterm),alpha2(nterm),alpha3(nterm)
 real(dp)::alpha4(nterm),alpha5(nterm)
 real(dp)::aux1(nx),aux2(nx)
 real(dp)::cos1(nx,nterm),sin1(nx,nterm)
 real(dp)::term1(nx,nterm),term2(nx,nterm)

! *********************************************************************

!
!Initialize
!
 sep=1.1d0
 y_out(:)=0.d0

!Get parameters from param array:
!
 alpha1(:)=param(1:nterm)
 alpha2(:)=param(nterm+1:2*nterm)
 alpha3(:)=param(2*nterm+1:3*nterm)
 alpha4(:)=param(3*nterm+1:4*nterm)
!
!
 raux=(2.d0*pi)/real(nterm,dp)
 do ii=1,nterm
   alpha5(ii)=sep**(ii)
!  alpha5(ii)=raux*real(ii-1,dp)
 end do
!
!calculate useful quantities
!
 do iexp=1,nterm
   aux1(:)=-alpha2(iexp)*x(:)**2
   term1(:,iexp)=alpha1(iexp)*exp(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha5(iexp)*x(:)**2
!
   sin1(:,iexp)=sin(aux1(:))
   cos1(:,iexp)=cos(aux1(:))
 end do
!
 do iexp=1,nterm
   aux1(:)=alpha3(iexp)*sin1(:,iexp)
   aux2(:)=alpha4(iexp)*cos1(:,iexp)
   term2(:,iexp)=aux1(:)+aux2(:)
   y_out(:)=y_out(:)+term1(:,iexp)*term2(:,iexp)
 end do
!
!Calculate derivatives:
!
 if(opt==2) then
!
!  alpha1
!
   do iexp=1,nterm
     aux1(:)=term1(:,iexp)/alpha1(iexp)
     aux2(:)=aux1(:)*term2(:,iexp)
     deriv(:,iexp)=aux2(:)
   end do
!
!  alpha2
!
   do iexp=1,nterm
     ii=nterm+iexp
     aux1(:)=-term1(:,iexp)*term2(:,iexp)
     aux2(:)=aux1(:)*x(:)**2
     deriv(:,ii)=aux2(:)
   end do
!
!  alpha3
!
   do iexp=1,nterm
     ii=2*nterm+iexp
     aux1(:)=term1(:,iexp)*sin1(:,iexp)
     deriv(:,ii)=aux1(:)
   end do
!
!  alpha4
!
   do iexp=1,nterm
     ii=3*nterm+iexp
     aux1(:)=term1(:,iexp)*cos1(:,iexp)
     deriv(:,ii)=aux1(:)
   end do
 end if

end subroutine gaussfit_calc_deriv_c4
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_rlsf
!! NAME
!!  gaussfit_rlsf
!!
!! FUNCTION
!!  Fits a given function to a sum of Gaussians.
!!  Uses the Levenberg-Marquardt algorithm.
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2020 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  The original Levemberg Marquardt routines were written by Armando Sole
!!  These were modified for the ARPUS spectra in the BigDFT code by A. Mirone.
!!  These were re-writen in Fortran and further modified in ABINIT for our particular needs.
!!
!! INPUTS
!!  option=1 fit to a1 cos(a2 x^2)+ a3 sin( a4 x^2)
!!         2 fit to a1 exp(-a2 x^2)*(a3 cos (a4 x^2) + a5 sin (a6 x^2) )
!!         3 fit to a1 cos (k x^2) + a2 sin (k x^2)
!!         4 fit to a1 exp(-a2 x^2)* (a3 cos(k x^2)+ a4 sin (k x^2))
!!  if(option==1)mparam=nterm_bounds(2)*4
!!  if(option==2)mparam=nterm_bounds(2)*6
!!  if(option==3)mparam=nterm_bounds(2)*2
!!  if(option==4)mparam=nterm_bounds(2)*4
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_rlsf(&
&chisq,constrains,limit,maxiter,&
&nterm,nparam,nx,option,parameters,&
&verbosity,weight,x,y)

!Arguments -------------------------------
 real(dp),parameter::deltachi=tol10
 integer, intent(in) ::maxiter,nparam,nterm,nx
 integer, intent(in) ::option ,verbosity
 integer, intent(in) ::constrains(nparam)
 real(dp),intent(out)::chisq
 !arrays
 real(dp),intent(in)::limit(nparam),weight(nparam)
 real(dp),intent(inout)::parameters(nparam)
 real(dp),intent(in)::x(nx),y(nx)

!Local variables-------------------------------
 integer::flag,ii,info,iter,jj,niter
 real(dp):: deltax
 real(dp)::chisq0,flambda,eta,lastdeltachi
 integer::ipvt(nparam)
 real(dp)::alpha(nparam,nparam)
 real(dp)::alpha0(nparam,nparam),beta(nparam)
 real(dp)::deltapar(nparam)
 real(dp)::tmp1(nparam,nparam)
 real(dp)::work(nparam)
 real(dp)::workpar(nparam)
 real(dp)::yfit(nx)
 character(len=500) :: msg

! *********************************************************************

!
 !flambda=1e-6
 flambda=1e-7
!iter=maxiter !later it is changed
 niter=0
 deltax=x(2)-x(1) !we assume this is a linear grid
!
 iter_loop: do iter=1,maxiter
!
   call gaussfit_chisq_alpha_beta(alpha0,beta,chisq0,&
&   nparam,nterm,nx,option,parameters,x,y)
!
   flag=0
   lastdeltachi=chisq0
!
!
   while_flag: do
     if(flag .ne. 0) exit while_flag
!
     tmp1=0.d0
     do ii=1,nparam
       tmp1(ii,ii)=1.d0*flambda  !identity matrix * flambda
     end do
     alpha=alpha0+tmp1*alpha0
!    Invert alpha matrix
     tmp1=alpha
     call dgetrf(nparam,nparam,tmp1,nparam,ipvt,info)
     if (.not.info==0) then
       if(verbosity>1) then
         write(msg,'(a)')'Matrix is singular'
         call wrtout(std_out,msg,'COLL')
       end if
       chisq=-1.d0
       exit iter_loop
     end if
     call dgetri(nparam,tmp1,nparam,ipvt,work,nparam,info)
     deltapar=0.d0
     if (.not.info==0) then
       if(verbosity>2) then
         write(msg,'(a)')'Matrix is singular'
         call wrtout(std_out,msg,'COLL')
       end if
       chisq=-1.d0
       exit iter_loop
     end if
!
     if(tmp1(1,1) .ne. tmp1(1,1)) then !If is NaN
       chisq=-1.d0
       exit iter_loop
     end if
     if(abs(tmp1(1,1)) == tmp1(1,1)*tmp1(1,1)) then !If is infinity
       chisq=-1.d0
       exit iter_loop
     end if
!
     do ii=1,nparam
       do jj=1,nparam
         deltapar(ii)=deltapar(ii)+beta(jj)*tmp1(jj,ii)
       end do
     end do
!    apply constrains
     workpar(1:nparam)=parameters(1:nparam)+deltapar(1:nparam)*weight(1:nparam)
     call gaussfit_apply_constrains(constrains,limit,nparam,workpar)
!
     if(option==1) then
       call gaussfit_calc_deriv_c2(nparam,nterm,nx,1,workpar,x,yfit)
     elseif(option==2) then
       call gaussfit_calc_deriv_c(nparam,nterm,nx,1,workpar,x,yfit)
     elseif(option==3) then
       call gaussfit_calc_deriv_c3(nparam,nterm,nx,1,workpar,x,yfit)
     elseif(option==4) then
       call gaussfit_calc_deriv_c4(nparam,nterm,nx,1,workpar,x,yfit)
     end if
     chisq=0.d0
     do ii=1,nx
       chisq=chisq + ((y(ii)-yfit(ii)))**2
     end do
     chisq=chisq*deltax
!
!    write(*,'("chisq ",f12.5,"  chisq0 ",f12.5)')chisq,chisq0
!
     if(chisq > chisq0) then
       flambda=flambda*2.0d0
       if( flambda > 1000.d0) then
         flag=1
!        iter=0
         if(verbosity>2) then
           write(msg,'(a)')'flambda > 1000.d0'
           call wrtout(std_out,msg,'COLL')
         end if
         exit iter_loop
       end if
     else
       flag=1
       parameters=workpar
       eta=0.d0
       lastdeltachi=(chisq0-chisq) !/(chisq0+eta)
       if(lastdeltachi<deltachi) cycle
       chisq0=chisq
       flambda=flambda/2.d0
       if(verbosity>2) then
         write(msg,'("iter = ",i4," chisq = ",e15.6)')iter,chisq
         call wrtout(std_out,msg,'COLL')
       end if
     end if
   end do while_flag
 end do iter_loop

end subroutine gaussfit_rlsf
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_chisq_alpha_beta
!! NAME
!!  gaussfit_chisq_alpha_beta
!!
!! FUNCTION
!!  Finds chisq, alpha and beta parameters for LSF using the Levenberg-Marquardt algorithm.
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2020 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  The original Levemberg Marquardt routines were written by Armando Sole
!!  These were modified for the ARPUS spectra in the BigDFT code by A. Mirone.
!!  These were re-writen in Fortran and further modified in ABINIT for our particular needs.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_chisq_alpha_beta(alpha,beta,chisq,&
& nparam,nterm,nx,option,parameters,x,y)

!Arguments -------------------------------
 integer,intent(in)::nparam,nterm,nx
 integer,intent(in)::option
 real(dp),intent(out)::chisq
 real(dp),intent(in)::parameters(nparam),x(nx),y(nx)
 real(dp),intent(out)::alpha(nparam,nparam),beta(nparam)

!Local variables-------------------------------
 integer::ii,jj,kk
 real(dp)::deltax,help1
 !arrays
 real(dp)::deltay(nx),deriv(nx,nparam),derivi(nx)
 real(dp)::yfit(nx)
 real(dp)::help0(nx),help2(nx),help3(nparam)

! *********************************************************************

 deltax=x(2)-x(1) !we assume a linear grid
!
 if(option==1) then
   call gaussfit_calc_deriv_c2(nparam,nterm,nx,2,parameters,x,yfit,deriv)
 elseif(option==2) then
   call gaussfit_calc_deriv_c(nparam,nterm,nx,2,parameters,x,yfit,deriv)
 elseif(option==3) then
   call gaussfit_calc_deriv_c3(nparam,nterm,nx,2,parameters,x,yfit,deriv)
 elseif(option==4) then
   call gaussfit_calc_deriv_c4(nparam,nterm,nx,2,parameters,x,yfit,deriv)
 end if
 deltay=y-yfit
 help0=deltay
!
 do ii=1,nparam
   derivi(:)=deriv(:,ii)
   help1=0.d0
   do jj=1,nx
     help1=help1+help0(jj)*derivi(jj)
   end do
   beta(ii)=help1
!  help1 = innerproduct(deriv,weight*derivi)
!  below I use help3 instead for the array dimenstions
   help3=0.d0
   do kk=1,nparam
     do jj=1,nx
       help3(kk)=help3(kk)+deriv(jj,kk)*derivi(jj)
     end do
   end do
!  !
   alpha(:,ii)=help3(:)
 end do
!
 help2(:)=help0(:)*deltay(:)
 chisq=sum(help2)
 chisq=chisq*deltax

end subroutine gaussfit_chisq_alpha_beta
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_set_param1
!! NAME
!!  gaussfit_set_param1
!!
!! FUNCTION
!!  Sets parameters for LSF
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_set_param1(nterm,nparam,nx,param,sep,x,y)

!Arguments -------------------------------
 integer,intent(in)::nterm,nparam,nx
 real(dp),intent(in)::sep
 real(dp),intent(in)::x(nx),y(nx)
 real(dp),intent(out)::param(nparam)

!Local variables-------------------------------
 integer::ii,jj
 real(dp)::raux

! *********************************************************************

!
 param(:)=1.0d0
!exps=1.0/(x(nx)**2)
!
!alpha1

!raux=maxval( y(:),nx )
!maxval gives problems in some architectures:
 raux=-9999999
 do ii=1,nx
   if(raux<y(ii)) raux=y(ii)
 end do
 param(1:nterm)=raux
!
!alpha2
!
!y(r_c)=e^{-\alpha1 r_c}
 raux=-log(abs(y(nx))+tol10)/(x(nx)**2)
 param(nterm+1:nterm+2)=raux
!
!raux=0.5d0*pi/real(nterm,dp)
 do jj=1,nterm
   ii=jj+3*nterm
   param(ii)=sep**(jj)
!  param(ii)=raux*real(jj,dp)
 end do
!
 do jj=1,nterm
   ii=jj+5*nterm
   param(ii)=sep**(jj)
!  param(ii)=raux*real(jj,dp)
 end do

end subroutine gaussfit_set_param1
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_set_param2
!! NAME
!!  gaussfit_set_param2
!!
!! FUNCTION
!!  Sets parameters for LSF
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_set_param2(nterm,nparam,nx,param,rpaw,x,y)

!Arguments -------------------------------
 integer,intent(in)::nterm,nparam,nx
 real(dp),intent(in)::rpaw
 real(dp),intent(in)::x(nx),y(nx)
 real(dp),intent(out)::param(nparam)

!Local variables-------------------------------
 integer::ii,jj
 real(dp)::exps,raux,sig,step

! *************************************************************************

 step=rpaw/real(nterm-1,dp)
!exps=1.0/(rpaw/(real(nterm,dp)/2.d0))**2
 exps=1.0/(step*1.0d0)**2
!alpha2 (width of gaussians)
!Set to exps*real(nterm,dp)**2, for a good guess
 param(nterm+1:2*nterm)=exps !*real(nterm,dp)**2
!
 do jj=1,nterm
!  alpha3
!  set to constant values of x
   ii=jj+2*nterm
   raux=step*real(jj-1,dp)
   param(ii)=raux
!  alpha1
!  set to the value of y at that point
   call gaussfit_param2_findsign()
   param(jj)=sig
 end do
!
 contains
!!***

!!****f* gaussfit_set_param2/gaussfit_param2_findsign
!! NAME
!!  gaussfit_param2_findsign
!!
!! FUNCTION
!!  Finds the value of y at a given point
!!  This was taken out of gaussfit_set_param2 to make the code more
!!  readable.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

 subroutine gaussfit_param2_findsign()

!Arguments -------------------------------
!Local variables-------------------------------
 integer::ix,minx
 real(dp)::dist,mindist,xx,yy

! *********************************************************************

   mindist=rpaw
   do ix=1,nx
     xx=x(ix)
     dist=abs(raux-xx)
     if(dist<mindist) then
       mindist=dist
       minx=ix
     end if
   end do
   yy=y(minx)
   sig=yy

 end subroutine gaussfit_param2_findsign

end subroutine gaussfit_set_param2
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_set_param3
!! NAME
!!  gaussfit_set_param3
!!
!! FUNCTION
!!  Sets parameters for LSF
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_set_param3(nterm,nparam,param,sep)

!Arguments -------------------------------
 integer,intent(in)::nterm,nparam
 real(dp),intent(in)::sep
 !real(dp),intent(in)::x(nx),y(nx)
 real(dp),intent(out)::param(nparam)

!Local variables-------------------------------
 integer::ii,jj

! *********************************************************************

 param(:)=1.0d0
!
 do jj=1,nterm
   ii=jj+nterm
   param(ii)=sep**(jj)
!  param(ii)=raux*real(i,dp)
 end do
!
 do jj=1,nterm
   ii=jj+3*nterm
   param(ii)=sep**(jj)
!  param(ii)=raux*real(i,dp)
 end do

end subroutine gaussfit_set_param3
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_set_param4
!! NAME
!!  gaussfit_set_param4
!!
!! FUNCTION
!!  Sets parameters for LSF
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_set_param4(nparam,param)

!Arguments -------------------------------
 integer,intent(in)::nparam
 real(dp),intent(out)::param(nparam)

!Local variables-------------------------------

! *********************************************************************

 param(:)=1.0d0

end subroutine gaussfit_set_param4
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_set_param5
!! NAME
!!  gaussfit_set_param5
!!
!! FUNCTION
!!  Sets parameters for LSF
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_set_param5(nterm,nparam,nx,param,rpaw,y)

!Arguments -------------------------------
 integer,intent(in)::nterm,nparam,nx
 real(dp),intent(in)::rpaw
 real(dp),intent(in)::y(nx)
 real(dp),intent(out)::param(nparam)

!Local variables-------------------------------
 integer::ix
 real(dp)::raux,a1,r_c,m

! *********************************************************************

 param(:)=1.0d0
!
!alpha1
!a1=maxval( y(:),nx )
!maxval gives problems in some architectures:
 a1=-9999999
 do ix=1,nx
   if(a1<y(ix)) a1=y(ix)
 end do
 param(1:nterm)=a1
!
!alpha2
!
 r_c=rpaw+0.5d0 !paw sphere + a bit more.
!this is not arbitrary since it is an initial guess
 m=0.01d0
 raux=log(a1/m)/r_c**2
 param(nterm+1:nterm*2)=raux

end subroutine gaussfit_set_param5
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_constrains_init
!! NAME
!!  gaussfit_constrains_init
!!
!! FUNCTION
!!  Initialise constrains for LSF.
!!  It will constrain the Gaussians width
!!  It will also constraint the Delta use in the LSF algorithm, to jump slowly at each step.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_constrains_init(cons1,cons2,limit,nparam,nterm,nx,option,rpaw,y)

!Arguments -------------------------------
 integer,intent(in)::nterm,option,nparam,nx
 integer,intent(out)::cons2(nparam)
 real(dp),intent(in)::rpaw
 real(dp),intent(in)::y(nx)
 real(dp),intent(out)::cons1(nparam),limit(nparam)

!Local variables-------------------------------
 integer :: ix
 real(dp)::rc,a1,mm,raux

! *********************************************************************

!
!DEFAULT: no weight
 cons1(:)=1.d0
 limit(:)=0.d0
 cons2=1
!
 if(option==4) then
     cons1(1:nterm)=0.2d0
     cons1(nterm+1:nterm*2)=0.3d0
 end if
!
 if(option==4) then
!  parameters

!  a1=maxval( y(:),nx)/real(nterm,dp)
!  maxval gives problems in some architectures:
   a1=-9999999
   do ix=1,nx
     if(a1<y(ix)) a1=y(ix)
   end do
   a1=a1/real(nterm,dp)

   mm=0.01
!
   rc=1.7d0*rpaw
   raux=log(a1/mm)/rc**2
!
!  Constraint exponential of gaussians to be positive (multiplied by -1),
!  so that it decays to zero.
!  Constraint as well its value, so that it does not get too big
!  and it decays soon,
!  This prevents that a gaussian grows at a very large x value.
   cons2(nterm+1:nterm*2)=restricted_and_positive
   limit(nterm+1:nterm*2)=raux
 end if

end subroutine gaussfit_constrains_init
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_apply_constrains
!! NAME
!!  gaussfit_apply_constrains
!!
!! FUNCTION
!!  Apply constrains to get new set of parameters
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_paw_gaussfit
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_apply_constrains(const,limit,nparam,ioparams)

!Arguments -------------------------------
 integer,intent(in):: nparam
 integer,intent(in):: const(nparam)
 real(dp),intent(in):: limit(nparam)
 real(dp),intent(inout):: ioparams(nparam)

!Local variables-------------------------------
 integer::ii

! *********************************************************************

 do ii=1,nparam
   if(const(ii)==restricted .or. const(ii)==restricted_and_positive) then
     if(ioparams(ii)<limit(ii)) ioparams(ii)=limit(ii)
   end if
   if(const(ii)==positive .or. const(ii)==restricted_and_positive ) ioparams(ii)=abs(ioparams(ii))

 end do

end subroutine gaussfit_apply_constrains
!!***

!----------------------------------------------------------------------

!!****f* m_paw_gaussfit/gaussfit_projector
!! NAME
!!  gaussfit_projector
!!
!! FUNCTION
!! Fit tproj to Gaussians
!!
!! INPUTS
!!  basis_size= size of the PAW basis
!!  orbitals= indicates the l quantum number for all orbitals.
!!  rpaw= PAW radius
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  tproj= projectors
!!  maxterm= maximum number of terms used to fit the projectors.
!!  mparam= maximum number of parameters (Gaussian coefficients and factors) used.
!!
!! OUTPUT
!!  nparam_array= number of parameters found.
!!  param = parameters found  (Gaussian coefficients and factors).
!!
!! NOTES
!! chisq=accuracy_p= sum_x abs(f(x)-y(x))/nx.
!!    nx is the number of points, f(x) and y(x) are the fitted and original functions.
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      gaussfit_main,paw_spline,paw_splint,pawrad_deducer0,pawrad_free
!!      pawrad_init,wrtout
!!
!! SOURCE

subroutine gaussfit_projector(basis_size,mparam,nparam_array,nterm_bounds,orbitals,param,pawrad,&
& rpaw,tproj,comm_mpi)

!Arguments ------------------------------------
 integer,intent(in) :: basis_size
 integer,intent(in) :: orbitals(basis_size)
 integer, optional,intent(in) :: comm_mpi
 real(dp),intent(in) :: rpaw
 type(pawrad_type),intent(in) :: pawrad
 real(dp),intent(in) :: tproj(:,:)
 integer,intent(in) :: mparam,nterm_bounds(2)
 integer,intent(out) :: nparam_array(basis_size)
 real(dp),intent(out) :: param(mparam,basis_size)
 type(pawrad_type)::mesh_tmp

!Local variables ------------------------------
 integer :: ibasis,ierr,il,ir
 integer :: msz1,msz2,option
 real(dp) :: raux(1),rr(1)
 real(dp),allocatable :: d2(:),tproj_tmp1(:),tproj_tmp2(:)
 character(len=500) :: msg
 character(80) :: outfile
 !debug: uncomment
 !integer::i,nterm ,unitp
 !real(dp),allocatable::y(:)
 !end debug

!************************************************************************

 if(size(tproj,2)<basis_size) then
   msg = 'wrong size for tproj in gaussfit_projector!'
   MSG_BUG(msg)
 end if

 option=4  !see gaussfit_main
 nparam_array(:)=0
 msz1=min(pawrad_ifromr(pawrad,rpaw)+2,size(tproj,1))
 !msz1=pawrad%mesh_size

!Augment the mesh size
!this is to make the Gaussians go to zero after paw_radius
!This is done by creating a new pawrad objet: mesh_tmp
!Change to a linear grid:
 mesh_tmp%mesh_type=1 !linear grid
 mesh_tmp%rstep=0.0005 !very fine grid
 msz2=ceiling(pawrad%rmax*two/mesh_tmp%rstep)
 mesh_tmp%lstep=zero !only needed for log grids
 call pawrad_init(mesh_tmp,mesh_size=msz2,mesh_type=mesh_tmp%mesh_type,&
& rstep=mesh_tmp%rstep,lstep=mesh_tmp%lstep)

 LIBPAW_ALLOCATE(tproj_tmp1,(msz1))
 LIBPAW_ALLOCATE(d2,(msz1))
 LIBPAW_ALLOCATE(tproj_tmp2,(msz2))

 do ibasis=1,basis_size

   write(msg,'(a," - Fitting wfn ",i4," to Gaussians")')ch10,ibasis
   call wrtout(std_out,  msg,'COLL')

   tproj_tmp1=zero; d2=zero; tproj_tmp2=zero

!  take out r^il factor:
!  il=psps%indlmn(1,ilmn,itypat)
   il=orbitals(ibasis)

   tproj_tmp1(2:msz1)=tproj(2:msz1,ibasis)/((pawrad%rad(2:msz1)+tol8)**(il))

!  take out 1/r factor from eq.(3) of M. Torrent CMS 42, 337 (2008)
!  since: <phi|proj>=1 from atompaw, and phi=phi*r, alors proj=proj/r

   tproj_tmp1(2:msz1)=tproj_tmp1(2:msz1)/(pawrad%rad(2:msz1))
   call pawrad_deducer0(tproj_tmp1(1:msz1),msz1,pawrad)

!  splint to a different mesh:
!  get second derivative of tproj and store it
   call paw_spline(pawrad%rad,tproj_tmp1(:),msz1,&
&   zero,zero,d2)

   do ir=2,msz2
     rr=mesh_tmp%rad(ir)
     if( rr(1)-rpaw > tol8 ) then
       !after rpaw projectors are zero
       raux=zero
     else
       call paw_splint(msz1,pawrad%rad,&
&       tproj_tmp1(:),d2(:),&
&       1,rr,raux,ierr=ierr)
     end if
     tproj_tmp2(ir)=raux(1)
   end do

!  Obtain the name for the output file
   if(ibasis<10) then
     write(outfile,'("wfn",i1,".fit")')ibasis
   elseif(ibasis<100) then
     write(outfile,'("wfn",i2,".fit")')ibasis
     write(msg,'(a,a,a,a)')ch10,&
&     "ib (basis index) is too big!",ch10,&
&     "Action: check your pseudopotentials"
     MSG_BUG(msg)
   end if

   if(present(comm_mpi)) then
     call gaussfit_main(mparam,nparam_array(ibasis),nterm_bounds,msz2,&
&     param(:,ibasis),mesh_tmp,option,outfile,rpaw,tproj_tmp2,comm_mpi)
   else
     call gaussfit_main(mparam,nparam_array(ibasis),nterm_bounds,msz2,&
&     param(:,ibasis),mesh_tmp,option,outfile,rpaw,tproj_tmp2)
   end if

!  check
!  LIBPAW_ALLOCATE(y,(mesh_tmp%mesh_size))
!  nterm=nparam_array(ibasis)/4
!  call calcgaussc4(nparam_array(ibasis),nterm,mesh_tmp%mesh_size,1,param(:,ibasis),&
!  &  mesh_tmp%rad,y)
!  ! unitp=600+ibasis
!  ! do ir=1,mesh_tmp%mesh_size
!  !  write(unitp,'(2(f16.7,x,f16.7))')mesh_tmp%rad(ir),y(ir)
!  ! end do
!  LIBPAW_DEALLOCATE(y)

 end do

!Deallocate
 call pawrad_free(mesh_tmp)
 LIBPAW_DEALLOCATE(tproj_tmp1)
 LIBPAW_DEALLOCATE(tproj_tmp2)
 LIBPAW_DEALLOCATE(d2)

end subroutine gaussfit_projector
!!***

!----------------------------------------------------------------------

end module m_paw_gaussfit
!!***
