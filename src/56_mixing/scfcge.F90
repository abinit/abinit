!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcge
!!
!! NAME
!! scfcge
!!
!! FUNCTION
!! Compute the next vtrial of the SCF cycle.
!! Uses a conjugate gradient minimization of the total energy
!! Can move only the trial potential (if moved_atm_inside==0), or
!! move the trial atomic positions as well (if moved_atm_inside==1).
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex= if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  dtn_pc(3,natom)=preconditioned change of atomic position, in reduced
!!    coordinates. Will be quickly transferred to f_atm(:,:,i_vrespc(1))
!!  etotal=the actual total energy
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  iscf =5 => SCF cycle, CG based on estimation of energy gradient
!!       =6 => SCF cycle, CG based on true minimization of the energy
!!  isecur=level of security of the computation
!!  istep= number of the step in the SCF cycle
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpicomm=the mpi communicator used for the summation
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  natom=number of atoms
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftot=total number of FFT grid points
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see i_vresid, ivrespc, i_rhor...)
!!  opt_denpot= 0 vtrial (and also f_fftgr) really contains the trial potential
!!              1 vtrial (and also f_fftgr) actually contains the trial density
!!  response= if 0, GS calculation, if 1, RF calculation, intrinsically harmonic !
!!  rhor(cplex*nfft,nspden)=actual density
!!  ucvol=unit cell volume in bohr**3
!!
!! OUTPUT
!! dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential that gave
!!       the input residual of the potential and Hellman-Feynman forces
!!                       at output, it is the new trial potential .
!!  xred(3,natom)=(needed if moved_atm_inside==1)
!!      reduced dimensionless atomic coordinates
!!      at input, those that generated the input residual of the potential
!!      and Hellman-Feynman forces, at output, these are the new ones.
!!  f_fftgr(cplex*nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!   The input vtrial is transferred, at output, in f_fftgr(:,:,1).
!!   The input f_fftgr(:,:,i_vresid(1)) contains the last residual.
!!     the value of i_vresid(1) is transferred to i_vresid(2) at output.
!!   The input f_fftgr(:,:,i_vresid(2)) contains the old residual.
!!     the value of i_vresid(2) is transferred to i_vresid(3) at output.
!!   The input f_fftgr(:,:,i_vresid(3)) contains the previous last residual.
!!   For the preconditioned potential residual, the same logic as for the
!!     the potential residual is used, with i_vrespc replacing i_vresid.
!!   The input rhor is transferred, at output, in f_fft(:,:,i_rhor(2)).
!!   The old density is input in f_fft(:,:,i_rhor(2)), and the value of
!!      i_rhor(2) is transferred to i_rhor(3) before the end of the routine.
!!   The input/output search vector is stored in f_fftgr(:,:,6)
!!  f_atm(3,natom,n_fftgr)=different functions defined for each atom :
!!   The input xred is transferred, at output, in f_atm(:,:,1).
!!   The input f_atm(:,:,i_vresid(1)) contains minus the HF forces.
!!     the value of i_vresid(1) is transferred to i_vresid(2) at output.
!!   The input f_atm(:,:,i_vresid(2)) contains minus the old HF forces.
!!     the value of i_vresid(2) is transferred to i_vresid(3) at output.
!!   The input f_atm(:,:,i_vresid(3)) contains minus the previous old HF forces.
!!   For the preconditioned change of atomic positions, the same logic as for the
!!     the potential residual is used, with i_vrespc replacing i_vresid.
!!   The input/output search vector is stored in f_atm(:,:,6)
!!  i_rhor(2:3)=index of the density (past and previous past) in the array f_fftgr
!!  i_vresid(3)=index of the residual potentials (present, past and previous
!!   past) in the array f_fftgr; also similar index for minus Hellman-Feynman
!!   forces in the array f_atm .
!!  i_vrespc(3)=index of the preconditioned residual potentials
!!                  (present, past and previous past) in the array f_fftgr ;
!!   also similar index for the preconditioned change of atomic position (dtn_pc).
!!
!! TODO
!! This routine is much too difficult to read ! Should be rewritten ...
!! Maybe make separate subroutines for line search and CG step ?!
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      aprxdr,findminscf,sqnormm_v,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
& f_fftgr,initialized,iscf,isecur,istep,&
& i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpicomm,mpi_summarize,&
& natom,nfft,nfftot,nspden,n_fftgr,n_index,opt_denpot,response,rhor,ucvol,vtrial,xred,errid,errmess)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfcge'
 use interfaces_14_hidewrite
 use interfaces_56_mixing, except_this_one => scfcge
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,initialized,iscf,isecur,istep,moved_atm_inside,mpicomm
 integer,intent(in) :: n_fftgr,n_index,natom,nfft,nfftot,nspden,opt_denpot,response
 integer,intent(out) :: dbl_nnsclo, errid
 character(len = 500), intent(out) :: errmess
 logical, intent(in) :: mpi_summarize
 real(dp),intent(in) :: etotal,ucvol
!arrays
 integer,intent(inout) :: i_rhor(n_index),i_vresid(n_index),i_vrespc(n_index)
 real(dp),intent(in) :: dtn_pc(3,natom),rhor(cplex*nfft,nspden)
 real(dp),intent(inout) :: f_atm(3,natom,n_fftgr)
 real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden),xred(3,natom)

!Local variables-------------------------------
!mlinmin gives the maximum number of steps in the line minimization
!   after which the algorithm is restarted (with a decrease of the
!   adaptative trial step length). This number should not be large,
!   since if the potential landscape is harmonic, the number of
!   search steps should be small. If it is large, we are not in the
!   harmonic region, and the CG algorithm will not be really useful,
!   so one can just restart the algorithm ...
!scalars
 integer,parameter :: mlinmin=5
 integer,save :: end_linmin,iline_cge,ilinear,ilinmin,isecur_eff,nlinear
 integer,save :: number_of_restart,status
 integer :: choice,iatom,idir,ifft,iline_cge_input,ilinmin_input,isp
 integer :: testcg,tmp,errid_
 real(dp),save :: d2edv2_old2,d_lambda_old2,dedv_old2,etotal_old
 real(dp),save :: etotal_previous,lambda_adapt,lambda_new,lambda_old,resid_old
 real(dp) :: d2e11,d2e12,d2e22,d2edv2_new,d2edv2_old
 real(dp) :: d2edv2_predict,d_lambda,de1,de2,dedv_mix
 real(dp) :: dedv_new,dedv_old,dedv_predict,determ,etotal_input
 real(dp) :: etotal_predict,gamma,lambda_input,lambda_predict2
 real(dp) :: lambda_predict=1.0_dp,ratio,reduction
 real(dp) :: resid_input,temp
 character(len=500) :: message
!arrays
 real(dp) :: resid_new(1)
 real(dp), allocatable :: tmp_fft1(:,:)

! *************************************************************************

 errid = AB7_NO_ERROR
 dbl_nnsclo = 0

!reduction gives the level of reduction of the error in
!the line minimization to be reached for the minimization to be
!considered successfull
 reduction=0.1_dp

!nlinear increases with the number of times the 2D minimization succeded
!to reach the true minimum directly. It is a measure of the
!degree of parabolicity of the problem, and is used to
!skip some steps by performing extrapolation.
 if(istep==1)then

!  Skipping some steps is sometimes unsecure, so it is possible
!  to make nlinear start at a negative value - if isecur is positive
   isecur_eff=isecur
   nlinear=min(-isecur_eff,0)
   ilinear=0

!  Response function calculation are intrinsically harmonic, so one
!  can shift isecur (by -2), and start with a positive nlinear
   if(response==1)then
     isecur_eff=isecur-2
     nlinear=-isecur_eff
     ilinear=nlinear
   end if

   iline_cge=0
   ilinmin=0
 end if

!Compute actual residual resid_new (residual of f_fftgr(:,:,i_vrespc(1))
 call sqnormm_v(cplex,i_vrespc(1),mpicomm,mpi_summarize,1,nfft,resid_new,n_fftgr,nspden,opt_denpot,f_fftgr)

!Save input residual and ilinmin for final printing
 resid_input=resid_new(1)
 etotal_input=etotal
 ilinmin_input=ilinmin
 iline_cge_input=iline_cge
!Transfer dtn_pc in f_atm
 if(moved_atm_inside==1)then
   f_atm(:,:,i_vrespc(1))=dtn_pc(:,:)
 end if

!=======================================================================
!Now the routine is decomposed in three mutually exclusive parts :
!if(istep==1)then initialize the algorithm
!else if(ilinmin>0)then perform the line minimisation
!else if(ilinmin==0)then determine the new search direction (CG step)
!=======================================================================


!--------------------------------------
!Here initialize the algorithm
 if(istep==1)then

!  At the beginning of each gstate run, lambda_adapt is forced to have the
!  same value, that is 1.0_dp. In the other cases when istep=1 (at different
!  broyden steps, for example), the previously obtained
!  adaptive value is kept.
   if(initialized==0)lambda_adapt=1.0_dp
   lambda_old=0.0_dp
   lambda_input=0.0_dp
   number_of_restart=0
   lambda_new=lambda_adapt

   f_fftgr(:,:,1)=vtrial(:,:)
   f_fftgr(:,:,i_rhor(2))=rhor(:,:)

!  This copy must be written in F77, because of stack problems on the DECs
   do isp=1,nspden
     do ifft=1,cplex*nfft
       f_fftgr(ifft,isp,6)=f_fftgr(ifft,isp,i_vrespc(1))
     end do
   end do
   vtrial(:,:)=f_fftgr(:,:,1)+(lambda_new-lambda_old)*f_fftgr(:,:,6)
   if(moved_atm_inside==1)then
     f_atm(:,:,1)=xred(:,:)
     f_atm(:,:,i_rhor(2))=xred(:,:)
!    There shouldn t be problems with the stack size for this small array.
     f_atm(:,:,6)=f_atm(:,:,i_vrespc(1))
     xred(:,:)=f_atm(:,:,1)+(lambda_new-lambda_old)*f_atm(:,:,6)
   end if
   tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
   tmp=i_vresid(2) ; i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
   ilinmin=1
   resid_old=resid_new(1)
   etotal_old=etotal

   status=0

!  --------------------------------------

!  Here performs the line minimisation
 else if(ilinmin>0)then

   lambda_input=lambda_new

!  The choice with the Brent algorithm has been abandoned in version 1.6.m

!  Compute the approximate energy derivatives dedv_new and dedv_old,
!  from vresid and vresid_old
   choice=2
   call aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&   f_atm,f_fftgr,i_rhor(2),i_vresid,moved_atm_inside,mpicomm,mpi_summarize,&
&   natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)
   d_lambda=lambda_new-lambda_old
   dedv_old=dedv_old/d_lambda
   dedv_new=dedv_new/d_lambda

!  DEBUG
!  write(std_out,'(a,4es12.4,i3)' )' scfcge:lold,lnew,dold,dnew,status',  &
!  &  lambda_old,lambda_new,dedv_old,dedv_new,status
!  ENDDEBUG

   if(status==0 .or. status==3)then
!    
!    Then, compute a predicted point along the line
!    The value of choice determines the minimization algorithm
!    choice=1 uses the two values of the derivative of the energy
!    choice=2 uses the two values of the energy, and and estimate of the
!    second derivative at the mid-point.

     choice=1
     if(iscf==6)choice=2
     call findminscf(choice,dedv_new,dedv_old,dedv_predict,&
&     d2edv2_new,d2edv2_old,d2edv2_predict,&
&     etotal,etotal_old,etotal_predict,&
&     lambda_new,lambda_old,lambda_predict,errid_,message)
     if (errid_ /= AB7_NO_ERROR) then
       call wrtout(std_out,message,'COLL')
     end if

!    Suppress the next line for debugging  (there is another such line)
     status=0

!    DEBUG
!    Keep this debugging feature : it gives access to the investigation of lines
!    in a different approach
!    if(response==1 .and. istep>8)then
!    lambda_predict=1.2d-2
!    if(istep>=15)lambda_predict=lambda_predict-0.002
!    if(istep>=14)stop
!    status=3
!    end if
!    ENDDEBUG

   else
     if(status/=-1)then
       status=-1
       lambda_predict=-2.5_dp
     else
       lambda_predict=lambda_predict+0.1_dp
     end if
   end if

!  If the predicted point is very close to the most recent
!  computed point, while this is the first trial on this line,
!  then we are in the linear regime :
!  nlinear is increased by one unit. For the time being, do this even when
!  moved_atm_inside==1 (the code still works when it is done, but it
!  seems to be a bit unstable). The maximal value of nlinear is 1, except
!  when isecur_eff is a negative number, less than -1.
   if( abs(lambda_predict-lambda_new)/&
&   (abs(lambda_predict)+abs(lambda_new)) < 0.01 .and. ilinmin==1  ) then
!    if(moved_atm_inside==0 .and. nlinear<max(1,-isecur_eff) )nlinear=nlinear+1
     if(nlinear<max(1,-isecur_eff))nlinear=nlinear+1
     ilinear=nlinear
   end if

!  If the predicted point is close to the most recent computed point,
!  or the previous one, set on the flag of end of line minization
   end_linmin=0
   if(abs(lambda_new-lambda_predict)*2.0_dp&
&   /(abs(lambda_predict)+abs(lambda_new)) <reduction) end_linmin=1
   if(abs(lambda_old-lambda_predict)*2.0_dp&
&   /(abs(lambda_predict)+abs(lambda_new)) <reduction) end_linmin=1

   if(status/=0)end_linmin=0

!  Save the closest old lambda, if needed,
!  also examine the reduction of the interval, and eventual stop
!  the present line minimisation, because of convergence (end_linmin=1)
!  Also treat the case in which the predicted value of lambda is negative,
!  or definitely too small in which case the algorithm has to be restarted
!  (not a very good solution, though ...)
!  Finally also treat the case where insufficiently converged
!  density at lambda=0.0_dp happens, which screws up the line minimisation.

!  Here restart the algorithm with the best vtrial.
!  Also make reduction in lambda_adapt
!  DEBUG
!  write(std_out,*)' scfcge : status=',status
!  ENDDEBUG
   if( end_linmin==0 .and. status==0 .and.                               &
&   (  (lambda_predict<0.005_dp*lambda_adapt .and. iscf==5)     .or.  &
&   (abs(lambda_predict)<0.005_dp*lambda_adapt .and. iscf==6).or.  &
&   ilinmin==mlinmin                                      )     )then
     if(number_of_restart>12)then
       errid = AB7_ERROR_MIXING_CONVERGENCE
       write(errmess,'(a,a,i0,a,a,a,a,a)')&
&       'Potential-based CG line minimization not',' converged after ',number_of_restart,' restarts. ',ch10,&
&       'Action : read the eventual warnings about lack of convergence.',ch10,&
&       'Some might be relevant. Otherwise, raise nband. Returning'
       MSG_WARNING(errmess)
       return
     end if
!    Make reduction in lambda_adapt (kind of steepest descent...)
     write(message,'(a,a,a)')&
&     'Potential-based CG line minimization has trouble to converge.',ch10,&
&     'The algorithm is restarted with more secure parameters.'
     MSG_WARNING(message)
     number_of_restart=number_of_restart+1
!    At the second restart, double the number of non-self consistent loops.
     if(number_of_restart>=2)dbl_nnsclo=1
     lambda_adapt=lambda_adapt*0.7_dp
     lambda_new=lambda_adapt
!    If the last energy is better than the old one, transfer the data.
!    Otherwise, no transfer must occur (very simple to code...)
     if(etotal<etotal_old .or. abs(lambda_old)<1.0d-8)then
       f_fftgr(:,:,1)=vtrial(:,:)
       f_fftgr(:,:,i_rhor(2))=rhor(:,:)
       do isp=1,nspden
         do ifft=1,cplex*nfft
           f_fftgr(ifft,isp,6)=f_fftgr(ifft,isp,i_vrespc(1))
         end do
       end do
       if(moved_atm_inside==1)then
         f_atm(:,:,1)=xred(:,:)
         f_atm(:,:,i_rhor(2))=xred(:,:)
         f_atm(:,:,6)=f_atm(:,:,i_vrespc(1))
       end if
       tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
       tmp=i_vresid(2) ; i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
       resid_old=resid_new(1)
       etotal_old=etotal
     end if
     lambda_old=0.0_dp
     ilinmin=1
!    Putting the flag to -1 avoids the usual actions taken with end_linmin=1
     end_linmin=-1
!    Also put ilinear and nlinear to 0
     ilinear=0
     nlinear=0

!    Here lambda_new is the closest to lambda_predict,
!    or lambda_old is still 0.0_dp, while the energy shows that the minimum
!    is away from 0.0_dp (insufficiently converged density at lambda=0.0_dp).
   else if( abs(lambda_new-lambda_predict)<abs(lambda_old-lambda_predict) &
&     .or.                                                           &
&     ( abs(lambda_old)<1.0d-6 .and.                               &
&     ilinmin>1              .and.                               &
&     etotal>etotal_previous         )                           &
&     )then
     f_fftgr(:,:,1)=vtrial(:,:)
     tmp=i_rhor(3) ; i_rhor(3)=i_rhor(2) ; i_rhor(2)=tmp
     f_fftgr(:,:,i_rhor(2))=rhor(:,:)
     tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(2)
     i_vrespc(2)=i_vrespc(1); i_vrespc(1)=tmp;
     tmp=i_vresid(3); i_vresid(3)=i_vresid(2)
     i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
     if(moved_atm_inside==1)then
       f_atm(:,:,1)=xred(:,:)
       f_atm(:,:,i_rhor(2))=xred(:,:)
     end if
     d_lambda_old2=lambda_old-lambda_new
     lambda_old=lambda_new
     etotal_old=etotal
     resid_old=resid_new(1)
     d2edv2_old2=d2edv2_new
     dedv_old=dedv_new
     dedv_old2=dedv_new
!    if(abs(lambda_new-lambda_predict)*2.0_dp&
!    &    /abs(lambda_new+lambda_predict)        <reduction) end_linmin=1

!    Here lambda_old is the closest to lambda_predict (except for avoiding
!    lambda_old==0.0_dp)
   else
     tmp=i_vresid(3) ; i_vresid(3)=i_vresid(1) ; i_vresid(1)=tmp
     f_fftgr(:,:,i_rhor(3))=rhor(:,:)
     if(moved_atm_inside==1) f_atm(:,:,i_rhor(3))=xred(:,:)
     tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(1) ; i_vrespc(1)=tmp
     d_lambda_old2=lambda_new-lambda_old
     etotal_previous=etotal
     d2edv2_old2=d2edv2_old
     dedv_old2=dedv_old
!    if(abs(lambda_old-lambda_predict)*2.0_dp&
!    &    /abs(lambda_old+lambda_predict)        <reduction) end_linmin=1
   end if

!  If the interval has not yet been sufficiently reduced,
!  continue the search
   if(end_linmin==0)then
     lambda_new=lambda_predict

!    DEBUG
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : continue search, lambda_old,lambda_new=',lambda_old,lambda_new
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : f_fftgr(3:4,1,1)=',f_fftgr(3:4,1,1)
!    write(std_out,'(a,2es16.6)' )&
!    &   ' scfcge : f_fftgr(3:4,1,6)=',f_fftgr(3:4,1,6)
!    ENDDEBUG

     vtrial(:,:)=f_fftgr(:,:,1)+(lambda_new-lambda_old)*f_fftgr(:,:,6)
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+(lambda_new-lambda_old)*f_atm(:,:,6)
     end if

     ilinmin=ilinmin+1
!    
!    Here generates a starting point for next line search
   else
     iline_cge=iline_cge+1
     if(end_linmin==1)ilinmin=0
     lambda_old=0.0_dp

!    In order to generate the new step, take into account previous
!    optimal lambdas (including those of previous ion moves),
!    and the selected new one, if it is positive.
!    However, wait iline_cge>1 to select new ones.
!    lambda_adapt has been initialized at 1.0_dp
     if(iline_cge>1 .and. lambda_new>0.0_dp )then
!      Actually compute a geometric mean
       lambda_adapt= ( lambda_adapt**(dble(iline_cge-1)) * abs(lambda_new)) &
&       **(1.0_dp/dble(iline_cge))
!      In order to recover the previous algorithm, it is enough
!      to decomment the next line
!      lambda_adapt=1.0_dp
     end if
     lambda_new=lambda_adapt

     vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
     end if

!    End choice between continue line minim and determine new direction
   end if

!  
!  -------------------------------

!  Here perform the CG step

 else if(ilinmin==0)then

!  Compute the approximate energy derivatives dedv_mix,dedv_new,dedv_old
   choice=3
   call aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&
&   f_atm,f_fftgr,i_rhor(2),i_vresid,moved_atm_inside,mpicomm,mpi_summarize,&
&   natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred)

   dedv_mix=dedv_mix/lambda_new
   dedv_new=dedv_new/lambda_new
   dedv_old=dedv_old/lambda_new

!  DEBUG
!  write(message, '(a,3es12.4)' )' scfcge: lambda_adapt',&
!  &     lambda_adapt
!  call wrtout(std_out,message,'COLL')

!  write(message, '(a,3es12.4)' )' scfcge: dedv_old,dedv_new,dedv_mix',&
!  &     dedv_old,dedv_new,dedv_mix
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG

!  Then, compute a predicted point, either along the line,
!  or in a 2D plane
   testcg=1
   if(testcg==0)then
!    This part corresponds to steepest descent,
!    in which the line minimisation can be done
!    using different algorithms, varying with the value of choice
     choice=1
     if(iscf==6)choice=2
     call findminscf(choice,dedv_new,dedv_old,dedv_predict,&
&     d2edv2_new,d2edv2_old,d2edv2_predict,&
&     etotal,etotal_old,etotal_predict,&
&     lambda_new,lambda_old,lambda_predict,errid_,message)
     if (errid_ /= AB7_NO_ERROR) then
       call wrtout(std_out,message,'COLL')
     end if
     lambda_predict2=0.0_dp
!    Suppress the next line for debugging (there is another such line)
     status=0
   else
!    This part corresponds to conjugate gradient
!    A 2D minimisation is performed
!    oldest direction is labelled 2
!    newest direction is labelled 1
     de1=dedv_old ;  de2=dedv_old2
     d2e11=(dedv_new-dedv_old)/lambda_new
     d2e22=d2edv2_old2
     d2e12=(dedv_mix-dedv_old)/d_lambda_old2
!    The system to be solved is
!    0 = de1 + lambda1 d2e11 + lambda2 d2d12
!    0 = de2 + lambda1 d2e12 + lambda2 d2d22
     determ=d2e11*d2e22-d2e12*d2e12
     lambda_predict=-(de1*d2e22-de2*d2e12)/determ
     lambda_predict2=(de1*d2e12-de2*d2e11)/determ
     d2edv2_new=d2e11 ;  d2edv2_old=d2e11
   end if

!  DEBUG
!  write(message, '(a,5es11.3)' )' scfcge: de1,de2,d2e11,d2e22,d2e12',&
!  &               de1,de2,d2e11,d2e22,d2e12
!  call wrtout(std_out,message,'COLL')
!  write(std_out,'(a,2es12.4)' )' scfcge: la_predict,la_predict2',&
!  &               lambda_predict,lambda_predict2
!  -----
!  write(std_out,*)'residues ',
!  !$       de1+lambda_predict*d2e11+lambda_predict2*d2e12,
!  !$       de2+lambda_predict*d2e12+lambda_predict2*d2e22
!  if(.true.)stop
!  ENDDEBUG
!  

!  Determine the region of the 2D search space
!  in which the predicted point is located,
!  or use linear indicator to decide interpolation
!  and advance to next 2D search.
   end_linmin=0
   write(message, '(a,2i3)' )' nlinear, ilinear',nlinear,ilinear
   call wrtout(std_out,message,'COLL')
   if(lambda_predict<0.0_dp)then
!    Something is going wrong. Just take a reasonable step
!    along the steepest descent direction (Region III).
!    Actually, Region I and region III are treated in the same way later.
!    In effect, this corresponds to restart the algorithm
     end_linmin=3
!    Also put ilinear and nlinear to 0
     ilinear=0
     nlinear=0
!    Decrease the adaptive step to predict next direction
     lambda_adapt=lambda_adapt*0.7_dp
   else if(ilinear>=1) then
!    Region IV : will do an interpolation
     end_linmin=4
     ilinear=ilinear-1
   else if(abs(lambda_predict2)>reduction          .or.&
&     lambda_predict<0.5_dp                .or.&
&     lambda_predict>2.5_dp                .or.&
&     lambda_predict-abs(lambda_predict2)/reduction <0.0_dp  ) then
!    Region II : lambda_predict is not too good, and not too bad.
     end_linmin=2
   else if (abs(1.0_dp-lambda_predict)<reduction)then
!    Region I, the out-of-line point is OK.
     end_linmin=1
   else
!    If everything fails, then region II.
     end_linmin=2
   end if

!  DEBUG
!  write(message, '(a,2es12.4,i2)' )&
!  &     ' scfcge : la_predict, la_predict2, region',&
!  &       lambda_predict,lambda_predict2,end_linmin
!  call wrtout(std_out,message,'COLL')
!  ENDDEBUG

!  Treat region I, in the same way as region III
   if(end_linmin==1 .or. end_linmin==3)then

!    In region I, the line search is
!    along vtrial-vtrial_old.
!    The closest point is the new point
!    thus to be transfered in the "old" locations

     do isp=1,nspden
       do ifft=1,cplex*nfft
         f_fftgr(ifft,isp,6)=(vtrial(ifft,isp)-f_fftgr(ifft,isp,1))/lambda_new
       end do
     end do
     f_fftgr(:,:,1)=vtrial(:,:)
     f_fftgr(:,:,i_rhor(2))=rhor(:,:)
     if(moved_atm_inside==1)then
       f_atm(:,:,6)=(xred(:,:)-f_atm(:,:,1))/lambda_new
       f_atm(:,:,1)=xred(:,:)
       f_atm(:,:,i_rhor(2))=xred(:,:)
     end if
     tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
     tmp=i_vresid(3) ; i_vresid(3)=i_vresid(2)
     i_vresid(2)=i_vresid(1) ; i_vresid(1)=tmp
     d_lambda_old2=-lambda_new
     lambda_old=lambda_new
     etotal_old=etotal
     resid_old=resid_new(1)
     d2edv2_old=d2edv2_new
     dedv_old=dedv_new

!    Region I or III : one is close of the 2D minimum,
!    or lambda_predict was negative (indicate a problem of convergence)
!    Compute next trial potential along the
!    PC residual and not along this search direction.
     ilinmin=0
!    Question : isn t it here that one should prevent region I to called
!    itself more than 1 time ???
!    Here the small difference between region I and region III
     if(end_linmin==3)ilinmin=1
     lambda_old=0.0_dp
     lambda_new=lambda_adapt

     vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
     if(moved_atm_inside==1)then
       xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
     end if
!    The new vtrial has been generated

   else

!    Here region II or IV
     ilinmin=1
     if (lambda_predict==0._dp) then
       gamma=zero
     else
       gamma=lambda_predict2/lambda_predict
     end if
!    Compute new search direction and trial potential
     write(message,*)' compute new search direction '
     call wrtout(std_out,message,'COLL')
     do isp=1,nspden
       do ifft=1,cplex*nfft
         f_fftgr(ifft,isp,6)=(vtrial(ifft,isp)-f_fftgr(ifft,isp,1))/lambda_new+ &
&         gamma*f_fftgr(ifft,isp,6)
       end do
     end do
     vtrial(:,:)=f_fftgr(:,:,1)+ lambda_predict*f_fftgr(:,:,6)
     if(moved_atm_inside==1)then
       f_atm(:,:,6)=(xred(:,:)-f_atm(:,:,1))/lambda_new+ gamma*f_atm(:,:,6)
       xred(:,:)=f_atm(:,:,1)+ lambda_predict*f_atm(:,:,6)
     end if

!    If end_linmin==2, then this vtrial is the good one

     if(end_linmin==2)then

       lambda_old=0.0_dp
       lambda_new=lambda_predict

     else if(end_linmin==4)then

!      predict the result of the computation at the trial potential
!      defined in the end_linmin==2 case
       gamma=lambda_predict2/d_lambda_old2
       ratio=lambda_predict/lambda_new

!      Take care of vtrial
       f_fftgr(:,:,1)=vtrial(:,:)

       ABI_ALLOCATE(tmp_fft1,(cplex*nfft,nspden))
!      Take care of vresid
       tmp_fft1(:,:)=f_fftgr(:,:,i_vresid(2))
       f_fftgr(:,:,i_vresid(2))=tmp_fft1(:,:)&
&       +ratio*(f_fftgr(:,:,i_vresid(1))-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_vresid(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_vresid(3))=tmp_fft1(:,:)

!      Take care of rhor
       tmp_fft1(:,:)=f_fftgr(:,:,i_rhor(2))
       f_fftgr(:,:,i_rhor(2))=tmp_fft1(:,:)&
&       +ratio*(rhor(:,:)-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_rhor(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_rhor(3))=tmp_fft1(:,:)

!      Take care of vrespc
       tmp_fft1(:,:)=f_fftgr(:,:,i_vrespc(2))
       f_fftgr(:,:,i_vrespc(2))=tmp_fft1(:,:)&
&       +ratio*(f_fftgr(:,:,i_vrespc(1))-tmp_fft1(:,:))&
&       +gamma*(f_fftgr(:,:,i_vrespc(3))-tmp_fft1(:,:))
       f_fftgr(:,:,i_vrespc(3))=tmp_fft1(:,:)
       ABI_DEALLOCATE(tmp_fft1)

       if(moved_atm_inside==1)then
         do idir=1,3
           do iatom=1,natom

!            Take care of xred
             f_atm(idir,iatom,1)=xred(idir,iatom)

!            Take care of -HF forces
             temp=f_atm(idir,iatom,i_vresid(2))
             f_atm(idir,iatom,i_vresid(2))=f_atm(idir,iatom,i_vresid(2))&
&             +ratio*(f_atm(idir,iatom,i_vresid(1))-f_atm(idir,iatom,i_vresid(2)))&
&             +gamma*(f_atm(idir,iatom,i_vresid(3))-f_atm(idir,iatom,i_vresid(2)))
             f_atm(idir,iatom,i_vresid(3))=temp

!            Take care of old xreds
             temp=f_atm(idir,iatom,i_rhor(2))
             f_atm(idir,iatom,i_rhor(2))=f_atm(idir,iatom,i_rhor(2))&
&             +ratio*(   xred(idir,iatom)          -f_atm(idir,iatom,i_rhor(2)))&
&             +gamma*(f_atm(idir,iatom,i_rhor(3))-f_atm(idir,iatom,i_rhor(2)))
             f_atm(idir,iatom,i_rhor(3))=temp

!            Take care of preconditioned changes of atomic positions
             temp=f_atm(idir,iatom,i_vrespc(2))
             f_atm(idir,iatom,i_vrespc(2))=f_atm(idir,iatom,i_vrespc(2))&
&             +ratio*(f_atm(idir,iatom,i_vrespc(1))-f_atm(idir,iatom,i_vrespc(2)))&
&             +gamma*(f_atm(idir,iatom,i_vrespc(3))-f_atm(idir,iatom,i_vrespc(2)))
             f_atm(idir,iatom,i_vrespc(3))=temp

           end do
         end do
       end if

!      Since we are at the 2D minimum, the derivative is supposed
!      to vanish. Note that dedv_old should not change, by contrast.
       dedv_old2=0.0_dp
       d_lambda_old2=-lambda_predict
       d2edv2_old2=-dedv_old/lambda_predict
       lambda_old=lambda_predict
       ilinmin=0

!      So, jump to the next line
       iline_cge=iline_cge+1
       write(message,*)' energy CG update : after 2D interpolation,'
       call wrtout(std_out,message,'COLL')
       write(message,*)'    computation in the next plane '
       call wrtout(std_out,message,'COLL')
       write(message,*)
       call wrtout(std_out,message,'COLL')
       lambda_old=0.0_dp
       lambda_new=lambda_adapt

       vtrial(:,:)=f_fftgr(:,:,1)+lambda_new*f_fftgr(:,:,i_vrespc(2))
       if(moved_atm_inside==1)then
         xred(:,:)=f_atm(:,:,1)+lambda_new*f_atm(:,:,i_vrespc(2))
       end if

!      The new trial potential is now generated

!      End the specific treatment of region IV
     end if
!    
!    End the choice between treatment of region I, II, or IV
   end if

!  End of choice between initialisation or more developed parts of the CG algorithm
 else
   errid = AB7_ERROR_MIXING_ARG
   errmess = 'scfcge : BUG You should not be here ! '
   return
 end if

!--------------------------------------

!Write information : it will be easy to read by typing  grep scfcge logfile

 if(istep==1)then
   write(message,'(a,a,a)') ' scfcge:',ch10,' scfcge:istep-iline_cge-ilinmin lambda      etot             resid '
   call wrtout(std_out,message,'COLL')
 end if

 if(ilinmin_input/=0 .or. istep==1)then
!  Usual line minimisation step

   if(iline_cge_input<10)then
     write(message, '(a,i4,a,i1,a,i1,es13.4,es20.12,es12.4)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-',ilinmin_input,lambda_input,etotal_input,resid_input
   else
     write(message, '(a,i3,a,i2,a,i1,es13.4,es20.12,es12.4)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-',ilinmin_input,lambda_input,etotal_input,resid_input
   end if
   call wrtout(std_out,message,'COLL')

   if( (end_linmin==1.or.end_linmin==-1) .and. istep/=1 )then

     if(end_linmin==1)then
       write(message, '(a,es13.4,a,i2,a,a)' )&
&       ' scfcge: predict         ',lambda_predict,&
&       ' suff. close => next line, ilinear=',ilinear,ch10,&
&       ' scfcge:'
     else if(end_linmin==-1)then
       write(message, '(a,es13.4,a,a,a)' )&
&       ' scfcge: predict         ',lambda_predict,&
&       ' restart the algorithm ',ch10,&
&       ' scfcge:'
     end if
     call wrtout(std_out,message,'COLL')

     if(iline_cge_input<9)then
       write(message, '(a,i4,a,i1,a,i1,es13.4,es20.12,es12.4)' ) &
&       ' scfcge: start   ',istep,'-',iline_cge,'-',0,0.0,etotal_old,resid_old
     else
       write(message, '(a,i3,a,i2,a,i1,es13.4,es20.12,es12.4)' ) &
&       ' scfcge: start   ',istep,'-',iline_cge,'-',0,0.0,etotal_old,resid_old
     end if
     call wrtout(std_out,message,'COLL')

   else if(istep/=1) then
     write(message, '(a,es13.4,a)' )&
&     ' scfcge: predict         ',lambda_predict,&
&     ' not close enough => continue minim.'
     call wrtout(std_out,message,'COLL')
   end if

 else
!  CG prediction
   if(iline_cge_input<10)then
     write(message, '(a,i4,a,i1,a,es11.4,es20.12,es12.4,a,i1)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-off',&
&     lambda_adapt,etotal_input,resid_input,', end=',end_linmin
   else
     write(message, '(a,i3,a,i2,a,es11.4,es20.12,es12.4,a,i1)' )&
&     ' scfcge: actual  ',istep,'-',iline_cge_input,'-off',&
&     lambda_adapt,etotal_input,resid_input,', end=',end_linmin
   end if
   call wrtout(std_out,message,'COLL')

   if(end_linmin==4)then
     write(message, '(a)' ) ' scfcge:'
     call wrtout(std_out,message,'COLL')
   end if

 end if

end subroutine scfcge
!!***
