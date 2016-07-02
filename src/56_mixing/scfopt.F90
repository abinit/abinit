!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfopt
!!
!! NAME
!! scfopt
!!
!! FUNCTION
!! Compute the next vtrial of the SCF cycle.
!! Possible algorithms are : simple mixing, Anderson (order 1 or 2), Pulay
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex= if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  iscf= 2 => simple mixing
!!      = 3,4 => Anderson mixing
!!      = 7 => Pulay mixing
!!  istep= number of the step in the SCF cycle
!!  mpicomm=the mpi communicator used for the summation
!!  comm_atom=the mpi communicator over atoms ; PAW only (optional argument)
!!  mpi_summarize=set it to .true. if parallelisation is done over FFT
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see ivrespc, i_vtrial...)
!!  opt_denpot= 0 vtrial (and also f_fftgr) really contains the trial potential
!!              1 vtrial (and also f_fftgr) actually contains the trial density
!!  pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential that gave
!!     the input preconditioned residual potential
!!     at output, it is the new trial potential .
!!  f_fftgr(cplex*nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!   The input vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(1)).
!!   The old vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(2)).
!!   The input preconditioned residual potential is in f_fftgr(:,:,i_vrespc(1))
!!   Two input old preconditioned residual potentials in f_fftgr(:,:,i_vrespc(2)) and f_fftgr(:,:,i_vrespc(3))
!!    Before output a permutation of i_vrespc(1), i_vrespc(2) and i_vrespc(3) occurs, without
!!    actually copying all the data (change of pointer).
!!  i_vrespc(n_index)=index of the preconditioned residual potentials (present and past) in the array f_fftgr
!!  i_vtrial(n_index)  =indices of the potential (present and past) in the array f_fftgr
!!  ==== if usepaw==1
!!    f_paw(npawmix,n_fftgr*mffmem*usepaw)=different functions used for PAW
!!                                           (same as f_fftgr but for spherical part)
!!    vpaw(npawmix*usepaw)=at input, the aug. occupancies (rhoij) that gave
!!                               the input preconditioned residual potential
!!                           at output, it is the new aug. occupancies.
!!
!! PARENTS
!!      m_ab7_mixing
!!
!! CHILDREN
!!      dgetrf,dgetri,dotprodm_v,sqnormm_v,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine scfopt(cplex,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,&
& mpicomm,mpi_summarize,nfft,npawmix,nspden,n_fftgr,&
& n_index,opt_denpot,pawoptmix,usepaw,vpaw,vresid,vtrial,errid,errmess, &
& comm_atom) ! optional

 use defs_basis
 use m_linalg_interfaces
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'scfopt'
 use interfaces_14_hidewrite
 use interfaces_56_mixing, except_this_one => scfopt
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,iscf,istep,n_fftgr,n_index,nfft
 integer,intent(in) :: npawmix,nspden,opt_denpot,pawoptmix,usepaw,mpicomm
 integer, intent(in),optional :: comm_atom
 integer,intent(out) :: errid
 character(len = 500), intent(out) :: errmess
 logical, intent(in) :: mpi_summarize
 real(dp), intent(out) :: vresid
!arrays
 integer,intent(inout) :: i_vrespc(n_index),i_vtrial(n_index)
 real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*usepaw),vpaw(npawmix*usepaw)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
!Local variables-------------------------------
!scalars
 integer,parameter :: npulaymax=50
 integer :: i_vstore,ierr,ifft,ii,index,isp,jj,comm_atom_,niter,npulay,tmp
 real(dp),save :: prod_resid_old,resid_old,resid_old2
 real(dp) :: aa1,aa2,bb,cc1,cc2,current,det,lambda,lambda2,resid_best
 character(len=500) :: message
!arrays
 integer,allocatable :: ipiv(:)
 real(dp),save :: amat(npulaymax+1,npulaymax+1)
 real(dp) :: mpibuff(2),prod_resid(1),prod_resid2(1),resid_new(1)
 real(dp),allocatable :: alpha(:),amatinv(:,:),amat_paw(:),rwork(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' scfopt : enter ; istep,iscf ',istep,iscf
!ENDDEBUG

 errid = AB7_NO_ERROR

 comm_atom_=xmpi_comm_self; if(present(comm_atom)) comm_atom_=comm_atom

 i_vstore=i_vtrial(1)
 if (iscf==4) i_vstore=i_vtrial(2)
 if (iscf==7) then
   if (modulo(n_fftgr, 2) == 0 ) then
     npulay=(n_fftgr-2)/2
   else
     npulay=(n_fftgr-1)/2
   end if
   i_vstore=i_vtrial(npulay)
 else
   npulay=0
 end if

!Compute the new residual resid_new, from f_fftgr/f_paw(:,:,i_vrespc(1))
 call sqnormm_v(cplex,i_vrespc(1),mpicomm,mpi_summarize,1,nfft,resid_new,n_fftgr,nspden,opt_denpot,f_fftgr)
 if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
     resid_new(1)=resid_new(1)+f_paw(index,i_vrespc(1))**2
   end do
   call xmpi_sum(resid_new(1),comm_atom_,ierr)
 end if
 vresid = resid_new(1)

!_______________________________________________________________
!Here use only the preconditioning, or initialize the other algorithms

 if(istep==1 .or. iscf==2)then

   write(message,'(2a)') ch10,'Simple mixing update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential :',resid_new(1)
   call wrtout(std_out,message,'COLL')

!  Store information for later use
   if (iscf==3.or.iscf==4) resid_old=resid_new(1)
   if (iscf==7) then
     amat(:,:)=zero
     amat(1,1)=resid_new(1)
   end if

!  Compute new vtrial (and new rhoij if PAW)
   if (iscf/=2) f_fftgr(:,:,i_vstore)=vtrial(:,:)
   vtrial(:,:)=vtrial(:,:)+f_fftgr(:,:,i_vrespc(1))
   if (usepaw==1) then
     if (iscf/=2) f_paw(:,i_vstore)=vpaw(:)
     vpaw(:)=vpaw(:)+f_paw(:,i_vrespc(1))
   end if

!  _______________________________________________________________
!  Here Anderson algorithm using one previous iteration
 else if((istep==2 .or. iscf==3).and.iscf/=7)then

   write(message,'(2a)') ch10,'Anderson update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential: ',resid_new(1)
   call wrtout(std_out,message,'COLL')

!  Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
   call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
     end do
     call xmpi_sum(prod_resid(1),comm_atom_,ierr)
   end if

!  Compute mixing factor
   lambda=(resid_new(1)-prod_resid(1))/(resid_new(1)+resid_old-2*prod_resid(1))
   write(message,*)' mixing of old trial potential    :',lambda
   call wrtout(std_out,message,'COLL')

!  Evaluate best residual square on the line
   resid_best=(1.0_dp-lambda)*(1.0_dp-lambda)*resid_new(1)&
&   +(1.0_dp-lambda)*lambda        *2*prod_resid(1)&
&   +lambda        *lambda        *resid_old
   write(message,*)' predicted best residual square on the line: ',resid_best
   call wrtout(std_out,message,'COLL')

!  Store information for later use
   if (iscf==4) then
     prod_resid_old=prod_resid(1)
     resid_old2=resid_old
   end if
   resid_old=resid_new(1)

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=(one-lambda)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&       +lambda      *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=(one-lambda)*(current                 +f_paw(index,i_vrespc(1)))&
&     +lambda      *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))
     f_paw(index,i_vstore)=current
   end do

!  _______________________________________________________________
!  Here Anderson algorithm using two previous iterations
 else if(iscf==4.and.iscf/=7)then

   write(message,'(2a)') ch10,'Anderson (order 2) update:'
   call wrtout(std_out,message,'COLL')

   write(message,*)' residual square of the potential :',resid_new(1)
   call wrtout(std_out,message,'COLL')

!  Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
   call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
     end do
   end if

!  Compute prod_resid2 from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(3))
   call dotprodm_v(cplex,1,prod_resid2,i_vrespc(1),i_vrespc(3),mpicomm,mpi_summarize,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
     do index=1,npawmix
       prod_resid2(1)=prod_resid2(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(3))
     end do
!    MPI reduction
     mpibuff(1)=prod_resid(1);mpibuff(2)=prod_resid2(1)
     call xmpi_sum(mpibuff,comm_atom_,ierr)
     prod_resid(1)=mpibuff(1);prod_resid2(1)=mpibuff(2)
   end if

!  Compute mixing factors
   aa1=resid_new(1)+resid_old -two*prod_resid (1)
   aa2=resid_new(1)+resid_old2-two*prod_resid2(1)
   bb =resid_new(1)+prod_resid_old-prod_resid(1)-prod_resid2(1)
   cc1=resid_new(1)-prod_resid (1)
   cc2=resid_new(1)-prod_resid2(1)
   det=aa1*aa2-bb*bb
   lambda =(aa2*cc1-bb*cc2)/det
   lambda2=(aa1*cc2-bb*cc1)/det
   write(message,*)' mixing of old trial potentials   :',lambda,lambda2
   call wrtout(std_out,message,'COLL')

!  Store information for later use
   prod_resid_old=prod_resid(1)
   resid_old2=resid_old
   resid_old=resid_new(1)

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=&
&       (one-lambda-lambda2)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&       +lambda             *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))&
&       +lambda2            *(f_fftgr(ifft,isp,i_vtrial(2))+f_fftgr(ifft,isp,i_vrespc(3)))
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=&
&     (one-lambda-lambda2)*(current                 +f_paw(index,i_vrespc(1)))&
&     +lambda             *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))&
&     +lambda2            *(f_paw(index,i_vtrial(2))+f_paw(index,i_vrespc(3)))
     f_paw(index,i_vstore)=current
   end do

!  _______________________________________________________________
!  Here Pulay algorithm
 else if(iscf==7)then

   niter=min(istep,npulay+1)

   write(message,'(2a,i2,a)') ch10,' Pulay update with ',niter-1,' previous iterations:'
   call wrtout(std_out,message,'COLL')

   if (npulay>npulaymax) then
     errid = AB7_ERROR_MIXING_CONVERGENCE
     write(errmess, '(4a)' ) ch10,&
&     ' scfopt : ERROR - ',ch10,&
&     '  Too much iterations required for Pulay algorithm (<50) !'
     return
   end if

!  Compute "A" matrix
   if (istep>npulay+1) then
     do jj=1,niter-1
       do ii=1,niter-1
         amat(ii,jj)=amat(ii+1,jj+1)
       end do
     end do
   end if
   if (usepaw==1.and.pawoptmix==1) then
     ABI_ALLOCATE(amat_paw,(niter))
     amat_paw(:)=zero
     do ii=1,niter
       do index=1,npawmix
         amat_paw(ii)=amat_paw(ii)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(1+niter-ii))
       end do
     end do
     call xmpi_sum(amat_paw,comm_atom_,ierr)
   end if
   do ii=1,niter
     call dotprodm_v(cplex,1,amat(ii,niter),i_vrespc(1),i_vrespc(1+niter-ii),mpicomm,mpi_summarize,1,1,&
&     nfft,n_fftgr,n_fftgr,nspden,opt_denpot,f_fftgr,f_fftgr)
     if (usepaw==1.and.pawoptmix==1) amat(ii,niter)=amat(ii,niter)+amat_paw(ii)
     if (ii<niter) amat(niter,ii)=amat(ii,niter)
   end do
   if (usepaw==1.and.pawoptmix==1)then
     ABI_DEALLOCATE(amat_paw)
   end if

!  Invert "A" matrix
   ABI_ALLOCATE(amatinv,(niter,niter))
   amatinv(1:niter,1:niter)=amat(1:niter,1:niter)
   ABI_ALLOCATE(ipiv,(niter))
   ABI_ALLOCATE(rwork,(niter))
   call dgetrf(niter,niter,amatinv,niter,ipiv,ierr)
   call dgetri(niter,amatinv,niter,ipiv,rwork,niter,ierr)
   ABI_DEALLOCATE(ipiv)
   ABI_DEALLOCATE(rwork)

!  Compute "alpha" factors
   ABI_ALLOCATE(alpha,(niter))
   alpha=zero
   det=zero
   do ii=1,niter
     do jj=1,niter
       alpha(ii)=alpha(ii)+amatinv(jj,ii)
       det=det+amatinv(jj,ii)
     end do
   end do
   alpha(:)=alpha(:)/det
   ABI_DEALLOCATE(amatinv)
   write(message,'(a,5(1x,g10.3))')' mixing of old trial potential : alpha(m:m-4)=',(alpha(ii),ii=niter,max(1,niter-4),-1)
   call wrtout(std_out,message,'COLL')

!  Save latest trial potential and compute new trial potential
   do isp=1,nspden
     do ifft=1,cplex*nfft
       current=vtrial(ifft,isp)
       vtrial(ifft,isp)=alpha(niter)*(current+f_fftgr(ifft,isp,i_vrespc(1)))
       do ii=niter-1,1,-1
         vtrial(ifft,isp)=vtrial(ifft,isp)+alpha(ii) &
&         *(f_fftgr(ifft,isp,i_vtrial(niter-ii))+f_fftgr(ifft,isp,i_vrespc(1+niter-ii)))
       end do
       f_fftgr(ifft,isp,i_vstore)=current
     end do
   end do

!  PAW: save latest rhoij and compute new rhoij
   do index=1,npawmix
     current=vpaw(index)
     vpaw(index)=alpha(niter)*(current+f_paw(index,i_vrespc(1)))
     do ii=niter-1,1,-1
       vpaw(index)=vpaw(index)+alpha(ii) &
&       *(f_paw(index,i_vtrial(niter-ii))+f_paw(index,i_vrespc(1+niter-ii)))
     end do
     f_paw(index,i_vstore)=current
   end do

   ABI_DEALLOCATE(alpha)
!  _______________________________________________________________
!  End of choice of optimization method
 end if

!Permute potential indices
 if (iscf==3) then
   tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
 else if (iscf==4) then
   tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
   tmp=i_vtrial(2) ; i_vtrial(2)=i_vtrial(1) ; i_vtrial(1)=tmp
 else if (iscf==7) then
   tmp=i_vtrial(  npulay)
   do ii=  npulay,2,-1
     i_vtrial(ii)=i_vtrial(ii-1)
   end do
   i_vtrial(1)=tmp
   tmp=i_vrespc(1+npulay)
   do ii=1+npulay,2,-1
     i_vrespc(ii)=i_vrespc(ii-1)
   end do
   i_vrespc(1)=tmp
 end if

end subroutine scfopt
!!***
