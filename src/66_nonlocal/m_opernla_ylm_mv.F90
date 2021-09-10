!!****m* ABINIT/m_opernla_ylm_mv
!! NAME
!!  m_opernla_ylm_mv
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2021 ABINIT group (LB,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_opernla_ylm_mv

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
#if defined HAVE_OPENMP
 use OMP_LIB
#endif

 use defs_abitypes, only : MPI_type
 use m_time,        only : timab

 implicit none

 private
!!***

 public :: opernla_ylm_mv
 integer,public,save :: opernla_mv_counter = -1
 integer,public,save :: opernla_mv_dgemv_counter = -1
!!***
!!***

contains
!!***

!!****f* ABINIT/opernla_ylm_mv
!! NAME
!! opernla_ylm_mv
!!
!! FUNCTION
!! "matrix-vector" alternative implementation of "opernla_ylm".
!!
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
!!
!! Here this is done in 3 steps:
!! (1) compute for every g   : scal(g)   = c(g).exp(2pi.i.g.R)
!! (2) compute for every lmn : scal(lmn) = Sum_g[scal(g).f_nl(g).Y_lm(g)]
!! (3) compute for every lmn : <p_lmn|c> = 4pi/sqrt(vol).(i)^l.scal(lmn)
!!
!! Step (2) is a real-matrix/complex-vector multiplication, here two options are possible:
!! - case nloalg(1)=2 : compute the real and imaginary parts separately using two calls of DGMEV
!! - case nloalg(1)=3 : in order to read the matrix only once, we compute both real and imaginary parts at the same time "by hand"
!!
!! Depending on the achitecture and the available blas library, one option could be more interesting than an other...
!!
!! INPUTS
!!  choice=chooses possible output:
!!         if choice>=0: compute projected scalars
!!         if choice<0: same as choice>0 but use already computed projected scalars
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  ia3=gives the number of the first atom in the subset presently treated
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  matblk=dimension of the array ph3d
!!  mpi_enreg=information about MPI parallelization
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  ucvol=unit cell volume (bohr^3)
!!  vect(2,npw*my_nspinor)=starting vector in reciprocal space
!!
!! OUTPUT
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars
!!
!! SIDE EFFECTS
!!
!! NOTES
!! 1-Not available yet for openMP
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!! 3-projector derivatives (abs(choice)>1) are not implemented yet
!!
!! PARENTS
!!      m_cgprj,m_nonlop_ylm
!!
!! CHILDREN
!!      dgemv,timab,xmpi_sum
!!
!! SOURCE

subroutine opernla_ylm_mv(choice,cplex,dimffnl,ffnl,gx,&
&       ia3,indlmn,istwf_k,matblk,mpi_enreg,nincat,nlmn,&
&       nloalg,npw,nspinor,ph3d,ucvol,vect)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,dimffnl,ia3,istwf_k,matblk
 integer,intent(in) :: nincat,nlmn,npw,nspinor
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,nlmn),nloalg(3)
 real(dp),intent(in),target :: ffnl(npw,dimffnl,nlmn)
 real(dp),intent(in) :: ph3d(2,npw,matblk)
 real(dp),intent(in) :: vect(:,:)
 real(dp),intent(out) :: gx(cplex,nlmn,nincat,nspinor)

!Local variables-------------------------------
!scalars
 logical :: use_dgemv
 integer :: ia,iaph3d,ierr,il,ilmn,ipw,ipw0,ipwshft,ispinor,jpw
 real(dp) :: wt
!arrays
 real(dp) :: buffer_r,buffer_i,tsec(2)
 real(dp),pointer :: ffnl_loc(:,:)
 real(dp),allocatable :: scali(:),scalr(:)
 real(dp),allocatable :: scalr_lmn(:),scali_lmn(:)
 complex(dp) :: ctmp,cil(4)
! *************************************************************************

 if (choice==-1) return

!Some checks
 if (abs(choice)>1) then
   ABI_ERROR('Only abs(choice)<=1 is available for now.')
 end if
 if (nloalg(1)<2.or.nloalg(1)>10) then
   ABI_ERROR('nloalg(1) should be between 2 and 10.')
 end if
! nthreads=1
!#if defined HAVE_OPENMP
! nthreads=OMP_GET_NUM_THREADS()
!#endif
! if (nthreads>1) then
!   ABI_ERROR('Only nthreads=1 is available for now.')
! end if

 use_dgemv = nloalg(1)==2.or.nloalg(1)==5.or.nloalg(1)==7
 if (choice>=0.or.abs(choice)>1) then
   if (use_dgemv) then
     if(opernla_mv_dgemv_counter>=0) opernla_mv_dgemv_counter = opernla_mv_dgemv_counter + 1
   else
     if(opernla_mv_counter>=0) opernla_mv_counter = opernla_mv_counter + 1
   end if
 end if

!Useful variables
 wt=four_pi/sqrt(ucvol);if (cplex==1) wt=2.d0*wt
 ipw0=1;if (istwf_k==2.and.mpi_enreg%me_g0==1) ipw0=2

!Allocate work space
 ABI_MALLOC(scalr,(npw))
 ABI_MALLOC(scali,(npw))
 ABI_MALLOC(scalr_lmn,(nlmn))
 ABI_MALLOC(scali_lmn,(nlmn))

 ffnl_loc => ffnl(:,1,:)

! i^l
 cil(1) = ( 1.0_DP, 0.0_DP) * wt
 cil(2) = ( 0.0_DP, 1.0_DP) * wt
 cil(3) = (-1.0_DP, 0.0_DP) * wt
 cil(4) = ( 0.0_DP,-1.0_DP) * wt

!$OMP PARALLEL PRIVATE(il,ilmn,ipw,jpw), &
!$OMP PRIVATE(ispinor,ipwshft,ia,iaph3d)

!Loop on spinorial components
 do ispinor =1,nspinor
   ipwshft=(ispinor-1)*npw

!  Loop on atoms (blocking)
   do ia=1,nincat
     iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1
!    Step (1) : Compute scal(g) = c(g).exp(2pi.i.g.R)
!$OMP DO
     do ipw=ipw0,npw
       jpw=ipw+ipwshft
       scalr(ipw)=(vect(1,jpw)*ph3d(1,ipw,iaph3d)-vect(2,jpw)*ph3d(2,ipw,iaph3d))
       scali(ipw)=(vect(2,jpw)*ph3d(1,ipw,iaph3d)+vect(1,jpw)*ph3d(2,ipw,iaph3d))
     end do
!$OMP END DO

!$OMP SINGLE
     if (ipw0==2) then
       scalr(1)=half*vect(1,1+ipwshft)*ph3d(1,1,iaph3d)
       scali(1)=half*vect(1,1+ipwshft)*ph3d(2,1,iaph3d)
     end if
!$OMP END SINGLE

!    --------------------------------------------------------------------
!    ALL CHOICES:
!    Accumulate Gx
!    --------------------------------------------------------------------

     if (choice>=0) then

!      Step (2) : Compute scal(lmn) = Sum_g[scal(g).f_nl(g).Y_lm(g)]
       if (use_dgemv) then
         call DGEMV('T',npw,nlmn,1.0_DP,ffnl_loc,npw,scalr,1,0.0_DP,scalr_lmn,1)
         call DGEMV('T',npw,nlmn,1.0_DP,ffnl_loc,npw,scali,1,0.0_DP,scali_lmn,1)
       else
         do ilmn=1,nlmn
!           do ipw=1,npw
!             scalr_lmn(ilmn) = scalr_lmn(ilmn) + scalr(ipw) * ffnl_loc(ipw,ilmn)
!             scali_lmn(ilmn) = scali_lmn(ilmn) + scali(ipw) * ffnl_loc(ipw,ilmn)
!           end do
!$OMP SINGLE
           buffer_r = 0.0_DP
           buffer_i = 0.0_DP
!$OMP END SINGLE
!$OMP DO REDUCTION(+:buffer_r,buffer_i)
           do ipw=1,npw
             buffer_r = buffer_r + scalr(ipw) * ffnl_loc(ipw,ilmn)
             buffer_i = buffer_i + scali(ipw) * ffnl_loc(ipw,ilmn)
           end do
!$OMP END DO
!$OMP SINGLE
           scalr_lmn(ilmn) = buffer_r
           scali_lmn(ilmn) = buffer_i
!$OMP END SINGLE
         end do
       end if
!      Step (3) : Compute gx(lmn) = 4pi/sqrt(vol) (i)^l scal(lmn)
!$OMP SINGLE
       if (cplex==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * cmplx(scalr_lmn(ilmn),scali_lmn(ilmn),kind=DP)
           gx(1,ilmn,ia,ispinor) = real(ctmp)
           gx(2,ilmn,ia,ispinor) = aimag(ctmp)
         end do
       else
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * cmplx(scalr_lmn(ilmn),scali_lmn(ilmn),kind=DP)
           gx(1,ilmn,ia,ispinor) = real(ctmp)
         end do
       end if
!$OMP END SINGLE

     end if

   end do ! End loop on atoms

 end do !  End loop on spinorial components
!$OMP END PARALLEL

!Deallocate temporary space
 ABI_FREE(scalr)
 ABI_FREE(scali)
 ABI_FREE(scalr_lmn)
 ABI_FREE(scali_lmn)

!Has to reduce arrays in case of FFT parallelization
 if (mpi_enreg%nproc_fft>1) then
   call timab(48,1,tsec)
   if (choice>=0) then
     call xmpi_sum(gx,mpi_enreg%comm_fft,ierr)
   end if
   call timab(48,2,tsec)
 end if

end subroutine opernla_ylm_mv
!!***

end module m_opernla_ylm_mv
!!***
