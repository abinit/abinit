!!****m* ABINIT/m_opernla_ylm_blas
!! NAME
!!  m_opernla_ylm_blas
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MT)
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

module m_opernla_ylm_blas

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

 public :: opernla_ylm_blas
!!***

contains
!!***

!!****f* ABINIT/opernla_ylm_blas
!! NAME
!! opernla_ylm_blas
!!
!! FUNCTION
!! For a given wave-function |c>, get all projected scalars
!! <p_lmn|c> where |p_lmn> are non-local projectors
!!   With:
!!   <p_lmn|c>=4pi/sqrt(vol) (i)^l Sum_g[c(g).f_nl(g).Y_lm(g).exp(2pi.i.g.R)]
!!
!! INPUTS
!!  choice=chooses possible output:
!!         if choice>=0: compute projected scalars
!!         if choice<0: same as choice>0 but use already computed projected scalars
!!         if ABS(choice)>1, then compute additional quantities:
!!           2: compute projected scalars and derivatives wrt atm pos.
!!           3: compute projected scalars and derivatives wrt strains
!!           22: compute projected scalars and 2nd derivatives wrt atm pos. and q-vector.
!!           23: compute projected scalars, derivatives wrt atm pos. and derivatives wrt strains
!!           25: compute projected scalars and 3rd derivatives wrt atm pos. and two q-vectors.
!!           4, 24: compute projected scalars, derivatives wrt atm pos.
!!                  and 2nd derivatives wrt atm pos.
!!           33: compute projected scalars and 2nd derivatives wrt strain and q-vector.
!!           5,51,52: compute projected scalars and derivatives wrt wave vector k
!!           53: compute projected scalars and derivatives wrt wave vector k in direction idir+1 and idir-1
!!           54: compute projected scalars, deriv. wrt atm pos., deriv. wrt wave vector k
!!               and 2nd derivatives wrt right wave vector k and atm pos.
!!           55: compute projected scalars, deriv. strains, deriv. wrt wave vector k
!!               and 2nd derivatives wrt right wave vector k and strain
!!           6: compute projected scalars, derivatives wrt atm pos., derivatives wrt strains,
!!              2nd derivatives wrt 2 strains and derivatives wrt strain and atm pos.
!!           7: not available
!!           8: compute projected scalars, derivatives wrt wave vector k
!!              and 2nd derivatives wrt 2 wave vectors k
!!  cplex=1 if <p_lmn|c> scalars are real or pure imaginary (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2) or (choice=22,signs=2)
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!                        - strain component (1:9) in the case (choice=33,signs=2) 
!!                        - (1:9) components to specify the atom to be moved and the second q-gradient 
!!                          direction in the case (choice=25,signs=2)
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  istwf_k=option parameter that describes the storage of wfs
!!  kpg(npw,nkpg)=(k+G) components          for ikpg=1...3   (if nkpg=3 or 9)
!!       [(k+G)_a].[(k+G)_b] quantities for ikpg=4...9   (if nkpg=9)
!!       (k+G) Cartesian components for choice==33
!!  matblk=dimension of the array ph3d
!!  mpi_enreg=information about MPI parallelization
!!  ndgxdt=second dimension of dgxdt
!!  nd2gxdt=second dimension of d2gxdt
!!  nincat=number of atoms in the subset here treated
!!  nkpg=second dimension of array kpg (0, 3 or 9)
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  [qdir]= optional, direction of the q-gradient (only for choice=22 choice=25 and choice=33) 
!!  signs=chooses possible output:
!!   signs=1: compute derivatives in all directions
!!   signs=2: compute derivative in direction IDIR only
!!            compatible only with 1st-order derivatives and "single" derivatives
!!  ucvol=unit cell volume (bohr^3)
!!  vect(2,npw*my_nspinor)=starting vector in reciprocal space
!!
!! OUTPUT
!!  if (choice>1) dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)=
!!     gradients of projected scalars wrt coords  (choice=2, 23, 4, 54, 6)
!!                                    wrt strains (choice=3, 23, 55)
!!                                    wrt k wave vect. (choice=5, 51, 52, 53, 54, 55, 8)
!!                                    wrt coords and q vect (choice=22)
!!                                    wrt coords and two q vects (choice=25)
!!                                    wrt strains and q vect (choice=33)
!!  if (choice=4, 24, 33, 54, 55, 6, 8) d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)=
!!     2nd grads of projected scalars wrt 2 coords (choice=4 or 24)
!!                                    wrt coords & k wave vect. (choice=54)
!!                                    wrt strains & k wave vect. (choice=55)
!!                                    wrt coords & strains (choice=6)
!!                                    wrt 2 strains (choice=6)
!!                                    wrt 2 k wave vect. (choice=8)
!!                                    wrt strains and q vect (choice=33)
!!     only compatible with signs=1
!!  cplex_dgxdt(ndgxdt) = used only when cplex = 1
!!             cplex_dgxdt(i) = 1 if dgxdt(1,i,:,:)   is real, 2 if it is pure imaginary
!!  cplex_d2gxdt(nd2gxdt) = used only when cplex = 1
!!             cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
!!
!! SIDE EFFECTS
!!  gx(cplex,nlmn,nincat,nspinor)= projected scalars - input if choice<0, output if choice>=0
!!
!! NOTES
!! 1-The openMP version is different from the standard version:
!!   the standard version is more effifient on one CPU core.
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!!
!! PARENTS
!!      getcprj,nonlop_ylm
!!
!! CHILDREN
!!      timab,xmpi_sum
!!
!! SOURCE

subroutine opernla_ylm_blas(choice,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnl,d2gxdt,dgxdt,ffnl,gx,&
&       ia3,idir,indlmn,istwf_k,kpg,matblk,mpi_enreg,nd2gxdt,ndgxdt,nincat,nkpg,nlmn,&
&       nloalg,npw,nspinor,ph3d,signs,ucvol,vect,qdir)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,dimffnl,ia3,idir,istwf_k,matblk,nd2gxdt
 integer,intent(in) :: ndgxdt,nincat,nkpg,nlmn,npw,nspinor,signs
 integer,intent(in),optional :: qdir
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,nlmn),nloalg(3)
 integer,intent(out) :: cplex_dgxdt(ndgxdt),cplex_d2gxdt(nd2gxdt)
 real(dp),intent(in) :: ffnl(npw,dimffnl,nlmn),kpg(npw,nkpg),ph3d(2,npw,matblk)
 real(dp),intent(in) :: vect(:,:)
 real(dp),intent(out) :: d2gxdt(cplex,nd2gxdt,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor),gx(cplex,nlmn,nincat,nspinor)

!Local variables-------------------------------
!scalars
 integer :: choice_,gama,gamb,gamc,gamd,ia,iaph3d,ialpha,ibeta,idelta,idelgam
 integer :: ierr,igamma,il,ilmn,ipw,ipw0,ipwshft,ispinor,jpw,mu,mu0
 integer :: mua,mub,ndir,nthreads,nua1,nua2,nub1,nub2
 real(dp), parameter :: two_pi2=two_pi*two_pi
 real(dp) :: aux_i,aux_i2,aux_i3,aux_i4
 real(dp) :: aux_r,aux_r2,aux_r3,aux_r4
 real(dp) :: buffer_i1,buffer_i2,buffer_i3,buffer_i4,buffer_i5,buffer_i6
 real(dp) :: buffer_ia,buffer_ib,buffer_ic,buffer_id,buffer_ie,buffer_if
 real(dp) :: buffer_r1,buffer_r2,buffer_r3,buffer_r4,buffer_r5,buffer_r6
 real(dp) :: buffer_ra,buffer_rb,buffer_rc,buffer_rd,buffer_re,buffer_rf
 real(dp) :: kpga,kpgb,kpgc,kpgd,scale,scale2,wt
 logical :: check,parity
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: nalpha(9)=(/1,2,3,3,3,2,2,1,1/),nbeta(9)=(/1,2,3,2,1,1,3,3,2/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 integer,parameter :: ffnl_dir_dat(6)=(/2,3,3,1,1,2/)
 integer,parameter :: idir1(9)=(/1,1,1,2,2,2,3,3,3/),idir2(9)=(/1,2,3,1,2,3,1,2,3/)
 integer :: ffnl_dir(2)
 real(dp) :: tsec(2)
! real(dp),allocatable :: kpg(:,:)
 real(dp),allocatable :: scali(:),scalr(:),scalari(:,:),scalarr(:,:),scalr_lmn(:),scali_lmn(:)
 complex(dp) :: ctmp,cil(4)
! *************************************************************************

 if (choice==-1) return

 if (abs(choice)>1) then
   MSG_ERROR('Only abs(choice)<=0 is available for now.')
 end if
 if (cplex/=2) then
   MSG_ERROR('Only cplex=2 is available for now.')
 end if
 if (istwf_k/=1) then
   MSG_ERROR('Only istwf_k=1 is available for now.')
 end if

!Useful variables
 choice_=abs(choice)
 wt=four_pi/sqrt(ucvol);if (cplex==1) wt=2.d0*wt
 ipw0=1;if (istwf_k==2.and.mpi_enreg%me_g0==1) ipw0=2
 cplex_dgxdt(:)  = 0 ; if (cplex == 1) cplex_dgxdt(:)  = 1
 cplex_d2gxdt(:) = 0 ; if (cplex == 1) cplex_d2gxdt(:) = 1
 nthreads=1

!Check compatibility of options
 if (signs==1) then
!  signs=1, almost all choices
   check=(choice_==0 .or.choice_==1 .or.choice_==2 .or.choice_==3 .or.choice_==4 .or.&
&   choice_==23.or.choice_==24.or.choice_==5 .or.choice_==51.or.choice_==52.or.&
&   choice_==53.or.choice_==54.or.choice_==55.or.&
&   choice_==6 .or.choice_==8 .or.choice_==81)
   ABI_CHECK(check,'BUG: choice not compatible (for signs=1)')
 end if
 if (signs==2) then
!  signs=2,less choices
   check=(choice_== 0.or.choice_== 1.or.choice_== 2.or.choice_==3.or.choice_==5.or.&
&  choice_==22.or.choice_==23.or.choice==25.or.choice_==33.or.choice_==51.or.choice_==52.or.choice_==54.or.choice_==55.or.&
&   choice_== 8.or.choice_==81)
   ABI_CHECK(check,'BUG: signs=2 not compatible with this choice')
 end if
 if (choice_==3.and.signs==2) then
!  1<=idir<=6 is  required when choice=3 and signs=2
   check=(idir>=1.and.idir<=6)
   ABI_CHECK(check,'BUG: choice=3 and signs=2 requires 1<=idir<=6')
 else if (choice_==25.and.signs==2) then
!  1<=idir<=9 is  required when choice=25 and signs=2
   check=(idir>=1.and.idir<=9)
   ABI_CHECK(check,'BUG: choice=25 and signs=2 requires 1<=idir<=9')
 else if (choice_==33.and.signs==2) then
!  1<=idir<=9 is  required when choice=33 and signs=2
   check=(idir>=1.and.idir<=9)
   ABI_CHECK(check,'BUG: choice=33 and signs=2 requires 1<=idir<=9')
 else if ((choice_==8.or.choice_==81.or.choice_==54).and.signs==2) then
!  1<=idir<=9 is required when choice=8 and signs=2
   check=(idir>=1.and.idir<=9)
   ABI_CHECK(check,'BUG: choice=8/81 and signs=2 requires 1<=idir<=9')
 else
!  signs=2 requires 1<=idir<=3 when choice>1
   check=(signs/=2.or.choice<=1.or.(idir>=1.and.idir<=3))
   ABI_CHECK(check,'BUG: signs=2 requires 1<=idir<=3')
 end if

#if defined HAVE_OPENMP
 nthreads=OMP_GET_NUM_THREADS()
#endif

 if (nthreads>1) then
   MSG_ERROR('Only nthreads=1 is available for now.')
 end if

!Allocate work space
 ABI_ALLOCATE(scalr,(npw))
 ABI_ALLOCATE(scali,(npw))
 ABI_ALLOCATE(scalr_lmn,(nlmn))
 ABI_ALLOCATE(scali_lmn,(nlmn))

! i^l
 cil(1) = ( 1.0_DP, 0.0_DP) * wt
 cil(2) = ( 0.0_DP, 1.0_DP) * wt
 cil(3) = (-1.0_DP, 0.0_DP) * wt
 cil(4) = ( 0.0_DP,-1.0_DP) * wt

!Loop on spinorial components
 do ispinor =1,nspinor
   ipwshft=(ispinor-1)*npw

!  Loop on atoms (blocking)
   do ia=1,nincat
     iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1
!    Compute c(g).exp(2pi.i.g.R)
     do ipw=ipw0,npw
       jpw=ipw+ipwshft
       scalr(ipw)=(vect(1,jpw)*ph3d(1,ipw,iaph3d)-vect(2,jpw)*ph3d(2,ipw,iaph3d))
       scali(ipw)=(vect(2,jpw)*ph3d(1,ipw,iaph3d)+vect(1,jpw)*ph3d(2,ipw,iaph3d))
     end do

     if (ipw0==2) then
       scalr(1)=half*vect(1,1+ipwshft)*ph3d(1,1,iaph3d)
       scali(1)=half*vect(1,1+ipwshft)*ph3d(2,1,iaph3d)
     end if

!    --------------------------------------------------------------------
!    ALL CHOICES:
!    Accumulate Gx
!    --------------------------------------------------------------------

     if (choice>=0) then

       call DGEMV('T',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,scalr,1,0.0_DP,scalr_lmn,1)
       call DGEMV('T',npw,nlmn,1.0_DP,ffnl(:,1,:),npw,scali,1,0.0_DP,scali_lmn,1)
!       scalr_lmn(:)=0.0_DP
!       scali_lmn(:)=0.0_DP
!       do ilmn=1,nmln
!         do ipw=1,npw
!           scalr_lmn(ilmn) = scalr_lmn(ilmn) + scalr(ipw) * ffnl(ipw,1,ilmn)
!           scali_lmn(ilmn) = scali_lmn(ilmn) + scali(ipw) * ffnl(ipw,1,ilmn)
!         end do                  
!       end do
       do ilmn=1,nlmn
         il=mod(indlmn(1,ilmn),4)+1
         ctmp = cil(il) * cmplx(scalr_lmn(ilmn),scali_lmn(ilmn),kind=DP)
         gx(1,ilmn,ia,ispinor) = real(ctmp)
         gx(2,ilmn,ia,ispinor) = aimag(ctmp)
       end do

     end if

   end do ! End loop on atoms

 end do !  End loop on spinorial components

!Deallocate temporary space
 ABI_DEALLOCATE(scalr)
 ABI_DEALLOCATE(scali)
 ABI_DEALLOCATE(scalr_lmn)
 ABI_DEALLOCATE(scali_lmn)

!Has to reduce arrays in case of FFT parallelization
 if (mpi_enreg%nproc_fft>1) then
   call timab(48,1,tsec)
   if (choice>=0) then
     call xmpi_sum(gx,mpi_enreg%comm_fft,ierr)
   end if
   if (choice_>1) then
     if (ndgxdt>0) then
       call xmpi_sum(dgxdt,mpi_enreg%comm_fft,ierr)
     end if
     if (nd2gxdt>0) then
       call xmpi_sum(d2gxdt,mpi_enreg%comm_fft,ierr)
     end if
   end if
   call timab(48,2,tsec)
 end if

end subroutine opernla_ylm_blas
!!***

end module m_opernla_ylm_blas
!!***
