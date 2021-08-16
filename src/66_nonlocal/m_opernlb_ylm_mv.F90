!!****m* ABINIT/m_opernlb_ylm_mv
!! NAME
!!  m_opernlb_ylm_mv
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

module m_opernlb_ylm_mv

 use defs_basis
 use m_abicore
 use m_errors
#if defined HAVE_OPENMP
 use OMP_LIB
#endif

 implicit none

 private
!!***

 public :: opernlb_ylm_mv
 integer,public,save :: opernlb_mv_counter = -1
 integer,public,save :: opernlb_mv_dgemv_counter = -1
!!***

contains
!!***

!!****f* ABINIT/opernlb_ylm_mv
!! NAME
!! opernlb_ylm_mv
!!
!! FUNCTION
!! "matrix-vector" alternative implementation of "opernlrb_ylm".
!!
!! * Operate with the non-local part of the hamiltonian,
!!   from projected scalars to reciprocal space.
!! * Operate with the non-local projectors and the overlap matrix,
!!   from projected scalars to reciprocal space.
!!
!!   The input is gxfac (gxfac_sij):
!!   gxfac(lmn) = Sum_l'm'n' D_l'm'n'.<p_l'm'n|c> (or S_l'm'n' for gxfac_sij)
!!   and here we compute :
!!   Sum_lmn <g|p_lmn> gxfac(lmn) = 4pi/sqrt(vol) exp(-2pi.i.g.R) Sum_lmn (-i)^l f_nl(g).Y_lm(g) gxfac(lmn)
!!   Here this is done in 3 steps:
!!   (1) compute for every lmn : gxfac_(lmn) = 4pi/sqrt(vol).(-i)^l.gxfac(lmn)
!!   (2) compute for every g   : scal(g)     = Sum_lmn f_nl(g).Y_lm(g).gxfac_(lmn)
!!   (3) compute for every g   : vect(g)     = exp(-2pi.i.g.R).scal(g)
!!
!!   Step (2) is a real-matrix/complex-vector multiplication, here two options are possible:
!!   - case nloalg(1)=2 : compute the real and imaginary parts separately using two calls of DGMEV
!!   - case nloalg(1)=3 : in order to read the matrix only once, we compute both real and imaginary parts at the same time "by hand"
!!
!!   Depending on the achitecture and the available blas library, one option could be more interesting than an other...
!!
!! INPUTS
!!  choice=chooses possible output (see below)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dimffnl=second dimension of ffnl
!!  ffnl(npw,dimffnl,nlmn)= nonlocal quantities containing nonlocal form factors
!!  gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))= reduced projected scalars related to Sij (overlap)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  matblk=dimension of the array ph3d
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw=number of plane waves in reciprocal space
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! --if (paw_opt=0)
!!    vectout(2,npwout*my_nspinor*ndat)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.
!! --if (paw_opt=0, 1 or 4)
!!    vect(2,npwout*nspinor)=result of the aplication of the concerned operator
!!                or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|V_nonlocal|vect_in>
!!  if (paw_opt=2)
!!    vect(2,npwout*nspinor)=final vector in reciprocal space:
!!      if (choice=1)  <G|V_nonlocal-lamdba.(I+S)|vect_in> (note: not including <G|I|c>)
!! --if (paw_opt=3 or 4)
!!    svect(2,npwout*nspinor)=result of the aplication of Sij (overlap matrix)
!!                  or one of its derivatives to the input vect.:
!!      if (choice=1)  <G|I+S|vect_in> (note: not including <G|I|c>)
!!
!! NOTES
!! 1-No openMP available for now
!! 2-Operate for one type of atom, and within this given type of atom,
!!   for a subset of at most nincat atoms.
!! 3-projector derivatives (abs(choice)>1) are not implemented yet
!!
!! PARENTS
!!      m_nonlop_ylm
!!
!! CHILDREN
!!      dgemv
!!
!! SOURCE

subroutine opernlb_ylm_mv(choice,cplex,cplex_fac,&
&                      dimffnl,ffnl,gxfac,gxfac_sij,&
&                      ia3,indlmn,matblk,nincat,nlmn,nloalg,npw,&
&                      nspinor,paw_opt,ph3d,svect,ucvol,vect)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,cplex,cplex_fac,dimffnl,ia3,matblk,nincat
 integer,intent(in) :: nlmn,npw,nspinor,paw_opt
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) ::  indlmn(6,nlmn),nloalg(3)
 real(dp),intent(in),target :: ffnl(npw,dimffnl,nlmn)
 real(dp),intent(in) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(in) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(in) :: ph3d(2,npw,matblk)
 real(dp),intent(inout) :: svect(:,:),vect(:,:)
!Local variables-------------------------------
!Arrays
!scalars
 logical :: use_dgemv
 integer :: ia,iaph3d
 integer :: il,ilmn,ipw,jpw,ipwshft,ispinor
 real(dp) :: buffer_r,buffer_i,wt
!arrays
! real(dp) :: tsec(2)
 real(dp) :: gxfac_(nlmn,2),gxfacs_(nlmn,2)
 real(dp),allocatable :: scalr(:),scali(:)
 real(dp),pointer :: ffnl_loc(:,:)
 complex(dp) :: ctmp, cil(4)

! *************************************************************************

 DBG_ENTER("COLL")

!Some checks
! nthreads=1
!#if defined HAVE_OPENMP
! nthreads=OMP_GET_NUM_THREADS()
!#endif
! if (nthreads>1) then
!   ABI_ERROR('Only nthreads=1 is available for now.')
! end if

 if (abs(choice)>1) then
   ABI_ERROR('Only abs(choice)<=1 is available for now.')
 end if
 if (nloalg(1)<2.or.nloalg(1)>10) then
   ABI_ERROR('nloalg(1) should be between 2 or 10.')
 end if

 use_dgemv = nloalg(1)==2.or.nloalg(1)==6.or.nloalg(1)==10
 if (use_dgemv) then
   if(opernlb_mv_dgemv_counter>=0) then
     opernlb_mv_dgemv_counter = opernlb_mv_dgemv_counter + 1
     if (paw_opt==4) opernlb_mv_dgemv_counter = opernlb_mv_dgemv_counter + 1
   end if
 else
   if(opernlb_mv_counter>=0) then
     opernlb_mv_counter = opernlb_mv_counter + 1
     if (paw_opt==4) opernlb_mv_counter = opernlb_mv_counter + 1
   end if
 end if

!Inits
 wt=four_pi/sqrt(ucvol)

 ffnl_loc => ffnl(:,1,:)


 ABI_MALLOC(scalr,(npw))
 ABI_MALLOC(scali,(npw))

!$OMP PARALLEL PRIVATE(buffer_r,buffer_i,cil,il,ilmn,ipw,jpw,ctmp), &
!$OMP PRIVATE(ispinor,ipwshft,ia,iaph3d,gxfac_,gxfacs_)

! (-i)^l
 cil(1) = ( 1.0_DP, 0.0_DP) * wt
 cil(2) = ( 0.0_DP,-1.0_DP) * wt
 cil(3) = (-1.0_DP, 0.0_DP) * wt
 cil(4) = ( 0.0_DP, 1.0_DP) * wt

! if (paw_opt/=3) then
!   ABI_MALLOC(gxfac_,(nlmn,2))
! end if
! if (paw_opt>=3) then
!   ABI_MALLOC(gxfacs_,(nlmn,2))
! end if

!Loop on spinorial components
 do ispinor=1,nspinor
   ipwshft=(ispinor-1)*npw

!  Loop on atoms (blocking)
   do ia=1,nincat
     iaph3d=ia;if (nloalg(2)>0) iaph3d=ia+ia3-1
!    Step (1) : scale gxfac with 4pi/sqr(omega).(-i)^l
     if (paw_opt/=3) then
       if (cplex_fac==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * cmplx( gxfac(1,ilmn,ia,ispinor), gxfac(2,ilmn,ia,ispinor), kind=DP )
           gxfac_(ilmn,1) =  real(ctmp)
           gxfac_(ilmn,2) = aimag(ctmp)
         end do
       else if (cplex_fac==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * gxfac(1,ilmn,ia,ispinor)
           gxfac_(ilmn,1) =  real(ctmp)
           gxfac_(ilmn,2) = aimag(ctmp)
         end do
       else
         ABI_BUG('Error : should not be possible to be here')
       end if
     end if

!    Step (1) bis: Scale gxfac_sij with 4pi/sqr(omega).(-i)^l
     if (paw_opt>=3) then
       if (cplex==2) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * cmplx( gxfac_sij(1,ilmn,ia,ispinor), gxfac_sij(2,ilmn,ia,ispinor), kind=DP )
           gxfacs_(ilmn,1) =  real(ctmp)
           gxfacs_(ilmn,2) = aimag(ctmp)
         end do
       else if (cplex==1) then
         do ilmn=1,nlmn
           il=mod(indlmn(1,ilmn),4)+1
           ctmp = cil(il) * gxfac_sij(1,ilmn,ia,ispinor)
           gxfacs_(ilmn,1) =  real(ctmp)
           gxfacs_(ilmn,2) = aimag(ctmp)
         end do
       else
         ABI_BUG('Error : should not be possible to be here')
       end if
     end if

!    Compute <g|Vnl|c> (or derivatives) for each plane wave:
     if (paw_opt/=3) then

!      Step (2) scal(g) = Sum_lmn f_nl(g).Y_lm(g).gxfac_(lmn)
       if (use_dgemv) then
         call DGEMV('N',npw,nlmn,1.0_DP,ffnl_loc,npw,gxfac_(:,1),1,0.0_DP,scalr,1)
         call DGEMV('N',npw,nlmn,1.0_DP,ffnl_loc,npw,gxfac_(:,2),1,0.0_DP,scali,1)
       else
!         scalr(:) = zero
!         scali(:) = zero
!         do ilmn=1,nlmn
!           do ipw=1,npw
!             scalr(ipw) = scalr(ipw) + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,1)
!             scali(ipw) = scali(ipw) + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,2)
!           end do
!         end do
!$OMP DO
         do ipw=1,npw
           buffer_r = zero
           buffer_i = zero
           do ilmn=1,nlmn
             buffer_r = buffer_r + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,1)
             buffer_i = buffer_i + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,2)
!             scalr(ipw) = scalr(ipw) + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,1)
!             scali(ipw) = scali(ipw) + ffnl_loc(ipw,ilmn) * gxfac_(ilmn,2)
           end do
           scalr(ipw) = buffer_r
           scali(ipw) = buffer_i
         end do
!$OMP END DO
       end if

!      Step (3) : vect(g) = exp(-2pi.i.g.R).scal(g)
!$OMP DO
       do ipw=1,npw
         jpw=ipw+ipwshft
         vect(1,jpw)=vect(1,jpw)+scalr(ipw)*ph3d(1,ipw,iaph3d)+scali(ipw)*ph3d(2,ipw,iaph3d)
         vect(2,jpw)=vect(2,jpw)-scalr(ipw)*ph3d(2,ipw,iaph3d)+scali(ipw)*ph3d(1,ipw,iaph3d)
       end do
!$OMP END DO

     end if

!    Compute <g|S|c> (or derivatives) for each plane wave:
     if (paw_opt>=3) then

!      Step (2) (bis) scal(g) = Sum_lmn f_nl(g).Y_lm(g).gxfacs_(lmn)
       if (nloalg(1)==3) then
!         scalr(:) = zero
!         scali(:) = zero
!         do ilmn=1,nlmn
!           do ipw=1,npw
!             scalr(ipw) = scalr(ipw) + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,1)
!             scali(ipw) = scali(ipw) + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,2)
!           end do
!         end do
!$OMP DO
         do ipw=1,npw
           buffer_r = zero
           buffer_i = zero
           do ilmn=1,nlmn
             buffer_r = buffer_r + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,1)
             buffer_i = buffer_i + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,2)
!             scalr(ipw) = scalr(ipw) + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,1)
!             scali(ipw) = scali(ipw) + ffnl_loc(ipw,ilmn) * gxfacs_(ilmn,2)
           end do
           scalr(ipw) = buffer_r
           scali(ipw) = buffer_i
         end do
!$OMP END DO
       else if (nloalg(1)==2) then
         call DGEMV('N',npw,nlmn,1.0_DP,ffnl_loc,npw,gxfacs_(:,1),1,0.0_DP,scalr,1)
         call DGEMV('N',npw,nlmn,1.0_DP,ffnl_loc,npw,gxfacs_(:,2),1,0.0_DP,scali,1)
       end if

!      Step (3) (bis) : svect(g) = exp(-2pi.i.g.R).scal(g)
!$OMP DO
       do ipw=1,npw
         jpw=ipw+ipwshft
         svect(1,jpw)=svect(1,jpw)+scalr(ipw)*ph3d(1,ipw,iaph3d)+scali(ipw)*ph3d(2,ipw,iaph3d)
         svect(2,jpw)=svect(2,jpw)-scalr(ipw)*ph3d(2,ipw,iaph3d)+scali(ipw)*ph3d(1,ipw,iaph3d)
       end do
!$OMP END DO

     end if

!    End loop on atoms
   end do
 end do !  End loop on spinors
! if (paw_opt/=3) then
!   ABI_FREE(gxfac_)
! end if
! if (paw_opt>=3) then
!   ABI_FREE(gxfacs_)
! end if
!$OMP END PARALLEL

 ABI_FREE(scalr)
 ABI_FREE(scali)


 DBG_EXIT("COLL")

end subroutine opernlb_ylm_mv
!!***

end module m_opernlb_ylm_mv
!!***
