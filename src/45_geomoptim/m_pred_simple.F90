!!****m* ABINIT/m_pred_simple
!! NAME
!!  m_pred_simple
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, JCC, SE)
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

module m_pred_simple

 use defs_basis
 use m_abicore
 use m_abimover
 use m_abihist

 use m_geometry,  only : fcart2fred, xred2xcart

 implicit none

 private
!!***

 public :: pred_simple
 public :: prec_simple
!!***

contains
!!***

!!****f* ABINIT/pred_simple
!! NAME
!! pred_simple
!!
!! FUNCTION
!! Ionmov predictors (4 & 6) Internal to SCFV
!!
!! IONMOV 4:
!! Conjugate gradient algorithm for simultaneous optimization
!! of potential and ionic degrees of freedom. It can be used with
!! iscf=2 and iscf=5 or 6
!!
!! IONMOV 5:
!! Simple relaxation of ionic positions according to (converged)
!! forces. Equivalent to ionmov=1 with zero masses, albeit the
!! relaxation coefficient is not vis, but iprcfc.
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preditor
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine pred_simple(ab_mover,hist,iexit)

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: iexit

!Local variables-------------------------------
!scalars
 integer  :: ihist_next,kk,jj

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   return
 end if

!All the operations are internal to scfcv.F90

!XRED, FCART and VEL
 ihist_next = abihist_findIndex(hist,+1)
 do kk=1,ab_mover%natom
   do jj=1,3
     hist%xred(jj,kk, ihist_next)=hist%xred (jj,kk,hist%ihist)
     hist%fcart(jj,kk,ihist_next)=hist%fcart(jj,kk,hist%ihist)
     hist%vel(jj,kk,ihist_next)=hist%vel(jj,kk,hist%ihist)
   end do
 end do

!ACELL
 do jj=1,3
   hist%acell(jj,ihist_next)=hist%acell(jj,hist%ihist)
 end do

!RPRIMD
 do kk=1,3
   do jj=1,3
     hist%rprimd(jj,kk,ihist_next)=hist%rprimd(jj,kk,hist%ihist)
   end do
 end do

 hist%ihist=ihist_next

end subroutine pred_simple
!!***

!!****f* ABINIT/prec_simple
!! NAME
!! prec_simple
!!
!! FUNCTION
!! Simple preconditioner, compute the force constant matrix
!! using the Badger's rule:
!!
!!                F=A/(r-B)^3
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information needed by the preconditioner
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces acell, rprimd, stresses
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      bonds_free,dsyev,dsysv,fcart2fred,make_bonds_new,xred2xcart
!!
!! SOURCE

subroutine prec_simple(ab_mover,forstr,hist,icycle,itime,iexit)

 use m_linalg_interfaces
 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,itime,icycle
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(in) :: hist
 type(abiforstr),intent(inout) :: forstr

!Local variables-------------------------------
!scalars
 integer :: period,ii,jj,index,kk,ksub,jsub
 integer :: info,lwork,new_order_forces
 real(dp) :: Z,badgerfactor,lambda,sigma,val_rms
 integer,save :: order_forces
 logical :: Compute_Matrix
!arrays
 type(go_bonds) :: bonds
 integer,allocatable :: periods(:,:)
 integer,allocatable :: iatoms(:,:)
 integer  :: ipiv(3*ab_mover%natom)
 real(dp) :: xcart(3,ab_mover%natom)
 real(dp) :: fcart(3,ab_mover%natom)
 real(dp) :: B(3*ab_mover%natom)
 real(dp) :: rprimd(3,3)
 real(dp) :: w(3*ab_mover%natom)
 real(dp),allocatable :: matrix_tmp(:,:)
 real(dp),allocatable :: work(:)
 real(dp) :: badger(6,6)
 real(dp),allocatable,save :: matrix(:,:)
 character(len=18)   :: fmt

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if (iexit/=0)then
   if(allocated(matrix))then
     ABI_DEALLOCATE(matrix)
   endif
   return
 end if

!##########################################################
!### 01. Show the Precondition parameters, set the badger
!###     matrix.

 write(std_out,*) 'Precondition option',ab_mover%goprecon
 write(std_out,*) 'Precondition parameters',ab_mover%goprecprm
 lambda=ab_mover%goprecprm(1)

 badger(:,:)=reshape( (/ -0.2573, 0.3401, 0.6937, 0.7126, 0.8335, 0.9491,&
& 0.3401, 0.9652, 1.2843, 1.4725, 1.6549, 1.7190,&
& 0.6937, 1.2843, 1.6925, 1.8238, 2.1164, 2.3185,&
& 0.7126, 1.4725, 1.8238, 2.0203, 2.2137, 2.5206,&
& 0.8335, 1.6549, 2.1164, 2.2137, 2.3718, 2.5110,&
& 0.9491, 1.7190, 2.3185, 2.5206, 2.5110, 0.0000 /), (/ 6, 6/) )

 write(fmt,'(a1,i4,a5)') '(',3*ab_mover%natom,'f8.3)'

!##########################################################
!### 02. Take the coordinates and cell parameters from HIST

 rprimd(:,:)=hist%rprimd(:,:,hist%ihist)
 fcart(:,:)=hist%fcart(:,:,hist%ihist)
 call xred2xcart(ab_mover%natom,rprimd,xcart,hist%xred(:,:,hist%ihist))

!##########################################################
!### 03. Decide based on kind of preconditioner if
!###     a new matrix should be computed

 if (ab_mover%goprecon==2)then

   val_rms=0.0
   do kk=1,ab_mover%natom
     do jj=1,3
       val_rms=val_rms+fcart(jj,kk)**2
     end do
   end do
   val_rms=sqrt(val_rms/dble(ab_mover%natom))
   new_order_forces=int(log(val_rms)/log(10.0))
 end if

 if (itime==1.and.icycle==1)then
   Compute_Matrix=.TRUE.
   order_forces=new_order_forces
   if (allocated(matrix))  then
     ABI_DEALLOCATE(matrix)
   end if

   ABI_ALLOCATE(matrix,(3*ab_mover%natom,3*ab_mover%natom))
 else
   Compute_Matrix=.FALSE.
   if ((ab_mover%goprecon==2).and.(order_forces.gt.new_order_forces)) then
     Compute_Matrix=.TRUE.
     order_forces=new_order_forces
   end if
   if (ab_mover%goprecon==3) Compute_Matrix=.TRUE.
 end if

!##########################################################
!### 04. Compute a new precondition matrix if required

 if (Compute_Matrix)then

!  Fix the tolerance for create a bond
   bonds%tolerance=1.35
   bonds%nbonds=1

!  Allocate the arrays with exactly the rigth nbonds
   ABI_ALLOCATE(bonds%bond_vect,(3,bonds%nbonds))
   ABI_ALLOCATE(bonds%bond_length,(bonds%nbonds))
   ABI_ALLOCATE(bonds%indexi,(ab_mover%natom,bonds%nbonds))
   ABI_ALLOCATE(bonds%nbondi,(ab_mover%natom))

!  Compute the bonds
   call make_bonds_new(bonds,ab_mover%natom,ab_mover%ntypat,rprimd,&
&   ab_mover%typat,xcart,ab_mover%znucl)

   write(std_out,'(a,a,63a,a)') ch10,'---PRECONDITIONER',('-',kk=1,63),ch10

!  For all bonds detect wich atoms are involved
!  and wich period they coprrespond in the periodic table
   if (bonds%nbonds>0)then

     ABI_ALLOCATE(periods,(2,bonds%nbonds))
     ABI_ALLOCATE(iatoms,(2,bonds%nbonds))
     periods(:,:)=0
     iatoms(:,:)=0

     write(std_out,'(a)') 'Bond of Atom | Bond Number | Index'

     do ii=1,ab_mover%natom
       Z=ab_mover%znucl(ab_mover%typat(ii))
       if (Z==1)then
         period=1
       elseif ((Z>1).and.(Z<10))then
         period=2
       elseif ((Z>10).and.(Z<18))then
         period=3
       elseif ((Z>18).and.(Z<36))then
         period=4
       elseif ((Z>36).and.(Z<54))then
         period=5
       elseif ((Z>55).and.(Z<86))then
         period=6
       else
!        Here are the cases for atoms larger than Fr(87) and
!        All the noble gases He-Rn
         period=-1
       end if

       do jj=1,bonds%nbondi(ii)
         index=bonds%indexi(ii,jj)

         write(std_out,'(i6,a,i6,a,i4)') ii,'       |',jj,'       |',index

!        The first atom should have index=0
!        To make easy fill the matrix using its
!        index

         if (index>0)then
           periods(1,index)=period
           iatoms(1,index)=ii
         elseif (index<0) then
           periods(2,-index)=period
           iatoms(2,-index)=ii
         end if
       end do
     end do

     write(std_out,'(a)') ch10

   end if

!  For all bonds compute the 3x3 matrix and fill also the big matrix

   matrix(:,:)=0.0_dp
   do ii=1,bonds%nbonds

     write(std_out,*) 'Bond number:',ii
     if (iatoms(1,ii)>0 .and. iatoms(2,ii)>0) then
       write(std_out,*) 'Between atoms:',iatoms(1,ii),' and ',iatoms(2,ii)
       badgerfactor=badger(periods(1,ii),periods(2,ii))
       write(std_out,*) 'Periods of atoms:',periods(1,ii),' and ',periods(2,ii)
       write(std_out,*) 'Badger factor:',badgerfactor

!      Compute the diadic product and
!      Insert the matrix into the big one
       do jj=1,3
         do kk=1,3
!          The non diagonal elements
           jsub=3*(iatoms(1,ii)-1)+jj
           ksub=3*(iatoms(2,ii)-1)+kk
           matrix(jsub,ksub)=matrix(jsub,ksub)-&
&           badgerfactor*bonds%bond_vect(jj,ii)*bonds%bond_vect(kk,ii)

           jsub=3*(iatoms(2,ii)-1)+jj
           ksub=3*(iatoms(1,ii)-1)+kk
           matrix(jsub,ksub)=matrix(jsub,ksub)-&
&           badgerfactor*bonds%bond_vect(jj,ii)*bonds%bond_vect(kk,ii)

!          The diagonal blocks
           jsub=3*(iatoms(1,ii)-1)+jj
           ksub=3*(iatoms(1,ii)-1)+kk
           matrix(jsub,ksub)=matrix(jsub,ksub)+&
&           badgerfactor*bonds%bond_vect(jj,ii)*bonds%bond_vect(kk,ii)

           jsub=3*(iatoms(2,ii)-1)+jj
           ksub=3*(iatoms(2,ii)-1)+kk
           matrix(jsub,ksub)=matrix(jsub,ksub)+&
&           badgerfactor*bonds%bond_vect(jj,ii)*bonds%bond_vect(kk,ii)

         end do !do kk=1,3
       end do !do jj=1,3

     end if

   end do

   if (bonds%nbonds>0)then
     ABI_DEALLOCATE(periods)
     ABI_DEALLOCATE(iatoms)
   end if

   call bonds_free(bonds)

   if (3*ab_mover%natom<100)then
!    Visualize the matrix
     do jj=1,3*ab_mover%natom
       write (std_out,fmt) matrix(jj,:)
     end do
   end if

   ABI_ALLOCATE(matrix_tmp,(3*ab_mover%natom,3*ab_mover%natom))
   
   matrix_tmp(:,:)=matrix(:,:)
   !write(*,*)"matrix_tmp",matrix_tmp

   ABI_ALLOCATE(work,(1))
   lwork=-1
   call DSYEV('V', 'U', 3*ab_mover%natom, matrix_tmp, 3*ab_mover%natom, w , work, lwork, info )
   lwork=work(1)
   write(std_out,*) '[DSYEV] Recommended lwork=',lwork
   ABI_DEALLOCATE(work)
   ABI_ALLOCATE(work,(lwork))
   call DSYEV('V', 'U', 3*ab_mover%natom, matrix_tmp, 3*ab_mover%natom, w , work, lwork, info )
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(matrix_tmp)

   write(std_out,*) 'DSYEV info=',info
   write(std_out,*) 'Eigenvalues:'
   write(std_out,fmt) w(:)

   sigma=0
   do jj=1,3*ab_mover%natom
     sigma=max(w(jj),sigma)
   end do

   matrix=lambda*matrix

   write(std_out,*) ch10
   do ii=1,3*ab_mover%natom
     matrix(ii,ii)=matrix(ii,ii)+(1-lambda)*sigma
   end do

 end if ! if (Compute_Matrix)

!##########################################################
!### 05. Use the precondition matrix to compute new residuals

 B=reshape(fcart,(/ 3*ab_mover%natom /))

 if (3*ab_mover%natom<100)then
!  Visualize the matrix
   do jj=1,3*ab_mover%natom
     write (std_out,fmt) matrix(jj,:)
   end do
 end if

!call dsysv( uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info )
!MGNAG FIXME: This call causes a floating point exception if NAG+MKL
 ABI_ALLOCATE(work,(1))
 lwork=-1
 call DSYSV( 'U', 3*ab_mover%natom, 1, matrix,&
& 3*ab_mover%natom, ipiv, B, 3*ab_mover%natom, work, lwork, info )

 lwork=work(1)
 write(std_out,*) '[DSYSV] Recomended lwork=',lwork
 ABI_DEALLOCATE(work)
 ABI_ALLOCATE(work,(lwork))
 call DSYSV( 'U', 3*ab_mover%natom, 1, matrix,&
& 3*ab_mover%natom, ipiv, B, 3*ab_mover%natom, work, lwork, info )
 ABI_DEALLOCATE(work)

 write(std_out,*) 'DSYSV info=',info
 write(std_out,*) 'Solution:'
 write(std_out,fmt) B(:)

 forstr%fcart=reshape(B,(/ 3, ab_mover%natom /) )
 call fcart2fred(forstr%fcart,forstr%fred,rprimd,ab_mover%natom)

end subroutine prec_simple
!!***

end module m_pred_simple
!!***
