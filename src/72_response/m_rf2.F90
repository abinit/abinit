!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rf2
!! NAME
!! m_rf2
!!
!! FUNCTION
!! This module defines structures and provides procedures used to compute the 2nd order Sternheimer 
!! equation with respect to the wave vector perturbation.
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2016 ABINIT group (LB)
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

MODULE m_rf2

 use defs_basis
 use m_profiling_abi

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_rf2/rf2_t
!! NAME
!!  rf2_t
!!
!! FUNCTION
!!  Datatype gathering all needed data
!!
!! SOURCE

 type,public :: rf2_t

!scalars
  integer :: ndir ! number of directions to consider
  integer :: nband_k ! number of bands
  integer :: size_wf ! number of coeffs in a wavefunction

!arrays
  integer :: iperts(2) ! perturbations to apply
   ! ipert = natom + 10 :
   !   iperts(1) = iperts(2) = natom+1
   !
   ! ipert = natom + 11 :
   !   iperts(1) = natom+1
   !   iperts(2) = natom+2

  integer :: idirs(2) ! directions of the perturbations (ndir=1 : idirs(1)=idirs(2) , ndir=2 : idirs(1)/=idirs(2))

  real(dp),allocatable :: RHS_Stern(:,:)
   ! Right-hand side of the 2nd order Sternheimer equation, for every bands.
   ! Namely, for a band "n" :
   ! |(RHS_Stern)_n> = (H^(2)-epsilon_n S^(2)) |u^(0)_n> + 2(H^(1)-epsilon_n S^(1))|u^(1)_n>
   !                       - 2 sum_m ( lambda_mn^(1) ( S^(1) |u^(0)_m> + S^(0) |u^(1)_m> )
   ! ( /!\ : in the sum, the index m runs over occupied bands )
   ! where :
   !  - epsilon_n = eigenvalue of the GS Hamiltonian
   !  - lambda_mn^(1) = <u^(0)_m| (H^(1)-(epsilon_n+epsilon_m)/2 S^(1) |u^(0)_n> (1st order Lagrange multiplier)
   !  - X^(2) = d^2X/(dk_dir1 dk_dir2)
   !  - 2X^(1)Y^(1) = dX/dk_dir1 dY/dk_dir2 + dX/dk_dir2 dY/dk_dir1
   ! Note : P_c^* will be apply to |(RHS_Stern)_n> in dfpt_cgwf.F90
   ! **
   ! Computed in "rf2_init"

  real(dp),allocatable :: amn(:,:)
   ! Scalar needed for the "orhtonormalization condition", see "dcwavef(:,:)"
   ! Namely :
   ! A_mn = <u^(0)_m| S^(2) |u^(0)_n> + 2 <u^(1)_m| S^(0) |u^(1)_n>
   !    + 2 <u^(1)_m| S^(1) |u^(0)_n> + 2 <u^(0)_m| S^(1) |u^(1)_n>
   ! **
   ! Computed in "rf2_init", stored only for testing purposes

  real(dp),allocatable :: dcwavef(:,:)
   ! Vector needed to enforce the "orhtonormalization condition" on 2nd order wave functions.
   ! Namely :
   ! |dcwavef_n> = -1/2 sum_m A_mn |u^(0)_m>
   ! **
   ! Computed in "rf2_init"

  real(dp),allocatable :: lambda_mn(:,:)
   ! 2nd order Lagrange multiplier.
   ! Namely :
   ! lambda_mn =  <u^(0)_m| H^(2) |u^(0)_n> + 2 <u^(1)_m| H^(0) |u^(1)_n>
   !          + 2 <u^(1)_m| H^(1) |u^(0)_n> + 2 <u^(0)_m| H^(1) |u^(1)_n>
   !          - A_mn * (epsilon_m + epsilon_n) / 2 
   ! **
   ! Computed in "rf2_init"

 end type rf2_t

 public :: rf2_getidirs
 public :: rf2_accumulate_bands
 public :: rf2_apply_hamiltonian
 public :: rf2_destroy

!!***

!----------------------------------------------------------------------

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_rf2/rf2_getidirs
!! NAME
!! rf2_getidirs
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_looppert,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_getidirs(idir,idir1,idir2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rf2_getidirs'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir
 integer,intent(out) :: idir1,idir2

! *************************************************************************

 select case(idir)
!  Diagonal terms :
   case(1)
     idir1 = 1
     idir2 = 1
   case(2)
     idir1 = 2
     idir2 = 2
   case(3)
     idir1 = 3
     idir2 = 3
!  Upper triangular terms :
   case(4)
     idir1 = 2
     idir2 = 3
   case(5)
     idir1 = 1
     idir2 = 3
   case(6)
     idir1 = 1
     idir2 = 2
!  Lower triangular terms :
   case(7)
     idir1 = 3
     idir2 = 2
   case(8)
     idir1 = 3
     idir2 = 1
   case(9)
     idir1 = 2
     idir2 = 1
 end select

end subroutine rf2_getidirs
!!***

!----------------------------------------------------------------------

!!****f* m_rf2/rf2_accumulate_bands
!! NAME
!! rf2_init
!!
!! FUNCTION
!! Compute the scalar product < vi | v1j > and add it to rf2%lambda_mn
!! If necessary, compute < vi | v2j > and add it to rf2%amn.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_accumulate_bands(rf2,choice,gs_hamkq,mpi_enreg,iband,idir1,idir2,ipert1,ipert2,&
                                 jband,print_info,vi,v1j,v2j)

 use defs_basis
 use defs_abitypes
 use m_hamiltonian
 use m_cgtools

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rf2_accumulate_bands'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,iband,idir1,idir2,ipert1,ipert2,jband,print_info
 type(rf2_t),intent(inout) :: rf2
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in) :: vi(2,rf2%size_wf),v1j(2,rf2%size_wf),v2j(2,rf2%size_wf)
 
!Local variables ---------------------------------------
!scalars
 integer :: nband_k,size_wf
 real(dp) :: dotr,doti,factor
 character(len=500) :: msg
 character(len=15) :: bra_i,ket_j,op1,op2
 character(len=2) :: pert1,pert2

! *************************************************************************

 nband_k = rf2%nband_k
 size_wf = rf2%size_wf
 factor = one
 if(rf2%ndir==1 .and. choice /= 3) factor = two ! in order to not compute same terms twice

 call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v1j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

 if(print_info/=0) then
   if (ipert1 == gs_hamkq%natom+1) then
     pert1 = "dk"
   else
     pert1 = "dE"
   end if
   if (ipert2 == gs_hamkq%natom+1) then
     pert2 = "dk"
   else
     pert2 = "dE"
   end if
   select case (choice)
    case(1)
      write(bra_i,'(2a,2(i1,a))') ' < du/',pert2,idir2,'(',iband,') | '
      write(ket_j,'(2a,2(i1,a))') ' | du/',pert1,idir1,'(',jband,') > '
      write(op1,'(a)')            '     H^(0)     '
      write(op2,'(a)')            '     S^(0)     '
    case(2)
      write(bra_i,'(2a,2(i1,a))') ' < du/',pert2,idir2,'(',iband,') | '
      write(ket_j,'(a,i1,a)')     ' | u^(0) (',jband,') > '
      write(op1,'(2a,i1,a)')      '     dH/',pert1,idir1,'    '
      write(op2,'(2a,i1,a)')      '     dS/',pert1,idir1,'    '
    case(3)
      write(bra_i,'(a,i1,a)')     ' < u^(0) (',iband,') | '
      write(ket_j,'(a,i1,a)')     ' | u^(0) (',jband,') > '
      write(op1,'(2(2a,i1),a)')   'd^2H/(',pert1,idir1,' ',pert2,idir2,')'
      write(op2,'(2(2a,i1),a)')   'd^2S/(',pert1,idir1,' ',pert2,idir2,')'
    case(4)
      write(bra_i,'(a,i1,a)')     ' < u^(0) (',iband,') | '
      write(ket_j,'(2a,2(i1,a))') ' | du/',pert1,idir1,'(',jband,') > '
      write(op1,'(2a,i1,a)')      '     dH/',pert2,idir2,'    '
      write(op2,'(2a,i1,a)')      '     dS/',pert2,idir2,'    '
   end select
   write(msg,'(3a,2(a,es22.13E3))') bra_i,op1,ket_j,' = ',dotr,',',doti
   call wrtout(std_out,msg)
 end if

 rf2%lambda_mn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + rf2%lambda_mn(:,iband+(jband-1)*nband_k)

 if (choice == 1 .or. gs_hamkq%usepaw==1) then
   call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v2j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   if(print_info/=0) then
     write(msg,'(3a,2(a,es22.13E3))') bra_i,op2,ket_j,' = ',dotr,',',doti
     call wrtout(std_out,msg)
   end if

   rf2%amn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + rf2%amn(:,iband+(jband-1)*nband_k)

 end if ! end choice

end subroutine rf2_accumulate_bands
!!***

!----------------------------------------------------------------------

!!****f* m_rf2/rf2_apply_hamiltonian
!! NAME
!! rf2_apply_hamiltonian
!!
!! FUNCTION
!! Apply the KS Hamiltonian (or derivative) to an input wave function.
!! If asked, it also does some checks.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_apply_hamiltonian(cg_jband,rf2,eig0_jband,eig1_k_jband,jband,gs_hamkq,gvnl1,&
                                  idir,ipert,mpi_enreg,print_info,prtvol,rf_hamk_idir,&
                                  work1,work2,work3)

 use defs_basis
 use defs_abitypes
 use m_hamiltonian
 use m_cgtools

 use m_pawcprj,  only : pawcprj_type,pawcprj_set_zero

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rf2_apply_hamiltonian'
 use interfaces_14_hidewrite
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: print_info,prtvol,idir,ipert,jband
 type(rf2_t),intent(in) :: rf2
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamk_idir
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 real(dp),intent(in) :: cg_jband(2,rf2%size_wf*print_info*rf2%nband_k,2)
 real(dp),intent(in) :: eig0_jband(rf2%nband_k),eig1_k_jband(2*rf2%nband_k)
 real(dp),intent(inout) :: gvnl1(2,rf2%size_wf)
 real(dp),intent(inout) :: work1(2,rf2%size_wf),work2(2,rf2%size_wf),work3(2,rf2%size_wf)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: tim_getghc=1,tim_getgh1c=1,tim_getgh2c=1 ! to change
! GETGH1C/GETGH2C OPTIONS :
 integer,parameter :: berryopt=0 ! no berry phase
! integer,parameter :: optlocal=ipert-gs_hamkq%natom-1
 integer,parameter :: optnl=2 ! non local terms must be computed here
! integer,parameter :: usevnl=0 ! non-local potential is not computed yet, to be computed here
! integer,parameter :: opt_gvnl1=ipert-gs_hamkq%natom-1 ! gvnl1 is used as ouput only
! ****************
 integer :: cpopt,iband,natom,nband_k,sij_opt,size_wf,optlocal,opt_gvnl1,usevnl
 real(dp) :: dotr,doti,dotr2,doti2,tol_test
 character(len=500) :: msg
 
!arrays
 type(pawcprj_type) :: dummy_cwaveprj(1,1)
 real(dp) :: dum1(1,1)
 real(dp),allocatable :: gvnlc(:,:),cgj(:,:),iddk(:,:)
 
! *********************************************************************

 usevnl     = 0
 optlocal   = 0
 opt_gvnl1  = 0
 if(ipert==gs_hamkq%natom+2) then
   usevnl = 1
   optlocal = 1
   opt_gvnl1 = 1
 end if
 call pawcprj_set_zero(dummy_cwaveprj)
 dum1 = zero
 
 natom = gs_hamkq%natom
 nband_k = rf2%nband_k
 size_wf = rf2%size_wf
 sij_opt=1;if (gs_hamkq%usepaw==0) sij_opt=0
 tol_test=tol8 !; if (gs_hamkq%usepaw==1) tol_test=tol8

 if (ipert == 0) then
   cpopt=-1 ! nonlop option : <p_lmn|in> (and derivatives) are computed here (and not saved)
   ABI_ALLOCATE(gvnlc,(2,size_wf))
   gvnlc(:,:) = zero
!  Test if < u^(0) | H^(0) | u^(0) > = eps_0
   if(print_info/=0) then
     ABI_ALLOCATE(cgj,(2,size_wf))
     cgj(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     call getghc(cpopt,cgj,dummy_cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,mpi_enreg,1,prtvol,&
                 sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgj,work2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!     write(msg,'(2(a,es22.13E3))') 'RF2 TEST GETGHC (before) dotr = ',dotr,  '    doti = ',doti
!     call wrtout(std_out,msg)
!     write(msg,'(a,es22.13E3)') 'RF2 TEST GETGHC Eig = ',eig0_jband
!     call wrtout(std_out,msg)
     dotr = dotr - eig0_jband(jband)
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol_test) then
       write(msg,'(a,es22.13E3)') 'RF2 TEST GETGHC : NOT PASSED dotr = ',dotr
       call wrtout(std_out,msg)
     end if

   end if ! end tests

   call getghc(cpopt,work1,dummy_cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,mpi_enreg,1,prtvol,&
               sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   ABI_DEALLOCATE(gvnlc)

 else if (ipert == natom+1 .or. ipert == natom+2) then

!  Test if < u^(0) | ( H^(1) - eps^(0) S^(1) ) | u^(0) > = eps^(1)
   if(print_info/=0) then
!     write(msg,'(2(a,i2))') 'RF2 TEST GETGH1 : ipert= ',ipert-natom,' idir = ',idir
!     call wrtout(std_out,msg)
     ABI_ALLOCATE(cgj,(2,size_wf))
     ABI_ALLOCATE(iddk,(2,size_wf))
     cgj(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     iddk(:,:) = zero
     if (ipert == natom+2) iddk(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,2)
     call getgh1c(berryopt,cgj,dummy_cwaveprj,work2,dum1,work3,gs_hamkq,iddk,idir,ipert,zero,&
                  mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)
     do iband=1,nband_k
       cgj(:,:) = cg_jband(:,1+(iband-1)*size_wf:iband*size_wf,1)
       call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgj,work2,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
       if (gs_hamkq%usepaw==1) then
         call dotprod_g(dotr2,doti2,gs_hamkq%istwf_k,size_wf,2,cgj,work3,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
         dotr = dotr - (eig0_jband(iband)+eig0_jband(jband))*dotr2/two
         doti = doti - (eig0_jband(iband)+eig0_jband(jband))*doti2/two
       end if
!       if (ipert == natom+2) then
!         write(msg,'(2(a,es22.13E3))') 'RF2 TEST GETGH1 (dotprod) : dots = ',dotr,',',doti
!         call wrtout(std_out,msg)
!       end if
       dotr = dotr - eig1_k_jband(1+2*(iband-1))
       doti = doti - eig1_k_jband(2+2*(iband-1))
!       if (ipert == natom+2) then
!         write(msg,'(2(a,i1),2(a,es22.13E3))') 'RF2 TEST GETGH1 : Eig1(',iband,',',jband,') = ',&
!                                         eig1_k_jband(1+2*(iband-1)),',',eig1_k_jband(2+2*(iband-1))
!         call wrtout(std_out,msg)
!       end if
       dotr = sqrt(dotr**2+doti**2)
       if (dotr > tol_test) then
         write(msg,'(2(a,i2),a,es22.13E3)') 'RF2 TEST GETGH1 : jband=',jband,' iband=',iband,&
                                                                       ' NOT PASSED dotr = ',dotr
         call wrtout(std_out,msg)
       end if
     end do ! end iband
     ABI_DEALLOCATE(cgj)
     ABI_DEALLOCATE(iddk)
   end if ! end tests

   call getgh1c(berryopt,work1,dummy_cwaveprj,work2,dum1,work3,gs_hamkq,gvnl1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)

 else if (ipert == natom+10) then
!    Compute  : d^2H/(dk_dir1 dk_dir2)|u^(0)>  (in work2)
!    and      : d^2S/(dk_dir1 dk_dir2)|u^(0)>  (in work3)
     call getgh2c(work1,dummy_cwaveprj,work2,work3,gs_hamkq,dum1,idir,ipert,zero,&
                   mpi_enreg,optlocal,optnl,rf_hamk_idir,sij_opt,tim_getgh2c,usevnl)
 
 end if

end subroutine rf2_apply_hamiltonian
!!***

!----------------------------------------------------------------------

!!****f* m_rf2/rf2_destroy
!! NAME
!! rf2_init
!!
!! FUNCTION
!!  Free all allocated arrays in a rf2_t object
!!
!! INPUTS
!! rf2 : the rf2_t object to free
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_destroy(rf2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rf2_destroy'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(rf2_t),intent(inout) :: rf2
!arrays

!Local variables ---------------------------------------
!scalars

! *************************************************************************

 if (allocated(rf2%RHS_Stern)) then
   ABI_DEALLOCATE(rf2%RHS_Stern)
 end if
 if (allocated(rf2%dcwavef)) then
   ABI_DEALLOCATE(rf2%dcwavef)
 end if
 if (allocated(rf2%amn)) then
   ABI_DEALLOCATE(rf2%amn)
 end if
 if (allocated(rf2%lambda_mn)) then
   ABI_DEALLOCATE(rf2%lambda_mn)
 end if

end subroutine rf2_destroy
!!***

END MODULE m_rf2
