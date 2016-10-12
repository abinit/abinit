!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_rf2
!! NAME
!! m_rf2
!!
!! FUNCTION
!! This module defines structures and provides procedures used to compute the 2nd order Sternheimer 
!! equation.
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
!!  Datatype gathering all needed data for the computation of the 2nd order Sternheimer equation.
!!
!! SOURCE

 type,public :: rf2_t

!scalars
  integer :: ndir ! number of directions to consider
  integer :: nband_k ! number of bands
  integer :: size_wf ! number of coeffs in a wavefunction

!arrays
  integer :: iperts(2) ! perturbations to compute
   ! ipert = natom + 10 :
   !   iperts(1) = iperts(2) = natom+1 (wavevector k)
   !
   ! ipert = natom + 11 :
   !   iperts(1) = natom+1 (wavevector k)
   !   iperts(2) = natom+2 (electric field)

  integer :: idirs(2) ! directions of the perturbations (ndir=1 : idirs(1)=idirs(2) , ndir=2 : idirs(1)/=idirs(2))

  real(dp),allocatable :: RHS_Stern(:,:)
   ! Right-hand side of the 2nd order Sternheimer equation, for every bands.
   ! Namely, for a band "n" :
   ! |(RHS_Stern)_n> = (H^(2)-epsilon_n S^(2)) |u^(0)_n> + 2(H^(1)-epsilon_n S^(1))|u^(1)_n>
   !                   - 2 sum_m ( lambda_mn^(1) ( S^(1) |u^(0)_m> + S^(0) |u^(1)_m> )
   !                      ( /!\ : in the sum, the index m runs over occupied bands )
   ! where :
   !  - epsilon_n = eigenvalue of the GS Hamiltonian
   !  - lambda_mn^(1) = <u^(0)_m| (H^(1)-(epsilon_n+epsilon_m)/2 S^(1) |u^(0)_n> (1st order Lagrange multiplier)
   !  - X^(2) = d^2X/(dk_dir1 dlambda_dir2)
   !  - 2X^(1)Y^(1) = dX/dk_dir1 dY/dlambda_dir2 + dX/dlambda_dir2 dY/dk_dir1
   !  lambbda_dir2 = k_dir2 (ipert==natom+10) , E_dir2 (ipert==natom+11)
   ! Note : P_c^* will be apply to |(RHS_Stern)_n> in dfpt_cgwf.F90
   ! **
   ! Computed in "rf2_init"

  real(dp),allocatable :: amn(:,:)
   ! Scalar needed for the "orhtonormalization condition", see "dcwavef(:,:)" below
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
!!  Get directions of the 1st and 2nd perturbations according to the input idir
!!
!! INPUTS
!!  idir : integer between 1 and 9
!!
!! OUTPUT
!!  idir1 : index of the 1st direction (1<=idir1<=3)
!!  idir2 : index of the 2nd direction (1<=idir2<=3)
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_looppert,dfpt_scfcv,rf2_init
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
!!  rf2 : rf2_t object containing all rf2 data
!!  choice : 4 possible values, see "select case (choice)" below
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  mpi_enreg=information about MPI parallelization
!!  iband : band index for vi
!!  idir1  (used only if print_info/=0)  : direction of the 1st perturbation
!!  idir2  (used only if print_info/=0)  : direction of the 2nd perturbation
!!  ipert1 (used only if print_info/=0)  : 1st perturbation
!!  ipert2 (used only if print_info/=0)  : 2nd perturbation
!!  jband : band index for v1j and v2j
!!  print_info : if /=0 : all < vi | v1j > and < vi | v2j > are printed in std_out
!!  vi,v1j,v2j : input vectors
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
    case(1) ! < u^(pert2) | H^(0) | u^(pert1) >
      write(bra_i,'(2a,2(i1,a))') ' < du/',pert2,idir2,'(',iband,') | '
      write(ket_j,'(2a,2(i1,a))') ' | du/',pert1,idir1,'(',jband,') > '
      write(op1,'(a)')            '     H^(0)     '
      write(op2,'(a)')            '     S^(0)     '
    case(2) ! < u^(pert2) | H^(pert1) | u^(0) >
      write(bra_i,'(2a,2(i1,a))') ' < du/',pert2,idir2,'(',iband,') | '
      write(ket_j,'(a,i1,a)')     ' | u^(0) (',jband,') > '
      write(op1,'(2a,i1,a)')      '     dH/',pert1,idir1,'    '
      write(op2,'(2a,i1,a)')      '     dS/',pert1,idir1,'    '
    case(3) ! < u^(0) | H^(pert1pert2) | u^(0) >
      write(bra_i,'(a,i1,a)')     ' < u^(0) (',iband,') | '
      write(ket_j,'(a,i1,a)')     ' | u^(0) (',jband,') > '
      write(op1,'(2(2a,i1),a)')   'd^2H/(',pert1,idir1,' ',pert2,idir2,')'
      write(op2,'(2(2a,i1),a)')   'd^2S/(',pert1,idir1,' ',pert2,idir2,')'
    case(4) ! < u^(0) | H^(pert2) | u^(pert1) >
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
!!  rf2 : rf2_t object containing all rf2 data
!!  cg_jband (used only if print_info/=0) : array containing |u^(0)(jband)> (see NOTES below)
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk> (see NOTES below)
!!  copt=0 if work1 is | u^(0) >, cpopt=1 if work1 is | u^(1) >
!!  dtfil <type(datafiles_type)>=variables related to files
!!  eig0 : 0-order eigenvalue for the present wavefunction at k
!!  eig1_k_jband : first-order lagrange multipliers for the band j (Lambda^(1)_ji for all i)
!!  jband : band index of the input wavefunction
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  gvnl1(2,npw1*nspinor)=  part of <G|K^(1)+Vnl^(1)|C> not depending on VHxc^(1)              (sij_opt/=-1)
!!                       or part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1) (sij_opt==-1)
!!  ibg=shift to be applied on the location of data in the array cprj
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation of the Hamiltonian :
!!     ipert = 0           : GS calculation, call of getghc
!!     ipert = natom+1,2   : 1st order calculation, call of getgh1c
!!     ipert = natom+10,11 : 2nd order calculation, call of getgh2c
!!  ikpt=number of the k-point
!!  isppol=1 index of current spin component
!!  mband=maximum number of bands
!!  mkmem =number of k points trated by this node (GS data).
!!  mpi_enreg=information about MPI parallelization
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  print_info : if /=0 : some tests are done (see NOTES below). Wrong results are printed in std_out
!!  prtvol=control print volume and debugging output (for getghc)
!!  rf_hamk_idir <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q (here q=0)
!!  work1 : input wave function
!!
!! OUTPUT
!!  work2 : array containing H^(ipert)|work1>
!!  work3 : array containing S^(ipert)|work1>
!!
!! NOTES
!!  * Tests are done if print_info/=0. In that case : cg_jband(:,jband,1) = |u^(0)(jband)>
!!    For ipert==natom+2 (electric field) : cg_jband(:,jband,2) = i * |du/dk_idir(jband)>
!!    According to ipert, we check that : 
!!    -- ipert = 0 :
!!      < u^(0)(jband)| H^(0) | u^(0)(jband)> = eig0(jband)
!!    -- ipert = natom+1 or 2 : 
!!      < u^(0)(jband)| ( H^(1) - (eig0(jband)+eig0(iband))/2 * S^(1) ) | u^(0)(iband)> = Lambda^(1)(jband,iband)
!!    --ipert = natom+10 or 11 :
!!      nothing is done...
!!  * Use of cprj :
!!    The cprj array is computed in dfpt_looppert.F90. In the context of rf2 calculation, it contains
!!    <Proj_i^(0)|u^(0)> and <Proj_i^(1)|u^(0)> (dir=1,2 and 3) for all GS wavefunctions u^(0).
!!    If |work1> is a GS wavefunction (copt=0), computed cwaveprj can be used.
!!    If |work1> is not a GS wavefunction (copt=1), cwaveprj must be computed on the fly.
!!    Note that <Proj_i^(2)|u^(0)> is always computed on the fly.
!!
!! PARENTS
!!      rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_apply_hamiltonian(rf2,cg_jband,cprj,copt,dtfil,eig0,eig1_k_jband,jband,gs_hamkq,&
                                 gvnl1,ibg,idir,ipert,ikpt,isppol,mband,mkmem,mpi_enreg,nsppol,&
                                 print_info,prtvol,rf_hamk_idir,work1,work2,work3)

 use defs_basis
 use defs_abitypes
 use m_hamiltonian
 use m_cgtools

 use m_pawcprj,      only : pawcprj_type,pawcprj_alloc,pawcprj_get,pawcprj_free

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
 integer,intent(in) :: copt,ibg,idir,ipert,ikpt,isppol,jband,mband,mkmem,nsppol,print_info,prtvol
 type(datafiles_type),intent(in) :: dtfil
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf2_t),intent(in) :: rf2
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamk_idir
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in) :: cg_jband(2,rf2%size_wf*print_info*rf2%nband_k,2)
 real(dp),intent(in) :: eig0(rf2%nband_k),eig1_k_jband(2*rf2%nband_k)
 real(dp),intent(inout) :: gvnl1(2,rf2%size_wf)
 real(dp),intent(inout) :: work1(2,rf2%size_wf),work2(2,rf2%size_wf),work3(2,rf2%size_wf)
 type(pawcprj_type),intent(in) :: cprj(gs_hamkq%natom,gs_hamkq%nspinor*mband*mkmem*nsppol*gs_hamkq%usecprj)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: berryopt=0,iorder_cprj=0,optnl=2,tim_getghc=1,tim_getgh1c=1,tim_getgh2c=1 ! to change
 integer :: cpopt,iband,natom,nband_k,ncpgr_loc,sij_opt,size_wf,optlocal,opt_gvnl1,usevnl
 real(dp) :: dotr,doti,dotr2,doti2,tol_test
 character(len=500) :: msg
 
!arrays
 real(dp) :: dum1(1,1)
 real(dp),allocatable :: gvnlc(:,:),cgj(:,:),iddk(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 usevnl     = 0
 optlocal   = 0
 opt_gvnl1  = 0
 if(ipert==gs_hamkq%natom+2.or.(ipert==gs_hamkq%natom+11.and.gs_hamkq%usepaw==1)) then
   usevnl = 1
   opt_gvnl1 = 1
   if (ipert==gs_hamkq%natom+2) then
     optlocal = 1
     if (gs_hamkq%usepaw==1) opt_gvnl1 = 2
   end if
 end if
 dum1 = zero
 
 natom = gs_hamkq%natom
 nband_k = rf2%nband_k
 size_wf = rf2%size_wf
 sij_opt=1;if (gs_hamkq%usepaw==0) sij_opt=0
 tol_test=tol8

! In the PAW case : get cwaveprj
 ncpgr_loc = 0
 ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   ABI_DATATYPE_DEALLOCATE(cwaveprj)
   ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,gs_hamkq%nspinor))
   if (ipert==natom+1.or.ipert==natom+2) ncpgr_loc = 1
   if (ipert==natom+10.or.ipert==natom+11) ncpgr_loc = 3
   call pawcprj_alloc(cwaveprj,ncpgr_loc,gs_hamkq%dimcprj)
   if ((copt==0.or.print_info/=0).and.(ipert==natom+1.or.ipert==natom+2)) then ! get only the gradient along idir
     call pawcprj_get(gs_hamkq%atindx1,cwaveprj,cprj,natom,jband,ibg,ikpt,iorder_cprj,&
&    isppol,mband,mkmem,natom,1,nband_k,gs_hamkq%nspinor,nsppol,dtfil%unpaw,&
&    mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,ncpgr=3,icpgr=idir)
   else if ( (copt==0.and.(ipert==0.or.ipert==natom+10.or.ipert==natom+11)).or.&
             (print_info/=0.and.ipert==0) ) then
     call pawcprj_get(gs_hamkq%atindx1,cwaveprj,cprj,natom,jband,ibg,ikpt,iorder_cprj,&
&    isppol,mband,mkmem,natom,1,nband_k,gs_hamkq%nspinor,nsppol,dtfil%unpaw,&
&    mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
   end if
 end if

! *******************************************************************************************
! apply H^(0)
! *******************************************************************************************
 if (ipert == 0) then

   ABI_ALLOCATE(gvnlc,(2,size_wf))
   gvnlc(:,:) = zero
!  Test if < u^(0) | H^(0) | u^(0) > = eig0(jband)
   if(print_info/=0) then
     ABI_ALLOCATE(cgj,(2,size_wf))
     cgj(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     cpopt = -1+3*gs_hamkq%usecprj*gs_hamkq%usepaw
     call getghc(cpopt,cgj,cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,mpi_enreg,1,prtvol,&
                 sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgj,work2,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     dotr = dotr - eig0(jband)
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol_test) then
       write(msg,'(a,es22.13E3)') 'RF2 TEST GETGHC : NOT PASSED dotr = ',dotr
       call wrtout(std_out,msg)
     end if
     ABI_DEALLOCATE(cgj)

   end if ! end tests

   if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!    copt=0 : |work1> is a GS wavefunction, cwaveprj are already in memory
!    copt=1 : |work1> is not a GS wavefunction, cwaveprj are computed here
     cpopt=2-3*copt
   else
     cpopt = -1
   end if
   call getghc(cpopt,work1,cwaveprj,work2,work3,gs_hamkq,gvnlc,zero,mpi_enreg,1,prtvol,&
               sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   ABI_DEALLOCATE(gvnlc)

! *******************************************************************************************
! apply H^(1)
! *******************************************************************************************
 else if (ipert==natom+1.or.ipert==natom+2) then

!  Test if < u^(0) | ( H^(1) - eps^(0) S^(1) ) | u^(0) > = eig^(1)
   if(print_info/=0) then
     ABI_ALLOCATE(cgj,(2,size_wf))
     ABI_ALLOCATE(iddk,(2,size_wf))
     cgj(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     iddk(:,:) = zero
     if (ipert==natom+2) iddk(:,:) = cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,2)
     call getgh1c(berryopt,0,cgj,cwaveprj,work2,dum1,work3,gs_hamkq,iddk,idir,ipert,zero,&
                  mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)
     do iband=1,nband_k
       cgj(:,:) = cg_jband(:,1+(iband-1)*size_wf:iband*size_wf,1)
       call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgj,work2,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
       if (gs_hamkq%usepaw==1.and.ipert/=natom+2) then ! S^(1) is zero for ipert=natom+2
         call dotprod_g(dotr2,doti2,gs_hamkq%istwf_k,size_wf,2,cgj,work3,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
         dotr = dotr - (eig0(iband)+eig0(jband))*dotr2/two
         doti = doti - (eig0(iband)+eig0(jband))*doti2/two
       end if
       dotr = dotr - eig1_k_jband(1+2*(iband-1))
       doti = doti - eig1_k_jband(2+2*(iband-1))
       dotr = sqrt(dotr**2+doti**2)
       if (dotr > tol_test) then
         write(msg,'(4(a,i2),a,es22.13E3)') 'RF2 TEST GETGH1 : ipert=',ipert-natom,' idir=',idir,&
                                            ' jband=',jband,' iband=',iband,' NOT PASSED dotr = ',dotr
         call wrtout(std_out,msg)
       end if
     end do ! end iband
     ABI_DEALLOCATE(cgj)
     ABI_DEALLOCATE(iddk)
   end if ! end tests

   call getgh1c(berryopt,copt,work1,cwaveprj,work2,dum1,work3,gs_hamkq,gvnl1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl)

! *******************************************************************************************
! apply H^(2)
! *******************************************************************************************
 else if (ipert==natom+10.or.ipert==natom+11) then

     call getgh2c(work1,cwaveprj,work2,work3,gs_hamkq,gvnl1,idir,ipert,zero,&
                   mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamk_idir,sij_opt,tim_getgh2c,usevnl)

 end if

 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   call pawcprj_free(cwaveprj)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

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
!! rf2 : the rf2_t object to destroy
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
