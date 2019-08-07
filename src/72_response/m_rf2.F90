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
!!  Copyright (C) 2015-2019 ABINIT group (LB,MT)
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
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_hamiltonian
 use m_cgtools

 use m_getgh1c,     only : getgh1c
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_copy, pawcprj_free, pawcprj_output
 use m_getghc,      only : getghc
 use m_nonlop,      only : nonlop
 use m_getgh2c,     only : getgh2c

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
  integer :: size_wf ! number of components in a wavefunction
  integer :: size_cprj ! number of components of a cprj variable (i.e. <p_i|WF>)

!arrays
  integer :: iperts(2) ! perturbations to compute
   ! ipert = natom + 10 :
   !   iperts(1) = iperts(2) = natom+1 (wavevector k)
   !
   ! ipert = natom + 11 :
   !   iperts(1) = natom+1 (wavevector k)
   !   iperts(2) = natom+2 (electric field)

  integer :: idirs(2) ! directions of the perturbations (ndir=1 : idirs(1)=idirs(2) , ndir=2 : idirs(1)/=idirs(2))

  real(dp), ABI_CONTIGUOUS pointer :: RHS_Stern(:,:) => null() ! need a pointer because is a ptr target
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
   !  lambda_dir2 = k_dir2 (ipert==natom+10) , E_dir2 (ipert==natom+11)
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

 public :: rf2_getidir
 public :: rf2_getidirs
 public :: rf2_accumulate_bands
 public :: rf2_apply_hamiltonian
 public :: rf2_destroy

!!***

!----------------------------------------------------------------------

CONTAINS  !==============================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_rf2/rf2_getidir
!! NAME
!! rf2_getidir
!!
!! FUNCTION
!!  Get the direction of the 2nd order perturbation according to inputs idir1 and idir2
!!
!! INPUTS
!!  idir1 : index of the 1st direction (1<=idir1<=3)
!!  idir2 : index of the 2nd direction (1<=idir2<=3)

!! OUTPUT
!!  idir : integer between 1 and 9
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_looppert,dfpt_scfcv,rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_getidir(idir1,idir2,idir)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: idir1,idir2
 integer,intent(out) :: idir
 integer :: dir_from_dirs(9) = (/1,6,5,9,2,4,8,7,3/)

! *************************************************************************

 idir = dir_from_dirs(3*(idir1-1)+idir2)

end subroutine rf2_getidir
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
!!  idir1  (used only if debug_mode/=0)  : direction of the 1st perturbation
!!  idir2  (used only if debug_mode/=0)  : direction of the 2nd perturbation
!!  ipert1 (used only if debug_mode/=0)  : 1st perturbation
!!  ipert2 (used only if debug_mode/=0)  : 2nd perturbation
!!  jband : band index for v1j and v2j
!!  debug_mode : if /=0 : all < vi | v1j > and < vi | v2j > are printed in std_out
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
                                 jband,debug_mode,vi,v1j,v2j)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,iband,idir1,idir2,ipert1,ipert2,jband,debug_mode
 type(rf2_t),intent(inout) :: rf2
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in) :: vi(2,rf2%size_wf),v1j(2,rf2%size_wf),v2j(2,rf2%size_wf)

!Local variables ---------------------------------------
!scalars
 integer :: nband_k,size_wf
 real(dp) :: dotr,dot2r,doti,dot2i,factor
 character(len=500) :: msg
 character(len=15) :: bra_i,ket_j,op1,op2
 character(len=2) :: pert1,pert2

! *************************************************************************

 DBG_ENTER("COLL")

 nband_k = rf2%nband_k
 size_wf = rf2%size_wf
 factor = one
 if(rf2%ndir==1 .and. choice /= 3) factor = two ! in order to not compute same terms twice

 call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v1j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

 if(debug_mode/=0) then
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
   dot2r = dotr ; if (abs(dot2r)<tol9) dot2r = zero ! in order to hide the numerical noise
   dot2i = doti ; if (abs(dot2i)<tol9) dot2i = zero ! in order to hide the numerical noise
   write(msg,'(3a,2(a,es17.8E3))') bra_i,op1,ket_j,' = ',dot2r,',',dot2i
   call wrtout(std_out,msg)
 end if

 rf2%lambda_mn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + rf2%lambda_mn(:,iband+(jband-1)*nband_k)

 if (choice == 1 .or. gs_hamkq%usepaw==1) then
   call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,vi,v2j,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

   if(debug_mode/=0) then
     dot2r = dotr ; if (abs(dot2r)<tol9) dot2r = zero ! in order to hide the numerical noise
     dot2i = doti ; if (abs(dot2i)<tol9) dot2i = zero ! in order to hide the numerical noise
     write(msg,'(3a,2(a,es17.8E3))') bra_i,op2,ket_j,' = ',dot2r,',',dot2i
     call wrtout(std_out,msg)
   end if

   rf2%amn(:,iband+(jband-1)*nband_k) = factor*(/dotr,doti/) + rf2%amn(:,iband+(jband-1)*nband_k)

 end if ! end choice

 DBG_EXIT("COLL")

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
!!  cg_jband (used only if debug_mode/=0) : array containing |u^(0)(jband)> for all bands
!!  cprj_jband(natom,nspinor*usecprj)= u^(0) wave functions for all bands
!!              projected with non-local projectors: cprj_jband=<p_i|u^(0)(jband)>
!!  cwave(2,size_wf) : input wave function |u>
!!  cwaveprj(natom,nspinor) : input wave function |u>
!!              projected with non-local projectors <p_i|u>
!!  eig0 : 0-order eigenvalue for the present wavefunction at k
!!  [enl]=optional (if not present, use hamk%ekb); non-local coeffs connecting projectors
!!  eig1_k_jband : first-order lagrange multipliers for the band j (Lambda^(1)_ji for all i)
!!  [ffnl1]=nonlocal form factors (needed for tests of ipert>=natom+10)
!!  [ffnl1_test]=nonlocal form factors used for tests (needed for tests of ipert>=natom+10)
!!  jband : band index of the input wavefunction
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  gvnlx1(2,npw1*nspinor)=  part of <G|K^(1)+Vnl^(1)|C> not depending on VHxc^(1)              (sij_opt/=-1)
!!                       or part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1) (sij_opt==-1)
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation of the Hamiltonian :
!!     ipert = 0           : GS calculation, call of getghc
!!     ipert = natom+1,2   : 1st order calculation, call of getgh1c
!!     ipert = natom+10,11 : 2nd order calculation, call of getgh2c
!!  ikpt=number of the k-point
!!  isppol=1 index of current spin component
!!  mkmem =number of k points trated by this node (GS data).
!!  mpi_enreg=information about MPI parallelization
!!  nband_k=number of bands at this k point for that spin polarization
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  debug_mode : if /=0 : some tests are done (see NOTES below). Wrong results are printed in std_out
!!  prtvol=control print volume and debugging output (for getghc)
!!  rf_hamk_idir <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q (here q=0)
!!  size_cprj=size of a cprj array (=gs_hamkq%nspinor)
!!  size_wf=size of a wavefunction array (=gs_hamkq%npw_k*gs_hamkq%nspinor)
!!
!! OUTPUT
!!  h_cwave(2,size_wf) : array containing H^(ipert)|cwave>
!!  s_cwave(2,size_wf) : array containing S^(ipert)|cwave>
!!
!! NOTES
!!  * Tests are done if debug_mode/=0. In that case : cg_jband(:,jband,1) = |u^(0)(jband)>
!!    For ipert==natom+2 (electric field) : cg_jband(:,jband,2) = i * |du/dk_idir(jband)>
!!    According to ipert, we check that :
!!    -- ipert = 0 :
!!      < u^(0)(jband)| H^(0) | u^(0)(jband)> = eig0(jband)
!!    -- ipert = natom+1 or 2 :
!!      < u^(0)(jband)| ( H^(1) - (eig0(jband)+eig0(iband))/2 * S^(1) ) | u^(0)(iband)> = Lambda^(1)(jband,iband)
!!    --ipert >= natom+10 :
!!      < u^(0) | H^(2) | u^(0) > from getgh2c is equal to nonlop with signs=1
!!  * Use of cprj :
!!    The cprj array is computed in dfpt_looppert.F90. In the context of rf2 calculation, it contains
!!    <Proj_i^(0)|u^(0)> and <Proj_i^(1)|u^(0)> (dir=1,2 and 3) for all wavefunctions u^(0) or
!!    <Proj_i^(0)|u^(1)> and <Proj_i^(1)|u^(1)> (dir=1,2 and 3) for all 1st-order wavefunctions u^(1).
!!    Note that <Proj_i^(2)|u^(0)> is always computed on the fly.
!!
!! PARENTS
!!      rf2_init
!!
!! CHILDREN
!!
!! SOURCE

subroutine rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave,cwaveprj,h_cwave,s_cwave,eig0,eig1_k_jband,&
&                                jband,gs_hamkq,gvnlx1,idir,ipert,ikpt,isppol,mkmem,mpi_enreg,nband_k,nsppol,&
&                                debug_mode,prtvol,rf_hamk_idir,size_cprj,size_wf,&
&                                conj,enl,ffnl1,ffnl1_test) ! optional

!Arguments ---------------------------------------------
!scalars
 logical,intent(in),optional :: conj
 integer,intent(in) :: idir,ipert,ikpt,isppol,jband,mkmem,nband_k,nsppol,debug_mode,prtvol,size_wf,size_cprj
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
! type(rf2_t),intent(in) :: rf2
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamk_idir
 type(MPI_type),intent(in) :: mpi_enreg

!arrays
 real(dp),intent(in),target :: cg_jband(2,size_wf*debug_mode*nband_k,2)
 real(dp),intent(in),optional,target :: enl(gs_hamkq%dimekb1,gs_hamkq%dimekb2,gs_hamkq%nspinor**2,gs_hamkq%dimekbq)
 real(dp),intent(in),optional :: ffnl1(:,:,:,:),ffnl1_test(:,:,:,:)
 real(dp),intent(in) :: eig0(nband_k),eig1_k_jband(2*nband_k)
 real(dp),intent(inout) :: gvnlx1(2,size_wf)
 real(dp),intent(inout) :: cwave(2,size_wf),h_cwave(2,size_wf),s_cwave(2,size_wf)
 type(pawcprj_type),intent(inout),target :: cprj_jband(:,:),cwaveprj(:,:)

!Local variables ---------------------------------------
!scalars
 integer,parameter :: berryopt=0,tim_getghc=0,tim_getgh1c=0,tim_getgh2c=0,tim_nonlop=0
 integer,parameter :: alpha(9)=(/1,2,3,2,1,1,3,3,2/),beta(9)=(/1,2,3,3,3,2,2,1,1/)
 integer :: choice,cpopt,iatom,iband,idirc,idir1,idir2,idir_dum,natom,nnlout,paw_opt,signs,sij_opt
 integer :: opt_gvnlx1,opt_gvnl2,optlocal,optnl,usevnl
 logical :: compute_conjugate,has_cprj_jband,has_cwaveprj,pert_phon_elfd
 real(dp) :: dotr,doti,dotr2,doti2,enlout(18*gs_hamkq%natom),tol_test
 real(dp), pointer :: enl_ptr(:,:,:,:)
 real(dp),allocatable,target :: enl_temp(:,:,:,:)
 character(len=500) :: msg

!arrays
 real(dp),target :: cwave_empty(0,0),svectout_dum(0,0)
 real(dp),allocatable :: gvnlxc(:,:),iddk(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: cwave_i(:,:),cwave_j(:,:)
 type(pawcprj_type),target :: cprj_empty(0,0)
 type(pawcprj_type),pointer :: cprj_j(:,:)
! *********************************************************************

 DBG_ENTER("COLL")

 compute_conjugate = .false.
 if(present(conj)) compute_conjugate = conj

 ABI_UNUSED(ikpt)
 ABI_UNUSED(isppol)
 ABI_UNUSED(mkmem)
 ABI_UNUSED(nsppol)

!Check sizes
 if (size(cprj_jband)/=0) then
   if (size(cprj_jband)/=nband_k*gs_hamkq%natom*gs_hamkq%nspinor*gs_hamkq%usecprj) then
     write(msg,'(2(a,i10))') &
     'Wrong cprj size! actual size = ',size(cprj_jband),&
     ' good size = ',nband_k*gs_hamkq%natom*gs_hamkq%nspinor*gs_hamkq%usecprj
     MSG_BUG(msg)
   end if
 end if
 if (size(cwaveprj)/=0) then
   if (size(cwaveprj)/=gs_hamkq%natom*gs_hamkq%nspinor*gs_hamkq%usecprj) then
     write(msg,'(2(a,i10))') &
     'Wrong cwaveprj size! actual size = ',size(cwaveprj),&
     ' good size = ',gs_hamkq%natom*gs_hamkq%nspinor*gs_hamkq%usecprj
     MSG_BUG(msg)
   end if
 end if

 natom = gs_hamkq%natom

 pert_phon_elfd = .false.
 if (ipert>natom+11.and.ipert<=2*natom+11) pert_phon_elfd = .true.
 usevnl     = 0
 opt_gvnlx1  = 0
 opt_gvnl2  = 0
 optnl      = 2
 sij_opt=1;if (gs_hamkq%usepaw==0) sij_opt=0
 if (ipert==natom+2.or.(ipert==natom+11.and.gs_hamkq%usepaw==1).or.pert_phon_elfd) then
   usevnl = 1
   opt_gvnlx1 = 2
   opt_gvnl2 = 1
 end if
 optlocal = 1
 if (pert_phon_elfd) then
   optlocal = 0
   sij_opt = 0
   optnl = 1
 end if
 tol_test=tol8

! In the PAW case : manage cprj_jband, cwaveprj
 has_cprj_jband=(gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1.and.size(cprj_jband)/=0)
 has_cwaveprj=(gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1.and.size(cwaveprj)/=0)
 cprj_j => cprj_empty

 if (debug_mode/=0) then
   if (ipert/=0) then
     write(msg,'(5(a,i4))') 'RF2 TEST rf2_apply_hamiltonian ipert = ',ipert,' idir =',idir,' isppol =',isppol,&
     & ' ikpt = ',ikpt,' jband = ',jband
   else
     write(msg,'(3(a,i4))') 'RF2 TEST rf2_apply_hamiltonian ipert =    0 idir =    0 isppol =',isppol,&
     & ' ikpt = ',ikpt,' jband = ',jband
   end if
   call wrtout(std_out,msg)
 end if

! *******************************************************************************************
! apply H^(0)
! *******************************************************************************************
 if (ipert == 0) then

   ABI_ALLOCATE(gvnlxc,(2,size_wf))
   gvnlxc(:,:) = zero

!  Test if < u^(0) | H^(0) | u^(0) > = eig0(jband)
   if(debug_mode/=0) then
     cwave_j => cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     if (has_cprj_jband) cprj_j => cprj_jband(:,1+(jband-1)*size_cprj:jband*size_cprj)
     cpopt = -1+3*gs_hamkq%usecprj*gs_hamkq%usepaw
     call getghc(cpopt,cwave_j,cprj_j,h_cwave,s_cwave,gs_hamkq,gvnlxc,zero,mpi_enreg,1,prtvol,&
                 sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_j,h_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     dotr = dotr - eig0(jband)
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > tol_test) then
       write(msg,'(a,es17.8E3)') 'RF2 TEST GETGHC : NOT PASSED dotr = ',dotr
       call wrtout(ab_out,msg)
       call wrtout(std_out,msg)
     end if
   end if ! end tests

   cpopt=-1;if (has_cwaveprj) cpopt=2
   call getghc(cpopt,cwave,cwaveprj,h_cwave,s_cwave,gs_hamkq,gvnlxc,zero,mpi_enreg,1,prtvol,&
               sij_opt,tim_getghc,0,select_k=KPRIME_H_KPRIME)
   ABI_DEALLOCATE(gvnlxc)

! *******************************************************************************************
! apply H^(1)
! *******************************************************************************************
 else if (ipert<=natom+2) then

!  Test if < u^(0) | ( H^(1) - eps^(0) S^(1) ) | u^(0) > = eig^(1)
   if(debug_mode/=0) then
     ABI_ALLOCATE(iddk,(2,size_wf))
     cwave_j => cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)
     if (has_cprj_jband) cprj_j => cprj_jband(:,1+(jband-1)*size_cprj:jband*size_cprj)
     iddk(:,:) = zero;if (ipert==natom+2) iddk(:,:)=cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,2)
     call getgh1c(berryopt,cwave_j,cprj_j,h_cwave,cwave_empty,s_cwave,gs_hamkq,iddk,idir,ipert,zero,&
                  mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl,conj=compute_conjugate)
     do iband=1,nband_k
       cwave_i => cg_jband(:,1+(iband-1)*size_wf:iband*size_wf,1)
       call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_i,h_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       if (gs_hamkq%usepaw==1.and.ipert/=natom+2) then ! S^(1) is zero for ipert=natom+2
         call dotprod_g(dotr2,doti2,gs_hamkq%istwf_k,size_wf,2,cwave_i,s_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         dotr = dotr - (eig0(iband)+eig0(jband))*dotr2/two
         doti = doti - (eig0(iband)+eig0(jband))*doti2/two
       end if
       dotr = dotr - eig1_k_jband(1+2*(iband-1))
       doti = doti - eig1_k_jband(2+2*(iband-1))
       dotr = sqrt(dotr**2+doti**2)
       if (dotr > tol_test) then
         write(msg,'(4(a,i2),a,es17.8E3)') 'RF2 TEST GETGH1 : ipert=',ipert,' idir=',idir,&
                                            ' jband=',jband,' iband=',iband,' NOT PASSED dotr = ',dotr
         call wrtout(ab_out,msg)
         call wrtout(std_out,msg)
       end if
     end do ! end iband
     ABI_DEALLOCATE(iddk)
   end if ! end tests

   call getgh1c(berryopt,cwave,cwaveprj,h_cwave,cwave_empty,s_cwave,gs_hamkq,gvnlx1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamk_idir,sij_opt,tim_getgh1c,usevnl,conj=compute_conjugate)

! *******************************************************************************************
! apply H^(2)
! *******************************************************************************************
 else if (ipert==natom+10.or.ipert==natom+11.or.pert_phon_elfd) then

! *******************************************************************************************
!  Test if < u^(0) | H^(2) | u^(0) > from getgh2c is equal to nonlop with signs=1
   if(debug_mode/=0.and.present(ffnl1).and.present(ffnl1_test)) then

     cwave_j => cg_jband(:,1+(jband-1)*size_wf:jband*size_wf,1)

     ABI_ALLOCATE(iddk,(2,size_wf))
     iddk(:,:) = cwave_j(:,:)

     call getgh2c(cwave_j,cprj_empty,h_cwave,s_cwave,gs_hamkq,iddk,idir,ipert,zero,&
                  mpi_enreg,optlocal,optnl,opt_gvnl2,rf_hamk_idir,sij_opt,tim_getgh2c,usevnl,&
                  conj=compute_conjugate,optkin=0,enl=enl)
     ABI_DEALLOCATE(iddk)

     call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_j,h_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

     idir1=alpha(idir);idir2=beta(idir)
     idirc=3*(idir1-1)+idir2

     enlout = zero

!    Change the pointer ffnl_k to ffnl1_test (idir_ffnl=0, for nonlop with signs=1)
     call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnl1_test)

     signs=1
     dotr2 = 0
     doti2 = 0

!    ************************
!    IPERT == NATOM+10
!    ************************
     if (ipert==natom+10) then

       cpopt=-1; choice=8; paw_opt=1; if (gs_hamkq%usepaw==0) paw_opt=0 ; nnlout = 6

       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idirc,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum)
       dotr2 = enlout(idir)

     end if ! IPERT == NATOM+10

!    ************************
!    IPERT == NATOM+11
!    ************************
     if (ipert==natom+11) then

       if (present(enl)) then
         enl_ptr => enl
       else if (associated(rf_hamk_idir%e1kbfr).and.associated(rf_hamk_idir%e1kbsc).and.optnl==2) then
         ABI_ALLOCATE(enl_temp,(gs_hamkq%dimekb1,gs_hamkq%dimekb2,gs_hamkq%nspinor**2,rf_hamk_idir%cplex))
         enl_temp(:,:,:,:) = rf_hamk_idir%e1kbfr(:,:,:,:) + rf_hamk_idir%e1kbsc(:,:,:,:)
         enl_ptr => enl_temp
       else if (associated(rf_hamk_idir%e1kbfr)) then
         enl_ptr => rf_hamk_idir%e1kbfr
       end if

!      Compute application of dS/dk1 to cwavef
!      sum_{i,j} s_ij  < psi^(0) | d(|p_i><p_j|)/dk(idir1) | psi^(0) >
       cpopt=-1 ; choice=5 ; paw_opt=3 ; nnlout=3
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum)
       dotr2 = enlout(idir1)

!      Compute part of H^(2) due to derivative of projectors (idir1) and derivative of Dij (idir2)
!      sum_{i,j} chi_ij(idir2) < psi^(0) | d(|p_i><p_j|)/dk(idir1) | psi^(0) >
       cpopt=-1 ; choice=5 ; paw_opt=1 ; nnlout=3
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum,enl=enl_ptr)
       dotr2 = dotr2 + enlout(idir1)

!      Compute derivatives due to projectors |d^2[p_i]/dk1dk2>,|d[p_i]/dk1>,|d[p_i]/dk2>
!      i * sum_{i,j} < psi^(0) | (d(|p_i><dp_j/dk(idir2)|)/dk(idir1) | psi^(0) >
       cpopt=-1 ; choice=81 ; paw_opt=3 ; nnlout=18
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum)
       dotr2 = dotr2 - enlout(2*idirc  )
       doti2 = doti2 + enlout(2*idirc-1)

     end if ! IPERT == NATOM+11

!    ************************
!    PERT_PHON_ELFD
!    ************************
     if (pert_phon_elfd) then

       if (present(enl)) then
         enl_ptr => enl
       else if (associated(rf_hamk_idir%e1kbfr).and.associated(rf_hamk_idir%e1kbsc).and.optnl==2) then
         ABI_ALLOCATE(enl_temp,(gs_hamkq%dimekb1,gs_hamkq%dimekb2,gs_hamkq%nspinor**2,rf_hamk_idir%cplex))
         enl_temp(:,:,:,:) = rf_hamk_idir%e1kbfr(:,:,:,:) + rf_hamk_idir%e1kbsc(:,:,:,:)
         enl_ptr => enl_temp
       else if (associated(rf_hamk_idir%e1kbfr)) then
         enl_ptr => rf_hamk_idir%e1kbfr
       end if

       iatom = gs_hamkq%atindx(ipert - natom - 11)  ! Atoms are type-sorted in enlout()

!      Compute application of dS/dtau1 to i*d[cwavef]/dk2
!      sum_{i,j} s_ij < psi^(0) | d(|p_i><p_j|)/dtau(idir1) | i*psi^(k(idir2)) >
       cpopt=-1 ; choice=2 ; paw_opt=3 ; nnlout=3*natom
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum)
       dotr2 = enlout(3*(iatom-1)+idir1)

!      Compute part of H^(2) due to derivative of projectors (idir1) and derivative of Dij (idir2)
!      sum_{i,j} chi_ij(idir2) < psi^(0) | d(|p_i><p_j|)/dtau(idir1) | psi^(0) >
       cpopt=-1 ; choice=2 ; paw_opt=1 ; nnlout=3*natom
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum,enl=enl_ptr)
       dotr2 = dotr2 + enlout(3*(iatom-1)+idir1)

!      Compute derivatives due to projectors |d^2[p_i]/dtau1dk2>,|d[p_i]/dtau1>,|d[p_i]/dk2>
!      i * sum_{i,j} < psi^(0) | (d(|p_i><dp_j/dk(idir2)|)/dtau(idir1) | psi^(0) >
       cpopt=-1 ; choice=54 ; paw_opt=3 ; nnlout=18*natom
       call nonlop(choice,cpopt,cprj_empty,enlout,gs_hamkq,idir_dum,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave_j,svectout_dum)
       dotr2 = dotr2 - enlout(18*(iatom-1)+2*idirc  )
       doti2 = doti2 + enlout(18*(iatom-1)+2*idirc-1)

     end if ! PERT_PHON_ELFD

     dotr = dotr - dotr2
     doti = doti - doti2
     dotr = sqrt(dotr**2+doti**2)
     if (dotr > 10*tol_test) then
       write(msg,'(a,es17.8E3)') 'RF2 TEST GETGH2C : NOT PASSED dotr = ',dotr
       call wrtout(ab_out,msg)
       call wrtout(std_out,msg)
     end if

!    Change the pointer ffnl_k back to ffnl1 (idir_ffnl=4, for nonlop with signs=2)
     call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnl1)

     if (associated(enl_ptr)) then
       nullify(enl_ptr)
     end if
     if (allocated(enl_temp)) then
       ABI_DEALLOCATE(enl_temp)
     end if

   end if ! end tests
! *******************************************************************************************

   call getgh2c(cwave,cwaveprj,h_cwave,s_cwave,gs_hamkq,gvnlx1,idir,ipert,zero,&
                mpi_enreg,optlocal,optnl,opt_gvnl2,rf_hamk_idir,sij_opt,tim_getgh2c,usevnl,&
                conj=compute_conjugate,enl=enl)

 else

   h_cwave = zero
   MSG_ERROR(" rf2_apply_hamiltonian can be used only for : 0<=ipert<=natom+2 and natom+10<=ipert<=2*natom+11.")
   return

 end if

 DBG_EXIT("COLL")

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

!Arguments ---------------------------------------------
!scalars
 type(rf2_t),intent(inout) :: rf2
!arrays

!Local variables ---------------------------------------
!scalars

! *************************************************************************

 if (associated(rf2%RHS_Stern)) then
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
