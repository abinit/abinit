!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_getgh2c
!! NAME
!!  m_getgh2c
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2019 ABINIT group (MT,JLJ)
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

module m_getgh2c

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors

 use m_pawcprj,     only : pawcprj_type,pawcprj_alloc,pawcprj_free
 use m_hamiltonian, only : gs_hamiltonian_type,rf_hamiltonian_type
 use m_nonlop,      only : nonlop

 implicit none

 private
!!***

 public :: getgh2c
!!***

contains
!!***

!!****f* ABINIT/getgh2c
!! NAME
!! getgh2c
!!
!! FUNCTION
!! Compute <G|H^(2)|C> (or <G|H^(2)-Eps.S^(2)|C>) for input vector |C> expressed in reciprocal space.
!! (H^(2) is the 2nd-order pertubed Hamiltonian, S^(2) is the 2nd-order perturbed overlap operator).
!! Result is put in array gh2c.
!! If required, part of <G|K(2)+Vnonlocal^(2)|C> not depending on VHxc^(2) is also returned in gvnl2.
!! If required, <G|S^(2)|C> is returned in gs2c (S=overlap - PAW only)
!! Available for the following cases :
!!  ipert = natom+10 (dkdk)   :       2nd derivative w.r.t wavevector
!!          natom+11 (dkdE)   : mixed 2nd derivative w.r.t wavector     and eletric field
!!  also if natom+12<=ipert<=2*natom+11 :
!!                   (dtaudE) : mixed 2nd derivative w.r.t atom. displ. and eletric field (nonlocal only)
!!
!! INPUTS
!!  cwavef(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C>
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  lambda=real use to apply H^(2)-lambda.S^(2)
!!  mpi_enreg=information about MPI parallelization
!!  optlocal=0: local part of H^(2) is not computed in gh2c=<G|H^(2)|C>
!!           1: local part of H^(2) is computed in gh2c=<G|H^(2)|C>
!!  optnl=0: non-local part of H^(2) is not computed in gh2c=<G|H^(2)|C>
!!        1: non-local part of H^(2) depending on VHxc^(2) is not computed in gh2c=<G|H^(2)|C>
!!        2: non-local part of H^(2) is totally computed in gh2c=<G|H^(2)|C>
!!  opt_gvnl2=option controlling the use of gvnl2 array:
!!            0: not used
!!            1: used as input:    - used only for PAW and ipert=natom+11/+12
!!               At input: contains the derivative w.r.t wavevector of cwavef (times i)
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 2nd-order Hamiltonian at k,k+q
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H^(2)|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S^(2)|C> have to be computed in gs2c in addition to gh2c
!!                       if -1, matrix elements <G|H^(2)-lambda.S^(2)|C> have to be computed in gh2c (gs2c not used)
!!  tim_getgh2c=timing code of the calling subroutine (can be set to 0 if not attributed)
!!  usevnl=1 if gvnl2=(part of <G|K^(2)+Vnl^(2)-lambda.S^(2)|C> not depending on VHxc^(2)) has to be input/output
!!
!! OUTPUT
!! gh2c(2,npw1*nspinor)= <G|H^(2)|C> or  <G|H^(2)-lambda.S^(2)|C>
!!                     (only kinetic+non-local parts if optlocal=0)
!! if (usevnl==1)
!!  gvnl2(2,npw1*nspinor)=  part of <G|K^(2)+Vnl^(2)|C> not depending on VHxc^(2)              (sij_opt/=-1)
!!                       or part of <G|K^(2)+Vnl^(2)-lambda.S^(2)|C> not depending on VHxc^(2) (sij_opt==-1)
!! if (sij_opt=1)
!!  gs2c(2,npw1*nspinor)=<G|S^(2)|C> (S=overlap).
!!
!! PARENTS
!!      m_rf2
!!
!! CHILDREN
!!      nonlop,pawcprj_alloc,pawcprj_free
!!
!! SOURCE

subroutine getgh2c(cwavef,cwaveprj,gh2c,gs2c,gs_hamkq,gvnl2,idir,ipert,lambda,&
&                  mpi_enreg,optlocal,optnl,opt_gvnl2,rf_hamkq,sij_opt,tim_getgh2c,usevnl,conj,enl,optkin)

!Arguments ------------------------------------
!scalars
 logical,intent(in),optional :: conj
 integer,intent(in) :: idir,ipert,optlocal,optnl,opt_gvnl2,sij_opt,tim_getgh2c,usevnl
 integer,intent(in),optional :: optkin
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
!arrays
 real(dp),intent(in),optional,target :: enl(gs_hamkq%dimekb1,gs_hamkq%dimekb2,gs_hamkq%nspinor**2,gs_hamkq%dimekbq)
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(inout),target :: gvnl2(:,:)
 real(dp),intent(out) :: gh2c(:,:),gs2c(:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_nonlop=0
 integer :: choice,cpopt,iatm,idir1,idir2,idirc,ipw,ipws,ispinor,my_nspinor
 integer :: natom,ncpgr,nnlout=1,npw,npw1,paw_opt,signs,usecprj
 logical :: compute_conjugate,has_kin,has_vnl,pert_phon_elfd
 real(dp) :: enlout_dum(1)
 character(len=500) :: msg
!arrays
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(9)=(/1,2,3,2,1,1,3,3,2/),beta(9)=(/1,2,3,3,3,2,2,1,1/)
!real(dp) :: tsec(2)
 real(dp) :: svectout_dum(1,1),vectout_dum(1,1)
 real(dp),allocatable :: nonlop_out(:,:)
 real(dp), pointer :: gvnl2_(:,:)
 real(dp), pointer :: ddkinpw(:),kinpw1(:),enl_ptr(:,:,:,:)
 real(dp),allocatable,target :: enl_temp(:,:,:,:)
 type(pawcprj_type),allocatable,target :: cwaveprj_tmp(:,:)
 type(pawcprj_type),pointer :: cwaveprj_ptr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")
 ABI_UNUSED(tim_getgh2c)

!Keep track of total time spent in getgh2c
!call timab(196+tim_getgh2c,1,tsec)

!======================================================================
!== Initialisations and compatibility tests
!======================================================================

 npw  =gs_hamkq%npw_k
 npw1 =gs_hamkq%npw_kp
 natom=gs_hamkq%natom

!Compatibility tests
 if(ipert/=natom+10.and.ipert/=natom+11.and.ipert>2*natom+11)then
   msg='only ipert==natom+10/+11 and natom+11<=ipert<=2*natom+11 implemented!'
   MSG_BUG(msg)
 end if
 pert_phon_elfd = .false.
 if (ipert>natom+11.and.ipert<=2*natom+11) pert_phon_elfd = .true.
 if (mpi_enreg%paral_spinor==1) then
   msg='Not compatible with parallelization over spinorial components!'
   MSG_BUG(msg)
 end if
 if (gs_hamkq%nvloc>1) then
   msg='Not compatible with nvloc=4 (non-coll. magnetism)!'
   MSG_BUG(msg)
 end if
 if((ipert==natom+11.or.pert_phon_elfd).and.gs_hamkq%usepaw==1.and.optnl>=1) then
   if (gs_hamkq%nvloc>1) then
     msg='Not compatible with nvloc=4 (non-coll. magnetism)!'
     MSG_BUG(msg)
   end if
   if (present(enl)) then
     enl_ptr => enl
   else if (associated(rf_hamkq%e1kbfr).and.associated(rf_hamkq%e1kbsc).and.optnl==2) then
     ABI_CHECK(size(rf_hamkq%e1kbfr,4)==1,'BUG in getgh2c: qphase>1!')
     ABI_CHECK(size(rf_hamkq%e1kbsc,4)==1,'BUG in getgh2c: qphase>1!')
     ABI_ALLOCATE(enl_temp,(gs_hamkq%dimekb1,gs_hamkq%dimekb2,gs_hamkq%nspinor**2,gs_hamkq%dimekbq))
     enl_temp(:,:,:,:) = rf_hamkq%e1kbfr(:,:,:,:) + rf_hamkq%e1kbsc(:,:,:,:)
     enl_ptr => enl_temp
   else if (associated(rf_hamkq%e1kbfr)) then
     ABI_CHECK(size(rf_hamkq%e1kbfr,4)==1,'BUG in getgh2c: qphase>1!')
     enl_ptr => rf_hamkq%e1kbfr
   else
     msg='For ipert=natom+11/pert_phon_elfd : e1kbfr and/or e1kbsc must be associated or enl optional input must be present.'
     MSG_BUG(msg)
   end if
   if (usevnl==0) then
     msg='gvnl2 must be allocated for ipert=natom+11/pert_phon_elfd !'
     MSG_BUG(msg)
   end if
   if(opt_gvnl2==0) then
     msg='opt_gvnl2=0 not compatible with ipert=natom+11/pert_phon_elfd !'
     MSG_BUG(msg)
   end if
 end if

!Check sizes
 my_nspinor=max(1,gs_hamkq%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwavef)<2*npw*my_nspinor) then
   msg='wrong size for cwavef!'
   MSG_BUG(msg)
 end if
 if (size(gh2c)<2*npw1*my_nspinor) then
   msg='wrong size for gh2c!'
   MSG_BUG(msg)
 end if
 if (usevnl/=0) then
   if (size(gvnl2)<2*npw1*my_nspinor) then
     msg='wrong size for gvnl2!'
     MSG_BUG(msg)
   end if
 end if
 if (sij_opt==1) then
   if (size(gs2c)<2*npw1*my_nspinor) then
     msg='wrong size for gs2c!'
     MSG_BUG(msg)
   end if
 end if

!PAW: specific treatment for usecprj input arg
!     force it to zero if cwaveprj is not allocated
 usecprj=gs_hamkq%usecprj ; ncpgr=0
 if(gs_hamkq%usepaw==1) then
   if (size(cwaveprj)==0) usecprj=0
   if (usecprj/=0) then
     ncpgr=cwaveprj(1,1)%ncpgr
     if (size(cwaveprj)<gs_hamkq%natom*my_nspinor) then
       msg='wrong size for cwaveprj!'
       MSG_BUG(msg)
     end if
   end if
 else
   if(usecprj==1)then
     msg='usecprj==1 not allowed for NC psps !'
     MSG_BUG(msg)
   end if
 end if

! tim_nonlop=8
! if (tim_getgh2c==1.and.ipert<=natom) tim_nonlop=7
! if (tim_getgh2c==2.and.ipert<=natom) tim_nonlop=5
! if (tim_getgh2c==1.and.ipert> natom) tim_nonlop=8
! if (tim_getgh2c==2.and.ipert> natom) tim_nonlop=5
! if (tim_getgh2c==3                 ) tim_nonlop=0

 idir1=alpha(idir);idir2=beta(idir)

 compute_conjugate = .false.
 if(present(conj)) compute_conjugate = conj

!======================================================================
!== Apply the 2nd-order local potential to the wavefunction
!======================================================================

 if (ipert/=natom+10.and.ipert/=natom+11.and.optlocal>0) then
   msg='local part not implemented'
   MSG_BUG(msg)
 else
!  In the case of ddk operator, no local contribution (also because no self-consistency)
!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gh2c(:,ipw)=zero
   end do

 end if

!======================================================================
!== Apply the 2st-order non-local potential to the wavefunction
!======================================================================

 has_vnl=(ipert==natom+10.or.ipert==natom+11.or.pert_phon_elfd)

!Use of gvnl2 depends on usevnl
 if (usevnl==1) then
   gvnl2_ => gvnl2
 else
   ABI_ALLOCATE(gvnl2_,(2,npw1*my_nspinor))
 end if

 if (has_vnl.and.(optnl>0.or.sij_opt/=0)) then

   idirc=3*(idir1-1)+idir2 !xx=1, xy=2, xz=3, yx=4, yy=5, yz=6, zx=7, zy=8, zz=9, (xyz,xyz)=(idir1,idir2)

! d^2[H_nl]/dk1dk2
!  -------------------------------------------
   if (ipert==natom+10) then
     if (gs_hamkq%usepaw==1) then
       if (usecprj==1) then
         cwaveprj_ptr => cwaveprj
       else
         ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,my_nspinor))
         call pawcprj_alloc(cwaveprj_tmp,0,gs_hamkq%dimcprj)
         cwaveprj_ptr => cwaveprj_tmp
       end if
       cpopt=-1+5*usecprj
       choice=8; signs=2; paw_opt=1; if (sij_opt/=0) paw_opt=sij_opt+3
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idirc,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,gs2c,tim_nonlop,cwavef,gvnl2_)
       if (usecprj==0) then
         call pawcprj_free(cwaveprj_tmp)
         ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
       end if
       nullify(cwaveprj_ptr)
     else
       choice=8; signs=2; cpopt=-1 ; paw_opt=0
       call nonlop(choice,cpopt,cwaveprj,enlout_dum,gs_hamkq,idirc,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwavef,gvnl2_)
     end if

! d^2[H_nl]/dk1dE2 : Non-zero only in PAW
!  -------------------------------------------
   else if (ipert==natom+11.and.gs_hamkq%usepaw==1) then

     ABI_ALLOCATE(nonlop_out,(2,npw1*my_nspinor))

     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,2,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if

     if (opt_gvnl2==1.and.optnl>=1) then

!      Compute application of dS/dk1 to i*d[cwavef]/dk2
!      sum_{i,j} s_ij d(|p_i><p_j|)/dk(idir1) | i*psi^(k(idir2)) >
       cpopt=-1 ; choice=5 ; paw_opt=3 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idir1,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,nonlop_out,tim_nonlop,gvnl2_,vectout_dum)

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=nonlop_out(:,ipw)
       end do

!      Compute part of H^(2) due to derivative of projectors (idir1) and derivative of Dij (idir2)
!      sum_{i,j} chi_ij(idir2) d(|p_i><p_j|)/dk(idir1) | psi^(0) >
       cpopt=4*usecprj ; choice=5 ; paw_opt=1 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idir1,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwavef,nonlop_out,enl=enl_ptr)

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=gvnl2_(:,ipw)+nonlop_out(:,ipw)
       end do

     else

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=zero
       end do

     end if ! opt_gvnl2==1

!    Compute derivatives due to projectors |d^2[p_i]/dk1dk2>,|d[p_i]/dk1>,|d[p_i]/dk2>
!    i * sum_{i,j} (d(|p_i><dp_j/dk(idir2)|)/dk(idir1) | psi^(0) >
     cpopt=-1+5*usecprj ; choice=81 ; paw_opt=3 ; signs=2
     call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idirc,(/zero/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,nonlop_out,tim_nonlop,cwavef,vectout_dum)

     if(compute_conjugate) then
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor ! Note the multiplication by -i
         gvnl2_(1,ipw)=gvnl2_(1,ipw)+nonlop_out(2,ipw)
         gvnl2_(2,ipw)=gvnl2_(2,ipw)-nonlop_out(1,ipw)
       end do
     else
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor ! Note the multiplication by i
         gvnl2_(1,ipw)=gvnl2_(1,ipw)-nonlop_out(2,ipw)
         gvnl2_(2,ipw)=gvnl2_(2,ipw)+nonlop_out(1,ipw)
       end do
     end if

     ABI_DEALLOCATE(nonlop_out)
     if (sij_opt==1) gs2c=zero
     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)

! d^2[H_nl]/dtau1dE2 : Non-zero only in PAW
!  -------------------------------------------
   else if (pert_phon_elfd.and.gs_hamkq%usepaw==1) then

     iatm = ipert-(natom+11)
     if (iatm<1.or.iatm>natom) then
       MSG_BUG(" iatm must be between 1 and natom")
     end if

     ABI_ALLOCATE(nonlop_out,(2,npw1*my_nspinor))

     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_DATATYPE_ALLOCATE(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,2,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if

     if (opt_gvnl2==1) then

!      Compute application of dS/dtau1 to i*d[cwavef]/dk2
!      sum_{i,j} s_ij d(|p_i><p_j|)/dtau(idir1) | i*psi^(k(idir2)) >
       cpopt=-1 ; choice=2 ; paw_opt=3 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idir1,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,nonlop_out,tim_nonlop,gvnl2_,vectout_dum,iatom_only=iatm)

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=nonlop_out(:,ipw)
       end do

!      Compute part of H^(2) due to derivative of projectors (idir1) and derivative of Dij (idir2)
!      sum_{i,j} chi_ij(idir2) d(|p_i><p_j|)/dtau(idir1) | psi^(0) >
       cpopt=4*usecprj ; choice=2 ; paw_opt=1 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idir1,(/zero/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwavef,nonlop_out,enl=enl_ptr,iatom_only=iatm)

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=gvnl2_(:,ipw)+nonlop_out(:,ipw)
       end do

     else

!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnl2_(:,ipw)=zero
       end do

     end if ! opt_gvnl2==1

!    Compute derivatives due to projectors |d^2[p_i]/dtau1dk2>,|d[p_i]/dtau1>,|d[p_i]/dk2>
!    i * sum_{i,j} (d(|p_i><dp_j/dk(idir2)|)/dtau(idir1) | psi^(0) >
     cpopt=-1+5*usecprj ; choice=54 ; paw_opt=3 ; signs=2
     call nonlop(choice,cpopt,cwaveprj_ptr,enlout_dum,gs_hamkq,idirc,(/zero/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,nonlop_out,tim_nonlop,cwavef,vectout_dum,iatom_only=iatm)

     if(compute_conjugate) then
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor ! Note the multiplication by -i
         gvnl2_(1,ipw)=gvnl2_(1,ipw)+nonlop_out(2,ipw)
         gvnl2_(2,ipw)=gvnl2_(2,ipw)-nonlop_out(1,ipw)
       end do
     else
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor ! Note the multiplication by i
         gvnl2_(1,ipw)=gvnl2_(1,ipw)-nonlop_out(2,ipw)
         gvnl2_(2,ipw)=gvnl2_(2,ipw)+nonlop_out(1,ipw)
       end do
     end if

     ABI_DEALLOCATE(nonlop_out)
     if (sij_opt==1) gs2c=zero
     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_DATATYPE_DEALLOCATE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)

   end if

!No non-local part
!-------------------------------------------
 else

   if (optnl>=1) then
 !$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gvnl2_(:,ipw)=zero
     end do
   end if
   if (sij_opt/=0) then
 !$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gs2c(:,ipw)=zero
     end do
   end if

 end if

 if (associated(enl_ptr)) then
   nullify(enl_ptr)
 end if
 if (allocated(enl_temp)) then
   ABI_DEALLOCATE(enl_temp)
 end if

!======================================================================
!== Apply the 2nd-order kinetic operator to the wavefunction
!======================================================================

 if (present(optkin)) then
   has_kin=(optkin/=0.and.ipert==natom+10)
 else
   has_kin=(ipert==natom+10)
 end if

!k-point perturbation
!-------------------------------------------
 if (associated(gs_hamkq%kinpw_kp)) then
   kinpw1 => gs_hamkq%kinpw_kp
 else if (optnl>=1.or.has_kin) then
   msg='need kinpw1 allocated!'
   MSG_BUG(msg)
 end if
 if (associated(rf_hamkq%ddkinpw_k)) then
   ddkinpw => rf_hamkq%ddkinpw_k
 else if (has_kin) then
   msg='need ddkinpw allocated!'
   MSG_BUG(msg)
 end if

 if (has_kin) then
   do ispinor=1,my_nspinor
 !$OMP PARALLEL DO PRIVATE(ipw,ipws) SHARED(cwavef,ispinor,gvnl2_,ddkinpw,kinpw1,npw,my_nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(zero)*1.d-11)then
         gvnl2_(1,ipws)=gvnl2_(1,ipws)+ddkinpw(ipw)*cwavef(1,ipws)
         gvnl2_(2,ipws)=gvnl2_(2,ipws)+ddkinpw(ipw)*cwavef(2,ipws)
       else
         gvnl2_(1,ipws)=zero
         gvnl2_(2,ipws)=zero
       end if
     end do
   end do
 end if

!======================================================================
!== Sum contributions to get the application of H^(2) to the wf
!======================================================================
!Also filter the wavefunctions for large modified kinetic energy

!Add non-local+kinetic to local part
 if (optnl>=1.or.has_kin) then
   do ispinor=1,my_nspinor
     ipws=(ispinor-1)*npw1
 !$OMP PARALLEL DO PRIVATE(ipw) SHARED(gh2c,gvnl2_,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh2c(1,ipw)=gh2c(1,ipw)+gvnl2_(1,ipw)
         gh2c(2,ipw)=gh2c(2,ipw)+gvnl2_(2,ipw)
       else
         gh2c(1,ipw)=zero
         gh2c(2,ipw)=zero
       end if
     end do
   end do
 end if

 if (usevnl==1) then
   nullify(gvnl2_)
 else
   ABI_DEALLOCATE(gvnl2_)
 end if

!call timab(196+tim_getgh2c,2,tsec)

 DBG_EXIT("COLL")

end subroutine getgh2c
!!***

end module m_getgh2c
!!***
