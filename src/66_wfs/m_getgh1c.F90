!!****m* ABINIT/m_getgh1c
!! NAME
!!  m_getgh1c
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (XG, DRH, MT, SPr)
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

module m_getgh1c

 use defs_basis
 use m_abicore
 use m_errors
 use m_dtset

 use defs_abitypes, only : MPI_type
 use defs_datatypes, only : pseudopotential_type
 use m_time,        only : timab, cwtime, cwtime_report
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, pawcprj_lincom, pawcprj_axpby
 use m_kg,          only : kpgstr, mkkin, mkkpg, mkkin_metdqdq
 use m_mkffnl,      only : mkffnl
 use m_pawfgr,      only : pawfgr_type
 use m_fft,         only : fftpac, fourwf
 use m_hamiltonian, only : gs_hamiltonian_type, rf_hamiltonian_type
 use m_cgtools,          only : projbd
 use m_nonlop,           only : nonlop
 use m_fourier_interpol, only : transgrid

 implicit none

 private
!!***

 public :: getgh1c
 public :: rf_transgrid_and_pack
 public :: getgh1c_setup
 public :: getdc1
 public :: getgh1dqc
 public :: getgh1dqc_setup
!!***

contains
!!***

!!****f* ABINIT/getgh1c
!!
!! NAME
!! getgh1c
!!
!! FUNCTION
!! Compute <G|H^(1)|C> (or <G|H^(1)-lambda.S^(1)|C>) for input vector |C> expressed in reciprocal space.
!! (H^(1) is the 1st-order pertubed Hamiltonian, S^(1) is the 1st-order perturbed overlap operator).
!! Result is put in array gh1c.
!! If required, part of <G|K(1)+Vnonlocal^(1)|C> not depending on VHxc^(1) is also returned in gvnlx1c.
!! If required, <G|S^(1)|C> is returned in gs1c (S=overlap - PAW only)
!!
!! INPUTS
!!  berryopt=option for Berry phase
!!  cwave(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!     if not allocated or size=0, they are locally computed (and not sorted)
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
!!  grad_berry(2,npw1*nspinor*(berryopt/4))= the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  lambda=real use to apply H^(1)-lambda.S^(1)
!!  mpi_enreg=information about MPI parallelization
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  optlocal=0: local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
!!           1: local part of H^(1) is computed in gh1c=<G|H^(1)|C>
!!  optnl=0: non-local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
!!        1: non-local part of H^(1) depending on VHxc^(1) is not computed in gh1c=<G|H^(1)|C>
!!        2: non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
!!  opt_gvnlx1=option controlling the use of gvnlx1 array:
!!            0: used as an output
!!            1: used as an input:   (only for ipert=natom+2)
!!                 NCPP: contains the ddk 1-st order WF
!!                 PAW: contains frozen part of 1st-order hamiltonian
!!            2: used as input/ouput:    - used only for PAW and ipert=natom+2
!!                 At input: contains the ddk 1-st order WF (times i)
!!                 At output: contains frozen part of 1st-order hamiltonian
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,k+q
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H^(1)|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S^(1)|C> have to be computed in gs1c in addition to gh1c
!!                       if -1, matrix elements <G|H^(1)-lambda.S^(1)|C> have to be computed in gh1c (gs1c not used)
!!  tim_getgh1c=timing code of the calling subroutine (can be set to 0 if not attributed)
!!  usevnl=1 if gvnlx1=(part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1)) has to be input/output
!!
!! OUTPUT
!! gh1c(2,npw1*nspinor)= <G|H^(1)|C> or  <G|H^(1)-lambda.S^(1)|C> on the k+q sphere
!!                     (only kinetic+non-local parts if optlocal=0)
!! if (usevnl==1)
!!  gvnlx1(2,npw1*nspinor)=  part of <G|K^(1)+Vnl^(1)|C> not depending on VHxc^(1)              (sij_opt/=-1)
!!                       or part of <G|K^(1)+Vnl^(1)-lambda.S^(1)|C> not depending on VHxc^(1) (sij_opt==-1)
!! if (sij_opt=1)
!!  gs1c(2,npw1*nspinor)=<G|S^(1)|C> (S=overlap) on the k+q sphere.
!!
!! PARENTS
!!      dfpt_cgwf,dfpt_nstpaw,dfpt_nstwf,dfpt_wfkfermi,m_gkk,m_phgamma,m_phpi
!!      m_rf2,m_sigmaph
!!
!! CHILDREN
!!      kpgstr,load_k_hamiltonian,load_k_rf_hamiltonian,load_kprime_hamiltonian
!!      mkffnl,mkkin,mkkpg
!!
!! SOURCE

subroutine getgh1c(berryopt,cwave,cwaveprj,gh1c,grad_berry,gs1c,gs_hamkq,&
&          gvnlx1,idir,ipert,lambda,mpi_enreg,optlocal,optnl,opt_gvnlx1,&
&          rf_hamkq,sij_opt,tim_getgh1c,usevnl,conj)

!Arguments ------------------------------------
!scalars
 logical,intent(in),optional :: conj
 integer,intent(in) :: berryopt,idir,ipert,optlocal,optnl,opt_gvnlx1,sij_opt,tim_getgh1c,usevnl
 real(dp),intent(in) :: lambda
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq
!arrays
 real(dp),intent(in) :: grad_berry(:,:)
 real(dp),intent(inout) :: cwave(2,gs_hamkq%npw_k*gs_hamkq%nspinor)
 real(dp),intent(out) :: gh1c(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 real(dp),intent(out) :: gs1c(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 real(dp),intent(inout),target :: gvnlx1(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=16
 integer :: choice,cplex1,cpopt,ipw,ipws,ispinor,istr,i1,i2,i3
 integer :: my_nspinor,natom,ncpgr,nnlout=1,npw,npw1,paw_opt,signs
 integer :: tim_fourwf,tim_nonlop,usecprj
 logical :: compute_conjugate,has_kin,usevnl2
 real(dp) :: weight !, cpu, wall, gflops
 !character(len=500) :: msg
!arrays
 real(dp) :: enlout(1),tsec(2),svectout_dum(1,1),vectout_dum(1,1)
 real(dp),allocatable :: cwave_sp(:,:),cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: gh1c_sp(:,:),gh1c1(:,:),gh1c2(:,:),gh1c3(:,:),gh1c4(:,:),gvnl2(:,:)
 real(dp),allocatable :: nonlop_out(:,:),vlocal1_tmp(:,:,:),work(:,:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: gvnlx1_(:,:)
 real(dp),pointer :: dkinpw(:),kinpw1(:)
 type(pawcprj_type),allocatable,target :: cwaveprj_tmp(:,:)
 type(pawcprj_type),pointer :: cwaveprj_ptr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 ! Keep track of total time spent in getgh1c
 call timab(196+tim_getgh1c,1,tsec)

 !call cwtime(cpu, wall, gflops, "start")
 !call cwtime_report(" getgh1c", cpu, wall, gflops)

!======================================================================
!== Initialisations and compatibility tests
!======================================================================

 npw  =gs_hamkq%npw_k
 npw1 =gs_hamkq%npw_kp
 natom=gs_hamkq%natom

 ! Compatibility tests
 if(gs_hamkq%usepaw==1.and.(ipert>=0.and.(ipert<=natom.or.ipert==natom+3.or.ipert==natom+4))) then
   if ((optnl>=1.and.(.not.associated(rf_hamkq%e1kbfr))) .or. &
       (optnl==2.and.(.not.associated(rf_hamkq%e1kbsc)))) then
     MSG_BUG('ekb derivatives must be allocated for ipert<=natom or natom+3/4 !')
   end if
 end if
 if(gs_hamkq%usepaw==1.and.(ipert==natom+2)) then
   if ((optnl>=1.and.(.not.associated(rf_hamkq%e1kbfr))) .or. &
       (optnl==2.and.(.not.associated(rf_hamkq%e1kbsc)))) then
     MSG_BUG('ekb derivatives must be allocated for ipert=natom+2 !')
   end if
   if (usevnl==0) then
     MSG_BUG('gvnlx1 must be allocated for ipert=natom+2 !')
   end if
 end if
 if(ipert==natom+2.and.opt_gvnlx1==0) then
   MSG_BUG('opt_gvnlx1=0 not compatible with ipert=natom+2 !')
 end if
 if (mpi_enreg%paral_spinor==1) then
   MSG_BUG('Not compatible with parallelization over spinorial components !')
 end if

 ! Check sizes
 my_nspinor=max(1,gs_hamkq%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwave)<2*npw*my_nspinor) then
   MSG_BUG('wrong size for cwave!')
 end if
 if (size(gh1c)<2*npw1*my_nspinor) then
   MSG_BUG('wrong size for gh1c!')
 end if
 if (usevnl/=0) then
   if (size(gvnlx1)<2*npw1*my_nspinor) then
     MSG_BUG('wrong size for gvnlx1!')
   end if
 end if
 if (sij_opt==1) then
   if (size(gs1c)<2*npw1*my_nspinor) then
     MSG_BUG('wrong size for gs1c!')
   end if
 end if
 if (berryopt>=4) then
   if (size(grad_berry)<2*npw1*my_nspinor) then
     MSG_BUG('wrong size for grad_berry!')
   end if
 end if

 ! PAW: specific treatment for usecprj input arg
 !      force it to zero if cwaveprj is not allocated
 usecprj=gs_hamkq%usecprj ; ncpgr=0
 if(gs_hamkq%usepaw==1) then
   if (size(cwaveprj)==0) usecprj=0
   if (usecprj/=0) then
     ncpgr=cwaveprj(1,1)%ncpgr
     if (size(cwaveprj)<gs_hamkq%natom*my_nspinor) then
       MSG_BUG('wrong size for cwaveprj!')
     end if
     if(gs_hamkq%usepaw==1.and.(ipert>=0.and.(ipert<=natom.or.ipert==natom+3.or.ipert==natom+4))) then
       if (ncpgr/=1)then
         MSG_BUG('Projected WFs (cprj) derivatives are not correctly stored !')
       end if
     end if
   end if
 else
   if(usecprj==1)then
     MSG_BUG('usecprj==1 not allowed for NC psps !')
   end if
 end if

 tim_nonlop=8
 if (tim_getgh1c==1.and.ipert<=natom) tim_nonlop=7
 if (tim_getgh1c==2.and.ipert<=natom) tim_nonlop=5
 if (tim_getgh1c==1.and.ipert> natom) tim_nonlop=8
 if (tim_getgh1c==2.and.ipert> natom) tim_nonlop=5
 if (tim_getgh1c==3                 ) tim_nonlop=0

 compute_conjugate = .false.
 if(present(conj)) compute_conjugate = conj

!======================================================================
!== Apply the 1st-order local potential to the wavefunction
!======================================================================

!Phonon perturbation
!or Electric field perturbation
!or Strain perturbation
!-------------------------------------------
 if (ipert<=natom+5.and.ipert/=natom+1.and.optlocal>0) then !SPr deb

   ABI_ALLOCATE(work,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))

   if (gs_hamkq%nvloc==1) then

     weight=one ; tim_fourwf=4
     call fourwf(rf_hamkq%cplex,rf_hamkq%vlocal1,cwave,gh1c,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&     gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&     use_gpu_cuda=gs_hamkq%use_gpu_cuda)
     if(gs_hamkq%nspinor==2)then
       ABI_ALLOCATE(cwave_sp,(2,npw))
       ABI_ALLOCATE(gh1c_sp,(2,npw1))
!$OMP PARALLEL DO
       do ipw=1,npw
         cwave_sp(1,ipw)=cwave(1,ipw+npw)
         cwave_sp(2,ipw)=cwave(2,ipw+npw)
       end do
       call fourwf(rf_hamkq%cplex,rf_hamkq%vlocal1,cwave_sp,gh1c_sp,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!$OMP PARALLEL DO
       do ipw=1,npw1
         gh1c(1,ipw+npw1)=gh1c_sp(1,ipw)
         gh1c(2,ipw+npw1)=gh1c_sp(2,ipw)
       end do
       ABI_FREE(cwave_sp)
       ABI_FREE(gh1c_sp)
     end if
   else ! Non-Collinear magnetism for nvloc=4
     if (gs_hamkq%nspinor==2) then
       weight=one ; tim_fourwf=4
       ABI_ALLOCATE(gh1c1,(2,npw1))
       ABI_ALLOCATE(gh1c2,(2,npw1))
       ABI_ALLOCATE(gh1c3,(2,npw1))
       ABI_ALLOCATE(gh1c4,(2,npw1))
       gh1c1(:,:)=zero; gh1c2(:,:)=zero; gh1c3(:,:)=zero ;  gh1c4(:,:)=zero
       ABI_ALLOCATE(vlocal1_tmp,(rf_hamkq%cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6)) !SPr: notation/dimension corrected vlocal_tmp -> vlocal1_tmp
       ABI_ALLOCATE(cwavef1,(2,npw))
       ABI_ALLOCATE(cwavef2,(2,npw))
       do ipw=1,npw
         cwavef1(1:2,ipw)=cwave(1:2,ipw)
         cwavef2(1:2,ipw)=cwave(1:2,ipw+npw)
       end do
!      gh1c1=v11*phi1
       vlocal1_tmp(:,:,:)=rf_hamkq%vlocal1(:,:,:,1)
       call fourwf(rf_hamkq%cplex,vlocal1_tmp,cwavef1,gh1c1,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      gh1c2=v22*phi2
       vlocal1_tmp(:,:,:)=rf_hamkq%vlocal1(:,:,:,2)
       call fourwf(rf_hamkq%cplex,vlocal1_tmp,cwavef2,gh1c2,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       ABI_FREE(vlocal1_tmp)
       cplex1=2
       ABI_ALLOCATE(vlocal1_tmp,(cplex1*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))
!      gh1c3=(re(v12)-im(v12))*phi1 => v^21*phi1
       if(rf_hamkq%cplex==1) then
         do i3=1,gs_hamkq%n6
           do i2=1,gs_hamkq%n5
             do i1=1,gs_hamkq%n4
               vlocal1_tmp(2*i1-1,i2,i3)= rf_hamkq%vlocal1(i1,i2,i3,3)
               vlocal1_tmp(2*i1  ,i2,i3)=-rf_hamkq%vlocal1(i1,i2,i3,4)
             end do
           end do
         end do
       else
       !SPr: modified definition of local potential components for cplex=2 (see dotprod_vn)
       !also, v21==v12* not always holds (e.g. magnetic field perturbation)
         do i3=1,gs_hamkq%n6
           do i2=1,gs_hamkq%n5
             do i1=1,gs_hamkq%n4
               vlocal1_tmp(2*i1-1,i2,i3)= rf_hamkq%vlocal1(2*i1  ,i2,i3,4)
               vlocal1_tmp(2*i1  ,i2,i3)=-rf_hamkq%vlocal1(2*i1-1,i2,i3,4)
             end do
           end do
         end do
       end if
       call fourwf(cplex1,vlocal1_tmp,cwavef1,gh1c3,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
!      gh1c4=(re(v12)+im(v12))*phi2 => v^12*phi2
       if(rf_hamkq%cplex==1) then
         do i3=1,gs_hamkq%n6
           do i2=1,gs_hamkq%n5
             do i1=1,gs_hamkq%n4
               vlocal1_tmp(2*i1,i2,i3)=-vlocal1_tmp(2*i1,i2,i3)
             end do
           end do
         end do
       else
         !for cplex=2 and time-reversal breaking perturbations,v21/=v12*
         do i3=1,gs_hamkq%n6
           do i2=1,gs_hamkq%n5
             do i1=1,gs_hamkq%n4
               vlocal1_tmp(2*i1-1,i2,i3)= rf_hamkq%vlocal1(2*i1-1,i2,i3,3)
               vlocal1_tmp(2*i1  ,i2,i3)= rf_hamkq%vlocal1(2*i1  ,i2,i3,3)
             end do
           end do
         end do
       end if
       call fourwf(cplex1,vlocal1_tmp,cwavef2,gh1c4,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
&       gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
&       use_gpu_cuda=gs_hamkq%use_gpu_cuda)
       ABI_FREE(vlocal1_tmp)
!      Build gh1c from pieces
!      gh1c_1 = (v11, v12) (psi1) matrix vector product
!      gh1c_2 = (v12*,v22) (psi2)
       do ipw=1,npw1
         gh1c(1:2,ipw)     =gh1c1(1:2,ipw)+gh1c4(1:2,ipw)
         gh1c(1:2,ipw+npw1)=gh1c3(1:2,ipw)+gh1c2(1:2,ipw)
       end do
       ABI_FREE(gh1c1)
       ABI_FREE(gh1c2)
       ABI_FREE(gh1c3)
       ABI_FREE(gh1c4)
       ABI_FREE(cwavef1)
       ABI_FREE(cwavef2)
     else
       MSG_BUG('nspinor/=1 for Non-collinear calculations!')
     end if
   end if ! nvloc

   ABI_FREE(work)

!  k-point perturbation (or no local part, i.e. optlocal=0)
!  -------------------------------------------
 else if (ipert==natom+1.or.optlocal==0) then

!  In the case of ddk operator, no local contribution (also because no self-consistency)
!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gh1c(:,ipw)=zero
   end do

 end if

!======================================================================
!== Apply the 1st-order non-local potential to the wavefunction
!======================================================================

!Use of gvnlx1 depends on usevnl
 if (usevnl==1) then
   gvnlx1_ => gvnlx1
 else
   ABI_ALLOCATE(gvnlx1_,(2,npw1*my_nspinor))
 end if

!Phonon perturbation
!-------------------------------------------
 if (ipert<=natom.and.(optnl>0.or.sij_opt/=0)) then

!  PAW:
   if (gs_hamkq%usepaw==1) then

     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_MALLOC(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,1,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if

!    1- Compute derivatives due to projectors |p_i>^(1)
!    Only displaced atom contributes
     cpopt=-1+5*usecprj ; choice=2 ; signs=2
     paw_opt=1;if (sij_opt/=0) paw_opt=sij_opt+3
     call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,gs1c,tim_nonlop,cwave,gvnlx1_,iatom_only=ipert)

!    2- Compute derivatives due to frozen part of D_ij^(1) (independent of VHxc^(1))
!    All atoms contribute
     if (optnl>=1) then
       ABI_ALLOCATE(nonlop_out,(2,npw1*my_nspinor))
       cpopt=1+3*usecprj ; choice=1 ; signs=2 ; paw_opt=1
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,nonlop_out,enl=rf_hamkq%e1kbfr)
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnlx1_(:,ipw)=gvnlx1_(:,ipw)+nonlop_out(:,ipw)
       end do
       ABI_FREE(nonlop_out)
     end if

!    3- Compute derivatives due to self-consistent part of D_ij^(1) (depending on VHxc^(1))
!    All atoms contribute
     if (optnl==2) then
       ABI_ALLOCATE(gvnl2,(2,npw1*my_nspinor))
       cpopt=4 ; choice=1 ; signs=2 ; paw_opt=1
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl2,enl=rf_hamkq%e1kbsc)
     end if

     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_FREE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)

!  Norm-conserving psps:
   else
!    Compute only derivatives due to projectors |p_i>^(1)
     cpopt=-1 ; choice=2 ; signs=2 ; paw_opt=0
     call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnlx1_,iatom_only=ipert)
     if (sij_opt==1) then
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gs1c(:,ipw)=zero
       end do
     end if
   end if

!  k-point perturbation
!  -------------------------------------------
 else if (ipert==natom+1.and.(optnl>0.or.sij_opt/=0)) then

   tim_nonlop=8 ; signs=2 ; choice=5
   if (gs_hamkq%usepaw==1) then
     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_MALLOC(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,1,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if
     cpopt=-1+5*usecprj; paw_opt=1; if (sij_opt/=0) paw_opt=sij_opt+3
!    JLJ: BUG (wrong result) of H^(1) if stored cprj are used in PAW DDKs with nspinor==2 (==1 works fine).
!    To be debugged, if someone has time...
     if(gs_hamkq%nspinor==2) cpopt=-1
     call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,gs1c,tim_nonlop,cwave,gvnlx1_)
     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_FREE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)
   else
     cpopt=-1 ; paw_opt=0
     call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnlx1_)
   end if

!  Electric field perturbation without Berry phase
!  -------------------------------------------
 else if(ipert==natom+2 .and. &
&   (berryopt/=4 .and. berryopt/=6 .and. berryopt/=7 .and. &
&   berryopt/=14 .and. berryopt/=16 .and. berryopt/=17) .and.(optnl>0.or.sij_opt/=0))then
!  gvnlx1 was already initialized in the calling routine, by reading a ddk file
!  It contains |i du^(0)/dk_band>

   if (gs_hamkq%usepaw==1) then
     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_MALLOC(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,1,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if
     if (opt_gvnlx1==2.and.optnl>=1) then

!      PAW: Compute application of S^(0) to ddk WF
       cpopt=-1 ; choice=1 ; paw_opt=3 ; signs=2
       ABI_ALLOCATE(nonlop_out,(2,npw1*my_nspinor))
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,0,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,nonlop_out,tim_nonlop,gvnlx1_,vectout_dum)
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnlx1_(:,ipw)=nonlop_out(:,ipw)
       end do

!      PAW: Compute part of H^(1) due to derivative of S
       cpopt=4*usecprj ; choice=51 ; paw_opt=3 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,nonlop_out,tim_nonlop,cwave,vectout_dum)
       if(compute_conjugate) then
!$OMP PARALLEL DO
         do ipw=1,npw1*my_nspinor ! Note the multiplication by -i
           gvnlx1_(1,ipw)=gvnlx1_(1,ipw)+nonlop_out(2,ipw)
           gvnlx1_(2,ipw)=gvnlx1_(2,ipw)-nonlop_out(1,ipw)
         end do
       else
!$OMP PARALLEL DO
         do ipw=1,npw1*my_nspinor ! Note the multiplication by i
           gvnlx1_(1,ipw)=gvnlx1_(1,ipw)-nonlop_out(2,ipw)
           gvnlx1_(2,ipw)=gvnlx1_(2,ipw)+nonlop_out(1,ipw)
         end do
       end if

!      PAW: Compute part of H^(1) due to derivative of electric field part of Dij
       cpopt=2 ; choice=1 ; paw_opt=1 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,0,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,nonlop_out,enl=rf_hamkq%e1kbfr)
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnlx1_(:,ipw)=gvnlx1_(:,ipw)+nonlop_out(:,ipw)
       end do
       ABI_FREE(nonlop_out)

     end if ! opt_gvnlx1==2

!    PAW: Compute derivatives due to part of D_ij^(1) depending on VHxc^(1)
     if (optnl>=2) then
       ABI_ALLOCATE(gvnl2,(2,npw1*my_nspinor))
       cpopt=-1+3*usecprj;if (opt_gvnlx1==2) cpopt=2
       choice=1 ; paw_opt=1 ; signs=2
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,0,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl2,enl=rf_hamkq%e1kbsc)
     end if

     if (sij_opt==1) then
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gs1c(:,ipw)=zero
       end do
     end if
     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_FREE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)
   end if  ! PAW

!  Electric field perturbation with Berry phase
!  -------------------------------------------
 else if(ipert==natom+2 .and. &
&   (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. &
&   berryopt==14 .or. berryopt==16 .or. berryopt==17 ) .and.(optnl>0.or.sij_opt/=0))then

   if (optnl>=1) then
     do ipw=1,npw1*my_nspinor
       gvnlx1_(1,ipw)=-grad_berry(2,ipw)
       gvnlx1_(2,ipw)= grad_berry(1,ipw)
     end do
   end if
   if (sij_opt==1) then
!$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gs1c(:,ipw)=zero
     end do
   end if

!  Strain perturbation
!  -------------------------------------------
 else if ((ipert==natom+3.or.ipert==natom+4).and.(optnl>0.or.sij_opt/=0)) then

   istr=idir;if(ipert==natom+4) istr=istr+3

!  PAW:
   if (gs_hamkq%usepaw==1) then

     if (usecprj==1) then
       cwaveprj_ptr => cwaveprj
     else
       ABI_MALLOC(cwaveprj_tmp,(natom,my_nspinor))
       call pawcprj_alloc(cwaveprj_tmp,1,gs_hamkq%dimcprj)
       cwaveprj_ptr => cwaveprj_tmp
     end if

!    1- Compute derivatives due to projectors |p_i>^(1)
!    All atoms contribute
     cpopt=-1+5*usecprj ; choice=3 ; signs=2
     paw_opt=1;if (sij_opt/=0) paw_opt=sij_opt+3
     call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,istr,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,gs1c,tim_nonlop,cwave,gvnlx1_)

!    2- Compute derivatives due to frozen part of D_ij^(1) (independent of VHxc^(1))
!    All atoms contribute
     if (optnl>=1) then
       ABI_ALLOCATE(nonlop_out,(2,npw1*my_nspinor))
       cpopt=1+3*usecprj ; choice=1 ; signs=2 ; paw_opt=1
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,istr,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,nonlop_out,enl=rf_hamkq%e1kbfr)
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gvnlx1_(:,ipw)=gvnlx1_(:,ipw)+nonlop_out(:,ipw)
       end do
       ABI_FREE(nonlop_out)
     end if

!    3- Compute derivatives due to part of D_ij^(1) depending on VHxc^(1)
!    All atoms contribute
     if (optnl>=2) then
       ABI_ALLOCATE(gvnl2,(2,npw1*my_nspinor))
       cpopt=4 ; choice=1 ; signs=2 ; paw_opt=1
       call nonlop(choice,cpopt,cwaveprj_ptr,enlout,gs_hamkq,istr,(/lambda/),mpi_enreg,1,nnlout,&
&       paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl2,enl=rf_hamkq%e1kbsc)
     end if

     if (usecprj==0) then
       call pawcprj_free(cwaveprj_tmp)
       ABI_FREE(cwaveprj_tmp)
     end if
     nullify(cwaveprj_ptr)

!    Norm-conserving psps:
   else
!    Compute only derivatives due to projectors |p_i>^(1)
     choice=3 ; cpopt=-1 ; signs=2 ; paw_opt=0
     call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,istr,(/lambda/),mpi_enreg,1,nnlout,&
&     paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnlx1_)
     if (sij_opt==1) then
!$OMP PARALLEL DO
       do ipw=1,npw1*my_nspinor
         gs1c(:,ipw)=zero
       end do
     end if
   end if

!  No non-local part
!  -------------------------------------------
 else if (usevnl>0.or.(sij_opt/=0)) then

   if (optnl>=1) then
!$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gvnlx1_(:,ipw)=zero
     end do
   end if
   if (sij_opt/=0) then
!$OMP PARALLEL DO
     do ipw=1,npw1*my_nspinor
       gs1c(:,ipw)=zero
     end do
   end if

 end if

!======================================================================
!== Apply the 1st-order kinetic operator to the wavefunction
!== (add it to nl contribution)
!======================================================================

!Phonon perturbation or Electric field perturbation
!-------------------------------------------
!No kinetic contribution

!k-point perturbation or Strain perturbation
!-------------------------------------------

 usevnl2=allocated(gvnl2)
 has_kin=(ipert==natom+1.or.ipert==natom+3.or.ipert==natom+4)
 if (associated(gs_hamkq%kinpw_kp)) then
   kinpw1 => gs_hamkq%kinpw_kp
 else if (optnl>=1.or.usevnl2.or.has_kin) then
   MSG_BUG('need kinpw1 allocated!')
 end if
 if (associated(rf_hamkq%dkinpw_k)) then
   dkinpw => rf_hamkq%dkinpw_k
 else if (has_kin) then
   MSG_BUG('need dkinpw allocated!')
 end if

 if (has_kin) then
!  Remember that npw=npw1 for ddk perturbation
   do ispinor=1,my_nspinor
!$OMP PARALLEL DO PRIVATE(ipw,ipws) SHARED(cwave,ispinor,gvnlx1_,dkinpw,kinpw1,npw,my_nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(zero)*1.d-11)then
         gvnlx1_(1,ipws)=gvnlx1_(1,ipws)+dkinpw(ipw)*cwave(1,ipws)
         gvnlx1_(2,ipws)=gvnlx1_(2,ipws)+dkinpw(ipw)*cwave(2,ipws)
       else
         gvnlx1_(1,ipws)=zero
         gvnlx1_(2,ipws)=zero
       end if
     end do
   end do
 end if

!======================================================================
!== Sum contributions to get the application of H^(1) to the wf
!======================================================================
!Also filter the wavefunctions for large modified kinetic energy

!Add non-local+kinetic to local part
 if (optnl>=1.or.has_kin) then
   do ispinor=1,my_nspinor
     ipws=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gh1c,gvnlx1_,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh1c(1,ipw)=gh1c(1,ipw)+gvnlx1_(1,ipw)
         gh1c(2,ipw)=gh1c(2,ipw)+gvnlx1_(2,ipw)
       else
         gh1c(1,ipw)=zero
         gh1c(2,ipw)=zero
       end if
     end do
   end do
 end if

!PAW: add non-local part due to first order change of VHxc
 if (usevnl2) then
   do ispinor=1,my_nspinor
     ipws=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gh1c,gvnl2,kinpw1,ipws,npw1)
     do ipw=1+ipws,npw1+ipws
       if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
         gh1c(1,ipw)=gh1c(1,ipw)+gvnl2(1,ipw)
         gh1c(2,ipw)=gh1c(2,ipw)+gvnl2(2,ipw)
       end if
     end do
   end do
   ABI_FREE(gvnl2)
 end if

 if (usevnl==1) then
   nullify(gvnlx1_)
 else
   ABI_FREE(gvnlx1_)
 end if

 call timab(196+tim_getgh1c,2,tsec)

 DBG_EXIT("COLL")

end subroutine getgh1c
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/rf_transgrid_and_pack
!! NAME
!!  rf_transgrid_and_pack
!!
!! FUNCTION
!! Set up local potential vlocal1 with proper dimensioning, from vtrial1
!! taking into account the spin. Same thing for vlocal from vtrial.
!!
!! INPUTS
!!  isppol=Spin index.
!!  nspden=Number of density components
!!  usepaw=1 if PAW, 0 for NC.
!!  cplex=1 if DFPT potential is real, 2 for complex
!!  nfftf=Number of FFT points on the FINE grid treated by this processor
!!  nfft=Number of FFT points on the COARSE grid treated by this processor
!!  ngfft(18)=Info on the coarse grid.
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  mpi_enreg=information about MPI parallelization
!!  vtrial(nfftf,nspden)=GS Vtrial(r) on the DENSE mesh
!!  vtrial1(cplex*nfftf,nspden)=INPUT RF Vtrial(r) on the DENSE mesh
!!
!! OUTPUT
!!  vlocal(n4,n5,n6,nvloc)= GS local potential in real space, on the augmented coarse fft grid
!!  vlocal1(cplex*n4,n5,n6,nvloc)= RF local potential in real space, on the augmented coarse fft grid
!!
!! PARENTS
!!      dfpt_vtorho,m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!      kpgstr,load_k_hamiltonian,load_k_rf_hamiltonian,load_kprime_hamiltonian
!!      mkffnl,mkkin,mkkpg
!!
!! SOURCE

subroutine rf_transgrid_and_pack(isppol,nspden,usepaw,cplex,nfftf,nfft,ngfft,nvloc,&
&                                pawfgr,mpi_enreg,vtrial,vtrial1,vlocal,vlocal1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: isppol,nspden,usepaw,cplex,nfftf,nfft,nvloc
 type(pawfgr_type),intent(in) :: pawfgr
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in),target :: vtrial(nfftf,nspden)
 real(dp),intent(inout),target :: vtrial1(cplex*nfftf,nspden)
 real(dp),intent(out) :: vlocal(ngfft(4),ngfft(5),ngfft(6),nvloc)
 real(dp),intent(out) :: vlocal1(cplex*ngfft(4),ngfft(5),ngfft(6),nvloc)

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,n4,n5,n6,paral_kgb,ispden
!arrays
 real(dp) :: rhodum(1), tsec(2)
 real(dp), ABI_CONTIGUOUS pointer :: vtrial_ptr(:,:),vtrial1_ptr(:,:)
 real(dp),allocatable :: cgrvtrial(:,:),cgrvtrial1(:,:),vlocal_tmp(:,:,:),vlocal1_tmp(:,:,:)

! *************************************************************************

 !call timab(1904, 1, tsec)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)
 paral_kgb = mpi_enreg%paral_kgb

 if (nspden/=4) then
   vtrial_ptr => vtrial
   if (usepaw==0.or.pawfgr%usefinegrid==0) then
     call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfft,vtrial_ptr,vlocal(:,:,:,1),2)
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,vtrial1,vlocal1(:,:,:,1),2)
   else
     ABI_ALLOCATE(cgrvtrial,(nfft,nspden))
     call transgrid(1,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial_ptr)
     call fftpac(isppol,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfft,cgrvtrial,vlocal(:,:,:,1),2)
     ABI_FREE(cgrvtrial)
     ABI_ALLOCATE(cgrvtrial,(cplex*nfft,nspden))
     call transgrid(cplex,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial1)
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,cgrvtrial,vlocal1(:,:,:,1),2)
     ABI_FREE(cgrvtrial)
   end if
   nullify(vtrial_ptr)
 else
   ! nspden==4 non-collinear magnetism
   vtrial_ptr => vtrial
   vtrial1_ptr => vtrial1
   ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
   ABI_ALLOCATE(vlocal1_tmp,(cplex*n4,n5,n6))
   if (usepaw==0.or.pawfgr%usefinegrid==0) then
     do ispden=1,nspden
       call fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfft,vtrial_ptr,vlocal_tmp,2)
       vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       call fftpac(ispden,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,ngfft,vtrial1_ptr,vlocal1_tmp,2)
       vlocal1(:,:,:,ispden)=vlocal1_tmp(:,:,:)
     end do
   else
     ! TODO FR EB check the correctness of the following lines for PAW calculations
     ABI_ALLOCATE(cgrvtrial,(nfft,nspden))
     ABI_ALLOCATE(cgrvtrial1,(nfft,nspden))
     call transgrid(cplex,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial_ptr)
     call transgrid(cplex,mpi_enreg,nspden,-1,0,0,paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial1,vtrial1_ptr)
     do ispden=1,nspden
       call fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfft,vtrial_ptr,vlocal_tmp,2)
       vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
       call fftpac(ispden,mpi_enreg,nspden,n1,n2,n3,n4,n5,n6,ngfft,vtrial1_ptr,vlocal1_tmp,2)
       vlocal1(:,:,:,ispden)=vlocal1_tmp(:,:,:)
     end do
     ABI_FREE(cgrvtrial)
   end if
   ABI_FREE(vlocal_tmp)
   ABI_FREE(vlocal1_tmp)
 end if !nspden

 !call timab(1904, 2, tsec)

end subroutine rf_transgrid_and_pack
!!***

!----------------------------------------------------------------------

!!****f* m_hamiltonian/getgh1c_setup
!! NAME
!!  getgh1c_setup
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_vtorho,m_gkk,m_phgamma,m_phpi,m_sigmaph
!!
!! CHILDREN
!!      kpgstr,load_k_hamiltonian,load_k_rf_hamiltonian,load_kprime_hamiltonian
!!      mkffnl,mkkin,mkkpg
!!
!! SOURCE

subroutine getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kpoint,kpq,idir,ipert,&          ! In
                natom,rmet,gprimd,gmet,istwf_k,npw_k,npw1_k,&                          ! In
                useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                           ! In
                dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&         ! Out
                ddkinpw,dkinpw2,rf_hamk_dir2,ffnl1_test)                               ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,istwf_k,npw_k,npw1_k,natom,useylmgr1
 integer,intent(out) :: nkpg,nkpg1
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
 type(rf_hamiltonian_type),intent(inout),optional :: rf_hamk_dir2
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k),kg1_k(3,npw1_k)
 real(dp),intent(in) :: kpoint(3),kpq(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1_k(npw1_k,3+6*((ipert-natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(in) :: ylm1_k(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),allocatable,intent(out) :: dkinpw(:),kinpw1(:),ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable,intent(out),optional :: dkinpw2(:),ddkinpw(:),ffnl1_test(:,:,:,:)
 real(dp),allocatable,intent(out) :: kpg_k(:,:),kpg1_k(:,:),ph3d(:,:,:),ph3d1(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: dimffnl1,dimffnlk,ider,idir0,idir1,idir2,istr,ntypat,print_info
 logical :: qne0
 !real(dp) :: cpu, wall, gflops
!arrays
 real(dp) :: ylmgr_dum(1,1,1), tsec(2)

! *************************************************************************

 ! MG: This routine is called **many times** in the EPH code for phonon and DDK perturbations
 ! Please, be extremely careful when adding extra stuff that may affect performance.

 ! Keep track of total time spent in getgh1c_setup (use 195 slot)
 call timab(195, 1, tsec)
 !call cwtime(cpu, wall, gflops, "start")

 if(.not.present(ddkinpw) .and. ipert==natom+10) then
   MSG_BUG("ddkinpw is not optional for ipert=natom+10.")
 end if
 if(.not.present(dkinpw2) .and. ipert==natom+10 .and. idir>3) then
   MSG_BUG("dkinpw2 is not optional for ipert=natom+10 and idir>3.")
 end if
 if(.not.present(rf_hamk_dir2) .and. ((ipert==natom+10 .and. idir>3) .or. ipert==natom+11)) then
   MSG_BUG("rf_hamk_dir2 is not optional for ipert=natom+10 (with idir>3) or ipert=natom+11.")
 end if

 ntypat = psps%ntypat
 qne0=((kpq(1)-kpoint(1))**2+(kpq(2)-kpoint(2))**2+(kpq(3)-kpoint(3))**2>=tol14)

 ! Compute k+G vectors
 nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(3)
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

 ! Compute k+q+G vectors
 nkpg1=0;if(ipert>=1.and.ipert<=natom) nkpg1=3*dtset%nloalg(3)
 ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
 if (nkpg1>0) call mkkpg(kg1_k,kpg1_k,kpq(:),nkpg1,npw1_k)

 ! ===== Preparation of the non-local contributions

 dimffnlk=0;if (ipert<=natom) dimffnlk=1
 ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,ntypat))

 ! Compute nonlocal form factors ffnlk at (k+G)
 ! (only for atomic displacement perturbation)
 if (ipert<=natom) then
   ider=0;idir0=0
   call mkffnl(psps%dimekb,dimffnlk,psps%ekb,ffnlk,psps%ffspl,&
     gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat,&
     psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
 end if

!Compute nonlocal form factors ffnl1 at (k+q+G)
 !-- Atomic displacement perturbation
 if (ipert<=natom) then
   ider=0;idir0=0
 !-- k-point perturbation (1st-derivative)
 else if (ipert==natom+1) then
   ider=1;idir0=idir
 !-- k-point perturbation (2nd-derivative)
 else if (ipert==natom+10.or.ipert==natom+11) then
   ider=2;idir0=4
 !-- Electric field perturbation
 else if (ipert==natom+2) then
   if (psps%usepaw==1) then
     ider=1;idir0=idir
   else
     ider=0;idir0=0
   end if
 !-- Strain perturbation
 else if (ipert==natom+3.or.ipert==natom+4) then
   if (ipert==natom+3) istr=idir
   if (ipert==natom+4) istr=idir+3
   ider=1;idir0=-istr
 !-- Magnetic field perturbation ( SPr, Zeeman )
 else if(ipert==natom+5)then
   ider=0;idir0=0
 end if

 ! Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
 dimffnl1=1+ider
 if (ider==1.and.idir0==0) dimffnl1=2+2*psps%useylm
 if (ider==2.and.idir0==4) dimffnl1=3+7*psps%useylm
 ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,ntypat))

 call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
   psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
   npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

 ! Compute ffnl for nonlop with signs = 1
 print_info = 0
 if (dtset%prtvol==-19.or.dtset%prtvol==-20.or.dtset%prtvol==-21.or.dtset%nonlinear_info>=3) then
   print_info = 1
 end if
 if (present(ffnl1_test).and.print_info/=0.and.(ipert==natom+10.or.ipert==natom+11)) then
   ABI_ALLOCATE(ffnl1_test,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
   idir0 = 0 ! for nonlop with signs = 1
   call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1_test,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,ider,idir0,&
     psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
     npw1_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)
 end if

 !call cwtime_report(" getgh1c_setup_mkffnl", cpu, wall, gflops)

 !===== Preparation of the kinetic contributions

 ! Note that not all these arrays should be allocated in the general case when wtk_k vanishes

 ! Compute (1/2) (2 Pi)**2 (k+q+G)**2:
 ABI_ALLOCATE(kinpw1,(npw1_k))
 kinpw1(:)=zero
 call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k,0,0)

 ABI_ALLOCATE(dkinpw,(npw_k)) ! 1st derivative (1st direction)
 dkinpw(:)=zero
 if(ipert==natom+10 .and. idir>3) then
   ABI_ALLOCATE(dkinpw2,(npw_k)) ! 1st derivative (2nd directions)
   dkinpw2(:)=zero
 end if
 if(ipert==natom+10) then
   ABI_ALLOCATE(ddkinpw,(npw_k)) ! 2nd derivative
   ddkinpw(:)=zero
 end if

 ! -- k-point perturbation (1st-derivative)
 if (ipert==natom+1) then
   ! Compute the derivative of the kinetic operator vs k
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,dkinpw,kpoint,npw_k,idir,0) ! 1st derivative
 end if

 !-- k-point perturbation (2nd-derivative)
 if (ipert==natom+10.or.ipert==natom+11) then
   ! Compute the derivative of the kinetic operator vs k in kinpw, second and first orders
   if(ipert==natom+10 .and. idir<=3) then
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,dkinpw,kpoint,npw_k,idir,0) ! 1st derivative
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,ddkinpw,kpoint,npw_k,idir,idir) ! 2nd derivative
   else
     select case(idir)
     ! Diagonal terms:
     case(1)
       idir1 = 1
       idir2 = 1
     case(2)
       idir1 = 2
       idir2 = 2
     case(3)
       idir1 = 3
       idir2 = 3
     ! Upper triangular terms:
     case(4)
       idir1 = 2
       idir2 = 3
     case(5)
       idir1 = 1
       idir2 = 3
     case(6)
       idir1 = 1
       idir2 = 2
     ! Lower triangular terms:
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
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,dkinpw,kpoint,npw_k,idir1,0) !  1st derivative, idir1
     if(ipert==natom+10) then
       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,dkinpw2,kpoint,npw_k,idir2,0) ! 1st derivative, idir2
       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,ddkinpw,kpoint,npw_k,idir1,idir2) ! 2nd derivative
     end if
   end if
 end if

 !-- Strain perturbation
 if (ipert==natom+3.or.ipert==natom+4) then
   if (ipert==natom+3) istr=idir
   if (ipert==natom+4) istr=idir+3
   ! Compute the derivative of the kinetic operator vs strain
   call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,gprimd,istr,kg_k,kpoint,npw_k)
 end if

 !call cwtime_report(" getgh1c_setup_mkkin", cpu, wall, gflops)

 !===== Load the k/k+q dependent parts of the Hamiltonian
 ! Load k-dependent part in the Hamiltonian datastructure
 ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
 call gs_hamkq%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
                      ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
 if (size(ffnlk)>0) then
   call gs_hamkq%load_k(ffnl_k=ffnlk)
 else
   call gs_hamkq%load_k(ffnl_k=ffnl1)
 end if

  ! Load k+q-dependent part in the Hamiltonian datastructure
  ! Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
 call gs_hamkq%load_kprime(kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
   kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,compute_gbound=.true.)
 if (qne0) then
   ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
   call gs_hamkq%load_kprime(ph3d_kp=ph3d1,compute_ph3d=.true.)
 end if

 !Load k-dependent part in the 1st-order Hamiltonian datastructure
 call rf_hamkq%load_k(npw_k=npw_k,dkinpw_k=dkinpw)
 if (ipert==natom+10) then
   call rf_hamkq%load_k(ddkinpw_k=ddkinpw)
   if (idir>3) call rf_hamk_dir2%load_k(dkinpw_k=dkinpw2,ddkinpw_k=ddkinpw)
 end if

 call timab(195, 2, tsec)
 !call cwtime_report(" getgh1c_setup_hams", cpu, wall, gflops)

end subroutine getgh1c_setup
!!***

!!****f* ABINIT/getdc1
!!
!! NAME
!! getdc1
!!
!! FUNCTION
!! Compute |delta_C^(1)> from one wave function C - PAW ONLY
!! Compute <G|delta_C^(1)> and eventually <P_i| delta_C^(1)> (P_i= non-local projector)
!! delta_C^(1) is the variation of wavefunction only due to variation of overlap operator S.
!! delta_C^(1)=-1/2.Sum_j [ <C_j|S^(1)|C>.C_j
!!         see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42)
!!
!! INPUTS
!!  cgq(2,mcgq)=wavefunction coefficients for ALL bands at k+Q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors: cprjq=<P_i|Cnk+q>
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  istwfk=option parameter that describes the storage of wfs
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=information about MPI parallelization
!!  natom= number of atoms in cell
!!  nband=number of bands
!!  npw1=number of planewaves in basis sphere at k+Q
!!  nspinor=number of spinorial components of the wavefunctions
!!  opt_cprj=flag governing the computation of <P_i|delta_C^(1)> (P_i= non-local projector)
!!  s1cwave0(2,npw1*nspinor)=<G|S^(1)|C> where S^(1) is the first-order overlap operator
!!
!! OUTPUT
!!  dcwavef(2,npw1*nspinor)=change of wavefunction due to change of overlap PROJECTED ON PLANE-WAVES:
!!         dcwavef is delta_C(1)=-1/2.Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!!  === if optcprj=1 ===
!!  dcwaveprj(natom,nspinor*optcprj)=change of wavefunction due to change of overlap PROJECTED ON NL-PROJECTORS:
!!
!! PARENTS
!!      dfpt_cgwf,dfpt_nstpaw
!!
!! CHILDREN
!!      pawcprj_axpby,pawcprj_lincom,projbd
!!
!! SOURCE

subroutine getdc1(cgq,cprjq,dcwavef,dcwaveprj,ibgq,icgq,istwfk,mcgq,mcprjq,&
&                 mpi_enreg,natom,nband,npw1,nspinor,optcprj,s1cwave0)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ibgq,icgq,istwfk,mcgq,mcprjq,natom,nband,npw1,nspinor,optcprj
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cgq(2,mcgq),s1cwave0(2,npw1*nspinor)
 real(dp),intent(out) :: dcwavef(2,npw1*nspinor)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq)
 type(pawcprj_type),intent(inout) :: dcwaveprj(natom,nspinor*optcprj)

!Local variables-------------------------------
!scalars
 integer, parameter :: tim_projbd=0
 integer :: ipw
 real(dp),parameter :: scal=-half
!arrays
 real(dp), allocatable :: dummy(:,:),scprod(:,:)
 type(pawcprj_type),allocatable :: tmpcprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!$OMP PARALLEL DO
 do ipw=1,npw1*nspinor
   dcwavef(1:2,ipw)=s1cwave0(1:2,ipw)
 end do

 ABI_ALLOCATE(dummy,(0,0))
 ABI_ALLOCATE(scprod,(2,nband))

!=== 1- COMPUTE: <G|S^(1)|C_k> - Sum_j [<C_k+q,j|S^(1)|C_k>.<G|C_k+q,j>]
!!               using the projb routine
!Note the subtlety: projbd is called with useoverlap=0 and s1cwave0
!in order to get Sum[<cgq|s1|c>|cgq>]=Sum[<cgq|gs1>|cgq>]
 call projbd(cgq,dcwavef,-1,icgq,0,istwfk,mcgq,0,nband,npw1,nspinor,&
& dummy,scprod,0,tim_projbd,0,mpi_enreg%me_g0,mpi_enreg%comm_fft)

!=== 2- COMPUTE: <G|delta_C^(1)> = -1/2.Sum_j [<C_k+q,j|S^(1)|C_k>.<G|C_k+q,j>] by substraction
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(dcwavef,s1cwave0,npw1,nspinor)
 do ipw=1,npw1*nspinor
   dcwavef(1:2,ipw)=scal*(s1cwave0(1:2,ipw)-dcwavef(1:2,ipw))
 end do

!=== 3- COMPUTE: <P_i|delta_C^(1)> = -1/2.Sum_j [<C_k+q,j|S^(1)|C_k>.<P_i|C_k+q,j>]
 if (optcprj==1.and.mcprjq>0) then
   ABI_MALLOC(tmpcprj,(natom,nspinor))
   call pawcprj_lincom(scprod,cprjq(:,ibgq+1:ibgq+nspinor*nband),dcwaveprj,nband)
   call pawcprj_axpby(zero,scal,tmpcprj,dcwaveprj)
   ABI_FREE(tmpcprj)
 end if

 ABI_FREE(dummy)
 ABI_FREE(scprod)

 DBG_EXIT("COLL")

end subroutine getdc1
!!***

!!****f* ABINIT/getgh1dqc
!! NAME
!!  getgh1dqc
!!
!! FUNCTION
!! Computes <G|dH^(1)/dq_{gamma}|C> or <G|d^2H^(1)/dq_{gamma}dq_{delta}|C>
!! for input vector |C> expressed in reciprocal space.
!! dH^(1)/dq_{gamma} and d^2H^(1)/dq_{gamma}dq_{delta} are the first
!! and second q-gradient (at q=0) of the 1st-order perturbed Hamiltonian.
!! The first (second) derivative direction is inferred from idir (qdir1).
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (MR,MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cwave(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!     if not allocated or size=0, they are locally computed (and not sotred)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  idir=first index of the perturbation
!!  ipert=type of the perturbation
!!  mpi_enreg=information about MPI parallelization
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  optlocal=0: local part of H^(1) is not computed
!!           1: local part of H^(1) is computed in gvloc1dqc
!!  optnl=0: non-local part of H^(1) is not computed
!!        1: non-local part of H^(1) depending on VHxc^(1) is not computed in gvloc1dqc
!!        2: non-local part of H^(1) is totally computed in gvloc1dqc
!!  qdir1= direction of the 1st q-gradient
!!  rf_hamkq <type(rf_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,k+q
!!  qdir2= (optional) direction of the 2nd q-gradient
!!
!! OUTPUT
!!  gh1dqc(2,npw1*nspinor)= <G|dH^(1)/dq_{\gamma}|C> on the k+q sphere
!!  gvloc1dqc(2,npw1*nspinor)= local potential part of gh1dqc
!!  gvnl1dqc(2,npw1*nspinor)= non local potential part of gh1dqc
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  Currently two Hamiltonian gradients at (q=0) are implemented:
!!     ipert<=natom -> 		    first q-derivative along reduced coordinates directions
!!                     		    of the atomic displacement perturbation hamiltonian
!!     ipert==natom+3 or natom+4 -> second q-derivative along cartesian coordinates
!!                                  of the metric perturbation hamiltonian.
!! 				    Which is equivalent (except for an i factor) to the first
!!                                  q-derivative along cartesian coordinates of the strain
!!                                  perturbation hamiltonian.
!!
!! PARENTS
!!      dfpt_qdrpwf
!!
!! CHILDREN
!!      fourwf,nonlopdq
!!
!! SOURCE

subroutine getgh1dqc(cwave,cwaveprj,gh1dqc,gvloc1dqc,gvnl1dqc,gs_hamkq,&
&          idir,ipert,mpi_enreg,optlocal,optnl,qdir1,rf_hamkq,&
&          qdir2)                                                        !optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,optlocal,optnl,qdir1
 integer,intent(in),optional :: qdir2
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout),target :: rf_hamkq

!arrays
 real(dp),intent(inout) :: cwave(2,gs_hamkq%npw_k*gs_hamkq%nspinor)
 real(dp),intent(out) :: gh1dqc(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 real(dp),intent(out) :: gvloc1dqc(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 real(dp),intent(out) :: gvnl1dqc(2,gs_hamkq%npw_kp*gs_hamkq%nspinor)
 real(dp),pointer :: dqdqkinpw(:),kinpw1(:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,iidir,ipw,ipws,ispinor,my_nspinor,natom,nnlout
 integer :: npw,npw1,paw_opt,signs,tim_fourwf,tim_nonlop
 logical :: has_kin
 character(len=500) :: msg
 real(dp) :: lambda,weight

!arrays
 integer,parameter :: ngamma(3,3)=reshape((/1,6,5,9,2,4,8,7,3/),(/3,3/))
 real(dp) :: enlout(1),svectout_dum(1,1)
 real(dp),ABI_CONTIGUOUS pointer :: gvnl1dqc_(:,:)
 real(dp), allocatable :: work(:,:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!======================================================================
!== Initialisations and compatibility tests
!======================================================================

 npw  =gs_hamkq%npw_k
 npw1 =gs_hamkq%npw_kp
 natom=gs_hamkq%natom

!Compatibility tests
 if (mpi_enreg%paral_spinor==1) then
   msg='Not compatible with parallelization over spinorial components !'
   MSG_BUG(msg)
 end if

!Check sizes
 my_nspinor=max(1,gs_hamkq%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwave)<2*npw*my_nspinor) then
   msg='wrong size for cwave!'
   MSG_BUG(msg)
 end if
 if (size(gh1dqc)<2*npw1*my_nspinor) then
   msg='wrong size for gh1dqc!'
   MSG_BUG(msg)
 end if

!=============================================================================
!== Apply the q-gradients of the 1st-order local potential to the wavefunction
!=============================================================================

!Phonon and metric (strain) perturbation
 if (ipert<=natom+5.and.ipert/=natom+1.and.ipert/=natom+2.and.optlocal>0) then

   ABI_ALLOCATE(work,(2,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6))

   weight=one ; tim_fourwf=4
   call fourwf(rf_hamkq%cplex,rf_hamkq%vlocal1,cwave,gvloc1dqc,work,gs_hamkq%gbound_k,gs_hamkq%gbound_kp,&
 & gs_hamkq%istwf_k,gs_hamkq%kg_k,gs_hamkq%kg_kp,gs_hamkq%mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
 & npw,npw1,gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,2,tim_fourwf,weight,weight,&
 & use_gpu_cuda=gs_hamkq%use_gpu_cuda)

   ABI_FREE(work)

 else

!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gvloc1dqc(:,ipw)=zero
   end do

 end if

!================================================================================
!== Apply the q-gradients of the 1st-order non-local potential to the wavefunction
!================================================================================

!Initializations
lambda=zero
nnlout=1
tim_nonlop=0

!Allocations
ABI_ALLOCATE(gvnl1dqc_,(2,npw1*my_nspinor))

!Phonon perturbation
!-------------------------------------------
 !1st q-gradient
 if (ipert<=natom.and..not.present(qdir2).and.optnl>0) then
   cpopt=-1 ; choice=22 ; signs=2 ; paw_opt=0
   call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&  paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl1dqc_,iatom_only=ipert,qdir=qdir1)

!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gvnl1dqc(1,ipw)=gvnl1dqc_(1,ipw)
     gvnl1dqc(2,ipw)=gvnl1dqc_(2,ipw)
   end do

 !2nd q-gradient
 else if (ipert<=natom.and.present(qdir2).and.optnl>0) then
   iidir=ngamma(idir,qdir2)
   cpopt=-1 ; choice=25 ; signs=2 ; paw_opt=0
   call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,iidir,(/lambda/),mpi_enreg,1,nnlout,&
&  paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl1dqc_,iatom_only=ipert,qdir=qdir1)

!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gvnl1dqc(1,ipw)=gvnl1dqc_(1,ipw)
     gvnl1dqc(2,ipw)=gvnl1dqc_(2,ipw)
   end do

!Metric (strain) perturbation
!-------------------------------------------
 else if ((ipert==natom+3.or.ipert==natom+4).and.optnl>0) then
   cpopt=-1 ; choice=33 ; signs=2 ; paw_opt=0
   call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamkq,idir,(/lambda/),mpi_enreg,1,nnlout,&
&  paw_opt,signs,svectout_dum,tim_nonlop,cwave,gvnl1dqc_,qdir=qdir1)

!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gvnl1dqc(1,ipw)=gvnl1dqc_(1,ipw)
     gvnl1dqc(2,ipw)=gvnl1dqc_(2,ipw)
   end do

 else

!$OMP PARALLEL DO
   do ipw=1,npw1*my_nspinor
     gvnl1dqc(:,ipw)=zero
   end do

 end if

!==============================================================================
!== Apply the q-gradients of the 1st-order kinetic operator to the wavefunction
!== (add it to nl contribution)
!==============================================================================

!Strain (metric) perturbation
!-------------------------------------------
 has_kin=(ipert==natom+3.or.ipert==natom+4)
 if (associated(gs_hamkq%kinpw_kp)) then
   kinpw1 => gs_hamkq%kinpw_kp
 else if (has_kin) then
   msg='need kinpw1 allocated!'
   MSG_BUG(msg)
 end if
 if (associated(rf_hamkq%dkinpw_k)) then
   dqdqkinpw => rf_hamkq%dkinpw_k
 else if (has_kin) then
   msg='need dqdqkinpw allocated!'
   MSG_BUG(msg)
 end if

 if (has_kin) then
!  Remember that npw=npw1
   do ispinor=1,my_nspinor
!$OMP PARALLEL DO PRIVATE(ipw,ipws) SHARED(cwave,ispinor,gvnl1dqc,dqdqkinpw,kinpw1,npw,my_nspinor)
     do ipw=1,npw
       ipws=ipw+npw*(ispinor-1)
       if(kinpw1(ipw)<huge(zero)*1.d-11)then
         gvnl1dqc(1,ipws)=gvnl1dqc(1,ipws)+dqdqkinpw(ipw)*cwave(1,ipws)
         gvnl1dqc(2,ipws)=gvnl1dqc(2,ipws)+dqdqkinpw(ipw)*cwave(2,ipws)
       else
         gvnl1dqc(1,ipws)=zero
         gvnl1dqc(2,ipws)=zero
       end if
     end do
   end do
 end if

!===================================================================================
!== Sum contributions to get the application of dH^(1)/dq or d^2H^(1)/dqdq to the wf
!===================================================================================

 do ispinor=1,my_nspinor
   ipws=(ispinor-1)*npw1
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gh1dqc,gvnl1dqc,kinpw1,ipws,npw1)
   do ipw=1+ipws,npw1+ipws
     if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
       gh1dqc(1,ipw)=gvloc1dqc(1,ipw)+gvnl1dqc(1,ipw)
       gh1dqc(2,ipw)=gvloc1dqc(2,ipw)+gvnl1dqc(2,ipw)
     else
       gh1dqc(1,ipw)=zero
       gh1dqc(2,ipw)=zero
     end if
   end do
 end do

 ABI_FREE(gvnl1dqc_)
 DBG_EXIT("COLL")

end subroutine getgh1dqc
!!***

!!****f* m_hamiltonian/getgh1dqc_setup
!! NAME
!!  getgh1dqc_setup
!!
!! FUNCTION
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_qdrpwf
!!
!! CHILDREN
!!      kpgstr,load_k_hamiltonian,load_k_rf_hamiltonian,load_kprime_hamiltonian
!!      mkffnl,mkkin,mkkpg
!!
!! SOURCE

subroutine getgh1dqc_setup(gs_hamkq,rf_hamkq,dtset,psps,kpoint,kpq,idir,ipert,qdir1,&    ! In
&                natom,rmet,gprimd,gmet,istwf_k,npw_k,npw1_k,nylmgr,&                    ! In
&                useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                            ! In
&                nkpg,nkpg1,kpg_k,kpg1_k,dqdqkinpw,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&       ! Out
&                qdir2)                                                                  ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,istwf_k,natom,npw_k,npw1_k,nylmgr,qdir1,useylmgr1
 integer,intent(in),optional :: qdir2
 integer,intent(out) :: nkpg,nkpg1
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k),kg1_k(3,npw1_k)
 real(dp),intent(in) :: kpoint(3),kpq(3),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),intent(in) :: ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
! real(dp),intent(in) :: ylmgr1_k(npw1_k,3+6*((ipert-natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(in) :: ylmgr1_k(npw1_k,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(in) :: ylm1_k(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),allocatable,intent(out) :: dqdqkinpw(:),kinpw1(:)
 real(dp),allocatable,intent(out) :: ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable,intent(out) :: kpg_k(:,:),kpg1_k(:,:),ph3d(:,:,:),ph3d1(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: dimffnl1,dimffnlk,ider,idir0,ig,mu,mua,mub,ntypat
 integer :: nu,nua,nub
 logical :: qne0
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: ffnl1_tmp(:,:,:,:)


! *************************************************************************

 ntypat = psps%ntypat
 qne0=((kpq(1)-kpoint(1))**2+(kpq(2)-kpoint(2))**2+(kpq(3)-kpoint(3))**2>=tol14)

!Compute (k+G) vectors
 nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(3)
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
 end if

!Compute (k+q+G) vectors
 nkpg1=0;if(ipert>=1.and.ipert<=natom) nkpg1=3*dtset%nloalg(3)
 ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
 if (nkpg1>0) then
   call mkkpg(kg1_k,kpg1_k,kpq(:),nkpg1,npw1_k)
 end if

!===== Preparation of the non-local contributions

 dimffnlk=0;if (ipert<=natom) dimffnlk=1
 ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,ntypat))

!Compute nonlocal form factors ffnlk at (k+G)
if (ipert<=natom) then
 ider=0;idir0=0
 call mkffnl(psps%dimekb,dimffnlk,psps%ekb,ffnlk,psps%ffspl,&
& gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
& psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat,&
& psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
end if

!Compute nonlocal form factors ffnl1 at (k+q+G)
!TODO: For the second order gradients, this routine is called for each 3 directions of the
!derivative and every time it calculates all the form factors derivatives. This could be
!done just once.
 !-- 1st q-grad of atomic displacement perturbation
 if (ipert<=natom.and..not.present(qdir2)) then
   ider=1;idir0=qdir1
 !-- 2nd q-grad of atomic displacement perturbation
 else if (ipert<=natom.and.present(qdir2)) then
   ider=2;idir0=4
 !-- 2nd q-grad of metric (1st q-grad of strain) perturbation
 else if (ipert==natom+3.or.ipert==natom+4) then
   ider=2;idir0=0
 end if

!Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
 dimffnl1=1+ider
 if (ider==2.and.(idir0==0.or.idir0==4)) dimffnl1=3+7*psps%useylm
 ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,ntypat))
 call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
& psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
& npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

!Convert nonlocal form factors to cartesian coordinates.
!For metric (strain) perturbation only.
 if (ipert==natom+3.or.ipert==natom+4) then
   ABI_ALLOCATE(ffnl1_tmp,(npw1_k,dimffnl1,psps%lmnmax,ntypat))
   ffnl1_tmp=ffnl1

   !First q-derivative
   ffnl1(:,2:4,:,:)=zero
   do mu=1,3
     do ig=1,npw1_k
       do nu=1,3
         ffnl1(ig,1+mu,:,:)=ffnl1(ig,1+mu,:,:)+ffnl1_tmp(ig,1+nu,:,:)*gprimd(mu,nu)
       end do
     end do
   end do

   !Second q-derivative
   ffnl1(:,5:10,:,:)=zero
   do mu=1,6
     mua=alpha(mu);mub=beta(mu)
     do ig=1,npw1_k
       do nua=1,3
         do nub=1,3
           nu=gamma(nua,nub)
           ffnl1(ig,4+mu,:,:)=ffnl1(ig,4+mu,:,:)+ &
         & ffnl1_tmp(ig,4+nu,:,:)*gprimd(mua,nua)*gprimd(mub,nub)
         end do
       end do
     end do
   end do

   ABI_FREE(ffnl1_tmp)
 end if

!===== Preparation of the kinetic contributions
!Compute (1/2) (2 Pi)**2 (k+q+G)**2:
 ABI_ALLOCATE(kinpw1,(npw1_k))
 kinpw1(:)=zero
 call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k,0,0)

ABI_ALLOCATE(dqdqkinpw,(npw_k))
 !-- Metric (strain) perturbation
 if (ipert==natom+3.or.ipert==natom+4) then
   call mkkin_metdqdq(dqdqkinpw,dtset%effmass_free,gprimd,idir,kg_k,kpoint,npw_k,qdir1)
 else
   dqdqkinpw(:)=zero
 end if

!===== Load the k/k+q dependent parts of the Hamiltonian

!Load k-dependent part in the Hamiltonian datastructure
 ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
 call gs_hamkq%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
& ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
 if (size(ffnlk)>0) then
   call gs_hamkq%load_k(ffnl_k=ffnlk)
 else
   call gs_hamkq%load_k(ffnl_k=ffnl1)
 end if

!Load k+q-dependent part in the Hamiltonian datastructure
!    Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
 call gs_hamkq%load_kprime(kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
& kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
& compute_gbound=.true.)
 if (qne0) then
   ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
   call gs_hamkq%load_kprime(ph3d_kp=ph3d1,compute_ph3d=.true.)
 end if

!Load k-dependent part in the 1st-order Hamiltonian datastructure
 call rf_hamkq%load_k(npw_k=npw_k,dkinpw_k=dqdqkinpw)

end subroutine getgh1dqc_setup
!!***

end module m_getgh1c
!!***
