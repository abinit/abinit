!!****m* ABINIT/m_getchc
!! NAME
!!  m_getchc
!!
!! FUNCTION
!! Compute <G|H|C> for input vector |C> expressed in reciprocal space;
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, LSI, MT)
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

module m_getchc

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi

 use defs_abitypes, only : mpi_type
 use m_time,        only : timab
 use m_pawcprj,     only : pawcprj_type!, pawcprj_alloc, pawcprj_free, pawcprj_getdim
! use m_bandfft_kpt, only : bandfft_kpt, bandfft_kpt_get_ikpt
 use m_hamiltonian, only : gs_hamiltonian_type, KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME
 use m_nonlop,      only : nonlop
! use m_fock,        only : fock_common_type, fock_get_getghc_call
! use m_fock_getghc, only : fock_getghc, fock_ACE_getghc
 use m_cgtools,     only : dotprod_g

 implicit none

 private
!!***

 public :: getchc
 public :: getcsc
!!***

contains
!!***

!!****f* ABINIT/getchc
!!
!! NAME
!! getchc
!!
!! FUNCTION
!! Compute <C_left|H|C> for input vectors |C> and |C_left>.
!! Note that |C_left> can be an array of "ndat" wavefunctions if only the non-local part is computed
!!
!! INPUTS
!! cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!       (same meaning as in nonlop.F90 routine)
!!       if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!       if cpopt= 0, <p_lmn|in> are computed here and saved
!!       if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!       if cpopt= 2  <p_lmn|in> are already in memory;
!!       if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!       if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!! cwavef(2,npw*my_nspinor)=planewave coefficients of wavefunction.
!! cwavef_left(2,npw*my_nspinor)=planewave coefficients of wavefunction left.
!! cwaveprj(natom,my_nspinor*(1+cpopt))= wave function projected on nl projectors
!! cwaveprj_left(natom,my_nspinor*(1+cpopt))= wave function projected on nl projectors (for left WF)
!! cwavef_r(2,n4,n5,n6,nspinor) = wave function in real space
!! cwavef_left_r(2,n4,n5,n6,nspinor) = wave function in real space (for left WF)
!! gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!! lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!        Typically lambda is the eigenvalue (or its guess)
!! mpi_enreg=information about MPI parallelization
!! ndat=number of left wavefunctions
!! sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!    (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                      if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!! type_calc= option governing which part of Hamitonian is to be applied:
!             0: whole Hamiltonian
!!            1: local part only
!!            2: non-local+Fock+kinetic only (added to the existing Hamiltonian)
!!            3: local + kinetic only (added to the existing Hamiltonian)
!! ===== Optional inputs =====
!!   [kg_fft_k(3,:)]=optional, (k+G) vector coordinates to be used for the FFT tranformation
!!                   instead of the one contained in gs_ham datastructure.
!!                   Typically used for real WF (in parallel) which are FFT-transformed 2 by 2.
!!   [kg_fft_kp(3,:)]=optional, (k^prime+G) vector coordinates to be used for the FFT tranformation
!!   [select_k]=optional, option governing the choice of k points to be used.
!!             gs_ham datastructure contains quantities needed to apply Hamiltonian
!!             in reciprocal space between 2 kpoints, k and k^prime (equal in most cases);
!!             if select_k=1, <k^prime|H|k>       is applied [default]
!!             if select_k=2, <k|H|k^prime>       is applied
!!             if select_k=3, <k|H|k>             is applied
!!             if select_k=4, <k^prime|H|k^prime> is applied
!!
!! OUTPUT
!!   chc(2*ndat)=matrix elements <C_left|H|C> (if sij_opt>=0)
!!                           or <C_left|H-lambda.S|C> (if sij_opt=-1)
!! SIDE EFFECTS
!!
!! PARENTS
!!   cgwf_paw
!!
!! NOTE
!!   Case k_left/=k is not implemented yet
!!
!! CHILDREN
!!
!! SOURCE

subroutine getchc(chc,cpopt,cwavef,cwavef_left,cwaveprj,cwaveprj_left,cwavef_r,cwavef_left_r,&
&                 gs_ham,lambda,mpi_enreg,ndat,&
&                 sij_opt,type_calc,&
&                 kg_fft_k,kg_fft_kp,select_k) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat
 integer,intent(in) :: sij_opt,type_calc
 integer,intent(in),optional :: select_k
 real(dp),intent(in) :: lambda
 real(dp),intent(inout) :: chc(2*ndat)
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout),target :: gs_ham
!arrays
 integer,intent(in),optional,target :: kg_fft_k(:,:),kg_fft_kp(:,:)
 real(dp),intent(inout) :: cwavef(:,:),cwavef_left(:,:),cwavef_r(:,:,:,:,:),cwavef_left_r(:,:,:,:,:)
 type(pawcprj_type),intent(inout),target :: cwaveprj(:,:),cwaveprj_left(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: re=1,im=2
 integer :: choice,cpopt_here,i1,i2,i3,idat,idir
 integer :: ig,igspinor,ispinor,ispinor_left,my_nspinor
 integer :: nnlout,nffttot,npw,npw_k1,npw_k2,nspinortot,n1,n2,n3
 integer :: paw_opt,select_k_,shift1,shift2,signs,tim_nonlop
 logical :: k1_eq_k2,has_fock
 logical :: nspinor1TreatedByThisProc,nspinor2TreatedByThisProc
 character(len=500) :: msg
!arrays
 integer, pointer :: gbound_k1(:,:),gbound_k2(:,:),kg_k1(:,:),kg_k2(:,:)
 real(dp) :: enlout(ndat),enlout_im(ndat),lambda_ndat(ndat),tsec(2),z_tmp(2)
 real(dp),allocatable :: gsc(:,:),gvnlxc(:,:)
 real(dp), pointer :: kinpw_k1(:),kinpw_k2(:),kpt_k1(:),kpt_k2(:)

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in getchc:
 call timab(1370,1,tsec)

!Select k-dependent objects according to select_k input parameter
 select_k_=1;if (present(select_k)) select_k_=select_k
 if (select_k_==KPRIME_H_K) then
!  <k^prime|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_kp
 else if (select_k_==K_H_KPRIME) then
!  <k|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_k
 else if (select_k_==K_H_K) then
!  <k|H|k>
   npw_k1    =  gs_ham%npw_fft_k ; npw_k2    =  gs_ham%npw_fft_k
   kpt_k1    => gs_ham%kpt_k     ; kpt_k2    => gs_ham%kpt_k
   kg_k1     => gs_ham%kg_k      ; kg_k2     => gs_ham%kg_k
   gbound_k1 => gs_ham%gbound_k  ; gbound_k2 => gs_ham%gbound_k
   kinpw_k1  => gs_ham%kinpw_k   ; kinpw_k2  => gs_ham%kinpw_k
 else if (select_k_==KPRIME_H_KPRIME) then
!  <k^prime|H|k^prime>
   npw_k1    =  gs_ham%npw_fft_kp; npw_k2    =  gs_ham%npw_fft_kp
   kpt_k1    => gs_ham%kpt_kp    ; kpt_k2    => gs_ham%kpt_kp
   kg_k1     => gs_ham%kg_kp     ; kg_k2     => gs_ham%kg_kp
   gbound_k1 => gs_ham%gbound_kp ; gbound_k2 => gs_ham%gbound_kp
   kinpw_k1  => gs_ham%kinpw_kp  ; kinpw_k2  => gs_ham%kinpw_kp
 end if
 k1_eq_k2=(all(abs(kpt_k1(:)-kpt_k2(:))<tol8))
 if (.not.k1_eq_k2) then
   MSG_ERROR('getchc is not implemented yet for k1/=k2')
 end if

!Check sizes
 my_nspinor=max(1,gs_ham%nspinor/mpi_enreg%nproc_spinor)
 if (size(cwavef)<2*npw_k1*my_nspinor) then
   msg='wrong size for cwavef!'
   MSG_BUG(msg)
 end if
 if (gs_ham%usepaw==1.and.cpopt>=0) then
   if (size(cwaveprj)<gs_ham%natom*my_nspinor) then
     msg='wrong size for cwaveprj!'
     MSG_BUG(msg)
   end if
 end if
 if (gs_ham%usepaw==1) then
   if (size(cwaveprj_left)<gs_ham%natom*my_nspinor*ndat) then
     msg='wrong size for cwaveprj_left!'
     MSG_BUG(msg)
   end if
 end if

!Eventually overwrite plane waves data for FFT
 if (present(kg_fft_k)) then
   kg_k1 => kg_fft_k ; kg_k2 => kg_fft_k
   npw_k1=size(kg_k1,2) ; npw_k2=size(kg_k2,2)
 end if
 if (present(kg_fft_kp)) then
   kg_k2 => kg_fft_kp ; npw_k2=size(kg_k2,2)
 end if

!paral_kgb constraint
 if (mpi_enreg%paral_kgb==1.and.(.not.k1_eq_k2)) then
   msg='paral_kgb=1 not allowed for k/=k_^prime!'
   MSG_BUG(msg)
 end if

!Do we add Fock exchange term ?
 has_fock=(associated(gs_ham%fockcommon))
 if (has_fock) then
   MSG_BUG('Fock not implemented yet')
 end if
! if (has_fock) fock => gs_ham%fockcommon

!Parallelization over spinors management
 if (mpi_enreg%paral_spinor==0) then
   shift1=npw_k1;shift2=npw_k2
   nspinor1TreatedByThisProc=.true.
   nspinor2TreatedByThisProc=(nspinortot==2)
 else
   shift1=0;shift2=0
   nspinor1TreatedByThisProc=(mpi_enreg%me_spinor==0)
   nspinor2TreatedByThisProc=(mpi_enreg%me_spinor==1)
 end if

 npw=gs_ham%npw_k
 nspinortot=gs_ham%nspinor

!============================================================
! Application of the local potential
!============================================================

 if ((type_calc==0).or.(type_calc==1).or.(type_calc==3)) then

   call timab(1371,1,tsec)
!  Need a Vlocal
   if (.not.associated(gs_ham%vlocal)) then
     MSG_BUG("We need vlocal in gs_ham!")
   end if
   if (ndat>1) then
     MSG_ERROR("ndat should be 1 for the local part")
   end if

   n1=gs_ham%ngfft(1)
   n2=gs_ham%ngfft(2)
   n3=gs_ham%ngfft(3)
   nffttot=n1*n2*n3
   chc = zero
!  Treat scalar local potentials
   if (gs_ham%nvloc==1) then
     if (gs_ham%istwf_k==2) then
       do i3=1,n3
         do i2=1,n2
           do i1=1,n1
             chc(1) = chc(1) + gs_ham%vlocal(i1,i2,i3,1)*cwavef_r(1,i1,i2,i3,1)*cwavef_left_r(1,i1,i2,i3,1)
           end do
         end do
       end do
       chc(2)=zero
     else
       do ispinor=1,my_nspinor
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               z_tmp(1) = cwavef_r(1,i1,i2,i3,ispinor)*cwavef_left_r(1,i1,i2,i3,ispinor) &
&                 +cwavef_r(2,i1,i2,i3,ispinor)*cwavef_left_r(2,i1,i2,i3,ispinor)
               z_tmp(2) = cwavef_r(2,i1,i2,i3,ispinor)*cwavef_left_r(1,i1,i2,i3,ispinor) &
&                 -cwavef_r(1,i1,i2,i3,ispinor)*cwavef_left_r(2,i1,i2,i3,ispinor)
               chc(1) = chc(1) + gs_ham%vlocal(i1,i2,i3,1)*z_tmp(1)
               chc(2) = chc(2) + gs_ham%vlocal(i1,i2,i3,1)*z_tmp(2)
             end do
           end do
         end do
       end do
     end if
   else ! nvloc = 4
     do ispinor=1,my_nspinor
       do ispinor_left=1,my_nspinor
         do i3=1,n3
           do i2=1,n2
             do i1=1,n1
               z_tmp(1) = cwavef_r(1,i1,i2,i3,ispinor)*cwavef_left_r(1,i1,i2,i3,ispinor_left) &
&                        +cwavef_r(2,i1,i2,i3,ispinor)*cwavef_left_r(2,i1,i2,i3,ispinor_left)
               z_tmp(2) = cwavef_r(2,i1,i2,i3,ispinor)*cwavef_left_r(1,i1,i2,i3,ispinor_left) &
&                        -cwavef_r(1,i1,i2,i3,ispinor)*cwavef_left_r(2,i1,i2,i3,ispinor_left)
               if (ispinor==ispinor_left) then
                 ! Then vloc is real : vloc_uu = vloc(1) and vloc_dd = vloc(2)
                 chc(1) = chc(1) + gs_ham%vlocal(i1,i2,i3,ispinor)*z_tmp(1)
                 chc(2) = chc(2) + gs_ham%vlocal(i1,i2,i3,ispinor)*z_tmp(2)
               else if (ispinor==1.and.ispinor_left==2) then ! Psi(left)_d^* Psi_u vloc_ud
                 ! Otherwise vloc is complex Re(vloc_ud) = vloc(3)
                 !                           Im(vloc_ud) = vloc(4)
                 !                               vloc_du = (vloc_ud)^*
                 chc(1) = chc(1) + gs_ham%vlocal(i1,i2,i3,3)*z_tmp(1) + gs_ham%vlocal(i1,i2,i3,4)*z_tmp(2)
                 chc(2) = chc(2) + gs_ham%vlocal(i1,i2,i3,3)*z_tmp(1) - gs_ham%vlocal(i1,i2,i3,4)*z_tmp(2)
               else if (ispinor==2.and.ispinor_left==1) then ! Psi(left)_u^* Psi_d vloc_ud
                 chc(1) = chc(1) + gs_ham%vlocal(i1,i2,i3,3)*z_tmp(1) - gs_ham%vlocal(i1,i2,i3,4)*z_tmp(2)
                 chc(2) = chc(2) + gs_ham%vlocal(i1,i2,i3,3)*z_tmp(1) + gs_ham%vlocal(i1,i2,i3,4)*z_tmp(2)
               end if
             end do
           end do
         end do
       end do
     end do
   end if
   chc = chc / dble(nffttot)

   call timab(1371,2,tsec)

 end if ! type_calc

 if ((type_calc==0).or.(type_calc==2).or.(type_calc==3).or.(type_calc==4)) then

!============================================================
! Application of the non-local potential and the Fock potential
!============================================================

  if ((type_calc==0).or.(type_calc==2).or.(type_calc==4)) then

     signs=1 ; choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=15
     cpopt_here=-1;if (gs_ham%usepaw==1) cpopt_here=cpopt
     ABI_ALLOCATE(gvnlxc,(0,0))
     ABI_ALLOCATE(gsc,(0,0))
!     if (has_fock) then
!       if (gs_ham%usepaw==1) then
!         cpopt_here=max(cpopt,0)
!         if (cpopt<2) then
!           ABI_DATATYPE_ALLOCATE(cwaveprj_fock,(gs_ham%natom,my_nspinor*ndat))
!           ABI_ALLOCATE(dimcprj,(gs_ham%natom))
!           call pawcprj_getdim(dimcprj,gs_ham%natom,gs_ham%nattyp,gs_ham%ntypat,&
!&           gs_ham%typat,fock%pawtab,'O')
!           call pawcprj_alloc(cwaveprj_fock,0,dimcprj)
!           ABI_DEALLOCATE(dimcprj)
!         else
!           cwaveprj_fock=>cwaveprj
!         end if
!         cwaveprj_nonlop=>cwaveprj_fock
!       else
!         cwaveprj_nonlop=>cwaveprj
!         cwaveprj_fock=>cwaveprj
!       end if
!     else
!     end if
     paw_opt=gs_ham%usepaw ; if (sij_opt/=0) paw_opt=sij_opt+3
     lambda_ndat = lambda

     enlout=zero
     enlout_im=zero
     call nonlop(choice,cpopt_here,cwaveprj,enlout,gs_ham,idir,lambda_ndat,mpi_enreg,1,&
&     nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_,&
&     cprjin_left=cwaveprj_left,enlout_im=enlout_im,ndat_left=ndat)

     do idat=1,ndat
       chc(2*idat-1) = chc(2*idat-1) + enlout(idat)
       chc(2*idat  ) = chc(2*idat  ) + enlout_im(idat)
     end do

     ABI_DEALLOCATE(gvnlxc)
     ABI_DEALLOCATE(gsc)

   end if ! if(type_calc...

!============================================================
! Assemble kinetic, local, nonlocal and Fock contributions
!============================================================

   if (type_calc==0.or.type_calc==2.or.type_calc==3) then

     if (ndat>1) then
       MSG_ERROR("ndat should be 1 for the kinetic part")
     end if

     call timab(1372,1,tsec)
!    Add modified kinetic contributions
     !  to <CP|H|C(n,k)>.
     if (gs_ham%istwf_k==1) then
!      !!$OMP PARALLEL DO PRIVATE(igspinor) COLLAPSE(2)
       do ispinor=1,my_nspinor
         do ig=1,npw_k2
           igspinor=ig+npw_k2*(ispinor-1)
           if (kinpw_k2(ig)<huge(zero)*1.d-11) then
             chc(1) = chc(1) +  kinpw_k2(ig)*cwavef(re,igspinor)*cwavef_left(re,igspinor)
             chc(1) = chc(1) +  kinpw_k2(ig)*cwavef(im,igspinor)*cwavef_left(im,igspinor)
             chc(2) = chc(2) +  kinpw_k2(ig)*cwavef(im,igspinor)*cwavef_left(re,igspinor)
             chc(2) = chc(2) -  kinpw_k2(ig)*cwavef(re,igspinor)*cwavef_left(im,igspinor)
           end if
         end do ! ig
       end do ! ispinor
     else if (gs_ham%istwf_k==2.and.mpi_enreg%me_g0==1) then
       if (kinpw_k2(1)<huge(zero)*1.d-11) then
         chc(1) = chc(1) + kinpw_k2(1)*cwavef(re,1)*cwavef_left(re,1)
       end if
       do ig=2,npw_k2
         if (kinpw_k2(ig)<huge(zero)*1.d-11) then
           chc(1) = chc(1) + 2*kinpw_k2(ig)*cwavef(re,ig)*cwavef_left(re,ig)
           chc(1) = chc(1) + 2*kinpw_k2(ig)*cwavef(im,ig)*cwavef_left(im,ig)
         end if
       end do
     else
       do ig=1,npw_k2
         if (kinpw_k2(ig)<huge(zero)*1.d-11) then
           chc(1) = chc(1) + 2*kinpw_k2(ig)*cwavef(re,ig)*cwavef_left(re,ig)
           chc(1) = chc(1) + 2*kinpw_k2(ig)*cwavef(im,ig)*cwavef_left(im,ig)
         end if
       end do ! ig
     end if
!    Special case of PAW + Fock : only return Fock operator contribution in gvnlxc
!     if (gs_ham%usepaw==1 .and. has_fock)then
!       gvnlxc=gvnlxc-gvnlc
!       ABI_DEALLOCATE(gvnlc)
!     endif
!
!     if ((type_calc==0).or.(type_calc==2)) then
!       if (has_fock.and.gs_ham%usepaw==1.and.cpopt<2) then
!         call pawcprj_free(cwaveprj_fock)
!         ABI_DATATYPE_DEALLOCATE(cwaveprj_fock)
!       end if
!     end if
     call timab(1372,2,tsec)
   end if

 end if ! type_calc

 call timab(1370,2,tsec)

 DBG_EXIT("COLL")

end subroutine getchc
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/getcsc
!! NAME
!! getcsc
!!
!! FUNCTION
!! Compute <C_left|S|C> for input vectors |C> and |C_left>.
!! Note that |C_left> can be an array of "ndat" wavefunctions
!!
!! INPUTS
!!  cpopt=flag defining the status of cwaveprj%cp(:)=<Proj_i|Cnk> scalars (PAW only)
!!        (same meaning as in nonlop.F90 routine)
!!        if cpopt=-1, <p_lmn|in> (and derivatives) are computed here (and not saved)
!!        if cpopt= 0, <p_lmn|in> are computed here and saved
!!        if cpopt= 1, <p_lmn|in> and first derivatives are computed here and saved
!!        if cpopt= 2  <p_lmn|in> are already in memory;
!!        if cpopt= 3  <p_lmn|in> are already in memory; first derivatives are computed here and saved
!!        if cpopt= 4  <p_lmn|in> and first derivatives are already in memory;
!!  cwavef(2,npw*my_nspinor)=planewave coefficients of wavefunction.
!!  cwavef_left(2,npw*my_nspinor)=planewave coefficients of wavefunction left.
!!  cprj(natom,my_nspinor*(1+cpopt))= wave function projected on nl projectors
!!  cprj_left(natom,my_nspinor*(1+cpopt))= wave function projected on nl projectors (for left WF)
!!  gs_ham <type(gs_hamiltonian_type)>=all data for the Hamiltonian to be applied
!!  lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!         Typically lambda is the eigenvalue (or its guess)
!!  mpi_enreg=information about MPI parallelization
!!  ndat=number of left wavefunctions
!!  mpi_enreg=information about MPI parallelization
!!
!! OUTPUT
!!   csc(2*ndat)=matrix elements <C_left|S|C>
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine getcsc(csc,cpopt,cwavef,cwavef_left,cprj,cprj_left,gs_ham,mpi_enreg,ndat,&
&                 tim_getcsc,&
&                 select_k) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpopt,ndat,tim_getcsc
 integer,intent(in),optional :: select_k
 real(dp),intent(out) :: csc(2*ndat)
 type(MPI_type),intent(in) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_ham
!arrays
 real(dp),intent(inout) :: cwavef(:,:)
 real(dp),intent(inout),target :: cwavef_left(:,:)
 type(pawcprj_type),intent(inout) :: cprj(:,:)
 type(pawcprj_type),intent(inout),target :: cprj_left(:,:)

!Local variables-------------------------------
!scalars
 integer :: choice,idat,idir,istwf_k
 integer :: npw,nspinor,paw_opt,select_k_,signs,tim_nonlop,nnlout
 !character(len=500) :: msg
!arrays
 real(dp) :: tsec(2),real_csc(ndat)
 real(dp),allocatable :: gsc(:,:),gvnlxc(:,:)
 real(dp),allocatable :: enlout(:),enlout_im(:)
! real(dp) :: dum
! real(dp),pointer :: cwavef_left_idat(:,:)
! *********************************************************************

 DBG_ENTER("COLL")

 call timab(1360+tim_getcsc,1,tsec)

 istwf_k = gs_ham%istwf_k
 npw = gs_ham%npw_k
 nspinor = gs_ham%nspinor

 if (size(cwavef,2)/=npw*nspinor) then
   MSG_BUG('Wrong size for cwavef')
 end if
 if (size(cwavef_left,2)/=npw*nspinor*ndat) then
   MSG_BUG('Wrong size for cwavef_left')
 end if
 if (size(cprj,2)/=nspinor) then
   MSG_BUG('Wrong size for cprj')
 end if
 if (size(cprj_left,2)/=nspinor*ndat) then
   MSG_BUG('Wrong size for cprj_left')
 end if

 call timab(1361,1,tsec)
 if (istwf_k==1) then
   call zgemv('C',npw*nspinor,ndat,cone,cwavef_left,npw*nspinor,cwavef,1,czero,csc,1)
 else ! nspinor==1 in that case
   call dgemv('C',2*npw,ndat,one,cwavef_left,2*npw,cwavef,1,zero,real_csc,1)
   do idat=1,ndat
     csc(2*idat-1) = two*real_csc(idat)
     csc(2*idat  ) = zero
   end do
   if (istwf_k==2 .and. mpi_enreg%me_g0==1) then ! Gamma k-point and I have G=0
     do idat=1,ndat
       csc(2*idat-1) = csc(2*idat-1) - cwavef_left(1,1+npw*(idat-1))*cwavef(1,1)
     end do
   end if
 end if
 call timab(1361,2,tsec)


 if (gs_ham%usepaw==1) then
   select_k_=1;if (present(select_k)) select_k_=select_k
   choice=1 ; nnlout=1 ; idir=0 ; tim_nonlop=16 ; paw_opt=3
   ABI_ALLOCATE(gsc,(0,0))
   ABI_ALLOCATE(gvnlxc,(0,0))
   ABI_ALLOCATE(enlout   ,(ndat))
   ABI_ALLOCATE(enlout_im,(ndat))
   enlout=zero
   enlout_im=zero
   signs=1
   call nonlop(choice,cpopt,cprj,enlout,gs_ham,idir,(/zero/),mpi_enreg,1,&
&   nnlout,paw_opt,signs,gsc,tim_nonlop,cwavef,gvnlxc,select_k=select_k_,&
&   cprjin_left=cprj_left,enlout_im=enlout_im,ndat_left=ndat)
   do idat=1,ndat
     csc(2*idat-1) = csc(2*idat-1) + enlout(idat)
     csc(2*idat  ) = csc(2*idat  ) + enlout_im(idat)
   end do
   ABI_DEALLOCATE(gsc)
   ABI_DEALLOCATE(gvnlxc)
   ABI_DEALLOCATE(enlout   )
   ABI_DEALLOCATE(enlout_im)
 end if

 call timab(1360+tim_getcsc,2,tsec)

 DBG_EXIT("COLL")

end subroutine getcsc
!!***

end module m_getchc
!!***
