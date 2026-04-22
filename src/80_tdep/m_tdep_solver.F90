!!****m* ABINIT/m_tdep_solver
!! NAME
!!  m_tdep_solver
!!
!! FUNCTION
!!  This module contains the TDEP Solver data type
!!  which fits the IFC on the forces and displacements of the sampling.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2026 ABINIT group (GA,FB,JB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_solver

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools
 use m_numeric_tools
 use m_linalg_interfaces
 use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
 use m_tdep_shell,       only : Shell_type
 use m_tdep_sym,         only : Symmetries_type
 use m_tdep_sampling,    only : tdep_Sampling_type
 use m_tdep_model,       only : tdep_Model_type
 use m_tdep_constraints, only : Constraints_type, tdep_calc_orthonorm

 implicit none

 type tdep_Solver_type

   integer :: order
   ! Maximum IFC order (2, 3, or 4)

   integer :: natom
   ! Number of atoms in the supercell

   integer :: natom_unitcell
   ! Number of atoms in the unitcell

   integer :: my_nstep
   ! Number of MD steps

   integer :: ntotcoeff
   ! Total number of coefficients

   integer :: ntotconst
   ! Total number of constraints

   integer :: ncoeff1st
   ! Number of 1st order coefficients

   integer :: ncoeff2nd
   ! Number of 2nd order coefficients

   integer :: ncoeff3rd
   ! Number of 3rd order coefficients

   integer :: ncoeff4th
   ! Number of 4th order coefficients

   integer :: nconst_1st
   ! Number of 1st order constraints

   integer :: nconst_2nd
   ! Number of 2nd order constraints

   integer :: nconst_3rd
   ! Number of 3rd order constraints

   integer :: nconst_4th
   ! Number of 4th order constraints

   integer :: nconst_rot2nd
   integer :: nconst_huang
   integer :: nconst_dynmat
   integer :: nconst_rot3rd
   integer :: nconst_asr3rd
   integer :: nconst_rot4th
   integer :: nconst_asr4th

   double precision, allocatable :: fcoeff(:,:)
   ! fcoeff(3*natom*my_nstep,ntotcoeff)
   ! All the cartesian displacement matrices
   ! at every order (u, u*u, u*u*u, u*u*u*u).

   double precision, allocatable :: const(:,:)
   ! const(ntotconst,ntotcoeff)
   ! The constraint matrices.

   double precision, allocatable :: theta(:)
   ! theta(ntotcoeff)
   ! The IFC coefficients at all orders, as a flat array. 

   double precision, allocatable :: Forces(:)
   ! Forces(3*natom*my_nstep)
   ! The cartesian forces for all configurations, as a flat array.
   ! These are weighted by the number of configurations.
   ! In case of dipole-dipole interaction, the long-range part of the
   ! forces should be removed.

 end type tdep_Solver_type

 public :: tdep_solver_init
 public :: tdep_solver_free
 public :: tdep_solver_set_residual_forces
 public :: tdep_calc_phi1fcoeff
 public :: tdep_calc_phi2fcoeff
 public :: tdep_calc_phi3fcoeff
 public :: tdep_calc_phi4fcoeff
 public :: tdep_calc_MoorePenrose
 public :: tdep_calc_constraints

contains

!=====================================================================================================

 subroutine tdep_solver_init(Solver, Invar, Shell1at, Shell2at, Shell3at, Shell4at)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell1at, Shell2at, Shell3at, Shell4at

  Solver%order = Invar%order
  Solver%natom = Invar%natom
  Solver%natom_unitcell = Invar%natom_unitcell
  Solver%my_nstep = Invar%my_nstep

  !Rotational invariances (1st order)
  !    constraints = 3
  Solver%nconst_1st = 3**2
  
  !Rotational invariances (2nd order) + Symetry of the Dynamical Matrix + Huang invariances
  !    constraints = natom*3**2 + (3*natom_unitcell)**2 + 3**4
  Solver%nconst_rot2nd = 3**3*Solver%natom_unitcell
  Solver%nconst_dynmat = (3*Solver%natom_unitcell)**2
  Solver%nconst_huang  = 3**4
  Solver%nconst_2nd = Solver%nconst_rot2nd + Solver%nconst_dynmat + Solver%nconst_huang
  
  !Rotational invariances (3rd order) + acoustic sum rules (3rd order)
  Solver%nconst_3rd=0
  if (Solver%order.ge.3) then
  !    constraints = natom_unitcell*natom*3**4 + natom_unitcell*natom*3**3
    Solver%nconst_rot3rd = 3**4 * Solver%natom_unitcell * Solver%natom
    Solver%nconst_asr3rd = 3**3 * Solver%natom_unitcell * Solver%natom
    Solver%nconst_3rd = Solver%nconst_rot3rd + Solver%nconst_asr3rd
  end if
  
  !Rotational invariances (4th order) + acoustic sum rules (4th order)
  Solver%nconst_4th=0
  if (Solver%order.ge.4) then
  !    constraints = natom_unitcell*natom**2*3**5 + natom_unitcell*natom**2*3**4
  !FB   Solver%nconst_rot4th = 3**5*natom_unitcell*natom**2
  !FB4TH   Solver%nconst_asr4th = 3**4*natom_unitcell*natom**2
    Solver%nconst_rot4th = 0
    Solver%nconst_asr4th = 0
    Solver%nconst_4th = Solver%nconst_rot4th + Solver%nconst_asr4th
  end if
 
  Solver%ncoeff1st = Shell1at%ntotcoeff
  Solver%ncoeff2nd = Shell2at%ntotcoeff
  Solver%ncoeff3rd = 0
  Solver%ncoeff4th = 0
  if (Solver%order.ge.3) Solver%ncoeff3rd = Shell3at%ntotcoeff
  if (Solver%order.ge.4) Solver%ncoeff4th = Shell4at%ntotcoeff

  Solver%ntotcoeff = Solver%ncoeff1st  + Solver%ncoeff2nd  + Solver%ncoeff3rd  + Solver%ncoeff4th
  Solver%ntotconst = Solver%nconst_1st + Solver%nconst_2nd + Solver%nconst_3rd + Solver%nconst_4th

  ABI_CALLOC(Solver%fcoeff, (3*Solver%natom*Solver%my_nstep,Solver%ntotcoeff))
  ABI_CALLOC(Solver%const, (Solver%ntotconst,Solver%ntotcoeff))
  ABI_CALLOC(Solver%Forces, (3*Solver%natom*Solver%my_nstep))
  ABI_CALLOC(Solver%theta, (Solver%ntotcoeff))

 end subroutine tdep_solver_init

!=====================================================================================================

 subroutine tdep_solver_free(Solver)

  type(tdep_Solver_type), intent(inout) :: Solver

  ABI_FREE(Solver%fcoeff)
  ABI_FREE(Solver%const)
  ABI_FREE(Solver%Forces)
  ABI_FREE(Solver%theta)

 end subroutine tdep_solver_free

!=====================================================================================================

!!****f* ABINIT/m_tdep_solver/tdep_solver_set_residual_forces
!! NAME
!!  tdep_solver_set_residual_forces
!!
!! FUNCTION
!! Compute the residual forces into the solver, that is, the differences
!! between the MD forces and the ones from the TDEP model.
!! The model forces may be zero, but they are non-zero when we remove
!! the long-range part of the IFC forces, or when we solve each IFC order
!! successively, and we remove the forces from the previously computed orders.
!!
!! INPUTS
!!  Solver = TDEP Solver object that will compute the coefficients.
!!  MD = TDEP Sampling object containing the input positions and forces.
!!  Model = TDEP Model object containing the IFC and corresponding forces.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  The following quantities in MD are computed:
!!  
!! NOTES
!!
!! SOURCE

 subroutine tdep_solver_set_residual_forces(Solver, MD, Model)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(tdep_Sampling_type), intent(in) :: MD
  type(tdep_Model_type), intent(in) :: Model

  integer :: ii,jj,istep,iatom

  do istep=1,MD%my_nstep
   do iatom=1,MD%natom
     do ii=1,3
      jj = ii + 3*(iatom-1) + 3*MD%natom*(istep-1)
      Solver%Forces(jj) = (MD%Forces(jj) - Model%Forces(jj)) * MD%weights(istep)
     end do
   end do
  end do

 end subroutine tdep_solver_set_residual_forces

!====================================================================================================

subroutine tdep_calc_phi1fcoeff(Solver,Invar,Shell1at,Sym)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell1at
  type(Symmetries_type),intent(in) :: Sym

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,iatshell,iat_mod
  integer :: icoeff,isym,mu,iatref
  double precision :: terme

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '############## Fill the matrices used in the pseudo-inverse #################'
  write(Invar%stdout,*) '#############################################################################'

  write(Invar%stdout,*) ' Compute the coefficients (at the 1st order) used in the Moore-Penrose...'
  do ishell=1,Shell1at%nshell
    if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
    do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
      iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell)
      iat_mod=mod(iatom+Invar%natom_unitcell-1,Invar%natom_unitcell)+1
      if (iat_mod==1) cycle
      iatref=Shell1at%iatref(ishell)
      isym=Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
      ncoeff     =Shell1at%ncoeff(ishell)
      ncoeff_prev=Shell1at%ncoeff_prev(ishell)
      do mu=1,3
        do icoeff=1,ncoeff
          terme=sum(Sym%S_ref(mu,:,isym,1)*Shell1at%proj(:,icoeff,ishell))
          do istep=1,Invar%my_nstep
            Solver%fcoeff(mu+3*(iatom-1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)= &
&           Solver%fcoeff(mu+3*(iatom-1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)+terme
!           Add all the other contributions, when iat_mod==1 (due to ASR)
            Solver%fcoeff(mu+3*(iatom-iat_mod+1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)= &
&           Solver%fcoeff(mu+3*(iatom-iat_mod+1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)-terme
          end do !istep
        end do
      end do
    end do !iatshell
  end do !ishell
  write(Invar%stdout,*) ' ------- achieved'

end subroutine tdep_calc_phi1fcoeff

!====================================================================================================

subroutine tdep_calc_phi2fcoeff(Solver,Invar,Shell2at,Sym,MD)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell2at
  type(Symmetries_type),intent(in) :: Sym
  type(tdep_Sampling_type), intent(in) :: MD

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,iatshell
  integer :: icoeff,isym
  integer :: mu,nu,alpha,beta,itrans
  double precision :: terme,temp
  double precision :: udiff(3),SSu(3,9)
  double precision, allocatable :: SS_ref(:,:,:,:,:)

! For each couple of atoms, transform the Phi2 (3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  ABI_MALLOC(SS_ref,(3,9,3,Sym%nsym,2)); SS_ref(:,:,:,:,:)=zero
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)
            SS_ref(mu,beta+(alpha-1)*3,nu,isym,1)=temp
            SS_ref(mu,alpha+(beta-1)*3,nu,isym,2)=temp
          end do
        end do
      end do
    end do
  end do

  write(Invar%stdout,*) ' Compute the coefficients (at the 2nd order) used in the Moore-Penrose...'
  do ishell=1,Shell2at%nshell
    do iatom=1,Invar%natom
      if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        if (iatom==jatom) cycle
        isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell2at%ncoeff(ishell)
        ncoeff_prev=Shell2at%ncoeff_prev(ishell)+Solver%ncoeff1st

        do istep=1,Invar%my_nstep
!         In order to impose the acoustic sum rule we use (u(j)-u(i))==u_j^\nu
          udiff(1)=(MD%ucart(1,jatom,istep)-MD%ucart(1,iatom,istep))*MD%weights(istep)
          udiff(2)=(MD%ucart(2,jatom,istep)-MD%ucart(2,iatom,istep))*MD%weights(istep)
          udiff(3)=(MD%ucart(3,jatom,istep)-MD%ucart(3,iatom,istep))*MD%weights(istep)

!         F_i^\mu(t)=\sum_{\alpha\beta,j,\nu}S^{\mu\alpha}.S^{\nu\beta}.\Phi_{ij}^{\alpha\beta}.u_j^\nu(t)
          SSu(:,:)=zero
          do nu=1,3
            SSu(:,:)=SSu(:,:)+SS_ref(:,:,nu,isym,itrans)*udiff(nu)
          end do
          do mu=1,3
            do icoeff=1,ncoeff
              terme=sum(SSu(mu,:)*Shell2at%proj(:,icoeff,ishell))
!FB              write(Invar%stdlog,*) 'indices=', mu+3*(iatom-1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev
              Solver%fcoeff(mu+3*(iatom-1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)= &
&             Solver%fcoeff(mu+3*(iatom-1)+3*Invar%natom*(istep-1),icoeff+ncoeff_prev)+terme
            end do
          end do

        end do !istep
      end do !iatshell
    end do !iatom
  end do !ishell
  write(Invar%stdout,*) ' ------- achieved'
  ABI_FREE(SS_ref)

end subroutine tdep_calc_phi2fcoeff

!====================================================================================================

subroutine tdep_calc_phi3fcoeff(Solver,Invar,Shell3at,Sym,MD)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell3at
  type(Symmetries_type),intent(in) :: Sym
  type(tdep_Sampling_type), intent(in) :: MD

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom
  integer :: icoeff,isym,itrans,iatshell
  integer :: mu,nu,xi,alpha,beta,gama,iindex
  double precision :: temp
  double precision :: udiff_ki(3),udiff_ji(3)
  double precision, allocatable :: SSS_proj(:,:,:,:)
  double precision :: SSS_tmp(27), proj_tmp(27)
  type(Constraints_type) :: Const

  ABI_MALLOC(Const%Sprod,(Sym%nsym,6))
  do isym=1,Sym%nsym
    do itrans=1,6
      ABI_MALLOC(Const%Sprod(isym,itrans)%SSS,(3,27,3,3)); Const%Sprod(isym,itrans)%SSS(:,:,:,:)=zero
    end do
  end do

! For each couple of atoms, transform the Phi3 (3x3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            do xi=1,3
              do gama=1,3
                temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)*Sym%S_ref(xi,gama,isym,1)
                Const%Sprod(isym,1)%SSS(mu,gama+(beta-1)*3+(alpha-1)*9,nu,xi)=temp !\Phi3_efg
                Const%Sprod(isym,2)%SSS(mu,gama+(beta-1)*3+(alpha-1)*9,xi,nu)=temp !\Phi3_egf
                Const%Sprod(isym,3)%SSS(nu,gama+(beta-1)*3+(alpha-1)*9,mu,xi)=temp !\Phi3_feg
                Const%Sprod(isym,4)%SSS(nu,gama+(beta-1)*3+(alpha-1)*9,xi,mu)=temp !\Phi3_fge
                Const%Sprod(isym,5)%SSS(xi,gama+(beta-1)*3+(alpha-1)*9,mu,nu)=temp !\Phi3_gef
                Const%Sprod(isym,6)%SSS(xi,gama+(beta-1)*3+(alpha-1)*9,nu,mu)=temp !\Phi3_gfe
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  write(Invar%stdout,*) ' Compute the coefficients (at the 3rd order) used in the Moore-Penrose...'
  do ishell=1,Shell3at%nshell
    do iatom=1,Invar%natom
      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
!FB        if (iatom==jatom.or.iatom==katom) cycle
        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell3at%ncoeff(ishell)
        ncoeff_prev=Shell3at%ncoeff_prev(ishell)+Solver%ncoeff2nd+Solver%ncoeff1st

        ABI_MALLOC(SSS_proj,(3,3,3,ncoeff)) ; SSS_proj(:,:,:,:)=zero
        do mu=1,3
          do nu=1,3
            do xi=1,3
              SSS_tmp(:)=Const%Sprod(isym,itrans)%SSS(mu,:,nu,xi)
              do icoeff=1,ncoeff
                proj_tmp(:)=Shell3at%proj(:,icoeff,ishell)
                SSS_proj(mu,nu,xi,icoeff)=DDOT(27,SSS_tmp,1,proj_tmp,1)
              end do
            end do
          end do
        end do
        do istep=1,Invar%my_nstep
          iindex=3*(iatom-1)+3*Invar%natom*(istep-1)
!         In order to impose the acoustic sum rule we use :
!FB          udiff_ji(:)=MD%ucart(:,jatom,istep)-MD%ucart(:,iatom,istep)
!FB          udiff_ki(:)=MD%ucart(:,katom,istep)-MD%ucart(:,iatom,istep)
          udiff_ji(:)=MD%ucart(:,jatom,istep)
          udiff_ki(:)=MD%ucart(:,katom,istep)
!         F_i^{\mu}(t)=\sum_{\alpha\beta\gamma,jk,\nu\xi} S^{\mu\alpha}.S^{\nu\beta}.S^{\xi\gamma}.
!                      \Phi3_{ijk}^{\alpha\beta\gamma}.udiff_k^\xi(t).udiff_j^\nu(t)
          do nu=1,3
            do xi=1,3
              Solver%fcoeff(iindex+1:iindex+3,ncoeff_prev+1:ncoeff_prev+ncoeff)= &
&             Solver%fcoeff(iindex+1:iindex+3,ncoeff_prev+1:ncoeff_prev+ncoeff)+&
&             SSS_proj(1:3,nu,xi,1:ncoeff)*udiff_ji(nu)*udiff_ki(xi)/2.d0*MD%weights(istep)
            end do
          end do
        end do !istep
        ABI_FREE(SSS_proj)
      end do !iatshell
    end do !iatom
  end do !ishell
  write(Invar%stdout,*) ' ------- achieved'
  do isym=1,Sym%nsym
    do itrans=1,6
      ABI_FREE(Const%Sprod(isym,itrans)%SSS)
    end do
  end do
  ABI_FREE(Const%Sprod)

end subroutine tdep_calc_phi3fcoeff

!====================================================================================================

subroutine tdep_calc_phi4fcoeff(Solver,Invar,Shell4at,Sym,MD)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell4at
  type(Symmetries_type),intent(in) :: Sym
  type(tdep_Sampling_type), intent(inout) :: MD

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom,latom
  integer :: icoeff,isym,iatshell,itrans,counter
  integer :: mu,nu,xi,zeta,alpha,beta,gama,delta,iindex_l,iindex_h
  integer :: ncoeff_prev_l,ncoeff_prev_h
  double precision :: temp,SSSS_tmp(81),proj_tmp(81)
  double precision, allocatable :: SSSS_proj(:,:,:,:,:)
  type(Constraints_type) :: Const

  ABI_MALLOC(Const%Sprod,(Sym%nsym,24))
  do isym=1,Sym%nsym
    do itrans=1,24
      ABI_MALLOC(Const%Sprod(isym,itrans)%SSSS,(3,81,3,3,3)); Const%Sprod(isym,itrans)%SSSS(:,:,:,:,:)=zero
    end do
  end do

! For each couple of atoms, transform the Phi4 (3x3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            do xi=1,3
              do gama=1,3
                do zeta=1,3
                  do delta=1,3
#if defined FC_NVHPC
                    if (itrans == -1) write(std_out, *)"NVHPC freezes here that is fixed by this print statement."
#endif

                    counter=delta+(gama-1)*3+(beta-1)*9+(alpha-1)*27
                    temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu  ,beta ,isym,1)*&
&                        Sym%S_ref(xi,gama ,isym,1)*Sym%S_ref(zeta,delta,isym,1)
                    Const%Sprod(isym,1 )%SSSS(mu,counter,nu,xi,zeta)=temp !\Phi4_efgh
                    Const%Sprod(isym,2 )%SSSS(mu,counter,xi,nu,zeta)=temp !\Phi4_egfh
                    Const%Sprod(isym,3 )%SSSS(nu,counter,mu,xi,zeta)=temp !\Phi4_fegh
                    Const%Sprod(isym,4 )%SSSS(nu,counter,xi,mu,zeta)=temp !\Phi4_fgeh
                    Const%Sprod(isym,5 )%SSSS(xi,counter,mu,nu,zeta)=temp !\Phi4_gefh
                    Const%Sprod(isym,6 )%SSSS(xi,counter,nu,mu,zeta)=temp !\Phi4_gfeh

                    Const%Sprod(isym,7 )%SSSS(mu,counter,nu,zeta,xi)=temp !\Phi4_efhg
                    Const%Sprod(isym,8 )%SSSS(mu,counter,xi,zeta,nu)=temp !\Phi4_eghf
                    Const%Sprod(isym,9 )%SSSS(nu,counter,mu,zeta,xi)=temp !\Phi4_fehg
                    Const%Sprod(isym,10)%SSSS(nu,counter,xi,zeta,mu)=temp !\Phi4_fghe
                    Const%Sprod(isym,11)%SSSS(xi,counter,mu,zeta,nu)=temp !\Phi4_gehf
                    Const%Sprod(isym,12)%SSSS(xi,counter,nu,zeta,mu)=temp !\Phi4_gfhe

                    Const%Sprod(isym,13)%SSSS(mu,counter,zeta,nu,xi)=temp !\Phi4_ehfg
                    Const%Sprod(isym,14)%SSSS(mu,counter,zeta,xi,nu)=temp !\Phi4_ehgf
                    Const%Sprod(isym,15)%SSSS(nu,counter,zeta,mu,xi)=temp !\Phi4_fheg
                    Const%Sprod(isym,16)%SSSS(nu,counter,zeta,xi,mu)=temp !\Phi4_fhge
                    Const%Sprod(isym,17)%SSSS(xi,counter,zeta,mu,nu)=temp !\Phi4_ghef
                    Const%Sprod(isym,18)%SSSS(xi,counter,zeta,nu,mu)=temp !\Phi4_ghfe

                    Const%Sprod(isym,19)%SSSS(zeta,counter,mu,nu,xi)=temp !\Phi4_hefg
                    Const%Sprod(isym,20)%SSSS(zeta,counter,mu,xi,nu)=temp !\Phi4_hegf
                    Const%Sprod(isym,21)%SSSS(zeta,counter,nu,mu,xi)=temp !\Phi4_hfeg
                    Const%Sprod(isym,22)%SSSS(zeta,counter,nu,xi,mu)=temp !\Phi4_hfge
                    Const%Sprod(isym,23)%SSSS(zeta,counter,xi,mu,nu)=temp !\Phi4_hgef
                    Const%Sprod(isym,24)%SSSS(zeta,counter,xi,nu,mu)=temp !\Phi4_hgfe

                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  write(Invar%stdout,*) ' Compute the coefficients (at the 4th order) used in the Moore-Penrose...'
  do ishell=1,Shell4at%nshell
    do iatom=1,Invar%natom
      if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
        latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
        isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell4at%ncoeff(ishell)
        ncoeff_prev=Shell4at%ncoeff_prev(ishell)+Solver%ncoeff3rd+Solver%ncoeff2nd+Solver%ncoeff1st
        ncoeff_prev_l=ncoeff_prev+1
        ncoeff_prev_h=ncoeff_prev+ncoeff
#if defined FC_NVHPC
        if (itrans == -1) write(std_out, *)"NVHPC freezes here that is fixed by this print statement."
#endif
        ABI_MALLOC(SSSS_proj,(3,3,3,3,ncoeff)) ; SSSS_proj(:,:,:,:,:)=zero
        do mu=1,3
          do nu=1,3
            do xi=1,3
              do zeta=1,3
                SSSS_tmp(:)=Const%Sprod(isym,itrans)%SSSS(mu,:,nu,xi,zeta)
                do icoeff=1,ncoeff
                  proj_tmp(:)=Shell4at%proj(:,icoeff,ishell)
!                 SSSS_proj(mu,nu,xi,zeta,icoeff)=DDOT(81,Const%Sprod(isym,itrans)%SSSS(mu,:,nu,xi,zeta),1,proj(:,icoeff,ishell),1)
                  SSSS_proj(mu,nu,xi,zeta,icoeff)=DDOT(81,SSSS_tmp,1,proj_tmp,1)
                end do
              end do
            end do
          end do
        end do
        do istep=1,Invar%my_nstep
          iindex_l=3*(iatom-1)+3*Invar%natom*(istep-1)+1
          iindex_h=3*(iatom-1)+3*Invar%natom*(istep-1)+3
!         F_i^{\mu}(t)=\sum_{\alpha\beta\gamma\delta,jkl,\nu\xi\zeta} S^{\mu\alpha}.S^{\nu\beta}.S^{\xi\gamma}.S^{\zeta\delta}.
!                      \Phi4_{ijkl}^{\alpha\beta\gamma\delta}.u_l^\zeta(t).u_k^\xi(t).u_j^\nu(t)
          do nu=1,3
            do xi=1,3
              do zeta=1,3
                Solver%fcoeff(iindex_l:iindex_h,ncoeff_prev_l:ncoeff_prev_h)= &
&               Solver%fcoeff(iindex_l:iindex_h,ncoeff_prev_l:ncoeff_prev_h)+&
&               SSSS_proj(1:3,nu,xi,zeta,1:ncoeff)*MD%ucart(nu,jatom,istep)*MD%ucart(xi,katom,istep)*MD%ucart(zeta,latom,istep)/6.d0 *&
&               MD%weights(istep)
              end do
            end do
          end do
        end do !istep
        ABI_FREE(SSSS_proj)
      end do !iatshell
    end do !iatom
  end do !ishell
  write(Invar%stdout,*) ' ------- achieved'
  do isym=1,Sym%nsym
    do itrans=1,24
      ABI_FREE(Const%Sprod(isym,itrans)%SSSS)
    end do
  end do
  ABI_FREE(Const%Sprod)

end subroutine tdep_calc_phi4fcoeff

!=====================================================================================================

 subroutine tdep_calc_MoorePenrose(Solver,simult,Invar,MPIdata)

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata
  integer, intent(in) :: simult

  integer :: INFO,ntotcoeff,ntotconst
  integer :: natnstep,nconcoef,ierr,ncoeff_prev,nconst_prev,iconst,icoeff
  integer, allocatable :: IPIV(:)
  double precision, allocatable :: WORK(:)
  double precision, allocatable :: ffcoeff_tmp(:,:),fforces_tmp(:),b_const(:)
  double precision, allocatable :: A_tot(:,:),A_inv(:,:),b_tot(:),x_tot(:)

  write(Invar%stdout,*) '################### And compute the pseudo-inverse ##########################'
  write(Invar%stdout,*) '#############################################################################'

  natnstep = 3 * Solver%natom * Solver%my_nstep

  if (simult.eq.0) then
!   Simultaneously (Invar%together=1)
    ncoeff_prev=0
    nconst_prev=0
    ntotcoeff=Solver%ntotcoeff
    ntotconst=Solver%ntotconst
  else if (simult.eq.1) then
!   Successively (Invar%together=0 and Invar%order=2)
    ncoeff_prev=0
    nconst_prev=0
    ntotcoeff=Solver%ncoeff1st +Solver%ncoeff2nd
    ntotconst=Solver%nconst_1st+Solver%nconst_2nd
  else if (simult.eq.2) then
!   Successively (Invar%together=0 and Invar%order=3)
    ncoeff_prev=Solver%ncoeff1st +Solver%ncoeff2nd
    nconst_prev=Solver%nconst_1st+Solver%nconst_2nd
    ntotcoeff=Solver%ncoeff3rd
    ntotconst=Solver%nconst_3rd
  else if (simult.eq.3) then
!   Successively (Invar%together=0 and Invar%order=4)
    ncoeff_prev=Solver%ncoeff1st +Solver%ncoeff2nd +Solver%ncoeff3rd
    nconst_prev=Solver%nconst_1st+Solver%nconst_2nd+Solver%nconst_3rd
    ntotcoeff=Solver%ncoeff4th
    ntotconst=Solver%nconst_4th
  end if
  nconcoef=ntotcoeff+ntotconst
  if ((ntotconst.gt.0).and.(simult.ge.2)) then
    ABI_MALLOC(b_const,(ntotconst)) ; b_const(:)=0.d0
    do iconst=1,ntotconst
      do icoeff=1,ncoeff_prev
        b_const(iconst)=b_const(iconst)+&
&         Solver%const(nconst_prev+iconst,icoeff)*Solver%theta(icoeff)
      end do
    end do
  end if

  ABI_CALLOC(ffcoeff_tmp,(ntotcoeff,ntotcoeff))
  ABI_CALLOC(fforces_tmp,(ntotcoeff))
  ABI_CALLOC(A_tot,(nconcoef,nconcoef))
  ABI_CALLOC(A_inv,(nconcoef,nconcoef))
  ABI_CALLOC(b_tot,(nconcoef))
  ABI_CALLOC(x_tot,(nconcoef))
  call DGEMM('T','N',ntotcoeff,ntotcoeff,natnstep,2.d0,&
&            Solver%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            Solver%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            0.d0,ffcoeff_tmp,ntotcoeff)
! NOTE, we have to solve F_ij = -\sum_j \Phi_ij u_j, so we add a minus sign
  call DGEMV('T',natnstep,ntotcoeff,-2.d0,&
&            Solver%fcoeff(:,ncoeff_prev+1:ncoeff_prev+ntotcoeff),natnstep,&
&            Solver%Forces,1,0.d0,fforces_tmp,1)
  call xmpi_sum(ffcoeff_tmp,MPIdata%comm_step,ierr)
  call xmpi_sum(fforces_tmp,MPIdata%comm_step,ierr)

  A_tot(1:ntotcoeff,1:ntotcoeff)=ffcoeff_tmp(1:ntotcoeff,1:ntotcoeff)
  ABI_FREE(ffcoeff_tmp)
  if (ntotconst.gt.0) then
    A_tot(ntotcoeff+1:nconcoef,1:ntotcoeff)=&
&                     Solver%const(nconst_prev+1:nconst_prev+ntotconst,ncoeff_prev+1:ncoeff_prev+ntotcoeff)
    A_tot(1:ntotcoeff,ntotcoeff+1:nconcoef)=&
&                     transpose(Solver%const(nconst_prev+1:nconst_prev+ntotconst,ncoeff_prev+1:ncoeff_prev+ntotcoeff))
!FB    ABI_FREE(Solver%const)
  end if
  b_tot(1:ntotcoeff)=fforces_tmp(:)
  if ((ntotconst.gt.0).and.(simult.ge.2)) then
    b_tot(ntotcoeff+1:nconcoef)=-b_const(1:ntotconst)
    ABI_FREE(b_const)
  end if
  ABI_FREE(fforces_tmp)

  ABI_MALLOC(WORK, (5 * nconcoef)); WORK(:) = 0.d0
  ABI_MALLOC(IPIV, (nconcoef)); IPIV(:) = 0
  A_inv(:,:) = A_tot(:,:)
  !BEGIN DEBUG
  !write(Invar%stdout,*) ' '
  !write(Invar%stdout,*) ' The matrix A_inv is (before DGETRF):'
  !do icoeff=1,nconcoef
  !  write(Invar%stdout,*) (A_inv(icoeff,iconst), iconst=1, nconcoef)
  !end do
  !END DEBUG

  ! Check for small pivot elements
  do icoeff=1,nconcoef
    if (abs(A_inv(icoeff, icoeff)) < tol12) then
      write(Invar%stdlog,*) ' WARNING: Small pivot value at index ', icoeff, ' : ', A_inv(icoeff, icoeff)
!      A_inv(icoeff, icoeff) = tol14 ! Regularization to avoid numerical issues
      A_inv(icoeff, icoeff) = max(EPSILON(1.0_dp) * maxval(abs(A_inv)), tol12)
    end if
  end do

  ! Perform LU factorization
  call DGETRF(nconcoef, nconcoef, A_inv, nconcoef, IPIV, INFO)
  if (INFO.ne.0) then
    write(Invar%stdout,*) 'ERROR: Singular matrix detected in DGETRF. INFO=', INFO
    stop
  end if

  ! Check for small pivot elements
  do icoeff=1,nconcoef
    if (abs(A_inv(icoeff, icoeff)) < tol12) then
      write(Invar%stdlog,*) ' WARNING: Small pivot value at index ', icoeff, ' : ', A_inv(icoeff, icoeff)
!      A_inv(icoeff, icoeff) = tol14 ! Regularization to avoid numerical issues
      A_inv(icoeff, icoeff) = max(EPSILON(1.0_dp) * maxval(abs(A_inv)), tol12)
    end if
  end do

  ! Compute matrix inverse using LU decomposition
  call DGETRI(nconcoef, A_inv, nconcoef, IPIV, WORK, 5 * nconcoef, INFO)
  if (INFO.ne.0) then
    write(Invar%stdout,*) 'ERROR: Matrix inversion failed in DGETRI. INFO=', INFO
    stop
  end if
  ! BEGIN DEBUG
  !write(Invar%stdout,*) ' '
  !write(Invar%stdout,*) ' The inverse matrix is (after DGETRI):'
  !do icoeff=1,nconcoef
  !  write(Invar%stdout,*) (A_inv(icoeff,iconst), iconst=1, nconcoef)
  !end do
  ! END DEBUG

  ABI_FREE(WORK)
  ABI_FREE(IPIV)

  call DGEMV('N',nconcoef,nconcoef,1.d0,A_inv,nconcoef,b_tot,1,0.d0,x_tot,1)
  write(Invar%stdout,*) ' The problem is solved'
  write(Invar%stdout,*) ' '
  !BEGIN DEBUG
  !write(Invar%stdout,*) ' The solutions are:'
  !do icoeff=1,nconcoef
  !  write(Invar%stdout,'(1x,i4,1x,f15.10)') icoeff,x_tot(icoeff)
  !end do
  !write(Invar%stdout,'(a,1x,f15.10)')'  condition number=',maxval(x_tot(:))/minval(x_tot(:))
  !END DEBUG

  Solver%theta(ncoeff_prev+1:ncoeff_prev+ntotcoeff)=x_tot(1:ntotcoeff)
  ABI_FREE(A_tot)
  ABI_FREE(A_inv)
  ABI_FREE(b_tot)
  ABI_FREE(x_tot)

 end subroutine tdep_calc_MoorePenrose


!====================================================================================================

subroutine tdep_calc_constraints(Solver,distance,Invar,MPIdata,Sym,&
&                                Shell1at,Shell2at,Shell3at,Shell4at) 

  type(tdep_Solver_type), intent(inout) :: Solver
  type(atdep_dataset_type),intent(in) :: Invar
  type(Symmetries_type),intent(in) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Shell_type),intent(in) :: Shell1at
  type(Shell_type),intent(in) :: Shell2at
  type(Shell_type),optional,intent(in) :: Shell3at
  type(Shell_type),optional,intent(in) :: Shell4at
  double precision, intent(in) :: distance(Invar%natom,Invar%natom,4)

  integer :: ishell,ncoeff,ncoeff_prev,iatom,jatom,katom,latom,iatshell,counter
  integer :: icoeff,iconst,nconst_loc,iconst_loc,iconst_new,isym,itrans,ntotcoeff,iat_mod
  integer :: mu,nu,xi,zeta,alpha,beta,gama,delta,lambda,natom_unitcell,natom,ii
  double precision :: terme,temp,terme1,terme2,terme3,terme4
  double precision, allocatable :: SS_ref(:,:,:,:,:)
  double precision, allocatable :: vect(:,:)
  double precision, allocatable :: const_rot1st(:,:,:)
  double precision, allocatable :: const_rot2nd(:,:,:,:,:)
  double precision, allocatable :: const_dynmat(:,:,:,:,:)
  double precision, allocatable :: const_huang(:,:,:,:,:)
!FB  double precision, allocatable :: const_asr4th(:,:,:,:,:,:)
!FB  double precision, allocatable :: const_rot4th(:,:,:,:,:,:,:)
  type(Constraints_type) :: Const3,Const4
  logical :: order2,order3,order4

  !TODO Move parts of this routine into m_tdep_constraints

  natom_unitcell=Invar%natom_unitcell
  natom         =Invar%natom

  order2 = .false.
  order3 = .false.
  order4 = .false.
  if (Invar%order.ge.2) order2=.true.
  if (Invar%order.ge.3) order3=.true.
!FB4th  if (Invar%order.ge.4) order4=.true.
  if (Invar%order.ge.4) order4=.false.

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '###################### Compute the constraints ##############################'

! For each couple of atoms, transform the Phi2 (3x3) ifc matrix using the symetry operation (S)
  if (order2.or.order3) then
    ABI_CALLOC(SS_ref,(3,9,3,Sym%nsym,2))
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta,nu,isym,1)
              SS_ref(alpha,nu+(mu-1)*3,beta,isym,1)=temp
              SS_ref(alpha,mu+(nu-1)*3,beta,isym,2)=temp
            end do
          end do
        end do
      end do
    end do
  end if
! For each couple of atoms, transform the Phi3 (3x3x3) ifc matrix using the symetry operation (S)
  if (order3.or.order4) then
    ABI_MALLOC(Const3%Sprod,(Sym%nsym,6))
    do isym=1,Sym%nsym
      do itrans=1,6
        ABI_CALLOC(Const3%Sprod(isym,itrans)%SSS,(3,27,3,3))
      end do  
    end do  
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              do gama=1,3
                do xi=1,3
                  temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta,nu,isym,1)*Sym%S_ref(gama,xi,isym,1)
                  Const3%Sprod(isym,1)%SSS(alpha,xi+(nu-1)*3+(mu-1)*9,beta ,gama) =temp !\Phi3_efg
                  Const3%Sprod(isym,2)%SSS(alpha,xi+(nu-1)*3+(mu-1)*9,gama ,beta) =temp !\Phi3_egf
                  Const3%Sprod(isym,3)%SSS(beta ,xi+(nu-1)*3+(mu-1)*9,alpha,gama) =temp !\Phi3_feg
                  Const3%Sprod(isym,4)%SSS(beta ,xi+(nu-1)*3+(mu-1)*9,gama ,alpha)=temp !\Phi3_fge
                  Const3%Sprod(isym,5)%SSS(gama ,xi+(nu-1)*3+(mu-1)*9,alpha,beta) =temp !\Phi3_gef
                  Const3%Sprod(isym,6)%SSS(gama ,xi+(nu-1)*3+(mu-1)*9,beta ,alpha)=temp !\Phi3_gfe
                end do
              end do
            end do
            end do
        end do
      end do
    end do
  end if

! For each couple of atoms, transform the Phi4 (3x3x3x3) ifc matrix using the symetry operation (S)
  if (order4) then
    ABI_MALLOC(Const4%Sprod,(Sym%nsym,24))
    do isym=1,Sym%nsym
      do itrans=1,24
        ABI_CALLOC(Const4%Sprod(isym,itrans)%SSSS,(3,81,3,3,3))
      end do  
    end do  
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              do gama=1,3
                do xi=1,3
                  do delta=1,3
                    do zeta=1,3
                      counter=zeta+(xi-1)*3+(nu-1)*9+(mu-1)*27
                      temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta  ,nu ,isym,1)*&
&                          Sym%S_ref(gama,xi ,isym,1)*Sym%S_ref(delta,zeta,isym,1)
                      Const4%Sprod(isym,1 )%SSSS(alpha,counter,beta,gama,delta)=temp !\Phi4_efgh
                      Const4%Sprod(isym,2 )%SSSS(alpha,counter,gama,beta,delta)=temp !\Phi4_egfh
                      Const4%Sprod(isym,3 )%SSSS(beta,counter,alpha,gama,delta)=temp !\Phi4_fegh
                      Const4%Sprod(isym,4 )%SSSS(beta,counter,gama,alpha,delta)=temp !\Phi4_fgeh
                      Const4%Sprod(isym,5 )%SSSS(gama,counter,alpha,beta,delta)=temp !\Phi4_gefh
                      Const4%Sprod(isym,6 )%SSSS(gama,counter,beta,alpha,delta)=temp !\Phi4_gfeh
  
                      Const4%Sprod(isym,7 )%SSSS(alpha,counter,beta,delta,gama)=temp !\Phi4_efhg
                      Const4%Sprod(isym,8 )%SSSS(alpha,counter,gama,delta,beta)=temp !\Phi4_eghf
                      Const4%Sprod(isym,9 )%SSSS(beta,counter,alpha,delta,gama)=temp !\Phi4_fehg
                      Const4%Sprod(isym,10)%SSSS(beta,counter,gama,delta,alpha)=temp !\Phi4_fghe
                      Const4%Sprod(isym,11)%SSSS(gama,counter,alpha,delta,beta)=temp !\Phi4_gehf
                      Const4%Sprod(isym,12)%SSSS(gama,counter,beta,delta,alpha)=temp !\Phi4_gfhe
  
                      Const4%Sprod(isym,13)%SSSS(alpha,counter,delta,beta,gama)=temp !\Phi4_ehfg
                      Const4%Sprod(isym,14)%SSSS(alpha,counter,delta,gama,beta)=temp !\Phi4_ehgf
                      Const4%Sprod(isym,15)%SSSS(beta,counter,delta,alpha,gama)=temp !\Phi4_fheg
                      Const4%Sprod(isym,16)%SSSS(beta,counter,delta,gama,alpha)=temp !\Phi4_fhge
                      Const4%Sprod(isym,17)%SSSS(gama,counter,delta,alpha,beta)=temp !\Phi4_ghef
                      Const4%Sprod(isym,18)%SSSS(gama,counter,delta,beta,alpha)=temp !\Phi4_ghfe
  
                      Const4%Sprod(isym,19)%SSSS(delta,counter,alpha,beta,gama)=temp !\Phi4_hefg
                      Const4%Sprod(isym,20)%SSSS(delta,counter,alpha,gama,beta)=temp !\Phi4_hegf
                      Const4%Sprod(isym,21)%SSSS(delta,counter,beta,alpha,gama)=temp !\Phi4_hfeg
                      Const4%Sprod(isym,22)%SSSS(delta,counter,beta,gama,alpha)=temp !\Phi4_hfge
                      Const4%Sprod(isym,23)%SSSS(delta,counter,gama,alpha,beta)=temp !\Phi4_hgef
                      Const4%Sprod(isym,24)%SSSS(delta,counter,gama,beta,alpha)=temp !\Phi4_hgfe
  
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do  
  end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute the constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ntotcoeff=Solver%ntotcoeff
! First order only
  write(Invar%stdout,*) '########################## At the 1st order #################################'
  if (order2) then
    ABI_CALLOC(const_rot1st, (3,3,ntotcoeff))
    ABI_CALLOC(const_rot2nd, (3,3,3,natom_unitcell,ntotcoeff))
    ABI_CALLOC(const_dynmat, (3,3,natom_unitcell,natom_unitcell,ntotcoeff))
    ABI_CALLOC(const_huang, (3,3,3,3,ntotcoeff))
    do ishell=1,Shell1at%nshell
      if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
        iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell) 
        if (iatom.ge.natom_unitcell) cycle
        if (iatom.eq.1) cycle
        isym=Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
        ncoeff     =Shell1at%ncoeff(ishell)
        ncoeff_prev=Shell1at%ncoeff_prev(ishell)
        do alpha=1,3
          do beta=1,3
            do icoeff=1,ncoeff
!             1/ Rotational invariances (1st order)
              terme1=sum(Sym%S_ref(alpha,:,isym,1)*Shell1at%proj(:,icoeff,ishell))*distance(1,iatom,beta +1)
              terme2=sum(Sym%S_ref(beta ,:,isym,1)*Shell1at%proj(:,icoeff,ishell))*distance(1,iatom,alpha+1)
              const_rot1st(alpha,beta,icoeff+ncoeff_prev)= &
&             const_rot1st(alpha,beta,icoeff+ncoeff_prev)+terme1-terme2
  
!             2/ Rotational invariances (for the 2nd order)
              do gama=1,3
                terme1=zero ; terme2=zero 
                if (alpha.eq.gama) terme1=sum(Sym%S_ref(beta,:,isym,1)*Shell1at%proj(:,icoeff,ishell))
                if (alpha.eq.beta) terme2=sum(Sym%S_ref(gama,:,isym,1)*Shell1at%proj(:,icoeff,ishell))
                const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)=&
&               const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)+terme1-terme2
                const_rot2nd(alpha,beta,gama,1,icoeff+ncoeff_prev)=&
&               const_rot2nd(alpha,beta,gama,1,icoeff+ncoeff_prev)-terme1+terme2
              end do
            end do    
          end do    
        end do  
      end do !iatshell
    end do !ishell

!   First + second order
    write(Invar%stdout,*) '########################## At the 2nd order #################################'
    do ishell=1,Shell2at%nshell
      do iatom=1,natom_unitcell
        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          if (iatom==jatom) cycle
          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell2at%ncoeff(ishell)
          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+Solver%ncoeff1st
          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!         1/ Rotational invariances (2nd order). Number of constraints = natom_unitcell*3**2
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do icoeff=1,ncoeff
                  terme1=sum(SS_ref(alpha,:,beta,isym,itrans)*Shell2at%proj(:,icoeff,ishell))*distance(iatom,jatom,gama+1)
                  terme2=sum(SS_ref(alpha,:,gama,isym,itrans)*Shell2at%proj(:,icoeff,ishell))*distance(iatom,jatom,beta+1)
                  const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)=&
&                 const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)+terme1-terme2                  
                end do
              end do
            end do
          end do
!         2/ Enforce the symetry of the dynamical matrix. Number of constraints = (3*natom_unitcell)**2
!            Note that we are unable to enforce the symetry when iatom=jatom (We have to write the equations)
          do alpha=1,3
            do beta=1,3
              do icoeff=1,ncoeff
                terme=sum(SS_ref(alpha,:,beta,isym,itrans)*Shell2at%proj(:,icoeff,ishell))-&
&                     sum(SS_ref(beta,:,alpha,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
                const_dynmat(alpha,beta,iatom,iat_mod,icoeff+ncoeff_prev)=&
&               const_dynmat(alpha,beta,iatom,iat_mod,icoeff+ncoeff_prev)+terme
              end do
            end do
          end do
!         3/ Huang invariances. Number of constraints = 3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  do icoeff=1,ncoeff
                    terme=sum(SS_ref(alpha,:,beta,isym,itrans)*Shell2at%proj(:,icoeff,ishell))*&
&                             distance(iatom,jatom,gama+1)*&
&                             distance(iatom,jatom,lambda+1)-&
&                         sum(SS_ref(gama,:,lambda,isym,itrans)*Shell2at%proj(:,icoeff,ishell))*&
&                             distance(iatom,jatom,alpha+1)*&
&                             distance(iatom,jatom,beta+1)
                    const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)=&
&                   const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)+terme
                  end do
                end do
              end do
            end do
          end do
        end do !iatshell
      end do !iatom
    end do !ishell
  end if !order=1,2

! Third order
  if (order3) then
    write(Invar%stdout,*) '########################## At the 3rd order #################################'
    ABI_MALLOC(Const3%AsrRot3,(natom_unitcell,natom,ntotcoeff))
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do icoeff=1,ntotcoeff
          ABI_MALLOC(Const3%AsrRot3(iatom,jatom,icoeff)%ABG,   (3,3,3)); Const3%AsrRot3(iatom,jatom,icoeff)%ABG(:,:,:)   =zero
          ABI_MALLOC(Const3%AsrRot3(iatom,jatom,icoeff)%ABGD,(3,3,3,3)); Const3%AsrRot3(iatom,jatom,icoeff)%ABGD(:,:,:,:)=zero
        end do  
      end do  
    end do  
    do ishell=1,Shell2at%nshell
      do iatom=1,natom_unitcell
        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
          if (iatom==jatom) cycle
          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell2at%ncoeff(ishell)
          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+Solver%ncoeff1st
          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!         1/ Rotational invariances (coming from the 2nd order). Number of constraints = natom_unitcell*natom*3**3
          if (Invar%order.ge.3) then
            do alpha=1,3
              do beta=1,3
                do gama=1,3
                  do lambda=1,3
                    do icoeff=1,ncoeff
                      terme1=zero ; terme2=zero ; terme3=zero ; terme4=zero ;
                      if (alpha.eq.lambda) terme1=sum(SS_ref(gama  ,:,beta  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
                      if (beta.eq.lambda)  terme2=sum(SS_ref(alpha ,:,gama  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
                      if (alpha.eq.gama)   terme3=sum(SS_ref(lambda,:,beta  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
                      if (beta.eq.gama)    terme4=sum(SS_ref(alpha ,:,lambda,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
                      if (distance(iatom,jatom,1).lt.Invar%rcut3) then
                        Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                       Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1+terme2-terme3-terme4
                      end if
                      Const3%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                     Const3%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)-terme1-terme2+terme3+terme4
                    end do
                  end do
                end do
              end do    
            end do  
          end if !proj3rd
        end do !iatshell
      end do !iatom
    end do !ishell
    do ishell=1,Shell3at%nshell
      do iatom=1,natom_unitcell
        if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell3at%ncoeff(ishell)
          ncoeff_prev=Shell3at%ncoeff_prev(ishell)+Solver%ncoeff2nd+Solver%ncoeff1st
!         2/ Acoustic sum rules (3rd order). Number of constraints = natom_unitcell*natom*3**3
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do icoeff=1,ncoeff
                  terme =sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,gama)*Shell3at%proj(:,icoeff,ishell))
                    Const3%AsrRot3(iatom,katom,icoeff+ncoeff_prev)%ABG(alpha,beta,gama)=&
&                   Const3%AsrRot3(iatom,katom,icoeff+ncoeff_prev)%ABG(alpha,beta,gama)+terme
                end do
              end do
            end do
          end do
!         2/ Rotational invariances (coming from the 3rd order). Number of constraints = natom_unitcell*natom*3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  do icoeff=1,ncoeff
                    terme1=sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,gama  )&
&                          *Shell3at%proj(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
                    terme2=sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,lambda)&
&                          *Shell3at%proj(:,icoeff,ishell))*distance(iatom,katom,gama+1)
                    Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                   Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1-terme2
                  end do
                end do
              end do
            end do    
          end do  
        end do !iatshell
      end do !iatom   
    end do !ishell   
  end if !order=3

! Fourth order
  if (order4) then
    write(Invar%stdout,*) '########################## At the 4th order #################################'
    ABI_MALLOC(Const4%AsrRot4,(natom_unitcell,natom,natom,ntotcoeff))
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do icoeff=1,ntotcoeff
            ABI_MALLOC(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD,   (3,3,3,3))
                       Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD(:,:,:,:)   =zero
!FB            ABI_MALLOC(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE,(3,3,3,3,3))
!FB                       Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE(:,:,:,:,:)=zero
          end do  
        end do  
      end do  
    end do  
!FB    do ishell=1,Shell2at%nshell
!FB      do iatom=1,natom_unitcell
!FB        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!FB        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
!FB          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
!FB          if (iatom==jatom) cycle
!FB          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
!FB          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
!FB          ncoeff     =Shell2at%ncoeff(ishell)
!FB          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+Solver%ncoeff1st
!FB          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!FB!         1/ Rotational invariances (coming from the 2nd order). Number of constraints = natom_unitcell*natom*3**3
!FB          if (Invar%order.ge.3) then
!FB            do alpha=1,3
!FB              do beta=1,3
!FB                do gama=1,3
!FB                  do lambda=1,3
!FB                    do icoeff=1,ncoeff
!FB                      terme1=zero ; terme2=zero ; terme3=zero ; terme4=zero ;
!FB                      if (alpha.eq.lambda) terme1=sum(SS_ref(gama  ,:,beta  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
!FB                      if (beta.eq.lambda)  terme2=sum(SS_ref(alpha ,:,gama  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
!FB                      if (alpha.eq.gama)   terme3=sum(SS_ref(lambda,:,beta  ,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
!FB                      if (beta.eq.gama)    terme4=sum(SS_ref(alpha ,:,lambda,isym,itrans)*Shell2at%proj(:,icoeff,ishell))
!FB                      if (distance(iatom,jatom,1).lt.Invar%rcut3) then
!FB                        Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                       Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1+terme2-terme3-terme4
!FB                      end if
!FB                      Const4%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                     Const4%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)-terme1-terme2+terme3+terme4
!FB                    end do
!FB                  end do
!FB                end do
!FB              end do    
!FB            end do  
!FB          end if !proj3rd
!FB        end do !iatshell
!FB      end do !iatom
!FB    end do !ishell
    do ishell=1,Shell4at%nshell
      do iatom=1,natom_unitcell
        if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
          isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell4at%ncoeff(ishell)
          ncoeff_prev=Shell4at%ncoeff_prev(ishell)+Solver%ncoeff3rd+Solver%ncoeff2nd+Solver%ncoeff1st
!         2/ Acoustic sum rules (4th order). Number of constraints = natom_unitcell*natom**2*3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do delta=1,3
                  do icoeff=1,ncoeff
                    terme =sum(Const4%Sprod(isym,itrans)%SSSS(alpha,:,beta,gama,delta)*Shell4at%proj(:,icoeff,ishell))
                      Const4%AsrRot4(iatom,katom,latom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,delta)=&
&                     Const4%AsrRot4(iatom,katom,latom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,delta)+terme
                  end do
                end do
              end do
            end do
          end do
!FB!         2/ Rotational invariances (coming from the 3rd order). Number of constraints = natom_unitcell*natom*3**4
!FB          do alpha=1,3
!FB            do beta=1,3
!FB              do gama=1,3
!FB                do lambda=1,3
!FB                  do icoeff=1,ncoeff
!FB                    terme1=sum(Const4%Sprod(isym,itrans)%SSS(alpha,:,beta,gama  )&
!FB                               &*Shell3at%proj(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
!FB                    terme2=sum(Const4%Sprod(isym,itrans)%SSS(alpha,:,beta,lambda)&
!FB                               &*Shell3at%proj(:,icoeff,ishell))*distance(iatom,katom,gama+1)
!FB                    Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                   Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1-terme2
!FB                  end do
!FB                end do
!FB              end do
!FB            end do    
!FB          end do  
        end do !iatshell
      end do !iatom   
    end do !ishell   
  end if !order=4

  if (order2.or.order3) then
    ABI_FREE(SS_ref)
  end if  
  if (order3.or.order4) then
    do isym=1,Sym%nsym
      do itrans=1,6
        ABI_FREE(Const3%Sprod(isym,itrans)%SSS)
      end do
    end do
    ABI_FREE(Const3%Sprod)
  end if  
  if (order4) then
    do isym=1,Sym%nsym
      do itrans=1,24
        ABI_FREE(Const4%Sprod(isym,itrans)%SSSS)
      end do
    end do
    ABI_FREE(Const4%Sprod)
  end if  

! Reduce the number of constraints by selecting the non-zero equations
  write(Invar%stdout,*) '################## Reduce the number of constraints #########################'
  iconst_new=0
  if (order2) then
!   1/ For Rotational invariances (1st order)
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,Solver%nconst_1st)) ; vect(:,:)=zero
    do alpha=1,3
      do beta=1,3
        iconst=iconst+1 
        vect(:,iconst)=const_rot1st(alpha,beta,:)
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_1st,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do
    end if
    ABI_FREE(vect)
    ABI_FREE(const_rot1st)
    Solver%nconst_1st=nconst_loc

!   2/ For Rotational invariances (2nd order)
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,Solver%nconst_rot2nd)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do alpha=1,3
        do beta=1,3
          do gama=1,3
            iconst=iconst+1  
            vect(:,iconst)=const_rot2nd(alpha,beta,gama,iatom,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_rot2nd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do
    end if
    ABI_FREE(vect)
    ABI_FREE(const_rot2nd)
    Solver%nconst_rot2nd=nconst_loc

!   3/ For symetry of the dynamical matrix
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,Solver%nconst_dynmat)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom_unitcell
        do alpha=1,3
          do beta=1,3
            iconst=iconst+1
            vect(:,iconst)=const_dynmat(alpha,beta,iatom,jatom,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_dynmat,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do
    end if
    ABI_FREE(vect)
    ABI_FREE(const_dynmat)
    Solver%nconst_dynmat=nconst_loc

!   4/ For Huang invariances
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,Solver%nconst_huang)) ; vect(:,:)=zero
    do alpha=1,3
      do beta=1,3
        do gama=1,3
          do lambda=1,3
            iconst=iconst+1
            vect(:,iconst)=const_huang(alpha,beta,gama,lambda,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_huang,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do
    end if
    ABI_FREE(vect)
    ABI_FREE(const_huang)
    Solver%nconst_huang=nconst_loc
    Solver%nconst_2nd=Solver%nconst_rot2nd+Solver%nconst_dynmat+Solver%nconst_huang
  end if

  if (order3) then
!   1/ For acoustic sum rules (3rd order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,Solver%nconst_asr3rd)) ; vect (:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3 
              iconst=iconst+1
              do ii=1,ntotcoeff
                vect(ii,iconst)=Const3%AsrRot3(iatom,jatom,ii)%ABG(alpha,beta,gama)
              end do  
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_asr3rd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    Solver%nconst_asr3rd=nconst_loc

!   2/ For Rotational invariances (3rd order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,Solver%nconst_rot3rd)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                iconst=iconst+1
                do ii=1,ntotcoeff
                  vect(ii,iconst)=Const3%AsrRot3(iatom,jatom,ii)%ABGD(alpha,beta,gama,lambda)
                  end do  
              end do
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_rot3rd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do icoeff=1,ntotcoeff
          ABI_FREE(Const3%AsrRot3(iatom,jatom,icoeff)%ABG)
          ABI_FREE(Const3%AsrRot3(iatom,jatom,icoeff)%ABGD)
        end do  
      end do  
    end do  
    ABI_FREE(Const3%AsrRot3)
    Solver%nconst_rot3rd=nconst_loc
    Solver%nconst_3rd=Solver%nconst_asr3rd+Solver%nconst_rot3rd
  end if

  if (order4) then
!   1/ For acoustic sum rules (4th order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,Solver%nconst_asr4th)) ; vect (:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do alpha=1,3
            do beta=1,3
              do gama=1,3 
                do delta=1,3 
                  iconst=iconst+1
                  do ii=1,ntotcoeff
                    vect(ii,iconst)=Const4%AsrRot4(iatom,jatom,katom,ii)%ABGD(alpha,beta,gama,delta)
                  end do  
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_asr4th,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        Solver%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    Solver%nconst_asr4th=nconst_loc

!FB!   2/ For Rotational invariances (3rd order)
!FB    iconst=0
!FB    ABI_MALLOC(vect ,(ntotcoeff,Solver%nconst_rot3rd)) ; vect(:,:)=zero
!FB    do iatom=1,natom_unitcell
!FB      do jatom=1,natom
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                iconst=iconst+1
!FB                vect(:,iconst)=Const4%AsrRot3(iatom,jatom,:)%ABGD(alpha,beta,gama,lambda)
!FB              end do
!FB            end do
!FB          end do
!FB        end do
!FB      end do
!FB    end do
!FB    call tdep_calc_orthonorm(ntotcoeff,Solver%nconst_rot3rd,nconst_loc,vect)
!FB    if (nconst_loc.ne.0) then
!FB      do iconst_loc=1,nconst_loc
!FB        iconst_new=iconst_new+1
!FB        Solver%const(iconst_new,:)=vect(:,iconst_loc)
!FB      end do  
!FB    end if
!FB    ABI_FREE(vect)
!FB    Solver%nconst_rot3rd=nconst_loc
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do icoeff=1,ntotcoeff
            ABI_FREE(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD)
!FB            ABI_FREE(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE)
          end do  
        end do  
      end do  
    end do  
    ABI_FREE(Const4%AsrRot4)
    Solver%nconst_rot4th=0
    Solver%nconst_4th=Solver%nconst_asr4th+Solver%nconst_rot4th
  end if

! Finalize the orthonormalization   
!FB  nconst=iconst_new
!FB  ABI_MALLOC(vect,(ntotcoeff,nconst)) ; vect(:,:)=zero
!FB  do iconst=1,nconst
!FB    vect(:,iconst)=Solver%const(iconst,:)
!FB  end do
!FB  ABI_FREE(Solver%const)
!FB  call tdep_calc_orthonorm(ntotcoeff,nconst,nconst_loc,vect)
!FB  Solver%ntotconst=nconst_loc
!FB  if (nconst_loc.ne.0) then
!FB    ABI_MALLOC(Solver%const ,(Solver%ntotconst,ntotcoeff)); Solver%const (:,:)=0.d0 
!FB    do iconst_loc=1,nconst_loc
!FB      Solver%const(iconst_loc,:)=vect(:,iconst_loc)
!FB    end do  
!FB  end if
!FB  ABI_FREE(vect)
  Solver%ntotconst=iconst_new
  if (MPIdata%iam_master) then
    open(unit=16,file=trim(Invar%output_prefix)//'_constraints.dat')
    write(16,*) ' ======== Constraints at the 1st order (Rotational Invariances) ========'
    write(16,*) ' Number of constraints =',Solver%nconst_1st
    write(16,*) ' ======== Constraints at the 2nd order (Rotational Invariances) ========'
    write(16,*) ' Number of constraints =',Solver%nconst_rot2nd
    write(16,*) ' ======== Constraints at the 2nd order (Dynamical Matrix) =============='
    write(16,*) ' Number of constraints =',Solver%nconst_dynmat
    write(16,*) ' ======== Constraints at the 2nd order (Huang) ========================='
    write(16,*) ' Number of constraints =',Solver%nconst_huang
    if (Invar%order.ge.3) then
      write(16,*) ' ======== Constraints at the 3rd order (Acoustic sum rules) ============'
      write(16,*) ' Number of constraints =',Solver%nconst_asr3rd
      write(16,*) ' ======== Constraints at the 3rd order (Rotational Invariances) ========'
      write(16,*) ' Number of constraints =',Solver%nconst_rot3rd
    end if
    if (Invar%order.ge.4) then
      write(16,*) ' ======== Constraints at the 4th order (Acoustic sum rules) ============'
      write(16,*) ' Number of constraints =',Solver%nconst_asr4th
      write(16,*) ' ======== Constraints at the 4th order (Rotational Invariances) ========'
      write(16,*) ' Number of constraints =',Solver%nconst_rot4th
    end if
    write(16,*) ' ======================================================================='
    write(16,*) ' Total number of constraints =',Solver%ntotconst
    close(16)
  end if  

end subroutine tdep_calc_constraints

!====================================================================================================

end module m_tdep_solver
!!***
