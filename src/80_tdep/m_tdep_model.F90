!!****m* ABINIT/m_tdep_model
!! NAME
!!  m_tdep_model
!!
!! FUNCTION
!!  This module contains the TDEP Model data type
!!  which holds the IFC at all orders and derived quantities.
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

module m_tdep_model

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools
 use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
 use m_tdep_sampling,    only : tdep_Sampling_type
 use m_tdep_shell,       only : Shell_type

 implicit none

!----------------------------------------------------------------------

!!****t* m_tdep_model/Phi2_type
!! NAME
!! Phi2_type
!!
!! FUNCTION
!!  Second-order IFC
!!
!! SOURCE

 type Phi2_type

   double precision, allocatable :: SR(:,:)
   ! SR(3*natom, 3*natom)
   ! Short-range part of the second-order IFC

   double precision, allocatable :: LR(:,:)
   ! LR(3*natom, 3*natom)
   ! Long-range part of the second-order IFC

   double precision, allocatable :: Tot(:,:)
   ! Tot(3*natom, 3*natom)
   ! Second-order IFC (SR + LR)


 end type Phi2_type


!----------------------------------------------------------------------

!!****t* m_tdep_model/tdep_Model_type
!! NAME
!! tdep_Model_type
!!
!! FUNCTION
!!  Holds the result of the TDEP fit, that is, the IFCs and resulting forces
!!  and energy contributions.
!!
!! SOURCE

 type,public :: tdep_Model_type

   integer :: order
   ! Maximum IFC order (2, 3, or 4)

   integer :: natom
   ! Number of atoms

   integer :: natom_unitcell
   ! Number of atoms in the unitcell

   integer :: my_nstep
   ! Number of MD steps held locally

   integer :: nshell3rd
   ! Number of shells for 3rd order IFC

   integer :: nshell4th
   ! Number of shells for 4th order IFC

   double precision :: U0
   ! Constant term contribution to energy

   double precision, allocatable :: Forces(:)
   ! Forces(3*natom*my_nstep)
   ! Forces computed from the fitted IFC and the displacements

   double precision, allocatable :: Phi1(:)
   ! Phi1(3*natom)
   ! First-order IFC

   type(Phi2_type) :: Phi2
   ! Second-order IFC

   double precision, allocatable :: Phi3(:,:,:,:)
   ! Phi3(3,3,3,Shell3at%nshell)
   ! Third-order IFC

   double precision, allocatable :: Phi4(:,:,:,:,:)
   ! Phi3(3,3,3,3,Shell4at%nshell)
   ! Fourth-order IFC

   double precision, allocatable :: Phi1Ui(:)
   ! Phi1Ui(my_nstep)
   ! First-order IFC contribution to energy

   double precision, allocatable :: Phi2UiUj(:)
   ! Phi2UiUj(my_nstep)
   ! Second-order IFC contribution to energy

   double precision, allocatable :: Phi3UiUjUk(:)
   ! Phi3UiUjUk(my_nstep)
   ! Third-order IFC contribution to energy

   double precision, allocatable :: Phi4UiUjUkUl(:)
   ! Phi4UiUjUkUl(my_nstep)
   ! Fourth-order IFC contribution to energy

 end type tdep_Model_type

 public :: tdep_init_phi2
 public :: tdep_destroy_phi2
 public :: tdep_model_init
 public :: tdep_model_free
 public :: tdep_calc_model

contains

!=====================================================================================================

subroutine tdep_init_phi2(Phi2,dipdip,natom)

  type(Phi2_type),intent(out) :: Phi2
  logical, intent(in) :: dipdip
  integer, intent(in) :: natom

  ABI_CALLOC(Phi2%SR ,(3*natom,3*natom))
  if (dipdip) then
    ABI_CALLOC(Phi2%Tot,(3*natom,3*natom))
    ABI_CALLOC(Phi2%LR ,(3*natom,3*natom))
  end if

end subroutine tdep_init_phi2

!=====================================================================================================

subroutine tdep_destroy_phi2(Phi2)

  type(Phi2_type),intent(inout) :: Phi2

  ABI_FREE(Phi2%SR)
  ABI_SFREE(Phi2%Tot)
  ABI_SFREE(Phi2%LR)

end subroutine tdep_destroy_phi2

!=====================================================================================================

subroutine tdep_model_init(Model, Invar, Shell3at, Shell4at)

 type(tdep_Model_type),intent(out) :: Model
 type(atdep_dataset_type),intent(in) :: Invar
 type(Shell_type),intent(in) :: Shell3at, Shell4at

 Model%order = Invar%order
 Model%natom = Invar%natom
 Model%natom_unitcell = Invar%natom_unitcell
 Model%my_nstep = Invar%my_nstep
 Model%nshell3rd = 1
 Model%nshell4th = 1

 Model%U0 = zero

 ABI_CALLOC(Model%Forces, (3*Model%natom*Model%my_nstep))
 ABI_CALLOC(Model%Phi1Ui      ,(Model%my_nstep))
 ABI_CALLOC(Model%Phi2UiUj    ,(Model%my_nstep))
 ABI_CALLOC(Model%Phi3UiUjUk  ,(Model%my_nstep))
 ABI_CALLOC(Model%Phi4UiUjUkUl,(Model%my_nstep))

 ABI_CALLOC(Model%Phi1,(3*Model%natom))

 call tdep_init_phi2(Model%Phi2,Invar%loto,Invar%natom)

 if (Model%order.ge.3) then
   Model%nshell3rd = Shell3at%nshell
   ABI_CALLOC(Model%Phi3,(3,3,3,Model%nshell3rd))
 end if

 if (Model%order.eq.4) then
   Model%nshell4th = Shell4at%nshell
   ABI_CALLOC(Model%Phi4,(3,3,3,3,Model%nshell4th))
 end if

end subroutine tdep_model_init

!=====================================================================================================

subroutine tdep_model_free(Model)

 type(tdep_Model_type),intent(inout) :: Model

 ABI_FREE(Model%Forces)
 ABI_FREE(Model%Phi1Ui)
 ABI_FREE(Model%Phi2UiUj)
 ABI_FREE(Model%Phi3UiUjUk)
 ABI_FREE(Model%Phi4UiUjUkUl)

 ABI_FREE(Model%Phi1)

 call tdep_destroy_phi2(Model%Phi2)

 if (Model%order.ge.3) then
   ABI_FREE(Model%Phi3)
 end if

 if (Model%order.eq.4) then
   ABI_FREE(Model%Phi4)
 end if

end subroutine tdep_model_free

!====================================================================================================

 subroutine tdep_calc_model(Model,MD,Invar,MPIdata)

  type(tdep_Model_type),intent(inout) :: Model
  type(tdep_Sampling_type),intent(in) :: MD
  type(atdep_dataset_type),intent(in) :: Invar
  type(MPI_enreg_type), intent(in) :: MPIdata

  integer :: ii,jj,istep,iatom
  double precision :: Delta_F2,Delta_U,Delta_U2
  double precision :: sigma,U_1,U_2,U_3,U_4,UMD
  double precision, allocatable :: tmp(:),Phi_tot(:)
  double precision, allocatable :: U_MD(:),U_TDEP(:),weights_tot(:)
  integer :: ierr

  ! Compute the different contributions to total energy from the model

  write(Invar%stdout,*)' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '######################### Energies, errors,...  #############################'
  write(Invar%stdout,*) '#############################################################################'

! Compute U0, U_TDEP, Delta_U and write them in the output file
  write(Invar%stdout,'(a)') ' Thermodynamic quantities and convergence parameters of THE MODEL,'
  write(Invar%stdout,'(a)') '      as a function of the step number (energies in eV/atom and forces in Ha/bohr) :'
  if (Invar%order.eq.4) then
    write(Invar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2 + U_3 + U_4'
    write(Invar%stdout,'(2a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj ',&
&                             '- 1/6 sum_ijk Phi3 ui uj uk - 1/24 sum_ijkl Phi4 ui uj uk ul >'
    write(Invar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(Invar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
    write(Invar%stdout,'(a)') '        and U_3 = < 1/6  sum_ijk  Phi3 ui uj uk >'
    write(Invar%stdout,'(a)') '        and U_4 = < 1/24 sum_ijkl Phi4 ui uj uk ul >'
  else if (Invar%order.eq.3) then
    write(Invar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2 + U_3'
    write(Invar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj - 1/6 sum_ijk Phi3 ui uj uk >'
    write(Invar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(Invar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
    write(Invar%stdout,'(a)') '        and U_3 = < 1/6  sum_ijk  Phi3 ui uj uk >'
  else
    write(Invar%stdout,'(a)') ' <U_TDEP> = U_0 + U_1 + U_2'
    write(Invar%stdout,'(a)') '       with U_0 = < U_MD - sum_i Phi1 ui - 1/2 sum_ij Phi2 ui uj >'
    write(Invar%stdout,'(a)') '        and U_1 = <      sum_i    Phi1 ui >'
    write(Invar%stdout,'(a)') '        and U_2 = < 1/2  sum_ij   Phi2 ui uj >'
  end if
  write(Invar%stdout,'(a)') '  Delta_U =   < U_MD - U_TDEP > '
  write(Invar%stdout,'(a)') '  Delta_U2= (< (U_MD - U_TDEP)^2 >)**0.5 '
  write(Invar%stdout,'(a)') '  Delta_F2= (< (F_MD - F_TDEP)^2 >)**0.5 '
  write(Invar%stdout,'(a)') '  Sigma   = (< (F_MD - F_TDEP)^2 >/<F_MD**2>)**0.5 '
  if (Invar%order.eq.4) then
    write(Invar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '            U_3              U_4            Delta_U          Delta_U2          Delta_F2          Sigma'
  else if (Invar%order.eq.3) then
    write(Invar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '            U_3            Delta_U          Delta_U2          Delta_F2          Sigma'
  else
    write(Invar%stdout,'(2a)') '     <U_MD>            U_0              U_1              U_2  ',&
&     '          Delta_U          Delta_U2          Delta_F2          Sigma'
  end if

! Compute eucledian distance for forces
  ABI_MALLOC(tmp,(11))                  ; tmp(:)   =0.d0
  do istep=1,MD%my_nstep
    do iatom=1,MD%natom
      do ii=1,3
        jj = ii + 3*(iatom-1) + 3*MD%natom*(istep-1)
        tmp(4)=tmp(4)+(MD%Forces(jj)-Model%Forces(jj))**2*MD%weights(istep)
        tmp(5)=tmp(5)+MD%Forces(jj)**2*MD%weights(istep)
      end do
    end do
  end do
! Compute energies
  ABI_MALLOC(U_TDEP,     (MD%nstep_tot)) ; U_TDEP(:)=0.d0
  ABI_MALLOC(U_MD,       (MD%nstep_tot)) ; U_MD(:)  =0.d0
  ABI_MALLOC(weights_tot,(MD%nstep_tot)) ; weights_tot(:)=0.d0
  ABI_MALLOC(Phi_tot,    (MPIdata%my_nstep)); Phi_tot(:)=0.d0
  do istep=1,MD%my_nstep
    tmp(7) =tmp(7) +MD%etot(istep)*MD%weights(istep)
    tmp(10)=tmp(10)+Model%Phi1Ui(istep)*MD%weights(istep)
    tmp(6) =tmp(6) +Model%Phi2UiUj(istep)*MD%weights(istep)
    tmp(8) =tmp(8) +Model%Phi3UiUjUk(istep)*MD%weights(istep)
    tmp(11)=tmp(11)+Model%Phi4UiUjUkUl(istep)*MD%weights(istep)
  end do
  call xmpi_sum(tmp,MPIdata%comm_step,ierr)
  tmp(1) = tmp(7)-tmp(10)-tmp(6)-tmp(8)-tmp(11)
  Phi_tot(:)=tmp(1)+Model%Phi1Ui(:)+Model%Phi2UiUj(:)+Model%Phi3UiUjUk(:)+Model%Phi4UiUjUkUl(:)
  call xmpi_gatherv(Phi_tot,MD%my_nstep,U_TDEP,MPIdata%nstep_all,MPIdata%shft_step,&
&                   MPIdata%master,MPIdata%comm_step,ierr)
  call xmpi_gatherv(MD%etot,MD%my_nstep,U_MD,MPIdata%nstep_all,MPIdata%shft_step,&
&                   MPIdata%master,MPIdata%comm_step,ierr)
  call xmpi_gatherv(MD%weights,MD%my_nstep,weights_tot,MPIdata%nstep_all,MPIdata%shft_step,&
&                   MPIdata%master,MPIdata%comm_step,ierr)
  do istep=1,MD%nstep_tot
    tmp(2) =tmp(2) + (U_MD(istep)-U_TDEP(istep)) * weights_tot(istep)
    tmp(9) =tmp(9) + (U_MD(istep)-U_TDEP(istep))**2 * weights_tot(istep)
  end do
  Model%U0       =tmp(1) /real(MD%natom)
  UMD      =tmp(7) /real(MD%natom)
  U_1      =tmp(10)/real(MD%natom)
  U_2      =tmp(6) /real(MD%natom)
  U_3      =tmp(8) /real(MD%natom)
  U_4      =tmp(11)/real(MD%natom)
  Delta_U  =tmp(2) /real(MD%natom)
  Delta_U2 =tmp(9) /real(MD%natom)
  Delta_F2 =tmp(4) /real(MD%natom*3)
  if (tmp(5).eq.0.d0) then
    sigma = 0.d0
  else
    sigma = dsqrt(tmp(4)/tmp(5))
  end if
  if (Invar%order.eq.4) then
    write(Invar%stdout,'(10(f12.5,5x))') UMD*Ha_eV,Model%U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,U_3*Ha_eV,U_4*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  else if (Invar%order.eq.3) then
    write(Invar%stdout,'(9(f12.5,5x))') UMD*Ha_eV,Model%U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,U_3*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  else
    write(Invar%stdout,'(8(f12.5,5x))') UMD*Ha_eV,Model%U0*Ha_eV,U_1*Ha_eV,U_2*Ha_eV,&
&     Delta_U*Ha_eV,Delta_U2**0.5*Ha_eV,Delta_F2**0.5,sigma
  endif
  ABI_FREE(tmp)
  write(Invar%stdout,'(a,1x,f12.5)') ' NOTE : in the harmonic and classical limit (T>>T_Debye), U_2=3/2*kB*T=',&
&   3.d0/2.d0*kb_HaK*Ha_eV*Invar%temperature

! Write : i) (U_TDEP vs U_MD) in etotMDvsTDEP.dat
!        ii) (Model%Forces vs MD%Forces) in fcartMDvsTDEP.dat
  write(Invar%stdout,'(a)') ' '
  write(Invar%stdout,'(a)') ' See the etotMDvsTDEP.dat & fcartMDvsTDEP.dat files'
  if (MPIdata%iam_master) then
    open(unit=32,file=trim(Invar%output_prefix)//'_etotMDvsTDEP.dat')
    open(unit=33,file=trim(Invar%output_prefix)//'_fcartMDvsTDEP.dat')
    write(32,'(a)') '#   Istep      U_MD(Ha)         U_TDEP(Ha)'
    write(33,'(a)') '# Forces_MD(Ha/bohr) Forces_TDEP(Ha/bohr)'
    do istep=1,MD%nstep_tot
      write(32,'(i6,1x,2(f17.6,1x))') istep,U_MD(istep),U_TDEP(istep)
    end do
    do istep=1,MD%my_nstep
      do iatom=1,MD%natom
        do ii=1,3
          write(33,'(2(f17.10,1x))') MD%Forces  (ii+3*(iatom-1)+3*MD%natom*(istep-1)),&
&                                    Model%Forces(ii+3*(iatom-1)+3*MD%natom*(istep-1))
        end do
      end do
    end do
    close(32)
    close(33)
  end if
  ABI_FREE(U_MD)
  ABI_FREE(U_TDEP)
  ABI_FREE(Phi_tot)
  ABI_FREE(weights_tot)

 end subroutine tdep_calc_model

!====================================================================================================

end module m_tdep_model
!!***
