!!****m* ABINIT/m_rttddft_tdef
!! NAME
!!  m_rttddft_tdef
!!
!! FUNCTION
!!  Contains definition of the tdef type
!!  related to time-dependent electric field
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2023 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft_tdef

 use defs_basis
 use defs_abitypes,   only: MPI_type

 use m_dtset,         only: dataset_type
 use m_errors,        only: msg_hndl
 use m_initylmg,      only: initylmg
 use m_profiling_abi, only: abimem_record
 use m_xmpi,          only: xmpi_bcast

 implicit none

 private
!!***

!! NAME
!! tdef_type: Time Dependent Electric Field type
!! Object containing the TD Electric field
!! for RT-TDDFT runs
!!
!! SOURCE
 type,public :: tdef_type

   integer  :: ef_type       !type of TD electric feld (Dirac or sin^2 pulse)
   real(dp) :: efield(3)     !TD external elec. field perturbation
   real(dp) :: ef_ezero(3)   !E_0 = |E_0|*polarization
   real(dp) :: vecpot(3)     !total vector potential
   real(dp) :: vecpot_ext(3) !external vector potential
   real(dp) :: vecpot_ind(3,2) !induced vector potential at t and t-dt
   real(dp) :: vecpot_red(3) !total vector potential in reduced coord. (in reciprocal space)
   real(dp) :: ef_tzero      !time at which elec field is switched on
   real(dp) :: ef_omega      !angular freq. of TD elec field
   real(dp) :: ef_tau        !time width of the pulse
   real(dp) :: ef_sin_a      !useful constant for sin^2 pulse
   real(dp) :: ef_sin_b      !useful constant for sin^2 pulse
   real(dp), allocatable :: kpa(:,:) ! contains kpts + A
                                     ! (in reduced coordinates in reciprocal space)
   logical  :: induced_vecpot ! Add the vector potential induced by the current density
                              ! in the Hamiltonian

   contains

   procedure :: init => tdef_init
   procedure :: update => tdef_update
   procedure :: restart => tdef_restart

 end type tdef_type
!!***

contains
!!***

!!****f* m_rttddft/tdef_init
!!
!! NAME
!!  tdef_init
!!
!! FUNCTION
!!  Update value of electric field and vector potential at time t
!!
!! INPUTS
!!  [tdef = tdef structure to update]
!!  td_ef_type = type of electric field (Dirac or sin^2 pulse)
!!  td_ef_pol = polarization
!!  td_ef_ezero = Amplitude (E_0 = |E_0|*polarization)
!!  td_ef_tzero = time at which the pulse is switched on
!!  td_ef_lambda = wavelength (for sin^2 pulse)
!!  td_ef_tau = time width of the pulse (for sin^2 pulse)
!!  time = propagation time
!!  nkpt = number of kpoints
!!  kpts = kpoints array
!!
!! OUTPUT
!!  [tdef = updated tdef structure]
!!
!! SOURCE
subroutine tdef_init(tdef,td_ef_type,td_ef_pol,td_ef_ezero,td_ef_tzero,td_ef_lambda,td_ef_tau,td_ef_induced_vecpot,nkpt,kpts)

 !Arguments ------------------------------------
 !scalars
 class(tdef_type), intent(inout) :: tdef
 integer,          intent(in)    :: td_ef_type, td_ef_induced_vecpot, nkpt
 real(dp),         intent(in)    :: td_ef_ezero, td_ef_tzero, td_ef_lambda, td_ef_tau
 !arrays
 real(dp),         intent(in)    :: kpts(:,:)

 !Local variables-------------------------------
 real(dp),         intent(in)    :: td_ef_pol(3)

! ***********************************************************************

 if (td_ef_type > 1 .or. td_ef_type < 0) then
   ABI_ERROR("Wrong value of td_ef_type")
 end if

 tdef%ef_type  = td_ef_type
 tdef%ef_ezero = td_ef_pol*td_ef_ezero
 tdef%ef_tau   = td_ef_tau
 tdef%ef_omega = 2.0_dp*pi*Sp_Lt/td_ef_lambda !2*pi*f=2*pi*c/lambda
 tdef%ef_tzero = td_ef_tzero
 tdef%ef_sin_a = 2.0_dp*pi/td_ef_tau + tdef%ef_omega
 tdef%ef_sin_b = 2.0_dp*pi/td_ef_tau - tdef%ef_omega
 if (td_ef_induced_vecpot == 0) then
   tdef%induced_vecpot = .false.
 else if (td_ef_induced_vecpot == 1) then
   tdef%induced_vecpot = .true.
 else
    ABI_ERROR("Wrong value of td_ef_induced_vecpot")
 end if

 tdef%efield = zero
 tdef%vecpot = zero
 tdef%vecpot_ext = zero
 tdef%vecpot_ind = zero
 tdef%vecpot_red = zero

 ABI_MALLOC(tdef%kpa,(3,nkpt))
 tdef%kpa = kpts

end subroutine tdef_init
!!***

!!****f* m_rttddft/tdef_update
!!
!! NAME
!!  tdef_update
!!
!! FUNCTION
!!  Update value of electric field and vector potential at time t
!!
!! INPUTS
!!  [tdef = tdef structure to update]
!!  dtset = dataset structure
!!  mpi_enreg = MPI communicators structure
!!  time = propagation time
!!  rprimd = cell vectors (direct space)
!!  gprimd = cell vectors (reciprocal space)
!!  kg = kpoints in reciprocal space
!!  mpsang = 1+maximum angular momentum for nonlocal pseudopotential (required by initylmg)
!!  npwarr = array holding npw for each k point
!!  ylm = real spherical harmonics for each G and k point
!!  ylmgr = gradient of real spherical harmonics for each G and k point
!!  current = total current density
!!
!! OUTPUT
!!  [tdef = updated tdef structure]
!!
!! SOURCE
subroutine tdef_update(tdef,dtset,mpi_enreg,time,rprimd,gprimd,kg,mpsang,npwarr,ylm,ylmgr,current,update_vecpot_ind)

 !Arguments ------------------------------------
 !scalars
 class(tdef_type),   intent(inout) :: tdef
 type(dataset_type), intent(inout) :: dtset
 type(MPI_type),     intent(inout) :: mpi_enreg
 real(dp),           intent(in)    :: time
 real(dp),           intent(in)    :: rprimd(3,3)
 real(dp),           intent(in)    :: gprimd(3,3)
 integer,            intent(in)    :: kg(:,:)
 integer,            intent(in)    :: mpsang
 integer,            intent(in)    :: npwarr(:)
 real(dp),           intent(in)    :: current(:,:)
 real(dp),           intent(out)   :: ylm(:,:)
 real(dp),           intent(out)   :: ylmgr(:,:,:)
 logical, optional,  intent(in)    :: update_vecpot_ind

 !Local variables-------------------------------
 character(len=500) :: msg
 integer            :: i
 logical            :: lvecpot_ind
 real(dp)           :: tmp(3)

! ***********************************************************************

 ! Update potential vector if time-dependent electric field perturbation is present
 select case (tdef%ef_type)
   !No external field
   case (0)
   !Dirac pulse: vector potential is just an Heaviside function
   case (1)
      if (abs(time-tdef%ef_tzero)<tol16) then
         tdef%efield(:) = tdef%ef_ezero(:)
      else
         tdef%efield(:) = zero
      end if
      if (time >= tdef%ef_tzero) then
         tdef%vecpot_ext(:) = -tdef%ef_ezero(:)
      end if
   !Pulse with sin^2 shape:
   !E(t) = E0*cos(w*(t-t0))*sin^2(pi*(t-t0)/tau)
   !A(t) = -(E0/2w)*sin(w*(t-t0))+E0/(4*(2pi/tau+w))*sin((2pi/tau+w)*(t-t0))+E0/(4(2pi/taur-w))*sin((2pi/tau-w)*(t-t0))
!  case(2)
!     if (time >= tdef%ef_tzero+tdef%ef_tau) then
!        tdef%efield(:) = zero
!     else if (time >= tdef%ef_tzero) then
!        t = time-tdef%ef_tzero
!        tdef%efield(:) = tdef%ef_ezero*cos(tdef%ef_omega*t)*sin(pi*t/tdef%ef_tau)**2
!        tdef%vecpot_ext(:) = tdef%ef_ezero*(-sin(tdef%ef_omega*t)/(two*tdef%ef_omega) &
!                                          & +sin(tdef%ef_sin_a*t)/(four*tdef%ef_sin_a) &
!                                          & +sin(tdef%ef_sin_b*t)/(four*tdef%ef_sin_b))
!     end if
   case default
      write(msg,"(a)") "Unknown electric field type - check the value of td_ef_type"
      ABI_ERROR(msg)
 end select

 lvecpot_ind = .true.
 if (present(update_vecpot_ind)) lvecpot_ind = update_vecpot_ind

 !Induced vector potential
 !Should deal with sppol?! How?
 !d^2A_ind/dt^2 = 4piJ(t)
 !A_ind(t+dt) = 2*A_ind(t) - A_ind(t-dt) + 4*pi*dt**2*J(t)
 if (tdef%induced_vecpot) then
    if (lvecpot_ind) then
      tmp = tdef%vecpot_ind(:,2) ! t - 2dt
      tdef%vecpot_ind(:,2) = tdef%vecpot_ind(:,1) ! t - dt
      tdef%vecpot_ind(:,1) = 2*tdef%vecpot_ind(:,2) - tmp(:) + four*pi*(dtset%dtele**2)*current(:,1) ! t
    end if
    tdef%vecpot = tdef%vecpot_ext + tdef%vecpot_ind(:,1)
 else
   tdef%vecpot = tdef%vecpot_ext
 end if

 tdef%vecpot_red = matmul(transpose(rprimd),tdef%vecpot)

 if (tdef%ef_type /= 0) then
   !Update the k+A grid used in cprojs
   !Divide by 2pi here seems necessary
   do i = 1,3
      tdef%kpa(i,:) = dtset%kptns(i,:) + tdef%vecpot_red(i)/(two*pi)
   end do
   ! update the spherical harmonics (computed at k+G+A)
   call initylmg(gprimd,kg,tdef%kpa,dtset%mkmem,mpi_enreg,mpsang,dtset%mpw,dtset%nband, &
               & dtset%nkpt,npwarr,dtset%nsppol,0,rprimd,ylm,ylmgr)
 end if

end subroutine tdef_update
!!***

!!****f* m_rttddft/tdef_restart
!!
!! NAME
!!  tdef_restart
!!
!! FUNCTION
!!  Update some values for restart of calculation with TD electric field
!!  Essentially needed if we account for induced vector potential
!!
!! INPUTS
!!  [tdef = tdef structure to update]
!!  mpi_enreg = MPI communicators structure
!!  restart_unit = unit of restart file to read
!!
!! OUTPUT
!!  [tdef = updated tdef structure]
!!
!! SOURCE
subroutine tdef_restart(tdef,mpi_enreg,restart_unit)

 !Arguments ------------------------------------
 !scalars
 class(tdef_type),   intent(inout) :: tdef
 type(MPI_type),     intent(inout) :: mpi_enreg
 integer,            intent(in)    :: restart_unit

 !Local variables-------------------------------
 integer :: ierr

! ***********************************************************************

 if (tdef%induced_vecpot) then
      if (mpi_enreg%me == 0) then
         read(restart_unit,*) tdef%vecpot_ind(:,2)
         read(restart_unit,*) tdef%vecpot_ind(:,1)
      end if
 end if
 !Send to all procs
 call xmpi_bcast(tdef%vecpot_ind,0,mpi_enreg%comm_world,ierr)

end subroutine tdef_restart
!!***

end module m_rttddft_tdef
!!***
