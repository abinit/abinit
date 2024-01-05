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
 use defs_abitypes,  only: MPI_type

 use m_dtset,        only: dataset_type
 use m_errors,       only: msg_hndl
 use m_initylmg,     only: initylmg

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
   real(dp) :: vecpot(3)     !associated vector potential (in velocity gauge)
   real(dp) :: vecpot_red(3) !vector potential in reduced coord.
   real(dp) :: ef_tzero      !time at which elec field is switched on
   real(dp) :: ef_omega      !angular freq. of TD elec field
   real(dp) :: ef_tau        !time width of the pulse
   real(dp) :: ef_sin_a      !useful constant for sin^2 pulse
   real(dp) :: ef_sin_b      !useful constant for sin^2 pulse
   real(dp),allocatable :: kpa(:,:) ! contains kpts + A 
                                    ! (in reduced coordinates in reciprocal space)

   contains

   procedure :: init => tdef_init
   procedure :: update => tdef_update

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
!!
!! OUTPUT
!!  [tdef = updated tdef structure]
!!
!! SOURCE
subroutine tdef_init(tdef, td_ef_type, td_ef_pol, td_ef_ezero, td_ef_tzero, td_ef_lambda, td_ef_tau, nkpt, kpts)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdef_type), intent(inout) :: tdef
 integer,          intent(in)    :: td_ef_type, nkpt
 real(dp),         intent(in)    :: td_ef_ezero, td_ef_tzero, td_ef_lambda, td_ef_tau
 !arrays
 real(dp),         intent(in)    :: kpts(:,:)

 !Local variables-------------------------------
 real(dp),         intent(in)    :: td_ef_pol(3)

! ***********************************************************************

 
 tdef%ef_type = td_ef_type
 tdef%ef_ezero(:) = td_ef_pol(:)*td_ef_ezero
 tdef%ef_tau      = td_ef_tau
 tdef%ef_omega    = 2.0_dp*pi*Sp_Lt/td_ef_lambda !2*pi*f=2*pi*c/lambda
 tdef%ef_tzero    = td_ef_tzero
 tdef%ef_sin_a    = 2.0_dp*pi/td_ef_tau + tdef%ef_omega
 tdef%ef_sin_b    = 2.0_dp*pi/td_ef_tau - tdef%ef_omega
 
 tdef%efield = 0.0_dp
 tdef%vecpot = 0.0_dp
 tdef%vecpot_red = 0.0_dp
 
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
!!  time = propagation time
!!  rprimd = primitive vectors
!!
!! OUTPUT
!!  [tdef = updated tdef structure]
!!
!! SOURCE
subroutine tdef_update(tdef, dtset, mpi_enreg, time, rprimd, gprimd, kg, mpsang, npwarr, ylm, ylmgr)

 implicit none

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
 real(dp),           intent(out)   :: ylm(:,:)
 real(dp),           intent(out)   :: ylmgr(:,:,:)

 !Local variables-------------------------------
 character(len=500) :: msg
 real(dp)           :: t
 integer            :: i

! ***********************************************************************

 ! Update potential vector if time-dependent electric field perturbation is present
 select case (tdef%ef_type)
   !No external field
   case (0)
   !Dirac pulse: vector potential is just an Heaviside function
   case (1)
      if (time >= tdef%ef_tzero) then
         tdef%efield(:) = tdef%ef_ezero(:)
         tdef%vecpot(:) = -tdef%ef_ezero(:)
      end if
      tdef%vecpot_red = matmul(transpose(rprimd),tdef%vecpot)/(2.0_dp*pi)
   !Pulse with sin^2 shape:
   !E(t) = E0*cos(w*(t-t0))*sin^2(pi*(t-t0)/tau)
   !A(t) = -(E0/2w)*sin(w*(t-t0))+E0/(4*(2pi/tau+w))*sin((2pi/tau+w)*(t-t0))+E0/(4(2pi/taur-w))*sin((2pi/tau-w)*(t-t0))
   case(2)
      if (time >= tdef%ef_tzero+tdef%ef_tau) then
         tdef%efield(:) = 0.0_dp
      else if (time >= tdef%ef_tzero) then
         t = time-tdef%ef_tzero
         tdef%efield(:) = tdef%ef_ezero*cos(tdef%ef_omega*t)*sin(pi*t/tdef%ef_tau)**2
         tdef%vecpot(:) = tdef%ef_ezero*(-sin(tdef%ef_omega*t)/(2.0_dp*tdef%ef_omega) &
                                       & +sin(tdef%ef_sin_a*t)/(4.0_dp*tdef%ef_sin_a) &
                                       & +sin(tdef%ef_sin_b*t)/(4.0_dp*tdef%ef_sin_b))
      end if
      tdef%vecpot_red(:) = matmul(transpose(rprimd),tdef%vecpot)/(2.0_dp*pi)
   case default
      write(msg,"(a)") "Unknown electric field type - check the value of td_ef_type"
      ABI_ERROR(msg)
 end select
 if (tdef%ef_type /= 0) then
   !Update the k+A grid used in cprojs
   do i = 1,3
      tdef%kpa(i,:) = dtset%kptns(i,:) + tdef%vecpot_red(i)
   end do
      ! update the spherical harmonics (computed at k+G+A)
      call initylmg(gprimd,kg,tdef%kpa,dtset%mkmem,mpi_enreg,mpsang,dtset%mpw,dtset%nband, &
                  & dtset%nkpt,npwarr,dtset%nsppol,0,rprimd,ylm,ylmgr)
 end if

end subroutine tdef_update
!!***

end module m_rttddft_tdef
!!***
