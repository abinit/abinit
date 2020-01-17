!!****m* ABINIT/glue_lotf
!! NAME
!! glue_lotf
!!
!! FUNCTION
!!  Contains the GLUE procedure and parameters for Lotf
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2019 ABINIT group (MMancini)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Parameters for the glue potential
!!  ref: Ercolessi et al, Phil Mag A 58(1), 213 (1988)
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

module glue_lotf

 use defs_basis
 use m_errors

 implicit none
 private

 !--Parameters for the glue potential (Ercolessi et al, Phil Mag A 58(1), 213 (1988))

 !--Function phi_r
 real(dp),public :: dphi,rcphi
 real(dp) :: a0I,a1I,a2I,a3I,a4I
 real(dp) :: a0II,a1II,a2II,a3II,a4II,a5II,a6II

 !--Function rho_r
 real(dp),public :: rmrho
 real(dp) :: rbrho
 real(dp) :: b0I,b1I,b2I,b3I
 real(dp) :: b0II,b1II,b2II,b3II,b0III,b1III,b2III,b3III

 !--Function U_n
 real(dp) :: n0U,nsU
 real(dp) :: c0I,c1I,c2I,c3I,c4I
 real(dp) :: c0II,c1II,c2II,c3II,c4II,c0III,c1III,c2III,c3III

 !--Auxiliary variables to change glue parameters
 real(dp),parameter :: ddphi = 0.028779246_dp
 real(dp),parameter :: scalefactor_phi = one


 public ::             &
   glue_init,          &
   glue_pair_devs,     &
   glue_pair,          &
   calc_coord,         &
   calc_rhop,          &
   calc_coord_new_d,   &
   rho_devs,           &
   rhop_value,         &
   eval_U_n,           &
   eval_Upp_n

contains !===========================================================
 !!***

!!****f* glue_lotf/glue_init
!! NAME
!! glue_init
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine glue_init()

  !--Function phi_r
! *************************************************************************
  dphi =  (0.2878207442141723d+01 + ddphi) / 0.529177d0
  rcphi = (0.3700000000000000d+01 + ddphi) / 0.529177d0

  a0I = -0.8000000000000000d-01 / (13.6057d0 * 2.d0) *&
    scalefactor_phi
  a1I =  0.0000000000000000d+00 / (13.6057d0 * 2.d0) * (0.529177d0) *&
    scalefactor_phi
  a2I =  0.7619231375231362d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**2 *&
    scalefactor_phi
  a3I = -0.8333333333333333d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**3 *&
    scalefactor_phi
  a4I = -0.1211483464993159d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**4 *&
    scalefactor_phi

  !        a0II = -0.8000000000000000d-01 / (13.6057d0 * 2.d0)
  a0II = a0I
  !        a1II =  0.0000000000000000d+00 / (13.6057d0 * 2.d0) * (0.529177d0)
  a1II = a1I
  !        a2II =  0.7619231375231362d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**2
  a2II = a2I
  a3II = -0.8333333333333333d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**3 *&
    scalefactor_phi
  a4II = -0.1096009851140349d+01 / (13.6057d0 * 2.d0) * (0.529177d0)**4 *&
    scalefactor_phi
  a5II =  0.2158417178555998d+01 / (13.6057d0 * 2.d0) * (0.529177d0)**5 *&
    scalefactor_phi
  a6II = -0.9128915709636862d+00 / (13.6057d0 * 2.d0) * (0.529177d0)**6 *&
    scalefactor_phi

  !--Function rho_r
  rbrho = 0.3500000000000000d+01 / 0.529177d0
  rmrho = (0.3900000000000000d+01 + ddphi) / 0.529177d0

  b0I =  0.1000000000000000d+01
  b1I = -0.6800000000000000d+00 * (0.529177d0)
  b2I =  0.7500000000000000d+00 * (0.529177d0)**2
  b3I = -0.1333333333333333d+01 * (0.529177d0)**3

  !        b0II =  0.1000000000000000d+01
  b0II = b0I
  !        b1II = -0.6800000000000000d+00 * (0.529177d0)
  b1II = b1I
  !        b2II =  0.7500000000000000d+00 * (0.529177d0)**2
  b2II = b2I
  b3II = -0.1527241171296038d+01 * (0.529177d0)**3

  b0III =  0.0000000000000000d+00
  b1III =  0.0000000000000000d+00 * (0.529177d0)
  b2III =  0.5578188675490974d+01 * (0.529177d0)**2
  b3III =  0.6132971688727435d+01 * (0.529177d0)**3

  !--Function U_n
  n0U =  0.1200000000000000d+02
  nsU =  0.9358157767784574d+01

  c0I = -0.2793388616771698d+01 / (13.6057d0 * 2.d0) * scalefactor_phi
  c1I = -0.3419999999999999d+00 / (13.6057d0 * 2.d0) * scalefactor_phi
  c2I =  0.3902327808424106d-01 / (13.6057d0 * 2.d0) * scalefactor_phi
  c3I =  0.7558829951858879d-02 / (13.6057d0 * 2.d0) * scalefactor_phi
  c4I =  0.3090472511796849d-03 / (13.6057d0 * 2.d0) * scalefactor_phi

  c0II = -0.3300000000000000d+01 / (13.6057d0 * 2.d0) * scalefactor_phi
  c1II =  0.0000000000000000d+00 / (13.6057d0 * 2.d0) * scalefactor_phi
  c2II =  0.8618226772941980d-01 / (13.6057d0 * 2.d0) * scalefactor_phi
  c3II =  0.4341701445034724d-02 / (13.6057d0 * 2.d0) * scalefactor_phi
  c4II = -0.3044398779375916d-03 / (13.6057d0 * 2.d0) * scalefactor_phi

  !        c0III = -0.3300000000000000d+01 / (13.6057d0 * 2.d0)
  c0III = c0II
  !        c1III =  0.0000000000000000d+00 / (13.6057d0 * 2.d0)
  c1III = c1II
  !        c2III =  0.8618226772941980d-01 / (13.6057d0 * 2.d0)
  c2III = c2II
  c3III =  0.4325981467602070d-02 / (13.6057d0 * 2.d0) * scalefactor_phi

 end subroutine glue_init
 !!***


!!****f* glue_lotf/glue_pair_devs
!! NAME
!! glue_pair_devs
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine glue_pair_devs(alpha_dum,RD,r_au,epot_2,fdum,dfdum)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) ::  alpha_dum(3),RD(3)
  real(dp),intent(in) ::  r_au
  real(dp),intent(out) ::  epot_2
  real(dp),intent(out) ::  fdum(3),dfdum(3,2)
  !Local ---------------------------
  real(dp) :: rst, r, r2, r3, r4, r5, r6
  real(dp) :: drcphi, dphi2, rcphi2, scale_factor
  ! *********************************************************************

!    --The calculation for a given pair
     drcphi = rcphi - dphi
     dphi2 = alpha_dum(1)

!    --Rewritten
     rcphi2 = dphi2 + drcphi

     scale_factor = alpha_dum(2)

     rst = sqrt(r_au)

     r = rst - dphi2
     r2 = r*r
     r3 = r2*r
     r4 = r3*r

!    --Energy
     epot_2 = zero

     if (rst < dphi2) then
       epot_2 = scale_factor * (a4I*r4 + a3I*r3 + a2I*r2 + a1I*r + a0I)
     elseif (rst < rcphi2) then
       r5 = r4*r
       r6 = r5*r
       epot_2 = scale_factor * (a6II*r6 + a5II*r5 + a4II*r4 + a3II*r3 + a2II*r2 +&
       a1II*r + a0II)
     end if

!    --Forces (without the minus sign, because we use rji) and derivatives w.r.t. parameter d
     fdum(:) = zero
     dfdum(:,:) = zero

     if (rst < dphi2) then
       fdum(:) = scale_factor * (4.d0*a4I*r3 + 3.d0*a3I*r2 + 2.d0*a2I*r + a1I) * RD(:)/rst
       dfdum(:,1) = -scale_factor * (12.d0*a4I*r2 + 6.d0*a3I*r + 2.d0*a2I) * RD(:)/rst
       dfdum(:,2) = fdum(:) / scale_factor

     elseif (rst < rcphi2) then
       fdum(:) = (scale_factor) * &
&       (6.d0*a6II*r5 + 5.d0*a5II*r4 + 4.d0*a4II*r3 + 3.d0*a3II*r2 + 2.d0*a2II*r + a1II) * &
&       RD(:)/rst

       dfdum(:,1) = -scale_factor * &
&       (30.d0*a6II*r4 + 20.d0*a5II*r3 + 12.d0*a4II*r2 + 6.d0*a3II*r + 2.d0*a2II) * RD(:)/rst

       dfdum(:,2) = fdum(:) / scale_factor

     end if
   end subroutine glue_pair_devs
!!***


!!****f* glue_lotf/glue_pair
!! NAME
!! glue_pair
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine glue_pair(RD,r_au,epot_2,fdum)

  implicit none

  !Arguments ------------------------
  real(dp) ::  RD(3), fdum(3), r_au, epot_2
  !Local ---------------------------
  real(dp) ::  rst, r, r2, r3, r4, r5, r6
  integer ::  k
  ! *********************************************************************
!    --The calculation for a given pair
     rst = sqrt(r_au)
     r = rst - dphi
     r2 = r*r
     r3 = r2*r
     r4 = r3*r

!    --Energy
     epot_2 = zero
     if (rst < dphi) then
       epot_2 = a4I*r4 + a3I*r3 + a2I*r2 + a1I*r + a0I
     elseif (rst < rcphi) then
       r5 = r4*r
       r6 = r5*r
       epot_2 = a6II*r6 + a5II*r5 + a4II*r4 + a3II*r3 + a2II*r2 +&
       a1II*r + a0II
     end if

!    --Forces (without the minus sign, because we use rji)
     fdum(:) = zero

     do k =1,3
       if (rst < dphi) then
         fdum(k) =  ( 4.d0*a4I*r3 + 3.d0*a3I*r2 + 2.d0*a2I*r + a1I ) * &
         RD(k)/rst
       elseif (rst < rcphi) then
         fdum(k) =  ( 6.d0*a6II*r5 + 5.d0*a5II*r4 + 4.d0*a4II*r3 + &
         3.d0*a3II*r2 + 2.d0*a2II*r + a1II ) * &
         RD(k)/rst
       end if
     end do

   end subroutine glue_pair
!!***

!!****f* glue_lotf/calc_coord
!! NAME
!! calc_coord
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine calc_coord(r_au,coordatom_dum)

  implicit none

  !Arguments ------------------------
  real(dp) :: r_au,coordatom_dum
  !Local ---------------------------
  real(dp) :: rst,r,rho_neigh
  ! *********************************************************************
!    --Adds the density of a single given neighbour

     rst = sqrt(r_au)
     rho_neigh = zero

     if (rst < dphi) then
       r = rst - dphi
       rho_neigh = b3I*r*r*r + b2I*r*r + b1I*r + b0I
     elseif (rst < rbrho) then
       r = rst - dphi
       rho_neigh = b3II*r*r*r + b2II*r*r + b1II*r + b0II
     elseif (rst < rmrho) then
       r = rst - rmrho
       rho_neigh = b3III*r*r*r + b2III*r*r + b1III*r + b0III
     end if

     coordatom_dum = coordatom_dum + rho_neigh

   end subroutine calc_coord
!!***

!!****f* glue_lotf/calc_rhop
!! NAME
!! calc_rhop
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!!
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

     subroutine calc_rhop(r_st,rhop_dum)

     implicit none

!    Arguments ------------------------
     real(dp),intent(in) :: r_st
     real(dp),intent(out) :: rhop_dum
!    Local ---------------------------
     real(dp) :: r,r2
!    *********************************************************************

     rhop_dum = zero

     if (r_st < dphi) then
       r = r_st - dphi
       r2 = r*r
       rhop_dum = (3.d0*b3I*r2 + 2.d0*b2I*r + b1I)
     elseif (r_st < rbrho) then
       r = r_st - dphi
       r2 = r*r
       rhop_dum = (3.d0*b3II*r2 + 2.d0*b2II*r + b1II)
     elseif (r_st < rmrho) then
       r = r_st - rmrho
       r2 = r*r
       rhop_dum = (3.d0*b3III*r2 + 2.d0*b2III*r + b1III)
     end if
   end subroutine calc_rhop
!!***

!!****f* glue_lotf/calc_coord_new_d
!! NAME
!! calc_coord_new_d
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine calc_coord_new_d(r_au,alpha_d,coordatom_dum)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: r_au,alpha_d
  real(dp),intent(out) ::coordatom_dum
  !Local ---------------------------
  real(dp) :: rst,r,rho_neigh
  ! *********************************************************************

!    --Adds the density of a single given neighbour
     rst = sqrt(r_au)

     rho_neigh = zero

     if (rst < alpha_d) then
       r = rst - alpha_d
       rho_neigh = b3I*r*r*r + b2I*r*r + b1I*r + b0I
     elseif (rst < rbrho) then
       r = rst - alpha_d
       rho_neigh = b3II*r*r*r + b2II*r*r + b1II*r + b0II

!      --Restriction on rmrho to keep rho and rhop continuous
     elseif (rst < (rmrho+alpha_d-dphi)) then
       r = rst - (rmrho+alpha_d-dphi)
       rho_neigh = b3III*r*r*r + b2III*r*r + b1III*r + b0III
     end if

     coordatom_dum = coordatom_dum + rho_neigh
   end subroutine calc_coord_new_d
!!***

!!****f* glue_lotf/rho_devs
!! NAME
!! rho_devs
!!
!! FUNCTION
!!
!! INPUTS
!!
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

     subroutine rho_devs(r_au,alpha_d,rho_neigh_d,rho_neigh_p,rho_neigh_pd)

     real(dp),intent(in) :: r_au,alpha_d
     real(dp),intent(out) :: rho_neigh_d,rho_neigh_p,rho_neigh_pd
!    Local ---------------------------
     real(dp) :: rst,r
!    *********************************************************************

!    --Adds the density of a single given neighbour
     rst = sqrt(r_au)

     rho_neigh_p = zero
     rho_neigh_d = zero
     rho_neigh_pd = zero

     if (rst < alpha_d) then
       r = rst - alpha_d
       rho_neigh_d = - ( 3.d0*b3I*r*r + 2.d0*b2I*r + b1I )
       rho_neigh_p = - rho_neigh_d
       rho_neigh_pd = - ( 6.d0*b3I*r + 2.d0*b2I )

     elseif (rst < rbrho) then
       r = rst - alpha_d
       rho_neigh_d = - ( 3.d0*b3II*r*r + 2.d0*b2II*r + b1II )
       rho_neigh_p = - rho_neigh_d
       rho_neigh_pd = - ( 6.d0*b3II*r + 2.d0*b2II )

!      --Restriction on rmrho to keep rho and rhop continuous
     elseif (rst < (rmrho+alpha_d-dphi)) then
       r = rst - (rmrho+alpha_d-dphi)
       rho_neigh_d = zero
       rho_neigh_p = 3.d0*b3III*r*r + 2.d0*b2III*r + b1III
       rho_neigh_pd = zero

     end if
   end subroutine rho_devs
!!***

!!****f* glue_lotf/rhop_value
!! NAME
!! rhop_value
!!
!! FUNCTION
!!
!!
!! INPUTS
!!
!! PARENTS
!!      m_eval_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine rhop_value(rst,alpha_d,rhop_dum)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: rst,alpha_d
  real(dp),intent(out) :: rhop_dum
  !Local ---------------------------
  real(dp) :: r
  ! *********************************************************************

!    --Adds the density of a single given neighbour
     rhop_dum = zero

     if (rst < alpha_d) then
       r = rst - alpha_d
       rhop_dum = 3.d0*b3I*r*r + 2.d0*b2I*r + b1I
     elseif (rst < rbrho) then
       r = rst - alpha_d
       rhop_dum = 3.d0*b3II*r*r + 2.d0*b2II*r + b1II
!      --Restriction on rmrho to keep rho and rhop continuous
     elseif (rst < (rmrho+alpha_d-dphi)) then
       r = rst - (rmrho+alpha_d-dphi)
       rhop_dum = 3.d0*b3III*r*r + 2.d0*b2III*r + b1III
     end if
   end subroutine rhop_value
!!***


!!****f* glue_lotf/eval_U_n
!! NAME
!! eval_U_n
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine eval_U_n(coordatom_i,epot_dum2,up_dum)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: coordatom_i
  real(dp),intent(out) :: epot_dum2,up_dum
  !Local ---------------------------
  real(dp) :: cns,cns2,cns3,cns4,cn0,cn02,cn03,cn04
  character(len=500) :: msg

  !--Evaluate U and Up for atom i
  epot_dum2 = zero
  up_dum = zero

  if (coordatom_i < zero) then

    write(msg,'(3a,i8,2a)')&
      &    'ERROR: negative coordination!!! ',ch10,&
      &    'coord = ', coordatom_i,ch10,&
      &    'The coordination cannot be negative!'
    MSG_ERROR(msg)

  elseif (coordatom_i < nsU) then

    cns = coordatom_i - nsU
    cns2 = cns * cns
    cns3 = cns2 * cns
    cns4 = cns3 * cns

    epot_dum2 = c4I*cns4 + c3I*cns3 + c2I*cns2 + c1I*cns + c0I
    up_dum = 4.d0*c4I*cns3 + 3.d0*c3I*cns2 + 2.d0*c2I*cns + c1I

  elseif (coordatom_i < n0U) then

    cn0 = coordatom_i - n0U
    cn02 = cn0 * cn0
    cn03 = cn02 * cn0
    cn04 = cn03 * cn0

    epot_dum2 = c4II*cn04 + c3II*cn03 + c2II*cn02 + c1II*cn0 + c0II
    up_dum = 4.d0*c4II*cn03 + 3.d0*c3II*cn02 + 2.d0*c2II*cn0 + c1II

  elseif (coordatom_i >= n0U) then

    cn0 = coordatom_i - n0U
    cn02 = cn0 * cn0
    cn03 = cn02 * cn0

    epot_dum2 = c3III*cn03 + c2III*cn02 + c1III*cn0 + c0III
    up_dum = 3.d0*c3III*cn02 + 2.d0*c2III*cn0 + c1III

  endif
 end subroutine eval_U_n
 !!***


!!****f* eval_lotf/eval_Upp_n
!! NAME
!! eval_Upp_n
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine eval_Upp_n(coordatom_i,up_dum,upp_dum)

  implicit none

  !Arguments ------------------------
  real(dp),intent(in) :: coordatom_i
  real(dp),intent(out) :: up_dum,upp_dum
  !Local ---------------------------
  real(dp) :: cns,cns2,cns3,cn0,cn02,cn03
  character(len=500) :: msg

! *************************************************************************

  !--Evaluate U and Up for atom i
  up_dum = zero
  upp_dum = zero

  if (coordatom_i < zero) then

    write(msg,'(3a,i8,2a)')&
      &    'ERROR: negative coordination!!! ',ch10,&
      &    'coord = ', coordatom_i,ch10,&
      &    'The coordination cannot be negative!'
    MSG_ERROR(msg)

  elseif (coordatom_i < nsU) then

    cns = coordatom_i - nsU
    cns2 = cns * cns
    cns3 = cns2 * cns

    up_dum = 4.d0*c4I*cns3 + 3.d0*c3I*cns2 + 2.d0*c2I*cns + c1I
    upp_dum = 12.d0*c4I*cns2 + 6.d0*c3I*cns + 2.d0*c2I

  elseif (coordatom_i < n0U) then

    cn0 = coordatom_i - n0U
    cn02 = cn0 * cn0
    cn03 = cn02 * cn0

    up_dum = 4.d0*c4II*cn03 + 3.d0*c3II*cn02 + 2.d0*c2II*cn0 + c1II
    upp_dum = 12.d0*c4II*cn02 + 6.d0*c3II*cn0 + 2.d0*c2II

  elseif (coordatom_i >= n0U) then

    cn0 = coordatom_i - n0U
    cn02 = cn0 * cn0

    up_dum = 3.d0*c3III*cn02 + 2.d0*c2III*cn0 + c1III
    upp_dum = 6.d0*c3III*cn0 + 2.d0*c2III

  endif
 end subroutine eval_Upp_n

end module glue_lotf
!!***
