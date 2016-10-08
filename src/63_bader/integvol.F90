!{\src2tex{textfont=tt}}
!!****f* ABINIT/integvol
!! NAME
!! integvol
!!
!! FUNCTION
!! This routine integrates the volume of the Bader atom
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (see side effects)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works on the data contained in the aimfields and aimprom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      coeffs_gausslegint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine integvol()

 use defs_basis
 use defs_aimfields
 use defs_aimprom
 use defs_parameters
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,   only : coeffs_gausslegint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integvol'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables ------------------------------
!scalars
 integer :: batom,ii,jj,nph,nth
 real(dp) :: chgint,ct1,ct2,nsphe,phimax,phimin
 real(dp) :: rsmax,rsmin,themax,themin
 logical :: gaus,weit
!arrays
 real(dp) :: shift(3)
 real(dp),allocatable :: rdint(:,:)
 real(dp),allocatable :: wgrs(:,:)

! *********************************************************************

 tpi=two_pi
 gaus=.true.
 weit=.true.


 rewind(unts)
 read(unts,*) batom,shift
 read(unts,*) nth,themin,themax
 read(unts,*) nph,phimin,phimax

 write(std_out,*) 'NTH NPH ',nth,nph

 ABI_ALLOCATE(wgrs,(nth,nph))
 ABI_ALLOCATE(rdint,(nth,nph))

 do ii=1,nth
   do jj=1,nph
     if (weit) then
       read(unts,*) th(ii),ph(jj),rs(ii,jj),wgrs(ii,jj)
     else
       read(unts,*) th(ii),ph(jj),rs(ii,jj)
     end if
   end do
 end do
 read(unts,*) rsmin,rsmax


 if (gaus) then
   ct1=cos(themin)
   ct2=cos(themax)
   call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
   call coeffs_gausslegint(phimin,phimax,ph,wph,nph)
 end if

 do ii=1,nth
   do jj=1,nph
     if (.not.weit) then
       if (gaus) then
         wgrs(ii,jj)=wcth(ii)*wph(jj)
       else
         wgrs(ii,jj)=1._dp
       end if
     end if
   end do
 end do

 nsphe=0._dp
 do ii=1,nth
   do jj=1,nph
     nsphe=nsphe+rs(ii,jj)**3/3._dp*wgrs(ii,jj)
   end do
 end do
 if (gaus.or.weit) then
   nsphe=nsphe*(pi/(themin-themax))*(tpi/(phimax-phimin))
 else
   nsphe=nsphe/(nth*nph)*2.0*tpi
 end if
 chgint=nsphe

 write(std_out,*) ':VOLTOT ',batom,chgint
 write(untout,'("Volume of the Bader atom: ", I6, F16.8)') batom,chgint

end subroutine integvol
!!***
