!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_fs_bands
!!
!! NAME
!! get_fs_bands
!!
!! FUNCTION
!! This routine determines the bands which contribute to the Fermi surface
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  eigenGS = ground state eigenvalues
!!  hdr = header from input GS file
!!  ep_b_min, ep_b_max=A non-zero value is used to impose certain bands.
!!  fermie=Fermi level.
!!  eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)=Energies.
!!
!! OUTPUT
!!  minFSband,maxFSband=Minimun and maximum index for the bands that cross the Fermi level
!!  nkptirr=Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!
!! TODO
!!  1) Indeces and dimensions should should be spin dependent.
!!  2) In the present status of the code, all the k-points in the IBZ are used!
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_fs_bands(eigenGS,hdr,fermie,ep_b_min,ep_b_max,minFSband,maxFSband,nkptirr)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_fs_bands'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ep_b_min, ep_b_max
 integer,intent(out) :: minFSband,maxFSband,nkptirr
 real(dp),intent(in) :: fermie
 type(hdr_type),intent(in) :: hdr
!arrays
 real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isppol,nband
 real(dp) :: epsFS,gausstol,gaussig
 character(len=500) :: message
 integer :: kpt_phonflag(hdr%nkpt)

! *************************************************************************

!supposes nband is equal for all kpts
 nband = hdr%nband(1)

!gausstol = minimum weight value for integration weights on FS
!should be set to reproduce DOS at Ef (Ref. PRB 34, 5065 p. 5067)
 gausstol = 1.0d-10

!use same band indices in both spin channels
 maxFSband=1
 minFSband=nband

!window of states around fermi Energy is contained in +/- epsFS
!should be adjusted to take into account a minimal but sufficient
!fraction of the kpoints: see the loop below.
!The 1000 is purely empirical!!!
!Should also take into account the density of kpoints.
!gaussig = width of gaussian energy window around fermi energy
!needed to get a good fraction of kpoints contributing to the FS
 
 gaussig = (maxval(eigenGS)-minval(eigenGS))/1000.0_dp
 
 write (message,'(a,f11.8,2a)')' get_fs_bands : initial energy window = ',gaussig,ch10,&
& ' The window energy will be increased until the full k-grid is inside the range'
 call wrtout(std_out,message,'COLL')

!NOTE: could loop back to here and change gaussig until we have
!a certain fraction of the kpoints in the FS region...
 nkptirr = 0
 
!Do not use restricted fermi surface: include all kpts -> one
 do while (nkptirr < hdr%nkpt)
   gaussig = gaussig*1.05_dp

!  we must take into account kpoints with states within epsFS:
   epsFS = gaussig*sqrt(log(one/(gaussig*sqrt(pi)*gausstol)))
   
!  check if there are eigenvalues close to the Fermi surface
!  (less than epsFS from it)
   kpt_phonflag(:) = 0

!  do for each sppol channel
   do isppol=1,hdr%nsppol
     do ikpt=1,hdr%nkpt
       do iband=1,nband
         if (abs(eigenGS(iband,ikpt,isppol) - fermie) < epsFS) then
           kpt_phonflag(ikpt) = 1
           if (iband > maxFSband) maxFSband = iband
           if (iband < minFSband) minFSband = iband
         end if
       end do
     end do
   end do ! isppol
   
!  if user imposed certain bands for e-p, make sure they are kept
   if (ep_b_min /= 0 .and. ep_b_min < minFSband) then
     minFSband = ep_b_min
   end if
   if (ep_b_max /= 0 .and. ep_b_max > maxFSband) then
     maxFSband = ep_b_max
   end if
   
!  number of irreducible kpoints (by all sym) contributing to the Fermi surface (to be completed by symops).
   nkptirr = sum(kpt_phonflag(:))
 end do

 write(std_out,*) ' Energy window around Fermi level= ',epsFS,' nkptirr= ',nkptirr

end subroutine get_fs_bands
!!***
