!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_all_gkk2
!!
!! NAME
!! get_all_gkk2
!!
!! FUNCTION
!! This routine determines where to store gkk2 matrix elements (disk or RAM)
!!   and calls interpolate_gkk to calculate them.
!!   This is the most time consuming step.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   acell = lengths of unit cell vectors
!!   amu = masses of atoms
!!   atmfrc = atomic force constants
!!   dielt = dielectric tensor
!!   dipdip = dipole-dipole contribution flag
!!   dyewq0 =
!!   elph_ds = datastructure for elphon data and dimensions
!!   kptirr_phon = irreducible set of fermi-surface kpoints
!!   kpt_phon = full set of fermi-surface kpoints
!!   ftwghtgkk = weights for FT of matrix elements
!!   gmet = metric in reciprocal space
!!   indsym = indirect mapping of atoms under symops
!!   mpert = maximum number of perturbations
!!   msym = maximum number of symmetries (usually nsym)
!!   nsym = number of symmetries
!!   ntypat = number of types of atoms
!!   onegkksize = size of one gkk record, in bytes
!!   rmet = real-space metric
!!   rprim = unit cell lattice vectors (dimensionless)
!!   rprimd = real-space unit-cell lattice vectors
!!   rpt = points in real space for FT, in canonical coordinates
!!   symrel = symmetry operations in reduced real space
!!   trans = Atomic translations : xred = rcan + trans
!!   typat = array of types of atoms
!!   ucvol = unit cell volume
!!   xred = reduced coordinates of atoms
!!   zeff = Born effective charges
!!
!! OUTPUT
!!   elph_ds = calculated |gkk|^2 are in elph_ds%gkk2
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      interpolate_gkk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine get_all_gkk2(crystal,ifc,elph_ds,kptirr_phon,kpt_phon)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

 use m_crystal,   only : crystal_t
 use m_ifc,       only : ifc_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_all_gkk2'
 use interfaces_77_ddb, except_this_one => get_all_gkk2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
 real(dp),intent(in) :: kptirr_phon(3,elph_ds%k_phon%nkptirr)

!Local variables-------------------------------
!scalars
 integer :: iost,onediaggkksize,sz1,sz2,sz3,sz4
 real(dp) :: realdp_ex
 !character(len=500) :: msg

! *************************************************************************

 if (elph_ds%nsppol /= 1) then
   MSG_ERROR('get_all_gkk2: nsppol>1 not coded yet!')
 end if

 onediaggkksize = elph_ds%nbranch*elph_ds%k_phon%nkpt*kind(realdp_ex)

 elph_ds%unit_gkk2 = 37
 if (elph_ds%gkk2write == 0) then
   write(std_out,*) 'get_all_gkk2 : keep gkk2 in memory. Size = ',&
&   4.0*dble(elph_ds%k_phon%nkpt)*dble(onediaggkksize)/&
&   1024.0_dp/1024.0_dp, " Mb"
   sz1=elph_ds%nbranch
   sz2=elph_ds%ngkkband
   sz3=elph_ds%ngkkband
   sz4=elph_ds%k_phon%nkpt
   ABI_ALLOCATE(elph_ds%gkk2,(sz1,sz2,sz3,sz4,elph_ds%k_phon%nkpt,1))
   elph_ds%gkk2(:,:,:,:,:,:) = zero

 else if (elph_ds%gkk2write == 1) then
   write(std_out,*) 'get_all_gkk2 : About to open gkk2 file : '
   write(std_out,*) elph_ds%unit_gkk2,onediaggkksize
   open (unit=elph_ds%unit_gkk2,file='gkk2file',access='direct',&
&   recl=onediaggkksize,form='unformatted',status='new',iostat=iost)
   if (iost /= 0) then
     MSG_ERROR('error opening gkk2file as new')
   end if
!  rewind (elph_ds%unit_gkk2)
   write(std_out,*) 'get_all_gkk2 : disk file with gkk^2 created'
   write(std_out,*) '  calculate from real space gkk and phonon modes'
   write(std_out,*) '  gkk2write = 1 is forced: can take a lot of time! '
   write(std_out,*) ' size = ', 4.0*dble(onediaggkksize)*dble(elph_ds%k_phon%nkpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'
 else
   MSG_ERROR('bad value of gkk2write')
 end if

!here do the actual calculation of |g_kk|^2
 MSG_ERROR("MGNOTE: interpolate_gkk is broken")
 ABI_UNUSED(kptirr_phon(1,1))
 call interpolate_gkk (crystal,ifc,elph_ds,kpt_phon)

 !MG: This was the old coding in version 7.6.2:

! call interpolate_gkk (elph_ds,kptirr_phon,kpt_phon,natom,nrpt,phon_ds,rcan,wghatm)
!
! and interpolate_gkk had the prototype:
!
!subroutine interpolate_gkk(elph_ds,kpt_phon,gprim,natom,nrpt,phon_ds,rpt,wghatm)

! hence we were associating kpt_phon to gprim!

end subroutine get_all_gkk2
!!***
