!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setBoxGeometry
!! NAME
!! wvl_setBoxGeometry
!!
!! FUNCTION
!! When wavelets are used, the box definition needs to be changed.
!! The box size is recomputed knowing some psp informations such as
!! the radius for coarse and fine grid. Then, the atoms are translated
!! to be included in the new box. Finally the FFT grid is computed using
!! the fine wavelet mesh and a buffer characteristic of used wavelets plus
!! a buffer used to be multiple of 2, 3 or 5.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  radii= the radii for each type of atoms, giving the fine and the coarse
!!         grid.
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! SIDE EFFECTS
!!  wvl <type(wvl_internal_type)>=internal variables used by wavelets, describing
!!                             the box are set.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      gstate,wvl_memory,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      mkrdim,nullify_locreg_descriptors,system_size,wrtout,xcart2xred
!!      xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_setBoxGeometry(prtvol, radii, rprimd, xred, wvl, wvl_crmult, wvl_frmult)

 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors

#if defined HAVE_BIGDFT
 use BigDFT_API, only: system_size,nullify_locreg_descriptors
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_setBoxGeometry'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtvol
 real(dp), intent(in) :: wvl_crmult, wvl_frmult
 type(wvl_internal_type), intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: radii(:,:)
 real(dp),intent(inout) :: rprimd(3,3),xred(:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ii
 logical,parameter :: OCLconv=.false.
 character(len=500) :: message
!arrays
 real(dp) :: rprim(3,3),acell(3)
 real(dp),allocatable :: xcart(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT
 if (prtvol == 0) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' wvl_setBoxGeometry : Changing the box for wavelets computation.'
   call wrtout(std_out,message,'COLL')
 end if

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, wvl%atoms%astruct%nat))
 call xred2xcart(wvl%atoms%astruct%nat, rprimd, xcart, xred)

 call nullify_locreg_descriptors(wvl%Glr)
 call system_size(wvl%atoms, xcart, radii, wvl_crmult, &
& wvl_frmult, wvl%h(1), wvl%h(2), wvl%h(3), OCLconv, wvl%Glr, wvl%shift)

 acell(:) = wvl%atoms%astruct%cell_dim(:) 

 if (prtvol == 0) then
   write(message, '(a,3F12.6)' ) &
&   '  | acell is now:         ', acell
   call wrtout(std_out,message,'COLL')
   write(message, '(a,2I5,a,a,2I5,a,a,2I5)' ) &
&   '  | nfl1, nfu1:           ', wvl%Glr%d%nfl1, wvl%Glr%d%nfu1, ch10, &
&   '  | nfl2, nfu2:           ', wvl%Glr%d%nfl2, wvl%Glr%d%nfu2, ch10, &
&   '  | nfl3, nfu3:           ', wvl%Glr%d%nfl3, wvl%Glr%d%nfu3
   call wrtout(std_out,message,'COLL')
 end if

!Change the metric to orthogonal one
 rprim(:, :) = real(0., dp)
 do ii = 1, 3, 1
   rprim(ii,ii) = real(1., dp)
 end do
 call mkrdim(acell, rprim, rprimd)

!Save shifted atom positions into xred
 call xcart2xred(wvl%atoms%astruct%nat, rprimd, xcart, xred)
 ABI_DEALLOCATE(xcart)

 if (prtvol == 0) then
   write(message, '(a,3I12)' ) &
&   '  | box size for datas:   ', wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i
   call wrtout(std_out,message,'COLL')
   write(message, '(a,3I12)' ) &
&   '  | box size for wavelets:', wvl%Glr%d%n1, wvl%Glr%d%n2, wvl%Glr%d%n3
   call wrtout(std_out,message,'COLL')
 end if
 
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  prtvol,wvl_crmult,wvl_frmult,wvl%h(1),&
& radii(1,1),rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_setBoxGeometry
!!***
