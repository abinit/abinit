!!****f* defs_wvltypes/wvl_descr_atoms_set_sym
!!
!! NAME
!! wvl_descr_atoms_set_sym
!!
!! FUNCTION
!! Add symmetry information to  wvl%atoms data structure.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=unit cell length scales (bohr)
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>= wavelet type
!!                 | nat      =  number of atoms
!!                 | ntypes   =  number of species
!!                 | alat1    =  acell(1)
!!                 | alat2    =  acell(2)
!!                 | alat3    =  acell(3)
!!                 | iatype   =  types for atoms
!!                 | lfrztyp  =  flag for the movement of atoms.
!!                 | natpol   =  integer related to polarisation at the first step
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_atoms_set_sym(wvl, efield, irrzon, nsppol, nsym, phnons, &
           & symafm, symrel, tnons, tolsym)

 use m_profiling_abi
 use m_errors

 use defs_basis
 use defs_wvltypes
 use m_ab7_symmetry
#if defined HAVE_BIGDFT
 use BigDFT_API, only: astruct_set_symmetries
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_atoms_set_sym'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: nsppol,nsym
  real(dp), intent(in)                   :: tolsym
  type(wvl_internal_type), intent(inout) :: wvl
  !arrays
  integer, intent(in)                    :: symafm(3,3,nsym), symrel(3,3,nsym)
  integer, target, intent(in)            :: irrzon(:,:,:)
  real(dp), target, intent(in)           :: phnons(:,:,:)
  real(dp), intent(in)                   :: efield(3), tnons(3,nsym)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
  integer :: errno
 character(len=500) :: message
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 write(message, '(a,a)' ) ch10,&
& ' wvl_descr_atoms_set_sym: Create symmetries for the wvl%atoms object.'
 call wrtout(std_out,message,'COLL')

 wvl%atoms%astruct%sym%symObj = -1
 nullify(wvl%atoms%astruct%sym%irrzon)
 nullify(wvl%atoms%astruct%sym%phnons)
 call astruct_set_symmetries(wvl%atoms%astruct, (nsym <= 1), tolsym, efield, nsppol)
 if (nsym > 1) then
   call symmetry_set_n_sym(wvl%atoms%astruct%sym%symObj, nsym, symrel, tnons, symafm, errno)
 end if
 wvl%atoms%astruct%sym%irrzon => irrzon
 wvl%atoms%astruct%sym%phnons => phnons

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) nsppol,nsym,tolsym,wvl%h(1),symafm(1,1,1),symrel(1,1,1),&
& irrzon(1,1,1),phnons(1,1,1),efield(1),tnons(1,1)
#endif  

end subroutine wvl_descr_atoms_set_sym
!!***
