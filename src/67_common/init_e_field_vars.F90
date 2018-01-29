!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_e_field_vars
!! NAME
!! init_e_field_vars
!!
!! FUNCTION
!! Initialization of variables and data structures used in polarization
!! calculations
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group 
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3) = primitive translations in recip space
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  mpi_enreg=information about MPI parallelization
!!  npwarr(nkpt) = number of planewaves in basis and boundary at this k point
!!  occ(mband*nkpt*nsppol) = occup number for each band at each k point
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive vectors
!!  symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!    reciprocal space primitive translations
!!  xred(3,natom) = location of atoms in reduced units
!!
!! OUTPUT
!!  dtefield <type(efield_type)> :: initialized polarization variables
!!  pwind(pwind_alloc,2,3) = array used to compute the overlap matrix smat
!!                         between k-points k and k +- dk where dk is
!!                         parallel to the direction idir
!!  pwind_alloc = first dimension of pwind and pwnsfac
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!
!! SIDE EFFECTS
!!
!! TO DO
!!
!! NOTES
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      initberry
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine init_e_field_vars(dtefield,dtset,gmet,gprimd,kg,&
&              mpi_enreg,npwarr,occ,pawang,pawrad,pawtab,psps,&
&              pwind,pwind_alloc,pwnsfac,rprimd,symrec,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_efield

 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init_e_field_vars'
 use interfaces_67_common, except_this_one => init_e_field_vars
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: pwind_alloc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield !vz_i needs efield2
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 integer,intent(in) :: symrec(3,3,dtset%nsym)
 integer,pointer :: pwind(:,:,:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,dtset%natom)
 real(dp),pointer :: pwnsfac(:,:)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
 logical :: initfield
!scalars

! *************************************************************************

 initfield = .false.

!initialization
 dtefield%has_qijb = 0
 dtefield%has_epawf3 = 0
 dtefield%has_epaws3 = 0
 dtefield%has_expibi = 0
 dtefield%has_rij = 0
 dtefield%usecprj = 0
 dtefield%berryopt = 0

 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4) .or. (dtset%berryopt == 6) .or.(dtset%berryopt == 7) .or. &
& (dtset%berryopt == 14) .or.(dtset%berryopt == 16) .or.(dtset%berryopt == 17)) then 
   nullify(pwind,pwnsfac)
   call initberry(dtefield,dtset,gmet,gprimd,kg,&
&   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
&   dtset%natom,dtset%nkpt,npwarr,dtset%nsppol,&
&   dtset%nsym,dtset%ntypat,occ,pawang,pawrad,pawtab,&
&   psps,pwind,pwind_alloc,pwnsfac,rprimd,symrec,&
&   dtset%typat,psps%usepaw,xred)
   initfield = .true.
 end if

 if (.not. initfield) then
   pwind_alloc = 1
   ABI_ALLOCATE(pwind,(pwind_alloc,2,3))
   ABI_ALLOCATE(pwnsfac,(2,pwind_alloc))
   pwind(:,:,:)=0
   pwnsfac(:,:)=zero
 end if

end subroutine init_e_field_vars
!!***
