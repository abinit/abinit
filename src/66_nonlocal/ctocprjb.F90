!{\src2tex{textfont=tt}}
!!****f* ABINIT/ctocprjb
!! NAME
!! ctocprjb
!!
!! FUNCTION
!!  Compute <Proj_i,k|Cn,k+b> for set of bands |Cn,k+b> expressed in reciprocal space.
!!  The projector is at  kpoint k, while the bands are at kpoint k+b. This construction
!!  is needed for orbital magnetization, in the computation of the <u_n,k1|H_k2k2|u_m,k3>
!!  matrix elements.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!
!! OUTPUTS
!!  cwaveprj(natom,nspinor) <type(pawcprj_type)>=projected input wave function <Proj_i,k|Cn,k+b> with all NL projectors
!!
!! TODO
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

subroutine ctocprjb(atindx1,cg,cprj,dtorbmag,icgb,idir,ifor,ikgb,&
     & ikptb,mcg,mkmem,mpsang,mpw,natom,nband_k,ncpgr,npwb,nspinor,ntypat,&
     & pawang,pawrad,pawtab,typat,ylm)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_pawang,     only : pawang_type
 use m_pawcprj,    only : pawcprj_type, pawcprj_put, pawcprj_alloc, pawcprj_free, pawcprj_getdim
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_orbmag, only : orbmag_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ctocprjb'
 use interfaces_66_nonlocal, except_this_one => ctocprjb
!End of the abilint section

 implicit none

!Arguments -------------------------------
 !scalars
 integer,intent(in) :: icgb,idir,ifor,ikgb,ikptb,mcg,mkmem,mpsang,mpw
 integer,intent(in) :: natom,nband_k,ncpgr,npwb,nspinor,ntypat
 type(pawang_type),intent(in) :: pawang
 !arrays
 integer,intent(in) :: atindx1(natom),typat(natom)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
 type(orbmag_type),intent(inout) :: dtorbmag
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(pawcprj_type),intent(out) ::  cprj(natom,nband_k)

!Local variables-------------------------------
 !scalars
 integer :: ilm, iband
 !arrays
 integer,allocatable :: dimlmn(:),nattyp_dum(:)
 real(dp),allocatable :: cwavef(:,:),ylm_kb(:,:)
 type(pawcprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

 ABI_ALLOCATE(cwavef,(2,npwb))

 ABI_ALLOCATE(dimlmn,(natom))
 call pawcprj_getdim(dimlmn,natom,nattyp_dum,ntypat,typat,pawtab,'R')
 ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,1))
 call pawcprj_alloc(cwaveprj,ncpgr,dimlmn)

 ABI_ALLOCATE(ylm_kb,(npwb,mpsang*mpsang))
 do ilm=1,mpsang*mpsang
   ylm_kb(1:npwb,ilm)=ylm(1+ikgb:npwb+ikgb,ilm)
 end do
 
 do iband = 1, nband_k

   cwavef(1:2,1:npwb) = cg(1:2,icgb+(iband-1)*npwb+1:icgb+iband*npwb)
   
   call getcprjb(cwavef,cwaveprj,dtorbmag,idir,ifor,ikptb,&
&   mpsang,natom,npwb,nspinor,ntypat,pawang,pawrad,pawtab,typat,ylm_kb)

    ! hard-coded for nspinor 1, nsppol 1
   call pawcprj_put(atindx1,cwaveprj,cprj,natom,iband,0,0,0,1,&
&   nband_k,mkmem,natom,1,nband_k,dimlmn,1,1,0)

 end do

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(ylm_kb)
 ABI_DEALLOCATE(dimlmn)

 call pawcprj_free(cwaveprj)
 ABI_DATATYPE_DEALLOCATE(cwaveprj)

 DBG_EXIT('COLL')

 end subroutine ctocprjb
!!***
