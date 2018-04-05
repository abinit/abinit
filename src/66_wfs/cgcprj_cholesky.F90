!{\src2tex{textfont=tt}}
!!****m* ABINIT/cgcprj_cholesky
!! NAME
!!  cgcprj_cholesky
!!
!! FUNCTION
!! Cholesky orthonormalization of the vectors stored in cg+cprj mode.
!!
!! This implementation is NOT band-parallelized
!! Also, it is far of being optimal at the level of linear algebra
!!
!! COPYRIGHT
!! Copyright (C) 2017-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimcprj(natom)=number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom
!!  icg=shift in cg array to locate current k-point and spinpol 
!!  ikpt=current k point index
!!  isppol=current spin polarization index
!!  istwf=input option parameter that describes the storage of wfs
!!  mcg=second dimension of cg array (mpw*nspinor*mband*mkmem*nsppol)
!!  mcprj=second dimension of cprj_k array
!!  mkmem=number of k points which can fit in memory
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nband=number of bands
!!  npw=number of planewaves in basis at this k point
!!  nspinor=number of spinor components
!!  nsppol=number of spin polarizations
!!  ntypat=number of types of atoms
!!  pawtab(dtset%ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  usepaw=1 if PAW is activated
!!
!! SIDE EFFECTS
!!  cg(2,mcg)= plane wave wavefunction coefficients for the set of input wavefunctions (all k points and spinpol)
!!  cprj_k(natom,mcprj) <type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk> with NL projectors for the specific k point and spinpol
!!
!! PARENTS
!!      wf_mixing
!!
!! CHILDREN
!!      dotprod_set_cgcprj,lincom_cgcprj,zpotrf,ztrsm
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine cgcprj_cholesky(atindx1,cg,cprj_k,dimcprj,icg,ikpt,isppol,istwf,mcg,mcprj,mkmem,&
&  mpi_enreg,natom,nattyp,nband,npw,nspinor,nsppol,ntypat,pawtab,usepaw)

 use defs_basis
 use defs_abitypes
 use m_pawtab, only : pawtab_type
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cgcprj_cholesky'
 use interfaces_66_wfs, except_this_one => cgcprj_cholesky
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,ikpt,isppol,istwf,mcg,mcprj,mkmem
 integer,intent(in) :: natom,nband,npw,nspinor,nsppol,ntypat,usepaw
!arrays
 integer, intent(in) :: atindx1(natom),dimcprj(natom),nattyp(ntypat)
 real(dp), intent(inout) :: cg(2,mcg) 
 type(pawcprj_type),intent(inout) :: cprj_k(natom,mcprj)
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer :: hermitian,ierr,ii,inplace
!arrays
 real(dp), allocatable :: dmn(:,:,:),smn(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(smn,(2,nband,nband))
 ABI_ALLOCATE(dmn,(2,nband,nband))

 hermitian=1
 call dotprod_set_cgcprj(atindx1,cg,cg,cprj_k,cprj_k,dimcprj,hermitian,&
& 0,0,icg,icg,ikpt,isppol,istwf,nband,mcg,mcg,mcprj,mcprj,mkmem,&
& mpi_enreg,natom,nattyp,nband,nband,npw,nspinor,nsppol,ntypat,pawtab,smn,usepaw)

!Cholesky factorization: O = U^H U with U upper triangle matrix.
 call ZPOTRF('U',nband,smn,nband,ierr)

!Solve X U = 1. 
 dmn=zero
 do ii=1,nband
   dmn(1,ii,ii)=one 
 end do
 call ZTRSM('Right','Upper','Normal','Normal',nband,nband,cone,smn,nband,dmn,nband)

 inplace=1
!This call does not take into account the fact that X=dmn is an upper triangular matrix...
!The number of operations might be divided by two.
 call lincom_cgcprj(dmn,cg,cprj_k,dimcprj,&
& icg,inplace,mcg,mcprj,natom,nband,nband,npw,nspinor,usepaw)

 ABI_DEALLOCATE(smn)
 ABI_DEALLOCATE(dmn)

 end subroutine cgcprj_cholesky
!!***
