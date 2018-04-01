!{\src2tex{textfont=tt}}
!!****f* ABINIT/setsym
!! NAME
!! setsym
!!
!! FUNCTION
!! Set up irreducible zone in  G space by direct calculation.
!! Do not call this routine if nsym=1 (only identity symmetry).
!! Only indsym and symrec get returned if iscf=0.
!! symrec needed to symmetrize coordinate gradients in sygrad.
!! (symrec is redundant and could be removed later in favor of symrel)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! iscf=(<= 0  =>non-SCF), >0 => SCF
!! natom=number of atoms in unit cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetries in space group (at least 1)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in terms of real space
!! primitive translations
!! tnons(3,nsym)=nonsymmorphic translations of space group in terms
!! of real space primitive translations (may be 0)
!! typat(natom)=atom type (integer) for each atom
!! xred(3,natom)=atomic coordinates in terms of real space primitive translations
!!
!! OUTPUT
!! indsym(4,nsym,natom)=indirect indexing of atom labels--see subroutine
!!   symatm for definition (if nsym>1)
!! irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!! phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!! symrec(3,3,nsym)=symmetry operations in terms of reciprocal
!!   space primitive translations (if nsym>1)
!!
!! NOTES
!! nsppol and nspden are needed in case of (anti)ferromagnetic symmetry operations
!!
!! PARENTS
!!      dfpt_looppert,gstate,m_dvdb,nonlinear,respfn,scfcv
!!
!! CHILDREN
!!      chkgrp,irrzg,mati3inv,symatm,symdet,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setsym(indsym,irrzon,iscf,natom,nfft,ngfft,nspden,nsppol,nsym,phnons,&
& symafm,symrec,symrel,tnons,typat,xred)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_time,          only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setsym'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace, except_this_one => setsym
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,natom,nfft,nspden,nsppol,nsym
!arrays
 integer,intent(in) :: ngfft(18),symafm(nsym),symrel(3,3,nsym),typat(natom)
 integer,intent(out) :: indsym(4,nsym,natom)
 integer,intent(inout) :: irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4)) !vz_i
 integer,intent(out) :: symrec(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)
 real(dp),intent(out) :: phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))

!Local variables-------------------------------
!scalars
 integer :: isym,ierr
 real(dp) :: tolsym8
!arrays
 integer,allocatable :: determinant(:)
 real(dp) :: tsec(2)

! *************************************************************************

 call timab(6,1,tsec)

!Check that symmetries have unity determinant
 ABI_ALLOCATE(determinant,(nsym))
 call symdet(determinant,nsym,symrel)
 ABI_DEALLOCATE(determinant)


!Get the symmetry matrices in terms of reciprocal basis
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!Check for group closure
 call chkgrp(nsym,symafm,symrel,ierr)
 ABI_CHECK(ierr==0,"Error in group closure")

 call chkgrp(nsym,symafm,symrec,ierr)
 ABI_CHECK(ierr==0,"Error in group closure")

!Obtain a list of rotated atom labels:
 tolsym8=tol8
 call symatm(indsym,natom,nsym,symrec,tnons,tolsym8,typat,xred)

!If non-SCF calculation, or nsym==1, do not need IBZ data
 if ( (iscf>0 .or. iscf==-3) .and. nsym>1 ) then
!  Locate irreducible zone in reciprocal space for symmetrization:
   call irrzg(irrzon,nspden,nsppol,nsym,ngfft(1),ngfft(2),ngfft(3),phnons,symafm,symrel,tnons)
 end if

 call timab(6,2,tsec)

end subroutine setsym
!!***
