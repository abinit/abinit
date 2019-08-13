!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lobpcg
!! NAME
!!  m_lobpcg
!!
!! FUNCTION
!!  This module provides the procedures used in the LOBPCGWF routine.
!!  They permit to hide the complex/real form of the WFs.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (FBottin,CS,FDahm,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_lobpcg

 use defs_basis
 use m_abicore
 use m_errors
 use defs_abitypes
 use m_wfutils
 use m_abi_linalg
 use m_cgtools
 use m_dtset

 use m_time,              only : timab

 implicit none

 private

!public procedures.
 public :: xprecon
!!***

CONTAINS
!----------------------------------------------------------------------


!! NAME
!!  xprecon
!!
!! FUNCTION
!!  precondition $<G|(H-e_{n,k})|C_{n,k}>$
!!  for a block of band (band-FFT parallelisation)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (FBottin,CS)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of blocks of bands
!!  $cg(vectsize,blocksize)=<G|C_{n,k}> for a block of bands$.
!!  $eval(blocksize,blocksize)=current block of bands eigenvalues=<C_{n,k}|H|C_{n,k}>$.
!!  $ghc(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  iterationnumber=number of iterative minimizations in LOBPCG
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  mpi_enreg=informations about MPI parallelization
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  $vect(vectsize,blocksize)=<G|H|C_{n,k}> for a block of bands$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996) [[cite:Kresse1996]]
!!           0 otherwise
!!  optpcon= 0 the TPA preconditionning matrix does not depend on band
!!           1 the TPA preconditionning matrix (not modified)
!!           2 the TPA preconditionning matrix is independant of iterationnumber
!!  vectsize= size of vectors
!!
!! OUTPUT
!!  vect(2,npw)=<g|(h-eval)|c_{n,k}>*(polynomial ratio)
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!      abi_xcopy,cg_precon_block,cg_zprecon_block,timab
!!
!! SOURCE

subroutine xprecon(cg,eval,blocksize,iterationnumber,kinpw,&
&  mpi_enreg,npw,nspinor,optekin,optpcon,pcon,ghc,vect,vectsize,&
&  timopt,tim_xprecon) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,iterationnumber,npw,nspinor,optekin
 integer,intent(in) :: optpcon,vectsize
 integer, intent(in), optional :: timopt,tim_xprecon
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: cg(vectsize,blocksize),eval(blocksize,blocksize)
 real(dp),intent(in) :: kinpw(npw)
 real(dp),intent(inout) :: ghc(vectsize,blocksize)
 real(dp),intent(inout) :: pcon(npw,blocksize),vect(vectsize,blocksize)

!Local variables-------------------------------
 complex(dpc),dimension(:,:),allocatable :: z_cg,z_eval,z_ghc,z_vect
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_xprecon).and.present(timopt)) then
   if(abs(timopt)==3) then
    call timab(tim_xprecon,1,tsec)
   end if
 end if

 if ( x_cplx == 1 ) then
   call cg_precon_block(cg,eval,blocksize,iterationnumber,kinpw,&
&    npw,nspinor,mpi_enreg%me_g0,optekin,optpcon,pcon,ghc,vect,vectsize,mpi_enreg%comm_bandspinorfft)
 else
   ABI_ALLOCATE(z_cg,(vectsize,blocksize))
   ABI_ALLOCATE(z_eval,(blocksize,blocksize))
   ABI_ALLOCATE(z_ghc,(vectsize,blocksize))
   ABI_ALLOCATE(z_vect,(vectsize,blocksize))

   call abi_xcopy(x_cplx*vectsize*blocksize,cg,1,z_cg,1)
   call abi_xcopy(x_cplx*vectsize*blocksize,ghc,1,z_ghc,1)
   call abi_xcopy(x_cplx*vectsize*blocksize,vect,1,z_vect,1)
   call abi_xcopy(x_cplx*blocksize*blocksize,eval,1,z_eval,1)

   call cg_zprecon_block(z_cg,z_eval,blocksize,iterationnumber,kinpw,&
&    npw,nspinor,optekin,optpcon,pcon,z_ghc,z_vect,vectsize,mpi_enreg%comm_bandspinorfft)

   call abi_xcopy(x_cplx*vectsize*blocksize,z_cg,1,cg,1)
   call abi_xcopy(x_cplx*vectsize*blocksize,z_ghc,1,ghc,1)
   call abi_xcopy(x_cplx*vectsize*blocksize,z_vect,1,vect,1)
   call abi_xcopy(x_cplx*blocksize*blocksize,z_eval,1,eval,1)

   ABI_DEALLOCATE(z_cg)
   ABI_DEALLOCATE(z_eval)
   ABI_DEALLOCATE(z_ghc)
   ABI_DEALLOCATE(z_vect)
 endif

 if (present(tim_xprecon).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xprecon,2,tsec)
   end if
 end if

end subroutine xprecon
!!***

!----------------------------------------------------------------------

END MODULE m_lobpcg
!!***
