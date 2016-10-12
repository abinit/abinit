!{\src2tex{textfont=tt}}
!!****f* ABINIT/corrmetalwf1
!!
!! NAME
!! corrmetalwf1
!!
!! FUNCTION
!! Response function calculation only:
!! Correct 1st-order wave-function, taking into account "metallic" occupations.
!! 1st-order WF orthogonal to C_n,k+q, restore the "active space" content of the first-order WF.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcgq)=planewave coefficients of wavefunctions at k+q
!!  cprjq(natom,mcprjq)= wave functions at k+q projected with non-local projectors
!!  cwavef(2,npw1*nspinor)= 1st-order wave-function before correction
!!  cwaveprj(natom,nspinor)= 1st-order wave-function before correction
!!                           projected on NL projectors (PAW)
!!  eig1(2*nband**2)=first-order eigenvalues (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  ghc(2,npw1*nspinor)=<G|H0-eig0_k.I|C1 band,k> (NCPP) or <G|H0-eig0_k.S0|C1 band,k> (PAW)
!!                      (C1 before correction)
!!  iband=index of current band
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  istwf_k=option parameter that describes the storage of wfs
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  npw1=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ(nband)=occupation number for each band for each k.
!!  rocceig(nband,nband)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this ratio has been attributed to the band n (second argument), zero otherwise
!!  timcount=index used to accumulate timing (0 from dfpt_vtowfk, 1 from dfpt_nstwf)
!!  usepaw=flag for PAW
!!
!! OUTPUT
!!  cwave1(2,npw1*nspinor)= 1st-order wave-function after correction
!!  cwaveprj1(natom,nspinor)= 1st-order wave-function after correction
!!                            projected on NL projectors (PAW)
!!  edocc(nband)=correction to 2nd-order total energy coming from changes of occupations
!!  wf_corrected=flag put to 1 if input cwave1 is effectively different from output cwavef
!!
!! NOTES
!!  Was part of dfpt_vtowfk before.
!!
!! PARENTS
!!      dfpt_vtowfk
!!
!! CHILDREN
!!      cg_zcopy,dotprod_g,pawcprj_copy,pawcprj_zaxpby,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine corrmetalwf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,edocc,eig1,fermie1,ghc,iband, &
&          ibgq,icgq,istwf_k,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,nspinor,occ,rocceig,timcount,&
&          usepaw,wf_corrected)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_cgtools

 use m_pawcprj, only : pawcprj_type, pawcprj_copy, pawcprj_zaxpby

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'corrmetalwf1'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ibgq,icgq,istwf_k,mcgq,mcprjq,natom,nband,npw1,nspinor,timcount,usepaw
 integer,intent(out) :: wf_corrected
 real(dp),intent(in) :: fermie1
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cgq(2,mcgq),cwavef(2,npw1*nspinor)
 real(dp),intent(in) :: eig1(2*nband**2),ghc(2,npw1*nspinor),occ(nband),rocceig(nband,nband)
 real(dp),intent(out) :: cwave1(2,npw1*nspinor),edocc(nband)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq),cwaveprj(natom,nspinor*usepaw)
 type(pawcprj_type),intent(inout) :: cwaveprj1(natom,nspinor*usepaw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ibandkq,index_cgq,index_cprjq,index_eig1,ii
 real(dp) :: facti,factr,invocc
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwcorr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(214+timcount,1,tsec)

!At this stage, the 1st order function cwavef is orthogonal to cgq (unlike when it is input to dfpt_cgwf).
!Here, restore the "active space" content of the 1st-order wavefunction, to give cwave1 .

!First copy input WF into output WF
 wf_corrected=0
 call cg_zcopy(npw1*nspinor,cwavef,cwave1)

 if (usepaw==1) then
   call pawcprj_copy(cwaveprj,cwaveprj1)
 end if

!Correct WF only for occupied states
 if (abs(occ(iband)) > tol8) then
   invocc=one/occ(iband)

   edocc(iband)=zero

!  Loop over WF at k+q subspace
   do ibandkq=1,nband

!    Select bands with variable occupation
     if (abs(rocceig(ibandkq,iband))>tol8) then

       wf_corrected=1

       index_eig1=2*ibandkq-1+(iband-1)*2*nband
       index_cgq=npw1*nspinor*(ibandkq-1)+icgq

       if(ibandkq==iband) then
         factr=rocceig(ibandkq,iband)*invocc*(eig1(index_eig1)-fermie1)
       else
         factr=rocceig(ibandkq,iband)*invocc*eig1(index_eig1)
       end if
       facti= rocceig(ibandkq,iband)*invocc*eig1(index_eig1+1)

!      Apply correction to 1st-order WF
!$OMP PARALLEL DO PRIVATE(ii) SHARED(cgq,cwave1,facti,factr,index_cgq,npw1,nspinor)
       do ii=1,npw1*nspinor
         cwave1(1,ii)=cwave1(1,ii)+(factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq))
         cwave1(2,ii)=cwave1(2,ii)+(facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq))
       end do

!      In the PAW case, also apply correction to projected WF
       if (usepaw==1) then
         index_cprjq=nspinor*(ibandkq-1)+ibgq
         call pawcprj_zaxpby((/factr,facti/),(/one,zero/),cprjq(:,index_cprjq+1:index_cprjq+nspinor),cwaveprj1)
       end if

!      The factor of two is needed because we compute the 2DTE, and not E(2)
       edocc(iband)=edocc(iband)-two*(factr*eig1(index_eig1)+facti*eig1(index_eig1+1))

     end if ! Variable occupations
   end do ! Loop over k+q subspace
 end if ! occupied states

!In the PAW case, compute <Psi^(1)_ortho|H-Eig0_k.S|Psi^(1)_parallel> contribution to 2DTE
 if (usepaw==1.and.wf_corrected==1) then
   ABI_ALLOCATE(cwcorr,(2,npw1*nspinor))
!$OMP WORKSHARE
   cwcorr(:,:)=cwave1(:,:)-cwavef(:,:)
!$OMP END WORKSHARE
   call dotprod_g(factr,facti,istwf_k,npw1*nspinor,1,cwcorr,ghc,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
   edocc(iband)=edocc(iband)+four*factr
   ABI_DEALLOCATE(cwcorr)
 end if

 call timab(214+timcount,2,tsec)

 DBG_EXIT("COLL")

end subroutine corrmetalwf1
!!***
