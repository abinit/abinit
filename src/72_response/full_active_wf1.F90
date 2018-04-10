!{\src2tex{textfont=tt}}
!!****f* ABINIT/full_active_wf1
!!
!! NAME
!! full_active_wf1
!!
!! FUNCTION
!! Response function calculation only:
!! Restore the full "active space" contribution to the 1st-order wavefunctions.
!! The 1st-order WF corrected in this way will no longer be othogonal to the other occupied states.
!! This routine will be only used in a non self-consistent calculation of the
!! 1st-order WF for post-processing purposes. Therefore, it does not compute
!! the contribution of the 2DTE coming from the change of occupations.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (GKA)
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
!!  eig0nk=energy of the band at k being corrected
!!  eig0_kq(nband)=energies of the bands at k+q
!!  elph2_imagden=imaginary parameter to broaden the energy denominators
!!  iband=index of current band
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icgq=shift to be applied on the location of data in the array cgq
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  npw1=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  timcount=index used to accumulate timing (0 from dfpt_vtowfk, 1 from dfpt_nstwf)
!!  usepaw=flag for PAW
!!
!! OUTPUT
!!  cwave1(2,npw1*nspinor)= 1st-order wave-function after correction
!!  cwaveprj1(natom,nspinor)= 1st-order wave-function after correction
!!                            projected on NL projectors (PAW)
!!
!! NOTES
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

subroutine full_active_wf1(cgq,cprjq,cwavef,cwave1,cwaveprj,cwaveprj1,eig1,&
&               fermie1,eig0nk,eig0_kq,elph2_imagden,&
&               iband,ibgq,icgq,mcgq,mcprjq,mpi_enreg,natom,nband,npw1,&
&               nspinor,timcount,usepaw)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_cgtools
 use m_time, only : timab

 use m_pawcprj, only : pawcprj_type, pawcprj_copy, pawcprj_zaxpby

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'full_active_wf1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,ibgq,icgq,mcgq,mcprjq,natom,nband,npw1,nspinor,timcount,usepaw
 real(dp),intent(in) :: fermie1, eig0nk
 real(dp),intent(in) :: elph2_imagden
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: cgq(2,mcgq),cwavef(2,npw1*nspinor)
 real(dp),intent(in) :: eig0_kq(nband) 
 real(dp),intent(in) :: eig1(2*nband**2)
 real(dp),intent(out) :: cwave1(2,npw1*nspinor)
 type(pawcprj_type),intent(in) :: cprjq(natom,mcprjq),cwaveprj(natom,nspinor*usepaw)
 type(pawcprj_type),intent(inout) :: cwaveprj1(natom,nspinor*usepaw) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ibandkq,index_cgq,index_cprjq,index_eig1,ii
 real(dp) :: facti,factr,invocc,eta,delta_E,inv_delta_E,gkkr
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwcorr(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(214+timcount,1,tsec)

!At this stage, the 1st order function cwavef is orthogonal to cgq (unlike when it is input to dfpt_cgwf).
!Here, restore the "active space" content of the 1st-order wavefunction, to give cwave1 .

!First copy input WF into output WF
 call cg_zcopy(npw1*nspinor,cwavef,cwave1)

 if (usepaw==1) then
   call pawcprj_copy(cwaveprj,cwaveprj1)
 end if

 eta = elph2_imagden

!Loop over WF at k+q subspace
 do ibandkq=1,nband

   delta_E = eig0nk - eig0_kq(ibandkq)
   inv_delta_E = delta_E / ( delta_E ** 2 + eta ** 2)

   index_eig1=2*ibandkq-1+(iband-1)*2*nband
   index_cgq=npw1*nspinor*(ibandkq-1)+icgq

   if(ibandkq==iband) then
     gkkr = eig1(index_eig1) - fermie1
   else
     gkkr = eig1(index_eig1)
   end if
   factr = inv_delta_E * gkkr
   facti = inv_delta_E * eig1(index_eig1+1)

!  Apply correction to 1st-order WF
!$OMP PARALLEL DO PRIVATE(ii) SHARED(cgq,cwave1,facti,factr,index_cgq,npw1,nspinor)
   do ii=1,npw1*nspinor
     cwave1(1,ii)=cwave1(1,ii)+(factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq))
     cwave1(2,ii)=cwave1(2,ii)+(facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq))
   end do

!  In the PAW case, also apply correction to projected WF
   if (usepaw==1) then
     index_cprjq=nspinor*(ibandkq-1)+ibgq
     call pawcprj_zaxpby((/factr,facti/),(/one,zero/),cprjq(:,index_cprjq+1:index_cprjq+nspinor),cwaveprj1)
   end if

 end do ! Loop over k+q subspace

 call timab(214+timcount,2,tsec)

 DBG_EXIT("COLL")

end subroutine full_active_wf1
!!***
