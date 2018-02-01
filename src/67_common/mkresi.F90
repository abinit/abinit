!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkresi
!! NAME
!! mkresi
!!
!! FUNCTION
!! Make residuals from knowledge of wf in G space and application of Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cg(2,mcg)=<G|Cnk>=Fourier coefficients of wavefunction
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  ikpt=index of k-point
!!  isppol=index of spin
!!  mcg=second dimension of the cg array
!!  mpi_enreg=information about MPI parallelization
!!  nband=number of bands involved in subspace matrix.
!!  npw=number of planewaves in basis sphere at this k point.
!!  prtvol=control print volume and debugging output
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  eig_k(nband)$= \langle C_n \mid H \mid C_n \rangle $ for each band.
!!  resid_k(nband)=residual for each band
!!   $= \langle C_n \mid H H \mid C_n \rangle- \langle C_n \mid H \mid C_n \rangle^2 $.
!!
!! PARENTS
!!      energy
!!
!! CHILDREN
!!      dotprod_g,getghc,prep_getghc,sqnorm_g,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkresi(cg,eig_k,gs_hamk,icg,ikpt,isppol,mcg,mpi_enreg,nband,prtvol,resid_k)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_cgtools

 use m_pawcprj,     only : pawcprj_type
 use m_hamiltonian, only : gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkresi'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,ikpt,isppol,mcg,nband,prtvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(inout) :: gs_hamk
!arrays
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: eig_k(nband),resid_k(nband)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_getghc=3
 integer :: blocksize,cpopt,iband,iband_last,iblock,iblocksize,ipw,ipw_shift
 integer :: my_nspinor,nblockbd,npw_k
 real(dp) :: doti,dotr
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable,target :: cwavef(:,:),ghc(:,:),gsc(:,:),gvnlc(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: cwavef_ptr(:,:),ghc_ptr(:,:),gsc_ptr(:,:)
 type(pawcprj_type) :: cwaveprj(0,0)

! *************************************************************************

!Keep track of total time spent in mkresi
 call timab(13,1,tsec)

!Parallelism over FFT and/or bands: define sizes and tabs
 my_nspinor=max(1,gs_hamk%nspinor/mpi_enreg%nproc_spinor)
 if (mpi_enreg%paral_kgb==1) then
   nblockbd=nband/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
 else
   nblockbd=nband/mpi_enreg%nproc_fft
   if (nband/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
 end if
 blocksize=nband/nblockbd

 npw_k=gs_hamk%npw_k
 ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(ghc,(2,npw_k*my_nspinor))
 ABI_ALLOCATE(gvnlc,(2,npw_k*my_nspinor))
 if (gs_hamk%usepaw==1)  then
   ABI_ALLOCATE(gsc,(2,npw_k*my_nspinor))
 else
   ABI_ALLOCATE(gsc,(0,0))
 end if

!Loop over (blocks of) bands
 do iblock=1,nblockbd
   iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband)
   if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,mpi_enreg%me_kpt)) cycle

!  Load |Cn>
   ipw_shift=(iblock-1)*npw_k*my_nspinor*blocksize+icg
!$OMP PARALLEL DO
   do ipw=1,npw_k*my_nspinor*blocksize
     cwavef(1,ipw)=cg(1,ipw+ipw_shift)
     cwavef(2,ipw)=cg(2,ipw+ipw_shift)
   end do

!  Compute H|Cn>
   cpopt=-1
   if (mpi_enreg%paral_kgb==0) then
     call getghc(cpopt,cwavef,cwaveprj,ghc,gsc,gs_hamk,gvnlc,zero,mpi_enreg,1,&
&     prtvol,gs_hamk%usepaw,tim_getghc,0)
   else
     call prep_getghc(cwavef,gs_hamk,gvnlc,ghc,gsc,zero,nband,mpi_enreg,&
&     prtvol,gs_hamk%usepaw,cpopt,cwaveprj,&
&     already_transposed=.false.)
   end if

!  Compute the residual, <Cn|(H-<Cn|H|Cn>)**2|Cn>:
   do iblocksize=1,blocksize
     iband=(iblock-1)*blocksize+iblocksize
     ipw_shift=(iblocksize-1)*npw_k*my_nspinor
     cwavef_ptr => cwavef(:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)
     ghc_ptr    => ghc   (:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)

!    First get eigenvalue <Cn|H|Cn>:
     call dotprod_g(dotr,doti,gs_hamk%istwf_k,npw_k*my_nspinor,1,cwavef_ptr,ghc_ptr,&
&     mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
     eig_k(iband)=dotr

!    Next need <G|(H-S<Cn|H|Cn>)|Cn> (in ghc):
     if (gs_hamk%usepaw==0) then
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(cwavef_ptr,ghc_ptr,eig_k,iband,npw_k,my_nspinor)
       do ipw=1,npw_k*my_nspinor
         ghc_ptr(1,ipw)=ghc_ptr(1,ipw)-eig_k(iband)*cwavef_ptr(1,ipw)
         ghc_ptr(2,ipw)=ghc_ptr(2,ipw)-eig_k(iband)*cwavef_ptr(2,ipw)
       end do
     else
       gsc_ptr => gsc(:,1+ipw_shift:npw_k*my_nspinor+ipw_shift)
!$OMP PARALLEL DO PRIVATE(ipw) SHARED(gsc_ptr,ghc_ptr,eig_k,iband,npw_k,my_nspinor)
       do ipw=1,npw_k*my_nspinor
         ghc_ptr(1,ipw)=ghc_ptr(1,ipw)-eig_k(iband)*gsc_ptr(1,ipw)
         ghc_ptr(2,ipw)=ghc_ptr(2,ipw)-eig_k(iband)*gsc_ptr(2,ipw)
       end do
     end if

!    Then simply square the result:
     call sqnorm_g(dotr,gs_hamk%istwf_k,npw_k*my_nspinor,ghc_ptr,&
&     mpi_enreg%me_g0,mpi_enreg%comm_fft)
     resid_k(iband)=dotr

   end do ! iblocksize

 end do ! iblock

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlc)
 ABI_DEALLOCATE(gsc)

 call timab(13,2,tsec)

end subroutine mkresi
!!***
