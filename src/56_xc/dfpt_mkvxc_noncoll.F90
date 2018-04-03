!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkvxc_noncoll
!! NAME
!! dfpt_mkvxc_noncoll
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to atomic displacement for non-collinear spins: assemble the first-order
!! density change with the frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (FR, EB, SPr)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!         if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhotoxc.F90)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- GS compensation density
!!  nhatdim= -PAW only- 1 if nhat array is used ; 0 otherwise
!!  nhat1(cplex*nfft,nspden*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!!  option=if 0, work only with the XC core-correction,
!!         if 1, treat both density change and XC core correction
!!         if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(nfft,nspden)=GS electron density in real space
!!  rhor1(cplex*nfft,nspden)=1st-order electron density in real space
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  vxc(nfft,nspden)=GS XC potential
!!
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_nstdy,dfpt_nstpaw,dfpt_rhotov,nres2vres
!!
!! CHILDREN
!!      dfpt_mkvxc,rotate_back_mag,rotate_back_mag_dfpt,rotate_mag,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat,nhatdim,nhat1,nhat1dim,&
&          nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,optnc,option,paral_kgb,qphon,rhor,rhor1,&
&          rprimd,usexcnhat,vxc,vxc1,xccc3d1,ixcrot)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xc_noncoll

 use m_time,      only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxc_noncoll'
 use interfaces_56_xc, except_this_one => dfpt_mkvxc_noncoll
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ixc,n3xccc,nfft,nhatdim,nhat1dim,nhat1grdim,optnc
 integer,intent(in) :: nkxc,nspden,option,paral_kgb,usexcnhat
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: vxc(nfft,nspden)
 real(dp),intent(in) :: nhat(nfft,nspden*nhatdim),nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in),target :: rhor(nfft,nspden),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: qphon(3),rprimd(3,3),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 integer,optional,intent(in) :: ixcrot
!Local variables-------------------------------
!scalars
!arrays
 real(dp) :: nhat1_zero(0,0),nhat1gr_zero(0,0,0),tsec(2)
 real(dp),allocatable :: m_norm(:),rhor1_diag(:,:),vxc1_diag(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: mag(:,:),rhor_(:,:),rhor1_(:,:)
! *************************************************************************

!  Non-collinear magnetism
!  Has to locally "rotate" rho(r)^(1) (according to magnetization),
!  Compute Vxc(r)^(1) in the spin frame aligned with \vec{m} and rotate it back

 DBG_ENTER("COLL")
 ABI_UNUSED(nhat1gr)

 call timab(181,1,tsec)

 if(nspden/=4) then
   MSG_BUG('only for nspden=4!')
 end if

 if(nkxc/=2*min(nspden,2)-1) then
   MSG_BUG('nspden=4 works only with LSDA.')
 end if

!Special case: no XC applied
 if (ixc==0.or.nkxc==0) then
   MSG_WARNING('Note that no xc is applied (ixc=0)')
   vxc1(:,:)=zero
   return
 end if



!Treat first LDA
 if(nkxc==1.or.nkxc==3)then

   vxc1(:,:)=zero

!  PAW: possibly substract compensation density
   if (usexcnhat==0.and.nhatdim==1) then
     ABI_ALLOCATE(rhor_,(nfft,nspden))
     rhor_(:,:) =rhor(:,:)-nhat(:,:)
   else
     rhor_ => rhor
   end if
   if (usexcnhat==0.and.nhat1dim==1) then
     ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
     rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
   else
     rhor1_ => rhor1
   end if

!  Magnetization
   mag => rhor_(:,2:4)
   ABI_ALLOCATE(rhor1_diag,(cplex*nfft,2))
   ABI_ALLOCATE(vxc1_diag,(cplex*nfft,2))
   ABI_ALLOCATE(m_norm,(nfft))

!  -- Rotate rho(r)^(1)
!  SPr: for option=0 the rhor is not used, only core density xccc3d1
!       rotate_mag is only to compute the m_norm
   call rotate_mag(rhor1_,rhor1_diag,mag,nfft,cplex,mag_norm_out=m_norm,&
&   rho_out_format=2)

!  -- Compute Vxc(r)^(1)=Kxc(r).rho(r)^(1)_rotated
!  Note for PAW: nhat has already been substracted; don't use it in dfpt_mkvxc
!                 (put all nhat options to zero).
!  The collinear routine dfpt_mkvxc wants a general density built as (tr[rho],rho_upup)
   call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1_zero,0,nhat1gr_zero,0,&
&   nkxc,2,n3xccc,option,paral_kgb,qphon,rhor1_diag,rprimd,0,vxc1_diag,xccc3d1)

   !call test_rotations(0,1)

!  -- Rotate back Vxc(r)^(1)
   if (optnc==1) then
     if(present(ixcrot)) then
       call rotate_back_mag_dfpt(option,vxc1_diag,vxc1,vxc,kxc,rhor1_,mag,nfft,cplex,&
&       mag_norm_in=m_norm,rot_method=ixcrot)
     else
       call rotate_back_mag_dfpt(option,vxc1_diag,vxc1,vxc,kxc,rhor1_,mag,nfft,cplex,&
&       mag_norm_in=m_norm)
     end if
   else
     call rotate_back_mag(vxc1_diag,vxc1,mag,nfft,mag_norm_in=m_norm)
     vxc1(:,3:4)=zero
   end if

   ABI_DEALLOCATE(rhor1_diag)
   ABI_DEALLOCATE(vxc1_diag)
   ABI_DEALLOCATE(m_norm)
   if (usexcnhat==0.and.nhatdim==1) then
     ABI_DEALLOCATE(rhor_)
   end if
   if (usexcnhat==0.and.nhat1dim==1) then
     ABI_DEALLOCATE(rhor1_)
   end if

 end if ! nkxc=1 or nkxc=3

 call timab(181,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxc_noncoll
!!***
