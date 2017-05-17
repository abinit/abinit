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
!! Copyright (C) 2001-2017 ABINIT group (FR, EB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!         if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat1(cplex*nfft,2nspden*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,nspden,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!!  optxc=0 if LDA part of XC kernel has only to be taken into account (even for GGA)
!!       1 if XC kernel has to be fully taken into
!!      -1 if XC kernel does not have to be taken into account
!!  option=if 0, work only with the XC core-correction,
!!         if 1, treat both density change and XC core correction
!!         if 2, treat only density change
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(nfft,nspden)=electron density in real space
!!  rhor1(cplex*nfft,nspden)=array for electron residual density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
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
!!      dfpt_mkvxc,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&
&          nkxc,nkxc_cur,nspden,n3xccc,optnc,option,optxc,paral_kgb,qphon,rhor,rhor1,rprimd,usexcnhat,vxc1,xccc3d1)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxc_noncoll'
 use interfaces_18_timing
 use interfaces_56_xc, except_this_one => dfpt_mkvxc_noncoll
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ixc,n3xccc,nfft,nhat1dim,nhat1grdim,optnc,optxc
 integer,intent(in) :: nkxc,nspden,option,paral_kgb,usexcnhat,nkxc_cur
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: kxc(nfft,nkxc),rhor(cplex*nfft,nspden)
 real(dp),intent(in),target :: rhor1(cplex*nfft,nspden),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)

!Local variables-------------------------------
!scalars
 integer :: ifft, ir
 real(dp),parameter :: m_norm_min=1.d-8
 real(dp) :: dum,dvdn,dvdz,fact,m_dot_m1
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: m_norm(:)
 real(dp),allocatable :: rhor1_diag(:,:)
 real(dp),allocatable :: vxc1_diag(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(181,1,tsec)
 
 if(nspden/=4) then
   MSG_BUG('only for nspden=4!')
 end if

 if(nkxc/=2*min(nspden,2)-1) then
   MSG_BUG('nspden=4 works only with LDA.')
 end if

 if(cplex==2) then
   MSG_BUG('nspden=4 not yet implemented with cplex=2 - that is qphon/=gamma.')
 end if

!Treat first LDA
 if(nkxc==1.or.nkxc==3)then

!FR EB If option=0 (i.e., for XC core-correction only) we apply the correction only on
! the diagonal elements of the potential which are vxc1(:,1:2) since XC core correction
! acts only on the electronic density (i.e., NOT on the magnetization density).
! Then, the corrections on vxc1(:,3:4) are ZERO.

   dvdn=zero; dvdz=zero; dum=zero; fact=zero; vxc1(:,:)=zero
   if (option==0) then
     if (n3xccc==0) then
       vxc1(:,:)=zero
     else !n3xccc/=0
       if(cplex==1)then
         do ir=1,nfft
           vxc1(ir,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(ir)*half
           vxc1(ir,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(ir)*half
         end do
       else
         do ir=1,nfft
           vxc1(2*ir-1,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir-1)*half
           vxc1(2*ir  ,1)=(kxc(ir,1)+kxc(ir,2))*xccc3d1(2*ir  )*half
           vxc1(2*ir-1,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir-1)*half
           vxc1(2*ir  ,2)=(kxc(ir,2)+kxc(ir,3))*xccc3d1(2*ir  )*half
         end do
       end if ! cplex==1
       vxc1(:,3:4)=zero
     end if

   else ! (option=1 or 2)
!  Non-collinear magnetism
!  Has to locally "rotate" rho(r)^(1) (according to magnetization),
!  compute Vxc(r)^(1) and rotate it back
!FR  The collinear routine dfpt_mkvxc wants a general density built as (tr[rho],rho_upup)
!    and the notation has to be consistent with dfpt_accrho and symrhg.

!    PAW: eventually substract compensation density
!    if (usexcnhat==0.and.nhat1dim==1) then
!      ABI_ALLOCATE(rhor1_,(cplex*nfft,nspden))
!      rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
!    else
!      rhor1_ => rhor1
!    end if

     ABI_ALLOCATE(rhor1_diag,(nfft,2))
     ABI_ALLOCATE(vxc1_diag,(nfft,2))
     ABI_ALLOCATE(m_norm,(nfft))

!      -- Rotate rho(r)^(1)
     do ifft=1,nfft
       rhor1_diag(ifft,1)=rhor1(ifft,1) !FR it is already the tr[rhor1] see symrhg.F90
       m_norm(ifft)=sqrt(rhor(ifft,2)**2+rhor(ifft,3)**2+rhor(ifft,4)**2)
       m_dot_m1=rhor(ifft,2)*rhor1(ifft,2)+rhor(ifft,3)*rhor1(ifft,3) &
&       +rhor(ifft,4)*rhor1(ifft,4)
       if (optxc /= -1) then 
         if(m_norm(ifft)>m_norm_min)then
           rhor1_diag(ifft,2)=half*(rhor1_diag(ifft,1)+m_dot_m1/m_norm(ifft)) !rhor1_upup
         else
           rhor1_diag(ifft,2)=half*rhor1_diag(ifft,1)
         end if
       else if (nkxc/=nkxc_cur.and.optxc/=-1) then
         rhor1_diag(ifft,2)=half*(rhor1_diag(ifft,1)+m_dot_m1/m_norm(ifft))
       end if
     end do

!      -- Compute Kxc(r).n^res(r)_rotated
     call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,nhat1,nhat1dim,nhat1gr,nhat1grdim,&
&     nkxc,2,n3xccc,option,paral_kgb,qphon,rhor1_diag,rprimd,usexcnhat,vxc1_diag,xccc3d1)

!    -- Rotate back Vxc(r)^(1)
! FR EB TODO: update the routine to cplex=2
!!   cplex=1:
!!     V is stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     N is stored as : n, m_x, m_y, m_z               (real)
!!   cplex=2:
!!     V is stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     N is stored as : n, m_x, m_y, mZ          (complex)
     if (optnc==1) then
       do ifft=1,nfft
         dvdn=(vxc1_diag(ifft,1)+vxc1_diag(ifft,2))*half
         dvdz=(vxc1_diag(ifft,1)-vxc1_diag(ifft,2))*half
         if(m_norm(ifft)>m_norm_min)then
           fact=dvdz/m_norm(ifft)
           dum=rhor(ifft,4)*fact
           vxc1(ifft,1)=dvdn+dum
           vxc1(ifft,2)=dvdn-dum
           vxc1(ifft,3)= rhor(ifft,2)*fact
           vxc1(ifft,4)=-rhor(ifft,3)*fact
!
         else
           vxc1(ifft,1:2)=dvdn
           vxc1(ifft,3:4)=zero
         end if
       end do
     else
       do ifft=1,nfft
         dvdn=(vxc1_diag(ifft,1)+vxc1_diag(ifft,2))*half
         dvdz=(vxc1_diag(ifft,1)-vxc1_diag(ifft,2))*half
         if(m_norm(ifft)>m_norm_min)then
           dum=dvdz*rhor(ifft,4)/m_norm(ifft)
           vxc1(ifft,1)=dvdn+dum
           vxc1(ifft,2)=dvdn-dum
         else
           vxc1(ifft,1:2)=dvdn
         end if
         vxc1(ifft,3:4)=zero
       end do
     end if ! optnc==1

     ABI_DEALLOCATE(rhor1_diag)
     ABI_DEALLOCATE(vxc1_diag)
     ABI_DEALLOCATE(m_norm)

!    PAW
!   if (option==0.or.(usexcnhat==0.and.nhat1dim==1)) then
!      ABI_DEALLOCATE(rhor1_)
!    end if

   end if ! option==1 or 2
 end if ! nkxc=1 or nkxc=3

 call timab(181,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxc_noncoll
!!***
