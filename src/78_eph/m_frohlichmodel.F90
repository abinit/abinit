!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_frohlichmodel
!! NAME
!!  m_frohlichmodel
!!
!! FUNCTION
!!  Compute ZPR, temperature-dependent electronic structure, and other properties
!!  using the Frohlich model 
!!
!! COPYRIGHT
!!  Copyright (C) 2018 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_frohlichmodel

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 private

 public :: frohlichmodel

contains
!!***

!!****f* m_frohlichmodel/frohlichmodel
!! NAME
!!  frohlichmodel
!!
!! FUNCTION
!! Main routine to compute properties based on the Frohlich model
!!
!! INPUTS
!! cryst<crystal_t>=Structure defining the unit cell (geometry, atomic positions and symmetry operations in real and reciprocal space)
!! dtfil<datafiles_type>=Variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=Electronic band structure information
!! efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band degeneracy at each k point
!! efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!   efmasval(:,:)%eig2_diag band curvature double tensor
!! ifc<ifc_type>=contains the dynamical matrix and the IFCs.
!!
!! PARENTS
!! eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine frohlichmodel(cryst,dtfil,dtset,ebands,efmasdeg,efmasval,ifc)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
!use m_xmpi
!use m_xomp
 use m_errors
!use m_hdr
 use m_crystal
 use m_crystal_io
 use m_ebands
 use m_efmas_defs
 use m_ifc

 use m_gaussian_quadrature, only : cgqf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'frohlichmodel'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
!arrays
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables ------------------------------
!scalars
 logical :: sign_warn
 integer :: deg_dim,iband,ideg,ikpt,info,iphi,itheta
 integer :: jband,lwork,nphi,ntheta
 real(dp) :: angle_phi,cosph,costh,sinph,sinth,weight,weight_phi
 character(len=500) :: msg
!arrays
 logical, allocatable :: saddle_warn(:), start_eigf3d_pos(:)
 real(dp) :: kpt(3),unit_r(3)
 real(dp), allocatable :: eigenval(:), rwork(:)
 real(dp), allocatable :: m_avg(:), m_avg_frohlich(:)
 real(dp), allocatable :: gq_points_th(:),gq_points_costh(:),gq_points_sinth(:),gq_weights_th(:)
 real(dp), allocatable :: gq_points_cosph(:),gq_points_sinph(:)
 complex(dpc), allocatable :: eigenvec(:,:), work(:)
 complex(dpc), allocatable :: eig2_diag_cart(:,:,:,:)
 complex(dpc), allocatable :: f3d(:,:)


!************************************************************************

 !!! Initialization of integrals 
 ntheta   = dtset%efmas_ntheta
 nphi     = 2*ntheta
 ABI_ALLOCATE(gq_points_th,(ntheta))
 ABI_ALLOCATE(gq_points_costh,(ntheta))
 ABI_ALLOCATE(gq_points_sinth,(ntheta))
 ABI_ALLOCATE(gq_weights_th,(ntheta))
 ABI_ALLOCATE(gq_points_cosph,(nphi))
 ABI_ALLOCATE(gq_points_sinph,(nphi))
 call cgqf(ntheta,1,zero,zero,zero,pi,gq_points_th,gq_weights_th)
 do itheta=1,ntheta
   gq_points_costh(itheta)=cos(gq_points_th(itheta))
   gq_points_sinth(itheta)=sin(gq_points_th(itheta))
 enddo
 weight_phi=two*pi/real(nphi,dp)
 do iphi=1,nphi
   angle_phi=weight_phi*(iphi-1)
   gq_points_cosph(iphi)=cos(angle_phi)
   gq_points_sinph(iphi)=sin(angle_phi)
 enddo

 do ikpt=1,dtset%nkpt

   kpt(:)=dtset%kptns(:,ikpt)
   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)

     deg_dim    = efmasdeg(ikpt)%degs_bounds(2,ideg) - efmasdeg(ikpt)%degs_bounds(1,ideg) + 1

     ABI_ALLOCATE(eig2_diag_cart,(3,3,deg_dim,deg_dim))

     !Convert eig2_diag to cartesian coordinates
     do iband=1,deg_dim
        do jband=1,deg_dim
          eig2_diag_cart(:,:,iband,jband)=efmasval(ideg,ikpt)%eig2_diag(:,:,iband,jband)
          eig2_diag_cart(:,:,iband,jband) = matmul(matmul(cryst%rprimd,eig2_diag_cart(:,:,iband,jband)),transpose(cryst%rprimd))/two_pi**2
        enddo
     enddo

     ABI_ALLOCATE(f3d,(deg_dim,deg_dim))
     ABI_ALLOCATE(m_avg,(deg_dim))
     ABI_ALLOCATE(m_avg_frohlich,(deg_dim))
     ABI_ALLOCATE(eigenval,(deg_dim))
     ABI_ALLOCATE(saddle_warn,(deg_dim))
     ABI_ALLOCATE(start_eigf3d_pos,(deg_dim))

     m_avg=zero
     m_avg_frohlich=zero
     saddle_warn=.false.

     !Initializations for the diagonalization routine 
     if(deg_dim>1)then

       ABI_ALLOCATE(eigenvec,(deg_dim,deg_dim))
       lwork=-1
       ABI_ALLOCATE(rwork,(3*deg_dim-2))
       ABI_ALLOCATE(work,(1))
       call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
       lwork=int(work(1))
       ABI_DEALLOCATE(work)
       ABI_ALLOCATE(work,(lwork))

     endif

     !Perform the integral over the sphere
     do itheta=1,ntheta
       costh=gq_points_costh(itheta) ; sinth=gq_points_sinth(itheta)
       weight=gq_weights_th(itheta)*weight_phi
       do iphi=1,nphi
         cosph=gq_points_cosph(iphi) ; sinph=gq_points_sinph(iphi)

         unit_r(1)=sinth*cosph
         unit_r(2)=sinth*sinph
         unit_r(3)=costh

         do iband=1,deg_dim
           do jband=1,deg_dim
             f3d(iband,jband)=DOT_PRODUCT(unit_r,MATMUL(eig2_diag_cart(:,:,iband,jband),unit_r))
           enddo
         enddo
 
         if(deg_dim==1)then
           eigenval(1)=f3d(1,1)
         else
           eigenvec = f3d ; eigenval = zero
           work=zero      ; rwork=zero
           call zheev('V','U',deg_dim,eigenvec,deg_dim,eigenval,work,lwork,rwork,info)
         endif

         m_avg = m_avg + weight*sinth*eigenval
         m_avg_frohlich = m_avg_frohlich + weight*sinth/(abs(eigenval))**half

         if(itheta==1 .and. iphi==1) start_eigf3d_pos = (eigenval > 0)
         do iband=1,deg_dim
           if(start_eigf3d_pos(iband) .neqv. (eigenval(iband)>0)) then
             saddle_warn(iband)=.true.
           end if
         end do

       enddo
     enddo

     if(deg_dim>1)then
       ABI_DEALLOCATE(eigenvec)
       ABI_DEALLOCATE(rwork)
       ABI_DEALLOCATE(work)
     endif

     m_avg = quarter/pi*m_avg
     m_avg = one/m_avg

     m_avg_frohlich = quarter/pi*m_avg_frohlich
     m_avg_frohlich = m_avg_frohlich**2

     if(deg_dim==1)then
       write(ab_out,'(2a,3(f6.3,a),i5)')ch10,&
&        ' - At k-point (',kpt(1),',',kpt(2),',',kpt(3),'), band ',&
&        efmasdeg(ikpt)%degs_bounds(1,ideg)
     else
       write(ab_out,'(2a,3(f6.3,a),i5,a,i5)')ch10,&
&        ' - At k-point (',kpt(1),',',kpt(2),',',kpt(3),'), bands ',&
&        efmasdeg(ikpt)%degs_bounds(1,ideg),' through ',efmasdeg(ikpt)%degs_bounds(2,ideg)
     endif

     sign_warn=.false.
     do iband=1,deg_dim
       if(saddle_warn(iband)) then
         write(ab_out,'(a,i5,a)') ' Band ',efmasdeg(ikpt)%degs_bounds(1,ideg)+iband-1,&
&          ' SADDLE POINT - Frohlich effective mass cannot be defined. '
         sign_warn=.true.
       else
         m_avg_frohlich(iband) = DSIGN(m_avg_frohlich(iband),m_avg(iband))
         write(ab_out,'(a,i5,a,f14.10)') &
&          ' Band ',efmasdeg(ikpt)%degs_bounds(1,ideg)+iband-1,&
&          ' Angular average effective mass for Frohlich model (<m**0.5>)**2= ',m_avg_frohlich(iband)
       endif
       if(start_eigf3d_pos(iband) .neqv. start_eigf3d_pos(1))then
         sign_warn=.true.
       endif
     enddo

     if(sign_warn .eqv. .false.)then
       write(ab_out,'(2a)')&
&       ' Angular and band average effective mass for Frohlich model.'
       write(ab_out,'(a,es16.6)') &
&       ' Value of     (<<m**0.5>>)**2 = ',(sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim)**2
       write(ab_out,'(a,es16.6)') &
&       ' Absolute Value of <<m**0.5>> = ', sum(abs(m_avg_frohlich(1:deg_dim))**0.5)/deg_dim
     else
       write(ab_out,'(a)')& 
&          ' Angular and band average effective mass for Frohlich model cannot be defined because of a sign problem.'
     endif

     ABI_DEALLOCATE(eig2_diag_cart)
     ABI_DEALLOCATE(f3d)
     ABI_DEALLOCATE(m_avg)
     ABI_DEALLOCATE(m_avg_frohlich)
     ABI_DEALLOCATE(eigenval)
     ABI_DEALLOCATE(saddle_warn)
     ABI_DEALLOCATE(start_eigf3d_pos)

   enddo ! ideg
 enddo ! ikpt

 ABI_DEALLOCATE(gq_points_th)
 ABI_DEALLOCATE(gq_points_costh)
 ABI_DEALLOCATE(gq_points_sinth)
 ABI_DEALLOCATE(gq_weights_th)
 ABI_DEALLOCATE(gq_points_cosph)
 ABI_DEALLOCATE(gq_points_sinph)

 end subroutine frohlichmodel

end module m_frohlichmodel
!!***
