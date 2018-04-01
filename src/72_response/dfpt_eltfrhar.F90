!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_eltfrhar
!! NAME
!! dfpt_eltfrhar
!!
!! FUNCTION
!! Compute the frozen-wavefunction hartree enegy contribution to the
!! elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, DCA, XG, GM, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  gsqcut =Fourier cutoff on G^2 for "large sphere" of radius double
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials
!!  mpi_enreg=informations about MPI parallelization
!!  nfft =(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=total electron density in G space
!!
!! OUTPUT
!!  eltfrhar(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!! NOTES
!! *based largely on hartre.f
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_eltfrhar(eltfrhar,rprimd,gsqcut,mpi_enreg,nfft,ngfft,rhog)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_time,     only : timab
 use m_geometry, only : metric
 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_eltfrhar'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(out) :: eltfrhar(6,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,id2,id3,ierr,ig,ig2,ig3,ii,ii1,ing,istr1,istr2,jj
 integer :: ka,kb,kd,kg,me_fft,n1,n2,n3,nproc_fft
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,d2eacc,d2etot,d2gs,deacc01,deacc10,dgs01,dgs10,eacc,fact,gs
 real(dp) :: term,ucvol
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: d2gm(3,3),dgm01(3,3),dgm10(3,3),gmet(3,3),gprimd(3,3),gqr(3)
 real(dp) :: rmet(3,3),tsec(2)
 real(dp),allocatable :: gq(:,:)

! *************************************************************************

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 eltfrhar(:,:)=0.0_dp

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 
!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Initialize a few quantities
 fact=0.5_dp*ucvol/pi
 cutoff=gsqcut*tolfix

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig
   end do
 end do

!Loop over 2nd strain index
 do istr2=1,6
!  Loop over 1st strain index, upper triangle only
   do istr1=1,istr2

     ka=idx(2*istr1-1);kb=idx(2*istr1);kg=idx(2*istr2-1);kd=idx(2*istr2)

     do ii = 1,3
       dgm01(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
       dgm10(:,ii)=-(gprimd(kg,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kg,ii))
     end do

     d2gm(:,:)=0._dp
     do ii = 1,3
       if(ka==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kb,ii)
       if(ka==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(kb,ii)
       if(kb==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(ka,ii)
       if(kb==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(ka,ii)
     end do
     d2gm(:,:)=0.5_dp*d2gm(:,:)

!    initialize energy accumulator
     eacc=0._dp
     deacc01=0._dp
     deacc10=0._dp
     d2eacc=0._dp

     id2=n2/2+2
     id3=n3/2+2
!    Triple loop on each dimension
     do i3=1,n3
       ig3=i3-(i3/id3)*n3-1
       gqr(3)=gq(3,i3)
       do i2=1,n2
         if (fftn2_distrib(i2)==me_fft) then
           gqr(2)=gq(2,i2)
           ig2=i2-(i2/id2)*n2-1
           i23=n1*(ffti2_local(i2)-1 +(n2/nproc_fft)*(i3-1))
!          Do the test that eliminates the Gamma point outside
!          of the inner loop
           ii1=1
           if(i23==0 .and. ig2==0 .and. ig3==0)then
             ii1=2
           end if

!          Final inner loop on the first dimension
!          (note the lower limit)
           do i1=ii1,n1
             gqr(1)=gq(1,i1)
             gs=(gmet(1,1)*gqr(1)*gqr(1)+gmet(2,2)*gqr(2)*gqr(2)+&
&             gmet(3,3)*gqr(3)*gqr(3)+2._dp*&
&             (gmet(1,2)*gqr(1)*gqr(2) + gmet(1,3)*gqr(1)*gqr(3)+&
&             gmet(2,3)*gqr(2)*gqr(3)) )
             ii=i1+i23
             if(gs<=cutoff)then
               dgs01=(dgm01(1,1)*gqr(1)*gqr(1)+dgm01(2,2)*gqr(2)*gqr(2)+&
&               dgm01(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (dgm01(1,2)*gqr(1)*gqr(2) + dgm01(1,3)*gqr(1)*gqr(3)+&
&               dgm01(2,3)*gqr(2)*gqr(3)) )
               dgs10=(dgm10(1,1)*gqr(1)*gqr(1)+dgm10(2,2)*gqr(2)*gqr(2)+&
&               dgm10(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (dgm10(1,2)*gqr(1)*gqr(2) + dgm10(1,3)*gqr(1)*gqr(3)+&
&               dgm10(2,3)*gqr(2)*gqr(3)) )
               d2gs =(d2gm(1,1)*gqr(1)*gqr(1)+d2gm(2,2)*gqr(2)*gqr(2)+&
&               d2gm(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (d2gm(1,2)*gqr(1)*gqr(2) + d2gm(1,3)*gqr(1)*gqr(3)+&
&               d2gm(2,3)*gqr(2)*gqr(3)) )

               term=(rhog(re,ii)**2+rhog(im,ii)**2)/gs
               eacc=eacc+term
               deacc01=deacc01+dgs01*term/gs
               deacc10=deacc10+dgs10*term/gs
               d2eacc=d2eacc+(-d2gs+2._dp*dgs01*dgs10/gs)*term/gs
             end if

!            End loop on i1
           end do
         end if
!        End loop on i2
       end do
!      End loop on i3
     end do

!    Add contributions taking account diagonal strain terms (from ucvol
!    derivatives)
     d2etot=d2eacc
     if(istr1<=3) d2etot=d2etot+deacc10
     if(istr2<=3) d2etot=d2etot+deacc01
     if(istr1<=3 .and. istr2<=3) d2etot=d2etot+eacc

     eltfrhar(istr1,istr2)=fact*d2etot

!    End loop on istr1
   end do
!  End loop in istr2
 end do

 ABI_DEALLOCATE(gq)

!Init mpi_comm
 call timab(48,1,tsec)
 call xmpi_sum(eltfrhar,mpi_enreg%comm_fft,ierr)
 call timab(48,2,tsec)
 
!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrhar(jj,ii)=eltfrhar(ii,jj)
   end do
 end do
end subroutine dfpt_eltfrhar
!!***
