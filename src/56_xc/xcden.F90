!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcden
!! NAME
!! xcden
!!
!! FUNCTION
!! Prepare density data before calling local or semi-local xc evaluators.
!!
!! NOTES
!! Can take into account a shift of the grid, for purpose of better accuracy.
!! Can also compute the gradient of the density, for use with GGAs.
!! Also eliminate eventual negative densities.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  ishift : if ==0, do not shift the xc grid (usual case); if ==1, shift the xc grid
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngrad : =1, only compute the density ; =2 also compute the
!!      gradient of the density. Note : ngrad**2 is also used to dimension rhonow
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(cplex*nfft,nspden)=real space electron density in electrons/bohr**3, on the
!!   unshifted grid (total in first half and spin-up in second half if nspden=2)
!!
!! OUTPUT
!!  rhonow(cplex*nfft,nspden,ngrad*ngrad)=electron (spin)-density in real space and
!!     eventually its gradient, either on the unshifted grid (if ishift==0,
!!     then equal to rhor),or on the shifted grid
!!    rhonow(:,:,1)=electron density in electrons/bohr**3
!!    if ngrad==2 : rhonow(:,:,2:4)=gradient of electron density in electrons/bohr**4
!!  OPTIONAL OUTPUT
!!  lrhonow(cplex*nfft,nspden)=Laplacian of the electron (spin)-density in real space
!!    in electrons/bohr**5 (in case of meta GGA)
!!
!! PARENTS
!!      afterscfloop,dfpt_mkvxcgga,dfpt_mkvxcstrgga,gammapositron_fft
!!      rhohxcpositron,rhotoxc
!!
!! CHILDREN
!!      fourdp,phase,ptabs_fourdp,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor,rhonow, & !Mandatory arguments
&  lrhonow)              !Optional arguments

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_time,     only : timab
 use m_mpinfo,   only : ptabs_fourdp
 use m_fft_mesh, only : phase

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcden'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ishift,nfft,ngrad,nspden,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
 real(dp),intent(out),optional :: lrhonow(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,idir,ifft,ig1,ig2,ig3
 integer :: ispden,n1,n2,n3,qeq0
 real(dp) :: gc23_idir,gcart_idir,ph123i,ph123r,ph1i,ph1r,ph23i,ph23r,ph2i,ph2r
 real(dp) :: ph3i,ph3r,work_im,work_re
 character(len=500) :: message
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: gcart1(:),gcart2(:),gcart3(:),ph1(:),ph2(:),ph3(:)
 real(dp),allocatable :: wkcmpx(:,:),work(:),workgr(:,:),worklp(:,:)

! *************************************************************************

!DEBUG
!write(std_out,*)' xcden : enter '
!ENDDEBUG

 if (ishift/=0 .and. ishift/=1) then
   write(message, '(a,i0)' )'ishift must be 0 or 1 ; input was',ishift
   MSG_BUG(message)
 end if

 if (ngrad/=1 .and. ngrad/=2) then
   write(message, '(a,i0)' )'ngrad must be 1 or 2 ; input was',ngrad
   MSG_BUG(message)
 end if

!Keep local copy of fft dimensions
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Initialize computation of G in cartesian coordinates
 id1=n1/2+2  ; id2=n2/2+2  ; id3=n3/2+2

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Check whether q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

 if(ishift==0)then

!  Copy the input rhor in the new location. Will check on negative values later

   do ispden=1,nspden
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,nfft,rhonow,rhor)
     do ifft=1,cplex*nfft
       rhonow(ifft,ispden,1)=rhor(ifft,ispden)
     end do
   end do

 end if

 if(ishift==1 .or. ngrad==2)then

   ABI_ALLOCATE(wkcmpx,(2,nfft))
   ABI_ALLOCATE(work,(cplex*nfft))

   if(ishift==1)then
!    Precompute phases (The phases correspond to a shift of density on real space
!    grid from center at 0 0 0 to (1/2)*(1/n1,1/n2,1/n3).)
     ABI_ALLOCATE(ph1,(2*n1))
     ABI_ALLOCATE(ph2,(2*n2))
     ABI_ALLOCATE(ph3,(2*n3))
     call phase(n1,ph1)
     call phase(n2,ph2)
     call phase(n3,ph3)
   end if

   do ispden=1,nspden

!    Obtain rho(G) in wkcmpx from input rho(r)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,nfft,rhor,work)
     do ifft=1,cplex*nfft
       work(ifft)=rhor(ifft,ispden)
     end do

     call timab(82,1,tsec)
     call fourdp(cplex,wkcmpx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call timab(82,2,tsec)

!    If shift is required, multiply now rho(G) by phase, then generate rho(r+delta)
     if(ishift==1)then
!$OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3,ph1i,ph1r,ph123i,ph123r,ph2i,ph2r,ph23i,ph23r,ph3i,ph3r,work_im,work_re) &
!$OMP&SHARED(n1,n2,n3,ph1,ph2,ph3,wkcmpx,mpi_enreg,fftn2_distrib)
       do i3=1,n3
         ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
         ph3r=ph3(2*i3-1)
         ph3i=ph3(2*i3  )
         do i2=1,n2
           ph2r=ph2(2*i2-1)
           ph2i=ph2(2*i2  )
           ph23r=ph2r*ph3r-ph2i*ph3i
           ph23i=ph2i*ph3r+ph2r*ph3i
           if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
             do i1=1,n1
               ifft=ifft+1
               ph1r=ph1(2*i1-1)
               ph1i=ph1(2*i1  )
               ph123r=ph1r*ph23r-ph1i*ph23i
               ph123i=ph1i*ph23r+ph1r*ph23i
!              Must use intermediate variables !
               work_re=ph123r*wkcmpx(1,ifft)-ph123i*wkcmpx(2,ifft)
               work_im=ph123i*wkcmpx(1,ifft)+ph123r*wkcmpx(2,ifft)
               wkcmpx(1,ifft)=work_re
               wkcmpx(2,ifft)=work_im
             end do
           end if
         end do
       end do
       call timab(82,1,tsec)
       call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       call timab(82,2,tsec)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,rhonow,work)
       do ifft=1,cplex*nfft
         rhonow(ifft,ispden,1)=work(ifft)
       end do
     end if

!    If gradient of the density is required, take care of the three components now
!    Note : this operation is applied on the eventually shifted rho(G)
     if(ngrad==2)then
       ABI_ALLOCATE(gcart1,(n1))
       ABI_ALLOCATE(gcart2,(n2))
       ABI_ALLOCATE(gcart3,(n3))
       ABI_ALLOCATE(workgr,(2,nfft))
       if (present(lrhonow)) then
         ABI_ALLOCATE(worklp,(2,nfft))
         lrhonow(:,ispden)=zero
       end if
       do idir=1,3

         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           gcart1(i1)=gprimd(idir,1)*two_pi*(dble(ig1)+qphon(1))
         end do
!        Note that the G <-> -G symmetry must be maintained
         if(mod(n1,2)==0 .and. qeq0==1)gcart1(n1/2+1)=zero
         do i2=1,n2
           ig2=i2-(i2/id2)*n2-1
           gcart2(i2)=gprimd(idir,2)*two_pi*(dble(ig2)+qphon(2))
         end do
         if(mod(n2,2)==0 .and. qeq0==1)gcart2(n2/2+1)=zero
         do i3=1,n3
           ig3=i3-(i3/id3)*n3-1
           gcart3(i3)=gprimd(idir,3)*two_pi*(dble(ig3)+qphon(3))
         end do
         if(mod(n3,2)==0 .and. qeq0==1)gcart3(n3/2+1)=zero

!        MG: Be careful here with OMP due to ifft. Disabled for the time being.
!        !$OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3,gcart_idir,gc23_idir) &
!        !$OMP&SHARED(gcart1,gcart2,gcart3,n1,n2,n3,wkcmpx,workgr)
         ifft = 0
         do i3=1,n3
           do i2=1,n2
             gc23_idir=gcart2(i2)+gcart3(i3)
             if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
               do i1=1,n1
                 ifft=ifft+1
                 gcart_idir=gc23_idir+gcart1(i1)
!                Multiply by i 2pi G(idir)
                 workgr(2,ifft)= gcart_idir*wkcmpx(1,ifft)
                 workgr(1,ifft)=-gcart_idir*wkcmpx(2,ifft)
!                Do the same to the gradient in order to get the laplacian
                 if (present(lrhonow)) then
                   worklp(2,ifft)= gcart_idir*workgr(1,ifft)
                   worklp(1,ifft)=-gcart_idir*workgr(2,ifft)
                 end if
               end do
             end if
           end do
         end do
         call timab(82,1,tsec)
         call fourdp(cplex,workgr,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
         call timab(82,2,tsec)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(idir,ispden,cplex,nfft,rhonow,work)
         do ifft=1,cplex*nfft
           rhonow(ifft,ispden,1+idir)=work(ifft)
         end do

         if (present(lrhonow)) then
           call fourdp(cplex,worklp,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
           do ifft=1,cplex*nfft
             lrhonow(ifft,ispden)=lrhonow(ifft,ispden)+work(ifft)
           end do
         end if

       end do
       ABI_DEALLOCATE(gcart1)
       ABI_DEALLOCATE(gcart2)
       ABI_DEALLOCATE(gcart3)
       ABI_DEALLOCATE(workgr)
       if (allocated(worklp))  then
         ABI_DEALLOCATE(worklp)
       end if
     end if

   end do  ! End loop on spins

   ABI_DEALLOCATE(wkcmpx)
   ABI_DEALLOCATE(work)
   if(ishift==1) then
     ABI_DEALLOCATE(ph1)
     ABI_DEALLOCATE(ph2)
     ABI_DEALLOCATE(ph3)
   end if

 end if  ! End condition on ishift and ngrad

end subroutine xcden
!!***
