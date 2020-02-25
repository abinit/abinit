!!****m* ABINIT/m_xctk
!! NAME
!!  m_xctk
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, DRH)
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

module m_xctk

 use defs_basis
 use m_abicore
 use m_errors

 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_mpinfo,   only : ptabs_fourdp
 use m_fft_mesh, only : phase
 use m_fft,      only : fourdp

 implicit none

 private
!!***

 public :: xcden
 public :: xcpot
!!***

contains
!!***

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

subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,qphon,rhor,rhonow, & !Mandatory arguments
&  lrhonow)              !Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ishift,nfft,ngrad,nspden
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
     call fourdp(cplex,wkcmpx,work,-1,mpi_enreg,nfft,1,ngfft,0)
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
       call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,1,ngfft,0)
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
         call fourdp(cplex,workgr,work,1,mpi_enreg,nfft,1,ngfft,0)
         call timab(82,2,tsec)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(idir,ispden,cplex,nfft,rhonow,work)
         do ifft=1,cplex*nfft
           rhonow(ifft,ispden,1+idir)=work(ifft)
         end do

         if (present(lrhonow)) then
           call fourdp(cplex,worklp,work,1,mpi_enreg,nfft,1,ngfft,0)
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
   if (allocated(wkcmpx))  then
     ABI_DEALLOCATE(wkcmpx)
   end if
   ABI_DEALLOCATE(work)
   if(ishift==1) then
     ABI_DEALLOCATE(ph1)
     ABI_DEALLOCATE(ph2)
     ABI_DEALLOCATE(ph3)
   end if

 end if  ! End condition on ishift and ngrad

end subroutine xcden
!!***

!!****f* ABINIT/xcpot
!! NAME
!! xcpot
!!
!! FUNCTION
!! Process data after calling local or semi-local xc evaluators
!! The derivative of Exc with respect to the density, gradient of density
!! in case of GGAs, and Laplacian of density in case of meta-GGA
!! are contained in depsxc(:,:).
!! In case of GGAs (and meta-GGAs) the gradient of the density contained
!! in rhonow(:,:,2:4) is already multiplied by the local partial derivative
!! of the XC functional.
!! Can take into account a shift of the grid, for purpose of better accuracy
!!
!! INPUTS
!!  cplex=if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  [depsxc(cplex*nfft,nspgrad)]=derivative of Exc with respect to the (spin-)density,
!!    or to the norm of the gradient of the (spin-)density,
!!    further divided by the norm of the gradient of the (spin-)density
!!   The different components of depsxc will be
!!   for nspden=1,             depsxc(:,1)=d(rho.exc)/d(rho)
!!     and if ngrad=2,         depsxc(:,2)=1/2*1/|grad rho_up|*d(n.exc)/d(|grad rho_up|)
!!                                     +1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!     and if use_laplacian=1, depsxc(:,3)=d(rho.exc)/d(lapl rho)
!!   for nspden>=2,            depsxc(:,1)=d(rho.exc)/d(rho_up)
!!                             depsxc(:,2)=d(rho.exc)/d(rho_down)
!!     and if ngrad=2,         depsxc(:,3)=1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!!                             depsxc(:,4)=1/|grad rho_down|*d(rho.exc)/d(|grad rho_down|)
!!                             depsxc(:,5)=1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!     and if use_laplacian=1, depsxc(:,6)=d(rho.exc)/d(lapl rho_up)
!!                             depsxc(:,7)=d(rho.exc)/d(lapl rho_down)
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  ishift : if ==0, do not shift the xc grid (usual case);
!!           if ==1, shift the xc grid
!!  use_laplacian : 1 if we use a  functional depending on the laplacian of the density
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngrad : =1, only take into account derivative wrt the density ;
!!          =2, also take into account derivative wrt the gradient of the density.
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  [rhonow(cplex*nfft,nspden,ngrad*ngrad)]=electron (spin)-density in real space and
!!     eventually its gradient already multiplied by the local partial derivative
!!     of the XC functional, either on the unshifted grid (if ishift==0,
!!     then equal to rhor), or on the shifted grid
!!    rhonow(:,:,1)=electron density in electrons/bohr**3
!!    if ngrad==2 : rhonow(:,:,2:4)=gradient of electron density in el./bohr**4,
!!     times local partial derivative of the functional, as required by the GGA
!!    In this routine, rhonow is used only in the GGA case (ngrad=2).
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output (all optional:
!!  [vxc(cplex*nfft,nspden)]=xc potential (spin up in first half and spin down in
!!   second half if nspden>=2). Contribution from the present shifted
!!   or unshifted grid is ADDED to the input vxc data.
!!  [vxctau(cplex*nfft,nspden,4)]=derivative of XC energy density with respect to
!!   kinetic energy density (depsxcdtau). The arrays vxctau(nfft,nspden,4) contains also
!!   the gradient of vxctau (gvxctau) which will be computed here in vxctau(:,:,2:4).
!!
!! PARENTS
!!      dfpt_mkvxcgga,dfpt_mkvxcstrgga,rhotoxc
!!
!! CHILDREN
!!      fourdp,phase,ptabs_fourdp,timab
!!
!! SOURCE

subroutine xcpot (cplex,gprimd,ishift,use_laplacian,mpi_enreg,nfft,ngfft,ngrad,nspden,&
&                 nspgrad,qphon,&
&                 depsxc,rhonow,vxc,vxctau) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ishift,nfft,ngrad,nspden,nspgrad,use_laplacian
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in),optional :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
 real(dp),intent(in),optional :: depsxc(cplex*nfft,nspgrad),gprimd(3,3),qphon(3)
 real(dp),intent(inout),optional :: vxc(cplex*nfft,nspden)
 real(dp),intent(inout),optional :: vxctau(cplex*nfft,nspden,4)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,idir,ifft,ig1,ig2,ig3,ispden,n1,n2,n3,qeq0
 real(dp),parameter :: lowden=1.d-14,precis=1.d-15
 real(dp) :: gc23_idir,gcart_idir,ph123i,ph123r,ph1i,ph1r,ph23i,ph23r,ph2i,ph2r
 real(dp) :: ph3i,ph3r,work_im,work_re
 character(len=500) :: message
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 logical :: with_vxc,with_vxctau
 real(dp) :: tsec(2)
 real(dp),allocatable :: gcart1(:),gcart2(:),gcart3(:),ph1(:),ph2(:),ph3(:)
 real(dp),allocatable :: wkcmpx(:,:),wkcmpxtau(:,:)
 real(dp),allocatable :: work(:),workgr(:,:),worklp(:,:),worktau(:,:)

! *************************************************************************

 if (ishift/=0 .and. ishift/=1) then
   write(message, '(a,i0)' )' ishift must be 0 or 1 ; input was',ishift
   MSG_BUG(message)
 end if

 if (ngrad/=1 .and. ngrad/=2 ) then
   write(message, '(a,i0)' )' ngrad must be 1 or 2 ; input was',ngrad
   MSG_BUG(message)
 end if

 with_vxc=present(vxc) ; with_vxctau=present(vxctau)
 if (with_vxc) then
   if ((.not.present(rhonow)).or.(.not.present(depsxc))) then
     message='need rhonow or depsxc!'
     MSG_BUG(message)
   end if
 end if
 
!Keep local copy of fft dimensions
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Initialize computation of G in cartesian coordinates
 id1=n1/2+2  ; id2=n2/2+2  ; id3=n3/2+2

 !Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Check whether q=0
 qeq0=0;if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

 if(with_vxc.and.ishift==0)then ! Add the value of depsxc to vxc
   do ispden=1,min(nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,depsxc,nfft,vxc,ispden)
     do ifft=1,cplex*nfft
       vxc(ifft,ispden)=vxc(ifft,ispden)+depsxc(ifft,ispden)
     end do
   end do
 end if

!If the grid is shifted, or if gradient corrections are present, there must be FFTs.
 if(ishift==1 .or. ngrad==2)then

   if (with_vxc) then
     ABI_ALLOCATE(work,(cplex*nfft))
     ABI_ALLOCATE(wkcmpx,(2,nfft))
   end if

   if(ishift==1)then
     ABI_ALLOCATE(ph1,(2*n1))
     ABI_ALLOCATE(ph2,(2*n2))
     ABI_ALLOCATE(ph3,(2*n3))
!    Precompute phases (The phases correspond to a shift of density on real space
!    grid from center at 0 0 0 to (1/2)*(1/n1,1/n2,1/n3).)
     call phase(n1,ph1)
     call phase(n2,ph2)
     call phase(n3,ph3)
   end if

   do ispden=1,min(nspden,2)

!    Initialize wkcmpx either to 0 or to the shifted vxc value
     if (with_vxc) then
       if(ishift==0)then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(nfft,wkcmpx)
         do ifft=1,nfft
           wkcmpx(:,ifft)=zero
         end do
       else
!      Obtain depsxc(G)*phase in wkcmpx from input depsxc(r+delta)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,depsxc,ispden,nfft,work)
         do ifft=1,cplex*nfft
           work(ifft)=depsxc(ifft,ispden)
         end do
         call timab(82,1,tsec)
         call fourdp(cplex,wkcmpx,work,-1,mpi_enreg,nfft,1,ngfft,0)
         call timab(82,2,tsec)
       end if
     end if

!    If gradient correction is present, take care of the three components now
!    Note : this operation is done on the eventually shifted grid
     if (ngrad==2) then
       ABI_ALLOCATE(gcart1,(n1))
       ABI_ALLOCATE(gcart2,(n2))
       ABI_ALLOCATE(gcart3,(n3))
       if (with_vxc) then
         ABI_ALLOCATE(workgr,(2,nfft))
         if (use_laplacian==1) then
           ABI_ALLOCATE(worklp,(2,nfft))
         end if
      end if
      if  (with_vxctau)  then
        ABI_ALLOCATE(worktau,(2,nfft))
        ABI_ALLOCATE(wkcmpxtau,(2,nfft))
      end if

       do idir=1,3
       
         if (with_vxc) then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,idir,ispden,nfft,rhonow,work)
           do ifft=1,cplex*nfft
             work(ifft)=rhonow(ifft,ispden,1+idir)
           end do
           call timab(82,1,tsec)
           call fourdp(cplex,workgr,work,-1,mpi_enreg,nfft,1,ngfft,0)
           call timab(82,2,tsec)

!          IF Meta-GGA then take care of the laplacian term involved.
!          Note : this operation is done on the eventually shifted grid
           if(use_laplacian==1)then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,idir,ispden,nspden,nfft,depsxc,work)
             do ifft=1,cplex*nfft
               if(nspden==1)then
                 work(ifft)=depsxc(ifft,2+ispden)
               else if(nspden==2)then
                 work(ifft)=depsxc(ifft,5+ispden)
               end if
             end do
             call timab(82,1,tsec)
             call fourdp(cplex,worklp,work,-1,mpi_enreg,nfft,1,ngfft,0)
             call timab(82,2,tsec)
           end if
         end if

         if(with_vxctau)then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ispden,nfft,vxctau,work)
           do ifft=1,cplex*nfft
             work(ifft)=vxctau(ifft,ispden,1)
           end do
           call timab(82,1,tsec)
           call fourdp(cplex,worktau,work,-1,mpi_enreg,nfft,1,ngfft,0)
           call timab(82,2,tsec)
         end if ! present vxctau

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

!        !$OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3,gc23_idir,gcart_idir) &
!        !$OMP&SHARED(gcart1,gcart2,gcart3,n1,n2,n3,wkcmpx,workgr)
         ifft = 0
         do i3=1,n3
           do i2=1,n2
             gc23_idir=gcart2(i2)+gcart3(i3)
             if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
               do i1=1,n1
                 ifft=ifft+1
                 gcart_idir=gc23_idir+gcart1(i1)
                 if(with_vxc)then
!                  Multiply by - i 2pi G(idir) and accumulate in wkcmpx
                   wkcmpx(1,ifft)=wkcmpx(1,ifft)+gcart_idir*workgr(2,ifft)
                   wkcmpx(2,ifft)=wkcmpx(2,ifft)-gcart_idir*workgr(1,ifft)
                   if(use_laplacian==1)then
!                    Multiply by - i 2pi G(idir) and accumulate in wkcmpx
                     wkcmpx(1,ifft)=wkcmpx(1,ifft)-gcart_idir**2*worklp(1,ifft)
                     wkcmpx(2,ifft)=wkcmpx(2,ifft)-gcart_idir**2*worklp(2,ifft)
                   end if
                 end if
                 if(with_vxctau)then
                   wkcmpxtau(1,ifft)= gcart_idir*worktau(2,ifft)
                   wkcmpxtau(2,ifft)=-gcart_idir*worktau(1,ifft)
                 end if
               end do
             end if
           end do
         end do

         if (with_vxctau) then
           call timab(82,1,tsec)
           call fourdp(cplex,wkcmpxtau,work,1,mpi_enreg,nfft,1,ngfft,0)
           call timab(82,2,tsec)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ispden,nfft,vxctau,work)
           do ifft=1,cplex*nfft
             vxctau(ifft,ispden,1+idir)=work(ifft)
           end do
         end if

       end do ! enddo idir

       ABI_DEALLOCATE(gcart1)
       ABI_DEALLOCATE(gcart2)
       ABI_DEALLOCATE(gcart3)
       if (with_vxc) then
         ABI_DEALLOCATE(workgr)
         if (use_laplacian==1) then
           ABI_DEALLOCATE(worklp)
         end if
       end if
       if (with_vxctau) then
         ABI_DEALLOCATE(worktau)
         ABI_DEALLOCATE(wkcmpxtau)
       end if

     end if

!    wkcmpx(:,:) contains now the full exchange-correlation potential, but
!    eventually for the shifted grid

     if (with_vxc) then
       if(ishift==1)then
!        Take away the phase to get depsxc(G)
         ifft=0
         do i3=1,n3
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
!                Multiply by phase.  Must use intermediate variables !
                 work_re= ph123r*wkcmpx(1,ifft)+ph123i*wkcmpx(2,ifft)
                 work_im=-ph123i*wkcmpx(1,ifft)+ph123r*wkcmpx(2,ifft)
                 wkcmpx(1,ifft)=work_re
                 wkcmpx(2,ifft)=work_im
               end do
             end if
           end do
         end do
       end if

       call timab(82,1,tsec)
       call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,1,ngfft,0)
       call timab(82,2,tsec)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(cplex,ispden,nfft,vxc,work)
       do ifft=1,cplex*nfft
         vxc(ifft,ispden)=vxc(ifft,ispden)+work(ifft)
       end do
     end if

   end do ! End loop on spins

   if(ishift==1)  then
     ABI_DEALLOCATE(ph1)
     ABI_DEALLOCATE(ph2)
     ABI_DEALLOCATE(ph3)
   end if
   if(with_vxc) then
     ABI_DEALLOCATE(wkcmpx)
     ABI_DEALLOCATE(work)
   end if

 end if ! End condition on ishift/ngrad

end subroutine xcpot
!!***

end module m_xctk
!!***
