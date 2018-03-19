!{\src2tex{textfont=tt}}
!!****f* ABINIT/strhar
!!
!! NAME
!! strhar
!!
!! FUNCTION
!! Compute Hartree energy contribution to stress tensor (Cartesian coordinates).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ehart=Hartree energy (hartree)
!!  gsqcut=cutoff value on $G^2$ for (large) sphere inside fft box.
!!  $gsqcut=(boxcut^2)*ecut/(2._dp*(\pi^2))$
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhog(2,nfft)= optional argument: Fourier transform of a second charge density (bohr^-3)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  harstr(6)=components of Hartree part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/bohr^3
!!   Definition of symmetric tensor storage: store 6 unique components
!!   in the order 11, 22, 33, 32, 31, 21 (suggested by Xavier Gonze).
!!
!! PARENTS
!!      stress
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine strhar(ehart,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,rprimd,&
&                 rhog2) ! optional argument

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_mpinfo,     only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strhar'
 use interfaces_18_timing
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: ehart,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),rhog(2,nfft)
 real(dp),intent(in),optional :: rhog2(2,nfft)
 real(dp),intent(out) :: harstr(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,id1,id2,id3,ierr,ig1,ig2,ig3,ii,irho2,me_fft,n1,n2,n3,nproc_fft
 real(dp) :: cutoff,gsquar,rhogsq,tolfix=1.000000001_dp,ucvol
!arrays
 real(dp) :: gcart(3),gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 call timab(568,1,tsec)

 harstr(:)=zero
!ehtest=0.0_dp (used for testing)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 irho2=0;if (present(rhog2)) irho2=1

!Conduct looping over all fft grid points to find G vecs inside gsqcut
!Include G**2 on surface of cutoff sphere as well as inside:
 cutoff=gsqcut*tolfix
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2
 ii=0

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2)==me_fft) then
       do i1=1,n1
         ig1=i1-(i1/id1)*n1-1
!        ii=ii+1
         ii=i1+n1*(ffti2_local(i2)-1+(n2/nproc_fft)*(i3-1))
!        **     GET RID OF THIS IF STATEMENT LATER for speed if needed
!        Avoid G=0:
!        if (ii>1) then
         if (ig1==0 .and. ig2==0 .and. ig3==0) cycle
!        Compute cartesian components of G
         gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+gprimd(1,3)*dble(ig3)
         gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+gprimd(2,3)*dble(ig3)
         gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+gprimd(3,3)*dble(ig3)
!        Compute |G|^2
         gsquar=gcart(1)**2+gcart(2)**2+gcart(3)**2

!        Keep only G**2 inside larger cutoff (not sure this is needed):
         if (gsquar<=cutoff) then
!          take |rho(G)|^2 for complex rhog
           if (irho2==0) then
             rhogsq=rhog(re,ii)**2+rhog(im,ii)**2
           else
             rhogsq=rhog(re,ii)*rhog2(re,ii)+rhog(im,ii)*rhog2(im,ii)
           end if
           harstr(1)=harstr(1)+(rhogsq/gsquar**2)*gcart(1)*gcart(1)
           harstr(2)=harstr(2)+(rhogsq/gsquar**2)*gcart(2)*gcart(2)
           harstr(3)=harstr(3)+(rhogsq/gsquar**2)*gcart(3)*gcart(3)
           harstr(4)=harstr(4)+(rhogsq/gsquar**2)*gcart(3)*gcart(2)
           harstr(5)=harstr(5)+(rhogsq/gsquar**2)*gcart(3)*gcart(1) 
           harstr(6)=harstr(6)+(rhogsq/gsquar**2)*gcart(2)*gcart(1)
         end if
!        end if
       end do
     end if
   end do
 end do

!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
#ifdef FC_IBM
 write(std_out,*)' strhar : before mpi_comm, harstr=',harstr
#endif

!Init mpi_comm
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(harstr,mpi_enreg%comm_fft ,ierr)
   call timab(48,2,tsec)
 end if

#ifdef FC_IBM
!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
 write(std_out,*)' strhar : after mpi_comm, harstr=',harstr
 write(std_out,*)' strhar : ehart,ucvol=',ehart,ucvol
#endif

!Normalize and add term -ehart/ucvol on diagonal
 harstr(1)=harstr(1)/pi-ehart/ucvol
 harstr(2)=harstr(2)/pi-ehart/ucvol
 harstr(3)=harstr(3)/pi-ehart/ucvol
 harstr(4)=harstr(4)/pi
 harstr(5)=harstr(5)/pi
 harstr(6)=harstr(6)/pi

 call timab(568,2,tsec)

end subroutine strhar
!!***
