!!****m* ABINIT/m_kxc
!! NAME
!! m_kxc
!!
!! FUNCTION
!! Helper functions to compute the XC kernel in reciprocal space.
!! WARNING: At present (10/01/14) these routines are not tested
!! since the ACFD code has been disabled.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, MF, XG, GMR, LSI, YMN, Rhaltaf, MS)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES:
!!  libxc_functionals.F90 uses a global structure (funcs) to store the XC parameters.
!!  This structure is initialized in driver with the value of ixc specified by the user in the input file.
!!  In order to change the value of ixc at run-time, we have to reinitialize the global structure
!!  with the new value of ixc before computing XC quantities.
!!  Moreover one has to reinstate the old functional before returning so that the other routines
!!  will continue to used the previous ixc. This task can be accomplished with the following pseudocode
!!
!!   ! Reinitialize the libxc module with the overriden values
!!   if (old_ixc<0)  call libxc_functionals_end()
!!   if (new_ixc<0) call libxc_functionals_init(new_ixc,nspden)
!!   ! Compute XC stuff here.
!!   ! Revert libxc module to the original settings
!!   if (new_ixc<0) call libxc_functionals_end()
!!   if (old_ixc<0) call libxc_functionals_init(old_ixc,nspden)
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_kxc

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_crystal
 use m_distribfft
 use m_xcdata
 use libxc_functionals
 use m_dtset

 use defs_abitypes,   only : MPI_type
 use m_io_tools,      only : open_file
 use m_pptools,       only : printxsf
 use m_numeric_tools, only : hermitianize
 use m_fft_mesh,      only : g2ifft
 use m_fft,           only : fourdp_6d, fourdp
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_spacepar,      only : hartre
 use m_rhotoxc,       only : rhotoxc
 use m_dfpt_mkvxc,    only : dfpt_mkvxc

 implicit none

 private
!!***

 public :: kxc_rpa         ! Hartree kernel
 public :: kxc_local       ! Compute local xc kernel in G space.
 public :: kxc_alda        ! AL(S)DA kernel in reciprocal space, on the FFT grid.
 public :: kxc_pgg         ! Compute the PGG-exchange kernel in reciprocal space (Phys. Rev. Lett. 76, 1212 (1996) [[cite:Petersilka1996]]).
 public :: kxc_eok         ! linear or non-linear (ixceok = 2) energy optimized kernel of Dobson and Wang.
 public :: kxc_driver      ! Driver routine (TODO)
 public :: kxc_ADA         ! Adiabatic density approximation


CONTAINS  !=========================================================================================================================
!!***

!!****f* m_kxc/kxc_rpa
!! NAME
!! kxc_rpa
!!
!! FUNCTION
!! Return the Hartree kernel:
!!  If option = 0, the bare Hartree kernel:
!!   krpa(ipw) = 4.0*pi/gsq(ipw) if gsq(ipw) /= 0.,
!!   krpa(ipw) = 0.0             if gsq(ipw) == 0. (1 <= ipw <= npw).
!!  If option /= 0, the Hartree kernel with a cut-off in real space beyond rcut_coulomb:
!!   krpa(ipw) = (4.0*pi/gsq(ipw))*(1.0-cos(sqrt(gsq(ipw))*rcut_coulomb)) if gsq(ipw) /= 0.,
!!   krpa(ipw) =  2.0*pi*rcut_coulomb**2                                  if gsq(ipw) == 0.
!!
!! INPUTS
!!  gsq(npw) = the squared norm of the planewaves.
!!  npw = number of planewaves in the gsq array.
!!  option = 0 for the bare Hartree kernel, /=0 for the cut-off Hartree kernel.
!!  rcut_coulomb = real space cut-off radius for the Coulomb interaction in Bohr.
!!
!! OUTPUT
!!  krpa(npw) = the Hartree kernel.
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_rpa(gsq,krpa,npw,option,rcut_coulomb)

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: npw,option
 real(dp),intent(in) :: rcut_coulomb
!arrays
 real(dp),intent(in) :: gsq(npw)
 real(dp),intent(out) :: krpa(npw)

!Local variables -------------------------------------------------------
!scalars
 integer :: ipw

!***********************************************************************

 if (option == 0) then

!  Compute the bare Hartree kernel.

   do ipw = 1,npw

     if (gsq(ipw) > tol12) then
       krpa(ipw) = four_pi/gsq(ipw)
     else
       krpa(ipw) = 0._dp
     end if

   end do

 else

!  Compute the Hartree kernel with a cut-off in real space beyond rcut_coulomb:

   do ipw = 1,npw

     if (gsq(ipw) > tol12) then
       krpa(ipw) = (four_pi/gsq(ipw))*(1._dp-cos(sqrt(gsq(ipw))*rcut_coulomb))
     else
       krpa(ipw) = two_pi*rcut_coulomb**2
     end if

   end do

 end if

end subroutine kxc_rpa
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_local
!! NAME
!! kxc_local
!!
!! FUNCTION
!! In a planewave basis set, the matrix of a local xc kernel:
!!  $f_{\rm xc}(\vec{r},\vec{r}') = f(\vec{r})\delta(\vec{r}-\vec{r}')$
!! is just:
!!  $f_{\rm xc}(\vec{G},\vec{G}') = f(\vec{G}-\vec{G}')$.
!! This subroutine calculates the matrix of such a local xc kernel
!! given $f(\vec{G})$ on the FFT grid.
!!
!! INPUTS
!!  ispxc = 1 for the up-up spin channel.
!!        = 2 for the up-down (and down-up) spin channels.
!!        = 3 for the down-down spin channel.
!!  ispxc must be 1 if nspden = 1.
!!  kg_diel(3,npwdiel) = reduced planewave coordinates for the kxc matrix.
!!  kxcg(2,nfft) = $f(\vec{G})$ on the FFT grid.
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  npwdiel = number of planewaves for the susceptibility matrix.
!!  nspden = number of spin-density components.
!!  option = 0 do not compute the first row and column of the matrix of the
!!             xc kernel (which we assume to the G = 0 row and column).
!!        /= 0 compute the full matrix of the xc kernel.
!!
!! OUTPUT
!!  kxc(2,npwdiel,nspden,npwdiel,nspden) = the matrix of the xc kernel.
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_local(ispxc,kg_diel,kxc,kxcg,nfft,ngfft,npwdiel,nspden,option)

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ispxc,nfft,npwdiel,nspden,option
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),ngfft(18)
 real(dp),intent(in) :: kxcg(2,nfft)
 real(dp),intent(out) :: kxc(2,npwdiel,nspden,npwdiel,nspden)

!Local variables -------------------------------------------------------
!For debbuging purposes:
!real(dp) :: c1,c2,c3
!scalars
 integer :: i1,i2,i3,ifft,ipw1,ipw2,ipwstart,isp1,isp2,j1,j2,j3,k1,k2,k3,n1,n2
 integer :: n3
 logical :: ok
 character(len=500) :: message

!***********************************************************************

!Check input parameters.

 if (nspden > 2) then
   message =' kxc_local does not work yet for nspden > 2.'
   MSG_ERROR(message)
 end if

 isp1 = 1
 isp2 = 1
 ok = .true.
 if (nspden == 1) then
   select case (ispxc)
     case (1)
       isp1 = 1
       isp2 = 1
     case default
       ok = .false.
   end select
 else
   select case (ispxc)
     case (1)
       isp1 = 1
       isp2 = 1
     case (2)
       isp1 = 1
       isp2 = 2
     case (3)
       isp1 = 2
       isp2 = 2
     case default
       ok = .false.
   end select
 end if

 if (.not.ok) then
   write (message,'(2(a,i0))')'  The input ispxc = ',ispxc,' is not compatible with nspden = ',nspden
   MSG_BUG(message)
 end if

 if (option == 0) then
   ipwstart = 2
   kxc(:,1,isp1,:,isp2) = 0._dp
   kxc(:,:,isp1,1,isp2) = 0._dp
 else
   ipwstart = 1
 end if

!Calculate the xc matrix.

 n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)

 do ipw2 = ipwstart,npwdiel

   j1 = kg_diel(1,ipw2) ; j2 = kg_diel(2,ipw2) ; j3 = kg_diel(3,ipw2)

!  Fill the diagonal.

   kxc(:,ipw2,isp1,ipw2,isp2) = kxcg(:,1)

!  Fill the off-diagonal elements.

   do ipw1 = ipw2+1,npwdiel

     i1 = kg_diel(1,ipw1) ; i2 = kg_diel(2,ipw1) ; i3 = kg_diel(3,ipw1)

!    Compute the difference between G vectors.
!    The use of two mod calls handles both i1-j1 >= n1 AND i1-j1 < 0.

     k1 = mod(n1+mod(i1-j1,n1),n1)
     k2 = mod(n2+mod(i2-j2,n2),n2)
     k3 = mod(n3+mod(i3-j3,n3),n3)

     ifft = k1+n1*(k2+n2*k3)+1

     kxc(1,ipw1,isp1,ipw2,isp2) =  kxcg(1,ifft)
     kxc(2,ipw1,isp1,ipw2,isp2) =  kxcg(2,ifft)

     kxc(1,ipw2,isp1,ipw1,isp2) =  kxcg(1,ifft)
     kxc(2,ipw2,isp1,ipw1,isp2) = -kxcg(2,ifft)

   end do

 end do

!If needed, copy the up-down to the down-up spin channel.

 if (ispxc == 2) then
   do ipw2 = 1,npwdiel
     do ipw1 = 1,npwdiel
       kxc(1,ipw2,isp2,ipw1,isp1) =  kxc(1,ipw1,isp1,ipw2,isp2)
       kxc(2,ipw2,isp2,ipw1,isp1) = -kxc(2,ipw1,isp1,ipw2,isp2)
     end do
   end do
 end if

!DEBUG
!See kxc_alda.f, "test kernel" DEBUG section.
!do ipw2 = 1,npwdiel
!j1 = kg_diel(1,ipw2) ; j2 = kg_diel(2,ipw2) ; j3 = kg_diel(3,ipw2)
!do ipw1 = ipw2+1,npwdiel
!i1 = kg_diel(1,ipw1) ; i2 = kg_diel(2,ipw1) ; i3 = kg_diel(3,ipw1)
!k1 = mod(n1+mod(i1-j1,n1),n1)
!k2 = mod(n2+mod(i2-j2,n2),n2)
!k3 = mod(n3+mod(i3-j3,n3),n3)
!ifft = k1+n1*(k2+n2*k3)+1
!c1 = 0._dp ; c2 = 0._dp ; c3 = 0._dp
!if (i1-j1 ==  0) c1 = c1+0.0_dp
!if (i2-j2 ==  0) c2 = c2+0.0_dp
!if (i3-j3 ==  0) c3 = c3+0.0_dp
!if (i1-j1 ==  1) c1 = c1+0.5_dp
!if (i2-j2 ==  2) c2 = c2+0.5_dp
!if (i3-j3 ==  3) c3 = c3+0.5_dp
!if (i1-j1 == -1) c1 = c1+0.5_dp
!if (i2-j2 == -2) c2 = c2+0.5_dp
!if (i3-j3 == -3) c3 = c3+0.5_dp
!if ((abs(kxcg(1,ifft)-c1*c2*c3) > tol10).or.(abs(kxcg(2,ifft)) > tol10)) then
!write (std_out,*) ' i1 i2 i3 ifft: ',i1,i2,i3,ifft
!write (std_out,*) ' exp.: ',c1*c2*c3,' got: ',kxcg(:,ifft)
!end if
!end do
!end do
!ENDDEBUG

end subroutine kxc_local
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_alda
!! NAME
!! kxc_alda
!!
!! FUNCTION
!! If option = 1:
!!  Compute the AL(S)DA kernel in reciprocal space, on the FFT grid.
!! If option = 2:
!!  Only computes the up-down channel of the AL(S)DA kernel, on the
!!  FFT grid, for use in the BPG kernel.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  ixc = choice of exchange-correlation functional.
!!  mpi_enreg=information about MPI parallelization
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  nspden = number of spin-density components.
!!  option = 1 compute the AL(S)DA kernel in reciprocal space.
!!         = 2 only computes the up-down channel of the AL(S)DA kernel,
!!             for use in the BPG kernel.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rhocut = cut-off density for the local kernels (ALDA, EOK),
!!           relative to max(rhor(:,:)).
!!  rprimd(3,3) = dimensional primitive translations for real space in Bohr.
!!
!! OUTPUT
!!  kxcg(2,nfft,*) = the AL(S)DA kernel in reciprocal space, on the FFT grid
!!   (the third dimension is 2*nspden-1 if option = 1, and 1 if option = 2).
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! Current restrictions are:
!!  a - Spin-polarized case not tested.
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_alda(dtset,ixc,kxcg,mpi_enreg,nfft,ngfft,nspden,option,rhor,rhocut,rprimd)

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ixc,nfft,nspden,option
 real(dp),intent(in) :: rhocut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,2*nspden-1),rprimd(3,3)
 real(dp),intent(out) :: kxcg(2,nfft,*)

!Local variables -------------------------------------------------------
!No improved xc quadrature.
!No core correction.
!Dummy here.
!For debugging purposes (see tests below):
!integer :: i1,i2,i3,k1,n1,n2,n3
!real(dp) :: kx,rho,rhomax,ftest
!scalars
 integer :: ifft,ikxc,isp,n3xccc,ncut,nk3xc,nkxc,optionrhoxc,tim_fourdp
 logical :: non_magnetic_xc
 real(dp),parameter :: gsqcut=1._dp
 real(dp) :: enxc,rhocuttot,rhomin,vxcavg
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 real(dp) :: strsxc(6)
 real(dp) :: dum(0)
 real(dp),parameter   :: dummyvgeo(3)=zero
 real(dp),allocatable :: kxcr(:,:),rhog(:,:),rhorcut(:,:),vhartree(:)
 real(dp),allocatable :: vxc(:,:),xccc3d(:)

!***********************************************************************
!For debugging purposes (see tests below):

!ftest(i1,n1,k1) = 0._dp+1._dp*cos(k1*two_pi*float(i1)/float(n1))

!***********************************************************************

!Check input parameters.

 if (nspden > 2) then
   message = ' kxc_alda does not work yet for nspden > 2.'
   MSG_ERROR(message)
 end if

!Allocate memory.

 ABI_MALLOC(rhorcut,(nfft,nspden))
 ABI_MALLOC(rhog,(2,nfft))
 ABI_MALLOC(vhartree,(nfft))
 ABI_MALLOC(vxc,(nfft,nspden))

 call xcdata_init(xcdata,dtset=dtset,intxc=0,ixc=ixc,nspden=nspden)

 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

 ! Reinitialize the libxc module with the overriden values
 if (dtset%ixc<0) then
   call libxc_functionals_end()
 end if
 if (ixc<0) then
   call libxc_functionals_init(ixc,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)
 end if

!to be adjusted for the call to rhotoxc
 nk3xc=1

!Cut-off the density.

 rhorcut(:,:) = rhor(:,:)

 do isp = 1,nspden

   rhomin = maxval(rhorcut(:,isp))*rhocut

   ncut = 0
   rhocuttot = 0._dp

   do ifft = 1,nfft
     if (rhorcut(ifft,isp) < rhomin) then
       ncut = ncut+1
       rhocuttot = rhocuttot+rhorcut(ifft,isp)
       rhorcut(ifft,isp) = rhomin
     end if
   end do

   if (ncut > 0) then
     write (message,'(a,es10.3,3a,i1,a,i6,a,f6.3,3a,f6.3,a)') &
&     'rhocut = ',rhocut,'.',ch10,&
&     'For isp = ',isp,' the density was cut-off at ',ncut,' (',100._dp*float(ncut)/float(ifft),'%) grid points.',ch10,&
&     'These points account for ',100._dp*rhocuttot/sum(rhor(:,isp)),'% of the total density.'
     MSG_WARNING(message)
   end if

 end do

!Calculate the AL(S)DA kernel in real space.

 rhog(:,:) = 0._dp !We do not need the Hartree potential.
 tim_fourdp=0

 if ((option == 1).or.((option == 2).and.(nspden == 2))) then

   nkxc = 2*nspden-1
   n3xccc=0
   ABI_MALLOC(kxcr,(nfft,nkxc))
   ABI_MALLOC(xccc3d,(n3xccc))

   optionrhoxc = 2 !See rhotoxc.f

   call hartre(1,gsqcut,3,0,mpi_enreg,nfft,ngfft,1,zero,rhog,rprimd,dummyvgeo,vhartree)
   call rhotoxc(enxc,kxcr,mpi_enreg,nfft,ngfft,dum,0,dum,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&   optionrhoxc,rhorcut,rprimd,strsxc,1,vxc,vxcavg,xccc3d,xcdata,vhartr=vhartree)

!  DEBUG
!  fx for tests.
!  write (std_out,'(a)') ' kxc_alda: Using exchange-only kernel for tests.'
!  rhomin = minval(rhor(:,1))
!  rhomax = maxval(rhor(:,1))
!  write (std_out,'(a,es12.5,a,es12.5)') ' kxc_alda: rhomin = ',rhomin,' rhomax = ',rhomax
!  write (std_out,'(a)') ' kxc_alda: loping below 0.2*rhomax.'
!  kx = (3._dp/4._dp)*((3._dp/pi)**(1._dp/3._dp))
!  do ifft = 1,nfft
!  rho = rhor(ifft,1)
!  rho = max(rho,0.2_dp*rhomax)
!  kxcr(ifft,1) = -(4._dp/9._dp)*kx*(rho**(-2._dp/3._dp))
!  write (std_out,'(i4,a,es12.5)') ifft,': ',kxcr(ifft,1)
!  end do
!  write (std_out,'(a,es12.5)') 'kxcrmin: ',minval(kxcr(:,1))
!  write (std_out,'(a,es12.5)') 'kxcrmax: ',maxval(kxcr(:,1))
!  ENDDEBUG

!  DEBUG
!  test kernel.
!  write(std_out,'(a)') ' kxc_alda: Using test kernel for tests.'
!  n1 = ngfft(1) ; n2 = ngfft(2) ; n3 = ngfft(3)
!  do i3 = 0,n3-1
!  do i2 = 0,n2-1
!  do i1 = 0,n1-1
!  ifft = i1+n1*(i2+n2*i3)+1
!  kxcr(ifft,1) = ftest(i1,n1,1)*ftest(i2,n2,2)*ftest(i3,n3,3)
!  end do
!  end do
!  end do
!  ENDDEBUG

!  Calculate the Fourier transform of the AL(S)DA kernel.

   if (option == 1) then
     do ikxc = 1,nkxc
       call fourdp(1,kxcg(:,:,ikxc),kxcr(:,ikxc),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     end do
   else
     call fourdp(1,kxcg(:,:,1),kxcr(:,2),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
   end if

 else if ((option == 2).and.(nspden == 1)) then

   nkxc = 2
   n3xccc=0
   ABI_MALLOC(kxcr,(nfft,nkxc))
   ABI_MALLOC(xccc3d,(n3xccc))

   optionrhoxc = -2 !See rhotoxc.f

   call hartre(1,gsqcut,3,0,mpi_enreg,nfft,ngfft,1,zero,rhog,rprimd,dummyvgeo,vhartree)
   call rhotoxc(enxc,kxcr,mpi_enreg,nfft,ngfft,dum,0,dum,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,&
&   optionrhoxc,rhorcut,rprimd,strsxc,1,vxc,vxcavg,xccc3d,xcdata,vhartr=vhartree)

   kxcr(:,2) = 0.5_dp*kxcr(:,2)

!  Calculate the Fourier transform of the up-down channel of the AL(S)DA kernel.

   call fourdp(1,kxcg(:,:,1),kxcr(:,2),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)

 else
   write (message,'(4a,i0)')'  Invalid option = ',option
   MSG_ERROR(message)
 end if

!DEBUG
!write(std_out,*) ' kxc_alda:  Exc  = ',enxc
!write(std_out,*) ' kxc_alda: <Vxc> = ',vxcavg
!ENDDEBUG

! Revert libxc module to the original settings
 if (ixc<0) then
   call libxc_functionals_end()
 end if
 if (dtset%ixc<0) then
   call libxc_functionals_init(dtset%ixc,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)
 end if

!Free memory.
 ABI_FREE(rhorcut)
 ABI_FREE(rhog)
 ABI_FREE(vhartree)
 ABI_FREE(vxc)
 ABI_FREE(kxcr)
 ABI_FREE(xccc3d)

end subroutine kxc_alda
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_pgg
!! NAME
!! kxc_pgg
!!
!! FUNCTION
!! Compute the PGG-exchange kernel in reciprocal space
!! (Phys. Rev. Lett. 76, 1212 (1996) [[cite:Petersilka1996]]).
!!
!! INPUTS
!!  gmet=reciprocal space metrix (bohr**-2)
!!  npw=number of plane waves
!!  rcut_coulomb=real space cutoff radius for Coulomb interaction (bohr)
!!  susmat(2,npw,npw)=density weighted squared density matrix in reciprocal space
!!  ucvol=unit cell volume (bohr**3)
!!
!! OUTPUT
!!  khxcg(2,npwdiel,nspden,npwdiel,nspden)=PGG-exhange kernel in G space, at
!!       full interaction strength
!!
!! NOTES
!! The density weighted squared density matrix (actually the reduced=summed-over-spin
!! density matrix) is convolved with the spherically cutoff Coulomb interaction.
!!
!! WARNINGS
!! a - 'rcut_coulomb' should be chosen consistently with cutoffs elsewhere,
!!     for instance dieltcel8.f
!! b - applicable for spin-polarized case as well, through input 'susmat',
!!     but this has not been checked
!!
!! TODO
!! If simply the squared density matrix is input through 'susmat' the
!! exchange energy is recovered as the zero-G component of the resulting 'khxcg'
!! (then not the kernel of course). This could help to check convergence
!! with respect to 'npw'. See +ex_pgg comment.
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_pgg(gmet,kg,khxcg,npw,rcut_coulomb,susmat,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 real(dp),intent(in) :: rcut_coulomb,ucvol
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),susmat(2,npw,npw)
 real(dp),intent(out) :: khxcg(2,npw,npw)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ig,ii,ikg11,ikg12,ikg13,ikg21,ikg22,ikg23,ipw1,ipw2
 integer :: j1,j2,j3,jg,jj
 real(dp),parameter :: diffgsq=1.d-2
 real(dp) :: gred1,gred2,gred3,gsquar,tpisq
!arrays
 integer :: kgmax(3)
 integer,allocatable :: index_g_inv(:,:,:),jgarr(:)
 real(dp),allocatable :: gsq(:),sumg(:),vcoul(:)

! *************************************************************************

!DEBUG
!write(std_out,*) '%kxc_pgg: enter'
!write(std_out,*) 'npw=',npw
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 kgmax(:)=0
 do ipw1=1,npw
   do jj=1,3
     kgmax(jj)=max( kg(jj,ipw1), kgmax(jj) )
   end do
 end do

!DEBUG
!write(std_out,*) 'kgmax:',kgmax(1:3)
!ENDDEBUG

!Perform allocations
 ABI_MALLOC(index_g_inv,(-2*kgmax(1):2*kgmax(1),-2*kgmax(2):2*kgmax(2),-2*kgmax(3):2*kgmax(3)))
 ABI_MALLOC(jgarr,(npw))
 ABI_MALLOC(gsq,(npw))
 ABI_MALLOC(sumg,(2))
 ABI_MALLOC(vcoul,(npw))

!DEBUG
!write(std_out,*) '%kxc_pg: creating plane wave index and coulomb potential'
!ENDDEBUG
 index_g_inv(:,:,:)=0
 do ipw1=1,npw
   index_g_inv(kg(1,ipw1),kg(2,ipw1),kg(3,ipw1))=ipw1

!  DEBUG
!  write(std_out,'(i5,2x,3i3,2x,i4)') ipw1,kg(1,ipw1),kg(2,ipw1),kg(3,ipw1)
!  ENDDEBUG

   gred1=dble(kg(1,ipw1))
   gred2=dble(kg(2,ipw1))
   gred3=dble(kg(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +2.0_dp*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3) )
!  Distinguish G=0 from other elements
   if(gsquar > 1.0d-12)then
     vcoul(ipw1)=four_pi/gsquar*(1._dp-cos(sqrt(gsquar)*rcut_coulomb))
   else
     vcoul(ipw1)=four_pi*0.5_dp*rcut_coulomb**2
   end if

 end do

!DEBUG
!write(std_out,*) '%kxc_pg: starting convolution integral'
!ENDDEBUG

!loop over G1,G2 components of the density matrix
 do ipw2=1,npw
   ikg21=kg(1,ipw2)
   ikg22=kg(2,ipw2)
   ikg23=kg(3,ipw2)

   do ii=1,npw
     j1=ikg21-kg(1,ii)
     j2=ikg22-kg(2,ii)
     j3=ikg23-kg(3,ii)
     jgarr(ii)=index_g_inv(j1,j2,j3)
   end do

   do ipw1=1,ipw2
     ikg11=kg(1,ipw1)
     ikg12=kg(2,ipw1)
     ikg13=kg(3,ipw1)

!    do the convolution integral
     sumg(:)=0._dp
     do ii=1,npw

       if( jgarr(ii) /= 0 ) then

         i1=ikg11-kg(1,ii)
         i2=ikg12-kg(2,ii)
         i3=ikg13-kg(3,ii)

!        j1=ikg21-kg(1,ii)
!        j2=ikg22-kg(2,ii)
!        j3=ikg23-kg(3,ii)

         ig=index_g_inv(i1,i2,i3)
!        jg=index_g_inv(j1,j2,j3)

         if( ig /= 0 ) then

           jg=jgarr(ii)

!          DEBUG
!          write(std_out,'(i5,2x,3i3,1x,3i3,2x,2i4)') ii,i1,i2,i3,&
!          &                                             kg(1,jg),kg(2,jg),kg(3,jg),&
!          &                                             ig,jg
!          ENDDEBUG

           sumg(1)=sumg(1)+susmat(1,ig,jg)*vcoul(ii)
           sumg(2)=sumg(2)+susmat(2,ig,jg)*vcoul(ii)

         end if

       end if

     end do
     khxcg(:,ipw1,ipw2)=-sumg(:)*ucvol

!    DEBUG
!    if(ipw1==ipw2) write(std_out,'(2i4,2(1x,es14.6))') ipw1,ipw2,khxcg(1,ipw1,ipw1),vcoul(ipw1)
!    write(std_out,'(2i4,3(1x,es14.6))') ipw1,ipw2,khxcg(1:2,ipw1,ipw2),vcoul(ipw1)
!    ENDDEBUG

   end do
 end do

!DEBUG
!verify hermiticity, note: ipw1 loop above must end at npw
!write(std_out,*) '%kxc_pgg: check hermiticity of pgg kernel'
!do ipw2=1,npw,max(2,npw/10)
!do ipw1=ipw2,npw,max(2,npw/10)
!write(std_out,'(2i4,2(1x,es14.6))') ipw1,ipw2,&
!&   khxcg(1,ipw1,ipw2)-khxcg(1,ipw2,ipw1),&
!&   khxcg(2,ipw1,ipw2)+khxcg(2,ipw2,ipw1)
!end do
!end do
!ENDDEBUG

!Impose hermiticity
 write(std_out,*) '%kxc_pg: imposing hermiticity'
 do ipw2=1,npw
   do ipw1=ipw2+1,npw
     khxcg(1,ipw1,ipw2)= khxcg(1,ipw2,ipw1)
     khxcg(2,ipw1,ipw2)=-khxcg(2,ipw2,ipw1)
   end do
 end do

!DEBUG
!write(std_out,'(a10,2(1x,es20.12))') '+ex_pgg? ', 0.5_dp*khxcg(1,1,1)/ucvol
!ENDDEBUG

 ABI_FREE(index_g_inv)
 ABI_FREE(jgarr)
 ABI_FREE(gsq)
 ABI_FREE(sumg)
 ABI_FREE(vcoul)

!DEBUG
!write(std_out,*) '%kxc_pgg: done'
!ENDDEBUG

end subroutine kxc_pgg
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_eok
!! NAME
!! kxc_eok
!!
!! FUNCTION
!!  Compute the linear (ixceok = 1) or non-linear (ixceok = 2)
!!  energy optimized kernel of Dobson and Wang, in reciprocal
!!  space, on the FFT grid.
!!  See J. Dobson and J. Wang, Phys. Rev. B 62, 10038 (2000) [[cite:Dobson2000]].
!!
!! INPUTS
!!  ixceok = 1 linear energy optimized kernel.
!!         = 2 non-linear energy optimized kernel.
!!  mpi_enreg=information about MPI parallelization
!!  nfft = number of fft grid points.
!!  ngfft(1:3) = integer fft box dimensions, see getng for ngfft(4:8).
!!  nspden = number of spin-density components.
!!  rhor(nfft,nspden) = electron density in real space in electrons/bohr**3
!!   (total in first half and spin-up in second half if nspden = 2).
!!  rhocut = cut-off density for the local kernels (ALDA, EOK),
!!           relative to max(rhor(:,:)).
!! OUTPUT
!!  kxcg(2,nfft,2*nspden-1) = the EOK kernel in reciprocal space, on the FFT grid.
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!!
!! PARENTS
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_eok(ixceok,kxcg,mpi_enreg,nfft,ngfft,nspden,rhor,rhocut)

!Arguments -------------------------------------------------------------
!scalars
 integer,intent(in) :: ixceok,nfft,nspden
 real(dp),intent(in) :: rhocut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,2*nspden-1)
 real(dp),intent(out) :: kxcg(2,nfft,2*nspden-1)

!Local variables -------------------------------------------------------
!Maximum value allowed for rs.
!scalars
 integer :: ifft,ikxc,ncut,nkxc,nlop,tim_fourdp
 real(dp),parameter :: rslim=50._dp,dummyvgeo(3)=zero
 real(dp) :: a2,a3,a4,rho,rhocuttot,rhomin,rs
 character(len=500) :: message
!arrays
 real(dp),allocatable :: kxcr(:,:)

!***********************************************************************

!Check input parameters.

 if (nspden > 1) then
   message = ' kxc_eok does not work yet for nspden > 1.'
   MSG_ERROR(message)
 end if

!Values of a2, a3 and a4 for case 1
 a2 = 0.0_dp
 a3 = 0.0_dp
 a4 = 0.0_dp

 select case (ixceok)
   case (1)
     a2 = -0.51887_dp
     a3 =  4.9359d-03
     a4 = -5.9603d-05
   case (2)
     a2 = -0.50044_dp
     a3 =  4.9653d-03
     a4 = -3.3660d-05
   case default
     message =  ' kxc_eok: ixceok /= 1 (linear EOK) or 2 (non-linear EOK).'
     MSG_ERROR(message)
 end select

!Allocate memory.

 nkxc = 2*nspden-1

 ABI_MALLOC(kxcr,(nfft,nkxc))

!Calculate the energy optimized kernel in real space.

 nlop = 0

 rhomin = rhocut*maxval(rhor(:,:))

 ncut = 0
 rhocuttot = 0._dp

 do ifft = 1,nfft

   rho = rhor(ifft,1)

   if (rho < rhomin) then
     ncut = ncut+1
     rhocuttot = rhocuttot+rho
     rho = rhomin
   end if

   rs = (3._dp/(4._dp*pi*rho))**(1._dp/3._dp)

   if (rs > rslim) then
     rs = rslim
     nlop = nlop+1
   end if

   kxcr(ifft,1) = a2*rs**2+a3*rs**3+a4*rs**4

 end do

 if (ncut > 0) then
   write (message,'(a,es10.3,3a,i1,a,i6,a,f6.3,3a,f6.3,a)') &
&   'rhocut = ',rhocut,'.',ch10,&
&   'For isp = ',1,' the density was cut-off at ',ncut,' (',100._dp*float(ncut)/float(ifft),'%) grid points.',ch10,&
&   'These points account for ',100._dp*rhocuttot/sum(rhor(:,1)),'% of the total density.'
   MSG_WARNING(message)
 end if

 if (nlop > 0) then
   write (message,'(a,f6.2,a,i6,a,f6.3,a)') &
&   'rs still exceeds ',rslim,' Bohr at ',nlop,' (',100._dp*float(nlop)/float(ifft),'%) grid points (after cut-off).'
   MSG_WARNING(message)
 end if

!Calculate the Fourier transform of the energy optimized kernel.
 tim_fourdp=0
 do ikxc = 1,nkxc
   call fourdp(1,kxcg(:,:,ikxc),kxcr(:,ikxc),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
 end do

!Free memory.

 ABI_FREE(kxcr)

end subroutine kxc_eok
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_driver
!! NAME
!! kxc_driver
!!
!! FUNCTION
!! Calculate the exchange-correlation kernel in reciprocal space.
!! Require density in real space on the *full* FFT mesh
!!
!! INPUTS
!! Dtset<dataset_type>=all input variables in this dataset
!! Cryst<crystal_t>=Info on the crystal structure.
!! ixc = choice for the exchange-correlation potential.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/variables/vargs.htm#ngfft
!! nfft_tot = Total number of points on the FFT grid.
!! nspden=Number of independent spin densities.
!! rhor(nfft_tot,nspden) = the charge density on the full FFT grid.
!!  (total in first half and spin-up in second half if nspden=2)
!! npw: the size of kernel matrix
!! dim_kxcg=dimension of the kernel.
!! comm=MPI communicator.
!! [dbg_mode]=Optional flag used to switch on the debug mode.
!!
!! OUTPUT
!!  FIXME: Why are we using nfft_tot instead of the G-sphere
!!  kxcg(nfft_tot,dim_kxcg) = the exchange-correlation potential on the FFT grid.
!!  warning: the kernel is not divided by the unit cell volume
!!
!! NOTES
!!  No xc quadrature
!!  No nl core correction
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_driver(Dtset,Cryst,ixc,ngfft,nfft_tot,nspden,rhor,npw,dim_kxcg,kxcg,gvec,comm,dbg_mode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,npw,nfft_tot,nspden,dim_kxcg,comm
 logical,optional,intent(in) :: dbg_mode
 type(crystal_t),intent(in) :: Cryst
 type(dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: gvec(3,npw),ngfft(18)
 real(dp),intent(in) :: rhor(nfft_tot,nspden)
 complex(gwpc),intent(out) :: kxcg(nfft_tot,dim_kxcg)

!Local variables ------------------------------
!scalars
 integer :: cplex,i1,i2,i3,ig,igp,iq,ir,n3xccc,ngfft1,ngfft2,izero
 integer :: ngfft3,nkxc,option,ikxc,nk3xc,my_rank,master,unt_dmp
 logical :: non_magnetic_xc
 real(dp) :: enxc,expo,gpqx,gpqy,gpqz,gsqcut,vxcavg
 character(len=500) :: msg,fname
 type(xcdata_type) :: xcdata
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: qphon(3),strsxc(6),dum(0)
 real(dp),parameter   :: dummyvgeo(3)=zero
 real(dp),allocatable :: kxcpw_g(:,:),kxcr(:,:),phas(:,:,:)
 real(dp),allocatable :: rhog(:,:),vhartr(:),kxcpw_r(:,:),vxclda(:,:)
 real(dp),allocatable :: xccc3d(:),xx(:,:)
 real(dp),allocatable :: my_kxcg(:,:)

!************************************************************************

 ABI_CHECK(Dtset%nsppol==1,'nsppol/=1 not coded')
 ABI_CHECK(Dtset%nspden==1,'nsppol/=1 not coded')
 ABI_CHECK(nfft_tot==PRODUCT(ngfft(1:3)),"mismatch in nfftot")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')
 my_rank = xmpi_comm_rank(comm)
 master=0

 write(msg,'(a,i3)') ' kxc_driver: calculating exchange-correlation kernel using ixc = ',ixc
 call wrtout(std_out,msg,'COLL')

 call xcdata_init(xcdata,dtset=Dtset,intxc=0,ixc=ixc,nspden=nspden)

 if (ALL(xcdata%xclevel/=(/1,2/))) then
   write(msg,'(a,i0)')"Unsupported xclevel = ",xcdata%xclevel
   MSG_ERROR(msg)
 end if

 ngfft1=ngfft(1)
 ngfft2=ngfft(2)
 ngfft3=ngfft(3)

 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

 if (ixc>=1.and.ixc<11) then ! LDA case
   nkxc= 2*min(nspden,2)-1   ! 1 or 3
 else                        ! GGA case
   nkxc=12*min(nspden,2)-5   ! 7 or 19
   ABI_CHECK(dtset%xclevel==2,"Functional should be GGA")
   MSG_ERROR("GGA functional not tested")
 end if

 ABI_MALLOC(kxcr,(nfft_tot,nkxc))

!gsqcut and rhog are zeroed because they are not used by rhotoxc if 1<=ixc<=16 and option=0
 gsqcut=zero

 ABI_MALLOC(rhog,(2,nfft_tot))
 ABI_MALLOC(vhartr,(nfft_tot))
 rhog(:,:)=zero
!MG FIXME this is the 3D core electron density for XC core correction (bohr^-3)
!should implement the non linear core correction
 n3xccc=0
 ABI_MALLOC(xccc3d,(n3xccc))
 ABI_MALLOC(vxclda,(nfft_tot,nspden))

 option=2 ! 2 for Hxc and kxcr (no paramagnetic part if nspden=1)
 qphon =zero

!to be adjusted for the call to rhotoxc
 nk3xc=1
 izero=0

 ! Reinitialize the libxc module with the overriden values
 if (dtset%ixc<0) then
   call libxc_functionals_end()
 end if
 if (ixc<0) then
   call libxc_functionals_init(ixc,Dtset%nspden,xc_tb09_c=Dtset%xc_tb09_c)
 end if

 call hartre(1,gsqcut,3,izero,MPI_enreg_seq,nfft_tot,ngfft,1,zero,rhog,Cryst%rprimd,dummyvgeo,vhartr)

!Compute the kernel.
 call rhotoxc(enxc,kxcr,MPI_enreg_seq,nfft_tot,ngfft,&
& dum,0,dum,0,nkxc,nk3xc,non_magnetic_xc,&
& n3xccc,option,rhor,Cryst%rprimd,&
& strsxc,1,vxclda,vxcavg,xccc3d,xcdata,vhartr=vhartr)

 ABI_FREE(rhog)
 ABI_FREE(vhartr)
!DEBUG print Kxc
 if (present(dbg_mode)) then
   if (dbg_mode .and. my_rank==master) then
     fname = 'xc_Kxc.xsf'
     if (open_file(fname,msg,newunit=unt_dmp,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
     call printxsf(ngfft1,ngfft2,ngfft3,kxcr(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 ABI_FREE(xccc3d)
 ABI_FREE(vxclda)

 ABI_MALLOC(my_kxcg,(2,nfft_tot))

 do ikxc=1,nkxc
   call fourdp(1,my_kxcg,kxcr(:,ikxc),-1,MPI_enreg_seq,nfft_tot,1,ngfft,0)
   kxcg(:,ikxc)=CMPLX(my_kxcg(1,:),my_kxcg(2,:))
 end do

!write(std_out,*)"kxcr(r=0)",kxcr(1,1)
!write(std_out,*)"my_kxg(G=0)",my_kxcg(:,1)
!write(std_out,*)"SUM kxcr/nfft_tot ",SUM(kxcr(:,1))/nfft_tot
!write(std_out,*)"SUM my_kxg ",SUM(kxcg(:,1))

 ABI_FREE(my_kxcg)

!MG this part is never executed, but one should use dfpt_mkvxc for the GGA kernel.
 if (xcdata%xclevel==2) then
   MSG_ERROR("check GGA implementation")
   cplex=2
   ABI_MALLOC(phas,(cplex*nfft_tot,npw,nspden))
   ABI_MALLOC(kxcpw_r,(cplex*nfft_tot,nspden))
   ABI_MALLOC(xx,(3,nfft_tot))
   ABI_MALLOC(kxcpw_g,(2,nfft_tot))

   kxcg = czero

!  find the coordinates for all r in the FFT grid
   ir=0
   do i3=1,ngfft3
     do i2=1,ngfft2
       do i1=1,ngfft1
         ir=ir+1
         xx(1,ir)=dble((i1-1))/ngfft1
         xx(2,ir)=dble((i2-1))/ngfft2
         xx(3,ir)=dble((i3-1))/ngfft3
       end do
     end do
   end do

   do iq=1,1

!    Calculate at once exp(i(G+q).r), for all possible q,G,r
     do ig=1,npw
       gpqx=dble(gvec(1,ig))
       gpqy=dble(gvec(2,ig))
       gpqz=dble(gvec(3,ig))
       do ir=1,nfft_tot
         expo=gpqx*xx(1,ir)+gpqy*xx(2,ir)+gpqz*xx(3,ir)
         phas(2*ir-1,ig,1)= cos(two_pi*expo)
         phas(2*ir,ig,1) =  sin(two_pi*expo)
       end do
     end do

!    Calculate $K(G,G'',q)=\frac{1}{nfft_tot}\sum_{r} exp(-i(q+G_{2}).r_{2} kxcr(r_{1}r_{2}) exp(i(q+G_{1}).r_{1} $

     do igp=1,npw

       kxcpw_r(:,:)=zero

       call dfpt_mkvxc(cplex,ixc,kxcr,MPI_enreg_seq,nfft_tot,ngfft,dum,0,dum,0,nkxc,non_magnetic_xc,&
&       nspden,n3xccc,option,qphon(:),phas(:,igp,:),Cryst%rprimd,1,kxcpw_r,xccc3d)

!      FFT the first index to --> to G space
       call fourdp(cplex,kxcpw_g(:,:),kxcpw_r(:,1),-1,MPI_enreg_seq,nfft_tot,1,ngfft,0)

!      kxcg(:,igp,iq)=CMPLX(kxcpw_g(1,igfft(:)),kxcpw_g(2,igfft(:)))
!      kxcg(:,igp)=CMPLX(kxcpw_g(1,igfft(:)),kxcpw_g(2,igfft(:)))

     end do ! igp
   end do ! iq

   ABI_FREE(phas)
   ABI_FREE(kxcpw_r)
   ABI_FREE(xx)
   ABI_FREE(kxcpw_g)
 end if !xclevel==2

! Revert libxc module to the original settings
 if (ixc<0) then
   call libxc_functionals_end()
 end if
 if (dtset%ixc<0) then
   call libxc_functionals_init(dtset%ixc,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)
 end if

 call destroy_mpi_enreg(MPI_enreg_seq)
 ABI_FREE(kxcr)

end subroutine kxc_driver
!!***

!----------------------------------------------------------------------

!!****f* m_kxc/kxc_ADA
!! NAME
!! kxc_ADA
!!
!! FUNCTION
!! Calculate exchange-correlation kernel in reciprocal space
!!
!! INPUTS
!! Dtset <type(dataset_type)>=all input variables in this dataset
!! Cryst<crystal_t>=Info on the unit cell.
!! ixc = choice for the exchange-correlation potential.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/variables/vargs.htm#ngfft
!! nfft = total number of points on the FFT grid.
!! rhor(nfft,nspden) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! npw: the size of kernel matrix
!! dim_kxcg=dimension of the kernel.
!! comm=MPI communicator.
!! [dbg_mode]=Set it to .TRUE. to switch on the debug mode.
!!
!! OUTPUT
!!  kxcg(nfft,dim_kxcg) = the exchange-correlation potential on the FFT grid.
!!  warning: the kernel is not divided by unit cell volume
!!
!! NOTES
!!  No xc quadrature
!!  No nl core correction
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,fourdp_6d,hartre,initmpi_seq
!!      libxc_functionals_end,libxc_functionals_init,printxsf,rhotoxc,wrtout
!!      xcdata_init
!!
!! SOURCE

subroutine kxc_ADA(Dtset,Cryst,ixc,ngfft,nfft,nspden,rhor,&
&                  npw,nqibz,qibz,fxc_ADA,gvec,comm,kappa_init,dbg_mode)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,nfft,nspden,npw,comm
 real(dp),intent(in),optional :: kappa_init
 logical,optional,intent(in) :: dbg_mode
 type(crystal_t),intent(in) :: Cryst
 type(dataset_type),intent(in) :: Dtset
!arrays
 integer,intent(in) :: gvec(3,npw),ngfft(18)
 integer,intent(in) :: nqibz
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: qibz(3,nqibz)
 complex(gwpc),intent(out) :: fxc_ADA(npw,npw,nqibz)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,ig,igp,ir,irp,n3xccc,ngfft1,ngfft2,izero !,isp
 integer :: ngfft3,nkxc,option,ikxc,ierr,nproc
 integer :: nk3xc,igrid,iqbz,my_rank,master,unt_dmp,gmgp_idx
 logical :: non_magnetic_xc
 real(dp) :: enxc,gsqcut,ucvol !,rs,Kx,Kc
 real(dp) :: vxcavg,kappa,abs_qpg_sq,abs_qpgp_sq
 real(dp) :: difx,dify,difz,inv_kappa_sq
 character(len=500) :: msg,fname
 type(MPI_type) :: MPI_enreg_seq
 type(xcdata_type) :: xcdata
!arrays
 real(dp) :: qpg(3),qpgp(3),qphon(3),strsxc(6),q_point(3),dum(0)
 real(dp),parameter   :: dummyvgeo(3)=zero
 real(dp),allocatable :: kxcr(:,:)
 real(dp),allocatable :: rhog(:,:),vhartr(:),vxclda(:,:)
 real(dp),allocatable :: xccc3d(:),my_rhor(:,:)
 real(dp),allocatable :: my_kxcg(:,:)
 real(dp),allocatable :: rhotilder(:,:)
 complex(gwpc),allocatable :: my_fxc_ADA_ggpq(:,:,:)
 complex(gwpc),allocatable :: FT_fxc_ADA_ggpq(:,:,:),dummy(:,:)
 real(dp),allocatable :: rvec(:,:),my_fxc_ADA_rrp(:,:)
 real(dp) :: rmrp(3),abs_rmrp
 integer :: n1,n2,n3,ig_idx_fft(npw)

! ************************************************************************

 ABI_CHECK(Dtset%nsppol==1,'nsppol/=1 not coded')
 ABI_CHECK(nspden==1,'nsppol/=1 not coded')
 ABI_CHECK(nfft==PRODUCT(ngfft(1:3)),"mismatch in nfftot")

!Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)

 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm)
 master=0

 write(msg,'(a,i3)') ' kxc_ADA: calculating exchange-correlation kernel using ixc = ',ixc
 call wrtout(std_out,msg,'COLL')
 call wrtout(std_out,' kxc_ADA: using smeared density','COLL')

 if (.not.present(kappa_init)) then
   kappa = 2.1_dp
 else
   kappa = kappa_init
 end if
 write(msg,'(a,F10.3)') ' kxc_ADA: inverse smearing length, kappa = ',kappa
 call wrtout(std_out,msg,'COLL')
 inv_kappa_sq = one/(kappa*kappa)

 call xcdata_init(xcdata,dtset=dtset,intxc=0,ixc=ixc,nspden=nspden)

 if (ALL(xcdata%xclevel/=(/1,2/))) then
   write(msg,'(a,i0)')"Unsupported xclevel = ",xcdata%xclevel
   MSG_ERROR(msg)
 end if

 ngfft1=ngfft(1)
 ngfft2=ngfft(2)
 ngfft3=ngfft(3)

 non_magnetic_xc=(abs(dtset%usepawu)==4.or.dtset%usepawu==14)

 if (ixc>=1.and.ixc<11) then      ! LDA case
   nkxc= 2*min(Dtset%nspden,2)-1  ! 1 or 3
 else                             ! GGA case
   nkxc=12*min(Dtset%nspden,2)-5  ! 7 or 19
   ABI_CHECK(dtset%xclevel==2,"Functional should be GGA")
   MSG_ERROR("GGA functional not implemented for ADA vertex")
 end if

 ABI_MALLOC(kxcr,(nfft,nkxc))

!gsqcut and rhog are zeroed because they are not used by rhotoxc if 1<=ixc<=16 and option=0
 gsqcut=zero

 ABI_MALLOC(rhog,(2,nfft))
 ABI_MALLOC(vhartr,(nfft))
 rhog(:,:)=zero
!MG FIXME this is the 3D core electron density for XC core correction (bohr^-3)
!should implement the non linear core correction
 n3xccc=0
 ABI_MALLOC(xccc3d,(n3xccc))
 ABI_MALLOC(vxclda,(nfft,nspden))

 option=2 ! 2 for Hxc and kxcr (no paramagnetic part if nspden=1)
 qphon(:)=zero

!to be adjusted for the call to rhotoxc
 nk3xc=1

!Compute the kernel.
 izero=0

!DEBUG print density
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     fname = 'xc_ADA_den.xsf'
     if (open_file(fname,msg,newunit=unt_dmp,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
     call printxsf(ngfft1,ngfft2,ngfft3,rhor(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

!Calculate the smeared density
 ABI_MALLOC(my_rhor,(nfft,nspden))
 ABI_MALLOC(rhotilder,(nfft,nspden))
 ucvol = Cryst%ucvol
 my_rhor = rhor
!do isp = 1,nsppol
!call calc_smeared_density(my_rhor(:,isp),1,rhotilder(:,isp),nfft,ngfft,npw,&
!&   gvec,Cryst%gprimd,Cryst%ucvol,MPI_enreg_seq,paral_kgb0,kappa_in=kappa)
!my_rhor(:,isp) = rhotilder(:,isp)
!end do

!DEBUG print smeared density
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     fname = 'xc_ADA_smeared_den.xsf'
     if (open_file(fname,msg,newunit=unt_dmp,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
     call printxsf(ngfft1,ngfft2,ngfft3,my_rhor(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 ! Reinitialize the libxc module with the overriden values
 if (dtset%ixc<0) then
   call libxc_functionals_end()
 end if
 if (ixc<0) then
   call libxc_functionals_init(ixc,Dtset%nspden,xc_tb09_c=Dtset%xc_tb09_c)
 end if

 call hartre(1,gsqcut,3,izero,MPI_enreg_seq,nfft,ngfft,1,zero,rhog,Cryst%rprimd,dummyvgeo,vhartr)
 call rhotoxc(enxc,kxcr,MPI_enreg_seq,nfft,ngfft,&
& dum,0,dum,0,nkxc,nk3xc,non_magnetic_xc,&
& n3xccc,option,my_rhor,Cryst%rprimd,&
& strsxc,1,vxclda,vxcavg,xccc3d,xcdata,vhartr=vhartr)

!Check for extreme (NaN) values
!do ir=1,nfft
!if (isnan(kxcr(ir,1))) kxcr(ir,1) = HUGE(kxcr(ir,1))
!end do

!DEBUG test with direct way of calculating Kxc
!do i1=1,nfft
!rs = (three/(four_pi*my_rhor(i1,1)))**third
!Kx = 16._dp/27._dp*0.3141592653589793e1_dp*(rs**2)*(-0.4581652_dp)
!
!Kc =  -0.4e1_dp / 0.9e1_dp * 0.3141592654e1_dp * rs ** 4 &
!* (0.207271333333333333333333333333e-1_dp * &
!(-0.177442658629204480000000e3_dp * rs - 0.17565190511219200000000e2_dp &
!* sqrt(rs) - 0.1332650665120000e2_dp * rs ** 2 &
!- 0.51031691247948928000000e2_dp * rs ** (0.3e1_dp / 0.2e1_dp)) &
!* rs ** (-0.3e1_dp / 0.2e1_dp) / (rs + 0.37274400e1_dp * sqrt(rs) &
!+ 0.129352000e2_dp) ** 2 / (-sqrt(rs) - 0.1049800_dp) &
!+ 0.518178333333333333333333333333e-2_dp * rs ** (-0.3e1_dp / 0.2e1_dp) &
!* (0.617071835390850041282140897280e3_dp * sqrt(rs) &
!+ 0.659369347307557491857191871552e5_dp * rs ** 2 + &
!0.700403648491298930017835369562e5_dp * rs ** (0.3e1_dp / 0.2e1_dp) &
!+ 0.398437532951539263722720308167e5_dp * rs ** (0.5e1_dp / 0.2e1_dp) &
!+ 0.368852071032531998953472000000e4_dp * rs ** (0.7e1_dp / 0.2e1_dp) &
!+ 0.5330602660480000e2_dp * rs ** (0.9e1_dp / 0.2e1_dp) &
!+ 0.143783940386264738593799346176e5_dp * rs ** 3 &
!+ 0.124672564145568409213848436081e5_dp * rs &
!+ 0.557398029956167136000000e3_dp * rs ** 4) &
!/ (rs + 0.37274400e1_dp * sqrt(rs) + 0.129352000e2_dp) ** 4 &
!/ (-sqrt(rs) - 0.1049800_dp) ** 2)
!kxcr(i1,1) = Kx + Kc
!end do
!END DEBUG

!DEBUG print Kxc
 if (present(dbg_mode)) then
   if (dbg_mode.and.my_rank==master) then
     fname = 'xc_ADA_Kxc.xsf'
     if (open_file(fname,msg,newunit=unt_dmp,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
     call printxsf(ngfft1,ngfft2,ngfft3,kxcr(:,1),Cryst%rprimd,(/zero,zero,zero/),&
&     Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xcart,Cryst%znucl,unt_dmp,0)
     close(unt_dmp)
   end if
 end if
!DEBUG

 ABI_FREE(xccc3d)
 ABI_FREE(vxclda)
 ABI_FREE(vhartr)

 ABI_MALLOC(my_kxcg,(2,nfft))

 do ikxc=1,nkxc
   call fourdp(1,my_kxcg,kxcr(:,ikxc),-1,MPI_enreg_seq,nfft,1,ngfft,0)
!  kxcg(:,ikxc)=CMPLX(my_kxcg(1,:),my_kxcg(2,:))
 end do
!TODO Check symmetry of kxcg

!set up ADA vertex
 ABI_MALLOC(my_fxc_ADA_ggpq,(npw,npw,nqibz))
 my_fxc_ADA_ggpq = czero
!Calculate f_xc(R,R')=(kappa^2/2)K_xc[\tilde{n}](G-G')
!x(1/(kappa^2+|q+G|^2) + 1/(kappa^2+|q+G'|^2
!First get G vectors and indices

 ierr=0
 do iqbz=1,nqibz
   q_point(:) = qibz(:,iqbz)
!
   do ig=1,npw
     do igp=1,npw
!      Calculate |q+G| and |q+G'|
       qpg(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:)+gvec(:,ig))
       qpgp(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:)+gvec(:,igp))
       abs_qpg_sq = 1.0_dp/(1.0_dp+dot_product(qpg,qpg)*inv_kappa_sq)
       abs_qpgp_sq = 1.0_dp/(1.0_dp+dot_product(qpgp,qpgp)*inv_kappa_sq)

       gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfft)
       if (gmgp_idx>0) then
         my_fxc_ADA_ggpq(ig,igp,iqbz) = half*CMPLX(my_kxcg(1,gmgp_idx), my_kxcg(2,gmgp_idx))*(abs_qpg_sq+abs_qpgp_sq)
       else
         ierr=ierr+1
         my_fxc_ADA_ggpq(ig,igp,iqbz) = czero
       end if
     end do
   end do
   if (ierr/=0) then
     write(msg,'(a,i4,3a)')&
&     ' Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&     ' Enlarge the FFT mesh to get rid of this problem. '
     MSG_WARNING(msg)
   end if
 end do

 fxc_ADA = my_fxc_ADA_ggpq

!do iqbz=1,nqibz
!call hermitianize(my_fxc_ADA_ggpq(:,:,iqbz),"All")
!end do


!DEBUG check symmetry
 if (.FALSE.) then
!  do iqbz=1,nkptgw
!  do ig=1,npw
!  do igp=ig,npw
!  if (ABS(REAL(fxc_ADA(ig,igp,iqbz))-REAL(fxc_ADA(igp,ig,iqbz)))>tol15.OR.&
!  ABS(AIMAG(fxc_ADA(ig,igp,iqbz))-AIMAG(-fxc_ADA(igp,ig,iqbz)))>tol15) then
!  write(std_out,*) 'Elements:'
!  write(std_out,*) 'fxc_ADA(ig,igp,iqbz):',ig,igp,iqbz,fxc_ADA(ig,igp,iqbz)
!  write(std_out,*) 'fxc_ADA(igp,ig,iqbz):',igp,ig,iqbz,fxc_ADA(igp,ig,iqbz)
!  MSG_ERROR('fxc_ADA not symmetric')
!  end if
!  end do
!  end do
!  end do

!  write(std_out,*)"kxcr(r=0)",kxcr(1,1)
!  write(std_out,*)"my_kxg(G=0)",my_kxcg(:,1)
!  write(std_out,*)"SUM kxcr/nfft ",SUM(kxcr(:,1))/nfft
!  write(std_out,*)"SUM my_kxg ",SUM(kxcg(:,1))

!  DEBUG Check FT to real space
!  The real-space expression is:
!  f_xc(R,R')=(1/2)(kappa^2/(4*Pi))
!  \{K_xc[\tilde{n(R)}]+K_xc[\tilde{n(R')}]\}
!  x exp(-kappa|R-R'|)/|R-R'|
   ABI_MALLOC(my_fxc_ADA_rrp,(nfft,nfft))
   ABI_MALLOC(FT_fxc_ADA_ggpq,(npw,npw,nqibz))
   ABI_MALLOC(rvec,(3,nfft))
   ABI_MALLOC(dummy,(nfft,nfft))
   my_fxc_ADA_rrp=zero; FT_fxc_ADA_ggpq=czero; dummy=czero; rvec=zero

!  First find coordinates of real-space fft points
   igrid = 0
   ngfft1 = ngfft(1)
   ngfft2 = ngfft(2)
   ngfft3 = ngfft(3)
   do i3=0,ngfft3-1
     difz=dble(i3)/dble(ngfft3)
     do i2=0,ngfft2-1
       dify=dble(i2)/dble(ngfft2)
       do i1=0,ngfft1-1
         difx=dble(i1)/dble(ngfft1)
         igrid = igrid + 1
         rvec(1,igrid)=difx*Cryst%rprimd(1,1)+dify*Cryst%rprimd(1,2)+difz*Cryst%rprimd(1,3)
         rvec(2,igrid)=difx*Cryst%rprimd(2,1)+dify*Cryst%rprimd(2,2)+difz*Cryst%rprimd(2,3)
         rvec(3,igrid)=difx*Cryst%rprimd(3,1)+dify*Cryst%rprimd(3,2)+difz*Cryst%rprimd(3,3)
       end do
     end do
   end do
   if (igrid/=nfft) then
     MSG_ERROR('kxc_ADA: igrid not equal to nfft')
   end if

!  Construct kernel in real space
   do ir=1,nfft
     do irp=ir,nfft
       rmrp(:) = rvec(:,ir)-rvec(:,irp)
       abs_rmrp = sqrt(dot_product(rmrp,rmrp))
       my_fxc_ADA_rrp(ir,irp) = eighth*kappa*kappa*piinv* &
       (kxcr(ir,1)+kxcr(irp,1))* &
       EXP(-kappa*abs_rmrp)/(abs_rmrp+1.e-3_dp)
!      (a small convergence factor is introduced
!      to avoid a singularity)
       my_fxc_ADA_rrp(irp,ir) = my_fxc_ADA_rrp(ir,irp)
     end do
   end do

!  Find FFT index for all G
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
!  Use the following indexing (N means ngfft of the adequate direction)
!  0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
!  1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
   do ig=1,npw
     i1=modulo(gvec(1,ig),n1)
     i2=modulo(gvec(2,ig),n2)
     i3=modulo(gvec(3,ig),n3)
     ig_idx_fft(ig)=i1+1+n1*(i2+n2*i3)
   end do
!  FT kernel to reciprocal space for each q
   do iqbz=1,nqibz
     dummy = CMPLX(my_fxc_ADA_rrp,0.0_dp)
!    Multiply with q-point phase factors exp(-iq.r)*f_xc(r,r')*exp(iq.r')
     do ir=1,nfft
       do irp=1,nfft
!        Calculate q (variables defined for other purposes
!        are being re-used as dummy variables)
         q_point(:) = qibz(:,iqbz)
         qpg(:) = two_pi*MATMUL(Cryst%gprimd,q_point(:))
         abs_qpg_sq = dot_product(qpg(:),rvec(:,ir))
         abs_qpgp_sq = dot_product(qpg(:),rvec(:,irp))
         dummy(ir,irp) = EXP(-j_dpc*abs_qpg_sq)* &
&         dummy(ir,irp)* &
&         EXP(j_dpc*abs_qpgp_sq)
       end do
     end do
     call fourdp_6d(2,dummy,-1,MPI_enreg_seq,nfft,ngfft, 0)
     do ig=1,npw
       do igp=1,npw
         FT_fxc_ADA_ggpq(ig,igp,iqbz) = dummy(ig_idx_fft(ig),ig_idx_fft(igp))
       end do
     end do

!    Output
     msg=''
     if (iqbz<10) write(msg,'(a,i1,a)') './debug_fxc_ADA_q',iqbz,'.dat'
     if ((iqbz>9).and.(iqbz<100)) write(msg,'(a,i2,a)') './debug_fxc_ADA_q',iqbz,'.dat'
     if ((iqbz>99).and.(iqbz<1000)) write(msg,'(a,i3,a)') './debug_fxc_ADA_q',iqbz,'.dat'

     !open(777,file=TRIM(msg),STATUS='REPLACE')
     !do igp=1,npw
     !  do ig=1,npw
     !    write(777,*) ig,igp,REAL(my_fxc_ADA_ggpq(ig,igp,iqbz)),AIMAG(my_fxc_ADA_ggpq(ig,igp,iqbz)), &
     !    REAL(FT_fxc_ADA_ggpq(ig,igp,iqbz)),AIMAG(FT_fxc_ADA_ggpq(ig,igp,iqbz)), &
     !    ABS(ABS(my_fxc_ADA_ggpq(ig,igp,iqbz))-ABS(FT_fxc_ADA_ggpq(ig,igp,iqbz)))
     !  end do
     !  write(777,*) ''
     !end do
     !close(777)

   end do ! iqbz

   MSG_ERROR('Stopping in kxc_ADA for debugging')

   ABI_FREE(rvec)
   ABI_FREE(my_fxc_ADA_rrp)
   ABI_FREE(FT_fxc_ADA_ggpq)

   if (xcdata%xclevel==2) then
     MSG_ERROR(" GGA not implemented for kxc_ADA")
   end if !xclevel==2

 end if ! Debugging section

! Revert libxc module to the original settings
 if (ixc<0) then
   call libxc_functionals_end()
 end if
 if (dtset%ixc<0) then
   call libxc_functionals_init(dtset%ixc,dtset%nspden,xc_tb09_c=dtset%xc_tb09_c)
 end if

 call destroy_mpi_enreg(MPI_enreg_seq)
 ABI_FREE(my_kxcg)
 ABI_FREE(my_rhor)
 ABI_FREE(rhotilder)
 ABI_FREE(rhog)
 ABI_FREE(kxcr)

end subroutine kxc_ADA
!!***

!----------------------------------------------------------------------

end MODULE m_kxc
