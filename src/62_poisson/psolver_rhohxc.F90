!{\src2tex{textfont=tt}}
!!****f* ABINIT/psolver_rhohxc
!! NAME
!! psolver_rhohxc
!!
!! FUNCTION
!! Given rho(r), compute Hartree potential considering the system as
!! an isolated one. This potential is obtained from the convolution
!! of 1/r and rho(r), treated in Fourier space. This method is a wrapper around
!! Psolver() developped for BigDFT.
!! It can compute the xc energy and potential if required. This computation is
!! built on the drivexc() routine of ABINIT but access it directly from real
!! space. The present routine is a real space counter part to rhotoxc().
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR,TRangel).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=MPI-parallelisation information.
!!  rhor(nfft,nspden)=electron density in real space in electrons/bohr**3
!!
!! OUTPUT
!!  enhartr=returned Hartree energy (hartree).
!!  enxc=returned exchange and correlation energy (hartree).
!!  envxc=returned energy of the Vxc potential (hartree).
!!  vhartr(nfft)=Hartree potential.
!!  vxc(nfft,nspden)=xc potential
!!  vxcavg=<Vxc>=unit cell average of Vxc = (1/ucvol) Int [Vxc(r) d^3 r].
!!
!! NOTE
!!  In psolver, with nspden == 2, rhor(:,1) = density up and
!!                                rhor(:,2) = density down.
!!  But in ABINIT (dtset%usewvl != 1) rhor(:,1) = total density and
!!                                    rhor(:,2) = density up .
!!  In ABINIT (dtset%usewvl != 1), the same convention is used as in psolver.
!!
!! PARENTS
!!      energy,rhotov,scfcv,setvtr
!!
!! CHILDREN
!!      h_potential,mean_fftr,metric,mkdenpos,psolver_kernel,wrtout
!!      wvl_rhov_abi2big,xc_potential
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psolver_rhohxc(enhartr, enxc, envxc, icoulomb, ixc, &
& mpi_enreg, nfft, ngfft, nhat,nhatdim,&
& nscforder, nspden, n3xccc, rhor, rprimd,&
& usexcnhat,usepaw,usewvl,vhartr, vxc, vxcavg, wvl,wvl_den,wvl_e,&
& xccc3d,xclevel,xc_denpos)

! use defs_basis,only: std_out,std_out_default
 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_abi2big
 use m_cgtools

 use m_xmpi, only: xmpi_comm_rank,xmpi_comm_size,xmpi_sum
 use m_geometry, only : metric

#if defined HAVE_BIGDFT
 use BigDFT_API, only : XC_potential,ELECTRONIC_DENSITY,coulomb_operator
 use poisson_solver, only : H_potential
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psolver_rhohxc'
 use interfaces_14_hidewrite
 use interfaces_41_xc_lowlevel
 use interfaces_62_poisson, except_this_one => psolver_rhohxc
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  !scalars
  integer, intent(in)           :: nhatdim,nspden,n3xccc
  integer, intent(in)           :: nfft, icoulomb, ixc, nscforder, usewvl
  integer,intent(in)            :: usexcnhat,usepaw,xclevel
  real(dp),intent(in)           :: rprimd(3,3)
  real(dp), intent(in)          :: xc_denpos
  real(dp), intent(out)         :: enxc, envxc, enhartr, vxcavg
  type(mpi_type), intent(in) :: mpi_enreg
  type(wvl_internal_type), intent(in) :: wvl
  type(wvl_denspot_type), intent(inout) :: wvl_den
  type(wvl_energy_terms), intent(inout) :: wvl_e
  !arrays
  integer, intent(in)    :: ngfft(18)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(in) :: nhat(nfft,nspden*nhatdim)
  real(dp),intent(inout) :: rhor(nfft, nspden)
  real(dp),intent(out)   :: vhartr(nfft)
  real(dp),intent(out)   :: vxc(nfft, nspden)

  !Local variables-------------------------------
#if defined HAVE_BIGDFT
! n_c and \hat{n} can be added/rested inside bigdft by passing
! them as pointers (rhocore and rhohat):
  logical, parameter :: add_n_c_here=.true.  !Add n_c here or inside bigdft
  logical, parameter :: rest_hat_n_here=.true.  !Rest \hat{n} here or inside bigdft
  !scalars
  integer :: me,nproc,comm
  integer :: ifft,ispin
  integer :: iwarn, opt_mkdenpos
  integer :: nfftot,ngrad
  integer :: n1i,n2i,n3d,n3i
  real(dp) :: tmpDown, tmpUp, tmpPot,ucvol,ucvol_local
  logical :: sumpion,test_nhat,use_psolver=.false.
  character(len=500) :: message
  character(len = 1) :: datacode, bndcode
  !arrays
  real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
  real(dp) :: hgrid(3)
  real(dp) :: vxcmean(1)
  real(dp), pointer :: rhocore(:,:,:,:),rhohat(:,:,:,:)
  real(dp), pointer :: pot_ion(:,:,:,:),rhonow(:,:)
  real(dp), dimension(6) :: xcstr
  type(coulomb_operator) ::  kernel
#endif

! *********************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 nfftot=PRODUCT(ngfft(1:3))
 comm=mpi_enreg%comm_fft
 if(usewvl==1) comm=mpi_enreg%comm_wvl
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 if(n3xccc>0) then
   if(nfft .ne. n3xccc)then
     write(message,'(a,a,a,2(i0,1x))')&
&     'nfft and n3xccc should be equal,',ch10,&
&     'however, nfft and n3xccc=',nfft,n3xccc
     MSG_BUG(message)
   end if
 end if
 if(nspden==4) then
   MSG_ERROR('nspden==4 not coded yet')
 end if

 if (ixc==0) then
   vxcavg=zero
   test_nhat=.false.

!  No xc at all is applied (usually for testing)
   MSG_WARNING('Note that no xc is applied (ixc=0).')

 else if (ixc/=20) then

!  ngrad=1 is for LDAs or LSDs, ngrad=2 is for GGAs
   ngrad=1;if(xclevel==2)ngrad=2
!  ixc 31 to 34 are for mgga test purpose only (fake functionals based on LDA but need the gradients too)
   if(ixc>=31 .and. ixc<=34)ngrad=2
!  Test: has a compensation density to be added/substracted (PAW) ?
!  test_nhat=((nhatdim==1).and.(usexcnhat==0.or.(ngrad==2.and.nhatgrdim==1)))
   test_nhat=((nhatdim==1).and.(usexcnhat==0))
 end if


!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (icoulomb == 0) then
!  The kernel is built with 'P'eriodic boundary counditions.
   bndcode = 'P'
 else if (icoulomb == 1) then
!  The kernel is built with 'F'ree boundary counditions.
   bndcode = 'F'
 else if (icoulomb == 2) then
!  The kernel is built with 'S'urface boundary counditions.
   bndcode = 'S'
 end if

!This makes the tests fail?
!For NC and n_c=0, call psolver, which uses less memory:
!if(usepaw==0 .and. n3xccc==0) use_psolver=.true.

 if(nspden > 2)then
   write(message, '(a,a,a,i0)' )&
&   'Only non-spin-polarised or collinear spin is allowed,',ch10,&
&   'while the argument nspden = ', nspden
   MSG_ERROR(message)
 end if

!We do the computation.
 write(message, "(A,A,A,3I6)") "psolver_rhohxc(): compute potentials (Vhartree and Vxc)...", ch10, &
& " | dimension:", ngfft(1:3)
 call wrtout(std_out, message,'COLL')

 if(usewvl==1) then
   hgrid=(/wvl_den%denspot%dpbox%hgrids(1),wvl_den%denspot%dpbox%hgrids(2),wvl_den%denspot%dpbox%hgrids(3)/)
 else
   hgrid=(/ rprimd(1,1) / ngfft(1), rprimd(2,2) / ngfft(2), rprimd(3,3) / ngfft(3) /)
 end if

 if (usewvl == 0) then
!  We get the kernel.
   call psolver_kernel( hgrid, 2, icoulomb, me, kernel, comm, ngfft, nproc, nscforder)
 elseif(usewvl==1) then
!  In this case, the kernel is already computed.
!  We just shallow copy it.
   kernel = wvl_den%denspot%pkernel
 end if

 if(usewvl==1) then
   if(wvl_den%denspot%rhov_is .ne. ELECTRONIC_DENSITY) then
     message= "psolver_rhohxc: rhov should contain the electronic density"
     MSG_ERROR(message)
   end if
 end if

 if(usewvl==1) then
   n1i=wvl%Glr%d%n1i; n2i=wvl%Glr%d%n2i; n3i=wvl%Glr%d%n3i
   n3d=wvl_den%denspot%dpbox%n3d
 else
   n1i=ngfft(1); n2i=ngfft(2) ; n3i=ngfft(3)
   n3d=ngfft(13)
 end if

 if (usewvl == 0) then
!  ucvol_local=product(hgrid)*half**3*real(n1i*n2i*n3i,dp)
!  write(*,*)'hgrid, n1i,n2i,n3i',hgrid,ngfft(1:3)
!  write(*,*)'ucvol_local',ucvol_local
   ucvol_local = ucvol
!  write(*,*)'ucvol_local',ucvol_local
 else
!  We need to tune the volume when wavelets are used because, not
!  all FFT points are used.
!  ucvol_local = (half * dtset%wvl_hgrid) ** 3 * ngfft(1)*ngfft(2)*ngfft(3)
   ucvol_local = product(wvl_den%denspot%dpbox%hgrids) * real(product(wvl_den%denspot%dpbox%ndims), dp)
 end if

!Core density array
 if(n3xccc==0 .or. add_n_c_here) then
   nullify(rhocore)
!  Pending, next line should follow the same logic that the rest
   if(usewvl==1 .and. usepaw==0)  rhocore=> wvl_den%denspot%rho_C
 else
   if(usepaw==1) then
     ABI_ALLOCATE(rhocore,(n1i,n2i,n3d,1)) !not spin dependent
     call wvl_rhov_abi2big(1,xccc3d,rhocore)

!    Make rhocore positive to avoid numerical instabilities in V_xc
     iwarn=0 ; opt_mkdenpos=0
     call mkdenpos(iwarn, nfft, nspden, opt_mkdenpos, rhocore, tol20 )
   end if
 end if

!write(*,*)'psolver_rhohxc, erase me, set rhocore=0'
!if( associated(wvl_den%denspot%rho_C))wvl_den%denspot%rho_C=zero
!if(associated(rhocore))rhocore=zero

!Rhohat array:
 if(test_nhat .and. .not. rest_hat_n_here) then
!  rhohat => nhat !do not know how to point 4 index to 2 index
!  here we have to copy since convention for spin changes.
   ABI_ALLOCATE(rhohat,(n1i,n2i,n3d,nspden))
   call wvl_rhov_abi2big(1,nhat,rhohat)
 else
   nullify(rhohat)
 end if

!Data are always distributed when using the wavelets, even if nproc = 1.
!The size given is the complete size of the box, not the distributed size
!stored in ngfft.
 if (nproc > 1 .or. usewvl > 0) then
   datacode = 'D'
 else
   datacode = 'G'
 end if

!If usewvl=1, vpsp(or v_ext) will be summed to vhartree
 if(usewvl==1) then
   pot_ion=>wvl_den%denspot%V_ext
   sumpion=.false.
!  Note:
!  if sumpion==.true.
!  call wvl_newvtr in setvtr and rhotov
!  if sumpion==.false.
!  modify setvtr and rhotov to not use wvl_newvtr and follow the normal ABINIT flow.
 else
!  This is not allowed
!  pot_ion=>vxc !this is a dummy variable here
   sumpion=.false.
 end if


!To make this work, make sure that xc_init has been called
!in gstate.
 if(.not. use_psolver) then
!  T.Rangel:
!  Use this approach for PAW and sometimes for NC since
!  in psolver() the core density is not added.
!
!  PAW case:
!  It is important to call H_potential before XC_potential:
!  In XC_potential, if test_nhat, we do:
!  1) rhor=rhor-rhohat,
!  2) makepositive(rhor,tol20)
!  3) after Psolver, we do rhor=rhor+rhohat,
!  I found that rhor at input and output are slightly different,
!  These differences lead to a difference of ~0.01 hartree in V_hartree.
!  If PAW, substract compensation density from effective density:
!  - if GGA, because nhat gradients are computed separately
!  - if nhat does not have to be included in XC

!  save rhor in rhonow to avoid modifying it.
   ABI_ALLOCATE(rhonow,(nfft,nspden))
!  copy rhor into rhonow:
!  ABINIT convention is followed: (ispin=1: for spin up + spin down)
   do ispin=1,nspden
     do ifft=1,nfft
       rhonow(ifft,ispin)=abs(rhor(ifft,ispin))+1.0d-20
     end do
   end do

   if(usewvl==1) then
     call H_potential(datacode,&
&     kernel,rhonow,pot_ion,enhartr,&
&     zero,sumpion)
   else
!    Vxc is passed as a dummy argument
     call H_potential(datacode,&
&     kernel,rhonow,vxc,enhartr,&
&     zero,sumpion)
   end if
!
   do ifft=1,nfft
     vhartr(ifft)=rhonow(ifft,1)
   end do
!  write(*,*)'erase me psolver_rhohxc l350, set vhartr=0'
!  vhartr=zero ; enhartr=zero
!
!  Since rhonow was modified inside H_potential:
!  copy rhor again into rhonow following the BigDFT convention:
   call wvl_rhov_abi2big(1,rhor,rhonow)

!  Add n_c here:
   if(n3xccc>0 .and. add_n_c_here) then
     do ispin=1,nspden
       do ifft=1,nfft
         rhonow(ifft,ispin)=rhonow(ifft,ispin)+xccc3d(ifft)
       end do
     end do
   end if
!  Remove \hat{n} here:
   if(test_nhat .and. rest_hat_n_here) then
     do ispin=1,nspden
       do ifft=1,nfft
         rhonow(ifft,ispin)=rhonow(ifft,ispin)-nhat(ifft,ispin)
       end do
     end do
   end if

!  Make the density positive everywhere (but do not care about gradients)
   iwarn=0 ; opt_mkdenpos=0
   call mkdenpos(iwarn, nfft, nspden, opt_mkdenpos, rhonow, xc_denpos)
!  do ispin=1,nspden
!  do ifft=1,nfft
!  rhonow(ifft,ispin)=abs(rhonow(ifft,ispin))+1.0d-20
!  end do
!  end do

!  If PAW, substract compensation density from effective density:
!  - if GGA, because nhat gradients are computed separately
!  - if nhat does not have to be included in XC
   if (test_nhat .and. .not. rest_hat_n_here) then

     call XC_potential(bndcode,datacode,me,nproc,comm,&
&     n1i,n2i,n3i,&
&     wvl_den%denspot%xc,hgrid(1),hgrid(2),hgrid(3),&
&     rhonow,enxc,envxc,nspden,rhocore,&
&     vxc,xcstr,rhohat=rhohat)

   else

     call XC_potential(bndcode,datacode,me,nproc,comm,&
&     n1i,n2i,n3i,&
&     wvl_den%denspot%xc,hgrid(1),hgrid(2),hgrid(3),&
&     rhonow,enxc,envxc,nspden,rhocore,&
&     vxc,xcstr)

   end if

!  write(*,*)'psolver_rhohxc: erase me, set vxc=0'
!  vxc=zero
!  enxc=zero
!  envxc=zero

!  deallocate temporary array
   ABI_DEALLOCATE(rhonow)

 else
!  NC case: here we optimize memory, and we reuse vhartree to store rhor:

!  We save total rhor in vhartr
   do ifft=1,nfft
     vhartr(ifft)  = rhor(ifft, 1)
   end do

!  In non-wavelet case, we change the rhor values.
   if (nspden == 2) then
     do ifft = 1, nfft
!      We change rhor for psolver call.
       tmpDown = rhor(ifft, 1) - rhonow(ifft, 2)
       tmpUp   = rhor(ifft, 2)
       rhor(ifft, 1) = tmpUp
       rhor(ifft, 2) = tmpDown
     end do
   end if
!  Make the density positive everywhere (but do not care about gradients)
   iwarn=0 ; opt_mkdenpos=0
   call mkdenpos(iwarn, nfft, nspden, opt_mkdenpos, rhor, xc_denpos)
!  do ispin=1,nspden
!  do ifft=1,nfft
!  rhor(ifft,ispin)=abs(rhor(ifft,ispin))+1.0d-20
!  end do
!  end do

!  Call Poisson solver, here rhor(:,1) will contain Vhartree at output
!  This does not compile, check mklocl_realspace where it do work.
!   call psolver(bndcode, datacode, me, nproc, n1i, &
!&   n2i,n3i, ixc, hgrid(1), hgrid(2), hgrid(3), &
!&   rhor, kernel, vxc, enhartr, enxc, envxc, 0.d0, .false., nspden)

!  PSolver work in place, we set back the rhor values.
   do ifft = 1, nfft, 1
     tmpPot     = rhor(ifft, 1)
!    Rhor total was saved in vhartr and current rhor(:,2) is down spin
     rhor(ifft, 1) = vhartr(ifft)
     if (nspden == 2) rhor(ifft, 2) = rhor(ifft, 1) - rhor(ifft, 2)
     vhartr(ifft)  = tmpPot
   end do
 end if

!Pass vhartr and vxc to BigDFT objects (useless?)
!if(usewvl==1) then
!  write(message, '(a,a,a,a)' ) ch10, ' rhotoxc_wvlpaw : but why are you copying me :..o('
! call wvl_vhartr_abi2big(1,vhartr,wvl_den)
!  (this can be commented out, since we do not use denspot%v_xc
! call wvl_vxc_abi2big(1,vxc,wvl_den)
!end if

!Compute vxcavg
 call mean_fftr(vxc, vxcmean, nfft, nfftot, nspden,mpi_comm_sphgrid=comm)
 vxcavg = vxcmean(1)

!Pass energies to wvl object:
 if(usewvl==1) then
   wvl_e%energs%eh  = enhartr
   wvl_e%energs%exc = enxc
   wvl_e%energs%evxc= envxc
 end if

!Nullify pointers and deallocate arrays
 if(test_nhat .and. .not. rest_hat_n_here) then
!  if(nspden==2) ABI_DEALLOCATE(rhohat)
   ABI_DEALLOCATE(rhohat)
   if(associated(rhohat)) nullify(rhohat)
 end if
 if( n3xccc>0 .and. .not. add_n_c_here) then
   if(usepaw==1) then
     ABI_DEALLOCATE(rhocore)
   end if
 end if
 if(associated(rhocore))  nullify(rhocore)
 if(associated(pot_ion)) nullify(pot_ion)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) nhatdim,nspden,n3xccc,nfft,icoulomb,ixc,nscforder,usewvl,&
& usexcnhat,usepaw,xclevel,rprimd(1,1),xc_denpos,enxc,envxc,enhartr,vxcavg,mpi_enreg%nproc,&
& wvl%h(1),wvl_den%symObj,wvl_e%energs,ngfft(1),xccc3d(1),nhat(1,1),rhor(1,1),vhartr(1),vxc(1,1)
#endif

 DBG_EXIT("COLL")

end subroutine psolver_rhohxc
!!***
