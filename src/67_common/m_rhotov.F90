!!****m* ABINIT/m_rhotov
!! NAME
!!  m_rhotov
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (XG, GMR, MT, EB)
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

module m_rhotov

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_ab7_mixing
 use m_abi2big
 use m_xmpi
 use m_cgtools
 use m_xcdata
 use m_dtset

 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : xred2xcart
 use m_energies,         only : energies_type
 use m_electronpositron, only : electronpositron_type, electronpositron_calctype, rhohxcpositron
 use libxc_functionals,  only : libxc_functionals_is_hybrid
 use m_spacepar,         only : hartre
 use m_dens,             only : constrained_dft_t,mag_penalty,constrained_residual
 use m_rhotoxc,          only : rhotoxc
 use m_xchybrid,         only : xchybrid_ncpp_cc
 use m_psolver,          only : psolver_rhohxc
 use m_wvl_psi,          only : wvl_psitohpsi

 implicit none

 private
!!***

 public :: rhotov
!!***

contains
!!***

!!****f* ABINIT/rhotov
!! NAME
!! rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given total density
!! the trial (local) potential and the residual potential.
!!
!! INPUTS
!!  [add_tfw]=flag controling the addition of Weiszacker gradient correction to Thomas-Fermi kin energy
!!  constrained_dft <type(constrained_dft_t>=data for constrained dft calculations
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | spinmagntarget=input variable that governs fixed moment calculation
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | ntypat=number of types of atoms in unit cell.
!!   | occopt=option for occupancies
!!   | typat(natom)=type (integer) for each atom
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc, see rhotoxc.F90 for a description
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  optxc=option to be used for the call to rhotoxc
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!   | definition for spin components:
!!   | case of nspden = 2
!!   |      rhor(:,1) => rho_up + rho_dwn
!!   |      rhor(:,2) => rho_up
!!   | case of nspden = 4
!!   |      rhor(:,1)   => rho_upup + rho_dwndwn
!!   |      rhor(:,2:4) => {m_x,m_y,m_z}
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  [taur(nfftf,nspden*dtset%usekden)]=array for kinetic energy density
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfft)=array for holding local psp
!!  [vxc_hybcomp(nfft,nspden)= compensation xc potential (Hartree) in case of hybrids] Optional output
!!       i.e. difference between the hybrid Vxc at fixed density and the auxiliary Vxc at fixed density
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  ==== if optres==0
!!    vtrial(nfft,nspden)= old value of trial potential
!!
!! OUTPUT
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_hartree=Hartree part of total energy (hartree units)
!!   | e_xc=exchange-correlation energy (hartree)
!!   | In case of hybrid compensation algorithm:
!!   | e_hybcomp_v=self-consistent potential compensation term for the exchange-correlation energy (hartree)
!!  ==== if optene==0.or.2
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if optene==1.or.2
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!  intgres(nspden,ngrcondft)=integrated residuals from constrained DFT. They are also Lagrange parameters, or gradients with respect to constraints.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if optxc==2.
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxc(nfft,nspden)=Vxc(r) (already computed above; gets recomputed below too)
!!  vxcavg=mean of the vxc potential
!!  ==== if optres==0
!!    vresidnew(nfft,nspden)=potential residual
!!    vnew_mean(nspden)=mean of the potential formed from vpsp, vhartr and vxc, might be spin-dependent
!!    vres_mean(nspden)=mean of the potential residual, might be spin-dependent
!!    vres2=square of the norm of the residual
!!    [vxctau(nfftf,dtset%nspden,4*dtset%usekden)]=derivative of XC energy density with respect to
!!      kinetic energy density (metaGGA cases) (optional output)
!!    [vtauresid(nfft,nspden*dtset%usekden)]=array for vxctau residue (see vtau below))
!!
!! SIDE EFFECTS
!! Input/Output:
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  vhartr(nfft)=array for holding Hartree potential
!!  ==== if optres==1
!!    vtrial(nfft,nspden)= new value of trial potential
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      dotprod_vn,hartre,mag_penalty,mean_fftr,psolver_rhohxc,rhohxcpositron
!!      rhotoxc,sqnorm_v,timab,wvl_psitohpsi,wvl_vtrial_abi2big,xcdata_init
!!      xchybrid_ncpp_cc,xred2xcart
!!
!! SOURCE

subroutine rhotov(constrained_dft,dtset,energies,gprimd,gsqcut,intgres,istep,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&
&  rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&
&  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,wvl,xccc3d,xred,&
&  electronpositron,taur,vxc_hybcomp,vxctau,vtauresid,add_tfw,xcctau3d) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfft,nhatgrdim,nkxc,optene,optres,optxc,usepaw,istep
 integer,intent(in) :: usexcnhat
 logical,intent(in),optional :: add_tfw
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vres2,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(constrained_dft_t),intent(inout) :: constrained_dft
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(wvl_data), intent(inout) :: wvl
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: intgres(:,:) ! (nspden,ngrcondft) ngrcondft=natom when condft is activated
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vnew_mean(dtset%nspden)
 real(dp),intent(out) :: vres_mean(dtset%nspden),vresidnew(nfft,dtset%nspden)
 real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(inout),optional :: vtauresid(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(out),optional,target :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
 real(dp),intent(out),optional :: vxc_hybcomp(:,:) ! (nfft,nspden)
 real(dp),intent(out),optional :: xcctau3d(n3xccc)

!Local variables-------------------------------
!scalars
 integer :: nk3xc,ifft,ipositron,ispden,nfftot,offset
 integer :: mpi_comm_sphgrid,ixc_current
!integer :: ii,jj,kk,ipt,nx,ny,nz           !SPr: debug
!real(dp):: rx,ry,rz                        !SPr: debug
 real(dp) :: doti,e_xcdc_vxctau
 logical :: add_tfw_,calc_xcdc,non_magnetic_xc,with_vxctau
 logical :: is_hybrid_ncpp,wvlbigdft=.false.
 type(xcdata_type) :: xcdata
!arrays
 real(dp) :: evxc,tsec(2),vmean(dtset%nspden),vzeeman(dtset%nspden)
 real(dp),target :: vxctau_dum(0,0,0)
 real(dp),allocatable :: rhowk(:,:),v_constr_dft_r(:,:),vnew(:,:),xcart(:,:)
 real(dp),pointer :: vxctau_(:,:,:)
!real(dp),allocatable :: vzeemanHarm(:,:)   !SPr: debug Zeeman field q/=0 real space

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(940,1,tsec)

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = (present(vxctau).and.present(taur).and.(dtset%usekden/=0))
 vxctau_ => vxctau_dum ; if (with_vxctau) vxctau_ => vxctau
 if (with_vxctau.and.optres==0.and.(.not.present(vtauresid))) then
   MSG_BUG('need vtauresid!')
 end if

!Check if we're in hybrid norm conserving pseudopotential with a core correction
 is_hybrid_ncpp=(usepaw==0 .and. n3xccc/=0 .and. &
& (dtset%ixc==41.or.dtset%ixc==42.or.libxc_functionals_is_hybrid()))

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

!mpi communicator for spherical grid
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) mpi_comm_sphgrid=mpi_enreg%comm_wvl

!Get size of FFT grid
 nfftot=PRODUCT(ngfft(1:3))

 ipositron=0;if (present(electronpositron)) ipositron=electronpositron_calctype(electronpositron)
 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw

!------Compute Hartree and xc potentials----------------------------------

!allocate vnew here.
!In wvl: vnew is used at call to wvl_psitohpsi
 if (optres==0) then
   ABI_ALLOCATE(vnew,(nfft,dtset%nspden))
   vmean(:)=zero ; vnew_mean(:)=zero
 end if

 if (ipositron/=1) then
!  if metaGGA, save current value of vxctau potential
   if (with_vxctau) vtauresid(:,:)=vxctau(:,:,1)
!  Compute xc potential (separate up and down if spin-polarized)
   if (dtset%icoulomb == 0 .and. dtset%usewvl == 0) then
     call hartre(1,gsqcut,usepaw,mpi_enreg,nfft,ngfft,rhog,rprimd,vhartr)
     !Use the proper exchange_correlation energy : either the origin one, or the auxiliary one
     ixc_current=dtset%ixc
     if(mod(dtset%fockoptmix,100)==11)ixc_current=dtset%auxc_ixc
     call xcdata_init(xcdata,dtset=dtset,ixc=ixc_current)
     non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

!    Use the periodic solver to compute Hxc.
     call timab(941,1,tsec)
     nk3xc=1
     if (ipositron==0) then
       if(.not.is_hybrid_ncpp .or. mod(dtset%fockoptmix,100)==11)then
         call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
&         nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,optxc,&
&         rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,&
&         taur=taur,vhartr=vhartr,vxctau=vxctau_,add_tfw=add_tfw_,xcctau3d=xcctau3d)
         if(mod(dtset%fockoptmix,100)==11)then
           energies%e_xc=energies%e_xc*dtset%auxc_scal
           vxc(:,:)=vxc(:,:)*dtset%auxc_scal
         end if
       else
         call xchybrid_ncpp_cc(dtset,energies%e_xc,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,&
&         strsxc,vxcavg,xccc3d,vxc=vxc)
       end if
     else
       call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
&       nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,optxc,&
&       rhor,rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,&
&       taur=taur,vhartr=vhartr,vxctau=vxctau_,add_tfw=add_tfw_,&
&       electronpositron=electronpositron,xcctau3d=xcctau3d)
     end if
     call timab(941,2,tsec)
   elseif (.not. wvlbigdft) then
!    Use the free boundary solver.
     call timab(943,1,tsec)
     call psolver_rhohxc(energies%e_hartree, energies%e_xc, evxc, &
&     dtset%icoulomb, dtset%ixc, &
&     mpi_enreg, nfft, &
&     ngfft, nhat,usepaw,&
&     dtset%nscforder, dtset%nspden, n3xccc, rhor,rprimd,&
&     usexcnhat,dtset%usepaw,dtset%usewvl,vhartr, vxc, vxcavg,&
&     wvl%descr,wvl%den,wvl%e,&
&     xccc3d,dtset%xclevel,dtset%xc_denpos)
     call timab(943,2,tsec)
   end if
!  For icoulomb==0 and usewvl Ehartree is calculated in psolver_rhohxc().
!  For PAW we recalculate this since nhat was not taken into account
!  in psolver_rhohxc: E_H= int v_H (n+nhat) dr
   if(.not. wvlbigdft .and. (dtset%icoulomb==0 .or. dtset%usepaw==1 ) ) then
     call timab(942,1,tsec)
     call dotprod_vn(1,rhor,energies%e_hartree,doti,nfft,nfftot,1,1,vhartr,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
     energies%e_hartree=half*energies%e_hartree
     call timab(942,2,tsec)
   end if
 else
   call timab(944,1,tsec)
   energies%e_hartree=zero;energies%e_xc=zero
   call rhohxcpositron(electronpositron,gprimd,kxc,mpi_enreg,nfft,ngfft,nhat,nkxc,dtset%nspden,n3xccc,&
&   dtset%paral_kgb,rhor,strsxc,ucvol,usexcnhat,usepaw,vhartr,vxc,vxcavg,xccc3d,dtset%xc_denpos)
   call timab(944,2,tsec)
 end if

 call timab(945,1,tsec)
 if (ipositron/=0) then
   call dotprod_vn(1,rhor,electronpositron%e_hartree,doti,&
&   nfft,nfftot,1,1,electronpositron%vha_ep,ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
   vhartr=vhartr+electronpositron%vha_ep
 end if

!------Compute parts of total energy depending on potentials--------

 if ( (optene==0.or.optene==2 ).and. .not. wvlbigdft) then
!  Compute local psp energy energies%e_localpsp
   call dotprod_vn(1,rhor,energies%e_localpsp,doti,nfft,nfftot,1,1,vpsp,ucvol,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)
 end if

 if(mod(dtset%fockoptmix,100)==11)then
   if (.not. wvlbigdft) then
!    Compute second compensation energy for hybrid functionals
     call dotprod_vn(1,rhor,energies%e_hybcomp_v,doti,nfft,nfftot,1,1,vxc_hybcomp,ucvol,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
   end if
 end if

 calc_xcdc=.false.
 if (optene==1.or.optene==2) calc_xcdc=.true.
 if (dtset%usewvl==1.and.dtset%nnsclo>0) calc_xcdc=.true.
 if (wvlbigdft) calc_xcdc=.false.
 if (dtset%usefock==1) calc_xcdc=.true.

 if (calc_xcdc) then

!  Compute double-counting XC energy energies%e_xcdc
   if (ipositron/=1) then
     if (usepaw==0.or.usexcnhat/=0) then
       call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol,&
&       mpi_comm_sphgrid=mpi_comm_sphgrid)
     else
       ABI_ALLOCATE(rhowk,(nfft,dtset%nspden))
       rhowk=rhor-nhat
       call dotprod_vn(1,rhowk,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol,&
&       mpi_comm_sphgrid=mpi_comm_sphgrid)
       ABI_DEALLOCATE(rhowk)
     end if
     if (with_vxctau)then
       call dotprod_vn(1,taur,e_xcdc_vxctau,doti,nfft,nfftot,dtset%nspden,1,vxctau(:,:,1),&
&       ucvol,mpi_comm_sphgrid=mpi_comm_sphgrid)
       energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
     end if
     if (ipositron==2) energies%e_xcdc=energies%e_xcdc-electronpositron%e_xcdc
   else
     energies%e_xcdc=zero
   end if

 end if

!------Produce residual vector and square norm of it-------------
!(only if requested ; if optres==0)

!Set up array for Zeeman field
!EB vzeeman(:) = factor*( Hz, Hx+iHy; Hx-iHy, -Hz)
!EB factor = -g/2 * mu_B * mu_0 = -1/2*B in a.u.
!EB <-- vzeeman might have to be allocated correctly --> to be checked
! vzeeman = 1/2 ( -B_z, -B_x + iB_y ; -B_x - iB_y , B_z)
 vzeeman(:) = zero
! ABI_ALLOCATE(vzeemanHarm,(nfft,dtset%nspden))  ! SPr: debug stuff
! vzeemanHarm(:,:) = zero                        !
 if (any(abs(dtset%zeemanfield(:))>tol8)) then
   if(dtset%nspden==2)then
!    EB The collinear case has to be checked :
!    EB Is it vzeeman(1) or (2) that has to be added here? to be checked in setvtr and energy as well
!    SPr: the density components are: rhor(1) => n_upup + n_dwndwn
!                                     rhor(2) => n_upup
!         the convention for the potential components is different:
!                                     v(1)    => v_upup
!                                     v(2)    => v_dndn
!         verified by comparing collinear and non-collinear calculations

     vzeeman(1) =-half*dtset%zeemanfield(3)  ! v_upup
     vzeeman(2) = half*dtset%zeemanfield(3)  ! v_dndn

     !vzeeman(1) = zero  ! v_upup
     !vzeeman(2) = zero  ! v_dndn

     !nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
     !do kk=0,nz-1
     !  do jj=0,ny-1
     !    do ii=0,nx-1
     !      ipt=1+ii+nx*(jj+ny*kk)
     !      !rx=(dble(ii)/nx)*rprimd(1,1)+(dble(jj)/ny)*rprimd(1,2)+(dble(kk)/nz)*rprimd(1,3)
     !      !ry=(dble(ii)/nx)*rprimd(2,1)+(dble(jj)/ny)*rprimd(2,2)+(dble(kk)/nz)*rprimd(2,3)
     !      !rz=(dble(ii)/nx)*rprimd(3,1)+(dble(jj)/ny)*rprimd(3,2)+(dble(kk)/nz)*rprimd(3,3)
     !      vzeemanHarm(ipt,1)= -half*dtset%zeemanfield(3)*cos(2*PI*(dble(ii)/dble(nx)))
     !      vzeemanHarm(ipt,2)=  half*dtset%zeemanfield(3)*cos(2*PI*(dble(ii)/dble(nx)))
     !    end do
     !  end do
     !end do

   else if(dtset%nspden==4)then

     vzeeman(1)=-half*dtset%zeemanfield(3)    ! v_upup
     vzeeman(2)= half*dtset%zeemanfield(3)    ! v_dndn
     vzeeman(3)=-half*dtset%zeemanfield(1)    ! Re(v_updn)
     vzeeman(4)= half*dtset%zeemanfield(2)    ! Im(v_updn)

     !vzeeman(1)=0.0
     !vzeeman(2)=0.0
     !vzeeman(3)=0.0
     !vzeeman(4)=0.0

     !nx=ngfft(1); ny=ngfft(2); nz=ngfft(3)
     !do kk=0,nz-1
     !  do jj=0,ny-1
     !    do ii=0,nx-1
     !      ipt=1+ii+nx*(jj+ny*kk)
     !      !rx=(dble(ii)/nx)*rprimd(1,1)+(dble(jj)/ny)*rprimd(1,2)+(dble(kk)/nz)*rprimd(1,3)
     !      !ry=(dble(ii)/nx)*rprimd(2,1)+(dble(jj)/ny)*rprimd(2,2)+(dble(kk)/nz)*rprimd(2,3)
     !      !rz=(dble(ii)/nx)*rprimd(3,1)+(dble(jj)/ny)*rprimd(3,2)+(dble(kk)/nz)*rprimd(3,3)
     !      vzeemanHarm(ipt,1)= -half*dtset%zeemanfield(3)*cos(2*PI*(dble(ii)/dble(nx)))
     !      vzeemanHarm(ipt,2)=  half*dtset%zeemanfield(3)*cos(2*PI*(dble(ii)/dble(nx)))
     !      vzeemanHarm(ipt,3)= -half*dtset%zeemanfield(1)*cos(2*PI*(dble(ii)/dble(nx)))
     !      vzeemanHarm(ipt,4)=  half*dtset%zeemanfield(2)*cos(2*PI*(dble(ii)/dble(nx)))
     !    end do
     !  end do
     !end do

   end if
 end if

!Compute the constrained potential for the magnetic moments
 ABI_ALLOCATE(v_constr_dft_r, (nfft,dtset%nspden))
 v_constr_dft_r = zero
 if (dtset%magconon==1.or.dtset%magconon==2) then
   call mag_penalty(constrained_dft,mpi_enreg,rhor,v_constr_dft_r,xred)
 end if

 if (optres==0) then


!  ------ Compute potential residual -------------

   if (.not. wvlbigdft) then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=1,min(dtset%nspden,2)
       do ifft=1,nfft
         vnew(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)+v_constr_dft_r(ifft,ispden)
         !vnew(ifft,ispden)=vnew(ifft,ispden)+vzeemanHarm(ifft,ispden)
         if(mod(dtset%fockoptmix,100)==11)vnew(ifft,ispden)=vnew(ifft,ispden)+vxc_hybcomp(ifft,ispden)
         vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
       end do
     end do
     if(dtset%nspden==4)then
!$OMP PARALLEL DO COLLAPSE(2)
       do ispden=3,4
         do ifft=1,nfft
           vnew(ifft,ispden)=vxc(ifft,ispden)+vzeeman(ispden)+v_constr_dft_r(ifft,ispden)
           !vnew(ifft,ispden)=vnew(ifft,ispden)+vzeemanHarm(ifft,ispden)
           if(mod(dtset%fockoptmix,100)==11)vnew(ifft,ispden)=vnew(ifft,ispden)+vxc_hybcomp(ifft,ispden)
           vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
         end do
       end do
     end if

     !If constrained_dft, must take into account the constraints, and recompute the residual and the new potential
     if( any(dtset%constraint_kind(:)/=0))then
       call constrained_residual(constrained_dft,energies%e_constrained_dft,intgres,mpi_enreg,rhor,vresidnew,xred)
       vnew(:,1:dtset%nspden)=vtrial(:,1:dtset%nspden)+vresidnew(:,1:dtset%nspden)
     endif

     offset   = 0

     if (dtset%iscf==0) vtrial=vnew

!    Pass vtrial to BigDFT object
     if(dtset%usewvl==1) then
       call wvl_vtrial_abi2big(1,vnew,wvl%den)
!      call wvl_vtrial_abi2big(1,vtrial,wvl%den)
     end if

   else
!    Compute with covering comms the different part of the potential.
!    only for wvlbigdft
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     call wvl_psitohpsi(dtset%diemix,energies%e_exactX, energies%e_xc, energies%e_hartree, &
&     energies%e_kinetic, energies%e_localpsp, energies%e_nlpsp_vfock, energies%e_sicdc, &
&     istep + 1, 1, dtset%iscf, mpi_enreg%me_wvl, dtset%natom, dtset%nfft,&
&     mpi_enreg%nproc_wvl, dtset%nspden, &
&     vres2, .true., energies%e_xcdc, wvl,&
&     wvlbigdft, xcart, strsxc,vtrial=vnew,vxc=vxc)
     ABI_DEALLOCATE(xcart)

     vresidnew = vnew - vtrial
     vtrial = vnew

     call mean_fftr(vxc, vmean(1:1),  nfft, nfftot, dtset%nspden,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     vxcavg = vmean(1)
     offset = 0
   end if

!  Compute mean values of potential and residual
   call mean_fftr(vnew(1+offset, 1),vnew_mean,nfft,nfftot,dtset%nspden,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)
   call mean_fftr(vresidnew(1+offset, 1),vmean,nfft,nfftot,dtset%nspden,&
&   mpi_comm_sphgrid=mpi_comm_sphgrid)

   ABI_DEALLOCATE(vnew)

!  Subtract the mean of the residual
!  Must take into account fixed occupation number in case of spin-polarized
   do ispden=1,dtset%nspden
     if (dtset%nspden==2.and.dtset%occopt>=3.and. abs(dtset%spinmagntarget+99.99_dp)<1.0d-10)then
       vres_mean(ispden)=(vmean(1)+vmean(2))*half
     else
       vres_mean(ispden)=vmean(ispden)
     end if
   end do

!$OMP PARALLEL DO COLLAPSE(2)
   do ispden=1,dtset%nspden
     do ifft=1,nfft
       vresidnew(ifft,ispden)=vresidnew(ifft,ispden)-vres_mean(ispden)
     end do
   end do

!  Compute square norm vres2 of potential residual vresid
   call sqnorm_v(1,nfft,vres2,dtset%nspden,optres,vresidnew(1+offset, 1),mpi_comm_sphgrid=mpi_comm_sphgrid)

!  Now take care of Vxctau residual (metaGGA)
   if (with_vxctau) then
     if (ipositron/=1) vtauresid(:,:)=vxctau(:,:,1)-vtauresid(:,:)
     if (ipositron==1) vtauresid(:,:)=zero
   end if

 else ! optres/=0

!  ------Produce new value of trial potential-------------

   if (.not. wvlbigdft) then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=1,min(dtset%nspden,2)
       do ifft=1,nfft
         vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)+vzeeman(ispden)+v_constr_dft_r(ifft,ispden)
         !vtrial(ifft,ispden)=vtrial(ifft,ispden)+vzeemanHarm(ifft,ispden)
         if(mod(dtset%fockoptmix,100)==11)vtrial(ifft,ispden)=vtrial(ifft,ispden)+vxc_hybcomp(ifft,ispden)
       end do
     end do
     if(dtset%nspden==4) then
!$OMP PARALLEL DO
       do ifft=1,nfft
         vtrial(ifft,3:4)=vxc(ifft,3:4)+vzeeman(3:4)+v_constr_dft_r(ifft,3:4)
         !vtrial(ifft,3:4)=vtrial(ifft,3:4)+vzeemanHarm(ifft,3:4)
         if(mod(dtset%fockoptmix,100)==11)vtrial(ifft,3:4)=vtrial(ifft,3:4)+vxc_hybcomp(ifft,3:4)
       end do
     end if
!    Pass vtrial to BigDFT object
     if(dtset%usewvl==1) then
       call wvl_vtrial_abi2big(1,vtrial,wvl%den)
     end if
   else
!    Compute with covering comms the different part of the potential.
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     call wvl_psitohpsi(dtset%diemix,energies%e_exactX, energies%e_xc, energies%e_hartree, &
&     energies%e_kinetic, energies%e_localpsp, energies%e_nlpsp_vfock, energies%e_sicdc, &
&     istep + 1, 1, dtset%iscf, mpi_enreg%me_wvl, &
&     dtset%natom, dtset%nfft, mpi_enreg%nproc_wvl,&
&     dtset%nspden,vres2, .true.,energies%e_xcdc,  wvl,&
&     wvlbigdft, xcart, strsxc, vtrial, vxc)
     ABI_DEALLOCATE(xcart)
!    Compute vxcavg
     call mean_fftr(vxc, vmean(1:1), nfft, nfftot, dtset%nspden,&
&     mpi_comm_sphgrid=mpi_comm_sphgrid)
     vxcavg = vmean(1)
   end if

 end if

 ABI_DEALLOCATE(v_constr_dft_r)
 !ABI_DEALLOCATE(vzeemanHarm) !SPr: debug for q/=0 magnetic field

 call timab(945,2,tsec)
 call timab(940,2,tsec)

 DBG_EXIT("COLL")

end subroutine rhotov
!!***

end module m_rhotov
!!***
