!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_newvtr
!! NAME
!!  m_newvtr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, MT)
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

module m_newvtr

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_errors
 use m_abi2big
 use m_ab7_mixing
 use m_cgtools
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,     only : MPI_type
 use m_time,     only : timab
 use m_geometry, only : metric
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type,pawrhoij_filter
 use m_prcref,   only : prcref_PMA
 use m_wvl_rho,  only : wvl_prcref
 use m_fft,      only : fourdp
 use m_xctk,     only : xcpot

 implicit none

 private
!!***

 public :: newvtr
!!***

contains
!!***

!!****f* ABINIT/newvtr
!! NAME
!! newvtr
!!
!! FUNCTION
!! Compute new trial potential by mixing new and old values.
!! Call prcref to compute preconditioned residual potential and forces,
!! Then, call one of the self-consistency drivers,
!! then update vtrial.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | spinmagntarget=input variable that governs fixed moment calculation
!!   | intxc=control xc quadrature
!!   | densfor_pred= governs the preconditioning of the atomic charges
!!   | iprcel= governs the preconditioning of the potential residual
!!   | iprcfc=governs the preconditioning of the forces
!!   | iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   |  iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!   |  iscf =2 => SCF cycle, simple mixing
!!   |  iscf =3 => SCF cycle, Anderson mixing
!!   |  iscf =4 => SCF cycle, Anderson mixing (order 2)
!!   |  iscf =5 => SCF cycle, CG based on the minimization of the energy
!!   |  iscf =7 => SCF cycle, Pulay mixing
!!   | isecur=level of security of the computation
!!   | ixc=exchange-correlation choice parameter.
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | paral_kgb=option for (kpt,g vectors,bands) parallelism
!!   | pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  etotal=the total energy obtained from the input vtrial
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only for electronic!
!     dielectric matrix
!!  mgfft=maximum size of 1D FFTs
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  nstep=number of steps expected in iterations.
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  [vtauresid(nfft,nspden*dtset%usekden)]=array for vtau residue (see vtau below))
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vnew_mean(nspden)=constrained mean value of the future trial potential (might be
!!    spin-polarized
!!  vpsp(nfft)=array for holding local psp
!!  vresid(nfft,nspden)=array for the residual of the potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  [vtau(nfftf,dtset%nspden,4*dtset%usekden)]=derivative of XC energy density
!!      with respect to kinetic energy density (metaGGA cases) (optional)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  vtrial(nfft,nspden)= at input, it is the "in" trial potential that gave vresid=(v_out-v_in)
!!       at output, it is an updated "mixed" trial potential
!!  ===== if usekden==1 =====
!!  [mix_mgga<type(ab7_mixing_object)>]=all data defining the mixing algorithm for
!!    the kinetic energy potential
!!  ===== if densfor_pred==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ==== if usepaw==1
!!    pawrhoij(natom)%nrhoijsel,rhoijselect,rhoijp= several arrays
!!                containing new values of rhoij (augmentation occupancies)
!!
!! WARNINGS
!! depending on the value of densfor_pred and moved_atm_inside,
!! the xc potential or the Hxc potential may have been subtracted from vtrial !
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developpers have to be careful when introducing others arrays:
!!    they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!!  Subtility in PAW and non-collinear magnetism:
!!    Potentials are stored in (up-up,dn-dn,Re[up-dn],Im[up-dn]) format
!!    On-site occupancies (rhoij) are stored in (n,mx,my,mz)
!!    This is compatible provided that the mixing factors for n and m are identical
!!    and that the residual is not a combination of V_res and rhoij_res (pawoptmix=0).
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ab7_mixing_copy_current_step,ab7_mixing_eval,ab7_mixing_eval_allocate
!!      ab7_mixing_eval_deallocate,ab7_mixing_use_moving_atoms,fourdp,mean_fftr
!!      metric,prcref_pma,timab,wvl_prcref,wvl_vtrial_abi2big
!!
!! SOURCE

subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
     &  dtn_pc,dtset,etotal,fcart,ffttomix,&
     &  gmet,grhf,gsqcut,&
     &  initialized,ispmix,&
     &  istep,&
     &  kg_diel,kxc,mgfft,mix,mixtofft,&
     &  moved_atm_inside,mpi_enreg,my_natom,nattyp,nfft,nfftmix,&
     &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&
     &  nstep,ntypat,n1xccc,&
     &  pawrhoij,&
     &  ph1d,&
     &  psps,rhor,rprimd,susmat,usepaw,&
     &  vhartr,vnew_mean,vpsp,vresid,&
     &  vtrial,vxc,xred,&
     &  nfftf,&
     &  pawtab,&
     &  rhog,&
     &  wvl,&
     &  mix_mgga,vtau,vtauresid) ! Optional arguments

!Arguments-------------------------------
  ! WARNING
  ! BEWARE THERE IS TWO DIFFERENT SIZE DECLARED FOR ARRAY NHAT IN RHOTOV AND RHOHXC
  ! THIS MIGHT RESULT IN A BUG
!scalars
 integer,intent(in) :: dielstrt,initialized,ispmix,istep,mgfft
 integer,intent(in) :: moved_atm_inside,my_natom,n1xccc,nfft
 integer,intent(in) :: nfftf,nfftmix,nkxc,npawmix,npwdiel,nstep
 integer,intent(in) :: ntypat,usepaw
 integer,intent(inout) :: dbl_nnsclo
 real(dp),intent(in) :: etotal,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(ab7_mixing_object),intent(inout) :: mix
 type(ab7_mixing_object),intent(inout),optional :: mix_mgga
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data), intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: kg_diel(3,npwdiel)
 integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft)),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),ngfftmix(18)
 real(dp),intent(in) :: dielar(7)
 real(dp),intent(in) :: fcart(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vnew_mean(dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
 real(dp),intent(inout) :: gmet(3,3)
 real(dp),intent(inout) :: kxc(nfft,nkxc),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhog(2,nfftf),vpsp(nfft)
 real(dp),intent(inout), target :: rhor(nfft,dtset%nspden)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),vtrial(nfft,dtset%nspden)
 real(dp),intent(inout), target :: xred(3,dtset%natom)
 real(dp),intent(inout),optional :: vtau(nfft,dtset%nspden,4*dtset%usekden)
 real(dp),intent(inout),optional :: vtauresid(nfft,dtset%nspden*dtset%usekden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,dplex,i_vresid1,i_vrespc1
! integer :: i1,i2,i3,ifft2,ifft3,ifft4,ifft5,ii1,ii2,ii3,ii4,ii5
 integer :: errid,iatom,ifft,indx,iq,iq0,irhoij,ispden,jfft,jrhoij,klmn,kklmn,kmix
 integer :: mpicomm,mpi_comm_sphgrid,n1,n2,n3,nfftot,qphase,tim_fourdp
 logical :: mpi_summarize,reset
 real(dp) :: dielng,diemix,fact,ucvol,ucvol_local,vme
 character(len=500) :: message
!arrays
 real(dp),parameter :: identity(4)=(/one,one,zero,zero/)
 real(dp) :: gprimd(3,3),rmet(3,3),tsec(2),vmean(dtset%nspden)
 real(dp),allocatable :: rhoijrespc(:)
 real(dp),allocatable :: rhoijtmp(:,:)
 real(dp),allocatable :: vresid0(:,:),vrespc(:,:),vreswk(:,:),vtrialg(:,:,:)
 real(dp),allocatable :: vtauresid0(:,:),vtaurespc(:,:),vtaug(:,:,:),vtau0(:,:)
 real(dp),pointer :: vtrial0(:,:),vpaw(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' newvtr : enter '
!write(std_out,*)' newvtr : ispmix,nfft,nfftmix=',ispmix,nfft,nfftmix
!ENDDEBUG

 call timab(93,1,tsec)
 call timab(901,1,tsec)
 tim_fourdp=8

!mpicomm over spherical grid:
 mpi_comm_sphgrid=mpi_enreg%comm_fft
 if(dtset%usewvl==1) mpi_comm_sphgrid=mpi_enreg%comm_wvl

!Compatibility tests
 if(nfftmix>nfft) then
   MSG_BUG('  nfftmix>nfft not allowed !')
 end if

 if(ispmix/=2.and.nfftmix/=nfft) then
   message = '  nfftmix/=nfft allowed only when ispmix=2 !'
   MSG_BUG(message)
 end if

 if (dtset%usekden==1) then
   if ((.not.present(vtauresid)).or.(.not.present(vtau)).or.(.not.present(mix_mgga))) then
      message='Several arrays are mising!'
      MSG_BUG(message)
   end if
   if (mix_mgga%iscf==AB7_MIXING_CG_ENERGY.or.mix_mgga%iscf==AB7_MIXING_CG_ENERGY_2.or.&
&      mix_mgga%iscf==AB7_MIXING_EIG) then
     message='kinetic energy potential cannot be mixed with the selected mixing algorithm!'
     MSG_ERROR(message)
   end if
 end if

 if(dtset%usewvl==1) then
   if(dtset%wvl_bigdft_comp==1) then
     message = 'newvtr: usewvl == 1 and wvl_bigdft_comp==1 not allowed (use wvl_newtr() instead)!'
     MSG_BUG(message)
   end if
   if(ispmix/=1 .or. nfftmix/=nfft) then
     MSG_BUG('newvtr: nfftmix/=nfft, ispmix/=1 not allowed for wavelets')
   end if
 end if

 if(usepaw==1.and.dtset%nspden==4.and.dtset%pawoptmix==1) then
   message = ' pawoptmix=1 is not compatible with nspden=4 !'
   MSG_ERROR(message)
 end if

 dielng=dielar(2)
 diemix=dielar(4)
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 if (usepaw==1.and.my_natom>0) then
   cplex=pawrhoij(1)%cplex_rhoij;dplex=cplex-1
   qphase=pawrhoij(1)%qphase
 else
   cplex=0;dplex=0 ; qphase=0
 end if

!Get size of FFT grid
 nfftot=PRODUCT(ngfft(1:3))

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if(dtset%usewvl==0) then
   ucvol_local=ucvol
#if defined HAVE_BIGDFT
 else
   ucvol_local = product(wvl%den%denspot%dpbox%hgrids) * real(nfftot, dp)
#endif
 end if

!------Treat the mean of potentiel residual

!Special care must be taken with components of the
!potential that are associated with NO density change.
!In general, only the global mean of the potential has
!such an anomalous feature. However, in the spin
!polarized case with fixed occupancies, also the
!mean of each spin-potential (independently of the other)
!has such a behaviour. The trick is to remove these
!variables before going in the predictive routines,
!then to put them back

!Compute the mean of the old vtrial
 call mean_fftr(vtrial,vmean,nfft,nfftot,dtset%nspden,mpi_comm_sphgrid)

!When (collinear) spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%spinmagntarget+99.99_dp)<1.0d-10)then
     vme=(vmean(1)+vmean(2))*half
   else
     vme=vmean(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)-vme
 end do

 call timab(901,2,tsec)

!Select components of potential to be mixed
 ABI_ALLOCATE(vtrial0,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vresid0,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vtau0,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 ABI_ALLOCATE(vtauresid0,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial0=vtrial;vresid0=vresid
   if (dtset%usekden==1) then
     vtau0=vtau(:,:,1);vtauresid0=vtauresid
   end if
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     call fourdp(1,vresid0(:,ispden),vresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
   end do
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       call fourdp(1,vtau0(:,ispden),vtau(:,ispden,1),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
       call fourdp(1,vtauresid0(:,ispden),vtauresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     end do
   end if
 else
   ABI_ALLOCATE(vtrialg,(2,nfft,dtset%nspden))
   ABI_ALLOCATE(vreswk,(2,nfft))
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=dielar(7)
     call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     call fourdp(1,vreswk,vresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     do ifft=1,nfft
       if (ffttomix(ifft)>0) then
         jfft=2*ffttomix(ifft)
         vtrial0(jfft-1,ispden)=vtrialg(1,ifft,ispden)
         vtrial0(jfft  ,ispden)=vtrialg(2,ifft,ispden)
         vresid0(jfft-1,ispden)=vreswk(1,ifft)
         vresid0(jfft  ,ispden)=vreswk(2,ifft)
       else
         vtrialg(:,ifft,ispden)=vtrialg(:,ifft,ispden)+fact*vreswk(:,ifft)
       end if
     end do
   end do
   if (dtset%usekden==1) then
     ABI_ALLOCATE(vtaug,(2,nfft,dtset%nspden))
     do ispden=1,dtset%nspden
       fact=dielar(4);if (ispden>1) fact=dielar(7)
       call fourdp(1,vtaug(:,:,ispden),vtau(:,ispden,1),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
       call fourdp(1,vreswk,vtauresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
       do ifft=1,nfft
         if (ffttomix(ifft)>0) then
           jfft=2*ffttomix(ifft)
           vtau0(jfft-1,ispden)=vtaug(1,ifft,ispden)
           vtau0(jfft  ,ispden)=vtaug(2,ifft,ispden)
           vtauresid0(jfft-1,ispden)=vreswk(1,ifft)
           vtauresid0(jfft  ,ispden)=vreswk(2,ifft)
         else
           vtaug(:,ifft,ispden)=vtaug(:,ifft,ispden)+fact*vreswk(:,ifft)
         end if
       end do
     end do
   end if
   ABI_DEALLOCATE(vreswk)
 end if

 call timab(902,1,tsec)

!Choice of preconditioner governed by iprcel, densfor_pred and iprcfc
 ABI_ALLOCATE(vrespc,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vtaurespc,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 ABI_ALLOCATE(vpaw,(npawmix*usepaw))
 if (usepaw==1)  then
   ABI_ALLOCATE(rhoijrespc,(npawmix))
 else
   ABI_ALLOCATE(rhoijrespc,(0))
 end if

 call timab(902,2,tsec)
 call timab(903,1,tsec)

 if(dtset%usewvl==0) then
   call prcref_PMA(atindx,dielar,dielinv,dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
&   istep,kg_diel,kxc,mgfft,moved_atm_inside,mpi_enreg,my_natom,&
&   nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&   ispmix,0,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,susmat,&
&   vhartr,vpsp,vresid0,vrespc,vxc,xred,&
&   etotal,pawtab,wvl)
 else
   call wvl_prcref(dielar,dtset%iprcel,my_natom,nfftmix,npawmix,dtset%nspden,pawrhoij,&
&   rhoijrespc,psps%usepaw,vresid0,vrespc)
 end if
!At present, only a simple precoditionning for vtau
! (is Kerker mixing valid for vtau?)
 if (dtset%usekden==1) then
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=abs(dielar(7))
     vtaurespc(1:ispmix*nfftmix,ispden)=fact*vtauresid0(1:ispmix*nfftmix,ispden)
   end do
 end if

 call timab(903,2,tsec)
 call timab(904,1,tsec)

!------Compute new vtrial and eventual new atomic positions

 if (mix%n_fftgr>0) then
   i_vresid1=mix%i_vresid(1)
   i_vrespc1=mix%i_vrespc(1)
 end if

!Initialise working arrays for the mixing object.
 if (moved_atm_inside == 1) then
   call ab7_mixing_use_moving_atoms(mix, dtset%natom, xred, dtn_pc)
 end if
 call ab7_mixing_eval_allocate(mix, istep)
!Copy current step arrays.
 if (moved_atm_inside == 1) then
   call ab7_mixing_copy_current_step(mix, vresid0, errid, message, &
&   arr_respc = vrespc, arr_atm = grhf)
 else
   call ab7_mixing_copy_current_step(mix, vresid0, errid, message, &
&   arr_respc = vrespc)
 end if
 if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(message)
 end if
 if (dtset%usekden==1) then
   call ab7_mixing_eval_allocate(mix_mgga, istep)
   call ab7_mixing_copy_current_step(mix_mgga, vtauresid0, errid, message, &
&        arr_respc = vtaurespc)
   if (errid /= AB7_NO_ERROR) then
     MSG_ERROR(message)
   end if
 end if
 ABI_DEALLOCATE(vresid0)
 ABI_DEALLOCATE(vrespc)
 ABI_DEALLOCATE(vtauresid0)
 ABI_DEALLOCATE(vtaurespc)

!PAW: either use the array f_paw or the array f_paw_disk
 if (usepaw==1) then
   indx=-dplex
   do iatom=1,my_natom
     ABI_ALLOCATE(rhoijtmp,(cplex*pawrhoij(iatom)%lmn2_size,1))
     do iq=1,qphase
       iq0=merge(0,cplex*pawrhoij(iatom)%lmn2_size,iq==1)
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero ; jrhoij=iq0+1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,1)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex
         end do
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+cplex;klmn=cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex; kklmn=iq0+klmn
           vpaw(indx:indx+dplex)=rhoijtmp(klmn:klmn+dplex,1)-pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vresid1)=pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vrespc1)=rhoijrespc(indx:indx+dplex)
         end do
       end do
     end do
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if

!------Prediction of the components of the potential associated with a density change

!Init mpicomm
 if(mpi_enreg%paral_kgb==1)then
   mpicomm=mpi_enreg%comm_fft
   mpi_summarize=.true.
 else
   mpicomm=0
   mpi_summarize=.false.
 end if

 reset = .false.
 if (initialized == 0) reset = .true.
 call ab7_mixing_eval(mix, vtrial0, istep, nfftot, ucvol_local, &
& mpicomm, mpi_summarize, errid, message, &
& reset = reset, isecur = dtset%isecur, &
& pawopt = dtset%pawoptmix, pawarr = vpaw, etotal = etotal, potden = rhor, &
& comm_atom=mpi_enreg%comm_atom)
 if (errid == AB7_ERROR_MIXING_INC_NNSLOOP) then
   dbl_nnsclo = 1
 else if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(message)
 end if
 if (dtset%usekden==1) then
   call ab7_mixing_eval(mix_mgga, vtau0, istep, nfftot, ucvol_local, &
&   mpicomm, mpi_summarize, errid, message, reset = reset)
   if (errid /= AB7_NO_ERROR) then
     MSG_ERROR(message)
   end if
 end if

!PAW: apply a simple mixing to rhoij (this is temporary)
 if(dtset%iscf==5 .or. dtset%iscf==6)then
   if (usepaw==1) then
     indx=-dplex
     do iatom=1,my_natom
       ABI_ALLOCATE(rhoijtmp,(cplex*qphase*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       do iq=1,qphase
         iq0=merge(0,cplex*pawrhoij(iatom)%lmn2_size,iq==1)
         if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
           do ispden=1,pawrhoij(iatom)%nspden
             do kmix=1,pawrhoij(iatom)%lmnmix_sz
               indx=indx+cplex;klmn=iq0+cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
               rhoijtmp(klmn:klmn+dplex,ispden)=rhoijrespc(indx:indx+dplex) &
&               -pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
             end do
           end do
         end if
         do ispden=1,pawrhoij(iatom)%nspden
           jrhoij=iq0+1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=iq0+cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&             +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex
           end do
         end do
       end do
       call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&           pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&           pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if
 end if

 ABI_DEALLOCATE(rhoijrespc)

!PAW: restore rhoij from compact storage
 if (usepaw==1.and.dtset%iscf/=5.and.dtset%iscf/=6) then
   indx=-dplex
   do iatom=1,my_natom
     ABI_ALLOCATE(rhoijtmp,(cplex*qphase*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
     rhoijtmp=zero
     do iq=1,qphase
       iq0=merge(0,cplex*pawrhoij(iatom)%lmn2_size,iq==1)
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           jrhoij=iq0+1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=iq0+cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+cplex;klmn=iq0+cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
           rhoijtmp(klmn:klmn+dplex,ispden)=vpaw(indx:indx+dplex)
         end do
       end do
     end do
     call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&         pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&         pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if
 ABI_DEALLOCATE(vpaw)

!Eventually write the data on disk and deallocate f_fftgr_disk
 call ab7_mixing_eval_deallocate(mix)
 if (dtset%usekden==1) call ab7_mixing_eval_deallocate(mix_mgga)

 call timab(904,2,tsec)

!Restore potential
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial=vtrial0
   if (dtset%usekden==1) vtau(:,:,1)=vtau0(:,:)
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
   end do
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       call fourdp(1,vtau0(:,ispden),vtau(:,ispden,1),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     end do
   end if
 else
   do ispden=1,dtset%nspden
     do ifft=1,nfftmix
       jfft=mixtofft(ifft)
       vtrialg(1,jfft,ispden)=vtrial0(2*ifft-1,ispden)
       vtrialg(2,jfft,ispden)=vtrial0(2*ifft  ,ispden)
     end do
     call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
   end do
   ABI_DEALLOCATE(vtrialg)
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       do ifft=1,nfftmix
         jfft=mixtofft(ifft)
         vtaug(1,jfft,ispden)=vtau0(2*ifft-1,ispden)
         vtaug(2,jfft,ispden)=vtau0(2*ifft  ,ispden)
       end do
       call fourdp(1,vtaug(:,:,ispden),vtau(:,ispden,1),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp)
     end do
     ABI_DEALLOCATE(vtaug)
   end if
 end if
 ABI_DEALLOCATE(vtrial0)
 ABI_DEALLOCATE(vtau0)

!In case of metaGGA, re-compute vtau gradient
 if (dtset%usekden==1.and.mix_mgga%iscf/=AB7_MIXING_NONE) then
   call xcpot(1,gprimd,0,0,mpi_enreg,nfft,ngfft,2,dtset%nspden,0,[zero,zero,zero],vxctau=vtau)
 end if

 call timab(905,1,tsec)

!------Treat the mean of the potential

!Compute the mean of the new vtrial
 call mean_fftr(vtrial,vmean,nfft,nfftot,dtset%nspden,mpi_comm_sphgrid)

!Reset the mean of the new vtrial, to the value vnew_mean
!When spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%spinmagntarget+99.99_dp)<1.0d-10)then
     vme=(vnew_mean(1)+vnew_mean(2)-vmean(1)-vmean(2))*half
   else
     vme=vnew_mean(ispden)-vmean(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)+vme
 end do

 if(moved_atm_inside==1 .and. istep/=nstep )then
   if(abs(dtset%densfor_pred)==1.or.abs(dtset%densfor_pred)==4)then
!    Subtract current local psp, but also vxc (for core charges)
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vtrial(:,ispden)-vpsp(:)*identity(ispden)-vxc(:,ispden)
     end do
   else if(abs(dtset%densfor_pred)==2.or.abs(dtset%densfor_pred)==5.or.abs(dtset%densfor_pred)==6)then
!    Subtract current vpsp+Hxc from vtrial. This should be rationalized later
     do ispden=1,dtset%nspden
       vtrial(:,ispden)=vtrial(:,ispden)-(vpsp(:)+vhartr(:))*identity(ispden)-vxc(:,ispden)
     end do
   end if
 end if

!In WVL: copy vtrial to BigDFT object:
 if(dtset%usewvl==1) then
   call wvl_vtrial_abi2big(1,vtrial,wvl%den)
 end if

 call timab(905,2,tsec)
 call timab(93,2,tsec)

end subroutine newvtr
!!***

end module m_newvtr
!!***
