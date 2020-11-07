!!****m* ABINIT/m_newrho
!! NAME
!!  m_newrho
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2020 ABINIT group (MT).
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

module m_newrho

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_ab7_mixing
 use m_abi2big
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes,     only : MPI_type
 use m_time,     only : timab
 use m_geometry, only : metric
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type,pawrhoij_filter
 use m_prcref,   only : prcref
 use m_wvl_rho, only : wvl_prcref
 use m_fft,     only : fourdp

 implicit none

 private
!!***

 public :: newrho
!!***

contains
!!***

!!****f* ABINIT/newrho
!! NAME
!! newrho
!!
!! FUNCTION
!! Compute new trial density by mixing new and old values.
!! Call prcref to compute preconditioned residual density and forces,
!! Then, call one of the self-consistency drivers, then update density.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | densfor_pred= governs the preconditioning of the atomic charges
!!   | iprcel= governs the preconditioning of the density residual
!!   | iprcfc= governs the preconditioning of the forces
!!   | iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   |  iscf =11 => determination of the largest eigenvalue of the SCF cycle
!!   |  iscf =12 => SCF cycle, simple mixing
!!   |  iscf =13 => SCF cycle, Anderson mixing
!!   |  iscf =14 => SCF cycle, Anderson mixing (order 2)
!!   |  iscf =15 => SCF cycle, CG based on the minimization of the energy
!!   |  iscf =17 => SCF cycle, Pulay mixing
!!   | isecur=level of security of the computation
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | pawoptmix=-PAW- 1 if the computed residuals include the PAW (rhoij) part
!!   | prtvol=control print volume and debugging
!!  etotal=the total energy obtained from the input density
!!  fnametmp_fft=name of _FFT file
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only for electronic
!!     dielectric matrix
!!  mgfft=maximum size of 1D FFTs
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned density residual must be computed;
!!    otherwise, compute only the preconditioned density residual.
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
!!  nresid(nfft,nspden)=array for the residual of the density
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vtrial(nfft,nspden)=the trial potential that gave vresid.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  tauresid(nfft,nspden*dtset%usekden)=array for kinetic energy density residue (out - in)
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  mix<type(ab7_mixing_object)>=all data defining the mixing algorithm for the density
!!  rhor(nfft,nspden)= at input, it is the "out" trial density that gave nresid=(rho_out-rho_in)
!!                     at output, it is an updated "mixed" trial density
!!  rhog(2,nfft)= Fourier transform of the new trial density
!!  ===== if usekden==1 =====
!!  [mix_mgga<type(ab7_mixing_object)>]=all data defining the mixing algorithm
!!     for the kinetic energy density
!!  ===== if densfor_pred==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ==== if usepaw==1
!!    pawrhoij(natom)%nrhoijsel=number of non-zero values of rhoij
!!    pawrhoij(natom)%rhoijp(cplex_rhoij*lmn2_size,nspden)= new (mixed) value of rhoij quantities in PACKED STORAGE
!!    pawrhoij(natom)%rhoijselect(lmn2_size)=select the non-zero values of rhoij
!!  taug(2,nfft*dtset%usekden)=array for Fourier transform of kinetic
!!     energy density
!!  taur(nfft,nspden*dtset%usekden)=array for kinetic energy density
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid (except f_fftgr).
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!      ab7_mixing_copy_current_step,ab7_mixing_eval,ab7_mixing_eval_allocate
!!      ab7_mixing_eval_deallocate,ab7_mixing_use_moving_atoms,fourdp,metric
!!      pawrhoij_filter,prcref,timab,wvl_prcref,wvl_rho_abi2big
!!
!! SOURCE

subroutine newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,&
&  gmet,grhf,gsqcut,initialized,ispmix,istep,kg_diel,kxc,mgfft,mix,mixtofft,&
&  moved_atm_inside,mpi_enreg,my_natom,nattyp,nfft,&
&  nfftmix,nfftmix_per_nfft,ngfft,ngfftmix,nkxc,npawmix,npwdiel,&
&  nresid,ntypat,n1xccc,pawrhoij,pawtab,&
&  ph1d,psps,rhog,rhor,rprimd,susmat,usepaw,vtrial,wvl,wvl_den,xred,&
&  mix_mgga,taug,taur,tauresid)

!Arguments-------------------------------
!scalars
 integer,intent(in) :: dielstrt,initialized,ispmix,istep,my_natom,mgfft
 integer,intent(in) :: moved_atm_inside,n1xccc,nfft
 integer,intent(in) :: nfftmix,nfftmix_per_nfft
 integer,intent(in) :: nkxc,npawmix,npwdiel,ntypat,usepaw
 integer,intent(inout) :: dbl_nnsclo
 real(dp),intent(in) :: etotal,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(ab7_mixing_object), intent(inout) :: mix
 type(ab7_mixing_object), intent(inout),optional :: mix_mgga
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(in) :: ffttomix(nfft*(nfftmix_per_nfft))
 integer,intent(in) :: kg_diel(3,npwdiel)
 integer,intent(in) :: mixtofft(nfftmix*nfftmix_per_nfft)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),ngfftmix(18)
 real(dp),intent(in) :: dielar(7),fcart(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(inout) :: rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in), target :: vtrial(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout), target :: dtn_pc(3,dtset%natom)
 real(dp),intent(inout) :: gmet(3,3)
!TODO: nresid appears to be only intent in here.
 real(dp),intent(inout) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden)
 real(dp),intent(inout), target :: xred(3,dtset%natom)
 real(dp),intent(inout) :: rhog(2,nfft)
 real(dp),intent(inout), optional :: taug(2,nfft*dtset%usekden)
 real(dp),intent(inout), optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(inout), optional :: tauresid(nfft,dtset%nspden*dtset%usekden)

 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourdp9=9
 integer :: cplex,dplex,errid,i_vresid1,i_vrespc1,iatom,ifft,indx,iq,iq0,irhoij,ispden,jfft
 integer :: jrhoij,klmn,kklmn,kmix,mpicomm,nfftot,qphase
 logical :: mpi_summarize,reset
 real(dp) :: fact,ucvol,ucvol_local
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3),rmet(3,3),ro(2),tsec(2),vhartr_dum(1),vpsp_dum(1)
 real(dp) :: vxc_dum(1,1)
 real(dp),allocatable :: magng(:,:,:),magntaug(:,:,:)
 real(dp),allocatable :: nresid0(:,:),nrespc(:,:),nreswk(:,:,:)
 real(dp),allocatable :: rhoijrespc(:),rhoijtmp(:,:)
! TODO : these should be allocatables not pointers: is there some reason to
!  keep them this way, eg an interface somewhere?
 real(dp), pointer :: rhomag(:,:), npaw(:)
 real(dp),allocatable :: tauresid0(:,:),taurespc(:,:)
 real(dp),allocatable :: taumag(:,:)

! *************************************************************************

 DBG_ENTER("COLL")
 call timab(94,1,tsec)

 nfftot=PRODUCT(ngfft(1:3))

!Compatibility tests
 if(nfftmix>nfft) then
   message='nfftmix>nfft not allowed!'
   MSG_BUG(message)
 end if

 if(dtset%usewvl==1) then
   if( (ispmix/=1 .or. nfftmix/=nfft)) then
     message='nfftmix/=nfft, ispmix/=1 not allowed for wavelets!'
     MSG_BUG(message)
   end if
   if(dtset%wvl_bigdft_comp==1) then
     message='usewvl == 1 and wvl_bigdft_comp==1 not allowed!'
     MSG_BUG(message)
   end if
 end if

 if(ispmix/=2.and.nfftmix/=nfft) then
   message='nfftmix/=nfft allowed only when ispmix=2!'
   MSG_BUG(message)
 end if

 if (dtset%usekden==1) then
   if ((.not.present(tauresid)).or.(.not.present(taug)).or. &
&      (.not.present(taur)).or.(.not.present(mix_mgga))) then
     message='Several arrays are missing!'
     MSG_BUG(message)
   end if
   if (mix_mgga%iscf==AB7_MIXING_CG_ENERGY.or.mix_mgga%iscf==AB7_MIXING_CG_ENERGY_2.or.&
&      mix_mgga%iscf==AB7_MIXING_EIG) then
     message='kinetic energy density cannot be mixed with the selected mixing algorithm!'
     MSG_ERROR(message)
   end if
 end if

 if (usepaw==1.and.my_natom>0) then
   cplex=pawrhoij(1)%cplex_rhoij;dplex=cplex-1
   qphase=pawrhoij(1)%qphase
 else
   cplex = 0;dplex = 0 ; qphase=0
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if(dtset%usewvl==0) then
   ucvol_local=ucvol
#if defined HAVE_BIGDFT
 else
   ucvol_local = product(wvl_den%denspot%dpbox%hgrids) * real(nfftot, dp)
#endif
 end if

!Select components of density to be mixed
 ABI_ALLOCATE(rhomag,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(nresid0,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(taumag,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 ABI_ALLOCATE(tauresid0,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 ! real space and all fft points are here
 if (ispmix==1.and.nfft==nfftmix) then
   rhomag(:,1:dtset%nspden)=rhor(:,1:dtset%nspden)
   nresid0(:,1:dtset%nspden)=nresid(:,1:dtset%nspden)
   if (dtset%usekden==1) then
     taumag(:,1:dtset%nspden)=taur(:,1:dtset%nspden)
     tauresid0(:,1:dtset%nspden)=tauresid(:,1:dtset%nspden)
   end if
 ! recip space and all fft points are here
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(1,nresid0(:,ispden),nresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
   end do
   rhomag(:,1)=reshape(rhog,(/2*nfft/))
   if (dtset%nspden>1) then
     do ispden=2,dtset%nspden
       call fourdp(1,rhomag(:,ispden),rhor(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end do
   end if
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       call fourdp(1,tauresid0(:,ispden),tauresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end do
     taumag(:,1)=reshape(taug,(/2*nfft/))
     if (dtset%nspden>1) then
       do ispden=2,dtset%nspden
         call fourdp(1,taumag(:,ispden),taur(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
       end do
     end if
   end if
 ! not all fft points are here - presumes recip space
 else
   fact=dielar(4)-1._dp
   ABI_ALLOCATE(nreswk,(2,nfft,dtset%nspden))
   do ispden=1,dtset%nspden
     call fourdp(1,nreswk(:,:,ispden),nresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
   end do
   do ifft=1,nfft
     if (ffttomix(ifft)>0) then
       jfft=2*ffttomix(ifft)
       rhomag (jfft-1:jfft,1)=rhog(1:2,ifft)
       nresid0(jfft-1:jfft,1)=nreswk(1:2,ifft,1)
     else
       rhog(:,ifft)=rhog(:,ifft)+fact*nreswk(:,ifft,1)
     end if
   end do
   if (dtset%nspden>1) then
     ABI_ALLOCATE(magng,(2,nfft,dtset%nspden-1))
     do ispden=2,dtset%nspden
       call fourdp(1,magng(:,:,ispden-1),rhor(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
       do ifft=1,nfft
         if (ffttomix(ifft)>0) then
           jfft=2*ffttomix(ifft)
           rhomag (jfft-1:jfft,ispden)=magng (1:2,ifft,ispden-1)
           nresid0(jfft-1:jfft,ispden)=nreswk(1:2,ifft,ispden)
         else
           magng(:,ifft,ispden-1)=magng(:,ifft,ispden-1)+fact*nreswk(:,ifft,ispden)
           if (dtset%nspden==2) magng(:,ifft,1)=two*magng(:,ifft,1)-rhog(:,ifft)
         end if
       end do
     end do
   end if
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       call fourdp(1,nreswk(:,:,ispden),tauresid(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end do
     do ifft=1,nfft
       if (ffttomix(ifft)>0) then
         jfft=2*ffttomix(ifft)
         taumag (jfft-1:jfft,1)=taug(1:2,ifft)
         tauresid0(jfft-1:jfft,1)=nreswk(1:2,ifft,1)
       else
         taug(:,ifft)=taug(:,ifft)+fact*nreswk(:,ifft,1)
       end if
     end do
     if (dtset%nspden>1) then
       ABI_ALLOCATE(magntaug,(2,nfft,dtset%nspden-1))
       do ispden=2,dtset%nspden
         call fourdp(1,magntaug(:,:,ispden-1),taur(:,ispden),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
         do ifft=1,nfft
           if (ffttomix(ifft)>0) then
             jfft=2*ffttomix(ifft)
             taumag (jfft-1:jfft,ispden)=magntaug(1:2,ifft,ispden-1)
             tauresid0(jfft-1:jfft,ispden)=nreswk(1:2,ifft,ispden)
           else
             magntaug(:,ifft,ispden-1)=magntaug(:,ifft,ispden-1)+fact*nreswk(:,ifft,ispden)
             if (dtset%nspden==2) magntaug(:,ifft,1)=two*magntaug(:,ifft,1)-taug(:,ifft)
           end if
         end do
       end do
     end if
   end if
   ABI_DEALLOCATE(nreswk)
 end if

!Retrieve "input" density from "output" density and density residual
 rhomag(:,1:dtset%nspden)=rhomag(:,1:dtset%nspden)-nresid0(:,1:dtset%nspden)
 if (dtset%usekden==1) then
   taumag(:,1:dtset%nspden)=taumag(:,1:dtset%nspden)-tauresid0(:,1:dtset%nspden)
 end if

!If nspden==2, separate density and magnetization
 if (dtset%nspden==2) then
   rhomag (:,2)=two*rhomag (:,2)-rhomag (:,1)
   nresid0(:,2)=two*nresid0(:,2)-nresid0(:,1)
   if (dtset%usekden==1) then
     taumag (:,2)=two*taumag (:,2)-taumag (:,1)
     tauresid0(:,2)=two*tauresid0(:,2)-tauresid0(:,1)
   end if
 end if

!If PAW, handle occupancy matrix
 if (usepaw==1.and.my_natom>0) then
   if (pawrhoij(1)%nspden==2) then
     do iatom=1,my_natom
       do iq=1,qphase
         iq0=merge(0,cplex*pawrhoij(iatom)%lmn2_size,iq==1)
         jrhoij=1+iq0
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           ro(1:1+dplex)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)
           pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)=ro(1:1+dplex)+pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)
           pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)=ro(1:1+dplex)-pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)
           jrhoij=jrhoij+cplex
         end do
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           klmn=iq0+cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex
           ro(1:1+dplex)=pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,1)
           pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,1)=ro(1:1+dplex)+pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,2)
           pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,2)=ro(1:1+dplex)-pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,2)
         end do
       end do
     end do
   end if
 end if

!Choice of preconditioner governed by iprcel, densfor_pred and iprcfc
 ABI_ALLOCATE(nrespc,(ispmix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(taurespc,(ispmix*nfftmix,dtset%nspden*dtset%usekden))
 ABI_ALLOCATE(npaw,(npawmix*usepaw))
 if (usepaw==1)  then
   ABI_ALLOCATE(rhoijrespc,(npawmix))
 else
   ABI_ALLOCATE(rhoijrespc,(0))
 end if
 if(dtset%usewvl==0) then
   call prcref(atindx,dielar,dielinv,&
&   dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&
&   istep,kg_diel,kxc,&
&   mgfft,moved_atm_inside,mpi_enreg,my_natom,&
&   nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&   ispmix,1,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&
&   susmat,vhartr_dum,vpsp_dum,nresid0,nrespc,vxc_dum,wvl,wvl_den,xred)
 else
   call wvl_prcref(dielar,dtset%iprcel,my_natom,nfftmix,npawmix,dtset%nspden,pawrhoij,&
&   rhoijrespc,psps%usepaw,nresid0,nrespc)
 end if
!At present, only a simple precoditionning for the kinetic energy density
! (is Kerker mixing valid for tau?)
 if (dtset%usekden==1) then
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=abs(dielar(7))
     taurespc(1:ispmix*nfftmix,ispden)=fact*tauresid0(1:ispmix*nfftmix,ispden)
   end do
 end if

!------Compute new trial density and eventual new atomic positions

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
   call ab7_mixing_copy_current_step(mix, nresid0, errid, message, arr_respc = nrespc, arr_atm = grhf)
 else
   call ab7_mixing_copy_current_step(mix, nresid0, errid, message, arr_respc = nrespc)
 end if
 if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(message)
 end if

!Same treatment for the kinetic energy density
 if (dtset%usekden==1) then
   call ab7_mixing_eval_allocate(mix_mgga, istep)
   call ab7_mixing_copy_current_step(mix_mgga, tauresid0, errid, message, arr_respc = taurespc)
   if (errid /= AB7_NO_ERROR) then
     MSG_ERROR(message)
   end if
 end if

 ABI_DEALLOCATE(nresid0)
 ABI_DEALLOCATE(nrespc)
 ABI_DEALLOCATE(tauresid0)
 ABI_DEALLOCATE(taurespc)

!PAW: either use the array f_paw or the array f_paw_disk
 if (usepaw==1) then
   indx=-dplex
   do iatom=1,my_natom
     ABI_ALLOCATE(rhoijtmp,(cplex*pawrhoij(iatom)%lmn2_size,1))
     do iq=1,qphase
       iq0=merge(0,cplex*pawrhoij(iatom)%lmn2_size,iq==1)
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero ; jrhoij=1+iq0
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,1)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex
         end do
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+cplex;klmn=cplex*pawrhoij(iatom)%kpawmix(kmix)-dplex ; kklmn=klmn+iq0
           npaw(indx:indx+dplex)=rhoijtmp(klmn:klmn+dplex,1)-pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vresid1)=pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vrespc1)=rhoijrespc(indx:indx+dplex)
         end do
       end do
     end do
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if

!------Prediction of the components of the density

!Init mpicomm
 if(mpi_enreg%paral_kgb==1)then
   mpicomm=mpi_enreg%comm_fft
   mpi_summarize=.true.
 else
   mpicomm=0
   mpi_summarize=.false.
 end if
 if(dtset%usewvl==1) then
   mpicomm=mpi_enreg%comm_wvl
   mpi_summarize=(mpi_enreg%nproc_wvl > 1)
 end if

 reset = .false.
 if (initialized == 0) reset = .true.

!Electronic density mixing
 call ab7_mixing_eval(mix, rhomag, istep, nfftot, ucvol_local, &
& mpicomm, mpi_summarize, errid, message, &
& reset = reset, isecur = dtset%isecur,&
& pawopt = dtset%pawoptmix, pawarr = npaw, &
& etotal = etotal, potden = vtrial, &
& comm_atom=mpi_enreg%comm_atom)
 if (errid == AB7_ERROR_MIXING_INC_NNSLOOP) then
   dbl_nnsclo = 1
 else if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(message)
 end if
!Kinetic energy density mixing (if any)
 if (dtset%usekden==1) then
   call ab7_mixing_eval(mix_mgga, taumag, istep, nfftot, ucvol_local, &
&   mpicomm, mpi_summarize, errid, message, reset = reset)
   if (errid /= AB7_NO_ERROR) then
     MSG_ERROR(message)
   end if
 end if

!PAW: apply a simple mixing to rhoij (this is temporary)
 if(dtset%iscf==15 .or. dtset%iscf==16)then
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
         if (pawrhoij(iatom)%nspden/=2) then
           do ispden=1,pawrhoij(iatom)%nspden
             jrhoij=iq0+1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=iq0+cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
               rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&               +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
               jrhoij=jrhoij+cplex
             end do
           end do
         else
           jrhoij=iq0+1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=iq0+cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
             ro(1:1+dplex)=rhoijtmp(klmn:klmn+dplex,1)
             rhoijtmp(klmn:klmn+dplex,1)=half*(ro(1:1+dplex)+rhoijtmp(klmn:klmn+dplex,2)) &
&             +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,1)
             rhoijtmp(klmn:klmn+dplex,2)=half*(ro(1:1+dplex)-rhoijtmp(klmn:klmn+dplex,2)) &
&             +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,2)
             jrhoij=jrhoij+cplex
           end do
         end if
       end do
       call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&           pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&           pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if
 end if

 !if (usepaw==1)  then
 ABI_DEALLOCATE(rhoijrespc)
 !end if

!PAW: restore rhoij from compact storage
 if (usepaw==1.and.dtset%iscf/=15.and.dtset%iscf/=16) then
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
           rhoijtmp(klmn:klmn+dplex,ispden)=npaw(indx:indx+dplex)
         end do
       end do
       if (pawrhoij(iatom)%nspden==2) then
         jrhoij=iq0+1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=iq0+cplex*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           ro(1:1+dplex)=rhoijtmp(klmn:klmn+dplex,1)
           rhoijtmp(klmn:klmn+dplex,1)=half*(ro(1:1+dplex)+rhoijtmp(klmn:klmn+dplex,2))
           rhoijtmp(klmn:klmn+dplex,2)=half*(ro(1:1+dplex)-rhoijtmp(klmn:klmn+dplex,2))
           jrhoij=jrhoij+cplex
         end do
       end if
     end do
     call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&         pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&         pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if   ! usepaw==1.and.dtset%iscf/=15.and.dtset%iscf/=16
 ABI_DEALLOCATE(npaw)

!Eventually write the data on disk and deallocate f_fftgr_disk
 call ab7_mixing_eval_deallocate(mix)
 if (dtset%usekden==1) call ab7_mixing_eval_deallocate(mix_mgga)

!Fourier transform the density
 if (ispmix==1.and.nfft==nfftmix) then
   !Real space mixing, no need to transform rhomag
   rhor(:,1:dtset%nspden)=rhomag(:,1:dtset%nspden)
   if(dtset%usewvl==0) then
     !Get rhog from rhor(:,1)
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
   end if
   if (dtset%usekden==1) then
     taur(:,1:dtset%nspden)=taumag(:,1:dtset%nspden)
     if(dtset%usewvl==0) then
       call fourdp(1,taug,taur(:,1),-1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end if
   end if
 else if (nfft==nfftmix) then
   !Reciprocal mixing space mixing, need to generate rhor in real space from rhomag in reciprocal space
   do ispden=1,dtset%nspden
     call fourdp(1,rhomag(:,ispden),rhor(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
   end do
   rhog(:,:)=reshape(rhomag(:,1),(/2,nfft/))
   if (dtset%usekden==1) then
     do ispden=1,dtset%nspden
       call fourdp(1,taumag(:,ispden),taur(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end do
     taug(:,:)=reshape(taumag(:,1),(/2,nfft/))
   end if
 else
   do ifft=1,nfftmix
     jfft=mixtofft(ifft)
     rhog(1:2,jfft)=rhomag(2*ifft-1:2*ifft,1)
   end do
   call fourdp(1,rhog,rhor(:,1),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
   if (dtset%nspden>1) then
     do ispden=2,dtset%nspden
       do ifft=1,nfftmix
         jfft=mixtofft(ifft)
         magng(1:2,jfft,ispden-1)=rhomag(2*ifft-1:2*ifft,ispden)
       end do
       call fourdp(1,magng(:,:,ispden-1),rhor(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     end do
     ABI_DEALLOCATE(magng)
   end if
   if (dtset%usekden==1) then
     do ifft=1,nfftmix
       jfft=mixtofft(ifft)
       taug(1:2,jfft)=taumag(2*ifft-1:2*ifft,1)
     end do
     call fourdp(1,taug,taur(:,1),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
     if (dtset%nspden>1) then
       do ispden=2,dtset%nspden
         do ifft=1,nfftmix
           jfft=mixtofft(ifft)
           magntaug(1:2,jfft,ispden-1)=taumag(2*ifft-1:2*ifft,ispden)
         end do
         call fourdp(1,magntaug(:,:,ispden-1),taur(:,ispden),+1,mpi_enreg,nfft,1,ngfft,tim_fourdp9)
       end do
       ABI_DEALLOCATE(magntaug)
     end if
   end if
 end if
 ABI_DEALLOCATE(rhomag)
 ABI_DEALLOCATE(taumag)

!Set back rho in (up+dn,up) form if nspden=2
 if (dtset%nspden==2) then
   rhor(:,2)=half*(rhor(:,1)+rhor(:,2))
   if (dtset%usekden==1) taur(:,2)=half*(taur(:,1)+taur(:,2))
 end if

!In WVL: copy density to BigDFT object:
 if(dtset%usewvl==1) then
   call wvl_rho_abi2big(1,rhor,wvl_den)
 end if

 call timab(94,2,tsec)

 DBG_EXIT("COLL")

end subroutine newrho
!!***

end module m_newrho
!!***
