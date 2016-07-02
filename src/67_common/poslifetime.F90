!{\src2tex{textfont=tt}}
!!****f* ABINIT/poslifetime
!! NAME
!! poslifetime 
!!
!! FUNCTION
!! Calculate the positron lifetime
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (GJ,MT,JW)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | nspden=number of spin-density components
!!   | ntypat=number of atom types
!!   | paral_kgb=flag controlling (k,g,bands) parallelization
!!   | pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!   | usepaw=flag for PAW
!!  gprimd(3,3)= dimensional reciprocal space primitive translations
!!  mpi_enreg= informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc= dimension of the xccc3d array (0 or nfft).
!!  nfft= number of FFT grid points
!!  ngfft(18)= contain all needed information about 3D FFT
!!  nhat(nfft,nspden)=charge compensation density (content depends on electronpositron%particle)
!!  option= if 1, calculate positron lifetime for whole density
!!          if 2, calculate positron lifetime for given state
!!          if 3, calculate positron lifetime for given state with IPM
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  rhor(nfft,nspden)=total electron/positron density (content depends on electronpositron%particle)
!!  ucvol=unit cell volume in bohr**3.
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  ===== Optional arguments, used only if option>1 =====
!!  pawrhoij_dop_el(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies of one state
!!  rhor_dop_el(nfft)=electron density of given state for the state dependent scheme
!!  ===== Optional argument =====
!!  pawrhoij_ep(my_natom*usepaw) <type(pawrhoij_type)>= atomic occupancies to be used in place of
!!                                                      electronpositron%pawrhoij_ep
!!
!! OUTPUT
!!  rate= annihilation rate of a given state needed for state dependent scheme for doppler broadening
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      outscfcv,posdoppler
!!
!! CHILDREN
!!      gammapositron,gammapositron_fft,mkdenpos,nderiv_gen,pawdensities
!!      pawxcsum,simp_gen,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine poslifetime(dtset,electronpositron,gprimd,my_natom,mpi_enreg,n3xccc,nfft,ngfft,nhat,&
&                      option,pawang,pawrad,pawrhoij,pawtab,rate,rate_paw,rhor,ucvol,xccc3d,&
&                      rhor_dop_el,pawrhoij_dop_el,pawrhoij_ep) ! optional arguments


 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_electronpositron

 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, simp_gen, nderiv_gen
 use m_pawtab, only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawxc, only: pawxcsum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'poslifetime'
 use interfaces_14_hidewrite
 use interfaces_41_xc_lowlevel
 use interfaces_56_xc
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,n3xccc,nfft,option
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: rate,rate_paw
 type(dataset_type), intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,dtset%nspden*dtset%usepaw),xccc3d(n3xccc)
 real(dp),intent(in),target :: rhor(nfft,dtset%nspden)
 real(dp),optional,intent(in) :: rhor_dop_el(nfft)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
 type(pawrhoij_type),optional,intent(in) :: pawrhoij_dop_el(my_natom*dtset%usepaw)
 type(pawrhoij_type),optional,target,intent(in) :: pawrhoij_ep(my_natom*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,ierr,ifft,igam,ii,ilm,ilm1,ilm2,iloop,ipt,ir,isel
 integer :: itypat,iwarn,iwarnj,iwarnp,lm_size,lmn2_size,mesh_size
 integer :: nfftot,ngamma,ngr,ngrad,nspden_ep,opt_dens,usecore
 logical,parameter :: include_nhat_in_gamma=.false.
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,fact2,intg
 real(dp) :: lambda_core    ,lambda_core_ipm    ,lambda    ,lambda_ipm
 real(dp) :: lambda_core_paw,lambda_core_paw_ipm,lambda_paw,lambda_paw_ipm
 real(dp) :: lifetime,lifetime_ipm,nbec,nbev,nbp,rdum,sqfpi,units
 character(len=500) :: msg
!arrays
 integer,allocatable :: igamma(:)
 logical,allocatable :: lmselect(:),lmselect_ep(:),lmselect_dum(:)
 real(dp) :: mpibuf(4)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/),lsign(2)=(/one,-one/)
 real(dp),allocatable :: d1gam(:,:,:),d2gam(:,:,:),ff(:),gam_(:,:,:),gamma(:,:),gammam(:,:,:),gg(:,:)
 real(dp),allocatable :: grhocore2(:),grhocor2_(:),grhoe2(:),grho2_(:)
 real(dp),allocatable :: nhat1(:,:,:),nhat1_ep(:,:,:),nhat1_j(:,:,:)
 real(dp),allocatable :: rho_(:),rho_ep_(:),rho1(:,:,:),rho1_ep(:,:,:),rho1_j(:,:,:)
 real(dp),allocatable :: rhoarr1(:),rhoarr1_ep(:),rhoarr1_j(:),rhoarr2(:)
 real(dp),allocatable :: rhocore(:),rhocor_(:),rhoe(:,:),rhop(:,:),rhor_dop_el_(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),rhosph_j(:),rhotot(:,:),rhotot_ep(:,:)
 real(dp),allocatable :: rhotot_j(:,:),trho1(:,:,:),trho1_ep(:,:,:),trho1_j(:,:,:)
 real(dp),allocatable :: v1sum(:,:),v2sum(:,:,:)
 real(dp),pointer :: rhor_(:,:),rhor_ep_(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_ep_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Tests for developers
 if (.not.associated(electronpositron)) then
   msg='electronpositron variable must be associated!'
   MSG_BUG(msg)
 end if
 if (option/=1) then
   if ((.not.present(rhor_dop_el)).or.(.not.present(pawrhoij_dop_el))) then
     msg='when option/=1, rhor_dop_el and pawrhoij_dop_el must be present!'
     MSG_BUG(msg)
   end if
 end if

!Constants
 fact=0.0
 cplex=1;nspden_ep=1
 usecore=n3xccc/nfft
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 iwarn=0;iwarnj=0;iwarnp=1
 sqfpi=sqrt(four_pi)

!Compatibility tests
 if (electronpositron%particle==EP_NOTHING) then
   msg='Not valid for electronpositron%particle=NOTHING!'
   MSG_BUG(msg)
 end if
 if (electronpositron%nfft/=nfft) then
   msg='nfft/=electronpositron%nfft!'
   MSG_BUG(msg)
 end if
 if (dtset%usepaw==1) then
   if(dtset%pawxcdev==0.and.ngrad==2) then
     msg='GGA is not implemented for pawxcdev=0 (use dtset%pawxcdev/=0)!'
     MSG_BUG(msg)
   end if
 end if

!Select type(s) of enhancement factor
 if ((electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3).and.option==1) then
   ngamma=2
   ABI_ALLOCATE(igamma,(ngamma))
   igamma(1)=1;igamma(2)=2
 else
   ngamma=1
   ABI_ALLOCATE(igamma,(ngamma))
   if (electronpositron%ixcpositron==-1) igamma(1)=0
   if (electronpositron%ixcpositron== 2) igamma(1)=4
   if (electronpositron%ixcpositron==11.or.electronpositron%ixcpositron==31) igamma(1)=3
   if (electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3) igamma(1)=2
 end if

!Select density according to nhat choice
 if (dtset%usepaw==0.or.include_nhat_in_gamma) then
   rhor_ => rhor
   rhor_ep_ => electronpositron%rhor_ep
 else
   ABI_ALLOCATE(rhor_,(nfft,dtset%nspden))
   ABI_ALLOCATE(rhor_ep_,(nfft,dtset%nspden))
   rhor_=rhor-nhat
   rhor_ep_=electronpositron%rhor_ep-electronpositron%nhat_ep
 end if

!Eventually overwrite electronpositron%pawrhoij_ep
 if (present(pawrhoij_ep)) then
   pawrhoij_ep_ => pawrhoij_ep
 else
   pawrhoij_ep_ => electronpositron%pawrhoij_ep
 end if

!Loop on different enhancement factors
 do igam=1,ngamma

!  Compute electron-positron annihilation rate using pseudo densities (plane waves)
!  ----------------------------------------------------------------------------------------

!  Select the densities and make them positive
   ABI_ALLOCATE(rhoe,(nfft,nspden_ep))
   ABI_ALLOCATE(rhop,(nfft,nspden_ep))
   if (electronpositron%particle==EP_ELECTRON) then
     rhoe(:,1)=rhor_ep_(:,1);rhop(:,1)=rhor_(:,1)
   else if (electronpositron%particle==EP_POSITRON) then
     rhoe(:,1)=rhor_(:,1);rhop(:,1)=rhor_ep_(:,1)
   end if
   call mkdenpos(iwarn ,nfft,nspden_ep,1,rhoe,dtset%xc_denpos)
   call mkdenpos(iwarnp,nfft,nspden_ep,1,rhop,dtset%xc_denpos)
   if (option/=1) then
     ABI_ALLOCATE(rhor_dop_el_,(nfft))
     rhor_dop_el_(:)=rhor_dop_el(:)
     call mkdenpos(iwarnp,nfft,1,1,rhor_dop_el_,dtset%xc_denpos)
   end if

!  Compute enhancement factor at each FFT grid point
!  gamma(:,1): using total   electronic density
!  gamma(:,2): using valence electronic density
   ABI_ALLOCATE(gamma,(nfft,2))
   if (option==1.or.option==2) then
     call gammapositron_fft(electronpositron,gamma,gprimd,igamma(igam),mpi_enreg,&
&     n3xccc,nfft,ngfft,rhoe,rhop,xccc3d)
   else
     gamma=one
   end if

!  Compute positron annihilation rates
   lambda     =zero;lambda_ipm     =zero
   lambda_core=zero;lambda_core_ipm=zero
   if (option==1) then
     do ifft=1,nfft
       lambda    =lambda    +rhop(ifft,1)*rhoe(ifft,1)*gamma(ifft,1)
       lambda_ipm=lambda_ipm+rhop(ifft,1)*rhoe(ifft,1)*gamma(ifft,2)
     end do
   else
     do ifft=1,nfft
       lambda    =lambda    +rhop(ifft,1)*rhor_dop_el_(ifft)*gamma(ifft,1)
       lambda_ipm=lambda_ipm+rhop(ifft,1)*rhor_dop_el_(ifft)*gamma(ifft,2)
     end do
   end if
   if (usecore==1) then
     do ifft=1,nfft
       lambda_core    =lambda_core    +rhop(ifft,1)*xccc3d(ifft)*gamma(ifft,1)
       lambda_core_ipm=lambda_core_ipm+rhop(ifft,1)*xccc3d(ifft)
     end do
   end if
   lambda         =lambda         *ucvol/dble(nfftot)
   lambda_ipm     =lambda_ipm     *ucvol/dble(nfftot)
   lambda_core    =lambda_core    *ucvol/dble(nfftot)
   lambda_core_ipm=lambda_core_ipm*ucvol/dble(nfftot)
   ABI_DEALLOCATE(gamma)
   ABI_DEALLOCATE(rhoe)
   ABI_DEALLOCATE(rhop)
   if (option/=1) then
     ABI_DEALLOCATE(rhor_dop_el_)
   end if
!  NC pseudopotential: check electrons/positron number
   if (dtset%usepaw==0.and.igam==ngamma) then
     nbec=zero;nbev=zero;nbp=zero
     if (electronpositron%particle==EP_ELECTRON) then
       do ifft=1,nfft
         nbec=nbec+xccc3d(ifft)
         nbev=nbev+electronpositron%rhor_ep(ifft,1)
         nbp =nbp +rhor(ifft,1)
       end do
     else
       do ifft=1,nfft
         nbec=nbec+xccc3d(ifft)
         nbev=nbev+rhor(ifft,1)
         nbp =nbp +electronpositron%rhor_ep(ifft,1)
       end do
     end if
     nbec=nbec*ucvol/dble(nfftot)
     nbev=nbev*ucvol/dble(nfftot)
     nbp =nbp *ucvol/dble(nfftot)
   end if

!  MPI parallelization
   if(mpi_enreg%nproc_fft>1)then
     call xmpi_sum(lambda    ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_ipm,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_core    ,mpi_enreg%comm_fft,ierr)
     call xmpi_sum(lambda_core_ipm,mpi_enreg%comm_fft,ierr)
     if (dtset%usepaw==0.and.igam==ngamma) then
       call xmpi_sum(nbec,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(nbev,mpi_enreg%comm_fft,ierr)
       call xmpi_sum(nbp ,mpi_enreg%comm_fft,ierr)
     end if
   end if


!  PAW: add on-site contributions to electron-positron annihilation rate
!  ----------------------------------------------------------------------------------------
   if (dtset%usepaw==1) then

     lambda_paw     =zero;lambda_paw_ipm     =zero
     lambda_core_paw=zero;lambda_core_paw_ipm=zero

!    Loop on atoms
     do iatom=1,my_natom

       itypat=pawrhoij(iatom)%itypat
       lmn2_size=pawtab(itypat)%lmn2_size
       mesh_size=pawtab(itypat)%mesh_size
       lm_size=pawtab(itypat)%lcut_size**2
       cplex=1
       ngr=0;if (ngrad==2) ngr=mesh_size

!      Allocations of "on-site" densities
       ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden_ep))
       ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden_ep))
       if (option/=1) then
         ABI_ALLOCATE(rho1_j ,(cplex*mesh_size,lm_size,nspden_ep))
         ABI_ALLOCATE(trho1_j,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(rho1_j ,(0,0,0))
         ABI_ALLOCATE(trho1_j,(0,0,0))
       end if
       if (include_nhat_in_gamma) then
         ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden_ep))
         ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(nhat1,(0,0,0))
         ABI_ALLOCATE(nhat1_ep,(0,0,0))
       end if
       if (include_nhat_in_gamma.and.option/=1) then
         ABI_ALLOCATE(nhat1_j,(cplex*mesh_size,lm_size,nspden_ep))
       else
         ABI_ALLOCATE(nhat1_j,(0,0,0))
       end if
       ABI_ALLOCATE(lmselect,(lm_size))
       ABI_ALLOCATE(lmselect_ep,(lm_size))
       ABI_ALLOCATE(lmselect_dum,(lm_size))

!      Compute "on-site" densities (n1, ntild1, nhat1) for electron and positron =====
       lmselect(:)=.true.
       opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
       call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1,nspden_ep,1,&
&       0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij(iatom),&
&       pawtab(itypat),rho1,trho1)
       lmselect_ep(:)=.true.
       call pawdensities(rdum,cplex,iatom,lmselect_ep,lmselect_dum,lm_size,nhat1_ep,nspden_ep,1,&
&       0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_ep_(iatom),&
&       pawtab(itypat),rho1_ep,trho1_ep)

!      For state dependent scheme in Doppler                                       =====
!      Compute "on-site" densities (n1, ntild1, nhat1) for a given electron state j=====
       if (option/=1) then
         opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
         call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1_j,nspden_ep,1,&
&         0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_dop_el(iatom),&
&         pawtab(itypat),rho1_j,trho1_j)
       end if

!      Compute contribution to annihilation rate:
!      Loop: first step: compute all-electron contribution (from n^1, n_c)
!      2nd   step: compute pseudo contribution (from tild_n^1, hat_n^1, tild_n_c)
       do iloop=1,2
         if (iloop==1) usecore=1
         if (iloop==2) usecore=pawtab(itypat)%usetcore
         ABI_ALLOCATE(rhocore,(mesh_size))

!        First formalism: use densities on r,theta,phi
         if (dtset%pawxcdev==0) then

           ABI_ALLOCATE(gamma,(mesh_size,2))
           ABI_ALLOCATE(rhoarr1,(mesh_size))
           ABI_ALLOCATE(rhoarr1_ep,(mesh_size))
           if (option/=1) then
             ABI_ALLOCATE(rhoarr1_j,(mesh_size))
           end if
!          Loop on the angular part
           do ipt=1,pawang%angl_size
!            Build densities
             rhoarr1=zero;rhoarr1_ep=zero;rhocore=zero
             if (option/=1) rhoarr1_j=zero
             if (iloop==1) then
               do ilm=1,lm_size
                 if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+rho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
               end do
               if (option/=1) then
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1_j(:)+rho1_j(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
               end if
               do ilm=1,lm_size
                 if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+rho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
               end do
               if (usecore==1) rhocore(:)=pawtab(itypat)%coredens(:)
             else
               if (include_nhat_in_gamma) then
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+(trho1(:,ilm,1)+nhat1(:,ilm,1))*pawang%ylmr(ilm,ipt)
                 end do
                 if (option/=1) then
                   do ilm=1,lm_size
                     if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1_j(:)+(trho1_j(:,ilm,1)+nhat1_j(:,ilm,1))*pawang%ylmr(ilm,ipt)
                   end do
                 end if
                 do ilm=1,lm_size
                   if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+(trho1_ep(:,ilm,1)+nhat1_ep(:,ilm,1))*pawang%ylmr(ilm,ipt)
                 end do
               else
                 do ilm=1,lm_size
                   if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+trho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
                 if (option/=1) then
                   do ilm=1,lm_size
                     if (lmselect(ilm)) rhoarr1_j(:)=rhoarr1(:)+trho1_j(:,ilm,1)*pawang%ylmr(ilm,ipt)
                   end do
                 end if
                 do ilm=1,lm_size
                   if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+trho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
                 end do
               end if
               if (usecore==1) rhocore(:)=pawtab(itypat)%tcoredens(:,1)
             end if
!            Make the densities positive
             if (electronpositron%particle==EP_ELECTRON) then
               call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
               call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
             else if (electronpositron%particle==EP_POSITRON) then
               call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
               call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
               if (option/=1) then
                 call mkdenpos(iwarnj,mesh_size,1,1,rhoarr1_j,dtset%xc_denpos)
               end if
             end if
!            Compute Gamma
             ABI_ALLOCATE(grhoe2,(ngr))
             ABI_ALLOCATE(grhocore2,(ngr))
             if (option==1.or.option==2) then
               if (electronpositron%particle==EP_ELECTRON) then
                 call gammapositron(gamma,grhocore2,grhoe2,igamma(igam),ngr,mesh_size,&
&                 rhocore,rhoarr1_ep,rhoarr1,usecore)
               else if (electronpositron%particle==EP_POSITRON) then
                 call gammapositron(gamma,grhocore2,grhoe2,igamma(igam),ngr,mesh_size,&
&                 rhocore,rhoarr1,rhoarr1_ep,usecore)
               end if
             else
               gamma(:,:)=one
             end if
             ABI_DEALLOCATE(grhoe2)
             ABI_DEALLOCATE(grhocore2)
!            Compute contribution to annihilation rates
             ABI_ALLOCATE(ff,(mesh_size))
             if (option/=1) rhoarr1(:)=rhoarr1_j(:)
             do ii=1,4
               if (ii==1) ff(:)=rhoarr1(:)*rhoarr1_ep(:)*gamma(:,1)*pawrad(itypat)%rad(:)**2
               if (ii==2) ff(:)=rhoarr1(:)*rhoarr1_ep(:)*gamma(:,2)*pawrad(itypat)%rad(:)**2
               if (electronpositron%particle==EP_ELECTRON) then
                 if (ii==3) ff(:)=rhoarr1   (:)*rhocore(:)*gamma(:,1)*pawrad(itypat)%rad(:)**2
                 if (ii==4) ff(:)=rhoarr1   (:)*rhocore(:)           *pawrad(itypat)%rad(:)**2
               else
                 if (ii==3) ff(:)=rhoarr1_ep(:)*rhocore(:)*gamma(:,1)*pawrad(itypat)%rad(:)**2
                 if (ii==4) ff(:)=rhoarr1_ep(:)*rhocore(:)           *pawrad(itypat)%rad(:)**2
               end if
               call simp_gen(intg,ff,pawrad(itypat))
               intg=intg*pawang%angwgth(ipt)*four_pi
               if (ii==1) lambda_paw         =lambda_paw         +lsign(iloop)*intg
               if (ii==2) lambda_paw_ipm     =lambda_paw_ipm     +lsign(iloop)*intg
               if (ii==3) lambda_core_paw    =lambda_core_paw    +lsign(iloop)*intg
               if (ii==4) lambda_core_paw_ipm=lambda_core_paw_ipm+lsign(iloop)*intg
             end do
             ABI_DEALLOCATE(ff)
           end do ! ipt
           ABI_DEALLOCATE(gamma)
           ABI_DEALLOCATE(rhoarr1)
           ABI_DEALLOCATE(rhoarr1_ep)
           if (option/=1) then
             ABI_DEALLOCATE(rhoarr1_j)
           end if
           
!          Second formalism: use (l,m) moments for densities
         else if (dtset%pawxcdev/=0) then

!          Build densities
           ABI_ALLOCATE(gammam,(mesh_size,2,lm_size))
           ABI_ALLOCATE(rhotot,(mesh_size,lm_size))
           ABI_ALLOCATE(rhotot_ep,(mesh_size,lm_size))
           ABI_ALLOCATE(rhosph,(mesh_size))
           ABI_ALLOCATE(rhosph_ep,(mesh_size))
           if (option/=1) then
             ABI_ALLOCATE(rhotot_j,(mesh_size,lm_size))
             ABI_ALLOCATE(rhosph_j,(mesh_size))
           end if
           if (usecore==0) rhocore(:)=zero
           if (iloop==1) then
             rhotot   (:,:)=rho1   (:,:,1)
             rhotot_ep(:,:)=rho1_ep(:,:,1)
             if (option/=1) rhotot_j (:,:)=rho1_j (:,:,1)
             if (usecore==1) rhocore(:)=pawtab(itypat)%coredens(:)
           else
             if (include_nhat_in_gamma) then
               rhotot   (:,:)=trho1   (:,:,1)+nhat1   (:,:,1)
               rhotot_ep(:,:)=trho1_ep(:,:,1)+nhat1_ep(:,:,1)
               if (option/=1) rhotot_j (:,:)=trho1_j (:,:,1)+nhat1_j (:,:,1)
             else
               rhotot   (:,:)=trho1   (:,:,1)
               rhotot_ep(:,:)=trho1_ep(:,:,1)
               if (option/=1) rhotot_j (:,:)=trho1_j (:,:,1)
             end if
             if (usecore==1) rhocore(:)=pawtab(itypat)%tcoredens(:,1)
           end if
           rhosph   (:)=rhotot   (:,1)/sqfpi
           rhosph_ep(:)=rhotot_ep(:,1)/sqfpi
           if (option/=1) rhosph_j (:)=rhotot_j(:,1)/sqfpi
!          Make spherical densities positive
           if (electronpositron%particle==EP_ELECTRON) then
             call mkdenpos(iwarnp,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
             call mkdenpos(iwarn ,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
           else if (electronpositron%particle==EP_POSITRON) then
             call mkdenpos(iwarn ,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
             call mkdenpos(iwarnp,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
             if (option/=1) then
               call mkdenpos(iwarnp,mesh_size,1,1,rhosph_j,dtset%xc_denpos)
             end if
           end if
!          Need gradients of electronic densities for GGA
           ABI_ALLOCATE(grhoe2,(ngr))
           ABI_ALLOCATE(grhocore2,(ngr))
           if (ngr>0) then
             if (electronpositron%particle==EP_ELECTRON) then
               call nderiv_gen(grhoe2,rhosph_ep,pawrad(itypat))
             else if (electronpositron%particle==EP_POSITRON) then
               call nderiv_gen(grhoe2,rhosph,pawrad(itypat))
             end if
             grhoe2(:)=grhoe2(:)**2
             if (usecore==1) then
               call nderiv_gen(grhocore2,rhocore,pawrad(itypat))
               grhocore2(:)=grhocore2(:)**2
             end if
           end if
!          Compute Gamma for (rho-,rho+),
!          (rho- +drho-,rho+), (rho- -drho-,rho+),
!          (rho-,rho+ +drho+), (rho-,rho+ -drho+),
!          (rho- +drho-,rho+ +drho+), (rho- -drho-,rho+ -drho+)
!          Do a seven steps loop
           ABI_ALLOCATE(gam_,(mesh_size,2,7))
           ABI_ALLOCATE(rho_,(mesh_size))
           ABI_ALLOCATE(rho_ep_,(mesh_size))
           ABI_ALLOCATE(rhocor_,(mesh_size))
           ABI_ALLOCATE(grho2_,(ngr))
           ABI_ALLOCATE(grhocor2_,(ngr))
           do ii=1,7
!            Apply delta to get perturbed densities
             rho_(:)=rhosph(:);rho_ep_(:)=rhosph_ep(:);if (usecore==1) rhocor_(:)=rhocore(:)
             if (ngr>0) grho2_(:)=grhoe2(:)
             if (ngr>0) grhocor2_(:)=grhocore2(:)
             if (ii==2.or.ii==4.or.ii==6) fact=(one+delta)
             if (ii==3.or.ii==5.or.ii==7) fact=(one-delta)
             fact2=fact**2
             if (ii==2.or.ii==3.or.ii==6.or.ii==7) then
               rho_(:)=fact*rho_(:)
               if (electronpositron%particle==EP_POSITRON) then
                 if (ngr>0) grho2_(:)=fact2*grho2_(:)
                 if (usecore==1)rhocor_(:)=fact*rhocor_(:)
                 if (ngr>0.and.usecore==1) grhocor2_(:)=fact2*grhocor2_(:)
               end if
             end if
             if (ii==4.or.ii==5.or.ii==6.or.ii==7) then
               rho_ep_(:)=fact*rho_ep_(:)
               if (electronpositron%particle==EP_ELECTRON) then
                 if (ngr>0) grho2_(:)=fact2*grho2_(:)
                 if (usecore==1)rhocor_(:)=fact*rhocor_(:)
                 if (ngr>0.and.usecore==1) grhocor2_(:)=fact2*grhocor2_(:)
               end if
             end if
!            Compute gamma for these perturbed densities
             if (option==1.or.option==2) then
               if (electronpositron%particle==EP_ELECTRON) then
                 call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma(igam),ngr,mesh_size,rhocor_,rho_ep_,rho_,usecore)
               else if (electronpositron%particle==EP_POSITRON) then
                 call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma(igam),ngr,mesh_size,rhocor_,rho_,rho_ep_,usecore)
               end if
             else
               gam_(:,:,:)=one
             end if
           end do ! end loop ii=1,7

           ABI_DEALLOCATE(rhocor_)
           ABI_DEALLOCATE(grho2_)
           ABI_DEALLOCATE(grhocor2_)
           ABI_DEALLOCATE(grhoe2)
           ABI_DEALLOCATE(grhocore2)
           rho_   (:)=rhosph   (:);if (electronpositron%particle==EP_POSITRON.and.usecore==1) rho_   (:)=rho_   (:)+rhocore(:)
           rho_ep_(:)=rhosph_ep(:);if (electronpositron%particle==EP_ELECTRON.and.usecore==1) rho_ep_(:)=rho_ep_(:)+rhocore(:)
!          Compute numerical first and second derivatives of Gamma
!          d1gam(1) = dgam/drho+ (particle=ELECTRON), dgam/drho- (particle=POSITRON)
!          d1gam(2) = dgam/drho- (particle=ELECTRON), dgam/drho+ (particle=POSITRON)
           ABI_ALLOCATE(d1gam,(mesh_size,2,2))
           d1gam(:,:,:)=zero
           do ir=1,mesh_size
             if (rho_     (ir)>tol14) d1gam(ir,1,1)=(gam_(ir,1,2)-gam_(ir,1,3))*half/(delta*rho_     (ir))
             if (rhosph   (ir)>tol14) d1gam(ir,2,1)=(gam_(ir,2,2)-gam_(ir,2,3))*half/(delta*rhosph   (ir))
             if (rho_ep_  (ir)>tol14) d1gam(ir,1,2)=(gam_(ir,1,4)-gam_(ir,1,5))*half/(delta*rho_ep_  (ir))
             if (rhosph_ep(ir)>tol14) d1gam(ir,2,2)=(gam_(ir,2,4)-gam_(ir,2,5))*half/(delta*rhosph_ep(ir))
           end do

!          d2gam(1) = d2gam/drho+_drho+ (particle=ELECTRON), dgam/drho-_drho- (particle=POSITRON)
!          d2gam(2) = d2gam/drho-_drho+ (particle=ELECTRON), dgam/drho+_drho- (particle=POSITRON)
!          d2gam(3) = d2gam/drho-_drho- (particle=ELECTRON), dgam/drho+_drho+ (particle=POSITRON)
           ABI_ALLOCATE(d2gam,(mesh_size,2,3))
           d2gam(:,:,:)=zero
           do ir=1,mesh_size
             if (rho_  (ir)>tol14) d2gam(ir,1,1)=(gam_(ir,1,2)+gam_(ir,1,3)-two*gam_(ir,1,1))/(delta*rho_  (ir))**2
             if (rhosph(ir)>tol14) d2gam(ir,2,1)=(gam_(ir,2,2)+gam_(ir,2,3)-two*gam_(ir,2,1))/(delta*rhosph(ir))**2
             if (rho_ep_(ir)>tol14) then
               d2gam(ir,1,3)=(gam_(ir,1,4)+gam_(ir,1,5)-two*gam_(ir,1,1))/(delta*rho_ep_(ir))**2
               if (rho_(ir)>tol14) then
                 d2gam(ir,1,2)=(gam_(ir,1,6)+gam_(ir,1,7)+two*gam_(ir,1,1) &
&                 -gam_(ir,1,2)-gam_(ir,1,3)-gam_(ir,1,4)-gam_(ir,1,5)) &
&                 *half/(delta*rho_(ir))/(delta*rho_ep_(ir))
               end if
             end if
             if (rhosph_ep(ir)>tol14) then
               d2gam(ir,2,3)=(gam_(ir,2,4)+gam_(ir,2,5)-two*gam_(ir,2,1))/(delta*rhosph_ep(ir))**2
               if (rhosph(ir)>tol14) then
                 d2gam(ir,2,2)=(gam_(ir,2,6)+gam_(ir,2,7)+two*gam_(ir,2,1) &
&                 -gam_(ir,2,2)-gam_(ir,2,3)-gam_(ir,2,4)-gam_(ir,2,5)) &
&                 *half/(delta*rhosph(ir))/(delta*rhosph_ep(ir))
               end if
             end if
           end do
           ABI_DEALLOCATE(rho_)
           ABI_DEALLOCATE(rho_ep_)
!          Compute useful sums of densities
           ABI_ALLOCATE(v1sum,(mesh_size,3))
           if ( dtset%pawxcdev>=2)  then
             ABI_ALLOCATE(v2sum,(mesh_size,lm_size,3))
           else
             ABI_ALLOCATE(v2sum,(0,0,0))
           end if
           rhotot(:,1)=sqfpi*rhosph(:);rhotot_ep(:,1)=sqfpi*rhosph_ep(:)
           call pawxcsum(1,1,1,lmselect,lmselect_ep,lm_size,mesh_size,3,dtset%pawxcdev,&
&           pawang,rhotot,rhotot_ep,v1sum,v2sum)
!          Compute final developpment of gamma moments
           gammam(:,:,:)=zero
           gammam(:,:,1)=gam_(:,:,1)*sqfpi
           gammam(:,1,1)=gammam(:,1,1)+(d2gam(:,1,2)*v1sum(:,2) &
&           +half*(d2gam(:,1,1)*v1sum(:,1)+d2gam(:,1,3)*v1sum(:,3)))/sqfpi
           gammam(:,2,1)=gammam(:,2,1)+(d2gam(:,2,2)*v1sum(:,2) &
&           +half*(d2gam(:,2,1)*v1sum(:,1)+d2gam(:,2,3)*v1sum(:,3)))/sqfpi
           do ilm=2,lm_size
             if (lmselect(ilm)) then
               gammam(:,1,ilm)=gammam(:,1,ilm)+d1gam(:,1,1)*rhotot(:,ilm)
               gammam(:,2,ilm)=gammam(:,2,ilm)+d1gam(:,2,1)*rhotot(:,ilm)
             end if
             if (lmselect_ep(ilm)) then
               gammam(:,1,ilm)=gammam(:,1,ilm)+d1gam(:,1,2)*rhotot_ep(:,ilm)
               gammam(:,2,ilm)=gammam(:,2,ilm)+d1gam(:,2,2)*rhotot_ep(:,ilm)
             end if
           end do
           if (dtset%pawxcdev>1) then
             do ilm=2,lm_size
               gammam(:,1,ilm)=gammam(:,1,ilm)+d2gam(:,1,2)*v2sum(:,ilm,2) &
&               +half*(d2gam(:,1,1)*v2sum(:,ilm,1)+d2gam(:,1,3)*v2sum(:,ilm,3))
               gammam(:,2,ilm)=gammam(:,2,ilm)+d2gam(:,2,2)*v2sum(:,ilm,2) &
&               +half*(d2gam(:,2,1)*v2sum(:,ilm,1)+d2gam(:,2,3)*v2sum(:,ilm,3))
             end do
           end if
           ABI_DEALLOCATE(gam_)
           ABI_DEALLOCATE(d1gam)
           ABI_DEALLOCATE(d2gam)
           ABI_DEALLOCATE(v1sum)
           ABI_DEALLOCATE(v2sum)

!          Compute contribution to annihilation rate

!          In state dependent scheme for Doppler, replace electronic density with one state
           if (option/=1) then
             rhotot(:,:) = rhotot_j(:,:)
             rhosph  (:) = rhosph_j  (:)
           end if

           ABI_ALLOCATE(gg,(mesh_size,4))
           gg=zero
           ABI_ALLOCATE(rhoarr1,(mesh_size))
           ABI_ALLOCATE(rhoarr2,(mesh_size))
           do ilm=1,lm_size
             do ilm1=1,lm_size
               if (lmselect(ilm1)) then
                 if (ilm1==1) rhoarr1(:)=sqfpi*rhosph(:)
                 if (ilm1/=1) rhoarr1(:)=rhotot(:,ilm1)
                 do ilm2=1,lm_size
                   if (lmselect_ep(ilm2)) then
                     if (ilm2==1) rhoarr2(:)=sqfpi*rhosph_ep(:)
                     if (ilm2/=1) rhoarr2(:)=rhotot_ep(:,ilm2)
                     if (ilm1>=ilm2) then
                       isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                     else
                       isel=pawang%gntselect(ilm,ilm1+ilm2*(ilm2-1)/2)
                     end if
                     if (isel>0) then
                       fact=pawang%realgnt(isel)
                       gg(:,1)=gg(:,1)+fact*rhoarr1(:)*rhoarr2(:)*gammam(:,1,ilm)
                       gg(:,2)=gg(:,2)+fact*rhoarr1(:)*rhoarr2(:)*gammam(:,2,ilm)
                     end if
                   end if
                 end do
               end if
             end do
           end do
           ABI_DEALLOCATE(rhoarr1)
           ABI_DEALLOCATE(rhoarr2)
           if (electronpositron%particle==EP_ELECTRON) then
             do ilm=1,lm_size
               if (lmselect(ilm)) gg(:,3)=gg(:,3)+rhotot(:,ilm)*rhocore(:)*gammam(:,1,ilm)
             end do
             gg(:,4)=sqfpi*rhotot(:,1)*rhocore(:)
           else if (electronpositron%particle==EP_POSITRON) then
             do ilm=1,lm_size
               if (lmselect_ep(ilm)) gg(:,3)=gg(:,3)+rhotot_ep(:,ilm)*rhocore(:)*gammam(:,1,ilm)
             end do
             gg(:,4)=sqfpi*rhotot_ep(:,1)*rhocore(:)
           end if
           do ii=1,4
             gg(:,ii)=gg(:,ii)*pawrad(itypat)%rad(:)**2
             call simp_gen(intg,gg(:,ii),pawrad(itypat))
             if (ii==1) lambda_paw         =lambda_paw         +lsign(iloop)*intg
             if (ii==2) lambda_paw_ipm     =lambda_paw_ipm     +lsign(iloop)*intg
             if (ii==3) lambda_core_paw    =lambda_core_paw    +lsign(iloop)*intg
             if (ii==4) lambda_core_paw_ipm=lambda_core_paw_ipm+lsign(iloop)*intg
           end do
           ABI_DEALLOCATE(gg)
           ABI_DEALLOCATE(gammam)
           ABI_DEALLOCATE(rhotot)
           ABI_DEALLOCATE(rhotot_ep)
           ABI_DEALLOCATE(rhosph)
           ABI_DEALLOCATE(rhosph_ep)
           if (option/=1) then
             ABI_DEALLOCATE(rhotot_j)
             ABI_DEALLOCATE(rhosph_j)
           end if

         end if ! dtset%pawxcdev

         ABI_DEALLOCATE(rhocore)

       end do ! iloop

       ABI_DEALLOCATE(rho1)
       ABI_DEALLOCATE(trho1)
       ABI_DEALLOCATE(rho1_ep)
       ABI_DEALLOCATE(trho1_ep)
       ABI_DEALLOCATE(rho1_j)
       ABI_DEALLOCATE(trho1_j)
       ABI_DEALLOCATE(nhat1)
       ABI_DEALLOCATE(nhat1_ep)
       ABI_DEALLOCATE(nhat1_j)
       ABI_DEALLOCATE(lmselect)
       ABI_DEALLOCATE(lmselect_ep)
       ABI_DEALLOCATE(lmselect_dum)

     end do ! iatom

!    Reduction in case of distribution over atomic sites
     if (mpi_enreg%nproc_atom>1) then
       mpibuf(1)=lambda_paw     ;mpibuf(2)=lambda_paw_ipm
       mpibuf(3)=lambda_core_paw;mpibuf(4)=lambda_core_paw_ipm
       call xmpi_sum(mpibuf,mpi_enreg%comm_atom,ierr)
       lambda_paw=mpibuf(1)     ;lambda_paw_ipm=mpibuf(2)
       lambda_core_paw=mpibuf(3);lambda_core_paw_ipm=mpibuf(4)
     end if

!    Add plane-wave and PAW contributions to annihilation rates
     lambda         =lambda         +lambda_paw
     lambda_ipm     =lambda_ipm     +lambda_paw_ipm
     lambda_core    =lambda_core    +lambda_core_paw
     lambda_core_ipm=lambda_core_ipm+lambda_core_paw_ipm
   end if ! dtset%usepaw


!  Convert into proper units and print
!  ---------------------------------------------------------------------------------------

!  Sum valence and core contributions to annihilation rates
   lambda        =lambda        +lambda_core
   lambda_ipm    =lambda_ipm    +lambda_core_ipm
   if (dtset%usepaw==1) then
     lambda_paw    =lambda_paw    +lambda_core_paw
     lambda_paw_ipm=lambda_paw_ipm+lambda_core_paw_ipm
   end if

!  Set annihilation rate in proper unit (picosec.)
   units=pi*(one/InvFineStruct)**3/Time_Sec/1.e12_dp/electronpositron%posocc

   lambda         =lambda         *units
   lambda_ipm     =lambda_ipm     *units
   lambda_core    =lambda_core    *units
   lambda_core_ipm=lambda_core_ipm*units
   lifetime    =one/lambda
   lifetime_ipm=one/lambda_ipm
   electronpositron%lambda=lambda
   electronpositron%lifetime=lifetime
   if (dtset%usepaw==1) then
     lambda_paw         =lambda_paw         *units
     lambda_paw_ipm     =lambda_paw_ipm     *units
     lambda_core_paw    =lambda_core_paw    *units
     lambda_core_paw_ipm=lambda_core_paw_ipm*units
   end if
   rate=lambda-lambda_core-lambda_paw+lambda_core_paw
   rate_paw=lambda_paw-lambda_core_paw
!  Print life time and additional information
   if (option==1) then
     if (igam==1) then
       write(msg,'(a,80("-"),2a)') ch10,ch10,' Results for electron-positron annihilation:'
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (ngamma>1.and.igam==1) then
       write(msg,'(a,i1,a)') ch10,ngamma,&
&       ' computations of positron lifetime have been performed (with different enhancement factors).'
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (ngamma>1) then
       write(msg,'(2a,i1)') ch10,"########## Lifetime computation ",igam
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (abs(electronpositron%ixcpositron)==1) then
       write(msg,'(4a)') ch10,' # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)'
     else if (electronpositron%ixcpositron==11) then
       write(msg,'(4a)') ch10,' # Zero-positron density limit of Arponen and Pajanne fitted by Sterne & Kaiser',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)'
     else if (electronpositron%ixcpositron==2) then
       write(msg,'(4a)') ch10,' # Electron-positron correlation provided by Puska, Seitsonen, and Nieminen',&
&       ch10,'   Ref: M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994)'
     else if (electronpositron%ixcpositron==3) then
       write(msg,'(8a)') ch10,' # Zero-positron density limit of Arponen and Pajanne provided by Boronski & Nieminen',&
&       ch10,'   + GGA corrections',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)',&
&       ch10,'         B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1994)'
     else if (electronpositron%ixcpositron==31) then
       write(msg,'(8a)') ch10,' # Zero-positron density limit of Arponen and Pajanne fitted by Sterne & Kaiser',&
&       ch10,'   + GGA corrections',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)',&
&       ch10,'         B. Barbiellini, M.J. Puska, T. Torsti and R.M. Nieminen, Phys. Rev. B 51, 7341 (1994)'
     end if
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,  msg,'COLL')
     if (igamma(igam)==0) then
       write(msg,'(a)')       ' # Enhancement factor set to one (test)'
     else if (igamma(igam)==1) then
       write(msg,'(3a)')      ' # Enhancement factor of Boronski & Nieminen',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)'
     else if (igamma(igam)==2) then
       write(msg,'(3a)')      ' # Enhancement factor of Boronski & Nieminen IN THE RPA LIMIT',&
&       ch10,'   Ref.: Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986)'
     else if (igamma(igam)==3) then
       write(msg,'(3a)')      ' # Enhancement factor of Sterne & Kaiser',&
&       ch10,'   Ref.: P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991)'
     else if (igamma(igam)==4) then
       write(msg,'(3a)')      ' # Enhancement factor of Puska, Seitsonen, and Nieminen',&
&       ch10,'   Ref.: M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994)'
     end if
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg, '(4(2a,es16.8))' ) ch10,&
&     ' Positron lifetime                         (ps)   =',lifetime    ,ch10,&
&     ' Positron lifetime with IPM for core elec. (ps)   =',lifetime_ipm,ch10,&
&     ' Annihilation rate                         (ns-1) =',lambda    *1000._dp,ch10,&
&     ' Annihilation rate with IPM for core elec. (ns-1) =',lambda_ipm*1000._dp
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     write(msg,'(2a,5(2a,es16.8))' ) ch10,&
&     ' Annihilation rate core/valence decomposition:',ch10,&
&     '   Core    contribution to ann.rate          (ns-1) =', lambda_core                 *1000._dp,ch10,&
&     '   Valence contribution to ann.rate          (ns-1) =',(lambda-lambda_core)         *1000._dp,ch10,&
&     '   Core    contribution to ann.rate with IPM (ns-1) =', lambda_core_ipm             *1000._dp,ch10,&
&     '   Valence contribution to ann.rate with IPM (ns-1) =',(lambda_ipm-lambda_core_ipm) *1000._dp
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
     if (dtset%usepaw==1) then
       write(msg, '(2a,6(2a,es16.8))' ) ch10,&
&       ' Annihilation rate PAW decomposition:',ch10,&
&       '   Plane-wave contribution to ann.rate          (ns-1) =',(lambda-lambda_paw)*1000._dp,ch10,&
&       '   Plane-wave valence contribution to ann.rate  (ns-1) =',(lambda-lambda_paw-lambda_core+lambda_core_paw)*1000._dp,ch10,&
&       '   On-site core contribution to ann.rate        (ns-1) =', lambda_core_paw*1000._dp,ch10,&
&       '   On-site valence contribution to ann.rate     (ns-1) =',(lambda_paw-lambda_core_paw)*1000._dp,ch10,&
&       '   Plane-wave contribution to ann.rate with IPM (ns-1) =',(lambda_ipm-lambda_paw_ipm)*1000._dp,ch10,&
&       '   Plane-wave core contrb. to ann.rate with IPM (ns-1) =',(lambda_core_ipm-lambda_core_paw_ipm)*1000._dp
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
     if (dtset%usepaw==0.and.igam==ngamma) then ! These tests are not relevant with PAW
       write(msg, '(2a,3(2a,es16.8))' ) ch10,&
&       ' ########## Some checks, for testing purpose:',ch10,&
&       '   Number of core electrons      =',nbec,ch10,&
&       '   Number of valence electrons   =',nbev,ch10,&
&       '   Number of positrons           =',nbp
       call wrtout(ab_out,msg,'COLL')
       call wrtout(std_out,msg,'COLL')
     end if
   end if !end if option
 end do ! Big loop on igam

 if (option==1) then
   write(msg, '(3a)' ) ch10,'      (*) IPM=Independent particle Model',ch10
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if !end if option

!Deallocate memory
 ABI_DEALLOCATE(igamma)
 if (dtset%usepaw==1.and.(.not.include_nhat_in_gamma)) then
   ABI_DEALLOCATE(rhor_)
   ABI_DEALLOCATE(rhor_ep_)
 end if

 DBG_EXIT("COLL")

end subroutine poslifetime
!!***
