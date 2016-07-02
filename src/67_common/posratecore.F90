!{\src2tex{textfont=tt}}
!!****f* ABINIT/posratecore
!! NAME
!! posratecore
!!
!! FUNCTION
!! Calculate the annihilataion rate of a given core state
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (JW,GJ,MT)
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
!!  iatom= index of the current atom in posdoppler
!!  mesh_sizej= size of the radial mesh for the current atom in posdoppler
!!  mpi_enreg= informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  option= if 1, use gamma
!!          if 2, use IPM (gamma=1)
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  rate= annihilation rate of a given core state needed for state dependent scheme for doppler broadening
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      posdoppler
!!
!! CHILDREN
!!      gammapositron,mkdenpos,nderiv_gen,pawdensities,pawxcsum,simp_gen
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine posratecore(dtset,electronpositron,iatom,my_natom,mesh_sizej,mpi_enreg,&
&                      option,pawang,pawrad,pawrhoij,pawrhoij_ep,&
&                      pawtab,rate,rhocorej)

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
#define ABI_FUNC 'posratecore'
 use interfaces_41_xc_lowlevel
 use interfaces_56_xc
 use interfaces_65_paw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,my_natom,option,mesh_sizej
 real(dp),intent(out) :: rate
 type(dataset_type), intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
!arrays
 real(dp),intent(in) :: rhocorej(mesh_sizej)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*dtset%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*dtset%usepaw)
 type(pawrhoij_type),intent(in),target :: pawrhoij_ep(my_natom*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,ierr,igamma,ii,ilm,ipt,ir
 integer :: itypat,iwarn,iwarnj,iwarnp,lm_size,lmn2_size,mesh_size
 integer :: ngr,ngrad,nspden_ep,opt_dens
 logical,parameter :: include_nhat_in_gamma=.false.
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,fact2,intg
 real(dp) :: mpibuf,rdum,sqfpi
 character(len=500) :: msg
!arrays
 logical,allocatable :: lmselect(:),lmselect_ep(:),lmselect_dum(:)
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/),lsign(2)=(/one,-one/)
 real(dp),allocatable :: d1gam(:,:),d2gam(:,:),ff(:),gam_(:,:,:),gamma(:,:),gammam(:,:),gg(:,:)
 real(dp),allocatable :: grhocore2(:),grhocor2_(:),grhoe2(:),grho2_(:)
 real(dp),allocatable :: nhat1(:,:,:),nhat1_ep(:,:,:)
 real(dp),allocatable :: rho_(:),rho_ep_(:),rho1(:,:,:),rho1_ep(:,:,:)
 real(dp),allocatable :: rhoarr1(:),rhoarr1_ep(:),rhoarr2(:)
 real(dp),allocatable :: rhocore(:),rhocor_(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),rhotot(:,:),rhotot_ep(:,:)
 real(dp),allocatable :: trho1(:,:,:),trho1_ep(:,:,:)
 real(dp),allocatable :: v1sum(:,:),v2sum(:,:,:)
 type(pawrhoij_type),pointer :: pawrhoij_ep_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Tests for developers
 if (.not.associated(electronpositron)) then
   msg='electronpositron variable must be associated!'
   MSG_BUG(msg)
 end if
!Constants
 fact=0.0
 cplex=1;nspden_ep=1
 ngrad=1;if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 iwarn=0;iwarnj=0;iwarnp=1
 sqfpi=sqrt(four_pi)

!Compatibility tests
 if (electronpositron%particle==EP_NOTHING) then
   msg='Not valid for electronpositron%particle=NOTHING!'
   MSG_BUG(msg)
 end if

 if (dtset%usepaw==1) then
   if(dtset%pawxcdev==0.and.ngrad==2) then
     msg='GGA is not implemented for pawxcdev=0 (use dtset%pawxcdev/=0)!'
     MSG_BUG(msg)
   end if
 end if

!Select type(s) of enhancement factor
 if (electronpositron%ixcpositron==-1) igamma=0
 if (electronpositron%ixcpositron== 2) igamma=4
 if (electronpositron%ixcpositron==11.or.electronpositron%ixcpositron==31) igamma=3
 if (electronpositron%ixcpositron==1.or.electronpositron%ixcpositron==3) igamma=2
 if (option==2) igamma=0

 pawrhoij_ep_ => pawrhoij_ep

 rate=zero

 itypat=pawrhoij(iatom)%itypat
 lmn2_size=pawtab(itypat)%lmn2_size
 mesh_size=pawtab(itypat)%mesh_size
 lm_size=pawtab(itypat)%lcut_size**2
 cplex=1
 ngr=0;if (ngrad==2) ngr=mesh_size

!Allocations of "on-site" densities
 ABI_ALLOCATE(rho1 ,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(trho1,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(rho1_ep ,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(trho1_ep,(cplex*mesh_size,lm_size,nspden_ep))
 ABI_ALLOCATE(lmselect,(lm_size))
 ABI_ALLOCATE(lmselect_ep,(lm_size))
 ABI_ALLOCATE(lmselect_dum,(lm_size))
 if (include_nhat_in_gamma) then
   ABI_ALLOCATE(nhat1,(cplex*mesh_size,lm_size,nspden_ep))
   ABI_ALLOCATE(nhat1_ep,(cplex*mesh_size,lm_size,nspden_ep))
 else
   ABI_ALLOCATE(nhat1,(0,0,0))
   ABI_ALLOCATE(nhat1_ep,(0,0,0))
 end if

!Compute "on-site" densities (n1, ntild1, nhat1) for electron and positron =====
 lmselect(:)=.true.
 opt_dens=1;if (include_nhat_in_gamma) opt_dens=0
 call pawdensities(rdum,cplex,iatom,lmselect,lmselect_dum,lm_size,nhat1,nspden_ep,1,&
& 0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij(iatom),&
& pawtab(itypat),rho1,trho1)
 lmselect_ep(:)=.true.
 call pawdensities(rdum,cplex,iatom,lmselect_ep,lmselect_dum,lm_size,nhat1_ep,nspden_ep,1,&
& 0,opt_dens,-1,0,pawang,0,pawrad(itypat),pawrhoij_ep_(iatom),&
& pawtab(itypat),rho1_ep,trho1_ep)
!Compute contribution to annihilation rate

 ABI_ALLOCATE(rhocore,(mesh_size))

!First formalism: use densities on r,theta,phi
 if (dtset%pawxcdev==0) then

   ABI_ALLOCATE(gamma,(mesh_size,2))
   ABI_ALLOCATE(rhoarr1,(mesh_size))
   ABI_ALLOCATE(rhoarr1_ep,(mesh_size))

!  Loop on the angular part
   do ipt=1,pawang%angl_size
!    Build densities
     rhoarr1=zero;rhoarr1_ep=zero;rhocore=zero
     do ilm=1,lm_size
       if (lmselect(ilm)) rhoarr1(:)=rhoarr1(:)+rho1(:,ilm,1)*pawang%ylmr(ilm,ipt)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) rhoarr1_ep(:)=rhoarr1_ep(:)+rho1_ep(:,ilm,1)*pawang%ylmr(ilm,ipt)
     end do
     rhocore(:)=pawtab(itypat)%coredens(:)
!    Make the densities positive
     if (electronpositron%particle==EP_ELECTRON) then
       call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
       call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
     else if (electronpositron%particle==EP_POSITRON) then
       call mkdenpos(iwarn ,mesh_size,1,1,rhoarr1   ,dtset%xc_denpos)
       call mkdenpos(iwarnp,mesh_size,1,1,rhoarr1_ep,dtset%xc_denpos)
     end if
!    Compute Gamma
     ABI_ALLOCATE(grhoe2,(ngr))
     ABI_ALLOCATE(grhocore2,(ngr))
     if (electronpositron%particle==EP_ELECTRON) then
       call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,mesh_size,&
&       rhocore,rhoarr1_ep,rhoarr1,1)
     else if (electronpositron%particle==EP_POSITRON) then
       call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,mesh_size,&
&       rhocore,rhoarr1,rhoarr1_ep,1)
     end if
     ABI_DEALLOCATE(grhoe2)
     ABI_DEALLOCATE(grhocore2)
!    Compute contribution to annihilation rates

     ABI_ALLOCATE(ff,(mesh_size))
     ff(:)=rhoarr1_ep(:)*rhocorej(:)*gamma(:,1)*pawrad(itypat)%rad(:)**2
     call simp_gen(intg,ff,pawrad(itypat))
     intg=intg*pawang%angwgth(ipt)*four_pi
     rate         =rate         +intg
     ABI_DEALLOCATE(ff)
   end do ! ipt
   ABI_DEALLOCATE(gamma)
   ABI_DEALLOCATE(rhoarr1)
   ABI_DEALLOCATE(rhoarr1_ep)

!Second formalism: use (l,m) moments for densities
 else if (dtset%pawxcdev/=0) then

!  Build densities
   ABI_ALLOCATE(gammam,(mesh_size,lm_size))
   ABI_ALLOCATE(rhotot,(mesh_size,lm_size))
   ABI_ALLOCATE(rhotot_ep,(mesh_size,lm_size))
   ABI_ALLOCATE(rhosph,(mesh_size))
   ABI_ALLOCATE(rhosph_ep,(mesh_size))

   rhotot   (:,:)=rho1   (:,:,1)
   rhotot_ep(:,:)=rho1_ep(:,:,1)
   rhocore(:)=pawtab(itypat)%coredens(:)
   rhosph   (:)=rhotot   (:,1)/sqfpi
   rhosph_ep(:)=rhotot_ep(:,1)/sqfpi
!  Make spherical densities positive
   if (electronpositron%particle==EP_ELECTRON) then
     call mkdenpos(iwarnp,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
     call mkdenpos(iwarn ,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
   else if (electronpositron%particle==EP_POSITRON) then
     call mkdenpos(iwarn ,mesh_size,1,1,rhosph   ,dtset%xc_denpos)
     call mkdenpos(iwarnp,mesh_size,1,1,rhosph_ep,dtset%xc_denpos)
   end if

!  Need gradients of electronic densities for GGA
   ABI_ALLOCATE(grhoe2,(ngr))
   ABI_ALLOCATE(grhocore2,(ngr))
   if (ngr>0) then
     if (electronpositron%particle==EP_ELECTRON) then
       call nderiv_gen(grhoe2,rhosph_ep,pawrad(itypat))
     else if (electronpositron%particle==EP_POSITRON) then
       call nderiv_gen(grhoe2,rhosph,pawrad(itypat))
     end if
     grhoe2(:)=grhoe2(:)**2
     call nderiv_gen(grhocore2,rhocore,pawrad(itypat))
     grhocore2(:)=grhocore2(:)**2
   end if
!  Compute Gamma for (rho-,rho+),
!  (rho- +drho-,rho+), (rho- -drho-,rho+),
!  (rho-,rho+ +drho+), (rho-,rho+ -drho+),
!  (rho- +drho-,rho+ +drho+), (rho- -drho-,rho+ -drho+)
!  Do a seven steps loop
   ABI_ALLOCATE(gam_,(mesh_size,2,7))
   ABI_ALLOCATE(rho_,(mesh_size))
   ABI_ALLOCATE(rho_ep_,(mesh_size))
   ABI_ALLOCATE(rhocor_,(mesh_size))
   ABI_ALLOCATE(grho2_,(ngr))
   ABI_ALLOCATE(grhocor2_,(ngr))

   do ii=1,7
!    Apply delta to get perturbed densities
     rho_(:)=rhosph(:);rho_ep_(:)=rhosph_ep(:);rhocor_(:)=rhocore(:)
     if (ngr>0) grho2_(:)=grhoe2(:)
     if (ngr>0) grhocor2_(:)=grhocore2(:)
     if (ii==2.or.ii==4.or.ii==6) fact=(one+delta)
     if (ii==3.or.ii==5.or.ii==7) fact=(one-delta)
     fact2=fact**2
     if (ii==2.or.ii==3.or.ii==6.or.ii==7) then
       rho_(:)=fact*rho_(:)
       if (electronpositron%particle==EP_POSITRON) then
         if (ngr>0) grho2_(:)=fact2*grho2_(:)
         rhocor_(:)=fact*rhocor_(:)
         if (ngr>0) grhocor2_(:)=fact2*grhocor2_(:)
       end if
     end if

     if (ii==4.or.ii==5.or.ii==6.or.ii==7) then
       rho_ep_(:)=fact*rho_ep_(:)
       if (electronpositron%particle==EP_ELECTRON) then
         if (ngr>0) grho2_(:)=fact2*grho2_(:)
         rhocor_(:)=fact*rhocor_(:)
         if (ngr>0) grhocor2_(:)=fact2*grhocor2_(:)
       end if
     end if
!    Compute gamma for these perturbed densities
     if (electronpositron%particle==EP_ELECTRON) then
       call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma,ngr,mesh_size,rhocor_,rho_ep_,rho_,1)
     else if (electronpositron%particle==EP_POSITRON) then
       call gammapositron(gam_(:,:,ii),grhocor2_,grho2_,igamma,ngr,mesh_size,rhocor_,rho_,rho_ep_,1)
     end if

   end do ! end loop ii=1,7

   ABI_DEALLOCATE(rhocor_)
   ABI_DEALLOCATE(grho2_)
   ABI_DEALLOCATE(grhocor2_)
   ABI_DEALLOCATE(grhoe2)
   ABI_DEALLOCATE(grhocore2)
   rho_   (:)=rhosph   (:);if (electronpositron%particle==EP_POSITRON) rho_   (:)=rho_   (:)+rhocore(:)
   rho_ep_(:)=rhosph_ep(:);if (electronpositron%particle==EP_ELECTRON) rho_ep_(:)=rho_ep_(:)+rhocore(:)
!  Compute numerical first and second derivatives of Gamma
!  d1gam(1) = dgam/drho+ (particle=ELECTRON), dgam/drho- (particle=POSITRON)
!  d1gam(2) = dgam/drho- (particle=ELECTRON), dgam/drho+ (particle=POSITRON)
   ABI_ALLOCATE(d1gam,(mesh_size,2))
   d1gam(:,:)=zero
   do ir=1,mesh_size
     if (rho_     (ir)>tol14) d1gam(ir,1)=(gam_(ir,1,2)-gam_(ir,1,3))*half/(delta*rho_     (ir))
     if (rho_ep_  (ir)>tol14) d1gam(ir,2)=(gam_(ir,1,4)-gam_(ir,1,5))*half/(delta*rho_ep_  (ir))
   end do

!  d2gam(1) = d2gam/drho+_drho+ (particle=ELECTRON), dgam/drho-_drho- (particle=POSITRON)
!  d2gam(2) = d2gam/drho-_drho+ (particle=ELECTRON), dgam/drho+_drho- (particle=POSITRON)
!  d2gam(3) = d2gam/drho-_drho- (particle=ELECTRON), dgam/drho+_drho+ (particle=POSITRON)
   ABI_ALLOCATE(d2gam,(mesh_size,3))
   d2gam(:,:)=zero
   do ir=1,mesh_size
     if (rho_  (ir)>tol14) d2gam(ir,1)=(gam_(ir,1,2)+gam_(ir,1,3)-two*gam_(ir,1,1))/(delta*rho_  (ir))**2
     if (rho_ep_(ir)>tol14) then
       d2gam(ir,3)=(gam_(ir,1,4)+gam_(ir,1,5)-two*gam_(ir,1,1))/(delta*rho_ep_(ir))**2
       if (rho_(ir)>tol14) then
         d2gam(ir,2)=(gam_(ir,1,6)+gam_(ir,1,7)+two*gam_(ir,1,1) &
&         -gam_(ir,1,2)-gam_(ir,1,3)-gam_(ir,1,4)-gam_(ir,1,5)) &
&         *half/(delta*rho_(ir))/(delta*rho_ep_(ir))
       end if
     end if
   end do

   ABI_DEALLOCATE(rho_)
   ABI_DEALLOCATE(rho_ep_)
!  Compute useful sums of densities
   ABI_ALLOCATE(v1sum,(mesh_size,3))
   if ( dtset%pawxcdev>=2)  then
     ABI_ALLOCATE(v2sum,(mesh_size,lm_size,3))
   else
     ABI_ALLOCATE(v2sum,(0,0,0))
   end if
   rhotot(:,1)=sqfpi*rhosph(:);rhotot_ep(:,1)=sqfpi*rhosph_ep(:)
   call pawxcsum(1,1,1,lmselect,lmselect_ep,lm_size,mesh_size,3,dtset%pawxcdev,&
&   pawang,rhotot,rhotot_ep,v1sum,v2sum)
!  Compute final developpment of gamma moments
   gammam(:,:)=zero
   gammam(:,1)=gam_(:,1,1)*sqfpi
   gammam(:,1)=gammam(:,1)+(d2gam(:,2)*v1sum(:,2) &
&   +half*(d2gam(:,1)*v1sum(:,1)+d2gam(:,3)*v1sum(:,3)))/sqfpi
   do ilm=2,lm_size
     if (lmselect(ilm)) then
       gammam(:,ilm)=gammam(:,ilm)+d1gam(:,1)*rhotot(:,ilm)
     end if
     if (lmselect_ep(ilm)) then
       gammam(:,ilm)=gammam(:,ilm)+d1gam(:,2)*rhotot_ep(:,ilm)
     end if
   end do
   if (dtset%pawxcdev>1) then
     do ilm=2,lm_size
       gammam(:,ilm)=gammam(:,ilm)+d2gam(:,2)*v2sum(:,ilm,2) &
&       +half*(d2gam(:,1)*v2sum(:,ilm,1)+d2gam(:,3)*v2sum(:,ilm,3))
     end do
   end if

   ABI_DEALLOCATE(gam_)
   ABI_DEALLOCATE(d1gam)
   ABI_DEALLOCATE(d2gam)
   ABI_DEALLOCATE(v1sum)
   ABI_DEALLOCATE(v2sum)
!  Compute contribution to annihilation rate
   ABI_ALLOCATE(gg,(mesh_size,4))
   gg=zero
   ABI_ALLOCATE(rhoarr2,(mesh_size))
   do ilm=1,lm_size
     if (lmselect_ep(ilm)) gg(:,1)=gg(:,1)+rhotot_ep(:,ilm)*rhocorej(:)*gammam(:,ilm)
   end do
   ABI_DEALLOCATE(rhoarr2)
   gg(:,1)=gg(:,1)*pawrad(itypat)%rad(:)**2
   call simp_gen(intg,gg(:,1),pawrad(itypat))
   rate         =rate         +intg
   ABI_DEALLOCATE(gg)
   ABI_DEALLOCATE(gammam)
   ABI_DEALLOCATE(rhotot)
   ABI_DEALLOCATE(rhotot_ep)
   ABI_DEALLOCATE(rhosph)
   ABI_DEALLOCATE(rhosph_ep)

 end if ! dtset%pawxcdev
 ABI_DEALLOCATE(rhocore)

 ABI_DEALLOCATE(rho1)
 ABI_DEALLOCATE(trho1)
 ABI_DEALLOCATE(rho1_ep)
 ABI_DEALLOCATE(trho1_ep)
 ABI_DEALLOCATE(lmselect)
 ABI_DEALLOCATE(lmselect_ep)
 ABI_DEALLOCATE(lmselect_dum)
 ABI_DEALLOCATE(nhat1)
 ABI_DEALLOCATE(nhat1_ep)

!Reduction in case of distribution over atomic sites
 if (mpi_enreg%nproc_atom>1) then
   mpibuf=rate
   call xmpi_sum(mpibuf,mpi_enreg%comm_atom,ierr)
   rate=mpibuf
 end if



 DBG_EXIT("COLL")

end subroutine posratecore
!!***
