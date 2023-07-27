!!****m* m_xc_tb09/m_xc_tb09
!! NAME
!!  m_xc_tb09
!!
!! FUNCTION
!!  This module contains routine(s) related to exchange-correlation Tran-Blaha 2009
!!    functional (modified Becke-Johnson)
!!
!! COPYRIGHT
!! Copyright (C) 2023-2023 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_xc_tb09

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_time,         only : timab
 use defs_abitypes,  only : MPI_type

 use m_xctk,         only : xcden
 use m_drivexc,      only : mkdenpos
 use m_geometry,     only : metric

 use m_pawang,       only : pawang_type
 use m_pawrad,       only : pawrad_type,pawrad_deducer0,nderiv_gen,simp_gen
 use m_pawtab,       only : pawtab_type
 use m_paw_an,       only : paw_an_type
 use m_pawrhoij,     only : pawrhoij_type
 use m_paw_denpot,   only : pawdensities
 use m_paral_atom,   only : get_my_atmtab,free_my_atmtab

 use libxc_functionals, only : libxc_functional_type,libxc_functionals_set_c_tb09, &
&                              libxc_functionals_is_tb09,libxc_functionals_ixc

 implicit none

 private

!public procedure(s)
 public :: xc_tb09_update_c ! Compute the C parameter for the TB09 functional (NCPP+PAW)

!Private variables
 real(dp),save :: grho_over_rho_pw    ! Contribution to TB09 c from the plane-wave pseudo-density
 real(dp),save :: grho_over_rho_paw   ! Contribution to TB09 c from the on-site PAW densities

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_xc_tb09/xc_tb09_update_c
!! NAME
!! xc_tb09_update_c
!!
!! FUNCTION
!! This routine computes from the density the C parameter of the Tran-Blaha 2009
!!  XC functional (modified Becke-Johnson)
!!  C = 1/V Int[|Grad(Rho)|/Rho] + Sum_PAW_aug_regions[ Natom/V [|Grad(Rho)|/Rho] ]
!!
!! INPUTS
!!  intxc = 1 if the XC functional has to be interpolated on a more refined mesh
!!  ixc= choice of exchange-correlation scheme
!!  mpi_enreg= information about MPI parallelization
!!  my_natom= number of atoms treated by current processor
!!  natom= total number of atoms in cell
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  ngfft(18)= contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhatdim= -PAW only- 0 if nhat array is not used ; 1 otherwise
!!  nhatgr(nfft,xcdata%nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nspden= number of spin-density components
!!  ntypat= number of types of atoms in unit cell.
!!  n3xccc= dimension of the xccc3d array (0 or nfft or cplx*nfft).
!!  paw_an(natom) <type(paw_an_type)>= paw arrays given on angular mesh
!!  pawang <type(pawang_type)> =paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)> =paw radial mesh and related data
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data (for the current atom)
!!  pawtab(ntypat) <type(pawtab_type)>= paw tabulated starting data
!!  pawxcdev= choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  rhor(nfft,nspden)= electron density in real space in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usepaw = 1 for PAW computation, 0 else
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xc_denpos= lowest allowed density
!!
!!  ----- Optional arguments -----
!!  [computation_type]= 'only_pw' : update only plane-wave contribution
!!                      'only_paw': update only PAW on-site contribution
!!                      'all' (default): compute all contributions  
!!  [mpi_atmtab(:)]= indexes of the atoms treated by current proc
!!  [comm_atom]= MPI communicator over atoms
!!  [xc_funcs(2)]= <type(libxc_functional_type)>=libxc XC functionals (if not the global one)
!!
!! OUTPUT
!!  No output
!!  The C parameter is directly set in the xc_functional datastructure
!!
!! SOURCE

subroutine xc_tb09_update_c(intxc,ixc,mpi_enreg,my_natom,natom,nfft,ngfft,nhat,nhatdim, &
&                           nhatgr,nhatgrdim,nspden,ntypat,n3xccc,paw_an,pawang,pawrad, &
&                           pawrhoij,pawtab,pawxcdev,rhor,rprimd,usepaw,xccc3d,xc_denpos, &
&                           computation_type,mpi_atmtab,comm_atom,xc_funcs) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intxc,ixc,n3xccc,nfft,my_natom,natom,nhatdim,nhatgrdim
 integer,intent(in) :: nspden,ntypat,pawxcdev,usepaw
 integer,intent(in),optional :: comm_atom
 real(dp),intent(in),optional :: xc_denpos
 character(len=*),intent(in),optional :: computation_type
 type(pawang_type),intent(in) :: pawang
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in),optional,target :: mpi_atmtab(:)
 real(dp),intent(in) :: nhat(nfft,nspden*nhatdim),nhatgr(nfft,nspden,3*nhatgrdim)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3),xccc3d(n3xccc)
 type(libxc_functional_type),intent(inout),optional :: xc_funcs(2)
 type(paw_an_type),intent(in) :: paw_an(my_natom*usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in),target :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,iatom_tot,ierr,ii,ilm,iloop,ir,itypat,ixc_from_lib,ipts
 integer :: ishift,iwarn,lm_size,lm_size_eff,mesh_size,my_comm_atom,ngrad,nfftot
 integer :: npts,nspden1,nzlmopt,option,opt_compch,opt_dens,pawprtvol,usecore
 integer :: usenhat,usexcnhat
 real(dp) :: compch_dum,factor,grho_over_rho,grho2,rho,sumg,ucvol,xc_c_tb09
 logical :: my_atmtab_allocated,paral_atom,test_tb09,use_exact_nhat_gradient
 logical :: compute_pw,compute_paw
 character(len=500) :: msg
!arrays
 integer,pointer :: my_atmtab(:)
 logical,allocatable :: lmselect(:),lmselect_tmp(:)
 real(dp) :: gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: dylmdr(:,:,:),ff(:)
 real(dp),allocatable :: rhotot(:),grhotot(:,:),rhonow(:,:,:)
 real(dp),allocatable :: rhoxc(:,:),drhoxc(:),dcorexc(:)
 real(dp),allocatable :: rho1(:,:,:),trho1(:,:,:),nhat1(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: corexc(:)

! *************************************************************************

!call timab(xx,1,tsec)

!Type of calculation
 compute_pw=.true. ; compute_paw=.true.
 if (present(computation_type)) then
   if (computation_type=='only_pw') compute_paw=.false.
   if (computation_type=='only_paw') compute_pw=.false.
 end if
 if ((.not.compute_pw).and.(.not.compute_paw)) return

!Nothing to do if XC functional is not TB09
 if (present(xc_funcs)) then
   test_tb09=libxc_functionals_is_tb09(xc_functionals=xc_funcs)
 else
   test_tb09=libxc_functionals_is_tb09()
 end if
 if (.not.test_tb09) return

!Check options
 if(ixc<0) then
   if (present(xc_funcs)) then
     ixc_from_lib=libxc_functionals_ixc(xc_functionals=xc_funcs)
   else
     ixc_from_lib=libxc_functionals_ixc()
   end if
   if (ixc /= ixc_from_lib) then
     write(msg, '(a,i0,a,a,i0)')&
&     'The value of ixc specified in input, ixc = ',ixc,ch10,&
&     'differs from the one used to initialize the functional ',ixc_from_lib
     ABI_BUG(msg)
   end if
 end if
 if (usepaw==1.and.compute_paw) then
   if (pawxcdev/=0) then
     msg='TB09 (mBJ) XC functional only available with pawxcdev=0!'
     ABI_BUG(msg)
   end if
   if(pawang%angl_size==0) then
     msg='pawang%angl_size=0!'
     ABI_BUG(msg)
   end if
   if(.not.allocated(pawang%ylmr)) then
     msg='pawang%ylmr must be allocated!'
     ABI_BUG(msg)
   end if
   if(.not.allocated(pawang%ylmrgr)) then
     msg='pawang%ylmrgr must be allocated!'
     ABI_BUG(msg)
   end if
 end if

!Initialize cell data 
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Other initializations
 usexcnhat=0; if (usepaw==1) usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)

!Initialization of Int[|Grad(Rho)|/Rho] contributions (plane-wave and on-site if PAW)
 if (compute_pw) grho_over_rho_pw=zero
 if (compute_paw) grho_over_rho_paw=zero
 
! ========================================================================
! =========== Plane wave contribution to Int[|Grad(Rho)|/Rho] ============

 if (compute_pw) then

!  Total density used in XC (with/without core density and compensation density)
   ABI_MALLOC(rhotot,(nfft))
   rhotot(1:nfft)=rhor(1:nfft,1)
   if (usexcnhat==0.and.nhatdim==1) rhotot(1:nfft)=rhotot(1:nfft)-nhat(1:nfft,1)
   if (n3xccc>0) rhotot(1:nfft)=rhotot(1:nfft)+xccc3d(1:nfft)

!  When possible, it is better to use the exact (analytical) gradient of nhat
!  In that case, we substract nhat from the density now
!  If not possible, the gradient will be computed in Fourier space (see call to xcden)
   use_exact_nhat_gradient=(nhatdim==1.and.nhatgrdim==1.and.usexcnhat==1.and.intxc==0)
   if (use_exact_nhat_gradient) rhotot(1:nfft)=rhotot(1:nfft)-nhat(1:nfft,1)
 
!  Loop on unshifted or shifted grids
   do ishift=0,intxc

!    Set up density on unshifted or shifted grid (will be in rhonow(:,:,1)),
!    as well as the gradient of the density, also on the unshifted
!    or shifted grid (will be in rhonow(:,:,2:4)), if needed.
     ABI_MALLOC(rhonow,(nfft,1,4))
     nspden1=1 ; cplex=1 ; ngrad=2 ; qphon(:)=zero
     call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden1,qphon,rhotot,rhonow)

!    PAW: if they were previously substracted, re-add nhat and its "exact" gradient
     if (use_exact_nhat_gradient) then
       rhonow(1:nfft,1,1)=rhonow(1:nfft,1,1)+nhat(1:nfft,1)
       do ii=1,3
         rhonow(1:nfft,1,ii+1)=rhonow(1:nfft,1,ii+1)+nhatgr(1:nfft,1,ii)
       end do
     end if

!    Make the density positive everywhere (but do not care about gradients)
     nspden1=1 ; iwarn=1  ; option=1
     call mkdenpos(iwarn,nfft,nspden1,option,rhonow(:,1,1),xc_denpos)

!    Accumulate integral of |Grad_rho|/Rho (to be used for TB09 XC)
     do ipts=1,nfft
       rho=rhonow(ipts,1,1)
       if (abs(rho)>tol10) then
         grho2=rhonow(ipts,1,2)**2+rhonow(ipts,1,3)**2+rhonow(ipts,1,4)**2
         grho_over_rho_pw=grho_over_rho_pw+sqrt(grho2)/rho
       end if
     end do

     ABI_FREE(rhonow)

!  End loop on unshifted or shifted grids
   end do

!  Normalize the integral
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
   grho_over_rho_pw=grho_over_rho_pw*ucvol/dble(nfftot)/dble(intxc+1)

!  Reduction in case of FFT distribution
   if (mpi_enreg%nproc_fft>1) then
     call xmpi_sum(grho_over_rho_pw,mpi_enreg%comm_fft,ierr)
   end if

   ABI_FREE(rhotot)

 end if ! compute_pw
 
! ========================================================================
! ========== PAW on-site contributions to Int[|Grad(Rho)|/Rho] ===========

 if (usepaw==1.and.compute_paw) then

!  Set up parallelism over atoms
   paral_atom=(present(comm_atom).and.(my_natom/=natom))
   nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
   my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
   call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

   usecore=1 ; nzlmopt=-1
   npts=pawang%angl_size

!  Compute Ylm gradients
!    dYlm/dr_i = { dYlm/dr_i^hat - r_i^hat * Sum_j[dYlm/dr_j^hat r_j^hat] } * (1/r)
   ABI_MALLOC(dylmdr,(3,npts,pawang%ylm_size))
   do ilm=1,pawang%ylm_size
     do ipts=1,npts
       factor=sum(pawang%ylmrgr(1:3,ilm,ipts)*pawang%anginit(1:3,ipts))
       dylmdr(1:3,ipts,ilm)=pawang%ylmrgr(1:3,ilm,ipts)-factor*pawang%anginit(1:3,ipts)
     end do
   end do

!  Loop on atom sites
   do iatom=1,my_natom
     iatom_tot=iatom;if (paral_atom) iatom_tot=my_atmtab(iatom)

     itypat=pawrhoij(iatom)%itypat
     mesh_size=pawtab(itypat)%mesh_size
     lm_size=paw_an(iatom)%lm_size
     lm_size_eff=min(lm_size,pawang%ylm_size)
     ABI_MALLOC(ff,(mesh_size))

!    Compute on-site densities
     ABI_MALLOC(rho1 ,(mesh_size,lm_size,1))
     ABI_MALLOC(trho1,(mesh_size,lm_size,1))
     ABI_MALLOC(nhat1,(mesh_size,lm_size,usexcnhat))
     ABI_MALLOC(lmselect,(lm_size))
     ABI_MALLOC(lmselect_tmp,(lm_size))
     lmselect_tmp(:)=.true. ; if (nzlmopt==1) lmselect_tmp(:)=paw_an(iatom)%lmselect(:)
     cplex=1 ; nspden1=1 ; opt_compch=0 ; opt_dens=1-usexcnhat ; pawprtvol=0
     call pawdensities(compch_dum,cplex,iatom_tot,lmselect_tmp,lmselect,&
&     lm_size,nhat1,nspden1,nzlmopt,opt_compch,opt_dens,-1,0,pawang,pawprtvol,&
&     pawrad(itypat),pawrhoij(iatom),pawtab(itypat),rho1,trho1)
     ABI_FREE(lmselect_tmp)

!    Allocation of temporary densities
     ABI_MALLOC(rhotot,(mesh_size))
     ABI_MALLOC(grhotot,(mesh_size,3))
     ABI_MALLOC(rhoxc,(mesh_size,lm_size))
     ABI_MALLOC(drhoxc,(mesh_size))
     ABI_MALLOC(dcorexc,(mesh_size))
     
!    Loop over AE and PS contributions
     do iloop=1,2
       if (iloop==1) then
         rhoxc(:,:)=rho1(:,:,1)
         corexc => pawtab(itypat)%coredens(:)
         usenhat=0
       else
         rhoxc(:,:)=trho1(:,:,1)
         corexc => pawtab(itypat)%tcoredens(:,1)
         usenhat=usexcnhat
       end if   

!      Add  hat density if needed
       if (usenhat>0) rhoxc(:,:)=rhoxc(:,:)+nhat1(:,:,1)

!      Need derivative of core density
       if (usecore==1) then
         call nderiv_gen(dcorexc,corexc,pawrad(itypat))
       end if

!      Do loop on the angular part (theta,phi)
       do ipts=1,npts

!        Compute the density and its gradient for this (theta,phi)
         rhotot(:)=zero ; grhotot(:,:)=zero
         do ilm=1,lm_size_eff
           if (lmselect(ilm)) then
             !Density
             rhotot(:)=rhotot(:)+rhoxc(:,ilm)*pawang%ylmr(ilm,ipts)
             !Gradient
             ff(1:mesh_size)=rhoxc(1:mesh_size,ilm)
             call nderiv_gen(drhoxc,ff,pawrad(itypat))
             ff(2:mesh_size)=ff(2:mesh_size)/pawrad(itypat)%rad(2:mesh_size)
             call pawrad_deducer0(ff,mesh_size,pawrad(itypat))
             do ii=1,3
               grhotot(1:mesh_size,ii)=grhotot(1:mesh_size,ii) &
&                +drhoxc(1:mesh_size)*pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts) &
&                +ff(1:mesh_size)*dylmdr(ii,ipts,ilm)
             end do
           end if
         end do
         if (usecore==1) then
           rhotot(:)=rhotot(:)+corexc(:)
           do ii=1,3
             grhotot(:,ii)=grhotot(:,ii)+dcorexc(:)*pawang%anginit(ii,ipts)
           end do
         end if

!        Make the density positive everywhere (but do not care about gradients)
         nspden1=1 ; iwarn=1
         call mkdenpos(iwarn,mesh_size,nspden1,0,rhotot,xc_denpos)

!        Accumulate integral of |Grad_rho|/Rho
         do ir=1,mesh_size
           rho=rhotot(ir)
           if (abs(rho)>tol10) then
             grho2=grhotot(ir,1)**2+grhotot(ir,2)**2+grhotot(ir,3)**2
             ff(ir)=sqrt(grho2)/rho
           end if
         end do
         ff(1:mesh_size)=ff(1:mesh_size)*pawrad(itypat)%rad(1:mesh_size)**2
         call simp_gen(sumg,ff,pawrad(itypat))

!        Add the contribution of the atom
         if (iloop==1) then
           grho_over_rho_paw=grho_over_rho_paw+sumg*four_pi*pawang%angwgth(ipts)
         else
           grho_over_rho_paw=grho_over_rho_paw-sumg*four_pi*pawang%angwgth(ipts)
         end if

       end do ! npts
     end do ! AE/PS loop

!    Deallocate temporary memory
     ABI_FREE(ff)
     ABI_FREE(rho1)
     ABI_FREE(trho1)
     ABI_FREE(nhat1)
     ABI_FREE(lmselect)
     ABI_FREE(rhotot)
     ABI_FREE(grhotot)
     ABI_FREE(rhoxc)
     ABI_FREE(drhoxc)
     ABI_FREE(dcorexc)

   end do ! atom sites

!  Deallocate temporary memory
   ABI_FREE(dylmdr)

!  Reduction in case of parallelism
   if (paral_atom) then
     call xmpi_sum(grho_over_rho_paw,my_comm_atom,ierr)
   end if

!  Destroy atom table used for parallelism
   call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 end if ! PAW and compute_paw
 
! ========================================================================
! ============ Assemble contributions and update mBJ C value  ============

 grho_over_rho = (grho_over_rho_pw + grho_over_rho_paw)/ucvol
 xc_c_tb09 = -0.012_dp + 1.023_dp * sqrt(grho_over_rho)
 if (xc_c_tb09<one) xc_c_tb09=one
 if (xc_c_tb09>two) xc_c_tb09=two

 write(msg,'(2a,f10.5,a)') ch10, &
&  ' In the mGGA functional TB09, c (computed value) = ',xc_c_tb09,ch10
 call wrtout(std_out,msg,'COLL')

 if (present(xc_funcs)) then
   call libxc_functionals_set_c_tb09(xc_c_tb09,xc_functionals=xc_funcs)
 else
   call libxc_functionals_set_c_tb09(xc_c_tb09)
 end if

!call timab(xx,2,tsec)

end subroutine xc_tb09_update_c
!!***

!----------------------------------------------------------------------

END MODULE m_xc_tb09
!!***
