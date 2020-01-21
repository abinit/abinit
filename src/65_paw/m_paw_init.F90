!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_init
!! NAME
!!  m_paw_init
!!
!! FUNCTION
!!  This module contains routines related tp PAW calculations initialization.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2019 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_init

 use defs_basis
 use m_errors
 use m_abicore
 use m_splines
 use m_dtset

 use m_time,    only : timab
 use m_pawpsp,  only : pawpsp_nl
 use m_paw_atom,only : atompaw_shpfun
 use m_pawang,  only : pawang_type, pawang_init, pawang_free
 use m_pawrad,  only : pawrad_type, simp_gen, nderiv_gen, poisson, pawrad_deducer0
 use m_pawtab,  only : pawtab_type
 use m_paw_numeric, only : paw_derfc

 implicit none

 private

!public procedures.
 public :: pawinit     ! Initialize some tabulated data for PAW calculations
 public :: paw_gencond ! Test whether we have to call pawinit to regenerate tabulated data.

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_init/pawinit
!! NAME
!! pawinit
!!
!! FUNCTION
!! Initialize some starting values of several arrays used in PAW calculations.
!!
!! 1-Initialize data related to angular mesh
!! 2-Tabulate normalized shape function g(r)
!! 3-Compute indklmn indexes giving some l,m,n,lp,mp,np info
!!                           from klmn=[(l,m,n),(lp,mp,np)]
!! 4-Compute various factors/sizes (depending on (l,m,n))
!! 5-Compute $q_ijL=\displaystyle
!!                  \int_{0}^{r_c}{(\phi_i\phi_j-\widetilde{\phi_i}\widetilde{\phi_j}) r^l\,dr}
!!                   Gaunt(l_i m_i,l_j m_j,l m))$
!!           $S_ij=\displaystyle \sqrt{4 \pi} q_ij0$
!! 6-Compute $e_ijkl= vh1_ijkl - Vhatijkl - Bijkl - Cijkl$
!!     With:
!!       $vh1_ijkl =\sum_{L,m} {vh1*Gaunt(i,j,Lm)*Gaunt(k,l,Lm)}$
!!       $Vhat_ijkl=\sum_{L,m} {vhatijL*Gaunt(i,j,Lm)*q_klL}$
!!       $B_ijkl   =\sum_{L,m} {vhatijL*Gaunt(k,l,Lm)*q_ijL}$
!!       $C_ijkl   =\sum_{L,m} {intvhatL*q_ijL*q_klL}$
!!     and:
!!       vh1 according to eq. (A17) in Holzwarth et al., PRB 55, 2005 (1997) [[cite:Holzwarth1997]]
!! 7-Compute Ex-correlation energy for the core density
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!             also determine the size of these pointers
!!  gsqcut_shp=effective cut-off to determine shape functions in reciprocal space
!!  hyb_range_fock=range coefficient for screened hybrid XC functionals
!!  lcutdens=max. l for densities/potentials moments computations
!!  lmix=max. l for which spherical terms will be mixed durinf SCF cycle
!!  mpsang=1+maximum angular momentum
!!  nphi="phi" dimension of paw angular mesh
!!  nsym=Number of symmetry elements in space group
!!  ntheta="theta" dimension of paw angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)=
!!                     selection rules for Gaunt coefficients
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!     %radfact(mesh_size)=Factor used to compute radial integrals
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %basis_size=Number of elements for the PAW nl basis
!!     %l_size=Maximum value of l+1 leading to non zero Gaunt coeffs
!!     %lmn_size=Number of (l,m,n) elements for the PAW basis
!!     %lmn2_size=lmn_size*(lmn_size+1)/2
!!     %dltij(lmn2_size)=factors used to compute sums over (ilmn,jlmn)
!!     %phi(mesh_size,basis_size)=PAW all electron wavefunctions
!!     %rshp=shape function radius (radius for compensation charge)
!!     %shape_type=Radial shape function type
!!     %shape_alpha=Alpha parameters in Bessel shape function
!!     %shape_lambda=Lambda parameter in gaussian shape function
!!     %shape_q=Q parameters in Bessel shape function
!!     %shape_sigma=Sigma parameter in gaussian shape function
!!     %tphi(mesh_size,basis_size)=PAW atomic pseudowavefunctions
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1=dev. on moments)
!!  usekden= 1 is kinetic energy density has to be computed, 0 otherwise
!!  xclevel=XC functional level (1=LDA, 2=GGA)
!!
!! OUTPUT
!!  pawang
!!     %gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)=selection rules for Gaunt coefficients
!!     %l_max=maximum value of angular momentum l+1
!!     %l_size_max=maximum value of angular momentum l_size=2*l_max-1
!!     %nsym=number of symmetry elements in space group
!!     %ngnt=number of non-zero Gaunt coefficients
!!     %realgnt(pawang%ngnt)=non-zero real Gaunt coefficients
!!     === only if pawxcdev==1 ==
!!       %anginit(3,angl_size)=for each point of the angular mesh, gives the coordinates
!!                             of the corresponding point on an unitary sphere
!!       %angl_size=dimension of paw angular mesh (angl_size=ntheta*nphi)
!!       %angwgth(angl_size)=for each point of the angular mesh, gives the weight
!!                           of the corresponding point on an unitary sphere
!!       %ntheta, nphi=dimensions of paw angular mesh
!!       %ylmr(l_size_max**2,angl_size)=real Ylm calculated in real space
!!       %ylmrgr(1:3,l_size_max**2,angl_size)=first gradients of real Ylm calculated in real space
!!     === only if pawspnorb==1 ==
!!     %ls_ylm(2,l_max**2,l_max**2,4)=LS operator in the real spherical harmonics basis
!!     %use_ls_ylm=flag activated if ls_ylm is allocated
!!     %usespnorb=flag activated for spin-orbit coupling
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated data read at start:
!!     %lcut_size_=max. value of l+1 leading to non zero Gaunt coeffs modified by lcutdens
!!     %lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
!!     %mqgrid_shp=number of points in reciprocal space for shape function
!!     %indklmn(8,lmn2_size)=array giving klm, kln, abs(il-jl) and (il+jl), ilmn and jlmn for each klmn=(ilmn,jlmn)
!!     %dshpfunc(mesh_size,l_size,4)=derivatives of shape function (used only for numerical shape functions)
!!     %eijkl(lmn2_size,lmn2_size)=part of the Dij that depends only from the projected occupation coeffs
!!     %exccore=Exchange-correlation energy for the core density
!!     %gnorm(l_size)=normalization factor of radial shape function
!!     %phiphj(:,:)=useful product Phi(:,i)*Phi(:,j)
!!     %qgrid_shp(mqgrid_shp)=points in reciprocal space for shape function
!!     %qijl(l_size**2,lmn2_size)=moments of the difference charge density between AE and PS partial wave
!!     %rad_for_spline(mesh_size)=radial grid used for spline (copy of pawrad%rad)
!!     %shapefunc(mesh_size,l_size)=normalized radial shape function
!!     %shapefncg(mqgrid_shp,l_size)=normalized radial shape function in reciprocal space
!!     %sij(lmn2_size)=nonlocal part of the overlap operator
!!     %tphitphj(:,:)=useful product tPhi(:,i)*tPhi(:,j)
!!
!! PARENTS
!!      bethe_salpeter,gstate,respfn,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawinit(effmass_free,gnt_option,gsqcut_eff,hyb_range_fock,lcutdens,lmix,mpsang,&
&                  nphi,nsym,ntheta,pawang,pawrad,pawspnorb,pawtab,pawxcdev,xclevel,&
&                  usekden,usepotzero)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: gnt_option,lcutdens,lmix,mpsang,nphi,nsym,ntheta,usekden
 integer,intent(in) :: pawspnorb,pawxcdev,xclevel,usepotzero
 real(dp),intent(in) :: effmass_free,gsqcut_eff,hyb_range_fock
 type(pawang_type),intent(inout) :: pawang
!arrays
 type(pawrad_type),intent(in) :: pawrad(:)
 type(pawtab_type),target,intent(inout) :: pawtab(:)

!Local variables ------------------------------
!scalars
 integer,parameter :: mqgrid_shp_default=300
 integer :: basis_size,i0lm,i0ln,ij_size,il,ilm,ilmn,iln,iloop,iq,isel,isel1
 integer :: itypat,j0lm,j0lmn,j0ln,jl,jlm,jlmn,jln,klm,klm1
 integer :: klmn,klmn1,kln,kln1,l_size,ll,lm0,lmax,lmax1,lmin,lmin1,lmn2_size
 integer :: lmn_size,lmnmix,mesh_size,meshsz,mm,nabgnt_option,ntypat,pw_mesh_size,usexcnhat,use_ls_ylm,use_ylm
 real(dp) :: dq,gnrm,intg,ql,ql1,rg,rg1,vh1,yp1,ypn
 character(len=500) :: message
!arrays
 integer,allocatable :: indl(:,:),klm_diag(:),kmix_tmp(:)
 integer, ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: der(:),ff(:),gg(:),hh(:),indklmn_(:,:),intvhatl(:)
 real(dp),allocatable :: rad(:),rgl(:,:),vhatijl(:,:),vhatl(:),work(:)
 real(dp),pointer :: eijkl(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 call timab(553,1,tsec)
! if kinetic energy density is used, set nabgnt_option to 1 for nablagaunt computation
 nabgnt_option=0 ; if (usekden>=1) nabgnt_option=1
 ntypat=size(pawtab)
 if (size(pawrad)/=ntypat) then
   MSG_BUG('pawrad and pawtab should have the same size!')
 end if

 ! Immediately set the value of usepotzero
 ! it will be used later on in this subroutine
 pawtab%usepotzero=usepotzero

!==================================================
!1- INITIALIZE DATA RELATED TO ANGULAR MESH
!* ANGULAR GRID
!* REAL SPHERICAL HARMONICS
!* REAL GAUNT COEFFICIENTS

 use_ylm=0;if (pawxcdev==0) use_ylm=1
 use_ls_ylm=0;if (pawspnorb>0) use_ls_ylm=1
 call pawang_free(pawang)
 call pawang_init(pawang,gnt_option,nabgnt_option,usekden,mpsang-1,nphi,nsym,ntheta,pawxcdev,use_ls_ylm,use_ylm,xclevel)
 usexcnhat=maxval(pawtab(1:ntypat)%usexcnhat)

!*******************
!Loop on atom types
!*******************
 do itypat=1,ntypat
   mesh_size=pawtab(itypat)%mesh_size
   l_size=pawtab(itypat)%l_size
   lmn_size=pawtab(itypat)%lmn_size
   lmn2_size=pawtab(itypat)%lmn2_size
   basis_size=pawtab(itypat)%basis_size
   ij_size=pawtab(itypat)%ij_size
   indlmn => pawtab(itypat)%indlmn(:,:)
   ABI_ALLOCATE(indklmn_,(8,lmn2_size))
   ABI_ALLOCATE(klm_diag,(lmn2_size))
   ABI_ALLOCATE(ff,(mesh_size))
   ABI_ALLOCATE(gg,(mesh_size))
   ABI_ALLOCATE(hh,(mesh_size))
   ABI_ALLOCATE(rad,(mesh_size))
   rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)

   if (pawtab(itypat)%usexcnhat/=usexcnhat) then
     write(message, '(7a)' )&
&     'You cannot simultaneously use atomic data with different',ch10,&
&     'formulation of XC [using compensation charge in XC or not] !',ch10,&
&     'Action: change at least one of your atomic data (psp) file',ch10,&
&     '        or use usexcnhat keyword in input file.'
     MSG_ERROR(message)
   end if

!  ==================================================
!  2- TABULATE SHAPE FUNCTION

!  Allocated shape function
   if (pawtab(itypat)%shape_type/=-1) then
     if (allocated(pawtab(itypat)%shapefunc))  then
       ABI_DEALLOCATE(pawtab(itypat)%shapefunc)
     end if
     ABI_ALLOCATE(pawtab(itypat)%shapefunc,(mesh_size,l_size))
   else if (.not.allocated(pawtab(itypat)%shapefunc))  then
     message='shapefunc should be allocated with shape_type=-1'
     MSG_ERROR(message)
   end if
   if (allocated(pawtab(itypat)%gnorm))  then
     ABI_DEALLOCATE(pawtab(itypat)%gnorm)
   end if
   ABI_ALLOCATE(pawtab(itypat)%gnorm,(l_size))

!  Compute shape function
   do il=1,l_size
     ll=il-1
     call atompaw_shpfun(ll,pawrad(itypat),gnrm,pawtab(itypat),ff)
     pawtab(itypat)%shapefunc(1:mesh_size,il)=ff(1:mesh_size)
     pawtab(itypat)%gnorm(il)=gnrm
   end do
!  In case of numerical shape function, compute some derivatives
   if (pawtab(itypat)%shape_type==-1) then
     if (allocated(pawtab(itypat)%dshpfunc))  then
       ABI_DEALLOCATE(pawtab(itypat)%dshpfunc)
     end if
     ABI_ALLOCATE(pawtab(itypat)%dshpfunc,(mesh_size,l_size,4))
     ABI_ALLOCATE(work,(mesh_size))
     do il=1,l_size
       call nderiv_gen(pawtab(itypat)%dshpfunc(:,il,1),pawtab(itypat)%shapefunc(:,il),pawrad(itypat))
       yp1=pawtab(itypat)%dshpfunc(1,il,1);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,1)
       call spline(rad,pawtab(itypat)%shapefunc(:,il),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,2))
       yp1=pawtab(itypat)%dshpfunc(1,il,2);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,2)
       call spline(rad,pawtab(itypat)%dshpfunc(:,il,1),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,3))
       yp1=pawtab(itypat)%dshpfunc(1,il,3);ypn=pawtab(itypat)%dshpfunc(mesh_size,il,3)
       call spline(rad,pawtab(itypat)%dshpfunc(:,il,2),mesh_size,yp1,ypn,pawtab(itypat)%dshpfunc(:,il,4))
     end do
     ABI_DEALLOCATE(work)
   end if

!  In some cases, has to store radial mesh for shape function in pawtab variable
   if (pawtab(itypat)%shape_type==-1) then
     if (allocated(pawtab(itypat)%rad_for_spline))  then
       ABI_DEALLOCATE(pawtab(itypat)%rad_for_spline)
     end if
     ABI_ALLOCATE(pawtab(itypat)%rad_for_spline,(mesh_size))
     pawtab(itypat)%rad_for_spline(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)
   end if

!  In some cases, has to store shape function in reciprocal space
   if (pawtab(itypat)%has_shapefncg>0) then
     if (gsqcut_eff<tol8) then
       message='Computation of shapefncg only possible when gsqcut>0!'
       MSG_BUG(message)
     end if
     pawtab(itypat)%mqgrid_shp=mqgrid_shp_default
     if (allocated(pawtab(itypat)%shapefncg))  then
       ABI_DEALLOCATE(pawtab(itypat)%shapefncg)
     end if
     if (allocated(pawtab(itypat)%qgrid_shp))  then
       ABI_DEALLOCATE(pawtab(itypat)%qgrid_shp)
     end if
     ABI_ALLOCATE(pawtab(itypat)%shapefncg,(pawtab(itypat)%mqgrid_shp,2,l_size))
     ABI_ALLOCATE(pawtab(itypat)%qgrid_shp,(pawtab(itypat)%mqgrid_shp))
     dq=1.1_dp*sqrt(gsqcut_eff)/dble(pawtab(itypat)%mqgrid_shp-1)
     do iq=1,pawtab(itypat)%mqgrid_shp
       pawtab(itypat)%qgrid_shp(iq)=dble(iq-1)*dq
     end do
     ABI_ALLOCATE(indl,(6,l_size))
     ABI_ALLOCATE(rgl,(mesh_size,il))
     do il=1,l_size
       indl(:,il)=0;indl(1,il)=il-1;indl(5,il)=il
       rgl(1:mesh_size,il)=rad(1:mesh_size)*pawtab(itypat)%shapefunc(1:mesh_size,il)
     end do
     call pawpsp_nl(pawtab(itypat)%shapefncg,indl,l_size,l_size,&
&     pawtab(itypat)%mqgrid_shp,pawtab(itypat)%qgrid_shp,pawrad(itypat),rgl)
     pawtab(itypat)%shapefncg=four_pi*pawtab(itypat)%shapefncg
     ABI_DEALLOCATE(indl)
     ABI_DEALLOCATE(rgl)
   else
     pawtab(itypat)%mqgrid_shp=0
   end if

!  ==================================================
!  3- COMPUTE indklmn INDEXES GIVING klm, kln, abs(il-jl) and (il+jl), ilmn and jlmn
!  for each klmn=(ilmn,jlmn)

   if (allocated(pawtab(itypat)%indklmn))  then
     ABI_DEALLOCATE(pawtab(itypat)%indklmn)
   end if
   ABI_ALLOCATE(pawtab(itypat)%indklmn,(8,lmn2_size))

   klm_diag=0
   do jlmn=1,lmn_size
     jl= indlmn(1,jlmn);jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2
     j0lm =jlm *(jlm -1)/2
     j0ln =jln *(jln -1)/2
     do ilmn=1,jlmn
       il= indlmn(1,ilmn);ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
       klmn=j0lmn+ilmn
       if (ilm<=jlm) then
         indklmn_(1,klmn)=j0lm+ilm
       else
         i0lm=ilm*(ilm-1)/2
         indklmn_(1,klmn)=i0lm+jlm
       end if
       if (iln<=jln) then
         indklmn_(2,klmn)=j0ln+iln
       else
         i0ln=iln*(iln-1)/2
         indklmn_(2,klmn)=i0ln+jln
       end if
       indklmn_(3,klmn)=min(abs(il-jl),lcutdens)
       indklmn_(4,klmn)=min(il+jl,lcutdens)
       indklmn_(5,klmn)=ilm
       indklmn_(6,klmn)=jlm
       indklmn_(7,klmn)=ilmn
       indklmn_(8,klmn)=jlmn
       pawtab(itypat)%indklmn(:,klmn)=indklmn_(:,klmn)
       if (ilm==jlm) klm_diag(klmn)=1
     end do
   end do

!  ==================================================
!  4- COMPUTE various FACTORS/SIZES (depending on (l,m,n))

   pawtab(itypat)%usespnorb=pawspnorb
   pawtab(itypat)%lcut_size=min(l_size,lcutdens+1)

   if (allocated(pawtab(itypat)%dltij))  then
     ABI_DEALLOCATE(pawtab(itypat)%dltij)
   end if
   ABI_ALLOCATE(pawtab(itypat)%dltij,(lmn2_size))
   pawtab(itypat)%dltij(:)=two
   do ilmn=1,lmn_size
     pawtab(itypat)%dltij(ilmn*(ilmn+1)/2)=one
   end do

   lmnmix=zero
   ABI_ALLOCATE(kmix_tmp,(lmn2_size))
   do jlmn=1,lmn_size
     jl=indlmn(1,jlmn)
     if (jl<=lmix) then
       j0lmn=jlmn*(jlmn-1)/2
       do ilmn=1,jlmn
         il=indlmn(1,ilmn)
         if (il<=lmix) then
           lmnmix=lmnmix+1
           kmix_tmp(lmnmix)=j0lmn+ilmn
         end if
       end do
     end if
   end do
   if (allocated(pawtab(itypat)%kmix))  then
     ABI_DEALLOCATE(pawtab(itypat)%kmix)
   end if
   ABI_ALLOCATE(pawtab(itypat)%kmix,(lmnmix))
   pawtab(itypat)%lmnmix_sz=lmnmix
   pawtab(itypat)%kmix(1:lmnmix)=kmix_tmp(1:lmnmix)
   ABI_DEALLOCATE(kmix_tmp)

!  ==================================================
!  5- STORE SOME USEFUL QUANTITIES FROM PARTIAL WAVES

   if (allocated(pawtab(itypat)%phiphj))  then
     ABI_DEALLOCATE(pawtab(itypat)%phiphj)
   end if
   if (allocated(pawtab(itypat)%tphitphj))  then
     ABI_DEALLOCATE(pawtab(itypat)%tphitphj)
   end if
   ABI_ALLOCATE(pawtab(itypat)%phiphj,(mesh_size,ij_size))
   ABI_ALLOCATE(pawtab(itypat)%tphitphj,(mesh_size,ij_size))
   do jln=1,basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       pawtab(itypat)%phiphj(1:mesh_size,kln)=pawtab(itypat)%phi(1:mesh_size,iln)&
&                                            *pawtab(itypat)%phi(1:mesh_size,jln)
       pawtab(itypat)%tphitphj(1:mesh_size,kln)=pawtab(itypat)%tphi(1:mesh_size,iln)&
&                                              *pawtab(itypat)%tphi(1:mesh_size,jln)
     end do
   end do

   if (usekden==1)  then
     pw_mesh_size=pawtab(itypat)%partialwave_mesh_size
     if (allocated(pawtab(itypat)%nablaphi)) then
       ABI_DEALLOCATE(pawtab(itypat)%nablaphi)
     end if
     ABI_ALLOCATE(pawtab(itypat)%nablaphi,(pw_mesh_size,basis_size))
     if (allocated(pawtab(itypat)%tnablaphi)) then
       ABI_DEALLOCATE(pawtab(itypat)%tnablaphi)
     end if
     ABI_ALLOCATE(pawtab(itypat)%tnablaphi,(pw_mesh_size,basis_size))
     ABI_ALLOCATE(der,(pw_mesh_size))
     do iln=1,basis_size
       call nderiv_gen(der,pawtab(itypat)%phi(1:pw_mesh_size,iln),pawrad(itypat))
       pawtab(itypat)%nablaphi(2:pw_mesh_size,iln)=der(2:pw_mesh_size) &
&          -pawtab(itypat)%phi(2:pw_mesh_size,iln)/pawrad(itypat)%rad(2:pw_mesh_size)
       call nderiv_gen(der,pawtab(itypat)%tphi(1:pw_mesh_size,iln),pawrad(itypat))
       pawtab(itypat)%tnablaphi(2:pw_mesh_size,iln)=der(2:pw_mesh_size) &
&          -pawtab(itypat)%tphi(2:pw_mesh_size,iln)/pawrad(itypat)%rad(2:pw_mesh_size)
       call pawrad_deducer0(pawtab(itypat)%nablaphi(1:pw_mesh_size,iln),pw_mesh_size,pawrad(itypat))
       call pawrad_deducer0(pawtab(itypat)%tnablaphi(1:pw_mesh_size,iln),pw_mesh_size,pawrad(itypat))
     end do
     ABI_DEALLOCATE(der)
     pawtab(itypat)%has_nablaphi=2
   end if

!  ==================================================
!  6- COMPUTE Qijl TERMS AND Sij MATRIX

!  Store some usefull quantities
   if (allocated(pawtab(itypat)%phiphj))  then
     ABI_DEALLOCATE(pawtab(itypat)%phiphj)
   end if
   if (allocated(pawtab(itypat)%tphitphj))  then
     ABI_DEALLOCATE(pawtab(itypat)%tphitphj)
   end if
   ABI_ALLOCATE(pawtab(itypat)%phiphj,(mesh_size,ij_size))
   ABI_ALLOCATE(pawtab(itypat)%tphitphj,(mesh_size,ij_size))
   do jln=1,basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       pawtab(itypat)%phiphj  (1:mesh_size,kln)=pawtab(itypat)%phi (1:mesh_size,iln)&
&       *pawtab(itypat)%phi (1:mesh_size,jln)
       pawtab(itypat)%tphitphj(1:mesh_size,kln)=pawtab(itypat)%tphi(1:mesh_size,iln)&
&       *pawtab(itypat)%tphi(1:mesh_size,jln)
     end do
   end do

!  Compute q_ijL and S_ij=q_ij0
   if (allocated(pawtab(itypat)%qijl))  then
     ABI_DEALLOCATE(pawtab(itypat)%qijl)
   end if
   if (allocated(pawtab(itypat)%sij))  then
     ABI_DEALLOCATE(pawtab(itypat)%sij)
   end if
   ABI_ALLOCATE(pawtab(itypat)%qijl,(l_size*l_size,lmn2_size))
   ABI_ALLOCATE(pawtab(itypat)%sij,(lmn2_size))
   pawtab(itypat)%qijl=zero
   pawtab(itypat)%sij=zero
   do klmn=1,lmn2_size
     klm=indklmn_(1,klmn);kln=indklmn_(2,klmn)
     lmin=indklmn_(3,klmn);lmax=indklmn_(4,klmn)
     do ll=lmin,lmax,2
       lm0=ll*ll+ll+1;ff(1)=zero
       ff(2:mesh_size)=(pawtab(itypat)%phiphj  (2:mesh_size,kln)&
&       -pawtab(itypat)%tphitphj(2:mesh_size,kln))&
&       *rad(2:mesh_size)**ll
       call simp_gen(intg,ff,pawrad(itypat))
       do mm=-ll,ll
         isel=pawang%gntselect(lm0+mm,klm)
         if (isel>0) pawtab(itypat)%qijl(lm0+mm,klmn)=intg*pawang%realgnt(isel)
       end do
     end do
     if (klm_diag(klmn)==1) pawtab(itypat)%sij(klmn)= &
&     pawtab(itypat)%qijl(1,klmn)*sqrt(four_pi)
   end do

!  ==================================================
!  6- COMPUTE Eijkl TERMS (Hartree)
!     Compute eventually short-range screened version of Eijkl (Fock)

   if (allocated(pawtab(itypat)%eijkl))  then
     ABI_DEALLOCATE(pawtab(itypat)%eijkl)
   end if
   ABI_ALLOCATE(pawtab(itypat)%eijkl,(lmn2_size,lmn2_size))
   if (abs(hyb_range_fock)>tol8) then
     if (allocated(pawtab(itypat)%eijkl_sr))  then
       ABI_DEALLOCATE(pawtab(itypat)%eijkl_sr)
     end if
     ABI_ALLOCATE(pawtab(itypat)%eijkl_sr,(lmn2_size,lmn2_size))
   end if

!  First loop is for eijkl (Hartree)
!  2nd loop is for eijkl_sr (short-range screened Fock exchange)
   do iloop=1,2
     if (iloop==2.and.abs(hyb_range_fock)<=tol8) cycle
     if (iloop==1) eijkl => pawtab(itypat)%eijkl
     if (iloop==2) eijkl => pawtab(itypat)%eijkl_sr

!    Compute:
!    vhatL(r) according to eq. (A14) in Holzwarth et al., PRB 55, 2005 (1997) [[cite:Holzwarth1997]]
!    intvhatL=$\int_{0}^{r_c}{vhatL(r) shapefunc_L(r) r^2\,dr}$
!    vhatijL =$\int_{0}^{r_c}{vhatL(r) \tilde{\phi}_i \tilde{\phi}_j \,dr}$
!    -----------------------------------------------------------------
     ABI_ALLOCATE(vhatl,(mesh_size))
     ABI_ALLOCATE(vhatijl,(lmn2_size,l_size))
     ABI_ALLOCATE(intvhatl,(l_size))
     intvhatl(:)=zero;vhatl(:)=zero;vhatijl(:,:)=zero
     do il=1,l_size
       vhatl(1)=zero;ff(1)=zero
       ff(2:mesh_size)=pawtab(itypat)%shapefunc(2:mesh_size,il)*rad(2:mesh_size)**2
       if (iloop==1) call poisson(ff,il-1,pawrad(itypat),vhatl)
       if (iloop==2) call poisson(ff,il-1,pawrad(itypat),vhatl,screened_sr_separation=hyb_range_fock)
       vhatl(2:mesh_size)=two*vhatl(2:mesh_size)/rad(2:mesh_size)
       gg(1:mesh_size)=vhatl(1:mesh_size)*ff(1:mesh_size)
       call simp_gen(intvhatl(il),gg,pawrad(itypat))
       do klmn=1,lmn2_size
         kln=indklmn_(2,klmn)
         hh(1:mesh_size)=vhatl(1:mesh_size)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
         call simp_gen(vhatijl(klmn,il),hh,pawrad(itypat))
       end do
     end do
     ABI_DEALLOCATE(vhatl)

!    Compute:
!    eijkl=$ vh1_ijkl - Vhatijkl - Bijkl - Cijkl$
!    With:
!          $vh1_ijkl =\sum_{L,m} {vh1*Gaunt(i,j,Lm)*Gaunt(k,l,Lm)}$
!          $Vhat_ijkl=\sum_{L,m} {vhatijL*Gaunt(i,j,Lm)*q_klL}$
!          $B_ijkl   =\sum_{L,m} {vhatijL*Gaunt(k,l,Lm)*q_ijL}$
!          $C_ijkl   =\sum_{L,m} {intvhatL*q_ijL*q_klL}$
!    and:
!      vh1 according to eq. (A17) in Holzwarth et al., PRB 55, 2005 (1997) [[cite:Holzwarth1997]]
!    Warning: compute only eijkl for (i,j)<=(k,l)
!    -----------------------------------------------------------------
     eijkl(:,:)=zero
     meshsz=pawrad(itypat)%int_meshsz;if (mesh_size>meshsz) ff(meshsz+1:mesh_size)=zero
     do klmn=1,lmn2_size
       klm=indklmn_(1,klmn);kln=indklmn_(2,klmn)
       lmin=indklmn_(3,klmn);lmax=indklmn_(4,klmn)
       do ll=lmin,lmax,2
         lm0=ll*ll+ll+1
         ff(1:meshsz)=pawtab(itypat)%phiphj  (1:meshsz,kln)
         if (iloop==1) call poisson(ff,ll,pawrad(itypat),gg)
         if (iloop==2) call poisson(ff,ll,pawrad(itypat),gg,screened_sr_separation=hyb_range_fock)
         ff(1:meshsz)=pawtab(itypat)%tphitphj(1:meshsz,kln)
         if (iloop==1) call poisson(ff,ll,pawrad(itypat),hh)
         if (iloop==2) call poisson(ff,ll,pawrad(itypat),hh,screened_sr_separation=hyb_range_fock)
         do klmn1=klmn,lmn2_size
           klm1=indklmn_(1,klmn1);kln1=indklmn_(2,klmn1)
           lmin1=indklmn_(3,klmn1);lmax1=indklmn_(4,klmn1)
           vh1=zero
           if ((ll.ge.lmin1).and.(ll.le.lmax1)) then
             ff(1)=zero
             ff(2:meshsz)=(pawtab(itypat)%phiphj  (2:meshsz,kln1)*gg(2:meshsz)&
&             -pawtab(itypat)%tphitphj(2:meshsz,kln1)*hh(2:meshsz))&
&             *two/rad(2:meshsz)
             call simp_gen(vh1,ff,pawrad(itypat))
           end if
           do mm=-ll,ll
             isel =pawang%gntselect(lm0+mm,klm)
             isel1=pawang%gntselect(lm0+mm,klm1)
             if (isel>0.and.isel1>0) then
               rg =pawang%realgnt(isel)
               rg1=pawang%realgnt(isel1)
               ql =pawtab(itypat)%qijl(lm0+mm,klmn)
               ql1=pawtab(itypat)%qijl(lm0+mm,klmn1)
               eijkl(klmn,klmn1)=eijkl(klmn,klmn1)&
&               +(   vh1                *rg *rg1&      ! vh1_ijkl
&              -    vhatijl(klmn ,ll+1)*rg *ql1&     ! Vhat_ijkl
&              -    vhatijl(klmn1,ll+1)*rg1*ql &     ! B_ijkl
&              -    intvhatl(ll+1)     *ql *ql1&     ! C_ijkl
&              )*two_pi
             end if
           end do
         end do
       end do
     end do
     ABI_DEALLOCATE(vhatijl)
     ABI_DEALLOCATE(intvhatl)
   end do ! iloop

!  ==================================================
!  7- COMPUTE gamma_ij TERMS
!  Corrections to get the background right

   if (pawtab(itypat)%usepotzero==1) then
     if (allocated(pawtab(itypat)%gammaij))  then
       ABI_DEALLOCATE(pawtab(itypat)%gammaij)
     end if
     ABI_ALLOCATE(pawtab(itypat)%gammaij,(lmn2_size))
     ABI_ALLOCATE(work,(mesh_size))
     do klmn=1,lmn2_size
       if (klm_diag(klmn)==1) then
         kln=indklmn_(2,klmn)
         ff(1)=zero
         ff(2:mesh_size)=pawtab(itypat)%phiphj(2:mesh_size,kln)-pawtab(itypat)%tphitphj(2:mesh_size,kln)
         ! First, compute q_ij^00
         call simp_gen(intg,ff,pawrad(itypat))
         ! Second, compute phi^2 - tphi^2 - 4pi*r^2 q_ij^00 g_0(r)
         ff(2:mesh_size)= ff(2:mesh_size) - intg*pawtab(itypat)%shapefunc(2:mesh_size,1)*rad(2:mesh_size)**2
         call poisson(ff,0,pawrad(itypat),work)
         ! work is r*vh; should be then multiplied by r to prepare the integration in the sphere
         work(1)=zero ; work(2:mesh_size)=work(2:mesh_size)*rad(2:mesh_size)
         ! Third, average it over the sphere
         call simp_gen(intg,work,pawrad(itypat))
         ! Finally, store it in pawtab%gammaij
         pawtab(itypat)%gammaij(klmn)=intg*four_pi
       else
         pawtab(itypat)%gammaij(klmn)=zero
       end if
     end do
     ABI_DEALLOCATE(work)
   end if

!  ==================================================
!  8- TAKE into account a modified effective mass for the electrons

   if (abs(effmass_free-one)>tol8) then
     if (pawtab(itypat)%has_kij/=2) then
       message='we need kij and has_kij/=2!'
       MSG_BUG(message)
     end if
     if (allocated(pawtab(itypat)%dij0)) then
       pawtab(itypat)%dij0(1:lmn2_size)=pawtab(itypat)%dij0(1:lmn2_size)-pawtab(itypat)%kij(1:lmn2_size)
     end if
     pawtab(itypat)%kij(1:lmn2_size)=pawtab(itypat)%kij(1:lmn2_size)/effmass_free
     if (allocated(pawtab(itypat)%dij0)) then
       pawtab(itypat)%dij0(1:lmn2_size)=pawtab(itypat)%dij0(1:lmn2_size)+pawtab(itypat)%kij(1:lmn2_size)
     end if
   end if

!  ***********************
!  End Loop on atom types
!  ***********************
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(gg)
   ABI_DEALLOCATE(hh)
   ABI_DEALLOCATE(indklmn_)
   ABI_DEALLOCATE(klm_diag)
   ABI_DEALLOCATE(rad)
 end do

 call timab(553,2,tsec)

 DBG_EXIT("COLL")

end subroutine pawinit
!!***

!----------------------------------------------------------------------

!!****f* m_paw_init/paw_gencond
!! NAME
!!   paw_gencond
!!
!! FUNCTION
!!   This routine tests whether we have to call pawinit in one of the optdriver
!!   routines since important values have changed wrt to the previous dataset.
!!   The function uses an internal array to store of the previous values
!!
!!   Usage example:
!!
!!   call paw_gencond(Dtset,gnt_option,"test",call_pawinit)
!!
!!   if (psp_gencond == 1 .or. call_pawinit) then
!!       call pawinit(...)
!!       call paw_gencond(Dtset, gnt_option, "save", call_pawinit)
!!   end if
!!
!!  where psp_gencond is the value returned by pspini.
!!
!! INPUT
!!   Dtset<type(dataset_type)>=all input variables for this dataset
!!   gnt_option=flag activated if pawang%gntselect and pawang%realgnt have to be allocated
!!              also determine the size of these pointers
!!   mode= "test" to test if pawinit must be called
!!         "save" to update the internal variables.
!!         "reset" to reset the internal variables
!!
!! OUTPUT
!!  call_pawinit=True if pawinit must be called. Meaninfull only if mode=="test"
!!
!! SIDE EFFECTS
!!  mode=="save" updates the internal variables.
!!        "reset" reset the internal variables to -1
!!
!! PARENTS
!!      bethe_salpeter,gstate,respfn,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine paw_gencond(Dtset,gnt_option,mode,call_pawinit)

!Arguments ------------------------------------
 integer,intent(in) :: gnt_option
 logical,intent(out) :: call_pawinit
 character(len=*),intent(in) :: mode
 type(dataset_type),intent(in) :: Dtset

!Local variables-------------------------------
!scalars
 integer,save :: gencond(10)=(/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/)

! *********************************************************************

 call_pawinit = .False.
 select case (mode)
 case ("test")

   if (gencond(1)/=Dtset%pawlcutd  .or.gencond(2) /=Dtset%pawlmix  .or.&
&      gencond(3)/=Dtset%pawnphi   .or.gencond(4) /=Dtset%pawntheta.or.&
&      gencond(5)/=Dtset%pawspnorb .or.gencond(6) /=Dtset%pawxcdev .or.&
&      gencond(7)/=Dtset%nsym      .or.gencond(8) /=gnt_option     .or.&
&      gencond(9)/=Dtset%usepotzero.or.gencond(10)/=Dtset%usekden) call_pawinit = .True.

 case ("save")
    ! Update internal values
   gencond(1)=Dtset%pawlcutd  ; gencond(2) =Dtset%pawlmix
   gencond(3)=Dtset%pawnphi   ; gencond(4) =Dtset%pawntheta
   gencond(5)=Dtset%pawspnorb ; gencond(6) =Dtset%pawxcdev
   gencond(7)=Dtset%nsym      ; gencond(8) =gnt_option
   gencond(9)=Dtset%usepotzero; gencond(10)=Dtset%usekden

 case ("reset")
   gencond = -1

 case default
   MSG_BUG("Wrong value for mode: "//trim(mode))
 end select

end subroutine paw_gencond
!!***

!----------------------------------------------------------------------

END MODULE m_paw_init
!!***
