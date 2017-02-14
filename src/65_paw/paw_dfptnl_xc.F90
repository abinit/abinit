!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_dfptnl_xc
!! NAME
!! paw_dfptnl_xc
!!
!! FUNCTION
!! Compute first-order change of XC potential and contribution to
!! 2nd-order change of XC energy inside a PAW sphere.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! INPUTS
!!  corexc1(cplex_den*nrad)=first-order change of core density on radial grid
!!  cplex_den= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  cplex_vxc= if 1, 1st-order XC potential is complex, if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nrad,pawang%angl_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor1
!!  nhat1(cplex_den*nrad,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1(cplex_den*nrad,lm_size,nspden)=first-order change of density
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!             2 if compensation density (nhat) has to be used in d2Exc and Vxc1
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option=0 or 2 or 3 ==
!!    d2enxc   =returned exchange-cor. contribution to 2nd-order XC energy
!!    d2enxc_im=returned IMAGINARY PART of exchange-cor. contribution to 2nd-order XC energy
!!              (optional argument)
!!
!! SIDE EFFECTS
!!    vxc1(cplex_vxc*nrad,pawang%angl_size,nspden)=1st-order XC potential
!!      Output if option==0 or 1
!!      Unused if option==2
!!      Input  if option==3
!!
!! PARENTS
!!      pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      libxc_functionals_getvxc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw_dfptnl_xc(cplex_1,cplex_2,cplex_3,d3exc1_iat,ixc,kxc,lm_size,lmselect1,lmselect2,lmselect3,&
&                 nhat1,nhat2,nhat3,nkxc,nrad,nspden,pawang,pawrad,rhor1,rhor2,rhor3,usexcnhat)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi, only : xmpi_comm_self,xmpi_sum

 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type,simp_gen
! use m_pawtab,     only : pawtab_type
! use m_paw_an,     only : paw_an_type
! use m_paw_ij,     only : paw_ij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw_dfptnl_xc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_1,cplex_2,cplex_3,ixc,lm_size,nkxc,nrad,nspden,usexcnhat
! real(dp),intent(out),optional :: d2enxc_im
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect1(lm_size),lmselect2(lm_size),lmselect3(lm_size)
! real(dp),intent(in) :: corexc1(cplex_den*nrad)
 real(dp),intent(out) :: d3exc1_iat(2)
 real(dp),intent(in) :: kxc(nrad,pawang%angl_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_1*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat2(cplex_2*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat3(cplex_3*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_1*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor2(cplex_2*nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor3(cplex_3*nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ipts,ir,ispden,jr,lm_size_eff,npts
 logical :: need_impart
 real(dp) :: d3exc1_int
 real(dp) :: rho_dn,rho_up,rhoim_dn,rhoim_up,ro11i,ro11r,ro12i,ro12r,ro21i,ro21r,ro22i,ro22r
 real(dp) :: v11i,v11r,v12i,v12r,v21i,v21r,v22i,v22r
 character(len=500) :: msg
!arrays
! real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),rho1arr(:,:),rho2arr(:,:),rho3arr(:,:),vxc1_(:,:)

! *************************************************************************

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

! if(option<0.or.option>3) then
!   msg='wrong option!'
!   MSG_BUG(msg)
! end if
!!if(xclevel==2) then
! if (.true.) then
!   if (present(d2enxc_im)) then
!     call pawxc3_gga(corexc1,cplex_den,cplex_vxc,d2enxc,ixc,kxc,lm_size,lmselect,nhat1,nkxc,nrad,nspden,&
!&     option,pawang,pawrad,rhor1,usecore,usexcnhat,vxc1,xclevel,&
!&     d2enxc_im=d2enxc_im)
!   else
!     call pawxc3_gga(corexc1,cplex_den,cplex_vxc,d2enxc,ixc,kxc,lm_size,lmselect,nhat1,nkxc,nrad,nspden,&
!&     option,pawang,pawrad,rhor1,usecore,usexcnhat,vxc1,xclevel)
!   end if
!   return
!   msg='GGA is not implemented!'
!   MSG_ERROR(msg)
! end if
! if(option/=3.and.nkxc/=2*min(nspden,2)-1) then
!   msg='nkxc must be 1 or 3!'
!   MSG_BUG(msg)
! end if
! if(nspden==4.and.option/=3) then
!   msg='nspden=4 not implemented (for vxc)!'
!   MSG_ERROR(msg)
! end if
! if(pawang%angl_size==0) then
!   msg='pawang%angl_size=0!'
!   MSG_BUG(msg)
! end if
! if(.not.allocated(pawang%ylmr)) then
!   msg='pawang%ylmr must be allocated!'
!   MSG_BUG(msg)
! end if
! if (option/=1) then
!   if (nrad<pawrad%int_meshsz) then
!     msg='When option=0,2, nrad must be greater than pawrad%int_meshsz!'
!     MSG_BUG(msg)
!   end if
! end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

 npts=pawang%angl_size
 lm_size_eff=min(lm_size,pawang%ylm_size)
! need_impart=present(d2enxc_im)
! if (option/=1) then
!   d2enxc=zero
!   if (need_impart) d2enxc_im=zero
! end if
! if (option<=1) vxc1(:,:,:)=zero

!Special case: no XC applied
 if (ixc==0) then
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

 ABI_ALLOCATE(rho1arr,(cplex_1*nrad,nspden))
 ABI_ALLOCATE(rho2arr,(cplex_2*nrad,nspden))
 ABI_ALLOCATE(rho3arr,(cplex_3*nrad,nspden))
! LIBPAW_ALLOCATE(vxc1_,(cplex_vxc*nrad,nspden))

!Restriction : all cplex must be 1
 if (cplex_1/=1.or.cplex_2/=1.or.cplex_3/=1) then
   msg='All cplex must be one (for the moment...)'
   MSG_BUG(msg)
 end if
!Restriction : all cplex must be 1
 if (nkxc>1) then
   msg='nkxc must be one (<=> nspden=1) (for the moment...)'
   MSG_BUG(msg)
 end if

 d3exc1_iat(:) = zero
 ABI_ALLOCATE(ff,(nrad))

!!----------------------------------------------------------------------
!!----- Loop on the angular part and inits
!!----------------------------------------------------------------------

!Do loop on the angular part (theta,phi)
 do ipts=1,npts

!  Copy the input 1st-order density for this (theta,phi) - PERT1
   rho1arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +rhor1(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect1(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&         +(rhor1(:,ilm,ispden)+nhat1(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT2
   rho2arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +rhor2(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect2(ilm)) rho2arr(:,ispden)=rho2arr(:,ispden) &
&         +(rhor2(:,ilm,ispden)+nhat2(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!  Copy the input 1st-order density for this (theta,phi) - PERT3
   rho3arr(:,:)=zero
   if (usexcnhat==0) then
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +rhor3(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
   else
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect3(ilm)) rho3arr(:,ispden)=rho3arr(:,ispden) &
&         +(rhor3(:,ilm,ispden)+nhat3(:,ilm,ispden))*pawang%ylmr(ilm,ipts)
       end do
     end do
   end if

!   if (usecore==1) then
!     rho1arr(:,1)=rho1arr(:,1)+corexc1(:)
!     if (nspden==2) rho1arr(:,2)=rho1arr(:,2)+half*corexc1(:)
!   end if

!!  
!!  ----------------------------------------------------------------------
!!  ----- Accumulate and store 1st-order change of XC potential
!!  ----------------------------------------------------------------------

!   if (option/=3) then
!!    Non-spin-polarized
!     if(nspden==1)then
!       if (cplex_vxc==1) then
!         if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
!           vxc1_(1:nrad,1)=kxc(1:nrad,ipts,1)*rho1arr(1:nrad,1)
!         else                    ! cplex_vxc==1 and cplex_den==2
!           do ir=1,nrad
!             vxc1_(ir,1)=kxc(ir,ipts,1)*rho1arr(2*ir-1,1)
!           end do
!         end if
!       else
!         if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
!           do ir=1,nrad
!             vxc1_(2*ir-1,1)=kxc(ir,ipts,1)*rho1arr(ir,1)
!             vxc1_(2*ir  ,1)=zero
!           end do
!         else                    ! cplex_vxc==2 and cplex_den==2
!           do ir=1,nrad
!             vxc1_(2*ir-1,1)=kxc(ir,ipts,1)*rho1arr(2*ir-1,1)
!             vxc1_(2*ir  ,1)=kxc(ir,ipts,1)*rho1arr(2*ir  ,1)
!           end do
!         end if
!       end if
!!      Spin-polarized
!     else
!       if (cplex_vxc==1) then
!         if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
!           do ir=1,nrad
!             rho_up=rho1arr(ir,2);rho_dn=rho1arr(ir,1)-rho_up
!             vxc1_(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!             vxc1_(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!           end do
!         else                    ! cplex_vxc==1 and cplex_den==2
!           do ir=1,nrad
!             jr=2*ir-1
!             rho_up=rho1arr(jr,2);rho_dn=rho1arr(jr,1)-rho_up
!             vxc1_(ir,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!             vxc1_(ir,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!           end do
!         end if
!       else
!         if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
!           do ir=1,nrad
!             jr=2*ir-1
!             rho_up=rho1arr(ir,2);rho_dn=rho1arr(ir,1)-rho_up
!             vxc1_(jr,1)=kxc(ir,ipts,1)*rho_up+kxc(ir,ipts,2)*rho_dn
!             vxc1_(jr,2)=kxc(ir,ipts,2)*rho_up+kxc(ir,ipts,3)*rho_dn
!           end do
!         else                    ! cplex_vxc==2 and cplex_den==2
!           do ir=1,nrad
!             jr=2*ir
!             rho_up  =rho1arr(jr-1,2);rho_dn  =rho1arr(jr-1,1)-rho_up
!             rhoim_up=rho1arr(jr  ,2);rhoim_dn=rho1arr(jr  ,1)-rhoim_up
!             vxc1_(jr-1,1)=kxc(ir,ipts,1)*rho_up  +kxc(ir,ipts,2)*rho_dn
!             vxc1_(jr  ,1)=kxc(ir,ipts,1)*rhoim_up+kxc(ir,ipts,2)*rhoim_dn
!             vxc1_(jr-1,2)=kxc(ir,ipts,2)*rho_up  +kxc(ir,ipts,3)*rho_dn
!             vxc1_(jr  ,2)=kxc(ir,ipts,2)*rhoim_up+kxc(ir,ipts,3)*rhoim_dn
!           end do
!         end if
!       end if
!     end if

!     if (option<=1) then
!       vxc1(1:cplex_vxc*nrad,ipts,1:nspden)=vxc1_(1:cplex_vxc*nrad,1:nspden)
!     end if

!   else  ! option==3
!     vxc1_(1:cplex_vxc*nrad,1:nspden)=vxc1(1:cplex_vxc*nrad,ipts,1:nspden)
!   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store 3nd-order change of XC energy
!  ----------------------------------------------------------------------

!    ----- Calculate d2Exc=Int[kxc(r).n^(p1)(r).n^(p2)(r).n^(p3)(r).dr]
!     if (need_impart) then
!       ABI_ALLOCATE(gg,(nrad))
!     end if

!!!    COLLINEAR MAGNETISM
   if (cplex_1==1.and.cplex_2==1.and.cplex_3==1) then ! all cplex are 1 :
     ff(:)=kxc(:,ipts,nspden)*rho1arr(:,nspden)*rho2arr(:,nspden)*rho3arr(:,nspden)
!     if (nspden==2) ff(:)=ff(:)+vxc1_(:,2)*(rho1arr(:,1)-rho1arr(:,2))
!     if (need_impart) gg(:)=zero
   end if
!   else if (cplex_vxc==2.and.cplex_den==2) then  ! cplex_vxc==2 and cplex_den==2
!     if (.not.need_impart) then      ! Real part only
!       do ir=1,nrad
!         jr=2*ir;v11r=vxc1_(jr-1,1);v11i=vxc1_(jr,1)
!         ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
!         ff(ir)=v11r*ro11r+v11i*ro11i
!       end do
!       if (nspden==2) then
!         do ir=1,nrad
!           jr=2*ir;v22r=vxc1_(jr-1,2);v22i=vxc1_(jr,2)
!           ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
!           ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
!           ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!         end do
!           end if
!         else
!           do ir=1,nrad                  ! Real and imaginary parts
!             jr=2*ir;v11r=vxc1_(jr-1,1);v11i=vxc1_(jr,1)
!             ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
!             ff(ir)=v11r*ro11r+v11i*ro11i
!             gg(ir)=v11r*ro11i-v11i*ro11r
!           end do
!           if (nspden==2) then
!             do ir=1,nrad
!               jr=2*ir;v22r=vxc1_(jr-1,2);v22i=vxc1_(jr,2)
!               ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
!               ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
!               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!               gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
!             end do
!           end if
!         end if
!       else                                          ! other cases for cplex_vxc and cplex_den
!         v11i=zero;ro11i=zero
!         do ir=1,nrad
!           jr=cplex_vxc*(ir-1)+1
!           v11r=vxc1_(jr,1);if (cplex_vxc==2) v11i=vxc1_(jr+1,1)
!           jr=cplex_den*(ir-1)+1
!           ro11r=rho1arr(jr,nspden);if (cplex_den==2) ro11i=rho1arr(jr+1,nspden)
!           ff(ir)=v11r*ro11r+v11i*ro11i
!           if (need_impart) gg(ir)=v11r*ro11i-v11i*ro11r
!         end do
!         if (nspden==2) then
!           v22i=zero;ro22i=zero
!           do ir=1,nrad
!             jr=cplex_vxc*(ir-1)+1
!             v22r=vxc1_(jr,2);if (cplex_vxc==2) v22i=vxc1_(jr+1,2)
!             jr=cplex_den*(ir-1)+1
!             ro22r=rho1arr(jr,1)-rho1arr(jr,2)
!             if (cplex_den==2) ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,2)
!             ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
!             gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
!           end do
!         end if
!       end if ! cplex_vxc and cplex_den

!!      NON-COLLINEAR MAGNETISM
!     else
!       if (cplex_vxc==1.and.cplex_den==1) then   ! cplex_vxc==1 and cplex_den==1
!         ff(:)=half*(vxc1_(:,1)*(rho1arr(:,1)+rho1arr(:,4)) &
!&         +vxc1_(:,2)*(rho1arr(:,1)-rho1arr(:,4))) &
!&         +vxc1_(:,3)*rho1arr(:,2) &
!&         -vxc1_(:,4)*rho1arr(:,3)
!         if (need_impart) gg(:)=zero
!       else                                      ! other cases for cplex_vxc and cplex_den

!!        V is stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!!        N is stored as : n, m_x, m_y, mZ          (each are complex)
!         do ir=1,nrad
!           jr=cplex_vxc*(ir-1)+1
!           v11r= vxc1_(jr,1);v22r= vxc1_(jr,2)
!           v12r= vxc1_(jr,3);v21i=-vxc1_(jr,1)
!           if (cplex_vxc==2) then
!             v11i= vxc1_(jr+1,1);v22i= vxc1_(jr+1,2)
!             v12i= vxc1_(jr+1,3);v21r= vxc1_(jr+1,1)
!           else
!             v11i=zero;v22i=zero
!             v12i=zero;v21i=zero
!           end if
!           jr=cplex_den*(ir-1)+1
!           ro11r= rho1arr(jr,1)+rho1arr(jr,4)
!           ro22r= rho1arr(jr,1)-rho1arr(jr,4)
!           ro12r= rho1arr(jr,2);ro12i=-rho1arr(jr,3)
!           ro21r= rho1arr(jr,2);ro21i= rho1arr(jr,3)
!           if (cplex_den==2) then
!             ro11i=rho1arr(jr+1,1)+rho1arr(jr+1,4)
!             ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,4)
!             ro12r=ro12r+rho1arr(jr+1,3);ro12i=ro12i+rho1arr(jr+1,2)
!             ro21r=ro21r-rho1arr(jr+1,3);ro21i=ro21i+rho1arr(jr+1,2)
!           else
!             ro11i=zero;ro22i=zero
!           end if
!!          Real part
!           ff(ir)=half*(v11r*ro11r+v11i*ro11i+v22r*ro22r+v22i*ro22i &
!&           +v12r*ro12r+v12i*ro12i+v21r*ro21r+v21i*ro21i)
!!          Imaginary part
!           if (need_impart) gg(ir)=half*(v11r*ro11i-v11i*ro11r+v22r*ro22i-v22i*ro22r &
!&           +v12r*ro12i-v12i*ro12r+v21r*ro21i-v21i*ro21r)
!         end do
!       end if ! cplex_vxc and cplex_den
!     end if ! nspden

     ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
     call simp_gen(d3exc1_int,ff,pawrad)
     d3exc1_iat(1)=d3exc1_iat(1)+d3exc1_int*pawang%angwgth(ipts)

!     if (need_impart) then
!       gg(:)=gg(:)*pawrad%rad(:)**2
!       call simp_gen(vxcrho,gg,pawrad)
!       d2enxc_im=d2enxc_im+vxcrho*pawang%angwgth(ipts)
!       ABI_DEALLOCATE(gg)
!     end if

!   end if

!  ----- End of the loop on npts (angular part)
 end do

 d3exc1_iat = d3exc1_iat*four_pi

!!Add the four*pi factor of the angular integration
! if (option/=1) then
!   d2enxc=d2enxc*four_pi
!   if (need_impart) d2enxc_im=d2enxc_im*four_pi
! end if

 ABI_DEALLOCATE(ff)
 ABI_DEALLOCATE(rho1arr)
 ABI_DEALLOCATE(rho2arr)
 ABI_DEALLOCATE(rho3arr)
! LIBPAW_DEALLOCATE(vxc1_)

end subroutine paw_dfptnl_xc
!!***
