!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdensities
!! NAME
!! pawdensities
!!
!! FUNCTION
!! Compute PAW on-site densities (all-electron ,pseudo and compensation) for a given atom
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, on-site densities are REAL, if 2, COMPLEX (response function only)
!!  iatom=index of current atom (note: this is the absolute index, not the index on current proc)
!!  lm_size=number of (l,m) moments
!!  lmselectin(lm_size)=flags selecting the non-zero LM-moments of on-site densities
!!                      (value of these flags at input; must be .TRUE. for nzlmopt/=1)
!!  nspden=number of spin-density components
!!  nzlmopt=if -1, compute all LM-moments of densities (lmselectin=.true. forced)
!!                 initialize "lmselectout" (index of non-zero LM-moments of densities)
!!          if  0, compute all LM-moments of densities (lmselectin=.true. forced)
!!                 force "lmselectout" to .true. (index of non-zero LM-moments of densities)
!!          if  1, compute only non-zero LM-moments of densities (stored before in "lmselectin")
!!  one_over_rad2(mesh_size)= contains 1/r**2 for each point of the radial grid -optional argument-
!!  opt_compch=flag controlling the accumulation of compensation charge density moments
!!             inside PAW spheres (compch_sph)
!!  opt_dens=flag controlling which on-site density(ies) is (are) computed
!!           0: all on-site densities (all-electron, pseudo and compensation)
!!           1: all-electron and pseudo densities (no compensation)
!!           2: only all-electron density
!!  opt_l=controls which l-moment(s) contribute to the density:
!!        <0 : all l contribute
!!        >=0: only l=opt_l contributes
!!        Note: opt_l>=0 is only compatible with opt_dens=2
!!  opt_print=1 if the densities moments have to be printed out (only if pawprtvol>=2)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data (for the current atom type)
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data (for the current atom)
!!  pawtab <type(pawtab_type)>=paw tabulated starting data (for the current atom type)
!!
!! OUTPUT
!!  nhat1(cplex*mesh_size,lm_size,nspden)= compensation charge on-site density for current atom
!!  rho1(cplex*mesh_size,lm_size,nspden)= all electron on-site density for current atom
!!  trho1(cplex*mesh_size,lm_size,nspden)= pseudo on-site density for current atom
!!  ==== if nzlmopt/=1
!!    lmselectout(lm_size)=flags selecting the non-zero LM-moments of on-site densities
!!                         (value of these flags at output if updated, i.e. if nzlmopt<1)
!!
!!  SIDE EFFECTS
!!  ==== if opt_compch==1
!!    compch_sph=compensation charge integral inside spheres computed over spherical meshes
!!               updated with the contribution of current atom
!!
!! PARENTS
!!      make_efg_onsite,pawdenpot,pawdfptenergy,poslifetime,posratecore
!!
!! CHILDREN
!!      pawrad_deducer0,simp_gen,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawdensities(compch_sph,cplex,iatom,lmselectin,lmselectout,lm_size,nhat1,nspden,nzlmopt,&
&          opt_compch,opt_dens,opt_l,opt_print,pawang,pawprtvol,pawrad,pawrhoij,pawtab,rho1,trho1,&
&          one_over_rad2) ! optional

 use m_profiling_abi
 use defs_basis
 use m_errors
 use m_pawang, only : pawang_type
 use m_pawrad, only : pawrad_type, pawrad_deducer0, simp_gen
 use m_pawtab, only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawdensities'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,iatom,lm_size,nspden,nzlmopt,opt_compch,opt_dens,opt_l,opt_print,pawprtvol
! jmb  real(dp),intent(out) :: compch_sph
 real(dp),intent(inout) :: compch_sph
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab
!arrays
 logical,intent(in) :: lmselectin(lm_size)
 logical,intent(inout) :: lmselectout(lm_size)
 real(dp),intent(in),target,optional :: one_over_rad2(pawtab%mesh_size)
 real(dp),intent(out) :: nhat1(cplex*pawtab%mesh_size,lm_size,nspden*(1-((opt_dens+1)/2)))
 real(dp),intent(out) ::  rho1(cplex*pawtab%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: trho1(cplex*pawtab%mesh_size,lm_size,nspden*(1-(opt_dens/2)))

!Local variables ---------------------------------------
!scalars
 integer :: dplex,ii,ilm,iplex,ir,irhoij,isel,ispden,jrhoij
 integer :: klm,klmn,kln,ll,lmax,lmin,mesh_size
 real(dp) :: m1,mt1,rdum
 character(len=500) :: msg
!arrays
 real(dp) :: compchspha(cplex),compchsphb(cplex),ro(cplex),ro_ql(cplex),ro_rg(cplex)
 real(dp),allocatable :: aa(:),bb(:)
 real(dp),pointer :: one_over_rad2_(:)

! *************************************************************************

 DBG_ENTER("COLL")

!Compatibility tests
 if (opt_dens/=2.and.opt_l>=0) then
   msg='  opt_dens/=2 incompatible with opt_l>=0 !'
   MSG_BUG(msg)
 end if
 if(nzlmopt/=0.and.nzlmopt/=1.and.nzlmopt/=-1) then
   msg='  invalid value for variable "nzlmopt".'
   MSG_BUG(msg)
 end if
 if(nspden>pawrhoij%nspden) then
   msg='  nspden must be <= pawrhoij%nspden !'
   MSG_BUG(msg)
 end if
 if (cplex>pawrhoij%cplex) then
   msg='  cplex must be <= pawrhoij%cplex !'
   MSG_BUG(msg)
 end if
 if (nzlmopt/=1) then
   if (any(.not.lmselectin(1:lm_size))) then
     msg='  With nzlmopt/=1, lmselectin must be true !'
     MSG_BUG(msg)
   end if
 end if
 if (pawang%gnt_option==0) then
   msg='  pawang%gnt_option=0 !'
   MSG_BUG(msg)
 end if

!Various inits
 rho1=zero
 if (opt_dens <2) trho1=zero
 if (opt_dens==0) nhat1=zero
 mesh_size=pawtab%mesh_size;dplex=cplex-1
 if (nzlmopt<1) lmselectout(1:lm_size)=.true.
 if (present(one_over_rad2)) then
   one_over_rad2_ => one_over_rad2
 else
   ABI_ALLOCATE(one_over_rad2_,(mesh_size))
   one_over_rad2_(1)=zero
   one_over_rad2_(2:mesh_size)=one/pawrad%rad(2:mesh_size)**2
 end if

!===== Compute "on-site" densities (n1, ntild1, nhat1) =====
!==========================================================

 do ispden=1,nspden

!  -- Loop over ij channels (basis components)
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     klm =pawtab%indklmn(1,klmn)
     kln =pawtab%indklmn(2,klmn)
     lmin=pawtab%indklmn(3,klmn)
     lmax=pawtab%indklmn(4,klmn)


!    Retrieve rhoij
     if (pawrhoij%nspden/=2) then
       ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,ispden)
     else
       if (ispden==1) then
         ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,1)&
&         +pawrhoij%rhoijp(jrhoij:jrhoij+dplex,2)
       else if (ispden==2) then
         ro(1:cplex)=pawrhoij%rhoijp(jrhoij:jrhoij+dplex,1)
       end if
     end if
     ro(1:cplex)=pawtab%dltij(klmn)*ro(1:cplex)

!    First option: all on-site densities are computed (opt_dens==0)
!    --------------------------------------------------------------
     if (opt_dens==0) then
       do ll=lmin,lmax,2
         do ilm=ll**2+1,(ll+1)**2
           if (lmselectin(ilm)) then
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               ro_ql(1:cplex)=ro(1:cplex)*pawtab%qijl(ilm,klmn)
               ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!              == nhat1(r=0)
               nhat1(1:cplex,ilm,ispden)=nhat1(1:cplex,ilm,ispden) &
&               +ro_ql(1:cplex)*pawtab%shapefunc(1,ll+1) 
!              == rho1(r>0), trho1(r>0), nhat1(r>0)
               do ir=2,mesh_size
                 rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%phiphj  (ir,kln)*one_over_rad2_(ir)
                 trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)=trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%tphitphj(ir,kln)*one_over_rad2_(ir)
                 nhat1(cplex*ir-dplex:ir*cplex,ilm,ispden)=nhat1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_ql(1:cplex)*pawtab%shapefunc(ir,ll+1)
               end do
             end if
           end if
         end do  ! End loops over ll,lm
       end do

!      2nd option: AE and pseudo densities are computed (opt_dens==1)
!      --------------------------------------------------------------
     else if (opt_dens==1) then
       do ll=lmin,lmax,2
         do ilm=ll**2+1,(ll+1)**2
           if (lmselectin(ilm)) then
             isel=pawang%gntselect(ilm,klm)
             if (isel>0) then
               ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!              == rho1(r>0), trho1(r>0)
               do ir=2,mesh_size
                 rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%phiphj  (ir,kln)*one_over_rad2_(ir)
                 trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)=trho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                 +ro_rg(1:cplex)*pawtab%tphitphj(ir,kln)*one_over_rad2_(ir)
               end do
             end if
           end if
         end do  ! End loops over ll,lm
       end do

!      3rd option: only all-electron on-site density is computed (opt_dens==2)
!      -----------------------------------------------------------------------
     else if (opt_dens==2) then
       if (opt_l<0.or.(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*opt_l)) then
         do ll=lmin,lmax,2
           do ilm=ll**2+1,(ll+1)**2
             if (lmselectin(ilm)) then
               isel=pawang%gntselect(ilm,klm)
               if (isel>0) then
                 ro_rg(1:cplex)=ro(1:cplex)*pawang%realgnt(isel)
!                == rho1(r>0)
                 do ir=2,mesh_size
                   rho1(cplex*ir-dplex:ir*cplex,ilm,ispden) =rho1(cplex*ir-dplex:ir*cplex,ilm,ispden)&
&                   +ro_rg(1:cplex)*pawtab%phiphj(ir,kln)*one_over_rad2_(ir)
                 end do
               end if
             end if
           end do  ! End loops over ll, lm
         end do
       end if
     end if

!    -- End loop over ij channels
     jrhoij=jrhoij+pawrhoij%cplex
   end do

!  Scale densities with 1/r**2 and compute rho1(r=0) and trho1(r=0)
   if (cplex==2)  then
     ABI_ALLOCATE(aa,(5))
     ABI_ALLOCATE(bb,(5))
   end if
   if (opt_dens==0.or.opt_dens==1) then
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           if (cplex==1) then
             call pawrad_deducer0(rho1 (:,ilm,ispden),mesh_size,pawrad)
             call pawrad_deducer0(trho1(:,ilm,ispden),mesh_size,pawrad)
           else
             do ii=0,1
               do ir=2,5
                 aa(ir)=rho1 (2*ir-ii,ilm,ispden)
                 bb(ir)=trho1(2*ir-ii,ilm,ispden)
               end do
               call pawrad_deducer0(aa,5,pawrad)
               call pawrad_deducer0(bb,5,pawrad)
               rho1 (2-ii,ilm,ispden)=aa(1)
               trho1(2-ii,ilm,ispden)=bb(1)
             end do
           end if
         end if
       end do
     end do
   else
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           if (cplex==1) then
             call pawrad_deducer0(rho1(:,ilm,ispden),mesh_size,pawrad)
           else
             do ii=0,1
               do ir=2,5
                 aa(ir)=rho1 (2*ir-ii,ilm,ispden)
               end do
               call pawrad_deducer0(aa,5,pawrad)
               rho1(2-ii,ilm,ispden)=aa(1)
             end do
           end if
         end if
       end do
     end do
   end if
   if (cplex==2)  then
     ABI_DEALLOCATE(aa)
     ABI_DEALLOCATE(bb)
   end if

!  -- Test moments of densities and store non-zero ones
   if (nzlmopt==-1) then
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         m1=zero;mt1=zero
         if (cplex==1) then
           m1=maxval(abs(rho1 (1:mesh_size,ilm,ispden)))
           if (opt_dens<2) mt1=maxval(abs(trho1(1:mesh_size,ilm,ispden)))
         else
           do ir=1,mesh_size
             rdum=sqrt(rho1(2*ir-1,ilm,ispden)**2+rho1(2*ir,ilm,ispden)**2)
             m1=max(m1,rdum)
           end do
           if (opt_dens<2) then
             do ir=1,mesh_size
               rdum=sqrt(trho1(2*ir-1,ilm,ispden)**2+trho1(2*ir,ilm,ispden)**2)
               mt1=max(mt1,rdum)
             end do
           end if
         end if
         if (ispden==1) then
           if ((ilm>1).and.(m1<tol16).and.(mt1<tol16)) then
             lmselectout(ilm)=.false.
           end if
         else if (.not.(lmselectout(ilm))) then
           lmselectout(ilm)=((m1>=tol16).or.(mt1>=tol16))
         end if
       end do
     end do
   end if

!  -- Compute integral of (n1-tn1) inside spheres
   if (opt_compch==1.and.ispden==1.and.opt_dens<2) then
     ABI_ALLOCATE(aa,(mesh_size))
     aa(1:mesh_size)=(rho1(1:mesh_size,1,1)-trho1(1:mesh_size,1,1)) &
&     *pawrad%rad(1:mesh_size)**2
! jmb 
!    compchspha(1)=zero
     call simp_gen(compchspha(1),aa,pawrad)
     compch_sph=compch_sph+compchspha(1)*sqrt(four_pi)
     ABI_DEALLOCATE(aa)
   end if

!  -- Print out moments of densities (if requested)
   if (abs(pawprtvol)>=2.and.opt_print==1.and.opt_dens<2) then
     ABI_ALLOCATE(aa,(cplex*mesh_size))
     ABI_ALLOCATE(bb,(cplex*mesh_size))
     if (opt_dens==0) then
       write(msg,'(2a,i3,a,i1,3a)') ch10, &
&       ' Atom ',iatom,' (ispden=',ispden,'):',ch10,&
&       '  ******* Moment of (n1-tn1)  ** Moment of (n1-tn1-nhat1)'
     else
       write(msg,'(2a,i3,a,i1,3a)') ch10, &
&       ' Atom ',iatom,' (ispden=',ispden,'):',ch10,&
&       '  ******* Moment of (n1-tn1)'
     end if
     call wrtout(std_out,msg,'PERS')
     do ll=0,pawtab%lcut_size-1
       do ilm=ll**2+1,(ll+1)**2
         if (lmselectin(ilm)) then
           do iplex=1,cplex
             if (opt_dens==0) then
               do ir=1,mesh_size
                 ii=cplex*(ir-1)+iplex
                 ro(1)=pawrad%rad(ir)**(2+ll)
                 aa(ir)=ro(1)*(rho1(ii,ilm,ispden)-trho1(ii,ilm,ispden))
                 bb(ir)=ro(1)*nhat1(ii,ilm,ispden)
               end do
               call simp_gen(compchspha(iplex),aa,pawrad)
               call simp_gen(compchsphb(iplex),bb,pawrad)
             else
               do ir=1,mesh_size
                 ii=cplex*(ir-1)+iplex
                 ro(1)=pawrad%rad(ir)**(2+ll)
                 aa(ir)=ro(1)*(rho1(ii,ilm,ispden)-trho1(ii,ilm,ispden))
               end do
               call simp_gen(compchspha(iplex),aa,pawrad)
             end if
           end do
           if (opt_dens==0) then
             if (cplex==1) then
               write(msg,'(3x,a,2i2,2(a,es14.7))') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1),&
&               ' **    M=',compchspha(1)-compchsphb(1)
             else
               write(msg,'(3x,a,2i2,2(a,2es14.7))') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1:2),&
&               ' **    M=',compchspha(1:2)-compchsphb(1:2)
             end if
           else
             if (cplex==1) then
               write(msg,'(3x,a,2i2,a,es14.7)') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1)
             else
               write(msg,'(3x,a,2i2,a,2es14.7)') &
&               'l,m=',ll,ilm-(ll**2+ll+1),': M=',compchspha(1:2)
             end if
           end if
           call wrtout(std_out,msg,'PERS')
         end if
       end do
     end do
     ABI_DEALLOCATE(aa)
     ABI_DEALLOCATE(bb)
   end if

!  ----- End loop over spin components
 end do

 if (.not.present(one_over_rad2))  then
   ABI_DEALLOCATE(one_over_rad2_)
 end if

 DBG_EXIT("COLL")

end subroutine pawdensities
!!***
