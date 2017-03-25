!{\src2tex{textfont=tt}}
!!****f* ABINIT/aim_follow
!! NAME
!! aim_follow
!!
!! FUNCTION
!! This routine follows the gradient line starting from the point
!! vv. It stop when it arrives to the atom (nearer than rminl(iat))
!! or - if srch=true - also if it arrives under the already known
!! part of Bader surface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! iatinit,iposinit= indexes of initial atom
!! npmax= maximum number of division in each step
!!
!! OUTPUT
!! iat,ipos= index of final atom
!! nstep= returns the number of step needed
!!
!! SIDE EFFECTS
!! srch=  (true/false) check if the line is outside or
!!             inside the atomic surface.
!! vv(3)= initial point in orthogonal coordinates
!!
!! PARENTS
!!      cpdrv,drvaim,rsurf
!!
!! CHILDREN
!!      critic,onestep,timein,vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine aim_follow(aim_dtset,vv,npmax,srch,iatinit,iposinit,iat,ipos,nstep)

 use defs_basis
 use defs_aimprom
 use defs_parameters
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'aim_follow'
 use interfaces_18_timing
 use interfaces_63_bader, except_this_one => aim_follow
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatinit,iposinit,npmax
 integer,intent(out) :: iat,ipos,nstep
 logical,intent(inout) :: srch
 type(aim_dataset_type),intent(in) :: aim_dtset
!arrays
 real(dp),intent(inout) :: vv(3)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,ii,iph,ires,ith,jj,kk,nit,np,nph,nsi,nth
 real(dp) :: deltar,dg,dist,dph,dth,fac2,facf,h0old,hh,hold,rho,rr,rsmed
 real(dp) :: t1,t2,t3,vcth,vph,vth,wall,xy,xyz
 logical :: fin,ldebold,srchold,stemp,stemp2
 character(len=50) :: formpc
 character(len=500) :: msg
!arrays
 real(dp) :: ev(3),grho(3),hrho(3,3),pom(3),vold(3),vt(3),vt1(3)
 real(dp) :: zz(3,3)

!************************************************************************
 formpc='(":CP",2I5,3F12.8,3E12.4,I4,2E12.4)'


 fin=.false.

 srchold=srch
 ldebold=ldeb
 h0old=h0

 nth=aim_dtset%nth
 nph=aim_dtset%nph

 if (slc==0) then
   rminl(:)=aim_dtset%rmin
 end if

 if (deb) then
   ldeb=.true.
 end if

 call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)

!Initial tests

 if (iat/=0) then
   if (rr<rminl(iat)) then
     fin=.true.
     write(std_out,*) 'rr < rmin iat=',iat,' ipos=',ipos
   elseif (rho<aim_rhomin) then
     fin=.true.
     write(std_out,*) 'CHARGE LT rhomin ',rho,' < ',aim_rhomin
     if (rho<zero) then 
       MSG_ERROR('RHO < 0 !!!')
     end if
   end if
 end if

 facf=aim_fac0
 hh=aim_hmax

 call timein(t1,wall)
 nstep=0
 nsi=0

!the principal cycle

 madw : do while(.not.fin)
   hold=hh

   dg=vnorm(grho,0)
   if (ldeb.or.deb) write(std_out,*) 'dg= ',dg

!  the time test

   call timein(t3,wall)
   t2=t3-t1
   if (t2>300.0) then
     write(std_out,*) 'TIME EXCEEDED 5 min IN FOLLOW'
     write(std_out,*) 'h0 =',h0,'  h =',hh,'  h0old =',h0old,'  dg =',dg
     write(std_out,*) 'facf =',facf
     msg =  'TIME EXCEEDED 5 min IN FOLLOW'
     MSG_ERROR(msg)
   end if

   if (dg<aim_dgmin) then
     write(std_out,*) 'gradient < dgmin ',dg,' < ',aim_dgmin
     fin=.true.
     iat=0
     ipos=0
!    testing for the CP
     if (npc>0) then
       call critic(aim_dtset,vv,ev,zz,aim_dmaxcrit,ires,0)
       if (ires==0) then
         do jj=1,npc
           pom(:)=pc(:,jj)-vv(:)+xatm(:,aim_dtset%batom)
           dist=vnorm(pom,0)
           if (dist<aim_tiny) cycle madw
         end do
         write(std_out,*) 'C.P. found !!'
         npc=npc+1
         do jj=1,3
           pc(jj,npc)=vv(jj)
           evpc(jj,npc)=ev(jj)
           do kk=1,3
             zpc(kk,jj,npc)=zz(kk,jj)
           end do
         end do
         i1=ev(1)/abs(ev(1))
         i2=ev(2)/abs(ev(2))
         i3=ev(3)/abs(ev(3))
         icpc(npc)=i1+i2+i3
         if (icpc(npc)==-3) then           ! pseudoatom handling
           npcm3=npcm3+1
           write(std_out,*) 'Pseudo-atom found !!'
         end if

         call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)
         write(22,formpc) 0,0,(pcrb(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),&
&         ev(1)+ev(2)+ev(3),rho
         write(std_out,formpc) 0,0,(pcrb(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),&
&         ev(1)+ev(2)+ev(3),rho
       else
         write(std_out,*) 'C.P. not found !!'
       end if
     end if

     cycle madw
   end if

   hh=h0/dg
   if (ldeb.or.deb) write(std_out,*) 'h= ',hh,' h0= ',h0,' dg= ',dg
   if (hh>aim_hmax) hh=aim_hmax
!  step modifications

   hh=hh*facf
   if (hh>(hold*aim_hmult)) then
     hh=hold*aim_hmult
   end if

   do ii=1,3
     vold(ii)=vv(ii)
   end do

   nit=0
   hold=hh

!  one step following the gradient line
!  
   call onestep(vv,rho,grho,hh,np,npmax,deltar)
   do while (((np>npmax).or.(deltar>aim_stmax)).and.(deltar>aim_dmin))
     nit=nit+1
     if (nit>5) then
       if (deltar>aim_stmax) then
         write(std_out,*) 'nit > 5 and deltar > stmax   nit=',nit
       else
         write(std_out,*) 'nit > 5 and np > npmax   nit=',nit
       end if
     end if
     do ii=1,3
       vv(ii)=vold(ii)
     end do
     hh=hh*0.3
     call onestep(vv,rho,grho,hh,np,npmax,deltar)
   end do


   nstep=nstep+1
   if (ldeb.or.deb) write(std_out,*) 'h= ',hh

   fac2=hh/hold
   if (fac2>=1._dp) then
     facf=facf*1.2
   else
     if (fac2>=aim_facmin) then
       facf=fac2
     else
       facf=aim_facmin
     end if
   end if

   if (deb.or.ldeb) then
     write(std_out,*) ':POS ',vv
     write(std_out,*) ':RBPOS ',vt1
     write(std_out,*) ':GRAD ',grho
   end if

   call vgh_rho(vv,rho,grho,hrho,rr,iat,ipos,slc)
   dg=vnorm(grho,0)
   pom(:)=vv(:)-xatm(:,iatinit)-atp(:,iposinit)

   if (iat /= 0) then
     fin=.true.
     write(std_out,*) 'r < rmin iat=',iat,' ipos=',ipos
     cycle madw
   end if

   if (rho<aim_rhomin) then
     fin=.true.
     write(std_out,*) 'charge < rhomin ',rho,' < ',aim_rhomin
     if (rho<zero) then 
       MSG_ERROR('RHO < 0 !!!')
     end if
     iat=0
     ipos=0
     cycle madw
   end if

   if (npcm3>0) then
     do jj=1,npc
       if (icpc(jj)==(-3)) then
         pom(:)=pc(:,jj)-vv(:)+xatm(:,aim_dtset%batom)
         dist=vnorm(pom,0)
         if (dist<(aim_dtset%rmin**2*0.1)) then
           iat=0
           ipos=0
           fin=.true.
           write(std_out,*) 'We are inside a pseudo-atom'
           cycle madw
         end if
       end if
     end do
   end if

   nsi=nsi+1

!  surface checking

   if (srch.and.(nsi>=nsimax)) then
     nsi=0
     ith=0
     iph=0
     do ii=1,3
       vt(ii)=vv(ii)-xatm(ii,iatinit)
     end do
     xy=vt(1)*vt(1)+vt(2)*vt(2)
     xyz=xy+vt(3)*vt(3)
     xyz=sqrt(xyz)
     if (xy<aim_snull) then
       vcth=1._dp
       if (vt(3)<0._dp) vcth=-vcth
       vph=0._dp
     else
       vcth=vt(3)/xyz
       vph=atan2(vt(2),vt(1))
     end if
     vth=acos(vcth)
     if (vth<th(1)) then
       ith=0
     else
       if (vth>th(nth)) then
         ith=nth
       else
         do ii=2,nth
           if (vth<th(ii)) then
             ith=ii-1
             exit
           end if
         end do
       end if
     end if

     if (vph<ph(1)) then
       iph=0
     else
       if (vph>ph(nph)) then
         iph=nph
       else
         do ii=2,nph
           if (vph<ph(ii)) then
             iph=ii-1
             exit
           end if
         end do
       end if
     end if

     stemp=(iph>0).and.(iph<nph)
     stemp=stemp.and.(ith>0).and.(ith<nth)

     if (stemp) then
       stemp2=rs(ith,iph)>0._dp
       stemp2=stemp2.and.(rs(ith+1,iph)>0._dp)
       stemp2=stemp2.and.(rs(ith+1,iph+1)>0._dp)
       stemp2=stemp2.and.(rs(ith,iph+1)>0._dp)
       if (stemp2) then
         dth=th(ith+1)-th(ith)
         dph=ph(iph+1)-ph(iph)
         rsmed=rs(ith,iph)*(th(ith+1)-vth)/dth*(ph(iph+1)-vph)/dph
         rsmed=rsmed+rs(ith+1,iph)*(vth-th(ith))/dth*(ph(iph+1)-vph)/dph
         rsmed=rsmed+rs(ith+1,iph+1)*(vth-th(ith))/dth*(vph-ph(iph))/dph
         rsmed=rsmed+rs(ith,iph+1)*(th(ith+1)-vth)/dth*(vph-ph(iph))/dph
         if (rsmed>xyz) then
           write(std_out,*) 'We are inside the surface'
           iat=iatinit
           ipos=iposinit
         else
           write(std_out,*) 'We are outside the surface'
           iat=0
           ipos=0
         end if
         fin=.true.
         cycle madw
       end if
     end if
   end if

 end do madw


 srch=srchold
 ldeb=ldebold
 h0=h0old


end subroutine aim_follow
!!***
