!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrho
!! NAME
!! integrho
!!
!! FUNCTION
!! This routine integrates the electron density inside the
!! atomic surface already calculated - it reads the file *.surf
!! The radial integration is always performed with splines and
!! the two angular integrations with Gauss quadrature
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! aim_dtset = the structured entity containing all input variables
!! znucl_batom=the nuclear charge of the Bader atom
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routine works primarily on the data contained in the aimfields and aimprom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      bschg1,coeffs_gausslegint,spline,vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine integrho(aim_dtset,znucl_batom)

 use defs_basis
 use defs_aimfields
 use defs_aimprom
 use defs_parameters
 use defs_abitypes
 use m_profiling_abi
 use m_splines
 use m_errors
 
 use m_numeric_tools, only : coeffs_gausslegint

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrho'
 use interfaces_63_bader, except_this_one => integrho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: batom,chs,iat,ii,inx,inxf,ipos,jj,kk,ll,nn,nph,nth
 real(dp) :: chg,chgint,cintr,ct1,ct2,lder,nsphe,phimax,phimin,rder
 real(dp) :: rsmax,rsmin,ss,stp,themax,themin,uu
 real(dp) :: znucl_batom,zz
 logical :: gaus,weit
!arrays
 real(dp) :: grho(3),hrho(3,3),shift(3),unvec(3),vv(3)
 real(dp),allocatable :: ncrho(:),nsp2(:),nsp3(:),nsp4(:),rdint(:,:),rr(:)
 real(dp),allocatable :: vdd(:),vrho(:),wgrs(:,:),work(:)

! *********************************************************************

 gaus=.true.
 weit=.true.

 write(std_out,*) 'npt = ',aim_dtset%npt

 rewind(unts)
 read(unts,*) batom,shift  ! Warning : batom is read, instead of coming from aim_dtset
 read(unts,*) nth,themin,themax ! Warning : these numbers are read, instead of coming from aim_dtset
 read(unts,*) nph,phimin,phimax ! Warning : these numbers are read, instead of coming from aim_dtset

 write(std_out,*) 'NTH NPH ',nth,nph

 ABI_ALLOCATE(wgrs,(nth,nph))
 ABI_ALLOCATE(rdint,(nth,nph))

 do ii=1,nth
   do jj=1,nph
     if (weit) then
       read(unts,*) th(ii),ph(jj),rs(ii,jj),wgrs(ii,jj)
     else
       read(unts,*) th(ii),ph(jj),rs(ii,jj)
     end if
   end do
 end do
 read(unts,*) rsmin,rsmax


 if (gaus) then
   ct1=cos(themin)
   ct2=cos(themax)
   call coeffs_gausslegint(ct1,ct2,cth,wcth,nth)
   call coeffs_gausslegint(phimin,phimax,ph,wph,nph)
 end if

 do ii=1,nth
   do jj=1,nph
     if (.not.weit) then
       if (gaus) then
         wgrs(ii,jj)=wcth(ii)*wph(jj)
       else
         wgrs(ii,jj)=1._dp
       end if
     end if
   end do
 end do


 do ii=1,nth
   do jj=1,nph
     if (rs(ii,jj) < rsmin) rsmin=rs(ii,jj)
   end do
 end do


!INTEGRATION OF THE CORE DENSITY

 nn=typat(batom)
 kk=ndat(nn)


!spherical integration of the core density in the sphere
!of the minimal Bader radius

!COEF. FOR SPHERICAL INTEGRATION

 ABI_ALLOCATE(nsp2,(kk))
 ABI_ALLOCATE(nsp3,(kk))
 ABI_ALLOCATE(nsp4,(kk))
 ABI_ALLOCATE(ncrho,(kk))

 do ii=1,kk
   ncrho(ii)=crho(ii,nn)*4._dp*pi*rrad(ii,nn)*rrad(ii,nn)
   nsp3(ii)=4._dp*pi*(2._dp*crho(ii,nn)+2._dp*rrad(ii,nn)*sp2(ii,nn)+&
&   rrad(ii,nn)*rrad(ii,nn)*sp3(ii,nn))
 end do

 if (rsmin < rrad(ndat(nn),nn)) then        ! search index
   inx=0
   if (rsmin < rrad(1,nn)) then
     MSG_ERROR('absurd')
   elseif (rsmin > rrad(ndat(nn),nn)) then
     inx=ndat(nn)
   else
     do while (rsmin >= rrad(inx+1,nn))
       inx=inx+1
     end do
   end if
 else
   inx=ndat(nn)
 end if

 cintr=4._dp/3._dp*pi*rrad(1,nn)**3*crho(1,nn)

!spline integration

 do ii=1,inx-1
   uu=rrad(ii+1,nn)-rrad(ii,nn)
   cintr=cintr+(ncrho(ii)+ncrho(ii+1))*uu/2._dp-uu*uu*uu/2.4d1*(nsp3(ii)+nsp3(ii+1))
 end do
 if (inx/=ndat(nn)) then
   uu=rsmin-rrad(inx,nn)
   zz=rrad(inx+1,nn)-rsmin
   ss=rrad(inx+1,nn)-rrad(inx,nn)
   cintr=cintr+ncrho(inx)/2._dp*(ss-zz*zz/ss)+ncrho(inx+1)/2._dp*uu*uu/ss+&
   nsp3(inx)/1.2d1*(zz*zz*ss-zz*zz*zz*zz/2._dp/ss-ss*ss*ss/2._dp)+&
   nsp3(inx+1)/1.2d1*(uu*uu*uu*uu/2._dp/ss-uu*uu*ss)
 end if


!INTEGRATION OF THE REST OF THE CORE DENSITY
!(for gauss quadrature)
!For the Gauss quadrature it is added
!to the radial integrated valence density

 rdint(:,:)=0._dp
 nsphe=0._dp
 do ii=1,nth
   do jj=1,nph
     if (inx==ndat(nn)) cycle
     inxf=inx
     if (rs(ii,jj) < rsmin) then
       write(std_out,*) rs(ii,jj),rsmin
       MSG_ERROR('in surface')
     elseif (rs(ii,jj) > rrad(ndat(nn),nn)) then
       inxf=ndat(nn)
     else
       do while (rs(ii,jj) >= rrad(inxf+1,nn))
         inxf=inxf+1
       end do
     end if

     if (inxf==inx) then
       uu=rrad(inx+1,nn)-rs(ii,jj)
       zz=rrad(inx+1,nn)-rsmin
       ss=rrad(inx+1,nn)-rrad(inx,nn)

       rdint(ii,jj)=(ncrho(inx)/2._dp/ss-nsp3(inx)/1.2d1*ss)*(zz*zz-uu*uu)+&
       nsp3(inx)/2.4d1/ss*(zz**4-uu**4)
       uu=rs(ii,jj)-rrad(inx,nn)
       zz=rsmin-rrad(inx,nn)
       rdint(ii,jj)=rdint(ii,jj)+(uu*uu-zz*zz)*(ncrho(inx+1)/2._dp/ss-nsp3(inx+1)/1.2d1*ss)+&
       nsp3(inx+1)/2.4d1/ss*(uu**4-zz**4)
     else
       uu=rrad(inx+1,nn)-rsmin
       zz=rsmin-rrad(inx,nn)

       rdint(ii,jj)=ncrho(inx)/2._dp/ss*uu*uu+ncrho(inx+1)/2._dp*(ss-zz*zz/ss)+&
       nsp3(inx)/1.2d1*(uu**4/2._dp/ss-uu*uu*ss)+nsp3(inx+1)/1.2d1*(zz*zz*ss-ss**3/2._dp-zz**4/2._dp/ss)
       if (inxf > inx+1) then
         do kk=inx+1,inxf-1
           uu=rrad(kk+1,nn)-rrad(kk,nn)
           rdint(ii,jj)=rdint(ii,jj)+(ncrho(kk)+ncrho(kk+1))*uu/2._dp-uu*uu*uu/2.4d1*(nsp3(kk)+nsp3(kk+1))
         end do
       end if

       if (inxf/=ndat(nn)) then
         uu=rs(ii,jj)-rrad(inxf,nn)
         zz=rrad(inxf+1,nn)-rs(ii,jj)
         ss=rrad(inxf+1,nn)-rrad(inxf,nn)
         rdint(ii,jj)=rdint(ii,jj)+ncrho(inxf)/2._dp*(ss-zz*zz/ss)+ncrho(inxf+1)/2._dp*uu*uu/ss+&
         nsp3(inxf)/1.2d1*(zz*zz*ss-zz*zz*zz*zz/2._dp/ss-ss*ss*ss/2._dp)+&
         nsp3(inxf+1)/1.2d1*(uu*uu*uu*uu/2._dp/ss-uu*uu*ss)
       end if
     end if
     rdint(ii,jj)=rdint(ii,jj)/4._dp/pi
     nsphe=nsphe+rdint(ii,jj)*wgrs(ii,jj)
   end do
 end do
 nsphe=nsphe*(pi/(themin-themax))*(two_pi/(phimax-phimin))

 write(untout,*)
 write(untout,*) "CHARGE INTEGRATION"
 write(untout,*) "=================="
 write(untout,'(" Core density contribution: ",/,/,"    ",F16.8)') cintr+nsphe

 write(std_out,*) ':INTECOR ', cintr+nsphe

 ABI_DEALLOCATE(ncrho)
 ABI_DEALLOCATE(nsp2)
 ABI_DEALLOCATE(nsp3)
 ABI_DEALLOCATE(nsp4)

!INTEGRATION OF THE VALENCE DENSITY

 ABI_ALLOCATE(rr,(aim_dtset%npt+1))
 ABI_ALLOCATE(vrho,(aim_dtset%npt+1))
 ABI_ALLOCATE(vdd,(aim_dtset%npt+1))

!in the case of the only irho appelation

 nn=0
 do ii=-3,3
   do jj=-3,3
     do kk=-3,3
       nn=nn+1
       atp(1,nn)=ii*1._dp
       atp(2,nn)=jj*1._dp
       atp(3,nn)=kk*1._dp
       call bschg1(atp(:,nn),1)
       if ((ii==0).and.(jj==0).and.(kk==0)) ipos=nn
     end do
   end do
 end do
 nnpos=nn
 iat=batom

!XG020629 There is a problem with this routine
!(or vgh_rho), when one uses the PGI compiler :
!The following line is needed, otherwise, iat and ipos
!are set to 0 inside vgh_now. Why ????
 write(std_out,*)' integrho : iat,ipos=',iat,ipos
!

 nsphe=0._dp
 ABI_ALLOCATE(work,(aim_dtset%npt+1))
 do ii=1,nth
   do jj=1,nph

     stp=rs(ii,jj)/aim_dtset%npt
     unvec(1)=sin(th(ii))*cos(ph(jj))
     unvec(2)=sin(th(ii))*sin(ph(jj))
     unvec(3)=cos(th(ii))
     do kk=0,aim_dtset%npt
       rr(kk+1)=kk*stp
       vv(:)=xatm(:,batom)+kk*stp*unvec(:)
       chs=-2
       call vgh_rho(vv,chg,grho,hrho,uu,iat,ipos,chs)
       vrho(kk+1)=chg*rr(kk+1)*rr(kk+1)
       if (kk==aim_dtset%npt) then
         rder=0._dp
         do ll=1,3
           rder=rder+grho(ll)*unvec(ll)
         end do
         rder=rder*rr(kk+1)*rr(kk+1)+2._dp*rr(kk+1)*chg
       end if
     end do
     lder=0._dp
     kk=aim_dtset%npt+1
     call spline(rr,vrho,kk,lder,rder,vdd)

!    INTEGRATION

     do kk=1,aim_dtset%npt
       rdint(ii,jj)=rdint(ii,jj)+stp/2._dp*(vrho(kk)+vrho(kk+1))&
&       -stp*stp*stp/24._dp*(vdd(kk)+vdd(kk+1))
     end do
     nsphe=nsphe+rdint(ii,jj)*wgrs(ii,jj)
   end do
 end do
 ABI_DEALLOCATE(work)

 if (gaus.or.weit) then
   nsphe=nsphe*(pi/(themin-themax))*(two_pi/(phimax-phimin))
 else
   nsphe=nsphe/(nth*nph)*2.0*two_pi
 end if
 chgint=cintr+nsphe

 write(untout,'(/," Different density contributions: Core (only spherical part) and the rest ",/,/,"      ",2F16.8)') &
& cintr, nsphe
 write(untout,'(/,a,i4,a,f14.8)') ' For atom number ',batom,', the number of electrons in the Bader volume is ',chgint
 write(untout,'(a,f15.7,a,f17.8)') ' The nuclear charge is',znucl_batom,', so that the Bader charge is ',znucl_batom-chgint
 write(untout,*)
 write(std_out,*) ':INTEPAR ', cintr, nsphe
 write(std_out,*) ':RHOTOT ',batom,chgint

end subroutine integrho
!!***
