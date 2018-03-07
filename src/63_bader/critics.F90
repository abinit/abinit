!{\src2tex{textfont=tt}}
!!****f* ABINIT/critics
!! NAME
!! critics
!!
!! FUNCTION
!! Search for critical points starting between
!!    atom inxat and its neighbors.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG,MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  aim_dtset= the structured entity containing all input variables
!!  dstmax=maximum distance to search for neighbors
!!  stwo, sthree, sfour: logical switches (TRUE/FALSE) indicating
!!                          to search CP starting in the middle point
!!                          of two, three or four atoms. One of these
!!                          atoms is inxat.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  This routines acts primarily on the data contained in the aim_prom module
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      bschg1,critic,sort_dp,vgh_rho
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine  critics(aim_dtset,inxat,stwo,sthree,sfour,dstmax)

 use defs_basis
 use defs_aimprom
 use defs_parameters
 use defs_abitypes
 use m_sort
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'critics'
 use interfaces_63_bader, except_this_one => critics
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: inxat
 real(dp),intent(in) :: dstmax
 logical,intent(in) :: sfour,sthree,stwo
!no_abirules
 type(aim_dataset_type), intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,iat,ii,ipos,ires,jii,jj,kjj,kk,ll,n1,n2,n3,nb,nc
 integer :: nshell
 real(dp) :: chg,dif1,dif2,diff,dist,olddist,rr
! real(dp) :: ss,uu
 logical :: found,inter
!arrays
 integer :: ibat(nnpos*natom),inat(nnpos*natom),ipibat(nnpos*natom)
 integer :: nnat(nnpos*natom),nr(nnpos*natom)
 real(dp) :: dif(3),dists(nnpos*natom),ev(3),grho(3),hrho(3,3)
 real(dp) :: pom(3),v1(3),v2(3),v3(3),v4(3),vi(3),vt(3),zz(3,3)

!************************************************************************
 vi(:)=xatm(:,inxat)

 nc=0
 do jii=1,nnpos
   do kjj=1,natom
     dist=0._dp
     dif(:)=xatm(:,inxat)-xatm(:,kjj)-atp(:,jii)

!    do ii=1,3
!    dif(ii)=xatm(ii,inxat)-xatm(ii,kjj)-atp(ii,jii)
!    end do
     dist=vnorm(dif,0)
     if (.not.((dist>dstmax).or.(dist<0.001))) then
       nc=nc+1
       dists(nc)=dist
       nnat(nc)=kjj
       inat(nc)=jii
     end if
   end do
 end do
 do n1=1,nc
   nr(n1)=n1
 end do
 call sort_dp(nc,dists,nr,tol14)
 nb=0
 olddist=0._dp
 nshell=0
!write(std_out,*) ':ORIAT ', (xatm(ii,inxat),ii=1,3)
 do n1=1,nc
   n2=nr(n1)
   n3=nnat(n2)
   if (dists(n1)<(2*dists(1))) then
     if ((dists(n1)-olddist)>aim_dlimit) then
       nshell=nshell+1
       olddist=dists(n1)
       if (nshell==5) exit
     end if
     nb=nb+1
     ibat(nb)=n3
     ipibat(nb)=inat(n2)
     write(std_out,*) ':NEIG ',inxat,n3,inat(n2),dists(n1)
!    write(std_out,*) ':POSAT',(xatm(ii,ibat(nb))+atp(ii,ipibat(nb)),ii=1,3)
   else
     exit
   end if
 end do

 npc=0
 npcm3=0

!
!.....SEARCH BETWEEN EACH PAIR OF ATOMS
!

 if (stwo) then
   do jii=1,nb
     do ii=1,3
       v1(ii)=xatm(ii,inxat)
       v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
       vt(ii)=(v1(ii)+v2(ii))/2._dp
     end do
     inter=.true.
     diff=0._dp
     pom(:)=vt(:)
     pom(:)=pom(:)-vi(:)
     diff=vnorm(pom,0)
     if (diff > maxcpdst) inter=.false.
     if (inter) then
       call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
       if (ires==0) then
         found=.false.
         if (npc > 0) then
           do jj=1,npc
             pom(:)=vt(:)-pc(:,jj)
             dist=vnorm(pom,0)
             if (dist < aim_dtset%dpclim) found=.true.
           end do
         end if
         if (.not.found) then
           pom(:)=vt(:)
           call bschg1(pom,-1)
           pcrb(:,npc+1)=pom(:)
           pom(:)=pom(:)-vi(:)
           diff=vnorm(pom,0)
           if (abs(diff) > maxcpdst) found=.true.
         end if
         if (.not.found) then
           npc=npc+1
           do jj=1,3
             pc(jj,npc)=vt(jj)
             evpc(jj,npc)=ev(jj)
             do kk=1,3
               zpc(kk,jj,npc)=zz(kk,jj)
             end do
           end do
           i1=ev(1)/abs(ev(1))
           i2=ev(2)/abs(ev(2))
           i3=ev(3)/abs(ev(3))
           icpc(npc)=i1+i2+i3
           if (icpc(npc)==-3) then
             npcm3=npcm3+1
           end if
           write(std_out,*) 'New critical point found'
           write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
           write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
           write(std_out,'("AUTOVAL: ",3F16.8)') ev
           write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&           ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
           call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
           write(22,'(":PC2",3F10.6,3E12.4,I4,2E12.4)') &
&           (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           write(std_out,'(":PC2",3F10.6,3E12.4,I4,2E12.4)')  &
&           (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           pom(:)=vt(:)-v1(:)
           dif1=vnorm(pom,0)
           pom(:)=vt(:)-v2(:)
           dif2=vnorm(pom,0)
           write(std_out,'(":DISPC ",2F12.8)') dif1,dif2
         end if
       end if
     end if
   end do
 end if
!
!.....SEARCH BETWEEN EACH THREE ATOMS
!
 if(sthree) then
   do jii=1,nb
     do kjj=jii+1,nb
       do ii=1,3
         v1(ii)=xatm(ii,inxat)
         v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
         v3(ii)=xatm(ii,ibat(kjj))+atp(ii,ipibat(kjj))
         vt(ii)=(v1(ii)+v2(ii)+v3(ii))/3._dp
       end do
       inter=.true.
       pom(:)=vt(:)
       pom(:)=pom(:)-vi(:)
       dist=vnorm(pom,0)
       if (abs(diff)>maxcpdst) then
         inter=.false.
         exit
       end if
       if (inter) then
         do jj=1,npc
           pom(:)=pc(:,jj)-vt(:)
           diff=vnorm(pom,0)
           if (diff<aim_dpc0) then
             inter=.false.
             exit
           end if
         end do
       end if
       if (inter) then
         call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
         if (ires==0) then
           found=.false.
           if (npc>0) then
             do jj=1,npc
               pom(:)=vt(:)-pc(:,jj)
               dist=vnorm(pom,0)
               if (dist<aim_dtset%dpclim) then
                 found=.true.
                 exit
               end if
             end do
           end if
           if (.not.found) then
             pom(:)=vt(:)
             call bschg1(pom,-1)
             pcrb(:,npc+1)=pom(:)
             pom(:)=pom(:)-vi(:)
             diff=vnorm(pom,0)
             if (abs(diff)>maxcpdst) found=.true.
           end if
           if (.not.found) then
             npc=npc+1
             do jj=1,3
               pc(jj,npc)=vt(jj)
               evpc(jj,npc)=ev(jj)
               do kk=1,3
                 zpc(kk,jj,npc)=zz(kk,jj)
               end do
             end do
             i1=ev(1)/abs(ev(1))
             i2=ev(2)/abs(ev(2))
             i3=ev(3)/abs(ev(3))
             icpc(npc)=i1+i2+i3
             if (icpc(npc)==-3) then
               npcm3=npcm3+1
             end if
             write(std_out,*) 'New critical point found'
             write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
             write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
             write(std_out,'("AUTOVAL: ",3F16.8)') ev
             write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&             ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
             call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
             write(22,'(":PC3",3F10.6,3E12.4,I4,2E12.4)') &
&             (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
             write(std_out,'(":PC3",3F10.6,3E12.4,I4,2E12.4)') &
&             (pc(jj,npc),jj=1,3),(ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
           end if
         end if
       end if
     end do
   end do
 end if

!
!.....SEARCH BETWEEN EACH FOUR ATOMS
!
 if (sfour) then
   do jii=1,nb
     do kjj=jii+1,nb
       do ll=jii+1,nb
         do ii=1,3
           v1(ii)=xatm(ii,inxat)
           v2(ii)=xatm(ii,ibat(jii))+atp(ii,ipibat(jii))
           v3(ii)=xatm(ii,ibat(kjj))+atp(ii,ipibat(kjj))
           v4(ii)=xatm(ii,ibat(ll))+atp(ii,ipibat(ll))
           vt(ii)=(v1(ii)+v2(ii)+v3(ii)+v4(ii))/4._dp
         end do
         inter=.true.
         pom(:)=vt(:)
         pom(:)=pom(:)-vi(:)
         diff=vnorm(pom,0)
         if (abs(diff)>maxcpdst) then
           inter=.false.
           exit
         end if
         if (inter) then
           do jj=1,npc
             pom(:)=pc(:,jj)-vt(:)
             diff=vnorm(pom,0)
             if (diff < aim_dpc0) then
               inter=.false.
               exit
             end if
           end do
         end if
         if (inter) then
           call critic(aim_dtset,vt,ev,zz,aim_dmaxcs,ires,0)
           if (ires==0) then
             found=.false.
             if (npc>0) then
               do jj=1,npc
                 pom(:)=vt(:)-pc(:,jj)
                 dist=vnorm(pom,0)
                 if (dist < aim_dtset%dpclim) found=.true.
               end do
             end if
             if (.not.found) then
               pom(:)=vt(:)
               pcrb(:,npc+1)=pom(:)
               pom(:)=pom(:)-vi(:)
               diff=vnorm(pom,0)
               if (abs(diff)>maxcpdst) found=.true.
             end if
             if (.not.found) then
               npc=npc+1
               do jj=1,3
                 pc(jj,npc)=vt(jj)
                 evpc(jj,npc)=ev(jj)
                 do kk=1,3
                   zpc(kk,jj,npc)=zz(kk,jj)
                 end do
               end do
               i1=ev(1)/abs(ev(1))
               i2=ev(2)/abs(ev(2))
               i3=ev(3)/abs(ev(3))
               icpc(npc)=i1+i2+i3
               if (icpc(npc)==-3) then
                 npcm3=npcm3+1
               end if
               write(std_out,*) 'New critical point found'
               write(std_out,'("POS: ",3F16.8)') (pc(ii,npc),ii=1,3)
               write(std_out,'("POS in base: ",3F16.8)') (pcrb(ii,npc),ii=1,3)
               write(std_out,'("AUTOVAL: ",3F16.8)') ev
               write(std_out,'("AUTOVEC: ",/,3F16.8,/,3F16.8,/,3F16.8)') &
&               ((zpc(ii,jj,npc),ii=1,3),jj=1,3)
               call vgh_rho(vt,chg,grho,hrho,rr,iat,ipos,0)
               write(22,'(":PC4",3F10.6,3E12.4,I4,2E12.4)') &
&               (pc(jj,npc),jj=1,3), (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
               write(std_out,'(":PC4",3F10.6,3E12.4,I4,2E12.4)') &
&               (pc(jj,npc),jj=1,3), (ev(jj),jj=1,3),icpc(npc),ev(1)+ev(2)+ev(3),chg
             end if
           end if
         end if
       end do
     end do
   end do
 end if

 write(std_out,*) npc
end subroutine critics
!!***
