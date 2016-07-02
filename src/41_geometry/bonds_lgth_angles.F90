!{\src2tex{textfont=tt}}
!!****f* ABINIT/bonds_lgth_angles
!! NAME
!! bonds_lgth_angles
!!
!! FUNCTION
!! From list of coordinates and primitive translations, output
!! a list of bonds lengths and bond angles.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  coordn = maximum coordination number to be taken into account
!!  fnameabo_app_geo=name of file for _GEO data
!!  natom  = number of atoms in unit cell
!!  ntypat = number of types of atoms in unit cell.
!!  rprimd(3,3)  = real space dimensional primitive translations (bohr)
!!  typat(natom) = type integer for each atom in cell
!!  znucl(ntypat)= real(dp), atomic number of atom type
!!  xred(3,natom)= reduced coordinates of atoms
!!
!! OUTPUT
!! data written in file fnameabo_app_geo
!!
!! NOTES
!!  The tolerance tol8 aims at giving a machine-independent ordering.
!!  (this trick is used in bonds.f, listkk.f, prtrhomxmn.f and rsiaf9.f)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine bonds_lgth_angles(coordn,fnameabo_app_geo,natom,ntypat,rprimd,typat,xred,znucl)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_atomdata
 use m_errors

 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bonds_lgth_angles'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => bonds_lgth_angles
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: coordn,natom,ntypat
 character(len=*),intent(in) :: fnameabo_app_geo
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rprimd(3,3),znucl(ntypat)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: done,ia,ib,ic,ii,ineighb,jneighb,mneighb,mu,ndig,nu,t1,t2,t3,tmax,temp_unit
 real(dp) :: adotb,asq,bsq,co,length,sq,thdeg
!real(dp)u1,u2,u3,v1,v2,v3
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: list_neighb(:,:,:)
 real(dp) :: bab(3),bac(3),dif(3),rmet(3,3)
 real(dp),allocatable :: sqrlength(:),xangst(:,:),xcart(:,:)
 character(len=8),allocatable :: iden(:)

! *************************************************************************

! Statement functions are obsolete
! dot(u1,u2,u3,v1,v2,v3)=rmet(1,1)*u1*v1+rmet(2,1)*u2*v1+&
!& rmet(3,1)*u3*v1+rmet(1,2)*u1*v2+rmet(2,2)*u2*v2+&
!& rmet(3,2)*u3*v2+rmet(1,3)*u1*v3+rmet(2,3)*u2*v3+rmet(3,3)*u3*v3

!Initialize the file
 write(message, '(a,a,a)' )' bonds_lgth_angles : about to open file ',trim(fnameabo_app_geo),ch10
 call wrtout(std_out,message,'COLL'); call wrtout(ab_out,message,'COLL')

 if (open_file(fnameabo_app_geo,message,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
   MSG_ERROR(message)
 end if
 rewind(temp_unit)

 write(message, '(a,a)' ) ch10,' ABINIT package : GEO file '
 call wrtout(temp_unit,message,'COLL')

!Compute maximum number of neighbors is the neighbor list,
!from the indicative coordination number
!Note : the following formula includes next nearest neighbors, but not others
 mneighb=1+coordn+coordn*(coordn-1)

 write(message, '(a,a,i2,a,a,i4,a,a,a,i4,a)' ) ch10,&
& ' Maximal coordination number, as estimated by the user : ',coordn,ch10,&
& '  giving a maximum of ',coordn*coordn,&
& ' nearest neighbors and next nearest neighbors, ',ch10,&
& '                  and ',(coordn*(coordn-1))/2,&
& ' distinct angles between nearest neighbors'
 call wrtout(temp_unit,message,'COLL')

!Compute metric tensor in real space rmet
 do nu=1,3
   do mu=1,3
     rmet(mu,nu)=rprimd(1,mu)*rprimd(1,nu)+&
&     rprimd(2,mu)*rprimd(2,nu)+&
&     rprimd(3,mu)*rprimd(3,nu)
   end do
 end do

 write(message, '(a,a)' )ch10,' Primitive vectors of the periodic cell (bohr)'
 call wrtout(temp_unit,message,'COLL')
 do nu=1,3
   write(message, '(1x,a,i1,a,3f10.5)' ) '  R(',nu,')=',rprimd(:,nu)
   call wrtout(temp_unit,message,'COLL')
 end do

 write(message, '(a,a)' ) ch10,&
& ' Atom list        Reduced coordinates          Cartesian coordinates (bohr)'
 call wrtout(temp_unit,message,'COLL')

!Set up a list of character identifiers for all atoms : iden(ia)
 ABI_ALLOCATE(iden,(natom))
 iden(:)='        '
 do ia=1,natom
   ndig=int(log10(dble(ia)+0.5d0))+1
   call atomdata_from_znucl(atom,znucl(typat(ia)))
   if(ndig==1) write(iden(ia), '(a,a,i1,a)' )  atom%symbol,'(',ia,')   '
   if(ndig==2) write(iden(ia), '(a,a,i2,a)' )  atom%symbol,'(',ia,')  '
   if(ndig==3) write(iden(ia), '(a,a,i3,a)' )  atom%symbol,'(',ia,') '
   if(ndig==4) write(iden(ia), '(a,a,i4,a)' )  atom%symbol,'(',ia,')'
   if(ndig>4)then
     close(temp_unit)
     write(message, '(a,i8,a,a)' )&
&     'bonds_lgth_angles cannot handle more than 9999 atoms, while natom=',natom,ch10,&
&     'Action : decrease natom, or contact ABINIT group.'
     MSG_BUG(message)
   end if
 end do

!Compute cartesian coordinates, and print reduced and cartesian coordinates
!then print coordinates in angstrom, with the format neede for xmol
 ABI_ALLOCATE(xangst,(3,natom))
 ABI_ALLOCATE(xcart,(3,natom))
 call xred2xcart(natom,rprimd,xcart,xred)
 xangst(:,:)=xcart(:,:)*Bohr_Ang

 do ia=1,natom
   write(message, '(a,a,3f10.5,a,3f10.5)' ) &
&   '   ',iden(ia),(xred(ii,ia)+tol10,ii=1,3),&
&   '    ',(xcart(ii,ia)+tol10,ii=1,3)
   call wrtout(temp_unit,message,'COLL')
 end do

 write(message, '(a,a,a,a,i4,a)' )ch10,&
& ' XMOL data : natom, followed by cartesian coordinates in Angstrom',&
& ch10,ch10,natom,ch10
 call wrtout(temp_unit,message,'COLL')

 do ia=1,natom
   call atomdata_from_znucl(atom,znucl(typat(ia)))
   write(message, '(a,a,3f10.5)' )'   ',atom%symbol,xangst(1:3,ia)
   call wrtout(temp_unit,message,'COLL')
 end do

 ABI_DEALLOCATE(xangst)
 ABI_DEALLOCATE(xcart)

 ABI_ALLOCATE(list_neighb,(0:mneighb+1,4,2))
 ABI_ALLOCATE(sqrlength,(0:mneighb+1))

!Compute list of neighbors
 do ia=1,natom

   write(message, '(a,a,a,a,a,a,a,a,a)' ) ch10,'===========',&
&   '=====================================================================',&
&   ch10,' ',iden(ia),ch10,ch10,' Bond lengths '
   call wrtout(temp_unit,message,'COLL')

!  Search other atoms for bonds, but must proceed
!  in such a way to consider a search box sufficiently large,
!  so increase the size of the search box until the
!  final bond length list do not change
   do tmax=0,5

!    Set initial list of neighbors to zero,
!    and initial square of bond lengths to a very large number.
!    Note that the dimension is larger than neighb to ease
!    the later sorting : neighbors 0 and neighb+1 are non-existent, while
!    neighbor 1 will be the atom itself ...
     list_neighb(0:mneighb+1,1:4,1)=0
     sqrlength(1:mneighb+1)=huge(0.0d0)
     sqrlength(0)=-1.0d0

!    Here search on all atoms inside the box defined by tmax
     do ib=1,natom
       do t3=-tmax,tmax
         do t2=-tmax,tmax
           do t1=-tmax,tmax
             dif(1)=xred(1,ia)-(xred(1,ib)+dble(t1))
             dif(2)=xred(2,ia)-(xred(2,ib)+dble(t2))
             dif(3)=xred(3,ia)-(xred(3,ib)+dble(t3))
             sq=rsdot(dif(1),dif(2),dif(3),dif(1),dif(2),dif(3),rmet)

!            Insert the atom at the proper place in the neighbor list.
             do ineighb=mneighb,0,-1
!              Note the tolerance
               if(sq+tol8>sqrlength(ineighb))then
                 sqrlength(ineighb+1)=sq
                 list_neighb(ineighb+1,1,1)=ib
                 list_neighb(ineighb+1,2,1)=t1
                 list_neighb(ineighb+1,3,1)=t2
                 list_neighb(ineighb+1,4,1)=t3
!                DEBUG
!                if(ineighb/=mneighb)then
!                write(std_out,*)' '
!                do ii=1,mneighb
!                write(std_out,*)ii,sqrlength(ii)
!                end do
!                end if
!                ENDDEBUG
                 exit
               else
                 sqrlength(ineighb+1)=sqrlength(ineighb)
                 list_neighb(ineighb+1,1:4,1)=list_neighb(ineighb,1:4,1)
               end if
             end do

           end do
         end do
       end do
!      end ib loop:
     end do

!    Now, check that the box defined by tmax was large enough :
!    require the present and old lists to be the same
     done=0

     if(tmax>0)then
       done=1
       do ineighb=1,mneighb
!        DEBUG
!        write(std_out,'(5i5,f12.5)' )ineighb,list_neighb(ineighb,1:4,1),&
!        &                                    sqrlength(ineighb)
!        write(std_out,'(5i5)' )ineighb,list_neighb(ineighb,1:4,2)
!        ENDDEBUG
         if( list_neighb(ineighb,1,1)/=list_neighb(ineighb,1,2) .or. &
&         list_neighb(ineighb,2,1)/=list_neighb(ineighb,2,2) .or. &
&         list_neighb(ineighb,3,1)/=list_neighb(ineighb,3,2) .or. &
&         list_neighb(ineighb,4,1)/=list_neighb(ineighb,4,2)       )then
           done=0
         end if
       end do
     end if

!    If done==1, then one can exit the loop : the correct list of
!    neighbors is contained in list_neighb(1:neighb,1:4,1),
!    with the first neighbor being the atom itself
     if(done==1)exit

!    If the work is not done, while tmax==5, then there is a problem .
     if(tmax==5)then
       close(temp_unit)
       write(message, '(2a)' )&
&       'Did not succeed to generate a reliable list of bonds ',&
&       'since tmax is exceeded.'
       MSG_BUG(message)
     end if

!    Copy the new list into the old list.
     list_neighb(1:mneighb,1:4,2)=list_neighb(1:mneighb,1:4,1)

!    Loop on tmax (note that there are exit instruction inside the loop)
   end do



!  Output the bond list
   do ineighb=2,mneighb
     ib=list_neighb(ineighb,1,1)
     length=sqrt(sqrlength(ineighb))
     write(message, '(a,a,a,a,3i2,t27,a,f10.5,a,f9.5,a)' )&
&     '  ',trim(iden(ia)),' - ',trim(iden(ib)),&
&     list_neighb(ineighb,2:4,1),'bond length is ',&
&     length,' bohr  ( or ',Bohr_Ang*length,' Angst.)'
     call wrtout(temp_unit,message,'COLL')
   end do

!  Output the angle list
   if(coordn>1)then

     write(message, '(a,a)' ) ch10,' Bond angles '
     call wrtout(temp_unit,message,'COLL')

     do ineighb=2,coordn
       do jneighb=ineighb+1,coordn+1

         ib=list_neighb(ineighb,1,1)
         ic=list_neighb(jneighb,1,1)
         do mu=1,3
           bab(mu)=xred(mu,ib)+dble(list_neighb(ineighb,1+mu,1))-xred(mu,ia)
           bac(mu)=xred(mu,ic)+dble(list_neighb(jneighb,1+mu,1))-xred(mu,ia)
         end do
         asq=rsdot(bab(1),bab(2),bab(3),bab(1),bab(2),bab(3),rmet)
         bsq=rsdot(bac(1),bac(2),bac(3),bac(1),bac(2),bac(3),rmet)
         adotb=rsdot(bab(1),bab(2),bab(3),bac(1),bac(2),bac(3),rmet)
         co=adotb/sqrt(asq*bsq)
         if( abs(co)-1.0d0 >= 0.0d0 )then
           if( abs(co)-1.0d0 <= 1.0d-12 )then
!            Allows for a small numerical inaccuracy
             thdeg=0.0d0
             if(co < 0.0d0) thdeg=180.0d0
           else
             MSG_BUG('the evaluation of the angle is wrong.')
           end if
         else
           thdeg=acos(co)*180.d0*piinv
         end if

         write(message, '(a,a,3i2,a,a,a,a,3i2,t44,a,f13.5,a)' )&
&         '  ',trim(iden(ib)),list_neighb(ineighb,2:4,1),' - ',&
&         trim(iden(ia)),' - ',trim(iden(ic)),&
&         list_neighb(jneighb,2:4,1),'bond angle is ',thdeg,' degrees '
         call wrtout(temp_unit,message,'COLL')
       end do
     end do

   end if
 end do !  End big ia loop:

 ABI_DEALLOCATE(iden)
 ABI_DEALLOCATE(list_neighb)
 ABI_DEALLOCATE(sqrlength)

 close(temp_unit)

 contains

   function rsdot(u1,u2,u3,v1,v2,v3,rmet)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rsdot'
!End of the abilint section

   real(dp) :: rsdot
   real(dp) :: u1,u2,u3,v1,v2,v3
   real(dp) :: rmet(3,3)
   rsdot=rmet(1,1)*u1*v1+rmet(2,1)*u2*v1+&
&   rmet(3,1)*u3*v1+rmet(1,2)*u1*v2+rmet(2,2)*u2*v2+&
&   rmet(3,2)*u3*v2+rmet(1,3)*u1*v3+rmet(2,3)*u2*v3+rmet(3,3)*u3*v3
 end function rsdot

end subroutine bonds_lgth_angles
!!***
