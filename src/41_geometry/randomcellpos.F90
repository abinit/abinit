!{\src2tex{textfont=tt}}
!!****f* ABINIT/randomcellpos
!! NAME
!!  randomcellpos
!!
!! FUNCTION
!!  FIXME: This subroutine creates a unit cell with random atomic positions. It is
!!         assumed that the cell parameters are given and fixed. Several methods are
!!        used to generate the cell.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2018 ABINIT group (AHR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! natom=number of atoms
!! npsp=number of pseudopotentials (needed for the dimension of znucl)
!! ntypat=number of type of atoms
!! random_atpos=input variable
!!   0 no generation of random atomic potision
!!   1 completely random atomic potisions
!!   2 random atomic positions, avoiding too close atoms (prevent coming closer than a fraction of the sum of covalent radii)
!!   3 same than 2 but also generates the rprim and acell randomly within some given ranges (angles between 50 and 130)
!! ratsph(1:ntypat)=radius of the atomic sphere
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! typat(1:natom)= input variable giving the type of each atom
!! znucl(1:npsp)=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      atomdata_from_znucl,mkrdim,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine randomcellpos(natom,npsp,ntypat,random_atpos,ratsph,rprim,rprimd,typat,xred,znucl,acell)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_atomdata

 use m_numeric_tools,  only : uniformrandom

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'randomcellpos'
 use interfaces_41_geometry, except_this_one => randomcellpos
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,npsp,ntypat,random_atpos
!arrays
 integer, intent(in)   :: typat(natom)
 real(dp),intent(in)   :: ratsph(ntypat)
 real(dp), intent(inout)  :: rprim(3,3)
 real(dp), intent(inout)  :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(in) :: znucl(npsp)
 real(dp), intent(inout) :: acell(3)

!Local variables-------------------------------
 integer ::   iatom=0,ii,idum=-20
 real(dp) ::  rij(3), rijd(3), radiuscovi, radiuscovj, dist, rati, ratj, angdeg(3)
 real(dp) ::  cosang,aa,cc,a2
 character(len=500) :: message
 type(atomdata_t) :: atom

! *************************************************************************

!DEBUG
!For the time being, print rprimd to keep it as an argument, in spite of abirule checking.
!write (std_out,*) ' randomcellpos : enter'
!write(std_out,*)' rprimd=',rprimd
!write(std_out,*)' znucl=',znucl
!write(std_out,*)' typat=',typat
!write(std_out,*)' random_atpos=',random_atpos
!ENDDEBUG

 if(random_atpos==2 .and. npsp/=ntypat)then
   write(message, '(a,i5,2a,i5,a,i5,4a)' )&
&   'Input variable random_atpos= ',random_atpos,ch10,&
&   'However, the number of pseudopotentials ',npsp,', is not equal to the number of type of atoms ',ntypat,ch10,&
&   'The use of alchemical mixing cannot be combined with the constraint based on the mixing of covalent radii.',ch10,&
&   'Action : switch to another value of random_atpos.'
   MSG_ERROR(message)
 end if

!random_atpos = 0   Default value, no random initialisation
!random_atpos = 1   Fully random (Is it really useful ???)
!random_atpos = 2   Random, but the sum of the two covalent radii is
!less than the interatomic distance
!random_atpos = 3   Random, but the sum of the two (other type of)
!radii is less than the interatomic distance
!random_atpos = 4   Random, but the sum of the two pseudopotential
!radii is less than the interatomic distance
!random_atpos = 5   Random, but the interatomic distance must be bigger
!than the sum of
!some input variable (well, instead of defining a new variable, why
!not use ratsph ?)
!Right now we are not using a factor for the tested distance.. something to be done, after a new variable has been defined

 if (random_atpos /= 0) then
   select case (random_atpos)
   case (1)
     do ii=1,natom
       xred(1,ii)=uniformrandom(idum)
       xred(2,ii)=uniformrandom(idum)
       xred(3,ii)=uniformrandom(idum)
     end do
   case (2)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       call atomdata_from_znucl(atom,znucl(typat(iatom)))
       radiuscovi = atom%rcov
       do ii=1,iatom-1
         rij=xred(:,iatom)-xred(:,ii)
!          periodic boundary conditions
         rij = rij - 0.5
         rij = rij - anint (rij)
!          coming back to cube between (0,1)
         rij = rij + 0.5
!          convert reduced coordinates to cartesian coordinates
         call xred2xcart(1,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         call atomdata_from_znucl(atom,znucl(typat(ii)))
         radiuscovj = atom%rcov
         if (dist<(radiuscovj+radiuscovi)) then
           iatom = iatom -1
           EXIT
         end if
       end do
       if (iatom>=natom) EXIT
     end do
   case(3)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       call atomdata_from_znucl(atom,znucl(typat(iatom)))
       radiuscovi = atom%rcov
       do ii=1,iatom-1
         rij=xred(:,iatom)-xred(:,ii)
!          periodic boundary conditions
         rij = rij - 0.5
         rij = rij - anint (rij)
!          coming back to cube between (0,1)
         rij = rij + 0.5
!          convert reduced coordinates to cartesian coordinates
         call xred2xcart(1,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         call atomdata_from_znucl(atom,znucl(typat(ii)))
         radiuscovj = atom%rcov
         if (dist<(radiuscovj+radiuscovi)) then
           iatom = iatom -1
           EXIT
         end if
       end do
       if (iatom>=natom) EXIT
     end do
     do ii=1,3
!        generates cells with angles between 60 and 120 degrees
       angdeg(ii)=60_dp+uniformrandom(idum)*60.0_dp
     end do
     if (angdeg(1)+angdeg(2)+angdeg(3)>360._dp) then
       angdeg(3)=360._dp-angdeg(1)-angdeg(2)
     end if
!      check if angles are between the limits and create rprim
     if( abs(angdeg(1)-angdeg(2))<tol12 .and. &
&     abs(angdeg(2)-angdeg(3))<tol12 .and. &
&     abs(angdeg(1)-90._dp)+abs(angdeg(2)-90._dp)+abs(angdeg(3)-90._dp)>tol12 )then
!        Treat the case of equal angles (except all right angles) :
!        generates trigonal symmetry wrt third axis
       cosang=cos(pi*angdeg(1)/180.0_dp)
       a2=2.0_dp/3.0_dp*(1.0_dp-cosang)
       aa=sqrt(a2)
       cc=sqrt(1.0_dp-a2)
       rprim(1,1)=aa        ; rprim(2,1)=0.0_dp                 ; rprim(3,1)=cc
       rprim(1,2)=-0.5_dp*aa ; rprim(2,2)= sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,2)=cc
       rprim(1,3)=-0.5_dp*aa ; rprim(2,3)=-sqrt(3.0_dp)*0.5_dp*aa ; rprim(3,3)=cc
!        DEBUG
!        write(std_out,*)' ingeo : angdeg=',angdeg(1:3)
!        write(std_out,*)' ingeo : aa,cc=',aa,cc
!        ENDDEBUG
     else
!        Treat all the other cases
       rprim(:,:)=0.0_dp
       rprim(1,1)=1.0_dp
       rprim(1,2)=cos(pi*angdeg(3)/180.0_dp)
       rprim(2,2)=sin(pi*angdeg(3)/180.0_dp)
       rprim(1,3)=cos(pi*angdeg(2)/180.0_dp)
       rprim(2,3)=(cos(pi*angdeg(1)/180.0_dp)-rprim(1,2)*rprim(1,3))/rprim(2,2)
       rprim(3,3)=sqrt(1.0_dp-rprim(1,3)**2-rprim(2,3)**2)
     end if
!      generate acell
     aa=zero
     do ii=1,npsp
       aa=znucl(ii)
     end do
     do ii=1,3
       acell(ii)=aa+uniformrandom(idum)*4.0
     end do
     call mkrdim(acell,rprim,rprimd)
   case(4)
     write(std_out,*) 'Not implemented yet'
   case(5)
     iatom=0
     do
       iatom=iatom+1
       xred(1,iatom)=uniformrandom(idum)
       xred(2,iatom)=uniformrandom(idum)
       xred(3,iatom)=uniformrandom(idum)
       rati=ratsph(typat(iatom))
       do ii=1,iatom-1
         ratj=ratsph(typat(ii))
!          apply periodic boundary conditions
         rij=(xred(:,iatom)-xred(:,ii))-0.5
         rij = rij - ANINT ( rij )
         rij = rij + 0.5
         call xred2xcart(natom,rprimd,rijd,rij)
         dist=dot_product(rijd,rijd)
         if (dist<(rati+ratj)) EXIT
       end do
       if (iatom==natom) EXIT
       if (ii<(iatom-1)) iatom=iatom-1
     end do
   end select
 end if

end subroutine randomcellpos
!!***
