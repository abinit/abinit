!{\src2tex{textfont=tt}}
!!****f* ABINIT/ddb_hybrid
!!
!! NAME
!! ddb_hybrid
!!
!! FUNCTION
!! Modify selected elements in atmfrc, dielt, zeff
!! as specified in an input file named: "modifs.ddb"
!! If this file does not exist, return immediately
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (PhG,XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell = lengths of lattice vectors
!! asr = flag to impose the acoustic sum rule
!! atmfrc(2,3,natom,3,natom,nrpt)= Analytical part of the
!!     Interatonmic Force Constants in real space.
!! dielt(3,3)=dielectric tensor
!! dipdip = flag to include dipole-dipole contribution
!! dyew = full Ewald matrix
!! dyewq0(3,3,natom)=contraction of the Ewald matrix at q=0
!!  modifications to the IFCs.
!! gmet = reciprocal space metric
!! gprim = reciprocal lattice vectors
!! iout = unit number for main output file. If -1 we are in a child mpi process which should not write
!! natom=number of atom in unit cell
!! nrpt = number of points in real space for the FT on the dynamical matrices
!! rcan = canonical positions of atoms
!! rmet = real-space metric
!! rprim = unit cell vectors
!! rpt = positions of points in real space for the FT on the dynamical matrices
!! ucvol = unit cell volume
!! wghatm = wieghts for pairs of atoms, in FT interpolations of dyn mat
!! xred = reduced positions of atoms
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement.
!!
!! OUTPUT
!! atmfrc(2,3,natom,3,natom,nrpt)= modified
!! dielt(3,3)=modified
!! zeff(3,3,natom)=modified
!!
!! PARENTS
!!
!! CHILDREN
!!      asrif9,canct9,ewald9,matr3inv,q0dy3_calc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ddb_hybrid(acell,asr,atmfrc,dielt,dipdip,dyew,dyewq0,&
& gmet,gprim,iout,natom,nrpt,rcan,rmet,&
& rprim,rpt,ucvol,wghatm,xred,zeff)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools, only : open_file
 use m_dynmat,   only : q0dy3_calc, asrif9, canct9
 use m_ewald,    only : ewald9

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hybrid'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: asr,dipdip,iout,natom,nrpt
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: acell(3),gmet(3,3),gprim(3,3),rcan(3,natom),rmet(3,3)
 real(dp),intent(in) :: rprim(3,3),rpt(3,nrpt),wghatm(natom,natom,nrpt)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt),dielt(3,3)
 real(dp),intent(inout) :: dyew(2,3,natom,3,natom),dyewq0(3,3,natom)
 real(dp),intent(inout) :: zeff(3,3,natom)

!Local variables-------------------------------
!Allocate should be used, instead of these fixed values for msrd and mddd
!scalars
 integer,parameter :: hyun=12,mddd=2,msrd=5000
 integer :: chk,d1,d2,flgr,genat1,genat2,ia,ib,id1,index,irpt,isrd,jj,kk,mu
 integer :: nbr,nddd,nsrd,nu,option,sumg0
 real(dp) :: dd,dd2,detdlt,detdlt2,normr2,rsq,rsq2
 logical :: ex
 character(len=500) :: msg
!arrays
 integer :: at1(msrd),at2(msrd),sta(natom,natom,nrpt)
 real(dp) :: del(3),delta(3),delta2(3),dielt2(3,3),difr(3),dyew2q0(3,3,natom)
 real(dp) :: ewab(3,3),ewab2(3,3),ewiaf(3,3),ewiaf2(3,3),ifcsr(3,3,msrd)
 real(dp) :: invdlt(3,3),invdlt2(3,3),qpt(3),ra(3),rb(3)
 real(dp) :: rpt2(3,msrd),xreda(3),zeff2(3,3,natom)

! ******************************************************************

#if defined HAVE_OS_MACOSX
 return ! macOSX seem to have problem with the inquire statement below
#endif

!Do the modifications only if the 'modifs.ddb' file exist.
 inquire(file='modifs.ddb',exist=ex)
 if(.not.ex)return

!PART 1: read input file for modifications
!+++++++++++++++++++++++++++++++++++++++++

 if (iout > 0) then
   write(iout,*)
   write(iout,*)
   write(iout,*)' **WARNING** some interatomic force constants'
   write(iout,*)'             will be artificially modified'
 end if

 if (open_file('modifs.ddb',msg,unit=hyun,status='old',form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 read(hyun,*) nsrd
 if (nsrd<0)then
   write(std_out,'(a,i8,a)' )&
&   'ddb_hybrid : nsrd is',nsrd,', which is lower than 0 => stop'
   MSG_ERROR("Aborting now")
 end if
 if (nsrd>msrd)then
   write(std_out,'(a,i8,a,a,i8,a)' )&
&   ' ddb_hybrid : nsrd is',nsrd,', which is larger than',ch10,&
&   msrd,' , the maximum number allowed => stop'
   MSG_ERROR("Aborting now")
 end if
 if (nsrd/=0)then
   isrd=0
   do ia=1,natom
     read(hyun,*) genat1
     do ib=1,(nsrd/natom)
       isrd=isrd+1
       at1(isrd)=genat1
       read(hyun,*) genat2
       at2(isrd)=genat2
       read(hyun,*) rpt2(1:3,isrd)
       read(hyun,*) ifcsr(1:3,1,isrd)
       read(hyun,*) ifcsr(1:3,2,isrd)
       read(hyun,*) ifcsr(1:3,3,isrd)
     end do
   end do
 end if

 read(hyun,*) nddd
 if (nddd<0)then
   write(std_out,'(a,i8,a)' )&
&   ' ddb_hybrid : nddd is',nddd,', which is lower than 0 => stop'
   MSG_ERROR("Aborting now")
 end if
 if (nddd>mddd)then
   write(std_out,'(a,i8,a,a,a)' )&
&   ' ddb_hybrid : nddd is',nddd,', which is not',ch10,&
&   '  one of the allowed values (0-1-2) => stop'
   MSG_ERROR("Aborting now")
 end if
 if (nddd/=0)then
   do ia=1,natom
     do id1=1,3
       read(hyun,*) zeff2(id1,1:3,ia)
     end do
   end do
   do id1=1,3
     read(hyun,*) dielt2(id1,1),dielt2(id1,2),dielt2(id1,3)
   end do
 end if

 close(hyun)

!PART 2: include modifications in the SR part
!++++++++++++++++++++++++++++++++++++++++++++

 if (nsrd/=0) then

!  write(iout,*)' List of SR modifications:'
!  write(iout,*)' dir, at, dir, at, R-vec, old-value, new-value'
!  write(iout,*)' old-value, new-value, weight'

   chk=0
   nbr=0

   do ia=1,natom
     do ib=1,natom
       do irpt=1,nrpt
         sta(ia,ib,irpt)=0
         if (wghatm(ia,ib,irpt)/=0) sta(ia,ib,irpt)=1
       end do
     end do
   end do

!  BIG loop on the N SR-modifications
   do isrd=1,nsrd

!    identify the R-vector
     flgr=0
     do irpt=1,nrpt
       difr(1:3)= abs(rpt2(1:3,isrd)-rpt(1:3,irpt))
       if (difr(1)<1.d-4) then
         if (difr(2)<1.d-4) then
           if (difr(3)<1.d-4) then
             flgr=irpt
           end if
         end if
       end if
     end do

     if (flgr==0)then
       write(std_out,'(a,a,a,i8,a)' )&
&       ' ddb_hybrid : not able to identify the',ch10,&
&       ' R-vector associated to SR-data',isrd,'=> stop'
       cycle
     end if

!    include the modifications in atmfrc
     if (wghatm(at1(isrd),at2(isrd),flgr)/=0.0_dp) then
       do d1=1,3
         do d2=1,3
           atmfrc(1,d1,at1(isrd),d2,at2(isrd),flgr) &
&           = ifcsr(d1,d2,isrd)/wghatm(at1(isrd),at2(isrd),flgr)
         end do
       end do
       sta(at1(isrd),at2(isrd),flgr)=0
       nbr=nbr+1
     end if

     chk=chk+1

   end do

!  end if nsrd><0
 end if

 if (iout > 0) then
   write(iout,*)' Number of SR modifications:',chk
   write(iout,*)' Modifications really included:',nbr

   write(iout,*)' Eventual ifc missing:'
   write(iout,*)' (ia,ib,irpt -- rpt1,rpt2,rpt3)'
   chk=0
   do ia=1,natom
     do ib=1,natom
       do irpt=1,nrpt
         if ((wghatm(ia,ib,irpt)/=0).and.(sta(ia,ib,irpt)==1)) then
           write(iout,*) ia,ib,irpt
           write(iout,*) rpt(1,irpt),rpt(2,irpt),rpt(2,irpt)
           chk=1
         end if
       end do
     end do
   end do
   if (chk==0) write(iout,*)' -no problem detected-'
 end if


!PART 3: include modifications in the DD part
!++++++++++++++++++++++++++++++++++++++++++++

 if (nddd/=0) then

   if (iout > 0) then
     write(iout,*)' The DD interaction has also been modified:'
     write(iout,*)' The Born effective chages are now:'
     do ia=1,natom
       write(iout,'(a,i4)' )' atom',ia
       do id1=1,3
         write(iout,*)zeff2(id1,1:3,ia)
       end do
     end do
     write(iout,*)' The dielectric tensor is now:'
     do id1=1,3
       write(iout,*) dielt2(id1,1:3)
     end do
   end if

!  The former values are replaced by the new ones in part 4

!  Modify dyew2q0 accordingly to zeff2 and dielt2 (if dip-dip non-zero)
   if (dipdip==1) then
     sumg0=0
     qpt(1)=0.0_dp
     qpt(2)=0.0_dp
     qpt(3)=0.0_dp
     call ewald9(acell,dielt2,dyew,gmet,gprim,natom,&
&     qpt,rmet,rprim,sumg0,ucvol,xred,zeff2)
     if (asr==1.or.asr==2) then
       option=asr
       call q0dy3_calc(natom,dyew2q0,dyew,option)
     else if (asr==0) then
       dyew2q0(:,:,:)=0.0_dp
     end if
   end if

!  end if nddd/=0
 end if

!Eventually keep the total ifc within the box unchanged
!This basically corresponds to modify the SR part accordingly
!to the change of the DD part within the box...

 if (nddd==2) then

!  Store the interatomic distances
!  call dist9(acell,dist,gprim,natom,nrpt,rcan,rprim,rpt)

!  calculating the inverse (transpose) of the dielectric tensor
   call matr3inv(dielt,invdlt)
   call matr3inv(dielt2,invdlt2)
!  calculating the determinant of the dielectric tensor
   detdlt=dielt(1,1)*dielt(2,2)*dielt(3,3)+dielt(1,3)*dielt(2,1)*&
&   dielt(3,2)+dielt(1,2)*dielt(2,3)*dielt(3,1)-dielt(1,3)*&
&   dielt(2,2)*dielt(3,1)-dielt(1,1)*dielt(2,3)*dielt(3,2)-&
&   dielt(1,2)*dielt(2,1)*dielt(3,3)
   detdlt2=dielt2(1,1)*dielt2(2,2)*dielt2(3,3)+dielt2(1,3)*&
&   dielt2(2,1)*&
&   dielt2(3,2)+dielt2(1,2)*dielt2(2,3)*dielt2(3,1)-&
&   dielt2(1,3)*&
&   dielt2(2,2)*dielt2(3,1)-dielt2(1,1)*dielt2(2,3)*&
&   dielt2(3,2)-&
&   dielt2(1,2)*dielt2(2,1)*dielt2(3,3)

!  Big loop on first atom ia
   do ia=1,natom

!    First transform canonical coordinates to reduced coordinates
     xreda(:)=gprim(1,:)*rcan(1,ia)+gprim(2,:)*rcan(2,ia)&
&     +gprim(3,:)*rcan(3,ia)

!    Then to cartesian coordinates
     ra(:)=xreda(1)*acell(1)*rprim(:,1)+&
&     xreda(2)*acell(2)*rprim(:,2)+&
&     xreda(3)*acell(3)*rprim(:,3)

!    Big intra-loop on the atoms in the box ib
     do index=1,(natom*nrpt)

       call canct9(acell,gprim,ib,index,irpt,natom,nrpt,&
&       rcan,rb,rprim,rpt)

       if (wghatm(ia,ib,irpt)/=0)then

         del(:)=ra(:)-rb(:)
         rsq=0.0_dp
         rsq2=0.0_dp
         delta(:)=0.0_dp
         delta2(:)=0.0_dp
         do jj=1,3
           do kk=1,3
             ewab(jj,kk)=0.0_dp
             ewab2(jj,kk)=0
             rsq=rsq+del(jj)*invdlt(kk,jj)*del(kk)
             rsq2=rsq2+del(jj)*invdlt2(kk,jj)*del(kk)
             delta(kk)=delta(kk)+invdlt(kk,jj)*del(jj)
             delta2(kk)=delta2(kk)+invdlt2(kk,jj)*del(jj)
           end do
         end do
         dd=sqrt(rsq)
         dd2=sqrt(rsq2)

!        Avoid zero denominators in 'term':
         if (sqrt(rsq)>=1.0d-12) then
           do mu=1,3
             do nu=1,3
               ewab(mu,nu)=(-3*delta(nu)*delta(mu)+invdlt(nu,mu)*dd**2)&
&               /dd**5/sqrt(detdlt)
             end do
           end do
         else
           if (ia/=ib)then
             write(std_out,*)' ddb_hybrid : interatomic distance vanishes ',&
&             ' Check the input for the following atoms :'
             write(std_out,*)ia,ib
             MSG_ERROR("Aborting now")
           end if
         end if

!        Avoid zero denominators in 'term':
         if (sqrt(rsq)>=1.0d-12) then
           do mu=1,3
             do nu=1,3
               ewab2(mu,nu)=&
&               (-3*delta2(nu)*delta2(mu)+invdlt2(nu,mu)*dd2**2)&
&               /dd2**5/sqrt(detdlt2)
             end do
           end do
         else
           if (ia/=ib)then
             write(std_out,*)' ddb_hybrid : inter-atomic distance vanishes ',&
&             ' Check the input for the following atoms :'
             write(std_out,*)ia,ib
             MSG_ERROR("Aborting now")
           end if
         end if

!        Take into account the effective charge tensor
         normr2=rpt(1,irpt)**2+rpt(2,irpt)**2+rpt(3,irpt)**2
         do mu=1,3
           do nu=1,3
             ewiaf(mu,nu)=0.0_dp
             ewiaf2(mu,nu)=0.0_dp
             if((ia==ib).and.(normr2<1.d-7))then
               ewiaf(mu,nu)=-dyewq0(mu,nu,ia)
               ewiaf2(mu,nu)=-dyew2q0(mu,nu,ia)
             end if
             do jj=1,3
               do kk=1,3
                 ewiaf(mu,nu)=ewiaf(mu,nu)&
&                 +zeff(jj,mu,ia)*zeff(kk,nu,ib)*&
&                 ewab(jj,kk)
                 ewiaf2(mu,nu)=ewiaf2(mu,nu)&
&                 +zeff2(jj,mu,ia)*zeff2(kk,nu,ib)*&
&                 ewab2(jj,kk)
               end do
             end do
           end do
         end do

!        add DD(old)-DD(new) to the SR part...
         do mu=1,3
           do nu=1,3
             atmfrc(1,mu,ia,nu,ib,irpt)=atmfrc(1,mu,ia,nu,ib,irpt)&
&             +(ewiaf(mu,nu)-ewiaf2(mu,nu))/wghatm(ia,ib,irpt)
           end do
         end do

!        end if wghatm><0
       end if

!      End loop ia
     end do

!    End loop index
   end do

!  end if nddd=2
 end if

!PART 4: actualize the matrices before exiting
!+++++++++++++++++++++++++++++++++++++++++++++

 if (nddd>0) then

!  Update of dielt and zeff
   zeff(:,:,:)=zeff2(:,:,:)
   dielt(:,:)=dielt2(:,:)

!  Update of dyewq0
   dyewq0(:,:,:)=dyew2q0(:,:,:)

 end if

 if ((nsrd>0).or.(nddd==2)) then
!  Reimpose the ASR on the ifc that have been modified...
   if(asr>0)then
     write(std_out,*)' ddb_hybrid : enter asrif9 '
     call asrif9(asr,atmfrc,natom,nrpt,rpt,wghatm)
     write(std_out,*)' ddb_hybrid : exit asrif9 '
   end if
 end if

end subroutine ddb_hybrid
!!***
