!!****m* ABINIT/m_xfpack
!! NAME
!!  m_xfpack
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2025 ABINIT group (XG, MJV, DCA, GMR, JCC, SE)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xfpack

 use defs_basis
 use m_errors
 use m_abicore
 use m_abimover

 use m_matrix,     only : matr3inv
 use m_geometry,   only : mkradim, mkrdim, metric, strainsym
 use m_results_gs , only : results_gs_type
 use m_bfgs,        only : hessupdt

 implicit none

 private
!!***

 public :: xfpack_vin2x
 public :: xfpack_x2vin
 public :: xfpack_f2vout
 public :: xfh_recover_new
 public :: xfh_update
!!***

contains
!!***

!!****f* ABINIT/xfpack_vin2x
!! NAME
!! xfpack_vin2x
!!
!! FUNCTION
!! Old option=2, transfer vin  to xred, acell and rprim
!!
!! INPUTS
!! acell0(3)=reference length scales of primitive translations (bohr), needed for some values of optcell.
!! natom=number of atoms in cell
!! ndim=dimension of vin array
!! nsym=order of group.
!! rprimd0(3,3)=reference real space primitive translations,
!!   needed for some values of optcell.
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of acell and rprim
!!  are contained in vin.
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!! ucvol0=reference unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! acell(3)=length scales of primitive translations (bohr)
!! rprim(3,3)=dimensionless real space primitive translations
!! vin(ndim)=vector that contains xred and some quantity derived
!!   from acell and rprim, depending on the value of optcell.
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SOURCE

subroutine xfpack_vin2x(acell,acell0,natom,ndim,nsym,optcell,&
& rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,nsym,optcell
 real(dp),intent(in) :: ucvol0
 real(dp),intent(out) :: ucvol
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: acell0(3),rprimd0(3,3)
 real(dp),intent(inout) :: acell(3),rprim(3,3)
 real(dp),intent(in) :: vin(ndim)
 real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: scale
 character(len=500) :: message
 logical :: equal=.TRUE.
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: rprimd_symm(3,3),scaling(3,3)

! *************************************************************************

!!DEBUG
!write(ab_out,*) ''
!write(ab_out,*) 'xfpack_vin2x'
!write(ab_out,*) 'natom=',natom
!write(ab_out,*) 'ndim=',ndim
!write(ab_out,*) 'nsym=',nsym
!write(ab_out,*) 'optcell=',optcell
!write(ab_out,*) 'ucvol=',ucvol
!write(ab_out,*) 'xred='
!do kk=1,natom
!write(ab_out,*) xred(:,kk)
!end do
!write(ab_out,*) 'VECTOR INPUT (vin) xfpack_vin2x INPUT'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!!DEBUG


!##########################################################
!### 1. Test for compatible ndim

 if(optcell==0 .and. ndim/=3*natom)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=0, ndim MUST be equal to 3*natom,',ch10,&
&   '  while ndim=',ndim,' and 3*natom=',3*natom,'.'
   ABI_BUG(messagE)
 end if

 if( (optcell==1) .and. ndim/=3*natom+1)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=1 ndim MUST be equal to 3*natom+1,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+1=',3*natom+1,'.'
   ABI_BUG(message)
 end if
 
 if( (optcell==2 .or. optcell==3) &
& .and. ndim/=3*natom+6)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=2 or 3, ndim MUST be equal to 3*natom+6,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+6=',3*natom+6,'.'
   ABI_BUG(message)
 end if

 if( optcell>=4 .and. ndim/=3*natom+3)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=4,5,6,7,8 or 9, ndim MUST be equal to 3*natom+3,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+3=',3*natom+3,'.'
   ABI_BUG(message)
 end if

!##########################################################
!### 3. option=2, transfer vin  to xred, acell and rprim

!Get xred, and eventually acell and rprim from vin
 xred(:,:)=reshape( vin(1:3*natom), (/3,natom/) )

 if(optcell==1)then

!  acell(:)=acell0(:)*vin(3*natom+1)/(ucvol0**third)
   acell(:)=acell0(:)*vin(3*natom+1)

 else if (optcell>=2)then

   scaling(:,:)=0.0_dp
   scaling(1,1)=1.0_dp ; scaling(2,2)=1.0_dp ; scaling(3,3)=1.0_dp

   if(optcell==2 .or. optcell==3)then
     scaling(1,1)=vin(3*natom+1)
     scaling(2,2)=vin(3*natom+2)
     scaling(3,3)=vin(3*natom+3)
     scaling(2,3)=vin(3*natom+4) ; scaling(3,2)=vin(3*natom+4)
     scaling(1,3)=vin(3*natom+5) ; scaling(3,1)=vin(3*natom+5)
     scaling(1,2)=vin(3*natom+6) ; scaling(2,1)=vin(3*natom+6)
   else if(optcell==4)then
     scaling(1,1)=vin(3*natom+1)
     if (abs(scaling(1,1) - 1.0_dp) <1.E-14) scaling(1,1)=1.0_dp
     scaling(3,1)=vin(3*natom+2)  
     if (abs(scaling(3,1)) <1.E-14) scaling(3,1)=0.0_dp
     scaling(2,1)=vin(3*natom+3)
     if (abs(scaling(2,1)) <1.E-14) scaling(2,1)=0.0_dp
   else if(optcell==5)then
     scaling(2,2)=vin(3*natom+1)
     if (abs(scaling(2,2) - 1.0_dp) <1.E-14) scaling(2,2)=1.0_dp
     scaling(3,2)=vin(3*natom+2)  
     if (abs(scaling(3,2)) <1.E-14) scaling(3,2)=0.0_dp
     scaling(1,2)=vin(3*natom+3)
     if (abs(scaling(1,2)) <1.E-14) scaling(1,2)=0.0_dp
   else if(optcell==6)then
     scaling(3,3)=vin(3*natom+1)
     if (abs(scaling(3,3) - 1.0_dp) <1.E-14) scaling(3,3)=1.0_dp
     scaling(2,3)=vin(3*natom+2)  
     if (abs(scaling(2,3)) <1.E-14) scaling(2,3)=0.0_dp
     scaling(1,3)=vin(3*natom+3)
     if (abs(scaling(1,3)) <1.E-14) scaling(1,3)=0.0_dp
   else if(optcell==7)then
     scaling(2,2)=vin(3*natom+2) ; scaling(3,3)=vin(3*natom+3)
     scaling(2,3)=vin(3*natom+1) ; scaling(3,2)=vin(3*natom+1)
   else if(optcell==8)then
     scaling(1,1)=vin(3*natom+1) ; scaling(3,3)=vin(3*natom+3)
     scaling(1,3)=vin(3*natom+2) ; scaling(3,1)=vin(3*natom+2)
   else if(optcell==9)then
     scaling(1,1)=vin(3*natom+1) ; scaling(2,2)=vin(3*natom+2)
     scaling(1,2)=vin(3*natom+3) ; scaling(2,1)=vin(3*natom+3)
   end if
   if(optcell<=3 .or. optcell>=7)then
     do ii=1,3
       do jj=1,3
         rprimd(ii,jj)=0.0_dp
         do kk=1,3
           rprimd(ii,jj)=rprimd(ii,jj)+scaling(ii,kk)*rprimd0(kk,jj)
         end do
       end do
     end do
   ! for optcell=4,5,6, implementing search for all 3 components of the vector to be relaxed according to Eq.10 of J. Chem. Phys.
   ! 136, 074103 (2012), i.e. search direction given by rprimd0 * stress
   else if(optcell==4)then
     rprimd(:,2) = rprimd0(:,2)
     rprimd(:,3) = rprimd0(:,3)
     rprimd(:,1) = 0.0_dp
     do ii=1,3
       do kk=1,3
          rprimd(ii,1) = rprimd(ii,1) + scaling(kk,1)*rprimd0(ii,kk)
       end do
     end do
   else if(optcell==5)then
     rprimd(:,1) = rprimd0(:,1)
     rprimd(:,3) = rprimd0(:,3)
     rprimd(:,2) = 0.0_dp
     do ii=1,3
       do kk=1,3
          rprimd(ii,2) = rprimd(ii,2) + scaling(kk,2)*rprimd0(ii,kk)
       end do
     end do
   else if(optcell==6)then
     rprimd(:,1) = rprimd0(:,1)
     rprimd(:,2) = rprimd0(:,2)
     rprimd(:,3) = 0.0_dp
     do ii=1,3
       do kk=1,3
          rprimd(ii,3) = rprimd(ii,3) + scaling(kk,3)*rprimd0(ii,kk)
       end do
     end do
   end if

!  Rescale if the volume must be preserved
   if(optcell==3)then
     call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
     scale=(ucvol0/ucvol)**third
     rprimd(:,:)=scale*rprimd(:,:)
   end if
   call strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
   do jj=1,3
     do ii=1,3
!      write(ab_out,*) 'DIFF',ii,jj,abs(rprimd0(ii,jj)-rprimd_symm(ii,jj))
       if (abs(rprimd0(ii,jj)-rprimd_symm(ii,jj))>1.E-14)&
&       equal=.FALSE.
     end do
   end do

   if (equal)then
     acell(:)=acell0(:)
     rprimd(:,:)=rprimd0(:,:)
   else
!    Use a representation based on normalised rprim vectors
     call mkradim(acell,rprim,rprimd_symm)
   end if

 end if

end subroutine xfpack_vin2x
!!***

!!****f* ABINIT/xfpack_x2vin
!! NAME
!! xfpack_x2vin
!!
!! FUNCTION
!! Old option=1, transfer xred, acell, and rprim to vin
!!
!! INPUTS
!! acell0(3)=reference length scales of primitive translations (bohr), needed for some values of optcell.
!! natom=number of atoms in cell
!! ndim=dimension of vin arrays
!! nsym=order of group.
!! rprimd0(3,3)=reference real space primitive translations,
!!   needed for some values of optcell.
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of acell and rprim
!!  are contained in vin.
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!! ucvol0=reference unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! acell(3)=length scales of primitive translations (bohr)
!! rprim(3,3)=dimensionless real space primitive translations
!! vin(ndim)=vector that contains xred and some quantity derived
!!   from acell and rprim, depending on the value of optcell.
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! SOURCE

subroutine xfpack_x2vin(acell,natom,ndim,nsym,optcell,&
  & rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,nsym,optcell
 real(dp),intent(in) :: ucvol0
 real(dp),intent(inout) :: ucvol !vz_i
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd0(3,3)
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vin(ndim)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: scale
 character(len=500) :: message
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),gprimd0(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: rprimd_symm(3,3),scaling(3,3)

! *************************************************************************

!!DEBUG
!write(ab_out,*) ''
!write(ab_out,*) 'xfpack_x2vin'
!write(ab_out,*) 'natom=',natom
!write(ab_out,*) 'ndim=',ndim
!write(ab_out,*) 'nsym=',nsym
!write(ab_out,*) 'optcell=',optcell
!write(ab_out,*) 'ucvol=',ucvol
!write(ab_out,*) 'xred='
!do kk=1,natom
!write(ab_out,*) xred(:,kk)
!end do
!write(ab_out,*) 'VECTOR INPUT (vin) xfpack_x2vin INPUT'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!!DEBUG


!##########################################################
!### 1. Test for compatible ndim

 if(optcell==0 .and. ndim/=3*natom)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=0, ndim MUST be equal to 3*natom,',ch10,&
&   '  while ndim=',ndim,' and 3*natom=',3*natom,'.'
   ABI_BUG(message)
 end if

 if( optcell==1 .and. ndim/=3*natom+1)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=1, ndim MUST be equal to 3*natom+1,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+1=',3*natom+1,'.'
   ABI_BUG(message)
 end if

 if( (optcell==2 .or. optcell==3) &
& .and. ndim/=3*natom+6)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=2,3,4,5,6, ndim MUST be equal to 3*natom+6,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+6=',3*natom+6,'.'
   ABI_BUG(message)
 end if

 if( optcell>=4 .and. ndim/=3*natom+3)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=4,5,6,7,8 or 9, ndim MUST be equal to 3*natom+3,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+3=',3*natom+3,'.'
   ABI_BUG(message)
 end if

!##########################################################
!### 2. option=1, transfer xred, acell, and rprim to vin

!Get vin from xred, acell, and rprim
 vin(1:3*natom)= reshape(xred(:,:), (/3*natom/) )

 if(optcell/=0)then
   call mkrdim(acell,rprim,rprimd)
   call strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
   call metric(gmet,gprimd,-1,rmet,rprimd_symm,ucvol)

   if(optcell==1)then

!    vin(3*natom+1)=ucvol**third
     vin(3*natom+1)=(ucvol/ucvol0)**third

   else if(optcell>=2)then

!    Generates gprimd0
     call matr3inv(rprimd0,gprimd0)
     if (optcell==2 .or. optcell==3 .or. optcell>=7)then
       do ii=1,3
         do jj=1,3
           scaling(ii,jj)=0.0_dp
           do kk=1,3
             scaling(ii,jj)=scaling(ii,jj)+rprimd_symm(ii,kk)*gprimd0(jj,kk)
           end do
         end do
       end do
!    Rescale if the volume must be preserved
       if(optcell==3)then
         scale=(ucvol0/ucvol)**third
         scaling(:,:)=scale*scaling(:,:)
       end if
       if(optcell==2 .or. optcell==3)then
         vin(3*natom+1)=scaling(1,1) ; vin(3*natom+4)=(scaling(2,3)+scaling(3,2))*0.5_dp
         vin(3*natom+2)=scaling(2,2) ; vin(3*natom+5)=(scaling(1,3)+scaling(3,1))*0.5_dp
         vin(3*natom+3)=scaling(3,3) ; vin(3*natom+6)=(scaling(1,2)+scaling(2,1))*0.5_dp
       else if(optcell>=7)then
         vin(3*natom+1)=scaling(1,1)
         vin(3*natom+2)=scaling(2,2)
         vin(3*natom+3)=scaling(3,3)
         if(optcell==7)vin(3*natom+1)=(scaling(2,3)+scaling(3,2))*0.5_dp
         if(optcell==8)vin(3*natom+2)=(scaling(1,3)+scaling(3,1))*0.5_dp
         if(optcell==9)vin(3*natom+3)=(scaling(1,2)+scaling(2,1))*0.5_dp
       end if
     end if

     if (optcell==4)then
       scaling(:,:) = 0.0_dp
       do ii=1,3
          do kk=1,3
            scaling(ii,1) = scaling(ii,1) + gprimd0(kk,ii)*rprimd_symm(kk,1)
          end do
       end do
       vin(3*natom+1) = scaling(1,1)
       vin(3*natom+2) = scaling(3,1)
       vin(3*natom+3) = scaling(2,1)
     else if (optcell==5)then
       scaling(:,:) = 0.0_dp
       do ii=1,3
          do kk=1,3
            scaling(ii,2) = scaling(ii,2) + gprimd0(kk,ii)*rprimd_symm(kk,2)
          end do
       end do
       vin(3*natom+1) = scaling(2,2)
       vin(3*natom+2) = scaling(3,2)
       vin(3*natom+3) = scaling(1,2)
     else if (optcell==6)then
       scaling(:,:) = 0.0_dp
       do ii=1,3
          do kk=1,3
            scaling(ii,3) = scaling(ii,3) + gprimd0(kk,ii)*rprimd_symm(kk,3)
          end do
       end do
       vin(3*natom+1) = scaling(3,3)
       vin(3*natom+2) = scaling(2,3)
       vin(3*natom+3) = scaling(1,3)
     end if

   end if

 end if

end subroutine xfpack_x2vin
!!***

!!****f* ABINIT/xfpack_f2vout
!! NAME
!! xfpack_f2vout
!!
!! FUNCTION
!! Old option=3, transfer gred and strten to vout
!!
!! INPUTS
!! natom=number of atoms in cell
!! ndim=dimension of vout arrays
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of strten
!!  are contained in vout.
!! strtarget(6)=target stresses ; they will be subtracted from strten when vout
!!  is computed.
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! gred(3,natom)=grads of Etot wrt reduced coordinates (hartree)
!! strten(6)=components of the stress tensor (hartree/bohr^3)
!! vout(ndim)=vector that contains gred and some quantity derived from
!!   strten, depending on the value of optcell, and taking care ot strtarget
!!
!! SOURCE

subroutine xfpack_f2vout(gred,natom,ndim,optcell,strtarget,strten,ucvol,vout)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,optcell
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: strtarget(6)
 real(dp),intent(in) :: gred(3,natom),strten(6)
 real(dp),intent(out) :: vout(ndim)

!Local variables-------------------------------
!scalars
 real(dp) :: strdiag
 character(len=500) :: message
!arrays
 real(dp) :: dstr(6)

! *************************************************************************

!!DEBUG
!write(ab_out,*) ''
!write(ab_out,*) 'xfpack_f2vout'
!write(ab_out,*) 'natom=',natom
!write(ab_out,*) 'ndim=',ndim
!write(ab_out,*) 'optcell=',optcell
!write(ab_out,*) 'ucvol=',ucvol
!!DEBUG


!##########################################################
!### 1. Test for compatible ndim

 if(optcell==0 .and. ndim/=3*natom)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=0, ndim MUST be equal to 3*natom,',ch10,&
&   '  while ndim=',ndim,' and 3*natom=',3*natom,'.'
   ABI_BUG(message)
 end if

 if( optcell==1 .and. ndim/=3*natom+1)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=1, ndim MUST be equal to 3*natom+1,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+1=',3*natom+1,'.'
   ABI_BUG(message)
 end if

 if( (optcell==2 .or. optcell==3) &
& .and. ndim/=3*natom+6)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=2 or 3, ndim MUST be equal to 3*natom+6,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+6=',3*natom+6,'.'
   ABI_BUG(message)
 end if

 if( optcell>=4 .and. ndim/=3*natom+3)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=4,5,6,7,8 or 9, ndim MUST be equal to 3*natom+3,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+3=',3*natom+3,'.'
   ABI_BUG(message)
 end if
!
!Get vout from gred and strten
!
 vout(1:3*natom)= reshape(gred(:,:), (/3*natom/) )
 dstr(:)=strten(:)-strtarget(:)

 if(optcell==1)then

   vout(3*natom+1)=( dstr(1)+dstr(2)+dstr(3))*ucvol

 else if(optcell>=2)then
!  Eventually take away the trace
   strdiag=0.0_dp
   if(optcell==3) strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   if(optcell==2 .or. optcell==3)then
     vout(3*natom+1:3*natom+3)=(dstr(1:3)-strdiag)*ucvol
!    For non-diagonal derivatives, must take into account
!    that eps(i,j) AND eps(j,i) are varied at the same time. Thus, derivative
!    is twice larger
     vout(3*natom+4:3*natom+6)=dstr(4:6)*ucvol*2.0_dp
   else if(optcell==7 .or. optcell==8 .or. optcell==9)then
!    Similar to case optcell==2 or optcell==3, but in 2 dimensions.
     vout(3*natom+1:3*natom+3)=dstr(1:3)*ucvol
     vout(3*natom+optcell-6)  =dstr(optcell-3)*ucvol*2.0_dp
   else if (optcell==4)then
     vout(3*natom+1) = dstr(1)*ucvol
     vout(3*natom+2) = dstr(5)*ucvol
     vout(3*natom+3) = dstr(6)*ucvol
   else if (optcell==5)then
     vout(3*natom+1) = dstr(2)*ucvol
     vout(3*natom+2) = dstr(4)*ucvol
     vout(3*natom+3) = dstr(6)*ucvol
   else if (optcell==6)then
     vout(3*natom+1) = dstr(3)*ucvol
     vout(3*natom+2) = dstr(4)*ucvol
     vout(3*natom+3) = dstr(5)*ucvol
   end if

 end if

end subroutine xfpack_f2vout
!!***


!!****f* ABINIT/xfh_recover_new
!! NAME
!! xfh_recover_new
!!
!! FUNCTION
!! Update the contents of the history xfhist taking values
!! from xred, acell, rprim, gred_corrected and strten
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE


subroutine xfh_recover_new(ab_xfh,ab_mover,acell,cycl_main,gred,&
& hessin,ndim,rprim,rprimd0,strten,ucvol,ucvol0,vin,vin_prev,vout,&
& vout_prev,xred)

!Arguments ------------------------------------
!scalars

integer,intent(in) :: ndim
integer,intent(out) :: cycl_main
real(dp),intent(inout) :: ucvol,ucvol0
type(ab_xfh_type),intent(inout) :: ab_xfh
type(abimover),intent(in) :: ab_mover


!arrays
real(dp),intent(inout) :: acell(3)
real(dp),intent(inout) :: hessin(:,:)
real(dp),intent(inout) :: xred(3,ab_mover%natom)
real(dp),intent(inout) :: rprim(3,3)
real(dp),intent(inout) :: rprimd0(3,3)
real(dp),intent(inout) :: gred(3,ab_mover%natom)
real(dp),intent(inout) :: strten(6)
real(dp),intent(inout) :: vin(:)
real(dp),intent(inout) :: vin_prev(:)
real(dp),intent(inout) :: vout(:)
real(dp),intent(inout) :: vout_prev(:)

!Local variables-------------------------------
!scalars
integer :: ixfh ! kk,jj

!*********************************************************************

 if(ab_xfh%nxfh/=0)then
!  Loop over previous time steps
   do ixfh=1,ab_xfh%nxfh

!    For that time step, get new (x,f) from xfhist
     xred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom        ,1,ixfh)
     rprim(1:3,1:3)=ab_xfh%xfhist(:,ab_mover%natom+2:ab_mover%natom+4,1,ixfh)
     acell(:)      =ab_xfh%xfhist(:,ab_mover%natom+1,1,ixfh)
     gred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom,2,ixfh)
!    This use of results_gs is unusual
     strten(1:3)   =ab_xfh%xfhist(:,ab_mover%natom+2,2,ixfh)
     strten(4:6)   =ab_xfh%xfhist(:,ab_mover%natom+3,2,ixfh)

!    !DEBUG
!    write (ab_out,*) '---READED FROM XFHIST---'

!    write (ab_out,*) 'XRED'
!    do kk=1,ab_mover%natom
!    write (ab_out,*) xred(:,kk)
!    end do
!    write (ab_out,*) 'FRED'
!    do kk=1,ab_mover%natom
!    write (ab_out,*) gred(:,kk)
!    end do
!    write(ab_out,*) 'RPRIM'
!    do kk=1,3
!    write(ab_out,*) rprim(:,kk)
!    end do
!    write(ab_out,*) 'ACELL'
!    write(ab_out,*) acell(:)
!    !DEBUG

!    Transfer it in vin, vout
     call xfpack_x2vin(acell,ab_mover%natom,&
&     ndim,ab_mover%nsym,ab_mover%optcell,rprim,rprimd0,&
&     ab_mover%symrel,ucvol,ucvol0,vin,xred)
     call xfpack_f2vout(gred,ab_mover%natom,&
&     ndim,ab_mover%optcell,ab_mover%strtarget,strten,&
&     ucvol,vout)
!    Get old time step, if any, and update inverse hessian
     if(ixfh/=1)then
       xred(:,:)     =ab_xfh%xfhist(:,1:ab_mover%natom,1,ixfh-1)
       rprim(1:3,1:3)=&
&       ab_xfh%xfhist(:,ab_mover%natom+2:ab_mover%natom+4,1,ixfh-1)
       acell(:)=ab_xfh%xfhist(:,ab_mover%natom+1,1,ixfh-1)
       gred(:,:)=ab_xfh%xfhist(:,1:ab_mover%natom,2,ixfh-1)
!      This use of results_gs is unusual
       strten(1:3)=ab_xfh%xfhist(:,ab_mover%natom+2,2,ixfh-1)
       strten(4:6)=ab_xfh%xfhist(:,ab_mover%natom+3,2,ixfh-1)
!      Tranfer it in vin_prev, vout_prev
       call xfpack_x2vin(acell,ab_mover%natom,&
&       ndim,ab_mover%nsym,ab_mover%optcell,rprim,rprimd0,&
&       ab_mover%symrel,ucvol,ucvol0,vin_prev,xred)
       call xfpack_f2vout(gred,ab_mover%natom,&
&       ndim,ab_mover%optcell,ab_mover%strtarget,strten,&
&       ucvol,vout_prev)

!      write(ab_out,*) 'Hessian matrix before update',ndim,'x',ndim
!      write(ab_out,*) 'ixfh=',ixfh
!      do kk=1,ndim
!      do jj=1,ndim,3
!      if (jj+2<=ndim)then
!      write(ab_out,*) jj,hessin(jj:jj+2,kk)
!      else
!      write(ab_out,*) jj,hessin(jj:ndim,kk)
!      end if
!      end do
!      end do

       call hessupdt(hessin,ab_mover%iatfix,ab_mover%natom,ndim,&
&       vin,vin_prev,vout,vout_prev)

!      !DEBUG
!      write(ab_out,*) 'Hessian matrix after update',ndim,'x',ndim
!      do kk=1,ndim
!      do jj=1,ndim,3
!      if (jj+2<=ndim)then
!      write(ab_out,*) jj,hessin(jj:jj+2,kk)
!      else
!      write(ab_out,*) jj,hessin(jj:ndim,kk)
!      end if
!      end do
!      end do
!      !DEBUG

     end if !if(ab_xfh%nxfh/=0)
   end do ! End loop over previous time steps

!  The hessian has been generated,
!  as well as the latest vin and vout
!  so will cycle the main loop
   cycl_main=1
 end if

end subroutine xfh_recover_new
!!***

!!****f* ABINIT/xfh_update
!! NAME
!! xfh_update
!!
!! FUNCTION
!! Update the contents of the history xfhist taking values
!! from xred, acell, rprim, gred_corrected and strten
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine xfh_update(ab_xfh,acell,gred_corrected,natom,rprim,strten,xred)

!Arguments ------------------------------------
!scalars
type(ab_xfh_type),intent(inout) :: ab_xfh
integer,intent(in) :: natom

!arrays
real(dp),intent(in) :: acell(3)
real(dp),intent(in) :: xred(3,natom)
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: gred_corrected(3,natom)
real(dp),intent(in) :: strten(6)

!Local variables-------------------------------
!scalars
!integer :: kk

!*********************************************************************

!DEBUG
!write (ab_out,*) '---WROTE TO XFHIST---'

!write (ab_out,*) 'XRED'
!do kk=1,natom
!write (ab_out,*) xred(:,kk)
!end do
!write (ab_out,*) 'FRED'
!do kk=1,natom
!write (ab_out,*) gred_corrected(:,kk)
!end do
!write(ab_out,*) 'RPRIM'
!do kk=1,3
!write(ab_out,*) rprim(:,kk)
!end do
!write(ab_out,*) 'ACELL'
!write(ab_out,*) acell(:)
!DEBUG

 ab_xfh%nxfh=ab_xfh%nxfh+1

 ab_xfh%xfhist(:,1:natom,1,ab_xfh%nxfh)=xred(:,:)
 ab_xfh%xfhist(:,natom+1,1,ab_xfh%nxfh)=acell(:)
 ab_xfh%xfhist(:,natom+2:natom+4,1,ab_xfh%nxfh)=rprim(:,:)
 ab_xfh%xfhist(:,1:natom,2,ab_xfh%nxfh)=gred_corrected(:,:)
 ab_xfh%xfhist(:,natom+2,2,ab_xfh%nxfh)=strten(1:3)
 ab_xfh%xfhist(:,natom+3,2,ab_xfh%nxfh)=strten(4:6)

end subroutine xfh_update
!!***

end module m_xfpack
!!***
