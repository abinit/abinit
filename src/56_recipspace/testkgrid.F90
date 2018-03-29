!{\src2tex{textfont=tt}}
!!****f* ABINIT/testkgrid
!!
!! NAME
!! testkgrid
!!
!! FUNCTION
!! Test different grids of k points. The algorithm used is based on the idea of testing different
!! one-dimensional sets of possible k point grids. It is not exhaustive (other families could be included),
!! but should do a respectable job in all cases. The Monkhorst-Pack set of grids (defined with respect to
!! symmetry axes, and not primitive axes) is always tested.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  bravais(11): bravais(1)=iholohedry
!!               bravais(2)=center
!!               bravais(3:11)=coordinates of rprim in the axes of the conventional bravais lattice (*2 if center/=0)
!!  iout=unit number for echoed output
!!  msym=default maximal number of symmetries
!!  nsym=number of symetries
!!  prtkpt=if non-zero, will write the characteristics of k grids, then stop
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry operations in real space in terms of primitive translations
!!  vacuum(3)=for each direction, 0 if no vacuum, 1 if vacuum
!!
!! OUTPUT
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk=number of k-point shifts in shiftk (always 1 from this routine)
!!  shiftk(3,210)=shift vectors for k point generation
!!
!! SIDE EFFECTS
!!  kptrlen=length of the smallest real space supercell vector associated with the lattice of k points.
!!
!! NOTES
!! Note that nkpt can be computed by calling this routine with input value nkpt=0
!! Note that kptopt is always =1 in this routine.
!!
!! PARENTS
!!      inkpts,m_ab7_kpoints
!!
!! CHILDREN
!!      getkgrid,leave_new,matr3inv,metric,smallprim,wrtout,xmpi_abort
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine testkgrid(bravais,iout,kptrlatt,kptrlen,&
& msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,vacuum)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_geometry,     only : metric

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'testkgrid'
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace, except_this_one => testkgrid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,msym,nsym,prtkpt
 integer,intent(out) :: nshiftk
 real(dp),intent(inout) :: kptrlen
!arrays
 integer,intent(in) :: bravais(11),symafm(msym),symrel(3,3,msym),vacuum(3)
 integer,intent(out) :: kptrlatt(3,3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: shiftk(3,210) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: kptopt=1,mkpt_list=100000
 integer :: ang90,center,dirvacuum,equal,igrid,igrid_current,iholohedry,ii,init_mult,iscale,iscf
 integer :: iset,mult1,mult2,mult3,ndims,nkpt,nkpt_current,nkpt_trial,nset
 real(dp) :: buffer_scale,determinant,fact,factor,kptrlen_current,kptrlen_max,kptrlen_target
 real(dp) :: kptrlen_trial,length1,length2,length3,length_axis1,length_axis2
 real(dp) :: length_axis3,merit_factor,mult1h,mult2h,mult3h,reduceda,reducedb
 real(dp) :: sca,scb,scc,surface,ucvol
 character(len=500) :: message
!arrays
 integer :: kptrlatt_current(3,3),kptrlatt_trial(3,3)
 integer,allocatable :: grid_list(:)
 real(dp) :: axes(3,3),gmet(3,3),gprimd(3,3),matrix1(3,3),matrix2(3,3)
 real(dp) :: metmin(3,3),minim(3,3),r2d(3,3),rmet(3,3),rsuper(3,3)
 real(dp) :: shiftk_current(3,210),shiftk_trial(3,210)
 real(dp),allocatable :: kpt(:,:),kptrlen_list(:),wtk(:)

! *************************************************************************

 kptrlen_target=kptrlen

!The vacuum array must be made of 0 or 1
 do ii=1,3
   if(vacuum(ii)/=0 .and. vacuum(ii)/=1)then
     write(message,'(a,a,a,i1,a,i3,a,a)')&
&     'The values of vacuum must be 0 or 1.',ch10,&
&     'However, the input vacuum(',ii,') is',vacuum(ii),ch10,&
&     'Action : correct vacuum in your input file.'
     MSG_ERROR(message)
   end if
 end do

!Specific preparation for 2-dimensional system
 if(sum(vacuum(:))==1)then

!  Make the non-active vector orthogonal to the active vectors,
!  and take it along the z direction
   if(vacuum(1)==1)then
     r2d(1,3)=rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3)
     r2d(2,3)=rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3)
     r2d(3,3)=rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3)
     r2d(:,1)=rprimd(:,2)
     r2d(:,2)=rprimd(:,3)
     dirvacuum=1
   else if(vacuum(2)==1)then
     r2d(1,3)=rprimd(2,3)*rprimd(3,1)-rprimd(3,3)*rprimd(2,1)
     r2d(2,3)=rprimd(3,3)*rprimd(1,1)-rprimd(1,3)*rprimd(3,1)
     r2d(3,3)=rprimd(1,3)*rprimd(2,1)-rprimd(2,3)*rprimd(1,1)
     r2d(:,1)=rprimd(:,3)
     r2d(:,2)=rprimd(:,1)
     dirvacuum=2
   else if(vacuum(3)==1)then
     r2d(1,3)=rprimd(2,1)*rprimd(3,2)-rprimd(3,1)*rprimd(2,2)
     r2d(2,3)=rprimd(3,1)*rprimd(1,2)-rprimd(1,1)*rprimd(3,2)
     r2d(3,3)=rprimd(1,1)*rprimd(2,2)-rprimd(2,1)*rprimd(1,2)
     r2d(:,1)=rprimd(:,1)
     r2d(:,2)=rprimd(:,2)
     dirvacuum=3
   end if
   surface=sqrt(sum(r2d(:,3)**2))
!  Identify the 2-D Bravais lattice
!  DEBUG
!  write(std_out,*)' r2d=',r2d(:,:)
!  ENDDEBUG
   call metric(gmet,gprimd,-1,rmet,r2d,ucvol)
   call smallprim(metmin,minim,r2d)
!  DEBUG
!  write(std_out,*)' minim=',minim(:,:)
!  ENDDEBUG
   ang90=0 ; equal=0 ; center=0
   axes(:,:)=minim(:,:)
   if(abs(metmin(1,2))<tol8)ang90=1
   if(abs(metmin(1,1)-metmin(2,2))<tol8)equal=1
   if(ang90==1)then
     if(equal==1)iholohedry=4
     if(equal==0)iholohedry=2
   else if(equal==1)then
     reduceda=metmin(1,2)/metmin(1,1)
     if(abs(reduceda+0.5_dp)<tol8)then
       iholohedry=3
     else if(abs(reduceda-0.5_dp)<tol8)then
       iholohedry=3
!      Use conventional axes
       axes(:,2)=minim(:,2)-minim(:,1)
     else
       iholohedry=2 ; center=1
       axes(:,1)=minim(:,1)+minim(:,2)
       axes(:,2)=minim(:,2)-minim(:,1)
     end if
   else
     reduceda=metmin(1,2)/metmin(1,1)
     reducedb=metmin(1,2)/metmin(2,2)
     if(abs(reduceda+0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,2)=2.0_dp*minim(:,2)+minim(:,1)
     else if(abs(reduceda-0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,2)=2.0_dp*minim(:,2)-minim(:,1)
     else if(abs(reducedb+0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,1)=2.0_dp*minim(:,1)+minim(:,2)
     else if(abs(reducedb-0.5_dp)<tol8)then
       iholohedry=2 ; center=1
       axes(:,1)=2.0_dp*minim(:,1)-minim(:,2)
     else
       iholohedry=1
     end if
   end if
!  Make sure that axes form a right-handed coordinate system
   determinant=axes(1,1)*axes(2,2)*axes(3,3) &
&   +axes(1,2)*axes(2,3)*axes(3,1) &
&   +axes(1,3)*axes(3,2)*axes(2,1) &
&   -axes(1,1)*axes(3,2)*axes(2,3) &
&   -axes(1,3)*axes(2,2)*axes(3,1) &
&   -axes(1,2)*axes(2,1)*axes(3,3)
   if(determinant<zero)then
     axes(:,1)=-axes(:,1)
   end if
!  Prefer symmetry axes on the same side as the primitive axes
   sca=axes(1,1)*r2d(1,1)+axes(2,1)*r2d(2,1)+axes(3,1)*r2d(3,1)
   scb=axes(1,2)*r2d(1,2)+axes(2,2)*r2d(2,2)+axes(3,2)*r2d(3,2)
   scc=axes(1,3)*rprimd(1,dirvacuum)&
&   +axes(2,3)*rprimd(2,dirvacuum)&
&   +axes(3,3)*rprimd(3,dirvacuum)
   if(sca<-tol8 .and. scb<-tol8)then
     axes(:,1)=-axes(:,1) ; sca=-sca
     axes(:,2)=-axes(:,2) ; scb=-scb
   end if
!  Doing this might change the angle between vectors, so that
!  the cell is not conventional anymore
!  if(sca<-tol8 .and. scc<-tol8)then
!  axes(:,1)=-axes(:,1) ; sca=-sca
!  axes(:,3)=-axes(:,3) ; scc=-scc
!  end if
!  if(scb<-tol8 .and. scc<-tol8)then
!  axes(:,2)=-axes(:,2) ; scb=-scb
!  axes(:,3)=-axes(:,3) ; scc=-scc
!  end if
   length_axis1=sqrt(axes(1,1)**2+axes(2,1)**2+axes(3,1)**2)
   length_axis2=sqrt(axes(1,2)**2+axes(2,2)**2+axes(3,2)**2)

!  DEBUG
!  write(std_out,*)' testkgrid : iholohedry, center =',iholohedry,center
!  write(std_out,*)' testkgrid : axis 1=',axes(:,1)
!  write(std_out,*)' testkgrid : axis 2=',axes(:,2)
!  write(std_out,*)' testkgrid : axis 3=',axes(:,3)
!  write(std_out,*)' testkgrid : length_axis=',length_axis1,length_axis2
!  ENDDEBUG

!  End special treatment of 2-D case
 end if

!3-dimensional system
 if(sum(vacuum(:))==0)then
   iholohedry=bravais(1)
   center=bravais(2)
   fact=1.0_dp
   if(center/=0)fact=0.5_dp
   matrix1(:,1)=bravais(3:5)*fact
   matrix1(:,2)=bravais(6:8)*fact
   matrix1(:,3)=bravais(9:11)*fact
   call matr3inv(matrix1,matrix2)
   do ii=1,3
     axes(:,ii)=rprimd(:,1)*matrix2(ii,1)+rprimd(:,2)*matrix2(ii,2)+rprimd(:,3)*matrix2(ii,3)
   end do
   length_axis1=sqrt(axes(1,1)**2+axes(2,1)**2+axes(3,1)**2)
   length_axis2=sqrt(axes(1,2)**2+axes(2,2)**2+axes(3,2)**2)
   length_axis3=sqrt(axes(1,3)**2+axes(2,3)**2+axes(3,3)**2)
!  DEBUG
!  write(std_out,*)' testkgrid : axes=',axes(:,:)
!  write(std_out,*)' length_axis=',length_axis1,length_axis2,length_axis3
!  ENDDEBUG
 end if

!This routine examine only primitive k lattices.
 nshiftk=1

!If prtkpt/=0, will examine more grids than strictly needed
 buffer_scale=one
 if(prtkpt/=0)buffer_scale=two

 if(prtkpt/=0)then
   write(message,'(a,a,a,a,a,a,a,a)' )ch10,&
&   ' testkgrid : will perform the analysis of a series of k-grids.',ch10,&
&   '  Note that kptopt=1 in this analysis, irrespective of its input value.',ch10,ch10,&
&   ' Grid#    kptrlatt         shiftk         kptrlen       nkpt  iset',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
   ABI_ALLOCATE(grid_list,(mkpt_list))
   ABI_ALLOCATE(kptrlen_list,(mkpt_list))
   grid_list(:)=0
   kptrlen_list(:)=0.0_dp
 end if

 if(sum(vacuum(:))==3)then

   kptrlatt(:,:)=0
   kptrlatt(1,1)=1
   kptrlatt(2,2)=1
   kptrlatt(3,3)=1
   shiftk(:,1)=0.0_dp
   kptrlen=1000.0_dp
   nkpt_current=1
   igrid_current=1

   if(prtkpt/=0)then
     write(message,&
&     '(a,3i4,a,es14.4,a,es14.4,i8,i6,a,a,3i4,a,es14.4,a,a,3i4,a,es14.4,a)' )&
&     '    1  ',kptrlatt(:,1),'  ',shiftk(1,1),'  ',kptrlen,1,1,ch10,&
&     '       ',kptrlatt(:,2),'  ',shiftk(2,1),ch10,&
&     '       ',kptrlatt(:,3),'  ',shiftk(3,1),ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
!    The unit cell volume is fake
     ucvol=kptrlen**3
   end if

 else

   nkpt=0 ; nkpt_current=0 ; iscf=1 ; iset=1
   kptrlen_current=0.0_dp
   mult1=0 ; mult2=0 ; mult3=0 ; init_mult=1
   ABI_ALLOCATE(kpt,(3,nkpt))
   ABI_ALLOCATE(wtk,(nkpt))
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!  Loop on different grids, the upper limit is only to avoid an infinite loop
   do igrid=1,1000

     kptrlatt_trial(:,:)=0
     kptrlatt_trial(1,1)=1
     kptrlatt_trial(2,2)=1
     kptrlatt_trial(3,3)=1
     shiftk_trial(:,1)=0.0_dp

!    1-dimensional system
     if(sum(vacuum(:))==2)then
       if(vacuum(1)==0)then
         kptrlatt_trial(1,1)=2*igrid ; shiftk_trial(1,1)=0.5_dp
       else if(vacuum(2)==0)then
         kptrlatt_trial(2,2)=2*igrid ; shiftk_trial(2,1)=0.5_dp
       else if(vacuum(3)==0)then
         kptrlatt_trial(3,3)=2*igrid ; shiftk_trial(3,1)=0.5_dp
       end if
     end if

!    2-dimensional system
     if(sum(vacuum(:))==1)then

!      Treat hexagonal holohedries separately
       if(iholohedry==3)then

!        write(std_out,*)' testkgrid : 2D, hexagonal'

         mult1=mult1+1
         nset=4
         if(iset==1)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           shiftk_trial(:,1)=0.0_dp
         else if(iset==2)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           shiftk_trial(1,1)=1.0_dp/3.0_dp
           shiftk_trial(2,1)=1.0_dp/3.0_dp
         else if(iset==3)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           shiftk_trial(:,1)=0.0_dp
         else if(iset==4)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           shiftk_trial(1,1)=0.5_dp
           shiftk_trial(2,1)=0.5_dp
         end if

       else
!        Now treat all other holohedries
         length1=length_axis1*mult1
         length2=length_axis2*mult2
!        DEBUG
!        write(std_out,*)' testkgrid : (2d) length=',length1,length2
!        ENDDEBUG
         if(abs(length1-length2)<tol8)then
           mult1=mult1+1
           mult2=mult2+1
         else if(length1>length2)then
           mult2=mult2+1
         else if(length2>length1)then
           mult1=mult1+1
         end if
         nset=4
!        iset==5 and 6 are allowed only for centered lattice
         if(center==1)nset=6
         if(iset==1 .or. iset==2)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult2
         else if(iset==3 .or. iset==4)then
           rsuper(:,1)=axes(:,1)*mult1-axes(:,2)*mult2
           rsuper(:,2)=axes(:,1)*mult1+axes(:,2)*mult2
         else if(iset==5 .or. iset==6)then
           rsuper(:,1)=axes(:,1)*(mult1-0.5_dp)-axes(:,2)*(mult2-0.5_dp)
           rsuper(:,2)=axes(:,1)*(mult1-0.5_dp)+axes(:,2)*(mult2-0.5_dp)
         end if
!        This was the easiest way to code all even mult1 and mult2 pairs :
!        make separate series for this possibility.
         if(iset==2 .or. iset==4 .or. iset==6)then
           rsuper(:,1)=2.0_dp*rsuper(:,1)
           rsuper(:,2)=2.0_dp*rsuper(:,2)
         end if
         shiftk_trial(1,1)=0.5_dp
         shiftk_trial(2,1)=0.5_dp

       end if

!      Put back the inactive direction
       if(dirvacuum==1)then
         rsuper(:,3)=rsuper(:,1)
         shiftk_trial(3,1)=shiftk_trial(1,1)
         rsuper(:,1)=rprimd(:,1)
         shiftk_trial(1,1)=0.0_dp
       else if(dirvacuum==2)then
         rsuper(:,3)=rsuper(:,1)
         shiftk_trial(3,1)=shiftk_trial(1,1)
         rsuper(:,1)=rsuper(:,2)
         shiftk_trial(1,1)=shiftk_trial(2,1)
         rsuper(:,2)=rprimd(:,2)
         shiftk_trial(2,1)=0.0_dp
       else if(dirvacuum==3)then
         rsuper(:,3)=rprimd(:,3)
         shiftk_trial(3,1)=0.0_dp
       end if

!      The supercell and the corresponding shift have been generated !
!      Convert cartesian coordinates into kptrlatt_trial
       do ii=1,3
         kptrlatt_trial(:,ii)=nint( gprimd(1,:)*rsuper(1,ii)+&
&         gprimd(2,:)*rsuper(2,ii)+&
&         gprimd(3,:)*rsuper(3,ii)  )
       end do

!      End of 2-dimensional system
     end if

!    3-dimensional system
     if(sum(vacuum(:))==0)then
!      Treat hexagonal holohedries separately
       if(iholohedry==6)then
         length1=length_axis1*mult1
         length3=length_axis3*mult3
!        DEBUG
!        write(std_out,*)' testkgrid : (hex) lengths=',length1,length2
!        ENDDEBUG
         if(abs(length1-length3)<tol8)then
           mult1=mult1+1
           mult3=mult3+1
         else if(length1>length3)then
           mult3=mult3+1
         else if(length3>length1)then
           mult1=mult1+1
         end if
         nset=4
         if(iset==1)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==2)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(1,1)=1.0_dp/3.0_dp
           shiftk_trial(2,1)=1.0_dp/3.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==3)then
           rsuper(:,1)=(axes(:,1)-axes(:,2))  *mult1
           rsuper(:,2)=(axes(:,1)+2*axes(:,2))*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.0_dp
           shiftk_trial(3,1)=0.5_dp
         else if(iset==4)then
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult1
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         end if

       else
!        Now treat all other holohedries
         length1=length_axis1*mult1
         length2=length_axis2*mult2
         length3=length_axis3*mult3
!        DEBUG
!        write(std_out,*)' testkgrid : length=',length1,length2,length3
!        ENDDEBUG
         if(length2>length1+tol8 .and. length3>length1+tol8)then
           mult1=mult1+1
         else if(length1>length2+tol8 .and. length3>length2+tol8)then
           mult2=mult2+1
         else if(length1>length3+tol8 .and. length2>length3+tol8)then
           mult3=mult3+1
         else if(abs(length2-length3)<tol8 .and. &
&           abs(length1-length3)<tol8 .and. &
&           abs(length1-length2)<tol8        )then
           mult1=mult1+1 ; mult2=mult2+1 ; mult3=mult3+1
         else if(abs(length1-length2)<tol8)then
           mult1=mult1+1 ; mult2=mult2+1
         else if(abs(length1-length3)<tol8)then
           mult1=mult1+1 ; mult3=mult3+1
         else if(abs(length2-length3)<tol8)then
           mult2=mult2+1 ; mult3=mult3+1
         end if
         nset=6
         if(center==-1 .or. center==-3)nset=8
         if(iset==1 .or. iset==2)then
!          Simple lattice of k points
           rsuper(:,1)=axes(:,1)*mult1
           rsuper(:,2)=axes(:,2)*mult2
           rsuper(:,3)=axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         else if(iset==3 .or. iset==4)then
!          FCC lattice of k points = BCC lattice in real space
           rsuper(:,1)=-axes(:,1)*mult1+axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,2)= axes(:,1)*mult1-axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,3)= axes(:,1)*mult1+axes(:,2)*mult2-axes(:,3)*mult3
           shiftk_trial(:,1)=0.5_dp
         else if(iset==5 .or. iset==6)then
!          BCC lattice of k points = FCC lattice in real space
           rsuper(:,1)=                 axes(:,2)*mult2+axes(:,3)*mult3
           rsuper(:,2)= axes(:,1)*mult1                +axes(:,3)*mult3
           rsuper(:,3)= axes(:,1)*mult1+axes(:,2)*mult2
!          The BCC lattice has no empty site with full symmetry
           shiftk_trial(:,1)=0.0_dp
         else if(iset==7 .or. iset==8)then
!          iset==7 and 8 are allowed only for centered lattice
           mult1h=mult1-0.5_dp
           mult2h=mult2-0.5_dp
           mult3h=mult3-0.5_dp
           if(center==-1)then
!            FCC lattice of k points = BCC lattice in real space
             rsuper(:,1)=-axes(:,1)*mult1h+axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,2)= axes(:,1)*mult1h-axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,3)= axes(:,1)*mult1h+axes(:,2)*mult2h-axes(:,3)*mult3h
             shiftk_trial(:,1)=0.5_dp
           else if(center==-3)then
!            BCC lattice of k points = FCC lattice in real space
             rsuper(:,1)=                  axes(:,2)*mult2h+axes(:,3)*mult3h
             rsuper(:,2)= axes(:,1)*mult1h                 +axes(:,3)*mult3h
             rsuper(:,3)= axes(:,1)*mult1h+axes(:,2)*mult2h
!            The BCC lattice has no empty site with full symmetry
             shiftk_trial(:,1)=0.0_dp
           end if
         end if
!        This was the easiest way to code all even mult1, mult2, mult3 triplets :
!        make separate series for this possibility.
         if(2*(iset/2)==iset)then
           rsuper(:,1)=2.0_dp*rsuper(:,1)
           rsuper(:,2)=2.0_dp*rsuper(:,2)
           rsuper(:,3)=2.0_dp*rsuper(:,3)
         end if
       end if

!      DEBUG
!      write(std_out,*)' testkgrid : gprimd=',gprimd(:,:)
!      write(std_out,*)' testkgrid : rsuper=',rsuper(:,:)
!      write(std_out,*)' testkgrid : iset  =',iset
!      ENDDEBUG


!      The supercell and the corresponding shift have been generated !
!      Convert cartesian coordinates into kptrlatt_trial
       do ii=1,3
         kptrlatt_trial(:,ii)=nint( gprimd(1,:)*rsuper(1,ii)+&
&         gprimd(2,:)*rsuper(2,ii)+&
&         gprimd(3,:)*rsuper(3,ii)  )
       end do

!      End of 3-dimensional system
     end if

!    DEBUG
!    write(std_out,*)' testkgrid : before getkgrid'
!    write(std_out,*)' testkgrid : rprimd=',rprimd(:,:)
!    write(std_out,*)' testkgrid : kptrlatt_trial=',kptrlatt_trial(:,:)
!    ENDDEBUG

     call getkgrid(0,0,iscf,kpt,&
&     kptopt,kptrlatt_trial,kptrlen_trial,&
&     msym,nkpt,nkpt_trial,nshiftk,nsym,rprimd,&
&     shiftk_trial,symafm,symrel,vacuum,wtk)

!    DEBUG
!    write(std_out,*)' testkgrid : after getkgrid'
!    ENDDEBUG

!    In case one does not need the full list of grids, will take a shortcut, and go to one of the last grids of the series,
!    that generates a kptrlen_trial that is just below kptrlen.
     if(prtkpt==0 .and. init_mult==1 .and. kptrlen_trial<(half-tol8)*kptrlen )then
       iscale=int((one-tol8)*kptrlen/kptrlen_trial)
       mult1=mult1*iscale
       mult2=mult2*iscale
       mult3=mult3*iscale
       init_mult=0
!       DEBUG
!       write(std_out,*)' testkgrid : iscale=',iscale
!       ENDDEBUG
       kptrlatt_trial(:,:)=kptrlatt_trial(:,:)*iscale
       call getkgrid(0,0,iscf,kpt,&
&       kptopt,kptrlatt_trial,kptrlen_trial,&
&       msym,nkpt,nkpt_trial,nshiftk,nsym,rprimd,&
&       shiftk_trial,symafm,symrel,vacuum,wtk)
     end if

     if( (kptrlen_trial+tol8>kptrlen*(1.0_dp+tol8) .and. nkpt_current==0) .or. &
&     (kptrlen_trial+tol8>kptrlen*(1.0_dp+tol8) .and. nkpt_trial<nkpt_current) .or. &
&     (nkpt_trial==nkpt_current  .and. kptrlen_trial>kptrlen_current*(1.0_dp+tol8)))then

       kptrlatt_current(:,:)=kptrlatt_trial(:,:)
       nkpt_current=nkpt_trial
       shiftk_current(:,:)=shiftk_trial(:,:)
       kptrlen_current=kptrlen_trial
       igrid_current=igrid
     end if

     if(prtkpt/=0)then
       write(message,'(i5,a,3i4,a,es14.4,a,es14.4,i8,i6,a,a,3i4,a,es14.4,a,a,3i4,a,es14.4,a)' )&
&       igrid,'  ',kptrlatt_trial(:,1),'  ',shiftk_trial(1,1),&
&       '  ',kptrlen_trial,nkpt_trial,iset,ch10,&
&       '       ',kptrlatt_trial(:,2),'  ',shiftk_trial(2,1),ch10,&
&       '       ',kptrlatt_trial(:,3),'  ',shiftk_trial(3,1),ch10
       call wrtout(std_out,message,'COLL')
       call wrtout(iout,message,'COLL')

!      Keep track of this grid, if it is worth
       if(kptrlen_trial > kptrlen_list(nkpt_trial)*(1.0_dp+tol8))then
         grid_list(nkpt_trial)=igrid
         kptrlen_list(nkpt_trial)=kptrlen_trial
       end if
     end if

!    Treat 1-D case
     if( sum(vacuum(:))==2 .and. kptrlen_trial>buffer_scale*(1.0_dp+tol8)*kptrlen )exit

!    Treat 2-D case or 3-D case
     if( sum(vacuum(:))<=1 .and. kptrlen_trial>buffer_scale*(1.0_dp+tol8)*kptrlen )then
!      The present set of sets of k points is finished :
!      either it was the last, or one has to go to the next one
       if(iset==nset)exit
       iset=iset+1
       mult1=0 ; mult2=0 ; mult3=0 ; init_mult=1
     end if

   end do ! igrid=1,1000

   ABI_DEALLOCATE(kpt)
   ABI_DEALLOCATE(wtk)

   kptrlatt(:,:)=kptrlatt_current(:,:)
   shiftk(:,:)=shiftk_current(:,:)
   kptrlen=kptrlen_current

 end if ! test on the number of dimensions

 if(prtkpt/=0)then

!  sqrt(1/2) comes from the FCC packing, the best one
   factor=sqrt(0.5_dp)/ucvol/dble(nsym)
   ndims=3
   if(sum(vacuum(:))/=0)then
     if(sum(vacuum(:))==1)then
!      sqrt(3/4) comes from the hex packing, the best one
!      one multiplies by 2 because nsym is likely twice the number
!      of symmetries that can be effectively used in 2D
       ndims=2 ; factor=sqrt(0.75_dp)/surface/dble(nsym)*2
       write(message,'(2a)' )ch10,' Note that the system is bi-dimensional.'
     else if(sum(vacuum(:))==2)then
       ndims=1 ; factor=1/ucvol
       write(message,'(2a)' )ch10,' Note that the system is uni-dimensional.'
     else if(sum(vacuum(:))==3)then
       ndims=0
       write(message,'(2a)' )ch10,' Note that the system is zero-dimensional.'
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')
   end if

!  The asymptotic value of the merit factor is determined
!  by the set of symmetries : in 3D, if it includes the
!  inversion symmetry, the limit will be 1, if not, it
!  will be two. In 2D, if it includes the inversion symmetry
!  and an operation that maps z on -z, it will tend to one,
!  while if only one of these operations is present,
!  it will tend to two, and if none is present, it will tend to four.
   write(message,'(11a)' )ch10,&
&   ' List of best grids, ordered by nkpt.',ch10,&
&   '  (stop at a value of kptrlen 20% larger than the target value).',ch10,&
&   '  (the merit factor will tend to one or two in 3 dimensions)',ch10,&
&   '  (and to one, two or four in 2 dimensions)',ch10,ch10,&
&   '    nkpt   kptrlen    grid#  merit_factor'
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

   kptrlen_max=0.0_dp
   do ii=1,mkpt_list
     if(kptrlen_list(ii)>kptrlen_max*(1.0_dp+tol8))then
       kptrlen_max=kptrlen_list(ii)
       merit_factor=kptrlen_max**ndims/dble(ii)*factor
       write(message, '(i6,es14.4,i6,f12.4)' )ii,kptrlen_max,grid_list(ii),merit_factor
       call wrtout(std_out,message,'COLL')
       call wrtout(iout,message,'COLL')
     end if
     if(kptrlen_max>1.2_dp*(1.0_dp-tol8)*kptrlen_target)exit
   end do

   write(message,'(a,a,es14.4,a,a,i6,a,a,a,es14.4,a,i6)' )ch10,&
&   ' For target kptrlen=',kptrlen_target,',',&
&   ' the selected grid is number',igrid_current,',',ch10,&
&   '     giving kptrlen=',kptrlen_current,' with nkpt=',nkpt_current
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')

   write(message,'(a,a,a,a)' )ch10,&
&   ' testkgrid : stop after analysis of a series of k-grids.',ch10,&
&   '  For usual production runs, set prtkpt back to 0 (the default).'
   call wrtout(std_out,message,'COLL',do_flush=.True.)
   call wrtout(iout,message,'COLL',do_flush=.True.)

   call xmpi_abort()
   call leave_new('PERS',exit_status=0,print_config=.false.)

 end if

end subroutine testkgrid
!!***
