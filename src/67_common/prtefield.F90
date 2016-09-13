!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtefield
!!
!! NAME
!! prtefield
!!
!! FUNCTION
!! Print components of electric field, displacement field and polarization in nice format
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, LBoeri, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt
!!   | efield 
!!   | dfield
!!   | red_efield 
!!   | red_efieldbar 
!!   | red_dfield 
!!  dtefield <type(efield_type)> 
!!   | efield2 
!!   | red_ptot1
!!  iunit = unit number to which the data is printed
!!  rprimd
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      gstate,update_e_field_vars
!!
!! CHILDREN
!!      metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtefield(dtset,dtefield,iunit,rprimd)

! use m_profiling_abi
 use defs_basis
 use defs_abitypes
 use m_efield

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtefield'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer :: iunit
 real(dp),intent(in) :: rprimd(3,3)
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(inout) :: dtset

!Local variables-------------------------------
! Do not modify the length of this string
!scalars
 integer :: idir,ii
 character(len=1500) :: message
 character(len=7)   :: flag_field(3)

 real(dp) ::    ucvol 
! arrays
 real(dp) :: ptot_cart(3),gmet(3,3),gprimd(3,3),rmet(3,3),red_pbar(3),red_dbar(3),red_dfieldbar(3)
 real(dp) :: red_efieldbar_lc(3),red_efield_lc(3) 


! *************************************************************************

!DEBUG
!write(iout, '(a)') ' prtefield : enter '
!ENDDEBUG

!write here

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)


 ptot_cart(:)=zero
 do idir = 1,3
   ptot_cart(idir)=rprimd(idir,1)*dtefield%red_ptot1(1) + rprimd(idir,2)*dtefield%red_ptot1(2) + &
&   rprimd(idir,3)*dtefield%red_ptot1(3)  
 end do
 ptot_cart(:)=ptot_cart(:)/ucvol

 if (dtset%berryopt == 4) then

!  to calculate e Eq.(25)  
   do idir=1,3
     dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
   end do

!  to calculate ebar Eq.(25)  
   do idir=1,3
     dtset%red_efieldbar(idir)  =dot_product(dtset%efield(:),rprimd(:,idir))
   end do


!  to calculate pbar
   do idir=1,3
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
   end do

!MGNAG: This message is to long and causes
! Runtime Error: wrtout_cpp.f90, line 893: Buffer overflow on output
! with NAG in test seq_tsv6_125 where we write to std_out!
! I cannot change the RECLEN of std_out!

   write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
&   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
&   '       P:  ', (ptot_cart(ii), ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
&   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   !!HONG need to change
&  '    pbar:  ', (red_pbar(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
&   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
&   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10,&
&   ' (S.I.), that is V/m for E, and C/m^2 for P', ch10, &  
&   '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3),ch10
   call wrtout(iunit,message,'COLL')
   
 end if ! berryopt ==4


 if (dtset%berryopt == 6) then 
   do idir=1,3
     dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
   end do

!  to calculate ebar   !! Need to be changed 
   do idir=1,3
     dtset%red_efieldbar(idir)  = dot_product(dtset%efield(:),rprimd(:,idir))
   end do

!  to calculate red_pbar
   do idir=1,3
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
   end do

!  to calculate red_dbar
   do idir=1,3
     red_dbar(idir)  = dot_product(dtset%dfield(:),rprimd(:,idir))
   end do

!  to calculate d 
   do idir=1,3
     dtset%red_dfield(idir)  =(ucvol/(4*pi))*dot_product(dtset%dfield(:),gprimd(:,idir))
   end do

   write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
&   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
&   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
&   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
&   ' e  +  p:  ', (dtset%red_efield(ii)+dtefield%red_ptot1(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10,&
&   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   !!HONG need to change
&  '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
&   '    dbar:  ', (red_dbar(ii),ii=1,3), ch10, &
&   ' eba+pba:  ', (dtset%red_efieldbar(ii)+red_pbar(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
&   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
&   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
&   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10,&
&   ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &  
&   '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10, &
&   '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3), ch10, &
&   'eps0*E+P:  ', (dtset%efield(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10  ! eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))=8.854187817620*5.14220652*1d-1
   call wrtout(iunit,message,'COLL')
   
   !MGNAG Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output

 end if  ! berryopt ==6


 if (dtset%berryopt == 14)  then

   do idir=1,3   ! ebar local
     red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir)) 
   end do


   do idir=1,3
     dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
   end do

!  to calculate pbar
   do idir=1,3
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
   end do

   write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))')  ' (a. u.)', ch10,&
&   '   ebar0:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   
&   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &
&   '    pbar:  ', (red_pbar(ii),ii=1,3) 
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
&   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       E:  ', (dtefield%efield2(ii), ii=1,3), ch10, &
&   '       P:  ', (ptot_cart(ii), ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10, &
&   ' (S.I.), that is V/m for E, and C/m^2 for P', ch10, &  
&   '-      E:  ', (dtefield%efield2(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3),ch10 
   call wrtout(iunit,message,'COLL')

   
 end if  ! berryopt ==14 


 if (dtset%berryopt == 16) then

!  to calculate e Eq.(25) 
   do idir=1,3
     dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtset%efield(:),gprimd(:,idir))
   end do

!  to calculate ebar
   do idir=1,3
     dtset%red_efieldbar(idir)  = dot_product(dtset%efield(:),rprimd(:,idir))
   end do

!  to calculate pbar
   do idir=1,3
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
   end do

!  to calculate dbar
   do idir=1,3
     red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
   end do

!  to calculate D
   do idir=1,3
     dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
   end do

   write(message,'(a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ' (a. u.)', ch10,&
&   '       e:  ', (dtset%red_efield(ii),ii=1,3), ch10, &
&   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
&   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
&   ' e  +  p:  ', (dtset%red_efield(ii)+dtefield%red_ptot1(ii),ii=1,3)  
    call wrtout(iunit,message,'COLL')

   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '    ebar:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   
&   '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
&   '    dbar:  ', (red_dfieldbar(ii),ii=1,3), ch10, &
&   ' eba+pba:  ', (dtset%red_efieldbar(ii)+red_pbar(ii),ii=1,3)
    call wrtout(iunit,message,'COLL')

   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       E:  ', (dtset%efield(ii), ii=1,3), ch10, &
&   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
&   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
&   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)') ch10, &
&   ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &  
&   '-      E:  ', (dtset%efield(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10, &
&   '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3), ch10, &
&   'eps0*E+P:  ', (dtset%efield(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10
   call wrtout(iunit,message,'COLL')

 end if  ! berryopt ==16

 if (dtset%berryopt == 17) then

   do idir=1,3   ! ebar local
     red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir)) 
   end do

   do idir=1,3
     dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
   end do

!  to calculate pbar
   do idir=1,3
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
   end do


!  do idir=1,3
!  if (dtset%rfdir(idir)==1) then   ! fixed ebar
!  red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))  ! local efieldbar
!  dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
!  red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
!  dtset%dfield(idir)=dtefield%efield2(idir)+4*pi*ptot_cart(idir)
!  dtset%red_dfield(idir)=dtset%red_efield+dtefield%red_ptot1(idir)  
!  dtset%red_dfieldbar(idir)=red_efieldbar_lc(idir)+red_pbar(idir)
!  E_lc(idir)=dtefield%efield2(idir)
!  e_lc(idir)=red_efieldbar_lc(idir)
!  ebar_lc(idir)=dtset%red_efieldbar(idir) 
!  else if (dtset%rfdir(idir)==2) then ! fixed d
!  dtset%red_efield(idir)  =(ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir))
!  dtset%red_efieldbar(idir)  = dot_product(dtefield%efield2(:),rprimd(:,idir))
!  red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
!  red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
!  dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))
!  E_lc(idir)=dtefield%efield2(idir)
!  e_lc(idir)=dtset%red_efield(idir)
!  ebar_lc(idir)=dtset%red_efieldbar(idir) 
!  end if
!  enddo

   

   do idir=1,3
     red_efield_lc(idir)= (ucvol/(4*pi))*dot_product(dtefield%efield2(:),gprimd(:,idir)) 
     red_efieldbar_lc(idir)=dot_product(dtefield%efield2(:),rprimd(:,idir))  ! local efieldbar
     red_pbar(idir)  = (4*pi/ucvol)*dot_product(dtefield%red_ptot1(:),rmet(:,idir))
     red_dfieldbar(idir)  = (4*pi/ucvol)*dot_product(dtset%red_dfield(:),rmet(:,idir))
     dtset%dfield(idir)  =(4*pi/ucvol)*dot_product(dtset%red_dfield(:),rprimd(:,idir))

   end do

   do idir=1,3
     if(dtset%jfielddir(idir)==1) then
       flag_field(idir)="E-field"
     else
       flag_field(idir)="D-field"
     end if
   end do

   write(message,'(a,a,a,6x,a,11x,a,11x,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') &
&   ' (a. u.)', ch10,&
&   '           ', (flag_field(ii),ii=1,3),ch10, & 
&   '   ebar0:  ', (dtset%red_efieldbar(ii),ii=1,3), ch10, &   
&   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &
&   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
&   ' e  +  p:  ', (red_efield_lc(ii)+dtefield%red_ptot1(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')

   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       e:  ', (red_efield_lc(ii),ii=1,3), ch10, &
&   '       p:  ', (dtefield%red_ptot1(ii), ii=1,3), ch10, &
&   '       d:  ', (dtset%red_dfield(ii),ii = 1, 3), ch10, &
&   ' e  +  p:  ', (red_efield_lc(ii)+dtefield%red_ptot1(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')
   
   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, & 
&   '    ebar:  ', (red_efieldbar_lc(ii),ii=1,3), ch10, &   
&   '    pbar:  ', (red_pbar(ii),ii=1,3), ch10, &
&   '    dbar:  ', (red_dfieldbar(ii),ii=1,3), ch10, &
&   ' eba+pba:  ', (red_efieldbar_lc(ii)+red_pbar(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')

   write(message,'(a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x))') ch10, &
&   '       E:  ', (dtefield%efield2(ii), ii=1,3), ch10, &
&   '       P:  ', (ptot_cart(ii), ii=1,3), ch10, &
&   '       D:  ', (dtset%dfield(ii),ii = 1, 3), ch10, &
&   'E+4*pi*P:  ', (dtset%efield(ii)+4.0d0*pi*ptot_cart(ii),ii=1,3)
   call wrtout(iunit,message,'COLL')

   write(message,'(a,a,a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a,a,3(es16.9,2x),a)')  ch10, &
&   ' (S.I.), that is V/m for E, and C/m^2 for P and D', ch10, &  
&   '       E:  ', (dtefield%efield2(ii)*(Ha_J/(e_Cb*Bohr_Ang*1d-10)), ii=1,3), ch10, &    !(Ha_J/(e_Cb*Bohr_Ang*1d-10))= 5.14220652*1d+11
&  '       P:  ', (ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2, ii=1,3), ch10, &
&   '       D:  ', ((1.0d0/(4*pi))*dtset%dfield(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii = 1, 3), ch10, &
&   'eps0*E+P:  ', (dtefield%efield2(ii)*eps0*(Ha_J/(e_Cb*Bohr_Ang*1d-10))+ptot_cart(ii)*(e_Cb)/(Bohr_Ang*1d-10)**2,ii=1,3),ch10
   call wrtout(iunit,message,'COLL')

 end if  ! berryopt ==17



end subroutine prtefield
!!***
