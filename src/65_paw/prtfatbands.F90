!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtfatbands
!! NAME
!! prtfatbands
!!
!! FUNCTION
!! Print dos_fractions_m in order to plot easily fatbands
!! if pawfatbnd=1  1 : fatbands are resolved in L.
!! if pawfatbnd=1  2 : fatbands are resolved in L and M.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
!!               = m-resolved projected dos inside PAW sphere.
!!  dtset
!!  pawfatbnd    = keyword for fatbands
!!  fermie       = Fermi energy
!!  eigen        = eigenvalues
!!  mbesslang    =maximum angular momentum for Bessel function expansion
!!  m_dos_flag   =option for the m-contributions to the partial DOS
!!  ndosfraction =natsph*mbesslang
!!
!! OUTPUT
!! (only writing)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,int2char4,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prtfatbands(dos_fractions_m,dtset,fildata,fermie,eigen,&
&  mbesslang,m_dos_flag,ndosfraction,pawfatbnd,pawtab)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_atomdata

 use m_fstrings,  only : int2char4
 use m_io_tools,  only : open_file
 use m_pawtab,    only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtfatbands'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m_dos_flag,mbesslang,ndosfraction,pawfatbnd
 type(dataset_type),intent(in) :: dtset
 real(dp),intent(in) :: fermie
 character(len=fnlen),intent(in) :: fildata
!arrays
 real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iall,il,iat,natsph,inbfatbands,iband,mband,ixfat,isppol,nkpt,lmax,ll,mm
 integer :: band_index,ikpt,nband_k
 real(dp) :: xfatband
 character(len=1) :: tag_l,tag_1m,tag_is
 character(len=2) :: tag_2m
 character(len=10) :: tag_il,tag_at,tag_grace
 character(len=500) :: message
 character(len=fnlen) :: tmpfil
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: unitfatbands_arr(:,:)
 real(dp),allocatable :: eigenvalues(:,:,:)
 character(len=2) :: symbol

!*************************************************************************

 DBG_ENTER("COLL")
 
 
 if(m_dos_flag.ne.0) then
   write(message,'(3a)')&
&   'm decomposed dos is activated',ch10, &
&   'Action: deactivate it with prtdosm=0 !'
   MSG_ERROR(message)
 end if

 if(dtset%nspinor==2) then
   message = "Fatbands are not yet available in the case nspinor==2 !"
   MSG_WARNING(message)
 end if

 natsph=dtset%natsph
 nkpt=dtset%nkpt
 mband=dtset%mband

 if(natsph>1000) then
   write(message,'(3a)')&
&   ' Too big number of fat bands !',ch10, &
&   ' Action: decrease natsph in input file !'
   MSG_ERROR(message)
 end if

!--------------  PRINTING IN LOG
 write(message,'(a,a,a,a,i5,a,a,1000i5)') ch10," ***** Print of fatbands activated ****** ",ch10,&
& "  Number of atom: natsph = ",natsph,ch10, &
& "  atoms  are             = ",(dtset%iatsph(iat),iat=1,natsph)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 iall=0;inbfatbands=0

 if(pawfatbnd==1) then
   inbfatbands=mbesslang-1
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L)",ch10
 else if(pawfatbnd==2) then
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L and M)",ch10
   inbfatbands=(mbesslang-1)**2
 end if
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(a,e12.5,a,e12.5,a)') "  Fermi energy is ",fermie*Ha_eV," eV = ",fermie," Ha"
 call wrtout(std_out,message,'COLL')

!--------------  OPEN AND NAME FILES FOR FATBANDS
 ABI_ALLOCATE(unitfatbands_arr,(natsph*inbfatbands,dtset%nsppol))
 unitfatbands_arr = -3

 do iat=1,natsph
   lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
   call int2char4(dtset%iatsph(iat),tag_at)
   ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
   call atomdata_from_znucl(atom,dtset%znucl(dtset%typat(dtset%iatsph(iat))))
   symbol = atom%symbol
   do il=1,inbfatbands
     iall=iall+1
     ll=int(sqrt(float(il-1)))  ! compute l
     if(ll.le.lmax) then  ! print only angular momentum included in the PAW data
       do isppol=1,dtset%nsppol
         write(tag_is,'(i1)')isppol
         if(pawfatbnd==1) then
           call int2char4(il-1,tag_il)
           ABI_CHECK((tag_il(1:1)/='#'),'Bug: string length too short!')
           tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//trim(tag_il)
         else if (pawfatbnd==2) then
           write(tag_l,'(i1)') ll
           mm=il-(ll**2+ll+1)      ! compute m
           if(mm<0) write(tag_2m,'(i2)') mm
           if(mm>=0) write(tag_1m,'(i1)') mm
           if(mm<0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m'//tag_2m
           if(mm>=0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m+'//tag_1m
         end if
         !unitfatbands_arr(iall,isppol)=tmp_unit+100+iall-1+(natsph*inbfatbands)*(isppol-1)
         !open (unit=unitfatbands_arr(iall,isppol),file=trim(tmpfil),status='unknown',form='formatted')
         if (open_file(tmpfil, message, newunit=unitfatbands_arr(iall,isppol), status='unknown',form='formatted') /= 0) then
           MSG_ERROR(message)
         end if

         write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitfatbands_arr(iall,isppol)
         call wrtout(std_out,message,'COLL')
         write(message,'(9a)') "# ",ch10,"# ABINIT package : FATBAND file ", ch10,&
&         "# It contains, for each band: the eigenvalues in eV (and the character of the band) as a function of the k-point",&
&         ch10,"# This file can be read with xmgrace (http://plasma-gate.weizmann.ac.il/Grace/)  ",ch10,"#  "
         call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
         do iband=1,mband
           call int2char4(iband-1,tag_grace)
           ABI_CHECK((tag_grace(1:1)/='#'),'Bug: string length too short!')
           write(message,'(16a)') ch10,"@    s",trim(tag_grace)," line color 1",&
&           ch10,"@    s",trim(tag_grace)," errorbar color 2",&
&           ch10,"@    s",trim(tag_grace)," errorbar riser linewidth 5.0", &
&           ch10,"@    s",trim(tag_grace)," errorbar linestyle 0"
           call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
         end do  !iband
         write(message,'(a,a)') ch10,'@type xydy'
         call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
       end do   ! isppol
     end if ! ll=<lmax
   end do   ! il
 end do  ! iat

 if(iall.ne.(natsph*inbfatbands)) then
   MSG_ERROR("error1 ")
 end if


 
!--------------  WRITE FATBANDS IN FILES
 if (pawfatbnd>0) then
   ABI_ALLOCATE(eigenvalues,(nkpt,mband,dtset%nsppol))
   band_index=0.d0
   do isppol=1,dtset%nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       do iband=1,mband
         eigenvalues(ikpt,iband,isppol)= eigen(iband+band_index)-fermie
       end do
       band_index=band_index+nband_k
     end do
   end do
   iall=0
   do iat=1,natsph
     lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
     do il=1,inbfatbands
       iall=iall+1
       ll=int(sqrt(float(il-1)))
       if(ll.le.lmax) then
         do isppol=1,dtset%nsppol
           do iband=1,mband
             write(message,'(a,a,i8)') ch10,"# BAND number :",iband
             call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
             do ikpt=1,nkpt
               if(pawfatbnd==1) then
                 xfatband=0.d0
                 do ixfat=(il-1)**2+1,il**2
                   xfatband=xfatband+dos_fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+ixfat)
                 end do ! ixfat
               else if (pawfatbnd==2) then
                 xfatband=dos_fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+il)
               end if
               write(message,'(i5,e20.5,e20.5)') ikpt-1,eigenvalues(ikpt,iband,isppol)*Ha_eV,xfatband
               call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
             end do ! ikpt
           end do  !iband
           write(message,'(a)') '&'
           call wrtout(unitfatbands_arr(iall,isppol),message,'COLL')
           !close(unitfatbands_arr(iall,isppol))
         end do  !isppol
       end if
     end do ! il
   end do ! iat
   ABI_DEALLOCATE(eigenvalues)
 end if

 do isppol=1,size(unitfatbands_arr, dim=2)
   do iat=1,size(unitfatbands_arr, dim=1)
     if (unitfatbands_arr(iat, isppol) /= -3) close (unitfatbands_arr(iat, isppol))
   end do 
 end do

 ABI_DEALLOCATE(unitfatbands_arr)

 DBG_EXIT("COLL")

end subroutine prtfatbands

!!***
