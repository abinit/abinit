!{\src2tex{textfont=tt}}
!!****f* ABINIT/adini
!! NAME
!! adini
!!
!! FUNCTION
!! Analysis of the input string "inpstr" (the content of input file)
!! and setting of the corresponding input variables
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  inpstr=character string containing the input data, to be treated
!!  lenstr=actual length of the string contained in inpstr
!!
!! OUTPUT
!!  aim_dtset=the structured entity containing all input variables
!!
!! WARNING
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!      consist,inread
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine adini(aim_dtset,inpstr,lenstr)

 use defs_basis
 use defs_abitypes
 use defs_aimprom
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'adini'
 use interfaces_42_parser
 use interfaces_63_bader, except_this_one => adini
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr
 character(len=*),intent(in) :: inpstr
!no_abirules
 type(aim_dataset_type), intent(inout) :: aim_dtset !vz_i

!Local variables ------------------------------
!scalars
 integer :: errcod,ii,inxh,ipos,jj,lenc,ll,outi,tstngr=0,tstvpt=0 !vz_z
 real(dp) :: outr
 logical :: nbtst,try
 character(len=20) :: cmot

! *********************************************************************

 if (iachar(inpstr(1:1)) < 32) then
   ipos=2
 else
   ipos=1
 end if

 write(std_out,*) 'ECHO of the INPUT'
 write(std_out,*) '************************'
 write(untout,*) 'ECHO of the INPUT'
 write(untout,*) '************************'

 mread:  do ii=1,lenstr
   try=.false.
   nbtst=.true.
   inxh=index(inpstr(ipos:lenstr),' ')
   if ((ipos >= lenstr)) exit
   if ((inxh==2).or.(inxh==1)) then
     ipos=ipos+inxh
     cycle
   end if
   lenc=inxh-1
   cmot(1:lenc)=inpstr(ipos:ipos+inxh-2)
   ipos=ipos+inxh
!  write(std_out,*) cmot(1:lenc), lenc

   select case (cmot(1:lenc))

!    DRIVER SPECIFICATIONS

   case ('SURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if ((inxh /= 2).and.(inpstr(ipos:ipos)/='-')) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%isurf=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%isurf
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%isurf
     ipos=ipos+inxh

   case ('CRIT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if ((inxh /= 2).and.(inpstr(ipos:ipos)/='-')) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%crit=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%crit
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%crit
     ipos=ipos+inxh

   case ('RSURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%irsur=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%irsur
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%irsur
     ipos=ipos+inxh

   case ('FOLLOW')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%foll=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%foll
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%foll
     ipos=ipos+inxh

   case ('IRHO')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%irho=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%irho
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%irho
     ipos=ipos+inxh

   case ('PLDEN')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%plden=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%plden
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%plden
     ipos=ipos+inxh


   case ('IVOL')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%ivol=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%ivol
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%ivol
     ipos=ipos+inxh

   case ('DENOUT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%denout=outi
     if ((aim_dtset%denout < -1).or.(aim_dtset%denout>3)) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%denout
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%denout
     ipos=ipos+inxh

   case ('LAPOUT')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%lapout=outi
     if ((aim_dtset%lapout < -1).or.(aim_dtset%lapout>3)) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%lapout
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%lapout
     ipos=ipos+inxh

   case ('DLTYP')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%dltyp=outi
     write(std_out,*) cmot(1:lenc),'     ', aim_dtset%dltyp
     write(untout,*) cmot(1:lenc),'     ', aim_dtset%dltyp
     ipos=ipos+inxh

   case ('GPSURF')
     inxh=index(inpstr(ipos:lenstr),' ')
     if (inxh /= 2) then
       write(std_out,*) 'ERROR in specif. of ', cmot(1:lenc)
       MSG_ERROR("Aborting now")
     end if
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%gpsurf=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%gpsurf
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%gpsurf
     ipos=ipos+inxh


!      END OF THE DRIVER SPECIFICATIONS

   case ('ATOM')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%batom=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%batom
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%batom
     ipos=ipos+inxh

   case ('NSA')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsa=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsa
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsa
     ipos=ipos+inxh

   case ('NSB')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsb=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsb
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsb
     ipos=ipos+inxh

   case ('NSC')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nsc=outi
     write(std_out,*) cmot(1:lenc),'       ', aim_dtset%nsc
     write(untout,*) cmot(1:lenc),'       ', aim_dtset%nsc
     ipos=ipos+inxh

   case ('INPT')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%npt=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%npt
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%npt
     ipos=ipos+inxh

   case ('NTHETA')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nth=outi
     write(std_out,*) cmot(1:lenc),'    ', aim_dtset%nth
     write(untout,*) cmot(1:lenc),'    ', aim_dtset%nth
     ipos=ipos+inxh

   case ('NPHI')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
     aim_dtset%nph=outi
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%nph
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%nph
     ipos=ipos+inxh

   case ('THETAMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%themin=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themin
     ipos=ipos+inxh

   case ('THETAMAX')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%themax=outr
     write(std_out, '(1x,a,a,es17.10)' ) cmot(1:lenc),'  ', aim_dtset%themax
     write(untout,'(1x,a,a,es17.10)') cmot(1:lenc),'  ', aim_dtset%themax
     ipos=ipos+inxh

   case ('PHIMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phimin=outr
     write(std_out, '(1x,a,a,es17.10)' ) cmot(1:lenc),'    ', aim_dtset%phimin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimin
     ipos=ipos+inxh

   case ('PHIMAX')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phimax=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimax
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%phimax
     ipos=ipos+inxh

   case ('ATRAD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%atrad=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'     ', aim_dtset%atrad
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'     ', aim_dtset%atrad
     ipos=ipos+inxh

   case ('RADSTP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%dr0=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dr0
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dr0
     ipos=ipos+inxh

   case ('FOLSTP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%folstp=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%folstp
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%folstp
     ipos=ipos+inxh

   case ('RATMIN')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%rmin=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%rmin
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%rmin
     ipos=ipos+inxh

   case ('COFF1')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%coff1=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff1
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff1
     ipos=ipos+inxh

   case ('COFF2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%coff2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%coff2
     ipos=ipos+inxh

   case ('DPCLIM')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%dpclim=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dpclim
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%dpclim
     ipos=ipos+inxh

   case ('LGRAD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lgrad=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad
     ipos=ipos+inxh

   case ('LGRAD2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lgrad2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lgrad2
     ipos=ipos+inxh

   case ('LSTEP')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lstep=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep
     ipos=ipos+inxh

   case ('LSTEP2')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%lstep2=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep2
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%lstep2
     ipos=ipos+inxh

   case ('RSURDIR')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%th0=outr
     ipos=ipos+inxh
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%phi0=outr
     ipos=ipos+inxh
     write(std_out, '(1x,a,a,2es17.10)') cmot(1:lenc),'   ', aim_dtset%th0, aim_dtset%phi0
     write(untout, '(1x,a,a,2es17.10)') cmot(1:lenc),'   ', aim_dtset%th0, aim_dtset%phi0

   case ('FOLDEP')
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
       aim_dtset%foldep(jj)=outr
       ipos=ipos+inxh
     end do
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%foldep
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%foldep

   case ('SCAL')
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
       aim_dtset%scal(jj)=outr
       ipos=ipos+inxh
     end do
     write(std_out,*) cmot(1:lenc),'      ', aim_dtset%scal
     write(untout,*) cmot(1:lenc),'      ', aim_dtset%scal

   case ('NGRID')
     try=.true.
     do jj=1,3
       inxh=index(inpstr(ipos:lenstr),' ')
       call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"INT",outi,outr,errcod)
       aim_dtset%ngrid(jj)=outi
       if (.not.nbtst) then
         tstngr=jj-1
         cycle mread
       end if
       if (inxh==0) then
         tstvpt=jj
         exit mread
       end if
       ipos=ipos+inxh
       if (ipos==lenstr-1) then
         tstngr=jj
         exit mread
       end if
     end do
!      Why no echo ?? XG 030218
     tstngr=3

   case ('VPTS')
     do jj=1,4
       do ll=1,3
         try=.true.
         inxh=index(inpstr(ipos:lenstr),' ')
         call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
         aim_dtset%vpts(ll,jj)=outr
         if (.not.nbtst) then
           tstvpt=jj-1
           cycle mread
         end if
         ipos=ipos+inxh
         if (ipos>=lenstr) then
           tstvpt=jj
           exit mread
         end if
       end do
     end do
!      Why no echo ?? XG 030218
     tstvpt=4

   case ('MAXATD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%maxatd=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxatd
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxatd
     ipos=ipos+inxh

   case ('MAXCPD')
     inxh=index(inpstr(ipos:lenstr),' ')
     call inread(inpstr(ipos:ipos+inxh-2),inxh-1,"DPR",outi,outr,errcod)
     aim_dtset%maxcpd=outr
     write(std_out, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxcpd
     write(untout, '(1x,a,a,es17.10)') cmot(1:lenc),'    ', aim_dtset%maxcpd
     ipos=ipos+inxh

   case default
     write(std_out,*) 'ERROR Bad key word ! ',cmot(1:lenc)
   end select
 end do mread

 write(std_out,*) '************************'

 call consist(aim_dtset,tstngr,tstvpt)

end subroutine adini
!!***
