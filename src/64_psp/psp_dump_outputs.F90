!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp_dump_outputs
!! NAME
!! psp_dump_outputs
!!
!! FUNCTION
!! (To be described ...)
!!
!! COPYRIGHT
!! Copyright (C) 2017-2018 ABINIT group (YP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!! (to be filled)
!!
!! OUTPUT
!! (to be filled)
!!
!! SIDE EFFECTS
!! (to be filled)
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!
!! SOURCE

!      call psp_dump_outputs(psps%lmnmax,psps%mpssoang,psps%lnmax, &
!&      psps%mqgrid_ff,psps%n1xccc,mmax,epsatm,qchrg,xcccrc,nctab, &
!&      indlmn,nproj,ekb,ffspl,vlspl,xccc1d)

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp_dump_outputs(pfx,pspcod,lmnmax,lnmax,mpssoang, &
&      mqgrid,n1xccc,mmax,maxrad,epsatm,qchrg,xcccrc,nctab, &
&      indlmn,nproj,ekb,ffspl,vlspl,xccc1d)

 use defs_basis
 use defs_datatypes, only : nctab_t
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp_dump_outputs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*), intent(in) :: pfx
 integer,intent(in) :: pspcod,lmnmax,lnmax,mpssoang,mqgrid,n1xccc
 integer,intent(in) :: mmax
 real(dp),intent(in) :: maxrad,epsatm,qchrg,xcccrc
 type(nctab_t),intent(in) :: nctab
!arrays
 integer,intent(in) :: indlmn(6,lmnmax),nproj(mpssoang)
 real(dp),intent(in) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(in) :: xccc1d(n1xccc,6)

!Local variables ------------------------------
!scalars
 integer, parameter :: dump = 64
 integer :: ierr, i, j ,k
 character(len=500) :: msg

 ! *********************************************************************

 open(unit=dump, file=trim(pfx)//"_psp_info.yaml", status='REPLACE', err=10, iostat=ierr)

 write(dump,'(3a)') "%YAML 1.2", ch10, "---"

 write(dump, '(2a)') ch10, "# Pseudopotential info"
 write(dump, '(a,1x,i8)') "pspcod:", pspcod

 write(dump, '(2a)') ch10, "# Array dimensions"
 write(dump, '(a)') "dims:"
 write(dump, '(4x,a,1x,i8)') "lmnmax:", lmnmax
 write(dump, '(4x,a,1x,i8)') "lnmax:", lnmax
 write(dump, '(4x,a,1x,i8)') "mpssoang:", mpssoang
 write(dump, '(4x,a,1x,i8)') "mqgrid:", mqgrid
 write(dump, '(4x,a,1x,i8)') "n1xccc:", n1xccc
 write(dump, '(4x,a,1x,i8)') "mmax:", mmax

 write(dump, '(2a)') ch10, "# Quantities"
 write(dump, '(a,1x,e12.5)') "maxrad:", maxrad
 write(dump, '(a,1x,e12.5)') "epsatm:", epsatm
 write(dump, '(a,1x,e12.5)') "qchrg:", qchrg
 write(dump, '(a,1x,e12.5)') "xcccrc:", xcccrc

 write(dump, '(2a)') ch10, "# Structure: nctab"
 write(dump, '(a)') "nctab:"
 write(dump,'(4x,a,":",1x,i4)') "mqgrid_vl", nctab%mqgrid_vl
 write(dump,'(4x,a,":",1x,l4)') "has_tvale", nctab%has_tvale
 write(dump,'(4x,a,":",1x,l4)') "has_tcore", nctab%has_tcore
 write(dump,'(4x,a,":",1x,e12.5)') "dncdq0", nctab%dncdq0
 write(dump,'(4x,a,":",1x,e12.5)') "d2ncdq0", nctab%d2ncdq0
 write(dump,'(4x,a,":",1x,e12.5)') "dnvdq0", nctab%dnvdq0

 if ( nctab%has_tvale ) then
   write(dump, '(2a)') ch10, "# Array: nctab_tvalespl(mqgrid_vl,2)"
   write(dump, '(a)') "nctab_tvalespl:"
   do j=1,2
     do i=1,nctab%mqgrid_vl
       if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "- -", nctab%tvalespl(i,j)
       else
         write(dump,'(4x,a,1x,e12.5)') "  -", nctab%tvalespl(i,j)
       end if
     end do
   end do
 end if

 if ( nctab%has_tcore ) then
   write(dump, '(2a)') ch10, "# Array: nctab_tcorespl(mqgrid_vl,2)"
   write(dump, '(a)') "nctab_tcorespl:"
   do j=1,2
     do i=1,nctab%mqgrid_vl
       if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "- -", nctab%tcorespl(i,j)
       else
         write(dump,'(4x,a,1x,e12.5)') "  -", nctab%tcorespl(i,j)
       end if
     end do
   end do
 end if

 write(dump, '(2a)') ch10, "# Array: integer indlmn(6,lmnmax)"
 write(dump, '(a)') "indlmn:"
 do i=1,lmnmax
   write(dump,'(4x,a,i4,5(",",i4),a)') "- [", indlmn(:,i), "]"
 end do

 write(dump, '(2a)') ch10, "# Array: integer nproj(mpssoang)"
 write(dump, '(a)') "nproj:"
 do i=1,mpssoang
   write(dump,'(4x,"-",1x,i4)') nproj(i)
 end do

 write(dump, '(2a)') ch10, "# Array: double ekb(lnmax)"
 write(dump, '(a)') "ekb:"
 do i=1,lnmax
   write(dump,'(4x,"-",1x,e12.5)') ekb(i)
 end do

 write(dump, '(2a)') ch10, "# Array: ffspl(mqgrid,2,lnmax)"
 write(dump, '(a)') "ffspl:"
 do k=1,lnmax
   do j=1,2
     do i=1,mqgrid
       if ( (i == 1) .and. (j == 1) ) then
         write(dump,'(4x,a,1x,e12.5)') "- - -", ffspl(i,j,k)
       else if ( i == 1 ) then
         write(dump,'(4x,a,1x,e12.5)') "  - -", ffspl(i,j,k)
       else
         write(dump,'(4x,a,1x,e12.5)') "    -", ffspl(i,j,k)
       end if
     end do
   end do
 end do

 write(dump, '(2a)') ch10, "# Array: vlspl(mqgrid,2)"
 write(dump, '(a)') "vlspl:"
 do j=1,2
   do i=1,mqgrid
     if ( i == 1 ) then
       write(dump,'(4x,a,1x,e12.5)') "- -", vlspl(i,j)
     else
       write(dump,'(4x,a,1x,e12.5)') "  -", vlspl(i,j)
     end if
   end do
 end do

 write(dump, '(2a)') ch10, "# Array: xccc1d(n1xccc,6)"
 write(dump, '(a)') "xccc1d:"
 do j=1,6
   do i=1,mqgrid
     if ( i == 1 ) then
       write(dump,'(4x,a,1x,e12.5)') "- -", xccc1d(i,j)
     else
       write(dump,'(4x,a,1x,e12.5)') "  -", xccc1d(i,j)
     end if
   end do
 end do

 write (dump,'(2a)') ch10, "..."

 close(dump)

 return

 10 continue

 if ( ierr /= 0 ) then
   write(msg,'(a,a,a,i8)') "Error writing pseudopotential information", &
&   ch10, "IOSTAT=", ierr
   MSG_WARNING(msg)
 end if

end subroutine psp_dump_outputs
!!***
