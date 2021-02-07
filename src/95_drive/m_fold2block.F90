!!****m* ABINIT/m_fold2block
!! NAME
!!  m_fold2block
!!
!! FUNCTION
!!  This module contains basic tools to operate on vectors expressed in reduced coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (AB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_fold2block

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 public
 !private

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_fold2block/SortC
!! NAME
!! SortC
!!
!! FUNCTION
!! Part of fold2Bloch that unfolds the wavefunction and calculates
!! the weight of each band.
!!
!! INPUTS
!! FX, FY, FZ= Number of folds in each dimention
!! vector= Array of energy coefficient position vectors
!! coefc= Array of energy coefficients of in a certain band
!! NV= Number of vectors/coefficients
!!
!! OUTPUT
!! wegihts= Calculated weights of a band in an unfolded state.
!! Depends on the number of folds the WFK file was structured with.
!!
!! PARENTS
!!      fold2Bloch
!!
!! CHILDREN
!!
!! SOURCE


 subroutine SortC(FX, FY, FZ, Vector, CoefC, NV, Weights)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: FX, FY, FZ, NV
!arrays
 integer, intent(in) :: Vector(3,NV)
 real(dp), intent(in) :: CoefC(2,NV)
 real(dp), intent(inout) ::  Weights(FX*FY*FZ)

!Local variables-------------------------------
!scalars
 real :: sumtot
 integer :: remainder_x, remainder_y, remainder_z, jj, kk, ll, pp, el
!arrays
 real(dp), allocatable :: Sums(:)
 real(dp), allocatable :: coefsqr(:),TGroupC(:,:,:,:)
 integer, allocatable :: counter(:,:,:)

! *************************************************************************

 ABI_MALLOC(TGroupC,(FX, FY, FZ, NV))
 ABI_MALLOC(Sums,(FX*FY*FZ))
 ABI_MALLOC(counter,(FX,FY,FZ))
 ABI_MALLOC(coefsqr,(NV))

 !Convert to sum of squares
 do jj=1, NV
   coefsqr(jj)=coefc(1,jj)*coefc(1,jj)+coefc(2,jj)*coefc(2,jj)
 end do

 !Initiates the counter and TGroupC elements at 0
 do jj=1,FX
   do kk=1,FY
     do ll=1,FZ
       counter(jj,kk,ll)=0
       do pp=1,NV
         TGroupC(jj,kk,ll,pp)=0.0
       end do
     end do
   end do
 end do

 !Sorts the energy coeeficients
 do jj=1, (NV)
   remainder_x=MODULO(Vector(1,jj), FX)
   remainder_y=MODULO(Vector(2,jj), FY)
   remainder_z=MODULO(Vector(3,jj), FZ)
   counter(remainder_x+1, remainder_y+1, remainder_z+1)=counter(remainder_x+1, remainder_y+1, remainder_z+1)+1
   TGroupC(remainder_x+1, remainder_y+1, remainder_z+1, counter(remainder_x+1, remainder_y+1, remainder_z+1))=Coefsqr(jj)
 end do

 !Sums the squares  of all coefficients per group
 el=1
 do jj=1, FX
   do kk=1, FY
     do ll=1, FZ
       if (counter(jj, kk, ll)>0) then
         Sums(el)=SUM(TGroupC(jj, kk, ll,1:counter(jj, kk, ll)))
         el=el+1
       else
         Sums(el)=0.0
         el=el+1
       end if
     end do
   end do
 end do

 !Assign weights to an array
 sumtot=SUM(Sums)
 do jj=1, (FX*FY*FZ)
   Weights(jj)=Sums(jj)/sumtot
 end do
 ABI_FREE(TGroupC)
 ABI_FREE(Sums)
 ABI_FREE(counter)
 ABI_FREE(coefsqr)

 end subroutine SortC
!!***

!!****f* m_fold2block/NewK
!! NAME
!! NewK
!!
!! FUNCTION
!! Part of fold2Bloch that determines the unfolded
!! positions of each K point.
!!
!! INPUTS
!! FX, FY, FZ= Number of folds in each dimention
!! XX, YY, ZZ= Original K point coordinates
!!
!! OUTPUT
!! nkval= Array of determined coordinates of unfolded K point locaations.
!! Depends on the number of folds the WFK file was structured with.
!!
!! PARENTS
!!      fold2Bloch
!!
!! CHILDREN
!!
!! SOURCE

 subroutine NewK(XX, YY, ZZ, FX, FY, FZ, NKVal)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: FX, FY, FZ
 real(dp), intent(in) :: XX, YY, ZZ
 real(dp), intent(inout) :: NKVal (3,FX*FY*FZ)
!arrays

!Local variables-------------------------------
!scalars
 real(dp) :: field_x, field_y, field_z
 real(dp) :: TX, TY, TZ
 integer :: loop, ii, jj, kk, size
!arrays

! *************************************************************************

 size=FX*FY*FZ

 !Brillouin zone range
 field_x=0.5*FX
 field_y=0.5*FY
 field_z=0.5*FZ
 loop=1

 !Determine new K point states
 do ii=0, (FX-1)
   if ((XX+ii)>(field_x)) then !if falls outside of field, loop it around
     TX=XX-(field_x*2)+ii
   else
     TX=XX+ii
   end if
   do jj=0, (FY-1)
     if ((YY+jj)>(field_y)) then
       TY=YY-(field_y*2)+jj
     else
       TY=YY+jj
     end if
     do kk=0,(FZ-1)
       if ((ZZ+kk)>(field_z)) then
         TZ=ZZ-(field_z*2)+kk
       else
         TZ=ZZ+kk
       end if
       NKVal(1,loop)=TX !Assign to an output array
       NKVal(2,loop)=TY
       NKVal(3,loop)=TZ
       loop=loop+1
     end do
   end do
 end do

 !reduce the values to fit in Brillouin zone
 do ii=1, FX*FY*FZ
   nkval(1,ii)=nkval(1,ii)/FX
   nkval(2,ii)=nkval(2,ii)/FY
   nkval(3,ii)=nkval(3,ii)/FZ
 end do

 END SUBROUTINE NewK
!!***

!!****f* m_fold2block/getargs
!! NAME
!! getargs
!!
!! FUNCTION
!! Part of fold2Bloch that interprets the command line
!! arguments and checks for their corectness.
!!
!! INPUTS
!! fname: Name of input file
!! folds: Number of folds in x, y, and z directions.
!!
!! OUTPUT
!! folds: Number of folds in x, y, and z directions.
!! Returns the folds array filled with each element.
!!
!! PARENTS
!!      fold2Bloch
!!
!! CHILDREN
!!
!! SOURCE

subroutine getargs(folds, fname)

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen), intent(inout) :: fname
!arrays
 integer, intent(inout) :: folds(3)

!Local variables-------------------------------
!scalars
 integer :: num_args, argcount, ii, ios
 character(len=fnlen) :: argfolds
 logical :: dir
!arrays
 character(len=fnlen), allocatable :: args(:)

! *************************************************************************

 !get number of args
 num_args=command_argument_count()
 ABI_MALLOC(args,(num_args))

 !Check for errors in arguments
 if (num_args>2) then
   write(std_out,*) "     Too many arguments."
   write(std_out,*) "     Usage: $fold2Bloch file_WFK x:y:z (folds)"
   ABI_ERROR("Aborting now")
 elseif (num_args<2) then
   if (num_args==1) then
     call get_command_argument(1,args(1))
     if ((args(1)=="-h").or.(args(1)=="--help")) then
       write(std_out,*) "     Usage: fold2Bloch file_WFK x:y:z (folds)"
       write(std_out,*) "     file_WFK is the WFK file name (ex. Ex_WFK)"
       write(std_out,*) "     x:y:z integers, greater than 0 that represent a multiplicity"
       write(std_out,*) "     in the corresponding directions used when constructing the supercell."
       ABI_ERROR("Aborting now")
     else
       write(std_out,*) "     Not all arguments are present."
       write(std_out,*) "     Make sure that file name and number of folds are indicated."
       write(std_out,*) "     Usage: $fold2Bloch file_WFK x:y:z (folds)"
       ABI_ERROR("Aborting now")
     end if
   else
     write(std_out,*) "     Not all arguments are present."
     write(std_out,*) "     Make sure that file name and number of folds are indicated."
     write(std_out,*) "     Usage: $fold2Bloch file_WFK x:y:z (folds)"
     ABI_ERROR("Aborting now")
   end if
 else
   do argcount=1, num_args
     call get_command_argument(argcount,args(argcount))
     if (argcount==1) then
       read(args(argcount), "(a)") fname !read first argument
     else
       read(args(argcount), "(a)") argfolds !read second argument
     end if
   end do
 end if

 ! Does intput file exist?
 inquire(file=fname, exist=dir)
 if (.not.(dir)) then
   write(std_out,*) "     Case file not found: ", trim(fname)
   write(std_out,*) "     Usage: $fold2Bloch file_WFK x:y:z (folds)"
   ABI_ERROR("Aborting now")
 end if

 !Was the number of folds entered in correct format?
 ii=0
 ii=INDEX(argfolds, ':') !look for first ":" to differentiate between axis
 if (ii==0) then
   write(std_out,*) "     Unknown number of folds. See below or type:"
   write(std_out,*) "     fold2Bloch <-h> or fold2Bloch <--help> for more information."
   write(std_out,*) '     Usage: $fold2Bloch file_WFK x:y:z (folds)'
   ABI_ERROR("Aborting now")
 end if
 read (argfolds(1:ii-1), *, iostat=ios) folds(1) !read X folds
 if ((ios/=0).or.(folds(1)<=0)) then
   ABI_ERROR('Number of folds has to be a positive integer greater than 0')
 end if
 argfolds=argfolds(ii+1:) !Start argfolds from the first ":"
 ii=0
 ii=INDEX(argfolds, ':') !look for second ":"
 if (ii==0) then
   write(std_out,*) '     Unknown number of folds. See below or type:'
   write(std_out,*) "     fold2Bloch <-h> or fold2Bloch <--help> for more information."
   write(std_out,*) '     Usage: $fold2Bloch file_WFK x:y:z (folds)'
   ABI_ERROR("Aborting now")
 end if
 read (argfolds(1:ii-1),*, iostat=ios) folds(2) !read Y folds
 if ((ios/=0).or.(folds(2)<=0)) then
   ABI_ERROR('Number of folds has to be a positive integer greater than 0')
 end if
 read(argfolds(ii+1:),*, iostat=ios) Folds(3) !read Z folds
 if ((ios/=0).or.(folds(3)<=0)) then
   ABI_ERROR('Number of folds has to be a positive integer greater than 0')
 end if

 ABI_FREE(args)

end subroutine getargs
!!***

!!****f* m_fold2block/progress
!! NAME
!! progress
!!
!! FUNCTION
!! Part of fold2Bloch that unfolds the wavefunction and calculates
!! the weight of each band.
!!
!! INPUTS
!! ikpt: K point index
!! nkpt: Total number of K points in the function
!! kpt: Current K point being procesed.
!!
!! OUTPUT
!! Write to screen K point being processed and percent complete.
!!
!! PARENTS
!!      fold2Bloch
!!
!! CHILDREN
!!
!! SOURCE

 subroutine progress(ikpt, nkpt, kpt)

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ikpt, nkpt
!arrays
 real(dp), intent(in) :: kpt(3)

! *************************************************************************

 write(std_out, '(a,i10, a1)', advance='no') ''//achar(27)//'[97m',100*ikpt/nkpt,'%'//achar(27)//'[0m'
 write(std_out, '(a25,3f12.6,a)') ''//achar(27)//'[32m Processing K point: ', kpt,''//achar(27)//'[0m'

end subroutine progress
!!***

end module  m_fold2block
