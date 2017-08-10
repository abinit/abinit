!{\src2tex{textfont=tt}}
!!****f* ABINIT/mblktyp1
!!
!! NAME
!! mblktyp1
!!
!! FUNCTION
!! This routine merges the derivative databases of type 0-4:
!! Total energy, (2nd derivatives (non-stat.),2nd derivatives (stationary),
!! 3rd derivatives, 1st derivatives
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG,MT,SP)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! The heading of the database is read, then the heading
!! of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in central memory.
!!
!! INPUTS
!!     chkopt=option for consistency checks between DDB files
!!     ddbun=define input and output unit numbers
!!     dscrpt=description of the output file
!!     filnam=name of input or output file
!!     mddb=maximum number of databases (cannot be made dynamic)
!!     nddb=number of input DDBs
!!     vrsddb=current version of the DDB
!!
!! OUTPUT
!!     msym=maximum number of symmetry elements in space group
!!     Merge the file
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!      ddb_free,ddb_malloc,ddb_write_blok
!!      read_blok8,wrtout
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mblktyp1(chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_ddb
 use m_ddb_hdr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mblktyp1'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chkopt,ddbun,mddb,nddb,vrsddb
 integer,intent(out) :: msym
 character(len=fnlen),intent(in) :: dscrpt
 character(len=fnlen),intent(in) :: filnam(mddb+1)

!Local variables -------------------------
!scalars
!Define input and output unit numbers:
 integer :: choice,dimekb,iblok,iblok1,iblok2
 integer :: iddb,ii,lmnmax,matom
 integer :: mband,mblktyp,mblok,mkpt,mpert,msize,mtypat
 integer :: nblok,nblokt,nq
 integer :: tmerge,usepaw
 integer :: ngfft(18),ngfft8(18)
 integer,allocatable :: mgblok(:)!,lloc(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff
 type(ddb_type) :: ddb
 type(ddb_hdr_type) :: ddb_hdr, ddb_hdr8
!arrays
 character(len=500) :: message

! *********************************************************************

 ! Make sure there is more than one ddb to be read
 if(nddb==1)then

   write(message, '(a,a,a,a,a)' )&
&   'The initialisation mode of MRGDDB, that uses nddb=1,',&
&   'has been disabled in version 2.1 of ABINIT.',&
&   'Action : you should use DDBs that include the symmetry',&
&   'information (and that can be used and merged without',&
&   'initialisation), or you should use ABINITv2.0.'
   MSG_ERROR(message)

 end if

!Evaluate the maximal dimensions of arrays
 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0
 msize=0  ; mtypat=0 ; lmnmax=0 ; usepaw=0 ; mblktyp = 1
 msym=192

 do iddb=1,nddb

   call ddb_hdr_open_read(ddb_hdr, filnam(iddb+1), ddbun, vrsddb,&
&                         dimonly=1)

   mblok=mblok+ddb_hdr%nblok
   mblktyp=max(mblktyp,ddb_hdr%mblktyp)
   matom=max(matom,ddb_hdr%matom)
   mkpt=max(mkpt,ddb_hdr%mkpt)
   mtypat=max(mtypat,ddb_hdr%mtypat)
   msym=max(msym,ddb_hdr%msym)
   mband=max(mband,ddb_hdr%mband)
   dimekb=max(dimekb,ddb_hdr%psps%dimekb)
   lmnmax=max(lmnmax,ddb_hdr%psps%lmnmax)
   usepaw=max(usepaw,ddb_hdr%usepaw)

   call ddb_hdr_free(ddb_hdr)

 end do

 mpert=matom+6
 msize=3*mpert*3*mpert
 if(mblktyp==3)msize=msize*3*mpert

 call ddb_malloc(ddb,msize,mblok,matom,mtypat)

!Allocate arrays
 ABI_ALLOCATE(mgblok,(mblok))

!**********************************************************************

!Read the first database

 write(std_out,*)' read the input derivative database information'
 call ddb_hdr_open_read(ddb_hdr, filnam(2), ddbun, vrsddb, &
&            matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
&            msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

 if(ddb_hdr%nblok>=1)then
!  Read the blocks from the input database.
   write(message, '(a,i5,a)' ) ' read ',ddb_hdr%nblok, &
&                              ' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   do iblok=1,ddb_hdr%nblok
     call read_blok8(ddb,iblok,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
!    Setup merged indicator
     mgblok(iblok)=0
   end do
 else
   write(message, '(a)' )' No bloks in the first ddb '
   call wrtout(std_out,message,'COLL')
 end if
!Close the first ddb
 close(ddbun)

!*********************************************

 nblok = ddb_hdr%nblok
!In case of merging of DDBs, iterate the reading
 do iddb=2,nddb

!  Open the corresponding input DDB,
!  and read the database file informations
   write(message, '(a,a,i6)' )ch10,&
&       ' read the input derivative database number',iddb
   call wrtout(std_out,message,'COLL')

   call ddb_hdr_open_read(ddb_hdr8, filnam(iddb+1), ddbun, vrsddb, &
&            matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
&            msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

   if (chkopt==1)then
!    Compare the current DDB and input DDB information.
!    In case of an inconsistency, halt the execution.
     write(message, '(a)' )' compare the current and input DDB information'
     call wrtout(std_out,message,'COLL')

     call ddb_hdr_compare(ddb_hdr, ddb_hdr8)

   else if(chkopt==0)then
!    No comparison between the current DDB and input DDB information.
     write(message, '(a)' )' no comparison between the current and input DDB information'
     call wrtout(std_out,message,'COLL')
     write(message, '(a,a,a)' )&
&     'No comparison/check is performed for the current and input DDB information ',&
&     'because argument --nostrict was passed to the command line. ',&
&     'Use at your own risk !'
     MSG_COMMENT(message)
   end if

   call wrtout(std_out,' Will try to merge this input DDB with the current one.','COLL')

!  First estimate of the total number of bloks, and error
!  message if too large
   write(message, '(a,i5)' ) ' Current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i5,a)' )' Will read ',ddb_hdr8%nblok,' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   nblokt=nblok+ddb_hdr8%nblok
   if(nblokt>mblok)then
     write(message, '(a,i5,a,a,a,i5,a)' )&
&     'The expected number of blocks',nblokt,' is larger than',ch10,&
&     'the maximum number of blocks',mblok,'.'
     MSG_ERROR(message)
   end if

!  Read the bloks from the temporary database, and close it.
!  Also setup the merging indicator
   do iblok=nblok+1,nblokt
     call read_blok8(ddb,iblok,ddb_hdr8%nband(1),mpert,msize,ddb_hdr8%nkpt,ddbun)
     mgblok(iblok)=0
   end do
   close(ddbun)

   nblok=nblokt
   write(message, '(a,i5)' ) ' Now, current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')

   ! In certain cases, the different DDB will have different information
   ! on the pseudos (depending on fullinit)
   ! Here, we copy the information of the last DDB file,
   ! only to make the tests pass...
   ddb_hdr%psps%indlmn(:,:,:) = ddb_hdr8%psps%indlmn(:,:,:)
   ddb_hdr%psps%pspso(:) = ddb_hdr8%psps%pspso(:)
   ddb_hdr%psps%ekb(:,:) = ddb_hdr8%psps%ekb(:,:)

   call ddb_hdr_free(ddb_hdr8)

 end do


 call wrtout(std_out,' All DDBs have been read ','COLL')

!*********************************************************

!Check the equality of blocks, and eventually merge them

 if(nblok>=1)then
   call wrtout(std_out,' check the equality of blocks, and eventually merge ','COLL')
   do iblok2=2,nblok
     do iblok1=1,iblok2-1
       tmerge=0

!      Check the block type identity
       if(ddb%typ(iblok1)==ddb%typ(iblok2))then

!        Check the wavevector identities
         tmerge=1
         if(ddb%typ(iblok1)==1.or.ddb%typ(iblok1)==2)then
           nq=1
         else if(ddb%typ(iblok1)==3)then
!          Note : do not merge permutation related elements ....
           nq=3
         else if(ddb%typ(iblok1)==4 .or. ddb%typ(iblok1)==0)then
           nq=0
         end if
         if(nq/=0)then
           do ii=1,nq
             diff=ddb%qpt(1+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1)&
&             -ddb%qpt(1+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if(abs(diff)>qtol)tmerge=0
             diff=ddb%qpt(2+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1)&
&             -ddb%qpt(2+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if(abs(diff)>qtol)tmerge=0
             diff=ddb%qpt(3+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1)&
&             -ddb%qpt(3+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if(abs(diff)>qtol)tmerge=0
           end do ! ii
         end if

!        Now merges,
         if(tmerge==1)then
           write(message, '(a,i5,a,i5)' )' merge block #',iblok2,' to block #',iblok1
           call wrtout(std_out,message,'COLL')
           mgblok(iblok2)=1
           do ii=1,msize
             if(ddb%flg(ii,iblok2)==1)then
               ddb%flg(ii,iblok1)=1
               ddb%val(1,ii,iblok1)=ddb%val(1,ii,iblok2)
               ddb%val(2,ii,iblok1)=ddb%val(2,ii,iblok2)
             end if
           end do
         end if

       end if
     end do
   end do

!  Count the final number of bloks
   tmerge=0
   do ii=1,nblok
     if(mgblok(ii)==1)tmerge=tmerge+1
   end do
   nblok=nblok-tmerge

!  Summarize the merging phase
   write(message, '(i6,a,i6,a)' )&
&   tmerge,' blocks are merged; the new DDB will have ',nblok,' blocks.'
   call wrtout(std_out,message,'COLL')

!  End the condition on existence of more than one blok in current DDB
 end if

!**********************************************************************

 write(message, '(a,a)' )' open the output database, write the',' preliminary information '
 call wrtout(std_out,message,'COLL')

 ddb_hdr%dscrpt = trim(dscrpt)
 ddb_hdr%nblok = nblok
 ddb_hdr%mblktyp = mblktyp

 call ddb_hdr_open_write(ddb_hdr, filnam(1), ddbun, fullinit=1)

 if(nddb>1)then

!  Write the whole database
   call wrtout(std_out,' write the DDB ','COLL')
   choice=2
   do iblok=1,nblok+tmerge
     if(mgblok(iblok)==0)then
       write(std_out,'(a,i4)' ) ' Write bloc number',iblok
       call ddb_write_blok(ddb,iblok,choice,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
     else
       write(message, '(a,i4,a)' )&
&       ' Bloc number',iblok,' was merged, so do not write it'
       call wrtout(std_out,message,'COLL')
     end if
   end do

!  Also write summary of bloks at the end
   write(ddbun, '(/,a)' )' List of bloks and their characteristics '
   choice=3
   do iblok=1,nblok+tmerge
     if(mgblok(iblok)==0)then
       call ddb_write_blok(ddb,iblok,choice,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
     end if
   end do

 end if

 close (ddbun)

!*********************************************************************

!Deallocate arrays

 ABI_DEALLOCATE(mgblok)

 call ddb_hdr_free(ddb_hdr)
 call ddb_free(ddb)

end subroutine mblktyp1
!!***
