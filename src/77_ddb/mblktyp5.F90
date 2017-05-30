!{\src2tex{textfont=tt}}
!!****f* ABINIT/mblktyp5
!!
!! NAME
!! mblktyp5
!!
!! FUNCTION
!! This routine merges the derivative databases of type 5:
!! second-order eigenvalue derivatives
!!   why is this separate from mblktyp1? Should be merged at some point for consistency
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
!!
!! INPUTS
!!     chkopt=option for consistency checks between DDB files
!!     codename=MRGDDB
!!     ddbun=define input and output unit numbers
!!     dscrpt=description of the output file
!!     filnam=name of input or output file
!!     mddb=maximum number of databases (cannot be made dynamic)
!!     nddb=number of input DDBs
!!
!! OUTPUT
!!     msym=maximum number of symmetry elements in space group
!!     Merge the file
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!      ddb_compare,ddb_free,ddb_getdims,ddb_io_out,ddb_malloc,ddb_write_blok
!!      ioddb8_in,pawtab_free,pawtab_nullify,psddb8,read_blok8,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mblktyp5 (chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)

 use defs_basis
 use m_xmpi
 use m_errors
 use m_profiling_abi
 use m_ddb

 use m_pawtab, only : pawtab_type,pawtab_nullify,pawtab_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mblktyp5'
 use interfaces_14_hidewrite
 use interfaces_72_response
 use interfaces_77_ddb, except_this_one => mblktyp5
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ddbun,mddb,nddb,vrsddb
 integer,intent(out) :: msym
 character(len=fnlen),intent(in) :: dscrpt,filnam(mddb+1)

!Local variables -------------------------
!scalars
!Define input and output unit numbers:
 integer,parameter :: ddbuntmp=3
 integer :: chkopt,choice,dimekb,dimekb_tmp,fullinit,fullmrgddb_init,iblok,iblok1,iblok2
 integer :: iddb,ii,intxc,intxc8,iscf,iscf8,ixc,ixc8,lmnmax,lmnmax_tmp,matom
 integer :: mband,mband_tmp,mblktyp,mblok,mkpt,mpert,msize,mtypat,natom
 integer :: natom8,nblok,nblok8,nblokt,nkpt,nkpt8,nspden,nspden8
 integer :: nspinor,nspinor8,nsppo8,nsppol,nsym,nsym8,ntypat,ntypat8,nunit
 integer :: occop8,occopt,temp,tmerge,usepaw,usepaw_tmp,useylm
 integer :: msym_tmp
 integer :: ngfft(18),ngfft8(18)
 integer, allocatable :: symafm(:),symafm8(:)
 integer, allocatable :: symre8(:,:,:),symrel(:,:,:)
 integer,allocatable :: indlmn(:,:,:),lloc(:),mgblok(:)
 integer,allocatable :: nband(:),nband8(:),pspso(:),typat(:),typat8(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff,dilatmx,dilatmx8,ecut,ecut8,ecutsm,ecutsm8,kptnr8,kptnrm
 real(dp) :: pawecutdg,pawecutdg8,dfpt_sciss,dfpt_sciss8,tolwf8,tolwfr,tphysel
 real(dp) :: tphysel8,tsmear,tsmear8
 type(ddb_type) :: ddb
!arrays
 real(dp) :: acell(3),acell8(3),rprim(3,3),rprim8(3,3)
 real(dp), allocatable :: tnons(:,:)
 real(dp), allocatable :: tnons8(:,:)
 real(dp),allocatable :: amu(:),amu8(:)
 real(dp),allocatable :: blkval2(:,:,:,:)
 real(dp),allocatable :: ekb(:,:),ekb8(:,:),kpt(:,:),kpt8(:,:),kpnt(:,:,:),occ(:),occ8(:)
 real(dp),allocatable :: spinat(:,:),spinat8(:,:),vel(:,:),wtk(:),wtk8(:)
 real(dp),allocatable :: xcart(:,:),xred(:,:),xred8(:,:),zion(:),zion8(:)
 real(dp),allocatable :: znucl(:),znucl8(:)
 character(len=500) :: message
 type(pawtab_type),allocatable :: pawtab(:),pawtab8(:)

! *********************************************************************

!Evaluate the maximal dimensions of arrays
 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0
 msize=0  ; mtypat=0 ; lmnmax=0 ; usepaw=0
 msym = 192
 do iddb=1,nddb
   call ddb_getdims(dimekb_tmp,filnam(iddb+1),lmnmax_tmp,mband_tmp,mblktyp,&
&   msym_tmp,natom,nblok,nkpt,ntypat,ddbun,usepaw_tmp,vrsddb,xmpi_comm_self)

   dimekb=max(dimekb,dimekb_tmp)
   lmnmax=max(lmnmax,lmnmax_tmp)
   matom=max(matom,natom)         ! MG Why this! I dont' understand why we always like to complicate things!
   mband=max(mband,mband_tmp)
   mblok=mblok+nblok
   mkpt=max(mkpt,nkpt)
   mtypat=max(mtypat,ntypat)
   usepaw=max(usepaw,usepaw_tmp)
   msym = max(msym_tmp, msym)
 end do

 ABI_ALLOCATE(symafm,(msym))
 ABI_ALLOCATE(symafm8,(msym))
 ABI_ALLOCATE(symre8,(3,3,msym))
 ABI_ALLOCATE(symrel,(3,3,msym))
 ABI_ALLOCATE(tnons,(3,msym))
 ABI_ALLOCATE(tnons8,(3,msym))

 mpert=matom+6
 msize=3*mpert*3*mpert
 call ddb_malloc(ddb,msize,mblok,matom,mtypat)

!Allocate arrays
 ABI_ALLOCATE(lloc,(mtypat))
 ABI_ALLOCATE(mgblok,(mblok))
 ABI_ALLOCATE(nband,(mkpt))
 ABI_ALLOCATE(nband8,(mkpt))
 ABI_ALLOCATE(typat,(matom))
 ABI_ALLOCATE(typat8,(matom))
 ABI_ALLOCATE(amu,(mtypat))
 ABI_ALLOCATE(amu8,(mtypat))
 ABI_ALLOCATE(ekb,(dimekb,mtypat))
 ABI_ALLOCATE(ekb8,(dimekb,mtypat))
 ABI_ALLOCATE(kpt,(3,mkpt))
 ABI_ALLOCATE(kpt8,(3,mkpt))
 ABI_ALLOCATE(indlmn,(6,lmnmax,mtypat))
 ABI_ALLOCATE(pspso,(mtypat))
 ABI_ALLOCATE(occ,(mband*mkpt))
 ABI_ALLOCATE(occ8,(mband*mkpt))
 ABI_ALLOCATE(spinat,(3,matom))
 ABI_ALLOCATE(spinat8,(3,matom))
 ABI_ALLOCATE(vel,(3,matom))
 ABI_ALLOCATE(wtk,(mkpt))
 ABI_ALLOCATE(wtk8,(mkpt))
 ABI_ALLOCATE(xcart,(3,matom))
 ABI_ALLOCATE(xred8,(3,matom))
 ABI_ALLOCATE(xred,(3,matom))
 ABI_ALLOCATE(znucl,(mtypat))
 ABI_ALLOCATE(znucl8,(mtypat))
 ABI_ALLOCATE(zion,(mtypat))
 ABI_ALLOCATE(zion8,(mtypat))
 ABI_DATATYPE_ALLOCATE(pawtab,(mtypat*usepaw))
 ABI_DATATYPE_ALLOCATE(pawtab8,(mtypat*usepaw))

 call pawtab_nullify(pawtab)
 call pawtab_nullify(pawtab8)

!This is needed to read the DDBs in the old format
 symafm(:)=1 ; symafm8(:)=1
 if(mtypat>=1)then
   pspso(:)=0
   znucl(:)=zero ; znucl8(:)=zero
   ekb(:,:)=zero ; ekb8(:,:)=zero
 end if
 if(matom>=1)then
   spinat(:,:)=zero ; spinat8(:,:)=zero
 end if

!**********************************************************************

!Read the first database
 if (nddb==1) then

   write(message, '(a,a,a,a,a)' )&
&   'The initialisation mode of MRGDDB, that uses nddb=1,',&
&   'has been disabled in version 2.1 of ABINIT.',&
&   'Action : you should use DDBs that include the symmetry',&
&   'information (and that can be used and merged without',&
&   'initialisation), or you should use ABINITv2.0.'
   MSG_ERROR(message)

 else if(nddb>=2)then

!  Open the first derivative database file
!  and read the preliminary information
   write(std_out,*)' read the input derivative database information'
   nunit=ddbun
   call ioddb8_in (filnam(2),matom,mband,&
&   mkpt,msym,mtypat,nunit,vrsddb,&
&   acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&   natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&   pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,&
&   tphysel,tsmear,typat,usepaw,wtk,xred,zion,znucl)

   ABI_ALLOCATE(blkval2,(2,msize,nband(1),mkpt))
   ABI_ALLOCATE(kpnt,(3,mkpt,mblok))


!  Read the psp information of the input DDB
   useylm=usepaw;choice=1
   call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&   nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

   if(nblok>=1)then
!    Read the blocks from the input database.
     write(message, '(a,i5,a)' ) ' read ',nblok,' blocks from the input DDB '
     call wrtout(std_out,message,'COLL')
     choice=1
     do iblok=1,nblok
       call read_blok8(ddb,iblok,nband(1),mpert,&
&       msize,nkpt,nunit,blkval2(1,1,1,1),kpnt(1,1,iblok))
!      Setup merged indicator
       mgblok(iblok)=1
     end do
   else
     call wrtout(std_out,' No bloks in the first ddb ','COLL')
   end if
!  Close the first ddb
   close(ddbun)

 end if

!*********************************************

!In case of merging of DDBs, iterate the reading
 do iddb=2,nddb

!  Open the corresponding input DDB, and read the database file informations
   write(message, '(a,a,i6)' )ch10,&
&   ' read the input derivative database number',iddb
   call wrtout(std_out,message,'COLL')
   nunit=ddbun
   call ioddb8_in (filnam(iddb+1),matom,mband,&
&   mkpt,msym,mtypat,nunit,vrsddb,&
&   acell8,amu8,dilatmx8,ecut8,ecutsm8,intxc8,iscf8,ixc8,kpt8,kptnr8,&
&   natom8,nband8,ngfft8,nkpt8,nspden8,nspinor8,nsppo8,nsym8,ntypat8,occ8,occop8,&
&   pawecutdg8,rprim8,dfpt_sciss8,spinat8,symafm8,symre8,tnons8,tolwf8,&
&   tphysel8,tsmear8,typat8,usepaw,wtk8,xred8,zion8,znucl8)

!  Read the psp information of the input DDB
   choice=1
   call psddb8 (choice,dimekb,ekb8,fullmrgddb_init,indlmn,lmnmax,&
&   nblok8,ntypat8,nunit,pawtab8,pspso,usepaw,useylm,vrsddb)

   if (chkopt==1)then
!    Compare the current DDB and input DDB information.
!    In case of an inconsistency, halt the execution.
     call wrtout(std_out,' compare the current and input DDB information','COLL')

!    Should also compare indlmn and pspso ... but suppose that
!    the checking of ekb is enough for the psps.
!    Should also compare many other variables ... this is still to be done ...
     call ddb_compare (acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&
&     fullinit,fullmrgddb_init,iscf,iscf8,ixc,ixc8,kpt,kpt8,kptnrm,kptnr8,&
&     natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&
&     nsppol,nsppo8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&
&     occopt,occop8,pawecutdg,pawecutdg8,pawtab,pawtab8,&
&     rprim,rprim8,dfpt_sciss,dfpt_sciss8,symrel,symre8,&
&     tnons,tnons8,tolwfr,tolwf8,typat,typat8,&
&     usepaw,wtk,wtk8,xred,xred8,zion,zion8)
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
   write(message, '(a,i5,a)' )' Will read ',nblok8,' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   nblokt=nblok+nblok8
   if(nblokt>mblok)then
     write(message, '(a,i0,a,a,a,i0,a)' )&
&     'The expected number of blocks',nblokt,' is larger than',ch10,&
&     'the maximum number of blocks',mblok,'.'
     MSG_ERROR(message)
   end if

!  Read the bloks from the temporary database, and close it.
!  Also setup the merging indicator
   choice=1
   nunit=ddbun
   do iblok=nblok+1,nblokt
     call read_blok8(ddb,iblok,nband(1),mpert,&
&     msize,nkpt8,nunit,blkval2(1,1,1,1),kpnt(1,1,iblok))
     mgblok(iblok)=1
   end do
   close(ddbun)

   nblok=nblokt
   write(message, '(a,i5)' ) ' Now, current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')
 end do

 call wrtout(std_out,' All DDBs have been read ','COLL')

!*********************************************************

!Check the equality of blocks, and eventually merge them

 if(nblok>=1)then
   call wrtout(std_out,' check the equality of blocks, and eventually merge ','COLL')
   do iblok2=2,nblok
     do iblok1=1,iblok2-1
!      Check the block type identity
       if(ddb%typ(iblok1)==ddb%typ(iblok2))then
!        Check the wavevector identities
         diff=abs(ddb%qpt(1,iblok1)-ddb%qpt(1,iblok2))
         diff=diff+abs(ddb%qpt(2,iblok1)-ddb%qpt(2,iblok2))
         diff=diff+abs(ddb%qpt(3,iblok1)-ddb%qpt(3,iblok2))
!        DEBUG
!        write(std_out,*),'ddb%qpt(1,iblok1)=',ddb%qpt(1,iblok1)
!        write(std_out,*),'ddb%qpt(1,iblok2)=',ddb%qpt(1,iblok2)
!        write(std_out,*),'ddb%qpt(2,iblok1)=',ddb%qpt(2,iblok1)
!        write(std_out,*),'ddb%qpt(2,iblok2)=',ddb%qpt(2,iblok2)
!        write(std_out,*),'ddb%qpt(3,iblok1)=',ddb%qpt(3,iblok1)
!        write(std_out,*),'ddb%qpt(3,iblok2)=',ddb%qpt(3,iblok2)
!        ENDDEBUG
         if(abs(diff)<qtol)mgblok(iblok2)=0
       end if
     end do
   end do
!  DEBUG
!  write(std_out,*)'ddb%qpt',ddb%qpt
!  write(std_out,*)'mgblok=',mgblok
!  ENDDEBUG

!  Count the final number of bloks
   tmerge=0
   do ii=1,nblok
     if(mgblok(ii)==1)tmerge=tmerge+1
   end do
   temp = nblok-tmerge
   nblok=tmerge

!  Summarize the merging phase
   write(message, '(i6,a,i6,a)' )&
&   temp,' blocks are merged; the new DDB will have ',nblok,' blocks.'
   call wrtout(std_out,message,'COLL')
 end if

!**********************************************************************

!Open the output database, then
!Write the preliminary informations
 write(message, '(a,a)' )' open the output database, write the',' preliminary information '
 call wrtout(std_out,message,'COLL')

 nunit=ddbun
 call ddb_io_out (dscrpt,filnam(1),matom,mband,&
& mkpt,msym,mtypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,&
& tphysel,tsmear,typat,usepaw,wtk,xred,zion,znucl)


!Write the psp information in the output DDB
!as well as the value of the number of blocks.
 call wrtout(std_out,' write the psp information ','COLL')
 fullinit=1 ; choice=2
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
& nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!DEBUG
!write(std_out,*)' mgblok= ',mgblok
!ENDDEBUG

 if(nddb>1)then

!  Write the whole database
   call wrtout(std_out,' write the DDB ','COLL')
   ii = 1 !unit indicator of what will be merged
!  Create a temporary file to decrease memory need.
   do iddb=1,nddb
     choice=1
     call ioddb8_in (filnam(iddb+1),matom,mband,mkpt,msym,mtypat,&
&     ddbuntmp,vrsddb,acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,&
&     kptnrm,natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,&
&     occ,occopt,pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,&
&     tphysel,tsmear,typat,usepaw,wtk,xred,zion,znucl)
     call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&     nblok,ntypat,ddbuntmp,pawtab,pspso,usepaw,useylm,vrsddb)
     do iblok=1,nblok
       if(mgblok(ii)==1) then
         call read_blok8(ddb,ii,nband(1),mpert,&
&         msize,nkpt,ddbuntmp,blkval2(:,:,:,:),kpnt(:,:,ii))
         choice=2
         call ddb_write_blok(ddb,ii,choice,nband(1),mpert,&
&         msize,nkpt,nunit,blkval2(:,:,:,:),kpnt(:,:,ii))
       else
         write(message, '(a,i4,a,i4,a)' )&
&         ' Bloc number',iblok,' of DDB ',iddb,&
&         ' was merged, so do not write it'
         call wrtout(std_out,message,'COLL')
       end if
       ii = ii+1
     end do
     close(ddbuntmp)
   end do !iddb=1,nddb

!  Also write summary of bloks at the end
   write(ddbun, '(/,a)' )' List of bloks and their characteristics '
   choice=3
   nunit=ddbun
   do iblok=1,nblokt
     if(mgblok(iblok)==1)then
       call ddb_write_blok(ddb,iblok,choice,nband(1),mpert,&
&       msize,nkpt,nunit)
     end if
   end do

 end if

 close (ddbun)

!*********************************************************************

!Deallocate arrays
 ABI_DEALLOCATE(lloc)
 ABI_DEALLOCATE(mgblok)
 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(nband8)
 ABI_DEALLOCATE(typat)
 ABI_DEALLOCATE(typat8)
 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(amu8)
 ABI_DEALLOCATE(blkval2)
 ABI_DEALLOCATE(ekb)
 ABI_DEALLOCATE(ekb8)
 ABI_DEALLOCATE(kpt)
 ABI_DEALLOCATE(kpt8)
 ABI_DEALLOCATE(kpnt)
 ABI_DEALLOCATE(indlmn)
 ABI_DEALLOCATE(pspso)
 ABI_DEALLOCATE(occ)
 ABI_DEALLOCATE(occ8)
 ABI_DEALLOCATE(spinat)
 ABI_DEALLOCATE(spinat8)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(wtk)
 ABI_DEALLOCATE(wtk8)
 ABI_DEALLOCATE(xcart)
 ABI_DEALLOCATE(xred)
 ABI_DEALLOCATE(xred8)
 ABI_DEALLOCATE(znucl)
 ABI_DEALLOCATE(znucl8)
 ABI_DEALLOCATE(zion)
 ABI_DEALLOCATE(zion8)
 ABI_DEALLOCATE(symafm)
 ABI_DEALLOCATE(symafm8)
 ABI_DEALLOCATE(symre8)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)
 ABI_DEALLOCATE(tnons8)

 call pawtab_free(pawtab)
 ABI_DATATYPE_DEALLOCATE(pawtab)
 call pawtab_free(pawtab8)
 ABI_DATATYPE_DEALLOCATE(pawtab8)
 call ddb_free(ddb)

end subroutine mblktyp5
!!***
