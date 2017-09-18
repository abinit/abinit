!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpspheads
!! NAME
!! inpspheads
!!
!! FUNCTION
!! Read the pseudopotential header of each psp file, in order to initialize pspheads(1:npsp).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, FrD, AF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!  pseudopotential file headers, as well as the psp file names
!!  ecut_tmp(3,2,npsp)= possible ecut values as read in psp files
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      atomic_info,getecutfromxml,paw_setup_copy,paw_setup_free,pawpsxml2ab
!!      psp_from_data,psxml2abheader,rdpawpsxml,upfheader2abi,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine inpspheads(filnam,npsp,pspheads,ecut_tmp)

 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_errors
 use m_hash_md5

 use m_io_tools, only : open_file
 use m_fstrings, only : basename, lstrip, sjoin, startswith
 use m_pawxmlps, only : paw_setup, paw_setuploc, npsp_pawxml, ipsp2xml, rdpawpsxml, &
&                       paw_setup_copy, paw_setup_free, getecutfromxml
 use m_psxml2ab

#if defined HAVE_PSML
 use m_psml
#endif
#if defined HAVE_BIGDFT
  use BigDFT_API, only: atomic_info,psp_from_data
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inpspheads'
 use interfaces_14_hidewrite
 use interfaces_57_iopsp_parser, except_this_one => inpspheads
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
!arrays
 real(dp),intent(inout) :: ecut_tmp(3,2,10)
 character(len=fnlen), intent(in) :: filnam(npsp)
 type(pspheader_type),intent(inout) :: pspheads(npsp) !vz_i

!Local variables-------------------------------
!In case a xc core correction is to be taken into account,
!the n1xccc value will be given by n1xccc_default. Otherwise it is set to 0.
!scalars
 integer,parameter :: n1xccc_default=2501
 integer :: extension_switch
 integer :: idum,ii,ilmax,ipsp,ipsp_pawxml,lang,lmax,mmax,mpsang,n1xccc,nmesh,npsp_pawxml0
 integer :: pspcod,pspso,test_paw,usexml,unt,useupf !,pspxc
 real(dp) :: al,e990,e999,fchrg,qchrg,r1,rchrg,rr ! ,rp,rs
 character(len=3) :: testxc
 character(len=500) :: message,errmsg
 character(len=70) :: testxml
 character(len=80) :: pspline
!arrays
 integer :: nproj(0:3),nprojso(1:3)
 integer,allocatable :: orb(:)
 real(dp) :: hdum(3)
#if defined HAVE_BIGDFT
 !new variables for wvl+paw
 character(len=2) :: symbol
 integer :: iasctype,nzatom, nelpsp, npspcode_,ixc_
 real(dp) :: rcov,ehomo
 real(dp) :: psppar(0:4,0:6)
 logical :: exists
#endif
#if defined HAVE_PSML
 character(len=3) :: atmsymb
 character(len=30) :: creator
#endif

!*************************************************************************

 test_paw=0
 npsp_pawxml0=0
 ipsp_pawxml=0

 do ipsp=1,npsp
   if (open_file(filnam(ipsp),message,newunit=unt,form="formatted",status="old") /= 0) then
     MSG_ERROR(message)
   end if

   rewind(unit=unt, err=10, iomsg=errmsg)
   read(unt,*, err=10, iomsg=errmsg) testxml
   if (testxml(1:5)=='<?xml') then
     read(unt,*, err=10, iomsg=errmsg) testxml
     if(testxml(1:4)=='<paw')npsp_pawxml0 =npsp_pawxml0+ 1
   end if
   close(unt, err=10, iomsg=errmsg)
 end do

 npsp_pawxml=npsp_pawxml0
 ABI_DATATYPE_ALLOCATE(paw_setup,(npsp_pawxml))
 ABI_ALLOCATE(ipsp2xml,(npsp))
 ipsp2xml=0

 do ipsp=1,npsp

!  Check if the file is written in XML
   pspheads(ipsp)%filpsp=trim(filnam(ipsp))

   usexml = 0
   if (open_file(filnam(ipsp),message,newunit=unt,form="formatted",status="old") /= 0) then
     MSG_ERROR(message)
   end if

   rewind(unit=unt, err=10, iomsg=errmsg)
   read(unt,*, err=10, iomsg=errmsg) testxml
   if(testxml(1:5)=='<?xml')then
     usexml = 1
     read(unt,*, err=10, iomsg=errmsg) testxml
     if(testxml(1:4)=='<paw')then
       test_paw = 1
     else
       test_paw = 0
     end if
   else
     usexml = 0
   end if

   close(unit=unt, err=10, iomsg=errmsg)

!  Check if pseudopotential file is a Q-espresso UPF file
   useupf = 0
   if (open_file(filnam(ipsp),message,newunit=unt,form="formatted",status="old") /= 0) then
     MSG_ERROR(message)
   end if

   rewind(unit=unt,err=10,iomsg=errmsg)
   read(unt,*,err=10,iomsg=errmsg) testxml ! just a string, no relation to xml.
   if (testxml(1:9)=='<PP_INFO>') then
     useupf = 1
   else
     useupf = 0
   end if
   close(unit=unt,err=10,iomsg=errmsg)

!  Read the header of the pseudopotential file
   if (usexml /= 1 .and. useupf /= 1) then
!    Open the psp file and read a normal abinit style header
     if (open_file(filnam(ipsp), message, newunit=unt, form='formatted', status='old') /= 0) then
       MSG_ERROR(message)
     end if

     rewind (unit=unt, err=10, iomsg=errmsg)

!    Read the three first lines
     read(unt, '(a)', err=10, iomsg=errmsg) pspheads(ipsp)%title
     read(unt,*, err=10, iomsg=errmsg)pspheads(ipsp)%znuclpsp,pspheads(ipsp)%zionpsp,pspheads(ipsp)%pspdat
     read(unt,*, err=10, iomsg=errmsg)pspheads(ipsp)%pspcod,pspheads(ipsp)%pspxc,pspheads(ipsp)%lmax,idum,mmax

     pspcod=pspheads(ipsp)%pspcod
     lmax=pspheads(ipsp)%lmax
     write(message,'(a,f5.1,a,i4,a,i4)')'  read the values zionpsp=',pspheads(ipsp)%zionpsp,' , pspcod=',pspcod,' , lmax=',lmax
     call wrtout(std_out,message,'PERS')

     nproj(0:3)=0 ; nprojso(1:3)=0

     pspheads(ipsp)%xccc=0
     pspheads(ipsp)%pspso=0

   else if (usexml==1 .and. test_paw==0) then
#if defined HAVE_PSML
     write(message,'(a,a)')  &
&     '- inpspheads : Reading pseudopotential header in XML form from ', trim(filnam(ipsp))
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')

! could pass pspheads(ipsp) directly and fill all of it in psxml2ab
     call psxml2abheader( filnam(ipsp), pspheads(ipsp), atmsymb, creator, 1)

! save some stuff locally for this ipsp
     pspcod = pspheads(ipsp)%pspcod
     lmax   = pspheads(ipsp)%lmax
     nproj = pspheads(ipsp)%nproj
     nprojso = pspheads(ipsp)%nprojso

#else
     write(message, '(a,a)') "XML norm-conserving pseudopotential has been input,", &
&     " but abinit is not compiled with libPSML support. Reconfigure and recompile."
     MSG_ERROR(message)
#endif
     
   else if(usexml==1.and.test_paw==1)then

     ipsp_pawxml=ipsp_pawxml+1
     ipsp2xml(ipsp)=ipsp_pawxml
     write(message,'(a,a)')  &
&     '- inpspheads : Reading pseudopotential header in XML form from ', trim(filnam(ipsp))
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     call getecutfromxml(filnam(ipsp),ecut_tmp(:,:,ipsp),unt)
     call rdpawpsxml(filnam(ipsp),paw_setuploc,unt)
     call pawpsxml2ab( paw_setuploc, pspheads(ipsp),1 )
     call paw_setup_copy(paw_setuploc,paw_setup(ipsp_pawxml))
     pspcod=17
     call paw_setup_free(paw_setuploc)

   else if (useupf == 1) then
     pspheads(ipsp)%pspcod = 11

     pspheads(ipsp)%xccc  = n1xccc_default ! will be set to 0 if no nlcc
!    call upfoctheader2abi (filnam(ipsp),  &
     call upfheader2abi (filnam(ipsp),  &
&     pspheads(ipsp)%znuclpsp, &
&     pspheads(ipsp)%zionpsp,  &
&     pspheads(ipsp)%pspxc,    &
&     pspheads(ipsp)%lmax,     &
&     pspheads(ipsp)%xccc,     &
&     nproj, nprojso)

     pspcod = pspheads(ipsp)%pspcod
     lmax   = pspheads(ipsp)%lmax
!    FIXME : generalize for SO pseudos
     pspheads(ipsp)%pspso = 0
   end if

!  DEBUG
!  write(std_out,*) pspheads(ipsp)%znuclpsp
!  write(std_out,*) pspheads(ipsp)%zionpsp
!  write(std_out,*) pspheads(ipsp)%pspcod
!  write(std_out,*) pspheads(ipsp)%pspxc
!  write(std_out,*) pspheads(ipsp)%lmax
!  stop
!  ENDDEBUG

!  Initialize nproj, nprojso, pspso, as well as xccc, for each type of psp
   pspheads(ipsp)%GTHradii = zero

   if(pspcod==1 .or. pspcod==4)then

!    Teter format
     do ilmax=0,lmax
       read(unt,*, err=10, iomsg=errmsg) lang,e990,e999,nproj(ilmax)
       read(unt,*, err=10, iomsg=errmsg)
     end do
     read(unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==2)then

!    GTH pseudopotentials
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(0) !rloc
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(1),hdum(1),hdum(2)
     if(abs(hdum(1))>1.d-9) nproj(0)=1
     if(abs(hdum(2))>1.d-9) nproj(0)=2
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(2),hdum(3)
     if(abs(hdum(3))>1.d-9) nproj(1)=1

   else if(pspcod==3)then

!    HGH pseudopotentials
     read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(0) !rloc
     do ilmax=0,lmax
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(ilmax + 1),hdum(1),hdum(2),hdum(3)
       if (abs(hdum(1))>1.d-9)nproj(ilmax)=1
       if (abs(hdum(2))>1.d-9)nproj(ilmax)=2
       if (abs(hdum(3))>1.d-9)nproj(ilmax)=3
       if (ilmax>0.and.ilmax<3) then
         read (unt,*, err=10, iomsg=errmsg) hdum(1),hdum(2),hdum(3)
         if (abs(hdum(1))>1.d-9)nprojso(ilmax)=1
         if (abs(hdum(2))>1.d-9)nprojso(ilmax)=2
         if (abs(hdum(3))>1.d-9)nprojso(ilmax)=3
         if(nprojso(ilmax)>0)pspheads(ipsp)%pspso=2
       end if
       if (ilmax==3) then
         read (unt,*, err=10, iomsg=errmsg) hdum(1)
         if (abs(hdum(1))>1.d-9)nprojso(3)=1
         if(nprojso(3)>0)pspheads(ipsp)%pspso=2
       end if
     end do

   else if(pspcod==5)then

!    PHONEY pseudopotentials
!    read parameter for Hamman grid
     pspso=1
     read (unt,fmt=*,err=50,end=50) r1,al,pspso
     50 continue
     do ilmax=0,lmax
       read (unt,*, err=10, iomsg=errmsg) lang,e990,e999,nproj(ilmax)
       read (unt,*, err=10, iomsg=errmsg)
       if (ilmax>0.and.pspso/=1) then
         read (unt,*, err=10, iomsg=errmsg) lang,e990,e999,nprojso(ilmax)
         read (unt,*, err=10, iomsg=errmsg)
         pspheads(ipsp)%pspso=pspso
!        Meaning of pspso internally to ABINIT has been changed in v5.4
!        So : file must contain pspso 1 , but ABINIT will have pspso 0 .
         if(pspso==1)pspheads(ipsp)%pspso=0
       end if
     end do
     read (unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==6)then

!    FHI pseudopotentials
     read (unt, '(a3)') testxc
!    Note : prior to version 2.2, this 4th line started with  4--  ,
!    and no core-correction was available.
     if(testxc/='4--')then
       backspace(unt, err=10, iomsg=errmsg)
       read (unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     else
       fchrg=0.0_dp
     end if
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
!    XG020728 : Should take lloc into account ??
     do ilmax=0,lmax
       nproj(ilmax)=1
     end do

   else if(pspcod==7)then

!    PAW pseudopotentials
     test_paw=1;pspheads(ipsp)%pawheader%pawver=1
     read (unt,'(a80)', err=10, iomsg=errmsg) pspline 
     pspline=adjustl(pspline)
     if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
&     read(unit=pspline(4:80),fmt=*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%pawver
     if (pspheads(ipsp)%pawheader%pawver==1) then   ! Compatibility with Abinit v4.2.x
       read (unit=pspline,fmt=*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%basis_size,&
&       pspheads(ipsp)%pawheader%lmn_size
       ABI_ALLOCATE(orb,(pspheads(ipsp)%pawheader%basis_size))
       orb(:)=0
       read (unt,*, err=10, iomsg=errmsg) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
       read (unt,*, err=10, iomsg=errmsg)
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%rpaw
       pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%mesh_size
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%shape_type
       if (pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
     else
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%basis_size,&
&       pspheads(ipsp)%pawheader%lmn_size
       ABI_ALLOCATE(orb,(pspheads(ipsp)%pawheader%basis_size))
       orb(:)=0
       read (unt,*, err=10, iomsg=errmsg) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
       pspheads(ipsp)%pawheader%mesh_size=mmax
       read (unt,*, err=10, iomsg=errmsg) nmesh
       do ii=1,nmesh
         read(unt,*, err=10, iomsg=errmsg)
       end do
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%rpaw
       pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
       read (unt,'(a80)', err=10, iomsg=errmsg) pspline 
       pspline=adjustl(pspline); write(std_out,*) pspline
       read(unit=pspline,fmt=*) pspheads(ipsp)%pawheader%shape_type
       if (pspheads(ipsp)%pawheader%pawver==2.and.&
&       pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
       if (pspheads(ipsp)%pawheader%pawver>=3.and.pspheads(ipsp)%pawheader%shape_type==-1) then
         rr=zero;read(unit=pspline,fmt=*,err=20,end=20) ii,rr
         20       continue
         if (rr>=tol8) pspheads(ipsp)%pawheader%rshp=rr
       end if
     end if
     do ilmax=0,lmax
       do ii=1,pspheads(ipsp)%pawheader%basis_size
         if(orb(ii)==ilmax) nproj(ilmax)=nproj(ilmax)+1
       end do
     end do
     pspheads(ipsp)%pawheader%l_size=2*maxval(orb)+1
     pspheads(ipsp)%xccc=1  ! We suppose apriori that cc is used (but n1xccc is not used in PAW)
     ABI_DEALLOCATE(orb)

!    WVL+PAW case, need to define GTHradii
#if defined HAVE_BIGDFT
     if(pspheads(ipsp)%usewvl==1) then
!      Obtain the HGH parameters by default from BigDFT

       call atomic_info(int(pspheads(ipsp)%znuclpsp), int(pspheads(ipsp)%zionpsp), &
&       symbol = symbol, ehomo = ehomo, rcov = rcov, nsccode = iasctype)

!      I use the XC: Perdew, Burke & Ernzerhof  as default, since
!      other XC potentials may not be in the BigDFT table.
       ixc_=1 
       call psp_from_data(symbol, nzatom, nelpsp, npspcode_, ixc_, psppar, exists)
       if(.not. exists) then
         write(message,'(4a)')ch10,&
&         "Chemical element not found in BigDFT table",ch10,&
&         "Action: upgrade BigDFT table"
         MSG_BUG(message)
       end if
!      
!      pspheads(ipsp)%pawheader%rpaw/4.0d0
       pspheads(ipsp)%GTHradii(0)=psppar(0,0) !rloc
       pspheads(ipsp)%GTHradii(1)=psppar(1,0) !rrs
       pspheads(ipsp)%GTHradii(2)=psppar(2,0) !rrp
!      pspheads(ipsp)%GTHradii(1) = one / sqrt(abs(two * ehomo))
!      write(*,*)pspheads(ipsp)%GTHradii(:)
     end if
#endif

   else if(pspcod==8)then

!    DRH pseudopotentials
     read(unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
     read(unt,*, err=10, iomsg=errmsg) nproj(0:lmax)
     read(unt,*, err=10, iomsg=errmsg) extension_switch
     if(any(extension_switch == [2, 3])) then
       pspso=2
       read(unt,*,err=10,iomsg=errmsg) nprojso(1:lmax)
     else
       pspso=0
     end if
     pspheads(ipsp)%pspso=pspso

   else if(pspcod==9)then
! placeholder: nothing to do everything is read above

   else if(pspcod==10)then

!    HGH pseudopotentials, full h/k matrices
     read (unt,*,err=10,iomsg=errmsg) pspheads(ipsp)%GTHradii(0) !rloc
     read (unt,*,err=10,iomsg=errmsg) idum
     if(idum-1/=lmax) then
       MSG_ERROR("in inpspheads: nnonloc-1 /= lmax")
     end if
     do ilmax=0,lmax
       read (unt,*,err=10,iomsg=errmsg) &
       pspheads(ipsp)%GTHradii(ilmax + 1),nproj(ilmax),(hdum(idum),idum=1,nproj(ilmax))
       do idum=2,nproj(ilmax) !skip the rest of h_ij
         read (unt,*, err=10, iomsg=errmsg)
       end do
       if (ilmax==0) cycle
       nprojso(ilmax)=nproj(ilmax)
       if(nprojso(ilmax)>0)then
         pspheads(ipsp)%pspso=2
         do idum=1,nprojso(ilmax) !skip the rest of k_ij
           read (unt,*, err=10, iomsg=errmsg)
         end do
       end if
     end do

   else if (pspcod == 11.or.pspcod == 17) then
!    already done above
   else

     write(message, '(a,i0,4a)' )&
&     'The pseudopotential code (pspcod) read from file is ',pspcod,ch10,&
&     'This value is not allowed (should be between 1 and 10). ',ch10,&
&     'Action: use a correct pseudopotential file.'
     MSG_ERROR(message)
   end if ! pspcod=...

!  Store in pspheads
   if (pspcod /= 17) then
     pspheads(ipsp)%nproj(0:3)=nproj(0:3)
     pspheads(ipsp)%nprojso(1:3)=nprojso(1:3)
   end if
!  DEBUG
!   write(message,'(a,10i6)') 'nproj = ', pspheads(ipsp)%nproj(:)
!   call wrtout(std_out,message,'PERS')
!   write(message,'(a,10i6)') 'nprojso = ', pspheads(ipsp)%nprojso(:)
!   call wrtout(std_out,message,'PERS')
!  ENDDEBUG

   close(unt)

   ! Compute md5 checksum
   pspheads(ipsp)%md5_checksum = md5_sum_from_file(filnam(ipsp))
 end do ! ipsp=1,npsp

!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! Likely troubles with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do

 write(message,'(2a,i0,a,i0,a)')ch10,&
& ' inpspheads: deduce mpsang = ',mpsang,', n1xccc = ',n1xccc,'.'
 call wrtout(std_out,message,'PERS')

!Test: if one psp is PAW, all must be
 if (test_paw==1) then
   do ipsp=1,npsp
     if (all(pspheads(ipsp)%pspcod /= [7, 17])) then
       write(message, '(5a)' )&
&       'One pseudopotential is PAW (pspcod=7 or 17) !',ch10,&
&       'All pseudopotentials must be PAW (this is not the case here) !',ch10,&
&       'Action: use only PAW pseudopotential files.'
       MSG_ERROR(message)
     end if
   end do
 end if

 return

 ! Handle IO error
 10 continue
 MSG_ERROR(errmsg)

end subroutine inpspheads
!!***
