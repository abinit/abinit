!!****m* ABINIT/m_pspheads
!! NAME
!! m_pspheads
!!
!! FUNCTION
!!  Functions used to read the pseudopotential header of each psp file,
!!  in order to initialize pspheads(1:npsp).
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2022 ABINIT group (DCA, XG, GMR, FrD, AF, MT, FJ, MJV, MG, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_pspheads

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_atomdata
 use m_hash_md5
 use m_psxml2ab
#if defined HAVE_LIBPSML
 use m_psml
#endif
#if defined HAVE_BIGDFT
 use BigDFT_API, only: atomic_info, psp_from_data
#endif

 use defs_datatypes, only : pspheader_type
 use m_time,         only : timab
 use m_io_tools,     only : open_file
 use m_numeric_tools,only : simpson
 use m_fstrings,     only : basename, lstrip, sjoin, startswith, atoi, itoa, ftoa
 use m_pawpsp,       only : pawpsp_read_header_xml,pawpsp_read_pawheader
 use m_pawxmlps,     only : rdpawpsxml,rdpawpsxml_header, paw_setup_free,paw_setuploc
 use pseudo_types,    only : pseudo_upf, deallocate_pseudo_upf !, pseudo_config
 use read_upf_new_module, only : read_upf_new

 implicit none

 private
!!***

 public :: inpspheads      ! Initialize pspheads(1:npsp).
 public :: pspheads_comm   ! Communicate pspheads to all processors
 public :: pawpsxml2ab

 public :: upf2_jl2srso
 public :: upfxc2abi       ! UPF XcC to Abinit pspxc (DEPRECATED. Only used for UPF1)
 public :: upfdft_to_ixc   ! UPF2 dft to Abinit pspxc.

contains
!!***

!!****f* m_pspheads/inpspheads
!! NAME
!! inpspheads
!!
!! FUNCTION
!! Read the pseudopotential header of each psp file, in order to initialize pspheads(1:npsp).
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!
!! OUTPUT
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!  pseudopotential file headers, as well as the psp file names
!!  ecut_tmp(3,2,npsp)= possible ecut values as read in psp files
!!
!! SOURCE

subroutine inpspheads(filnam, npsp, pspheads, ecut_tmp)

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
 integer :: idum,ii,ilmax,ipsp,lang,lmax,mmax,mpsang,n1xccc,nmesh
 integer :: pspcod,pspso,test_paw,usexml,unt,useupf
 real(dp) :: al,e990,e999,fchrg,qchrg,r1,rchrg,rr
 character(len=3) :: testxc
 character(len=500) :: msg,errmsg
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
#if defined HAVE_LIBPSML
 character(len=3) :: atmsymb
 character(len=30) :: creator
#endif

!*************************************************************************

 test_paw=0

 do ipsp=1,npsp

   pspheads(ipsp)%filpsp=trim(filnam(ipsp))

   ! Check if the file is written in XML
   usexml = 0
   if (open_file(filnam(ipsp), msg, newunit=unt, form="formatted", status="old") /= 0) then
     ABI_ERROR(msg)
   end if

   rewind(unit=unt, err=10, iomsg=errmsg)
   read(unt, "(a)", err=10, iomsg=errmsg) testxml

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

   ! Check if pseudopotential file is a QE UPF2 file
   ! "<UPF version="2.0.1">
   useupf = 0
   if (testxml(1:4) == '<UPF') then
     ii = index(testxml, '"')
     if (ii /= 0) then
       useupf = atoi(testxml(ii+1:ii+1))
     else
       ABI_ERROR(sjoin("Cannot find version attribute in UPF2 file:", filnam(ipsp)))
     end if
   end if

   close(unit=unt, err=10, iomsg=errmsg)

   ! Check if pseudopotential file is a QE UPF1 file
   if (useupf == 0) then
     if (open_file(filnam(ipsp), msg, newunit=unt, form="formatted", status="old") /= 0) then
       ABI_ERROR(msg)
     end if

     rewind(unit=unt, err=10, iomsg=errmsg)
     read(unt,*,err=10,iomsg=errmsg) testxml ! just a string, no relation to xml.
     if (testxml(1:9)=='<PP_INFO>') then
       useupf = 1
     else
       useupf = 0
     end if
     close(unit=unt,err=10,iomsg=errmsg)
   end if

   ! Read the header of the pseudopotential file
   if (usexml /= 1 .and. useupf == 0) then
     ! Open the psp file and read a normal abinit style header
     if (open_file(filnam(ipsp), msg, newunit=unt, form='formatted', status='old') /= 0) then
       ABI_ERROR(msg)
     end if
     rewind (unit=unt, err=10, iomsg=errmsg)

     ! Read the three first lines
     read(unt, '(a)', err=10, iomsg=errmsg) pspheads(ipsp)%title
     read(unt,*, err=10, iomsg=errmsg)pspheads(ipsp)%znuclpsp,pspheads(ipsp)%zionpsp,pspheads(ipsp)%pspdat
     read(unt,*, err=10, iomsg=errmsg)pspheads(ipsp)%pspcod,pspheads(ipsp)%pspxc,pspheads(ipsp)%lmax,idum,mmax

     pspcod=pspheads(ipsp)%pspcod
     lmax=pspheads(ipsp)%lmax
     write(msg,'(a,f5.1,a,i4,a,i4)')'  read the values zionpsp=',pspheads(ipsp)%zionpsp,' , pspcod=',pspcod,' , lmax=',lmax
     call wrtout(std_out,msg,'PERS')

     nproj(0:3)=0 ; nprojso(1:3)=0

     pspheads(ipsp)%xccc=0
     pspheads(ipsp)%pspso=0

   else if (usexml==1 .and. test_paw==0) then
#if defined HAVE_LIBPSML
     write(msg,'(4a)')  &
       '- inpspheads : Reading pseudopotential header in XML form from ',ch10,&
&      '-   ',trim(filnam(ipsp))
     call wrtout([std_out, ab_out], msg)

     ! could pass pspheads(ipsp) directly and fill all of it in psxml2ab
     call psxml2abheader( filnam(ipsp), pspheads(ipsp), atmsymb, creator, 1)

     ! save some stuff locally for this ipsp
     pspcod = pspheads(ipsp)%pspcod
     lmax   = pspheads(ipsp)%lmax
     nproj = pspheads(ipsp)%nproj
     nprojso = pspheads(ipsp)%nprojso

#else
     write(msg, '(2a)') "XML norm-conserving pseudopotential has been input,", &
       " but abinit is not compiled with libPSML support. Reconfigure and recompile."
     ABI_ERROR(msg)
#endif

   else if(usexml==1.and.test_paw==1)then

     write(msg,'(4a)')  &
       '- inpspheads : Reading pseudopotential header in XML form from ',ch10,&
       '-   ',trim(filnam(ipsp))
     call wrtout([std_out, ab_out], msg)

     call pawpsxml2ab(filnam(ipsp),ecut_tmp(:,:,ipsp), pspheads(ipsp),1)
     pspcod=17; pspheads(ipsp)%pspcod=pspcod

   else if (useupf > 0) then
     pspheads(ipsp)%xccc  = n1xccc_default ! will be set to 0 if no nlcc

     if (useupf == 1) then
       pspheads(ipsp)%pspcod = 11
       call upf1_to_psphead(filnam(ipsp), pspheads(ipsp)%znuclpsp, pspheads(ipsp)%zionpsp, pspheads(ipsp)%pspxc, &
         pspheads(ipsp)%lmax, pspheads(ipsp)%xccc, nproj, nprojso)

       ! FIXME: generalize for SO pseudos
       pspheads(ipsp)%pspso = 0

     else
       ! UPF2 format
       pspheads(ipsp)%pspcod = 12
       call upf2_to_psphead(filnam(ipsp), pspheads(ipsp)%znuclpsp, pspheads(ipsp)%zionpsp, pspheads(ipsp)%pspxc, &
         pspheads(ipsp)%lmax, pspheads(ipsp)%xccc, nproj, nprojso)

       pspheads(ipsp)%pspso = merge(2, 0, any(nprojso > 0))
     end if

     pspcod = pspheads(ipsp)%pspcod
     lmax   = pspheads(ipsp)%lmax
   end if

   !write(std_out,*) "pspheads(ipsp)%znuclpsp", pspheads(ipsp)%znuclpsp
   !write(std_out,*) "pspheads(ipsp)%zionpsp",  pspheads(ipsp)%zionpsp
   !write(std_out,*) "pspheads(ipsp)%pspcod",   pspheads(ipsp)%pspcod
   !write(std_out,*) "pspheads(ipsp)%pspxc",    pspheads(ipsp)%pspxc
   !write(std_out,*) "pspheads(ipsp)%lmax",     pspheads(ipsp)%lmax

   ! Initialize nproj, nprojso, pspso, as well as xccc, for each type of psp
   pspheads(ipsp)%GTHradii = zero

   if(pspcod==1 .or. pspcod==4)then

     ! Teter format
     do ilmax=0,lmax
       read(unt,*, err=10, iomsg=errmsg) lang,e990,e999,nproj(ilmax)
       read(unt,*, err=10, iomsg=errmsg)
     end do
     read(unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==2)then

     ! GTH pseudopotentials
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(0) !rloc
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(1),hdum(1),hdum(2)
     if(abs(hdum(1))>1.d-9) nproj(0)=1
     if(abs(hdum(2))>1.d-9) nproj(0)=2
     read(unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%GTHradii(2),hdum(3)
     if(abs(hdum(3))>1.d-9) nproj(1)=1

   else if(pspcod==3)then

     ! HGH pseudopotentials
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

     ! PHONEY pseudopotentials
     ! read parameter for Hamman grid
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
         ! Meaning of pspso internally to ABINIT has been changed in v5.4
         ! So file must contain pspso 1, but ABINIT will have pspso 0.
         if(pspso==1)pspheads(ipsp)%pspso=0
       end if
     end do
     read (unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default

   else if(pspcod==6)then

     ! FHI pseudopotentials
     read (unt, '(a3)') testxc
     ! Note: prior to version 2.2, this 4th line started with  4--  ,
     ! and no core-correction was available.
     if(testxc/='4--')then
       backspace(unt, err=10, iomsg=errmsg)
       read (unt,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
     else
       fchrg=0.0_dp
     end if
     if (fchrg>1.d-15) pspheads(ipsp)%xccc=n1xccc_default
     ! XG020728 : Should take lloc into account ??
     do ilmax=0,lmax
       nproj(ilmax)=1
     end do

   else if(pspcod==7)then

     ! PAW pseudopotentials
     test_paw=1;pspheads(ipsp)%pawheader%pawver=1
     read (unt,'(a80)', err=10, iomsg=errmsg) pspline
     pspline=adjustl(pspline)
     if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
&      read(unit=pspline(4:80),fmt=*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%pawver
     if (pspheads(ipsp)%pawheader%pawver==1) then   ! Compatibility with Abinit v4.2.x
       read (unit=pspline,fmt=*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%basis_size,&
         pspheads(ipsp)%pawheader%lmn_size
       ABI_MALLOC(orb,(pspheads(ipsp)%pawheader%basis_size))
       orb(:)=0
       read (unt,*, err=10, iomsg=errmsg) (orb(ii), ii=1,pspheads(ipsp)%pawheader%basis_size)
       read (unt,*, err=10, iomsg=errmsg)
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%rpaw
       pspheads(ipsp)%pawheader%rshp=pspheads(ipsp)%pawheader%rpaw
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%mesh_size
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%shape_type
       if (pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
     else
       read (unt,*, err=10, iomsg=errmsg) pspheads(ipsp)%pawheader%basis_size,pspheads(ipsp)%pawheader%lmn_size
       ABI_MALLOC(orb,(pspheads(ipsp)%pawheader%basis_size))
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
           pspheads(ipsp)%pawheader%shape_type==3) pspheads(ipsp)%pawheader%shape_type=-1
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
     ABI_FREE(orb)

#if defined HAVE_BIGDFT
     ! WVL+PAW case, need to define GTHradii
     if(pspheads(ipsp)%usewvl==1) then
       ! Obtain the HGH parameters by default from BigDFT

       call atomic_info(int(pspheads(ipsp)%znuclpsp), int(pspheads(ipsp)%zionpsp), &
         symbol = symbol, ehomo = ehomo, rcov = rcov, nsccode = iasctype)

       ! I use the XC: Perdew, Burke & Ernzerhof  as default, since
       ! other XC potentials may not be in the BigDFT table.
       ixc_=1
       call psp_from_data(symbol, nzatom, nelpsp, npspcode_, ixc_, psppar, exists)
       if(.not. exists) then
         write(msg,'(4a)')ch10,&
          "Chemical element not found in BigDFT table",ch10,&
          "Action: upgrade BigDFT table"
         ABI_BUG(msg)
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

     ! DRH pseudopotentials
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

     ! HGH pseudopotentials, full h/k matrices
     read (unt,*,err=10,iomsg=errmsg) pspheads(ipsp)%GTHradii(0) !rloc
     read (unt,*,err=10,iomsg=errmsg) idum
     if(idum-1/=lmax) then
       ABI_ERROR("in inpspheads: nnonloc-1 /= lmax")
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

   else if (any(pspcod == [11, 12, 17])) then
     ! already done above

   else
     write(msg, '(a,i0,4a)' )&
       'The pseudopotential code (pspcod) read from file is ',pspcod,ch10,&
       'This value is not allowed.',ch10,&
       'Action: use a correct pseudopotential file.'
     ABI_ERROR(msg)
   end if ! pspcod

   ! Store in pspheads
   if (pspcod /= 17) then
     pspheads(ipsp)%nproj(0:3)=nproj(0:3)
     pspheads(ipsp)%nprojso(1:3)=nprojso(1:3)
   end if
   !write(std_out,'(a,*(i0,1x))') 'nproj = ', pspheads(ipsp)%nproj(:)
   !write(std_out,'(a,*(i0,1x))') 'nprojso = ', pspheads(ipsp)%nprojso(:)

   close(unt)

   ! Compute md5 checksum
   pspheads(ipsp)%md5_checksum = md5_sum_from_file(filnam(ipsp))
 end do ! ipsp=1,npsp

 ! Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do

 write(msg,'(2a,i0,a,i0,a)')ch10,' inpspheads: deduce mpsang = ',mpsang,', n1xccc = ',n1xccc,'.'
 call wrtout(std_out,msg,'PERS')

 ! Test: if one psp is PAW, all must be
 if (test_paw==1) then
   do ipsp=1,npsp
     if (all(pspheads(ipsp)%pspcod /= [7, 17])) then
       write(msg, '(5a)' )&
        'One pseudopotential is PAW (pspcod=7 or 17) !',ch10,&
        'All pseudopotentials must be PAW (this is not the case here) !',ch10,&
        'Action: use only PAW pseudopotential files.'
       ABI_ERROR(msg)
     end if
   end do
 end if

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine inpspheads
!!***

!!****f* ABINIT/pspheads_comm
!! NAME
!! pspheads_comm
!!
!! FUNCTION
!! Communicate pspheads to all processors
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!  test_paw=0 if no PAW, 1 if PAW
!!
!! SIDE EFFECTS
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names. On one processor at input,
!!   on all processors at output
!!
!! SOURCE

subroutine pspheads_comm(npsp,pspheads,test_paw)

!Arguments ------------------------------------
 integer,intent(in) :: npsp
 integer,intent(inout) :: test_paw
 type(pspheader_type),intent(inout) :: pspheads(npsp)

!Local variables-------------------------------
#if defined HAVE_MPI
!scalars
 integer,parameter :: master=0
 integer :: ierr,comm
!arrays
 integer,allocatable :: list_int(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: list_dpr(:)
 character(len=fnlen),allocatable :: list_char(:)
#endif

!*************************************************************************

#if defined HAVE_MPI
 call timab(48,1,tsec)

 comm = xmpi_world

 ! Broadcast the characters (file names and titles)
 ABI_MALLOC(list_char,(3*npsp))
 list_char(1:npsp)=pspheads(1:npsp)%filpsp
 list_char(npsp+1:2*npsp)=pspheads(1:npsp)%title
 list_char(2*npsp+1:3*npsp)=pspheads(1:npsp)%md5_checksum

 call xmpi_bcast(list_char,master,comm,ierr)

 pspheads(1:npsp)%filpsp=list_char(1:npsp)
 pspheads(1:npsp)%title=list_char(npsp+1:2*npsp)
 pspheads(1:npsp)%md5_checksum=list_char(2*npsp+1:3*npsp)(1:md5_slen)
 ABI_FREE(list_char)

 ! Brodcast the integers
 ABI_MALLOC(list_int,(1+13*npsp))
 list_int(1        :   npsp) = pspheads(1:npsp)%nproj(0)
 list_int(1+   npsp: 2*npsp) = pspheads(1:npsp)%nproj(1)
 list_int(1+ 2*npsp: 3*npsp) = pspheads(1:npsp)%nproj(2)
 list_int(1+ 3*npsp: 4*npsp) = pspheads(1:npsp)%nproj(3)
 list_int(1+ 4*npsp: 5*npsp) = pspheads(1:npsp)%lmax
 list_int(1+ 5*npsp: 6*npsp) = pspheads(1:npsp)%xccc
 list_int(1+ 6*npsp: 7*npsp) = pspheads(1:npsp)%pspxc
 list_int(1+ 7*npsp: 8*npsp) = pspheads(1:npsp)%pspdat
 list_int(1+ 8*npsp: 9*npsp) = pspheads(1:npsp)%pspcod
 list_int(1+ 9*npsp:10*npsp) = pspheads(1:npsp)%pspso
 list_int(1+10*npsp:11*npsp) = pspheads(1:npsp)%nprojso(1)
 list_int(1+11*npsp:12*npsp) = pspheads(1:npsp)%nprojso(2)
 list_int(1+12*npsp:13*npsp) = pspheads(1:npsp)%nprojso(3)
 list_int(1+13*npsp)         = test_paw

 call xmpi_bcast(list_int,master,comm,ierr)

 pspheads(1:npsp)%nproj(0)   = list_int(1        :   npsp)
 pspheads(1:npsp)%nproj(1)   = list_int(1+   npsp: 2*npsp)
 pspheads(1:npsp)%nproj(2)   = list_int(1+ 2*npsp: 3*npsp)
 pspheads(1:npsp)%nproj(3)   = list_int(1+ 3*npsp: 4*npsp)
 pspheads(1:npsp)%lmax       = list_int(1+ 4*npsp: 5*npsp)
 pspheads(1:npsp)%xccc       = list_int(1+ 5*npsp: 6*npsp)
 pspheads(1:npsp)%pspxc      = list_int(1+ 6*npsp: 7*npsp)
 pspheads(1:npsp)%pspdat     = list_int(1+ 7*npsp: 8*npsp)
 pspheads(1:npsp)%pspcod     = list_int(1+ 8*npsp: 9*npsp)
 pspheads(1:npsp)%pspso      = list_int(1+ 9*npsp:10*npsp)
 pspheads(1:npsp)%nprojso(1) = list_int(1+10*npsp:11*npsp)
 pspheads(1:npsp)%nprojso(2) = list_int(1+11*npsp:12*npsp)
 pspheads(1:npsp)%nprojso(3) = list_int(1+12*npsp:13*npsp)
 test_paw                    = list_int(1+13*npsp)
 ABI_FREE(list_int)

 ! Unbeliveable, this cannot be sent with the others, for woopy
 ABI_MALLOC(list_int,(npsp))
 list_int(1:npsp) = pspheads(1:npsp)%usewvl
 call xmpi_bcast(list_int,master,comm,ierr)
 pspheads(1:npsp)%usewvl     = list_int(1:npsp)
 ABI_FREE(list_int)

 ! Broadcast zionpsp and znuclpsp
 ABI_MALLOC(list_dpr,(7*npsp))
 list_dpr(1       :  npsp) = pspheads(1:npsp)%zionpsp
 list_dpr(1+  npsp:2*npsp) = pspheads(1:npsp)%znuclpsp
 list_dpr(1+2*npsp:3*npsp) = pspheads(1:npsp)%GTHradii(0)
 list_dpr(1+3*npsp:4*npsp) = pspheads(1:npsp)%GTHradii(1)
 list_dpr(1+4*npsp:5*npsp) = pspheads(1:npsp)%GTHradii(2)
 list_dpr(1+5*npsp:6*npsp) = pspheads(1:npsp)%GTHradii(3)
 list_dpr(1+6*npsp:7*npsp) = pspheads(1:npsp)%GTHradii(4)

 call xmpi_bcast(list_dpr,master,comm,ierr)

 pspheads(1:npsp)%zionpsp     = list_dpr(1       :  npsp)
 pspheads(1:npsp)%znuclpsp    = list_dpr(1+  npsp:2*npsp)
 pspheads(1:npsp)%GTHradii(0) = list_dpr(1+2*npsp:3*npsp)
 pspheads(1:npsp)%GTHradii(1) = list_dpr(1+3*npsp:4*npsp)
 pspheads(1:npsp)%GTHradii(2) = list_dpr(1+4*npsp:5*npsp)
 pspheads(1:npsp)%GTHradii(3) = list_dpr(1+5*npsp:6*npsp)
 pspheads(1:npsp)%GTHradii(4) = list_dpr(1+6*npsp:7*npsp)
 ABI_FREE(list_dpr)

 ! Broadcast additional integers for PAW psps (testpaw was sent, previously)
 if (test_paw==1) then
   ABI_MALLOC(list_int,(6*npsp))
   list_int(1       :  npsp)=pspheads(1:npsp)%pawheader%basis_size
   list_int(1+  npsp:2*npsp)=pspheads(1:npsp)%pawheader%l_size
   list_int(1+2*npsp:3*npsp)=pspheads(1:npsp)%pawheader%lmn_size
   list_int(1+3*npsp:4*npsp)=pspheads(1:npsp)%pawheader%mesh_size
   list_int(1+4*npsp:5*npsp)=pspheads(1:npsp)%pawheader%pawver
   list_int(1+5*npsp:6*npsp)=pspheads(1:npsp)%pawheader%shape_type

   call xmpi_bcast(list_int,master,comm,ierr)

   pspheads(1:npsp)%pawheader%basis_size=list_int(1       :  npsp)
   pspheads(1:npsp)%pawheader%l_size    =list_int(1+  npsp:2*npsp)
   pspheads(1:npsp)%pawheader%lmn_size  =list_int(1+2*npsp:3*npsp)
   pspheads(1:npsp)%pawheader%mesh_size =list_int(1+3*npsp:4*npsp)
   pspheads(1:npsp)%pawheader%pawver    =list_int(1+4*npsp:5*npsp)
   pspheads(1:npsp)%pawheader%shape_type=list_int(1+5*npsp:6*npsp)
   ABI_FREE(list_int)

   ! broadcast rpaw values
   ABI_MALLOC(list_dpr,(2*npsp))

   list_dpr(1       :  npsp) = pspheads(1:npsp)%pawheader%rpaw
   list_dpr(1+1*npsp:2*npsp) = pspheads(1:npsp)%pawheader%rshp

   call xmpi_bcast(list_dpr,master,comm,ierr)

   pspheads(1:npsp)%pawheader%rpaw = list_dpr(1       :  npsp)
   pspheads(1:npsp)%pawheader%rshp = list_dpr(1+  npsp:2*npsp)

   ABI_FREE(list_dpr)
 end if

 call timab(48,2,tsec)

#else
 ! Code to use unused dummy arguments
 if(pspheads(1)%lmax == -10) pspheads(1)%lmax=-10
 if(test_paw == -1) test_paw = -1
#endif

end subroutine pspheads_comm
!!***

!!****f* m_pspheads/pawpsxml2ab
!! NAME
!! pawpsxml2ab
!!
!! FUNCTION
!!  From a XML format pseudopotential file which has already been read in,
!!  convert to abinit internal datastructures.
!!
!! INPUTS
!!  ecut_tmp(3,2)= possible ecut values as read in psp files
!!  filenam= input file name (atomicdata XML)
!!  option= 1 if header only is read; 0 if the whole data are read
!!
!! OUTPUT
!! pspheads data structure is filled
!!
!! SOURCE

subroutine pawpsxml2ab(filnam, ecut_tmp, pspheads, option)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: option
 character(len=fnlen), intent(in) :: filnam
 type(pspheader_type),intent(inout) :: pspheads !vz_i
!arrays
 real(dp),intent(inout) :: ecut_tmp(3,2)

!Local variables-------------------------------
 integer :: ii,il,lloc,lmax,pspcod,pspxc
 real(dp) :: r2well,zionpsp,znuclpsp
! character(len=100) :: xclibxc, msg

! *********************************************************************

 if (option==1) then
   call rdpawpsxml_header(ecut_tmp,filnam,paw_setuploc)
   paw_setuploc%idgrid= paw_setuploc%radial_grid(1)%id
 else
   call rdpawpsxml(filnam,paw_setuploc)
 end if

 call pawpsp_read_header_xml(lloc,lmax,pspcod, pspxc,paw_setuploc,r2well,zionpsp,znuclpsp)

 pspheads%lmax=lmax
 pspheads%pspxc=pspxc
 pspheads%zionpsp=zionpsp
 pspheads%znuclpsp=znuclpsp

 call pawpsp_read_pawheader(pspheads%pawheader%basis_size,&
   pspheads%lmax,pspheads%pawheader%lmn_size,&
   pspheads%pawheader%l_size,pspheads%pawheader%mesh_size,&
   pspheads%pawheader%pawver,paw_setuploc,pspheads%pawheader%rpaw,&
   pspheads%pawheader%rshp,pspheads%pawheader%shape_type)

 pspheads%nproj=0
 do il=0,pspheads%lmax
   do ii=1,pspheads%pawheader%basis_size
     if(paw_setuploc%valence_states%state(ii)%ll==il) pspheads%nproj(il)=pspheads%nproj(il)+1
   end do
 end do

 pspheads%nprojso=0
 pspheads%pspdat=27061961
 pspheads%pspso=1
 pspheads%xccc=1
 pspheads%title=paw_setuploc%atom%symbol

 if (option==1) call paw_setup_free(paw_setuploc)

end subroutine pawpsxml2ab
!!***

!!****f* m_pspheads/upf1_to_psphead
!! NAME
!! upf1_to_psphead
!!
!! FUNCTION
!!  This routine wraps a call to a PWSCF module, which reads in
!!  a UPF1 (PWSCF / Espresso) format pseudopotential, then transfers
!!  data for the HEADER of abinit psps only!
!!
!! INPUTS
!!  filpsp = name of file with UPF1 data
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!  lmax_ = maximal angular momentum
!!  znucl = charge of species nucleus
!!  zion = valence charge
!!  n1xccc = default number of points. Set to 0 if no nlcc is present
!!  nproj_l= number of projectors for each channel
!!  nprojso_l= number of projectors for each channel for SO correction projectors
!!
!! SOURCE

subroutine upf1_to_psphead(filpsp, znucl, zion, pspxc, lmax_, n1xccc, nproj_l, nprojso_l)

 use m_read_upf_pwscf,  only : read_pseudo
 use pseudo_pwscf ! pwscf module with all data explicit!

!Arguments -------------------------------
 character(len=fnlen), intent(in) :: filpsp
 integer,intent(inout) :: n1xccc
 integer,intent(out) :: pspxc, lmax_
 real(dp),intent(out) :: znucl, zion
!arrays
 integer,intent(out) :: nproj_l(0:3)
 integer,intent(out) :: nprojso_l(1:3)

!Local variables -------------------------
 integer :: iproj, ll, iunit
 character(len=500) :: msg
 type(atomdata_t) :: atom

! *********************************************************************

 if (open_file(filpsp, msg, newunit=iunit, status='old',form='formatted') /= 0) then
   ABI_ERROR(msg)
 end if

 ! read in psp data to static data in pseudo module, for ipsx == 1
 call read_pseudo(1,iunit)
 close (iunit)

 ! copy over to abinit internal arrays and vars
 ! FIXME: The API is broken. It does not recognize PBEsol
 ! should use upfdft_to_ixc
 call upfxc2abi(dft(1), pspxc)
 lmax_ = lmax(1)
 call atomdata_from_symbol(atom,psd(1))
 znucl = atom%znucl
 zion = zp(1)

 nproj_l = 0
 do iproj = 1, nbeta(1)
   ll = lll(iproj,1)
   nproj_l(ll) = nproj_l(ll) + 1
 end do

 nprojso_l = 0 !FIXME deal with so
 !do iproj = 1, nbeta(1)
 !nprojso_l(ll+1) = nprojso_l(ll+1) + 1
 !end do

 if (.not. nlcc(1)) n1xccc = 0

end subroutine upf1_to_psphead
!!***

!!****f* m_pspheads/upf2_to_psphead
!! NAME
!! upf2_to_psphead
!!
!! FUNCTION
!!  This routine wraps a call to a PWSCF module, which reads in
!!  a UPF2 (PWSCF / Espresso) format pseudopotential, then transfers
!!  data for the HEADER of abinit psps only!
!!
!! INPUTS
!!  filpsp = name of file with UPF1 data
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!  lmax_ = maximal angular momentum
!!  znucl = charge of species nucleus
!!  zion = valence charge
!!  n1xccc = default number of points. Set to 0 if no nlcc is present
!!  nproj_l= number of projectors for each channel
!!  nprojso_l= number of projectors for each channel for SO correction projectors
!!
!! SOURCE

subroutine upf2_to_psphead(filpsp, znucl, zion, pspxc, lmax, n1xccc, nproj_l, nprojso_l)

!Arguments -------------------------------
 character(len=fnlen), intent(in) :: filpsp
 integer,intent(inout) :: n1xccc
 integer,intent(out) :: pspxc, lmax
 real(dp),intent(out) :: znucl, zion
!arrays
 integer,intent(out) :: nproj_l(0:3), nprojso_l(1:3)

!Local variables -------------------------
 integer :: ierr , iprj, ll, mmax, irad
 real(dp) :: amesh, damesh
 character(len=500) :: msg
 logical :: linear_mesh
 type(pseudo_upf) :: upf
 type(atomdata_t) :: atom
! arrays
 real(dp),allocatable :: vsr(:,:,:), esr(:,:), vso(:,:,:), eso(:,:)

! *********************************************************************

 ! See also https://github.com/QEF/qeschemas/blob/master/UPF/qe_pp-0.99.xsd
 call read_upf_new(filpsp, upf, ierr)
 ABI_CHECK(ierr == 0, sjoin("read_upf_new returned ierr:", itoa(ierr)))

 call atomdata_from_symbol(atom, upf%psd)
 znucl = atom%znucl
 zion = upf%zp
 lmax = upf%lmax
 mmax = upf%mesh

 ! Consistency check
 ABI_CHECK(upf%typ == "NC", sjoin("Only NC pseudos in UPF2 format are supported while type is:", upf%typ))
 ABI_CHECK(upfdft_to_ixc(upf%dft, pspxc, msg) == 0, msg)
 if (.not. upf%nlcc) n1xccc = 0

 ! Check that rad grid is linear starting at zero
 linear_mesh = .True.
 amesh = upf%r(2) - upf%r(1); damesh = zero
 do irad=2,mmax-1
   damesh = max(damesh, abs(upf%r(irad)+amesh-upf%r(irad+1)))
 end do
 linear_mesh = damesh < tol8

 if (.not. linear_mesh .or. abs(upf%r(1)) > tol16) then
   write(msg,'(3a)')&
   'Assuming pseudized valence charge given on linear radial mesh starting at zero.',ch10,&
   'Action: check your pseudopotential file.'
   ABI_ERROR(msg)
 end if

 nproj_l = 0; nprojso_l = 0

 if (.not. upf%has_so) then
   ! Scalar case
   do iprj=1,upf%nbeta
     ll = upf%lll(iprj)
     nproj_l(ll) = nproj_l(ll) + 1
   end do

 else
   ! Pseudo in j = l + s representation.
   call upf2_jl2srso(upf, nproj_l, nprojso_l, vsr, esr, vso, eso)

   ABI_FREE(vsr)
   ABI_FREE(esr)
   ABI_FREE(vso)
   ABI_FREE(eso)
 end if

 call deallocate_pseudo_upf(upf)

end subroutine upf2_to_psphead
!!***

!!****f* m_pspheads/upf2_jl2srso
!! NAME
!! upf2_jl2srso
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine upf2_jl2srso(upf, nproj_l, nprojso_l, vsr, esr, vso, eso)

!Arguments -------------------------------
 type(pseudo_upf),intent(in) :: upf
!arrays
 integer,intent(out) :: nproj_l(0:3), nprojso_l(1:3)
 real(dp),allocatable,intent(out) :: vsr(:,:,:), esr(:,:), vso(:,:,:), eso(:,:)

!Local variables -------------------------
 integer :: iprj, ii, ll, l1, il, ik, lmax, mmax, mxprj
 real(dp) :: jtot, eprmin !eps_srso,
 !character(len=500) :: msg
! arrays
integer :: irc6(6),nproj6(6), done_ilk(6,2)
 real(dp),allocatable :: vkb(:,:,:,:), evkb(:,:,:)

! *********************************************************************

 lmax = upf%lmax; mmax = upf%mesh
 nproj_l = 0; nprojso_l = 0

 ! Pseudo in j = l + s representation.
 irc6 = zero; nproj6 = zero
 do iprj=1,upf%nbeta
   ll = upf%lll(iprj)
   nproj6(ll+1) = nproj6(ll+1) + 1
   !irc6(ll+1) = max(upf%kbeta(iprj), irc6(ll+1))
   irc6(ll+1) = mmax
 end do

 ! Divide by two for l > 0 as this is sr_so_r expects.
 nproj6(2:) = nproj6(2:) / 2
 mxprj = maxval(nproj6)

 ABI_MALLOC(vkb, (mmax,mxprj,4,2))
 ABI_MALLOC(evkb, (mxprj,4,2))
 ABI_MALLOC(vsr, (mmax,2*mxprj,4))
 ABI_MALLOC(esr, (2*mxprj,4))
 ABI_MALLOC(vso, (mmax,2*mxprj,4))
 ABI_MALLOC(eso, (2*mxprj,4))

 done_ilk = 0
 do iprj=1,upf%nbeta
   jtot = upf%jjj(iprj)
   ll = upf%lll(iprj)
   il = ll + 1
   if (ll == 0) then
     ik = 1
   else
     ! l+1/2 --> ik 1, l-1/2 --> ik 2
     if (abs(jtot - (ll + half)) < tol6) then
       ik = 1
     else if (abs(jtot - (ll - half)) < tol6) then
       ik = 2
     else
       ABI_ERROR(sjoin("Cannot detect ik index from jtot:", ftoa(jtot)))
     end if
   end if

   done_ilk(il, ik) = done_ilk(il, ik) + 1
   ii = done_ilk(il, ik)
   evkb(ii,il,ik) = upf%dion(iprj,iprj) * half  ! convert from Rydberg to Ha
   vkb(:,ii,il,ik) = upf%beta(:,iprj)
 end do

 call sr_so_r(lmax, irc6, nproj6, upf%r, mmax, mxprj, evkb, vkb, vsr, esr, vso, eso)

 ! MG: This is done in oncvpsp 3.3 but not in oncvpsp4
 ! drop sr, so orthonormal projectors with neglibible coefficients
 ! modify cutoff if desired

 eprmin=2.0d-5
 write(std_out,'(/a,1p,e10.2,a)') 'Orthonormal projectors with coefficients <', &
     eprmin,' Ha will be dropped'

 do l1=1,lmax+1
  if(abs(esr(3,l1))<eprmin) esr(3,l1)=0.0d0
  if(abs(esr(4,l1))<eprmin) esr(4,l1)=0.0d0
  if(abs(eso(3,l1))<eprmin) eso(3,l1)=0.0d0
  if(abs(eso(4,l1))<eprmin) eso(4,l1)=0.0d0
 end do

#if 0
 ! MG: This is done in oncvpsp 4 but not in oncvpsp 3.3
 ! set smallest components to zero (following the approach used in oncvpsp)
 eps_srso=1.0d-3
 do l1=1,lmax+1
   if (nproj6(l1) > 0) then
     do iprj=2,2*nproj6(l1)
       if (abs(esr(iprj,l1)) < eps_srso*abs(esr(1,l1))) esr(iprj,l1) = 0.0d0
       if (l1 == 1) cycle
       if (abs(eso(iprj,l1)) < eps_srso*abs(eso(1,l1))) eso(iprj,l1) = 0.0d0
     end do
   end if
 end do
#endif

 ! set up projector number for sr_so calculations based on non-zero coefficients
 ! note that energies and projectors have been sorted sr_so_r
 ! so the relevant projectors are packed in the first positions.
 do l1=1,lmax+1
   ll = l1 - 1
   do ii=1,2*nproj6(l1)
    if (abs(esr(ii,l1)) > 0.0d0) nproj_l(ll) = nproj_l(ll) + 1
    if (abs(eso(ii,l1)) > 0.0d0) nprojso_l(ll) = nprojso_l(ll) + 1
   end do
   !write(std_out, '(a,3(i0,1x))')' ll, nproj_l, nprojso_l',ll, nproj_l(ll), nprojso_l(ll)
 end do

 ABI_FREE(vkb)
 ABI_FREE(evkb)

end subroutine upf2_jl2srso
!!***

!!****f* m_pspheads/upfxc2abi
!! NAME
!! upfxc2abi
!!
!! FUNCTION
!!  This routine wraps a call to an OCTOPUS module, which reformats
!!  a UPF (PWSCF / Espresso) string describing XC functionals,
!!  and returns the abinit internal code pspxc
!!
!! INPUTS
!!  dft = string with x/c functionals from PWSCF format
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!
!! NOTES
!!   FIXME: extend to more functionals with libxc
!!   Could be included in separate module, eg read_upf_pwscf or funct_pwscf
!!   Left without defs_basis or calls to abinit routines ON PURPOSE
!!
!! SOURCE

subroutine upfxc2abi(dft, pspxc)

 use funct_pwscf  ! pwscf module for naming xc functionals

!Arguments -------------------------------
 character(len=*), intent(in) :: dft
 integer, intent(out) :: pspxc

!Local variables -------------------------
 integer :: iexch,icorr,igcx,igcc
 integer :: totalindex, offset

! *********************************************************************
!extract from char*20 :: dft(:)
!###  The following has been copied from pwscf src/Modules/upf_to_internal.f90:
!workaround for rrkj format - it contains the indices, not the name
 if ( dft(1:6)=='INDEX:') then
   read( dft(7:10), '(4i1)') iexch,icorr,igcx,igcc
   call set_dft_from_indices(iexch,icorr,igcx,igcc)
 else
   call set_dft_from_name( dft )
   iexch = get_iexch()
   icorr = get_icorr()
   igcx = get_igcx()
   igcc = get_igcc()
 end if

 !reset dft string to avoid stray spaces
 call set_dft_from_indices(iexch,icorr,igcx,igcc)
 write(std_out,'(a)') ' upf2abinit: XC string from pseudopotential is :'
 write(std_out,'(3a)') '>', dft, '<'
 ABI_WARNING("upfxc2abi is not guaranteed to return the right ixc from QE XC string e.g. PBEsol. Please crosscheck!")

 offset = 100
 totalindex = offset*offset*offset*iexch + offset*offset*icorr + offset*igcx + igcc
 select case (totalindex)
 case (00000000)  !(" NOX  NOC NOGX NOGC") ! no xc
   pspxc = 0
 case (01010000)  !(" SLA   PZ NOGX NOGC") ! slater exchange + Perdew Zunger
   pspxc = 2
 case (01050000)  !(" SLA  WIG NOGX NOGC") ! slater exchange + Wigner corr
   pspxc = 4
 case (01060000)  !(" SLA   HL NOGX NOGC") ! Hedin + Lundqvist
   pspxc = 5
 case (02000000)  !(" SL1  NOC NOGX NOGC") ! full slater exchange
   pspxc = 6
 case (01040000)  !(" SLA   PW NOGX NOGC") ! slater exchange + Perdew Wang
   pspxc = 7
 case (01000000)  !(" SLA  NOC NOGX NOGC") ! Perdew Wang + no corr
   pspxc = 8
 case (01040304)  !(" SLA   PW  PBX  PBC") ! LDA + PBE GGA
   pspxc = 11 ! PBE
 case (01000300)  !(" SLA  NOC  PBX NOGC") ! exchange part of PBE GGA
   pspxc = 12
 case (01040404)  !(" SLA   PW  RPB  PBC") ! rev PBE
   pspxc = 14
 case (00000505)  !(" NOX  NOC HTCH HTCH") ! HTCH 120
   pspxc = 17
 case (01030103)  !(" SLA  LYP  B88 BLYP") ! BLYP
   pspxc = -106131
 case (01040101)  !(" SLA   PW  B88  P86") ! BP86
   pspxc = -106132
 case (00030603)  !(" NOX  LYP OPTX BLYP") ! OLYP
   pspxc = -110131
!    FIXME: important cases left to be patched with libxc:
!    vosko wilkins nusair
!    ortiz ballone
!    pbe0
!    Gunnarson-Lunqvist
!    make general approach: check gradient parts first, then lda.
!    event. check if they are consistent.
 case default
   ABI_ERROR('upf2abinit: XC functional not recognized')
 end select

end subroutine upfxc2abi
!!***

!!****f* m_pspheads/updft_to_ixc
!! NAME
!! updft_to_ixc
!!
!! FUNCTION
!!  Returns the abinit internal `ixc` from `dft` string with XC functional in QE format.
!!
!! SOURCE

integer function upfdft_to_ixc(dft, ixc, msg) result(ierr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: dft
 character(len=*),intent(out) :: msg
 integer,intent(out) :: ixc

!*************************************************************************

 ! This list taken from oncvpsp/src/upfout.f90
 ! It should be OK as long as the UPF2 NC pseudos are generated with oncvpsp
 ! but it does not cover all QE possibilities.
 ierr = 0; msg = ""
 select case (dft)
 case ("PZ")
   ixc = -001009
 case ("PBE")
   ixc = -101130 !; ixc = 11
 case ("PW91")
   ixc = -109134
 case ("PBESOL")
   ixc = -116133
 case ("REVPBE")
   ixc = -102130
 case ("BP")
   ixc = -106132
 case ("BLYP")
   ixc = -106131
 case ("WC")
   ixc = -118130
 case default
   ierr = 1
   write(msg, "(4a)") &
     "Cannot find ABINIT ixc value corresponding to QE dft:", trim(dft), ch10, &
     "Please update mapping in m_pspheads/upfdft_to_ixc."
 end select

end function upfdft_to_ixc
!!***
!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
 subroutine sr_so_r(lmax,irc,nproj,rr,mmax,mxprj,evkb,vkb, &
&                   vsr,esr,vso,eso)

! reformulates non-local potentials based on j = l +/- 1/2 to scalar-
! relativistic and L dot S projectors
! uses relationship <L dot S> = (J^2 - L^2 - S^2)/2
! so L dot S = +/- l/2 for j = l +/- 1/2

!lmax  maximum angular momentum
!irc  core radii indices
!nproj  number of projectors for each l
!rr  log radial grid
!mmax  size of radial grid
!mmax  dimension of log grid
!mxprj  dimension of number of projectors
!vkb  vkb projectors
!evkb  coefficients of BKB projectors
!vsr  normalized scalar projectors
!esr  energy  coefficients of vscal
!vso  normalized spin-orbig projectors
!esol  energy  coefficients of vso

 !implicit none
 !integer, parameter :: dp=kind(1.0d0)

!Input variables
 integer,intent(in) :: lmax,mmax,mxprj
 integer,intent(in) :: irc(6),nproj(6)
 real(dp),intent(in) :: rr(mmax),vkb(mmax,mxprj,4,2),evkb(mxprj,4,2)

!Output variables
 real(dp),intent(out) :: vsr(mmax,2*mxprj,4),esr(2*mxprj,4)
 real(dp),intent(out) :: vso(mmax,2*mxprj,4),eso(2*mxprj,4)

!Local variables
 integer :: ii,jj,kk,ik1,ik2,ip1,ip2,ipk,ll,l1,info,nn
 real(dp) :: amesh
 real(dp) :: apk,tt
 real(dp) :: sovl(2*mxprj,2*mxprj),sovlev(2*mxprj),ascl(2*mxprj,2*mxprj),aso(2*mxprj,2*mxprj)
 real(dp) :: asclst(2*mxprj,2*mxprj),wsclst(2*mxprj),asost(2*mxprj,2*mxprj),wsost(2*mxprj)
 real(dp) :: asclt(2*mxprj,2*mxprj),asot(2*mxprj,2*mxprj)
 real(dp) :: sphalf(2*mxprj,2*mxprj),smhalf(2*mxprj,2*mxprj)
 real(dp) :: fscl(mxprj),fso(mxprj),work(10*mxprj)
 real(dp), allocatable :: vkbt(:,:),vkbst(:,:)
 logical :: sorted
 character(len=500) :: msg

 ABI_UNUSED(irc(1))

 ! Check that rad grid is linear starting at zero
 !linear_mesh = .True.
 amesh = rr(2) - rr(1) !; damesh = zero
 !do irad=2,mmax-1
 !  damesh = max(damesh, abs(upf%r(irad)+amesh-upf%r(irad+1)))
 !end do
 !linear_mesh = damesh < tol8

 !if (.not. linear_mesh .or. abs(upf%r(1)) > tol16) then
 !  write(msg,'(3a)')&
 !  'Assuming pseudized valence charge given on linear radial mesh starting at zero.',ch10,&
 !  'Action: check your pseudopotential file.'
 !  ABI_ERROR(msg)
 !end if

 !allocate(vkbt(mmax,2*mxprj),vkbst(mmax,2*mxprj))
 ABI_MALLOC(vkbt, (mmax,2*mxprj))
 ABI_MALLOC(vkbst, (mmax,2*mxprj))

 do l1=1,lmax+1
  ll=l1-1

  if(ll==0) then
   vsr(:,:,l1)=0.0d0
   vso(:,:,l1)=0.0d0
   esr(:,l1)=0.0d0
   eso(:,l1)=0.0d0
   if(nproj(l1)>=1) then
    do ii=1,nproj(l1)
     vsr(:,ii,l1)=vkb(:,ii,l1,1)
     esr(ii,l1)=evkb(ii,l1,1)
    end do
   end if
   cycle
  end if

  nn=2*nproj(l1)

  fscl(1)=(ll+1)/dble(2*ll+1)
  fscl(2)=ll/dble(2*ll+1)
  fso(1)=2/dble(2*ll+1)
  fso(2)=-2/dble(2*ll+1)

! construct overlap matrix and diagonal energy matrices

   sovl(:,:)=0.0d0
   ascl(:,:)=0.0d0
   aso(:,:)=0.0d0
   vkbt(:,:)=0.0d0

   do ik1=1,2
    do ip1=1,nproj(l1)
     ii=ip1+(ik1-1)*nproj(l1)

     ascl(ii,ii)=fscl(ik1)*evkb(ip1,l1,ik1)
     aso(ii,ii)=fso(ik1)*evkb(ip1,l1,ik1)

     vkbt(:,ii)=vkb(:,ip1,l1,ik1)

     do ik2=1,2
      do ip2=1,nproj(l1)
       jj=ip2+(ik2-1)*nproj(l1)

       ! MG: This routine cannot be used as it assumes log mesh.
       !call vpinteg(vkb(1,ip1,l1,ik1),vkb(1,ip2,l1,ik2),irc(l1),2*l1, &
       !             sovl(ii,jj),rr)

       sovl(ii,jj) = simpson(amesh, vkb(:,ip1,l1,ik1) * vkb(:,ip2,l1,ik2))
       !print *, "vkb(1,ip1,l1,ik1),vkb(1,ip2,l1,ik2), irc(l1)"
       !print *, vkb(1,ip1,l1,ik1),vkb(1,ip2,l1,ik2), irc(l1)
       !print *, sovl(ii, jj), "ll", ll

      end do
     end do
    end do
   end do

   call dsyev( 'V', 'U', nn, sovl, 2*mxprj, sovlev, work, 10*mxprj, info )

   if(info .ne. 0) then
    write(msg,'(a,i4)') 'sr_so_r: S matrix eigenvalue ERROR, info=',info
    ABI_ERROR(msg)
   end if

! construct S^(-1/2) AND s^(1/2)

   do jj=1,nn
    tt=sqrt(sovlev(jj))
    do ii=1,nn
     sphalf(ii,jj)=tt*sovl(ii,jj)
     smhalf(ii,jj)=sovl(ii,jj)/tt
    end do
   end do

! take linear combinations to form orthonormal basis functions

   vkbst(:,:)=0.0d0

   do jj=1,nn
    do ii=1,nn
     vkbst(:,jj)=vkbst(:,jj) + smhalf(ii,jj)*vkbt(:,ii)
    end do
   end do

! construct A^(-1)* = S^(1/2)^T A^(-1) S^(1/2)

     asclt(:,:)=0.0d0
     asclst(:,:)=0.0d0
     asot(:,:)=0.0d0
     asost(:,:)=0.0d0

     do ii=1,nn
      do jj=1,nn
       do kk=1,nn
        asclt(ii,jj)=asclt(ii,jj)+ascl(ii,kk)*sphalf(kk,jj)
        asot(ii,jj) =asot(ii,jj) + aso(ii,kk)*sphalf(kk,jj)
       end do
      end do
     end do

     do ii=1,nn
      do jj=1,nn
       do kk=1,nn
        asclst(ii,jj)=asclst(ii,jj)+asclt(kk,jj)*sphalf(kk,ii)
        asost(ii,jj) =asost(ii,jj) + asot(kk,jj)*sphalf(kk,ii)
       end do
      end do
     end do

! find eigenvalues and eigenvectors of the A* matrices

!      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

     call dsyev( 'V', 'U', nn, asclst, 2*mxprj, wsclst, work, 10*mxprj, info )

     if(info .ne. 0) then
      write(msg,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
      ABI_ERROR(msg)
     end if

     call dsyev( 'V', 'U', nn,  asost, 2*mxprj,  wsost, work, 10*mxprj, info )


     if(info .ne. 0) then
      write(msg,'(a,i4)') 'sr_so_r: A* matrix eigenvalue ERROR, info=',info
      ABI_ERROR(msg)
     end if

! take linear combinations to form orthonormal projectors

     vsr(:,:,l1)=0.0d0
     vso(:,:,l1)=0.0d0
     esr(:,l1)=0.0d0
     eso(:,l1)=0.0d0

     do ii=1,nn
      esr(ii,l1)=wsclst(ii)
      eso(ii,l1)=  wsost(ii)
      do jj=1,nn
       vsr(:,ii,l1)=vsr(:,ii,l1) + asclst(jj,ii)*vkbst(:,jj)
       vso(:,ii,l1)= vso(:,ii,l1) +   asost(jj,ii)*vkbst(:,jj)
      end do
     end do

! bubble-sort on coefficient magnitudes for scalar and then s-o
! (Yes, I know bubble-sort is the least-efficient sorting algorithm.)

     do ii=1,100
      sorted=.true.
      do jj=2,nn
       if(abs(esr(jj-1,l1))<abs(esr(jj,l1))) then
        tt=esr(jj,l1)
        vkbt(:,1)=vsr(:,jj,l1)
        esr(jj,l1)=esr(jj-1,l1)
        vsr(:,jj,l1)=vsr(:,jj-1,l1)
        esr(jj-1,l1)=tt
        vsr(:,jj-1,l1)=vkbt(:,1)
        sorted=.false.
       end if
      end do
      if(sorted) exit
     end do

     do ii=1,100
      sorted=.true.
      do jj=2,nn
       if(abs(eso(jj-1,l1))<abs(eso(jj,l1))) then
        tt=eso(jj,l1)
        vkbt(:,1)=vso(:,jj,l1)
        eso(jj,l1)=eso(jj-1,l1)
        vso(:,jj,l1)=vso(:,jj-1,l1)
        eso(jj-1,l1)=tt
        vso(:,jj-1,l1)=vkbt(:,1)
        sorted=.false.
       end if
      end do
      if(sorted) exit
     end do

     write(std_out,'(/a,i2)') &
&         ' Orthonormal scalar projector coefficients, l = ',ll
     write(std_out,'(1p,6e12.4)') (esr(jj,l1),jj=1,nn)
     write(std_out,'(/a,i2)') &
&         ' Orthonormal spin-orbit projector coefficients, l = ',ll
     write(std_out,'(1p,6e12.4)') (eso(jj,l1),jj=1,nn)

! Set sign of projectors (physically irrelevant) so that they are positive
! at their peak (needed for compaisons apparently)

     do jj=1,nn
       apk=0.0d0
       do ii=1,mmax
         if(abs(vso(ii,jj,l1))>apk) then
           apk=abs(vso(ii,jj,l1))
           ipk=ii
         end if
       end do
       if(vso(ipk,jj,l1)<0.0d0) then
         vso(:,jj,l1)=-vso(:,jj,l1)
       end if
       apk=0.0d0
       do ii=1,mmax
         if(abs(vsr(ii,jj,l1))>apk) then
           apk=abs(vsr(ii,jj,l1))
           ipk=ii
         end if
       end do
       if(vsr(ipk,jj,l1)<0.0d0) then
         vsr(:,jj,l1)=-vsr(:,jj,l1)
       end if
     end do

 end do ! l1

 ABI_FREE(vkbt)
 ABI_FREE(vkbst)
 return
end subroutine sr_so_r
!!***

!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! subroutine vpinteg(gg,hh,nn,mm,ss,rr)
!
!! integrals that go into construction of Vanderbilt separable pseudopotential
!
!! product of functions gg*hh goes like rr**mm at rr -> 0
!! integral on usual log mesh from 1 to nn
!
!!Input variables
! integer,intent(in) :: nn,mm
! real(dp),intent(in) :: gg(nn),hh(nn),rr(nn)
!
!!Output variable
! real(dp),intent(out) :: ss
!
!!Local variables
! real(dp) :: r0,amesh,al
! integer :: ii
!
! al = 0.01d0 * dlog(rr(101)/rr(1))
! amesh = exp(al)
!
! r0=rr(1)/dsqrt(amesh)
! ss=r0**(mm+1)*(gg(1)*hh(1)/rr(1)**mm)/dfloat(mm+1)
!
! do ii = 4, nn - 3
!   ss =  ss + al*gg(ii)*hh(ii)*rr(ii)
! end do
!
! ss=ss + al*(23.d0*rr(nn-2)*gg(nn-2)*hh(nn-2) &
!&        + 28.d0*rr(nn-1)*gg(nn-1)*hh(nn-1) &
!&        +  9.d0*rr(nn  )*gg(nn  )*hh(nn  ))/24.d0
!
!
! return
! end subroutine vpinteg

end module m_pspheads
!!***
