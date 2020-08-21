!!****m* ABINIT/m_psps
!! NAME
!!  m_psps
!!
!! FUNCTION
!!  This module provides method to allocate/free/initialize the
!!  pseudopotential_type object.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (XG,DC,MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_psps

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_nctk
 use m_copy
 use m_dtset
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,      only : itoa, sjoin, yesno
 use m_io_tools,      only : open_file
 use m_symtk,         only : matr3inv
 use defs_datatypes,  only : pspheader_type, pseudopotential_type, pseudopotential_gth_type, nctab_t
 use m_paw_numeric,   only : paw_spline
 use m_pawrad,        only : pawrad_type, pawrad_init, pawrad_free, simp_gen
 use m_pawpsp,        only : pawpsp_cg
 use m_parser,        only : chkint_eq
 use m_memeval,       only : getdim_nloc, setmqgrid

 implicit none

 private

 ! Helper functions
 public :: test_xml_xmlpaw_upf     ! Test if a pseudo potential file is in XML, XML-PAW or in UPF format.

 public :: psps_init_global        ! Allocate and init all part of psps structure that are independent of a given dataset.
 public :: psps_init_from_dtset    ! Allocate and init all part of psps structure that are dependent of a given dataset.
 public :: psps_free               ! Deallocate all memory of psps structure.
 public :: psps_copy               ! Copy the psps structure.
 public :: psps_print              ! Print info on the pseudopotentials.
 public :: psps_ncwrite            ! Write psps data in netcdf format.

 public :: nctab_init              ! Create the object.
 public :: nctab_free              ! Free memory.
 public :: nctab_copy              ! Copy the object.
 public :: nctab_eval_tvalespl     ! Evaluate spline-fit of the atomic pseudo valence charge in reciprocal space.
 public :: nctab_eval_tcorespl     ! Evalute spline-fit of the model core charge in reciprocal space.
 public :: nctab_mixalch           ! Mix the pseudopotential tables. Used for alchemical mixing.
!!***

contains

!!****f* m_psps/test_xml_xmlpaw_upf
!! NAME
!!  test_xml_xmlpaw_upf
!!
!! FUNCTION
!!  Test if a pseudo potential file is in XML, XML-PAW or in UPF format.
!!
!! INPUTS
!!  path=Pseudopotential file
!!
!! OUTPUT
!!  usexml=1 if XML file
!!  xmlpaw=1 if PAW file in XML format
!!  useupf=1 if UPF file.
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine test_xml_xmlpaw_upf(path, usexml, xmlpaw, useupf)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 integer,intent(out) :: usexml, xmlpaw, useupf

!Local variables-------------------------------
!scalars
 integer :: temp_unit
 character(len=500) :: msg,errmsg
 character(len=70) :: testxml

! *************************************************************************

!  Check if the file pseudopotential file is written in XML
 usexml = 0; xmlpaw = 0; useupf = 0

 if (open_file(path,msg,newunit=temp_unit,form='formatted',status='old') /= 0) then
   MSG_ERROR(msg)
 end if
 rewind (unit=temp_unit,err=10,iomsg=errmsg)

 read(temp_unit,*,err=10,iomsg=errmsg) testxml
 if(testxml(1:5)=='<?xml')then
   usexml = 1
   read(temp_unit,*,err=10,iomsg=errmsg) testxml
   if(testxml(1:4)=='<paw') xmlpaw = 1
 else
   usexml = 0
 end if

 !Check if pseudopotential file is a Q-espresso UPF file
 rewind (unit=temp_unit,err=10,iomsg=errmsg)
 read(temp_unit,*,err=10,iomsg=errmsg) testxml ! just a string, no relation to xml.
 if(testxml(1:9)=='<PP_INFO>')then
   useupf = 1
 else
   useupf = 0
 end if

 close(unit=temp_unit,err=10,iomsg=errmsg)

 return

 ! Handle IO error
10 continue
 MSG_ERROR(errmsg)

end subroutine test_xml_xmlpaw_upf
!!***

!!****f* m_psps/psps_init_global
!! NAME
!! psps_init_global
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are independent of a given dataset.
!!
!! INPUTS
!! npsp=the number of read pseudo files.
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      m_driver
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_init_global(mtypalch, npsp, psps, pspheads)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtypalch,npsp
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer :: ii, mpsang, n1xccc

! *************************************************************************

!Allocation of some arrays independent of the dataset
 ABI_ALLOCATE(psps%filpsp,(npsp))
 ABI_ALLOCATE(psps%pspcod,(npsp))
 ABI_ALLOCATE(psps%pspdat,(npsp))
 ABI_ALLOCATE(psps%pspso,(npsp))
 ABI_ALLOCATE(psps%pspxc,(npsp))
 ABI_ALLOCATE(psps%title,(npsp))
 ABI_ALLOCATE(psps%zionpsp,(npsp))
 ABI_ALLOCATE(psps%znuclpsp,(npsp))
 call psp2params_init(psps%gth_params, npsp)

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

 ! Transfer md5 checksum
 ABI_ALLOCATE(psps%md5_pseudos, (npsp))
 psps%md5_pseudos = pspheads(1:npsp)%md5_checksum
!Set values independant from dtset
 psps%npsp   = npsp
!Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
!n1xccc=maxval(pspheads(1:npsp)%xccc)
 mpsang=1
 n1xccc=pspheads(1)%xccc
 do ii=1,psps%npsp
   mpsang=max(pspheads(ii)%lmax+1,mpsang)
   n1xccc=max(pspheads(ii)%xccc,n1xccc)
 end do
 psps%mpsang = mpsang
 psps%n1xccc = n1xccc
! Determine here whether the calculation is PAW
! If paw, all pspcod necessarily are 7 or 17 (see iofn2)
 psps%usepaw  =0
 if (pspheads(1)%pspcod==7.or.pspheads(1)%pspcod==17) psps%usepaw=1
 psps%mtypalch = mtypalch

end subroutine psps_init_global
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_init_from_dtset
!! NAME
!! psps_init_from_dtset
!!
!! FUNCTION
!! Allocate and initialise all part of psps structure that are dependent of a given dataset.
!!
!! INPUTS
!! dtset=<type dataset_type>a given dataset
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      m_driver
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_init_from_dtset(dtset, idtset, psps, pspheads)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 type(pspheader_type),intent(in) :: pspheads(psps%npsp)

!Local variables-------------------------------
!scalars
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgridff_old=0
 integer,save :: mqgridvl_old=0,ntypat_old=-1,usepaw_old=-1
 integer :: ipsp,lmnmax,lmnmaxso,lnmax,lnmaxso,newmqgrid,newmqgriddg,nptsgvec
 integer :: changed,ii,itypat
 real(dp) :: gprimd_orig(3,3)

! *************************************************************************

 psps%optnlxccc   = dtset%optnlxccc
!Determine the number of points needed in reciprocal space to represent the
!pseudopotentials (either set by hand from input variable or set automatically by abinit)
 nptsgvec         = 200 !This has to be chosen one and for all or else ??
 newmqgrid        = dtset%mqgrid
 newmqgriddg      = dtset%mqgriddg
 !JB:Which image to use ? I guess 1 always works
 call matr3inv(dtset%rprimd_orig(:,:,1),gprimd_orig)
 if ( dtset%usewvl == 0) then
   call setmqgrid(newmqgrid,newmqgriddg,dtset%ecut*dtset%dilatmx**2,&
&       dtset%pawecutdg*dtset%dilatmx**2,gprimd_orig,nptsgvec,psps%usepaw)
 else
   call setmqgrid(newmqgrid,newmqgriddg,one,one,gprimd_orig,nptsgvec,psps%usepaw)
 end if
 psps%mqgrid_ff   = newmqgrid
 if (psps%usepaw == 1) then
   psps%mqgrid_vl = newmqgriddg
 else
   psps%mqgrid_vl = newmqgrid
 end if

!Determine the maximum number of projectors, for the set of pseudo atom
 call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch_orig,dtset%nimage,psps%npsp,dtset%npspalch,&
& dtset%ntypat,dtset%ntypalch,pspheads)

 psps%npspalch = dtset%npspalch
 psps%ntypat   = dtset%ntypat
 psps%ntypalch = dtset%ntypalch
 psps%ntyppure = dtset%ntyppure

!Set the flag for reciprocal space or real space calculations
 psps%vlspl_recipSpace = (dtset%icoulomb /= 1)
!changed by RShaltaf
 psps%positron = dtset%positron
 psps%useylm   = dtset%useylm

!Added by T. Rangel for WVL+PAW
 psps%usewvl   = dtset%usewvl

! Define treatment of the model core density for NC pseudos.
 psps%nc_xccc_gspace = dtset%nc_xccc_gspace

 if (idtset > 1) then
   if (allocated(psps%algalch))  then
     ABI_DEALLOCATE(psps%algalch)
   end if
   if (allocated(psps%mixalch))  then
     ABI_DEALLOCATE(psps%mixalch)
   end if
 end if
 ABI_ALLOCATE(psps%algalch,(psps%ntypalch))
 ABI_ALLOCATE(psps%mixalch,(psps%npspalch,psps%ntypalch))
 psps%algalch(1:psps%ntypalch)=dtset%algalch(1:psps%ntypalch)
!This value will be overwritten elsewhere in case there are different images ...
 psps%mixalch(1:psps%npspalch,1:psps%ntypalch)=dtset%mixalch_orig(1:psps%npspalch,1:psps%ntypalch,1)

!Set mpspso and psps%pspso
!Warning: mpspso might be different for each dataset.
!         mpspso not relevant in case of PAW.
 psps%mpspso=1
 do ipsp=1,dtset%npsp
   if(dtset%nspinor==1)then
     psps%pspso(ipsp)=0
     ! This is needed to treate SOC perturbatively in sigma.
     !if (dtset%optdriver == RUNL_SIGMA .and. dtset%so_psp(ipsp) /= 0) then
     !  MSG_WARNING("Setting pspso to 2 although nspinor == 1")
     !  psps%pspso(ipsp) = 2
     !end if
!    Ideally the following line should not exist,
!      but at present, the space has to be booked
     if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
   else if (psps%usepaw==0) then
     if(dtset%so_psp(ipsp)/=1)then
       psps%pspso(ipsp)=dtset%so_psp(ipsp)
     else
       psps%pspso(ipsp)=pspheads(ipsp)%pspso
     end if
     if(psps%pspso(ipsp)/=0)psps%mpspso=2
     if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
   else
     psps%pspso(ipsp)=1+dtset%pawspnorb
   end if
 end do

!Set mpssoang, lmnmax, lnmax
 if(psps%mpspso==1)then
   psps%mpssoang=psps%mpsang
   psps%lmnmax  =lmnmax
   psps%lnmax   =lnmax
 else
   psps%mpssoang=2*psps%mpsang-1
   psps%lmnmax=lmnmaxso
   psps%lnmax=lnmaxso
 end if

!T. Rangel: for wvl + paw do not change psps%lmnmax
 if (psps%useylm==0 .and. psps%usepaw/=1 ) then
   psps%lmnmax=psps%lnmax
 end if

!Set dimekb
 if (psps%usepaw==0) then
   psps%dimekb=psps%lnmax
 else
   psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
 end if

!The following arrays are often not deallocated before the end of the dtset loop
!and might keep their content from one dataset to the other, if the conditions are fulfilled
 changed = 0

 if(dimekb_old/=psps%dimekb .or. ntypat_old/=dtset%ntypat .or. usepaw_old/=psps%usepaw) then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%ekb))  then
       ABI_DEALLOCATE(psps%ekb)
     end if
   end if
   ABI_ALLOCATE(psps%ekb,(psps%dimekb,dtset%ntypat*(1-psps%usepaw)))
   dimekb_old=psps%dimekb
 end if

 if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=dtset%ntypat)then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%indlmn))  then
       ABI_DEALLOCATE(psps%indlmn)
     end if
   end if
   ABI_ALLOCATE(psps%indlmn,(6,psps%lmnmax,dtset%ntypat))
   lmnmax_old=psps%lmnmax
 end if

 if(mqgridff_old/=psps%mqgrid_ff .or. lnmax_old/=psps%lnmax .or. ntypat_old/=dtset%ntypat)then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%ffspl))  then
       ABI_DEALLOCATE(psps%ffspl)
     end if
     if (allocated(psps%qgrid_ff))  then
       ABI_DEALLOCATE(psps%qgrid_ff)
     end if
   end if
   ABI_ALLOCATE(psps%ffspl,(psps%mqgrid_ff,2,psps%lnmax,dtset%ntypat))
   ABI_ALLOCATE(psps%qgrid_ff,(psps%mqgrid_ff))
   mqgridff_old=psps%mqgrid_ff
   lnmax_old=psps%lnmax
 end if

 if(mqgridvl_old/=psps%mqgrid_vl .or. ntypat_old/=dtset%ntypat)then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%qgrid_vl))  then
       ABI_DEALLOCATE(psps%qgrid_vl)
     end if
     if (allocated(psps%vlspl))  then
       ABI_DEALLOCATE(psps%vlspl)
     end if
     if (allocated(psps%nctab)) then
       do ii=1,size(psps%nctab)
         call nctab_free(psps%nctab(ii))
       end do
       ABI_FREE(psps%nctab)
     end if
   end if
   if (idtset/=1 .and. .not.psps%vlspl_recipSpace) then
     if (allocated(psps%dvlspl))  then
       ABI_DEALLOCATE(psps%dvlspl)
     end if
   end if

   ABI_ALLOCATE(psps%qgrid_vl,(psps%mqgrid_vl))
   ABI_ALLOCATE(psps%vlspl,(psps%mqgrid_vl,2,dtset%ntypat))

   if (psps%usepaw == 0) then
     ! If you change usepaw in the input, you will get what you deserve!
     ABI_MALLOC(psps%nctab, (dtset%ntypat))
     do itypat=1,dtset%ntypat
       call nctab_init(psps%nctab(itypat), psps%mqgrid_vl, .False., .False.)
     end do
   end if

   if (.not.psps%vlspl_recipSpace) then
     ABI_ALLOCATE(psps%dvlspl,(psps%mqgrid_vl,2,dtset%ntypat))
   end if
   mqgridvl_old=psps%mqgrid_vl
 end if

 if(ntypat_old/=dtset%ntypat.or. usepaw_old/=psps%usepaw)then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%xccc1d))  then
       ABI_DEALLOCATE(psps%xccc1d)
     end if
   end if
   ABI_ALLOCATE(psps%xccc1d,(psps%n1xccc*(1-psps%usepaw),6,dtset%ntypat))
   usepaw_old=psps%usepaw
 end if

 if(ntypat_old/=dtset%ntypat)then
   changed = changed + 1
   if(idtset/=1) then
     if (allocated(psps%xcccrc))  then
       ABI_DEALLOCATE(psps%xcccrc)
     end if
     if (allocated(psps%ziontypat))  then
       ABI_DEALLOCATE(psps%ziontypat)
     end if
     if (allocated(psps%znucltypat))  then
       ABI_DEALLOCATE(psps%znucltypat)
     end if
   end if
   ABI_ALLOCATE(psps%xcccrc,(dtset%ntypat))
   ABI_ALLOCATE(psps%znucltypat,(dtset%ntypat))
   ABI_ALLOCATE(psps%ziontypat,(dtset%ntypat))
   ntypat_old=dtset%ntypat
 end if

 psps%ziontypat(:)=dtset%ziontypat(:)

end subroutine psps_init_from_dtset
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_free
!! NAME
!! psps_free
!!
!! FUNCTION
!! Deallocate all memory of psps structure.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! psps=<type pseudopotential_type>the pseudopotentials description
!!
!! PARENTS
!!      m_ddb_hdr,m_driver
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_free(psps)

!Arguments ------------------------------------
!scalars
 type(pseudopotential_type),intent(inout) :: psps

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

!Allocation of some arrays independent of the dataset
 ABI_SFREE(psps%filpsp)
 ABI_SFREE(psps%pspcod)
 ABI_SFREE(psps%pspdat)
 ABI_SFREE(psps%pspso)
 ABI_SFREE(psps%pspxc)
 ABI_SFREE(psps%title)
 ABI_SFREE(psps%zionpsp)
 ABI_SFREE(psps%znuclpsp)
 ABI_SFREE(psps%algalch)
 ABI_SFREE(psps%mixalch)
 ABI_SFREE(psps%ekb)
 ABI_SFREE(psps%indlmn)
 ABI_SFREE(psps%ffspl)
 ABI_SFREE(psps%qgrid_ff)
 ABI_SFREE(psps%qgrid_vl)
 ABI_SFREE(psps%vlspl)
 ABI_SFREE(psps%dvlspl)
 ABI_SFREE(psps%xccc1d)
 ABI_SFREE(psps%xcccrc)
 ABI_SFREE(psps%ziontypat)
 ABI_SFREE(psps%znucltypat)
 ABI_SFREE(psps%md5_pseudos)

 ! Free types.
 call psp2params_free(psps%gth_params)

 if (allocated(psps%nctab)) then
   do ii=1,size(psps%nctab)
     call nctab_free(psps%nctab(ii))
   end do
   ABI_FREE(psps%nctab)
 end if

end subroutine psps_free
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_copy
!! NAME
!! psps_copy
!!
!! FUNCTION
!! Copy the psps structure.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_copy(pspsin, pspsout)

!Arguments ------------------------------------
!scalars
 type(pseudopotential_type),intent(in) :: pspsin
 type(pseudopotential_type),intent(out) :: pspsout

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

 ! integer
 pspsout%dimekb         = pspsin%dimekb
 pspsout%lmnmax         = pspsin%lmnmax
 pspsout%lnmax          = pspsin%lnmax
 pspsout%mproj          = pspsin%mproj
 pspsout%mpsang         = pspsin%mpsang
 pspsout%mpspso         = pspsin%mpspso
 pspsout%mpssoang       = pspsin%mpssoang
 pspsout%mqgrid_ff      = pspsin%mqgrid_ff
 pspsout%mqgrid_vl      = pspsin%mqgrid_vl
 pspsout%mtypalch       = pspsin%mtypalch
 pspsout%npsp           = pspsin%npsp
 pspsout%npspalch       = pspsin%npspalch
 pspsout%ntypat         = pspsin%ntypat
 pspsout%ntypalch       = pspsin%ntypalch
 pspsout%ntyppure       = pspsin%ntyppure
 pspsout%n1xccc         = pspsin%n1xccc
 pspsout%optnlxccc      = pspsin%optnlxccc
 pspsout%positron       = pspsin%positron
 pspsout%usepaw         = pspsin%usepaw
 pspsout%usewvl         = pspsin%usewvl
 pspsout%useylm         = pspsin%useylm
 pspsout%nc_xccc_gspace = pspsin%nc_xccc_gspace

 ! logical
 pspsout%vlspl_recipSpace = pspsin%vlspl_recipSpace

 ! integer allocatable
 if (allocated(pspsin%algalch)) then
   call alloc_copy(pspsin%algalch, pspsout%algalch)
 end if
 if (allocated(pspsin%indlmn)) then
   call alloc_copy(pspsin%indlmn, pspsout%indlmn)
 end if
 if (allocated(pspsin%pspdat)) then
   call alloc_copy(pspsin%pspdat, pspsout%pspdat)
 end if
 if (allocated(pspsin%pspcod)) then
   call alloc_copy(pspsin%pspcod, pspsout%pspcod)
 end if
 if (allocated(pspsin%pspso)) then
   call alloc_copy(pspsin%pspso, pspsout%pspso)
 end if
 if (allocated(pspsin%pspxc)) then
   call alloc_copy(pspsin%pspxc, pspsout%pspxc)
 end if

 ! real allocatable
 if (allocated(pspsin%ekb)) then
   call alloc_copy( pspsin%ekb, pspsout%ekb)
 end if
 if (allocated(pspsin%ffspl)) then
   call alloc_copy( pspsin%ffspl, pspsout%ffspl)
 end if
 if (allocated(pspsin%mixalch)) then
   call alloc_copy(pspsin%mixalch, pspsout%mixalch)
 end if
 if (allocated(pspsin%qgrid_ff)) then
   call alloc_copy(pspsin%qgrid_ff, pspsout%qgrid_ff)
 end if
 if (allocated(pspsin%qgrid_vl)) then
   call alloc_copy(pspsin%qgrid_vl, pspsout%qgrid_vl)
 end if
 if (allocated(pspsin%vlspl)) then
   call alloc_copy(pspsin%vlspl, pspsout%vlspl)
 end if
 if (allocated(pspsin%dvlspl)) then
   call alloc_copy(pspsin%dvlspl, pspsout%dvlspl)
 end if
 if (allocated(pspsin%xcccrc)) then
   call alloc_copy(pspsin%xcccrc, pspsout%xcccrc)
 end if
 if (allocated(pspsin%xccc1d)) then
   call alloc_copy(pspsin%xccc1d, pspsout%xccc1d)
 end if
 if (allocated(pspsin%zionpsp)) then
   call alloc_copy(pspsin%zionpsp, pspsout%zionpsp)
 end if
 if (allocated(pspsin%ziontypat)) then
   call alloc_copy(pspsin%ziontypat, pspsout%ziontypat)
 end if
 if (allocated(pspsin%znuclpsp)) then
   call alloc_copy(pspsin%znuclpsp, pspsout%znuclpsp)
 end if
 if (allocated(pspsin%znucltypat)) then
   call alloc_copy(pspsin%znucltypat, pspsout%znucltypat)
 end if

 ! allocate and copy character strings
 ABI_ALLOCATE(pspsout%filpsp,(pspsout%npsp))
 ABI_ALLOCATE(pspsout%title,(pspsout%npsp))
 ABI_ALLOCATE(pspsout%md5_pseudos,(pspsout%npsp))
 do ii=1,pspsout%npsp
   pspsout%filpsp(ii) = pspsin%filpsp(ii)
   pspsout%title(ii) = pspsin%title(ii)
   pspsout%md5_pseudos(ii) = pspsin%md5_pseudos(ii)
 end do

 ! allocate and copy objects
 if (allocated(pspsin%nctab)) then
   ABI_DATATYPE_ALLOCATE(pspsout%nctab,(pspsout%ntypat))
   do ii=1,pspsout%ntypat
     call nctab_copy(pspsin%nctab(ii), pspsout%nctab(ii))
   end do
 end if

 call psp2params_copy(pspsin%gth_params, pspsout%gth_params)

end subroutine psps_copy
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_print
!! NAME
!! psps_print
!!
!! FUNCTION
!!  Print the content of a pseudopotential_type derived type
!!
!! INPUTS
!!  psps=<type pseudopotential_type>=Info on the pseudopotentials.
!!  unit(optional)=unit number for output
!!  prtvol(optional)=verbosity level
!!  mode_paral(optional): either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_print(psps,unit,prtvol,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ierr,ips,ipsp_alch,ityp_alch,itypat,unt,my_prtvol
 character(len=4) :: mode
 character(len=500) :: msg
!arrays
 integer :: cond_values(4)
 character(len=9) :: cond_string(4)

! *************************************************************************

 ! Provide defaults
 my_prtvol=0; if (present(prtvol)) my_prtvol=prtvol
 unt=std_out; if (present(unit)) unt=unit
 mode='COLL'; if (present(mode_paral)) mode=mode_paral
 ierr=0; cond_string(1:4)=' '; cond_values(:)=0

 ! General info including spin-orbit
 call wrtout(unt,' ==== Info on pseudopotentials ==== ', mode)

 SELECT CASE (psps%usepaw)
 CASE (0)
   call wrtout(unt,'  Norm-conserving pseudopotentials ', mode)
   !call wrtout(unt, sjoin('  Max number of Kleinman-Bylander energies ', itoa(psps%dimekb)), mode)
   !do itypat=1,psps%ntypat
   ! write(msg,'(a,i4,a,f9.4)')' Type ',itypat,' K-B energies ',(psps%ekb(ikbe,itypat),ikbe=1,psps%dimekb)
   !end do
 CASE (1)
   write(msg,'(a)')
   call wrtout(unt,'  PAW calculation', mode)
   !call wrtout(unt,sjoin('  Max number of D_ij coefficients ', itoa(psps%dimekb)), mode)
 CASE DEFAULT
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usepaw',psps%usepaw,2,[0,1],unt)
 END SELECT

 !SELECT CASE (psps%positron)
 !CASE (0)
 !  call wrtout(unt, '  Standard Electron Calculation ', mode)
 !CASE (1,2)
 !  write(msg,'(a,i0)')'  Positron Calculation with positron .. ',psps%positron
 !  call wrtout(unt,msg,mode)
 !CASE DEFAULT
 !  call chkint_eq(0,0,cond_string,cond_values,ierr,'positron',psps%positron,3,[0,1,2],unt)
 !END SELECT

 write(msg,'(a,i4,2a,i4)')&
  '  Number of pseudopotentials .. ',psps%npsp,ch10,&
  '  Number of types of atoms   .. ',psps%ntypat
 call wrtout(unt,msg,mode)

 if (psps%usepaw==0) then
   SELECT CASE (psps%mpspso)
   CASE (1)
     call wrtout(unt,'  Scalar calculation (no spin-orbit term) ',mode)
   CASE (2)
     write(msg,'(3a,i3)')&
      '  Calculation with spin-orbit coupling ',ch10,&
      '  Max number of channels (spin-orbit included) ',psps%mpssoang
     call wrtout(unt,msg,mode)
     do itypat=1,psps%ntypat
       if (psps%pspso(itypat) /= 1) then
         write(msg,'(a,i4,a,i2,a)')&
          '  - Atom type ',itypat,' has spin-orbit characteristics (pspso= ',psps%pspso(itypat),")"
         call wrtout(unt,msg,mode)
       end if
     end do
   CASE DEFAULT
     call chkint_eq(0,0,cond_string,cond_values,ierr,'mpspso',psps%mpspso,2,[1,2],unt)
   END SELECT
 else
   SELECT CASE (maxval(psps%pspso))
   CASE (0,1)
     msg='  Scalar calculation (no spin-orbit term) '
   CASE (2)
     msg='  Calculation with spin-orbit coupling '
   END SELECT
   call wrtout(unt,msg,mode)
 end if

 ! Info on nonlocal part
 SELECT CASE (psps%useylm)
 CASE (0)
   msg = '  Nonlocal part applied using Legendre polynomials '
 CASE (1)
   msg = '  Nonlocal part applied using real spherical harmonics '
 CASE DEFAULT
   call chkint_eq(0,0,cond_string,cond_values,ierr,'psps%useylm',psps%useylm,2,(/0,1/),unt)
 END SELECT
 call wrtout(unt,msg,mode)

 write(msg,'(a,i3)')'  Max number of non-local projectors over l and type ',psps%mproj
 call wrtout(unt,msg,mode)

 write(msg,'(a,i3,2a,i3,2a,i3)')&
 '  Highest angular momentum +1 ....... ',psps%mpsang,ch10,&
 '  Max number of (l,n)   components .. ',psps%lnmax, ch10,&
 '  Max number of (l,m,n) components .. ',psps%lmnmax
 call wrtout(unt,msg,mode)

 !FIXME for paw n1xccc==1
 ! Non-linear Core correction
 if (psps%n1xccc/=0) then
   write(msg,'(3a,2(a,i4,a),2a)')ch10,&
    ' Pseudo-Core Charge Info: ',ch10,&
    '   Number of radial points for pseudo-core charge .. ',psps%n1xccc,ch10,&
    '   XC core-correction treatment (optnlxccc) ........ ',psps%optnlxccc,ch10,&
    '   Radius for pseudo-core charge for each type ..... ',ch10
   call wrtout(unt,msg,mode)
   do itypat=1,psps%ntypat
     write(msg,'(a,i4,a,f7.4)')'  - Atom type ',itypat,' has pseudo-core radius .. ',psps%xcccrc(itypat)
     call wrtout(unt,msg,mode)
   end do
 end if

 ! Alchemical mixing
 if (psps%mtypalch/=0) then
   write(msg,'(3a,3(a,i4,a))')ch10,&
    ' Calculation with alchemical mixing:',ch10,&
    '   Number of pure pseudoatoms .... ',psps%ntyppure,ch10,&
    '   Number of pseudos for mixing .. ',psps%npspalch,ch10,&
    '   Alchemical pseudoatoms ........ ',psps%ntypalch,ch10
   call wrtout(unt,msg,mode)
   do ipsp_alch=1,psps%npspalch
     do ityp_alch=1,psps%ntypalch
       write(std_out,*)' mixalch ',psps%mixalch(ipsp_alch,ityp_alch)
     end do
   end do
   do ityp_alch=1,psps%ntypalch
     write(msg,'(a,i4,a,i4)')' For alchemical atom no. ',ityp_alch,' algalch is .. ',psps%algalch(ityp_alch)
     call wrtout(unt,msg,mode)
   end do
 end if

 ! Info in Q-grid for spline of form factors
 write(msg,'(3a,a,i6,a,a,i6)')ch10,&
  ' Info on the Q-grid used for form factors in spline form: ',ch10,&
  '   Number of q-points for radial functions ffspl .. ',psps%mqgrid_ff,ch10,&
  '   Number of q-points for vlspl ................... ',psps%mqgrid_vl
 call wrtout(unt,msg,mode)

 if (psps%vlspl_recipSpace) then
   call wrtout(unt,'   vloc is computed in Reciprocal Space ',mode)
 else
   call wrtout(unt,'   vloc is computed in Real Space ',mode)
 end if
 if (psps%usepaw == 0) then
   if (psps%nc_xccc_gspace == 0) call wrtout(unt,'   model core charge treated in real-space', mode)
   if (psps%nc_xccc_gspace == 1) call wrtout(unt,'   model core charge treated in G-space', mode)
 end if

 !TODO additional stuff that might be printed
 call wrtout(unt, "", mode)
 do itypat=1,psps%ntypat
   write(msg,'(a,i0,a,i0)')'  XC functional for type ',itypat,' is ',psps%pspxc(itypat)
   call wrtout(unt,msg,mode)
   !write(std_out,*)psps%ziontypat(itypat),psps%znucltypat(itypat)
   if (psps%usepaw == 0) then
     call wrtout(unt, sjoin("  Pseudo valence available: ", yesno(psps%nctab(itypat)%has_tvale)), mode)
   end if
 end do

 !integer, pointer :: pspxc(:)
 ! pspxc(ntypat)
 ! For each type of psp, the XC functional that was used to generate it, as given by the psp file
 if (my_prtvol>=3) then
   do ips=1,psps%npsp
     write(std_out,*)' Pseudo number   ',ips,' read from ',trim(psps%filpsp(ips))
     write(std_out,*)' Format or code  ',psps%pspcod(ips)
     write(std_out,*)' Generation date ',psps%pspdat(ips)
     write(std_out,*)' Content of first line: ', trim(psps%title(ips))
   end do
 end if

 call wrtout(unt, "", mode)

end subroutine psps_print
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psps_ncwrite
!! NAME
!! psps_ncwrite
!!
!! FUNCTION
!!  Writes on file the most important arrays defined in the pseudopotential_type
!!  for futher post-processing. This function should be called by master node only.
!!
!! INPUTS
!!   path=File name.
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psps_ncwrite(psps, path)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
!scalars
 integer :: ipsp,itypat,ncid,ncerr

! *************************************************************************

#ifdef HAVE_NETCDF
 NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))

 ! Define dimensions
 ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("fnlen", fnlen+1), nctkdim_t("md5_slen", md5_slen+1), nctkdim_t("ntypat", psps%ntypat), &
     nctkdim_t("npsp", psps%npsp), nctkdim_t("lnmax", psps%lnmax), &
     nctkdim_t("lmnmax", psps%lnmax), nctkdim_t("dimekb", psps%dimekb), &
     nctkdim_t("mqgrid_vl", psps%mqgrid_vl), nctkdim_t("mqgrid_ff", psps%mqgrid_ff) &
 ])
 NCF_CHECK(ncerr)

 if (psps%n1xccc /= 0) then ! 0 means unlimited!
   NCF_CHECK(nctk_def_dims(ncid, nctkdim_t("n1xccc", psps%n1xccc)))
 end if

 ! Define variables
 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "usepaw", "useylm"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("ziontypat", "dp", "ntypat"), &
   nctkarr_t("znucltypat", "dp", "ntypat"), &
   nctkarr_t("qgrid_vl", "dp", "mqgrid_vl"), &
   nctkarr_t("qgrid_ff", "dp", "mqgrid_ff"), &
   nctkarr_t("vlspl", "dp", "mqgrid_vl, two, ntypat"), &
   nctkarr_t("xcccrc", "dp", "ntypat"), &
   nctkarr_t("indlmn", "int", "six, lmnmax, ntypat"), &
   nctkarr_t("ffspl", "dp", "mqgrid_ff, two, lnmax, ntypat"), &
   nctkarr_t("filpsp", "char", "fnlen, npsp"), &
   nctkarr_t("md5_pseudos", "char", "md5_slen, npsp") &
 ])
 NCF_CHECK(ncerr)

 if (psps%usepaw == 0) then
   NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("ekb", "dp", "dimekb, ntypat")))
   if (psps%n1xccc /= 0) then
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("xccc1d", "dp", "n1xccc, six, ntypat")))
   end if
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("nc_tvalespl", "dp", "mqgrid_vl, two, ntypat"), &
     nctkarr_t("nc_tcorespl", "dp", "mqgrid_vl, two, ntypat")  &
   ])
   NCF_CHECK(ncerr)
 end if

 ! Write data
 NCF_CHECK(nf90_enddef(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("usepaw"), psps%usepaw))
 NCF_CHECK(nf90_put_var(ncid, vid("useylm"), psps%useylm))
 NCF_CHECK(nf90_put_var(ncid, vid("ziontypat"), psps%ziontypat))
 NCF_CHECK(nf90_put_var(ncid, vid("znucltypat"), psps%znucltypat))
 do ipsp=1,psps%npsp
   NCF_CHECK(nf90_put_var(ncid, vid("filpsp"), trim(psps%filpsp(ipsp)), start=[1, ipsp]))
   NCF_CHECK(nf90_put_var(ncid, vid("md5_pseudos"), trim(psps%md5_pseudos(ipsp)), start=[1, ipsp]))
 end do
 NCF_CHECK(nf90_put_var(ncid, vid("qgrid_vl"), psps%qgrid_vl))
 NCF_CHECK(nf90_put_var(ncid, vid("qgrid_ff"), psps%qgrid_ff))
 NCF_CHECK(nf90_put_var(ncid, vid("indlmn"), psps%indlmn))

 ! Local part in q-space and second derivative
 if (allocated(psps%vlspl)) then
   NCF_CHECK(nf90_put_var(ncid, vid("vlspl"), psps%vlspl))
 end if

 ! Form factors for each type of atom
 ! for each type and each (l,n) channel, ffnl(q) and second derivative
 if (allocated(psps%ffspl)) then
   NCF_CHECK(nf90_put_var(ncid, vid("ffspl"), psps%ffspl))
 end if

 ! Pseudo-core charge for each type of atom, on the real-space radial
 NCF_CHECK(nf90_put_var(ncid, vid("xcccrc"), psps%xcccrc))
 if (psps%usepaw == 0 .and. allocated(psps%xccc1d)) then
   NCF_CHECK(nf90_put_var(ncid, vid("xccc1d"), psps%xccc1d))
 end if

 ! NC-only: add tcore_spl and tvalespl in q-space
 if (psps%usepaw == 0) then
   NCF_CHECK(nf90_put_var(ncid, vid("ekb"), psps%ekb))
   do itypat=1,psps%ntypat
     if (psps%nctab(itypat)%has_tvale) then
       ncerr = nf90_put_var(ncid, vid("nc_tvalespl"), psps%nctab(itypat)%tvalespl, start=[1,1,itypat])
       NCF_CHECK(ncerr)
     end if
     if (psps%nctab(itypat)%has_tcore) then
       ncerr = nf90_put_var(ncid, vid("nc_tcorespl"), psps%nctab(itypat)%tcorespl, start=[1,1,itypat])
       NCF_CHECK(ncerr)
     end if
   end do
 end if

 NCF_CHECK(nf90_close(ncid))

#else
 MSG_WARNING("netcdf support not activated. psps file cannot be created!")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine psps_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psp2params_init
!! NAME
!! psp2params_init
!!
!! FUNCTION
!! Allocate and initialise the data structure holding parameters for the GTH
!! pseudo-potentials.
!!
!!  MJV note: this should be renamed: psp2 suggests it relates to pspcod 2,
!!     whereas it is actually 3
!!    the parameters would also be better off separated into C and h arrays
!!
!! INPUTS
!!  npsp=number of true pseudo used (not alchemy).
!!
!! OUTPUT
!!  gth_params <type (pseudopotential_gth_type)>=the values to allocate and initialise.
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psp2params_init(gth_params, npsp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
 type(pseudopotential_gth_type),intent(out) :: gth_params

! *********************************************************************

!Check array, no params are currently set.
 ABI_ALLOCATE(gth_params%set,(npsp))
 gth_params%set(:) = .false.

!Check array, have geometric information been filled?
 ABI_ALLOCATE(gth_params%hasGeometry,(npsp))
 gth_params%hasGeometry(:) = .false.

!Coefficients for local part and projectors
 ABI_ALLOCATE(gth_params%psppar,(0:4, 0:6, npsp))
 gth_params%psppar = zero

!Coefficients for spin orbit part
 ABI_ALLOCATE(gth_params%psp_k_par,(1:4, 1:3, npsp))
 gth_params%psp_k_par = zero

!Different radii
 ABI_ALLOCATE(gth_params%radii_cf,(npsp, 3))

end subroutine psp2params_init
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psp2params_copy
!! NAME
!! psp2params_copy
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psp2params_copy(gth_paramsin, gth_paramsout)

!Arguments ------------------------------------
!scalars
 type(pseudopotential_gth_type),intent(in) :: gth_paramsin
 type(pseudopotential_gth_type),intent(out) :: gth_paramsout

! *********************************************************************

 if (allocated(gth_paramsin%psppar)) then
   call alloc_copy( gth_paramsin%psppar, gth_paramsout%psppar)
 end if
 if (allocated(gth_paramsin%radii_cf)) then
   call alloc_copy( gth_paramsin%radii_cf, gth_paramsout%radii_cf)
 end if
 if (allocated(gth_paramsin%psp_k_par)) then
   call alloc_copy( gth_paramsin%psp_k_par, gth_paramsout%psp_k_par)
 end if
 if (allocated(gth_paramsin%hasGeometry)) then
   call alloc_copy( gth_paramsin%hasGeometry, gth_paramsout%hasGeometry)
 end if
 if (allocated(gth_paramsin%set)) then
   call alloc_copy( gth_paramsin%set, gth_paramsout%set)
 end if

end subroutine psp2params_copy
!!***

!----------------------------------------------------------------------

!!****f* m_psps/psp2params_free
!! NAME
!! psp2params_free
!!
!! FUNCTION
!! Deallocate a previously allocated data structure for storage of GTH parameters.
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  gth_params <type (pseudopotential_gth_type)>=the values to deallocate.
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine psp2params_free(gth_params)

!Arguments ------------------------------------
!scalars
 type(pseudopotential_gth_type),intent(inout) :: gth_params

! *********************************************************************

 ABI_SFREE(gth_params%set)
 ABI_SFREE(gth_params%hasGeometry)

!Coefficients for local part and projectors
 ABI_SFREE(gth_params%psppar)

!Coefficients for spin orbit part
 ABI_SFREE(gth_params%psp_k_par)

!Different radii
 ABI_SFREE(gth_params%radii_cf)

end subroutine psp2params_free
!!***

!!****f* m_psps/nctab_init
!! NAME
!!  nctab_init
!!
!! FUNCTION
!!  Create nctab_t.
!!
!! INPUTS
!!  mqgrid_vl=Number of q-points
!!  has_tcore=True if the pseudo has NLCC.
!!  has_tvale=True if the atomic valence density is available.
!!
!! PARENTS
!!      m_pspini,m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_init(nctab, mqgrid_vl, has_tcore, has_tvale)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mqgrid_vl
 logical,intent(in) :: has_tcore,has_tvale
 type(nctab_t),intent(inout) :: nctab

! *************************************************************************

 nctab%mqgrid_vl = mqgrid_vl

 ! The array for the model core charge is always allocated and initialized with zeros.
 ! This approach is similar to the one used in the PAW code.
 ! has_tcore tells us whether the model core charge is present or not.
 nctab%has_tcore = has_tcore
 nctab%dncdq0 = zero; nctab%d2ncdq0 = zero
 ABI_CALLOC(nctab%tcorespl, (mqgrid_vl, 2))
 nctab%tcorespl = zero

 ! tvalespl is allocated only if available.
 nctab%has_tvale = has_tvale
 nctab%dnvdq0 = zero
 if (has_tvale) then
   ABI_MALLOC(nctab%tvalespl, (mqgrid_vl, 2))
   nctab%tvalespl = zero
 end if

end subroutine nctab_init
!!***

!!****f* m_psps/nctab_free
!! NAME
!!  nctab_free
!!
!! FUNCTION
!! Free memory allocated in nctab_t
!!
!! PARENTS
!!      m_pspini,m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_free(nctab)

!Arguments ------------------------------------
!scalars
 type(nctab_t),intent(inout) :: nctab

! *************************************************************************

 nctab%mqgrid_vl = 0
 nctab%has_tvale = .False.
 ABI_SFREE(nctab%tvalespl)
 nctab%has_tcore = .False.
 ABI_SFREE(nctab%tcorespl)

end subroutine nctab_free
!!***

!!****f* m_psps/nctab_copy
!! NAME
!!  nctab_copy
!!
!! FUNCTION
!!
!! PARENTS
!!      m_psps
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_copy(nctabin, nctabout)

!Arguments ------------------------------------
!scalars
 type(nctab_t),intent(in) :: nctabin
 type(nctab_t),intent(out) :: nctabout

! *************************************************************************

 nctabout%mqgrid_vl  = nctabin%mqgrid_vl
 nctabout%has_tvale  = nctabin%has_tvale
 nctabout%has_tcore  = nctabin%has_tcore
 nctabout%dncdq0     = nctabin%dncdq0
 nctabout%d2ncdq0    = nctabin%d2ncdq0
 nctabout%dnvdq0     = nctabin%dnvdq0

 if (allocated(nctabin%tvalespl)) then
   call alloc_copy(nctabin%tvalespl, nctabout%tvalespl)
 end if
 if (allocated(nctabin%tcorespl)) then
   call alloc_copy(nctabin%tcorespl, nctabout%tcorespl)
 end if

end subroutine nctab_copy
!!***

!!****f* m_psps/nctab_eval_tvalespl
!! NAME
!!  nctab_eval_tvalespl
!!
!! FUNCTION
!!  Evalute spline-fit of the atomic pseudo valence charge in reciprocal space.
!!
!! INPUTS
!!  zion=nominal valence of atom as specified in psp file. Used to rescale the f(q=0) component
!!  mesh<pawrad_type>Radial mesh (r-space) used for the valence denity.
!!  valr(mesh%mesh_size)=Valence density in real space.
!!  mqgrid_vl=Number of points in the reciprocal space grid
!!  qgrid_vl(mqgrid_vl)=The coordinates of all the points of the radial q-grid
!!
!! SIDE EFFECTS
!!  nctabl%tvalspl(mqgrid_vl,2)
!!  nctab%d2ncdq0
!!
!! PARENTS
!!      m_psp8,m_psp9
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_eval_tvalespl(nctab, zion, mesh, valr, mqgrid_vl, qgrid_vl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mqgrid_vl
 real(dp),intent(in) :: zion
 type(nctab_t),intent(inout) :: nctab
 type(pawrad_type),intent(in) :: mesh
!arrays
 real(dp),intent(in) :: valr(mesh%mesh_size),qgrid_vl(mqgrid_vl)

!Local variables-------------------------------
!scalars
 real(dp) :: fact,yp1,ypn,d2nvdq0

! *************************************************************************

 nctab%has_tvale = .True.
 if (.not. allocated(nctab%tvalespl)) then
   ABI_MALLOC(nctab%tvalespl, (mqgrid_vl, 2))
 else
   ABI_CHECK(size(nctab%tvalespl, dim=1) == mqgrid_vl, "wrong mqgrid_vl")
 end if

 call pawpsp_cg(nctab%dnvdq0,d2nvdq0,mqgrid_vl,qgrid_vl,nctab%tvalespl(:,1),mesh,valr,yp1,ypn)
 call simp_gen(yp1,mesh%rad**2 * valr, mesh); write(std_out,*)" valence charge integrates to: ",four_pi*yp1

 ! Rescale the integral to have the correct number of valence electrons.
 ! In some cases, indeed, the radial mesh is not large enough and some valence charge is missing
 ! pawpsp_cg extrapolates the integrand beyond rmax but this is not enough.
 ! Remember that tvalespl is used to build an initial guess for rhor hence it's very important
 ! to have the correct electrostatic.
 fact = zion / nctab%tvalespl(1,1)
 nctab%tvalespl(:,1) = nctab%tvalespl(:,1) * fact

 ! Compute second derivative of tvalespl(q)
 call paw_spline(qgrid_vl,nctab%tvalespl(:,1),mqgrid_vl,yp1,ypn,nctab%tvalespl(:,2))

end subroutine nctab_eval_tvalespl
!!***

!!****f* m_psps/nctab_eval_tcorespl
!! NAME
!!  nctab_eval_tcorespl
!!
!! FUNCTION
!!  Evalute spline-fit of the model core charge in reciprocal space.
!!
!! INPUTS
!!  xcccrc=maximum radius of the pseudo-core charge
!!  n1xccc=Number of radial points for the description of the pseudo-core charge
!!     (in the framework of the non-linear XC core correction)
!!  mqgrid_vl=Number of points in the reciprocal space grid
!!  qgrid_vl(mqgrid_vl)=The coordinates of all the points of the radial q-grid
!!  xccc1d(n1xccc,6)= The component xccc1d(n1xccc,1) is the pseudo-core charge
!!   on the radial grid. The components xccc1d(n1xccc,ideriv) give the ideriv-th derivative of the
!!   pseudo-core charge with respect to the radial distance.
!!
!! SIDE EFFECTS
!!  nctabl%tcorespl(mqgrid_vl,2)
!!  nctab%d2ncdq0
!!  nctab%dncdq0
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_eval_tcorespl(nctab, n1xccc, xcccrc, xccc1d, mqgrid_vl, qgrid_vl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1xccc,mqgrid_vl
 real(dp),intent(in) :: xcccrc
 type(nctab_t),intent(inout) :: nctab
!arrays
 real(dp),intent(in) :: xccc1d(n1xccc,6),qgrid_vl(mqgrid_vl)

!Local variables-------------------------------
!scalars
 real(dp) :: amesh,yp1,ypn
 type(pawrad_type) :: core_mesh

! *************************************************************************

 ABI_CHECK(mqgrid_vl == nctab%mqgrid_vl, "wrong mqgrid_vl")

 if (.not. allocated(nctab%tcorespl)) then
   ABI_CALLOC(nctab%tcorespl, (mqgrid_vl, 2))
 else
   ABI_CHECK(size(nctab%tcorespl, dim=1) == mqgrid_vl, "wrong mqgrid_vl")
 end if

 ! Skip loop if this atom has no core charge
 if (abs(xcccrc) < tol16) then
   nctab%has_tcore = .False.
   return
 end if

 nctab%has_tcore = .True.
 ! XCCC is given on a linear mesh.
 amesh = xcccrc / dble(n1xccc-1)
 call pawrad_init(core_mesh, mesh_size=n1xccc, mesh_type=1, rstep=amesh)

 ! Compute 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 n(r))dr].
 !write(std_out,*)"xccc1d: amesh, min, max, minloc ",amesh,maxval(xccc1d(:,1)),minval(xccc1d(:,1)),minloc(xccc1d(:,1))
 call pawpsp_cg(nctab%dncdq0, nctab%d2ncdq0, mqgrid_vl, qgrid_vl, nctab%tcorespl(:,1), &
                core_mesh, xccc1d(:,1), yp1, ypn)

 ! Compute second derivative of tcorespl(q)
 call paw_spline(qgrid_vl, nctab%tcorespl(:,1), mqgrid_vl, yp1, ypn, nctab%tcorespl(:,2))
 call pawrad_free(core_mesh)

end subroutine nctab_eval_tcorespl
!!***

!!****f* m_psps/nctab_mixalch
!! NAME
!!  nctab_mixalch
!!
!! FUNCTION
!!  Mix the pseudopotential tables. Used for alchemical mixing.
!!
!! INPUTS
!! nctabs(npspalch)=NC tables to be mixed
!! npspalch=Number of alchemical pseudos.
!! ntypalch=Number of types of alchemical pseudoatoms
!! algalch(ntypalch)=For each type of pseudo atom, the algorithm to mix the pseudopotentials
!! mixalch(npspalch,ntypalch)=Mixing coefficients to generate alchemical pseudo atoms
!!
!! OUTPUT
!! mixtabs(ntypalch)=NC tables describing the alchemical pseudos
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      nctab_free,nctab_init
!!
!! SOURCE

subroutine nctab_mixalch(nctabs, npspalch, ntypalch, algalch, mixalch, mixtabs)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npspalch,ntypalch
!arrays
 integer,intent(in) :: algalch(ntypalch)
 real(dp),intent(in) :: mixalch(npspalch, ntypalch)
 type(nctab_t),intent(in) :: nctabs(npspalch)
 type(nctab_t),target,intent(inout) :: mixtabs(ntypalch)

!Local variables-------------------------------
!scalars
 integer :: ipspalch,itypalch
 logical :: has_tcore, has_tvale
 real(dp) :: mc
 type(nctab_t),pointer :: mix

! *************************************************************************

 ABI_CHECK(all(nctabs(:)%mqgrid_vl == nctabs(1)%mqgrid_vl), "Wrong mqgrid_vl")
 ABI_CHECK(all(algalch == 1), "algalch /= 1 not implemented")

 do itypalch=1,ntypalch

   ! has_tcore is true is at least one pseudo has nlcc.
   ! has_tvale is true if *all* mixed pseudos have the PS valence charge.
   has_tcore = .False.; has_tvale = .True.
   do ipspalch=1,npspalch
     if (abs(mixalch(ipspalch,itypalch)) < tol6) cycle
     if (nctabs(ipspalch)%has_tcore) has_tcore = .True.
     if (.not. nctabs(ipspalch)%has_tvale) has_tvale = .False.
   end do
   !write(std_out,*)has_tvale, has_tcore

   call nctab_free(mixtabs(itypalch))
   call nctab_init(mixtabs(itypalch), nctabs(1)%mqgrid_vl, has_tcore, has_tvale)
   mix => mixtabs(itypalch)

   do ipspalch=1,npspalch
     mc = mixalch(ipspalch,itypalch)
     if (abs(mc) < tol6) cycle
     ! Linear combination of the quantities
     ! Mix core for NLCC
     if (has_tcore) then
       mix%tcorespl = mix%tcorespl + mc * nctabs(ipspalch)%tcorespl
       mix%dncdq0 = mix%dncdq0 + mc * nctabs(ipspalch)%dncdq0
       mix%d2ncdq0 = mix%d2ncdq0 + mc * nctabs(ipspalch)%d2ncdq0
     end if
     ! Mix pseudo valence charge.
     if (has_tvale) then
       mix%tvalespl = mix%tvalespl + mc * nctabs(ipspalch)%tvalespl
       mix%dnvdq0 = mix%dnvdq0 + mc * nctabs(ipspalch)%dnvdq0
     end if
   end do

 end do

end subroutine nctab_mixalch
!!***

end module m_psps
!!***
