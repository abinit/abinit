!{\Src2tex{textfont=tt}}
!!****f* ABINIT/invars1
!! NAME
!! invars1
!!
!! FUNCTION
!! Initialize the dimensions needed to allocate the input arrays
!! for one dataset characterized by jdtset, by
!! taking from string the necessary data.
!! Perform some preliminary checks and echo these dimensions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, AR, MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iout=unit number of output file
!!  jdtset=number of the dataset looked for
!!  lenstr=actual length of string
!!  msym=default maximal number of symmetries
!!  npsp1= number of pseudopotential files
!!  zionpsp(npsp1)= valence charge over all psps
!!
!! OUTPUT
!!  mband_upper=estimation of the maximum number of bands for any k-point
!!
!! SIDE EFFECTS
!! Input/Output (the default value is given in the calling routine)
!!  dtset=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized, while some others will still be initialized later.
!!   The list of records of dtset initialized in the present routine is:
!!   acell_orig,densty,iatfix,kptopt,kptrlatt,
!!   mkmem,mkqmem,mk1mem,natsph,natvshift,nconeq,nkpt,nkptgw,nkpthf,
!!   nqptdm,nshiftk,nucdipmom,nzchempot,optdriver,
!!   rprim_orig,rprimd_orig,shiftk,
!!   spgroup,spinat,typat,vel_orig,vel_cell_orig,xred_orig
!!  bravais(11)=characteristics of Bravais lattice (see symlatt.F90)
!!  symafm(1:msym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,1:msym)=symmetry operations in real space in terms
!!   of primitive translations
!!  tnons(3,1:msym)=nonsymmorphic translations for symmetry operations
!!  string*(*)=string of characters containing all input variables and data
!!
!! NOTES
!! Must set up the geometry of the system, needed to compute
!! k point grids in an automatic fashion.
!! Treat separately mband_upper, since
!! fband, charge and zionpsp must be known for being able to initialize it.
!!
!! Defaults are provided in the calling routine.
!! Defaults are also provided here for the following variables :
!! mband_upper, occopt, fband, charge
!! They should be kept consistent with defaults of the same variables
!! provided to the invars routines.
!!
!! PARENTS
!!      invars1m
!!
!! CHILDREN
!!      atomdata_from_znucl,chkint_ge,ingeo,inkpts,inqpt,intagm,inupper
!!      invacuum,mkrdim,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine invars1(bravais,dtset,iout,jdtset,lenstr,mband_upper,msym,npsp1,&
& string,symafm,symrel,tnons,zionpsp)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_atomdata

 use m_fstrings, only : inupper

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars1'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => invars1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,lenstr,msym,npsp1
 integer,intent(out) :: mband_upper
 character(len=*),intent(inout) :: string
 type(dataset_type),intent(inout) :: dtset
!arrays
 integer,intent(inout) :: bravais(11),symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)
 real(dp),intent(in) :: zionpsp(npsp1)

!Local variables-------------------------------
 character :: blank=' ',string1
!scalars
 integer :: chksymbreak,found,ierr,iatom,ii,ikpt,iimage,index_blank,index_lower
 integer :: index_typsymb,index_upper,ipsp,iscf,intimage,itypat,leave,marr
 integer :: natom,nkpt,nkpthf,npsp,npspalch
 integer :: nqpt,nspinor,nsppol,ntypat,ntypalch,ntyppure,occopt,response
 integer :: rfddk,rfelfd,rfphon,rfstrs,rfuser,rf2_dkdk,rf2_dkde,rfmagn
 integer :: tfband,tnband,tread,tread_alt
 real(dp) :: charge,fband,kptnrm,kptrlen,zelect,zval
 character(len=2) :: string2,symbol
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 integer :: cond_values(4),vacuum(3)
 integer,allocatable :: iatfix(:,:),intarr(:),istwfk(:),nband(:),typat(:)
 real(dp) :: acell(3),rprim(3,3)
!real(dp) :: field(3)
 real(dp),allocatable :: amu(:),dprarr(:),kpt(:,:),kpthf(:,:),mixalch(:,:),nucdipmom(:,:)
 real(dp),allocatable :: ratsph(:),reaalloc(:),spinat(:,:)
 real(dp),allocatable :: vel(:,:),vel_cell(:,:),wtk(:),xred(:,:),znucl(:)
 character(len=32) :: cond_string(4)


!************************************************************************

!Some initialisations
 ierr=0
 cond_string(1:4)=' '
 cond_values(1:4)=(/0,0,0,0/)

!Read parameters
 marr=dtset%npsp;if (dtset%npsp<3) marr=3
 marr=max(marr,dtset%nimage)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!---------------------------------------------------------------------------

 rfddk=0; rfelfd=0; rfphon=0; rfmagn=0; rfstrs=0; rfuser=0; rf2_dkdk=0; rf2_dkde=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfddk',tread,'INT')
 if(tread==1) rfddk=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfelfd',tread,'INT')
 if(tread==1) rfelfd=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmagn',tread,'INT')
 if(tread==1) rfmagn=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfphon',tread,'INT')
 if(tread==1) rfphon=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfstrs',tread,'INT')
 if(tread==1) rfstrs=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfuser',tread,'INT')
 if(tread==1) rfuser=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkdk',tread,'INT')
 if(tread==1) rf2_dkdk=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rf2_dkde',tread,'INT')
 if(tread==1) rf2_dkde=intarr(1)

 response=0
 if(rfddk/=0.or.rf2_dkdk/=0.or.rf2_dkde/=0.or.rfelfd/=0.or.rfphon/=0.or.rfstrs/=0.or.rfuser/=0.or.rfmagn/=0)response=1

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optdriver',tread,'INT')
 if (tread==1) then
   dtset%optdriver=intarr(1)
 else
!  If optdriver was not read, while response=1, set optdriver to 1
   if(response==1)dtset%optdriver=1
 end if

!---------------------------------------------------------------------------
!For now, waiting express parallelisation for recursion
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tfkinfunc',tread,'INT')
 if(tread==1) dtset%tfkinfunc=intarr(1)

!---------------------------------------------------------------------------
! wvl_bigdft_comp, done here since default values of nline, nwfshist and iscf
! depend on its value (see indefo)
 if(dtset%usewvl==1) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wvl_bigdft_comp',tread,'INT')
   if(tread==1) dtset%wvl_bigdft_comp=intarr(1)
 end if

!---------------------------------------------------------------------------

 natom=dtset%natom
 npsp=dtset%npsp
 ntypat=dtset%ntypat

!No default value for znucl
 call intagm(dprarr,intarr,jdtset,marr,dtset%npsp,string(1:lenstr),'znucl',tread,'DPR')
 if(tread==1)then
   dtset%znucl(1:dtset%npsp)=dprarr(1:dtset%npsp)
 end if
 if(tread/=1)then
   write(message, '(3a)' )&
&   'The array znucl MUST be initialized in the input file while this is not done.',ch10,&
&   'Action: initialize znucl in your input file.'
   MSG_ERROR(message)
 end if

!The default for ratsph has already been initialized
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'ratsph',tread,'LEN')
 if(tread==1)then
   do ii=1,dtset%ntypat
     dtset%ratsph(ii)=dprarr(ii)
   end do
 end if
 ABI_ALLOCATE(ratsph,(dtset%ntypat))
 do ii=1,dtset%ntypat
   ratsph(ii)=dtset%ratsph(ii)
 end do

!Special treatment of _TYPAX (from a XYZ file), taking into account
!the fact that znucl does NOT depend on the dataset
!Examine all occurences of '_TYPAX'

 do
   index_typsymb=index(string(1:lenstr),'_TYPAX')
   if(index_typsymb==0)exit
!  Replace '_TYPAX' by '_TYPAT'
   string(index_typsymb:index_typsymb+5)='_TYPAT'
   index_upper=index_typsymb+5
!  Must start from the first blank after the tag (including possible dtset_char)
   index_upper=index(string(index_upper:lenstr),blank)+index_upper-1
   index_lower=index_upper

!  Examine all atoms (the end of the symbol string is delimited by a XX )
   do
     index_blank=index(string(index_upper:lenstr),blank)+index_upper-1
     string2=string(index_blank+1:index_blank+2)
     if(string2=="XX")exit
     found=0
!    Find the matching symbol
     do ipsp=1,dtset%npsp
       call atomdata_from_znucl(atom,dtset%znucl(ipsp))
       symbol = atom%symbol
       call inupper(symbol)
       call inupper(string2)

!      DEBUG
!      write(std_out,'(a)')' invars1 : before test, trim(adjustl(symbol)),trim(adjustl(string2))'
!      write(std_out,'(5a)' )'"',trim(adjustl(symbol)),'","',trim(adjustl(string2)),'"'
!      ENDDEBUG

       if(trim(adjustl(symbol))==trim(adjustl(string2)))then
         found=1
         index_upper=index_blank+1
!        Cannot deal properly with more that 9 psps
         if(ipsp>=10)then
           message = 'Need to use a pseudopotential with number larger than 9. Not allowed yet.'
           MSG_ERROR(message)
         end if

!        DEBUG
!        write(std_out,*)' invars1 : found ipsp=',ipsp
!        ENDDEBUG

         write(string1,'(i1)')ipsp
         string(index_lower:index_lower+1)=blank//string1
         index_lower=index_lower+2
       end if
     end do ! ipsp
!    if not found ...
     if(found==0)then
       write(message,'(6a)' )&
&       'Did not find matching pseudopotential for XYZ atomic symbol,',ch10,&
&       'with value ',string2,ch10,&
&       'Action: check that the atoms required by the XYZ file correspond to one psp file.'
       MSG_ERROR(message)
     end if
   end do ! Loop on atoms
!  One should find blanks after the last significant type value
   string(index_lower:index_blank+2)=blank
 end do ! loop to identify _TYPAX

!---------------------------------------------------------------------------

!Here, set up quantities that are related to geometrical description
!of the system (acell,rprim,xred), as well as
!initial velocity(vel), and spin of atoms (spinat), nuclear dipole moments
! of atoms (nucdipmom),
!the symmetries (symrel,symafm, and tnons)
!and the list of fixed atoms (iatfix,iatfixx,iatfixy,iatfixz).
!Arrays have already been
!dimensioned thanks to the knowledge of msym and mxnatom

!ji: We need to read the electric field before calling ingeo
!****** Temporary ******

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'berryopt',tread,'INT')
 if(tread==1) dtset%berryopt=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'berrysav',tread,'INT')
 if(tread==1) dtset%berrysav=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'bfield',tread,'DPR')
 if (tread==1) dtset%bfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'dfield',tread,'DPR')
 if (tread==1) dtset%dfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'efield',tread,'DPR')
 if (tread==1) dtset%efield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_dfield',tread,'DPR')
 if (tread==1) dtset%red_dfield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_efield',tread,'DPR')
 if (tread==1) dtset%red_efield(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'red_efieldbar',tread,'DPR')
 if (tread==1) dtset%red_efieldbar(1:3) = dprarr(1:3)

 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'jfielddir',tread,'INT')
 if(tread==1) dtset%jfielddir(1:3)=intarr(1:3)

!We need to know nsppol/nspinor/nspden before calling ingeo
 nsppol=dtset%nsppol
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsppol',tread,'INT')
 if(tread==1) nsppol=intarr(1)

!Alternate SIESTA definition of nsppol
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'SpinPolarized',tread_alt,'LOG')
 if(tread_alt==1)then
   if(tread==1)then
     write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&     ' invars1: ERROR -',ch10,&
&     '  nsppol and SpinPolarized cannot be specified simultaneously',ch10,&
&     '  for the same dataset.',ch10,&
&     '  Action : check the input file.'
     call wrtout(std_out,  message,'COLL')
     leave=1
   else
!    Note that SpinPolarized is a logical input variable
     nsppol=1
     if(intarr(1)==1)nsppol=2
     tread=1
   end if
 end if
 dtset%nsppol=nsppol

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nspinor',tread,'INT')
 if(tread==1) dtset%nspinor=intarr(1)

!Has to read pawspnorb now, in order to adjust nspinor
!Also, if nspinor=1, turn on spin-orbit coupling by default, here for the PAW case. NC case is treated elsewhere.
 if (dtset%usepaw>0)then
!  Change the default value
   if(dtset%nspinor==2)dtset%pawspnorb=1
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawspnorb',tread,'INT')
   if(tread==1)then
     dtset%pawspnorb=intarr(1)
     if(dtset%pawspnorb>0) dtset%nspinor=2
   else
     if(dtset%nspinor==2)then
       write(message, '(4a)' ) ch10,&
&       ' invars1: COMMENT -',ch10,&
&       '  With nspinor=2 and usepaw=1, pawspnorb=1 has been switched on by default.'
       call wrtout(iout,  message,'COLL')
     end if
   end if
 end if
 nspinor=dtset%nspinor

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nspden',tread,'INT')
 if(tread==1) then
   dtset%nspden=intarr(1)
 else
   dtset%nspden=dtset%nsppol
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntypalch',tread,'INT')
 if(tread==1) dtset%ntypalch=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nzchempot',tread,'INT')
 if(tread==1) dtset%nzchempot=intarr(1)

 ntypalch=dtset%ntypalch
 if(ntypalch>ntypat)then
   write(message, '(3a,i0,a,i0,a,a)' )&
&   'The input variable ntypalch must be smaller than ntypat, while it is',ch10,&
&   'ntypalch=',dtset%ntypalch,', and ntypat=',ntypat,ch10,&
&   'Action: check ntypalch vs ntypat in your input file.'
   MSG_ERROR(message)
 end if

 ntyppure=ntypat-ntypalch
 dtset%ntyppure=ntyppure
 npspalch=npsp-ntyppure
 dtset%npspalch=npspalch
 if(npspalch<0)then
   write(message, '(a,i0,2a,i0,a,a)' )&
&   'The number of available pseudopotentials, npsp=',npsp,ch10,&
&   'is smaller than the requested number of types of pure atoms, ntyppure=',ntyppure,ch10,&
&   'Action: check ntypalch versus ntypat and npsp in your input file.'
   MSG_ERROR(message)
 end if

 if(ntypalch>0)then
   call intagm(dprarr,intarr,jdtset,marr,ntypalch,string(1:lenstr),'algalch',tread,'INT')
   if(tread==1) dtset%algalch(1:ntypalch)=intarr(1:ntypalch)
 end if

!Read the Zeeman field
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'zeemanfield',tread,'BFI')
 if(tread==1) then
   if(dtset%nspden == 2)then
     write(message,'(7a)')&
&     'A Zeeman field has been specified without noncollinear spins.',ch10,&
&     'Only the z-component of the magnetic field will be used.'
     MSG_WARNING(message)
   else if (dtset%nspden == 1)then
     write(message, '(a,a,a)' )&
&     'A Zeeman field has been specified for a non-spin-polarized calculation.',ch10,&
&     'Action: check the input file.'
     MSG_ERROR(message)
   end if

   dtset%zeemanfield(1:3) = dprarr(1:3)
 end if


 ABI_ALLOCATE(amu,(ntypat))
 ABI_ALLOCATE(mixalch,(npspalch,ntypalch))
 ABI_ALLOCATE(vel,(3,natom))
 ABI_ALLOCATE(vel_cell,(3,3))
 ABI_ALLOCATE(xred,(3,natom))
 intimage=2 ; if(dtset%nimage==1)intimage=1
 do ii=1,dtset%nimage+1
   iimage=ii
   if(dtset%nimage==1 .and. ii==2)exit
   if(dtset%nimage==2 .and. ii==3)exit
   if(dtset%nimage> 2 .and. ii==intimage)cycle ! Will do the intermediate reference image at the last reading
   if(dtset%nimage>=2 .and. ii==dtset%nimage+1)iimage=intimage

   write(message,'(a,i0)')' invars1 : treat image number: ',iimage
   call wrtout(std_out,message,'COLL')

!  Need to reset nsym to default value for each image
   dtset%nsym=0

!  Call ingeo for each image in turn, with the possible default values
   acell=dtset%acell_orig(1:3,iimage)
   amu=dtset%amu_orig(1:ntypat,iimage)
   mixalch=dtset%mixalch_orig(1:npspalch,1:ntypalch,iimage)
   rprim=dtset%rprim_orig(1:3,1:3,iimage)
   vel=dtset%vel_orig(1:3,1:natom,iimage)
   vel_cell=dtset%vel_cell_orig(1:3,1:3,iimage)
   xred=dtset%xred_orig(1:3,1:natom,iimage)
   ABI_ALLOCATE(iatfix,(3,natom))
   ABI_ALLOCATE(nucdipmom,(3,natom))
   ABI_ALLOCATE(spinat,(3,natom))
   ABI_ALLOCATE(typat,(natom))
   ABI_ALLOCATE(znucl,(dtset%npsp))
   nucdipmom(1:3,1:natom)=dtset%nucdipmom(1:3,1:natom)
   spinat(1:3,1:natom)=dtset%spinat(1:3,1:natom)
   znucl(1:dtset%npsp)=dtset%znucl(1:dtset%npsp)

   call ingeo(acell,amu,dtset,bravais,dtset%genafm(1:3),iatfix,&
&   dtset%icoulomb,iimage,iout,jdtset,dtset%jellslab,lenstr,mixalch,&
&   msym,natom,dtset%nimage,dtset%npsp,npspalch,dtset%nspden,dtset%nsppol,&
&   dtset%nsym,ntypalch,dtset%ntypat,nucdipmom,dtset%nzchempot,&
&   dtset%pawspnorb,dtset%ptgroupma,ratsph,&
&   rprim,dtset%slabzbeg,dtset%slabzend,dtset%spgroup,spinat,&
&   string,symafm,dtset%symmorphi,symrel,tnons,dtset%tolsym,typat,vel,vel_cell,xred,znucl)
   dtset%iatfix(1:3,1:natom)=iatfix(1:3,1:natom)
   dtset%nucdipmom(1:3,1:natom)=nucdipmom(1:3,1:natom)
   dtset%spinat(1:3,1:natom)=spinat(1:3,1:natom)
   dtset%typat(1:natom)=typat(1:natom)
   ABI_DEALLOCATE(iatfix)
   ABI_DEALLOCATE(nucdipmom)
   ABI_DEALLOCATE(spinat)
   ABI_DEALLOCATE(typat)
   ABI_DEALLOCATE(znucl)
   dtset%acell_orig(1:3,iimage)=acell
   dtset%amu_orig(1:ntypat,iimage)=amu
   dtset%mixalch_orig(1:npspalch,1:ntypalch,iimage)=mixalch
   dtset%rprim_orig(1:3,1:3,iimage)=rprim
   dtset%vel_orig(1:3,1:natom,iimage)=vel
   dtset%vel_cell_orig(1:3,1:3,iimage)=vel_cell
   dtset%xred_orig(1:3,1:natom,iimage)=xred
   call mkrdim(dtset%acell_orig(1:3,iimage),dtset%rprim_orig(1:3,1:3,iimage),dtset%rprimd_orig(1:3,1:3,iimage))
 end do
 ABI_DEALLOCATE(amu)
 ABI_DEALLOCATE(mixalch)
 ABI_DEALLOCATE(vel)
 ABI_DEALLOCATE(vel_cell)
 ABI_DEALLOCATE(xred)

!Examine whether there is some vacuum space in the unit cell
 call invacuum(jdtset,lenstr,natom,dtset%rprimd_orig(1:3,1:3,intimage),string,vacuum,&
& dtset%xred_orig(1:3,1:natom,intimage))

!DEBUG
!write(std_out,*)' invars1: before inkpts, dtset%mixalch_orig(1:npspalch,1:ntypalch,:)=',&
!dtset%mixalch_orig(1:npspalch,1:ntypalch,1:dtset%nimage)
!stop
!ENDDEBUG

!---------------------------------------------------------------------------

!Set up k point grid number
!First, get additional information
 dtset%kptopt=1
 if(dtset%nspden==4)dtset%kptopt=4
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'kptopt',tread,'INT')
 if(tread==1) dtset%kptopt=intarr(1)

 dtset%qptopt=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qptopt',tread,'INT')
 if(tread==1) dtset%qptopt=intarr(1)

 iscf=5
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iscf',tread,'INT')
 if(tread==1) iscf=intarr(1)


 dtset%natsph=dtset%natom
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natsph',tread,'INT')
 if(tread==1) dtset%natsph=intarr(1)

 dtset%natsph_extra=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natsph_extra',tread,'INT')
 if(tread==1) dtset%natsph_extra=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natvshift',tread,'INT')
 if(tread==1) dtset%natvshift=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nconeq',tread,'INT')
 if(tread==1) dtset%nconeq=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nkptgw',tread,'INT')
 if(tread==1) dtset%nkptgw=intarr(1)
 if (dtset%nkptgw<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input nkptgw must be >= 0, but was ',dtset%nkptgw,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

!Number of points for long wavelength limit.
!Default is dtset%gw_nqlwl=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_nqlwl',tread,'INT')
 if(tread==1) dtset%gw_nqlwl=intarr(1)
 if (dtset%gw_nqlwl<0) then
   write(message, '(a,i12,a,a,a,a)' )&
&   'Input gw_nqlwl must be > 0, but was ',dtset%gw_nqlwl,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

!Read number of k-points (if specified)
 nkpt=0
 if(dtset%kptopt==0)nkpt=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nkpt',tread,'INT')
 if(tread==1) nkpt=intarr(1)
 dtset%nkpt=nkpt

 call chkint_ge(0,0,cond_string,cond_values,ierr,'nkpt',nkpt,0,iout)
 if(dtset%kptopt==0)then
   cond_string(1)='kptopt' ; cond_values(1)=0
   call chkint_ge(1,1,cond_string,cond_values,ierr,'nkpt',nkpt,1,iout)
 end if

 nkpthf=nkpt
 dtset%nkpthf=nkpt

!Will compute the actual value of nkpt, if needed. Otherwise,
!test that the value of nkpt is OK, if kptopt/=0
!Set up dummy arrays istwfk, kpt, wtk

 if(nkpt/=0 .or. dtset%kptopt/=0)then
   ABI_ALLOCATE(istwfk,(nkpt))
   ABI_ALLOCATE(kpt,(3,nkpt))
   ABI_ALLOCATE(kpthf,(3,nkpthf))
   ABI_ALLOCATE(wtk,(nkpt))
!  Here, occopt is also a dummy argument
   occopt=1 ; dtset%nshiftk=1 ; dtset%kptrlatt(:,:)=0

   kptrlen=20.0_dp ; wtk(:)=1.0_dp
   dtset%shiftk(:,:)=half

   nqpt=0
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqpt',tread,'INT')
   if(tread==1) nqpt=intarr(1)

   chksymbreak=1
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chksymbreak',tread,'INT')
   if(tread==1) chksymbreak=intarr(1)


!  Use the first image to predict k and/or q points,
!  except if an intermediate image is available
   intimage=1 ; if(dtset%nimage>2)intimage=(1+dtset%nimage)/2

!  Find the q-point, if any.
   if(nqpt==1)then
     call inqpt(chksymbreak,std_out,jdtset,lenstr,msym,natom,dtset%qptn,dtset%wtq,&
&     dtset%rprimd_orig(1:3,1:3,intimage),dtset%spinat,string,dtset%typat,&
&     vacuum,dtset%xred_orig(1:3,1:natom,intimage),dtset%qptrlatt)
   end if

!  Find the k point grid
   call inkpts(bravais,chksymbreak,dtset%fockdownsampling,iout,iscf,istwfk,jdtset,&
&   kpt,kpthf,dtset%kptopt,kptnrm,dtset%kptrlatt_orig,dtset%kptrlatt,kptrlen,lenstr,msym,&
&   nkpt,nkpthf,nqpt,dtset%ngkpt,dtset%nshiftk,dtset%nshiftk_orig,dtset%shiftk_orig,dtset%nsym,&
&   occopt,dtset%qptn,response,dtset%rprimd_orig(1:3,1:3,intimage),dtset%shiftk,&
&   string,symafm,symrel,vacuum,wtk)

   ABI_DEALLOCATE(istwfk)
   ABI_DEALLOCATE(kpt)
   ABI_DEALLOCATE(kpthf)
   ABI_DEALLOCATE(wtk)
!  nkpt and nkpthf have been computed, as well as the k point grid, if needed
   dtset%nkpt=nkpt
   dtset%nkpthf=nkpthf
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqptdm',tread,'INT')
 if(tread==1) dtset%nqptdm=intarr(1)

 if (dtset%nqptdm<-1) then
   write(message, '(a,i12,a,a,a,a)' )&
&   'Input nqptdm must be >= 0, but was ',dtset%nqptdm,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nzchempot',tread,'INT')
 if(tread==1) dtset%nzchempot=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cd_customnimfrqs',tread,'INT')
 if(tread==1) dtset%cd_customnimfrqs=intarr(1)

 if (dtset%cd_customnimfrqs<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input cd_customnimfrqs must be >= 0, but was ',dtset%cd_customnimfrqs,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gw_customnfreqsp',tread,'INT')
 if(tread==1) dtset%gw_customnfreqsp=intarr(1)

 if (dtset%gw_customnfreqsp<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input gw_customnfreqsp must be >= 0, but was ',dtset%gw_customnfreqsp,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gwls_n_proj_freq',tread,'INT')
 if(tread==1) dtset%gwls_n_proj_freq=intarr(1)

 if (dtset%gwls_n_proj_freq<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input gwls_n_proj_freq must be >= 0, but was ',dtset%gwls_n_proj_freq,ch10,&
&   'This is not allowed.',ch10,&
&   'Action : check the input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_calc_dirs',tread,'INT')
 if(tread==1) dtset%efmas_calc_dirs=intarr(1)

 if (ABS(dtset%efmas_calc_dirs)>3) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input efmas_calc_dirs must be between -3 and 3, but was ',dtset%efmas_calc_dirs,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'efmas_n_dirs',tread,'INT')
 if(tread==1) dtset%efmas_n_dirs=intarr(1)

 if (dtset%efmas_n_dirs<0) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'Input efmas_n_dirs must be >= 0, but was ',dtset%efmas_n_dirs,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if
!---------------------------------------------------------------------------

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnos',tread,'INT')
 if(tread==1) dtset%nnos=intarr(1)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ga_n_rules',tread,'INT')
 if(tread==1) dtset%ga_n_rules=intarr(1)

!Perform the first checks

 leave=0

!Check that nkpt is greater than 0
 if (nkpt<=0) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  After inkpts, nkpt must be > 0, but was ',nkpt,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,  message,'COLL')
   leave=1
 end if

!Check that nsppol is 1 or 2
 if (nsppol/=1 .and. nsppol/=2) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nsppol must be 1 or 2, but was ',nsppol,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Check that nspinor is 1 or 2
 if (nspinor/=1 .and. nspinor/=2) then
   write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  Input nspinor must be 1 or 2, but was ',nspinor,ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Check that nspinor and nsppol are not 2 together
 if (nsppol==2 .and. nspinor==2) then
   write(message, '(8a)' ) ch10,&
&   ' invars1: ERROR -',ch10,&
&   '  nspinor and nsappol cannot be 2 together !',ch10,&
&   '  This is not allowed.',ch10,&
&   '  Action : check the input file.'
   call wrtout(std_out,message,'COLL')
   leave=1
 end if

!Here, leave if an error has been detected earlier
 if(leave==1) then
   message = ' Other errors might be present in the input file. '
   MSG_ERROR(message)
 end if


!Now, take care of mband_upper

 mband_upper=1
 occopt=1
 fband=0.5_dp

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'occopt',tread,'INT')
 if(tread==1) occopt=intarr(1)

!Also read fband, that is an alternative to nband. The default
!is different for occopt==1 and for metallic occupations.
 if(occopt==1)fband=0.125_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fband',tfband,'DPR')
 if(tfband==1)fband=dprarr(1)

!fband cannot be used when occopt==0 or occopt==2
 if(tfband==1 .and. (occopt==0 .or. occopt==2) )then
   write(message, '(3a)' )&
&   'fband cannot be used if occopt==0 or occopt==2 ',ch10,&
&   'Action: correct your input file, suppress fband, or change occopt.'
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(nband,(nkpt*nsppol))
 tnband=0

!Compute ziontypat
!When the pseudo-atom is pure, simple copy
 if(ntyppure>0)then
   do itypat=1,ntyppure
     dtset%ziontypat(itypat)=zionpsp(itypat)
   end do
 end if
!When the pseudo-atom is alchemical, must make mixing
 if(ntypalch>0)then
   do itypat=ntyppure+1,ntypat
     dtset%ziontypat(itypat)=zero
     do ipsp=ntyppure+1,npsp
       dtset%ziontypat(itypat)=dtset%ziontypat(itypat) &
&       +dtset%mixalch_orig(ipsp-ntyppure,itypat-ntyppure,1)*zionpsp(ipsp)
     end do
   end do
 end if



 if (occopt==0 .or. occopt==1 .or. (occopt>=3 .and. occopt<=8) ) then
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nband',tnband,'INT')
!  Note : mband_upper is initialized, not nband
   if(tnband==1) mband_upper=intarr(1)

   if(tfband==1 .and. tnband==1)then
     write(message, '(3a)' )&
&     'fband and nband cannot be used together. ',ch10,&
&     'Action: correct your input file, suppress either fband or nband.'
     MSG_ERROR(message)
   end if

!  In case nband was not read, use fband, either read, or the default,
!  to provide an upper limit for mband_upper
   if(tnband==0)then

     charge=0.0_dp
     call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'charge',tread,'DPR')
     if(tread==1) charge=dprarr(1)

!    Only take into account negative charge, to compute maximum number of bands
     if(charge > 0.0_dp)charge=0.0_dp

!     mband_upper=nspinor*((nint(zion_max)*natom+1)/2 - floor(charge/2.0_dp)&
!&     + ceiling(fband*natom-1.0d-10))
     zval=0.0_dp
     do iatom=1,natom
       zval=zval+dtset%ziontypat(dtset%typat(iatom))
     end do
     zelect=zval-charge
     mband_upper=nspinor * ((ceiling(zelect-1.0d-10)+1)/2 + ceiling( fband*natom - 1.0d-10 ))

     nband(:)=mband_upper

!    DEBUG
!    write(std_out,*)' invars1 : zion_max,natom,fband,mband_upper '
!    write(std_out,*)zion_max,natom,fband,mband_upper
!    ENDDEBUG

   end if

   nband(:)=mband_upper

 else if (occopt==2) then
   ABI_ALLOCATE(reaalloc,(nkpt*nsppol))
   call intagm(reaalloc,nband,jdtset,nkpt*nsppol,nkpt*nsppol,string(1:lenstr),'nband',tnband,'INT')
   if(tnband==1)then
     do ikpt=1,nkpt*nsppol
       if (nband(ikpt)>mband_upper) mband_upper=nband(ikpt)
     end do
   end if
   ABI_DEALLOCATE(reaalloc)
 else
   write(message, '(a,i0,3a)' )&
&   'occopt=',occopt,' is not an allowed value.',ch10,&
&   'Action: correct your input file.'
   MSG_ERROR(message)
 end if

!Check that mband_upper is greater than 0
 if (mband_upper<=0) then
   write(message, '(a,i0,4a)' )&
&   'Maximal nband must be > 0, but was ',mband_upper,ch10,&
&   'This is not allowed.',ch10,&
&   'Action: check the input file.'
   MSG_ERROR(message)
 end if

!The following 3 values are needed to dimension the parallelism over images
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'imgmov',tread,'INT')
 if(tread==1) dtset%imgmov=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntimimage',tread,'INT')
 if(tread==1) dtset%ntimimage=intarr(1)
 call intagm(dprarr,intarr,jdtset,marr,dtset%nimage,string(1:lenstr),'dynimage',tread,'INT')
 if(tread==1)then
   dtset%dynimage(1:dtset%nimage)=intarr(1:dtset%nimage)
 else if (dtset%imgmov==2.or.dtset%imgmov==5) then
   dtset%dynimage(1)=0;dtset%dynimage(dtset%nimage)=0
 end if
 dtset%ndynimage=count(dtset%dynimage(1:dtset%nimage)/=0)

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wfoptalg',tread,'INT')
 if(tread==1) then
   dtset%wfoptalg=intarr(1)
 else
   if (dtset%usepaw==0)    dtset%wfoptalg=0
   if (dtset%usepaw/=0)    dtset%wfoptalg=10
   if (dtset%optdriver==RUNL_GSTATE) then
     if (dtset%paral_kgb/=0) dtset%wfoptalg=14
   end if
 end if

!---------------------------------------------------------------------------
!Some PAW+DMFT keywords
 dtset%usedmft=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usedmft',tread,'INT')
 if(tread==1) dtset%usedmft=intarr(1)

!Some ucrpa keywords
 dtset%ucrpa=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ucrpa',tread,'INT')
 if(tread==1) dtset%ucrpa=intarr(1)

 if (dtset%ucrpa > 0 .and. dtset%usedmft > 0) then
   write(message, '(7a)' )&
&   'usedmft and ucrpa are both activated in the input file ',ch10,&
&   'In the following, abinit assume you are doing a ucrpa calculation and ',ch10,&
&   'you define Wannier functions as in DFT+DMFT calculation',ch10,&
&   'If instead, you want to do a full dft+dmft calculation and not only the Wannier construction, use ucrpa=0'
   MSG_WARNING(message)
 end if

!Some PAW+U keywords
 dtset%usepawu=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usepawu',tread,'INT')
 if(tread==1) dtset%usepawu=intarr(1)

 if (dtset%usepawu > 0 .and. dtset%usedmft > 0) then
   write(message, '(7a)' )&
&   'usedmft and usepawu are both activated ',ch10,&
&   'This is not an usual calculation:',ch10,&
&   'usepawu will be put to 10:',ch10,&
&   'LDA+U potential and energy will be put to zero'
   MSG_WARNING(message)
 end if
 if ( dtset%usedmft > 0 .and. dtset%usepawu >= 0 ) dtset%usepawu = 1

 dtset%usedmatpu=0
 dtset%lpawu(1:dtset%ntypat)=-1
 if (dtset%usepawu>0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'lpawu',tread,'INT')
   if(tread==1) dtset%lpawu(1:dtset%ntypat)=intarr(1:dtset%ntypat)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usedmatpu',tread,'INT')
   if(tread==1) dtset%usedmatpu=intarr(1)

 end if

!Some PAW+Exact exchange keywords
 dtset%useexexch=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'useexexch',tread,'INT')
 if(tread==1) dtset%useexexch=intarr(1)


 dtset%lexexch(1:dtset%ntypat)=-1

 if (dtset%useexexch>0) then
   call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'lexexch',tread,'INT')
   if(tread==1) dtset%lexexch(1:dtset%ntypat)=intarr(1:dtset%ntypat)
 end if
! LDA minus half keyword
 call intagm(dprarr,intarr,jdtset,marr,dtset%ntypat,string(1:lenstr),'ldaminushalf',tread,'INT')
 if(tread==1) dtset%ldaminushalf(1:dtset%ntypat)=intarr(1:dtset%ntypat)
!Some plowan data
 dtset%plowan_natom=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_natom',tread,'INT')
 if(tread==1) dtset%plowan_natom=intarr(1)

 dtset%plowan_nt=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'plowan_nt',tread,'INT')
 if(tread==1) dtset%plowan_natom=intarr(1)


!PAW potential zero keyword
 dtset%usepotzero=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'usepotzero',tread,'INT')
 if(tread==1) dtset%usepotzero=intarr(1)

!Macro_uj (determination of U in PAW+U), governs also allocation of atvshift
 dtset%macro_uj = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'macro_uj',tread,'INT')
 if(tread==1) dtset%macro_uj=intarr(1)

 ABI_DEALLOCATE(nband)
 ABI_DEALLOCATE(ratsph)
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)

end subroutine invars1
!!***
