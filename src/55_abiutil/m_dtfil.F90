!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dtfil
!! NAME
!!  m_dtfil
!!
!! FUNCTION
!!   object and procedures dealing with input/output filenames
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (XG, MT)
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

module m_dtfil

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_xmpi
 use m_build_info

 use m_clib,         only : clib_rename
 use m_fstrings,     only : int2char4
 use m_io_tools,     only : open_file
 use m_libpaw_tools, only : libpaw_log_flag_set

 implicit none

 private
!!***

 public :: dtfil_init
 public :: dtfil_init_img
 public :: dtfil_init_time
 public :: mkfilename
 public :: isfile
 public :: iofn1
 ! MG FIXME: Deprecated
 public :: status
!!***

contains
!!***

!!****f* m_dtfil/dtfil_init
!!
!! NAME
!! dtfil_init
!!
!! FUNCTION
!! Initialize most of the dtfil structured variable
!! (what is left should be initialized inside the itimimage,
!! iimage and itime loops).
!!
!! INPUTS
!! dtset=<type datasets_type>contain all input variables for the current dataset
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! idtset=number of the dataset
!! image_index= -optional argument- index of image to be used when appending
!!             "_IMGxxx" string to file names. To be used only when an algorithm
!!             using images of the cell is activated
!! jdtset_(0:ndtset)=actual index of the datasets
!! mpi_enreg=information about MPI parallelization
!! ndtset=number of datasets
!!
!! OUTPUT
!! dtfil=<type datafiles_type>infos about file names, file unit numbers
!!  (part of which were initialized previously)
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! PARENTS
!!      driver,gstateimg
!!
!! CHILDREN
!!      appdig,int2char4,mkfilename
!!
!! SOURCE

subroutine dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset,&
&                      image_index) ! optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: idtset,ndtset
 integer, optional, intent(in) :: image_index
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil !vz_i
!arrays
 integer :: jdtset_(0:ndtset)
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
!scalars
! Define input and output unit numbers (do not forget, unit 5 and 6 are standard input and output)
! Also, unit number 21, 22 and 23 are used in dfpt_nstdy, for the 3 dot wavefunctions.
! Unit 50,51,52 and 53 are used in dfpt_looppert (for ipert=natom+2, ipert=natom+10 and ipert=natom+11).
! Others unit numbers will be used in the case of the variational and 2n+1 expressions.
! In defs_basis, one defines :
!   std_in=5, ab_in=5, std_out=6, ab_out=7, tmp_unit=9, tmp_unit2=10
 integer,parameter :: unchi0=42,unddb=16,unddk=50,undkdk=54,undkde=55,unkg1=19,unkg=17,unkgq=18
 integer,parameter :: unpaw=26,unpaw1=27,unpawq=28,unpos=30
 integer,parameter :: unwff1=1,unwff2=2,unwff3=8,unwffgs=3,unwfkq=4,unwft1=11
 integer,parameter :: unwft2=12,unwft3=15,unwftgs=13,unwftkq=14,unylm=24,unylm1=25
 integer,parameter :: unkss=40,unscr=41,unqps=43
 integer :: ii,iimage,ireadden,ireadkden,ireadwf,ixx,jdtset,will_read
 character(len=10) :: appen,tag
 character(len=9) :: stringvar
 character(len=15) :: stringfile
 character(len=500) :: message
 character(len=fnlen) :: filsus,filddbsin,fildens1in,fildensin,filpawdensin,filkdensin,filqps,filscr,fil_efmas
 character(len=fnlen) :: fnamewff1,fnamewffddk,fnamewffdelfd,fnamewffdkdk,fnamewffdkde,fnamewffk,fnamewffq
 character(len=fnlen) :: filbseig,filfft,filhaydock,fil_bsreso,fil_bscoup
 character(len=fnlen) :: filwfkfine
 character(len=fnlen) :: filnam_ds(5)
 character(len=fnlen) :: tmpfil(14)
 integer :: idtmpfil(14)

!******************************************************************

 DBG_ENTER("COLL")

 iimage=0;if (present(image_index)) iimage=image_index

 dtfil%unchi0 =unchi0
 dtfil%unddb  =unddb
 dtfil%unddk  =unddk
 dtfil%undkde =undkde
 dtfil%undkdk =undkdk
 dtfil%unkg   =unkg
 dtfil%unkgq  =unkgq
 dtfil%unkg1  =unkg1
 dtfil%unkss  =unkss
 dtfil%unqps  =unqps
 dtfil%unscr  =unscr
 dtfil%unwff1 =unwff1
 dtfil%unwff2 =unwff2
 dtfil%unwff3 =unwff3
 dtfil%unwffgs=unwffgs
 dtfil%unwffkq=unwfkq
 dtfil%unwft1 =unwft1
 dtfil%unwft2 =unwft2
 dtfil%unwft3 =unwft3
 dtfil%unwftgs=unwftgs
 dtfil%unwftkq=unwftkq
 dtfil%unylm  =unylm
 dtfil%unylm1 =unylm1
 dtfil%unpaw  =unpaw
 dtfil%unpaw1 =unpaw1
 dtfil%unpawq =unpawq
 dtfil%unpos  =unpos
 filnam_ds(1:5)=filnam(1:5)
 jdtset=dtset%jdtset

!If multi dataset mode, special treatment of filenames 3 and 4 (density and
!wavefunctions input and output, as well as other output files)
 if(ndtset>0)then
   call appdig(jdtset,'',appen)
   filnam_ds(3)=trim(filnam(3))//'_DS'//trim(appen)
   filnam_ds(4)=trim(filnam(4))//'_DS'//trim(appen)
 end if

!If multi image mode (nimage>1), special treatment of filenames 4 and 5
 if(iimage>0)then
   call appdig(iimage,'',appen)
   filnam_ds(4)=trim(filnam_ds(4))//'_IMG'//trim(appen)
   filnam_ds(5)=trim(filnam_ds(5))//'_IMG'//trim(appen)
 end if

!According to getwfk and irdwfk, build _WFK file name, referred as fnamewffk
 if (iimage>0.and.dtfil%getwfk_from_image/=0) then
   if (dtfil%getwfk_from_image==-1) then
     call appdig(iimage,'',appen)
   else
     call appdig(dtfil%getwfk_from_image,'',appen)
   end if
   stringfile='_IMG'//trim(appen)//'_WFK'
 else
   stringfile='_WFK'
 end if
 stringvar='wfk'
 call mkfilename(filnam,fnamewffk,dtset%getwfk,idtset,dtset%irdwfk,jdtset_,&
& ndtset,stringfile,stringvar,will_read)

 if(dtset%optdriver/=RUNL_RESPFN)ireadwf=will_read
 if(ndtset/=0 .and. dtset%optdriver==RUNL_RESPFN .and. will_read==0)then
   write(message, '(5a,i3,3a,i3,a,i3,3a)' )&
&   'At least one of the input variables irdwfk and getwfk ',ch10,&
&   'must refer to a valid _WFK file, in the response function',ch10,&
&   'case, while for idtset=',idtset,',',ch10,&
&   'they are irdwfk=',dtset%irdwfk,', and getwfk=',dtset%getwfk,'.',ch10,&
&   'Action: correct irdwfk or getwfk in your input file.'
   MSG_ERROR(message)
 end if

!Treatment of the other get wavefunction variable, if response function case or nonlinear case
 if ( ANY(dtset%optdriver == (/RUNL_RESPFN, RUNL_NONLINEAR, RUNL_EPH/)) ) then

!  According to getwfq and irdwfq, build _WFQ file name, referred as fnamewffq
   stringfile='_WFQ' ; stringvar='wfq'
   call mkfilename(filnam,fnamewffq,dtset%getwfq,idtset,dtset%irdwfq,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
!  If fnamewffq is not initialized thanks to getwfq or irdwfq, use fnamewffk
   if(will_read==0)fnamewffq=fnamewffk

!  According to get1wf and ird1wf, build _1WF file name, referred as fnamewff1
   stringfile='_1WF' ; stringvar='1wf'
   call mkfilename(filnam,fnamewff1,dtset%get1wf,idtset,dtset%ird1wf,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   ireadwf=will_read

!  According to getddk and irdddk, build _1WF file name, referred as fnamewffddk
   stringfile='_1WF' ; stringvar='ddk'
   call mkfilename(filnam,fnamewffddk,dtset%getddk,idtset,dtset%irdddk,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

!  According to getdelfd, build _1WF file name, referred as fnamewffdelfd
   stringfile='_1WF' ; stringvar='delfd'
   call mkfilename(filnam,fnamewffdelfd,dtset%getdelfd,idtset,0,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

!  According to getdkdk, build _1WF file name, referred as fnamewffdkdk
   stringfile='_1WF' ; stringvar='dkdk'
   call mkfilename(filnam,fnamewffdkdk,dtset%getdkdk,idtset,0,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

!  According to getdkde, build _1WF file name, referred as fnamewffdkde
   stringfile='_1WF' ; stringvar='dkde'
   call mkfilename(filnam,fnamewffdkde,dtset%getdkde,idtset,0,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)

 end if

!-------------------------------------------------------------------------------------------
!Build name of files from dtfil%filnam_ds(3)

!SP :According to getddb, build _DDB file name, referred as filddbsin
 stringfile='_DDB'
 stringvar='ddb'
 call mkfilename(filnam,filddbsin,dtset%getddb,idtset,dtset%irdddb,jdtset_,&
& ndtset,stringfile,stringvar,will_read)

!According to getden, build _DEN file name, referred as fildensin
!A default is available if getden is 0
 if (iimage>0.and.dtfil%getden_from_image/=0) then
   if (dtfil%getden_from_image==-1) then
     call appdig(iimage,'',appen)
   else
     call appdig(dtfil%getden_from_image,'',appen)
   end if
   stringfile='_IMG'//trim(appen)//'_DEN'
 else
   stringfile='_DEN'
 end if
 stringvar='den'
 call mkfilename(filnam,fildensin,dtset%getden,idtset,dtset%irdden,jdtset_,&
& ndtset,stringfile,stringvar,will_read)

 if(will_read==0)fildensin=trim(filnam_ds(3))//'_DEN'
 ireadden=will_read

 if ((dtset%optdriver==RUNL_GWLS.or.dtset%optdriver==RUNL_GSTATE) &
& .and.dtset%iscf<0) ireadden=1
!if (optdriver==RUNL_GSTATE.and.ireadwf/=0) ireadden=0

!According to getpawden, build _PAWDEN file name, referred as filpawdensin
!A default is available if getden is 0
 if (iimage>0.and.dtfil%getpawden_from_image/=0) then
   if (dtfil%getpawden_from_image==-1) then
     call appdig(iimage,'',appen)
   else
     call appdig(dtfil%getpawden_from_image,'',appen)
   end if
   stringfile='_IMG'//trim(appen)//'_PAWDEN'
 else
   stringfile='_PAWDEN'
 end if
 stringvar='pawden'
 call mkfilename(filnam,filpawdensin,dtset%getpawden,idtset,dtset%irdden,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filpawdensin=trim(filnam_ds(3))//'_PAWDEN'

!According to getden and usekden, build _KDEN file name, referred as filkdensin
!A default is available if getden is 0
 if(dtset%usekden==1)then
   if (iimage>0.and.dtfil%getden_from_image/=0) then
     if (dtfil%getden_from_image==-1) then
       call appdig(iimage,'',appen)
     else
       call appdig(dtfil%getden_from_image,'',appen)
     end if
     stringfile='_IMG'//trim(appen)//'_KDEN'
   else
     stringfile='_KDEN'
   end if
   stringvar='kden'
   call mkfilename(filnam,filkdensin,dtset%getden,idtset,dtset%irdden,jdtset_,&
&   ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filkdensin=trim(filnam_ds(3))//'_KDEN'
   ireadkden=will_read
   if ((dtset%optdriver==RUNL_GSTATE.or.dtset%optdriver==RUNL_GWLS).and.dtset%iscf<0) ireadkden=1
 else
   ireadkden=0
 end if

!According to get1den, build _DEN file name, referred as fildens1in
!A default is available if get1den is 0
 stringfile='_DEN' ; stringvar='1den'
 call mkfilename(filnam,fildens1in,dtset%get1den,idtset,dtset%ird1den,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fildens1in=trim(filnam_ds(3))//'_DEN'

!According to getefmas and irdefmas, build _EFMAS file name, referred as fil_efmas
!A default is available if getefmas is 0
 stringfile='_EFMAS.nc' ; stringvar='efmas'
 call mkfilename(filnam,fil_efmas,dtset%getefmas,idtset,dtset%irdefmas,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fil_efmas=trim(filnam_ds(3))//'_EFMAS.nc'

!According to getscr and irdscr, build _SCR file name, referred as filscr
!A default is available if getscr is 0
 stringfile='_SCR' ; stringvar='scr'
 call mkfilename(filnam,filscr,dtset%getscr,idtset,dtset%irdscr,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filscr=trim(filnam_ds(3))//'_SCR'

!According to getsuscep and irdsuscep, build _SUS file name, referred as filsus
!A default is available if getsuscep is 0
 stringfile='_SUS' ; stringvar='sus'
 call mkfilename(filnam,filsus,dtset%getsuscep,idtset,dtset%irdsuscep,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filsus=TRIM(filnam_ds(3))//'_SUS'

!According to getqps and irdqps, build _QPS file name, referred as filqps
!A default is available if getqps is 0
 stringfile='_QPS' ; stringvar='qps'
 call mkfilename(filnam,filqps,dtset%getqps,idtset,dtset%irdqps,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filqps=trim(filnam_ds(3))//'_QPS'

!According to getbseig and irdbseig, build _BSEIG file name, referred as filbseig
!A default is available if getbseig is 0
 stringfile='_BSEIG' ; stringvar='bseig'
 call mkfilename(filnam,filbseig,dtset%getbseig,idtset,dtset%irdbseig,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filbseig=trim(filnam_ds(3))//'_BSEIG'

!According to gethaydock and irdhaydock, build _HAYD file name, referred as filhaydock.
!A default is available if gethaydock is 0
 stringfile='_HAYDR_SAVE' ; stringvar='haydock'
 call mkfilename(filnam,filhaydock,dtset%gethaydock,idtset,dtset%irdhaydock,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filhaydock=trim(filnam_ds(3))//'_HAYDR_SAVE'

!According to getbsr and irdbsr, build _BSR file name, referred as fil_bsreso
!A default is available if getbsr is 0
 stringfile='_BSR' ; stringvar='bsreso'
 call mkfilename(filnam,fil_bsreso,dtset%getbsreso,idtset,dtset%irdbsreso,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0) fil_bsreso=trim(filnam_ds(3))//'_BSR'

!According to getbsc and irdbsc, build _BSC file name, referred as fil_bscoup
!A default is available if getbsc is 0
 stringfile='_BSC' ; stringvar='bscoup'
 call mkfilename(filnam,fil_bscoup,dtset%getbscoup,idtset,dtset%irdbscoup,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fil_bscoup=trim(filnam_ds(3))//'_BSC'

!According to getwfkfine and irdwfkfine, build _WFK file name, referred as filwfkfine
!A default is avaible if getwfkfine is 0
 stringfile='_WFK' ; stringvar='wfkfine'
 call mkfilename(filnam,filwfkfine,dtset%getwfkfine,idtset,dtset%irdwfkfine,jdtset_,&
& ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filwfkfine=trim(filnam_ds(3))//'_WFK'

 dtfil%ireadden      =ireadden
 dtfil%ireadkden     =ireadkden
 dtfil%ireadwf       =ireadwf
 dtfil%filnam_ds(1:5)=filnam_ds(1:5)

 dtfil%fnameabi_bsham_reso=fil_bsreso
 dtfil%fnameabi_bsham_coup=fil_bscoup
 dtfil%fnameabi_bseig=filbseig
 dtfil%fnameabi_haydock=filhaydock
 dtfil%fnameabi_sus  =filsus
 dtfil%fnameabi_qps  =filqps
 dtfil%fnameabi_scr  =filscr
 dtfil%fnameabi_efmas=fil_efmas
 dtfil%filddbsin     =filddbsin
 dtfil%fildensin     =fildensin
 dtfil%fildens1in    =fildens1in
 dtfil%filkdensin    =filkdensin
 dtfil%filpawdensin  =filpawdensin
 dtfil%fnameabi_wfkfine = filwfkfine
 dtfil%filstat       =filstat
 dtfil%fnamewffk     =fnamewffk
 dtfil%fnamewffq     =fnamewffq
 dtfil%fnamewffddk   =fnamewffddk
 dtfil%fnamewffdelfd =fnamewffdelfd
 dtfil%fnamewffdkdk  =fnamewffdkdk
 dtfil%fnamewffdkde  =fnamewffdkde
 dtfil%fnamewff1     =fnamewff1

 dtfil%fnameabi_hes=trim(dtfil%filnam_ds(3))//'_HES'
 dtfil%fnameabi_phfrq=trim(dtfil%filnam_ds(3))//'_PHFRQ'
 dtfil%fnameabi_phvec=trim(dtfil%filnam_ds(3))//'_PHVEC'

!-------------------------------------------------------------------------------------------
!Build name of files from dtfil%filnam_ds(4)

 dtfil%fnameabo_ddb=trim(dtfil%filnam_ds(4))//'_DDB'
 dtfil%fnameabo_den=trim(dtfil%filnam_ds(4))//'_DEN'
 dtfil%fnameabo_dos=trim(dtfil%filnam_ds(4))//'_DOS'
 dtfil%fnameabo_eelf=trim(dtfil%filnam_ds(4))//'_EELF'
 dtfil%fnameabo_eig=trim(dtfil%filnam_ds(4))//'_EIG'
 dtfil%fnameabo_eigi2d=trim(dtfil%filnam_ds(4))//'_EIGI2D'
 dtfil%fnameabo_eigr2d=trim(dtfil%filnam_ds(4))//'_EIGR2D'
 dtfil%fnameabo_em1=trim(dtfil%filnam_ds(4))//'_EM1'
 dtfil%fnameabo_em1_lf=trim(dtfil%filnam_ds(4))//'_EM1_LF'
 dtfil%fnameabo_em1_nlf=trim(dtfil%filnam_ds(4))//'_EM1_NLF'
 dtfil%fnameabo_fan=trim(dtfil%filnam_ds(4))//'_FAN'
 dtfil%fnameabo_gkk=trim(dtfil%filnam_ds(4))//'_GKK'
 dtfil%fnameabo_gw=trim(dtfil%filnam_ds(4))//'_GW' ! TODO change name
 dtfil%fnameabo_gwdiag=trim(dtfil%filnam_ds(4))//'_GWDIAG'
 dtfil%fnameabo_kss=trim(dtfil%filnam_ds(4))//'_KSS'
 dtfil%fnameabo_moldyn=trim(dtfil%filnam_ds(4))//'_MOLDYN'
 dtfil%fnameabo_pot=trim(dtfil%filnam_ds(4))//'_POT'
 dtfil%fnameabo_qps=trim(dtfil%filnam_ds(4))//'_QPS'
 dtfil%fnameabo_qp_den=trim(dtfil%filnam_ds(4))//'_QP_DEN'
 dtfil%fnameabo_qp_pawden=trim(dtfil%filnam_ds(4))//'_QP_PAWDEN'
 dtfil%fnameabo_qp_dos=trim(dtfil%filnam_ds(4))//'_QP_DOS'
 dtfil%fnameabo_qp_eig=trim(dtfil%filnam_ds(4))//'_QP_DB.nc' ! TODO change name
 dtfil%fnameabo_rpa=trim(dtfil%filnam_ds(4))//'_RPA'
 dtfil%fnameabo_scr=trim(dtfil%filnam_ds(4))//'_SCR'
 dtfil%fnameabo_sgm=trim(dtfil%filnam_ds(4))//'_SGM'
 dtfil%fnameabo_sgr=trim(dtfil%filnam_ds(4))//'_SGR'
 dtfil%fnameabo_sig=trim(dtfil%filnam_ds(4))//'_SIG'
 dtfil%fnameabo_spcur=trim(dtfil%filnam_ds(4))//'_SPCUR'
 dtfil%fnameabo_sus=trim(dtfil%filnam_ds(4))//'_SUS'
 dtfil%fnameabo_vha=trim(dtfil%filnam_ds(4))//'_VHA'
 dtfil%fnameabo_vpsp=trim(dtfil%filnam_ds(4))//'_VPSP'
 dtfil%fnameabo_vso=trim(dtfil%filnam_ds(4))//'_VSO'
 dtfil%fnameabo_vxc=trim(dtfil%filnam_ds(4))//'_VXC'
 dtfil%fnameabo_wan=trim(dtfil%filnam_ds(4))//'_WAN'
 dtfil%fnameabo_wfk=trim(dtfil%filnam_ds(4))//'_WFK'
 dtfil%fnameabo_wfq=trim(dtfil%filnam_ds(4))//'_WFQ'
 dtfil%fnameabo_w90=trim(dtfil%filnam_ds(4))//'_w90'
 dtfil%fnameabo_1wf=trim(dtfil%filnam_ds(4))//'_1WF'
 dtfil%fnameabo_nlcc_derivs=trim(dtfil%filnam_ds(4))//'_nlcc_derivs_'
 dtfil%fnameabo_pspdata=trim(dtfil%filnam_ds(4))//'_pspdata_'

!-------------------------------------------------------------------------------------------
!Build name of files from dtfil%filnam_ds(5)

 dtfil%fnametmp_eig=trim(dtfil%filnam_ds(5))//'_EIG'
 dtfil%fnametmp_1wf1_eig=trim(dtfil%filnam_ds(5))//'_1WF1_EIG' ! This appendix should be changed !
 dtfil%fnametmp_kgs=trim(dtfil%filnam_ds(5))//'_KGS'
 dtfil%fnametmp_sustr=trim(dtfil%filnam_ds(5))//'_SUSTR'
 dtfil%fnametmp_tdexcit=trim(dtfil%filnam_ds(5))//'_TDEXCIT'
 dtfil%fnametmp_tdwf=trim(dtfil%filnam_ds(5))//'_TDWF'
 dtfil%fnametmp_cg=trim(dtfil%filnam_ds(5))//'_cg'
 dtfil%fnametmp_cprj=trim(dtfil%filnam_ds(5))//'_cprj'

!'_WF1' -> dtfil%unwft1
!'_WF2' -> dtfil%unwft2
!'_KG' ->  dtfil%unkg
!'_DUM' -> tmp_unit (real dummy name)
!'_YLM' -> dtfil%unylm
!'_PAW' -> dtfil%unpaw

 tmpfil(1)=trim(dtfil%filnam_ds(5))//'_WF1'  ! tmpfil(1)
 tmpfil(2)=trim(dtfil%filnam_ds(5))//'_WF2'  ! tmpfil(2)

 tmpfil(3)=trim(dtfil%filnam_ds(5))//'_KG'   ! tmpfil(3)
 tmpfil(4)=trim(dtfil%filnam_ds(5))//'_DUM'  ! tmpfil(4)
 tmpfil(5)=' '  ! to avoid Valgrind complain.
 tmpfil(6)=trim(dtfil%filnam_ds(5))//'_YLM'  ! tmpfil(6)
 tmpfil(7)=trim(dtfil%filnam_ds(5))//'_PAW'  ! tmpfil(7)

 if(xmpi_paral==1)then ! parallel case : the index of the processor must be appended
   call int2char4(mpi_enreg%me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   ixx=1
   if (xmpi_mpiio == 1 .and. dtset%iomode == IO_MODE_MPI ) ixx=3
   do ii=ixx,7
     tmpfil(ii)=trim(tmpfil(ii))//'_P-'//trim(tag)
   end do
 end if

 dtfil%fnametmp_wf1=trim(tmpfil(1))
 dtfil%fnametmp_wf2=trim(tmpfil(2))

 dtfil%fnametmp_kg =trim(tmpfil(3))
 dtfil%fnametmp_dum=trim(tmpfil(4))
 dtfil%fnametmp_ylm=trim(tmpfil(6))
 dtfil%fnametmp_paw=trim(tmpfil(7))

!Create names for the temporary files based on dtfil%filnam_ds(5)
!by appending adequate string.
!'_1WF1' -> dtfil%unwft1
!'_1WF2' -> dtfil%unwft2
!'_KG'   -> dtfil%unkg
!'_KGQ'  -> dtfil%unkgq (not used for the time being)
!'_KG1'  -> dtfil%unkg1
!'_DUM'  -> tmp_unit (real dummy name)
!'_WFGS' -> dtfil%unwftgs
!'_WFKQ' -> dtfil%unwftkq
!'_YLM'  -> dtfil%unylm
!'_YLM1' -> dtfil%unylm1
!'_PAW'  -> dtfil%unpaw
!'_PAW1' -> dtfil%unpaw1
!'_PAWQ' -> dtfil%unpawq
 tmpfil(1) =trim(dtfil%filnam_ds(5))//'_1WF1'
 tmpfil(2) =trim(dtfil%filnam_ds(5))//'_1WF2'
 tmpfil(3) =trim(dtfil%filnam_ds(5))//'_KG'
 tmpfil(4) =trim(dtfil%filnam_ds(5))//'_KGQ'
 tmpfil(5) =trim(dtfil%filnam_ds(5))//'_KG1'
 tmpfil(6) =trim(dtfil%filnam_ds(5))//'_DUM'
 tmpfil(7) =trim(dtfil%filnam_ds(5))//'_WFGS'
 tmpfil(8) =trim(dtfil%filnam_ds(5))//'_WFKQ'
 tmpfil(9) =' ' ! for Valgrind, to avoid uninitialized
 tmpfil(10)=trim(dtfil%filnam_ds(5))//'_YLM'
 tmpfil(11)=trim(dtfil%filnam_ds(5))//'_YLM1'
 tmpfil(12)=trim(dtfil%filnam_ds(5))//'_PAW'
 tmpfil(13)=trim(dtfil%filnam_ds(5))//'_PAW1'
 tmpfil(14)=trim(dtfil%filnam_ds(5))//'_PAWQ'

 if(xmpi_paral==1) then
   idtmpfil(:)=0
   do ii=1,14
     idtmpfil(ii)=ii
   end do
   if (xmpi_mpiio==1.and.dtset%iomode==IO_MODE_MPI)then
     idtmpfil(1)=0              !_1wf1
     idtmpfil(2)=0              ! s1wf2
     idtmpfil(7)=0              !  WFGS
     idtmpfil(8)=0              !  WFKQ
   end if
   call int2char4(mpi_enreg%me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   do ii=1,14
     if(idtmpfil(ii) /= 0) tmpfil(ii)=trim(tmpfil(ii))//'_P-'//trim(tag)
   end do
 end if

 dtfil%fnametmp_1wf1=trim(tmpfil(1))
 dtfil%fnametmp_1wf2=trim(tmpfil(2))
 dtfil%fnametmp_kg  =trim(tmpfil(3))
 dtfil%fnametmp_kgq =trim(tmpfil(4))
 dtfil%fnametmp_kg1 =trim(tmpfil(5))
 dtfil%fnametmp_dum =trim(tmpfil(6))
 dtfil%fnametmp_wfgs=trim(tmpfil(7))
 dtfil%fnametmp_wfkq=trim(tmpfil(8))
 dtfil%fnametmp_ylm =trim(tmpfil(10))
 dtfil%fnametmp_ylm1=trim(tmpfil(11))
 dtfil%fnametmp_paw =trim(tmpfil(12))
 dtfil%fnametmp_paw1=trim(tmpfil(13))
 dtfil%fnametmp_pawq=trim(tmpfil(14))

!Prepare the name of the _FFT file
 filfft=trim(dtfil%filnam_ds(5))//'_FFT'
!There is a definite problem in the treatment of // by CPP ...
 if(xmpi_paral==1 .or. mpi_enreg%paral_kgb==1)then
   call int2char4(mpi_enreg%me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   filfft=trim(filfft)//'_P-'//trim(tag)
 end if
 dtfil%fnametmp_fft=filfft

!These keywords are only used in algorithms using images of the cell
 if (iimage==0) then
   dtfil%getwfk_from_image   =0
   dtfil%getden_from_image   =0
   dtfil%getpawden_from_image=0
 end if

 DBG_EXIT("COLL")

end subroutine dtfil_init
!!***

!!****f* m_dtfil/dtfil_init_time
!!
!! NAME
!! dtfil_init_time
!!
!! FUNCTION
!! Inside the itimimage, iimage and itime loops (this is only needed for optdriver=0),
!! initialize the remaining parts of dtfil.
!!
!! INPUTS
!! iapp=indicates the eventual suffix to be appended to the generic output root
!!         if 0 : no suffix to be appended (called directly from gstate)
!!         if positive : append "_TIM//iapp" (called from move or brdmin)
!!         if -1 : append "_TIM0" (called from brdmin)
!!         if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtfil=<type datafiles_type>infos about file names, file unit numbers
!!  (part of which were initialized previously)
!!
!! PARENTS
!!      gstate,mover
!!
!! CHILDREN
!!      fappnd
!!
!! SOURCE

subroutine dtfil_init_time(dtfil,iapp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init_time'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: iapp
 type(datafiles_type),intent(inout) :: dtfil

!Local variables-------------------------------
!scalars
 character(len=fnlen) :: filapp,filprot

!******************************************************************

 DBG_ENTER("COLL")

!--------------------------------------------------------
!Names based on dtfil%filnam_ds(4)+iapp

!Prepare the name of the auxiliary files DOS, EIG...
 call fappnd(filapp,dtfil%filnam_ds(4),iapp)
 dtfil%fnameabo_app=trim(filapp)
 dtfil%fnameabo_app_atmden_core=trim(filapp)//'_ATMDEN_CORE'
 dtfil%fnameabo_app_atmden_val=trim(filapp)//'_ATMDEN_VAL'
 dtfil%fnameabo_app_atmden_full=trim(filapp)//'_ATMDEN_FULL'
 dtfil%fnameabo_app_n_tilde=trim(filapp)//'_N_TILDE'
 dtfil%fnameabo_app_n_one=trim(filapp)//'_N_ONE'
 dtfil%fnameabo_app_nt_one=trim(filapp)//'_NT_ONE'
 dtfil%fnameabo_app_bxsf=trim(filapp)//'_BXSF'
 dtfil%fnameabo_app_cif=trim(filapp)//'.cif'
 dtfil%fnameabo_app_den=trim(filapp)//'_DEN'
 dtfil%fnameabo_app_dos=trim(filapp)//'_DOS'
 dtfil%fnameabo_app_eig=trim(filapp)//'_EIG'
 dtfil%fnameabo_app_elf=trim(filapp)//'_ELF'
 dtfil%fnameabo_app_elf_down=trim(filapp)//'_ELF_DOWN'
 dtfil%fnameabo_app_elf_up=trim(filapp)//'_ELF_UP'
 dtfil%fnameabo_app_fatbands=trim(filapp)//'_FATBANDS'
 dtfil%fnameabo_app_gden1=trim(filapp)//'_GDEN1'
 dtfil%fnameabo_app_gden2=trim(filapp)//'_GDEN2'
 dtfil%fnameabo_app_gden3=trim(filapp)//'_GDEN3'
 dtfil%fnameabo_app_geo=trim(filapp)//'_GEO'
 dtfil%fnameabo_app_kden=trim(filapp)//'_KDEN'
 dtfil%fnameabo_app_lden=trim(filapp)//'_LDEN'
 dtfil%fnameabo_app_nesting=trim(filapp)//'_NEST'
 dtfil%fnameabo_app_opt=trim(filapp)//'_OPT'
 dtfil%fnameabo_app_opt2=trim(filapp)//'_OPT2'
 dtfil%fnameabo_app_pawden=trim(filapp)//'_PAWDEN'
 dtfil%fnameabo_app_pot=trim(filapp)//'_POT'
 dtfil%fnameabo_app_stm=trim(filapp)//'_STM'
 dtfil%fnameabo_app_vclmb=trim(filapp)//'_VCLMB'
 dtfil%fnameabo_app_vha=trim(filapp)//'_VHA'
 dtfil%fnameabo_app_vhxc=trim(filapp)//'_VHXC'
 dtfil%fnameabo_app_vpsp=trim(filapp)//'_VPSP'
 dtfil%fnameabo_app_vxc=trim(filapp)//'_VXC'
 dtfil%fnameabo_app_wfk=trim(filapp)//'_WFK'
 dtfil%fnameabo_app_vha_1dm=trim(filapp)//'_VHA_1DM'
 dtfil%fnameabo_app_vclmb_1dm=trim(filapp)//'_VCLMB_1DM'
 dtfil%fnameabo_app_1dm=trim(filapp)//'_1DM'

!--------------------------------------------------------
!Names based on dtfil%filnam_ds(5)+iapp

!Prepare the name of the auxiliary files for protection
 call fappnd(filprot,dtfil%filnam_ds(5),iapp)
 dtfil%fnametmp_app_den=trim(filprot)//'_DEN'
 dtfil%fnametmp_app_kden=trim(filprot)//'_KDEN'

 DBG_EXIT("COLL")

end subroutine dtfil_init_time
!!***

!!****f* m_dtfil/fappnd
!!
!! NAME
!! fappnd
!!
!! FUNCTION
!! Create the modified root name to be used for output of density, potential,
!! and geometry files. See the description of the iapp input variable.
!!
!! INPUTS
!! filnam= generic output root name
!! iapp=indicates the eventual suffix to be appended to the generic output root
!!      (the suffixe depends on the presence of the suff (optional) argument.
!!        if 0 : no suffix to be appended (called directly from gstate)
!!        if positive : append "_SUF//iapp" (called from move or brdmin)
!!        if -1 : append "_SUF0" (called from brdmin)
!!        if -2, -3, -4, -5: append "_SUFA", ... ,"_SUFD", (called from move)
!!      SUF can be TIM (default) or IMG
!! suff= --optional argument--indicates the suffixe to be appended:
!!       SUF=TIM (default) or SUF=IMG or ...
!!
!! OUTPUT
!! filapp= filename with appended string
!!
!! PARENTS
!!      dtfil_init_time
!!
!! CHILDREN
!!
!! SOURCE

subroutine fappnd(filapp,filnam,iapp,&
&                 suff) ! optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fappnd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp
 character(len=fnlen),intent(in) :: filnam
 character(len=fnlen),intent(out) :: filapp
 character(len=3),optional,intent(in) :: suff

!Local variables-------------------------------
!scalars
 integer :: ndig
 character(len=3) :: suffixe
 character(len=8) :: nchar
 character(len=500) :: msg

! *************************************************************************

 if(iapp==0)then
   filapp=trim(filnam)
 else
   suffixe="TIM"
   if (present(suff)) suffixe=trim(suff(1:3))
   if(iapp>0)then
!    Create character string for filename. Determine the number of digits in iapp.
     ndig=int(log10(dble(iapp)+0.5_dp))+1
!    Make integer format field of exact size (internal write)
!    for assumed nchar string of 8 characters
     write(nchar, '(i8)' ) iapp
     if (ndig>8) then
       write(msg,'(5a,i0,2a,i0,2a)')&
&       'Requested file name extension has more than the allowed 8 digits.',ch10,&
&       'Action: resubmit the job with smaller value for ntime.',ch10,&
&       'Value computed here was ndig=',ndig,ch10,&
&       'iapp= ',iapp,' filnam=',TRIM(filnam)
       MSG_ERROR(msg)
     end if
!    Concatenate into character string, picking off exact number of digits
!    The potential or density label will be appended in ioarr
     filapp=trim(filnam)//'_'//suffixe(1:3)//nchar(9-ndig:8)
   else if(iapp==-1)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'0'
   else if(iapp==-2)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'A'
   else if(iapp==-3)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'B'
   else if(iapp==-4)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'C'
   else if(iapp==-5)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'D'
   end if
 end if

end subroutine fappnd
!!***

!!****f* m_dtfil/dtfil_init_img
!! NAME
!! dtfil_init_img
!!
!! FUNCTION
!! Initialize few scalars in the dtfil structured variable
!! when an alogrithm using image of the cell is selected.
!! (initialize index of images from which read files)
!!
!! INPUTS
!!  dtset=<type datasets_type>=input variables for the current dataset
!!  dtsets(0:ndtset_alloc)=<type datasets_type>=input variables for all datasets
!!  idtset=number of the dataset
!!  jdtset(0:ndtset)=actual index of the datasets
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtfil=<type datafiles_type>= only getxxx_from_image flags are modified
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset,ndtset,ndtset_alloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: idtset,ndtset,ndtset_alloc
 type(datafiles_type),intent(out) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 integer :: jdtset(0:ndtset)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------
!scalars
 integer :: iget
!arrays

! *********************************************************************

 DBG_ENTER("COLL")

!Default values
 dtfil%getwfk_from_image   =0 ! Get standard WFK from previous dataset
 dtfil%getden_from_image   =0 ! Get standard DEN from previous dataset
 dtfil%getpawden_from_image=0 ! Get standard PAWDEN from previous dataset

 if (dtset%optdriver==RUNL_GSTATE.and.dtset%nimage>1) then

!  Define getwfk_from_image
   if (dtset%getwfk/=0.or.dtset%irdwfk/=0) then
     iget=-1
     if(dtset%getwfk<0) iget=jdtset(idtset+dtset%getwfk)
     if(dtset%getwfk>0) iget=dtset%getwfk
     if(dtset%irdwfk>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getwfk_from_image=-1     ! Get WFK from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getwfk_from_image=1      ! Get WFK from the first image of previous dataset
       end if
     end if
   end if

!  Define getden_from_image
   if (dtset%getden/=0.or.dtset%irdden/=0) then
     iget=-1
     if(dtset%getden<0) iget=jdtset(idtset+dtset%getden)
     if(dtset%getden>0) iget=dtset%getden
     if(dtset%irdden>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getden_from_image=-1     ! Get DEN from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getden_from_image=1      ! Get DEN from the first image of previous dataset
       end if
     end if
   end if

!  Define getpawden_from_image
   if (dtset%getpawden/=0.or.dtset%irdpawden/=0) then
     iget=-1
     if(dtset%getpawden<0) iget=jdtset(idtset+dtset%getpawden)
     if(dtset%getpawden>0) iget=dtset%getpawden
     if(dtset%irdpawden>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getpawden_from_image=-1     ! Get PAWDEN from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getpawden_from_image=1      ! Get PAWDEN from the first image of previous dataset
       end if
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine dtfil_init_img
!!***

!!****f* m_dtfil/mkfilename
!!
!! NAME
!! mkfilename
!!
!! FUNCTION
!! From the root (input or output) file names, produce a real file name.
!!
!! INPUTS
!! character(len=fnlen):: filnam(5)=the root file names
!!  (only filnam(3) and filnam(4) are really needed)
!! get=input 'get variable', if 1, must get the file from another dataset
!! idtset=number of the dataset
!! ird=input 'iread variable', if 1, must get the file from the input root
!! jdtset_(0:ndtset)=actual index of the dataset
!! ndtset=number of datasets
!! stringfil character(len=*)=the string of characters to be appended e.g. '_WFK' or '_DEN'
!! stringvar tcharacter(len=*)=the string of characters to be appended
!!   that defines the 'get' or 'ird' variables, e.g. 'wfk' or 'ddk'
!!
!! OUTPUT
!! character(len=fnlen):: filnam_out=the new file name
!! will_read=1 if the file must be read ; 0 otherwise (ird and get were zero)
!!
!! PARENTS
!!      dtfil_init,finddistrproc
!!
!! CHILDREN
!!      appdig,wrtout
!!
!! SOURCE

subroutine mkfilename(filnam,filnam_out,get,idtset,ird,jdtset_,ndtset,stringfil,stringvar,will_read)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkfilename'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: get,idtset,ird,ndtset
 integer,intent(out) :: will_read
 character(len=*),intent(in) :: stringfil
 character(len=*),intent(in) :: stringvar
 character(len=fnlen),intent(out) :: filnam_out
!arrays
 integer,intent(in) :: jdtset_(0:ndtset)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
!scalars
 integer :: jdtset,jget
 character(len=4) :: appen
 character(len=500) :: message
 character(len=fnlen) :: filnam_appen

! *************************************************************************

!Here, defaults if no get variable
 will_read=ird

 filnam_appen=trim(filnam(3))
 if(ndtset>0)then
   jdtset=jdtset_(idtset)
   call appdig(jdtset,'',appen)
   filnam_appen=trim(filnam_appen)//'_DS'//appen
 end if
 filnam_out=trim(filnam_appen)//trim(stringfil)

!Treatment of the multi-dataset case  (get is not relevant otherwise)
 if(ndtset/=0)then

   if(ndtset==1.and.get<0.and.(jdtset_(1)+get>0))then
     write(message, '(7a,i3,a,i3,5a)' )&
&     'You cannot use a negative value of get',trim(stringvar),' with only 1 dataset!',ch10, &
&     ' If you want to refer to a previously computed dataset,',ch10, &
&     ' you should give the absolute index of it (i.e. ', &
&     jdtset_(idtset)+get,' instead of ',get,').',ch10, &
&     'Action: correct get',trim(stringvar),' in your input file.'
     MSG_ERROR(message)
   end if

   if(idtset+get<0)then
     write(message, '(a,a,a,a,a,i3,a,a,a,i3,a,a,a,a)' )&
&     'The sum of idtset and get',trim(stringvar),' cannot be negative,',ch10,&
&     'while they are idtset=',idtset,', and get',trim(stringvar),'=',get,ch10,&
&     'Action: correct get',trim(stringvar),' in your input file.'
     MSG_ERROR(message)
   end if

   if(get>0 .or. (get<0 .and. idtset+get>0) )then

     if(ird/=0 .and. get/=0)then
       write(message, '(a,a,a,a,a,a,a,a,a,a,a,i3,a,i3,a,a,a,a,a,a,a)' )&
&       'The input variables ird',trim(stringvar),' and get',trim(stringvar),' cannot be',ch10,&
&       'simultaneously non-zero, while for idtset=',idtset,',',ch10,&
&       'they are ',ird,', and ',get,'.',ch10,&
&       'Action: correct ird',trim(stringvar),' or get',trim(stringvar),' in your input file.'
       MSG_ERROR(message)
     end if

     will_read=1

!    Compute the dataset from which to take the file, and the corresponding index
     if(get<0 .and. idtset+get>0) jget=jdtset_(idtset+get)
     if(get>0) jget=get
     call appdig(jget,'',appen)

!    Note use of output filename (filnam(4))
     filnam_out=trim(filnam(4))//'_DS'//trim(appen)//trim(stringfil)

     if(jdtset>=100)then
       write(message, '(a,a,a,a,a,i5,a,a)' )&
&       ' mkfilename : get',trim(stringvar) ,'/=0, take file ',trim(stringfil),&
&       ' from output of DATASET ',jget,'.',ch10
     else
       write(message, '(a,a,a,a,a,i3,a,a)' )&
&       ' mkfilename : get',trim(stringvar) ,'/=0, take file ',trim(stringfil),&
&       ' from output of DATASET ',jget,'.',ch10
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if ! conditions on get and idtset

 end if ! ndtset/=0

end subroutine mkfilename
!!***


!!****f* m_dtfil/isfile
!! NAME
!! isfile
!!
!! FUNCTION
!! Inquire Status of FILE
!! Checks that for status =
!! 'old': file already exists
!! 'new': file does not exist; if file exists,
!! filnam is modified to filnam.A or filnam.B,....
!!
!! INPUTS
!! filnam=character string to specify filename
!! status='old' or 'new'
!!
!! OUTPUT
!! stops processing if old file does not exist; changes name
!! and returns new name in redefined filnam if new file already exists.
!!
!! PARENTS
!!      anaddb,iofn1,m_effective_potential,m_polynomial_coeff,m_vcoul
!!      multibinit,ujdet
!!
!! CHILDREN
!!      clib_rename,int2char4
!!
!! SOURCE

subroutine isfile(filnam, status)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isfile'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=3),intent(in) :: status
 character(len=fnlen),intent(inout) :: filnam

!Local variables-------------------------------
!scalars
 logical :: ex,found
 integer :: ii,ios, ioserr
 character(len=500) :: message
 character(len=fnlen) :: filnam_tmp
 character(len=fnlen) :: trialnam

! *************************************************************************

 filnam_tmp=filnam

 if (status=='old') then !  Check that old file exists
   inquire(file=filnam,iostat=ios,exist=ex)

   if (ios/=0) then
     write(message,'(4a,i0,2a)')&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE statement returns error code',ios,ch10,&
&     'Action: identify which problem appears with this file.'
     MSG_ERROR(message)
   else if (.not.ex) then
     write(message, '(5a)' )&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE finds file does not exist.',&
&     'Action: check file name and re-run.'
     MSG_ERROR(message)
   end if

 else if (status=='new') then

   ! Check that new output file does NOT exist
   ioserr = 0
   trialnam = filnam
   ii = 0
   inquire(file=trim(trialnam),iostat=ios,exist=ex)
   if ( ios /= 0 ) then
     write(message,'(4a)') 'Something is wrong with permissions for ', &
&     'reading/writing on this filesystem.',ch10,&
&     'Action : Check permissions.'
     MSG_ERROR(message)
   end if

   if ( ex .eqv. .true. ) then
     write(message,'(4a)')'Output file ',trim(trialnam),ch10,' already exists.'
     MSG_COMMENT(message)
     found=.false.

     ii=1
     do while ( (found .eqv. .false.) .and. (ii < 10000) )
       call int2char4(ii,message)
       trialnam=trim(trim(filnam_tmp)//message)
       inquire(file=trim(trialnam),iostat=ios,exist=ex)
       if ( (ex .eqv. .false.) .and. (ios == 0)) then
         found  = .true.
       end if
       if ( ios /= 0 )  ioserr=ioserr+1
       if ( ioserr > 10 ) then
!        There is a problem => stop
         write(message, '(2a,i0,2a)' )&
&         'Check for permissions of reading/writing files on the filesystem', &
&         '10 INQUIRE statements returned an error code like ',ios,ch10,&
&         'Action: Check permissions'
         MSG_ERROR(message)
       end if
       ii=ii+1
     end do
     if ( found .eqv. .true. ) then
       write(message,'(4a)') 'Renaming old ',trim(filnam),' to ',trim(trialnam)
       MSG_COMMENT(message)
       call clib_rename(filnam,trialnam,ioserr)
       if ( ioserr /= 0 ) then
         write(message,'(4a)') 'Failed to rename file ', trim(filnam),' to ',trim(trialnam)
         MSG_ERROR(message)
       end if
     else
       write(message,'(3a)')&
&       'Have used all names of the form filenameXXXX, X in [0-9]',ch10,&
&       'Action: clean up your directory and start over.'
       MSG_ERROR(message)
     end if
   end if

   ! if ii > 0 we iterated so rename abi_out to abi_outXXXX
   ! and just write to abi_out
 else ! status not recognized
   write(message,'(3a)')'  Input status= ',status,' not recognized.'
   MSG_BUG(message)
 end if

end subroutine isfile
!!***

!!****f* m_dtfil/iofn1
!! NAME
!! iofn1
!!
!! FUNCTION
!! Begin by eventual redefinition of unit std_in and std_out
!! Then, print greetings for interactive user.
!! Next, Read filenames from unit std_in, AND check that new
!! output file does not already exist.
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  character(len=fnlen) :: filnam(5)=character strings giving file names
!!  character(len=fnlen) :: filstat=character strings giving name of status file
!!
!! NOTES
!! If it does exist, isfile will create a new name to avoid overwriting the output file.
!! Also create name of status file
!!
!! File names refer to following files, in order:
!!  (1) Formatted input file  (std_in)
!!  (2) Formatted output file (std_out)
!!  (3) Root name for generic input files (wavefunctions, potential, density ...)
!!  (4) Root name for generic output files (wavefunctions, potential, density,
!!                                          DOS, hessian ...)
!!  (5) Root name for generic temporary files (wftmp1,wftmp2,kgunit,status ...)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,isfile,libpaw_log_flag_set,xmpi_barrier
!!      xmpi_bcast
!!
!! SOURCE

subroutine iofn1(filnam,filstat,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'iofn1'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 character(len=fnlen), intent(out) :: filstat
 character(len=fnlen), intent(out) :: filnam(5)

!Local variables-------------------------------
 character(len=1) :: blank
 integer,parameter :: master=0
 integer :: me,ios,nproc,ierr
 logical :: ex
 character(len=fnlen) :: fillog,tmpfil
 character(len=10) :: tag
 character(len=500) :: message,errmsg

!*************************************************************************

 ! NOTE: In this routine it's very important to perform tests
 ! on possible IO errors (err=10, iomsg) because we are initializing the IO stuff
 ! It there's some problem with the hardware or some misconfiguration,
 ! it's very likely that the code will crash here and we should try to give useful error messages.

 blank = ' '; tmpfil = ''

!Determine who I am in comm
 me = xmpi_comm_rank(comm)
 nproc = xmpi_comm_size(comm)

!Define values of do_write_log and do_write_status parameters
!if a _NOLOG file exists no LOG file and no STATUS file are created for each cpu core
!if a _LOG file exists, a LOG file and a STATUS file are created for each cpu core
!if the #_of_cpu_core>NPROC_NO_EXTRA_LOG OR presence of ABI_MAIN_LOG_FILE, LOG file is only created for master proc
!if the #_of_cpu_core>NPROC_NO_EXTRA_STATUS OR presence of ABI_MAIN_LOG_FILE, STATUS file is only created for master proc
 inquire(file=ABI_NO_LOG_FILE,iostat=ios,exist=ex)
 if (ios/=0) ex=.false.
 if (ex) then
   do_write_log=.false. ; do_write_status=.false.
   call abi_log_status_state(new_do_write_log=.false.,new_do_write_status=.false.)
   call libpaw_log_flag_set(.false.)
 else
   inquire(file=ABI_ENFORCE_LOG_FILE,iostat=ios,exist=ex)
   if (ios/=0) ex=.false.
   if (ex) then
     do_write_log=.true. ; do_write_status=.true.
     call abi_log_status_state(new_do_write_log=.true.,new_do_write_status=.true.)
     call libpaw_log_flag_set(.true.)
   else
     inquire(file=ABI_MAIN_LOG_FILE,iostat=ios,exist=ex)
     if (ios/=0) ex=.false.
     if (ex.and.me/=0) then
       do_write_log=.false. ; do_write_status=.false.
       call abi_log_status_state(new_do_write_log=.false.,new_do_write_status=.false.)
       call libpaw_log_flag_set(.false.)
     else
       if (me/=0) then
         do_write_log= (nproc<NPROC_NO_EXTRA_LOG)
         call abi_log_status_state(new_do_write_log=(nproc<NPROC_NO_EXTRA_LOG))
         call libpaw_log_flag_set((nproc<NPROC_NO_EXTRA_LOG))
         do_write_status= (nproc<NPROC_NO_EXTRA_STATUS)
         call abi_log_status_state(new_do_write_status=(nproc<NPROC_NO_EXTRA_STATUS))
       end if
     end if
   end if
 end if

 if (me==master) then

!  Eventually redefine standard input and standard output

   if (do_write_log) then
#if defined READ_FROM_FILE
!    Take care of the output file
     tmpfil(1:fnlen)=blank
     tmpfil(1:3)='log'
     call isfile(tmpfil,'new')
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(tmpfil,message,unit=std_out,form='formatted',status='new',action="write") /= 0) then
       MSG_ERROR(message)
     end if
#endif
   else
!    Redirect standard output to null
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(NULL_FILE,message,unit=std_out,action="write") /= 0) then
       MSG_ERROR(message)
     end if
   end if

#if defined READ_FROM_FILE
!  Now take care of the "files" file
   tmpfil(1:fnlen)=blank
   tmpfil(1:9)='ab.files'
   write(message, '(4a)' )&
&   'Because of cpp option READ_FROM_FILE,',ch10,&
&   'read file "ab.files" instead of standard input ' ,ch10
   MSG_COMMENT(message)
   call isfile(tmpfil,'old')
   close(std_in, err=10, iomsg=errmsg)
   if (open_file(tmpfil,message,unit=std_in,form='formatted',status='old',action="read") /= 0) then
     MSG_ERROR(message)
   end if
#endif

!  Print greetings for interactive user
   write(std_out,*,err=10,iomsg=errmsg)' ABINIT ',trim(abinit_version)
   write(std_out,*,err=10,iomsg=errmsg)' '

!  Read name of input file (std_in):
   write(std_out,*,err=10,iomsg=errmsg)' Give name for formatted input file: '
   read(std_in, '(a)',err=10,iomsg=errmsg ) filnam(1)
   write(std_out, '(a)',err=10,iomsg=errmsg ) trim(filnam(1))
   write(std_out,*)' Give name for formatted output file:'
   read (std_in, '(a)',err=10,iomsg=errmsg ) filnam(2)
   write (std_out, '(a)',err=10,iomsg=errmsg ) trim(filnam(2))
   write(std_out,*)' Give root name for generic input files:'
   read (std_in, '(a)',err=10,iomsg=errmsg ) filnam(3)
   write (std_out, '(a)',err=10,iomsg=errmsg ) trim(filnam(3))
   write(std_out,*, err=10, iomsg=errmsg )' Give root name for generic output files:'
   read (std_in, '(a)', err=10, iomsg=errmsg ) filnam(4)
   write (std_out, '(a)', err=10, iomsg=errmsg ) trim(filnam(4))
   write(std_out,*, err=10, iomsg=errmsg)' Give root name for generic temporary files:'
   read (std_in, '(a)', err=10, iomsg=errmsg ) filnam(5)
   write (std_out, '(a)', err=10, iomsg=errmsg ) trim(filnam(5))

!  Check that old input file exists
   call isfile(filnam(1),'old')

!  Check that new output file does NOT exist
   call isfile(filnam(2),'new')

!  Check that root name for generic input and output differ
   if ( trim(filnam(3))==trim(filnam(4)) ) then
     write(message, '(a,a,a)' )&
&     'Root name for generic input and output files must differ ',ch10,&
&     'Action: correct your "file" file.'
     MSG_ERROR(message)
   end if

!  Check that root names are at least 20 characters less than fnlen
   if ( len_trim(filnam(3)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic input files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(3)),ch10,&
&     'Action : correct your "file" file.'
     MSG_ERROR(message)
   end if
   if ( len_trim(filnam(4)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic output files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(4)),ch10,&
&     'Action: correct your "file" file.'
     MSG_ERROR(message)
   end if
   if ( len_trim(filnam(5)) >= (fnlen-20) ) then
     write(message, '(a,a,a,a,a,i4,a,i4,a,a)' )&
&     'Root name for generic temporary files is too long. ',ch10,&
&     'It must be 20 characters less than the maximal allowed ',ch10,&
&     'length of names, that is ',fnlen,', while it is ',len_trim(filnam(5)),ch10,&
&     'Action: correct your "file" file.'
     MSG_ERROR(message)
   end if

 end if ! master only

!Communicate filenames to all processors
 call xmpi_bcast(filnam,master,comm,ierr)

!Check
!Create a name for the status file, based on filnam(5)
 filstat=trim(filnam(5))//'_STATUS'

!Redefine the log unit if not the master
 if (me/=master) then
   call int2char4(me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   filstat=trim(filstat)//'_P-'//trim(tag)
   if (do_write_log) then
     fillog=trim(filnam(5))//'_LOG_'//trim(tag)
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(fillog,message,unit=std_out,status='unknown',action="write") /= 0) then
       MSG_ERROR(message)
     end if
   else
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(NULL_FILE,message,unit=std_out,action="write") /= 0) then
       MSG_ERROR(message)
     end if
   end if
 end if

 call xmpi_barrier(comm)
 return

! Handle possibe IO errors
 10 continue
 MSG_ERROR(errmsg)

end subroutine iofn1
!!***


!!****f* ABINIT/status
!! NAME
!! status
!!
!! FUNCTION
!! Routine for description of the status of the calculation
!! Eventually open the status file, write different information,
!! and close the file. The output rate and shift are governed by istat
!!
!! INPUTS
!!  counter=value of the loop counter at that level
!!  file (optional argument)=name of the status file
!!  istat=gives the rate or shift of output. The status file will be opened
!!     and written only once every "istatr" calls.
!!     This variable is saved at the first call (just after the first
!!     call to invars0. The shift "istatshft" is saved at the second call.
!!     In subsequent calls, istat has no meaning.
!!  level=number of the level of the calling subroutine (see the description later)
!!  routine=string of 14 characters indicating the status inside the level
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Warning : The string "routine" can have any size but
!! it is truncated to a size of 14.
!! because of the behaviour of some compilers, the string
!! "routine" should always have 14 characters in the calling subroutine
!!
!! PARENTS
!!      abinit,dfpt_accrho,dfpt_looppert,dfpt_rhofermi,dfpt_scfcv,dfpt_vtowfk
!!      dfpt_wfkfermi,dfptnl_loop,dfptnl_mv,dfptnl_resp,driver,gstate,gstateimg
!!      m_ab7_invars_f90,mover,nonlinear,respfn,scfcv,vtorhotf
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

subroutine status(counter,filstat,istat,level,routine)

 use m_time,       only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'status'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,istat,level
 character(len=*),intent(in) :: routine
 character(len=*),intent(in) :: filstat

!Local variables-------------------------------
!scalars
 integer,parameter :: mcounter=2,mlevel=120
 integer,save :: output_rate=1,shift_rate=1,statnu=0
 integer :: ilevel,temp_unit
 character(len=12) :: headwr
 character(len=500) :: message
!arrays
 integer,save :: active(mlevel),actual_counter(mlevel,mcounter)
 integer,save :: ncounter(mlevel)
 integer,save :: list_level(29)=&
&  (/1,2,3,100,101,102,110,111,112,113,114,10,11,12,13,14,15,16,17,18,20,30,31,40,41,50,51,52,90/)
 real(dp) :: tsec(2)
 character(len=20),save :: nm_levels(mlevel),nm_routine(mlevel)
 character(len=12),save :: nm_counter(mlevel,mcounter)

!***********************************************************************

 if (.not.do_write_status .or. output_rate==0) return

 call timab(73,1,tsec)

!Note : all processors have their own file, so no special
!attention must be paid to the parallel case.
!Initialisation
 if(statnu==0)then
   nm_routine(:)='                  '
   active(:)=0
   actual_counter(:,:)=0

!  List of names for each level
!  Numbers from 1 to 9 are for abinit and driver
!  Numbers from 100 to 120 are for optdriver=0 routines (GS)
!  Numbers between 10 and 19 are for optdriver=1 routines (RF)
!  Numbers between 20 and 29 are for optdriver=2 routines (suscep)
!  Numbers between 30 and 39 are for optdriver=3 routines (screening)
!  Numbers between 40 and 49 are for optdriver=4 routines (sigma)
!  Numbers between 50 and 59 are for optdriver=5 routines (nonlinear)
!  When you add a level number, or modify one, do not forget to change list_level

   nm_levels(1)   ='abinit              '
   ncounter(1)=0
   nm_counter(1,1)='            '

   nm_levels(2)   ='driver              '
   ncounter(2)=1
   nm_counter(2,1)='jdtset     ='

   nm_levels(3)   ='ab7_invars_load     '
   ncounter(3)=0
   nm_counter(3,1)='            '

!  Optdriver=0
   nm_levels(100)   ='gstateimg           '
   ncounter(100)=2
   nm_counter(100,1)='idynimage  ='
   nm_counter(100,2)='itimimage  ='

   nm_levels(101)   ='gstate              '
   ncounter(101)=1
   nm_counter(101,1)='itime      ='

   nm_levels(102)   ='mover               '
   ncounter(102)=2
   nm_counter(102,1)='icycle     ='
   nm_counter(102,2)='itime      ='

   nm_levels(110)   ='scfcv               '
   ncounter(110)=1
   nm_counter(110,1)='istep      ='

   nm_levels(111)   ='vtorho(tf)          '
   ncounter(111)=2
   nm_counter(111,1)='isppol     ='
   nm_counter(111,2)='ikpt       ='

   nm_levels(112)   ='vtowfk              '
   ncounter(112)=2
   nm_counter(112,1)='inonsc     ='
   nm_counter(112,2)='iband      ='

   nm_levels(113)   ='cgwf                '
   ncounter(113)=1
   nm_counter(113,1)='iline      ='

   nm_levels(114)   ='getghc              '
   ncounter(114)=0
   nm_counter(114,1)='            '


!  Optdriver=1
   nm_levels(10)   ='respfn              '
   ncounter(10)=0
   nm_counter(10,1)='            '

   nm_levels(11)   ='dfpt_looppert       '
   ncounter(11)=1
   nm_counter(11,1)='respcase   ='

   nm_levels(12)   ='dfpt_scfcv          '
   ncounter(12)=1
   nm_counter(12,1)='istep      ='

   nm_levels(13)   ='dfpt_vtorho         '
   ncounter(13)=2
   nm_counter(13,1)='isppol     ='
   nm_counter(13,2)='ikpt       ='

   nm_levels(14)   ='dfpt_vtowfk         '
   ncounter(14)=2
   nm_counter(14,1)='inonsc     ='
   nm_counter(14,2)='iband      ='

   nm_levels(15)   ='dfpt_cgwf           '
   ncounter(15)=1
   nm_counter(15,1)='iline      ='

   nm_levels(16)   ='getgh1c             '
   ncounter(16)=0
   nm_counter(16,1)='            '

   nm_levels(17)   ='dfpt_rhofermi       '
   ncounter(17)=2
   nm_counter(17,1)='isppol     ='
   nm_counter(17,2)='ikpt       ='

   nm_levels(18)   ='dfpt_wfkfermi       '
   ncounter(18)=2
   nm_counter(18,1)='inonsc     ='
   nm_counter(18,2)='iband      ='


!  Optdriver=2
   nm_levels(20)   ='suscep              '
   ncounter(20)=0
   nm_counter(20,1)='            '


!  Optdriver=3
   nm_levels(30)   ='screening           '
   ncounter(30)=1
   nm_counter(30,1)='iqpt       ='

!  Optdriver=4
   nm_levels(40)   ='sigma               '
   ncounter(40)=1
   nm_counter(40,1)='ikpt_gw    ='

!  Optdriver=5
   nm_levels(50)   ='nonlinear           '
   ncounter(50)=0
   nm_counter(50,1)='            '

   nm_levels(51)   ='dfptnl_loop         '
   ncounter(51)=2
   nm_counter(51,1)='pert1case  ='
   nm_counter(51,2)='pert3case  ='

   nm_levels(52)   ='mv_/dfptnl_resp     '
   ncounter(52)=2
   nm_counter(52,2)='ikpt       ='

!  Optdriver=9
   nm_levels(90)   ='bethe_salpether     '
   ncounter(90)=0
   nm_counter(90,1)='            '

   if(istat<0)then
     write(message, '(a,i7,a,a,a,a,a)' )&
&     'the value of the input variable istatr is',istat,' ,',ch10,&
&     'while it must be a positive number.',ch10,&
&     'Action : change istatr in your input file.'
     MSG_ERROR(message)
   end if
   output_rate=istat
 end if

!The true input variable "shift_rate" is only available at the second call
 if(statnu==1)then
   if(istat<0 .or. istat>output_rate)then
     write(message, '(a,i7,3a,i7,2a)' )&
&     'the value of the input variable istatshft is',istat,' ,',ch10,&
&     'while it must be a positive number smaller or equal to istatr=',output_rate,ch10,&
&     'Action: change istatshft in your input file.'
     MSG_ERROR(message)
   end if
   shift_rate=istat
   if(shift_rate==output_rate)shift_rate=0

!  At this second call, also feed information that the abinit routine called ab7_invars_load
   write(unit=nm_routine(1),fmt='(a20)') 'call ab7_invars_load'
   active(1)=1
 end if

!Check the value of level
 if( minval(abs(list_level(:)-level)) /= 0)then
   write(message, '(a,i5,3a)' )&
&   '  The value of level in the calling routine is',level,' ,',ch10,&
&   '  which is not an allowed value.'
   MSG_BUG(message)
 end if

!Assign the info about the actual routine
 write(unit=nm_routine(level),fmt='(a20)') trim(adjustl(routine))
 if(trim(adjustl(nm_routine(level)))=='exit')then
!  The value 2 will be changed to 0 at the end of the routine.
   active(level)=2
 else if(trim(adjustl(nm_routine(level)))=='')then
   active(level)=0
 else
   active(level)=1
 end if

!Assign the info about the actual counter
 if(counter>=0)then
   if(ncounter(level)==1)then
     actual_counter(level,1)=counter
   else if(ncounter(level)==2)then
     actual_counter(level,2)=counter/100
!    The counter number 1 does not allow more than 99 passes
     actual_counter(level,1)=counter-(counter/100)*100
   end if
 end if

!============================================================

!After treatment of present information, output of the status
 statnu=statnu+1

!DEBUG
! write(std_out,*)' status : statnu, output_rate, shift_rate=',statnu,output_rate, shift_rate
! write(std_out,*)'level,routine=',level,routine
! write(std_out,*)'active(level)=',active(level)
! write(std_out,*)'counter,actual_counter(level,1:2)=',counter,actual_counter(level,1:2)
! write(std_out,*)'List of active levels :'
! do ilevel=1,mlevel
!   if(active(ilevel)/=0)write(std_out,*)' Active level number=',ilevel
! end do
!ENDDEBUG

 if(statnu > 2 )then
   if( mod(statnu,output_rate)==shift_rate )then

     if (open_file(filstat,message,newunit=temp_unit,form='formatted',status='unknown') /= 0) then
       MSG_ERROR(message)
     end if

     rewind temp_unit
     write(temp_unit,*)

     headwr='(a,i4,a,i6 )'
     if(statnu>=100000)   headwr='(a,i4,a,i9 )'
     if(output_rate>=1000)headwr='(a,i6,a,i6 )'
     if(statnu>=100000 .and. output_rate>=1000)   headwr='(a,i6,a,i9 )'
     if(statnu>=100000000)headwr='(a,i6,a,i12)'
     write(temp_unit,headwr)' Status file, with repetition rate',output_rate,', status number',statnu
     write(temp_unit,*)

!    Treat every level one after the other
     do ilevel=1,mlevel
!      This level must be activated in order to have a corresponding output
       if(active(ilevel)==1 .or. active(ilevel)==2)then

         write(temp_unit,'(4a)')'  Level ',nm_levels(ilevel),' : ',adjustl(nm_routine(ilevel))

!        Is there a counter for this level ?
         if(ncounter(ilevel)>=1)then

           if(actual_counter(ilevel,1)>0)then
             write(temp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,1),actual_counter(ilevel,1)
           end if
           if(ncounter(ilevel)==2)then
             if(actual_counter(ilevel,2)>0)then
               write(temp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,2),actual_counter(ilevel,2)
             end if
           end if

         end if
       end if ! End of the check on activation of the level
     end do ! End of the loop on the levels

     close (temp_unit)
   end if ! End of the repetition rate check
 end if ! statnu > 2

 if (active(level)==2)then
   active(level)=0
   nm_routine(level)='              '
 end if

 call timab(73,2,tsec)

end subroutine status
!!***

end module m_dtfil
!!***
