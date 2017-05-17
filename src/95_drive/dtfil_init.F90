!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtfil_init
!!
!! NAME
!! dtfil_init
!!
!! FUNCTION
!! Initialize most of the dtfil structured variable
!! (what is left should be initialized inside the itimimage,
!! iimage and itime loops).
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset,&
&                      image_index) ! optional argument

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fstrings, only : int2char4

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init'
 use interfaces_32_util
 use interfaces_54_abiutil
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
 integer,parameter :: unwff1=1,unwff2=2,unwffgs=3,unwfkq=4,unwft1=11
 integer,parameter :: unwft2=12,unwftgs=13,unwftkq=14,unylm=24,unylm1=25
 integer,parameter :: unkss=40,unscr=41,unqps=43
 integer :: ii,iimage,ireadden,ireadkden,ireadwf,ixx,jdtset,will_read
 character(len=10) :: appen,tag
 character(len=9) :: stringvar
 character(len=15) :: stringfile
 character(len=500) :: message
 character(len=fnlen) :: filsus,filddbsin,fildens1in,fildensin,filpawdensin,filkdensin,filqps,filscr
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
 dtfil%unwffgs=unwffgs
 dtfil%unwffkq=unwfkq
 dtfil%unwft1 =unwft1
 dtfil%unwft2 =unwft2
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
