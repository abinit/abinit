!!****m* ABINIT/m_dtfil
!! NAME
!!  m_dtfil
!!
!! FUNCTION
!!   object and procedures dealing with input/output filenames
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (XG, MT)
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
 use m_abicore
 use m_errors
 use m_xmpi
 use m_build_info
 use m_dtset

 use defs_abitypes,  only : MPI_type
 use m_clib,         only : clib_rename
 use m_fstrings,     only : int2char4, rmquotes, sjoin, strcat, basename
 use m_io_tools,     only : open_file, file_exists
 use m_libpaw_tools, only : libpaw_log_flag_set
 use m_parser,       only : parsefile, intagm

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_dtfil/datafiles_type
!! NAME
!! datafiles_type
!!
!! FUNCTION
!! The datafiles_type structures datatype gather all the variables
!! related to files, such as filename, and file units.
!! For one dataset, it is initialized in 95_drive/dtfil_init1.F90,
!! and will not change at all during the treatment of the dataset.
!!
!! SOURCE

 type, public :: datafiles_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! These keywords are only used in algorithms using images of the cell
  integer :: getwfk_from_image
   ! index of image from which read WFK file (0 if standard WFK)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getden_from_image
   ! index of image from which read DEN file (0 if standard DEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: getpawden_from_image
   ! index of image from which read PAWDEN file (0 if standard PAWDEN)
   !    -1: the same image as current one
   !     0: no image
   !    >0: index of an image

  integer :: ireadddb
   ! ireadddb non-zero  if the ddb file must be read

  integer :: ireadden
   ! ireadden non-zero  if the den file must be read

  integer :: ireadkden
   ! ireadkden non-zero  if the kden file must be read

  integer :: ireadwf
   ! if(optdriver/=1), that is, no response-function computation,
   !   ireadwf non-zero  if the wffk file must be read
   !   (if irdwfk non-zero or getwfk non-zero)
   ! if(optdriver==1), that is, response-function computation,
   !   ireadwf non-zero  if the wff1 file must be read
   !   (if ird1wf non-zero or get1wf non-zero)

  integer :: unchi0  ! unit number for chi0 files
  integer :: unddb   ! unit number for Derivative DataBase
  integer :: unddk   ! unit number for ddk 1WF file
  integer :: undkdk  ! unit number for 2WF file (dkdk)
  integer :: undkde  ! unit number for 2WF file (dkde)
  integer :: unkg    ! unit number for k+G data
  integer :: unkgq   ! unit number for k+G+q data
  integer :: unkg1   ! unit number for first-order k+G+q data
  integer :: unkss   ! unit number for KSS file
  integer :: unqps   ! unit number for QPS file
  integer :: unscr   ! unit number for SCR file
  integer :: unwff1  ! unit number for wavefunctions, number one
  integer :: unwff2  ! unit number for wavefunctions, number two
  integer :: unwff3  ! unit number for wavefunctions, number three
  integer :: unwffgs ! unit number for ground-state wavefunctions
  integer :: unwffkq ! unit number for k+q ground-state wavefunctions
  integer :: unwft1  ! unit number for wavefunctions, temporary one
  integer :: unwft2  ! unit number for wavefunctions, temporary two
  integer :: unwft3  ! unit number for wavefunctions, temporary three
  integer :: unwftgs ! unit number for ground-state wavefunctions, temporary
  integer :: unwftkq ! unit number for k+q ground-state wavefunctions, temporary
  integer :: unylm   ! unit number for Ylm(k) data
  integer :: unylm1  ! unit number for first-order Ylm(k+q) data
  integer :: unpaw   ! unit number for temporary PAW data (for ex. rhoij_nk) (Paw only)
  integer :: unpaw1  ! unit number for temporary PAW first-order cprj1=<c1_k,q|p>(1) data
  integer :: unpawq  ! unit number for temporary PAW cprjq=<c+_k+q|p> at k+qdata
  integer :: unpos   ! unit number for restart molecular dynamics

  ! TODO: All this strings should be initialized with ABI_NOFILE
  ! so that we can easily test for path /= ABI_NOFILE instead of getwfk /= 0 or irdwfk /= 0

  character(len=fnlen) :: filnam_ds(5)
   ! if no dataset mode, the five names from the standard input:
   !   ab_in, ab_out, abi, abo, tmp
   ! if dataset mode, the same 5 filenames, appended with //'_DS'//trim(jdtset)

  character(len=fnlen) :: filddbsin
   ! if no dataset mode             : abi//'DDB'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DDB'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DDB'

  character(len=fnlen) :: fildensin
   ! if no dataset mode             : abi//'DEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'DEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'DEN'

  character(len=fnlen) :: fildvdbin
   ! if no dataset mode             : abi//'DVDB'
   ! if dataset mode, and getdvdb==0 : abi//'_DS'//trim(jdtset)//'DVDB'
   ! if dataset mode, and getdvdb/=0 : abo//'_DS'//trim(jgetden)//'DVDB'

  character(len=fnlen) :: filpotin
   ! Filename used to read POT file.
   ! Initialize via getpot_filepath

  character(len=fnlen) :: filkdensin
   ! if no dataset mode             : abi//'KDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'KDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'KDEN'

  character(len=fnlen) :: filpawdensin
   ! if no dataset mode             : abi//'PAWDEN'
   ! if dataset mode, and getden==0 : abi//'_DS'//trim(jdtset)//'PAWDEN'
   ! if dataset mode, and getden/=0 : abo//'_DS'//trim(jgetden)//'PAWDEN'

  character(len=fnlen) :: filsigephin
   ! Filename used to read SIGEPH file.
   ! Initialize via getsigeph_filepath

  character(len=fnlen) :: filstat
   ! tmp//'_STATUS'

  character(len=fnlen) :: fnamewffk
   ! the name of the ground-state wavefunction file to be read (see driver.F90)

  character(len=fnlen) :: fnamewffq
   ! the name of the k+q ground-state wavefunction file to be read (see driver.F90)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffddk
   ! the generic name of the ddk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdelfd
   ! the generic name of the electric field response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkdk
   ! the generic name of the 2nd order dkdk response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewffdkde
   ! the generic name of the 2nd order dkde response wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fnamewff1
   ! the generic name of the first-order wavefunction file(s) to be read (see driver.F90)
   ! (the final name is formed by appending the number of the perturbation)
   ! only useful in the response-function case

  character(len=fnlen) :: fildens1in   ! to be described by MVeithen
  character(len=fnlen) :: fname_tdwf
  character(len=fnlen) :: fname_w90

  character(len=fnlen) :: fnametmp_wf1
  character(len=fnlen) :: fnametmp_wf2
  character(len=fnlen) :: fnametmp_1wf1
  character(len=fnlen) :: fnametmp_1wf2
  character(len=fnlen) :: fnametmp_wfgs
  character(len=fnlen) :: fnametmp_wfkq
   ! Set of filenames formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _WF1, _WF2 ...
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_kg
  character(len=fnlen) :: fnametmp_kgq
  character(len=fnlen) :: fnametmp_kg1
  character(len=fnlen) :: fnametmp_dum
  character(len=fnlen) :: fnametmp_ylm
  character(len=fnlen) :: fnametmp_ylm1
  character(len=fnlen) :: fnametmp_paw
  character(len=fnlen) :: fnametmp_paw1
  character(len=fnlen) :: fnametmp_pawq
   ! Set of filenames formed from trim(dtfil%filnam_ds(5))//APPEN where APPEN is _KG, _DUM, followed
   ! by the index of the processor.
   ! See dtfil_init

  character(len=fnlen) :: fnametmp_cg
  character(len=fnlen) :: fnametmp_cprj
  character(len=fnlen) :: fnametmp_eig
  character(len=fnlen) :: fnametmp_1wf1_eig
  character(len=fnlen) :: fnametmp_fft
  character(len=fnlen) :: fnametmp_fft_mgga
  character(len=fnlen) :: fnametmp_kgs
  character(len=fnlen) :: fnametmp_sustr
  character(len=fnlen) :: fnametmp_tdexcit
  character(len=fnlen) :: fnametmp_tdwf

!@Bethe-Salpeter
! New files introduced for the Bethe-Salpeter part.

   character(len=fnlen) :: fnameabi_bsham_reso
    ! if no dataset mode             : abi//'BSR'
    ! if dataset mode, and getbsreso==0 : abi//'_DS'//trim(jdtset)//'BSR'
    ! if dataset mode, and getbsreso/=0 : abo//'_DS'//trim(jget_reso_bsham)//'BSR'

   character(len=fnlen) :: fnameabi_bsham_coup
    ! if no dataset mode             : abi//'BSC'
    ! if dataset mode, and getbscoup==0 : abi//'_DS'//trim(jdtset)//'BSC'
    ! if dataset mode, and getbscoup/=0 : abo//'_DS'//trim(jget_coup_bsham)//'BSC'

  character(len=fnlen) :: fnameabi_bseig
   ! The name of the file containing the eigenstates and eigenvalues of the Bethe-Salpeter Hamiltonian
   ! if no dataset mode             : abi//'BS_EIG'
   ! if dataset mode, and getbseig==0 : abi//'_DS'//trim(jdtset)//'BS_EIG'
   ! if dataset mode, and getbseig/=0 : abo//'_DS'//trim(jget_bseig)//'BS_EIG'

   character(len=fnlen) :: fnameabi_haydock
   ! The prefix used to construct the names of the files containing the coefficients of the
   ! continued fractions produced by the Haydock iterative algorithm.
   ! if no dataset mode             : abi//'HAYDOCK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'HAYDOCK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'HAYDOCK'

   character(len=fnlen) :: fnameabi_wfkfine
   ! The name of the file containing the wavefunctions on a fine grid
   ! if no dataset mode             : abi//'WFK'
   ! if dataset mode, and gethaydock==0 : abi//'_DS'//trim(jdtset)//'WFK'
   ! if dataset mode, and gethaydock/=0 : abo//'_DS'//trim(jget_bseig)//'WFK'

!END @BEthe-Salpeter

!The following filenames do not depend on itimimage, iimage and itime loops.
!Note the following convention:
!  fnameabo_* are filenames used for ouput results (writing)
!  fnameabi_* are filenames used for data that should be read by the code.
!  fnametmp_* are filenames used for temporary files that should be erased at the end of each dataset.
!
!If a file does not have the corresponding "abi" or the corresponding "abo" name, that means that
!that particular file is only used for writing or for reading results, respectively.

  character(len=fnlen) :: fnameabi_efmas
  character(len=fnlen) :: fnameabi_hes
  character(len=fnlen) :: fnameabi_phfrq
  character(len=fnlen) :: fnameabi_phvec
  character(len=fnlen) :: fnameabi_qps
  character(len=fnlen) :: fnameabi_scr            ! SCReening file (symmetrized inverse dielectric matrix)
  character(len=fnlen) :: fnameabi_sus            ! KS independent-particle polarizability file
  character(len=fnlen) :: fnameabo_ddb
  character(len=fnlen) :: fnameabo_den
  character(len=fnlen) :: fnameabo_dos
  character(len=fnlen) :: fnameabo_dvdb
  character(len=fnlen) :: fnameabo_eelf
  character(len=fnlen) :: fnameabo_eig
  character(len=fnlen) :: fnameabo_eigi2d
  character(len=fnlen) :: fnameabo_eigr2d
  character(len=fnlen) :: fnameabo_em1
  character(len=fnlen) :: fnameabo_em1_lf
  character(len=fnlen) :: fnameabo_em1_nlf
  character(len=fnlen) :: fnameabo_fan
  character(len=fnlen) :: fnameabo_gkk
  character(len=fnlen) :: fnameabo_gw
  character(len=fnlen) :: fnameabo_gw_nlf_mdf
  character(len=fnlen) :: fnameabo_kss
  character(len=fnlen) :: fnameabo_moldyn
  character(len=fnlen) :: fnameabo_pot
  character(len=fnlen) :: fnameabo_qps            ! Quasi-Particle band structure file.
  character(len=fnlen) :: fnameabo_qp_den
  character(len=fnlen) :: fnameabo_qp_pawden      ! Full QP density
  character(len=fnlen) :: fnameabo_qp_dos
  character(len=fnlen) :: fnameabo_qp_eig
  character(len=fnlen) :: fnameabo_rpa
  character(len=fnlen) :: fnameabo_rpa_nlf_mdf
  character(len=fnlen) :: fnameabo_scr
  character(len=fnlen) :: fnameabo_sgm
  character(len=fnlen) :: fnameabo_sgr
  character(len=fnlen) :: fnameabo_sig
  character(len=fnlen) :: fnameabo_spcur
  character(len=fnlen) :: fnameabo_sus
  character(len=fnlen) :: fnameabo_vha
  character(len=fnlen) :: fnameabo_vpsp
  character(len=fnlen) :: fnameabo_vso
  character(len=fnlen) :: fnameabo_vxc
  character(len=fnlen) :: fnameabo_wan
  character(len=fnlen) :: fnameabo_wfk
  character(len=fnlen) :: fnameabo_wfq
  character(len=fnlen) :: fnameabo_w90
  character(len=fnlen) :: fnameabo_1wf
  character(len=fnlen) :: fnameabo_gwdiag
  character(len=fnlen) :: fnameabo_nlcc_derivs
  character(len=fnlen) :: fnameabo_pspdata

!The following filenames are initialized only iniside itimimage, iimage and itime loops,
!and are appended with the adequate specifier 'app'.

  character(len=fnlen) :: fnameabo_app
  character(len=fnlen) :: fnameabo_app_atmden_core
  character(len=fnlen) :: fnameabo_app_atmden_full
  character(len=fnlen) :: fnameabo_app_atmden_val
  character(len=fnlen) :: fnameabo_app_n_tilde
  character(len=fnlen) :: fnameabo_app_n_one
  character(len=fnlen) :: fnameabo_app_nt_one
  character(len=fnlen) :: fnameabo_app_bxsf
  character(len=fnlen) :: fnameabo_app_cif
  character(len=fnlen) :: fnameabo_app_den
  character(len=fnlen) :: fnameabo_app_dos
  character(len=fnlen) :: fnameabo_app_elf
  character(len=fnlen) :: fnameabo_app_elf_down
  character(len=fnlen) :: fnameabo_app_elf_up
  character(len=fnlen) :: fnameabo_app_eig
  character(len=fnlen) :: fnameabo_app_fatbands
  character(len=fnlen) :: fnameabo_app_gden1
  character(len=fnlen) :: fnameabo_app_gden2
  character(len=fnlen) :: fnameabo_app_gden3
  character(len=fnlen) :: fnameabo_app_geo
  character(len=fnlen) :: fnameabo_app_kden
  character(len=fnlen) :: fnameabo_app_lden
  character(len=fnlen) :: fnameabo_app_nesting
  character(len=fnlen) :: fnameabo_app_pawden
  character(len=fnlen) :: fnameabo_app_pot
  character(len=fnlen) :: fnameabo_app_opt
  character(len=fnlen) :: fnameabo_app_opt2
  character(len=fnlen) :: fnameabo_app_stm
  character(len=fnlen) :: fnameabo_app_vclmb
  character(len=fnlen) :: fnameabo_app_vha
  character(len=fnlen) :: fnameabo_app_vhxc
  character(len=fnlen) :: fnameabo_app_vhpsp
  character(len=fnlen) :: fnameabo_app_vpsp
  character(len=fnlen) :: fnameabo_app_vxc
  character(len=fnlen) :: fnameabo_app_wfk
  character(len=fnlen) :: fnameabo_app_1dm
  character(len=fnlen) :: fnameabo_app_vha_1dm
  character(len=fnlen) :: fnameabo_app_vclmb_1dm
  character(len=fnlen) :: fnametmp_app_den
  character(len=fnlen) :: fnametmp_app_kden

 end type datafiles_type
!!***

 public :: dtfil_init
 public :: dtfil_init_img
 public :: dtfil_init_time
 public :: mkfilename
 public :: isfile
 public :: iofn1
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
!! jdtset_(0:ndtset)=actual index of the datasets
!! mpi_enreg=information about MPI parallelization
!! ndtset=number of datasets
!! [image_index]= index of image to be used when appending
!!             "_IMGxxx" string to file names. To be used only when an algorithm
!!             using images of the cell is activated
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
!!      m_driver,m_gstateimg
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine dtfil_init(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset,&
&                      image_index) ! optional argument

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
! TODO: Remove all these units and use get_unit API
 integer,parameter :: unchi0=42,unddb=16,unddk=50,undkdk=54,undkde=55,unkg1=19,unkg=17,unkgq=18
 integer,parameter :: unpaw=26,unpaw1=27,unpawq=28,unpos=30
 integer,parameter :: unwff1=1,unwff2=2,unwff3=8,unwffgs=3,unwfkq=4,unwft1=11
 integer,parameter :: unwft2=12,unwft3=15,unwftgs=13,unwftkq=14,unylm=24,unylm1=25
 integer,parameter :: unkss=40,unscr=41,unqps=43
 integer :: ii,iimage,ireadden,ireadkden,ireadwf,ixx,jdtset,will_read
 character(len=10) :: appen,tag
 character(len=9) :: stringvar
 character(len=15) :: stringfile
 character(len=500) :: msg
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

 ! If multi dataset mode, special treatment of filenames 3 and 4 (density and
 ! wavefunctions input and output, as well as other output files)
 if(ndtset>0)then
   call appdig(jdtset,'',appen)
   filnam_ds(3)=trim(filnam(3))//'_DS'//trim(appen)
   filnam_ds(4)=trim(filnam(4))//'_DS'//trim(appen)
 end if

 ! If multi image mode (nimage>1), special treatment of filenames 4 and 5
 if(iimage>0)then
   call appdig(iimage,'',appen)
   filnam_ds(4)=trim(filnam_ds(4))//'_IMG'//trim(appen)
   filnam_ds(5)=trim(filnam_ds(5))//'_IMG'//trim(appen)
 end if

 ! According to getwfk and irdwfk, build _WFK file name, referred as fnamewffk
 if (iimage >0 .and. dtfil%getwfk_from_image /= 0) then
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
 call mkfilename(filnam,fnamewffk,dtset%getwfk,idtset,dtset%irdwfk,jdtset_,ndtset,stringfile,stringvar,will_read, &
                 getpath=dtset%getwfk_filepath)

 if (dtset%optdriver /= RUNL_RESPFN) ireadwf = will_read
 if(ndtset/=0 .and. dtset%optdriver==RUNL_RESPFN .and. will_read==0)then
   write(msg, '(5a,i0,3a,i0,a,i0,3a)' )&
   'At least one of the input variables irdwfk and getwfk ',ch10,&
   'must refer to a valid _WFK file, in the response function',ch10,&
   'case, while for idtset = ',idtset,',',ch10,&
   'they are irdwfk= ',dtset%irdwfk,', and getwfk= ',dtset%getwfk,'.',ch10,&
   'Action: correct irdwfk or getwfk in your input file.'
   MSG_ERROR(msg)
 end if

!Treatment of the other get wavefunction variable, if response function case or nonlinear case
 if (ANY(dtset%optdriver == [RUNL_RESPFN, RUNL_NONLINEAR, RUNL_EPH, RUNL_LONGWAVE])) then

   ! According to getwfq and irdwfq, build _WFQ file name, referred as fnamewffq
   stringfile='_WFQ' ; stringvar='wfq'
   call mkfilename(filnam,fnamewffq,dtset%getwfq,idtset,dtset%irdwfq,jdtset_,ndtset,stringfile,stringvar,will_read, &
                   getpath=dtset%getwfq_filepath)
   ! If fnamewffq is not initialized thanks to getwfq or irdwfq, use fnamewffk
   if(will_read==0) fnamewffq = fnamewffk

   ! According to get1wf and ird1wf, build _1WF file name, referred as fnamewff1
   stringfile='_1WF' ; stringvar='1wf'
   call mkfilename(filnam,fnamewff1,dtset%get1wf,idtset,dtset%ird1wf,jdtset_,ndtset,stringfile,stringvar,will_read)
   ireadwf=will_read

   ! According to getddk and irdddk, build _1WF file name, referred as fnamewffddk
   stringfile='_1WF' ; stringvar='ddk'
   call mkfilename(filnam,fnamewffddk,dtset%getddk,idtset,dtset%irdddk,jdtset_,ndtset,stringfile,stringvar,will_read)

   ! According to getdelfd, build _1WF file name, referred as fnamewffdelfd
   stringfile='_1WF' ; stringvar='delfd'
   call mkfilename(filnam,fnamewffdelfd,dtset%getdelfd,idtset,0,jdtset_,ndtset,stringfile,stringvar,will_read)

   ! According to getdkdk, build _1WF file name, referred as fnamewffdkdk
   stringfile='_1WF' ; stringvar='dkdk'
   call mkfilename(filnam,fnamewffdkdk,dtset%getdkdk,idtset,0,jdtset_,ndtset,stringfile,stringvar,will_read)

   ! According to getdkde, build _1WF file name, referred as fnamewffdkde
   stringfile='_1WF' ; stringvar='dkde'
   call mkfilename(filnam,fnamewffdkde,dtset%getdkde,idtset,0,jdtset_,ndtset,stringfile,stringvar,will_read)
 end if

!-------------------------------------------------------------------------------------------
 ! Build name of files from dtfil%filnam_ds(3)

 ! According to getddb, build _DDB file name, referred as filddbsin
 stringfile='_DDB'; stringvar='ddb'
 call mkfilename(filnam,filddbsin,dtset%getddb,idtset,dtset%irdddb,jdtset_,ndtset,stringfile,stringvar,will_read, &
                  getpath=dtset%getddb_filepath)

 ! According to getpot, build _POT file name
 stringfile='_POT'; stringvar='pot'
 call mkfilename(filnam, dtfil%filpotin, 0, idtset, 0, jdtset_, ndtset, stringfile, stringvar, will_read, &
                  getpath=dtset%getpot_filepath)

 ! According to getdvdb, build _DVDB file name
 stringfile='_DVDB'; stringvar='dvdb'
 call mkfilename(filnam,dtfil%fildvdbin,dtset%getdvdb,idtset,dtset%irddvdb,jdtset_,ndtset,stringfile,stringvar,will_read, &
                  getpath=dtset%getdvdb_filepath)
 if (will_read == 0) dtfil%fildvdbin = ABI_NOFILE

 ! According to getsigeph_filepath, build _SIGEPH file name
 stringfile='_SIGEPH.nc'; stringvar='sigeph'
 call mkfilename(filnam, dtfil%filsigephin, 0, idtset, 0, jdtset_, ndtset, stringfile, stringvar, will_read, &
                  getpath=dtset%getsigeph_filepath)
 ! If getsigeph_filepath is not used, will read the output as assumed in the transport driver when called after sigeph
 if (will_read == 0) dtfil%filsigephin = strcat(filnam_ds(4), "_SIGEPH.nc")

 ! According to getden, build _DEN file name, referred as fildensin
 ! A default is available if getden is 0
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
 call mkfilename(filnam,fildensin,dtset%getden,idtset,dtset%irdden,jdtset_,ndtset,stringfile,stringvar, will_read, &
                 getpath=dtset%getden_filepath)

 if(will_read==0)fildensin=trim(filnam_ds(3))//'_DEN'
 ireadden=will_read

 if ((dtset%optdriver==RUNL_GWLS.or.dtset%optdriver==RUNL_GSTATE) .and.dtset%iscf<0) ireadden=1

 ! According to getpawden, build _PAWDEN file name, referred as filpawdensin
 ! A default is available if getden is 0
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
 call mkfilename(filnam,filpawdensin,dtset%getpawden,idtset,dtset%irdden,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filpawdensin=trim(filnam_ds(3))//'_PAWDEN'

 ! According to getden and usekden, build _KDEN file name, referred as filkdensin
 ! A default is available if getden is 0
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
   call mkfilename(filnam,filkdensin,dtset%getden,idtset,dtset%irdden,jdtset_,ndtset,stringfile,stringvar,will_read)
   if(will_read==0)filkdensin=trim(filnam_ds(3))//'_KDEN'
   ireadkden=will_read
   if ((dtset%optdriver==RUNL_GSTATE.or.dtset%optdriver==RUNL_GWLS).and.dtset%iscf<0) ireadkden=1
 else
   ireadkden=0
 end if

 ! According to get1den, build _DEN file name, referred as fildens1in
 ! A default is available if get1den is 0
 stringfile='_DEN' ; stringvar='1den'
 call mkfilename(filnam,fildens1in,dtset%get1den,idtset,dtset%ird1den,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fildens1in=trim(filnam_ds(3))//'_DEN'

 ! According to getefmas and irdefmas, build _EFMAS file name, referred as fil_efmas
 ! A default is available if getefmas is 0
 stringfile='_EFMAS.nc' ; stringvar='efmas'
 call mkfilename(filnam,fil_efmas,dtset%getefmas,idtset,dtset%irdefmas,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fil_efmas=trim(filnam_ds(3))//'_EFMAS.nc'

 ! According to getscr and irdscr, build _SCR file name, referred as filscr
 ! A default is available if getscr is 0
 stringfile='_SCR' ; stringvar='scr'
 call mkfilename(filnam,filscr,dtset%getscr,idtset,dtset%irdscr,jdtset_,ndtset,stringfile,stringvar,will_read, &
                 getpath=dtset%getscr_filepath)
 if(will_read==0)filscr=trim(filnam_ds(3))//'_SCR'

 ! According to getsuscep and irdsuscep, build _SUS file name, referred as filsus
 ! A default is available if getsuscep is 0
 stringfile='_SUS' ; stringvar='sus'
 call mkfilename(filnam,filsus,dtset%getsuscep,idtset,dtset%irdsuscep,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filsus=TRIM(filnam_ds(3))//'_SUS'

 ! According to getqps and irdqps, build _QPS file name, referred as filqps
 ! A default is available if getqps is 0
 stringfile='_QPS' ; stringvar='qps'
 call mkfilename(filnam,filqps,dtset%getqps,idtset,dtset%irdqps,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filqps=trim(filnam_ds(3))//'_QPS'

 ! According to getbseig and irdbseig, build _BSEIG file name, referred as filbseig
 ! A default is available if getbseig is 0
 stringfile='_BSEIG' ; stringvar='bseig'
 call mkfilename(filnam,filbseig,dtset%getbseig,idtset,dtset%irdbseig,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filbseig=trim(filnam_ds(3))//'_BSEIG'

 ! According to gethaydock and irdhaydock, build _HAYD file name, referred as filhaydock.
 ! A default is available if gethaydock is 0
 stringfile='_HAYDR_SAVE' ; stringvar='haydock'
 call mkfilename(filnam,filhaydock,dtset%gethaydock,idtset,dtset%irdhaydock,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)filhaydock=trim(filnam_ds(3))//'_HAYDR_SAVE'

 ! According to getbsr and irdbsr, build _BSR file name, referred as fil_bsreso
 ! A default is available if getbsr is 0
 stringfile='_BSR' ; stringvar='bsreso'
 call mkfilename(filnam,fil_bsreso,dtset%getbsreso,idtset,dtset%irdbsreso,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0) fil_bsreso=trim(filnam_ds(3))//'_BSR'

 ! According to getbsc and irdbsc, build _BSC file name, referred as fil_bscoup
 ! A default is available if getbsc is 0
 stringfile='_BSC' ; stringvar='bscoup'
 call mkfilename(filnam,fil_bscoup,dtset%getbscoup,idtset,dtset%irdbscoup,jdtset_,ndtset,stringfile,stringvar,will_read)
 if(will_read==0)fil_bscoup=trim(filnam_ds(3))//'_BSC'

 ! According to getwfkfine and irdwfkfine, build _WFK file name, referred as filwfkfine
 ! A default is avaible if getwfkfine is 0
 stringfile='_WFK' ; stringvar='wfkfine'
 call mkfilename(filnam,filwfkfine,dtset%getwfkfine,idtset,dtset%irdwfkfine,jdtset_,ndtset,stringfile,stringvar,will_read, &
                 getpath=dtset%getwfkfine_filepath)
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
 ! Build name of files from dtfil%filnam_ds(4)

 dtfil%fnameabo_ddb=trim(dtfil%filnam_ds(4))//'_DDB'
 dtfil%fnameabo_den=trim(dtfil%filnam_ds(4))//'_DEN'
 dtfil%fnameabo_dos=trim(dtfil%filnam_ds(4))//'_DOS'
 dtfil%fnameabo_dvdb=trim(dtfil%filnam_ds(4))//'_DVDB'
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
 ! Build name of files from dtfil%filnam_ds(5)
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

 if(xmpi_paral==1)then ! parallel case: the index of the processor must be appended
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

 ! Prepare the name of the _FFT file
 filfft=trim(dtfil%filnam_ds(5))//'_FFT'
 if(xmpi_paral==1 .or. mpi_enreg%paral_kgb==1)then
   call int2char4(mpi_enreg%me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   filfft=trim(filfft)//'_P-'//trim(tag)
 end if
 dtfil%fnametmp_fft=filfft
 dtfil%fnametmp_fft_mgga=trim(filfft)//'_MGGA'

 ! These keywords are only used in algorithms using images of the cell
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
!!      m_gstate,m_mover
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine dtfil_init_time(dtfil,iapp)

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
!! [suff]= --optional argument--indicates the suffixe to be appended:
!!         SUF=TIM (default) or SUF=IMG or ...
!!
!! OUTPUT
!! filapp= filename with appended string
!!
!! PARENTS
!!      m_dtfil
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine fappnd(filapp,filnam,iapp,&
&                 suff) ! optional argument

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
        'Requested file name extension has more than the allowed 8 digits.',ch10,&
        'Action: resubmit the job with smaller value for ntime.',ch10,&
        'Value computed here was ndig=',ndig,ch10,&
        'iapp= ',iapp,' filnam= ',trim(filnam)
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
!!      m_driver
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset,ndtset,ndtset_alloc)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: idtset,ndtset,ndtset_alloc
 type(datafiles_type),intent(out) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: jdtset(0:ndtset)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------
!scalars
 integer :: iget

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
!! filnam(5)=the root file names (only filnam(3) and filnam(4) are really needed)
!! get=input 'get variable', if 1, must get the file from another dataset
!! idtset=number of the dataset
!! ird=input 'iread variable', if 1, must get the file from the input root
!! jdtset_(0:ndtset)=actual index of the dataset
!! ndtset=number of datasets
!! stringfil=the string of characters to be appended e.g. '_WFK' or '_DEN'
!! stringvar=the string of characters to be appended
!!   that defines the 'get' or 'ird' variables, e.g. 'wfk' or 'ddk'
!! [getpath]=String with filename to be used as input, exclude get and ird option.
!!
!! OUTPUT
!! filnam_out=the new file name
!! will_read=1 if the file must be read ; 0 otherwise (ird and get were zero)
!!
!! PARENTS
!!      m_dtfil,m_mpi_setup
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine mkfilename(filnam,filnam_out,get,idtset,ird,jdtset_,ndtset,stringfil,stringvar,will_read, &
                      getpath) ! Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: get,idtset,ird,ndtset
 integer,intent(out) :: will_read
 character(len=*),intent(in) :: stringfil
 character(len=*),intent(in) :: stringvar
 character(len=fnlen),intent(out) :: filnam_out
 character(len=fnlen),optional,intent(in) :: getpath
!arrays
 integer,intent(in) :: jdtset_(0:ndtset)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
!scalars
 integer :: jdtset,jget
 character(len=4) :: appen
 character(len=500) :: msg
 character(len=fnlen) :: filnam_appen

! *************************************************************************

 ! Here, defaults if no get variable
 will_read=ird

 filnam_appen = trim(filnam(3))
 if (ndtset > 0) then
   jdtset = jdtset_(idtset)
   call appdig(jdtset, '', appen)
   filnam_appen = trim(filnam_appen)//'_DS'//appen
 end if
 filnam_out = trim(filnam_appen)//trim(stringfil)

 if (present(getpath)) then
   if (getpath /= ABI_NOFILE) then
     if (ird /= 0 .or. get /= 0)then
       write(msg, '(11a,i0,3a,i0,a,i0,7a)' ) &
       'When the input variable: ', trim(getpath), ' is used ',ch10, &
       'the input variables ird',trim(stringvar),' and get',trim(stringvar),' cannot be',ch10,&
       'simultaneously non-zero, while for idtset = ',idtset,',',ch10,&
       'they are ',ird,', and ',get,'.',ch10,&
       'Action: correct ird',trim(stringvar),' or get',trim(stringvar),' in your input file.'
       MSG_ERROR(msg)
     end if
     filnam_out = rmquotes(getpath)
     write(msg, '(5a)' )' mkfilename: get',trim(stringvar) ," from: ",trim(filnam_out), ch10
     call wrtout([std_out, ab_out], msg)
     ! Check whether file exists taking into account a possible NC file extension.
     if (xmpi_comm_rank(xmpi_world) == 0) then
       if (.not. file_exists(filnam_out) .and. .not. file_exists(strcat(filnam_out, ".nc"))) then
         MSG_ERROR(sjoin("Cannot find file:", filnam_out, "(with or without .nc extension)"))
       end if
     end if
     will_read = 1; return
   end if
 end if

 ! Treatment of the multi-dataset case (get is not relevant otherwise)
 if (ndtset /= 0) then

   if(ndtset==1 .and. get<0 .and. (jdtset_(1)+get>0)) then
     write(msg, '(7a,i0,a,i0,5a)' )&
     'You cannot use a negative value of get',trim(stringvar),' with only 1 dataset!',ch10, &
     'If you want to refer to a previously computed dataset,',ch10, &
     'you should give the absolute index of it (i.e. ', jdtset_(idtset)+get,' instead of ',get,').',ch10, &
     'Action: correct get',trim(stringvar),' in your input file.'
     MSG_ERROR(msg)
   end if

   if (idtset + get < 0) then
     write(msg, '(5a,i0,3a,i0,4a)' )&
     'The sum of idtset and get',trim(stringvar),' cannot be negative,',ch10,&
     'while they are idtset = ',idtset,', and get',trim(stringvar),' = ',get,ch10,&
     'Action: correct get',trim(stringvar),' in your input file.'
     MSG_ERROR(msg)
   end if

   if(get>0 .or. (get<0 .and. idtset+get>0) )then

     if(ird/=0 .and. get/=0)then
       write(msg, '(7a,i0,3a,i0,a,i0,7a)' )&
       'The input variables ird',trim(stringvar),' and get',trim(stringvar),' cannot be',ch10,&
       'simultaneously non-zero, while for idtset = ',idtset,',',ch10,&
       'they are ',ird,', and ',get,'.',ch10,&
       'Action: correct ird',trim(stringvar),' or get',trim(stringvar),' in your input file.'
       MSG_ERROR(msg)
     end if

     will_read=1

     ! Compute the dataset from which to take the file, and the corresponding index
     if(get<0 .and. idtset+get>0) jget=jdtset_(idtset+get)
     if(get>0) jget=get
     call appdig(jget,'',appen)

     ! Note use of output filename (filnam(4))
     filnam_out=trim(filnam(4))//'_DS'//trim(appen)//trim(stringfil)

     if(jdtset>=100)then
       write(msg, '(5a,i5,2a)' )&
       ' mkfilename : get',trim(stringvar) ,'/=0, take file ',trim(stringfil),' from output of DATASET ',jget,'.',ch10
     else
       write(msg, '(5a,i3,2a)' )&
       ' mkfilename : get',trim(stringvar) ,'/=0, take file ',trim(stringfil),' from output of DATASET ',jget,'.',ch10
     end if
     call wrtout([std_out, ab_out], msg)
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
!!      'old': file already exists
!!      'new': file does not exist; if file exists,
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
!!      anaddb,m_dtfil,m_effective_potential,m_polynomial_coeff,m_vcoul
!!      multibinit,ujdet
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine isfile(filnam, status)

!Arguments ------------------------------------
!scalars
 character(len=3),intent(in) :: status
 character(len=fnlen),intent(inout) :: filnam

!Local variables-------------------------------
!scalars
 logical :: ex,found
 integer :: ii,ios, ioserr
 character(len=500) :: msg
 character(len=fnlen) :: filnam_tmp, trialnam

! *************************************************************************

 filnam_tmp=filnam

 if (status=='old') then
   ! Check that old file exists
   inquire(file=filnam,iostat=ios,exist=ex)

   if (ios/=0) then
     write(msg,'(4a,i0,2a)')&
     'Checks for existence of file: ',trim(filnam),ch10,&
     'but INQUIRE statement returns error code',ios,ch10,&
     'Action: identify which problem appears with this file.'
     MSG_ERROR(msg)
   else if (.not.ex) then
     write(msg, '(5a)' )&
     'Checks for existence of file: ',trim(filnam),ch10,&
     'but INQUIRE finds file does not exist.',&
     'Action: check file name and re-run.'
     MSG_ERROR(msg)
   end if

 else if (status=='new') then

   ! Check that new output file does NOT exist
   ioserr = 0
   trialnam = filnam
   ii = 0
   inquire(file=trim(trialnam),iostat=ios,exist=ex)
   if (ios /= 0) then
     write(msg,'(3a)') &
     'Something is wrong with permissions for reading/writing on this filesystem.',ch10,&
     'Action: Check permissions.'
     MSG_ERROR(msg)
   end if

   if ( ex .eqv. .true. ) then
     write(msg,'(3a)')'Output file: ',trim(trialnam),' already exists.'
     MSG_COMMENT(msg)
     found=.false.

     ii=1
     do while ( (found .eqv. .false.) .and. (ii < 10000) )
       call int2char4(ii,msg)
       trialnam=trim(trim(filnam_tmp)//msg)
       inquire(file=trim(trialnam),iostat=ios,exist=ex)
       if ( (ex .eqv. .false.) .and. (ios == 0)) then
         found  = .true.
       end if
       if ( ios /= 0 )  ioserr=ioserr+1
       if ( ioserr > 10 ) then
         ! There is a problem => stop
         write(msg, '(2a,i0,2a)' )&
         'Check for permissions of reading/writing files on the filesystem', &
         '10 INQUIRE statements returned an error code like ',ios,ch10,&
         'Action: Check permissions'
         MSG_ERROR(msg)
       end if
       ii=ii+1
     end do
     if ( found .eqv. .true. ) then
       write(msg,'(4a)') 'Renaming old: ',trim(filnam),' to: ',trim(trialnam)
       MSG_COMMENT(msg)
       ioserr = clib_rename(filnam, trialnam)
       if ( ioserr /= 0 ) then
         write(msg,'(4a)') 'Failed to rename file: ', trim(filnam),' to: ',trim(trialnam)
         MSG_ERROR(msg)
       end if
     else
       write(msg,'(3a)')&
       'Have used all names of the form filenameXXXX, X in [0-9]',ch10,&
       'Action: clean up your directory and start over.'
       MSG_ERROR(msg)
     end if
   end if
   ! if ii > 0 we iterated so rename abi_out to abi_outXXXX and just write to abi_out
 else
   MSG_BUG(sjoin('Input status:', status, ' not recognized.'))
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
!!  input_path: String with input file path. Empty string activates files file legacy mode.
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
!!  (4) Root name for generic output files (wavefunctions, potential, density, DOS, hessian ...)
!!  (5) Root name for generic temporary files (wftmp1,wftmp2,kgunit,status ...)
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_log_status_state,int2char4,intagm,isfile,libpaw_log_flag_set
!!      parsefile,xmpi_barrier,xmpi_bcast
!!
!! SOURCE

subroutine iofn1(input_path, filnam, filstat, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 character(len=fnlen), intent(in) :: input_path
 character(len=fnlen), intent(out) :: filstat
 character(len=fnlen), intent(out) :: filnam(5)

!Local variables-------------------------------
 character(len=1) :: blank
 integer,parameter :: master = 0
 integer :: me, ios, nproc, ierr, ndtset, lenstr, marr, jdtset, tread, i1,i2
 logical :: ex
 character(len=fnlen) :: fillog, tmpfil, fname
 character(len=10) :: tag
 character(len=500) :: msg, errmsg
 character(len=strlen) :: string
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)
!*************************************************************************

 ! NOTE: In this routine it's very important to perform tests
 ! on possible IO errors (err=10, iomsg) because we are initializing the IO stuff
 ! It there's some problem with the hardware or some misconfiguration,
 ! it's very likely that the code will crash here and we should try to give useful error messages.

 blank = ' '; tmpfil = ''

 ! Determine who I am in comm
 me = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Define values of do_write_log and do_write_status parameters
 ! if a _NOLOG file exists no LOG file and no STATUS file are created for each cpu core
 ! if a _LOG file exists, a LOG file and a STATUS file are created for each cpu core
 ! if the #_of_cpu_core>NPROC_NO_EXTRA_LOG OR presence of ABI_MAIN_LOG_FILE, LOG file is only created for master proc
 ! if the #_of_cpu_core>NPROC_NO_EXTRA_STATUS OR presence of ABI_MAIN_LOG_FILE, STATUS file is only created for master proc
 inquire(file=ABI_NO_LOG_FILE, iostat=ios, exist=ex)
 if (ios /= 0) ex=.false.
 if (ex) then
   do_write_log=.false. ; do_write_status=.false.
   call abi_log_status_state(new_do_write_log=.false.,new_do_write_status=.false.)
   call libpaw_log_flag_set(.false.)
 else
   inquire(file=ABI_ENFORCE_LOG_FILE, iostat=ios, exist=ex)
   if (ios/=0) ex=.false.
   if (ex) then
     do_write_log=.true. ; do_write_status=.true.
     call abi_log_status_state(new_do_write_log=.true.,new_do_write_status=.true.)
     call libpaw_log_flag_set(.true.)
   else
     inquire(file=ABI_MAIN_LOG_FILE, iostat=ios, exist=ex)
     if (ios /= 0) ex=.false.
     if (ex .and. me /= master) then
       do_write_log=.false. ; do_write_status=.false.
       call abi_log_status_state(new_do_write_log=.false.,new_do_write_status=.false.)
       call libpaw_log_flag_set(.false.)
     else
       if (me /= master) then
         do_write_log= (nproc<NPROC_NO_EXTRA_LOG)
         call abi_log_status_state(new_do_write_log=(nproc<NPROC_NO_EXTRA_LOG))
         call libpaw_log_flag_set((nproc<NPROC_NO_EXTRA_LOG))
         do_write_status= (nproc<NPROC_NO_EXTRA_STATUS)
         call abi_log_status_state(new_do_write_status=(nproc<NPROC_NO_EXTRA_STATUS))
       end if
     end if
   end if
 end if
 !do_write_log = .True.

 if (me == master) then
   !  Eventually redefine standard input and standard output
   if (do_write_log) then
#if defined READ_FROM_FILE
     ! Take care of the output file
     tmpfil(1:fnlen)=blank
     tmpfil(1:3)='log'
     call isfile(tmpfil,'new')
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(tmpfil,msg,unit=std_out,form='formatted',status='new',action="write") /= 0) then
       MSG_ERROR(msg)
     end if
#endif
   else
     ! Redirect standard output to null
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(NULL_FILE, msg, unit=std_out, action="write") /= 0) then
       MSG_ERROR(msg)
     end if
   end if

#if defined READ_FROM_FILE
   ! Now take care of the "files" file
   tmpfil(1:fnlen)=blank
   tmpfil(1:9)='ab.files'
   write(msg, '(4a)' )&
    'Because of CPP option READ_FROM_FILE,',ch10,&
    'read file "ab.files" instead of standard input ' ,ch10
   MSG_COMMENT(msg)
   call isfile(tmpfil,'old')
   close(std_in, err=10, iomsg=errmsg)
   if (open_file(tmpfil,msg,unit=std_in,form='formatted',status='old',action="read") /= 0) then
     MSG_ERROR(msg)
   end if
#endif

   ! Print greetings for interactive user
   write(std_out,*,err=10,iomsg=errmsg)' ABINIT ',trim(abinit_version)
   write(std_out,*,err=10,iomsg=errmsg)' '

   if (len_trim(input_path) == 0) then
     ! Legacy Files file mode.
     write(std_out, "(2a)")" DeprecationWarning: ",ch10
     write(std_out, "(a)") "     The files file has been deprecated in Abinit9 and will be removed in Abinit10."
     write(std_out, "(2a)")"     Use the syntax `abinit t01.abi` where t01.abi is an input with pseudopotentials e.g.",ch10
     write(std_out, "(3a)")'            pseudos = "al.psp8, as.psp8"',ch10,ch10

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

   else
     ! Get prefix from input file. Default values are provided
     filnam(1) = input_path
     filnam(2) = trim(input_path)//".abo"
     filnam(3) = "i"
     filnam(4) = "o"
     filnam(5) = "t"

     fname = basename(input_path)
     i1 = index(fname, ".")
     !if (i1 /= 0) then
     if (i1 > 1) then
       ! file ext is present --> use prefix to initialize filnam
       i2 = index(input_path, ".", back=.True.)
       filnam(2) = input_path(:i2) // "abo"
       filnam(3) = fname(:i1-1) // "i"
       filnam(4) = fname(:i1-1) // "o"
       filnam(5) = fname(:i1-1) // "t"
     end if

     ! Read the file, stringify it and return the number of datasets.
     call parsefile(input_path, lenstr, ndtset, string, xmpi_comm_self)

     marr = max(1, ndtset)
     ABI_MALLOC(dprarr, (marr))
     ABI_MALLOC(intarr, (marr))
     jdtset = 0

     ! Allow user to change default values
     call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "output_file", tread, 'KEY', key_value=filnam(2))
     call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "indata_prefix", tread, 'KEY', key_value=filnam(3))
     call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "outdata_prefix", tread, 'KEY', key_value=filnam(4))
     call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "tmpdata_prefix", tread, 'KEY', key_value=filnam(5))

     ABI_FREE(dprarr)
     ABI_FREE(intarr)
   end if

   ! Check that old input file exists
   call isfile(filnam(1), 'old')
   ! Check that new output file does NOT exist
   call isfile(filnam(2), 'new')

   ! Check that root name for generic input and output differ
   if ( trim(filnam(3)) == trim(filnam(4)) ) then
     write(msg, '(3a)' )&
     'Root name for generic input and output files must differ ',ch10,&
     'Action: correct your "file" file.'
     MSG_ERROR(msg)
   end if

   ! Check that root names are at least 20 characters less than fnlen
   if ( len_trim(filnam(3)) >= (fnlen-20) ) then
     write(msg, '(a,a,a,a,a,i0,a,i0,a,a)' )&
     'Root name for generic input files is too long. ',ch10,&
     'It must be 20 characters less than the maximal allowed ',ch10,&
     'length of names, that is ',fnlen,', while it is: ',len_trim(filnam(3)),ch10,&
     'Action: correct your "file" file.'
     MSG_ERROR(msg)
   end if
   if ( len_trim(filnam(4)) >= (fnlen-20) ) then
     write(msg, '(a,a,a,a,a,i0,a,i0,a,a)' )&
     'Root name for generic output files is too long. ',ch10,&
     'It must be 20 characters less than the maximal allowed ',ch10,&
     'length of names, that is ',fnlen,', while it is: ',len_trim(filnam(4)),ch10,&
     'Action: correct your "file" file.'
     MSG_ERROR(msg)
   end if
   if ( len_trim(filnam(5)) >= (fnlen-20) ) then
     write(msg, '(a,a,a,a,a,i0,a,i0,a,a)' )&
     'Root name for generic temporary files is too long. ',ch10,&
     'It must be 20 characters less than the maximal allowed ',ch10,&
     'length of names, that is ',fnlen,', while it is: ',len_trim(filnam(5)),ch10,&
     'Action: correct your "file" file.'
     MSG_ERROR(msg)
   end if

 end if ! master only

 ! Communicate filenames to all processors
 call xmpi_bcast(filnam,master,comm,ierr)

 ! Create a name for the status file, based on filnam(5)
 filstat=trim(filnam(5))//'_STATUS'

 ! Redefine the log unit if not the master
 if (me /= master) then
   call int2char4(me,tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   filstat=trim(filstat)//'_P-'//trim(tag)
   if (do_write_log) then
     fillog=trim(filnam(5))//'_LOG_'//trim(tag)
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(fillog,msg,unit=std_out,status='unknown',action="write") /= 0) then
       MSG_ERROR(msg)
     end if
   else
     close(std_out, err=10, iomsg=errmsg)
     if (open_file(NULL_FILE,msg,unit=std_out,action="write") /= 0) then
       MSG_ERROR(msg)
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

end module m_dtfil
!!***
