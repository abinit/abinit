!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_chkinp
!! NAME
!!  m_chkinp
!!
!! FUNCTION
!! Check consistency of Abinit input data against itself.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, MKV, DRH, MVer)
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

module m_chkinp

 use defs_basis
 use m_gwdefs
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xomp
 use libxc_functionals
 use m_dtset

 use defs_datatypes,   only : pspheader_type
 use defs_abitypes,    only : MPI_type
 use m_numeric_tools,  only : iseven
 use m_symtk,          only : chkgrp, chkorthsy
 use m_geometry,       only : metric
 use m_fftcore,        only : fftalg_has_mpi
 use m_exit,           only : get_timelimit
 use m_parser,         only : chkdpr, chkint, chkint_eq, chkint_ge, chkint_le, chkint_ne

 implicit none

 private
!!***

 public :: chkinp
!!***

contains
!!***

!!****f* ABINIT/chkinp
!! NAME
!! chkinp
!!
!! FUNCTION
!! Check consistency of input data against itself.
!! Please: use the alphabetic order
!! Please: use the routines chkint_eq, chkint_ne, chkint_ge, chkint_le, and chkdpr
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables
!!  iout=unit number for output file
!!  mpi_enregs(0:ndtset_alloc)=informations about MPI parallelization
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set.
!!  npsp=number of pseudopotentials
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!  comm: MPI communicator (MPI_COMM_WORLD)
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdpr,chkgrp,chkint,chkint_eq,chkint_ge,chkint_le,chkint_ne,chkorthsy
!!      dtset_copy,dtset_free,metric,wrtout,xmpi_sum
!!
!! SOURCE

subroutine chkinp(dtsets,iout,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,ndtset,ndtset_alloc,npsp, comm
 type(MPI_type),intent(in) :: mpi_enregs(0:ndtset_alloc)
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 logical :: twvl,allow
 logical :: wvlbigdft=.false.
 integer :: ttoldfe,ttoldff,ttolrff,ttolvrs,ttolwfr
 integer :: bantot,ia,iatom,ib,iband,idtset,ierr,iexit,ii,iimage,ikpt,ilang,intimage,ierrgrp
 integer :: ipsp,isppol,isym,itypat,iz,jdtset,jj,kk,maxiatsph,maxidyn,minplowan_iatom,maxplowan_iatom
 integer :: mband,mgga,miniatsph,minidyn,mod10,mpierr,all_nprocs
 integer :: mu,natom,nfft,nfftdg,nkpt,nloc_mem,nlpawu,nproc,nspden,nspinor,nsppol,optdriver,response,usepaw,usewvl
 integer :: fftalg !,fftalga,fftalgc,
 real(dp) :: delta,dz,sumalch,summix,sumocc,ucvol,wvl_hgrid,zatom
 character(len=1000) :: msg
 type(dataset_type) :: dt
!arrays
 integer :: cond_values(4),nprojmax(0:3)
 integer :: gpu_devices(5)=(/-2,-2,-2,-2,-2/)
 integer,allocatable :: ierr_dtset(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: frac(:,:)
 character(len=32) :: cond_string(4)
 character(len=32) :: input_name

! *************************************************************************

 DBG_ENTER("COLL")

 all_nprocs = xmpi_comm_size(comm)

!Some initialisations
 cond_string(1:4)='#####'
 cond_values(1:4)=(/0,0,0,0/)
 ABI_ALLOCATE(ierr_dtset,(ndtset_alloc))
 ierr_dtset=0

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc
   if(mpi_enregs(idtset)%me<0) cycle
   jdtset=dtsets(idtset)%jdtset
   if(ndtset==0)jdtset=0
   ierr=0

   if(jdtset/=0)then
     write(msg, '(a,a,a,i4,a)' ) ch10,' chkinp: Checking input parameters for consistency,',' jdtset=',jdtset,'.'
   else
     write(msg, '(a,a)' ) ch10,' chkinp: Checking input parameters for consistency.'
   end if
   call wrtout(iout,msg,'COLL')
   call wrtout(std_out,msg,'COLL')

   ! Will test directly on the dataset "dt"
   dt = dtsets(idtset)%copy()

   ! Copy or initialize locally a few input dataset values
   fftalg   =dt%ngfft(7)
   !fftalga  =fftalg/100; fftalgc=mod(fftalg,10)
   natom    =dt%natom
   nkpt     =dt%nkpt
   nspden   =dt%nspden
   nspinor  =dt%nspinor
   nsppol   =dt%nsppol
   optdriver=dt%optdriver
   usepaw   =dt%usepaw
   usewvl   =dt%usewvl
   intimage=1 ; if(dtsets(idtset)%nimage>2)intimage=2
   rprimd(:,:)=dtsets(idtset)%rprimd_orig(:,:,intimage)    ! For the purpose of checking symmetries
   response=0
   if(dt%rfelfd/=0.or.dt%rfphon/=0.or.dt%rfstrs/=0.or.dt%rfddk/=0.or.dt%rfuser/=0 &
&   .or.dt%rf2_dkdk/=0.or.dt%rf2_dkde/=0.or.dt%rfmagn/=0.or.dt%d3e_pert1_elfd/=0 &
&   .or.dt%d3e_pert2_elfd/=0.or.dt%d3e_pert3_elfd/=0.or.dt%d3e_pert1_phon/=0 &
&   .or.dt%d3e_pert2_phon/=0.or.dt%d3e_pert3_phon/=0) response=1
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
   nproc=mpi_enregs(idtset)%nproc
   mgga=0;if(dt%ixc>=31.and.dt%ixc<=34)mgga=1
   if (dt%ixc<0.and.libxc_functionals_ismgga()) mgga=1

!  =====================================================================================================
!  Check the values of variables, using alphabetical order
!  PLEASE: use the routines chkint_eq, chkint_ne, chkint_ge, chkint_le, chkdpr

    !  iomode Must be one of 0, 1, 3
   call chkint_eq(0,0,cond_string,cond_values,ierr,'iomode',dt%iomode,3, [IO_MODE_FORTRAN,IO_MODE_MPI,IO_MODE_ETSF],iout)
!  However, if mpi_io is not enabled, must be one of 0, 3.
   if(xmpi_mpiio==0)then
     cond_string(1)='enable_mpi_io' ;  cond_values(1)=0
!    Make sure that iomode is 0 or 3
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iomode',dt%iomode,2,[IO_MODE_FORTRAN,IO_MODE_ETSF],iout)
   end if
   if (dt%iomode == IO_MODE_NETCDF .and. dt%npspinor == 2) then
     MSG_ERROR_NOSTOP("npspinor == 2 not compatible with netcdf", ierr)
   end if

!  accuracy
   call chkint_eq(0,0,cond_string,cond_values,ierr,'accuracy',dt%accuracy,7,(/0,1,2,3,4,5,6/),iout)

!  adpimd
   call chkint_eq(0,0,cond_string,cond_values,ierr,'accuracy',dt%adpimd,2,(/0,1/),iout)

!  amu Check that atomic masses are > 0 if ionmov = 1
   do iimage=1,dt%nimage
     if (dt%ionmov==1) then
       do itypat=1,dt%ntypat
         cond_string(1)='ionmov' ; cond_values(1)=1
         write(input_name,'(a4,i2,a1,i2,a1)')'amu(',itypat,',',iimage,')'
         call chkdpr(1,1,cond_string,cond_values,ierr,input_name,dt%amu_orig(itypat,iimage),1,tol8,iout)
       end do
     end if
   end do

!  autoparal
   call chkint_eq(0,0,cond_string,cond_values,ierr,'autoparal',dt%autoparal,5,(/0,1,2,3,4/),iout)

!  auxc_scal
   call chkdpr(0,0,cond_string,cond_values,ierr,'auxc_scal',dt%auxc_scal,1,0.0_dp,iout)

!  bdberry
   if(dt%berryopt>0.and.dt%nberry>0.and.&
&   dt%berryopt/= 4.and.dt%berryopt/= 6.and.dt%berryopt/= 7.and.&
&   dt%berryopt/=14.and.dt%berryopt/=16.and.dt%berryopt/=17) then
     do ii=1,2*nsppol
       cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
       cond_string(2)='nberry'   ; cond_values(2)=dt%nberry
       write(input_name,'(a4,i1,a1)')'bdberry(',ii,')'
       call chkint_ge(2,2,cond_string,cond_values,ierr,input_name,dt%bdberry(ii),1,iout)
     end do
!    bdberry(2) must be greater than bdberry(1)
     cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
     cond_string(2)='nberry'   ; cond_values(2)=dt%nberry
     cond_string(3)='bdberry(1)'   ; cond_values(3)=dt%bdberry(1)
     call chkint_ge(3,3,cond_string,cond_values,ierr,'bdberry(2)',dt%bdberry(2),dt%bdberry(1),iout)
     if(nsppol==2)then
!      bdberry(4) must be greater than bdberry(3)
       cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
       cond_string(2)='nberry'   ; cond_values(2)=dt%nberry
       cond_string(3)='bdberry(3)'   ; cond_values(3)=dt%bdberry(3)
       call chkint_ge(3,3,cond_string,cond_values,ierr,'bdberry(4)',dt%bdberry(4),dt%bdberry(3),iout)
     end if
!    Make sure all nband(nkpt) are >= bdberry
     do isppol=1,nsppol
       do ikpt=1,nkpt
         if (dt%nband(ikpt+(isppol-1)*nkpt)<=dt%bdberry(2*isppol)) then
           cond_string(1)='ikpt'   ; cond_values(1)=ikpt
           cond_string(2)='isppol' ; cond_values(2)=isppol
           cond_string(3)='nband'  ; cond_values(3)=dt%nband(ikpt+(isppol-1)*nkpt)
           call chkint_le(0,3,cond_string,cond_values,ierr,&
&           'bdberry',dt%bdberry(2*isppol),dt%nband(ikpt+(isppol-1)*nkpt),iout)
           if(ierr==1)exit
         end if
       end do
     end do
   end if

!  berryopt
!  berryopt must be between -3 to +4, 6,7,14,16,17
   call chkint_eq(0,0,cond_string,cond_values,ierr,&
&   'berryopt',dt%berryopt,13,(/-3,-2,-1,0,1,2,3,4,6,7,14,16,17/),iout)
!  berryopt must be positive when mkmem==0
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_ge(1,1,cond_string,cond_values,ierr,'berryopt',dt%berryopt,0,iout)
   end if
!  berryopt must be positive when occopt/=1
   if(dt%occopt/=1)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_ge(1,1,cond_string,cond_values,ierr,'berryopt',dt%berryopt,0,iout)
   end if
!  berryopt cannot be 4,6,7,14,16,17 when toldfe, tolvrs, toldff and tolrff are zero or negative
   if ((dt%toldfe < tiny(one)).and.(dt%tolvrs < tiny(one)).and.&
&   (dt%toldff < tiny(one)).and.(dt%tolrff < tiny(one))) then
     cond_string(1)='toldfe' ; cond_values(1)=dt%toldfe
     cond_string(2)='toldff' ; cond_values(2)=dt%toldff
     cond_string(3)='tolrff' ; cond_values(3)=dt%tolrff
     cond_string(4)='tolvrs' ; cond_values(4)=dt%tolvrs
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/4/),iout)
     cond_string(1)='toldfe' ; cond_string(2)='toldff' ; cond_string(3)='tolrff' ; cond_string(4)='tolvrs'
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/6/),iout)
     cond_string(1)='toldfe' ; cond_string(2)='toldff' ; cond_string(3)='tolrff' ; cond_string(4)='tolvrs'
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/7/),iout)
     cond_string(1)='toldfe' ; cond_string(2)='toldff' ; cond_string(3)='tolrff' ; cond_string(4)='tolvrs'
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/14/),iout)
     cond_string(1)='toldfe' ; cond_string(2)='toldff' ; cond_string(3)='tolrff' ; cond_string(4)='tolvrs'
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/16/),iout)
     cond_string(1)='toldfe' ; cond_string(2)='toldff' ; cond_string(3)='tolrff' ; cond_string(4)='tolvrs'
     call chkint_ne(4,4,cond_string,cond_values,ierr,'berryopt',dt%berryopt,1,(/17/),iout)
   end if
!  Non-zero berryopt and usepaw==1 cannot be done unless response==0
   if (usepaw==1.and.dt%berryopt/=0) then
     cond_string(1)='usepaw' ; cond_values(1)=1
     cond_string(2)='berryopt' ; cond_values(2)=dt%berryopt
     call chkint_eq(1,2,cond_string,cond_values,ierr,'response',response,1,(/0/),iout)
   end if
!  Non-zero berryopt and usepaw==1 and kptopt/=3 cannot be done unless symmorphi=0
!  (that is, nonsymmorphic symmetries do not work yet
!  Update MT 2017-05-31: nonsymmorphic symmetries seem also to be an issue for NCPP
   if (usepaw==1.and.dt%berryopt/=0.and.dt%kptopt/=3) then
  !if (dt%berryopt/=0.and.dt%kptopt/=3) then
     cond_string(1)='usepaw'; cond_values(1)=1
     cond_string(2)='berryopt'; cond_values(2)=dt%berryopt
     cond_string(3)='kptopt'; cond_values(3)=dt%kptopt
     call chkint_eq(1,3,cond_string,cond_values,ierr,'symmorphi',dt%symmorphi,1,(/0/),iout)
   end if
!  Restrictions for berryopt=4,6,7,14,16,17
   if (usepaw==1.and.&
      (dt%berryopt== 4.or.dt%berryopt== 6.or.dt%berryopt== 7.or.&
       dt%berryopt==14.or.dt%berryopt==16.or.dt%berryopt==17)) then
!     if (nsppol==2.and.nproc>1) then
!       write(msg,'(3a)') &
!&       'For berryopt = 4,6,7,14,16,17 and nsppol=2, nproc must = 1 ',ch10,&
!&       'Action: change number of processes'
!       MSG_ERROR_NOSTOP(msg, ierr)
!     end if
   end if
!  Non-zero berryopt and usepaw==1 and kpt // requires nproc to be a divisor of nkpt
   if (usepaw==1.and.dt%berryopt/=0.and.nproc>1.and.mod(dt%nkpt,nproc)/=0) then
     write(msg, '(3a)' )&
      'For berryopt /= 0 with PAW in parallel, nproc must be a divisor of nkpt ',ch10,&
      'Action: change number of processes or kpts such that nproc divides nkpt evenly'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  Finite electric/displacement field checks
   if (dt%berryopt==4) then
     if (maxval(abs(dt%dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efieldbar(1:3)))>tiny(0.0_dp)) then
       write(msg,'(5a)' ) &
        'When berryopt==4, only efield is needed, other input field',ch10,&
        '(dfield,red_dfield,red_efield,red_efieldbar) should be zero.',ch10,&
        'Action: delete unneeded field in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%berryopt==14) then
     if (maxval(abs(dt%dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%efield(1:3)))>tiny(0.0_dp)) then
       write(msg,'(5a)') &
        'When berryopt==14, only red_efieldbar is needed, other input field',ch10,&
        '(dfield,red_dfield,efield,red_efield) should be zero.',ch10,&
        'Action: delete unneeded field in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%berryopt==6) then
     if (maxval(abs(dt%red_dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efieldbar(1:3)))>tiny(0.0_dp)) then
       write(msg,'(5a)') &
         'When berryopt==6, only dfield and efield are needed, other input field',ch10,&
         '(red_dfield,red_efield,red_efieldbar) should be zero.',ch10,&
         'Action: delete unneeded field in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%berryopt==16) then
     if (maxval(abs(dt%dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%efield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efieldbar(1:3)))>tiny(0.0_dp)) then
       write(msg,'(5a)')  &
         'When berryopt==16, only red_dfield and red_efield are needed, other input field',ch10,&
         '(dfield,efield,red_efieldbar) should be zero.',ch10,&
         'Action: delete unneeded field in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%berryopt==17) then
     if (maxval(abs(dt%dfield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%efield(1:3)))>tiny(0.0_dp).or.&
         maxval(abs(dt%red_efield(1:3)))>tiny(0.0_dp)) then
       write(msg,'(5a)') &
        'When berryopt==17, only red_dfield and red_efieldbar are needed, other input field',ch10,&
        '(dfield,efield,red_efield) should be zero.',ch10,&
        'Action: delete unneeded field in input file.'
        MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if ((dt%jfielddir(1)/=1.and.dt%jfielddir(1)/=2).or.&
         (dt%jfielddir(2)/=1.and.dt%jfielddir(2)/=2).or.&
         (dt%jfielddir(3)/=1 .and.dt%jfielddir(3)/=2)) then
       write(msg,'(5a)') &
        'When berryopt==17, jfielddir can only be 1 or 2 to controls whether reduced electric field',ch10,&
        '(jfielddir=1) or reduced electric displacement field (jfielddir=2) is chosen to be fixed', ch10,&
        'Action: change jfielddir to be 1 or 2 in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  berrystep
   call chkint_ge(0,0,cond_string,cond_values,ierr,'berrystep',dt%berrystep,1,iout)
   if(nproc>1)then
     cond_string(1)='nproc'; cond_values(1)=mpi_enregs(idtset)%nproc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'berrystep',dt%berrystep,1,(/1/),iout)
   end if

!  boxcutmin
   call chkdpr(0,0,cond_string,cond_values,ierr,'dilatmx',dt%dilatmx,1,zero,iout)
!  if(response==1)then
!    cond_string(1)='response' ; cond_values(1)=1
!    call chkdpr(1,1,cond_string,cond_values,ierr,'boxcutmin',dt%boxcutmin,1,two,iout)
!  end if

!  builtintest
   call chkint_eq(0,0,cond_string,cond_values,ierr,'builtintest',dt%builtintest,8,(/0,1,2,3,4,5,6,7/),iout)

!  chkdilatmx
   call chkint_eq(0,0,cond_string,cond_values,ierr,'chkdilatmx',dt%chkdilatmx,2,(/0,1/),iout)

!  chksymbreak
   call chkint_eq(0,0,cond_string,cond_values,ierr,'chksymbreak',dt%chksymbreak,3,(/0,1,-1/),iout)
   if(dt%chksymbreak==1)then
!    Check the values of tnons
     do isym=1,dt%nsym
       do ii=1,3
         delta=dt%tnons(ii,isym)*eight
         if(abs(delta-nint(delta))>tol6)then
           delta=dt%tnons(ii,isym)*three*four
           if(abs(delta-nint(delta))>tol6)then
             write(msg, '(8a,i4,2a,9i3,2a,3es16.6,4a)' ) ch10,&
&             ' chkinp: WARNING -',ch10,&
&             '   Chksymbreak=1 . Found potentially symmetry-breaking value of tnons, ', ch10,&
&             '   which is neither a rational fraction in 1/8th nor in 1/12th :', ch10,&
&             '   for the symmetry number ',isym,ch10,&
&             '   symrel is ',dt%symrel(1:3,1:3,isym),ch10,&
&             '   tnons is ',dt%tnons(1:3,isym),ch10,&
&             '   Please, read the description of the input variable chksymbreak,',ch10,&
&             '   then, if you feel confident, you might switch it to zero, or consult with the forum.'
             call wrtout(iout,msg,'COLL')
             call wrtout(std_out,msg,'COLL')
             !ierr=ierr+1 ! moved this to a warning: for slab geometries arbitrary tnons can appear along the vacuum direction
           end if
         end if
       end do
     end do
   end if

!  constraint_kind
   do itypat=1,dt%ntypat
     cond_string(1)='itypat';cond_values(1)=itypat
     call chkint_eq(0,1,cond_string,cond_values,ierr,'constraint_kind',dt%constraint_kind(itypat),8,(/0,1,2,3,10,11,12,13/),iout)
     !Only potential self-consistency is currently allowed with constrained_dft
     if (dt%iscf>10) then
       cond_string(2)='iscf';cond_values(2)=dt%iscf
       call chkint_eq(2,2,cond_string,cond_values,ierr,'constraint_kind',dt%constraint_kind(itypat),1,(/0/),iout)
     endif
     if (dt%ionmov==4) then
       cond_string(2)='ionmov';cond_values(2)=dt%ionmov
       call chkint_eq(2,2,cond_string,cond_values,ierr,'constraint_kind',dt%constraint_kind(itypat),1,(/0/),iout)
     endif

   enddo

!  densfor_pred
   if(dt%iscf>0)then
     cond_string(1)='iscf';cond_values(1)=dt%iscf
     call chkint_le(0,1,cond_string,cond_values,ierr,'densfor_pred',dt%densfor_pred,6,iout)
     call chkint_ge(0,1,cond_string,cond_values,ierr,'densfor_pred',dt%densfor_pred,-6,iout)
     if (dt%densfor_pred<0.and.mod(dt%iprcel,100)>=61.and.(dt%iprcel<71.or.dt%iprcel>79)) then
       cond_string(1)='iscf';cond_values(1)=dt%iscf
       cond_string(2)='iprcel';cond_values(2)=dt%iprcel
       call chkint_ge(1,2,cond_string,cond_values,ierr,'densfor_pred',dt%densfor_pred,0,iout)
     end if
   end if

!  diecut
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
     cond_string(2)='4*ecut' ; cond_values(1)=4*dt%ecut
!    Checks that presently diecut is 4*ecut
     call chkdpr(1,1,cond_string,cond_values,ierr,'diecut',dt%diecut,0,4*dt%ecut,iout)
   end if

!  diemac
   call chkdpr(0,0,cond_string,cond_values,ierr,'diemac',dt%diemac,1,0.01_dp,iout)

!  dilatmx
   call chkdpr(0,0,cond_string,cond_values,ierr,'dilatmx',dt%dilatmx,1,zero,iout)
   if(dt%chkdilatmx==1)then
     cond_string(1)='chkdilatmx' ; cond_values(1)=1
!    Checks that presently chkdilatmx is smaller than 1.15
     call chkdpr(1,1,cond_string,cond_values,ierr,'dilatmx',dt%dilatmx,-1,1.15_dp,iout)
   end if

!  dmatpuopt
   if (dt%usepawu>0) then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmatpuopt',dt%dmatpuopt,10,(/1,2,3,4,5,6,7,8,9,10/),iout)
   end if

!  dmatudiag
   if (dt%usepawu>0) then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmatudiag',dt%dmatudiag,3,(/0,1,2/),iout)
   end if


!  dmftbandi, dmftbandf
   if (dt%usedmft>0) then
     call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftcheck',dt%dmftcheck,4,(/-1,0,1,2/),iout)
     if(dt%dmftcheck/=-1) then
       cond_string(1)='usedmft' ; cond_values(1)=dt%usedmft
       call chkint_eq(1,1,cond_string,cond_values,ierr,'occopt',dt%occopt,1,(/3/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftbandi',dt%dmftbandi,1,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftbandf',dt%dmftbandf,dt%dmftbandi,iout)
       cond_string(1)='mband' ; cond_values(1)=dt%mband
       call chkint_le(0,1,cond_string,cond_values,ierr,'dmftbandi',dt%dmftbandi,dt%mband,iout)
       cond_string(1)='mband' ; cond_values(1)=dt%mband
       call chkint_le(0,1,cond_string,cond_values,ierr,'dmftbandf',dt%dmftbandf,dt%mband,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_entropy',dt%dmft_entropy,0,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_iter',dt%dmft_iter,0,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_kspectralfunc',dt%dmft_kspectralfunc,2,(/0,1/),iout)
       if(dt%dmft_kspectralfunc==1) then
         cond_string(1)='dmft_kspectralfunc' ; cond_values(1)=1
         call chkint_eq(0,1,cond_string,cond_values,ierr,'iscf',dt%iscf,2,(/-2,-3/),iout)
       endif
       if((dt%dmft_solv<6.or.dt%dmft_solv>7).and.dt%ucrpa==0) then
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nwlo',dt%dmft_nwlo,1,iout)
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nwli',dt%dmft_nwli,1,iout)
       end if
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_read_occnd',dt%dmft_read_occnd,3,(/0,1,2/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_occnd_imag',dt%dmft_occnd_imag,2,(/0,1/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_rslf',dt%dmft_rslf,3,(/-1,0,1/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_mxsf',dt%dmft_mxsf,1,zero,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_mxsf',dt%dmft_mxsf,-1,one,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_solv',dt%dmft_solv,9,(/-2,-1,0,1,2,5,6,7,8/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_tolfreq',dt%dmft_tolfreq,-1,0.01_dp,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_tollc',dt%dmft_tollc,-1,tol5,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_charge_prec',dt%dmft_charge_prec,-1,tol4,iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkdpr(0,1,cond_string,cond_values,ierr,'dmft_charge_prec',dt%dmft_charge_prec,1,tol20,iout)
       if(dt%usepawu==14) then
         cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_dc',dt%dmft_dc,1,(/5/),iout)
       endif
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_dc',dt%dmft_dc,4,(/0,1,2,5/),iout)
       if(dt%getwfk==0.and.dt%irdwfk==0.and.dt%irdden==0.and.dt%getden==0.and.dt%ucrpa==0) then
         write(msg,'(3a,i3,a,i3,a,i3,a,i3,a)' )&
&         'When usedmft==1, A WFC file or a DEN file have to be read. In the current calculation:',ch10, &
&         '  getwfk =',dt%getwfk, &
&         '  irdwfk =',dt%irdwfk, &
&         '  getden =',dt%getden, &
&         '  irdden =',dt%irdden, &
&         '  Action: use a restart density or wfc file'
         if(dt%iscf>0) MSG_ERROR(msg)
       end if
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_t2g',dt%dmft_t2g,2,(/0,1/),iout)
       cond_string(1)='usedmft' ; cond_values(1)=1
       call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_x2my2d',dt%dmft_x2my2d,2,(/0,1/),iout)
       if (dt%dmft_solv>=4.and.dt%ucrpa==0) then
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftqmc_l',dt%dmftqmc_l,1,iout)
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkdpr(0,1,cond_string,cond_values,ierr,'dmftqmc_n',dt%dmftqmc_n,1,one,iout)
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftqmc_seed',dt%dmftqmc_seed,0,iout)
         cond_string(1)='usedmft' ; cond_values(1)=1
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftqmc_therm',dt%dmftqmc_therm,1,iout)
       end if
       if (dt%dmft_solv>=5) then
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftctqmc_basis',dt%dmftctqmc_basis,3,(/0,1,2/),iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftctqmc_check',dt%dmftctqmc_check,4,(/0,1,2,3/),iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftctqmc_gmove',dt%dmftctqmc_gmove,0,iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftctqmc_meas',dt%dmftctqmc_meas,1,iout)
#if defined HAVE_TRIQS_v2_0 || defined HAVE_TRIQS_v1_4
         if (dt%dmft_solv==6.or.dt%dmft_solv==7) then
           cond_string(1)='dmft_solv' ; cond_values(1)=5
           call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftqmc_l',dt%dmftqmc_l,2*dt%dmft_nwli+1,iout)
           cond_string(1)='usedmft' ; cond_values(1)=1
           call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftctqmc_triqs_nleg',dt%dmftctqmc_triqs_nleg,1,iout)
         end if
#endif
       end if
       if (dt%dmft_solv==5) then
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftctqmc_correl',dt%dmftctqmc_correl,2,(/0,1/),iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftctqmc_grnns',dt%dmftctqmc_grnns,2,(/0,1/),iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftctqmc_mrka',dt%dmftctqmc_mrka,0,iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmftctqmc_mov',dt%dmftctqmc_mov,2,(/0,1/),iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmftctqmc_order',dt%dmftctqmc_order,0,iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nwlo',dt%dmft_nwlo,2*dt%dmftqmc_l,iout)
       end if
       if (dt%dmft_entropy>=1) then
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_ge(0,1,cond_string,cond_values,ierr,'dmft_nlambda',dt%dmft_nlambda,3,iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_le(0,1,cond_string,cond_values,ierr,'dmft_entropy',dt%dmft_entropy,dt%dmft_nlambda,iout)
         cond_string(1)='dmft_solv' ; cond_values(1)=5
         call chkint_eq(0,1,cond_string,cond_values,ierr,'dmft_dc',dt%dmft_dc,1,(/1/),iout)
         if (dt%dmft_solv /= 5 ) then
           write(msg,'(3a,i3,a,i3,a,i3,a,i3,a)' )&
&           'When dmft_entropy>=1, the impurity solver has to be currently  dmft_solv=5:',ch10, &
&           'Action: change your dmft_solv input'
           MSG_ERROR(msg)
         end if
       end if
     end if
   end if

#if !defined HAVE_TRIQS_v2_0 && !defined HAVE_TRIQS_v1_4
   if(dt%dmft_solv>=6.and.dt%dmft_solv<=7) then
     write(msg, '(a,a,a)' )&
&     ' dmft_solv=6, or 7 is only relevant if the TRIQS library is linked',ch10,&
&     ' Action: check compilation options'
     MSG_ERROR(msg)
   end if
#endif

!  dosdeltae
   call chkdpr(0,0,cond_string,cond_values,ierr,'dosdeltae',dt%dosdeltae,1,0.0_dp,iout)

!  dynimage between 0 and 1
   maxidyn=maxval(dt%dynimage(:))
   minidyn=minval(dt%dynimage(:))
   call chkint_ge(0,0,cond_string,cond_values,ierr,'dynimage',minidyn,0,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'dynimage',maxidyn,1,iout)

!  ecut
!  With planewaves, one must use positive ecut
   if(usewvl==0)then
     if (abs(dt%ecut+1._dp)<tol8) then
       write(msg, '(3a)' )&
        'The input keyword "ecut" is compulsory !',ch10,&
        'Action: add a value for "ecut" in the input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     else
       cond_string(1)='usewvl' ; cond_values(1)=usewvl
       call chkdpr(1,1,cond_string,cond_values,ierr,'ecut',dt%ecut,1,tol8,iout)
     end if
   end if

!  pawecutdg (placed here to stop before ngfftdg)
   if (usepaw==1) then
     if(usewvl==0) then
       call chkdpr(1,0,cond_string,cond_values,ierr,'pawecutdg',dt%pawecutdg,1,tol8,iout)
       cond_string(1)='ecut' ; cond_values(1)=dt%ecut
       call chkdpr(1,1,cond_string,cond_values,ierr,'pawecutdg',dt%pawecutdg,1,dt%ecut,iout)
     else
       if(dt%pawecutdg > 0.d0) then
         MSG_ERROR('In PAW+WVL do not use pawecutdg')
       end if
     end if
   end if

!  ecuteps
   if( ANY(optdriver == [RUNL_SCREENING]) )then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecuteps',dt%ecuteps,1,0.0_dp,iout)
     if (dt%ecuteps <= 0) then
       MSG_ERROR_NOSTOP("ecuteps must be > 0 if optdriver == 3", ierr)
     end if
     if(dt%fftgw<20 .and. dt%fftgw/=0)then
       if(dt%ecutwfn<dt%ecuteps-tol8)then
         write(msg,'(a,es16.6,a,es16.6,a,6a)')&
          'The values of ecutwfn and ecuteps are ', dt%ecutwfn,' and ',dt%ecuteps,ch10,&
          'With fftgw lower than 20, one expect ecuteps to be smaller or equal to ecutwfn.',ch10,&
          'Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
          'Action: use another value of fftgw (e.g. 21), or adjust ecutwfn with ecuteps.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end if
   end if

!  ecutsigx
!  @MG FIXME reinstate this check, after having rewritten FFT treatment in GW
   if( ANY( optdriver==(/RUNL_SIGMA/) ) .and..FALSE.)then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsigx',dt%ecutsigx,1,0.0_dp,iout)
     if(dt%fftgw<20)then
       if(dt%ecutwfn<dt%ecutsigx-tol8)then
         write(msg,'(a,es16.6,a,es16.6,a,6a)')&
          'The values of ecutwfn and ecutsigx are ', dt%ecutwfn,' and ',dt%ecutsigx,ch10,&
          'With fftgw lower than 20, one expect ecutsigx to be smaller or equal to ecutwfn.',ch10,&
          'Indeed, one is wasting memory without gaining CPU time or accuracy.',ch10,&
          'Action: use another value of fftgw (e.g. 21), or adjust ecutwfn with ecutsigx.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end if
   end if

   if ( optdriver==RUNL_BSE) then
     ! Check for BSE calculations that are not implemented.
     if (dt%nspinor == 2) then
       MSG_ERROR_NOSTOP("BSE with nspinor 2 not implemented", ierr)
     end if
   end if

   ! Check for GW calculations that are not implemented.
   if (ANY(optdriver == [RUNL_SCREENING, RUNL_SIGMA])) then
     if (dt%nspinor == 2) then
       if (dt%usepaw == 1) then
         MSG_ERROR_NOSTOP("GW with PAW and nspinor 2 not implemented", ierr)
       end if
       !if (optdriver == RUNL_SCREENING .and. dt%symchi == 1) then
       !  MSG_ERROR_NOSTOP("Screening with symchi 1 and nspinor 2 not implemented", ierr)
       !end if
       !if (optdriver == RUNL_SIGMA .and. dt%symsigma == 1) then
       !  MSG_ERROR_NOSTOP("Self-energy with symsigma 1 and nspinor 2 not implemented", ierr)
       !end if
       if (optdriver == RUNL_SIGMA .and. &
       any(mod(dt%gwcalctyp, 10) == [SIG_GW_AC, SIG_QPGW_PPM, SIG_QPGW_CD])) then
         MSG_ERROR_NOSTOP("analytic-continuation, model GW with nspinor 2 are not implemented", ierr)
       end if
       !if (optdriver == RUNL_SIGMA .and. mod(dt%gwcalctyp, 100) >= 10) then
       !  MSG_ERROR_NOSTOP("Self-consistent GW with nspinor == 2 not implemented", ierr)
       !end if
       if (dt%gwcomp /= 0) then
         MSG_ERROR_NOSTOP("gwcomp /= 0 with nspinor 2 not implemented", ierr)
       end if
     end if ! nspinor 2

     if (maxval(abs(dt%istwfk(1:nkpt))) > 1 .and. mod(dt%gwcalctyp, 100) >= 20) then
       write(msg, "(3a)")"Self-consistent GW with istwfk > 1 not supported.",ch10, &
       "Please regenerate your WFK file with istwfk *1"
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

     ! Avoid wasting CPUs if nsppol==2.
     if (dt%nsppol==2 .and. .not. iseven(nproc) .and. nproc > 1) then
       write(msg,'(3a)') "Spin-polarized GW calculations should be run with an even number of processors ",ch10,&
        " for achieving an optimal distribution of memory and CPU load. Change the number of processors."
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  ecutsm
   call chkdpr(0,0,cond_string,cond_values,ierr,'ecutsm',dt%ecutsm,1,0.0_dp,iout)
!  With non-zero optcell, one must use non-zero ecutsm
   if(dt%optcell/=0 )then
     cond_string(1)='optcell' ; cond_values(1)=dt%optcell
     call chkdpr(1,1,cond_string,cond_values,ierr,'ecutsm',dt%ecutsm,1,tol8,iout)
   end if

!  ecutwfn <= ecut. This is also needed for the correct evaluation
!  of the Kleynman-Bylander form factors as the spline in Psps% is done with ecut
!  while we need |q+G| up to ecut. enlargement due to the q is already
!  taken into account by enlarging the spline mesh by around 20%.
   if ( ANY(optdriver == [RUNL_SCREENING, RUNL_SIGMA, RUNL_BSE]) ) then
     call chkdpr(0,0,cond_string,cond_values,ierr,'ecutwfn',dt%ecuteps,1,0.0_dp,iout)
     if(dt%ecut<dt%ecutwfn-tol8)then
       write(msg,'(a,es16.6,a,es16.6,a,6a)')&
        'The values of ecut and ecutwfn are ', dt%ecut,' and ',dt%ecutwfn,ch10,&
        'One expects ecutwfn to be smaller or equal to ecut.',ch10,&
        'Action: adjust ecutwfn with ecut.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

   ! Check variables used to specify k-points in self-energy.
   if (dt%nkptgw /= 0 .and. (any(dt%sigma_erange > zero .or. dt%gw_qprange /= 0))) then
     MSG_ERROR_NOSTOP("nkptw cannot be used with sigma_erange or gw_qprange", ierr)
   end if
   if (any(dt%sigma_erange > zero) .and. dt%gw_qprange /= 0) then
     MSG_ERROR_NOSTOP("sigma_erange and gw_qprange are mutually exclusive", ierr)
   end if

!  efmas
   if(optdriver==RUNL_RESPFN) then !.and.usepaw==1)then
     cond_string(1)='optdriver' ; cond_values(1)=1
     cond_string(2)='usepaw'    ; cond_values(2)=0 !usepaw
     cond_string(3)='ieig2rf'   ; cond_values(3)=1
     cond_string(4)='nsym'      ; cond_values(4)=1
     !cond_string(5)='useylm'    ; cond_values(5)=1
     call chkint_eq(1,4,cond_string,cond_values,ierr,'efmas',dt%efmas,2,(/0,1/),iout)
     if (dt%paral_rf==1) then
       cond_string(1)='paral_rf' ; cond_values(1)=1
       call chkint_eq(1,1,cond_string,cond_values,ierr,'efmas',dt%efmas,1,(/0/),iout)
     end if
   end if

!  efmas_calc_dirs
   if(dt%efmas==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'efmas_calc_dirs',dt%efmas_calc_dirs,7,(/-3,-2,-1,0,1,2,3/),iout)
   end if

!  efmas_deg
   if(dt%efmas==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'efmas_deg',dt%efmas_deg,2,(/0,1/),iout)
   end if

!  efmas_deg_tol
   if(dt%efmas==1) then
     call chkdpr(0,0,cond_string,cond_values,ierr,'efmas_deg_tol',dt%efmas_deg_tol,1,0.0_dp,iout)
   end if

!  efmas_dim
   if(dt%efmas==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'efmas_dim',dt%efmas_dim,3,(/1,2,3/),iout)
   end if

!  efmas_n_dirs
   if(dt%efmas==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'efmas_n_dirs',dt%efmas_n_dirs,0,iout)
   end if

!  efmas_ntheta
   if(dt%efmas==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'efmas_ntheta',dt%efmas_ntheta,1,iout)
   end if

!  enable_mpi_io
   if(dt%iomode==IO_MODE_MPI) then
     cond_string(1)='iomode' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'enable_mpi_io',xmpi_mpiio,1,(/1/),iout)
   end if

   !  eph variables
   if (optdriver == RUNL_EPH) then
     cond_string(1)='optdriver'; cond_values(1)=RUNL_EPH
     call chkint_eq(1,1,cond_string,cond_values,ierr,'eph_task',dt%eph_task, &
       14, [0, 1, 2, -2, 3, 4, -4, 5, -5, 6, 7, 15, -15, 16], iout)

     if (any(dt%ddb_ngqpt <= 0)) then
       MSG_ERROR_NOSTOP("ddb_ngqpt must be specified when performing EPH calculations.", ierr)
     end if
     if (dt%eph_task == 1 .and. dt%ph_nqpath <= 0) then
       MSG_ERROR("When eph_task == 1, the q-path for the linewidth must be specified via ph_nqpath and ph_qpath")
     end if

     if (dt%eph_task == 2 .and. dt%irdwfq == 0 .and. dt%getwfq == 0) then
       MSG_ERROR_NOSTOP('Either getwfq or irdwfq must be non-zero in order to compute the gkk', ierr)
     end if
     if (any(dt%eph_task == [-5, -14])) then
       ABI_CHECK(dt%ph_nqpath > 0, "ph_nqpath must be specified when eph_task in [-5, -14]")
     end if

     cond_string(1)='optdriver' ; cond_values(1)=RUNL_EPH
     call chkint_eq(1,1,cond_string,cond_values,ierr,'eph_frohlichm',dt%eph_frohlichm,2,[0,1],iout)

     if (dt%eph_stern /= 0) then
       ! Check requirements for Sternheimer.
       MSG_ERROR_NOSTOP_IF(dt%tolwfr == zero, "tolwfr must be specified when eph_stern /= 0", ierr)
       if (dt%getpot_path == ABI_NOFILE) then
         MSG_ERROR_NOSTOP(" getpot_path is required when eph_stern /= 0", ierr)
       end if
       if (all(dt%sigma_bsum_range /= 0)) then
         MSG_ERROR_NOSTOP("sigma_bsum_range cannot be used when eph_stern /= 0", ierr)
       end if
     end if
   end if

   if (any(dt%eph_np_pqbks /= 0)) then
     ! Perform basic consistency check for MPI grid.
     ! (q-points and k-points will be computed at runtime so cannot perform checks at this level.
     if (product(dt%eph_np_pqbks) /= all_nprocs) then
       write(msg, "(a,i0,3a, 6(a,1x,i0))") &
         "Cannot create 5d Cartesian grid with nprocs: ", all_nprocs, ch10, &
         "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
         "nprocs_pert (", dt%eph_np_pqbks(1), ") x nprocs_qpt (", dt%eph_np_pqbks(2), &
         ") x nprocs_bsum (", dt%eph_np_pqbks(3), ") x nprocs_kcalc (", dt%eph_np_pqbks(4), &
         ") x nprocs_spin (", dt%eph_np_pqbks(5), ") != ", all_nprocs
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

     ! Check spin
     if (dt%eph_np_pqbks(5) /= 0) then
       if (dt%nspinor == 2) then
         MSG_ERROR_NOSTOP("Spin parallelism cannot be used when nspinor == 2", ierr)
       else if (dt%nspinor == 1 .and. dt%eph_np_pqbks(5) > dt%nsppol) then
         MSG_ERROR_NOSTOP("Number of procs for spin cannot be greater than nsppol", ierr)
       end if
     end if

    if (dt%eph_np_pqbks(1) /= 0 .and. dt%eph_np_pqbks(1) > 3 * dt%natom ) then
      MSG_ERROR_NOSTOP("Number of procs for perturbation parallelism  cannot be greater than 3 * natom", ierr)
    end if
   end if

   ! exchmix
   call chkdpr(0,0,cond_string,cond_values,ierr,'exchmix',dt%exchmix,1,0.0_dp,iout)

   ! extrapwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'extrapwf',dt%extrapwf,3,(/0,1,2/),iout)
   if (dt%extrapwf>0.and.dt%densfor_pred<5) then
     write(msg,'(3a)')&
     'extrapwf keyword (extrapolation of WF) is only compatible with',ch10,&
     'densfor_pred=5 or 6; please change densfor_pred value.'
     MSG_ERROR_NOSTOP(msg,ierr)
     ! MT oct 14: Should use chkint_eq but the msg is not clear enough
   end if

   ! fermie_nest
   call chkdpr(0,0,cond_string,cond_values,ierr,'fermie_nest',dt%fermie_nest,1,0.0_dp,iout)

   ! fftgw
   call chkint_eq(0,0,cond_string,cond_values,ierr,'fftgw',dt%fftgw,8, [00,01,10,11,20,21,30,31],iout)

   !  fockoptmix
   call chkint_eq(0,0,cond_string,cond_values,ierr,'fockoptmix',&
     dt%fockoptmix,12,(/0,1,11,201,211,301,401,501,601,701,801,901/),iout)
   if(dt%paral_kgb/=0)then
     cond_string(1)='paral_kgb' ; cond_values(1)=dt%paral_kgb
     ! Make sure that dt%fockoptmix is 0, 1 or 11 (wfmixalg==0)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'fockoptmix',dt%fockoptmix,3,(/0,1,11/),iout)
   end if

   ! frzfermi
   call chkint_eq(0,0,cond_string,cond_values,ierr,'frzfermi',dt%frzfermi,2,(/0,1/),iout)

   ! fxcartfactor
   call chkdpr(0,0,cond_string,cond_values,ierr,'fxcartfactor',dt%fxcartfactor,1,zero,iout)

   ! ga_algor
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ga_algor',dt%ga_algor,3,(/1,2,3/),iout)

   ! ga_fitness
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ga_fitness',dt%ga_fitness,3,(/1,2,3/),iout)

   ! ga_opt_percent
   call chkdpr(0,0,cond_string,cond_values,ierr,'ga_opt_percent',dt%ga_opt_percent,1,tol8,iout)

   ! getxred
   if(dt%getxcart/=0)then
     cond_string(1)='getxcart' ; cond_values(1)=dt%getxcart
     ! Make sure that dt%getxred is 0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'getxred',dt%getxred,1,(/0/),iout)
   end if

   ! goprecon
   call chkint_eq(0,0,cond_string,cond_values,ierr,'goprecon',dt%goprecon,4,(/0,1,2,3/),iout)

   ! gpu_devices
   if (dt%use_gpu_cuda==1) then
     if (all(gpu_devices(:)==-2)) then
       gpu_devices(:)=dt%gpu_devices(:)
     else if (any(dt%gpu_devices(:)/=gpu_devices(:))) then
       write(msg,'(3a)')&
        'GPU device(s) selection cannot be different from one dataset to another!',ch10,&
        'Action: change gpu_devices in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

   ! gw_invalid_freq
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gw_invalid_freq',dt%gw_invalid_freq,3,(/0,1,2/),iout)

   ! gw_sctype
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gw_sctype',dt%gw_sctype,&
     4,(/GWSC_one_shot,GWSC_only_W,GWSC_only_G,GWSC_both_G_and_W/),iout)

   ! gw_sigxcore
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gw_sigxcore',dt%gw_sigxcore,2,[0,1],iout)

   ! gwcomp
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwcomp',dt%gwcomp,2,[0,1],iout)
   if (dt%gwcomp/=0) then
     if (optdriver==RUNL_SCREENING .and. (dt%awtr /=1 .or. dt%spmeth /=0)) then
       write(msg,'(3a)' )&
        'When gwcomp /= 0, the Adler-Wiser formula with time-reversal should be used',ch10,&
         'Action: set awtr to 1 or/and spmeth to 0'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

     ! Extrapolar trick with HF, SEX and COHSEX is meaningless for Sigma
     if(optdriver == RUNL_SIGMA) then
       mod10=MOD(dt%gwcalctyp,10)
       if ( ANY(mod10 == [SIG_HF, SIG_SEX, SIG_COHSEX]) ) then
         write(msg,'(3a)' )&
         'gwcomp/=0, is meaningless in the case of HF, SEX or COHSEX calculations. ',ch10,&
         'Action: set gwcomp to 0 or change gwcalctyp'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end if
     if (optdriver==RUNL_SIGMA .and. ALL( dt%ppmodel /= [0,1,2] )) then
       write(msg,'(a,i0,a)')&
         'The completeness trick cannot be used when ppmodel is ',dt%ppmodel,'. It should be set to 0, 1 or 2. '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwmem',dt%gwmem,4,[0,1,10,11],iout)
   call chkint_eq(0,0,cond_string,cond_values,ierr,'gwpara',dt%gwpara,3,[0,1,2],iout)

   ! gwrpacorr
   if(dt%gwrpacorr>0) then
     mod10=MOD(dt%gwcalctyp,10)
     if (optdriver /= RUNL_SCREENING) then
       write(msg,'(3a)' )&
       'gwrpacorr>0 can only be used when calculating the screening',ch10,&
       'Action: set gwrpacorr to 0 or optdriver to 3'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if( mod10 /= SIG_GW_AC ) then
       write(msg,'(3a)' )&
       'gwrpacorr > 0 can only be used with purely imaginary frequencies',ch10,&
       'Action: set gwrpacorr to 0 or change gwcalctyp'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

   ! gwls_stern_kmax
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_stern_kmax',dt%gwls_stern_kmax,1,iout)

   ! gwls_npt_gauss_quad
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_npt_gauss_quad',dt%gwls_npt_gauss_quad,1,iout)

   ! gwls_diel_model
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_diel_model',dt%gwls_diel_model,1,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'gwls_diel_model',dt%gwls_diel_model,3,iout)

   ! gwls_model_parameter
   call chkdpr(0,0,cond_string,cond_values,ierr,'gwls_model_parameter',dt%gwls_model_parameter,1,zero,iout)

   ! gwls_print_debug
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_print_debug',dt%gwls_print_debug,0,iout)

   ! gwls_nseeds
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_nseeds',dt%gwls_nseeds,1,iout)

   ! gwls_kmax_complement
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_kmax_complement',dt%gwls_kmax_complement,0,iout)

   ! gwls_kmax_poles
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_kmax_poles',dt%gwls_kmax_poles,0,iout)

   ! gwls_kmax_analytic
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_kmax_analytic',dt%gwls_kmax_analytic,0,iout)

   ! gwls_kmax_numeric
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_kmax_numeric',dt%gwls_kmax_numeric,0,iout)

   ! gwls_band_index
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_band_index',dt%gwls_band_index,1,iout)

   ! gwls_exchange
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_exchange',dt%gwls_exchange,0,iout)

   ! gwls_correlation
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_correlation',dt%gwls_correlation,0,iout)

   ! gwls_first_seed
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_first_seed',dt%gwls_first_seed,1,iout)

   ! gwls_recycle
   call chkint_ge(0,0,cond_string,cond_values,ierr,'gwls_recycle',dt%gwls_recycle,0,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'gwls_recycle',dt%gwls_recycle,2,iout)

   ! iatsph must between 1 and natom
   maxiatsph=maxval(dt%iatsph(1:dt%natsph))
   miniatsph=minval(dt%iatsph(1:dt%natsph))
   call chkint_ge(0,0,cond_string,cond_values,ierr,'iatsph',miniatsph,1,iout)
   call chkint_le(0,0,cond_string,cond_values,ierr,'iatsph',maxiatsph,natom,iout)

   ! icoulomb
   call chkint_eq(0,0,cond_string,cond_values,ierr,'icoulomb',dt%icoulomb,3,(/0,1,2/),iout)
   if (dt%nspden > 2) then
     cond_string(1)='nspden' ; cond_values(1)=nspden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'icoulomb',dt%icoulomb,1,(/0/),iout)
   end if

   ! ieig2rf
   if(optdriver==RUNL_RESPFN.and.usepaw==1)then
     cond_string(1)='optdriver' ; cond_values(1)=1
     cond_string(2)='usepaw'    ; cond_values(2)=usepaw
     call chkint_eq(1,2,cond_string,cond_values,ierr,'ieig2rf',dt%ieig2rf,1,(/0/),iout)
   end if
   if(optdriver==RUNL_RESPFN.and.dt%paral_rf==1)then
     cond_string(1)='paral_rf' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ieig2rf',dt%ieig2rf,1,(/0/),iout)
   end if

   ! imgmov
   call chkint_eq(0,0,cond_string,cond_values,ierr,'imgmov',dt%imgmov,9,(/0,1,2,4,5,6,9,10,13/),iout)
   if (dt%imgmov>0 .and. dt%imgmov/=6) then ! when imgmov>0, except imgmov==6, allow only ionmov0 and optcell 0 (temporary)
     cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ionmov',dt%ionmov,1,(/0/),iout)
     if (dt%imgmov==9.or.dt%imgmov==10.or.dt%imgmov==13) then
       cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
       !Temporarily deactivate NPT algorithms (not yet usable)
       call chkint_eq(1,1,cond_string,cond_values,ierr,'optcell',dt%optcell,1,(/0/),iout)
     else
       cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
       call chkint_eq(1,1,cond_string,cond_values,ierr,'optcell',dt%optcell,1,(/0/),iout)
     end if
   end if

   ! imgwfstor
   call chkint_eq(0,0,cond_string,cond_values,ierr,'imgwfstor',dt%imgwfstor,2,(/0,1/),iout)
   if (dt%extrapwf/=0) then ! extrapwf/=0 not allowed presently with imgwfstor
     cond_string(1)='extrapwf' ; cond_values(1)=dt%extrapwf
     call chkint_eq(1,1,cond_string,cond_values,ierr,'imgwfstor',dt%imgwfstor,1,(/0/),iout)
   endif
   if (dt%ntimimage<=1) then ! imgwfstor activate only when there is more than one time step for images
     cond_string(1)='ntimimage' ; cond_values(1)=dt%ntimimage
     call chkint_eq(1,1,cond_string,cond_values,ierr,'imgwfstor',dt%imgwfstor,1,(/0/),iout)
   endif

!  intxc
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
     ! Make sure that dt%intxc is 0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'intxc',dt%intxc,1,(/0/),iout)
   end if
!  TEMPORARY
   if(optdriver==RUNL_RESPFN)then ! Make sure that dt%intxc is 0
     cond_string(1)='optdriver' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'intxc',dt%intxc,1,(/0/),iout)
   end if

   ! ionmov
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ionmov',&
     dt%ionmov,24,(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25,26,27/),iout)

   ! When optcell/=0, ionmov must be 2, 3, 13, 22 or 25 (except if imgmov>0)
   if(dt%optcell/=0)then
     if (dt%imgmov==0) then
       cond_string(1)='optcell' ; cond_values(1)=dt%optcell
       call chkint_eq(1,1,cond_string,cond_values,ierr,'ionmov',dt%ionmov,6,(/2,3,13,15,22,25/),iout)
     else
       cond_string(1)='optcell' ; cond_values(1)=dt%optcell
       call chkint_eq(1,1,cond_string,cond_values,ierr,'ionmov',dt%ionmov,1,(/0/),iout)
     end if
   end if
   if (dt%ionmov == 13) then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
!    Make sure that nnos is not null
     call chkint_ge(1,1,cond_string,cond_values,ierr,'nnos',dt%nnos,1,iout)
   end if

!  iprcel
   call chkint(0,0,cond_string,cond_values,ierr,'iprcel',dt%iprcel,1,(/0/),1,21,iout)   !  0 or superior to 21
   if(nsppol==2 .and. (dt%occopt>=3 .and. dt%occopt<=8).and.mod(dt%iprcel,10)>49 )then
     write(msg,'(5a)')&
     'For spin-polarized metallic systems (occopt>3),',ch10,&
     'only RPA dielectric matrix can be evaluated) !',ch10,&
     'Action: change iprcel value in input file (mod(iprcel,100)<50) !'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if(dt%npspinor>1.and.dt%iprcel>0)then
     write(msg,'(5a)')&
     'When parallelization over spinorial components is activated (npspinor>1),',ch10,&
     'only model dielectric function is allowed (iprcel=0) !',ch10,&
     'Action: change iprcel value in input file !'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

   ! irandom
   call chkint_eq(0,0,cond_string,cond_values,ierr,'irandom',dt%irandom,3, [1,2,3], iout)

   ! iscf
   if (usewvl ==0) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,&
      'iscf',dt%iscf,18, [-3,-2,-1,1,2,3,4,5,6,7,11,12,13,14,15,16,17,22], iout)
   else
!    If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
     wvlbigdft=(dt%usewvl==1.and.dt%wvl_bigdft_comp==1)
     cond_string(1)='wvl_bigdft_comp' ; cond_values(1)=dt%wvl_bigdft_comp
     if(wvlbigdft) then
       call chkint_eq(1,1,cond_string,cond_values,ierr,&
&       'iscf',dt%iscf,15,(/0,1,2,3,4,5,6,7,11,12,13,14,15,16,17/),iout)
     else
       call chkint_eq(1,1,cond_string,cond_values,ierr,&
&       'iscf',dt%iscf,18,(/-3,-2,-1,1,2,3,4,5,6,7,11,12,13,14,15,16,17,22/),iout)
     end if
!    If wvl+metal, iscf cannot be 0
     if (dt%occopt>2) then
       cond_string(1)='occopt' ; cond_values(1)=dt%occopt
       call chkint_eq(1,1,cond_string,cond_values,ierr,&
&       'iscf',dt%iscf,18,(/-3,-2,-1,1,2,3,4,5,6,7,11,12,13,14,15,16,17,22/),iout)
     end if
   end if
!  If ionmov==4, iscf must be 2, 12, 5 or 6.
   if(dt%ionmov==4)then
     cond_string(1)='ionmov' ; cond_values(1)=4
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,4,(/2,12,5,6/),iout)
   end if
!  If PAW, iscf cannot be -1, 11
   if (usepaw==1 .and. usewvl==0) then
     cond_string(1)='PAW' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,11,(/-3,-2,2,3,4,7,12,13,14,17,22/),iout)
   end if
!  Mixing on density is only allowed for GS calculations or for drivers where it is not used.
   if(optdriver /= RUNL_GSTATE .and. all(optdriver/=[RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE,RUNL_EPH,&
&     RUNL_WFK,RUNL_NONLINEAR])) then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_le(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,9,iout)
   end if
!  When pawoptmix=1 and nspden=4, iscf must be >=10
   if(dt%pawoptmix/=0.and.nspden==4)then
     cond_string(1)='nspden'    ; cond_values(1)=nspden
     cond_string(2)='pawoptmix' ; cond_values(2)=dt%pawoptmix
     call chkint_ge(2,2,cond_string,cond_values,ierr,'iscf',dt%iscf,10,iout)
   end if

!  istatimg
   call chkint_eq(0,0,cond_string,cond_values,ierr,'istatimg',dt%istatimg,2,(/0,1/),iout)
   if (dt%string_algo==2) then
     cond_string(1)='string_algo' ; cond_values(1)=dt%string_algo
     call chkint_eq(1,1,cond_string,cond_values,ierr,'istatimg',dt%istatimg,1,(/1/),iout)
   end if

!  istwfk
   if(dt%usefock==1 .and. dt%optdriver/=RUNL_SIGMA .and. mod(dt%wfoptalg,10)/=5 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
     write(msg,'(3a)' )&
      'When usefock==1, unless sigma calculation, all the components of istwfk must be 1.',ch10,&
      'Action: set istwfk to 1 for all k-points'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

   if(dt%usewvl==1 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
     write(msg,'(3a)' )&
      'When usewvl==1, all the components of istwfk must be 1.',ch10,&
      'Action: set istwfk to 1 for all k-points'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

   if(response==1 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
     ! Force istwfk to be 1 for RF calculations
     ! Other choices cannot be realized yet, because of the ddk perturbation.
     write(msg,'(5a)' )&
     'When response==1, all the components of istwfk must be 1.',ch10,&
     'Not yet programmed for time-reversal symmetry.',ch10,&
     'Action: set istwfk to 1 for all k-points'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if(dt%nbandkss/=0 .and. dt%kssform/=3 .and. maxval( abs(dt%istwfk(1:nkpt)-1) ) >0)then
     write(msg,'(5a)' )&
     'When nbandkss/=0 and kssform/=3 all the components of istwfk must be 1.',ch10,&
     'Not yet programmed for time-reversal symmetry.',ch10,&
     'Action: set istwfk to 1 for all k-points'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if(dt%berryopt/=0 .and. maxval(dt%istwfk(:))/=1)then
     write(msg,'(5a)' )&
     'When berryopt/=0, all the components of istwfk must be 1.',ch10,&
     'Not yet programmed for time-reversal symmetry.',ch10,&
     'Action: set istwfk to 1 for all k-points'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if (dt%optdriver==RUNL_GSTATE) then
     if ((dt%wfoptalg==4.or.dt%wfoptalg==14.or.dt%wfoptalg==114).and.maxval(dt%istwfk(:)-2)>0) then
       write(msg, '(a,a,a,a,a)' )&
       'Only the gamma point can use time-reversal and wfoptalg=4 or 14',ch10,&
       'Action: put istwfk to 1 or remove k points with half integer coordinates ',ch10,&
       'Also contact ABINIT group to say that you need that option.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
!     if ((dt%wfoptalg==4.or.dt%wfoptalg==14).and.any(dt%istwfk(:)==2) .and.dt%paral_kgb==1.and.fftalg/=401.and.fftalg/=312) then
!       write(msg, '(a,i3,a,a,a)' )&
!&       ' For istwfk=2, the value fftalg= ',fftalg, &
!&       ' is not allowed in case of wfoptalg=4 or 14 !', ch10,&
!&       ' Change if to fftalg=401.'
!       MSG_ERROR_NOSTOP(msg, ierr)
!     end if
   end if

!  ixc
   call chkint(0,0,cond_string,cond_values,ierr,&
&   'ixc',dt%ixc,33,(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,21,22,23,24,26,27,31,32,33,34,40,41,42,50/),-1,0,iout) ! One of the values, or negative
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
!    Make sure that ixc is 1, 7, 8, 20, 21 or 22 or negative
     call chkint(1,1,cond_string,cond_values,ierr,'ixc',dt%ixc,6,(/1,7,8,20,21,22/),-1,0,iout)
   end if
   if(response==1)then
     cond_string(1)='response' ; cond_values(1)=1
!    Make sure that ixc is between 0 and 9, or 11, 12, 14, 15, 23 or 24 or negative
     call chkint(1,1,cond_string,cond_values,ierr,&
&     'ixc',dt%ixc,16,(/0,1,2,3,4,5,6,7,8,9,11,12,14,15,23,24/),-1,0,iout)
   end if
   if(nspden/=1)then
     cond_string(1)='nspden' ; cond_values(1)=nspden
!    Make sure that ixc is 0, 1 , the gga, or Fermi-Amaldi, or negative
     call chkint(1,1,cond_string,cond_values,ierr,&
&     'ixc',dt%ixc,24,(/0,1,7,8,9,11,12,13,14,15,16,17,20,23,24,26,27,31,32,33,34,40,41,42/),-1,0,iout)
   end if
   if(dt%usepaw>0.and.dt%ixc<0) then
     if (libxc_functionals_is_hybrid()) then
       MSG_ERROR_NOSTOP('Meta-GGA functionals are not compatible with PAW!', ierr)
     end if
   end if
   if (dt%usepaw>0.and.(dt%ixc==-427.or.dt%ixc==-428)) then
     MSG_WARNING('Range-separated Hybrid Functionals have not been extensively tested in PAW!!!')
   end if
   allow=(dt%ixc > 0).and.(dt%ixc /= 3).and.(dt%ixc /= 7).and.(dt%ixc /= 8)
   if(.not.allow)then
     allow=(dt%ixc < 0).and.(libxc_functionals_is_hybrid().or.libxc_functionals_ismgga())
   end if
   if(allow)then
     cond_string(1)='ixc' ; cond_values(1)=dt%ixc
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if

!  ixcpositron
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ixcpositron',dt%ixcpositron,8,(/0,-1,1,11,2,3,31,4/),iout)

!  ixcrot
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ixcrot',dt%ixcrot,3,(/1,2,3/),iout)

!  tim1rev
   call chkint_eq(0,0,cond_string,cond_values,ierr,'tim1rev',dt%tim1rev,2,(/0,1/),iout)

!  kptnrm and kpt
!  Coordinates components must be between -1 and 1.
   if(dt%kptnrm<1.0-1.0d-10)then
     write(msg, '(a,es22.14,a,a,a)' )&
      'The input variable kptnrm is',dt%kptnrm,' while it must be >=1.0_dp.',ch10,&
      'Action: change the input variable kptnrm.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   do ikpt=1,nkpt
     do mu=1,3
       if ( abs(dt%kpt(mu,ikpt))> dt%kptnrm*1.0000001_dp ) then
         write(msg, '(a,i5,a,a,a,a,3es22.14,a,a,a,a)' )&
          'For k point number',ikpt,'  the reduced coordinates',ch10,&
          'generated by the input variables kpt and kptnrm are',ch10,&
          dt%kpt(1,ikpt)/dt%kptnrm,dt%kpt(2,ikpt)/dt%kptnrm,dt%kpt(3,ikpt)/dt%kptnrm,ch10,&
          'while they must be between -1.0_dp and 1.0_dp (included).',ch10,&
          'Action: check kpt and kptnrm in the input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end do

!  jellslab
   call chkint_eq(0,0,cond_string,cond_values,ierr,'jellslab',dt%jellslab,2,(/0,1/),iout)

   if (dt%jellslab==1) then
     if(dt%nimage>1)then
       cond_string(1)='nimage' ; cond_values(1)=dt%nimage
       call chkint_eq(1,1,cond_string,cond_values,ierr,'jellslab',dt%jellslab,1,(/0/),iout)
     end if
!    slabwsrad must be positive
     cond_string(1)='jellslab' ; cond_values(1)=dt%jellslab
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabwsrad',dt%slabwsrad,1,zero,iout)
!    slabzbeg must be positive
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabzbeg',dt%slabzbeg,1,zero,iout)
!    slabzend must be bigger than slabzbeg
     call chkdpr(1,0,cond_string,cond_values,ierr,'slabzend',dt%slabzend,1,dt%slabzbeg,iout)
!    rprimd(3,3) must be bigger than slabzend
     call chkdpr(1,0,cond_string,cond_values,ierr,'rprimd33',rprimd(3,3),1,dt%slabzend,iout)
!    Third real space primitive translation has to be orthogonal to the other ones,
!    actually, for convenience it is useful that rprimd is something like:
!    a  b  0
!    c  d  0
!    0  0  e
     if(abs(rprimd(1,3))+abs(rprimd(2,3))+abs(rprimd(3,1))+abs(rprimd(3,2))>tol12) then
       write(msg,'(3a)')&
        'Third real space vector is not orthogonal to the other ones,',ch10,&
        'this is needed to use jellium'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    Atoms have to be placed in the vacuum space
     do iatom=1,natom
       zatom=(dt%xred_orig(3,iatom,intimage)-anint(dt%xred_orig(3,iatom,intimage)-half+tol6))*rprimd(3,3)
       if(abs(zatom-dt%slabzbeg)<tol8 .or. abs(zatom-dt%slabzend)<tol8) then
         if(dt%znucl(dt%typat(iatom))>tol6) then
           write(msg,'(a,i0,a)')'atom number=',iatom,' lies precisely on the jellium edge !'
           MSG_WARNING(msg)
         end if
         cycle
       end if
       if(zatom>dt%slabzbeg .and. zatom<dt%slabzend) then
         write(msg,'(a,i0,a)')' atom number=',iatom,' is inside the jellium slab.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end if

!  kssform
   call chkint_eq(0,0,cond_string,cond_values,ierr,'kssform',dt%kssform,3,(/0,1,3/),iout)

   if (dt%kssform/=0 .and. dt%nbandkss/=0) then ! Check for outkss limitations.
     call wrtout(std_out," Checking if input is consistent with KSS generation",'COLL')
     call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,1,(/0/),iout)
     call chkint_eq(0,0,cond_string,cond_values,ierr,'iomode',dt%iomode,2,(/IO_MODE_FORTRAN,IO_MODE_ETSF/),iout)
   end if

!  localrdwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,2,(/0,1/),iout)
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if
   if(dt%mkqmem==0)then
     cond_string(1)='mkqmem' ; cond_values(1)=dt%mkqmem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if
   if(dt%mk1mem==0)then
     cond_string(1)='mk1mem' ; cond_values(1)=dt%mk1mem
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if
   if(dt%iomode==IO_MODE_MPI)then
     cond_string(1)='iomode' ; cond_values(1)=dt%iomode
     call chkint_eq(1,1,cond_string,cond_values,ierr,'localrdwf',dt%localrdwf,1,(/1/),iout)
   end if


!  LOTF
#if defined HAVE_LOTF
   if (dt%ionmov==23) then
     write(msg, '(a,a)' ) ch10,&
&     '=== LOTF METHOD ================================================================'
     call wrtout(ab_out,msg,'COLL')
     cond_string(1)='ionmov' ; cond_values(1)=23
     call chkint_eq(0,1,cond_string,cond_values,ierr,'lotf_classic',dt%lotf_classic,1,(/5/),iout)
     cond_string(1)='ionmov' ; cond_values(1)=23
     call chkint_ge(0,1,cond_string,cond_values,ierr,'lotf_nitex',dt%lotf_nitex,1,iout)
     cond_string(1)='ionmov' ; cond_values(1)=23
     call chkint_ge(0,1,cond_string,cond_values,ierr,'lotf_nneigx',dt%lotf_nneigx,2,iout)
     cond_string(1)='ionmov' ; cond_values(1)=23
     call chkint_eq(0,1,cond_string,cond_values,ierr,'lotf_version',dt%lotf_version,1,(/2/),iout)
   end if
#endif

!  magconon
   call chkint_eq(0,0,cond_string,cond_values,ierr,'magconon',dt%magconon,3,(/0,1,2/),iout)
!!  impose nspden 4 for the moment and spinors
!   if (dt%magconon == 1) then
!     if (dt%nspinor /= 2 .or. dt%nspden /= 4) then
!       write (msg, '(4a)') &
!&       ' magnetization direction constraint is only compatible with non-collinear calculations', ch10,&
!&       ' Action: set nspinor 2 and nspden 4 in the input file.'
!       MSG_ERROR_NOSTOP(msg,ierr)
!     end if
!   end if

!  macro_uj
   if(dt%macro_uj/=0) then
     if (dt%ionmov/=0) then
       write(msg, '(3a,i2,2a,i2,3a)' )&
        'Determination of U can not be combined with ionic movements.',ch10,&
        'Here  ionmov= ',dt%ionmov,ch10,&
        'and macro_uj=',dt%macro_uj,'.',ch10,&
        'Action: change ionmov in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     else if (dt%nstep<3) then
       write(msg, '(3a,i1,2a,i2,3a)' )&
        'Determination of U needs at least 3 scf steps:',ch10,&
        ' nstep = ',dt%nstep,ch10,&
        ' and macro_uj=',dt%macro_uj,'.',ch10,&
        'Action: increase nstep in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  mep_solver
   call chkint_eq(0,0,cond_string,cond_values,ierr,'mep_solver',dt%mep_solver,5,(/0,1,2,3,4/),iout)
!  String method
   if(dt%imgmov==2) then
     cond_string(1)='imgmov'      ; cond_values(1)=dt%imgmov
!    Some restriction for the solver
     if(dt%string_algo==0)then
       cond_string(2)='string_algo' ; cond_values(2)=dt%string_algo
       call chkint_eq(1,1,cond_string,cond_values,ierr,'mep_solver',dt%mep_solver,1,(/0/),iout)
     end if
     if(dt%string_algo==1.or.dt%string_algo==2)then
       cond_string(2)='string_algo' ; cond_values(2)=dt%string_algo
       call chkint_eq(1,1,cond_string,cond_values,ierr,'mep_solver',dt%mep_solver,2,(/0,4/),iout)
     end if
   end if
!  NEB
   if(dt%imgmov==5)then
     cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
!    Some restriction for the solver
     call chkint_eq(1,1,cond_string,cond_values,ierr,'mep_solver',dt%mep_solver,4,(/0,1,2,3/),iout)
!    Static image energy is needed if spring constant is variable
     if (abs(dt%neb_spring(1)-dt%neb_spring(2))>=tol8.and.dt%istatimg==0) then
       write(msg, '(7a)' )&
        'When using variable NEB spring constants (which is the default for CI-NEB),',ch10,&
        'all the energies of the cell images are needed (including static images!).',ch10,&
        'You cannot use istatimg=0!',ch10,&
        'Action: put istatimg=1 in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
!    Static image energy is needed for CI-NEB or improved tangent
     if ((dt%neb_algo==1.or.dt%neb_algo==2).and.dt%istatimg==0) then
       write(msg, '(7a)' )&
        'When using Improved-tangent-NEB or CI-NEB,',ch10,&
        'all the energies of the cell images are needed (including static images!).',ch10,&
        'You cannot use istatimg=0!',ch10,&
        'Action: put istatimg=1 in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  mffmem
   call chkint_eq(0,0,cond_string,cond_values,ierr,'mffmem',dt%mffmem,2,(/0,1/),iout)

!  mixalch_orig
!  For each type of atom, the sum of the psp components must be one.
   do iimage=1,dt%nimage
     if(dt%ntypalch>0)then
       do itypat=1,dt%ntypalch
         sumalch=sum(dt%mixalch_orig(:,itypat,iimage))
         if(abs(sumalch-one)>tol10)then
           if(dt%npspalch<=6)then
             write(msg, '(2a,6es12.4)' )ch10,' chkinp : mixalch(:,itypat,iimage)=',dt%mixalch_orig(:,itypat,iimage)
           end if
           call wrtout(iout,msg,'COLL')
           call wrtout(std_out,  msg,'COLL')
           write(msg, '(a,i4,2a,i4,2a,f8.2,4a)' )&
            'For the alchemical atom number',itypat,ch10,&
            'image number',iimage,ch10,&
            'the sum of the pseudopotential coefficients is',sumalch,ch10,&
            'while it should be one.',ch10,&
            'Action: check the content of the input variable mixalch.'
           MSG_ERROR_NOSTOP(msg, ierr)
         end if
       end do
     end if
   end do

!  mixesimgf
!  The sum of the mixing image factors must be one
   if(dt%imgmov==6)then
     summix=sum(dt%mixesimgf(1:dt%nimage))
     if(abs(summix-one)>tol10)then
       write(msg, '(2a,20es12.4)' )ch10,' chkinp : mixesimgf(1:dt%nimage)=',dt%mixesimgf(1:dt%nimage)
       call wrtout(iout,msg,'COLL')
       call wrtout(std_out,  msg,'COLL')
       write(msg, '(a,es12.4,4a)' )&
        'The sum of the mixing image factors is',summix,ch10,&
        'while it should be one.',ch10,&
        'Action: check the content of the input variable mixesimgf.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  natom
   if(dt%prtgeo>0)then
     cond_string(1)='prtgeo' ; cond_values(1)=dt%prtgeo
     call chkint_le(1,1,cond_string,cond_values,ierr,'natom',natom,9999,iout)
   end if

!  nband
!  Make sure all nband(nkpt) are > 0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       if (dt%nband(ikpt+(isppol-1)*nkpt)<=0) then
         cond_string(1)='ikpt' ; cond_values(1)=ikpt
         cond_string(2)='isppol' ; cond_values(2)=isppol
         call chkint_ge(0,2,cond_string,cond_values,ierr,'nband',dt%nband(ikpt+(isppol-1)*nkpt),1,iout)
       end if
     end do
   end do
   if(nproc/=1.and.nsppol==2.and.usewvl==0)then
     do ikpt=1,nkpt
       if (dt%nband(ikpt)/=dt%nband(ikpt+nkpt)) then
         write(msg, '(5a,i4,a,2i5,a)' )&
          'the number of bands in the spin up case must be equal to',ch10,&
          'the number of bands in the spin down case.',ch10,&
          'This is not the case for the k point number :',ikpt,&
          'The number of bands spin up and down are :',dt%nband(ikpt),dt%nband(ikpt+nkpt),&
          'Action: change nband, or use the sequential version of ABINIT.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end if

!  nbandkss
!  Must be greater or equal to -1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nbandkss',dt%nbandkss,-1,iout)
!  When ionmov/=0
   if(dt%ionmov/=0 .and. dt%nbandkss/=0)then
     write(msg,'(11a)')&
      'Ions (or cell) are allowed to move (ionmov/=0),',ch10,&
      'and a _KSS file is requested (nbandkss/=0).',ch10,&
      'A _KSS file will be created at each geometry-optimisation step.',ch10,&
      'Note that this is time consuming !',ch10,&
      'Action: use datasets (one for geometry optimisation,',ch10,&
      '        one for states output).'
     MSG_WARNING(msg)
   end if

!  nbdblock
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,iout)
!  When wfoptalg==0, nbdblock must be 1
   if(mod(dt%wfoptalg,10)==0)then
     cond_string(1)='wfoptalg' ; cond_values(1)=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
   end if
!  When wfoptalg==2, nbdblock must be 1
   if(dt%wfoptalg==2)then
     cond_string(1)='wfoptalg' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
   end if
!  When wfoptalg==3, nbdblock must be 1, and iscf must be -2
   if(dt%wfoptalg==3)then
     cond_string(1)='wfoptalg' ; cond_values(1)=3
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nbdblock',dt%nbdblock,1,(/1/),iout)
     call chkint_eq(1,1,cond_string,cond_values,ierr,'iscf',dt%iscf,1,(/-2/),iout)
   end if
!  When wfoptalg==4, nbdblock must be a divisor of nband
   if(mod(dt%wfoptalg,10)==4.and.dt%optdriver==RUNL_GSTATE)then
     do isppol=1,nsppol
       do ikpt=1,nkpt
         if(mod(dt%nband(ikpt+(isppol-1)*nkpt),dt%nbdblock)/=0) then
           write(msg, '(5a)' )&
            'For the moment, when wfoptalg=4,',ch10,&
            'nband must be a multiple of nbdblock.',ch10,&
            'Action: check the value of the input variable nbdblock.'
           MSG_ERROR_NOSTOP(msg, ierr)
         end if
       end do
     end do
   end if

!  nberry
!  must be between 0 and 20
   if(dt%berryopt/=0)then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'nberry',dt%nberry,0,iout)
     call chkint_le(0,0,cond_string,cond_values,ierr,'nberry',dt%nberry,20,iout)
     if(xmpi_paral==1)then
!      MPI Parallel case
       if (dt%nberry/=0.and.dt%berryopt>0.and.&
           dt%berryopt/= 4.and.dt%berryopt/= 5.and.dt%berryopt/= 6.and.dt%berryopt/= 7.and.&
           dt%berryopt/=14.and.dt%berryopt/=15.and.dt%berryopt/=16.and.dt%berryopt/=17) then
         write(msg,'(a,a,a,a,a,i4,a,a,a)')&
          'Berry phase calculation of polarisation with positive berryopt is not',ch10,&
          'allowed in the parallel version of ABINIT.',ch10,&
          'So, the value of nberry=',dt%nberry,' is not allowed,',ch10,&
          'Action: change berryopt to negative values or change nberry, or use the sequential version.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end if
   end if

   if (dt%optcell /=0 .and. dt%berryopt == 4)  then
     write(msg,'(a,a,a,a,a,a,a,a,a,a,a,a,a)') ch10,&
      ' chkinp : WARNING -',ch10,&
      '  Constant unreduced E calculation with relaxation of cell parameters is allowed.',ch10,&
      '  But we strongly recommend users to use reduced ebar calculation (berryopt=14)',ch10,&
      '  with the relaxation of cell parameters, for internal consistency purpose.',ch10, &
      '  For more information, please refer to "M. Stengel, N.A. Spaldin and D.Vanderbilt,', ch10, &
      '  Nat. Phys., 5, 304,(2009)" and its supplementary notes.', ch10 ! [[cite:Stengel2009]]
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
   end if

   if (dt%optcell /=0 .and. (dt%berryopt == 6 ))  then
     write(msg,'(12a)') ch10,&
&     ' chkinp : WARNING -',ch10,&
&     '  Constant unreduced D calculation with relaxation of cell parameters is allowed.',ch10,&
&     '  But we strongly recommend users to use reduced d calculation (berryopt=16)',ch10,&
&     '  with the relaxation of cell parameters, for internal consistency purpose.',ch10, &
&     '  For more information, please refer to "M. Stengel, N.A. Spaldin and D.Vanderbilt,', ch10, &
&     '  Nat. Phys., 5, 304,(2009)" and its supplementary notes.' ! [[cite:Stengel2009]]
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,msg,'COLL')
   end if

!  ndynimage
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ndynimage',dt%ndynimage,1,iout)

!  neb_algo
   call chkint_eq(0,0,cond_string,cond_values,ierr,'neb_algo',dt%neb_algo,4,(/0,1,2,3/),iout)

!  nfft and nfftdg
!  Must have nfft<=nfftdg
   if (usepaw==1) then
     nfft  =dt%ngfft(1)  *dt%ngfft(2)  *dt%ngfft(3)
     nfftdg=dt%ngfftdg(1)*dt%ngfftdg(2)*dt%ngfftdg(3)
     cond_string(1)='nfft' ; cond_values(1)=nfft
     call chkint(1,1,cond_string,cond_values,ierr,'nfftdg',nfftdg,1,(/0/),1,nfft,iout) ! Must be 0 or nfft
   end if

!  diismemory
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'diismemory',dt%diismemory,1,iout)

!  nimage
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nimage',dt%nimage,1,iout)
   if (usewvl==1) then
     cond_string(1)='usewvl' ; cond_values(1)=usewvl
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if
   if (optdriver/=RUNL_GSTATE) then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if
   if (dt%tfkinfunc==2) then
     cond_string(1)='tfkinfunc' ; cond_values(1)=dt%tfkinfunc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if
   if (dt%prtxml==1) then
     cond_string(1)='prtxml' ; cond_values(1)=dt%prtxml
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nimage',dt%nimage,1,(/1/),iout)
   end if
   if (dt%imgmov==9.or.dt%imgmov==13) then
     if (dt%pitransform==1.and.(mod(dt%nimage,2)/=0)) then
       write(msg,'(6a)')ch10,&
        'Path-Integral Molecular Dynamics (imgmov=9,13)',ch10,&
        'in normal mode tranformation (pitransform=1).',ch10,&
        'requires nimage to be even!'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%imgmov==10.and.dt%pitransform>0) then
     write(msg,'(4a)')ch10,&
      'Path-Integral Molecular Dynamics (imgmov=10) with QTB',ch10,&
      'requires primitive coordinates (pitransform=0).'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  nkpt
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nkpt',nkpt,1,iout)
!  If prtdos>=2, nkpt must be greater or equal to 2
   if(dt%prtdos>=2)then
     cond_string(1)='prtdos' ; cond_values(1)=dt%prtdos
     call chkint_ge(1,1,cond_string,cond_values,ierr,'nkpt',nkpt,2,iout)
   end if
!  Must be smaller than 50 if iscf=-2 (band structure)
!  while prteig=0 and prtvol<2, except if kptopt>0
   if(dt%iscf==-2 .and. dt%prteig==0 .and. dt%prtvol<2 .and. dt%kptopt<=0)then
     cond_string(1)='iscf'   ; cond_values(1)=dt%iscf
     cond_string(2)='prteig' ; cond_values(2)=dt%prteig
     cond_string(3)='prtvol' ; cond_values(3)=dt%prtvol
     call chkint_le(1,3,cond_string,cond_values,ierr,'nkpt',nkpt,50,iout)
   end if

!  nloalg(1)= nloc_alg
!  Must be 2, 3, 4
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nloc_alg',dt%nloalg(1),3,(/2,3,4/),iout)

!  nloc_mem= nloalg(2)*(nloalg(3)+1)
!  nloalg(2) must be -1 or 1 ; nloalg(3) is 0 or 1.
   nloc_mem=dt%nloalg(2)*(dt%nloalg(3)+1)
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nloc_mem',nloc_mem,4,(/-2,-1,1,2/),iout)

!  npband
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npband',dt%npband,1,iout)

!  npfft
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npfft',dt%npfft,1,iout)
!  If usepaw==1 and pawmixdg==0, npfft must be equal to 1
   if(usepaw==1 .and. dt%pawmixdg==0)then
     cond_string(1)='usepaw  ' ; cond_values(1)=usepaw
     cond_string(2)='pawmixdg' ; cond_values(2)=dt%pawmixdg
     call chkint_eq(1,2,cond_string,cond_values,ierr,'npfft',dt%npfft,1,(/1/),iout)
   end if
#ifdef HAVE_OPENMP
   if (dt%wfoptalg==114) then
     if ( xomp_get_num_threads(.true.) > 1 .and. dt%npfft > 1 ) then
       write(msg,'(4a,i4,a,i4,a)') "When compilied with OpenMP, the FFT parallelization is not ",&
        "compatible with multiple threads.",ch10,"Please set npfft to 1 (currently npfft=",&
        dt%npfft, ") or export OMP_NUM_THREADS=1 (currently ",xomp_get_num_threads(.true.),")"
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
#endif

!  npimage
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npimage',dt%npimage,1,iout)
!  At present, parallelism over images is not coded ...
!  call chkint_eq(0,0,cond_string,cond_values,ierr,'npimage',dt%npimage,1,(/1/),iout)

!  npkpt
!  Must be greater or equal to 1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npkpt',dt%npkpt,1,iout)

!  nppert
   cond_string(1)='paral_rf' ; cond_values(1)=1
   call chkint_ge(1,1,cond_string,cond_values,ierr,'nppert',dt%nppert,1,iout)

!  nproc
   if (response==1.and.nsppol==2.and.nproc>1.and.modulo(nproc,2)>0) then
     write(msg,'(3a)' ) &
      'For DFPT parallel calculations on spin-polarized systems (nsppol=2),',ch10,&
      'the number of processors must be even !'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  nproj
!  If there is more than one projector for some angular momentum
!  channel of some pseudopotential
   do ilang=0,3
!    nprojmax(ilang)=maxval(pspheads(1:npsp)%nproj(ilang)) ! Likely problems with HP compiler
     nprojmax(ilang)=pspheads(1)%nproj(ilang)
     if(npsp>=2)then

       do ii=2,npsp
         nprojmax(ilang)=max(pspheads(ii)%nproj(ilang),nprojmax(ilang))
       end do
     end if
   end do

!  npspinor
!  Must be equal to 1 or 2
   call chkint_eq(0,0,cond_string,cond_values,ierr,'npspinor',dt%npspinor,2,(/1,2/),iout)
!  If nspinor==1, npspinor must be equal to 1
   if(dt%nspinor==1 )then
     cond_string(1)='nspinor' ; cond_values(1)=dt%nspinor
     call chkint_eq(0,1,cond_string,cond_values,ierr,'npspinor',dt%npspinor,1,(/1/),iout)
   end if

!  npimage

!  npvel (must be positive)
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npvel',dt%npvel,0,iout)

!  npwkss
!  Must be greater or equal to -1
   call chkint_ge(0,0,cond_string,cond_values,ierr,'npwkss',dt%npwkss,-1,iout)

!  nqpt
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nqpt',dt%nqpt,2,(/0,1/),iout)

!  nscforder
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nscforder',dt%nscforder,10,(/8,14,16,20,24,30,40,50,60,100/),iout)

!  nspden
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nspden',nspden,3,(/1,2,4/),iout)

   if(nsppol==2)then  !  When nsppol=2, nspden must be 2
     cond_string(1)='nsppol' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,1,(/2/),iout)
   end if
   if(nspden==2 .and. nsppol==1 .and. response==1)then
     write(msg,'(13a)')&
      'nspden==2 together with nsppol==1 is not allowed',ch10,&
      'for response function calculations.',ch10,&
      'For antiferromagnetic materials, use nspden==2 and nsppol=2.',ch10,&
      'In this case, Shubnikov symmetries will be used to decrease',ch10,&
      'the number of perturbations. In a future version, it will also be',ch10,&
      'used to decrease the number of spin components (to be coded).',ch10,&
      'Action: change nsppol to 1, or check nspden.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if(nspden==4.and.response==1)then
     write(msg,'(3a)')&
      'nspden==4 allowed in response formalism.',ch10,&
      'BUT Non collinear magnetism under development in perturbative treatment.'
     MSG_WARNING(msg)
   end if
!  TR symmetry not allowed for NC magnetism, in the present version
!  (to be investigated further)
   if (nspden==4.and.(dt%kptopt==1.or.dt%kptopt==2)) then
     write(msg, '(8a)' ) ch10,&
      'When non-collinear magnetism is activated (nspden=4),',ch10,&
      'time-reversal symmetry cannot be used in the present',ch10,&
      'state of the code (to be checked and validated).',ch10,&
      'Action: choose kptopt different from 1 or 2.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
!  When densfor_pred<0 or 3, nspden must be 1 or 2
   if(dt%densfor_pred<0.or.dt%densfor_pred==3)then
     cond_string(1)='densfor_pred' ; cond_values(1)=dt%densfor_pred
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if
!  When ionmov=4 and iscf>10, nspden must be 1 or 2
   if(dt%ionmov==4.and.dt%iscf>10)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     cond_string(1)='iscf' ; cond_values(1)=dt%iscf
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if
!  When iprcel>49, nspden must be 1 or 2
   if(mod(dt%iprcel,100)>49)then
     cond_string(1)='iprcel' ; cond_values(1)=dt%iprcel
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspden',nspden,2,(/1,2/),iout)
   end if
   if(mgga==1.and.nspden==4)then
     write(msg, '(3a)' )&
&     'The meta-GGA XC kernel is not yet implemented for non-colinear magnetism case',ch10, &
&     'Please use "nspden=1 or 2".'
     MSG_ERROR(msg)
   end if

!  nspinor
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nspinor',nspinor,2,(/1,2/),iout)
   if(nspden==2)then !  When nspden=2, nspinor must be 1
     cond_string(1)='nspden' ; cond_values(1)=2
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if

   if(nspden==4)then  !  When nspden=4, nspinor must be 2
     cond_string(1)='nspden' ; cond_values(1)=4
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/2/),iout)
   end if
!  When iscf=-1, nspinor must be 1
   if(dt%iscf==-1)then
     cond_string(1)='iscf' ; cond_values(1)=-1
!    Make sure that nsppol is 1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if
!  spin-orbit is not implemented for the strain perturbation
   if(dt%rfstrs/=0)then
     cond_string(1)='rfstrs' ; cond_values(1)=dt%rfstrs
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if
!  When usepawu=2 or -2, nspinor must be 1
   if(abs(dt%usepawu)==2)then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
!    Make sure that nspinor is 1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',nspinor,1,(/1/),iout)
   end if

!  nsppol
   call chkint_eq(0,0,cond_string,cond_values,ierr,'nsppol',nsppol,2,(/1,2/),iout)

!  nstep
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nstep',dt%nstep,0,iout)
   if(dt%nstep==0)then
!    nstep==0 computation of energy not yet implemented with Fock term, see m_energy.F90
     cond_string(1)='usefock' ; cond_values(1)=dt%usefock
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usefock',dt%usefock,1,(/0/),iout)
   endif

!  nsym
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nsym',dt%nsym,1,iout)
!  check if nsym=1 in phonon calculation in finite electric field
   if (response==1.and.&
&   (dt%berryopt== 4.or.dt%berryopt== 6.or.dt%berryopt== 7.or.&
&   dt%berryopt==14.or.dt%berryopt==16.or.dt%berryopt==17)) then
     cond_string(1)='response' ; cond_values(1)=1
     cond_string(2)='berryopt' ; cond_values(2)=4
     call chkint_eq(1,2,cond_string,cond_values,ierr,'nsym',dt%nsym,1,(/1/),iout)
   end if

!  ntime
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ntime',dt%ntime,0,iout)

!  ntimimage
   call chkint_ge(0,0,cond_string,cond_values,ierr,'ntimimage',dt%ntimimage,1,iout)

!  ntypalch
   if (usepaw==1) then
     cond_string(1)='pspcod' ; cond_values(1)=7;cond_values(2)=17
     call chkint_eq(1,2,cond_string,cond_values,ierr,'ntypalch',dt%ntypalch,1,(/0/),iout)
   end if

!  nucdipmom

   if (any(abs(dt%nucdipmom)>tol8)) then

!    nucdipmom requires PAW
     if(usepaw/=1)then
       write(msg, '(3a)' )&
        'Nuclear dipole moments (variable nucdipmom) input as nonzero but PAW not activated => stop',ch10,&
        'Action: re-run with PAW '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    nucdipmom requires complex rhoij
     if(dt%pawcpxocc/=2)then
       write(msg, '(3a)' )&
       'Nuclear dipole moments (variable nucdipmom) require complex rhoij => stop',ch10,&
       'Action: re-run with pawcpxocc = 2 '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    nucdipmom requires no force or stress calculation
     if(dt%optforces/=0 .OR. dt%optstress/=0)then
       write(msg, '(3a)' )&
       'Nuclear dipole moments (variable nucdipmom) cannot be used with force or stress calculations => stop',ch10,&
       'Action: re-run with optforces = 0 and optstress = 0 '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    nucdipmom requires kptopt > 2
     if(dt%kptopt<=2) then
       write(msg, '(a,i4,a,a,a)' )&
       ' Nuclear dipole moments (variable nucdipmom) break time reveral symmetry but kptopt = ',dt%kptopt,&
       ' => stop ',ch10,&
       'Action: re-run with kptopt greater than 2 '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

     ! nucdipmom is not currently compatible with spinat (this is necessary because both are used in symfind)
     if( any(abs(dt%spinat) > tol8) ) then
       write(msg, '(3a)' )&
        ' Nuclear dipole moments (variable nucdipmom) input as nonzero but spinat is also nonzero => stop',ch10,&
        'Action: re-run with spinat zero '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

   end if

!  nzchempot
   call chkint_ge(0,0,cond_string,cond_values,ierr,'nzchempot',dt%nzchempot,0,iout)
!  Cannot be used with response functions at present
   if (response==1) then
     cond_string(1)='response' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nzchempot',dt%nzchempot,1,(/0/),iout)
   end if
   if(dt%nzchempot>0)then
     do itypat=1,dt%ntypat
       do iz=2,dt%nzchempot
         dz=dt%chempot(1,iz,itypat)-dt%chempot(1,iz-1,itypat)
         if(dz<-tol12)then
           write(msg, '(a,2i6,a,a,d17.10,a,a, a,d17.10,a,a, a,a,a)' )&
            ' For izchempot,itypat=',iz,itypat,ch10,&
            ' chempot(1,izchempot-1,itypat) = ',dt%chempot(1,iz-1,itypat),' and', ch10,&
            ' chempot(1,izchempot  ,itypat) = ',dt%chempot(1,iz  ,itypat),',',ch10,&
            ' while they should be ordered in increasing values =>stop',ch10,&
            'Action: correct chempot(1,*,itypat) in input file.'
           MSG_ERROR_NOSTOP(msg, ierr)
         end if
       end do
       dz=dt%chempot(1,dt%nzchempot,itypat)-dt%chempot(1,1,itypat)
       if(dz>one)then
         write(msg, '(a,2i6,a,a,d17.10,a,a, a,d17.10,a,a, a,a,a)' )&
          ' For nzchempot,itypat=',dt%nzchempot,itypat,ch10,&
          ' chempot(1,1,itypat) = ',dt%chempot(1,1,itypat),' and', ch10,&
          ' chempot(1,nzchempot  ,itypat) = ',dt%chempot(1,dt%nzchempot,itypat),'.',ch10,&
          ' However, the latter should, at most, be one more than the former =>stop',ch10,&
          'Action: correct chempot(1,nzchempot,itypat) in input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end if

!  occ
!  Do following tests only for occopt==0 or 2, when occupation numbers are needed
   if ((dt%iscf>0.or.dt%iscf==-1.or.dt%iscf==-3) .and. (dt%occopt==0 .or. dt%occopt==2) ) then
     do iimage=1,dt%nimage
!      make sure occupation numbers (occ(n)) were defined:
       sumocc=zero
       bantot=0
       do isppol=1,nsppol
         do ikpt=1,nkpt
           do iband=1,dt%nband(ikpt+(isppol-1)*nkpt)
             bantot=bantot+1
             sumocc=sumocc+dt%occ_orig(bantot,iimage)
             if (dt%occ_orig(bantot,iimage)<-tol8) then
               write(msg, '(a,3i6,a,e20.10,a,a,a)' )&
                'iband,ikpt,iimage=',iband,ikpt,iimage,' has negative occ=',dt%occ_orig(bantot,iimage),' =>stop',ch10,&
                'Action: correct this occupation number in input file.'
               MSG_ERROR_NOSTOP(msg, ierr)
             end if
           end do
         end do
       end do
       if (sumocc<=1.0d-8) then
         write(msg, '(a,1p,e20.10,a,a,a)')&
          'Sum of occ=',sumocc, ' =>occ not defined => stop',ch10,&
          'Action: correct the array occ in input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     enddo
   end if

!  occopt
   call chkint_eq(0,0,cond_string,cond_values,ierr,'occopt',dt%occopt,9,(/0,1,2,3,4,5,6,7,8/),iout)
!  When prtdos==1 or 4, occopt must be between 3 and 8
   if(dt%prtdos==1.or.dt%prtdos==4)then
     write(cond_string(1), "(A)") 'prtdos'
     cond_values(1)=dt%prtdos
!    Make sure that occopt is 3,4,5,6,7, or 8
     call chkint_eq(1,1,cond_string,cond_values,ierr,'occopt',dt%occopt,6,(/3,4,5,6,7,8/),iout)
   end if
!  When nsppol==2 and spinmagntarget is the default value (-99.99d0), occopt cannot be 1.
   if(nsppol==2.and.dt%occopt==1.and.abs(dt%spinmagntarget+99.99d0)<tol8)then
     if(natom/=1 .or. abs(dt%znucl(dt%typat(1))-one)>tol8)then
       write(msg,'(a,i3,2a,i3,4a,f7.2,7a)' )&
        'This is a calculation with spin-up and spin-down wavefunctions,         ... nsppol=',nsppol,ch10,&
        'in which the occupation numbers are to be determined automatically.     ... occopt=',dt%occopt,ch10,&
        'However, in this case, the target total spin magnetization',ch10,&
        'must be specified, while the default value is observed.                 ... spinmagntarget=',dt%spinmagntarget,ch10,&
        'Action: if you are doing an antiferromagnetic calculation, please use nsppol=1 with nspden=2 ;',ch10,&
        'on the other hand, if you are doing a ferromagnetic calculation, either specify your own spinmagntarget,',ch10,&
        'or let the code determine the total spin-polarization, by using a metallic value for occopt (e.g. 7 or 4 ...).'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  optcell
   call chkint_eq(0,0,cond_string,cond_values,ierr,'optcell',dt%optcell,10,(/0,1,2,3,4,5,6,7,8,9/),iout)
!  With dt%berryopt=4, one must have optcell==0
!  if(dt%berryopt==4)then
!  cond_string(1)='berryopt' ; cond_values(1)=dt%berryopt
!  call chkint_eq(1,1,cond_string,cond_values,ierr,'optcell',dt%optcell,1,(/0/),iout)
!  end if

!  optdriver
   call chkint_eq(0,0,cond_string,cond_values,ierr,'optdriver',optdriver,9,&
&   [RUNL_GSTATE,RUNL_RESPFN,RUNL_SCREENING,RUNL_SIGMA,RUNL_NONLINEAR,RUNL_BSE, RUNL_GWLS, RUNL_WFK,RUNL_EPH],iout)
   if (response==1.and.all(dt%optdriver/=[RUNL_RESPFN,RUNL_NONLINEAR])) then
     write(msg,'(a,i3,3a,14(a,i2),4a)' )&
&     'The input variable optdriver=',dt%optdriver,ch10,&
&     'This is in conflict with the values of the other input variables,',ch10,&
&     'rfphon=',dt%rfphon,' rfddk=',dt%rfddk,' rf2_dkdk=',dt%rf2_dkdk,' rf2_dkde=',dt%rf2_dkde,&
&     ' rfelfd=',dt%rfelfd,'  rfmagn=',dt%rfmagn,' rfstrs=',dt%rfstrs,' rfuser=',dt%rfuser,&
&     ' d3e_pert1_elfd=',dt%d3e_pert1_elfd,' d3e_pert2_elfd=',dt%d3e_pert2_elfd,' d3e_pert3_elfd=',dt%d3e_pert3_elfd,&
&     ' d3e_pert1_phon=',dt%d3e_pert1_phon,' d3e_pert2_phon=',dt%d3e_pert2_phon,' d3e_pert3_phon=',dt%d3e_pert3_phon,ch10,&
&     'Action: check the values of optdriver, rfphon, rfddk, rf2dkdk, rf2dkde, rfelfd, rfmagn, rfstrs, rfuser,',ch10,&
&     'd3e_pert1_elfd, d3e_pert2_elfd, d3e_pert3_elfd, d3e_pert1_phon, d3e_pert2_phon, and d3e_pert3_phon in your input file.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if(usepaw==1)then
     ! Is optdriver compatible with PAW?
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(1,1,cond_string,cond_values,ierr,&
&     'optdriver',optdriver,7,[RUNL_GSTATE,RUNL_RESPFN,RUNL_SCREENING,RUNL_SIGMA,RUNL_BSE, RUNL_WFK,RUNL_NONLINEAR],iout)
   end if

!  Linear and Non-linear response calculations
   !Non-linear response not compatible with spinors
   if(nspinor/=1)then
     cond_string(1)='nspinor' ; cond_values(1)=nspinor
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   !Non-linear response only for insulators
   if(dt%occopt/=1 .and. dt%occopt/=2)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   !Non-linear response not compatible with mkmem=0
   if(dt%mkmem==0)then
     cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if
   !Response function need all k-points
   if(dt%kptopt==1 .or. dt%kptopt==4) then
     cond_string(1)='kptopt' ; cond_values(1)=dt%kptopt
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,2,(/RUNL_RESPFN,RUNL_NONLINEAR/),iout)
   end if
   !dkdk and dkde non-linear response only for occopt=1 (insulators)
   if (dt%rf2_dkdk==1) then
     cond_string(1)='rf2_dkdk' ; cond_values(1)=dt%rf2_dkdk
     call chkint_eq(1,1,cond_string,cond_values,ierr,'occopt',dt%occopt,1,(/1/),iout)
   end if
   if (dt%rf2_dkde==1) then
     cond_string(1)='rf2_dkde' ; cond_values(1)=dt%rf2_dkde
     call chkint_eq(1,1,cond_string,cond_values,ierr,'occopt',dt%occopt,1,(/1/),iout)
   end if
   !PEAD non-linear response only for occopt=1 (insulators)
   if(dt%usepead==0.and.dt%optdriver==RUNL_NONLINEAR)then
     cond_string(1)='usepead'   ; cond_values(1)=dt%usepead
     cond_string(2)='optdriver' ; cond_values(2)=dt%optdriver
     call chkint_eq(1,2,cond_string,cond_values,ierr,'occopt',dt%occopt,1,(/1/),iout)
   end if
   !PAW non-linear response only with DFPT (PEAD not allowed)
   if(usepaw==1.and.dt%optdriver==RUNL_NONLINEAR)then
     cond_string(1)='usepaw'    ; cond_values(1)=usepaw
     cond_string(2)='optdriver' ; cond_values(2)=dt%optdriver
     call chkint_eq(1,2,cond_string,cond_values,ierr,'usepead',dt%usepead,1,(/0/),iout)
     cond_string(1)='usepaw'    ; cond_values(1)=usepaw
     cond_string(2)='optdriver' ; cond_values(2)=dt%optdriver
     call chkint_eq(1,2,cond_string,cond_values,ierr,'pawxcdev',dt%pawxcdev,1,(/0/),iout)
   end if
   !Non-linear response not compatible with autoparal
   if(dt%optdriver==RUNL_NONLINEAR)then
     cond_string(1)='optdriver' ; cond_values(1)=dt%optdriver
     call chkint_eq(1,1,cond_string,cond_values,ierr,'autoparal',dt%autoparal,1,(/0/),iout)
   end if
   !Linear Response function only for LDA/GGA
   allow=(dt%ixc>0).and.((dt%ixc>16.and.dt%ixc<23).or.(dt%ixc>=40))
   if(.not.allow) allow=(dt%ixc<0).and.libxc_functionals_ismgga()
   if(allow)then
     cond_string(1)='ixc' ; cond_values(1)=dt%ixc
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_RESPFN/),iout)
   end if
   !PAW+Linear Response+GGA function restricted to pawxcdev=0
   if (dt%usepaw==1.and.dt%optdriver==RUNL_RESPFN) then
     allow=(dt%ixc>0).and.((dt%ixc>=11.and.dt%ixc<=16).or.(dt%ixc>=23.and.dt%ixc<=39))
     if(.not.allow) allow=(dt%ixc<0).and.libxc_functionals_isgga()
     if(allow) then
       if (dt%pawxcdev/=0)then
         write(msg,'(7a)' )&
&         'You are performing a DFPT+PAW calculation using a GGA XC functional:',ch10,&
&         '  This is restricted to pawxcdev=0!',ch10,&
&         '  Action: change pawxcdev value in your input file!',ch10,&
&         '    and be careful to run the preparatory Ground-State calculations also with pawxcdev=0!'
         MSG_ERROR_NOSTOP(msg, ierr)
       else
         write(msg,'(5a)' )&
&         'You are performing a DFPT+PAW calculation using a GGA XC functional:',ch10,&
&         '  - This is restricted to pawxcdev=0!',ch10,&
&         '  - Be careful to run the preparatory Ground-State calculations also with pawxcdev=0!'
         MSG_WARNING(msg)
       end if
     end if
   end if
   !Non linear Response function only for LDA (restricted to ixc=3/7/8)
   allow=(dt%ixc>0).and.(dt%ixc/=3.and.dt%ixc/=7.and.dt%ixc/=8)
   if(.not.allow) allow=(dt%ixc<0).and.(libxc_functionals_isgga().or.libxc_functionals_ismgga())
   if(allow)then
     cond_string(1)='ixc' ; cond_values(1)=dt%ixc
     call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',dt%optdriver,1,(/RUNL_NONLINEAR/),iout)
   end if

!  optforces
!  When ionmov>0, optforces must be >0
   if(dt%ionmov>0)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optforces',dt%optforces,2,(/1,2/),iout)
   end if
!  When imgmov>0, optforces must be >0
   if(dt%imgmov>0)then
     cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optforces',dt%optforces,2,(/1,2/),iout)
   end if
!  When iscf=22, optforces must be 0 or 2
   if(dt%iscf==22)then
     cond_string(1)='iscf' ; cond_values(1)=dt%iscf
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optforces',dt%optforces,2,(/0,2/),iout)
   end if

!  optstress
!  When optcell>0, optstress must be >0
   if(dt%optcell>0)then
     cond_string(1)='optcell' ; cond_values(1)=dt%optcell
     call chkint_eq(1,1,cond_string,cond_values,ierr,'optstress',dt%optstress,1,(/1/),iout)
   end if

  !  orbmag
  ! only values of 0 (default) 1, 2, 3 are allowed
  call chkint_eq(0,0,cond_string,cond_values,ierr,'orbmag',dt%orbmag,4,(/0,1,2,3/),iout)
  ! when orbmag /= 0, symmorphi must be 0 (no tnons)
  if(dt%orbmag .NE. 0) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'symmorphi',dt%symmorphi,1,(/0/),iout)
  end if
  ! only kptopt 4 and 3 are allowed for Chern number
  if((dt%orbmag.EQ.1).OR.(dt%orbmag.EQ.3)) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'kptopt',dt%kptopt,2,(/3,4/),iout)
  end if
  ! only kptopt 3 is allowed for orbmag
  if(dt%orbmag.GT.1) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'kptopt',dt%kptopt,1,(/3/),iout)
  end if
  ! only kpt parallelism is allowed at present
  if(dt%orbmag .GT. 0) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'paral_atom',dt%paral_atom,1,(/0/),iout)
  end if
  if(dt%orbmag .GT. 0) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,1,(/0/),iout)
  end if
  ! require usexcnhat 0
  if(dt%orbmag .NE. 0) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usexcnhat',dt%usexcnhat_orig,1,(/0/),iout)
  end if
  ! require PAW
  if(dt%orbmag .NE. 0) then
     cond_string(1)='orbmag';cond_values(1)=dt%orbmag
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',dt%usepaw,1,(/1/),iout)
  end if

!  paral_atom
   call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_atom',dt%paral_atom,2,(/0,1/),iout)
   if (dt%paral_atom/=0) then
     if (dt%optdriver/=RUNL_GSTATE.and.dt%optdriver/=RUNL_RESPFN) then
       write(msg, '(5a)' )&
&       'Parallelisation over atoms is only compatible with',ch10,&
&       'ground-state or response function calculations !',ch10,&
&       'Action: change paral_atom in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (dt%optdriver==RUNL_NONLINEAR) then
       cond_string(1)='optdriver' ; cond_values(1)=dt%optdriver
       call chkint_eq(1,1,cond_string,cond_values,ierr,'paral_atom',dt%paral_atom,1,(/0/),iout)
     end if
     if (dt%usedmft==1) then
       cond_string(1)='usedmft' ; cond_values(1)=dt%usedmft
       call chkint_eq(1,1,cond_string,cond_values,ierr,'paral_atom',dt%paral_atom,1,(/0/),iout)
     end if
     if (dt%prtden>1.and.dt%paral_kgb==0) then
       cond_string(1)='paral_kgb' ; cond_values(1)=dt%paral_kgb
       cond_string(2)='prtden' ; cond_values(2)=dt%prtden-1
       call chkint_eq(1,2,cond_string,cond_values,ierr,'paral_atom',dt%paral_atom,1,(/0/),iout)
     end if
   end if

!  paral_kgb
   call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,2,(/0,1/),iout)
!  Warning
   if(dt%paral_kgb==1.and.dt%iomode/=IO_MODE_MPI) then
     write(msg,'(11a)' )&
&     'When k-points/bands/FFT parallelism is activated',ch10,&
&     '(paral_kgb=1), only MPI-IO input/output is allowed !',ch10,&
&     'iomode/=1 in your input file',ch10,&
&     'You will not be able to perform input/output !'
     MSG_WARNING(msg)
   end if
   if(dt%paral_kgb==1.and.dt%nstep==0) then
     msg='When k-points/bands/FFT parallelism is activated, nstep=0 is not allowed!'
     MSG_ERROR_NOSTOP(msg,ierr)
   end if
   if(dt%paral_kgb==1.and.dt%usefock>0) then
     msg='Hartree-Fock or Hybrid Functionals are not compatible with bands/FFT parallelism!'
     MSG_ERROR_NOSTOP(msg,ierr)
   end if


!  paral_rf
   if ((response==0).and.(dt%paral_rf/=0)) then
     write(msg,'(a,i3,3a,14(a,i2),4a)' )&
&     'The input variable optdriver=',dt%optdriver,ch10,&
&     'This is in conflict with the values of the other input variables,',ch10,&
&     'rfphon=',dt%rfphon,' rfddk=',dt%rfddk,' rf2_dkdk=',dt%rf2_dkdk,' rf2_dkde=',dt%rf2_dkde,&
&     ' rfelfd=',dt%rfelfd,'  rfmagn=',dt%rfmagn,' rfstrs=',dt%rfstrs,' rfuser=',dt%rfuser,&
&     ' d3e_pert1_elfd=',dt%d3e_pert1_elfd,' d3e_pert2_elfd=',dt%d3e_pert2_elfd,' d3e_pert3_elfd=',dt%d3e_pert3_elfd,&
&     ' d3e_pert1_phon=',dt%d3e_pert1_phon,' d3e_pert2_phon=',dt%d3e_pert2_phon,' d3e_pert3_phon=',dt%d3e_pert3_phon,ch10,&
&     'Action: check the values of optdriver, rfphon, rfddk, rf2dkdk, rf2dkde, rfelfd, rfmagn, rfstrs, rfuser,',ch10,&
&     'd3e_pert1_elfd, d3e_pert2_elfd, d3e_pert3_elfd, d3e_pert1_phon, d3e_pert2_phon, and d3e_pert3_phon in your input file.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  pawcpxocc
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawcpxocc',dt%pawcpxocc,2,(/1,2/),iout)
     if (dt%usepawu/=0.and.nspinor==2.and.dt%pawcpxocc==1) then
       write(msg, '(5a)' )&
&       'When non-collinear magnetism is activated ,',ch10,&
&       'and LDA+U activated ',ch10,&
&       'PAW occupancies must be complex !'
       MSG_ERROR_NOSTOP(msg, ierr)
     else if (dt%pawspnorb==1.and.(dt%kptopt==0.or.dt%kptopt>=3).and.dt%pawcpxocc==1) then
       if (optdriver==RUNL_GSTATE.and.dt%iscf<10) then
         write(msg, '(11a)' )&
&         'When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         'and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         'PAW occupancies are complex !',ch10,&
&         'Their imaginary part is used to evaluate total energy by direct',ch10,&
&         'scheme, needed here because SCF potential mixing has been chosen (iscf<10).',ch10,&
&         'Action: put pawcpxocc=2 in input file, or choose SCF density mixing (iscf>=10).'
         MSG_ERROR_NOSTOP(msg, ierr)
       else if (optdriver==RUNL_GSTATE.and.dt%iscf>=10) then
         write(msg, '(11a)' )&
&         'When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         'and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         'PAW occupancies are complex !',ch10,&
&         'By setting pawcpxocc=1 in input file, their imaginary part',ch10,&
&         'is not computed. As a consequence, total energy computed',ch10,&
&         'is not available. Put pawcpxocc=2 in input file if you want it.'
         MSG_WARNING(msg)
       else
         write(msg, '(11a)' )&
&         'When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&         'and time-reversal symmetry is broken (kptopt/=1 and kptopt/=2)',ch10,&
&         'PAW occupancies are complex !',ch10,&
&         'Action: put pawcpxocc=2 in input file to compute their imaginary part.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end if
     if (dt%pawspnorb==1.and.dt%kptopt==0) then
       write(msg, '(7a)' )&
&       'When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
&       'time-reversal symmetry might be broken.',ch10,&
&       'Using kptopt=0 might be risky: if (kx,ky,kz) is present in k-points list,',ch10,&
&       '(-kx,-ky,-kz) (or equivalent) should also be present.'
       MSG_WARNING(msg)
     end if
   end if

!  pawcross
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawcross',dt%pawcross,2,(/0,1/),iout)
   end if

!  pawfatbnd
   call chkint_eq(0,0,cond_string,cond_values,ierr,'pawfatbnd',dt%pawfatbnd,3,(/0,1,2/),iout)
   if(usepaw/=1.and.dt%pawfatbnd>1) then
     MSG_ERROR_NOSTOP('pawfatbnd without PAW is not possible', ierr)
   end if
   if(dt%prtdosm==1.and.dt%pawfatbnd>0)then
     MSG_ERROR_NOSTOP('pawfatbnd>0  and prtdosm=1 are not compatible', ierr)
   end if
!  for the moment pawfatbnd is not compatible with fft or band parallelization
   !if (dt%pawfatbnd > 0 .and. (dt%npfft > 1 .or. dt%npband > 1)) then
   !  msg = 'pawfatbnd and FFT or band parallelization are not compatible yet. Set pawfatbnd to 0  '
   !  MSG_ERROR_NOSTOP(msg,ierr)
   !end if

!  pawlcutd
   if (usepaw==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'pawlcutd',dt%pawlcutd,0,iout)
   end if

!  pawlmix
   if (usepaw==1) then
     call chkint_ge(0,0,cond_string,cond_values,ierr,'pawlmix',dt%pawlmix,0,iout)
   end if

!  pawmixdg
   if (usepaw==1) then
     if(dt%ionmov==4)then
       cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawmixdg',dt%pawmixdg,1,(/1/),iout)
     end if
     if(dt%iscf==5.or.dt%iscf==6.or.dt%iscf==15.or.dt%iscf==16)then
       cond_string(1)='iscf' ; cond_values(1)=dt%iscf
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawmixdg',dt%pawmixdg,1,(/1/),iout)
     end if
     if(usewvl==1)then
       cond_string(1)='usewvl' ; cond_values(1)=1
       call chkint_eq(1,1,cond_string,cond_values,ierr,'pawmixdg',dt%pawmixdg,1,(/1/),iout)
     end if
   end if

!  pawnhatxc
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawnhatxc',dt%pawnhatxc,2,(/0,1/),iout)
   end if

!  pawnzlm
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawnzlm',dt%pawnzlm,2,(/0,1/),iout)
   end if

!  pawoptmix
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawoptmix',dt%pawoptmix,2,(/0,1/),iout)
   end if

!  pawprtdos
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawprtdos',dt%pawprtdos,3,(/0,1,2/),iout)
   end if

!  pawprtvol
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawprtvol',dt%pawprtvol,7,(/-3,-2,-1,0,1,2,3/),iout)
   end if

!  pawspnorb
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawspnorb',dt%pawspnorb,2,(/0,1/),iout)
     if (dt%pawspnorb==1.and.(dt%kptopt==1.or.dt%kptopt==2)) then
       write(msg, '(7a)' )&
        'When spin-orbit coupling is activated (pawspnorb=1),',ch10,&
        'time-reversal symmetry is broken; k-points cannot',ch10,&
        'be generated using TR-symmetry.',ch10,&
        'Action: choose kptopt different from 1 or 2.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  pawstgylm, pawsushat
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawstgylm',dt%pawstgylm,2,(/0,1/),iout)
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawsushat',dt%pawstgylm,2,(/0,1/),iout)
   end if

!  pawusecp
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,2,(/0,1/),iout)
!      if (dt%mkmem/=0)then
!        cond_string(1)='mkmem' ; cond_values(1)=dt%mkmem
!        call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
!      end if
!      if (dt%mk1mem/=0)then
!        cond_string(1)='mk1mem' ; cond_values(1)=dt%mk1mem
!        call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
!      end if
!      if (dt%mkqmem/=0)then
!        cond_string(1)='mkqmem' ; cond_values(1)=dt%mkqmem
!        call chkint_eq(1,1,cond_string,cond_values,ierr,'pawusecp',dt%pawusecp,1,(/1/),iout)
!      end if
   end if

!  pawxcdev
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'pawxcdev',dt%pawxcdev,3,(/0,1,2/),iout)
   end if

!  pimass
!  Check that masses are > 0 if imgmov = 9 or 13
   if (dt%imgmov==9.or.dt%imgmov==13) then
     do itypat=1,dt%ntypat
       cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
       write(input_name,'(a4,i1,a1)')'pimass(',itypat,')'
       call chkdpr(1,1,cond_string,cond_values,ierr,input_name,dt%pimass(itypat),1,tol8,iout)
     end do
   end if

!  pimd_constraint
   call chkint_eq(0,0,cond_string,cond_values,ierr,'pimd_constraint',dt%pimd_constraint,2,(/0,1/),iout)
   if(dt%pimd_constraint==1.and.dt%nconeq>1 )then
     cond_string(1)='pimd_constraint' ; cond_values(1)=dt%pimd_constraint
!    Make sure that nconeq=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nconeq',dt%nconeq,1,(/1/),iout)
   end if

!  pitransform
   call chkint_eq(0,0,cond_string,cond_values,ierr,'pitransform',dt%pitransform,3,(/0,1,2/),iout)
!  When imgmov is not one of 9 or 13, pitransform must be 0
   if(dt%imgmov/=9 .and. dt%imgmov/=13 )then
     cond_string(1)='imgmov' ; cond_values(1)=dt%imgmov
!    Make sure that pitransform=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'pitransform',dt%pitransform,1,(/0/),iout)
   end if
   if(dt%pimd_constraint/=0 )then
     cond_string(1)='pimd_constraint' ; cond_values(1)=dt%pimd_constraint
!    Make sure that pitransform=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'pitransform',dt%pitransform,1,(/0/),iout)
   end if

!  plowan_compute
   cond_string(1)='usepaw' ; cond_values(1)=1
   call chkint_eq(1,1,cond_string,cond_values,ierr,'plowan_compute',dt%plowan_compute,4,(/0,1,2,10/),iout)
   if(dt%plowan_compute>0) then
!    plowan_bandi/plowan_bandf
     !call chkint_ge(0,0,cond_string,cond_values,ierr,'plowan_bandi',dt%plowan_bandi,              1,iout)
     !call chkint_ge(0,0,cond_string,cond_values,ierr,'plowan_bandf',dt%plowan_bandf,dt%plowan_bandi,iout)

     !call chkint_le(0,0,cond_string,cond_values,ierr,'plowan_bandi',dt%plowan_bandi,dt%plowan_bandf,iout)
     !call chkint_le(0,0,cond_string,cond_values,ierr,'plowan_bandi',dt%plowan_bandf,dt%mband       ,iout)

     call chkint_ge(0,0,cond_string,cond_values,ierr,'plowan_natom',dt%plowan_natom,              0,iout)

     maxplowan_iatom=maxval(dt%plowan_iatom(1:dt%plowan_natom))
     minplowan_iatom=minval(dt%plowan_iatom(1:dt%plowan_natom))
     call chkint_ge(0,0,cond_string,cond_values,ierr,'plowan_iatom',minplowan_iatom,              1,iout)
     call chkint_le(0,0,cond_string,cond_values,ierr,'plowan_iatom',maxplowan_iatom,          natom,iout)

     kk=0
     do jj = 1, dt%plowan_natom
       do ii = 1, dt%plowan_nbl(jj)
         kk=kk+1
         cond_string(1)='usepaw' ; cond_values(1)=1
         call chkint_eq(1,1,cond_string,cond_values,ierr,'plowan_lcalc',dt%plowan_lcalc(kk),4,(/0,1,2,3/),iout)
       end do
     end do

     call chkint_ge(0,0,cond_string,cond_values,ierr,'plowan_nt'   ,dt%plowan_nt,                 0,iout)
     cond_string(1)='usepaw' ; cond_values(1)=1
     call chkint_eq(1,1,cond_string,cond_values,ierr,'plowan_realspace',dt%plowan_realspace,3,(/0,1,2/),iout)
   end if

!  posdoppler
   call chkint_eq(0,0,cond_string,cond_values,ierr,'posdoppler',dt%posdoppler,2,(/0,1/),iout)

!  positron
   call chkint_eq(0,0,cond_string,cond_values,ierr,'positron',dt%positron,7,(/-20,-10,-2,-1,0,1,2/),iout)
   if ((dt%positron==2.or.dt%positron<0).and.(dt%ixcpositron==3.or.dt%ixcpositron==31)) then
     if ((dt%ixc<11.or.dt%ixc>17).and.dt%ixc/=23.and.dt%ixc/=26.and.dt%ixc/=27) then
       write(msg, '(7a)' )&
&       'For the electronic ground-state calculation in presence of a positron,',ch10,&
&       'when GGA is selected for electron-positron correlation (ixcpositron=3 or 31),',ch10,&
&       'electron-electron XC must also be GGA !',ch10,&
&       'Action: choose another psp file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
   if (dt%positron/=0.and.dt%ionmov==5) then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'positron',dt%positron,1,(/0/),iout)
   end if
   if (dt%positron<0.and.usepaw==0) then
     write(msg, '(5a)' )&
&     'You cannot use positron<0 (automatic two-component DFT)',ch10,&
&     'with norm-conserving pseudopotentials !',ch10,&
&     'Action: choose PAW.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if
   if ((dt%positron==1.or.dt%positron<0).and.dt%iscf<10.and.dt%tolvrs>tiny(one)) then
     write(msg, '(7a)' )&
&     'You cannot perform a positronic ground-state calculation (positron=1 or <0)',ch10,&
&     'using SCF potential mixing (iscf<10) and tolvrs !',ch10,&
&     '(in that case, the potential is constant)',ch10,&
&     'Action: change iscf or select another convergence criterion.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  posocc
   call chkdpr(0,0,cond_string,cond_values,ierr,'posocc',dt%posocc,-1,one,iout)

!  postoldfe, postoldff
   call chkdpr(0,0,cond_string,cond_values,ierr,'postoldff',dt%postoldff,1,zero,iout)
   if (dt%positron<0) then
     if ( (abs(dt%postoldfe)> tiny(0.0_dp).and.abs(dt%postoldff)> tiny(0.0_dp)).or.&
&     (abs(dt%postoldfe)<=tiny(0.0_dp).and.abs(dt%postoldff)<=tiny(0.0_dp))) then
       write(msg,'(5a)' )&
&       'One and only one of the input tolerance criteria postoldfe or postoldff',ch10,&
&       'must differ from zero !',ch10,&
&       'Action: change postoldfe or postldff in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (abs(dt%postoldff)>tiny(0.0_dp).and.dt%optforces/=1)then
       write(msg,'(3a)' )&
&       'When postoldff is set to a non-zero value, optforces must be set to 1 !',ch10,&
&       'Action: change your input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if

!  prepanl
!  Must have prtden=1 to prepare a nonlinear calculation
   if (dt%prtden /= 1) then
     cond_string(1)='prtden' ; cond_values(1)=dt%prtden
     call chkint_ne(1,1,cond_string,cond_values,ierr,'prepanl',dt%prepanl,1,(/1/),iout)
   end if

!  prtbbb
!  Not allowed for PAW
   if(usepaw==1.and.dt%prtbbb==1)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtbbb',dt%prtbbb,1,(/0/),iout)
   end if

!  prtden
   if (usepaw==1) then
     call chkint_le(0,0,cond_string,cond_values,ierr,'prtden',dt%prtden,7,iout)
   else
     call chkint_le(0,0,cond_string,cond_values,ierr,'prtden',dt%prtden,1,iout)
   end if

!  prtdensph
   if (usepaw==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdensph',dt%prtdensph,2,(/0,1/),iout)
   end if

!  prtdos
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdos',dt%prtdos,6,(/0,1,2,3,4,5/),iout)

! for the moment prtdos 3,4,5 are not compatible with fft or band parallelization
   if (dt%prtdos > 3 .and. (dt%npfft > 1 .or. dt%npband > 1)) then
     MSG_ERROR_NOSTOP('prtdos>3 and FFT or band parallelization are not compatible yet. Set prtdos <= 2', ierr)
   end if

! prtdos 5 only makes sense for nspinor == 2. Otherwise reset to prtdos 2
   if (dt%prtdos == 5 .and. dt%nspinor /= 2) then
     dt%prtdos = 2
     MSG_WARNING('prtdos==5 is only useful for nspinor 2. Has been reset to 2')
   end if
   if (dt%prtdos == 5 .and. dt%npspinor /= 1) then
     MSG_ERROR_NOSTOP('prtdos==5 not available with npspinor==2', ierr)
   end if
   ! Consistency check for prtdos 5 with PAW
   if (dt%prtdos == 5 .and. dt%usepaw == 1) then
     if (dt%pawprtdos == 2) then
       MSG_ERROR_NOSTOP('prtdos==5 is not compatible with pawprtdos 2', ierr)
     end if
     MSG_ERROR_NOSTOP('prtdos==5 is not available with PAW', ierr)
   end if

!  prtdosm
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtdosm',dt%prtdosm,3,(/0,1,2/),iout)
   if(usepaw==1.and.dt%pawprtdos==1)then
     cond_string(1)='pawprtdos' ; cond_values(1)=dt%pawprtdos
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtdosm',dt%prtdosm,1,(/0/),iout)
   end if
   if(usepaw==1.and.dt%prtdosm>=1)then
     cond_string(1)='prtdosm' ; cond_values(1)=dt%prtdosm
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtdos',dt%prtdos,1,(/3/),iout)
   end if
   if(dt%prtdosm==2.and.dt%pawprtdos/=2)then
     MSG_ERROR(' pawprtdos/=2  and prtdosm=2 are not compatible')
   end if

!  prtefmas
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtefmas',dt%prtefmas,2,(/0,1/),iout)
   if(optdriver/=RUNL_RESPFN)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtefmas',dt%prtefmas,1,(/0/),iout)
   end if

!  prtelf
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtelf',dt%prtkden,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtelf',dt%prtelf,1,(/0/),iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtelf',dt%prtelf,1,(/0/),iout)
   end if

!  prtfsurf only one shift allowed (gamma)
   if (dt%prtfsurf == 1) then

     if (abs(dt%kptrlatt(1,2))+abs(dt%kptrlatt(1,3))+abs(dt%kptrlatt(2,3))+&
         abs(dt%kptrlatt(2,1))+abs(dt%kptrlatt(3,1))+abs(dt%kptrlatt(3,2)) /= 0 ) then
       write(msg,'(4a)')ch10,&
         'prtfsurf does not work with non-diagonal kptrlatt ', ch10,&
         'Action: set nshift 1 and shiftk 0 0 0'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (dt%nshiftk > 1) then
       write(msg,'(4a)') ch10,&
        'prtfsurf does not work with multiple kpt shifts ', ch10, &
        'Action: set nshift 1 and shiftk 0 0 0'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (sum(abs(dt%shiftk(:,1:dt%nshiftk))) > tol8) then
       write(msg,'(4a)')ch10,&
        'prtfsurf does not work with non-zero kpt shift ',ch10,&
        'Action: set nshift 1 and shiftk 0 0 0'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    Occcupations, Fermi level and k weights have to be calculated correctly.
     if (.not.(dt%iscf>1.or.dt%iscf==-3)) then
       write(msg,'(4a)')ch10,&
        'prtfsurf==1 requires either iscf>1 or iscf==-3 ',ch10,&
        'Action: change iscf in the input file. '
       MSG_ERROR_NOSTOP(msg, ierr)
     end if

!    Make sure all nband are equal (well it is always enforced for metals)
     if (any(dt%nband(1:nkpt*nsppol) /= maxval(dt%nband(1:nkpt*nsppol)) )) then
       write(msg,'(3a)')&
        'The number of bands has to be constant for the output of the Fermi surface.',ch10,&
        'Action: set all the nbands to the same value in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end if ! prtfsurf==1

!  prtgden
   call chkint(0,0,cond_string,cond_values,ierr,'prtgden',dt%prtgden,1,(/0/),1,0,iout)

   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint(0,1,cond_string,cond_values,ierr,'prtgden',dt%prtgden,1,(/0/),0,0,iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint(0,1,cond_string,cond_values,ierr,'prtgden',dt%prtgden,1,(/0/),0,0,iout)
   end if

!  prtkden
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtkden',dt%prtkden,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
   end if

!  prtlden
   call chkint(0,0,cond_string,cond_values,ierr,'prtlden',dt%prtlden,1,(/0/),1,0,iout)

   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint(0,1,cond_string,cond_values,ierr,'prtlden',dt%prtlden,1,(/0/),0,0,iout)
   end if
   if(usepaw/=0)then
     cond_string(1)='usepaw' ; cond_values(1)=usepaw
     call chkint(0,1,cond_string,cond_values,ierr,'prtlden',dt%prtlden,1,(/0/),0,0,iout)
   end if

!  prtstm
   call chkint_ge(0,0,cond_string,cond_values,ierr,'prtstm',dt%prtstm,0,iout)
   if(optdriver/=RUNL_GSTATE)then
     cond_string(1)='optdriver' ; cond_values(1)=optdriver
     call chkint_eq(0,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%occopt/=7)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%nstep/=1)then
     cond_string(1)='nstep' ; cond_values(1)=dt%nstep
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%ionmov/=0)then
     cond_string(1)='ionmov' ; cond_values(1)=dt%ionmov
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
!  tolwfr must be 0 to make a problem (another tol variable is used). Here, check that it is very very small.
   if(abs(dt%tolwfr)<tol16*tol16)then
     cond_string(1)='tolwfr' ; cond_values(1)=0
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtden/=0)then
     cond_string(1)='prtden' ; cond_values(1)=dt%prtden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtnabla>0)then
     cond_string(1)='prtnabla' ; cond_values(1)=dt%prtnabla
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvxc>0)then
     cond_string(1)='prtvxc' ; cond_values(1)=dt%prtvxc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvha>0)then
     cond_string(1)='prtvha' ; cond_values(1)=dt%prtvha
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if
   if(dt%prtvhxc>0)then
     cond_string(1)='prtvhxc' ; cond_values(1)=dt%prtvhxc
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtstm',dt%prtstm,1,(/0/),iout)
   end if

!  prtvclmb - needs prtvha as well
   if(dt%prtvclmb > 0)then
     cond_string(1)='prtvclmb' ; cond_values(1)=dt%prtvclmb
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtvha',dt%prtvha,1,(/1/),iout)
   end if

!  prtvolimg
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtvolimg',dt%prtvolimg,3,(/0,1,2/),iout)

!  prtwant
   if (dt%prtwant/=0) then
     cond_string(1)='prtwant' ; cond_values(1)=dt%prtwant
     call chkint_eq(0,0,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,1,(/0/),iout)
   end if
#if !defined HAVE_WANNIER90
   if(dt%prtwant==2) then
     write(msg, '(a,a,a)' )&
&     ' prtwant==2 is only relevant if wannier90 library is linked',ch10,&
&     ' Action: check compilation options'
     MSG_ERROR_NOSTOP(msg,ierr)
   end if
#endif

!  prtwf
   call chkint_eq(0,0,cond_string,cond_values,ierr,'prtwf',dt%prtwf,5,[-1,0,1,2,3],iout)

   if (dt%prtkbff == 1 .and. dt%useylm /= 0) then
     MSG_ERROR_NOSTOP("prtkbff == 1 requires useylm == 0", ierr)
   end if

!  random_atpos
   call chkint_eq(0,0,cond_string,cond_values,ierr,'random_atpos',dt%random_atpos,5,(/0,1,2,3,4/),iout)

!  ratsph
!  If PAW and (prtdos==3 or dt%prtdensph==1), must be greater than PAW radius
   if(usepaw==1.and.(dt%prtdos==3.or.dt%prtdensph==1))then
     do itypat=1,dt%ntypat
       if (pspheads(itypat)%pawheader%rpaw>dt%ratsph(itypat)) then
         write(msg, '(7a,i2,a,f15.12,3a)' )&
&         'Projected DOS/density is required in the framework of PAW !',ch10,&
&         'The radius of spheres in which DOS/density has to be projected',ch10,&
&         'must be greater or equal than the (max.) PAW radius !',ch10,&
&         'Rpaw(atom_type ',itypat,')= ',pspheads(itypat)%pawheader%rpaw,' au',ch10,&
&         'Action: modify value of ratsph in input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end if


!  recgratio
   if (dt%tfkinfunc==2) then
     write(msg, '(a,a)' ) ch10,&
&     '=== RECURSION METHOD ==========================================================='
     call wrtout(ab_out,msg,'COLL')
     cond_string(1)='tfkinfunc' ; cond_values(1)=2
     call chkint_ge(0,1,cond_string,cond_values,ierr,'recgratio',dt%recgratio,1,iout)
     if(dt%recgratio>1) then
       write(msg, '(a,a)' )&
&       '=== Coarse Grid is used in recursion ==========================================='
       call wrtout(ab_out,msg,'COLL')
       write(msg, '(a,i3,a,a,i3,a,i3,a,i3)' ) 'grid ratio =',dt%recgratio,&
&       ch10,'fine grid =   ',dt%ngfft(1),' ',dt%ngfft(2),' ',dt%ngfft(3)
       call wrtout(ab_out,msg,'COLL')
       write(msg, '(a,i2,a,i2,a,i2)' ) 'coarse grid = ',&
&       dt%ngfft(1)/dt%recgratio,' ',dt%ngfft(2)/dt%recgratio,' ',dt%ngfft(3)/dt%recgratio
       call wrtout(ab_out,msg,'COLL')
     else
       write(msg, '(a,i2,a,i2,a,i2)' ) 'fine grid =   ',dt%ngfft(1),' ',dt%ngfft(2),' ',dt%ngfft(3)
       call wrtout(ab_out,msg,'COLL')

     end if
   end if


!  rfatpol
   call chkint_ge(0,0,cond_string,cond_values,ierr,'rfatpol(1)',dt%rfatpol(1),1,iout)
   cond_string(1)='natom' ; cond_values(1)=natom
   call chkint_le(1,1,cond_string,cond_values,ierr,'rfatpol(2)',dt%rfatpol(2),natom,iout)

!  rfmeth
   call chkint_eq(0,0,cond_string,cond_values,ierr,'rfmeth',dt%rfmeth,6,(/-3,-2,-1,1,2,3/),iout)

!  rprimd
!  With optcell beyond 4, one has constraints on rprimd.
   if(dt%optcell==4 .or. dt%optcell==7 )then
     cond_string(1)='optcell' ; cond_values(1)=4
     if(dt%optcell==7)cond_values(1)=7
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
   else if(dt%optcell==5 .or. dt%optcell==8 )then
     cond_string(1)='optcell' ; cond_values(1)=5
     if(dt%optcell==8)cond_values(1)=8
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,1)',rprimd(2,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,2)',rprimd(1,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
   else if(dt%optcell==6 .or. dt%optcell==9 )then
     cond_string(1)='optcell' ; cond_values(1)=6
     if(dt%optcell==9)cond_values(1)=9
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,1)',rprimd(3,1),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(3,2)',rprimd(3,2),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(1,3)',rprimd(1,3),0,0.0_dp,iout)
     call chkdpr(1,1,cond_string,cond_values,ierr,'rprimd(2,3)',rprimd(2,3),0,0.0_dp,iout)
   end if

!  so_psp
   if(usepaw==0)then
     do ipsp=1,npsp
!      Check that so_psp is between 0 and 3
       if ( dt%so_psp(ipsp)<0 .or. dt%so_psp(ipsp)>3 ) then
         write(msg, '(a,i3,a,i3,a,a,a,a,a)' )&
&         'so_psp(',ipsp,' ) was input as ',dt%so_psp(ipsp),' .',ch10,&
&         'Input value must be 0, 1, 2, or 3.',ch10,&
&         'Action: modify value of so_psp (old name : so_typat) in input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
!      If nspinor=1, the spin-orbit contribution cannot be taken into account
       if ( nspinor==1 .and. (dt%so_psp(ipsp)==2 .or. dt%so_psp(ipsp)==3) ) then
         write(msg, '(a,i2,a,i3,a,a,a,a,a)' )&
&         'so_psp(',ipsp,') was input as ',dt%so_psp(ipsp),', with nspinor=1 and usepaw=0.',ch10,&
&         'When nspinor=1, so_psp cannot be required to be 2 or 3.',ch10,&
&         'Action: modify value of so_psp (old name : so_typat) or nspinor in input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
!      If nspinor=2, the spin-orbit contribution should be present in the pseudopotentials,
!      unless the user explicitly allows not to treat it.
       if ( nspinor==2 .and. dt%so_psp(ipsp)/=0 .and. pspheads(ipsp)%pspso==0 ) then
         write(msg, '(a,i2,a,i3,9a)' )&
&         'so_psp(',ipsp,') was input as ',dt%so_psp(ipsp),', with nspinor=2 and usepaw=0.',ch10,&
&         'This requires a treatment of the spin-orbit interaction. However, it has been detected ',ch10,&
&         'that the pseudopotential that you want to use does not specify the spin-orbit coupling.',ch10,&
&         'Action: choose a pseudopotential that contains information about the spin-orbit interaction,',ch10,&
&         ' or deliberately switch off the spin-orbit interaction by setting so_psp=0 for that pseudopotential in the input file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do ! ipsp
   end if ! usepaw==0

!  spinmagntarget
   if(abs(dt%spinmagntarget+99.99d0)>tol8 .and. abs(dt%spinmagntarget)>tol8)then
     if(nsppol==1)then
       write(msg, '(a,f8.2,4a)' )&
&       'spinmagntarget was input as ',dt%spinmagntarget,ch10,&
&       'When nsppol=1, spinmagntarget is required to be 0.0d0 or the default value.',ch10,&
&       'Action: modify value spinmagntarget or nsppol in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if(optdriver==RUNL_RESPFN)then
       write(msg, '(a,f8.2,4a)' )&
&       'spinmagntarget was input as ',dt%spinmagntarget,ch10,&
&       'For a response function run, spinmagntarget is required to be 0.0d0 or the default value.',ch10,&
&       'Action: modify value spinmagntarget or nsppol in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if(dt%prtdos==1)then
       write(msg, '(a,f8.2,4a)' )&
&       'spinmagntarget was input as ',dt%spinmagntarget,ch10,&
&       'When prtdos==1, spinmagntarget is required to be 0.0d0 or the default value.',ch10,&
&       'Action: modify value spinmagntarget or nsppol in input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if
!  If nsppol==2 and spinmagntarget==0.0 , suggest to use anti-ferromagnetic capability of ABINIT.
   if(nsppol==2.and.abs(dt%spinmagntarget)<tol8)then
     write(msg,'(a,i3,2a,f7.2,6a)' )&
&     ' This is a calculation with spin-up and spin-down wavefunctions,         ... nsppol=',nsppol,ch10,&
&     ' in which the target spin-polarization is zero.                  ... spinmagntarget=',dt%spinmagntarget,ch10,&
&     ' Tip ... It might be possible that the ground state is either non-spin-polarized, or antiferromagnetic.',ch10,&
&     ' In the former case, it is advantageous to use nsppol=1 and nspden=1,',ch10,&
&     ' while in the latter  case, it is advantageous to use nsppol=1 and nspden=2.'
     call wrtout(iout,msg,'COLL')
   end if

!  stmbias
   cond_string(1)='prtstm' ; cond_values(1)=dt%prtstm
   if(dt%prtstm/=0)then
!    If non-zero prtstm, stmbias cannot be zero : test is positive or zero
     if(dt%stmbias>-tol10)then
!      Then, enforce positive
       call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',dt%stmbias,1,2*tol10,iout)
     end if
   else
     call chkdpr(1,1,cond_string,cond_values,ierr,'stmbias',dt%stmbias,0,zero,iout)
   end if

!  string_algo
   call chkint_eq(0,0,cond_string,cond_values,ierr,'string_algo',dt%string_algo,2,(/1,2/),iout)

!  symafm
   if(nsppol==1 .and. nspden==2)then
!    At least one of the symmetry operations must be antiferromagnetic
     if(minval(dt%symafm(1:dt%nsym))/=-1)then
       write(msg, '(5a)' )&
&       'When nsppol==1 and nspden==2, at least one of the symmetry operations',ch10,&
&       'must be anti-ferromagnetic (symafm=-1), in order to deduce the spin-down density',ch10,&
&       'from the spin-up density.'
       call wrtout(iout,msg,'COLL')
       call wrtout(std_out,  msg,'COLL')
       write(msg, '(7a)' ) &
&       'However, it is observed that none of the symmetry operations is anti-ferromagnetic.',ch10,&
&       'Action: Check the atomic positions, the input variables spinat, symrel, tnons, symafm.',ch10,&
&       '        In case your system is not antiferromagnetic (it might be ferrimagnetic ...),',ch10,&
&       '        you must use nsppol=2 with nspden=2 (the latter being the default when nsppol=2).'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end if

!  symrel and tnons
!  Check the point group closure (TODO should check the spatial group closure !!)
   call chkgrp(dt%nsym,dt%symafm,dt%symrel,ierrgrp)
   if (ierrgrp==1) ierr=ierr+1

!  Check the orthogonality of the symmetry operations
!  (lengths and absolute values of scalar products should be preserved)
   iexit=0

   call chkorthsy(gprimd,iexit,dt%nsym,rmet,rprimd,dt%symrel)

!  symchi
   if (all(dt%symchi /= [0, 1])) then
     write(msg, '(a,i0,2a)' )'symchi was input as ',dt%symchi,ch10,'Input value must be 0, 1.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  symsigma
   if (all(dt%symsigma /= [0, 1, -1])) then
     write(msg, '(a,i0,a,a,a,a)' )&
      'symsigma was input as ',dt%symsigma,ch10,&
      'Input value must be 0, 1, or -1.',ch10,&
      'Action: modify value of symsigma in input file.'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  MG now it is possible to perform a GW calculation with non-symmorphic operations if required by the user
!  tnons
   if (dt%symmorphi==0) then
     if(dt%nbandkss/=0)then
       do isym=1,dt%nsym
         if(sum(dt%tnons(:,isym)**2)>tol6)then
           write(msg, '(3a,i3,a,3f8.4,3a)' )&
&           'When nbandkss/=0, all the components of tnons must be zero.',ch10,&
&           'However, for the symmetry operation number ',isym,', tnons =',dt%tnons(:,isym),'.',ch10,&
&           'Action: use the symmetry finder (nsym=0) with symmorphi==0.'
           MSG_ERROR_NOSTOP(msg,ierr)
         end if
       end do
     end if
     if (ANY(optdriver ==[RUNL_SCREENING,RUNL_SIGMA])) then
       do isym=1,dt%nsym
         if (sum(dt%tnons(:,isym)**2)>tol6) then
           write(msg,'(3a,i3,a,3f8.4,3a)')&
&           'When optdriver==RUNL_SCREENING or RUNL_SIGMA, all the components of tnons must be zero.',ch10,&
&           'However, for the symmetry operation number ',isym,', tnons =',dt%tnons(:,isym),'.',ch10,&
&           'Action: use the symmetry finder (nsym=0) with symmorphi==0.'
           MSG_ERROR_NOSTOP(msg, ierr)
         end if
       end do
     end if
   end if !of symmorphi

!  tfwkinfunc
   call chkint_eq(0,0,cond_string,cond_values,ierr,'tfwkinfunc',dt%tfkinfunc,5,(/0,1,2,11,12/),iout)
   if(dt%ionmov==4)then
     cond_string(1)='ionmov' ; cond_values(1)=4
     call chkint_eq(1,1,cond_string,cond_values,ierr,'tkinfunc',dt%tfkinfunc,1,(/0/),iout)
   end if

!  tolmxde
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolmxde',dt%tolmxde,1,zero,iout)

!  toldff
   call chkdpr(0,0,cond_string,cond_values,ierr,'toldff',dt%toldff,1,zero,iout)

!  tolimg
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolimg',dt%tolimg,1,zero,iout)

!  tolrde
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolrde',dt%tolrde,1,zero,iout)

!  tolrff
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolrff',dt%tolrff,1,zero,iout)

!  tolwfr
   call chkdpr(0,0,cond_string,cond_values,ierr,'tolwfr',dt%tolwfr,1,zero,iout)

!  tsmear
   call chkdpr(0,0,cond_string,cond_values,ierr,'tsmear',dt%tsmear,1,zero,iout)
!  Check that tsmear is non-zero positive for metallic occupation functions
   if(3<=dt%occopt .and. dt%occopt<=8)then
     cond_string(1)='occopt' ; cond_values(1)=dt%occopt
     call chkdpr(1,1,cond_string,cond_values,ierr,'tsmear',dt%tsmear,1,tol8,iout)
   end if

!  ucrpa
   call chkint_eq(0,0,cond_string,cond_values,ierr,'ucrpa',dt%ucrpa,5,(/0,1,2,3,4/),iout)
   if (dt%ucrpa>=1) then
     cond_string(1)='ucrpa' ; cond_values(1)=dt%ucrpa
     call chkint_eq(1,1,cond_string,cond_values,ierr,'nspinor',dt%nspinor,1,(/1/),iout)
   end if

!  usedmatpu
   if (usepaw==1.and.dt%usepawu>0) then
!    abs(dt%usedmatpu)<=nstep
     cond_string(1)='nstep' ; cond_values(1)=dt%nstep
     call chkint_le(1,1,cond_string,cond_values,ierr,'abs(usedmatpu)',abs(dt%usedmatpu),dt%nstep,iout)
!    lpawu must be constant or -1
     if (dt%usedmatpu/=0) then
       do itypat=1,dt%ntypat
         if (dt%lpawu(itypat)/=-1.and.dt%lpawu(itypat)/=maxval(dt%lpawu(:))) then
           write(msg, '(3a)' )&
&           'When usedmatpu/=0 (use of an initial density matrix for LDA+U),',ch10,&
&           'lpawu must be equal for all types of atoms on which +U is applied !'
           MSG_ERROR_NOSTOP(msg,ierr)
         end if
       end do
     end if
   end if
   !usedmatpu not allowed with experimental usepawu<0
   if (usepaw==1.and.dt%usepawu<0) then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usedmatpu',dt%usedmatpu,1,(/0/),iout)
   end if
   !usedmatpu only allowed for GS
   if (usepaw==1.and.dt%usedmatpu/=0.and.dt%optdriver/=RUNL_GSTATE) then
     MSG_ERROR_NOSTOP('usedmatpu/=0 is only allowed for Ground-State calculations!', ierr)
   end if

!  usedmft
   if (dt%usedmft>0) then
     cond_string(1)='usedmft' ; cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'usedmft',dt%usedmft,2,(/0,1/),iout)
     if (dt%paral_kgb>0) then
       cond_string(1)='usedmft' ; cond_values(1)=dt%usedmft
       call chkint_eq(1,1,cond_string,cond_values,ierr,'npspinor',dt%npspinor,1,(/1/),iout)
     end if
!    call chkint_eq(1,1,cond_string,cond_values,ierr,'paral_kgb',dt%paral_kgb,1,(/0/),iout)
   end if

!  useexexch and lexexch
!  Local exact-exchange and restrictions
   if(dt%useexexch/=0)then
     cond_string(1)='useexexch' ; cond_values(1)=dt%useexexch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'useexexch',dt%useexexch,1,(/1/),iout)
     cond_string(1)='useexexch' ; cond_values(1)=dt%useexexch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/1/),iout)
     cond_string(1)='useexexch' ; cond_values(1)=dt%useexexch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'pawxcdev',dt%pawxcdev,2,(/1,2/),iout)
     cond_string(1)='useexexch' ; cond_values(1)=dt%useexexch
     call chkint_eq(1,1,cond_string,cond_values,ierr,'ixc',dt%ixc,2,(/11,23/),iout)
     do itypat=1,dt%ntypat
       cond_string(1)='lexexch' ; cond_values(1)=dt%lexexch(itypat)
       call chkint_eq(1,1,cond_string,cond_values,ierr,'lexexch',dt%lexexch(itypat),5,(/-1,0,1,2,3/),iout)
     end do
   end if

!  usefock and restrictions
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usefock',dt%usefock,2,(/0,1/),iout)

!  usekden
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usekden',dt%usekden,2,(/0,1/),iout)
   if(dt%usekden==0)then
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'prtkden',dt%prtkden,1,(/0/),iout)
     if(mgga==1)then
       write(msg, '(3a)' )&
&       'The functional is a MGGA, but the kinetic energy density',ch10, &
&       'is not present. Please set "usekden 1" in the input file.'
       MSG_ERROR(msg)
     end if
   else if(dt%usekden/=0)then
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usewvl',usewvl,1,(/0/),iout)
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/0/),iout)
     cond_string(1)='usekden' ; cond_values(1)=dt%usekden
     call chkint_eq(1,1,cond_string,cond_values,ierr,'intxc',dt%intxc,1,(/0/),iout)
     do ipsp=1,npsp
!      Check that xccc is zero (metaGGAs cannot be used at present with non-linear core corrections)
       if ( pspheads(ipsp)%xccc/=0 ) then
         write(msg, '(3a,i0,3a)' )&
&         'When usekden is non-zero, it is not possible to use pseudopotentials with a non-linear core correction.',ch10,&
&         'However, for pseudopotential number ',ipsp,', there is such a core correction.',ch10,&
&         'Action: either set usekden=0 in input file, or change this pseudopotential file.'
         MSG_ERROR_NOSTOP(msg, ierr)
       end if
     end do
   end if

!  usepawu and lpawu
!  PAW+U and restrictions
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usepawu',dt%usepawu,11,(/-4,-3,-2,-1,0,1,2,3,4,10,14/),iout)
   if(dt%usepawu/=0)then
     cond_string(1)='usepawu' ; cond_values(1)=dt%usepawu
     call chkint_eq(1,1,cond_string,cond_values,ierr,'usepaw',usepaw,1,(/1/),iout)
     do itypat=1,dt%ntypat
       cond_string(1)='lpawu' ; cond_values(1)=dt%lpawu(itypat)
       call chkint_eq(1,1,cond_string,cond_values,ierr,'lpawu',dt%lpawu(itypat),5,(/-1,0,1,2,3/),iout)
     end do
     if(dt%pawspnorb>0) then
       write(msg,'(3a)' ) &
&       '  LDA+U+SpinOrbit is still on test ',ch10,&
&       '  (not yet in production)'
       MSG_WARNING(msg)
     end if
   end if

!  useexexch AND usepawu
!  Restriction when use together
   if(dt%useexexch/=0.and.dt%usepawu/=0)then
     do itypat=1,dt%ntypat
       if (dt%lpawu(itypat)/=dt%lexexch(itypat).and.&
&       dt%lpawu(itypat)/=-1.and.dt%lexexch(itypat)/=-1) then
         write(msg, '(5a,i2,3a)' )&
&         'When PAW+U (usepawu/=0) and local exact-exchange (useexexch/=0)',ch10,&
&         'are selected together, they must apply on the same',ch10,&
&         'angular momentum (lpawu/=lexexch forbidden, here for typat=',itypat,') !',ch10,&
&         'Action: correct your input file.'
         MSG_ERROR_NOSTOP(msg,ierr)
       end if
     end do
   end if

!  usedmft/usepawu and lpawu
!  Restriction when use together
   if(dt%usedmft>0.or.dt%usepawu/=0)then
     nlpawu=0
     do itypat=1,dt%ntypat
       if (dt%lpawu(itypat)/=-1) then
         nlpawu=nlpawu+1
       end if
     end do
     if(nlpawu==0) then
       write(msg, '(6a)' )&
&       'When DFT+U or DFT+DMFT is used',ch10,&
&       'at least one value of lpawu should be different from -1',ch10,&
&       'Action: correct your input file.'
       MSG_ERROR(msg)
     end if
   end if

!  usepotzero
   if(dt%usepotzero/=0)then
     if(dt%iscf<10) then
       write(msg, '(3a)' )&
&       'usepotzero can only be used with density mixing (not implemented yet)',ch10,&
&       'Action: choose iscf > 10 '
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end if

!  usexcnhat
   call chkint_eq(0,0,cond_string,cond_values,ierr,'usexcnhat',dt%usexcnhat_orig,3,(/-1,0,1/),iout)

!  useylm
   call chkint_eq(0,0,cond_string,cond_values,ierr,'useylm',dt%useylm,2,(/0,1/),iout)
   if (usepaw==1) then
     if(usewvl==0) then
       write(cond_string(1), "(A)") 'pspcod'
       cond_values(1)=7;cond_values(2)=17
       call chkint_eq(1,2,cond_string,cond_values,ierr,'useylm',dt%useylm,1,(/1/),iout)
     else
       write(cond_string(1), "(A)") 'usewvl'
       cond_values(1)=1
       call chkint_eq(1,1,cond_string,cond_values,ierr,'useylm',dt%useylm,1,(/0/),iout)
     end if
   end if

!  use_gpu_cuda
   call chkint_eq(0,0,cond_string,cond_values,ierr,'use_gpu_cuda',dt%use_gpu_cuda,2,(/0,1/),iout)
   if (dt%use_gpu_cuda==1) then
     if (dt%nspinor==2) then
       write(msg,'(3a)')&
&       'Use of GPU is not allowed when nspinor==2 !',ch10,&
&       'Action: impose use_gpu_cuda=0 in your input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
!    if (dt%optdriver==RUNL_GSTATE.and.mod(dt%wfoptalg,10)/=4) then
!    write(msg,'(6a)') ch10,&
!    &       ' chkinp : ERROR -',ch10,&
!    &       '  When GPU is in use (use_gpu_cuda=1), wfoptalg must be 4 or 14 !',ch10,&
!    &       '  Action: change wfoptalg in your input file.'
!    call wrtout(std_out,msg,'COLL')
!    ierr=ierr+1
!    end if
     if (dt%useylm==0) then
       write(msg,'(3a)')&
&       'Use of GPU is not allowed when useylm==0 !',ch10,&
&       'Action: impose uselym=1 in your input file.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (dt%tfkinfunc>0) then
       write(msg,'(5a)')&
&       'use_gpu_cuda=1 is not allowed when tfkinfunc>0 !',ch10,&
&       'Action: suppress use_gpu_cuda=0 from your input file',ch10,&
&       '        (GPU will be used but with another mechanism)'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     if (dt%ngfft(4)/=dt%ngfft(1).or.dt%ngfft(5)/=dt%ngfft(2).or.dt%ngfft(6)/=dt%ngfft(3)) then
       write(msg,'(3a)')&
&       'When GPU is in use (use_gpu_cuda=1), ngfft(4:6) must be equal to ngfft(1:3) !',ch10,&
&       'Action: suppress ngfft in input file or change it.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
#ifndef HAVE_GPU_CUDA
     write(msg,'(6a)') ch10,&
&     ' invars0: ERROR -',ch10,&
&     '   Input variables use_gpu_cuda is on but abinit hasn''t been built with gpu mode enabled !',ch10,&
&     '   Action: change the input variable use_gpu_cuda or re-compile ABINIT with Cuda enabled.'
     call wrtout(std_out,msg,'COLL')
     ierr=ierr+1
#endif
#ifndef HAVE_GPU_CUDA_DP
     write(msg,'(10a)') ch10,&
&     ' invars0: ERROR -',ch10,&
&     '   Input variables use_gpu_cuda is on but abinit hasn''t been built',ch10,&
&     '   with gpu mode in DOUBLE PRECISION enabled !',ch10,&
&     '   Action: change the input variable use_gpu_cuda',ch10,&
&     '   or re-compile ABINIT with double precision Cuda enabled.'
     call wrtout(std_out,msg,'COLL')
     ierr=ierr+1
#endif
   end if

!  use_slk
   if (dt%paral_kgb==1) then
     call chkint_eq(0,0,cond_string,cond_values,ierr,'use_slk',dt%use_slk,2,(/0,1/),iout)
   end if

!  use_yaml
   dt%use_yaml = 1
   call chkint_eq(0,0,cond_string,cond_values,ierr,'use_yaml',dt%use_yaml,2,(/0,1/),iout)

!  vdw_xc
   call chkint_eq(0,1,cond_string,cond_values,ierr,'vdw_xc',dt%vdw_xc,9,(/0,1,2,5,6,7,10,11,14/),iout)
   if (dt%usepaw==1.and.(.not.(dt%vdw_xc==0.or.dt%vdw_xc==5.or.dt%vdw_xc==6.or.dt%vdw_xc==7))) then
     write(msg,'(a,i2,a)')'vdw_xc=',dt%vdw_xc,' is not yet available with Projector Augmented-Wave (PAW) formalism!'
     MSG_ERROR_NOSTOP(msg, ierr)
   end if

!  vdw DFT-D2
   if (dt%vdw_xc==5.or.dt%vdw_xc==6.or.dt%vdw_xc==7) then
!    Only for GS or RF calculations
     if(optdriver/=RUNL_GSTATE.and.optdriver/=RUNL_RESPFN) then
       cond_string(1)='vdw_xc' ; cond_values(1)=dt%vdw_xc
       call chkint_ne(1,1,cond_string,cond_values,ierr,'optdriver',optdriver,2,(/RUNL_GSTATE,RUNL_RESPFN/),iout)
     end if
!    Restriction for DFT-D2
     if (dt%vdw_xc==5) then
!    Only with PBE, BP86  or BLYP GGA XC
       if(dt%ixc/=11.and.dt%ixc/=-101130.and.dt%ixc/=-130101.and. &
&       dt%ixc/=18.and.dt%ixc/=-106131.and.dt%ixc/=-131106.and. &
&       dt%ixc/=19.and.dt%ixc/=-106132.and.dt%ixc/=-132106.and. &
&       dt%ixc/=-202231.and.dt%ixc/=-231202) then
         write(msg,'(8a)') ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  Van der Waals DFT-D2 correction (vdw_xc=5) only available for the following XC functionals:',ch10,&
&         '      GGA-PBE, GGA-BLYP, GGA-BP86, mGGA-TPSS',ch10,&
&         '  Action: change your pseudopotential file.'
         call wrtout(std_out,msg,'COLL')
         ierr=ierr+1
       end if
!       Only for the first 5 lines of the periodic table
       do itypat=1,dt%ntypat
         if (dt%znucl(itypat)<0.or.dt%znucl(itypat)>54) then
           write(msg,'(4a,f5.1,a)') ch10,&
&           ' chkinp : ERROR -',ch10,&
&           '  Van der Waals DFT-D2 correction (vdw_xc=5) not available for atom type Z=',dt%znucl(itypat),' !'
           call wrtout(std_out,msg,'COLL')
           ierr=ierr+1
         end if
       end do
     end if
!    Restriction for DFT-D3/DFT-D3(BJ)
     if (dt%vdw_xc==6.or.dt%vdw_xc==7) then
!    Only with PBE, BP86  or BLYP GGA XC
       if(dt%ixc/=11.and.dt%ixc/=-101130.and.dt%ixc/=-130101.and. &
&       dt%ixc/=18.and.dt%ixc/=-106131.and.dt%ixc/=-131106.and. &
&       dt%ixc/=19.and.dt%ixc/=-106132.and.dt%ixc/=-132106.and. &
&       dt%ixc/=-202231.and.dt%ixc/=-231202.and.&
&       dt%ixc/=14.and.dt%ixc/=-102130.and.dt%ixc/=-130102.and. &
&       dt%ixc/=-170.and.dt%ixc/=41.and.dt%ixc/=-406) then
         write(msg,'(4a,i2,5a)') ch10,&
&         ' chkinp : ERROR -',ch10,&
&         '  Van der Waals DFT-D correction (vdw_xc=',dt%vdw_xc,') only available for the following XC functionals:',ch10,&
&         '      GGA-PBE, GGA-BLYP, GGA-BP86, mGGA-TPSS, GGA-RevPBE, PBE0',ch10,&
&         '  Action: change your pseudopotential file.'
         call wrtout(std_out,msg,'COLL')
         ierr=ierr+1
       end if
!       Only up to chemical species 96
       do itypat=1,dt%ntypat
         if (dt%znucl(itypat)<0.or.dt%znucl(itypat)>96) then
           write(msg,'(4a,i2,1a,f5.1,a)') ch10,&
&           ' chkinp : ERROR -',ch10,&
&           '  Van der Waals DFT-D correction (vdw_xc=',dt%vdw_xc,') not available for atom type Z=',dt%znucl(itypat),' !'
           call wrtout(std_out,msg,'COLL')
           ierr=ierr+1
         end if
       end do
     end if
   end if

!  wfoptalg
!  Must be greater or equal to 0
   call chkint_ge(0,0,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,0,iout)
!  wfoptalg==0,1,4,10,14 or 114 if PAW
   if (usepaw==1) then
     write(cond_string(1), "(A)") 'usepaw'
     cond_values(1)=1
     call chkint_eq(0,1,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,6,(/0,1,4,10,14,114/),iout)
   end if
!  wfoptalg/=114 if PAW+Fock
   if (usepaw==1 .and. dt%usefock==1) then
     write(cond_string(1), "(A)") 'usepaw'
     write(cond_string(2), "(A)") 'usefock'
     cond_values(1)=1
     cond_values(1)=2
     call chkint_ne(1,2,cond_string,cond_values,ierr,'wfoptalg',dt%wfoptalg,1,(/114/),iout)
   end if

   ! Check if FFT library supports MPI-FFT.
   if (dt%npfft > 1 .and..not. fftalg_has_mpi(fftalg)) then
     write(msg,"(a,i0,a)")"fftalg: ",fftalg," cannot be used in MPI-FFT mode (npfft > 1)"
     MSG_ERROR_NOSTOP(msg,ierr)
   end if

   ! Chebyshev
   if(dt%wfoptalg == 1) then
     if(dt%nspinor > 1) then
       msg='Nspinor > 1 not yet compatible with wfoptalg 1'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     !! TODO obsolete?
     if(dt%usefock > 0) then
       MSG_ERROR_NOSTOP('Fock not yet compatible with wfoptalg 1 (use Fock-level parallelism)', ierr)
     end if
     if(maxval(abs(dt%istwfk(1:nkpt))) > 2) then
       MSG_ERROR_NOSTOP('Istwfk > 2 not compatible with wfoptalg 1', ierr)
     end if
     if(dt%ecutsm > 0) then
       MSG_ERROR_NOSTOP('Ecutsm > 0 not yet compatible with wfoptalg 1', ierr)
     end if
     !! TODO check bandpp instead of overwriting it
   end if

!  np_slk
   call chkint_ge(0,0,cond_string,cond_values,ierr,'np_slk',dt%np_slk,0,iout)
   if (dt%np_slk>0 .and. dt%wfoptalg /= 114 ) then
     if(dt%np_slk <= dt%npfft*dt%npband*dt%npspinor .and. MOD(dt%npfft*dt%npband*dt%npspinor, dt%np_slk) /= 0) then
       MSG_ERROR_NOSTOP('np_slk must divide npfft*npband*npspinor.',ierr)
     end if
   end if
!  slk_rankpp
   call chkint_ge(0,0,cond_string,cond_values,ierr,'slk_rankpp',dt%slk_rankpp,0,iout)


!  wtk
!  Check that no k point weight is < 0:
   do ikpt=1,nkpt
     if (dt%wtk(ikpt)< -tiny(0.0_dp) ) then
       write(msg, '(a,i5,a,1p,e12.4,a,a,a)' )&
&       'At k point number',ikpt,'  wtk=',dt%wtk(ikpt),' <0.',ch10,&
&       'Action: check wtk in input file. Each wtk must be >=0.'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end do

!  xc_denpos
   call chkdpr(0,0,cond_string,cond_values,ierr,'xc_denpos',dt%xc_denpos,1,tiny(one),iout)

!  xc_tb09_c
   call chkdpr(0,0,cond_string,cond_values,ierr,'xc_tb09_c',dt%xc_tb09_c,1,0.0_dp,iout)

!  xred
!  Check that two atoms are not on top of each other
   do iimage=1,dt%nimage
     if(natom>1)then
       ABI_ALLOCATE(frac,(3,natom))
       do ia=1,natom
!        Map reduced coordinate xred(mu,ia) into [0,1)
         frac(1,ia)=dt%xred_orig(1,ia,iimage)-aint(dt%xred_orig(1,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(1,ia,iimage))
         frac(2,ia)=dt%xred_orig(2,ia,iimage)-aint(dt%xred_orig(2,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(2,ia,iimage))
         frac(3,ia)=dt%xred_orig(3,ia,iimage)-aint(dt%xred_orig(3,ia,iimage))+0.5_dp-sign(0.5_dp,dt%xred_orig(3,ia,iimage))
       end do
       do ia=1,natom-1
         do ib=ia+1,natom
           if( abs(frac(1,ia)-frac(1,ib))<1.0d-6 .and. &
&           abs(frac(2,ia)-frac(2,ib))<1.0d-6 .and. &
&           abs(frac(3,ia)-frac(3,ib))<1.0d-6         ) then
             if(iimage>1)then
               write(msg,'(2a,i5)') ch10,' The following was observed for image=',iimage
               call wrtout(iout,msg,'COLL')
               call wrtout(std_out,msg,'COLL')
             end if
             write(msg, '(a,i4,a,i4,a,a,a,a,a,a)' )&
&             'Atoms number',ia,' and',ib,' are located at the same point',' of the unit cell',ch10,&
&             '(periodic images are taken into account).',ch10,&
&             'Action: change the coordinate of one of these atoms in the input file.'
             MSG_ERROR_NOSTOP(msg,ierr)
           end if
         end do
       end do
       ABI_DEALLOCATE(frac)
     end if
   end do

!  znucl
!  Check that znucl and znuclpsp agree
   do ipsp=1,npsp
     if (abs(dt%znucl(ipsp)-pspheads(ipsp)%znuclpsp)> tol12 ) then
       write(msg, '(a,i0,a,es12.4,a,a,es12.4,2a)' )&
&       'For pseudopotential ',ipsp,'  znucl from user input file= ',dt%znucl(ipsp),ch10,&
&       'while znucl from pseudopotential file=',pspheads(ipsp)%znuclpsp,ch10,&
&       'Action: check znucl in input file, or check psp file. They must agree.'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end do

!  bandFFT
   if(dt%paral_kgb==1.and.dt%optdriver==RUNL_GSTATE) then
     if (mod(dt%wfoptalg,10) /= 4 .and. dt%wfoptalg /= 1) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of wfoptalg is found to be ',dt%wfoptalg,ch10,&
&       'This is not allowed in the case of band-FFT parallelization.',ch10,&
&       'Action: put wfoptalg = 4, 14 or 114 in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
!    Make sure all nband are equal
     if (any(dt%nband(1:nkpt*nsppol) /= maxval(dt%nband(1:nkpt*nsppol)) )) then
       write(msg,'(a,a,a)')&
&       'The number of bands have to remain constant in the case of band-FFT parallelization.',ch10,&
&       'Action: set all the nbands to the same value in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if(maxval(abs(dt%istwfk(1:nkpt)-1)) > 1)then
       write(msg,'(5a)' )&
&       'One of the components of istwfk is not equal to 1 or 2.',ch10,&
&       'Time-reversal symmetry is not yet programmed in the case of band-FFT parallelization.',ch10,&
&       'Action: set istwfk to 1 or 2 for all k-points'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if (dt%mkmem == 0) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of mkmem is found to be ',dt%mkmem,ch10,&
&       'An out-of-core solution can''t be used in the case of band-FFT parallelization.',ch10,&
&       'Action: put mkmem = nkpt in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end if

!  WVL - wavelets checks and limitations
   if(dt%usewvl == 1) then
     if (dt%wvl_hgrid <= 0) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of wvl_hgrid is found to be ',dt%wvl_hgrid,ch10,&
&       'This value is mandatory and must be positive.',ch10,&
&       'Action: put wvl_hgrid to a positive value in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if (dt%nsym /= 1 .and. dt%icoulomb == 1) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of nsym is found to be ',dt%nsym,ch10,&
&       'No symetry operations are allowed for isolated systems.',ch10,&
&       'Action: put nsym = 1 in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if (dt%optstress > 0) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of optstress is found to be ', dt%optstress, ch10,&
&       'There is no stress computation available with the wavelet code.',ch10,&
&       'Action: put optstress = 0 in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if (usepaw == 1) then
       MSG_WARNING('WVL+PAW is under development')
     end if
     if (dt%nspden > 2) then
       write(msg,'(a,i0,a,a,a,a)')&
&       'The value of nspden is found to be ', dt%nspden, ch10, &
&       'The wavelet computation is not allowed with non-colinear spin.',ch10,&
&       'Action: put nspden = 1 or 2 in your input file'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
     if (dt%nspden /= dt%nsppol) then
       write(msg,'(a,i0,a,a,i0,a,a)')&
&       'The value of nspden is found to be ', dt%nspden, ch10, &
&       'and the one of nsppol is found to be ', dt%nsppol, ch10, &
&       'In wavelet computation, nspden and nsppol must be equal.'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
!    We check the consistency of occupation, empty bands are not allowed.
     if (dt%nsppol == 2) then
       mband = dt%nelect
     else
       mband = dt%mband
     end if
     do iimage=1,dt%nimage
       do ii = 1, mband, 1
         if (dt%occ_orig(ii,iimage) < tol8 .and. dt%iscf == 0) then
           write(msg,'(a,f7.4,a,a,a,a,a,a)')&
&           'One value of occ is found to be ', dt%occ_orig(ii,iimage), ch10, &
&           'The direct minimization is not allowed with empty bands.',ch10,&
&           'Action: use occopt = 1 for automatic band filling or', ch10, &
&           'change occ value in your input file'
           MSG_ERROR_NOSTOP(msg,ierr)
         end if
       end do
     enddo
     if (npsp /= dt%ntypat) then
       write(msg, '(a,a,a,a,I0,a,I0,a,a,a)' ) ch10,&
&       'wvl_wfs_set:  consistency checks failed,', ch10, &
&       'dtset%npsp (', npsp, ') /= dtset%ntypat (', dt%ntypat, ').', ch10, &
&       'No alchemy pseudo are allowed with wavelets.'
       MSG_ERROR_NOSTOP(msg,ierr)
     end if
   end if

   ! Test on tolerances (similar tests are performed in scprqt, so keep the two versions in synch)
   if (any(optdriver == [RUNL_GSTATE, RUNL_RESPFN])) then
     ttolwfr=0 ; ttoldff=0 ; ttoldfe=0 ; ttolvrs=0; ttolrff=0
     if(abs(dt%tolwfr)>tiny(zero))ttolwfr=1
     if(abs(dt%toldff)>tiny(zero))ttoldff=1
     if(abs(dt%tolrff)>tiny(zero))ttolrff=1
     if(abs(dt%toldfe)>tiny(zero))ttoldfe=1
     if(abs(dt%tolvrs)>tiny(zero))ttolvrs=1

     ! If non-scf calculations, tolwfr must be defined
     if(ttolwfr /= 1 .and. ((dt%iscf<0 .and. dt%iscf/=-3) .or. dt%rf2_dkdk/=0 .or. dt%rf2_dkde/=0))then
       write(msg,'(a,a,a,es14.6,a,a)')&
        'when iscf <0 and /= -3, or when rf2_dkdk/=0 or rf2_dkde/=0, tolwfr must be strictly',ch10,&
        'positive, while it is ',dt%tolwfr,ch10,&
        'Action: change tolwfr in your input file and resubmit the job.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
     !  toldff only allowed when prtfor==1
     !if((ttoldff == 1 .or. ttolrff == 1) .and. dt%prtfor==0 )then
     !  MSG_ERROR_NOSTOP('toldff only allowed when prtfor=1!', ierr)
     !end if

     ! If SCF calculations, one and only one of these can differ from zero
     ! FIXME: this test should be done on input, not during calculation
     if(ttolwfr+ttoldff+ttoldfe+ttolvrs+ttolrff /= 1 .and. (dt%iscf>0 .or. dt%iscf==-3))then
       write(msg,'(6a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,a,es14.6,a,a,a)' )&
&       'For the SCF case, one and only one of the input tolerance criteria ',ch10,&
&       'tolwfr, toldff, tolrff, toldfe or tolvrs ','must differ from zero, while they are',ch10,&
&       'tolwfr=',dt%tolwfr,', toldff=',dt%toldff,', tolrff=',dt%tolrff,', toldfe=',dt%toldfe,ch10,&
&       'and tolvrs=',dt%tolvrs,' .',ch10,&
&       'Action: change your input file and resubmit the job.'
       MSG_ERROR_NOSTOP(msg, ierr)
     end if
   end if


!  If molecular dynamics or structural optimization is being done
!  (dt%ionmov>0), make sure not all atoms are fixed
!  if (dt%ionmov > 0) then
!  if (natfix == natom) then
!  write(msg, '(a,a,a,a,i4,a,i5,a,a,i5,a,a,a,a,a,a)' ) ch10,&
!  &   ' setup1: ERROR -',ch10,&
!  &   '  ionmov is ',dt%ionmov,' and number of fixed atoms is ',natfix,ch10,&
!  &   '  while number of atoms natom is ',natom,'.',ch10,&
!  &   '  Thus all atoms are fixed and option ionmov to move atoms',&
!  &           ' is inconsistent.',ch10,&
!  &   '  Action: change ionmov or natfix and iatfix in input file and resubmit.'
!  call wrtout(std_out,msg,"COLL")
!  ierr = ierr + 1
!  end if
!  end if

!  Should check that the symmetry operations are consistent with iatfixx,
!  iatfixy and iatfixz (diagonal symmetry operations)

!  Should check values of fftalg

!  rfasr=2 possible only when electric field response is computed.

!  Must have nqpt=1 for rfphon=1

!  ** Here ends the checking section **************************************

   call dt%free()
   ierr_dtset(idtset)=ierr

 end do !  End do loop on idtset

 if (maxval(dtsets(:)%usewvl) > 0) then
   write(msg,'(4a)') ch10,&
&   ' Comparison between wvl_hgrid and ecut',ch10,&
&   '  real-space mesh | eq. Ec around atoms | eq. Ec further from atoms'
   MSG_COMMENT(msg)
   wvl_hgrid = zero
   twvl = .false.
   do idtset=1,ndtset_alloc
!    Give an indication to the equivalent ecut corresponding to
!    given hgrid.
     if (dtsets(idtset)%usewvl == 1 .and. &
&     wvl_hgrid /= dtsets(idtset)%wvl_hgrid) then
       write(msg,'(F11.3,A,F16.1,A,F16.1,A)') &
&       dtsets(idtset)%wvl_hgrid, " bohr  |", &
&       two * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht  | ", &
&       half * pi * pi / (dtsets(idtset)%wvl_hgrid ** 2), " Ht"
       call wrtout(std_out,msg,'COLL')
       wvl_hgrid = dtsets(idtset)%wvl_hgrid
     end if
     twvl = twvl .or. (dtsets(idtset)%usewvl == 1 .and. dtsets(idtset)%iomode /= IO_MODE_ETSF)
   end do
   if (twvl) then
     write(msg,'(5a)') &
&     'Restart files from wavelets in non ETSF format does not follow', ch10, &
&     'the ABINIT standards.', ch10, &
&     'Put iomode to 3 to use ETSF retart files.'
     MSG_WARNING(msg)
   end if
 end if

 ! If there was a problem, then stop.
 call xmpi_sum(ierr_dtset, comm, mpierr)
 ierr=sum(ierr_dtset(1:ndtset_alloc)/mpi_enregs(1:ndtset_alloc)%nproc)

 if (ierr==1) then
   write(msg,'(a,i0,3a)')&
   'Checking consistency of input data against itself gave ',ierr,' inconsistency.',ch10,&
   'The details of the problem can be found above.'
   MSG_ERROR(msg)
 end if
 if (ierr>1) then
   write(msg,'(a,i0,3a)')&
   'Checking consistency of input data against itself gave ',ierr,' inconsistencies.',ch10,&
   'The details of the problems can be found above.'
   MSG_ERROR(msg)
 end if

 ABI_DEALLOCATE(ierr_dtset)

 if (ndtset_alloc /= 1 .and. get_timelimit() > zero) then
   MSG_ERROR("--timelimit option cannot be used when ndtset > 1")
 end if

 DBG_EXIT("COLL")

end subroutine chkinp
!!***

end module m_chkinp
!!***
