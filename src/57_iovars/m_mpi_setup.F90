!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mpi_setup
!! NAME
!! m_mpi_setup
!!
!! FUNCTION
!!  Initialize MPI parameters and datastructures for parallel execution
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (FJ, MT, FD)
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

module m_mpi_setup

 use defs_basis
 use defs_abitypes
 use m_distribfft
 use m_xmpi
 use m_xomp
 use m_hdr
 use m_sort
 use m_errors
 use m_abicore

 use m_time,         only : abi_wtime
 use m_io_tools,     only : flush_unit
 use m_parser,       only : intagm
 use m_geometry,     only : mkrdim, metric
 use m_fftcore,      only : fftalg_for_npfft, getng,  kpgcount
 use m_mpinfo,       only : init_mpi_enreg, mpi_distrib_is_ok, initmpi_atom, proc_distrb_cycle, &
                            initmpi_grid, initmpi_pert, initmpi_img, distrb2, distrb2_hf, initmpi_world
 use m_libpaw_tools, only : libpaw_write_comm_set
 use m_dtset,        only : get_npert_rbz
 use m_kg,           only : getmpw
 use m_dtfil,        only : mkfilename

 implicit none

 private
!!***

 public :: mpi_setup
!!***

contains
!!***

!!****f* ABINIT/mpi_setup
!! NAME
!! mpi_setup
!!
!! FUNCTION
!! Big loop on the datasets:
!! - compute mgfft,mpw,nfft,... for this data set;
!! - fill mpi_enreg
!!  *** At the output of this routine, all the dtsets input variables are known ***
!! The content of dtsets should not be modified anymore afterwards.
!!
!! INPUTS
!!  filnam(5)=character strings giving file names
!!  ndtset= number of datasets to be read; if 0, no multi-dataset mode
!!  ndtset_alloc=number of datasets, corrected for allocation of at least
!!      one data set.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are initialized here, while other were already
!!   initialized previously.
!!
!! SIDE EFFECTS
!!   mpi_enregs=information about MPI parallelization
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      abi_io_redirect,distrb2,distrb2_hf,finddistrproc,get_npert_rbz,getmpw
!!      getng,init_distribfft,init_mpi_enreg,initmpi_atom,initmpi_grid
!!      initmpi_img,initmpi_pert,intagm,libpaw_write_comm_set,metric,mkrdim
!!      wrtout
!!
!! SOURCE

subroutine mpi_setup(dtsets,filnam,lenstr,mpi_enregs,ndtset,ndtset_alloc,string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,ndtset,ndtset_alloc
 type(MPI_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
 character(len=*),intent(in) :: string
!arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------------
!scalars
 integer :: blocksize,exchn2n3d,iband,idtset,iexit,ii,iikpt,iikpt_modulo, prtvol
 integer :: isppol,jdtset,marr,mband_lower,mband_upper
 integer :: me_fft,mgfft,mgfftdg,mkmem,mpw,mpw_k,optdriver
 integer :: nfft,nfftdg,nkpt,nkpt_me,npert,nproc,nproc_fft,nqpt
 integer :: nspink,nsppol,nsym,paral_fft,response,tnband,tread0,usepaw,vectsize
 integer :: fftalg,fftalga,fftalgc
#ifdef HAVE_LINALG_ELPA
 integer :: icol,irow,np
#endif
 logical :: fftalg_read,ortalg_read,wfoptalg_read,do_check
 real(dp) :: dilatmx,ecut,ecut_eff,ecutdg_eff,ucvol
 character(len=500) :: message
!arrays
 integer :: ngfft(18),ngfftdg(18),ngfftc(3),tread(12)
 integer,allocatable :: intarr(:),istwfk(:),symrel(:,:,:)
 integer,pointer :: nkpt_rbz(:)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: gmet(3,3),gprimd(3,3),kpt(3),qphon(3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: dprarr(:),kpt_with_shift(:,:)
 real(dp),pointer :: nband_rbz(:,:)
 character(len=6) :: nm_mkmem(3)

!*************************************************************************

 DBG_ENTER("COLL")

 iexit=0;mpw_k=0

 call init_mpi_enreg(mpi_enregs(0))
 call initmpi_img(dtsets(0),mpi_enregs(0),-1)

 do idtset=1,ndtset_alloc
   call init_mpi_enreg(mpi_enregs(idtset))

   ! Handy read-only variables.
   optdriver = dtsets(idtset)%optdriver
   prtvol = dtsets(idtset)%prtvol

!  Read parallel input parameters
   marr=max(5,dtsets(idtset)%npsp,dtsets(idtset)%nimage)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
   nkpt  =dtsets(idtset)%nkpt
   nsppol=dtsets(idtset)%nsppol
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   usepaw=dtsets(idtset)%usepaw
   mband_upper=maxval(dtsets(idtset)%nband(1:nkpt*nsppol))
   mband_lower=minval(dtsets(idtset)%nband(1:nkpt*nsppol))

!  Compute metric for this dataset
   call mkrdim(dtsets(idtset)%acell_orig(1:3,1),dtsets(idtset)%rprim_orig(1:3,1:3,1),rprimd)
   call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

   ! Read paral_kgb and disable it if not supported in optdriver.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'paral_kgb',tread(1),'INT')
   if (tread(1)==1) dtsets(idtset)%paral_kgb=intarr(1)

   if(xmpi_paral==0.and.dtsets(idtset)%paral_kgb==1)then
     dtsets(idtset)%paral_kgb=0
     write(message, '(5a)' ) &
&     'When ABINIT is compiled without MPI flag,',ch10,&
&     'setting paral_kgb/=0 is useless. paral_kgb has been reset to 0.',ch10,&
&     'Action: modify compilation option or paral_kgb in the input file.'
     MSG_WARNING(message)
   end if

   if ( ALL(optdriver /= [RUNL_GSTATE, RUNL_GWLS]) .and. dtsets(idtset)%paral_kgb/=0) then
     dtsets(idtset)%paral_kgb=0
     write(message, '(a,i0,a)') &
&     "paral_kgb != 0 is not available in optdriver ",optdriver,". Setting paral_kgb to 0"
     MSG_COMMENT(message)
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'max_ncpus',tread0,'INT')
   if (tread0==1) dtsets(idtset)%max_ncpus=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'paral_atom',tread0,'INT')
   if(tread0==1) dtsets(idtset)%paral_atom=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'paral_rf',tread0,'INT')
   if (tread0==1.and.any(optdriver==[RUNL_RESPFN, RUNL_NONLINEAR])) dtsets(idtset)%paral_rf=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npimage',tread(2),'INT')
   if(tread(2)==1) dtsets(idtset)%npimage=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nppert',tread(3),'INT')
   if (tread(3)==1.and.optdriver==RUNL_RESPFN) dtsets(idtset)%nppert=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npkpt',tread(4),'INT')
   if(tread(4)==1) dtsets(idtset)%npkpt=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npspinor',tread(5),'INT')
   if(tread(5)==1) dtsets(idtset)%npspinor=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npfft',tread(6),'INT')
   if(tread(6)==1) dtsets(idtset)%npfft=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'npband',tread(7),'INT')
   if(tread(7)==1) dtsets(idtset)%npband=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bandpp',tread(8),'INT')
   if(tread(8)==1) dtsets(idtset)%bandpp=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_slk',tread(9),'INT')
   if(tread(9)==1) dtsets(idtset)%use_slk=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'np_slk',tread(10),'INT')
   if(tread(10)==1) dtsets(idtset)%np_slk=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'slk_rankpp',tread(12),'INT')
   if(tread(12)==1) dtsets(idtset)%slk_rankpp=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pw_unbal_thresh',tread0,'DPR')
   if(tread0==1) dtsets(idtset)%pw_unbal_thresh=dprarr(1)
   mpi_enregs(idtset)%pw_unbal_thresh=dtsets(idtset)%pw_unbal_thresh

   call intagm(dprarr,intarr,jdtset,marr,5,string(1:lenstr),'gpu_devices',tread0,'INT')
   if(tread0==1) dtsets(idtset)%gpu_devices(1:5)=intarr(1:5)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gpu_linalg_limit',tread(11),'INT')
   if(tread(11)==1) dtsets(idtset)%gpu_linalg_limit=intarr(1)


   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nphf',tread0,'INT')
   if(tread0==1) dtsets(idtset)%nphf=intarr(1)

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'autoparal',tread0,'INT')
   if(tread0==1) dtsets(idtset)%autoparal=intarr(1)

   wfoptalg_read=.false.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wfoptalg',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%wfoptalg=intarr(1)
     wfoptalg_read=.true.
   else
     if (dtsets(idtset)%usepaw==0) dtsets(idtset)%wfoptalg=0
     if (dtsets(idtset)%usepaw/=0) dtsets(idtset)%wfoptalg=10
     if ((optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS).and.dtsets(idtset)%paral_kgb/=0) dtsets(idtset)%wfoptalg=114
   end if

   ! Dump the list of irreducible perturbations and exit.
   if (dtsets(idtset)%paral_rf==-1.and.optdriver/=RUNL_NONLINEAR) then
     call get_npert_rbz(dtsets(idtset),nband_rbz,nkpt_rbz,npert)
     ABI_DEALLOCATE(nband_rbz)
     ABI_DEALLOCATE(nkpt_rbz)
     iexit = iexit + 1
   end if

!  From total number of procs, compute all possible distributions
!  Ignore exit flag if GW/EPH calculations because autoparal section is performed in screening/sigma/bethe_salpeter/eph
   call finddistrproc(dtsets,filnam,idtset,iexit,mband_upper,mpi_enregs(idtset),ndtset_alloc,tread)
   if (any(optdriver == [RUNL_SCREENING, RUNL_SIGMA, RUNL_BSE, RUNL_EPH, RUNL_NONLINEAR])) iexit = 0

   if ((optdriver/=RUNL_GSTATE.and.optdriver/=RUNL_GWLS).and. &
&   (dtsets(idtset)%npkpt/=1   .or.dtsets(idtset)%npband/=1.or.dtsets(idtset)%npfft/=1.or. &
&   dtsets(idtset)%npspinor/=1.or.dtsets(idtset)%bandpp/=1)) then
!&   .or.(dtsets(idtset)%iscf<0)) then
     dtsets(idtset)%npkpt=1 ; dtsets(idtset)%npspinor=1 ; dtsets(idtset)%npfft=1
     dtsets(idtset)%npband=1; dtsets(idtset)%bandpp=1  ; dtsets(idtset)%nphf=1
     dtsets(idtset)%paral_kgb=0
     MSG_COMMENT('For non ground state calculation, set bandpp, npfft, npband, npspinor npkpt and nphf to 1')
   end if

!  Take into account a possible change of paral_kgb (change of thwe default algorithm)
   if (.not.wfoptalg_read) then
     if (dtsets(idtset)%usepaw==0) dtsets(idtset)%wfoptalg=0
     if (dtsets(idtset)%usepaw/=0) dtsets(idtset)%wfoptalg=10
     if ((optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS).and.dtsets(idtset)%paral_kgb/=0) dtsets(idtset)%wfoptalg=114
     if (mod(dtsets(idtset)%wfoptalg,10)==4) then
       do iikpt=1,dtsets(idtset)%nkpt
         if (any(abs(dtsets(idtset)%kpt(:,iikpt))>tol8)) dtsets(idtset)%istwfk(iikpt)=1
       end do
     end if
   end if

   dtsets(idtset)%densfor_pred=2
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'densfor_pred',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%densfor_pred=intarr(1)
   else
     if (dtsets(idtset)%paral_kgb==1) dtsets(idtset)%densfor_pred=6
   end if
   if((dtsets(idtset)%iscf==5.or.dtsets(idtset)%iscf==6) &
&   .and. dtsets(idtset)%ionmov==4 .and. dtsets(idtset)%densfor_pred/=3 )then
     dtsets(idtset)%densfor_pred=3
     write(message, '(a,a,a)' )&
&     'When ionmov==4 and iscf==5 or 6, densfor_pred must be 3.',ch10,&
&     'Set densfor_pred to 3.'
     MSG_COMMENT(message)
   end if

#ifdef HAVE_LOTF
!  LOTF need densfor_pred=2
   if(dtsets(idtset)%ionmov==23) dtsets(idtset)%densfor_pred=2
#endif

   if (usepaw==0) then
     dtsets(idtset)%ortalg=2
   else
     dtsets(idtset)%ortalg=-2
   end if
   ortalg_read=.false.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ortalg',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%ortalg=intarr(1)
     ortalg_read=.true.
   else if (dtsets(idtset)%wfoptalg>=10 .and. dtsets(idtset)%ortalg>0) then
     dtsets(idtset)%ortalg=-dtsets(idtset)%ortalg
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'iomode',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%iomode=intarr(1)
   else
     if ((xmpi_mpiio==1).and.(dtsets(idtset)%paral_kgb==1)) dtsets(idtset)%iomode=IO_MODE_MPI
#ifdef HAVE_NETCDF_DEFAULT
     dtsets(idtset)%iomode=IO_MODE_ETSF
#endif
   end if

   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'pawmixdg',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%pawmixdg=intarr(1)
   else if (dtsets(idtset)%npfft>1.and.usepaw==1) then
     dtsets(idtset)%pawmixdg=1
   end if

   call initmpi_img(dtsets(idtset),mpi_enregs(idtset),-1)
   nproc=mpi_enregs(idtset)%nproc_cell

!  Cycle if the processor is not used
   if (mpi_enregs(idtset)%me<0.or.iexit>0) then
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     cycle
   end if

   response=0
   if (dtsets(idtset)%rfddk/=0 .or. dtsets(idtset)%rf2_dkdk/=0 .or. dtsets(idtset)%rf2_dkde/=0 .or. &
&   dtsets(idtset)%rfelfd/=0 .or. dtsets(idtset)%rfphon/=0 .or. dtsets(idtset)%rfstrs/=0 .or. &
&   dtsets(idtset)%rfuser/=0 .or. dtsets(idtset)%rfmagn/=0) response=1

!  If no MPI, set all npxxx variables to 1
   if (nproc==1) then
     dtsets(idtset)%npkpt    = 1 ; dtsets(idtset)%npband   = 1
     dtsets(idtset)%npfft    = 1 ; dtsets(idtset)%npspinor = 1
     dtsets(idtset)%nphf     = 1
   end if

!    --IF CUDA AND RECURSION:ONLY BAND PARALLELISATION
   if(dtsets(idtset)%tfkinfunc==2 .and. nproc/=1)then
     dtsets(idtset)%npband = dtsets(idtset)%npband*dtsets(idtset)%npkpt*dtsets(idtset)%npspinor*dtsets(idtset)%npfft
     dtsets(idtset)%npkpt = 1
     dtsets(idtset)%npfft = 1
     dtsets(idtset)%npspinor = 1
     write(message, '(5a,i6,a)' )&
&     'If HAVE_GPU_CUDA and recursion are used ',ch10,&
&     'only the band parallelisation is active, we set:',ch10,&
&     'npfft= 1, npkpt= 1, npband=',dtsets(idtset)%npband,' .'
     MSG_WARNING(message)
   end if

   if (dtsets(idtset)%npspinor>=2.and.dtsets(idtset)%nspinor==1) then
     dtsets(idtset)%npspinor=1
     dtsets(idtset)%npfft=2*dtsets(idtset)%npfft
     write(message,'(3a)')&
&     'npspinor is bigger than nspinor !',ch10,&
&     'We set npspinor to 1 ; we set npfft to 2*npfft'
     MSG_WARNING(message)
   end if

!  Some checks on parallelization data
   if(dtsets(idtset)%paral_kgb < 0 ) then
     cycle
   else if(dtsets(idtset)%paral_kgb/=0.and.(dtsets(idtset)%bandpp/=1.or.dtsets(idtset)%npband/=1.or.&
&     dtsets(idtset)%npfft/=1.or.dtsets(idtset)%npkpt/=1.or.dtsets(idtset)%npspinor/=1))then
     if(dtsets(idtset)%npkpt*dtsets(idtset)%npfft*dtsets(idtset)%npband*dtsets(idtset)%npspinor > nproc )then
       write(message,'(7a)')&
&       'The product of npkpt, npfft, npband and npspinor is bigger than the number of processors.',ch10,&
&       'The user-defined values of npkpt, npfft, npband or npspinor will be modified,',ch10,&
&       'in order to bring this product below nproc .',ch10,&
&       'At present, only a very simple algorithm is used ...'
       MSG_WARNING(message)

       if(dtsets(idtset)%npkpt*dtsets(idtset)%npband*dtsets(idtset)%npspinor <= nproc) then
         dtsets(idtset)%npfft=1
         MSG_WARNING('Set npfft to 1')
       else if(dtsets(idtset)%npkpt*dtsets(idtset)%npspinor <= nproc)then
         dtsets(idtset)%npfft=1
         dtsets(idtset)%npband=1
         MSG_WARNING('Set npfft and npband to 1')
       else if(dtsets(idtset)%npkpt <= nproc)then
         dtsets(idtset)%npfft=1
         dtsets(idtset)%npband=1
         dtsets(idtset)%npspinor=1
         MSG_WARNING('Set npfft ,npband and npspinor to 1')
       else
         dtsets(idtset)%npfft=1
         dtsets(idtset)%npband=1
         dtsets(idtset)%npkpt=1
         dtsets(idtset)%npspinor=1
         MSG_WARNING('Set npfft, npband, nspinor and npkpt to 1')
       end if
     else if(dtsets(idtset)%npkpt*dtsets(idtset)%npfft*dtsets(idtset)%npband*dtsets(idtset)%npspinor < nproc)then
       write(message,'(2a)')&
&       'The number of processor must not be greater than npfft*npband*npkpt*npsinor ',&
&       'when npfft or npkpt or npband or nspinor are chosen manually in the input file.'
       MSG_ERROR(message)
     end if
   end if

!  LOBPCG and ChebFi need paral_kgb=1 in parallel
   if ((dtsets(idtset)%npband*dtsets(idtset)%npfft>1).and. &
&   (mod(dtsets(idtset)%wfoptalg,10)==1.or.mod(dtsets(idtset)%wfoptalg,10)==4)) then
     dtsets(idtset)%paral_kgb=1
   end if

!  Check size of Scalapack communicator
#ifdef HAVE_LINALG_ELPA
   if(dtsets(idtset)%paral_kgb>0.and.dtsets(idtset)%np_slk>0) then
     np=min(dtsets(idtset)%np_slk,dtsets(idtset)%npband*dtsets(idtset)%npfft*dtsets(idtset)%npspinor)
     irow=int(sqrt(float(np)))
     do while(mod(np,irow)/=0)
       irow=irow-1
     end do
     icol=nproc/irow
     if (icol>mband_lower) then
       do while(icol>mband_lower)
         icol=icol-1
         do while(mod(np,icol)/=0)
           icol=icol-1
         end do
       end do
       dtsets(idtset)%np_slk=icol
       write(message,'(5a,i6,a)')&
&       'The number of band*fft*spinor processors was not consistent with',ch10,&
&       'the size of communicator used for ELPA library (np_slk).',ch10,&
&       'np_slk value has been adjusted to ',dtsets(idtset)%np_slk,'.'
       MSG_COMMENT(message)
     end if
   end if
#endif

   !Additional check in case of a parallelized Hartree-Fock calculation
   !   %usefock == option to perform Fock exchange calculation
   !   %nphf   == number of processors for Fock exchange calculation
   if ((dtsets(idtset)%usefock==1).and.(dtsets(idtset)%nphf/=1)) then

     if ((dtsets(idtset)%nphf<0).or.(dtsets(idtset)%nphf==0)) then
       MSG_ERROR('The value of variable nphf should be a non negative integer.')
     end if
     if (dtsets(idtset)%paral_kgb/=0) then
       message = 'Option paral_kgb should be turned off (value 0) for a parallelized Hartree-Fock calculation.'
       MSG_ERROR(message)
     end if
     if (response/=0) then
       message = 'A response function calculation is not yet possible with a parallelized Hartree-Fock calculation.'
       MSG_ERROR(message)
     end if
     if (dtsets(idtset)%npspinor>1) then
       message = 'The parallelism on spinors is not supported by a parallelized Hartree-Fock calculation.'
       MSG_ERROR(message)
     end if
     if (dtsets(idtset)%npkpt*dtsets(idtset)%nphf > nproc )then
       write(message,'(a,3(a,i0))') ch10,&
&       'The product of variables npkpt and nphf is bigger than the number of processors: nkpt= ',&
&       dtsets(idtset)%npkpt,' nphf= ',dtsets(idtset)%nphf  ,' and nproc= ', nproc
       MSG_ERROR(message)
     end if
   end if ! Fock

   !When using chebfi, the number of blocks is equal to the number of processors
   if((dtsets(idtset)%wfoptalg == 1)) then
     !Nband might have different values for different kpoint, but not bandpp.
     !In this case, we just use the largest nband, andthe input will probably fail
     !at the bandpp check later on
     dtsets(idtset)%bandpp = mband_upper / dtsets(idtset)%npband
   end if

!  Set mpi_enreg
   mpi_enregs(idtset)%paral_kgb=dtsets(idtset)%paral_kgb
   if(dtsets(idtset)%paral_kgb/=0)then
     mpi_enregs(idtset)%nproc_kpt=dtsets(idtset)%npkpt
     mpi_enregs(idtset)%nproc_fft=dtsets(idtset)%npfft
     mpi_enregs(idtset)%nproc_band=dtsets(idtset)%npband
     mpi_enregs(idtset)%nproc_spinor=min(dtsets(idtset)%npspinor,dtsets(idtset)%nspinor)
     mpi_enregs(idtset)%bandpp=dtsets(idtset)%bandpp
!    Additional setting in case of hybrid functional calculation => not yet tested (CMartins)
!     if (dtsets(idtset)%usefock==1) then
!       mpi_enregs(idtset)%nproc_hf = dtsets(idtset)%nphf
!       if (dtsets(idtset)%nphf>1) mpi_enregs(idtset)%paral_hf=1
!     end if
   else
     mpi_enregs(idtset)%bandpp = dtsets(idtset)%bandpp
!    Additional setting in case of a Fock exchange of PBE0 calculation
     if (dtsets(idtset)%usefock==1) then
       if (dtsets(idtset)%nphf>1) mpi_enregs(idtset)%paral_hf=1
       mpi_enregs(idtset)%nproc_hf = dtsets(idtset)%nphf
       if (dtsets(idtset)%npkpt/=1) then
         mpi_enregs(idtset)%nproc_kpt = dtsets(idtset)%npkpt
       else
         mpi_enregs(idtset)%nproc_kpt = mpi_enregs(idtset)%nproc_cell/mpi_enregs(idtset)%nproc_hf
       end if
     else
       mpi_enregs(idtset)%nproc_kpt = mpi_enregs(idtset)%nproc_cell
     end if
   end if

   if(dtsets(idtset)%paral_kgb>=0) then

!    Compute processor distribution over perturbations
     mpi_enregs(idtset)%paral_pert=dtsets(idtset)%paral_rf
     if (mpi_enregs(idtset)%paral_pert==1) then
       dtsets(idtset)%nppert=max(1,dtsets(idtset)%nppert)
       if(dtsets(idtset)%nppert>mpi_enregs(idtset)%nproc) then
         message=' The number of processors must not be smaller than nppert !'
         MSG_ERROR(message)
       end if
       call initmpi_pert(dtsets(idtset),mpi_enregs(idtset))
       mpi_enregs(idtset)%nproc_kpt = mpi_enregs(idtset)%nproc_cell
       nproc=mpi_enregs(idtset)%nproc_cell
     end if
!    Cycle if the processor is not used
     if (mpi_enregs(idtset)%me<0) then
       ABI_DEALLOCATE(intarr)
       ABI_DEALLOCATE(dprarr)
       cycle
     end if

!    Compute processor distribution over kpt (and eventually band-fft)
     call initmpi_grid(mpi_enregs(idtset))
     if(dtsets(idtset)%usewvl==1) mpi_enregs(idtset)%comm_fft=mpi_enregs(idtset)%comm_cell

!    Initialize tabs used for k/spin parallelism (with sequential-type values)
     ABI_ALLOCATE(mpi_enregs(idtset)%proc_distrb,(nkpt,mband_upper,nsppol))
     ABI_ALLOCATE(mpi_enregs(idtset)%my_kpttab,(nkpt))
     mpi_enregs(idtset)%proc_distrb(:,:,:)=0
     mpi_enregs(idtset)%my_kpttab(:)=(/(ii,ii=1,nkpt)/)
     mpi_enregs(idtset)%my_isppoltab(:)=1;if (dtsets(idtset)%nsppol==1) mpi_enregs(idtset)%my_isppoltab(2)=0

!    HF or hybrid calculation : initialization of the array distrb_hf
     if (dtsets(idtset)%usefock==1) then
       ABI_ALLOCATE(mpi_enregs(idtset)%distrb_hf,(dtsets(idtset)%nkpthf,dtsets(idtset)%nbandhf,1))
!      The dimension of distrb_hf are given by %nkpthf and %nbandhf.
!      We assume that there will be no dependence in spinpol for all the occupied states.
       mpi_enregs(idtset)%distrb_hf=0
     end if

!    Define k-points distribution (determine who I am)
!    Note that nkpt_me may differ from processor to processor
!    This fact will NOT be taken into account when
!    the memory needs will be evaluated in the subroutine memory.
!    Also, the reduction of k points due to symmetry in RF calculations
!    is NOT taken into account. This should be changed later ...
     nkpt_me=nkpt
     if(xmpi_paral==1 .and. dtsets(idtset)%usewvl == 0) then
       nkpt_me=0
       if(response==0 .or. (response==1 .and. dtsets(idtset)%efmas==1))then
         mpi_enregs(idtset)%paralbd=0
         call distrb2(mband_upper,dtsets(idtset)%nband,nkpt,nproc,nsppol,mpi_enregs(idtset))
         do iikpt=1,nkpt
           if(.not.(proc_distrb_cycle(mpi_enregs(idtset)%proc_distrb,iikpt,1,1,-1,mpi_enregs(idtset)%me_kpt)))&
&           nkpt_me=nkpt_me+1
         end do ! ikpt=1,nkpt
!        HF or hybrid calculation : define the occupied states distribution (in array distrb_hf)
         if (dtsets(idtset)%usefock==1) then
           call distrb2_hf(dtsets(idtset)%nbandhf,dtsets(idtset)%nkpthf,nproc,nsppol,mpi_enregs(idtset))
         end if
       else ! response==1
!        Wrongly assumes that the number of elements of the
!        k-point sets of the two spin polarizations is the maximal
!        value of one of these k-point sets ...
!        This is to be corrected when RF is implemented
!        for spin-polarized case.
         mpi_enregs(idtset)%paralbd=1
!        nproc=mpi_enregs(idtset)%nproc_cell*mpi_enregs(idtset)%nproc_pert
         call distrb2(mband_upper,dtsets(idtset)%nband,nkpt,nproc,nsppol,mpi_enregs(idtset))
         do isppol=1,nsppol
           nspink=0
           do iikpt=1,nkpt
             do iband=1,dtsets(idtset)%nband(iikpt+(isppol-1)*nkpt)
               if(mpi_enregs(idtset)%proc_distrb(iikpt,iband,isppol)==mpi_enregs(idtset)%me_cell)then
                 nspink=nspink+1
                 exit
               end if
             end do ! iband
           end do ! iikpt
           if(nspink>nkpt_me)nkpt_me=nspink
         end do ! isppol
!        Is nband present in input file or automatically estimated ?
         tnband=0
         call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nband',tnband,'INT')
!        If the number of bands was estimated, there might be a side effect
!        when the definitive number of bands is known. k points
!        might be attributed to different processors than the present
!        proc_distrb describes. At most, the number of k points could increase by 1 ...
         if(tnband==0)nkpt_me=nkpt_me+1
!        In any case, the maximal number of k points is nkpt
         if(nkpt_me>nkpt)nkpt_me=nkpt
       end if
     end if
   end if

!  Take care of mkmems. Use the generic name -mkmem- for mkmem as well as mkqmem
!  and mk1mem.
   nm_mkmem(1)='mkmem '
   nm_mkmem(2)='mkqmem'
   nm_mkmem(3)='mk1mem'

   do ii=1,3

!    Read in mkmem here if it is in the input file
!    TODO: mkmem is not supported any longer. These variables can be removed.
     if(ii==1)then
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mkmem',tread0,'INT')
     else if(ii==2)then
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mkqmem',tread0,'INT')
     else if(ii==3)then
       call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mk1mem',tread0,'INT')
     end if

!    Note that mkmem is used as a dummy variable, representing mkmem as well
!    as mkqmem, and mk1mem.
     if(tread0==1) then
       mkmem=intarr(1)
       if (mkmem<0) then
!        mkmem is unreasonable; must be zero or positive
         write(message, '(4a,i0,4a)')&
&         nm_mkmem(ii),' must be positive or null but ',nm_mkmem(ii),' =',mkmem,ch10,&
&         'Use default ',nm_mkmem(ii),' = nkpt .'
         MSG_WARNING(message)
         mkmem=nkpt
       end if

     else

       !  mkmem was not set in the input file so default to incore solution
       !write(message,'(6a)') &
       !'mpi_setup: ',nm_mkmem(ii),' undefined in the input file.','Use default ',nm_mkmem(ii),' = nkpt'
       !call wrtout(std_out,message,'COLL')
       mkmem=nkpt
     end if

!    Check whether nkpt distributed on the processors <= mkmem;
!    if so then may run entirely in core,
!    avoiding i/o to disk for wavefunctions and kg data.
!    mkmem/=0 to avoid i/o; mkmem==0 to use disk i/o for nkpt>=1.
     if (nkpt_me<=mkmem .and. mkmem/=0 ) then
       write(message, '(a,i0,a,a,a,i0,a)' ) &
        ' mpi_setup: With nkpt_me=',nkpt_me,' and ',nm_mkmem(ii),' = ',mkmem,', ground state wf handled in core.'
       if (prtvol > 0) call wrtout(std_out,message)
       if(nkpt_me<mkmem .and. nkpt_me/=0)then
         write(message,'(3a)')' Resetting ',nm_mkmem(ii),' to nkpt_me to save memory space.'
         mkmem=nkpt_me
         if (prtvol > 0) call wrtout(std_out,message)
       end if
     else if(mkmem/=0)then
       write(message, '(a,i0,3a,i0,5a)' ) &
       ' mpi_setup: With nkpt_me=',nkpt_me,'and ',nm_mkmem(ii),' = ',mkmem,&
       ' ground state wf require disk i/o.',ch10,&
       ' Resetting ',nm_mkmem(ii),' to zero to save memory space.'
       mkmem=0
       if (prtvol > 0) call wrtout(std_out,message)
     end if
     if(dtsets(idtset)%usewvl == 0 .or. dtsets(idtset)%usepaw==1)then
       if(ii==1)dtsets(idtset)%mkmem=mkmem
     end if
     if(ii==2)dtsets(idtset)%mkqmem=mkmem
     if(ii==3)dtsets(idtset)%mk1mem=mkmem

     if(dtsets(idtset)%usewvl == 1 .and. dtsets(idtset)%usepaw==1 )then
       if(dtsets(idtset)%mkmem .ne. dtsets(idtset)%nkpt) then
         MSG_ERROR("mkmem is not allowed for WVL+PAW")
       end if
     end if

   end do  ! End the loop on the three possiblities mkmem, mkqmem, mk1mem.

   if(dtsets(idtset)%paral_kgb==1) mpi_enregs(idtset)%paralbd=0

!  Check if some MPI processes are empty (MBPT code uses a complete different MPI algorithm)
   do_check = all(optdriver /= [RUNL_SCREENING, RUNL_SIGMA, RUNL_BSE, RUNL_EPH])
   if (dtsets(idtset)%usewvl == 0 .and. do_check) then
     if (.not.mpi_distrib_is_ok(mpi_enregs(idtset),mband_upper,&
&     dtsets(idtset)%nkpt,dtsets(idtset)%mkmem,nsppol,msg=message)) then
       write(message,'(5a)') trim(message),ch10,&
&       'YOU ARE STRONGLY ADVICED TO ACTIVATE AUTOMATIC PARALLELIZATION!',ch10,&
&       'PUT "AUTOPARAL=1" IN THE INPUT FILE.'
       MSG_WARNING(message)
     end if
   end if

!  call mpi_setup1(dtsets(idtset),jdtset,lenstr,mband_upper,mpi_enregs(idtset),string)
!  Printing of processor distribution
!  MPIWF : here, set up the complete ngfft, containing the information
!  for the parallelisation of the FFT
   call abi_io_redirect(new_io_comm=mpi_enregs(idtset)%comm_world)
   call libpaw_write_comm_set(mpi_enregs(idtset)%comm_world)

!  Default values for sequential case
   paral_fft=0; nproc_fft=1; me_fft=0

   if(dtsets(idtset)%usewvl == 0)then
     if(optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS) then
       paral_fft=1           ! parallelisation over FFT
       if (mpi_enregs(idtset)%nproc_cell>0) then
         if(mpi_enregs(idtset)%paral_kgb == 1) then

           if((dtsets(idtset)%use_gpu_cuda==1).and.(mpi_enregs(idtset)%nproc_fft/=1))then
             write(message,'(3a,i0)') &
&             'When use_gpu_cuda is on, the number of FFT processors, npfft, must be 1',ch10,&
&             'However, npfft=',mpi_enregs(idtset)%nproc_fft
             MSG_ERROR(message)
           end if

           if(modulo(dtsets(idtset)%ngfft(2),mpi_enregs(idtset)%nproc_fft)/=0)then
             write(message,'(3a,i0,a,i0)') &
&             'The number of FFT processors, npfft, should be a multiple of ngfft(2).',ch10,&
&             'However, npfft=',mpi_enregs(idtset)%nproc_fft,' and ngfft(2)=',dtsets(idtset)%ngfft(2)
             MSG_BUG(message)
           end if

           do iikpt=1,nkpt*nsppol
             iikpt_modulo = modulo(iikpt,nkpt)+1
             if ((dtsets(idtset)%istwfk(iikpt_modulo)==2)) then !.and.(dtsets(idtset)%ngfft(7)==401)) then
               if ((mpi_enregs(idtset)%bandpp==0).or. &
               ((mpi_enregs(idtset)%bandpp/=1).and.(modulo(mpi_enregs(idtset)%bandpp,2)/=0))) then
                 write(message,'(3a,i0)') &
&                 'The number bandpp should be 1 or a multiple of 2',ch10,&
&                 'However, bandpp=',mpi_enregs(idtset)%bandpp
                 MSG_BUG(message)
               end if
               if(modulo(dtsets(idtset)%nband(iikpt),mpi_enregs(idtset)%nproc_band*mpi_enregs(idtset)%bandpp)/=0)then
                 write(message,'(3a,i0,a,i0)') &
&                 'The number of bands for the k-point, nband_k, should be a multiple of nproc_band*bandpp.',ch10,&
&                 'However, nband_k=',dtsets(idtset)%nband(iikpt),' and nproc_band*bandpp=', &
&                 mpi_enregs(idtset)%nproc_band* mpi_enregs(idtset)%bandpp
                 MSG_BUG(message)
               end if
             else if ((dtsets(idtset)%istwfk(iikpt_modulo)==2) .and. (dtsets(idtset)%ngfft(7)==400)) then
               MSG_BUG('The fftalg=400 with istwfk=2 is not valid')
             else
               if(modulo(dtsets(idtset)%nband(iikpt),mpi_enregs(idtset)%nproc_band*mpi_enregs(idtset)%bandpp)/=0)then
                 write(message,'(3a,i0,a,i0)') &
&                 'The number of band for the k-point, nband_k, should be a multiple of nproc_band*bandpp.',ch10,&
&                 'However, nband_k=',dtsets(idtset)%nband(iikpt),' and nproc_band*bandpp=', &
&                 mpi_enregs(idtset)%nproc_band* mpi_enregs(idtset)%bandpp
                 MSG_BUG(message)
               end if
               if ((mpi_enregs(idtset)%bandpp==0)) then
                 write(message,'(a,i0,2a,i0,2a,i0)')&
&                 'The number bandpp should not be 0 with fftalg=',dtsets(idtset)%ngfft(7),ch10,&
&                 'and istwfk=',dtsets(idtset)%istwfk(iikpt_modulo),ch10,&
&                 'However, bandpp=',mpi_enregs(idtset)%bandpp
                 MSG_BUG(message)
               end if
             end if
           end do

           if (xmpi_paral==1) then
             if(modulo(nkpt*nsppol,mpi_enregs(idtset)%nproc_kpt)/=0)then
               write(message,'(3a,i0,a,i0)') &
&               'The number of KPT processors, npkpt, should be a multiple of nkpt*nsppol.',ch10,&
&               'However, npkpt=',mpi_enregs(idtset)%nproc_kpt,' and nkpt*nsppol=',nkpt*nsppol
               MSG_WARNING(message)
             end if
           end if
         else
           do iikpt=1,nkpt*nsppol
             iikpt_modulo = modulo(iikpt,nkpt)+1
             if(modulo(dtsets(idtset)%nband(iikpt),mpi_enregs(idtset)%nproc_band*mpi_enregs(idtset)%bandpp)/=0)then
               write(message,'(3a,i0,a,i0)') &
&               'The number of band for the k-point, nband_k, should be a multiple of npband*bandpp.',ch10,&
&               'However, nband_k=',dtsets(idtset)%nband(iikpt),' and npband*bandpp=', &
&               mpi_enregs(idtset)%nproc_band* mpi_enregs(idtset)%bandpp
               MSG_BUG(message)
             end if
           end do
         end if
       end if
       nproc_fft=mpi_enregs(idtset)%nproc_fft
       me_fft=mpi_enregs(idtset)%me_fft
     end if
   end if

!  Compute mgfft,mpw,nfft for this data set (it is dependent of mpi_enreg)
   ABI_ALLOCATE(istwfk,(nkpt))
   ABI_ALLOCATE(kpt_with_shift,(3,nkpt))

   ! Set the default value of fftalg for given npfft but allow the user to override it.
   ! Warning: If you need to change npfft, **DO IT** before this point so that here we get the correct fftalg
   dtsets(idtset)%ngfft(7) = fftalg_for_npfft(dtsets(idtset)%npfft)
   dtsets(idtset)%ngfftdg(7) = fftalg_for_npfft(dtsets(idtset)%npfft)

   fftalg_read=.false.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fftalg',tread0,'INT')

   if (tread0==1) then
     dtsets(idtset)%ngfft(7)=intarr(1)
     if (usepaw==1) dtsets(idtset)%ngfftdg(7)=intarr(1)
     fftalg_read=.true.
   end if

   ecut     =dtsets(idtset)%ecut
   dilatmx  =dtsets(idtset)%dilatmx
   ngfft(:) =dtsets(idtset)%ngfft(:)
   istwfk(:)=dtsets(idtset)%istwfk(1:nkpt)
   nsym     =dtsets(idtset)%nsym

   nqpt=dtsets(idtset)%nqpt
   qphon(:)=zero;if(nqpt/=0) qphon(:)=dtsets(idtset)%qptn(:)

   ABI_ALLOCATE(symrel,(3,3,nsym))
   symrel(:,:,1:nsym)=dtsets(idtset)%symrel(:,:,1:nsym)
   ecut_eff=ecut*dilatmx**2

   if (usepaw==1) then
     call wrtout(std_out,'getng is called for the coarse grid:','COLL')
   end if
   kpt=k0; if (response==1.and.usepaw==1) kpt=qphon ! this is temporary

   call getng(dtsets(idtset)%boxcutmin,ecut_eff,gmet,kpt,me_fft,mgfft,nfft,&
&   ngfft,nproc_fft,nsym,paral_fft,symrel,&
&   use_gpu_cuda=dtsets(idtset)%use_gpu_cuda)

   dtsets(idtset)%ngfft(:)=ngfft(:)
   dtsets(idtset)%mgfft=mgfft
   dtsets(idtset)%nfft=nfft
   kpt_with_shift(:,:)=dtsets(idtset)%kpt(:,1:nkpt)/dtsets(idtset)%kptnrm

   exchn2n3d=dtsets(idtset)%exchn2n3d
   nproc_fft=ngfft(10) ; me_fft=ngfft(11)
   fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=mod(fftalg,10)

   ! Initialize tables for MPI-FFT.
   call init_distribfft(mpi_enregs(idtset)%distribfft,'c',mpi_enregs(idtset)%nproc_fft,ngfft(2),ngfft(3))

   if(response/=0)then
!    This value of mpw is used in the first part of respfn.f
     call getmpw(ecut_eff,exchn2n3d,gmet,istwfk,kpt_with_shift,mpi_enregs(idtset),mpw_k,nkpt)
   end if
   if(nqpt/=0)then
     kpt_with_shift(1,:)=kpt_with_shift(1,:)+qphon(1)
     kpt_with_shift(2,:)=kpt_with_shift(2,:)+qphon(2)
     kpt_with_shift(3,:)=kpt_with_shift(3,:)+qphon(3)
   end if
   if (dtsets(idtset)%usewvl == 0) then
     call getmpw(ecut_eff,exchn2n3d,gmet,istwfk,kpt_with_shift,mpi_enregs(idtset),mpw,nkpt)
     ! Allocate tables for parallel IO of the wavefunctions.
     if( xmpi_mpiio==1 .and. mpi_enregs(idtset)%paral_kgb == 1 .and. &
&     any(dtsets(idtset)%iomode == [IO_MODE_MPI, IO_MODE_ETSF])) then
       ABI_ALLOCATE(mpi_enregs(idtset)%my_kgtab,(mpw,dtsets(idtset)%mkmem))
     end if
   else
     mpw = 0
   end if

!  The dimensioning, in the RF case, should be done only with mpw,
!  but mpw is used in the first part of respfn.f, and should at least
!  be equal to mpw_k . The chosen way to code is not optimal, only convenient :
!  it leads to a small waste of memory.
   if(response/=0 .and. mpw_k>mpw)mpw=mpw_k
   dtsets(idtset)%ngfft(:)=ngfft(:)

!  Initialize ngfftc to the initial guess for the coarse mesh
   ngfftc(:) = 2

!  In case of PAW, compute fine FFT parameters
   if (usepaw==1) then
     ecutdg_eff=dtsets(idtset)%pawecutdg*dtsets(idtset)%dilatmx**2
     ngfftdg(:)=dtsets(idtset)%ngfftdg(:)
     call wrtout(std_out,'getng is called for the fine grid:','COLL')
!    Start with the coarse mesh as an initial guess for the fine mesh
!    This ensures that the fine mesh will not be any coarser than the coarse mesh in each dimension
     ngfftc(:) = ngfft(1:3)
     kpt=k0; if (response==1.and.usepaw==1) kpt=qphon  ! this is temporary

     call getng(dtsets(idtset)%bxctmindg,ecutdg_eff,gmet,kpt,me_fft,mgfftdg,&
&     nfftdg,ngfftdg,nproc_fft,nsym,paral_fft,symrel,ngfftc,&
&     use_gpu_cuda=dtsets(idtset)%use_gpu_cuda)

     dtsets(idtset)%ngfftdg(:)=ngfftdg(:)
     dtsets(idtset)%mgfftdg=mgfftdg
     dtsets(idtset)%nfftdg=nfftdg
!    Compute fft distribution for fine grid
     fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=mod(fftalg,10)
     call init_distribfft(mpi_enregs(idtset)%distribfft,'f', mpi_enregs(idtset)%nproc_fft,ngfftdg(2),ngfftdg(3))
   end if

   dtsets(idtset)%mpw=mpw
   ABI_DEALLOCATE(symrel)
   ABI_DEALLOCATE(istwfk)
   ABI_DEALLOCATE(kpt_with_shift)
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)

!  Initialize data for the parallelization over atomic sites (PAW)
   if (dtsets(idtset)%natom==1) dtsets(idtset)%paral_atom=0
   if (dtsets(idtset)%usepaw==0) dtsets(idtset)%paral_atom=0
   if (dtsets(idtset)%usewvl/=0) dtsets(idtset)%paral_atom=0
   if (dtsets(idtset)%usedmft==1) dtsets(idtset)%paral_atom=0
   if (optdriver/=RUNL_GSTATE.and.optdriver/=RUNL_RESPFN.and.optdriver/=RUNL_GWLS) dtsets(idtset)%paral_atom=0
   if (dtsets(idtset)%macro_uj/=0) dtsets(idtset)%paral_atom=0

   call initmpi_atom(dtsets(idtset),mpi_enregs(idtset))

!  In case of the use of a GPU (Cuda), some defaults can change
!  according to a threshold on matrix sizes
   if (dtsets(idtset)%use_gpu_cuda==1.or.dtsets(idtset)%use_gpu_cuda==-1) then
     if (optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS) then
       vectsize=dtsets(idtset)%mpw*dtsets(idtset)%nspinor/dtsets(idtset)%npspinor
       if (all(dtsets(idtset)%istwfk(:)==2)) vectsize=2*vectsize
       blocksize=dtsets(idtset)%npband*dtsets(idtset)%bandpp
       if (dtsets(idtset)%paral_kgb==0) blocksize=dtsets(idtset)%npfft
       if ((vectsize*blocksize**2)>=dtsets(idtset)%gpu_linalg_limit) then
         if (.not.wfoptalg_read) then
           dtsets(idtset)%wfoptalg=14
           if (.not.fftalg_read) then
             dtsets(idtset)%ngfft(7) = fftalg_for_npfft(dtsets(idtset)%npfft)
             if (usepaw==1) dtsets(idtset)%ngfftdg(7) = fftalg_for_npfft(dtsets(idtset)%npfft)
           end if
           if (.not.ortalg_read) dtsets(idtset)%ortalg=-abs(dtsets(idtset)%ortalg)
         end if
       end if
     end if
   end if

!  initialize data for the parallelization for WVL:
   if(dtsets(idtset)%usewvl==1) then
     mpi_enregs(idtset)%comm_wvl=mpi_enregs(idtset)%comm_cell
     mpi_enregs(idtset)%nproc_wvl=xmpi_comm_size(mpi_enregs(idtset)%comm_wvl)
     mpi_enregs(idtset)%me_wvl=xmpi_comm_rank(mpi_enregs(idtset)%comm_wvl)
   end if

 end do

!This is not a very clean exit in case of paral_kgb<0
 if (iexit/=0)then
   message="Stopping now!"
   MSG_STOP(message)
 end if

 DBG_EXIT("COLL")

end subroutine mpi_setup
!!***

!!****f* ABINIT/finddistrproc
!! NAME
!! finddistrproc
!!
!! FUNCTION
!!   Given a total number of processors, find a suitable distribution
!!   that fill all the different levels of parallelization
!!   (npimage, nppert, npkpt, npspinor, npband, npfft, bandpp)
!!   Also determine parameters of parallel Linear Algebra routines
!!   (use_slk, np_slk, gpu_linalg_limit)
!!
!! INPUTS
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   for all datasets; at this stage only datasets with index lower than
!!   idtset are already initalized
!!  filnam(5)=character strings giving file names
!!  idtset=number of the current dataset
!!  mpi_enreg=information about MPI parallelization
!!  mband=maximum number of bands.
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set
!!  tread(11)=flags indicating wether parallel input parameters were read from input file
!!            tread(1)  : paral_kgb      tread(6) : npfft
!!            tread(2)  : npimage        tread(7) : npband
!!            tread(3)  : nppert         tread(8) : bandpp
!!            tread(4)  : npkpt          tread(9) : use_slk
!!            tread(5)  : nspinor        tread(10): np_slk
!!            tread(11) : gpu_linalg_limit
!!
!! SIDE EFFECTS
!!  iexit= if incremented, an exit is required
!!  dtset%paral_kgb= flag for band-fft parallelism
!!  dtset%npimage  = number of processors for parallelisation over image
!!  dtset%nppert   = number of processors for parallelisation over perturbations
!!  dtset%npspinor = number of processors for parallelisation on k points
!!  dtset%npkpt    = number of processors for parallelisation on k points
!!  dtset%npfft    = number of processors for parallelisation on fft grid
!!  dtset%npband   = number of processors for parallelisation on bands
!!  dtset%nphf     = number of processors for parallelisation on occupied states for fock exchange
!!  dtset%bandpp   = internal parameter for lobpcg parallelisation algorithm
!!  dtset%use_slk  = flag for ScalaPAck use
!!  dtset%np_slk   = number of processors used in ScaLapack routines
!!  dtset%gpu_linalg_limit=threshold activating Linear Algebra on GPU
!!
!! PARENTS
!!      mpi_setup
!!
!! CHILDREN
!!      compute_kgb_indicator,get_npert_rbz,hdr_free,hdr_read_from_fname
!!      initmpi_world,kpgcount,metric,mkfilename,mkrdim,sort_dp,wrtout
!!      xmpi_bcast
!!
!! SOURCE

 subroutine finddistrproc(dtsets,filnam,idtset,iexit,mband,mpi_enreg,ndtset_alloc,tread)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset,mband,ndtset_alloc
 integer,intent(inout) :: iexit
 type(dataset_type),intent(inout),target :: dtsets(0:ndtset_alloc)
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: tread(11)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
!scalars
!128 should be a reasonable maximum for npfft (scaling is very poor for npfft>20)
 integer,parameter :: ALGO_NOT_SET=-1, ALGO_DEFAULT_PAR=2
 integer,parameter :: ALGO_CG=0, ALGO_LOBPCG_OLD=1, ALGO_LOBPCG_NEW=2, ALGO_CHEBFI=3
 integer,parameter :: NPFMAX=128,BLOCKSIZE_MAX=3000,MAXBAND_PRINT=10
 integer,parameter :: MAXCOUNT=250,MAXPRINT=10,MAXBENCH=25,MAXABIPY=5,NPF_CUTOFF=20
 real(dp),parameter :: relative_nband_range=0.025
 integer :: wf_algo,wf_algo_global,bpp,bpp_max,bpp_min,optdriver,autoparal,nblocks,blocksize
 integer :: npi_max,npi_min,npc,npc_max,npc_min
 integer :: npk,npk_max,npk_min,npp_max,npp_min
 integer :: nps,nps_max,nps_min,npf,npf_max,npf_min
 integer :: npb,npb_max,npb_min,max_ncpus,ount,paral_kgb
 integer :: work_size,nks_per_proc,tot_ncpus
 integer :: ib1,ib2,ibest,icount,ii,imin,jj,kk,mcount,mcount_eff,mpw
 integer :: n2,n3,ncell_eff,ncount,nimage_eff,nkpt_eff,npert_eff
 integer :: nproc,nproc1,nprocmin,np_slk,nthreads,use_linalg_gpu,omp_ncpus
 logical :: dtset_found,file_found,first_bpp,iam_master
 logical :: with_image,with_pert,with_kpt,with_spinor,with_fft,with_band,with_bandpp,with_thread
 real(dp):: acc_c,acc_k,acc_kgb,acc_kgb_0,acc_s,ecut_eff,eff,ucvol,weight0
 character(len=9) :: suffix
 character(len=20) :: strg
 character(len=500) :: msg,msgttl
 character(len=fnlen) :: filden
 type(hdr_type) :: hdr0
!arrays
 integer :: idum(1),idum3(3),ngmax(3),ngmin(3)
 integer,allocatable :: nband_best(:),isort(:),jdtset_(:)
 integer,allocatable :: my_algo(:),my_distp(:,:),nproc_best(:)
 integer,pointer :: nkpt_rbz(:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),allocatable :: weight(:)
 real(dp),pointer :: nband_rbz(:,:)
 type(dataset_type),pointer :: dtset

!******************************************************************

 DBG_ENTER("COLL")

!Select current dataset
 dtset => dtsets(idtset)

!Is automatic parallelization activated?
 autoparal = dtset%autoparal
 if (autoparal==0) return

!Is it available
 if ((dtset%usefock==1).AND.(dtset%nphf/=1)) then
   msg="autoparal>0 not available for Hartree-Fock or hybrid XC calculations!"
   MSG_ERROR(msg)
 end if
 if ((autoparal>1).and.dtset%wfoptalg/=4.and.dtset%wfoptalg/=14) then
   msg="autoparal>1 only available for the old LOBPCG algorithm (wfoptalg=4/14)!"
   MSG_ERROR(msg)
 end if

!Unit number used for outputting the autoparal sections
 ount = ab_out

!Handy local variables
 iam_master = (mpi_enreg%me==0)
 optdriver = dtset%optdriver
 max_ncpus = dtset%max_ncpus ; if (dtset%paral_kgb<0) max_ncpus=abs(dtset%paral_kgb)
 nthreads=xomp_get_max_threads()
 nproc=mpi_enreg%nproc
 if (max_ncpus>0) nproc = dtset%max_ncpus/nthreads
 if (xmpi_paral==0.and.max_ncpus<=0) nproc=1

 nprocmin=2
 if (xmpi_paral==1.and.max_ncpus<=0) nprocmin=max(2,nproc-100)
 if (max_ncpus>0.and.autoparal/=0) nprocmin=1

 wf_algo_global=ALGO_NOT_SET
 if (dtset%wfoptalg==0.and.tread(1)==1) wf_algo_global=ALGO_CG
 if (dtset%wfoptalg==4.or.dtset%wfoptalg==14) wf_algo_global=ALGO_LOBPCG_OLD
 if (dtset%wfoptalg==114) wf_algo_global=ALGO_LOBPCG_NEW
 if (dtset%wfoptalg==1) wf_algo_global=ALGO_CHEBFI

!Some peculiar cases (with direct exit)
 if (max_ncpus<=0) then
   if (nproc==1.and.max_ncpus<=0) then
     if (tread(1)==0.or.xmpi_paral==0) dtset%paral_kgb= 0
     if (tread(2)==0.or.xmpi_paral==0) dtset%npimage  = 1
     if (tread(3)==0.or.xmpi_paral==0) dtset%nppert   = 1
     if (tread(4)==0.or.xmpi_paral==0) dtset%npspinor = 1
     if (tread(5)==0.or.xmpi_paral==0) dtset%npkpt    = 1
     if (tread(6)==0.or.xmpi_paral==0) dtset%npfft    = 1
     if (tread(7)==0.or.xmpi_paral==0) dtset%npband   = 1
     if (tread(8)==0.or.xmpi_paral==0) dtset%bandpp   = 1
     if (tread(9)==0.or.xmpi_paral==0) dtset%use_slk  = 0
     if (tread(10)==0.or.xmpi_paral==0) dtset%np_slk  = 1000000
     return
   end if
   if ((dtset%optdriver/=RUNL_GSTATE.and.dtset%optdriver/=RUNL_RESPFN.and.dtset%optdriver/=RUNL_GWLS).or. &
&    (dtset%optdriver==RUNL_GSTATE.and.dtset%usewvl==1)) then
     dtset%paral_kgb= 0
     dtset%npimage  = max(1,dtset%npimage)
     dtset%nppert   = max(1,dtset%nppert)
     dtset%npspinor = max(1,dtset%npspinor)
     dtset%npkpt    = max(1,dtset%npkpt)
     dtset%npfft    = max(1,dtset%npfft)
     dtset%npband   = max(1,dtset%npband)
     dtset%bandpp   = max(1,dtset%bandpp)
     return
   end if
 end if

!Need the metric tensor
 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Determine some quantities related to plane waves
!  - Crude estimation of the number of PW
!  - Number of G vectors in each direction
 mpw=0;ngmin=0;ngmax=0
 if (optdriver==RUNL_GSTATE) then
   ecut_eff = dtset%ecut*dtset%dilatmx**2
   mpw = nint(ucvol*((two*ecut_eff)**1.5_dp)/(six*pi**2)) ! Crude estimation
   if (all(dtset%istwfk(1:dtset%nkpt)>1)) mpw=mpw/2+1
   call kpgcount(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,dtset%kpt,ngmax,ngmin,dtset%nkpt)
   write(msg,'(a,i8)') ' getmpw sequential formula gave: ',mpw
   call wrtout(std_out,msg,'COLL')
 end if

!Parallelization over images
 npi_min=1;npi_max=1;nimage_eff=1
 if (optdriver==RUNL_GSTATE) then
   nimage_eff=dtset%ndynimage
   if (dtset%ntimimage<=1) nimage_eff=dtset%nimage
   npi_min=max(1,dtset%npimage)
   npi_max=min(nproc,nimage_eff)
   if (tread(2)==1) npi_max=dtset%npimage
 end if

!Parallelization over k-points and spin components (GS)
 npk_min=1;npk_max=1;nkpt_eff=0
 if (optdriver==RUNL_GSTATE) then
   nkpt_eff=dtset%nkpt*dtset%nsppol
   npk_min=max(1,dtset%npkpt)
   npk_max=min(nproc,nkpt_eff)
   if (tread(4)==1) npk_max=dtset%npkpt
 end if

!Parallelization over perturbations, k-points and spin components (DFPT)
 npp_min=1;npp_max=1;npert_eff=1
 if (optdriver==RUNL_RESPFN) then
   if (dtset%paral_rf==1) then
     call get_npert_rbz(dtset,nband_rbz,nkpt_rbz,npert_eff)
     do jj=1,npert_eff
       ii=dtset%nsppol*nkpt_rbz(jj)*maxval(nband_rbz(:,jj))
       nkpt_eff=max(nkpt_eff,ii)
     end do
     npp_min=max(1,dtset%nppert)
     npp_max=min(nproc,npert_eff)
     if (tread(3)==1) then
       npp_max=dtset%nppert
       if (npp_max>npert_eff) then
         npp_min=npert_eff;npp_max=npert_eff
         msg='nppert is bigger than npert; we set nppert=npert'
         MSG_WARNING(msg)
       end if
     end if
     npk_min=1
     npk_max=min(nproc,nkpt_eff)
     ABI_DEALLOCATE(nkpt_rbz)
     ABI_DEALLOCATE(nband_rbz)
   else
     nkpt_eff=nproc
     npk_min=nproc-5
     npk_max=nproc
   end if
 end if

!Parallelization over spinorial components
 nps_min=1;nps_max=1
 if (optdriver==RUNL_GSTATE) then
   nps_min=max(1,dtset%npspinor)
   nps_max=min(nproc,dtset%nspinor)
   if (tread(5)==1) nps_max=dtset%npspinor
 end if

!KGB Parallelization

 npf_min=1;npf_max=1
 npb_min=1;npb_max=1
 bpp_min=1;bpp_max=1
 n2=0;n3=0
 if (dtset%optdriver==RUNL_GSTATE) then

!  >> FFT level
   npf_min=max(1,dtset%npfft)
   npf_min=min(npf_min,ngmin(2))
   npf_max=min(nproc,NPFMAX)
   if (tread(6)==1) then
     npf_max=dtset%npfft
     if (npf_max>ngmin(2)) then
       write(msg,'(3a)') &
&       "Value of npfft given in input file is too high for the FFT grid!",ch10,&
&       "Action: decrease npfft or increase FFT grid (ecut, ngfft, ...)."
       MSG_ERROR(msg)
     end if
   end if
   npf_max=min(npf_max,ngmin(2))
   !Deactivate MPI FFT parallelism for GPU
   if (dtset%use_gpu_cuda==1) then
     npf_min=1;npf_max=1
   end if
   !Deactivate MPI FFT parallelism for GPU
   if (tread(1)==1.and.dtset%paral_kgb==0) then
     npf_min=1;npf_max=1
   end if
   !Deactivate MPI FFT parallelism for multi-threaded LOBPCG / CHEBFI
   if ((wf_algo_global==ALGO_LOBPCG_NEW.or.wf_algo_global==ALGO_CHEBFI).and.nthreads>1) then
     npf_min=1;npf_max=1
   end if

!  Number of FFT procs has to be a multiple of FFT grid sizes
!  In case of a restart from a density file, it has to be
!  compatible with the FFT grid used for the density
   n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
   if (n2==0.and.n3==0.and.(dtset%getden/=0.or.dtset%irdden/=0.or.dtset%iscf<0)) then
     dtset_found=.false.;file_found=.false.
     !1-Try to find ngfft from previous dataset
     if (dtset%getden/=0) then
       do ii=1,ndtset_alloc
         jj=dtset%getden;if (jj<0) jj=dtset%jdtset+jj
         if (dtsets(ii)%jdtset==jj) then
           dtset_found=.true.
           n2=dtsets(ii)%ngfftdg(2);n3=dtsets(ii)%ngfftdg(3)
         end if
       end do
     end if
     !2-If not found, try to extract ngfft from density file
     if (.not.dtset_found) then
       !Retrieve file name
       suffix='_DEN';if (dtset%nimage>1) suffix='_IMG1_DEN'
       ABI_ALLOCATE(jdtset_,(0:ndtset_alloc))
       jdtset_=0;if(ndtset_alloc/=0) jdtset_(0:ndtset_alloc)=dtsets(0:ndtset_alloc)%jdtset
       call mkfilename(filnam,filden,dtset%getden,idtset,dtset%irdden,jdtset_,ndtset_alloc,suffix,'den',ii)
       ABI_DEALLOCATE(jdtset_)
       !Retrieve ngfft from file header
       idum3=0
       if (mpi_enreg%me==0) then
         inquire(file=trim(filden),exist=file_found)
         if (file_found) then
           call hdr_read_from_fname(hdr0,filden,ii,xmpi_comm_self)
           idum3(1:2)=hdr0%ngfft(2:3);if (file_found) idum3(3)=1
           call hdr_free(hdr0)
           MSG_WARNING("Cannot find filden"//filden)
         end if
       end if
       call xmpi_bcast(idum3,0,mpi_enreg%comm_world,ii)
       n2=idum3(1);n3=idum3(2);file_found=(idum3(3)/=0)
     end if
   end if

!  >> Band level
   npb_min=max(1,dtset%npband)
   npb_max=min(nproc,mband)
   if (tread(7)==1) npb_max=dtset%npband
   if (tread(1)==1.and.dtset%paral_kgb==0) then
     npb_min=1;npb_max=1
   end if

!  >> banddp level
   bpp_min=max(1,dtset%bandpp)
   bpp_max=mband
   if (wf_algo_global==ALGO_LOBPCG_OLD) bpp_max=max(4,nint(mband/10.)) ! reasonnable bandpp max
   if (tread(8)==1) bpp_max=dtset%bandpp
   if (wf_algo_global==ALGO_CHEBFI) bpp_min=1 ! bandpp not used with ChebFi
   if (wf_algo_global==ALGO_CHEBFI) bpp_max=1

 end if ! RUNL_GSTATE

!Disable KGB parallelisation in some cases:
!  - no GS
!  - paral_kgb=0 present in input file
!  - nstep=0
!  - Hartree-Fock or hybrid calculation (for now on)
 if ( (optdriver/=RUNL_GSTATE).or.(dtset%paral_kgb==0.and.tread(1)==1).or. &
&     (dtset%nstep==0).or.(dtset%usefock==1)) then
   nps_min=1; nps_max=1
   npf_min=1; npf_max=1
   npb_min=1; npb_max=1
   bpp_min=1; bpp_max=1
 end if

!Which levels of parallelism do we have?
 with_image =(npi_min/=1.or.npi_max/=1)
 with_pert  =(npp_min/=1.or.npp_max/=1)
 with_kpt   =(npk_min/=1.or.npk_max/=1)
 with_spinor=(nps_min/=1.or.nps_max/=1)
 with_fft   =(npf_min/=1.or.npf_max/=1)
 with_band  =(npb_min/=1.or.npb_max/=1)
 with_bandpp=(bpp_min/=1.or.bpp_max/=1)
 with_thread=(nthreads>1)

!Allocate lists
 ABI_ALLOCATE(my_distp,(10,MAXCOUNT))
 ABI_ALLOCATE(weight,(MAXCOUNT))
 ABI_ALLOCATE(my_algo,(MAXCOUNT))
 my_distp(1:7,:)=0;weight(:)=zero
 my_distp(8,:)=dtset%use_slk
 my_distp(9,:)=dtset%np_slk
 my_distp(10,:)=dtset%gpu_linalg_limit
 my_algo(:)=wf_algo_global
 icount=0;imin=1

!Cells= images or perturbations
 npc_min=1;npc_max=1;ncell_eff=1
 if (optdriver==RUNL_GSTATE) then
   ncell_eff=nimage_eff;npc_min=npi_min;npc_max=npi_max
 end if
 if (optdriver==RUNL_RESPFN) then
   ncell_eff=npert_eff;npc_min=npp_min;npc_max=npp_max
 end if

!Loop over all possibilities
!Computation of weight~"estimated acceleration"
!================================================================

!Cells= images or perturbations
 npc_min=1;npc_max=1;ncell_eff=1
 if (optdriver==RUNL_GSTATE) then
   ncell_eff=nimage_eff;npc_min=npi_min;npc_max=npi_max
 end if
 if (optdriver==RUNL_RESPFN) then
   ncell_eff=npert_eff;npc_min=npp_min;npc_max=npp_max
 end if

!>>>>> CELLS
 do npc=npc_min,npc_max
   acc_c=one;if (npc>1) acc_c=0.99_dp*speedup_fdp(ncell_eff,npc)

!  >>>>> K-POINTS
   do npk=npk_min,npk_max
!    -> for DFPT runs, impose that nsppol divide npk
     if (optdriver==RUNL_RESPFN.and.modulo(npk,dtset%nsppol)>0.and.npk>1) cycle
     acc_k=one;if (npk>1) acc_k=0.96_dp*speedup_fdp(nkpt_eff,npk)

!    >>>>> SPINORS
     do nps=nps_min,nps_max
       acc_s=one;if (nps>1) acc_s=0.85_dp*speedup_fdp(dtset%nspinor,nps)

!      >>>>> FFT
       do npf=npf_min,npf_max
!        -> npf should divide ngfft if set (if unset, ngfft=0 so the modulo test is ok)
         if((modulo(n2,npf)>0).or.(modulo(n3,npf)>0)) cycle
!        -> npf should be only divisible by 2, 3 or 5
         ii=npf
         do while (modulo(ii,2)==0)
           ii=ii/2
         end do
         do while (modulo(ii,3)==0)
           ii=ii/3
         end do
         do while (modulo(ii,5)==0)
           ii=ii/5
         end do
         if(ii/=1) cycle

!        Change algo if npfft>1
         wf_algo=wf_algo_global
         if (optdriver==RUNL_GSTATE.and.npf>1.and. &
&            wf_algo_global==ALGO_NOT_SET) wf_algo=ALGO_DEFAULT_PAR

!        FFT parallelism not compatible with multithreading
         if (wf_algo==ALGO_LOBPCG_NEW.or.wf_algo==ALGO_CHEBFI) then
           if (nthreads>1.and.npf>1) cycle
         end if

!        >>>>> BANDS
         do npb=npb_min,npb_max
           nproc1=npc*npk*nps*npf*npb
           if (nproc1<nprocmin)     cycle
           if (nproc1>nproc)        cycle
           if (modulo(mband,npb)>0) cycle

!          Change algo if npband>1
           if (optdriver==RUNL_GSTATE.and.npb>1.and. &
&              wf_algo_global==ALGO_NOT_SET) wf_algo=ALGO_DEFAULT_PAR

!          Base speedup
           acc_kgb_0=one;if (npb*npf*nthreads>1) acc_kgb_0=0.7_dp*speedup_fdp(mpw,(npb*npf*nthreads))

           if (npb*npf>4.and.wf_algo==ALGO_LOBPCG_OLD) then
!            Promote npb=npf
             acc_kgb_0=acc_kgb_0*min((one*npf)/(one*npb),(one*npb)/(one*npf))
!            Promote npf<=20
             if (npf>20)then
               acc_kgb_0=acc_kgb_0* &
&                 0.2_dp+(one-0.2_dp)*(sin((pi*(npf-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))) &
&                 /((pi*(npf-NPF_CUTOFF))/(one*(NPFMAX-NPF_CUTOFF))))**2
             end if
           end if

           first_bpp=.true.
           do bpp=bpp_min,bpp_max

             if (wf_algo==ALGO_LOBPCG_NEW) then
               blocksize=npb*bpp;nblocks=mband/blocksize
               if (modulo(bpp,nthreads)>0) cycle
               if ((bpp>1).and.(modulo(bpp,2)>0)) cycle
               if (modulo(mband,npb*bpp)>0) cycle
             else if (wf_algo==ALGO_LOBPCG_OLD) then
               blocksize=npb*bpp;nblocks=mband/blocksize
               if (modulo(mband/npb,bpp)>0) cycle
               if ((bpp>1).and.(modulo(bpp,2)>0)) cycle
               if (one*npb*bpp >max(1.,mband/3.).and.(mband>30)) cycle
               if (npb*npf<=4.and.(.not.first_bpp)) cycle
             else if (wf_algo==ALGO_CHEBFI) then
               !Nothing
             else
               if (bpp/=1.or.npb/=1) cycle
             end if

             first_bpp=.false.

             acc_kgb=acc_kgb_0
!            OLD LOBPCG: promote bpp*npb>mband/3
             if (wf_algo==ALGO_LOBPCG_OLD) then
               if (npb*npf>4.and.mband>30) acc_kgb=acc_kgb*(one-(three*bpp*npb)/(one*mband))
             end if
!            NEW LOBPCG: promote minimal number of blocks
!                        promote block size <= BLOCKSIZE_MAX
             if (wf_algo==ALGO_LOBPCG_NEW) then
               acc_kgb=acc_kgb*(one-0.9_dp*dble(nblocks-1)/dble(mband-1))
               if (blocksize>BLOCKSIZE_MAX) acc_kgb=acc_kgb*max(0.1_dp,one-dble(blocksize)/dble(10*BLOCKSIZE_MAX))
               if (nthreads==1) then
!                Promote npband vs bandpp & npfft
                 if (blocksize>1) acc_kgb=acc_kgb*(0.1_dp*bpp+0.9_dp-blocksize)/(one-blocksize)
                 if (npb*npf>4.and.mband>100) acc_kgb=acc_kgb*(one-0.8_dp*((three*bpp*npb)/(one*mband)-one)**2)
                 tot_ncpus=max(npb,npf);if (tot_ncpus==2) tot_ncpus=0
                 acc_kgb=acc_kgb*(one-0.8_dp*((dble(npb)/dble(npf))-2_dp)**2/(tot_ncpus-2_dp)**2)
                 eff=max(npf,20);acc_kgb=acc_kgb*(one-0.8_dp*min(one,(eff-20)**2))
               end if
             end if

!            CHEBFI: promote npfft=npband and nband>=npfft
             if (wf_algo==ALGO_CHEBFI) then
               if (npf>1) then
                 if (npb>npf) then
                   acc_kgb=acc_kgb*(one-0.8_dp*0.25_dp*((dble(npb)/dble(npf))-one)**2/(nproc1-one)**2)
                 else
                   acc_kgb=acc_kgb*(one-0.8_dp*nproc1**2*((dble(npb)/dble(npf))-one)**2/(nproc1-one)**2)
                 end if
               end if
             end if

!            Resulting "weight"
!            weight0=acc_c*acc_k*acc_s*acc_kgb
             weight0=nproc1*(acc_c+acc_k+acc_s+acc_kgb)/(npc+npk+nps+(npf*npb))

!            Store data
             icount=icount+1
             if (icount<=MAXCOUNT) then
               my_algo(icount)=merge(ALGO_CG,wf_algo,wf_algo==ALGO_NOT_SET)
               my_distp(1:7,icount)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
               weight(icount)=weight0
               if (weight0<weight(imin)) imin=icount
             else
               if (weight0>weight(imin)) then
                 my_algo(imin)=merge(ALGO_CG,wf_algo,wf_algo==ALGO_NOT_SET)
                 my_distp(1:7,imin)=(/npc,npk,nps,npf,npb,bpp,nproc1/)
                 weight(imin)=weight0
                 idum=minloc(weight);imin=idum(1)
               end if
             end if

           end do ! bpp
         end do ! npb
       end do ! npf
     end do ! nps
   end do ! npk
 end do ! npc

!Compute number of selected distributions
 mcount_eff=icount
 mcount=min(mcount_eff,MAXCOUNT)

!Stop if no solution found
 if (mcount==0) then
!  Override here the 0 default value changed in indefo1
   dtset%npimage  = max(1,dtset%npimage)
   dtset%nppert   = max(1,dtset%nppert)
   dtset%npkpt    = max(1,dtset%npkpt)
   dtset%npspinor = max(1,dtset%npspinor)
   dtset%npfft    = max(1,dtset%npfft)
   dtset%npband   = max(1,dtset%npband)
   dtset%bandpp   = max(1,dtset%bandpp)
   write(msg,'(a,i0,2a,i0,a)')  &
&  'Your input dataset does not let Abinit find an appropriate process distribution with nCPUs=',nproc*nthreads,ch10, &
&  'Try to comment all the np* vars and set max_ncpus=',nthreads*nproc,' to have advices on process distribution.'
   MSG_WARNING(msg)
   if (max_ncpus>0) then
     call wrtout(ab_out,msg,'COLL')
     call flush_unit(ab_out)
   end if
   iexit=iexit+1
 end if

!Sort data by increasing weight
 if (mcount>0) then
   ABI_ALLOCATE(isort,(mcount))
   isort=(/(ii,ii=1,mcount)/)
   call sort_dp(mcount,weight,isort,tol6)
   ncount=min(mcount,MAXPRINT)
 end if

!Deduce a global value for paral_kgb
 paral_kgb=dtset%paral_kgb
 if (tread(1)==0) then
   if (any(my_algo(:)/=ALGO_CG)) paral_kgb=1
 end if

!Print output for abipy
 if (iam_master.and.max_ncpus>0.and. &
&    (mcount>0.or.wf_algo_global==ALGO_CG)) then
   write(ount,'(2a)')ch10,"--- !Autoparal"
   if (optdriver==RUNL_GSTATE.and.paral_kgb==0) then
     write(ount,"(a)")"# Autoparal section for GS run (band-by-band CG method)"
   else if (optdriver==RUNL_GSTATE) then
     write(ount,'(a)')'#Autoparal section for GS calculations with paral_kgb'
   else if (optdriver==RUNL_RESPFN) then
     write(ount,'(a)')'#Autoparal section for DFPT calculations'
   else
    msg='Unsupported optdriver'
     MSG_ERROR(msg)
   end if
   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",autoparal
   write(ount,"(a,i0)")"    paral_kgb: ",paral_kgb
   write(ount,"(a,i0)")"    max_ncpus: ",max_ncpus
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ount,"(a,i0)")"    mband: ",mband
   write(ount,"(a)")"configurations:"
   if (optdriver==RUNL_GSTATE.and.paral_kgb==0) then
     work_size = dtset%nkpt * dtset%nsppol
     do ii=1,max_ncpus
       if (ii > work_size) cycle
       do omp_ncpus=1,nthreads
         nks_per_proc = work_size / ii
         nks_per_proc = nks_per_proc + MOD(work_size, ii)
         eff = (one * work_size) / (ii * nks_per_proc)
         write(ount,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
         write(ount,"(a,i0)")"      mpi_ncpus: ",ii
         write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
         write(ount,"(a,f12.9)")"      efficiency: ",eff
         !write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
       end do
     end do
   else if (optdriver==RUNL_GSTATE) then
     omp_ncpus=nthreads
     do jj=mcount,mcount-min(ncount,MAXABIPY)+1,-1
       ii=isort(jj)
       tot_ncpus = my_distp(7,ii)
       eff = weight(jj) / tot_ncpus
       write(ount,'(a,i0)')'    - tot_ncpus: ',tot_ncpus
       write(ount,'(a,i0)')'      mpi_ncpus: ',tot_ncpus
       write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
       write(ount,'(a,f12.9)')'      efficiency: ',eff
       !write(ount,'(a,f12.2)')'      mem_per_cpu: ',mempercpu_mb
       write(ount,'(a)'   )'      vars: {'
       write(ount,'(a,i0,a)')'            npimage: ',my_distp(1,ii),','
       write(ount,'(a,i0,a)')'            npkpt: ', my_distp(2,ii),','
       write(ount,'(a,i0,a)')'            npspinor: ',my_distp(3,ii),','
       write(ount,'(a,i0,a)')'            npfft: ', my_distp(4,ii),','
       write(ount,'(a,i0,a)')'            npband: ',my_distp(5,ii),','
       write(ount,'(a,i0,a)')'            bandpp: ',my_distp(6,ii),','
       write(ount,'(a)')   '            }'
     end do
   else if (optdriver==RUNL_RESPFN) then
     do jj=mcount,mcount-min(ncount,MAXABIPY)+1,-1
       ii=isort(jj)
       tot_ncpus = my_distp(7,ii)
       eff = weight(jj) / tot_ncpus
       write(ount,'(a,i0)')'    - tot_ncpus: ',tot_ncpus
       write(ount,'(a,i0)')'      mpi_ncpus: ',tot_ncpus
       !write(ount,'(a,i0)')'      omp_ncpus: ',omp_ncpus !OMP not supported  (yet)
       write(ount,'(a,f12.9)')'      efficiency: ',eff
       !write(ount,'(a,f12.2)')'      mem_per_cpu: ',mempercpu_mb
       write(ount,'(a)'   )'      vars: {'
       write(ount,'(a,i0,a)')'             nppert: ', my_distp(1,ii),','
       write(ount,'(a,i0,a)')'             npkpt: ', my_distp(2,ii),','
       write(ount,'(a)')   '            }'
      end do
   end if
   write(ount,'(a)')"..."
 end if

!Print out tab with selected choices
 if (mcount>0.and.iam_master) then
   if (nthreads==1) then
     write(msg,'(a,1x,100("="),2a,i0,2a)') ch10,ch10,&
&     ' Searching for all possible proc distributions for this input with #CPUs<=',nthreads*nproc,':',ch10
   else
     write(msg,'(a,1x,100("="),2a,i0,a,i0,2a)')  ch10,ch10,&
&     ' Searching for all possible proc distributions for this input with #CPUs<=',nthreads*nproc,&
&     ' and ',nthreads,' openMP threads:',ch10
   end if
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   !Titles of columns
   msgttl='~'
   if (with_image)  msgttl=trim(msgttl)//'~~~~~~~~~~~'
   if (with_pert)   msgttl=trim(msgttl)//'~~~~~~~~~~~'
   msgttl=trim(msgttl)//'~~~~~~~~~~~~~' ! kpt
   if (with_spinor) msgttl=trim(msgttl)//'~~~~~~~~~~'
   if (with_fft)    msgttl=trim(msgttl)//'~~~~~~~~~~~~~'
   if (with_band)   msgttl=trim(msgttl)//'~~~~~~~~~~~~~'
   if (with_bandpp) msgttl=trim(msgttl)//'~~~~~~~~~~~~~'
   if (with_thread) msgttl=trim(msgttl)//'~~~~~~~~~~'
   msgttl=trim(msgttl)//'~~~~~~~~~~~~~' ! nproc
   if (with_thread) msgttl=trim(msgttl)//'~~~~~~~~~~~~~'
   msgttl=trim(msgttl)//'~~~~~~~~~~~'   ! CPUs
   msgttl=' '//trim(msgttl)
   call wrtout(std_out,msgttl,'COLL');if(max_ncpus>0) call wrtout(ab_out,msgttl,'COLL')
   msg='|'
   if (with_image)  msg=trim(msg)//'   npimage|'
   if (with_pert)   msg=trim(msg)//'    nppert|'
   msg=trim(msg)//'       npkpt|'
   if (with_spinor) msg=trim(msg)//' npspinor|'
   if (with_fft)    msg=trim(msg)//'       npfft|'
   if (with_band)   msg=trim(msg)//'      npband|'
   if (with_bandpp) msg=trim(msg)//'      bandpp|'
   if (with_thread) msg=trim(msg)//' #Threads|'
   msg=trim(msg)//'  #MPI(proc)|'
   if (with_thread) msg=trim(msg)//'       #CPUs|'
   msg=trim(msg)//'    WEIGHT|'
   msg=' '//trim(msg)
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   msg='|'
   write(strg,'(i4,a,i4,a)') npi_min,'<<',npi_max,'|';if (with_image)  msg=trim(msg)//trim(strg)
   write(strg,'(i4,a,i4,a)') npp_min,'<<',npp_max,'|';if (with_pert)   msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i5,a)') npk_min,'<<',npk_max,'|';                 msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i2,a)') nps_min,'<<',nps_max,'|';if (with_spinor) msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i5,a)') npf_min,'<<',npf_max,'|';if (with_fft)    msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i5,a)') npb_min,'<<',npb_max,'|';if (with_band)   msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i5,a)') bpp_min,'<<',bpp_max,'|';if (with_bandpp) msg=trim(msg)//trim(strg)
   write(strg,'(i9,a)'     ) nthreads            ,'|';if (with_thread) msg=trim(msg)//trim(strg)
   write(strg,'(i5,a,i5,a)') 1      ,'<<',nproc  ,'|';                 msg=trim(msg)//trim(strg)
   write(strg,'(i4,a,i6,a)') nthreads,'<<',nthreads*nproc,'|';if (with_thread) msg=trim(msg)//trim(strg)
   write(strg,'(a,i6,a)')   '  <=',nthreads*nproc,'|';                 msg=trim(msg)//trim(strg)
   msg=' '//trim(msg)
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,msgttl,'COLL');if(max_ncpus>0) call wrtout(ab_out,msgttl,'COLL')
   !Loop over selected choices
   do jj=mcount,mcount-ncount+1,-1
     ii=isort(jj)
     msg='|'
     write(strg,'(i10,a)') my_distp(1,ii),'|';if (with_image)  msg=trim(msg)//trim(strg)
     write(strg,'(i10,a)') my_distp(1,ii),'|';if (with_pert)   msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') my_distp(2,ii),'|';                 msg=trim(msg)//trim(strg)
     write(strg,'(i9,a)')  my_distp(3,ii),'|';if (with_spinor) msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') my_distp(4,ii),'|';if (with_fft)    msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') my_distp(5,ii),'|';if (with_band)   msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') my_distp(6,ii),'|';if (with_bandpp) msg=trim(msg)//trim(strg)
     write(strg,'(i9,a)')  nthreads      ,'|';if (with_thread) msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') my_distp(7,ii),'|';                 msg=trim(msg)//trim(strg)
     write(strg,'(i12,a)') nthreads*my_distp(7,ii),'|';if (with_thread) msg=trim(msg)//trim(strg)
     write(strg,'(f10.3,a)') weight(jj)  ,'|';                 msg=trim(msg)//trim(strg)
     msg=' '//trim(msg)
     call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   end do
   !End of tab
   call wrtout(std_out,msgttl,'COLL');if(max_ncpus>0) call wrtout(ab_out,msgttl,'COLL')
   write(msg,'(a,i6,a,i6,a)')' Only the best possible choices for nproc are printed...'
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
 end if ! mcount>0

!Determine an optimal number of bands
 if (optdriver==RUNL_GSTATE.and. &
&    (any(my_algo(1:mcount)==ALGO_LOBPCG_OLD.or. &
&         my_algo(1:mcount)==ALGO_LOBPCG_NEW.or. &
&         my_algo(1:mcount)==ALGO_CHEBFI))) then
   if (mcount>0) then
     icount=isort(mcount)
     npc=my_distp(1,icount);npk=my_distp(2,icount)
     nps=my_distp(3,icount);npf=my_distp(4,icount)
   else
     npc=1;if (with_image ) npc=npi_min
     npk=1;if (with_kpt   ) npk=npk_min
     nps=1;if (with_spinor) nps=nps_min
     npf=1;if (with_fft   ) npf=npf_min
   end if
   nproc1=npc*npk*nps*npf
   msg=ch10//' >>> Possible (best) choices for the number of bands (nband) are:'
   if (with_image.or.with_kpt.or.with_spinor.or.with_fft) msg=trim(msg)//ch10//'     with:'
   write(strg,'(a,i0)') ' npimage=' ,npc;if (with_image)  msg=trim(msg)//trim(strg)
   write(strg,'(a,i0)') ' npkpt='   ,npk;if (with_kpt)    msg=trim(msg)//trim(strg)
   write(strg,'(a,i0)') ' npspinor=',nps;if (with_spinor) msg=trim(msg)//trim(strg)
   write(strg,'(a,i0)') ' npfft='   ,npf;if (with_fft)    msg=trim(msg)//trim(strg)
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   ib1=mband-int(mband*relative_nband_range);if (my_algo(icount)==ALGO_CHEBFI) ib1=mband
   ib2=mband+int(mband*relative_nband_range)
   ABI_ALLOCATE(nproc_best,(1+ib2-ib1))
   ABI_ALLOCATE(nband_best,(1+ib2-ib1))
   nproc_best(:)=1
   nband_best=(/(ii,ii=ib1,ib2)/)
   bpp=merge(1,nthreads,my_algo(icount)==ALGO_CHEBFI)
   do ii=ib1,ib2
     do jj=1,nproc/nproc1
       ibest=1
       do kk=1,jj
         if (mod(jj,kk)/=0) cycle
         if (mod(ii,kk*bpp)==0) ibest=max(ibest,kk)
       end do
       nproc_best(1+ii-ib1)=max(nproc_best(1+ii-ib1),ibest)
     end do
   end do
   call sort_int(1+ib2-ib1,nproc_best,nband_best)
   kk=-1
   do ii=1+ib2-ib1,max(ib2-ib1-MAXBAND_PRINT,1),-1
     write(msg,'(3(a,i6),a,i3,a,i5,a)') '     nband=',nband_best(ii),' using ',nproc1*nproc_best(ii)*nthreads,&
&        ' CPUs =',nproc1*nproc_best(ii),' MPI x',nthreads,' threads (npband=',nproc_best(ii),')'
     call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
     if (nband_best(ii)==mband) kk=nproc_best(ii)
   end do
   if (kk==maxval(nproc_best(:))) then
     if (my_algo(icount)/=ALGO_CHEBFI) then
       write(msg,'(a,i6,a)') ' >>> The present nband value (',mband,') seems to be the best choice!'
     end if
     if (my_algo(icount)==ALGO_CHEBFI) then
       write(msg,'(a,i6,a)') ' >>> The present nband value (',mband,') seems to be a good choice!'
     end if
     call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
   end if
   ABI_DEALLOCATE(nproc_best)
   ABI_DEALLOCATE(nband_best)
 end if

 if (optdriver==RUNL_GSTATE.and.(any(my_algo(1:mcount)==ALGO_CHEBFI))) then
   write(msg,'(5a)') &
&   ' >>> Note that with the "Chebyshev Filtering" algorithm, it is often',ch10,&
&   '     better to increase the number of bands (10% more or a few tens more).',ch10,&
&   '     Advice: increase nband and put nbdbuf input variable to (nband_new-nband_old).'
   call wrtout(std_out,msg,'COLL');if(max_ncpus>0) call wrtout(ab_out,msg,'COLL')
 end if

!Refinement of the process distribution by mean of a LinAlg routines benchmarking
 if (mcount>0.and.optdriver==RUNL_GSTATE.and.autoparal/=1) then
   icount=isort(mcount)
   if (autoparal/=3) then
     if (autoparal==2) then
       write(msg,'(5a,9(a10,a1))') ch10, &
&       ' Values below have been tested with respect to Linear Algebra performance;',ch10,&
&       ' Weights below are corrected according:',ch10,&
&       'npimage','|','npkpt' ,'|','npspinor'  ,'|','npfft'     ,'|','npband','|',' bandpp ' ,'|',&
&       'nproc'  ,'|','weight','|','new weight','|'
     else
       write(msg,'(5a,11(a10,a1))') ch10, &
&       ' Values below have been tested with respect to Linear Algebra performance;',ch10,&
&       ' Weights below are corrected according:',ch10,&
&       'npimage','|','npkpt' ,'|','npspinor'  ,'|','npfft'     ,'|','npband','|',' bandpp ' ,'|',&
&       'nproc'  ,'|','weight','|','new weight','|','best npslk','|','linalggpu' ,'|'
     end if
     call wrtout(std_out,msg,'COLL');if (max_ncpus > 0) call wrtout(ab_out,msg,'COLL')
   end if
   acc_k=zero
   ncount=min(MAXBENCH,mcount);if (autoparal==3) ncount=1
   do jj=mcount,mcount-ncount+1,-1
     ii=isort(jj)
     npf=my_distp(4,ii);npb=my_distp(5,ii);bpp=my_distp(6,ii)
     if ((npb*npf*bpp>1).and.(npf*npb<=mpi_enreg%nproc)) then
       use_linalg_gpu=dtset%use_gpu_cuda
       call compute_kgb_indicator(acc_kgb,bpp,xmpi_world,mband,mpw,npb,npf,np_slk,use_linalg_gpu)
       if (autoparal/=2) then
         my_distp(9,ii)=np_slk
         if (np_slk>0) my_distp(8,ii)=1
!        * gpu_linalg_limit:
!        No use of GPU: htgspw_01.outuge value ~2  *vectsize*blocksize**2 tested
!        Use of GPU:    tiny value ~0.5*vectsize*blocksize**2 tested
         my_distp(10,ii)=2*dtset%mpw*(npb*bpp)**2/npf
         if (use_linalg_gpu==1) my_distp(10,ii)=my_distp(10,ii)/4
       end if
       if (abs(acc_k)<=tol12) acc_k=acc_kgb ! Ref value : the first one computed
!      * Weight (corrected by 10% of the computed ratio)
       weight0=weight(jj)*(one + 0.1_dp*acc_k/acc_kgb)
       if (autoparal==2) then
         write(msg, '(7(i10,a1),f9.2,a2,f9.5,a2)') &
&         my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(3,ii),'|',my_distp(4,ii),'|',&
&         my_distp(5,ii),'|',my_distp(6,ii),'|',my_distp(7,ii),'|',weight(jj),'=>', weight0,' |'
       else if (autoparal==3) then
         write(msg,'(a,5(a,i3))') ch10,' For npband=',npb,', npfft=',npf,' and bandpp=',bpp, &
&         ', compute_kgb_indicator recommends you to set np_slk=',my_distp(9,ii),&
&         ' and use_linalg_gpu=',use_linalg_gpu
       else
         write(msg, '(7(i10,a1),f9.2,a2,f9.5,a2,2(i10,a1))') &
&         my_distp(1,ii),'|',my_distp(2,ii),'|',my_distp(3,ii),'|',my_distp(4,ii),'|',&
&         my_distp(5,ii),'|',my_distp(6,ii),'|',my_distp(7,ii),'|',weight(jj),'=>', weight0,' |',&
&         my_distp(9,ii),'|',use_linalg_gpu,'|'
       end if
       call wrtout(std_out,msg,'COLL');if (max_ncpus>0) call wrtout(ab_out,msg,'COLL')
!      We store the best value in weight(mcount) and keep icount
       if (weight0 > weight(mcount)) then
         icount=ii;weight(mcount)=weight0
       end if
     end if
   end do
 end if

!Store new process distribution
 if (mcount>0.and.max_ncpus<=0) then
   icount=isort(mcount)
   nproc1=my_distp(7,icount)
!  Work load distribution
   if (optdriver==RUNL_GSTATE) then
     dtset%npimage= my_distp(1,icount)
     dtset%nppert = 1
     dtset%npkpt  = my_distp(2,icount)
   end if
   if (optdriver==RUNL_RESPFN) then
     dtset%npimage= 1
     dtset%nppert = my_distp(1,icount)
     dtset%npkpt  = 1
   end if
   dtset%npspinor = my_distp(3,icount)
   dtset%npfft    = my_distp(4,icount)
   dtset%npband   = my_distp(5,icount)
   dtset%bandpp   = my_distp(6,icount)
   if (tread(1)==0)  dtset%paral_kgb= merge(0,1,my_algo(icount)==ALGO_CG)
!  The following lines are mandatory : the DFT+DMFT must use ALL the
!  available procs specified by the user. So nproc1=nproc.
!  Works only if paral_kgb is not activated??
   if (dtset%usedmft/=0.and.optdriver==RUNL_GSTATE) then
     if (dtset%paral_kgb==0) then
       dtset%npspinor = 1 ; dtset%npfft    = 1
       dtset%npband   = 1 ; dtset%bandpp   = 1
       dtset%npimage  = 1
     end if
     nproc1 = nproc
   end if
   if (dtset%npband*dtset%npfft*dtset%bandpp>1) dtset%paral_kgb=1
!  LinAlg parameters: we change values only if they are not present in input file
   if (dtset%paral_kgb==1) then
     if (tread(9)==0) dtset%use_slk=my_distp(8,icount)
     if (tread(10)==0) dtset%np_slk=my_distp(9,icount)
     if (tread(11)==0) dtset%gpu_linalg_limit=my_distp(10,icount)
   end if
!  New definition of "world" MPI communicator
   if (optdriver==RUNL_RESPFN.and.dtset%paral_rf==1) then
     nproc1=max(nproc1,dtset%nsppol*dtset%nkpt) ! Take into account the code in respfn.F90
     nproc1=min(nproc1,nproc)
     nproc1=(nproc1/dtset%nppert)*dtset%nppert
   end if
   call initmpi_world(mpi_enreg,nproc1)
 end if

!Final advice in case max_ncpus > 0
 if (max_ncpus>0.and.mcount>0) then
   write(msg,'(6a)') ch10,&
&   ' Launch a parallel version of ABINIT with a distribution of processors among the above list,',ch10,&
&   ' and the associated input variables (npkpt, npband, npfft, bandpp, etc.).',ch10,&
&   ' The higher weight should be better.'
   call wrtout(std_out,msg,'COLL');if (max_ncpus>0) call wrtout(ab_out,msg,'COLL')
 end if

 if (mcount>0) then
   ABI_DEALLOCATE(isort)
 end if
 ABI_DEALLOCATE(my_distp)
 ABI_DEALLOCATE(my_algo)
 ABI_DEALLOCATE(weight)

!Final line
 write(msg,'(a,100("="),2a)') " ",ch10,ch10
 call wrtout(std_out,msg,'COLL');if (max_ncpus>0) call wrtout(ab_out,msg,'COLL')

!max_ncpus requires a stop
 if (max_ncpus>0) then
   iexit = iexit + 1 ! will stop in the parent.
 end if

 DBG_EXIT("COLL")

 contains

   function speedup_fdp(nn,mm)
   !Expected linear speedup for a nn-sized problem and mm processes
   real(dp) :: speedup_fdp
   integer,intent(in) :: nn,mm
   speedup_fdp=(one*nn)/(one*((nn/mm)+merge(0,1,mod(nn,mm)==0)))
 end function speedup_fdp

end subroutine finddistrproc
!!***

!!****f* ABINIT/compute_kgb_indicator
!! NAME
!! compute_kgb_indicator
!!
!! FUNCTION
!! Only for "KGB" parallelism (LOBPCG algorithm for Ground-state):
!!  Give an indicator of performance for a given distribution of processors
!!  (npband, npfft and bandpp).
!!  Determine best choice of parameters for Scalapack and/or Magma Linear Algebra routines.
!!
!! INPUTS
!!  bandpp=internal lobpcg optimization variable
!!  glb_comm=communicator for global MPI communications
!!  mband=maximum number of bands.
!!  mband=maximum number of plane waves
!!  npband=number of processor 'band'
!!  npfft = number of processor 'fft'
!!  uselinalggpu=indicate if we also test the gpu linear algebra
!!
!! OUTPUT
!!  acc_kgb = indicator of performance
!!  npslk = number of process to used in communicators
!!
!! SIDE EFFECTS
!! This routine can be used to find an indicator in order to refine automatic process distribution.
!!   This indicator is returned in acc_kgb
!! This routine can be used to find the optimal values of np_slk parameter (ScaLapack)
!!   and wheter or not we should use Magma for Linear Algebra in lobpcgwf
!!
!! PARENTS
!!      finddistrproc
!!
!! CHILDREN
!!      abi_linalg_finalize,abi_linalg_init,abi_xhegv,abi_xorthonormalize
!!      wrtout,xmpi_bcast,xmpi_comm_free
!!
!! SOURCE

subroutine compute_kgb_indicator(acc_kgb,bandpp,glb_comm,mband,mpw,npband,npfft,npslk,uselinalggpu)

 use m_abi_linalg

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,glb_comm,mband,mpw,npband,npfft
 integer,intent(inout) :: npslk,uselinalggpu
 real(dp),intent(inout) :: acc_kgb

!Local variables-------------------------------
!scalars
 integer,parameter :: max_number_of_npslk=10,max_number_of_iter=10
 integer :: blocksize,bigorder,ierr,ii,islk,islk1,iter,jj,keep_gpu
 integer :: kgb_comm,my_rank,np_slk,np_slk_max,np_slk_best,nranks
 integer :: use_lapack_gpu,use_slk,vectsize,wfoptalg
 real(dp) :: min_eigen,min_ortho,time_xeigen,time_xortho
 character(len=500) :: message
!arrays
 integer,allocatable :: ranks(:),val_npslk(:)
 real(dp),allocatable :: eigen(:),grama(:,:),gramb(:,:)
 complex(dpc),allocatable :: blockvectorbx(:,:),blockvectorx(:,:),sqgram(:,:)

!******************************************************************

 DBG_ENTER("COLL")

#ifdef DEBUG_MODE
 write(message,'(a,3i3)') 'compute_kgb_indicator : (bpp,npb,npf) = ', bandpp, npband, npfft
 call wrtout(std_out,message,'PERS')
#endif

!Create local communicator for test
 if (xmpi_paral==1) then
   nranks=npfft*npband
   ABI_ALLOCATE(ranks,(nranks))
   ranks=(/((my_rank-1),my_rank=1,nranks)/)
   kgb_comm=xmpi_subcomm(glb_comm,nranks,ranks,my_rank_in_group=my_rank)
   ABI_DEALLOCATE(ranks)
 else
   kgb_comm=xmpi_comm_self
   my_rank=0
 end if

!Only for process in the new subgroup
 if (my_rank/=xmpi_undefined) then

!  We enforce vectsize >=blocksize  (This is not true in lobpcg but
!  these are rare cases and this simplify the matrix constructions below...)
   blocksize=npband*bandpp
   vectsize=max(1+mpw/(npband*npfft),blocksize)
   bigorder=3*blocksize

   ABI_ALLOCATE(blockvectorx,(vectsize,blocksize))
   ABI_ALLOCATE(blockvectorbx,(vectsize,blocksize))
   ABI_ALLOCATE(sqgram,(blocksize,blocksize))
   ABI_ALLOCATE(grama,(2*bigorder,bigorder))
   ABI_ALLOCATE(gramb,(2*bigorder,bigorder))
   ABI_ALLOCATE(eigen,(bigorder))
   ABI_ALLOCATE(val_npslk,(max_number_of_npslk)) ! not too much values tested

   min_eigen=greatest_real
   min_ortho=greatest_real
   np_slk_best=-1 ; np_slk_max=0
#ifdef HAVE_LINALG_SCALAPACK
   np_slk_max=min(max_number_of_npslk,npband*npfft)
#endif

!  Preselect a range of available np_slk values
   val_npslk(1:)=0 ; val_npslk(2)=1
   do islk=3,np_slk_max
     np_slk=val_npslk(islk-1)*2
     do while ((modulo(npband*npfft,np_slk)>0).and.(np_slk<(npband*npfft)))
       np_slk=np_slk+1
     end do
     if(np_slk>(npband*npfft).or.np_slk>mband) exit
     val_npslk(islk)=np_slk
   end do
   np_slk_max=islk-1

!  Loop over np_slk values
   islk1=1
#ifdef HAVE_LINALG_MAGMA
   islk1=1-uselinalggpu
#endif
   do islk=islk1,np_slk_max

     time_xortho=zero ; time_xeigen=zero

     use_slk=0
     if (islk==0) then
!      This is the test for the GPU
       use_lapack_gpu=1 ; np_slk=0
     else
       use_lapack_gpu=0 ; np_slk=val_npslk(islk)
       if (np_slk>0) use_slk=1
     end if

!    Initialize linalg parameters for this np_slk value
!    For the first np_slk value, everything is initialized
!    For the following np_slk values, only Scalapack parameters are updated
     wfoptalg=14 ! Simulate use of LOBPCG
     call abi_linalg_init(bigorder,RUNL_GSTATE,wfoptalg,1,&
&                         use_lapack_gpu,use_slk,np_slk,kgb_comm)

!    We could do mband/blocksize iter as in lobpcg but it's too long
     do iter=1,max_number_of_iter

!      Build matrixes
       do ii=1,vectsize
         do jj=1,blocksize
           if (ii>jj) then
             blockvectorx(ii,jj) =czero
             blockvectorbx(ii,jj)=czero
           else
             blockvectorx(ii,jj) =cone
             blockvectorbx(ii,jj)=cone
           end if
         end do
       end do
       grama=zero;gramb=zero
       do jj=1,bigorder
         do ii=jj,bigorder
           if (ii==jj) then
             grama(2*ii-1,jj)=one
             gramb(2*ii-1,jj)=one
           else
             grama(2*ii-1:2*ii,jj)=one
             grama(2*jj-1,ii)= one
             grama(2*jj  ,ii)=-one
           end if
         end do
       end do

!      Call to abi_xorthonormalize
       time_xortho=time_xortho-abi_wtime()
       call abi_xorthonormalize(blockvectorx,blockvectorbx,blocksize,kgb_comm,sqgram,vectsize)
       time_xortho = time_xortho + abi_wtime()

!      Call to abi_xhegv
       time_xeigen=time_xeigen-abi_wtime()
       call abi_xhegv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,&
&       x_cplx=2,use_slk=use_slk,use_gpu=use_lapack_gpu)
       time_xeigen=time_xeigen+abi_wtime()

     end do ! iter

!    Finalize linalg parameters for this np_slk value
!    For the last np_slk value, everything is finalized
!    For the previous np_slk values, only Scalapack parameters are updated
     call abi_linalg_finalize()

     time_xortho= time_xortho*mband/blocksize
     time_xeigen= time_xeigen*mband/blocksize
     if (time_xortho<min_ortho) min_ortho=time_xortho
     if (time_xeigen<min_eigen) then
       min_eigen=time_xeigen
       np_slk_best=np_slk
       keep_gpu=use_lapack_gpu
     end if

   end do ! np_slk

#ifdef DEBUG_MODE
   write(message,'(2(a,es15.3),a,i3)') ' In the best case, xortho took ',min_ortho,&
&   ' and xeigen took ',min_eigen,' for np_slk=',np_slk_best
   call wrtout(std_out,message,'PERS')
#endif

!  Final values to be sent to others process
   acc_kgb=min_ortho+four*min_eigen
   npslk=max(np_slk_best,1)
   uselinalggpu=keep_gpu

   ABI_DEALLOCATE(blockvectorx)
   ABI_DEALLOCATE(blockvectorbx)
   ABI_DEALLOCATE(sqgram)
   ABI_DEALLOCATE(grama)
   ABI_DEALLOCATE(gramb)
   ABI_DEALLOCATE(eigen)
   ABI_DEALLOCATE(val_npslk)

 end if ! my_rank in group

!Free local MPI communicator
 call xmpi_comm_free(kgb_comm)

!Broadcast of results to be sure every process has them
 call xmpi_bcast(acc_kgb,0,glb_comm,ierr)
 call xmpi_bcast(npslk,0,glb_comm,ierr)
 call xmpi_bcast(uselinalggpu,0,glb_comm,ierr)

#ifndef DEBUG_MODE
 ABI_UNUSED(message)
#endif

 DBG_EXIT("COLL")

end subroutine compute_kgb_indicator
!!***

end module m_mpi_setup
!!***
