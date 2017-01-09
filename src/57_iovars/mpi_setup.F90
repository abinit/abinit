!{\src2tex{textfont=tt}}
!!****f* ABINIT/mpi_setup
!! NAME
!! mpi_setup
!!
!! FUNCTION
!! Big loop on the datasets :
!! - compute mgfft,mpw,nfft,... for this data set ;
!! - fill mpi_enreg
!!  *** At the output of this routine, all the dtsets input variables are known ***
!! The content of dtsets should not be modified anymore afterwards.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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
!!   mpi_enregs=informations about MPI parallelization
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mpi_setup(dtsets,filnam,lenstr,mpi_enregs,ndtset,ndtset_alloc,string)

 use defs_basis
 use defs_abitypes
 use defs_parameters
 use m_distribfft
 use m_xmpi
 use m_errors
 use m_profiling_abi

 use m_fftcore,      only : fftalg_for_npfft
 use m_mpinfo,       only : init_mpi_enreg,mpi_distrib_is_ok
 use m_libpaw_tools, only : libpaw_write_comm_set

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mpi_setup'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_42_parser
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
 use interfaces_56_recipspace
 use interfaces_57_iovars, except_this_one => mpi_setup
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lenstr,ndtset,ndtset_alloc
 type(MPI_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
 character(len=*),intent(inout) :: string
!arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------------
!scalars
 integer :: blocksize,exchn2n3d,iband,idtset,iexit,ii,iikpt,iikpt_modulo
 integer :: isppol,jdtset,marr,mband_upper
 integer :: me_fft,mgfft,mgfftdg,mkmem,mpw,mpw_k,optdriver
 integer :: nfft,nfftdg,nkpt,nkpt_me,npert,nproc,nproc_fft,nqpt
 integer :: nspink,nsppol,nsym,paral_fft,response,tnband,tread0,usepaw,vectsize
 integer :: fftalg,fftalga,fftalgc
 logical :: fftalg_read,ortalg_read,wfoptalg_read,do_check
 real(dp) :: dilatmx,ecut,ecut_eff,ecutdg_eff,ucvol
 character(len=500) :: message
!arrays
 integer :: ngfft(18),ngfftdg(18),ngfftc(3),tread(11)
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

!  Read parallel input parameters
   marr=max(5,dtsets(idtset)%npsp,dtsets(idtset)%nimage)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
   nkpt  =dtsets(idtset)%nkpt
   nsppol=dtsets(idtset)%nsppol
   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=0
   usepaw=dtsets(idtset)%usepaw
   mband_upper=maxval(dtsets(idtset)%nband(1:nkpt*nsppol))

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
&     'Action : modify compilation option or paral_kgb in the input file.'
     MSG_WARNING(message)
   end if

   if ( ALL(optdriver /= [RUNL_GSTATE, RUNL_RESPFN, RUNL_GWLS]) .and. dtsets(idtset)%paral_kgb/=0) then
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
   if (tread0==1.and.optdriver==RUNL_RESPFN) dtsets(idtset)%paral_rf=intarr(1)

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

   ! Dump the list of irreducible perturbations and exit.
   if (dtsets(idtset)%paral_rf==-1) then
     call get_npert_rbz(dtsets(idtset),nband_rbz,nkpt_rbz,npert)
     ABI_DEALLOCATE(nband_rbz)
     ABI_DEALLOCATE(nkpt_rbz)
     iexit = iexit + 1
   end if

!  From total number of procs, compute all possible distributions
!  Ignore exit flag if GW calculations because autoparal section is performed in screening/sigma/bethe_salpeter
   call finddistrproc(dtsets,filnam,idtset,iexit,mband_upper,mpi_enregs(idtset),ndtset_alloc,tread)
   if (any(optdriver == [RUNL_SCREENING, RUNL_SIGMA, RUNL_BSE])) iexit = 0

   if ((optdriver/=RUNL_GSTATE.and.optdriver/=RUNL_GWLS).and. &
&   (dtsets(idtset)%npkpt/=1   .or.dtsets(idtset)%npband/=1.or.dtsets(idtset)%npfft/=1.or. &
&   dtsets(idtset)%npspinor/=1.or.dtsets(idtset)%bandpp/=1)) then
!&   .or.(dtsets(idtset)%iscf<0)) then
     dtsets(idtset)%npkpt=1 ; dtsets(idtset)%npspinor=1 ; dtsets(idtset)%npfft=1
     dtsets(idtset)%npband=1; dtsets(idtset)%bandpp=1
     dtsets(idtset)%paral_kgb=0
     message = 'For non ground state calculation, set bandpp, npfft, npband, npspinor and npkpt to 1'
     MSG_WARNING(message)
   end if

!  Read again some input data to take into account a possible change of paral_kgb
   wfoptalg_read=.false.
   call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'wfoptalg',tread0,'INT')
   if(tread0==1) then
     dtsets(idtset)%wfoptalg=intarr(1)
     wfoptalg_read=.true.
   else
     if (dtsets(idtset)%usepaw==0) dtsets(idtset)%wfoptalg=0
     if (dtsets(idtset)%usepaw/=0) dtsets(idtset)%wfoptalg=10
     if ((optdriver==RUNL_GSTATE.or.optdriver==RUNL_GWLS).and.dtsets(idtset)%paral_kgb/=0) dtsets(idtset)%wfoptalg=14
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

   mpi_enregs(idtset)%paral_kgb=dtsets(idtset)%paral_kgb

   call initmpi_img(dtsets(idtset),mpi_enregs(idtset),-1)

!  Cycle if the processor is not used
   if (mpi_enregs(idtset)%me<0) then
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     cycle
   end if

   response=0
   if (dtsets(idtset)%rfddk/=0 .or. dtsets(idtset)%rf2_dkdk/=0 .or. dtsets(idtset)%rf2_dkde/=0 .or. &
&   dtsets(idtset)%rfelfd/=0 .or. dtsets(idtset)%rfphon/=0 .or. dtsets(idtset)%rfstrs/=0 .or. &
&   dtsets(idtset)%rfuser/=0) response=1

   nproc=mpi_enregs(idtset)%nproc_cell

!  --IF CUDA AND RECURSION:ONLY BAND PARALLELISATION
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

!      mkmem was not set in the input file so default to incore solution
       write(message,'(6a)') &
&       'mpi_setup: ',nm_mkmem(ii),' undefined in the input file.','Use default ',nm_mkmem(ii),' = nkpt'
       call wrtout(std_out,message,'COLL')
       mkmem=nkpt
     end if

!    Check whether nkpt distributed on the processors <= mkmem;
!    if so then may run entirely in core,
!    avoiding i/o to disk for wavefunctions and kg data.
!    mkmem/=0 to avoid i/o; mkmem==0 to use disk i/o for nkpt>=1.
     if (nkpt_me<=mkmem .and. mkmem/=0 ) then
       write(message, '(a,i0,a,a,a,i0,a)' ) &
&       ' mpi_setup: With nkpt_me=',nkpt_me,' and ',nm_mkmem(ii),' = ',mkmem,', ground state wf handled in core.'
       call wrtout(std_out,message,'COLL')
       if(nkpt_me<mkmem .and. nkpt_me/=0)then
         write(message,'(3a)')' Resetting ',nm_mkmem(ii),' to nkpt_me to save memory space.'
         mkmem=nkpt_me
         call wrtout(std_out,message,'COLL')
       end if
     else if(mkmem/=0)then
       write(message, '(a,i0,3a,i0,5a)' ) &
&       ' mpi_setup: With nkpt_me=',nkpt_me,'and ',nm_mkmem(ii),' = ',mkmem,&
&       ' ground state wf require disk i/o.',ch10,&
&       ' Resetting ',nm_mkmem(ii),' to zero to save memory space.'
       mkmem=0
       call wrtout(std_out,message,'COLL')
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
   do_check = all(optdriver /= [RUNL_SCREENING, RUNL_SIGMA, RUNL_BSE])
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
   MSG_ERROR_NODUMP("aborting now")
 end if

 DBG_EXIT("COLL")

end subroutine mpi_setup
!!***
